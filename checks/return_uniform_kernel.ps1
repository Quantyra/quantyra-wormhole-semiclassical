param([string]$Output = "evidence/return-uniform-kernel.json")
$ErrorActionPreference = 'Stop'
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path

function OneMinusExpNegative([double]$x) {
    if ($x -lt 1e-5) {
        $term = $x
        $sum = $term
        for ($n = 2; $n -le 12; $n++) {
            $term *= -$x / $n
            $sum += $term
        }
        return $sum
    }
    return 1.0 - [math]::Exp(-$x)
}

$tValues = @(0.25, 0.1, 0.01, 1e-4 | ForEach-Object { 2.0 * [math]::Acos($_) })
$dValues = @(-8, -6, -4, -2, 0, 2, 4, 6, 8 | ForEach-Object { [math]::Pow(10.0, $_) })
$omegaValues = @(-10, -8, -6, -4, -2, 0, 2, 4, 6, 8 | ForEach-Object { [math]::Pow(10.0, $_) })
$fractions = @(0..20 | ForEach-Object { $_ / 20.0 })
$checked = 0
$maxRatio = 0.0
$maxDrop = 0.0
$tMin = 2.0 * [math]::Acos(0.25)
$c0 = 0.5 * [math]::Cosh($tMin / 2.0) / [math]::Sinh($tMin / 2.0)
$diagonalCases = 0
$maxDiagonalRatio = 0.0
$remainderCases = 0
$profileCases = 0
$maxVertexRatio = 0.0
$maxProfileIntegralRatio = 0.0

foreach ($throat in $tValues) {
    foreach ($exterior in $dValues) {
        $length = $throat + $exterior
        foreach ($omega in $omegaValues) {
            if ($omega -ge 1.0) {
                $tail = [math]::Exp(-$omega * $length)
                $freeDiag = 1.0 / (2.0 * $omega) + $tail / ($omega * (1.0 - $tail))
                $diagonalRatio = $omega * $freeDiag / $c0
                if ($diagonalRatio -gt (1.0 + 2e-12)) {
                    throw "Free diagonal bound failed: omega=$omega L=$length ratio=$diagonalRatio"
                }
                $maxDiagonalRatio = [math]::Max($maxDiagonalRatio, $diagonalRatio)
                $freeRemainder = $tail / ($omega * (1.0 - $tail))
                $uniformRemainder = [math]::Exp(-$omega * $tMin) / ($omega * (1.0 - [math]::Exp(-$tMin)))
                if ($freeRemainder -gt ($uniformRemainder * (1.0 + 2e-12))) {
                    throw "Free return remainder bound failed: omega=$omega L=$length"
                }
                $diagonalCases++
                $remainderCases++
            }
            foreach ($fraction in $fractions) {
                $separation = $throat * $fraction
                $circular = [math]::Min($separation, $length - $separation)
                $n1 = OneMinusExpNegative ($omega * $circular)
                $n2 = OneMinusExpNegative ($omega * ($length - $circular))
                $den = 2.0 * $omega * (OneMinusExpNegative ($omega * $length))
                $drop = $n1 * $n2 / $den
                $bound = $circular / 2.0
                if ($drop -gt ($bound * (1.0 + 3e-12) + 1e-300)) {
                    throw "Kernel bound failed: omega=$omega L=$length d=$circular drop=$drop bound=$bound"
                }
                if ($bound -gt 0) { $maxRatio = [math]::Max($maxRatio, $drop / $bound) }
                $maxDrop = [math]::Max($maxDrop, $drop)
                $checked++
            }
        }
    }
}

foreach ($epsilon in @(0.25, 0.2, 0.1, 0.01, 1e-4, 1e-8)) {
    $A = [math]::Acos($epsilon)
    $b = [math]::Sqrt(1.0 - $epsilon * $epsilon)
    $kappa = -$epsilon / ($A * $b)
    for ($j = 0; $j -le 100; $j++) {
        $z = -$A + 2.0 * $A * $j / 100.0
        $g = 1.0 + $kappa + 2.0 * $kappa * $z * [math]::Tan($z)
        $ratio = [math]::Abs($g) / (4.0 / 3.0)
        if ($ratio -gt (1.0 + 2e-12)) { throw "Moving-mouth vertex envelope failed: epsilon=$epsilon z=$z g=$g" }
        $maxVertexRatio = [math]::Max($maxVertexRatio, $ratio)
        $profileCases++
    }
    $exactSec4 = 2.0 * ($b / $epsilon + $b * $b * $b / (3.0 * $epsilon * $epsilon * $epsilon))
    $profileBound = (16.0 / 9.0) * $exactSec4
    $asymptoticBound = (38.0 / 27.0) / ($epsilon * $epsilon * $epsilon)
    $integralRatio = $profileBound / $asymptoticBound
    if ($integralRatio -gt (1.0 + 2e-12)) { throw "Moving-mouth L2 bound failed: epsilon=$epsilon ratio=$integralRatio" }
    $maxProfileIntegralRatio = [math]::Max($maxProfileIntegralRatio, $integralRatio)
}

function Get-NormalizedSha256([string]$path) {
    $bytes = [System.IO.File]::ReadAllBytes($path)
    $text = [System.Text.Encoding]::UTF8.GetString($bytes).Replace("`r`n", "`n")
    $normalized = [System.Text.Encoding]::UTF8.GetBytes($text)
    $sha = [System.Security.Cryptography.SHA256]::Create()
    try { return ([BitConverter]::ToString($sha.ComputeHash($normalized))).Replace('-', '').ToLowerInvariant() }
    finally { $sha.Dispose() }
}

$notePath = Join-Path $root 'notes/return-uniform-low-frequency-kernel.md'
$scriptPath = Join-Path $root 'checks/return_uniform_kernel.ps1'
$result = [ordered]@{
    date = '2026-10-12'
    status = 'passed'
    script_sha256 = Get-NormalizedSha256 $scriptPath
    input_sha256 = @{ 'notes/return-uniform-low-frequency-kernel.md' = Get-NormalizedSha256 $notePath }
    powershell = $PSVersionTable.PSVersion.ToString()
    checked_cases = $checked
    throat_lengths = $tValues
    exterior_returns = $dValues
    euclidean_frequencies = $omegaValues
    separation_grid_points_per_geometry = $fractions.Count
    max_drop_over_circular_separation_half = $maxRatio
    max_absolute_drop = $maxDrop
    free_diagonal_bound_cases = $diagonalCases
    max_omega_times_free_diagonal_over_C0 = $maxDiagonalRatio
    free_return_remainder_bound_cases = $remainderCases
    C0 = $c0
    moving_mouth_vertex_profile_cases = $profileCases
    max_vertex_envelope_ratio = $maxVertexRatio
    max_profile_integral_bound_ratio = $maxProfileIntegralRatio
    analytic_lipschitz_bound = '|dG/dd| <= 1/2'
    kernel_operator_bound = '||Kbar|| <= T^2/2 by Schur test'
    interacting_remainder_positive_proved = $false
    potential_seagull_low_frequency_uniform = $true
    high_frequency_interaction_correction_uniform = $true
    far_ultraviolet_potential_bubble_tail_bound = '38 C a^2 epsilon / 27, for fixed external frequency and epsilon^-2 >= 2|w|+1'
    moving_mouth_bubble_midband_controlled = $false
    full_response_shrinking_cutoff_uniformity_proved = $false
    shrinking_cutoff_uniformity_proved = $false
    physical_matching_proved = $false
    lean_invocations = 0
}
$outputPath = Join-Path $root $Output
$result | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $outputPath -Encoding utf8
Write-Output (ConvertTo-Json @{ status = 'passed'; checked_cases = $checked; max_ratio = $maxRatio })
