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

foreach ($throat in $tValues) {
    foreach ($exterior in $dValues) {
        $length = $throat + $exterior
        foreach ($omega in $omegaValues) {
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
    analytic_lipschitz_bound = '|dG/dd| <= 1/2'
    kernel_operator_bound = '||Kbar|| <= T^2/2 by Schur test'
    interacting_remainder_positive_proved = $false
    shrinking_cutoff_uniformity_proved = $false
    physical_matching_proved = $false
    lean_invocations = 0
}
$outputPath = Join-Path $root $Output
$result | ConvertTo-Json -Depth 8 | Set-Content -LiteralPath $outputPath -Encoding utf8
Write-Output (ConvertTo-Json @{ status = 'passed'; checked_cases = $checked; max_ratio = $maxRatio })

