$ErrorActionPreference = 'Stop'
$rng = [System.Random]::new(20261012)
$count = 10000
$maxResidual = 0.0
for ($i = 0; $i -lt $count; $i++) {
    $q = 1.0 + 999.0 * $rng.NextDouble()
    $G = 1e-4 + 0.0999 * $rng.NextDouble()
    $re = 0.5 + 4.5 * $rng.NextDouble()
    $A = $re * $re * $re / $G
    $alpha = $q * $G / (4.0 * [Math]::PI * $re * $re)
    $ell = 16.0 * $re * $re * $re / ($q * $G)
    $xStar = 1.0 / $ell

    $slopeFromMMP = [Math]::PI * $alpha / 2.0
    $slopeFromRN = 2.0 * $re / $ell
    $massFromJT = -$A * $xStar * $xStar / 2.0
    $massFromParameters = -$q * $q * $G / (512.0 * $re * $re * $re)
    $residual = [Math]::Max([Math]::Abs($slopeFromMMP - $slopeFromRN), [Math]::Abs($massFromJT - $massFromParameters))
    $scale = [Math]::Max(1.0, [Math]::Max([Math]::Abs($massFromJT), [Math]::Abs($slopeFromMMP)))
    $relativeResidual = $residual / $scale
    if ($relativeResidual -gt $maxResidual) { $maxResidual = $relativeResidual }
    if ($relativeResidual -gt 1e-12) { throw "RN/JT match residual $relativeResidual at sample ${i}" }
}
Write-Output "PASS: $count random parameter samples; max scaled residual $maxResidual."
