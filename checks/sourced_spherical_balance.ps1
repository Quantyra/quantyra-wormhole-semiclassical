$ErrorActionPreference = 'Stop'

# Check the trace-reversed reduced source against the differentiated charged
# Misner-Sharp mass for unrestricted scalar 2D jets. This is a numerical
# consistency check of the identity derived in the accompanying note.
$rng = [System.Random]::new(20261012)
$maxResidual = 0.0
$count = 10000
for ($i = 0; $i -lt $count; $i++) {
    $R = 0.5 + 4.0 * $rng.NextDouble()
    $a = 0.1 + 3.0 * $rng.NextDouble()
    $Y = -2.0 + 4.0 * $rng.NextDouble()
    $boxR = -3.0 + 6.0 * $rng.NextDouble()
    $Ra = -2.0 + 4.0 * $rng.NextDouble()
    $Hcontract = -3.0 + 6.0 * $rng.NextDouble()
    $C = $Y - 1.0 + $a / ($R * $R)

    $massGradientTimes2G = $Ra * (1.0 - $Y - $a / ($R * $R)) - 2.0 * $R * $Hcontract
    $traceReversedSource = 2.0 * $R * ($boxR * $Ra - $Hcontract) + $C * $Ra - (2.0 * $R * $boxR + 2.0 * $C) * $Ra
    $residual = [Math]::Abs($massGradientTimes2G - $traceReversedSource)
    if ($residual -gt $maxResidual) { $maxResidual = $residual }
    if ($residual -gt 1e-10) { throw "Identity residual $residual at sample $i" }
}

Write-Output "PASS: $count arbitrary-jet samples; max absolute residual $maxResidual."
