$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.Numerics

$rng = [System.Random]::new(20261012)
$count = 10000
$maxResidual = 0.0
for ($i = 0; $i -lt $count; $i++) {
    $x = 0.5 + 2.0 * $rng.NextDouble()
    $epsilon = 0.01 + 0.2 * $rng.NextDouble()
    $arg = $epsilon * $x
    if ($arg -ge 0.25) { $epsilon = 0.24 / $x; $arg = 0.24 }
    $a = [Math]::Acos($arg)
    $b = [Math]::Sqrt(1.0 - $arg * $arg)
    $r = $epsilon / $b
    $T = 2.0 * $a / $x
    $D = 0.01 + 5.0 * $rng.NextDouble()
    $L = $T + $D
    $omega = [System.Numerics.Complex]::new(-2.0 + 4.0 * $rng.NextDouble(), 0.02 + 2.0 * $rng.NextDouble())
    $imag = [System.Numerics.Complex]::ImaginaryOne
    $E = [System.Numerics.Complex]::Exp($imag * $omega * $T)
    $Z = [System.Numerics.Complex]::Exp($imag * $omega * $L)
    $one = [System.Numerics.Complex]::One
    $psi = $one
    $phi = -$imag * $omega * $psi
    $deltaP = $psi * ($one - $E) + $r * $phi * ($one + $E)
    $deltaR = $deltaP / ($one - $Z)
    $etaDirect = $x * ($deltaR - $psi - $r * $phi)
    $etaTransfer = $x * ($imag / $omega * ($Z - $E) / ($one - $Z) + $r * ($E + $Z) / ($one - $Z)) * $phi
    $residual = [System.Numerics.Complex]::Abs($etaDirect - $etaTransfer)
    if ($residual -gt $maxResidual) { $maxResidual = $residual }
    if ($residual -gt 1e-9) { throw "Boundary-map transfer residual $residual at sample $i" }
}

Write-Output "PASS: $count finite-cut complex-frequency transfer samples; max absolute residual $maxResidual."
