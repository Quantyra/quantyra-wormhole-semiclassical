$ErrorActionPreference = 'Stop'
$rng = [System.Random]::new(20261012)
$count = 10000
$maxSpatial = 0.0
$maxMomentum = 0.0
$maxCasimir = 0.0
for ($i = 0; $i -lt $count; $i++) {
    $x = -1.2 + 2.4 * $rng.NextDouble()
    $A = -2.0 + 4.0 * $rng.NextDouble()
    $B = -2.0 + 4.0 * $rng.NextDouble()
    $C = -2.0 + 4.0 * $rng.NextDouble()
    $D = -2.0 + 4.0 * $rng.NextDouble()
    $E = -2.0 + 4.0 * $rng.NextDouble()
    $F = -2.0 + 4.0 * $rng.NextDouble()
    $C1 = -1.0 + 2.0 * $rng.NextDouble()
    $C2 = -1.0 + 2.0 * $rng.NextDouble()
    $C3 = -1.0 + 2.0 * $rng.NextDouble()
    $sinx = [Math]::Sin($x); $cosx = [Math]::Cos($x)
    $tanx = $sinx / $cosx; $secx = 1.0 / $cosx

    # H=A+B sin(x)+C cos(x), with its two integrals from 0 to x.
    $H = $A + $B * $sinx + $C * $cosx
    $P = $A * $x + $B * (1.0 - $cosx) + $C * $sinx
    $Q = $A * (1.0 - $cosx) + $B * ($x / 2.0 - [Math]::Sin(2.0 * $x) / 4.0) + $C * $sinx * $sinx / 2.0
    $phi = -$tanx * $P + $secx * $Q + $C1 * $tanx + $C2 * $secx
    $phiPrime = -$secx * $secx * $P + $secx * $tanx * $Q + $C1 * $secx * $secx + $C2 * $secx * $tanx
    $phiSecond = -2.0 * $secx * $secx * $tanx * $P - $H + $secx * ($tanx * $tanx + $secx * $secx) * $Q + 2.0 * $C1 * $secx * $secx * $tanx + $C2 * ($secx * $tanx * $tanx + $secx * $secx * $secx)
    $spatialResidual = $phiSecond - $tanx * $phiPrime - $secx * $secx * $phi + $H

    # J=D+E sin(x)+F cos(x); differentiate the integrated first-order solution.
    $J = $D + $E * $sinx + $F * $cosx
    $W = $D * $sinx + $E * $sinx * $sinx / 2.0 + $F * ($x / 2.0 + [Math]::Sin(2.0 * $x) / 4.0)
    $pi = $secx * ($C3 + $W)
    $piPrime = $secx * $tanx * ($C3 + $W) + $secx * $cosx * $J
    $momentumResidual = $piPrime - $tanx * $pi - $J

    $maxSpatial = [Math]::Max($maxSpatial, [Math]::Abs($spatialResidual))
    $maxMomentum = [Math]::Max($maxMomentum, [Math]::Abs($momentumResidual))

    $phiHomogeneous = $C1 * $tanx + $C2 * $secx
    $piHomogeneous = $C3 * $secx
    $phiPrimeHomogeneous = $C1 * $secx * $secx + $C2 * $secx * $tanx
    $casimir = $phiHomogeneous * $phiHomogeneous + $cosx * $cosx * $piHomogeneous * $piHomogeneous - $cosx * $cosx * $phiPrimeHomogeneous * $phiPrimeHomogeneous
    $expectedCasimir = $C2 * $C2 + $C3 * $C3 - $C1 * $C1
    $maxCasimir = [Math]::Max($maxCasimir, [Math]::Abs($casimir - $expectedCasimir))
    if ([Math]::Abs($spatialResidual) -gt 1e-9 -or [Math]::Abs($momentumResidual) -gt 1e-9) {
        throw "JT Cauchy residual at sample ${i}: spatial=$spatialResidual momentum=$momentumResidual"
    }
}
if ($maxCasimir -gt 1e-10) { throw "JT Casimir residual $maxCasimir" }
Write-Output "PASS: $count source-profile samples; max spatial residual $maxSpatial; max momentum residual $maxMomentum; max homogeneous Casimir residual $maxCasimir."
