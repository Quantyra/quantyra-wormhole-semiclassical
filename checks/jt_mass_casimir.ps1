$ErrorActionPreference = 'Stop'
$rng = [System.Random]::new(51012)
$maxIdentity = 0.0
$maxBalanceBridge = 0.0
$maxExpansionScaled = 0.0
$samples = 10000

for ($i = 0; $i -lt $samples; $i++) {
    $phi = 2.0 * $rng.NextDouble() - 1.0
    $p0 = 2.0 * $rng.NextDouble() - 1.0
    $p1 = 2.0 * $rng.NextDouble() - 1.0
    $h00 = 2.0 * $rng.NextDouble() - 1.0
    $h01 = 2.0 * $rng.NextDouble() - 1.0
    $h11 = 2.0 * $rng.NextDouble() - 1.0
    $box = -$h00 + $h11
    $trace = $box - 2.0 * $phi
    $S00 = -$box - $h00 + $phi
    $S01 = -$h01
    $S11 = $box - $h11 - $phi
    $rhs0 = 2.0 * (($S00 + $trace) * (-$p0) + $S01 * $p1)
    $rhs1 = 2.0 * ($S01 * (-$p0) + ($S11 - $trace) * $p1)
    $lhs0 = 2.0 * ($phi * $p0 - ($h00 * (-$p0) + $h01 * $p1))
    $lhs1 = 2.0 * ($phi * $p1 - ($h01 * (-$p0) + $h11 * $p1))
    $maxIdentity = [Math]::Max($maxIdentity, [Math]::Abs($lhs0 - $rhs0))
    $maxIdentity = [Math]::Max($maxIdentity, [Math]::Abs($lhs1 - $rhs1))

    # Expand the exact spherical mass balance and compare it with the JT
    # Casimir flux using S_ab=(2G/r_e^2) T_ab.
    $re = 0.1 + 10.0 * $rng.NextDouble()
    $grav = 0.01 + 2.0 * $rng.NextDouble()
    $t00 = 2.0 * $rng.NextDouble() - 1.0
    $t01 = 2.0 * $rng.NextDouble() - 1.0
    $t11 = 2.0 * $rng.NextDouble() - 1.0
    $ttrace = -$t00 + $t11
    $s00 = 2.0 * $grav / ($re * $re) * $t00
    $s01 = 2.0 * $grav / ($re * $re) * $t01
    $s11 = 2.0 * $grav / ($re * $re) * $t11
    $strace = -$s00 + $s11
    $casimir0 = $re / (4.0 * $grav) * (($s00 + $strace) * (-$p0) + $s01 * $p1)
    $casimir1 = $re / (4.0 * $grav) * ($s01 * (-$p0) + ($s11 - $strace) * $p1)
    $spherical0 = 1.0 / (2.0 * $re) * (($t00 + $ttrace) * (-$p0) + $t01 * $p1)
    $spherical1 = 1.0 / (2.0 * $re) * ($t01 * (-$p0) + ($t11 - $ttrace) * $p1)
    $maxBalanceBridge = [Math]::Max($maxBalanceBridge, [Math]::Abs($casimir0 - $spherical0))
    $maxBalanceBridge = [Math]::Max($maxBalanceBridge, [Math]::Abs($casimir1 - $spherical1))

    $eps = 0.001
    $x = $eps * (2.0 * $rng.NextDouble() - 1.0)
    $q0 = $eps * (2.0 * $rng.NextDouble() - 1.0)
    $q1 = $eps * (2.0 * $rng.NextDouble() - 1.0)
    $grad2 = -$q0 * $q0 + $q1 * $q1
    $exact = 0.5 * ([Math]::Sqrt(1.0 + $x) + 1.0 / [Math]::Sqrt(1.0 + $x) - $grad2 / (4.0 * [Math]::Sqrt(1.0 + $x)))
    $quadratic = 1.0 + ($x * $x - $grad2) / 8.0
    $scaled = [Math]::Abs($exact - $quadratic) / ($eps * $eps * $eps)
    $maxExpansionScaled = [Math]::Max($maxExpansionScaled, $scaled)
}

$maxNormalization = 0.0
for ($i = 0; $i -lt 1000; $i++) {
    $re = 0.1 + 10.0 * $rng.NextDouble()
    $G = 0.01 + 2.0 * $rng.NextDouble()
    $ell = 0.2 + 8.0 * $rng.NextDouble()
    $beta = 2.0 * $re / $ell
    $lhs = $re / (8.0 * $G) * (-$beta * $beta)
    $rhs = -$re * $re * $re / (2.0 * $G * $ell * $ell)
    $maxNormalization = [Math]::Max($maxNormalization, [Math]::Abs($lhs - $rhs) / [Math]::Max(1.0, [Math]::Abs($rhs)))
}

if ($maxIdentity -gt 1e-12 -or $maxBalanceBridge -gt 1e-12 -or $maxExpansionScaled -gt 1.0 -or $maxNormalization -gt 1e-12) {
    throw "JT Casimir check failed: identity=$maxIdentity balance_bridge=$maxBalanceBridge expansion_scaled=$maxExpansionScaled normalization=$maxNormalization"
}
"PASS samples=$samples identity_max=$maxIdentity balance_bridge_max=$maxBalanceBridge expansion_scaled_max=$maxExpansionScaled normalization_scaled_max=$maxNormalization"
