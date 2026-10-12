$ErrorActionPreference = 'Stop'
$rng = [System.Random]::new(20261012)
$count = 10000
$maxTrace = 0.0
$maxConservation = 0.0
$maxMassU = 0.0
$maxMassV = 0.0
for ($i = 0; $i -lt $count; $i++) {
    $re = 0.5 + 4.5 * $rng.NextDouble()
    $mu = 0.001 + 0.3 * $re * $rng.NextDouble()
    $Mbar = $re - $mu / 2.0
    $r = $re * (1.1 + 4.0 * $rng.NextDouble())
    $c = 1.0 + 100.0 * $rng.NextDouble()
    $tu = -1.0 + 2.0 * $rng.NextDouble()
    $tv = -1.0 + 2.0 * $rng.NextDouble()

    $F = 1.0 - 2.0 * $Mbar / $r + $re * $re / ($r * $r)
    $Fp = 2.0 * $Mbar / ($r * $r) - 2.0 * $re * $re / ($r * $r * $r)
    $Fpp = -4.0 * $Mbar / [Math]::Pow($r, 3) + 6.0 * $re * $re / [Math]::Pow($r, 4)
    $Fppp = 12.0 * $Mbar / [Math]::Pow($r, 4) - 24.0 * $re * $re / [Math]::Pow($r, 5)

    $A = $c * (2.0 * $F * $Fpp - $Fp * $Fp) / (192.0 * [Math]::PI)
    $Aprime = $c * $F * $Fppp / (96.0 * [Math]::PI)
    $B = $c * $F * $Fpp / (96.0 * [Math]::PI)
    $Bprime = $c * ($Fp * $Fpp + $F * $Fppp) / (96.0 * [Math]::PI)
    $Tuu = $tu + $A; $Tvv = $tv + $A; $Tuv = $B

    $traceResidual = -4.0 * $Tuv / $F + $c * $Fpp / (24.0 * [Math]::PI)
    $conservationResidual = $F * ($Aprime - $Bprime) + $Fp * $B

    # Project the trace-reversed reduced source onto grad(r) in RN null coordinates.
    $trace = -4.0 * $Tuv / $F
    $T_u_u = -2.0 * $Tuv / $F; $T_u_v = -2.0 * $Tuu / $F
    $T_v_u = -2.0 * $Tvv / $F; $T_v_v = -2.0 * $Tuv / $F
    $r_u = -$F / 2.0; $r_v = $F / 2.0
    $massU = ($T_u_u - $trace) * $r_u + $T_u_v * $r_v
    $massV = $T_v_u * $r_u + ($T_v_v - $trace) * $r_v
    $massUResidual = $massU + $Tuu + $Tuv
    $massVResidual = $massV - $Tvv - $Tuv

    $speed = -0.9 + 1.8 * $rng.NextDouble()
    $du = 1.0 - $speed; $dv = 1.0 + $speed
    $mapResidual = [Math]::Max([Math]::Max(0.0, -$du), [Math]::Max(0.0, -$dv))

    $maxTrace = [Math]::Max($maxTrace, [Math]::Abs($traceResidual))
    $maxConservation = [Math]::Max($maxConservation, [Math]::Abs($conservationResidual))
    $maxMassU = [Math]::Max($maxMassU, [Math]::Abs($massUResidual))
    $maxMassV = [Math]::Max($maxMassV, [Math]::Abs($massVResidual))
    if ([Math]::Abs($traceResidual) -gt 1e-10 -or [Math]::Abs($conservationResidual) -gt 1e-10 -or [Math]::Abs($massUResidual) -gt 1e-10 -or [Math]::Abs($massVResidual) -gt 1e-10 -or $mapResidual -gt 0.0) {
        throw "RN collar identity residual at sample ${i}"
    }
}
Write-Output "PASS: $count RN samples; max trace $maxTrace; conservation $maxConservation; mass-U $maxMassU; mass-V $maxMassV."
