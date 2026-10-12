$ErrorActionPreference = 'Stop'

function Assert-Near([string]$Label, [double]$Value, [double]$Tolerance = 1e-9) {
    if ([Math]::Abs($Value) -gt $Tolerance) {
        throw "$Label residual $Value exceeds $Tolerance"
    }
}

$checks = 0
$maxResidual = 0.0
foreach ($rho in @(-3.0, -1.25, -0.2, 0.4, 1.7, 4.0)) {
    $f = 1.0 + $rho * $rho
    foreach ($tau in @(0.0, 0.37, 1.2, 2.8)) {
        # Mode rho: (Phi, Phi_tau, Phi_rho, Phi_tautau, Phi_taurho, Phi_rhorho)
        $modes = @(
            @{ n='rho'; p=$rho; pt=0.0; pr=1.0; ptt=0.0; ptr=0.0; prr=0.0 },
            @{ n='cos'; p=[Math]::Sqrt($f)*[Math]::Cos($tau); pt=-[Math]::Sqrt($f)*[Math]::Sin($tau); pr=$rho/[Math]::Sqrt($f)*[Math]::Cos($tau); ptt=-[Math]::Sqrt($f)*[Math]::Cos($tau); ptr=-$rho/[Math]::Sqrt($f)*[Math]::Sin($tau); prr=[Math]::Cos($tau)/[Math]::Pow($f,1.5) },
            @{ n='sin'; p=[Math]::Sqrt($f)*[Math]::Sin($tau); pt=[Math]::Sqrt($f)*[Math]::Cos($tau); pr=$rho/[Math]::Sqrt($f)*[Math]::Sin($tau); ptt=-[Math]::Sqrt($f)*[Math]::Sin($tau); ptr=$rho/[Math]::Sqrt($f)*[Math]::Cos($tau); prr=[Math]::Sin($tau)/[Math]::Pow($f,1.5) }
        )
        foreach ($m in $modes) {
            $e1 = -$m.ptt/$f + $rho*$m.pr - $m.p
            $e2 = -$m.ptr + $rho*$m.pt/$f
            $e3 = -$f*$m.prr - $rho*$m.pr + $m.p
            foreach ($v in @($e1,$e2,$e3)) {
                $maxResidual = [Math]::Max($maxResidual,[Math]::Abs($v))
                Assert-Near "$($m.n) homogeneous mode" $v
                $checks++
            }
        }
    }
}

# MMP static solution Phi=alpha(1+rho atan rho), checking rho Phi'-Phi=-alpha/f.
$alpha = 0.73
foreach ($rho in @(-12.0,-4.0,-1.0,-0.1,0.0,0.25,2.0,8.0,20.0)) {
    $f = 1.0 + $rho*$rho
    $phi = $alpha*(1.0 + $rho*[Math]::Atan($rho))
    $phiPrime = $alpha*([Math]::Atan($rho) + $rho/$f)
    $residual = $rho*$phiPrime - $phi + $alpha/$f
    $maxResidual = [Math]::Max($maxResidual,[Math]::Abs($residual))
    Assert-Near 'MMP static radial identity' $residual 1e-8
    $checks++
}

[pscustomobject]@{ status='passed'; checks=$checks; max_residual=$maxResidual; grid='6 rho x 4 tau x 3 modes x 3 components; 9 static radii' }
