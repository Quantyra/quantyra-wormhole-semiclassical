$ErrorActionPreference = 'Stop'
$maxScaledResidual = 0.0
$cases = 0
foreach ($rho in @(-2.1,-0.6,0.0,0.9,2.7)) {
    $f=1.0+$rho*$rho
    foreach ($jet in @(
        @{pt=0.7;pr=-0.4;ptt=0.3;ptr=-0.2;prr=0.5},
        @{pt=-0.1;pr=0.8;ptt=-0.6;ptr=0.9;prr=-0.3},
        @{pt=1.1;pr=0.2;ptt=-0.4;ptr=-0.7;prr=0.6}
    )) {
        [double]$eps=1e-5; [double]$X=1.0+$eps*0.23
        [double]$Xt=$eps*$jet.pt; [double]$Xr=$eps*$jet.pr
        [double]$Xtt=$eps*$jet.ptt; [double]$Xtr=$eps*$jet.ptr; [double]$Xrr=$eps*$jet.prr
        # Background Christoffels for ds^2=-(1+rho^2)dt^2+drho^2/(1+rho^2).
        $Htt=$Xtt-$f*$rho*$Xr
        $Htr=$Xtr-$rho/$f*$Xr
        $Hrr=$Xrr+$rho/$f*$Xr
        $box=-$Htt/$f+$f*$Hrr
        $grad2=-$Xt*$Xt/$f+$f*$Xr*$Xr
        $V=2.0-2.0/$X
        $exact0 = -$f*$box-$Htt+($Xt*$Xt-(-$f)*$grad2/2.0)/(2.0*$X)-(-$f)*$V/2.0
        $exact1 = -$Htr+($Xt*$Xr-0.0*$grad2/2.0)/(2.0*$X)
        $exact2 = $box/$f-$Hrr+($Xr*$Xr-(1.0/$f)*$grad2/2.0)/(2.0*$X)-(1.0/$f)*$V/2.0
        # Linearized equation E_ab/eps = g_ab Box(phi)-Hess_ab(phi)-g_ab phi.
        $phitt=$jet.ptt-$f*$rho*$jet.pr
        $phitr=$jet.ptr-$rho/$f*$jet.pr
        $phirr=$jet.prr+$rho/$f*$jet.pr
        $phibox=-$phitt/$f+$f*$phirr
        $linear0 = -$f*$phibox-$phitt+$f*0.23
        $linear1 = -$phitr
        $linear2 = $phibox/$f-$phirr-0.23/$f
        foreach ($scaled in @(
            [Math]::Abs(($exact0-$eps*$linear0)/($eps*$eps)),
            [Math]::Abs(($exact1-$eps*$linear1)/($eps*$eps)),
            [Math]::Abs(($exact2-$eps*$linear2)/($eps*$eps))
        )) { $maxScaledResidual=[Math]::Max($maxScaledResidual,$scaled); $cases++ }
    }
}
if ($maxScaledResidual -gt 20.0) { throw "Linearization residual too large: $maxScaledResidual" }
[pscustomobject]@{status='passed';components_checked=$cases;max_residual_over_epsilon_squared=$maxScaledResidual;eps=1e-5;grid='5 rho x 3 arbitrary jets x 3 metric components'}
