$ErrorActionPreference = 'Stop'
$checks=0
$maxResidual=0.0
foreach ($q in @(1.0,7.0,31.0,100.0)) {
    $k=$q/(8.0*[Math]::PI)
    foreach ($rho in @(-9.0,-3.0,-1.0,-0.1,0.0,0.2,1.4,4.0,11.0)) {
        $f=1.0+$rho*$rho
        # Constant conformal-coordinate state difference, pulled back as a covariant tensor.
        $Ttt=-$k
        $Txx=-$k
        $Trr=-$k/($f*$f)
        $trace=(-1.0/$f)*$Ttt+$f*$Trr
        $mixedTau=(-1.0/$f)*$Ttt
        $mixedRho=$f*$Trr
        # For a static diagonal tensor: nabla_a T^a_rho = d_rho T^rho_rho
        # + Gamma^a_(a rho) T^rho_rho - Gamma^c_(a rho) T^a_c.
        $dMixedRho=2.0*$k*$rho/($f*$f)
        $gammaTrace=0.0 # Gamma^tau_(tau rho)+Gamma^rho_(rho rho)=rho/f-rho/f
        $connection=($rho/$f)*$mixedTau+(-$rho/$f)*$mixedRho
        $divRho=$dMixedRho+$gammaTrace*$mixedRho-$connection
        # Coordinate transformation: d x/d rho=1/f; covariant T_rhorho=T_xx/f^2.
        $pullbackResidual=$Trr-($Txx/($f*$f))
        foreach ($residual in @($trace,$divRho,$pullbackResidual)) {
            $maxResidual=[Math]::Max($maxResidual,[Math]::Abs($residual))
            if ([Math]::Abs($residual) -gt 1e-12) { throw "Ward/pullback residual $residual at q=$q rho=$rho" }
            $checks++
        }
    }
}
[pscustomobject]@{status='passed';checks=$checks;max_absolute_residual=$maxResidual;grid='4 fluxes x 9 radii x trace/divergence/pullback';scope='fixed-throat state difference only'}
