# Full trace component of the finite-seed Newton defect

2026-10-22, E003/S169. Ordinary argument before frozen verification.
This evaluates a component of the selected full-defect calculation;
the signed coupled bolt pairing remains uncomputed. The closed feasibility
audit and the full stable semiclassical wormhole objective are unchanged.

## 1. The actual defect, before estimating it

Use the current seed, unchanged physical action and first reference
correction (f,a,s_R) of fixed-clock-seed-bolt-feasibility.md. Its reference
constraints use G, not the effective Einstein prefactor used to prove the
full Hessian estimate. Neither split changes the physical action.
Write R4 for scalar curvature, R for the bolt radius, D=D_R, and

    T=R4+12lambda0+epsilon[7C4^2/12+11J^2/6-20I-11R4^2/72]-D Delta R4,
    epsilon=G/(60pi), T0=T[seed], I=q^2/(2r^4),
    J^a_b=diag(-tau,-tau,tau,tau), c=F'/r+R4/6.

Thus C4^2=12c^2. All variations below first keep the areal coordinate r
fixed. Let Q=G(rho-p), for the COMPLETE evaluated joint seed source.
The reference Einstein constraints and angular Ward equation imply

    delta R4=-T0,
    delta Ric^t_t=-2G rho/r^2-T0/2,
    delta Ric^r_r= 2G p/r^2-T0/2,
    delta Ric^theta_theta=Q/r^2.

The sphere sectional curvature is (1-F)/r^2. Therefore, exactly,

    delta J^2=2tau T0+8tau Q/r^2,
    delta C4^2=-24c[(f+Q)/r^2+T0/6], delta I=0,
    delta(Delta)R4=f(R4''+R4'/r)-Q R4'/r.                 (1)

In particular the nonzero R4' terms cannot be dropped on this seed.
Substituting into T0+delta T gives the ENTIRE areal trace defect

    e_T=D Delta T0+epsilon K_T T0+A_f f/r^2+A_Q Q/r^2,    (2)
    K_T=11R4/36-7c/3+11tau/3,
    A_f=-14epsilon c-D(r^2 R4''+r R4'),
    A_Q=epsilon(44tau/3-14c)+D r R4'.

No traceless quantum variation has been inferred from the trace anomaly.
The anomaly determines the full trace derivative; the evaluated rho,p
still enter (2), and the independent shape variation remains nonlocal.
At infinity (2) recovers the earlier nonzero r^-4 coefficient defect.
Here (2) applies on the whole finite exterior, including the regular bolt.

## 2. Quantify the current reference correction, without an old transfer

The earlier current-source test transferred the quantum lapse bound and
Theta' bound K=1e15, but NOT the complete energy derivatives. Check those
derivatives now, including the large local coefficient. Put

    Q_rho=(G/lambda) sup |rho'|/r, Q_p=(G/lambda) sup |p'|/r.

The same current rational/paired-domain hypotheses used in that test
give the quantum vacuum-subtracted derivative bounds
(31K+1e6)Gq^2/(lambda B^4) and (29K+1e6)Gq^2/(lambda B^4).
These are the explicit Theta'/Ward formulas (27)-(31) of
stationary-centered-exterior-source.md. The local Weyl/LLL subtotal fits
inside 1e6 without the old reference term. The previously checked fixed
mu0 restoration is retained. No second derivative of Theta is assumed.

Here is the separate current local accounting. Set x=lambda B^2,
t=r/B, z=1/t, u=1/x, beta=delta/x, and write

    F/x=z^-2+u-(4+2u+beta)z+(3+u+du)z^2-2du z^3+du z^4,
    Rbar=B^2 R4=-12x+d(4z^5-6z^6).

For H_ab=2R4 Ric_ab-g_ab R4^2/2+2(g_ab Delta-Hess_ab)R4,
the added curvature-square source is rho_R2=8pi c_R r^2 H^t_t,
p_R2=-8pi c_R r^2 H^r_r. In dimensionless form

    B^4 H^t_t=-2Rbar taubar
                 +2[F Rbar_tt+(F_t/2+2F/t)Rbar_t],
    B^4 H^r_r=-2Rbar taubar+2(F_t/2+2F/t)Rbar_t.         (3)

The absolute polynomial coefficient sums for
partial_t(t^2 B^4 H^a_a)/(t x^2) are below 200 for both rows.
Use u<=1/(8e9), |d|<=8.6e-19 and 0<=beta<=7h_max^2,
h_max=2000^-10. The resulting derivative contribution is
<200(8pi)c_R G lambda.

The original-G reference density is explicitly

    rho_ref=[3(lambda0-lambda)r^2+(Q_seed-Q_phys)/r^2
                         -4d B^3/r^3+3d B^4/r^4]/(2G),
    p_ref=-rho_ref,
    Q_seed=B^2(1+3x+d), Q_phys=piG q^2/g_R^2=F_l Z.

Its vacuum cancels the quantum vacuum exactly. The remaining derivative
bound is at most

    |Delta Q|/(lambda B^4)+(1-F_l) Z/(lambda B^4)+12|d|/x
       <3.4e-22+1.179e-8(3.001)+12(8.6e-19)/(8e9).

Using Gq^2/(lambda B^4)<1.001/2000^8 and G lambda<2000^-10,
the sum of these three contributions proves

    Q_rho<1e-7, Q_p<1e-7.                              (4)

This deliberately loose bound retains the local pieces separately.
Their cancellation could sharpen (4), but is unnecessary for this trace
test. It does not replace the sharper already proved bolt value or lapse.

Use |a|<3.3e-10+9.182e-37, |s_R|/(Bh)<3.33e-10,
|G rho(R)/(kappa R)|<2e-12, and the unchanged rational inequalities
1.9lambda<=F''<=12.1lambda, 5<kappa/(lambda Bh)<7.
The exact reference reconstruction and subtraction at the moving bolt
in stationary-centered-first-order-geometry.md, applied to (4), give

    fbar=f+s_R F', |fbar/F|<3.004e-7,
    H=a+fbar/F, |H|<3.01e-7,
    xi=s_R-F integral_r^infinity H/F,
    |xi|/(r-B)<3.02e-6, |xi'|<8e-5.                    (5)

For specificity, the fbar bound is |a(R)|+3Q_rho+80|s_R|/B.
The displacement uses (F/(r-B)) integral_r^infinity 1/F<10.
On y=r-R<=Bh, the subtraction in the same reconstruction bounds
|H|<1e-7 y/(Bh): the terms are at most (Q_rho+Q_rho+Q_p)/10,
.7|a(R)| and a term <1e-25, in units y/(Bh).
Then |J|<[1e-7+2.8 sup|H|]/(2kappa), and F'/(2kappa)<2.21 there.
On y>=Bh use F' integral_r^infinity 1/F<260.
The identity xi'=H-F'J proves the stated global derivative bound.
Only a value bound on a, not a claimed total-variation bound, is used.

## 3. Full radial trace and transport bounds

Keep the earlier physical weight, with no renormalization of the test:

    w(t)=((t-1+h)/t)^(2/5), ||v||_w=sup |v|/w.

The product trace relation in product-centered-scalar-cold-response.md
gives an explicit convenient polynomial representation. Put
S=epsilon q^2/(lambda B^4), eta=epsilon lambda, d_l=D lambda,
r(z)=4z^5-6z^6, and L p=z^4 partial_z[(F/x)partial_z p]. Then

    T0/lambda=S A+eta Bp-d_l du Lr,
    A=-10z^4-20z^5+30z^6,
    Bp=P/x^2+(P_*/x^2)r/2,
    P/x^2=7(c/lambda)^2+(22/3)(tau/lambda)^2
                                      -11[(R4/lambda)^2-144]/72,
    P_*/x^2=2(6+du)^2+(22/3)(6+du)u+(22/3)u^2+22.       (6)

Every zero-temperature Bp and Lp used below vanishes at z=1.
For p=sum c_n z^n with sum c_n=0, choose an integer j and write
p=sum c_n(z^n-z^j); then ||p||_w<=sum |n-j||c_n|.
Choose the smallest certified sum over 0<=j<=degree(p); this preserves
the exact zero before taking absolute values. For A the choice is j=5,
and for its leading Laplacian it is j=9.
For (t-1)partial_t p use sum n|c_n|. Expand the thermal terms
in beta. Since w>=h_min, beta/w<=7h_max^2/h_min,
h_min=h_max/8, their coefficient sums are also rational upper bounds.
Thus this is a whole-exterior coefficient bound, not a radial sample.
The frozen checker evaluates these sums, proving

    ||T0||_w/lambda<9e-28,
    ||(t-1)partial_t T0||_w/lambda<7e-27,
    ||Delta T0||_w/lambda^2<4e-25.                     (7)

The elementary F/x polynomial also gives F/(lambda r^2)<=7(1-z)^2.
Equation (5) and the horizon source with (4) imply

    ||f/r^2||_w/lambda<2.11e-6,
    ||Q/r^2||_w/lambda<2.01e-7.                       (8)

For the latter, integrate |rho'|,|p'|<=lambda r/(G) times
their Q bounds from R, using p(R)=-rho(R). The result is bounded
by 1e-7(1-R^2/r^2)+2.8e-11 h R^2/r^2 before division by w.
The horizon term in (8) is negligible even with w>=h^(2/5).

Use |c|<35lambda, |tau|<3.01lambda and |R4|<12.001lambda.
Then |K_T|<100lambda. The exact R4 polynomial gives

    |A_f|<=490eta+316d_l |d|/x,
    |A_Q|<535eta+56d_l |d|/x.

Equations (2),(7),(8) supply a numerical bound on the areal trace
defect without evaluating a quantum Hessian kernel.

## 4. Moving density and bolt pullback are part of the defect

The physical solver uses a fixed conformal disk. In that chart the
areal pullback displacement is xi in (5). Varying sqrt(g)T before
dividing by the SEED volume gives the exact density-correct trace defect

    e_disk=e_T+a T0+(1/r^2)partial_r(r^2 xi T0)
          =e_T+(a+xi'+2xi/r)T0+xi T0'.                (9)

In particular the seed is off shell: neither extra term is zero.
Equation (5) gives |a+xi'+2xi/r|<8.7e-5. Equations (7) and
|xi|/(r-B)<3.02e-6 then bound the complete disk trace row by

    ||e_disk||_w/lambda<1.1e-31.                      (10)

This row is in the physical fixed-disk density equations of
recentered-stationary-constraints.md. Any subsequent York coordinate
change must carry its off-shell contacts too. Equation (10) is NOT
the full residual's action-dual norm, a bound on the shape row, or a
quantitative identification of the full weighted inverse norm.

The already proved finite-seed conformal row estimate is

    ||A^-1 F||_cold,01 < (4/3)e18 ||F||_w/lambda.

Consequently the trace defect with ZERO independent shape-error input
has cold value/first-disk-derivative response <1.467e-13. This closes
that scalar-only defect allocation. It does not close the coupled
bolt-difference target <1.667e-9, nor raw C5 or total Y0<=5e-9.
If e_xi is the unknown shape error, the actual scalar error equation
also contains B e_xi. A sufficient scalar-row target is

    ||B e_xi||_w/lambda<=1.25014e-27,
    (4/3)e18(1.1e-31+1.25014e-27)=1.667e-9.           (11)

The bound in (11) is UNKNOWN, not a failing measured upper bound.
It concerns two scalar cold components, not the physical radius
functional, which also contains the independent sphere/shape response.

## 5. Decision and limits

The executed full trace calculation is favorable and removes this
explicit local trace defect as the dominant loss. No additional trace
source optimization is warranted now. The decisive missing quantity
is the signed SHAPE contribution and its coupled bolt return through
the actual finite-domain Hessian, with a certified weighted remainder.
Its evaluation, including the tadpole, bubble, common contacts and
actual state/domain, remains UNEXECUTED. Anomaly identities, form
coercivity, a product multiplier or supporting checker passes cannot
provide its value. Keep the construction and the chosen direct defect
method, with the same <1.667e-9 cold and <=5e-9 total targets.

The previous inverse realization and physical linear phase remain in
force, within their stated sector and hypotheses. This calculation
does not establish a quantum-source C2 map on a preserved neighborhood,
a nonlinear solution, physical matching, semiclassical validity or
all-sector retarded stability. No S169 task checkbox is completed.

Remaining to-do list: full shape-defect/adjoint-bolt pairing and weighted
error; preserved nonlinear solution, matching, validity, all-sector
stability and the remaining full stable-wormhole requirements.
