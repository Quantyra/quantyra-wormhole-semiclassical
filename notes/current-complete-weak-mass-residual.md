# Complete weak residual for the correlated current mass adjoint

2026-10-25 HST, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: progress (domain test and original output comparison).
The full stable semiclassical wormhole objective remains incomplete.
The closed feasibility audit and original physical targets are unchanged.

CLAIM TO VERIFY: the COMPLETE current adjoint residual satisfies
R=||j Psi_ref||E*<2.1e-5<3.27169e-5. Together with the already closed
source-pairing allowance, this certifies the mass of the FULL COUPLED
LINEAR correction at the seed. It does not certify its original graph,
the mass-column operator norm kappa_m, a nonlinear solution or stability.

## 1. Use the complete weak operator, not another pointwise estimate

The previous attempt left a pointwise outer derivative uncomputed.
For this ONE adjoint residual, the retained complete weak kernel is
sufficient. Its constant probe already contains the density contact;
its gluing is already the actual finite-domain return. Extract its
absolute bilinear estimates BEFORE the later outer spectral comparison
and IMS localization, which are not needed for the present input.

Keep Q, z_h, gamma, beta=lambda B^2/C, the paired mass/log physical
test domain, common subtraction and finite coefficients. The physical
mass lift and its .1/.2 transition, original .05/.1 graph cutoff and
regular bolt are unchanged. Use the pole-removed correlated trial
Psi_ref=(n,b) of current-reference-clock-adjoint-trial.md:

    p=beta m_T=gamma^2/(z_h D_A), D_A=z_h^2 Q''(z_h)/2+3gamma z_h,
    A=I_Q+2p J, I_Q=int_0^z 1/Q, J=int_0^z t^3/Q^2,
    n=-QA/z, b=1+Q'A/2+p z^3/Q.                         (1)

Both reference equations vanish, with their actual bolt matching.
The two separate shape profiles are NOT used as global directions.
They have opposite logarithmic bolt singularities; the sum in (1)
is regular. No zero physical fixed-clock mode is asserted.

Let alpha=2GC/(B^2 beta)=2eta_G/beta. Normalize orbit integrals by
removing their common factor 2pi N/gamma. Thus the reduced area is dz,
and the reduced unit-tangent measure has total volume 2pi z_h. Multiplying
a canonical quadratic form with coefficient C by alpha converts it to
the same h=I''/(|W|beta) used in the mass certificate. This is accounting,
not an improvement by rescaling.

## 2. Close the actual joint profile bounds, including the bolt

Write delta=z_h-z and

    Q=delta U, U=2gamma+a delta+b3 delta^2+b4 delta^3+b5 delta^4+b6 delta^5.

The actual intervals give a in (.999999,1.000001), b3 in (-1.334,-1.333),
b4 in (.5,.501), |b5|,|b6|<1e-27, 0<gamma<1e-30 and z_h>.99.
Consequently P=(U-2gamma)/delta lies between .16 and 1.001 on the real
interval, Q>=delta^2/8, and Q<.17. All are statements about the actual
polynomial, with the positive thermal term retained.

For z<=z_c=1-1e-8 compare integrals with the explicit polynomial
Qstar=(1-z)^2(1+2z+3z^2)/6. This is ONLY an inequality for the exact
primitive, not a replacement metric or quantum domain. Direct coefficient
bounds give |Q/Qstar-1|<1e-8 and |Q'-Qstar'|<1e-9(1-z). Its primitive is

    Istar=1/(1-z)-1-(4/3)log(1-z)+(2/3)log(1+2z+3z^2)
       +(7sqrt(2)/6)[atan((3z+1)/sqrt(2))-atan(1/sqrt(2))].       (2)

Since Qstar'=-2z^2(1-z), bstar=1-z^2(1-z)Istar. To bound it on the
WHOLE interval, let H=1.5/[z^2(1-z)]-Istar. The sign of H' is the sign
of 5z^3-z-2. This polynomial has one positive zero, lying in [.825,.829].
H tends to positive infinity at both ends; interval evaluation over that
root bracket gives H>3. Thus -.5<=bstar<=1. The integral comparison
Istar/(1+1e-8)<=I_Q<=Istar/(1-1e-8), Istar<=6/(1-z), and
J<=64/[3(z_h-z)^3] give |b-bstar|<1e-6 on this whole middle/outer interval.
No sampled maximum is used.

At the bolt use the exact pole-removed identity. Set c=z_h^2/(2D_A),

    V=[delta U+2p(z_h-delta)^3-cU^2]/delta^2,
    W=[-c(U+delta U_delta)U/2+p(z_h-delta)^3]/delta,
    A=c/delta+Areg, Areg_z=V/U^2,
    K=V(0)/a, Rreg=Areg-K/U.

Both divisions are exact polynomial divisions. Let N1=(V-K U_delta)/delta
and Dnum=U-K(U+delta U_delta)/2+W. Then

    b=Dnum/U-(U+delta U_delta)Rreg/2,
    Rreg_z=delta N1/U^2.                              (3)

The constant and linear coefficients of Dnum have factors gamma^2 and
gamma, respectively. The checker verifies these identities before interval
substitution. The resulting coefficient bounds on |delta|<=2e-8 give
|V|<100, |W|<100(gamma+|delta|), |N1|<1e4 and
|Dnum|<5000(gamma+|delta|)^2. Real U>=gamma+delta/2 therefore gives
|Dnum/U|<1e4(gamma+delta). The exact anchor from (2), retaining its
comparison error and the tiny J term, gives |Rreg(z_c)|<1e6.

Integrate (3), using delta/U^2<=4/(2gamma+delta), to obtain

    |Rreg(delta)|<=1e6+4e4 log[(2gamma+delta_c)/(2gamma+delta)].

Here delta_c=z_h-z_c<1e-8. Also U+delta U_delta<=3(gamma+delta).
The elementary maximum x log(A/x)<=A/e consequently gives

    |b|<=1e4(gamma+1e-8)
       +1.5[(gamma+1e-8)1e6+4e4(2gamma+1e-8)/e]<.02.       (4)

This bounds the whole arbitrarily small bolt interval without dividing
independently rounded roots. Equations (2)-(4) prove

    ||b||infinity<=B0=1.000001, ||n||infinity<100.          (5)

For the latter, on z<=.1 use .16<Q<.17 in the integral in (1). On
z>=.1 outside the cold interval, A<9/delta and
Q<=1.001 delta^2+2gamma delta. In the cold interval use
n=-U(c+delta Areg)/z and |Areg|<=1e9+400/(gamma+delta), as below.

We need deliberately loose ONE-SIDED geodesic derivative bounds, not
the earlier unknown whole-chart S145 norm. These can be proved uniformly
from the actual polynomial. In the holomorphic normal Goursat chart,

    z(w,0)=z0+sqrt(Q(z0))w/2, z(0,v)=z0+sqrt(Q(z0))v/2,
    z_wv=Q'(z)exp(2sigma)/4, sigma_wv=Q''(z)exp(2sigma)/8.

On |w|,|v|<=1/32 the integral equations preserve

    |delta-delta0|<=min((gamma+delta0)/8,.02), |sigma|<.026.  (6)

Indeed |Q'|<6(gamma+delta0), also |Q'|<.7, |Q''|<20 and |Q'''|<30
in this cone. The two displacement estimates are <.047(gamma+delta0)
and <.014, and the sigma estimate is <.003. Their scaled Picard row sums
are <.02, proving existence and analytic dependence on this disk.
This includes the regular-bolt chart. At the outer seam it is the
analytic continuation of the ONE-SIDED coefficients used solely to
bound their jets. The physical even double remains C2,1 throughout.

For delta0>=delta_c, (6) gives |Q|>delta0^2/32 and |Q'|<6delta0.
Continue A from its real value <9/delta0 along the straight segment;
its increment is <5/delta0. Hence |b|<50 on that complex chart.
For delta0<=delta_c, the pole-removed polynomials instead give
|U|>.8(gamma+delta0), |Areg|<1e9+420/(gamma+delta0), and |W/U|<142.
Thus |b|<800. Retain the generous uniform complex bound 1e6.

Cauchy's formula gives straight-ray jets through order four bounded by
2^k k! 32^k 1e6. At a normal center gamma_geod''=0; the retained seed
chart bounds give |gamma_geod'''|<=15 and |gamma_geod''''|<=4e6.
The chain rule then proves, on each physical side and at the regular bolt,

    |Xb|<=K1=1e8, |X^2b|<=K2=1e11,
    |X^3b|<=K3=1e14, |X^4b|<=K4=1e16.                (7)

No analytic assumption is made about general mass/log test fields.
These bounds apply only to the explicitly analytic reference trial.

## 3. Complete principal form with its actual seam and all angular modes

Use Qgeo and P=iX on the actual complete geodesic flow of the paired
C2,1 double, exactly as in stationary-curved-quantum-tt.md. In reduced
unit-tangent measure its spectral multiplier obeys

    0<=q_geo(k)<=2C(k^2+|k|^(17/4)).                  (8)

Here is an all-mode proof, avoiding an unbounded comparison of Qgeo
with the physical energy near the seam. The retained Bessel moments
and 0<=g(s)<=min(s^6/720,s^4/24) give, per massive level,

    P_M(k)<=min(k^6/(525pi^2 M^4), k^4/(90pi^2 M^2)).

Use min(x,y)<=x^(1/8)y^(7/8) and the complete counting bound
sum d/M^(9/4)<=36C 18^(-1/4). Its coefficient is <.02C.
The common fourth restoration has coefficient |F_eff|C/(90pi^2)<.113C.
Since k^4<=k^2+|k|^(17/4), (8) follows. No angular, thermal, geodesic-time
or finite-domain truncation is introduced. The Bessel integral theorem
is used on positive masses; the paired domain and global flow are supplied
by their separate ordinary arguments.

The trial is C2 across the actual seam and piecewise C4. Let f=X^3b.
It has |f|<=K3, one-sided |Xf|<=K4 and jump at most 2K3. In reduced
Liouville measure the set of trajectories crossing the seam in time
t<=t0=1/32 has measure <32t. This follows by integrating the normal
flux: the normalized seam length is sqrt(beta), and
int_0^(2pi)|cos theta|dtheta=4; 32 is a conservative reserve.
There is at most one crossing in this time. Between two crossings the
normal coordinate would satisfy Y''=c(t)Y with |c|<=2, since
|a'/a|<=2|Y| on |Y|<=1/32. Poincare on the interval between its two
zeros would require t>=pi/sqrt(2), a contradiction.

Thus

    ||(exp(-itP)-1)f||2^2<=4pi K4^2 t^2+256K3^2 t, t<=t0,
    ||(exp(-itP)-1)f||2^2<=8pi K3^2, t>=t0.

Spectral calculus and Tonelli give the fractional-difference identity.
Its coefficient for power 1/4 is <1 (integrating
1-cos t>=2t^2/pi^2 on [0,pi] already proves this). Therefore

    ||P^(13/4)b||2^2
      <=(8pi/3)K4^2 t0^(3/2)+512K3^2 sqrt(t0)
                                      +16pi K3^2/sqrt(t0)=Bfrac^2. (9)

The continuity through second derivatives ensures b is in Dom(P^3);
the fundamental theorem along the flow proves this successively.
Equation (9) then supplies the extra quarter derivative. It does not
pretend the cubic seam is H4 or smooth.

## 4. Physical test bounds and the entire massive residual

Keep the fixed-seed physical energy and output comparison of the previous
attempt: ||T/z^2||2<17.6 E and ||nu/z^2||2<12.65 E. The retained TT
gradient bound C J1/E0<=2e5 eta_G and l>=.99 F_l E0 also give

    ||XT||2<=T1 E, T1^2=pi*1e5 beta/(.99F_l).          (10)

All X norms in (7)-(10) are in reduced unit-tangent measure. The retained
joint positive reserve l>=.99F_l E0+C S2/(100pi), together with
J2<=2S2+2*500^4 J0,outer, gives

    sqrt(alpha)||X^2T||2<39 E.                        (11)

For an explicit substitution use j0=1000eta_G/beta, jo=(16*1e-8/25)j0,
j1=2e5eta_G. The square of the left constant is at most
150pi^2+[(3pi/2)500^4 jo+(pi/2)j1]/(.99F_l)<39^2.
This uses the joint positive reserve; it does not bound Qgeo[T] by E.

Now extract the COMPLETE actual massive kernel decomposition from the
ordinary curved-form proof, sections 3-6, with its current-seed transfer:

* Principal part: Qgeo with exactly its common fourth restoration.
* Row-subtracted actual gluing: L2 operator norm <=6e24 C. Its complete
  constant row is in the next item, not counted twice or omitted.
* Actual constant contact: |Theta_p0|<=1e20q^2(1+1/Y), including tadpole,
  bubble, density variation, common finite parts and the actual-to-cap
  mass derivative. It is not a derivative of a source-value envelope.
* Remaining gain-two/contact form: for epsilon=1e-24 its absolute
  diagonal bound is (epsilon C/3)||X^2T||2^2+1e15C||XT||2^2.
  The second contact and curvature contact are included before this bound.
  Polarization gives the PRODUCT of the combined norms on the H2 core.
  It does not give the sum of the two matching derivative products.

Account explicitly for the full geodesic tail outside the local pair
radius 1/128. After its already fixed polynomial contacts are retained,
its row-subtracted kernel is bounded by the second tail moment
|J2,tail|<=36C*128^2. Unitarity and
||exp(-irP)T-T||<=r||PT|| put this contribution inside an additional
1e7 C||XT||2^2 reserve, including angular factors. Its second contact
was fixed by the whole-flow quadratic probe; its fourth contact is the
same common restoration, not an extra free coefficient. The finite-
difference gradient coefficient plus 1e13 for the second contact and
this 1e7 tail still totals less than 1e15. Thus no long geodesic return
is silently replaced by a compact cap or omitted from (12).

Specifically use the retained finite-difference estimate with
K=1e13, R0=1/400 and split radius bsplit=sqrt(epsilon/(4K)). Its
gradient coefficient K[1/4+2log(R0/bsplit)]+1e13+1e7 remains below
1e15 even for epsilon=1e-24. This is a change in the proved physical
estimate, with the same norms and contacts. For epsilon=1e-6 the valid
combined-norm product is too large; the earlier separate-product draft
was invalid. Both passing draft arithmetic snapshots are preserved and
superseded, not accepted as proof. The corrected product is evaluated
below before acceptance. No sign or cross-derivative cancellation is
assumed.

Only this last form is polarized on b. The full p0 form need not have
a finite diagonal value on a root-zero field. Its MIXED pairing with
physical tests is estimated separately using (5) and Y>=sqrt(5)z.
The seam Delta R contact pairs b(0)T(0)=0; integrations producing a
first derivative also use the physical T'(0)=0. The log phase obeys
these conditions. No fourth derivative of the reflected metric is used.

Consequently the FULL canonical massive Hessian acting on b has E* norm
bounded by the sum

    Rgeo =2alpha [sqrt(2pi)K1+Bfrac] T1,
    Rgain=alpha sqrt[2pi((epsilon/3)K2^2+1e15 K1^2)]
                  *sqrt[(epsilon/3)*39^2/alpha+1e15 T1^2],
    epsilon=1e-24,
    Rglue=17.6 alpha*6e24 B0,
    Rp0=35.2 L_q*1e20 B0 sqrt(1/5+1/(2sqrt(5))+1/15).   (12)

The principal mixed bound follows from (8), spectral Cauchy--Schwarz
with derivative orders (1,1) and (13/4,1), and (9)-(10). The gluing
uses ||b||L2(dz)<=B0. The p0 square-root is the EXACT integral of
[z^2+z/sqrt(5)]^2 over [0,1], enlarged from [0,z_h]. Here
L_q=Gq^2/(lambda B^4) is unchanged. No separate density cost is added:
the actual constant contact already contains it.

For the gain-two form, the diagonal estimate makes it a contraction
relative to G(T)^2=(epsilon/3)||X^2T||2^2+1e15||XT||2^2 on the
gradient seminorm completion (constants lie in its kernel).
Polarization and the self-adjoint Riesz representative therefore give
|B(b,T)|<=C G(b)G(T), exactly the product used in (12). The fixed
gradient contacts and tail are included in the same G bound.

This also addresses the needed continuum passage. On compact physical
tests the restored actual kernel/constant-contact identity is the original
common-regulator Hessian. Its all-mode bounds above act on the JOINT
bounded, regular trial and dominate every displayed mixed term. Extend
in the physical E norm using (9)-(12). This equals the retained auxiliary
q_P-hv_ref by uniqueness of that extension. It is NOT a C2 theorem for
an unrestricted open-clock nonlinear family. The common regulator is
removed on compact tests before the auxiliary outer cutoff, as before.
No joint rate, smooth surrogate or finite selection of modes is assumed.

## 5. Both complete local rows, with the boundary representative retained

Use the actual affine canonical action generator and compact tests.
Set phi=-log z+nu, shape u, Delta=(Q partial_z)', and

    L=-Q''-2Delta u, M=Delta phi, N=Q phi'^2, D=L-6M-6N,
    Wdim=exp(-2phi)[2+C exp(-2u)D]+12C beta,
    alpha_R=8pi c_R g/(C^2 beta), H=g/(2pi C^2 beta), g=G lambda.

The ENTIRE local difference from F_l I0, excluding the LLL evaluated
separately below, has density

    alpha_R exp(4phi+2u) Wdim^2
      +[k_vac exp(4phi)-F_l epsilon(4exp(-phi)-3exp(-2phi))+a_area]exp(2u)
      +H{C^2 exp(-2u)[phi L^2/60+LM/30-M^2/10+LN/30-7MN/90-7N^2/180]
                  -C phi L/18+(1/15-q^2/3)phi exp(2u)-CN/18},       (13)

where k_vac=11g beta/(120pi)=H*11C^2 beta^2/60. The area coefficient is
F_l Delta Q/(CB^2), so |a_area|<3.4e-22 beta<6e-23 by the current
reference-potential bound. The exact relation C u_seed=1 is imposed
before estimating; otherwise the curvature-square cancellation is lost.
No product-gradient subtraction, finite coupling or vacuum is retuned.

Compute both Euler densities f_nu,f_b of (13) at the actual seed. They
are finite Laurent/log polynomials. Their compact-test variation along
the joint trial has a particularly useful exact form. The radial
generator is xi^z=QA=-nz, with (xi^z)'+1=2b; the reference family also
changes Q by 2p z^3. If a3 denotes the independent coefficient in
Q=beta+u_seed z^2-a3 z^3+b4 z^4-2epsilon z^5+epsilon z^6, covariance gives

    H_local Psi_ref: f_i ->2b f_i-nz f_i'-2p partial_(a3) f_i,
    r_i=-(2b f_i-nz f_i'-2p partial_(a3) f_i)/beta.       (14)

This is a LOCAL compact-test Ward identity, not separate global mass
or clock domain variations. The time-density and off-shell terms are
retained. The seed is not declared stationary for the full action.

The code differentiates (13) by the Euler formula through order two,
imposes C u_seed=1 and the exact k_vac relation, then estimates the
resulting whole polynomials. For a monomial z^j log^k z the exact
weighted L2 norm is sqrt((2k)!/(2j+5)^(2k+1)), j>-5/2.
All combined terms satisfy this condition. With (5), (14) and (6) of
the preceding output comparison, the ENTIRE two-row local residual
is <1e-18. Its computed norm ledger is retained, including the
uncertain a_area interval rather than choosing its midpoint.

The reduced Weyl density differs from the raw one by the known boundary
flux. That flux need not vanish on an unrestricted open-clock action
variation. Here (14) is derived on compact physical tests; its weighted
L2 bound proves a unique E* extension. The original boundary completion
and normal mass flux remain in h and K_m. A noncontinuous mass-boundary
functional cannot be appended to this extension. Thus no additional
uncomputed boundary term is discarded or counted twice. This argument
does not redefine the physical boundary charge or action domain.

The exact fixed-domain LLL Hessian is -q/(12pi)int grad b.grad T.
Using (7),(10) gives its normalized cost

    R_LLL <=alpha(q/C)K1 T1/(12pi sqrt(pi))<1e-17.       (15)

It retains the prepared period and boundary realization of the accepted
LLL formula; no separately adjustable thermal state is introduced.

## 6. Execute the whole residual and state exactly what it closes

Freeze the ordinary proof, checker and hypothesis review before execution.
Substitute all current Arb seed intervals in (2)-(15), including the
latest source-pairing bound. The COMPLETE residual ledger is

    R <=Rgeo+Rgain+Rglue+Rp0+1e-18+1e-17 <2.1e-5.      (16)

The expected dominant contributions are about 7.42501e-6 (actual gluing),
9.19378e-6 (constant contact) and 4.05326e-6 (gain-two gradient bound).
All three are upper bounds, not observed signed values. No cancellation
between unknown signed quantities is assumed. The calculation preserves
all angular/thermal modes, coincidence contacts and actual returns.

Apply the unchanged full coupled weak inverse ||h^-1||E*->E<=1,
||F||E*<1.1e-10 and the complete normal flux K_m>1.9999999764:

    |m_linear| <=[I_new+1.1e-10 R]/1.9999999764,
    I_new<9.811394479e-16.                              (17)

If (16) passes, (17) closes the selected actual LINEAR mass target
2.29e-15. This is the mass of the complete two-row linear correction,
including its scalar/quartic, correlated slow, fast and bolt couplings
through h. It is not m_T, an independent mass-column allowance or a
nonlinear solution's mass. No finite trial space replaces that inverse.

The original norm conversion remains 2133<||d_chi||H3<2174. Compute
the sufficient remaining lifted allowance 1e-11-2174|m_linear| from
(17), but do not assert that the actual lifted graph meets it. In
particular kappa_m<=5e-5 remains UNKNOWN: one affine-source mass
observable does not control an independent phase column.

The full graph still must retain the scalar quartic row 20beta Fstar c4
and the exact slow pair z^4-z^(4+e_star), e_star about 1.865097e-34.
The conditional kappa_m costs and residual budgets remain

    affine 2.14765e-12, direct graph .0024731,
    Y_other<=5.34234999994958e-12,
    eta_other+(.0195239+4.583e-13)gamma_E<.0012069.

No other entry is set to zero. The .135973<.1599164 screen keeps only
its two-fast-input family scope. Original physical-neighborhood input
comparisons, nonlinear Y0/Z0/Z2, exact compatibility, matching,
validity/noise and all-sector stability are unproved. The preserved-ball
and contraction conditions stay Y0+1e-8 Z0+5e-17 Z2<=1e-8 and
Z0+1e-8 Z2<1. This is not a controlled stationary wormhole theorem.

RESEARCH DECISION: if (16)-(17) pass, retain the construction and the
complete weak adjoint method for this observable. Stop optimizing its
source pairing or residual. The mass conversion obstacle is then removed
for the linear correction; the ORIGINAL coupled graph is the decisive
remaining quantitative task. If a supporting bound fails, retain that
failure and do not publish the conditional mass conclusion as established.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED HERE): certify the original
coupled fast graph of the SAME full linear correction, with scalar quartic,
correlated slow pair, physical lift transition and finite-domain returns;
test it against the remaining allowance from (17), or directly against
the original 1e-11 graph target. Do not substitute the fast-input quantum
screen or infer a neighborhood estimate from this single-source result.

Remaining to-do list: original coupled graph and independent phase columns;
original physical-norm/neighborhood comparisons; nonlinear admission and
contraction; exact compatibility; physical matching; semiclassical validity
and noise; all-sector stability; charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
