# Quantitative higher York maps with the physical mass column

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary proof
and verification contract before execution. **Linear seed theorem;
no nonlinear candidate regularity or physical wormhole conclusion.**

The genuine York extraction and its fixed return now have numerical
bounds in the global fourth/log norms. The mass comparison is made
against an explicit scalar cubic, including its induced fast part.
The radial fourth-power return tail is retained throughout.

## 1. Fixed scalar norms and statement

Use the analytic seed, physical pairing, clock, temperature, state
and domain of [the global disk bound](collar-patch-global-disk-bound.md).
Keep its epsilon=10^-160, r=epsilon/32, L=4r, t=10^-170,
delta=7/2, chi=c_r, charts chi_i, A0 and gamma0. The letter r
without a subscript in these scales is a cutoff radius. Write R_phys
for the physical areal radius, u=B/R_phys, k=sqrt(lambda B^2/C).

For a scalar F define its analysis map

    JF=(Y^-delta chi F, (A0^-1/2 chi_i F)_i).

Its ambient spaces are H^s(R) direct-summed with the chart H^s(R2)
spaces. Let A_s(F)=||JF|| in that Hilbert sum, and A_4log use
the same multiplier (1+|xi|^2)^2 log(2+|xi|) in every summand.
At the bolt these are the retained rotational disk norms.

Let (tau,c) be the genuine York extraction of (S,T); c is the
sphere component of h_TT, and Rcal(tau,c)=(tau+S_c,T_c) is the
exact inverse return. The preceding disk mass coefficient is
m=m_Y=-3k^3 m_c. Choose the scalar mass representative

    c_mass(Y)=c_L(Y)u(Y)^3,
    c_f=c-m_c c_mass.

No new tail or independent mass is selected. Define

    Ynorm(tau,c)^2=A_4(tau)^2+A_4log(c_f)^2+abs(m_c)^2. (1)

G(S,T) is exactly the original-variable mass-extended disk norm
in the preceding theorem, including its specified e_m column.
Then

    Ynorm(Rcal^-1(S,T)) <=10^25000 G(S,T),
    G(Rcal(tau,c)) <=10^40000 Ynorm(tau,c).             (2)

The two different mass representatives are quantitatively compared
below; (2) is not a definition of one norm by pulling back the other.
For the actual compatible inverse on ker L0, the global disk result
consequently gives

    Ynorm(Rcal^-1 A_m^-1 f)
       <=10^412000 lambda^-2 gamma0^-1/2 Fdisk(f).     (3)

The sharper independent bound on abs(m_c) remains available.
Equation (3) quantifies the fast fourth/log TT remainder, the
scalar fourth norm, the mass coefficient and the fixed return.
It does not evaluate the actual comparison residual after inversion
or prove a useful Newton test. All estimates are at the analytic
seed, not at an arbitrary merely C5 candidate.

## 2. Weighted projection energy, including the bolt

Use physical proper radius rho and
A=a'/a, R=R_phys'/R_phys, w=A+2R, W=u^-2. The radial vector
projection is exactly

    P Yv=(3/2)(Yv''+wYv')-VYv=T'+4RT,
    V>=1.49lambda+(3/2)A^2+3R^2.

These are [the retained projection and mass guards](tt-mass-numerical-bound.md).
The [weighted return polynomial](weighted-york-source-conversion.md)
also gives A>=.995R>=0 on the whole actual exterior.
Put Z0=WYv. Integration by parts in the physical radial measure gives

    E=Re <W Yv,W(-P Yv)>
      =(3/2)||Z0'||^2+integral(V-6R^2)|Z0|^2
      =(3/2)||Z0'+(3/2)RZ0||^2
         +integral[V+ (9/2)AR-(57/8)R^2]|Z0|^2
      >=(3/2)||Z0'+(3/2)RZ0||^2
        +1.49lambda||Z0||^2+(3/2)||AZ0||^2+.35||RZ0||^2. (4)

The identity uses R'=AR-R^2. The last R coefficient is .3525
before rounding. The source side is exactly

    E=Re <WT,Z0'+AZ0>.

Weighted Cauchy-Schwarz in the three terms of (4) costs
2/3+2/3+(9/4)/.35=163/21<8. Consequently

    ||sqrt(lambda)WYv||<3||WT||,
    ||WYv'||<17||WT||, ||WAYv||<3||WT||,
    ||WRYv||<5||WT||.                                (5)

The weighted coercive form constructs the solution in the fast
weighted energy space. Its critical realization agrees with the
original projection by uniqueness. The regular bolt endpoint and
decaying outer endpoint are preserved. Equivalently obtain (4) on
the smooth core and pass by weighted form closure; no value of the
singular polar coefficient A is bounded by a constant.

For the inverse return, the already proved weighted integral
estimates give, for Yv=-eta,

    ||sqrt(lambda)WYv||<16||Wc||,
    ||WYv'||<=22||Wc||, ||WAYv||<=11||Wc||,
    ||WRYv||<35.2||Wc||.                             (6)

The bolt-to-radius TT integral determining its radial component,
and the outer-fixed eta integral, are unchanged. Their r_phys^-4
tail is not discarded. Bounds (5),(6) supply the low norms for
the following higher estimates independently of those estimates.

## 3. Smooth geometric coefficient bounds without a polar singularity

Use the dimensionless metric gbar4=lambda g4 and its covariant
derivative nabla. All L2 norms still use the fixed physical volume;
the constant volume rescaling cancels in operator bounds. Put
qvec=grad log R_phys, P2 the orbit projector and Q2=I-P2.
Their tensor identities are

    nabla u=-u qvec,
    nabla qvec=a P2+b Q2-qvec tensor qvec,
    nabla_a(P2)_bc=(Q2)_ab(qvec)_c+(Q2)_ac(qvec)_b,
    nabla Q2=-nabla P2,
    a=F_r/(2lambda R_phys), b=F/(lambda R_phys^2).

These tensors are smooth at the bolt. The curvature has the eight
projector-product terms with sectional coefficients

    K1=-F_rr/(2lambda), K2=-a, K3=u^2/xi-b,
    xi=lambda B^2.

Writing b=f(u)=1+nu u^2-2alpha u^3+beta u^4-2e u^5+e u^6,
the accepted guards are nu<.001, alpha<2.002, beta<3.002,
abs(e)<1e-15. Thus a=f-u f'/2 and
K1=-f+u f'-(u^2/2)f''. Every polynomial above and its first six
(u partial_u) derivatives has coefficient sum below B0=10^9.

The displayed tensor differential system and Leibniz imply, by
induction, jet bounds B0 K0^n n!, K0=10^10, for qvec,P2,Q2
and the required scalar polynomial families through n=6.
Indeed a derivative of a sum of at most three products costs
3(n+1)B0^2 K0^n n!, bounded by B0 K0^(n+1)(n+1)!.
The curvature products cost at most
16 binom(n+2,2) B0^3 K0^n n!, including Ricci contraction.
It follows that

    ||nabla^j qvec||_infinity,
    ||nabla^j Rm||_infinity, ||nabla^j Ric||_infinity <10^110,
                                                       0<=j<=6. (7)

This derives global seed coefficients from smooth tensors; no
derivative of a singular radial frame at the bolt is used.

## 4. First-order systems and higher weighted estimates

View Yv partial_rho as a vector field, and write Z=WYv as a
geometric vector (its gbar4 norm is sqrt(lambda)W abs(Yv)).
Let D=div Yv. The exact radial York identity is

    P Yv=(3/2)grad D+2Ric(Yv).

For extraction, set U=WT and Psi=W(D-2T/3). Then

    grad Psi=2qvec Psi+(8/3)qvec U-(4/3)Ric(Z),
    div Z=(2/3)U+Psi+2<qvec,Z>, curl Z=0,
    Wc=-(2/3)U+Psi/2-2<qvec,Z>,
    W(tau-theta)=-Psi/4, theta=S+T/3.                (8)

For reconstruction, set U=Wc and Psi=W(D+2RYv+c). The exact
inverse return instead gives

    grad Psi=4qvec Psi-6qvec U-2Ric(Z)
                                      -12qvec<qvec,Z>,
    div Z=Psi-U, curl Z=0,
    Wtheta_c=Psi/6,
    WT_c=-(3/2)U+Psi/2-3<qvec,Z>.                    (9)

These identities exhibit the leading constants explicitly. They
follow by substitution from the projection and return, not from
a principal-symbol replacement. Equations (5),(6) imply in both
cases ||Z||+||Psi||<500||U||. Multiplication by the radial W keeps
Z curl-free, including at the regular bolt.

Here is a quantitative global Hodge estimate used to differentiate
(8),(9). For a curl-free vector Z and 0<=j<=5,

    ||nabla^(j+1)Z|| <=2||nabla^j div Z||
                          +10^120 sum_(l=0)^j||nabla^l Z||. (10)

For j=0 this is the integrated Bochner identity. For j>=1 apply
the following integration-by-parts estimate to U0=nabla^(j-1)Z:

    ||nabla^2 U0|| <=||Delta U0||
          +8(j+1)(1+K)(||nabla U0||+||U0||), K=10^110.

In components [Delta,nabla] on a rank-l tensor has contraction
bound 16(l+1)[K|nabla U0|+K|U0|]. Integration once gives the
displayed second-derivative estimate. Commute Delta past j-1
derivatives and use Delta Z=grad div Z+Ric(Z). The remaining
coefficient is bounded by
2^(j-1)+16(j+1)(j-1)2^max(j-2,0)+16(j+1), times 1+K,
well below 10^120 for j<=5. This proves (10). Cutoff limits
justify the integrations on the complete exterior; induction
controls their errors. The bolt is an interior point.

All coefficient jets through five in (8),(9) are <10^230 by
(7), including the quadratic qvec term. If N_j is the sum of
Z and Psi derivative norms through j, and I_j the corresponding
sum for U, product rule and (10) give

    N_(j+1)<=10^235(N_j+I_j), 0<=j<=5,
    N_n<=10^(240n+3) I_max(n-1,0), 0<=n<=6.          (11)

Substituting in the output formulas proves bounds <10^2000 at
each integer Sobolev order s=3,5 for c and tau-theta from T,
and theta_c,T_c from c. Here the weighted scalar norm is the
square sum of covariant derivatives of W times the scalar;
dividing both sides by sqrt(gamma0) leaves the bounds unchanged.
The low norm, not an unproved higher inverse, starts the induction.

## 5. Transfer to the exact fourth/log atlas norms

Integer endpoint bounds must be transferred before applying a
logarithmic interpolation argument. Let H_s be the ambient direct
sum defining J. There is a radial extension E such that EJ F=F:

    Eh=average_time[ zeta_0 Y^delta h_outer(log Y)
                    +sqrt(A0) sum_i f(s)chi_i h_i(z_i) ].

Use exactly the partition f(s),zeta_0 in the global disk proof.
The first term is defined only on its outer support; the other
terms are zero extended through their chart boundaries. Time
averaging uses the actual period and is contractive in Hcov^s,
because it averages isometries preserving W. At the bolt it is
the smooth rotation average. For a stationary F the partition
identity and chi=1 on supp zeta_0 give EJ F=F exactly.

For s=3,5 both endpoint conversions have bounds

    ||JF||_H_s <=10^10000 ||WF||_Hcov^s/sqrt(gamma0),
    ||W Eh||_Hcov^s/sqrt(gamma0)<=10^10000 ||h||_H_s. (12)

In this display Hcov denotes the unnormalized covariant derivative
sum, so the factors sqrt(gamma0) are shown explicitly. Here is a
finite reserve for (12). On the interior cover u>r/24>10^-164;
the canonical/physical conformal factor is k/u, with .4<k<.5.
The chart jets, partition jets, and their coordinate transitions
are those already quantified in the global disk proof. In the
normalized metric the orbit metric is (k/u)^2 gbar and the sphere
radius is sqrt(xi)/u. Use u^-1<10^164, t^-1=10^170, the previous
partition/transition chain base 10^550 and the smooth tensor jet
bound 10^110. Their product is 10^994<10^1000. The conformal
connection formula and the P2,Q2 identities express each next
scalar derivative as a derivative of one coefficient or a product
with one connection coefficient. Their Leibniz sums cost at most
10(j+1) at order j<=5. Thus every derivative-conversion coefficient
through five
is bounded by 10^1000 j!(10J0)^j, J0=10^1000. Product/Bell
sums through five fit the same reserve after enlarging 10 to 100.
The squared weighted physical measure is exactly

    gamma0^-1 dmu4 W^2=(k/A0)u^-8 dvol_gbar.

Its square-root conversion on a patch is <10^1000. Forward and
reverse patch overlap is <10^6. Thus a sufficient endpoint budget
is 10^1000*4sqrt(10^6)*3^s times the maximum of those finite
derivative sums, strictly below 10^10000. This also covers the
sphere derivative terms through the smooth qvec,P2 identities.
At the end, u=kY w_seed, j=a/a(0) and dvol_gbar=j dY dy give
the exact weighted density k^-7 w_seed^-8 j dx. Its conormal jets
and the relation partial_rhohat=-w_seed partial_x satisfy the same
reserve. No period bound or polar connection bound is assumed.

For any of the scalar maps in section 4 apply J K E to the
ambient spaces. Its H3 and H5 bounds are <10^22000. The
dimension-independent dyadic Schur proof from the global disk
estimate gives H4log bound 28672 times this, and H4 bound 1536
times this. Enlarge to 10^23000. Since EJ is the identity on the
actual radial space, these are bounds on the actual maps in A_4
and A_4log. The mixed scalar/shape and b=2T conversions fit
10^24000 in each direction on the FAST spaces. This does not
infer an H4-to-H4log estimate for an arbitrary scalar input.

## 6. Explicit comparison of the two mass representatives

The physical unit column is theta_e=-m_0^-1 c_LY^3/12,
T_e=m_0^-1 c_LY^3/2, S_e=-T_e/2, m_0=sqrt(j).
Choose the explicit comparison projection vector, in outward
dimensionless proper-radius component, q_e=-T_e/6. It vanishes
near the bolt. The actual projection differs by a vector H solving

    Pbar H=E_e=T_e'+4Rbar T_e-Pbar q_e,

where primes now denote outward dimensionless proper derivatives.
For D=Y partial_Y, w_seed=u/(kY), Rbar=j,
Abar=j-w_seed Dj/j, and Vbar the dimensionless projection potential,

    E_e=[w_seed^2 D^2/4
       +(w_seed Dw_seed/4-(1+(Abar+2j)/4)w_seed)D
       +4j-Vbar/6]T_e.                              (13)

The normal operator is (D-3)(D-4)/4. In particular the uncut
mass cubic cancels before taking a fast norm.

Here are explicit reserves for the source and local mass error:

    ||W E_e||_Hcov^4/sqrt(gamma0)<10^1000 L^-1/2,
    ||W tau_*||_Hcov^5/sqrt(gamma0)
      +||W(c_*-c_0)||_Hcov^5/sqrt(gamma0)<10^1000 L^-1/2, (14)
    tau_*=-(div q_e)/4,
    c_*=-T_e+(div q_e)/2-2Rbar q_e,
    c_0=-c_L u^3/(3k^3).

To verify (14), the retained analytic collar has |j-1|,|w_seed-1|
<.001 on |Y|<.1. Cauchy on the smaller disks bounds j',w_seed'
by .02 and j'' by .8 on |Y|<=.02. The displayed Abar has modulus
<1.002 (in particular <3) and derivative <.1 there;
Vbar=(partial_rhohat(Abar+2j))/2
+2Abar^2+4j^2 has modulus <10. On |Y|<.01 all coefficient
differences from their normal values vanish through degree one.
On the required smaller interval 0<Y<=2L<.001, their derivatives
through D^6 are <10^30 Y^2, using the convergent majorant
10*100^2 sum_(n>=2)n^6(.1)^(n-2). The same bound holds for
m_0^-1-1 and w_seed^3-1.

The beta cutoff has log jets through six <10^14. Product rule
in (13) bounds all D^j E_e, j<=4, by 10^100 Y^5 below L and
10^100 Y^3 on L<Y<2L. The same bounds through j=5 hold for
tau_* and c_*-c_0; the latter uses w_seed^3-m_0^-1=O(Y^2).
The sixth piecewise derivative of the C5 cutoff is sufficient.
The collar's unit radial vector n=qvec/j and -w_seed n have
covariant jets <=j!10^(60(j+1)) through five, by (7)'s sharper
recurrence and these analytic scalar bounds. In differentiating
the weighted fields use nabla(WF)=W(nabla F+2qvec F).
The finite chain rule, including these weight derivatives, uses
products of at most five geometric jets with total derivative
order at most five. The displayed jet bound costs at most
10^600 times the finite factorial/Bell sums, below 10^700.
Finally the weighted density is
<10^5 Y^-8 dY: the two integrals are L^3/3 and 1/(2L).
These bounds prove (14) with reserve.

For H the coercive form (4), now with vector forcing E_e, gives
||Z_H||+||Psi_H||<500||WE_e||, where
Z_H=WH and Psi_H=W div H. Their exact equations are

    grad Psi_H=2qvec Psi_H+(2/3)WE_e-(4/3)Ric(Z_H),
    div Z_H=Psi_H+2<qvec,Z_H>, curl Z_H=0.

The same induction through input order four bounds both fields
through five. The corrections to (tau_*,c_*) are
-Psi_H/(4W) and [Psi_H/2-2<qvec,Z_H>]/W. Equations (11),(12),
(14), and H4log<=2H5 in the ambient spaces consequently give

    A_4(tau_e)+A_4log(c_e-c_0)<10^13000 L^-1/2.      (15)

This is a bound on the actual projected physical column, not an
arbitrary induced norm. Its outer mass is -1/(3k^3), exactly
as required by m=-3k^3 m_c. Both columns have their original
regular-bolt realizations; the projection's return tail is kept.

## 7. Assembly, verification and remaining work

Extraction applies the fast bound 10^24000 to the disk fast
remainder, then adds m times (15) and the mass coefficient.
Use L^-1/2<10^82 and abs(m_c)<6abs(m). This gives the first
inequality of (2). For the reverse map, the mass discrepancy is

    Rcal(0,c_mass)+3k^3 e_m
           =3k^3 Rcal(tau_e,c_e-c_0).

Its fast norm is bounded by 10^24000 times (15). Since 3k^3<1,
the resulting reverse mass constant is below 10^37082; all terms
fit 10^40000. The preceding global disk inverse then proves (3).

Freeze proof and checker before execution. Verify the weighted
square completion and coercivity, both exact first-order systems,
global polynomial/tensor jet budgets, commutator/Hodge coefficient
counts, higher induction, endpoint and logarithmic transfers,
normal mass cancellation, cutoff/source/integral reserves and
final two-way mass conversion. Pin the source receipts and hashes;
preserve every run and failed source snapshot. Finite checks support
the ordinary global operator, geometry and closure arguments.
They do not formally verify those arguments or a nonlinear solution.

Remaining to-do list: actual residual/inverse composition and S166
completion audit; nonlinear family and signed compatibility;
exterior/state, charge and return; opening, stability, transport,
preparation and separate creation requirements.
