# Near-diagonal restoration and a quantified integration remainder

2026-10-23, E003/S169. Ordinary argument before execution.
The full stable semiclassical wormhole objective remains incomplete.

This advances the actual differentiated finite-domain solve, not the
parked cold-value estimate chain. It identifies and evaluates the common
fourth-order local restoration on the current finite seed, and bounds
the short-distance remainder of a Taylor-subtracted integral acting on
the actual reference direction. It also tests the regularity and mode
resolution needed before that integral can supply a differentiated
metric solve. It does not calculate all local moments, the complete
renormalized Hessian, or the physical correction.

## 1. Fixed action, region and direction

Keep G=1 only as the unit, k=2000, q=k^4, c_R=1e22, the original
mu0=1/b0, preparation z_i=.265, exact product parameter enclosure and
paired finite-seed realization of
[the executed full separated sum](current-seed-separated-quantum-response.md).
Its canonical proper coordinate is rho, z=B/r and

    Q=(1-z)^2[beta(1+2z+3z^2)+u z^2+epsilon z^4]-tau z^3,
    partial_rho=-sqrt(Q) partial_z, a=sqrt(Q)/gamma,
    R2=-Q'', Delta=Q partial_z^2+Q' partial_z.

Primes on Q mean z derivatives. Work at observation z in [1/2,7/8]
and proper separations |t|<=.01. Since Q<.2, the connecting intervals
stay inside J=[.49,.88]. On J direct polynomial bounds give Q>.01.
The inner cut z=15/16 is separated from J by proper length at least
log(.12/.0625)>.25. All directions are tested with compact interior
cutoffs when interpreting distributions. No statement about closed
endpoint C5 regularity, unequal mouths or the retarded domain follows.

The actual reference shape direction u_h is the u component of h1 in
[the cold action defect](current-seed-cold-action-defect.md). Distinguish
it from the constant u=1/C in Q. Its retained fixed seed chart bounds are

    |u_h|<1.2e-8, |partial_x u_h|<8e-6,
    |partial_x^j u_h|<=1.3e-4 j!100^j, 2<=j<=5.

These follow from [the current path transfer](current-seed-directional-cold-response.md)
and [its normalizer calculation](comparison-path-interior-response.md).
The analytic seed normal chart satisfies sigma(0)=sigma'(0)=0 and
|sigma^(j)(0)|<=.027(j+1)j!32^j. The radial axis can be chosen real;
proper differentiation along it is D=exp(-sigma) partial_x. At the center,

    D2 u_h=u_h'', D3 u_h=u_h'''-sigma'' u_h',
    D4 u_h=u_h''''-4sigma''u_h''-sigma'''u_h',
    D5 u_h=u_h'''''-10sigma''u_h'''-5sigma'''u_h''
                        +(7sigma''^2-sigma'''')u_h'.       (1)

The common normal cover permits these bounds at every point in J. The
executable substitutes the stated constants rather than identifying
divided chart jets with raw proper derivatives. In particular H5=2e8
is a sufficient upper bound for sup_J |D5 u_h|. This is a bound on a
specified direction, not on the unknown correction or its nonlinear ball.

## 2. The prescribed common fourth-order contact

The retained common angular finite part, not a separate mode scheme, is

    F_eff=4gamma_E-2H_q+4log(mu0 B),
    alpha=-C F_eff/(960pi),
    Gamma_a2[b]=alpha integral dvol_gamma exp(-2b)(R2-2Delta b)^2. (2)

This is the same restoration as
[the first-order vertex calculation](flat-variable-dirac-three-point.md).
Use H_q=psi(q+1)+gamma_E and the actual B/b0 interval. The c_R term,
Weyl transgression, zero-momentum vacuum restoration and boundary action
remain additional retained terms. None is reset by (2).

The affine shape Hessian of (2), per seed area, is exactly

    H_a2 v=8alpha[Delta^2 v+R2 Delta v+Delta(R2 v)+R2^2 v/2]. (3)

This includes the density, curvature and gradient contacts. With the
unchanged row normalization -2G/B^2, set

    c4=GC F_eff/(60pi B^2),
    H_a2,norm v=c4 [Q^2 v''''+4QQ'v'''
      +(QQ''+2Q'^2)v''-(QQ'''+Q'Q'')v'
      +(-QQ''''-Q'Q'''+Q''^2/2)v].                    (4)

For a bound in proper derivatives put H=a'/a=-Q'/(2sqrt(Q)). The
bracket in (3) is sum_(j=0)^4 c_j D^j with

    c4_local=1, c3_local=2H,
    c2_local=2DH+H^2+2R2,
    c1_local=D2H+H DH+2R2 H+2DR2,
    c0_local=Delta R2+R2^2/2.                         (5)

Do not confuse c4_local=1 with the overall c4 in (4). Differentiate
(5) once, using D=-sqrt(Q) partial_z. Coefficient enclosures on all of
J and (1) then bound both H_a2,norm u_h and D H_a2,norm u_h. This
operation needs only five direction jets. It is not differentiation of
a previously asserted source-value bound.

## 3. Leading singularity and the information missing from separated data

For a local flat chart at fixed M and proper thermal momentum w, the
direct first-order bubble is -w^2/(M^2+w^2) exp(-2sqrt(M^2+w^2)|t|).
The ultraviolet angular density is exactly obtained by continuous
interpolation of M_n^2=n(n+q)/C and d_n=2(q+2n):

    d_n dn=2C d(M^2)=4C M dM.

The two thermal signs and per-area normalization give the continuum
leading kernel

    -(4C/pi) integral_0^infinity integral_0^infinity
       [M^3 w^2/(M^2+w^2)] exp(-2|t|sqrt(M^2+w^2)) dM dw
       =-2C/(5pi |t|^5).                             (6)

Polar integration uses integral sin^3(theta)cos^2(theta)dtheta=2/15
and integral k^4 exp(-2|t|k)dk=3/(4|t|^5). This is the LEADING local
symbol, not a replacement for the finite-seed kernel. The ordinary
interior heat localization in
[the logarithmic-response proof](stationary-shape-logarithmic-response.md)
justifies this leading transfer; lower state/curvature contributions
are retained. The image term is smooth on this compact interior.

After the row normalization, the coefficient is

    kappa5=4GC/(5pi B^2)>0,
    K_norm(t)~kappa5/|t|^5.

Since partial_t^4(1/|t|)=24/|t|^5 off zero and the Fourier finite part
of 1/|t| has logarithmic coefficient -2log|k|, this gives the principal
multiplier -GC k^4 log|k|/(15pi B^2). Independently, differentiating
(4)'s c4 with respect to log(mu0) gives +GC/(15pi B^2), the required
local scale counterpart. No finite scale is changed in this test.

The off-diagonal distribution has scaling degree five in one radial
dimension. Its extensions differ by local derivatives of delta up to
order four, restricted by covariance and symmetry. Thus even exact
separated kernels cannot determine all local moments. Adding a local
fourth-order operator changes the metric equation without changing any
of those separated values. Equation (2) fixes the common a2 allocation;
the other finite local/state-dependent moments still need evaluation
in the same prescription. Setting the missing moments to zero is not
a permissible numerical convenience.

## 4. An all-mode near-diagonal majorant with the correct power

The preceding first-order comparison supplies, for all massive modes
and both direct/image paths in J,

    |B_direct+B_image|<=8 exp[-c |t|sqrt(M_n^2+5p^2)], c=1.98.

Unlike the earlier far-separated tail, do NOT bound M_n only by
18sqrt(n) here: that bound loses the correct ultraviolet spectral
density. Put a0=c/sqrt(2), b0_t=c sqrt(5/2). For fixed p the angular
summand

    f(x)=2sqrt(1+4C M(x)^2/q^2) M(x)^2
                         exp[-c|t|sqrt(M(x)^2+5p^2)]

is nonnegative and unimodal. Its logarithmic derivative in M is
2/M+dM/(1+dM^2)-c|t| M/sqrt(M^2+5p^2), d=4C/q^2; this derivative
is strictly decreasing. The first two terms have negative combined
derivative, and the last term is decreasing. Therefore
sum_(n>=1) f(n)<=integral_0^infinity f(x)dx+2sup f.
The extra factor two is a harmless variation allowance.

Thermal monotonicity gives gamma sum_(j>=0) g(gamma(j+1/2))
<=integral_0^infinity g(p)dp+gamma g(0). Keep the actual gamma.
The main joint integral, with the exact angular Jacobian, is
64C/[q sqrt(5)(c|t|)^5]. For the angular-sum error use
sqrt(M^2+5p^2)>=(M+sqrt(5)p)/sqrt(2) and
sqrt(1+4CM^2/q^2)<=1+2sqrt(C)M/q. The maxima of
M^r exp(-aM) are (r/(ae))^r. These prove the complete majorant

    |K_norm(A,A+t)|<=sum_(j=2)^5 kappa_j/|t|^j,      (7)
    D0=16Gq/(pi B^2 sqrt(Q_min)), Q_min=.01,
    kappa_5=D0 64C/[q sqrt(5)c^5],
    kappa_4=D0 [216sqrt(C)e^-3/(q a0^3 b0_t)
                                      +gamma 24C/(q c^4)],
    kappa_3=D0 [16e^-2/(a0^2 b0_t)
                                      +gamma 216sqrt(C)e^-3/(q c^3)],
    kappa_2=D0 gamma 16e^-2/c^2.

Here the majorant kappa_5 is DIFFERENT from the exact leading kappa5
in section 3. It includes comparison losses and the image term. Every
angular level and thermal number is covered in (7); there is no
zero-temperature replacement or assumed ultraviolet cutoff.

Let T4,A v be the degree-four proper Taylor polynomial at A. Then
|v(A+t)-T4,A v(A+t)|<=H5 |t|^5/120. The omitted interval of the
Taylor-subtracted integral therefore has the rigorous VALUE allowance

    E_delta <= H5/60 sum_(j=2)^5 kappa_j delta^(6-j)/(6-j). (8)

This is not the entire renormalized response: its polynomial moments,
cutoff terms and local restoration must be put back in the prescribed
scheme. Nor can a derivative be taken through (8) without additional
regularity. Quantify (8) at delta=.01,.0001,.00001,.000001.

## 5. The differentiated obstruction and its quantitative target

A C5 label alone does not bound the full response in C1. Interior
packets v_N=N^-5 chi exp(iN rho) have bounded C5 norms, whereas
D H v_N has the nonzero principal growth const*log N. Ordinary
order-four contacts, order-three-log remainders and boundary smoothing
cannot remove that leading divergence. This does not say the actual
smooth u_h lacks the additional regularity; its needed quantitative
modulus has not been supplied by a C5 bound.
This is a mapping obstruction for the retained uncut operator, not a
claim that arbitrarily large packet frequencies lie in a proved
semiclassical validity band. A restricted band would need its own
explicit physical and operator error control.

For the principal singular part alone, a fifth-derivative Holder
seminorm H_(5,alpha) gives the centered local remainder bound

    E_derivative,principal <= kappa5 H_(5,alpha) delta^alpha/(12alpha).

For alpha=1/2 and the diagnostic tolerance epsilon_num=1e-28, the
required seminorm threshold is

    H_(5,1/2)<=6 epsilon_num/(kappa5 sqrt(delta)).     (9)

Calculate that threshold at the same cutoffs. The seminorm remains
UNKNOWN. Lower symbols, moment derivatives and endpoints require their
own bounds; (9) is not a sufficient full-operator or metric certificate.
A Dini modulus is a possible sharper alternative: its integral
integral_0^delta omega_5(t)dt/t must be controlled. This is the regularity
question to settle before continuing differentiated quadrature.

There IS an older higher-jet allowance; it must not be called unknown.
[The original finite-middle source proof](finite-middle-source-four.md)
uses |z0^j nabla^j u|<=4e24 j!100^j through j=9, from
[the original high-jet ledger](seed-nine-metric-eleven.md). At z0>=.49,
its sixth proper radial derivative allowance is

    H6_old=4e24*6!*100^6/(.49)^6.

It belongs to the original comparison. The current c_R reference-path
transfer cited in section 1 explicitly supplies C5; it does not by
itself transfer this higher ledger. Conditional on such a transfer,
the elementary interpolation modulus min(2H5,H6_old |t|) would give
H_(5,1/2)<=sqrt(2H5 H6_old), with H5=2e8. Calculate its ratio to
(9) at delta=1e-5. That known coarse allowance fails to certify the
selected principal quadrature tolerance by many orders, even if it
transfers. It is not a lower bound on the actual seminorm or an
obstruction to the actual solution. A sharper direct derivative or
modulus calculation is required; merely labeling the path smooth or
reusing this huge reserve does not finish it.

For a concrete resolution diagnostic also compute the angular index
and proper thermal momentum at 32 exponential attenuation units:

    M_32=32/(a0 delta),
    n_32=2C M_32^2/[q+sqrt(q^2+4C M_32^2)],
    p_32=32/(b0_t delta), j_32=p_32/gamma-1/2.        (10)

These are scale counts, NOT an asserted truncation-error tolerance.
They assess whether extending the few-level separated evaluator
directly to the diagonal is practical. Exact binning can group indices;
the mode counts do not prove that every mode must be evaluated alone.
The local parametrix must capture the ultraviolet part and match the
resolved finite-geometry modes with an error bound.

## 6. Execution contract and decision boundary

Use the unchanged parameter box and python-flint 0.9.0 / Arb, 256 bits.
Enclose geometry coefficients on 64 and 128 cells covering J. Require
all Q>.01, agreement of contact and proper-jet algebra, and agreement
of the ultraviolet polar and scale coefficients. Check H5<2e8, the
common-contact value <1e-27, its first derivative <1e-24, and (8) at
delta=1e-5 <1e-28. These are predeclared component precision/size
tests, not inferred physical solution thresholds. Keep every run,
failure, numerical result, source snapshot and dependency hash.

After the first successful component run, the source audit located the
older higher-jet ledger above. Preserve that successful snapshot, add
the explicit conditional comparison and its provenance, and replay the
whole final package. This is a documented scientific addition after a
successful run, not an unrecorded edit to its certification scope.

The [Arb API](https://python-flint.readthedocs.io/en/latest/arb.html)
supplies enclosed elementary functions, digamma and Euler's constant.
The continuum argument above is separate from finite symbolic checks.
No Lean or new external compute is selected.

Continue the actual H e=-r1 route, but prioritize a matched local
parametrix with the prescribed moments and a quantified fifth-derivative
modulus of the actual direction. Do not extend a few-level separated
quadrature toward the diagonal without this subtraction and regularity
step. This execution evaluates the common fourth-order contact and
short-distance VALUE remainder. It does not execute the remaining local
moment calculation, full differentiated response or metric inverse.
All closed-audit nonlinear, compatibility, boundary, matching, validity
and stability targets remain unchanged. No oversized allowance or
regularity counterexample proves that the actual wormhole cannot exist.

Remaining to-do list: fixed-scheme local moments and matched parametrix;
quantified actual-direction derivative modulus and full weighted
response; actual correction and nonlinear admission; compatibility,
matching, semiclassical validity, causal realization, all-sector
stability and every remaining full stable-wormhole requirement.
