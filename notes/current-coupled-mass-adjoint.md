# Current coupled mass adjoint and its quantitative residual target

2026-10-25 HST, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: NO PROGRESS. Full stable-wormhole objective incomplete;
the closed feasibility audit and original physical targets are unchanged.

OUTCOME: the requested actual coupled return is still UNRESOLVED. This
attempt identifies a current-theory mass certificate that requires weak
primal/adjoint residuals, rather than an unknown energy-to-H3 comparison.
It establishes the fixed-seed auxiliary domain and the complete mass-flux
normalization. It does NOT evaluate the signed pairing or either residual.
The numerical targets below are sufficient allowances, not computed errors.

## 1. What is being attempted, and what is not transferred

The full conversion calculation gives 2133<||d_chi||H3<2174 and
g_old=g_lifted+m d_chi. A sufficient CONVERSION-ONLY split is

    ||g_lifted||H3<=5e-12, actual |m|<=2.29e-15.             (1)

Neither input is known. The older scalar-TT mass-adjoint certificate is
for a different operator and cannot be inserted into the current solve.
Here use the COMPLETE current H, including both original canonical rows,
the c_R=1e22 action, all angular/thermal terms, local contacts, off-shell
York/density returns, fixed clock, actual paired mass/log domain and
regular bolt. The reference mass lift and its separate .1/.2 cutoff are
unchanged. Its local cancellation is not a global zero-mode assertion.

The strong common-regulated core is still the missing step for the FULL
graph return. For the scalar observable m, first test whether the current
weak inverse and exact boundary flux avoid that step. They do, with an
additional affine source pairing that must actually be evaluated. There
is no resulting bound for mass on the unrestricted energy space alone.

## 2. The original units and current two-row boundary flux

Let H_c denote the canonical density Hessian, A=z^4 H_c, and beta=x/C.
Use the physical positive-diagonal energy E of
current-global-return-energy-attempt.md. The exact conversion in
current-coupled-input-energy-bound.md, section 2, is

    E(u)^2=(a+l)(u,u)/(|W| beta),
    h(u,v)=I''(u,v)/(|W| beta)
          =-beta^-1 integral u^T H_c v dz.                  (2)

Indeed W=-N pi B^2/(G gamma), beta_T=2pi B/(C gamma), and
|W|/(N 4pi beta_T)=CB/(8piG). Both forms and sources use (2).
Write h u=F for the ACTUAL full linear equation, with its original sign.
The current coupled saddle inverse satisfies ||h^-1||_(E*->E)<=1.
The source calculation gives ||F||E*<1.1e-10; its proof bounds the source
dual norm before applying the inverse, not just the resulting solution.

Keep the exact full normal matrix from current-coupled-outer-return-test.md:

    A0(s)=[[Jc, Jc a(s)],
           [a(3-s)Jc, a(3-s)Jc a(s)+4D0 N(s)]],
    a(s)=(s+3)/(3(s+1)),
    Jc=-6beta[C_v-D_R lambda s(s-3)](s-4)(s+1),
    D0=-8pi g beta, g=G lambda,
    N(s)=d0 q(s)(q(s)+2)H(s)-A_T q(s), q=s(s-3),
    d0=1/(480pi^2), A_T=Fstar/(48pi g),
    Fstar=F_lambda+23g/(180pi),
    H(3)=1-gamma_E+kappa_M,
    kappa_M=-6/5-log(mu0/sqrt(lambda)).                    (3)

The outer root-zero column is p=(-1,1), and the original unit mass is
v_m=(-1/2,1) at weight three. The correct dual identities are

    A0(0)p=0, p^T A0(3)=0, A0(3)v_m=0,
    p^T A0'(3)v_m=4D0 N'(3).                             (4)

The scalar and mixed terms cancel in the last identity. A left-null
test at weight ZERO would be the wrong test: the adjoint reflection is
s -> 3-s. This corrects an exploratory test in the preceding goal turn;
no failed left-null identity is used here.

With the outward lower-end sign and the minus sign in (2), the positive
mass flux is

    K_m=4D0 N'(3)/beta
       =2Fstar-(2g/(5pi))(1-gamma_E+kappa_M)
       =2 Z_c,
    Z_c=F_lambda+(g/(5pi))[
               log(mu0/sqrt(lambda))+gamma_E+151/180].     (5)

The last expression is EXACTLY the already retained physical boundary
charge in recentered-mass-boundary.md, section 5. It is a separate check
of the complete contact and physical normalization, not a new choice of
charge or counterterm. The scalar quartic row 20beta Fstar is still in H;
it is not removed by taking the differentiated shape equation.

## 3. Auxiliary domain and common-regulated Green limit

Take P(z)=p on the physical canonical disk. It is regular at the bolt:
nu+b=0 leaves the orbit metric unchanged and the constant sphere variation
is smooth there. P is NOT in E at the outer end. It is an auxiliary test,
not an admitted physical correction or a new free boundary source.

Define q_P=h P initially against compact physical-half tests, using the
actual fixed-seed kernel and its retained contacts. The current weighted
kernel argument applies here at a weight just below zero: the normal
symbol is holomorphic there, the paired one-sided coefficients are analytic,
and the separated kernel weights require -1<delta<4. Choose delta=-1/8.
P belongs to that auxiliary weighted graph, and A0(0)p=0. The two-power
actual-minus-normal gain therefore makes z^4 H_c P=O_b(z^2 log^N z).
Any separated finite-domain moment starts at weight four and is kept.
Thus the physical source represented by q_P is in critical physical L2:
its outer squared integrand is O(log^(2N)z), which is integrable in dz.
The genuine York transpose is bounded at this critical weight, and the
positive scalar/TT gaps embed this source in E*. At the regular bolt the
fixed-seed compact/disk kernel argument supplies the same membership.
This proves finiteness, not a useful numerical value for ||q_P||E*.

These are exactly the fixed-seed direct/reflected kernel and complete
contact hypotheses of current-seed-full-hessian-realization.md sections
3 and 5, recentered-reference-inverse.md sections 1-3, and the weighted
Green argument in stationary-weighted-tt-inverse.md section 4. No older
inverse constant, smooth doubled seam, new boundary realization or
neighborhood regularity estimate is transferred. The normal family is
the CURRENT one in (3). The proof uses the already established two-power
gain, not the unknown numerical higher differentiated core.

For clarity, do not interpret q_P as a metric derivative of the action on
an enlarged physical nonlinear ball. Use P_e=chi_e P, where chi_e is zero
for z<=e and one for z>=2e. Each P_e is in the original fixed-source form
core (after ordinary smoothing), with the original regular bolt and paired
projectors. There is no wall at e. At positive common proper-time cutoff
keep the full operator and all local restorations; remove that regulator
for each fixed e using the existing core identity. Only then take e->0.
No unproved simultaneous-regulator rate is used or needed for this identity.

The actual source F has the retained weight-four finite-log asymptotic
expansion. Hence F(P)=lim F(P_e) exists as an ordinary integral of the
COMBINED canonical source; it is NOT controlled by ||F||E* alone. The
weak solution u has

    u=m v_m z^3+O_b(z^4 log^N z),                         (6)

with its correlated scalar/slow terms and regular bolt. In the cutoff
pairing, the actual-minus-normal commutator contribution tends to zero:
the normal mass pairing has weight zero and its coefficient error gains
two powers, giving O(e^2 log^N e); the fast term gains at least one power.
Separated tails obey the same weighted bounds used in the retained Green
argument. Their finite moments are in q_P, not assigned zero.

The surviving normal commutator is K_m m. One way to compute it is to
apply the complete multiplier A0(D) to chi_e p: integral chi_e' dx=1
for x=log z, so its weight-zero pairing is v_m^T A0'(0)p.
A0(3-s)^T=A0(s) changes this to -p^T A0'(3)v_m; (2) gives (5).
The multiplier is holomorphic at the two roots and has polynomial/log
growth; the divided multiplier acting on the smooth compact chi_e'
justifies this zero-frequency identity without a mode cutoff. Contacts
are part of the multiplier. Bolt boundary work vanishes by regular-disk
matching; no punctured-bolt extension is introduced.

Consequently the actual two-row Green identity is

    K_m m=F(P)-q_P(u).                                   (7)

The canonical mass/contact distribution has not been discarded against
the nonzero trace of P. The cutoff commutator retains it, and (5) matches
its independently derived open-clock charge. This is why a bare constant
test with its boundary term set to zero would give a wrong certificate.

## 4. Current coupled adjoint and computable error formula

Set v=h^-1 q_P in E and Psi=P-v. This is the unique affine-domain
solution h Psi=0 with Psi-P in E. It uses the complete current coupled
inverse, not an uncoupled TT inverse. The auxiliary v can contain a
forced weight-two term; it is not proposed as a nonlinear metric change.
Equation (7) and symmetry on E give

    m=F(P-v)/K_m.                                        (8)

Take any bolt-regular approximate forward field u_a and adjoint correction
v_a in the actual energy/form domain for which the full residuals are
certified. Define

    d=h u_a-F, r=h v_a-q_P,
    I_a=F(P-v_a)+r(u_a)
       =F(P)-q_P(u_a)+d(v_a).                            (9)

All pairings in (9) have just been defined: the only non-E input is the
separate convergent source moment F(P). Since h is symmetric and its
coupled inverse has norm at most one, exactly

    m-I_a/K_m = -r(h^-1 d)/K_m,
    |m-I_a/K_m|<=||r||E* ||d||E*/K_m.                    (10)

Thus a numerical enclosure I_hat with certified absolute error E_I gives

    |m| <= (|I_hat|+E_I+delta_r delta_d)/K_m.             (11)

No H3 output norm occurs in this mass error. Every scalar, slow, mass,
interior and bolt coupling is inside h,d,r. This does not allow their
omission from assembly or prove any graph output estimate. Finite trial
representations are permitted only if their errors use the FULL operator,
not a finite angular/thermal selection.

## 5. Quantitative test and limits of what was executed

Evaluate (5) with the current enclosed g,F_lambda and complete contact
range 39/5<kappa_M<54/5. The conservative lower target for verification is

    K_m>1.9999999764.                                    (12)

For the sufficient mass half-budget in (1), the following explicit targets
would suffice:

    |I_hat|<=3.3e-15, E_I<=1e-16,
    delta_r<=1e-5, delta_d<=1.1e-10.                     (13)

The last allowance is already available with u_a=0, because d=-F. The
other THREE quantities have NOT been evaluated. Their full numerator
is <=4.5e-15; (11)-(12) then give |m|<2.250000027e-15<2.29e-15.
With that mass bound AND the separately unproved lifted graph <=5e-12,
the previous upper conversion 2174 gives original graph <1e-11.

More generally, after reserving 3.4e-15 for |I_hat|+E_I, the allowable
product is K_m*2.29e-15-3.4e-15. With delta_d=1.1e-10, this permits an
adjoint residual slightly above 1.0727e-5. This is a conditional error
budget, not an observed residual, an adjoint norm or a mass estimate.
It is less demanding in derivative order than a global H3 certificate;
attainability of the signed full pairing at absolute error 1e-16 remains
undecided. No numerical claim of realistic attainability follows from
these tolerances alone.

EXECUTED: the current coupled normal/charge consistency, physical energy
conversion and resulting residual thresholds, supported by exact symbolic
identities and outward Arb arithmetic. The ordinary cutoff/domain argument
above is required independently of those checks. NOT EXECUTED: the complete
I_hat, r or common-regulated core application. In particular (13) has not
passed as an actual feasibility test. The preceding exploratory wrong-
weight null test is recorded in the review; its corrected version is used.

The original kappa_m<=5e-5 remains unknown. Its conditional remaining
budgets remain Y_other<=5.34234999994958e-12 and
eta_other+(.0195239+4.583e-13)gamma_E<.0012069. The screen
.028913+1e-11*1.0706e10=.135973<.1599164 still has only its proved
two-fast-input admitted-family scope. Neither (11) nor (13) establishes
independent phase derivatives, physical-norm embedding or a nonlinear ball.

RESEARCH DECISION: retain the construction; sharpen the signed admission
method using this CURRENT COUPLED mass adjoint. The first mass decision can
use a weak residual certificate without waiting for an energy-to-H3 theorem.
Do not optimize the already controlled cold/scalar or L2 angular tails.
The main coupled graph milestone remains open, not replaced by this method.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): apply the complete restored
current operator to a bolt-regular adjoint trial v_a, and enclose its
residual r together with the SAME signed source pairing F(P-v_a), testing
(11) against 2.29e-15. The source moment, actual contacts and all finite-
domain returns must be assembled with their common regulator. A failed
upper estimate would not prove the actual mass too large or rule out the
construction; an enclosure excluding the target would reject only this
sufficient conversion split. The original graph still needs its own bound.

Remaining to-do list: executed coupled signed mass/graph certificate;
original physical-norm/neighborhood comparisons; nonlinear admission and
contraction; exact compatibility; physical matching; semiclassical validity
and noise; all-sector stability; charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
