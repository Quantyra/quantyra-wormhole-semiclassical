# Complete coupled application residual: forcing cancellation and original rows

2026-10-26 HST, E003/S169. Ordinary argument and execution contract before
arithmetic. Previous goal turn: PROGRESS (complete bolt map bounded).
The full stable semiclassical wormhole objective remains active and the
feasibility audit remains closed. No new physical source is chosen.

## 1. The first gap in the application estimate

Keep the actual fixed-seed Hessian H, paired mass/log domain and adjoint,
fixed temperature, original finite-part equations, mass cutoff transitions
and regular bolt. All angular and thermal sums, ultraviolet logarithm,
coincidence/density contacts and common regulator removal are unchanged.
Write its physical canonical rows and frozen York congruence as
 H=[-A,-C;-C*,D], R(s,xi)=(s+Uxi,Vxi),
 R*HR=[-A,b;b*,ell], B=b*A^-1b, N=|W_action| beta,
 E(s,xi)^2=[a(s,s)+ell(xi,xi)]/N.
The accepted full inverse has energy-dual to energy norm <=1.

The preceding P32 application diagnostic granted an F factor on the error.
That hypothesis is not supplied by a numerical inverse. For an actual
coupled residual H e=(r_nu,r_T), the exact equations instead give
 e_nu=-A^-1 C e_T-A^-1 r_nu,
 e_nu=(C*)^-1(D e_T-r_T).                                  (1)
The SECOND forcing term in (1) must be retained. The homogeneous output
map R_S T=((C*)^-1 D T,T) does not by itself reconstruct an application
error. In York coordinates the complete error is simply R(s,xi); it need
not have F or F^3 factors. No operator is commuted in this identity.

The first unsupported LOW-ORDER implication was a numerical output bound
for this unfiltered coupled error. Section 2 supplies it. The first
unsupported STRONG implication remains the source/graph/physical mapping
of the actual application residual. Sections 3-4 specify its complete
two-row content and test what the known estimates actually certify.

## 2. An unconditional physical bolt bound for a full coupled residual

Reuse the ordinary actual-domain endpoint proof from
current-complete-bolt-energy-bound.md, NOT just its receipt. It gives
 |s(z_h)|<=Cs E_s(s),
 |xi(z_h)|<=Kxi Cinf E_xi(xi),
 Cs=S0/sqrt(.99)+4(S2+2sqrt(.17)S1)/gamma,
 S0=1/sqrt(23.99), S1=sqrt(beta/5.99), S2=sqrt(beta/(6d)),
 Kxi=1+4/[3(alpha-2)], alpha=(1+sqrt(21))/2, Cinf=7.477e50.
The scalar bound uses the fourth-order scalar energy. The York trace
uses its original auxiliary Poisson profile with integrable endpoint
kernel and regular-bolt flux. Both extend by the original form closure.
They require neither a smooth quantum double nor any new quantum theorem.

Since the physical bolt functional is L_b R(s,xi)=ell_b[s+xi/2](z_h),
ell_b=(1+h)/h, Cauchy--Schwarz in the TWO positive energies gives
 |L_b e|<=Kapp E(e),
 Kapp=ell_b sqrt(Cs^2+(Kxi Cinf/2)^2).                      (2)
This includes arbitrary scalar residuals and every admitted energy
direction. It does not require r_nu=0, an F factor, or a homogeneous
solution. This is not a new bound on R_S V:E_xi->R, which is a different
map. All physical normalization and positive-temperature costs remain.

For a conforming trial u~, residual r=f-Hu~, uniqueness on the actual
domain and the full saddle estimate prove
 |L_b(u-u~)|<=Kapp ||r||E*.                                (3)
This is a complete two-row backward-error estimate for the bolt, NOT a
graph/C5 estimate or an evaluation of r. Its sufficient half-bolt-margin
target is ||r||E*<=(1.667e-9/2)/Kapp. Execute it.

## 3. Recover the undifferentiated source, including its missing moment

The original normalized action is W_action int F dz, W_action<0.
For canonical normalized residual densities F=(F_nu,F_T), set
 r=z^4 F, D_z=z partial_z, sigma=4.005, delta=sigma-4=.005,
 G=(D_z-4)r_T-(D_z-6)r_nu/3.                              (4)
The normalized energy-dual pairing is EXACTLY
 -(1/beta) int r dot R(s,xi) z^-4 dz.
Thus no action factor disappears in the following comparison.

First the prescribed affine source powers/logarithms must cancel exactly
between source and trial. An unmatched residual power below sigma makes
the original fast graph infinite, even if its coefficient is tiny. In
particular use the exact implicit slow root, not a rounded exponent.
When the residual satisfies this phase condition, write on the physical
half-line x=log z<log z_h
 r_nu=z^sigma f(x), G=z^sigma g(x), f,g in H3(dx).
These are the original scalar and differentiated shape rows. Define their
ACTUAL undifferentiated finite part
 c4=lim_(z->0) z^-4[r_T-r_nu/3].
It is an output of the original row, not an adjustable new source datum.
With K_delta h(x)=int_(-infinity)^x exp[-delta(x-y)]h(y)dy,
the exact recovery, including all transition and bolt values, is
 r_T=c4 z^4+z^sigma[f/3+K_delta(g-2f/3)].                   (5)
Indeed (D_z-4)(r_T-r_nu/3)=G-2r_nu/3. The integration constant in (5)
is mandatory. No choice at the upper endpoint is imposed.

Zero-extension is used only to prove L2 multiplier bounds for these
one-sided data; it changes no physical field or boundary condition.
Young/Plancherel gives ||K_delta||L2<=200 and
 ||(I-2K_delta)/3||L2<=133,
since sup_t |(delta-2+it)/(3(delta+it))|=133. Consequently, with
 H_mu=L2(z^-4 dz), M=[[1,0],[133,200]],
 Rfast=(||f||H3^2+||g||H3^2)^(1/2),
 ||r||_(H_mu^2)<=Krow Rfast+|c4|/sqrt(5), Krow=||M||2.     (6)
Here z_h<1 and sigma>3/2 justify the physical-weight embedding; the
quartic integral is z_h^5/5<=1/5. The same estimates hold for complex
fields and then for the real physical subspace.

For an explicit residual norm conversion, the actual metric-return bounds
give ||R(s,xi)||_(H_mu^2)<=MR E(s,xi), where
 MR^2=[16/23.99+1225/(.99 F_lambda)]/[4(3-sqrt(5))].
This follows from the complete tensor Gram matrix [[16,8],[8,8]], scalar
tensor norm 4||s|| and outer-fixed TT disk-return norm <=35||xi||. Hence
 ||r||E*<=MR/beta [Krow Rfast+|c4|/sqrt(5)].                (7)
Equations (3),(7) are a COMPLETE source-to-bolt bound for residuals whose
global original rows and actual finite part are enclosed. The original
localization cutoff is not moved; all its transition equations are inside
f,g. This global residual quantity is stronger than the recorded outer
observations alone. It is NOT asserted equal to the original Y norm.
In particular it does not assign I_Y, bound the infinite matching quotient,
or supply H3 output / C5 physical comparisons.

The missing moment is not optional: r_nu=0,r_T=z^4 has f=g=0 but c4=1
and ||r||H_mu=z_h^(5/2)/sqrt(5)>0. This is an actual smooth physical-side
canonical source, normalizable at infinity and regular at the bolt, with
the original paired symmetry. The complete form inverse applies to it.
This observation rejects certifying an original two-row residual using
only the differentiated pair. It does not prove that the actual quantum
tail has this particular residual. No independent moment is set to zero.

Execute the coefficients of Rfast and |c4| in (3),(7), and their individual
half-margin ceilings (each with the OTHER residual debit zero). They are
not simultaneous allocations and are not values of actual residuals.

## 4. Apply the map to the COMPLETE retained transport approximation

The actual already-defined full retained operator is H~=H-DeltaH, with
all low words at EVERY angular/thermal label, all original local contacts,
common finite restorations and finite returns. Its source is f~=f-Deltaf.
No independently renormalized truncated determinant is substituted.
This is its established FORM realization. Membership of H~^-1 f~ in the
original strong mass/log graph, and its numerical phase constants, have
not been certified. The bolt functional extends to the form completion,
so (8)-(9) do not require that missing assertion. A numerical trial used
in the strong residual test (4)-(7) must itself retain the original phase.
The established form comparison is
 epsilon=||DeltaH||E->E*<=1.020e-25,
 e_f=||Deltaf||E*<=1.365e-78,
 ||H~^-1||<=1/(1-epsilon), E(H^-1 f)<=1.1e-10.
For the exact mathematical responses u=H^-1f, u~=H~^-1 f~, subtract
BOTH equations before estimating:
 H~(u-u~)=Deltaf-DeltaH u,
 E(u-u~)<=(e_f+epsilon*1.1e-10)/(1-epsilon).                (8)
Thus the actual complete approximation has the unconditional bolt-error
upper certificate
 Uapp=Kapp (e_f+epsilon*1.1e-10)/(1-epsilon).               (9)
The sign of the source/Hessian terms is retained in the identity; their
norms are separated only in the final specified estimator. This includes
the source error omitted in the earlier illustrative 1.122e-35 expression.
Its source-only and Hessian-only debits must both be reported. These are
errors between defined full responses, NOT numerically executed inverse
applications or signed trial centers.

Execute (9) against half the original bolt allowance, and the separate
improvements needed for e_f and epsilon if each were the sole debit.
Increasing the Chebyshev degree does not decrease (8). P76's already
bounded exact-polynomial tail therefore does not certify this application.
Do not optimize the transport tail: the retained strong operator and
its original finite-part residual are still unevaluated, even if this
tail were known exactly.

All remaining outputs require their OWN backward-error maps. The known
actual-domain transition counterexample in
current-strong-quotient-reconstruction-test.md rules out a finite generic
energy-to-source-H3 comparison, even for Schur-harmonic quotient fields.
Our point trace does not repair that. Graph, independent phase and C5
comparisons remain unknown; no conditional ceiling is substituted for one.
The original P32 source requirement Cdef3 IY<=1.0965284986e32, graph and
physical ceilings with their remaining margins, actual kappa_m<=5e-5,
every signed coupled entry and the restricted quantum-fast screen remain
exactly as in the preceding checkpoint. Neither a graph defect nor the
mass target has been numerically certified here.

Outcome to test: a finite UNCONDITIONAL bolt-application comparison, an
excessive specified full approximation estimator, and an unresolved complete
strong application map with its precise source/boundary obligations. No
obstruction to the physical construction follows.

Verification contract: freeze argument/checker, verify the exact forcing
and row-recovery algebra and all original-unit coefficients at 512 bits,
retain outcomes, ordinary review, unchanged acceptance. Do not call a
checker a full inverse application or a proof of the missing strong maps.

Remaining to-do list: execute the complete original strong application
residual and signed coupled solve; nonlinear admission/contraction; exact
compatibility; physical matching; semiclassical validity/noise; all-sector
stability; opening, transport, preparation, separate creation and the full
stable semiclassical wormhole objective.
