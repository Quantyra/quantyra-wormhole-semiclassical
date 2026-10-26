# Signed strong local core for the full coupled operator

2026-10-25 HST, E003/S169. Ordinary argument before frozen execution.
The full stable-wormhole objective is active and incomplete. The closed
feasibility audit and original physical norms remain unchanged.

This is implementation of both differential rows needed by the complete
signed solve, not a completed coupled bound. It does not optimize a local
source or substitute a local matrix for the full operator. The canonical
quantum block still needs its actual strong application. No kappa_m,
nonlinear admission, contraction or physical solution is certified here.

## 1. Original action, coordinates and exact strong operator

Use the affine canonical chart of current-finite-local-action-assembly.md:
phi=-log(z)+nu, shape T, Delta v=(Q v')', D=z partial_z, and the ORIGINAL
dz action measure. The fixed factor is W=-N pi B^2/(G gamma), with G=1
as before. Use the exact finite-seed Q and original paired mass/log and
regular-disk domains. There is no new inner wall or chosen mass condition.

The implementation differentiates precisely that note's Einstein/Maxwell,
specified R2 and four-dimensional Weyl densities. Each depends on the
two fields and their Euler jets through order two. Its second weak
variation has the exact form

 sum_(a,b=nu,T) sum_(i,j=0)^2 int
                  (D^i v_a) C_(ai,bj)(z) (D^j w_b) dz,
 C_(ai,bj)=partial^2 L/[partial(D^i u_a)partial(D^j u_b)] at u=0.

Since D*=-D-1 in dz, the ORIGINAL strong density Hessian is

 H_ab=sum_(i,j=0)^2 (-D-1)^i C_(ai,bj) D^j
     =sum_(k=0)^4 H_ab,k(z) D^k,
 H_ab,j+l += (-1)^l binom(i,l)(-D-1)^(i-l) C_(ai,bj).

This is an identity on arbitrary direction jets. It does not select a
metric trial space or physical quantum modes. Multiplication by z^4
returns A=z^4 H_c, the original operator in the graph argument. Neither
the physical factor nor the original source/admission norms are changed.
All cubic and quartic variable coefficients remain in this signed core;
none is estimated by the failed independent-positive matrix.

For a coefficient depending on z and the finite seed jets, exact total
differentiation is

 D a=z partial_z a+sum_(j>=0) z Q^(j+1) partial_(Q^j) a.

The generated differential expressions need Q through order four. The
interface also retains the higher exact seed jets needed to differentiate
coefficients. It raises an error if the declared jet order is exceeded.
The actual Q is never replaced by its extremal or product approximation.

## 2. Endpoints and the retained Weyl boundary representative

Before imposing physical endpoint data the first-variation boundary
expression for an arbitrary local density is EXACTLY

 [ z{(L_u1-(D+1)L_u2) v+L_u2 Dv} ]_lower^upper.

The second-variation expression is obtained by varying every L_ui and
retaining both independent direction jets. These expressions are outputs
of the generator; they are not silently set to zero. On the original
mass/fast tangent class the retained local Hessian has the endpoint
cancellations proved in the finite-local-action note. At the regular
bolt Q(z_h)=0 must be used as an exact correlated root identity. Merely
evaluating the degree-six polynomial with independently rounded bolt data
would not justify that cancellation.

The reduced Weyl density differs from the unreduced action by

 W_endpoint= C Q/(60 pi B^2 z)
       {phi D[exp(-2T)L]-(D phi)exp(-2T)L},
 L=-Q''-2Delta T.

This endpoint representative is retained explicitly. Its first variation
is not declared zero. The original common boundary subtraction and
completion are still required for physical source pairings. Its second
variation vanishes on the original fixed-source mass/fast tangents, as
already proved. A bulk strong source is not the full boundary functional.

## 3. Independent trace-row identity and reference mass cancellation

The scalar and mixed rows of the complete physical action are local:
the canonical determinant has no nu slot in this affine chart. The
independent Euler trace identities in current-finite-local-action-assembly.md
therefore check the ENTIRE first Hessian row, allocation by allocation,
for arbitrary formal field jets and arbitrary physical-side seed jets.
This checks the Euclidean sign, density contacts and fourth derivatives,
not merely the symmetry of a sampled matrix.

For a separate comparison, the reference density is

 F_lambda{ exp(2nu)Dcurv/(2z^2)+exp(2nu+2T)/(Cz^2)
     +3beta exp(4nu+2T)/z^4-(1/C+3beta+eps)exp(2T)
     +4eps z exp(-nu+2T)-3eps z^2 exp(-2nu+2T) }.

Here Dcurv=L-6Delta phi-6Q phi'^2. This reference is NOT added to the
physical local sum. Its two bulk seed rows vanish for

 Q=beta+z^2/C-(4beta+2/C+tau)z^3
                         +(3beta+1/C+eps)z^4-2eps z^5+eps z^6.

The exact uncut mass tangent has J'=z^3/Q^2,
n_m=-2beta QJ/z, b_m=beta z^3/Q+beta Q'J. Its areal return is
f=2xz, alpha=0, so (rf)'=alpha'=0 and both reference equations cancel.
The implementation must reproduce this identity with BOTH strong rows.
The arbitrary original mass cutoff remains explicit: derivatives of
chi_m multiply these exact profiles by the ordinary product rule. For
chi_m=1 through .1 and zero by .2, retain

 f_cut=2xz chi_m-2F n_m Dchi_m, alpha_cut=n_m Dchi_m.

The cut lift vanishes near the regular bolt, has zero reference bolt
shift, and its entire transition carries reference mass -xB. It is not
a global fixed-temperature zero mode. No numerical cutoff derivative
constant is assigned to the previously unspecified profile.

## 4. Complete assembly obligation and quantitative role

The physical Hessian in this chart is exactly

 H_full=H_Einstein+H_Maxwell+H_R2+H_Weyl
                                  +diag(0,Gamma_hat''[T]).

Gamma_hat includes the separately declared LLL and the complete massive
canonical determinant with its common regulator, all angular/thermal
modes, coincidence contacts, paired boundary and finite-domain returns.
The angular-resummed action identity supplies a representation, not an
evaluated strong application. It must replace its entire previous massive
allocation once; it must not duplicate the older local restorations.

The generator's local differential representation is exact, so it adds
zero symbolic local truncation error to a future full residual. Evaluation
roundoff, interpolation of a candidate field, the quantum application,
physical boundary equations and any inverse estimate still need their
own certificates. In particular local exactness does not reduce either
global budget without the complete composed inverse and output maps.

Original mass conversion 2133<||d_chi||H3<2174 and the known LINEAR
|m|<1.627532e-15 retain debit <=3.538254568e-12. At CONDITIONAL
kappa_m<=5e-5, all remaining entries must still obey

 Y_other+4.084e-22 A_o+Y_phase<=5.34235e-12,
 eta_other+.0195239 gamma_E+3.712e-12 A_o J_E+eta_phase<.0012069.

No actual value of these unknown maps is supplied. The fast screen
.135973<.1599164 retains its two-fast/admitted-family scope only.
Physical-X/neighborhood comparisons and all original nonlinear gates
remain unproved. Exact compatibility, physical matching, validity/noise
and all-sector stability remain open.

## 5. Attempted complete return and the decision boundary

The attempted reduction starts with BOTH original rows, in the same chart:

 H_nunu nu+H_nuT T=f_nu,
 H_Tnu nu+(H_TT,local+Gamma_hat'')T=f_T,

together with the original phase/interface equations and regular bolt.
Even if the already retained weak scalar inverse is used, its exact
Schur expression is

 S_T=H_TT,local+Gamma_hat''-H_Tnu H_nunu^-1 H_nuT.

Thus eliminating the scalar does NOT remove or estimate Gamma_hat''.
The local generator accounts for the entire differential part of these
equations; it does not give a strong inverse for H_nunu or S_T. Nor does
reference annihilation of the uncut mass profile imply annihilation by
this physical operator. The original cutoff adds the explicitly nonzero
two-row commutator sum H_ab,k binom(k,j)(D^j chi_m)D^(k-j)Phi_m,b,
j=1,...,k. Its numerical derivative constants are not specified in the
retained choice of chi_m. They cannot be assigned zero or guessed from
the different .05/.1 graph cutoff.

For the full commonly regulated massive term, insertion of the previous
all-angular representation gives the metric Hessian of ALL terms of
equation (9) in current-angular-resummed-quantum-operator.md. The exact
weak Duhamel formula has both the bubble and the second metric vertex.
The required next operation is to regard that combined functional as a
source in the original near/transition graph norm and apply the original
border. This is the FIRST UNSUPPORTED CONTINUUM STEP. Reinspection of
the actual-domain proof gives an admitted-family C2 scalar functional
and bounds with BOTH inputs in the admitted phase norm. It does not
give a bound with the testing input in the dual of the desired strong
source norm. In particular integration by parts in a localized HEAT TRACE
is not an estimate for spatial derivatives of its variational gradient.
The positive angular bubble bound between unbounded metric vertices
does not supply that estimate either.

No theorem for a smooth closed metric, a fixed square-operator domain,
or a finite mode selection is used to cross this gap. The first-order
paired domain, mass/log weighted C5 hypotheses, common regulator,
completed coincidence contacts, all thermal images, actual finite disk
and gluing defect stay as in the retained weak argument. These are the
proved scope of that argument, not a strong source theorem. No source
field or coefficient for Gamma_hat'' is numerically interpreted here.

Consequently this attempt does NOT execute a complete coupled numerical
return. Its executed part, once verified, is the arbitrary-function
two-row differential application, its reference mass cancellation and
its nonzero cutoff commutator. The combined quantum application and
conforming strong solve remain unexecuted. This is an UNRESOLVED
estimate, not a verified excessive bound for the actual coupled column,
and not an obstruction to this construction. It closes NO additional
admission, feedback or physical-solution inequality. The requested
complete coupled feasibility milestone remains open.

For clarity, the residual budgets give the individual ceilings
A_o<=1.30811704211557e10, A_o J_E<3.25134698275862e8 and
gamma_E<.061816542801387 ONLY when every other debit in the respective
inequality is zero. These are conditional allowances, not estimates of
those maps, and they cannot all be independently spent. The previously
proved linear energy upper bound E(h_lin)<1.1e-10 alone also fails to
certify the original neighborhood requirement E_center<=1e-11, even
before any positive-radius debit. It is an insufficient known upper
bound, not evidence that the actual center is outside that neighborhood.
The required original nonlinear inequalities remain

 E_center+1e-8 alpha_embed<=1e-11,
 R_X alpha_embed*.028913+Z0_other<=.25,
 R_X alpha_embed^2*1.0706e10+Z2_other<=5e7,
 Y0+1e-8 Z0+5e-17 Z2<=1e-8, Z0+1e-8 Z2<1.

The physical-X comparison constants and independent boundary-phase
columns remain unknown. A change of boundary basis changes none of
these physical requirements.

RESEARCH DECISION: retain the construction; sharpen the method at the
combined quantum STRONG APPLICATION, with the signed differential core
now explicit. Do not extend the estimate chain or optimize controlled
cold/scalar terms. The failure to obtain a strong map from weak trace
control is the reason for this decision, not a large actual mass return.

ONE NEXT CALCULATION, PROPOSED NOT EXECUTED: evaluate the full regulated
geometric Hessian as an original strong source operator on a conforming
full coupled parametrix P, retaining both original equations and its
phase/interface border. Establish regulator convergence in THAT source
norm and bound delta=||I-S P||<1 (sufficient target <=1/2) and the original
output maps M_i=||T_i P||. The common calculation must include mass,
quartic, correlated slow and fast columns; its column error is then
M_i||(I-S P)b||/(1-delta). This supplies an actual kappa_m and debits
only after those domain and inverse assertions are proved. It is not
authorized by the executed local algebra or angular scalar weights.

Remaining to-do list: complete actual quantum strong application and
coupled parametrix/residual; original norm comparisons; nonlinear
admission/contraction; compatibility, physical matching, validity/noise,
all-sector stability and the full stable semiclassical wormhole objective.
