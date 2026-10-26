# Complete remaining source: exact square return before taking norms

2026-10-26 HST, E003/S169. Ordinary argument and outcome contract BEFORE
arithmetic. Previous turn: progress (original weighted bolt output factor).
Preserve accepted 1ee884b, complete P32, the closed audit and the full
stable semiclassical wormhole objective. The COMPLETE estimate is unresolved.

## 1. Attempt and first unsupported step

The current target is M_rem=||W(Dbar-Dbar_ref)V F^3||_(E->H_mu),
H_mu=L2(z^-4 dz), with the ORIGINAL weight s=7/5 and .1/.2 transition.
E^2=ell/(|W_action| beta), F=(ell+B)^-1 B, B=b*A^-1b,
||F||E<=f=.338/1.338. These are the full operators on the accepted paired
mass/log realization, common regulator, fixed clock/temperature, and
regular bolt. The canonical quantum determinant, exact ultraviolet log,
all angular/thermal sums, contacts and finite returns stay in Dbar.
The physical row scaling is rho=8piG A_r^2 for BOTH C and D.

The complete source is not just the nonlocal determinant. It includes the
curvature square in I_full=I_ref-I_pot+I_sq+Gamma_q, where
I_sq=-c_R integral (R+12lambda)^2 dvol. Applying this fourth-order row
and then bounding Cbar*^-1 independently loses an exact cancellation.
The attempted two-row substitution
 Dbar V F^3=(V*)^-1 Bbar F^2-Cbar* Abar^-1 Cbar V F^3
does not by itself compute the required weighted derivatives. Nor does
the form comparison .0022 bound a strong source. No full M_rem value was
obtained by those substitutions. This attempt instead retains the square
and the mixed inverse in their actual order and tests the resulting WHOLE
output inequality below. No numerical quantum inverse is applied.

## 2. Exact square Hessian, including its off-shell density terms

Keep the affine canonical shape chart g2->exp(2tT)g2, with the physical
sphere radius fixed. Put w=R+12lambda, V0=2/r^2+12lambda and
P0=-Delta2+V0. The exact warped-curvature identity is
 R(t)+12lambda=exp(-2tT)(w-V0-2t Delta2 T)+V0.
The volume density is exp(2tT). Differentiating the PRODUCT twice gives

 (exp(2tT) w(t)^2)''|0
   =8[(P0-w)T]^2+4(2V0 w-w^2)T^2.                    (1)

The last term is mandatory on this nonstationary seed. It is not inferred
from the square's value or dropped because w is small. Polarization in
the original physical measure supplies the adjoint of P0-w.

Let t0bar=A_r^2(R+12lambda0), wbar=A_r^2 w,
V0bar=A_r^2 V0, Pbar=A_r^2(-Delta2+2/r^2+12lambda0),
J=Pbar-t0bar=A_r^2(P0-w), m=2V0bar wbar-wbar^2, and
dc=96piG c_R/A_r^2. Then the EXACT normalized shape square is

 Dbar_sq=-(2dc/3) J*J-(dc/3)m.                       (2)

Here dc differs from d=D_R/A_r^2 by G/(60pi A_r^2); they are not equated.
The fixed reference potential and all quantum/local restorations remain
in D_Q=Dbar-Dbar_ref-Dbar_sq. In particular D_Q is NOT a choice of stress
or a determinant renormalized independently from its contacts.

## 3. Actual-domain cancellation, before any numerical interpretation

Use the accepted mixed factorization in its original order:
 (Cbar*)^-1=(1/2)L_D^-1 (I+R*)^-1(Pbar*)^-1,
 L_D=1-d A_r^2 Delta4, ||W R* W^-1||<r_mix=1.556e-23.

For T=V F^3 u on the common core, J T is regular at the disk center and
has the original one-sided mass/slow expansion at infinity. Its mass
power is at least z^3; the weighted physical L2 norm is integrable.
The retained three-factor regularity gives local derivatives through
fourteen, so J T is in the actual normalizable Pbar* domain. Endpoint
Green products vanish for these powers and the original regular bolt.
No quartic/slow coefficient is prescribed away. Thus
 (Pbar*)^-1 J* J T=J T-(Pbar*)^-1 t0bar J T.          (3)
This is an inverse identity on its domain, not a commutation of factors.
The actual paired quantum domain is unchanged; the mixed auxiliaries
use their already established physical ODE domains. No smooth-double
Hessian theorem is invoked. Bounds extend only where their displayed
norms are finite. The original-Y extension is still a separate gap.

Put Z=(I+R*)^-1. Equations (2)-(3) give exactly

 (Cbar*)^-1 Dbar_sq T
  =-(dc/3)L_D^-1 J T
    +(dc/3)L_D^-1 Z[R*+(Pbar*)^-1 t0bar]J T
    -(dc/6)L_D^-1 Z(Pbar*)^-1 m T.                  (4)

The first term has differential order ZERO: L_D^-1 has order minus two
and J has order two. The fourth-order square has not survived as a
strong source requirement in that term. This genuine composition gain
does not remove its low-frequency contribution or evaluate it. The
independent quantum logarithm is retained in the combined output below.

## 4. Bound the COMPLETE mixed correction without an unknown A2

The initial frozen development verified the A2-dependent estimate below.
Before acceptance we sharpen its entire mixed correction, not the reference
or cold source. Preserve that development snapshot. Both weight signs
are essential for taking the following adjoints.

Write mD=1-d s^2>.99999996536. For the conjugated L_D let x=L_D,W v.
Its real form gives Re<x,v>>=mD ||v||^2+d||grad v||^2. Therefore
 ||x-v||^2<=||x||^2+(1-2mD)||v||^2<=||x||^2.
Also maximizing ||x||||v||-mD||v||^2 gives
 ||W A_r grad L_D^-1 f||<=gD||Wf||,
 gD=1/(2 sqrt(d mD))+s/mD.
Thus ||W(I-L_D^-1)W^-1||<=1; the two resolvents are not bounded
separately. This is an elementary coercivity calculation on the actual
closed domain. The same calculation holds for W^-1.

Use -A_r^2 Delta4 L_D^-1=(I-L_D^-1)/d and the actual adjoint expression
 J*=-A_r^2 Delta4-2A_r^2 grad(log r).grad
                      +Vbar-2A_r^2 Delta4 log r-t0bar.
The retained |A_r grad log r|<=1 and |A_r^2 Delta4 log r|<=27 beta give

 KJL=||W J* L_D^-1 W^-1||
   <=1/d+2gD+(2.001+54beta+tmax)/mD.                 (4a)

The earlier mixed proof gives, for u=Pbar^-1 f and BOTH weights,
 ||Wu||<=26||Wf||, ||W A_r grad u||<=66||Wf||,
 ||W A_r^2 Hess u||<=950||Wf||.
Put t1=1.4e-26 beta^(3/2), t2=8.6e-25 beta^2. These follow from the
actual gradient/Hessian bounds on w; the constant vacuum shift has zero
derivatives. In four dimensions |Delta u|<=2|Hess u|. Consequently

 Dt=1900*tmax+132*t1+26*t2
   bounds ||W A_r^2 Delta4(t0bar Pbar^-1 f)||/||Wf||,
 EJ=2d Dt+[G/(60pi A_r^2)]1e12
                   +40beta(4e-27)*26/(60pi)
   bounds ||W Ebar Pbar^-1 W^-1||.                  (4b)

Ebar is the ENTIRE mixed correction: the first term includes derivatives
of t0, the second the complete curvature anomaly, the third its magnetic
density contact. None is a nonlocal shape approximation. Since
R=(1/2)L_D^-1 Ebar Pbar^-1, taking the adjoint of the opposite-weight
estimate yields ||W R* J W^-1||<=KJL EJ/2.

Similarly J* t0bar Pbar^-1 is bounded by the product rule:
 KJP=1900+132+26(2.001+54beta+tmax),
 Kt=tmax KJP+132*t1+26(t2+2t1).
Its opposite-weight adjoint is (Pbar*)^-1 t0bar J. Thus

 ||W[R*+(Pbar*)^-1 t0bar]J W^-1||<=Ksum,
 Ksum=KJL EJ/2+Kt.                                 (4c)

These closed composed maps carry their original domains and finite
returns. No J is commuted through R or a multiplier. Formula (4c)
removes the A2 input from the ENTIRE second term of (4), at a computed
cost Ksum. It does not bound the main quantum map or the actual correction.

## 5. Complete original bolt budget and known/unknown factors

Define the JOINT main output, keeping its signed sum before the norm,

 J_b=||ell_b ev_b[-(dc/3)L_D^-1 J+(Cbar*)^-1 D_Q]V F^3||_(E->R).

It includes the square's leading return and ALL the remaining quantum,
potential, contact and boundary returns. Its actual norm is UNKNOWN.
The initial A2=||W J V F^3||_(E->H_mu) remains UNKNOWN, but is no longer
needed to bound the mixed correction in (4). The previously proved
||Vv/z^2||2<17.6 E(v) and W<=h^-7/5 do give the finite, conservative

 A0=||W V F^3||_(E->H_mu)<=17.6 f^3 h^-7/5.          (5)

This full-weight cost is paid explicitly, not improved by normalization.
The existing point comparison for L_D^-1 is 2e6/sqrt(gamma), and the
Pbar* inverse in weighted L2 is <26. Put
 Lphys=ell_b (2.5h)^(7/5)/sqrt(gamma),
 wmax=1.1e-27 beta,
 tmax=beta[1.1e-27+12*11 G lambda/(360pi)],
 mmax=2(2.001)wmax+wmax^2.
The true signs of w and m need not be assumed. The exact inequalities
V0bar<2.001 and |t0bar|<=tmax use the retained vacuum relation.
The initial WHOLE remainder estimate in (4) had original bolt norm at most

 R2 A2+R0 A0,
 R2=Lphys*(2e6)*dc/3*(r_mix+26*tmax)/(1-r_mix),
 R0=Lphys*(2e6)*dc/6*26*mmax/(1-r_mix).              (6)

The sharper executed bound from (4c) replaces R2 A2 by Rmix0 A0, where

 Rmix0=Lphys*(2e6)*dc/3*Ksum/(1-r_mix).               (6a)

All compositions use the same weight; no factor is moved through a
variable coefficient. The error bounds contain the actual complete
mixed remainder, not only its principal symbol.

Consequently the complete original bolt test is

 p_aff[J_b+(Rmix0+R0)A0+C_b3_ref]<=B_remaining,
 p_aff=(1.1e-10)h3.                                 (7)

Here C_b3_ref and h3 are the unchanged accepted P32 quantities. Execute
(4a)-(7), freezing p_aff R0 A0<1e-18 and the COMPLETE mixed-correction
debit p_aff Rmix0 A0<1e-9. The initial full-margin A2 and J_b sensitivity
slices remain in the receipt as historical comparisons, not simultaneous
allocations. The revised J_b allowance subtracts the computed mixed debit
FIRST. No actual A2, J_b, source norm or signed center is substituted by
its allowance. Other trial/error debits must still be deducted.

For orientation also evaluate the ceiling for the UNMULTIPLIED filtered
square point map ell_b ev_b L_D^-1 J V F^3: it is three times the joint
remaining ceiling divided by dc, only if all other main contributions
are zero. This is not that map's norm and is not an allocation adopted
in (7), which retains the combined signed main map.

The first unsupported quantitative step is now explicit: bound the joint
map J_b in (7) on the actual form inputs. Its entire mixed/off-shell
correction has a numerical bound without the uncomputed A2.
Qualitative regularity of F^3 does not provide those constants. The
10^4000 C interior and 10^74739.22 generic absorption coefficient are
absent from (4)-(7); neither has been replaced by an assumed small bound.
The retained M_rem source norm itself remains unevaluated, not refuted.
No claim that its existing sufficient target fails is made.

## 6. Decision and scope

If the exact algebra and both whole-remainder debits pass, retain the
construction and P32 and SHARPEN THE OUTPUT COMPOSITION using (4),(7).
The reason is an actual two-derivative cancellation on the retained
domain, not an isolated contact improvement or a smaller normalization.
The square and quantum terms must remain together in J_b. This attempt
does NOT establish that the joint target is realistically attainable.
It supplies an exact alternative and an evaluated complete-budget
interface, with the decisive norm still missing. Further optimization
of R0, R2, the cold source or the already bounded reference row is not
the next task.

The original defect still requires C_def,3 I_Y<=1.0965284986e32 and a
bounded original-Y extension. Graph/physical ceilings 1.0650980044e31
and 9.9684408960e33 retain their margin hypotheses. Every coupled entry
requires |c_hat_ij|+e_ij+C_i3 I_j h2. Actual kappa_m<=5e-5 is untested;
its conditional costs 2.14765e-12, .0024731 and .0195239 do not establish
the remaining entries. The original joint affine/feedback budgets,
mass conversion 2133--2174, scalar quartic, correlated slow pair and
infinite matching quotient remain. The .135973<.1599164 quantum-fast
screen retains only its admitted two-fast scope. Neither an independent
phase nor the nonlinear neighborhood inherits that screen.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): enclose the JOINT complete
main output J_b against the remaining allowance in (7), using
the actual three Schur factors and common-restored quantum operator.
That calculation must also supply the actual remaining margin; another
local square/contact calculation is not its execution.

Verification contract: freeze this ordinary argument and exact checker;
verify (1)-(4) symbolically, the square coefficient against the original
local action, actual interval guards and all displayed physical debits
at 512 bits. Retain failures, review, then run unchanged acceptance.
Supporting checks do not apply F numerically or prove its missing norms.

Remaining to-do list: complete joint output/source maps and every coupled
entry; nonlinear admission/contraction; exact compatibility; physical
matching; semiclassical validity/noise; all-sector stability; opening,
transport, preparation, separate creation and the full stable semiclassical
wormhole objective.
