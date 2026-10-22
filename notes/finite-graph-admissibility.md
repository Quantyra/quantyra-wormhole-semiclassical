# Actual source regularity and admissible positive graph transfers

2026-10-21, E003/S167. Ordinary argument before supporting verification.
Keep the original finite comparison, full physical Hessian, quantum
state/domain, common subtraction, mass column and fixed return.

The actual comparison residual and its Schur source have every finite
weighted Sobolev order. The actual positive solves and their transfers
therefore belong to the original fourth-order certificate domain.
This validates the exact transferred-moment identity on these sources.
The higher norms are finite but unevaluated; the numerical cancellation
target, signed zero and nonlinear solution remain unproved.

## 1. Precise statement and unchanged operators

Use the [original comparison regularity](finite-comparison-conormal-regularity.md),
[strong Schur source](finite-strong-schur-source.md) and
[weighted positive inverse](finite-positive-weighted-inverse.md).
Write z=B/r_seed, W=z^-2, lambda=lambda_e, g=G lambda and
gradbar=lambda^-1/2 grad for the physical finite metric g_ell.
For a radial scalar or tensor define

    ||v||_(H_W^k)=sum_(j=0)^k ||W gradbar^j v||_L2(dmu_ell),
    H_W^infinity=intersection_(k>=0) H_W^k.

Use smooth radial disk fields at the bolt, not arbitrary interval data.
The outer norms are equivalent to z^(7/2) H_b^k with measure dz/z;
unweighted physical H^k corresponds to the exponent 3/2. Derivatives of
log W are bounded at every finite order, so ||Wv||_H^k is equivalent
to this norm. All constants below are finite for the fixed benchmark
and positive temperature. Higher numerical values are not asserted.

The actual residual is f_ell=Jvol F_ell(g_ell) for the bookkeeping
action I_ell=I0+ell Gamma_rest. This action is distinct from the FULL
Hessian used in the following operators unless ell=1. Set

    (f_tau,f_c)=lambda^-2 Rcal_ell* (Jvol^-1 f_ell),
    p=80pi^2(f_tau-Bbar Dbar^-1 f_c),
    T=-Box_ell/lambda, H=T+4, D0=T(T+4),
    E=60pi/g+11/3,
    Q=(80pi^2/E) Bbar Dbar^-1 Bbar*,
    K=-divbar((1+alpha)gradbar)+4+b+Q.              (1)

Here Dbar is the established positive critical TT realization, Bbar
is its full off-shell mixed row, and alpha,b are the unchanged local
coefficients of the positive-inverse proof. In particular |alpha|<1e-32.

For every ell in [0,1]:

    f_ell in (H_W^infinity)^2,
    f_tau in H_W^infinity, f_c in H^infinity,
    p in H_W^infinity.                              (2)

Starting with q0=p, the ACTUAL unmodified positive inverse gives

    u_j=K^-1 q_j in H_W^infinity intersect D_cert,
    q_(j+1)=D0 u_j in H_W^infinity, j>=0.           (3)

Only finitely many orders are used for any fixed iteration count.
This is a regularity and domain result; it supplies no small high-order
constant or useful bound on N(q3). At ell=0 the residual and all these
iterates are zero, as in the preceding source theorem.

## 2. Differentiating the existing weighted kernel argument

Let A_ell be the FULL affine disk Hessian in (S,Tdisk) coordinates.
Use X_delta^m,Y_delta and the one mass vector e_m=chi z^3(1,-2)
from [the actual weighted map](stationary-global-weighted-map.md).
For integer k>=0, X_delta^(k) has k extra b-Sobolev derivatives in
both columns, retaining the shape logarithm; Y_delta^(k)=z^delta H_b^k.
Append the SAME e_m to define X_delta^(k),m. At the bolt use disk norms.
We now prove, for each finite k and -1<delta<4,

    A_ell:X_delta^(k),m -> (Y_delta^(k))^2.          (4)

The earlier proof supplied k=0, finite diagonal symbol estimates and
separated-kernel bounds with two powers of actual-minus-normal gain.
The new input is that all rescaled coefficient jets needed for a chosen
k exist, uniformly along the original finite path. This follows from
the preceding conormal theorem, not from a general C5 hypothesis.

Here is why nonsmooth compactification does not prevent this extension.
On z,y comparable to a dyadic scale h, first take interior balls of
radius proportional to h, staying a fixed fraction of h from z=0.
All coefficients there are smooth. The scaled derivatives of the
actual-minus-flat coefficients have bounds C_j h^2(1+|log h|)^N_j
for every specified finite j. A z^4 log z coefficient obeys these
bounds despite lacking an ordinary fourth derivative at z=0.

Near the direct diagonal, construct the interior Dirac parametrix to
the additional finite order required by k. Its coefficients and
remainders use only those scaled interior jets. The boundary/global
remainder solves the homogeneous elliptic equation away from the
diagonal. Interior elliptic estimates on the same h-balls upgrade the
existing zeroth/low derivative remainder bounds to any specified
scaled observation derivative. Apply the equation in each variable
for joint derivatives. Subtract enough direct parametrix terms first
when the observation points coincide. This also treats the finite
Green contacts: their inhomogeneous local remainders are explicitly
determined by the same finite jets.

Thus no smooth-boundary theorem is applied at infinite ordinary order
to z^4 log z. The established boundary realization and low-order
boundary estimates are retained; higher differentiation is performed
on interior scaled balls. On separated scales the same argument in
each observation variable preserves the existing power of z+y.
Separated compact-to-end terms retain their physical z^4 times finite
logarithms. All global regular parts, thermal terms and contacts remain.

In x=log z, let B=partial_x on the outer collar, smoothly cut off
before the bolt. For a kernel written against dy/y the commutator is

    k_comm(z,y)=(z partial_z+y partial_y) k_A(z,y),   (5)

with the usual cutoff terms on a compact annulus. Formula (5) is a
distribution identity, including the local finite-part contacts.
The complete normal operator is translation invariant in x: the flat
logarithmic scaling term and its Weyl contact must cancel before this
assertion. Repeated commutators of A-N retain the two-power gain.
In particular they obey the same separated weighted Schur strip

    -1<delta<4,
    exp[-(4-delta)|x-x'|], exp[-(1+delta)|x-x'|].    (6)

Every finite number of log factors remains integrable against these
strict exponential margins. Near the diagonal, extract the unchanged
fourth-log shape symbol before estimating the lower-order finite part.
Normalize by <p>^4 log(2+|p|) in that column and by <p>^4 in the local
column. The differentiated symbols have the same finite L2 seminorm
bounds, with new constants, because the required extra coefficient
jets were just established. One may use [Grubb, Theorem 2.3](https://web.math.ku.dk/~grubb/CPDE14.pdf)
at output L2 for each differentiated normalized symbol. Its hypotheses
involve finitely many frequency derivatives and spatial Holder control.
We obtain higher orders by commutation, not by applying an unrestricted
Sobolev theorem to coefficients known only to be C5.

The smooth local reflected model is consistent with
[Gaarde--Grubb, Theorem 2.4 and Proposition 2.9](https://arxiv.org/pdf/math/0703878).
Those results do not themselves assert (4) for the present logarithmic
compactification; the scaled interior/remainder argument above supplies
the needed extension. No general polyhomogeneous heat theorem is assumed.

For compact inputs, Leibniz and (5) now give

    B^k A v=sum_(j=0)^k binom(k,j) ad_B^j(A) B^(k-j)v.

Use the just-established k=0 bounds on every summand, then complete.
Compact cutoffs are dense in the fast remainder norms. For the mass
column retain N e_m=0; it is not estimated as a generic z^3 source.
The earlier mass-column calculation differentiates by the same scaled
estimates: comparable scales give z^5 times finite logarithms, while
the far integral z^4 integral_(2z)^epsilon y^-3 y^3 log^N(y)dy and
separated terms give z^4 times finite logarithms. Every b derivative
has the same admissible rate. This proves (4), including the mass column.

Constants in (4) can be uniform in ell on [0,1]: the finite required
jets are uniformly bounded there, the paired realization has its
retained positive gap, and the fixed compact interior/bolt estimates
are uniform. No uniformity as temperature tends to zero is inferred.

## 3. The actual quantum source along the finite path

Do not infer source differentiability from (4) alone. Identify it first
on compact test variations by the actual determinant formula. Along
the original path the orbit Dirac operator in fixed half-density
coordinates is a differentiable family with the original complementary
domains. The path variations and their first normal derivatives vanish
at the outer boundary. Its low coefficient norms are bounded uniformly;
the positive square therefore has the common form domain, differentiable
resolvent, and the retained spectral gap. For compact test variations,
differentiate the resolvent/heat formula with its two vertices AND
second contact vertex. At short time the test collar is smooth and the
common four-dimensional subtraction controls the differentiated integral;
at long time the fixed gap controls it. The full angular trace is summed
before subtraction. Boundary-supported local terms have zero action on
these compact tests, and the existing boundary prescription is retained.

This yields the distributional derivative of the actual bulk source
with respect to the path parameter. It is not a new Banach-space C1
claim on an arbitrary neighborhood of metrics.

Use the affine disk coordinates of
[the full operator](general-background-newton-operator.md),
Phi=log(r/B), Vdisk=the product-orbit conformal factor. For
t=r_seed/B, tau=(t-1)/t, k_ell=exp(ell V), T_ell=1+tau(k_ell-1),
the original path tangent is exactly

    hdot_ell=(S_ell,Tdisk_ell)
            =(tau k_ell V/T_ell, u).              (7)

The fixed disk coordinate changes contribute no ell derivative.
By the preceding metric theorem, hdot_ell is a fixed correlated mass
cubic plus an all-conormal order-four logarithmic remainder. It is
therefore in X_delta^(k),m for every finite k and delta<4. Its norms
are bounded along this one path. In particular its leading coefficients
obey 4S_3+2Tdisk_3=0. The positive fixed-volume ratio

    Jvol=T_ell^4 exp(2ell u)

and its inverse have bounded b derivatives of all finite orders.

Let G_rest(ell) be the gradient of Gamma_rest at g_ell per FIXED
seed volume. The compact-test differentiation gives exactly

    d_ell G_rest(ell)=Jvol A_rest,ell hdot_ell,
    A_rest,ell=A_full,ell-A_local,ell.              (8)

The full density contact is included by this fixed-volume derivative.
The local action I0 has the same higher-order mapping and mass-column
property by direct differentiation of its finite-order local formulas:
the fixed vacuum jets cancel, and the mass cubic is the asymptotic
mass variation with unchanged curvature/clock. Thus (4) applies to
A_rest,ell as well. Reference-potential and Weyl contributions have
not been dropped or counted twice.

The right side of (8) is uniformly bounded in each Y_(7/2)^(k).
Its distributional continuity and separability give strong measurability
there (pair first with the dense compact tests). Hence it has a Bochner
integral in each such Hilbert space. The fixed-seed source belongs to
all of them: use the actual seed source expansion, massive/Weyl
cancellation and exact Ward reconstruction from the preceding theorem.
Integrating (8), and using its equality on compact tests to identify
the result, proves the same higher weighted regularity for the ACTUAL
finite-metric quantum source.

Add the explicit local I0 gradient and multiply the quantum gradient
by ell as specified in I_ell. This proves the first assertion in (2).
The argument does not require, or assert, a complete pointwise
polyhomogeneous expansion of the finite-metric quantum source.
It also does not supply numerical higher norms or their ell^2 constants.
The already proved low-norm ell^2 bounds remain unchanged.

## 4. Critical TT inverse and the mixed row at higher orders

Keep the actual scalar-c realization and fixed return from
[the finite TT theorem](finite-comparison-tt-realization-and-range.md).
Its York equation, return equation and their physical transposes now
map every finite Sobolev order: differentiate their regular elliptic
equations on the disk and the uniformly controlled outer annuli.
The strict strip (6) controls the unchanged tails. No homogeneous return
mode is added. The exact cancellation

    d(K_Y Tdisk)=(2/3)Tdisk dlog r-(1/3)Ric(Y,.),
    N Tdisk=-2Tdisk/3+2K_Y Tdisk-2RY

therefore gives N+(2/3)I of order minus one at every finite order.
The inverse return L has the analogous gain L+(3/2)I. In particular
N and L preserve the fourth-log spaces: their improving remainder maps
H^4 into H^5, which embeds in H^(4,log). This identifies the ORIGINAL
scalar-c domain with the corresponding fourth-log domain in this
more regular family; it does not replace the old realization.

On the critical physical line delta=3/2 the existing TT inverse is
positive, unique and bounded. Its scalar normal family is invertible
on that line, and its interior principal symbol is elliptic of order
four with one logarithm. Apply the same finite differentiated interior
parametrix and tail-small normal construction used for (4). For each
integer s>=0 this gives

    ||Dbar^-1 f||_(H^(s+4,log)) <= C_s ||f||_H^s.  (9)

To see that this estimates the existing inverse, start with its L2
solution. The elliptic estimate has an additional C_s||c||_L2 term;
the already established critical inverse bounds that term. Smooth
input approximation and closedness identify the limits with that
same solution. Equivalently use b-translation difference quotients
on the outer end: the differentiated commutators have the bounds in
section 2, while the critical normal line has no mass root. Interior
disk regularity treats the other end. This is a regularity extension,
not an application of positivity alone or a change of outer data.

The mixed row remains exactly

    Bmix c=-(8piG)^-1[Ttrace h_c-div(Etr eta_c)].    (10)

It has at most two derivatives of h_c, plus its fixed return. The
coefficients WJ, WC4, WZ and WEtr, and every specified physical
derivative needed in (10), are bounded. This follows from their actual
local curvature expressions and the now differentiated conormal
metric expansions, including the correlated mass cubic. It is not
obtained by differentiating their previous unsigned numerical bounds.
Their weighted values were already bounded in the strong-source proof;
their known expansions preserve those powers under b differentiation.
Leibniz and the higher return maps therefore prove

    Bbar:H^(s+2) -> H_W^s, s>=0,
    Bbar*:H^s -> H^(s-2), s>=2.                   (11)

All pairings and transposes are the original physical scalar-c ones.
The Weyl-invariant determinant is still in Dbar; it does not add a
second conformal trace to Bbar.

The York transpose in (1) preserves ordinary finite Sobolev orders.
Its conformal column is the original scalar conformal direction, so
f_tau retains the weighted regularity; f_c at least retains unweighted
regularity at every order. Combining this with (9)-(11) proves the
Schur-source assertion p in H_W^infinity. It also proves

    Q:H_W^s -> H_W^s, s>=2,                       (12)

by the sequence H_W^s -> H^s -> H^(s-2) -> H^(s+2,log) -> H_W^s.
For s=0 use the existing strong bound Q=AA*, ||WA||<1500 and W>=1;
no negative-Sobolev extension of Dbar^-1 is required in this proof.

## 5. Positive inverse bootstrap and the certificate core

The positive-inverse theorem gives u=K^-1q in H_W^2 for every q in
H_W^0. Keep its numerical constants and original realization.
The local part Lloc=K-Q is uniformly elliptic, with smooth bounded
physical coefficient derivatives at every finite order in this family.
For k>=0 the radial annular and regular-disk elliptic estimate is

    ||v||_(H_W^(k+2))
       <= C_k (||Lloc v||_(H_W^k)+||v||_(H_W^0)).  (13)

Conjugate by W to obtain a uniformly elliptic scalar expression whose
extra coefficients involve bounded derivatives of log W. On the outer
end use x=log z, freeze its positive leading coefficient on uniform
small intervals, absorb that principal error, and differentiate the
remaining bounded coefficients. Patch with bounded overlap. On the
compact bolt use the smooth two-dimensional disk elliptic estimate.
This proves (13) also for distributional solutions by approximation.
There is no ordinary finite-boundary extension hidden in this step.

For q in H_W^2, u is already in H_W^2, so (12) puts Qu in H_W^2.
Equation Lloc u=q-Qu and (13) give u in H_W^4. For q with further
derivatives repeat: Q preserves the Sobolev order already obtained,
and (13) gains two. This proves (3)'s regularity by induction, since
D0 is a fourth-order local operator with bounded coefficient derivatives.

It remains to prove the ORIGINAL core-limit condition. Let chi_R be
a smooth radial cutoff equal to one through physical distance R and
zero beyond R+1, placed only in the asymptotic end. Its physical
derivatives through every fixed order are uniformly bounded. For
u in H_W^infinity, u_R=chi_R u is smooth, even at the bolt and compact
at the outer end, so u_R is in C_reg. Weighted H^4 integrability gives

    u_R -> u in H^1,
    W D0(u_R-u) -> 0 in L2.

The local part of K converges using weighted H^2. For the nonlocal
part use the ALREADY PROVED L2_W boundedness of Q:

    ||W Q(u_R-u)|| <= ||WQW^-1|| ||W(u_R-u)|| -> 0.

Do not try to localize the returned TT tail or drop Q after weighting.
Thus WK(u_R-u) also tends to zero. These are precisely the three
terms ||H^(1/2)(u_R-u)||, ||WK(u_R-u)|| and ||WD0(u_R-u)|| defining
D_cert. This proves actual membership, not only distributional existence
of D0u. No physical scalar, tensor or matter boundary condition changes.

## 6. The actual three-transfer identity and what remains numerical

Apply the [unchanged moment certificate](coercive-signed-moment-certificate.md)
to (3). The exact solves have r_j=q_j-Ku_j=0, so

    beta_ell=M(p)=E^-3 M(q3),
    |beta_ell| < N(q3)/(25 E^3).                   (14)

This identity now applies to the actual original residual, including
its full TT response. It is not a numerical evaluation of either side.
For the [sufficient cold target](cold-projection-value-upper-bound.md),
three exact solves still need N(q3)<=1e57 epsilon_c. Domain membership
and a finite norm do not imply that inequality.

An explicit finite-order estimate to quantify next is obtained by
letting R_s be an operator bound for

    D0 K^-1:H_W^(s+2) -> H_W^s, s=0,2,4.

The proof above gives finite such constants. Then, with
N6(p)=||p||_(H_W^6)/sqrt(gamma0),

    N(q3)<=R_0 R_2 R_4 N6(p).                     (15)

None of R_0,R_2,R_4 or N6(p) is numerically bounded here. These are
the quantities needed for a sufficient three-step estimate, rather
than an unspecified infinite smoothness assumption. Even when K=H,
D0K^-1=T, and three transfers cost six derivatives. High-frequency
data show why finite input L2 size and positivity cannot replace (15).
Resonant data can retain their signed moment exactly under (14).

The checker tests the commutator/weight and density identities, exact
tail integrals and a genuinely coupled finite-Fourier continuum model
with positive noncommuting tensor return. Its adjoint moment is unchanged
through three exact transfers, providing an explicit no-suppression
control. It does not evaluate the actual quantum kernel or high norms,
or formally verify the continuum regularity argument.

Remaining to-do list: numerical R_0,R_2,R_4 and actual weighted H6
source bounds, or a sharper actual transfer estimate; cancellation and
signed compatibility; cold/higher response, nonlinear preservation and
existence; physical boundary/state and validity, stability, matching,
charge/return, opening, transport, preparation, separate creation and audit.
