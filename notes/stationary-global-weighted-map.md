# The actual weighted stationary map and its remaining compatibility

2026-10-17 HST, E003/S139 supporting S074/S075. Ordinary argument for
independent review before focused verification. **PROGRESS.**

The complete stationary Hessian, including its actual reflecting Dirac
response, is Fredholm on the mixed weighted remainder space below.
Its index is -2. Adding the one physical mass datum raises the index to
-1, so at least one actual compatibility condition remains. This proves
a global linear mapping result; it does not prove that the preparation
source meets that condition or that a nonlinear exterior exists.

## 1. Fixed comparison, density and domains

Use the retained nonextremal comparison, positive period, fixed clock,
magnetic flux, paired kernel-free Dirac realization and subtraction
scheme of the [fast phase](stationary-fast-phase-variation.md). The
comparison is smooth on the open exterior and at its regular bolt, with
the five conormal outer jets stated there. In particular its compactified
metric has the fixed constant and quadratic jets, a mass cubic and an
O(z^4 log^p z) remainder. This includes the accepted smooth seed.
Nothing here upgrades a general C5 comparison to an all-orders smooth
compactification, or proves nonlinear C1 dependence in the spaces below.

Use the [global disk](stationary-disk-and-normal-inverse.md), with
tangents u=(S,T)=(delta Phi,delta V). Near infinity take z proportional
to -log xi. The fixed comparison's physical radial volume is

```
dmu=gamma(z) z^-4 dz=gamma(z) z^-3 dz/z,
gamma(z)=gamma0[1+O(z^2)], gamma0>0.                 (1)
```

At the bolt use smooth rotationally invariant two-dimensional disk
norms, not an interval norm with a freely chosen boundary condition.
Put theta=S+T/3 and b=2T. In the outer collar define

```
X_delta: z^-delta theta in H_b^4,
         z^-delta b in H_b^(4,log),
Y_delta: z^-delta (f1,f2) in L_b^2 x L_b^2,
dvol_b=dz/z, H_b^(4,log): <D_logz>^4 log(2+|D_logz|).
                                                               (2)
```

Patch these to the corresponding disk/interior norms by fixed smooth
cutoffs. Finite changes of those cutoffs give equivalent norms. Use
the dual transpose change for output variables. Compact smooth radial
directions, smooth at the bolt and vanishing near infinity, are dense.
At delta=7/2 the physical L2 weight is z^-2, since the physical measure
already contains z^-3 relative to dz/z. It is not z^-7/2 in physical L2.

Define A by the actual second action variation at this fixed comparison:

```
I''(u,v)=integral v dot A u dmu.                     (3)
```

This definition keeps every contact, state/domain derivative and
off-shell defect. Equivalently differentiate the action gradient
expressed per FIXED comparison volume. If F denotes its gradient per
moving physical volume, then in the disk coordinates
A u=DF[u]+(4S+2T)F at the base. The second term cannot be silently
dropped when invoking symmetry. Both systems have the same nonlinear
zeros. At the centered fast comparison F=O(z^2 log^N z) at worst;
this zero-order density correction is compact X_delta to Y_delta.

The mass-extended domain is

```
X_delta^m=X_delta plus span{e_m},
e_m=chi(z) z^3 (1,-2), chi=1 near zero,             (4)
```

in the original (S,T) variables. A nonzero fixed normalization converts
its coefficient to delta m. This is one correlated cubic, not two
independent parameters. Its bounded image is proved in section 4.

## 2. The actual boundary kernel: two powers of improvement

Let K denote the coordinate-density kernel of the Hessian before the
physical factor gamma^-1 z^4. Separate local contact distributions from
the kernel at distinct points. The nonlocal part is obtained by twice
differentiating the actual determinant; equivalently use the exact
two-vertex heat formula and contact vertex in the
[shape response](stationary-shape-logarithmic-response.md).
No separately subtracted angular modes are introduced.

Here is the finite parametrix estimate needed for the mapping result.
At a boundary point choose boundary normal coordinates, tangential
normal coordinates, a parallel spin frame and radial gauge for the
magnetic connection. The first metric jets and extrinsic curvature
vanish. The first-order Dirac expression has the form

```
D=D0+E, E=sum a^j(X) partial_j+c(X),
a^j=O(R^2), c=O(R), R=|X|.                         (5)
```

The paired boundary projectors are fixed in this frame. Their induced
square domain is the actual complementary mixed realization; its Robin
coefficient and variations are included by differentiating this
first-order problem. D0 has the direct-plus-reflected Green kernel S0
of degree -3. The first correction -S0 E S0 has degree -1: either the
R^2 coefficient cancels two powers or the R coefficient accompanies one
fewer derivative. Further insertions improve that degree. Finite Taylor
remainders and the elliptic remainder solve give the same bounds after
up to two derivatives of the Green kernel. The fixed positive gap
controls the large-time/global remainder. Its constants need not be
uniform as the prepared temperature tends to zero.

There are two scales, which must not be conflated. On a collar annulus
z,y comparable to h, with separation r much smaller than h, the Green
difference decomposes into terms bounded by

```
C h^2 r^-3, C h r^-2, C r^-1(1+|log r|)^N,         (6)
```

and their differentiated versions. The first is the change of frozen
principal symbol. At r comparable to h, all three gain two powers
against S0. Near the diagonal they remain distributions of their stated
orders, not integrable majorants. After rescaling the annulus by h,
all finite coefficient seminorms used in this calculation are bounded
by C h^2(1+|log h|)^N. The five conormal jets control this finite
expansion; it does not require arbitrary derivatives at z=0.

Two stress insertions use at most two additional propagator derivatives.
Thus away from the radial diagonal the four-dimensional difference has
degree at worst -6, with finite logarithms, against the leading degree
-8. Integrating the three tangential differences gives

```
|Delta K(z,y)| <= C (z+y)^-3(1+|log(z+y)|)^N
                         when z/y outside [1/2,2]. (7)
```

The same bound holds for the finite log-coordinate derivatives needed
below. For example integral_R3 (a^2+|v|^2)^-3 dv is a finite constant
times a^-3. A regular global Green remainder paired with a singular
propagator has degree at worst -5 before tangential integration and
hence gains at least three powers. Distant tangential regions and
thermal windings have bounded kernels here; they are retained.

At z comparable to y, (6), finite-part extension and the common heat
contacts give an annular operator bound

```
||Delta A_h||_(H4 x H4log -> L2 x L2)
                  <= C h^2(1+|log h|)^N.           (8)
```

The normalization includes gamma^-1 z^4 and dyadic rescaling. The
full-order logarithmic coefficient lies only in the shape block:
the complete conformal Ward identity makes the conformal and mixed
blocks local. Its coefficient difference is O(h^2); lower symbol
terms in (6) have the same rescaled improvement. To justify the
operator norm from these kernel bounds, first extract the universal
direct-diagonal degree -8 stress term. Its radial finite part is the
explicit fourth-order logarithmic operator, with coefficient change
O(h^2); the Ward cancellation puts this term only in the shape block.
The remaining radial diagonal kernel has order at most three, bounded
by a finite part of |t|^-4 times finite logarithms, plus local contacts
of order at most four with O(h^2) coefficients. Subtract the input's
Taylor polynomial through degree three in this finite part. On the
fixed rescaled interval H4 embeds in C^(3,alpha), for 0<alpha<1/2;
the remaining absolute integral is bounded by
integral_0^1 r^(alpha-1)|log r|^N dr<infinity. The finite-part Taylor
terms are bounded derivatives through order three, and the local
fourth derivative is bounded H4 to L2. The extracted logarithmic term
is bounded H4log to L2. The reflected Green factor is smooth on these
comparable annuli; its products with a singular direct factor are
included in the lower-order diagonal remainder. This proves (8) using the same h^2 coefficient
bounds without treating a finite part as an absolute kernel.

An alternative symbol norm formulation is supplied by
[Grubb, Theorem 2.3](https://web.math.ku.dk/~grubb/CPDE14.pdf). It gives
this bound from C^tau spatial symbol control, tau>0, and two frequency
derivatives in dimension one. Divide by <p>^4 or <p>^4 log(2+|p|)
first. The five metric jets supply, for example, tau=1/2 for the
normalized finite symbols, including the local fourth-jet contacts.
This source supplies a norm theorem; the finite expansion and the
diagonal extraction above supply the required symbol bounds.
At the bolt the corresponding argument is on a smooth disk.

Contact terms are essential to this estimate. The universal leading
contact is the one already present in the normal family. The logarithm
of the rescaling in the flat determinant cancels its Weyl partner
before taking the limit. All differences in the common local degree-four
contacts contain a positive geometric jet from (5), and retain the h^2
factor, possibly with log h. Globally dependent contacts of order at
most two are also retained. Their universal reflected singular pieces
are included in the leading parametrix; their remaining coefficients
obey the same rescaled improvement. Explicitly, the inverse-square
Green kernel in Q12 has a leading image coefficient O(z^-2). Its
curvature remainder is O(log^N z); if a derivatives land on the
direction and 2-a on that remainder, its coefficient is
O(z^(-(2-a)) log^N z), for a=0,1,2. Consequently multiplication by
z^4 and conversion of partial_z^a to b-derivatives leave the factor
O(z^2 log^N z). This retains those contacts without calling them
smoothing. Additional boundary-supported
finite-cut distributions have zero action on compact test directions;
the extension below is defined by completion of those directions.
It does not impose a different finite-cut Dirichlet problem.

For the analytic framework, [Gaarde--Grubb, Theorem 2.4 and
Proposition 2.9](https://arxiv.org/pdf/math/0703878) give the local
reflected symbol-kernel estimates for differential boundary problems
and explain why logarithmic boundary operators require more than the
ordinary small calculus. We use the finite expansion (5)-(8), not an
unqualified global classical-symbol theorem. [Falomir et al., section
III](https://arxiv.org/pdf/hep-th/9608101) explain the local contact
correction between Green-kernel regular parts and determinant finite
parts. Our contacts are fixed by the retained heat prescription.

## 3. Weighted bounds and a global Fredholm parametrix

Write x=log z. The normal off-diagonal radial kernel is bounded by
C(z+y)^-5 before multiplying by z^4. In the conjugated L2(dx) map its
kernel is bounded by z^(4-delta)y^(1+delta)(z+y)^-5. Hence its two
off-diagonal decays are

```
exp[-(4-delta)|x-x'|] when z<y/2,
exp[-(1+delta)|x-x'|] when y<z/2.                   (9)
```

They are integrable exactly for -1<delta<4. Equation (7) supplies an
extra (z+y)^2 log^N(z+y). On 0<z,y<epsilon its weighted Schur norm is
at most C epsilon^2(1+|log epsilon|)^N, with a changed constant on any
closed subinterval of these weights. Near the diagonal use (8) and a
bounded-overlap dyadic partition. The finite-part logarithm is handled
in H4log, not by taking its absolute value.

Separated compact-to-end kernels give physical O(z^4 log^N z), which
belongs to Y_delta for delta<4. End-to-compact bounds require delta>-1.
They are compact: truncate the end, use the smooth separated kernel
on the resulting compact sets, and let the weighted tail tend to zero.
The same argument applies to adjoints. We have therefore proved the
actual bounded map A:X_delta -> Y_delta, and a tail-small difference
from its normal family, for every -1<delta<4.

Take a weight whose coupled normal family is invertible. Use its bounded
Mellin inverse in the end and the accepted exact local parametrix in the
compact interior and regular disk. For a sufficiently short end, (8)-(9)
make the normal-inverse error smaller than one in operator norm.
Invert that tail error by a Neumann series. Cutoff commutators gain a
derivative, up to one logarithm, and are supported in a compact annulus
or have the separated-kernel tails just bounded. They are compact by
Rellich and tail truncation. This constructs left and right parametrices
modulo compact operators, so the actual map is Fredholm.

A full-order coefficient error is NOT declared compact merely because
it decays at the end. On a fixed compact annulus high-frequency packets
would disprove that statement. Tail smallness and the exact interior
parametrix are the two different steps used above.

These arguments are at output L2 and the indicated mixed domain. They
do not assert all Sobolev orders for arbitrary C5 data. In particular
they neither establish a nonlinear Banach-space action theorem nor
provide a numerical inverse constant on the preparation interval.

## 4. The correlated mass column is bounded

The [normal theorem](stationary-disk-and-normal-inverse.md) gives
N[z^3(1,-2)]=0: the pure TT mass cubic's flat and Weyl singularities
cancel. For its cut-off extension, comparable dyadic scales in (8)
give O(z^5 log^N z). At y>2z use (7):

```
z^4 integral_(2z)^epsilon y^-3 y^3 log^N y dy
                                      =O(z^4),    (10)
```

and y<z/2 gives O(z^5 log^N z). The cut-off normal commutator and the
far actual kernel give O(z^4 log^N z). Local contact differences obey
the comparable-scale bound. Therefore

```
A e_m=O(z^4 log^N z) in the end, A e_m in Y_(7/2). (11)
```

For this explicit smooth cut-off cubic the same annular Taylor
subtraction gives the displayed pointwise bound. General directions
use only the L2 bound (8); the mass mapping needs only its dyadic L2
version.

The extension is smooth on the compact region and can vanish near the
bolt. Changing its cutoff changes it by an element of X_(7/2), hence
does not change its class in the cokernel. This proves boundedness of
the mass-extended linear map, not that mass removes every obstruction.

## 5. Root count in the whole relevant strip

Use the normal notation q=s(s-3), Q=q+2=(s-1)(s-2),
F=digamma(s-1)+(pi/2)cot(pi s/2), H=F+kappa_M, and
L=A_T/d0. The retained benchmark gives
39/5<kappa_M<54/5<11 and L>30*2000^10>400.
The TT bracket is QH-L. We now exclude its zeros throughout
-1/2<=Re s<=7/2, not just on the previously proved endpoint contour.

By F(3-s)=F(s), it suffices to take s=delta+it with
3/2<=delta<=7/2. For |t|<=2, put a=delta-1 in [1/2,5/2]. The series
for digamma gives

```
|digamma(a+it)-digamma(a)|
 <= |t| sum_(n>=0) 1/(n+a)^2 <=12,
|digamma(a)|<3, |digamma(s-1)|<15, |Q|<9.          (12)
```

Here the sum is at most 1/a^2+1/a<=6. To bound the removable cot
product put v=pi(s-2)/2. Then |v|<=5pi/4,
|v-pi|>=pi/4 and |v+pi|>=3pi/4. The partial fraction identity
v cot v=1+sum_(n>=1)2v^2/(v^2-n^2 pi^2) gives

```
|v cot v|<1+50/3+200/39<23,
|Q (pi/2)cot(pi s/2)|<75,
|QH|<9(15+11)+75<400.                             (13)
```

For n>=2 use n^2-25/16 >=39n^2/64 and
sum_(n>=2)n^-2<1. This estimate includes s=2 by its removable limit.

For t>=2 the decreasing summand in
Im digamma(a+it)=sum_(n>=0)t/[(n+a)^2+t^2] has integral
pi/2-arctan(a/t). Comparing sum and integral proves

```
-a <= t[Im digamma(a+it)-pi/2] <=1.
```

The explicit cotangent formula gives
|t[Im((pi/2)cot(pi s/2))+pi/2]|<5t exp(-pi t)<5/32.
Thus |t Im F|<3. Also Re digamma(a+it)>=digamma(a)>-2,
and the real cot term has absolute value below 1/2, so Re H>53/10>5.
For these frequencies Re Q=delta^2-3delta+2-t^2<0 and
0<=Im Q/t=2delta-3<=4. It follows that

```
Re(QH)<12<L.                                      (14)
```

Complex conjugation treats negative t. Equations (13)-(14), and
reflection to the left half, exclude every zero of the TT bracket.
The only TT zeros in the strip are consequently the simple roots 0
and 3. At s=1,2 the product P4 F is -2, not zero; its pole is removable
in N_T. The genuine Mellin poles -1 and 4 lie outside the strip.

The scalar normal block has exactly the simple roots
3/2+/-i omega in the strip. Its roots -1 and 4 are outside. The
decomposition C and its reflected adjoint are holomorphic invertible
there. Hence det H_disk=4 J_c N_T has exactly four zeros, counted
with multiplicity, in -1/2<Re s<7/2. None lies on the two edges.

## 6. The global index and actual adjoint moments

For clarity, the weight-change argument also applies to this logarithmic
operator. After the mixed order reduction its normal inverse is
meromorphic on the open strip, with the four finite-rank simple poles
just counted, and has the required large-frequency symbol bounds.
One can choose that reduction holomorphically: for
p(s)=100-(s-3/2)^2, Re p>=96+(Im s)^2 on this closed strip. The
diagonal normal reductions p(s)^2 and p(s)^2 log p(s), using the
principal logarithm, have no zeros there and have the required orders.
No analytic continuation of the absolute value in (2) is assumed.

For the finite asymptotic step, retain coefficients as well as bounds
in the construction of section 2. Taylor-expand the Dirac, projector
and contact coefficients through relative degree three, and the regular
global Green remainder through its constant term at the boundary.
The latter can contribute a degree-three term: its two-vertex cross
kernel has degree -5 before tangential integration. Its coefficient
can depend on the actual state and global Green regular part. Retain
it. Define this global term using cutoff sources separated from the
boundary diagonal, keeping local singular parametrix terms separate.
It solves a homogeneous elliptic boundary problem in the observation
collar. The compactified coefficients are at least C^(2,alpha), with
the induced mixed boundary data of sufficient C^(1,alpha) regularity,
for any fixed 0<alpha<1. Local boundary regularity gives a C1 Taylor
bound for this separated term. Subtracting its boundary constant
therefore gains one power. This uses no globally smooth symbol
expansion across the z^4 log^p boundary jets; those jets stay in the
finite local remainder of degree at least four, with finite logarithms.
The globally dependent inverse-square Green contact coefficients must
also be retained: their constant and linear Taylor terms can contribute
at relative degrees two and three. The separated inverse-square Green
remainder solves the homogeneous mixed second-order problem. The
compactified metric is C^(3,alpha), and its induced coefficients and
boundary data give a C2 Taylor bound for this remainder. After the
constant and linear terms, its quadratic remainder contributes at
relative degree at least four. Thus both the nonlocal Dirac term and
every finite contact coefficient are extracted; none is dropped on
grounds of global dependence. The five conormal metric jets suffice
for this finite extraction.
On cut-off powers used below this yields

```
A=N+z^2 sum_j (log z)^j A_(2j)
   +z^3 sum_j (log z)^j A_(3j)+R.                  (17)
```

The sums are finite homogeneous Mellin operators at the indicated
degrees; both A_(2j) and A_(3j) include every contact coefficient,
including globally dependent finite Green coefficients. A_(3j) also
includes the nonlocal global Dirac coefficient just described. At the
degree-two or degree-three far-integration thresholds their meromorphic
finite parts produce O(z^4 log^N z), which is beyond the final weight.
The remainder R has this same output order or the local input exponent
plus four, whichever is smaller. In particular it sends the relevant
leading powers of real parts 0, 3/2 and 3 strictly beyond 7/2. This
coefficient extraction, not the big-O inequalities alone, supplies the
finite expansion used in the index argument.

Moving a Mellin contour across a pole adds its residue mode. The
remainder estimate gains two powers, up to finite logarithms. Iterate
the normal correction on finitely many overlapping weight intervals
of length less than two; away from roots it improves decay, and at a
root its finite residue records the added asymptotic datum. Compact
remainders in the constructed parametrices contribute the same finite
kernel/range terms on both sides of the relative-index exact sequence.
Thus increasing the decay weight from -1/2 to 7/2 decreases the index
by four.

The representatives include their forced corrections. The complete
degree-two coefficient, including contacts, sends each leading scalar mode z^(3/2+/-i omega) to a
source of degree 7/2+/-i omega, exactly at the final weight. The normal
matrix there is invertible; subtract its particular solution before
claiming a remainder in X_(7/2). The next error has strictly larger
real exponent. The TT source mode at s=0 similarly needs its forced
degree-two correction. The degree-three background can require a
degree-three logarithmic correction at the mass root; meromorphic
inversion supplies that fixed particular term, while its homogeneous
mass coefficient is the independently counted s=3 datum. These
corrections use every degree-two and degree-three coefficient retained
in (17); the O(z^4 log^p z) remainder is
already beyond the required strip. This
finite construction prevents leading powers alone from being mistaken
for representatives in the range-jump sequence.

This is the contour and boundary-pairing proof of the relative index
formula, as in [Melrose, sections 5.17 and 6.2](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/melrose.pdf).
The finite-strip kernel estimates and logarithmic order reduction above
are supplied explicitly rather than assuming membership in a classical
calculus. No continuation through the poles at -1 or 4 is used.

By (1) and the exact symmetry (3), the complementary weight is 3-delta.
The weight -1/2 space is an auxiliary linear completion of compact
tests, not a new physical nonlinear phase with growing boundary data.
The cokernel at delta identifies with actual distributional solutions
of A chi=0 in the complementary L2 weight; the two-sided parametrix
upgrades them to the complementary mixed domain, with regular disk
behavior. This uses a finite bootstrap through four derivatives in the
graded symbol expansion: the leading coefficients have five bounded
conormal derivatives, while lower-order coefficients need fewer. It
does not infer all-Sobolev regularity from the C^(1/2) order-zero norm
bound or require all-orders boundary smoothness. Consequently the
indices at 7/2 and -1/2 are opposites.
Writing i_delta for the index gives

```
i_(7/2)-i_(-1/2)=-4, i_(-1/2)=-i_(7/2),
index(A:X_(7/2)->Y_(7/2))=-2,
index(A_m:X_(7/2)^m->Y_(7/2))=-1.                 (15)
```

The last step is the finite-dimensional domain extension (4), whose
bounded image was proved in (11). It does not require transversality.
If K_- is the actual regular-bolt null space of A at weight -1/2, then
the cokernel of the mass-extended map is precisely

```
C_m={chi in K_- : integral chi dot A e_m dmu=0}.
```

The exact range condition for A_m u=f is

```
integral chi dot f dmu=0 for every chi in C_m.      (16)
```

All integrals are convergent by the complementary weights. If the
mass-extended kernel has dimension k, then dim C_m=k+1. In particular
there is at least one compatibility condition even with mass free.
These are the actual coupled adjoint fields, not the earlier scalar
conformal moments or a normal-power substitution. The cutoff choice
in (4) does not affect (16). A bounded inverse exists on chosen closed
kernel and range complements, but its norm is not quantified here.

## 7. Consequence for the solver and the physical goal

The full-action route now has an actual global Fredholm map and a
bounded physical mass column. It cannot be used as an unrestricted
inverse onto arbitrary fast sources. The next decision requires the
actual coupled kernel, adjoint fields and projection of the preparation
defect, followed by a quantitative inverse on the compatible range.
An extra freely adjustable datum cannot be invented to remove (16).

The weight 7/2 is an auxiliary remainder domain. It does not itself
establish the physical O(z^4 log^p z) phase, choose its logarithmic
boundary data, or prove nonlinear differentiability and closure there.
A negative index is not nonexistence: the actual source may satisfy
the remaining condition, and a nonlinear solution still needs proof.
No candidate equation defect is converted into a distance to a solution.
The physical action/state, Hamiltonian, S138 transfer and S074-S080
requirements remain open.

Focused verification may check the finite kernel scaling, weighted
integrals, strip-bound rational arithmetic, root multiplicities and
index/domain-extension bookkeeping. The Green-parametrix estimates,
finite-part contact analysis, Fredholm construction and actual range
theorem remain ordinary arguments. No root numerics or inverse-norm
fit can replace them. Review the complete argument before execution.

Remaining to-do list: determine the actual coupled kernel and adjoint
compatibility projections; quantify the inverse on the compatible range
and close the nonlinear physical logarithmic phase, or justify controlled
spatial reduction; match exterior/state and charge, transfer S138, and
complete driven stress, radius, controller, capture, tail, spatial area
and S074-S080.
