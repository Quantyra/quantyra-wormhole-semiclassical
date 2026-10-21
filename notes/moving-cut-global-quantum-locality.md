# Moving-cut global quantum locality

2026-10-21, E003/S168. Ordinary argument before focused verification.
For the fixed charged seed, paired domain and prescribed Euclidean state,
the massive two-vertex response coming from separated regions or from
the difference between two identical local geometries has an ordinary
finite moving-cut limit. Its global diagonal terms have the same limit.
This isolates the unresolved comparison to the local curved boundary
kernel and its ultraviolet contacts. It does not evaluate that comparison.

Keep the original angular tower, common subtraction, finite coefficients
and exponential mass chart. Use [the chart-contact result](energy-domain-and-chart-contact.md)
and [the retained actual comparison](stationary-curved-quantum-tt.md).
The statement concerns the fixed Euclidean seed. It does not choose the
physical gravitational boundary phase or establish a causal state.

## 1. A fixed smooth half-domain, rather than a smooth double

Let Y be inward canonical proper distance from the limiting outer end.
The actual half-metric dY^2+a(Y)^2 dtau^2 is smooth through Y=0 from
the right, a(0)>0 and a'(0)=0. A cut at Y=d has extrinsic curvature
of order d. Its even double contains a'(d)|x|, which is only Lipschitz.
Thus the C2,1 seam comparison at d=0 is not uniform evidence at d>0.

Fix a smooth rho which equals one near Y=0 and zero well before the
bolt. For 0<=d<=d0 small, use

    F_d(x)=x+d rho(x),  1/2<F_d'(x)<3/2.                 (1)

This maps the limiting half-domain onto the cut domain and is identity
near the bolt. Pull back spinors with the volume half-density and the
normal spin-frame identification. The time period, spin structure and
magnetic bundle remain fixed. The pulled-back operators A_d=Dsharp_d
then have a common domain

    D={u in H1: Pi_- u|boundary=0},
    A_d^*=-A_d,  ||(A_d+z)^-1||<=1/Re z, Re z>0.        (2)

Here Dsharp=-J D2 is the actual paired realization of
[the cold-bolt proof](stationary-reference-cold-bolt.md), not the
unpaired D2. Its normal symbol anticommutes with chi and its tangential
symbol commutes with chi. The normal-times-tangential symbol therefore
exchanges the chi eigenspaces, giving the elliptic boundary condition;
the normal Green form pairs those eigenspaces, giving self-adjointness
after multiplication by i. These statements are unaffected by d.

The general boundary regularity and self-adjointness results used here
are [Baer--Ballmann, Theorems 3.9 and 3.11, Corollary 3.17](https://arxiv.org/pdf/1307.3021v1).
They are applied to smooth half-manifolds. Uniformity below follows
from the explicit small smooth perturbation in (1), not from a theorem
about a smooth doubled metric.

On this common domain, for each fixed finite k,

    ||A_d-A_0||_(H^(k+1)->H^k) <= C_k d,
    ||u||_H1 <= C(||A_d u||_2+||u||_2).                (3)

The seed is fixed and smooth on its compact canonical disk, with a
regular bolt. A finite smooth cover therefore supplies finite constants.
The H1 estimate at d=0 persists by absorbing C d ||u||_H1. Local
higher estimates persist the same way with nested cutoffs and the same
boundary projector. This does not assert bounds uniform as temperature
tends to zero, or use more quantitative seed-jet bounds than were proved.

With R_d(z)=(A_d+z)^-1, the resolvent identity and (3) give, for real M>0,

    ||R_d(M)-R_0(M)|| <= C d (1+1/M)/M.                (4)

Indeed A_0 R_0=1-M R_0 bounds the H1 norm of R_0 by C(2+1/M).
The same argument works on |z-M|<=M/4 with a fixed larger constant.
The pulled-back coefficients depend smoothly on d; away from a source,
local regularity also gives convergence of resolvent kernels and their
finite jets up to the smooth half-boundary.

## 2. Separation gives a uniform summable estimate

Multiplication by a real smooth weight preserves D. If |grad f|<=1,

    exp(alpha f)(A_d+z)exp(-alpha f)
                         = A_d+z-alpha sigma_d(df).

For alpha<=Re z/(2 c_0), where ||sigma_d(df)||<=c_0, a Neumann series
and (2) bound its inverse by 2/Re z. Weights with fixed separated
plateaux therefore give an exponentially small L2 resolvent between
sets separated by a fixed positive distance. Constants are uniform
for small d; the metric equivalence in (1) preserves the separation.

If xi and eta have separated supports, apply local boundary or interior
regularity to (A_d+z)R_d(z)eta=0 near supp xi. Nested cutoffs, all with
positive source separation, give

    ||xi R_d(z) eta||_(L2->H^s)
                    <= C_s (1+M)^s M^-1 exp(-c M),     (5)

for any fixed integer s and |z-M|<=M/4. Shrink c when choosing the
nested buffers. The positive lower mass M_1 absorbs negative powers.
The adjoint has the same estimates. The constants and c are finite
for this fixed seed, cover and period; no numerical calibration is made.

On a compact two-dimensional patch the inclusion H^s->L2 is trace
class for s>2, and Hilbert--Schmidt for s>1. One can check this on
rectangles by summing (1+m^2+n^2)^(-s/2), then use bounded extension
from half-rectangles and a finite partition. Consequently (5) gives
trace-norm estimates using s=4 and Hilbert--Schmidt estimates using s=2.
No pointwise square of the singular Dirac kernel is integrated here.

## 3. The actual-to-cap difference, including both vertices

For each patch choose an auxiliary complete smooth half-geometry or
interior geometry identical to the actual pulled-back metric and domain
on a larger patch. Its coefficients form a smooth family in d. Complete
it outside that patch with fixed smooth collars; it may be a cylinder
with compact boundary and a fixed end. This preserves the paired
elliptic realization and (2). Cutoffs and completions are independent
of M. Let B_d(z) be its resolvent, chi a cutoff supported in the common
patch, and xi a smaller cutoff with chi=1 near supp xi. Set Q_d=[A_d,chi].
The two geometric resolvent identities give

    E_d(z)=xi(R_d(z)-B_d(z))xi
                    =-xi B_d(z) Q_d R_d(z) Q_d B_d(z)xi. (6)

The intervening single-commutator term is
(chi^2 B_d-B_d chi^2)/2 and vanishes between the inner cutoffs.
Every outer leg in (6) has a positive separation from its commutator.
Equations (2), (5) and the compact Sobolev inclusions imply

    ||E_d(z)||_1 <= K(1+M)^12 exp(-c M),               (7)

uniformly in d and the mass circle. The degree 12 is a conservative
reserve: four regularity steps on either end and at most four additional
mass powers suffice for all the expressions below. The same envelope,
with a larger K and smaller c, bounds the finite diagonal jets needed
below. For fixed M these expressions converge in trace norm and in
those diagonal jets as d->0. To see trace-norm convergence, use (4)
and local regularity on the separated factors in (6), with an extra
regularity step if necessary, and factor through H4->L2.

Let T be supported where xi=1. Multiplication by T is bounded, and
write Rloc=xi R xi, Bloc=xi B xi. The local massive bubble difference is
defined without taking a divergent trace of either bubble separately:

    M^2 Tr(T Rloc T Rloc-T Bloc T Bloc)
     =M^2 Tr(T Bloc T E+T E T Bloc+T E T E).             (8)

Each right-hand term is trace class, bounded by an envelope of the
form K(1+M)^12 exp(-cM)||T||_infinity^2. Exponentials absorb any
extra polynomial reserve after shrinking c. The mixed version uses
||T1||_infinity ||T2||_infinity. The global part of the exponential
second mass vertex is also retained:

    -M Tr(T^2 E).                                      (9)

All local measure and common-subtraction terms cancel in this
actual-to-cap difference because the geometries, domains and variation
jets agree on the support. Equations (8)-(9) therefore identify a finite
part of the ORIGINAL Hessian. They do not renormalize a bare tadpole
or replace the original exponential chart by an affine chart.

For a pair of separated supports, the two resolvents in the actual
bubble are Hilbert--Schmidt by (5); their product is trace class with
the same exponential envelope. A finite cover separates the complete
quadratic response into these pieces and local cap responses: use a
partition psi_j(x), and for each j a cutoff chi_j(y) equal to one near
supp psi_j. The remaining psi_j(x)(1-chi_j(y)) pairs are separated.
The local pieces use (8) with polarized insertions. The cover can include
the outer boundary, the smooth interior and the bolt.

## 4. Global diagonal terms are not omitted

Put F_d(z;x)=tr(R_d(z)-B_d(z))(x,x) on the inner patch. Equation (6)
is holomorphic on the mass circle and its separated factors control
the diagonal. Cauchy's formula gives

    |M^2 partial_M F_d(M;x)+M F_d(M;x)|
                         <=5M sup_|z-M|=M/4 |F_d(z;x)|. (10)

The left expression, with a minus sign, is precisely the difference
of the complete constant conformal probes:

    M^2 tr(R_M^2)-M tr R_M
                           =-M^2 partial_M tr R_M-M tr R_M.

Keeping only the chart contact would lose the mass derivative, which
contains the two-vertex row sum. Keeping only a row-subtracted bubble
would likewise lose its compensating diagonal. Both contributions
converge here. For a localized row subtraction use the row over the
same localization region, not the full constant-probe row. The raw
trace formula (8) already avoids that bookkeeping ambiguity.

This argument differentiates a holomorphic identity with uniform
complex-mass bounds. It does not differentiate an unproved remainder
estimate. Local ultraviolet diagonal contacts, which are shared by
actual and cap, have not been evaluated by (10).

## 5. Angular summation and the moving-cut limit

Retain M_n^2=n(n+q)/C and d_n=2(q+2n). The exact counting bound is
N(E)=sum_(M_n<=E) d_n<=4 C E^2. Define

    A(p,a)=sup_(x>=0) (1+x)^p exp(-a x)
          =1 if p<=a, otherwise (p/a)^p exp(-p+a).

Stieltjes integration of the counting function gives

    sum d_n exp(-b M_n) <=8 C/b^2,
    sum d_n (1+M_n)^p exp(-c M_n)
                                 <=32 C A(p,c/2)/c^2. (11)

A useful uniform tail bound above mass L is

    sum_(M_n>L) d_n (1+M_n)^p exp(-c M_n)
                   <=128 C A(p,c/2)c^-2 exp(-c L/4).   (12)

The fixed mass gap is positive. Thus (7)-(10) and the separated
Hilbert--Schmidt products can be summed by dominated convergence,
including the diagonal mass derivative. No angular truncation defines
the theory, and no common ultraviolet counterterm is altered.

For an energy-class outer shape, T_d(x,tau)=v(F_d(x))H(tau), with a
fixed smooth finite-frequency H, converges in H2 on a collar and hence
uniformly. This follows first for smooth v, then by density and the
uniform H2 bounds for composition by (1). Use the corresponding
partition for the rest of a regular tangent. Consequently the complete
global massive response described above has an ORDINARY finite limit.
In particular it has the translated logarithmic weak limit required
by the preceding obstruction test, and supplies no persistent oscillation.

The dependence on the prescribed thermal identification, bolt and
remote attachment remains in the actual resolvents. It is not replaced
by a vacuum kernel or discarded. The theorem controls the global
completion contribution for this fixed Euclidean state; it does not
construct a physical exterior/state or prove real-time stability.

## 6. The remaining decisive comparison

The earlier normal-model and chart-contact limits, together with this
global locality result, leave the local curved finite-cut cap Hessian
minus the normal representative. Its near-diagonal and reflected kernels,
fixed local contacts, curvature and all Weyl/LLL/classical shape partners
must still be combined. The positive-cut first-derivative jump is exactly
why the limiting C2,1 geodesic proof alone does not settle this step.
No full pure-shape cutoff limit, branch exclusion or nonlinear solution
is claimed. S168.3-S168.4 and all physical gates remain unchecked.

Verification tests the exact two-sided gluing and mass-derivative
identities, moving-operator continuity, the separated trace decomposition
and angular counting/tail bounds. Finite matrices are diagnostics, not
evaluations of the actual charged kernel or formal continuum proofs.
In a finite matrix model Q=[A,chi] need not commute with chi. Its exact
version of (6) includes +xi B[Q,chi]B xi/2. For a first-order differential
operator Q is multiplication by sigma(dchi), so this defect is zero.
The initial diagnostic omitted the finite-matrix defect and failed;
its frozen source and run are retained. The corrected diagnostic keeps
the defect and separately checks the differential multiplication identity.

Remaining to-do list: local curved finite-cut quantum comparison with
all local shape partners; physical boundary phase and source compatibility;
continuous causal constraints, exact nonlinear existence and all-sector
stability; exterior/state, charge/return, opening, transport, preparation,
separate creation and full audit.
