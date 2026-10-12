# Electric Gauss law selects a singular mouth boundary in the current model

2026-10-11, E003/S053. The regulated annulus determinant result does not
determine the global gauge modes. In the prescribed spherical Maxwell/LLL
current reduction, a coupling that is small at the throat center nevertheless
produces a strong charged-current barrier near each mouth. This note derives
a uniform resolvent estimate for replacing that barrier by two collective
boundary constraints, retaining the unequal exterior returns. It also checks
the bare electric inertia of the matched classical regions.

This is a result for the stated quadratic current model. Compact flux and
spin sectors, the actual exterior photon, angular interactions, the quantum
state and renormalized stress remain open. In particular, neither an isolated
rotor nor a fixed annulus covariance determines the physical vacuum. The
Schwinger mechanism and large-coupling boundary limits are established methods;
no physical damping, stability or consequential novelty is accepted.

## Physical normalization and a nonuniform weak-coupling limit

Use inner dimensionless time t=t_phys/ell, r=tan(s), and the leading product
metric P^2[-(1+r^2)dt^2+dr^2/(1+r^2)+dOmega^2]. Write lambda=P/ell.
The unit-charge Maxwell action is -F^2/(4g_c^2). Its spherical electric
coefficient in (t,r) is 4pi/g_c^2. In optical coordinate s it is

```
A_e(s)=4pi cos^2(s)/g_c^2,
e_2(s)^2=1/A_e(s),
V(s)=a sec^2(s),                 a=q g_c^2/(4pi^2).          (1)
```

Indeed the canonically normalized current chi=q^(-1/2) sum_j phi_j
couples as sqrt(q/pi) chi F_ts. Gauss law fixes the displacement
A_e F_ts+sqrt(q/pi) chi in a sector. Completing its electric square gives
(e_2^2/2)(n-sqrt(q/pi) chi)^2. The Hessian of fluctuations about a fixed
background displacement is V chi^2/2. This calculation fixes neither the
allowed n nor its wave function, and supplies no missing exterior action.
It agrees with [MMP, exact v3, Appendix A](https://arxiv.org/html/1807.04726v3#A1)
and the existing [singlet reduction](schwinger-return-feedback.md).

Truncate at r=+/-R, put epsilon=arctan(1/R), and T=pi-2epsilon. The
common optical interval is [-T/2,T/2]. For each finite R the potential
is smooth and bounded; its integrated strength is exactly

```
I=int V(s) ds=2a R.                                         (2)
```

In the retained family g_c^2=q^(-2), lambda=(16pi q^3)^(-1), a matching
cut R=q^2 is within the classical overlap: R tends to infinity and
lambda R=1/(16pi q) tends to zero. Yet I=q/(2pi^2), which grows. Thus
the small-integrated-potential hypothesis of
[the earlier static-curvature bound](singlet-static-curvature.md), equation
(3), fails at this deeper cut. That earlier theorem used epsilon=q^(-1/4)
and remains valid there. Extending it to the physical mouth regions requires
an additional calculation. A resolvent estimate below is not a bound on
its ultraviolet-subtracted vacuum energy or clock derivatives.

## Bare global electric inertia

For a prescribed source-free spherical electric displacement n, the inner
energy in throat units is n^2 g_c^2 R/(4pi). Add two RN radial regions
from dimensionless areal radius r_s=1+lambda R to b>r_s. RN time is
t_phys/P=t/lambda. Their combined energy in throat units is
n^2 g_c^2(1/r_s-1/b)/(4pi lambda). Consequently this explicitly matched
classical comparator gives

```
E_e=n^2/(2 I_e),
I_e=2pi lambda/[g_c^2 F],
F=lambda R+1/(1+lambda R)-1/b.                              (3)
```

This follows either by integrating the Maxwell Hamiltonian or by integrating
the electrostatic potential difference and using E=n Delta A_t/2. It is
the bare capacitance of these prescribed regions, not the fermion-screened
capacitance or a completed one-exterior binary. The remaining far region is
not included. For b tending to infinity, F=1+O((lambda R)^2), and (3)
gives I_e~1/(8q), rather than an inertia growing with inverse microscopic
coupling. The redshift and the two mouths are essential to this conclusion.

The analogous *bare* angular comparison is already encoded by the paired
slow Kerr-Newman energy: ell E_spin=16 S^2/q, giving I_rot=q/32 for a
relative angular coordinate conjugate to S. This is consistent with the
existing relative-frame response Delta w_t=2J/lambda and S=KJ. It does
not identify that electrovacuum charge with the matter-coupled conserved
Q=K J_geom+sqrt(T_q/pi) chi, or quantize its three noncommuting axes.

Static free-fermion holonomy curvatures on equal circles of length L are
q/(pi L) for the electric angle and T_q/(pi L) for a rotation about one
axis. Combining them with these bare inertias gives nominal oscillator
frequencies of order q, above the order-one throat frequencies. Therefore
no slow-variable justification for an instantaneous-vacuum rotor follows.
Furthermore, large gauge transformations relate fermion sea sectors; one
cannot replace their spectral-flow matching by periodically minimizing a
single sea's parabola. [Azakov, exact v1, section 2.3](https://arxiv.org/html/hep-th/0511116v1#S2.SS3)
explicitly retains that matching in the circle Schwinger model. These tests
reject this shortcut to the state; they do not predict a physical instability.

## The boundary replacement theorem with unequal returns

Let q>=1 circles have a common coordinate interval of length T and exterior
lengths D_j>0. Their fields remain distinct. At each mouth impose the
original channelwise continuity. On the direct-sum L2 space define

```
v=(1,...,1)/sqrt(q),
h_0[f]=sum_j int_circle_j |f_j'|^2,
h_a[f]=h_0[f]+int_0^T a sec^2(x-T/2) |v^* f(x)|^2 dx,
B f=(v^* f(0), v^* f(T)).                                   (4)
```

H_a is the operator of h_a on the common H1 circle domain. H_D is the
operator of h_0 restricted to Bf=0. This restriction is dense in L2 and
closed in the form norm. It constrains two collective endpoint values,
not every channel separately. The neutral endpoint values still match.
Unequal exterior lengths need not preserve the neutral subspace in transit.

Choose 0<h<=T/2 and define

```
beta=(a/2)(1/epsilon-1/(epsilon+h)),
gamma=a log(1+h/epsilon),
C=2[1/T+avg_j(1/D_j)]+[T+avg_j D_j]/2,
c=pi^2 a.
```

If gamma<1, the following explicit estimate holds in the same L2 space:

```
||(H_a+1)^(-1)-(H_D+1)^(-1)|| <= E,
E=min(1,max(c/(1+c), gamma/(1-gamma)+C/(beta+C))).             (5)
```

It is independent of the largest individual return length except through
the displayed mean lengths. No uniform free gap, small integrated mass,
fixed photon annulus, or deletion of a low-frequency interval is used.

### Boundary-layer lower bound

Put z=v^*f on the throat. Near its left end cos(x-T/2)=sin(epsilon+x)
and hence V>=a/(epsilon+x)^2. The fundamental theorem of calculus gives

```
|z(x)|^2 >= |z(0)|^2/2-x int_0^x |z'(y)|^2 dy.
```

Integrate on [0,h]. The positive term is beta |z(0)|^2; the negative
term has magnitude at most gamma int_0^h |z'|^2. Repeat at the right
end. The layers are disjoint and |v|=1, so

```
h_a >= (1-gamma) h_0 + beta |Bf|^2.                         (6)
```

### Upper bound on the constrained domain

For d(x)=min(x,T-x), sin(epsilon+d)>=2d/pi. When Bf=0, the interval
Hardy inequality gives int |z|^2/d^2<=4 int |z'|^2. It follows that

```
h_a[f] <= (1+c) h_0[f],                 f in ker B.           (7)
```

For completeness the Hardy inequality on each half interval follows by
integrating (|z|^2/x)' and applying Cauchy-Schwarz; the outer endpoint
term has the favorable sign. Density extends it from smooth functions
with vanishing endpoint trace to H1 functions. No singular-endpoint
choice has been inserted into H_a.

### Penalty resolvent and the uniform constant

Let R_0=(H_0+1)^(-1), G=B R_0 B^*, and let H_beta have form
h_0+beta|Bf|^2. Lagrange multipliers, or direct completion of the
positive form square, give

```
R_beta=R_0-R_0 B^*(beta^(-1)+G)^(-1) B R_0,
R_D   =R_0-R_0 B^*G^(-1) B R_0.                             (8)
```

The adjoints in these formulas are interpreted between the form domain
and its dual. In particular G is a finite positive 2 by 2 matrix.
To bound its inverse, interpolate f=v b linearly between b_0,b_1 on the
throat and between b_1,b_0 on each exterior interval. This is a right
inverse of B, with

```
h_0[f]+||f||^2 <= C (|b_0|^2+|b_1|^2).
```

The minimum of the left side at Bf=b is b^*G^(-1)b, so G>=C^(-1)I.
Writing W=R_0^(1/2) B^*G^(-1/2) gives W^*W=I and

```
R_beta-R_D=R_0^(1/2) W (I+beta G)^(-1) W^* R_0^(1/2),
0<=R_beta-R_D<= [C/(beta+C)] I.                              (9)
```

Equations (6) and (7), and the variational characterization of inverses,
sandwich R_a between the inverse of (1+c)H_D+1 and that of
(1-gamma)H_0+beta B^*B+1. Their distances below and above R_D are at
most c/(1+c) and gamma/(1-gamma)+C/(beta+C), respectively. This proves
(5). Positivity and 0<=R_a,R_D<=I supply the final cap by one. This is
a form/penalty argument of an established type, not a new general
large-coupling theorem.

## Consequence at the proposed damping resolution

Take R=q^2, a=1/(4pi^2 q), h=T/4, and the midpoint point-dipole returns
D_j=delta f((j+1/2)/q), delta=log^(-2)q. The existing geometric bounds
give D_min>=c_0 delta and avg D_j=O(delta log q). Then

```
beta~q/(8pi^2), gamma=O(log q/q), C=O(log^2 q),
E=O(log^2 q/q).                                             (10)
```

For zeta=(omega+i eta)^2 with omega in [1/2,3/2], eta>0, the resolvent
identity and nonnegativity give

```
||R_a(zeta)-R_D(zeta)||
 <= [1+|1+zeta|/(2 omega eta)]^2 E.                          (11)
```

Here R_H(zeta)=(H-zeta)^(-1). In particular eta=delta yields
O(log^6 q/q), tending to zero. The two operators in each comparison
use the same T and D_j; no identification of different-q Hilbert spaces
is assumed. This transfers bounded resolvent readouts at that resolution,
not stress tensors, coincident derivatives, or a real-axis uniform limit.

For equal exterior lengths the comparison is especially transparent:
H_D consists of q-1 free circles, one charged throat Dirichlet interval,
and one charged exterior Dirichlet interval. This differs from a large
*bulk* mass that removes the charged throat field altogether. For unequal
returns retain the constraints (4), rather than substituting q-1 loops.

The uniform charged test field has Rayleigh quotient 2aR/(T+D), which
diverges at R=q^2. Nevertheless a singlet test function cos(pi s/T)
on the throat, zero outside, has Rayleigh quotient at most
(pi/T)^2+2a. For D<T, (5) makes the lowest charged frequency tend to one.
The diverging constant-field estimate therefore is not the physical gap
even in this quadratic model. The spatial current response is necessary.

## Verification, prior work and decision

[The checker](../checks/gauss_mouth_boundary.py) verifies the Maxwell
coefficients, independent energy/potential integrals, the matrix resolvent
identity, the continuum Green-function lift bound, graded finite-element
boundary convergence, an independently integrated charged-mode ODE, and
retained-family constants. The ODE comparison initially exceeded the
declared absolute eigenvalue tolerance at 64 mouth cells; refinement to
128 cells addresses that error without loosening the tolerance. Its
[receipt](../evidence/gauss-mouth-boundary.json) pins exact inputs and
refinements. Numerics supplement the ordinary proof; no Lean was run.

[BelHadjAli, Ben Amor and Brasche, December 2009 preprint](https://www.mathematik.tu-clausthal.de/fileadmin/Sites/Mathematik/Documents/Mathematik_Berichte/2009_10.pdf)
reviews large-coupling convergence and operator formulas; those methods
are prior work. The [source record](../evidence/gauss-mouth-source-review.json)
discloses the bounded proof reading and searches. No comprehensive novelty
finding follows from this continuation.

Continue physical long-return damping as primary and matched finite-time
transfer as fallback. Use (5) for the quadratic electric mouth limit.
Do not extrapolate the small-I static estimate through the rest of the
mouth or replace the full state by an instantaneous-vacuum rotor. The
angular/global state and the renormalized mouth stress now determine the
next physical test. No separate graph-formalization campaign is selected.

Remaining to-do list: derive the angular and compact gauge sectors with
the physical exterior, control the coupled state and renormalized mouth
stress, complete smooth matching and the physical damping/transfer result,
refresh novelty, independently verify, GCP-certify suitable mathematics
and prepare the reproducible manuscript.
