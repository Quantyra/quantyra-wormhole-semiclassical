# The moving mouth adds pressure to the clock vertex

2026-10-11, E003/S051-S053. The finite product-throat Gaussian model has
an explicit linear moving-mouth vertex. Including the momentum coupling
from the moving coordinate system cancels the apparent bulk mass-profile
correction: on shell the vertex is the throat Hamiltonian plus the sum of
the two mouth pressures, multiplied by epsilon/sqrt(1-epsilon^2).

This identifies a real limitation of the previous canonical-clock bound.
Already in the free CFT, the mouth-pressure spectral response requires
more than a static subtraction in the sharp-mouth limit. We compute its
divergent local inertia and prove a uniform sufficient smoothing bound.
For optical smoothing width h comparable to epsilon, the free spectral
correction is small at the retained damping resolution. The physical
join has not been shown to supply that smoothing or the matching local
terms. This is progress on the prescribed product model, not acceptance
of physical wormhole damping, stability or novelty.

## 1. Pull back the moving interval before varying the action

Use throat-radius units, the same Gaussian oscillator-sector convention
as [the finite determinant](finite-mouth-determinant.md), q real channel
fields f, and the fixed normalized charged vector v. The rank-one
potential is a sec^2(sigma) vv*, with a>=0. Exterior optical return
lengths D_j>0 stay fixed. Fields are continuous at the identified mouth
traces; no reflecting boundary condition is imposed on the throat.
The background clock is x=exp(phi)=1. Write

```text
A=acos(epsilon), b=sqrt(1-epsilon^2), T=2A,
r=epsilon/b, kappa=-r/A,       0<epsilon<=1/4.
```

The [proper-cut condition](proper-clock-anomaly.md) is
cos^2(s)=epsilon^2(x^2-sdot^2). To linear order s=A+kappa A phi.
The quadratic velocity correction to s is irrelevant to the linear
vertex but is relevant to the full local quadratic action.
Set sigma=J(t) z, J=s/A, on the fixed interval -A<=z<=A, and
u=Jdot z/J. The pulled-back throat action is

```text
S=1/2 integral dt dz [(J/x)|fdot-u f_z|^2
                    -(x/J)|f_z|^2-a x J sec^2(Jz)|v* f|^2].   (1)
pi=(J/x)(fdot-u f_z),
H=1/2 integral dz [(x/J)(|pi|^2+|f_z|^2)
                   +a x J sec^2(Jz)|v* f|^2]
    +1/2 integral dz u {pi,f_z}.                            (2)
```

The braces denote symmetric ordering, summed over channels. The exterior
Hamiltonian is unchanged in this prescribed coordinate family. At the
static background let V0=a sec^2(z), K=integral(|pi|^2+|f_z|^2)/2,
and H_th=K+integral V0|v*f|^2/2. Expanding (2) gives

```text
H=H_0+phi V+phidot G+O(phi^2,phi phidot,phidot^2),
V=(1-kappa)K
  +1/2 integral V0[1+kappa+2kappa z tan(z)]|v*f|^2,
G=kappa/2 integral z {pi,f_z}.                             (3)
```

In particular, varying only the static potential misses G.

## 2. The Ward identity reduces the transition vertex to pressure

At the background, fdot=pi and pidot=f_zz-V0 vv* f. Define the inside
mouth pressure

```text
T_zz=(|pi|^2+|f_z|^2-V0|v*f|^2)/2.
```

Differentiate G and integrate once by parts in z:

```text
Gdot=kappa { [z T_zz]_-A^A -K
             +1/2 integral (V0+z V0')|v*f|^2 }.            (4)
```

Since z V0'=2z tan(z)V0, the effective vertex in the time-integrated
linear perturbation is exactly

```text
O=V-Gdot=H_th-kappa [z T_zz]_-A^A
        =H_th+r[T_zz(A)+T_zz(-A)].                         (5)
```

Time perturbations vanish at their initial and final endpoints. Equation
(5) is an identity for classical solutions and separated, on-shell
transition matrix elements. For a drive exp(-i w t), its matrix element
is V-i w G at a transition energy w. It does not identify all equal-time
contact terms after time ordering. Quantum anomaly, the second-order
proper-cut displacement, and the exterior gluing action still contribute
to the physical conservative response. The potential contribution in
T_zz is essential; it cannot be discarded at the deep cut.

The free vacuum supplies a separate sign check. Its ground energy is
-pi c/(6L), with L=D+2acos(epsilon x)/x. At x=1 its clock derivative is
-pi c(T+2r)/(6L^2), exactly the vacuum expectation of (5).

An independent local test uses any real solutions
f_i''=(V0-omega_i^2)f_i, i=1,2. Define

```text
E=omega_1+omega_2,
K12=integral(f_1' f_2'-omega_1 omega_2 f_1 f_2),
I12=integral z(omega_1 f_1 f_2'+omega_2 f_2 f_1'),
B12=[z(f_1' f_2'-(V0+omega_1 omega_2)f_1 f_2)]_-A^A.
```

Differentiating the bracket's integrand directly gives

```text
-K12+integral(V0+z V0')f_1 f_2+E I12=-B12.                (6)
```

Dividing the pair coefficients by 4sqrt(omega_1 omega_2) recovers (5).
This local identity does not presume that arbitrary test solutions are
eigenmodes of a globally matched return problem.

## 3. Exact free-channel weights for a specified smooth mouth

Now set a=0 and take one nonchiral CFT with equal chiral charges c, its
circle vacuum, and L=T+D. For a complex Dirac channel c=1. Here T_zz=T_00.
Convolve the distribution

```text
F=1_[-T/2,T/2]+r(delta_-T/2+delta_T/2)
```

with the periodic normalized Gaussian of variance h^2, h>0. This defines
an observable O_h, not a derivation of a four-dimensional mouth profile.
With k_n=2pi n/L and the same stress normalization as the
[canonical collar calculation](collar-conservative-response.md),

```text
(F_h)_n=2/L [sin(k_n T/2)/k_n+r cos(k_n T/2)] exp(-h^2 k_n^2/2),
W_n=2c(n^2-1)/(3L^2 n)
       [sin(k_n T/2)+r k_n cos(k_n T/2)]^2 exp(-h^2 k_n^2),
Delta chi_(r,h)(w)=sum_(n>=2) 2w^2 W_n/[k_n(w^2-k_n^2)].   (7)
```

The zero harmonic is a conserved energy and n=1 has zero vacuum weight.
For Im w>0 (7) is the complete statically-subtracted spectral response of
this observable. It is not the full metric influence action.
The form factor agrees independently with the finite-cut ray calculation:

```text
1-exp(iwT)-iw r[1+exp(iwT)]
 =-2i exp(iwT/2)[sin(wT/2)+r w cos(wT/2)].                 (8)
```

Thus a finite low-frequency ray response and a local spectral divergence
are compatible. The [transparent flux equation](transparent-clock-flux.md)
is a different, already specified observable with its own local balance
law; its finite-cut bound is not automatically a bound on (7).

## 4. A static subtraction does not remove the sharp pressure inertia

Fix T,D,r>0 and w in the upper half-plane with Re w>0. Let nu=2pi/L,
s_n=sin(k_nT/2), c_n=cos(k_nT/2), and E_h(k)=exp(-h^2 k^2).
Expanding (7), including the identically zero n=1 term, separates
the canonical interval response, a cross term, and an endpoint term:

```text
Cross_h=2c r w^2/(3pi L) sum (k^2-nu^2)sin(kT)E_h(k)/[k(w^2-k^2)],
End_h=2c r^2 w^2/(3pi L) sum (k^2-nu^2)c_n^2 E_h(k)/(w^2-k^2).
```

Use the exact decompositions

```text
(k^2-nu^2)/[k(w^2-k^2)]=-1/k+(w^2-nu^2)/[k(w^2-k^2)],
(k^2-nu^2)/(w^2-k^2)=-1+(w^2-nu^2)/(w^2-k^2).             (9)
```

The endpoint term is -w^2 C_h plus an absolutely convergent residual,
where

```text
C_h=2c r^2/(3pi L) S_h,
S_h=sum_(n>=1) cos^2(k_n T/2) exp(-h^2 k_n^2)
   =[Theta_h(0)+Theta_h(T)-2]/4,
Theta_h(t)=L/(2sqrt(pi)h) sum_(m in Z) exp[-(t-mL)^2/(4h^2)]. (10)
```

For fixed positive T,D, only the zero image in Theta_h(0) survives as
h->0, so S_h=L/(8sqrt(pi)h)-1/2+exponentially small terms.
The cross term's -sin(kT)/k part tends to the ordinary Fourier sawtooth;
its remaining summand is absolutely summable. The canonical interval's
subtracted response converges by the previously proved collar estimate.
Consequently a finite F(w;T,D,r) exists and

```text
Delta chi_(r,h)(w)=-c r^2 w^2/(12pi^(3/2)h)+F(w;T,D,r)+o(1). (11)
```

This is a proved obstruction to importing the canonical sharp-mouth
spectral kernel using static subtraction alone. The divergent term is
local in time and has the form of inertia. It is not a demonstrated
physical divergence or instability: the prescribed observable has not
included all local geometric and gluing variations. Those terms must be
derived from the same physical action, rather than chosen to cancel (11).
Nor is the fixed-D coefficient uniform as D/h->0: nearby mouth images
then overlap, which (10) retains exactly.

For reproducible evaluation of the finite part, standard circle sums
give, with U=wL/2 and Z=w(L/2-T),

```text
Q1=sum sin(kT)/k=(L-2T)/4,
Q2=sum sin(kT)/[k(w^2-k^2)]=[Q1-(L/4)sin(Z)/sin(U)]/w^2,
Q3=sum cos^2(kT/2)/(w^2-k^2)
  =L[cos(U)+cos(Z)]/[8w sin(U)]-1/(2w^2),
Delta chi_(0,0)=2c/(3pi L) [-nu^2 TD/8
  -(w^2-nu^2)L sin(wT/2)sin(wD/2)/(4w sin(U))],
F=Delta chi_(0,0)+2c r w^2[-Q1+(w^2-nu^2)Q2]/(3pi L)
  +c r^2 w^2/(3pi L)+2c r^2 w^2(w^2-nu^2)Q3/(3pi L).     (12)
```

These follow by integrating/differentiating the circle resolvent Fourier
series; the subtracted canonical sum uses sum sin^2(kT/2)/k^2=TD/8.

## 5. A uniform sufficient smoothing condition

Let T range over a fixed compact positive interval [T0,T1], D>0 be
arbitrary, 0<h<=T0/4, Re w in [m,M] with m>0, and 0<eta=Im w<=1.
There is a C depending only on T0,T1,m,M such that

```text
|Delta chi_(r,h)(w)-Delta chi_(0,0)(w)|
 <= C c [h+h^2/eta+r(1+eta^-1)
                  +r^2(h^-1+1+eta^-1)].                  (13)
```

Here r>=0 need not be tied to epsilon for the inequality.
The bound is uniform in every return length. Proof:

1. The canonical smoothing difference is bounded by equation (11) of
   the collar note, giving Cc(h+h^2/eta).
2. The Gaussian-damped sum sin(kT)/k is heat convolution of the periodic
   sawtooth. Its absolute value is at most L/4, since the normalized
   heat kernel is positive. In the residual cross sum use
   |sin(kT)|/k<=T1 below a fixed large frequency R and
   |w^2-k^2|>=eta Re w. The number of those modes is at most RL/(2pi).
   Above R>2|w| the summands are O(k^-3); the decreasing lattice tail
   divided by L is bounded uniformly. This gives Ccr(1+eta^-1).
3. Since cos^2<=1 and the Gaussian decreases,
   S_h<=L/(4sqrt(pi)h), hence C_h<=c r^2/(6pi^(3/2)h).
   The residual endpoint sum is bounded by splitting at the same R:
   counting below R gives CL/eta, and the k^-2 tail gives CL.
   This proves the last term in (13).

For q independent channels this bound adds without any upper bound on
D_j. At eta comparable to delta, choose h comparable to r comparable to
epsilon. The correction divided by q delta is O(epsilon/delta^2).
It vanishes for the retained epsilon=q^-2, delta=log^-2(q).
The comparison is to the canonical response at the same finite T;
neither a change of T nor a physical constraint equation is implicit.

By contrast, taking h to zero first at fixed epsilon fails by (11).
For example h=epsilon^4 makes r^2/h grow as epsilon^-2, at fixed positive
D and a fixed nonzero w. Static determinant control alone does not
select the safe joint limit. Equation (13) supplies a sufficient free
spectral condition for the next physical matching calculation.

## 6. Verification, predecessors and decision

Run `python checks/moving_mouth_clock_vertex.py`. The
[receipt](../evidence/moving-mouth-clock-vertex.json) binds this note and
its dependencies. Symbolic checks derive the pulled-back Hamiltonian,
linear coefficients and local Ward identity. Independent ODE integrals
test (6), including large mouth potentials. Spatial heat-kernel Fourier
integrals test (7); direct spectral sums test (9)-(12), refine Gaussian
tails, approach the divergent coefficient and finite part, and sample
the uniform smoothing estimate. Finite tests supplement the ordinary
proof; they do not establish its all-return uniformity on their own.

The [source supplement](../evidence/moving-mouth-source-review.json)
records the targeted search and exact cached primary sources.
[Jaekel--Reynaud, Motional Casimir force](https://arxiv.org/abs/quant-ph/0101069v1),
the motional-force section and appendix B, derives boundary-force response and its
relation to moving scattering. Their
[single-mirror paper](https://arxiv.org/abs/quant-ph/0101068v1), equation
(32) and appendix A, likewise connects displacement, pressure and
retarded response, with high-frequency transparency part of the model.
These are established mechanisms, published in 1992 and posted to arXiv
in 2001. Our mouths are transmitting joins, not those reflecting mirrors;
their response formulas are not imported as a wormhole calculation.
The CFT stress spectrum and canonical collar estimates are established
dependencies already reviewed in the linked notes.

**Decision:** keep physical long-return damping provisional, with matched
finite-duration transfer as fallback. The product Gaussian linear vertex
is now derived, but its interacting pressure response has not inherited
the canonical singlet bound: the charged throat Hamiltonian commutes
with an isolated charged throat, whereas its endpoint pressure does not.
The next decisive task is the covariant mouth/exterior matching, including
local inertia and the massive pressure response, together with corrected
background and angular/compact-state control. Neither a source-complete
novelty verdict nor a consequential physical endpoint is accepted.
No Lean/Lake or GCP run occurred; a supporting algebra identity does not
warrant a separate certification campaign before that endpoint exists.

Remaining to-do list: match the full local and nonlocal moving-mouth
response, resolve corrected geometry and coupled state/exterior, establish
the physical endpoint and novelty, independently verify, certify suitable
mathematics on GCP, and prepare the reproducible manuscript.
