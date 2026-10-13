# Electric thermodynamic matching and finite optical cuts

S072, 2026-10-13. **The classical fixed-temperature RN charge weight
matches the capacitance at the actual horizon radius, without an extra
extremal-entropy factor. The radial electric midpoint bound also admits
a finite-cut estimate uniform in the number of fermion channels.**
These remove two specific ambiguities in the preceding state-matching
target. They do not identify the full interacting black-hole density
operator or its charged matrix elements, and do not compute driven stress.

## Thermodynamic charge is the total asymptotic charge

Use G, AdS radius L_A, fixed magnetic geometric charge P, Maxwell coupling
g and integer total electric charge Q. Write alpha=G g^2/(4pi). The
classical dyonic RN family is

```
f(r)=1-2GM/r+(P^2+alpha Q^2)/r^2+r^2/L_A^2,
M(r,Q)=r/(2G)[1+r^2/L_A^2+(P^2+alpha Q^2)/r^2],
S(r)=pi r^2/G,
T(r,Q)=[1+3r^2/L_A^2-(P^2+alpha Q^2)/r^2]/(4pi r).       (1)
```

The primary thermodynamic and ensemble conventions are reviewed in
[Iliesiu and Turiaci, section 2.1](https://arxiv.org/html/2003.02860v3).
Electric charge is varied here at fixed magnetic P; their expansion
along a purely electric family is not copied. Their sections 3.1--3.2
also distinguish the classical free energy from quantum prefactors.
Their partition function does not include the present light charged
fermion coupling; see their section 5.1. It cannot be appended to the
already gauged radial fermions as an independent factor.

Fix T>0 and let R=r(0) be on the large stable branch. Put w=Q^2 and
C_T(r)=1+6r^2/L_A^2-6pi T r. Assume C_T(R)>0 and
R>=pi T L_A^2/2, so C_T is nondecreasing for r>=R. The equation and
its charge derivative are

```
3r^4/L_A^2-4pi T r^3+r^2-P^2-alpha w=0,
dr/dw=alpha/[2r C_T(r)].                                  (2)
```

The canonical classical free energy is F(T,Q)=M(r(Q),Q)-T S(r(Q)).
On shell M_r=T S_r. All implicit radius/entropy terms therefore cancel
in the first charge derivative:

```
dF/dw = g^2/(8pi r(w)),
F(T,Q)-F(T,0)=g^2 Q^2/(8pi R)+O(Q^4),
I_bare(T)=4pi R/g^2.                                      (3)
```

This is the full classical on-shell free energy, not just the mass.
It equals the electrostatic response coefficient at that temperature.
Let b be the extremal magnetic radius, C_b=1+6b^2/L_A^2. Expanding
R=b+2pi b^2 T/C_b+O(T^2) gives

```
beta g^2/(8pi R)=beta g^2/(8pi b)-g^2/(4 C_b)+O(T),
partial_(Q^2) S_ext = g^2/(4 C_b).                        (4)
```

Thus the familiar extremal-entropy term is recovered if one expands
about b. Multiplying exp[-beta Q^2/(2 I_bare(T))] by that entropy
factor a second time is inconsistent. This observation does not revoke
the earlier entropy calculation, which used the extremal coefficient.
It fixes how to combine the two at a common nonzero temperature.

## A controlled charge window, with a tail bound

The quadratic approximation can be bounded without subtracting large
black-hole masses numerically. Define

```
a=beta g^2/(8pi R), d=alpha/[2R^2 C_T(R)], p=a/d=beta R C_T(R)/G,
I(Q)=beta [F(T,Q)-F(T,0)].
```

Equation (2) implies 0<=r(w)-R<=R d w. Integrating (3) gives

```
(a/d) log(1+d Q^2) <= I(Q) <= a Q^2,
0<=a Q^2-I(Q)<=a d Q^4/2.                                (5)
```

The second bound follows by integrating
0<=1-R/r(u)<=d u. It is valid at every charge on this large branch;
no near-extremal extrapolation of the Schwarzian formula is used for
the far tail. The exact classical RN family is the comparator there.

For 0<a<=1 and p>=100, choose Q_*=ceil(sqrt(80/a)). On |Q|<=Q_*,
Q_*^2<100/a, hence the log-weight error is at most 5000/p.
For the classical probability proportional to exp[-I(Q)], the tail
outside this window is less than exp(-39). Here is an explicit proof.
The normalization is at least sum_Q exp(-a Q^2)>=sqrt(pi/a).
For |Q|<=1/sqrt(d), the first bound in (5) gives I(Q)>=a Q^2/2.
For |Q|>=1/sqrt(d), write u=sqrt(d)|Q|>=1 and use
(1+u^2)^(-p)<=(2u)^(-p). Monotone sum/integral bounds therefore give

```
Prob(|Q|>Q_*) <= 2 sqrt(a/pi) exp(-a Q_*^2/2)
 +sqrt(2) erfc(sqrt(a/2) Q_*)
 +2^(1-p)/sqrt(pi) [sqrt(a)+sqrt(p)/(p-1)] < exp(-39).     (6)
```

The estimate deliberately overcounts the two ranges. For the final
inequality use Q_*>=sqrt(80/a), erfc(v)<=exp(-v^2)/(sqrt(pi)v)
and p>=100. The same tail bound covers the Gaussian comparator.
The total-variation distance of these two *classical charge
probabilities* is at most exp(5000/p)-1+2 exp(-39), by conditioning
on the common window and comparing its positive weights. This is not
an operator-norm or fermion-correlation bound.

For the candidate g=k^-4, L_A/l_p=k^5, P^2=pi G k^16,
X=b^2/L_A^2, X(1+3X)=pi k^6, N=2k^4, eta=1/k and
T=N eta z0/(2pi phi_r), phi_r=b^3/(G C_b), the earlier temperature
bound gives R>=b, C_T(R)>=C_b/2 and a<1/(12k).
Since X>=1/3 and X^2>=pi k^6/6,

```
p >= pi^2 k^23/(12 z0) >= pi^2 k^23/[12(7/25)].           (7)
```

For clarity, C_T(R)>=C_b/2 follows directly: the temperature equation
has its root between b and 2b, because 4pi T b/C_b is much smaller
than one in this family; on that interval
6pi T R<=12pi T b=6 N eta z0 G C_b/b^2<C_b/2. The same inequality
ensures R>=pi T L_A^2/2. All these comparisons are analytic for
k>=2000, and the checker also substitutes the endpoints.
At k=2000, the on-window classical nonquadratic log-weight error
5000/p is below 2.04e-73. The deliberately conservative tail bound,
not this tiny central error alone, controls normalized probabilities.

Quantum determinants, angular branches, matter-dependent entropy and
other gravitational fluctuations are not included in (5)--(7). In
particular a small classical distribution error cannot be substituted
for a bound on the matrix elements of the charged boundary field.

## What the screened susceptibility matches

For any consistently defined Euclidean generating functional at this
fixed boundary metric, let its grand potential have expansion
Omega(mu)=Omega(0)-I_stat mu^2/2+O(mu^4). Legendre transformation gives
F(Q)=F(0)+Q^2/(2 I_stat)+O(Q^4). This is a statement about the
thermodynamic Hessian. When the charged matter is retained in the same
functional, its polarization belongs in I_stat once.

The prior radial reduction gives I_0=lambda/K(0) from that projected
electric action. Equations (3)--(4) show which bare gravitational
reference this coefficient must replace at the actual temperature.
They do not license multiplying the full RN charge partition by the
radial charge partition, or interpreting a full asymptotic-charge weight
as an independent horizon-cut weight. At a cut, Gauss law instead reads
Q_infinity=h+Q_exterior. The two charges coincide only when there is no
exterior charged matter. This distinction is also explicit in
[Harlow, section 2](https://arxiv.org/html/1510.07911v2).

An exact total-charge probability is insufficient to fix the TFD charged
matrix elements. The physical initial state still needs a joint
gravity/matter construction with the same boundary conditions. Fixing
each side's charge separately kills a charge-one two-sided expectation;
fixing only the pair's total charge permits a charge-correlated TFD.
The latter is the intended preparation target, not an assumption that
an arbitrary entangled pair has the required correlation.

## Finite optical cutoff: a uniform radial midpoint bound

The preceding [radial correlator](radial-electric-correlator.md) used
the optical-volume limit at fixed species count to simplify its theta
function. That limit is unnecessary for a useful bound. Retain its
neutral-inner-cut ensemble on any finite optical interval of length L_o,
with n channels, lambda=n/pi, nonnegative V and the original integer
species lattice. Define

```
c=pi beta/(2L_o), x=beta K_L(0)/(2lambda), y=x-c/n>=0,
A=c 1_n+y 11^T,
D=D_L(1/2), 0<=D<=y/4.
```

The radial formula is equivalently

```
G_V,L(beta/2)/G_0,L(beta/2)
 = [C_A,m/C_e(c)] exp(y/4-D),                            (8)
```

where C_A,m is its actual species-charge trace. The square bracket
is *exactly* the finite-circle constant-charging ratio already bounded
in [electric-reservoir-midpoint.md](electric-reservoir-midpoint.md).
That theorem retains the complete fermion sea and radial occupations;
it is uniform in channel count and circumference. Applying it with
charging parameter y, then using 0<=D<=y/4, proves

```
exp(-y/2) <= G_V,L/G_0,L <= exp(y/2)/cosh(y/2).           (9)
```

This reuse does not replace the radial interaction by y Q^2: its
nontrivial oscillator factor remains in (8). The comparison theorem
is applied only to the identical charge trace.

For a spherical exterior beginning at areal radius r_c, the zero-
frequency trial z=(x+L_o)/L_o gives

```
K_L(0)-1/L_o <= integral V dx
 =lambda integral_(r_c)^infinity dr/kappa
 =lambda g^2/(4pi r_c),
y<=beta g^2/(8pi r_c).                                  (10)
```

Since r_c>=R>=b, the candidate has y<1/(12k), uniformly in n and
L_o. Thus the relative loss against its *same finite-cut free
reference* is below 1-exp[-1/(24k)], or 2.08332e-5 at k=2000.
No artificial reflecting cut has thereby been identified with a
transparent physical horizon. The finite reference G_0,L must also
be matched to the free JT preparation before using this in capture.
The stronger infinite-volume result remains valid in its earlier scope.

## Verification and next action

[Checker](../checks/electric_thermodynamic_matching.py),
[receipt](../evidence/electric-thermodynamic-matching.json),
and [source record](../evidence/electric-thermodynamic-source-review.json).
The checks compare the on-shell mass/entropy difference with an
independent charge integral, test the classical probability tails and
window error, and verify the finite-cut formula on both constant and
nonconstant radial potentials. High-precision arithmetic retains the
very small candidate radius and action corrections without catastrophic
mass subtraction. The ordinary inequalities establish the entire
parameter family; numerical examples are supplementary.

The initial numerical run failed when an integration node so close to
zero charge made the secant solver's two starting points coincide at
working precision. Newton iteration with the analytic positive derivative
resolved that numerical failure. The final candidate calculation uses
400 decimal digits to resolve its tiny remainder after subtraction of
large free energies. The original source and failure log remain retained;
the physical formulas and acceptance inequalities were not changed.

Continue with one joint finite-temperature gravity/matter functional
and the actual angular/coset branching. Use (3) as the bare matching
condition, not a second charge weight; retain the boundary ladder and
initial covariance, and derive the driven stress with one subtraction.
The unresolved physical horizon/gluing and quantum-sector conditions
remain explicit. This is progress on matching, not completion of the
four-dimensional state or stable-wormhole goal.

Remaining to-do list: joint gravitational/angular initial-state and
operator matching; driven covariance and complete stress with massive
modes and geometry; preparation, quantum fluctuations, sufficient
general stability and finite-energy transport.
