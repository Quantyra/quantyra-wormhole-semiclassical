# A joint weak-coupling window for the response investigation

2026-10-10, S051/S052/S053. **Narrow the next physical-response calculation
to a stated joint limit**, rather than use small four-dimensional coupling
alone to discard the nonspherical photons. Two ordinary results support
this decision: an exact norm for the projected angular exchange kernel,
and a coupling-independent bound for the prescribed dipole electric
potential, including its accumulated effect on the proposed damping time.
Neither result closes the interacting quantum stress response.

The provisional primary remains physical long-return damping of the
symmetric throat clock. The fallback remains finite-duration, band-resolved
energy transfer with physical matching. This is not a new matter sector or
an accepted stability/novelty result. The limit varies a renormalized
coupling as well as flux; it is **not** a large-flux theorem in one fixed
theory. Earlier fixed-coupling results retain their stated scope.

## 1. Exact angular sum and operator norm

Use the sphere LLL projector and its established form factor from the
[finite-flux audit](finite-flux-response-audit.md). The angular Hilbert
space has dimension q and spin (q-1)/2. Let Y_la be a real orthonormal
basis of degree-l harmonics, a=1,...,2l+1, and set

```
T_la = P_q Y_la P_q,                 T_la*=T_la,
beta_l(q) = (q-1)! q! / [(q-1-l)! (q+l)!],  0<=l<q,
beta_l(q) = 0,                     l>=q.
```

The trace identity is Tr(T_la^2)=q beta_l/(4 pi). Rotational covariance
and irreducibility imply

```
sum_a T_la^2 = (2l+1) beta_l/(4 pi) I_q.                 (1)
```

Write H_q=sum_(k=1)^q 1/k. Then, including q=1,

```
S_q = sum_(l=1)^(q-1) (2l+1) beta_l/[l(l+1)] = H_q-1,
sum_(l>=1,a) T_la^2/[l(l+1)] = c_q I_q,
c_q = (H_q-1)/(4 pi).                                    (2)
```

**Proof of the sum.** Legendre expansion of ((1+x)/2)^(q-1), evaluated
at x=1, gives sum_l(2l+1) beta_l(q)=q. The factorial ratios give

```
beta_l(q+1)-beta_l(q) = l(l+1) beta_l(q+1)/[q(q+1)],
                                                       1<=l<=q.
```

Consequently S_(q+1)-S_q=1/(q+1). Since S_1=0, (2) follows. This is a
finite identity, without a large-q interchange of sums.

On H_q tensor its conjugate space define

```
K_q = sum_(l>=1,a) T_la tensor conjugate(T_la) / [l(l+1)].
```

The row/column operator Cauchy-Schwarz inequality and (2) give
||K_q||<=c_q. For the normalized vectorized identity v_q,
K_q v_q=c_q v_q, since vectorization turns each summand into T_la^2.
Thus **||K_q||=c_q**. The unconjugated two-particle kernel obeys the
same upper bound; its norm need not equal c_q. Any N-particle sum of
such angular pair kernels has norm at most N(N-1)c_q/2 by the triangle
inequality. No independent boson has been assigned to each generator.
The matrices for different angular harmonics generally do not commute.

Here l=0 is excluded: its dynamical Schwinger sector was retained in the
[uniform singlet calculation](singlet-stress-continuum.md). The factor
1/[l(l+1)] is the inverse angular photon mass operator at zero Euclidean
longitudinal momentum. A common additional nonnegative scalar in these
denominators preserves the upper bound. This does not permit replacing
the Lorentzian photon propagator near its poles by an instantaneous
interaction, or bounding the renormalized many-body stress by an angular
matrix norm. Radial dependence, current normalization, the filled sea,
vertex corrections and regulator removal still need their own argument.

For an angular pair Hamiltonian with coefficient g^2 and N=O(q), the
result motivates the conservative per-particle proxy

```
eta_ang = g^2 q (H_q-1)/(4 pi) = O(g^2 q log q).            (3)
```

This is a precisely identified proxy, not a proved QED response error.
At fixed g it is not small uniformly as q grows. The free-current
polarization parameter g^2 q beta_l/(4 pi^2) at small l makes the same
basic distinction between small g^2 and small collective response.

## 2. The dipole potential remains bounded when screening is weak

Keep exactly the point-dipole, local-LLL equation of the
[screening baseline](rotating-electric-screening.md), with d>0,

```
r_+ = |x-d ez/2|,       r_- = |x+d ez/2|,
B = q/2 [(x-d ez/2)/r_+^3-(x+d ez/2)/r_-^3],
V = g^2 |B|/(2 pi^2),
L A_t := (-Delta+V) A_t = -2 Omega ey.B.                   (4)
```

Use the unique D^(1,2)(R^3) solution, with its L6 representative and no
added constant. The earlier proof establishes this solution for every
g>0 by Hardy, Sobolev and the coercive variational form. For g=0 use
the same variational argument for -Delta. Define the explicit unscreened
solution

```
A_0 = -Omega q/2 (y/r_+ - y/r_-).                          (5)
```

Indeed -Delta(y/r_i)=2y/r_i^3. There is no point delta term: the
integrated flux of its locally O(1/r_i) gradient through a small sphere
vanishes. At infinity A_0=O(1/r), grad A_0=O(1/r^2), and near either
charge its gradient is O(1/r_i). Thus A_0 belongs to D^(1,2).

For Omega>=0, on Q={y>0,z>0}, r_+<r_-, hence A_0<=0 and the source
in (4) is nonpositive. V is invariant under reflection in either y=0
or z=0, while that source changes sign. Uniqueness makes A_t odd under
both reflections. Both A_t and A_0 have zero Sobolev trace on the two
planes. Their restrictions lie in D_0^(1,2)(Q); the isolated charge on
the boundary has zero capacity and causes no extra trace datum.

Test the weak equation by (A_t)_+ to obtain A_t<=0 on Q. Next
L(A_0-A_t)=V A_0<=0; testing with (A_0-A_t)_+ gives A_0<=A_t.
The tests are legitimate Lipschitz truncations in D_0^(1,2), with V
controlled by the translated Hardy inequality. Reflect into the other
quadrants and reverse signs for Omega<0. Consequently

```
|A_t(x)| <= |A_0(x)| <= |Omega| q/2,
||A_t||_infinity <= min(|Omega|q/2, 4 pi^2|Omega|/g^2).     (6)
```

The two terms y/r_+ and y/r_- have the same sign and absolute values
at most one, so their difference is at most one, not two. The constant
in the unscreened bound is sharp: approach either charge along ey.
The second bound in (6) is the earlier maximum-principle result.
Both hold for arbitrary d>0. No strong-screening expansion was used.

For any self-adjoint reference Dirac Hamiltonian to which this bounded
real multiplication term is added, Duhamel gives

```
||exp[-it(H_0-A_t)]-exp[-it H_0]|| <= |t| |Omega| q/2.     (7)
```

This controls the electric term in that specified operator on the entire
time interval. It does not construct the reference operator for the full
rotating spacetime. Finite-mouth electric data, finite-flux nonlocal
screening, higher rotation orders and the quantum state are not fixed by
(4). The equation on all of space remains an elliptic baseline, not a
global corotating vacuum.

## 3. A compatible window for the quantities actually controlled

Retain the leading MMP scales, with g evaluated at the magnetic scale:

```
r_e^2=pi G q^2/g^2,        ell=16 r_e^3/(q G),
ell/r_e=16 pi q/g^2,       delta=d/ell,
Omega=2 sqrt(r_e/d^3),     t_d=ell/delta.
```

t_d is the time scale suggested by the provisional model width; no
physical decay law is assumed. Set delta=q^(-a), g^2=q^(-b). Then
d/r_e=16 pi q^(1+b-a). The following table gives powers of q, omitting
fixed constants and the displayed logarithms. The two propagation/collar
entries remain nominal requirements from earlier work, not newly proved
physical errors.

| Quantity | Power | Meaning |
| --- | --- | --- |
| 1/(q delta^5) | -1+5a | proved Gaussian singlet absorption bound |
| delta log(2+q delta^2) | -a, times log q | proved soft-pair bound |
| log(2+q delta)/(q delta^2) | -1+2a, times log q | finite-channel quadrature bound |
| 1/[delta sqrt(q delta)] | -1/2+3a/2 | hypothetical accumulated full-return bound |
| sqrt(r_e ell)/d | -(1+b)/2+a | nominal collar error divided by delta |
| Omega t_d | -(1+b)/2+5a/2 | orbital change over t_d |
| q Omega t_d/2 | (1-b)/2+5a/2 | proved electric-term bound (7) |
| eta_ang/delta | 1-b+a, times log q | angular interaction proxy only |
| (M ell)^2/delta | 1-b+a | central throat singlet mass parameter |

Here M ell=g sqrt(q)/(2 pi) refers to the central optical mass, not
the maximum of the mass profile near a mouth. A sufficient open window for every listed
quantity to vanish is

```
0<a<1/5,               b>1+5a.                            (8)
```

The remaining Gaussian discretization term 1/(q delta) and the frequency
smoothing error delta^2 also vanish. Take, for example, a=1/16 and b=2.
Then d/r_e scales as q^(47/16), the electric bound as q^(-11/32), and
eta_ang/delta as q^(-15/16) log q. Choosing the collar at
epsilon=sqrt(r_e ell) makes both nominal relative collar errors scale
as q^(-23/16). The Gaussian rank bound is q^(-11/16).

This limit has g^2 q -> 0. The local approximation A_t approximately
-(4 pi^2 Omega/g^2) ey.B/|B|, which required g^2 q -> infinity, must
**not** be imported into it. Equation (6) replaces that approximation
for the bounded electric term. Small orbital change on the causal size
t_d does not produce a state beyond the light cylinder. Nor does (8)
prove an interacting stress estimate: possible small Lorentzian
denominators must be bounded before assigning their error a power of q.

## 4. Mass generation is an established issue, with a limited comparison

[MMP v3](https://arxiv.org/html/1807.04726v3), Appendix A, discusses the
Schwinger singlet, the residual massless sector and massive KK vectors.
Its leading throat construction and no-shortcut conclusion are
predecessors. We do not claim that photon dynamics was omitted there.

[Gusynin, Miransky and Shovkovy, hep-ph/9908320v1](https://arxiv.org/abs/hep-ph/9908320v1),
equations (2), (47)-(49), analyzes constant-field flat QED using a
specified improved Schwinger-Dyson truncation. For N physical fermion
flavors it obtains

```
m_dyn approximately C sqrt(|eB|) (N alpha)^(1/3)
          exp[-pi/{alpha log(C_1/(N alpha))}],
alpha=g^2/(4 pi),             C_1 approximately 1.82.
```

The paper's argument for the approximation is gauge-specific; its
Appendix C checks the Ward identity to the retained infrared order.
This is not a rigorous curved-wormhole gap estimate. Here N=1, not q:
flux degeneracy is not the number of four-dimensional flavors.
[Maldacena, 2004.06084v2](https://arxiv.org/html/2004.06084v2), section 3
and reference [21], already calls out this smaller nonperturbative mass
in the magnetic-black-hole discussion. The citation chain therefore
rules out presenting the concern itself as a new discovery.

As a comparison only, insert |eB|=q/(2r_e^2) into the flat formula.
The throat infrared diagnostic is m_dyn r_e, also the central optical
mass times ell, because the central lapse is r_e/ell. Using m_dyn ell
without this redshift would compare different clocks. At g^2=q^(-b),

```
log[(m_dyn r_e)/delta]
 = (1/2+a-b/3) log q
   - 4 pi^2 q^b/[b log q+log(4 pi C_1)] + O(1).             (9)
```

It tends to minus infinity faster than any multiple of log q for b>0.
This is an algebraic consequence of the *flat* formula, not an upper
bound on the curved theory's gap. In particular the planar infrared
approximation can fail once its correlation length exceeds r_e. A
formal fixed-g extrapolation is not uniform in q; running at a fixed
reference coupling is a different limit again. None of these observations
establishes loss of wormhole support at a finite physical charge.

## 5. Verification, novelty boundary and next decisive test

Run `python checks/weak_coupling_window.py`. The script checks the finite
sum with exact rationals; independently constructs T_la from sphere
quadrature; checks the sum of squares, conjugate kernel norm and
noncommutation; verifies the unscreened Poisson identity symbolically
and by finite differences; tests discrete comparison on bounded boxes;
and derives the scale powers from the dimensional definitions. Results
and exact source hashes are in
[the receipt](../evidence/weak-coupling-window.json). The box checks are
discretization checks, not a substitute for the whole-space proof. The
mass-formula evaluations test only the stated comparison (9).

The [source supplement](../evidence/weak-coupling-source-review.json)
records ten search queries and scoped primary
reads. Broad query misses do not establish novelty. The central source
register's preexisting edits are preserved; this supplement does not
claim completion of the full proof/citation review. The finite projection
and standard comparison principles are reused mathematics. The physical
damping or pulse-transfer result, if obtained, still needs a direct
comparison with the closest prior work. No Lean/Lake invocation occurred.

**Decision:** continue in the explicit joint window (8) as a candidate
physical approximation regime. The next decisive calculation must derive
the renormalized real response and a bound on the interacting stress at
width delta, with the actual clock vertex and a finite-time state. The
exact angular norm alone is insufficient for this step. A failure to
close it is not a stability or no-go theorem; reassess the primary and
fallback before adding a further formalization campaign.

Remaining to-do list: control the renormalized physical stress response,
complete endpoint/throat matching and the finite-time state, decide
primary/fallback novelty, obtain the consequential result, independently
verify and GCP-certify suitable endpoints, and prepare the manuscript.
