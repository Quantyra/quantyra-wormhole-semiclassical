# The lapse collar has a large conservative response

2026-10-10, S051-S053. **The small-collar absorption estimate does not
control the clock stiffness.** For the specified CFT Hamiltonian family,
the exact ground energy contains a negative spatial-gradient term. It
grows as q/h when the optical collar width h shrinks. At the same time,
the complete statically subtracted response has a uniform sharp-collar
limit at damping-scale frequency resolution. These are compatible facts.

The consequence is an explicit failure of the proposed inference from
small collar matrix-element errors to a physical resonance. The closed
static-plus-dynamic model in [the clock note](collective-clock-response.md)
remains a defined model; its assumed static completion is not derived by
taking the sharp limit of the canonical lapse Hamiltonian. The full MMP
geometry has not been shown unstable. The physical primary and fallback
remain open, and this established CFT mechanism is not accepted novelty.

## 1. Exact static energy in the canonical CFT family

Take a unitary nonchiral CFT with equal chiral central charges c and an
NS vacuum on a circle of length L. A complex Dirac field has c=1; its
total vacuum energy is -pi c/(6L). Normal order the stress once in this
fixed circle theory, retaining that Casimir constant. For a positive
smooth periodic function v(s), define

```
H[v] = integral_0^L v(s) [T_+(s)+T_-(s)] ds,
J[v] = integral_0^L ds/v(s),       v_bar=L/J[v].
```

This specifies the operator, including its scalar part. Its exact lowest
energy in the vacuum sector is

```
E[v] = -pi c/(6 J[v]) - c/(24 pi) integral_0^L v'(s)^2/v(s) ds.   (1)
```

This is established inhomogeneous-CFT theory. For clarity, derive it from
the standard stress transformation. Set y'(s)=v_bar/v(s), so y(s+L)=y(s)+L.
The two chiral stress operators transform as

```
U T_+/- (s) U^-1 = y'(s)^2 T_+/- (y(s)) - c {y,s}/(24 pi).
```

Hence U H[v] U^-1=v_bar H[1]-(c/(12 pi)) integral v{y,s} ds.
Direct differentiation gives
`v{y,s}=-v''+v'^2/(2v)`. The v'' integral vanishes by periodicity.
Unitary equivalence and the circle ground energy give (1).
The required representation and this Hamiltonian transformation are
[Moosavi v3](https://arxiv.org/html/1912.04821v3), sections 2.1 and 5.1,
Lemma 5.3 and its proof, especially (5.19)-(5.20). This is not a newly
derived CFT vacuum mechanism.

At fixed positive collar width every term in (1) is finite. The issue
below is a singular family of spatial profiles, not an uncanceled
ultraviolet divergence of each smooth Hamiltonian. Changing its scalar
part changes the defined family. It cannot be called an arbitrary
renormalization after the canonical operator has been fixed.

## 2. Contact and susceptibility agree with the static energy

Let f be real, smooth and periodic and define

```
F=integral f,   Q=integral f^2,   I=integral f'^2,
f_n=L^-1 integral f(s) exp(-2 pi i n s/L) ds,
O_f=integral f(s) T_00(s) ds.
```

For the affine perturbation v=1+a f, differentiating (1) gives

```
E_a'(0) = -pi c F/(6 L^2),
E_a''(0) = -c I/(12 pi) + pi c [Q-F^2/L]/(3 L^2).          (2)
```

Independently, the already reproduced two-chirality stress weights are

```
W_n = 2 pi^2 c n(n^2-1)|f_n|^2/(3 L^2),
E_n = 2 pi n/L,                          n>=2.
```

With the retarded convention of [the channel calculation](long-return-response.md),
`chi_f(w)=sum W_n[(w-E_n)^-1-(w+E_n)^-1]`. Thus

```
chi_f(0)=-2 sum W_n/E_n
        =-c I/(12 pi)+pi c [Q-F^2/L]/(3 L^2).               (3)
```

Here Parseval gives `I=(8 pi^2/L) sum_(n>=1) n^2|f_n|^2`
and `Q-F^2/L=2L sum_(n>=1)|f_n|^2`. This proves (2) by an
independent spectral route, including the absent first harmonic.

The actual clock family previously specified is
`v_phi=1+(exp(phi)-1)f`, so H''(0)=O_f. Its matter stiffness is

```
K_m = <O_f>+chi_f(0)
    = -pi c F/(6 L^2)
      +pi c [Q-F^2/L]/(3 L^2)-c I/(12 pi).                 (4)
```

The first term is the contact term. It does not cancel the last term.
A uniform f has I=0 and no spectral response. A first-harmonic affine
profile has chi_f(0)=0: its gradient and variance terms cancel exactly.
For example, v=1+a cos(2 pi s/L), |a|<1, has E[v]=-pi c/(6L).
These checks fix both the factor of two and the Schwarzian sign.

## 3. Gaussian collars: a uniform negative term

Use precisely the periodic Gaussian smoothing from the earlier drive:
f_h is the convolution of the interval indicator 1_[0,T] with a
normalized Gaussian of variance h^2, on L=T+D. Then F=T, 0<=f_h<=1,
and for n!=0,

```
(f_h)_n=(1-exp(-i k_n T))/(i k_n L) exp(-h^2 k_n^2/2),
k_n=2 pi n/L.
```

Parseval followed by Poisson summation gives the exact identity

```
I_h = 1/(sqrt(pi) h) sum_(m in Z)
      [exp(-(mL)^2/(4h^2))-exp(-(T+mL)^2/(4h^2))].         (5)
```

For a_0=min(T,D) and h<=a_0/2,

```
|sqrt(pi) h I_h-1| <= 8 exp[-a_0^2/(4h^2)].                (6)
```

To see this, bound the nonzero first-image terms by
`2 exp[-L^2/(4h^2)]/[1-exp(-L^2/(4h^2))]`. The second-image sum is
at most the same geometric denominator times
`exp[-T^2/(4h^2)]+exp[-D^2/(4h^2)]`. Use L>=2a_0 and h<=a_0/2.
The constants are uniform over arbitrarily long D.

Consequently, for q channels with D_j>=d, h/min(T,d)->0,

```
sum_j K_m,j = -q c/(12 pi^(3/2) h) + O(q c/T)
              + O(q c exp[-min(T,d)^2/(4h^2)]/h).          (7)
```

Indeed the contact is bounded by pi c/(6T), and the variance term is
bounded by pi c/(12T), since Q<=T and T D/L^2<=1/4. Thus (7) needs
no finite maximum return length. For compact transition collars an even
more elementary obstruction is `I>=1/h_left+1/h_right`, by
Cauchy-Schwarz on each unit-height transition. Small L1 error from the
interval profile places no upper bound on I.

The leading clock stiffness used in the closed MMP reduction is
`K_0=q/(8 ell)`, with c=1 and T=pi ell. The magnitude of the leading
extra term in (7), relative to K_0, is

```
2 ell/(3 pi^(3/2) h).                                     (8)
```

It grows when h/ell shrinks. Weak gauge coupling cannot suppress this
free-CFT contribution. It is present even when all gauge interactions
are removed. This does not identify the canonical collar with the
metric perturbation solving the four-dimensional constraints.

## 4. A complete diagnostic: this lapse-only completion is unstable

For a finite number of smooth channels, let the clock have positive
kinetic mass M and a remaining potential U(phi). Fix its first derivative
so phi=0 is stationary; this can be done by a linear term, without
altering the second derivative. Assume `U''(0)=O(q/ell)` uniformly in h.
For sufficiently small h/min(T,d), (7) makes

```
K_stat=U''(0)+sum_j K_m,j < 0.
```

The exact retarded *linear-response* inverse propagator in this model is

```
K(w)=-M w^2+K_stat+sum_j [chi_j(w)-chi_j(0)].               (9)
```

For gamma>=0 it obeys

```
K(i gamma)=M gamma^2+K_stat
    +2 gamma^2 sum_(j,n) W_jn/[E_jn(E_jn^2+gamma^2)].      (10)
```

The smooth profiles make these sums convergent. This continuous real
function starts negative, tends to positive infinity, and is strictly
increasing for gamma>0: each summand has positive derivative. There is
therefore exactly one positive gamma with K(i gamma)=0. Its derivative
is nonzero. The retarded clock has an exponentially growing linear mode
in this defined completion. Positive spectral absorption does not prevent it.

A concrete diagnostic uses ell=1, T=pi, c=1, q identical returns D,
M=q/8 and
`U(phi)=(q/16)exp(2phi)+(q/24)exp(phi)+lambda phi`.
The first two terms are the leading gravity and throat-anomaly potential
scales used before; lambda only enforces stationarity at the chosen
finite D. In particular U''(0)=7q/24. The numerical check locates the
unique positive root of (10) at three collar widths. This is a
well-defined lapse-only diagnostic, **not** the complete MMP perturbation
problem. It does not prove instability, collapse or loss of traversability
for an actual wormhole.

## 5. The subtracted dynamic kernel still converges

Let Delta chi_h(w)=chi_h(w)-chi_h(0), and let Delta chi_0 be the sharp,
statically subtracted kernel in the clock note. For Re w>0, Im w=eta>0,

```
|Delta chi_h(w)-Delta chi_0(w)|
 <= 8 c |w|^2 h/(9 pi^2)
    +8 c |w|^5 h^2/[3 pi^2 eta Re w].                    (11)
```

**Proof.** The sharp weight satisfies `W_n<=c E_n/(3 pi L)` and the
smooth weight equals it times exp(-h^2 E_n^2). From the subtracted
spectral formula, the difference is bounded by

```
2 c |w|^2/(3 pi L) sum_(n>=2)
       (1-exp(-h^2 E_n^2))/|E_n^2-w^2|.
```

For E_n>=2|w| the denominator is at least 3E_n^2/4. The decreasing
function min(h^2,E^-2) has integral 2h, so its lattice sum is at most
Lh/pi. This gives the first term in (11). For E_n<2|w|, there are at
most |w|L/pi terms. Their numerator is at most 4h^2|w|^2 and their
denominator is at least eta Re w, giving the second term. No bound on
the largest return length entered. The proof also justifies the sharp
subtracted sum's convergence at each such w.

At |w|~1/ell, eta~delta/ell, summing q channels and dividing by the
model damping-kernel scale q delta/ell makes (11)

```
O(h/(ell delta) + h^2/(ell^2 delta^2)).                    (12)
```

It vanishes if h/ell=o(delta), even while the static ratio (8) diverges.
For example retain `delta=q^(-1/16)`, `g^2=q^(-2)` and choose
`h/ell=q^(-1/4)`. This respects h<<d, r_e<<h and even the nominal
gradient condition 1/h<<sqrt(q)/ell. It makes (12) small but (8) large.
Thus the problem is not repaired by the
[joint weak-coupling window](weak-coupling-window.md).

This theorem concerns the free-channel susceptibility. It does not
extend the positive-band singlet theorem to the full interacting real
response, or control a Lorentzian photon pole.

## 6. Equal optical propagation does not fix the conservative completion

The scalar functional

```
C_h(phi)=c/(24 pi) integral v_phi'(s)^2/v_phi(s) ds          (13)
```

would exactly cancel the spatial-gradient piece of (1), and would add
`c I_h/(12 pi)` to the stiffness while leaving the positive-frequency
transition weights unchanged. Equation (13) states the size of the
missing completion; it does **not** establish that it is an admissible
four-dimensional counterterm or that the matched geometry supplies it.
The total geometric vertex can also contain spatial-metric, shift,
moving-endpoint and state contributions absent from H[v_phi].

There is an explicit two-dimensional comparison that makes this issue
concrete. On the same circle prescribe two metric families,

```
g_A=-v(t,s)^2 dt^2+ds^2,
g_B=v(t,s)^(-2) g_A=-dt^2+v(t,s)^(-2) ds^2.                (14)
```

Both have exactly the same coordinate null curves ds/dt=+/-v, and agree
at v=1. Identify the initial CFT states when v=1, and use the conformal
renormalization fixed by the stated circle vacuum, with no additional
dimensionful local metric terms. Their quantum matter actions differ by
the local Weyl-anomaly functional. With
`rho=-log v`, `R_A=-2 v_ss/v`, and trace convention
`<T^a_a>=cR/(24pi)`, integration of the trace anomaly gives

```
W[g_B]-W[g_A]
 = c/(24 pi) integral sqrt(-g_A) [rho R_A+(grad rho)^2]
 = -c/(24 pi) integral dt ds [v_s^2/v+v_t^2/v^3].          (15)
```

For the second equality integrate `2 log(v) v_ss` by parts in s.
Time perturbations vanish at the endpoints, or use periodic time for
this identity. The normalization and Lorentzian sign follow from
[Mazur--Mottola v1](https://arxiv.org/html/hep-th/0106151v1), section III,
the Wess-Zumino functional and trace anomaly (HTML equations 17, 26, 32).
Only the matter central charge is used; their optional gravity/ghost
contributions are not included. Equivalently integrate
`sqrt(-g_lambda) R_lambda=sqrt(-g_A)(R_A-2 lambda Box_A rho)`
from lambda=0 to 1. This supplies (15) directly.

For a static v, g_B is a flat cylinder of circumference J[v], so its
ground energy is -pi c/(6J[v]). Equation (15), with W=-time*E in the
static limit, reproduces exactly the positive completion (13). For
v=v_phi, the time-dependent term additionally changes the clock mass by

```
M_B-M_A = -c Q/(12 pi).                                   (16)
```

The nonlocal response and positive-frequency transition weights are
unchanged by this local action difference. Thus matching null travel
times and absorption alone fixes neither stiffness nor inertia. The
mass change can be of leading order, not merely of collar order.

For the diagnostic of section 4 with q channels, restore ell explicitly:
M_A=q ell/8 and U''=7q/(24ell). Since Q_j<=T=pi ell, the B family obeys

```
M_B >= q ell/24 >0,
K_B = U''+sum_j [contact_j+pi(Q_j-T^2/L_j)/(3L_j^2)]
    >= q/(8ell) >0.                                      (17)
```

Here each variance term is nonnegative and each contact is at least
-1/(6ell). With these coefficients there is **no upper-half-plane
linear-response pole**. To prove it, let z=w^2 and divide its kernel by z:

```
K_B(w)/z = K_B/z-M_B-2 sum W/[E(E^2-z)],
Im[K_B(w)/z]
 = -Im(z) [K_B/|z|^2+2 sum W/(E|E^2-z|^2)].               (18)
```

It cannot vanish when Im z!=0. When w=i gamma, (10) is strictly positive
with K_stat=K_B and mass M_B. These cases exhaust Im w>0. This proves
the stated absence of growing poles, not decay or nonlinear stability.
The A family has the growing pole already proved. Both mechanical
diagnostics use the same bare clock potential and identical spectral
weights. The Weyl-related matter completions explain the difference.

These are two-dimensional prescribed metric families attached to a
mechanical clock, not alternative solutions of the same four-dimensional
Einstein constraints. In particular the physical MMP perturbation is
not free to choose A or B while keeping its gravitational action fixed.
The comparison establishes the information missing from an optical
channel model; it does not select a physical completion.

Nor can a time-dependent reference subtraction be specified only for
the instantaneous ground energy while its causal response is borrowed
unchanged from another state. The [winding subtraction](warped-winding-bound.md)
uses a common local geometry and removes its local terms in both
functionals; it is a different observable. The covariant-response
framework of [Anderson, Molina-Paris and Mottola](https://arxiv.org/abs/gr-qc/0209075v1),
section III, (3.4)-(3.10), likewise treats local terms and the retarded
commutator as variations of the same renormalized effective action.
Its scalar flat-space stability result is not transferred to this model.

**Decision:** keep the physical primary and fallback provisional and
withhold any pole/stability inference from the absorption coefficient.
The immediate physical task is now explicit: derive the common-observer,
constraint-compatible metric family and its conservative collar/state
contribution, and decide whether its residual stiffness is positive on
the required scale. If it cannot be controlled, reassess the primary
before investing in a formalization or manuscript claim. The prescribed
pulse-transfer calculation remains valid in its stated scope, but is not
silently substituted for the requested physical result.

## 7. Verification and source scope

Run `python checks/collar_conservative_response.py`; the
[receipt](../evidence/collar-conservative-response.json) binds the exact
note, script and prior dependencies. Checks compare the transformed
ground energy with independent filled-fermion matrices at increasing
cutoffs; differentiate it against the stress spectrum and contact term;
compare the heat-kernel and Fourier gradient sums; verify the complex
response bound at finite resolution; and locate the growing roots of the
defined completion. These are supplementary computations, not a new
four-dimensional solution or formal certification.

The [source supplement](../evidence/collar-response-source-review.json)
records scoped primary reads and seven additional search queries. The
finite-circle transformation is established; the Fewster-Hollands sharp
energy inequality is a related predecessor, not a new wormhole no-go
theorem. Full source dependencies and the final novelty review remain
incomplete. No Lean/Lake or cloud run occurred. Historical notes and
their receipts retain their exact identities and model qualifications.

Remaining to-do list: derive and control the matched conservative
response and finite-time state, close interacting stress and physical
channel matching, decide primary/fallback novelty, obtain the consequential
physical result, independently verify and GCP-certify suitable endpoints,
and prepare the reproducible manuscript.
