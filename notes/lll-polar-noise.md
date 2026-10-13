# Free LLL stress noise and the constrained polar metric

S072, 2026-10-13. The free lowest-Landau-level (LLL) reference supplies
an exact ordered angular stress kernel. Its traceless stress does not
directly force either propagating polar Einstein-Maxwell master on the
round product, but it does perturb the remaining metric components.
The resulting time-smeared metric variance has a convergent sum over all
LLL angular harmonics, including the light branch's angular range.

This is the free Dirac comparison, with the same fixed sphere projector
at both insertions. It is not the interacting state satisfying the local
electric and angular Gauss constraints. Induced anomaly/Hall response,
quadratic dipole sources, massive levels and physical boundaries remain
required. The small comparison variance below is not a capture bound.

## Fixed-background operator projection

Normalize the q sphere spinor zero modes by
`integral dOmega eta_i^dagger eta_j=delta_ij`, and set
`Psi=b^(-1) sum_i chi_i eta_i`. For a real orthonormal scalar harmonic
Y_lm, define the Hermitian q by q matrix

```
A_lm = integral dOmega Y_lm eta^dagger eta,
Tr(A_lm A_l'm') = tau_l delta_ll' delta_mm',
tau_l = q beta_l/(4pi),
beta_l = (q-1)! q! / [(q-1-l)! (q+l)!], 0<=l<q;
beta_l = 0, l>=q.                                      (1)
```

These are the already established finite-flux vertices, not independent
field lines. Their exact projection and metric vertex were checked in
[the finite-flux calculation](../checks/lll_projection.py); their
normalization is retained from [the angular current note](angular-current-renormalization.md).
For l>=1, Tr A_lm=0. The base four-dimensional stress coefficient is

```
delta T4_ab,lm = b^(-2) t_ab[A_lm],
D^a t_ab[A_lm]=0,       g2^ab t_ab[A_lm]=0.              (2)
```

Here t[A] is the normal-ordered free two-dimensional Dirac stress with
matrix A inserted in the flavour bilinear. These statements concern
fixed-background operator fluctuations, away from coincident contact
terms. The local trace anomaly is proportional to the flavour identity;
it has no l>=1 fluctuating component. Common c-number subtraction leaves
connected correlations unchanged. This does not remove the anomaly's
variation when the geometry changes.

The angular chirality of all zero modes is the same. Since an angular
gamma matrix reverses it, the direct projected operators J_i and T_ij
vanish. Mixed stress needs more care. In a local orthonormal frame choose
the supported chirality so that `D_1 eta=-i D_2 eta`. Then

```
eta_i^dagger <->D_A eta_j
  = -i epsilon_A^B D_B(eta_i^dagger eta_j),
T4_aA = (1/4) epsilon_A^B D_B J4_a.                    (3)
```

Reversing magnetic orientation reverses this sign. Equation (3) follows
from the symmetrized Dirac stress: its gamma_A term has a vanishing LLL
matrix element, and its gamma_a term carries i/4. The derivative includes
the spin and monopole connections before the neutral product is formed.
This mixed stress is axial on the sphere. Its scalar-gradient projection
vanishes, while its axial projection generally does not. Thus a bare
statement that all mixed stress vanishes would discard current noise.
The accepted axial metric vertex `l(l+1)/4` is the integrated version of
(3). Base currents and this mixed stress belong to the other parity
sector for the original magnetic background; they are retained there.

## Complete two-point kernel of this free reference

Take the flavour-invariant Neveu-Schwarz (antiperiodic) vacuum on a
conformal circle of circumference ell_c. The q complex Dirac fields have
two independent chiral copies. This is a specified comparison domain,
not reflecting strip or physical mouth matching. For either chirality,
filled negative half-integer fermion modes give

```
J_n[A] = sum_(r in Z+1/2) :c_r^dagger A c_(r+n):,
L_n[A] = sum_r (r+n/2) :c_r^dagger A c_(r+n):,
<J_n[A] J_-n[B]> = n Tr(AB),
<L_n[A] L_-n[B]> = n(n^2-1) Tr(AB)/12,
<L_n[A] J_-n[B]> = 0,                       n>0.       (4)
```

The first two formulas are simply sums over the n particle-hole pairs
crossing the filled-sea edge. Their weights are respectively 1 and
`(r+n/2)^2`; the mixed sum of `r+n/2` vanishes. They hold for arbitrary
Hermitian A,B, including noncommuting matrices. They do not assert a
closed independent Virasoro algebra for each angular vertex or a
Gaussian distribution of the composite stress. Higher moments require
the fermion Wick expansion and do not factor into stress pairs.

Let d=2pi/ell_c and let u be the relevant chiral coordinate. In physical
Dirac density normalization the fluctuating fields are

```
j_A(u) = ell_c^(-1) sum_n J_n[A] exp(-i d n u),
t_A(u) = (2pi/ell_c^2) sum_n L_n[A] exp(-i d n u).
```

The l>=1 means vanish. With `z=exp[-i d(u-v-i epsilon)]`, epsilon>0,

```
<j_A(u) j_B(v)> = Tr(AB) z/[ell_c^2 (1-z)^2],
<t_A(u) t_B(v)> = (2pi/ell_c^2)^2 Tr(AB)
                            z^2/[2(1-z)^4],
<t_A(u) j_B(v)> = 0.                                  (5)
```

The distributional boundary value epsilon->0 is understood. Equivalently
the ordered stress spectrum is the positive measure in (4), at positive
frequencies d n for n>=2. The reverse ordering is its negative-frequency
partner. Their difference supplies the nonlocal retarded commutator;
local variations of the stress definition are additional response
contacts. On the conformal product the fluctuating t_++ and t_-- remain
these fields; the common Weyl-anomaly term is a c-number. Equation (1),
the tensor lift (2), and the axial derivative (3) therefore fix every
connected two-point function of the direct free LLL stress and current.

For a normalized real time window w at a fixed spatial point, put
`W_n=integral dt w(t) exp(-i d n t)`. Its chiral stress variance is

```
Var(t_A[w]) = (2pi/ell_c^2)^2 Tr(A^2) S_w/12,
S_w = sum_(n>=2) n(n^2-1) |W_n|^2.                    (6)
```

Finite S_w is a precise regularity condition. Point evaluation has
W_n=1 and divergent variance. A Gaussian window of standard deviation
sigma has `|W_n|^2=exp[-d^2 sigma^2 n^2]`. Smearing is part of the
observable, not a regulator whose removal is claimed to be harmless.
The current analogue needs `sum n |W_n|^2<infinity`.

## Vanishing direct master source, nonzero metric response

Use the variables of [the sourced product response](quadrupole-product-response.md).
For the free prescribed polar stress in (2), S_T=0 and S^a_a=0. Its
angular electric current is zero. After duality this means Bcal=0 in
this polar sector; the original base electric current belongs to the
other sector and is not set to zero. There is no polar Jtilde source.
Consequently

```
r_F = 0, r_A = 0,
chi=0, Acal=0,
h_ab,lm = [2 kappa^2 / l(l+1)] t_ab[A_lm], l>=2.       (7)
```

The zero master fields select no added homogeneous wave. The remaining
metric follows from the sourced constraints. It is generally nonzero.
In static coordinates `ds2^2=-B dT^2+d xi^2/B`, introduce the optical
coordinate x by dx/dxi=1/B. Every free traceless conserved stress is
locally a sum of two null profiles,

```
t_TT=F(T-x)+G(T+x),
t_Txi=[-F(T-x)+G(T+x)]/B,
t_xixi=[F(T-x)+G(T+x)]/B^2.                            (8)
```

Substitution into all linearized four-dimensional Einstein-Maxwell
components verifies (7), with source T4_ab=t_ab/b^2 and no polar Maxwell
perturbation. In particular one must include the metric variation of
the background Maxwell stress even when delta F=0. The independent
checker does this for l=2,3, arbitrary frequency including zero and
symbolic C=b^2/a^2. The general l result follows from the sourced
constraint reconstruction and the harmonic eigenvalue L=l(l+1).

This is the linear response to a prescribed free source on the classical
product. Once the induced stress, shifted reference background or
nonlinear residuals are included, (7) is not the full semiclassical
equation and zero master forcing cannot be assumed. In particular the
quadratic dipole electric stress from the all-axis calculation is not
traceless in this two-dimensional sense.

## Summing the constrained metric noise

Use `ds2^2=a^2 sec^2(x)(-d tau^2+dx^2)`. At the center x=0 define the
fractional metric component `H=delta g_tau tau/a^2`, in the same scalar
gauge and background coordinates as (7). This is a specified linear
metric observable, not a diffeomorphism-invariant capture probability.
Let H_w be its time smear, retaining only l>=2. Orthogonality, the two
independent chiral copies and (1),(6),(7) give

```
V_w = (1/(4pi)) integral dOmega <H_w(Omega)^2>
    = kappa^4 q S_w A_q / (6 a^4 ell_c^4),
A_q = sum_(l=2)^(q-1) (2l+1) beta_l/[l^2(l+1)^2]
    <= 1/4,                                           (9)
```

where kappa^4 means `(kappa^2)^2`. The inequality is a complete angular
sum: `0<=beta_l<=1` and `(2l+1)/[l^2(l+1)^2]=1/l^2-1/(l+1)^2`.
It does not count harmonics as independent excited particles. More
precisely, for q>=2,

```
0 <= 1/4-A_q <= [H_(q-1)+H_q-5/2]/q+1/q^2.             (10)
```

Use the established bound `1-beta_l<=l(l+1)/q` in the finite sum and
retain the omitted tail. Thus there is no hidden angular ultraviolet
divergence in this smeared *constrained polar component*. Other parity,
derivative observables and nonlinear sources have not been bounded by
(9). In particular the time regularity in (6) is still necessary.

For the candidate `q=k^4`, `L4/l_p=k^5`, `C=sqrt(1+12pi k^6)`, the
product relation gives `a^2/l_p^2=k^10(C-1)/(6C)`. At ell_c=pi,

```
sqrt(V_w) <= sqrt(96 S_w)/pi * C/(C-1) * k^(-8).        (11)
```

The receipt evaluates a Gaussian window sigma=1/4, as well as two
nearby widths, by a positive spectral sum with an explicit tail bound.
At k=2000, sigma=1/4 gives `sqrt(V_w)<=3.007696e-26`.
The Gaussian has support on all times. It is a comparison measurement,
not a causal finite switching pulse or an all-future stability test.

Finite-q resolution is also retained in the soft angular sector of the
previous response. For every l through its l_max,
`1-beta_l<=l_max(l_max+1)/q`. At k=2000 this is below 0.007412; for the
quadrupole alone the exact fractional loss is
`6q/[(q+1)(q+2)]<3.75e-13`. These are relative changes in the free
two-point spectral weight, not relative errors in the interacting source
or capture energy. They cannot justify discarding other corrections of
similar angular derivative order.

## Induced response, provenance and next physical calculation

[Sadhukhan, sections 5.2, 5.4-5.5 and 6](https://arxiv.org/html/2609.33511v1)
treats an angle-parametric LLL expectation response, with the stated
restriction l much smaller than sqrt(q). Its perturbed lift includes
anomaly pressure and an induced Hall current required by four-dimensional
force balance. Section 5.4 also describes free chiral radiation that
changes constrained metric components while leaving its polar masters
unchanged at that order. These are precedents, not a novelty claim here.
The paper does not supply the finite-q ordered noise kernel or the
constrained driven state used by the present program. Its proposed
effective-action representative and Euler coefficient carry qualifications
in section 6.4; no first-principles determinant matching is imported.

The exact mixed operator (3) and the vanishing vacuum mean of that
operator must not be confused with induced Hall response. Changing the
metric or gauge field varies the projector and the operator definition
and can involve virtual massive levels and local contacts. A fixed
projector Wick calculation alone does not compute those terms. Likewise,
zero free stress-current two-point function in (5) does not make all
mixed correlations vanish in an interacting or driven state.

The [checker](../checks/lll_polar_noise.py) verifies the direct tensor
response, actual fermion particle-hole spectral weights with noncommuting
vertices, independent zero-mode matrix variances, all-angular sum bounds,
and high-precision smeared candidate estimates. Its
[receipt](../evidence/lll-polar-noise.json) records the source and dependency
hashes. This is a free reference calculation with an ordinary derivation,
not Lean certification or verification of an evolving wormhole.

The next step remains the ordered source kernel in the joint Gauss-law
state, including the quadratic dipole/metric residuals, induced currents
and massive/contact terms. The free kernel above provides a normalized
baseline. It does not justify assigning this vacuum to the dressed,
boundary-driven physical system or dropping the light propagating branch.

Remaining to-do list: interacting constrained source and induced response;
physical state/operator and boundary matching; evolving renormalized
stress and geometry; preparation, sufficient stability and finite-energy
transport, with creation distinguished from opening an existing bridge.
