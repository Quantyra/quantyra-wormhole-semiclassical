# Reflected angular modes: a throat-load bound and the remaining energy terms

2026-10-10. S051/S052/S054. A prescribed, static spherical throat gives a
quantitative bound on its imaginary-frequency scattering loss using only
two short mouth collars. Combined with the exterior current identity, this
gives an error estimate with every higher angular mode retained. Its boundary
error is an angularly weighted norm; the previous unweighted mouth-trial
estimate is insufficient to evaluate it. This is an ordinary conditional
matching result, not an accepted physical MMP return or novelty claim.

Two further requirements follow. A massive channel can reflect almost
perfectly despite exponentially small transmission. Its nonzero transmission
also has the throat front time, which need not equal the exterior front.
Eliminating these channels from a Gaussian determinant produces both a
separate determinant and a length-dependent effective ground-channel return.
Neither term can be discarded from an energy derivative by passivity alone.

## Specified throat and the known radial reduction

Prescribe a smooth static spherical segment, with positive lapse and radius,
one free massless charged Dirac field, magnetic flux integer `q>=1`, and no
electric potential. There are no photon interactions, gravitational field
equations, rotating mouths or state backreaction in this calculation. Write

```text
ds^2=N(x)^2(-dt^2+dx^2)+R(x)^2 dOmega^2,  0<=x<=L.
```

Here `x` is optical distance and `t` the common static coordinate time.
The unit-sphere magnetic Dirac eigenvalues are zero (`q` modes) and
`+/-lambda_n`, `lambda_n=sqrt(n(n+q))`, with degeneracy `q+2n` per sign.
After the spinor rescaling `N^(-1/2)/R`, each nonzero angular eigenmode
has optical mass `m_n(x)=lambda_n N(x)/R(x)`. This is the established warped
spherical reduction; a constant unitary transformation handles negative
eigenvalues. A static radial pure gauge can be removed on this segment;
its eventual closed-loop holonomy is not fixed by that removal.

Use `H=-i sigma_3 partial_x-m sigma_2`. At `H psi=i xi psi`, `xi>0`,
the radial equation and current are

```text
psi'=A psi, A=[[-xi,m],[m,xi]],
j=|psi_+|^2-|psi_-|^2,       j'=-2xi |psi|^2.                 (1)
```

Define incoming data `a=(psi_+(0),psi_-(L))`, outgoing data
`b=(psi_-(0),psi_+(L))`, and `b=K_n(xi)a`. This includes both reflection
and transmission. Existence follows from the finite-dimensional ODE:
the linear map from initial values to incoming data has trivial kernel,
since integrating (1) with `a=0` forces `psi=0`; it is therefore invertible.
For complex incoming data the same real-coefficient ODE applies. Exactly,

```text
||a||^2-||K_n a||^2=2xi integral_0^L |psi|^2 dx.              (2)
```

The direct sum over all modes is a bounded contraction on flux `ell2`.
Its summed bulk norm is at most `||a||^2/(2xi)`. Partial sums therefore
construct the corresponding `L2` weak radial solution, retaining the
angular degeneracies; no uniform angular gap in the exterior is assumed.

## Collar control of the loss, including arbitrarily high modes

Suppose there are disjoint optical collars `[0,ell]` and `[L-ell,L]`,
`2ell<=L`, with `N/R<=beta` on both. Neither `beta` nor `ell` may be
inferred from an arbitrary coordinate radius; they must be checked on the
specified throat. Let `c_n=xi+lambda_n beta`. In either collar,
`||A||=sqrt(xi^2+m_n^2)<=c_n`. Gronwall applied in either direction gives

```text
|psi(x)| >= exp(-c_n x)|psi(0)|,
|psi(L-x)| >= exp(-c_n x)|psi(L)|,     0<=x<=ell.
```

Integrate these two bounds in (2). The endpoint full norms include the
prescribed incoming component at each end. Thus for every angular mode,

```text
I-K_n^*K_n >= delta_n I,
delta_n=xi/(xi+lambda_n beta)
              [1-exp(-2(xi+lambda_n beta)ell)] >0.             (3)
```

Because each block is square, its singular values give the same bound
for `I-K_n K_n^*`. This argument is independent of all interior variations
between the two collars and of their optical separation, provided the
collars remain disjoint. Equal signs of the interior mass are not needed
for (3). Unit phases on port bases preserve the singular-value statement.

There is no mode-independent positive lower bound in (3). For constant
mass `m>0`, put `k=sqrt(xi^2+m^2)`. The exact solution gives

```text
K=[[-r,t],[t,r]],
r=m sinh(kL)/(k cosh(kL)+xi sinh(kL)),
t=k/(k cosh(kL)+xi sinh(kL)),
I-K^*K = [2xi sinh(kL)/(k cosh(kL)+xi sinh(kL))] I.             (4)
```

As `m->infinity` at fixed `xi,L>0`, reflection tends to one,
transmission is `O(exp(-mL))`, and the defect is asymptotic to `2xi/m`.
Consequently the direct-sum load can have norm exactly one while every
finite block is strictly contractive. Exponentially small transmission
does not justify replacing this load by absorption.

The transfer formula (4) is a change of basis and analytic continuation
of the known rectangular Dirac mass-barrier solution; it is a reproduced
baseline, not a new scattering effect.

## A boundary error estimate after the load is attached

This part is conditional on an actual exterior solution with the following
trace/current identity. It does not construct a matched nonspherical
four-dimensional throat/exterior operator.

Split incoming and outgoing flux spaces into spherical ground channels
`G` and all higher channels `H`. Feed prescribed data `g` into `G`, and
attach the spherical throat load only on `H`:

```text
in_G psi=g,                 in_H psi=K out_H psi.             (5)
```

All higher-mode reflections and cross-mouth transmissions are in `K`.
The ground throat propagation is left open at this step. Suppose the
exterior error `e=psi-Phi` against a smooth trial has forcing `f` and
obeys

```text
||out e||^2+2xi||e||^2=||in e||^2+2 Re <e,f>.
```

Define the actual loaded boundary mismatches

```text
b_G=g-in_G Phi,
b_H=K out_H Phi-in_H Phi,
in_H e=K out_H e+b_H.                                        (6)
```

For one finite angular block set `D_n=I-K_n^*K_n`,
`E_n=I-K_n K_n^*`. Completing the square, for any vector `y`, gives

```text
||y||^2-||K_n y+b_n||^2
 = ||D_n^(1/2)(y-D_n^(-1)K_n^*b_n)||^2
       - <b_n,E_n^(-1)b_n>.                                 (7)
```

The identity uses `I+K_n D_n^(-1)K_n^*=E_n^(-1)` and permits nonnormal
complex blocks. Its inverse-defect cost is sharp for arbitrary `y`.
Using (3) in (7), summing and completing the bulk square yields

```text
||out_G e||^2 <= ||b_G||^2
       + sum_n delta_n^(-1)||b_n||^2 + ||f||^2/(2xi).          (8)
```

The index in this sum includes both angular signs and all degeneracies.
For the infinite sum, require its right side finite and `L2` traces for
the exact solution. Apply (7) to partial block sums and take limits:
the differences of squared flux norms converge by `ell2`, and the cost
sum converges by hypothesis. This is an estimate for a solution with the
stated traces, not a proof of loaded boundary existence from (8) alone.

For an elementary sufficient boundary norm, set

```text
d_0=1-exp(-2(xi+sqrt(q+1) beta)ell).
sum_n delta_n^(-1)||b_n||^2
 <= d_0^(-1) [||b_H||^2+(beta/xi)||Lambda^(1/2)b_H||^2],       (9)
```

where `Lambda` has eigenvalues `lambda_n` on the high angular subspace.
This is half an angular derivative in the magnetic spectral norm, with
explicit flux/frequency/collar dependence. Ordinary `L2` convergence is
insufficient: a single coefficient `lambda_n^(-1/4)` tends to zero in
`L2`, while its cost in (8) grows as `sqrt(lambda_n)` at fixed collars
and frequency. That example shows insufficiency of the available norm,
not failure of the actual dipole trial.

This estimate is unweighted in spatial distance. If `P` is the nominated
ground receiving projection, `R_0` the reference return, and
`epsilon_out=||P out_G Phi-R_0||`, then for unit inputs

```text
||R_loaded-R_0|| <= epsilon_out
 + sqrt[||b_G||^2+sum_n delta_n^(-1)||b_n||^2+||f||^2/(2xi)].   (10)
```

All terms must be bounded in operator norm uniformly in those inputs.
On `0.05<=xi D<=4`, a sufficient *absolute* target for the right side is
`0.001 exp(-4)=1.83156388887e-5`, which would imply the earlier desired
`0.001 exp(-xi D)` accuracy. No calculation here establishes this target.

The spherical angular basis in (5) is essential. The exact deformed
dipole port profiles of the existing mouth trial are not the isolated
spherical ground modes. Their discrepancy and its half-derivative cost
must be included in (6)-(9), or the actual deformed-throat load must be
derived. Previous fixed-port `O(q^(-1/2))` bulk and `L2` trace limits do
not prove a joint large-flux/small-port bound in this norm.

One component of that matching cost can be bounded now. Let `g_epsilon`
be an exact normalized deformed angular zero mode at `epsilon=rho/d<1`,
`L_0` the isolated unit-sphere magnetic Dirac operator, and `P_0` its
ground projection in the same azimuthal sector. The port formulas give

```text
L_epsilon=L_0+V,    ||V|| <= M=q epsilon^2/[4(1-epsilon)^2].
```

Indeed `|delta A_phi/sin(theta)|=q epsilon^2 sin(theta)/(2D_theta)`
and `D_theta>=2(1-epsilon)^2`. Set `h=(I-P_0)g_epsilon`.
Since `L_0 h=-V g_epsilon` and the first nonzero eigenvalue has modulus
`sqrt(q+1)`, spectral calculus and Cauchy-Schwarz give the all-mode bounds

```text
||h||^2 <= M^2/(q+1),
||Lambda^(1/2)h||^2 <= ||h|| ||L_0 h|| <= M^2/sqrt(q+1),
sum_n delta_n^(-1)||h_n||^2
 <= M^2/d_0 [1/(q+1)+beta/(xi sqrt(q+1))].                    (10a)
```

Orthogonality of azimuthal sectors makes this uniform for arbitrary
normalized ground-sector superpositions. It needs no angular truncation.
It bounds this particular basis error, not the full `b_H` of (6).
When `beta<=C/rho`, `xi d>=y_0>0`, and `d_0` is bounded below, its cost
is `O(q epsilon^4+q^(3/2)epsilon^3)`. Thus
`epsilon=o(q^(-1/2))` is a sufficient joint scale for this component to
vanish. It is not asserted necessary. At fixed nonzero port fraction the
bound gives no vanishing estimate. Bulk residual constants, wrong-spin
traces, collar geometry, weak gravity and the other boundary components
must still be controlled simultaneously before this is a physical regime.

## Two requirements for the energy and travel-time argument

First, an attached throat is an additional propagation path. In (4), at
fixed finite `m,L`, `t(xi) exp(xi L)->1` as `xi->infinity`. If `L<D`,
then `t(xi) exp(xi D)` is unbounded. Thus a spherical higher-mode load
does not satisfy an exterior cross-mouth front bound `exp(-xi D)`.
This is already true for the constant-mass example; mass suppresses
low-frequency transport without changing the wavefront speed. An exterior
weight that differs between mouths need not preserve passivity of `K`.
Any claim of the exterior front for the effective return needs a new
proof with these paths included. This does not establish a ground-channel
shortcut or a self-consistent short wormhole.

Second, fix a finite angular cutoff temporarily and write exterior
scattering in blocks `S_GG,S_GH,S_HG,S_HH`. Let `Q(T)` be ground throat
propagation (including the specified spin/gauge phase), and `K(T)` the
high-mode throat load. Wherever the inverses exist, block elimination
gives the exact identity

```text
det[I-S diag(Q,K)]
 = det[I-S_HH K] det[I-R_eff Q],
R_eff=S_GG+S_GH K(I-S_HH K)^(-1)S_HG.                         (11)
```

For fixed `S`, a length variation has

```text
R_eff'=S_GH(I-K S_HH)^(-1)K'(I-S_HH K)^(-1)S_HG,
partial_T log det[I-S diag(Q,K)]
 = -tr[(I-S_HH K)^(-1)S_HH K']
   -tr[(I-R_eff Q)^(-1)(R_eff' Q+R_eff Q')].                  (12)
```

The term containing `Q'` alone is incomplete unless the other two terms
vanish or are bounded under the actual metric variation. Even if that
variation leaves a remote collar fixed and makes some derivatives small,
the required bound must be established for that family. Passivity at
each `T` gives no derivative bound. If `S` varies, its derivatives also
enter. Extending (11)-(12) to the physical determinant needs trace-class
or renormalized relative-determinant control; the finite identities alone
do not prescribe a four-dimensional vacuum subtraction.

## Verification, predecessors and decision

Run `python checks/throat_load_feedback.py --output
evidence/throat-load-feedback.json` with the pinned requirements.
The [receipt](../evidence/throat-load-feedback.json) contains four exact
algebra checks; 27 two-port ODE solutions with refinement and independent
integrated-current checks; 36 nonnormal square-completion and 36 directly
closed passive-network tests; five high-mass limits; three front-bound
counterexamples; four angular-tail examples; and 18 determinant/length-
derivative comparisons against direct determinants and finite differences.
Another 27 exact-profile cases use Jacobi spectral projections and two
angular quadrature orders to check (10a). Their partial 28-mode sums are
lower sums; the upper bound comes from the ordinary spectral argument,
not from assuming that a numerical tail vanishes.
The constant profiles also agree with the analytic barrier formula.
These are supplementary tests of the argument, not a dipole PDE solution,
a finite-flux certificate or an infinite-dimensional determinant check.

[Freivogel, Fumagalli and Tomasevic v2](https://arxiv.org/html/2606.12528v2#A3),
Appendix C, especially (271)-(294), supplies the warped spherical reduction.
The printed angular operator (286) has a misplaced gauge factor relative
to the covariant derivative in (289)-(290); use the latter and the direct
unitary radial normalization above. Optical mass includes the lapse;
the local proper mass and coordinate-frequency mass differ.
[Sundberg and Jaffe v1](https://arxiv.org/abs/hep-th/0308010v1), section 2,
especially (2.11)-(2.21), supplies the mass-barrier transfer and reflecting
limit. This section was reread from the previously retained source;
its live HTML retrieval failed, while live arXiv metadata remained available.
The loss/current argument, Gronwall estimate and Schur complement are
standard methods. No new general method is claimed.

Four targeted searches on 2026-10-10 covered Dirac mass-barrier imaginary
frequency reflection, MMP higher Landau-level feedback, dissipative Dirac
boundary defects, and Casimir Schur-complement elimination. They do not
complete a citation/edition audit or establish novelty. Their exact queries,
reading identities and impact are recorded in the planning milestone.

Decision: continue the existing physical matching question, using the
weighted boundary norm (8), and do not promote the absorbing-port energy
functional to a physical no-shortcut theorem. This supplies a concrete
criterion for reflected-mode control and identifies missing energy/front
hypotheses. It is not a separate paper target or a reason to start Lean
formalization before the physical consequence and novelty comparison close.

Remaining to-do list: establish the joint physical mouth/throat scaling
and the weighted boundary errors, construct the matched loaded operator,
control higher-mode front and determinant derivative terms plus full
vacuum/interaction corrections, complete the common-observer and novelty
comparison, then select and GCP-certify suitable endpoints and prepare
the reproducible manuscript.
