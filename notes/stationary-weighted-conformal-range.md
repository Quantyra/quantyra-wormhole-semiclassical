# Weighted conformal range and the remaining coupled equations

2026-10-17, S139/S074. Complete ordinary argument before verification.

The [coupled unweighted obstruction](stationary-coupled-asymptotic-obstruction.md)
does not exclude a weighted boundary problem. This note constructs the
weighted range of the full trace restricted to conformal variations on
an actual fast comparison geometry. It also gives an exact conformal/shape
coordinate system, retains the radial constraint missing from trace plus
anisotropy, and states the resulting coupled linear reduction.

This is progress toward the coupled solver, not that solver's existence
theorem. In particular no nonlocal shape inverse, actual shooting
Jacobian, uniform inverse constant or physical branch is established.

## Geometry and the exact conformal operator

Keep the retained action, magnetic flux, positive prepared temperature,
clock, paired Dirac domain and subtraction of the preceding notes.
Let g be a regular static spherical comparison with a nonextremal bolt,
smooth on the finite exterior and at the Euclidean bolt, with

```
F=lambda r^2+1-2m/r+O(r^-2 log^p r),
alpha=O(r^-4 log^p r), lambda=lambda_e>0,
```

and five differentiated weighted bounds. A genuine solution is NOT
assumed in this note. Use the positive Euclidean section for the radial
operator and norms. Its static scalar Laplacian is the same radial
expression as the Lorentzian stationary one. Write

```
epsilon=G/(60pi), a=1+11epsilon lambda/3,
lambda0=lambda+11epsilon lambda^2/6,
R=R4, J_ab=Ric_ab-(R/4)g_ab,
M=F_em^2, Z=R+12lambda, w=exp(alpha)r^2, p_r=wF,
E=E_trace=(a+epsilon Box)Z-S[g]=-g^ab E_ab.
```

The [exact retained trace](stationary-trace-range-condition.md) is
equivalently

```
E=R+12lambda0+epsilon[Box R-11R^2/72
                          +7C4^2/12+11J^2/6-20M].       (1)
```

For a conformal variation delta g=2sigma g, hold the magnetic two-form
fixed in the current coordinates. This is not a fixed-areal-radius
variation. In four dimensions the exact first variations are

```
delta R=-6Box sigma-2R sigma,
delta(Box R)=Box(delta R)-2sigma Box R+2 grad sigma.grad R,
delta C4^2=-4sigma C4^2,
delta J^2=-4sigma J^2-4J^ab nabla_a nabla_b sigma,
delta M=-4sigma M.                                    (2)
```

The last equality differs from delta M=0 in fixed-q AREAL coordinates;
the radial coordinate restoration accounts for the difference. No
independent quantum state variation is used. The complete bulk trace
is local even though the other tensor equations retain their nonlocal
Dirac response.

Substituting (2) into (1) and using
`nabla_a J^ab=(1/4)nabla^b R` gives the exact operator at the comparison:

```
L_g sigma := D_conf E[sigma]
 =-6epsilon Box^2 sigma+div(A grad sigma)+V sigma,
A^ab=(-6-epsilon R/6)g^ab-(22epsilon/3)J^ab,
div A=-2epsilon grad R,
V=-2R+11epsilon R^2/18-4epsilon Box R
                      +80epsilon M-7epsilon C4^2/3-22epsilon J^2/3
 =2R+48lambda0-4E.                                    (3)
```

Thus L is formally symmetric in the measure w dr, also at a comparison
with nonzero equation defect. The conformal variation of the trace
density is `(L+4E)sigma`; this adds a real multiplication term and does
not change the formal symmetry of L. We use the scalar derivative L,
not the density derivative, in the range theorem.

In radial form, set `b=-6-epsilon R/6-(22epsilon/3)J^r_r` and
`D=Box=w^-1 partial_r(p_r partial_r)`. Then

```
L=-6epsilon D^2+w^-1 partial_r(p_r b partial_r)+V.       (4)
```

All coefficients in (4) are determined by the actual comparison
geometry. Four metric derivatives enter V; a fifth gives the first
conormal coefficient bounds. No sixth asymptotic derivative or negative
scalar-curvature assumption is needed for the theorem below.

## A noncritical weighted range theorem

Fix `delta=3` and `0<eta<1`. On the end use x=log(r/B) and the norm
obtained from the uniformly local C^(j,eta) norm of exp(delta x)u on
unit x intervals. On a fixed compact region, including the bolt, use
the usual C^(j,eta) norm on the smooth Euclidean geometry restricted
to static spherical functions. Define X_delta by j=4 and Y_delta by
j=0. At the bolt these are rotation-invariant smooth-disk regularity
conditions, not arbitrary boundedness in the singular areal coordinate.
The theorem also holds for every fixed 3/2<delta<4.

The operator L:X_delta -> Y_delta is bounded Fredholm and

```
index L=-1,
k=dim ker L is either 0 or 1,
dim coker L=k+1.                                     (5)
```

Its range is exactly the common kernel of k+1 moment functionals
defined below. These statements concern the conformal TRACE block;
they assign no index to the complete nonlocal tensor system.

Here is a direct proof, including the endpoint data.

At the bolt set v=D sigma. A regular solution of L sigma=f obeys

```
p_r sigma'(r)=integral_rh^r w v ds,
p_r v'(r)=(1/(6epsilon)) integral_rh^r
                    w[b v+F b' sigma'+V sigma-f]ds.   (6)
```

Because p_r has a simple zero and F sigma' is regular, these are
regular Volterra equations with two initial values sigma_h and v_h.
They give a two-dimensional regular homogeneous space H and an affine
regular space for each forcing. A contraction on a sufficiently short
bolt interval followed by ordinary continuation proves existence and
uniqueness for those data. Smooth forcing gives smooth disk solutions;
the regular radial elliptic estimate gives the stated Holder mapping.
The local interval and constants may depend on the fixed positive
temperature. No uniform extremal-limit statement is made.

At infinity the coefficients converge, with exponentially decaying
conormal differences in x, to the constant-coefficient polynomial

```
I(s)=-6 lambda [a+epsilon lambda s(s+3)] [s(s+3)-4].     (7)
```

Its four distinct roots are

```
1, -4, -3/2+i omega, -3/2-i omega,
omega^2=a/(epsilon lambda)-9/4>0.                     (8)
```

For completeness, rewrite the fourth-order equation as a first-order
four-component system in x and conjugate by exp(delta x). The limiting
matrix has one eigenvalue with negative real part and three with
positive real part. Its stable/unstable integral kernels decay
exponentially. On a sufficiently remote end the coefficient difference
has small sup norm, so their integral operator is inverted by a
contraction. This constructs a bounded particular solution for every
weighted forcing, plus a one-dimensional weighted homogeneous space S.
It also supplies the weighted derivative estimates. Extension to a
fixed interior section is ordinary ODE continuation. This argument
does not use an assumed Fredholm theorem for the nonlocal response.

At that section identify solutions by their four Cauchy data. Matching
the bolt particular solution plus H to the outer particular solution
plus S is a finite-dimensional map

```
H direct-sum S -> R^4, (h,s) -> h-s.                  (9)
```

Its index is 2+1-4=-1. Its kernel is H intersect S and has dimension
at most one. Its cokernel has dimension k+1. The particular solution
maps are bounded, so its finite-dimensional range condition is closed
and supplies a bounded solution operator on a complement of the kernel.
Every matching obstruction is reached by a smooth compact forcing:
equivalently, its nonzero adjoint solution pairs nontrivially with a
bump where that solution is nonzero. This identifies the cokernel of
L, rather than merely giving an upper bound for it. Thus (5) follows.

The ordinary edge-operator theory in
[Mazzeo--Vertman, section 2](https://arxiv.org/html/1307.2266)
explains the role of nonindicial weights and asymptotic Cauchy data.
Its general theorems have additional normal-operator and regularity
hypotheses. Here the radial Volterra/matching proof supplies the result
directly, including the retained logarithmic coefficient falloffs;
that reference is not a substitute for a full tensor boundary theorem.

## Actual adjoint tests and exact solvability

The Lagrange concomitant for two radial scalars is

```
W(chi,sigma)=-6epsilon p_r[chi(D sigma)'-chi'D sigma
                                 +(D chi)sigma'-(D chi)'sigma]
                      +p_r b(chi sigma'-chi' sigma),
W'=w(chi L sigma-sigma L chi).                        (10)
```

The regular bolt space H is Lagrangian for this nondegenerate
four-dimensional skew pairing: its dimension is two and W vanishes
there. The outer fast line S has a three-dimensional annihilator.
This annihilator is exactly the homogeneous solutions with no growing
r^1 coefficient. To see the nonzero pairing, the r^1/r^-4 pair has
second-order Wronskian -5lambda in that order; its multiplier in (10)
is -6(a+4epsilon lambda), which is nonzero. The two oscillatory branches
and the fast branch pair to zero with the fast line. Decaying coefficient
corrections preserve these limiting pairings.

Let N be the space of solutions L chi=0 regular at the bolt and with
no growing r^1 coefficient. Finite-dimensional symplectic linear algebra
gives `dim N=k+1`. Its elements have oscillatory r^-3/2 asymptotics
and/or the r^-4 branch. Choose any real basis chi_1,...,chi_(k+1).
Then the exact range condition and solution statement are

```
L sigma=f has sigma in X_delta
    iff ell_j(f):=integral_rh^infinity w chi_j f dr=0 for every j.
The solution is unique modulo ker L.                 (11)
```

These moments converge continuously on Y_delta since delta>3/2.
The Green identity proves necessity. Sufficiency follows by taking
the annihilator of (9), which is precisely N; its pairing with the
particular-solution mismatch is (11). The functionals are independent
because compactly supported forcings distinguish linearly independent
continuous chi_j. This also proves the cokernel assertion used above.

These chi_j solve the full fourth-order conformal derivative on the
ACTUAL comparison. They are not generally the previous second-order
test phi_g satisfying (a+epsilon Box)phi_g=0. Replacing chi by phi,
or assuming k=0 without controlling the matching matrix, would lose
the new range information. A finite inverse on a chosen complement
does not supply its numerical size or uniformity over preparation.

The auxiliary X_delta space contains tails outside the fixed fast
physical phase. Its application therefore needs the following bootstrap;
the weighted Fredholm statement alone is not phase preservation.
The auxiliary weight is deliberately not the critical root delta=4.
If a forcing has r^-4 log^p falloff, the simple -4 root generally
produces r^-4 log^(p+1) in sigma. For delta=3 the coefficient difference
acting on an X_delta solution is O(r^-5 log^q r); variation of constants
then proves this improved asymptotic order and excludes the other roots.
The coefficient of the first resonant logarithm for forcing c r^-4 is
`c/I'(-4)`, with

```
I'(-4)=30lambda(a+4epsilon lambda)>0.                 (12)
```

One additional conormal derivative follows from the first-order system
when the forcing has one and the five metric jets are controlled.
Smooth finite-region data retain smooth bolt regularity. This recovers
fast asymptotics with the stated logarithmic gain. It does NOT prove
closure of a Newton iteration in a phase with a preassigned fixed log
power: the coupled asymptotic coefficients must still be matched.

## An exact conformal/shape coordinate system at fixed temperature

The split need not introduce a guessed independent scalar horizon datum.
For any metric above define

```
1/u(r)=integral_r^infinity exp(alpha(s))/s^2 ds,
sigma(u)=log(r/u),
f(u)=exp(2alpha(r))F(r)u^2/r^2.                       (13)
```

The integral is positive, du/dr=exp(alpha)u^2/r^2>0, and u/r tends to
one. A direct substitution gives the exact global representation

```
g=exp(2sigma(u))[-f(u)dU^2+du^2/f(u)+u^2 dOmega^2],
F=f(1+u sigma_u)^2,
alpha=sigma-log(1+u sigma_u).                         (14)
```

The Euclidean sign is used at the bolt. The fast lapse implies
sigma=O(u^-4 log^p u), while f has the same lambda, constant term 1
and mass coefficient m as F. Thus this split preserves the fixed
boundary clock and fast data rather than changing them by definition.
The integral normalization removes its residual radial ambiguity.

Let u_h be the conformal representative's bolt. Its regularity data are

```
f(u_h)=0, f'(u_h)=2kappa,
r_h=exp(sigma_h)u_h, beta=2pi/kappa.                  (15)
```

The second equality follows directly from (13): f'_h=exp(alpha_h)F'_h.
A smooth conformal variation preserves the period; sigma_h need not
vanish. The two regular conformal ODE amplitudes in (6) are therefore
legitimate tangent data at fixed temperature, but they are not free
parameters of solutions to all field equations. Shape data include
u_h and f, subject to (15), the outer data and the remaining equations.
At fixed u, varying (15) gives
delta f_h=-2kappa delta u_h and
delta f'_h=-f''_h delta u_h. For actual differentiation use a smooth
pullback to a common bolt manifold; singular areal-coordinate components
at a moving horizon are not legitimate extra degrees of freedom.

## The radial constraint and the coupled linear reduction

Write E_t,E_r,E_theta for the mixed diagonal components of the FULL
retained tensor equation. Covariance with fixed monopole gives

```
E_r'+(alpha'+F'/(2F))(E_r-E_t)+(2/r)(E_r-E_theta)=0.    (16)
```

The vanishing full trace and vanishing anisotropy E_r-E_t do not alone
give E_ab=0. Instead (16) integrates exactly to

```
E_r=E_t=C/r^4, E_theta=-C/r^4.                       (17)
```

This residual is compatible with the allowed falloffs. A fixed magnetic
flux does not set C=0 by itself; the field equations must do so. Retain
one radial constraint, for example K=r0^4 E_r(r0)=0 at a fixed interior
areal section. Once trace and anisotropy vanish it equals C. Equivalently,
for regular data one can impose r_h^4 E_r(r_h)=0. Thus the exact reduced
equations are trace, anisotropy and this one constraint. These are two
functional equations plus one integration constraint, not three
independent tensor functions. At a smooth bolt the anisotropy is O(F),
so its natural regular representative is (E_r-E_t)/F.

Let T denote the trace residual and R_low the anisotropy/constraint
residual in these coordinates. On a common smooth manifold the
linearization has blocks

```
[ L  B ] [s] = -[T    ],
[ C  D ] [v]   [R_low].                              (18)
```

Here s is the conformal variation and v denotes all allowed shape
variations, including the bolt-position datum and its constraints.
B,C,D are derivatives of the actual full residual, not the AdS limiting
operator. At a comparison with nonzero defect they include its pullback,
endpoint, state/domain and coefficient variations.

Let Pi collect the k+1 independent moments (11), choose a fixed
complement of ker L, and let G be its bounded range inverse, extended
by a fixed projection off the cokernel. The trace block is equivalent to

```
Pi(Bv+T)=0,
s=-G(T+Bv)+z, z in ker L.
```

Consequently the full linear step is equivalent to

```
Pi(Bv+T)=0,
(D-CGB)v+Cz=-R_low+CGT, z in ker L.                  (19)
```

This is an algebraic reduction on smooth tangents; using it as a bounded
coupled solve requires appropriate spaces and bounds for the NONLOCAL
shape blocks. No such bounds or invertibility are inferred from (5).
If k=0 there is one projected shape equation. If k=1 there are two and
one conformal kernel amplitude. Neither Pi B's rank nor the invertibility
of the lower system follows from counting. In particular a horizon
shift is not automatically transverse to the compatibility condition.

Equation (19) uses the actual candidate derivative. It does not drop
the defect-weight or moving-endpoint terms in the previous nonlinear
compatibility derivative. The new moment basis is different and has
its own dependence on the comparison geometry. No residual size is
converted here into an error on an exact exterior or its charge.

There is also a concrete reason not to assign the full system index -1
by assuming its cross-coupling is compact. Already on the AdS limiting
geometry, the accepted full conformal tensor factorization gives

```
(DE_r-DE_t)/F=-2 partial_r^2[(a+epsilon Box)sigma].
```

Its natural r^2-rescaled outer symbol on sigma~r^s is
`-2s(s-1)[a+epsilon lambda s(s+3)]`. At the admissible trace root
s=-4 this is `-40(a+4epsilon lambda)`, which is nonzero. Thus an
allowed homogeneous trace tail generally has a leading shape residual.
The actual lower blocks and their Schur normal operator must be
analyzed together. Neither full index -1 nor full non-surjectivity
follows from the conformal trace theorem.

The practical gain is an OFF-SHELL range formulation for the accepted
candidate: its actual L and chi_j legitimately project the Newton trace
equation. This avoids substituting the genuine-solution compatibility
formula at a candidate, without claiming that the new projected equation
is the old nonlinear moment's derivative. The remaining shape response
and its projection are explicit obligations in (19).

## Verification contract and next decision

The complete ordinary argument above precedes a focused exact checker.
It may verify conformal-variation coefficients, the divergence form,
the density/defect term, radial Green identity, indicial polynomial and
nonzero pairings, resonance coefficient, coordinate transformation,
temperature, Ward residual and block elimination. Every symbolic
residual must vanish exactly; there is no shooting or numerical tolerance.
Endpoint existence, weighted estimates, matching/Fredholm conclusions,
asymptotic bootstrap and full functional-domain claims remain ordinary
arguments. Freeze exact input bytes, retain every run, and preserve
failed input snapshots. No predecessor checker or Lean run is selected.

The next coupled construction must control the actual matching matrix,
physical logarithmic asymptotics and nonlocal shape response in (19),
including the radial constraint. This replaces a presumed single scalar
shooting knob with an explicit conditional range problem. It neither
selects a physical boundary branch nor closes S139.1-S139.2.

Remaining to-do list: verify this weighted conformal range and coupled
reduction; control the actual shape/constraint map and physical branch
or justify controlled spatial reduction; establish exterior/state and
charge matching, transfer S138, and complete S074-S080.
