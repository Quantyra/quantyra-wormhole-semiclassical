# Electric screening with the Laplacian retained

2026-10-10. S051/S052. This is a controlled baseline for the **prescribed
local LLL screening equation**, not a new screening mechanism or a physical
wormhole solution. It proves existence, a global potential bound and a
local large-flux limit for the point-dipole equation. It also separates the
bounded electric multiplication term from the still uncontrolled rotating
Dirac operator. The [joint matching argument](joint-mouth-matching.md)
remains conditional on its loaded operator and physical interface.

## Equation and state boundary

Use separation `d=1` initially. Put magnetic charges at `z=+/-1/2` and set

```text
bvec(x) = 1/2 [(x-ez/2)/|x-ez/2|^3-(x+ez/2)/|x+ez/2|^3],
b=|bvec|, B=q bvec, n=ey.bvec/b, h^2=2 pi^2/(g^2 q),
a=-g^2 A_t/(4 pi^2 Omega).
```

Here `q>0`, `g>0`, and rotation is about `ey`. The linear-in-rotation
equation under examination is

```text
(-Delta+g^2 |B|/(2 pi^2)) A_t = -2 Omega ey.B,
             (-h^2 Delta+b) a = b n.                         (1)
```

This is [MMP v3](https://arxiv.org/html/1807.04726v3), Appendix G,
(G.35)-(G.39), including its specified bosonized static current. MMP
already predicts the screened potential and suppressed drift. Appendix A
also includes dynamical photons; the support state is not a completely
free collection of charged fermions. The result below justifies the
elliptic part of this baseline. It does not derive the current at finite
flux, fix charge/holonomy sectors on the closed wormhole, or control the
higher orders in rotation. A constant potential and the choice of electric
state cannot be changed independently when interpreting (1).

There is also a global limitation: the linear rotating-coordinate expansion
is not uniform where `|Omega| |x|` reaches one. The corotating Killing field
becomes spacelike outside the light cylinder, as MMP section 6.1.3 discusses.
Solving (1) on all of space is a well-defined elliptic baseline, not a
construction of a global corotating ground state. A physical application
must instead match a finite timelike region to an appropriate exterior
radiative state and control the omitted rotation terms.

## Global weak solution and exact amplitude bound

Let `D^{1,2}(R^3)` be the completion of compactly supported smooth real
functions in the gradient norm, represented in `L6`. For every `h>0`,
(1) has a unique solution `a` in this space and

```text
|a| <= 1 almost everywhere,
|A_t| <= 4 pi^2 |Omega|/g^2.                                (2)
```

**Proof.** The dipole field has no zero at a finite nonsingular point.
Off the axis, a zero of its cylindrical radial component requires `z=0`,
where its axial component is nonzero. On the axis the axial component is
nonzero away from the two charges. Near each charge `b=O(r_i^-2)`; at
infinity `b=O(r^-3)`, with a direction-uniform positive leading dipole
magnitude. In particular `b<=C(r_+^-2+r_-^-2)`. The translated Hardy
inequality makes `integral b u v` a bounded bilinear form on `D^{1,2}`.
The source `bn=ey.bvec` belongs to `L^{6/5}`: its local inverse-square
singularities and its inverse-cube far tail are both integrable to that
power. Sobolev and Cauchy-Schwarz therefore bound the source functional.
The form `h^2 integral grad u.grad v+integral b u v` is coercive in the
gradient norm. Its elementary Hilbert-space variational construction
gives existence and uniqueness.

Test the weak equation with `(a-1)_+`, which belongs to `D^{1,2}` by
Lipschitz truncation and is supported where `a>1`. On that set `n<=1`,
so both `h^2 integral |grad(a-1)_+|^2` and
`integral b(a-n)(a-1)_+` are nonnegative and their sum is zero.
Thus `a<=1`. Apply the same argument to `-a`. Smoothness away from the
charges follows by local Poisson regularity and bootstrapping. The `L6`
condition excludes an added constant electric potential. On balls of
radius proportional to `|x|` at large `|x|`, the scaled potential and
source are `O(|x|^-1)`; the interior Poisson estimate gives
`|a(x)|<=C |x|^-1/2 ||a||_{L6(ball)}+C/|x|`, with constants allowed to
depend on fixed `h`. Thus `a` tends to zero at infinity. No boundary
condition at a fictitious material matching shell was used. QED.

The same truncation proves `|a|<=max(1,M)` in a finite-mouth exterior
with Dirichlet data bounded by `M`. This is a *conditional* extension:
the actual throat electrostatic data have not been derived. The full-space
point-charge solution is not asserted to be those data.

## Local approximation, including shrinking collars

For a cube `Q=x0+[-rho,rho]^3` away from the charges, suppose `b>=b0>0`
there and `M2=sup_Q |Delta n|`. Put `w=a-n`. Then

```text
(-h^2 Delta+b)w=h^2 Delta n,
|w| <= h^2 M2/b0 + 12 exp[-sqrt(b0) rho/(2h)]
                  on x0+[-rho/2,rho/2]^3.                  (3)
```

To verify (3), the nonnegative barrier

```text
v=sum_{j=1}^3 cosh(sqrt(b0)(x_j-x0_j)/h)
                       /cosh(sqrt(b0)rho/h)
```

is at least one on the boundary and satisfies
`(-h^2 Delta+b)v=(b-b0)v>=0`. Since `|w|<=2`, compare `+/-w` with
`h^2 M2/b0+2v`. In the inner cube each cosh ratio is at most
`2 exp[-sqrt(b0)rho/(2h)]`. This proves the displayed bound, with no
assumption that screening remains strong at infinity.

Fix a compact set of trimmed field-line labels, as in the joint trial.
Let `r_*` be distance to the nearer charge, capped by one. At every
point of its finite tube one can choose `rho=c r_*` so that

```text
b0>=c_b/r_*^2, M2<=C/r_*^2, |grad n|<=C/r_*.
```

These follow from the nonsingular dipole field on the compact rescaled
near-mouth regions and the compact middle region. Thus (3) gives

```text
|a-n| <= C h^2 + 12 exp(-c/h),
|grad a| <= C/r_*                                             (4)
```

uniformly down to any positive shrinking radius, in particular `p=q^-3/4`.
For the gradient claim use a further interior ball of radius `h r_*`.
The Poisson gradient estimate for `w` bounds its gradient by
`C (h r_*)^-1 ||w||_infinity + C h r_* ||Delta w||_infinity`.
Here `Delta w=(b/h^2)w-Delta n`, and (3) on a slightly larger cube
gives `|grad w|<=C h/r_*` for sufficiently small `h`; the exponential
term is absorbed. Add the displayed bound on `grad n`.

An exact check at a single monopole is especially useful: if
`b=c/r^2`, `n=y/r`, then `Delta n=-2n/r^2` and
`a=c n/(c+2h^2)` solves (1). Its error is uniformly `O(h^2)` even as
`r` shrinks. A fixed-distance compactness statement alone would not have
supplied this uniformity.

On fixed nested compact sets the next approximation is
`a2=n+h^2 Delta n/b`. Its equation has residual
`-h^4 Delta(Delta n/b)`, so the same comparison gives a local `O(h^4)`
error, plus an exponentially small boundary term with a bounded prefactor.
This refinement is tested below but is not needed for (2) or (4).

## What this controls in the return calculation

Restore `d` and keep the fixed positive frequency window `xi d=Xi`.
In a specified stationary flat Dirac model an electric potential enters
as a real scalar multiplication `V`; its sign is immaterial here.
Consequently it does not change the current identity used for the
[throat-load estimate](throat-load-feedback.md). For the same trial `Phi`
it adds a residual of size at most `||V||_infinity ||Phi||`.
The existing fixed-window trial obeys `||Phi||<=C sqrt(d)`.
The current estimate weights a residual by `1/sqrt(2xi)`. Therefore,
**if the same loaded problem exists and the electric boundary data obey
the bound above**, the electric term alone adds at most

```text
C d ||V||_infinity <= C |Omega| d/g^2                      (5)
```

to the return error. No norm of a derivative of a large-flux wavefunction
was used. For fixed small `g` and `d/r_e=q^nu`, MMP's orbital scaling
gives `Omega d=2q^(-nu/2)`. With `3/4<nu<1`, (5) is
`o(q^-3/8)` and thus compatible with the previous joint error rate.
This is a scaling consequence in the prescribed current model, not a
finite parameter certificate for a fixed interacting theory. If `g`
varies, both `g^2 q -> infinity` and `Omega d/g^2 ->0` must be checked.

The same local estimates yield a bound on the *electric contribution*
to the usual guiding-center drift: `|grad A_t|/|B|<=C |Omega| r_*/(g^2 q)`
and integrated angular displacement at most `C |Omega| d/(g^2 q)` on
a trimmed line. This statement uses the reduced guiding-center equation;
it does not prove a quantum return holonomy or sum the inertial terms.

In particular the shift/Coriolis part of the rotating Dirac operator is
a differential operator. A bare bound on its coefficient does not imply
the multiplication estimate (5); acting on the unmodified static Gaussian
can introduce growing transverse derivatives. The real rotating metric,
spin frame, induced magnetic/current response and throat boundary data
still have to be treated together. Neither (2) nor (5) bounds the load's
length derivative, its separate determinant, or the total vacuum energy.

## Verification and research decision

Run `python checks/rotating_electric_screening.py --output
evidence/rotating-electric-screening.json` with the pinned requirements.
The check independently solves the azimuthal `m=1` elliptic equation on
finite cylindrical boxes using a sparse finite-difference operator,
compares three meshes and two outer boundaries, checks the discrete
maximum bound and linear residual, and compares both local asymptotic
orders. Exact algebra checks the monopole solution, cosh barrier, source
normalization and joint exponents. Finite-box numerics illustrate the
ordinary proof; they do not certify the infinite-domain error or the
physical local-current closure. The source state, finite mouths and
rotating Dirac problem are not silently replaced by the numerical box.

The cached [2023 accepted manuscript](https://www.osti.gov/servlets/purl/2419919)
was also compared directly for Appendix A, Appendix G and section 6.1.3.
Its screening equation is (G.35), its leading potential is (G.36), and
its drift scaling is (G.37); these correspond to (G.38)-(G.40) in v3.
The equation and assumptions used here agree in those two editions.
This is a scoped accepted-manuscript comparison, not a claim to have
verified every derivation or the publisher's final typeset version.

All six exact identities, four rational exponents and 16 finite-box solves
pass. The 64 point samples use four locations. Halving the mesh changes
them by approximately one quarter as much on the second refinement; the
largest final change is `9.01e-6`. The outer-box change is `0.006762` at
`h=0.2`, falling to `6.02e-11` at `h=0.025`; the coarse screening case
is retained. These are observed differences, not certified continuum
error bars. The scaled leading and next errors approach finite values
over the tested range, with no inferred uniform numerical constant.
An initial symbolic check left an exact zero in unevaluated radical form
and failed its assertion. Clearing the positive monopole denominators
before exact polynomial factorization resolves it, without a tolerance
or a change to the equation. The failed input and execution record are
retained outside the repository.

Decision: retain this justified baseline and its conditional electric
return budget. MMP already supplied the physical mechanism and leading
scaling; standard elliptic comparison supplies the proof method. There
is no standalone novelty acceptance and no reason to formalize this
baseline before a complete physical consequence survives comparison.
The next decisive calculation is the rotating metric/gauge Dirac
matching with a specified electric state, not more static-port tests.

Remaining to-do list: derive that rotating physical interface and loaded
operator, control full energy and its length variation, finish the
common-observer and primary/fallback novelty review, then independently
verify, GCP-certify suitable complete endpoints and prepare the manuscript.
