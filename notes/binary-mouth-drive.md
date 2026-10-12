# Circular binary mouth data and the additional electric drive

S053 continuation, 2026-10-11. This calculation identifies the leading
companion data for the [RN dipole transfer](dynamic-rn-dipole.md). It does
not construct a binary wormhole. The circular orbit and retarded-field
methods are established results, not a new physical endpoint.

The useful decision is specific: retain the force-balanced magnetic dipole,
but include the electric dipole in the other perturbation sector and an
explicit angular gluing map. Orbital acceleration cannot be inserted
directly as the acceleration parameter of the local magnetic C-metric.

## Model and conventions

Use geometric Einstein-Maxwell units, with the normalization of
[magnetic mouth force](magnetic-mouth-force.md). Each mass is `P>0`, the
geometric magnetic charges are `P_+=P`, `P_-=-P`, and `d` is the binary
separation in asymptotically inertial harmonic coordinates. Put

```text
x_+(t)= (d/2)n(t),       x_-(t)=-(d/2)n(t),
n=(cos(Omega t),sin(Omega t),0),
tangent=(-sin(Omega t),cos(Omega t),0),
v=Omega d/2,            chi=P/d << 1.
```

We use the conservative, nonspinning charged point-body approximation for
the orbit. A common electromagnetic duality rotation maps the purely
electric two-body problem to this purely magnetic one. Finite mouth size,
the wormhole stress and radiation reaction require separate matching.
The field calculation below uses a *prescribed circular orbit in flat
spacetime*. Its exact retardation is an independent check of leading local
data; its higher coefficients are not complete post-Newtonian GR fields.

## Orbit, local acceleration and clocks

The Newtonian action has potential term `(m1 m2-Q1 Q2)/d`; here the two
attractions are equal. Thus

```text
Omega^2=4P/d^3 [1+O(chi)],       v^2=chi [1+O(chi)],
x_+''=-2P n/d^2 [1+O(chi)],    E_orbit=-P^2/(G d) [1+O(chi)].       (1)
```

In the regular external gravitational background, `g00=-1+2U+...`,
`U=P/|x-x_-|`, and `Gamma^i_00=-partial_i U`. At the positive body,
`grad U=-P n/d^2`. Consequently its covariant acceleration relative to
that background, to this order, is

```text
a_regular=x_+''-grad U=-P n/d^2=(P_+/P) B_companion,
B_companion=-P n/d^2.                                             (2)
```

The negative body's covariant acceleration has the opposite sign. These
are effective center worldlines in the regular background, not clocks or
observers placed on a black-hole horizon. Equation (2) is the proper
force that enters the local accelerating mouth solution. In the signed
C-metric convention of the earlier note, `A_C=-a_regular` along the
chosen axis and `P A_C+P_+ B=0`. Using the full coordinate acceleration
would double this electromagnetic acceleration and spoil axis balance.

The regular-background clock rate of either center is

```text
d tau_regular/dt=1-U-v^2/2+O(chi^2)=1-3chi/2+O(chi^2).              (3)
```

Their equality follows from the binary symmetry at this order. Absorbing
the constant companion potential in a local time coordinate does not
remove (3) from a comparison with the exterior clock. The RN self-field
and the quantum-open throat introduce additional lapse and matching data.

For reference, the established 1PN Lagrangian gives a sharper orbital
check. With `G=c=1`, its velocity-independent 1PN terms are

```text
-m1 m2(m1+m2)/(2d^2) + Q1 Q2(m1+m2)/d^2
                         -(m1 Q2^2+m2 Q1^2)/(2d^2).              (4)
```

Specializing the full Lagrangian to equal opposite charges and a circular
orbit, omitting the constant rest masses, gives

```text
Lcirc=P d^2 Omega^2/4+2P^2/d
      +P d^4 Omega^4/64+(7/4)P^2 d Omega^2-4P^3/d^2,
Omega^2=(4P/d^3)[1-8chi+O(chi^2)],
G E_orbit=-P^2/d+(15/4)P^3/d^2+O(P^4/d^3).                        (5)
```

Differentiate `Lcirc` with respect to `d` at fixed `Omega`, then compute
`Omega partial_Omega Lcirc-Lcirc`. A separate Cartesian variation of the
two-body action, performed before imposing circular motion, gives the
same radial correction. These are conservative PN expansions; no bound
on the omitted GR remainder is proved here.

The adopted source is [Khalil et al., v2, (3.1), appendix A](https://arxiv.org/pdf/1809.03109v2).
Their zero-dilaton specialization has the positive mixed term in (4).
[Tang, Liu and Han, v2, (11), PDF p.4](https://arxiv.org/pdf/2512.04806v2)
prints the opposite sign; visual inspection confirms the PDF/HTML agree.
The [published equation (11)](https://link.springer.com/article/10.1140/epjc/s10052-026-15381-0)
also displays that sign in its current HTML.
We reject that printed term for the following two independently
reproducible reasons, without inferring that every result of that paper
is affected:

* Setting `Q_i=m_i` and velocities to zero cancels (4), as extremal static
  balance requires. The alternate printed sign instead leaves
  `-2m1 m2(m1+m2)/d^2`, a separation-dependent static force.
* In exact RN, let areal radius `R=rho+M`, with `rho` harmonic. The static
  charged test-particle action is `-m sqrt(1-2M/R+Q^2/R^2)-eQ/R`.
  Its expansion is `-m+(mM-eQ)/rho +[-m(M^2+Q^2)/2+eQM]/rho^2`.
  Equation (4) has precisely this test-body limit; the alternate sign
  does not. This test needs no binary PN field solution.

## Retarded companion fields in the instantaneous rest frame

At `t=0`, let `theta=Omega Rret`. The light cone of the positive mouth and
the retarded negative companion satisfies

```text
theta=2v cos(theta/2),    Rret=d cos(theta/2),
nret=(cos(theta/2),-sin(theta/2),0),
beta_source=(-v sin theta,-v cos theta,0),
dot_beta_source=(2v^2/d)(cos theta,-sin theta,0),
k=1-nret.beta_source=1+v sin(theta/2).                              (6)
```

For `0<=v<1`, the root in `[0,2v]` is unique: the defining function has
derivative `1+v sin(theta/2)>0`. It has the expansion
`theta=2v-v^3+(13/12)v^5+O(v^7)` by direct substitution.

The dual Lienard-Wiechert fields of magnetic charge `Q=-P` are

```text
B=Q [(1-beta^2)(nret-beta)/(k^3 Rret^2)
       +nret cross ((nret-beta) cross dot_beta)/(k^3 Rret)],
E=-nret cross B.                                                   (7)
```

These follow by electromagnetic duality from the retarded potentials;
[Fitzpatrick's derivation](https://farside.ph.utexas.edu/teaching/jk1/Electromagnetism/node155.html)
keeps both source-time derivatives and the acceleration term. For an
observer moving with `u=v tangent`, use

```text
E'=gamma(E+u cross B)-gamma^2/(gamma+1) u(u.E),
B'=gamma(B-u cross E)-gamma^2/(gamma+1) u(u.B).                     (8)
```

In the axes obtained by the instantaneous rotation-free boost from the
inertial frame, (6)-(8) imply at the positive mouth

```text
B'/(P/d^2)=-(1+5v^4+O(v^6))n
           +[-(8/3)v^3+(112/15)v^5+O(v^7)]tangent,
E'/(P/d^2)=(2v+2v^5+O(v^7)) e_z.                                 (9)
```

The cancellation of the `v^2` radial magnetic coefficient and `v^3`
electric coefficient is a property of this flat prescribed-orbit model.
Curvature corrections enter already at relative `chi~v^2`, and the
source-only retarded field excludes the body's own regular radiation
field. Equation (9) therefore supplies neither physical GR 1PN field
coefficients nor a radiation-reaction force. Its defensible leading data
are

```text
B'_+=B'_-=-P n/d^2 [1+O(chi)] + subleading components,
E'_+=+2vP e_z/d^2 [1+O(chi)],
E'_-=-2vP e_z/d^2 [1+O(chi)].                                    (10)
```

The indicated orders are PN ordering statements, not rigorous uniform
error bounds for a constructed binary metric. The electric field contains
one contribution from the moving source and one from the mouth boost.
Omitting either loses a factor of two. It is a uniform electric `l=1`
field about each mouth, so it is absent from the previously reconstructed
magnetic scalar branch. In a magnetic RN background it also couples to
the axial metric sector. At leading order its direction is the fixed
orbital normal: it is a stationary drive in these axes, whereas the
magnetic coefficient rotates at `Omega`. The missing calculation is
therefore initially a stationary electric/axial dipole, including its
angular-momentum integration data. Duality would rotate the background too and
cannot simply turn this new perturbation into the previous solution
while holding the magnetic background fixed.

Instantaneously boosted axes and Fermi-Walker transported axes are not
identical along an accelerating orbit. Their difference starts beyond
the leading field data. This calculation has not fixed the transported
mouth tetrads, their spin lift, or the fermion return holonomy.

## Finite ports and the magnetic master datum

Take a port radius `a=p d`, with `P<<a<<d`. For the quasistatic flat
companion, put `y=(x-x_+)/d`. The dimensionless vector field is
`(n+y)/|n+y|^3`. Its derivative has operator norm at most
`2/(1-p)^3` on `|y|<=p`; hence

```text
|B(x)-B(x_+)|/(P/d^2) <= 2p/(1-p)^3.                              (11)
```

This is an actual bound for that static field. Moreover its magnetic
scalar potential is harmonic inside the port. Its degree-one solid
harmonic is exactly its gradient at the center, by the convergent
Legendre expansion of `1/|d n+x|`. Projecting onto `l=1` therefore gives
the uniform-field coefficient exactly, at every `a<d`; the gradients
populate other multipoles. Neither this exact projection nor (11)
eliminates those multipoles from the coupled physical problem.

For the *regular external background*, accelerated Fermi coordinates
have `g00=-(1+a_regular.x)^2-R0i0j x^i x^j+O(|x|^3)` and the corresponding
quadratic `g0i,gij` terms; see [Poisson, Pound and Vega, section 9, (9.14)-(9.16)](https://arxiv.org/pdf/1102.0529v3).
The leading companion curvature is `O(P/d^3)`. Thus the port acceleration
and gravitational tidal scales are respectively `chi p` and `chi p^2`.
This is a Taylor/PN scale assessment, not a uniform bound on a full
RN-plus-binary junction. A rotation of the axes adds its own connection
terms. The strong self-field remains in the RN inner problem.

In the notation of the earlier master calculation, the local magnetic
drive is the scalar `epsilon_i(omega)=P_i B'_i.omega`. Thus in two mouth
frames aligned initially with the common inertial axes,

```text
epsilon_+=-chi^2 n,        epsilon_-=+chi^2 n                       (12)
```

as coefficient vectors. The static regular-axis solution has
`W=Z/Z0=(3/2)epsilon_i`. This identifies the leading *magnetic* outer
datum of the finite-annulus transfer, rather than deriving all boundary
data of a smooth throat.

Let the negative sphere be glued into the positive throat angular chart
by `omega_-=O omega`, where `O` is an orthogonal matrix with determinant
`-1`. This makes the opposite oriented monopole fluxes consistent. The
scalar coefficient pulled back from that mouth is

```text
epsilon_-^common=+chi^2 O^T n.                                   (13)
```

For an explicitly chosen antipodal gluing `O=-I`, (12)-(13) agree. The
previous prescribed master junction then yields, componentwise in `l=1`,

```text
D=-(3/2)lambda^2 chi^2 [1+O(s+R^-2+eta)] n                       (14)
```

at both ends, to leading magnetic order. This uses the already proved
bound for that junction and its definitions `s=lambda R` and `eta`;
binary corrections and the other sectors are additional data, not
included in this error expression. For general `O`, the two drives differ
as (13) specifies. The local force law does not select `O`, its time
dependence or the spinor phase. The antipodal choice is a transparent
model datum, not a derived quantum holonomy.

## Consequence for the retained hierarchy and verification

For the existing family `lambda=1/(16pi q^3)`, `d/ell=(log q)^-2`,
`p=q^-3/2`, all the following necessary scales tend to zero:

```text
chi=(log q)^2/(16pi q^3),      |epsilon|=chi^2,
|E'|/|B'|~2v=(log q)/(2sqrt(pi) q^1.5),
Omega ell=(log q)^3/(2sqrt(pi) q^1.5),
clock correction~(3/2)chi,     acceleration port scale=chi p,
tidal port scale=chi p^2.                                        (15)
```

A prescribed rigid rotation of the LLL spin `j=(q-1)/2` has generator
norm `j|Omega|`. Over a time `ell`, its bound is
`j Omega ell=O(q^-1/2(log q)^3)`, which also vanishes but need not be
small at modest flux. This reuses the representation in
[LLL return holonomy](lll-return-holonomy.md). It is only a candidate
rigid-frame estimate: it is not the norm of the actual many-channel
return, does not control the full Dirac spectrum, and proves no adiabatic
vacuum or Casimir-energy bound. Choosing a rotating basis alone cannot
determine a physical holonomy.

[The checker](../checks/binary_mouth_drive.py) verifies (4)-(9)
symbolically, including the independent Cartesian orbit and exact RN
test limit. Separately it solves the retarded light cone at displaced
spacetime points, differentiates the scalar/vector potentials with a
fourth-order finite difference, and compares those fields with (7)-(8)
at several speeds and phases. Three step sizes test convergence. It
also checks the Lorentz invariants, light-cone root, degree-one port
projection and the hierarchy. The [receipt](../evidence/binary-mouth-drive.json)
pins the script and input notes. Numerical samples do not certify the
PN remainder or replace the ordinary bounds above.

The [source review](../evidence/binary-mouth-drive-source-review.json)
records exact editions, read scope, cache hashes and the rejected sign.
No novelty, full junction, finite-disturbance stability, state response or
shortcut conclusion follows from these component checks. No Lean or
cloud run is justified by this established classical calculation alone.

Remaining to-do list: solve the stationary electric/axial response and
remaining multipoles through a smooth quantum-open junction; derive
transported frames, gluing and physical return holonomy; control the
time-dependent quantum state, renormalized source and exterior response;
finish the consequential endpoint and novelty review, independently
verify, GCP-certify suitable mathematics and prepare the manuscript.
