# Retarded LLL stress and mass flux in a fixed RN collar

2026-10-12. This derives the complete two-dimensional conformal stress
tensor of the massless spherical LLL sector on a fixed, source-free RN
collar, and composes its state map with a moving throat/RN matching sphere.
It gives an exact collar mass-flux law for the state-dependent chiral
response. It does not solve the moving RN geometry, determine the boundary
state from the full exterior, or include massive angular modes.

## Conformal LLL on the RN collar

Take a magnetic RN mouth with fixed flux `q`, radius `r_e`, and

```text
ds4^2=-F(r)dt^2+F(r)^(-1)dr^2+r^2dOmega_2^2,
F(r)=1-2GM/r+r_e^2/r^2,
dr_*/dr=1/F(r),  u=t-r_*,  v=t+r_*.
```

In the source-free MMP collar the mass is constant and the chosen below-
extremal mass shift makes `F>0`; there is no horizon in this patch. After
the spherical spinor rescaling, the `q` LLL channels are massless two-
dimensional Dirac fields. Their action is Weyl invariant on
`ds2^2=-F du dv`; in this restricted sector the static RN collar has no
scattering potential and transports the two chiral state functions
independently. The (O(r_e/ell)) matching corrections, gauge polarization,
and all nonzero angular modes are excluded here.

Let `R_+(u)` and `R_-(v)` be the monotone state maps fixed by the past
two-ended state and propagated through the collar. For total central charge
`c=q` and zero background electric potential, the chiral state functions
are

```text
t_+(u)=R_+'(u)^2 t_+^0(R_+(u))-c/(24pi) Sch(R_+,u),
t_-(v)=R_-'(v)^2 t_-^0(R_-(v))-c/(24pi) Sch(R_-,v).        (1)
```

The full covariantly renormalized stress in the fixed RN metric is

```text
T_uu=t_+(u)+c[2F F''-(F')^2]/(192pi),
T_vv=t_-(v)+c[2F F''-(F')^2]/(192pi),
T_uv=c F F''/(96pi).                                      (2)
```

Primes on `F` mean `d/dr`. These terms follow from
`sigma=(1/2)log F` in the standard two-dimensional conformal-anomaly
stress. They give
`T^a_a=-c F''/(24pi)=c R2/(24pi)` because `R2=-F''`, and obey covariant
conservation. Equation (1) is the state-dependent retarded response;
the local terms in (2) are fixed by the RN background. They must be
included in the static reference stress and not counted again as a
frequency-dependent matter load.

For a vector current, the two covariant null components are independently
chiral and transform with weight one,

```text
j_u(u)=R_+'(u) j_+^0(R_+(u)),
j_v(v)=R_-'(v) j_-^0(R_-(v)).                              (3)
```

There is no Schwarzian term in the anomaly-free vector current. The
antiperiodic, zero-holonomy-current reference has `j_+^0=j_-^0=0`; a pure
metric perturbation then has no vector current in this fixed-gauge LLL
model. If a radial gauge perturbation or nonzero current is present, (3)
must instead be coupled to the two-dimensional Maxwell Gauss constraint;
setting the electric flux to zero is not consistent in that case.

## Matching the retarded state map at a moving sphere

Let the matching sphere have areal radius `R_b(t)` in the RN Killing time.
The RN null coordinates at the interface are

```text
u_b(t)=t-r_*(R_b(t)),   v_b(t)=t+r_*(R_b(t)),
u_b'=1-R_b'/F(R_b),   v_b'=1+R_b'/F(R_b).                  (4)
```

Timelikeness is `|R_b'|<F(R_b)`, so both maps remain increasing. If the
throat outgoing null phase on the same interface is `Y_+(U_b(t))`,
continuity of the LLL characteristic gives
`R_+(u_b(t))=Y_+(U_b(t))`; the other side similarly gives
`R_-(v_b(t))=Y_-(V_b(t))`. Therefore, for example,

```text
R_+'(u_b)=Y_+'(U_b) U_b'/[1-R_b'/F(R_b)].                  (5)
```

The Schwarzian chain rule supplies the corresponding second/third
derivative contact in (1). More explicitly, with
`Z(lambda)=Y_+(U_b(lambda))` and `u_b=u_b(lambda)`,

```text
Sch(Z as a function of u_b)
  = (U_b'/u_b')^2 Sch(Y_+,U_b)
    + [Sch(U_b,lambda)-Sch(u_b,lambda)]/(u_b')^2.           (6)
```

Thus the RN chiral state stress at the cut is

```text
t_+^RN(u_b)=(U_b'/u_b')^2 t_+^th(U_b)
 - c/(24pi) [Sch(U_b,lambda)-Sch(u_b,lambda)]/(u_b')^2,     (7)
```

where `t_+^th=Y_+'^2 t_+^0-c Sch(Y_+,U_b)/(24pi)` is the throat state
function. The corresponding vector current has only its weight-one factor,
`j_u^RN=(U_b'/u_b')j_U^th`. For a static matching sphere, equation (4)
is a constant tortoise translation and the stress acquires only the
constant redshift factor `x_*^2`. A moving sphere contributes through both
the redshift and the Schwarzian contact in (7), once `R_b(lambda)` and the
induced time matching have been solved from the coupled Einstein equations.
Equations (6)-(7) determine the interface transform but do not specify that
dynamics or the local renormalized join action.

## RN mass function sourced by the conformal stress

The mass aspect in the spherical balance note obeys
`nabla_a m=(T_a^b-delta_a^b T)nabla_b r`. In the static RN null chart,
`r_u=-F/2`, `r_v=F/2`, and `g_uv=-F/2`. Substitution gives the exact
collar identities

```text
partial_u m=-(T_uu+T_uv),
partial_v m= +(T_vv+T_uv).                                  (8)
```

For the finite state difference on the same fixed metric, the geometric
terms cancel, `Delta T_uv=0`, and

```text
partial_u Delta m=-Delta t_+(u),
partial_v Delta m=+Delta t_-(v),
Delta m=C-int^u Delta t_+(s)ds+int^v Delta t_-(s)ds.       (9)
```

Here `C` is fixed by the RN mass data on an initial cut. Equation (7)
identifies how retarded LLL energy changes the collar mass aspect; it is
the missing sourced boundary datum that vacuum Birkhoff cannot provide.
It remains to compose the two mouth maps with the actual field-line return
and solve for `R_b(t)` so that (5), (6), the throat constraints and the
asymptotic clock normalizations hold simultaneously.

Christensen--Fulling's anomaly/conservation method is the established
two-dimensional stress construction; see their primary
[1977 paper](https://doi.org/10.1103/PhysRevD.15.2088). The RN specialization
and the characteristic composition above are ordinary applications to the
MMP LLL sector. The baseline mouth mass and static clock normalization are
reproduced in `rn-jt-static-boundary-charge-match.md` and MMP
[equations (5.39)-(5.42)](https://arxiv.org/html/1807.04726v3#S5.SS3).

Run `checks/rn_collar_lll_flux.ps1` for deterministic checks of the RN
anomaly trace/conservation identities, null mass-balance projections, the
timelike interface-map signs, and the composed Schwarzian formula (6).
These checks test the formulas on the fixed RN background; they do not
verify the four-dimensional renormalized state, moving-geometry response,
or stability.

Remaining to-do list: determine the matched initial state maps and moving
RN sphere from the coupled boundary equations; include electric
polarization and massive modes; solve the full two-mouth mass/clock
response and reassess the stability target.
