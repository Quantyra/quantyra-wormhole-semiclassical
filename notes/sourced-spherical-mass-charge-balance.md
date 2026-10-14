# Sourced spherical mass and charge balance

2026-10-12. This records the exact balance laws that connect the spherical
Einstein-Maxwell constraints to matter flux in the throat and mouths. It uses
the normalization of `exact-spherical-em-reduction.md` and `dirac-spherical-
source-reduction.md`. The result is an identity for the reduced equations,
not an evaluation of the two-ended Dirac state or a stability result.

## Mass aspect

Write `X=R^2`, fix magnetic charge so that the magnetic length is `r_e`, and
define the charged Misner-Sharp mass aspect

```text
m(x) = R/(2G) [1 - (nabla R)^2 + r_e^2/R^2]
     = sqrt(X)/(2G) [1 - (nabla X)^2/(4X) + r_e^2/X].
```

In a source-free RN collar, `(nabla R)^2=1-2GM/R+r_e^2/R^2`, so `m=M`.
The fixed magnetic field is already included in this definition; its energy
must not be counted again as matter. Let `T_ab^(2)` denote all remaining
matter stress after spherical reduction, including the renormalized Dirac
source in a specified scheme, and write `T^(2)=g^ab T_ab^(2)`. The exact
spherical Einstein equations imply

```text
nabla_a m = (T_a^(2)b - delta_a^b T^(2)) nabla_b R.          (1)
```

This is the two-dimensional form of the spherical mass-flux identity. The
normalization follows because `T_ab^(2)` is the sphere-integrated orbit-space
stress (`4 pi R^2` times the four-dimensional orbit-space stress). Equation
(1) follows by differentiating `m`, replacing the second derivatives of `R`
with the metric equation from the exact reduced action and its metric
trace. The angular equation is not needed for this algebraic identity,
as checked in the [dated constraint correction](spherical-anomaly-sources.md).
It is valid before a near-throat
expansion. Diffeomorphism invariance requires the reduced stress and dilaton
source to obey their Ward identity; treating a prescribed stress history
that violates it as a solution would invalidate (1).

For a source-free collar, (1) gives `nabla_a m=0`, recovering the constant
RN mass. With matter present it gives the actual energy transfer that must
be matched to the two asymptotic masses. In particular, the asymptotic mass
change between two cuts is the integral of the corresponding flux, with
orientation and cut normal fixed by the chosen time coordinate. It is not
an independent radiative spherical gravitational mode.

## Electric Gauss balance

If charged Dirac matter is retained, set `F_ab=f epsilon_ab` with
`epsilon_uv=+sqrt(-g)` and define the electric flux
`Q_e=(4 pi/g_c^2) X f`. The Maxwell equation in the source convention of the
mode-reduction note gives

```text
nabla_a(X F^ab) = -(g_c^2/(4 pi)) j^b,
partial_a Q_e = epsilon_ab j^b,                            (2)
```

where the second equality uses the corresponding volume-form convention;
reversing orientation reverses both `f` and the displayed flux sign. Thus
electric flux is constant only in regions with zero radial current. The
fixed magnetic flux remains topological and obeys its Bianchi identity.
The [2026-10-13 correction](spherical-anomaly-sources.md) distinguishes
this scaled charge from `D=X f`, whose derivative has the factor
`g_c^2/(4 pi)`. The predecessor equation attached that factor to Q_e.

## Ward identities required of the retarded source

For a gauge- and diffeomorphism-invariant renormalized effective action
`Gamma[g,X,A]` with the source definitions in the mode-reduction note,
infinitesimal gauge variation and diffeomorphism variation give

```text
nabla_a j^a = 0,
nabla^a T_ab^(2) = (J_X/(4G)) partial_b X + j^a F_ba.      (3)
```

For the second identity, use
`delta g^ab=-2 nabla^(a xi^b)`, `delta X=xi^a partial_a X`, and the
gauge-covariant form of the vector variation
`delta A_a=xi^b F_ba + partial_a(xi^b A_b)`. Integrating by parts and using
the gauge identity removes the pure-gradient term and yields (3). The sign
of the Lorentz-force term follows the convention
`j^a=(1/sqrt(-g)) delta Gamma/delta A_a`; changing the sign in the Maxwell
source convention changes that term consistently. In the fixed zero-electric
sector with no charged source, set `j=0`; one must not do so when the charged
Dirac modes are retained.

The retarded response functions cannot be checked by testing only a
stress-stress correlator. Their metric, area and gauge variations must
satisfy the linearization of (3), including contact terms from varying the
connection, inverse metric, `F_ba`, and reference state. The in-in
functional version assumes the initial state and boundary conditions are
transformed along with the sources; a state held fixed under a physical
boundary deformation requires its explicit variation to appear in the
response. Equations (1)-(3) therefore provide coupled consistency conditions
for the yet-unknown `delta T_ab`, `delta J_X`, and `delta j^a`, rather than
independent adjustable loads.

## Matching consequence and scope

Equations (1)-(2) separate the data needed for a physical calculation:
the throat constraints determine the mass aspect and dilaton; the fermion
state determines stress and current; and integrating their fluxes through
each collar determines the RN mass and electric-flux data at the two ends.
The boundary clocks must then be normalized against those matched exterior
solutions. The vacuum Birkhoff result alone supplies none of the sourced
frequency-dependent kernel. Nor does this balance law determine that kernel:
one still needs the renormalized retarded variations of `T_ab^(2)`, `J_X`,
and `j^a` in the two-ended state, including the massive angular tower,
holonomy, boundary conditions and local counterterms.

The identity is the ordinary spherical Einstein constraint written as a
mass balance, not a new theorem. MMP section 5.3.2 supplies the static
charged-mouth matching baseline; the exact reduced action and source
conventions used here are recorded in the two satellite derivations linked
above.

Remaining to-do list: independently verify (1) from the reduced equations
with unrestricted two-dimensional metric data; evaluate the retarded matter
stress/current response; integrate the balances to both asymptotic mouths;
then reassess whether a physical radial-stability target survives.
