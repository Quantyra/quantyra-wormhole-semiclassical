# Spherical electrovac mouths: the l=0 matching data

2026-10-12, consequence of the S051 constraint and competitor audits. In
the source-free Reissner-Nordstrom (RN) mouth approximation used by MMP,
the spherically symmetric gravitational/electromagnetic exterior has no
local radiative degree of freedom. Its matching data are conserved mass,
fixed magnetic charge (and, if allowed, a conserved electric charge), plus
coordinate gauge. This removes a proposed vacuum gravitational-wave load
from the monopole problem; it does not remove the dynamical fermion
transport or its energy flux.

## Local exterior solution space

In an open mouth region neglect the quantum source to the same order as
the MMP classical near-extremal matching and impose spherical symmetry.
The magnetic flux is fixed and quantized. The source-free spherical
Einstein-Maxwell equations imply the RN form

```text
ds^2=-F(r)dt^2+F(r)^(-1)dr^2+r^2 dOmega_2^2,
F(r)=1-2GM/r+Q_m^2/r^2,
```

up to radial/time coordinate transformations. The integration constant
`M` is constant in the vacuum region by the spherical Einstein constraints;
the magnetic charge `Q_m` is fixed by the Bianchi identity and flux sector.
A spherically symmetric electric field is another two-dimensional Maxwell
constraint, not a propagating photon; if charged currents are admitted, its
Gauss law must be retained, and otherwise its conserved charge is fixed.

Linearizing within this exact family gives
`delta F=-2G delta M/r + 2 Q_m delta Q_m/r^2`. In the fixed magnetic and
zero-electric-charge sector, `delta Q_m=0` and `delta M` is constant in
time. A lapse perturbation left by a time-dependent radial gauge choice is
removed by a reparametrization of `t`; it is not an independent local
degree of freedom. Thus for Fourier frequency `omega != 0`, the source-free
mouth carries no spherical gravitational perturbation. At `omega=0`, a
mass shift is allowed and is fixed by the global energy/matching condition.

This statement applies only to the classical source-free mouth region. If
the Dirac stress or current is retained there, the mass/charge functions
obey sourced constraints and the fermion propagation must be included; one
cannot use Birkhoff's theorem to discard that matter load. In a matched
calculation the throat clock cannot be extended as an arbitrary radiative
gravitational wave through a vacuum RN collar. It must be supported by a
constraint-compatible throat solution and by the stress/current/energy
transfer that changes or matches the exterior conserved charges.

## Consequence for the target and prior work

MMP's [equations (5.32)-(5.33)](https://arxiv.org/html/1807.04726v3#S5.SS3)
use the exact static RN mouth and retain its mass parameter for throat
matching. Sadhukhan's [2026 perturbation analysis](https://arxiv.org/html/2609.33511v1)
states explicitly that its `l=0` mouth perturbation reduces by generalized
Birkhoff to mass plus gauge, with quantized magnetic charge unchanged; it
does not treat the resulting throat JT boundary problem. These are
established facts, not a new no-go.

The primary target should therefore distinguish three pieces: (i) the
nonpropagating throat dilaton constraints, (ii) asymptotic mass/time
matching data in each vacuum RN mouth, and (iii) the retarded charged-fermion
state and energy/current flux through the throat and the ambient exterior.
Only (iii) can produce a frequency-dependent matter load in a vacuum-matched
`l=0` calculation. A projected absorbing port or a prescribed metric clock
does not establish that load. The leading transparent-channel model remains
a comparator, not the full exterior solution.

Remaining to-do list: derive the sourced spherical charge/mass balance
through the nonvacuum throat and fermion-filled collars; match that balance
to the two asymptotic RN masses and normalized clocks; then insert the
resulting state-dependent matter kernel into the JT boundary equation and
reassess pole stability.
