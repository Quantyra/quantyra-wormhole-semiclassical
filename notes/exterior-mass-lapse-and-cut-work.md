# Exterior mass, lapse and moving-cut work

2026-10-15, E003/S137.3, supporting S074. The spherical equations give
a predictive exterior system once its stress and initial/boundary data
are supplied. They do not determine the finite-cut defect beta from an
asymptotic mass history alone. This note gives the exact system and an
explicit leading free-LLL response on a matched static baseline. Classical
RN-AdS supplies a quantified comparison, not the missing quantum baseline.

The fixed laboratory waveform, eta=1/2000 and original preparation
rectangle remain unchanged. S136 did not certify their capture targets.
No numerical experiment below supplies a physical state or beta bound.

## Exact exterior constraints in the physical frame

On an untrapped exterior patch with monotone area radius R, choose
polar-areal coordinates

```text
ds4^2=-exp(2alpha)F dU^2+dR^2/F+R^2 dOmega^2, F>0,
F=1+R^2/L4^2-2Gm(U,R)/R+r_m^2/R^2,
alpha(U,infinity)=0.
```

U is the fixed asymptotic laboratory time. The fixed magnetic energy
and cosmological term are already included in F. T below is the complete
remaining sphere-integrated orbit stress, T_ab=4pi R^2 T_ab^(4).
Electric energy, if present, is part of this remaining stress unless a
separately consistent charge-work split is made. The exact constraints are

```text
m_R=-T_U^U=:rho2,       m_U=T_U^R,
alpha_R=G*(T_R^R-T_U^U)/(R*F)=G*(rho2+p2)/(R*F).       (1)
```

They follow directly from the four-dimensional Einstein components.
Equivalently the AdS/magnetic-subtracted mass satisfies

```text
m=R/(2G)*[1+R^2/L4^2+r_m^2/R^2-(nabla R)^2],
nabla_a m=(T_a^b-delta_a^b*T_c^c)*nabla_b R.            (2)
```

This is the AdS extension of the retained
[mass balance](sourced-spherical-mass-charge-balance.md), with the same
source normalization. Its energy-supply and work structure is the
standard spherical first law: [Hayward](https://arxiv.org/pdf/gr-qc/9710089),
equations (2.7), (8.9) and (8.14). No new mass theorem is claimed.

Given suitable falloff and consistent sources, the radial reconstruction is

```text
m(U,R)=M_infinity(U)-integral_R^infinity rho2(U,r)dr,
alpha(U,R)=-G*integral_R^infinity [rho2+p2]/(r*F)dr.     (3)
```

The mass evolution in (1), the source Ward identities and initial
constraints must also hold; arbitrary independent functions in (3) need
not solve the full equations. The angular/radius source is still required.
The radial lapse integral is a constraint in this slicing, not by itself
a retarded boundary kernel. Causal prediction additionally requires the
initial state and characteristic matter evolution, followed by null-ray
tracing dR/dU=+/-exp(alpha)F and the appropriate half-density factors.
The static [two-leg LLL graph](exterior-lll-delayed-boundary.md) is its
fixed-geometry limit. Prescribing M(U) inside an otherwise vacuum RN metric
does not provide that driven solution.

## Stored exterior energy and moving-cut pressure work

Let w be a future unit tangent and n the outward unit normal of a timelike
cut. Write dot R=wR, R_n=nR, j=T(w,n), p=T(n,n). Contracting (2) gives

```text
dm_cut/dtau=j*R_n-p*dot R.                             (4)
```

Positive j has the inward-energy-flux convention. The pressure term
cannot be omitted for a moving cut. In polar coordinates, for R=Rc(U),

```text
dm_cut/dU=T_U^R|c+rho2|c*Rc',
E_ext=M_infinity-m_cut=integral_Rc^infinity rho2 dr,
dE_ext/dU=(T_U^R)_infinity-(T_U^R)_c-rho2_c*Rc'.        (5)
```

E_ext is the gravitating exterior mass difference; it is not generally
the proper-volume integral of a local energy density. The limiting boundary
flux includes the actual paired-source work and its prescribed charge
convention. Work at the physical boundary and flux through the internal
cut occur at different times; their difference is stored in the exterior.
Reflection unitarity does not set this storage or source work to zero.

For a fixed areal cut its proper lapse is exp(alpha_c)*sqrt(F_c). For a
moving one it is sqrt(exp(2alpha_c)*F_c-Rc'^2/F_c). Neither is determined
by m alone when the radial pressure is unspecified.

Indeed keep F=F_RN and m constant but choose alpha=epsilon*A(U)*chi(R),
where chi equals one near the cut and vanishes outside a compact exterior
shell. Defining the remaining stress by Einstein's equation gives
rho2=T_U^R=0 but p2=(R*F/G)*epsilon*A*chi'. The angular stress follows
from the same metric and the resulting source is conserved by Bianchi.
The cut lapse changes although all mass histories/fluxes are unchanged.
This is a nonuniqueness witness for unrestricted stresses, not a free
fermion state, causal matter model or physical solution claim. A known
trace or constitutive law can remove this freedom; the free traceless
LLL difference below is one such restricted case.

## Leading free LLL response on a matched static exterior

Suppose a static reference metric/state has actually been matched, with
alpha0=0, F0>0, and optical coordinate
x(R)=-integral_R^infinity dr/F0, so x=0 at infinity.
This restricts the static comparison class: a corrected quantum exterior
can have nonconstant alpha0(R). Its redshift weights and the metric-induced
variation of its reference stress would require a separate calculation;
stationarity of a product throat does not impose alpha0=0 outside it.
Retain only the restricted baseline's
state-dependent, traceless free LLL difference; known reference anomaly,
radius and local terms are separate. At leading order in gravitational
response, canonical propagation gives

```text
H(U,x)=J_plus(U-x)+J_minus(U+x),
delta T_UU=delta T_xx=H,
delta T_Ux=-J_plus(U-x)+J_minus(U+x),
delta rho2=delta p2=H/F0.                              (6)
```

J_plus and J_minus must be derived from the same initial covariance and
paired laboratory source, with the actual propagation phases and both
legs retained. They are not independent reservoirs. Their functions are
not supplied by theta_b alone. An outgoing function with a future
boundary-arrival label U-x is determined by an earlier cut departure;
the label is not permission to use future boundary input.

Equations (1) integrate explicitly at this order:

```text
delta m(U,R)=delta M_infinity(U)-integral_x(R)^0 H(U,xi)dxi,
delta M_infinity'=-J_plus(U)+J_minus(U),
delta alpha(U,R)=-2G*integral_x(R)^0 H(U,xi)/[r(xi)F0(r(xi))]dxi.
                                                               (7)
```

Differentiation gives m_x=H and m_U=-J_plus(U-x)+J_minus(U+x),
so the radial and temporal constraints agree. A boundary action/state
calculation must supply the flux in the second line with its work/contact
normalization; that line is not an extra independent source.

Assume a whole-history chirality envelope |J_plus(v)|+|J_minus(w)|<=Hstar
for every pair of arguments encountered throughout the drive and its
exterior history. A bound on the waveform bandwidth or on H at one
boundary time is insufficient. Then

```text
|delta E_ext|<=Hstar*T_ext,
|delta alpha_c|<=2G*Hstar*I_alpha,
I_alpha=integral_Rc^infinity dr/[r*F0(r)^2],
|delta M_infinity(U)-delta M_infinity(U0)|<=Hstar*|U-U0|.
                                                               (8)
```

These are leading conditional estimates. They do not bound the nonlinear
metric/state feedback, higher angular channels or the omitted contact
terms. The standard dynamic Schwarzian/JT mass relation is already
retained in [S135's boundary energy](jt-prescribed-source-constraint.md);
it is not rederived or counted as new exterior matching here.

## Explicit RN-AdS comparison and the finite-cut warning

For the classical extremal RN-AdS comparison of
[S134](rn-ads-finite-cut-matching.md), put xi=b^2/L4^2, C=1+6xi,
s=(R-b)/b, delta=(Rc-b)/b in (0,1], and

```text
F0=s^2*J(s)/(1+s)^2, J(s)=C+4xi*s+xi*s^2.
I_alpha=integral_delta^infinity (1+s)^3/[s^4*J(s)^2]ds,
P(delta)=(delta^-3-1)/3+(3/2)*(delta^-2-1)
                             +3*(delta^-1-1)+log(1/delta),
I_alpha<=P(delta)/C^2+1/(4xi^2).                        (9)
```

For delta<=s<=1 use J>=C and integrate each term of
(1+s)^3/s^4. For s>=1 use J>=xi*(1+s)^2 and integrate s^-5/xi^2.
This preserves the sharp leading 1/(3C^2 delta^3) coefficient. The
travel-time bound in S134 similarly supplies

```text
T_ext<=(b/C)*[delta^-1+2log(1/delta)-delta+C/xi].        (10)
```

For a transparent dimensionless comparison, let
B=Neta, mu=2GB/b^2, lambda0=kappa0*b/C=mu*z_i/2,
e=lambda0/delta and assume Hstar=A*B*kappa0^2. Equation (9) gives the
explicit bound

```text
|delta alpha_c|<=A*mu*lambda0^2*
                       [P(delta)+C^2/(4xi^2)].         (11)
```

Its leading term is 2A*e^3/(3z_i). At the classical k=2000 comparison,
e=1e-5 and 1/4<=z_i<=7/25, (11) is below 2.667e-15*A.
A is an assumed envelope for the entire stress history, not an assigned
physical coefficient or a consequence of the prescribed source alone.

This small lapse comparison is not a beta bound. Near an ideal JT end,
Phi is of order B/e and the two terms in Phi_sigma+cot(sigma)Phi are
individually of order B/e^2. A relative mass correction O(e^2) contributes
the leading clock charge. An O(e^3) exterior correction can contribute
O(e) to the normalized cut charge. At e=1e-5 that order alone exceeds
the proposed 2e-6 combined ledger. The exact known missing free-collar
energy E_miss in [S136](physical-clock-remainder.md) must cancel before
absolute-value estimation. No proof here leaves an O(e^2) beta remainder;
claiming one from (11) would discard the required charge cancellation.

Moreover the quantum-corrected product radius is not a corrected global
RN exterior. Its displacement can exceed the classical seam delta.
The companion [reference-centering calculation](quantum-reference-clock-centering.md)
quantifies that comparison in its specified one-loop scheme.
It does not determine F0, alpha0, the radial reference stress/radius source,
asymptotic mass subtraction, clock map or initial covariance. Replacing b
in the classical RN formulas without those data is not a valid transfer.
Known constant anomaly shifts must be transferred consistently between
stress and surface charge; raw physical area and reference-centered JT
fields cannot be identified in that step.

## Frozen verification

Before the first checker run, require exact symbolic four-dimensional
Einstein checks of all three constraints in (1), the moving-cut work
identity and the pressure-only nonclosure witness. Use 60 decimal digits.
Independent RN integral variables and the explicit upper bound (9) use
relative tolerance 1e-35 at (xi,delta)=(1,1/10),(17,1/50).
Two smooth chiral verification profiles are checked by quadrature versus
closed primitives, fourth-order finite differences and integrated moving-cut
flux/work; absolute tolerance 1e-35, difference step 1e-10. These profiles
verify formulas and are not claimed quantum states. The candidate comparison
uses exact rational parameter inequalities and (9), with no physical
chirality-envelope assignment. No long nested quadratures, historical test
reruns or Lean execution are needed. Preserve every failed attempt and
record source hashes and the final command/cwd/exit status in UTF8/LF.

Remaining to-do list: obtain the matched quantum exterior/reference state
and boundary source work; solve its coupled characteristic/stress response;
prove the finite-cut charge cancellation and remainder needed for beta;
retain S136's unresolved capture and the physical S074-S079 gates.
