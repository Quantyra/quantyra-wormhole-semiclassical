# Magnetic mouth force balance and orbital hierarchy

2026-10-11, S053. The leading uniform magnetic perturbation of a classical
RN mouth admits an explicit coupled Einstein-Maxwell response. In the
weak-field Ernst/C-metric branch, regularity of both symmetry axes fixes
an acceleration-field balance. Setting that acceleration to zero leaves
a string or strut. Thus the static spherical background is an auxiliary
reference; a physical binary needs force-balanced motion or specified
support, as already recognized by MMP.

The first-order local solution is verified directly below. The selected
weak-coupling family also passes the leading circular-orbit, unscreened
drift and radiation hierarchy checks. This is a reproduction and a model
consistency decision, not a new support, stability or transport theorem.
It does not construct the two-mouth geometry or its quantum state.

## 1. Normalization and the coupled local solution

Use a geometrically normalized field `F_geo=sqrt(kappa/(2g^2)) F_phys`,
where `kappa=8pi G_N`. Then

    G_mu nu = 2(F_mu a F_nu^a - g_mu nu F^2/4),
    P^2 = pi G_N q^2/g^2,        F_0=-P dc wedge dphi,
    f=1-2M/r+P^2/r^2,            c=cos(theta).             (1)

P is the geometric magnetic charge and M the geometric mass. Neither is
the auxiliary port fraction p. The unperturbed metric is
`diag(-f,1/f,r^2/(1-c^2),r^2(1-c^2))`. The physical mass is M/G_N.
For an extremal mouth M=P=r_e>0. The geometric external magnetic field
B has inverse-length units. A remote equal-magnitude monopole at distance
d supplies `|B|=P/d^2` to leading order, with a sign fixed by orientation.

Introduce one bookkeeping parameter multiplying both B and the signed
acceleration parameter A. On a finite classical mouth patch define the
first-order diagonal relative metric corrections by

    h_tt/g0_tt = h_rr/g0_rr = -2c(A r+B P),
    h_cc/g0_cc   = -2c(A r+A M+B P),
    h_phiphi/g0_phiphi = -2c(A r-A M-B P),
    delta a_phi = (B/2)(r^2+P^2)(1-c^2).                 (2)

An irrelevant constant in a_phi has been removed within each monopole
gauge patch; the perturbation in (2) is a globally smooth zero-flux
one-form. Direct construction of the connection, Ricci tensor, inverse
metric, Maxwell stress and Maxwell divergence verifies every background
and first-order Einstein-Maxwell component for arbitrary M,P,A,B.
This is a local electrovacuum calculation away from the axes and horizon
coordinate singularity. It neither supplies a distributional axial source
nor decides its absence until elementary flatness is checked.

The formulas are the weak-field limit of the established magnetic Ernst
response and accelerated charged C-metric, with our signed A convention.
They retain the metric effect on Maxwell. In particular the coupled radial
potential has r^2+P^2, whereas the regular test field on the fixed RN
metric has r^2-P^2. In the frozen metric, the potential in (2) gives

    d_r(f d_r[(B/2)(r^2+P^2)])
           -2[(B/2)(r^2+P^2)]/r^2 = -2B P^2/r^2.       (3)

The metric contribution cancels this residual in the full first-order
Maxwell equation. A small gravitational perturbation cannot simply be
omitted in a near-extremal magnetic mouth merely because the imposed field
is weak: the monopole multiplies that metric perturbation at the same order.
The prior [fixed-metric stress bound](local-maxwell-stress.md) keeps its
stated scope, but is not the solution of this coupled problem.

## 2. Axis regularity and force

Let phi have period 2pi C. Circumference divided by 2pi times proper
distance from an axis is, to first order,

    K_N=C[1+2(A M+B P)],
    K_S=C[1-2(A M+B P)].                                (4)

Allow C=1+delta C. Requiring both K_N=K_S=1 gives

    delta C=0,                 A M+B P=0.                (5)

The sign of A depends on the chosen accelerated coordinate patch. Its
magnitude obeys the established weak-field force law M|A|=|P B|.
At A=0 and BP nonzero, one common azimuthal period cannot regularize both
axes in this branch. Regularizing one produces a deficit or excess at the
other of magnitude `8pi |BP|+O((BP)^2)`. Its tension scale is
`|BP|/G_N`, as expected from the magnetic force. This is not a uniqueness
theorem excluding every static quantum-supported wormhole or additional
external support. It excludes treating this particular unsupported
magnetized mouth branch as smooth and static.

The usual nonrotating magnetic Harrison completion makes the cone check
explicit beyond the linear approximation. Write

    O=1+A r c,     H=1+2MA c+P^2 A^2 c^2,
    Lambda=(1-BPc/2)^2+B^2 r^2 H(1-c^2)/(4O^2).

Its angular metric has
`g_theta theta=Lambda^2 r^2/(O^2 H)` and
`g_phi phi=r^2 H sin(theta)^2/(O^2 Lambda^2)`. Therefore

    K_N=C(1+2MA+P^2 A^2)/(1-BP/2)^4,
    K_S=C(1-2MA+P^2 A^2)/(1+BP/2)^4.                    (6)

For the extremal seed M=P and the small-field branch |BP|<2, equality
of these two cone factors gives exactly

    A=-B/[1+(BP/2)^2],       C=[1+(BP/2)^2]^2.           (7)

These are identities for the stated Ernst-family cone factors. Our
independent bulk verification is through first order, not an all-orders
reconstruction of the Harrison theorem. The period and the physical
magnetic flux also change at second order: in this parameterization the
geometric flux charge is `C P/[1-(BP/2)^2]`. Thus an all-orders fixed-q
comparison must adjust the seed P; at first order it is unchanged.

There is a concrete source correction. Brenner, Giribet and Montecchio,
[published Section VI](https://doi.org/10.1103/PhysRevD.103.124006), display
the static metric with `Lambda=1-BPc+B^2[r^2(1-c^2)+P^2 c^2]/4`.
It gives `C_N=(1-BP/2)^4` and `C_S=(1+BP/2)^4` directly. The expanded
C expressions printed in the following paragraph do not equal that
displayed metric. The PDF page was visually inspected; we use the
direct metric limit, not the inconsistent expansion. The source's
qualitative conclusion that a static monopole needs an axial force agrees.

## 3. Regular extremal response and horizon data

Set M=P and impose the first-order balance A=-B. Equation (2) becomes

    h_tt/g0_tt = h_rr/g0_rr = 2Bc(r-P),
    h_cc/g0_cc = h_phiphi/g0_phiphi = 2Brc.              (8)

In ingoing coordinates `v=t+int dr/f`, the time/radial part is
`[1+2Bc(r-P)](-f dv^2+2dv dr)` through first order. It is regular at
the extremal horizon. The intrinsic horizon metric and magnetic flux
changes are precisely the angular relabeling generated by
`xi^c=-BP(1-c^2)`: their Lie derivatives reproduce both perturbations.
There is no first-order intrinsic horizon distortion from this l=1 field.
Normal derivatives and the exterior field are not removed by that statement.

More explicitly, the relative angular area and magnetic flux densities
at fixed r have coefficients `2Br c` and `B(r^2+P^2)c/P`. Their difference is

    (B/P)(r-P)^2 c.                                     (9)

In the orthonormal frame normal to these r=constant surfaces,

    B_hat(r)=P/r^2+B(1-P/r)^2 cos(theta)+O(B^2),
    B_hat(theta)=-B(1-P/r)sin(theta)+O(B^2).              (10)

The radial perturbation now vanishes quadratically at the extremal
horizon. These formulas refer to the displayed radius gauge and frame;
the difference (9) is invariant under angular relabeling, not under
arbitrary radial gauge changes away from the product limit. Its quadratic
near-horizon power is consistent with the regular mass-six scalar on
Poincare AdS2, for which `(z^2 u')'-6u=0` has solutions z^2 and z^-3.
This consistency check does not reconstruct the full gauge-invariant
master field on a quantum-open throat.

On a bounded central RN neck, the fractional metric correction is
`O(|BP|)=O((P/d)^2)`. On a classical patch extending to r_port=p d it is
at most `C |B| r_port=C(P/d)p`. With p=q^-3/2 and the selected family
below, this is `O(q^-9/2 (log q)^2)`. There is no inverse redshift loss
in (8), but the statement is for the local classical solution. It does
not verify the geometric hypotheses of the earlier stress theorem through
the semiclassical core or the two-mouth exterior.

## 4. The selected family permits the leading orbital hierarchy

MMP v3 [Section 6.1](https://arxiv.org/html/1807.04726v3#S6.SS1) already
uses orbital motion to balance the attraction of opposite mouths. For
two equal extremal mouths in a Newtonian circular orbit, gravitational and
magnetic attraction contribute equally. The relative orbit has

    Omega=2 sqrt(P/d^3),
    |E_orbit|=P^2/(G_N d).                              (11)

The acceleration parameter in (5) accounts for the electromagnetic
force in the local accelerating-mouth reference; it is not the total
coordinate acceleration of the binary. The remote gravitational field
and the matching coordinates also enter that binary problem.

Keep our previously selected weak-coupling family, rather than copying
fixed-g powers of q from a source:

    g^2=q^-2,    ell/P=16pi q^3,    delta=d/ell=(log q)^-2.

Then

    Omega ell=(log q)^3/[2 sqrt(pi) q^(3/2)] -> 0,
    Omega d sqrt(q)=(log q)/[2 sqrt(pi) q] -> 0.          (12)

The second is MMP's naive unscreened orbital displacement relative to the
LLL spacing. This family passes even that leading diagnostic, without
invoking the extra screening suppression discussed in its Appendix G.
It is not a derived bound on the full quantum return or holonomy.

The leading circular-orbit radiation powers, in the same conventions, are

    P_EM=(4pi/g^2)(2/3)(q d/2)^2 Omega^4,
    P_grav=(2/15) P^2 d^4 Omega^6/G_N,
    P_grav/P_EM=(4/5)P/d.                               (13)

In units of one throat crossing and its characteristic excitation energy,

    P_EM ell^2=(log q)^8/(24pi q^2) -> 0,
    P_EM ell/(q/ell)=(log q)^8/(24pi q^3) -> 0,
    (|E_orbit|/P_EM)/ell=24pi^2 q^6/(log q)^6 -> infinity.
                                                               (14)

This passes MMP's stronger total-emission diagnostic as well as its
comparison with the order-q/ell collapse energy. The last expression is
the leading energy-loss timescale; it is not the exact inspiral time.
The radiation estimate is leading classical circular-orbit radiation,
not the energy actually deposited into the throat. No quantitative
nonlinear stability or collapse bound is inferred from it.

The asymptotic inequalities follow directly from powers beating logarithms.
The corresponding first two ratios in (12) decrease for log q>2 and
log q>1; the first ratio in (14) decreases for log q>4. Supplementary
evaluations at q=10^3,10^4,10^6,10^12 retain the constants. At q=10^6,
`Omega ell=7.44e-7`, `Omega d sqrt(q)=3.90e-6`, and
`P_EM ell^2=1.77e-5` after rounding upwards.

These are necessary scale comparisons, not an adiabatic theorem. The
exterior continuum has no uniform spectral gap. Slow orbital motion alone
does not control an evolving vacuum, particle creation, resonant channel
mixing, the renormalized stress, or an unbounded rotation generator in the
full Dirac operator. A co-rotating coordinate choice also does not make a
radiating asymptotically flat binary globally stationary.

## 5. Verification and research decision

Run `checks/magnetic_mouth_force.py --output evidence/magnetic-mouth-force.json`.
The script checks the direct linear tensors, axis limits, extremal balance,
horizon angular relabeling, radiation normalizations and scale identities.
Twelve cone cases solve the balance condition for extremal and nonextremal
seeds, then independently integrate proper angular distance in the stated
nonlinear angular metric at two cutoffs. Their circumference ratios converge
to regular axes. These numerical checks support the angular calculation;
they do not verify the full nonlinear spacetime equations.

The first development run passed all bulk tensor equations and then found
a factor-of-two error in the expected frozen-metric residual. Equation (3)
was independently differentiated and corrected. No tolerance was changed.
The [source review](../evidence/magnetic-mouth-force-source-review.json)
records the exact versions, reviewed derivations and the printed cone
expansion inconsistency. The force law, accelerating response and orbital
proposal are established prior work; no new physical endpoint is accepted.

Continue physical long-return damping, with matched finite-duration band
transfer as fallback. Carry force-balanced moving mouth data in the physical
model; retain the static construction only as a reference whose error must
be bounded. The local l=1 classical balance is now explicit. The decisive
next step is matching it, the companion gravitational/tidal field and higher
multipoles through the semiclassical open throat, followed by a controlled
time-dependent quantum response. This does not justify an unrestricted
inverse estimate based on positive bulk masses.

Remaining to-do list: match the force-balanced binary mouth/metric data
through the quantum-open throat, control time-dependent full Dirac/state
and renormalized response including exterior determinant and discarded tail;
finish literature/proof/novelty review and the consequential endpoint;
independently verify it, certify suitable mathematics exclusively on GCP,
and prepare the reproducible manuscript.
