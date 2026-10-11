# Classical dipole matching into the leading throat

2026-10-11, S053. The force-balanced RN dipole can be put in a gauge that
matches the physical mass-six scalar of the leading AdS2 x S2 throat.
The required radial gauge displacement changes the matching coefficient:
removing the angular area perturbation alone is insufficient. The complete
four-dimensional classical metric/Maxwell reconstruction and the throat's
driven response are derived below and checked independently.

This closes a specific leading classical scalar matching calculation. It
does not match the full moving binary, determine relative boundary frames
or holonomy, include the semiclassical source response, or construct the
quantum state. The original consequential research goal remains incomplete.

## 1. A common bulk gauge for the force-balanced mouth

Use the geometric normalization of [the force calculation](magnetic-mouth-force.md).
Take M=P=r_e, signed external magnetic field B, and A=-B through first
order. Let

    z=(r-P)/P,     epsilon=B P,     f=z^2/(1+z)^2.

Factor epsilon out of all perturbations and temporarily set P=1. In
coordinates (t,z,c,phi), c=cos(theta), the known coupled mouth response is

    h_tt/g0_tt=h_zz/g0_zz=2 z c,
    h_cc/g0_cc=h_phiphi/g0_phiphi=2(1+z)c,
    delta A_phi=[((1+z)^2+1)/2](1-c^2).                  (1)

The background gauge potential is A0_phi=-c. Seek a coordinate change
`h_new=h-L_xi g0`, `delta A_new=delta A-L_xi A0` with

    xi^z=U(z)c,
    xi^c=[U(z)/(1+z)-(1+z)](1-c^2),    xi^t=xi^phi=0.   (2)

The angular metric perturbation vanishes exactly. Removing the z-c
component as well requires

    U'+(1+2z)U/[z^2(1+z)]=1+z.                         (3)

The regular solution at z=0 is

    I(z)=int_0^infinity exp(-s)/(1+zs)^3 ds,
    U=z^2(1+z) I(z),
    z^2 I'+(1+3z)I=1.                                 (4)

For example, the integrating factor of (3) is
`z exp(-1/z)/(1+z)`. Its homogeneous solution grows as
`(1+z)exp(1/z)/z`; regularity excludes it. The integral representation
in (4) gives the unique regular solution, rather than just a formal series.
It follows by changing variables in
`U=(1+z)exp(1/z)/z int_0^z v exp(-1/v) dv`.

The new dimensionless magnetic amplitude and relative time/radial metric
amplitudes, with epsilon restored, are

    a_RN=epsilon z^2[1/2+I(z)],
    H_t=2epsilon z[1-I(z)],
    H_z=2epsilon[(1+3z)I(z)-1],
    delta A_new=P a_RN sin(theta)^2 dphi.                (5)

Thus, as z tends to zero,

    a_RN=(3/2)epsilon z^2[1+O(z)],
    H_t=H_z=6epsilon z^2[1+O(z)].                       (6)

These estimates have direct remainders. Taylor's theorem inside the
positive integral gives, for z>=0,

    0<I<=1,   I>=1-3z,
    |I-(1-3z)|<=12z^2,
    |I-(1-3z+12z^2)|<=60z^3.                           (7)

They control (6) uniformly on a fixed small interval. Keeping just the
angular shift would instead give a coefficient 1/2 in a_RN and leave a
mixed metric component. The radial displacement in (2) contributes the
other unit to the coefficient 3/2.

This is a bulk gauge choice with nonzero boundary values. The coordinate
change and monopole gauge transition must be carried into any boundary
matching of clocks, frames and charged spinors. We do not declare those
boundary data physically irrelevant.

## 2. Direct reconstruction on the leading product throat

Scale out P and use

    g0=-(1+rho^2)d tau^2+d rho^2/(1+rho^2)+dOmega^2,
    F0=-dc wedge dphi,           f=1+rho^2.

For the dipole Y=c define the perturbation by

    h_ab=4 a(tau,rho) Y g0_ab,       a,b in {tau,rho},
    h_aA=h_AB=0,
    delta A=a(tau,rho) sin(theta)^2 dphi.                (8)

Every linear Einstein-Maxwell equation is satisfied if

    (Box_AdS2-6)a=0.                                   (9)

This is checked by direct four-dimensional tensor contraction, including
mixed time/angular and radial/angular components. The time dependence
therefore includes the induced electric field; it is not obtained by
changing a static magnetic amplitude while freezing Maxwell's equations.
For `a=exp(-i sigma tau) u(rho)`, (9) reads

    d_rho(f u')+[sigma^2/f-6]u=0.                        (10)

Rotational covariance gives the other two l=1 harmonics, and linear
superposition describes a rotating dipole in this classical scalar branch.
The other parity branch and boundary vector degrees are not included.

The physical mass six is established in Michelson-Spradlin,
[hep-th/9906056v3](https://arxiv.org/abs/hep-th/9906056v3), Section 3.7.
In its normalization, our bulk gauge has pi=0, b=2a, H_ab=4a g_ab;
its scalar pi-2b is -4a. We supply the explicit reconstruction (8), rather
than assigning physical matching to a mass eigenvalue alone. The source's
extra boundary modes remain outside (8).

The RN dipole and its metric/Maxwell reconstruction are also established
in Burko's [polar dipole paper](https://arxiv.org/abs/gr-qc/9507012v1).
Its equations (40), (41) and (87) give, for the physical dipole master,

    (d_rstar^2+omega^2) Z=V Z,
    V=Delta/r^5 [2M-2Q^2/r+6M W],
    W=3M Delta/(r varpi^2)+M/varpi,
    varpi=3M-2Q^2/r,  Delta=r^2-2Mr+Q^2.

At M=Q=P, the limit of V/f at r=P is 6/P^2, with static throat
exponents 2 and -3. The checker verifies this independent predecessor
limit. The paper uses an electric background; electromagnetic duality
relates the magnetic sector, but its full reconstruction variables and
boundary gauges have not been identified with ours. Equations (8)-(10)
are checked directly rather than importing that uncompleted identification.

## 3. Exact low-frequency response to both end coefficients

Put x=atan(rho). Equation (10) is

    [-d_x^2+6 sec(x)^2]u=sigma^2 u.                     (11)

Two Darboux steps applied to exp(i sigma x) give the exact complex solution

    [1-sigma^2+3rho^2+3i sigma rho] exp(i sigma x).

Take its real and imaginary parts,

    E_sigma=(1-sigma^2+3rho^2)cos(sigma x)
                                    -3sigma rho sin(sigma x),
    O_sigma=(1-sigma^2+3rho^2)sin(sigma x)
                                    +3sigma rho cos(sigma x).    (12)

Specify the two leading growing coefficients
`u(rho)/rho^2 -> D_+` at +infinity and `-> D_-` at -infinity.
For real |sigma|<=1/2, the unique solution in this frequency sector is

    u=(D_++D_-)/(6 cos(pi sigma/2)) E_sigma
           +(D_+-D_-)/(6 sin(pi sigma/2)) O_sigma,      (13)

with the continuous limit at sigma=0. That limit is

    u_0=(D_++D_-)(1+3rho^2)/6
      +(D_+-D_-)[(1+3rho^2)atan(rho)+3rho]/(3pi).       (14)

In particular,

    u(0)=(1-sigma^2)(D_++D_-)/(6 cos(pi sigma/2)),
    u'(0)=sigma(4-sigma^2)(D_+-D_-)/(6 sin(pi sigma/2)). (15)

There is no small-frequency pole. For any fixed core |rho|<=R0,
the amplitude and its first two rho derivatives are bounded by
`C_R0 (|D_+|+|D_-|)`, uniformly in this frequency interval. Their difference
from (14), at the same instantaneous end amplitudes, is bounded by
`C_R0 sigma^2 (|D_+|+|D_-|)`. These statements follow directly from the
even analytic quotients in (13), using O_sigma/sigma at zero.

Uniqueness does not rely only on the explicit denominators. For a
normalizable difference of two solutions, the operator in (11) factors as

    H=(d_x+3tan x)^*(d_x+3tan x)+9.                     (16)

The boundary term vanishes on its normalizable form domain, and
`cos(x)^3` attains the lower bound nine. Hence no nonzero normalizable
difference exists at |sigma|<=1/2. The corresponding normalizable
resolvent norm is at most `1/(9-sigma^2)`. This is the classical dipole
operator's gap, not a gap of the full exterior Dirac continuum.

Equations (13)-(15) specify the harmonic particular response with no added
free oscillations. Time-dependent boundary values alone do not eliminate
homogeneous initial excitations. We have not proved that a physical
binary selects this solution by adiabatic preparation or a retarded vacuum.

## 4. Matching coefficient and the physical scales

Let lambda=P/ell. The established RN/throat overlap relation is
`z=lambda |rho|`, with `1<<|rho|<<lambda^-1`. For each mouth, (6) supplies
the leading static coefficient

    D_+=(3/2)epsilon_+ lambda^2,
    D_-=(3/2)epsilon_- lambda^2.                        (17)

The signs and harmonic frames are part of the supplied mouth data;
opposite magnetic charges do not by themselves fix D_+=D_-. With equal
coefficients epsilon_+=epsilon_-=epsilon, the static response is
`u_0=(epsilon lambda^2/2)(1+3rho^2)`.

There is also a finite-overlap check. Solve (10) on [-R,R] with
`u(+R)=D_+ R^2` and `u(-R)=D_- R^2`. The leading r terms cancel in
(12), so its exact infinite-end solution differs from these boundary
values by O(|D_+|+|D_-|). The even/odd boundary matrix has entries of
order R^2 with nonzero coefficients uniformly for |sigma|<=1/2.
Consequently its core solution differs from (13) by
`O((|D_+|+|D_-|)/R^2)`, also for two core derivatives.
The static RN coefficient (17) has a relative remainder O(lambda R).
Choosing `R=lambda^-1/3` makes these two classical comparison errors
O(lambda^(2/3)). This is not a bound on the omitted semiclassical equations.

For a leading external field |B|=P/d^2, `|epsilon|=(P/d)^2`.
The supplied classical dipole data therefore produce a core metric
perturbation of order

    (P/d)^2(P/ell)^2=lambda^4/delta^2,
    delta=d/ell.                                      (18)

In the selected family `lambda=1/(16pi q^3)`,
`delta=(log q)^-2`, this is `O(q^-12 (log q)^4)`. Relative to the
specified leading Casimir opening scale `alpha=g^2/(4pi^2 q)`, it is
`O(q^-9 (log q)^4)`. For a harmonic at the orbital frequency, the
product-throat transfer correction to the same instantaneous static
response is relatively `O((Omega ell)^2)=O(q^-3 (log q)^6)`.

The last statement uses the prescribed end coefficients in (17) as driven
data. Our RN gauge calculation is static. Finite-frequency corrections
from the actual moving RN exterior, its companion tidal field and
radiative matching have not been bounded. In the formal overlap
`|rho|>>|sigma|`, the local RN frequency parameter Omega P is small,
but that observation is not a proof of the missing dynamic matching.

Nor does (18) bound boundary coordinate shifts, monopole-frame rotations,
charged-spinor phases, other l=1 modes, l>=2 tidal perturbations, or
second-order mean stress. It is a small classical response in one
specified physical scalar branch. The core radius and quantum source
response still require their own coupled calculation.

## 5. Verification and contribution decision

Run `checks/dipole_throat_matching.py --output evidence/dipole-throat-matching.json`.
The script constructs all background and linear Einstein-Maxwell tensors
for (8), verifies the gauge transformation and the regular coefficient,
and checks the Darboux solutions, gap factorization, central transfer,
scales and Burko's extremal potential limit. Independent
quadratures and a stiff initial-value solve check (4). Refined boundary-value
solves on three overlap domains compare the finite-end solutions with
(13) and the predicted inverse-square core error. These are supplementary
component numerics, not a four-dimensional wormhole simulation.

The [source review](../evidence/dipole-throat-source-review.json) distinguishes
the known mass-six reconstruction from the present matching calculation.
Sadhukhan's [polar-sector paper](https://arxiv.org/html/2609.33511v1), Section
2.2, explicitly excludes a gauge-fixed l=1 calculation; Porfyriadis's
[near-extreme matching paper](https://doi.org/10.1140/epjc/s10052-019-7347-6),
Section 2, likewise treats l>=2. Those exclusions do not establish novelty;
Burko already treats l=1 and reconstructs its metric and Maxwell fields.
Its complete main derivation has now been read, including the special
dipole constraint and reconstruction. Only the potential limit above has
been independently checked from its formulas; its cited earlier proofs
have not all been reviewed. The classical dipole physics is prior work.
The matched binary and quantum response still need derivation and a
completed contribution comparison before accepting a new physical result.

Continue physical long-return damping with matched finite-duration band
transfer as fallback. The leading classical scalar boundary coefficient
and its regular product-throat response are now explicit. Retain the
coordinate/frame information from the mouth and complete the remaining
binary, boundary and quantum equations before interpreting this as a
physical transport or stability result. No new consequential endpoint,
full-source review, or Lean certificate is claimed.

Remaining to-do list: derive finite-frequency binary/mouth matching and
the remaining boundary-frame, vector and tidal data; include semiclassical
metric/current response and the actual time-dependent quantum state,
exterior determinant and discarded tail; finish the literature/proof/novelty
review and consequential endpoint; independently verify it, certify suitable
mathematics exclusively on GCP, and prepare the reproducible manuscript.
