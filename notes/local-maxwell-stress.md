# Local Maxwell stress and reference throat focusing

2026-10-11, S053. The source-free Maxwell completion has a controlled local
stress in a regular spherical neck/core satisfying an explicit geometric
inequality. The field correction is O(q/d^2), independently of the auxiliary
port radius. Its radial null stress is positive and quadratic in the
transverse correction; the monopole cross term cancels exactly. Its direct
contribution to the leading spherical focusing equation is subleading to
the specified LLL Casimir opening in the selected scale family.

This is a test-field source and reference-focusing result. It does not solve
the nonspherical Einstein-Maxwell system, derive the quantum state, or
control the complete second-order monopole response. In particular, angular
metric corrections and boundary data cannot be discarded using a bulk
Kaluza-Klein mass gap alone.

## 1. A comparison function tied to the Einstein equations

Keep the source-free field, gauges and prescribed spherical throat of
[coupled matching](coupled-maxwell-matching.md). Initially d=1. Write

    ds^2=N(x)^2(-dt^2+dx^2)+R(x)^2 dOmega^2,  w=N/R,
    a_T=sin(theta)^2 b(x,theta),
    b_xx+w^2(Delta_(S^4)b-2b)=0.                         (1)

The sphere in (1) is the mathematical axis-regular lift used in the
[Dirac comparison](maxwell-dirac-comparison.md). Its maximum-principle
argument, without an optical cap or attenuation hypothesis, already gives

    |b|<=C q p^2                                         (2)

everywhere on the throat. The constant is uniform in its middle profiles.

Take an inner segment J whose two endpoints have areal radii at least
c p, for a fixed c>0. It contains the neck and core to be estimated.
On J require

    (R^2)''-2N^2<=0.                                    (3)

No such condition is imposed on an arbitrary outer connector to the
artificial Euclidean collars. This distinction matters: those collars
are not an exact Einstein-Maxwell description of a massive mouth. The
physical exterior and all junctions still need their own metric solution.

The positive function R^2/(c p)^2 is a supersolution of (1) and is at
least one on the two endpoints of J. Apply the maximum principle to
both signs of b using (2). It follows that

    |b(x,theta)|<=C q R(x)^2,          x in J.            (4)

The auxiliary cut p has disappeared. This is an all-angular-mode
comparison; no derivative bound on the port trace is assumed.

Condition (3) has a direct Einstein interpretation. Put phi=log N and
define rho=T_hat(t)hat(t), p_r=T_hat(x)hat(x) for the effective source of
the spherical reference, G=8pi G_N T. Direct tensor contraction gives

    8pi G_N rho = 1/R^2
       + [2phi' R'-2R''-(R')^2/R]/(N^2 R),
    8pi G_N p_r = -1/R^2
       + [2phi' R'+(R')^2/R]/(N^2 R),
    (R^2)''-2N^2 = 8pi G_N N^2 R^2 (p_r-rho).           (5)

Thus rho>=p_r suffices. This is not the null energy condition, which
would instead constrain rho+p_r. In particular a negative traceless
two-dimensional Casimir source has rho=p_r and does not itself violate
this comparison inequality. Nor does writing G=8pi G_N T prove that T
is the physical quantum source; for a prescribed metric it is just a
way to state the geometric condition. Formula (5) remains valid at an
areal minimum and does not use R as a coordinate.

On an extremal RN segment with areal coordinate r, N=1-r_e/r and
dx/dr=N^-2, direct substitution gives

    (R^2)''-2N^2=-2N^2 r_e^2/r^2<0.                    (6)

On a constant-radius AdS throat it is -2N^2. For the first-order
spherical reference `R^2=r_e^2 C(rho)`, `N^2=r_e^2 S(rho)/ell^2`,
`dx/drho=ell/S`, (3) is equivalently `S'C'+S C''<=2`.
The previous Casimir reference has, on its controlled overlap,
`C=1+alpha c`, `S=1+rho^2+alpha s`, with the left side O(alpha(1+|rho|))
and smaller higher-order terms when the specified corrections are small.
This verifies compatibility on the analytic pieces and their valid
overlap; it is not a proof for an arbitrary physical interpolation.

## 2. From potential to local magnetic field

The additional local regularity hypothesis is explicit. At each point
of the region being estimated, let R_0 be its areal radius. In proper
distance ds=N dx, require an interval of radius c_0 R_0 contained in J
on which R/R_0 and N/N_0 are bounded above and below by fixed positive
constants, with uniformly bounded scaled C^(2,alpha) norms, 0<alpha<1.
Constants are independent of q,p and the central redshift N_0. The
estimate is for this interior region; a physical junction is not included
without checking the same local hypotheses there.

In proper distance (1) becomes

    b_ss+(N_s/N)b_s+R^-2(Delta_(S^4)b-2b)=0.              (7)

Rescale the interval by R_0. The operator is uniformly elliptic with
controlled coefficients. An interior gradient estimate, applied to (4), gives

    |b_s|<=C q R_0,       |grad_(S^4)b|<=C q R_0^2.       (8)

For completeness this estimate follows from a local Bernstein argument.
In sigma=s/R_0, put r=R/R_0, n=N/N_0 and
`h=d sigma^2+r^2 dOmega_4^2`. Equation (7), multiplied by R_0^2, is
`D b=V b`, where `D=Delta_h+Z.grad_h`,
`Z=(n'/n-4r'/r) partial_sigma`, and `V=2/r^2>=0`.
The hypotheses bound the Ricci tensor of h, Z, its covariant derivative,
V and grad(V), with constants independent of N_0. Normalize b to a
function u with |u|<=1 using (4) and the local upper bound on r.
For W=|grad u|^2 the Bochner identity gives

    D W >= 2|Hess u|^2-C W-C,
    D(u^2)=2W+2V u^2>=2W.

Choose a smooth cutoff eta(sigma) supported away from the two interval
endpoints and equal to one near its center, with bounded first two
derivatives. For `F=eta^2 W+A u^2`, the product cross term is bounded
below by `-eta^2 |Hess u|^2-C W`. Therefore
`D F>=(2A-C)W-C`. Take a fixed sufficiently large A. At any interior
maximum of F, D F<=0, hence W<=C and F<=C; on the interval boundary
F=A u^2<=A. The compact sphere has no boundary. The maximum principle
thus gives W<=C where eta=1, proving (8) after rescaling.
This uses the coefficient bounds, rather than an energy-to-pointwise
inference. There is no angular-cutoff dependence or inverse power of N_0.
The exact RN and AdS pieces satisfy these proper-distance bounds on
fixed fractions of their areal scale. For example RN has
`R_s=N`, `R(N_s/N)=r_e/r<=1`; in the AdS core
`N=N_min cosh(s/r_e)` and R=r_e. Their higher scaled derivatives are
bounded as well. Smoothness without quantitative bounds would not suffice.

With B_0=q/(2R^2) in the radial direction, the corrected magnetic field is

    B_hat(r)=B_0+beta_r,
    beta_r=[2cos(theta)b+sin(theta)b_theta]/R^2,
    beta_theta=-sin(theta)b_s/R,       beta_phi=0.         (9)

Regularity on the lifted sphere makes (9) regular at both poles.
Equations (4),(8) prove, now restoring separation d,

    |beta|<=C q/d^2,      |beta|/|B_0|<=C(R/d)^2.        (10)

As a separate exact check, an applied uniform test field on extremal RN
has the regular l=1 potential
`a=(B_ext/2)(r^2-r_e^2) sin(theta)^2`. The source-free radial equation
is satisfied identically. Its components are
`beta_r=B_ext(1-r_e^2/r^2)cos(theta)` and
`beta_theta=-B_ext(1-r_e/r)sin(theta)`, so |beta|<=|B_ext|.
This known type of test-field behavior does not include its gravitational
response or determine the actual coupled wormhole boundary coefficients.

## 3. Stress, angular average and the null cancellation

For a static purely magnetic field in an orthonormal frame,

    T_tt=|B|^2/(2g^2),
    T_ij=(delta_ij |B|^2/2-B_i B_j)/g^2.                 (11)

Subtract the monopole stress at the same metric. Every component obeys

    |Delta T_ab|<=C q^2/g^2 [1/(R^2 d^2)+1/d^4].        (12)

For the specified unit-normalized radial null vector k=e_hat(t)+e_hat(r),

    Delta T_kk=T_kk=|B_transverse|^2/g^2
       =beta_theta^2/g^2>=0,
    Delta T_kk<=C q^2/(g^2 d^4).                         (13)

The reference monopole has T_kk=0. There is no term linear in beta in
(13), even pointwise. Other components do have linear terms, including
off-diagonal radial/angular stress, so (13) must not be used for them.

The spherical mean of each diagonal linear cross term also vanishes:
the zero-flux perturbation satisfies
`int beta_r dOmega = (2pi/R^2)[a(pi)-a(0)]=0`.
Equivalently its angular derivative integrates to zero. The first
nonzero mean Maxwell correction is quadratic, O(q^2/(g^2 d^4)).
This is an average on the fixed spherical reference, not the monopole
part of a fully backreacted nonspherical metric.

In a bounded central core where R is comparable to r_e, the specified
leading LLL Casimir source has magnitude O(q/r_e^4). Therefore

    |Delta T_ab|/(q/r_e^4)<=C (q/g^2)(r_e/d)^2,
    Delta T_kk/(q/r_e^4)<=C (q/g^2)(r_e/d)^4,             (14)

where r_e/d is small and its fourth-power term is absorbed in the first
bound. This is a comparison to that stated source, not a proof of the
source in the actual interacting exterior-connected state.

## 4. Signed contribution to the leading focusing equation

There is also an integrated consequence that keeps the null weighting.
On the leading AdS core use `f=1+rho^2` and
`ds^2=r_e^2[-f d tau^2+d rho^2/f+dOmega^2]`, with external time
`t=ell tau`. The exact spherical null
Einstein identity in optical coordinates is

    R''-2(N'/N)R'=-(8pi G_N/2) N^2 R T_kk.              (15)

Linearizing R about r_e and using dx/drho=ell/f gives

    partial_rho^2 Delta R=-(8pi G_N r_e^3/(2f)) Delta T_kk.
                                                               (16)

Here Delta T_kk may be replaced by its spherical mean to compute the
direct linear spherical response to this specified Maxwell source.
The trace of the invariant product reference does not enter (15).
For the established leading LLL source with fixed chi>0,

    T_kk,Q=-chi q/(16pi^2 r_e^4 f),
    partial_rho^2 R_Q=alpha chi r_e/f^2,
    alpha=8pi G_N q/(32pi^2 r_e^2)=g^2/(4pi^2 q).        (17)

Choose a reflection point at rho=0 and compare the changes of slope up
to any rho_b>0 in the controlled reference/overlap region. The direct
Maxwell term decreases that slope. Equations (13),(16)-(17) give

    0<=-Delta R'_M(rho_b)/R'_Q(rho_b)
       <=(32pi^2 C/chi)(q/g^2)(r_e/d)^4,                (18)

since

    I_1=int_0^rho_b d rho/f=atan(rho_b),
    I_2=int_0^rho_b d rho/f^2
       =[atan(rho_b)+rho_b/(1+rho_b^2)]/2,
    1<=I_1/I_2<=2.                                      (19)

The comparison is uniform in rho_b as long as (10) and the leading
reference approximation hold. It does not extrapolate a first-order
metric to infinity at fixed alpha. In an overlap require alpha rho_b<<1
and all other retained reference corrections small.

With `g^2=q^-2`, `delta=d/ell=(log q)^-2` and the leading MMP scale
`ell/r_e=16pi q^3`, the component ratio in (14) is
`O(q^-3 (log q)^4)`. The null ratio and direct focusing ratio are
`O(q^-9 (log q)^8)`. Thus this Maxwell source cannot by itself overwhelm
the prescribed leading Casimir opening in that limit. These estimates
do not solve the full support or stability problem.

## 5. Why this is not yet the Einstein response

The nonspherical linear Maxwell stress in (12) sources nonspherical
metric fields. Those fields change Maxwell's equation and its stress.
Their quadratic contribution to the spherical equations, and the
quantum response to them, are not included in the direct source term
(18). A small source alone does not prove a uniform inverse estimate
for a long-domain Einstein constraint problem.

The relevant bulk comparison is established in Michelson and Spradlin,
[hep-th/9906056v3](https://arxiv.org/abs/hep-th/9906056v3), Section 3.
For l>=2 its scalar mass matrix, with L=l(l+1), is

    M=[[L+2,-4L],[-1,L]].                                (20)

The similarity diag(1,2sqrt(L)) makes it symmetric, with eigenvalues
`l(l-1)` and `(l+1)(l+2)`, both at least two. The physical l=1 scalar
has mass squared six; an extra scalar boundary mode is not fixed by
that bulk field equation. The source also retains massless boundary
vectors and the l=0 boundary graviton. A homogeneous-boundary static
scalar inverse on the exact product can use its positive mass, but
the actual mouth data and constraints are not homogeneous boundaries.
We therefore do not use (20) to claim the full Einstein inverse exists
or is uniformly bounded. The displayed scalar/vector mass matrices and
the l=1 scalar elimination are checked algebraically in the supplement;
this is not an independent reconstruction of all four-dimensional mode
equations. The sourced global reconstruction remains open.

The LLL stress and leading focusing/matching comparator are the existing
[MMP v3](https://arxiv.org/html/1807.04726v3), Sections 5.2-5.3.2, and
[the reproduced reference](casimir-reference-matching.md). Their magnetic
mouths also require orbital/other support; no static global equilibrium
is supplied here. Neither the underlying Casimir construction nor
magnetic attenuation is accepted as a new consequential wormhole result.

## 6. Verification and decision

Run `checks/local_maxwell_stress.py --output
evidence/local-maxwell-stress.json`. The script constructs the Einstein
and Maxwell tensors directly, verifies the two radial identities, null
and angular cancellations, the RN test field and bulk mass matrices.
Independent numerical radial solutions, tensor contractions and focusing
quadratures test the relevant signs and scales, with refinement. These
are supplementary checks, not a solved global metric or quantum vacuum.
The source supplement records reviewed derivations and inaccessible
material without treating the latter as reviewed.

Continue physical long-return damping and matched finite-duration band
transfer. The local classical Maxwell source is now controlled under
explicit inner-geometry hypotheses. The next decisive issue is the
coupled metric and boundary response, followed by the state-dependent
renormalized quantum response; that requirement is not discharged by
the positive bulk mass spectrum or the direct focusing diagnostic.

Remaining to-do list: solve the coupled nonspherical metric/Maxwell
response with physical mouth data, verify the inner geometric and local
regularity conditions, derive the quantum state and full renormalized
response including exterior determinant and discarded tail; finish the
literature/proof/novelty review and consequential endpoint; independently
verify it, certify suitable mathematics exclusively on GCP, and prepare
the reproducible manuscript.
