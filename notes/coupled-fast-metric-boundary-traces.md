# Fast curvature traces in the reconstructed metric

2026-10-21, E003/S168. Ordinary argument before focused verification.
The full stable-wormhole goal remains incomplete.

The fast curvature coefficients remain visible in the conformal metric
amplitude after the full auxiliary thermal shape responds. The shape
response cannot cancel their leading term. Both coefficients nevertheless
leave the leading conformal boundary metric and its clock unchanged.
Thus those prescribed leading sources do not determine the extra
gravitational boundary condition.

For the strong auxiliary thermal system, retaining zero, one or two
independent fast trace amplitudes gives Fredholm index -1, 0 or +1,
respectively. The last realization has a nonzero homogeneous kernel at
every admitted thermal frequency. These are statements about specified
auxiliary domains, not constructed physical modes or a choice of a
gravitational boundary action.

## 1. What the retained preparation actually specifies

The [same-source charge calculation](stationary-same-source-boundary-charge.md)
fixes the leading conformal metric, clock, magnetic bundle and paired
Dirac prescription. It explicitly leaves the full gravitational boundary
completion as a matching obligation. The later
[stationary variational phase](stationary-fast-phase-variation.md) proves
differentiability on a chosen stationary class that excludes oscillatory
tails; it explicitly does not select that class as the physical dynamical
phase. The [local anomaly boundary calculation](spherical-anomaly-boundary.md)
also shows why leading metric/radius data alone do not determine all
higher-derivative boundary data.

None of those records supplies a time-dependent fast-trace phase. We do
not infer one from the matter reflector, the stationary weight, or the
availability of an auxiliary Friedrichs shape inverse.

## 2. An exact reconstruction identity on the maximal graph domain

Use the same charged seed, action, fixed state/domain and physical
H=L2(dy) as the [coupled boundary index](coupled-thermal-boundary-index.md).
Fix a positive Matsubara frequency above the retained thresholds. Put

    E=60pi/G+11lambda/3, d0=E+4lambda>0,
    Z=R+12lambda, H_R=2R+rR_r,
    K_gamma=-partial_y^2+gamma^2+F(F_r/r-H_R/6),
    Atilde=-partial_y^2+F(F_r/r-E+11Z/36-c_R).

Let D_max be the maximal H graph domain of Atilde, with its limit-point
horizon endpoint. Its outer trace pair is (A_c,A_s), with

    U(y)=sqrt(y)[A_c cos(k log(y/y0))+A_s sin(k log(y/y0))]
                                                   +o(y^(3/2)),
    k^2=E/lambda-9/4>0.                                  (1)

The corresponding differentiated remainder is o(y^(1/2)). These are
the retained scalar trace results for H forcing, including the bounded
seed remainder. They apply to the full reduced equation because its
feedback is bounded H->H.

Define the exact bounded coefficient

    h=c_R-(23Z+6rZ_r)/36.

Then

    K_gamma-(Atilde+gamma^2)=F(d0+h).                     (2)

At the outer end h=O(r^-3), so Fh=O(y); at the horizon Fh tends to
zero. In particular Fh is a bounded multiplier on the whole exterior.
For U in D_max, the physical reconstruction q_U=R_gamma U therefore
obeys the exact identity

    q_U=U/(6d0)+z_U,
    z_U=-(1/(6d0))K_gamma^-1[(Atilde+gamma^2)U+FhU].    (3)

The inverse on the right is the ordinary unweighted positive resolvent
on H. Its source belongs to H, so z_U belongs to D(K_gamma). To prove
(3), apply K_gamma in distributions and use K_gamma R_gamma U=FU/6.
Both sides are in H, and uniqueness of the limit-point K_gamma
realization identifies them. No application of the unbounded multiplier
F to an arbitrary H input is used to justify a resolvent norm.

The gain in (3) is useful at the boundary. Every z in D(K_gamma) has

    z=o(y^(3/2)), z_y=o(y^(1/2)).                         (4)

Here is a direct endpoint proof. Locally write
(-partial_y^2+6/y^2)z=g with g in H; the bounded potential remainder
and gamma^2 z are included in g. Excluding the non-L2 y^-2 solution,
variation of constants on (0,b) gives

    z=c y^3+(1/5)[y^-2 integral_0^y s^3 g(s)ds
                           +y^3 integral_y^b s^-2 g(s)ds]. (5)

The first term in brackets divided by y^(3/2) is bounded by a constant
times ||g||_(0,y), which tends to zero. Split the second integral at
eta: the part from y to eta has arbitrarily small coefficient when eta
is small, and the part from eta to b is O(y^3). This proves the first
claim; differentiating (5), with cancellation of its two g(y) terms,
proves the second. The proof needs no pointwise bound on g.

## 3. The full shape cannot erase the metric trace

For any w in the full closed shape energy space E_gamma, the retained
coercivity and reconstruction estimates give

    w in H1_0, Q_gamma w in H1_0,
    ||(Q_gamma w)_y||^2 <=10N_gamma[w]/9,
    N_gamma[w]<=33G a_gamma[w].                         (6)

The zero outer traces follow by continuous H1 convergence from the
compact shape core and its reconstructed tails. For any f in H1_0,

    |f(y)| <= sqrt(y)||f_y||_(0,y)=o(sqrt(y)).            (7)

This is a pointwise statement for the absolutely continuous
representative. Apply it to w and Q_gamma w. With r^-1=lambda y+O(y^3),
q=R_gamma U+Q_gamma w, p=q/r and v=w/r, equations (1)-(7) imply

    p(y)=[lambda/(6d0)] y^(3/2)
          [A_c cos(k log(y/y0))+A_s sin(k log(y/y0))]
                                                +o(y^(3/2)),
    v(y)=o(y^(3/2)).                                    (8)

The result applies in particular to w=w_ell-K_cross^dagger U after
solving the complete auxiliary thermal shape row, for every continuous
shape source ell. Its leading coefficient is independent of that
source. Equation (8) is not obtained by setting w to zero or cutting off
its conformal tail. It is stronger than the pure-vacuum conformal lift:
the seed is charged, its residual contacts remain, and shape is included.

In the retained optical chart use the fixed background defining factor
Omega=B_ref/r_seed. The angular relative variation of the compactified
metric is 2p; the orbit relative variation is 2(p+v). Their leading
oscillatory coefficients are therefore both

    [lambda/(3d0)](A_c,A_s).                             (9)

Every such coefficient pair leaves the leading boundary metric and
clock unchanged, since p,v tend to zero. Conversely, neither pair can
be hidden by the energy-class shape response in this chart.

If p has a C2 extension to this fixed compactification, then p(0)=0
and p=o(y) force p_y(0)=0, so p=O(y^2). Equation (8) then forces both
A_c and A_s to vanish: any nonzero sinusoid has nonvanishing values
along a sequence tending to zero. This proves necessity only. Zero fast
traces do not by themselves prove C2 regularity of every energy-class
shape response. Nor is C2 regularity in this chart imposed as a physical
axiom, or asserted invariant under nonsmooth changes of compactification.

## 4. How many trace conditions the auxiliary equation needs

For a fixed complex linear subspace L of the two-dimensional trace space,
let D_L be the subspace of D_max whose trace lies in L. Let d=dim L.
The trace map is continuous and onto; cutoffs of the two normalized
local homogeneous scalar solutions realize arbitrary pairs. Thus

    D_L=D_min plus d independent trace lifts.            (10)

The direct sum is topological in the graph norm. The bounded full feedback
B_gamma from the preceding result acts on all of H, so L_gamma extends
to each D_L. The index changes by the number of added domain coordinates:

| Admitted outer amplitudes | Domain choice | Index of full auxiliary pair |
| --- | --- | ---: |
| None | A_c=A_s=0 | -1 |
| One independent amplitude | One homogeneous trace relation | 0 |
| Both independently | No fast trace relation | +1 |

For completeness, view L_gamma on D_min as a bounded Fredholm graph-space
map. Extending its domain by C^d first adds d zero columns, raising its
index by d. Replacing those columns by the actual images of the trace
lifts is a finite-rank perturbation, so the index is unchanged. The
positive shape block remains an isomorphism after triangular elimination.
This proves the table for the complete auxiliary system.

Consequently the maximal realization has dim kernel=dim cokernel+1>=1
at every admitted thermal frequency. These homogeneous strong row
solutions reconstruct metrics satisfying the leading boundary-source
condition (9). They are not asserted to meet the remaining tensor
constraints, the physical horizon domain, the full boundary variational
problem or an exact nonlinear background. Their existence identifies
underdetermination of this auxiliary prescription, not physical instability.

A scalar no-flux phase cos(theta)A_c+sin(theta)A_s=0 is one possible
one-relation domain for the diagonal scalar Green form. Index zero does
not imply invertibility, positivity, a full gravitational no-flux law,
or stability. The full feedback need not be symmetric in physical H.
No particular theta, boundary action or source ensemble is chosen here.

## 5. Sources, verification and next decision

[Ishibashi--Wald, section 1 and section 3](https://arxiv.org/pdf/hep-th/0402184v4)
explain how dynamical boundary domains supplement the bulk equation in
AdS, and why a positive self-adjoint scalar realization is an additional
requirement. Their AdS scalar classification does not supply our charged
quantum domain. The finite-extension index argument uses the standard
compact-perturbation result in
[Haller, Corollary I.3.3](https://www.mat.univie.ac.at/~stefan/files/ASIT/ASIT.pdf).
The charged reconstruction identity and energy-space trace transfer are
derived here; generic boundary-extension theory is not claimed as new.

The checker verifies the exact charged operator difference, endpoint
coefficient cancellations, Green representation, trace normalization and
domain-index arithmetic. High-precision Green quadrature in the explicit
inverse-square reconstruction model tests (3), its leading coefficient
and a nonzero shape tail at two precisions. Exact finite-rank infinite
shift models test all three domain dimensions and the distinction between
index zero and invertibility. These fixtures are not an evaluation of
the actual quantum shape kernel or a charged spectrum. The continuum
boundary theorem remains an ordinary proof; no Lean campaign is selected.

The physical task is now more specific: match one dynamical gravitational
trace rule, or prove compatibility of the stronger stationary exclusion
with the actual source and initial-data class. Leading metric/clock data
alone do not supply it. The original matter domain, action, counterterms,
state and all physical completion gates remain unchanged and incomplete.

Remaining to-do list: matched gravitational boundary rule and source
compatibility; continuous causal response and constrained stability;
nonlinear existence and signed compatibility; exterior/state,
charge/return, opening, all-sector stability, transport, preparation,
separate creation and full audit.
