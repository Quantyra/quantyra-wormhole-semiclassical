# A quantitative canonical matter domain on the physical conformal ball

2026-10-22, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: progress, technical 6b1f0a0 and planning c379e9c.
The full stable semiclassical wormhole objective remains incomplete.

The actual canonical fermion realization is stable on a quantitative
open real-metric neighborhood of the product-centered seed. For

    N1(b)=max(||b||infinity, ||grad_gamma0 b||infinity/100)<=rho,
    rho=1e-8,

the massive first-order relative defect is <3.78e-8, uniformly over the
entire angular tower. The massless sector remains invertible by an exact
conformal inverse formula, without needing its unknown absolute spectral
gap to establish this neighborhood. The graph norms are equivalent with
factor <1.000000511. These statements quantify a matter-domain input of
the physical solve; they do NOT prove the renormalized metric Hessian
inverse, its nonlinear map, or a complete gravitational boundary domain.

The preceding feasibility audits also overlooked a stronger already
available outer admission theorem. Section 6 records the correction
without changing the frozen historical notes or importing its constants
into the new construction without a transfer proof.

## 1. Physical family, fixed units and existing conformal structure

Keep the S169 c_R=10^22 theory, full vectorlike Dirac field, magnetic
bundle, fixed q, mu0, positive temperature and original paired projector.
Let B and A=B/sqrt(C) be the current exact product radii. On the same
fixed smooth disk and thermal angle use

    r=B+(r0-B)exp(V),
    g4=(r/B)^2 [A^2 gamma_b+B^2 dOmega^2],
    gamma_b=exp(2b)gamma0.                            (1)

The canonical orbit gamma0 is the ACTUAL product-centered seed, not an
RN reference to which an unquantified transfer is silently applied.
All derivatives and operator masses below use its original unit A.
The norm N1 is the value/first-jet component of the retained j!100^j
units. The complete physical solution norm X must dominate this component
along with cold weights, interior raw C5/cap, outer weighted jets and
mass/return data. This note does not show that an unknown correction has
small X or N1. In particular, N1 cannot replace the missing higher norm.

Allow real radial b smooth on the open disk, regular and even at the
bolt, and C1 through its canonical outer boundary. The physical fast
tail b=m Y^3+O_(5)(Y^4 log^p Y) has this first-order regularity even
though its fifth ordinary derivative need not be bounded at Y=0. Fix
the clock and leading boundary data; the tail vanishes there. Higher
coefficient requirements remain those of the complete admission space.
Nonradial/time-dependent/electric perturbations and changing paired angles
are not included. Both copies have identical fields and fixed bundle map.

The existing exact common-subtraction identity is

    Gamma_D[g4]=Gamma_hat[gamma_b]+W4[log(r/B);gamma_b]. (2)

The canonical full tower has masses M_n=sqrt(n(n+q)/C), M_n>18sqrt(n),
and multiplicities d_n=2(q+2n), with q massless orbit fields. The fixed
sphere, scale, thermal angle and first-order boundary projectors imply
that Gamma_hat is independent of V in (1). Thus its V and mixed b,V
derivatives vanish in these exact physical coordinates. This was already
proved in [the canonical audit](recentered-canonical-tt-audit.md), and
is not a newly discovered cancellation. The local Weyl functional,
boundary cocycle, common restoration and off-shell chart/density contacts
remain. The four-dimensional radius source is NOT set to zero by (2).

At a finite physical cut the Weyl transformation and its projector are
exact; at the outer limit we use the specified canonical realization of
[the conformal-domain proof](stationary-exterior-conformal-domain.md).
The unbounded absolute factor r/B at infinity is not used as a bounded
similarity on physical L2. The present comparison is between canonical
metrics related by the BOUNDED relative factor exp(b). Physical action
limits still require the retained boundary transgression and fast phase.

## 2. Actual first-order domains and their graph comparison

The orbit operator D has the complementary adjoint domains

    Dom D: Pi_- psi=0, Dom Ddagger: Pi_+ psi=0,
    Pi_+=(1+chi)/2, chi=tau3 gamma_n.                  (3)

There is no boundary at the smooth bolt. In a fixed spinor L2 space,
two-dimensional conformal covariance with the unitary half-density
identification gives

    D_b=exp(-b/2)D_0 exp(-b/2).                       (4)

The bundle conformal formula is [Walpuski, Proposition 15.5,
printed page 79](https://walpu.ski/Teaching/SpinGeometry.pdf). Multiplying
it by the two-dimensional L2 density weight exp(b) gives (4).
That proposition alone is not a boundary-domain theorem. Here scalar
multipliers exp(+/-b/2) commute with Pi_+ and Pi_-, preserve H1, and
have bounded inverses. The trace condition in (3) is therefore the same
closed H1 domain as a SET, for both adjoints, throughout this C1 ball.
The usual local ellipticity can also be checked directly: the allowed
gamma_n eigenspace is transverse to the decaying tangential-symbol
subspace since gamma_n anticommutes with i gamma_n gamma_t.

Use the retained parallel involution J which exchanges these domains
and commutes with scalar multipliers. The paired identity
Ddagger=-J D J implies that

    A0=-J D0, A0*=-A0,
    A_b=exp(-b/2) A0 exp(-b/2), A_b*=-A_b             (5)

as CLOSED operators on the allowed domain. A bounded invertible real
multiplier M gives (M A0 M)*=M A0* M with precisely the transported
domain; the H1 multiplication argument identifies it with (3).
The Clifford symbol c0(db)=[A0,b] is Hermitian and has norm |db|.
Expanding (5) gives the useful exact identities

    A_b=exp(-b)[A0-c0(db)/2],
    A0=exp(b)A_b+c0(db)/2.                            (6)

For epsilon=||b||infinity, ell1=||db||infinity and graph norms
G_b(psi)=||psi||2+||A_b psi||2, these imply in BOTH directions

    G_b<=exp(epsilon)(1+ell1/2)G_0,
    G_0<=exp(epsilon)(1+ell1/2)G_b.                   (7)

At epsilon<=1e-8 and ell1<=1e-6 the factor is <1.000000511.
This is a comparison in fixed unscaled Hilbert and graph norms, not
a choice to declare a large unknown inverse equal to one. Only the
first spatial derivative of b is used. C1 approximation preserving the
radial bolt and outer trace passes (4)-(7) to weighted logarithmic tails;
no closed-boundary C5 or global spatial analyticity is being asserted.

The induced square domains are still obtained by composition: Pi_- psi=0
and Pi_+(nabla_n+K/2)psi=0, with complementary conditions for the partner.
Their Robin coefficients vary with the metric. They are not frozen, and
no artificial angular-mass Robin term is introduced. Equations (4)-(7)
concern first-order graphs and do not identify the two squares by a
scalar similarity or discard their anomaly.

## 3. Uniform transfer for every massive angular sector

For M=M_n>18, let R_b,M=(A_b+M)^-1. Skew-adjointness gives directly

    ||R_b,M||<=1/M,
    ||A_b R_b,M||<=1                                (8)

for every real b in this family. The second inequality follows from
|it/(M+it)|<=1 by the spectral theorem. It is sharper than the triangle
bound two and requires the actual skew-adjoint realization, not merely
an anti-Hermitian differential expression.

Equation (6), composed on the RIGHT with R_0,M, gives

    ||(A_b-A0)R_0,M||
       <=exp(epsilon)-1+exp(epsilon)ell1/(2M)
       <=eta(epsilon,ell1)
       :=exp(epsilon)-1+exp(epsilon)ell1/36.           (9)

At the target ball eta(1e-8,1e-6)<3.78e-8, far below .25. This is a
numeric graph-relative transfer for the canonical MATTER inverse.
It is not the .25 target for the gravitational/quantum Hessian; these
are different operators on different spaces. The resolvent identity
and (8) also give

    M||R_b,M-R_0,M||<3.78e-8.                         (10)

Multiplicity does not alter a direct-sum operator norm, so (9)-(10)
hold uniformly on the massive direct sum. They do not bound the trace
of that sum. Coincident kernels, ultraviolet subtractions and derivatives
of their finite parts still need independent angular/spatial majorants.

For completeness the exact multiplier formula also defines a complex
C1 parameter family on the same H1/projector domain. If N1(b)<=.01,
then (9) with epsilon=.01,ell1=1 bounds its graph perturbation by <.0382.
The complexified fields remain radial, so their Clifford gradient has
one complex scalar coefficient times a real unit Clifford direction;
no unaccounted conversion from a general complex vector norm is used.
Neumann inversion supplies ||R_b,M||<1.04/M on that complex tube.
It is a first-order operator tube only, NOT an enlargement of an actual
renormalized source theorem or permission to disregard its C5 contracts.

Useful real first and second directional bounds follow without any new
domain variation. At an admitted b, for U=N1(u),Vd=N1(v),

    A_b'[u]=-u A_b-c_b(du)/2,
    A_b''[u,v]=uv A_b+[u c_b(dv)+v c_b(du)]/2,
    c_b=exp(-b)c0.

Therefore ||A_b'[u]R_b,M||<=a1 U and
||A_b''[u,v]R_b,M||<=a2 U Vd, where

    a1=1+50 exp(1e-8)/18 <3.78,
    a2=1+100 exp(1e-8)/18 <6.56.

Differentiating the actual inverse on its fixed graph gives

    M||R_b,M'[u]||<3.78 U,
    M||R_b,M''[u,v]||<35.14 U Vd,                    (11)

since 2(3.78)^2+6.56=35.1368. These are operator response constants,
not diagonal stress or complete metric-map Lipschitz constants Z2.
In particular summing (11) with d_n alone is divergent and is invalid.

## 4. The LLL: exact inversion instead of a small-gap perturbation test

At the fixed positive temperature the actual canonical seed is a compact
smooth disk with a finite positive conformal factor, including at its
outer boundary. The retained kernel-free argument is short enough to
check directly. Map a zero mode to the flat disk by conformal covariance.
In a Cartesian spin frame its components are f(z) holomorphic and
g(bar z) antiholomorphic, regular at zero. On the circle the allowed
trace obeys f=+/-exp(-i phi)g. The two sides have Fourier powers n>=0
and -(n+1)<=-1 respectively, so both vanish. The complementary adjoint
uses the opposite sign and has the same conclusion. Ellipticity on the
compact disk and trivial kernel/cokernel give a bounded inverse R0=A0^-1.
Its absolute norm has not been numerically evaluated and is not set one.

The exact identity (5), with the actual domain, gives

    R_b=exp(b/2) R0 exp(b/2),
    ||R_b||<=exp(epsilon)||R0||,
    ||R_b-R0||<=(exp(epsilon)-1)||R0||.               (12)

Thus at epsilon=1e-8 the relative change is <1.00000001e-8. If lambda0
is the positive lowest singular-square eigenvalue of the reference,
then lambda_b>=exp(-2epsilon)lambda0. Kernel absence is preserved on
the entire ball, with a quantitative RELATIVE gap. No temperature-zero
limit or uniform absolute gap as h->0 is claimed.

The previous audit's ||delta D R0||<=1/2 was a sufficient generic test,
not a necessary condition for this exact conformal family. Its unknown
gap-dependent norm need not be computed to establish (12). For C1
directions, endpoint differentiation gives exactly

    R_b'[u]=(u R_b+R_b u)/2,
    R_b''[u,v]=(uv R_b+u R_b v+v R_b u+R_b uv)/4.      (13)

Consequently the first/second norms are <=exp(epsilon)||R0|| times
||u||infinity and ||u||infinity||v||infinity respectively. There is no
power of the unknown inverse gap from repeated generic insertions.
Both partner domains obey the same formula. The q angular zero modes
only supply their finite degeneracy. The exact HH stress/Polyakov
formula remains the appropriate way to evaluate this sector's canonical
stress; (12)-(13) alone do not bound the complete four-dimensional source.

## 5. What requirement is now closed, and what is still decisive

The stated real canonical conformal ball has a common paired first-order
domain, controlled graph equivalence, uniform massive inverses and
relative LLL invertibility. This removes the proposed need for a numerical
LLL absolute gap BEFORE asserting first-order domain stability in this
particular stationary family. It also supplies explicit first/second
massive resolvent bounds without varying the actual projector.

The full physical metric problem remains open. Its canonical nonlocal
shape Hessian is logarithmic order four; the local Weyl/radius action,
auxiliary scalar, off-shell contacts, constraints and mass return must
be included. Neither a matter resolvent graph nor an L2 gap identifies
that Hessian's boundary domain or bounds its inverse in the physical
admission norm. The finite-part ultraviolet subtraction, complete
source derivatives and outer physical tails cannot be estimated by
taking the trace of (11). There is no machine-certified continuum
metric solution in the supporting checks.

Keep the actual-correction budget from the preceding audit:

    1.632e-11 K_m+1.36e-9 K_c+Y_rest<=5e-9,
    rho=1e-8, Z0<=.25, Z2<=5e7.

The K_m,K_c are PHYSICAL metric return norms, not the numbers in (8).
Their sufficient allocations 61.2745098039 and 1.4705882352, and
Y_rest<=2e-9, remain unknown. Semiclassical omission/noise control,
all-sector retarded stability and full wormhole requirements are unchanged.

## 6. Dated correction: use the strongest available outer contract

The preceding [cold feasibility audit](cold-source-feasibility-audit.md)
and [physical return audit](auxiliary-physical-return-feasibility.md)
quoted the older outer background threshold 1e-7 and noted that the
upper allowance 1.01e-7 failed that test. That statement is true for
the OLDER [S146 contract](outer-finite-background-response.md), but it
is incomplete as a description of all retained admission results.

The later [actual-domain outer theorem](finite-outer-admission-and-response.md),
already acknowledged in the recentered canonical audit, proves the
stronger contract N<=1.5e-7. Its N includes the proper-collar mass and
weighted C5 remainder AND the actual global C1 metric norm G/1000.
It supplies regular first response <1.8e-7 and mixed response <.02,
with actual-to-cap separation and both domain inverses accounted for.

For the ORIGINAL finite comparison N<1.01e-7, the margin is therefore
>4.9e-8. Adding a correction with N<=1e-8 leaves N<1.11e-7<1.5e-7.
The first source increment is <1.8e-15 and second-order Taylor remainder
<1e-18 in that theorem's output norm. The correct conditional threshold
at radius 1e-8 is now

    b_o+1e-8 E_o<=1.5e-7;
    E_o<=1 implies b_o<=1.4e-7.                      (14)

If the old 1.01e-7 background allowance and all theorem hypotheses
transfer, E_o<=4.9 suffices (strictness at an open boundary retained).
Thus that numerical old allowance is NOT itself an outstanding outer
admission obstruction under the stronger theorem. The two norms and
their physical domains cannot be conflated by just changing a threshold.

No S169 transfer is claimed here: the new seed coefficients, cap
attachment/response hypotheses, full physical return and unknown
correction still require verification in the stronger N norm. The matter
domain result (9) supplies only part of those hypotheses. The old
complete interior H4 loss and the lack of an N2-to-C5 embedding are
unaffected. Preserve both earlier audits and this dated scope correction.

Remaining to-do list: complete physical source/Hessian graph and actual
coupled correction; outer theorem transfer and admission; nonlinear
preservation, matching, semiclassical validity, all-sector stability,
charge/return, opening, transport, preparation, separate creation and audit.
