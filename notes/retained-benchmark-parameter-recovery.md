# Recovering the retained radius, vacuum and temperature data

2026-10-19, E003/S141, supporting S139/S074. Ordinary argument recorded
before its focused checker. This corrects a dependency audit, not the
retained physical model.

## Comparison and acceptance scope

The comparison is between [the reference model, equation (1)](uniform-dirac-reference.md),
[the horizon equations](stationary-horizon-product-matching.md),
[the centered vacuum equation (1)](stationary-centered-covariant-seed.md),
and [the classical-radius audit](classical-radius-clock-width.md).
The independent source is
[Bintanja, Freivogel and Rolph, arXiv:2308.00871v2, equations (2.18), (2.24)](https://arxiv.org/html/2308.00871v2#S2.SS1.SSS4).
It supplies the established magnetic RN-AdS radius and charge convention;
no novelty is claimed for recovering that radius. The source was inspected
online on 2026-10-19.

The hardest bookkeeping issue is that `lambda0=L4^-2`, the corrected
`lambda_e=ell_outer^-2`, a bookkeeping homotopy, and the dimensionless
thermal gap were all called `lambda` in different notes. Here they are
distinct. Likewise the benchmark integer is `K=2000`, the scalar radial
wavenumber is `k_rad`, and neither is the Killing surface gravity `kappa`
or the outer normal-symbol constant `kappa_M`.

Acceptance requires direct horizon substitution, rational enclosures of
the derived constants, and a separate high-precision numerical root solve.
Before execution, require a `kappa_M` enclosure narrower than `10^-20`,
strict rejection of the old spectral point by the vacuum equation,
positive finite temperature bounds over the entire retained preparation
interval, and `18<M_1<19` for the retained product-radius bounds.
These are parameter checks. They do not test a variable Dirac vertex,
construct the exterior or decide its signed source condition.

## 1. The classical radius is already determined

Keep precisely the inputs of uniform-reference equation (1):

    G=l_p^2, K=2000, q=K^4, g_R(mu0)=K^-4,
    L4=sqrt(G) K^5, lambda0=1/(G K^10), eta=1/2000,
    Q=pi G q^2/g_R^2=pi G K^16.

Here Q is the metric charge coefficient, equal to the source paper's
`r_e^2`, not its dimensionless magnetic charge. For

    F0(r)=1+lambda0 r^2-2m0/r+Q/r^2,

the equations `F0(b0)=F0'(b0)=0` give

    Q=b0^2+3lambda0 b0^4,
    m0=b0(1+2lambda0 b0^2),
    C0=1+6lambda0 b0^2=sqrt(1+12pi K^6),
    xi0=(C0-1)/6,
    b0^2=G K^10 xi0, a0^2=b0^2/C0, mu0=1/b0.        (1)

The function `y+3lambda0 y^2` strictly increases on `y>=0` from
zero to infinity. Thus there is exactly one positive classical `b0`.
The expression in (1) is an exact numerical specification after fixing
length units, even though it contains pi and radicals. It was explicitly
present in uniform-dirac-reference.md before S140. No new measurement,
choice of charge, lower-endpoint substitution or corrected radius is needed.
The renormalized coupling at this reference scale is a retained model
input, not a measured coupling of a laboratory system.

Direct factorization gives the independent check

    F0(r)=(r-b0)^2/r^2
              [1+lambda0(r^2+2b0 r+3b0^2)],
    F0''(b0)/2=C0/b0^2=1/a0^2.                       (2)

In particular the exterior is positive for `r>b0`. This is the classical
reference; it is not a wormhole and is not the semiclassical exterior.

## 2. The corrected outer curvature is not a free sample in its bounds

The retained outer prescription fixes the positive root of

    lambda_e+A lambda_e^2=lambda0, A=11G/(360pi).

It has the cancellation-free exact expression

    lambda_e=2lambda0/[1+sqrt(1+4A lambda0)],
    delta_e=1-lambda_e/lambda0=A lambda_e^2/lambda0,
    0<delta_e<A lambda0=11/(360pi K^10).              (3)

Uniqueness follows because the left side increases strictly for positive
curvature. At K=2000 the fractional shift is below `10^-34`.
This specifies the vacuum curvature, conditional on the already selected
AdS-invariant outer prescription. It does not establish that a full
monopole HH solution approaches that vacuum.

The S140 point `G=1, lambda_test=(1999/2000)K^-10`
lies in a loose enclosing interval but fails this equation. Its normalized
residual is

    [lambda_test+A lambda_test^2-lambda0]/lambda0
       =-1/2000+[11/(360pi K^10)](1999/2000)^2.       (4)

The last term is less than `10^-34`, so the residual is strictly negative
and is not a rounding error. A guard interval is not permission to vary
an algebraically fixed input. Equation (4) rejects the spectral test
point as a realization of this retained model. Its abstract spectral
arithmetic can remain valid on its own larger parameter set.

## 3. Normal-symbol constant and the physical temperature family

With the same classical subtraction scale, equations (1)-(3) give

    kappa_M=-6/5-log(mu0 ell_outer)
           =-6/5+(1/2)log(lambda_e b0^2).            (5)

It is a single determined dimensionless constant. The old width 3 and
the intersection width `(3/2)log(2)` were widths of coarse enclosures;
they are not physical freedom in this model. A finite numerical enclosure
still has a rounding width. The mathematical parameter has no unspecified
degree of freedom. Neither observation equates `kappa_M` with `kappa`.

The retained cold-bolt equation (2) also fixes

    nu=2q eta G C0/b0^3,
    kappa(z_i)=z_i nu, 1/4<=z_i<=7/25,
    beta(z_i)=2pi/kappa(z_i).                         (6)

Consequently the prescribed temperature and Euclidean period are a
defined positive finite family. No new single temperature should be
chosen: S139 explicitly retains the full preparation interval and fixed
laboratory waveform. Its envelope is

    nu/4<=kappa<=7nu/25,
    50pi/(7nu)<=beta<=8pi/nu.                         (7)

G=1 in the receipt means Planck length units. Restoring dimensions,
`b0,a0,beta` scale as `sqrt(G)`, `mu0,nu,kappa` as `G^-1/2`,
and `lambda_e` as `G^-1`. The scalar phase constants are also fixed:

    E=60pi/G+11lambda_e/3,
    k_rad^2=E-9lambda_e/4=60pi/G+17lambda_e/12,
    c0=k_rad(E+4lambda_e).                           (8)

Thus `c0` scales as `G^-3/2`. Its factor-eight change between G=1
and G=4 at fixed dimensionless benchmark is the predicted scale
dependence; it does not show that there is no defined parameter domain.
This does not calculate the separate source field or its moment norm.

## 4. Restoring the proper angular gap on the retained product box

The reference has `M_n^2=n(n+q)/C0`. The corrected homogeneous product
has `B^2=z_star b0^2`, `A_star^2=u_star a0^2` and

    C_star=C0 z_star/u_star,
    M_n^2=n(n+q)/C_star,
    |z_star-1|<1.4e-22, |u_star-1|<4e-26.             (9)

These are the existing homogeneous theorem's bounds; the checker
propagates them and does not re-prove the determinant theorem.
Writing `dz=14/10^23`, `du=4/10^26`, one has

    C0(1-dz)/(1+du)<C_star<C0(1+dz)/(1-du).

The rational enclosure of these endpoints proves `18<M_1<19`.
Moreover `(n+q)/(1+q)>=1` for every integer n>=1, so
`M_n>=M_1 sqrt(n)>18 sqrt(n)` for the entire tower. This is a
proper-mass gap, not a gap for the optical Hamiltonian.

The prior witness `x=10^10,d=0` shows the weakness of the loose guards
`x>7.9e9, |d|<10^-12`. The actual product has
`x=lambda_e b0^2 z_star<8.2e9`, so that witness is outside (9).
It therefore does not refute the actual reference gap. Propagating (9)
also recovers `100<1.3e8/(60pi M_1^3)<124`; this formal scalar
heat-coefficient product remains greater than one and supplies no
bound for the variable Dirac three-vertex. The nonlinear obstruction
is not removed by recovering the proper parameter domain.

## 5. Consequence for the research task

S140 correctly found no literal numeral in its selected files and
computed the widths of their loose bounds. Its inference of an
unavailable classical radius was unsupported: its input list omitted
the linked note explicitly defining that radius, and its regular
expression could not recognize an exact algebraic dependency anyway.
The old files and receipts remain preserved as historical evidence.
The new result supersedes that physical-input classification.

The previous goal checkpoint made progress by retaining an auditable
calculation, but its conclusion does not survive the dependency check.
This turn changes the next action: no user-supplied radius or arbitrary
temperature selection is needed. The next task is the actual quantum
exterior response at these retained data, with the full preparation range.

The source moment `m_num`, the auxiliary source norm, and the seam trace
are not set to zero or enclosed here. In the scalar-phase note they are
derived from a fixed source and operator: absence of their computation
is not in itself evidence that an external physical measurement is
required. Their computation and the nonlinear mapping estimates remain
research obligations. Finite EFT coefficients outside the explicitly
selected reference prescription retain their separate matching status.

Remaining to-do list: construct and bound the actual exterior response,
its signed source condition, variable Dirac vertices and nonlinear
remainder; establish the boundary completion and charge; transfer the
driven state/stress and geometry estimates; prove the specified stability,
finite-energy transport and preparation results. The area budget remains
unevaluated. S139.1-S139.2 and S074.2-S074.6 remain open; no full
completion-plan outcome or stable semiclassical wormhole is proved.
