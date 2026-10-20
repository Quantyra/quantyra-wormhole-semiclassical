# Scalar reference resolvent on the actual finite comparison path

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary proof and
verification contract before execution. **The scalar reference
resolvent is controlled at the finite comparison; the full candidate
Hessian inverse and nonlinear solve remain unproved.**

This result uses the actual finite geometry, rather than identifying
its Laplacian with the seed operator. Keep the original physical
temperature, clock, flux, state/domain and comparison path of
[the metric completion](stationary-completed-metric-response-norm.md).
The path is

    sigma_ell=sigma0+ell u,
    r_ell=B+(r0-B)exp(ell V), 0<=ell<=1,
    t=r0/B, tau=(t-1)/t, nu=tau V,
    T=1+tau[exp(ell V)-1], S=log T, w=ell u+S.

Here T is a radius ratio, not a Laplacian. In seed proper radius rho,
with angular time theta of period 2pi and a0=sqrt(F0)/kappa,

    g_ell=e^(2w)(d rho^2+a0^2 d theta^2)+r0^2 e^(2S)dOmega^2,
    d rho_ell=e^w d rho, a_ell=e^w a0, r_ell=e^S r0.

The fixed physical period is beta=2pi/kappa. Regular even u,S and
the new proper radius give a_ell'(0)=1; no temperature is changed.
Use dimensionless proper coordinates y=sqrt(lambda)rho and
z=sqrt(lambda)rho_ell, so z(y)=integral_0^y exp(w(v))dv.

Let Tbar_ell=(-Box_(g_ell))/lambda on static spherical scalars,
with its actual regular bolt realization and physical L2 measure

    dmu_ell=8pi^2 a_ell r_ell^2 d rho_ell.

The unitary map to L2(dz) multiplies by
sqrt(8pi^2 a_ell r_ell^2/sqrt(lambda)). The following estimates hold
uniformly on the stated finite comparison path:

    U_ell Tbar_ell U_ell^-1
       =H_B0+9/4+q_ell,
    H_B0=-partial_z^2-1/(4z^2), regular critical realization,
    ||q_ell||_infinity<400,
    ||q_ell||_L1(dz)<72000.                          (1)

At the retained fast energy Ebar=60pi/g+11/3, g=G lambda, put
kbar=sqrt(Ebar-9/4)>4e17. Then

    ||q_ell||_1/kbar<2e-13,
    |G_ell(kbar;z,z')|<1.001/kbar.                  (2)

G_ell is the actual outgoing kernel of Tbar_ell-Ebar in these
unitary variables, not the full coupled Hessian inverse. It has the
regular bolt condition and no nonzero regular outgoing homogeneous
solution. Weighted kernels with L2 multiplication weights have the
same limiting-absorption and Hilbert--Schmidt conclusions as the
retained seed Bessel construction. Their bound is
`1.001||omega1||_2||omega2||_2/kbar`.

## 1. Low metric information used, and its coordinate conversion

Only two spatial derivatives of the actual path enter this proof.
The [completed metric estimates](stationary-completed-metric-response-norm.md)
give |u|,|V|<1.2e-8, |b|=|(t-1)V_t|<3e-7 and the independent
small low relative coefficients. On the cold and middle cover Y>=.05,
use the deliberately weakened radial reserves in canonical proper
distance zeta (in A units)

    |u_zeta|,|nu_zeta|<1e-4,
    |u_zetazeta|,|nu_zetazeta|<1.                   (3)

These follow from the explicit low radial estimates: the middle
second jets are below 8e-4 and 1e-3. On the cold patch the radius
two-jet allowance is at most `(1e-5)(.00331)(2!100^2)=.000662`;
the orbit Hessian is smaller. The first-two-jet radial/normal
conversion has factor below two, with its lower normalization
terms inside the separately retained .01 coefficient reserve.
Even `2(.000662+.01)<1` fits (3). No fifth-jet allowance is used
as a bound for these low derivatives. The outer guard z_seed/Y<.412
gives Y>1/(.412 t), so t<=40 is within Y>=.05 wherever the
outer collar is needed to check this cover.

Since d zeta/dy=sqrt(C/x)/t, with x=lambda B^2 and 5.9<C/x<6.1,
and |d_y log t|<=sqrt(5), (3) gives

    |u_y|,|nu_y|<.00025,
    |u_yy|,|nu_yy|<6.2.                             (4)

The actual scalar-radius completion is also controlled, not just its
tangent. Write k=exp(ell V), B_ell=k-1-ell V k. Exact differentiation
in y gives

    T_y=ell k nu_y+B_ell tau_y,
    T_yy=ell k nu_yy+B_ell tau_yy
                                  +tau k ell^2 V_y^2. (5)

The seed guards imply |tau_y|<=sqrt(5), |tau_yy|<22 and
|V_y|<=sqrt(7)|b| on t<=40; the last follows from F0/x<=7(t-1)^2.
Also |B_ell|<2e-16, |k|<1.001 and T>.999. Thus S_y=T_y/T and
S_yy=T_yy/T-(T_y/T)^2, with (4), imply the generous global core
bounds

    |w|,|S|<epsilon=3e-8,
    |w_y|,|S_y|<.001,
    |w_yy|,|S_yy|<100,       t<=40.                 (6)

Both first derivatives vanish at the regular bolt. In particular
|S_y(y)|<=100y, a fact needed to treat the polar coefficient.

Here is the outer information in the same fixed coordinates. The
common-coordinate conversion and finite endpoint ledger in the metric
proof, through order TWO only, allow for each f=ell u,S

    f=a3 Y^3+R, |a3|<1e-9,
    |R^(j)|<=1e-6 j!100^j Y^(4-j)L(Y), j=0,1,2,
    L(Y)=1+log(1/Y), Y<=.1.                         (7)

The retained common tangent remainder is below 2e-8 and its cubic
below 1e-11. For S the exact log-gap completion adds quadratic terms
starting at Y^6. Its denominators T and 1-1/t exceed .999 and .93
on this collar. Here is an explicit two-jet budget for that finite
step. Put q=ell nu and

    S=q+q^2 H(q,tau),
    H(q,tau)=integral_0^1 (1-v)
      (1-tau)exp(vq/tau)/[tau(1+tau(exp(vq/tau)-1))^2]dv.

The tangent bounds imply `|Y^j q^(j)|<Q_j Y^3`,
Q=(1e-8,1e-6,2e-4), j=0,1,2, using Y L(Y)<.34.
The baseline analytic guards give |Y tau'|<.05 and
|Y^2 tau''|<.001. On |q|<.01 and the complex .1-neighborhood
of .93<=tau<=1, |tau|>.83 and the displayed denominator factor
has modulus above .986. Thus |H|<1. Cauchy disks of radii
.005 and .05 about the real arguments bound

    |H_q|<200, |H_tau|<20,
    |H_qq|<80000, |H_qtau|<4000, |H_tautau|<800.

The two product/chain rules then give
`|partial_Y^j(q^2 H)|<1e-10 Y^(6-j)` for j=0,1,2.
This correction's N_(4,2) coefficient is below 1e-12, which
fits (7) together with the tangent's 2e-8 coefficient. Its leading
order is six, so the correlated cubic is unchanged. No moving
proper-coordinate map is differentiated.

For t>=40 the baseline outer identities give

    1/t<Y<3/t<.1, |Y_y|<1.1Y, |Y_yy|<3Y,
    L(Y)<=1+log t, (1+log t)/t<.12.

They follow from Y=sqrt(C) integral_0^(1/t) calB^-1/2 and the
strict outer guards in the metric proof. Applying the first/second
chain rules to (7), and adding the two fields for w, gives

    |w|,|w_y|,|w_yy|,|S_y|,|S_yy|<10t^-3, t>=40.  (8)

For example the largest numerator, for w_yy, is bounded by
32.52(1e-9)Y^3+49000(1e-6)Y^4 L(Y), which is below
`[32.52(1e-9)27+49000(1e-6)81(.12)]t^-3<10t^-3`.
The global low value bounds still give |w|,|S|<epsilon there.
Equations (6),(8) in particular give the global first/second reserves
.001 and 100. The core and outer covers overlap, so no seam is omitted.

## 2. Exact transformed potential, including the bolt cancellation

Set J0=d_y log(a0 r0^2) and V0=J0_y/2+J0^2/4. The retained seed
estimates give

    J0=1/y+j0, 0<=j0<8,
    V0=-1/(4y^2)+9/4+q0, |q0|<30,
    |V0-9/4|<=11/t^2 for t>=2.

For the actual metric J_ell=d_z log(a_ell r_ell^2)
=exp(-w)(J0+w_y+2S_y). Direct differentiation yields

    V_ell=exp(-2w)[V0+w_yy/2+S_yy+J0 S_y
                                      -w_y^2/4+S_y^2]. (9)

In particular the J0 w_y term and the mixed w_y S_y square cancel.
Estimating the pieces before that cancellation would give an incorrect
polar bound. Formula (9) uses the actual four-dimensional scalar
Laplacian; it is not the orbit two-dimensional Laplacian.

Let R(y)=z(y)/(y exp(w(y))). Global |w|<epsilon and w_y(0)=0,
|w_yy|<=K=100 give

    exp(-2epsilon)<=R<=exp(2epsilon),
    |R-1|<=exp(2epsilon)K y^2/3.

The second estimate follows by integrating
|w(v)-w(y)|<=K(y^2-v^2)/2 before exponentiation. Consequently

    |1/(4z^2)-exp(-2w)/(4y^2)|
        <=exp(10epsilon)K/6<17.                    (10)

This is a bound on the cancellation of the two singular terms,
not separate bounds on either of them. Also (6) and J0<=1/y+8 give
|J0 S_y|<=100+.008. Inserting these estimates in (9), with
q_ell=V_ell+1/(4z^2)-9/4, proves |q_ell|<400 on t<=40.
It has a bounded continuous bolt limit by the even metric jets.

The derivative conjugation has the same physical regular realization.
Indeed J_ell-1/z is bounded. Its first contribution is
exp(-w)(j0+w_y+2S_y), of size below 8.02. For the remaining
exp(-w)/y-1/z use both |R-1|<=exp(2epsilon)K y^2/3 and
|R-1|<=exp(2epsilon)-1<3epsilon. Their geometric mean bounds
that term by less than .002. Hence

    |J_ell/2-1/(2z)|<5.                             (11)

Thus the physical derivative differs from the critical Bessel
derivative by bounded multiplication. The scalar form domain is the
same critical regular domain; bounded q_ell identifies its operator
realization as well. In particular the log bolt branch is excluded
without imposing zero scalar value at the bolt.

## 3. Full potential integral and outgoing inversion

The retained seed proper lengths obey .2<y_2<172. Since
F0/x>=t^2/2 for t>=2, integration to t=40 gives y_40<177.
Consequently

    .199<z_40<178.

The core contribution to ||q_ell||_1 is below 400(178)=71200.
For the outer contribution, (8),(9) and J0<13 imply

    |V_ell-9/4|<=12/t^2+200/t^3,
    dz<=1.5 dt/t, t>=40.

The integral of the displayed potential bound is less than .01.
The remaining critical term has integral 1/(4z_40)<1.3. Thus

    ||q_ell||_1<71201.31<72000.

The same outer bounds also give |q_ell|<400 there. Restoring physical
units would give potential L1 norm below 72000sqrt(lambda); (1)
and (2) use the dimensionless proper coordinate consistently.

Reuse the exact critical Bessel kernel and its upper-half-plane
bound from [the outgoing proof](stationary-scalar-outgoing-inverse.md),

    |G_B0(zeta;z,z')|<=1/|zeta|, Im zeta>=0.

Factoring the real bounded potential q_ell into its square-root
absolute value and sign gives a sandwiched operator of norm at most
||q_ell||_1/|zeta|. The Neumann series therefore constructs the actual
outgoing resolvent near kbar, with kernel bound

    |G_ell|<=1/[|zeta|(1-||q_ell||_1/|zeta|)].       (12)

The retained kbar>4e17 proves (2). Dominated convergence with any two
L2 multiplication weights gives the operator-norm outgoing limit.
The same invertible sandwich excludes a regular outgoing homogeneous
solution. All angular quantum and TT terms remain outside this scalar
reference statement; none is discarded from the full physical action.

The exact reference fourth-order identity is also valid for this SAME
finite-comparison scalar operator:

    [(Tbar_ell+4)(Tbar_ell-Ebar)]_out^-1
      =[(Tbar_ell-Ebar-i0)^-1-(Tbar_ell+4)^-1]/(Ebar+4).

The positive elliptic term is retained. This is a reference resolvent
identity, not a factorization of the full candidate Schur complement.

## 4. Consequence and limits

The actual finite comparison's scalar reference block remains suitable
for the outgoing construction. Its bound is uniform along 0<=ell<=1
and retains the physical cold length and regular bolt. It does not
require spatial analyticity of the finite metric or a sixth derivative.

The actual TT positivity, mixed vertices, density contacts and their
quantitative transfer to the finite comparison still require proof.
The candidate's adjoint and coupled range must be treated with those
blocks, rather than silently replacing them by the seed quantities.
No full candidate Hessian inverse, fixed projection-direction norm,
cold correction bound, preserved nonlinear family or projected solution
is certified. The separate signed compatibility equation remains.

## Frozen verification contract

Freeze the note and checker before execution. Reconcile the completed
metric, seed outgoing and current S167 evidence/input identities.
Verify the physical coordinate conversion, exact finite radius jets,
Liouville potential cancellation, singular-term integral estimates,
all low-jet/outer/core budgets and the final small scattering ratio.
Use independent symbolic finite metric warps to check the transformed
potential and regular bolt expansion. They are geometry fixtures, not
solutions or evaluations of the physical quantum determinant.
Continuum form domains, kernel limits and inverse statements remain
ordinary proofs. Preserve every terminal run and failed source if any.
No Lean/cloud campaign is selected.

Remaining to-do list: finite-comparison coupled blocks, fixed projection
direction and cold/higher inverse/source bounds; nonlinear preservation
and projected full solve; signed compatibility; exterior/state, charge
and return; opening, stability, transport, preparation, separate creation
and full-result audit.
