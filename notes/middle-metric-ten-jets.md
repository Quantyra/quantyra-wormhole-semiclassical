# Middle metric reconstruction through ten derivatives

2026-10-21, E003/S167. Ordinary argument before supporting verification.
On the ORIGINAL common-disk finite comparison, write

    sigma_ell=sigma0+ell u,
    t_ell=1+(t-1)exp(ell V),
    nu=(1-1/t)V, w_ell=log(t_ell/t), 0<=ell<=1.

On the entire middle region s=t-1>=.001, Y>=.05, the result is

    |nabla^j(ell u)|, |nabla^j w_ell|
        <ell 1e17 j!100^j, 0<=j<=10.                  (1)

At ell=0 the fields vanish exactly. Derivatives and tensor norms in (1)
use the seed's Euclidean orbit metric in A units. These are full ordered
covariant tensors, including angular components; no symmetrization of
higher derivatives is made. Keep the original, much smaller C5 and
independent C0/C1/C2 contracts from the
[endpoint proof](stationary-completed-metric-response-norm.md), with the
[repaired source dependency](repaired-source-metric-transfer.md).
The large allowance in (1) is used only when the higher orders are needed.

This completes the middle region's two logarithmic metric fields. It
does not assert the cold or outer order-ten estimates, derivatives in
the endpoint-dependent proper gauge, a full finite-metric H6 source,
nonlinear existence, or stability. In particular the endpoint proof's
C6 auxiliary cutoff is not silently made C11. The fixed common disk
needs no new cutoff for this calculation.

## 1. Inputs and analytic seed coefficients

Use the original source, reference action, HH state, paired domain and
common subtraction of the [exterior proof](stationary-centered-exterior-source.md).
Set z=1/t, tau=1-z, e=d/C, gamma=Gq^2/(CB^2)<6.8e-28. With outward proper
radius rho, write p=sqrt(b), where

    b(z)=calB(z)/C
      =x/C+z^2/C-(2+4x+delta)z^3/C
        +(1+3x+d)z^4/C-2e z^5+e z^6,
    z'=-p, p'=-b_z/2, a=(log calB)'=-b_z/p,
    a'=b_zz-a^2/2, L=p/z.

Primes on radial fields mean rho derivatives. Subscripts on b mean z
derivatives. The retained bounds give t<50.1, p<.42, |a|<2.02,
1/C<1e-10, x/C<1/5.9, |e|<1e-24 and delta/C<1.001e-60.
The bound for a is twice the real logarithmic-warping bound in the
exterior proof. Also tau>=1/1001. On |z|<=1.02, absolute polynomial
coefficient sums give

    |b_z|<5, |b_zz|<12, |b_zzz|<17, |b_zzzz|<13.

Solve the analytic ODEs on |Delta rho|<=.005. The bootstrap
|p|<.45, .017<|z|<1.02, |a|<4 closes: their displacements are below
.0125, .00225 and .1, respectively. For a use its Riccati equation,
not a lower bound for the small real p. Uniqueness agrees with -b_z/p
where initially defined and continues it without a spurious throat
division. The bounds below therefore imply Cauchy bounds
|c^(k)|<=c_* k!200^k:

| Coefficient | c_* |
| --- | ---: |
| a/2 | 2 |
| L | 40 |
| p z | 1 |
| z^2 | 2 |
| b_zz | 12 |
| z b_zzz/2 | 10 |
| a/2+2L | 100 |
| 2 S, S=-b_z/z+2b/z^2 | 20000 |
| c_nu=z b_zzz/6+4(x/C)/z^2+2e z^3-4e z^4 | 20000 |

These analytic extensions concern explicit seed coefficients only.
No analytic continuation of the nonlocal quantum source is assumed.

The [complete middle source](seed-middle-eight-jets.md) bounds the massive
canonical trace through eight by 1e44 j!100^j q^2. Add the canonical LLL
trace -q C b_zz/(24pi), whose complex value on the above disk is below
q^2. For the combined normalized canonical trace theta, use

    theta_0<=1e20,
    |theta^(j)|/(j!100^j)<=2e44, 1<=j<=8.             (2)

The first entry follows from the old value/first-derivative data
|Theta(2)|<=1e11q^2, |Theta_t|<=1e15q^2 and t<50.1.
Radial derivatives of the massive trace are bounded by its full physical
tensors, since the radial curves are unit-speed geodesics. Thus (2) uses
the same physical normalization as its source theorem.

## 2. Retain the explicit local sources

Let rho_loc and (p-rho)_loc be the Weyl-minus-reference allocations;
the LLL is already included in theta and its canonical Ward pressure.
Define

    E_rho=B^2 rho_loc/(q^2 z^2),
    E_delta=B^2 (p-rho)_loc/(q^2 z^2),
    E_A=4pi B^4 t^2 Arest/q^2,
    E_u=E_A-3z^2 E_delta, E_nu=z^2 E_delta.

On the .005 complex proper disk, the exact local formulas imply

    |E_rho|, |E_delta|, |E_A|<1e16,
    |E_u|, |E_nu|<1e18.                             (3)

Here is a deliberately loose coefficient ledger establishing (3).
Use |t|<60, |z|<1.02, |log(t)+log(mu0 B)|<5, C<q and q>=1.
The Weyl compact Euler identity (20) of the exterior proof has quadratic
coefficient sum below 200, linear sum below 5, and the explicit magnetic
term. F has coefficient sum below 2C and powers -4 through 2. After
multiplying the density by t^2, its absolute Laurent/log coefficient
sum is below 1000q^2, with powers within [-10,4]. Consequently its
value divided by q^2 is below 1000*6*2*60^4<1.6e11, even before the
factor 1/(4pi). The explicit E_F formula and
p_W=-rho_W-F E_F/(2pi) cost less than another 1e12 in this normalization.
The exact reference profiles (35) of the
[horizon proof](stationary-horizon-source-defect.md), with
|B^4 p0A|<3300q^2, |B^4 p_v|<q^2 and |B^4 D_P|<3q^2, contribute less
than 1e10 after the same t^2 rescaling. Both local stress quantities fit
1e14, hence the stated 1e16 reserve.

For Arest use the retained c=1/20, a_an=11/360, b_an=1/30 and magnetic
anomaly, subtracting both reference pressures. The seed curvature
coefficient counts in the [reconstruction proof](stationary-centered-first-order-geometry.md)
give |coeff(B^2 R4)|<61C and derivative coefficient sums 366C and
2562C. In B^4 Box R4=F R4_tt+(F_t+2F/t)R4_t the coefficient sum is
below 9516C^2. The Weyl-square sum is below 1241C^2, and the Euler sum
below 1000C^2. Including the reference profiles, 4pi t^2 B^4 Arest
has coefficient sum below 1e6q^2 and powers within [-10,2]. Its value
divided by q^2 is below 2*1e6*60^2<1e10, hence also below 1e16.
All these are Laurent/log identities, not derivatives of unsigned
real-axis estimates. Cauchy now bounds the normalized kth derivative
of each entry in (3) by its displayed constant times 2^k.

## 3. Exact coupled equations and radial recurrence

Let P=B^2 p_can/(q^2 z^2) and F=f/C, with f the original Eulerian
constraint tangent. Canonical conservation, regularity at the bolt
and the original mass constraint give

    P'+a P=a theta/2,
    F'=-L F-2gamma p z(P-theta+E_rho).                (4)

The initial values at every middle observation obey |P|<1e20 and
|F|<1e-8. For P the unimodal calB Ward integral has total variation at
most three times its observation value; its entire interior interval
has t<50.1 and the old trace value bound. For F use the original
mass-integral bound B integral |rho_rest|dr<1e16q^2 and
|s_R|/(R-B)<3.4e-10: |F|<2gamma*1e16+1e-60<1e-8.

The common-coordinate curvature identities are exactly

    delta Rhat=-8piG r^2 Arest/C-6F
                  +6gamma z^2(theta+E_delta)+nu z b_zzz,
    A^2 t^2 delta R4=-8piG r^2 Arest/C
                  +nu[-20e z^3+36e z^4].

All orbit curvatures are in A units. Substituting these into the
conformal and scalar-radius equations of the endpoint proof yields

    u''=-(a/2)u'+b_zz u-(z b_zzz/2)nu
                          +3F-3gamma z^2 theta+gamma E_u,
    nu''=-(a/2+2L)nu'+2S u+c_nu nu
                          -F+gamma z^2 theta+gamma E_nu. (5)

The Arest terms cancel in the radius equation. This cancellation uses
the same common coordinate in BOTH curvature variations, including the
nonzero d-dependent Lie derivative of R40. The anomaly remains in the
u equation. No source allocation, physical equation or finite
subtraction has been removed.

Normalize scalar radial jets by j!100^j. Put T_0=1e20, T_j=2e44 for
1<=j<=8, g=6.8e-28 and E_j=1e16*2^j. Positive rational majorants for
(4) start at P_0=1e20, F_0=1e-8 and obey, 0<=n<=7,

    P_(n+1)=sum_(i=0)^n 4*2^i(P_(n-i)+T_(n-i)/2)/[100(n+1)],
    F_(n+1)={sum 40*2^i F_(n-i)
       +2g sum 2^i(P_(n-i)+T_(n-i)+E_(n-i))}/[100(n+1)]. (6)

For (5), retain U_j=N_j=1e-4 for 0<=j<=5, supplied by the original
radial ledger before its radial-to-normal conversion. In particular
this is not an inference that an arbitrary chart C5 bound controls
proper derivatives. Define for each n

    C_n(c,X)=sum_(i=0)^n c 2^i X_(n-i),
    D_n(c,X)=sum_(i=0)^n c 2^i 100(n-i+1) X_(n-i+1).

For n=4,...,8, use

    U_(n+2)={D_n(2,U)+C_n(12,U)+C_n(10,N)+3F_n
                      +3g C_n(2,T)+g 1e18 2^n}/[1e4(n+1)(n+2)],
    N_(n+2)={D_n(100,N)+C_n(20000,U)+C_n(20000,N)+F_n
                      +g C_n(2,T)+g 1e18 2^n}/[1e4(n+1)(n+2)]. (7)

Every term is an ordinary Leibniz term; factorial normalization is why
the convolutions have no binomial coefficient. The only new entries
in (7) involve indices at most n+1, so the recurrence is not circular.
Its exact rational evaluation gives upward-rounded bounds

| Order | U_j | N_j |
| --- | ---: | ---: |
| 6 | 4.086e13 | 1.362e13 |
| 7 | 6.042e13 | 2.205e13 |
| 8 | 9.222e13 | 3.825e13 |
| 9 | 1.448e14 | 6.642e13 |
| 10 | 2.329e14 | 1.161e14 |

Thus H=1e15 bounds every high radial tangent coefficient. The loss
from the source allowance 2e44 is offset by gamma and the two radial
integrations. The separate old low jets remain small.

## 4. The ORIGINAL finite radius completion

It is necessary to bound log(t_ell/t), not replace the finite family
by a convenient interpolation. First control tau^-1 on a larger seed
coefficient disk. Write

    b=tau^2 g0(z)-(delta/C)z^3,
    g0=[x(1+2z+3z^2)+z^2+d z^4]/C,
    alpha=p/tau, tau'=alpha tau,
    alpha'=-tau g0_z/2+(delta/C)[z^3/tau^2+3z^2/(2tau)].

On |Delta rho|<=.05, bootstrap |alpha|<3. Then
|tau|<=tau0 exp(.15)<1.2, |tau^-1|<1001exp(.15)<1200 and
|z|<=z0+tau0(exp(.15)-1)<=1. Since |g0_z|<2, |alpha'|<1.3;
the real |alpha0|<1.01 closes the bootstrap. The thermal terms are
below 1e-53. Analytic continuation here may cross z=0; only the
explicit seed coefficient tau^-1 is continued, not the physical domain.
Its normalized Taylor coefficient sum through ten is below

    sum_(j=0)^10 |(tau^-1)^(j)|/(j!100^j)
       <1200 sum_(j=0)^10 (1/5)^j<1500.              (8)

Work in the Taylor algebra through degree ten with variable
100 Delta rho and coefficient l1 norm. Split each field into degrees
0,...,5 and 6,...,10. The latter subspace is an ideal whose square is
zero at this truncation order. For nu the low part has norm <=.0006
and the high part <=5H. Consequently V=nu/tau has low norm <=.9 and
high norm <=7501H. This is a finite-jet statement; no Taylor series
convergence of nu is required.

Set h(v)=(exp(v)-1)/v=integral_0^1 exp(tv)dt. Exactly,

    t_ell/t=1+ell nu h(ell V).

This identity keeps the tau cancellation. Since exp(.9)<2.5, h and
its first derivative at the low V part have l1 norm <2.5. The low
part of X=ell nu h(ell V) has norm <=ell*.0015. Its high part has norm
below ell*(.0015+23.7515H): the terms are the low-only contribution,
5H*2.5 and .0006*2.5*7501H. Products of two high parts vanish through
ten. Thus, for X=X_low+X_high,

    log(1+X)=log(1+X_low)+(1+X_low)^-1 X_high

in this algebra. The inverse norm is at most 1/(1-.0015), proving
the full radial coefficient norm of w_ell is below ell*25H, hence
below ell*3e16. The orbit logarithm ell u has radial reserve ell*H.
The formulas retain ell even at its zero endpoint.

## 5. Full physical tensors and implication

In the seed orbit metric write n=d rho. Then

    nabla n=(a/2)(g-n tensor n), nabla g=0.

Expand each ordered derivative in tensors made from n and g, keeping
the slot order. Their Hilbert norms satisfy |n|=1, |g|=sqrt(2).
Differentiating one n costs at most (1+sqrt(2))|a|/2<3|a|/2. The
radial derivatives of this coefficient have majorants
Gamma_k=6 k!200^k. With initial radial reserve J=3e16, define

    B_(0,k)=J k!100^k,
    B_(r+1,k)=B_(r,k+1)
       +r sum_(i=0)^k binom(k,i) Gamma_i B_(r,k-i), r+k<=9. (9)

This bounds the sum of physical tensor norms times radial coefficient
derivatives. There are at most r factors n at rank r. The radial
derivative of the tensor basis is zero because radial geodesics are
unit speed. Hence (9) covers angular components and all connection
products, including Gamma_0. Its amplification
B_(j,0)/(J j!100^j) is below 3.12 for j<=10; the tenth is below 3.120.
Since 3.12*3e16<1e17, (1) follows with the factor ell retained.

The new quantitative consequence is a middle-region order-ten
coefficient bound for the actual original finite comparison. It does
not yet bound the quantum functional on those finite metrics: the
six-derivative full trace also needs the cold weighted and outer
conormal metric estimates, endpoint maps and correlated mass cubic.
Then the differentiated quantum, tensor-return and Schur-source
estimates are required before the transfer target
R0 R2 R4 N6(p)<=1e57 epsilon_c can be assessed. No such target, signed
compatibility condition or physical stability gate has passed.

## Verification scope

The focused checker verifies exact coefficient margins, curvature/Lie
and anomaly cancellation identities, rational radial and physical
tensor ledgers, the finite-jet algebra margins, independent coordinate
covariant tensors, and independent complex-contour evaluations of the
nonlinear finite radius map at 60/90 digits. It pins its predecessor
and source proofs. These checks support the displayed ordinary
continuum argument; they are not a machine verification of the proof
or an evaluation of the actual charged quantum source. Failed runs,
if any, are retained with frozen source snapshots and raw logs.

Remaining to-do list: cold and outer metric derivatives through ten,
endpoint maps and correlated mass cubic; finite-metric quantum/Schur
estimates; positive or sharper actual transfer; signed compatibility,
nonlinear existence, physical state/boundary and stability; remaining
matching, charge/return, opening, transport, preparation, creation and audit.
