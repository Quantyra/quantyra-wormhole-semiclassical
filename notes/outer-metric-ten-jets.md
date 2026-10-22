# Outer finite-comparison metric derivatives through ten

2026-10-21, E003/S167. Ordinary argument before supporting verification.
Keep the original seed, paired domain, HH state, physical action,
common subtraction and original finite comparison

    sigma_ell=sigma0+ell u, t_ell=1+(t-1)exp(ell V),
    nu=(1-z)V, z=1/t, w_ell=log(t_ell/t), 0<=ell<=1.

Use the SEED proper outer coordinate

    Y=sqrt(C) integral_0^z calB(v)^(-1/2)dv,
    L(Y)=1+log(1/Y), 0<Y<=.1.

There are the fixed mass coefficients u3,nu3, with u3+2nu3=0,
such that the original finite fields satisfy

    |nabla^j[ell u-ell u3Y^3]|,
    |nabla^j[w_ell-ell nu3Y^3]|
        <ell 1e14 j!100^j Y^(4-j)L(Y), j<=10.       (1)

These are full ordered physical covariant tensors in the seed orbit
metric in A units. At ell=0 the differences vanish. The cubic is
subtracted as the scalar Y^3 BEFORE taking the tensor derivative.
Keep the original smaller low C5 and C0/C1/C2 contracts.

Together with the [cold](cold-metric-ten-jets.md) and
[middle](middle-metric-ten-jets.md) results, (1) quantifies all three
regions of the original logarithmic metric fields through ten.
This note also constructs quantitative forward conormal C11 proper
maps and a C11 auxiliary global identification. It does not bound
their higher inverse/composition operators, the proper endpoint
chords through ten, finite-metric H6, or a nonlinear solution.

## 1. Keep the singular-source cancellation before integration

The [outer source theorem](seed-outer-eight-jets.md) gives ordinary
weighted derivatives through eight of Theta-c C^2/(5pi Y), with
allowance K q^2 L, K=1e44. Its smaller original C3 bounds remain
available. The [endpoint metric proof](stationary-completed-metric-response-norm.md)
supplies the exact analytic seed guards on |Y|<=.15:

    |calB/x-1|<.0011, |z/Y|<.412, |Y/z|<3,
    |z_Y|<.42, |Y z_Y/z|<1.002,
    |(calB_Y/calB)/Y|<1.

Consequently the normalized ordinary jets of an analytic multiplier
of modulus M on a disk |Y-Y0|<=Y0/2 are bounded by M*.02^j.
The full coefficient sum is at most M/(1-.02), independently of
the derivative order. Every multiplier factor is retained below.

Put L0=-a_seed x/(5pi), l0=L0/q^2 and R=Theta/q^2-l0/z.
The difference between the two principal profiles is analytic at
Y=0. The seed guards bound it by the local reserve 1e8; thus

    N_(0,8)(R)<3K,
    N_(m,J)(F)=max_(j<=J) sup |F^(j)|/[j!100^j Y^(m-j)L].

The actual Ward integral, with its fixed bolt constant, obeys

    I(Y)/q^2=I0/q^2+integral_0^Y calB_Z[l0/z+R]dZ,
    |I0|<4e15 x q^2,
    B^2 rho_can/q^2=z^2 I/(2calB q^2)-l0 z-z^2 R,
    B^2(p-rho)_can/q^2=l0 z+z^2 R.                 (2)

The Weyl density contains its exact vacuum term plus l0 z; its
pressure difference contains minus twice the vacuum term minus l0 z.
The reference cancels the vacuum terms. Therefore the l0 z terms
cancel exactly in BOTH complete physical sources in (2), before
forming absolute norms or integrating mass. No divergent sector
mass is integrated separately. This is the same identity as in the
[complete exterior source](stationary-centered-exterior-source.md).

The remaining local Weyl/reference/LLL terms are finite sums
Y^2[A(Y)+B(Y)log Y] with analytic A,B on |Y|<=.15. The RN Weyl
Laurent/log coefficient sums 2076 and 3730, the degree-six correction
bound below 200 per quadratic amplitude, and the reference bound
4pi(3301+12)<5e4 give analytic coefficient reserves below 1e6 after
the bounded z/Y, calB/x and log(z/Y) factors. The LLL is rational
there with the same nonvanishing calB denominator and C/q<.00352;
its HH constant is fixed and fits this reserve. Weighted logarithmic
derivatives have coefficients at most 1/(j100^j), j>=1. Together
with the analytic coefficient sum and the Y^2 polynomial factor,
these local terms fit 1e8 in N_(2,8). These are differentiated exact
Laurent/log formulas, not differentiated unsigned remainder bounds.

For completeness, the nonlocal terms have the same order-independent
product bounds as the old C3 calculation:

* z^2 R has N_(2,8)<.63K, since the multiplier coefficient sum
  is below .21 and R has allowance 3K.
* I_R=integral calB_Y R has N_(2,8)<2.5 C K q^2.
  The value integral uses integral_0^Y yL dy=Y^2(L/2+1/4);
  derivative j>=1 divides the preceding integrand jet by 100j.
* z^2 I_R/(2calB q^2) has N_(4,8)<1.6K and hence N_(2,8)<.016K.
* The constant I0 contribution is below 4.2e14 in N_(2,8).
  The principal I integral is analytic after its explicit l0/z
  factor and fits the local reserve already included.

It follows that the COMPLETE sources satisfy

    Rrho=B^2 rho_rest/q^2, Rdelta=B^2(p-rho)_rest/q^2,
    N_(2,8)(Rrho), N_(2,8)(Rdelta)<2K.              (3)

The trace anomaly minus reference has its exact outer vacuum
canceled as well. The remaining local polynomial begins at z^4.
The coefficient counts in the middle proof, now divided by Y^2,
give N_(2,8)(-4piG r^2 Arest/C)<gamma 1e8, with
gamma=Gq^2/(CB^2)<6.8e-28. Its analytic coefficient jets use .02^j.

## 2. Exact mass mode and the proper tangent

Write m=m1/(xB), where m1 is the original fixed geometric mass
variation, |m|<4e-11. Separate it BEFORE subtracting a cubic:

    f=-2xz[m+mu(Y)], mu(0)=0,
    mu_Y=-gamma_x z_Y Rrho/z^2, gamma_x<4e-27,
    |mu|<7e-11 YL.                                (4)

The last value bound is the old sharp source-value result and is
not replaced by K. From (3), the analytic multiplier
z_Y(Y/z)^2 and its normalized coefficient sum give

    |mu^(j)|/[j!100^j Y^(1-j)L]
          <=8 gamma_x K/(100j), 1<=j<=9.           (5)

Let j=f/F0 and Jp=integral_0^Y j. The exact linear proper-gauge
radius and inverse orbit-warping directions are

    w=-z_Y Jp/(2z),
    v=-a_lapse-j/2-calB_Y Jp/(4calB).              (6)

Use subscripts M and T for their exact mass-mode and tail parts.
The mass mode has a_lapse,M=0 and

    j_M=-2xz^3 m/calB, J_M=integral j_M,
    j3=-2(x/C)^(3/2)m, w_M,3=-j3/8, v_M,3=-j3/2. (7)

The tail w_T is obtained by the two integrals in (4),(6). For v_T
use the invariant curvature equation, retaining its proper-coordinate
Lie term, rather than differentiating the lapse a_lapse ten times:

    v_T''+(calB_Y/calB)v_T'
       =-3 f_T/C+3gamma Rdelta-4piG r^2 Arest/C
                                      +(z b_zzz/2)w_T,
    b=calB/C.                                    (8)

The exact mass mode satisfies its own homogeneous-source curvature
identity and has already been removed from (8). Here |z b_zzz/(2Y)|<4
on the analytic collar, and (calB_Y/calB)/Y has modulus below one.
The older proper tail ledger, before adding the analytic mass
remainder, is below 1e-9 and safely fits N_(4,5)<1e-8 for v_T,w_T.
We retain this low estimate separately throughout.

## 3. Explicit ordinary weighted jet ledger

All following arrays are positive normalized majorants. Put

    A(c)_i=c*.02^i,
    P(m)_i=|binom(m,i)|/100^i,
    (X*Z)_j=sum_(i=0)^j X_i Z_(j-i).

For m=-1, P(-1)_i=.01^i. Multiplication by an analytic coefficient
and by Y^m is therefore encoded by A(c)*P(m). If X bounds a field
of weight Y^m L, its primitive has weight Y^(m+1)L and majorant

    I_m(X)_0=X_0 (m+2)/(m+1)^2,
    I_m(X)_j=X_(j-1)/(100j), j>=1.                 (9)

The value factor is the exact integral of y^m L bounded using L>=1;
the high derivatives use the fundamental theorem of calculus, not a
false uniform two-order gain for value norms.

Take mu_0=7e-11 and mu_j from (5), through nine. The analytic
coefficients have moduli

    |-2(x/C)(z/Y)|<.14,
    |-2x(z/Y)^3/calB|<.15,
    |Yz_Y/(2z)|<.51.

Thus f_T/C, j_T, J_T and w_T are bounded respectively by

    F=A(.14)*P(1)*mu,
    Jsmall=A(.15)*P(3)*mu,
    J=I_4(Jsmall), W=A(.51)*P(-1)*J.               (10)

Only Jsmall through nine is used to obtain W through ten. Keep the
smaller old W_j<=1e-8 for j<=5. Set V_j=1e-8 for j<=5, for the
proper v_T majorant. For n=4,...,8, (8) gives

    V_(n+2)={3F_n+6gamma K+gamma 1e8*.02^n
       +.001[A(4)*P(1)*W]_n
       +.01[A(1)*P(1)*Dv]_n}/[1e4(n+1)(n+2)],
    Dv_j=100(j+1)V_(j+1).                         (11)

The factors .001 and .01 are Y^3 and Y^2, respectively, after
matching the weights of the Lie and drift terms to Y^2 L. This
recurrence only uses earlier V entries. Exact evaluation gives
max V<1.372e12 and max W<9.098e11.

## 4. Return to the original common fields and retain the cubic

Let U(Y)=integral_0^Y u, and G0=sqrt(calB/x). The linear gauge
identities, without differentiating the path in moving proper gauge,
are

    U=-G0 integral_0^Y v/G0,
    u=-v+(calB_Y/(2calB))U,
    nu=w-(z_Y/z)U, V=nu/(1-z).                    (12)

Both G0 and its inverse have modulus below 1.001 on the analytic
collar; (1-z)^-1 has modulus below 1.07. Apply (9) and convolution
to the tails in (12):

    Utail=A(1.001)*I_4[A(1.001)*Vproper],
    u_tail<=Vproper+.01 A(.5)*P(1)*Utail,
    nu_tail<=W+A(1.01)*P(-1)*Utail,
    V_tail<=A(1.07)*nu_tail.                       (13)

The maximum in these original-field tail arrays is below 1.372e12.
The corresponding exact analytic mass modes divided by Y^3 have
modulus below 1e-11 on |Y|<=.15. For example
|j_M/Y^3|<.15|m|, |J_M/Y^4|<.15|m|/4, and (6),(12) give this
bound for every field. Subtract their actual cubic coefficients.
Schwarz applied to F_M/Y^3-F_M,3 bounds the quotient by Y by
2e-11/.15; subsequent analytic coefficient and P(4) factors give
N_(4,10)<1e-9 for each analytic mass remainder.

Consequently, after subtracting their cubics, the original u,nu,V
have ordinary weighted reserve H=1e13 through ten. Their old
low remainders through five fit 1e-7. The exact cubic coefficients are

    u3=j3/2=-m(x/C)^(3/2),
    nu3=V3=-j3/4=m(x/C)^(3/2)/2, u3+2V3=0.       (14)

They are written in the seed proper Y, not the common disk defining
function. Changing defining functions would change their numeric
coefficients and requires its own transfer. The proper coefficients
in (7) have the different relation v_M,3=4w_M,3.

## 5. The original finite radius map and physical tensors

Use Taylor coefficients Y^j F^(j)/(j!100^j) at each real Y, through
degree ten. Low means degrees zero through five; the ideal of higher
degrees has square zero at this truncation order. The low coefficient
sum for each of u,nu,V is at most

    epsilon(Y)=1.1e-11 Y^3+6e-7 Y^4 L
               <2.1e-7 Y^3<1e-9,

and its high sum is at most 5H Y^4 L. The bounds YL<.331 and
Y^4L<.000331 hold throughout this collar. Exactly,

    t_ell/t=1+ell nu h(ell V), h(v)=(exp(v)-1)/v.

On the low part, h,h' and the logarithmic inverse factors have norms
below 2. A product containing two high parts vanishes. The nonlinear
remainder w_ell-ell nu therefore has coefficient norm bounded by

    ell[10 epsilon(Y)^2+100 epsilon(Y)H Y^4L].     (15)

The first term retains its Y^6 factor; it is not bounded by a constant
before division by Y^4L. Dividing (15) by Y^4L gives less than
ell*1e-6 H. Thus the full ordinary remainder of w_ell, and that of
ell u, fit ell*3e13 in N_(4,10). Their cubic coefficients are precisely
ell nu3 and ell u3; the physical lapse's cubic remains zero because
2u3+4V3=0. No auxiliary interpolation replaces the original family.

For the full physical tensors, the seed orbit metric has
n=dY and nabla n=beta(g-n tensor n), beta=calB_Y/(2calB).
On the radius-Y/2 analytic disk, |beta|<.75Y. The tensor-product
factor below three and Cauchy give coefficient bounds
(9/4)Y i!(2/Y)^i. After normalizing a rank-r-plus-k derivative by
Y^(4-r-k)L, each connection insertion gains Y^2<=.01. Hence the
ordered majorant recurrence is

    B_(0,k)=k!100^k,
    B_(r+1,k)=B_(r,k+1)
       +r sum_(i=0)^k binom(k,i)(9/400)i!2^i B_(r,k-i). (16)

Its amplification B_(j,0)/(j!100^j) is below 1.002 through ten.
The zeroth connection and all nonsymmetric tensor components are
included. Since 1.002*3e13<1e14, this proves (1).

## 6. Forward C11 proper maps and an explicit higher cutoff

For each fixed ell use the exact forward proper map

    Phi_ell(Y)=integral_0^Y exp(ell u(Z))dZ.

The same low/high jet algebra applied to exp(ell u) yields
N_(4,10)(exp(ell u)-1-ell u3Y^3)<ell*3e13. Integrating with (9)
therefore gives the explicit forward estimate

    N_(5,11)(Phi_ell-Y-ell u3Y^4/4)<ell*1e14.      (17)

Phi'>0, Phi(0)=0 and Phi'(0)=1. This is conormal C11 on the OPEN
collar, with the previously established C4 boundary extension; an
allowed Y^5 logY term prevents any inference of ordinary C11 at Y=0.
The old low forward-map bounds remain unchanged because Phi itself
has not changed.

To obtain a global forward identification at the needed order,
explicitly replace only the old auxiliary C6 transition by the C11 step

    eta(Y)=1 for Y<=.05, eta(Y)=0 for Y>=.1,
    eta(Y)=sum_(i=0)^11 binom(23,i)s^i(1-s)^(23-i),
    s=(Y-.05)/.05,
    Psi_ell(Y)=Y+eta(Y)[Phi_ell(Y)-Y].             (18)

Its derivative in s is -[23!/(11!)^2]s^11(1-s)^11. On a complex
radius-.05 neighborhood of [0,1], its derivative has modulus below 5
and eta has modulus below 1.25. Each local polynomial branch hence
has |eta^(j)(Y)|<=1.25 j!400^j through eleven; the endpoint jets
match the constant branches through eleven. No global analytic
cutoff is assumed. On the transition, the weighted coefficient sum
is below 1.25 sum_(j=0)^11 .4^j<2.084. This explicitly quantifies
the forward cutoff multiplication through eleven using (17) and
the known mass term: N_(4,11)(Psi_ell-Y)<ell*3e13 on the outer
collar and Psi_ell-Y=0 farther inward. The low derivative is bounded directly:
|eta_Y|<100, |Phi-Y|<1e-10Y and |Phi'-1|<2e-8 give
|Psi'-1|<2.1e-8. Thus Psi is monotone and is identity toward the
cold region. It agrees with the original proper map for Y<=.05.

This higher auxiliary identification preserves the endpoint physical
metric and its boundary data. It is not a modification of t_ell or
sigma_ell. Quantitative C11 inverse-map estimates and the resulting
C10 proper endpoint chords remain separate work. In particular (17)
is not silently used as an inverse or composition bound.

## Verification and next implication

The checker pins the regional source and metric inputs and verifies
exact cancellation, weighted product/integration margins, the mass-mode
curvature and cubic identities, rational reconstruction and full tensor
ledgers, the C11 cutoff and forward-map bounds, and independent
weighted-log integrals and nonlinear finite-map coefficients. Continuum
claims remain ordinary proofs, not machine verification or samples of
the actual charged source. Failed runs, if any, retain source snapshots
and raw logs.

All regional original-field order-ten bounds are now quantitative.
Complete the inverse endpoint pullbacks and then the differentiated
finite-metric quantum, tensor-return and H6 Schur-source estimates.
The sufficient target R0 R2 R4 N6(p)<=1e57 epsilon_c has not passed.
Exact signed compatibility, nonlinear mapping and existence, physical
state/boundary admissibility and stability remain unproved. No S167
acceptance task or physical gate is completed by (1) or (17).

Remaining to-do list: quantitative inverse endpoint maps and C10
proper pullbacks; finite-metric quantum/Schur estimates; positive or
sharper actual transfer; signed compatibility, nonlinear existence,
physical state/boundary and stability; matching, charge/return,
opening, transport, preparation, separate creation and audit.
