# Canonical cold response on the retained comparison path

2026-10-20, E003/S156 supporting S139/S074. Ordinary proof and
verification protocol before acceptance.

The actual retained comparison path has canonical cold first and
second responses, uniformly INCLUDING its endpoints. Its null
response has an integrable quarter-power bolt modulus. The cold
canonical DEVIATION contribution to the constraint Hessian is
therefore finite. Its allowance is coarse and is not an existence
or contraction estimate.

This removes the actual-path admissibility condition left open in
[S153](cold-response-from-real-endpoints.md) for this one path.
It does not certify an arbitrary solver family. It uses neither
real metric estimates outside the known parameter interval nor
parameter differentiation of moving proper coordinates.

## 1. Two different source bounds on one fixed path

Keep the fixed common disk, actual paired domain, physical clock,
scheme and parameters of [S155](comparison-path-interior-response.md):

    sigma_lambda=sigma0+lambda u,
    r_lambda=B+(r0-B)exp(lambda V), 0<=lambda<=1,
    P5_40(u)<1.3e-4.                               (1)

The cold observation region lies within S155's cover Y0>=.051.
The metric estimate's finite-endpoint proof applies for EACH
original parameter lambda in [0,1]. Its once-made endpoint
identification is identity throughout the cold region. The real
canonical endpoint proof consequently evaluates the same scalar
there as the common-chart family (1). The auxiliary proper chord
used to prove that endpoint difference need not be the original
path between its endpoints.

Use [the retained weights](cold-response-from-real-endpoints.md):

    eta_G<8.5e-33, exp(-80)<h<1e-30,
    L_h=1+log(.01/h)<81,
    a_x=s[1+log(.01/s)], h<=s<=.001,
    a_x>=h L_h>exp(-80), a_x<.00331,
    a_x<=2h L_h on chi<=1.

For S_lambda=(2G/(B^2 C))Theta_lambda and f_lambda=S_lambda-S_0,
the REAL endpoint theorem gives

    |f_lambda(x)|<=1.7e-15 lambda a_x,
    |f_lambda(chi)-f_lambda(0)|
        <=1.7e-14 lambda h L_h sqrt(chi), chi<=1. (2)

The proof uses the full angular tower and original actual domain.
No endpoint estimate has been relabeled as a second derivative.

Separately, the local remainder, contact and actual-domain estimates
in S154/S155 are uniform when the complex field's low and divided
high norms are at most .011. Their source bound depends on that
coefficient condition, not on the bidisk parameterization used to
obtain it. Keep the coordinates and cap fixed and put x=2lambda-1.
On the closed Bernstein ellipse E_rho, rho=256,

    lambda=1/2+(w+w^-1)/4, |w|=256,
    |lambda|<=.5+(256+1/256)/4=64.5009765625,
    P5_40(lambda u)<.008386<.011.                  (3)

The interior is also within this coefficient cone. The fixed-domain
inverse and uniformly convergent subtracted angular sum therefore
give parameter holomorphy there, with S155's nonconstant bound
K40<1.25e-11. Hence both f_lambda(x) and
g_lambda(chi)=f_lambda(chi)-f_lambda(0) have modulus at most

    B_c=5e-11                                     (4)

on that ellipse. Two source values suffice for f, four for g.
The fixed vacuum scalar cancels in these differences; it has not
been assigned zero value. No complex spatial coordinate or complex
extension of the real weighted estimate (2) is used.

## 2. Analytic smallness on a CLOSED real interval

Here is a self-contained interval lemma. Suppose F is holomorphic
on a neighborhood of the closed ellipse E_rho, |F|<=B there, and
|F|<=a on [-1,1], with 0<a<=B and log(rho)>4. Define

    k=max(1, log(B/a)/log(rho)).

On the ellipse minus the slit [-1,1], let

    p(z)=log|w(z)|,
    z=(w+w^-1)/2, |w|>1.

The inverse branch is w=z+sqrt(z^2-1), with sqrt asymptotic to z.
Its log modulus is single-valued and harmonic on the slit
complement. It is zero continuously on the slit and log(rho)
on the ellipse. This is the real part of the inverse hyperbolic
cosine; its logarithmic branch is recorded in
[DLMF 4.37](https://dlmf.nist.gov/4.37), inspected 2026-10-20.
The maximum principle for log|F| gives

    |F(z)|<=a exp[k p(z)].                         (5)

Zeros are handled by the subharmonic maximum principle. The
endpoints cause no singular boundary data because p tends to zero.
If only an open-ellipse bound is available, use smaller ellipses
and pass to the limit; (3) in fact has strict coefficient slack.

For |z-x|<=r, x in [-1,1], the ellipse-distance identity and
triangle inequality imply

    2cosh p(z)=|z-1|+|z+1|<=2+2r,
    p(z)<=sqrt(2r).                                (6)

For j=1,2 choose r_j=2j^2/k^2. Then p<=2j/k<=4<log(rho),
so the full Cauchy circle stays inside the ellipse, including
when x=1 or -1. Equations (5)-(6) give

    |F'(x)|<=a e^2 k^2/2,
    |F''(x)|<=a e^4 k^4/32.

For f(lambda)=F(2lambda-1) this becomes

    |f'(lambda)|<=a e^2 k^2,
    |f''(lambda)|<=a e^4 k^4/8, 0<=lambda<=1.      (7)

If a=0 the identity theorem gives zero. The derivatives at the
two interval endpoints are those of the unweighted holomorphic
family, bounded using only the real values within the interval.

The fourth logarithmic power in the second bound cannot in general
be replaced by a uniform second power at endpoints. Let

    f_n(lambda)=a_n T_n(2lambda-1),
    a_n=2B/(rho^n+rho^-n).

The Joukowski identity T_n((w+w^-1)/2)=(w^n+w^-n)/2 proves the
ellipse bound B and real bound a_n. Direct differentiation gives
f_n''(1)/a_n=4n^2(n^2-1)/3, whereas log(B/a_n) grows linearly
in n. These are manufactured analytic functions, not quantum
sources. They identify why a centered half-disk estimate alone
would leave a gap at lambda=0,1.

## 3. Cold scalar responses and their target space

For f in (2), take a=1.7e-15 a_x. Then
log(B_c/a)<91 and k<=91/log(256). For g use
a=1.7e-14 h L_h sqrt(chi), giving
k<=[91+.5log(1/chi)]/log(256). Equation (7) yields

    |partial_lambda S_lambda(x)|<3.5e-12 a_x,
    |partial_lambda^2 S_lambda(x)|<9e-10 a_x,     (8)
    |partial_lambda S_lambda(chi)
                      -partial_lambda S_lambda(0)|
        <3.5e-11 h L_h chi^(1/4),
    |partial_lambda^2 S_lambda(chi)
                      -partial_lambda^2 S_lambda(0)|
        <9e-9 h L_h chi^(1/4).                   (9)

Indeed for t>=0 and p=2,4,
exp(-t/4)(91+t/2)^p<=91^p: its logarithmic derivative
is -1/4+p/[2(91+t/2)]<0. The sharper estimates before
this simplification are sqrt(chi) times a quadratic or quartic
logarithm. Dividing by chi^(1/4) makes their tails tend uniformly
to zero at the bolt. Away from the bolt, uniform unweighted
holomorphy gives continuous parameter derivatives. Splitting
at chi=delta proves C2 into the cold weighted space with
quarter-power oscillation, with one-sided norm derivatives at
the endpoints. Integration along [0,1] proves the weighted
Taylor theorem there. No weighted holomorphy outside the real
interval or C2 into the stronger square-root space is asserted.

The unweighted S155 bounds remain available simultaneously.
The larger coefficients in (8) supply a different, essential
cold weight; they do not improve the unweighted interior bound.

## 4. HH pressure, bolt cancellation and cold constraints

Keep the frozen homogeneous scalar S_* at the same B,A,mu0:

    Q=S-S_*, R=P-S_*/2=F_lambda[Q],
    F_lambda[f](c)=[2W_lambda(c)]^-1 integral_0^c W_lambda' f,
    rho_c=R-Q, n_c=2R-Q.                           (10)

These are normalized canonical deviation profiles. The homogeneous
canonical density remains allocated to the full source and is NOT
included in the deviation estimate below. The null profile n_c
equals the full canonical null profile.

S153's fixed-seed value calculation plus (2) gives
|Q|<3e-7 a_x and |S(chi)-S(0)|<5e-4 h L_h chi^(1/4).
Use S155's already proved bounds on the SAME actual path:

    ||F_lambda||<=1.2,
    ||partial_lambda F_lambda||<.000253,
    ||partial_lambda^2 F_lambda||<2.5e-11.         (11)

All operator derivatives annihilate spatial constants. Apply the
full product rule to Q for the ordinary cold bound and to
S-S(0) for the bolt bound. If (z0,z1,z2) bounds that source,
the corresponding density/null profile bounds are

    3.4 z0,
    3.4 z1+2(.000253)z0,
    3.4 z2+4(.000253)z1+2(2.5e-11)z0.            (12)

The source weights are monotone on [0,c], so they can be taken
outside these operator bounds. This includes both changing-kernel
cross terms without requiring a varying spatial source gradient.
The resulting rounded table, uniform on 0<=lambda<=1, is

| Profile weight | Value | First lambda response | Second lambda response |
| --- | ---: | ---: | ---: |
| max of absolute rho_c,n_c / a_x | 1.1e-6 | 1.7e-10 | 3.1e-9 |
| absolute n_c / [h L_h chi^(1/4)], chi<=1 | .0017 | 2.6e-7 | 3.1e-8 |

The quarter-power bolt integral uses S153's explicit variant
of [S152](fixed-chart-global-constraint-transfer.md):
3L_h integral_0^1 chi^(-3/4)dchi=12L_h<972.
The remaining regional constants are unchanged. The actual
radius/shape direction size is d_g=3e-7. Thus define

    m2=1.003(3.1e-9)+4.02 d_g(1.7e-10)+5.54 d_g^2(1.1e-6),
    c2=1.006(3.1e-9)+10.08 d_g(1.7e-10)+26.5 d_g^2(1.1e-6),
    b2=1.006(3.1e-8)+10.08 d_g(2.6e-7)+26.5 d_g^2(.0017).

Normalization cancels one factor 2eta_G as in S153. For the
COLD canonical DEVIATION alone the second-response bounds are

    fixed-seed zero-datum mass primitive: 46m2/2 <7.2e-8,
    absolute mass density in BC units: .000002m2/2 <3.2e-15,
    absolute lapse density L1: [13122c2+972b2]/2 <3.6e-5. (13)

These are derivatives of maps into absolute-density spaces;
no derivative of a scalar absolute value is taken. The Taylor
allowances about zero are half (13) times lambda^2.

In particular the lapse allowance at lambda=1 is 1.8e-5.
It is finite but far above the retained ~1e-8 conformal displacement.
Those are different norms; the comparison signals that this bound
does not certify a small correction. No inverse, solver radius or
global physical error is inferred. S153's conditional mixed-family
theorem still has its own hypotheses; this result supplies only the
actual one-dimensional path and its endpoint response.

## 5. Predeclared verification and remaining work

Before acceptance, revalidate S155, S153 and S152 receipts and
their current input hashes. Check the ellipse/cone margins,
Green-function ellipse identity, Cauchy radii, parameter factors,
logarithmic weights, every Ward cross term and regional budget.
The maximum principle and infinite angular-source arguments are
ordinary proofs, not numerical claims.

Use 70 digits. Test Chebyshev fixtures of degrees 2,8,16,32 at
lambda=0,.3,1. Compare coefficient-polynomial derivatives with
64/128-node Cauchy contours of lambda radius .01, tolerance
1e-35(1+|target|), normalized by the real amplitude. Compare
first/second real differences at steps 1e-7,5e-8 with tolerance
1e-7(1+|target|), refinement ratio .35 unless both errors are
below 1e-40. At endpoints use second-order ONE-SIDED differences
entirely inside [0,1]. Check the endpoint n^4 identity and reject
a degree-independent n^2 second-response bound.

Independently integrate a manufactured Ward family with
W_lambda(c)=c^2 exp(2lambda gamma c^2), gamma=.03,
Q_lambda(c)=hL_h sqrt(c)(1+.2lambda+.1lambda^2), h=exp(-79).
Use lambda=0,.4,1 and cuts exp(-4),.1,1. Substitute c=y^2
in the pressure integral to remove the square-root quadrature
endpoint. Normalize by hL_h sqrt(c). Compare 32/64 Gauss nodes
and independent adaptive quadrature of the full second-response
formula to 1e-25. Use real steps 1e-4,5e-5 with the same
difference tolerance/refinement policy and endpoint one-sided
differences. Check that freezing the Ward kernel changes the
result. At gamma=0 also check the exact null profile and its
finite absolute bolt lapse integral. The deliberately larger
fixture gamma makes cross terms detectable; it is not the
retained metric or a continuum quantum determinant.

Preserve every run, raw output, failed source snapshot and
input identity before/after. Do not rewrite historical evidence.

Remaining to-do list: outer recombined and other action-sector
responses, including homogeneous allocations and the E0 residual;
a general admissible family and usable full higher-order inverse;
coupled exterior/state, charge, return and signed-source solve;
opening, stability, transport, preparation and separate creation.
