# Quantifying the flat angular-tower correction at the outer end

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary argument
and verification contract before execution. **Outer model component;
the actual curved weighted inverse remains open.**

The full flat angular tower differs from its logarithmic normal
operator by a bounded third-order Fourier multiplier. After the
physical fourth-power factor and weight 7/2, this difference is
explicitly small on a short outer collar. Its finite fourth contact
is determined, not chosen. This closes the spectrum-discreteness
part of the normal comparison; the actual curved, thermal and
global-domain differences remain separate obligations.

## 1. Fixed operators and theorem

Use the full masses and multiplicities of
[the compact massive proof](compact-massive-logarithmic-graph.md):

    M_n^2=n(n+q)/C, d_n=2(q+2n), n>=1,
    K(p)=C kappa4 p^4/(30pi)
       +sum_n d_n M_n^2/(2pi) integral_0^1 R(a p^2/M_n^2) d alpha,
    a=alpha(1-alpha), R(t)=log(1+t)-t+t^2/2,
    kappa4=H_q/2-gamma_E-log(mu0 B).

This is the translation-invariant, zero-winding GEODESIC flat
multiplier used in the preceding comparison, with every massive
angular level and its original fourth restoration. The actual
Hessian's potential/tadpole contact is not contained in K; it
belongs to the actual-minus-model comparison. K is not the
complete actual thermal/capped response. Define

    ell_o=B/sqrt(C),
    K_infty(p)=C p^4[log(p/(mu0 ell_o))-31/30]/(30pi), p>0,
    K_infty(0)=0, E(p)=K(p)-K_infty(p).

The coefficient -31/30 will be derived below. In particular no
ordinary fourth derivative may be added to E. On the retained box,

    |E(p)| <=1000 C (1+p^2)^(3/2), p>=0.                  (1)

Let E_ev act on the even extension from Y>0, with Fourier multiplier
E(|p|) on the real line. For 0<epsilon<=1/16 put

    R_epsilon v=1_(0<Y<epsilon) (30pi/C) Y^4 E_ev v,
    delta=7/2, x=log Y.

For smooth v compactly supported in (0,epsilon), with zero extension
of Y^-delta v to the whole x-line,

    ||Y^-delta R_epsilon v||_L2(dx)
       <=1e30 epsilon ||Y^-delta v||_H4log(dx).           (2)

The input norm has multiplier (1+p^2)^2 log(2+|p|).
The estimate extends to the stated weighted completion. The cutoff
is a localization of an operator difference, not a replacement
boundary condition for the physical problem.

Replacing ONLY K_infty by K in the flat normal physical density
and b=2T identification gives the shape-entry difference

    lambda^-2 Delta N_b=d0 R_epsilon, d0=1/(480pi^2).

Thus (2) and the accepted complete normal inverse give, at
epsilon=1e-60,

    ||N^-1 Delta N_b|| <8.7e49 epsilon <1e-9.             (3)

This is the composition budget for THIS ONE correction. It is not
the norm of the actual-minus-normal operator and does not establish
a global Neumann inverse or a nonlinear contraction.

## 2. A real-series estimate with its finite constant

Set t>0 and

    F(t)=-2[psi(1+nu_+)+psi(1+nu_-)],
    nu_++nu_-=q, nu_+ nu_-=t.

The [digamma series](https://dlmf.nist.gov/5.7.E6), inspected
2026-10-20, gives the real representation

    F(t)=lim_(N->infinity) [2 sum_(n=1)^N
                      (2n+q)/(n(n+q)+t)-4log N].

It avoids estimating separate complex roots nu_+/- near their
coalescence. Define

    phi_t(s)=(2s+q)/(s(s+q)+t),
    psi_t(s)=(2s+q)/(s(s+q)+t)^2=-partial_t phi_t(s).

Both functions on s>=0 are nonnegative and have at most one
turning point, a maximum. Also

    sup phi_t <=q/t+1/sqrt(t), sup psi_t <=q/t^2+t^-3/2.

Indeed 2s/(s^2+t)<=1/sqrt(t), while the q numerator is bounded
by q/t; divide by the denominator again for psi_t. Their total
variations are at most twice these suprema.

For any such function h, summing the fundamental theorem of calculus
over [n-1,n] gives

    |sum_(n=1)^N h(n)-integral_0^N h(s)ds|
       <=integral_0^N |h'(s)|ds.

Since integral_0^N phi_t=log[N(N+q)+t]-log t and
integral_0^infinity psi_t=1/t, passage to the limit proves

    |F(t)+2log t| <=4q/t+4/sqrt(t),
    |t F'(t)+2| <=4q/t+4/sqrt(t).                       (4)

Differentiation of the convergent positive series is uniform on
closed positive t intervals. The second line is a separate
derivative estimate, not differentiation of the first bound.

Write z=p^2 and k(z)=K(sqrt z). The retained exact identity is

    k''(z)=-C/(2pi) integral_0^1 a^2
                  [F(Caz)+Caz F'(Caz)+4log(mu0 B)]d alpha.

Insert -2log(Caz)-2 from (4). Using

    integral a^2=1/30, integral a^2 log a=-47/900,

gives

    k_infty''(z)=C[log z-2log(mu0 ell_o)-17/30]/(30pi).

Twice differentiating the K_infty defined above gives exactly this
expression. Its constant is therefore -31/30, not the lower-bound
constant -38/15 used in the earlier graph proof.

The sum of the two errors in (4), together with
integral a=1/6 and integral a^(3/2)=3pi/128, yields

    |(k-k_infty)''(z)| <=2q/(3pi z)+(3/32)sqrt(C)/sqrt(z). (5)

This proves the high-frequency comparison quantitatively, including
the ordinary fourth coefficient; an arbitrary fourth polynomial
would contradict (5).

For the two integration constants use the actual low-frequency
model, not an asymptotic prescription. The retained box gives

    C>1, 300<q/C<400, 0<kappa4<100, |log(mu0 ell_o)|<20.

For the upper kappa4 bound use H_q<=1+log q, q=2000^4<2^64,
gamma_E>=0 and |log(mu0B)|<1. For the logarithm bound use
C<7(2000)^3<2^36 and the same bound on log(mu0B).
For the corrected C, one may use C0<6.2(2000)^3 and
C/C0<1.01/.99, whose product is less than 7(2000)^3.
R(t)<=t^3/3 and sum d_n/M_n^4<4C^2/q give K(1)<2C.
The earlier inequality pK'(p)<=6K(p) gives k'(1)<6C.
The explicit K_infty and its derivative at z=1 now bound
|E(1)|+|partial_z E(sqrt z)|_(z=1) by 20C.

Integrating (5) from 1 to z>=1 gives, generously,

    |E(p)| <=100C p^2+2q p^2 log p+sqrt(C)p^3, p>=1.

Use log p<=p and q/C<400. For p<=1 use K(p)<=2C p^4
and p^4 |log p|<=1. These prove (1), including all finite
lower-order contacts. No angular cutoff is taken.

There is an independent normalization check. The Weyl principal
coefficient contains -1/6 in the same units, so

    -31/30-1/6=-6/5.

With u=sqrt(xi/C)Y in the normal geometry, xi=lambda B^2,
omega=-log u, and ell_A=B/sqrt(xi)=lambda^-1/2, the combined
logarithm is

    log|D_Y|-omega-log(mu0 ell_o)
       =log|D_Y|+log Y-log(mu0 ell_A).

This reproduces the constant and logarithmic cancellation of the
accepted full normal symbol. It does not calculate curvature errors.

## 3. A separated kernel estimate for the same difference

Use the exact angular counting function N(M)=sum_(M_n<=M)d_n.
On M_n<=M<M_(n+1), including n=0,

    N(M_n)=2C M_n^2+2n,
    |N(M)-2C M^2| <=4q+2sqrt(C)M.                       (6)

The flat two-dimensional off-diagonal bubble is

    B_M(r)=M^4[K0(Mr)^2-K1(Mr)^2]/(2pi^2).

The continuum measure d(2CM^2) gives
-16C/(15pi^2 r^6), the separated kernel of K_infty.
Its ordinary fourth constant has no off-diagonal contribution.
Stieltjes integration by parts, using (6), bounds the difference
by the integral of (4q+2sqrt(C)M)|partial_M B_M|.
The exact derivative is

    partial_M B_M=M^3[2K0(Mr)^2-K1(Mr)^2]/pi^2.

The needed moments are

    integral t^3 K0(t)^2=1/3, integral t^3 K1(t)^2=2/3,
    integral t^4 K0(t)^2=27pi^2/512,
    integral t^4 K1(t)^2=45pi^2/512.

They follow from the positive
[Bessel integral representation](https://dlmf.nist.gov/10.32.E9),
inspected 2026-10-20, by a beta integral. More generally the expression
is 2^(a-3) Gamma(a/2)^2 Gamma(a/2+nu) Gamma(a/2-nu)/Gamma(a)
for integral t^(a-1)K_nu(t)^2. The boundary terms in the Stieltjes
integration vanish at both zero and infinity. It follows that

    |K_E,2D(r)| <=10[q/r^4+sqrt(C)/r^5].

Integrating the flat time separation gives the stationary kernel

    |K_E,1D(s)| <=20[q/|s|^3+sqrt(C)/|s|^4], s!=0.      (7)

The even half-line operator has direct plus image kernels, with
no additional factor one half. Formula (1) controls its diagonal
extension and all contacts; (7) will be used only for separated
inputs. In particular (7) alone is not used to set local coefficients.

## 4. The complete weighted estimate on the half-line

Let h_j=epsilon 2^(-j-1) and I_j=(h_j,2h_j), j>=0.
Choose eta_j(Y)=rho(Y/h_j), with rho=1 on [1/2,4],
zero outside [1/4,8], and the beta(6,6) transitions between them.
Its first three derivatives are bounded by 1e6.
Every point lies in at most seven supports of eta_j.

For Y in I_j write E_ev v=E_ev(eta_j v)+E_ev[(1-eta_j)v].
At delta=7/2 there is a useful exact cancellation of measures:

    ||Y^-delta Y^4 f||_L2(dY/Y;I_j)=||f||_L2(dY;I_j).  (8)

Put w(x)=Y^-delta v(Y). The derivative identity

    partial_Y^k[Y^delta w(log Y)]
       =Y^(delta-k) product_(m=0)^(k-1)(partial_x+delta-m)w

and Leibniz give, for 0<h_j<=1,

    ||(eta_j v)_even||_H3(dY)
        <=1e15 h_j ||w||_(H3(dx);supp eta_j).             (9)

Here the local H3 norm is the sum of derivative squares with
binomial coefficients. To see an explicit reserve, all the polynomial
coefficient sums through k=3 are less than 40; the cutoff derivatives
are below 1e6; Y/h_j lies in [1/4,8]. Each differentiated term thus
costs at most 2^3(1e6)8^4(40)h_j in L2(dx). The binomial H3
sum and the even-extension factor fit inside 1e15.

Equations (1), (8)-(9), disjointness of I_j and the support overlap
bound make the NEAR part, after multiplication by 30pi/C,
smaller than 1e22 epsilon ||w||_H3(dx).

For the far part, (1-eta_j(y)) vanishes unless y<Y/2 or y>2Y.
Both direct and image kernels in (7) then have combined bound

    400[q/max(Y,y)^3+sqrt(C)/max(Y,y)^4].

After multiplying by (30pi/C)Y^4, conjugating the weight and
using dy/y, the kernel is bounded by

    12000pi[(q/C)epsilon^2+epsilon/sqrt(C)]
       times (Y/y)^(1/2) if Y<y/2,
             (y/Y)^(9/2) if y<Y/2.

The two exponential integrals in log coordinates have total at
most 2+2/9<3 in either Schur row. With q/C<400 and C>1 this
gives a far norm below 1e8 epsilon ||w||_L2(dx).
Finally ||w||_H3<2||w||_H4log. These bounds prove (2).
The proof sums infinitely many outer annuli; it is not a check on
a finite set of sample scales.

## 5. Physical normal identification and remaining actual errors

In the normal density, the b=2T massive entry divided by lambda^2
is Y^4 K_ev/(16pi C). Thus its difference is exactly d0 times
R_epsilon. The complete normal inverse in
[S163](weighted-normal-inverse-scale.md) is less than
8.7e19/lambda^2 in the same fourth/fourth-log normalization.
This proves (3). It does not replace that inverse by a shape-only
bound or vary the physical parameters.

Still unbounded quantitatively are the actual-minus-flat geometric
response, thermal winding and global-domain terms, their finite
contacts, and the outer cutoff/compact-annulus matching. These must
be added before any claim about the actual outer inverse. The
compact graph and the mass-functional estimate retain their scopes.

## 6. Frozen verification

Freeze this proof and its checker before execution. Pin the complete
flat-model and normal-inverse receipts and inputs. Verify the rational
series denominators and separate variation bounds, exact integration
of the logarithm and its -31/30 constant, low-frequency integration
constants, all-mass counting inequalities, Bessel derivative and
moments, continuum and stationary kernels, dyadic cutoff derivatives,
weight/measure and derivative identities, infinite overlap/Schur
bounds, and the physical normal composition budget.
Preserve every run and any failed source snapshots.

Finite algebraic and numerical enclosures support the ordinary
infinite-series, distribution and operator arguments; they do not
formally certify the continuum theorem or the unproved actual
curved/thermal problem. No wormhole solution or stability follows.

Remaining to-do list: actual curved and thermal outer response,
finite contacts and annulus matching; higher York maps and full
inverse/residual composition; nonlinear signed compatibility;
exterior/state, opening, stability, transport, preparation and creation.
