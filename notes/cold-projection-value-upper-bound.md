# Actual cold value of the full fixed projection direction

2026-10-21, E003/S167. Ordinary argument before focused verification.

The original two-component projection direction has a quantitative
pointwise upper bound throughout the cold region. This closes the missing
projection-direction factor in a SUFFICIENT cold-value subtraction test.
It does not evaluate the signed coefficient or certify the actual source
norms required by that test. The projection, action, state and physical
boundary data are unchanged.

## 1. Statement and retained normalization

Use the seed proper coordinate y=sqrt(lambda)rho, t=r_seed/B=1+s,
G0=F0/x, h=t_h-1, gamma0=4pi beta_period B^3 and g=G lambda. Let

    P=sqrt(gamma0) p_cut,
    w_c(s)=s[1+log(.01/s)],
    X_c(v)=sup_(h<=s<=.001) |v(s)|_R2/w_c(s).

This is a value norm of BOTH fixed-density disk-source components on the
entire cold region. It has no spatial derivatives and is not a metric
norm. The following bounds concern the actual full seed adjoint:

    |U chi0(y)|_R2<1.04, y>=.04,
    |P(s)|_R2<.8/sqrt(h), h<=s<=.001,
    3e43<X_c(P)<1e49.                                 (1)

The lower bound is the earlier [annular lower bound](cold-projection-source-lower-bound.md),
since its interval y in [1,2] lies in this cold region. The upper bounds
are proved below; an L2 error alone would not prove them.

The [fixed-cutoff theorem](quantile-cutoff-localization-and-jets.md) and
the lower-bound normalization give

    p_cut=theta t^-4 chi0/(alpha sigma0),
    sigma0=lambda^-1/4 N0, alpha N0>2.9,
    0<=theta<=1, theta=0 for y<=y(a_-), y(a_-)>.04,
    1/(8*2000^10)<=h<1/2000^10.

Here U=sqrt(4pi beta_period r^2 sqrt(F0)) is the PHYSICAL Liouville
multiplier. In particular dmu=U^2 d rho and

    sqrt(h) P=theta b(y) U chi0/(alpha N0),
    b=t^-5 (G0/h^2)^-1/4.                             (2)

The same scalar multiplier applies to both disk components. No period,
volume or lambda factor is dropped.

## 2. Separate the actual adjoint into three controlled corrections

All operators in this section are at the seed. Write T=-Box,
Q^2=T+4lambda, E=60pi/g+11/3, e=E+4 and kbar^2=E-9/4.
The [actual adjoint construction](finite-comparison-adjoint-and-fixed-projection.md)
and [two-channel estimates](on-shell-adjoint-error-and-trace-gradient.md)
give, with n_mix=sqrt(1+a_mix^2)>=1,

    psi=[phi_ref+z_out-z_el]/n_mix,
    z_out=R_PV VJg_R/(lambda e),
    z_el=Q^-2 VJg_R/(lambda e),
    c_psi=-D^-1/2(g_R)_T/[sqrt(80pi^2)n_mix],
    chi0=(psi,0)+Rcal(0,c_psi).                       (3)

Set Gamma_i=lambda^-3/4 ||(g_R)_i|| and A_i=sqrt(g)Gamma_i.
The seed trace gradient is below 3e-19. Therefore, retaining the scalar
channel and its cross-channel feedback,

    xi<30.1(3e-19)+1e-8<1.000001e-8,
    A_T<1.026xi+3e-7(2520)<.000757,
    A_s<4e-7xi+1.000002(2520)<2521.                   (4)

These are full coupled estimates, not estimates for a scalar-only
reference. The spectral sign J has norm at most one.

Define the scaled physical norm and energy

    ||f||_*=lambda^1/4 ||f||_L2(dmu)=||Uf||_L2(dy),
    H(f)^2=4||f||_*^2+||gradbar f||_*^2,
    gradbar=lambda^-1/2 grad.

For z_el the positive scalar form, by duality, gives

    H(z_el)<=sum_i [(a_i sqrt(g)+b_i/2)/(g e)] A_i,
    a=(4e6,4), b=(1000,2e-5).                         (5)

Indeed the physical vertex is div f_i+q_i with
||omega^-1 f_i||<=a_i sqrt(lambda)||input|| and
||omega^-1 q_i||<=b_i lambda g^-1/2||input||. Since omega<=1,
||Q^-1 div||<=1 and ||Q^-1||<=1/(2sqrt(lambda)), division by
lambda e and multiplication by lambda^-1/4 prove (5). These form
bounds apply to weak sources by continuity, without differentiating
the quantum source. With g e>180 and sqrt(g)<4e-17, (4)-(5) give

    H(z_el)<.00225,
    ||z_el||_*<.001125, ||gradbar z_el||_*<.00225.     (6)

For z_out use the physical Liouville kernel estimates
|G_PV|<1.001/k, |D_input G_PV|<3 and
integral omega^2 d rho<200/sqrt(lambda) from the
[mixed-form theorem](finite-comparison-weighted-mixed-form.md).
The input derivative follows from symmetry of the scalar Green kernel;
the real part preserves these upper bounds. Cauchy--Schwarz in the input
coordinate yields

    |U z_out|<=sqrt(200)/e
       sum_i [1.001 b_i/(kbar sqrt(g))+3a_i] Gamma_i
      <(15 sqrt(g)/180)
       sum_i [1.001 b_i/13+3a_i] A_i <2e-13.          (7)

Here kbar sqrt(g)>13. The correction has oscillatory tails and is NOT
assumed to be unweighted L2. Its pointwise kernel bound is used directly.
The source-divergence boundary term is absent on the regular compact
core; weighted continuity extends the estimate to the retained vertices.

For v_T=Rcal(0,c_psi), the spectral gap D>=.90 E_TT and the
[actual disk-return estimates](finite-comparison-disk-return-form.md), as
combined in [the finite response theorem](on-shell-finite-range-response.md),
give

    ||v_T||_*<5 A_T<.003785,
    ||gradbar v_T||_*<87 A_T<.065859.                 (8)

For clarity the physical tensor estimates before scaling are
||h||<=25sqrt(G E_TT/ lambda) and
||gradbar h||<=420sqrt(G E_TT/lambda). The return's derivative is
at most .8||gradbar h||+73||h||, hence at most
2161sqrt(G E_TT/lambda). Inserting (3) and the .90 spectral gap
gives the derivative factor 2161/sqrt(.90*80pi^2)<87.
Its value factor is 25/sqrt(9pi)<5. All estimates retain the returned
mass and shape fields; no higher logarithmic return bound is assumed.

## 3. A pointwise estimate without differentiating the fast reference

For any vector H1 field v on an interval I of length one,

    sup_I |v|^2<=||v||_L2(I)^2+2||v||_L2(I)||v'||_L2(I). (9)

Choose a point with |v|^2 no greater than its average, integrate
(|v|^2)'=2 Re<v,v'>, and apply Cauchy--Schwarz. Approximation proves
the assertion for H1 fields. This explicit interval estimate is standard
one-dimensional Sobolev control; general embedding is discussed in
[Teschl, Theorem 7.14](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe2.pdf),
printed page 196. The charged constants here come from (4)-(8).

For y>=.04 take I=[y,y+1]. On the whole interval the seed gives

    m=partial_y log U=R0+A0/2,
    R0<2.24, A0<=1/y+3.25, 0<=m<17.

Consequently, for any f with finite scaled value and gradient norms,

    ||(Uf)'||_L2(I)<=||gradbar f||_*+17||f||_*.

Using (6),(8) in (9) gives uniformly on y>=.04

    |U z_el|^2<(.001125)^2
                   +2(.001125)[.00225+17(.001125)]<(.0072)^2,
    |U v_T|^2<(.003785)^2
                   +2(.003785)[.065859+17(.003785)]<(.032)^2.

The reference satisfies |U phi_ref|<1.000001. Equations (3),(7) and
n_mix>=1 now imply

    |U chi0|_R2<1.000001+.0072+.032+2e-13<1.04.       (10)

The large fast derivative of phi_ref is never put into a Sobolev bound.
Only the small elliptic and returned TT corrections use (9).

## 4. Whole-cold projection bound and a sufficient certificate

Where theta is nonzero in the cold region, y>.04. The exact cold guards
give s/h>=cosh(sqrt(5.98)y) and G0>5.98(s^2-h^2). Therefore

    G0/h^2>5.98 sinh^2(sqrt(5.98)*.04)
            >=(5.98)^2(.04)^2>(2.1)^-4,
    b<2.1.

Equations (2),(10) prove |P|<2.1(1.04)/(2.9sqrt(h))<.8/sqrt(h).
At the remaining cold points theta=0, so the same bound holds.
The weight w_c is increasing for h<=s<=.001. Also
1+log(.01/h)>70: exp(1)<11/4 and (11/4)^69<2000^10/100 suffice.
It follows that

    X_c(P)<.8/(70 h^(3/2))
           <(.8/70)(23)(2000^15)<1e49.              (11)

This proves (1) on the full cold interval, including the region where the
cutoff is identically zero. It does not bound source derivatives.

For the original finite comparison and fixed-density residual, retain

    f_comp=f_ell-p_cut L_ell^0(f_ell)/sigma_ell,
    beta_ell=L_ell^0(f_ell)/(lambda^2 sqrt(gamma0) sigma_ell),
    lambda^-2(f_ell-f_comp)=beta_ell P.

Thus |beta_ell|<=1e-49 epsilon_c is now SUFFICIENT for
X_c(lambda^-2(f_ell-f_comp))<epsilon_c. The earlier
epsilon_c/(3e43) threshold remains only a necessary condition. No
contradiction is implied by the gap between the two bounds.

In the [coercive certificate](coercive-signed-moment-certificate.md),
keep the actual density-correct source, operator, domain and E>1.125e35.
For epsilon_c>0, the actual normalized bounds

    N(r0)<=2e-49 epsilon_c, N(r1)<=2e-14 epsilon_c,
    N(r2)<=2e21 epsilon_c, N(q3)<=2e56 epsilon_c       (12)

imply |beta_ell|<2.71e-50 epsilon_c and therefore a subtraction with
cold value below .271 epsilon_c. Three EXACT admissible solves can
instead use N(q3)<=1e57 epsilon_c; the cold subtraction is then below
.281 epsilon_c. These sufficient targets are one million times tighter
than the corresponding older necessary-budget targets. They have NOT
been verified for the actual residual. The positivity/resonance caveat,
full TT source, operator errors and endpoint membership all remain.

This controls the subtracted term conditional on (12). To bound f_comp
itself, add the original source bound in this same cold norm. It supplies
neither a cold metric inverse nor a nonlinear solution. A small beta
still does not imply the exact physical equation beta=0.

## 5. Verification and next action

Freeze the ordinary proof, checker and source review. Reconcile the
preceding certificate and all transitive hashes. Verify the new rational
budgets, physical scaling, vector interval inequality, thermal geometry
and sufficient targets. Independent smooth warp fixtures check the
Liouville norm/derivative and projection normalization at two nonunit
lambda values. Controls must detect omission of the Liouville factor
and failure of the old coefficient threshold to ensure the new sufficient
budget. No oscillatory physical source is sampled. These fixtures check
the identities and constants, not the charged quantum kernel or beta.

Next apply the sufficient certificate to the actual source, by weighted
residuals or justified graph estimates. Higher regularity and endpoint
conditions remain unresolved. Do not extend derivatives just to bound
the projection's VALUE: (1) already settles that factor. No Lean run is
selected; continuum estimates remain ordinary proofs.

Remaining to-do list: actual weighted source certificate and signed
compatibility; cold/higher metric response, nonlinear preservation and
existence; physical boundary/state and validity checks; stability,
exterior matching, charge/return, opening, transport, preparation,
separate creation and full audit.
