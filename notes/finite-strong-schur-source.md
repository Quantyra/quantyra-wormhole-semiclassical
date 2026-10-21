# A strong weighted Schur source at the actual finite comparison

2026-10-21, E003/S167. Ordinary argument before focused verification.

The actual density-correct residual has an ordinary weighted L2 Schur
source, with a numerical bound. This is the initial source of the
coercive moment certificate, not a bound for its three transferred
sources. No signed cancellation, nonlinear solution or stability result
is inferred from this unsigned estimate.

## 1. Result and unchanged source

Use the full retained physical Hessian at each original comparison
g_ell, the actual TT realization, fixed return and physical measure.
Let lambda denote the outer curvature, g=G lambda, t=r_seed/B, W=t^2,
gradbar=lambda^-1/2 grad and Dbar=D/lambda^2. With Jvol=dmu_ell/dmu0,

    (f_tau,f_c)=lambda^-2 Rcal_ell* (Jvol^-1 f_ell),
    F=||W f_tau||_ell+||f_c||_ell,
    p=80pi^2(f_tau-Bbar Dbar^-1 f_c),
    N(p)=||W p||_ell/sqrt(gamma0).

The density conversion precedes the York transpose. Keep the TT datum
f_c, its mass/return response and the full quantum determinant in D.
We prove

    ||W V_T y||_ell <2e4 g^-1/2 ||y||_ell,
    V_T=sqrt(80pi^2) Bbar Dbar^-1/2,
    ||W p||_ell <3e6 F,
    N(p)<2e28 ell^2.                                 (1)

The last inequality uses the actual residual, without projecting it or
assuming its compatibility. At ell=0 that residual is zero. The actual
signed moment remains unevaluated. The older direct moment bound is
sharper than applying Cauchy--Schwarz to (1); the purpose here is strong
weighted source control needed before graph transfers can be estimated.

## 2. Keep the positive second-derivative margin

The [finite disk form](finite-comparison-quantum-disk-tt.md), before
rounding away its positive terms, gives

    <c,Dc> >=.91 E_ell[c]+C/(100pi) sum_j J2[chi_j T_c]
                                      +C/(40pi) J2[S_c+T_c/3].

Here T_c is the actual disk shape, not the scalar TT amplitude. The
same partition identity used in the [seed graph proof](tt-form-and-scalar-fourth-graph.md)
is valid on the finite canonical orbit:

    Delta T=sum_j chi_j Delta(chi_j T)+Q T,
    Q=sum_j |grad chi_j|^2,
    sum_j J2[chi_j T]>=J2[T]/2-integral Q^2 T^2.

On its transition |theta'|<=500 and C J0,outer/E_ell<3.4e-37.
Thus the lost energy is below 1e-28 E_ell. Two-dimensional conformal
covariance, with the actual orbit metric g2, gives exactly

    C J2[T]=||Box_2 T||_L2(dmu_ell)^2/(4pi),
    Box_2=partial_rho^2+A partial_rho,
    A=a'/a, R=r'/r, Box_4=Box_2+2R partial_rho.

Consequently

    <c,Dc>>=.90 E_ell[c]+||Box_2 T_c||^2/(800pi^2).  (2)

This extends to the actual spectral form domain. First use the compact
shape to compact TT approximation of the [finite realization](finite-comparison-tt-realization-and-range.md),
section 4; it converges in the stronger N_ell norm, including J2.
Next approximate in the already established critical operator graph
domain. Applying (2) to differences makes Box_2 T_c Cauchy in physical
L2. The bounded return and distributional closedness identify its limit.
Finally spectral truncations extend (2) to Dom(D^(1/2)). No identification
with an absolute Weyl/model norm or scalar-c logarithmic norm is used.

For c=Dbar^-1/2 y, the retained low-energy estimates therefore give

    ||Box_2 T_c/lambda||<89||y||,
    (||c||,||gradbar c||,||h_c||,||gradbar h_c||,
       ||sqrt(lambda)eta||,||A eta||,||div eta||)
       <(5.3,50,27,443,27,20,170)sqrt(g)||y||.       (3)

These are unweighted physical norms. The constants follow from the .90
gap and the [finite return](finite-comparison-disk-return-form.md) and
[mixed-energy ledger](finite-comparison-weighted-mixed-form.md).
For example c'= (c'+2Rc)-2Rc gives a coefficient below 29, within 50.
The full tensor derivative is bounded by 420/sqrt(.90)<443, and
||A eta||<.713||h||. The source y need only be ordinary physical L2.

## 3. Two tensor derivatives without a seed-only identity

Primes denote actual proper-radius derivatives. Put K=a''/a,
E2=R'+AR and E1=K-E2. The actual low geometry satisfies

    |K|<10.501lambda, |E1|<52lambda, |E2|<22.001lambda,
    0<=AR<12lambda, R^2<5.1lambda, |R'|<35lambda,
    |Rm|<55lambda, |Ric|<15lambda,
    |grad Rm|<1e16 lambda^(3/2), |E1'|<4e16 lambda^(3/2). (4)

For AR use exp(-2w)(A0+w_y)(R0+S_y), retaining A0 S_y<3.45e-6;
no bound on singular A0 alone is used. The first-curvature component
ledger in the mixed-form proof gives the stated gradient reserve.
The radial orthonormal frame is parallel along the radial geodesic.
Expressing E1 through its three sectional curvatures and R^2 bounds
|E1'| by 3|grad Rm|+2|R R'|, which fits (4).

The exact TT and return identities are

    b'=-2A(b+c)-2R(b-c), eta'=A eta+b+c,
    T_c'= -3c'/2+2Rb-E1 eta,
    Box_2 T_c=-3Box_2 c/2
       +(2R'-2AR-4R^2-E1)b+(4R^2-4AR-E1)c
       -(E1'+2AE1)eta.                              (5)

In particular R'=AR-R^2 is NOT imposed at a finite comparison.
Equation (5) and (3)-(4), with ||b||<2.002||c||, imply

    ||Box_4 c/lambda||<100||y||.                      (6)

An explicit majorant is
(2/3){89+[166.4(11)+120.4(5.3)+104(20)]sqrt(g)
                  +(4e16)(27)sqrt(g)}+4.5(50)sqrt(g)<100.
The term A eta, rather than A itself, appears in this calculation.

The exact tensor rough-Laplacian entries are

    (Box h)_rr=-2K(b+c)-2(R'+2R^2-AR)(b-c)-2(A-R)c',
    (Box h)_ss=Box_4 c+2R^2(b-c),
    (Box h)_tt=-(Box h)_rr-2(Box h)_ss.               (7)

The scalar Bochner identity gives
||Hess c||<=||Box_4 c||+4sqrt(lambda)||c'||; its time component
controls A c'. Thus

    ||Box h/lambda||
      <=7||Box_4 c/lambda||+26||gradbar c||+1000||c||<701||y||.

The torsion-free Ricci commutator, [Semmelmann--Weingart, equation (7)](https://arxiv.org/pdf/1708.04775),
gives the same integrated rank-two estimate used in the seed proof:

    ||grad^2 h||^2<=||Box h||^2
       +(||Ric||inf+4||Rm||inf)||grad h||^2
       +4||grad Rm||inf ||h||||grad h||.

Inserting (3)-(4), g<(4e-17)^2, proves

    ||gradbar^2 h_c||<800||y||.                       (8)

The inequalities first hold on the regular core. The approximation used
for (2), applied to differences, extends (5)-(8) by closed weak derivatives
to the actual form inputs. The bolt is kept as a regular disk point.

## 4. Preserve the two outer powers in the local coefficients

The weak mixed estimate weakened several weights to t^-1/2. Its own
component proof retains stronger information. The following rounded
weighted bounds hold at every original comparison:

    W|J|, W|C4| <500lambda,
    W|grad J|, W|grad C4| <3e17 lambda^(3/2),
    W|grad^2 J| <1e25 lambda^2,
    W|Hess R4| <1e20 lambda^2,
    W|grad Box R4| <1e27 lambda^(5/2).              (9)

Here J is tracefree Ricci. For the values use the core differences
|Delta J|,|Delta C4|<.001lambda, and the outer differences
12000lambda/t^3,6000lambda/t^3. The seed values are 7lambda/t^4
and 40lambda/t^3. These give the first row with margin.

For the derivatives use the exact mixed-form component arrays
M=(30,3000,1e6,4e8), P3=1117000, P4=442510000. On the core t<=40,
physical conversion costs 2.5 per derivative. On the outer end the
metric-jet difference is .001 Mj/t^2; the telescoping product reserve
is 1e6. Thus the first-gradient outer budget is
1000(2.5^3)(1e6 P3)/t^2, while the scalar-Hessian budget is
1000(2.5^4)(1e6 P4)/t^2. The corresponding core budgets, multiplied
by W<=1600, also fit (9). The seed outer derivatives are added.

For completeness a full tensor Hessian count avoids inferring grad^2 J
from Hess R4. With the prior connection arrays G0,G1 and curvature
arrays Q0,Q1,Q2, a sufficient component norm is

    |grad^2 J|<=128{Q2+(8M1+36G0)Q1
                      +(4M2+16G1+144M1G0+320G0^2)Q0}.

There are 64 square-root component factors for a rank-six tensor and a
factor two for Ricci contraction. Differentiating the lowering metric,
four tensor slots and then five slots yields the displayed terms.
Its polynomial coefficients in (M4,M1M3,M2^2,M1^2M2,M1^4) are
(384,34560,19968,1241088,6236160), each below 1e7. The core and
outer conversions above are far below the 1e25 reserve in (9).
No scalar Hessian estimate is substituted for a tensor estimate.

The fifth-jet recurrence in the [trace-gradient proof](finite-trace-gradient-and-projection-transversality.md)
gives B5<2.5e20 before physical conversion. The core weighted budget
is 1600(2.5^5)B5; the outer weighted difference is
2e4(2.5^5)B5. Both, including the explicit seed term, are below 1e27.
Only the established five metric derivatives are used.

The full off-shell trace contact also needs its actual weight. Write
Ebar=Etr/lambda and Zbar=(R4+12lambda)/lambda. Then

    W|Zbar|<1.3e-10, W|gradbar Zbar|<.0017,
    W|Ebar|<2e-10, W|gradbar Ebar|<.002.            (10)

On t<=40 the first two bounds follow from 3.1e-12/t and 4e-5/t.
Here is the outer check, retaining the exact curvature remainder rather
than differentiating its value bound. In the notation of the preceding
trace-gradient proof, for t>=40 the old low jets give

    |u|<8e-9/t^2, |nu|,|V|<2e-9/t^2,
    |u_y|<3.2e-7/t, |nu_y|<1.5e-8/t^2,
    |w_y|<3.3e-7/t, |n|<1e-8/t^2,
    |V_y|<4.3e-8/t^2, |V_yy|<2.1e-6/t^2,
    |d|=|k/T-1|<2.01e-9/t^3, |d_y|<4.6e-8/t^3,
    |B|<2.1e-18/t^4, |B_y|<8.7e-17/t^4,
    |L4nu|<6.3e-7, |(L4nu)_y|<500,
    |L4tau|<22/t, |(L4tau)_y|<146/t.

Here n=u+nu, B=k-1-ell V k; the letter d in this display is only
the ratio difference, not the seed polynomial parameter. In the six
terms of the exact remainder (8) of that proof, the resulting bounds
are W|N|<1e-12 and W|N_y|<2e-7. The dominant differentiated term
is 6(1.001)|d|(500), whose W multiple is below 1.51e-7 at forty.
The first exponential derivative uses its double zero:
-2exp(-2w)[(S_y-ell nu_y)+2ell n w_y], with
S_y-ell nu_y=ell d nu_y+(B/T)tau_y. The other terms are smaller.
The seed and linear curvature tails have order at least t^-4 and fit
the displayed reserves. Conversion to actual proper radius costs less
than 1.000001. This proves the first line of (10).

Finally insert (9) in the exact complete trace

    Ebar=(1+11e/3)Zbar+e Box R4/lambda^2
       -11e Zbar^2/72+(7e/12)C4^2/lambda^2
                     +(11e/6)J^2/lambda^2-20epsilon M_em/lambda,
    e=g/(60pi)<g/180.

The Laplacian value costs at most 2e20 e; its gradient at most 1e27 e.
The differentiated curvature squares fit 1e20 e by the retained
unweighted values |C4|<50lambda, |J|<8lambda and the weighted first
gradients in (9); the Z square fits the same reserve by (10).
The magnetic terms use the fixed flux bound.
This proves the second line of (10), including the density/return contact.

## 5. Strong mixed vertex and the actual residual

Use the exact off-shell mixed row

    Bmix c=-(8piG)^-1[Ttrace h_c-div(Etr eta_c)].

Expand Box(J:h) as a tensor product, retaining Box J, its two gradient
terms and J:Box h. The exact Weyl variation contains at most two tensor
derivatives and has principal bound 6|C4||grad^2 h|. The Lichnerowicz,
curvature, magnetic and scalar-curvature terms remain. Equations (9)
give the sufficient full local estimate

    lambda^-1 ||W Ttrace h||
      <=501||h||+e[5000||gradbar^2 h||
                         +6e17||gradbar h||+4e25||h||]. (11)

For the second derivatives, Box h has norm at most 2|grad^2 h|:
2(500)(1+11/6)+(7/12)6(500)<5000. The gradient term comes from
2 grad J.grad h. The value terms include |Box J|<=2|grad^2 J|,
Hess R4 and every algebraic contraction; they fit 4e25. The R4
prefactor and fixed magnetic term fit the reserve from 500 to 501.
The Weyl-invariant determinant remains in D; its compact conformal trace
does not add another nonlocal term to (11).

The full return term is bounded separately by (10):

    lambda^-1 ||W div(Etr eta)||
         <=2e-10||div eta||+.002||sqrt(lambda)eta||. (12)

Since sqrt(80pi^2)/(8pi)=sqrt(5)/2<1.2, (3),(8),(11)-(12)
bound the coefficient of g^-1/2||y|| in ||W V_T y|| by

    1.2{501(27)+(sqrt(g)/180)5000(800)
             +(g/180)[6e17(443)+4e25(27)]
             +2e-10(170)+.002(27)} <2e4.           (13)

This proves the first assertion of (1). Core approximation and the
weighted coefficient bounds make the strong output Cauchy in weighted
L2; its distributional limit agrees with the original weak mixed map.
No new domain or outer condition has been introduced.

Set r=sqrt(80pi^2)Dbar^-1/2 f_c. The critical inverse gives
||r||<150sqrt(g)||f_c||. Then p=80pi^2 f_tau-V_T r and
80pi^2<800 prove ||Wp||<3e6 F. For the actual unprojected residual,
the [finite source conversion](on-shell-finite-range-response.md) and
[matched forcing estimate](projected-low-correction-and-mass.md) give

    F/sqrt(gamma0)<27(2.2e-13)ell^2/g,
    g>(3e-17)^2,
    N(p)<3e6(27)(2.2e-13)ell^2/(3e-17)^2<2e28 ell^2.

This is an actual unsigned starting-source bound. It supplies neither
membership of K^-1p in the certificate core-limit domain nor a bound
after three transfers D0 K^-1. The [sufficient cold certificate](cold-projection-value-upper-bound.md)
still requires its stated residual/terminal estimates and exact signed
compatibility remains a separate equation.

## 6. Verification and next decision

Freeze this argument, checker and source review. Pin the previous
projection-value receipt and transitive inputs. Check the general-warp
return identities and tensor rough Laplacian independently in four
coordinates, keeping R' independent of AR-R^2. Verify curvature counts,
outer remainder budgets, physical form conversion and all final constants.
An omission control must detect the seed-only R' substitution. The
checks do not evaluate the charged kernel or its signed moment; the
continuum estimates and domain extensions remain ordinary proofs.

The next source task begins with the strong bound (1), not merely a weak
distribution. Obtain actual weighted transfer estimates or admissible
trials with the required errors. Five metric derivatives do not justify
three iterations, so higher source control and endpoint approximation
must still be proved. No Lean run is selected.

Remaining to-do list: actual weighted transfer certificate and signed
compatibility; cold/higher metric response, nonlinear preservation and
existence; physical boundary/state and validity checks; stability,
exterior matching, charge/return, opening, transport, preparation,
separate creation and full audit.
