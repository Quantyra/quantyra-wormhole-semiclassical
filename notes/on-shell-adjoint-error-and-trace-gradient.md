# On-shell adjoint error and the remaining local trace-gradient bound

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary proof before
frozen verification. **The full seed adjoint error is below .007.
At finite comparisons the full error is bounded by one explicit local
trace-gradient norm. Its required finite bound remains unproved.**

Use the actual adjoint, fixed seed-density norm and scalar reference
of [the cold overlap theorem](cold-reference-overlap-and-adjoint-budget.md).
Write e_ell=||chi_ell-Phi_ell||_D. At each original comparison define

    z=sqrt(lambda)rho_ell,
    d_tr=||partial_z(Etr/lambda)||_L2(0,infinity;dz).

This is a trace-gradient norm. The complete trace includes the
anomaly; its off-shell density/gauge return is retained in the mixed
row below. The result is

    e_ell < 6.7[30.1 d_tr+1e-8]+.0064.             (1)

In particular d_tr<=3e-4 would give e_ell<.07, inside the preceding
sufficient .1 budget. At the analytic seed the explicit local trace
gives d_tr<3e-19, hence e0<.007 unconditionally. The finite-comparison
trace VALUE bound does not establish d_tr<=3e-4. No full projection
transversality or nonlinear solve is asserted here.

## 1. The reference equation cancels the leading local trace variation

Keep the full physical action, actual TT realization D, mixed row B,
return Y=eta partial_rho, fixed flux and paired quantum domain. Put
T=-Box, epsilon=G/(60pi), a=1+11epsilon lambda/3 and E=a/epsilon.
The exact trace can be written

    Etr=(a-epsilon T)Z+epsilon H,
    Z=R4+12lambda,
    H=7C4^2/12+11J^2/6-20M_em-11Z^2/72.

For a genuine TT input h, r_dot=-J:h and delta T Z=h:Hess Z.
Consequently

    delta Etr=(a-epsilon T)r_dot
                       +epsilon[-h:Hess Z+delta H]. (2)

The actual scalar reference obeys T phi_ref=E phi_ref. Pairing (2)
with it kills the ENTIRE first term, including its Einstein part.
This is a selfadjoint integration by parts on compact regular TT
inputs, extended by the bounds below. It does not set Etr to zero.
The returned mixed row is still

    Bc=-(delta Etr-div(Etr Y))/(8piG).              (3)

Its return term may now be integrated as

    <grad phi_ref,Etr Y>
       =-<phi_ref,Etr div Y+grad Etr dot Y>.        (4)

The bolt boundary vanishes by regular parity and physical radial
volume. Compact-TT inputs and their outer-fixed return vanish near
infinity. Subsequent energy approximation retains the actual z^4
tails. The finite C5 metric defines grad Etr; its weighted falloff
makes d_tr finite. A numerical bound for it is a separate task.

For reference, one integration by parts in delta C4^2 and delta J^2
using the exact formulas of the [TT inverse](stationary-weighted-tt-inverse.md)
gives the following majorant for the surviving local part of (2):

    epsilon integral [200lambda omega |grad phi_ref||grad h|
      +6e15 lambda^(3/2)omega |phi_ref||grad h|
      +(2e20lambda^2+20.02q_mag^2/B0^4)omega |phi_ref||h|]dmu_ell.

Here B0 is the radius B, q_mag the magnetic charge and omega=t^-1/2.
These are the [proved finite curvature coefficients](finite-comparison-weighted-mixed-form.md).
The 8lambda|h| term and the gradient of J:h in that earlier absolute
bound belonged to (a-epsilon T)r_dot and are removed together by (2).
The Weyl, Ricci, magnetic and Hessian-Z terms remain. No nonlocal
TT response is discarded: it remains inside D and its spectral maps.

## 2. A sharper bound for the reference TT vertex

Let g=V*phi_ref=(g_T,g_s), with
V_T=sqrt(80pi^2) B D^-1/2 and the same local scalar channel V_s.
Set gamma=G lambda and

    Gamma_T=lambda^-3/4 ||g_T||,
    Gamma_s=lambda^-3/4 ||g_s||,
    xi=sqrt(gamma) Gamma_T.

The finite spectral form has E_TT[D^-1/2 y]<=||y||^2/.90. In units
sqrt(G E_TT), the retained bounds for h,grad h,Y are respectively
25/sqrt(lambda), 420 and 25/lambda. Also

    div Y=(b+c)+2(A+R)eta,
    ||div Y||<160 sqrt(G E_TT/lambda).             (5)

For (5), use ||b+c||<=3.002||c||, ||A eta||<.713||h||,
R<2.25sqrt(lambda), ||eta||<4.83sqrt(8piG E_TT)/lambda,
and ||h||<4.8sqrt(8piG E_TT/lambda). This bounds the polar
product A eta, not the singular coefficient A separately.

The scalar reference satisfies

    |U phi_ref|<1.000001, |U grad phi_ref|<3k,
    integral omega^2 d rho_ell<200/sqrt(lambda),
    4e17<kbar=k/sqrt(lambda)<5e17,
    3e-17<sqrt(gamma)<4e-17, sqrt(gamma)kbar<14.     (6)

The derivative bound follows from D0 v_B=-k sqrt(pi k rho/2)J1,
whose modulus is below 2k, the differentiated Bessel resolvent
correction below 2(72000sqrt(lambda))/(1-2e-13), and
|D_rho-D0|<5sqrt(lambda). Thus it applies at the bolt as well.

Multiply the pairing bounds by sqrt(80pi^2)/(8piG) and use
sqrt(5)/(2sqrt(.90))<1.2. The gradient-return term in (4) gives

    sqrt(gamma)lambda^-3/4 |g_T,gradient|
       <=30 lambda^-5/4 ||phi_ref grad Etr||_L2
       <30.1 d_tr.                                (7)

The physical Liouville normalization is essential here:
||phi_ref grad Etr||<1.000001 lambda^(5/4)d_tr.
The remaining terms in (2),(4), using (5),(6), have total xi below
1e-8. The largest of them is the trace VALUE term, bounded by
1.2(160)(3.2e-12)(1.000001)sqrt(200)<9e-9.
The on-shell derivative term costs less than 2e-11; the other local
terms are smaller. Therefore

    xi<30.1 d_tr+1e-8.                             (8)

The local scalar reference vertex must be kept separate. Its gradient
is large at this energy. The established coefficients 4sqrt(lambda)
and 2e-5sqrt(lambda/G), together with (6), imply only

    Gamma_s<180 kbar.                             (9)

Using (9) as a bound for the TT vertex would lose the needed estimate.

## 3. Propagate the two channels through the actual positive elimination

Use A,L,N,M_R and J from the [real adjoint construction](finite-comparison-adjoint-and-fixed-projection.md).
All the following are block operator-norm bounds, not replacements
of operators by scalar matrices. With b0=2.6e-13, positivity of A
and its exact local Schur elimination give the entrywise majorant

    ||N_ij|| <= Nmaj_ij,
    Nmaj=[[1.000001,2.6e-7],[2.6e-7,1.000001]].

For example the off-diagonal bound is sqrt(b0)/(2(1-b0)); the TT
diagonal correction is at most b0/(4(1-b0)). A_TT may be large.
For vertex coefficients f=(4e6,4), g=(1000,2e-5), the actual
zero/one/two-derivative kernel bounds give

    ||L_ij|| <=201 g_i g_j/(180 kbar)
       +(600/180)(f_i g_j+g_i f_j)sqrt(gamma)
       +601 f_i f_j/kbar.

In particular a valid majorant is

    Lmaj=[[.02405,2.405e-8],[2.405e-8,2.405e-14]].

The same bounds hold for Re L. The positive matrix
(I-Nmaj Lmaj)^-1 Nmaj then bounds the two-channel real inverse.
Direct rational inequalities give, for g_R=M_R^-1g,

    Gamma_RT <=1.026 Gamma_T+3e-7 Gamma_s,
    Gamma_Rs <=4e-7 Gamma_T+1.000002 Gamma_s.       (10)

Thus the large local reference derivative enters the TT response
only through a small off-diagonal factor. This is a consequence of
the actual positive block elimination and resolvent estimates.

## 4. From auxiliary channels to the full disk adjoint

The exact TT reconstruction is

    c_psi=-D^-1/2 (g_R)_T/[sqrt(80pi^2)sqrt(1+a_mix^2)].

The critical gap and actual disk return give

    ||Rcal(0,c_psi)||_D <5 sqrt(gamma) Gamma_RT.    (11)

Indeed ||D^-1/2||<=sqrt(8piG/(.90lambda)) and
||Rcal(0,c)||_L2<25||c||: the finite disk return is below 4.6||h||
and its finite TT lift is below 4.8||c||. Conversion to seed volume
costs less than 1.000001, and t^-2<=1. The square-root constant
is 25.1/sqrt(9pi)<5. No higher logarithmic return bound is assumed.

For the scalar correction keep BOTH pieces
R4_PV=(R_PV-Q^-2)/(E+4lambda). The two multiplication weights have
integrals below 32/sqrt(lambda) and 200/sqrt(lambda), so their product
of L2 norms is below 80/sqrt(lambda). The actual kernel bounds and
the spectral inequalities

    ||Q^-2||<=1/(4lambda),
    ||Q^-2 div||<=1/(4sqrt(lambda))

therefore give, channel by channel,

    ||(R4_PV V_T y,0)||_D <1.5sqrt(gamma)lambda^-3/4||y||,
    ||(R4_PV V_s y,0)||_D <3e-8sqrt(gamma)lambda^-3/4||y||. (12)

For clarity, after dividing by sqrt(gamma), the elliptic constants
are at most 1.000001(g_i/4+f_i sqrt(gamma)/4)/180.
The outgoing constants add at most
1.000001[80(1.001)g_i/(13)+240 f_i]sqrt(gamma)/180.
These fit the rounded numbers in (12); the use of 13 follows from
sqrt(gamma)kbar>sqrt(180)>13. The factors of lambda cancel as
displayed, rather than being reset to one in a physical estimate.

Combine (9)-(12), sqrt(gamma)kbar<14 and the normalization error
below .0013 from the preceding overlap theorem. Since
sqrt(gamma)Gamma_s<2520, the result is

    ||chi_ell-Phi_ell||_D <6.7xi+.0064.            (13)

Equations (8),(13) prove (1). The common-density norm, returned mass
component, low-frequency subtraction and both channels are retained.
For d_tr<=3e-4 the right side is below .07. This is a sufficient
condition; it is not a numerical evaluation of the finite d_tr.

## 5. The seed satisfies its full-adjoint budget

At ell=0 use the explicit retained seed, not a finite derivative
inferred from a value bound. Its Z/lambda=(d/x)(4t^-5-6t^-6),
|d|/x<1.3e-22, gives

    |partial_y(Z/lambda)|<1.7e-20/t^5.

For the differentiated Laplacian, F0/x<=5t^2, F0,t/x<=22t and
|F0,tt|/x<=21 imply

    |partial_y(Box Z/lambda^2)|<9e4(|d|/x)t^-5.

To check this bound, the first three t-derivative coefficient sums
of Z/lambda are 56,372,2856 times |d|/x. Differentiating
(F0/x)Z_tt+(2F0/(xt)+F0,t/x)Z_t costs
22(372)+5(2856)+75(56)+32(372)=38568 before the sqrt(5) proper
conversion. The seed Weyl/Ricci derivative bounds are respectively
(40,1000) and (7,100) in their retained lambda units and powers of
t^-1. The magnetic derivative is bounded using R0<2.24 and
Gq_mag^2/(lambda B0^4)<4e-27. In the full trace formula these give

    |partial_y(Etr0/lambda)|<2e-20/t^4.

The same weight integral below 200 gives d_tr,0<3e-19. Thus

    e0=||chi0-Phi0||_D<.007.                       (14)

This proves the seed half of the cold-overlap budget for the FULL
coupled adjoint, including its interior correction and returned TT
field. It is stronger than a bound on the asymptotic phase alone.

## Verification and remaining decision

Freeze the ordinary argument and checker. Reconcile the cold-overlap,
mixed-form, finite-adjoint and seed-outgoing receipts and all hashes.
Check the on-shell trace cancellation and return integration by parts,
the two-channel positive majorants, scalar/TT reconstruction factors,
seed derivative budgets and conditional threshold. Use independent
finite block systems, an exact smooth compact integration fixture and
explicit scalar-polynomial seed derivatives. These tests do not certify
the unproved finite trace-gradient threshold. Preserve all scientific
runs and exact failed sources. Continuum domains and spectral bounds
remain ordinary arguments; no Lean run is selected.

Next prove d_tr<=3e-4 along the original finite comparison, or use
the more flexible bound (1) in the signed-overlap criterion. Its local
trace is known exactly, but differentiating it requires a justified
five-jet estimate. The earlier value bound alone is insufficient.
This local task now replaces an unquantified full-adjoint interior
estimate. If its upper bound is too large, that is not a zero-overlap
or nonexistence result. Pi_cut stays fixed, and every S167 acceptance
task remains open.

Remaining to-do list: finite trace-gradient bound and fixed-projection
transversality; useful full response and cold/higher control; nonlinear
preservation and solve; signed physical compatibility; exterior/state,
charge/return, opening, stability, transport, preparation, creation and audit.
