# Physical TT contacts and the fixed-density block at finite comparisons

2026-10-20, E003/S167 supporting S139/S074. Ordinary argument before
frozen verification. **The physical exponential and affine tensor TT
forms are positive. The fixed-volume returned disk block is identified;
mixed blocks, full operator domains/range and the nonlinear goal remain
open. No dynamical stability result is asserted.**

Let g_ell be the original finite comparison, 0<=ell<=1, and let
h=diag(b,-b-2c,c,c) be a genuine spherical TT input on its smooth,
bolt-regular, compact-outer core. For the full retained action I,

    Q_I,exp[h] >= .97 E_ell[c],
    H_I,aff[h,h] >= .97 E_ell[c].                    (1)

Here Q_exp differentiates g exp(s g^-1 h); H_aff differentiates
g+s h. Both are restricted TT forms, not assertions that an
unprojected tensor operator preserves TT. The result also holds for
I^(alpha)=I0+alpha Gamma_rest, 0<=alpha<=1, at every g_ell.
Thus it includes the actual bookkeeping path alpha=ell and the full
physical theory alpha=1, without identifying these two parameters.

For the original actual disk return R_ell c=(S_c,T_c), the fixed-volume
disk block already satisfies

    <c,D_ell c>_muell=I_disk''[R_ell c,R_ell c]>=.90 E_ell[c]. (2)

No extra physical tensor contact is needed to deduce (2). Its
definition, as in the [seed TT inverse](stationary-weighted-tt-inverse.md),
is the returned DISK Hessian itself. The prior disk theorem therefore
already supplied this form lower bound. The transfer in (1) establishes
the separate physical tensor statements. Neither (1) nor (2) identifies
the full operator realization or supplies the candidate coupled inverse.

## 1. A uniform bound for the complete remainder gradient

Keep the original action split I0=S_classical+Gamma_ref and
Gamma_rest=Gamma_Dirac-Gamma_ref, the prepared period, clock,
physical flux, state, paired realization and common coefficients.
The [complete source theorem](comparison-path-complete-source-response.md)
and its four regional components bound the actual value of
Gamma_rest's source on g_ell. No parameter or spatial derivative
of a varying source is used here.

Write t0=r0/B, t=r_ell/B, x=lambda B^2, C/x<6.1,
lambda=lambda_e, s=t0-1 and L_e(t0)=1+log t0. The normalized profiles are

    rho=2G r_ell^2 rho2/(B^2 C),
    N=2G r_ell^2(rho2+p2)/(B^2 C),
    A=8piG r_ell^4 pS/(B^2 C).

The symbol A in this section is the angular source, not the orbit
logarithmic derivative used later. Summing the already allocated
canonical, Weyl and LLL/reference values gives

| Fixed seed region | max(abs(rho),abs(N)) | abs(A) |
| --- | ---: | ---: |
| cold, s<=.001 | 1.2e-6 a_x | 1.7e-6 a_x |
| middle, 1.001<=t0<=40 | 1.8e-9 | 2.6e-9 |
| outer, t0>=40 | 3.2e-7 L_e^2 | 4.7e-7 L_e^2 |

Here a_x=s[1+log(.01/s)]<=.00331. For the cold orbit row, the
three contributions are 1.1e-6 a_x, 1e-15 a_x and 1e-23 a_x.
For the middle they are 1.7e-9, 1e-16 and 1e-18. At the outer
end they are 3.1e-7 L_e^2, 1e-9 L_e^2 and 1e-21 L_e^2.
The angular column is exactly the retained complete radius-source
ledger. It is not reconstructed from an assumed solution of the
orbit constraints.

In particular the outer row uses the actual correlated cubic and
vacuum cancellation of
[the outer recombination](comparison-path-outer-recombination.md).
The singular canonical term and the Weyl mass term are combined
before this norm; their pressure memory and fixed reference vacuum
are retained. Applying the bound to either singular sector separately
would be invalid. The cold row similarly combines the product source
with the frozen reference before applying its vanishing weight.

Define the action gradient by

    Gamma_rest'[v]=integral mathcalE_rest:v dmu_ell.

In the orthonormal frame its entries are, up to Euclidean sign,
one half of (-rho4,p_r,pS,pS), where rho4=rho2/(4pi r_ell^2).
Consequently, if the orbit and angular row bounds are R and A_b,

    (8piG/lambda)||mathcalE_rest||
      <=(C/x)/(2t^4) sqrt(5R^2+2A_b^2).              (3)

The factor one half belongs to action variation; the 4pi sphere
factor belongs to converting rho2 to rho4. Neither is omitted.
For the cold row use t>=1 and a_x<=.00331; (3) is below 4e-8.
For the middle it is below 2e-8. For the outer row, |log(t/t0)|<3e-8
gives t^-4<1.001t0^-4. The function (1+log t0)^2/t0^4 decreases
for t0>=1, and L_e(40)<5, so (3) is below 3e-11. Hence globally

    (8piG/lambda)||mathcalE_rest||_infinity <4e-8.    (4)

The [classical gradient bound](finite-comparison-classical-disk-contact.md)
has the same action normalization and is below 1.8e-5. Therefore
for I^(alpha), uniformly in alpha and ell,

    (8piG/lambda)||mathcalE_Ialpha||_infinity<1.81e-5. (5)

These are full four-component Frobenius norms. The outer physical
weight makes them bounded even though the normalized source ledger
contains logarithms. They are not bounds for a general nonlinear
neighborhood beyond the specified comparison path.

## 2. Exact gauge and chart contact, with its cancellation

Use the actual outer-fixed return vector Y=eta partial_rho and
write A_orb=a'/a, R=r'/r. Then

    eta'-A_orb eta=b+c, k=h-L_Yg,
    S_c=c/2-R eta, T_c=-(b+3c)/2-(A_orb-R)eta.

For any of the retained diffeomorphism-invariant functionals J,
the affine tensor identity and exponential accelerations give

    J_disk''[R_ell c,R_ell c]-Q_J,exp[h]=J'[Z],
    Z=L_Y(2h-L_Yg)+k^2-h^2.                         (6)

This holds off shell. No gauge direction is assumed to annihilate
the Hessian. Set s_i=(eta',A_orb eta,R eta,R eta); the relative
components of L_Yg are 2s_i. Since
(L_Yv)_i=eta v_i'+2s_i v_i for a diagonal covariant tensor,
the exact cancellation is

    Z_i=2eta(h_i'-s_i').                            (7)

The previously proved actual return estimates imply, with
e=8piG E_ell and Cc=||c||_muell,

    ||h||<4.8 Cc, ||h'||<75sqrt(e),
    ||eta||<4.83 Cc/sqrt(lambda), ||s'||<215sqrt(e),
    integral |Z| dmu_ell<3000 e/lambda,
    integral |h^2| dmu_ell<=||h||^2<24 Cc^2.         (8)

All singular polar products are grouped as in the classical contact
proof. Z=O(rho^2) at the bolt, while the outer-fixed return vanishes
outside the support of the core input. The vector has regular odd
parity. Thus no bolt displacement, extra entropy contact or omitted
outer boundary functional enters (6). The retained quantum state and
paired boundary domain are carried by the same admissible diffeomorphism.

Apply (6) to Gamma_rest only. Equations (4),(8) give

    |Gamma_rest'[Z]|<.00012 E_ell.                  (9)

This avoids spending the larger classical-gradient contact a second
time. The [classical exponential form](finite-comparison-classical-tt-form.md)
already has Q_I0,exp>=.99 E_ell. The
[complete quantum disk proof](finite-comparison-quantum-disk-tt.md)
gives (Gamma_hat,mass+W4-Gamma_ref),disk''>=-.01 E_ell,
and the canonical LLL costs at most 8e-24 E_ell. Thus

    Gamma_rest,disk''>=-(.01+8e-24)E_ell,
    Q_Ialpha,exp> [.99-alpha(.01+8e-24+.00012)]E_ell. (10)

The resulting lower margin exceeds .9798. Finally

    H_Ialpha,aff[h,h]=Q_Ialpha,exp[h]-Ialpha'[h^2],
    |Ialpha'[h^2]|<24(1.81e-5)E_ell=.0004344 E_ell. (11)

Both bounds in (1) follow with room to spare. Only the source VALUE
was used to control these contacts, not a derivative of the stress.

## 3. Fixed density, the actual disk block and the moving-volume map

The [general-background Newton identification](general-background-newton-operator.md)
uses disk coordinates (Phi,V), tangent u=(S,T), and gradient F per
moving physical four-volume. Let l(u)=4S+2T. At each fixed endpoint,

    A_disk u=DF[u]+F l(u),
    mathfrakF=Jvol F, D mathfrakF[u]=Jvol A_disk u,
    Jvol=dmu_ell/dmu0.                               (12)

These are exact identities, not additional model terms. With the
correct fixed-volume dual transpose,

    <R_ell d,D mathfrakF[R_ell c]>_mu0
       =I_disk''[R_ell d,R_ell c]
       =<d,D_ell c>_muell.                          (13)

Thus (2) is already the fixed-density TT-form bound required by the
definition D_ell=R_ell^* A_disk R_ell. The previous milestone's
unresolved physical-tensor contacts did not obstruct this algebraic
identification. A sharp full logarithmic map for R_ell and an actual
operator realization are separate requirements, still open.

For comparison, the density term itself has a useful quantitative
bound. The stronger orbit combination cancels R eta exactly:

    S_c+T_c=-(b+2c)/2-A_orb eta,
    ||S_c+T_c||<5.1 Cc, ||S_c||<12 Cc.

Here ||A_orb eta||<1.008||b+c|| and ||b+c||<3.002 Cc.
The relative returned tensor is
k=diag(2(S_c+T_c),2(S_c+T_c),2S_c,2S_c), so

    ||k||^2<1400 Cc^2, ||l(R_ell c)||<34.2 Cc,
    F dot R_ell c=mathcalE_Ialpha:k.

For two TT inputs c,d, Cauchy and (5) give the energy-dual bound

    |integral l(R_ell c) F dot R_ell d dmu_ell|
        <.025 sqrt(E_ell[c] E_ell[d]).              (14)

Indeed 34.2sqrt(1400)<1300 and 1300(1.81e-5)<.025.
Consequently the returned moving-volume derivative has quadratic
lower bound .87 E_ell, using (2). It remains generally nonsymmetric:
(F tensor l)^*=l^T tensor F^T. Do not identify its adjoint moment
with that of A_disk. For the intended fixed-density Newton operator
use (12)-(13), including the density term, and retain the .90 bound.
Do not substitute the .97 tensor bound for this differently defined D.

## 4. Form scope and verification

The contacts (9), (11) and (14) are bounded in the actual TT energy.
They therefore extend continuously to the explicit stronger core
completion in the quantum disk theorem. The physical tensor forms
inherit that completion and their lower inequalities there. This
does not identify the mass/fast operator domains, establish all outer
approximations, or imply that an arbitrary semibounded form equals
the retained nonlocal operator. No coupled kernel or source-range
theorem follows merely from these positive restricted forms.

The frozen checker reconciles the complete disk, complete source,
general Newton and seed block definitions with their exact inputs.
It verifies source assembly and physical normalization, all regional
and contact budgets, the gauge/chart cancellation, exact fixed-density
chain rule and output pairing, and the sharper orbit-return combination.
Omission controls detect the half-stress factor, a missing contact,
and the false off-shell symmetry of the moving-volume derivative.
Continuum source/domain arguments remain ordinary proofs. This is
not formal verification, a numerical solution or dynamical stability.

Pi_cut remains fixed at the seed. Every S167 acceptance task remains
open. The next task is the finite-background mixed block and full
logarithmic realization/range, followed by the remaining nonlinear
regularity, solve and signed compatibility obligations.

Remaining to-do list: mixed blocks; physical-direction membership and
full logarithmic domains/coupled range; fixed projection and cold/higher
response; preserved nonlinear family and projected solve; signed
compatibility; exterior/state, charge/return, opening, stability,
transport, preparation, creation and full-result audit.
