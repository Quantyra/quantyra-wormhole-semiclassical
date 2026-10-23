# The complete current-action product metric inverse

2026-10-22, E003/S169. Ordinary argument before supporting verification.
Progress; the full stable semiclassical wormhole objective is incomplete.

The exact homogeneous product of the CURRENT c_R=10^22 theory has a
complete stationary spherical metric Hessian inverse on complete-H2 L2,
with dimensionless norm <.51 in the physical pairing defined below.
This includes the full Dirac angular tower and the restored curvature
scalar. The two cold roots also remain close: the non-Ward root obeys
|p_root|<4e-8 and |nu-1|<1.34e-8. These are product conclusions, not a
finite-geometry inverse, cold C5 admission bound or retarded stability.

The old 16000 product inverse concerns a REDUCED zero-added-c action.
It is not an inverse of the present theory. Rather than transplant it,
we derive the new full matrix and retain an exact determinant cancellation.
Its norm is not a replacement for the still-unknown K_m or K_c.

## 1. Product, quantum input and its checked hypotheses

Use the exact product from [the current product enclosure](recentered-product-and-cold-residual.md):
X=A^-2, Y=B^-2, R=2(Y-X), C=X/Y, with the same G, q=2000^4,
g_R=2000^-4, lambda0=2000^-10/G, c_R=10^22 and mu0=1/b_classical.
Do not replace it by the nearby RN reference or the finite thermal seed.
The homogeneous pressures obey exactly

    E_A=-Y-3lambda0+piG q^2 Y^2/g_R^2
                    +32piG c_R(X^2-Y^2)-8piG p_A=0,
    E_S= X-3lambda0-piG q^2 Y^2/g_R^2
                    -32piG c_R(X^2-Y^2)-8piG p_S=0.       (1)

In particular tau_D=2(p_A+p_S)=(X-Y-6lambda0)/(4piG).
At the product the Box R anomaly vanishes. The exact product trace is

    tau_*=[6X^2+10XY+(6-30q^2)Y^2]/(1440pi^2).          (2)

We work on complete H2_A with its decaying massive inverse and the
regular origin. This is the invariant product realization used by the
retained bubble proof, not the actual finite outer reflecting domain.
For stationary spherical fields S,T write

    g=exp(2S)[A^2 exp(2T)g_H2+B^2 dOmega^2],
    ell=-Delta_H2, p=ell+2, gamma0=X/(8piG).

The pairing is physical four-volume, 4pi B^2 A^2 dvol_H2, with the
ordinary Euclidean norm of the two components S,T. On complete H2 the
L2 spectrum is ell>=1/4, hence p>=9/4. The source in the inverse below
is the physical variational source divided by gamma0.

The angular data still obey C/q<.00352, q/C<400, eta_G=GX<8.5e-33,
M_n^2=n(n+q)/C, d_n=2(q+2n), M_n>18sqrt(n), and
|log(mu0 B)|<1e-8. These follow for the exact product from its enclosure,
not a finite-seed kernel transfer. The H2 normal cap is unchanged.

Let q_red(p) denote the canonical quantum shape response with the small
anomaly reduction, BEFORE restoring c_* below. The exact bubble identity
and bounds in [the product contour proof](stationary-product-contour-inverse.md),
sections 2-5, apply to these angular parameters:

    q_red(p)=p alpha_q+p^2 H(p), |alpha_q|<1.7e-8,
    H(p)=beta_q+p int dmu(tau)/(tau+p),
    eta_G/5<beta_q<2eta_G/3, mu>=0,
    supp mu subset (1404,infinity),
    int dmu(tau)/tau<6e-5 eta_G.                         (3)

Here is the hypothesis transfer, including the finite terms. The pole
residues and their summable inverse-cubic correction use only M>18 and
the displayed angular moments. The common finite coefficient uses
F_eff=4gamma_E-2H_q+4log(mu0 B). The new logarithm allowance still gives
-100<F_eff<-51, so the original strict beta_q inequalities follow
unchanged. The loose slope estimate is the compact-probe VALUE proof in
[the reduced response](stationary-reduced-cold-response.md), section 9:
same H2 cap, probe norm <4, remainder 8e22 M^-3, summed <(16/3)e22 C,
restoration <8e13 C and tail <5e6 C. All remain below 1e24 C; its
restoration is far within that reserve with |log(mu0 B)|<1e-8. Multiplying
by 2G/B^2 gives |alpha_q|<1.7e-8. No finite exterior inverse is used.

For clarity, the underlying per-mode poles are

    tau_j=(2M+2j)(2M+2j+3),
    A_j=(2/pi)(2M+2j+3/2)
             W(j+1)W(2M+j+1)/W(M+j+1)^2>0,
    W(x)=Gamma(x+1/2)/Gamma(x),
    Pi_M(p)=M^2 p sum_j A_j/[tau_j(tau_j+p)]-p/(12pi),
    Q_red=-q p/(12pi)-C[F_eff/(120pi)+1/(180pi)]p^2
                      +sum_n d_n[Pi_Mn+p^2/(120pi M_n^2)],
    q_red=(2G/B^2)Q_red.                                (4)

One common four-dimensional subtraction, the LLL and the stationary
contact are present. The primary bubble formula, [Carmi--Di Pietro--Komatsu
(6.25), (6.29), (D.10)](https://arxiv.org/html/1810.04185v3), leaves an
additive constant ambiguous. The difference at p=0 and the product Ward
contact fix that constant here; their Gross-Neveu finite coupling is not
imported. The rank-two normalization bridge remains the retained one.

## 2. Derive the full matrix, including the large local coefficient

Put delta_c=1/(5760pi^2), c_*=c_R-delta_c, and define

    beta=64piG c_* X,
    c=6+6beta(1-Y/X)-11GY/(90pi),
    d=2(R+24lambda0)/X,
    b=1-G(6X+5Y)/(180pi),
    a(ell)=9beta ell^2+c ell+d,
    j(ell)=b+(3/2)beta ell.                            (5)

These beta,b,c,d differ from names in earlier reduced-action notes.
They are fixed current-action coefficients, not fitted or tunable inputs.
The exact complete physical Hessian is gamma0 M, where

    M(p)=[[-a(ell), -2p j(ell)],
          [-2p j(ell), q_red(p)-beta p^2]], ell=p-2.   (6)

Derive the entries directly from the full action. In proper product
frequency k^2=X ell, the Euclidean int R^2 Hessian entries are

    H_SS^R2=72k^4-24R k^2,
    H_ST^R2=24k^2(k^2+2X),
    H_TT^R2=8k^4+32Xk^2+16X^2+16Y^2.                 (7)

Their coefficient in the physical Euclidean action is -c_R. The Dirac
conformal column equals the anomaly-complete column plus delta_c times
(7), giving precisely -c_* for the fourth-order terms. Differentiating
the complete anomaly gives

    H_SS^Gamma*=11Y k^2/(720pi^2),
    H_ST^Gamma*=[X^2+(5q^2-1)Y^2]/(120pi^2)
                                    +(6X+5Y)k^2/(720pi^2).

For the first line, delta(sqrt(g)Euler)=8sqrt(g)G_ab Hess^ab S,
G_ab=-Yg_ab in orbit directions, and the Euler anomaly coefficient is
-11/(5760pi^2). The conformal derivative of the Euclidean determinant
is minus the trace. C^2 and F^2 densities are conformally invariant.
This fixes the sign without a high-frequency extrapolation.

The EH/cosmological entries are

    H_SS^EH=-3k^2/(4piG)-(R+24lambda0)/(4piG),
    H_ST^EH=-(k^2+2Y+12lambda0)/(4piG).

Using the TRACE equation (1)-(2) factors their full mixed sum by
k^2+2X and gives -2p j. For the TT row use the ORBIT equation (1):
its local Einstein/Maxwell constant is 4p_A-16c_R(X^2-Y^2).
The canonical quantum contact is Gamma_D,TT+4p_A. The same contact
minus 8delta_c X^2 p^2 is gamma0 q_red. Combining with (7) therefore
gives gamma0(q_red-beta p^2). In particular the negative local p^2
term in (6) must not be dropped or called a change of subtraction.

This calculation uses exact product stationarity for its contacts.
The same formula is NOT an off-shell finite-seed Hessian. The current
auxiliary rows from [the preceding audit](physical-metric-block-and-preconditioner.md)
remain the appropriate implementation there, with all density contacts.

The exact product enclosure gives ample rational guards

    0<beta<1.18e-8, 5.99<c<6.01,
    3.99<d<4.01, .999<b<1.001.                         (8)

Indeed X=6lambda+Y+d_geometry Y, |d_geometry|/C<1.71e-29,
Y/X<2.1e-11 and G lambda<2000^-10. Equations (5) then imply (8),
including lambda0 rather than silently replacing it by lambda.

## 3. Exact determinant cancellation and a complete product inverse

Define

    B(ell)=4b^2-beta d+beta ell(12b-c),
    D(p)=B(ell)+a(ell)[H(p)+alpha_q/p].                 (9)

The identity 4j^2-beta a=B cancels the beta^2 ell^2 terms EXACTLY,
and

    det M=-p^2 D,
    (M^-1)_SS=[beta-H-alpha_q/p]/D,
    (M^-1)_ST=-2j/(pD), (M^-1)_TT=a/(p^2 D).           (10)

This is the required coupling cancellation. Bounding the large R^2
entries separately would miss it. It does not eliminate a physical
scalar branch: the full a,j and determinant are retained.

For real ell>=1/4, H>=0. The rational guards imply

    B>3.99+5.97beta ell,
    a/p<=9beta ell+7.8,
    D>3.98+5.96beta ell+aH>0.                          (11)

In the middle inequality use ell/p<=1 and d/p<4.01/(9/4).
The small negative allowance |alpha_q|a/p costs less than
1.326e-7+1.53e-7 beta ell, leaving the strict reserve in (11).
Also a>=5.4875. Equations (10)-(11) yield the uniform entry bounds

    |(M^-1)_SS| <1/5.4875+1.18e-8/3.98
                                     +1.7e-8/[(9/4)3.98] <.183,
    |(M^-1)_ST| <.504/p <=.224,
    |(M^-1)_TT| <.284.                               (12)

For the last inequality divide by 3.98+5.96beta ell and regard the
ratio as a convex combination of its beta ell=0 and infinite limits.
The two required maxima over ell>=1/4 are

    sup (6.01ell+4.01)/[3.98(ell+2)^2]
                         =6.01^2/[4*3.98*(12.02-4.01)]<.284,
    sup 9ell/[5.96(ell+2)^2]=9/(8*5.96)<.189.

The maximum absolute row sum of the symmetric inverse is <.508<.51.
The complete H2 spectral transform makes M(P) a real symmetric matrix
multiplier. Its maximal domain is {h in L2: M(P)h in L2}. This is a
closed self-adjoint realization by the spectral theorem; bounded
measurable inverse (10) maps every L2 source into that domain. Hence

    ||h||_(L2,dvol4;R2)<.51 ||gamma0^-1 f||_(L2,dvol4;R2). (13)

This constructs the complete PRODUCT inverse, including arbitrary
frequencies and the full positive Stieltjes remainder, rather than
inverting only its principal symbol. It neither identifies the physical
finite outer domain nor controls growing cold weights or raw C5 jets.

## 4. The cold roots of the complete current matrix

The determinant has the exact Ward factor p. Its other factor is

    F(p)=p[B(p-2)+a(p-2)H(p)]+a(p-2)alpha_q.            (14)

On |p|<=4e-8, direct polynomial bounds from (8) give |a|<8.04,
|B|>3.99. The Stieltjes representation gives

    |H(p)|<[(2/3)+(|p|/(1-|p|/1404))*6e-5]eta_G
           <.667 eta_G.

On the circle the perturbation a alpha_q+p aH has magnitude
<8.04*1.7e-8+4e-8*8.04*.667*8.5e-33,
strictly less than |pB|>4e-8*3.99. Since B has no zero in this disk,
Rouche's theorem proves exactly one F zero, counted with multiplicity.
Real coefficients make that zero real. If alpha_q=0 it coincides with
the Ward root; no artificial splitting is assumed. With
nu=(-1+sqrt(9-4p))/2 on the branch nu(0)=1,

    |p_root|<4e-8, |nu-1|<1.34e-8,
    exp(84|nu-1|)<1.00000113.                           (15)

This controls the possible exponent shift over the retained finite neck.
It does not bound a residue, source amplitude or boundary matching.
cosh chi is not an L2 zero mode and is not declared an allowed variation
of the finite exterior's mass or temperature. The positive scalar branch
and any other poles outside this small disk are not classified by (15).

## 5. Execute the long-neck source-envelope comparison

The L2 result does not answer the decisive cold-value question. Test the
retained unsigned envelope in the CLASSICAL coupled product, in the same
physical normalization gamma0. Its matrix is

    [[-(6ell+4),-2p],[-2p,0]].

For a pure canonical shape source F_T, the radius equation is exactly
-2P2 S=F_T. Set s=h cosh chi, let the observation seam be s_c=.001,
L=arcosh(s_c/h), and use the actual envelope

    a_s=s[1+log(.01/s)], F_T=a_s,
    1/(8*2000^10)<=h<=2000^-10.

Take a regular bolt and S(L)=0. This is a specified finite-neck
CLASSICAL comparison, not an imposed boundary condition on the full
nonlocal metric operator or on the physical exterior. A nonzero matched
seam value would add a homogeneous term, whose sign is not known.

Writing S=-(h cosh chi)v/2 gives

    (cosh^2 chi sinh chi v')'
             =-cosh^2 chi sinh chi [L0-log cosh chi],
    L0=1+log(.01/h), v(L)=0.

Thus the exact bolt gain in |S|/(h cosh chi), per unit envelope, is

    K_env(h)=(1/2) int_0^L J(cosh u)/[cosh^2 u sinh u] du,
    J(x)=[x^3(L0-log x+1/3)-(L0+1/3)]/3.                (16)

The apparent endpoint singularity is removable. This is the positive
Green solution, not an estimate obtained by multiplying independent
operator norms. The actual linear physical return S/tau adds 1+s;
at the bolt that is 1+h, so the comparison is already slightly smaller
than that return. It also avoids blaming an unknown seam amplitude for
the large gain: the seam value is zero in this explicit test.

There is an elementary uniform rigorous lower bound. For u>=1,
J(cosh u)>=((cosh^3 u-1)/3)[L0-log cosh u]. Since sinh u<=cosh u,
log cosh u<=u, 1-sech(1)^3>.72, L>69 and L0>72 for the whole h range,

    K_env >(.72/6) int_1^69 (72-u) du=301.92.          (17)

Direct 70/110-digit quadrature of (16) gives the endpoint illustrations

    K_env(2000^-10)=443.44127776367363352937426969...,
    K_env(2000^-10/8)=469.11675291031201196687018080... .

These displayed quadratures are numerical illustrations. Inequality
(17), based on interval-checked elementary guards and positive integrals,
is the certificate that the comparison gain exceeds 300. The gain is
monotone as h decreases, since both L and L0 increase in (16)'s original
positive double integral; these endpoint values characterize the model's
range but are not interval quadrature certificates.

At the retained source allowance 1.36e-9, the envelope-saturating
comparison gives radius response about 6.03e-7 to 6.38e-7. To allocate
2e-9 to this channel using the same unsigned envelope would instead
require its normalized coefficient about 4.26e-12 to 4.51e-12, over
300 times smaller. The rigorous lower gain alone rules out K<=1.4705882352
for this comparison value problem, and gives response >4.106112e-7.

This is NOT a lower bound on the actual quantum correction, nor an
operator lower bound for the full finite geometry. The actual signed
source need not saturate the envelope, the full matrix differs, and
matched seam data may cancel part of the response. It demonstrates why
the old pointwise source bound plus a hoped-for unit cold inverse is not
a credible feasibility certificate. Small cold exponent displacement in
(15) does not remove this nearly resonant forced accumulation.

## 6. Close the audit and select one decisive calculation

The earlier principal-symbol decision is now supported by an exact
current-action product inverse: retain the coupled matrix, the quantum
logarithm and the determinant cancellation. No product L2 zero or large
cold exponent is forced by the restored local coefficient. A positive
Euclidean Schur denominator does not establish retarded stability.

The source-driven finite correction must still be computed with its
mass/clock boundary phase and physical radius return. In particular
S=delta log r= tau delta V at the finite seed; tau tends to the tiny
thermal gap at the bolt. The L2 norm in (13) does not control S/tau,
the growing cold source a_s, seam values, outer tails or their C5 jets.
Multiplying .51 by a pointwise canonical source allowance is invalid.

Keep the unchanged target

    1.632e-11 K_m+1.36e-9 K_c+Y_rest<=5e-9,
    K_m<=61.2745098039, K_c<=1.4705882352, Y_rest<=2e-9,
    rho=1e-8, Y0<=5e-9, Z0<=.25, Z2<=5e7,
    Y0+Z0 rho+(Z2/2)rho^2<=rho, Z0+Z2 rho<1.

These physical constants remain unknown. The strongest outer contract
remains b_o+1e-8 E_o<=1.5e-7 with its transfer hypotheses. Exact range
compatibility, full source vertices and actual boundary realization
remain prior to claiming a nonlinear solution.

The principal product conditions are realistically attainable: they are
now proved in their stated spectral space. The required cold physical
return target is not supported by the unsigned-envelope route; (16)
exhibits the loss already in its classical limit. A useful cancellation
must concern the actual signed residual and its boundary return, rather
than the positivity or size of a different norm. The product Ward factor
suggests a specific possibility: its leading growing orbit variation is
a diffeomorphism, so stationary contact terms cancel that component in
the COMPLETE product residual. It has not been proved that the finite
seed and its actual matching preserve the needed cancellation. Counting
the canonical stress difference alone can miss the matching contacts.

Prior obstacles are the definition and boundedness of that physical
boundary-to-bolt map: fixed clock/temperature, correlated mass and scalar
fast phase, bolt regularity, common quantum boundary realization and
weighted C5 tails. A mere N2 bound has no generic C5 upgrade. No usable
full finite inverse, projection compatibility theorem, or physical
nonlinear radius has been supplied. Their targets remain conditional;
the present result does not establish that their constants are attainable.

DECISION: retain the fixed construction, but park the generic unsigned
source-times-inverse estimate chain. Select ONE next calculation: enclose
the signed bolt-radius response functional of the COMPLETE matched
finite-seed residual, including the seam/mass/clock and scalar fast-phase
contribution and a weighted remainder. Use the full current product
kernel as the comparison, and test the joint response against the 2e-9
cold allocation. This tests whether Ward/contact and boundary cancellation
actually remove the accumulation in (16). If the required map cannot be
realized in the prescribed domain, that is the prior obstacle to resolve;
further source tightening should wait.

The product inverse, cold-root test and classical envelope comparison
have been executed here. The signed finite-seed response calculation has
NOT. Overall feasibility remains undecided; the audit closes with this
specified sharper method and no assertion of existence or nonexistence.
No further estimate chain is selected in this audit. The full stable
wormhole objective remains unchanged.

Remaining to-do list: actual weighted correction and residual/tails;
finite-geometry metric domain/inverse, boundary transfer, compatibility
and nonlinear admission; matching, semiclassical validity, all-sector
stability, charge/return, opening, transport, preparation, separate
creation and final audit.
