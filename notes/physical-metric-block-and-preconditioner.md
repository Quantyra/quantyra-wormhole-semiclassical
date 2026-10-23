# Physical metric equations and the required preconditioner

2026-10-22, E003/S169. Ordinary argument before supporting verification.
Progress on the full stable semiclassical wormhole objective, not a
controlled solution. The previous canonical matter neighborhood does not
bound the metric Hessian or its inverse in the physical admission norm.

## 1. Exact bulk equations in the current physical variables

Keep the same c_R=10^22 action, mu0=1/b0, flux, prepared positive
temperature, product B,A and paired boundary realization. Use the exact
[anomaly repartition](anomaly-complete-auxiliary-response.md):

    delta_c=1/(5760pi^2), c_*=c_R-delta_c,
    F_*=1-384piG c_*lambda, D=96piG c_*,
    Phi=F_*exp(2u), phi=(F_*exp(2u)-1)/(32piG c_*),
    Gamma_*ren=Gamma_D+delta_c int R^2-p_v Vol_g,
    p_v=-3(lambda0-lambda)/(8piG).

The last line is Lorentzian notation; the Euclidean R^2 addition is
-delta_c int R^2. Include the same repartition in the boundary action.
This is an exact organization of the existing action, not a new finite
coefficient or scale choice. The Einstein-frame Euclidean bulk action is

    I_E=F_*/(16piG) int sqrt(g_E)
          [-R_E-6lambda+6|du|_E^2+W(u)]
          +I_EM,E+Gamma_*ren[exp(-2u)g_E,A],
    W(u)=(1-exp(-2u))^2/(64piG c_*).

Use gamma0 for the actual product-centered seed orbit in units A:

    g_E=exp(2P)[A^2 exp(2b)gamma0+B^2 dOmega^2],
    p=P-u=log(r/B), r=B+(r0-B)exp(V), A^2/B^2=1/C.

The unperturbed values are b=0, V=0 and P=log(r0/B)+u_b. All derivatives
below use gamma0, and R0=R[gamma0]. Set

    K=F_*B^2/(4G), d=lambda A^2, w=A^2 W,
    Q=2piG A^2 q^2/(F_*g_R^2 B^4),
    U(P,u)=-(2/C)exp(2P)+(w(u)-6d)exp(4P)+Q,
    Psi=Gamma_*ren[exp(2(P-u))(A^2 exp(2b)gamma0+B^2 dOmega^2)]/K.

The magnetic convention is F=(q/2)sin(theta)dtheta wedge dphi, hence
F_ab F^ab=q^2/(2r^4). Its integrated Euclidean coefficient is precisely
Q above. Sphere integration and one integration by parts give I_E/K
as the integral of the following density against sqrt(gamma0), plus Psi:

    L=-exp(2P)R0-4exp(2P)dP.db-6exp(2P)|dP|^2
                         +6exp(2P)|du|^2+exp(2b)U(P,u).       (1)

For example R[g_E]=exp(-2P){A^-2 exp(-2b)
[R0-2Delta b-6Delta P-6|dP|^2]+2/B^2}. The Weyl Laplacian flux cancels
the corresponding transformed GHY flux. On a finite orbit boundary the
remaining local GHY term is -2K int exp(2P)(k0+n0 b) ds0. Keep the
quantum boundary cocycle, counterterms and physical asymptotic boundary
phase as well. Compact variations suffice for the following bulk rows;
they do not determine or change endpoint conditions.

Writing Psi_P etc. for functional derivatives with reference density,
the exact three rows are

    E_b=4 div(exp(2P)dP)+2exp(2b)U+Psi_b,
    E_P=exp(2P)[-2R0+4Delta b+12Delta P
                          +12|dP|^2+12|du|^2]+exp(2b)U_P+Psi_P,
    E_u=-12 div(exp(2P)du)+exp(2b)U_u+Psi_u.             (2)

These are off-shell equations. No seed residual, varying density or
coordinate-pullback contact has been set to zero. Psi includes the
complete canonical nonlocal determinant, local Weyl/radius functional,
vacuum subtraction and the local anomaly repartition.

Since Psi depends on P,u only through P-u, Psi_P+Psi_u=0 exactly.
More precisely, the sum of the last two rows is

    E_P+E_u=-2 A^2 exp(4P-2u+2b)(R[g]-phi).             (3)

To check it directly, use

    4W+W'-24lambda=2exp(-2u)phi,
    R[g]=exp(-2p){A^-2 exp(-2b)
                  [R0-2Delta b-6Delta p-6|dp|^2]+2/B^2}.

Thus the curvature constraint is retained in the actual metric equations,
with a nonzero coefficient. For estimates one may use the previously
proved massive scalar row instead, PROVIDED the exact traced-row
equivalence and its contacts are used: the constraint-recovery factor
is >1-2e-8 on |phi+12lambda|<=lambda, |R-phi|<=lambda. Those inequalities
remain separate admission requirements. Equation (3) is an explicit
implementation of the existing auxiliary equivalence, not an extra
independent equation to impose after solving an inconsistent system.

## 2. The retained cancellation and the mixed differential orders

Freeze a smooth compact interior background, use proper Einstein-frame
radial frequency k, and let

    theta=delta P+delta b/3, v=delta u.

The derivative square in (1) is exactly

    -6|dP|^2-4dP.db=-6|d(P+b/3)|^2+(2/3)|db|^2.

Consequently the local second-variation symbol, per physical Einstein
volume, has coefficients -a k^2, b_E k^2, a k^2 on theta,b,v, where

    a=3F_*/(4piG), b_E=F_*/(12piG).

The full Dirac interior symbol in these disk variables, as normalized in
[the retained shape calculation](stationary-shape-logarithmic-response.md), is

    k^4(theta-v)^2/(80pi^2)
      +k^4 b^2[log(k/mu_E)-6/5]/(120pi^2),
    mu_E=mu0 exp(-u_b).                                  (4)

Here the change of proper frequency and of volume cancels in the fourth
order coefficient. mu_E follows from the fixed physical mu0; it is not
a scale retuning. The Euclidean -delta_c int R^2 Hessian contributes
-72delta_c k^4(theta-v)^2. Because 72delta_c=1/(80pi^2), the conformal
fourth-order term in (4) cancels EXACTLY. The shape logarithm remains.
No part of the full quantum shape response is discarded by this identity.
Curvature-dependent second-order terms, state-dependent lower terms,
residual contacts and the nonlocal remainder still belong to the operator.

Moller, [Theorem 1 and Lemma 6.1](https://arxiv.org/html/0705.3857v3),
supplies a local Dirac zeta-Hessian symbol on closed smooth spin manifolds
with stable kernel dimension. Its global spectral theorem is not a
theorem for this boundary problem. The retained compact heat-parametrix
argument transfers the local leading symbol only; it explicitly retains
lower state/domain terms. The actual paired first-order domain stability
established in the preceding audit does not supply a renormalized metric
Hessian graph bound or justify summing divergent angular operator bounds.

A candidate frozen-frequency preconditioner therefore has model symbol

    diag(-a k^2,
         b_E k^2+k^4[log(k/mu_E)-6/5]/(120pi^2),
         a[k^2+m0^2]),
    m0^2=1/D.                                           (5)

The scalar mass in (5) is the local vacuum mass; the complete scalar row
has C_v/D and the already bounded seed corrections. A global realization
of (5) needs low-frequency completion, boundary conditions and weights.
In particular its extension of the ultraviolet expression to small k is
a MODEL choice, not an evaluation of the complete low-frequency Hessian.

## 3. Execute the frequency test before choosing the solver

Use t=Gk^2, z=Gmu_E^2. The shape quantum/Einstein principal ratio is

    r(t)=t[0.5log(t/z)-6/5]/(10pi F_*).                    (6)

With benchmark kappa_benchmark=2000 (distinct from frequency k),

    g0=2000^-10, g+11g^2/(360pi)=g0,
    C0=sqrt(1+12pi*2000^6), x0=(C0-1)/6,
    z0=Gmu0^2=g0/x0, z=z0 exp(-2u_b), |u_b|<5.4e-37,
    F_*=1-384pi(10^22-1/(5760pi^2))g.

Numerically at u_b=0,

    F_*=0.9999999882190275490382753557650873334...,
    z0=1.1928772263497732068596907614445158211e-43,
    m0^2G=3.3157279810811528285184117427599083007e-25,
    m0^2/lambda=339530545.2627100496402853624586... .

The unique positive r=1 crossover has the explicit expression

    t_*=20pi F_*/W0(20pi F_*/[z exp(12/5)]),

where W0 is the positive Lambert branch. Direct monotonic bracketing
(not reliance on the Lambert implementation) gives, including the whole
u_b interval,

    0.65442489137407894964960401158
          < t_* <0.65442489137407894964960401159,
    k_*sqrt(G)=0.8089653214904078628662897526...,
    r(1)=1.5348073818309983800810189182... .

Since r(t) diverges, preconditioning the quantum shape term with an
Einstein second-order inverse cannot give a bounded, let alone <1/4,
operator defect in the corresponding same-order Sobolev graph. A compact
oscillatory shape family supported away from the ends proves this: the
fourth-log term is proportional to n^4 log n whereas the classical
principal inverse contributes n^-2; lower symbol orders cannot cancel
the nonzero leading coefficient. This is an all-frequency statement in
the continuum model, not a claim about a physical excitation at the
Planck scale. The high-frequency ratio is negligible at the scalar mass,
but that does not make its supremum small.

The MODEL in (5) has no negative shape zero: the minimum of r over t>0
is at t=z exp(7/5) and equals -z exp(7/5)/(20pi F_*), with magnitude
<7.70e-45. Thus its shape multiplier has pointwise inverse allowance
at most 1/(1-7.70e-45) relative to b_E k^2. This elementary multiplier
fact gives no global C5 estimate, no boundary realization, and no
coercivity conclusion for the actual low-frequency quantum operator.

The crossover is near the Planck frequency, where semiclassical validity
has not been established. It is not evidence that the resulting field
contains appreciable power there. Conversely a cutoff cannot simply be
introduced into the specified equations: an EFT/order-reduction change
requires a documented theory and error analysis. Retarded stability is
not determined by Euclidean multiplier signs.

## 4. Physical return, feasibility ledger and research decision

The above diagonalization changes variables, not the admission norm.
With tau=(r0-B)/r0=s/t_radius, the exact physical radius return is

    delta V=(theta-v-b/3)/tau,
    V=log[1+(exp(zeta)-1)/tau],
    zeta=(P-P_b)-(u-u_b).

The first formula is the tangent formula; the second is nonlinear.
Preserving |V|<=rho is ensured by |zeta|/tau<=1-exp(-rho), and
monotonic radius additionally needs 1+s V_t>0. Taking separate norms
of theta,v,b before this combination can destroy the cold cancellation.
The complete physical metric inverse and this joint return are UNKNOWN.

The existing canonical middle/cold C3 allowances in original units are
6e15 q^2 and 5e17 q^2 a_s, with conversion <2.72e-27, hence normalized
1.632e-11 and 1.36e-9 a_s. These are not the complete physical sources.
The complete regional returns must meet

    1.632e-11 K_m+1.36e-9 K_c+Y_rest<=5e-9,
    K_m<=61.2745098039, K_c<=1.4705882352, Y_rest<=2e-9

as one sufficient allocation. They are unknown constants, not supplied
by the pointwise MODEL inverse above or by the matter inverse 1/M.
At rho=1e-8 the complete projected Newton map requires

    Y0<=5e-9, Z0<=.25, Z2<=5e7,
    Y0+Z0*rho+(Z2/2)*rho^2<=rho,
    Z0+Z2*rho<1.                                         (7)

Use a fixed physical projection with proved range/compatibility, or
identity until such a projection is derived for this action. No old
zero-c signed projection has been transferred. All domain constraints
and actual-boundary admission margins must hold along the whole ball:

    b45+1e-8 E45<=1e-3 (E45<=87000 if b45<=1.3e-4),
    b_c+1e-8 E_c<=5e-7 (b_c<=4.9e-7 if E_c<=1),
    b_o+1e-8 E_o<=1.5e-7 (b_o<=1.4e-7 if E_o<=1).

The last contract is the stronger retained outer theorem, still awaiting
S169 proper-collar/global-C1/actual-domain transfer. Its original
background allowance 1.01e-7 leaves 4.9e-8, not a negative margin.
The first-jet matter domain is controlled at |b|<=1e-8, |db|<=1e-6;
that is one ingredient, not the full admissible metric neighborhood.

For comparison the old complete interior H4 allowance remains
C_I=9.1591500000000000015e103, with global C_O=2.365283519093731664e186.
The old converted source allowances are C_F,I=4.141944868545865e118
and C_F=1.069626977892751e201. Even discarding the old low-order floor,
that certificate needs C_F<=2.38393e-25 at R=.97, C_D4=2: the interior
and global allowances miss by >1.7374e143 and >4.4868e225. Keeping the
floor instead requires R<3.90858e-8 while R>.97. These are failed
sufficient bounds in the original zero-c route, not lower bounds on
the new correction and not nonexistence theorems. The newer regional C3
estimates must not be substituted for a complete H4 norm.

DECISION: continue the fixed construction using the explicit equations
(2) and a mixed-order preconditioner retaining the full shape logarithm.
Reject an all-frequency second-order small-quantum perturbation test.
This decision follows from the executed frequency calculation, not from
new normalization or passing checks. Keep the scalar massive row and
measure the joint physical radius correction directly. The local
cancellation in section 2 makes this approach structurally plausible;
it does not yet show the numerical regional return targets are realistic.

The decisive remaining feasibility calculation is the complete physical
correction and residual on the actual finite geometry, including low
frequencies, contact vertices, common angular subtraction, mass/clock
boundary phase and outer tail. This note derives its explicit local
rows but does NOT claim to have executed that global solve. The next
calculation must give a certified weighted remainder and either meet
K_m,K_c,Y_rest and (7), or document the specific failed bound/mode.
An approximate correction with residual <=1.25e-9 would give error
<=5e-9 only after a preserved ball and contraction <=.75 are proved.

Compatibility, seed-to-finite matching and the physical boundary domain
are prior to adopting an inverse; neither a closed-manifold theorem nor
the canonical matter domain alone supplies them. Small G|Riem| on the
seed and a positive local scalar mass do not control induced fluctuations
or all-sector retarded dynamics. No additional theory switch or claim
of physical validity, stability, existence or nonexistence is selected.

Remaining to-do list: actual full physical correction and residual/tails;
mixed-order metric domain/inverse and outer transfer; compatibility and
nonlinear admission; matching, semiclassical validity, all-sector
stability, charge/return, opening, transport, preparation, separate
creation and the final audit.
