# Test the uniform bolt certificate against the actual trace source

2026-10-23, E003/S169. Ordinary argument before frozen execution.
Previous turn: progress. The full stable-wormhole objective is active and
incomplete; the feasibility audit stays closed.

The global pair has not been computed. Before doing that expensive solve,
test whether ANY pair could pass the currently selected full-V uniform
remainder certificate. A signed projection of the actual trace source
can bound the necessary norms from below. This is a test of a sufficient
certificate, not a lower bound on the actual physical bolt error.

## 1. A global displacement estimate with the actual thermal root

Use t=r/B, z=1/t, v=1+h, z_h=1/v, x=lambda B^2, u=1/x and d as in
product-centered-scalar-cold-response.md. The ACTUAL seed has

    F(t)/x=(t-1)^2[1+2/t+(3+u)/t^2+du/t^4]-delta/(xt),
    delta=v F_ext(v)>0, t>=v.

The polynomial in 2(t-1)F_t-3F, apart from its positive (t-1)^2
factor and the positive thermal term, is

    p(z)=1-2z-5z^2+12z^3+u(-3z^2+4z^3)+du(-7z^4+8z^5).

Exact Bernstein subdivision on [0,1], with u<1/(8e9), |d|<8.6e-19,
proves p>0. Thus -(2(t-1)/F)' >=1/F and the expression tends to
zero at infinity. Integration gives the global physical inequality

    F(r) integral_r^infinity dr'/F(r') <=2(r-B).       (1)

The denominator is r-B, NOT r-R_bolt: no false extremal comparison is
made at the nonextremal root. For the existing actual reference correction,

    xi=s_R-F integral_r^infinity H/F,
    |H|<1.08e-9, |a|<3.3e-10, |s_R|/B<3.33e-10 h,

from current-seed-directional-cold-response.md, hence

    |xi|/B <=3.33e-10 h+2(1.08e-9)(t-1).             (2)

This is a sharper estimate of the existing correction; it introduces no
new source-response derivative or normalization. The preserved paired
quantum-source and Ward arguments are inputs to these constants.

## 2. Compute a signed projection with all density and endpoint terms

Test the full defect r1 against the actual form vector v_s=(z^4,0).
It is smooth at the bolt, decays at the outer end and lies in V. It
adds no artificial wall or extra boundary condition. In the fixed affine
chart with the frozen York map, its scalar component is just nu=z^4.

Let T0/lambda=tau(z), epsilon=1/(60pi), D=96pi c_R-epsilon,
eta=epsilon lambda, d_l=D lambda, S=epsilon q^2/(lambda B^4).
Use the product identity and actual thermal polynomial from the full
trace argument, not a zero-temperature replacement:

    tau=S A+eta Bp-d_l du Lr,
    A=-10z^4-20z^5+30z^6=-10z^4(1-z)(1+3z),
    r=4z^5-6z^6, Lf=z^4 partial_z[(F/x)partial_z f],
    rn=-12+du r,
    tn=(3+u+du)z^4-du(5z^5-9z^6/2),
    cn=-z^3 partial_z(F/x)+rn/6,
    P=7cn^2+(22/3)tn^2-11(rn^2-144)/72,
    Ps=2(6+du)^2+(22/3)(6+du)u+(22/3)u^2+22,
    Bp=P+Ps r/2.

The exact product equation du=-5S+eta Ps/2 justifies this correlated
representation. The finite-temperature mass term stays in F, cn and L.
The complete known areal term has normalized integral

    I_known=integral_0^z_h [d_l Ltau+eta Kt tau]dz,
    Kt=11rn/36-7cn/3+11tn/3.                         (3)

The remaining areal terms are Af f/r^2+AQ Q/r^2, with Q=G(rho-p),
not the seed orbit polynomial. The full trace proof gives, since w<=1,

    err_areal <=(490eta+316d_l |d|/x)(2.11e-6)
                 +(535eta+56d_l |d|/x)(2.01e-7).    (4)

These are known upper bounds on unevaluated terms; their signs are NOT
assumed. Their source inputs retain the complete quantum tower.

The actual disk density is
e_disk=e_T+a T0+r^-2 partial_r(r^2 xi T0). After integration against
z^4 with the physical radial measure, the z^-4 measure cancels. Integrate
the pullback term once, keeping its nonzero bolt endpoint:

    integral_0^z_h e_disk/lambda dz
       = integral_0^z_h [e_T/lambda+a tau+4z(xi/B)tau]dz
                           -z_h^2(s_R/B)tau(z_h).    (5)

The outer term vanishes by the actual decay. No derivative of an unknown
xi is estimated. Write I0>=int|tau|, I1>=int(1-z)|tau| over [0,z_h].
Because A<=0 on [0,1], its exact contributions are S*22/21 and S*23/84;
bound the small remainder coefficient by coefficient on [0,1]. Then

    err_density <=(3.3e-10+4*3.33e-10 h) I0
           +8(1.08e-9) I1+3.33e-10 h |tau(z_h)|.    (6)

In the existing action units W=-N pi B^2/(G gamma), the corrected PLUS
sign for the normalized trace row gives

    <r1,v_s>=W I,
    I=(x/C) integral_0^z_h e_disk/lambda dz,
    I >=(x/C)[I_known-err_areal-err_density].         (7)

Freeze a positive lower target I>1e-38. The unknown shape source has no
pairing with this scalar test. This is an actual signed source projection,
not a norm inferred by differentiating a source value bound.

## 3. An obstruction to the current certificate for EVERY trial pair

The complete scalar quadratic form is a(v_s,v_s)=|W| A44, where A44
is the exact scalar_z4 entry of current-finite-local-action-assembly.md.
It includes Einstein, R2 and anomaly terms; the canonical shape
determinant has no scalar row. Also L_b(v_s)=ell_b z_h^4,
ell_b=(1+h)/h. Keep these physical factors unchanged.

In the positive diagonal V norm, the actual Hessian has normalized block
matrix [-I,Cmix;Cmix*,I], ||Cmix||^2<.338. Squaring it gives
||H||<=sqrt(1.338). For the numerical comparison Hc, the already proved
||H-Hc||<.0022 yields ||Hc||<M=1.16.

Let R0=|<r1,v_s>|/sqrt(a44), L0=|L_b(v_s)|/sqrt(a44). For ANY trial
pair with norms E,P and comparison residual norms D0,Db0,

    M E+D0>=R0, M P+Db0>=L0.

Therefore its sufficient remainder expression satisfies

    eta E P+(Db0+eta P)(D0+eta E)
        >=eta(1+eta) R0 L0/M^2,
    R0 L0/M^2 >= I_lower ell_b z_h^4/(M^2 A44_upper). (8)

Here eta=.0022 and M>1+eta. To prove the first line minimize in the
nonnegative residuals, then put E=R0 u/M and P=L0 v/M for u,v in
[0,1]. The expression is bilinear on that square. Its four corner values
are 1, eta/M, eta/M and eta(1+eta)/M^2, in units R0 L0. The last is
the minimum. Values outside the square can only increase the expression.
Thus this is not merely a restriction on exact model solves.

Freeze the test R0 L0/M^2>1e-5 and compare (8) with the ENTIRE target
1.667e-9, even before |Q0| is added. If it fails, no trial refinement
can make this particular full-V/.0022 certificate succeed. Calculate
the necessary uniform eta thresholds for the entire and half budgets.
Those thresholds are upper necessities, not proofs that such an operator
bound exists or would suffice with nonzero Q0 and other errors.

This does NOT prove a nonzero actual correction, nonexistence, instability,
or that ||H-Hc|| is as large as its certified upper bound. The actual
signed error can cancel and the operator difference acts only in TT.

## 4. Use that exact block structure before optimizing sources

Let Hc=[-A,B;B*,Dc] and Rc=H-Hc, with Rc only TT/TT in the frozen
chart. Eliminate the unchanged scalar block EXACTLY in its actual form
realization. Define

    Sc=Dc+B* A^-1 B,
    r_eff=r_T+B* A^-1 r_s,
    ell_eff=ell_T+B* A^-1 ell_s,
    q_s=<ell_s,A^-1 r_s>.

Then Sc e_T=-r_eff, Sc psi_T=ell_eff. The physical bolt correction is
q_s+<ell_eff,e_T>, including the effective point source and all boundary
returns. The actual Schur form is S=Sc+Rc and S>=l. The comparison is
positive with Sc>=(1-.0022)l. No new domain or zero-mass condition is
introduced; the point source remains in the scalar form dual.

The existing scalar-only inverse estimate gives an actual bolt allowance

    |q_s| <=(1+h)[1+log(.01/h)]*(4/3)e18*(1.1e-31). (9)

The cold norm is divided by a_s, so its previously quoted 1.4667e-13
is NOT itself the normalized bolt value. Equation (9) retains the extra
L_h factor. Require this bolt allowance <1.2e-11.

For trial TT fields, with all norms now in l and its dual, the correct
block-specific certificate is

    |q_s+Q_T|+eta E_T P_T
               +(Db_T+eta P_T)(D_T+eta E_T)<1.667e-9, (10)
    d_T=Sc e_T+r_eff, db_T=ell_eff-Sc psi_T,
    Q_T=ell_eff(e_T)-<psi_T,d_T>.

The proof is the original two-residual identity applied to S. This uses
the known zero scalar/mixed blocks of Rc; it is not a renamed norm.
The lower bound (8) on FULL V norms does not force a corresponding TT
lower bound. Do not assume the TT norms or cancellations are small.

Decision if the frozen obstruction passes: park the full-V uniform-error
certificate and pursue the exact scalar Schur reduction with directional
TT residuals. The single next calculation is the actual finite-geometry
Schur forward/bolt-adjoint pair, including r_eff, ell_eff and its two
certified TT residuals. This full physical solve remains UNEXECUTED.
The source projection/certificate feasibility test is executed only by
its retained receipt. If the TT certificate fails too, compute direct
signed Rc pairings before extending a generic estimate chain.

Remaining to-do list: actual Schur forward/adjoint pair and physical
bolt/exterior certificate; X/admission, nonlinear closure and exact
compatibility; matching and causal realization, semiclassical validity/noise,
all-sector stability, charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
