# Complete outer-input coupling, with the global York return

2026-10-27 HST, E003/S169. Ordinary argument before interval execution.
The preceding complement attempt supplied no new bound. This calculation
tests a specific proposed domain split for the complete S_Q application.
It is not a cold quantum-column optimization or an inverse implementation.
Its contribution to the complete residual is given in section 4.

## 1. Actual domain and the first missing implication

Keep E_t with ||u||_Et^2=ell_t[u,u]/N, N=|W_action| beta, and
H=L2(z^-4 dz), z=B_radius/r, beta=lambda B_radius^2/C. The accepted
constants are E_t^2>=c_H H^2, ||C_H||_(A->H)<=K_b and
B_H=C_H a^-1 b, where A is the scalar form space with norm a/N.
Both original equations, all density contacts and the outer-fixed York
return are retained. The seed paired mass/log realization, common
regulator, exact UV logarithm, all angular/thermal modes, fixed period,
cutoff transitions and regular bolt remain unchanged.

Fix alpha=1e-9 BEFORE execution. Let chi be any smooth multiplier with
0<=chi<=1, chi=1 for z<=alpha/2 and chi=0 for z>=alpha. Its transition
is included. On the actual core chi u has the original mass/log behavior
and is zero near the bolt. The goal is to bound B_H(chi u), including its
nonlocal returns, uniformly for the ENTIRE unit E_t input ball.

This does not impose a finite wall, remove a boundary phase, or restrict
the output domain. Multiplication is used only in the input of the
already established bounded H realization of B_H. Its extension from the
core is controlled in H, so no unproved norm for chi:E_t->E_t is needed.
The eventual output still solves on the full paired domain.

The first missing implication is localization of the mixed form. Its
York return is global even for an outer-supported input. A pointwise
local-coefficient estimate alone would miss that term. We bound it next.

## 2. Retain and estimate the actual global return

Write the genuine TT lift as diag(-b_rad-2xi,b_rad,xi,xi), where

    (r^2 F b_rad)'=-(r^2 F'-2rF)xi.

Bolt regularity fixes the integral from the bolt. If xi is supported
in z<=alpha, then b_rad is also supported there; the accepted bound
||b_rad||_H<=2||xi||_H still applies. The outer-fixed disk return obeys
eta'_rho-A_orbit eta=b_rad+xi, d rho=dr/sqrt(F). With f=b_rad+xi,

    sqrt(lambda) eta(z)
       =-sqrt(beta) sqrt(Q(z))/z int_0^z f(s)/Q(s) ds.       (1)

This is the ORIGINAL return, not a homogeneous solution set to zero.
In the interior its integral is generally nonzero. Its sqrt(Q) factor
has the correct regular-bolt behavior. On 0<=z<=alpha let q_min be a
proved positive lower bound for Q. Cauchy--Schwarz gives

    |int_0^z f/Q|^2 <= ||f||_H^2 min(z,alpha)^5/(5q_min^2).

Let I_Q be any upper bound on int_0^z_h z^2 Q(z)dz. Then, since z_h<1,

    ||z^4 sqrt(lambda) eta||_H
       <= K_eta alpha^(5/2) ||xi||_H,
    K_eta=3 sqrt[beta I_Q/(5q_min^2)].                      (2)

Use I_Q<=sum_k |Q_k|/(k+3), not an uncertain endpoint cancellation.
All transitions and the full interior return occur in (1)-(2).

The full local trace T0 has an exact z^4 factor. With x=lambda B^2,
d_seed=C Q6, e=lambda/(60pi), d=D_R lambda, and G=1, set

    wbar=(d_seed/x)(4z^5-6z^6),
    taubar=[(1+3x+d_seed)z^4-d_seed(5z^5-9z^6/2)]/x,
    wCbar=-(C Q3/(2x))z^3-(C Q4/x)z^4
                         +(10d_seed/(3x))z^5-(5d_seed/(2x))z^6,
    Deltabar p=(z^2/beta)[Qp''+(Q'-2Q/z)p'],
    v=q^2/(60pi lambda B^4), C_v=1+11e/3.

The original exact vacuum relation
12(lambda0/lambda-1)=22e removes the CONSTANT term, giving

    T0/lambda=C_v wbar-(11e/72)wbar^2+28e wCbar^2
                 +(22e/3)taubar^2-10v z^4-d Deltabar wbar. (3)

This is the complete anomaly/R2 trace, not the traceless quantum response.
Its use in the mixed return is exactly the accepted fixed-chart Ward
identity. The canonical massive determinant has no scalar/mixed derivative
in that chart; it remains in the TT operator and its complete comparison.
No additional renormalization or vacuum subtraction is selected here.
Execute the absolute polynomial coefficient sum K_T, proving
|T0|/lambda<=K_T z^4 throughout the actual geometry. Equations (2)-(3)
then give the full return contribution K_eta K_T h1 alpha^(5/2), with
h1=1/sqrt(5.99). No derivative of T0 or pointwise singular bolt connection
is required.

## 3. Both rows give the complete localized source bound

Retain the constants of current-structured-source-graph-test.md. Its
covariant transpose and projector estimates retain these spatial powers:
tau=O(z^4); C4,grad C4,Hess C4=O(z^3); J,grad J,Hess J=O(z^4);
w,Hess w=O(z^5); the magnetic term is O(z^4). These are quantitative
polynomial bounds on this same seed, not regularity imported from RN.
The coefficient sums for j/lambda and wC/lambda are checked against the
previous (3.6,20,108) and (5.1,40,153) guards. Their covariant tensor
derivatives follow the retained projector identities, including the bolt.

For the main Laplacian term, integrate onto the scalar test BEFORE
estimating; its factor is tau xi. For the remaining local transposes,
H_TT xi is supported in z<=alpha. Thus all cutoff derivatives remain in
the exact weak identity; none has been separately discarded. Combining
these terms and (2) gives

    |b(s,xi)|/N <= K_alpha ||s||_A ||xi||_H,
    K_alpha=K_main alpha^4+K_curv alpha^3+K_w alpha^5
             +K_Hessw alpha^5+K_EM alpha^4
             +K_eta K_T h1 alpha^(5/2).                    (4)

The w term could retain a higher power; the stated power is a safe bound.
K_main,...,K_EM are the ESTABLISHED complete transpose allocations. The
return allocation is replaced by its global weighted bound, not omitted.

Solve the COMPLETE scalar equation a(s_xi,v)=b(v,xi). Its variational
definition and (4) give ||s_xi||_A<=K_alpha||xi||_H. Apply the global
mixed bound in the other slot and the original H-to-E_t comparison:

    R_outer:=||B_N M_chi||_(Et->Et*)<=K_b K_alpha/c_H,
    B_N=B/N.                                               (5)

All scalar returns remain in a^-1. Equation (5) is a uniform operator
bound on actual structured sources, not one column or an arbitrary-source
surrogate. It extends to E_t by the bounded H source realization.

## 4. Explicit contribution to the requested complete bound

Let S_Q=(ell_Q+B)^-1 B and S_inner=S_Q M_(1-chi), both on the FULL output
domain. No numerical S_inner is claimed. The actual coercivity proves

    ||S_Q-S_inner||_(Et->Et)<=a_outer=R_outer/(1-e_Q).       (6)

The known 12-term accounting converts this attained omission error to
the complete residual debit

    debit_outer=(1+q_t) A12 a_outer.                        (7)

Here the stored A12 bounds the exact retained K powers and also allows
the previously specified K-application screen. No actual K implementation
is assigned zero error. Check additionally
||S_inner||<=q_t/(1-e_Q)+a_outer against the existing, slightly larger
q_t/(1-e0) caller bound, e0=(.0022+epsilon)/(1-epsilon). This ensures the
input split does not silently enlarge the separate K-error multiplier.

Execute (4)-(7) and subtract a_outer from the sufficient 5e-32 S-application
allowance. The resulting allowance for approximating S_inner remains
CONDITIONAL; its numerical inverse, quantum applications and arithmetic
errors are still unknown. K needs its separate application certificate.
If (6) fits, the input split is justified at the stated cost. It does
not justify truncating output tails or replacing the actual boundary
problem by a wall at alpha. B_H applied to an interior input is generally
global, and so is its inverse output.

The first complement-inverse problem therefore persists. This result can
close only the outer-input omission error in that proposed implementation,
not its main solve or the complete T_t-(I+T_t)F_n requirement. In particular
no mass coefficient, independent-phase output, original geometric norm,
nonlinear admission/contraction, exact compatibility, physical matching,
semiclassical validity/noise or all-sector stability claim follows.

Freeze this argument and the complete-source constants before execution.
The arithmetic verifies the ordinary bounds and their exact budget role;
no multiplier, transport, cold-column or polynomial-degree refinement is
performed. If the omission fits, stop optimizing it. The next decisive
calculation remains the full S_inner application with its original global
output and a certified complement, against the remaining S allowance.
That inverse/application calculation is PROPOSED, not executed here.

## 5. Executed result and scope

The initial 512-bit receipt is preserved. It gives

    K_eta <1.603568,
    K_T <1.217898e-27,
    K_alpha <2.458167e-36,
    R_outer <6.103621e-36,
    a_outer <6.116646e-36,
    debit_outer <8.202157e-36.

The main mixed term dominates. The FULL interior York-return contribution
to K_alpha is <2.523391e-50; it is small only after the weighted global
calculation (1)-(3). This is not a formal local cancellation. The actual
S-omission uses less than .012234 percent of its sufficient 5e-32 screen.
The remaining inner-application ceiling is 4.9993883354e-32; it is not an
attained inner error. The final receipt also checks the K-caller comparison.

What changes: the entire outer-input sector can be omitted at this
certified cost while retaining the full output and operator domain.
This is partial progress in a complete application-error budget, not a
certified numerical S/K solve. The principal inverse/complement estimate
remains unknown. No independent boundary-phase output or physical margin
is obtained. Do not further optimize this now affordable omission term.

Remaining to-do list: full S_inner/K application and complementary-space
errors; complete residual and signed G_t; all original coupled/admission
and physical/stability requirements; the full stable-wormhole objective.
