# Actual global scalar return in the reduced-source projection

2026-10-23, E003/S169. Ordinary argument before frozen verification.
Previous turn: progress, checkpoint reconciliation completed and pushed.
The full stable-wormhole objective remains active and incomplete. The
feasibility audit stays closed. This calculates the scalar component of
the selected reduced-source test, not the independent shape component.

## Exact quantity and a required distinction

For xi=z^4,J=0 the fixed York return is S=z^4/2,T=-z^4. Let
u=A^-1 B xi for the ACTUAL complete positive scalar block A. Then

    r_eff(xi)/W = R_nu(z^4/2+u) - R_b(z^4),
    R_nu(v)=(x/C) integral_0^z_h z^-4 (e_disk/lambda) v dz,
    W=-N pi B^2/(G gamma), z_h=1/(1+h), x=lambda B^2.    (1)

R_b is the complete canonical shape row of the Newton defect r1, not
the seed gradient or its trace anomaly. The preceding point-return
bound determines u(z_h), not the entire function needed in (1).
The scalar part of (1) and R_b are correlated through the same first
reference correction. Separate intervals need not preserve a useful
cancellation. No shape source is set to zero in this calculation.

## Actual scalar form and the global residual certificate

Use the unchanged action, seed, paired domain, fixed York chart and
physical mass/log phase. Put Q=z^2 F/C, R=B^2 R4, tau_c=B^2 tau_Ric,
epsilon=1/(60pi), D_R=96pi c_R, D=D_R-epsilon. For radial v define

    L v=(Q v')'-2Q v'/z,
    c4=6D C/B^2,
    c1=6+(-2D_R+epsilon/6)R/B^2-(22epsilon/3)tau_c/B^2,
    m=(2R+48lambda0 B^2)/C.

The full normalized positive scalar form is exactly

    a_bar(v,w)=integral [c4 L v L w+c1 Q z^-2 v'w'
                                      +m z^-4 vw] dz.              (2)

This is the independently checked scalar form of the local action,
including its anomaly and off-shell density terms. At the actual seed,
c4>0, c1>5.99 and m>3.99. In particular

    a_bar(v,v)>=m_* integral z^-4 v^2 dz, m_*=3.99.                 (3)

The source is f=(x/C)z^-4 g, with the complete mixed polynomial g
retained by current-schur-bolt-functional.json. Its minus sign relative
to the normalized mixed action follows from W<0. For a conforming trial
u_p, compute its COMPLETE strong residual d=A_bar u_p-f. Here

    L* v=(Qv')'+(2Qv/z)',
    A_bar v=c4 L*L v-(c1 Q z^-2 v')'+m z^-4 v.                     (4)

Regularity at the actual bolt and the z^4 outer decay eliminate the
integration-by-parts endpoints. The trial belongs to the scalar form
domain; no boundary equation or zero-mass constraint is replaced. The
actual form inverse and Cauchy--Schwarz give

    ||u-u_p||_(a_bar) <= ||z^2 d||_2/sqrt(m_*),
    ||z^-2(u-u_p)||_2 <= ||z^2 d||_2/m_*.                           (5)

This is an infinite-domain residual certificate. Matrix convergence
alone, or equality of a finite set of moments, would not imply (5).

## Trial choice and exact integration

The exploratory polynomial trials z^4,... converge slowly at the cold
end. The largest tested polynomial trial misses the independently bounded
bolt value by about .075. Its interior samples do not certify a global
inverse. Both the polynomial and fractional exploratory source snapshots
and outputs are preserved as diagnostics, not accepted solution evidence.

At the extremal cold end the light scalar equation has the positive
exponent nu satisfying nu(nu+1)=4x/C. Use a FIXED numerical midpoint
nu=(sqrt(1+16x/C)-1)/2, about .45742710775, only to select trials

    z^(4+j), z^(4+j)(1-z)^nu, j=0,...,10.                          (6)

Every trial is smooth at the ACTUAL positive-temperature bolt, where
1-z_h=h/(1+h)>0. At infinity the fractional factor is analytic and
the same z^4 decay holds. The trial choice does not replace the actual
geometry by its extremal limit. Any exponent-selection error is included
in the complete residual (4); the exponent is not asserted to solve the
finite-temperature equation. No approximation to the fast scalar branch
is declared exact.

For stable arithmetic retain the exact factorization

    Q=(1-z)^2[beta(1+2z+3z^2)+z^2/C+(d/C)z^4]-tau_th z^3,
    beta=x/C,
    tau_th=(1+h)h^2[beta(1+2/(1+h)+3/(1+h)^2)
                      +1/(C(1+h)^2)+(d/C)/(1+h)^4].              (7)

Expanding independent coefficients before multiplying very large cold
derivatives would destroy this correlation. Equation (7) retains it.
Integrate every term z^k(1-z)^(j+t nu), k>=0, by expanding z^k in
y=1-z and taking power primitives from y=h/(1+h) to 1. There is no
quadrature omission or endpoint excision. Negative powers are integrated
at the actual positive lower endpoint, not at zero.

The matrix uses midpoint entries solely to choose a trial. In the
certificate restore all parameter and mixed-source intervals. Split the
residual into its coefficient midpoints and uncertainties and use

    ||z^2 d||_2 <= ||z^2 d_mid||_2
        +sum radius(d_kjt) ||z^(k+2)y^(j+t nu)||_2.                 (8)

This triangle inequality avoids subtracting uncertain large coefficients
inside a squared norm. It changes neither the physical residual nor
its tolerance. The final run must have ||z^2 d||_2<1e-9.

The first frozen 18-dimensional run missed that target with 2.406e-9.
Independent entry midpoints need not be the matrix of one consistent
point operator. A replay selected coefficient midpoints BEFORE assembly,
including one fixed lower integration endpoint, and then restores every
actual interval for the residual. The failed snapshot is retained and
the target is unchanged. That replay gave the same residual to the
displayed accuracy, so matrix midpoint inconsistency did not explain the
loss. Retain its snapshot too and increase the trial dimension to 22.
Acceptance still depends on (4)-(5), not matrix conditioning assumptions.

## Scalar source pairing, including unknown terms and endpoints

Use the exact trace-defect identity and finite-seed polynomials already
read in current-uniform-bolt-certificate-obstruction.md. Write

    tau=T0/lambda, K=d_l L_lambda tau+eta K_t tau,
    e_T/lambda=K+A_f f/(lambda r^2)+A_Q G(rho-p)/(lambda r^2),
    eta=epsilon lambda, d_l=D lambda.

The complete fixed-disk trace is e_disk=e_T+a T0+
r^-2 partial_r(r^2 xi_disp T0). For v=z^4/2+u_p, integration gives

    R_nu(v)=(x/C){integral [z^-4 v(e_T/lambda+a tau)
                          +z^-2 v' (xi_disp/B)tau] dz
                 -v(z_h)z_h^-2(s_R/B)tau(z_h)}.                   (9)

The outer endpoint vanishes by the actual decay. The bolt endpoint is
retained. The previously proved bounds are

    |f/(lambda r^2)|<=2.11e-6, |G(rho-p)/(lambda r^2)|<=2.01e-7,
    |a|<3.3e-10, |s_R/B|<3.33e-10 h,
    |xi_disp/B|<=3.33e-10 h+2(1.08e-9)(1-z)/z.                    (10)

The run verifies v'>=0 on the entire domain by interval Bernstein
subdivision of its two polynomial factors; v(0)=0 then gives v>=0.
On interior cells a degree-six Taylor enclosure of (1-z)^(nu-1),
with its explicit seventh-derivative remainder at the nearer endpoint,
preserves cancellation between the factors before taking their Bernstein
minimum. A direct factor bound handles the final cold cell. This reduces
the number of certified cells without changing the inequality. The
earlier successful run's larger subdivision remains in its receipt.
If P_abs denotes the absolute-coefficient polynomial, use
tau_abs=-S A+(tau-S A)_abs, where A=-10z^4-20z^5+30z^6<=0. Thus
|tau|<=tau_abs without losing the leading exact product zero. Integrating
(10), tau_abs, |A_f|<=A_f,abs and |A_Q|<=A_Q,abs against v and v'
gives explicit symmetric allowances for all unevaluated areal, density
and transport terms in (9). They are upper bounds, not measured values.

For the trial error in (9), use the global trace identity directly. The
retained |a+xi_disp'+2xi_disp/r|<8.7e-5 and (10) give
|xi_disp/(r-B)|<2.493e-9. An explicit positive polynomial majorant is

    E=K_abs+2.11e-6 A_f,abs+2.01e-7 A_Q,abs
           +8.7e-5 tau_coeff_abs+2.493e-9[z(1-z)tau']_abs,
    |e_disk|/lambda<=E,
    C_s=(x/C) sqrt(integral z^-4 E^2 dz).

All these polynomials start at least at z^3, so this integral is finite.
Together with (5), the COMPLETE scalar pairing error is bounded by

    |R_nu(u-u_p)| <= C_s ||z^2 d||_2/m_*.                         (11)

The first 22-dimensional run reached the residual target but rejected
the outer majorant: independently interval-expanded leading curvature
constants left a zero-containing constant coefficient. Its absolute
value would give a nonintegrable allowance. Do not drop that coefficient.
Use the exact identity, with a=4+2/x+delta/x and b=3+1/x+d/x,

    F/x=z^-2+1/x-a z+b z^2-2(d/x)z^3+(d/x)z^4,
    c/lambda=a z^3-2b z^4+(20/3)(d/x)z^5-5(d/x)z^6.

The constant cancels algebraically before evaluation. This retains the
same finite-temperature geometry and physical vacuum cancellation.
The rejected run and its sources are retained; no positive coefficient
or uncertain physical parameter is discarded.

Require (11)<1e-40, so numerical scalar inversion consumes less than
the 2.6035780e-40 entire-budget necessary reduced-source ceiling. This
does not require the unknown signed source to be close to zero. Record
the actual scalar interval and the width of the unevaluated contact
allowances; a wide interval cannot certify cancellation with R_b.

For the known K pairing, integrate each fixed trial moment first and
then multiply by its shared uncertain physical coefficient. Expanding
those coefficient intervals through the ill-conditioned polynomial trial
before integrating gave a valid but much wider allowance in the first
successful run. The replay preserves the same coefficient uncertainty
and exact integral, without that avoidable interval dependency loss.
All post-success scientific changes require and receive a full replay.

## Decision and scope

This executes the GLOBAL scalar part needed by the reduced-source test.
The full signed reduced projection still requires R_b(z^4), with the
quantum Hessian acting on the actual first correction and the same
scalar transport/contact terms. If their uncertainty exceeds the target,
evaluate the correlated joint expression rather than refining an already
small scalar numerical residual. A passing scalar inverse check does
not determine the unknown shape source or certify the uniform TT bound.

The single next calculation remains that full signed reduced-source
projection, now using a controlled global scalar test function. Keep
the common quantum subtractions, differentiated paired domain, mass/log
phase and exact compatibility. Exterior admission, preserved nonlinear
neighborhood, causal realization, matching, semiclassical validity/noise
and all-sector stability remain indispensable. No actual controlled
solution or nonexistence result follows from this partial calculation.

Remaining to-do list: complete signed reduced-source projection and its
correlated quantum/contact errors; actual Schur forward/adjoint solve and
physical bolt/exterior certificate; X/admission, nonlinear closure and
compatibility; matching, causal realization, semiclassical validity/noise,
all-sector stability, charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
