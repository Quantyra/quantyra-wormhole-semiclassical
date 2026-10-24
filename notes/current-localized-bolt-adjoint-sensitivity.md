# The localized physical bolt functional after exact scalar elimination

2026-10-24, E003/S169. Ordinary argument before acceptance execution.
Previous goal turn: progress. The full stable-wormhole goal stays active
and incomplete. Preserve the closed audit and every existing physical
source, regulator, contact, endpoint and normalization.

The broad xi=z^4 screening trial does not upper-bound the adjoint norm.
Test that norm with a genuine TT direction localized at the actual bolt.
This is an adjoint sensitivity calculation, not refinement of the already
controlled physical cold source or scalar correction. It determines how
much the chosen uniform .0022 comparison could lose in a correction
certificate. No finite-mode physical problem replaces the actual operator.

## 1. A genuine localized test on the original finite seed

Use Q,z_h,h,gamma and the fixed York return of
current-physical-york-dual-assembly.md. Fix a=1e-4 and set

    t=(z_h-z)/h, v=t/a, z_a=z_h-ha,
    q(t)=Q(z_h-ht)/h^2,
    psi(t)=(1-t/a)^8 for 0<=t<=a, zero for t>=a,
    J=psi/h,
    A(t)=q''(t)/z^4+4h q'(t)/z^5,
    K(t)=integral_t^a A(s)psi(s)ds,
    xi=-(q psi)'-z^4 K,
    S=xi/2+h q psi/z,
    T=-xi+(q psi)'/2.                              (1)

Primes in this note mean t derivatives except where specified otherwise.
All fields are zero for z<z_a. In z coordinates,
K=integral_0^z (Q'/z^4)'J and xi=(QJ)'-z^4K, exactly as in the
preceding genuine-TT construction. Hence the original TT constraint,
J(0)=0 and (Q^2 J'/z^4)(z_h)=0 hold. The matched field is C6, belongs
to the existing C5 phase and is regular in bolt Cartesian coordinates:
t is a smooth function of the squared disk radius there. The nonzero
values at the bolt are allowed. It is not a puncture boundary condition.

Retain the exact thermal identities q(0)=0, q'(0)=2gamma/h. Expand
the remaining polynomial from the actual Q derivatives at z_h; never
estimate Q(z_h) by subtracting independent intervals. With A>=0,

    K>=0, K<=a sup A/9,
    xi(z_h)=-2gamma/h-z_h^4 K(0),
    |xi(z_h)|/2>=gamma/h.                           (2)

The test changes neither the reference metric nor the prescribed quantum
stress. Its amplitude is arbitrary because only a homogeneous norm ratio
will be used. The actual correction is not this test direction.

## 2. A global Poisson barrier for a compact cold scalar source

The scalar Schur return cannot be omitted. Let u=A_tr^-1(lambda g)
for the ACTUAL complete mixed row g induced by (1), as in
current-schur-bolt-functional.md. It is supported on 0<=t<=a.
To control u at the bolt, keep the exact scalar factorization

    A_tr=6 L_D L_4+Q_rem, L_D=-D Delta+C_v, L_4=-Delta+4lambda,
    L_D w>=m w, L_4 w>=.63lambda w,
    m=1-1.1e-8, ||(1/6)L_D^-1 Q_rem L_4^-1||w<epsilon=1e-26.

These are the existing actual finite-geometry estimates with w increasing
outward. No scalar Galerkin inverse or product geometry is substituted.
L_D and L_4 are resolvents of the SAME regular-bolt/zero-source Laplacian;
they commute, first on its spectral realization and then on their bounded
Green extension. This does not commute the perturbation Q_rem.

Here is an explicit whole-exterior Green upper bound. Write
P(z)=beta(1+2z+3z^2)+z^2/C+(d/C)z^4. The actual d/C is negative.
P>=beta and (P/z^3)'<0 on 0<z<=z_h. Using the exact root relation,

    Q(z) >= beta[(1-z)^2-(1-z_h)^2],
    q(t) >= beta t(t+kappa), kappa=2/(1+h).          (3)

This is a finite-temperature inequality on the entire geometry. In these
variables -Delta/(C/B^2)=-z^4 partial_t(q z^-2 partial_t).
The nonnegative Poisson solution to source 1_[0,a], regular at the
bolt and zero at infinity, has bolt value

    V_P=integral_0^(z_h/h) (z^2/q)
                           integral_0^min(t,a) z(s)^-4 ds dt
       <= [z_h^2/(beta z_a^4)] I(a),
    I(a)=log(1+a/kappa)+(a/kappa)log(1+kappa/a).     (4)

The integral formula follows by integrating the divergence equation from
the regular bolt. Its solution decreases outward. A positive mass only
decreases the Green solution by the maximum principle, and dividing by
the outward-increasing w puts the largest ratio at the bolt. The positive
Green/weighted bounds and the order of the perturbation series therefore
give, for |g|<=g_max on its true support,

    |u(z_h)| <= g_max [ (x/C) V_P/(6m)
                         +epsilon/(6*.63*m*(1-epsilon)) ]=:U_s. (5)

For the second term put y=L_4 u, expand
y=(1+(1/6)L_D^-1 Q_rem L_4^-1)^-1 (L_D^-1 lambda g)/6,
and then apply L_4^-1 to the remainder. The w(z_h) factors cancel
because the source support has w>=w(z_h). This proves (5) for the
actual scalar realization, including its allowed tail at the outer end.

## 3. Compute the complete mixed source before taking absolute values

Use the EXACT canonical trace-row variation of
current-finite-local-action-assembly.md with (nu,u_shape)=(S,T):

    g=-(C/x) z^4 delta(F_E+F_R2+F_W), x=lambda B^2.  (6)

F_W includes the complete Weyl/magnetic anomaly, and F_R2 the selected
finite curvature coupling. The canonical determinant has zero scalar
and mixed row; it is not approximated by omitting these rows. Off-shell
density and frozen-return contacts are in (6).

Compute with Delta=(q partial_v)'/a^2 and phi_v=ha/z at the seed.
Represent the first variation exactly as rational polynomials in v,z,K,
with z=z_h-ha v and
partial_v K=-a(q''/z^4+4h q'/z^5)psi. K enters linearly. Algebraically
combine the resulting source into

    g=z^-M [P_g(v)+K(v) R_g(v)]

before applying Bernstein bounds on [0,1]. This retains the TT and
curvature cancellations; bounding individual fourth derivatives first
would lose them. Bound the remaining K by (2), and z^-M by z_a^-M.
Carry every parameter interval and every zero-containing coefficient.
Freeze a diagnostic target g_max<100, U_s<.02. Then

    |ell_eff(xi)|/ell_b
       =|xi(z_h)/2+u(z_h)| >= gamma/h-U_s=:b_low>0,
    ell_b=(1+h)/h.                                  (7)

This is an actual scalar return estimate. A positive b_low demonstrates
noncancellation on this test; it does not evaluate the full point adjoint.

## 4. Bound the complete test energy, retaining the quantum spectrum

Use the full known affine local action, a0 interval, a2 and LLL in the
preceding weak-to-admission audit. Positive product bounds on the actual
S,T jets through order two control their exact diagonal. Use t derivatives,
the actual q and phi_t=h/z; integrate over the true support dz=h dt.
There is no factor h^-1 left from changing derivative coordinates.

The complete positive massive geodesic term obeys

    P_geo/|W| <= [2/(525pi B^2)] [4C^2/q_mag] Z,
    Z=h integral_0^a q[(5/8)q^2(T''')^2+(3/2)q T''' U+U^2]dt,
    U=(3/2)q'T''+(1/2)q''T'.                        (8)

Here q_mag=2000^4. Formula (8) is the exact coordinate conversion of
the full sixth-moment bound, not a truncation of angular, thermal or
flow frequencies. All fixed-scale contacts and the full finite-geometry
remainder enter via the unchanged TT-only comparison error eta=.0022.
If U_H is the resulting normalized absolute comparison diagonal,

    l(xi,xi)/|W| <= L_c=U_H/(1-eta),
    ||ell_eff||_(l*) >= ell_b b_low/sqrt(|W| L_c),
    ||Sc^-1 ell_eff||_l >= ell_b b_low/[M_S sqrt(|W| L_c)],
    M_S=1.341.                                      (9)

The physical normalization |W|=N pi B^2/gamma is retained for N=1,2.
Freeze the conservative diagnostic target for the last quantity >1e8.
These are LOWER bounds. No adjoint approximation or upper norm is claimed.

## 5. Recalculate the necessary projection precision

The old global xi_g=z^4 test is still valid, with its retained upper
L_g=l(xi_g)/|W|. Use DIFFERENT tests for the primal and dual norms:

    R0= sqrt(|W|) |I_g|/sqrt(L_g),
    L0=ell_b b_low/sqrt(|W| L_c), I_g=r_eff(xi_g)/W.

The all-trial minimization in
current-uniform-bolt-certificate-obstruction.md applies to these two
independent norm lower bounds. Every trial pair's uniform remainder
expression is at least eta(1+eta) R0 L0/M_S^2. Thus even giving that
remainder the ENTIRE 1.667e-9 bolt budget requires

    |I_g| < T_new,
    T_new=1.667e-9 M_S^2 sqrt(L_g L_c)
                            /[eta(1+eta) ell_b b_low].            (10)

This is necessary ONLY for that sufficient uniform-error certificate.
It is not necessary for a controlled wormhole. The earlier 2.60e-40
screen and 1.287e-40 quantum/contact allocation remain valid weaker
screens; passing them did not prove the full certificate would work.
Compare the actually certified reduced-source allowance 1.3079691e-40
and rest error 2.123e-42 to (10). Never compare them as measured nonzero
values: their signed centers are zero and the sign is still unknown.

Freeze T_new<1e-52 as the sensitivity decision target. If passed, the
uniform .0022 proof route needs at least twelve more orders of projected
source precision than its present complete allowance. This is a known
insufficient upper certificate, not a proved nonzero physical error or
a proof that the actual operator difference saturates .0022.

For any positive evaluated |I_g|, the necessary uniform eta would obey
eta(1+eta)<1.667e-9 M_S^2 sqrt(L_g L_c)/(ell_b b_low |I_g|).
Substituting an upper allowance for |I_g| gives a hypothetical design
target only, never a necessary eta for the unmeasured actual source.

DECISION CONTRACT: if the sensitivity target passes, prioritize a
directional/weighted bound on the actual TT operator remainder, retaining
the cold suppression and the exact scalar Schur block, over further
global uniform-norm or already passing source-screen optimization. The
test alone does not rule out the .0022 certificate for the actual source.
Do not claim a solved correction, new source sign, nonlinear admission or
stability. Update the checkpoint with the evaluated numerical threshold
and explicitly identify the actual calculation still unexecuted.

Remaining to-do list: actual coupled correction with directional operator
error and strong weighted return; X-ball admission and nonlinear source
remainder; exact compatibility; physical matching and causal realization;
semiclassical validity/noise and all-sector stability; charge/return,
opening, transport, preparation, separate creation and stable wormholes.
