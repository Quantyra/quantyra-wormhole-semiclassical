# Moving reflection gluing and the angular boundary load

2026-10-11, E003/S053. Reflection gluing preserves the leading untwisted
return at a frozen orbital phase. It is not time independent during an
orbit. Its angular velocity is twice the orbital velocity; in corotating
mouth charts this becomes an odd boundary angular potential. Including
that datum changes the specified classical junction's rotation parameter
from order e to order P Omega, larger by d/P than the electric-only result.

The projected classical source is still small over the proposed damping
time in the retained joint limit. This closes a necessary boundary-data
check; it does not construct the physical quantum state or its stress.
Known frame transport and gauge covariance underlie the result. No new
general wormhole theorem or consequential novelty is accepted.

## Moving angular identification and coordinate changes

Use the leading circular orbit, n(t)=(cos(Omega t),sin(Omega t),0),
orbital normal l=e_z and R(t)=rotation about l by Omega t. The reflected
identification of [the compatibility calculation](gluing-support-compatibility.md)
has

```text
F(t)=I-2n(t)n(t)^T=R(t) F0 R(t)^T,       F0=diag(-1,1,1),
dot(F) F^T=2 Omega L_l,                 L_l x=l cross x,
F(t)F0=R_l(2 Omega t).                                         (1)
```

Thus a constant F0 in fixed, inertially oriented mouth labels does not
remain the instantaneous exterior reflection. In the spin-j LLL space,
j=(q-1)/2, the relative proper lift can be chosen as

```text
U(t)=exp[-2i Omega t J_l],
||U(t)-I|| <= min(2,2j |Omega t|).                             (2)
```

The sign is an active-rotation convention. A common lift sign is separate
from this relative generator. These are kinematic identities for sphere
labels, not the full retarded fermion return along the orbit.

One can make F constant by using corotating axes at each mouth. This
also changes both the metric connection and the Maxwell time component.
In the magnetic RN convention

```text
ds_ang^2=r^2[dtheta^2+sin^2(theta)(dphi-w dt)^2],
A=-cos(theta)dphi+V cos(theta)dt,
phi_inertial=phi_rotating+Omega t
    => w_rot=w_inertial-Omega, V_rot=V_inertial-Omega.           (3)
```

Here the rotation axis is the orbital normal. The invariant V-w is
unchanged. In the negative-mouth chart the monopole sign is opposite;
after pulling back through F0 the common magnetic orientation is restored.
Axial rotation coefficients transform by det(F0)F0^T, which reverses l.
Hence inertially fixed outer axes give common-chart shifts -Omega and
+Omega at the upper and lower ends. The corresponding Maxwell gauge
shifts have the same parity. Keeping V unchanged in (3) would create a
spurious physical electric field. The full Dirac convention was checked
in [the rotating-field audit](rotating-dirac-and-state.md).

The boundary load below uses these corotating charts, constant F0, zero
radial angular connection and odd core potentials. It explicitly changes
the frozen auxiliary boundary condition omega_out=0 to omega_out=-P Omega
on the upper RN annulus. It is a specified stationary comparison to the
moving problem. Coordinate rotation by itself does not create a physical
spin: the relative seam, radial connection and boundary frames must all
be specified before interpreting a solved rotation charge.

## What actually determines a changing return holonomy

For clarity, let a Hermitian matrix connection have D_mu=partial_mu+iA_mu,
F_tx=partial_t A_x-partial_x A_t+i[A_t,A_x], and transport
partial_x W(x,a)=-i A_x W(x,a), W(a,a)=I. Direct differentiation gives

```text
partial_t W(b,a)+i A_t(b)W(b,a)-i W(b,a)A_t(a)
 = -i integral_a^b W(b,x)F_tx(x)W(x,a) dx.                      (4)
```

To prove it, differentiate W^-1 partial_t W in x, integrate, and
integrate the partial_x A_t term by parts; the commutator supplies the
derivatives of W and W^-1. Endpoint terms give the two A_t terms.
Equation (4) is covariant under changes of the frames at both endpoints.
For a time-dependent seam psi_-=U psi_+, the transformed Hamiltonian
contains i dot(U) U^-1 in addition to U H_+ U^-1. Neither term can be
silently discarded when expressing the same evolution in another frame.

This explains the limitation of (1): the angular identification's rate
alone is not the rate of gauge-invariant loop eigenphases. Parallel
transport, electromagnetic phases and the physical seam lift all enter.
A compensating radial connection can move time dependence between these
descriptions. It must be transformed along with the boundary data, rather
than being set to zero in both descriptions without justification.

## Transported mouth axes

The axes of an instantaneous rotation-free Lorentz boost are not
Fermi-Walker axes. For a prescribed circular orbit in Minkowski space,
speed v, laboratory time t and gamma=(1-v^2)^(-1/2), take signature -+++.
The equations for a spatial tetrad vector S are

```text
u=gamma(1,v tangent),     a=gamma^2(0,-v Omega n),
dS/dt=[u (a.S)-a (u.S)]/gamma,     u.S=0.                       (5)
```

If B(v(t)) is the standard rotation-free boost, substitution gives

```text
S_i(t)=B(v(t)) [0,R_l(theta_T(t))e_i],
theta_T(t)=-(gamma-1)Omega t.                                   (6)
```

The boost has B00=gamma, B0i=Bi0=gamma v_i and
Bij=delta_ij+gamma^2 v_i v_j/(gamma+1). Differentiating this matrix and
using v dot d_t(v)=0 verifies (5) for each column, while preserving
S_i.S_j=delta_ij and u.S_i=0. This supplies the ordinary derivation of
(6). [Bini et al., v1, sections 3.2, 4.1 and the Minkowski specialization](https://arxiv.org/html/1408.4598v1)
establish the frame relations and Thomas precession; this is prior work.

For the retained leading orbit v^2=P/d=chi_b, the correction to a
rotation-free boosted frame is O(chi_b Omega). Companion curvature and
gravitomagnetic transport also enter at that relative order. Therefore
(6) verifies a flat-space comparison, not the complete binary precession
coefficient. Those omitted terms do not replace the leading order-Omega
boundary datum. No strong-field transported tetrad through the throat is
inferred from this flat calculation.

## The prescribed classical junction with moving-frame data

Retain the RN/product equations and junction of the previous compatibility
note. Let lambda=P/ell, throat endpoints +/-R, s=lambda R, r_s=1+s,
r_out>=2, and omega_hat=P Omega. The upper electric coefficient is
e=P E=2v(P/d)^2, so e=chi_b omega_hat at leading binary order.

With the prior u0,u1,f and H_o=rho(1+rho^2), g_o=1+3rho^2,

```text
u_RN=e u0+b u1-2J/r,
V_th=k H_o-2J rho,       w_th=-2k H_o-2J rho,
b u1-k g_o=-e u0+2J(1/r_s-1),
b f u1'-lambda(1+R^2)k g_o'=-e f u0'-2J f/r_s^2.               (7)
```

All quantities in the last two equations are evaluated at the upper
interface. Their determinant is nonzero. Let (b_E,k_E) and (b_J,k_J)
be the solutions for unit e and unit J. Reconstruct

```text
w_out=lambda[-2k H_o(R)-2JR]
        +integral_(r_s)^r_out (-6J-2u_RN)/r^4 dr
     =e w_E+J w_J.
```

The moving-frame condition is now w_out=-omega_hat. Thus

```text
J=(-omega_hat-e w_E)/w_J,
J/omega_hat -> [1-chi_b(1-w0(r_out))]/[1-wKN(r_out)],
w0=2/r-(2r-1)/r^4,     wKN=(2r-1)/r^4,                        (8)
```

as s->0 and R->infinity. The denominator is bounded away from zero for
r_out>=2, and the overlap remainders are uniform as r_out increases:
the decaying-mode integrals converge, while the growing mode's w integral
has an integrable tail. For large r_out and small chi_b, J~P Omega.
The earlier fixed common-frame electric result instead had J~ -e.
That difference is a changed boundary problem, not a sign correction to
the earlier calculation.

The fields remain small in the following precise asymptotic sense. For
s sufficiently small, R sufficiently large and r_out>=2, (7) gives

```text
k_E=O(lambda^2),   b_E=O(s^5(s+R^-2)),
k_J=O(lambda/R),   b_J=O(s^4).
```

Using integral_(r_s)^infinity u1/r^2=O(s^-2), the reconstruction of V
and w then gives, with a constant independent of these matching scales,

```text
sup_RN |w| <= C(|J|+|e|),
sup_RN |V| <= C(|J|+|e| r_out),
lambda sup_throat (|V_th|+|w_th|) <= C(|J|s+|e|s^3).            (9)
```

For the last estimate, bound |H_o| by C R^3 and use lambda R=s;
the spin k term contributes O(J s^2), below O(J s). The first two
follow by adding the regular antiderivatives to the decaying-mode
integrals; the only growing potential is e r. These arguments prove
asymptotic uniform constants, not explicit finite-q thresholds.

At a physical port a=p d, r_out=a/P=p/chi_b. Thus e r_out=p omega_hat.
For p bounded by a fixed small constant and chi_b small, (8)-(9) bound the
projected classical multiplication source throughout throat and annuli:

```text
A_eff=-J_l (w+q V)/(q+1),
||A_eff,physical|| <= C j |Omega|.                             (10)
```

RN time is t/P and throat time is t/ell, so the respective conversions
are division by P and ell; (9) supplies exactly these factors. Equation
(10) concerns the free LLL vertex already derived in the current audit.
It does not control omitted Landau levels, the complete exterior
derivative operator or the quantum source's self-consistent response.

## Finite-time consequence and the unresolved quantum sector

Suppose a well-posed projected one-particle problem with a common domain
has Hamiltonians H0(t) and H0(t)+B(t), B bounded and self-adjoint. For
their unitary propagators and identical initial covariance 0<=C0<=I,
Duhamel's identity and unitarity give

```text
||U_B(t)-U_0(t)|| <= min(2, integral_0^t ||B(s)|| ds),
||C_B(t)-C_0(t)|| <= 2 min(2, integral_0^t ||B(s)|| ds).          (11)
```

The second formula follows by expanding U_B C0 U_B^*-U_0 C0 U_0^*.
A trace-class observable A has error at most ||A||_1 times this bound.
This controls bounded smeared finite-band observables once their domain,
state and actual source have been supplied. It cannot be taken through
the coincident-point derivatives and ultraviolet subtraction defining
the renormalized stress without further estimates.

Apply (11) only to the bounded component (10). With
g_c^2=q^-2, lambda=(16pi q^3)^-1, delta=d/ell=log^-2 q and
the proposed damping time t_d=ell/delta,

```text
Omega ell=(log q)^3/(2 sqrt(pi) q^(3/2)),
q Omega t_d=(log q)^5/(2 sqrt(pi) sqrt(q)) -> 0.                (12)
```

Thus this newly retained classical angular load does not invalidate the
finite-time *bounded-source* hierarchy. Convergence is slow; modest-flux
values need not be small. This is not an adiabatic theorem for every
channel. The point-dipole model has D_max of order qd, so its smallest
free-circle gap is of order 1/(qd). Comparing this gap to the crude
j Omega norm gives a ratio of order q^2 Omega d, growing like sqrt(q)
log q. That sufficient uniform-gap argument fails even though (12)
vanishes. It is not evidence of unavoidable excitation or instability.

The spin normalization exposes the remaining state issue. With K=P^2/G
and ell=16P^3/(Gq), the classical solution has

```text
KJ ~ K P Omega = q Omega ell/16 -> 0.                          (13)
```

This small continuous classical value does not specify a quantum spin
eigenstate, ensemble or boundary charge. The actual free-current Gauss
law conserves Q=K J_geom+sqrt(T_q/pi) chi_current, not K J_geom alone.
The spatial angular Wilson line, its compact sectors and the three
noncommuting angular currents must be treated in the same exterior/throat
state. A local positive quadratic energy does not determine them.
The global corotating ground state remains excluded by the already
reviewed light-cylinder argument; this work does not restore that option.

## Verification and research decision

[The checker](../checks/moving_reflection_load.py) verifies (1)-(3),
integrates the original four-vector Fermi-Walker equation at two
tolerances, and checks (4) on noncommuting matrix connections by separate
transport, time-variation and curvature integrals. Fifteen classical
junction cases, each at two tolerances, solve the original RN flux and
both potential equations with J unknown. The invariant V-w is checked
against its independent flux-derivative reconstruction.
[The receipt](../evidence/moving-reflection-load.json) pins exact inputs;
[the source record](../evidence/moving-reflection-source-review.json)
records editions and review boundaries. No sampled calculation replaces
the ordinary estimates or tests the physical interacting state.

Continue the physical long-return/finite-transfer program with the
moving-frame load retained. The next decisive calculation is the coupled
finite-time state and angular holonomy response, including the global
rotation sector. Frame algebra, Thomas precession and Duhamel continuity
are established ingredients, not a novel endpoint. No Lean or cloud
campaign, publication or outreach is undertaken for this component.

Remaining to-do list: determine the coupled angular/electric quantum state
and physical return, control renormalized stress and omitted response,
complete smooth semiclassical matching and the consequential endpoint,
refresh novelty, independently verify, GCP-certify suitable mathematics
and prepare the reproducible manuscript.
