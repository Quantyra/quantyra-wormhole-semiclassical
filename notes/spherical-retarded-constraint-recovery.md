# Exact spherical constraint transport and retarded recovery

2026-10-20, E003/S168. Ordinary argument before focused verification.
The full research goal remains incomplete.

The gravitational Ward block on a static spherical magnetic exterior
reduces exactly to two free transport equations after weighting the
equation residuals. With zero initial constraints and homogeneous
contractive incoming boundary data, the trace and orbit-trace equations
recover every spherical tensor equation. The corresponding constraint
resolvent exists throughout Re s>0 and has norm at most 1/Re s.

This is a constraint-recovery theorem, not an inverse for the metric
response. It applies on any exact background satisfying its stated
domain hypotheses, without an AdS approximation or an estimate of the
Dirac kernel. For the off-shell seed, explicit residual forcing and
action-density contacts survive. The actual gravitational boundary
condition and electric-current response are not supplied by this theorem.

## 1. Fixed model, chart and equation residuals

Retain the complete selected action, common subtraction, state and
paired fermion domain from [the causal-column audit](charged-conformal-response-and-constraints.md)
and [thermal response](thermal-shape-response-and-subtractions.md).
The [all-frequency normal inverse](curvature-free-all-frequency-bound.md)
is separate; it is not invoked as a causal inverse here.

On either zero-source exterior choose tortoise coordinate x with outer
endpoint x=0 and future horizon x=-infinity. Write

    ds^2=N(x)^2(-dt^2+dx^2)+r(x)^2 dOmega^2,
    n=log N, l=log r, U=N^2 r^2.

The result also holds on a regular subinterval with characteristic data
at its ends. Use the spherical conformal chart of the preceding notes:
delta Phi=p, delta V=v, so h^t_t=h^x_x=2(p+v),
h^theta_theta=h^phi_phi=2p, h_tx=0. Global availability of this gauge
with physical boundary data is not asserted merely from local coordinates.

Let calE be the full equation tensor defined in the causal-column note.
Its stationary mixed components are b_t,b_r,b_o,b_o. Set

    e^a_b=delta(calE^a_b),
    e_t=e^t_t, e_r=e^x_x, e_o=e^theta_theta,
    k=e^x_t=-e^t_x,
    tau=e_t+e_r+2e_o, o=e_t+e_r,
    E_c=U e_r, J_c=U k.                                 (1)

Here tau is minus that note's equation-trace variation. The symbols E_c
and J_c denote constraint residuals, not curvature energy or electric
current. Off shell the mixed variation includes the variation of the
inverse metric. The equality for k holds in the displayed diagonal
chart; an arbitrary off-shell h_tx would require its additional contact.

## 2. The magnetic Ward identity and the off-shell forcing

The vectorlike Dirac Ward identity and Maxwell identity give
div calE proportional to F_ba times the Maxwell/current residual.
At fixed purely magnetic F=(q/2)sin(theta)dtheta wedge dphi, its orbit
components vanish identically for any spherical metric. Consequently
the t,x Ward identities hold even before the orbit Maxwell equations
are solved. This does NOT prove the latter equations. If an electric
perturbation is included, the same linear statement holds when the
background Maxwell/current residual vanishes, as assumed for the
selected neutral stationary state; otherwise retain its force source.

On a nonzero stationary metric residual, differentiating the Ward
identity gives

    nabla_a e^a_b=W_b,
    W_b=sum_a (b_a-b_b) delta Gamma^a_ab,                 (2)
    W_t=(b_r-b_t)(p+v)_t+2(b_o-b_t)p_t,
    W_x=(b_t-b_r)(p+v)_x+2(b_o-b_r)p_x.

The sum in (2) includes both angular components. It follows directly
by varying nabla_a calE^a_b=0; the two connection variations have
opposite signs. Neither a small residual nor spherical symmetry makes
W vanish. On an exact solution b_t=b_r=b_o=0, W=0.

This use of complete renormalized response, its local contacts and
independent state sources is consistent with
[Anderson, Molina-Paris and Mottola, section I, equations (2)-(3)](https://arxiv.org/html/0907.0823v2).
That paper treats de Sitter response; it does not supply the charged
boundary condition or the transport theorem proved below. No removal
of short-scale branches from that source is adopted.

## 3. Exact reduction, including nonzero reduced rows

Direct four-dimensional covariant differentiation in this chart yields

    (e_t)_t+k_x+(2n'+2l')k=W_t,
    -k_t+(e_r)_x+n'(e_r-e_t)+2l'(e_r-e_o)=W_x.          (3)

Primes mean x derivatives. Substitute e_t=o-e_r,
e_o=(tau-o)/2 and U'/U=2n'+2l'. All geometric coefficients in the
homogeneous constraint block cancel, leaving

    (E_c)_t-(J_c)_x=F_c, F_c=U(o_t-W_t),
    (J_c)_t-(E_c)_x=G_c,
    G_c=-U[W_x+(n'-l')o+l' tau].                         (4)

These identities are exact on the charged exterior, not normal-limit
equations. They reconstruct the residual tensor from the two reduced
rows and its initial/boundary constraint data. In particular, on an
exact solution with tau=o=0, they are the free wave first-order system.
With C_+=E_c+J_c and C_-=E_c-J_c,

    (partial_t-partial_x)C_+=F_c+G_c,
    (partial_t+partial_x)C_-=F_c-G_c.                    (5)

C_+ travels into the exterior from x=0, while C_- travels toward x=0.
This fixes which datum is incoming and prevents an implicit reversal
of the horizon condition.

## 4. Retarded recovery and the full right-half-plane Ward resolvent

First impose, as an explicit constraint boundary condition,

    C_+(t,0)=rho C_-(t,0), |rho|<=1,                    (6)

with zero initial E_c,J_c and no incoming data from the past horizon.
For real-time real fields rho is real. A constant complex rho also
defines the complexified spectral problem. Reflection J_c(0)=0 is
rho=1; C_+(0)=0 is absorption. This family illustrates exactly the
needed constraint-domain hypothesis. It is not deduced from the
fermion projector or selected as the physical metric boundary condition.

The norm below is the unweighted L2(dx) norm of the pair (E_c,J_c).
For smooth finite-energy fields, integration by parts gives

    (1/2)partial_t ||(E_c,J_c)||^2
       =Re(conj(E_c)J_c)(0)+Re<(E_c,J_c),(F_c,G_c)>,
    Re(conj(E_c)J_c)(0)=(|C_+(0)|^2-|C_-(0)|^2)/4<=0.   (7)

Thus ||constraints(t)||<=||constraints(0)||+
integral_0^t ||(F_c,G_c)(u)||du. One may first prove this with
compact smooth sources, then extend by the characteristic formulas
and L2 approximation. Zero rows, zero residual background and zero
constraint data imply all e_t,e_r,e_o,k vanish. This is the promised
time-dependent recovery; static bolt regularity alone is unnecessary
and would not suffice in its place.

More explicitly, set gamma=Re s>0 and use exp(s t) modes or the
zero-initial-data Laplace transform. For f_+=F_c+G_c and f_-=F_c-G_c,

    C_-(x)=integral_-infinity^x exp[-s(x-y)] f_-(y)dy,
    C_+(x)=rho exp(s x) C_-(0)
                  +integral_x^0 exp[s(x-y)] f_+(y)dy.   (8)

The expressions belong to H1 for L2 sources, obey (6), and solve
(s-partial_x)C_+=f_+, (s+partial_x)C_-=f_-.
The horizon homogeneous solution exp(-s x) is excluded by L2.
The remaining homogeneous term is fixed by (6), proving uniqueness.
Apply (7) to the spectral equations to obtain

    ||(E_c,J_c)|| <= ||(F_c,G_c)||/gamma.                (9)

This proves existence and the resolvent bound on the closed domain
H1 pairs with (6). For a boundary defect C_+(0)-rho C_-(0)=b, add
C_+=b exp(s x), C_-=0. Its pair norm is exactly

    |b|/(2 sqrt(gamma)).                                (10)

The triangle inequality adds (10) to the right side of (9).
No norm bound for the metric, the stress kernel or the fast scalar
follows from (9). There are no right-half-plane poles of THIS Ward
transport realization, not a proof excluding physical metric poles.
The weighted residual norm and boundary trace still have to be shown
admissible for any proposed physical metric mode.

Boundary data cannot be omitted. On any such geometry,
C_+=b(t+x), C_-=0 solves the homogeneous reduced-row Ward problem.
The corresponding covariant tensor is

    delta calE=(b(t+x)/(2r^2))(dt+dx)^2,                (11)

with zero trace and orbit trace. For smooth b it is regular at the
future horizon in the ingoing coordinate t+x, yet nonzero. This is a
conserved tensor counterexample, not a constructed metric response.
For a causal b it is also a zero-initial-data boundary injection.
It shows why future-horizon regularity and two rows alone cannot
justify recovery. A reflecting matter boundary is not by itself an
assertion about (11)'s gravitational equation residual.

## 5. Off-shell coordinate Hessians retain another contact

There is a second off-shell distinction besides W. Let
tau_0=b_t+b_r+2b_o and o_0=b_t+b_r. Up to the common action
normalization and sphere factor, the chart gradients are 2U tau_0
and 2U o_0. Since delta log U=d=4p+2v, their linearizations divided
by 2U are

    row_a=tau+d tau_0, row_b=o+d o_0.                   (12)

Therefore the source terms for a chart Hessian calculation must be
read as

    F_c=U[partial_t row_b-o_0 d_t-W_t],
    G_c=-U[W_x+(n'-l')(row_b-d o_0)
                         +l'(row_a-d tau_0)].           (13)

For homogeneous chart Hessian rows row_a=row_b=0 these sources generally
remain nonzero. A finite off-shell seed kernel therefore cannot be
called a constrained mode of an exact solution. Both (2) and (12)
must be retained before any limit or quantitative residual estimate.
The cancellation of leading residual terms in the previous AdS normal
limit does not make the finite seed on shell.

## 6. Separate electric equations and the next physical step

For clarity the unresolved spherical Maxwell subsystem can also be
written without a gauge choice for its potential. Let f=delta F_tx,
keep magnetic q fixed, and define the electric flux

    Q_e=r^2 f/(g_c^2 N^2).

For the convention nabla_a F^ab/g_c^2=j^b its linearized equations are

    partial_x Q_e=U delta j^t,
    -partial_t Q_e=U delta j^x.                         (14)

They follow by differentiating the Maxwell density. Metric variations
alone produce no orbit Maxwell density from the magnetic background.
Integrability is partial_t(U delta j^t)+partial_x(U delta j^x)=0.
Also the Maxwell stress has zero first variation under a purely
electric f about a purely magnetic background: the electric/magnetic
cross contractions vanish. Neither fact sets the Dirac current to zero.
Its metric/current response, any electric feedback and its charge data
must be computed or excluded by a demonstrated symmetry of the actual
state and domain. No charge-conjugation symmetry of the magnetic state
is assumed without its bundle and boundary intertwiner.

The next constrained transfer now has a sharper target: solve the
actual two-row retarded metric problem with the residual contacts,
establish its admissible gravitational boundary and initial data, and
check (14). On an exact background and the domain (6), no additional
bulk tensor equation or Ward pole remains to be solved. The actual
state, physical boundary realization, variable-curvature mixing and
existence or justified exclusion of the fast branch remain open.
The finite-frequency normal inverse alone does not close these steps.

The checker independently computes four-dimensional Christoffel
symbols and divergences, the off-shell connection variation, density
contacts, magnetic/electric contractions and characteristic kernels.
Independent complex-frequency quadrature fixtures compare (8) against
closed exponential-source solutions and test (9)-(10). They do not
compute a metric spectrum. The continuum theorem is the ordinary
argument above; no Lean run or formal certification is claimed.

Remaining to-do list: actual two-row retarded stability transfer with
physical initial/boundary and current data; S167 signed coefficient,
cold/higher response, nonlinear solve and compatibility; exterior/state,
charge/return, opening, stability, transport, preparation, separate
creation and full audit.
