# Global optical gauge and the boundary clock

2026-10-21, E003/S168. A complete global gauge-range calculation exposes
an additional condition in using the two optical rows as a description
of arbitrary spherical perturbations. With an explicit gauge domain
fixing the outer endpoint AND its time parametrization, the gauge map
has a codimension-one range at each right-half-plane frequency. A
compact metric probe need not satisfy that condition. Allowing the
boundary time to transform restores the map, but its optical image can
leave the preceding curvature reconstruction space. This is a domain
decision, not an instability, a physical boundary selection or a new
nonlinear solution. The full goal remains incomplete.

Keep the action, state, flux, temperature and prepared clock of the
[charged shape calculation](charged-curvature-free-shape-form.md).
That calculation and its restricted positivity remain valid. Its
extension to a global physical reduction must include the issue below.
The prescribed boundary clock used here is an EXPLICIT test of one
possible gravitational gauge realization. The earlier work did not
derive a complete physical metric boundary domain; this note does not
silently select one or assert that every kinematically allowed probe
is an on-shell disturbance.

## 1. Exact gauge equations on the entire optical half-line

Use inward optical coordinate y, with the outer end at zero and horizon
at infinity. The background is

    g=F(y)(-dt^2+dy^2)+r(y)^2 dOmega^2,
    r_y=-F, n=(log F)/2, l=log r.

Write a general spherical perturbation as

    h_tt=-2F A, h_yy=2F B, h_ty=F C,
    h_theta_theta=2r^2 p_0,
    h_phi_phi=sin(theta)^2 h_theta_theta.

Here B is a metric amplitude, not the fixed reference radius used in
the seed formulas. For h_new=h-Lie_xi g, direct differentiation gives

    A_new=A-n_y xi^y-partial_t xi^t,
    B_new=B-n_y xi^y-partial_y xi^y,
    C_new=C+partial_y xi^t-partial_t xi^y,
    p=p_0-l_y xi^y.                                      (1)

The optical chart requires A_new=B_new=p+v and C_new=0. Define
d=A-B, zeta_+=xi^t+xi^y, zeta_-=xi^t-xi^y and f_+=d+C,
f_-=d-C. The exact gauge equations are

    (partial_t-partial_y)zeta_+=f_+,
    (partial_t+partial_y)zeta_-=f_-.                      (2)

The background conformal factor cancels, with no field equation or
asymptotic approximation. zeta_+ travels OUT to y=0; zeta_- enters the
exterior there. Preserving the endpoint requires xi^y(t,0)=0, hence
zeta_-(t,0)=zeta_+(t,0). Also preserving its clock requires xi^t(t,0)=0,
hence BOTH traces zero. The latter is an extra condition on the outgoing
solution, not a second freely specifiable incoming datum.

## 2. Exact frequency range and bounds

Let Re gamma=a>0. First state the domain precisely: zeta belongs to
H1((0,infinity);C^2), with the stated traces at zero. This is an optical
gauge-vector norm, not a gravitational energy or a replacement for
future-horizon regularity. Section 5 relates it to the latter for compact
sources at sufficiently large real frequency. Put

    L=diag(-partial_y,+partial_y),
    Theta_gamma(f)=integral_0^infinity exp(-gamma s) f_+(s) ds. (3)

The functional is complex linear, with norm 1/sqrt(2a) on the plus
component. The unique square-integrable outgoing solution is

    zeta_+(y)=integral_y^infinity exp[-gamma(s-y)] f_+(s) ds. (4)

Its trace is Theta_gamma(f). With endpoint preservation the other
component is

    zeta_-(y)=exp(-gamma y)Theta_gamma(f)
               +integral_0^y exp[-gamma(y-s)] f_-(s) ds. (5)

Young's inequality and the differential equations put (4)-(5) in H1
for f in L2. The explicit formulas prove existence and uniqueness.
Integration by parts gives

    Re <zeta,L zeta>=(|zeta_+(0)|^2-|zeta_-(0)|^2)/2.

Thus the endpoint-preserving realization has

    ||(gamma+L)^-1||_(L2->L2)<=1/a.                     (6)

For the smaller domain D_00 with BOTH traces zero, (4) instead proves

    Ran((gamma+L)|D_00)=ker Theta_gamma,
    ker((gamma+L)|D_00)={0}.                            (7)

The range is closed of codimension one, and the graph-domain map has
index -1. Its inverse on that range obeys (6); the derivative bound
||zeta_y||<=||f||+|gamma| ||zeta|| controls the graph norm. A nonzero
compact positive f_+ proves that the range is proper. The theorem is
global, with no artificial horizon wall and no uncomputed kernel.

There is no residual homogeneous L2 gauge after endpoint preservation:
the homogeneous solutions are (a_+ exp(gamma y),a_- exp(-gamma y));
integrability sets a_+=0 and the endpoint relation then sets a_-=0.
The obstruction is coverage of the gauge chart, not a surviving L2
gauge zero mode. For an allowed boundary-vanishing gauge eta,
f_+ changes by -(gamma-partial_y)eta_+. Integration by parts gives

    Theta_gamma(f_new)=Theta_gamma(f)-eta_+(0).          (8)

Consequently Theta is invariant under such small gauges. The necessary
endpoint decay holds for H1 eta and a>0. No physical charge or symplectic
norm is being assigned to this invariant.

## 3. Retarded data and a compact counterexample

For zero initial gauge vector the outgoing characteristic is

    zeta_+(t,y)=integral_0^t f_+(s,y+t-s) ds,
    c(t)=zeta_+(t,0)=integral_0^t f_+(s,t-s) ds.         (9)

The incoming characteristic has initial value zero and boundary value
c(t) for endpoint preservation. For y<t it is
c(t-y)+integral_(t-y)^t f_-(s,y-t+s)ds; for y>=t its integral starts
at zero. These formulas also prove the retarded estimate
||zeta(t)||<=integral_0^t ||f(s)||ds. Fixing the boundary clock requires
c(t)=0 for every t, equivalently (3) applied to the Laplace transformed
forcing whenever that transform exists. If an initial gauge vector is
prescribed and nonzero, retain it: the frequency source is
f_hat+zeta(0), and the condition is Theta_gamma(f_hat+zeta(0))=0.
It cannot be erased by silently choosing different initial coordinates.

Choose nonnegative, nonzero smooth bumps q(t) supported in (0,1) and
k(y) in (1,2). Let A=q(t)k(y), B=C=p_0=0. The perturbation is zero near
the initial slice and both spatial endpoints. Yet

    c(t)=integral_0^t q(s)k(t-s) ds

is positive for some t in (1,3). Hence no zero-initial, endpoint- and
clock-preserving solution of (2) exists for this compact probe. Its
amplitude can be arbitrarily small. This counterexample concerns the
kinematic space of metric probes. It does not establish that the full
retarded equations admit this probe as an unforced physical solution.

## 4. What the clock change does to the previous reconstruction space

For a frequency source compact in ell<y<b, let c=Theta_gamma(f).
Below ell, (4)-(5) give exactly

    xi^t=c cosh(gamma y), xi^y=c sinh(gamma y).          (10)

If the original h vanishes there, its optical image from (1) is

    p=(F/r)c sinh(gamma y),
    v=(F_r/2-F/r)c sinh(gamma y)-gamma c cosh(gamma y). (11)

The actual charged asymptotics r=1/(lambda y)+O(y),
F=1/(lambda y^2)+O(1) imply

    p=gamma c+O(y^2), v=-gamma c+O(y^2),
    p+v=O(y^2), r p=gamma c/(lambda y)+O(y).            (12)

This is a pure coordinate variation in the source-free outer region.
Its gauge-invariant curvature Psi=delta R-R_r delta r is zero exactly,
as can also be verified in the complete optical curvature identity.
Nevertheless r p is NOT in L2(dy) if c is nonzero. The regular L2
auxiliary reconstruction in the preceding work therefore cannot simply
be declared a complete representation of these transformed probes.
Allowing the clock to transform is not accomplished by dropping its
nonintegrable boundary contribution. If c=0 and the original source is
compact, both gauge components vanish below ell and this issue is absent.

Under the passive convention h_new=h-Lie_xi g, the original boundary
clock is represented in the new coordinates by t-c(t) to first order.
The boundary metric/frame and state data must be transported along with
it. Equation (12) displays the associated metric change; one cannot
hold the original boundary sources fixed while silently discarding c.
No boundary action, thermal-state response or physical admissibility of
this enlarged class has been calculated here.

## 5. Horizon regularity and its exceptional frequency

For the actual simple horizon let kappa_h=F_r(r_h)/2>0. In ingoing
coordinates v_h=t-y,r, the metric is -F dv_h^2+2dv_h dr+r^2dOmega^2.
Outside a compact forcing region the outgoing homogeneous gauge is
zeta_+=a_+ exp(gamma y), zeta_-=0, multiplied by exp(gamma t). Its
ingoing components are

    xi^(v_h)=0,
    xi^r=-(a_+/2)F exp(2gamma y) exp(gamma v_h).         (13)

Since y=-(log(r-r_h))/(2kappa_h)+smooth(r), its radial coefficient
scales as (r-r_h)^(1-gamma/kappa_h). For real gamma>kappa_h, continuity
at the future horizon excludes a_+. In this regime the L2 choice in
(4) agrees with the necessary horizon regularity condition for compact
forcing. The ingoing homogeneous solution exp(-gamma y) has finite
xi^(v_h) and xi^r proportional to F in these coordinates and is smooth.

At gamma=kappa_h the outgoing radial coefficient can instead be smooth
and nonzero. This horizon displacement is excluded by the L2 optical
gauge norm but not by smoothness alone. Therefore (7) MUST NOT be
asserted for every positive frequency under horizon regularity alone.
For example F=2kappa_h(r-r_h) gives exactly a constant radial coefficient
at that frequency. No exceptional physical metric mode is inferred.
Zero-initial retarded characteristic statement (9) is independent of
this separated-mode exception and retains its specified initial data.

## 6. Research decision and evidence scope

Before using a global inverse of the two optical rows to conclude
physical charged stability, do one of the following in the ORIGINAL
physical initial/boundary problem: prove its admitted perturbations
satisfy Theta=0, or retain the transported boundary clock/frame and
its metric/state boundary response, or use a metric gauge that keeps
those data explicit. The condition here is distinct from the fast
curvature forcing moment and from the Ward constraint equations.
Neither of those earlier identities proves this gauge condition.
An optical perturbation already in the chart has f_+=0; our result
does not invalidate the restricted annular shape theorem.

All formulas here are kinematic. On an off-shell seed, a coordinated
metric/state diffeomorphism gives delta calE=Lie_xi calE, not zero.
The background residual and state/boundary contacts remain necessary
in the full response. Thus Psi=0 in (11) is not a new null vector of
the complete action Hessian. No full physical perturbation domain,
global shape inverse, coupled moment value or charged mode is supplied.

[Harlow and Wu, introduction and section 4.3](https://arxiv.org/html/1906.08616v4)
explain why spatial boundary conditions and asymptotic domains are part
of the Hamiltonian problem. That established framework motivates keeping
the boundary data explicit. It does not supply (2)-(13), a gauge-range
theorem for this exterior, or the retained fermion/state boundary response.
The transport proof and metric consequences above are derived here.

Freeze this note, the checker and source review before acceptance.
Verify the Lie derivative directly in four dimensions, the characteristic
signs, exact range witnesses, gauge invariance of Theta, all boundary
jets, the pure-gauge curvature identity and the horizon exception.
Compare independently integrated gauge ODEs with Green quadratures for
real and complex frequencies, including compatible and incompatible
compact sources. Retain the complete incoming tail. Check actual seed
outer jets against radial quadrature at successive radii. The finite
fixtures are not physical modes, a Dirac-kernel computation or formal
verification. Preserve failed snapshots and outcomes; no Lean selected.

Remaining to-do list: physical gauge coverage or explicit boundary clock
response, global shape/cross-block control, total compatibility and charged
stability; S167 signed source, higher response, nonlinear existence and
compatibility; exterior/state, charge/return, opening, all-sector stability,
transport, preparation, separate creation and final audit.
