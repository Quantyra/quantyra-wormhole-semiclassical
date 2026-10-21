# Conserved sources and the zero-fast compatibility conditions

2026-10-21, E003/S168. Ordinary argument before focused verification.
The stable semiclassical wormhole goal remains incomplete.

Conservation, smooth thermal-horizon behavior and zero outer energy flux
do not by themselves enforce the source conditions of the auxiliary
zero-fast curvature domain. We prove this on the retained smooth charged
seed at each admitted thermal harmonic n>=2. We also pull those conditions
back to the original stress-source variables, including the shape source
induced by curvature reconstruction. Arbitrary conserved stresses are not
thereby realized as variations of the selected quantum state.

## 1. Original source work and its exact pullback

Retain the action, paired fermion domain, HH state, clock, magnetic charge,
subtraction and finite contacts of the
[coupled boundary index theorem](coupled-thermal-boundary-index.md).
Fix gamma=n kappa_h above both of its finite auxiliary thresholds, with
integer n>=2. The inward optical coordinate is y=-x, H=L2(dy), and the
shape energy space E_gamma embeds continuously and injectively into H.
Write that inclusion as j when it helps distinguish the two norms. Use
inner products conjugate linear in the first argument. Put

    epsilon=G/(60pi), alpha=epsilon/(2G),
    q=r p=R_gamma U+Q_gamma w, w=r v.

Define the independent normalized source pair (f,g) in H x H by its
original linear work, with time-mode normalization divided out,

    F_source(p,v)=alpha <r p,f>_H+<r v,g>_H.             (1)

The sign convention is fixed by equating the congruent Hessian to this
work. Its strong trace source is S=-epsilon f. Substitution of the full
reconstruction, without dropping its Q_gamma part, gives

    F_source(U,w)=alpha <R_gamma U,f>+ell_(f,g)(w),
    ell_(f,g)(w)=alpha <Q_gamma w,f>+<j w,g>.             (2)

Thus even a source with g=0 usually induces a shape functional. Let
a_gamma(w,w_ell)=ell_(f,g)(w). All retained off-shell density and radial
contacts remain in the operators below; adding an independent stress at
fixed metric does not remove them. The existing Schur elimination gives

    w=w_ell-K_cross^dagger U,
    L_gamma U=f+epsilon^-1 D_gamma w_ell,               (3)
    L_gamma:D_min -> H, ind L_gamma=-1.

D_min removes both leading fast traces A_c=A_s=0. It is a candidate
domain, not an adopted physical boundary prescription. The range of L
is closed. Its full adjoint nullspace has dimension dim ker L+1>=1;
neither its actual dimension nor its profiles have been evaluated.

For each z in ker L_gamma^*, define nu_z in E_gamma by

    a_gamma(w,nu_z)=epsilon^-1 <D_gamma w,z>_H,
    mu_z=z+alpha Q_gamma nu_z.                          (4)

Boundedness of D:E_gamma->H and the Riesz theorem justify (4). Then

    <z,f+epsilon^-1 D_gamma w_ell>
          =<mu_z,f>_H+<j nu_z,g>_H.                    (5)

Indeed the second term on the left equals a_gamma(nu_z,w_ell), hence
ell_(f,g)(nu_z). This also fixes the conjugations in (5). Solvability of
the complete auxiliary pair is exactly the vanishing of (5) for EVERY
adjoint null vector. This says nothing about uniqueness when ker L is
nonzero. The original-source compatibility map is finite dimensional,
but its numerical coefficients are not supplied here.

For nonzero z, the pair (mu_z,j nu_z) cannot vanish: j nu_z=0 implies
nu_z=0, hence mu_z=z. Smooth compact pairs are dense in H x H, so at
least one such pair has a nonzero overlap (5). Scaling it yields
arbitrarily small incompatible sources. This density argument uses the
full original pair and retains the induced shape source; it does not
replace physical sources by independent arbitrary Schur-row forces.

## 2. Every compact pair has a conserved stress completion

On the actual static seed, write

    ds^2=N(x)^2(-dt^2+dx^2)+r(x)^2 dOmega^2,
    n_0=log N, l=log r, W=N^2 r^2, -infinity<x<0.

The horizon is smooth and simple with surface gravity kappa_h, and the
outer end has the retained AdS asymptotics. This W is a positive weight,
not a Ward forcing or a quantum second vertex. For an independent
spherical mixed stress T, set

    tau=T^t_t+T^x_x+2T^theta_theta, o=T^t_t+T^x_x,
    k=T^x_t=-T^t_x, E=W T^x_x, J=W k.

All radial coefficients multiply exp(gamma t). The conformal metric test
has h^t_t=h^x_x=2(p+v), h^theta_theta=h^phi_phi=2p. Direct contraction
and angular integration give source work

    (1/2) integral sqrt(-g) T^{ab} h_ab
             =4pi integral dy W(tau p+o v).

For complex mode tests take conjugates of p,v in this pairing. Matching
(1) fixes, with its prescribed sign convention,

    tau=epsilon r f/(8pi G W), o=r g/(4pi W).           (6)

Any common paired-copy or time normalization multiplies both sides.
Smooth compact f,g give smooth compact tau,o, because the prefactors are
smooth and positive in the open exterior.

The exact four-dimensional Ward reduction from
[spherical constraint transport](spherical-retarded-constraint-recovery.md)
applied to an independent conserved stress at fixed metric is

    gamma E-J'=F, F=gamma W o,
    gamma J-E'=G_1, G_1=-W[(n_0'-l')o+l' tau].         (7)

Both F and G_1 have compact radial support. This uses no equation of
motion for the seed metric. For C_+=E+J, C_-=E-J, the unique L2 solution
with J(0)=0 and no exponentially growing horizon coefficient is

    C_-(x)=integral_-infinity^x exp[-gamma(x-s)]
                                      (F-G_1)(s) ds,
    C_+(x)=exp(gamma x) C_-(0)
               +integral_x^0 exp[gamma(x-s)] (F+G_1)(s) ds. (8)

Differentiate (8) to recover (7). The previous energy estimate gives
||(E,J)||_L2 <= ||(F,G_1)||_L2/gamma. Reconstruct the stress by

    T^x_x=E/W, T^t_t=o-E/W,
    T^theta_theta=T^phi_phi=(tau-o)/2, k=J/W.           (9)

Equations (7) imply all four covariant divergence components vanish;
the angular ones vanish by spherical symmetry. This completion is
linear and has no need to tune an overlap of tau,o. It is neutral:
no independent electric current or change of magnetic flux is added.
This does not determine the actual Dirac current response.

Below the compact source support, C_-=0 and C_+=A exp(gamma x).
The covariant stress there is exactly

    T=A exp[gamma(t+x)] (dt+dx)^2/(2r^2).

In regular ingoing v=t+x and Kruskal V=exp(kappa_h v), it becomes

    T=A V^(n-2) dV^2/(2 kappa_h^2 r^2).               (10)

The seed r is smooth in the Kruskal product. For integer n>=2, (10)
extends smoothly at the bifurcation surface and at the Euclidean bolt:
under thermal continuation V is a smooth complex Cartesian coordinate,
and V^(n-2)dV^2 is polynomial. Its real and imaginary Euclidean parts
are smooth. Only the time harmonic is continued; the compact radial
sources need not be analytic. The n=1 expression would generally have
a V^-1 pole. That harmonic is deliberately excluded, consistently with
the [first-harmonic Dirac flux result](first-thermal-horizon-flux.md).
There are arbitrarily high admitted integer harmonics above the finite
auxiliary thresholds.

Near the outer end, E,J are bounded and J(0)=0. Since W=O(r^4) with
positive leading coefficient, mixed orbit stresses are O(r^-4);
their covariant orbit components are O(r^-2). The conserved energy
flux through a sphere is proportional to J and vanishes at the outer
boundary. This is a stress-flux condition, not a selection of the metric's
fast phase. The construction can be copied to the paired spherical
sector; it does not prove the fermion domain realizes the copied stress.

Combining (5), density, and (6)-(10) proves the promised result:
conservation plus these regularity and stress-flux conditions alone
cannot make every compact original source pair zero-fast compatible.
It supplies arbitrarily small conserved incompatible sources, without
constructing a matter state or an exact semiclassical metric.

## 3. What the actual quantum preparation still has to supply

The fixed massless Dirac theory has a local state-independent trace
anomaly. At fixed metric, gauge field, renormalization and boundary
domain, the bulk trace of an independent admissible state variation is
therefore zero. In (6) that means f_state=0, not an arbitrary compact f.
The surviving compatibility conditions are precisely

    <j nu_z,g_state>=0 for every z in ker L_gamma^*.     (11)

If nu_z is nonzero, some compact traceless conserved source violates its
condition, by the same density argument. If nu_z=0, every traceless
source passes that particular condition. Neither alternative has been
established for the full actual adjoint profiles; no nonzero nu_z is
assumed. Even the former alternative would not prove that an admissible
Dirac state realizes the incompatible source.

For a strictly fixed HH preparation without an independent driver or
initial-state change, there is no freely adjustable g_state. The
metric-driven quantum response is already part of L_gamma and its shape
block; it must not be counted again as an independent source. Initial
metric/state data or an authorized preparation driver must first be
inserted in the original rows, including their initial reconstruction
and off-shell contacts, before (5) is applied. A homogeneous physical
mode problem likewise cannot be decided by an arbitrary-source example.

The distinction between metric-driven response and independent state
variation, and the local state-independent anomaly, agrees with
[Anderson, Molina-Paris and Mottola, sections I and III](https://arxiv.org/html/0907.0823v2).
Their de Sitter equations do not compute the charged preparation map.
No short-scale exclusion prescription from that paper is adopted.

The [conditional nonzero-fast family obstruction](combined-shape-cutoff-obstruction.md)
remains in force. The present result prevents conservation alone from
closing the zero-fast alternative. The physical next step is to identify
the admitted preparation/initial-data map into the original rows and
evaluate its adjoint overlaps, while retaining the separate static signed
condition. No new boundary condition, counterterm or order reduction is
selected. No continuous causal inverse, physical instability, nonlinear
solution or stable wormhole is established.

Focused checks derive the divergence directly from four-dimensional
Christoffel symbols, check source work and complex weighted adjoints,
compare the Green construction with an independent ODE integrator, and
check the thermal horizon powers. Numerical transport fixtures use an
explicit smooth model geometry, not an evaluated charged Dirac kernel.
The proof on the charged seed is the geometry-independent argument above.
No Lean campaign is selected. Freeze sources before acceptance and retain
every run and any failed source snapshot.

Remaining to-do list: actual preparation/initial-data image and coupled
adjoint overlaps; physical gravitational boundary domain and static signed
condition; continuous causal constraints, nonlinear existence and
all-sector stability; exterior/state, charge/return, opening, transport,
preparation, separate creation and full audit.
