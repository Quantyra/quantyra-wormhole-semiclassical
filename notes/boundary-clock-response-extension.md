# Removing the coordinate singularity from the boundary clock response

2026-10-21, E003/S168. The clock-associated curvature-free boundary
column is the sum of a transported coordinate column and a metric
remainder whose outer conformal amplitude belongs to the old L2 space.
The remainder's response against compact bulk readouts is the unique
limit of the existing interior response. This resolves the outer input
extension for that column. It does not evaluate the response numerically,
provide a boundary force law, select physical boundary data, or settle
the clock compatibility and charged stability problem.

Use the unchanged seed, action, prepared state, magnetic tower and common
subtraction. Keep the [ADM vertices](adm-dirac-boundary-response.md),
[transported-state identity](transported-clock-dirac-response.md), and
[boundary reconstruction](boundary-clock-curvature-reconstruction.md).
No boundary susceptibility is set to zero. In particular the whole
curvature-free column is not a coordinate variation.

## 1. A causal coordinate lift with the required boundary jets

Write n=log(F)/2, l=log(r), r_y=-F. Under a positive coordinate pullback
xi=a partial_t+b partial_y, its ADM logarithmic metric amplitudes are

    A_xi=n_y b+a_t, B_xi=n_y b+b_y, p_xi=l_y b,
    d_xi=A_xi-B_xi=a_t-b_y,
    w_xi=A_xi-p_xi=a_t+(n_y-l_y)b,
    s_xi=b_t-a_y.                                             (1)

Here w denotes the optical mass/shape amplitude, and s the shift.
Let rho be one near y=0 and zero beyond a fixed outer collar. For a
smooth clock displacement c(t), zero near the initial surface, choose

    a=rho[c+c'' y^2/2],
    b=rho[c' y+c''' y^3/6].                                  (2)

This uses only time derivatives at the current time. It preserves the
endpoint, has boundary clock a(t,0)=c(t), and vanishes toward the horizon.
Exactly,

    d_xi=-rho'[c' y+c''' y^3/6],
    s_xi=rho c'''' y^3/6-rho'[c+c'' y^2/2].                    (3)

In particular the boundary shift is zero. The geometric projector
there remains Pi_0 in this frame; the clock, frame, cutoff and state
are nevertheless transported as in the preceding theorem. Near zero,

    l_y=-1/y-2lambda y/3+3m lambda^2 y^2/2+O(y^3),
    n_y-l_y=lambda y-3m lambda^2 y^2+O(y^3).                  (4)

Thus p_xi=-c'+O(y^2), w_xi=c'+O(y^2), d_xi=0 and s_xi=O(y^3).
For a clock-associated optical boundary column p->-v_0, V->v_0,
the transported clock has c'=v_0. This uses the actual clock source;
an independently driven metric with a different clock is outside this
identification. In the earlier passive gauge convention the clock is
t-c_old, and v_0=-c_old', so the signs agree with c=-c_old here.

Subtracting (1)-(2) from this column leaves

    p_rem=O(y^2), w_rem=O(y^2), d_rem=0, s_rem=O(y^3),
    q_rem=r p_rem=O(y) in L2 near zero.                       (5)

The endpoint and clock of the remainder are unchanged. Initial clock
data cannot be discarded: in Laplace variables v_0=gamma c_hat-c(0).
Equation (2) is used with zero initial clock jets here. Independent
initial-state terms remain the ones in the transported-state theorem.

## 2. Exact optical lift for the fixed-frequency proof

For gamma with Re gamma>0, an alternative local lift is

    a=(v_0/gamma) rho cosh(gamma y),
    b=(v_0/gamma) rho sinh(gamma y).                           (6)

It has the same boundary clock. Where rho=1, d_xi=s_xi=0 exactly,
and its optical amplitudes are

    p_xi/v_0=(l_y/gamma)sinh(gamma y),
    V_xi/v_0=cosh(gamma y)+(n_y-l_y)sinh(gamma y)/gamma.        (7)

This is a fixed-frequency proof device, not a causal reconstruction
algorithm. The cosh lift corresponds to advanced as well as delayed
clock values in time; the local polynomial lift (2) avoids that issue.
Their difference vanishes as a=O(y^4), b=O(y^5) at the boundary and
has zero clock trace. Its full response is its known covariant coordinate
response, including off-shell contacts. No uniform large-frequency
bound or causal inverse is inferred from the cosh factors.

For the unit curvature-free column V=chi, with chi=1 near zero,

    q_B=-r chi-R_gamma u_chi,
    u_chi=4r F^-1(chi_yy-gamma^2 chi)-12chi_y
                                      +r(C_v+H_R)chi.         (8)

The existing real-frequency reconstruction extends to each fixed
Re gamma>0 without assuming a new weighted resolvent bound. Choose a
real Gamma above the established threshold. On L2 sources define

    R_gamma=[I+(Gamma^2-gamma^2)(K+gamma^2)^-1] R_Gamma.       (9)

The positive self-adjoint K makes the ordinary resolvent bounded and
holomorphic there. The resolvent identity gives the required equation
with singular multiplier F/6 in distributions. It agrees with the
previous real construction and is independent of Gamma. This is an
old-L2 statement, not a complex-frequency uniform bound on the larger Z
space. Since u_chi is in L2, (8) is defined by (9).

Put z=1/r. The analytic outer IVP is

    z_y=lambda+z^2-2m z^3+Qhat z^4+c3 z^5+c4 z^6, z(0)=0.

Consequently K=-partial_y^2+6/y^2+2lambda+O(y). In the region chi=1,

    F u_chi/6=2(lambda-gamma^2)/(3lambda y)
                                 +(2gamma^2+4lambda)y/9+O(y^2).

The L2 particular solution q_c=R_gamma u_chi has

    q_c=a_1 y+a_log y^3 log(y)+O(y^3),
    a_1=(lambda-gamma^2)/(9lambda),
    a_log=-(gamma^2+lambda)(gamma^2+2lambda)/(45lambda).       (10)

These coefficients follow from
(-partial_y^2+6/y^2)y=6/y and
(-partial_y^2+6/y^2)y^3 log(y)=-5y. The regular homogeneous
coefficient at y^3 is global data and is not set to zero. The singular
y^-2 homogeneous branch is excluded by L2. Recurrence, or variation of
parameters with the regular y^3 solution, controls the differentiated
remainders; in particular the resulting p_B is C3 at the endpoint.

For v_0=1, subtraction of (7) gives the explicit leading coefficients

    p_B=-1+(gamma^2-lambda)y^2/9+O(y^4 log(y)),
    p_rem=5(gamma^2+2lambda)y^2/18
                              -3m lambda^2 y^3/2+O(y^4 log(y)),
    V_rem=-(gamma^2+2lambda)y^2/2
                                  +3m lambda^2 y^3+O(y^4).    (11)

The ratio of the quadratic coefficients is -5/9, agreeing with the
leading curvature-free constraint 36 p_2+20 V_2=0. The invariant
curvature remains zero exactly, since both the original column and
the coordinate column have zero delta R-R_r delta r. The remainder
is generally nonzero; removing it would change the perturbation.

## 3. The remainder belongs to the previously controlled core closure

First fix gamma in a compact subset of Re gamma>0 and a small collar.
Use (6) so the remainder is optical near zero; its lapse/shift mismatch
outside that neighborhood is supported strictly in the interior and
is handled by the existing full ADM vertices. No horizon response is
approximated by an outer model. The conformal reconstruction tail stays
as in (8); global energy membership at the horizon is a separate issue.

Let theta_delta(y)=theta(y/delta), with theta=1 on [0,1] and zero
on [2,infinity). Remove this endpoint piece from p_rem,V_rem to make
interior probes. If f(y)=f_2 y^2+f_3 y^3+O(y^4 log(y)), with the
differentiated remainder from (10), its even tail has the scaling

    f_delta(y)=theta_delta f=delta^2 psi_delta(y/delta),
    sup_(0<delta<delta_0)||psi_delta||_H3(R)<infinity.

The |y|^3 term has bounded third weak derivative, so it causes no
failure of this H3 bound. In the physical-half logarithmic form norm

    ||f||_V^2=(1/2) integral (1+k^2)^2[1+log(2+k^2)]
                                               |f_even_hat(k)|^2 dk,

Fourier rescaling and log(2+k^2/delta^2)
<=2|log delta|+log(2+k^2) prove

    ||f_delta||_V^2 <= C delta(1+|log delta|).                (12)

The ordinary weighted local terms satisfy, by direct rescaling,

    integral_0^(2delta) [|p_delta''|^2+|p_delta'/y|^2
                                      +|p_delta/y^2|^2]dy <=C delta,
    integral_0^(2delta) |log y| |V_delta''|^2 dy
                                      <=C delta(1+|log delta|). (13)

These weights dominate the outer coefficients of the retained local
optical Hessian at fixed frequency. Lower terms obey the same or better
bounds. They include the fixed-density and off-shell contacts; none is
removed by imposing background field equations.

The existing uniform canonical comparison controls all interior outer
probes by C_gamma ||V||_V^2 at each nonzero Matsubara frequency. To see
the upper bound as well as the lower one, use its absolute comparison
to Kflat(sqrt(L)); Kflat has at most fourth-logarithmic growth, the
half-density coefficients have bounded weak derivatives through order
three, and its lower-order contact errors are bounded in H2. The same
four-positive-anchor spectral interpolation used in the annular theorem
then gives locally uniform bounded canonical forms at complex causal
frequencies on this core closure. No constant depends on the distance
of the support from zero. The constants remain uncalibrated.

Thus (12) makes the optical endpoint tail vanish in the actual canonical
response form. Polarization gives convergence against every member of
that optical form space. For any fixed smooth compact bulk readout of
the general metric, its overlap with the endpoint tail eventually
vanishes. The remaining off-diagonal retarded kernel is smooth in the
separated canonical spatial neighborhoods at fixed causal frequency;
the reflected diagonal is also separated. Its pairing with the tail
tends to zero. The full angular sum and common subtraction are those
of the retained response, not a sum of separately renormalized modes.

Equations (12)-(13) therefore define a unique response of the regular
remainder against compact bulk readouts, independent of theta. This is
an input-extension theorem, not a uniform bound on all global readouts,
a standalone boundary-force susceptibility, or a solved coupled inverse.
In particular the original V=1 boundary value cannot simply be cut off:
its endpoint tail has fourth-log norm of order delta^-3 |log delta|,
instead of (12). The coordinated subtraction is essential.

## 4. What is now determined, and what is still open

Let R denote the complete linearized equation tensor response, with
the state, frame and domain sources included. For the clock-associated
column H_B and a fixed smooth compact bulk readout J,

    <J,R(H_B)>=<J,Lie_xi calE>
                         +lim_(delta->0)<J,R(H_rem,delta)>.   (14)

The remainder probes use the unchanged prepared state and boundary
sources. The first term includes the off-shell seed residual and the
observable/density contacts already derived. It is not zero unless the
appropriate background equation vanishes. Changing rho or choosing
(2) in place of (6) transfers a boundary-clock-free coordinate column
between the two terms; covariance preserves their sum. No new state
parameter or thermal reset is introduced.

Equation (14) gives the bulk response as a controlled limit of the
existing interior operator, rather than a new arbitrary boundary
susceptibility. It does not supply numerical values of that operator,
the boundary equation that fixes the physical clock, a boundary action
Hessian evaluated on nonvanishing traces, or the global physical
metric domain. Those obligations remain before a charged mode or a
stability verdict. The zero-curvature boundary column has not vanished.

The local covariance framework is established in
[Sanders, section 4.3/Theorem 4.18](https://arxiv.org/pdf/0911.1304v2)
and [Zahn, sections 4.1-4.2](https://arxiv.org/pdf/1210.4031v3).
They do not supply the timelike paired boundary theorem. Its current
argument is pinned in the transported-clock note; the subtraction,
outer coefficients and closure estimates are derived above. The
canonical continuity input is the previously pinned uniform comparison
and retarded interpolation, with their original limitations intact.

The checker verifies exact coordinate jets and Frobenius coefficients,
charged reconstruction against the boundary expansion, and independent
Fourier quadrature of the cutoff convergence. These are supplementary
identity and numerical fixtures, not an evaluated charged Dirac kernel.
An exploratory unsimplified symbolic series calculation was interrupted
before producing output; the retained checker uses expanded Laurent
coefficients and does not depend on that unfinished computation.

Remaining to-do list: evaluate the regular bulk response and the physical
boundary/clock equation; global coupled charged stability; nonlinear
existence and compatibility; exterior/state, charge/return, opening,
all-sector stability, transport, preparation, separate creation and audit.
