# Causal horizon localization and the reflected zero-mode response

2026-10-21, E003/S168. **Finite-time compact exterior readouts have a
cutoff-independent causal response, even when a global static-slice
perturbation chart fails. The magnetic zero modes give an explicitly
evaluated reflected stress response.** The massive angular response,
physical gravitational domain and coupled stability remain unresolved.

Keep the [ingoing areal formulation](ingoing-areal-retarded-system.md),
the [full ADM vertices](adm-dirac-boundary-response.md), prepared HH state,
zero-source reflecting pair and common four-dimensional subtraction.
The zero-mode calculation uses the already declared
[Polyakov allocation](stationary-lll-quantum-exterior.md); it is one
component of the complete action, not a replacement for that action.

## 1. Exact causal localization of the full response

On a static exterior y>=0, radial null curves satisfy |dy/dt|=1. Let a
readout have compact support in 0<=t<=T and 0<y<=Y, and let a smooth metric
source vanish before t=0, with the initial state fixed and the physical
outer source/domain unchanged. Boundary corner compatibility is assumed.
Multiply the source by a smooth radial cutoff chi_L equal to one up to L
and zero beyond L+1, with L>Y+T. Its difference from the uncut source is
outside the causal past of the readout. Consequently

    delta <O>_ret[h]=delta <O>_ret[chi_L h].                (1)

This is equality, not an estimate in powers of 1/L. To establish it, the
Dirac equation's characteristic principal symbol gives speed at most one;
the angular mass is a zeroth-order term. The paired local reflecting
projector has zero normal current, so the same causal energy argument
holds at the boundary. Reflected paths cannot shorten the radial travel
time. Fermionic causal propagation gives spacelike commutativity of the
even stress operators, hence support of the retarded commutator in the
causal cone. Every local second vertex and renormalization contact has
diagonal support and agrees near the readout.

The statement holds for each angular channel and for the complete
renormalized distribution: the common subtraction is local, and a
distributional limit retains the closed causal support. It is not an
assertion that a finite spectral cutoff preserves exact causality. At
zero inter-exterior coupling, purification correlations do not transmit
an even-operator retarded signal between the two copies. The prepared
state is not rethermalized after applying the cutoff.

For each finite T and cutoff, a sufficiently small real source amplitude
defines a spacelike static-slice family on the relevant compact region.
No common amplitude neighborhood for all T is required to define this
linear response. This supplies a consistent family of local retarded
responses without a single global static-slice perturbation chart.
It does not give a uniform infinite-time estimate, a stress trace directly
on the horizon, or permission to exchange a Laplace transform with an
uncontrolled infinite-time limit.

## 2. What high-frequency ingoing regularity does give

For a simple future horizon let F_r(r_h)=2kappa_h>0 and
r-r_h~a exp(-2kappa_h y). A smooth ingoing mode has

    delta F=e^(gamma v) fbar(r), delta psi=e^(gamma v) zetabar(r),
    v=t-y, k_gamma=e^(-gamma y)fbar/F,
    zeta_gamma=e^(-gamma y)zetabar, d_gamma=zeta_gamma+k_gamma.

Thus for Re gamma>2kappa_h, all fixed finite y-derivative orders allowed
by the smooth horizon data obey

    |partial_y^j k_gamma|+|partial_y^j zeta_gamma|
         <=C_(j,gamma) exp[-(Re gamma-2kappa_h)y]          (2)

on a sufficiently long horizon tail. For example
k_gamma'=(F_r-gamma)k_gamma-exp(-gamma y) fbar_r.
Generic fbar(r_h) nonzero makes k_gamma unbounded below the threshold;
at equality its modulus need not vanish. This is a domain threshold,
not a pole or a condition excluding slower physical disturbances. If
fbar=O(F), the corresponding decay threshold is Re gamma>0.

There is a useful uniform angular graph estimate. For
h_nu=-i alpha partial_y+beta m_nu, m_nu=lambda_nu sqrt(F)/r,
put K=sup_y|F/r-F_r/2|, finite on the retained simple-horizon AdS seed.
Then |m_nu'|<=K|m_nu|. In the folded convention alpha=tau_3 sigma_2,
beta=I sigma_1 and chi=tau_3 sigma_1. The boundary form involving
i alpha beta vanishes on Pi_0=(I+chi)/2 because it anticommutes with chi.
Integration by parts and Young's inequality give

    ||phi'||^2+(1/2)||m_nu phi||^2
                   <=||h_nu phi||^2+(K^2/2)||phi||^2.     (3)

The constant is independent of angular number, mass sign and multiplicity.
It holds first on the regular operator core and then on its graph closure.
The first ingoing vertex is

    B_gamma=i(I-alpha)(d_gamma partial_y+d_gamma'/2)
                                        +beta m_nu w_gamma,
    w_gamma=zeta_gamma+k_gamma/2.

Equation (3) bounds B from the h graph domain to L2 by the bounded
coefficients d,d',w. Cutting off the horizon tail changes that graph map
by at most C_gamma exp[-(Re gamma-2kappa_h)L]. This also holds on the
direct sum of all angular channels. Fixed admissible outer traces are
still required; no outer gravitational domain is selected here.

This graph statement does not bound a renormalized stress trace, perform
the angular sum, prove a retarded gravitational inverse, or supply an
infinite-time estimate for (1). In particular it does not by itself
justify a Laplace transform of the full quantum response. It removes a
specific vertex-domain problem on the high-growth-frequency half-plane.

## 3. Exact canonical zero-mode susceptibility

There are q massless two-dimensional Dirac fields per exterior, with
central charge c=q and prepared beta_T=2pi/kappa_h. A constant reflecting
phase does not affect their stress. Unfold the outgoing component at y>0
to s=y>0 and the incoming component to s=-y<0. Reflection joins them into
one chiral field on the line, with h=i partial_s and

    t_0=<T(s)>=c pi/(12 beta_T^2)=c kappa_h^2/(48pi).

There is no extra factor of two from unfolding or the other exterior.
For an areal perturbation compactly supported away from y=0 and infinity,
the zero-mode Hamiltonian follows exactly from the ADM formula:

    h_in=-i partial_y,
    h_out=i[(2 speed-1)partial_y+speed'],
    speed=1/[2-exp(psi)(1+k)].

At first order H_int=integral g(s)T(s)ds, where g(s)=2d(s) for s>0
and zero for s<0. Here d=zeta+k. Extend these compact profiles smoothly
by zero; nonzero endpoint jets would require their boundary contacts.

For an independent derivation of the central and thermal coefficients,
let f_beta(E)=(1+exp(beta_T E))^-1. Momentum eigenvalues are -p and
the source matrix element is -(p+p')g(p-p')/2. Put ell=p-p'. The
linear covariance equation gives

    delta t_hat(ell,gamma)=-ic g_hat/(gamma-i ell)
       * integral dE/(2pi) E^2
           [f_beta(E+ell/2)-f_beta(E-ell/2)].             (4)

The integral is finite without discarding the filled sea. Since -f_beta'
has total mass one, zero mean and variance pi^2/(3 beta_T^2), integration
over shifts in [-ell/2,ell/2] gives exactly

    integral dE E^2[f_beta(E+ell/2)-f_beta(E-ell/2)]
                  =-ell pi^2/(3 beta_T^2)-ell^3/12.     (5)

For completeness the variance follows by setting u=exp(beta_T E):
integral_0^infinity u^a/(1+u)^2 du=Gamma(1+a)Gamma(1-a)
=pi a/sin(pi a). Its second derivative at zero is pi^2/3.

For negative ell the same odd formula holds. Equations (4)-(5) imply
the retarded stress equation, with zero initial state variation,

    (partial_t-partial_s)delta t=P[d],
    P[d]=4t_0 partial_s d-c partial_s^3 d/(12pi)
         =c(kappa_h^2 d'-d''')/(12pi).                  (6)

For nonzero independent initial data add delta t_initial(s+t).
There is no instantaneous change in the Fermi occupations. Equation (6)
is the thermal chiral stress response in the original zero-mode
allocation, not a new finite-counterterm choice.

The exact time-domain solution is

    delta t(t,s)=integral_0^t P[d](tau,s+t-tau)dtau.      (7)

For Re gamma>0 its outgoing resolvent is

    delta t_hat(s)=integral_s^infinity
                          exp[-gamma(u-s)]P[d_hat](u)du. (8)

Compact radial input makes this entire right-half-plane expression
unambiguous in the massless sector, independently of the general graph
threshold in section 2. The opposite characteristic direction from the
coordinate transport in the preceding note is essential: this is the
OUTGOING quantum stress, which reflects into the incoming component.

## 4. Evaluated reflection and physical off-source flux

For d_hat supported in [a,b] with 0<a<b, define

    D_gamma=integral_a^b exp(-gamma u)d_hat(u)du,
    R_gamma=c gamma(kappa_h^2-gamma^2)D_gamma/(12pi).     (9)

Three integrations by parts in (8) give

    delta t_hat(s)=exp(gamma s)R_gamma, s<0,
    delta t_hat(s)=c/(12pi)[d_hat''+gamma d_hat'
       +(gamma^2-kappa_h^2)d_hat
       +gamma(kappa_h^2-gamma^2)integral_s^infinity
                           exp[-gamma(u-s)]d_hat(u)du], s>0. (10)

This is an evaluated nonlocal response. At gamma=kappa_h the reflected
coefficient vanishes, but the bulk response is generally
c(d_hat''+kappa_h d_hat')/(12pi), not zero. This cancellation is neither
a pole nor a statement about stability of the coupled gravitational system.

In the unchanged metric region beyond the source, the outgoing variation
is zero and the incoming one is exp(-gamma y)R_gamma. Therefore in regular
ingoing coordinates its nonlocal zero-mode stress is

    delta T^(4)_ab dx^a dx^b
                    =exp(gamma v)R_gamma dv^2/(4pi r^2). (11)

It is finite at the future horizon. The local two-dimensional measure
restoration and four-dimensional Weyl/contact variations vanish at these
off-source readouts. The MASSIVE nonlocal stress does not vanish by this
argument and must still be added. The coefficient of the flux source in
the mass equation is 4pi r^2 delta T^r_v=exp(gamma v)R_gamma.
This reads the full tensor before imposing a metric gauge; differentiating
only the gauge-fixed two-field action could miss that flux equation.

For d(t,y)=A(t)D(y), with A supported in [0,T_A] and D in [a,b], the
incoming pulse at y arrives no earlier than t=y+a and ends by t=y+b+T_A.
Its integral in v is zero because integral P[D]dy=0. Thus the integrated
first-order zero-mode flux source returns to zero. This is not a full
mass-return result, a nonlinear energy balance or a second-order heating
calculation.

For canonical metric-force readouts within the source, also retain the
second vertex. In the affine psi,k chart, g_i=2d_i and
g_ij=2 speed_ij, with speed_ij from the preceding note. The response is

    delta F_i=integral_0^infinity 2d_i delta t_j dy
                              +integral_0^infinity 2t_0 speed_ij dy. (12)

The second term is an instantaneous contact. Physical covariant readouts
on the source additionally require the fixed measure/Weyl restoration;
(12) alone is not the full four-dimensional stress susceptibility.

## 5. Source review and next calculation

[Grosse and Murro, Remark 1.2 and Proposition 3.4](https://arxiv.org/pdf/1806.06544)
provide local reflecting-boundary finite propagation. The retained folded
projector's current cancellation verifies the relevant boundary hypothesis;
the source is not a charged-wormhole stability theorem.
[Bernard and Doyon, discussion of reflection and mean current](https://arxiv.org/pdf/1202.0239)
give chiral unfolding and the thermal normalization pi c/(12 beta_T^2).
The thermal chiral algebra is established mathematics. The evaluated
transfer here is its application to the retained ingoing metric vertex,
prepared state and spherical flux normalization, not a new CFT anomaly.

The focused checker independently integrates the Fermi occupation
difference, compares reflected transfer with retarded pulse quadrature,
checks the special-frequency cancellation without erasing the bulk term,
and checks the horizon graph coefficients and contacts. Numerical compact
profiles are piecewise C5, approximable by smooth profiles in the norms
used by these formulas. No Lean campaign is selected.

The next missing response is the massive angular sum with common local
subtraction, including physical outer traces and covariant horizon
matching. The zero-mode result supplies an exact normalization and causal
benchmark for that calculation; it cannot decide the sign or poles of the
total coupled operator. S168.3-S168.4 and every S167 task remain open.

Remaining to-do list: evaluate the massive response and full physical
domain; decide charged stability; establish nonlinear existence and
compatibility; complete exterior/state, charge/return, opening,
all-sector stability, transport, preparation, separate creation and audit.
