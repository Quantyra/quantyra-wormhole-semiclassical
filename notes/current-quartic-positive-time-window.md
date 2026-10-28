# Finite-window contract for the actual quartic quantum application

2026-10-27, E003/S169. The preceding turn made progress: it proved the
fixed quartic column's full energy-dual regulator estimate. The active
goal remains the COMPLETE residual T_t-(I+T_t)F_n, including numerical
application errors. This note reduces one remaining application to a
finite heat-time/angular window with explicit certified tails. It does
not evaluate that window or complete the goal.

## 1. Keep the same operator and application error

Use precisely the paired mass/log realization and original E_t, beta,
action factor, York return, fixed period, scale and regular bolt from
current-quartic-column-energy-dual.md. The trial is psi=z^4, J=0,
v=Vpsi=-z^4. The two original equations, all cutoff/phase conditions,
exact UV logarithm, coincidence contacts and finite returns remain in
the complete coupled operator. No seed-to-neighborhood transfer is made.

Let r_v(t,z) be the ORIGINAL canonical geometric heat-Hessian density
after the completed a0/a1/a2 subtraction. For a core h=Vw,

    D2 r(t)[h,v]=int_0^L h(z) r_v(t,z) dz,
    Q_v[w]=-p integral_0^infinity Theta(t) D2 r(t)[h,v] dt/t,
    p=G/(B^2 beta), G=1,
    Theta(t)=sum_(n>=1) d_n exp(-m_n^2 t),
    d_n=2(q+2n), m_n^2=n(n+q)/C, q=2000^4.

The fixed-column regulator proof defines this functional for every
w in E_t, not merely the selected core. The original local and LLL terms
are separate retained allocations. The product-centered heat remainder
is NOT silently substituted for r_v.

For t>=t0=.001, the proved full thermal/finite-domain heat estimate is

    |D2 r(t)[h,v]| <= [CH+C0/t+C2 t] ||w||E_t,
    CH=64 M CLt/3, C0=2 CLt/(3pi),
    C2=CLt ||P4[v]||L2(0,1)/(480pi), M<734.

The original energy-to-L2 constant CLt and polynomial P4[v] come from
the preceding full-column argument, including its contact integrations
by parts and actual boundary conditions. This bound covers the complete
quantum action on the indicated time region, not selected thermal modes.

## 2. Entire omitted heat-time and angular tails

Set a=(q+1)/C>324. Since m_n^2>=a for every n>=1, the complete counting
bound Theta(t0)<=4C/t0 gives, for T>=t0,

    d_time(T) = p (4C/t0) exp[-a(T-t0)]/a
                     [CH/T+C0/T^2+C2].                    (1)

This bounds the ENTIRE integral on t>=T in E_t*. The two decreasing
factors 1/t and 1/t^2 were bounded at T; the contact term is integrated
with the same exponential. Zero modes of the orbit operator cause no
infrared problem because the massive angular gap remains.

For a finite retained angular range 1<=n<=N, m_n^2>=a n, and the whole
discarded counting function obeys, for t>=t0,

    Theta_>N(t) <= A_N exp[-a(N+1)t],
    A_N=(2q+4(N+1))/(1-x)+4x/(1-x)^2, x=exp(-a t0).

This is the exact sum of a positive geometric-series majorant, including
the growing degeneracies. Integrating it gives

    d_ang(N)=p A_N exp[-a(N+1)t0]/[a(N+1)]
                         [CH/t0+C0/t0^2+C2].               (2)

The two estimates overlap on n>N,t>T; counting that region twice is
conservative. There is no discarded angular mode or omitted Matsubara
tail. The retained heat kernel still contains ALL thermal/radial modes,
the paired boundary and regular-bolt returns.

## 3. What an actual retained-window evaluator has to certify

Define Theta_N(t) by its exact finite sum with the actual m_n, and

    Z_window(z)=-p integral_(t0)^T Theta_N(t) r_v(t,z) dt/t.

If a validated numerical evaluator returns a spatial density Zhat with
||Z_window-Zhat||L2<=e_Z, its column error is at most CLt e_Z. This is
a condition on the error, not on the unknown signed value of Z_window.
Qualitative smoothing does not establish e_Z.

For a time-dependent density approximation, a uniform L2 kernel error
rho and a separate weighted time/quadrature/spatial error e_num give

    e_Z <= p W_N rho+e_num,
    W_N=integral_(t0)^T Theta_N(t)dt/t.                     (3)

To bound W_N numerically, take 256 geometric time cells and integrate
the decreasing positive integrand by its left-endpoint value times the
exact cell length, using interval arithmetic. This is a rigorous upper
sum for an elementary scalar weight, not a quadrature for r_v. Every
coefficient and mass in Theta_N is retained. No normalized error is
identified with a physical output tolerance without CLt and the original
coefficient/output factors.

For the SPECIFIED one-dimensional exact Ritz sensitivity, let alpha be
the proved coefficient envelope and d_cut the cap plus retained-transport
debit from the preceding assembly. The sufficient application test is

    d_cut+alpha[d_time+d_ang+CLt(p W_N rho+e_num)] <=1e-30.  (4)

All arithmetic in (4) is executed for T=.05,N=60. The purpose is to
quantify the remaining evaluator requirement. Neither rho nor e_num is
declared attained. In a complete multicolumn F_n, (4) is replaced by its
actual coefficient map and the sum of ALL columns, source, scalar-solve
and coefficient errors. Those comparisons remain unknown.

## 4. Attempted implementation and the first unsupported step

The existing current_subtracted_thermal_action.py does NOT return r_v or
Z_window. Inspection confirms that its directions are compactly supported
in [.5,.875], its initial directional sensitivities vanish, and its
certified mass range is 18<M<26 for angular levels n=1,2. Substituting
v=-z^4 breaks the zero initial-sensitivity argument. The actual v is
nonzero at the regular bolt. Using the same startup values would omit
boundary/normalization returns. It also evaluates a subtracted mass-vertex
characteristic, whose geometric Weyl conversion and common local
restoration have not been proved for this full-domain application.

The actual-disk parametrix proof is also not an implemented evaluator:
its boundary-compatible P contains the exact capped heat kernel, while
the code stores positive majorants for transport and propagation, not
signed kernel values. A Schur or pointwise remainder table cannot be used
as those values. Extending its t<=.001 constants to .05 would be another
unsupported implication. No such extension or smooth surrogate is used.

Thus the first unexecuted step in the direct method is a certified
geometric heat-Hessian density on [t0,T] in the stated L2 norm, preserving
the actual paired realization. The signed midpoint and its numerical
error remain UNKNOWN. The complete positive-time calculation was
attempted through both available implementations, but neither supplies
this map. There is no new bound on the full coupled residual this turn.

## 5. Decision and scope

Execute (1)-(4), retain the tried T,N values and report the certified
tail debit and numerical allowance. This establishes a finite-window
accuracy contract; it is not a completed signed-evaluation milestone.
It does not improve the accepted all-time zero-centered column bound
until the retained-window error is actually certified.

Retain the direct paired geometric method. Against the same obstacle,
the raw thermal ODE alternative additionally needs its geometric action
identity and nonzero-bolt variation normalization, so simply increasing
its angular range or mesh is not an adequate next step. The single next
calculation is the actual Z_window enclosure in (3), with a signed center
and proved spatial/time/kernel errors satisfying (4). This is proposed,
not executed. P32, the original norms, the closed audit and all accepted
results remain unchanged; no degree, cold or scalar optimization follows.

Remaining to-do list: actual positive-time application and full residual;
all source/solve/application errors, signed G_t and every original coupled
and admission budget; nonlinear existence/contraction, exact compatibility,
physical matching, semiclassical validity/noise and all-sector stability;
transport, preparation, separate creation and the full objective.
