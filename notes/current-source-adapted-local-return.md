# Complete local two-row return of the source-adapted trial

2026-10-28 HST, E003/S169. Ordinary argument before execution. Previous
turn: progress, 8dc41ff bounded the entire source approximation and proved
the fractional-flow coefficient. Preserve its actual-domain trial, all
accepted results, P32, closed audit, norms and full stable-wormhole goal.
The complete remaining test is Ereturn+Kgeo*N_9/2+Eapplication<4.313e32.
This calculation attempts the ENTIRE local Ereturn, including both rows
and the original York adjoint. It does not evaluate N_9/2 or the later
actual-minus-reference quantum applications. Report the complete inverse
unresolved if those contributions are not bounded.

## Original row norms and the boundary adjoint

Let y=QJ/z^2 for an arbitrary original E_t test xi. Direct substitution
in the accepted York equation gives

 -y''+Pxi y=(4xi/z-xi')/z^2, y(0)=y(L)=0,
 Pxi=6/z^2+(Q''-4Q'/z)/Q,
 Uxi=xi/2+z y, Vxi=-xi-z y-z^2 y'/2.                 (1)

On the exact finite seed Pxi>=4.5/z^2 and Q Pxi>=1. To check these without
an endpoint singular division, expand z^2 Q(Pxi-4.5/z^2) and z^2 QPxi.
The leading terms in the first are 1.5beta-4.5u z^2+18beta z^3-7.5beta
z^4, so z<=1 and the retained tiny thermal/eps terms leave it positive.
The second is 6beta+2(3beta+u+eps)z^4-12eps z^5+12eps z^6; division
by z^2 and z^-2+z^2>=2 give QPxi>1.9. The code checks the coefficient
guards. The singular endpoints use the ORIGINAL Friedrichs/form limit,
not a new wall: the mass/log rates and regular bolt make all fluxes vanish.

With cY=1+2/sqrt(4.5), testing (1) and integrating its xi' once gives

 Ey=(int[y'^2+Pxi y^2])^1/2 <=cY ||xi||Hmu <=cY Cxi ||xi||Et,
 ||Uxi||Hmu<=CU ||xi||Et, CU=Cxi(1/2+cY/sqrt(4.5)),
 ||Vxi||Hmu<=CV ||xi||Et, CV=Cxi(1+cY/sqrt(4.5)+cY/2), Cxi=1.01. (2)

These are fixed-seed row comparisons, not C5 output or independent phase
estimates. They retain the original z^-4 measure and global York return.

Let fn,fu be the classical comparison's complete canonical strong rows,
evaluated at nu0=xi0/2, u0=v0 from 8dc41ff. They include Einstein, Maxwell
and R2 terms with their exact parameters. Their model congruence is
fn/2+V0 fu=K(H)xi0, V0=-3/2+a R_(2a). Let rf=R_(2a)fu and define the
EXPLICIT original-domain adjoint trial

    Y=(z^2/2) Q2 partial_z rf.                       (3)

It vanishes at both original endpoints. At the bolt the flux has at
worst the permitted power/log form and Y has finite Dirichlet energy;
at zero Y=O(z^2), so the original boundary pairing vanishes. No independent
mass/log coefficient is prescribed. For the exact original y in (1),

 <fn,Uxi>+<fu,Vxi>-<K(H)xi0,xi>
   =<A,xi>+<R,Ytest>, Ytest=y,
 z^2 A=2z Q2 rf_z,
 R=z fn+2z h'-2h-(z^2/2)[(Q2/Q)(Q''-4Q'/z)-2a]rf_z,
 h=Q2 rf_z, h'=2a rf-fu.                             (4)

This identity uses BOTH rows before norms and the same original endpoints.
It includes the entire finite-domain adjoint return, not a local York
surrogate. Since ||y/sqrt(Q)||2<=Ey, its complete dual cost is bounded by
Cxi||z^2 A||2+cY Cxi||sqrt(Q)R||2.

On 0<s=L-z<L, the actual ratio .1<Q/Q2<2 follows from a in [.999,1.001],
the cubic coefficient in [-1.34,-1.32], the quartic in [.49,.51], and the
tiny fifth/sixth coefficients; split s<=.05 and s>=.05. These guards are
checked on exact root-based coefficients. Taylor bounds with M3=1000 give

 z^2|Q''-4Q'/z-2a Q/Q2| <=8gamma+(4M3/3+8a+2M3 L)s.

For a model profile f define Nf=||f||2, Mf=||s f||2, Hf=||Q2 f_s||2
on [0,L], and Gf=||sqrt(Q2)f_s||2. Then

 Wf=sqrt(2gamma)*sqrt(Nf Mf)+sqrt(a) Mf >=||sqrt(Q2)f||2,
 ||(gamma+s)sqrt(Q2)f_s||2<=2gamma Gf+2Hf/sqrt(a).

Use these in (4), Q<=.17, z<=1 and the ratio guards to bound the whole
adjoint contribution explicitly. No unknown adjoint norm is replaced
by a conditional ceiling.

## Actual field and exact strong local rows

Set x=s/gamma, q2=2x+a*x^2, q=Q/gamma^2, Rq=q/q2 and z=L-gamma*x.
Let f=v0, r=R_muY f, W(z)=int_0^z (tQ''-4Q')'r dt and
M=z q_xx+4gamma q_x. The actual field of 8dc41ff is EXACTLY

 nu=z^4/(3L^5)[-z Rq f+(z Rq muY-M)r
                              -(7gamma Rq-z Rq_x)h_r+W],
 u=z^4/L^5[z Rq f+(2M/3-z Rq muY)r
                              +(5gamma Rq-z Rq_x)h_r-2W/3],
 h_r=q2 r_x, W_x=(z q_xxx+3gamma q_xx)r.              (5)

There is no numerical York solve or cutoff interpolation in (5). Original
physical transitions, boundary conditions and bolt behavior are unchanged.
The nonzero rho^4 log(rho) term is retained. Four weak local derivatives
are permitted: the fourth strong rows have at most a logarithmic bolt
singularity in L2. Their integration-by-parts boundary fluxes vanish.
The failed fifth-derivative source hypothesis is not used.

The implementation differentiates the exact Euler rows of this generator,
which is algebraically equivalent to its Hessian and adjoints below. In
particular it uses the independent trace row in the accepted local-action
note. This retains the second gradient contact through the variation of
N_phi=Q phi'^2; it does not linearize that square before varying the row.
Its original finite epsilon coefficient is of order 1e-29. The fifth and
sixth root-based polynomial coefficients are bounded by 1e-25, not zero.
The first development receipt preserves a rejected, unnecessarily small
1e-50 guard for those coefficients; no result is claimed from that run.

Use the EXACT affine action generator of current-finite-local-action-
assembly.md. Evaluate its Hessian in the five variables
(nu,u,Delta_q nu,Delta_q u,h_nu), with h=q2 partial_x,
Delta_q=Rq Delta2+Rq_x h. At the background

 L0=-q_xx, M0=gamma q_x/z+gamma^2 q/z^2, N0=gamma^2 q/z^2,
 delta L=-2 Delta_q u, delta M=Delta_q nu,
 delta N=2gamma Rq h_nu/z.

Its second gradient contact is delta2 N=2Rq h_nu h_eta/q2. If SN is
the first derivative of the unvaried generator with respect to N, its
strong nu row is -2[(SN Rq)Delta2 nu+(SN Rq)'h_nu]. This contact MUST
be added to the five-variable Hessian. The remaining adjoints are
Delta_q*=Delta_q and h*=-h-q2'. Divide the action Hessian by -beta to
obtain the positive TT-energy convention. Source densities are integrated
against ORIGINAL dz=gamma dx.

The algebra closes on r,h_r,f,h_f,g,h_g,k,W, where g=Hf, k=H^2 f:
Delta2 r=muY r-f; Delta2 h_r=muY h_r-h_f+q2'(muY r-f);
Delta2 f=-g; Delta2 h_f=-h_g-q2'g; Delta2 g=-k.
For coefficients A, Delta2(A b)=(q2 A''+q2'A')b+2A'h_b+A Delta2 b.
Also Delta2 W=D h_r+(q2D'+q2'D)r, D=z q_xxx+3gamma q_xx.
Thus every retained strong row is an explicit finite coefficient vector;
no strong r derivative is bounded by dividing by the tiny bolt radius.

Subtract the comparison classical rows BEFORE taking norms. The actual
classical difference is bounded through (2) by the original L2 norms of
z^2(fn_actual-fn) and z^2(fu_actual-fu). Add the ENTIRE Weyl row, exact
LLL, common a2 and shared bounded a0 contributions. Their signs and
coincidence/density contacts follow the accepted local assembly. The
retained nonlocal geodesic term is separate and remains in the final budget.
No actual thermal/quantum difference from the reference is silently
removed; its later applications remain obligations.

## Full profile enclosures without near-coincident-pole losses

All comparison profiles are rational functions F(H)delta_b, multiplied
by At. The shape multiplier is

 Fv(p)=K0 [1+(Cv-1)lambda_n/(p+lambda_n)]
                         /[(p+muT)(p+muS)(p+lambda_l)].

Use its full product, including Cv-1. The r profile adds R_muY; g and k
multiply by p and p^2. The model canonical rows and rf use the exact
classical Hessian polynomials, with leading p^2 cancellations performed
algebraically. Their rational functions decay at least as 1/p; there is
no omitted polynomial/delta term. The poles near 2 are NOT separated
into coefficients divided by their tiny differences. Enclose all of them
by one circle |mu-2|=.25 and use Cauchy's residue formula on -F(-mu)Q_mu.
Its norm contribution is at most .25 sup_circle|F(-mu)| times that of Q_.9.
Indeed Re k(mu)>.9 on this circle, using the principal square root and
a in its proved range. The full positive Q integral gives |Q_k|<=Q_.9;
the differentiated integral bounds the flux by twice its real .9 counterpart.
The remaining poles muS, muY and lambda_n are isolated and retained with
their complete residue balls. This is a geometric profile calculation,
not a finite selection of quantum angular or thermal modes.

For each real comparison degree k>0, reuse the established complete Q
norm, moment and flux bounds. On a cell s in [s0,s1], use the entire
integrated logarithmic bound when a*s1/gamma<=1, and monotonicity with
Q_k<=Q_0 at the left endpoint otherwise; for a*s0/gamma>=1 also use
Q_k<=pi(a*s0/gamma)^(-k-1). The flux is <=gamma and has tail
<=pi gamma(k+1)(a*s0/gamma)^(-k). Multiply these full-cell bounds by
residue/contour weights. The spectral norm from 8dc41ff additionally
gives N(r),N(f),N(g),N(k)<=At sqrt(gamma) Nmodel, and the same
gradient bound for r,f,g. Use the smaller of these proven bounds on
each cell, retaining the actual gamma measure.

W has the entire bound |W(s)|<=Dmax int_s^L |r(t)|dt, with Dmax bounded
from the exact finite polynomial. For a Q_k term the integral is at most
gamma/[a k(k+1)] and, beyond a*s/gamma>=1, at most
pi gamma(a*s/gamma)^(-k)/(a*k). These retain the global return.

All real cells cover z in [.01,L]. Coefficients are enclosed on whole
cells by outward Taylor-jet interval arithmetic; they are not sampled.
For 0<=z<=.01, the weighted strong coefficients are holomorphic and
vanish to at least order two, except an explicit term linear in log z.
This follows by substituting nu,u=O(z^4) before taking the two local
Euler variations. The classical terms are analytic. In the Weyl term
the only remaining logarithm is the original phi0=-log z; its derivatives
are rational. Evaluate the constant and linear-log coefficients separately
on the entire circle |z|=.02, enclosing every arc, and use the maximum
principle. The circle is only an analytic bound for the original one-sided
coefficients; it does not smooth or change the paired operator. Integrate
the full outer logarithm using int_0^b log^2 z dz=b[(log b)^2-2log b+2].
No original boundary region or numerical application error is assigned zero.

The computed Ereturn is the complete classical row difference, full local
quantum rows and (4)'s adjoint return, all in E_t*. Its role in the complete
test is explicit. Report failure if the full expression exceeds 4.313e32,
and keep N_9/2 unknown unless separately bounded on this actual field.
A local pass alone is not a first inverse or original admission certificate.

Remaining to-do list: complete joint first-reference residual, later
actual quantum applications and G_t; every coupled/admission requirement;
existence, compatibility, matching, validity/noise and all-sector stability.
