# Explicit quantum contacts in the coupled graph, with their normal matching

2026-10-24, E003/S169. Ordinary argument before frozen acceptance.
The full stable-wormhole objective is ACTIVE and INCOMPLETE; audit CLOSED.

This calculates the entire EXPLICIT LOCAL block W4+a0+a2+LLL from the
retained action, after its normal subtraction. It does not bound every
contact inside the remaining subtracted massive heat Hessian. In particular
the local mass-vertex chart conversion and the finite diagonal extension
of its remainder must remain in that remainder until explicitly matched.

## 1. Domain and exact scope of the split

Keep the actual finite seed Q, paired cubic C2,1 seam, thermal state,
regular bolt, physical mu0 and all angular levels. The action chart is
phi=-log z+nu, with shape b. The normalized action is I/W, with the same
negative W as current-finite-local-action-assembly.md. Put

 f=GC/(pi B2), u=1/C, Delta v=(Qv')', R=-Q'',
 L=R-2Delta b, M=Delta phi, V=Q(phi')^2,
 A=phi L^2/60+LM/30-M^2/10+LV/30-7MV/90-7V^2/180,
 Fe=4 gamma_E-2H_q+4log(mu0 B), k_mag=1/15-q^2/3.

The explicit local density is exactly

 L_exp/f = exp(-2b) A/2-u phi L/36-u V/36
             +k_mag u^2 phi exp(2b)/2
             +Fe exp(-2b)L^2/480+q u Q(b')^2/12
             +(c0/f)exp(2b)/4.                              (1)

Here c0=2G J0/(pi C B2) is the ONE shared vacuum coefficient. Keep its
known enclosure |c0|<=U0<2.39121e-25, not a guessed signed value. The
canonical a1 variation vanishes in this phase. The LLL is its full
kernel-free paired cocycle. Equation (1) is an exact action allocation,
not a finite-mode approximation or an adjustable forcing.

The original boundary representatives are retained. The W4 integration
by parts has second-variation endpoint flux O(z^3 log z) on mass/fast
directions; it vanishes at the actual bolt because Q(z_h)=0. The a2
integrated paired completion keeps the curvature-divergence seam term;
b(0)=b'(0)=0 removes its limiting Hessian boundary contribution. The a0
and LLL pieces have the already established fixed-source realization.
Their first-variation boundary work is not erased from the global source.

These facts give the local Hessian on the admitted mass/log phase without
a smooth reflected surrogate. They do not prove uniform X-ball admission.
For the graph estimate below, start with compact smooth graph inputs and
extend by the displayed bounded-map estimates, retaining all original
global matching terms separately. No new physical wall is imposed.

Let H_exp be the Euler Hessian of (1), A_exp=z^4 H_exp, and A_exp,0 its
coefficient of z^0, RETAINING its explicit log z. Define the remainder of
the full normal operator by subtraction, so

 A_actual-A_normal = (A_class-A_class,0)
                   +(A_exp-A_exp,0)
                   +[A_mass,rem-(A_normal-A_class,0-A_exp,0)].    (2)

This is an exact partition. A_exp,0 alone is not translation invariant.
Its outer logarithm cancels with the canonical massive normal term in the
last bracket, not with an invented finite contact. The already bounded
log-frequency principal component, separated input and coincident-domain
transfer belong to that remaining bracket and are allocated only once.

## 2. Differentiate before discarding powers or taking absolute values

Use D=z partial_z, sigma=4.005, a(s)=(s+3)/(3(s+1)), theta=nu+a(D)b.
The row operation is T(F)=(D-4)F_b-(D-6)F_nu/3, BEFORE output cutoff.
The block to bound is

 chi z^-sigma [1,0;-(D-6)/3,D-4](A_exp-A_exp,0)
       [1,-a(D);0,1] z^sigma diag(Jc_sigma^-1,B_sigma^-1).       (3)

Jc and B are the COMPLETE retained normal multipliers. The cutoff is
the same chi at a=.01, with derivative guards (1,3,10,60). Neither the
normalization nor the collar is changed.

One constructive exact evaluator is as follows. Treat (nu,nu',nu'') and
(b,b',b'') as independent jets; take the second jet derivatives of (1).
For a direction z^s, replace its kth ordinary derivative by (s)_k/z^k,
and differentiate a coefficient with partial_z+s/z. The Euler sum is
sum_(r,k)(-partial_z)^r L_(i,r;j,k) partial_z^k, with the FULL product
rule. Multiply by z^4, perform T and the input return, then extract powers.

After subtraction every entry is a finite sum

 f sum_(p=2)^12 sum_(ell=0)^1 z^p (log z)^ell c_ij,p,ell(D).    (4)

The scalar-input orders in the two rows are at most four and three.
After division by D+1 the shape-input orders are at most two and five,
with a possible constant rational remainder. The missing orders are
exact cancellations, not small-coefficient assumptions. The full
polynomial coefficients and the unsimplified normal pieces are retained
in the receipt. a0 contributes at p=4 and is not set to zero.

## 3. Verify normal compatibility and the finite fourth coefficient

This check precedes the numerical bound. Write w=G lambda/pi,
r_R=96pi G c_R lambda and ell=lambda0 B2/C=beta(1+11w/360), using the
actual vacuum relation. Then f=w/beta and D_R lambda=r_R-w/60.
The explicit quantum normal pieces in the first row/scalar-input column
and the two mixed transformed entries are

 -f beta^2 s(s-3)(3s^2-9s-1)/30,
 -f 11beta^2(s-3)/[45(s+1)],
 -f 11beta^2 s/45.                                      (5)

Adding the independently retained CLASSICAL normal pieces cancels both
mixed entries exactly and gives

 Jc=-6beta[1+11w/180-(r_R-w/60)s(s-3)](s-4)(s+1).

Thus the quantum contacts do not leave an omitted scalar/shape normal
coupling. The scalar normal inverse used in (3) is the correct one.

For the highest shape order put P4(s)=s(s-1)(s-2)(s-3),
K=GC beta^2/(15pi B2). The local coefficient is

 K[Fe/4+1/6-log z](s-4)P4(s),                         (6)

up to lower differential orders. No independent fourth coefficient may
be added to it. The complete angular flat model has

 K_flat(p)=C p^4[log(p/(mu0 ell_o))-31/30]/(30pi)
                       + E_flat(p),  ell_o=B/sqrt(C),
 |E_flat(p)|<=1000C(1+p^2)^(3/2).

These are the fixed-constant, entire-spectrum result of
outer-flat-tower-normal-remainder.md, not an assumption based on finitely
many angular modes. Its digamma series and sector asymptotic are
[DLMF 5.7.6](https://dlmf.nist.gov/5.7.E6) and
[5.11.2](https://dlmf.nist.gov/5.11.E2), rechecked 2026-10-24. They supply
the spectral identity, not a rough-boundary or curved-operator theorem.

In the physical normalization, after removing the already allocated a2
piece, the massive high-frequency normal coefficient is

 -K[log|xi_x|-log z+log(beta)/2-log(mu0 ell_o)-31/30+Fe/4].

Adding (6) gives -K[log|xi_x|-log(mu0 ell_A)-6/5],
ell_A=ell_o/sqrt(beta)=lambda^-1/2. Both Fe and log z cancel, and
-31/30-1/6=-6/5 agrees with the retained full Mellin normal family.
This matches the finite highest-order contact. It DOES NOT estimate
the curved lower-order kernel or its remaining diagonal extension.

## 4. All-frequency bounds with the physical crossover retained

Write L_big=10pi Fstar/(G lambda), delta=.005,
r=D_R lambda, A_v=Cv-r sigma(sigma-3)>0. The already proved lower bounds
are |B|>=K L_big |s-4||s||s-3|/2 at all frequencies, and for t>=4,

 |B(sigma+it)|>=K |s-4||s||s-3|
                         [L_big+t^2 log(2+t)]/2.             (7)

For |t|<=1 use the low constants from current-coupled-local-graph-bound.md;
for all frequencies their order-zero through order-three bounds are

 Cb=[2/(K L_big delta sigma(sigma-3)),
     2/(K L_big sigma(sigma-3)), 2/(K L_big sigma), 2/(K L_big)].

For the two additional powers needed here, set
T=sqrt(L_big/log L_big)>4 and use

 Cb4=max(Cb0,8/(K L_big),1/[K sqrt(L_big log6)]),
 Cb5=max(Cb0,2T^2/(K L_big),2/[K log(2+T)]).           (8)

The t^4 bound uses t/(L_big+t^2 log6)<=1/(2sqrt(L_big log6)).
For t^5 split at T: the classical part controls |t|<=T, while the full
logarithmic graph controls |t|>=T. This is a bound on EVERY frequency,
not numerical sampling near the crossover. No physical modes are cut.

The scalar constants for powers zero through four are

 Cj=[1/(6beta A_v delta(sigma+1)),
     1/(6beta A_v(sigma+1)), 1/(6beta A_v),
     1/(12beta sqrt(A_v r)), 1/(6beta r)].             (9)

For each c in (4), first perform division by s+1 if appropriate and
expand at s=sigma+it BEFORE taking absolute values. Bound its polynomial
by sum_k |c_k| Cb_k or Cj_k, and its rational remainder by
|rem| Cb0/(sigma+1). Denote the resulting multiplier bound b_ij,p,ell.

Put A_a=log(1/a) and M(p)=p^3+12p^2+51p+100. The H3 multiplication
bound for chi z^p is a^p M(p). For chi z^p log z it is

 a^p[A_a M(p)+M'(p)].                                (10)

Indeed p A_a>1 for p>=2, so sup_(x<=log a) e^(px)|x|=a^p A_a;
differentiate the product through three orders and use (1,3,10,60).
The input/output Sobolev norms use binomial derivative weights, giving
the exact product constant sum_(j=0)^3 binom(3,j)||alpha^(j)||_infinity.
Equations (8)--(10) prove each entry bound by a positive finite sum.
Take its maximum row sum for the joint scalar/shape H3 norm.

Freeze the target eta_exp<.004 and the sharper reporting guard <.0033.
Record all four entries, the whole-frequency constants, vacuum enclosure,
normal matching and contribution by z power and log power. This uses no
new local regularity assertion about the unknown nonlocal heat Hessian.

## 5. Decision and complete remaining obligations

If the executed bound passes, subtract eta_exp from the preceding .18698
allocation and select the conservative combined UNKNOWN target <=.18368
only if available. This target includes the residual canonical massive
response, all its unallocated chart/diagonal contacts, angular discreteness,
thermal contributions and other original-domain/global feedback. The
bounded far-input and identical-collar transfer are not added again.

RESEARCH DECISION: retain the construction and the full normal graph.
Use the exact contact/normal matching to isolate the remaining canonical
geometric heat response. Do not optimize the passing explicit contacts,
cold/scalar source, or remote transfer while that response remains unknown.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): obtain an actual strong
H3 bound for the commonly subtracted canonical massive response on the
identical-collar cylinder, after its full normal and already allocated
log-frequency principal terms are removed. Include the mass-vertex chart
conversion and finite lower contacts, and prove the differentiated common
regulator limit on the cubic paired realization. Test the combined remaining
feedback <=.18368. Existing weak and lower-output-order envelopes cannot
be differentiated to establish this; they are not a certificate or an
obstruction. Feasibility of that estimate remains undecided.

Joint affine <=7.49e-12 and known-source <=1e-14, independent Mellin
Efar<=1e-10, original mass/slow coefficients and the actual correction
remain UNKNOWN. The earlier sufficient projection screen is unchanged
and is not recomputed here. At rho=1e-8 the unclosed requirements are
Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1, with allocations Y0<=5e-9,
Z0<=.25, Z2<=5e7. Uniform nonlinear mass/log admission, exact compatibility,
physical matching, semiclassical validity/noise and all-sector stability
remain open. No controlled solution or nonexistence theorem follows.

Remaining to-do list: remaining canonical quantum graph bound and global
assembly; joint sources, original moments and independent far tail; actual
controlled correction and nonlinear admission; compatibility, matching,
validity/noise and all-sector stability; charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
