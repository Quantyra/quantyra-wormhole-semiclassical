# A heat cutoff cannot discard the logarithmic principal form

2026-10-23, E003/S169. Ordinary argument and source review before execution.
Previous goal turn: progress. The stable semiclassical wormhole goal remains
incomplete and the quantitative feasibility audit stays closed.

The proposed differentiated heat bound must be tested before pursuing a
seam estimate. There is a prior obstruction: after subtracting only the
local a0/a1/a2 terms, the omitted small-time part still carries the entire
fourth-order logarithmic principal symbol. Its norm does NOT tend to zero
on the full physical form space as the heat cutoff tends to zero. This
holds on smooth compact interiors, independently of the outer seam.

This corrects the possible operator-norm use of equations (11)-(13) in
[the preceding assembly](current-canonical-heat-restoration-assembly.md).
Its local rows, corrected scale, finite-entry conditional bound and general
two-residual algebra remain valid. A fixed-direction remainder estimate is
not disproved. The stronger proposed Aop bound for the UNREMOVED complete
nonlocal principal part is impossible, not merely an unknown constant.

## 1. Positive heat representation of the complete angular principal model

Use exactly the current action, rank, angular multiplicities, fixed mu0,
actual B,C,q and physical phase. Retain the complete principal multiplier
K from [the angular calculation](current-all-angular-principal-operator.md):

    K(p)=C p^4 kappa(p)/(120pi), kappa(0)=kappa0=-F_eff>0,
    M_n^2=n(n+q)/C, d_n=2(q+2n), a=alpha(1-alpha),
    K_rem(p)=sum_n d_n M_n^2/(2pi) integral_0^1
                  R(a p^2/M_n^2) d alpha,
    R(x)=log(1+x)-x+x^2/2.

Define r(s)=1-exp(-s)-s+s^2/2. Since r(0)=r'(0)=r''(0)=0 and
r'''(s)=exp(-s)>0 for s>=0,

    0<=r(s)<=min(s^2/2,s^3/6),
    R(x)=integral_0^infinity exp(-u) r(ux) du/u.       (1)

The second identity follows by differentiating in x and using R(0)=0;
the derivative is 1/(1+x)-1+x. All terms are integrable. With
Theta1(t)=sum d_n M_n^2 exp(-M_n^2 t), Tonelli gives the EXACT positive
representation

    K_rem(p)=(1/(2pi)) integral_0^1 integral_0^infinity
                         Theta1(t) r(a p^2 t) dt/t d alpha. (2)

This is a representation of the full angular PRINCIPAL MODEL. It does
not assert that a raw fixed-metric mass determinant equals the full curved
geometric determinant. Its logarithmic symbol has the physical compact-
interior identification already proved in the retained heat-Hessian note.

For any fixed numerical heat cutoff T>0, let K_<T be (2) on 0<t<T, and
let K_keep,T contain the SAME local fourth contact and (2) on t>=T.
Thus K=K_<T+K_keep,T with both parts nonnegative. No physical cutoff or
change of finite coupling is introduced.

Using int_0^1 a^2 d alpha=1/30 and (1),

    K_keep,T(p)<=C p^4/(120pi) [kappa0+J(T)],
    J(T)=C^-1 sum_n d_n exp(-M_n^2 T)(T+M_n^-2)
        =C^-1[T Theta(T)+integral_T^infinity Theta(t)dt]. (3)

The last identity retains the complete angular tower. The fixed gap
M_n>=18 and the known counting bound N(lambda)<=4C lambda, with
N(lambda)=0 for lambda<324, imply

    Theta(t)<=4C exp(-324t)(324+1/t),
    J(T)<=Jbar(T)=4[(324T+2)exp(-324T)+E1(324T)],
    E1(x)=integral_x^infinity exp(-s)ds/s.             (4)

Derive the first line from Theta(t)=t integral N(lambda)exp(-t lambda)
d lambda. Equation (4) is a finite known allowance for every T>0.
It is used to give counterexamples, not optimized as a source bound.

## 2. The decisive uniform feasibility test

The retained complete-angular formula has

    kappa(p)=4log[p/(mu0 ell_o)]-62/15+o(1),
    ell_o=B/sqrt(C),
    kappa(p)>=4log[p/(mu0 ell_o)]-152/15.             (5)

Therefore for EVERY fixed T>0,

    0<=K_<T/K<=1,
    lim_(p->infinity) K_<T/K=1,
    sup_(p>0) K_<T/K=1.                              (6)

This is a lower-bound obstruction, not a failure of an oversized upper
estimate. The fourth local contact and any finite number of additional
fourth-order local terms cannot restore the missing log p.

For an explicit witness exceeding an omitted fraction eta in (0,1), put

    ceiling(T)=kappa0+Jbar(T),
    p_eta(T)=mu0 ell_o exp[38/15+ceiling(T)/(4(1-eta))]. (7)

At p>p_eta, (3)-(5) imply K_<T/K>eta. Execute eta=1/4 at
T=1e-2,1e-8,1e-16 using p=2p_eta and the existing all-angular interval
quadrature. Preserve all nodes and enclose E1. The samples illustrate the
all-T proof (6); they do not establish it by sampling. These finite p
values are for the principal model, not thresholds for the full Hessian
with its lower-order terms or claims of physical instability.

There is no contradiction with convergence at a fixed p. The other
bound in (1), together with Theta1(t)<=12C/t^2, gives

    K_<T(p)<=C p^6 T/(140pi),
    K_<T(p)/K(p)<=6 p^2 T/[7 kappa(p)].               (8)

Thus a fixed frequency band or extra regularity can give convergence.
The order of quantifiers matters: (8) is not uniform as p tends to
infinity. A physical ultraviolet band would need its own justification;
none is added to the current uncut theory or goal.

## 3. Why the actual paired domain cannot remove this obstruction

The actual physical TT diagonal form is l, and
||(s,xi)||_V^2=a(s,s)+l(xi,xi). Work in a smooth compact annulus strictly
inside the actual finite geometry. The retained compact-interior symbol
argument gives a nonzero positive p^4 log p term in l. Magnetic/state,
Weyl finite fourth contacts and local curvature terms do not change its
coefficient. The mixed scalar/TT terms have lower principal order.

The spectral-symbol source is Moller, Theorem 1 and equation (1.29),
[arXiv:0705.3857v3](https://arxiv.org/html/0705.3857v3).
It assumes a smooth closed manifold with locally stable kernel dimension.
Here it supplies only the local symbol: nested interior cutoffs and the
short-time direct parametrix identify that symbol on the actual annulus.
Separated boundary/state terms cannot alter its simultaneous diagonal
singularity. This is the retained ordinary transfer argument in
stationary-shape-logarithmic-response.md, not an application of a smooth
closed theorem across the C2,1 seam or at a moving boundary.

One can choose genuine admissible TT oscillatory fields without inventing
a boundary condition. In a physical radial coordinate take
xi_p=chi cos(p r)-c_p chi0, with both cutoffs compact in this annulus,

    w(r)=r^2 F'(r)-2rF(r),
    c_p=integral w chi cos(pr)dr / integral w chi0 dr,
    b_p=-(r^2F)^-1 integral_(r_bolt)^r w xi_p dr.      (9)

Choose a smaller annulus with w not identically zero and a chi0 with
nonzero denominator. Such a region exists for the actual rational F;
w is not identically zero. Repeated integration by parts gives
c_p=O(p^-N) for every N. The moment in (9) is exactly zero, so b_p and
xi_p vanish near BOTH endpoints and satisfy the TT constraint. Its
principal tensor has eigenvalues (-2xi_p,0,xi_p,xi_p), with nonzero
spin-two logarithmic coefficient. The O(1/p) radial component and
the negligible moment correction do not change the leading form.

For fixed T>0 the t>=T geometric Duhamel Hessian, including the second
vertex, has a finite bound C_T ||xi||_H2^2 on this annulus. Here is why
one must not assume it has the missing logarithm. Each first variation
L_i has at most two derivatives of its metric coefficient; the second
vertex has at most two total derivatives. In the two-vertex trace split
the simplex at 1/2. On each half one heat factor has time at least t/2.
Using a spectral expansion and Cauchy--Schwarz, the other heat factor
is bounded by one; summing its matrix elements gives ||L_i phi_m||^2.
Interior elliptic bounds on the smooth annulus and the long heat factor
bound this by C_t||xi||_H2^2. The constants are uniformly finite for
t>=T; the massive angular gap controls the remaining large-time integral.
The differentiated a0/a1/a2 subtractions are also bounded by H2 there.
This argument allows C_T to diverge as T goes to zero; no numerical
value or smooth-seam theorem is required to establish the obstruction.

Consequently the retained t>=T piece and local terms are O(p^4), while
l(xi_p,xi_p) has positive leading p^4 log p. If R_<T is the actual
omitted geometric short-time form after only those local subtractions,

    <xi_p,R_<T xi_p>/l(xi_p,xi_p) -> 1,
    ||R_<T||_(V to V*) >= 1                         (10)

whenever that bounded-operator norm is defined; if it is unbounded the
small-norm proposal fails more strongly. Physical action signs are used
in (10); changing to the negative normalized action multiplies numerator
and denominator consistently. No temperature or radius rescaling repairs
the ratio. This compact-core obstruction precedes the seam and adjoint
regularity questions. It does not require accepting an unverified global
smooth pseudodifferential calculus.

A putative bound ||D^2 r(t)||_(V to V*)<=Aop t^(1+alpha), alpha>0,
for all sufficiently small t, would imply ||R_<T||<=4GC Aop T^alpha/
(alpha B^2) ->0. This contradicts (10). Thus the corresponding Aop in
the preceding checkpoint cannot be finite for the complete unremoved
principal remainder. Its finite-entry prefactor 2.34375e-32 does not
overcome this loss. The old conditional finite-direction calculation is
unaffected; a uniform V-to-V* interpretation is ruled out.

## 4. Required change of method and remaining physical test

Keep the complete angular logarithmic principal operator analytically,
with its actual metric/York symbol and EXACT finite fourth contact.
The already executed all-angular quadrature provides a uniform relative
model error below 2e-5. It is still a model, not the physical remainder
bound. A bare scalar Laplace spectral multiplier is not silently equated
with the actual curved Hessian. The contact-complete bolt parametrix and
outer normal phase remain necessary.

Only after that subtraction should one seek a small curved/seam
remainder. The claimed two-order gain in the retained finite-kernel
construction still needs quantitative constants on the actual paired
mass/log phase. The single next calculation is that actual-minus-complete-
principal remainder in the physical form norm, including its paired seam
and the logarithmic bolt-adjoint extension. This is PROPOSED, not executed
by the present counterexample. Further source optimization remains parked.

Direct bounds on the chosen approximants can avoid an unnecessary
uniform small-remainder requirement. For symmetric R=H-H0, put
eta_e>=||R e_app||_V*, eta_psi>=||R psi_app||_V*, with the preceding
definitions E,P,D0,Db0,Q0. Then a sufficient physical bolt test is

    |Q0|+min(P eta_e,E eta_psi)
                    +(Db0+eta_psi)(D0+eta_e)<1.667e-9. (11)

Symmetry gives the two bounds on the same signed pairing. Actual
directional residual bounds, E,P,D0,Db0,Q0 remain unknown. This option
does not allow an unchecked high-frequency tail to be dropped, or an
extra smoothness assumption on the bolt adjoint.

Research decision: reject local-subtraction-plus-heat-cutoff as a
uniform small-form approximation. Continue with the complete logarithmic
principal part and the physical weak primal-dual residual method. This
uniform feasibility test IS EXECUTED; the full physical <1.667e-9 test,
exterior <6.67e-10, Y0<=5e-9, admission, nonlinear neighborhood and exact
compatibility remain UNEXECUTED. No wormhole nonexistence, physical
instability or completed construction is inferred.

Remaining to-do list: quantitative actual-minus-principal paired remainder
and physical forward/adjoint residuals; X/admission, nonlinear closure and
exact compatibility; matching, causal realization, semiclassical validity/
noise, all-sector stability, charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
