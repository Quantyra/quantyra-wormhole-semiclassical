# All-angular principal operator and a finite-cutoff obstruction

2026-10-23, E003/S169. Ordinary argument and source review before execution.
The full stable semiclassical wormhole objective remains incomplete.

This resumes the paused physical weak primal-dual calculation. A finite
angular approximation has a domain/order problem before its source error
can be optimized. We settle that problem and give an executable, uniformly
controlled replacement for the principal spectral model. The full curved
operator, its common remainder and the physical correction are not computed
by this package. The closed feasibility audit remains closed.

## 1. Same action, finite part and comparison operator

Keep c_R=1e22, G=1 as the unit, q=2000^4, the actual finite-seed C and B,
mu0=1/b0, paired state and mass/log boundary phase. No cutoff is introduced
into the physical theory. Use the rank-two normalization of
[the curved TT argument](stationary-curved-quantum-tt.md), section 7:

    A_n=n(n+q), M_n^2=A_n/C, d_n=2(q+2n), a=alpha(1-alpha),
    L=log(mu0 B), F0=4 gamma_E-2 H_q, F_eff=F0+4L,
    R(x)=log(1+x)-x+x^2/2,
    K(p)=-C F_eff p^4/(120 pi)
       +sum_n d_n M_n^2/(2pi) integral_0^1 R(a p^2/M_n^2) d alpha. (1)

K is a scalar comparison multiplier, not the complete Hessian on the
physical seed. Its a2 principal restoration uses the original finite part.
Vacuum, curvature contacts, Weyl, boundary and chart-return terms must
still be retained in the actual operator. Dividing by C p^4/(120pi) below
defines a coefficient for computation; it does not improve a physical norm.

The retained proof identifies the same positive p^4 log p principal
behavior on smooth compact interiors of the actual TT problem. Its
lower-order curved and remote terms do not erase that behavior. The
regular bolt is an interior point in Cartesian disk coordinates. No
time-Fourier cutoff at a punctured polar origin is used in this model.

## 2. Sum the complete angular principal part before integration

For t>=0, define

    nu_+/-=(q +/- sqrt(q^2-4t))/2,
    F(t)=-2[psi(1+nu_+)+psi(1+nu_-)]; F(0)=F0.          (2)

For t>q^2/4 the roots are conjugates with real part q/2. Below that
point both are nonnegative. Thus all digamma arguments have real part
at least one; there are no poles. Root exchange leaves F invariant.
Evaluate the smaller root as 2t/(q+sqrt(q^2-4t)) to avoid cancellation.
At the apparent coalescence F=-4 psi(1+q/2); it is not a physical branch
singularity. Complex balls enclose the real sum through this point.

Partial fractions of (2), using the convergent digamma difference, give

    F0-F(t)=sum_n d_n t/[A_n(A_n+t)] >=0,
    -F'(t)=sum_n d_n/(A_n+t)^2.                         (3)

The positive convergent sums justify differentiation on compact t ranges.
Twice differentiating (1) in z=p^2 and then integrating back with
K(0)=K_z(0)=0 gives

    K_z(z)=-C z/(2pi) integral_0^1 a^2[F(Caz)+4L]d alpha,
    K(p)=C p^4/(120pi) kappa(p),
    kappa(p)=(15/4) integral_0^1 y sqrt(1-y)
                         [-F(Cp^2 y/4)-4L]dy.           (4)

To obtain the single integral, write z integration as z v and substitute
x=a v. The length of {alpha:a>=x} is sqrt(1-4x); then y=4x.
The weight integrates to 4/15, hence kappa(0)=-F_eff. This fixes every
factor against the original fourth contact, including the spin count.

The source formulas are [DLMF 5.9.13 and 5.9.16](https://dlmf.nist.gov/5.9),
applied only at Re(1+nu)>0, and [5.11.2](https://dlmf.nist.gov/5.11#E2),
whose sector condition holds for these roots. They justify the digamma
identities/asymptotic, not the physical domain or any quantitative inverse.

## 3. Uniform quadrature bound, including arbitrarily large p

Let N(A)=sum_(A_n<=A)d_n. At A=A_n it is 2n(q+n+1)<=4A_n;
between levels N(A)<=4A. For g(s)=-F(exp(s))-4L, (3) implies

    0<=g'(s)=sum_n d_n t/(A_n+t)^2<=4,
    |g''(s)|<=sum_n d_n t/(A_n+t)^2<=4, t=exp(s).      (5)

The upper bound follows by Stieltjes integration:
sum d_n t/(A_n+t)^2 = integral_0^infinity 2t N(A)/(A+t)^3 dA
<=8t integral A/(A+t)^3 dA=4. Termwise the absolute logarithmic
derivative of t/(A+t)^2 is at most t/(A+t)^2. This establishes the
second bound without a smooth-seam assumption or a high-mass expansion.

Use y_j=2^(-j/16), j=0,...,320, so y_min=2^-20 and each log interval
has length Delta=log(2)/16. Interpolate g linearly in log y on each
interval. Its pointwise error is at most Delta^2/2 by (5). Integrate
this interpolant with EXACT weight moments W=int y sqrt(1-y)dy and
J=int y sqrt(1-y)log(y)dy. With r=sqrt(1-y), primitives are

    W(y)=-(2/3)r^3+(2/5)r^5,
    J(y)=-2(r^3/3-r^5/5)log y-4r^5/25
                         +8r^3/45+8r/15-8atanh(r)/15. (6)

Use W(0)=-4/15 and J(1)=W(1)=0 by limits. On [a,b], the two
interpolation weights are (log(b)W_ab-J_ab)/Delta and
(J_ab-log(a)W_ab)/Delta; both are nonnegative.

For the first interval [0,y_min], use g(log(Cp^2 y_min/4)) as the
upper constant. Equation (5) gives a one-sided normalized error at most

    E_tip=(15/4)y_min^2,
    E_bulk=(Delta^2/2)(15/4) integral_y_min^1 y sqrt(1-y)dy. (7)

Indeed int_0^b y log(b/y)dy=b^2/4 and sqrt(1-y)<=1. Therefore the
unrounded quadrature Q satisfies Q-E_tip-E_bulk<=kappa<=Q+E_bulk
for EVERY p>0, with p=0 handled exactly. These are uniform analytic
error bounds, not estimates inferred from a frequency sample. Store all
node values and moments, and add outward interval arithmetic error.

At the retained parameters, require -F_eff>50 and a final coefficient
width <.002 at the declared fixtures. The midpoint error bound from (7)
then gives <2e-5 relative error for the positive principal multiplier,
uniformly in p. This is a model estimate, not Z0 for the full metric map.
For any nonnegative self-adjoint scalar Laplacian L0, the same relative
quadratic-form bound transfers to K(sqrt(L0)) by its spectral measure.
The zero mode is unchanged. A closed Laplace form on the actual C2,1
double suffices for THIS statement; it supplies no heat expansion and
does not identify this comparison with the physical determinant Hessian.

## 4. Decisive test of a fixed angular truncation

Keep the common fourth contact in (1), but retain only n<=N in its
positive remainder; call the result K_N. For each fixed N,

    0<K_N(p)<=C p^4/(120pi) kappa_N,infinity,
    kappa_N,infinity=-F_eff+2(H_N+H_(N+q)-H_q),
    lim_(p->infinity) 120pi K_N/(C p^4)=kappa_N,infinity. (8)

This uses 0<=R(x)<=x^2/2 and the fixed finite sum; no tail theorem is
needed. In contrast (2), (4) and the digamma asymptotic give

    kappa(p)=4log[p/(mu0 ell_o)]-62/15+o(1), ell_o=B/sqrt(C). (9)

The constant follows from the normalized weight mean
E(log y)=2log2-31/15. The bound (5) controls the small-y endpoint
when taking the limit. Consequently

    lim_(p->infinity) (K-K_N)/K=1,
    sup_p (K-K_N)/K=1, for every finite N.                (10)

Thus there is NO uniform angular-tail bound below 1/4 relative to this
uncut principal form for a fixed N. Moreover K/K_N is unbounded: using
the inverse of the truncated fourth-order model as a same-space
preconditioner cannot repair the missing logarithm. Compact oscillatory
packets in any smooth interior carry this principal-order obstruction
to the actual uncut TT operator. Local Weyl/fourth contacts can change a
p^4 coefficient on that compact support, not restore the missing log p.
Lower-order local/mixed and boundary terms cannot cancel it either.

This is an ACTUAL obstruction to that uniform approximation, stronger
than an oversized upper bound. It is not a nonexistence or instability
theorem for wormholes. No physical instability is inferred at arbitrarily
high frequencies; a restricted effective-theory band would itself need
a justified regulator, physical matching and error budget. None is added.
Direct residual bounds on a specified smooth forward/adjoint approximation
can still use finite modes if their omitted contributions are certified;
(10) does not forbid function-specific convergence or a finite Galerkin
calculation accompanied by a correct infinite-dimensional tail treatment.

Quantify the obstruction with the retained all-frequency lower bound
kappa>=4log[p/(mu0 ell_o)]-152/15. For any 0<eta<1, set

    p_eta=mu0 ell_o exp[38/15+kappa_N,infinity/(4(1-eta))]. (11)

At p>p_eta the omitted principal fraction exceeds eta. Execute eta=1/4
at N=2,10^6,10^16, using p=2p_eta and the all-angular evaluator. These
are thresholds for the principal model, not already computed thresholds
for the full Hessian with its lower-order classical terms.

For an independent normalization check, the exact low-mode remainder is

    integral R(a x)d alpha
      =2[sqrt((4+x)/x)atanh(sqrt(x/(4+x)))-1]-x/6+x^2/60.

Its omitted positive tail obeys

    0<=K-K_N<=p^6/(840pi) S4(N),
    S4(N)=(2C^2/q)[zeta(2,N+1)-zeta(2,N+q+1)].         (12)

Use R(x)<=x^3/3, int a^3=1/140 and the telescoping identity
d_n/M_n^4=(2C^2/q)[n^-2-(n+q)^-2]. This is a bound at a specified
p or on a stronger regularity class, and is consistent with (10).

## 5. Boundary priority and research decision

The physical gravitational adjoint has both bolt delta sources and lies
in the action form space. The scalar rho^2 log rho term is not H3 at the
bolt; simply assuming a uniform H3 bound on the entire coupled pair
would be invalid. Keep its weak source and residual, as in
[the bolt-adjoint domain](current-bolt-adjoint-domain.md). The auxiliary
thermal-action adjoint is a different object. Its ell*Y/ell*Z cut terms
remain necessary. The accepted <2.92582e-11 per-state cap remains valid
only in its stated mass range and does not become an angular tail bound.

The actual outer mass cubic is reflected C2,1 and its log remainder is
one-sided. [Grubb, Theorem 2.3](https://web.math.ku.dk/~grubb/CPDE14.pdf)
requires symbol regularity and frequency derivatives; it does not provide
the missing symbol seminorms or X inverse. The smooth elliptic boundary
calculus of [Gaarde--Grubb](https://arxiv.org/pdf/math/0703878) cannot
by itself identify this reflected determinant Hessian. The source review
preserves the finite-kernel/contact arguments as additional hypotheses.

Decision: continue the physical weak primal-dual route with the complete
angular logarithmic principal operator retained analytically. Park a
fixed-angular-cutoff inverse with only local fourth restoration. This
replaces one infeasible uniform approximation without changing the theory,
normalization, physical target, or the closed audit. The quadrature closes
the all-frequency MODEL evaluation requirement. It does not yet close the
joint curved/bolt/boundary remainder, actual physical inverse/admission,
nonlinear compatibility, matching, semiclassical validity or stability.

The single decisive physical calculation remains the actual forward/bolt-
adjoint residual certificate, now using the complete logarithmic model in
its approximation and tail treatment. Its next unresolved operation is
the common-subtracted curved and boundary remainder on that trial pair,
including original Weyl/vacuum/a0/a1/a2, York return and both point sources.
The target is |Q_app|+||d_b||_V*||d||_V*<1.667e-9, exterior return
<6.67e-10 and every X/admission and nonlinear inequality. This physical
test is still PROPOSED, not executed by the principal-operator test.

## 6. Frozen execution contract

Use the unchanged finite-seed parameter box, python-flint 0.9.0 at 256
bits and mpmath at 60 digits as an independent diagnostic. Before execution
freeze this argument, checker and source review. Preserve every attempt,
source snapshot, raw record and failure. Check the moment/normalization
algebra, full digamma sum across root coalescence, coefficient width <.002,
uniform principal relative bound <2e-5, and the three eta=1/4 witnesses.
Compare independent real integrals with the enclosures and finite modes
plus their analytic tail at p=1,20. No Lean or full physical solve is selected
for this package. Passing checks support these arguments; they do not
substitute for a domain or continuum solution proof.

Remaining to-do list: common curved/bolt/boundary restoration and remainder;
actual coupled correction with certified dual residual and X/admission;
nonlinear closure and exact compatibility; physical matching, causal
realization, semiclassical validity/noise, all-sector stability, charge/return,
opening, transport, preparation, separate creation and stable wormholes.
