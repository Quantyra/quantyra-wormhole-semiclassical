# Exact compact shape response and coupled interior ellipticity

2026-10-17 HST. E003/S139, supporting S074/S075. **PROGRESS.**

The complete stationary action has a local parametrix for its coupled
conformal and shape Hessian on every fixed smooth compact radial interior.
The shape variable requires one logarithmic factor beyond four derivatives;
the conformal variable requires four derivatives. We derive the exact heat
Hessian, including its contact vertex, before identifying this operator.
The resulting semibounded energy estimate permits an L2 remainder. It is
neither a global boundary inverse nor a physical stability theorem.

This continues the [weighted conformal theorem](stationary-weighted-conformal-range.md)
and [fast stationary phase](stationary-fast-phase-variation.md), without
changing their equations, subtraction scheme, fixed temperature or gates.
The [unweighted obstruction](stationary-coupled-asymptotic-obstruction.md)
and the retained candidate equation defects keep their scopes.

## 1. Actual geometry, domain and exact compact variation

Use the accepted fixed-clock conformal/shape coordinates, after Euclidean
continuation of the prepared time:

```
g = exp(2 sigma)[f d tau^2 + du^2/f + u^2 dOmega^2],
beta = 2pi/kappa, f(u_h)=0, f'(u_h)=2kappa,
Omega = sigma + log(u/B),
g = exp(2 Omega) ghat,
ghat = (B^2/u^2)[f d tau^2 + du^2/f] + B^2 dOmega^2.       (1)
```

B is a fixed reference radius. The product orbit has a regular disk bolt
and the actual outer mixed Dirac domain. The angular operator includes
the fixed magnetic connection. Its full heat trace, including all angular
multiplicities, is denoted Theta_q(t/B^2). The positive square is
Q4=Q2(f)+D_S,q^2/B^2, with Q2=D2^dagger D2 on the complementary domains
already selected in the [domain construction](stationary-exterior-conformal-domain.md).
The [paired phase](stationary-paired-dirac-phase.md) remains treated as there.
First-order D2 need not be selfadjoint for this positive-square argument.

Work in the accepted no-kernel neighborhood at fixed positive temperature.
Fix a compact interval I strictly between bolt and outer end, with f>0.
Take smooth shape directions v_i=delta_i f supported in I, and keep sigma
fixed for this step. These variations leave period, magnetic flux, bolt,
outer geometry and the operator domains unchanged. The orbit density
B^2/u^2 du d tau is independent of f. Identify orthonormal spinors by the
diagonal frame, so adjoints use one Hilbert space and one domain.

Write b_i=v_i/f, and distinguish this b_i from magnetic charge q. With
constant Clifford matrices (and the same harmless overall Dirac phase),

```
A = u sqrt(f)/B, H = u/(B sqrt(f)), c = f'/(4f)-1/(2u),
D2 = gamma_r A(partial_u+c) + gamma_tau H partial_tau,
D_i = gamma_r A[(b_i/2)(partial_u+c)+b_i'/4]
      - gamma_tau H(b_i/2)partial_tau.                    (2)
```

For the AFFINE two-parameter family f+t1 v1+t2 v2 the mixed vertex is

```
D_12 = gamma_r A[-b1 b2/4 (partial_u+c)-(b1 b2)'/8]
       + gamma_tau H(3 b1 b2/4)partial_tau,
Q_i = D_i^dagger D2 + D2^dagger D_i,
Q_12 = D_12^dagger D2 + D2^dagger D_12
       + D_1^dagger D_2 + D_2^dagger D_1.                 (3)
```

In the last line D_1,D_2 are the two variation vertices, not the unvaried
operator denoted D2. Equation (3) follows by differentiating the square;
delta_12 A=-A b1 b2/4, delta_12 H=3H b1 b2/4 and
delta_12 c=-(b1 b2)'/4. In particular the second variation contains a
contact vertex even though f is affine. A pure two-resolvent bubble
would omit part of the physical response.

For trace-class heat operators, cyclicity and Duhamel's formula give

```
K_12(t) = -t Tr_2(Q_12 exp(-t Q2))
 + t^2 integral_0^1 Tr_2[Q_1 exp(-t(1-a)Q2)
                         Q_2 exp(-t a Q2)] da.           (4)
```

There is one simplex integral, with the displayed t^2 coefficient;
cyclicity combines the two time orderings. The exact modulus Hessian is

```
delta_12 Gamma_mod[ghat]
 = (1/2) partial_s|0 { mu0^(2s)/Gamma(s)
     integral_0^infinity t^(s-1) Theta_q(t/B^2) K_12(t) dt }.
                                                               (5)
Gamma_mod = -(1/2) log det Q4 = (1/2) zeta'_Q4(0).
```

The braces mean common four-dimensional Mellin continuation, with the
retained scale and local subtractions, not a convergent unrenormalized
integral at s=0. Sum the complete angular heat trace before continuation.
There is no separately renormalized mode Hessian or lowest-level truncation.
Smooth compact coefficient perturbations have a fixed domain, relatively
bounded operator derivatives and differentiable positive-time semigroups.
Their short-time mixed singularities are removed by differentiating the
same local heat subtractions; their large-time derivatives are integrable
in the fixed kernel-free neighborhood. This proves (5) for these directions.
It does not assert global twice differentiability for moving endpoint data.

For the physical metric the complete identity is
Gamma_mod[g]=Gamma_mod[ghat]+W_E[Omega;ghat], with the full retained bulk
and boundary Weyl transgression. Add the selected local gravitational and
magnetic terms. Their variations, and all contact terms, remain included.
Equivalently one can calculate the interior symbol directly at g; one must
not then add a second copy of the same Weyl contribution.

## 2. Interior symbol and its normalization

The source input is [Moller, arXiv:0705.3857v3](https://arxiv.org/html/0705.3857v3),
Theorem 1 and Lemma 6.1: the rank-four Dirac zeta Hessian has a local
order-four logarithmic symbol. The stated theorem assumes a closed smooth
spin manifold. We use its local symbol, not its global gauge factorization
or sphere stability result. [Okikiolu, arXiv:math/0211100](https://arxiv.org/abs/math/0211100)
separates a locally determined symbol from a possibly globally dependent
order-two contribution. Such lower contributions are retained here.

Here is the interior transfer argument for the actual boundary problem.
Choose nested compact cutoffs around I. The short-time interior parametrix
of the strongly elliptic square has the same near-diagonal coefficients
as an extension to a closed smooth manifold. In the two-vertex heat
integral, the simultaneous short-time region supplies the fourth-order
logarithm. Regions with one heat time bounded away from zero, and the
contact term, can retain global Green coefficients multiplying derivatives
of the variations, of order at most two. Off-diagonal separated supports
give smooth kernels. Boundary reflections cannot change the simultaneous
interior singularity; the fixed magnetic connection is lower order in the
Dirac symbol. Thus localization transfers the leading symbol while
retaining lower state/domain terms. We do NOT declare the entire boundary
or state-dependent remainder smoothing. This argument needs smoothness
on the fixed compact interior, not uniform derivative bounds at an end.

For clarity reproduce the scalar specialization of the source formula.
Let p=|xi|_g be proper cotangent frequency and Pi the projection onto
xi-perpendicular covectors. For the metric tangent K=g^-1 delta g set
T2=tr[(K Pi)^2], T1=tr(K Pi). In four dimensions the principal quadratic
symbol before differentiation is

```
W(s) = p^(4-2s)/(16pi^2) [(2s-3) T2 + T1^2],
Hess zeta'(0) = [W'(0)+eta2 W(0)]/60,
eta2 = 2 sum_{j=1}^5 1/j - sum_{j=1}^2 1/j = 46/15.     (6)
```

Take radial xi and physical tangents sigma variation s and shape b=v/f.
In the ordered (tau, radial, sphere, sphere) orthonormal frame,

```
K = diag(2s+b, 2s-b, 2s, 2s),
T2 = 12s^2+4s b+b^2, T1 = 6s+b,
3T2-T1^2 = 2b^2.                                      (7)
```

Equations (5)-(7), including the fermionic factor 1/2, give the principal
quadratic form per four-component Dirac field per exterior:

```
p^4/(960pi^2) [12s^2+4s b+b^2+2b^2(log p-eta2/2)].       (8)
```

Make log p dimensionless with the retained reference scale. Changing that
scale changes an ordinary order-four shape term, not the logarithmic
coefficient. In particular Phi K=Pi K Pi-tr(K Pi)Pi/3 has
K:Phi K=2b^2/3, and the shape logarithm is p^4 log p b^2/(480pi^2).
The sphere/period integrations multiply the radial density by positive
factors; they do not change this per-volume normalization. Fixed twisting
does not change the ultraviolet rank-four coefficient.

There is an independent sign check from the retained complete trace.
With epsilon=G/(60pi), Euclidean conventions give

```
delta_s I_E = -integral E_trace s dvol/(8pi G),
delta_s E_trace = -6epsilon Box^2 s + lower derivatives,
delta R principal = p^2(6s+b).
```

The density/contact corrections off shell are of order zero. Consequently
H_ss=6epsilon p^4/(8pi G)=p^4/(80pi^2) and
H_sb=epsilon p^4/(8pi G)=p^4/(480pi^2), exactly as in (8).
Here H_sb is the matrix entry, whose cross quadratic term is 2H_sb s b.
These signs describe a Euclidean action Hessian, not a Lorentzian energy.

The retained trace coefficient fixes the scalar-curvature-squared scheme.
Permitted retained local Weyl terms or a change of log reference can add
an ordinary order-four b^2 coefficient. They cannot cancel (8)'s shape
logarithm. An arbitrary change of the R^2 coefficient would change the
trace problem and is outside this theorem. Einstein terms are order two,
Maxwell terms order zero. A nonlinear coordinate change in metric space
adds a first-variation contact with no derivatives of s,b; it does not
alter the principal calculation at a comparison with nonzero defect.

## 3. Complete coupled interior theorem

Set theta=s+b/6. Completing the square in (8) gives

```
p^4 theta^2/(80pi^2)
 + p^4 b^2 [log p + 1/3-eta2/2]/(480pi^2).               (9)
```

Let H_I denote the localized complete Euclidean action Hessian in these
two radial variables. The output covectors transform with the transpose
of s=theta-b/6: F_theta=F_s, F_b,new=F_b,old-F_s/6. On I, the original
Euler gradients are invertible smooth density multiples of trace and
pressure anisotropy: F_s=-E_trace/(8pi G),
F_b,old=(E_tau-E_r)/(16pi G), per physical volume. Division by f on I
is harmless. Thus the variational system controls both local equations.
It does not impose the separate radial integration constraint.

After (9), the symbol has positive diagonal leading entries

```
a0 p^4, d0 p^4 log p,  a0=1/(80pi^2), d0=1/(480pi^2),
```

an ordinary order-four term only in the shape diagonal, and remainders
of order at most three with one logarithm. Lower global contact terms
are included in these remainders. Smooth density conjugation affects
lower symbols only. Constants below may depend on I, the comparison,
the fixed state/domain and retained coefficients.

Define H^{t,log} locally by Fourier multiplier
<p>^t log(2+|p|). The full local bounded map is

```
H_I: H^{t+4}_theta x H^{t+4,log}_b -> H^t x H^t.          (10)
```

At sufficiently large p include the entire ordinary b p^4 coefficient
in the diagonal symbol and invert it. The reciprocal entries have sizes
p^-4 and p^-4/log p; their derivatives satisfy the corresponding symbol
bounds. Composition with an order-three-log remainder loses one power
of p, up to logarithms. Iterating this error in the local symbol calculus
therefore constructs a two-sided parametrix modulo a smooth kernel,
between nested interior cutoffs. This proves mixed elliptic regularity
and the usual local estimate with a lower norm term. It is a parametrix
for a fixed smooth comparison; it asserts no global boundary condition,
nonlinear neighborhood estimate or quantitative inverse norm.

The same operator has the coupled energy bound, for compact smooth
directions supported in I,

```
H_I[(theta,b),(theta,b)]
 >= c (||theta||_H2^2 + ||b||_H^{2,sqrt(log)}^2)
    - C ||(theta,b)||_L2^2,        c>0, C finite.          (11)
```

To prove it, quantize square roots of the two positive principal entries.
Their adjoint products differ from the diagonal operator by order-three
logarithmic terms. The ordinary b p^4 coefficient is absorbed at large
p since 1/log p tends to zero. For a cross remainder one can use
p^3 log p |theta b| <= delta p^4 |theta|^2
+ C_delta p^2 log^2 p |b|^2; the latter is smaller than the shape energy
at large p. Diagonal lower terms obey the same interpolation bound.
Low frequencies and compact-cutoff commutators contribute at most an
arbitrarily small principal energy plus C_delta L2 norm. Absorb the small
terms to obtain (11). No smallness of the background equation defect is
used. Neither c nor C is uniform as temperature goes to zero or over the
preparation interval. Negative low-frequency directions remain possible.

## 4. Consequence for the selected global solver route

An ordinary four-derivative shape domain with no logarithmic factor is
insufficient. On a smaller interval choose a nonzero cutoff chi and
b_n=n^(-t-4) chi exp(i n y), in smooth proper radial coordinate y.
Its H^{t+4} norm stays bounded, but the shape response in H^t grows as
log n, by the nonzero principal coefficient. Lower symbols cannot
cancel that growth. The same dyadic construction with amplitude
2^(-(4+eta)j) proves failure of C^{4,eta} -> C^{0,eta}, 0<eta<1.
Endpoint decay weights do not help compactly supported packets.

A compatible local Holder scale uses C^{4,eta} for theta and
B^{4+eta,log}_{infinity,infinity} for b, whose dyadic norm is
sup_j (1+j)2^((4+eta)j)||Delta_j b||_infinity. A positive extra derivative
also bounds the logarithm, consistent with earlier C5-to-C0 estimates,
but using an unnecessarily stronger conformal domain does not establish
a same-scale Fredholm theorem. Equation (10) identifies the sharp local
mixed scale for the next global construction.

The accepted conformal inverse has local order -4, while its cross
variations have order four without a logarithm, by the local conformal
Ward identity. Accordingly the conformal Schur correction has ordinary
order four and cannot cancel the shape p^4 log p term. The square
completion exhibits this directly. Global projections, kernel directions,
moving-bolt data and endpoint constraints still require their own maps.

In particular (10)-(11) do not settle the normal operator at infinity,
the allowed physical r^-4 log tails, the actual conformal matching matrix,
or the shape response to a moving bolt. Trace and anisotropy still leave
the C/r^4 residual unless one radial constraint is imposed. The earlier
nonzero shape residual at the fast AdS trace root remains relevant.
Interior semiboundedness allows all of those global obstructions; it
does not supersede the unweighted escaping-sequence result or prove an
exact RN exterior, stability, nonexistence, physical charge matching or
any S074-S079 gate. The milestone is the exact compact response and a
complete interior inverse construction in the required mixed scale.

## 5. Focused verification contract

The accompanying checker is limited to exact finite algebra: product
density and Dirac vertices, the mixed contact, a noncommuting square
expansion, scalar Duhamel normalization, principal projection and zeta
coefficients, Euclidean trace signs, square completion and reciprocal
symbols. All symbolic residuals must vanish exactly; rational inequalities
must hold. There are no numerical tolerances or fitted spectral data.
Heat localization, mixed pseudodifferential estimates, semiboundedness,
packet asymptotics and all global domain issues remain ordinary arguments.
Freeze the reviewed note and checker before execution; append all runs,
retain failed input bytes, and hash raw UTF8/LF bytes without normalization.

Remaining to-do list: establish the global weighted shape/constraint map,
bolt and outer normal problems, physical logarithmic asymptotics and actual
conformal matching matrix; prove the projected coupled inverse and choose
the physical branch, or justify controlled spatial reduction; establish
exterior/state and charge matching, transfer S138, and complete driven
stress, radius, controller, capture, tail, spatial area and S074-S080.
