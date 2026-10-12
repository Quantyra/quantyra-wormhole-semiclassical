# Collective mouth constraints and the static Casimir response

2026-10-11, E003/S053. The electric mouth limit has a finite, explicitly
computable nonlocal vacuum response even with unequal exterior returns.
For a fixed rank r of collective Dirichlet constraints, its first two
throat-length derivatives are bounded independently of the number of
channels and every exterior length. In the retained scaling its static
clock correction is smaller than q delta. This conclusion concerns a
defined Gaussian graph, not its full dynamical response or the mouth
stress of the smooth semiclassical geometry.

There is a separate obstruction to using this limit as a physical static
estimate: in the previous optical subtraction prescription, the local
finite term c integral V has clock curvature of order q at the deeper
cut. The ideal-boundary determinant does not determine that term or its
four-dimensional matching. Norm-resolvent convergence cannot settle it.
The next required comparison must retain the local renormalization terms.

## Defined model and observable

Let q real massless scalar current fields have separate circles of lengths
L_j=T+D_j, T>0 and D_j>0. They share a throat coordinate [0,T]; they are
not identified as one field. Start with the usual channelwise continuity
at both mouths. Fix a real q by r matrix W with W^t W=I_r, 1<=r<=q,
independent of T. Impose W^t f(0)=W^t f(T)=0 in the closed free form.
Denote the free and constrained operators H_0 and H_W. Their domains
retain all neutral matching and all exterior intervals. Zero oscillators
contribute zero; compact charge and holonomy sectors are not selected.

For the electric limit derived in [the Gauss-law note](gauss-mouth-boundary.md),
r=1 and W_j=1/sqrt(q). The rank-r theorem below is mathematical. In
particular, the three noncommuting rotational currents are not three
additional independent bosons on these q circles. The
[rotational action](rotational-current-response.md) also has a dynamical
Maxwell dipole field; it has not been replaced by this boundary model.

Define E_W(T,D) as the difference of oscillator Casimir energies of H_W
and H_0 with scale-invariant vertex self energies normalized to zero.
Equivalently use the explicit convergent expression (2) below. Hold D
and W fixed when differentiating T. This is a global energy derivative,
not a pointwise renormalized stress or a physical observer comparison.

## Boundary determinant with every return retained

For k>0 put x=exp(-kT), y_j=exp(-kD_j) and

```text
C_+(k)=W^t diag[(1+y_j)/(1-x y_j)] W,
C_-(k)=W^t diag[(1-y_j)/(1-x y_j)] W,
F(k)=(1-x^2)^r det C_+(k) det C_-(k).                    (1)
E_W(T,D)=(1/(2pi)) integral_0^infinity log F(k) dk.       (2)
```

All matrices in (1) are positive. The free circle Green function gives
the endpoint covariance G=B(H_0+k^2)^(-1)B^t, where
Bf=(W^t f(0),W^t f(T)). Its diagonal and off-diagonal blocks are
W^t diag[(1+x y_j)/(1-x y_j)]W/(2k) and
W^t diag[(x+y_j)/(1-x y_j)]W/(2k). Passing to endpoint sum/difference
blocks gives 2k G_+=(1+x)C_+ and 2k G_-=(1-x)C_-.
Thus F=det(2kG). This uses a compression before the determinant;
replacing unequal returns by their mean is not the same operation.

An independent scattering construction fixes the normalization and sign.
At either mouth write P=W W^t, Q=I-P. Incoming and outgoing amplitudes
on the throat/exterior edges obey the orthogonal scattering matrix

```text
S_W=[[-P,Q],[Q,-P]],     L=diag(x I_q, y_1,...,y_q).
d_W(k)=det[I_(2q)-(S_W L)^2],
d_0(k)=product_j (1-x y_j)^2,
F(k)=d_W(k)/d_0(k).                                      (3)
```

Indeed continuity makes the shared endpoint value
c=Q(a_throat+a_exterior). Projected Kirchhoff matching is the natural
condition of the restricted form, giving the displayed S_W. Propagation
to the second identical mouth and back gives (S_W L)^2. For explicit
determinant elimination set U=(W;W), so S_W=S_0-U U^t. The determinant
lemma for I-S_W L gives a reduced factor
I+U^t L(I-S_0 L)^(-1)U=(1+x)C_+. For I+S_W L it gives
I-U^t L(I+S_0 L)^(-1)U=(1-x)C_-. Each free determinant is
product_j(1-x y_j). Multiplying proves (3).
The free energy from d_0 is -sum_j pi/(6L_j), the standard real scalar
circle result. Imaginary-frequency integration of the secular determinant
after removing the common length term gives (2). This is the established
quantum-graph spectral-determinant method, independently specialized here.

There is no generic-length or nondegeneracy requirement: for k>0 all the
matrices are positive/invertible, so equal lengths follow directly.
As k tends to zero, C_+ is k^(-1) times a positive definite matrix and
C_- tends to W^t diag[D_j/(T+D_j)]W. Therefore F tends to a finite
strictly positive value. At infinity log F decays exponentially for
each fixed positive T,D. This proves convergence of (2), including the
zero-frequency endpoint. Constants in that convergence statement need
not be uniform in D_min; the derivative bounds below are uniform.

## A response bound independent of the return distribution

For all q,r,W,T,D as specified,

```text
|partial_T E_W|   <= pi r/(6 T^2),
|partial_T^2 E_W| <= 2pi r/(3 T^3).                       (4)
```

Here is a direct proof. Let ell_+=log det C_+, ell_-=log det C_-.
Writing c_j=(1+/-y_j)/(1-x y_j), differentiation yields

```text
-k x/(1-x) <= partial_T log c_j <= 0,
partial_T^2 log c_j = k^2 x y_j/(1-x y_j)^2 >= 0.
```

Cauchy-Binet expresses det(W^t diag(c) W) as
sum_{|A|=r} det(W_A)^2 product_{j in A} c_j. These nonnegative
weights are not all zero. A log-sum of exponentials of convex functions
is convex (differentiate: mean second derivative plus a nonnegative
variance). Hence ell_+'' and ell_-'' are nonnegative. Also

```text
-r k x/(1-x) <= ell_+/-' <= 0,
0 <= ell_+/-''
  =Tr(C^-1 C''-C^-1 C' C^-1 C')
  <=Tr(C^-1 C'')
  <=r k^2 x(1+x)/(1-x)^2 <=2r k^2 x/(1-x)^2.             (5)
```

The penultimate inequality follows coefficientwise from
c_j''/c_j=k^2 z_j(1+z_j)/(1-z_j)^2 with z_j=x y_j<=x.
The subtracted trace is nonnegative since C^-1/2 C' C^-1/2 is symmetric.
The derivative and second derivative of r log(1-x^2) are respectively
2r k x^2/(1-x^2) and -4r k^2 x^2/(1-x^2)^2. Combining their signs
with (5), rather than adding both upper magnitudes, gives

```text
|partial_T log F| <=2r k x/(1-x),
|partial_T^2 log F| <=4r k^2 x/(1-x)^2.                   (6)
```

Both are integrable at zero and infinity. Their integrals use
integral k/(exp(kT)-1) dk=pi^2/(6T^2) and
integral k^2 exp(-kT)/(1-exp(-kT))^2 dk=pi^2/(3T^3).
Dominated differentiation on compact positive T intervals proves (4).
No maximum return, minimum free gap, or large-q approximation enters.

The exterior-only limit, which can itself be large, is

```text
Y_W(k)=W^t diag(y_j) W,
E_ext(D,W)=(1/(2pi)) integral log det(I_r-Y_W(k)^2) dk,
-pi r/(24 D_min) <= E_ext <=0,
|E_W(T,D)-E_ext(D,W)| <=pi r/(6T).                        (7)
```

It follows by taking T to infinity in (1) and integrating (4).
The eigenvalues of Y_W lie in [0,exp(-k D_min)]. Direct dominated
convergence for fixed D identifies the limit; (4) controls the integral
from finite T to infinity. Do not discard E_ext from a variation of the
exterior or the complete geometry. It is independent only of this T.

For the prescribed relative clock T(phi)=T_0 exp(-phi),

```text
|partial_phi E_W| <=pi r/(6T),
|partial_phi^2 E_W| <=5pi r/(6T).                         (8)
```

For fixed r, T bounded away from zero and delta=log^(-2)q, both
derivatives are o(q delta). This closes the static *nonlocal boundary*
comparison for the electric r=1 limit. It does not close the finite
singular-potential or four-dimensional static premise.

## Reproduced equal-return baseline

If every D_j=D, rotate the fields by a fixed orthogonal matrix. H_W
has q-r free circles and r pairs of Dirichlet intervals of lengths T,D.
Using the established interval energy -pi/(24L) gives exactly

```text
E_W=r[-pi/(24T)-pi/(24D)+pi/(6(T+D))].                    (9)
```

Equations (1)-(2) reproduce this by integrating
r log[(1-exp(-2kT))(1-exp(-2kD))/(1-exp(-k(T+D)))^2].
This checks the sign, the real-scalar normalization and both length
derivatives. The result is not equivalent to simply removing r circles:
the charged throat intervals remain, and unequal returns mix currents
outside the throat. These baseline energies and the determinant method
are prior work.

## Why the finite mouth still needs renormalization matching

The earlier [static optical prescription](singlet-static-curvature.md)
contains the explicitly fixed local term c I, I=integral V(s) ds. For
the deeper finite-mouth family V=a sec^2(s), a=1/(4pi^2 q), R=q^2,
the [Gauss-law calculation](gauss-mouth-boundary.md) gives

```text
I=2aR=q/(2pi^2).                                         (10)
```

For the clock H(phi)=exp(phi) H_throat+H_exterior, put the throat into
unit-speed optical coordinate u=exp(-phi)s. The quadratic potential is
V_phi(u)=exp(2phi)V(exp(phi)u); its optical integral is I_phi=exp(phi)I.
Two optical prescriptions differing only by a finite coefficient Delta c
therefore differ by

```text
Delta E(phi)=Delta c exp(phi) I,
partial_phi Delta E(0)=partial_phi^2 Delta E(0)
                       =Delta c q/(2pi^2).               (11)
```

For fixed nonzero Delta c the ratio to q delta grows as log^2 q. This
is a precise failure of a scheme-independent inference from (4) or from
norm-resolvent convergence. Changing scheme in a physical theory also
changes the matched renormalized gravitational/gauge coefficients; it
does not change the physics or authorize tuning c after inspection.
No assertion that (11) survives uncancelled in the full physical energy
is made. Instead, that cancellation/matching must be calculated before
claiming an o(q delta) *total* static error. A fixed c alone does not
bound the remaining finite-profile determinant either.

The potential in (10) reaches the movable endpoint and changes with the
clock. It does not satisfy a fixed compact-interior-potential hypothesis
that would make a graph-length force automatically free of local terms.
Harrison and Kirsten explicitly distinguish these cases in
[1201.3754v1, sections 4-5](https://arxiv.org/pdf/1201.3754v1).
The ideal-boundary finite energy cannot replace this physical calculation.

## Prior work and decision

[Harrison and Kirsten, 0911.2509v3, sections 7.3-7.4](https://arxiv.org/pdf/0911.2509v3)
derive graph forces from vertex matching and ultraviolet asymptotics.
Their [2012 calculation](https://arxiv.org/pdf/1201.3754v1) includes
potentials and identifies the extra renormalization dependence when the
compact-support restriction is relaxed. [MMP Appendix A](https://arxiv.org/html/1807.04726v3#A1)
already supplies the electric and rotational screening mechanism. The
present bound is a specialized comparison inside the research model;
it is not an accepted independent novel physics target.

Continue long-return damping as the provisional primary and matched
finite-duration transfer as fallback, conditional on S051's still-open
novelty review. Narrow the immediate static work to the finite-mouth
determinant and local four-dimensional coefficient matching. Do not add
angular graph constraints by counting currents or commission formalization
of this supporting lemma as the research endpoint.

[The checker](../checks/collective_boundary_casimir.py) tests independent
vertex scattering, equal-return spectra, numerical differentiation,
the distribution-independent inequalities and the clock counterterm.
[Its receipt](../evidence/collective-boundary-casimir.json) and
[source review](../evidence/collective-boundary-casimir-source-review.json)
pin the inputs, environment and reading scope. Floating-point checks
supplement the ordinary proof. No physical stability/shortcut or novel
result has been accepted; no Lean/Lake invocation occurred.

Remaining to-do list: match the finite-mouth determinant and local
renormalized coefficients; construct the coupled angular/compact state
with the physical exterior; control smooth matching and finite-disturbance
transport; complete novelty and independent verification, certify suitable
endpoints on GCP, and prepare the reproducible manuscript.
