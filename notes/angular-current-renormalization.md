# Angular current running and the limit of a contact replacement

2026-10-11, S051/S053. The nonspherical photon interaction requires more
than the previously proved angular pair norm. Its local current-current
approximation generates new angular coupling ratios already at one loop.
Keeping all rotation-invariant couplings gives a closed one-loop flow and
a channel-uniform bound on its size. In the declared weak-coupling family,
polynomial scale separations do not make that flow strong.

This does **not** close the physical interaction error. The contact
Hamiltonian has no cutoff-uniform relative form bound on the free vacuum,
and the lowest throat photon frequency remains comparable to the clock
frequency. The next physical response calculation must retain the photon
propagators or supply an equivalent renormalized matching argument. The
Gaussian [singlet decay result](singlet-static-curvature.md) is preserved.
The primary remains physical long-return damping, with finite-duration
energy transfer as fallback; no separate current-algebra project is selected.

## 1. Local interaction and normalization

Take `q>=2`; at `q=1` the interaction considered here is absent.
Let `T_la=P_q Y_la P_q` be the real Hermitian sphere-LLL vertices from
the [angular norm proof](weak-coupling-window.md), `1<=l<q`,
`a=1,...,2l+1`. With `beta_l=(q-1)!q!/[(q-1-l)!(q+l)!]`,

`Tr(T_la T_mb)=q beta_l delta_lm delta_ab/(4pi)`,
`sum_la T_la^2/[l(l+1)]=c_q I_q`,
`c_q=(H_q-1)/(4pi)`.                                    (1)

The identity harmonic is excluded because its Schwinger dynamics was
retained separately. Normalize `E_la=T_la/sqrt(q beta_l/(4pi))`.
These form a real orthonormal basis of traceless Hermitian matrices:
`Tr(E_A E_B)=delta_AB`, `A=(l,a)`.

In a locally flat massless current theory use

`J_X(z) J_Y(0) ~ Tr(XY)/z^2 + J_[X,Y](0)/z`,              (2)

and an independent antiholomorphic copy. In physical fermion-density
normalization `J_X=2pi :psi_R^dagger X psi_R:`. Define the local
Euclidean perturbation by

`S_int=(1/(2pi)) int d^2z sum_AB G_AB J_EA(z) Jbar_EB(zbar)`.

The zero-longitudinal-momentum exchange of a nonspherical vector has

`G_A B(0)=delta_AB sigma_l`,
`sigma_l=lambda q beta_l/[4pi l(l+1)]`, `lambda=g^2/pi`.   (3)

For this coefficient, the reduced vector action contains
`l(l+1) A_a A^a/(2g^2)` and source `i A_a j^a`. Its algebraic Gaussian
integral gives `g^2 j_a j^a/[2l(l+1)]`. With
`j_tau=rho_R+rho_L`, `j_s=i(rho_R-rho_L)`, this is
`2g^2 rho_R rho_L/[l(l+1)]`; inserting `J=2pi rho` gives (3).
This calculation neglects the vector kinetic term. It is a local contact
matching coefficient, not the full Lorentzian exchange kernel on the
wormhole. Curvature, joins, return boundary conditions and higher Landau
levels are not erased by writing (3).

## 2. One-loop closure needs all angular couplings

Put `[E_A,E_B]=i sum_C f_ABC E_C`, with real structure constants.
For logarithmic length scale `t`, the one-loop beta function in this
normalization is

`dG_CD/dt=(1/2) sum_ABEF G_AB G_EF f_AEC f_BFD`.           (4)

This is a statement about the one-loop truncated differential equation.
To derive it, multiply the two simple-pole terms in (2). The coefficient
of `J_C Jbar_D/|z|^2` is the negative of the sum in (4), because of the
two factors of `i`. Integrating a radial shell with measure `d^2z/(2pi)`
gives `dt`; the second-order exponential expansion contributes `1/2`.
Re-exponentiation therefore gives the positive sign in (4).

Under the spin `(q-1)/2` rotation representation, traceless matrices
decompose into spins `l=1,...,q-1`, each once. Equation (4) is equivariant
under rotations. Uniqueness of its finite-dimensional ODE consequently
preserves the rotation-invariant family: each angular block remains a
scalar `sigma_l(t)` and different blocks do not mix. Commutators are
traceless, so no identity-current coupling is generated at this order.
There are `q-1` running coefficients, not `q^2-1` independent bosons.

The fixed ratios in (3) are generally **not** preserved. An explicit
counterexample uses `q=3`. Write the spin-one and spin-two coefficients
as `a,b`. Direct commutators give

`a'=(a^2+5b^2)/2`, `b'=3ab`.                             (5)

The tree coefficients obey `a=15b>0`. At that point
`a'/a=23a/45`, whereas `b'/b=3a`. Thus running a single coefficient
`lambda(t)` multiplying the original angular kernel misses a generated
operator. This does not mean that four-dimensional QED acquires arbitrary
new parameters: the additional effective coefficients are determined by
matching and running. It means that the one-parameter contact ansatz is
not closed. For `q=2` there is just one angular block and `sigma'=2sigma^2`,
consistent with the usual `Tr(E_A E_B)=delta_AB` normalization.

## 3. A bound on the complete positive one-loop flow

For a real symmetric positive semidefinite `G`, factor
`G_AB=sum_p u_pA u_pB` and set `X_p=sum_A u_pA E_A`. Define

`Q(G)=sum_AB G_AB E_A E_B=sum_p X_p^2`, `kappa=||Q(G)||`.

The beta function in (4) is a positive Gram matrix of the Hermitian
commutators `C_pr=-i[X_p,X_r]`, with coefficient `1/2`. It is therefore
positive semidefinite. The positive cone is invariant; in forward time
`G(t)-G(0)>=0` wherever the solution exists. Moreover

`dQ/dt=(1/2) sum_pr C_pr^2 <= 2 kappa Q <=2 kappa^2 I`.    (6)

Indeed `(A-B)^dagger(A-B)<=2(A^dagger A+B^dagger B)` applied to
`A=X_p X_r`, `B=X_r X_p`, followed by summation, gives
`sum_pr C_pr^2<=4 sum_p X_p Q X_p<=4 kappa Q`.
The upper right derivative of the largest eigenvalue in (6) gives

`kappa(t)<=kappa_0/(1-2 kappa_0 t)`,
`||Q(G(t)-G(0))||<=2 kappa_0^2 t/(1-2 kappa_0 t)`,        (7)

for `0<=t<1/(2kappa_0)`. For the rotation-invariant solution `Q=kappa I`
by irreducibility, and `kappa=Tr(G)/q`. More generally
`Tr G=Tr Q<=q kappa`, so (7) bounds every coefficient and prevents
finite-time escape on the displayed interval. It therefore proves
existence on that interval as well as the estimate.

For positive `G`, the angular pair operator
`K(G)=sum_AB G_AB E_A tensor conjugate(E_B)` has norm at most `kappa`
by row/column Cauchy-Schwarz. In the rotation-invariant case, the
vectorized identity is an eigenvector with eigenvalue `kappa`, so equality
holds. The same argument applies to the positive increment of `G`.
Also `||G||<=Tr G=q kappa`. From (1)-(3),

`kappa_0=lambda c_q=g^2(H_q-1)/(4pi^2)`.

For `g^2=q^-2` and any fixed `A>0`, let `t<=A log q`. Then

`kappa(t)/kappa_0=1+O(q^-2 log^2 q)`,
`||G(t)||=O(q^-1 log q)`,
`||G(t)-G(0)||=O(q^-3 log^3 q)`.                          (8)

Both coefficient norms are smaller than `delta=q^-1/16`. This rules out
one-loop strong running caused solely by any polynomial ratio of scales
within this contact model. It is not an all-loop estimate, a bound on
stress correlation functions, or a proof that this local theory applies
through every one of those scales. The length of the actual longest
inverse-tail return is polynomial, so its logarithm alone supplies no
parametric one-loop obstruction.

## 4. Why a finite angular norm does not bound the filled sea

Consider a diagnostic cylinder of circumference `2pi`, with the free
left/right current vacuum `Omega`. Its currents satisfy

`[J_m^A,J_n^B]=i f_ABC J_(m+n)^C + m delta_AB delta_(m+n,0)`,
`J_n Omega=Jbar_n Omega=0` for `n>=0`,
`H_0 J_-n^A Jbar_-n^B Omega=2n J_-n^A Jbar_-n^B Omega`.

The spatially integrated contact interaction, cut at current mode `N`,
has vacuum column

`V_N Omega=sum_(n=1)^N sum_AB G_AB J_-n^A Jbar_-n^B Omega`.

At fixed `n`, the products have Gram matrix
`n^2 delta_AC delta_BD`; distinct `n` have distinct total energies.
Therefore exactly

`||V_N Omega||^2=||G||_F^2 N(N+1)(2N+1)/6`,
`||(H_0+1)^(-1/2)V_N Omega||^2`
`                  =||G||_F^2 sum_(n=1)^N n^2/(2n+1)`.    (9)

For every nonzero fixed `G`, the second quantity grows as
`||G||_F^2 N^2/4`. Thus no cutoff-independent constant can make
`|<psi,V_N Omega>|<=C ||(H_0+1)^(1/2)psi||` for all finite-energy
`psi`. In particular, a relative form bound on the free Hamiltonian,
uniform under removal of this contact regulator, is impossible. A
vacuum-energy constant does not change these off-vacuum matrix elements.
This is an obstruction to that proposed *bare-operator estimate*, not
an impossibility of defining a renormalized current theory or QED.

For the projected tree coefficients,

`||G(0)||_F^2=lambda^2 q^2/(16pi^2)`
`       *sum_(l=1)^(q-1) (2l+1) beta_l^2/[l^2(l+1)^2]`.

The sum tends to one. It is bounded by one since
`(2l+1)/[l^2(l+1)^2]=1/l^2-1/(l+1)^2`; dominated convergence and
`beta_l->1` at fixed `l` give the limit. Thus this vacuum effect is
not controlled by simply treating a finite number of angular particles.

A finite-regulator statement survives. With a further finite fermion
mode cutoff, the matrices are bounded, and Duhamel gives
`||(exp(-it(H_0+V_N))-exp(-it H_0))Omega||<=|t| ||V_N Omega||`.
For `g^2=q^-2`, `N=floor(sqrt(q))`, and `t=O(delta^-1)`, (9) makes
this upper bound `O(q^-3/16)`. This is a prepared-free-vacuum comparison
in a cutoff contact model. It does not construct the interacting ground
state, justify discarded higher frequencies, cover unequal returns, or
bound the clock stress response. The order of cutoff and weak-coupling
limits matters, and the finite statement cannot remove those premises.

## 5. The lowest throat photons do not become infinitely fast

There is a second reason to retain the propagator. In the fixed-metric
Maxwell sector on the exact constant-radius `AdS2 x S2` throat, the
electric-parity vector harmonic has
covariant mass squared `l(l+1)/r_e^2`. In two dimensions its field
strength dual obeys a scalar equation with that mass. In global throat
coordinates, `N=r_e x sec(sigma)`, it reduces to

`H_l=-partial_sigma^2+l(l+1)sec^2(sigma)`,
`sigma in (-pi/2,pi/2)`, `omega^2=x^2 eigenvalue(H_l)`.

For `a=l+1`, factorization gives
`H_l=(-partial_sigma+a tan sigma)(partial_sigma+a tan sigma)+a^2`.
The normalizable function `cos^a(sigma)` saturates the bound. Hence the
lowest frequency is exactly `(l+1)x`; the full Friedrichs spectrum is
`(n+l+1)x`, `n=0,1,...`, as substitution of
`cos^(l+1)(sigma) C_n^(l+1)(sin sigma)` verifies.

For `l=1`, the lowest frequency is `2x`, while the proposed symmetric
clock approaches frequency `x`. Their ratio is independent of `q` and
`g`. Large flux alone therefore supplies no small frequency ratio for
the lowest nonspherical photon. Small coupling can still make its total
effect small, but does not justify a uniformly accurate instantaneous
replacement mode by mode. The physical finite-cutoff mouth loads and
possible exterior radiation require separate matching; this spectrum is
the exact global-throat comparator, not the full wormhole spectrum.
Gravito-electromagnetic mixing is not included in this comparator. In
particular, MMP Appendix A.1 also identifies the sphere-isometry `SU(2)`
gauge fields, massless at the classical reduction level and with a
fermion-induced mass. Those are an additional sector to match; the
Maxwell gap above is not assigned to them.

## 6. Source review, verification and decision

The [script](../checks/angular_current_renormalization.py) and
[receipt](../evidence/angular-current-renormalization.json) check direct
LLL quadrature, current commutators, one-loop block closure and flow,
the exact `q=3` counterexample, explicit fermionic vacuum matrix elements,
and the throat photon eigenfunctions. Numerical flow checks verify the
finite ODE and its analytic bound; they do not sum higher loops.

[Gerganov, LeClair and Moriconi, hep-th/0011189v3](https://arxiv.org/abs/hep-th/0011189v3)
give the current OPE, anisotropic closure condition and one-loop term,
especially equations (1)-(2), (10)-(11) and the quadratic term of (17).
Their all-order expression is not imported. [Ludwig and Wiese,
cond-mat/0211531v1](https://arxiv.org/abs/cond-mat/0211531v1) report a
four-loop counterexample to the claimed general expression, explicitly
at level zero. That does not itself calculate a level-one wormhole
remainder. Their model and one-loop derivation are checked separately;
the full four-loop calculation is not newly certified here. The
[source audit](../evidence/angular-current-source-review.json) records
the read boundaries, versions, hashes and searches. The later
[bosonization discussion](https://arxiv.org/html/2606.12528v2#S4.SS3.SSS1)
also explicitly leaves additional quartic interactions and full flavor
current/boundary data outside the singlet simplification.

Continue the physical target. The one-loop local running is controlled,
but a cutoff-free many-body response estimate cannot be inferred from
the angular pair norm or the running couplings. The next decisive step is
the renormalized finite-frequency photon/current response with its physical
matching terms. This finding does not invalidate the Gaussian theorem,
prove physical damping, or establish novelty. No Lean, Lake, GCP,
publication or external-message action occurred.

Remaining to-do list: control the full nonspherical photon/current response
and physical background, return, collar/vertex and state matching; finish
novelty and the consequential physical result; independently verify,
GCP-certify suitable endpoints and prepare the reproducible manuscript.
