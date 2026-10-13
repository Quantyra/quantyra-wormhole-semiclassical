# Phase action, quantum noise and recoil at the boundary

S072, 2026-10-13. **The projected RN current sector determines an
off-shell quadratic boundary-phase action and a positive oscillator
representation of its fluctuations.** An exact quantum-current identity
then retains the bath backreaction and mixed noise terms during a
transfer interaction. This supplies terms missing from the earlier
mean-response calculation. The electric
oscillator contribution to the prepared midpoint correlation is small
in the selected finite-band comparison. However, a sharp one-mode
counterexample proves that a small correlation loss alone does not
bound energy deposited during switching. The physical interacting
fermion stress remains the next calculation, not a completed result.

## Boundary variation fixes the action

Start from [the constrained Gaussian model](causal-compact-response.md)
on the fixed RN exterior. Use tortoise coordinate `x<0`, with AdS
boundary `x=0`, and Euclidean time of period `beta`. The electric
sector has `lambda=N_f q/pi`, `kappa=4pi r^2/g^2` and

```
S_E = (1/2) integral d tau dx [
  (kappa/f)(partial_x a_tau-partial_tau a_x)^2
  + lambda ((a_tau-partial_tau chi)^2
            +(a_x-partial_x chi)^2)].                     (1)
```

This is the positive Euclidean version of the center-of-mass charged
boson in [MMP appendix A](https://arxiv.org/html/1807.04726v3).
Fix the physical Maxwell source at zero and prescribe
`chi(0,tau)=alpha(tau)`. In unitary gauge the boundary value is
`a_tau(0)=-partial_tau alpha`. Retain `a_x` while varying (1).
For a nonzero Matsubara frequency `s=|omega_n|`, define
`p=(kappa/f)(partial_x a_tau-partial_tau a_x)`. The constraints are

```
p_x=lambda a_tau,       p_tau=-lambda a_x,
[-partial_x^2+V+s^2]p=0,    V=f lambda/kappa.              (2)
```

Regularity gives decay toward the Euclidean horizon. Let `K(s)` be the
Dirichlet-to-Neumann coefficient of (2): the decaying solution with
unit boundary value has derivative `K(s)` at zero. Then
`p(0)=-i omega_n I(s) alpha_n`, where `I(s)=lambda/K(s)` is precisely
the earlier boundary-transfer response, including its normalization.
Integration by parts in (1) leaves `S_E=(1/2) integral a_tau p` at the
outer boundary. Thus, with `alpha(tau)=sum alpha_n exp(i omega_n tau)`,

```
S_boundary = (beta/2) sum_(n!=0)
             omega_n^2 I(|omega_n|) |alpha_n|^2.           (3)
```

This derives (3) from the boundary variational problem. A static
susceptibility by itself would not have justified it. The Gaussian
determinant with homogeneous boundary data is independent of `alpha`;
it does not add another source-dependent quadratic term.

For an independent check eliminate `a_x` first. The electric action is
`(lambda/2) integral [a_tau,x^2/(V+s^2)+a_tau^2]`. Its minimized
coefficient at `a_tau(0)=1` is `I(s)`. It is the primal variational
problem dual to the one defining `K(s)`.

The axial sector follows by retaining both potentials. Write the
previous constant mass matrix as `M=L L^T`, set `p=L Z` and
`b=L^T a_tau`. Its canonical potential is
`V=f L^T K_radial^(-1)L`. The source obeys
`b(0)=-i omega_n sqrt(lambda_J) e_2 alpha_n`, since
`L^T n=sqrt(lambda_J)e_2`. Elimination of the radial potentials gives

```
S_E = (1/2) integral [b_x^T(V+s^2 1)^(-1)b_x+b^T b].      (4)
```

The dual problem fixes `Z_2(0)=1` and leaves `Z_1(0)` free, producing
its natural zero-derivative boundary condition. Equations (3)--(4)
therefore apply to each linearized angular generator with `I=I_J` and
`lambda=lambda_J`. They do not construct the nonlinear compact SU(2)
action or its representation weights.

The derivation covers fluctuations in the topologically trivial sector.
It does not determine winding-sector weights at the Euclidean horizon,
the black-hole degeneracies, a microscopic compact charge lattice, or
the global nonlinear angular completion. Those are genuine additional
state data. In particular it does not remove the qualifications on the
[joint compact thermal comparison](quantum-junction-state.md).

## Positive spectral completion of the Gaussian phase

The phase propagator of (3) is

```
G_alpha(s)=1/[s^2 I(s)]
         =1/(I_0 s^2)+integral_0^infinity d sigma(w)/(s^2+w^2),
I_0=I(0),       d sigma >= 0.                             (5)
```

Here is a direct proof of the positivity, including the mixed axial
boundary condition. First put a reflecting cut at finite negative `x`.
Let `Z_0` minimize the zero-frequency quadratic form with unit light
boundary value. Perturb it by a function whose light boundary value
vanishes; its heavy boundary value remains free. Denote the resulting
positive operator by `H_D`. In its orthonormal eigenbasis write
`H_D u_a=w_a^2 u_a` and `c_a=<u_a,Z_0>`. Minimization at `z=s^2`
and Parseval give

```
K(s)=K(0)+sum_a [s^2 w_a^2/(s^2+w_a^2)] |c_a|^2,
d sigma(w)=sum_a w_a^2 |c_a|^2 delta_(w_a)(dw)/lambda.     (6)
```

The zero-boundary functions are dense in the spatial L2 space, so no
boundary point mass appears in Parseval. Passing the cut to the horizon
preserves the positive measure and the Stieltjes representation for
positive `s`. Equivalently these are positive variational resolvents
converging on the half-line. The possible linear-in-`s^2` remainder in
the Dirichlet-to-Neumann representation vanishes: the bounded RN
potential and the trial `Z=exp(sx)e_2` give
`s<=K(s)<=s+||V||_infinity/(2s)` at large `s`.

The earlier half-line bound supplies the quantitative estimate

```
0 <= K(s)-K(0) <= s/(1-epsilon^2),
integral s^2/(s^2+w^2) d sigma(w) <= B s,
B=1/[lambda(1-epsilon^2)],
Sigma(W):=integral_0^W d sigma <= 2 B W.                  (7)
```

For electricity `epsilon=0`; for one angular generator its value is
the mixing bound in the causal note. Equation (7) has no spectral-gap
or low-temperature premise. There is no extra zero-frequency bath atom,
since it would violate (7) as `s` tends to zero.

A finite element boundary basis has nonzero spatial norm. Consequently
its discrete version of (6) also contains a positive mesh contact term
`b_mesh s^2`. The checker retains this term explicitly, verifies the
resolvent identity with it, and checks that it decreases on refinement.
It is not interpreted as physical phase noise. The continuum proof
above, rather than a finite mesh, establishes (5).

## Hamiltonian, operator insertion and the required recoil term

At a finite frequency cutoff, (5) is realized by independent canonical
coordinates in the diagonal Hamiltonian

```
H_G=Q_0^2/(2I_0)+sum_a (P_a^2+w_a^2 X_a^2)/2,
alpha=alpha_0+sum_a sqrt(sigma_a) X_a,
[alpha_0,Q_0]=i,        [X_a,P_b]=i delta_ab.              (8)
```

Compactifying `alpha_0` is a possible completion, not a derivation of
the physical winding weights. The oscillators in (8) come from the
same radial/current system as (1); adding a second independent copy of
that system would double count them. The conserved shift charge when
there is no transfer is `Q_0`. A force `j` coupled to `alpha` gives
`Q_0(s)=j(s)/s` and `alpha(s)=G_alpha(s)j(s)`, so
`Q_0(s)=I(s)s alpha(s)` for zero initial perturbations. Homogeneous
oscillator fluctuations supply the noise omitted by a mean-response
relation. Conservation of the isolated charge is therefore consistent
with the earlier driven response.

The electric charge-one factor of a transfer insertion is `V=exp(i alpha)`.
It shifts both the zero-mode charge and oscillator momenta. In particular

```
V^dagger H_G V-H_G
  = Q_0/I_0+1/(2I_0)
    +sum_a sqrt(sigma_a) P_a+(1/2)sum_a sigma_a.           (9)
```

The last positive term is the quadratic recoil vertex; it cannot be
dropped while keeping the linear oscillator coupling. The mean energy
of a single instantaneous insertion in a centered thermal state is
`1/(2I_0)+Sigma(Lambda)/2`. A two-sided transfer has both sides'
contributions. This is the energy of an insertion, not the energy
produced by a smooth finite-time interaction: the latter weights
allowed transitions with the actual switch and fermion matrix elements.

For the oscillator factor `X=sum sqrt(sigma_a)X_a`, equilibrium quantum
mechanics gives

```
<exp(i X(t))exp(-i X(0))>=exp[J(t)],
J(t)=integral d sigma(w)/(2w) [
  coth(beta w/2)(cos(wt)-1)-i sin(wt)].                    (10)
```

The imaginary term is the commutator and is essential to recoil and
detailed balance. Its spectral transform is the energy-exchange
distribution of this Gaussian environment; its first moment is
`Sigma(Lambda)/2`, agreeing with (9). The general oscillator/Wick method
and its distinction from perturbative tunneling rates are reviewed in
[Ingold and Nazarov, sections 2--3.4](https://arxiv.org/html/cond-mat/0508728v1).
Their reservoir golden-rule calculation is not imported as the stress
of this entangled, evolving wormhole setup.

## Exact interacting bath identity

For `H=H_G+H_neutral+H_int(alpha,t)`, where the neutral operators commute
with the canonical bath variables, define the transfer current
`j=dot Q_0=-partial_alpha H_int`. No weak-coupling expansion is needed
for the oscillator equations:

```
dot P_a=-w_a^2 X_a+sqrt(sigma_a) j,
a_a(T)=exp(-i w_a T) a_a(0)+D_a,
D_a=i sqrt(sigma_a/(2w_a)) integral_0^T
        exp[-i w_a(T-t)] j(t) dt.                        (11)
```

The operators `j(t)` on the right are the full interacting Heisenberg
currents, including their dependence on the bath. Therefore (11) is an
identity, not a prescribed classical-current approximation. With
`A_a=exp(-i w_a T)a_a(0)`, the exact change of oscillator energy is

```
Delta <H_bath>=sum_a w_a [<D_a^dagger D_a>
                         +2 Re <A_a^dagger D_a>].         (12)
```

The first term is a positive current two-point functional. The second
is its mixed correlation with initial quantum noise. Even an initially
factorized neutral/bath state does not make the second term vanish:
the later current depends on the bath. Neither a mean current nor an
equilibrium vertex overlap determines (12). The corresponding identities
apply on both sides with their joint prepared state; cross-side
correlations must be retained in evaluating the currents.

The checker evolves a finite single-transfer block with the Gauss
constraint `Q=n_R` imposed, so it does not introduce an extra independent
charge rotor. Its Hamiltonian contains `sigma_+ exp(i X)+h.c.` and the
bath is fully quantum during a cubic switch. Independent unitary
evolution verifies (11)--(12), including a nonzero mixed term. A larger
oscillator truncation checks the retained result. This verifies the
operator identity beyond free evolution. It is not a discretization
of the many-channel RN fermion stress or a physical capture simulation.

## Prepared midpoint correlation and a uniform band bound

For an oscillator TFD, equal forward times `t` on the two sides give
the electric phase factor `exp[-W(t)]`, where

```
W(t)=integral_0^Lambda d sigma(w)/(2w) [
  coth(beta w/2)-cos(2wt)/sinh(beta w/2)],
W(0)=integral_0^Lambda d sigma(w) tanh(beta w/4)/(2w).      (13)
```

This is an equilibrium prepared-state factor before including the
transfer interaction. It is neither an ordinary uncorrelated-reservoir
correlator nor a prescription to reset the bath after each event.

Equation (7) implies the useful bounds

```
W(0) <= B F(beta Lambda),
F(v)=v/4                 for v<=4,
F(v)=1+log(v/4)          for v>=4,
0 <= W(t)-W(0) <= 2pi^2 B (t/beta)^2.                    (14)
```

For the first bound use `tanh(beta w/4)/(2w)` bounded by
`min(beta/8,1/(2w))` and integrate by parts against
`Sigma(w)<=2Bw`. For the second, `1-cos(2wt)<=2w^2t^2`, while
`w/sinh(beta w/2)` decreases to zero. Integration by parts and
`integral_0^infinity v/sinh(v) dv=pi^2/4` then give (14).
The free-evolution increment is ultraviolet finite even without a
cutoff. The initial bare vertex factor is not: the high-frequency
behavior `G_alpha(s)~1/(lambda s)` produces a logarithmic short-distance
singularity. A continuum charged vertex and its coupling require a
specified renormalization/matching scale.

For the selected `q=k^4`, `N_f=1`, the electric value is `B=pi/k^4`.
With `z_0=1/2`, the proposed ramp has `t/beta<=1/(8pi)`.
The flux equation and thermodynamic clock give
`beta/l_p<2 k^(27/2)`. For any retained Gaussian band with
`Lambda l_p<=1`, (14) yields, uniformly for integer `k>=2000`,

```
W(t) < (pi/k^4)[1-log(2)+(27/2)log(k)+1/32] < 2.03e-11.  (15)
```

The bracket divided by `k^4` decreases in this range. The cutoff is a
bound on the selected Gaussian band, not a claim that this truncation
remains the full four-dimensional theory up to the Planck scale.
Higher Landau levels, other harmonics and ultraviolet matching may
enter earlier. No nonlinear SU(2) vertex bound is inferred by replacing
an electric charge with a Casimir in (13).

## Why this is still not a switching-energy bound

The limitation can be proved, rather than merely cautioned about.
Take one oscillator with `d sigma=2B Omega delta_Omega`. It obeys the
full inequality in (7), since
`2Omega s/(s^2+Omega^2)<=1`. Its midpoint exponent is
`W(0)=B tanh(beta Omega/4)<=B`, but its mean recoil energy is
`B Omega`. It can therefore have arbitrarily small correlation loss
and arbitrarily large energy cost as the unresolved frequency grows.
This example satisfies the fluctuation bound; it is not an RN spectrum
or an obstruction to the chosen candidate. It proves that the bound
alone cannot supply the desired energy conclusion.

Consequently the next calculation must retain the charged vertex,
its commutator/noise kernel, the quadratic recoil term (9), and the
neutral fermion operator in the time-dependent interaction. It must
use the actual smooth switch to weight transitions and subtract the
same renormalized background energy as the JT capture calculation.
A classical average of the phase-velocity term in the prior matrix
flux would include unrenormalized zero-point energy and is not this
calculation. The angular nonlinear completion and the actual compact
state weights also remain required.

## Verification and decision

[Checker](../checks/phase_influence_action.py),
[receipt](../evidence/phase-influence-action.json), and
[source scopes](../evidence/phase-influence-source-review.json).
The checker compares independent primal and dual RN radial actions,
reconstructs their positive bath spectra including the finite-mesh
contact term, and refines the mesh. Independent finite Fock-space
traces verify (9), (10) and (13), including the Wightman commutator phase.
Three interacting quantum switches check (11)--(12), including their
nonzero mixed noise/current term and oscillator truncation refinement.
Quadratures check the measure bounds and the uniform parameter
comparison. These support the analytic arguments and do not calculate
the physical switching stress.

This turn is progress from a prescribed mean response to a derived
off-shell Gaussian action, fluctuation kernel and recoil Hamiltonian.
It also rules out the tempting inference from a small coherence error
to small switching energy. The broad formation goal remains active.
No publication, outreach, paid compute, novelty acceptance or Lean
certificate is claimed.

Remaining to-do list: calculate the renormalized interacting junction
stress with this noise/recoil kernel and the actual compact/neutral
state; complete nonlinear angular matching; bound evolving-geometry
and omitted-mode errors; establish preparation, sufficient stability
and finite-energy transport.
