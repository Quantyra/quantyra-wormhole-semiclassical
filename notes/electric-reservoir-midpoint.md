# Electric midpoint correlation with a radial fermion reservoir

S072, 2026-10-13. **A global electric charging Hamiltonian admits a
cutoff-independent bound on the physical boundary midpoint, with the
radial kinetic energy and actual fermion occupations retained.** Its
relative loss is at most `1-exp(-x/2)`, where `x=beta/(2 I)`. In a
specified large radial reservoir limit its exact factor becomes the
compact U(1) rotor factor. No independent rotor is added to the fermions.

This resolves the electric constant-charging comparison that the
[canonical first moment](boundary-packet-spectral-matching.md) could
not control. It does not replace the full radial Coulomb interaction
by a total-charge square, establish the screened value of `I` in this
Hamiltonian, or include the angular interaction. Those physical matching
conditions remain essential.

## Hamiltonian and the actual charged transitions

Start with finitely many particle and antiparticle pairs, with
nonnegative dimensionless radial energies `u_r=beta E_r`. Degenerate
species are separate values of `r`. The seed has zero electric charge.
Use the full fermion Fock space and

```
Hcal = beta H = sum_r u_r(n_r,+ + n_r,-) + x Q^2,
Q = sum_r(n_r,+ - n_r,-),          x>0.                    (1)
```

There is no constraint fixing an individual boundary charge. Opposite
charges can be paired in a globally neutral two-sided TFD. The midpoint
for a tagged particle annihilator is
`b_r=Tr[exp(-Hcal/2)c_r^dagger exp(-Hcal/2)c_r]/Z`.
Its free value is `b_r,0=1/[2 cosh(u_r/2)]`.

Remove the tagged particle AND its antiparticle from the reservoir.
In the free kinetic Gibbs state of the remaining pairs, let `P_Q` be
the probability of net charge `Q`, and set

```
z=exp(-u_r),        f(s)=sum_Q P_Q exp[-x(Q+s)^2],
D=(1+z^2)f(0)+2z f(1).
```

The four tagged occupations give the exact result

```
b_r = sqrt(z)(1+z) exp(-x/4) f(1/2)/D,
r_r = b_r/b_r,0 = (1+z)^2 exp(-x/4) f(1/2)/D.             (2)
```

This includes all other radial modes in their actual fermionic state.
It is a ratio of full thermal traces, not a fixed-charge approximation.
For an independent spectral check, let `a=0,1` be the tagged antiparticle
occupation. The positive anticommutator spectral measure of the tagged
particle has an atom at dimensionless frequency

```
v_Qa = u_r+x[2(Q-a)+1]
```

with weight

```
A_Qa = P_Q z^a {exp[-x(Q-a)^2]+z exp[-x(Q+1-a)^2]}/D.
sum_Qa A_Qa=1,       b_r=sum_Qa A_Qa/[2 cosh(v_Qa/2)].    (3)
```

Frequencies may have either sign. Equation (3) retains the whole
charged spectrum and uses the established fermionic thermal kernel.
It does not subtract inequalities for divergent first moments.

## Uniform bound without a reservoir limit

Write `p_s=1/(1+exp(u_s))` and `v_s=p_s(1-p_s)` for each remaining
pair. Its charge probabilities at `-1,0,1` are `v_s,1-2v_s,v_s`.
The reservoir characteristic function is

```
phi(alpha)=sum_Q P_Q exp(iQ alpha)
          =product_s [1-4v_s sin^2(alpha/2)] >=0.          (4)
```

In particular `P_Q=P_-Q`. This determinant/product structure follows
from the finite fermion trace; see [Klich, section 2, equations (8)--(12)](https://arxiv.org/pdf/cond-mat/0209642).
We use that finite trace identity, not a noninteracting scattering
formula for the interacting Hamiltonian (1).

The Gaussian Fourier transform and (4) give

```
f(s)=(1/sqrt(4pi x)) integral_R exp[-alpha^2/(4x)]
                                      phi(alpha) cos(s alpha) d alpha,
0<f(s)<=f(0).
```

Charge reflection also gives
`f(s)=exp(-x s^2) sum_Q P_Q exp(-xQ^2) cosh(2xQs)`, hence

```
exp(-x s^2) f(0) <= f(s) <= f(0).                         (5)
```

Use the lower bound for `s=1/2` and the upper bound for `s=1` in (2).
For the reverse estimate use the other two bounds and
`2z/(1+z)^2<=1/2`. The result is

```
exp(-x/2) <= r_r <= exp(x/4)/cosh(x/2) <= exp(x/4).        (6)
```

The bound holds at every mode energy, for every number of pairs and
at every finite regulator. In particular, it does not contain an error
multiplied by the divergent equal-point CAR normalization. The unpaired
or charge-biased reservoir case does not satisfy the proof's positivity
and reflection hypotheses; the checker includes a charged-seed example
that violates the lower bound if those hypotheses are discarded.

## Physical point field and both regulator limits

For a linear boundary field built from these radial annihilators and
antiparticle creators, distinct-mode midpoint terms vanish by the
occupation symmetries of (1). Its correlator is a sum of `b_r` with
nonnegative spatial mode weights. Equation (6) therefore bounds its
ratio to the identically smeared free correlator as well.

For an explicit massless chiral reservoir, take an antiperiodic optical
circle of length `ell`, with positive momenta
`E_n=2pi(n+1/2)/ell`, both charge signs, and a fixed finite number `q`
of species. The per-species field normalization is `1/sqrt(ell)`.
Spatial Fourier smearing by `exp(-epsilon E_n/2)` gives

```
b_ell,epsilon = (2/ell) sum_(n>=0) exp(-epsilon E_n) b_n,
b_0,point = lim_(ell->infinity,epsilon->0) b_ell,epsilon,0
          = integral_0^infinity dE/[2pi cosh(beta E/2)]
          = 1/(2beta).                                   (7)
```

The smearing acts on spatial momenta, not on the shifted transition
frequencies in (3). At fixed `ell`, the infinite free charge distribution
exists since `sum_n exp(-beta E_n)<infinity`; its finite-cutoff
distributions converge. Equations (2) and (6), and the integrable free
midpoint weights, then give the ultraviolet limit by domination.
The reservoir limit below establishes the optical-volume limit too.
The uniform comparison survives removing the spatial regulator in
either order. No density matrix on an unregulated infinite Dirac sea
is used.

Thus, even before calculating that limit, the initial electric comparison
gives the physical lower bound `b_point>=exp(-x/2)/(2beta)` whenever
the stated continuum limit is taken. Restoring `kappa=2pi/beta` yields
the same free normalization `kappa/(4pi)` as the strip checkpoint.

## When the compact factor follows from the fermions

Here the rotor factor is a limit of (2), not an assumed multiplicity.
Define for the reservoir with the tagged pair removed

```
V=sum_s v_s,
dmu(alpha)=phi(alpha)d alpha / integral_(-pi)^pi phi,
m2=integral_(-pi)^pi alpha^2 dmu(alpha).
```

The measure is positive by (4). Fourier inversion implies

```
0<=P_Q/P_0<=1,
0<=1-P_Q/P_0<=Q^2 m2/2.                                 (8)
```

For `V>=1/2` a fully explicit, conservative concentration bound is

```
m2 <= min[pi^2, c_*/V],      c_*=e sqrt(2) pi^(7/2)/32.   (9)
```

Indeed `phi(alpha)<=exp[-4V alpha^2/pi^2]` on `[-pi,pi]`.
On `|alpha|<=1/sqrt(2V)<=1`, each logarithm in (4) is at least
`-2v_s alpha^2`, so `phi>=exp(-1)`. The denominator normalized
by `2pi` is at least `exp(-1)/(pi sqrt(2V))`; integrating the
Gaussian upper bound for the second moment proves (9).

Let

```
Theta_s(x)=sum_(Q in Z) exp[-x(Q+s)^2],
T_0(x)=sum_Q Q^2 exp(-xQ^2),
T_half(x)=sum_Q (Q+1/2)^2 exp[-x(Q+1/2)^2],
C_e(x)=exp(-x/4) Theta_half(x)/Theta_0(x).
```

Using (8) in the three sums in (2) gives an explicit finite-reservoir
certificate, uniform in the tagged energy:

```
d_N=(m2/2)[T_half/Theta_half+1/4],
d_D=(m2/2)[T_0/Theta_0+1/2],
1-d_N <= r_r/C_e(x) <= 1/(1-d_D),        provided d_D<1.  (10)
```

For the denominator use
`sum_Q Q^2 exp[-x(Q+1)^2]=T_0+Theta_0` and
`2z/(1+z)^2<=1/2`. No continuum approximation to the charge sum
has entered the certificate. A negative lower endpoint is simply
uninformative; (6) still applies.

As `V` tends to infinity at fixed `x>0`, (9)--(10) prove `r_r->C_e(x)`
uniformly in mode energy. For the optical reservoir in (7),
`V=q ell/(4pi beta)+o(ell/beta)`; removing one pair changes `V`
by at most `1/4`. Consequently

```
b_point = C_e(x)/(2beta).                                (11)
```

This derives the U(1) heat-kernel midpoint from the actual fermion
charge weights, while retaining the physical point-field normalization.
The limit is not uniform as `x` tends to zero at unspecified reservoir
size: (10) requires the displayed charge sums and `V`, rather than
only a large number of modes. Bound (6) needs no such size condition.

Poisson summation also gives

```
C_e(x)=exp(-x/4) [sum_n (-1)^n exp(-pi^2 n^2/x)]
                         /[sum_n exp(-pi^2 n^2/x)],
0<=1-Theta_half/Theta_0
  <=4 exp(-pi^2/x)/[1-exp(-8pi^2/x)].                    (12)
```

The last upper bound need not be small at large `x`; the exact sums
are always used. In the weak charging regime the compact correction
to the Gaussian free-particle factor is exponentially small.

## Independent finite circle formula

The same optical-circle model also has an exact bosonized expression
before taking its volume limit. The constructive fermion correspondence
retains integer species charges and charge-ladder (Klein) operators;
see [von Delft and Schoeller, sections 4E--4F, 6--7 and appendix B](https://arxiv.org/pdf/cond-mat/9805275v3).
Completeness is proved there by matching the charged fermion product
with the oscillator and integer-charge partition sum. Antiperiodic
boundary conditions set their `delta_b=1`. In our normalization the
free dimensionless Hamiltonian is
`c sum_m Q_m^2+2c sum_(m,n>0) n N_mn`, with `c=pi beta/ell`.
At the spatial origin the fermion is its charge ladder times an
oscillator vertex. The oscillator factor is independent of `Q_m`.

Adding `x(sum_m Q_m)^2` in this exact correspondence changes only
the charge quadratic form to `A=c 1_q+x 11^T`. Define
`Theta_A(v)=sum_(Q in Z^q) exp[-(Q+v)^T A(Q+v)]`.
The midpoint charge shift is by the unit vector `e_m`, giving

```
C_A,m=exp[-(c+x)/4] Theta_A(e_m/2)/Theta_A(0),
b_ell,point,m = b_ell,point,0 C_A,m/C_e(c).              (13)
```

To derive it, complete the square in
`[Q^T A Q+(Q+e_m)^T A(Q+e_m)]/2` inside the charge trace.
The unchanged oscillator trace cancels in the ratio to the free field.
The charge-ladder sign cancels between the two insertions. The result
is after removal of the common spatial regulator at finite `ell`;
one should not impose this simple factor on a finite spatial smear.

Equation (13) preserves the full integer species-charge lattice. It
does not split that lattice into independent total-charge and neutral
rotors. The one-species case is the ratio `C_e(c+x)/C_e(c)`.
Independent lattice sums for one, two and three species agree with
the original fermion transition calculation. This supplies a second
derivation and a normalization check, and explains why some finite
circle corrections are already smaller than a numerical ultraviolet
truncation error.

## Candidate comparison and physical matching still needed

If the electric coefficient in (1) is matched to the candidate's
coefficient with `x<1/(48 k z0)`, `1/4<=z0<=7/25`, then (6) gives

```
1-b_point/b_0,point <= 1-exp(-x/2) < x/2
                   <1/(96 k z0)<=1/(24 k).              (14)
```

At `k=2000` the last value is `2.084e-5`. This is an electric
constant-charging initial-correlation bound. It is not a bound on
angular loss, switching energy or the complete four-dimensional state.

The earlier [Gaussian matching result](boundary-vertex-matching.md)
already includes the free-particle factor `exp[-beta/(8I0)]=exp(-x/4)`.
If a consistent joint decomposition later identifies (11) with that
compact sector, replace that factor by `C_e(x)`. Multiplying by both
counts the same electric zero mode twice. The present proof alone
does not establish the required decomposition or a microscopic finite
vertex normalization.

In particular the actual radial Gauss energy has the form
`(1/2) integral dx w(x) Q(x)^2`, with enclosed charge `Q(x)` and
positive radial weight. A charge-neutral separated pair can change
this energy while leaving total `Q` zero. Its replacement by `Q^2/(2I)`
therefore needs an integration/matching argument. Substituting an inertia
that already includes fermion polarization into (1) with the same free
fermions cannot be assumed to perform that integration consistently.

The next physical calculation must retain that radial charge kernel
and the boundary Wilson dressing, for example in the projected abelian
bosonized theory, and establish which part of (11) survives its oscillator
and compact constraints. Joint angular branching and the driven two-time
covariance remain required. Equilibrium midpoint control does not supply
the point-split switching stress.

Evidence: [checker](../checks/electric_reservoir_midpoint.py),
[receipt](../evidence/electric-reservoir-midpoint.json), and
[source record](../evidence/electric-reservoir-source-review.json).
Direct many-body transition sums, an independent Fourier evaluation,
finite-reservoir certificates and two regulator refinements check the
ordinary argument. This is a scoped mathematical matching result;
no physical wormhole creation or novelty acceptance follows from it.

Remaining to-do list: match the full radial electric Hamiltonian and
dressed boundary operator without repeating its charge mode; include
joint angular and massive sectors; control driven covariance, stress
and geometry; establish preparation, fluctuations, general stability
and finite-energy transport.
