# Radial Gauss Hamiltonian and the charged boundary correlator

S072, 2026-10-13. **The spherical electric LLL theory on a fixed radial
background admits an exact boundary thermal two-point formula that keeps
the radial Coulomb profile, the fermion charge lattice and one common
oscillator subtraction.** The screened charge energy follows from
diagonalizing this Hamiltonian; it is not an additional interaction imposed
on the same free fermions. In the specified neutral-inner-cut ensemble,
the infinite optical-volume midpoint is the earlier matched Gaussian
answer times a single theta-function correction. This is an initial-state
result in a precisely defined projected theory, not a four-dimensional
opening, preparation or stability theorem.

## Model, Gauss law and ensemble

Take n massless complex Dirac channels on an optical interval [-L,0],
with reflecting ends whose unfolding is an antiperiodic circle of length
2L. Use a fixed static spherical metric, one asymptotic time, zero electric
chemical potential, and no fermion mass or nonspherical gauge fields.
Set lambda=n/pi, kappa=4pi r^2/g^2, w=f/kappa and V=lambda w>=0. On each
finite interval assume V bounded and piecewise continuous. At the inner
cut fix electric flux to zero. Sum over all integer species charges
N=(N_1,...,N_n); the outer flux is Q=sum N_m. This is a specified matter
ensemble. It does not sum over charged black-hole interiors or derive
their degeneracy, nor insert a separate horizon rotor.

The spherical projection and charged/neutral separation are established
in [Maldacena, Milekhin and Popov, appendix A](https://arxiv.org/html/1807.04726v3#A1).
The charge lattice and boundary fermion ladder follow from constructive
bosonization on the unfolded circle; see [von Delft and Schoeller,
sections 4E--4F, 6--7 and appendix B](https://arxiv.org/html/cond-mat/9805275v3).
The following radial reduction is derived here from those inputs.

In radial gauge let varphi be the canonical charged density boson, so
rho=sqrt(lambda) varphi_x and varphi(-L)=0. Gauss law gives
p(x)=integral_-L^x rho=sqrt(lambda) varphi(x), with
varphi(0)=Q/sqrt(lambda). Before gauge fixing, use the fermion transported
by its equal-time Wilson line to the outer boundary. Gauge fixing removes
that line from its expression, not its action on electric flux.
The Maxwell energy is exactly (1/2) integral w(x) p(x)^2 dx. Thus

```
H_e = (1/2) integral_-L^0 [Pi^2+varphi_x^2+V varphi^2] dx.       (1)
```

This retains the enclosed-charge profile. The fermion kinetic energy
supplies the first two terms; no second copy of its electric current is
retained. Charge-independent vacuum constants cancel in the normalized
thermal correlator. The neutral bosons remain free in this electric-only
model. Their charge energy is pi/(2L) [sum N_m^2-Q^2/n].

Let z_s solve [-d_x^2+V+s^2]z_s=0, z_s(-L)=0, z_s(0)=1, and define
K_L(s)=z_s'(0). With z=z_0, write varphi=Q z/sqrt(lambda)+y, where y
vanishes at both ends. Integration by parts removes the mixed term.
For normalized Dirichlet modes H_D u_a=omega_a^2 u_a, define
c_a=integral u_a z. Then the exact diagonal Hamiltonian is

```
beta H = N^T A N + beta sum_a (P_a^2+omega_a^2 Y_a^2)/2
         + beta H_neutral,osc + constant,
c=pi beta/(2L), x=beta K_L(0)/(2lambda),
P=1_n-11^T/n, A=c P+x 11^T.                                  (2)
```

Here 1_n is the identity matrix and 1 is the column of ones. A is
positive definite. K_L(0)>=1/L, hence x>=c/n. The screened coefficient
I_L=lambda/K_L(0) is derived in (2). Its oscillators are simultaneously
changed from the free ones; adding x Q^2 to the original free Hamiltonian
would not give (2).

## The physical boundary ladder also translates the radial oscillators

A boundary creation in species m changes N to N+e_m. At equal time it
does not change enclosed charge at any fixed interior point x<0. In
the centered coordinates in (2), it therefore translates
Y_a -> Y_a-c_a/sqrt(lambda). Its electric oscillator factor is the Weyl
translation exp[i sum c_a P_a/sqrt(lambda)], with the sign reversed for
annihilation. The neutral translation and Klein signs are the usual
ones. The sign cancels between the two insertions in this two-point
function, but must remain in a multi-species driven Hamiltonian.

This gives a direct canonical interpretation of the boundary phase
spectrum from [phase-influence-action.md](phase-influence-action.md):

```
sigma_a=omega_a^2 c_a^2/lambda,
K_L(s)=K_L(0)+sum_a s^2 omega_a^2 c_a^2/(s^2+omega_a^2).          (3)
```

To prove (3), minimize the quadratic form at s after writing its trial
function as z+v with v Dirichlet. Mode minimization gives
K_L(s)=K_L(0)+s^2 ||z||^2-s^4 sum c_a^2/(s^2+omega_a^2);
Parseval gives (3). The dense continuum Dirichlet basis leaves no
boundary mass. A finite-element basis instead has a mesh contact b s^2.
The same mesh has the same b with and without V, so it cancels in their
difference. The checker retains it rather than attributing it to noise.

For 0<t<1, set tau=beta t. A single translated oscillator has thermal
two-point factor exp[-d^2 omega F_t(beta omega)/2], where d=c_a/sqrt(lambda)
and

```
F_t(v)=(1-exp(-t v))(1-exp(-(1-t)v))/(1-exp(-v)).                (4)
```

This follows by the Weyl multiplication rule and the thermal Gaussian
trace, keeping the commutator. At t=1/2, F_t(v)=tanh(v/4).
The entire charge trace, including the neutral charge constraints, is

```
Theta_A(v)=sum_(N in Z^n) exp[-(N+v)^T A(N+v)],
C_A,m(t)=exp[-t(1-t) A_mm] Theta_A(t e_m)/Theta_A(0).            (5)
```

Indeed (1-t)N^T A N+t(N+e_m)^T A(N+e_m)
=(N+t e_m)^T A(N+t e_m)+t(1-t)A_mm. These are the original fermion
sectors, not independent total-charge and neutral rotors.

## Exact finite-interval answer and common subtraction

Let G_V,L(tau)=<psi_m(tau) psi_m^dagger(0)> in this ensemble. Normalize
the endpoint operator by the same ultraviolet convention as G_0,L in
the free massless theory. With s_j=2pi j/beta, set

```
D_L(t)=2/(beta lambda) sum_(j>=1)
       (1-cos(2pi j t)) [K_L(s_j)-s_j coth(s_j L)]/s_j^2.

G_V,L(beta t)/G_0,L(beta t)
 = exp[-D_L(t)]
   [Theta_A(t e_m)/Theta_A(0)]
   /[Theta_(c 1_n)(t e_m)/Theta_(c 1_n)(0)].                   (6)
```

Derivation: multiply (4) and (5), and divide by their free counterparts.
The charge-energy difference is t(1-t)(x-c/n). The Matsubara identity
sum_(j>=1)(1-cos(2pi jt))/j^2=pi^2 t(1-t), together with (3), absorbs
exactly that difference into D_L. All neutral oscillator factors cancel.
One must not multiply (6) by another compact factor or free charged
oscillator. For n=1 the formula still applies, without neutral modes.

For completeness, define delta K_L(s)=K_L(s)-s coth(sL). The maximum
principle gives 0<z_s,V<=z_s,0 in the interior. Differentiating the
minimum with respect to s^2 gives dK_L/d(s^2)=integral z_s^2. Therefore
delta K_L is nonnegative and decreasing in s^2, and

```
0<=delta K_L(s)<=K_L(0)-1/L,
0<=D_L(t)<=t(1-t)(x-c/n).                                    (7)
```

The first inequality also follows by dropping V in the variational
problem. It proves absolute convergence of (6), including the common
point limit. The ultraviolet kinetic normalization remains the free
one. This fixes the reduced-theory operator; it does not determine a
finite matching coefficient from a more complete four-dimensional theory.

## Horizon limit, compact correction and candidate comparison

Now take nested intervals exhausting the fixed RN optical half-line,
with bounded V>=0 and finite K(0)>0. For RN, V=f lambda/kappa is bounded,
decays at the horizon and tends to a finite value at the AdS boundary.
The positive variational problems converge to the regular half-line
solution. In particular K_L(s)->K(s) for s>=0. At s=0 this statement is
convergence of infima: compactly supported cutoffs of the regular finite-
energy minimizer supply recovery sequences, even when the minimizer has
a nonzero limiting value at a nonextremal horizon. No finite reflecting
horizon is retained in the answer.

At fixed finite n and beta, c->0 and x->beta K(0)/(2lambda)>0. Poisson
summation proves the required lattice limit, without replacing species
charges by independent rotors. The dual exponent is

```
pi^2 p^T A^(-1)p
 = pi^2 [||P p||^2/c+(sum p_m)^2/(n^2 x)].
```

Every nonconstant integer vector p is suppressed. The surviving vectors
p=l 1 have exponent pi^2 l^2/x and phase exp(2pi i l t). Dominated
convergence is justified by c<=c_0 and x in a fixed compact positive
interval: the summands are bounded by a summable positive Gaussian on
Z^n. This is a limit at each fixed n, not a finite-L estimate uniform in n.

Writing Theta_t(x)=sum_(Q in Z) exp[-x(Q+t)^2], the free endpoint tends
to 1/[2beta sin(pi t)], and (6) gives

```
G_V(beta t)= exp[-D(t)] Theta_t(x)/Theta_0(x)
             /[2beta sin(pi t)],
D(t)=2/(beta lambda) sum_(j>=1)
       (1-cos(2pi jt))[K(s_j)-s_j]/s_j^2,
0<=D(t)<=x t(1-t),    x=beta/(2 I_0), I_0=lambda/K(0).          (8)
```

The theta ratio is at most one by its positive dual Fourier series.
At the midpoint,

```
G_V(beta/2)=(1/(2beta)) exp[-D(1/2)] Theta_(1/2)(x)/Theta_0(x),
0<=1-Theta_(1/2)/Theta_0
 <=4 exp(-pi^2/x)/(1-exp(-8pi^2/x)).                          (9)
```

The Gaussian expression in the prior boundary matching note is exactly
D(1/2). The only compact replacement is the theta ratio in (9); its
Gaussian zero-charge contribution is already included in D. The
constant-charging reservoir had no radial oscillator shift; it is a
separate comparator, not the Hamiltonian used to prove (8).

For the inherited candidate, I_0>=I_bare and
x<x_bare<1/(48 k z0)<=1/(12k) throughout 1/4<=z0<=7/25.
Thus the fractional initial electric midpoint loss in this ensemble is
bounded by

```
1-2beta G_V(beta/2)
 <=1-exp[-1/(48k)]
   +4 exp(-12pi^2 k)/(1-exp(-96pi^2 k)).                     (10)
```

At k=2000 the first term is 1.04166125e-5. The logarithm of the second
term is less than -236869; report this as a logarithmic bound, not a
floating-point zero. This controls initial electric coherence only.
The source assumptions behind I_0 and beta are the retained fixed-RN
ones. Angular projection, gravitational charge-sector weights and
microscopic matching are not included in (10).

## Verification, decision and next physical step

[Checker](../checks/radial_electric_correlator.py),
[receipt](../evidence/radial-electric-correlator.json), and
[source record](../evidence/radial-electric-source-review.json).
Independent shifted-oscillator thermal matrices test the physical
boundary ladder and charge trace. An analytic constant-potential spectrum
is compared with the Matsubara boundary formula at several Euclidean
times. Nonconstant finite-element Gauss Hamiltonians test the canonical
and boundary-action reductions, retaining mesh contacts and refinement.
Direct and Poisson-dual lattice sums test the charge selection and
half-line limit. The proof of (7)--(10), not small-species numerics,
supports the parameter-family statement.

Two development checks were corrected without changing the Hamiltonian
or analytic result. The first requested the optical-limit tolerance
before the two-species neutral lattice had converged; extending the
existing sequence from c=0.3 to c=0.15 resolved it. The second omitted
the finite mesh's constant high-frequency delta-K tail. The checker now
sums that tail exactly and bounds the remaining inverse-fourth-power
tail. Both original snapshots and failure logs are retained in the
source record. Neither failure was a continuum-theorem counterexample.

This completes the radial electric initial correlator for the stated
fixed-background, spherical LLL, neutral-inner-cut ensemble. Matching that
ensemble to the actual prepared black-hole pair remains essential: a
different horizon/seed charge measure or angular interaction changes
the trace. A Euclidean equilibrium correlator does not determine the
nonstationary state after coupling. The next calculation should retain
this charged boundary operator and the physical neutral/angular branching
in the driven Hamiltonian, with initial sector weights specified, and
calculate the two-time covariance and stress rather than substituting
the midpoint amplitude into the free-JT energy. No new publication,
paid compute, outreach, Lean run or novelty acceptance is selected.

Remaining to-do list: match the prepared gravitational/angular sector
weights and boundary operator; control driven state, complete stress,
massive effects and evolving geometry; establish preparation, quantum
fluctuations, sufficient general stability and finite-energy transport.
