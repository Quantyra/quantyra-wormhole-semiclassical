# Initial gravitational fluctuations and the preparation tolerances

S072, E003, 2026-10-13. **The disk-Schwarzian Hartle-Hawking state has
explicit initial length and momentum bounds compatible with the scales
of the existing classical preparation window.** At the central nominal
redshift in the candidate family, separate Born probabilities of leaving
the redshift and velocity windows are below `1.172e-8` and `2.602e-5`.
The midpoint mean correction is below `9.3e-12` throughout the nominal
redshift interval. These are ordinary analytic bounds, with independent
numerical checks below.

This is a pure-JT initial-state result, also applicable to the invariant
matter-vacuum sector specified below. It is neither a joint classical
capture probability nor quantum evolution through the switch. The
gravitational observable does not commute with its velocity. Adding the
charged matter state, its determinant and constraints is still necessary;
the answer must not be multiplied by the radial matter correlator as an
independent correction. No preparation from ordinary matter or topology
creation follows.

## State and normalization

Let `C>0` be the single-boundary Schwarzian coefficient in physical time
`U`, with action `-C integral Sch(tan(pi f/beta),U)dU`. Define

```
h=beta/(2C),  a=h/4=beta/(8C),
dmu(k)=2k sinh(2pi k)dk,
Z_h=(pi/h)^(3/2) exp(pi^2/h).
```

Use the full prefactor in [Mertens, Turiaci and Verlinde, equation
(4.10)](https://arxiv.org/html/1705.08408v3#S4.SS2), together with their
`C=1/2` convention and equations (1.8)-(1.12). Restoring `C` gives the
midpoint bilocal

```
G_Delta = 1/[2pi^2 Z_h (2C)^(2Delta)]
          integral dmu(k1)dmu(k2) exp[-h(k1^2+k2^2)/2]
          product_(signs) Gamma(Delta +/- i(k1 +/- k2))/Gamma(2Delta).
G_Delta,classical=(pi/beta)^(2Delta).                    (1)
```

The source drops overall constants in some later formulas; they cannot
be dropped in this comparison. We keep its convention removing the
constant Casimir zero-point energy. An overall topological factor cancels
in normalized disk expectation values. No sum over additional topologies
or microscopic completion is assumed.

A length representation makes the second moment unambiguous. Write
`p=-i partial_ell` on `L2(R,dell)`. The two-boundary Hamiltonian is

```
H_pair=(p^2+4 exp(-ell))/C,    H_L=H_R=H_pair/2,
psi_k(ell)=sqrt(2)/pi K_(2ik)(4 exp(-ell/2)),
<psi_k,psi_j>=delta(k-j)/[2k sinh(2pi k)],
Psi_beta(ell)=Z_h^(-1/2) integral dmu(k) exp(-h k^2/2) psi_k(ell).
                                                               (2)
```

The canonical length Hamiltonian and Bessel modes are those of
[Harlow and Jafferis, equations (39)-(44)](https://arxiv.org/html/1804.01081v2#S3.SS1),
with `C=2 phi_b`. Here the real mode normalization in (2) is fixed by
the Kontorovich-Lebedev identity. The disk norm is the norm of a
two-boundary state, not a thermal trace over a factorized pure-JT bulk
Hilbert space; see that source's section 4.

For completeness the normalization and insertion can be checked
directly. Put `x=4 exp(-ell/2)`, so `dell=-2dx/x`. Bessel orthogonality
gives the delta normalization in (2). The convergent Mellin integral is

```
integral_0^infinity x^(2Delta-1) K_(2ik)(x)K_(2ij)(x) dx
 = 2^(2Delta-3) product Gamma(Delta +/- i(k +/- j))/Gamma(2Delta).
```

Consequently multiplication by
`B_Delta=exp(-Delta ell)/C^(2Delta)` has matrix element equal to the
Gamma product divided by `2pi^2 (2C)^(2Delta) Gamma(2Delta)`. Inserting
(2) reproduces (1). This fixes the real relative phases of the state as
well as its norm. In this representation `B_(1/2)^2=B_1` is an actual
positive-operator identity, not a guessed four-point factorization or
a coincident product of matter fields. `Psi_beta` is real; its Gaussian
spectral decay and the exponential wall give finite moments used here
and zero mean momentum. The moment and momentum integrals below also
verify their finite quadratic forms.

## Exact one-dimensional midpoint formulas

Set `F=G_(1/2)/(pi/beta)` and `R=G_1/(pi/beta)^2`. Define two normalized
even densities and their Gaussian moments:

```
s(v)=sech(pi v),              w(v)=2v/sinh(pi v),  w(0)=2/pi,
S_n(a)=integral_R v^n exp(-a v^2)s(v)dv,
W_n(a)=integral_R v^n exp(-a v^2)w(v)dv.
```

Then the exact identities, for every `a>0`, are

```
F(a)=(1+2a/pi^2)S_0(a)-(4a^2/pi^2)S_2(a),
R(a)=(1+6a/pi^2)W_0(a)-(4a^2/pi^2)W_2(a).             (3)
```

To derive them rather than posit an asymptotic expansion, change
variables in (1) to `u=k1+k2`, `v=k1-k2`. The domain is `u>|v|`, the
Jacobian is `1/2`, and the Boltzmann factor is `exp[-a(u^2+v^2)]`.
The half-dimension Gamma product is
`pi^2/[cosh(pi u)cosh(pi v)]`; at dimension one it is
`pi^2 u v/[sinh(pi u)sinh(pi v)]`. Use respectively

```
sinh(pi(u+v))sinh(pi(u-v))
 = cosh(pi u)^2-cosh(pi v)^2
 = sinh(pi u)^2-sinh(pi v)^2.
```

The difference multiplied by `u^2-v^2` is symmetric under `u<->v`.
Evenness extends the wedge to a quarter of the plane. The two terms
give equal integrals after exchange of `u,v`. The unnormalized numerators
therefore become

```
N_half=(pi^2/4) integral_R2 (u^2-v^2)cosh(pi u)/cosh(pi v)
                              exp[-a(u^2+v^2)]du dv,
N_one =(pi^2/4) integral_R2 (u^2-v^2)u v sinh(pi u)/sinh(pi v)
                              exp[-a(u^2+v^2)]du dv.
```

Performing the Gaussian `u` moments and dividing by the exact norm and
classical midpoint gives (3). All steps converge absolutely for `a>0`.

## Bounds with a controlled remainder

The density moments are

```
integral s=1, integral v^2 s=1/4, integral v^4 s=5/16,
integral w=1, integral v^2 w=1/2, integral v^4 w=1.
```

They follow, for example, by differentiating the elementary Fourier
transform of `sech(pi v)` and integrating the positive exponential series
of `1/sinh(pi v)` on the positive half-line. The pointwise inequality
`0<=exp(-x)-1+x<=x^2/2` gives a remainder bound before integration.
Let `gamma=1/4-2/pi^2>0`, and write

```
F=1-gamma a+r_half,        R=1+(6/pi^2-1/2)a+r_one.
```

For `0<a<=1/10`, (3) and the positive density moments imply

```
-3a^2/(2pi^2) <= r_half <= (5a^2/32)(1+2a/pi^2),
-5a^2/pi^2    <= r_one  <= (a^2/2)(1+6a/pi^2).
```

For example, substitute `S_0=1-a/4+r_0`,
`0<=r_0<=5a^2/32`, and `0<=S_2<=1/4` in (3); the first remainder is
`(1+2a/pi^2)r_0-a^2/(2pi^2)-4a^2 S_2/pi^2`. The second follows in
the same way from `W_0=1-a/2+rho_0`, `0<=rho_0<=a^2/2`, and
`0<=W_2<=1/2`. In particular

```
|r_half|<=a^2/5,  |r_one|<=3a^2/5,
0<1-F<=gamma a+a^2/5,
0<=R-2F+1<=2a/pi^2+a^2.                              (4)
```

The last expression is the second moment centered on the *nominal*
classical value, which is the relevant center for the preparation window.
The actual variance is `R-F^2`, and differs by `(1-F)^2`.

The exact norm gives
`<k^2>=pi^2/(16a^2)+3/(8a)`. Since the potential moment in (2) is
`4<exp(-ell)>=pi^2 R/(16a^2)`, we also have

```
<p>=0,
<p^2>=3/(8a)+pi^2(1-R)/(16a^2)
     <=pi^2/(32a)+5/16.                               (5)
```

This subtracts two terms using an analytic remainder bound; a
double-precision subtraction at the candidate parameters would not
justify the result.

## Transfer to the existing preparation scales

Use the clock convention in [switching-radiation-capture.md](switching-radiation-capture.md),
`u=N U/C`, and set `M=N eta`. At equal global boundary times the
classical bilocal is `B_(1/2)=(N/C)y/2`, with `y=exp(phi)`.
Thus the canonical length variables corresponding to the initial
redshift and velocity in [returned-switch-capture.md](returned-switch-capture.md)
are

```
phi=-ell/2+log(2/N),       p_phi=-2p,
z_hat=2 exp(-ell/2)/M,     v_hat=-p/M,
[z_hat,v_hat]=i z_hat/(2M).                            (6)
```

Indeed `(C/N)H_pair=p_phi^2/(4N)+N exp(2phi)` before coupling, so
`phi_u=p_phi/(2N)` and `v=phi_u/eta` agree with the accepted classical
normalization. This identifies initial operators in the pure symmetric
sector; it does not quantize the prescribed global-time matter drive.

Choose a nominal `z_c` by `beta=2pi C/(M z_c)`. Then

```
a=pi/(4M z_c),
<z_hat>=z_c F,    <z_hat^2>=z_c^2 R,
<(z_hat-z_c)^2> <= z_c^2(2a/pi^2+a^2),
<v_hat^2> <= [pi^2/(32a)+5/16]/M^2.                    (7)
```

For `q=k^4`, `N=2q`, `eta=1/k`, integer `k>=2000`, the collective
scale is `M=2k^3`. All candidate nominal redshifts in `[1/4,7/25]`
have `a<=pi/(8 k^3(1/4))`; equation (4) gives the uniform mean loss
`1-F<9.3e-12`. This bound concerns this pure gravitational bilocal,
not an interacting fermionic vertex.

For the preparation-tail comparison choose the **interior** nominal
point `z_c=53/200`, at distance `d=3/200` from both ends of that
interval. At `k=2000`, `a=pi/(16e9*1.06)`, approximately `1.852e-10`.
The spectral Markov inequalities applied separately to the two
self-adjoint observables give, uniformly for `k>=2000`,

```
Prob_z(|z_hat-z_c|>=3/200) < 1.172e-8,
Prob_v(|v_hat|>=1/2000)   < 2.602e-5.                  (8)
```

Both bounds decrease with `k`, as follows directly from (7): their
terms scale as `1/M` and `1/M^2`. The interior choice matters. A packet
centered exactly on an edge of the classical interval has no positive
distance to that edge, so (8) cannot be asserted there.

Equation (8) refers to separate measurements on the same initial state.
There is no common spectral projection defining a sharp classical
rectangle for (6), and no justification here for applying a classical
union bound to obtain capture probability. In particular the initial
state spreads during evolution, the switch must be specified using
physical boundary clocks, and its quantum matter stress is an operator.
The accepted classical all-time trajectory bound does not answer these
quantum questions. The result removes an initial *scale* objection in
pure JT; the quantum capture theorem remains open.

## Physical state matching and verification

There is one justified matter extension. In minimally coupled JT plus
an AdS2 quantum field theory with invariant vacuum `Omega`, the reduced
Hilbert space and boundary Hamiltonians are given in
[Penington and Witten, sections 2.3 and 3.2](https://arxiv.org/html/2301.07257v2#S2.SS3).
In our variables their sum is

```
H_pair=[p^2+(J_3,matter)^2/4
        -2 J_1,matter exp(-ell/2)+4 exp(-ell)]/C.
```

Since all three matter generators annihilate `Omega`, this sector
reduces exactly to (2). Their half-disk state without insertions is
`Psi_beta tensor Omega` in this reduced representation. Thus (4)-(8)
also describe the initial gravitational observables of that state,
including the uncoupled free-field vacuum benchmark. This use of the
tensor product follows after imposing the gravitational constraints.
It requires an AdS2-invariant matter theory without direct dilaton
coupling and the stated vacuum. A boundary interaction or matter
insertion can leave this sector.

The four-dimensional charged-matter state is not supplied by the disk
norm. [Iliesiu and Turiaci, section 5.1](https://arxiv.org/html/2003.02860v3#S5.SS1)
leave additional charged scalar/fermion contributions outside their
analyzed near-extremal partition function. Here the magnetic LLL channels
are explicit dynamical matter. Their determinant, charged boundary
matrix elements, angular constraints and horizon gluing must be derived
in one common gravitational functional. Tensoring a pure-JT state with
an independent rotor or multiplying equilibrium suppression factors
would bypass that requirement. Equations (1)-(8) are a benchmark for
that matching, not its completion. In particular the radial coefficients
in [radial-axial-correlator.md](radial-axial-correlator.md) depend on the
four-dimensional radius, hence on the spherical area/dilaton when
gravity is dynamical. Applying the minimally coupled vacuum result
requires a controlled throat reduction and an error bound for that
dependence, along with its boundary conditions and flux sectors.

The [checker](../checks/schwarzian_initial_state.py),
[receipt](../evidence/schwarzian-initial-state.json) and
[source/dependency record](../evidence/schwarzian-initial-source-review.json)
compare three independent representations: the direct two-energy Gamma
integral, Bessel matrix elements, and the length wavefunction with
momentum computed by differentiating it. The length check first performs
the spectral Gaussian integral using the integral representation of
`K_(2ik)`, then integrates its square; it does not reuse (3) for the
wavefunction. Mesh refinement and 70-digit moment/candidate calculations
supplement the inequalities. Finite numerical agreement is not interval
certification or evidence for four-dimensional stability.

Next derive the charged matter/gravity state and its boundary operator
with common gauge constraints, then its evolution under a drive defined
in physical boundary time. The noncommuting angular and massive sectors
remain mandatory. A small initial pure-gravity mean correction alone
does not bound driven stress, positive area or finite-energy transport.

Remaining to-do list: full physical state/operator matching; noncommuting
angular and massive interactions; driven quantum capture, complete stress
and geometry; initial-pair preparation, sufficient stability and
finite-energy signal transport.
