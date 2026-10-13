# Smooth scalar backgrounds: static endpoint test

S069, 2026-10-12. **Numerical rejection for the two selected backgrounds.**
This is a stationary Einstein-scalar calculation, not a wormhole simulation.
No construction, formation rate, all-time stability theorem or novelty
acceptance follows. The parent planning lane is Quantyra Space, E003.

## Model and calculation

The models are S06A044 and S08A06 of
[Marks et al., 2504.17775v3](https://arxiv.org/html/2504.17775v3), equation (1)
and Table 1. That paper reports finite-duration stability evidence for the
unloaded stars. We reproduce the stationary configurations; we do not
reproduce their time evolutions or radial oscillation eigenvalues.

Set `G=c=hbar=mu_scalar=1`, write the lapse as `A=exp(nu)` and distinguish
it from the scalar amplitude `psi`. With
`ds^2=-A^2 dt^2+dr^2/g+r^2 dOmega^2`, `g=1-2m/r`, and
`Phi=psi(r) exp(i omega t)`, the action's factor of one-half gives

```text
V = psi^2 (1-2 psi^2/sigma^2)^2
V_x = dV/d(psi^2) = (1-2 psi^2/sigma^2)(1-6 psi^2/sigma^2)
T = omega^2/A^2
rho = (T psi^2 + g psi'^2 + V)/2
p_r = (T psi^2 + g psi'^2 - V)/2
p_t = (T psi^2 - g psi'^2 - V)/2
m' = 4 pi r^2 rho
nu' = (m+4 pi r^3 p_r)/(r^2 g)
psi'' = -(2/r+nu'+g'/(2g)) psi' + (V_x-T) psi/g.
```

These equations also agree with the GR reduction of
[Evstafyeva et al., 2310.05200v2](https://arxiv.org/html/2310.05200v2),
equation (13). Their appendices A/B explain the gravitational correction
to the scalar tail and the sensitivity of shooting. We use collocation
with `scipy.integrate.solve_bvp`; we do not integrate a tuned frequency
outward through the entire star.

The variables are `(m,nu,f=psi/sigma,h=f')`, with eigenfrequency represented
by a logistic parameter. Regular center series fix mass, amplitude and
amplitude derivative at `r0`. The outer clock is normalized by
`nu(R)=log(1-2m(R)/R)/2`. The decaying Robin condition is

```text
f'(R) + [kappa+(1+epsilon)/R] f(R) = 0
kappa = sqrt(1-omega^2)
epsilon = M(1-2 omega^2)/kappa.
```

This is the leading curved-space asymptotic condition, not an exact
finite-radius boundary condition. Domain changes test its effect.

## Endpoint criterion, derived before the scan

For a test string of tension `mu_string`, let `k=mu_string/m_endpoint` and
`I(r)=integral_0^r A/sqrt(g) dr`. The pair energy relative to the unbroken
string is `2U`, where `U=m_endpoint(A-kI)`. A zero-energy resting endpoint
requires both `U=0` and `U'=0`. Eliminating `k` gives

```text
D = I sqrt(g) A'/A^2 - 1 = 0,
k = A/I,
U'' = m_endpoint A^2 D'/(I sqrt(g))   at such a root.
```

Thus a strict minimum requires `D'>0`. This is not the light-ring condition
`r nu'=1`. For numerical integration use `F=I/A`, which satisfies
`F'=1/sqrt(g)-F nu'`, `F(r0)=r0-4 pi p_c r0^3/3+O(r0^5)`.
It avoids subtracting small redshifted quantities. We separately integrate
the original `I` for verification.

There is also a useful energy interpretation. At any stationary endpoint,
`k=sqrt(g) nu'` and `U/m_endpoint=-A D`. Negative `D` therefore means that
the stationary pair, if one chooses its supporting tension, has positive
energy relative to the unbroken string. It cannot be the selected
zero-energy resting endpoint. Finite-energy or thermal creation is a
different problem; it is not excluded by this test.

## Results and verification

All quantities below use `mu_scalar=1`; the particle number is
`N=4 pi integral r^2 omega psi^2/(A sqrt(g)) dr`.

| Model | Mass | r99 | N | max D | Radius of max D |
| --- | ---: | ---: | ---: | ---: | ---: |
| S06A044 | 10.428721891 | 31.747785993 | 65.883424432 | -0.4222203815 | 29.43999347 |
| S08A06 | 5.651749735 | 17.932829090 | 26.103616118 | -0.4316236298 | 15.72322144 |

The masses, radii and charges round to every displayed digit in Table 1.
Both solutions have two light rings: approximately `(28.0943,31.2382)` and
`(15.9429,16.5276)`, respectively. Neither has a zero of `D` in the numerical
profile or in its controlled asymptotic continuation.

The [receipt](../evidence/smooth-scalar-formation.json) retains four runs
per model: adaptive tolerances `1e-6`, `1e-9`, and `1e-10`, outer radii
70/90/110, inner radii `1e-3/1e-4/1e-5`, and initial meshes 901/1201/1801.
The largest observed change in `max D` is `1.68e-9`; changes in mass are
below `2.90e-7`. Each run searches all sampled local maxima and refines
them; each has one peak. The missed peak height on the 24,001-point search
grid is below `9e-8`. Independent Simpson integrals with three resolutions
agree with the `F` integration within `6.6e-13` in `D`.

We also evaluate the angular Einstein equation, which is not one of the
first-order BVP equations, using derivatives of the collocation polynomial:

```text
g (nu''+nu'^2+nu'/r) + (g'/2)(nu'+1/r) = 8 pi p_t.
```

Its largest absolute residual over all runs is below `3.76e-8` in scalar
mass units. The finest first-order pointwise residuals are below `4.6e-12`;
the finest boundary residuals are below `4.8e-17`. Residuals and convergence
are numerical evidence, not a rigorous a posteriori enclosure.

At large radius, the decaying scalar tail has rate `kappa>0.99` for these
models. For the tail estimate we allow the much looser envelope
`|psi|,|psi'| <= 1e-12 exp[-(r-R)/2]`, beginning at `R>=70`. The matching
values lie far below that envelope. The resulting mass-tail allowance is
below `1.6e-19`. The canonical scalar satisfies the radial null energy
condition, hence `(log(A/sqrt(g)))'>=0` and `I<=r`. With increasing exterior
lapse and the stated tail envelope,

```text
D+1 <= M_infinity/(R sqrt(g_R) A_R)
       + 4 pi sup_{r>=R}(r^2 p_r)/(sqrt(g_R) A_R).
```

The largest exterior upper estimate in the receipt remains below `-0.78`
for `D`. It is conditional on the asymptotic envelope and continuation on
the regular branch, rather than an interval-verified global theorem.
Near the regular center, `D=-1+O(r^2)`.

We adopt a conservative **numerical** allowance of `1e-4` on `D`, much
larger than the measured variations and integration errors. Both margins
exceed 0.42 under that allowance. We do not present the allowance as a
mathematically certified absolute error bound or a proof covering other
central amplitudes. The source-defined profile reproduction and large sign
margin justify rejecting these candidates at the numerical research gate.

## Reproduction and scope

Run from the satellite root with Python and NumPy/SciPy versions recorded
in the receipt:

```powershell
$env:OPENBLAS_NUM_THREADS='1'
python checks/smooth_scalar_formation.py --study --out evidence/smooth-scalar-formation.json
```

On the current workstation the functional interpreter is
`C:/Users/dfred/QuantyraTools/SpaceProofs/python/python.exe`.
The script solves both profiles from analytic guesses, performs the eight
comparisons, checks the rounded external benchmarks, and exports
[S06A044](../evidence/smooth-scalar-profiles/S06A044.csv) and
[S08A06](../evidence/smooth-scalar-profiles/S08A06.csv). The receipt pins the
code and profile bytes. [Source provenance](../evidence/smooth-scalar-source-review.json)
records the reading scope and cached-body hashes.

The radial compactness threshold from the rejected S068 membrane was not
used. No perfect-fluid warping theorem was applied to these smooth stars.
These hypothetical matter profiles also do not establish Standard Model
support, loaded stability, black-hole replacement, a quantum opening or a
formation rate. No universe-wide impossibility or novelty claim is made.

S069's selected candidates are closed at the endpoint gate. S070 continues
with the physical fluctuation content of the AdS formation saddle. The
broader stable-wormhole objective remains active.

Remaining to-do list: determine the admissible physical fluctuations of the
AdS tunnelling saddle, then establish a controlled formation channel and
finite-disturbance traversability if one survives.
