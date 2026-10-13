# Collective electric response and the higher-mode comparison

S072, 2026-10-13. **The electrically projected neutral kernels yield a
local damping term and an additional causal memory term before the
first return. Retaining both source derivatives cancels the apparent
switch-derivative contact.** Resumming this linear response exposes a
collective charge scale that is faster than the imposed opening ramp.
Under the current reference calibration it is not below the first
higher-Landau mass scale. The earlier slow prescribed-drive comparison
does not establish that those modes can be omitted here.

This is the next calculation after [angular transfer matching](angular-transfer-matching.md).
It uses the [electrically neutral reference](neutral-transfer-current.md)
with angular dynamics held fixed. It is second order in its normalized
coupling and linear in its phase, followed by a resummation of that
specified response. It does not bound the higher-order remainder,
complete the angular/coset projection or compute the full noisy quantum
stress. The candidate is not rejected by this frequency comparison.

## Large-channel kernels with a quantitative finite-band comparison

Write `q=n`, `mu=q/(2pi)`, `b=1/(4pi)` and use global strip time.
The inherited exact kernels have positive spectral weights `v_E(q)`
at even `E>=2` for `X`, and `w_E(q)` at odd `E>=3` for `Y`. At large
`q`, their per-channel limits are

```
v_E^infinity = q E/pi^2                  (positive even E),
w_E^infinity = q(E-1/E)/pi^2             (positive odd E).
```

The odd expression is zero at `E=1`. For every integer `q>=2` and
energy band `E<=L`, `L>=3`, a uniform coefficient bound is

```
0 <= 1-v_E(q)/v_E^infinity <= epsilon_L,
0 <= 1-w_E(q)/w_E^infinity <= epsilon_L   (odd E>=3),
epsilon_L = 9 H_L/(2q)+1/q^2
          <= 9[1+log(L)]/(2q)+1/q^2.                    (1)
```

Here `H_L` is the harmonic number. To prove this, use the inherited
coefficients `f_E(a)=[z^E]exp(a L(z))`, where
`L(z)=4 sum_(odd l) z^l/l` and `c=1/q`. Positivity gives

```
f_E(1)-f_E(1-c) <= c [L exp(L)]_E,
[L exp(L)]_E <= 16 E H_E + 4/E * 1_(E odd).
```

For odd `E`, positivity of the odd powers of `L` also gives
`f_E(c)-4c/E <= 4c^3(E-1/E)`. In the loss of `w_E`, the linear
term subtracts `4c/E`, canceling the last term of the derivative
bound. Divide by `4(E-1/E)` and use
`E^2/(E^2-1)<=9/8`. This proves the odd bound. For even `E`, simply
drop the nonnegative cross-species term subtracted from the loss;
the relative loss is at most `4H_E/q`, which implies (1).

Consequently any finite-band Wightman noise quadratic form, expressed
as a sum of positive weights times squared Fourier transforms, obeys
the same relative comparison. This is stronger than convergence at
each fixed energy. It is not a relative bound on a signed retarded
kernel near a cancellation, nor a permission to interchange an
unbounded ultraviolet limit with `q->infinity`. At the illustrative
band `L=k^2`, (1) is below `4.557e-12` at `q=2000^4`.

## A common Abel regulator and its finite retarded part

Use `chi_AA=i Theta(t)<[A(t),A(0)]>` and multiply both spectra by
the same Abel factor `exp(-epsilon E)`. For `Re(s)>0`,

```
chi_XX,epsilon(s) = (2q/pi^2) sum_even
                       E^2 exp(-epsilon E)/(s^2+E^2),
chi_YY,epsilon(s) = (2q/pi^2) sum_odd
                       (E^2-1) exp(-epsilon E)/(s^2+E^2).
```

These formulas refer to the displayed large-`q` comparison spectra.
The elementary sums are

```
sum_even 1/(s^2+E^2) = pi coth(pi s/2)/(4s)-1/(2s^2),
sum_odd  1/(s^2+E^2) = pi tanh(pi s/2)/(4s),
sum_even exp(-epsilon E) = 1/(2epsilon)-1/2+O(epsilon),
sum_odd  exp(-epsilon E) = 1/(2epsilon)+O(epsilon).
```

Subtract the *same* `C_epsilon=q/(pi^2 epsilon)` in both susceptibilities.
The finite parts are exactly

```
chi_XX^ren(s) = -mu s coth(pi s/2),
chi_YY^ren(s) = -mu (s+1/s) tanh(pi s/2).              (2)
```

The divergent common constant cancels in the combined phase response
below. A regulator-dependent *difference* of finite local subtractions
would instead change its quadratic restoring coefficient and needs
physical matching. We have not proved the finite-`q` continuum stress
counterterms by taking the large-`q` limit first.

The expansions `coth(pi s/2)=1+2 sum_(m>=1)e^(-m pi s)` and
`tanh(pi s/2)=1+2 sum_(m>=1)(-1)^m e^(-m pi s)` make the delayed
returns explicit. On sources supported after zero, before `t=pi`,

```
chi_XX^ren*f = -mu f',
chi_YY^ren*f = -mu [f'+integral_0^t f(s) ds].          (3)
```

Assume `h(0)=0` and sufficient smoothness; otherwise retain the
distributional onset terms. These are causal one-sided convolutions.
They are not pointwise evaluations of divergent sine series.

## The switch-derivative cancellation

The complete inherited source response is

```
<j> = -mu h alpha - h alpha (chi_XX*h)
                         + h [chi_YY*(h alpha)]
                         + higher orders.
```

Using (3) gives

```
<j(t)> = -mu [h(t) alpha(t) + h(t)^2 alpha'(t)
                     + h(t) integral_0^t h(s) alpha(s) ds]
           + higher orders,                    0<t<pi. (4)
```

In particular, `+mu h h' alpha` from the local restoring expectation
cancels `-mu h h' alpha` from the current-current response. Omitting
the local term leaves both the common divergent contact and an
incorrect derivative term. The memory term in (4) is the consequence
of the electrically neutral spectrum `E-1/E`; using the unprojected
free-fermion spectrum `E` would miss it. This reproduces the inherited
double-counting distinction at the level of a causal equation.

Noise has not disappeared: at this order its symmetrized covariance is

```
N_j(t,s)=h(t)h(s) Re C_YY(t-s),
C_YY(t)=q/pi^2 sum_(odd E>=3)(E-1/E) exp[-iE(t-i0)].    (5)
```

The radial electric phase noise from [phase-influence-action.md](phase-influence-action.md)
is another retained input. Neither can be replaced by a zero-point
average of the classical phase-velocity stress.

## Collective feedback, with the clocks stated explicitly

For a static-inertia comparison let `I_rel` be the inertia conjugate
to the *relative* electric phase. With the initial clock map
`t_dot=2pi/beta` and equal mouths of capacitance `K_e`,

```
x=beta/(2K_e),    I_rel=K_e t_dot/2=pi/(2x).
```

The factor two accounts for equal and opposite charge on the mouths.
Hold this initial clock and inertia fixed for this comparison. Feeding
(4) into `I_rel alpha''=<j>` sums repeated insertions of the quadratic
response without a Dyson expansion in `q h`. For constant `h`, its
local inertial/restoring/damping part has

```
Omega^2 = mu h/I_rel,
Gamma = mu h^2/I_rel = h Omega^2,
alpha''+Gamma alpha'+Omega^2 alpha
                  +Gamma integral_0^t alpha(s) ds=0.   (6)
```

For a variable switch retain the two separate `h` factors inside and
outside the integral in (4). Equation (6) omits the noises (5) and the
radial bath; it describes a homogeneous response component. The
first-return memory equation must not be continued to arbitrarily
late times while deleting the delayed returns in (2).

For comparison with the old free switch, calibrate the leading
negative reference energy as `mu h=N eta`, `N=2q`; thus `h=4pi eta`.
This is the existing reference normalization, not a derivation of the
renormalized four-dimensional coupling. Using the frozen initial RN
relation, with `C=1+6X`, `X=b^2/L_4^2`,

```
x = pi q/[8 eta z_0 C(1+3X)],
Omega^2 = q^2/[2 z_0 C(1+3X)],
m_1^2 R_AdS2^2 = (q+1)/C.                            (7)
```

The last expression is the first higher-Landau *mass parameter* in
the throat, not an exact global pair-production threshold. It follows
from the sphere Dirac eigenvalue `n(n+q)/b^2` at `n=1` and
`R_AdS2^2=b^2/C`; see [the magnetic spectrum, equations 2.29-2.30](https://arxiv.org/html/2603.15730v1).
Write this dimensionless mass as `M_1`. Then

```
Omega^2/M_1^4 = q^2 C/[2z_0(q+1)^2(1+3X)]
               -> 1/z_0                 (q,X large). (8)
```

Thus the bare collective scale grows as `M_1^2`, even when `M_1` is
large compared with an order-one prescribed ramp frequency. It is
independent of `eta` under this particular initial-temperature and
reference-energy calibration. This is a new comparison, not an
additional proof of control from a large Landau mass.

For the selected `q=k^4`, `eta=1/k`, `z_0=1/2`,

```
Omega ~ k/sqrt(6pi),     Gamma ~ 2k/3,
Gamma/(2Omega) -> sqrt(6pi)/3 > 1,
M_1 ~ sqrt(k/[6sqrt(pi/3)]).
```

The two local decay rates are
`r_+-=[Gamma+-sqrt(Gamma^2-4Omega^2)]/2`. At `k=2000`,

| Quantity in global strip units | Value |
| --- | ---: |
| `Omega` | 460.659 |
| `Gamma` | 1333.333 |
| Local slow / fast decay rate | 184.756 / 1148.577 |
| First higher-Landau mass `M_1` | 18.048 |
| Slow / fast rate divided by `M_1` | 10.237 / 63.640 |

These rates characterize a response kernel; they do not prove that a
given smooth source excites those rates with appreciable weight, that
heavy particles are produced or that the physical solution is unstable.
The memory and delayed-return terms are also absent from this two-pole
diagnostic. Nevertheless an argument that integrates out the higher
Landau levels solely because the *imposed* ramp is slow has not covered
the collective response.

The radial electric Gaussian comparison has a different scale. Its
small Laplace-frequency parameter, converted to global time, is
`a_rad=2pi I_rel/q`. At `s=Omega` it is `1.364e-5` for `k=2000`.
That positive-real-frequency comparison is not a bound at all complex
poles and says nothing about the omitted higher-Landau polarization.
Small static dressing of `K_e` does not remove the separation in the table.

## A self-consistent deterministic switch and its energy accounting

To check the feedback, take `a(t)=3(t/T)^2-2(t/T)^3` on `0<=t<=T`,
then one, with `T=0.25` and observation time `0.6<pi`; let `h=h_f a`.
This is a specified global-time comparison, not the already matched
physical-time ramp of the free JT theorem. Set `m'=a alpha` and solve

```
alpha''+Gamma_f a^2 alpha'+Omega_f^2 a alpha
                                 +Gamma_f a m=0.
```

For `E/I_rel=(alpha'^2+Omega_f^2 a alpha^2)/2`, direct differentiation
gives the exact bookkeeping identity within this equation,

```
E_f + R + W_mem = E_0 + W_rest,
R/I_rel = integral Gamma_f a^2 alpha'^2 dt,
W_mem/I_rel = integral Gamma_f a m alpha' dt,
W_rest/I_rel = integral (Omega_f^2 a'/2) alpha^2 dt.     (9)
```

The memory work is signed and must be retained. For a diagnostic
initial compact Gaussian component use
`Var(alpha_0)=x/2`, `Var(alpha'_0)=2x/pi^2` and zero symmetric cross
covariance. These are the local electric heat-kernel comparison
moments, not a derivation of the physical joint compact/radial state.

At `k=2000`, solving (9) gives
`R/(mu h_f)=2.581e-7`. The response to an initial displacement first
falls to half at `t=0.028924`, whose inverse is `1.916 M_1`.
The corresponding ratios are `1.923` and `1.928` at `k=4000,8000`:
even this smooth ramp produces homogeneous phase variation on the
higher-mode mass time scale in this comparison.

As a check against a tempting inconsistent calculation, omit damping
from the equation and then evaluate the same positive friction load on
that trajectory. It is approximately **580 times** the supplied
initial-plus-restoring energy at `k=2000`. That is not a physical
radiation prediction; it demonstrates why feedback cannot be added
only after computing the phase motion. With feedback retained, (9)
balances, including its nonzero memory term.

Nor is the damped homogeneous component a complete quantum covariance:
the checker explicitly finds its `det Cov(alpha,p)` below `1/4`.
The noise-driven contribution is needed to restore the canonical
uncertainty bound. Thus the small number quoted for `R` is not an
upper bound on the total quantum energy or on the four-dimensional
capture error. No late-time phase stability follows from this run.

## Evidence and next scientific decision

[Checker](../checks/collective_electric_response.py),
[receipt](../evidence/collective-electric-response.json), and
[source and dependency record](../evidence/collective-electric-source-review.json).
Six channel counts check (1) through energy 512. Direct Abel sums
converge to (2); independent analytic integration of the oscillatory
mode responses checks the cancellation in (4), while the bubble-only
answer diverges. Three coupled switches verify (9), and an ODE
refinement checks their numerical accuracy. The tests do not establish
an all-order interaction remainder or a continuum quantum stress.

The next response calculation must include or bound massive-Landau
and other omitted-mode load at the *collective* scale, alongside the
neutral/radial noise and common stress subtraction. An overdamped
alternative has local slow rate approximately `1/h`; making that
rate small compared with `M_1` would require `h M_1>>1`, together
with a bound on the fast transient and its spectral weight. No new
parameter family or physical endpoint is accepted from that observation.

Remaining to-do list: compute or bound omitted-mode response at the
collective scale; include quantum noise, mixed energy terms and
higher-order interaction errors; match physical clocks, vertices and
renormalized stress; control evolving geometry and establish preparation,
sufficient stability and finite-energy transport.
