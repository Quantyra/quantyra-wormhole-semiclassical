# Neutral fermion current after removing the electric charge mode

S072, 2026-10-13. **The free LLL strip determines explicit neutral
transfer-current kernels once its common charged boson is removed.**
The removal changes the current noise, rather than merely multiplying
the free-fermion answer by a coherence factor. The resulting neutral
influence functional and its second-order response include a local
term from the change of the transfer expectation value. Direct fermion
determinants and quantum evolution verify these statements.

This provides fermionic input for the [phase/noise calculation](phase-influence-action.md).
It is not yet the interacting four-dimensional stress: nonlinear
angular projection, physical vertex normalization, the coupled phase
evolution and continuum stress subtraction remain required.

## Reference state and the mode that must be removed

Use the same unfolded antiperiodic chiral circle of circumference `2pi`
and prepared strip vacuum as in the
[free switching calculation](switching-radiation-capture.md). There are
`n=N_f q` complex channels. The selected family has `N_f=1`, hence
`n=q`. The two endpoint fields are separated by `pi`. Define

```
T_m(t)=i psi_L,m^dagger(t) psi_R,m(t),
b=<T_m>=1/(4pi),
L(t)=4 sum_(E=1,3,5,...) exp[-iE(t-i0)]/E.
```

For positive Euclidean separation `tau`,
`L_E(tau)=2 log coth(tau/2)`. For real time retain the positive-frequency
series (or its analytic continuation with `|exp[-i(t-i0)]|<1`); arbitrary
principal powers on opposite sides of a cut would lose the commutator.

Abelian bosonization writes the reference transfer as a product of a
common charged vertex and a neutral vertex. If `phi_m` denotes the
relative endpoint boson, set

```
phi_c=(1/n)sum_m phi_m,       beta_m=phi_m-phi_c,
<phi_c(t)phi_c(t')>=L(t-t')/n,
<beta_m(t)beta_l(t')>=P_ml L(t-t'),
P_ml=delta_ml-1/n,           sum_m beta_m=0.               (1)
```

The electric gauge field couples to the common mode, leaving the
`SU(n)_1` sector, as described in
[MMP appendix A](https://arxiv.org/html/1807.04726v3). The radial electric
phase action already contains this common current degree of freedom.
An independent copy of the full free-fermion current kernel would count
it again. Equation (1) specifies the subtraction at the oscillator level.

This is only the electric subtraction. The remaining sector still
contains the angular `SU(2)` current algebra. Adding the three angular
Gaussian phases independently to all of these neutral fields would
repeat the same problem. The nonlinear angular/coset projection and
the compact charge-spin sector selection remain to be performed.

## Exact neutral multipoint functions

Normal order relative to the stated vacuum and normalize

```
B_m=b :exp(i beta_m):,       C_ref=:exp(i phi_c):,
T_m=B_m C_ref,              <B_m>=b, <C_ref>=1.           (2)
```

These are fixed normalization conventions, not a claim that a physical
unitary charged vertex has expectation one. For example, in a finite
Gaussian regulator `C_ref=Z_ref^(-1)exp(i phi_c)`, with
`Z_ref=exp[-L_reg(0)/(2n)]`. Replacing the reference charge mode by the
gauged compact/radial phase must retain and match this factor in the
vertex or in its coupling. Identifying both `C_ref=exp(i alpha)` and
the old coupling without that matching would be unjustified.

For any ordered product with `epsilon_a=+1` denoting `B_m` and `-1`
denoting `B_m^dagger`, Gaussian contraction gives

```
<product_(a=1)^r B_(m_a)^(epsilon_a)(t_a)>
 = b^r exp[-sum_(a<b) epsilon_a epsilon_b
                (delta_(m_a,m_b)-1/n) L(t_a-t_b)].        (3)
```

The order is the operator order, not an instruction to sort real times.
Contour-ordered products use the corresponding ordering and analytic
continuations. The charged reference factor is the same exponential
with the coefficient `1/n` in place of `delta-1/n`. Multiplying the two
restores the full coefficient `delta_(m_a,m_b)` and therefore the
original free-fermion correlator.

This can also be obtained directly by Wick's theorem. For each species,
put the creation endpoint at `x_a=0` and annihilation endpoint at
`y_a=pi` for `epsilon_a=+1`; exchange them for `-1`. At Euclidean times
the fermionic contraction matrix is

```
D_ab=-i/[4pi sin((y_b-x_a+i(tau_b-tau_a))/2)].
```

Its determinant, multiplied by `product_a i epsilon_a`, is the
ordered transfer correlator. The trigonometric Cauchy determinant
factorizes into (3) times the charged factor. This checks the sign,
the endpoint phases and the normalization without assuming that the
neutral vertices themselves are free fermions. They are not.

For a prescribed phase history, a finite-regulator Hamiltonian can now
be written as

```
H_int(t)=-h(t)[B(t)exp(i alpha(t))+B^dagger(t)exp(-i alpha(t))],
B=sum_m B_m.                                             (4)
```

Here `h` includes the chosen fixed vertex-normalization factor; it has
not been equated to the physical four-dimensional coupling or to the
earlier angle beyond a matched reference convention. Expanding the
forward/backward evolution of (4) and using (3) specifies every term
of its neutral influence functional. At finite mode cutoff the
normal-ordered vertices are bounded multiples of unitaries, so this
Dyson expansion exists for bounded finite-time `h`. Equal forward and
backward histories give one by unitarity. Removing the cutoff requires
the specified boundary renormalization; that limit is not asserted here.

The general use of source-dependent fermion actions is established;
[Snyman and Nazarov, sections I, II.6 and III](https://arxiv.org/html/0801.2293v1)
provide a scattering-action comparator. Their channel-diagonal reservoir
specialization is not used to replace this prepared strip covariance.

## Closed current and restoring-operator kernels

Define the Hermitian operators and mean

```
X=B+B^dagger,       Y=i(B-B^dagger),       mu_X=2nb.
```

For (4), the current conjugate to the phase is
`j=-partial_alpha H_int=h[Y cos(alpha)-X sin(alpha)]`.
Let `a=1-1/n`, `c=1/n`, and write `L=L(t-t')`. Equation (3) gives

```
C_XX(t-t')=<X(t)X(t')>-mu_X^2
  =4n b^2 [cosh(aL)+(n-1)cosh(cL)-n],
C_YY(t-t')=<Y(t)Y(t')>
  =4n b^2 [sinh(aL)-(n-1)sinh(cL)],
C_XY=C_YX=0.                                            (5)
```

For example
`<B(t)B^dagger(t')>=b^2 exp(-L/n)[n exp(L)+n(n-1)]`.
The second term contains the weak but numerous correlations between
different species after their common mode has been removed. Dropping
those terms would not implement the constraint in (1).

There is a positive discrete spectrum

```
C_XX(t)=sum_(E even>=2) v_E(n) exp[-iE(t-i0)],
C_YY(t)=sum_(E odd>=3)  w_E(n) exp[-iE(t-i0)],
0<=v_E(n)<=n E/pi^2,
0<=w_E(n)<=n(E-1/E)/pi^2.                                (6)
```

Proof. In positive powers of `L`, the coefficients of `C_XX` are
`[a^(2r)+(n-1)c^(2r)]/(2r)!` for `r>=1`. Those of `C_YY` are
`[a^(2r+1)-(n-1)c^(2r+1)]/(2r+1)!` for `r>=1`; the linear term
cancels exactly. For `n>=2` all are nonnegative and no larger than
the corresponding coefficients of `cosh(L)-1` and `sinh(L)-L`.
Since `L=4sum_(odd E)z^E/E` has positive coefficients, the comparison
also holds coefficient by coefficient in `z`. Finally

```
cosh(L)-1=4 sum_(even E>=2) E z^E,
sinh(L)-L=4 sum_(odd E>=3) (E-1/E) z^E.
```

The case `n=1` has no neutral field and both connected kernels vanish.
For `n=2`, `beta_2=-beta_1`, so `B=B^dagger` and `Y` vanishes exactly.
This does not eliminate physical charge noise: it is then carried by
the charged phase and its coupling to `X`. For every `n`, the neutral
current has no `E=1` oscillator. That oscillator is present in the
unprojected free current and belongs to the common charged sector.

An efficient exact coefficient recurrence is obtained from
`F_a(z)=((1+z)/(1-z))^(2a)=sum f_E(a)z^E`:

```
f_0=1, f_1=4a,
E f_E=4a f_(E-1)+(E-2)f_(E-2),
v_E=4n b^2[f_E(a)+(n-1)f_E(c)]  (E even>0),
w_E=4n b^2[f_E(a)-(n-1)f_E(c)]  (E odd>0).               (7)
```

For fixed `E` as `n` grows, the bounds in (6) are approached. In
contrast, the unprojected free current has weight `n E/pi^2` at every
positive odd `E`, including one. A small overall charge-coherence
deficit cannot substitute for this frequency-dependent subtraction.

## Response through second order, including the local term

Define the susceptibility with the sign appropriate to `H_pert=-f O`:

```
chi_AA(t)=i Theta(t)<[A(t),A(0)]>=-2 Theta(t) Im C_AA(t),
N_YY(t)=Re C_YY(t).                                      (8)
```

At linear order in an external phase and through second order in `h`,
the neutral part of the induced current is

```
<j(t)> = -h(t) mu_X alpha(t)
  -h(t) alpha(t) integral_0^t h(s) chi_XX(t-s) ds
  +h(t) integral_0^t h(s) chi_YY(t-s) alpha(s) ds
  +O(h^3 alpha, alpha^2),
N_j(t,s)=h(t)h(s) N_YY(t-s)+O(h^3).                      (9)
```

The first term is the prepared state's restoring expectation. The
second is local in the phase at time `t`, but includes the change of
`<X(t)>` caused by the switch `-h X`. The third is the response of
`Y` to `-h alpha Y`. Symmetry under `beta_m -> -beta_m` removes the
mixed `X,Y` correlators. This derivation keeps both source derivatives;
retaining only the `Y,Y` bubble misses the second term.

Equations (8)--(9) are distributions in the continuum. The local and
nonlocal contributions must share one regulator and subtraction.
They are ordinary convergent integrals in the retained finite-band
checks. Neither an isolated formally divergent local term nor a
zero-point phase-velocity average is a physical switching energy.

The checks evolve the neutral Hamiltonian (4) in a finite oscillator
band with `n=3`, varying a small external phase. The complete expression
(9) approaches the measured response with a third-order coupling
remainder. Omitting the local `X` response leaves a second-order error.
Increasing the oscillator truncation verifies the comparison. This
tests the interacting source-response structure, not just a Gaussian
two-point identity.

## Physical consequence and next calculation

The candidate is retained, but its electric phase must be coupled to
the projected neutral functional (3), not to an independent copy of
the complete free-fermion noise. The formulas now supply the restoring
term, current noise, retarded kernel and the local response needed to
start evaluating the exact bath-current energy identity in the prior
note. They do not yet evaluate that identity in the coupled physical
state. Small `h` per channel alone is not a uniform error estimate at
the candidate's large collective coupling; its phase evolution may
require resummation.

All frequencies in (5)--(9) refer to unit global strip time. They
cannot be compared directly with the RN bath's asymptotic-time
frequencies without the common boundary trajectory. The full reference
transfer bilinear has dimension one; its neutral and charge factors
carry weights `1-1/n` and `1/n`. Assigning the full redshift weight to
both would again count part of the field twice. Gauge-induced changes
of this matching, the compact sector weights, the angular projection,
and the continuum stress counterterms remain required.

[Checker](../checks/neutral_transfer_current.py),
[receipt](../evidence/neutral-transfer-current.json), and
[source scopes](../evidence/neutral-transfer-source-review.json).
The verification compares full fermion Wick determinants with the
factorized multipoint functions; proves and checks positive spectra;
tests non-time-ordered oscillator vertex products; and compares exact
finite-band quantum response with both the complete and incomplete
second-order expressions. The large-`n` decision uses the analytic
bounds, not an extrapolation from a small number of species.

The previous turn was progress on the phase action and recoil. This
turn is progress on the fermionic influence itself: it removes a
double-counting path and supplies explicit electrically neutral
current kernels. No full formation, novelty acceptance, publication,
outreach, paid compute or Lean certificate is claimed.

Remaining to-do list: use the projected influence in the coupled quantum
phase/stress calculation with matched clocks, vertex normalization and
counterterms; complete angular/coset and compact-state matching; bound
evolving-geometry and omitted-mode errors; establish preparation,
sufficient stability and finite-energy transport.
