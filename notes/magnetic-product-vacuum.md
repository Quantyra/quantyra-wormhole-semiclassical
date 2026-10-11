# Renormalized magnetic-product vacuum and the nearby stationary reference

2026-10-10. S051/S052 feasibility result. The missing two homogeneous
pressures can be calculated from an explicit all-level spectral determinant.
They reproduce the [trace-matching condition](throat-trace-anomaly.md) and
determine both first-order radius shifts in a stated scheme. A local
implicit-function argument supplies a nearby stationary product for
sufficiently small loop strength, at each fixed magnetic charge. This is
a homogeneous reference calculation, not a globally completed wormhole,
a stability theorem, or a certified novelty claim.

## Model and spectral input

Use one free massless unit-charge four-component Dirac field on
`AdS2(a) x S2(b)`, magnetic flux `q>=1`, no electric field, and the
Einstein-Maxwell normalization of the preceding note. Define the reference
state by the complete Euclidean `H2(a)` L2 spectral kernel and analytic
continuation. Keep the magnetic flux, subtraction scale and renormalized
couplings fixed when varying either radius. No finite-radius boundary
action, boundary coupling or wormhole-specific Casimir state is included.
The result must still be matched to the actual throat and exterior state.

The sphere eigenvalues and multiplicities are

```text
lambda_n^2 = n(n+q)/b^2,
d_0=q,        d_n=2(q+2n) for n>=1.
```

The positive heat traces per physical two-volume are

```text
K_S(s;b,q) = [q + 2 sum_(n>=1) (q+2n) exp(-s n(n+q)/b^2)]/(4 pi b^2),
K_H(s;a) = (1/(pi a^2)) integral_0^infinity
           lambda coth(pi lambda) exp(-s lambda^2/a^2) d lambda,
K_4 = K_H K_S.
```

These are literature inputs. [Castellano, Montella and Zatti,
arXiv:2603.15730v1](https://arxiv.org/html/2603.15730v1), (2.29)-(2.31),
(2.84)-(2.85) at zero electric field, and (3.40)-(3.42), give the magnetic
spectrum and product factorization; their `g=q/2`. The imported equations
also agree with the [published edition](https://doi.org/10.1016/j.nuclphysb.2026.117500),
pages 6, 12 and 17. [Camporesi and Higuchi,
v1](https://arxiv.org/html/gr-qc/9505009v1), (3.77)-(3.80) and section 4.2,
independently give the hyperbolic normalization after tracing two spin
components. In particular, its prefactor is `1/(pi a^2)`, not half this.
The 2026 work already establishes determinants in this geometry; a broad
first-determinant claim is unavailable. Its supersymmetric and electric
instability results are not being imported into this different model.

At small `s`, both two-dimensional kernels start with `1/(2 pi s)`.
At large `s`, the sphere kernel tends to `q/(4 pi b^2)` and

```text
K_H(s;a) ~ 1/(2 pi^(3/2) a sqrt(s)).
```

Thus `integral K_4 ds/s` is infrared finite at zero mass: the sphere zero
modes are paired with a continuous hyperbolic spectrum, not a normalizable
four-dimensional zero mode. Radius derivatives preserve this convergence.
This statement concerns the chosen Euclidean reference prescription.

## Subtraction and the two pressures

Let `v=s/b^2`, `eta=a^2/b^2`, and `k(v)=b^4 K_4(s)`. Its local expansion is

```text
k(v) = A0/v^2 + A1/v + A2 + O(v),
A0 = 1/(4 pi^2),
A1 = (1/eta-1)/(24 pi^2),
A2 = (15q^2-3-5/eta-3/eta^2)/(720 pi^2).
```

All three coefficients follow independently by Euler-Maclaurin summation
of the sphere spectrum and the thermal expansion of `coth`. In particular,
`-A2/b^4` equals the previously computed four-dimensional anomaly.

Fix zeta subtraction and set additional finite local terms to zero in that
scheme. With `L=log(mu^2 b^2)`, the Euclidean action density is
`w_E=b^-4 F(eta,L,q)`, where

```text
I = integral_0^1 [k(v)-A0/v^2-A1/v-A2] dv/v
    + integral_1^infinity k(v) dv/v,
F = (1/2) [I - A0/2 - A1 + A2 (gamma_E+L)].
```

To define the Mellin transform first introduce a positive mass regulator.
The unregulated massless integral has no common UV/IR convergence strip;
it must not be used as an ordinary convergent Mellin transform. After the
three explicit local subtractions, the massless remainder is integrable
near `z=0` (in the strip `-1<Re(z)<1/2`). Take the mass to zero there;
the local coefficients and convergent remainders have the stated limits.
The `A2/z` pole combines with
`1/Gamma(z)=z+gamma_E z^2+...`. For the Dirac determinant
`w_E=(zeta'(0)+log(mu^2) zeta(0))/2` per physical volume. This fixes both
the finite constants and the sign; it is not an unspecified subtraction.

Continue the bulk density to Lorentzian signature with `L_quantum=-w_E`.
For the invariant stress `T_ab=p_A g_ab`, `T_ij=p_S g_ij`, varying the
two radii independently gives

```text
p_A = -w_E - (a/2) partial_a w_E
    = -b^-4 (F+eta F_eta),
p_S = -w_E - (b/2) partial_b w_E
    =  b^-4 (F+eta F_eta-F_L),
F_L = A2/2,
2p_A+2p_S = -A2/b^4.
```

These are bulk homogeneous variations for the stated prescription; no
finite-boundary stress is inferred. The dependence on `mu` is retained
through differentiation. Changing scheme requires compensating changes in
the local gravitational and Maxwell couplings.

For an independent direct pressure integral, replace the hyperbolic factor
`h(v,eta)=b^2 K_H` by

```text
h_P = h + eta partial_eta h
    = v/(pi eta^2) integral_0^infinity
      lambda^3 coth(pi lambda) exp(-v lambda^2/eta) d lambda.
```

Use the same subtraction formula with coefficients `Pj=Aj+eta Aj'`:

```text
P0=A0,       P1=-1/(24 pi^2),
P2=(15q^2-3+3/eta^2)/(720 pi^2).
```

Its finite part is `F+eta F_eta`. The program also differentiates the
original action numerically in both radii at fixed `mu` to check that this
procedure did not inadvertently vary flux, scale or physical volume.

## Ultraviolet and angular-tail control

The following are analytic remainder bounds, separate from floating-point
and quadrature error estimates. The numerical table is not an interval
certificate of the whole calculation.

For `f(n)=(q+2n)exp[-v n(n+q)]`, Euler-Maclaurin gives the sphere trace
polynomial through `N` Bernoulli terms:

```text
S_N(v) = 1/(4 pi) [2/v
  + 2 sum_(k=1)^N B_(2k)/v
      sum_(j=0)^k (-v)^(2k-j) q^(2k-2j)/((2k-2j)! j!)].
```

The Euler-Maclaurin remainder obeys

```text
|S-S_N| <= C_N exp(q^2 v/4) v^(N-1),
C_N = zeta(2N) H_N / [pi (2pi)^(2N)],
H_N = sqrt(pi)/2 * sqrt(2^(2N+1) (2N+1)!).
```

Indeed, differentiating the shifted Gaussian expresses `f^(2N)` in the
physicists' Hermite polynomial `H_(2N+1)`. Its absolute integral on the
positive half-line is at most the displayed `H_N` by Cauchy-Schwarz and
Hermite orthogonality. Combine this with
`sup |B_(2N)(fractional part)|/(2N)! <= 2 zeta(2N)/(2pi)^(2N)`
and the factor two in the sphere degeneracy sum.

The hyperbolic expansion through `v^M` is

```text
h_M = 1/(2pi v) + sum_(k=0)^M c_k v^k,
c_k = B_(2k+2)/[(2k+2) k! pi eta^(k+1)].
```

Taylor's integral remainder for `exp(-x)`, applied to the positive thermal
measure `2 lambda/(exp(2pi lambda)-1)`, yields

```text
|h-h_M| <= D_M v^(M+1),
D_M = 2 (2M+3)! zeta(2M+4)
      /[(M+1)! (2pi)^(2M+4) pi eta^(M+2)].
```

For `h_P`, its coefficient is `-k c_k` and its error constant is
`(M+1)D_M`. On `0<v<=v0`, elementary positive-kernel bounds are

```text
h <= C_h/v,      C_h=1/(2pi)+v0/(12pi eta),
h_P <= C_P/v,    C_P=1/(2pi)+v0^2/(120pi eta^2),
S <= C_S/v,      C_S=[2+3q v0+2sqrt(2v0/e)]/(4pi).
```

The last follows by bounding a positive unimodal sum by its integral and
maximum. Put `E_N=C_N exp(q^2 v0/4)` and
`C_SN=C_S+E_N v0^N`. The integrated subtracted-product error on `[0,v0]`
is at most

```text
C_h E_N v0^(N-2)/(N-2) + C_SN D_M v0^M/M.
```

Replace `C_h,D_M` by `C_P,(M+1)D_M` for the pressure channel. The script
uses `N=10`, `M=9`; it integrates the finite Laurent products exactly on
this interval before applying the overall determinant factor `1/2`.

Above `v0`, truncate the sphere sum after an integer `n_max` with
`v n_max(n_max+q)>=80`. The terms are decreasing beyond this point, so

```text
0 <= S-S_cut <= exp[-v n_max(n_max+q)]/(2pi v).
```

Using `h <= 1/(2pi v)+1/(2pi^(3/2)sqrt(eta v))` gives the integrated
mode-tail bound recorded in the receipt. The same upper bound covers
`h_P`. The remaining spectral quadratures use exponential tails; their
adaptive error estimates and floating-point cancellations are checked
numerically, not promoted to rigorous end-to-end error bounds.

## Stationary reference and its limits

The invariant Einstein-Maxwell equations reduce to

```text
-1/b^2 = -kappa q^2/(8 g_c^2 b^4) + kappa p_A,
 1/a^2 =  kappa q^2/(8 g_c^2 b^4) + kappa p_S.
```

The invariant vector current vanishes: neither homogeneous factor admits
a nonzero vector fixed by its full isotropy group. The constant magnetic
flux therefore satisfies the source-free Maxwell equation in this state.

Introduce a loop bookkeeping parameter multiplying the determinant and
hold the renormalized couplings fixed. At zero loop strength,
`a^2=b^2=r_e^2=kappa q^2/(8g_c^2)`. In variables
`u=a^2/r_e^2`, `z=b^2/r_e^2`, the classical residuals are
`(-1/z+1/z^2, 1/u-1/z^2)`. Their Jacobian at `(1,1)` is

```text
[[ 0, -1 ],
 [ -1, 2 ]],       determinant = -1.
```

For each fixed `q` and scale, the subtracted integrals and radius
derivatives are smooth in a neighborhood of positive equal radii. The UV
subtractions above give locally uniform integrability, and the IR estimate
does the same for the massless sector. The ordinary implicit-function
theorem therefore gives a unique nearby stationary product branch for
sufficiently small loop parameter in this specified model. This is a
local existence argument with no quantitative uniform charge-dependent
radius of validity, and no dynamical stability content. Fixed finite
curvature-squared terms counted at loop order do not change the first-order
shifts on the conformally flat `R4=0` reference; a different classical
higher-curvature theory requires its own stationary equations.

With `P_A=r_e^4 p_A`, `P_S=r_e^4 p_S` evaluated at the classical product,
the first-order squared-radius shifts are

```text
(b^2-r_e^2)/r_e^2 = -(kappa/r_e^2) P_A,
(a^2-r_e^2)/r_e^2 = -(kappa/r_e^2) (2P_A+P_S),
(a^2/b^2-1)       = -(kappa/r_e^2) (P_A+P_S).
```

The last reproduces the scheme-independent trace condition. The separate
shifts depend on the fixed coupling/scale definition; changing `mu` while
holding a numerical coupling unchanged compares different physical inputs.

At `eta=1` and `mu r_e=1` in the stated zeta scheme, the numerical results
are as follows. The last two columns are coefficients of `g_c^2` in the
first-order squared-radius shifts, not finite-coupling existence bounds.

| q | r_e^4 p_A | r_e^4 p_S | sphere shift/g_c^2 | AdS shift/g_c^2 |
| --- | --- | --- | --- | --- |
| 3 | -0.0245236964 | 0.0157988167 | 0.0217988413 | 0.0295542899 |
| 8 | -0.0788783453 | 0.0121048708 | 0.0098597932 | 0.0182064775 |
| 16 | -0.0928255450 | -0.1765902968 | 0.0029007983 | 0.0113200433 |
| 64 | 4.9328845225 | -9.2551477100 | -0.0096345401 | -0.0011926198 |

No general sign claim is inferred from these samples. The scale-dependent
common shift can change sign while the anomaly fixes the relative shift.

## Evidence and next decision

[The script](../checks/magnetic_product_vacuum.py) and
[receipt](../evidence/magnetic-product-vacuum.json) verify seven exact
identities and the nonzero classical Jacobian. Four charge cases compare
two UV cutoffs; their largest discrepancy is `1.57e-11`. Independent
two-radius finite differences at `q=3` agree with the direct pressure
integrals within `1.35e-11`. These residuals are numerical cross-checks.
The ordinary integral formulas, UV bounds and local existence argument
are the mathematical statements; no Lean certification has been run.

This closes the homogeneous reference-pressure feasibility calculation in
one explicit scheme and state. It does not close the inhomogeneous
wormhole's stress response. In particular, the invariant product stress
has zero contraction with an AdS2 null vector and therefore does not itself
provide the state-dependent null stress that opens an MMP throat. The
global field-line boundary state, exterior matching, response to changing
geometry, finite disturbances and physical shortcut test remain required.
The next decisive calculation is how the Casimir-supported throat and its
coupled response change about this corrected reference. The determinant
construction is established prior work; any consequential contribution
must survive the remaining proof, edition and citation comparison.

Remaining to-do list: match the quantum reference to the inhomogeneous
Casimir state and global geometry, complete primary/fallback novelty
selection, then establish and verify a consequential physical result.
