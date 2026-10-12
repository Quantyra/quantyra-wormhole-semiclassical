# Interacting mouth pressure in a fixed frequency band

2026-10-11, E003/S051-S053. In the product-throat Gaussian model, the
actual moving-mouth pressure has a uniform finite-band estimate, including
arbitrarily long and unequal exterior returns. Its contribution to the
statically subtracted **band-filtered** response is small relative to
`q delta` in the retained scaling. The proof includes soft modes and the
large charged potential at the cut; it does not assume a spectral gap or
replace the channels by independent charged and neutral circles.

This closes the resonant-band pressure estimate. It does **not** bound the
renormalized high-frequency remainder or the mass-dependent local terms.
Consequently the full interacting response, physical damping and novelty
remain unresolved. The [complete prescribed free response](free-proper-mouth-response.md)
and [canonical interacting response](singlet-complex-response.md) retain
their previous scopes.

## 1. Model, observable and theorem

Use q scalar circles, each with the same throat `[-A,A]` and its own
return length `D_j>0`. The Hilbert norm is the sum of their ordinary
optical L2 norms. Matching at both mouths is continuity of field and
derivative in each channel. In the throat put

```text
v=(1,...,1)/sqrt(q), H=-partial_z^2+V(z) vv*,
V(z)=a sec^2 z, A=acos(epsilon), T=2A, b=sqrt(1-epsilon^2),
0<epsilon<=min(a,1/4), a<=1, M=a/epsilon>=1, r=epsilon/b.
```

There is no potential on the returns. Each fixed model has bounded V
and compact graph. Work in its Gaussian oscillator vacuum, fixing the
q-1 constant neutral zero modes as in the preceding notes. Let `f_n`
be real orthonormal positive eigenmodes, with frequencies `omega_n`.
No lower bound on the first positive frequency is assumed uniformly in
the family. The conserved quadratic form is

`Q(f)=integral |f'|^2 + integral_th V |v*f|^2`.

The [derived on-shell vertex](moving-mouth-clock-vertex.md) is `O=J+P`,
where `J=H_th` and

`P=r[T_zz(-A)+T_zz(A)]`,
`T_zz=(|pi|^2+|f'|^2-V|v*f|^2)/2`, with the inside mouth limits.       (1)

For any quadratic observable X let `x_nm` be its symmetric ordered-pair
amplitude, normalized so `sum_nm |x_nm|^2` is its vacuum pair weight.
In particular if `X=(pi B pi+f C f)/2`,

`x_nm=(C_nm-omega_n omega_m B_nm)/sqrt(8 omega_n omega_m)`.            (2)

This convention handles diagonal pairs as well as off-diagonal pairs.
Define `N_X(R)=sum_(omega_n+omega_m<=R) |x_nm|^2` and
`A_X(R)=sum_(omega_n+omega_m<=R) |x_nm|^2/(omega_n+omega_m)`.
Constants below depend only on the fixed allowed range of T. For all
q, all positive returns and all R>0,

```text
N_P(R) <= C R^2 B_R,              A_P(R) <= C R B_R,
B_R=1+q epsilon^2(1+R^2)^4.                                      (3)
N_J(R) <= C q R^2(1+Lbar R),     A_J(R) <= C q R(1+Lbar R),
Lbar=T+q^-1 sum_j D_j.                                            (4)
```

The polynomial in (3) is deliberately nonsharp. Its useful assertion
is uniformity at fixed R, including the strong mouth barrier `M`.

Let `0<=chi(E)<=1` be a fixed smooth pair-energy cutoff supported in
`E<=R`. For X define only the band-filtered subtracted response

`Delta_X^chi(w)=sum_nm 2w^2 chi(E)|x_nm|^2/[E(w^2-E^2)]`,
`E=omega_n+omega_m`, `w=u+i eta`, `0<c0<=u<=c1`, `0<eta<=1`.

Then

`|Delta_(J+P)^chi-Delta_J^chi|`
` <= C R(1+eta^-1)[B_R+sqrt(q(1+Lbar R)B_R)]`.                       (5)

The same bound holds for a sharp band when comparing J+P with J on
the same model. Smoothness is needed for the canonical mass comparison
used in section 6, not for (3)-(5).

## 2. The charged trace is suppressed by boundary capacity

Write `u=v*f` in the throat and `U=(u(-A),u(A))`. At inward distance
`0<=d<=epsilon` from a mouth,

`V=a/[epsilon cos d+b sin d]^2 >= a/(4 epsilon^2)`.

The one-dimensional trace inequality gives

`|u(0)|^2 <= (2/epsilon) integral_0^epsilon |u|^2`
`                  +2epsilon integral_0^epsilon |u'|^2`.

The two mouth layers are disjoint. Multiplying by `a/(8epsilon)` and
adding their contributions proves

`(M/8)|U|^2 <= Q(f)`.                                             (6)

It follows by form duality that `ell H^+ ell* <= (8/M) I_2`, where
ell is the charged two-mouth trace and H^+ is the inverse on positive
modes. The trace annihilates the exact zero modes. For `omega_n<=R`,

```text
S_U=sum U_n U_n*/omega_n <= (8R/M) I_2,
S_pi=sum omega_n U_n U_n* <= (8R^3/M) I_2.                          (7)
```

These use `1/omega<=R/omega^2` and `omega<=R^3/omega^2`; no gap is
inserted. Also, integration from either endpoint bounds
`integral_th |u|^2 <= C(|U|^2+integral_th |u'|^2) <= C Q(f)`.
Consequently `J_ch H^+ J_ch*<=C I` on throat L2, where `J_ch f=v*f`.
This is a local charged inverse bound, not a global graph gap.

## 3. Normal derivatives via a positive Dirichlet-to-Neumann map

Let Gamma be the Poisson operator for `-u''+(V+1)u=0` on the throat
with its two endpoint values, and let N be its outward
Dirichlet-to-Neumann matrix. The maximum principle on the two positive
endpoint solutions gives `||Gamma||<=C`. A linear interpolating trial
function in its variational principle gives

`0<=N<=C(1+integral_th V)I<=CM I`, since `integral_th V=2ab/epsilon`.

For a global eigenmode the charged throat restriction obeys
`(-partial_z^2+V+1)u_n=(omega_n^2+1)u_n`. Green's formula with
the two columns of Gamma gives the **exact** outward trace relation

`p_n=N U_n-(omega_n^2+1)Gamma* u_n`.                               (8)

The sign follows by testing a zero-boundary inhomogeneous solution:
its outward derivative is minus the pairing of the source with Gamma.
Combining (7), the local charged inverse bound and (8) yields

`S_p=sum_(omega_n<=R) p_n p_n*/omega_n`
`                 <= C R[M+(1+R^2)^2] I_2`.                       (9)

This estimate retains large boundary derivatives instead of assuming
they stay O(1). Applying (2) and the Hilbert-Schmidt covariance bound to
the three charged pressure terms proves

```text
N_Pch(R) <= C r^2 (||S_pi||_HS^2+||S_p||_HS^2+V_b^2||S_U||_HS^2),
V_b=a/epsilon^2,
N_Pch(R) <= C R^2[1+epsilon^2(1+R^2)^4].                          (10)
```

In particular `r V_b ||S_U||<=C R`: the large potential at the cut
does not imply a large low-frequency charged pressure weight. The
derivative contribution uses `r^2 M^2<=C a^2<=C`.

## 4. Neutral point covariances without independent-circle assumptions

Fix disjoint throat segments of length h=1/2 just inside both mouths,
using inward coordinate d. Choose a smooth real chi supported there
with integral one, and `eta_0(d)=integral_d^h chi(s)ds`. For a neutral
component y of any global eigenmode, `-y''=omega^2 y` on the throat.
Two integrations by parts give

```text
y(0)=integral chi y - integral eta_0 y',
y'(0)=integral chi y' + omega^2 integral eta_0 y.                   (11)
```

For any neutral test vector at the two mouths the corresponding
gradient functional g has `g H^+ g*<=C` by Cauchy-Schwarz in Q; it
annihilates the global zero modes. The L2 functionals l in (11) have
`sum_n |l(f_n)|^2<=C` by Bessel. Applying these operator statements to
(11), with `omega<=R`, proves

```text
sum omega_n Y_n Y_n* <= C(R+R^3) I_(2(q-1)),
sum Y'_n Y'_n*/omega_n <= C(R+R^3) I_(2(q-1)).                    (12)
```

Here Y collects all neutral endpoint values. The constants are
independent of q because the channel projection is orthogonal and
the two supporting segments are disjoint. Taking Hilbert-Schmidt
norms costs their dimension, O(q), rather than its square. Thus

`N_Pneutral(R)<=C q r^2 R^2(1+R^2)^2`.                            (13)

Unequal returns mix the charged and neutral modes globally. Their
pair amplitudes need not be orthogonal. The triangle inequality in
the pair Hilbert space, squared with a factor two, combines (10)
and (13) and proves the first half of (3). There is no assertion of
factorization of the vacuum into charged and neutral circles.

## 5. Infrared weights and canonical cross terms

For each graph there are finitely many modes below R. Equations (3)
and (4) imply their cumulative weights vanish quadratically at zero,
uniformly in the graph. Stieltjes integration by parts gives

`A_X(R)=N_X(R)/R+integral_0^R N_X(t)/t^2 dt`.

This proves their stated A bounds without dividing by a smallest
positive frequency. To prove (4), use the energy-space contraction
X for the throat Hamiltonian from the canonical-response note. Its
pair weight is `(omega_n omega_m/2)|<e_n^-,X e_m^+>|^2`.
Bessel bounds the sum over m by one, and `omega_n omega_m<=R^2`.
Min-max against the free disjoint circles, including their q zero
eigenvalues before removing the actual q-1 zeros, gives

`# {0<omega_n<=R} <= q+(R/pi)sum_j(T+D_j)`.

No inaccurate q-1 free-channel substitution is made.

Finally `|2w^2/(w^2-E^2)|<=C(1+eta^-1)` on the stated rectangle.
Use the positive pair measure `chi(E)/E` and Cauchy-Schwarz:

`|Delta_(J+P)^chi-Delta_J^chi|`
` <= C(1+eta^-1)[A_P(R)+2sqrt(A_J(R) A_P(R))]`.

This proves (5), including the real dispersive part of this finite band.

## 6. Consequence and the remaining ultraviolet question

Fix R and chi, with chi one on the resonance interval. The canonical
interacting-versus-free comparison remains
`C_chi(1+eta^-4)` with this smooth filter: insert chi into the scalar
test function in the dyadic functional-calculus proof of the
[canonical result](singlet-complex-response.md). On finitely many
nonzero annuli its first four derivatives only change the constants;
the infrared summation is unchanged (a filter smooth at zero suffices).
The high-frequency annuli vanish. That argument is independent of q,
the individual returns and the size of V in each fixed model.

If Lbar stays bounded, `epsilon=q^-2`, `a=1/(4pi^2 q)` and
`eta` is comparable to `delta=log^-2 q`, the difference between the
proper interacting **band response** and the canonical free band
response, divided by `q delta`, is bounded by

```text
C_chi [1/(sqrt(q) delta^2)+1/(q delta^5)] -> 0.                    (14)
```

The harmless lower terms are included in this upper bound. Its
two leading scales are `log^4(q)/sqrt(q)` and `log^10(q)/q`.
This is an asymptotic statement with unspecified constants, not a
usable finite-q threshold. The strong-barrier hypothesis holds for
all sufficiently large q in that family.

Equation (14) does not license removal of chi. In particular, (3)
has insufficient high-frequency decay to control the full dispersion
integral, even after subtracting the static value. The previous free
anomaly cancellation is not a proof of the corresponding interacting
local matching. The next mathematical target is a uniform estimate
for the renormalized ultraviolet difference together with its fixed
mass-dependent contact terms. The physical Einstein constraints,
angular/compact state and exterior matching remain additional tasks.

## 7. Checks, prior work and decision

The [verification script](../checks/interacting_mouth_pressure_band.py)
and [receipt](../evidence/interacting-mouth-pressure-band.json) check
the pair normalization, trace and Green identities, graph spectra,
pressure amplitudes, positive-measure cross bound and asymptotic rates.
Conforming graded finite elements resolve the large mouth potential;
independent scalar ODE calculations test the boundary map and a lowest
eigenfrequency. Mesh refinements are supplementary evidence, not proofs
of the uniform estimates. Exact zero modes are removed, and long
unequal returns are included.

Pressure fluctuation spectra, Gaussian two-frequency contractions and
their use in motional linear response are established; see Jaekel and
Reynaud, [Motional Casimir force](https://arxiv.org/abs/quant-ph/0101069),
especially equations (15)-(19). Their partially transmitting mirror
model is not the charged graph in this note. The proof above derives
the needed uniform capacity and pressure bounds directly. This is a
scoped technical extension supporting the provisional research route,
not an accepted consequential novelty claim. The [source supplement](../evidence/interacting-mouth-pressure-source-review.json)
records the exact reading scope and searches.

Continue physical long-return damping as provisional primary and
matched finite-duration transfer as fallback. No Lean or cloud builds
were run; formalization remains subordinate to closing the full
ordinary argument and its physical consequence.

Remaining to-do list: control the renormalized interacting ultraviolet
response and its local terms; close constraint-compatible geometry,
coupled state and exterior; establish the physical endpoint and novelty;
independently verify, GCP-certify suitable mathematics and prepare the
reproducible manuscript.
