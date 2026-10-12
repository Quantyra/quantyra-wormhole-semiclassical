# JT mass Casimir and the RN boundary charge

2026-10-12. The RN mass deficit first appears at quadratic order in the
near-product throat dilaton. Expanding the exact charged Misner-Sharp mass
gives a gauge-invariant JT Casimir whose source evolution is fixed by the
traceless stress. This reconciles the linear sourced JT constraints with the
quadratic static RN matching and identifies which order a stability
calculation must retain. The result is a reproduction of the standard JT
mass Casimir and MMP's static match, not a new mass theorem.

## Exact mass expansion

For fixed magnetic charge `r_e`, the spherical mass aspect is

```text
m=R/(2G)[1-(nabla R)^2+r_e^2/R^2],
R=r_e sqrt(1+Phi),  g_ab=r_e^2 gbar_ab.
```

Thus, using the unit-radius metric `gbar`,

```text
(nabla R)^2=(nabla_bar Phi)^2/[4(1+Phi)],
m=r_e/(2G)[sqrt(1+Phi)+1/sqrt(1+Phi)
            -(nabla_bar Phi)^2/(4 sqrt(1+Phi))].             (1)
```

Expanding at small `Phi` and its first derivatives gives

```text
m=r_e/G + r_e/(8G) C + O(Phi^3),
C=Phi^2-(nabla_bar Phi)^2.                                  (2)
```

There is no term linear in `Phi`: at the product throat, changing the
area dilaton infinitesimally does not change the RN mass. Mass matching
through first order in the JT equation is therefore insufficient. One must
retain the quadratic mass Casimir and match perturbations around the
corrected static background.

## Sourced Casimir evolution

Write the leading JT constraint in unit-radius `AdS2` as

```text
L_ab[Phi]=gbar_ab Box Phi-nabla_a nabla_b Phi-gbar_ab Phi=S_ab,
S=gbar^ab S_ab.
```

Taking the trace gives `Box Phi-2Phi=S`, so the Hessian is
`nabla_a nabla_b Phi=gbar_ab(Phi+S)-S_ab`. Differentiating `C` in (2)
then yields

```text
nabla_a C=2(S_ab-gbar_ab S)nabla^b Phi.                    (3)
```

For the state-dependent LLL source, `S=0`; its stress transfers the JT
Casimir according to `nabla_a C=2 S_ab nabla^b Phi`. In a source-free
overlap, `C` is constant. This is the JT form of the RN Birkhoff mass
constant. With a source, (3) is the local mass-work identity that must
agree with the full spherical mass-flux balance when the reduction
normalizations and background corrections are included.

At the leading product-throat order the normalization agreement can be made
explicit. From the exact reduction,
`S_ab=(2G/r_e^2) delta T_ab^(2)` and
`nabla_b R=(r_e/2) partial_b Phi`, while raising an orbit-space index in
the exact mass balance contributes `r_e^-2`. Therefore

```text
nabla_a m = r_e/(4G) (S_ab-gbar_ab S) nabla_bar^b Phi
           = 1/(2r_e) (delta T_ab^(2)-gbar_ab gbar^cd delta T_cd^(2))
                         nabla_bar^b Phi + higher throat orders.
```

This is precisely the small-dilaton expansion of
`nabla_a m=(T_a^b-delta_a^b T)nabla_b R`. It closes the local normalization
bridge between the sourced JT Casimir flux and the exact spherical mass
balance. It does not determine the retarded stress, moving-boundary contact,
or integration constants at either asymptotic mouth; those still require the
corrected matched solution.

The three homogeneous dilaton modes can be written

```text
Phi_h=a rho + sqrt(1+rho^2)[b cos(tau)+c sin(tau)].        (4)
```

Their Casimir is `C_h=b^2+c^2-a^2`. Thus the three source-free embedding
modes are not all mass-neutral; their invariant combination determines
the asymptotic mass data. A coordinate/gauge classification of each mode
is still separate from this Casimir statement.

## Reproduction of the MMP static mouth mass

The large-`rho` MMP throat profile has
`Phi ~= beta rho`, `beta=2r_e/ell`, from its static solution and overlap
matching. In a source-free RN overlap the leading homogeneous mode has
`C=-beta^2`, so (2) gives

```text
deltaM=m-r_e/G=-r_e beta^2/(8G)
                     =-r_e^3/(2G ell^2)=-A x_*^2/2,
A=r_e^3/G,  x_*=1/ell.                                    (5)
```

This reproduces MMP's per-mouth mass shift in
[equations (5.39)-(5.42)](https://arxiv.org/html/1807.04726v3#S5.SS3)
and the boundary relation in `rn-jt-static-boundary-charge-match.md`.
The derivation clarifies why the shift is quadratic in the leading throat
flare even though its static value is fixed by the first-order LLL stress.

For a linear stability perturbation around the supported background,
`delta C=2(Phi_0 deltaPhi-gbar^ab partial_a Phi_0 partial_b deltaPhi)`
at first order in the perturbation, plus variations of the background
metric and source. Dropping `Phi_0` and using only the product-background
linear equation would incorrectly set this RN mass response to zero. The
quadratic completion, corrected static metric, and state-dependent stress
must be treated together.

## Verification and scope

Run `checks/jt_mass_casimir.ps1`. It checks the local identity (3) on
deterministic Lorentzian jets, the leading normalization bridge from the
JT Casimir flux to the spherical Misner-Sharp balance, the exact expansion
(1)-(2), and the MMP normalization in (5). These checks verify the algebra
and perturbative order; they do not supply the corrected two-ended background
or finite-frequency boundary response.

Remaining to-do list: evaluate the Casimir variation on the corrected
MMP background; equate its source flux to the RN collar mass balance;
include the matched stress/current and join contacts; then solve the
physical retarded stability problem.
