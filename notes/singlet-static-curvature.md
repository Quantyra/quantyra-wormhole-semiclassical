# Renormalized vacuum curvature of the Gaussian singlet

2026-10-11, S053. The Gaussian spherical-photon model now has a bound on
its vacuum energy and first two clock derivatives, with the local
renormalization terms retained. In the working weak-coupling family this
closes the static premise of the [complex-response theorem](singlet-complex-response.md).
Together those results transfer the linear decay interval to the specified
canonical Gaussian model. Nonspherical interactions and the physical
four-dimensional transfer remain unresolved; this is not acceptance of
physical damping or novelty.

## 1. Model, normalization and statement

Use the [common-throat graph](schwinger-return-feedback.md): `q` circles
with a common throat coordinate, throat optical length `T(phi)`, and
fixed positive exterior lengths `D_j`. The scalar fields are separate
channel fields; only their unit singlet `v=(1,...,1)/sqrt(q)` receives a
nonnegative throat potential `m_phi(s)^2`. Each finite-model profile is
bounded. Charge sectors are fixed, exact zero oscillators contribute zero,
and a mode lifted to positive frequency is included in its oscillator
vacuum. This is the same convention as the previous response calculation;
it does not resolve compact-boson holonomy or fermion-state matching.

Set a fixed reference length and frequency to one. On `y=s/T` in `[0,1]`,
the positive Euclidean free form at `k>0` is

`A_k[f]=sum_j {T^-1 int_0^1 |f_j'|^2 dy`
`          + k^2 T int_0^1 |f_j|^2 dy`
`          + int_ext (|f_j'|^2+k^2|f_j|^2) ds}`.

The potential form is `V[f]=int_0^1 T m_phi(Ty)^2 |v*f(y)|^2 dy`.
The common form domain has the usual circle continuity at the two joins.
Let `D=partial_phi`. On a fixed clock interval assume

`0<T_min<=T<=T_max`, `-c_j A_k<=D^j A_k<=c_j A_k`,
`-c_j V<=D^j V<=c_j V`, `j=1,2`.                           (1)

The constants are uniform in `q,k,D_j` and the model sequence. Derivatives
and inequalities are in the form sense, with twice continuously varying
coefficients. Define

`I=int_0^T m_phi(s)^2 ds`, `Lbar=q^-1 sum_j (T+D_j)`.

For a fixed positive subtraction scale `mu` and fixed finite coefficient
`c`, the optical prescription is

`E_mu,c=(1/(2pi)) int_0^infinity [log det(1+B_k)`
`                         -I/(2sqrt(k^2+mu^2))] dk + c I`, (2)

where `B_k=V^(1/2) A_k^(-1) V^(1/2)` is understood on the single throat
space via the potential embedding. It is positive and trace class. The
determinant is the relative determinant of `A_k+V` and `A_k`, so varying
the kinetic measure does not introduce an uncancelled determinant.

There are `I_0>0,C<infinity`, depending only on the constants in (1), such
that, whenever `I<=I_0`, for `j=0,1,2`,

`|D^j E_mu,c| <= C [sqrt(I)+I Lbar+I^2`
`                              +I(1+|log mu|+|c|)]`.       (3)

The assumptions hold throughout the clock interval, so differentiation
under the renormalized integral is justified locally by the same bound.
Equation (3) is uniform in the number and individual lengths of returns.
It requires small integrated potential, not a uniform free spectral gap.

## 2. The low-frequency zero mode is retained

For points in the throat, the singlet-restricted free resolvent kernel is

`Gbar_k(s,t)=q^-1 sum_j [exp(-k|s-t|)+exp(-k(L_j-|s-t|))]`
`                             /[2k(1-exp(-k L_j))]`,
`L_j=T+D_j`.

Thus `Tr B_k=I avg_j coth(k L_j/2)/(2k)`. The zero-frequency part is
the constant kernel `avg_j(1/L_j)/k^2`. After multiplication by the mass
it has rank one and eigenvalue `a/k^2`, where
`a=I avg_j(1/L_j)<=I/T_min`. The remaining resolvent is positive and

`Tr B_k,reg <= I Lbar/12`.                                (4)

For example, expand each circle in Fourier modes. The nonzero modes at
`k=0` sum to `L_j/12`; increasing `k^2` decreases that positive operator.
Equivalently `coth u<=1/u+u/3` follows from the partial-fraction expansion.
For positive trace-class `X,Y`,
`log det(1+X+Y)<=log det(1+X)+Tr Y`, obtained by integrating the
resolvent along `X+tY`. Consequently, writing `L_k=log det(1+B_k)`,

`0<=L_k<=log(1+a/k^2)+I Lbar/12`,
`int_0^1 L_k dk<=pi sqrt(a)+I Lbar/12`.                    (5)

Here `int_0^infinity log(1+a/k^2) dk=pi sqrt(a)`.
The `sqrt(I)` contribution is the lifted zero oscillator, not a factor
of `q` zero oscillators. Subtracting the entire first Born term at zero
frequency would give an artificial infrared divergence. Equation (2)
subtracts a local ultraviolet term with an infrared-regular reference.

## 3. Two derivatives without a spectral gap

Write `R=A_k^-1`, `S=(A_k+V)^-1` and

`Z=Tr(S V)=Tr[(R-S) A_k]=Tr[B_k/(1+B_k)]<=L_k`.

All traces here are either positive finite traces or trace-class
differences. At a finite Galerkin cutoff, direct differentiation gives

`D L_k=Tr(S V')-Tr[(R-S) A_k']`,

`D^2 L_k=Tr(S V'')-Tr[(R-S) A_k'']`
` +Tr[R A_k' R A_k'-S A_k' S A_k']`
` -2Tr(S A_k' S V')-Tr(S V' S V')`.                       (6)

By (1), the first line of the second identity has magnitude at most
`2 c_2 Z`. Put `Delta=R-S>=0`. The difference in its next line equals
`Tr(Delta A_k' R A_k')+Tr(Delta A_k' S A_k')`. Both terms are
nonnegative and at most `c_1^2 Z`, since
`A_k' S A_k'<=A_k' R A_k'<=c_1^2 A_k`. For the cross term,
`||S^(1/2) A_k' S^(1/2)||<=c_1` and
`||S^(1/2) V' S^(1/2)||_1<=c_1 Z`. Finally
`P=S^(1/2)V S^(1/2)<=1`; the relative bound on `V'` gives
`Tr(S V' S V')<=c_1^2 Z`. These yield

`|D L_k|<=2c_1 Z`,
`|D^2 L_k|<=(2c_2+5c_1^2) Z`.                             (7)

The bounds pass to the form limit; locally the differentiated resolvents
converge with the trace factors supplied by `V`. Together with (5) they
control both infrared derivative integrals uniformly. This argument does
not assume that the derivative of the Birman-Schwinger operator is
bounded above and below by that operator. Such a claim need not hold
when its positive square-root factors vary noncommutatively.

## 4. Ultraviolet subtraction and its derivatives

For `k>=1`, the exact trace gives `b_k=Tr B_k<=C I/k`, uniformly in
the returns. Choose `I_0` so `b_k<=1/2`. Then

`|L_k-Tr B_k|<=b_k^2/2`.

The differentiated version, for `j<=2`, is

`|D^j(L_k-Tr B_k)|<=C_j b_k^2`.                           (8)

One direct proof expands `L_k-Tr B_k` in the absolutely convergent
series `sum_(n>=2) (-1)^(n+1) Tr[(R V)^n]/n`. Each derivative of `R`
is `-R A_k' R`, and
`R''=2R A_k' R A_k' R-R A_k'' R`. Insert the relative bounded forms
from (1). Schatten Holder applied to the `2n` factors
`R^(1/2)V^(1/2)` and their adjoints bounds each resulting trace by
`C_j n^j b_k^n`; the number of placements contributes at most another
quadratic polynomial in `n` for two derivatives. Summing that polynomial
against `b_k^n`, `b_k<=1/2`, proves (8). This proof also supplies
trace-norm convergence of the differentiated series.

Moreover, `Tr B_k-I/(2k)=I avg_j[1/(k(exp(kL_j)-1))]`.
Its first two clock derivatives are bounded by an integrable exponential
times `C I(1+k^2)`, since `D_j` is fixed and `T,T',T''` are bounded.
The remaining subtraction is local:

`int_1^infinity [1/k-1/sqrt(k^2+mu^2)] dk`
`                       =log[(1+sqrt(1+mu^2))/2]`,
`int_0^1 dk/sqrt(k^2+mu^2)=asinh(1/mu)`.

Both are bounded by `C(1+|log mu|)`. Equations (5), (7) and (8), with
`|D^j I|<=c_j I`, prove (3). The coefficient `c` stays fixed through
clock differentiation; it is not fitted to remove the curvature.

The reference-scale change is exactly

`E_mu2,c-E_mu1,c=I log(mu2/mu1)/(4pi)`.                    (9)

A physical change of prescription compensates this by changing the
specified coefficient `c`. This is not an arbitrary nonlocal subtraction
or a free adjustment of the clock stiffness.

## 5. Covariant prescription and the proper-cutoff profile

On a static two-dimensional metric `g=N(s)^2(-dt^2+ds^2)`, write
`m_opt^2=N^2 m_cov^2`. At separated points the covariant scalar Green
function equals the flat optical Green function with this potential.
The logarithmic coincident subtraction uses physical distance
`N|Delta s|`, so

`G_cov,ren-G_opt,ren=log N/(2pi)`.

Integrating `delta E/delta m_opt^2=G_ren/2` from zero potential gives
the mass-dependent local difference

`E_cov-E_opt=(1/(4pi)) int m_opt^2 log N ds`,               (10)

up to the already specified finite multiple of `I`. The mass-independent
Weyl anomaly cancels between the massive and massless models on the
same metric. Equation (10) is a local subtraction comparison; it does
not import a compact-surface determinant or a four-dimensional state.
For matrix potential `m_opt^2 vv*`, its channel trace is one times
`m_opt^2`, so no additional `q` appears. If
`|D log N|+|D^2 log N|<=C` and `||log N||_infinity<=H`, (10) and
its first two derivatives cost at most `C I(1+H)`.

For the truncated AdS throat, set `x=exp(phi)`, `z=epsilon x<=1/4`,
`b=sqrt(1-z^2)`, `a=acos z`, `rho=2y-1`. The spherical profile is

`T=2a/x`, `m_opt^2=(q g^2/(4pi^2)) x^2 sec^2(a rho)`,
`N=r_e x sec(a rho)`.

Its integrated potential and derivatives are exactly

`I=q g^2 b/(2pi^2 epsilon)`,
`D I=-q g^2 epsilon x^2/(2pi^2 b)`,
`D^2 I=-q g^2 epsilon [2x^2/b+epsilon^2 x^4/b^3]/(2pi^2)`.

Since `D a=-z/b`, `D^2 a=-z/b^3`,
`D log(T m_opt^2)=1-z/(ab)-2rho z tan(a rho)/b`.
Using `|tan(a rho)|<=b/z` and `sec^2(a rho)<=1/z^2`, its absolute
value is at most 3 and its second logarithmic derivative at most 5.
Thus `|D^2(T m_opt^2)|<=14 T m_opt^2`. The coefficients `T` and
`1/T` obey smaller relative bounds; (1) holds with `c_1=3,c_2=16`.
The first logarithmic derivative of `N` lies in `[0,1]`, and the second
has magnitude at most `2/b^2<3`. Hence (3) and (10) apply even with
these moving proper cutoffs. Each cutoff is finite; no singular-endpoint
determinant is asserted.

The canonical family used for the response theorem is simpler:
`T=T_0 exp(-phi)`, `m_phi(Ty)^2=exp(2phi)m_0(T_0 y)^2`.
Then `V_phi=exp(phi)V_0` in the fixed coordinate, precisely the
Hamiltonian vertex `H_phi=exp(phi)H_throat+H_ext`. It also satisfies
(1). A truncated spherical reference profile, rescaled to `T_0=pi`,
has `I=O(q g^2/epsilon)` uniformly on a fixed clock interval.
Its local covariant term can be included with `N_phi=exp(phi)N_0`.
The proper-cutoff family has an additional changing shape and vertex;
its static estimate alone does not identify its dynamical response with
this canonical one.

## 6. Consequence in the existing weak-coupling family

Take the already declared `g^2=q^-2`, `epsilon=q^-1/4`,
`delta=q^-1/16`, with all lengths in reference throat units. For the
inverse-tail returns `D_j=kappa delta/u_j`, `u_j=(j-1/2)/q`,
`Lbar=T+O(delta log q)=O(1)` and `I=O(q^-3/4)`.
The actual dipole upper bound `D_j<=sqrt(5) delta/u_j` gives the same
mean-length estimate, although its full dynamical transfer remains open.
For polynomial reference scales and lapses, and a fixed finite coefficient
`c` (or one growing at most logarithmically with `q`), the total static
energy correction and both derivatives obey

`R_q=O(q^-3/8+q^-3/4 log q)=O(q^-3/8)`.

Here the physical reference lapse obeys `r_e/ell_0=O(q^-3)`, so
`H=O(log q)` also at the proper cutoff. Renormalization coefficients
from a yet-unperformed four-dimensional matching are not assumed fixed
by this Gaussian calculation.

Near the free equilibrium the free total energy has second derivative
bounded below by `c q` and third derivative `O(q)`. The bound on the
first two derivatives of the correction therefore gives a unique nearby
minimum, displaced by `O(R_q/q)=O(q^-11/8)`, and changes its static
curvature by `O(R_q)`. This follows from monotonicity of the first
derivative in that neighborhood and the intermediate value theorem.
At the same reference clock,

`(K0_m-K0_free)/(q delta)=O(q^-21/16)`.

Recentring the free response at the displaced throat changes the kernel
by `O(q |Delta T| eta^-2)`, using the length-derivative estimate in
the [finite-channel proof](finite-channel-decay.md). At `eta~delta`,
its relative size is `O(R_q/(q delta^3))=O(q^-19/16)`.
The previous singlet dynamical comparison is `O(1/(q delta^5))`,
or `O(q^-11/16)`, relative to `q delta`. All these errors vanish.

Consequently the canonical Gaussian model, with (2) and its specified
local term, has the same scaled-kernel limit and finite-time linear decay
envelope as the free inverse-tail theorem. Its limiting width is
`gamma delta`, `gamma=4 kappa(8log 2-13pi^2/24)/(3pi^2)`.
The positive spectral realization applies because the resulting static
curvature stays positive and the massive subtracted response has the
summability already proved. Every finite system remains recurrent.
This is a conclusion about the derived Gaussian model, with its stated
zero-sector and vertex conventions, rather than a choice of a stabilizing
counterterm. It supplies no finite-disturbance survival threshold.

## 7. Verification, sources and decision

Run `python checks/singlet_static_curvature.py`; the
[numerical receipt](../evidence/singlet-static-curvature.json) and
[source audit](../evidence/singlet-static-source-review.json) retain the
checks and their scope. Independent periodic transfer matrices check the
Fredholm determinant and energy calculation; direct finite matrices check
both noncommuting derivative identities. The ordinary estimates above
establish uniformity, not the finite numerical sample.

[Graham et al., hep-th/0207120v2](https://arxiv.org/abs/hep-th/0207120v2),
sections 2.3-2.4, 3.2 and Appendix A, supply the prior framework for
imaginary-frequency vacuum energies and Born subtraction with fixed
counterterms. Their infinite-line hypotheses are not a substitute for
the compact graph and zero-mode proof here. [Bilal and Leduc,
1606.01901v2](https://arxiv.org/abs/1606.01901v2), equations (2.10)-(2.11)
and (3.11)-(3.14), fix the coincident subtraction and positive sign of
the local logarithm in (10). Their full compact-surface action is not
used as a Lorentzian wormhole action. These established techniques are
prior art; this intermediate estimate is not accepted as standalone
novelty or an expanded research target.

Continue the physical long-return damping primary and finite-duration
energy-transfer fallback. The Gaussian static and dynamical singlet
gaps are now controlled in the declared weak-coupling family. The
nonspherical current interactions, corrected physical background,
complete return map, collar response and four-dimensional state still
need quantitative bounds. No Lean, Lake, GCP or publication action occurred.

Remaining to-do list: control nonspherical interactions and physical
background, return, vertex and state matching; finish novelty and the
consequential physical result; independently verify, GCP-certify suitable
endpoints and prepare the reproducible manuscript.
