# Actual energy-domain transfer and the second-variation chart contact

2026-10-21, E003/S168. Ordinary argument before focused verification.
The complete stable-wormhole goal remains incomplete.

Two obligations left by the [normal moving-cut theorem](finite-cut-normal-shape-limit.md)
can be resolved without assuming pointwise higher boundary jets.
The full auxiliary shape energy controls that theorem's compensated
logarithmic domain on a fixed outer collar. Also, the actual fixed-seed
first variation applied to the squared shape has a finite weak cutoff
limit. This is precisely the contact relating exponential and affine
mass parametrizations of the complete renormalized functional.

Neither result computes the actual nonlocal finite-cut Hessian. The
remaining comparison must retain the two-vertex response and its fixed
local allocation, curved boundary kernels, thermal state and remote
geometry. No physical fast-branch exclusion or stability follows.

## 1. The actual closed energy controls the localized compensated form

Fix a positive Matsubara frequency above the retained global threshold.
Use the same charged seed, HH state, paired domain, full angular tower,
physical clock and subtraction. Let a_gamma be the full closed auxiliary
form, with v=w/r and s=Q_gamma w/r its curvature-free reconstruction.
Its [retained bounds](global-thermal-shape-inverse.md) give

    N_gamma<=33G a_gamma,
    C J2[v]<=200pi a_gamma,
    C J2[s+v/3]<=40pi a_gamma.                         (1)

Begin on its compact w core. Define

    F(v)=Qcan[v]-C/(30pi) integral(omega+19/60)(Delta v)^2,
    omega=log(r/B), c=C/(30pi).

All integrals in this definition use the retained normalized canonical
orbit metric and actual time-mode norm. The exact Weyl completion,
classical lower bound and residual budget in
[the global thermal proof](global-thermal-curvature-free-form.md),
equations (2),(9),(10), imply F(v)<=a_gamma[w]. This is an upper bound
on the compensated form extracted from the COMPLETE action; it does
not assume either logarithmic summand is separately finite in the closure.

Use the fixed partition chi_o^2+chi_i^2=1 from that proof. Its quantum
double commutator and weighted fourth-order localization are bounded
by K(J2+J1+J0). The same proof bounds F(chi_i v) below by a finite
negative multiple of J1+J0. Hence

    F(chi_o v)<=K_gamma a_gamma[w].                    (2)

The constants are finite for this seed and fixed frequency. J1 and J0
are controlled by J2 at positive frequency. Localization of J2 adds
only fixed lower-order terms, so (1) also bounds all the localized
H2 norms below. No numerical value of K_gamma or new threshold is claimed.

## 2. A two-sided comparison with the Mellin logarithmic domain

Use inward canonical proper distance Y, and rescale time by the fixed
outer lapse a0. On the same complete even auxiliary cylinder used in
the retained scalar comparison, write a(Y)=a0 exp(h(Y)), with h in
C^(2,1), bounded with its first three weak derivatives, and h(0)=h'(0)=0.
Put nu=gamma/a0>0 and, for T=chi_o v,

    u=sqrt(a/a0) T,
    H=-partial_Y^2+h''/2+h'^2/4+nu^2 exp(-2h),
    P=-partial_Y^2+nu^2, H=P+W.

W and W' are bounded. The factorization of H at nu=0 gives
H>=exp(-2||h||_infinity)nu^2; P>=nu^2. If m is the smaller positive
gap, the resolvent identity gives

    log H-log P=integral_0^infinity
                 (H+t)^-1 W(P+t)^-1 dt,
    ||log H-log P||<=||W||_infinity/m.                 (3)

This follows by an absolutely convergent operator-norm integral. It
does not require H and P to commute. Because multiplication by W maps
H1 to H1 and log sqrt(P) maps H1 to L2, expanding Hu=Pu+Wu proves

    |<Hu,log sqrt(H)Hu>-<Pu,log sqrt(P)Pu>|
                                                <=K_gamma ||u||_H2^2. (4)

For example the two cross terms use <Pu,log sqrt(P)Wu>; no logarithmic
norm of Pu is estimated by H2. The low-momentum difference between
p^4 log sqrt(p^2+nu^2) and p^4 log|p| is bounded by a constant times
1+p^2, while the high-momentum difference is O(nu^2 p^2).

Use the improved full-tower remainder from the preceding normal theorem:
it is O((1+p^2)log(e+p^2)), with the original fourth coefficient.
Its functional calculus for H is bounded as a quadratic form by a
finite multiple of ||u||_H2^2. Also omega=-log Y+constant+O(Y^2).
Replacing Hu by -u'' in the weighted local term costs only H2:
the cross terms use ||u''||_2 ||log(Y)Wu||_2, and u(0)=0 bounds
||log(Y)u||_2 by its H1 norm on this fixed collar. All remaining
constant contacts are retained in the same H2 error budget.

The actual/geodesic and geodesic/scalar TWO-SIDED comparisons from
[the uniform outer proof](uniform-logarithmic-frequency-comparison.md)
now yield

    |F(T)-a0 c P_normal(u)|<=K_gamma ||T||_H2^2,        (5)

where a0 retains the original time-mode normalization. P_normal is
the zero-cut compensated even form of the previous theorem. Its Mellin
symbol is Re psi(1/2+i xi)-pi/[2cosh(pi xi)]. The second term is bounded.
Consequently its form domain is exactly the logarithmic domain with
seminorm

    E_D(g)=(1/4) integral_R^2 k(|z-z'|)|g(z)-g(z')|^2 dz dz',
    k(r)=1/[2sinh(r/2)], z=-log Y, g=sqrt(Y)u''(Y).      (6)

Here the digamma identity is the same established
[DLMF integral](https://dlmf.nist.gov/5.9.E16). We are proving transfer
of the actual energy norm, not claiming a new Mellin formula.

The smooth density factor does not change the domain. Indeed if
m0=sqrt(a/a0), then

    sqrt(Y)u''=m0 g_T+ell,
    g_T=sqrt(Y)T'', ell=sqrt(Y)(2m0'T'+m0''T).

ell belongs to H1(dz), with norm bounded by ||T||_H2. Both m0 and
its inverse are bounded and Lipschitz as functions of z. Directly from
(6), for any bounded Lipschitz multiplier b,

    E_D(bg)<=2||b||_infinity^2 E_D(g)
                         +14 zeta(3) Lip(b)^2 ||g||_2^2,
    E_D(ell)<=7 zeta(3)||ell'||_2^2.                  (7)

The constants follow from integral_0^infinity r^2 k(r)dr=14zeta(3).
Apply (7) also to 1/m0. Combining (1),(2),(5)-(7) proves on the core

    ||g_T||_2^2+E_D(g_T)<=K_gamma a_gamma[w].           (8)

An a_gamma-Cauchy core sequence is therefore Cauchy in this logarithmic
space. Its L2 limit is the weak derivative g_T already determined by
J2 convergence. The closed Fourier form in (6) identifies the limits.
Thus (8) holds for the FULL auxiliary energy domain, not only the
earlier stronger model-plus-absolute-Weyl completion. The preceding
normal weak-cutoff theorem now applies to every such localized shape.
This is a continuous inclusion; no equality with a gravitational
physical domain or with a global logarithmic Sobolev space is claimed.

## 3. A weak moving-cut lemma for cubic tests

The remaining contact uses a fixed-seed FIRST variation. The earlier
[cubic-direction argument](coupled-fast-quadratic-boundary-obstruction.md)
provides its actual heat-kernel decomposition on a collar at cut d:
finite sums of fixed-scale finite parts FP[x^-p log^j x], p<=3,
contacts of order at most two, and uniformly L1 remainders. Their
coefficients and rescaled d derivatives are bounded and converge.
They are fixed-seed jet expressions. Off the collar the kernels converge
smoothly after the common subtraction; the retained gap controls large
heat time. The balanced paired spin trace has no weight-zero image
layer, and tr chi=0 removes the possible third smearing derivative in
[Vassilevich (5.33)](https://arxiv.org/html/hep-th/0306138v3).

Let U(Y) be continuous, smooth for Y>0 initially, with |U(Y)|<=M Y^3
on a fixed collar. No bound on U''' or a pointwise limit of U'' is
required. For every smooth compact phi, set

    psi_delta(d)=phi(log(delta/d))/d, delta=exp(-Z0).

Its support lies in [a delta,b delta] for fixed positive a,b. Averaging
the cut pairing against phi(Z-Z0), Z=-log d, is its pairing against
psi_delta(d) dd. For a finite part T_p, move the distribution onto this
smooth test before estimating:

    integral psi_delta(d)<T_p(x),U(d+x)> dd
                         =<U(Y),(T_p*psi_delta)(Y)>.    (9)

On Y<=2b delta, homogeneity and the finite-part scaling contact bound
the right kernel by C_phi delta^-p(1+|log delta|)^m, with a fixed finite
m also covering logarithmic partners. Multiplication by U=O(Y^3) and
integration over this interval therefore cost

    C_phi M delta^(4-p)(1+|log delta|)^m ->0.           (10)

On Y>2b delta, the kernel is ordinary. It tends pointwise to
(integral phi)Y^-p log^j Y and is dominated by a fixed multiple of
Y^-p(1+|log Y|)^j. This is integrable against U near zero. Dominated
convergence supplies its finite limit. Smooth cutoff factors in the
heat decomposition and convergent coefficient jets obey the same
argument. No cut-dependent subtraction length is introduced.

A contact U^(q)(d), q<=2, is handled by q integrations by parts onto
psi_delta times its bounded smooth coefficient. Its bound is
C delta^(3-q), with possible fixed logarithmic factors, and tends to
zero. Uniformly L1 remainder kernels converge against U(d+.) by first
removing a small collar, then using |U|<=M Y^3 there. Their smooth
off-collar convergence suffices; no differentiated remainder is assumed.
This proves a finite translated-logarithmic weak limit of the complete
fixed-seed first variation on U. It is not a pointwise cutoff theorem.

For the shape, the retained energy estimates give
v=Y^(3/2)V(-log Y) with V in H2 on the outer half-line. Thus V is
bounded and U=v^2 satisfies the cubic envelope. Polarization also
covers v1 v2. The weak limit is the original limiting first-variation
pairing, with its finite interior/state value retained. At the geodesic
limit all derivative-of-smearing boundary terms vanish, while the
remaining value contact pairs to U(0)=0. Off the outer collar, the
fixed-state first-variation kernels are smooth and convergent; the
retained compactified-orbit bounds give the ordinary contribution.

The absence of a pointwise hypothesis is substantive. Take a smooth
base v0(Y)=Y^2 near zero and disjoint bumps
b_n=d_n^13 psi((Y-d_n)/d_n^8), d_n=2^-n, with
psi(0)=psi'(0)=0, psi''(0)=1. Their H2 seminorm costs are O(d_n^2).
The Mellin derivative square costs O(d_n^-12), so the logarithmic
energy costs O(d_n^2(1+|log d_n|)) by Jensen's inequality for log.
The sum of square roots of these costs is finite. Nevertheless, for
v=v0+sum b_n, (v^2)''(d_n)=12d_n^2+2/d_n diverges. A compact smooth
bump suffices; the checker uses a C5 polynomial bump with the same
center jets. This is a regularity example, not an actual shape solve.

## 4. The actual second-variation chart contact

At every positive cut the original complete canonical functional is
Gamma_d[V], including its two-dimensional measure restoration and the
common four-dimensional subtraction. Compare two paths only as a
chain-rule identity:

    V_exp=V0+eta v H(tau),
    V_aff=V0+log(1+eta v H(tau)).

The latter makes the optical mass affine. Keeping the original
exponential path as the physical parametrization, exact differentiation
gives

    Gamma_d,exp''[v,v]-Gamma_d,aff''[v,v]
                                      =Gamma_d'[v^2 H^2]. (11)

The stationary seed sees only the mean of H^2, with its actual time
normalization. Section 3 proves that the right side has a finite weak
cutoff limit for the entire auxiliary energy class. It cannot supply
a persistent nonconstant cutoff oscillation. Its finite value is
retained, not set to zero.

In the optical mass representation, the nonlocal determinant part of
(11) is -Tr(S M v^2), alongside the two-vertex term Tr(S Mv S Mv).
The measure/restoration pieces must follow the same chain rule. The
statement concerns their COMPLETE allocation, not an independently
renormalized bare tadpole. The actual constant row also contains the
two-vertex row sum; (11) alone does not evaluate that row or the whole
finite-cut response. No alternative physical mass chart is adopted.

The proof uses smooth finite-cut paths and extends their quadratic
contact by the weak estimates. It is not a C2 theorem for finite
amplitude rough families at the limiting boundary. Nor does it
identify the affine Hessian with the retained geodesic model.

## 5. Remaining actual comparison

The normal theorem now applies to the actual auxiliary shape energy
domain, and the chart contact has a finite necessary cutoff limit.
What remains is the actual nonlocal finite-cut comparison after this
contact is allocated, including the two-vertex kernel and its fixed
diagonal terms, charged curvature, boundary geometry, thermal state
and remote attachment. The original local shape partners must remain
in the same cutoff action. The limiting geodesic-seam estimate alone
still does not bound that moving-domain difference.

Verification checks the noncommuting logarithm-resolvent identity,
Mellin localization, finite-part convolution and scale contacts,
concentrated squared-shape jets and the complete finite-regulator
parametrization identity. These support the ordinary transfer proofs;
they do not evaluate the actual charged Dirac Hessian or a physical
mode. No Lean run is selected. S168.3-S168.4 remain unchecked.

Remaining to-do list: actual nonlocal finite-cut comparison with all
local shape partners; physical boundary phase and source compatibility;
continuous causal constraints, exact nonlinear existence and all-sector
stability; exterior/state, charge/return, opening, transport, preparation,
separate creation and full audit.
