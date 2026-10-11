# Axial Casimir response with the Dirac curl vertex

2026-10-10, Space S051/S052. The magnetic-cylinder calculation now includes
both longitudinal components of the axial metric and gauge sources and their
Euclidean frequency dependence. The mixed lowest/higher-Landau-level response
is exact. The remaining massive/massive response and contact term have an
analytic all-level bound. A chirality-changing curl vertex, absent in the
static spatial test, reduces to a transverse derivative shift and a massive
tadpole. This supplies a calculable response beyond the prescribed drift.

The geometry remains a prescribed product cylinder, with one free Dirac
field. This is neither a closed semiclassical gravitational model nor a
physical stability result. In particular, these four source components do
not constitute a complete four-dimensional Ward-identity audit. The
[static sign theorem](casimir-cylinder-twist-response.md) is preserved and
reproduced as a limit.

## Observable and source path

Use the same magnetic flux, radius, circumference and antiperiodic spin
structure as the static note, with integer $q\ge3$. Work in Euclidean
coordinates $x=(\tau,z)$, with background

$$ds_E^2=d\tau^2+dz^2+R^2d\Omega^2,\qquad
A_\varphi=q\cos\theta/2.$$

The uncompactified Euclidean time implements zero physical temperature.
Put $u=\cos\theta$, $Y=P_2(u)$ and $Q=(\Omega,k)$, with real $\Omega$ and
$k=2\pi j/\mathcal L$, $j\ne0$. Let
$f(x)=\sqrt2\cos(\Omega\tau+kz)$ and choose real amplitude vectors $a,h$:

$$g_{a\varphi}(\epsilon)=-\epsilon h_a f(x)(1-u^2)Y'(u),\qquad
A_a(\epsilon)=\epsilon a_a f(x)Y(u),\qquad a\in\{\tau,z\}.$$

All other metric and gauge components are fixed through second order. The
source path is affine in the coordinate metric, including $g_{ab}$.
Subtract the infinite-line determinant with the same periodic sources,
regulator and local couplings. Define $\mathscr Q_E$ by the coefficient
of $\epsilon^2\mathcal T\mathcal L$ in the resulting effective action
$\Delta\Gamma_E=-\Delta\log\det D$, taking $\mathcal T\to\infty$.
The normalization of $f$ makes its spacetime mean square one.

This fixes the Euclidean observable. Lorentzian temporal sources require
their Wick-rotation factors; the Euclidean Hessian is not a gravitational
canonical energy or a stability matrix.

## Full linear vertex and angular transitions

The independent tetrad calculation now verifies, before projecting, the
full warped-background Lorentzian metric vertex

$$\delta D=-\frac iN\sum_{a=0,1}\gamma^a
\left[v_a(m-qu/2)I+\frac{\sin\theta\,\partial_\theta v_a}{4}\sigma_3\right]
+\frac{R\sin\theta}{4N^2}
(\partial_t v_z-\partial_z v_t)\gamma^0\gamma^1\gamma^3.$$

Here $v_a=-h_aY'/R^2$, and $\sigma_3$ acts on transverse chirality. Arbitrary
$N(t,z),R(t,z)$ are retained in this **linear operator** identity. The loop
calculation below specializes to constant radius and $N=1$; it does not
thereby derive a vacuum response on the warped background.

Such Pauli-type curl couplings are established dimensional-reduction
structure. Maraner and Pachos,
[0704.2076v2](https://arxiv.org/abs/0704.2076v2), section III.6.2, HTML
equation (58), explicitly retain an induced Pauli term in the reduced
Dirac operator. The new comparison target is its magnetic-cylinder vacuum
response, not the existence of that coupling.

In Euclidean signature the algebraic vertex is $-i\gamma^a V_a$, with

$$V_{a,\sigma}=a_aY-\frac{h_a}{R^2}T_\sigma,\qquad
T_\sigma=Y'(m-qu/2)-\frac\sigma4(1-u^2)Y''.$$

The additional vertex is
$-w(\partial_\tau h_z-\partial_z h_\tau)\gamma^\tau\gamma^z\gamma^3/(4R)$,
where $w=\sin\theta Y'$. If $D_S$ is the anti-Hermitian angular Dirac
operator, its commutator with $Y$ gives

$$w\sigma_2=-iR[D_S,Y]\sigma_3.$$

On a negative-chirality zero mode this is $iR D_SY$. Therefore the curl
has no zero-mode/zero-mode matrix element. It connects the zero mode only
to the same two angular levels as the centered quadratic harmonic. Choose
the massive chirality basis so that $D_S=iM_n\sigma_1$ in that block. Its
curl matrix element is then $-RM_n$ times the corresponding $Y$ overlap.
The known weights, masses and algebraic shifts are

$$S_1=\frac{6q(q-1)}{5(q+1)(q+3)},\quad
S_2=\frac{12q}{5(q+2)(q+3)},$$
$$M_1^2=\frac{q+1}{R^2},\quad M_2^2=\frac{2(q+2)}{R^2},\qquad
b_1=(q-2)/2,\quad b_2=(2q+1)/2.$$

For each level put $c=a+b_nh/R^2$ and $F=\Omega h_z-kh_\tau$.
The angular weights include the two massive angular eigenvalue signs;
they are not counted a second time.

## The pure LLL finite-size difference vanishes

The algebraic LLL projection is the two-dimensional vector potential
$\langle Y\rangle_m(a-3h/(2R^2))$. Its compact and infinite-line current
kernels agree at nonzero allowed $k$. One direct proof uses the chiral
occupation sum. For $k=2\pi j/\mathcal L>0$, right-moving energies
$r=2\pi(n+1/2)/\mathcal L$ give

$$\frac1{\mathcal L}\sum_n
\big[\mathbf1_{n<0}-\mathbf1_{n+j<0}\big]
=\frac j{\mathcal L}=\frac k{2\pi}.$$

Only $n=-j,\ldots,-1$ contribute. Their response denominator is the same
$i\Omega-k$ for every term. The infinite-line occupation integral gives
the identical $k/(2\pi)$ coefficient. The left-moving channel gives the
corresponding identity with $k\mapsto-k$. Common local regulator terms
cancel in the subtraction. The curl projection is zero, so the pure LLL
contribution to $\mathscr Q_E$ is exactly zero. This does not remove the
LLL current response from either unsubtracted theory, and it excludes the
global charge/holonomy mode $k=0$.

## Exact mixed response and curl reduction

Use two-dimensional Euclidean gamma matrices $(\sigma_1,\sigma_2)$ and
$\gamma_* = \sigma_3$. For loop momentum $P=(p,r)$, let $P'=P+Q$. In an
angular block ordered $(+,-)$, the massive propagator numerator is

$$-i\begin{pmatrix}\gamma\cdot P'&M\gamma_*\\
M\gamma_*&\gamma\cdot P'\end{pmatrix}.$$

The vertex row at $-Q$ and column at $Q$ are respectively
$[MF I/4,-i\gamma\cdot c]$ and $[-MF I/4,-i\gamma\cdot c]^T$.
These signs follow from $\partial_a\leftrightarrow iQ_a$ and the routing
$S_L(P)V(-Q)S_H(P+Q)V(Q)$. Tracing gives the real numerator

$$\mathcal N(P)=4(P\cdot c)(P'\cdot c)-2(P\cdot P')c^2
+M^2F(pc_z-rc_\tau)+\frac{M^2F^2}{8}P\cdot P'.$$

The denominator is $(r^2+p^2)[(r+k)^2+(p+\Omega)^2+M^2]$.
The mixed-sector factor $1/2$ in the quadratic determinant is canceled
by its two orientations. Thus no additional factor $1/2$ multiplies this
trace when weighted by $S_n$.

Define $\mathcal V_M(c;Q)$ as the compact-minus-line loop for this
numerator with $F=0$. Then the exact reduction is

$$\boxed{\mathcal M_M(c,h;Q)=
\mathcal V_M\!\left(c+\frac{Q^2h-Q(Q\cdot h)}4;Q\right)
-\frac{M^2F^2}{8}\Theta_M},$$

$$\Theta_M=\frac1\pi\int_0^\infty
\frac{n_F(\mathcal L\sqrt{p^2+M^2})}{\sqrt{p^2+M^2}}\,dp,
\qquad 0\le\Theta_M\le K_0(M\mathcal L)/\pi.$$

This is not a derivative expansion: the shift contains $Q^2$, while
$\mathcal V_M$ retains its full frequency and momentum dependence. The
massive tadpole is exponentially suppressed at long circumference.

For an explicit proof and numerical definition, expand any of the trace
numerators as $n_2r^2+n_1r+n_0$ and set

$$d=M^2+\Omega^2+k^2+2p\Omega,\qquad H=d^2+4k^2p^2,$$
$$A=\frac{dn_1-2k(n_0-n_2p^2)}H,\quad
B=\frac{d(n_0-n_2p^2)+2kp^2n_1}H,\quad C=n_2-kA-B.$$

The partial fractions are $(Ar+B)/(r^2+p^2)+
[-A(r+k)+C]/[(r+k)^2+(p+\Omega)^2+M^2]$. Since $k$ is a bosonic circle
momentum, both odd terms vanish in the sum and integral. The finite-size
difference is the convergent one-dimensional quadrature

$$\int_{-\infty}^{\infty}\frac{dp}{2\pi}
\left[-\frac{B\,n_F(\mathcal L|p|)}{|p|}
-\frac{C\,n_F(\mathcal L E)}E\right],\qquad
E=\sqrt{(p+\Omega)^2+M^2}.$$

The value at the isolated point $p=0$ is immaterial. In the difference
between the full numerator and its shifted-vector version, the residues are

$$B_{\rm diff}=\frac{F(F\Omega+8c_z)}8p,\qquad
C_{\rm diff}=\frac{M^2F^2}8-
\frac{F(F\Omega+8c_z)}8(p+\Omega).$$

The first integrates to zero by reflection in $p$; the second leaves only
the stated tadpole after reflection in $p+\Omega$. This proves the boxed
identity without an expansion or a Landau-level cutoff. The vector kernel
is a real symmetric, traceless two-by-two matrix, as follows directly by
summing the two vector-component trace numerators.

## Contact term and the controlled full-sector approximation

The affine metric path requires a quadratic coframe correction. The LLL
kinetic tadpole, with the same current/measure reasoning as in the static
note, gives

$$\mathscr Q_{E,\mathrm{contact},0}
=\frac{\pi W}{12R^2\mathcal L^2}(h_\tau^2-h_z^2),\qquad W=6q/5.$$

Its signs also follow by varying
$-\pi q\sqrt{g_{\tau\tau}}/(6\mathcal L\sqrt{g_{zz}})$ and canceling
the quadratic base metric induced by the first-order coframe. One-point
currents vanish in the prescribed zero-density vacuum. The residual
massive contact contribution is bounded below, not discarded.

The full four-source response therefore has the representation

$$\mathscr Q_E=
\frac{\pi W(h_\tau^2-h_z^2)}{12R^2\mathcal L^2}
+\sum_{n=1}^2S_n\mathcal M_{M_n}(a+b_nh/R^2,h;Q)
+\mathcal R_E.$$

For the entire remaining angular spectrum $M_n=\sqrt{n(n+q)}/R$, define

$$B_n=|a|+\frac{3(q+n)}{R^2}|h|+\frac{3|F|}{8R},$$

where vector norms are Euclidean. The last term uses
$|\sin\theta P'_2(\cos\theta)|\le3/2$. Then, for real $\Omega$ and $k\ne0$,

$$|\mathcal R_E|\le
\frac4{\pi|k|}\sum_{n\ge1}(q+2n)B_n^2M_nK_1(\mathcal LM_n)
+\frac{6|h|^2}{5\pi R^2}\sum_{n\ge1}(q+2n)M_n^2
[K_0(\mathcal LM_n)+K_2(\mathcal LM_n)].$$

Proof: rotate compact $z$ to Euclidean time. The external $\Omega$ becomes
a real spatial momentum transfer, whereas the external frequency $k$ is
nonzero. The real spectral denominator is bounded by $1/|k|$ even for
degenerate intermediate energies. Each change in occupation is bounded by
$e^{-\mathcal L E}$. Cauchy-Schwarz and completeness bound the two vertex
matrix elements by $B_n^2$. Shifting the continuous momentum removes
$\Omega$ from the occupation integral; the two angular and two energy
signs give the displayed conservative multiplicity. Rotational invariance
and the absolute Euclidean metric amplitude bound the contact as in the
static proof. The bounds on $M_n$ and $K_\nu$ from that note reduce these
infinite sums to polynomial geometric sums. The code evaluates these
closed bounds, with no angular cutoff.

This error estimate is proved on the real Euclidean frequency axis.
It is not a bound obtained by substituting a complex Lorentzian frequency
into its numerical value.

## Retarded domain below the first angular pair threshold

A separate vacuum spectral argument gives a useful Lorentzian domain. At
fixed nonzero allowed $k$, the pure LLL poles cancel in the finite-size
difference as proved above. Every remaining noncontact intermediate state
contains at least one massive angular excitation. On the infinite line,

$$|p|+\sqrt{(p+k)^2+M^2}\ \ge\ \sqrt{k^2+M^2}.$$

For $p>0$ the left side increases with $p$; for $p<0$ it decreases towards
$p=0$, proving the minimum. On the antiperiodic circle its minimum is
$\pi/\mathcal L+\sqrt{(|k|-\pi/\mathcal L)^2+M^2}$, which is strictly
larger. Two massive excitations cost at least
$\sqrt{k^2+(M_i+M_j)^2}$, by the usual energy-momentum inequality.

Consequently the retarded compact-minus-line response has no absorptive
spectral weight for

$$|\omega|<\omega_*\equiv\sqrt{k^2+(q+1)/R^2}.$$

After the common UV subtraction it is analytic in $\omega$ inside that
disk. To see why the infinite spectrum does not create a lower singularity,
use the vacuum Lehmann representation: all noncanceling spectral support
lies outside this interval; finitely many local UV subtractions make the
free-field spectral integral locally uniformly convergent, and the
subtraction/contact polynomials are entire. This is a gap/analyticity
statement, not a claimed new general spectral theorem or a quantitative
Lorentzian remainder estimate. Particular source combinations may have a
higher threshold. It supplies a domain for a subsequent controlled
low-frequency matching calculation, not a stability or front-speed result.

## Verification and remaining physical requirements

The [response script](../checks/casimir_cylinder_response.py) and
[receipt](../evidence/casimir-cylinder-response.json) retain the independent
spinor trace, general partial fractions and exact curl-residue identities.
Four direct Matsubara sums agree within $8.6\times10^{-15}$; four raw versus
reduced curl integrals agree within $7.9\times10^{-18}$. Three static limits
reproduce the previous integral within $3.5\times10^{-18}$. Five real
Euclidean frequency samples give four-source response matrices and
all-level error bounds. Their signs are numerical observations for those
samples, not a global frequency-dependent sign theorem. Floating quadrature
errors and bound evaluations are not interval certificates.

The [full vertex check](../checks/lll_axial_vertex.py) independently derives
the warped-background curl term from the tetrad and spin connection. The
new source and the static script it imports are hash-pinned in the receipt.
No Lean compiler was invoked.

The subsequent [Ward completion](casimir-cylinder-ward-completion.md)
restores the angular gauge and axial tensor partner sources and includes
background tadpole terms for this selected six-source block. Its finite
pullback and independent Dirac commutator checks also expose opposite
static matter curvatures on metric paths with the same linear tangent.
The next physical requirement is a controlled curved total response.
Other perturbation sectors, internal gauge interactions, renormalized infinite-line
response, gravity and global exterior matching also remain outside this
free-cylinder result. The static negative sign still cannot be interpreted
as a physical wormhole instability.

Remaining to-do list: derive the curved total response;
complete the closest-source contribution comparison and primary/fallback
selection; establish the physical consequence; then GCP-certify suitable
complete endpoints and prepare the manuscript.
