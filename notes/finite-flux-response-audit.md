# Finite magnetic flux: projected charge response and axial metric coupling

2026-10-10, Space S051/S052. The angularly local charge response is a
large-flux approximation. Direct projection gives an exact finite-flux form
factor and a reciprocal axial metric/charge vertex in the **free lowest
Landau level (LLL) sector**. These corrections can precede the retained
gravitational backreaction order. Including their reciprocal coupling in
the prescribed drift model preserves the positive-energy structure.

This advances the microscopic source audit; it does not derive the complete
Casimir-state response, integrate out higher Landau levels, or prove physical
MMP stability. The form factor is the established spherical Berezin
eigenvalue. Novelty of the response application and energy extension remains
under comparison. The research decision is to examine a joint finite-flux
and backreaction expansion before assigning physical precision to an
alpha-only correction.

## Scope of the projection

Start with one massless unit-charge Dirac field, magnetic flux $q>0$ an
integer, and

$$ds^2=N(t,z)^2(-dt^2+dz^2)+R(t,z)^2d\Omega^2,
\qquad A_\varphi=\tfrac q2\cos\theta.$$

Retain the $q$ sphere-Dirac zero modes, with
$j=(q-1)/2$ and $m=-j,\ldots,j$. The chiral angular spinors are those in
Maldacena, Milekhin and Popov,
[1807.04726v3](https://arxiv.org/html/1807.04726v3), Appendix C, especially
(C.18). Their normalized scalar factors in a spin frame are

$$\eta_m=\sqrt{\frac q{4\pi}{q-1\choose j-m}}
\left(\sin\frac\theta2\right)^{j-m}
\left(\cos\frac\theta2\right)^{j+m}e^{im\varphi}.$$

The full spinor includes the negative transverse-chirality basis vector.
Their norm uses $d\Omega$; phases on overlapping gauge/spin patches cancel
in the quantities below. The background separation is exact. Treating this
projected sector as the full interacting theory is an additional low-energy
approximation. We examine infinitesimal axisymmetric sources, so angular
vertices are diagonal in $m$ and ordinary Abelian bosonization suffices for
the current response. General angular sources require the matrix current
algebra; counting its generators as independent quantum scalar fields would
overcount the central charge.

## Known projection kernel, independently reproduced

For the zero-mode projector $P_q$, the squared kernel is

$$|P_q(\Omega,\Omega')|^2=
\left(\frac q{4\pi}\right)^2
\left(\frac{1+\hat n\cdot\hat n'}2\right)^{q-1}.$$

This follows by summing the normalized binomial basis. Define
$B_qf=(4\pi/q)\int|P_q(\Omega,\Omega')|^2f(\Omega')d\Omega'$.
Rotational invariance diagonalizes $B_q$ on harmonics. Its eigenvalue is

$$\beta_l(q)=\frac q2\int_{-1}^1
\left(\frac{1+u}2\right)^{q-1}P_l(u)\,du
=\begin{cases}
\displaystyle\frac{(q-1)!q!}{(q-1-l)!(q+l)!},&0\le l<q,\\
0,&l\ge q.
\end{cases}$$

Proof: insert Rodrigues' formula for $P_l$ and integrate by parts $l$ times.
For $l<q$ the remaining beta integral gives the factorial ratio. For $l\ge q$
the $l$th derivative of the degree-$(q-1)$ weight vanishes. Endpoint terms
vanish because the Rodrigues polynomial has zeros of order $l$ at both ends.

For a unit-normalized harmonic $Y_l$,

$$\operatorname{Tr}[(P_qY_lP_q)^\dagger(P_qY_lP_q)]
=\frac q{4\pi}\beta_l(q).$$

This formula is already known: Demni, Mouayn and Yaqine,
[2021 publisher edition](https://doi.org/10.15407/mag17.04.422), equation
(6.10), agrees with $2\nu=q-1$ and spectral label $l$; their Fubini-Study
eigenvalue is $4l(l+1)$. The end of section 5 of
[1610.02875v1](https://arxiv.org/html/1610.02875v1) contains the corresponding
gamma ratio. We use the ratio directly and do not import a Laplacian
normalization from an infinite-product expression. The general higher-level
formula and its full dependency chain are not needed or claimed reviewed.

For $l<q$, the product form gives

$$\beta_l(q)=\prod_{n=1}^l\frac{q-n}{q+n},\qquad
\beta_l(q)=1-\frac{l(l+1)}q+\frac{l^2(l+1)^2}{2q^2}+O(q^{-3})$$

at fixed $l$. Moreover $0\le\beta_l\le1$ and
$1-\beta_l\le\min(1,l(l+1)/q)$ for all $l$. The bound follows from
$1-\prod(1-a_n)\le\sum a_n$ with $a_n=2n/(q+n)$; above the cutoff the
right side is at least one. Thus the angularly local approximation is
controlled by $l(l+1)/q$, but that error need not be below the requested
backreaction precision.

The planar comparison is also explicit. For $l<q$, expanding the logarithm
of each product factor gives

$$0\le-\log\beta_l-\frac{L}{q}
\le\frac{L^2}{6q^3[1-(l/q)^2]},\qquad L=l(l+1).$$

Indeed the remainder of $2\operatorname{artanh}(x)$ after $2x$ is at most
$2x^3/[3(1-x^2)]$, and $\sum_{n=1}^l n^3=L^2/4$.
As $q\to\infty$, $L/q\to x<\infty$ gives $\beta_l\to e^{-x}$.
With $p_\perp^2=L/R^2$ and $|eB|=q/(2R^2)$, this is the known planar factor
$e^{-p_\perp^2/(2|eB|)}$. Hattori and Itakura,
[2205.04312v2](https://arxiv.org/html/2205.04312v2), equations (14)-(21) and
Appendix A.1, derive that factor by the transverse Gaussian loop integral.
The spherical formula is not the first finite-width correction to a
Schwinger response.

The free projected fermion loop has a common two-dimensional current kernel
multiplied by the matrix trace of its angular vertices. Consequently its
charge-response coefficient is the local-lift value times $\beta_l(q)$.
In the baseline notation,

$$s_l=s_0\beta_l(q),\qquad s_0=\frac{g_c^2q}{4\pi^2}.$$

This statement concerns the two-point charge kernel, with a consistent
choice of the two-dimensional vacuum/retarded response and treatment of
global zero modes. It does not compute all local contact terms in the
four-dimensional determinant. The vanishing at $l\ge q$ is a property of the
LLL projection, not vanishing of the complete QED response.

## Axial metric vertex from the Dirac operator

Use the first-order coframe

$$e^0=Ndt,\quad e^1=Ndz,\quad e^2=Rd\theta,\quad
e^3=R\sin\theta(d\varphi+\epsilon v_a dx^a),$$

and $D_\mu=\partial_\mu+\Omega_\mu-iA_\mu$. Thus
$h_{a\varphi}=R^2\sin^2\theta\,v_a$.
The [tetrad check](../checks/lll_axial_vertex.py) constructs the connection,
spin connection and Dirac operator, retaining $t,z$ derivatives of $N,R,v_a$.
Projection to negative transverse chirality gives a shift of the effective
two-dimensional gauge potential by

$$v_a\left(m-\frac q2\cos\theta\right)
-\frac{\sin\theta}{4}\partial_\theta v_a.$$

The last term is essential; dropping the spin connection gives a different
answer. For $h_{a\varphi}=-h_a(1-u^2)P_l'(u)$, $u=\cos\theta$, one has
$v_a=-h_aP_l'(u)/R^2$. Define
$Q_{lm}=\int |\eta_m|^2P_l(u)d\Omega$. The projected potential is

$$\mathcal A_{a,m}=Q_{lm}\left(a_a-\frac{l(l+1)}{4R^2}h_a\right).$$

To prove the angular integral, write the normalized density as
$\rho_m(u)\propto(1+u)^{j+m}(1-u)^{j-m}$. It satisfies
$(1-u^2)\rho_m'=2(m-ju)\rho_m$. Integration by parts and Legendre's
equation give

$$\left\langle P_l'(m-ju)\right\rangle_m
=\tfrac12 l(l+1)Q_{lm},$$

$$\left\langle P_l'(m-qu/2)+\tfrac14(1-u^2)P_l''\right\rangle_m
=\tfrac14 l(l+1)Q_{lm}.$$

For $l=1$, $Q_{1m}=2m/(q+1)$. A rigid rotation together with its dragged
monopole potential reduces to the expected pure-gradient phase with angular
momentum $m$; this independently fixes the relative sign of the metric term.
The stability theorem below still excludes $l=0,1$ from its domain.

This is a current vertex of the mixed stress $T_{a\varphi}$, not a claim that
the two-dimensional scalar stress $t_{ab}$ has a nonzero vacuum two-point
function with the vector current. The latter can vanish while the former
survives. Sadhukhan's
[2610.09847v1](https://arxiv.org/html/2610.09847v1), sections 5-6, uses the
leading angularly local lift and leaves gradient corrections outside its
scope. The calculation identifies one such correction; it does not refute
that stated leading truncation or establish a full four-dimensional response.

## Why the expansion order matters

Set $L=l(l+1)$, $C=R^2$, and $w_a=4a_a/q$. The projected source depends on

$$w_a^{\rm eff}=w_a-\delta h_a/C,\qquad \delta=L/q.$$

For $G=h_0'-\dot h_1-(C'/C)h_0$ and $E=\dot w_1-w_0'$, its electric
field is $E_{\rm eff}=E+\delta G/C$. Both the form factor and the reciprocal
vertex have leading corrections of order $1/q$ at fixed multipole.

The source parameters obey

$$\alpha=\frac{g_c^2}{4\pi^2q}=\frac{s_0}{q^2}.$$

Along the sequence $q\to\infty$ with fixed $s_0>0$,
$q^{-1}=\sqrt{\alpha/s_0}$. Therefore

$$s_l-s_0=-L\sqrt{s_0\alpha}+O(\alpha),\qquad
\delta=L\sqrt{\alpha/s_0}.$$

For nonzero multipoles these terms precede an $O(\alpha)$ backreaction
correction. At fixed weak $g_c$, the relative form-factor error divided by
$\alpha$ tends to $4\pi^2L/g_c^2$, rather than zero. Thus an alpha-only
truncation with unexpanded screening is not, by itself, a uniform physical
accuracy guarantee for the charge kernel. This is a precise hierarchy
obstruction, not an instability theorem. No cancellation or size of a full
observable's shift is inferred without calculating its other contributions.
Higher levels, renormalization and determinant contact terms still require
matching at the intended order.

## Reciprocal correction and positive energy

Retain the **assumed ground-state drift response** of the
[earlier theorem](global-drift-energy.md), and add the projected charge
channel with $s=s_l>0$ and constant $\delta=L/q$. Variation of the same charge
coupling must supply both the Maxwell current and the mixed Einstein source.
In radial coordinates with $ds_2^2=-Sdt^2+dr^2/S$, its auxiliary-field action
is proportional to

$$I_X=\frac{s\mathcal N_l}{\widetilde g^2}\int dt\,dr\,
\left[\frac{\dot X^2}{2S}-\frac S2(X')^2-XE_{\rm eff}\right],$$

where $\mathcal N_l=\int Y_l^2d\Omega$ and the normalized electromagnetic
coupling satisfies $\kappa/\widetilde g^2=1/2$. Its variations with respect to
$(w_0,w_1,h_0,h_1)$, after removing $\mathcal N_l/\widetilde g^2$, are
$(-sX',s\dot X,\delta sX'/C,-\delta s\dot X/C)$. The angular norm of the
axial vector harmonic is $L\mathcal N_l$; this supplies the factor $1/L$
when converting the metric variations to mixed tensor components.

In the harmonic-stripped equations, the extra terms on the Einstein left
sides are $\delta sS X'/(2LC)$ and $\delta s\dot X/(2LCS)$. Their mixed
stress is separately conserved. Keeping only a modified anomaly equation
would miss this reciprocal contribution.

Set

$$\Phi=CE-sX,\qquad F_g=CG+2\Phi/L-\delta sX/L,\qquad
H=F_g-2\Phi/L+\delta sX/L.$$

Reconstruction retains $h_0=SF_g'/k$, $h_1=\dot F_g/(kS)$ and the same
$w_a$ formulas in the earlier theorem, where $k=L-2$. The equations become

$$\frac{\ddot F_g}{SC}-\left(\frac S C F_g'\right)'+\frac{k}{C^2}H=0,$$

$$\frac{\ddot\Phi}{SA_r}-\left(\frac S{A_t}\Phi'\right)'
+\frac L C(\Phi+sX)-\frac{2L}{C^2}H=0,$$

$$\frac{\ddot X}S-(SX')'+\frac{\Phi+sX}C+\frac{\delta}{C^2}H=0.$$

Their positive kinetic and gradient matrices $K,B$ are unchanged, and the
potential becomes

$$z^TVz=\frac{(LF_g-2\Phi+\delta sX)^2}{C^2}
+\frac L C(\Phi+sX)^2.$$

The same energy/flux proof therefore excludes finite-energy growing modes
for $S,C,A_t,A_r>0$, $L>2$, $s>0$ and admissible no-incoming-energy boundaries.
It holds for any real constant $\delta$, so this particular correction does
not require a smallness estimate to preserve positivity. The constant kernel
is now $(-s(2+\delta)/L,-s,1)$; radial characteristic cones are unchanged.
For the physical projection this three-field assertion covers $2\le l<q$.
It does not apply a spurious positive scalar kinetic coefficient above the
LLL cutoff, where $s_l=0$ and the matter channel must be removed/reanalyzed.

This is a conditional theorem for the enlarged defined response class.
Projection of the linear Dirac vertex does not determine second-order
contact terms or the full vacuum metric response, so it cannot by itself
promote this enlarged model to complete semiclassical gravity.

## Source comparison and reproducibility

Gralla's [1811.08422v2](https://arxiv.org/html/1811.08422v2), sections 2-4,
constructs a covariant bosonized coherent-plasma action. Its coherent-field
restriction is relevant: quantizing one common scalar does not restore the
central charge of all $q$ independent fermions. Iqbal's
[1909.12609v1](https://arxiv.org/html/1909.12609v1), Appendix A, explicitly
discusses the omitted gapless modes and limits of that matching. MMP's
Appendix A already distinguishes the screened singlet from the remaining
massless sector. Neither plasma action alone establishes the required
Casimir-state determinant. These are prior structures to use and compare,
not a first-bosonization claim.

Karabali and Nair's
[1606.06405v1](https://arxiv.org/html/1606.06405v1), section 2, already combines
the electromagnetic potential and spin connection in a filled integer Hall
state. Its odd-dimensional topological action is a different problem from
the gapless four-dimensional Dirac vacuum here; neither its coefficient nor
its state can simply be imported. A broad claim to discovering metric/gauge
mixing is unavailable. The specific spherical axial coefficient and its
coupled dynamical consequence remain the comparison target.

Gorbar, Hashimoto and Miransky's
[hep-th/0501135v1](https://arxiv.org/html/hep-th/0501135v1), equations (6),
(18) and (24)-(26), give a further limitation on a general decoupling argument.
In their massive, large-flavor, flat constant-field calculation, the LLL
two-point vertex is transverse, while its three-point vertex needs a summed
higher-level contribution to restore the ordinary gauge Ward identity.
Their result does not invalidate the two-point projection computed here or
prove the same effect for a massless spherical Casimir state. It does show
why the size of an individual Landau gap alone is not a bound on the complete
omitted determinant. The next matching calculation must check the summed
response and the Ward identities at the claimed order.

The [projection check](../checks/lll_projection.py) compares the kernel with
explicit mode integrals for $q=1,\ldots,12$, including vanishing above the
cutoff, orthogonality and metric vertices. Its
[receipt](../evidence/lll-projection.json) records 102 kernel and 102 mode-norm
cases, 286 orthogonality cases and 806 metric-vertex integrals. The
[Dirac receipt](../evidence/lll-axial-vertex.json) verifies the independently
constructed tetrad/spin-connection calculation. The
[energy extension](../checks/finite_flux_energy.py) checks reconstruction,
reciprocity, the action variations, local energy law and kernel in 18 exact
identities, retained
in its [receipt](../evidence/finite-flux-energy.json). The proofs above give
the arbitrary-parameter arguments; finite cases are independent checks.
All run without Lean. No new GCP certification or manuscript is claimed.

## Bounded matching calculation

The subsequent [static response theorem](casimir-cylinder-twist-response.md)
completes the spatial sign component of this test, including contact and
an analytic all-level remainder. The later
[Euclidean response](casimir-cylinder-dynamical-response.md) includes temporal
and curl terms with a bound on the remaining massive sector. The full Ward
audit below remains required; the test as a whole is incomplete.

The next feasibility test isolates the state-dependent response on a
prescribed product cylinder, before attempting a globally completed
wormhole. Use one free massless unit-charge four-dimensional Dirac field on
$\mathbb R_t\times S^1_z\times S^2$, with circle circumference $\mathcal L$,
sphere radius $R_0$, flux $q>0$, zero temperature and antiperiodic fermions
around the circle. This is a source calculation on a fixed geometry, not a
self-consistent wormhole construction. Start with the axial/polar $l=2$
sector and $q\ge3$, with longitudinal frequencies and momenta below the
first nonzero angular Dirac eigenvalue. Internal photon and graviton loops
are outside this first free-fermion test.

Compute the difference between the compact-circle vacuum response and the
infinite-line vacuum response at the same local geometry and external
Fourier momentum. A nonzero allowed external spatial momentum avoids the
global charge/holonomy zero mode. Use the same regulator and renormalized
couplings in both terms. Local state-independent counterterms then cancel
in this difference, while a finite-size stress response remains. First
reproduce the LLL Casimir energy $-\pi q/(6\mathcal L)$ and the projected
charge kernel; derive the higher angular spectrum and matrix elements as
part of the calculation rather than replacing them by a gap estimate.

For $D=D_0+\epsilon D_1+\epsilon^2D_2+\cdots$, the quadratic coefficient
of $\log\det D$ is

$$\operatorname{Tr}(D_0^{-1}D_2)
-\tfrac12\operatorname{Tr}(D_0^{-1}D_1D_0^{-1}D_1).$$

The Euclidean fermion effective action $-\log\det D$ has the opposite sign.

Separate the LLL, mixed LLL/higher-level, and purely higher-level parts
before taking the finite-size difference. Equivalently the low-sector
Schur complement contains
$P D_2P-PD_1Q(QD_0Q)^{-1}QD_1P$, with $Q=1-P$; the higher-sector
determinant must also be retained. Fix the metric path
$g(\epsilon)=g_0+\epsilon h$: the coframe used for the linear vertex adds
quadratic $g_{ab}$ terms if used unchanged beyond first order, so it cannot
be reused for $D_2$ without that correction.

Acceptance for this test is the complete state-dependent quadratic kernel,
its charge/diffeomorphism Ward identities and a justified low-energy
remainder estimate. Compare its leading mixed stress and polarization
current with the prescribed drift relation, including any finite-flux
correction. A finite Landau-level cutoff or a Ward identity alone does not
close the test. If it fails to reproduce or control the assumed response,
narrow the stability claim or change the target before adding formalization.
Success would close this fixed-background comparison only; curved/global
matching and the final novelty decision would still be required.

Remaining to-do list: compare the microscopic vertex with its closest
predecessors; match the Casimir ground-state and higher-level/contact response
in a joint expansion; decide the primary/fallback and physical domain; then
GCP-certify the complete surviving endpoint and prepare the manuscript.
