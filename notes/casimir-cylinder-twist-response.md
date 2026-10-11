# Static Casimir response to an axial twist at finite magnetic flux

2026-10-10, Space S051/S052. A fixed-background free-Dirac calculation now
resolves a component of the previously assumed Casimir response. On a long
magnetic product cylinder, an affine axial metric twist has a strictly
negative quadratic **finite-size vacuum-energy** coefficient. The result
includes virtual LLL/higher-level transitions and an analytic bound on the
entire remaining massive sector. It is not an instability theorem for the
metric: classical gravity, electromagnetism and the infinite-line vacuum
response have not been included in this energy comparison.

The leading large-flux terms reproduce the spatial drift prescription. A
finite-flux correction survives and is not fixed by the number of zero modes
alone. Its novelty and its consequence in a curved, dynamical wormhole remain
under review. This advances the [specified cylinder test](finite-flux-response-audit.md)
without claiming that its full dynamical kernel or Ward-identity audit is
complete.

## Defined observable and result

Use one free massless four-dimensional Dirac field of unit charge on

$$ds^2=-dt^2+dz^2+R^2d\Omega^2,\qquad
z\sim z+\mathcal L,\qquad A_\varphi=\frac q2\cos\theta,$$

with integer $q\ge3$, antiperiodic fermions in $z$, zero temperature, zero
chemical potential and no added longitudinal holonomy. Photons and the
metric are prescribed external sources in this one-loop matter calculation.
Put $u=\cos\theta$, $Y=P_2(u)$, $k=2\pi/\mathcal L$ and
$f(z)=\sqrt2\cos(kz)$, so $\mathcal L^{-1}\int f^2dz=1$. The source path is

$$g_{z\varphi}(\epsilon)=-\epsilon h f(z)(1-u^2)Y'(u),\qquad
A_z(\epsilon)=\epsilon a f(z)Y(u),$$

with all other metric and gauge components fixed. In particular $g_{zz}=1$
is held fixed through second order; this matters for the contact term.
Sufficiently small $|\epsilon h|<\sqrt2R/3$ keeps the spatial metric positive.

Subtract the infinite-line vacuum energy per period for the same periodic
sources, regulator and local couplings from the compact-circle energy.
Write its quadratic coefficient as

$$\Delta E(\epsilon)-\Delta E(0)
=\epsilon^2\mathcal L\,\mathcal Q(a,h)+O(\epsilon^3).$$

Local state-independent counterterms cancel in this difference. The
comparison is between ground states at the specified sources, rather than
an energy budget for preparing or sustaining the geometry.

**Static response theorem.** For $\mathcal L/R\ge16$,

$$\mathcal Q(0,h)<-
\frac{127}{128}\frac{\pi D(q)}{12R^2\mathcal L^2}h^2
\quad(h\ne0),$$

$$D(q)=\frac{18q(2q^3+3q^2+q+3)}{5(q+1)^2(q+2)^2}>2.$$

The coefficient refers to this explicit affine source path. It is not the
total gravitational canonical energy, and a Hessian on a prescribed
nonstationary gravitational background cannot be promoted to a
gauge-independent metric-stability criterion. No physical MMP stability or
shortcut conclusion follows from this sign alone.

## Angular reduction and exact virtual-transition weights

The sphere Dirac spectrum has $q$ zero modes and paired nonzero eigenvalues
$\pm\sqrt{n(n+q)}/R$, each sign with multiplicity $q+2n$, $n\ge1$.
Greiter and Thomale,
[Landau level quantization of Dirac electrons on the sphere](https://arxiv.org/abs/1807.05816v2),
derive this spectrum in the HTML display equations (13)-(25), with
$2b_0=q$. The corresponding equations (9)-(17) in their
[publisher edition](https://doi.org/10.1016/j.aop.2018.04.020), pages 35-36,
have been compared directly. The spectrum itself is established
prior work.

For an LLL mode $\eta_m$ with $j=(q-1)/2$ and $-j\le m\le j$, the normalized
weight is $\rho_m(u)\propto(1-u)^{j-m}(1+u)^{j+m}$. Acting on $F(u)\eta_m$,
the negative-chirality squared angular Dirac operator is

$$H_mF=-(1-u^2)F''+[(q+1)u-2m]F'. $$

This follows by multiplying the first-order operators in MMP Appendix C;
it is self-adjoint for $\rho_mdu$. Its polynomial eigenfunctions are
$P_n^{(j-m,j+m)}(u)$ with eigenvalues $\lambda_n=n(n+q)$.
The centered quadratic $Y-\langle Y\rangle_m$ therefore has components
only at $n=1,2$.

The full static spatial Dirac vertex, checked in both transverse chiralities,
is $-i\gamma^zV$, where

$$V_\sigma=aY-\frac h{R^2}T_\sigma,\qquad
T_\sigma=Y'(m-qu/2)-\frac\sigma4(1-u^2)Y'',\quad\sigma=\pm1.$$

There is no longitudinal-curl vertex for this static spatial source.
On the LLL chirality $\sigma=-1$, Legendre's equation gives the useful
operator identity

$$T_-=(3/2-H_m/2)Y.$$

Thus a transition to level $n$ couples to $a+b_nh/R^2$, where

$$\lambda_1=q+1,\quad \lambda_2=2(q+2),\qquad
b_1=(q-2)/2,\quad b_2=(2q+1)/2.$$

Summing the squared normalized matrix elements over all LLL modes gives

$$S_1=\frac{6q(q-1)}{5(q+1)(q+3)},\qquad
S_2=\frac{12q}{5(q+2)(q+3)}.$$

One direct proof decomposes the centered $Y$ using the projectors
$(\lambda_2-H_m)/(\lambda_2-\lambda_1)$ and
$(H_m-\lambda_1)/(\lambda_2-\lambda_1)$, then integrates with $\rho_m$.
The moments obey

$$\langle u^r\rangle_m=
\frac{2m\langle u^{r-1}\rangle_m+(r-1)\langle u^{r-2}\rangle_m}{q+r}.$$

This proves the formulas for arbitrary $q$. Separate raw Jacobi-polynomial
integrals check them without using that projector decomposition.

Define $G_{TT}=\sum_{n=1}^2S_nb_n^2/\lambda_n$ and
$W=\sum_m\langle(1-u^2)(Y')^2\rangle_m=6q/5$. Then

$$G_{TT}=\frac{3q(q^4+4q^2+9q-5)}{10(q+1)^2(q+2)^2},\qquad
W-4G_{TT}=D(q).$$

For the other entries, $G_{YY}=\sum S_n/\lambda_n$ and
$G_{YT}=-\sum S_nb_n/\lambda_n$.

## Quadratic determinant and metric contact term

For $D=D_0+\epsilon D_1+\epsilon^2D_2+\cdots$, the Euclidean matter action
$-\log\det D$ has quadratic coefficient

$$-\operatorname{Tr}(D_0^{-1}D_2)
+\tfrac12\operatorname{Tr}(D_0^{-1}D_1D_0^{-1}D_1).$$

The first-order coframe $e^\varphi=R\sin\theta(d\varphi+\epsilon v dz)$,
$v=-hfY'/R^2$, must be supplemented by
$e^z=[1-\epsilon^2R^2\sin^2\theta\,v^2/2]dz$ to describe the affine metric
path. The inverse frame then has the opposite quadratic correction. The
Dirac kinetic contact term contributes

$$\mathcal Q_{\mathrm{contact},0}
=-\frac{\pi W}{12R^2\mathcal L^2}h^2.$$

Here the LLL Casimir energy per two-dimensional Dirac field is
$-\pi/(6\mathcal L)$, obtained from
$-(2/\pi)\int_0^\infty dp\log(1+e^{-\mathcal Lp})$.
At zero chemical potential the longitudinal vector and axial one-point
currents vanish by momentum reflection; rotational invariance excludes a
tangential one-point current on the sphere. The connection terms therefore
add no finite-size tadpole. Densitizing the
fermion produces a similarity transformation and state-independent measure
terms, which cancel in the subtraction. The massive contact contribution
is retained and bounded below, rather than identified with this LLL term.

Both vertices inside a purely LLL loop are longitudinal spatial gauge
potentials. At nonzero $k$ they are removable by a periodic phase for each
axisymmetric LLL species; the loop contributes zero in this static channel.
This does not assert that the general dynamical Schwinger kernel vanishes.

An independent coordinate check fixes the contact sign. For the auxiliary
$l=1$ choice $Y=u$, take $a=-qh/(2R^2)$. These linear sources are the pullback
of the product background under $\varphi\mapsto\varphi+\epsilon\xi(z)$,
with $\xi'=-hf/R^2$. The full vertex becomes $V=-hm/R^2$ in each angular
$m$ sector, so its nonzero-$k$ bubble vanishes. The exact pullback also has
$\delta^{(2)}g_{zz}=h^2f^2\sin^2\theta/R^2$. Its total energy change is zero
by coordinate invariance. Removing that quadratic metric change to recover
the affine path leaves minus its first variation against the unperturbed
stress. For the LLL, varying $-\pi/(6\mathcal L\sqrt{g_{zz}})$ gives precisely
the negative contact coefficient above, with the $l=1$ angular weight.
This verifies the path dependence and sign; it is not a replacement for the
full $l=2$ dynamical Ward-identity audit.

## Exact mixed loop and its sign bound

Rotate to Euclidean time and sum the antiperiodic momenta in the compact
$z$ direction. Set $M_n=\sqrt{\lambda_n}/R$. Including both orientations
of the mixed loop, its contribution is

$$\mathcal Q_{0H}(a,h)=\sum_{n=1}^2S_nJ(M_n,\mathcal L,k)
\left(a+\frac{b_nh}{R^2}\right)^2,$$

$$J(M,\mathcal L,k)=\frac1\pi\int_0^\infty dp\left[
\frac{4M^2p\,n_F(\mathcal Lp)}{(M^2+k^2)^2+4k^2p^2}
-\frac{2M^2(M^2+k^2+2p^2)\,n_F(\mathcal L\sqrt{p^2+M^2})}
{[(M^2+k^2)^2+4k^2p^2]\sqrt{p^2+M^2}}\right],$$

where $n_F(x)=(e^x+1)^{-1}$. This includes the massive thermal-pole term.
It is exact for the mixed sector, not an expansion in $k/M$.

To check the sum, take the time-component trace numerator
$2[p^2-r(r+k)]$ over
$(r^2+p^2)[(r+k)^2+p^2+M^2]$. Its partial fractions have even numerators

$$B=\frac{4M^2p^2}{\mathcal D},\qquad
F=-\frac{2M^2(M^2+k^2+2p^2)}{\mathcal D},\qquad
\mathcal D=(M^2+k^2)^2+4k^2p^2.$$

The odd terms sum to zero, since $k$ is a bosonic circle momentum. Using
$\mathcal L^{-1}\sum_r(r^2+E^2)^{-1}-(2E)^{-1}=-n_F(\mathcal LE)/E$
and the opposite sign of the spatial-component trace yields $J$ above.
Direct large Matsubara sums independently reproduce the partial-fraction
answer.

Every term in the second integrand is nonnegative and
$\mathcal D\ge M^4$. Consequently

$$J(M,\mathcal L,k)\le\frac{\pi}{3M^2\mathcal L^2}.$$

At fixed $q$ and $k=2\pi/\mathcal L$, this upper bound is the leading
long-cylinder coefficient. An explicit absolute error bound follows by
expanding $1/(1+x)$ with a bounded remainder:

$$0\le\frac{\pi}{3M^2\mathcal L^2}-J
\le\frac{\pi}{3M^2\mathcal L^2}
\left[\frac{2k^2}{M^2}+\frac{k^4}{M^4}
+\frac{14\pi^2k^2}{5M^4\mathcal L^2}\right]
+\frac2\pi\left[K_2(M\mathcal L)+\frac{k^2}{M^2}K_0(M\mathcal L)\right].$$

The Fermi integrals used here are $\pi^2/(12\mathcal L^2)$ and
$7\pi^4/(120\mathcal L^4)$. The $K_\nu$ terms bound the massive-pole
integral using $n_F(x)\le e^{-x}$.

## Bound on the entire massive sector

The remaining quadratic response consists of a loop with both angular
propagators outside the LLL and the massive contact tadpole. It is not set
to zero. In an angular level $n$, $|m|\le(q-1)/2+n$ and the explicit vertex
above gives

$$\|V\psi_n\|\le B_n\|\psi_n\|,\qquad
B_n=|a|+\frac{3(q+n)}{R^2}|h|.$$

After interpreting compact $z$ as Euclidean time, the unperturbed Hermitian
one-particle Hamiltonian has energies $\pm\sqrt{p^2+M_n^2}$. The real part
of a two-point spectral denominator is bounded by $1/|k|$. The thermal
occupation change has magnitude at most $e^{-\mathcal L\sqrt{p^2+M_n^2}}$.
Applying completeness to the second vertex and counting the two angular
signs and two energy signs gives the conservative bound

$$|\mathcal Q_{HH,\mathrm{bubble}}|
\le\frac4{\pi|k|}\sum_{n=1}^\infty(q+2n)B_n^2M_nK_1(\mathcal LM_n).$$

This argument does not assume a nonzero difference between two massive
energies; the nonzero external $k$ controls that denominator.
For the contact term, rotational invariance gives angular trace
$2(q+2n)\,6/5$ in level $n$. The compact-minus-line kinetic tadpole and
$\int_0^\infty E e^{-\mathcal LE}dp
=M^2[K_0(M\mathcal L)+K_2(M\mathcal L)]/2$ yield

$$|\mathcal Q_{H,\mathrm{contact}}|
\le\frac{6h^2}{5\pi R^2}\sum_{n=1}^\infty(q+2n)M_n^2
[K_0(\mathcal LM_n)+K_2(\mathcal LM_n)].$$

Both sums converge. They bound all higher levels, including their
degeneracy and growing metric vertex, and can be evaluated without a cutoff.
For explicit upper bounds use

$$M_n\ge M_1+(n-1)/R,\qquad M_n\le(n+q/2)/R,$$

$$K_\nu(x)\le e^{-x}\sqrt{\frac\pi{2x}}e^{\nu^2/(2x)}.$$

The last inequality follows directly from the positive integral for
$K_\nu$ and $\cosh t\ge1+t^2/2$. The remaining sums are polynomial
geometric series through degree four. The supplied code evaluates these
closed bounds, not a finite-level approximation to them.

## Proof of the uniform negative coefficient

For $a=0$, the exact mixed-loop upper bound and the LLL contact give

$$\mathcal Q_{0H}(0,h)+\mathcal Q_{\mathrm{contact},0}
\le-\frac{\pi D(q)}{12R^2\mathcal L^2}h^2.$$

Write $d=\mathcal L/R\ge16$, $x=\sqrt{q+1}\ge2$. The numerator of
$D(q)-2$, after $q=y+3$, is
$26y^4+306y^3+1238y^2+1908y+698$, over a positive denominator.
Thus the magnitude on the right is at least
$\pi h^2/(6R^4d^2)$.

For $\nu=0,1,2$ and $\mathcal LM_1\ge32$, the preceding Bessel bound is
less than $e^{-\mathcal LM_n}$. Also
$\sum_{n\ge1}n^4r^{n-1}\le2$ for $r=e^{-d}<1/64$.
Using $q+2n\le2(q+n)$, $q+n\le n(q+1)$ and $M_n\le(q+n)/R$, the ratios
of the massive bubble and contact bounds to that lower magnitude are at most

$$\frac{432}{\pi^3}d^3x^8e^{-dx},\qquad
\frac{288}{5\pi^2}d^2x^6e^{-dx},$$

respectively. Each decreases in both $d\ge16$ and $x\ge2$.
Using $\pi>3$ and $e>2$, their maxima are smaller than $1/256$ and
$1/32768$, respectively. Their sum is less than $1/128$.
This proves the stated strict negative bound for the **full** static
finite-size coefficient, including the infinite massive sector.

## Drift comparison and scientific limits

Prior geometric response calculations already include mixed electromagnetic
and metric kernels and quadratic metric contact terms. For example,
[Abanov and Gromov](https://arxiv.org/abs/1401.3703v1), sections II.1-II.2,
study a nonrelativistic two-dimensional gas at integer Landau filling.
That state and action differ from this four-dimensional Dirac vacuum.
Likewise, [Fischetti, Wallis and Wiseman](https://arxiv.org/abs/2003.09428v2),
sections 2.1, 3.2-3.3 and 3.5, find lower free energy under shape deformations of
an unmagnetized spatial two-sphere, using analytic and numerical evidence.
Their fixed-area geometry subtraction differs from our compact-minus-line
subtraction. These precedents rule out a broad claim that a negative
geometric Casimir response is new; the precise finite-flux coefficient,
uniform all-level bound and any dynamical consequence need their own
contribution comparison.

The leading derivative coefficient from the mixed sector is the Gram matrix
of $(Y,T_-)$ with $H_m^{-1}$ on its nonzero subspace. Up to a conformal term
that has zero contraction with the flat-cylinder Casimir stress, its
effective two-dimensional metric correction is

$$\sum_m\delta\gamma_{ab}
=4R^2G_{YY}a_aa_b-4G_{YT}(a_ah_b+h_aa_b)
+\frac{4G_{TT}-W}{R^2}h_ah_b.$$

The spatial component is derived above. This algebraic expression is not
a completed time-dependent determinant; longitudinal-curl vertices must be
included before using the temporal or mixed dynamical components.
At large $q$, $G_{YY}\sim W/q^2$ and $G_{YT}\sim-W/(2q)$, which reproduce
the leading spatial field-line drift response. But $4G_{TT}-W=-D(q)$,
with $D(q)\to36/5$. The finite-flux metric term is order one in the species
sum, while the background Casimir energy is order $q$.

This is a controlled counterexample to treating the angularly local spatial
drift response as an exact finite-flux Casimir response. It does not refute
W22's stated leading-gradient truncation. State-independent counterterms
cannot remove this compact-minus-line correction. Its sign alone neither
invalidates the earlier conditional positive-energy theorem nor establishes
that the complete physical system has a negative mode. The dynamical
response, total semiclassical action, curved background and global exterior
completion are separate requirements.

The [script](../checks/casimir_cylinder_twist.py) and
[receipt](../evidence/casimir-cylinder-twist.json) retain the exact identities,
52 independent Jacobi-mode checks, three independent Matsubara sums, five
mixed-loop error-bound checks, the Casimir baseline and eight numerical
evaluations of the all-level bound.
The [full static Dirac vertex check](../checks/lll_axial_vertex.py) verifies
both angular chiralities. Numerical quadrature is supplementary to the
ordinary sign proof; no Lean certificate exists yet.

Remaining to-do list: complete the closest-source novelty comparison;
derive the temporal/curl response and its Ward identities; assess the
curved-background and stability consequence; then select and GCP-certify
the complete surviving endpoint and prepare the manuscript.
