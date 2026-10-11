# Angular Ward completion and the meaning of the Casimir sign

2026-10-10. S051/S052. Ordinary calculation with independent symbolic and
numerical checks; no Lean certification. This completes charge and angular
diffeomorphism covariance in the selected axial-metric/polar-gauge block.
It does not compute the other perturbation blocks or the curved quantum
response, and does not establish wormhole stability.

There is also a concrete consequence: the negative static coefficient in
the [affine-source calculation](casimir-cylinder-twist-response.md) becomes
strictly positive on another specified metric family with the **same linear
perturbation**. Both signs include bounds on the entire massive sector.
The two families differ at second order and the matter background has a
nonzero stress tadpole. Thus a sign assigned to this isolated, off-shell
matter Hessian cannot by itself diagnose a gravitational instability.

## Sources and established method

Use the free unit-charge four-dimensional Dirac field and compact-minus-line
subtraction of the [Euclidean response](casimir-cylinder-dynamical-response.md).
The circle is antiperiodic, the state is the vacuum, and temperature,
chemical potential and longitudinal holonomy vanish. Initially the metric is
the Euclidean product of the longitudinal plane/circle and a round sphere.
Write $x^a=(\tau,z)$, $u=\cos\theta$, $s=1-u^2$, $C=R^2$ and $Y=P_l(u)$.
Take $l\ge2$; the computed quantum kernel uses $l=2$, integer $q\ge3$ and
nonzero allowed longitudinal momentum. The monopole potential in a local
patch is $A_\phi=qu/2$.

The affine six-source path is

$$\delta g_{a\phi}=-sY'h_a(x),\qquad
\delta g_{u\phi}=-sY''H(x),\qquad
\delta A_a=a_a(x)Y,\qquad \delta A_u=b(x)Y'.$$

Other first-order components and all second-order components of this path
vanish. $H$ has dimensions of length squared, $h_a$ of length, $a_a$ of
inverse length and $b$ is dimensionless. Smooth small sources are periodic
in $z$. Time averages can be taken over a period or in a suitable large-time
limit. Monopole patches are related by their usual gauge transition; the
calculation keeps flux fixed.

The invariant axial metric combination is established perturbation theory.
[Martel and Poisson](https://arxiv.org/html/gr-qc/0502028v1#S5.SS1),
equations (49)-(55), give $\tilde h_a=h_a-\partial_a h_2/2+r_a h_2/r$.
With $h_2=2H$ and $C=r^2$ this is the combination used below. Their
Schwarzschild perturbation dynamics are not imported. The selected section
has been read; its remaining derivation dependencies and publisher edition
remain unaudited. Neither this variable nor the general Ward method is a
novelty claim.

## Linear symmetries and exact removal of angular partners

For the active displacement $\xi^\phi=-\chi(x)Y'$, combine the Lie change
of the potential with the compensating gauge transformation so that
$\delta A=i_\xi F$. The linear changes are

$$h_a\mapsto h_a+C\partial_a\chi,\quad H\mapsto H+C\chi,
\quad b\mapsto b+q\chi/2,\quad a_a\mapsto a_a.$$

Ordinary gauge parameter $\lambda(x)Y$ adds $\partial_a\lambda$ to $a_a$
and $\lambda$ to $b$. Consequently, with $\eta=H/C$,

$$\widehat h_a=h_a-C\partial_a\eta,\qquad
\widehat a_a=a_a-\partial_a b+\frac q2\partial_a\eta$$

are invariant. These identities hold also for varying $C(x)$; replacing
$C\partial_a(H/C)$ by $\partial_a H$ on a warped background would be wrong.

To determine the second-order terms, perform the finite pullback
$\phi_{\rm old}=\phi_{\rm new}+\epsilon\eta Y'$. If
$\theta_i=\partial_i(\eta Y')$ and
$d_i=(-sY'h_a,-sY''H)$ for $i=(a,u)$, its quadratic metric term is

$$K_{ij}=\theta_i d_j+\theta_j d_i+Cs\theta_i\theta_j.$$

Explicitly,

$$K_{ab}=-s(Y')^2[h_a\partial_b\eta+h_b\partial_a\eta
                    -C\partial_a\eta\partial_b\eta],$$
$$K_{au}=-sY'Y''h_a\eta,\qquad
K_{uu}=-Cs(Y'')^2\eta^2,\qquad K_{i\phi}=K_{\phi\phi}=0.$$

The remaining linear metric is $\delta g_{a\phi}=-sY'\widehat h_a$;
the angular partner vanishes. Apply also

$$\Lambda=-\epsilon[bY+(q/2)\eta(uY'-Y)].$$

The resulting potential has $A_u=0$, $A_a=\epsilon\widehat a_aY$,
unchanged $A_\phi$, and **no quadratic gauge term**. These are exact
pullbacks of the specified affine family, not a truncation of a guessed
gauge-invariant action. The script verifies their full matrices, including
all derivatives of $C(x)$.

An independent operator check uses a general axisymmetric displacement
$\xi^\phi=-\psi(\tau,z,\theta)$. From the coframe

$$e^0=d\tau,\quad e^1=dz,\quad e^2=R\,d\theta,\quad
e^3=R\sin\theta(d\phi-\epsilon d\psi)$$

the script derives the Christoffel and spin connections, adds
$\delta A_\theta=-q\sin\theta\,\psi/2$, and verifies all 16 entries of

$$\delta D=[T,D_0],\qquad T=-i\psi(m-q\cos\theta/2).$$

This is independent of substituting invariants into the previously computed
kernel. In particular it retains the angular source omitted by that kernel.
Ordinary gauge covariance similarly follows from $\delta D=[i\lambda Y,D_0]$.
We use a covariant regulator for the vectorlike four-dimensional Dirac
determinant: there is no local vector-current or diffeomorphism anomaly.
The small periodic transformations preserve the spin structure. Both sides
of the finite-size subtraction undergo the same transformation.

## Ward identity with the background tadpole

Let $\Delta\Gamma_1[K]$ denote the first variation of the finite-size
effective action in the metric direction $K$, and let $\Delta\Gamma_{2,U}$
be the exact four-source quadratic functional with $H=b=0$. The convention
is that $\Gamma_2$ is the coefficient of $\epsilon^2$, not its second
derivative. Covariance gives

$$\boxed{\Delta\Gamma_{2,\rm affine}[a,b,h,H]
=\Delta\Gamma_{2,U}[\widehat a,\widehat h]+\Delta\Gamma_1[K].}$$

Proof: after the finite transformation the first-order sources are the
hatted ones and the only second-order source is $K$. Expand the transformed
effective action and use its equality to the original action. The first
variation in an axial or gauge direction vanishes by background symmetry;
its first variation in the diagonal $K$ directions does not. No stationarity
of the matter action with respect to the background metric was assumed.

For the homogeneous product background define its exact finite-size vacuum
energy $E(L,R,q)$ and

$$\mu_\tau=\frac{E}{2L},\qquad \mu_z=\frac{\partial_L E}{2},\qquad
\mu_\Omega=\frac{R\partial_R E}{4L}.$$

These signs follow by varying $\Gamma/T=\sqrt{g_{\tau\tau}}
E(L\sqrt{g_{zz}},R)$ and the sphere radius, avoiding stress-tensor sign
convention ambiguities. Rotational symmetry fixes the angular stress;
off-diagonal background stresses vanish. With angular average normalized
by $1/(4\pi)$,

$$\frac{\Delta\Gamma_1[K]}{TL}
=\langle\mu_\tau K_{\tau\tau}+\mu_zK_{zz}
             +\mu_\Omega g_S^{ij}K_{ij}\rangle_{x,S^2}.$$

The Legendre norms are

$$w_l=\langle s(Y')^2\rangle=\frac{l(l+1)}{2l+1},\qquad
v_l=\langle s^2(Y'')^2\rangle
 =\frac{l(l+1)(l-1)(l+2)}{2l+1}.$$

They follow from the norm of $P_l^1$ and $P_l^2$, respectively; the script
also integrates them directly for $l=2,\ldots,7$. Thus the required tadpole is

$$\frac{\Delta\Gamma_1[K]}{TL}
=-w_l\sum_{a=\tau,z}\mu_a
 \langle2h_a\partial_a\eta-C(\partial_a\eta)^2\rangle_x
 -v_l\mu_\Omega\langle\eta^2\rangle_x.$$

For $l=2$, $w_2=6/5$, $v_2=24/5$. The massless modes contribute
$\mu_\tau=-\pi q/(12L^2)$, $\mu_z=\pi q/(12L^2)$ and
$\mu_\Omega=0$. The massive modes are included below, not discarded.

For a pure linear gauge displacement
$(h_a,H,b,a_a)=(C\partial_a\chi,C\chi,q\chi/2,0)$ both invariants
vanish. Its **affine** quadratic coefficient is

$$-w_l C\sum_a\mu_a\langle(\partial_a\chi)^2\rangle_x
  -v_l\mu_\Omega\langle\chi^2\rangle_x.$$

The actual finite orbit has positive quadratic metric term
$K^{\rm orbit}_{ij}=Cs\,\partial_i(\chi Y')\partial_j(\chi Y')$.
Its tadpole cancels this expression exactly. This explains why imposing a
homogeneous null-vector Ward identity on the isolated affine matter Hessian
would fail: its Ward identity contains the background first variation.
On a solution of the **total** semiclassical equations, the total first
variation cancels among all sectors; that is a different assertion.

## Exact pressures and bounds on all omitted levels

Let $M_n=\sqrt{n(n+q)}/R$. One massive two-dimensional Dirac contributes

$$E_M=-\frac2\pi\int_0^\infty dp\,
 \log(1+e^{-L\sqrt{p^2+M^2}})
=-\frac{2M}{\pi}\sum_{w\ge1}\frac{(-1)^{w+1}}w K_1(wLM).$$

The first equality follows from the subtracted antiperiodic momentum sum,
or interchange of Euclidean space and time. Expanding the logarithm and
integrating gives the winding form and its massless limit $-\pi/(6L)$.
The two angular eigenvalue signs give

$$E=-\frac{\pi q}{6L}+2\sum_{n\ge1}(q+2n)E_{M_n}.$$

Differentiation under the exponentially convergent integrals yields, with
$e=\sqrt{p^2+M^2}$ and $n_F(x)=(e^x+1)^{-1}$,

$$\partial_LE_M=\frac2\pi\int_0^\infty e\,n_F(Le)\,dp>0,
\qquad R\partial_RE_M=-\frac{2LM^2}\pi\int_0^\infty\frac{n_F(Le)}e\,dp<0.$$

For levels $n>N$, absolute pressure tails are bounded by

$$|\mu_\tau^{>N}|\le\frac2{\pi L}\sum_{n>N}(q+2n)M_nK_1(LM_n),$$
$$|\mu_z^{>N}|\le\frac1\pi\sum_{n>N}(q+2n)M_n^2[K_0(LM_n)+K_2(LM_n)],$$
$$|\mu_\Omega^{>N}|\le\frac1\pi\sum_{n>N}(q+2n)M_n^2K_0(LM_n).$$

These follow from $n_F(x)\le e^{-x}$ and the positive Bessel integrals.
Set $n_0=N+1$. Using
$M_n\ge M_{n_0}+(n-n_0)/R$, $M_n\le(n+q/2)/R$ and
$K_\nu(x)\le e^{-x}\sqrt{\pi/(2x)}e^{\nu^2/(2x)}$
reduces every tail to polynomial geometric sums of degree at most three.
The script evaluates these infinite bounds, compares the integral and
winding forms, and checks derivatives by finite differences. Floating
evaluations are not interval certificates; the displayed inequalities are
the ordinary analytic bounds.

The identity above applies to the exact $\Gamma_{2,U}$. For the previously
computed mixed/LLL approximation, evaluate its massive remainder bound on
$(\widehat a,\widehat h)$ and retain the exact tadpole or these pressure
tail bounds. Completing a truncated kernel by invariants does **not** prove
the uncomputed massive kernel is zero. For complex Fourier amplitudes use
the Hermitian polarization of the real quadratic form; derivatives are
$iQ_a$, so angular and longitudinal sources need not have the same phase.
Longitudinal inversion makes the four-vector-source kernel even in $Q$;
reality then gives a real symmetric kernel before this substitution.

## Opposite static signs on paths with the same tangent

Retain $l=2$, $k=2\pi/L$, $f(z)=\sqrt2\cos kz$, $a=0$ and constant
amplitude $h\ne0$. Compare the earlier affine metric family to

$$ds^2=d\tau^2+dz^2+R^2d\theta^2
 +R^2\sin^2\theta\left[d\phi-\epsilon\frac{h f(z)Y'}{R^2}dz\right]^2,
\qquad A_\phi=qu/2.$$

This completed-shift family has exactly the same linear $g_{z\phi}$ and
gauge source. It adds
$\epsilon^2 K_{zz}^{\rm shift}=\epsilon^2 s(Y')^2h^2f^2/R^2$.
Therefore its static coefficient satisfies the **exact** relation

$$\mathcal Q_{\rm shift}(0,h)
=\mathcal Q_{\rm affine}(0,h)+\frac{w_2\mu_z}{R^2}h^2.$$

The added LLL pressure cancels the entire negative LLL contact.
The added massive pressure is positive. Write
$G_{TT}=\sum_{n=1,2}S_nb_n^2/\lambda_n$ in the notation of the static
note. Its rational formula satisfies $G_{TT}>1/4$ for every $q\ge3$:
the numerator of $G_{TT}-1/4$ has strictly positive coefficients after
$q=y+3$, with a positive denominator.

For $x=ML\ge32$ the previous mixed-loop error bound, divided by
$J_0=\pi/(3M^2L^2)$, is at most

$$\frac{8\pi^2}{x^2}+\frac{136\pi^4}{5x^4}
   +\frac{6x^2}{\pi^2}K_2(x)+24K_0(x)<\frac18.$$

To verify the last strict inequality without numerical quadrature, use
$\pi^2<10$, $\pi>3$ and $K_0,K_2<e^{-x}$ for $x\ge32$.
The power terms are less than $1/10$ at $32$ and decrease. The remaining
bound $(2x^2/3+24)e^{-x}$ decreases and is less than $1/1000$ at $32$,
using $e>2$. Hence $J>(7/8)J_0$ for both mixed levels when $L/R\ge16$.
The all-massive remainder proved in the static note is less than
$\pi h^2/(768R^2L^2)$, uniformly on this domain. Discarding the positive
added massive pressure gives

$$\mathcal Q_{\rm shift}(0,h)
>\left(\frac7{96}-\frac1{768}\right)\frac{\pi h^2}{R^2L^2}
=\frac{55\pi h^2}{768R^2L^2}
>\boxed{\frac{\pi h^2}{16R^2L^2}}>0.$$

Meanwhile the affine coefficient remains strictly negative by the previous
theorem. This proves path dependence of the isolated matter curvature, not
a disagreement between two descriptions of one finite geometry. The families
differ at second order. At a stationary point the **total** quadratic action
is independent of such second-order path choices; the missing Einstein,
Maxwell, infinite-line quantum and support contributions are essential to
that conclusion in the physical model.

## Acceptance and next physical question

The [script](../checks/casimir_cylinder_ward.py) and
[receipt](../evidence/casimir-cylinder-ward.json) retain the operator
commutator, exact finite pullback, angular norms, tadpole cancellation,
pressure checks and uniform sign-proof inequalities. No angular cutoff
enters either static sign theorem. The fixed product-background quantum
calculation and the kinematic identity valid for warped $C(x)$ are separate:
the latter does not supply a curved quantum stress tensor.

The next decision is whether the finite-flux response can be incorporated
into a controlled **total on-shell** throat/exterior model with a new
stability or obstruction consequence. More proofs about this off-shell sign
alone would not close that question. The source-proof/edition review and
primary/fallback selection remain open. No specialist review gate,
publication, GCP run or physical-shortcut claim is introduced here.

Remaining to-do list: complete novelty selection, derive the curved total
response and global physical consequence, verify suitable complete endpoints
on GCP and prepare the reproducible manuscript.
