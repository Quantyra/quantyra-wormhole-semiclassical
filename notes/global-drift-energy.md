# Exact energy and radial causal cones for a prescribed drift model

2026-10-10, Space S051/S052. A change of dynamical variables gives an exact
positive energy for the axial-metric/polar-gauge sector of a specified local
response model on a smooth static spherical background. The resulting
conditional theorem excludes finite-energy exponentially growing modes,
including oscillatory growth. It also relates one radial characteristic cone
to the source's radial null energy. **This is not a derivation of the response
from quantum field theory, an exact fermion-vacuum wormhole construction, or
a stability theorem for the physical MMP wormhole. Novelty remains unresolved.**

The starting response is the drift lift and Schwinger channel in Sadhukhan,
[arXiv:2610.09847v1](https://arxiv.org/html/2610.09847v1), sections 5-6.
Here we define its extension to arbitrary background coefficients as a
constitutive assumption and examine that mathematical model exactly. The
[tensor check](../checks/throat_geometry.py) independently reproduces the
source's corrected throat equations through first order in backreaction;
the [reduction check](../checks/throat_reduction.py) reproduces every radial
constraint and master definition at that order. Those reproduction results
are conditional on the source's quantum response.

## Defined model and assumptions

Use dimensionless magnetic-radius units, coordinates $(t,r,u,\varphi)$ with
$u=\cos\theta$, and

$$ds^2=-S(r)dt^2+S(r)^{-1}dr^2+C(r)d\Omega^2,\qquad C=R^2.$$

Take smooth $S,C>0$ on a connected interval, possibly all of $\mathbb R$.
Normalize the magnetic field to $F_{u\varphi}=2$ and
$\kappa/g_c^2=1/2$. This rescales perturbation amplitudes and fixes the magnetic
radius to one; it does not assert that physical magnetic charge four is in
the large-charge regime. Let $L=l(l+1)>2$, $k=L-2$, and let the Schwinger
coefficient $s$ be a positive constant. Results concern linear perturbation
amplitude in this parity sector only. Spherical symmetry permits using an
axisymmetric harmonic $Y(u)$ for the radial calculation.

Define the **scaled** residual source tensor by

$$\mathcal T_{\mu\nu}=G_{\mu\nu}-\tfrac12 T^{\rm EM}_{\mu\nu}
=\kappa T^Q_{\mu\nu}.$$

Thus the background solves Einstein's equation by definition. The existence
of a quantum state with this stress has not been established. Put
$P=\mathcal T_{\theta\theta}/C$ and
$K_a=\mathcal T_a{}^a-P$, with no sum for $a=t,r$. Prescribe

$$h_{a\varphi}=-h_a(1-u^2)Y_u,\qquad
\delta\mathcal T_{a\varphi}=P h_{a\varphi}+K_a U_a,\qquad
U_a=-(h_a+Cw_a/2)(1-u^2)Y_u.$$

Here $w_t=w_0,w_r=w_1$, and the perturbed gauge potential has components
$\delta A_t=w_0Y,\delta A_r=w_1Y$ in the chosen angular gauge.
All other independent source perturbations vanish. This excludes extra
axial source degrees of freedom and additional derivative or nonlocal
response terms.

For completeness the assumed conserved current is also explicit. Write
$D_a=g^{aa}K_a(h_a+Cw_a/2)$, without summing, and use the density
$\mathcal J^\mu=g_c^2\sqrt{-g}J^\mu$. Its polarization part is

$$\mathcal J^t_{\rm pol}=-CLD_tY,\quad
\mathcal J^r_{\rm pol}=-CLD_rY,\quad
\mathcal J^u_{\rm pol}=-(1-u^2)
\big[\partial_t(CD_t)+\partial_r(CD_r)\big]Y_u,\quad
\mathcal J^\varphi_{\rm pol}=0.$$

Add $\mathcal J^t_{\rm Sch}=sX'Y$,
$\mathcal J^r_{\rm Sch}=-s\dot XY$, and the anomaly equation

$$-\ddot X/S+(SX')'=E,\qquad E=\dot w_1-w_0'.$$

Both currents conserve charge separately. The prescribed stress/current
satisfy the linear source force-balance identity, as checked directly from
the tensors. These identities establish internal consistency; they do not
uniquely determine a quantum constitutive relation.

Define two dimensionless response factors:

$$A_t=1-\tfrac12 C^2K_t
=\frac{C(CS''-SC''+2)}4,$$

$$A_r=1-\tfrac12 C^2K_r
=\frac{C^2S''+CSC''+2C-S(C')^2}4.$$

The energy theorem additionally assumes **$A_t,A_r>0$ everywhere on its
domain**. No sign assumption is made about the source energy density itself.

## Exact reduction without dividing by frequency

Introduce the electric displacement and gravitational dual amplitude

$$\Phi=CE-sX,\quad
G=h_0'-\dot h_1-(C'/C)h_0,\quad
F_g=CG+2\Phi/L,\quad b_a=w_a+2h_a/C.$$

$F_g$ denotes an amplitude, not the electromagnetic field tensor. Stripping
the nonzero harmonic factors, the temporal and radial Einstein equations are

$$-\frac S{2C}(CG)'+\frac{k}{2C}h_0+\frac{A_t}{C}b_0=0,
\qquad
-\frac{\dot G}{2S}+\frac{k}{2C}h_1+\frac{A_r}{C}b_1=0.$$

The corresponding Maxwell equations are

$$\Phi'+\frac{LA_t}{S}b_0=0,\qquad
-\dot\Phi-LSA_r b_1=0.$$

They give the time-domain reconstruction

$$h_0=\frac S k F_g',\quad h_1=\frac{\dot F_g}{kS},\quad
w_0=-\frac{2h_0}C-\frac S{LA_t}\Phi',\quad
w_1=-\frac{2h_1}C-\frac{\dot\Phi}{LSA_r}.$$

The definitions of $G,E$ and the anomaly equation now close the system:

$$\frac{\ddot F_g}{SC}-\left(\frac S C F_g'\right)'
+\frac{k}{C^2}(F_g-2\Phi/L)=0,$$

$$\frac{\ddot\Phi}{SA_r}-\left(\frac S{A_t}\Phi'\right)'
+\left(\frac L C+\frac4{C^2}\right)\Phi
-\frac{2L}{C^2}F_g+\frac{Ls}C X=0,$$

$$\frac{\ddot X}S-(SX')'+\frac{\Phi+sX}C=0.$$

The angular Einstein equation $\dot h_0-S(Sh_1)'=0$ is automatic under
reconstruction. The remaining Maxwell equation follows from charge
conservation; the other tensor components vanish. Conversely, the displayed
system and reconstruction solve every radial equation and both field
definitions. The verification retains the off-equation residuals of both
Einstein equations, rather than checking only the Maxwell constraints.
No step divides by a temporal frequency. This removes the nonzero-frequency
restriction of the original master-variable reconstruction.

## Positive energy and its precise consequence

For $z=(F_g,\Phi,X)^T$ set

$$K=\operatorname{diag}\left(\frac{L^2}{kSC},\frac1{SA_r},\frac{Ls}S\right),
\qquad
B=\operatorname{diag}\left(\frac{L^2S}{kC},\frac S{A_t},LsS\right).$$

Multiplying the three equations by $L^2/k,1,Ls$, respectively, gives
$K\ddot z-(Bz')'+Vz=0$, where $V$ is symmetric and

$$z^TVz=\frac{(LF_g-2\Phi)^2}{C^2}
+\frac L C(\Phi+sX)^2.$$

Consequently

$$e=\tfrac12\left[\dot z^TK\dot z+z'^TBz'
+\frac{(LF_g-2\Phi)^2}{C^2}+\frac L C(\Phi+sX)^2\right],
\qquad j=-\dot z^TBz',$$

satisfy $\partial_t e+\partial_r j=0$. Proof: multiply the equation by
$\dot z^T$ and use that all coefficients are static and $K,B,V$ are symmetric.
The two potential squares prove $V\ge0$; the sign assumptions give
$K,B>0$. Equivalently these are the Euler-Lagrange equations of the quadratic
action with density $(\dot z^TK\dot z-z'^TBz'-z^TVz)/2$.

**Conditional theorem.** Consider sufficiently regular solutions with finite
total energy and boundary conditions for which the boundary flux exists and
$j(b)-j(a)\ge0$. Then the total energy is nonincreasing. It is conserved for
zero total flux, including reflecting boundaries and admissible vanishing
flux at both infinite ends. There is no nonzero finite-energy normal mode
$z=e^{pt}v(r)$ with $\operatorname{Re}p>0$ satisfying those same conditions.

For the mode assertion extend energy to complex fields using Hermitian
products and real flux. Its value is $e^{2\operatorname{Re}p\,t}$ times a
nonnegative constant. For a nonzero mode with $p\ne0$ the kinetic contribution
$|p|^2\int v^*Kv/2$ is positive, contradicting nonincrease. This covers
oscillatory growing modes as well as purely exponential ones. It presupposes
the stated regularity, integrability and boundary conditions; finite energy
alone is not a substitute for verifying a physical matching condition.

The argument works on a single smooth two-ended domain. If coefficients are
glued piecewise, continuity of $z$ and $Bz'$ would cancel interface flux, but
the physical junction conditions must establish these requirements. No such
matching for the actual MMP exterior is claimed.

There is no uniform spectral-gap or pointwise-boundedness claim. On domains
that permit it, the constant zero-energy kernel is
$z=c(-2s/L,-s,1)$. It reconstructs zero metric/gauge perturbations and zero
current, but a full microscopic equivalence of constant boson offsets is not
proved. Secular solutions and low multipoles $l=0,1$ require their own
analysis. This is a linear energy statement, not finite-amplitude stability.

## Radial cones and the null-energy relation

The principal symbol has coordinate radial speeds
$S,S\sqrt{A_r/A_t},S$. Relative to the metric light speed, these are
$1,\sqrt{A_r/A_t},1$. For the radial metric-null vector
$n=S^{-1/2}\partial_t+S^{1/2}\partial_r$ one finds

$$\mathcal N=\kappa T_Q(n,n)
=-\frac{S[2CC''-(C')^2]}{2C^2}=-\frac{2SR''}R,$$

$$A_r-A_t=-\tfrac12 C^2\mathcal N=SR^3R'',\qquad
v_\Phi^2-1=-\frac{C^2\mathcal N}{2A_t}.$$

The magnetic Maxwell background has zero radial null contraction. At a
strict smooth flare, $R'=0,R''>0$, positive $S,A_t$ therefore imply
$v_\Phi>1$ in this exact local model. The displacement channel has a nonzero
electric field; it is not a discarded constraint mode.

This gives a conditional incompatibility: a model with this exact response
at all frequencies, a strict flare and the positive-energy signs cannot
also have all radial characteristics inside the metric cone. It is **not**
a no-go theorem for semiclassical wormholes. Low-energy quantum response
can be dispersive and nonlocal, so its finite derivative approximation does
not determine the microscopic front velocity. No outside travel-time
advantage or controllable faster-than-light signal has been established.
Nor does cone widening alone imply a time machine: when the assumed static
$t$ is global, its gradient is timelike for all these radial cones. It is a
common time function, so the radial effective problem is stably causal.

These distinctions have important precedents. Liberati, Sonego and Visser,
[gr-qc/0107091v2](https://arxiv.org/html/gr-qc/0107091v2), section 3.2 and
Appendix B, already discuss widened effective cones, the null-energy sign
and stable causality in the Casimir setting. Shore,
[gr-qc/0302116v1](https://arxiv.org/html/gr-qc/0302116v1), section 3.2,
relates wavefronts to the high-frequency limit. Hollowood and Shore,
[1006.1238v1](https://arxiv.org/html/1006.1238v1), sections 1-4, explain why
low-frequency curvature corrections cannot simply be extrapolated to that
limit. Their QED calculation is not a derivation for the present massless,
magnetized wormhole source. The broad energy/cone idea is not our novelty claim.

## Nonempty hypotheses and baseline limits

A smooth mathematical example is, for any $a>0$,

$$v=\sqrt{r^2+a^2},\quad R=1+v,\quad C=R^2,\quad
S=(1-1/R)^2.$$

It has two asymptotic ends and a horizonless minimum radius $1+a$. Direct
substitution gives
$A_t=1+a^2/v>1$ and $A_r=1+a^2+2a^2/v>1$. Thus the exact hypotheses are
nonempty. Its residual source and drift response are prescribed, not obtained
from a fermion state. In fact the displacement's asymptotic relative speed
approaches $\sqrt{1+a^2}$, so even asymptotic Maxwell propagation should not be
silently attributed to this artificial constitutive example. No formation,
same-universe embedding, quantum support or novelty of this geometry is claimed.

For the electrovacuum mouth $C=r^2$, $S=1-2M/r+1/r^2$, both response factors
are exactly one. For the first-order MMP throat, with $f=1+r^2$ and the
background corrections of the [baseline note](axial-baseline-and-throat-energy.md),

$$A_t=1+\alpha\frac{2r^2-1}{6f}+O(\alpha^2),\qquad
A_r=1+\alpha\frac{2r^2+5}{6f}+O(\alpha^2).$$

These have positive leading margins on the domain of the expansion. The
source's throat expansion requires $|\alpha r|\ll1$, while its lowest-Landau-level
description restricts angular momenta below the magnetic scale
($l\ll\sqrt q$) and local frequencies below the Landau gap. The screening
parameter is retained without expansion. Mathematical statements for every
$L>2$ do not remove these physical cutoffs.

The
exact theorem controls arbitrary smooth backgrounds **within the defined
response class**. It does not control omitted response terms outside that
class or prove existence of a global quantum background. The earlier
first-order square completion and its indefinite second-order truncation
remain valid observations; they were not exact positivity proofs. The present
variables give an exact structure without making that truncation.

## Contribution gate and verification

The energy technique is established; compare Jansen, Rostworowski and
Rutkowski, [1909.04049v2](https://arxiv.org/abs/1909.04049v2), sections 2 and 4,
and Kimura and Tanaka, [1809.00795v3](https://arxiv.org/abs/1809.00795v3),
section II.1. The former's ordinary four-dimensional scalar has different
angular dynamics from the flux-line charge field here. Dutta Roy,
[2110.05019v3](https://arxiv.org/html/2110.05019v3), Appendices A-B, instead
uses vanishing axial source components in its tetrad and a different
effective matter model. Neither comparison by itself establishes novelty.

The candidate contribution is the full constrained reduction, positive
energy and explicit radial cone/null-energy relation for this specified
three-channel response, potentially followed by a controlled quantum
completion. General S-deformation positivity, vacuum cone widening and
global time-function arguments are already known. Whether the particular
result is new and physically consequential remains under review.

The [exact check](../checks/general_drift_energy.py) derives the background
and linear tensors, verifies all source/constraint identities, checks the
potential Hessian and local energy law, and reproduces the response factors
in both baseline limits and the smooth example. Its
[receipt](../evidence/general-drift-energy.json) pins code and dependency
hashes. Symbolic frequency identities are polynomial identities with static
coefficients and translate to time derivatives; the energy identity is
checked directly in the time domain. These are algebraic verification
receipts, not a PDE simulation or a Lean certificate. No Lean run has started.

Remaining to-do list: finish the close-source proof and edition comparison;
derive or reject the quantum drift approximation with its omitted terms;
assess a global physical background and novelty; then select the primary
and fallback before GCP certification and manuscript preparation.
