# Axial baseline reproduction and first order throat energy

2026-10-10, Space S051/S052. The classical mouth equations and their
master-variable reconstruction reproduce exactly. The corrected throat
matrix also admits a positive symmetrizing weight and an energy factorization
**through first order** in the background correction. These calculations
advance the global stability feasibility test; they do not establish global
stability, the quantum closure, or novelty.

The baseline is Sadhukhan,
[arXiv:2610.09847v1](https://arxiv.org/html/2610.09847v1), especially
sections 3, 5-7 and Appendices A-D. All its main text and appendices have
now been read. The cited quantum-response dependencies and global validity
argument remain under audit. The energy method, including a positive
frequency weight, is established in Kimura and Tanaka,
[arXiv:1809.00795v3](https://arxiv.org/html/1809.00795v3), section II.1
and equation (21). No new general stability method is claimed.

## Mouth reproduction

Use the notation of the [mouth energy argument](mouth-energy-factorization.md):
$e=r_e$, $M_0=GM$, $L=l(l+1)$, $k=L-2$, $s=m_S^2$, and
$f=1-2M_0/r+e^2/r^2$. For $l\ge2$, take time dependence $\exp(pt)$
and an axisymmetric spherical harmonic $Y(u)$, $u=\cos\theta$.
The geometry check starts from

$$g_{\mu\nu}=\operatorname{diag}(-f,f^{-1},r^2/(1-u^2),r^2(1-u^2)),
\quad h_{a\phi}=-h_a(1-u^2)\partial_uY,$$

with $a=t,r$. It constructs the connections, curvature, Maxwell stress
and Maxwell divergence directly. Setting the magnetic normalization to
$q=4$ makes the rescaled gauge amplitudes equal to the perturbation
amplitudes, $F_{u\phi}=2$ and $\kappa/g_c^2=e^2/2$.
This is an algebraic normalization test, not a claim that $q=4$ lies in
the large-flux physical regime.

The assumed quantum terms are zero linear quantum stress in the mouth and
the prescribed charge current density

$$g_c^2\sqrt{-g}\,J^t=sX'Y,\qquad
g_c^2\sqrt{-g}\,J^r=-spXY.$$

The six nontrivial Einstein-Maxwell components reproduce equations
(A.1)-(A.3), (A.6), (A.13)-(A.14), with explicit harmonic and denominator
factors retained in the [geometry receipt](../evidence/mouth-geometry.json).
The background Einstein and Maxwell equations, zero linear volume
variation, the other seven independent Einstein components, Maxwell's
azimuthal component and charge conservation also vanish exactly.

Separately, the [reduction check](../checks/mouth_reduction.py) reconstructs
the radial fields from $z=(g,d,\xi)$ solving the published master system:

$$h_1=\frac r f g,\quad h_0=\frac f p(fh_1)',\quad
E=\frac{L}{per^2}(d+s\xi),\quad X=\frac L{pe}\xi,$$

$$F_D=\frac L{pe}d,\qquad
w_1=-\frac{2h_1}{r^2}-\frac p{Lf}F_D,\qquad
w_0=-\frac{2h_0}{r^2}-\frac f L F_D'.$$

All eight residuals vanish: the six field equations, anomaly equation
(A.15), and electric-field definition. Here primes mean $d/dr$.
The inverse reconstruction requires $p\ne0$ and the stated nonzero
denominators. It therefore covers a putative exponentially growing mode,
but does not settle the static sector, secular solutions, or boundedness of
all original perturbations. The two checks verify classical geometry and
algebra conditional on the given source. Neither derives the quantum
effective action.

## First order throat structure

Here $\rho$ and $x$ are the source's dimensionless throat coordinates,
$f=1+\rho^2$ and $\partial_x=(f+\alpha\gamma)\partial_\rho$.
This $f$ differs from the mouth function. In the source's variables after
removing the first derivative, write

$$-z_{xx}+(fN+\alpha V_a)z
=\Omega^2(I+\alpha V_b)z+O(\alpha^2),$$

$$N=\begin{pmatrix}k&-2&0\\-2k&L+4&Ls\\0&1&s\end{pmatrix},
\qquad V_b=\operatorname{diag}(0,-1/f,0).$$

$V_a$ is Appendix D, equations (D.9)-(D.15). The new
[check](../checks/throat_energy.py) independently transcribes both the
pre-transformation coefficients and $V_a$, then reproduces their conversion,
including the matrix commutator in equation (7.7).

Define constant matrices

$$H_0=\operatorname{diag}(k,1,Ls),\qquad
H_1=\operatorname{diag}(-k/6,0,-Ls/6),\qquad H_\alpha=H_0+\alpha H_1.$$

Direct substitution gives

$$H_0V_a+H_1fN=(H_0V_a+H_1fN)^T.$$

Thus the matrix asymmetry is removable through first order. Both
$H_\alpha$ and
$H_\alpha(I+\alpha V_b)=\operatorname{diag}(k(1-\alpha/6),
1-\alpha/f,Ls(1-\alpha/6))$ are positive for $L>2$, $s>0$,
$0\le\alpha<1$. This does not say that $H_\alpha(fN+\alpha V_a)$
is exactly symmetric at finite $\alpha$: its second-order term has been
discarded in this assertion.

Choose the diagonal deformation $S=\alpha S_1$ with

$$S_1=\operatorname{diag}\left(
\frac{\rho+f\arctan\rho}{2},\frac{\rho}{2f},0\right).$$

The gravitational entry equals $R_x/R$ to first order for the source's
corrected radius. The electric entry is essential to the following
cancellation. Set

$$A_0=H_0fN,\qquad
A_1=H_0V_a+H_1fN+H_0f\partial_\rho S_1.$$

For the neutral vector $v_0=(-2s/k,-s,1)^T$,
$Nv_0=0$ and $v_0^TA_1v_0=0$. The check also reproduces the undeformed
neutral potential and frequency coefficient in equation (7.26).

A stronger identity exposes all three channels. Write

$$\begin{aligned}
T&=12\rho^3\arctan\rho+17\rho^2+24\rho\arctan\rho-6\log f+9,\\
a&=-T/(12f),\\
b&=-\big(18\rho^3\arctan\rho+22\rho^2+30\rho\arctan\rho
-6\log f+17\big)/(12f),\\
c&=a+\rho^2/(4f),\\
B_0&=\begin{pmatrix}k&-2&0\\0&1&s\end{pmatrix},\qquad
B_1=\begin{pmatrix}ka&-2b&0\\0&c&sa\end{pmatrix},\qquad J=\operatorname{diag}(1,L).
\end{aligned}$$

Then exactly

$$A_0=fB_0^TJB_0,\qquad
A_1=f(B_0^TJB_1+B_1^TJB_0).$$

Let $Q$ have derivative coefficient $H_\alpha$ and potential coefficient
$H_0fN+\alpha(H_0V_a+H_1fN)$. Consequently this first-order weighted energy,
including its endpoint term, has the formal expansion

$$Q[z]=\int\left[(z_x+Sz)^*H_\alpha(z_x+Sz)
+f\left|(B_0+\alpha B_1)z\right|_J^2\right]dx
-[z^*H_\alpha Sz]_a^b+O(\alpha^2).$$

This is equality of the zeroth and first coefficients on a fixed compact
throat interval. The remainder includes the square completion, $S^2$,
the corrected coordinate derivative, and the omitted physical operator.
No uniform remainder estimate on the matched global geometry is supplied.

## Why this is not yet a stability theorem

The squares cannot be retained while silently dropping their second-order
cost. In fact the exact determinant of the truncated deformed matrix at
$\rho=0,L=6,s=1$ is

$$\det(A_0+\alpha A_1)=\frac{256}{3}\alpha^2(3\alpha-2),$$

which is negative for $0<\alpha<2/3$. This is pointwise indefiniteness
at second order, not a proof of a negative spectral eigenvalue. It explicitly
rules out claiming exact positivity from the first-order identity alone.
The neutral channel has no uniform global spectral gap with which to
absorb an arbitrary error. The throat expansion also requires
$|\alpha\rho|\ll1$; it cannot be extended to either infinity.

The next scientific test is to derive a global energy domain and a relative
bound for interface and omitted terms, including the neutral channel.
Throat and mouth fields, weights, tortoise coordinates and endpoint signs
must be transformed consistently. Merely joining two positive boxes would
not establish these physical conditions. Reproducing the throat's tensors
with the stated quantum response remains part of that test.

## Quantum closure question for finite disturbances

The source's local drift stress is an additional matter-response argument,
not a consequence of charge conservation alone. Its matching test is
insensitive to alternative drift closures at the leading overlap order.
A first-principles Landau-level derivation and error hierarchy remain relevant.

One geometric distinction matters for a nonlinear extension. The source's
tilted invariant planes have tangents
$e_a=\partial_a+\epsilon v_a\partial_\phi$ with
$v_a=2\beta_aY'(\theta)/(q\sin\theta)$. For axisymmetric amplitudes,

$$[e_\tau,e_\rho]=
\frac{2\epsilon E Y'(\theta)}{q\sin\theta}\partial_\phi+O(\epsilon^2),
\qquad E=\partial_\tau\beta_\rho-\partial_\rho\beta_\tau.$$

For nonzero coefficient this vector is not in the two-plane: its $\tau$
and $\rho$ components vanish. Thus the perturbed planes need not form
two-dimensional surfaces. The field-sheet theorem in Gralla and Jacobson,
[arXiv:1401.6159v4](https://arxiv.org/html/1401.6159v4), sections 3.2.3-3.2.5,
requires a closed, degenerate electromagnetic two-form. The electric
perturbation need not satisfy that degeneracy condition.

This calculation does **not** refute the source's linear stress lift:
an algebraic local tensor need not define integrable sheets, and a linear
current can be evaluated on unperturbed sheets. It identifies a condition
that must be addressed before treating tilted planes as independent quantum
worldsheets at finite amplitude. It is not a novelty claim or a wormhole
no-go theorem.

Remaining to-do list: complete the quantum-response and novelty audit,
reproduce the corrected throat tensors, and establish or reject global
matching and a useful error bound before selecting formalization or a manuscript.
