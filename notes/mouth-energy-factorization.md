# Exact energy bound for the coupled mouth operator

2026-10-10. Research calculation under Space S051/S052. **Ordinary argument
and computational verification; novelty and physical-model audit pending.**

## Question and prior result

Sadhukhan's [axial-metric/polar-gauge paper, version 1](https://arxiv.org/html/2610.09847v1)
gives a three-channel mouth operator in equations (5.20)-(5.22), with Schwinger
screening parameter $m_S^2$. Section 10.2 reports no negative eigenvalues for
specified finite grids and parameters. A negative eigenvalue of the
pointwise potential does not imply a negative eigenvalue of the differential
operator. The following factorization addresses that distinction analytically.

The method is the established matrix S-deformation/energy method, not a new
proof technique. Compare Jansen, Rostworowski and Rutkowski,
[arXiv:1909.04049v2](https://arxiv.org/abs/1909.04049v2), section 4, and its
predecessors. Whether this application is already a direct instance of an
existing physical theorem is still being checked. A simple application alone
may be insufficient for the program's contribution requirement.

The 2019 paper's equations (40) and (43) already supply the energy method and
the gravitational deformation $f/r$ in four dimensions. Its action (2)
describes a conventional Einstein-Maxwell-scalar theory. The extra charge
field here instead comes from fermions confined to flux worldsheets, with
anomaly response and no ordinary four-dimensional scalar angular kinetic
term. This distinction prevents an immediate identification of the physical
models, but does not make the elementary factorization a new general method.

## Operator and assumptions

Write $M_0=GM$, $e=r_e>0$, $L=l(l+1)$, $k=L-2>0$, and $s=m_S^2>0$.
Here $s$ is the squared screening parameter, not the deformation function.
Let

$$f(r)=1-\frac{2M_0}{r}+\frac{e^2}{r^2},\qquad
\frac{dr}{dx}=f(r),\qquad C(r)=\frac{f(r)}{r^3}.$$

On an interval with $r>0$ and $f>0$, the source operator is
$\mathcal H=-\partial_x^2+V$, where

$$V=C\begin{pmatrix}
Lr-6M_0+4e^2/r&-2e&0\\
-2ek&Lr+4e^2/r&Lsr\\
0&r&sr
\end{pmatrix}.$$

Its components are $z=(g,d,\xi)^T=(\Psi_g,\Psi_D,\Xi)^T$ in the source
conventions. This $g$ is a master amplitude, not a gauge coupling; $d$ is a
master amplitude, not the mouth separation. The positive constant weight
$H=\operatorname{diag}(k,1,Ls)$ makes $HV$ symmetric. It is equivalent to the
source's constant symmetrizer. The zero-screening limit requires separate
variables because this weight degenerates at $s=0$.

For the physical multipoles $l\ge2$, $k>0$. The algebra below holds more
generally for $L>2$. It is exact for this operator at any $M_0$ on any such
interval. This does not promote that operator to an exact finite-backreaction
description of the physical MMP mouth. The source derives it within specified
lowest-Landau-level, mouth, angular and perturbative approximations.

## Factorization

Choose $S=\operatorname{diag}(f/r,0,0)$ and
$D=\partial_x+S$. Since $H$ commutes with $S$,

$$W=V+\partial_xS-S^2.$$

Only the gravitational diagonal changes. Direct differentiation gives

$$f\frac{d}{dr}\left(\frac f r\right)-\frac{f^2}{r^2}
=C\left(-2r+6M_0-\frac{4e^2}{r}\right).$$

Consequently

$$z^*HWz=C\left[\frac{|krg-2ed|^2}{r}
+Lr|d+s\xi|^2\right]\ge0.$$

This identity retains all three channels and every positive screening
parameter. It does not use elimination of a heavy mode, a large-screening
limit, extremality, an angular truncation of a scan, or positivity of $V$.
The latter can fail even though the differential operator is nonnegative.

For a continuously differentiable complex test field on $[a,b]$ in $x$,
integration by parts yields

$$\begin{aligned}
Q[z]&=\int_a^b\left(z'^*Hz'+z^*HVz\right)dx\\
&=\int_a^b\left((Dz)^*H(Dz)+
C\left[\frac{|krg-2ed|^2}{r}+Lr|d+s\xi|^2\right]\right)dx
-\left[\frac{kf}{r}|g|^2\right]_a^b.
\end{aligned}$$

Compact support or Dirichlet endpoints make the last term zero. The quadratic
form is then nonnegative. Its Friedrichs realization in the weighted Hilbert
space has no negative spectrum, hence no square-integrable exponentially
growing modes for the corresponding self-adjoint wave problem. This operator
consequence uses the usual quadratic-form construction, which is not yet
formalized here. Complex fields follow by applying the real identity to
their real and imaginary parts.

For a right-hand mouth $r\in[r_0,\infty)$, with the outer endpoint term
vanishing and $r_0>0$, $f(r_0)>0$, the retained term is **positive**:

$$Q[z]\ge\frac{kf(r_0)}{r_0}|g(r_0)|^2.$$

This statement concerns the energy integral; identifying it with
$\langle z,\mathcal H z\rangle$ at a free inner boundary also requires the
boundary term from $-\partial_x^2$ and a specified self-adjoint boundary or
matching condition. An arbitrary Robin condition can change the conclusion.
For a left-hand mouth use the global coordinate orientation consistently;
reversing $dr/dx$ reverses the deformation sign.

An abstract pair of mouths connected to a nonnegative throat energy can
inherit positivity if a common positive energy weight and compatible
matching remove the operator's interface terms. Verifying these hypotheses
for the physical throat is an unresolved step. Numerical stability of two
separate finite boxes would not by itself verify them.

## Verification and scope

[The check](../checks/mouth_energy.py) differentiates $f$, verifies the
weighted matrix and exact Gram identity symbolically, and checks the stated
determinant of the source potential. It then independently integrates both
energy expressions in the radial coordinate, including cases with a nonzero
inner boundary value. A negative pointwise-potential example confirms that
the check is addressing the indefinite-potential problem. Numerical agreement
checks transcription and boundary signs; the algebra supplies the universal
identity. [The retained output](../evidence/mouth-energy.json) records versions
and the exact check-script hash.

No claim is made here about low multipoles outside the source reduction,
nonlinear disturbances, quantum-state closure, corrected full-throat spectra,
the full three-dimensional exterior, formation, or a travel-time shortcut.
The [baseline audit](axial-baseline-and-throat-energy.md) now independently
reproduces the classical mouth tensors and the nonzero-frequency master
reconstruction, conditional on the prescribed quantum current. The physical
quantum approximation and corrected throat still require their own audit.
This bound also supplies no uniform spectral gap: a massless channel can
approach zero energy. It cannot automatically control arbitrary small omitted
terms. The observation that $M_0$ cancels from the deformed bracket only covers
mass changes inside the stated operator, not all corrections at that order.

Remaining to-do list: complete the novelty and physical-source audit; reproduce
the corrected throat calculation; establish or reject the matching and correction
control needed for a global stability result; certify a suitable complete
endpoint on GCP only after that argument and its payoff are clear.
