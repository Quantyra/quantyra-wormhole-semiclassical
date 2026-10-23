# The full flow multiplier, with preserved bolt contact

2026-10-23, E003/S169. Ordinary argument before frozen execution.
The full stable-wormhole goal remains active and incomplete. The feasibility
audit stays closed. This implements a prerequisite of its selected global
forward/adjoint calculation; that physical solve is not executed here.

## 1. Identify the operator before discretizing it

Keep every physical parameter, paired domain, common subtraction and fixed
York map of current-geodesic-two-sided-remainder.md. Write P=iX for the
actual complete geodesic-flow generator and L=log(mu0 B). With base lift J,
the comparison quantum form is J* A(P) J, where

    A(w)=sum_n d_n integral r[-B_Mn(r)]g(rw)dr
                         -C F_eff w^4/(90pi^2).

The SM angular measure is dtheta, of mass 2pi. In a flat chart ONLY,

    K(p)=integral_0^(2pi) A(p cos(theta))dtheta.

Consequently K(P) is not the desired operator. Derive the inverse, rather
than change normalization and assume that repairs its shape. Setting
k(z)=K(sqrt(z)), Abel inversion gives

    A(w)=w^2/pi integral_0^1 k'(w^2(1-t^2))dt.           (1)

Indeed K(sqrt(z))=2 integral_0^z A(sqrt(y))/sqrt(y(z-y))dy;
the Abel inverse, k(0)=0, and y=z(1-t^2) give (1). For the retained
finite part F from current-all-angular-principal-operator.md,

    k'(z)=-C z/(64pi) integral_0^1
                 y^2 [F(Czy/4)+4L]/sqrt(1-y)dy.

This follows either from the convergent mass sum or integration by parts
in its known kappa formula. Substituting into (1), putting v=y(1-t^2),
and doing the remaining positive integral yields the SINGLE integral

    A(w)=C w^4 chi(w)/(90pi^2),
    chi(w)=(45/32) integral_0^1 v E(1-v)
                         [-F(Cw^2v/4)-4L]dv.           (2)

Here E(m)=integral_0^(pi/2) sqrt(1-m sin^2(theta))dtheta uses the
PARAMETER m, not the modulus. The inner identity is
integral_0^(pi/2)(1-(1-v)sin^2(theta))^(-3/2)dtheta=E(1-v)/v.
Equation (2) retains every angular mass and all flow frequencies.

The Bessel moments in the existing geodesic argument and Abel inversion
identify (2) with its original flow multiplier, not a scalar Laplacian
replacement. The original Bessel form has polynomial-log growth, so the
Abel transform is injective in this class. The identity is a scalar
spectral identity; applying it to P retains the ACTUAL geodesics and J.
It does not provide a global matrix, geodesic propagation, or a solution.

## 2. Exact contacts and a normalization trap

Let W(s)=integral_0^1 v^(s+1) E(1-v)dv. Euler integration and Gauss
summation give, near s=0,

    W(s)=(pi/2) Gamma(s+2)Gamma(s+3)
                        /[Gamma(s+7/2)Gamma(s+5/2)],
    W(0)=32/45, W'(0)/W(0)=mu=4log(2)-97/30.           (3)

The conventions/derivatives are [DLMF 19.2](https://dlmf.nist.gov/19.2)
and [19.4](https://dlmf.nist.gov/19.4), after m=k^2; the summation is
[15.4.20](https://dlmf.nist.gov/15.4.E20). These special-function identities
do not establish the physical paired-domain theorem.

Thus chi(0)=-F_eff and the exact high-frequency constant is

    chi(w)=4log|w|+2log C-4L+4log2-97/15+o(1).          (4)

Angular averaging adds 4 times the cos^4-weighted mean of log|cos|,
namely 4(7/12-log2), and recovers the existing -62/15 constant in
kappa. Even after normalizing the leading fourth coefficient, using
the old flat coefficient as chi would leave the nonzero contact error

    Delta_contact=4log2-7/3.

This is not an improved estimate or an equivalent operator. It matters
for the point adjoint: a fixed wrong fourth contact produces the already
documented non-L2 residual. Retain exact contacts symbolically, with
intervals enclosing their coefficients rather than redefining them.

## 3. Positive full-tower quadrature and exact moment correction

Use v_j=2^(-j/256), j=0,...,5120, delta=log2/256, v_min=2^-20.
Cell masses are EXACT elliptic differences. If m=1-v, put

    U(v)=32/45-a(m)E(m)-b(m)K(m),
    a(m)=(22+28m-18m^2)/45,
    b(m)=(-22+28m-6m^2)/45.

Then U'=vE(1-v), U(0)=0 and U(1)=32/45. Use normalized masses
lambda_j=(45/32)[U(v_right)-U(v_left)], and the tip mass
lambda_tip=(45/32)U(v_min). Their sum is exactly one. Place ordinary
nodes at each log midpoint and the tip node at log(v_min).

For f(s)=-F(exp(s))-4L the retained counting proof gives |f'|<=4,
|f''|<=4, uniformly in all frequencies. In log coordinates the density
is w(s)=exp(2s)E(1-exp(s)); 2<=w'/w<=5/2. The upper bound follows
by bounding sin^2(theta)/(1-m sin^2(theta))<=1/(1-m) in its defining
positive integral. On a length-delta cell the weighted mean's distance
from its midpoint is <=(5/2)exp(5delta/2)delta^2/12. Taylor's theorem
therefore bounds the normalized midpoint error by

    eps_bulk=delta^2[1/2+(5/6)exp(5delta/2)].            (5)

The tip approximation has one-sided error <=eps_tip=(45pi/64)v_min^2,
using E<=pi/2 and integral_0^b v log(b/v)dv=b^2/4.

A naive quadrature has a small but NONZERO error in (4). Repair it
exactly without changing either physical contact: set

    sigma=mu-sum_j lambda_j s_j,
    corrected_node_j=exp(s_j+sigma),
    chi_num(w)=sum_j lambda_j[-F(Cw^2 corrected_node_j/4)-4L].

All weights and sigma denote their exact formulas, not floating point
midpoints. Verify the corrected nodes remain in (0,1). Mass one keeps
chi_num(0)=-F_eff. Its exact log moment mu keeps (4). The extra uniform
error is at most 4|sigma|, so

    |chi-chi_num|<=eps=eps_bulk+eps_tip+4|sigma|.        (6)

For fixed finite nodes, the digamma remainder is O(w^-2). The integral's
remainder is also O(w^-2): split v at 4q^2/(Cw^2), use the retained
digamma estimate on the upper part and the logarithmic derivative bound
below; the v weight makes both terms integrable. Hence matching (4)
removes the unwanted fourth-contact mismatch. This proves the appropriate
asymptotic order, not a numerical residual bound on the global adjoint.

## 4. Convert quadrature error to the actual form budget

Do not compare its relative error with Qgeo alone: the outer Weyl
cancellation prevents that shortcut. Equation (6) gives directly

    |Delta Q|<=C eps ||X^2T||_SM^2/(90pi^2)
               <=C eps[J2/(120pi)+J1/(180pi)].          (7)

Use the joint actual coercivity and partition from the previous result:
l>=.99F_lambda E0+C S2/(100pi), CJ2<=2C S2+2(500)^4 CJ0,outer.
With its actual j1 and jo,

    eta_quad=max(5eps/3,
       eps[(500)^4 jo/(60pi)+j1/(180pi)]/(.99F_lambda)). (8)

This is an absolute V-to-V* bound by polarization and form closure,
including the singular bolt adjoint. Freeze the target eta_quad<4e-5
and eta_previous+eta_quad<.0022. This consumes the old rounding reserve;
it does not improve the physical continuum bound by normalization.
The unchanged <1.667e-9 bolt test and ideal EP<3.78031966038352e-7
then remain usable for the quadrature-defined comparison, IF its actual
global solve and residuals are computed. They are still unknown.

## 5. Execution scope and research decision

Require exact algebra/contact checks, positive interval weights, exact
mass/log moments, a frequency-uniform analytic error bound, the actual
form conversion (8), and independent real-integral diagnostics at low,
transition and very high frequencies. Diagnostic quadrature is not the
proof of the uniform bound. Preserve every snapshot and failed attempt.

This closes the spectral-multiplier implementation and its numerical
form-error allowance if the frozen targets pass. It does NOT close
geodesic propagation on the actual double, the global matrix, the signed
source, E/P/D0/Db0/Q0, bolt/exterior admission, or a nonlinear neighborhood.
The earlier paired-domain argument remains an input; this calculation
does not independently prove the quantum-source or boundary estimates.

Decision: continue the actual weak forward/adjoint method using (2),
the original J and the contact-complete point source. The single next
calculation is the global comparison pair and its two certified residuals
on the actual finite geometry, evaluating (not replacing) the geodesic
flow. If those residuals/norms fail the existing target, test the signed
directional remainder before optimizing generic source constants.
The multiplier test is executed only by its retained receipt; the full
physical feasibility test is still proposed and unexecuted.

Remaining to-do list: global comparison pair and physical bolt/exterior
certificate; X/admission, nonlinear closure and exact compatibility;
physical matching and causal realization, semiclassical validity/noise,
all-sector stability, charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
