# Finite-frequency RN dipole and prescribed throat matching

2026-10-11, S053. The physical RN dipole master and a finite-frequency
metric/Maxwell reconstruction are derived directly in the magnetic gauge.
An explicit annulus estimate controls frequency corrections without
imposing a black-hole horizon condition. Value and flux matching to the
leading product throat preserves the previous coefficient 3/2 at leading
order for a prescribed two-mouth classical master problem.

This is progress toward physical matching. It does not determine the
actual binary's outer drive, construct a smooth four-dimensional junction,
fix boundary frames/holonomy, include the other modes, or derive the
semiclassical state and source response. Those limitations are distinct
from the finite-frequency classical calculation completed here.

## 1. A physical master in the magnetic RN gauge

Set the geometric extremal charge P=1 temporarily; physical frequency
Omega becomes w=Omega P. Let r be radius in units P, z=r-1,
v=3r-2=1+3z, f=z^2/r^2 and c=cos(theta). The background is

    g0=-f dt^2+dr^2/f+r^2 dOmega^2,    F0=-dc wedge dphi.

Use a common exp(-i w t) factor and the metric perturbation

    h_tt=-f T c,    h_rr=H c/f,    h_tr=J c,
    h_aA=h_AB=0,    delta A=a(r) sin(theta)^2 dphi.       (1)

Here H is a radial metric amplitude, not the overlap radius R below.
Direct Einstein constraints give

    J=i w r H,
    r z^2 H'+(2r^2-1)H=4a.                             (2)

Neither identity divides by w. Define

    Z=a-z^2 H/(2v),    b=a-Z.                           (3)

The residual coordinate freedom preserving the angular and mixed-angular
gauge has xi^r=chi c, xi^c=(chi/r)(1-c^2),
xi^t=-i w r chi c/f and chi'=(1-1/f)chi/r.
It changes a by chi/r and H by 2v chi/(r z^2), so Z is unchanged.
The residual homogeneous radial profile is proportional to
`exp(1/z)/z` in a. Boundary values of this coordinate change still matter
when comparing clocks and charged-spinor frames; (3) is a bulk invariant.

Equations (2) become the first-order reconstruction equation

    b'+(r/z^2)b=2 Z/(r v).                             (4)

Together with (4), reconstruct all the amplitudes by

    a=Z+b,
    H=2v b/z^2,
    J=i w r H,
    T=(w^2 r^4-z)H/v +4z^2 Z'/(r v)
                         +[4/v+8z^2/(r^2 v^2)] Z.     (5)

The Maxwell equation, before imposing an equation on Z, reduces exactly to

    (f Z')'+(w^2/f-V/f)Z=0,
    V=2z^2(9r^4-18r^2+16r-4)/(r^6 v^2).               (6)

Equivalently, `(d_rstar^2+w^2)Z=VZ`, with d_rstar=f d_r.
After (6), every linear Einstein-Maxwell component for (1) vanishes.
The direct contraction includes the time/radial component and induced
electric field. Thus (4)-(6) provide an actual finite-frequency classical
reconstruction, rather than a time-dependent replacement in a static field.

The potential is precisely the extremal specialization of Burko's
[polar dipole master](https://arxiv.org/abs/gr-qc/9507012v1), equations
(40), (41), (87). The magnetic calculation above independently obtains
the field identification and verifies the tensor equations; it does not
assume the paper's electric tetrad reconstruction in our boundary gauge.
The RN dipole dynamics themselves are established prior work.

The exact static growing solution of (6) is

    Z0=r z^2/v.                                        (7)

For the force-balanced mouth in the [previous note](dipole-throat-matching.md),
`a=epsilon z^2(1/2+I)`, `H=2epsilon(vI-1)` and epsilon=BP.
Substitution into (3) cancels I exactly:

    Z=(3/2)epsilon Z0.                                 (8)

Near z=0, Z0/z^2 tends to one. This identifies the coefficient 3/2
directly in the physical master; the earlier nonlocal gauge integral
is confined to reconstruction. On a finite annulus, (4) requires a
coordinate-gauge boundary datum. Nothing here sets that datum from the
binary or discards its boundary interpretation.

## 2. An explicit finite-annulus frequency bound

Take 0<s<L<infinity and write Z=Z0 W. Equation (6) is exactly

    (p W')'+w^2 k W=0,
    p=z^6/(1+3z)^2,    k=z^2(1+z)^4/(1+3z)^2.          (9)

For initial values W(s)=A and p(s)W'(s)=M, it is equivalent to

    W(z)=A+M J_s(z)-w^2 (T_s W)(z),
    J_s(z)=int_s^z du/p(u),
    (T_s W)(z)=int_s^z du/p(u) int_s^u k(v)W(v)dv.      (10)

Since k(v)<=v^2(1+v)^2, the positive integral operator obeys
`||T_s||_infinity <= K(s,L)`, where

    K=(s^-2-L^-2)/6 +(5/2)(s^-1-L^-1)
      +(31/5)log(L/s)+(57/10)(L-s)+(9/10)(L^2-s^2).    (11)

To verify (11), replace the inner integral by
`u^3/3+u^4/2+u^5/5` and multiply by p(u)^-1. The resulting integrand is
`1/(3u^3)+5/(2u^2)+31/(5u)+57/10+9u/5`, whose integral is (11).
This bound explicitly includes the redshift near the neck and the large
outer radius; small w alone is not a uniform estimate on an unbounded domain.

Let eta=w^2 K<1 and F=A+M J_s. The Neumann series in the supremum norm
gives the exact, finite-interval estimates

    ||W|| <= ||F||/(1-eta),
    ||W-F|| <= eta ||F||/(1-eta).                       (12)

For A=1, put zeta=|M| J_s(L). Then
`||W-1|| <= (zeta+eta)/(1-eta)`. If this last bound is beta<1,
normalizing the outer value to one changes the bound to
`||W/W(L)-1|| <= 2beta/(1-beta)`.

For the special zero-load case M=0 and real w, positivity improves this:
`1-eta <= W <= 1`, hence `||W/W(L)-1|| <= eta/(1-eta)`.
Indeed, up to a first zero, (9) makes W nonincreasing, and (10) would give
W>=1-eta>0 at that zero, a contradiction. We do not impose this special
zero load on the physical wormhole; the general estimate (12) retains M.

No horizon regularity or ingoing condition was used in (9)-(12).

## 3. Two-mouth matching in the prescribed classical master problem

Let lambda=P/ell, sigma=Omega ell and w=lambda sigma. Retain the exact
product-throat scalar equation and solution from the previous note,
with its two growing coefficients D=(D_plus,D_minus). At both ends use
z=lambda |rho| and an overlap radius R, so s=lambda R.

Define the piecewise classical master problem explicitly: use the product
equation on |rho|<=R, the RN equation on s<=z<=L at each mouth, and match
the master value and outward flux. In the common z coordinate, the
product and RN flux coefficients are respectively lambda^2+z^2 and f.
Thus at the positive end the initial RN data are

    Z(s)=u(R),
    f(s)Z'(s)=(lambda^2+s^2)u'(R)/lambda.               (13)

Use the outward derivative at the negative end. This is a specified
classical master junction, not a claimed smooth Einstein-Maxwell or
semiclassical four-dimensional junction. Its interpretation as the leading
physical overlap still requires control of the omitted background and sources.

For |sigma|<=1/2, the exact throat formulas give uniformly

    u(+/-R)=D_plus/minus R^2+O(||D||),
    R d_out u-2u=O(||D||).                             (14)

The changes of the remainders in (14) from sigma=0 are O(sigma^2 ||D||).
This follows by expanding the even analytic trigonometric quotients in
the exact solution; the terms linear in R cancel. It avoids division by
either end coefficient and covers zero or oppositely signed end data.

Since Z0(s)=s^2(1+s)/(1+3s), equations (13)-(14) imply

    A=lambda^-2[D_end+O((s+R^-2)||D||)],
    |M| J_s(L)=O(lambda^-2(s+R^-2)||D||).               (15)

For the second estimate, use
`s Z0'/Z0=2-2s/[(1+s)(1+3s)]`, the flux ratio
`(lambda^2+s^2)/f(s)=(1+s)^2(1+R^-2)`, and the bounded quantity

    p(s)J_s(L)/s <= [1/5+(3/2)s+3s^2]/(1+3s)^2.

Apply (12) to both mouths. The exact linear map from D to the normalized
outer master data takes the form

    lambda^2 (W_plus(L),W_minus(L)) = (Id+E_sigma) D,
    ||E_sigma||=O(s+R^-2+eta),
    ||E_sigma-E_0||=O(sigma^2/R^2+eta),                (16)

uniformly for small s, large R and eta<=eta0<1. The constants
do not grow with L: all outer-radius growth is included in eta.
For sufficiently small displayed error, (Id+E_sigma) is invertible.
Prescribing the static-normalized outer drive W_end(L)=(3/2)epsilon_end
therefore gives

    D=(3/2)lambda^2[Id+O(s+R^-2+eta)] epsilon.           (17)

The coefficient correction from finite frequency alone, at the same
outer drive, is relatively O(sigma^2/R^2+eta). The core throat response
itself still has the larger O(sigma^2) correction at fixed coefficients.
These are operator statements; their error norms multiply the norm of
the two-component drive, not a possibly vanishing individual component.

There is an additional conservative nonresonance check for this junction.
Writing (6) in optical time units gives potential V/lambda^2 on each
RN annulus and 6(1+rho^2) on the product core. The polynomial in V obeys

    9r^4-18r^2+16r-4-3r^4
                      =6z^4+24z^3+18z^2+4z >=0.

Consequently V >= (2/3)z^2/(1+z)^4. This last function has its only
interior extremum at z=1, a maximum. With zero outer Dirichlet data,
the self-adjoint master operator is bounded below by

    G^2=min{6, (2/(3lambda^2))
                    min[s^2/(1+s)^4,L^2/(1+L)^4]}.     (18)

Value and flux matching cancel the interface terms. Thus sigma^2<G^2
excludes a homogeneous resonance and bounds the normalizable resolvent.
This statement concerns this positive classical master operator only;
it is not a spectral gap of the full Dirac problem or a quantum adiabatic theorem.

## 4. The selected scale family and its unresolved input

Use the retained g^2=q^-2, lambda=1/(16pi q^3), delta=d/ell=(log q)^-2,
and choose R=lambda^-1/3. Then s=lambda^(2/3). Extend the local RN
annulus to the prior outer port radius r_out/P=q^-3/2 d/P, so
`L=16pi q^(3/2)/(log q)^2-1`. At the leading orbital frequency,

    sigma=(log q)^3/(2 sqrt(pi) q^(3/2)),
    eta=O(q^-5 (log q)^6)+O(q^-6 (log q)^2),
    s+R^-2=O(q^-2).                                   (19)

The intermediate terms of (11) are smaller than these leading terms.
The dynamic end-coefficient error in (16) vanishes as
`O(q^-5 (log q)^6)`; the core's relative harmonic correction is
`O(q^-3 (log q)^6)`. The conservative bound (18) eventually equals six.
These estimates close the finite-frequency classical scalar propagation
through the prescribed RN annuli; they do not compute the actual outer drive.

In particular, companion tides, radiative boundary data, relative
coordinate/monopole frames, charged-spinor phases, other dipole/vector
degrees, higher multipoles, and semiclassical metric/current response
remain outside (17). A low-frequency harmonic particular solution does
not specify preparation, free initial excitations or the actual quantum state.
The leading product geometry also requires a controlled comparison with
the quantum-open background across the overlap. No thin shell is asserted
to solve those omitted four-dimensional equations.

## 5. Verification and contribution decision

Run `checks/dynamic_rn_dipole.py --output evidence/dynamic-rn-dipole.json`.
The checker constructs the background and linear four-dimensional tensors,
obtains the physical master directly from Maxwell before substituting its
equation, and checks every reconstructed component. It verifies the static
solution/coefficient, residual gauge cancellation, Volterra primitives and
positive-potential bound. Independent explicit and implicit ODE solvers
compare 36 annulus cases with the inequalities and 12 two-mouth transfer
matrices with (16). Four scale evaluations supplement (19).

The [source audit](../evidence/dynamic-rn-dipole-source-review.json) records
the earlier dipole work and scoped current comparisons. Garcia-Saenz's
[dynamical-running paper](https://arxiv.org/html/2610.01152v1), Section 4.3
and Appendix B, treats axial RN modes and the special dipole normalization.
Its vanishing extremal logarithmic running does not set the full response
or determine a wormhole's boundary condition. Rai-Santoni's
[static Love analysis](https://arxiv.org/html/2404.06544v2), Section 5,
uses polar/axial duality and emphasizes source/response identification.
Neither partial comparison is a full proof review of those papers.

Continue physical long-return damping with matched finite-duration band
transfer as fallback. Retain this verified classical propagation result;
the actual binary, semiclassical junction and quantum response remain
decisive. The RN master and reconstruction are prior work. No consequential
new physical stability/transport result, actual shortcut or completed
program-wide novelty review is accepted here.

Remaining to-do list: derive and match the actual binary outer drive,
boundary frames/holonomy and remaining modes through the quantum-open
geometry; control the time-dependent quantum state, renormalized response,
exterior determinant and discarded tail; finish the consequential endpoint
and literature/proof/novelty review; independently verify, GCP-certify
suitable mathematics, and prepare the reproducible manuscript.
