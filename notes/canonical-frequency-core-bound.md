# Canonical frequency bounds at the centered seed

2026-10-20, E003/S168. Ordinary argument, followed by focused algebraic
verification. The full research goal remains incomplete.

The canonical Dirac shape Hessian at the retained centered seed admits
the following frequency-independent comparison on its compact smooth
Euclidean core:

    |Qcan[T]-Qgeo[T]| <= C[.001 J2+3e24 J1+3e26 J0],             (1)
    Qcan[T] >= C[J2/(10pi)-1e25 J1-1e27 J0].                   (2)

Qcan includes the massive angular tower, its prescribed local contacts
and the canonical LLL. It excludes the four-dimensional Weyl action and
classical action, whose coupled frequency polynomial was derived in
[the preceding note](local-coupled-frequency-polynomial.md). Qgeo is
the complete geodesic comparison from
[the curved seed calculation](stationary-curved-quantum-tt.md), section 2.
The norms are Jj from that note, with BOTH orbit derivatives included.
These are direct shape forms, not York TT forms or reduced physical energies.

In particular, if T is a nonzero bosonic thermal Fourier mode with
dimensionless frequency nu and a_* = sup a in
gbar=dY^2+a(Y)^2 d tau^2, then

    nu^2/a_*^2 >= 1e28  ==>  Qcan[T] >= C J2/(20pi).           (3)

This is a coarse sufficient threshold, not a predicted crossing. It is
uniform in the support's distance from the outer seam, on this fixed seed.
It does not overcome the negative Weyl logarithm by itself. A separate
result below extends the scalar model's logarithmic estimate to every
frequency, while leaving its comparison with the actual Dirac form open.

## 1. Core, state and the extension beyond static tangents

Keep the fixed centered seed, flux, clock, prepared HH period, paired
zero-source reflecting domain, full angular multiplicities and common
four-dimensional subtraction. Begin with real smooth spherical T(tau,Y),
periodic in tau and supported away from the outer seam and bolt. Finite
Fourier sums suffice. The even double extension is smooth at the seam
because T vanishes near it. No variation changes the endpoint projector,
spin structure or period. No exact semiclassical background is assumed.

In the fixed Hilbert representation the orbit Dirac variation is

    Dsharp(T)=exp(-T/2) Dsharp exp(-T/2).

It retains the paired adjoint/grading identities: multiplication by T
commutes with their constant fiber intertwiners. Let R be the reflection
lift used for the square and Gamma the two-dimensional chirality. At the
basepoint R Dsharp R=-Dsharp; even T preserves this relation. Gamma
anticommutes with Dsharp(T), and R anticommutes with Gamma. Consequently
Dsharp(T)^2 and its heat kernel commute with both R and Gamma; conjugation
by Gamma makes the reflected spin trace negative. That trace is zero
along the entire conformal path, including time-dependent paths. The
physical action remains half of the evenly varied double action. This
is the same image/phase allocation as the
[paired-domain audit](stationary-paired-dirac-phase.md) and
[boundary-layer identity](stationary-massive-boundary-layer.md), not a
replacement by an unpaired self-adjoint Dirac problem.

Sections 2-6 of the curved seed calculation estimate kernels and local
contacts on the full two-dimensional orbit. Their normal quadratic
probes q_v=(v.z)^2 range over every tangent direction. The geodesic flow,
angular Hessian identity, two-endpoint gluing and constant-contact
calculation do not set the time derivative of T to zero. In particular:

- the full two-dimensional separated error is bounded by 1e8 C/d^4;
- the ultraviolet second-contact allowance is 1e13 C;
- the row-sum-subtracted global gluing form costs at most 6e24 C J0;
- the actual constant contact has |Theta_p0|<=1e20 q^2(1+1/Y),
  with density p0=Theta_p0/C in gbar units;
- the curvature-gradient contact is the fixed equation (3) there.

These estimates therefore apply to the present T without a new
frequency-dependent constant. They are background kernel estimates;
we are not applying static source estimates to a time-dependent metric.
Only the comparison with the auxiliary radial scalar model in section 7,
the radial IMS allocation and the static TT return are NOT imported.

The LLL Polyakov cocycle, already fixed in
[the seed form calculation](stationary-seed-quantum-form.md), gives
Q_LLL=-q J1/(12pi) for these tangents as well. The mass vertex vanishes
for the LLL but this local term does not. All massive finite contacts
are retained in the same scheme. The literal logarithmic shape chart
is used; no off-shell equation removes its second contact.

At nonzero bosonic Matsubara frequencies the Euclidean second variation
equals the upper-half-plane retarded form in the fixed HH state, by the
[thermal spectral representation](thermal-shape-response-and-subtractions.md).
Equation (1) is thus an actual seed bound at those positive frequencies.
The static equilibrium form is not identified with the retarded zero
limit. Continuation to every non-Matsubara positive imaginary frequency
and to a nonlinear exterior is not supplied by this argument.

## 2. Explicit lower-order budget and the seam term

For transparency, the stronger intermediate comparison is

    |Qcan-Qgeo| <= .001 C J2+1e17 C J1+1e25 C J0
                     +2e25 C (J0+Jminus1),
    Jminus1=integral_gbar T^2/Y.                              (4)

Here is the numerical allocation. In the gain-two form estimate use
epsilon=1e-3, K=1e13, R=1/400 and b^2=epsilon/(4K). Then
2 log(R/b)<28, so its gradient coefficient is below 3e14 before the
unit-circle angular norm. The identities

    ||XT||_SM^2=pi J1,
    ||X^2T||_SM^2 <= (3pi/4)J2+(pi/2)J1

give fourth loss at most epsilon*pi/4<.001. Together with the finite
second contact, the explicit curvature-gradient contact, and the LLL
(q/C<400), the gradient reserve 1e17 suffices. The 6e24 gluing remainder
fits 1e25. For p0, (q/C)^2<160000 gives
1e20 q^2/C <1.6e25 C<2e25 C. The distributional seam contact pairs to
zero on this core; it is not erased for nonzero boundary traces.

Unlike the unbounded Weyl logarithm, this p0 weight is lower order. On
0<Y<1/10 the retained |h|<=Y^2, h=log(a/a(0)), implies
amax/amin<=exp(.02)<5/4. Since T(0,tau)=0,

    |T(Y,tau)|^2/Y <= integral_0^Y |partial_Y T(s,tau)|^2 ds.

Integrating this inequality over the collar and then tau gives at most
J1/8. On its complement 1/Y<=10. Therefore

    Jminus1 <= J1/8+10 J0.                                  (5)

Substitution in (4) yields (1):
1e17+2e25/8<3e24 and 1e25+22e25<3e26. This is an explicit
zero-trace estimate, not an assumption that the endpoint stress is bounded.

The benchmark has 13.9<kappa4<17. The lower bound is already proved;
the upper bound follows from H_q<=1+log q, log(2000)<8,
gamma_E>0 and |log(mu0 B)|<1e-21. Since |Rbar|<=2,

    Qgeo >= C kappa4 J2/(30pi)-C kappa4 J1/(45pi).

Using its positive tower remainder and (1) proves (2), with ample
reserve. No quartic finite coefficient has been chosen or changed.
The comparison (1) is two-sided, which matters if one later uses three
frequency anchors: above all three nodes the middle interpolation weight
is negative, so three lower bounds alone would not imply a lower bound.

## 3. The sufficient nonzero-frequency threshold

The seed's canonical warp has finite positive a_*; it vanishes only
at the bolt. On the nu Fourier sector the scalar orbit Laplacian is

    L_nu=-a^-1 partial_Y(a partial_Y)+nu^2/a^2 >= z,
    z=nu^2/a_*^2.

One may use its Friedrichs scalar realization merely to express these
norm inequalities; this is not a choice of gravitational boundary data.
For compact tangents, spectral calculus or Cauchy-Schwarz gives

    J1<=J2/z,  J0<=J2/z^2.

Thus (2) has coefficient C[1/(10pi)-1e25/z-1e27/z^2]. At z>=1e28
it exceeds C/(20pi), proving (3). Frequencies here are conjugate to tau
in gbar, so this threshold must not be inserted as an unconverted
laboratory or proper physical frequency. No Planck-scale validity is
inferred. There is no claim that every smaller frequency is negative.

For any fixed smooth compact radial profile and each nonzero thermal
frequency, (1) also gives an actual finite form enclosure in terms of
Qgeo and Jj. This is not a numerical evaluation of three exterior
operator anchors or a uniform norm bound on their Schur complement.

## 4. Frequency-uniform logarithmic estimate for the scalar model

Keep exactly the complete even auxiliary collar from section 7 of the
curved seed calculation. Put a_ext=a(0)exp(h_ext),

    H_0=-partial_Y^2+V,
    V=h_ext''/2+h_ext'^2/4, |V|<=22,
    H_nu=H_0+nu^2/a_ext^2 >=0.

H_0=(partial_Y-h_ext'/2)^*(partial_Y-h_ext'/2), and
H_nu+24>=-partial_Y^2 for every real nu. For a smooth even compact v
write w=H_nu v. Operator monotonicity of log, with positive
regularization followed by its removal, gives

    <w,log sqrt(H_nu)w>
      >= <w,log|D|w>-12<v,H_nu v>.

Indeed log(1+24/H_nu)<=24/H_nu on its spectral support. Multiplication
by w=H_nu v controls the lower endpoint without a spectral-gap assumption.
[Beckner, Corollary 2, equation (10)](https://arxiv.org/pdf/math/0612486v2)
with n=1,p=2 now gives

    <w,log sqrt(H_nu)w>
      >= integral -log|Y| |w|^2+D_H||w||^2-12<v,H_nu v>,     (6)
    D_H=log2+psi(1/4)=-gamma_E-2log2-pi/2 > -4.

The constant is independent of nu, despite the frequency potential.
The full-line proof applies to even v; divide by two for physical-half
norms and by the common thermal mode norm where appropriate. The
quarter-shift digamma constant and its normalization are established
literature, not a new uncertainty principle.

The accepted all-frequency multiplier bound applied to H_nu says

    <v,Kflat(sqrt(H_nu))v>
      >= C/(30pi)[<w,log sqrt(H_nu)w>
                        -(log(mu0 ell_o)+38/15)||w||^2].

Here Kflat(p) is the function of momentum p, so its spectral argument
is sqrt(H_nu). For support in |Y|<=.02, v=sqrt(a_ext)T and
H_nu v=-sqrt(a_ext) Delta_nu T. The seed inequality
omega<=-log|Y|-.5log(x/C)+.5log2 and
-log(mu0 ell_A)>9 therefore imply

    <v,Kflat(sqrt(H_nu))v>
       -C/(30pi) integral(omega+19/60)|H_nu v|^2
      >= C/(30pi)[(33/20)||H_nu v||^2-12<v,H_nu v>].         (7)

The rational margin follows from
9-38/15-19/60-4-(log2)/2>33/20, using log2<7/10.
This proves that the scalar-model logarithmic cancellation itself does
not fail at nonzero frequency. It does NOT extend the static TT estimate
of the rest of W4 to arbitrary time-dependent tangents.

## 5. What remains between these estimates and stability

Equations (1)-(3) concern the actual canonical seed form. Equations
(6)-(7) concern the auxiliary scalar model. They must not be conflated.
The earlier scalar comparison used a bounded radial potential V. In a
nonzero sector, subtracting only D_Y^2 adds nu^2/a_ext^2; its supremum
and derivative grow with nu. Substituting this in the old Feynman-Kac
estimate would replace its fixed gain-two constants by frequency-growing
ones. That does not establish the uniform comparison needed on escaping
annuli. Comparing against D_Y^2+nu^2/a(0)^2 still leaves a variable
potential of size O(nu^2 Y^2), requiring a new joint space/frequency
estimate. It has not been supplied here.

For example, on a compact annulus [d,2d], omega grows as log(1/d).
The lower bound (2) supplies a fixed positive J2 coefficient; subtracting
the known negative Weyl term therefore leaves a lower estimate tending
to minus infinity relative to J2 as d->0. A failed lower estimate is not
an instability or a negative value of the actual form. The scalar model
indicates the missing logarithmic compensation, but cannot replace its
proof for Qcan. The next specific task is that two-dimensional
Dirac/geodesic/scalar comparison with uniform contacts and returns.

Even after that comparison, the full local polynomial, frequency-
dependent conformal/shape coupling, initial/incoming Ward constraints,
Maxwell/current response, admitted boundary realization and transfer
away from the seed must be treated. Positivity of the direct canonical
shape form is neither stability nor exclusion of the fast scalar branch.
S168.3-S168.4 and all S167 tasks remain open. No action coefficient,
state, physical disturbance class or requirement is changed.

## Verification scope

The companion checker pins the predecessor evidence and scientific
sources, checks the rational budgets, the arbitrary-jet half-density
identity and the mode norm reductions, and tests the operator-log step
on independent noncommuting finite matrices. Those finite fixtures are
not the wormhole spectrum. Continuum estimates above are ordinary
arguments reusing the stated kernel bounds, not formal certificates.
No Lean invocation or authoritative numerical spectrum is selected.

Remaining to-do list: actual uniform logarithmic frequency comparison;
coupled retarded constraints and stability transfer or justified exclusion;
signed coefficient, cold/higher response, nonlinear solve and compatibility;
exterior/state, charge/return, opening, stability, transport, preparation,
separate creation and full audit.
