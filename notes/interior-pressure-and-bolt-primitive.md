# Interior canonical pressure and the bolt lapse primitive

2026-10-19, E003/S148 supporting S139/S074. Ordinary argument and
verification protocol before execution.

Bolt regularity determines the canonical pressure integration
constant. This note bounds that pressure functional at fixed
interior cuts, derives the transport to a cut of fixed outer proper
distance, and gives an exact signed bolt lapse primitive. The
moving cut's second response still needs one spatial trace bound.
The lapse identity is not a proof of full physical bolt smoothness.

## 1. Fixed disk, source normalization and hypotheses

Use the common disk of the centered reconstruction,

    ghat_2=A^2 exp(2sigma)(dchi^2+sinh^2 chi dvarphi^2),
    sigma=sigma0+b, 0<=chi<D0<82,
    W_b=exp(2sigma)sinh^2 chi.

The physical period, paired domain, B,A,C and mu0 stay fixed.
All variations are stationary and smooth radial fields at the bolt,
or the stated C5 radial limits. The full global shape perturbation
is conformal in this common disk, as required by S145. Its norm
is the supremum of S145's chart/cap/global norms over the interior
cover; b has norm at most 1e-3, and direction norms are U,V.
In particular |b|<=.001, |b_chi|<=.008 and
|h|<=U, |h_chi|<=8U, since exp(sigma0)<7.31<8.

Let c_* be the seed point with outer proper distance Y0=1/25.
The normal charts of S145 stay before the actual outer boundary
on 0<=chi<=c_*. This is a sufficient coefficient hypothesis,
not a demonstration that every previously constructed endpoint
lies in the raw C5 ball. The seed outer estimate gives r/B<80
there: Y0<=sqrt(8)log[1/(1-B/r)] implies
r/B<=1/[1-exp(-.04/sqrt(8))]<80.

Define normalized canonical scalar trace and pressure by

    S=(2G/(B^2 C))Theta,
    P=(2G/(B^2 C))(B^2 p_hat).

S is the Lorentzian trace. Its sign is opposite to S145's
canonical scalar ACTION source, whose mass term is -M tr R.
Their absolute derivative bounds agree; no coordinate-density
factor is inserted into this scalar Ward equation.

The seed satisfies |Theta0|<=4e15 q^2(r/B) and
Gq^2/(CB^2)<6.8e-28, so |S0|<4.352e-10 on this cover.
S145 also supplies a first-response bound, not just its second:
the nonconstant part on its parameter circle is bounded by

    K=2 eta_G[(2/3)1e17+7e11*24*256^3+1e17]<5e-12,
    eta_G<8.5e-33.

The terms are the subtracted local remainder, actual-domain
comparison and common curvature restoration. The fixed vacuum
constant has zero scalar variation, not zero value.
Cauchy on radius 1/(200U) gives |D_h S|<1e-9 U.
Together with the real background radius and S145,

    ||S||inf<5e-10,
    ||D_h S||inf<1e-9 U,
    ||D_h D_k S||inf<2e-7 U V.                      (1)

The complete angular remainder converges uniformly in these
norms on the interior cover. The individual renormalized mode
diagonals and local terms are continuous in the observation;
the limit and its two parameter derivatives are therefore
continuous at the bolt. This statement gives no quantitative
spatial derivative or modulus of continuity.

## 2. HH regularity fixes the pressure functional

The canonical Ward identity is

    (W_b P)'=W_b' S/2.

A nonzero homogeneous integration constant gives P~constant/chi^2
at the smooth Euclidean bolt. Regularity excludes it. Hence

    P_b(c)=F_b[S](c)
          =[2W_b(c)]^-1 integral_0^c W_b'(t)S(t)dt,
    P_b(0)=S(0)/2.                                  (2)

This is a prescription for the canonical HH pressure on every
specified geometry, not an independent state parameter. It does
not compute the unknown geometry or its trace. It includes the
full canonical massive trace, not the outer regular remainder
alone. S147's separated cubic allocation is made only after
the full pressure and its matching value have been reconstructed.

The seed W0 has one maximum, past/current ratio at most 1.01,
and total variation on [0,c] at most 1.02 W0(c), as proved in
the centered complete-response calculation. Consequently

    integral_0^c W_b(t)/W_b(c)dt <84,
    J_b(c)=integral_0^c |W_b'(t)|dt/W_b(c) <2.4.      (3)

For example the upper bounds are
exp(.004)*1.01*82 and
exp(.004)[1.02+2*.008*1.01*82].
No monotonicity of the perturbed warp is assumed.

With K_b(t,c)=W_b(t)/W_b(c), ell_b=W_b'/W_b and
d_h=h(t)-h(c),

    F_b[S]=(1/2)integral_0^c K_b ell_b S,
    D_h F_b[S]=integral_0^c K_b[d_h ell_b+h']S,
    D_h D_k F_b[S]
      =integral_0^c K_b[2d_h d_k ell_b
                         +2d_h k'+2d_k h']S.        (4)

Here the derivatives act on the operator, holding S fixed.
Since |d_h|<=2U and |h'|<=8U, (3) gives

    ||F_b||<=1.2,
    ||D_h F_b||<=700 U,
    ||D_h D_k F_b||<=5500 U V.                      (5)

The unsimplified last two constants are
2*2.4+8*84=676.8 and 8*2.4+64*84=5395.2.
Every operator derivative in (4) annihilates a spatial constant:
F_b[constant]=constant/2 for every geometry.

Applying the full product rule to (1), not freezing the Ward kernel,

    ||P_b||inf<6e-10,
    ||D_h P_b||inf<3.6e-7 U,
    ||D_h D_k P_b||inf<4.4e-6 U V.                  (6)

The mixed subtotal is 4.39e-6. The bolt values and parameter
derivatives follow by continuous extension of (2); no new
integration constant appears upon differentiating.

## 3. Transport to the outer proper cut

The cut in S147 is Y=1/20, not a fixed chi. In this common
conformal disk its location c_b is determined by

    Y_b(c)=integral_c^D0 E_b(t)dt=1/20,
    E_b=exp(sigma0+b).

The whole boundary collar participates in this integral.
The global bounds .39<exp(sigma0)<7.31 and |b|<=.001 imply
that the seed Y0 at c_b lies between .05 exp(-.001) and
.05 exp(.001). Thus c_b lies inside the cover ending at c_*.
The value estimate in (6) applies to P_b(c_b).

Let I_h=integral_c^D0 E_b h and similarly I_hk. Exact
implicit differentiation gives

    c_h=I_h/E_b(c),
    c_hk=I_hk/E_b(c)-h(c)c_k-k(c)c_h
                          -(log E_b)'(c)c_h c_k.    (7)

The outer seed curvature and geodesic boundary give
|a0,Y/a0|<=2Y0. At these cuts c>1 and
sigma0'=-E0(a0,Y/a0)-coth c, so
|(log E_b)'|<3. Therefore

    |c_h|<.13 U, |c_hk|<.45 U V.                   (8)

These are derived cut responses, not independent free quantities.
The fixed-proper-cut pressure is Pm(b)=P_b(c_b), with

    D_h Pm=P_h+P_chi c_h,
    D_h D_k Pm=P_hk+(P_h)_chi c_k+(P_k)_chi c_h
                       +P_chichi c_h c_k+P_chi c_hk. (9)

All right-hand quantities are evaluated at c_b. In particular
the fixed-cut bound (6) cannot be relabeled as the complete
moving-cut Hessian.

Ward differentiation yields, with ell=W'/W,

    P_chi=ell(S/2-P),
    (P_h)_chi=2h'(S/2-P)+ell(S_h/2-P_h),
    P_chichi=ell'(S/2-P)+ell(S_chi/2-P_chi).         (10)

On this cut |ell|<2 and |ell'|<120. For the seed these follow
from ell0=-2E0 a0,Y/a0, |E0|<7.31,
|sigma0'|<3 and |(a0,Y/a0)_Y|<=1+(2Y0)^2.
The b terms add 2b' and 2b''; at the observation
|b''_chi|<=100*1e-3 from its covariant C2 bound and connection.
These allowances give

    |P_chi|<2e-9,
    |(P_h)_chi|<7.4e-7 U,
    |P_chichi|<1.1e-7+|S_chi|.                     (11)

Consequently the value and first pressure responses at the actual
proper cut are bounded:

    |Pm|<6e-10, |D_h Pm|<3.7e-7 U.

IF a spatial bound |S_chi|<=J is supplied there, then

    |D_h D_k Pm|<(4.6e-6+.017 J) U V.              (12)

The inherited C5 value theorem does not supply J. Neither J nor
the last transport term is set to zero. If the trace is not
differentiable at the cut, (9)-(12) require a stronger regularity
space or a joint transport argument; this note does not claim
their unconditional C2 existence.

## 4. An exact signed bolt lapse primitive

There is a useful value statement even before a spatial trace
derivative is controlled. In the common disk put r=r(chi) and
use the exact geometric identity

    F=C exp(-2sigma)(r_chi/r)^2.

For the canonical sector, the lapse-constraint source is

    alpha_src,chi
      =(B^2 C/2) r_chi (2P-S)/(r^3 F).

This defines the source integral on a specified geometry; it is
not an assertion that that geometry solves the constraint.
On a bolt neighborhood where W'>0, Ward gives

    alpha_src,chi=-H_bolt P_chi,
    H_bolt=B^2 C r_chi W/(r^3 F W')
          =B^2 exp(2sigma) W/(r r_chi W').          (13)

For a regular nonextremal bolt,
r=r_h+r2 chi^2+O(chi^4), r2>0, and W=w2 chi^2+O(chi^4).
The coefficient extends smoothly and evenly:

    H_bolt(0)=B^2 exp(2sigma_h)/(4r_h r2)
             =B^2 C/[F'_h r_h^3].

This is finite at fixed positive temperature. It is generally
of order 1/h for the retained cold gap, not bounded by one.
An arbitrary value of this coefficient is not declared small.

Writing Z=P-P(0), integration by parts gives the exact primitive

    alpha_src(c)-alpha_src(0)
       =-H_bolt(c) Z(c)+integral_0^c H_bolt'(t)Z(t)dt. (14)

Continuity of S implies Z(0)=0; a C1 extension of H_bolt suffices
to take the lower endpoint limit. Thus the signed lapse primitive
is finite under the actual value/continuity hypotheses. It does
not require an integral of an absolute unsigned bound 1/chi.

For the norm B(H)=|H(c)|+integral_0^c |H'|, (6) gives

    ||Z||inf<8.5e-10,
    ||Z_h||inf<3.7e-7 U,
    ||Z_hk||inf<4.5e-6 U V.

If the geometry family has B(H)<=B0,
B(H_h)<=B1 U and B(H_hk)<=B2 U V, the mixed primitive obeys

    |D_h D_k alpha_src(c)|
      <[4.5e-6 B0+7.4e-7 B1+8.5e-10 B2] U V.      (15)

The BV coefficient bounds in (15) remain geometric inputs. They
are not replaced by one, and the inverse-gap size of B0 means
this estimate alone does not give a useful physical error budget.
Weighted cold response and/or the coupled horizon equations
are still needed to control that size.

## 5. What the primitive does and does not establish

Equation (14) is a signed primitive, not an absolute lapse-
variation bound or a smoothness theorem. This distinction is
necessary even for W=chi^2. Let t=log(1/chi)>=1 and

    S(chi)=sin(t^2)/t, S(0)=0.

It is continuous, bounded and smooth for chi>0. Its regular
pressure is P(t)=integral_0^infinity exp(-2s)S(t+s)ds.
One integration by parts in the phase proves |P(t)|<=1/t^2.
Yet integral_1^infinity |S(t)|dt diverges:
it equals (1/2)integral_1^infinity |sin v|dv/v.
Hence integral |2P-S|dchi/chi diverges, whereas the signed
primitive with H_bolt=1 equals -P(c)+P(0) and is finite.
This is a counterexample to inference from mere scalar continuity,
not a proposed Dirac trace or an admissible physical solution.

A sufficient condition for absolute control is a Dini bound on
the trace oscillation. With
omega_S(c)=sup_(0<=t<=c)|S(t)-S(c)|, the exact constant
cancellation in (2) gives

    |2P(c)-S(c)|<=J_b(c) omega_S(c)<=2.4 omega_S(c).

The parameter responses have the analogous bound
2.4 omega_hk+1400(U omega_k+V omega_h)+11000UV omega_S.
An integrable omega(c)/c is sufficient near a fixed nonextremal
bolt, and a power modulus is one way to supply it. This note
does not infer such a quantitative modulus from the pointwise
source estimate.

## 6. Verification and research consequence

The checker is to validate the exact Ward derivatives, moving-cut
chain rule and lapse primitive; all numerical constants and
normalizations; and the cold coefficient's inverse-gap behavior.
Independent quadratures and mixed finite differences use smooth
manufactured radial traces and nonzero regular pressures. Cut
locations are also solved independently before differencing.
Omitting a kernel or cut-transport term must be detected.
These manufactured traces are not physical quantum states.

The pre-execution protocol fixes 55-digit arithmetic for the
manufactured checks, 32 and 64 Gauss nodes, mixed-difference steps
.001 and .0005, relative-to-(1+|target|) error tolerance 2e-7,
and absolute quadrature/primitive tolerance 1e-11. Step refinement
must reduce the error by a factor below .35 unless both errors
are below 1e-35. The independent proper-cut solve uses adaptive
quadrature. Failed runs, if any, retain their source snapshots.

The retained model and ensemble originate in the
[Bintanja--Freivogel--Rolph construction](https://arxiv.org/html/2308.00871v2).
Its distinction between a semiclassical construction, phase choice
and preparation remains in force. The estimates here derive from
the repository's fixed-domain source calculations and the explicit
Ward identities, not from assuming the original paper already
proves this finite-background match.

The canonical pressure integration constant is no longer an
unassigned datum on a specified geometry: (2) fixes it. Its
fixed-cut response and the proper-cut value/first response are
bounded. A complete proper-cut Hessian still needs J, and the
bolt primitive still needs quantitative cold geometric/source
control for the physical error norm. The local/radius sectors,
global inverse composition and physical gates remain open.

Remaining to-do list: bound the spatial trace/cut transport and
cold bolt response in the actual solution space; complete the
other source mappings and global exterior/state solve; determine
charge and signed source; establish opening, stability, transport,
preparation and creation.
