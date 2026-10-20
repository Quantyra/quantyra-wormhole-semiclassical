# A strong outer bound for the actual massive-to-geodesic remainder

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary proof and
verification contract before execution. **Actual quantum remainder
controlled; the curved geodesic-to-flat comparison remains open.**

The first derivative contact can be estimated by a covariant probe,
without differentiating a bound for the second contact. This gives
a strong weighted estimate for the entire actual canonical massive
Hessian minus the retained geodesic model, including its constant
contact, global term, thermal contributions and derivative contacts.

## 1. Statement and fixed allocation

Keep the seed, paired domain, common subtraction and physical
normalizations of [the curved proof](stationary-curved-quantum-tt.md)
and [the preceding global/contact estimate](outer-global-and-zero-contact-bound.md).
Put xi=lambda B^2, u=B/r, delta=7/2, x=log Y, and

    E=H_mass-Q_geo,
    D_epsilon v=1_(0<Y<epsilon) C u^4/(16pi xi^2) E v.

H_mass is the fixed-canonical-area Hessian. Q_geo is the complete
positive geodesic comparison, with its common fourth restoration
and full-flow returns. For stationary smooth physical directions
v supported in 0<Y<epsilon and 0<epsilon<=1e-15,

    ||Y^-7/2 D_epsilon v||_L2(dx)
        <=1e64 epsilon ||Y^-7/2 v||_H4log(dx).           (1)

Use even extension across the seam and zero extension of the
weighted input on the log line. The estimate extends by completion.
For b=2T, D_epsilon is the dimensionless physical shape contribution
of E, with dual output f_b=(f_T-f_S/3)/2.
The actual density is retained.

Consequently the accepted complete normal inverse gives

    ||N^-1 lambda^2(0,D_epsilon)|| <8.7e83 epsilon,
    <1e-30 when epsilon=1e-120.                        (2)

This is a component of the full outer comparison. Q_geo minus the
flat model, its combination with the Weyl/principal density, and
annulus matching still require estimates. Equation (1) includes the
global/constant block bounded in the previous note; that block must
not be added to (1) a second time.

## 2. A global probe supplied by covariance

On the doubled canonical orbit write

    gbar=dY^2+a(Y)^2 dt^2, a0=a(0),
    j=a/a0, v_*=a'/a0, V=j partial_Y.

Here v_* is a probe, not an unknown field or a new physical datum.
The even double has two regular bolts. V vanishes smoothly there.
At the seam a is C2,1, so V is C2,1 and v_* is C1,1.
The elementary identities are

    L_V gbar=2v_* gbar,
    v_*'=-(Rbar/2)j,
    v_*''=-(Rbar'/2)j-(Rbar/2)v_*,
    h' v_*'=-(Rbar/2)v_*, h=log j.                    (3)

Let T_mass be the Euclidean scalar action source per canonical
area. Diffeomorphism covariance, followed by the fixed-area density
conversion, gives at every Y>0

    H_mass v_*=j T_mass'+2v_* T_mass.                  (4)

The Lorentzian trace Theta used in the source notes has the opposite
Euclidean action sign: T_mass=-Theta/C. Only the absolute bounds
|T_mass|=|Theta|/C and |T_mass'|=|Theta'|/C are used below.

For clarity, (4) is an auxiliary identity on the FULL doubled
Dirac problem, whose restriction to even conformal directions
gives the original physical Hessian. v_* is odd; it is not
introduced as an allowed odd perturbation of the physical pair.
The local coefficients of the full-double kernel are the ones
acting on even physical directions after folding.

The smooth covariance input is the unitary naturality of the
Dirac operator under spin diffeomorphisms; see
[Dabrowski--Dossena](https://arxiv.org/pdf/1209.2021), inspected
2026-10-20. That source concerns the smooth case. Here is the
additional passage needed at this seam. Smooth a only inside
|Y|<2eta, preserving positivity, evenness, uniform C2,1 bounds and
the metric outside that interval. Keep the bolts unchanged.
For each smoothed metric use its smooth vector (a_eta/a0)partial_Y.
The diffeomorphism identity holds first for finite angular and
proper-time cutoffs. The isotopy preserves the spin structure,
time period and bolt domain.

Fix an observation Y>0, then take eta<Y/8. The local metric and
all its jets at that observation agree with the original metric.
Resolvent identities localize the smoothing difference to distance
at least Y/2. Endpoint derivatives cost only finite powers of M;
the intervening massive kernels have exponential decay in M Y.
The counting bound N(M)<=4CM^2 makes every such fixed-Y sum
convergent. Thus the scalar source and its first observation
derivative converge locally, as do their responses to v_*, while
the common local heat contacts at the observation are unchanged.
One may also approximate v_* by smooth fields agreeing near
the observation: its separated kernel contribution converges
by dominated convergence. This proves (4) as a pointwise identity
away from the seam. It does not assume a bounded sixth seam jet,
move the physical boundary, or drop a seam-supported term.

The accepted source bounds in
[the conormal source proof](stationary-centered-conormal-source.md)
give, with L(Y)=1+log(1/Y),

    |Theta-c C^2/(5pi Y)|<=1e20 q^2 L(Y),
    |Theta'-(-c C^2/(5pi Y^2))|<=2e22 q^2 L(Y)/Y.

The seed coefficient is |c|=|a_seed|sqrt(xi)/C^(3/2)<1.
Use xi/C<1, a_seed/C<.334, q/C>300 and L(Y)<=1/Y.
For Y<1/10 this yields the generous bounds

    |Theta|<2e20 q^2/Y, |Theta'|<3e22 q^2/Y^2.          (5)

The global polynomial guards calB<=xi+1 and |calB'|<=5C give
j<2 and |v_*|<7. The curvature guards |Rbar|<=2,
|Rbar'|<=17 imply from (3)

    |v_*'|<=2, |v_*''|<=24, |h'v_*'|<=7,
    |Hess v_*|<100 a.e. on the complete double.        (6)

These include the regular bolts. Equations (4)-(6), q/C<400,
give

    |H_mass v_*(Y)|<1e29 C/Y^2.                        (7)

## 3. Local high jets and the nonzero probe derivative

Use the exact polynomial

    P(u)=calB(u)=xi+u^2-2a_seed u^3
                         +Q0 u^4+d u^4(1-u)^2,
    j=sqrt(P/xi), u'=sqrt(P/C),
    v_*=P'(u)/(2sqrt(C xi)).

The retained guards are xi>=1000, a_seed/xi<2.01,
Q0/xi<3.01, |d|/xi<1e-15, 5.9<C/xi<6.1 and C<2^36.
On the complex disk |u|<=1/20, direct coefficient sums give

    |P/xi-1|<.001,
    |(P/xi)^(-1/2)-1|<.002.

The proper-coordinate integral Y(u) differs from sqrt(C/xi)u
by less than .002sqrt(C/xi)|u|. Its normalized derivative has
positive real part on this convex disk, so it is univalent.
On its boundary |Y(u)|>.998sqrt(5.9)/20>.1. Rouche's theorem
then supplies the analytic inverse for |Y|<.1. On this disk,
|j|<2, |v_*|<1 and |h'|<2, by the same polynomial bounds.
This is the analytic one-sided continuation, not the even double.

At a center 0<Y<1e-15 Cauchy on radius 1/20 gives Taylor
coefficient bounds 20^k for v_* and 2*20^k for h'.
Along a unit geodesic,

    Y''=h'(Y)(1-Y'^2), |Y'(0)|<=1.

The positive formal majorant

    z(0)=0, z'(0)=1,
    z''=2(1-20z)^(-1)[1+(z')^2]                      (8)

therefore bounds its jets through six. The finite Bell sums
sum_(k=1)^n k!20^k Bell_(n,k)(z',...,z^(n-k+1))
are below 1e30 for n<=6. Hence |X^n v_*|<1e30 on flow
segments of length less than Y/2, n<=6. Such segments stay
on the positive analytic side, since proper Y is 1-Lipschitz.

Finally, |P'''|<=17C and u<=Y/sqrt(5)<Y imply

    |P''(u)-2|<=17 C Y<1, 0<Y<=1e-15.

Since j>1/2 there,

    v_*'=P''j/(2C)>1/(4C).                            (9)

The small derivative is divided by its actual explicit lower
bound; it is not replaced by one.

## 4. Bounding the geodesic probe and the fixed Taylor remainder

Q_geo can act pointwise on v_* at Y>0 although v_* is not in
its global fourth/log domain. Define this action by smoothing at
the seam and taking the separated-kernel limit just described.
Let RY=Y/2. For r<RY the even Taylor remainder in g(rP_flow)
is at most r^6 sup|X^6 v_*|/720. The common bounds

    sum d_n |B_Mn(r)|<=72C/r^6,
    sum d_n |B_Mn(r)|<=1152C r^-6 exp(-9r)

give a near action below 1e30 C Y^2 after the full 2pi angular
integral. For r>=RY, the global Hessian bound (6) gives

    |v_*-cos(rP_flow)v_*|<=50r^2,
    |(r^2/2)X^2v_*|<=50r^2.

Their radial integral is below 1e5 C/Y^2. The quartic polynomial
term costs at most
400 C |X^4v_*(Y)|[1+log(2/Y)], using the exponential envelope
and integral_RY^infinity exp(-9r)dr/r<=1+log(1/RY).
The ordinary fourth restoration costs less than
10 C |X^4v_*|, using kappa4<100. Since
1+log(2/Y)<3/Y^2 for Y<1/10, all terms fit

    |Q_geo v_*(Y)|<1e35 C/Y^2.                        (10)

This includes long-flow escape and return, not just a local cap.

Use exactly the fixed R=1/400 geodesic Taylor allocation of
[the strong interior proof](strong-interior-quantum-remainder.md):

    E v=p0 v+A1.grad v+A2:Hess v+T3 v+B0 v,            (11)
    T3 v(X)=integral_(d<R) K(X,Z)
       [v(Z)-v(X)-r grad_e v(X)-r^2 Hess_ee v(X)/2]dZ.

The double's gain-two kernel obeys |K|<=1e8 C/d^4;
the remaining bounded kernel has both rows below 4e24 C,
so ||B0||_2<=8e24 C and B0 1=0. Also p0=H_mass 1,
with its already proved bound. These facts hold on the seam
charts as well as in the interior; only the old FIRST-contact
probe estimate required the interior restriction.

The Taylor operator satisfies the same Liouville-flow estimate

    ||T3 v||_0<=2e6 C ||nabla^3 v||_0.                (12)

For v_* at observation Y, split its point action at r=Y/2.
Below that scale use the local third jet <1e30. Above it,
the global Hessian bound gives the Taylor remainder at most
100r^2. With polar Jacobian at most 2r, the two contributions
are bounded by 2e38 C Y and
1e12 C log(2R/Y), respectively. They imply

    |T3 v_*(Y)|<1e40 C/Y^2.                          (13)

## 5. Fixing the derivative contacts

First fix A2 by the centered normal probes
P_e=chi(z)(e.z)^2, |e|=1, of the curved proof's section 6.
They have zero value and gradient, Hess P_e(0)=2 e tensor e,
sup|P_e|<=1/1024 and coordinate jets
|D^j P_e|<=3 j!100^j through five.
That proof bounds their complete local finite second response
by 1e13 C, with the known curvature contact separately below
1e6 C. Include G0 using its row bound, rather than omitting it.
Thus |E P_e|<1e22 C.

For the allocation (11), the uniform normal-chart bounds
|Gamma|<=15, |D Gamma|<=6000 and unit coordinate speed <2
give |X^3 P_e|<1e12 by the explicit third chain rule.
Hence |T3 P_e|<2e18 C and
|B0 P_e|<=8e24 C/1024. Evaluating (11) and polarizing bounds
the Frobenius norm by the generous reserve

    |A2|<1e30 C.                                     (14)

This uses the existing second-probe calculation, not a
differentiated coefficient envelope or a form-to-operator inference.

Apply (11) to v_*. Its gradient is purely radial, so at Y>0

    (A1)_Y v_*'=H_mass v_*-Q_geo v_*-p0 v_*
                            -A2:Hess v_*-T3 v_*-B0 v_*.

Equations (6), (7), (10), (13), (14), the bound on p0 and
|B0 v_*|<=8e24 C sup|v_*| put the numerator below
2e40 C/Y^2. Divide using (9), then use C<2^36:

    |(A1)_Y|<8e40 C^2/Y^2 <1e60 C/Y^2.               (15)

The tangential component does not act on stationary inputs.
This is the missing strong first-contact estimate. It has been
extracted from the actual point action on a covariant probe.

## 6. The strong weighted estimate

The physical density in (1) is rho Y^4/(16pi C), with rho<3/2.
Use ||Y^-7/2 Y^4 F||_L2(dY/Y)=||F||_L2(dY),
and 1/2<exp(h)<2 on the collar. If w=Y^-7/2 v, then

    partial_Y^k v=Y^(7/2-k)
                    product_(j=0)^(k-1)(partial_x+7/2-j)w.

The coefficient sums for k<=3 are below 40. In particular

    ||Y^-2 v'||_L2(dY)<=5 epsilon ||w||_H1,
    ||v''||_L2(dY)<=40 epsilon^2 ||w||_H2.

For a stationary radial field,

    |nabla^3 v|<=|v'''|+3|h'||v''|
                                    +( |h''|+2h'^2 )|v'|.

The Riccati bounds |h'|<=2Y, |h''|<=2 therefore give
||nabla^3 v||_0<=1000 epsilon ||w||_H3.
The H3 and lower norms are at most twice H4log.
Using (12), (14), (15), and the same p0/B0 argument as in
the preceding global-contact proof, safe component reserves are

| Component after physical density and weight | H4log-to-L2 bound |
| --- | --- |
| First derivative contact | 1e61 epsilon |
| Second derivative contacts | 1e34 epsilon^2 |
| Taylor remainder | 1e11 epsilon |
| Constant contact and bounded B0 | 1e27 epsilon^3 |

Their sum is below 1e64 epsilon, proving (1). The estimates apply
to all scales down to zero, not a finite mesh of end annuli.
Equation (2) follows from ||N^-1||<8.7e19/lambda^2.

## 7. Verification and remaining comparison

Freeze this argument and checker before execution. Pin the curved,
strong-interior, global-contact, flat-correction and source receipts.
Verify the covariant/density identity, polynomial geometry and
complex-disk reserves, geodesic majorant and Bell bounds, second
probe extraction, local/far point integrals, division by the small
probe derivative, weighted derivative identities and full component
budget. Retain every run and any failed source snapshot.
Finite algebraic checks support, but do not formally certify,
the ordinary covariance, limiting and operator arguments.

This controls H_mass-Q_geo. It leaves Q_geo minus the even flat
model, the actual local Weyl/principal-density combination and
the compact annulus. It does not supply a global inverse constant,
a nonlinear solution, signed compatibility or stability.

Remaining to-do list: curved geodesic-to-flat comparison and local
Weyl/density assembly; annulus matching; higher York maps and full
inverse/residual composition; nonlinear signed compatibility;
exterior/state, charge, return, opening, stability, transport,
preparation and creation.
