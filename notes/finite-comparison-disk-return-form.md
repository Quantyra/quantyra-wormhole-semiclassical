# Finite-comparison disk return, logarithmic form mapping and LLL bound

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary argument and
verification contract before execution. **The actual disk return has
quantitative value and first-derivative bounds, and the retained LLL
disk Hessian has a small lower bound at every finite comparison. The
remaining full quantum TT form and coupled inverse are still open.**

The [York gap](finite-comparison-york-gap.md) supplies the actual global
decomposition. Here the return to the disk slice is controlled as well.
The form-level logarithmic map uses fewer derivatives than the full
H4-log operator-domain isomorphism. This distinction permits TT form
work now, without claiming that the latter isomorphism is proved.

## 1. Geometry, components and the outer-fixed return

Keep the original finite path and physical measure of the
[scalar reference proof](finite-comparison-scalar-resolvent.md):

    g_ell=exp(2w)(d rho^2+a0^2 dtheta^2)+r0^2 exp(2S)dOmega^2,
    w=ell u+S, S=log[1+((t-1)/t)(exp(ell V)-1)],
    d rho_ell=exp(w)d rho, y=sqrt(lambda)rho,
    dmu_ell=exp(2w+2S)dmu0,
    |w|,|S|<epsilon=3e-8, |w_y|,|S_y|<.001.        (1)

Use (rho_ell,theta,sphere,sphere) component order throughout this
note. A genuine radial TT tensor is

    h=diag(b,-b-2c,c,c),
    b'+2(A+R)b=-2(A-R)c,
    A=a_ell'/a_ell, R=r_ell'/r_ell,               (2)

where primes here and below mean physical rho_ell derivatives unless
a y subscript is present. In particular b is the RADIAL component.
The point norm has the useful exact expression

    |h|^2=2(b+c)^2+4c^2.

Thus |b+c|<=|h|/sqrt(2), |c|<=|h|/2,
|b|<=sqrt(3)|h|/2 and |b+3c|<=sqrt(3/2)|h|.
The same c estimate applies to radial component derivatives.

The disk return is fixed at the outer end:

    eta=-a_ell integral_(rho_ell)^infinity (b+c)/a_ell d s,
    eta'-A eta=b+c,
    S_c=c/2-R eta,
    T_c=-(b+3c)/2-(A-R)eta.                       (3)

Subtracting L_eta g_ell makes the two orbit components equal. The
remaining scalar/radius disk components are exactly (S_c,T_c).
There is no free multiple of a_ell: it grows at infinity and is
excluded from physical L2 and the chosen outer-fixed domain.
For smooth spherical tensors at the bolt, b+c=O(rho_ell^2), so (3)
has eta=O(rho_ell) with regular odd vector parity. The subtraction
preserves the bolt period, rather than moving a bolt by a singular
areal-radius vector.

## 2. Quantitative return inverse and values

Set beta_u=exp(w+S) and v=beta_u eta. This unitary identification
uses proper-frame vector components and actual physical measure.
Let A0=partial_y log a0, R0=partial_y log r0 and B0=partial_y-A0.
Direct differentiation gives the exact dimensionless return map

    B_ell=lambda^-1/2 beta_u (partial_(rho_ell)-A) beta_u^-1
          =exp(-w)[B0-(2w_y+S_y)].                (4)

The seed guards a0,yy/a0>=1/2, A0>=0 and R0>=0 imply

    ||B0 v||^2
      =||v_y||^2+integral [a0,yy/a0+A0^2+2A0 R0]|v|^2 dmu0
      >=||v_y||^2+||A0 v||^2+||v||^2/2.          (5)

The integration is first on the regular vector core and then its
graph closure. B0 is one-to-one and onto physical L2: its regular
homogeneous vector a0 is not L2 at infinity; the formal adjoint
homogeneous vector 1/(a0^2 r0^2) is not L2 at the bolt. The gap (5)
makes its range closed, so the absent adjoint kernel makes it onto.
Its inverse is the outer-fixed Volterra integral in (3) at the seed.

With d=1/200=.005, (1),(4),(5) give

    ||(B_ell-B0)v||
      <=[exp(epsilon)-1+exp(epsilon) .003 sqrt(2)]||B0 v||
      <d||B0 v||.                                (6)

Indeed B_ell B0^-1 is a bounded invertible perturbation of identity.
The graph domain is unchanged and the actual outer-fixed inverse
exists uniformly. For F=b+c and f=beta_u F,

    sqrt(lambda)||v||<=sqrt(2)/(1-d)||f||,
    sqrt(lambda)||B0 v||<=(1-d)^-1||f||,
    sqrt(lambda)||eta||_(dmu_ell)<1.006||h||_(dmu_ell). (7)

Strict estimates concern nonzero h; zero returns to zero.
The seed R0<=sqrt(5) and (1) imply R<2.25sqrt(lambda).
Also (5)-(7) give the regular polar multiplication estimate

    ||A eta||<1.008||F||<.713||h||.               (8)

This bounds the product A eta, not the singular coefficient A alone.
Using (3), the component inequalities following (2), and (7)-(8),

    ||S_c||<2.6||h||, ||T_c||<3.7||h||,
    ||(S_c,T_c)||<4.6||h||.                      (9)

All norms in (7)-(9) are actual physical L2 norms. In particular this
is a uniform candidate return bound, not just a seed statement.

## 3. First derivative and cancellation of polar coefficients

The low-jet ledger in the scalar reference proof sharpens the needed
second derivatives without any new source estimate. Its core bounds
give |u_y|<.00025, |u_yy|<6.2. Its exact finite-radius formula gives

    |S_yy|<6.3,

because (1.001*6.2+22*2e-16+1.001*7*(3e-7)^2)/.999
plus the square of the first-log derivative is below 6.3. The
retained outer jets are smaller than these reserves for t>=40.
Thus these three bounds are global. Evenness at the bolt implies
|ell u_y|<=6.2y and |S_y|<=6.3y. Since A0<=1/y+3.25,

    |A0 ell u_y|<6.2+.0008125,
    |A0 S_y|<6.3+.00325.                         (10)

Define E1=A'-R'+A(A-R) and E2=R'+AR. Their exact transformations are

    E1/lambda=exp(-2w)[E1,0/lambda+ell u_yy+A0 ell u_y],
    E2/lambda=exp(-2w)[E2,0/lambda+S_yy+A0 S_y].  (11)

The seed coefficient bounds give |E1,0|<39lambda and
|E2,0|<27lambda. For example E1,0=F0,rr/2-F0,r/r+F0/r^2
and E2,0=F0,r/r-F0/r^2. Equations (10)-(11) yield

    |E1|<52lambda, |E2|<41lambda.                (12)

Transversality in (2) cancels the dangerous A b and A c terms in
the differentiated return. The actual identities are

    T_c'=-(3/2)c'+2R b-E1 eta,
    S_c'=c'/2-R(b+c)-E2 eta.                     (13)

These do not bound separate singular polar terms. Applying (7),(12)
to (13), with ||h'|| the L2 norm of the radial derivatives of its
relative components, gives

    ||T_c'||<.75||h'||+57sqrt(lambda)||h||,
    ||S_c'||<.25||h'||+44sqrt(lambda)||h||,
    ||(S_c',T_c')||<.8||h'||+73sqrt(lambda)||h||. (14)

In particular ||T_c'||^2<6500[||h'||^2+lambda||h||^2].
The proof extends by the regular TT H1 closure. These are derivative
bounds for the return of an INPUT tensor, not a bound on the unknown
physical metric correction or its cold-weighted jets.

## 4. Form-level logarithm versus the full operator domain

There is also a finite, uniform form-level mapping statement. Use
the seed physical tensor/scalar Sobolev spaces H^k with regular disk
charts at the bolt, fixed compact charts, and physical rescaled end
charts, together with the explicit measure identifications in (1).
Consider either the identified York TT projector, or the disk-return
pair map (3) on tracefree radial tensors with regular tensor parity.
The latter formula is defined even before imposing transversality;
equation (13) and the derivative bound (14) use transversality.

For k=2,3 both maps have finite uniform bounds K_k from H^k to H^k.
Here uniform means over the stated finite comparison path, not over
arbitrary metrics. The reason is the retained C5 coefficient control:
the global York theorem gives the H4 tensor map with a regular H5
vector, while (4) gives the return vector one derivative more than its
source. Its derivative recurrence through four involves only metric
derivatives through five. Near the bolt use the regular polar tensor
parities and its Volterra expression; on outer rescaled charts use
the retained five conormal jets. The lower-order inverses have the
uniform gaps already proved. The integer Sobolev estimates therefore
have uniform finite constants; interpolation gives k=2,3. This does
not provide new numerical values for K_2,K_3 or make them small.

The logarithmic interpolation can be stated without assuming any
extra coefficient derivative. For the compatible pair (H2,H3), set

    K(t,f)^2=inf_(f=f0+f1)(||f0||_H2^2+t^2||f1||_H3^2),
    ||f||_formlog^2=||f||_H2^2+integral_0^1 K(t,f)^2 dt/t. (15)

For either map L, K(t,Lf)<=max(K_2,K_3)K(t,f), by applying L to
each decomposition. Thus (15) is preserved with bound
max(K_2,K_3). On a flat chart, with nu=sqrt(1+|xi|^2), minimization
and the exact t integral give

    ||f||_formlog^2
      =integral nu^4[1+(1/2)log(1+nu^2)]|fhat(xi)|^2 dxi,
    log(2+nu)<1+(1/2)log(1+nu^2)<2log(2+nu), nu>=1. (16)

The same local norm equivalence holds for fixed chart definitions.
This is H2 with a SQUARE-ROOT logarithmic multiplier, the frequency
regularity of a fourth-order logarithmic QUADRATIC FORM. It is not
H4 with a full logarithmic multiplier, which is needed for the
operator-domain isomorphism. The latter remains unproved for the
finite comparison. Nor does (15) identify the complete quantum form
domain: outer position-log weights, actual nonlocal kernels and their
contacts still require their own estimates. The form-level kinematic
map is available for that work without first asserting the stronger
operator-domain map.

## 5. Actual canonical LLL disk Hessian at the finite comparison

The retained canonical HH/reflecting Polyakov allocation has central
charge q, as in [the seed contact proof](stationary-seed-quantum-form.md).
For the exponential DISK conformal tangent T_c, its exact cocycle is

    Gamma_hat,LLL,disk''=-q/(12pi) integral |grad T_c|^2 dvol_ghat.

Two-dimensional conformal invariance and the actual four-volume
measure make this

    Gamma_hat,LLL,disk''
      =-q/(48pi^2) integral r_ell^-2 |T_c'|^2 dmu_ell
      >=-6500q/(48pi^2 B^2) E_h,
    E_h=||h'||^2+lambda||h||^2.                   (17)

Here r_ell>B on the original finite path. The compact-outer regular
core has no boundary cocycle contribution to this quadratic term;
the physical finite-energy extension is its H1 closure. The selected
state, temperature, boundary projector and normalization are fixed.

Since Gq/B^2=(q/C)(GC/B^2)<400(8.5e-33)=3.4e-30 and pi>3,

    Gamma_hat,LLL,disk''> -1.3e-27 E_h/(8pi G).    (18)

The numerical comparison in (18) is to the EXPLICIT radial H1 energy
E_h/(8pi G), not to an unproved candidate classical TT lower bound.
It controls this component of the actual returned disk Hessian.
W4, massive angular modes, reference subtraction, density and the
full action's off-shell affine/gauge contacts are not absorbed into
the LLL term or dropped. Thus (18) does not prove full TT positivity.

## Frozen verification and remaining scope

Freeze this note and checker before execution. Reconcile the York,
global TT and seed contact receipts and their common input identities.
Verify the proper-frame return transformation, exact norm identities,
first-order form integration, E1/E2 cancellation, differentiated TT
return and all value/gradient/LLL constants. Check the quadratic
K-functional minimization and logarithmic weight inequalities, and
independent regular-bolt return examples. Fixtures do not evaluate the
physical determinant, signed moment or nonlinear solution. Continuum
domains, interpolation and weak inverse statements remain ordinary
proofs. Preserve every run and failure. No Lean campaign is selected.

The full quantum TT form, mixed/contact blocks and actual coupled
range still require proof. The stronger H4-log map remains open;
the form-level statement does not close that task. Pi_cut remains
fixed at the seed, with previous correction bounds unchanged.

Remaining to-do list: massive/Weyl/reference TT and mixed/contact forms,
outer position-log control, full operator-domain mapping and coupled
range; fixed projection direction and cold/higher inverse/source bounds;
nonlinear preservation and projected solve; signed compatibility;
exterior/state, charge and return; opening, stability, transport,
preparation, separate creation and full-result audit.
