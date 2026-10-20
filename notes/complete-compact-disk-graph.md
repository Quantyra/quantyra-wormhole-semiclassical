# A compact graph estimate for the complete disk Hessian

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary argument
and verification contract before execution. **Complete local disk
operator; global TT inverse and physical solution remain open.**

The compact fourth/log estimate now includes the full retained
Einstein-Maxwell, Weyl, LLL and massive action. The reference
bookkeeping action cancels exactly, off shell. Both disk equations
are estimated together, with explicit physical source conversion.

## 1. Statement and fixed identification

Use the analytic seed and the charts of
[the compact massive graph](compact-massive-logarithmic-graph.md):
center Y>=1/4, input support V={|z|<1/256}, observation
O={|z|<1/128}, gbar=exp(2sigma)|dz|^2 and |sigma|<.027.
The smooth bolt may be the center. Put

    theta=S+T/3, X_4=||theta||_H4+||T||_H4log,
    U_j=||theta||_Hj+||T||_Hj, j=0,3,
    F_O=||L_theta(theta,T)||_L2(O)+||L_T(theta,T)||_L2(O).

All chart norms use Euclidean measure and unitary Fourier
multipliers (1+|xi|^2)^(j/2); H4log includes log(2+|xi|).
L is the Hessian of the COMPLETE Euclidean disk action per fixed
gbar area, in coordinates (theta,T). It includes the density and
exponential-coordinate contacts. In particular it is not the
derivative of a moving-volume gradient with that contact omitted.

For every smooth pair supported in V,

    X_4 <= (2000/C) F_O + 2e96 U_3,
    X_4 <= (5000/C) F_O + 1e400 U_0.                (1)

This is a local a priori estimate, not an inverse statement.
The large value constant is not a nonlinear contraction estimate.
Applying it to a cutoff of a global solution requires the full
commutator, which is not evaluated here.

## 2. Regular seed jets, including at the bolt

Write u=B/r, omega=-log u, C=B^2/ell_o^2, x=lambda B^2,
g=G lambda and H=1e6. On O the preceding graph proof gives
1/28<u<1 and 0<omega<4. Put

    J(u)=u^2 F/C=(x/C) f(u),
    f=1+nu u^2-2a u^3+b u^4-2e u^5+e u^6.

The retained bounds nu<.001, a<2.002, 0<b<3.002,
|e|<1e-15 and x/C<1 imply

    0<=J<9, |J'|,|J''|,|J'''|<1e4.

There is a useful exact two-dimensional identity,

    |nabla u|^2=J,  Hess_gbar u=(J'/2)gbar.        (2)

Indeed u_Y=sqrt(J) and the canonical time warp is proportional
to sqrt(J); its two Hessian entries agree. The identity extends
smoothly through the bolt. No estimate of the singular polar
connection is used there.

In the conformal chart,

    partial_ij u=(J'/2)exp(2sigma)delta_ij
                                  +Gamma^k_ij partial_k u.

The retained analytic chart gives exp(2|sigma|)<2,
|D^j sigma|<=j!1000^j through four, and component Christoffel
bounds 100,25600,1e7 through two. Differentiating this displayed
identity twice gives the following coordinate component bounds:

    |Du|<4, |D2u|<2e4, |D3u|<1e8, |D4u|<2e11.

For example the third-derivative budget is
4e4+2e7+204800+4e6<1e8. The fourth-derivative budget is

    1e4(16+2e4)+1e4(4)(4000)+1e4(1.6e7)/2
       +2[1e7(4)+2(25600)(2e4)+100(1e8)] <2e11.

Finite Faà di Bruno for -log u, using
|(-log u)^(k)|<=(k-1)!28^k, now proves, for every component,

    |D^j omega|,|D^j sigma|<=j! H^j, 1<=j<=4.     (3)

These are seed bounds, not regularity assumptions on an unknown
solution or bounds on a moving coordinate system.

## 3. Complete local allocation and fourth-order diagonal

Retain the exact split from
[the covariant seed](stationary-centered-covariant-seed.md):

    I_full=I0+W4+Gamma_LLL+Gamma_mass-Gamma_ref
          =I_classical+W4+Gamma_LLL+Gamma_mass.       (4)

The equality holds for all spherical fields before differentiation.
It therefore preserves all off-shell contacts. In particular the
potential in the remaining classical action is the PHYSICAL

    V_phys(X)=2+6lambda0 X-2Q/X,

with the original charge Q, not the seed comparison charge.

Here is an explicit classical formula useful for checking every
local coefficient. Set Z=exp(2omega), h=omega+S. In gbar units,
up to the compact-variation boundary integration,

    I_classical=-B^2/(4G) integral_gbar [
       exp(2h)(Rbar+6|nabla h|^2+4nabla h.nabla T)
       +(2/C)exp(2h+2T)
       +(6lambda0 B^2/C)exp(4h+2T)
       -(2Q/(B^2 C))exp(2T)].                       (5)

It follows from X R2+(nabla X)^2/(2X)+V_phys after the
two-dimensional Weyl transformation and one integration by parts.
The kinetic coefficient is six, and the mixed coefficient is four.
The second variation of the radius exponential is included.

For the Weyl quadratic form, use exactly section 8 of
[the curved TT proof](stationary-curved-quantum-tt.md), now with
dimensionless gbar derivatives. Its principal terms are

    [C/(20pi)] integral (Delta theta)^2
      -[C/(30pi)] integral (omega+1/6)(Delta T)^2.   (6)

For explicit coefficient verification, put
R=Rbar, U=Delta omega, V=|nabla omega|^2,
e0=1/C, k0=(1/15-q^2/3)/C^2,
d=Delta T, f=Delta theta, w=nabla omega.nabla S and
j=|nabla S|^2. The remaining quadratic form is

    -C/(2pi) integral [d D+f L+J0 j+W T w-7w^2/45+Z0],
    D=S(-R/15+e0/9)+T[(6omega+1)R/45+11V/135]-11w/135,
    L=T[-R/15+2U/5+7V/45]-7w/45,
    J0=R/30-7(U+V)/90-e0/18,
    W=-2R/15+14(U+V)/45,
    Z0=T^2[omega R^2/30+R(U+V)/15-U^2/5
             -7UV/45-7V^2/90+2k0 omega]+ST[-R^2/30+2k0].

Its highest total derivative order on the two tangents is three.
Consequently its Euler operator has order at most three; its
background coefficients use at most four derivatives.
Equation (6)'s divergence-form completion has lower-order
terms as well. For a scalar v and b0=omega+1/6,

    Delta_gbar^2 v
      =exp(-4sigma)[Delta_0^2 v-4Dsigma.D Delta_0v
                         +(4|Dsigma|^2-2Delta_0sigma)Delta_0v],

    Delta_gbar(b0 Delta_gbar v)
      =exp(-4sigma)[b0 Delta_0^2 v
          +(2Db0-4b0 Dsigma).D Delta_0v
          +(Delta_0b0-4Db0.Dsigma
                  +b0(4|Dsigma|^2-2Delta_0sigma))Delta_0v].

All those terms are retained in the coefficient ledger below.
The canonical LLL contribution is exactly
-q/(12pi) integral|nabla T|^2, hence its operator is
q Delta_gbar/(12pi). Its fixed HH domain is unchanged.

The only nonlocal entry is the complete H_mass in the T equation.
Equations (4)-(6) thus yield the EXACT decomposition

    L_theta= [C/(20pi)]exp(-4sigma)Delta_0^2 theta+R_theta,
    L_T=H_mass T-[C/(30pi)]b0 exp(-4sigma)Delta_0^2 T+R_T. (7)

The fourth-order cross entries vanish in (theta,T) coordinates.
R_theta and R_T contain every remaining local contribution.

## 4. Explicit lower-coefficient budget

The following generous ledger can be checked by finite polynomial
differentiation. Put E=exp(-2sigma), Z=exp(2omega), so
E_i=-2E sigma_i and Z_i=2Z omega_i. Use (3) and

    E<2, Z<28^2, |omega|<4, e0<1, |k0|<60000,
    lambda0 B^2/C<1, Q/(B^2 C)<2, 1/g<2e33.        (8)

For the charge bound, B^2=z_star b0^2, C=C0 z_star/u_star
and Q=b0^2+3lambda0 b0^4 give
Q/(B^2 C)<u_star/z_star^2<1.01/.99^2<2.
Also lambda0/lambda<2 and C/x>5.9.
The lower bound g>.999*2000^-10 proves the last inequality.
These are retained benchmark bounds, not free choices of
physical parameters.

Expand the quadratic densities in (5)-(6), set S=theta-T/3,
and apply one HALF of their Euler operator, then multiply by E
to convert Euclidean density to gbar density. This half is
required because the displayed quadratic forms are I''[v,v],
not the quadratic Taylor coefficient I''[v,v]/2.
For (5), differentiate twice before applying this convention.

Subtract the two principal entries in (7). The resulting
coefficient sum is bounded as follows. Each background jet of
order j is assigned the weight j!H^j, each tangent jet weight
one, and E,Z,omega and constants the bounds in (8). Sum the
absolute monomial coefficients over BOTH output rows.

| Remaining allocation | Coefficient-sum bound |
| --- | ---: |
| Complete lower Weyl operator, including (6)'s completion | 1e42 C |
| Physical classical Hessian, including the radius contacts | 1e59 C |
| Canonical LLL Hessian | 100 C |

The Weyl count has background total order at most four and
tangent order at most three. In the classical count the
background order is at most two and the tangent order at most
two; its factor B^2/G=(x/g)<C/g is bounded by (8).
These finite polynomial bounds are checked on the complete
operators, not on selected entries or sampled seed points.
For each tangent derivative of order at most three its
Euclidean L2 norm is bounded by that tangent's H3 norm. Thus

    ||R_theta||_L2(O)+||R_T||_L2(O) <=1e60 C U_3.   (9)

For source context, the retained charge convention is independently
given in [Bintanja--Freivogel--Rolph](https://arxiv.org/html/2308.00871v2#S2.SS1.SSS4).
The LLL cocycle normalization is the one already fixed using
[Bilal--de Lacroix--Erbin, equations (1.3)-(1.4)](https://arxiv.org/pdf/2109.03637).
Both were inspected 2026-10-20. Neither paper supplies (9).

## 5. Coupled graph estimate and physical source conversion

For supported theta, the elementary Fourier estimate
||theta||_H4<=2||Delta_0^2 theta||_2+2||theta||_2
and exp(4|sigma|)20pi<80 give, from (7)-(9),

    ||theta||_H4 <=(200/C)||L_theta||_L2(O)+1e64 U_3.

The preceding massive graph, in its intermediate form and with
b=b0, gives

    ||T||_H4log <=(1600/C)||L_T||_L2(O)
                         +1600e60 U_3+1e96||T||_H3.

Adding proves the first inequality of (1). For both Fourier
norms, the retained interpolation bound gives
U_3<=epsilon X_4+(2/epsilon)^3 U_0.
Take epsilon=1/(4e96), absorb one half, and enlarge the
forcing/value constants to 5000/C and 1e400. This proves (1).

Here is the output normalization in the original physical disk
variables. If f=(f_S,f_T) is the full fixed-volume Hessian output
paired against (S,T), then

    dmu_4=(4pi r^4/C)dvol_gbar,
    (L_theta,L_T)=(4pi r^4/C)(f_S,f_T-f_S/3).       (10)

On O, B<r<28B and exp(|sigma|)<1.1. The two-row L2 sum is
bounded by 5/3 times the original pair norm before the measure
conversion. Consequently

    F_O <=1e4 B^2 C^-1/2 ||f||_L2(dmu_4;O),
    U_0 <=C^1/2 B^-2 ||(S,T)||_L2(dmu_4).          (11)

The weighted global source norm Fdisk=||u^-2 f||_L2(dmu_4)
also bounds the first norm in (11), since u<1. For supported
directions (1) therefore implies

    X_4 <=5e7 B^2 C^-3/2 Fdisk
               +1e400 C^1/2 B^-2 ||(S,T)||_L2(dmu_4).

This is a source/value graph bound, not a global inverse.
For a general solution the cutoff commutator appears in f and
must be estimated; (10) does not remove it.

## 6. Frozen verification and remaining work

Before execution freeze this proof and its checker. Independently
differentiate the complete integrated Weyl density twice and
compare its quadratic polynomial with (6). Independently vary
the physical classical density and retain the second radius
variation. Apply the two-dimensional Euler operator to both;
verify their principal diagonal, order bounds and every weighted
coefficient. Verify the regular u/omega jet ledger, benchmark
conversions, Fourier absorption, density/source transpose and
physical norm factors. Pin predecessor receipts and exact inputs;
preserve every run and failure.

This completes a local graph estimate for both equations of the
full retained disk Hessian. It does not give its global kernel,
range, inverse norm, cutoff commutators, higher York transfer,
nonlinear solution or physical stability.

Remaining to-do list: patch commutators and higher York maps;
outer fast fourth/log estimate and full inverse/residual
composition; nonlinear signed compatibility and all exterior,
opening, stability, transport, preparation and creation requirements.
