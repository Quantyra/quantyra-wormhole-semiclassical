# Assemble the complete cold action defect before applying an inverse

2026-10-22, E003/S169. Ordinary argument before supporting verification.
The full stable semiclassical wormhole objective remains incomplete.
This executes the cold part of the selected physical defect calculation.
It does not solve the full finite-domain linear or nonlinear problem.

## 1. Fix the rows, chart and exact defect identity

Keep the current c_R=1e22 action, G, q=2000^4, g_R=2000^-4,
mu0=1/b0, prepared positive temperature and paired HH realization.
Use the ORIGINAL-G reference split of
[the first correction](fixed-clock-seed-bolt-feasibility.md):

    I=I0+Gamma_rest, I0'(g0)=0,
    Gamma_rest=Gamma_hat,mass+Gamma_hat,LLL+W4+I_R2-Gamma_ref.

All boundary actions are repartitioned consistently. This is not the
effective-F_l comparison used to prove the full Hessian's form bound.
On the fixed canonical seed disk write

    g=exp(2phi)[A^2 exp(2b)gamma0+B^2 dOmega^2],
    phi0=log(r0/B), h1=(nu,u), nu=xi/r0,
    r_ell=B+(r0-B)exp(ell V), nu=(r0-B)V/r0.

The Hessian H in this note is the fixed-density action Hessian in the
AFFINE variables (phi,b). A fixed seed York return can subsequently be
applied with its correct dual map. Differentiating gradient rows along
the displayed path uses only its tangent (nu,u). Differentiating the
action twice along that nonlinear radius path would also produce an
acceleration contact; it is not the operation used here.

The complete reference constraints, angular Ward equation, zero-current
identity and fixed bolt/clock data imply, in these same affine rows,

    I0'' h1=-Gamma_rest'(g0),
    r1=I'(g0)+I''(g0)h1=Gamma_rest''(g0)h1.           (1)

This identity includes the equation-to-density contacts. It does not
replace H by I0'', discard the off-shell Ward remainder, or authorize
a second reference-constraint iteration. The sought error is still
e=-H^-1 r1 in the ACTUAL physical realization.

Normalize each action row by -2G/B^2 per seed canonical disk area.
The same factor multiplies H and r1 in every inverse equation below;
the physical correction and its admission norms do not change.
For a source with sphere-integrated energy rho and pressure p, its
shape row is

    S=2G r^2(p-rho)/(B^2 C),
    F_b=exp(2b)S, C=B^2/A^2.                         (2)

For the massive canonical source (2) is exactly S=2G Theta/(B^2C).
For the geometric trace defect e_disk of current-seed-full-trace-defect.md,
the normalized conformal row is -r^4 e_disk/(B^2 C). Thus both rows
below have the same action normalization. The arbitrary sign of an
overall equation convention has no role in the absolute bounds.

## 2. The R2/reference shape column with all its cold contacts

Set R=R4, O=Ric^t_t+Ric^r_r and w=R+12lambda. Primes here mean
areal r derivatives. For ds^2=-exp(2alpha)F dt^2+dr^2/F+r^2dOmega^2,
the exact curvature-square tensor gives

    H_sum=H^t_t+H^r_r
      =2RO-R^2+2F R''+2(F'+F alpha'+4F/r)R'.         (3)

This follows from the metric R2 variation including both Hessian terms;
the standard variation is reviewed in [De Felice--Tsujikawa, section 2.1](https://arxiv.org/html/1002.4928v2#S2.SS1).
No boundary or quantitative source estimate is imported from that paper.

At fixed areal r the actual reference correction obeys

    delta R=-T0, delta O=2G(p-rho)/r^2-T0,
    f'=(-f-2G rho)/r, F a'=G(rho+p)/r.

Differentiating (3), before taking any norm, yields

    delta H_sum=-2O T0+4GR(p-rho)/r^2
       -2F T0''-2(F'+4F/r)T0'
       +2f w''+2(f'+F a'+4f/r)w'.                    (4)

No T0 derivative has been identified with a nonlocal shape variation.
All terms in (4) are fixed-seed local curvature or the evaluated full
source entering h1. In particular F a' stays grouped at the bolt.

Let S_L be the seed shape source of I_R2-Gamma_ref, INCLUDING the
reference vacuum. Its exact fixed-disk density derivative is

    r1,b,L=-2chi_R(R/lambda) S_rest+2u S_L+xi S_L'
       +(16piG c_R r^4/(B^2 C))
          [2O T0+2F T0''+2(F'+4F/r)T0'
                         -2f w''-2(f'+F a'+4f/r)w'], (5)
    chi_R=16piG c_R lambda,
    S_rest=2G r^2(p-rho)/(B^2 C).

Here 2u S_L is the canonical density contact and xi S_L' is the
ENTIRE radius/pullback contact, including the r^4 factor in the R2
source. Omitting either would change the action defect. The leading
term in (5) contains the complete reference source S_rest, not an
independently varied quantum stress. The reference potential has no
fixed-areal metric variation in its mixed stress. This is why its
large seed terms enter only the already combined S_L and its derivative.

For explicit substitution put t=r/B=1+s, z=1/t, x=lambda B^2,
chi_R<4.91e-10, |d|<8.6e-19, x>8e9 and 5.9<C/x<6.1.
Use the exact polynomials P_rho,P_p in
[the directional calculation](current-seed-directional-cold-response.md):

    S_L,nv=(2x/C)z^-4(P_p-P_rho),
    S_L,vac=6(x/C)[(lambda0-lambda)/lambda]t^4.

The leading charge term of S_L,nv is constant in t. Absolute polynomial
coefficient sums give

    |S_L,nv|<1.16e-22,
    |d_t S_L,nv|<4e-28, |d_t^2 S_L,nv|<1e-27.       (6)

The fixed vacuum relation gives |S_L,vac|<2e-35 on s<=.001;
its first two t derivatives are bounded by the corresponding factors
4/t and 12/t^2. It is retained, not reset to zero. Its growing outer
counterpart must still be recombined with the Dirac vacuum there.

The improved Qrho,Qp<1.3e-10 and the actual horizon bound
|G rho(R_bolt)/(kappa R_bolt)|<2e-12 imply on h<=s<=.001

    |G rho|/(lambda B^2), |G p|/(lambda B^2)<1.31e-10 s,
    |S_rest|<8.9e-11 s, |d_t S_rest|<9e-11.           (7)

To see the first inequality, integrate the actual ordinary derivative
from t=1+h, use p(R_bolt)=-rho(R_bolt), kappa/(lambda B h)<7,
and (t^2-(1+h)^2)/2<=1.001(s-h). The horizon coefficient
1.4014e-11 h is smaller than the derivative coefficient. Thus (7)
does not assume a new source-response derivative or a product replacement.

Write f_n=f/(lambda B^2), A_n=F a_t/x, D_n=f_n,t+A_n+4f_n/t.
The known fbar/F<7.5e-10 and |s_R|/(Bh)<3.33e-10, together with
F/x<=7s^2 and F_t/x<=14s, give

    |f_n|<1e-8 s^2, |f_n,t|<4e-10 s,
    |A_n|<2.62e-10 s, |D_n|<1.1e-9 s,
    |f_n,tt|<2.7e-10, |D_n,t|<5.5e-10.              (8)

The derivative identities following (4), with Qrho,Qp, prove (8)
without bounding a'/F or differentiating a singular bolt quotient.
The retained trace polynomial has degree at most 12 and exact
zero-temperature product zero. Its first coefficient-derivative sum
is <7e-27. Consequently, including the tiny thermal offset,

    |T0|/lambda<8e-27 s,
    |d_t T0|/lambda<7e-27,
    |d_t^2 T0|/lambda<1e-25, |d_t^3 T0|/lambda<1.3e-24. (9)

For example the second and third derivative sums cost at most 13 and
182 times the first sum. The exact w polynomial gives coefficients
10,56,372,2856 times |d|/x through order three. Also |O|/lambda<13
and |d_t O|/lambda<30. All these are spatial seed bounds.

Use the current h1 path bounds |u|<=1e-6 a_s,
|u_chi|<=8e-6 a_s, |V|<1.2e-8, |xi_chi|/B<.005a_s,
|t_chi|<=2s, and a_s=s[1+log(.01/s)]. The xi_chi bound follows
from xi/B=t nu and the retained fixed-chart radius first jet
|nu_chi|<.0048a_s. No inverse unknown occurs here.

In units lambda^2, the bracket in (5) divided by s is bounded by

    2[13(8e-27)+7(.001)(1e-25)+(14+.028)(7e-27)
      +(1e-8)(.001)372(|d|/x)+(1.1e-9)56(|d|/x)]<4.1e-25.

Its ordinary t derivative is <4e-25 using (8)-(9). The prefactor
converts it by chi_R(x/C)t^4. The principal term in (5) is at most
2(4.91e-10)(12.001)(8.9e-11)s. The two contacts contribute
<2.4e-28 a_s, including the vacuum. Hence

    |r1,b,L|<1.06e-18 a_s,
    |r1,b,L(chi)-r1,b,L(0)|<5e-18 hL_h sqrt(chi), chi<=1. (10)

For the second line differentiate the grouped expression (5) along
chi. The principal derivative is <2(4.91e-10)(12.001)(9e-11)2s,
plus the w' term. The bracket, density and pullback derivatives fit
the remaining margin in 5e-18 hL_h, since a_s<=2hL_h there.
Integrate from zero and use chi<=sqrt(chi). This is a shape ROW
bound with all stated contacts, not the earlier reference-bolt functional.

## 3. Add the local quantum rows without changing their allocation

For W4 use the explicit conformal stress polynomials (1),(7) of
[the complete source calculation](comparison-path-complete-source-response.md).
They follow from the retained [Weyl anomaly action](https://arxiv.org/html/hep-th/0106151v1),
not from identifying orbit and independent radius variations. Their
coefficient counts are <1e4 through total derivative order four and
<1e5 after one spatial derivative, including the degree factors for
varying exp(-2sigma). The magnetic monomials are counted separately.

At the current cold seed, differentiated background jets obey j!100^j,
and the fixed-chart actual tangent jets obey

    u: 4e-6 a_s j!100^j,
    nu: 4e-5 a_s j!100^j, j<=5.

These follow from the path transfer already executed, with the same
finite normal-coordinate conversion. The rational geometry, product
center and fixed-scale restoration satisfy that calculation's guards.
The polynomial identity and its coefficient counts do not require
the old reference potential coefficients; that reference is handled
in (5). Use eta_G=GC/B^2<8.5e-33 and q^2/C^2<160000.
For EACH Weyl density/null profile the direct first variation is bounded by

    (eta_G/pi) 2(1e4)(4e-5)(4!100^4)a_s <6e-24 a_s,

plus a magnetic contribution <1e-32 a_s. The derivative of that
variation uses 1e5 and 5!100^5 instead and is <3e-20 a_s.
The shape profile is N-2rho. Its canonical density contact is 2u S_W;
the retained seed bound |S_W|<3e-15 a_s supplies that contact. Thus

    |r1,b,W|<1e-22 a_s,
    bolt oscillation <1e-18 hL_h sqrt(chi).             (11)

For the last line use the order-five bound, proper-to-chi factor <1.2,
a_s<=2hL_h, and differentiate 2u S_W. The seed order-five polynomial
bound, obtained with the old 100a_s seed difference, is <2e-13 a_s
for the shape gradient. On chi<=1, hL_h<81*2000^-10 makes the
two contact-product contributions negligible inside (11).

The canonical LLL scalar is k_L R_gamma with
k_L=Gq/(12pi B^2)<9.1e-32. The EXACT two-dimensional identity
delta[exp(2ell u)R_(exp(2ell u)gamma0)]=-2Delta_gamma0 u
includes its density contact. The fixed-chart second and third jets
therefore give

    |r1,b,LLL|<3e-32 a_s,
    bolt oscillation <1e-27 hL_h sqrt(chi).             (12)

These are local row statements. The thermal LLL state and boundary
ensemble remain fixed; the trace identity is not a determination of
an independently prescribed pressure or a changed boundary state.

Finally the already executed actual-domain massive contribution,
including its fixed-density contact, is

    |r1,b,mass|<1.72e-15 a_s,
    bolt oscillation <1.72e-14 hL_h sqrt(chi).          (13)

The entire subtracted angular tower, actual finite paired domain and
common finite restoration are retained in (13). The parameter derivative
exists in the weaker quarter-power target; the derivative itself has
the displayed square-root modulus. No square-root norm convergence or
spatial C3 response is newly asserted.

## 4. The other physical row and the assembled cold defect

The full trace defect is already known exactly, including its off-shell
pullback. Its density factor agrees with the present chart because

    a+xi'+2xi/r=2u+4nu.

Use (8)-(9) directly in that trace identity, instead of dividing its
whole-exterior weight by a_s. With epsilon=G/(60pi), D_R lambda<2.946e-9,
|K_T|<100lambda and the retained bounds

    |A_f|<=490epsilon lambda+316D_R lambda |d|/x,
    |A_Q|<=535epsilon lambda+56D_R lambda |d|/x,

the AREAL trace defect is bounded by

    |e_T|/(lambda s)
      <=(2.946e-9)[7(.001)(1e-25)+(14+.014)(7e-27)]
        +100epsilon lambda(8e-27)
        +|A_f|(1e-8)(.001)+|A_Q|(2.62e-10).

Add the full density/transport contacts

    [(2e-6)(.00331)+(4.8e-8)(.001)](8e-27)
                                  +(1.2e-8)(7e-27).

Their sum is <5e-34. Since the normalized trace row has factor
(x/C)t^4, it follows that

    |r1,phi|<8.6e-35 s.                              (14)

One chi derivative of the SAME exact identity uses only (8)-(9),
rho',p', the retained u_chi,nu_chi and xi_chi. The dominant transport
term is (|xi_chi|/B)|d_t T0|/lambda<3.5e-29 a_s. All other
terms, including Delta T0's derivative and the local coefficient
derivatives, fit below 1e-28 a_s in e_disk/lambda. After multiplication
by (x/C)t^4 and integration on chi<=1, the trace-row oscillation is
<1e-27 hL_h sqrt(chi). This deliberately loose trace allocation is
already far below the shape row. No further trace optimization is selected.

Combining (10)-(14) proves for the COMPLETE two-row cold defect in the
fixed physical disk action coordinates

    max_i |r1,i|<1.73e-15 a_s,
    max_i |r1,i(chi)-r1,i(0)|<1.73e-14 hL_h sqrt(chi). (15)

This closes the cold VALUE/bolt-oscillation source assembly for the
actual full action, not only the canonical massive column. It includes
the local R2/reference, Weyl, LLL and canonical density/pullback contacts.
It does NOT identify these rows with a TT dual source without the fixed
York transpose, prove C3 response bounds, or discard the global tail.

## 5. Quantitative consequence for the selected physical solve

The remaining inverse is global even when its input is supported in
the cold region. To state a useful numerical target without relabeling
the old C3 constants, fix a smooth 0<=eta_c<=1 equal to one for
s<=.0005 and zero for s>=.001. Let r_c=eta_c r1 and r_e=(1-eta_c)r1.
At chi<=1 the cutoff is identically one. Define the source envelope
U_c by |f_i|<=a_s and bolt oscillation <=10hL_h sqrt(chi), with
support in this cold patch. Define G_c as the operator bound from
THIS envelope to the retained physical cold return of H^-1, including
the fixed disk/York conversion and physical boundary realization.
Then (15) gives the exact conditional inequality

    ||return_cold(H^-1 r1)||
       <=1.73e-15 G_c+||return_cold(H^-1 r_e)||.        (16)

One sufficient allocation of the <1.667e-9 difference target is

    G_c<=578034.6820809248,
    ||return_cold(H^-1 r_e)||<6.67e-10.               (17)

Both quantities in (17) are UNKNOWN, not measured upper bounds.
The ten in U_c defines this stated source envelope; changing it
would change G_c reciprocally. No improvement is obtained by that
normalization. G_c is NOT the earlier K_c for seed C3 sources, and
the classical envelope gain 443..469 is not its value or upper bound.
It is only a scale comparison explaining why (15) is worth carrying
into the actual solve. A numerical or analytic bound for that full
return is still required. The exterior/mass/clock contribution in
(16) can be decisive and has not been assigned zero.

Retain the full coefficient norm Y0<=5e-9, rho=1e-8, Z0<=.25,
Z2<=5e7, all raw/weighted C5 admissions and D_ell embedding <=.025.
The cold value result does not close those conditions. The stationary
form inverse and physical fast phase are available within their
ordinary-proof scope, but their quantitative weighted return is not.
Matter boundary projectors alone do not determine gravitational
constraint data or retarded stability. Compatibility, matching,
semiclassical validity and all-sector stability remain required.

Decision: keep the closed audit's direct full physical defect solve.
The cold rows (15) are now assembled; do not further optimize them.
The SINGLE next calculation remains H e=-r1 with the full actual
domain and a weighted remainder, now with (15) as a known input and
(17) as a new explicit conditional allocation. Compute the exterior
contribution and full return together; do not replace them by a
product inverse or a source-only reference iteration. That actual
full solve is still UNEXECUTED. S169.2-S169.5 remain unchecked and
the full stable-wormhole objective is not reduced to this milestone.

## 6. Predeclared supporting verification

Freeze note, checker and source review before execution. Independently
derive the static R2 shape tensor, its first variation, and the reference
constraint/density cancellation. Check the exact local polynomial and
trace coefficient bounds, the grouped bolt estimates, local Weyl/LLL
normalizations and all numerical sums. Use a manufactured smooth radial
family to compare the unsimplified derivative with (4)-(5), including
nonzero density and pullback contacts. The fixture is not an evaluated
quantum determinant. Retain every run, source snapshot and raw failure;
revalidate inherited hashes. No Lean is selected. Numerical checks do
not certify the infinite-tower continuum or actual inverse arguments.

Remaining to-do list: full finite-domain defect return and weighted
admission, including the exterior; nonlinear compatibility and solution;
physical matching, retarded boundary realization, semiclassical validity,
all-sector stability and every remaining stable-wormhole requirement.
