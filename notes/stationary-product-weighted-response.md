# Bolt-regular product reconstruction and weighted cold response

2026-10-17 HST, E003/S139 supporting S074/S075. Ordinary proof before
focused verification. **PROGRESS; actual exterior and nonlinear control open.**

The complete-H2 product inverse has an explicit invariant kernel,
including the whole near-zero pole cluster. For sources supported in a
cold ball of radius at most 84, its normalized cosh-weighted C0, L2 and
specified radial C2 norm is below 1.3e8. It retains the thermal-gap
factor for a forcing bounded by h cosh(chi).

This continues the [complete product contour theorem](stationary-product-contour-inverse.md).
All its actual parameters, exact matrix M(p), physical normalization
gamma0=(8pi G A^2)^-1 and positive Stieltjes function H(p) are retained.
The [reduced cold response](stationary-reduced-cold-response.md) supplies
60<L_c<84, at the fixed prepared temperature.

The result is a bolt-regular particular product response with its
complete-H2 decaying tail. It does not impose arbitrary boundary data
on a restricted nonlocal throat operator. The actual variable-radius
seam and fixed mass/temperature return constraints remain open.

## 1. Normalized matrix and the precise inputs used

Write P2=-Delta_H2+2 and f(p)=M(p)^-1, where

    M(p)=[[-a(p),-2b p],[-2b p,p g(p)]],
    a(p)=c(p-r0), g(p)=alpha+pH(p),
    D(p)=4b^2p+a(p)g(p), det M(p)=-pD(p).

The companion estimates are

    5.9<c<6.1, .999<b<1.001, 1.3<r0<1.4,
    |alpha|<1.7e-8, eta_G=GC/B^2<8.5e-33,
    H(p)=beta+p integral dmu(tau)/(tau+p),
    eta_G/5<beta<2eta_G/3,
    supp mu is in (1404,infinity), mu>=0,
    integral dmu(tau)/tau<6e-5 eta_G.                    (1)

The exact inverse entries are

    f_SS=-g/D, f_ST=f_TS=-2b/D, f_TT=a/(pD).            (2)

For w=.99 or 1.01 and xi real, put

    nu=w+i xi, p=2-nu(nu+1), r=|p|,
    W_r=1+rH(r).

The companion all-frequency estimates give

    |D(p)|>=r[3+.75rH(r)], r>=.0299>1/34,
    |f_SS|<=3/r, |f_ST|<=3/(rW_r),
    |f_TT|<=12(1+r)/(r^2W_r), ||f(p)||<16000.           (3)

For |xi|>=4 these same entry bounds imply the useful simpler estimate

    ||f(p)||<=16/r<=17/xi^2.                            (4)

Indeed r>=xi^2-.0301, so r>15.96 and the largest absolute row sum
is at most [3+12(1+1/r)]/r<16/r. We use the Euclidean matrix norm.

## 2. The ordinary spherical density and the resolvent kernel

For radial smooth compact u on H2 define

    uhat(zeta)=2pi integral_0^infinity u(chi)
                     P_(-1/2+i zeta)(cosh chi) sinh chi dchi.

The ordinary Mehler-Fock inversion is

    u(chi)=(1/(2pi)) integral_0^infinity zeta tanh(pi zeta)
                uhat(zeta) P_(-1/2+i zeta)(cosh chi) dzeta.       (5)

Thus the inverse invariant-kernel density is
zeta tanh(pi zeta)/(2pi), not plain d zeta. This is the mu=0 case
of DLMF 14.20.11-12, using
Gamma(1/2+i zeta)Gamma(1/2-i zeta)=pi/cosh(pi zeta):
https://dlmf.nist.gov/14.20#vi
It agrees with the unit-radius primary transform in the contour proof.

We reconstruct the shifted formula through the resolvent, which keeps
this density automatically and avoids presuming a shifted Plancherel
identity. For p off [9/4,infinity), choose

    nu(p)=(-1+sqrt(9-4p))/2, Re nu>-1/2.

The exact invariant kernel of (P2-p)^-1 is

    R_p(d)=Q_nu(p)(cosh d)/(2pi).                       (6)

Here Q is the ordinary Legendre Q, not Olver's Q divided by
Gamma(nu+1). The distinction is explicit in DLMF 14.3.10:
https://dlmf.nist.gov/14.3.E10
To check (6), Q solves the radial Legendre equation, decays at infinity,
and has coefficient -log d at d=0. The flux of -Delta applied to
-log d/(2pi) is one. Its decay selects the complete-H2 L2 resolvent;
there is no extra bolt extension parameter. This also derives the
normalization in (6) independently of a spectral-density convention.
See DLMF 14.3.7 and 14.8 for the equation normalization and endpoints:
https://dlmf.nist.gov/14.3.E7
https://dlmf.nist.gov/14.8

An especially useful integral, in this ordinary Q normalization, is

    Q_nu(cosh d)=e^(-(nu+1)d) integral_0^infinity
      e^(-(nu+1)s)/sqrt[(1-e^-s)(1-e^(-2d-s))] ds.       (7)

For Re nu>-1 and d>0 it follows by Euler's beta integral, or directly
by inserting v=e^-s: the right side is

    e^(-(nu+1)d) B(nu+1,1/2)
                   2F1(1/2,nu+1;nu+3/2;e^-2d).

The quadratic hypergeometric representation is the same solution as
DLMF 14.3.7. Equivalently, differentiation shows the radial Legendre
equation and the leading coefficient
sqrt(pi)Gamma(nu+1)/Gamma(nu+3/2) at infinity; these determine Q.
Equation (7) fixes every coefficient used below.

## 3. An absolutely integrable complex-Q estimate

Let .99<=w<=1.01 and nu=w+i xi. Denote
A_d=(1-e^-2d)^-1/2. The amplitude of the s integral in (7) is
positive and decreasing. At |xi|<=1 its absolute integral is bounded
by A_d B(w+1,1/2)<=2A_d.

For |xi|>=1 split at s0=1/|xi|. Since
1-e^-s>=s/(1+s), on 0<s<=s0<=1 the amplitude is at most
sqrt(2) A_d s^-1/2. This part integrates to
2sqrt(2) A_d/sqrt(|xi|). On [s0,infinity), integration by parts
against e^(-i xi s), using monotonicity of the amplitude, bounds
the tail by twice its value at s0 divided by |xi|. The same upper
bound 2sqrt(2) A_d/sqrt(|xi|) applies. Consequently

    |Q_(w+i xi)(cosh d)|
       <=9 e^(-(w+1)d)/sqrt[(1-e^-2d)(1+|xi|)].         (8)

This explicit half power in xi is essential. Merely bounding Q by
Q_w would leave a logarithmically divergent contour integral.

For the low-p cluster below, |p|=.001 implies Re nu(p)>.999.
Using (7) without oscillation then gives the simpler estimate

    |Q_nu(p)(cosh d)|<=2 e^-1.999d/sqrt(1-e^-2d).       (9)

Indeed nu-1=-2p/[3+sqrt(9-4p)], so |nu-1|<.000334.

## 4. Locating every pole crossed by the contour shift

Let Gamma_w be the parabola p=2-(w+i xi)(w+1+i xi), oriented with
xi increasing from -infinity to infinity. It runs from upper right
to lower right around its left vertex. Its right-hand region is
the interior when closed far to the right in the positive orientation.

The proof of (3) extends to the region to the right of Gamma_.99:
on its compact part |p|<21 and |p|>=.0299 the perturbation of
4b^2 in D/p is below 5e-6 plus a negligible eta_G term; on the
remaining part the same right-half-plane sector estimates used in
the companion contour proof apply, with smaller |Im p| at fixed Re p.
Thus f is analytic in that region, which contains [9/4,infinity).

For clarity the entire larger right-hand region is described by
x=Re p>=y^2/(3.02)^2-.0301, where y=Im p. If x>=15.9699, then
|y|<=3.02sqrt(x+.0301), so |y|/x<.76 and |y|/(x-r0)<.83.
The companion sector proof applies throughout this region, giving
Re[D(p)/p]>3+r|H(p)|. If -.0301<=x<=15.9699, then |y|<=12.08
and |p|<21. Thus H(p)<eta_G in modulus. Outside |p|<.001,

    |a(p)alpha/p|<=6.1(1+1.4/.001)(1.7e-8)<.00015,
    |a(p)H(p)|<137eta_G.

Hence D does not vanish there. On |p|=.001, Rouche comparison
with 4b^2p gives exactly one zero of D inside. Together with the
explicit factor p in det M, this is the complete low cluster, counted
with multiplicity. It is the same Ward/physical cluster as in the reduced cold-response theorem.
There are no other poles of f between the contours or to their right.

On the cluster circle, (1)-(2) give

    |a|<8.547, |D|>3.99|p|,
    ||f(p)||<2.2e6.                                   (10)

The cluster is treated by a single Cauchy integral. If alpha is nonzero,
the separate Ward residue of f_TT is 1/alpha, and its cancellation with
the nearby physical residue must not be bounded separately. If alpha
vanishes, a double pole can occur. The circle bound (10) covers both
cases without dividing by alpha or assuming the two roots are distinct.

## 5. Exact Dunford reconstruction, including the cluster

For each real lambda>=9/4, elementary Cauchy integration gives

    f(lambda)=-(1/(2pi i)) integral_Gamma_.99
                      f(p)/(lambda-p) dp.              (11)

For a fixed lambda one closes the contour at Re p=R>lambda and
lets R tend to infinity. Along that closing segment the sector bounds
give f(p)=O(R^-1), and its length is O(sqrt(R)), so its contribution
vanishes. This is first a scalar spectral identity. No operator
resolvent on a closing segment through the unbounded spectrum is used.

Along the open parabola the actual resolvent norm is bounded by its
distance from [9/4,infinity). Together with (4), this makes (11)
an operator-norm convergent integral after replacing lambda by P2.
The ordinary spectral theorem then proves the operator identity.
Shifting to Gamma_1.01 gives

    f(P2)=-(1/(2pi i)) integral_Gamma_1.01 f(p)(P2-p)^-1 dp
         +(1/(2pi i)) integral_(|p|=.001,CCW)
                                    f(p)(P2-p)^-1 dp. (12)

The sign follows from the pole of (lambda-p)^-1 having residue -1.
The second term is exactly the sum of the two crossed pole contributions,
including the derivative-resolvent term if they coincide. Thus (12)
contains the ordinary spectral density, all crossed residues and the
regular-origin domain. It is not the assertion that the inverse in
(5) can be shifted with no correction or changed to plain d xi.

Using (6), the invariant matrix kernel K(d) of f(P2) is the sum of
the two corresponding Q integrals in (12). They converge absolutely
for every d>0. To bound the long contour, use (3)-(4),
|dp/dxi|<=2|xi|+3.02, and (8). The elementary integral estimate is

    integral_R ||f(p(xi))|| (2|xi|+3.02)/sqrt(1+|xi|) dxi
     <2(16000)(4)(11.02)
            +2(17)(2.755) integral_4^infinity xi^-3/2 dxi
     <1.5e6.

The factors 1/(2pi) from Cauchy integration and 1/(2pi) from the
Green kernel, with pi>3, therefore give a bound 4e5 for that part.
For the circle use its length 2pi(.001), (9) and (10). Its coefficient
is at most (.001)(2.2e6)(2)/(2pi)<1000. Altogether,

    ||K(d)|| <= [4e5 e^-2.01d+1e3 e^-1.999d]
                                        /sqrt(1-e^-2d).          (13)

The majorant is locally integrable against the H2 area element, since
its singularity is only d^-1/2. It is also integrable at infinity.
Thus the absolutely convergent kernel in (13) is the actual bounded
complete-H2 inverse on compact sources, not just a formal symbol.

## 6. Explicit physical weighted value and L2 bounds

Fix 0<L<=84. Let B_L be the geodesic ball about the regular bolt,
and let omega(chi)=cosh chi. The logarithmic derivative of omega is
bounded by one. For any x,y in H2,

    omega(y)/omega(x)<=e^d(x,y),
    omega(x)/omega(y)<=e^d(x,y).                        (14)

Let F be a two-component source supported in B_L, extend it by zero
to complete H2, and define U=f(P2)F. For x,y in B_L, d(x,y)<=2L.
Equations (13)-(14) give a row and column bound for the conjugated
kernel omega^-1 K omega. Indeed

    2pi sinh d e^d e^-2d/sqrt(1-e^-2d)
         =pi sqrt(1-e^-2d)<=pi.

Therefore both its row and column integrals over B_L are bounded by

    pi[4e5 integral_0^(2L) e^(-.01d)dd
               +1e3 integral_0^(2L) e^(.001d)dd]
     <=pi[4e7+2000L e^(.002L)]<1.3e8.                 (15)

For the last inequality use L<=84, e^.168<1.2 and pi<22/7.
Schur's test and the same row bound give respectively

    ||U/omega||_(L2(B_L))<=1.3e8 ||F/omega||_(L2(B_L)),
    ||U/omega||_(L-infinity(B_L))
                      <=1.3e8 ||F/omega||_(L-infinity(B_L)).     (16)

The L2 measure is the actual H2 area element. Using the physical measure 4pi B^2 A^2 dvol_H2 in both L2 norms
does not change the norm bound. For the physical Hessian H_red=gamma0 M(P2), replace
F by gamma0^-1 times the physical source. Thus the physical response
constant in (16) is 1.3e8(8pi G A^2).

In particular a normalized coupled source with |F(chi)|<=a h cosh chi
and support chi<=L produces the concrete radius/shape estimate

    |S(chi)|,|T(chi)|<=1.3e8 a h cosh chi, 0<=chi<=L.   (17)

This keeps the thermal gap factor h. It is an estimate for both entries
of the exact quantum product inverse, not a scalar-alone throat barrier.
The constants are conservative; the finite near-root cancellation is
preserved by (12), although (15) does not optimize its neck-length cost.

## 7. Bolt-regular C2 spaces and derivative control

For a smooth radial vector V define the dimensionless radial jet norm

    ||V||_(X_L^2)=sup_(0<=chi<=L) [1/cosh chi]
         max{|V|,|V'|,|V''|,|coth chi V'|}.             (18)

At the bolt the last expression is its regular even limit. Define the
source space as the completion in (18) of smooth radial sources compactly
supported in B_L; their zero extensions are C2 and have no boundary
jump. The output space consists of regular radial C2 fields on B_L with
the same norm, without prescribing an outer boundary value.

The complete inverse f(P2) commutes with hyperbolic isometries. Its
commutation with infinitesimal generators holds first on smooth compact
sources, by differentiating the unitary group action and then using
the bounded inverse. In polar coordinates use the two smooth boost
Killing fields

    Z1=cos theta partial_chi-coth chi sin theta partial_theta,
    Z2=sin theta partial_chi+coth chi cos theta partial_theta.

For radial V, each Z_i V is bounded in magnitude by |V'|. Also

    Z1^2 V=cos^2 theta V''+sin^2 theta coth chi V',
    Z2^2 V=sin^2 theta V''+cos^2 theta coth chi V',
    Z1 Z2 V=sin theta cos theta [V''-coth chi V'].

Each displayed second derivative is bounded by the maximum of the
last two jets in (18); the same is true with indices reversed on a
radial function. These derivatives preserve the compact source support.
Apply (16) to F, Z_i F and Z_i Z_j F. At theta=0 the corresponding
output jets are U, U', U'' and coth chi U'. Hence

    ||f(P2)F restricted_to B_L||_(X_L^2)
                                <=1.3e8 ||F||_(X_L^2).          (19)

Smooth compact input gives a smooth radial output at the bolt. One
direct domain check uses the physical spectrum p>=9/4: section 8 of
the companion contour proof establishes the same entry bounds there as
in (3), and these imply ||p f(p)||<24. Thus f(P2) maps L2 to Dom P2;
applying the same
argument to powers of P2 shows that smooth compact sources yield
smooth output. Its radial expansion is regular and even; a logarithmic
bolt solution is excluded. By (19) and the complete-product L2 inverse,
the same statement extends to the completed C2 source space as an exact
distributional solution of M(P2)U=F with the specified complete-H2 tail.

The derivative estimate commutes actual geometric isometries, so it
does not differentiate a singular contour integral term by term. It
does require the indicated C2 source jets, including the jets of any
source cutoff used to produce a compact forcing.

## 8. The finite-neck and actual-seed boundary that remains

The map in (16),(19) is exactly

    restriction_to_B_L composed_with M(P2)^-1
                                    composed_with zero_extension.

It selects the complete-H2 decaying tail. It is not an independent
inverse of the restricted nonlocal Hessian on B_L. For a local second-
order equation, one might add a regular homogeneous solution to adjust
an outer datum. Here the full product quantum operator samples both
sides of the seam, so changing the exterior continuation changes the
interior operator as well. No arbitrary Dirichlet or Robin datum at
chi=L is supplied by the preceding construction.

The remaining actual calculation is now concrete: write the true
coupled cold/exterior return problem using this particular product
inverse and its exterior tail, estimate the difference of the actual
variable-radius quantum kernel and the product kernel in the weighted
jet/source norm, and solve the finite-dimensional mass/temperature
matching system with that kernel correction included. A source cutoff
must be accounted for through (18), not treated as a free operation.
The fixed reference radius R=B(1+h), laboratory temperature and physical
mass return constraints have not been changed or automatically solved.

Equations (13),(16),(19) close the product spatial reconstruction and
give explicit weighted physical norms. They do not close the actual
quantum seam inverse or the nonlinear source-to-observable theorem.
All full-action signed source-moment and E003/S074-S080 physical gates
remain open as in the frozen milestones. No scientific checker or
numerical contour calculation was executed for this proof.

The shared [focused checker](../checks/stationary_product_weighted_inverse.py)
audits finite identities and rational reserves for the contour and
weighted-response proofs. Its [JSON receipt](../evidence/stationary-product-weighted-inverse.json)
and [raw log](../evidence/stationary-product-weighted-inverse.log) record
the actual execution outcome and exact input hashes. Analytic
continuation, infinite pole sums, spectral domains, contour deformation
and spatial reconstruction remain ordinary proofs.

Remaining to-do list: actual variable-radius quantum kernel comparison and coupled exterior
return/mass/temperature matching; nonlinear source-to-observable bound.
