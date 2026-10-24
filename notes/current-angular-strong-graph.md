# Full angular discretization in the current strong graph

2026-10-24, E003/S169. Ordinary argument before acceptance.
Previous turn: progress on the regulator inventory, not nonlinear closure.
The full stable-wormhole objective is ACTIVE and INCOMPLETE; audit CLOSED.

The complete curved nonlinear response is still unbounded. This note
closes a previously unallocated SEED component in precisely the current
H3 source graph: the full flat angular tower minus its matched logarithmic
normal model, on output z<.01 and input z'<.02. The bound is <3.038e-9.
It includes every angular level, all frequencies, the reflected kernel
and the prescribed finite fourth contact. It does not replace the actual
curved/thermal operator by this comparison, or prove its metric variation.

## 1. Exact allocation, domain and normal check

Retain sigma=4.005, x=log z, D=z partial_z, the actual full B_sigma and
the original shape-row operation z^-sigma(D-4)z^4. Let a=.01, b=.02.
The output chi is the existing cutoff, with derivative guards (1,3,10,60).
Its H3 multiplication allowance is 1+9+30+60=100.

The already specified flat entire-tower multiplier is

 K(p)=C kappa4 p^4/(30pi)
      +sum_(n>=1) d_n M_n^2/(2pi) integral_0^1
            [log(1+t)-t+t^2/2] d alpha,
 t=alpha(1-alpha)p^2/M_n^2,
 M_n^2=n(n+q)/C, d_n=2(q+2n),
 K_infty(p)=C p^4[log(p/(mu0 ell_o))-31/30]/(30pi),
 E(p)=K(p)-K_infty(p).                              (1)

The complete-sum proof in outer-flat-tower-normal-remainder.md gives
|E(p)|<=1000C(1+p^2)^(3/2). This is the exact matched multiplier,
including its finite polynomial terms; an off-diagonal bound alone
would not define it. In z coordinates set E_beta(p)=E(sqrt(beta)p).
Since 0<beta<1, the same Fourier bound is valid for E_beta.

For w=B_sigma^-1 g, put v=z^sigma w on the positive half. The operator
bounded here is

 R_ang g=-(2G/B^2) chi z^-sigma z^5 partial_z
                          E_beta,ev[1_(0<z<b) v].   (2)

The even response has direct PLUS reflected kernels, with no extra
half. The factor -2G/B^2 is the existing canonical action normalization.
Indeed its logarithmic coefficient is -GC beta^2/(15pi B^2), the retained
K-normal coefficient. G=1 in the numerical units; it is not an adjustable
normalization. The identity (D-4)z^4=z^5 partial_z is exact.

The input restriction in (2) is the existing near/far decomposition of
an operator, not a boundary condition. Its discontinuity at b is separated
from every observation z<a. The local Fourier arguments below use input
cutoffs strictly INSIDE b. All actual inputs beyond b, actual-minus-cylinder
propagation, and seagull/domain returns retain their previous allocations.
The exact remaining near operator is, by subtraction, the actual cylinder
response minus this flat component and the already allocated principal
and conversion terms. Its thermal, curvature and diagonal terms remain
unknown; none is dropped or declared equal to the flat model.

Before applying the inverse, check the full normal contact, not just its
leading frequency order. With P4=s(s-1)(s-2)(s-3), the even Fourier formula
gives a flat contribution proportional to

 P4[psi(s-1)+(pi/2)cot(pi s/2)+kappa+O+1/6]
                                     -s(s-1)(2s-5). (3)

Differentiating the remaining log z after the inner two derivatives
produces the last polynomial. The retained Weyl contact cancels both
that term and O. This is the exact identity in sections 3-4 of
outer-complete-disk-normal-comparison.md and
stationary-disk-and-normal-inverse.md. At s=3 the flat -6/z and Weyl
+6/z cancel. There is no missing finite term to amend in B_sigma.
The real-frequency graph and the original mass/quartic moments are
unchanged. A differentiated row does not eliminate the need to impose
those moments in the full solution.

## 2. Obtain derivative bounds from the actual complete kernel

Do not differentiate the real envelope for E. Instead use its defining
angular counting integral. The retained exact discrepancy is
|N(M)-2CM^2|<=4q+2sqrt(C)M. For the radial two-dimensional bubble

 B_M(r)=M^4[K0(Mr)^2-K1(Mr)^2]/(2pi^2),
 partial_M B_M=M^3[2K0(Mr)^2-K1(Mr)^2]/pi^2.

The integral representation [DLMF 10.32.9](https://dlmf.nist.gov/10.32.E9)
gives |K_nu(Mr)|<=K_nu(M Re r) for nu=0,1 and Re r>0. Therefore the
SAME Stieltjes integration by parts, now in this complex half-plane,
gives an analytic separated kernel with

 |k_2(r)|<=10[q/(Re r)^4+sqrt(C)/(Re r)^5].          (4)

The constants follow from the retained squared-Bessel moments. For
example the respective sharper coefficients are 16/(3pi^2) and 99/256;
both are below 10. These moments follow by squaring the positive integral
and evaluating its beta integrals, as in the cited original proof. The
absolutely integrable bounds justify the angular sum, time integral and
complex continuation on compact subsets of Re r>0. This is an argument
for this translation-invariant comparison kernel only, not an analytic
extension of a rough physical metric.

For real s>0 and |zeta-s|<=s/4,
Re sqrt(zeta^2+t^2)>=sqrt(s^2/2+t^2). Integrating (4) in the entire
time-separation line gives

 |k_1(zeta)| <=20pi q/s^3+(160sqrt(2)/3)sqrt(C)/s^4
             <100[q/s^3+sqrt(C)/s^4].

Cauchy on this disk consequently proves, for r=0,...,5,

 |partial_s^r k_1(s)| <=A_r[q/s^(3+r)+sqrt(C)/s^(4+r)],
 A_r=100 r!4^r.                                    (5)

The beta rescaling multiplies the two terms by beta and beta^(3/2),
respectively, so (5) remains a conservative bound for E_beta. Negative
s are supplied by evenness. All finite contacts are handled by (1) in
the near estimate below; they are not inferred from (5).

## 3. Near interactions: keep the seventh derivative in L2

Partition (0,a) into I_j=(h_j,2h_j), h_j=a 2^(-j-1), for EVERY j>=0.
Choose eta(t) supported in [1/2,3], one on [3/4,5/2], with beta(9,9)
transitions on [1/2,3/4] and [5/2,3]. This is C8. Put eta_j(y)=eta(y/h_j).
The supports overlap at most four times, and lie below 3a/2<b.
Thus eta_j v has a valid even H7 extension, independently of the
sharp near/far partition at b. No endpoint of this cover is a wall.

Exact Bernstein coefficient bounds for the cutoff derivatives are

 Aeta_0,...,Aeta_7 =
 [1,68,4352,522240,43868160,4562288640,
                              364983091200,32118512025600]. (6)

These deliberately generous constants are allowed because derivatives
of the cutoff are paired with the CORRESPONDING lower derivative of w,
not with the largest inverse constant indiscriminately.

Write F_n(X)=product_(m=0)^(n-1)(sigma+X-m)=sum f_nm X^m.
For n<=7, Leibniz and dy=y d(log y) give

 ||partial_y^n(eta_j v)||_2
 <=h_j^(sigma+1/2-n) sum_(i=0)^n binomial(n,i) Aeta_i
       T_(sigma+1/2-n+i) sum_m |f_(n-i),m|
                         ||D^m w||_(L2(supp eta_j;dx)),       (7)
 T_p=max_(1/2<=t<=3)t^p =3^p if p>=0, (1/2)^p otherwise.

The whole-line even extension costs sqrt(2). By Plancherel and (1),
for j=0,...,3,

 ||partial_z^(1+j) E_beta(eta_j v)_ev||_2
 <=1000C sum_(k=0)^3 sqrt(binomial(3,k))
                             ||partial_z^(1+j+k)(eta_j v)_ev||_2. (8)

The symbol index j in (8) is independent of the annulus label.
This keeps the highest derivative in L2 and uses no W3,infinity tube.
Define c_lj by c_00=1 and
c_(l+1),j=(5+j-sigma)c_lj+c_l,(j-1), with absent entries zero. Then

 (D-sigma)^l[z^5 partial_z E v]
             =sum_(j=0)^l c_lj z^(5+j) partial_z^(1+j) E v. (9)

On I_j the conversion from dz to dx costs h_j^-1/2. Combining (7)-(9)
leaves h_j^(4-k), k<=3; all powers are strictly positive. Bound them
by (a/2)^(4-k), sum all annuli, and use overlap four. The H3 binomial
weights and the output cutoff give the explicit finite expression

 eta_near=(2G/B^2) 100 sqrt(8) 1000C
   sum_(l=0)^3 sum_(j=0)^l sum_(k=0)^3
    sqrt(binomial(3,l)) |c_lj| 2^(5+j-sigma)
    sqrt(binomial(3,k)) (a/2)^(4-k)
    sum_(i=0)^n binomial(n,i) Aeta_i T_(sigma+1/2-n+i)
                    sum_m |f_(n-i),m| W_m, n=1+j+k.        (10)

Here W_m bounds ||D^m B_sigma^-1 g||_2/||g||_H3. Use the EXISTING full
inverse, with K=GC beta^2/(15pi B^2), L=10pi Fstar/(G lambda), delta=.005:

 Cb=[2/(KL delta sigma(sigma-3)),2/(KL sigma(sigma-3)),
                                      2/(KL sigma),2/(KL)],
 Cb4=max(Cb0,8/(KL),1/[K sqrt(L log6)]),
 W=[Cb0,Cb1,Cb2,Cb3,Cb1,Cb2,Cb3,Cb4].                (11)

The last four entries follow from the unit third-derivative coefficient
in the H3 norm of D^j B_sigma^-1g. No interpolation or sampled inverse is
substituted. In particular the large Cb4 occurs only at input derivative
seven; multiplying every term of (10) by it would be a substantial
avoidable loss.

## 4. Remote interactions: retain the image cancellation

On I_j, the complementary input vanishes unless y<=3h_j/4 or y>=5h_j/2.
Thus |z-y|>=max(z,y)/5. For an odd derivative r=1,3, the high-input
combination is

 k_1^(r)(z-y)+k_1^(r)(z+y)
             =k_1^(r)(y+z)-k_1^(r)(y-z).

The mean-value theorem and (5) give an EXTRA factor z. Without this
reflection cancellation the following far bound would not be justified
at sigma>4. Even r=2,4 uses the ordinary sum bound.

Set low_r=2 A_r 5^(4+r), and
high_r=low_r for even r, 2 A_(r+1)5^(5+r) for odd r;
e_r=4+r+(r mod 2)-sigma. These e_r are 6-sigma or 8-sigma, positive.
After multiplying by z^(5+j-sigma), writing dy=y d(log y) and using
y,z<=b, the positive Schur majorant is proportional to

 [q b^2+sqrt(C)b] low_r (y/z)^(sigma+1) for y<z,
 [q b^2+sqrt(C)b] high_r (z/y)^e_r for z<y, r=1+j.

The full-line log-coordinate convolution norms are therefore bounded by
low_r/(sigma+1)+high_r/e_r. Truncating the physical input/output sets
only decreases this positive bound. The entire far contribution is

 eta_far=(2G/B^2)100[q b^2+sqrt(C)b] Cb0
   sum_(l=0)^3 sqrt(binomial(3,l)) sum_(j=0)^l |c_lj|
                       [low_(1+j)/(sigma+1)+high_(1+j)/e_(1+j)]. (12)

No derivative of the annulus partition is taken: each true output
derivative is first split into near and far integrals on its annulus.
Their sum is the original distributional derivative. The resulting
global estimates extend by completion from smooth graph data. The
original scalar/shape triangular return has the same shape component,
so this shape-only bound also bounds this component in the coupled graph.

## 5. Execute, classify and decide

Using the current correlated q,C,B,beta,lambda,F_lambda and the full
normal inverse, interval arithmetic gives

 eta_near=3.037932454286561...e-9,
 eta_far =1.575630476581287...e-16,
 ||R_ang||_(H3->H3)<=eta_near+eta_far<3.038e-9.       (13)

The near input-derivative-zero and derivative-seven terms contribute
1.241347e-9 and 1.796423e-9, respectively. Other derivatives together
cost less than 1.63e-13. These identify the cutoff allowance and genuine
highest-derivative inverse as the main losses in THIS certificate.
Neither warrants optimization: their entire sum is already negligible
against the .1599164 nonlocal target. The old 1e30 epsilon comparison
was in a different graph and was unusable at the current collar; (13)
is a new strong estimate in the unchanged current units, not a rescaling.

This closes ONLY the near-input SEED angular-discreteness component.
The prescribed fourth contact is included, but unknown actual diagonal
contacts are not. All angular levels are summed; zero-winding continuum
time is a comparison component, and actual thermal terms remain in the
unestimated difference. There is no new bound for the first metric
variation of that difference, uniform nonlinear admission, or the full
original X-to-graph and preconditioned return maps.

Allocate the conservative 4e-9 to (13). The remaining curved/thermal and
other unallocated nonlocal SEED feedback must then be <.159916396;
M_A<=3.18e27 and M_V<=2.66e30 remain unknown separate contact targets.
For a graph ball of radius 1e-11 the conditional combined requirement is

 eta_curved_seed+1e-11 C3 < .159916396.               (14)

C3 must bound the COMPLETE metric variation, including any variation
of this split; (13) is not its value. Neither eta_curved_seed nor C3 is
known. Relating (14) to Z2<=5e7 still needs the original norm comparisons
and all other sectors. An excessive old upper allowance is not a lower
bound or a nonexistence theorem; the remaining actual constants here are
UNKNOWN. There is no demonstrated physical obstruction or solution.

The prior sufficient E_rest<2.123e-42<1.287e-40 enclosure is unchanged;
zero remains its enclosure center, not a measured signed projection.
No new full signed projection or actual global correction is computed.
At rho=1e-8, Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1 remain
unclosed; the allocations Y0<=5e-9, Z0<=.25, Z2<=5e7 are unproved.
Compatibility, physical matching, semiclassical validity/noise and
all-sector stability remain explicit obligations.

RESEARCH DECISION: retain the construction and the direct strong graph
method. Park the now controlled flat angular component. Its proof uses
translation invariance to move output derivatives and proves kernel
derivatives from a complex integral; neither operation may be applied
to the actual curved kernel merely by differentiating its real envelope.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): bound the actual
curved/thermal remainder and its first metric variation in the coupled
H3 graph uniformly on the admitted mass/log neighborhood, combining the
complete common-regulated vertices before norms and retaining the actual
finite-domain return. Test (14) and the original nonlinear norm comparison.
This remains the decisive feasibility calculation; (13) closes one
required strong component and does not replace that calculation.

Remaining to-do list: curved/thermal strong response and its diagonal
contacts, uniform complete metric variation, actual global correction,
original X/graph comparison and nonlinear admission/contraction, mass/slow
matching, exact compatibility, physical matching/causal realization,
semiclassical validity/noise, all-sector stability, charge/return, opening,
transport, preparation, separate creation and the full stable-wormhole goal.
