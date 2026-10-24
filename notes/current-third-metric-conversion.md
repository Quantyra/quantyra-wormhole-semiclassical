# Third metric variation: common-regulator conversion and the remaining nonlocal test

2026-10-24, E003/S169. Ordinary argument before acceptance.
Previous turn: progress. Objective ACTIVE and INCOMPLETE; audit CLOSED.

The selected calculation is the metric variation of the subtracted
quantum Hessian in the actual coupled graph. Before estimating a radial
three-vertex loop, its conversion from the geometric action must also
be known to third order. The preceding conversion theorem only certified
two jets. This note closes the third seed jet, including its regulator
and boundary terms, and bounds its actual outer fast graph contribution.
It does not certify the complete three-vertex strong response.

## 1. Actual domain, mass representation and the common heat trace

Keep the original paired first-order D, its complementary adjoint domain,
grading, regular disk bolt, fixed thermal spin structure and entire
angular tower. All fields and variations are identical on both copies.
For b in the admitted cubic/mass/log phase put F=exp(b). On the fixed
seed Hilbert space the mass representation is D+M F, with Dom D unchanged.
Multiplication by F and F^-1 preserves the first-order domain. Its square
is defined by composition, not by a separately selected boundary problem.
For real b, accretivity and the adjoint estimate give inverse norm
<=exp(||b||infinity)/M. The same argument gives a positive massive gap.

The exact positive-cutoff identity from current-variable-contact-conversion.md
is retained for EVERY b and r in [0,1]:

 A_r=exp(-rb/2)(D+M exp(b))exp(-rb/2), P_r=A_r^* A_r,
 J_e(b)=Gamma_g,e(b)-Gamma_m,e(b)
       =sum_n d_n int_0^1 Tr_bar[b exp(-e P_r(b))] dr.   (1)

The grading supplies the partner-square equality. All squares have their
composed domains. Both finite geometry and thermal spectrum occur in (1).
This is not an identity between separately renormalized mode determinants.

There is a useful direct way to control the additional mass insertions.
In the actual proper collar coordinate let f=1/a be the fixed inverse
time radius. At frequency p the mass square, after radial half-density
conjugation, is

 -partial_y^2+p^2 f^2+M^2 F^2+p f' Sigma1+M F' Sigma2, (2)

up to fixed signs of the two matrices. The matrices are Hermitian,
anticommute, and square to one. Equation (2) follows by squaring the
first-order operator, including [MF,D]. It retains the mass-gradient
vertex. Omitting that term would invalidate the conversion.

For a variance-two bridge xi_s, h=sqrt(t), set A_f=int f(y+h xi)^2,
A_F=int F(y+h xi)^2. The exact matrix Feynman--Kac expression is

 (4pi t)^(-1/2) E{exp[-t(p^2 A_f+M^2 A_F)]
   tr Texp[-t int(p f'(y+h xi)Sigma1
                            +M F'(y+h xi)Sigma2) ds]}. (3)

T denotes ordered multiplication, not a scalar exponential of its
integral. Scalar positive potentials commute with every matrix and have
already been factored out. The Trotter product formula proves (3) for
bounded local coefficients; the same integrable matrix norm majorant
passes to the original W3,infinity coefficients. Auxiliary local
extensions retain their entire separated gluing defect. They do not
replace the paired seam or the actual bolt.

Here is the uniform ultraviolet argument needed below. Expand the ordered
matrix factor through four insertions. A word has zero trace unless each
of Sigma1 and Sigma2 occurs an even number of times. The retained terms
therefore contain even powers p^(2a) M^(2b), a+b<=2. The FULL Gaussian
frequency integral and angular sums become Gaussian moments and
(-partial_s)^b Theta(s), s=t A_F. Use the same angular heat expansion
Theta(s)=2C/s+H0+H1 s+O(s^2), and its differentiated versions through
order two. These expansions are uniform for A_F in a bounded positive
sector: the finite Euler--Maclaurin remainder for
2(q+2n)exp[-s n(n+q)/C] is bounded there by its Gaussian derivative
integral. Keeping eight derivatives suffices also after two s derivatives.
No angular level is omitted.

For a term with l matrix insertions, expand spatial translations through
degree 5-l in h. Its integral remainder has total spatial derivative
count six, including the l derivatives already on f,F. In a localized
smeared trace, integrate any excess above three on one factor onto the
other factors and the localization. There can be at most one such factor
at a time. This takes at most three integrations and preserves matrix
order and total derivative count. Denominator derivatives of A_f,A_F
produce further translated factors but no inverse power of h. The
result is controlled by W3,infinity norms and finite bridge moments.
Odd h terms vanish under bridge reflection. The ordered tail starts
at six insertions: its absolute norm is bounded by the corresponding
scalar exponential tail. Gaussian p moments and angular counting turn
its t^6(|p|+M)^6 factor into O(t^3) before the t^-2 trace prefactor.
Completing the Gaussian squares controls the exponential remainder.

Consequently the COMPLETE smeared trace has remainder O(t) after its
t^-2, t^-1 and constant terms, uniformly on bounded coefficient families
with positive ellipticity margins. This is a weak, integrated remainder.
It is not an H3 output estimate. All constants depend on the stated
coefficient, localization and finite-seed atlas bounds; no useful
numerical value for this generic constant is claimed.

For the actual antiperiodic thermal sum, Poisson summation of the retained
Gaussian polynomial moments includes every image. Its nonzero images
are exponentially small at small t, with their polynomial prefactors.
For the ordered tail use the entire discrete Gaussian-moment bound
directly; its extra lattice endpoint allowance is a higher power of t.
Thus no zero-temperature physical substitution is made. Regular bolt
charts use the fixed smooth seed spin connection and its heat transport;
the variable mass factor has the same ordered expansion. Local chart
defects are separated and exponentially bounded before summing modes.

The argument extends to a small complex W3,infinity tube in F with
Re F^2 bounded below. Ordered series, Gaussian majorants and the massive
gap are uniform on a smaller tube. Cauchy thus permits three MASS-variable
variations of this subtracted integrated trace and its large-time tail.
The renormalized mass action has a bounded weak third variation on
bounded admitted families. This qualitative statement supplies no
numerical strong norm or Z2 certificate.

For (1) also use f_r and F_r=exp((1-r)b) in the r-proper coordinate.
The coordinate return and the smeared insertions have the bounded first
two parameter derivatives into W3,infinity already proved on the
weighted C5 phase. At b=0 the THIRD variation of (1) requires only TWO
variations of its heat factor: the outside smearing b vanishes at zero.
The preceding O(e) endpoint remainder therefore still tends to zero,
uniformly for three bounded admitted directions. This proves common
regulator removal for the THIRD SEED conversion jet. It does not use
five dilation jets or claim a uniform strong geometric neighborhood limit.

## 2. Compute the finite conversion without a mode or boundary surrogate

For a geometric interpretation, angular reduction with the original
spinor density is an isometry. The metrics underlying A_r are the
four-dimensional warped products

 g_r=exp(2rb)g2+exp(-2(1-r)b)gS
                         =exp(2rb)[g2+exp(-2b)gS].    (4)

The usual warped spin connection contributes the gradient of the sphere
radius; conjugation by that radius cancels it and leaves precisely the
mass M exp((1-r)b). Scalar conjugations preserve both paired traces.
Thus the finite endpoint coefficient in (1) is the four-dimensional
Weyl transgression with the zero angular coefficient subtracted.
This uses the FULL heat trace. It does not assume that a separately
renormalized two-dimensional determinant has the four-dimensional anomaly.

The smooth local coefficient is the original one in
spherical-anomaly-action.md and current-finite-local-action-assembly.md;
the general local formulas are reviewed in
[Vassilevich, sections 4.1 and 7.1](https://arxiv.org/html/hep-th/0306138v3).
The weak remainder argument above, in addition to that smooth algebra,
supplies the continuation to this actual C2,1 paired seam.

At the seam b=b'=b''=0. In the smeared a4 boundary terms the first three
smearing jets therefore vanish; the balanced projector trace cancels the
remaining third-normal-smearing term. Curvature-divergence distributions
are retained before pairing. The ordinary original boundary completion
is unchanged. At the bolt there is no new radial boundary. The zero
angular subtraction here concerns the ENDPOINT HEAT coefficient only;
possible global LLL zero-mode measures cancel in that subtraction and
are not presumed absent from the separate physical LLL problem.

Let J be normalized by -2G/B^2, fq=GC/(pi B^2), R=-Q'',
d=Delta b=(Qb')', g=Q(b')^2, u=1/C, k=1/15-q^2/3. Its finite local
representative is

 J(b)/fq=int dz {
 e^(-2b)[b R^2/120-b R d/30+b d^2/30+R d/60
                         +d^2/60-R g/60-d g/180+7g^2/360]
 -u b R/36+u b d/18+u g/36+k u^2 b e^(2b)/2
                                      -q u(g+R b)/12 }. (5)

Equivalently (5) is minus the retained normalized W4 at phi=-b, shape=b,
minus the zero-angular endpoint cocycle fq q/(12C)int(g+Rb).
The latter is an ultraviolet subtraction identity, not an extra LLL
physical term. Every coupling and mu0 remains fixed. The endpoint trace
has no logarithmic divergence after the common power terms; its finite
part agrees with the original common Mellin/zeta difference. No finite
fourth contact is added to the already matched normal operator.

Equation (5) agrees with the preceding COMPLETE first and second jets.
In particular its first Euler row is

 fq[R^2/120+Delta R/60+H0 R/(12C)-H1/C],
 H0=-q-1/3, H1=(q^2/6-1/30)/C,

and its second row is exactly equation (10) of
current-variable-contact-conversion.md. Checking these identities is a
normalization and finite-contact test, not the proof of the third jet.
The proof is (1)-(4) with the common remainder. Formula (5) supplies a
local differentiable continuation; the quantitative result below concerns
its THIRD JET AT THE SEED, not a uniform nonlinear strong bound.

## 3. The cancellation in the actual third jet

The homogeneous cubic density of (5) is exactly

 L3=R b^2 Delta b/10+R b Q(b')^2/30+b^3 R^2/60
                            -(Delta b)Q(b')^2/180+k b^3/C^2. (6)

The b(Delta b)^2 terms CANCEL between the exponential density and its
vertices. This is an identity in the specified conversion; there is no
choice of a source or counterterm. In the stationary problem the last
derivative cubic also loses its second derivative after integration by
parts. Put

 A=-RQ/6, B=-QQ'/540,
 V=R^2/60+Delta R/30+k/C^2.

Then L3=A b(b')^2+B(b')^3+V b^3+partial_z Phi, where

 Phi=RQ b^2 b'/10-QR' b^3/30-Q^2(b')^3/540.           (7)

All three boundary fluxes vanish on the admitted mass/log directions;
at z=0 the slowest is O(z^6), up to fixed logs, and at the regular bolt
Q=0 with bounded regular radial jets. Keep the seam distribution in
Delta R before its b^3 pairing. Its trace vanishes. This is an identity
on the original paired domain, extended by its existing completion.

The derivative of the conversion Hessian, as a strong bilinear local
row at b=0, is fq times

 N(u,v)=-2A[u'v'+uv''+vu'']-2A'[uv'+vu']
       -6B[u'v''+v'u'']-6B'u'v'+6Vuv.                (8)

This contains every third-order conversion contact. The third action
variation is int w fq N(u,v), symmetric in all three arguments after
the retained boundary terms. It is not the nonlocal three-vertex loop.

## 4. Execute its full outer fast-graph norm

Use sigma=4.005, a=.01, x=log z and the ACTUAL full B_sigma inverse:
u=z^sigma B_sigma^-1 g1, v=z^sigma B_sigma^-1 g2.
Apply the original shape row chi z^-sigma(D-4)z^4 to fq N(u,v).
The scalar row of this canonical conversion is zero, so the original
coupled return introduces no new scalar block here.

To make every coefficient explicit, evaluate (8) on u=z^s,v=z^t:

 F(z;s,t)=-2A[st+s(s-1)+t(t-1)]/z^2-2A'(s+t)/z
             -6B st(s+t-2)/z^3-6B' st/z^2+6V,
 Rrow=z^4[(s+t)F+z partial_z F]=sum_p z^p c_p(s,t).   (9)

Substitute the COMPLETE correlated finite-seed Q polynomial. The powers
are p=2,...,12. Expand c_p(sigma+xi,sigma+eta)=sum c_pij xi^i eta^j.
Every i,j<=3 and i+j<=4. After three output derivatives each input
therefore needs at most six derivatives, and the lesser order is <=3.
No eighth derivative in L-infinity, no sampled frequency inverse and no
finite selection of physical modes enters the estimate.

The retained H6 graph constant is 1801. For j=0,...,3 the exact Sobolev
evaluation constants in H6 are c_j=sqrt((63,7,3,3)_j/512). Define

 P_ij=1801^2 sum_(l=0)^3 sqrt(binomial(3,l))
                  sum_(m=0)^l binomial(l,m)c_min(i+m,j+l-m),
 M(t)=t^3+12t^2+51t+100.

The integer Sobolev product rule bounds each product in H3 by P_ij.
The existing output cutoff has derivative guards (1,3,10,60), so the
coefficient multiplication bound is a^(sigma+p) M(sigma+p). Thus

 ||chi z^-sigma(D-4)z^4 fq N(u,v)||_H3
 <= Cconv,3 ||g1||_H3 ||g2||_H3,
 Cconv,3=fq sum_p a^(sigma+p)M(sigma+p)
                                         sum_ij |c_pij|P_ij. (10)

Interval arithmetic at the actual parameters gives

 Cconv,3=9.84082368339454...e-33 <9.841e-33.           (11)

This closes the SEED OUTER FAST conversion contribution to the selected
bilinear graph calculation. It is not the global or neighborhood Z2:
mass/slow coefficients, regular-chart transfer, the full inverse and
their relation to the original X norm remain separate obligations.
The already certified linear conversion bound <8.357937945e-5 is
unchanged; (11) is a different derivative and is not added to it.

## 5. What is still missing and the research decision

The mass action's third weak variation has the cyclic resolvent-word
structure

 -Tr(R V_u R V_v R V_w+R V_u R V_w R V_v)
 +Tr(R V_uv R V_w+R V_uw R V_v+R V_vw R V_u)
 -Tr(R V_uvw),  V_I=M exp(b) product_(i in I) i.      (12)

At a positive heat cutoff these words must be realized by their
differentiated heat/Duhamel integrals, with the regulator variation
contacts. Simply inserting heat factors into bare resolvent words is
not an exact cutoff identity. Equation (12) specifies the necessary
vertex inventory, not an already evaluated strong kernel. The last
four terms are essential coincidence contacts. All restorations remain
fixed by the same heat coefficients and (1)-(5). Individual traces or
separately summed vertices are not ultraviolet-finite expressions.
The weak argument in section 1 proves the combined functional's
integrability, not a strong estimate on its output derivatives.

UNKNOWN: the actual H3 norm of the moment-subtracted three-vertex response
and its uniform neighborhood control. No finite complete strong upper
bound has been established, so it is not justified to call that bound
excessive. The previous five-jet W3,infinity method obstruction remains
valid; (7)-(11) show how the actual conversion avoids that artificial
loss. They establish no obstruction to the physical operator or solution.

RESEARCH DECISION: retain the construction and the direct Sobolev method.
The conversion is now exact through the required third seed variation
and its outer fast contact contribution is negligible on the present
scale. Do not reoptimize this contact or the cold/scalar components.
Concentrate the estimate on the COMPLETE subtracted mass three-vertex
combination, preserving its cancellations before taking norms.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): bound (12), after its
matched common principal/local subtraction and (11), as a bilinear map
from the coupled fast graph to H3, uniformly in the admitted neighborhood.
Keep the highest derivatives in L2 and every angular/thermal tail,
reflected term, diagonal contact and actual finite-domain return. Test
whether a finite numerical bound can fit the remaining total Z2<=5e7
allocation after the original X/graph comparison; that comparison is
not presumed. This is the remaining decisive strong calculation.

No full signed projection or controlled solution was computed here.
The preserved sufficient E_rest<2.123e-42<1.287e-40 enclosure has zero
as an enclosure center, not a measured sign. Seed M_A<=3.18e27,
M_V<=2.66e30 and eta_nonlocal<.1599164 remain unknown requirements.
At rho=1e-8, Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1 are still
unclosed; Y0<=5e-9, Z0<=.25 and Z2<=5e7 are unproved allocations.
Exact compatibility, physical matching, validity/noise and all-sector
stability remain unresolved. The closed audit and full objective persist.

Remaining to-do list: complete nonlocal metric response and seed quantum
norms, actual global correction, mass/slow matching, full nonlinear
admission and contraction, compatibility, physical/causal matching,
validity/noise, all-sector stability, charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
