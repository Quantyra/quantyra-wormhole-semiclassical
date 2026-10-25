# Restored canonical quantum commutator and the remaining constrained return

2026-10-24, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: progress. Full stable-wormhole objective ACTIVE/INCOMPLETE.
The closed feasibility audit and original physical norm targets are unchanged.

OUTCOME: the COMPLETE COUPLED GLOBAL RETURN remains UNRESOLVED. The fixed-seed
canonical quantum commutator below satisfies its proposed form sub-budget.
The actual global gluing loss is bounded directly, with all angular/thermal
returns. Weak weighted constraint/gauge returns are also bounded. These do
not identify canonical multiplication with multiplication in the TT variable,
bound the original strong phase rows, or certify a post-elimination column.

## 1. Operator and scope

Use the actual paired mass/log domain, complementary adjoints and composed
square of the current finite seed. Keep the common regulator, the specified
state and temperature, both original equations and all finite-domain returns.
Use the SAME auxiliary W=exp omega of
[the previous attempt](current-weighted-global-return-attempt.md), s=7/5,
omega_z=s theta/(h+1-z), theta=I_(10z-1)(11,11) on [.1,.2], zero/one outside.
W=1 near the seam and is smooth at the regular bolt. No domain is smoothed.

For a symmetric canonical Hessian form q on v=(S,T), define

    c_W[q;v]=Re q(W^-1 v,W v)-q(v,v).                       (1)

This is the symmetrized similarity commutator as a quadratic form. For
q=Gamma_q'' it includes the complete massive determinant, canonical LLL,
Weyl completion and the fixed vacuum subtraction. The latter is a local
zeroth-order form and has zero commutator. No off-shell gradient is declared
zero; we first calculate in the actual canonical chart. Returning (1) to a
weighted TT coordinate is a further operation, displayed in section 6.

The first unsupported step in the full problem is still the strong
constrained interior/phase return. This calculation addresses the earlier
missing COMMON-RESTORED quantum estimate, without inferring that strong
return from a bounded canonical form. The exact cut reference mass lift,
its transition and fixed-temperature bolt matching remain those of
[the coupled calculation](current-coupled-global-return-attempt.md).

## 2. Uniform weight jets on the actual finite geometry

Put d=h+1-z. The retained exact polynomial is

    Q=(1-z)^2 P-tau z^3,
    P=beta(1+2z+3z^2)+u z^2+eps z^4,
    beta=(1-u-eps)/6, eps<0, tau>0.

Coefficient bounds give Q<=1.000001 d^2, |Q'|<=4d,
|Q''|<=9, |Q'''|<=15, |Q4|<=13, |Q5|<=1. In particular the thermal
term in |Q'|/d is <=3tau/h<1e-30, not divided by a vanishing Q.
For a unit-speed orbit geodesic, z'^2<=Q and z''=Q'/2. Repeated
differentiation gives bounds for |z^(j)|/d, j=1,...,6,

    (1.001, 2, 4.51, 16.6, 72, 400).                       (2)

For example z5=Q4 z'^3/2+3Q'''Q'z'/4+(Q'')^2 z'/4,
and z6=Q5 z'^4/2+3Q4 Q'z'^2/2+5Q''Q'''z'^2/4
+3Q'''(Q')^2/8+(Q'')^2Q'/8. The displayed reserves cover each term.
Across the paired seam W is constant, so no high seam jet is requested.

For theta, the Bernstein derivative bounds are T0=1 and
Tk=(21)!/(21-k)! 2^(k-1)10^k, 1<=k<=5. Thus

    |partial_z^j omega| d^j
      <=s sum_(k=0)^(j-1) binom(j-1,k) Tk (1.001)^k (j-1-k)!.

Faà di Bruno with (2), followed by the exponential Bell polynomial,
proves |W^(-a)X^j W^a|<=Aj for a=+/-1, with

    A1<1.401, A2<302, A3<121824,
    A4<4.774e7, A5<1.803e10, A6<6.623e12.                  (3)

The sharper gradient L=1.400000000001 bounds |grad omega|.
These are derivatives along EVERY complete geodesic, not selected modes.

## 3. Actual gluing: use the second weight difference

Retain the two-buffer identities and constants of
[the full curved form](stationary-curved-quantum-tt.md), section 4.
The local bubble-difference row is <=1e24 C, on distance <1/128.
The far row for one mass uses
R_M=9.8304e8 exp(-M/530), with a commutator supported within 1/128
of the observation point, and is <=2 R_M^2 before summation.
Adjoints and grading give the same reversed-row estimate.

Local multiplication by cosh(omega(x)-omega(y))-1 costs at most

    K_local/C <=1e24 (L/128)^2 exp(L/128)/2 <6.048e19.       (4)

For the far part choose, separately for EACH mass, t_M=M/(4L)>1.
The positive Taylor series gives cosh a-1 <=cosh(t_M a)/t_M^2.
Split the exponential equally between the two Dirac kernels in the bubble.
The actual weighted inverses then cost

    ||W^(+/-t_M/2)(Dsharp+M)^-1 W^(-/+t_M/2)||<=8/(7M).

This is the same elementary bounded-perturbation proof on the actual
paired domain; its perturbation has norm M/8. In the endpoint cap row,
the weight ratio costs at most exp(t_M L/(2*128)) per factor. The
two-buffer identity is retained. No weighted extension of the physical
operator to the cap is made. Consequently the absolute far commutator row
for one mass is bounded by

    2(9.8304e8)^2(8/7)^2 (16 L^2/M^2) exp(-a M),
    a=1/265-1/512>0.                                    (5)

The complete count N(E)<=4CE^2 and M1>18 imply, by Stieltjes integration,

    sum d_n exp(-aM_n)/M_n^2
      <=C[4 exp(-18a)+8 E1(18a)]
      <=C[4 exp(-18a)+8(log(1/(18a))+1)]<39.218 C.         (6)

The last elementary bound splits the E1 integral at one. Equations
(4)-(6) give K_far<3.105e21 C and total K_glue<3.166e21 C.
The row-sum subtraction in the original gluing form is a diagonal
multiplier: its commutator is EXACTLY zero. It is retained in the Hessian,
not charged again in (4)-(6). Schur's two-row argument therefore gives

    |c_W[q_glue;T]|<3.166e21 C J0[T].                      (7)

All thermal propagation and finite returns are contained in the actual
resolvent. The sum (6) is infinite and certified analytically. This replaces
the previous >3.56e42 upper-bound expression; it is not a normalization
change. Its old E0 conversion costs <1.647e-7 E0. The .001 gluing allowance
is closed for this canonical commutator, with considerable reserve.

## 4. Restore the ultraviolet cancellation before taking bounds

Use the retained actual/geodesic comparison, with the same finite contacts.
The complete geodesic term is sum d P_M+Q4, where

    P_M=int_0^infinity r dr [-B_M(r)] g(rP),
    g(x)=1-x^2/2+x^4/24-cos x, P=iX,
    Q4=-C F_eff ||X^2T||^2/(90pi^2), |F_eff|<70.

The polynomial and cosine pieces are combined BEFORE the angular limit.
The previous calculation, using the positive-real
[Bessel representation](https://dlmf.nist.gov/10.32.E9), gives

    sum d int r^7 exp(Lr)|B_M(r)|dr <.002671854 C.          (8)

The fixed constant probe for the commutator has sixth remainder bounded
by A6 exp(Lr)r^6/720. Its Q4 point action costs at most A4 C/2.
Including the full 2pi direction integral, its complete point value is

    |p0,geo|/C <=2pi(.002671854)A6/720+A4/2<1.783e8.       (9)

To fix the SECOND contact use a centered normal quadratic q_v with unit
v, equal to (v.x)^2 inside radius R/2, zero after R=1/400 and bounded
by R^2. Along each geodesic from its center it is exactly r^2(v.e)^2
inside R/2. There the sixth remainder is bounded using

    Aq6=A6 R^2/4+6A5 R+30A4.

For r>=R/2 keep all returns. Use sum d|B_M(r)|<=1152 C r^-6 exp(-9r).
The second polynomial commutator acts as zero on this quadratic at its
center; the fourth acts with modulus <=12 L^2. The cosine term is bounded
by R^2(cosh(Lr)-1). Thus the entire second-probe action, including Q4, is
bounded by

    C[2pi(.002671854)Aq6/720
      +2pi(1152L^2+576L^2[log(2/(9R))+1])+6L^2]<93062 C. (10)

No fourth-contact ambiguity is fitted by these probes: it was already
fixed by the common restoration. They determine the remaining zero and
second distributional contacts. They are used to bound the ENTIRE form,
not presented as an isolated contact result.

The actual-minus-geodesic local kernel is <=2e8 C/r^4, with second
contact reserve <=1e13 C, plus its explicit curvature contact. Its
zeroth multiplier commutes. Its symmetric weighted constant probe has
quadratic contribution bounded by 1e15 C L^2; the third Taylor remainder
costs <=1e6 C A3 at radius R. These give <3e15 C. The weighted second
probe has no local second-contact contribution at its center and is
absolutely integrable. Changing the finite-part anchor costs less than
1e17 C in either remaining contact. This last reserve also covers (9)-(10).
The actual constant/tadpole term and seam delta multiplier are retained
and commute; W=1 in their seam neighborhood.

On r<=R the principal separated kernel bound 1e4 C/r^6 times
cosh(delta omega)-1 gains two powers. The gain-two comparison error gains
two more powers. Together, allowing the angular/Jacobian factors as in
the retained normal-cap proof, both fit K=1e6 in the finite-part form
estimate (9) of stationary-curved-quantum-tt.md. The complete flow tail
after R is bounded by the same tilted exponential envelope and costs
less than 1e10 C J0. Long returns are not called local contacts.

For clarity the elementary form estimate uses X=||v'||, Z=||v''||:
split at b<=R, bound the Taylor-subtracted near integral by
KC[bXZ+b^2Z^2/8], and the rest by 2KC log(R/b)X^2.
Set b^2=epsilon/(4K), epsilon=1e-4. Young's inequality and the angular
Hessian/weak Bochner identity give a safe fourth coefficient pi epsilon/2
and gradient coefficient K[1/4+2log(R/b)]<1.268e7. Add the fixed contacts
just bounded. All these inequalities are uniform at a common positive
regulator and persist when the combined regulator is removed. Separating
the divergent second/fourth angular moments would invalidate this step.

Including (7), the full actual MASSIVE canonical commutator therefore obeys

    |c_W[Gamma_mass'';T]|
      <=(pi/20000) C J2[T]+1e19 C J1[T]+4e21 C J0[T].     (11)

The reserves include the finite local comparison and common curvature
contacts, not merely the principal symbol. The argument starts on the
actual strong core and extends by the retained mass/log form completion.
W is constant near the outer seam; the logarithmic boundary flux conditions
are unchanged. No smooth-boundary pseudodifferential theorem is invoked.

## 5. LLL and full Weyl completion in the same canonical chart

The specified canonical LLL form is -q J1[T]/(12pi). Its commutator is
+q int |grad omega|^2 T^2/(12pi), bounded by 21 C J0[T] since q/C<400.
This uses the fixed HH disk realization, including its boundary allocation;
it does not apply a massive gap to the LLL.

Use the EXACT Weyl polynomial in section 8 of stationary-curved-quantum-tt.md,
with its actual boundary representative. Denote its background log radius
by Omega=log(r/B), distinct from the auxiliary omega. A local commutator
vanishes where W is constant, so it suffices to use z>=.1. In normalized
orbit units, the original coefficient bounds imply

    |Omega|<3, |grad Omega|<11, |Rbar|<=2,
    |Delta Omega|<140, |grad Omega|^2<101, B^-2 ell_o^2<1.

For its D,L,J coefficients this gives |D_S|<.25, |D_T|<10,
|D_gradS|<1, |L_T|<72, |L_gradS|<2 and |J|<19.
The weight has |Delta omega|<600 and |Hess omega|<424 from section 2.
For the fourth form int a(Delta u)^2, the exact weighted difference is

    int a[2c u Delta u+c^2u^2-4(v.grad u)^2
                         -4b u(v.grad u)-b^2u^2],
    v=grad omega, b=Delta omega, c=|v|^2.                 (12)

Integrate the u Delta u term once; |grad c|<=2L|Hess omega|.
Only a=constant or a proportional to Omega+1/6 occurs in these fourth
terms, so |a|,|grad a|<1 is sufficient after their displayed coefficients.
For a mixed product Delta u A.grad v0, the symmetric difference is
c u A.grad v0-2(v.grad u)(A.v)v0-bu(A.v)v0.
For (Delta u)v0 it is c u v0. These identities bound all D,L terms.
The gradient squares have zeroth-order differences. First-derivative
bilinears and the full Z polynomial (including its q^2 coefficients)
commute after symmetric pairing. No large magnetic contact is discarded.

Expanding theta=S+T/3 and using 2|ab|<=a^2+b^2 gives, conservatively,

    |c_W[W4'';(S,T)]|
      <=1e8 C [J0[S]+J0[T]+J1[S]+J1[T]].                 (13)

For example every coefficient after (12) is <1e5; the displayed Weyl
polynomial, including the theta expansion, has fewer than 1000 resulting
products. This reserve is intentionally much larger than their sum.
Its boundary term has zero commutator because W=1 in a whole neighborhood
of the boundary; the bolt is regular. No punctured-bolt cutoff is used.

Equations (11)-(13), with ample room in their rounded constants, prove

    |c_W[Gamma_q'';(S,T)]|
      <=(pi/20000) C J2[T]
         +1e19 C(J1[S]+J1[T])+4e21 C(J0[S]+J0[T]).        (14)

On an UNWEIGHTED genuine TT tangent returned to the canonical chart, the
existing physical comparisons give C(J1[S]+J1[T])/E0<2.55e-27 and
C(J0[S]+J0[T])/E0<6.24e-29. Therefore its lower-order loss in (14) is
<2.752e-7 E0, while pi/20000<1/(200pi). This closes the proposed
CANONICAL quantum-form sub-budget, within the .006398999997 E0 margin.
For comparison with the localized fourth margin, use
J2[T]<=2 sum J2[chi_j T] plus the fixed cutoff's lower derivative terms.
The coefficient is then pi/10000<1/(200pi); the cutoff terms fit the
unused lower-order reserves above (their coefficients are below 1e14
before multiplication by pi/20000). Thus localization does not spend the
fourth reserve twice or require a changed physical norm.
It does not prove a strong norm or a constrained weighted inverse.

## 6. Exact constraint return and the remaining distinction

The weight decreases with areal r. The positive-kernel Schur argument for
the actual TT constraint therefore proves

    ||W b[xi]||_(L2(r^2dr))<=2||W xi||,
    ||b[Wxi]-W b[xi]||<=4||Wxi||.                         (15)

Indeed the original dominating kernel is U(r)^-1 U'(t)1_(t<=r),
U=r^2F, and its two Schur integrals are at most two. The extra ratio
W(r)/W(t)<=1 improves this domination. The exact correction is still
delta b=-(r^2F)^-1 int_R^r W' t^2F b dt, with its genuine quartic outer
tail. Its strong trace coefficient is not set to zero by (15).

There is also a finite weighted VALUE bound for the gauge vector eta,
eta'-A eta=b+xi in physical proper radius, A=F_r/(2sqrt F).
If f=partial_rho log W, then -sA<=f<=0. To prove this at the actual seed,

    d(2Q-zQ')-2zQ
      =(1-z)^3[2beta(1+z)-2eps z^4]
        +tau z^3(1+z)+h(2Q-zQ')>=0.                     (16)

The last term is positive because F_r>0. With physical measure
m=a_orb r^2, conjugate psi=Weta to phi=sqrt(m)psi. Its first-order
coefficient is B=3A/2+Rphysical+f>=(3/2-s)A.
Since F_rr>=lambda, F_r^2>=2lambda F and A>=sqrt(lambda/2).
The decaying-end Volterra kernel is bounded by
exp[-(3/2-s)sqrt(lambda/2)(t-rho)]. Young's inequality gives

    ||Weta||<=[3sqrt(2)/(3/2-s)]lambda^-1/2 ||Wxi||
              <42.427 lambda^-1/2 ||Wxi||.                (17)

This retains the actual bolt/outer realization; it is not a freely chosen
Green function. It does not bound A Weta, strong derivatives or point traces.

The exact difference between weighted TT return and canonical weighting
is particularly important. In z let J'= (b+xi)/Q, J(0)=0, and
S=xi/2+QJ/z, T=-(b+3xi+Q'J)/2. Set
delta J=J[Wxi]-WJ and delta b=b[Wxi]-Wb. Then

    delta b'+(Q'/Q-4/z)delta b=-W' b,
    delta J'=delta b/Q-W'J,
    D_S=Q delta J/z, D_T=-(delta b+Q'delta J)/2.            (18)

Their regular-bolt and original outer conditions must be retained. Thus
K(Wxi)=W K(xi)+D, for the actual canonical return K. The form (14) controls
q(W^-1Kxi,WKxi)-q(Kxi,Kxi); the constrained commutator additionally contains
the two cross terms with D_+,D_- and q(D_-,D_+). Equations (15),(17)
DO NOT control those terms in the required quantum graph/phase norm.
Replacing K(Wxi) by WKxi would silently omit them.

This is now the precise first missing estimate for this attempted
constrained return: a strong bound for (18) and its original finite-part
rows, followed by the signed gravitational/scalar/interior resolvent.
The current source/fast-input constants do not provide that estimate.

## 7. Research decision and unchanged global requirements

The global Schur columns L_xy-L_xI L_II^-1 L_Iy remain unevaluated for
mass, scalar quartic, correlated slow pair and fast remainder. Neither
kappa_m<=5e-5 nor any actual mass coefficient is established. At that
conditional target, the remaining affine/graph budgets are 5.34235e-12
and .0012069 (rounded), before independent energy and all other entries.
The difference-energy debit about .0195239 remains. None is assumed to fit.
Both original graph ceilings still require actual |m|<2.205539e-9; the
4.345e-8 allowance is not the measured mass. No boundary basis or weight
has changed this requirement.

The conditional fast screen remains .135973<.1599164 only in its proved
two-fast-input scope. Independent phase variations and neighborhood
comparisons are unproved. The original requirements remain
Ecenter+1e-8 alpha<=1e-11,
R_X alpha*.028913+Z0_other<=.25,
R_X alpha^2*1.0706e10+Z2_other<=5e7;
Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1.
No full coupled admission or contraction requirement is closed by (14).

RESEARCH DECISION: continue the signed weighted method, retaining the
construction. The complete canonical quantum commutator is quantitatively
attainable at the seed; its former divergent/contact and enormous gluing
estimates are replaced by (14). The full feasibility question remains open.
Do not spend another turn sharpening these already sufficient constants.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): bound the complete
constraint-return correction (18) in the original strong quantum graph
and phase norms, including its cross terms in the bordered interior
resolvent. Test the resulting actual coupled columns against the unchanged
source/feedback and original physical-norm requirements. A weak L2 estimate
or another canonical quantum calculation cannot substitute for this return.

Executed: full seed canonical quantum commutator (all modes, common contacts,
finite returns and Weyl/LLL allocations); exact constraint return equations
and their weak weighted bounds. Unexecuted: strong constrained return,
actual coupled solve and coefficients, neighborhood transfer and nonlinear
solution. No nonexistence or stability conclusion is drawn.

Remaining to-do list: strong constrained global/phase return; actual correction and original physical norm comparisons; nonlinear admission/contraction and exact compatibility; physical matching/causal realization; semiclassical validity/noise; all-sector stability; charge/return, opening, transport, preparation, separate creation and the full stable-wormhole objective.
