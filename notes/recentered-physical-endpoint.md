# The reference inverse reaches the physical stationary endpoint

2026-10-22, E003/S169. Ordinary argument before focused verification.
The full stable semiclassical wormhole objective remains incomplete.

The actual recentered reference inverse, applied to its complete static
Dirac residual, has the physical mass-plus-fast tail, including the five
areal jets required by the stationary boundary theorem. The proof retains
the actual remote kernel moment. A differentiated row annihilates that
moment for purposes of regularity; the original equation fixes it and is
never replaced by the differentiated equation. This is a reference linear
regularity result, not a numerical admission bound or a nonlinear solution.

## 1. Setting and the existing inverse

Keep the c_R=10^22 action, new RN reference, fixed magnetic flux, mu0,
positive prepared period and actual paired Dirac domain of
[the reference inverse](recentered-reference-inverse.md). The previous
[stationary constraint result](recentered-stationary-constraints.md) and
[mass-boundary result](recentered-mass-boundary.md) retain their scopes.
Use physical-volume rows, or equivalently the fixed-reference density
rows divided by its volume; keep their off-shell contacts throughout.

Let H u=f be the EXISTING bulk inverse equation. The disk tangent is
u=(S,T)=(delta Phi,delta V). After subtracting its correlated mass cubic,
the previous result gives u in the mixed weight-7/2 graph space. We do
not solve a new differentiated problem and assume that it is equivalent.
All normal rows below have been divided by the fixed lambda^2; restore
that factor in dimensional multiplier bounds. Put D=z partial_z.

The analytic one-sided RN compactification has the actual smooth mixed
boundary problem. Its reflected double need not be C-infinity: the mass
cubic is treated on the half-space. The fixed compact square is kernel
free, including the LLL as already proved. Thus its separated Green
kernel is smooth up to each boundary face away from the diagonal.
All constants below are finite at this fixed positive temperature.
No uniform zero-temperature bound or new numerical gap is asserted.

## 2. The actual remote moment and its exact removal

Work directly in disk conformal/shape variables. The conformal and mixed
rows are local by the full trace Ward identity; only the shape/shape row
has a nonlocal kernel. In its coordinate-density kernel subtract the
complete universal normal kernel K0, including the paired image. Local
contacts are kept separately. Call the nonlocal difference K_R(z,y).

The retained two-power boundary parametrix gives, for 0<=z<y/2,

    |partial_z^j K_R(z,y)|<=C_j y^(-3-j)L(y)^N,          (1)
    L(y)=1-log y, 0<y<=a<1.

Here j is any of the finitely many orders used below (six suffice).
This strengthens the previous conormal bound to the one-sided face
trace needed now. The justification is specific to this reference:
the Dirac coefficients and boundary operators are smooth to the face;
the direct/reflected resolvent parametrix can be differentiated there.
Its first coefficient error is quadratic, and its zeroth-order connection
error is linear, so the first Green correction gains two powers. Two
stress vertices and tangential integration give degree -3 instead of -5.
Each normal derivative consumes one power. Finite logarithmic terms
give L^N, with N increased if necessary. A separated global Green
remainder is smooth at this face and obeys a weaker instance of (1).

The differentiated differential-boundary resolvent estimates in
[Gaarde--Grubb, Theorem 2.4](https://arxiv.org/pdf/math/0703878) justify
the face differentiations. Proposition 2.9 also warns that the reflected
logarithm has a nonstandard singular term. We extract that term as K0;
we do not call the full logarithmic operator smoothing. These are local
parametrix inputs, not a ready-made determinant-Hessian theorem. The
two-vertex heat formula, its contact and the full common angular heat
sum remain those in the retained shape-response proof.

On comparable annuli retain the proved H4log-to-L2 finite-part estimate
with gain z^2, and its finite differentiated versions. Never majorize
the diagonal finite part by an absolute kernel. On the analytic reference,
commuting three conormal derivatives through the original elliptic
system upgrades a source with three conormal derivatives to u with
seven (and the shape logarithmic multiplier). Each commutator has the
same order and already controlled lower derivatives; local disk
regularity handles the bolt. This is a finite bootstrap at the fixed
reference, not an open-ball C2 quantum-source theorem.

For the fast remainder v in weight 7/2, the REMAINDER moment is the
ordinary convergent functional

    M_R[v]=integral_0^a K_R(0,y)v(y)dy
                         + integral_compact K_R(0,y)v(y)dy.         (2)

Include the separated cutoff and fixed-density factors in this definition.
If several matrix components occur, use their finite sum. The universal
K0 moment is not included in (2); it is generally divergent at this
weight and is handled exactly by the normal operator. For v=y^(7/2)w,
Cauchy--Schwarz gives the explicit collar factor

    |M_R[v]| <= C_0 sqrt(I_N(3,a)) ||w||_L2(dy/y)
                                      +C_sep ||v||_compact,
    I_N(b,a)=integral_0^a y^(b-1)L(y)^(2N)dy
      =a^b sum_(j=0)^(2N) binom(2N,j)L(a)^(2N-j) j!/b^(j+1).      (3)

The first normal Taylor remainder has coefficient
C_1 sqrt(I_N(1,a)) times z. For y<2z, the annular bound and the
discarded part of (2) are O(z^(3/2)L^N); for y>2z use (1) and Taylor's
formula. Repeating D derivatives gives the same powers: when the
higher derivative integral stops converging at zero, its lower endpoint
supplies z^(3/2)L^N. Consequently

    R_nonlocal v=z^4 M_R[v]+O_b^3(z^5 L^N).             (4)

The O_b statement here may be read in annular Sobolev norms sufficient
for the three output derivatives; it is not a pointwise diagonal-kernel
claim. Local remainder terms are O_b^3(z^(11/2)L^N). Smooth density
factors differing from their limit by O(z^2) do not change the leading
moment. The compact part contributes its actual boundary trace to (2).

The mass cubic is treated separately. Its remainder moment is also
integrable (y^-3 times y^3); its first Taylor remainder can be
O(z log^N z), still beyond weight four after the physical z^4 factor.
Local corrections to the mass column start at degree five. No mass
moment is assigned zero or promoted to a free parameter.

Thus the exact operation D-4 annihilates the leading moment in (4).
It also removes the universal remote pole: K0(0,y)=-24d0 y^-5,
d0=1/(480pi^2), and

    (D-4)[z^4 integral K0(z,y)v(y)dy]
                              =z^5 integral partial_z K0(z,y)v(y)dy

away from its retained local contacts. The improved remote kernel is
integrable at the upper Mellin endpoint 4, in fact on a strip extending
at least to 5. This derives the endpoint continuation from the actual
kernel operation. It is not a continuation of an unjustified formula
past the original strip. Contacts are differentiated as distributions.

## 3. Remove the scalar mixing without a singular York inverse

Let a(s)=(s+3)/(3(s+1)), q(s)=s(s-3), and use the previous symbols

    J(s)=3[C_v-D_R lambda q(s)][q(s)-4]/(4pi G lambda),
    N(s)=d0 q(s)[q(s)+2][F_M(s)+kappa_M]-A_T q(s).

D_R denotes the positive curvature coefficient, distinct from D=z d_z.
The disk normal matrix is C(3-s)^T diag(J,N)C(s), with
C(s)=[[1,a(s)],[0,2]]. Its two rows are

    f1=J theta, f2=a(3-D)J theta+4N T,
    theta=S+a(D)T.                                      (5)

The input operator a(D) is harmless here: (D+1)^-1 T equals
z^-1 integral_0^z T(y)dy; its growing homogeneous solution is excluded.
The output a(3-D) has a pole at 4 and must NOT be inverted there.
Instead apply the local row combination

    g2=(D-4)f2-(D-6)f1/3,
    g2=4 M(D)T, M(s)=(s-4)N(s).                        (6)

Indeed (s-4)a(3-s)=(s-6)/3 exactly. This keeps the scalar mixing and
avoids introducing a singular output York transformation. The same
combination applied to the actual equation has the higher-order
remainders of section 2; its only nonlocal leading moment cancels.
The conformal row is local, so it introduces no other remote moment.

## 4. A proved endpoint multiplier, not a sampled contour

Write L_T=A_T/d0>30*2000^10 and 39/5<kappa_M<54/5. On s=4+it,
t>0, the exact digamma reflection and recurrence identities give

    Re F_M=Re psi(3+it),
    -t Im F_M=pi t/(2sinh(pi t))+1/2
                              +t^2/(1+t^2)+t^2/(4+t^2).           (7)

Thus 1/2<-t Im F_M<3. Re H>0 for H=F_M+kappa_M, and Re H<13 for
0<t<=sqrt(6); it is <12 for t<=1/2. For these upper bounds use
Re psi(3+it)<=psi(3)+t^2 sum_(n>=0)(n+3)^-3 and the sum <5/54.

For t<=1/2, q(q+2) has real part 24-35t^2+t^4>=245/16 and imaginary
part 50t-10t^3. Therefore

    Re M >= (245/32)d0+(5A_T-600d0)t^2 >6d0.           (8)

For 1/2<=t<=sqrt(6), Re[(q+2)H]<93 and |t q|>=2, giving
|M|>=2(A_T-93d0)>6d0. For t>=sqrt(6), Re(q+2)<=0 and (7) gives
Re[(q+2)H]<15, again a stronger lower bound. At t=0 use
M(4)=24d0. Complex conjugation treats t<0. Hence on the WHOLE line

    |M(4+it)|>=6d0,
    ||(4 lambda^2 M)^-1||_L2->L2 <=20pi^2/lambda^2.    (9)

This numerical bound is for the transformed NORMAL shape multiplier.
It is not a global admission-norm inverse. At high frequency its symbol
has order five times log, so the usual finite differentiated multiplier
bounds give the corresponding derivative gain, with finite constants.

There are no zeros to cross from 7/2 to 4. To check the missing strip,
put e=s-4, Re e in [-1/2,0]. For |Im s|<=2,

    F_M=1/e+H_reg-kappa_M, |H_reg|<17, |q+2|<11,
    Re[(q+2)/e]=Re[6/e+5+e]<=5.

The cotangent partial fractions bound its pole-subtracted part by 3;
the digamma recurrence/series bound |psi(s-1)| by 3. Thus
Re[(q+2)(F_M+kappa_M)]<192<L_T away from the pole. For
2<=|t|<=sqrt(6) the bound is <48; for |t|>=sqrt(6) it is <20,
using |t Im F_M|<4, Re H>5 and the same elementary series. These
are continuum inequalities over the strip. The only pole at 4 has
already been removed in M, with its nonzero value in (9).

## 5. The actual source and the physical five-jet conclusion

On this analytic RN reference, the full residual rows have the finite
one-sided expansion

    f_i=z^4 P_i(log z)+O_b^3(z^(4+nu)), 0<nu<1,        (10)

with finite polynomials P_i. For the trace row P1 is constant; in the
shape row P2 has degree at most one. Coefficients include the actual
global state; no finite coefficient in P2 is set to zero.

Here is why (10) holds for this source, beyond a generic conormal bound.
The scalar row is the exact local polynomial trace already evaluated.
For the canonical shape source use the full four-dimensional smooth
boundary resolvent/heat parametrix, summed over the sphere before
renormalizing. Its finite expansion through the constant term gives
powers and a single logarithm in normal distance; its one-order remainder
and three conormal derivatives are O(z log^N z). The smooth separated
large-time part has a boundary Taylor expansion because the compact
square is kernel free. The known 1/z mass-cubic coefficient cancels its
Weyl partner, and the constant physical AdS vacuum cancels the exact
lambda recentering. These are the retained mass-layer cancellations,
valid for the new analytic RN coefficients without invoking an Einstein
equation. The remaining physical factor is z^4. The LLL regular HH
allocation and local R^2/Maxwell terms have ordinary finite expansions.
The fixed scale and all restored heat coefficients remain included.
This proves (10) at the reference, with finite but unquantified remainder
and global constant. It does not prove a uniform expansion for arbitrary
nonlinear C5 metric perturbations.

Subtract the known mass cubic from u and commute the three derivatives
as in section 2. The actual-minus-normal first row is beyond weight four;
the second row has (4). The combination (6) removes its moment, so its
right side is z^4 times the polynomial

    Q(log z)=P2'(log z)+(2P1(log z)-P1'(log z))/3        (11)

plus a remainder in weight four with two Sobolev derivatives (indeed
with some positive excess weight). Since M(4)!=0, the finite triangular
polynomial system 4M(4+partial_t)P_T(t)=Q(t) has a unique polynomial
solution. At this source Q is constant and

    P_T=Q/(96d0)=5pi^2 Q                              (12)

in the lambda^2-divided row convention. Subtract a cut-off z^4 P_T.
The normal equation for the remainder has weight-four H2 data. The
Mellin inverse on (9) and the zero-free strip give T_remainder in
z^4 H_b^(7,log). Cutoff commutators include their actual separated
tails; after (6) they are beyond weight four. Equivalence with the
old weight-7/2 solution follows by contour movement through this
zero-free strip, not by choosing a second inverse solution.

The first row is local regular-singular. Its simple root 4 supplies
z^4 P_theta(log z) of degree at most deg P1+1; its independent z^4
coefficient is fixed by the global solution. Its other roots are -1
and Delta_+/- outside the strip. Local contour/Frobenius reduction
puts the remainder beyond weight four. Finally S=theta-a(D)T;
a(D) is regular at 4 and preserves these finite logarithmic bounds.
One-dimensional Sobolev embedding H_b^7 into six bounded conormal
derivatives then gives

    u=m_c e_m+O_b^6(z^4(1+|log z|)^p), p finite.       (13)

For this source p=1 suffices. The smooth disk-to-areal formula involves
one derivative of Phi. Six disk jets therefore give the FIVE areal
jets of R_F and R_alpha required by the physical fast phase. The
correlated mass cubic is unchanged. Bolt regularity and fixed period
are inherited from the original inverse; no horizon condition is added.

The discarded moment has NOT become an additional equation or datum.
For a solution of the differentiated system alone, the original shape
residual could still equal c z^4. We began with H u=f, so that constant
is exactly zero. Equations (2)-(4) describe its actual global origin.
This also explains why the present regularity proof cannot be used as
an independent construction method without retaining the original row.

## 6. Quantitative consequence and the next calculation

The reference linear response is now an admissible tangent of the
stationary variational phase, and the zero-current and static recovery
identities apply. This closes the previously missing source-specific
physical-tail implication. It does not quantify its norm or prove an
open nonlinear ball. The moment coefficients C_j,C_sep, the source
coefficients and the global mass coefficient still need useful bounds.
For an allowed moment budget B_M, (3) requires

    C_0 sqrt(I_N(3,a)) ||z^-7/2 v||+C_sep||v||_compact <=B_M.       (14)

These constants are unknown, not upper bounds already demonstrated to
fail. The exact kernel moment is removed only from the REGULARITY
equation; dropping it from the original solver would be incorrect.

There is a serious quantitative trap in separately bounding the terms
of (13). The retained meromorphic zero is at 4+e_* with
e_* about 1.865097e-34. A very slow compensating tail can cancel a
large constant z^4 coefficient over an enormous log-radius interval.
The following exactly solvable normal comparison shows the right bound
to seek, without claiming it is the actual quantum response:

    -A w(t)-R integral_0^t w(s)ds=b+a t,
    A=4A_T, R=24d0, eta=R/A,
    w(t)=-(b/A)exp(-eta t)-(a/R)[1-exp(-eta t)].        (15)

Direct differentiation and w(0)=-b/A prove (15). Since
1-exp(-eta t)<=eta t,

    sup_(t>=0) |w(t)|/(1+t) <=max(|b|/A,|a|/A).        (16)

The separate constant coefficient is |a|/R; losing its cancellation
costs A/R=1/eta about 5.36e33. This is an actual cancellation in a
specified operator comparison, not a normalization improvement. It
does not assign that bound to the full kernel. A direct estimate of
the combined actual correction, retaining its remote moment and
remainder, is therefore the next quantitative step.

Continue the construction with this method. At rho=1e-8 the complete
nonlinear requirements remain Y0<=5e-9, Z0<=.25, Z2<=5e7 plus every
positive admission margin and exact compatibility. The reference energy
inverse 1 and normal endpoint bound (9) live in different stated norms;
neither is the missing complete admission inverse. The old huge uniform
allowances remain parked and are not nonexistence evidence. Matching
about the new full reference, semiclassical validity and all-sector
retarded stability remain unproved. S169.2 stays unchecked.

Freeze this argument, checker and review before execution. Check the
row combination, the whole-contour inequalities' algebra, moment factors
and source-polynomial inversion. Independent separated-kernel integrals
must retain a nonzero remote moment and test its exact annihilation;
they are fixtures, not samples of the actual Dirac kernel. Verify the
Volterra comparison against its integral equation, including a regime
where separate terms are much larger than their sum. Use 100/150 digits
for actual benchmark scales and retain every run. These checks do not
certify the continuum parametrix, source expansion or a nonlinear ball.
No Lean invocation is selected.

Remaining to-do list: quantify the combined actual linear correction,
mass and kernel moments in the physical admission norm; finite-state
transfer and nonlinear existence; compatibility, matching, validity,
all-sector stability, charge/return, opening, transport, preparation,
separate creation and full-result audit.
