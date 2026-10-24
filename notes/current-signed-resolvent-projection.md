# Actual paired resolvent projection: subtract before numerical evaluation

2026-10-23, E003/S169. Ordinary argument and protocol before execution.
Previous turn: progress. The full stable-wormhole goal remains active and
incomplete; the closed feasibility audit is preserved.

## 1. The actual functional, and the first implementation gap

The previous checkpoint provides a selected common-regulator response,
complete local contacts and an angular-tail certificate. It does not
provide an executable value of Q_(P64) or values of the actual h1. A finite
matrix or sampled frequencies would not fill that gap. The direct radial
resolvent below represents the SAME paired operator and its derivative.
Its approximation error must be proved independently; the old Q_(P64)
error cannot be transferred to a different approximation by name.

Keep the fixed-Hilbert anti-self-adjoint representative Dsharp, the actual
paired square domain, smooth thermal bolt and physical-half rank-two trace.
For R_b=(L_b+M^2)^-1, S_b,plus=(D_b+M)^-1 and its adjoint S_b,minus,

    R_b=(S_b,plus+S_b,minus)/(2M),
    S_b,plus=exp(b/2)[Dsharp+M exp(b)]^-1 exp(b/2).       (1)

Both identities hold as bounded resolvents on the actual domain. The
scalar multipliers preserve its first-order projector. They do not assume
self-adjointness of the original unpaired Euclidean D. In the canonical
proper frame, put f=sqrt(Q), dz/d rho=-f, W=p/f, p=gamma(j+1/2).
The retained left/right regular Green construction has coincident trace
T=(1-UV)/(1+UV). Thus at b=0

    M^2 delta_v tr R = M[v T+delta_v T].                 (2)

This derives the density vertex in addition to the mass variation. The
two frequency signs are included explicitly below. Formula (2) supplies
a resolvent identity, not a determinant-normalization guess. In particular
it avoids using the uncut characteristic determinant before matching its
finite thermal subtraction.

Subtract its local flat-frequency term before integration. The remaining
per-frequency trace is integrable at large p. Its common finite part is
identified by the finite-M conversion in current-geometric-mass-conversion:
the flat continuous density cancels the massive conversion contact in
the product-matched action; its finite-temperature difference remains.
The curvature conversion, geometric a2 and all existing finite local
restorations must still be retained. The difference after matching the
principal local term is trace class, so its trace is basis independent.
The already proved selected heat majorant supplies the common angular/
regulator limit. None of this proves a rough-family nonlinear C2 theorem.

## 2. Exact signed radial expression and a finite expansion with a residual

For m=M exp(b), K=sqrt(m^2+W^2), theta=atan(W/m), define the rotated
left/right variables X=tan(atan(U)-theta/2), Y=tan(atan(V)-theta/2).
With the proper derivative toward the seam and c=theta'/2,

    X'=-2K X-c(1+X^2),
    Y'= 2K Y-c(1+Y^2),       Y(seam)=X(seam).           (3)

The exact mass source minus its flat source, for BOTH frequency signs, is

    F_diff = [4m^2 X Y/K+2mW(X+Y)/K]/(1+XY).           (4)

This follows by writing T=cos(alpha+beta)/cos(alpha-beta), not by
estimating two large coincident traces separately. Its v variation,
paired with h=b_h1, is the required signed kernel after the stated local
conversion/subtractions. The exact paired boundary is present in (3).

At b=0 set r=p/sqrt(p^2+M^2 f^2), t=r^2, l=1-t, a=f'/f, a_i=D^i a,
v_i=D^i v, D=d/d rho. A finite radial expansion is

    X_N=sum_(j=1)^N M^-j r l^((j+1)/2) P_j(t),
    Y_N=sum_(j=1)^N (-1)^j M^-j r l^((j+1)/2) P_j(t).

It is an auxiliary approximation, not a different physical domain or a
truncation of physical modes. For a coefficient polynomial P define

    B_j P = D_geo P + a[(j+2)t-1]P -2a t(1-t) partial_t P,
    D_geo a_i=a_(i+1), D_geo v_i=v_(i+1).

The exact coefficient recursion is

    P_1=a/4,
    P_(j+1)=-(B_j P_j)/2
                 +(a/4)t(1-t) sum_(k+l=j) P_k P_l.    (5)

Define V_j by delta_v X_N=sum M^-j r l^((j+1)/2)V_j; the variation is at
FIXED physical p and fixed proper background coordinate. Then

    V_1=[v_1+a(3t-2)v_0]/4,
    V_(j+1)=-(B_j V_j)/2+(1-t)v_0(B_j P_j)/2
      +t(1-t)/4 { [v_1+a(3t-2)v_0] sum P_k P_l
                         +a sum(V_k P_l+P_k V_l)}.     (6)

The spatial derivative commutes with this mass variation. In particular
it is wrong to differentiate only M while freezing r or the metric
vertex. Verify (5)-(6) directly against the exact Riccati equation.

The ENTIRE residual after (5) is

    R_N=M^-N D[r l^((N+1)/2)P_N]
                +c sum_(k+l>=N) M^(-k-l) x_k x_l,     (7)

where x_j=r l^((j+1)/2)P_j. Differentiate (7) to retain the complete
v residual. Coefficient bounds are useful only after propagation by (3),
the true left regular solution and the actual seam matching. They are
not themselves a projection certificate.

## 3. Local cancellations, complete thermal counting and endpoint obligations

Insert the finite expansions into (4), expand only through M^-5, and
retain the exact rational remainder of 1/(1+X_N Y_N). The continuous
frequency integrals of every retained coefficient are exact beta
integrals: dp=M f(1-r^2)^(-3/2)dr. Independently check that the first
coefficient cancels the curvature conversion, the second is the SAME
geometric a2 Hessian, and the third is the retained a6/A3 coefficient.
No local finite coupling or product contact may be reset in this step.

Finite temperature is not set to zero. For the first 80 angular levels,
bound the difference between gamma sum_j F(gamma(j+1/2)) and int_0^infinity
F(p)dp by the complete total-variation bound gamma int|F'|/2. It follows
cell by cell from the fundamental theorem of calculus and includes ALL
frequency blocks and the infinite tail. For the rational flat profile
whose integral needs a finite part, use its proper-time representation
and the already retained flat thermal restoration before removing the
common regulator. The difference has the same convergent variation bound.
Do not assign independent finite parts to divergent frequency sums.

Use an observation split zc=1-2^-49 and an inner computational starting
point zL=1-2^-56. The interval data X(zL) in [-1,1] and its bounded metric
variation must enclose the TRUE regular bolt solution. They are not wall
data. The buffer suppresses their uncertainty on 0<=z<=zc. The original
positive-temperature Q, h and gamma remain correlated and unchanged.

On zc<z<zh bound the selected linear-in-h heat/resolvent observable with
the existing actual cold envelope |h|<=1e-6 s[1+log(.01/s)]. This splits
an integral, not the metric domain: h times an observation characteristic
is not asserted to be a new admissible metric variation. Keep its UV
local coefficient, positive-time part and actual bolt endpoint. Every
paired-seam reflection and differentiated boundary contribution from
Y(seam)=X(seam) must likewise be included. The outer v=-z^4 jets provide
suppression; they are not permission to omit the reflection.

## 4. Decision protocol

First execute the exact local coefficient identities and actual-geometry
residual sizes for N=9 and N=11. The geometry uses x=-log(1-z),
F=f/(1-z), D=-F partial_x and a=F-F_x. Since F^2 is a degree-four
polynomial in y=exp(-x) plus the RETAINED -tau(1-y)^3/y^2, complete
interval Taylor jets give every required derivative on each whole cell.
There is no finite Taylor replacement for the metric.

The first development execution passed all three local identities but its
interval-jet residual coefficients grew from approximately 1.516e8 at N=9
to 8.255e10 at N=11. These are known upper allowances, not values of the
physical response. Increasing the order did not improve that estimator.
Before accepting its verdict, preserve the exact geometry correlations:

    T=F^2, a_i=A_i/T^((i+1)/2), v_i=B_i/T^(i/2),
    A_0=T+y T_y/2, B_0=-(1-y)^4,
    A_(i+1)=y[T A_i'-(i+1)A_i T_y/2],
    B_(i+1)=y[T B_i'-i B_i T_y/2].                    (8)

Primes in (8) mean y derivatives. This is repeated application of the
EXACT D=F y partial_y. Every numerator is a finite Laurent polynomial,
including the actual thermal term. Complete Bernstein ranges on each
cell, plus the explicit negative-power terms, enclose the whole functions.
No model, endpoint or normalization is changed. Compare (8) with the
independent Taylor-jet computation. The second executed diagnostic did
NOT improve the residual allowance (about 2.510e8 at N=9): bounding the
jets separately still discards their correlations inside P_j and V_j.

The next evaluation therefore retains those correlations in the COMPLETE
coefficient. All P_j and V_j are homogeneous of proper derivative weight
j. Write P_j=N_j(y,t)/T^(j/2), V_j=Z_j(y,t)/T^(j/2), and define

    C_j N=y[T partial_y N-j N T_y/2]
              +A_0[(j+2)t-1]N-2A_0 t(1-t) partial_t N.

Then substitute D_geo in (5)-(6) by this numerator operation, a by A_0,
and v_1 by B_1, retaining the common denominator. The resulting finite
bivariate Laurent polynomials are bounded jointly in y and t BEFORE
dividing by the positive T power. Bernstein subdivision covers the
ENTIRE t in [0,1], hence every p>=0. Negative y powers keep their exact
thermal coefficients and actual positive lower endpoint. The positive
polynomial can be bounded on a larger interval reaching zero without
changing the physical endpoint of the full Laurent expression.

The first joint-enclosure execution stopped before a coefficient result:
forming a wide midpoint-radius ball for y in [2^-56,2^-12] rounded its
lower endpoint past zero, so inverting that ball produced NaN. This is
an enclosure implementation failure, not a zero of the actual T. For
negative powers bound the monotone function by its separately evaluated
POSITIVE endpoint values before taking the hull. The failed run remains
in the receipt. No physical interval or thermal term is removed.

No complete projection claim is permitted until the Riccati residual,
its metric variation, boundary reflection, inner buffer, cold observation,
complete thermal difference and angular tail have all been enclosed.
The old Q_(P64) propagation error is not an error for this new evaluator.
If a full signed interval centered at zero is narrow enough to pass
|P_hat_rest|+E_rest<1.287e-40, it certifies that sufficient screen without
determining the sign. Otherwise report the actual interval and the
dominant unknown or excessive allowance, not an existence obstruction.

No cold/scalar/local-contact refinement is selected. Nonlinear admission,
rough nonlinear action/domain, exact compatibility, physical matching,
causal realization, semiclassical validity/noise and all-sector stability
remain independent requirements for the full stable-wormhole objective.

## 5. Propagate the residual on the original domain

The accepted joint diagnostic bounds the entire N=9 residual coefficients
by R<3585 and Rv<9001, instead of the earlier 1.5e8 variation allowance.
The geometry guards are |a|<1.01, |Dv|<.5, f<.409, F<1.01 and bulk proper
length <35.023. The f guard follows from Q<=beta+u/16<.409^2; the negative
epsilon and positive thermal coefficient have their original signs.
These are upper bounds, not source values. The decisive remaining test is
the propagation and complete source accounting below.

Set A=1.01, c0=A/4, cv0=(.5+A)/4, M1=sqrt((q+1)/C), and
J=(1/2)sum B_j/M1^j, Jv=(1/2)sum D_j/M1^j, where B_j,D_j bound P_j,V_j.
Both are <.1. Write lmax=[1+(p/(M*.409))^2]^-1. The exact residual bounds
are |R_N|<=R lmax^((N+1)/2)/M^N and the analogous bound with Rv for its
variation. Define

    eta=c0 max(1+J,2J)/(2M1), b0=2(1-eta),
    Ca=[1+(cv0(1+J)+5c0)/(2M1)]/(1-eta),
    E=R lmax^((N+2)/2)/[2(1-eta) M^(N+1)],
    Ev=Rv lmax^((N+2)/2)/[2(1-eta) M^(N+1)].

For e=X-X_N the exact equations, before inequalities, are

    e'=-[2K+c(X+X_N)]e-R_N,
    ev'=-[2K+2c X_N]ev
         -[2Kv+cv(X+X_N)+2c Xv]e-delta R_N.

The true global bounds |X|<=1, |Xv|<=2.5 follow from 0<=U<=1,
|Uv|<=1 and the rotation. The damped variation equation therefore gives
left forced errors E,Ca E+Ev. The right equation, including the left
forced error transmitted at the seam, gives 2E,2Ca E+Ev. This is a
comparison for the actual variable coefficients, not a constant-flow
replacement. The boundary mismatch of the formal expansions is accounted
for separately below.

The true smallness used in the source Lipschitz bound is stronger. The
buffer length is at least 7 log(2)/1.01. On its first half, damping
2M-c0 brings |X| below B=.008; on its second half, damping 2K-2c0 B
and forcing <=2KB+cv0(1+B^2) bring |Xv| below Bv=.02. The explicit
exponential comparison includes initial bounds 1 and 2.5. These bounds
remain invariant through the observed bulk and the inward solution
initialized by the actual seam data.

For phi=(4sqrt(l)XY+2r(X+Y))/(1+XY), on |X|,|Y|<=.1 its first
derivatives are <=3, second derivatives <=6, and the coefficient-metric
variation of a first derivative is <=3|v|. For example
phi_X=[4sqrt(l)Y+2r(1-Y^2)]/(1+XY)^2. The mean value theorem gives

    |delta F-delta F_N|/M
      <=3(|evX|+|evY|)+(12Bv+3)(|eX|+|eY|).

Thus the forced residual coefficient is

    Kerr=[(9Ca+9+36Bv)R+6Rv]/[2(1-eta)].

For e>1/2 put J_e=sqrt(pi)Gamma(e-1/2)/(2Gamma(e)). Monotonicity gives
the COMPLETE thermal sum, including every half-integer frequency,

    gamma sum_(j>=0) lmax^e <=M*.409*J_e+gamma/2.       (9)

With H=1.2e-8, L the complete observed proper length and d_n=2(q+2n),
the first 80 angular levels have forced error

    E_rad <=G H L Kerr/(pi B^2 F_lambda)
       [ .409 J_((N+2)/2) sum_(n=1)^80 d_n/M_n^(N-1)
                            +(gamma/2) sum d_n/M_n^N ]. (10)

The angular tail is still present separately. No finite set of physical
thermal modes is evaluated or substituted in (9).

## 6. Entire rational remainder and thermal difference

Let A_+(u)=sum B_j u^j and V_+(u)=sum D_j u^j. In each product below
keep its even powers. Put Q_+=(A_+^2)_even, Qv_+=(2A_+ V_+)_even,
S_+=2sum_even B_j u^j, Sv_+=2sum_even D_j u^j and

    N_+=4Q_++2S_+, Nv_+=4Qv_++8Q_++2Sv_++2S_+,
    H_+=Nv_+/(1-Q_+/4)+N_+ Qv_+/[4(1-Q_+/4)^2].

This is a nonnegative coefficient majorant for the full rational source.
At u=1/M1 require Q_+/4<1. Subtract only its coefficients of u^2,u^4,u^6:

    Ktail=M1^8[H_+(1/M1)-sum_(k=1)^3 H_+,2k/M1^(2k)].

Positivity bounds the ENTIRE remaining series at u=sqrt(l)/M. Its source
is <=Ktail r^2 l^(9/2)/M^7. Integrating the profile uses
int_0^1 r^2(1-r^2)^3 dr=16/315. Its total variation in p is <=2,
so the full sum-to-integral error is <=gamma. Consequently

    E_app <=G Ktail/(pi B^2 F_lambda)
          [ H(16/315)sum d_n/M_n^6+gamma H L sum d_n/M_n^7 ]. (11)

The retained local coefficients are polynomials H_2k(t) of degree <=12.
The interval coefficient majorant bounds their sup norm. Markov's
inequality on [0,1] gives |H_2k'|<=288 sup|H_2k|. The product with
r^2 l^(k+1/2) has total variation <=290 H_+,2k. The midpoint error is
therefore <=145 gamma H_+,2k M^(1-2k), integrated against H L. The flat
metric-Hessian profile has total variation <=6M, so add 3gamma M.
These finite-temperature errors are summed for all 80 angular levels.
For the flat term take the common proper-time cutoff first; the variation
bound controls the convergent difference with its continuous integral.
Its matched local finite part is precisely the conversion already fixed
in section 1. No separate zeta assignments or zero-temperature replacement
are used. The complete higher angular contribution is bounded separately.

The first two continuous integrals cancel exactly the curvature conversion
and the same a2 allocation. The third is O6/4. The existing whole-domain
product-subtracted A3 bound gives the raw bound
Uraw=2.3454197e-13+6.704e-12/(63pi). Its normalized allowance is
G Uraw sum d_n/M_n^4/(B^2 F_lambda). Add the cold piece when using this
whole-domain bound for the bulk, and keep the full product contact <=2e-48.

## 7. Actual seam reflection and the regular-bolt uncertainty

Compute all seam jets by the original physical-side polynomial Q, with
Q(0)=beta, Q'(0)=0 and v=-z^4. Then a0=0, v0=v1=v2=v3=v5=0.
The formal mismatch is 2sum_odd x_j/M^j: its background starts at j=3
and its metric variation at j=7. These zeros are algebraic, not rounded
small numbers. All higher odd terms of the finite expansion are retained.

Let d be proper distance inward from the seam. Since z<=.409d,
|v|<=.409^4 d^4. The exact mass Riccati sensitivity comparison gives

    |Uv|<=.409^4[d^4+2d^3/M+3d^2/M^2+3d/M^3+1.5/M^4],
    |Vv|<=.409^4[d^4+1.5/M^4].

The first follows by integrating (d+s)^4 against 2M exp(-2Ms) from
the regular side; the second includes the transmitted seam value and
the inward forcing. Rotation gives both |Xv|,|Yv|<=.409^4 P(d), where
P(d)=2.5d^4+4d^3/M+6d^2/M^2+6d/M^3+3/M^4.
Also |cv|<=.409^4[d^3+A d^4/4]. Put

    Qb(d)=(B+J)[d^4/4+A d^5/20]
            +2c0[d^5/2+d^4/M+2d^3/M^2+3d^2/M^3+3d/M^4].

If b, bv are a formal boundary mismatch and its variation, the propagated
source error is at most

    M exp(-I(d)) {3|bv|+|b|*.409^4[
                  3d^4 I(d)/(1-eta)+3Qb(d)+12P(d)+3d^4]},
    I(d)>=2(1-eta)Kmin d, Kmin=M/sqrt(lmax).

This follows from the exact error equations in section 5, with the
inhomogeneous terms integrated before taking bounds. Integrate every
distance monomial to infinity using k!/[2(1-eta)Kmin]^(k+1).
For d^4 I exp(-I), use I exp(-I)<=(1+a)exp(-a), I>=a, giving
144/[2(1-eta)Kmin]^5. Apply (9) to each resulting frequency profile.
This includes the actual reflected boundary variation, not a new Robin
condition. The upper integration extension only enlarges an error bound.

At the inner computational point keep true unknown data |X|<=1 and
|Xv|<=2.5. With A_min=2(1-eta)M1*7log(2)/1.01>13, both exp(-I) and
I exp(-I) retain an lmax^6 majorant for I>=A_min Kmin/M1.
Initial errors are <=1+J and 2.5+Jv. A factor four bounds the transmitted
inner error and its subsequent inward propagation, including its
variation. Integrate this majorant by (9). This is suppression of true
regular-bolt uncertainty, not the imposition of a wall datum.

## 8. Complete cold observation, spectral tails and the final test

On z>=zc use Hc=1e-6 sc[1+log(.01/sc)], sc=2^-49/(1-2^-49).
This is a bound for the original observation h. The operator, domain and
metric variation are unchanged. The interior selected heat defect is
pointwise, so Bcold=1101 Hc, while
Ucold<=Hc*2^-49 sup|O6prod|/(4pi) and Ic<=Hc*2^-49.
Use the earlier complete positive-time tail formula with N=0, H=Hc,
these Ucold,Ic and the unchanged trace anchor 1521 at .00025.
The all-angular bound majorizes the first 80 levels. Add the cold UV
allowance 4GC[Ucold*.001+Bcold*.001^2/2]/(B^2 F_lambda).
The local kernels s1,s2 and the finite-time boundary term are retained.
Add Cprod Ic to convert this cold product-subtracted response back to
the cold actual response. This is accounting for an observation split;
it does not introduce a discontinuous admissible metric direction.

Assemble exactly once: entire local rest <1.875137e-42; n>80 UV subset
<1.5483364e-43; n>80 positive-time <4.624e-45; low-angular product contact
<2e-48; and equations (10)-(11), A3, thermal, seam, inner and cold errors.
The n>80 UV estimate uses the old complete sum only as a positive upper
bound on the OMITTED subset. Low-level UV is already included in the
radial calculation and is not added again.

Take P_hat_rest=0 and E_rest equal to this complete absolute allowance.
A pass is a certified signed INTERVAL containing the actual projection,
without measuring its sign. It establishes the user's sufficient screen
only. If it fails, the allowance is excessive; it is neither an unknown
number nor a lower bound obstructing a wormhole. The decisive execution
must report every component and the strict inequality, not just check counts.

Remaining to-do list: execute and certify the complete quantum/contact
projection, then the actual controlled correction, nonlinear domain and
admission, exact compatibility, matching, validity/noise and all-sector
stability; preserve the full opening/transport/preparation/creation goal.
