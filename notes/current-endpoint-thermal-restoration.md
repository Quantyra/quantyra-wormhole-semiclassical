# Endpoint normalization and the complete flat thermal restoration

2026-10-23, E003/S169. Ordinary argument before frozen execution.
The stable semiclassical wormhole objective remains incomplete.

The selected calculation is the physical weak H e=-r1 solve. The accepted
compact thermal difference cannot simply be assigned to bolt-reaching
directions: its finite-cut initial sensitivities were zero. This supplies
the missing endpoint identities and a reusable complete thermal restoration
for the SAME canonical mass action. It does not assemble the remaining
angular/vacuum/Weyl terms or compute the metric correction.

## 1. Normalize at the actual bolt, not at a finite cut

Keep the current seed, fixed period, paired seam, original scale and full
Dirac tower. In canonical proper distance, a(rho)=rho+O(rho^3),
W=omega/a, omega=j+1/2, and m=M exp(b), with b smooth and even at the
bolt. The mass-vertex variations here hold the canonical seed metric a
fixed. The Weyl/physical metric conversion remains a separate retained
part of the action. They are not a claim about all metric variations.

Use x'=m x+W y, y'=W x-m y, R=sqrt(x^2+y^2), U=y/x. Fix the regular
branch by rho^-omega R -> 1, independent of b. The equations for
r=x+y and s=x-y are r'=Wr+ms, s'=mr-Ws. A Volterra construction of
the regular branch (factor rho^omega in r, then integrate the s equation
with integrating factor rho^omega) gives, with m0=M exp(b(0)),

    s/r=m0 rho/(1+2omega)+O(rho^3),
    U=1-2m0 rho/(1+2omega)+O(rho^2),
    X=tan(atan(U)-atan(W/m)/2)
      =m0 rho/[2omega(1+2omega)]+O(rho^3),
    log R=omega log rho+O(rho^2).                         (1)

Factoring the regular powers gives integral equations with bounded smooth
coefficients on a small disk; iteration and parameter differentiation
give the displayed expansion and its first two mass variations. In
particular Y_i=D_i X and Z_ij=D_iD_j X have the same leading coefficient
as X, multiplied by u_i(0) and u_i(0)u_j(0). They tend to zero at the
true bolt but are generally NONZERO at a computational epsilon.

In the exact flat-polar comparator a=rho, constant m, let j=omega-1/2.
The regular solution is

    r=C sqrt(rho) I_j(m rho), s=C sqrt(rho) I_(j+1)(m rho),
    C=sqrt(2) j! (2/m)^j.                               (2)

The series and recurrence in [DLMF 10.25.2](https://dlmf.nist.gov/10.25#E2)
and [10.29.2](https://dlmf.nist.gov/10.29#E2) verify the normalization and
the differential equations. The coefficient C must be differentiated.
Omitting it loses +2j u(0) from the paired action's first variation.
Summing that error over j is divergent. This is not a finite-counterterm
choice or a reason to reject the correctly normalized construction.

Define the flat primitive with its mass-independent singular subtraction:
Phi(rho)-omega log rho -> 0, Phi'=K=sqrt(m^2+W^2). Thus

    A_diff=-2 log R+2 Phi=4 integral_0^L K X^2/(1+X^2) d rho. (3)

At the paired seam L this is the same characteristic variation as before.
Both log R and Phi have mass-direction derivatives tending to zero at
the bolt. No new inner contact is introduced. In (2), explicitly,

    A_diff(rho)=m^2 rho^2/[2omega(1+2omega)^2]+O(rho^4).

The zero in (3) refers to the normalized finite part, not a wall condition.

## 2. The auxiliary adjoint and finite-cut terms

Retain the exact sensitivities and auxiliary action adjoint of
[the thermal difference](current-subtracted-thermal-action.md):

    a_Y=2K+2cX, ell'=a_Y ell-4K g'(X), ell(L)=0,
    g(X)=X^2/(1+X^2), c=theta'/2.

At omega=1/2, (1) gives

    ell=rho[C_ell-2m0 log rho+o(1)].                    (4)

For omega>1/2 its leading particular term is
4m0 rho/[(1+2omega)(2omega-1)], with a homogeneous O(rho^(2omega))
term and higher terms. The local regular-singular integration proves
ell Y_i and ell Z_ij tend to zero, including the lowest-channel logarithm.
This is not the singular gravitational bolt adjoint or its delta source.

If numerical adjoint integration begins at epsilon>0, integration by parts
instead gives the REQUIRED terms

    grad_i = tip_grad_i + integral_epsilon^L G_i
                              +ell(epsilon)Y_i(epsilon),
    Hess_ij= tip_Hess_ij+ integral_epsilon^L H_ij
                              +ell(epsilon)Z_ij(epsilon).       (5)

G_i,H_ij are the already derived grouped integrands. Dropping the last
terms is valid for the earlier compact prefix, where sensitivities are
identically zero, but not for a physical bolt direction. A small omitted-tip
bound does not also bound these boundary terms. The returned initial
boxes and the backward ell enclosure must be used in (5).

At the outer seam, the allowed one-sided mass/log direction is
u=c3 z^3+z^4(c4L log z+c4)+a smoother weighted remainder. It and its
first derivative are continuous on the paired extension. The first-order
transfer and the mass-action formulas need only these derivatives; their
finite-mode integration by parts has no seam term since ell(L)=0.
The paired endpoint characteristic still varies. This establishes the
endpoint identity for such directions, not a smooth closed C5 domain or
the common angular heat expansion. The physical mass/log boundary action
and canonical-to-York dual map remain as in the closed audit.

## 3. A uniform all-thermal bound on the omitted bolt cap

At b=0, assume 0<rho<=epsilon, rho/2<=a<=2rho, |a'|<=2,
M epsilon<=1e-3, and |u_i|<=1, epsilon |u_i'|<=1. The derivative is
proper radial. Smooth bolt directions meet the last condition on a
sufficiently small cap; it must be checked for the actual basis.
Write Qw=omega(1+omega), delta=M epsilon. The exact coefficient formulas
with f replaced by a and p by omega give

    omega/(2rho)<=K<=3omega/rho,
    |c|<=M/omega, |c_i|<=2M/omega, |c_ij|<=3M/omega,
    |K_i|<=2M^2 rho/omega, |K_ij|<=4M^2 rho/omega.       (6)

For example |1-2lambda|<=1 and |1-8lambda+8lambda^2|<=1 on
lambda in [0,1] bound the differentiated c coefficients. Signed scalar
comparison from the regular initial limit yields

    |X|<=2M rho/Qw, |Y_i|<=10M rho/Qw,
    |Z_ij|<=14M rho/Qw.                                (7)

Details of the invariant barriers: |X|<=1 follows from 0<=U<=1 and
0<=tan(theta/2)<=1. The X forcing is at most 2M/omega, and 2K>=omega/rho.
The stated linear barrier solves that comparison exactly. Then
a_Y>=omega/(2rho), since 4delta^2/[omega^3(1+omega)]<1/2.
The Y forcing divided by M/omega is at most
4+8delta^2/Qw<5. The Z forcing has allowance
6+96delta^2/Qw+360delta^2/Qw^2<7. The regular zero initial limit and
the barriers 10M rho/Qw and 14M rho/Qw close the two comparisons.
The same mixed equation has these bounds when both direction norms are 1.

Using |g|<=X^2, |g'|<=2|X|, |g''|<=2 in the ORIGINAL variational
integrand gives gradient coefficient <=480+32delta^2/omega^2 and
Hessian coefficient <=3072+704delta^2/omega^2. Both are below 3100:

    |integrand_i|, |integrand_ij|
          <=3100 M^2 rho/[omega(1+omega)^2].            (8)

The thermal sum is exact:

    sum_(j>=0) 1/[omega(1+omega)^2]=6-pi^2/2<1.1.

Partial fractions telescope the first two terms and the odd reciprocal
squares supply the last one. Consequently the entire omitted thermal
cap in each gradient/Hessian entry has allowance

    E_diff<=1550(6-pi^2/2)(M epsilon)^2.                (9)

This is the unweighted paired-frequency action, per angular state, with
all j included. Multiply by the original gamma if using the previous
publication units. That conversion changes no physical error budget.
The formula applies to any M with the stated cap/geometry hypotheses;
it is not uniform in M at a fixed epsilon without those hypotheses.

On the ACTUAL finite seed use the already factored chi geometry, with
y=2sinh(chi/2)^2, t=1+h(1+y), Bchi=2K(v)v+yH1(t,v):

    a=sqrt(y Bchi)/(K(v)t^(3/2)),
    d rho/d chi=sqrt((y+2)/(t Bchi)).                  (10)

Enclose a/chi by sinh(chi/2)/(chi/2) in [1,cosh(.001/2)] and
differentiate the nonsingular factor sqrt(2)sinh(chi/2)sqrt(Bchi)/
(K(v)t^(3/2)). A whole interval chi in [0,.001] verifies the geometry
guards in (6), including a/rho through the integral bounds on d rho/d chi.
The current h box is inherited, not replaced by zero. This transfer is
checked before inserting the benchmark M in (9).

## 4. Restore the complete local thermal difference

The subtracted flat action must be put back. Its regulator-independent
thermal part is the nonzero-winding part of the SAME proper-time sum,
with the zero-temperature part retained for the common four-dimensional
renormalization. At a>0 and y=ma>0, define

    T(m,a)=-(4/a) integral_y^infinity
                   sqrt(t^2-y^2)/(exp(2pi t)+1) dt,
    S0(y)=integral_0^infinity [exp(2pi y cosh v)+1]^-1 dv,
    S1(y)=-y S0'(y).

Then the first two exponential-mass derivatives are exactly

    T_b=4m^2 a S0(y),
    T_bb=4m^2 a[2S0(y)-S1(y)].                         (11)

For two directions the local Hessian is T_bb u_i u_j, including the
exponential-mass tadpole. These are not estimates of a raw frequency sum.

Here is the regulator and summation justification. Apply
[DLMF 1.8.14](https://dlmf.nist.gov/1.8#E14) to the shifted Gaussian;
the Gaussian and its second derivative are integrable. This gives

    sum_(j>=0) exp[-s(j+1/2)^2]
      =sqrt(pi)/(2sqrt(s)) [1+2sum_(k>=1)(-1)^k exp(-pi^2 k^2/s)].

Subtract the displayed zero-winding term INSIDE the proper-time integral.
For y>0, the integrals of the absolute winding terms converge; this
justifies interchange. The positive-real K integrals in
[DLMF 10.32.9-10](https://dlmf.nist.gov/10.32#E9) give

    aT=-(2y/pi) sum_(k>=1)(-1)^(k+1) K1(2pi k y)/k,
    S0=sum_(k>=1)(-1)^(k+1) K0(2pi k y),
    S1=2pi y sum_(k>=1)(-1)^(k+1) k K1(2pi k y).      (12)

The Fermi integral follows by summing the exponential series under the
same convergent integrals. The derivative identities are those of
DLMF 10.29. The mass-independent limit T(0,a)=-1/(12a) is RETAINED;
subtracting it only when bounding mass variations does not erase it from
the metric/Weyl action. The proper-time definition avoids assigning an
unjustified value to a divergent hard-cutoff energy sum.

For an efficient rigorous implementation near y=0, the convergent
frequency finite part gives

    S0=-[EulerGamma+log(2y)]/2+sum_(k>=1) c_k y^(2k),
    c_k=-((-1)^k binom(2k,k)/4^k)(2^(2k+1)-1)zeta(2k+1)/2,
    S1=1/2-sum_(k>=1)2k c_k y^(2k),
    aT=-1/12+y^2[1/2-EulerGamma-log(2y)]
                          +sum_(k>=1)2c_k y^(2k+2)/(k+1). (13)

Indeed the half-integer harmonic sum minus integral tends to
EulerGamma+log(2y); the remaining binomial series is absolutely convergent
for 2y<1, so summing it yields (13). The Gaussian representation above
identifies this finite part with (12), including its constant.

For y<=1/8 put r=4y^2. After K terms, use
E0=1.203 r^(K+1)/(1-r) for S0,
E1=2.406 r^(K+1)[(K+1)/(1-r)+r/(1-r)^2] for S1, and
ET=2y^2 E0/(K+2) for aT. These follow from |binomial coefficient|<=1
and zeta(2k+1)<=zeta(3)<1.203. For y>=1/16, use (12) and an
alternating-series remainder bounded by the next term. K0(x), xK1(x)
and K1(x)/x are positive decreasing for x>0 by their positive integrals
and derivative identities. Both algorithms therefore have analytic tails.
Their overlap checks implementation, not the proof of summability.

## 5. The thermal restoration at the tip and angular limitation

The positive integrals imply S0(y)<=log_+(1/y)/2+1 and 0<S1(y)<=1.
For the first bound use cosh v>=exp(v)/2 and split its resulting
integral at 1. For the second, S1<=2pi y K1(2pi y)<=1; the substitution
t=2pi y sinh v bounds the latter integral by integral_0^infinity exp(-t)dt.
If 2M epsilon<=1 and the geometry guards hold, integration of (11) gives

    E_thermal,H<=4M^2 epsilon^2[log(2/(M epsilon))+7/2],
    E_thermal,G<=2M^2 epsilon^2[log(2/(M epsilon))+5/2]. (14)

For unit direction bounds the Hessian allowance also exceeds the gradient
allowance. Substitute the known M<=26 into (9),(14) at epsilon=5e-9.
Require their combined omitted-tip allowance <1e-10 in UNWEIGHTED action
units. This is an assembly precision, not the physical bolt-error target.
For a general M and desired tolerance, (9),(14) provide the explicit
cap threshold; epsilon must also meet the geometry and direction guards.

Neither (9) nor (14) can be summed absolutely over the angular tower at
a fixed cap. The known allowance contains M_n^2, while
d_n=2(q+2n), M_n^2=n(n+q)/C; hence d_n M_n^2 grows cubically.
The fixed-cap hypothesis itself eventually fails. This rejects that
certificate, not the actual commonly subtracted angular action.
The bolt and high-angular local terms must be combined in a matched
parametrix BEFORE taking a uniform angular remainder. A constant mass
direction tests its required endpoint cancellation; its finite common
part must be restored in the original scheme. No new finite counterterm
or normalization can be used to make this issue disappear.

## 6. Execution and decision

Freeze this argument, the implementation and source review before running.
Verify the regular-branch/finite-cut identities, current full-box geometry,
all-thermal cap constants and the <1e-10 combined cap target. Execute the
restoration evaluator using K=16 near zero and a certified next-term
tail <1e-22 in (12); require dimensionless T, gradient and Hessian widths
<1e-12 on the stated validation inputs. Use an independent 60-digit Fermi
integral and normalized flat-polar Bessel solution as diagnostics only.
The Fermi diagnostic compares within absolute 1e-45, not within arbitrarily
narrow Arb roundoff balls. Its inputs are y=1e-40,1e-8,.01,1/16,1/8,
1/4,1/2,1,4,20; both rigorous algorithms overlap at 1/16,.1,1/8.
The flat-polar endpoint diagnostic uses j=0,1,8,100, M=20 and rho=5e-9.
The first run exposed a diagnostic truncation error: stopping the Fermi
integral at exponent 100 leaves an action tail about 3.81e-43 at y=1e-40,
above the fixed 1e-45 comparison target. Its source and numerical result
are retained. Increase that diagnostic endpoint to exponent 160 and replay
the full package; no analytic enclosure or acceptance target is weakened.
The evaluator's coverage follows from the explicit analytic tails, not
from testing finitely many y or j. Retain every run and snapshot.

Continue the physical weak solve. This removes the compact-only endpoint
identity restriction for the canonical mass variations, implements the
full local thermal restoration, and supplies a usable all-frequency tip
bound. It does not execute the actual coupled correction or prove a
complete global operator. The next combined operator operation is the
matched bolt/angular subtraction with the original Weyl and common
vacuum/a0/a1/a2 restoration, on the physical mass/log phase. Keep both
gravitational point sources in the dual residual. No additional isolated
finite-mode sample or cold-value refinement is selected.

The <1.667e-9 physical bolt difference, <6.67e-10 exterior return,
Y0/Z0/Z2 and all X/admission inequalities remain uncomputed. The audit
stays closed. Matching, causal realization, semiclassical validity/noise
and all-sector stability remain requirements of the full objective.

Remaining to-do list: matched bolt/angular subtraction and original-scheme
restoration in the physical weak operator; actual coupled correction and
certified dual error/X admission; nonlinear closure and exact compatibility;
matching, causal realization, validity/noise, stability, charge/return,
opening, transport, preparation, separate creation and stable wormholes.
