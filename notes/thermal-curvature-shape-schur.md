# Thermal curvature-shape coupling and compact Schur feedback

2026-10-21, E003/S168. Ordinary argument before focused verification.
The full research goal remains incomplete.

The complete thermal shape inverse can be used with the ORIGINAL
physical curvature space H=L2(dy), U=r(delta R-r R_r p). Its mixed
curvature/shape form extends continuously to H times the full shape
energy space. At each sufficiently large positive Matsubara frequency,
eliminating the auxiliary shape equation produces a bounded, positive,
compact correction to the curvature form, with a minus sign in the
Schur complement. No additional curvature weight is imposed.

The reduced Euclidean form remains unbounded below along normalized
curvature profiles escaping to the outer end. The induced shape energy
tends to zero on this sequence. This is a coupled variational statement
for the seed and its auxiliary shape realization. It is not a Lorentzian
growing mode, a physical constraint solution, or an instability theorem
for a nonlinear wormhole. Negative Euclidean conformal directions alone
do not settle those questions.

## 1. Spaces, reconstruction and the exact mixed form

Keep the selected action, physical charge, HH state, clock, angular
tower, finite contacts and common subtraction. Fix a positive thermal
frequency gamma at least

    Gamma_* = max(Gamma_T, Gamma_R),

where Gamma_T is the threshold of the
[full thermal shape inverse](global-thermal-shape-inverse.md), and
Gamma_R is the real-frequency threshold in
[physical curvature reconstruction](physical-curvature-reconstruction.md).
Both are finite; neither is numerically calibrated here.

Use w=r v and q=r p. Write the two reconstruction maps distinctly as

    R_gamma U=(1/6)K_gamma^-1(F U),
    Q_gamma w=-(1/3)K_gamma^-1 B_gamma w,
    q=R_gamma U+Q_gamma w.

The first formula is its established weighted-kernel extension, not
an application of the unbounded multiplier F on all of H. The bounds are

    ||R_gamma||_(H->H)<=1/(18lambda), ||Q_gamma||_(H->H)<4/5.

Let E_gamma be the closed full shape form space, with squared norm
a_gamma[w]. Its retained estimates include

    N_gamma[w]<=33G a_gamma[w],
    C J2[v]<=200pi a_gamma[w],
    ||w||<=sqrt(40G) a_gamma[w]^(1/2)/gamma,
    ||(Q_gamma w)_y||<=sqrt(10N_gamma[w]/9).                  (1)

The last estimate first holds on the compact core by reconstruction
energy, then passes to E_gamma. Define b_gamma(U,w) initially for
compact U,w as the full Hessian pairing of (p_U,0) with (p_w,v),
where p_U=R_gamma U/r and p_w=Q_gamma w/r. The canonical determinant
is independent of the four-dimensional conformal variable, so only
the complete local classical/Weyl action enters this mixed pairing.
Its canonical contribution to a_gamma is nevertheless retained.

Put epsilon=G/(60pi) and let calT be the exact seed equation trace,
including its nonzero residual. The polarized density identity from
[the charged shape form](charged-curvature-free-shape-form.md) gives

    b_gamma(U,w)=-(1/(2G)) <R_gamma U,D_gamma w>,
    D_gamma w=epsilon g_w
              +F[(r calT_r+4calT)q_w+2calT w],
    q_w=Q_gamma w.                                         (2)

The pairings use dy and the actual time-mode norm. For complex fields
use the Hermitian polarization of this real thermal identity.
Every residual, chart and volume contact in (2) remains present.

The exact trace source from the physical reconstruction is

    g_w=A1 U_v(v)/6-B1(q_w)_y-B0 q_w-2A0 v,
    U_v(v)=r[2F^-1(gamma^2 v-v_yy)+C_v v],
    A1=rF d_R, A0=rF(d_R F_r+e_R F),
    B1=2rF e_R,
    B0=2e_R F^2+rF e_R F_r+F d_R L,
    L=F_r+rF_rr+r(2R+rR_r)/6.

Here d_R,e_R,C_v are the explicit unchanged seed coefficients in
[the metric source](explicit-curvature-metric-source.md). Rearranging,

    D_gamma w=-(epsilon r^2 d_R/3)(v_yy-gamma^2v)
               -epsilon B1(q_w)_y+M_w w+M_q q_w,
    M_w=epsilon(A1 C_v/6-2A0/r)+2F calT,
    M_q=F(r calT_r+4calT)-epsilon B0.                       (3)

This uses the exact trace before estimating the mixed Hessian.
Estimating its fourth-order terms separately would unnecessarily demand
extra derivatives or weights of U. Equation (2) needs only R_gamma U
in H, even when that conformal lift is outside the unweighted energy
domain of K_gamma. The mixed form for such U is the unique continuous
limit below; no individually divergent fourth-order summand is assigned
an independent value.

## 2. Bounded and compact coupling in physical curvature L2

Define finite nonnegative constants

    K2=(1/3)||r sqrt(F) d_R||_infinity,
    K1=sqrt(10/9)||B1||_infinity,
    K0=||M_w||_infinity+(4/5)||M_q||_infinity.

Since C J2[v]=integral (r^2/F)(v_yy-gamma^2v)^2dy, (1)-(3) give

    ||D_gamma w||<=epsilon K2 sqrt(C J2[v])
                   +epsilon K1 sqrt(N_gamma[w])+K0||w||,
    |b_gamma(U,w)|<=L_gamma ||U|| a_gamma[w]^(1/2),
    L_gamma=[epsilon K2 sqrt(200pi)+epsilon K1 sqrt(33G)
                         +K0 sqrt(40G)/gamma]/(36G lambda). (4)

These are finite definitions, not evaluated numerical error bars.
The bound is uniform along the stated large thermal ray after replacing
gamma in its last denominator by Gamma_*.

The endpoint behavior also makes the coupling compact at EACH fixed
admitted gamma. At the outer end

    d_R=O(r^-4), e_R=O(r^-6), C_v=O(1),
    B1=O(y^3), B0=O(y^2),
    calT=O(r^-4), r calT_r+4calT=O(r^-5).

Thus r sqrt(F)d_R, B1, M_w and M_q all tend to zero. At the
nonextremal horizon all these quantities vanish with a positive power
of F. These facts follow from the exact rational seed, with the vacuum
cosmological shift retained; they are not an assumption calT=0.
For cutoffs chi exhausting the open exterior, (3)-(4) therefore imply

    ||(1-chi)D_gamma||_(E_gamma->H)->0.                    (5)

The physical reconstruction R_gamma maps bounded H sets to bounded
H2 on every compact interior interval. Hence chi R_gamma:H->H is
compact by local Rellich compactness, and R_gamma^* chi is compact
by adjunction. Define

    K_gamma,cross=-(1/(2G))R_gamma^* D_gamma:E_gamma->H.

It is a norm limit of the compact maps
-(1/(2G))R_gamma^* chi D_gamma, by (5). Consequently it is compact
and has norm at most L_gamma. Its E_gamma-adjoint
K_gamma,cross^dagger:H->E_gamma is compact too. Compactness here is
not claimed uniformly in gamma, at a boundary trace, or at a nonlinear
background.

## 3. Exact auxiliary shape elimination

Let c_gamma[U] be the pure-conformal full-action quadratic form
evaluated on p_U=R_gamma U/r, initially for compact U. The transformed
two-row quadratic form is

    F_gamma[U,w]=c_gamma[U]+2 Re b_gamma(U,w)+a_gamma[w].

For every U in H the shape equation has the unique weak solution

    w_*(U)=-K_gamma,cross^dagger U,
    a_gamma[w_*]<=L_gamma^2||U||^2,
    ||w_*||<=sqrt(40G)L_gamma||U||/gamma.                   (6)

In particular this is defined for all original physical curvature L2
inputs. The value c_gamma[U] need not be finite for every such input;
its domain is not enlarged by the existence of the shape solution.
On the curvature core, completing the square gives exactly

    F_gamma[U,w]=c_gamma[U]-<U,M_gamma U>
                       +a_gamma[w+K_gamma,cross^dagger U],
    M_gamma=K_gamma,cross K_gamma,cross^dagger,
    0<=M_gamma, M_gamma compact, ||M_gamma||<=L_gamma^2.    (7)

Thus the auxiliary reduced curvature form is c_gamma-<U,M_gamma U>.
This is a Schur reduction of the full congruence, including the Q^*
trace-row contribution. It is not a substitution into only the original
v equation, and does not invert or delete the fast conformal factor.
Nonzero initial jets and independent state/boundary sources still require
the transformed source terms of the retained reconstruction.

The positive closed-form/Riesz framework is standard; see
[Ashbaugh et al., equations (2.18)-(2.21)](https://www.mat.univie.ac.at/~gerald/ftp/articles/KreinWeyl.pdf).
The new content is the exact charged mixed row, its physical-L2 bound,
and its endpoint compactness. The closest prior comparison controlled
conformal feedback with shape held fixed; (7) now includes the response
of the complete auxiliary shape block.

## 4. Negative curvature directions survive the shape response

This section concerns the Euclidean form, not Lorentzian stability.
It supplies a concrete sign consequence of (7) rather than an unevaluated
coupling criterion. Keep gamma fixed in the admitted thermal sequence.
Choose q_delta=delta^-1/2 chi(y/delta), with nonzero real smooth chi
compactly supported in (1,2), and set

    p_delta=q_delta/r,
    U_delta=6F^-1 K_gamma q_delta.

Both are compact. The exact reconstruction gives R_gamma U_delta=q_delta.
As delta tends to zero, with eta=y/delta,

    ||U_delta||^2 ->36lambda^2 Hchi,
    Hchi=integral_1^2(-eta^2 chi''+6chi)^2deta>0.           (8)

The full conformal form, including its density contact, is retained as

    c_gamma[U]=(1/(2G)) integral {
       (6epsilon/F)[q_yy-(gamma^2+F F_r/r)q]^2
       -k_o[(q_y+Fq/r)^2+gamma^2q^2]
       -F(2R+48lambda0)q^2}dy,
    k_o=6+epsilon R/6+(22epsilon/3)j,
    j=-F_rr/4+(F-1)/(2r^2).                               (9)

On the scaled support, F=1/(lambda y^2)+O(1), F/r=1/y+O(y),
F F_r/r=2/y^2+O(1), R=-12lambda+o(1), j=o(1), and
lambda0-lambda=11epsilon lambda^2/6. Define

    Dchi=integral chi'^2, Pchi=integral chi^2/eta^2,
    Achi=integral eta^2 chi''^2, Gchi=Dchi+6Pchi.

Substitution into (9), with all global functions kept until this limit,
gives the exact leading coefficient

    2G delta^2 c_gamma[U_delta]
       -> -6Gchi+epsilon lambda(6Achi+26Dchi-60Pchi).       (10)

The identity integral(eta chi''-2chi/eta)^2=Achi+4Dchi+4Pchi
fixes the coefficients in (10). For a concrete test take
chi=eta^2(eta-1)^6(2-eta)^6 on [1,2], zero outside. It is C5;
all the displayed integrals are rational and Achi<1e6 Gchi.
It can be approximated by C_c^infinity((1,2)) in H2 with the same
strict inequality. Fix one such smooth approximation. The benchmark
epsilon lambda<(4e-17)^2/(180)<1e-35 makes (10) less than -5Gchi.
For sufficiently small delta, (8)-(10) therefore imply

    c_gamma[U_delta/||U_delta||]
           <-[Gchi/(36G lambda^2 Hchi)] delta^-2.          (11)

All constants refer to the chosen smooth profile. The sufficiently
small delta is not numerically calibrated for the physical benchmark.
This is an asymptotic test theorem, not a certified operating scale.

The normalized U_delta have shrinking compact support near zero and
converge weakly to zero in H. Compactness in section 2 implies

    a_gamma[w_*(U_delta/||U_delta||)]->0,
    <U_delta/||U_delta||,M_gamma U_delta/||U_delta||>->0.    (12)

The reduced form in (7) is thus unbounded below even after solving the
complete auxiliary shape equation. No independent cutoff discards its
induced conformal tails. There is no assertion that these variational
directions solve the lapse/shift constraints, obey the physical boundary
ensemble, or yield a pole of the retarded metric system.

The interpretation is especially important here. Euclidean conformal
indefiniteness is established physics, and imposing physical constraints
can change the relevant degrees of freedom. For an explicit primary
comparison, [Griffin--Kosower, pages 2-6](https://lss.fnal.gov/archive/1989/pub/Pub-89-177-T.pdf)
distinguish the covariant conformal problem from their constrained FRW
calculation. Their result is not transferred to this charged higher-order
system, and no conformal contour rotation is adopted. Gauge-invariant
curvature alone does not make arbitrary variations constraint-admissible.

## 5. Verification and remaining problem

The checker verifies the mixed trace normalization and coefficient
rearrangement, endpoint limits, rational dilation integrals and leading
negative margin. Independent charged reconstruction solves compare (2)
with mixed derivatives of the original nonlinear optical classical/Weyl
densities, keeping the nonzero trace residual. Outer-dilation fixtures
check (8)-(10) against the full conformal action with refinement. They
do not evaluate M_gamma or the canonical kernel numerically; compactness,
form-dual extension and the infinite-dimensional Schur theorem remain
ordinary proofs. No Lean campaign or formal continuum certification is
claimed.

Finite-size convergence need not be monotone before the asymptotic regime.
Writing u=1/r, the exact seed gives

    u=lambda y+lambda^2 y^3/3-m lambda^3 y^4/2+O(y^5),
    F^-1=lambda y^2-lambda^2 y^4/3+m lambda^3 y^5+O(y^6),
    F_r/r-(2R+rR_r)/6=6lambda+2m lambda^3 y^3+O(y^4).

With h=-eta^2 chi''+6chi, the relative curvature norm correction is
C2 delta^2+C3 delta^3+O(delta^4), where

    C2=(2/Hchi) integral h[(lambda/3)eta^4 chi''+gamma^2 eta^2 chi],
    C3=(2m lambda^2/Hchi) integral h[-eta^5 chi''+2eta^3 chi].

In the third fixture C2 is about -6.2027 and C3 about 1197.752.
Their opposite signs explain why an initial test at delta=.005 can
have smaller error than delta=.0025. The retained failed run records
that failure. Refined fixtures use .00125, .000625 and .0003125 with
the same convergence and identity tolerances. This is a numerical
asymptotic check, not calibration of the physical benchmark.

The physical shape realization, continuous causal frequencies, initial
constraints, horizon/outer source domain and nonlinear transfer remain
unresolved. The curvature/shape mixed form is now controlled for this
thermal auxiliary realization; the constrained causal curvature equation
and physical stability decision are still open. No fast branch, finite
counterterm or higher-order action term is removed. S168.3-S168.4 and
all physical completion gates remain unchecked.

Remaining to-do list: identify the physical domain and constrained causal
curvature problem; control continuous frequencies and nonlinear transfer;
prove nonlinear existence and signed compatibility; complete exterior/state,
charge/return, opening, all-sector stability, transport, preparation,
separate creation and full audit.
