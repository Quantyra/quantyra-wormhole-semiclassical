# Continuous causal feedback and fixed-phase auxiliary resonances

2026-10-21, E003/S168. Ordinary argument before focused verification.
The stable-wormhole goal remains incomplete.

The actual eliminated strong response V_gamma of the retained charged
seed is real, uniformly bounded and locally operator-norm Lipschitz for
all sufficiently high real Laplace frequencies. The auxiliary closed
shape forms have a common form domain, with locally equivalent norms.
This proves the remaining continuity hypothesis of the
[boundary-phase theorem](coupled-boundary-phase-winding.md).

Consequently EVERY fixed separated real reflecting phase has an infinite
sequence of auxiliary homogeneous frequencies gamma_j tending to infinity,

    gamma_j/gamma_j^0=1+O((gamma_j^0)^-2),
    gamma_j^0=(2/y0) exp[(arg Gamma(ik)+pi/2-theta+j pi)/k],
    gamma_(j+1)/gamma_j -> exp(pi/k), k=sqrt(a-1/4)>0.           (1)

These are zeros of the actual reduced seed operator on the specified
auxiliary domain. No physical phase is selected. The original constraint,
boundary/action and initial/state admissibility problem is still open,
as are nonlinear existence and transfer to a physical geometry. The
frequency thresholds are finite but uncalibrated, with no controlled
semiclassical frequency window. No analytic-pole or simplicity assertion
is made, and (1) is not a charged-wormhole instability theorem.

## 1. Spaces and previously proved bounds

Write s=gamma^2. Distinguish the scalar canonical Hilbert space
H_c=L2(a(Y)dY), with a=sqrt(CF)/r, from physical H=L2(dy), w=r v.
Use the actual doubled orbit and even scalar seam condition. The
[global causal inverse proof](global-continuous-causal-shape-inverse.md)
gives, on the common compact gravitational core,

    H_s=A+s B, A=-a^-1 partial_Y(a partial_Y), B=a^-2,
    H_s>=s b0, J2_s[v]=||H_s v||_c^2,
    ||Bv||_c<=2||H_s v||_c/s,
    ||sqrt(B)v'||_c<=||H_s v||_c/sqrt(2s).                    (2)

The last bound is the positive gradient term in the same graph identity.
It extends to D(H_s) by the scalar graph core, which includes even
profiles nonzero at the seam. It does not assert that the gravitational
core, which vanishes near the seam, is a core for that scalar operator.
Take gamma/kappa_h>4 and all the previous thresholds. The scalar
limit-point endpoints and the fixed even seam condition give the
needed scalar graph core. At the seam a'=0, so the graph-identity
boundary terms vanish also for those nonzero even profiles.

Let Kflat(sqrt(x))=x^2 g(x), G_s=g(H_s). The previous proof gives

    g(x)=g0+integral [1-(1+t x)^-1] dnu(t), dnu>=0,
    0<=x g'(x)<=C/(40pi), ||G_s'||<=C/(40pi s),
    |Qsc'''(s)[v]|<=27C J2_s[v]/(20pi s^3).                  (3)

The positive measure includes the full angular sum. Enlarge the finite
threshold so that g(x)>0 and is bounded above and below by positive
constant multiples of log(1+x^2) for x>=s b0. This follows from the
retained lower bound g(x)>=C log(x)/(60pi)-constant and (3).
All constants below may depend on a compact positive s-interval I in
this enlarged ray. No uniform derivative bound as s tends to infinity
is needed. The already proved uniform bound on V_gamma is retained.

The actual canonical error e_s=D(s)-Qsc(s) obeys |e_s[v]|<=7C J2_s/800.
For the full scalar-comparison shape form a_sc,s and actual form a_s,

    a_s=a_sc,s+e_s, |e_s|<=eta a_sc,s, eta<3/5,
    a_sc,s>=N_s/(33G)+4C J2_s[v]/(75pi)
                                    +C J2_s[p_s+v/3]/(40pi). (4)

All contacts, conformal tails and the original HH response remain fixed.

## 2. The actual canonical error is locally Lipschitz in a graph norm

Fix s0 and a slightly larger compact high-frequency interval J around I.
On the compact core, H_t v=H_s0 v+(t-s0)Bv, so (2) makes all J2_t
norms equivalent on J. Let X be completion in J2_s0. The diagonal error
bound and polarization give bounded Hermitian forms e_s on X, uniformly
for s in J. This alone would not prove norm continuity.

Choose four fixed distinct positive anchors a1<a2<a3<a4 in J. They
need not be thermal because the preceding theorem already controls all
real s in J. Third divided differences give

    0<=D[a1,a2,a3,a4][v]
       =integral dmu_v(t)/product_j(t+a_j)<=K_I J2_s0[v].       (5)

The upper bound follows by writing D=Qsc+e: the divided difference of
Qsc is bounded by its third derivative (3), and the four error values
are bounded forms on X. No separate bound on a raw unsubtracted D is
used. For s in J, the quotient product_j(t+a_j)/(t+s)^4 is bounded
uniformly in t>=0. The actual positive spectral representation then gives

    0<=D'''(s)[v]<=K_I J2_s0[v],
    |e_s'''[v]|<=K_I J2_s0[v].                               (6)

If necessary first choose the anchors in a larger interval, so (6)
holds throughout J. Compact-probe differentiation is justified by the
retained convergent subtracted spectral integral. For a fixed h>0 with
s+-h in J whenever s in I, the central Taylor identity gives

    |e_s'[v]|<=[M_I/h+K_I h^2/6] J2_s0[v],                  (7)

where |e_t[v]|<=M_I J2_s0[v] on J. Integrate (7) and polarize.
Thus e_s is locally Lipschitz as an operator X->X*, uniformly in the
support of the compact probe. This argument removes the distinction
between compact-probe analyticity and the required uniform error bound.

## 3. Transport of the bolt multiplier through the logarithmic weight

We need more than (2) to differentiate a form with an unbounded logarithm.
Put K_P=||4a'^2+2aa''||_infinity, M_a=||a'||_infinity and

    P=-2(B'/B)partial_Y+[2(B'/B)^2-(Delta_rad B)/B].

The exact identities are

    (H_s-P)B=B H_s,
    P=4(a'/a)partial_Y+(4a'^2+2aa'')B.

Equation (2) gives

    ||Pz||<=theta_s ||H_s z||,
    theta_s=4M_a/sqrt(2s)+sqrt(2)K_P/s.                       (8)

Enlarge the finite threshold so theta_s<=1/2. Let T_s=B H_s^-1.
Then ||T_s||<=2/s and

    T_s D(H_s) subset D(H_s), ||H_s T_s z||<=4||H_s z||/s.     (9)

Here is a domain proof avoiding an assumed core for H_s squared. For
z in D(H_s), define q=(H_s-P)^-1 Bz; the inverse exists by (8) and
||H_s q||<=2||Bz||. Multiplication by B^-1=a^2 maps D(H_s) into itself:
its coefficients and first two covariant derivatives are bounded, its
derivative vanishes at the seam, and it preserves the limit-point ends.
The differential identity gives H_s(B^-1q)=z. Hence q=B H_s^-1z,
proving (9). Graph closure extends the earlier estimates used here.

Let L_s=log(1+H_s^2). Normalize U=T_s/(4/s). Equations (2),(9) imply
U*U<=I and U*H_s^2 U<=H_s^2 as forms. Contractive operator Jensen for
the concave function log(1+x), whose value at zero is zero, gives

    U* L_s U<=log(1+U*H_s^2 U)<=L_s.                         (10)

For the unbounded positive operator this follows by its variational
resolvent compression identity and integration, or monotone bounded
resolvent regularization. This is the concave version of
[Hansen--Pedersen, Corollary 2.3](https://arxiv.org/pdf/math/0204049).
Consequently T_s is bounded on the logarithmic form domain of H_s.
The equivalence following (3) gives

    ||G_s^(1/2) T_s z||<=K_I ||G_s^(1/2) z||.                (11)

No commutation of B with H_s or its logarithm has been assumed.

## 4. Localize the logarithmic energy before differentiating it

Use the same chi_o,chi_i with sum chi_j^2=1 and transition
.01<|Y|<.02. The commutators [G_s,chi_j] are bounded on H_c, locally
uniformly in s. To show this directly, for R_t=(I+tH_s)^-1,

    [1-R_t,chi]=t R_t[H_s,chi]R_t,
    ||[H_s,chi]z||<=K_I ||H_s^(1/2)z||.

The coefficients of the first-order commutator are confined to the
fixed smooth transition. Thus its resolvent integrand is bounded by
K_I sqrt(t)/(1+t s b0). Equation (3), evaluated at x=1/t0, gives
nu([t0,2t0])<=9C/(80pi), because x t/(1+x t)^2>=2/9 there.
Summing the resulting convergent dyadic bounds at t=0 and infinity
proves boundedness of the commutator. In particular

    E_s=G_s-sum chi_j G_s chi_j=sum chi_j[chi_j,G_s]

is bounded on H_c. This controls both signs of the localization error.

Set h=H_s v, W=omega+19/60, c0=C/(30pi), and

    F_s[v]=<h,(G_s-c0 W)h> >=4C J2_s[v]/(75pi).               (12)

The outer uncertainty bound applies to chi_o h itself. It gives a
positive lower bound for <chi_o h,(G_s-c0 W)chi_o h>.
On chi_i, W is bounded. The bounded E_s therefore implies

    ||G_s^(1/2)chi_i h||^2<=F_s[v]+K_I J2_s[v].              (13)

The outer multiplier chi_o B has bounded derivatives and vanishes near
the bolt. Multiplication by it maps D(H_s) to itself, locally uniformly
in s. Since |g(x)|<=K_I(1+x),

    ||G_s chi_o Bv||+||W chi_o Bv||<=K_I ||H_s v||.           (14)

For the W term use the one-dimensional Hardy/first-derivative estimate
near Y=0: compact gravitational inputs vanish there and
|v(Y)|<=sqrt(Y)||v'|| locally. The integrable weight Y log^2(Y)
then gives (14). No boundedness of W itself at the outer end is asserted.

For the interior, write

    chi_i Bv=T_s chi_i h+T_s[H_s,chi_i]H_s^-1 h.

The second input has bounded H_s^(1/2) norm in terms of ||h||: its
support is the transition, and local elliptic regularity for H_s^-1
controls the first two derivatives there. This can also be read directly
from the smooth radial ODE on that compact interval. Using (11) and
g(x)<=K_I(1+x) gives

    ||G_s^(1/2)chi_i Bv||
         <=K_I[||G_s^(1/2)chi_i h||+||h||].                  (15)

For each compact v the derivative of (12) is exactly

    F_s'[v]=<h,G_s' h>+2 Re<Bv,(G_s-c0 W)h>.

Localize the cross term with chi_o,chi_i and the bounded E_s. Use (14)
on the outer part, (13),(15) and bounded W on the interior, and (2),(3)
on the remainder. This proves

    |F_s'[v]|<=K_I(F_s[v]+J2_s[v])<=K_I' F_s[v].              (16)

The localization is essential: the full uncompensated logarithmic norm
of h is not bounded by the compensated energy near the outer seam.

## 5. Derivatives of the full local action and reconstruction

Return to physical y and write a_o=F/r, K_s=-partial_y^2+s+V_K,
B_s=K_s-2a_o partial_y+c, q_s=Q_s w=-K_s^-1 B_s w/3.
Use a_o here to avoid confusing it with the canonical warp a(Y).
For fixed compact w set u_s=q_s+w/3. The exact identities are

    q_s'=-K_s^-1 u_s, u_s'=-K_s^-1 u_s,
    ||Q_s||<4/5, ||u_s||<=17||w||/15.                       (17)

The prime on q_s,u_s in this section means an s derivative. The positive
reconstruction operator has V_K>=3lambda F. Therefore

    ||K_s^(1/2)q_s'||<=||u_s||/sqrt(s),
    integral F|q_s'|^2<=||u_s||^2/(3lambda s).                (18)

The actual q_s already has energy <10N_s/9. Thus q_s', its physical
conformal field p_s'=q_s'/r and their first derivatives satisfy the
same types of weighted estimates, with finite constants times N_s.

The derivative of the positive fourth-order conformal square can be
controlled without estimating a fourth derivative of K_s^-1. Define

    J_s=-partial_y^2-2a_o partial_y+s-(a_o)_y-a_o^2,
    C J2_s[p_s+v/3]=||F^-1/2 J_s u_s||^2.

Equation (17) gives

    partial_s(J_s u_s)=(K_s-J_s)K_s^-1 u_s,
    K_s-J_s=2a_o partial_y+F D0,
    D0=-H_R/6+2F/r^2.                                     (19)

D0 and F/r^2 are bounded. The operator F^-1/2(K_s-J_s) is bounded
from the K_s energy space into H, with bound at most
2sqrt(5lambda)+||D0||/sqrt(3lambda). Thus the derivative of
F^-1/2 J_s u_s has H norm <=K_I||w||. Cauchy--Schwarz controls the
derivative of its squared norm by K_I(N_s+C J2_s[p_s+v/3]).
Also |partial_s J2_s[v]|<=4J2_s[v]/s by (2).

For completeness, the other local terms have no hidden unbounded
frequency derivative. The exact classical formula is

    G S_C=3<q_s,K_s q_s>-2 integral F J w^2/r^2
          -integral F[(2T_C+rR_r/2)q_s^2+2T_C q_s w].

Differentiating its first term uses
<q_s',K_s q_s>=-<u_s,q_s>. The remaining coefficients are bounded
apart from their displayed F weight; (18) and the retained F-weighted
q_s,w estimates bound |partial_s S_C| by K_I N_s.

Use the exact unrestricted Weyl identity, equation (8) of the
[global form proof](global-thermal-curvature-free-form.md), before its
Young inequalities. Together with Qsc it reads

    a_sc,s=S_C+F_s+C J2_s[v]/(200pi)
                   +C J2_s[p_s+v/3]/(20pi)+R_loc,s.          (20)

R_loc contains the unchanged dT D, f L, J|grad p|^2, Wc v womega,
womega^2 and Z terms. Its absolute value and its s derivative obey

    |R_loc,s|+|partial_s R_loc,s|
       <=K_I[N_s+C J2_s[v]+C J2_s[p_s+v/3]].                 (21)

Here are the estimates needed for the derivative statement. Replace
p by p_s' in each linear coefficient D or in womega=-r p_y/B^2;
(18) controls their weighted squares exactly as for p. Explicit time
derivatives of gradient pairings produce p^2,pv,v^2, bounded by the
same first-derivative norms divided by inf(I)>0. The dT derivative
is Bv and is controlled by (2); the f derivative is (19). Pair each
of these L2 terms with its original partner by Cauchy--Schwarz and
Young. The background bounds |Rhat|<=13lambda, |U|<=22lambda(r/B)^2,
V<=5lambda(r/B)^2, omega<= (r/B)^2 and F/r^2<=5lambda are unchanged.
They are exactly the weights used in the preceding residual table;
no derivative of the outer logarithm or a new endpoint trace occurs.
This proves (21) for every term of that displayed identity, not just
for its favorable squares. The same estimates give |S_C|<=K_I N_s.

Equations (4),(20),(21) show F_s<=K_I a_sc,s. Combining (16)-(21)
therefore gives the full comparison-form derivative bound

    |partial_s a_sc,s[w]|<=K_I a_sc,s[w].                    (22)

The bounds are obtained on compact w with its global reconstructed q.
The s derivative of a decaying horizon tail may contain log(rho);
(17)-(19) control it by the resolvent estimates. It has not been
assumed to be a smooth periodic Euclidean harmonic or independently cut.

## 6. Common closed form space and norm-continuous strong feedback

Gronwall applied to (22), followed by integration and polarization,
gives locally equivalent comparison norms and

    |a_sc,s[w]-a_sc,t[w]|<=K_I |s-t| a_sc,s0[w].              (23)

The model closures are completions of the SAME compact gravitational
core. Since they embed continuously in physical H and their norms are
equivalent, their domains coincide as subsets of H. Let E denote this
common space locally, using a_sc,s0 as its fixed inner product. Covering
any compact subinterval by overlapping neighborhoods identifies the
same domain across the entire sufficiently high real ray; the norm
equivalence constants need only be local.

Equation (7), J2_s0<=K a_sc,s0 and (23) make the actual forms locally
Lipschitz as bounded maps E->E*. Equation (4) makes them uniformly
coercive on E on each compact parameter interval. If A_s denotes that
form map, its bounded inverse satisfies

    A_s^-1-A_t^-1=A_s^-1(A_t-A_s)A_t^-1.                    (24)

This proves local operator-norm Lipschitz continuity E*->E, and also
of the physical inverse H->H through the fixed embeddings. This is the
usual bounded-form resolvent construction, consistent with
[Teschl, Theorem 6.25](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe2.pdf).
The common-domain estimate (22), not that abstract theorem alone, was
the missing ingredient.

The mixed map from the retained exact trace row is

    D_s w=-(epsilon r^2 d_R/3)(v_yy-s v)
                   -epsilon B1(q_s)_y+M_w w+M_q q_s.

It is locally Lipschitz E->H. Indeed its direct derivative is controlled
by ||r/sqrt(F) v||=sqrt(C)||Bv||_c and (2), with bounded
r sqrt(F)d_R; its q_s derivative is controlled in H1 by (17),(18).
All the displayed background coefficients are fixed, bounded where
required by the retained mixed-row proof.

The physical reconstruction R_s=(1/6)K_s^-1 F is its weighted-kernel
extension, with ||R_s||<=C_R=1/(18lambda). On compact sources the
resolvent identity gives

    R_s-R_t=-(s-t)K_s^-1 R_t.                              (25)

It extends to all H by boundedness. The same identity with one derivative
shows local continuity of rho partial_y R_s, rho=ell y/(ell+y), since
||rho partial_y K_s^-1||<=ell/sqrt(s). Thus the full coefficient map
C_s=A1/6+(B1 partial_y+B0)R_s is locally Lipschitz H->H.

In the fixed E representation put b_s=-(1/(2G))R_s^*D_s:E->H.
Let b_s^*:H->E* be its dual map. The varying-energy adjoint is

    K_cross,s^dagger=A_s^-1 b_s^*.

Equations (24),(25) and the mixed-map bound make it locally Lipschitz
H->E. Consequently the ACTUAL strong remainder

    V_s=W_seed+C_s-epsilon^-1 P_seed R_s
                          +epsilon^-1 D_s K_cross,s^dagger   (26)

is locally operator-norm Lipschitz on H. This is the strong feedback,
not merely the positive variational Schur correction. It generally
is not self-adjoint. The previous uniform bound ||V_s||<=M_* holds on
the whole ray, and its curvature graph domain is unchanged.
Reciprocity of the real HH spectral form, real local coefficients and
real closures make (26) conjugation preserving for every real s here.

## 7. Fixed-phase zeros and the exact physical scope

All hypotheses of section 4 of the boundary-phase theorem now hold.
Its normalized maximal homogeneous vector is

    h_gamma=(I+Q_gamma V_gamma)^-1 b_gamma,
    b_gamma=gamma sqrt(y)K_(ik)(gamma y)/||sqrt(y)K_(ik)(y)||.

Here Q_gamma is the PURE compatible inverse adjoint used in that theorem,
not the metric reconstruction Q_s above. Its local norm continuity
follows from the explicit right-inverse formula already proved there.
The Neumann inverse makes h_gamma H-continuous. Its differential
equation then gives continuity in the fixed maximal graph norm and
hence in the two boundary traces. Their lifted phase is

    beta_gamma=arg Gamma(ik)-k log(gamma y0/2)+O(gamma^-2).

At a fixed real phase theta, the trace equation is
cos(theta)A_c+sin(theta)A_s=0. Continuity and the displayed error give
opposite signs on either side of every sufficiently high pure root,
at logarithmic distances O((gamma_j^0)^-2). The intermediate value
argument proves (1). The root intervals eventually separate for fixed
k; no uniform assertion as k tends to infinity is used.
At those roots the auxiliary homogeneous pair is

    U=h_gamma, w=-K_cross,gamma^dagger h_gamma,
    q=R_gamma U+Q_gamma,metric w.

The two original reduced rows hold weakly/strongly in their established
spaces with zero independent source. The homogeneous maximal kernel
and resonant cokernel each have dimension one. Removing both fast
traces instead leaves the previously proved injective minimal operator
with a range of codimension one; it does not give arbitrary-source
solvability. No projection of sources or selection of a phase is adopted.

The original Ward/Maxwell constraints, preparation data and physical
gravitational boundary/action requirements have not been imposed merely
by solving these two rows. The earlier conditional nonlinear fast-family
and literal-action obstruction remains in force with its stated scope.
The seed is off shell. These facts prevent promotion of (1) to a
physical instability, a nonlinear no-go theorem or a completed wormhole.
No excluded high-frequency branch, altered finite term or order
reduction is introduced. S167 and S168.3-S168.4 remain open.

Verification pins the global inverse evidence. Exact identities test the
transport commutator, weighted budgets, reconstruction derivative and
common-form inverse algebra. Independent positive-measure quadrature,
matrix logarithms and nonnormal coupled-row fixtures test the mechanisms
at two precisions. They do not evaluate the charged spectral kernel or
formally certify the continuum proof. No Lean run is selected.

Remaining to-do list: physical boundary/action, constraint and initial/state
admissibility of the auxiliary modes or justified exclusion; quantitative
thresholds and signed static compatibility; nonlinear existence and
all-sector stability; exterior/state, charge/return, opening, transport,
preparation, separate creation and full audit.
