# Explicit charged metric source and its horizon moment

2026-10-20, E003/S168. The exact curvature row admits a useful local
rearrangement: absorb its curvature-proportional coupling into the scalar
potential, leaving a first-spatial-derivative metric source. On the actual
seed this source is a bounded map H1(dy)^2 -> L2(dy). Its horizon moment
has explicit L2 adjoint kernels and is nontrivial on unconstrained compact
metric profiles. Its value on constrained physical solutions remains open.

This continues [the physical curvature and range theorem](physical-curvature-source-and-boundary-range.md)
and [the charged conformal response](charged-conformal-response-and-constraints.md).
The action, finite counterterms, magnetic flux and temperature are unchanged.
The seed is still off shell. No stable wormhole or charged growing mode is
constructed, and no physical metric boundary prescription is selected.

## 1. Exact local metric coupling

On the monotone open exterior use a polar-areal representative

    ds^2=-exp(2 alpha) F dt^2+F^-1 dr^2+r^2 dOmega^2,
    alpha=v a(t,r), F=F0(r)+v b(t,r).

Here v is a variation parameter, distinct from epsilon=G/(60pi).
Below F denotes F0, and primes denote r derivatives. Locally the remaining
time-coordinate freedom removes the mixed metric component. This does not
construct a gauge regular at either endpoint. A residual time change has
zero variation of static scalar invariants and does not change this row.

Let R2 be the two-dimensional orbit scalar, H=Box_2 r, and put

    Wc=-F''+2F'/r+2(1-F)/r^2,
    j=-F''/4+(F-1)/(2r^2),
    R=-F''-4F'/r-2(F-1)/r^2.

On this background C^2=Wc^2/3 and J^a_b=diag(j,j,-j,-j), where
J_ab=Ric_ab-R g_ab/4. The general warped-product curvature identities
are given in [Maeda--Nozawa, Appendix A, (A2)-(A5)](https://arxiv.org/pdf/0709.1199v4).
Only their geometric identities at n=4 and spherical sectional curvature
one are used; their Gauss--Bonnet field equations are not adopted.

For the above variation, direct differentiation gives

    delta H=b'+F a',
    delta R2=-b''-3F'a'-2F a''-b_tt/F^2,
    Psi=delta R=delta R2-2b/r^2-4 delta H/r,
    delta C^2=(2/3)Wc(delta R2+2 delta H/r-2b/r^2),
    delta J^2=2j(delta R2+2b/r^2),
    K_R=delta(Box)R=b R''+(b'+F a'+2b/r)R'.             (1)

For example the exact orbit scalar before linearization is

    R2=-F_rr-3F_r alpha_r-2F(alpha_rr+alpha_r^2)
       -exp(-2alpha)[F_tt/F^2-alpha_t F_t/F^2-2F_t^2/F^3].

It follows by writing the orbit metric as -A^2dt^2+B^2dr^2 and using
R2=-2/(AB)[partial_r(A_r/B)-partial_t(B_t/A)]. For delta J^2, the
angular mixed Ricci component is K=(1-F-rH)/r^2, so
delta J^2=2j(Psi-4 delta K). This also accounts for the inverse-metric
variation in the invariant. The b_tt term is retained in (1).

The preceding trace coupling is Bmet=K_R+7 delta C^2/12+11 delta J^2/6.
Eliminate delta R2 using Psi, without any field equation, to get EXACTLY

    Bmet=c_R Psi+d_R(b'+F a')+e_R b,
    c_R=7Wc/18+11j/3,
    d_R=7Wc/(3r)+44j/(3r)+R',
    e_R=44j/(3r^2)+R''+2R'/r.                          (2)

This identity holds for any static F in this chart. Fixed magnetic flux
has delta_A M=0, as in the preceding trace derivation. In a general gauge
Psi=delta R-R_r delta r and f_A=(DT[h]-T_r delta r)/epsilon; the off-shell
contact remains. Equation (2) refers to their polar-areal representative.

## 2. Rearranged wave operator and global source bound

With E=1/epsilon+11lambda/3 and Z=R+12lambda, the exact trace becomes

    (Box+E-11Z/36+c_R)Psi+d_R(b'+F a')+e_R b=f_A.

Set y=integral_r^infinity ds/F(s), u=r Psi, and beta=b/F. Thus
partial_y=-F partial_r. Then

    u_tt+Atilde u=g_mix-rF f_A,
    Atilde=-partial_y^2+Vtilde,
    Vtilde=F(F'/r-E+11Z/36-c_R),
    g_mix=-A1(a_y+beta_y)+A0 beta,
    A1=rF d_R, A0=rF(d_R F'+e_R F).                    (3)

Only the curvature multiple c_R Psi has moved to the diagonal. No
order reduction, term deletion, time integration or altered coefficient
has occurred. The physical curvature and trace-source weights stay
L2(r^2dy) and L2(r^2F^2dy). Initial u0,u1 still encode the required
higher metric jets; this rearrangement does not free those data.

For the retained rational seed

    F=lambda r^2+1-2m/r+Qhat/r^2+c3/r^3+c4/r^4,
    c3=-2d B^3, c4=d B^4, Qhat=B^2(1+3lambda B^2+d),

the positive-temperature outer horizon r_h>B is simple and m>0.
The three coefficients in (2) can be evaluated as polynomials. Set
z=B/r, x=lambda B^2, mu=m/B, q=Qhat/B^2 and

    f(z)=z^2 F(B/z)=x+z^2-2mu z^3+q z^4-2d z^5+d z^6,
    D(z)=28mu-(128/3)q z+(440/3)d z^2-100d z^3,
    H(z)=-(44/3)q+(460/3)d z-246d z^2.

Then d_R=B^-3 z^4 D and e_R=B^-4 z^6 H, while

    c_R=B^-2[(14/3)mu z^3-(25/3)q z^4
                            +(305/9)d z^5-(169/6)d z^6],
    A1=B^-2 z f D,
    A0=B^-3 f[D(2f-z f')+z H f].                      (4)

The prime in (4) is a z derivative. On 0<z<z_h=B/r_h<1 these are
bounded. If |P|_1 is the sum of absolute polynomial coefficients, explicit
valid upper bounds for C1=sup|A1| and C0=sup|A0| are

    C1 <= B^-2 |f|_1 |D|_1,
    C0 <= B^-3 |f|_1[|D|_1(2|f|_1+|z f'|_1)+|H|_1|f|_1].

Consequently every a,beta in H1(0,infinity;dy) satisfies

    ||g_mix||_2 <= C1(||a_y||_2+||beta_y||_2)+C0||beta||_2. (5)

This proves source-space membership for this stated class, without
assuming any coefficient is small. It is not a complete physical metric
norm and does not show all horizon-regular metric data lie in this class.

Writing h0=28lambda^2 m>0, the endpoint limits are

    A1=h0 y+O(y^2), A0=2h0+O(y), y->0,
    A1,A0=O(exp(-2kappa y)), y->infinity.               (6)

In particular the zeroth-order source coefficient does not vanish at
the outer boundary. Polynomial decay of the individual curvature
coefficients does not justify discarding their weighted source.
The added potential -F c_R tends to zero at both ends and is globally
bounded. Thus Vtilde has the SAME inverse-square coefficient and bounded
remainder limit (2lambda-E)/3 as before, and an integrable horizon tail.
The preceding scalar endpoint and range proofs apply to Atilde. Every
scalar self-adjoint phase is still unbounded below. This remains a
diagonal conclusion, not a coupled metric stability theorem.

## 3. The metric contribution as an explicit bounded functional

For Re gamma>0 construct wtilde_gamma by (Atilde+gamma^2)wtilde=0 and
wtilde exp(gamma y)->1 at the horizon. It is the NEW horizon solution;
the old w_gamma cannot be reused after changing the diagonal. Write
Mtilde_gamma(g)=integral wtilde_gamma g dy. The previous theorem gives
Ran(Atilde|D_min+gamma^2)=ker Mtilde_gamma and zero kernel, where D_min
has both fast outer traces zero. A bounded potential change preserves
this domain and the leading trace definitions.

Integration by parts gives, with w=wtilde_gamma,

    J_gamma[a,beta]=Mtilde_gamma(g_mix)
       =integral (K_a a+K_beta beta)dy-[w A1(a+beta)]_0^infinity,
    K_a=(w A1)_y, K_beta=(w A1)_y+w A0.                (7)

The integral is bilinear for complex gamma. For compact smooth profiles
the boundary term vanishes. It also vanishes for H1 profiles: they are
bounded, w A1=O(y^(3/2)) at zero and decays exponentially at infinity.

To check the kernels, let
S=A_c cos(k log(y/y0))+A_s sin(k log(y/y0)),
T=-A_c sin(k log(y/y0))+A_s cos(k log(y/y0)). The controlled outer
expansion w=sqrt(y)[S+O(y^2)] and (6) give

    K_a=h0 sqrt(y)[(3/2)S+kT]+O(y^(3/2)),
    K_beta=h0 sqrt(y)[(7/2)S+kT]+O(y^(3/2)).           (8)

At the horizon both kernels are O(exp(-(Re gamma+2kappa)y)). They
therefore belong to L2 globally. In particular (7) extends uniquely
from H1 pairs to a bounded functional on L2(dy)^2 with

    |J_gamma[a,beta]| <= (||K_a||_2^2+||K_beta||_2^2)^(1/2)
                                         (||a||_2^2+||beta||_2^2)^(1/2).

For a general L2 pair this defines the moment, not necessarily an L2
source g_mix. The bounds are locally uniform on compact subsets of
Re gamma>0, by the preceding Volterra bounds and smooth coefficient
dependence. No uniform bound up to the imaginary axis is claimed.

For real gamma>0 all quantities are real and K_a is not identically
zero on any sufficiently small outer interval. Otherwise w A1 would
be constant there, its zero endpoint limit would force it to vanish,
and A1=h0 y+O(y^2) would force the nonzero homogeneous w to vanish.
Choose a smooth nonnegative compact chi supported where K_a is nonzero.
Taking a=chi K_a, beta=0 gives J=integral chi K_a^2>0. Changing the
sign gives a negative value, and both profiles can be arbitrarily small.
Thus geometry alone does not make the isolated metric moment vanish.
These profiles have not been shown to satisfy the other metric row or
the constraints and are not physical instability counterexamples.

## 4. What this condition does and does not decide

For transformed curvature initial data u0,u1 the complete requirement is

    Mtilde_gamma[gamma u0+u1+g_mix_hat-rF f_A_hat]=0.    (9)

The remaining source contains no time derivatives, so its transform
uses the metric profiles directly. The initial curvature terms and the
off-shell contact in f_A must still be retained. The old and new
functionals are equivalent only as part of their respective full
operator equations; their values on arbitrary isolated sources differ.

Crucially, (9) is the trace equation projected onto its range condition,
not an independent additional physical equation after that trace has
already been solved. For a compact frequency-domain metric test profile,
form Psi and f_A from that SAME profile, with zero transform initial
contacts. Then U=r Psi is compact and (3), followed by Green's identity,
gives Mtilde_gamma(rF f_A)=J_gamma. A nonzero isolated J cannot prove
incompatibility when the corresponding f_A is generated as well.

The outstanding question is whether the full retarded shape equation,
physical initial constraints and admitted boundary data enforce the
required total moment for the prescribed physical source. Neither an
arbitrary independent scalar forcing nor the compact examples above
answer it. Neutral-current decoupling and Ward recovery are available;
their use does not supply that missing solution. The off-shell seed must
also be transferred to an exact background before physical stability.

## Verification scope

The checker compares (1)-(2) against independent four-dimensional
Christoffel/Riemann contractions and exact first variations at rational
metric jets, including time-dependent and charged backgrounds. Exact
symbolic checks verify (3)-(6) and the projected Green identity. Numerical
quadrature checks (7) on actual-form seed coefficients with explicitly
chosen test weights. Those weights are not horizon solutions, and those
fixtures are not a benchmark spectrum or evolution. The global functional
and range arguments remain ordinary mathematics; no Lean is selected.

Remaining to-do list: evaluate the total horizon moment with the retarded
shape row and admissible metric data; constrained stability transfer or
justified exclusion; S167 signed source, cold/higher response, nonlinear
solve and compatibility; exterior/state, charge/return, opening, full-sector
stability, transport, preparation, separate creation and final audit.
