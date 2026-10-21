# Complete shape cutoff limit and the conditional fast-family obstruction

2026-10-21, E003/S168. Ordinary argument before focused verification.
The full stable-wormhole goal remains incomplete.

The original pure-shape quadratic action has a finite translated-
logarithmic weak cutoff limit on the full retained auxiliary energy
space. Its canonical and local logarithmic terms must be combined
before taking that limit. Consequently it cannot provide the persistent
oscillation required by the coupled nonlinear trace expansion. This
closes that expansion's conditional boundary-action obstruction.

This does not exclude every nonlinear solution or physical fast mode.
The theorem concerns the stated second-order expansion and convergence
of its literal quadratic cutoff coefficient, with the original action,
state and paired domain. Physical source and boundary admissibility,
causal constraints, nonlinear existence and stability remain unresolved.

## 1. Original action and shape domain

Use the exact allocation

    I=S_classical+W4+Gamma_hat,LLL+Gamma_hat,rest.

The reference functional in I0+Gamma_rest cancels algebraically and is
not an additional term. The declared LLL is included once. If Qcan in
the [canonical comparison](smooth-local-quantum-comparison.md) includes
that continuous lower-order allocation, subtract it there before adding
it separately below; this leaves the comparison's conclusion unchanged.
Keep the prepared period, physical clock, flux, common subtraction,
finite coefficients, GHY, Einstein counterterms, loop counterterms and
the converted literal Weyl boundary functional. The cut is fixed in
the background coordinate before varying the metric.

Fix one positive Matsubara frequency above the retained auxiliary and
reconstruction thresholds. A pure shape tangent is (s,v), with
s=Q_gamma w/r, v=w/r, w in E_gamma. Near the outer end put Z=-log Y,
where Y is dimensionless inward canonical proper distance. Then

    s=Y^(3/2)S(Z), v=Y^(3/2)V(Z), S,V in H2((Z0,infinity)). (1)

This is the domain result already proved in the
[energy-shape trace argument](energy-shape-nonlinear-trace-forcing.md).
Its fixed smooth coordinate changes preserve (1). Begin with fields
smooth for Y>0; the estimates extend their quadratic cutoff distributions.
They do not differentiate arbitrary finite-amplitude rough metrics.

For either H=S or V define

    P0 H=H, P1 H=3H/2-H', P2 H=H''-2H'+3H/4.

Then h^(j)(Y)=Y^(3/2-j)Pj H for j=0,1,2. All Pj H belong to L2(dZ).
It follows directly that every bulk term

    b(Y) Y^(i+j-4) h^(i) k^(j) dY, 0<=i,j<=2,          (2)

is absolutely integrable when b is bounded: it becomes a bounded
coefficient times Pi H Pj K dZ. Every boundary term with coefficient
b(Y)Y^(i+j-3) is in L1(dZ) and hence tends to zero in translated local
L1. The same holds with a factor Y^epsilon log^m Y, epsilon>0.
Neither claim assumes pointwise decay of the second derivatives.

## 2. Extract the only critical logarithmic terms

Write the canonical metric as ell_o^2[dY^2+a(Y)^2 dtau^2],
C=B^2/ell_o^2, c=C/(30pi), beta=a'/a. The retained smooth seed has

    beta=O(Y), beta'=O(1),
    omega=log(r/B)=-log Y+omega0+O(Y^2),
    omega'= -1/Y+O(Y), omega''=1/Y^2+O(1).             (3)

All stated remainders have the differentiated bounds required here.
They follow from the retained rational areal seed and its proper-
coordinate expansion. Primes in this section mean inward Y derivatives.
Normalize by the actual positive time-mode norm.

The exact converted bulk density is the one in
[the finite-cut Weyl functional](stationary-exterior-boundary-weyl.md),
not the literal boundary heat path added to an already integrated bulk.
Its Euclidean Hessian contains the single critical logarithmic term

    QW,log,bulk=-c integral_d a omega (Delta_bar v)^2 dY.       (4)

All other explicit omega terms contain background canonical curvature,
sphere curvature or magnetic coefficients and at most two derivatives
on the shapes. They gain at least two Y powers relative to (2).
The exact completed-square bulk identity in
[the unrestricted Weyl calculation](global-thermal-curvature-free-form.md)
also displays (4). Its ordinary 1/6 contact is retained in the remainder.

For an independent boundary extraction, use the original outward
proper derivative. In dimensionless units under the pure shape path,

    A_eta=-exp(-eta v)(beta+eta v'),
    A_dot_eta=exp(-2eta v)[beta'+eta v''-eta v'(beta+eta v')],
    T_eta=-exp(-eta v)(omega'+eta s'),
    T_dot_eta=exp(-2eta v)[omega''+eta s''-eta v'(omega'+eta s')].

Here a dot is an outward proper derivative, not a Y derivative. Put
k_s=1/C. The original polynomial is

    G_eta=(omega+eta s)[12 A_eta A_dot_eta+44 A_eta^3/7+10 k_s A_eta]
       -6 A_dot_eta T_eta-24 A_eta^2 T_eta/7-6 A_eta T_dot_eta
       -96 A_eta T_eta^2/7-18 T_eta T_dot_eta-85 T_eta^3/21. (5)

The Euclidean boundary functional is +C a exp(eta v)G_eta/(360pi).
Twice differentiating it fixes the only critical logarithmic Hessian:

    QW,log,boundary=-2c a(d) omega(d) v'(d)v''(d).             (6)

The sign includes the outward-to-inward conversion. The remaining
omega terms in (5) contain beta or beta', and gain at least two Y
powers relative to the critical boundary weights in (2). Terms without
an explicit omega satisfy those weights with no logarithmic loss.
In particular s's'', omega' v s'', and omega'' v s' are allowed:
their rescaled products are in L1(dZ). There are no third shape
derivatives in this converted static boundary representative.

These formulas retain the standard anomaly coefficients. The underlying
bulk cocycle is [Mazur--Mottola, section IV](https://arxiv.org/html/hep-th/0106151v1);
the mixed boundary heat coefficient is [Vassilevich (5.33)](https://arxiv.org/html/hep-th/0306138v3).
Neither source alone fixes this finite-cut conversion; (5) is the
previously derived conversion with all its fluxes retained.

## 3. Match the density-conjugated normal compensator

Put a0=a(0), m=sqrt(a/a0), u=mv, and

    W=m''/m+gamma^2/a^2.

Since 2m'/m=beta, the fixed-frequency identity is exact:

    m Delta_bar v=u''-Wu,  m'=O(Y), m''=O(1), W=O(1).        (7)

The normal compensator supplied by the canonical theorem is

    N_d=a0 c [integral_d log Y (u'')^2 dY
                                      +2log(d)u'(d)u''(d)]. (8)

Subtract its bulk part from (4). The difference is exactly

    -a0 c integral_d [(omega+log Y)(u'')^2
                          -2omega W u u''+omega W^2 u^2]dY. (9)

It has an ordinary finite limit. The first coefficient is bounded by
(3); u'' is L2. Also omega u is L2, and W is bounded, so the cross
term is absolutely integrable. No separate logarithmic u'' norm is
required or inferred from H2.

The difference of (6) and the boundary part of (8) is

    -2a0 c [(omega+log d)m^2 v'v''
                          +log d (u'u''-m^2 v'v'')] at d.   (10)

The first term is an L1(dZ) product with bounded coefficient. For the
second use the exact expansion

    u'u''-m^2 v'v''=mm' v v''+2mm'(v')^2
                  +(mm''+2(m')^2)v v'+m'm''v^2.             (11)

Every term on the right gains at least two Y powers relative to its
critical boundary weight, so multiplication by log Y preserves L1(dZ).
Thus (10) tends to zero in translated local L1. This matches both the
bulk and the boundary logarithms of the ORIGINAL Weyl action; (8) is
not substituted for that action and no finite coefficient is changed.

## 4. All remaining original local terms

For clarity, the following is a complete weight classification of the
quadratic local terms, not a claim that they vanish identically.

| Part | Bulk coefficient for shape jets i,j | Boundary coefficient |
| --- | --- | --- |
| Weyl after (4),(6), no explicit omega | O(Y^(i+j-4)), i,j<=2 | O(Y^(i+j-3)), i,j<=2, i+j<=3 |
| Remaining explicit omega terms | O(Y^(i+j-2) log Y) or better | O(Y^(i+j-1) log Y) or better |
| Einstein, cosmological and magnetic bulk | O(Y^(i+j-4)), i,j<=2 | GHY: O(Y^(i+j-3)), i+j<=1 |
| Einstein and loop intrinsic counterterms | no new bulk term | O(Y^-3)h k plus lower weights |
| Fixed gamma^2 corrections | at least two powers better than the first row, possible log Y | at least two powers better, at most one normal derivative |
| Fixed gamma^4 corrections | four powers better, possible log Y | none of derivative weight three |
| Declared canonical LLL | -q J1[v]/(12pi) | affine cocycle; no additional critical quadratic term |

The Weyl rows follow by literal expansion of its bulk polynomial and
(5): the only unbounded background jets are omega', omega'', with
weights -1,-2. beta has weight +1 and beta' weight 0. After the
extracted logarithms, the checker audits every monomial and its weight.
The square Delta(s+v/3)^2 is a fourth-order term covered by (2), not
discarded as a first-order error.

For the classical rows, the unintegrated spherical Einstein density has
e^(2omega) multiplying curvature and two shape derivatives, and the
cosmological potential has e^(4omega). These have weights -2 and -4.
The magnetic potential is no worse. GHY contains e^(2omega) times one
normal jet; the intrinsic boundary volume and curvature counterterms
have e^(3omega) and e^omega. The fixed loop c0,c1 terms have those same
weights. Their actual coefficients remain unchanged. This also retains
the scalar-curvature vacuum shift in lambda0.

The nonzero-frequency bulk terms are the explicit polynomial in
[the local frequency calculation](local-coupled-frequency-polynomial.md).
Its gamma^2 coefficient has at most two radial derivatives and its
gamma^4 coefficient none. In the full spherical boundary heat coefficient
two tangential derivatives leave at most one normal derivative; periodic
time integration removes tangential divergences. Its coefficients obey
the displayed better weights using (3). Thus the static polynomial is
not asserted to be an exact dynamic boundary formula.

The declared LLL has the same fixed conformal domain and state. Its
Polyakov Hessian is the retained Dirichlet form, with affine boundary
cocycle; it has an ordinary limit here. No zero-mode normalization term
is introduced. Far from the outer end, the fixed smooth seed and the
same form-domain approximation give finite local integrals. The bolt
is an interior polar point and no punctured-bolt boundary is added.

Equations (1)-(2) now prove ordinary convergence of all remaining bulk
terms and vanishing translated local L1 boundary tails. Their finite
interior values are retained. Combining (8)-(11) with the canonical
normal theorem proves

    Qshape,d[s,v] has a finite translated-logarithmic weak limit. (12)

This is the literal pure-shape Hessian extended by the stated estimates.
It is not an ordinary pointwise cutoff theorem on arbitrary energy
profiles, and is not a nonlinear differentiability theorem.

## 5. Conditional obstruction and the next physical question

Apply (12) to exactly the family in section 1 of the retained
[nonlinear trace argument](energy-shape-nonlinear-trace-forcing.md):
its first-order tangent has a nonzero fast amplitude a_f and the full
auxiliary energy-class shape; its normalized second-order tails h,j
have the stated bounded jets and controlled distributional remainders.
It is required to satisfy the original equations through second order
and to have an ordinary limit of the literal quadratic cutoff coefficient.
These are explicit hypotheses, not properties proved for all solutions.

The nonlinear trace forces the persistent vector, in the original
Lorentzian Taylor-coefficient convention,

    Rc=-(78815t^3+74391t^2+16137t+81)
                         /[168(5t+3)(11t+3)],
    Rs=-k(6655t^3+111t^2-3015t-567)
                           /[84(5t+3)(11t+3)],
    t=G lambda/(60pi)>0, k^2=1/t+17/12.                    (13)

It multiplies S0 <H^2> a_f^2, S0=lambda B^3/G, and rotates with twice
the fast phase. The exact positive identity

    -Rc-3/56=t(78815t^2+73896t+15705)
                              /[168(5t+3)(11t+3)]>0         (14)

excludes cancellation by choosing a phase. The mixed first-order
action and canonical first variation on the second-order tail already
have finite necessary cutoff limits. The pure-shape Lorentzian Taylor
coefficient is minus one half of the Euclidean Hessian in (12), hence
also has a finite weak limit. It cannot cancel (13).

A nonzero sinusoid has no translated distributional limit to a constant:
choose a compact smooth test with nonzero Fourier transform at 2k and
translate through two different phases. The vanishing weak remainder
cannot change the resulting different limits. The smooth relation
log(r/B)=-log Y+constant+O(Y^2) preserves this conclusion. Therefore
no family satisfies all the hypotheses above with a_f nonzero.

This closes the specified regular perturbative family, not the research
goal. In particular existence of a finite action at each nonzero
amplitude does not by itself justify interchanging its second amplitude
derivative with the cutoff limit. Families outside the bounded-tail
expansion, a different physically justified boundary completion, and
the zero-fast-amplitude sector are not excluded by this theorem.
No new completion or counterterm has been selected.

The next decision is physical admissibility: determine the permitted
gravitational boundary data and whether the prepared source and causal
constraints admit the zero-fast sector, or require a different justified
family. A Euclidean action obstruction is not a Lorentzian instability
theorem. Exact nonlinear existence and all-sector stability remain open.

Remaining to-do list: physical boundary and source/initial-data
compatibility, including the zero-fast sector and static signed condition;
continuous causal constraints, exact nonlinear existence and all-sector
stability; exterior/state, charge/return, opening, transport, preparation,
separate creation and full audit.
