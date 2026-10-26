# An explicit full Schur candidate: small form defect, unbounded graph output

2026-10-25 HST, E003/S169. Ordinary argument before frozen arithmetic.
Full stable semiclassical wormhole objective ACTIVE/INCOMPLETE. Closed audit
unchanged. This is an executed analytic feasibility test, not an executed
numerical application of the quantum inverse or a completed strong certificate.

OUTCOME: the requested ORIGINAL-SPACE coupled estimate is UNRESOLVED. One
explicit candidate P has defect at most .338 in a specified form source
space. Its original fast-output norm on that space is INFINITE. Consequently
this candidate certificate cannot close admission, even though its residual
is below one half. This is an estimator failure, not a construction obstruction.
No new physical gate closes. The result must not be reported as certification
of delta<=1/2 in the original coupled source graph.

## 1. Actual operator, domain and the candidate's first domain question

Use the fixed finite seed, action c_R=1e22, original clock, scale and state.
The physical canonical Hessian has the two rows

    H=[[-A,-C],[-C*,D]], K=C^-1 A(C*)^-1,
    S=I+KD=K L, L=D+C* A^-1 C.                         (1)

These are physical ACTION blocks. Using the common rescaled blocks in
current-global-mixed-schur-preconditioner.md gives exactly the same S and P
below; no physical row conversion or admission norm is changed.

Use the original paired first-order Dirac domain and complementary adjoint,
the mass cubic and finite log phases on each physical side, and the regular
rotational disk at the finite-temperature bolt. Take the retained fixed York
return R, including its global radial tail and off-shell chart contacts:

    (nu,T)=R(s,xi)=(s+U xi,V xi),
    U xi=xi/2+Q J/z, V xi=-xi-(QJ)'/2,
    -(Q^2 J'/z^4)'=Q(4xi/z-xi')/z^4,
    J(0)=0, (Q^2 J'/z^4)(z_h)=0.                      (2)

The endpoint equations in (2) are the original outer-fixed gauge return
and regular-bolt condition, not additional physical boundary conditions.
The actual mass/log data are not replaced by smooth reflected fields. In
particular the reflected cubic is only C2,1. The fixed seed forms and their
own closures are those of current-seed-full-hessian-realization.md:

    R* H R=[[-a,b],[b*,ell]], a=A,
    a>0, ell>0, |b(s,xi)|^2<=q a(s,s)ell(xi,xi), q=.338. (3)

The full ell contains the exact ultraviolet logarithm, complete angular
and thermal sums, common-restored coincidence terms and finite returns.
It is NOT a local principal model, a finite-mode matrix or a separately
regulated sum. ell^-1 below means the Riesz inverse of this FULL positive
closed form. It has not been numerically evaluated. Its existence here
uses the retained fixed-seed realization, not a new smooth-boundary theorem
or transfer to a nonlinear neighborhood.

The first question in defining the candidate is which space its inverses
act on. Let V_ell be the form domain with norm ell(xi,xi)^(1/2), and let
V_ell* be its continuous dual. Let X_F=V(V_ell), with the transported form
norm, and define J_F=V* K^-1 on the common strong core. Define Y_F by
transport from V_ell*, with norm ||y||Y_F=||J_F y||V_ell*; equivalently take
the completion of that core under this norm. The inverse transports are
understood on the original fixed-gauge ranges. They are available from
the retained York realization and mixed inverse. This is a FORM completion,
not identification with the original strong mass/log source graph.

Define the single candidate

    P=V ell^-1 J_F : Y_F -> X_F.                      (4)

This is a precise infinite-dimensional operator, although evaluating its
full ell^-1 is still a nonlocal solve. It does not conceal S^-1 inside P.
It is not an implemented numerical preconditioner. On original strong data
one must separately prove the strong range/phase bounds of (4). Those
bounds are the FIRST UNSUPPORTED STEP for using (4) as the requested P.
The form construction does not supply them. Section 3 addresses, and
rejects, the tempting energy-to-output shortcut before using any numbers.

All original mass cutoff transitions and the regular-bolt returns remain
inside (1)-(4). The exact reference lift still has

    n_m=-2 beta Q J_m/z, b_m=beta z^3/Q+beta Q'J_m,
    J_m=int_0^z t^3/Q^2,
    f_cut=2xz chi_m-2F n_m D_z chi_m,
    alpha_cut=n_m D_z chi_m.

The uncut reference rows cancel; chi_m is one through .1 and zero by .2.
The cut lift has zero reference bolt displacement and is not a global
fixed-temperature zero mode. Neither the scalar quartic nor the exact
implicit-root correlated slow pair is projected out of (1). The scalar
quartic original row, 20 beta Fstar c4, remains approximately 3.333333294 c4.
The full matching rows must be recovered from a strong solution; the form
norm does not certify independent boundary-phase outputs.

## 2. Complete coupled residual in the declared form spaces

Block multiplication on the actual common core gives

    b=-AU-CV,
    V* L V=ell+b* a^-1 b=:ell+B,
    J_F S V=ell+B,
    E=I-SP=-J_F^-1 B ell^-1 J_F.                    (5)

These identities retain BOTH equations. In particular the cancellation
in the Schur logarithmic principal part has not been undone. Put
T_b=a^-1/2 b ell^-1/2. Equation (3) means ||T_b||^2<=q. Therefore

    ||B ell^-1||V_ell*->V_ell*<=q,
    ||E||Y_F->Y_F<=.338<1/2,
    ||P||Y_F->X_F=1, ||(I-E)^-1||<=1/.662.          (6)

The spectral square roots merely express the Riesz-form proof: boundedness
of b and coercivity of a,ell are its hypotheses, already specified in (3).
They do not require a smooth even double. Completion extends (5)-(6).
No additional global pseudodifferential theorem is invoked. The actual
common regulator is removed in the retained forms before this algebra;
raw-cutoff norm convergence is neither asserted nor needed.

Equation (6) is a complete coupled FORM estimate, not new evidence that
the original strong graph contraction is .338. The original norm is kept
in the output test below. Renaming Y_F as the physical source norm would
be precisely the invalid normalization shortcut this test excludes.

## 3. Executed actual-domain output test: M_graph is infinite

Choose a nonzero smooth envelope psi supported in .002<z<.003, x=log z,
and the genuine TT coordinate xi_k=z^sigma psi(x)exp(ikx), sigma=4.005.
Use (2), not the approximation of (2), to define T_k=V xi_k. The direction
vanishes near both ends before the gauge return; the exact return retains
the permitted tail and regular bolt. Smoothness and integration by parts
make that tail smaller than every fixed inverse power of k as k tends to
infinity at the FIXED positive temperature. No uniform-in-temperature
claim or finite-frequency onset is made. These are directions inside the
actual domain, not replacement boundary conditions or a finite mode list.

On the packet, (2) gives principal values T_k=-(3/2)xi_k+O(k^-1),
U xi_k=xi_k/2+O(k^-1). Thus nu=-T/3 to leading order. Choose the exact
source y_k=J_F^-1 ell xi_k. Then

    P y_k=T_k, ||y_k||Y_F^2=ell(xi_k,xi_k).          (7)

The known full local logarithmic symbol, in the ORIGINAL action units, is

    H_Q,TT ~ -K0 (Q/beta)^2 z^-4 D_z^4 log|D_z|,
    K0=G lambda beta/(15pi)>0,
    W=-N pi B^2/(G gamma).                          (8)

This uses the retained interior transfer in current-common-regulator-
strong-limit-test.md. Moller Theorem 1 and Lemma 6.1 concern a closed smooth
spin manifold and stable kernel. They are NOT applied to the rough global
domain. The compact-interior singular coefficient transfers through the
original heat kernel; separated boundary/bolt returns and local finite
contacts have lower logarithmic order there. All are retained in ell.
See https://arxiv.org/html/0705.3857v3, Theorem 1, for the checked scope.

The original normal source symbol has leading K0 k^5 log k, and its H3
norm adds three derivatives. Equation (8) and dz=z dx consequently give

    ||G_graph T_k||H3
       >=(3/2) K0 k^8 log k ||psi||L2 (1+o(1)),
    ||y_k||Y_F
       =(3/2) sqrt(|W|K0) k^2 sqrt(log k)
           [int (Q/beta)^2 z^(2sigma-3)|psi|^2 dx]^(1/2)(1+o(1)). (9)

Fourth-order local terms and the scalar/mixed returns are lower than the
fourth logarithm in the second line. The first line uses the unchanged
original .05/.1 localization, equal to one on this packet. Extra mass,
quartic and slow finite-part subtractions do not cancel its interior
high-frequency term. The first line can be read solely as the lower
bound from the third derivative on the packet interval. If an induced
phase already prevents membership in the original global graph, the
claimed output map fails earlier; no membership is inferred from the
form construction. Where original phase extraction makes that graph
finite, the interior lower bound still applies. The source y_k is
manufactured using the FULL ell,
not an evaluation or estimate of any actual forcing column.

Hence M_graph=||G_graph P||Y_F->H3 is infinite. Quantitatively,

    liminf_k ||G_graph P y_k||/(||y_k||Y_F k^6 sqrt(log k))
      >=sqrt[K0/(|W| q_hi^2 .003^(5.01))] >1e-50,  (10)

in the inherited G=1 action units, for either retained multiplicity N=1,2.
For the last conservative bound use q_hi=1.001, gamma>5e-34, and the
actual parameter enclosure for B^2=C beta/(G lambda). The checker encloses
that expression at 512 bits. It also verifies the whole-interval Q/beta
upper guard using the actual polynomial including tau. Scaling units
changes the finite coefficient, never the infinity or derivative deficit.

This is an asymptotic whole-operator test. It does not provide a certified
finite k threshold. Scaling each direction into any positive amplitude
ball leaves the norm ratio unchanged. Even a perfect residual inverse in
Y_F cannot produce a bounded full graph-output map on this source space.
This failure is not an unknown finite constant or merely a large bound.
It is also not a failure of the original stronger-data construction.

## 4. Bottleneck reset and one specified alternative

New evidence toward the complete estimate is limited but definite: (4)
is an explicit full-operator candidate, (5) evaluates its coupled defect
on its true spaces, and (7)-(10) disprove the output comparison required
to turn that certificate into original admission. The earlier raw-cutoff
failure only tested UV matching. The new candidate keeps the full UV term
and still fails the requested certificate for a different reason. No
signed center, actual mass column or original strong defect was computed.

Required gain: SIX radial derivatives and a square root of the logarithm
relative to the form data in (9), plus the independent phase and physical
bolt maps. Neither the small number .338 nor another scalar contact
estimate supplies this gain. The original stronger source norm can supply
derivatives only if its complete finite-domain quantum/source maps are
actually bounded; that remains unknown.

Compare the specified iterated candidate P_N=P sum_(j=0)^(N-1) E^j.
Its exact defect is E^N. For a specified original output O_i, if the FULL
composite M_i,2=||O_i P E^2||Y_F->Z_i is finite and bounded, then

    ||O_i(S^-1-P_N)y|| <= M_i,2 q^(N-2)||y||Y_F/(1-q), N>=2. (11)

On compact physical charts the exact York mixed row has order at most two,
a has order four, and ell has order four with its log. Thus B ell^-1 has
order -4/log. Two residual powers gain eight derivatives and two inverse
logs, enough for the LOCAL six-derivative deficit in (9). This explicitly
tests the available local cancellation with the FULL operator's principal
orders. It is not a proof of the global boundary/phase composite (11).
The mass/log endpoint, generated moments and regular-bolt derivative maps
must be estimated before declaring M_i,2 finite. Five physical areal jets
need a separate higher composite; do not infer them from this graph test.

The decisive alternative test is (11) for the original graph and every
independent phase row, with both original equations retained. Fix N=32
for its quantitative comparison. If signed centers leave the FULL remaining
affine and feedback margins, the conditional product ceilings are

    M_graph,2 U_aff <=5.34235e-12*.662/.338^30,
    M_graph,2 U_fb  < .0012069*.662/.338^30,          (12)

where U_aff and U_fb are the actual Y_F input norms. They are UNKNOWN
in these physical action units; the recorded rescaled energy bound is
not inserted without conversion. Positive trial debits reduce (12).
For a mass column with computed trial size k_trial<5e-5, replace the
right numerator by (5e-5-k_trial)*.662. No zero trial is assumed.
The same composite certificate must include all independent boundary
inputs. A graph defect bound for P_N additionally needs the original
source-graph composite and its input embedding. Neither is supplied by (6).

Equation (12) quantifies a specified sharper alternative, not an attained
constant or new norm convention. Local derivative counting favors testing
the twice-residual composite over attempting the disproved M_graph.
No claim is made that the global target is realistically attained until
the common-restored differentiated kernel and phase bounds are computed.

## 5. Complete budget: no missing entry assigned zero

The mass-lift conversion remains 2133<||d_chi||H3<2174. The known LINEAR
|m|<1.627532e-15 spends at most 3.538254568e-12; the known distant shape
piece is <1.24e-20, leaving >6.4617454e-12 for the original near shape.
These facts do not bound the nonlinear mass or the independent mass column.

Actual kappa_m is UNKNOWN. If the sufficient post-elimination kappa_m<=5e-5
were established, its debits would be 2.14765e-12 affine, .0024731 feedback
and .0195239 energy. Every other entry must still fit JOINTLY:

    Y_other+4.084e-22 A_o+Y_phase <=5.34235e-12,
    eta_other+.0195239 gamma_E+3.712e-12 A_o J_E+eta_phase <.0012069.

All displayed unknown maps remain unknown. The quartic, correlated slow,
fast and interface inputs have not been solved. Individual zero-other-cost
ceilings A_o<=1.30811704211557e10, A_o J_E<3.25134698275862e8 and
gamma_E<.061816542801387 are conditional allowances, not observed entries.

The retained two-fast/admitted-family screen is
.028913+1e-11(1.0706e10)=.135973<.1599164, reserve .0239434.
It gives no independent boundary-phase derivative or uniform nonlinear
admission. The original missing comparisons and conditions remain

    E_center+1e-8 alpha_embed<=1e-11,
    R_X alpha_embed*.028913+Z0_other<=.25,
    R_X alpha_embed^2*1.0706e10+Z2_other<=5e7,
    Y0+1e-8 Z0+5e-17 Z2<=1e-8, Z0+1e-8 Z2<1.

No nonlinear admission/contraction, exact compatibility, physical matching,
semiclassical validity/noise or all-sector stability is established.
No actual controlled wormhole follows. An estimator failure is not a
nonexistence theorem and does not justify changing construction here.

DECISION: retain the construction; reject the unbounded form-to-graph
output certificate for (4). Sharpen the METHOD to the full composite (11)
or the original strong-source estimate, retaining exact UV and all phases.
Do not pursue further weak residual improvements as feasibility progress.

ONE NEXT CALCULATION, PROPOSED NOT EXECUTED: bound the COMPLETE
O_graph/phase P E^2 composite on the actual paired domain, including
common-restored contacts, mass/quartic/slow rows and finite-bolt returns,
and evaluate its input products against (12). Its global domain and
numerical constants are the first missing argument. The local gain test
and the failure of M_graph ARE executed; the repaired global bound is not.

Remaining to-do list: complete original-space coupled residual and output
comparisons; actual signed columns and joint budgets; nonlinear admission/
contraction; exact compatibility; physical matching; validity/noise;
all-sector stability; opening, transport, preparation, separate creation
and the full stable semiclassical wormhole objective.
