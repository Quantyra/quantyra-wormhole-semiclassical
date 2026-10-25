# Complete fixed-seed input-energy conversion; coupled strong test remains open

2026-10-25 HST, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: no progress. Full stable-wormhole objective incomplete.
The closed feasibility audit and original physical targets are preserved.

This addresses the energy input factors in the complete return error, not
the missing strong output map. All quantum terms enter through the actual
two-sided finite-seed comparison. No physical Hessian application or signed
coupled solve is claimed. The requested quantitative feasibility milestone
remains incomplete.

## 1. Exact global inverse of the frozen York conversion

Use the current finite polynomial Q, exact positive-temperature root z_h,
and the original paired mass/log phase. In the notation of
current-physical-york-dual-assembly.md set X=QJ=z^2 y and
A=Q''-4Q'/z. Its two kinematic equations give exactly

    -y''+[6/z^2+2A/(3Q)]y=2T'/(3z^2)-8T/(3z^3),
    y(0)=y(z_h)=0,
    xi=-T-z y-z^2 y'/2,
    s=nu+T/2+z^2 y'/4-z y/2.                         (1)

These endpoint conditions express the original outer-fixed return and
finite J at the regular bolt. They are not new physical boundary data.
The potential's 1/(z_h-z) singularity has the regular solution y=O(z_h-z);
then J=z^2y/Q is finite. On the mass/log core y=O(z^2 log^k z) at zero.
At the bolt a smooth stationary T yields the regular one-sided solution
by substitution in (1). Completion uses this variational realization.

Let P=6/z^2+2A/(3Q). On the WHOLE physical interval P>=5/z^2.
Indeed Q+(2/3)z^2 A equals

    beta-3u z^2+3a z^3-(5/3)b4 z^4-2eps z^5+5eps z^6,

where a=4beta+2u+tau and b4=3beta+u+eps. It is bounded below by
beta-3u-5|eps|>0 since 3a-5b4/3>0. Keep Q>0 only on the actual interval;
no division by a separately rounded bolt value is performed.

Put l=1-log z and

    N(f)=||sqrt(l) f/z^2||2+||sqrt(l) f'/z||2
         +||sqrt(l) f''||2+||sqrt(l) Q^(3/2) f'''||2,          (2)

with measure dz on (0,z_h). This is an explicitly specified INPUT norm,
not a redefinition of the physical admission norm or fast output.
Testing (1) with l y and integrating before bounding gives

    E_y=int [l y'^2+(l P-1/(2z^2))y^2],
    E_y>=||sqrt(l)y'||2^2+(9/2)||sqrt(l)y/z||2^2,
    sqrt(E_y)<=c ||sqrt(l)T/z^2||2, c=2sqrt(17)/9.     (3)

The right pairing is -(2/3)int l(T/z^2)[y'+(2-1/l)y/z].
Cauchy--Schwarz in the two energy components proves (3). Truncation and
the displayed endpoint rates justify the logarithmic test; no unproved
boundary term is dropped.

An important cancellation occurs in the physical scalar Laplacian. Define
L_s=Qs''+(Q'-2Q/z)s'. Direct substitution, BEFORE division by Q, gives

    L_s=L_nu+(QT')'/3+(2Q'/(3z)-2Q/z^2)T
       +(Q+z^2Q''/6-2zQ'/3)y'
       +(z^3Q'''/6-2z^2Q''/3+5zQ'/3-2Q)y/z.          (4)

There is no inverse Q in these coefficients. The actual seed has
|Q|<.17, |Q'|<.3, |Q''|<2.1, |Q'''|<8.1 and |A|<3.
Also (A/Q)_-<1e-9: for z>=1e-9 use
A>=z(6a-4b4)-6u-6|eps|>0; below it Q>.16.
Equations (3)-(4) and

    s'=nu'+T'/3+2T/(3z)+(A/(6Q)+1/z^2)z^2y

therefore give, with N=N(nu)+N(T),

    ||sqrt(l)s/z^2||2<N,
    ||sqrt(l Q)s'/z||2<N,
    ||sqrt(l)L_s||2<3N,  ||sqrt(l)(Qs')'||2<4N.       (5)

For the apparent gradient denominator, keep the potential energy:
int l |A|y^2/Q <=(1.5+2e-9/(9/2)) E_y. Thus its contribution is at most
[sqrt(3(1.5+2e-9/(9/2)))/6+sqrt(.17/(9/2))]sqrt(E_y).
The remaining T/z^2 coefficient is 2sqrt(.17)/3; their sum is <1.
For (4) the corresponding coefficient is
.54+.72c+3.59c/sqrt(9/2)<3. These are global bounds, not cold asymptotics.

The TT scalar return S=nu-s depends only on T. The same bounds hold for
S with N(T), by setting nu=0 and reversing the sign. In particular BOTH
diagonal form inputs are controlled; no indefinite total-action bound
is incorrectly substituted for their positive sum.

## 2. Complete action-energy bound in the original units

Use the exact local generator in current-finite-local-action-assembly.md
and the same affine chart/frozen York congruence. Normalize N=1 temporarily
by homogeneity, not by changing an admission target. For EACH of (s,0)
and (S,T), (5) permits the following integral product majorants:

    nu0<=1, T0<=1, gradient<=1, L_nu<=3, Delta nu<=4,
    Delta T<=1, D0*z^2<=6, |delta D|_L2(l dz)<=20,
    integral |delta^2 D|/z^2<=12.

In the Weyl polynomial use seed/first/second majorants
phi=(1,1,0), L=(3,2,0), M=(1,4,0), Nphi=(1,1,2).
The seed phi entry represents |log z|<=l, not a uniform bound by one.
Seed M,Nphi,D carry z^-2; variation nu and T carry the compensating z^2
in (2). All second products are bounded by weighted L2 Cauchy--Schwarz.
There is at most one undifferentiated phi, so the single l weight suffices.
Exponential jets are (1,|k|,k^2). Multiplying these triples gives a
reproducible positive upper bound for every term of the exact generator.

Retain the shared trial vacuum interval, actual a2, exact LLL and complete
positive geodesic term. Its full sixth-moment bound is

    Pgeo/|W| <= [2/(525pi B^2)] [4C^2/q_flux] Z,
    Z=int Q[(5/8)Q^2 T'''^2+(3/2)Q T''' U+U^2]dz,
    U=(3/2)Q'T''+(1/2)Q''T',   Z<2 N(T)^2.          (6)

This includes the whole angular spectrum and every flow frequency. The
actual thermal, paired-boundary, coincidence and finite-domain differences
remain in |l-l0|<=.0022 l. The scalar diagonal is unchanged by that
comparison. No determinant is replaced by a finite list of angular modes.
The complete normalized comparison diagonal majorant is <78.000004;
using it for both diagonals gives

    [a(s)+l(xi)]/|W| < 2(78.000004) N^2/(1-.0022).

Restore the original radial energy, removing the SAME factor 4pi beta_T
per half from the action, beta_T=2pi B/(C gamma). Since
|W|/(N_halves 4pi beta_T)=CB/(8piG), the exact conversion is

    E(u)^2=[a(s)+l(xi)]/(|W| beta),
    E(u)<31 [N(nu)+N(T)].                            (7)

The evaluated upper expression is near 30.627827. This is an INPUT-to-
energy bound at the seed. It is not energy-to-graph, R_X, or an X-ball
estimate. The inverse direction is not asserted.

## 3. Exact contribution to complete-return errors

For R_s and H_s of current-short-travel-strong-domain.md, keep
eta_s=.0022000000000000032 and q_s=eta_s/(1-eta_s). Equations (7) and
the ACTUAL coupled inverse prove for EVERY input in (2), independently
of whether it is fast, a boundary phase or interior,

    E(H^-1 R_s v)<=31 eta_s N(v),
    E(H_s^-1 R_s v)<=31 q_s N(v).                    (8)

Thus for an actual source column F_j=R_s Phi_j its input factor in
the complete tail is U_j<=31 eta_s N(Phi_j), and

    error_ij<=M_i,p q_s^(Niter-p) 31 eta_s N(Phi_j).   (9)

This closes a numerical conversion for the FULL coupled energy response
to these specified columns, including scalar/TT mixing. It does not set
an unspecified post-elimination F_j equal to R_s Phi_j. Additional source
maps in the original global matching system still need their own norms.
Nor does (8) compute a signed output, M_i,p, or the graph feedback norm.

For the retained Niter=64, p=4 graph-tail allocation, a sufficient
rounded-down replacement for M_graph,4 U_feedback<=1.5062e156 is

    M_graph,4 N_feedback<=2.2085e157,                 (10)

ONLY when the complete feedback source under consideration is R_s times
an input map with the stated N_feedback norm. Use the upward guard
31 eta_s<.0682000000001 in this substitution. Neither N_feedback nor
M_graph,4 is assigned one. The original matching system may contain
additional input maps, so (10) is not a certificate for its full feedback.
The sharper established affine energy bound 1.1e-10 is unchanged.

The exact mass lift, its arbitrary previously fixed cutoff transition,
scalar quartic and exact correlated slow pair all belong to (2). Their
actual cutoff derivative norms must be inserted before quoting individual
numerical U_j; no unspecified smooth cutoff is assigned convenient jets.
For the correlated pair its jets must be evaluated using B0 e_* and
expm1, as in operators/exact_slow_phase.py, rather than separate B0 terms.
The present acceptance evaluates the universal constant (7)-(8), not those
individual column norms. This limitation is explicit.

## 4. Unresolved feasibility decision

The first missing STRONG argument remains numerical bounds for the complete
differentiated R_s kernel and contacts, composed with the actual inverse,
original phase/moment equations and finite-domain reconstruction. The
known H3->L2 interior estimate supplies no differentiated output constant.
No external smooth-boundary theorem is imported to fill that gap.

No kappa_m, signed coupled center, actual mass or complete feedback norm
was computed. At conditional kappa_m=5e-5, the remaining affine and graph
budgets are still 5.34235e-12 and .0012069, with the independent energy
debit .0195239. The conversion >.00906808164|m| in the ORIGINAL fast norm
is unchanged. The screen .135973<.1599164 retains its two-fast-input
admitted-family scope; (8) gives independent phase control only in energy.
No nonlinear neighborhood follows from these seed statements.

RESEARCH DECISION: retain the construction and the complete strong-return
method. Use (7)-(9) for its actual input-energy costs; do not insert them
for the missing strong output constants. The kinematic conversion does
not force a tiny-bolt denominator loss. Quantitative feasibility remains
undecided, with no construction obstruction established.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): bound the complete original-
output composition T_old(H_s^-1 R_s)^4, including common-restored contacts,
all angular/thermal contributions and original matching/bolt returns.
Combine its bound with (9) only for the identified source columns; the
remaining matching input maps and signed solve are still required.

Remaining to-do list: complete strong return and signed coupled response;
original physical-norm comparisons; nonlinear admission/contraction and
exact compatibility; physical matching; semiclassical validity/noise;
all-sector stability and the full stable semiclassical wormhole objective.
