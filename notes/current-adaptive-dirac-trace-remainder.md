# Complete-spectrum Dirac trace remainder; strong coupled return unresolved

2026-10-25 HST, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: PROGRESS (actual-domain inverse series). The full
stable semiclassical wormhole objective stays active and incomplete. The
closed feasibility audit, physical construction and admission targets stay.

OUTCOME: UNRESOLVED strong coupled estimate. This attempt closes a specific
continuum gap in the new propagator method: a mode-dependent truncation
has a trace-class remainder, summed over the ENTIRE angular and thermal
spectrum, with common-regulator removal and mass-field variations. Its
normalized Hessian error is <9.204e-77 in L-infinity -> L1, not in the
original H3 graph. No signed center, mass column or coupled inverse is
computed. A general conversion from this weak error to the strong graph
does not exist. The actual remainder needs its own differentiated bound.

## 1. Domain, exact representation and the first unsupported passage

Use current-global-dirac-transport-series.md on the ACTUAL paired mass/log
domain and its C2,1 identical-field double, with both original regular
bolts, all M_n^2=n(n+q_flux)/C, d_n=2(q_flux+2n), and
p=gamma(j+1/2), j in Z. Keep the exact finite Q and correlated h, gamma.
The fixed seed unitary factors conjugate the original first-order and
composed adjoint domains. There is no smooth seam replacement, new wall,
Robin condition, or finite physical mode selection. The massless LLL,
Weyl allocation, full geometric-to-mass conversion and original local
counterterms remain in their existing allocations. The conversion is the
common-cutoff identity of current-third-metric-conversion.md, not merely
a two-dimensional massless cocycle. Independent odd-mouth variations are
not covered.

In the frozen mass chart A_b=D+M exp(b), the rotated expression is

 B_b=B0+V_b, B0=partial_y-k sigma3, k=sqrt(M^2+p^2/Q),
 V_b=alpha J-M(exp(b)-1) O^T sigma3 O,
 alpha=(atan[p/(M sqrt(Q))])'/2.

The diagonal G0=B0^-1 is the full-interval transport inverse already
proved. Set K=k^(1/2) G0 k^(1/2), S_b=K V_b/k on ordinary L2(dy;C2).
The prior energy identity gives ||K||<=1 and

 ||S_b||<=.385+expm1(||b||infinity).

Thus ||S_b||<q_c=.3961 for complex ||b||infinity<=.011. This complex
extension is used only for the remainder. Real physical fields still
obey the admitted Lipschitz, mass/log, fixed-period and regular-disk
hypotheses. It does not admit an arbitrary physical L-infinity ball.

The exact inverse and its N-term approximation are

 R=k^(-1/2)(I+S_b)^-1 K k^(-1/2),
 R_N=k^(-1/2) sum_(r=0)^(N-1)(-S_b)^r K k^(-1/2),
 R-R_N=k^(-1/2)(-S_b)^N(I+S_b)^-1 K k^(-1/2).       (1)

The FIRST unsupported step in the preceding ordinary argument was taking
a coincident trace and full spectral sum of its operator-norm error.
A fixed positive q_c^N times sum d_n still diverges. Equation (1) alone
does not license that passage. Sections 2-5 address the passage in trace
norm; section 6 identifies precisely why the strong-output passage remains.

## 2. Hilbert--Schmidt control on the actual finite interval

Each scalar entry of G0 is exp(-integral k) on its directed triangle.
For either direction, integrate its square against k at the output:

 integral k(y) exp[-2 integral_y^s k] dy <=1/2.

The reversed triangle has the identical identity. The two matrices J
and O^T sigma3 O are Frobenius-orthogonal, also for complex b. Consequently

 ||S_b||HS^2 <= integral [alpha^2/k
                       +M^2 |exp(b)-1|^2/k] dy.     (2)

This is a genuine integral-kernel estimate, not a point evaluation of an
operator bound. The directed integrals run to the actual regular bolts;
the paired seam is traversed. There is no discarded return or endpoint.

The polynomial Q' has the root z=0 and degree at most five. On each half
there are at most five monotonicity intervals. On each interval theta=
atan[p/(M sqrt(Q))] varies by at most pi/2. Including both halves gives
integral |alpha|<=5pi/2. Using the prior all-mode |alpha|/k<.385 yields

 integral alpha^2/k <.385 (5pi/2).

Also M^2/k<=M. Write L for the doubled proper length. A quantitative
bound on L must use the original tiny positive temperature. In the exact
bolt factorization of current-coincident-dirac-action-assembly.md, let
t=r/B>=v=1+h. With A=3beta+u+eps,

 Q=(t-v)[2gamma v+(t-v)H1(t,v)]/t^3,
 H1=beta(t+2v)+A/(tv^2)-2eps(2t+v)/(t^2v^3)
                            +eps(3t^2+2tv+v^2)/(t^3v^4).

The enclosed seed has beta>.16, A>0, u>0, eps<0, |eps|<1e-28,
v>=1. Hence H1>=.16t-6|eps|>=.15t on the ENTIRE t>=v interval.
For s=t-v between 0 and 1 the proper integrand is bounded by
1/[v sqrt(s(2gamma+.15s))]. For t>=v+1 it is bounded by
1/[sqrt(.15)t(t-v)]. Integrating both pieces, then both halves, gives

 L <= 4/sqrt(.15) asinh(sqrt(.15/(2gamma)))+2/sqrt(.15)
   <394.854<395.                                     (3)

We used log(v+1)/v<=1, not an artificial finite outer endpoint. The
nonzero gamma is retained; the zero-temperature interval would be infinite.
The current enclosure is gamma=5.05779942482658936611...e-34.

Combining (2)-(3), M>18 and ||b||infinity<=.011 gives

 ||S_b||HS^2/M
 < .385(5pi/36)+395 expm1(.011)^2 <.217<1.             (4)

The loose final guard M is convenient for a complete spectral sum. It
does not change a metric norm, temperature or action normalization.

## 3. Trace-norm remainder and complete contact variations

For an allowed observation v let V_b'[v]=-M exp(b)v O^T sigma3 O.
The mass-source remainder is the trace of (R-R_N)V_b'[v]. In (1), the
outer factors satisfy

 ||k^(-1/2)||<=M^(-1/2),
 ||k^(-1/2)V_b'[v]||<=M^(1/2) exp(.011)||v||infinity.

They are bounded on the PHYSICAL unweighted Hilbert space for every mode.
There is no use of an unbounded similarity to assert trace invariance.
Two S factors are Hilbert--Schmidt. The elementary trace-ideal inequality
||AB||1<=||A||HS||B||HS and the bounded-factor ideal inequality prove

 ||(R-R_N)V_b'[v]||1
 <= exp(.011) ||v||infinity M q_c^(N-2)/(1-q_c).     (5)

Those inequalities follow directly by Cauchy--Schwarz in orthonormal
sequences; their hypotheses were established in (2)-(4). See also
[Melrose's Hilbert-space notes, section 23, Propositions 3.21 and 3.23](https://math.mit.edu/~rbm/18.102Sp20/Chapter3.pdf).
They require no smooth-boundary elliptic theorem.

For fixed N the tail is integrable as an action variation. Indeed the
absolutely convergent trace series

 t_N(b)=sum_(ell=N+1)^infinity (-1)^(ell+1) Tr(S_b^ell)/ell

has derivative Tr[(-S_b)^N(I+S_b)^-1 S_b'[v]], equal to (5)'s trace by
cycling bounded factors. The action -Tr log A has the negative of t_N.
Its derivative includes V_b''[u,v], and further derivatives include all
higher exp(b) vertices. Thus the tadpole is not removed by retaining only
the bubble. N below depends on seed spectral labels, not on b; no
derivative of a truncation index or changing thermal period is introduced.

The remainder has an L1 source density: cycle the scalar multiplication
by v to the end of the trace-class product. Restricting a trace-class
functional to multiplication operators is represented by its integrable
diagonal, with L1 norm at most the trace norm. This follows by a nuclear
sum of rank-one kernels and L2 Cauchy--Schwarz. It does not assert that an
arbitrary operator-norm remainder has a diagonal. The density in dz is
the proper-density pushforward; its L1 norm is coordinate invariant.

## 4. An index-dependent order, with no omitted modes

Choose, for EVERY n>=1 and EVERY signed thermal label,

 N(M,p)=128+ceil{[4 log(M/18)+2 log(1+|p|)]/log(1/q_c)}. (6)

This is an infinite-spectrum rule, not a finite selected-mode computation.
No angular or thermal level is dropped. For an implementation using
enclosed seed constants choose an integer satisfying the following
inequality for the whole enclosure; extra terms are harmless:

 q_c^(N-2)<=q_c^126 (18/M)^4 (1+|p|)^(-2).          (7)

The exact angular counting function is at most 4C M^2 and M_n>18, so
Stieltjes integration gives sum d_n/M_n^3<2C/3. For the entire positive
half-integer thermal lattice, monotonicity gives

 gamma sum_(j>=0) [1+gamma(j+1/2)]^-2 <=1+gamma.

Multiply by TWO for the negative frequencies. Equations (5)-(7) now
give an absolutely trace-norm summable remainder, including the full
degeneracies. No continuum-frequency replacement is made.

The exact identical paired physical action is half the doubled action.
Its original per-half normalization has magnitude pi B^2/(G gamma).
For a conservative bound we use the full doubled trace (factor two larger
than the half) times G gamma/(pi B^2). With G=1 in the existing units,
the normalized SOURCE tail is bounded by

 E1=exp(.011)/(1-q_c) q_c^126 18^4
                         *4C(1+gamma)/(3pi B^2)
   <9.204e-79.                                      (8)

Here C=49119682270.77026011517..., B^2=8.38309244070745772632...e42.
The temperature factor cancels against the exact thermal count; this
does not renormalize an admission target. The unnormalized trace error
is also finite. No spectral sum of retained low-word terms is yet evaluated.

On the real admitted subfamily ||b||infinity<=.001, use a complex mass
polydisk of total radius .01. Cauchy in each normalized direction gives
the following bounds for the same WHOLE remainder:

 ||source tail||L1                         <= E1 <9.204e-79,
 ||Hessian tail[u]||L1                    <=100 E1 ||u||infinity
                                                    <9.204e-77 ||u||infinity,
 ||third action tail[u,w]||L1             <=40000 E1 ||u||infinity||w||infinity
                                                    <3.682e-74 ||u||infinity||w||infinity. (9)

These include the differentiated exponential contact and the varied
paired/finite-domain propagators. Actual phase/input supremum norms must
be inserted. They are not assigned one; a coefficient allowance is not
an actual coefficient. The LLL and separate Weyl/local contacts are kept
exactly in the retained part, not bounded or altered by (9).

## 5. Common regulator and its metric derivatives

At positive common cutoff epsilon the exact mass-force formula contains
E_epsilon(b)=exp[-epsilon A_b^*A_b] multiplying A_b^-1 delta A_b
(with real part for the adjoint-square definition). Insert (1) into THIS
formula. Keep the same E_epsilon, actual domain, local heat subtractions,
seam/boundary completion and finite restorations on the retained part.
Do not assign a heat kernel to R_N or use a different subtraction there.
The geometric force also retains the derivative of the EXACT conversion
J_e(b)=sum_n d_n integral_0^1 Tr_bar[b exp(-e P_r(b))]dr,
P_r=A_r^*A_r, A_r=exp(-rb/2)A_b exp(-rb/2). This conversion, its boundary
completion and all finite restorations stay in the retained center.
The present tail split is applied only to the mass-action part; it does
not equate geometric and mass heat traces at positive cutoff.

For the trace-class tail E_epsilon converges strongly to I and is a
contraction. Equation (5), then the complete majorant (7), justify both
the regulator limit and both spectral sums. In particular no local
divergent counterterm belongs to this tail. All original subtractions
and finite contacts must remain in the retained center. This is an exact
error split of the original regulated force, not a new physical action.

Metric derivatives commute with this limit on the stated mass tube.
Here is the domain argument needed for that assertion. Let A_0=D+M on
the original domain. D is the already established skew realization of
the paired mass representation. Thus ||M psi||<=||A_0 psi||. On the
same closed form domain consider the holomorphic form associated with
A_(bar b)^* A_b. Its difference from the b=0 form is bounded by

 (2e+e^2)||A_0 psi||||A_0 phi||,
 e=expm1(.011), 2e+e^2=expm1(.022)<.023.

This estimate holds on the full direct sum as well. The form is coercive
and sectorial, with common domain and uniform angle. Its resolvent can
be obtained by a Neumann series between the unperturbed form isometries;
the resolvent is holomorphic and maximal accretive. The implicit-Euler
resolvent products give the contractive heat semigroup, holomorphic in
the bounded mass parameters for positive time and strongly tending to I.
For real b this is precisely the original adjoint square, not a square
with its operator domain frozen. This uses a fixed FORM domain only.

Apply Cauchy's formula to the trace of this heat operator times the
trace-class remainder on the parameter polydisks. The trace is dominated
by (5)-(7), uniformly in epsilon and on their boundaries. Dominated
convergence passes the regulator limit through those integrals. Terms
arising from derivatives of E_epsilon therefore have their proper limit;
they were not silently dropped. The result is the derivative of the
convergent unregulated tail action in section 3, with (9)'s bounds.

The field theory's remaining Weyl allocation and common restored local
action are unchanged. The prior qualitative admitted-family action proof
is still needed for that retained part. The present trace-ideal argument
does not provide its strong output jets or extend its physical domain.

## 6. Attempt the original coupled test; precise failure of a conversion

The new remainder is controlled in L-infinity -> L1, including every
mass/quartic/correlated-slow/fast/interior input that has its ACTUAL
supremum norm supplied. It is not a bound in the original graph and does
not evaluate any of the retained signed entries. In particular there
is NO valid generic constant converting this L1 error to an H3 forcing
norm. Compact oscillations supported strictly inside the original collar
have uniformly bounded L1 norm and unbounded differentiated H3 norm.
They satisfy the actual boundary/bolt conditions because they vanish
near both ends. Taking an increasingly large constant cannot repair that
space mismatch. This observation does not assert that the actual tail
has those oscillations or that its strong norm is large.

The original localized coupled rows are, without dropping either equation,

 V=chi[-F-A(m Phi_m)-R w]+[A0,chi]w=A0(chi w),
 G=(D-4)V_b-(D-6)V_nu/3,   w=u-m Phi_m.              (10)

An L1 estimate supplies neither the derivatives in (10) nor the phase
moments, the independent interior source norm, or the bolt return. The
mass determinant has canonical nu row zero; its shape row still couples
through the ORIGINAL two-row inverse/York return. No scalar equation or
scalar quartic is deleted. A sufficient local remainder norm for (10)
would control the relevant weighted W^(4,2) source jets, as well as the
separately specified phase moments and interior/bolt maps. The precise
bound on these ACTUAL remainder maps, not a nonexistent L1 embedding,
is the next missing argument. Boundary cutoffs and their differentiated
contacts must stay together in obtaining it.
For an explicit comparison take V_nu=0, V_b=psi(x)sin(Lambda x),
x=log z, with nonzero smooth psi supported strictly inside the collar.
Then ||V_b||L1 is bounded independently of Lambda, whereas the leading
term in partial_x^3[(partial_x-4)V_b] has size Lambda^4 psi sin(Lambda x)
up to phase and O(Lambda^3). Its L2 norm divided by Lambda^4 tends to
||psi||2/sqrt(2). Thus the particular L1-to-G-H3 comparison fails on
actual-domain compact sources, not just on a smooth surrogate interval.

No threshold on a strong multiplier is inferred from (9), because the
generic multiplier is unbounded. For an already justified weak observable
with bounded mass test v, however, the executed error is quantitatively
E1 ||v||infinity, and the Hessian pairing error is
<9.204e-77 ||u||infinity||v||infinity. This is its explicit contribution
to a future complete coupled certificate. It replaces the formerly
divergent raw trace error; it does not replace a missing strong estimator.

The retained low-word center in (6) is itself an infinite angular/thermal
sum with the ORIGINAL common subtraction and restoration. It is not
computed here. R_N is not a renormalized determinant on its own, and a
tiny tail says nothing about that signed center's size. Increasing N
alone will not close any remaining strong-output or center requirement.

## 7. Unchanged feasibility ledger and decision

The exact reference mass lift remains
J_m=int_0^z t^3/Q^2, n_m=-2beta QJ_m/z,
b_m=beta z^3/Q+beta Q'J_m. The uncut areal f=2xz, alpha=0 cancels BOTH
local reference rows. Its actual chi_m transition (.1 to .2) retains
f_cut=2xz chi_m-2F n_m Dchi_m and alpha_cut=n_m Dchi_m. The original
boundary conditions and regular-bolt matching stay; the reference cut
lift has zero bolt shift and net mass -xB per coefficient. It is not a
global fixed-temperature zero mode. The scalar quartic and the correlated
B0(z^4-z^(4+e_*))+Delta z^(4+e_*) remain exact and coupled.

The original graph conversion 2133<||d_chi||H3<2174 and established
LINEAR |m|<1.627532e-15 still debit <=3.538254568e-12. With the distant
shape bound <1.24e-20, near SHAPE has >6.4617454e-12 available. Scalar
graph and phase conditions have their own allowances.

kappa_m<=5e-5 remains UNKNOWN. If it holds, its conditional debits are
2.14765e-12 affine and .0024731 feedback. All OTHER entries must still obey

 Y_other<=5.34234999994958e-12,
 eta_other+(.0195239+4.583e-13)gamma_E<.0012069.

Even eta_other=0 requires gamma_E<.061816543. None of these unknown
coupled coefficients was computed or assigned its allowance this turn.

The proved two-fast admitted-family screen remains
.028913+1e-11(1.0706e10)=.135973<.1599164. Its margin is .0239434.
It supplies neither independent boundary-phase columns nor the missing
original physical-X comparisons. Nonlinear admission/contraction still
requires, in the original units and with ACTUAL constants,

 E_center+1e-8 alpha_embed<=1e-11,
 R_X alpha_embed*.028913+Z0_other<=.25,
 R_X alpha_embed^2*1.0706e10+Z2_other<=5e7,
 Y0+1e-8 Z0+5e-17 Z2<=1e-8,
 Z0+1e-8 Z2<1.

RESEARCH DECISION: retain the construction and sharpen the specified
Dirac method to the complete-spectrum adaptive trace split (6). Its
common-regulator weak tail is now quantitatively controlled. Stop
optimizing its small number. Quantitative attainability of the complete
coupled graph is still undecided; no controlled solution is established.
The generic weak-to-strong shortcut fails, not the construction.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): compute the common-regulated
signed bordered return block for the affine source and all mass, quartic,
correlated-slow, fast and interior columns, using (6) and certifying its
ACTUAL differentiated source remainder and original phase/bolt outputs.
The calculation must retain both rows and their full local restoration;
compare the resulting complete errors/centers with the above residual
budgets. No isolated contact or finite-mode proxy supplies that block.

EXECUTED: ordinary actual-domain Hilbert--Schmidt/trace-ideal proof,
complete spectral summation of the adaptive remainder, common-regulator
and mass-variation argument, outward bounds (3)-(9), and rejection of a
generic L1-to-H3 conversion. NOT EXECUTED: the retained signed center,
strong differentiated tail, full bordered solve or kappa_m. Supporting
checks verify arithmetic and identities; they are not the domain proof.

Remaining to-do list: signed complete coupled strong return and actual
phase/interior columns; original scalar graph and physical-X/neighborhood
comparisons; nonlinear admission/contraction; exact compatibility;
physical matching; semiclassical validity/noise; all-sector stability;
charge/return, opening, transport, preparation, separate creation and the
full stable semiclassical wormhole objective.
