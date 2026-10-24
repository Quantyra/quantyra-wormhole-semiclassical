# Geometric conversion before the quantum projection is summed

2026-10-23, E003/S169. Ordinary argument before frozen execution.
Previous turn: progress, certified magnetic and complete cold components.
The full stable-wormhole objective stays active and incomplete. The audit
stays closed. This addresses the representation/domain prerequisite of
the remaining quantum projection, not another cold-source optimization.

## 1. Identify the missing finite-mode functional on the actual domain

Use the same physical-half rank-two orbit Dirac realization Dsharp,
paired complementary domains, fixed period, fixed scale and positive
mass M as current-canonical-heat-restoration-assembly.md. On its fixed
spinor Hilbert space Dsharp is anti-self-adjoint. Identical scalar
multiplications preserve its first-order paired domain. The original
unpaired Euclidean Dirac is NOT being declared self-adjoint.

The existing ODE evaluator varies the mass of Dsharp+M exp(b) at fixed
orbit metric gamma0. The geometric operator is instead

    A_geo(b)=exp(-b/2) Dsharp exp(-b/2)+M.

For the determinant modulus, define a path between them:

    A_t=exp(-t b/2)[Dsharp+M exp(b)]exp(-t b/2),
    gamma_t=exp(2t b)gamma0, m_t=M exp((1-t)b),
    partial_t A_t=-{b,A_t}/2, 0<=t<=1.              (1)

Work first at a smooth finite outer cut with the ACTUAL paired projector,
and include its induced square boundary condition. For bounded real b,
Re<A_t psi,psi>>=min(m_t)||psi||^2 in the corresponding geometric
Hilbert space, and the adjoint has the same positive real part. Closed
range and the zero adjoint kernel give invertibility. There is no omitted
massive zero-mode normalization or phase crossing on this path. The
massless LLL and its declared kernel-free allocation remain separate.

Let P_t=A_t^dagger A_t and Pbar_t=A_t A_t^dagger. Trace cyclicity first
in the trace-class half-plane gives

    partial_t zeta_P(s)=s[Tr(b P_t^-s)+Tr(b Pbar_t^-s)]. (2)

For Gamma=+zeta_P'(0)/2+log(mu_bar)zeta_P(0), (2) yields the average
of the two smeared a2 coefficients. Here mu_bar=mu0 B/sqrt(C) stays
fixed. The scale term has zero t derivative along this homogeneous
path. This derivation retains both squares, not a guessed determinant
multiplicativity rule or a finite-dimensional determinant analogy.

The rank-two Lichnerowicz potential has trace
tr E=-R_t/2-2m_t^2; the extra Clifford derivative of m_t is traceless
and changes sign between the partners. Consequently both bulk
coefficients, and hence their average, are

    a2_bulk(b)=-integral dvol_t b[R_t/(24pi)+m_t^2/(2pi)]. (3)

At the cut, the balanced projectors have equal rank. The normal derivative
of the smearing has zero traced coefficient. With S=-K_t Pi_+/2 the
remaining boundary coefficient is -integral ds_t b K_t/(12pi).
No mass Robin condition is added: scalar m_t preserves the first-order
domain, so composition induces the same boundary subspace. The partner
exchanges projectors and gives the same trace. This is the local a2
calculation on a smooth physical side, not a heat theorem across a
smooth reflected seam.

Since dvol_t R_t=dvol_0(R-2t Delta b),
dvol_t m_t^2=dvol_0 M^2 exp(2b), and
ds_t K_t=ds_0(K+t partial_n b), integrating (3) and its boundary term
gives the exact finite-M conversion

    Gamma_geo,M[b]=Gamma_mass,M[b]+C_M[b],
    C_M[b]=-(1/(24pi)){integral [b R+|grad b|^2]
                              +2 integral_boundary b K}
                       -(M^2/(2pi)) integral b exp(2b). (4)

Both actions use the same modulus and auxiliary two-dimensional zeta
scale, and agree at b=0. The b partial_n b terms cancel between the
integrated bulk and boundary coefficients. Formula (4) holds for every
positive M; no small-mass expansion is invoked.

For the retained outer phase b=c3 z^3+O(z^4 log^p z), b and its first
normal derivative vanish at the limiting seam. K=O(z), so bK=O(z^4)
up to logarithms. The finite-cut boundary term tends to zero; its
Hessian is already zero in this fixed affine chart. The actual bolt is
a smooth disk point, not a boundary. Take the SAME finite-cut limit on
the two sides of (4): the explicitly computed difference and its mixed
derivative have this limit. This does not assert a smooth C5 double or
commute an infinite angular limit with metric differentiation. A separate
uncut norm theorem is still required for that infinite sum.

The ordinary source tools are the paired-square/domain proofs in
stationary-paired-dirac-phase.md, stationary-exterior-conformal-domain.md
and stationary-exterior-boundary-weyl.md. The heat conventions and
determinant-variation rule agree with [Vassilevich, sections 5.3 and 7](https://arxiv.org/html/hep-th/0306138v3).
[Namuduri--Bilal, section 5.1](https://arxiv.org/html/2308.05802v2#S5.SS1)
provides a separate smooth-surface check on the curvature and massive
contact structure (their Majorana action is half our Dirac action).
Their small-mass series supplies no bound for M1>18 or for our paired seam.

## 2. The required contacts in the unchanged physical normalization

Divide by W=-N pi B^2/(G gamma), so an orbit density is multiplied by
-2G/B^2. Let cL=G/(12pi B^2), cM=G M^2/(pi B^2), and use
Delta=(Q partial_z)' with R=-Q''. The conversion in I/W units is

    C_M/W=cL integral [R b+Q b'^2]+cM integral b exp(2b),
    Gconv[v]=integral (cL R+cM)v,
    Hconv[h,v]=2cL integral Q h'v'+4cM integral h v.  (5)

The boundary source term from (4) is retained until its physical limit;
it is not an extra Hessian in (5). The density and exponential tadpole
are essential. At the exact product R_*=-2, the product-gradient
subtraction of the preceding milestone changes this conversion to

    Hconv[h,v]-2gconv,* integral h v
       =2cL integral [Q h'v'+2h v]+2cM integral h v. (6)

The last term is NOT removed by the pure-area a0 identity. The functional
b exp(2b) is not a constant multiple of exp(2b). It must cancel against
the correlated part of the mass determinant before separate bounds are
taken. Formula (6) is a representation term, not a new finite coupling
and not an additional physical source on top of Gamma_geo.

There is an exact continuum flat-density test. At constant b, the
rank-two mass determinant per area is

    V_mass(b)=M^2 exp(2b)/(4pi)
                    [log(M^2/mu_bar^2)+2b-1].        (7)

Adding the mass term in (4) gives
V_geo(b)=exp(2b) M^2[log(M^2/mu_bar^2)-1]/(4pi).
Thus V_geo''(0)-2V_geo'(0)=0 exactly, whereas
V_mass''(0)-2V_mass'(0)=M^2/pi. The conversion contributes -M^2/pi.
This test checks the sign, factor and product contact in continuum
normalization; a finite matrix determinant would not reproduce the anomaly.

## 3. Quantify the loss if the conversion is bounded separately

For a finite angular cutoff N, retain the ACTUAL masses and multiplicities
M_n^2=n(n+q)/C, d_n=2(q+2n). Exact finite sums give

    sum_1^N d_n=2N(N+q+1),
    sum_1^N d_n n(n+q)=N(N+1)(N+q)(N+q+1).          (8)

The normalized absolute coefficient of the flat product-subtracted
conversion contact is therefore

    K_N=2G N(N+1)(N+q)(N+q+1)/(pi C B^2).           (9)

Its equal and opposite mass-determinant term cancels algebraically for
every N. Summing unsigned pieces first would diverge as N^4.

Use the actual first correction and the global shape test v=-z^4.
For the same bulk cut z_c=1-2^-49, write
h=b_h=a/2-z w-z^2 w'/2. Integration by parts gives

    integral_0^zc h v
       =integral_0^zc [a v/2+(z^2 v'/2)w]
                                  -[z^2 v w/2]_0^zc. (10)

Use the previously certified A0=3.768e-12, a_loc=4e-28,
f0=4.34e-11, s0=3.33e-10 h_bolt and
L=-log(1-z)-z:

    |a|<=A0 L+a_loc,
    |w|<=s0+2(1-z)(A0 L+f0+a_loc)/z.

The positive envelopes may be extended from z_c to 1 only inside the
integral upper bound. Exact moments are
integral z^4 L=29/100,
integral z^4(1-z)L=103/4200, integral z^4(1-z)=1/30.
Retaining the finite-cut endpoint, this gives

    |integral_0^zc h v| <= U_pair,
    U_pair=A0(29/200+4*103/4200)+4f0/30
            +a_loc(1/10+4/30)+s0/3+zc^6 w_bound(zc)/2. (11)

Freeze U_pair<6.704e-12. Substitute the actual interval B,C,h and
F_lambda into K_N U_pair/F_lambda for N=1,2,16,128. Compare these
KNOWN upper allowances with the remaining 1.287e-40 working target.
They are neither measured source values nor lower bounds; their purpose
is to test the invalid strategy of treating the two correlated contacts
as independent errors. The exact cancellation (7) costs zero at every N.

The same cancellation can be executed directly in the continuous thermal
integral, independently of the zeta density (7). At flat proper thermal
momentum p, the mass gradient is -M^2/sqrt(M^2+p^2). Differentiating
the exponential mass including its tadpole and subtracting twice the
gradient leaves M^4/(M^2+p^2)^(3/2). Both frequency signs give

    integral_0^Lambda dp M^4/[pi(M^2+p^2)^(3/2)]
                         =M^2 Lambda/[pi sqrt(Lambda^2+M^2)]. (11a)

The integral to infinity is M^2/pi, exactly the contact in (7).
Consequently the conversion cancels it before numerical summation.
If this one known contact were instead truncated at Lambda=kappa M,
its fractional omitted tail would be 1-kappa/sqrt(1+kappa^2), bounded
by 1/(2kappa^2). Evaluate the sufficient kappa required to bring the
separate allowance (9),(11) below 2.5e-43. Integrating (11a) analytically
avoids this otherwise unnecessary resolution cost. This is a local
continuum contact test, not a zero-temperature replacement for the
actual finite-temperature variable-geometry remainder.

## 4. Preserve the common four-dimensional prescription

Do NOT assign independent finite values to sum d_n or sum d_n M_n^2
in (5). Do NOT separately zeta-sum the conversion and mass determinant.
Use (4) at finite N, combine each mode's terms, then subtract its
GEOMETRIC heat polynomial before taking the common sum/finite part.

For one mass the geometric local allocation is

    P_M[b]=A0[b] M^2[log(M^2/mu_bar^2)-1]/(8pi)
               -A1[b] log(M^2/mu_bar^2)/(8pi)
               +A2[b]/(8pi M^2),
    R_M[b]=Gamma_mass,M[b]+C_M[b]-P_M[b].           (12)

Here A0=2 integral exp(2b), A1=-(1/6)integral(R-2Delta b), and
A2=-(1/120)integral exp(-2b)(R-2Delta b)^2, including the paired
divergence completion. A1 has zero variations on the physical phase.
Restore exactly the previously computed common Gamma_a0,Gamma_a2
and fixed scale ONCE. Their complete geometric heat remainder is
sum d_n R_M, provided the differentiated sum is justified. Equation
(12) introduces no new scheme or independent mode finite coefficient.

The bare transfer-ODE determinant, even at fixed n, has an uncut thermal
divergence. Its quadrature must be matched to the auxiliary mass modulus
in (4) with an explicit thermal subtraction before being combined with
(12). The old compact finite-frequency matrices do not perform this
matching. Neither (4) nor the flat test evaluates that finite remainder.

For a correctly admitted full-domain pair or a smooth observation
partition with all contacts retained, let the combined product-subtracted
reduced heat remainder obey |D_hD_v r(t)-product_contact(t)|
<=A_hv t^(1+alpha), alpha>0, for 0<t<=t0. The full counting bound
Theta(t)<=4C/t then gives a normalized omitted cap

    E_cap <=4GC A_hv t0^alpha/(alpha B^2 F_lambda). (13)

Allocate E_cap<=2.5e-43, one quarter of the selected 1e-42 total error.
Its required threshold is
A_hv t0^alpha/alpha<=2.5e-43 B^2 F_lambda/(4GC).
Evaluate this coefficient at the actual geometry, including the example
alpha=1,t0=.001. A_hv and alpha are UNKNOWN for this differentiated
paired remainder, not supplied by a checker or an undifferentiated
heat bound. The common regulator and angular limit must be removed
uniformly in the two directions. A sharp characteristic truncation of
v is not an admitted metric variation; the preceding observation/domain
warning remains in force.

For a separate full-time angular-tail bound on the time-reduced,
product-subtracted R_M pairing, suppose |R_M,hv|<=B_hv/M^4.
The EXACT counting derivative is
d_n/M_n^4=-2C^2 partial_n[1/(n(n+q))]. Its positive decreasing
continuation gives sum_(n>N) d_n/M_n^4<=2C^2/[N(N+q)]. Therefore

    E_angular<=4GC^2 B_hv/[B^2 F_lambda N(N+q)].     (14)

Assign it another 2.5e-43 and evaluate its B_hv threshold at the same
four cutoffs. B_hv is UNKNOWN. A short-time A_hv alone is not a bound
on this full-time coefficient. Reserve the remaining 5e-43 for thermal
quadrature/tail and spatial/rounding errors; neither is currently proved.
These are explicit sufficient allocations, not new measured errors.

## 5. Research decision

Use the geometrically converted, heat-subtracted mass representation
(4),(12) to evaluate the selected quantum/contact projection. This closes
the finite-mode representation/contact identification and identifies an
exact cancellation essential to its numerical feasibility. It does not
close the infinite thermal/angular remainder or certify a smaller full
source. Do not replace the selected projection by more compact-probe
matrix checks or assign zeros to the missing local moments.

The SINGLE next calculation remains the complete signed quantum/contact
projection with |P_hat_rest|+E_rest<1.287e-40, now requiring the actual
thermal subtraction and a differentiated common-tail certificate for
(12). An error <=1e-42 with |P_hat_rest|<1.277e-40 is a sufficient
allocation. The current execution tests the conversion/contact strategy;
the full projection and actual controlled solution remain UNEXECUTED.
The finite-cut boundary identity precedes its interpretation, and the
remaining infinite-limit regularity requirement precedes more source
optimization. The full nonlinear, physical and stability goal is retained.

Remaining to-do list: complete the geometrically converted thermal/angular
projection and differentiated common-tail bound; global Schur forward/adjoint
solution, physical bolt/exterior return and admission; nonlinear closure and
exact compatibility; physical matching, causal realization, validity/noise,
all-sector stability, charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
