# Explicit global Green input with its original equation residuals

2026-10-28 HST, E003/S169. The preceding filtered-output comparison was
progress. The full signed quantum certificate remains UNRESOLVED. This
calculation constructs and encloses the actual auxiliary Green INPUT,
including the global error measured by the original operators. It does
not evaluate T_t, F_t or Z_r. Preserve P32, the closed audit and all
accepted results and physical norms.

## 1. Certificate and the errors this calculation must supply

Use H_mu=L2(z^-4 dz), the actual normalizable/regular-bolt second-order
realizations, and the paired mass/log quantum form E_t. All clocks,
periods, spin, subtraction scale, local contacts, full angular/thermal
sums and original cutoffs remain unchanged. With s=L-z, the exact rows are

    L_D=1-d[z^2 Q partial_s^2+(z^2 Q_s+2zQ)partial_s],
    Pbar=-z^2(Q partial_s)' + v(z),
    v(z)=2z^2/C+12lambda0 A_radius^2.                         (1)

The source normalization is c=beta alpha/2, alpha=(1+h)/h. Let m_b be
evaluation at L in the H_mu dual, h_b=c L_D^-1 m_b and t_0=Pbar^-1 h_b.
The exact point singularity is defined by the original resolvent, not by
a new puncture boundary. The actual averaged source and point replacement
retain the errors already paid in current-filtered-mixed-green-residual.md.

For explicit h_hat,t_hat the remaining input errors there are

    e_D=||h_hat-h_b||H_mu,
    e_P=||Pbar t_hat-h_hat||H_mu.                              (2)

Their contribution to the COMPLETE three-filter bolt certificate is

    E_input <= E_green,0+K_green(1+r_mix)(e_D+e_P),
    E_green,0<7.415e-12, K_green<2.162e-42,
    r_mix=1.556e-23.                                        (3)

That comparison retains both original rows, their action normalization,
and the full quantum inverses in the filters. Its proposed one-percent
allocation is not an available margin. The actual full certificate is
still pref_t(c_a+c_tau+Zr_hat+e_Zr)+E32+Ereg+e_other<=1.667e-9.

## 2. One explicit trial, with exact singular coefficients

Set a=Q''(L)/2, Q2(s)=2gamma s+a s^2 and v0=v(L). Define

    Delta0=L^2(Q2 partial_s)',
    L_D0=1-d Delta0, P0=-Delta0+v0,
    ell(ell+1)=v0/(L^2 a), nu(nu+1)=1/(d L^2 a),
    w=1+a s/gamma,
    H(s)=A_h Q_nu(w), T(s)=A_t[Q_ell(w)-Q_nu(w)],
    A_h=c L^2/(d gamma), A_t=c L^2/[gamma(1-d v0)].           (4)

Q_k here is the standard real Legendre function on w>1. These degrees
belong to auxiliary geometric equations, not selected matter modes.
Use chi_S=1 for s<=S, zero for s>=2S, and, on S<s<2S,
chi_S=1-10x^3+15x^4-6x^5 with x=(s-S)/S. Then

    h_hat=chi_SD H, t_hat=chi_SP T,
    SD=1e-36, SP=1e-16.                                    (5)

These cutoffs belong to the computational trial. They impose no physical
boundary: the errors below use the ORIGINAL global inverses and include
every derivative on both cutoff transitions. The original mass/graph
cutoffs are untouched. Zero outer data for a trial do not set the actual
Green field's generated outer coefficients or phase returns to zero.

The first continuum check is the source and operator domain, before any
numerical interpretation. The Legendre equation implies, away from s=0,

    L_D0 H=0, P0 T=H.                                     (6)

The second identity uses BOTH terms of T before taking norms:
A_t(1/d-v0)=A_h. Their identical logarithmic singularities cancel.
The standard normalization gives Q_k(w)=-log(w-1)/2+
log(2)/2-EulerGamma-digamma(k+1)+O((w-1)log(w-1)). See
[DLMF 14.8.9](https://dlmf.nist.gov/14.8.E9) together with
[14.3.10](https://dlmf.nist.gov/14.3.E10) for the conversion from Olver's Q.
Thus T(0)=A_t[digamma(nu+1)-digamma(ell+1)] is finite and its
Q T_s flux vanishes at the bolt. No delta remains in Pbar t_hat.

For h_hat the exact original flux is
-d lim_(s->0) Q h_hat_s/z^2=d gamma A_h/L^2=c. Hence
L_D h_hat=c m_b+r_D distributionally. The displayed estimates show
r_D in H_mu. Cancellation of this exact source flux makes h_hat-h_b
the original L_D inverse of r_D; there is no arbitrary logarithmic or
homogeneous boundary coefficient. Therefore e_D<=||r_D||H_mu, since
||L_D^-1||H_mu<=1. The C2 transitions add no delta terms. Both trial
fields vanish near the paired outer seam and retain the original weak
bolt realization, with the specified source singularity.

The filtered map established previously acts on these H_mu arguments
through the bounded B source. It does not apply a C5 heat theorem to
the point field or promote it to a nonlinear admitted metric.

## 3. Whole-interval bounds and checked special-function hypotheses

For real k>=1, w=cosh(eta)>1, the positive integral is

    Q_k(w)=int_0^infinity [w+sqrt(w^2-1)cosh(u)]^(-k-1) du. (7)

This is [DLMF 14.12.6](https://dlmf.nist.gov/14.12.E6) at order zero
with the same standard-Q conversion. The hypothesis k>-1 holds; all
branches in (7) are positive real. The execution verifies ell>=1, nu>1.
The cutoff, domain and finite-geometry conclusions do not come from DLMF.

From (7), Q_k<=Q_1 for k>=1. Splitting x=a s/gamma at one, use
Q_1(1+x)<=log(3/x)/2 for x<=1 and <=1/(3x^2) for x>=1.
Integrating x^2 times these squares gives
[(log 3)^2/3+2log(3)/9+2/27]/4+1/9<1. Consequently

    ||s Q_k(1+a s/gamma)||L2(0,infinity)<(gamma/a)^(3/2).

The positive Legendre equation, its source flux and decay also give
0<=-Q2 partial_s Q_k<=gamma. For cutoff tails put

    U_k(S)=pi exp[-(k+1)eta]/[2 sqrt(1-exp(-2eta))],
    V_k(S)=pi gamma(k+1)exp(-k eta)/4,
    eta=acosh(1+a S/gamma).                                (8)

These bound Q_k and -Q2 partial_s Q_k on s>=S. To derive U, use
cosh(u)>=1+u^2/2 in (7), then bound the remaining power by its k=1
integral. Integrating this bound against sinh(eta) and using the
Legendre equation gives V. No sampled radial values replace these
whole-interval estimates.

Keep Q2 and the actual Q distinct. The exact polynomial gives
|Q'''|<=M=1000 on the whole interval, hence
R=Q-Q2 obeys |R|<=M s^3/6, |R_s|<=M s^2/2. On the complete support
z>.99, Q/Q2<2 and the H_mu/L2 norm factor is <2. Define

    A=2+M/(6a), B=5M/(6a), g=gamma/a.

Eliminating f_ss by the model equation for EACH Legendre term gives
|z^2Q/Q2-L^2|<=A s and
|z^2(R_s-R Q2_s/Q2)|/Q2<=B. The original L_D drift adds at most four
to B. Thus the entire interior residuals have upper bounds

    D_geom=2 A_h[A g^(3/2)/L^2
                        +d(B+4)gamma sqrt(2SD)],
    P_geom=2 A_t[(A(v0+1/d)/L^2+8/C)g^(3/2)
                        +2B gamma sqrt(2SP)].               (9)

The factor 2 is a proved measure bound, not a changed norm. The large
source normalization remains in A_h and A_t. Terms proportional to
Q''' and 1/d are retained, even when generous.

For a cutoff radius S set

    Qmax=4gamma S+4a S^2+(4M/3)S^3,
    Qsmax=2gamma+4aS+2MS^2.

Since |chi'|<=2/S, |chi''|<=15/S^2, the FULL cutoff bounds are

    D_cut=2 sqrt(SD) d A_h {
      [15Qmax/SD^2+2(Qsmax+2Qmax)/SD] U_nu(SD)
                         +8 V_nu(SD)/SD},
    P_cut=2 sqrt(SP) A_t {
      [15Qmax/SP^2+2Qsmax/SP][U_ell(SP)+U_nu(SP)]
                         +8[V_ell(SP)+V_nu(SP)]/SP}.

Use Qmax,Qsmax at the corresponding radius in each expression. The
remaining mismatch between the two different cutoffs is bounded by

    P_tail=2 sqrt(2SP) A_h U_nu(SD).

Equations (1),(6) now prove globally, including all transition sources
and the omitted finite-domain return,

    e_D<=D_geom+D_cut,
    e_P<=P_geom+P_cut+P_tail.                              (10)

The actual inverses, not a reference-domain inverse, convert these
residuals into (2)-(3). Neither the exterior nor the bolt boundary was
replaced by the quadratic reference.

## 4. Correlated parameters and executed interval evaluation

Do not index the sparse seed coefficient list as a dense polynomial:
its powers are 0,2,3,4,5,6. Use the exact factorization
Q=(1-z)^2 Pseed(z)-tau z^3, Pseed(1)=1. Write delta=h/(1+h),
L=1-delta, epsilon_seed<0, u=1/C, beta=(1-u-epsilon_seed)/6,
p1=8beta+2u+4epsilon_seed, p2=6beta+2u+12epsilon_seed, and

    PL=1-p1 delta+p2 delta^2/2-4epsilon_seed delta^3
                                      +epsilon_seed delta^4,
    tau=delta^2 PL/L^3,
    a=1+a_c,
    a_c=-3p1 delta+3p2 delta^2-40epsilon_seed delta^3
                  +15epsilon_seed delta^4-3delta^2 PL/L^2.

Keep lambda0=lambda+11lambda^2/(360pi) exactly. To certify ell>=1,
compute v0-2L^2a as
2delta(2-delta)(1-u)-2L^2 a_c-2epsilon_seed+
12beta*11lambda/(360pi), without subtracting independent rounded values
near two. Recover ell-1 by rationalizing its quadratic equation. This
preserves the accepted parameter correlations; it improves no physical
normalization or admission target.

The frozen 512-bit test evaluates (8)-(10) over the complete retained
parameter enclosures and tests e_D<1e27, e_P<1e29, e_D+e_P<4e30 and the
original full input debit <8e-12. It also encloses values of the explicit
trial via (7), including the finite bolt limit of T. A direct library
Legendre evaluation returned nonfinite intervals at the large degree;
it is not used as evidence. The implemented positive integral uses an
analytic logarithm check and adds its entire u>128 tail, bounded by
(2/sqrt(w^2-1))^(k+1) exp[-128(k+1)]/(k+1).

The analytic functions (4)-(5), with correlated exact seed parameters,
are the trial; no interpolation or midpoint selection defines a different
field. The all-domain error certificate comes from (10), not from the
finite table of values. Later integrations and quantum inverse applications
still require their own numerical error accounting.

## 5. What this changes and what remains unexecuted

This supplies explicit global computational input functions and proved
original-equation errors, where e_D and e_P were previously unknown for
the proposed approximation. It closes this input accuracy contract for
the analytic trial only. It does not establish an available remaining
bolt margin or evaluate the signed three-residual center.

The first remaining unsupported numerical implication is the complete
paired quantum solve needed for T_t chi_hat and the successive F_t
applications, including the energy-dual norm of their residuals. The
existing selected-mode evaluator and quartic column do not implement
that map on this input. Full spectra, ultraviolet logarithm, common
contacts and global returns must stay in those solves. Their actual
application errors remain UNKNOWN.

Retain the spatial signed method and this explicit input. No further
Green, degree or tail optimization is selected: its contribution is
already below the stated input screen. The one next calculation is the
complete signed transfer with certified quantum applications. That is
PROPOSED, not executed by this input calculation. No original-source,
graph, other physical output, mass-column or independent boundary-phase
requirement inherits this bolt comparison.

Remaining to-do list: signed quantum certificate and application errors;
every coupled budget; nonlinear admission/contraction; exact compatibility;
physical matching/opening; semiclassical validity/noise; all-sector
stability; the full stable semiclassical wormhole objective.
