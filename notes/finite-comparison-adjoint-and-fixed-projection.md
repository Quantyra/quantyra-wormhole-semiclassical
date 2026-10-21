# Finite adjoint phase and the fixed projection's exact overlap test

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary argument
before frozen verification. **The actual finite adjoint is constructed
with a proper-radius phase enclosure. Fixed-projection transversality
is reduced to a specified signed overlap; its value remains open.**

The finite TT/range theorem permits the real adjoint construction
previously available only at the seed. At every original comparison,
the actual unit-Liouville-amplitude adjoint has scalar phase

    dist_(R/pi Z)(delta_ell,-pi/4)<.025.              (1)

This is the phase in that metric's own bolt-anchored proper radius.
It is not an interior adjoint estimate or a seed-to-candidate phase
comparison in a common coordinate. The latter retains a proper-length
quadrature multiplied by the large physical frequency.

The Newton map used by the fixed seed projection has a different
output density from the endpoint Hessian. After the exact density
conversion, its projected derivative is invertible if and only if
one specified signed overlap sigma_ell is nonzero. We give its
convergent real resolvent representation and the exact inverse
formula. No lower bound for sigma_ell is asserted. In particular
neither (1) nor the preceding small outgoing contraction supplies it.

## 1. Retained operators and finite real construction

Use the [finite TT realization and range](finite-comparison-tt-realization-and-range.md),
including its pulled-back scalar-c domain and correlated mass column.
At a fixed ell, let T=-Box_(g_ell), D the actual scalar-c TT operator,
B its mixed row, E=60pi/G+11lambda/3, k^2=E-9lambda/4, and
Q^2=T+4lambda. All derivatives, adjoints and L2 pairings in this
section use the actual metric and dmu_ell. The regular Bessel bolt,
paired Dirac domain, clock, flux and finite coefficients remain fixed.

As in the [seed phase construction](stationary-scalar-phase-moment.md),
set V=(sqrt(80pi^2) B D^(-1/2), U_s*) and J=diag(I,-J_s).
Here V is an auxiliary vertex, unrelated to the metric radius tangent.
The complete scalar Schur operator is

    P=(T+4lambda)(T-E)-V J V*.

The source/response spaces now use omega=t^(-1/2):
P_omega consists of g+div f with omega^-1 g,omega^-1 f in L2;
V_omega has norm ||omega tau||+||omega grad tau||, with the
regular bolt realization. The finite mixed theorem gives bounded
paired maps V:Haux->P_omega and V*:V_omega->Haux. Keep

    R+=(T-E-i0)^-1, R4+=(R+-Q^-2)/(E+4lambda),
    A=V*Q^-2 V/(E+4lambda)>=0,
    L=V*R+ V/(E+4lambda), N=(I+A J)^-1,
    ||N||<beta=1.000002, ||L||<l0=.02405.           (2)

These estimates include the positive low-frequency TT block; that
block is not required to be small. The local scalar block of A is
below 2.6e-13, which is what makes its positive elimination uniform.
Conjugation gives R-, and R_PV=(R++R-)/2 is defined between the
weighted source/response spaces. No unweighted L2 principal-value
operator at the continuous-spectrum energy is used.

Let phi_ref be the actual real scalar reference field, normalized by

    U_ell phi_ref=cos(k rho_ell+delta_ref,ell)+o(1),
    U_ell tau=sqrt(8pi^2 a_ell r_ell^2) tau.

The exact regular/Jost Wronskian from the seed proof applies to the
finite potential: it has the same regular bolt and integrable outer
potential. Thus, with c0=k(E+4lambda),

    Im R+=(1/k)|phi_ref><phi_ref|,
    g=V*phi_ref, Im L=|g><g|/c0,
    ||g||^2/c0<=l0.                                (3)

The field belongs to V_omega because omega is L2 in proper radius;
the physical derivative is regular at the bolt. The weak pairings
and rank-one identity extend from compact sources by these bounds.

Put L_R=(L+L*)/2, K_R=N L_R J and

    M_R=I+A J-L_R J=(I+A J)(I-K_R),
    ||K_R||<q=.025,
    g_R=M_R^-1 g,
    a_mix=<g,J g_R>/c0,
    Psi=phi_ref+R4_PV V J g_R,
    psi=Psi/sqrt(1+a_mix^2).                        (4)

Here R4_PV=(R_PV-Q^-2)/(E+4lambda). The identity M_R*J=J M_R
shows a_mix is real, including when J has a kernel. The same exact
algebra as at the seed gives V*Psi=g_R and P Psi=0. Moreover

    ||M_R^-1||<=beta/(1-q)<1.026,
    |a_mix|<beta l0/(1-q)<.0247,
    U_ell psi=cos(k rho_ell+delta_ref,ell
                                      +arctan(a_mix))+o(1). (5)

The finite weak reconstruction upgrades psi to outer weight one
with three ordinary derivatives. The mixed transpose therefore
produces critical L2 TT forcing. The actual adjoint is

    c_psi=-D^-1 B*psi, chi_ell=Rcal_ell(psi,c_psi). (6)

The previously proved asymptotic theorem supplies its allowed
oscillatory/mass/fast expansion. Its nonzero scalar amplitude and
the physical kernel exclusion make it the unique actual adjoint
line. This construction has not changed the physical domain.

For the reference phase use dimensionless z=sqrt(lambda)rho_ell.
The [finite scalar comparison](finite-comparison-scalar-resolvent.md)
has ||q_geo,ell||_1<72000 and kbar=k/sqrt(lambda)>4e17.
The regular critical Bessel field and its bounded outgoing kernel
give, by the real rank-one construction with multiplication vertex,

    delta_ref,ell=-pi/4+arctan(a_geo),
    |a_geo|<=eta_geo/(1-eta_geo), eta_geo<2e-13.

The sign vertex is -sign(q_geo,ell). Combining this with (5)
proves (1), without a high-frequency replacement of Q^-2.

## 2. A signed moment with an explicit series remainder

For any real p in P_omega, define

    h_R(p)=V*R4_PV p,
    m_ell(p)=<phi_ref,p>+<g_R,J h_R(p)>.

Symmetry of the real inverse gives exactly

    <psi,p>=m_ell(p)/sqrt(1+a_mix^2),
    A_out(p)=i exp(i delta_ell) m_ell(p)
                       /[c0 sqrt(1+a_mix^2)].      (7)

The sign and normalization follow also from the rank-one resolvent
identity M=M_R-i|g><g|J/c0: M^-1g=g_R/(1-i a_mix).
The incoming/outgoing ratio is (1+i a_mix)/(1-i a_mix), of modulus
one. All quantities in (7) retain the full TT response and the
negative elliptic term in R4_PV.

For m>=0 take g_R^(m)=sum_(j=0)^m K_R^j N g. Then

    ||g_R-g_R^(m)|| <= beta q^(m+1)/(1-q) ||g||,
    |m_ell(p)-[<phi_ref,p>+<g_R^(m),J h_R(p)>]|
                 <=beta q^(m+1)/(1-q)||g||||h_R(p)||. (8)

For example m=9 makes the dimensionless coefficient in (8) less
than 1e-16. This bounds ONLY the series tail. It is not an error
bound for evaluating N, g, h_R, the actual kernels or a signed source.

## 3. Density-correct fixed projection test and exact inverse

Keep the [endpoint-preserving seed projection](endpoint-preserving-compatible-projection.md)
Pi_cut=I-p_cut eta0, eta0(p_cut)=1, unchanged. Its Hilbert space is

    H0: ||f||=||t^2 f||_L2(dmu0), ||p_cut||<1.006.

The fast output norm Y_(7/2) is equivalent to H0, with the usual
regular disk norm. Thus p_cut is a legitimate compact smooth output
and Pi_cut is bounded there. No higher derivative norm is inferred
from the 1.006 constant.

Write w=ell u+S, S=log[1+tau(exp(ell V_geom)-1)]. The exact
[fixed-density identity](finite-comparison-physical-tt-contacts.md) is

    Jvol=dmu_ell/dmu0=exp(2w+2S),
    mathfrakF=Jvol F, B_ell=D mathfrakF=Jvol A_ell. (9)

B_ell is the physical Newton derivative in the fixed seed output
density. A_ell is the endpoint action Hessian per dmu_ell. The
density derivative is already included in A_ell; (9) does not add
or remove its contact. In particular the range functional for B_ell is

    L_ell^0(f)=<chi_ell,Jvol^-1 f>_ell
             =integral chi_ell dot f dmu0,
    sigma_ell=L_ell^0(p_cut).                       (10)

The inverse density multiplies the source BEFORE the nonlocal York
transpose. It cannot be commuted through that transpose. The earlier
endpoint notation <chi_ell,p_cut>_ell is applicable to the fixed
projection only if its argument denotes Jvol^-1 p_cut. Equation (10)
specifies that conversion explicitly.

The operator Pi_cut B_ell:X_m->ker eta0 is an isomorphism exactly
when sigma_ell!=0. Indeed, if sigma_ell!=0, for y in ker eta0 its
unique inverse is

    (Pi_cut B_ell)^-1 y
      =(B_ell|range)^-1[y-p_cut L_ell^0(y)/sigma_ell]. (11)

The bracket is in ker L_ell^0=Ran B_ell and projects to y. A
homogeneous projected solution would have B_ell x=a p_cut, whose
L_ell^0 moment forces a=0, then x=0 by the physical kernel theorem.
Conversely, if sigma_ell=0, p_cut belongs to Ran B_ell, so its unique
preimage is a nonzero projected kernel vector. The projected range is
Ran B_ell intersect ker eta0, of codimension one in ker eta0. Thus
at zero overlap the projected kernel and cokernel both have dimension
one, although the unprojected physical kernel is still zero.

With C_ell the unquantified norm of the inverse on Ran B_ell,

    ||(Pi_cut B_ell)^-1||
       <=C_ell[1+1.006 ||L_ell^0||_(H0*)/|sigma_ell|]. (12)

This exhibits the exact missing quantitative ratio. Rescaling chi_ell
does not improve it. At the seed, choose the sign convention shared
with the unit-amplitude construction; sigma0=||t^-4 chi0||_H0>0.
This identity at ell=0 does not establish a uniform finite lower bound.

For a concrete evaluation of (10), form

    (p_tau,p_c)=Rcal_ell* (Jvol^-1 p_cut),
    p_eff=p_tau-B D^-1 p_c, p=80pi^2 p_eff.

Then the fully coupled, density-correct overlap is

    sigma_ell=m_ell(p)/[80pi^2 sqrt(1+a_mix^2)].     (13)

The source is admissible: the original p_cut is compact smooth;
York transposes retain only their permitted z^4 tails, and B D^-1
is a bounded weak-source map by the mixed energy estimate. Formula
(8) applies. A rigorous interval for the numerator excluding zero
would prove transversality. An interval containing zero proves neither
zero nor nonzero, and the projected singularity would not itself be
a physical obstruction or a failure of the unprojected range theorem.

## 4. Keep the geometric phase in a common coordinate

Let y=sqrt(lambda)rho0 and z_ell(y)=integral_0^y exp(w(v))dv.
The outer jets give the convergent exact shift

    d_ell=lim_(y->infinity)(z_ell(y)-y)
         =integral_0^infinity [exp(w(y))-1]dy.       (14)

The finite scalar proof gives |w|<3e-8 globally, y_40<177,
|w|<10t^-3 for t>=40 and dy<1.5dt/t there. Split at t=40
and 1000, and use |exp(w)-1|<=1.001|w|. Then

    |d_ell|<1.001[177(3e-8)+1.5(3e-8)log25+5e-9]
             <5.5e-6.                             (15)

The retained benchmark g=G lambda>.999(2000)^-10 gives kbar<5e17.
Thus the available upper allowance for |kbar d_ell| is 2.75e12;
this is not a measured phase shift and not a positive lower bound.
In the common seed proper coordinate the outer phase is exactly

    kbar y+delta_ell+kbar d_ell.

The existing small metric value bound therefore does not itself
certify a small phase difference modulo pi. The much smaller true
integral, cancellation, cold localization or oscillatory overlap
must be proved if used. Nor does the small asymptotic phase in (1)
bound the interior correction R4_PV VJg_R in the source-pairing norm.
The positive low-frequency block and elliptic subtraction remain in
that correction. We do not infer overlap positivity from either fact.

## Verification and next action

Freeze this argument and checker; reconcile the finite TT/range,
seed phase, fixed projection and density-contact receipts. Check the
real block factorization, rank-one elimination, phase and series
budgets, density conversion and projected inverse algebra. Independent
finite matrices must exhibit both the transverse and zero-overlap
cases and detect omission of the output density. A direct finite
resolvent sum must agree with the real-plus-rank-one formula, and
its truncation must satisfy (8). These fixtures do not evaluate the
physical signed overlap. Domain and limiting-absorption statements
remain ordinary proofs; no Lean or physical numerical solve is claimed.
Retain every run and any failure's exact source snapshot.

The first run stopped at the projected inverse fixture because its
symbolic matrix norm retained an absolute value around a rational
expression equal to zero. The failed receipt and exact sources are
retained. The corrected checker cancels each rational matrix entry
before taking that norm. The tested inverse formula is unchanged.

Next bound or evaluate (13) for the original fixed p_cut, with all
actual coordinate, density, TT and low-frequency contributions. The
small contraction supplies a fast auxiliary series, but the signed
terms and their source-pairing error remain the central missing data.
After transversality, useful full weighted response and nonlinear
preservation are still required. No S167 acceptance task closes here.

Remaining to-do list: signed fixed-projection overlap and useful full
weighted response; cold/higher bounds and nonlinear preservation;
projected solve and signed physical compatibility; exterior/state,
charge/return, opening, stability, transport, preparation, creation and audit.
