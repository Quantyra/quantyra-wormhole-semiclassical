# Actual canonical interior elimination and its original outer return

2026-10-25 HST, E003/S169. Ordinary argument before frozen execution.
Full stable semiclassical wormhole objective ACTIVE/INCOMPLETE. The closed
feasibility audit and original admission norms are unchanged.

OUTCOME: the complete coupled estimate is UNRESOLVED. This calculation
bounds the actual canonical interior compression in energy and its ENTIRE
return to the original outer graph and three phase observations. It does
not evaluate the signed complementary blocks, the finite border inverse,
the physical bolt return, or the post-elimination mass coefficient.

## 1. First unsupported step and the domain used to repair it

The full Hessian inverse bound does not pass automatically to an interior
compression: the Hessian is indefinite. For example diag(-1,1), compressed
to span(1,1), is zero. This elementary example refutes that inference,
not the actual interior operator. Nor may one impose compact support on
both York components: a canonical interior field has an exterior York
tail. We retain and bound that tail below.

Use the actual finite seed and frozen York congruence of
current-seed-full-hessian-realization.md and
current-physical-york-dual-assembly.md. Their closed form has

    H=[-a,b;b*,l], V^2=a(s)+l(xi), |b(s,xi)|^2<=.338 a(s)l(xi),
    l>=.99 E_eff, ||H||_(V->V*)<=sqrt(1.338).                  (1)

The last bound follows by squaring [-I,C;C*,I] in diagonal energy
coordinates; its square is diag(I+CC*,I+C*C). It is not positivity of H.
All normalizations below use the original energy
E^2=(a+l)/(|W_action| beta); the same factor defines its dual norm.

Let c=.2. Define W_I to be the kernel, in V, of canonical restriction
(nu,T)|_(0,c). The York-to-canonical map is continuous into physical L2
by current-coupled-schur-graph-attempt.md, section 2, so W_I is closed.
Use the actual global York return and original regular-bolt form domain.
The quantum propagators still run on the whole paired geometry, with the
original complementary domains, common regulator and finite returns.
This is a subspace of metric directions, not a new Dirac wall or state.

The exact energy elimination is H_II=H|_(W_I x W_I). The complement is
the FULL quotient V/W_I, including the transition and all matching data.
For any lift v of a quotient class, put

    u_I=P_I (H v)|_W_I, P_I=H_II^-1,
    h=v-u_I, H(h,w_I)=0 for every w_I in W_I.                  (2)

Changing the lift by w_I changes u_I by the same w_I. Thus h is
independent of that choice. The remaining quotient equations impose the
original interface matching and both original equations. We neither
discard their infinite-dimensional part nor assert it is solved by the
three phase observations. For a source f, the same elimination adds
P_I f|_W_I. This gives an exact weak block factorization; strong matching
and physical bolt traces still require their own estimates.

## 2. Exterior York tail and a bounded scalar extension

The exact kinematic inverse in current-coupled-input-energy-bound.md is

    -y''+[6/z^2+2A/(3Q)]y=2T'/(3z^2)-8T/(3z^3), A=Q''-4Q'/z,
    xi=-T-z y-z^2 y'/2, s=nu+T/2+z^2 y'/4-z y/2.             (3)

Its outer-fixed and regular-bolt conditions remain y(0)=y(z_h)=0.
If (s,xi) is in W_I, on (0,c) the unique energy-admissible homogeneous
branch is y=lambda_tail z^3 f(z), with

    f''+6f'/z=V0 f, V0=2A/(3Q), f(0)=1, f'(0)=0,
    p_s=(z^4 f+z^5 f')/4, p_xi=-(5z^4 f+z^5 f')/2,
    s=lambda_tail p_s, xi=lambda_tail p_xi.                  (4)

The other indicial branch y~z^-2 gives xi~z^-1 and is excluded by
the original energy, not by a new boundary condition. Distributional
equation (3) and local elliptic regularity on (0,c) give the same branch
for energy limits. No analyticity of an arbitrary mass/log input is used.

Let M_I be the projection of W_I onto its TT coordinate. Retain its
global kinematic constraint; it is NOT all TT fields whose local tail
looks like p_xi. For xi in M_I define E_s xi=lambda_tail p_ext, where
Delta=.0002, d=c+Delta, s_c=p_s(c), s_c'=p_s'(c), and

    p_ext=p_s,                                         0<z<=c,
    p_ext=s_c c/z+(s_c'+s_c/c)Delta t(1-t)^2,
                      t=(z-c)/Delta,                   c<z<d,
    p_ext=s_c c/z,                                     d<=z<=z_h. (5)

Values and first derivatives match at c,d. This C1 piecewise function
is in the scalar H2 form domain; second derivative jumps cause no delta
in its weak second derivative. It is regular in disk coordinates at the
actual bolt and has the original admissible outer decay. It is an
auxiliary linear form test, not a proposed C5 nonlinear correction.
Local H2 approximation across c,d proves membership without changing
the seed, the paired domain or the mass/log boundary fields.

Every W_I field has a unique representation (s,xi)=(t+E_s xi,xi),
with t scalar-form-admissible and zero on (0,c). Conversely that
representation with xi in M_I belongs to W_I: its canonical T vanishes
there by the ORIGINAL global York map, and its scalar cancels the
returned S. The bound proved below makes E_s continuous. If xi_j in
M_I converges in l, then (E_s xi_j,xi_j) converges in V and W_I is
closed, so its limit lies in W_I. Thus M_I is closed. The t space is
also a closed kernel of restriction. These are the complete Hilbert
spaces needed for the coercive signed test, including their constraints.

## 3. Executed whole-interval extension estimate

The exact scalar form and current guards give the upper bound

    a(s)<(8piG)^-1 int[6D_R(Delta4 s)^2
                          +6.0001 |grad s|^2+24.001 lambda s^2]dvol.

Indeed D_R lambda<2.946e-9, |R4/lambda+12|<1.1e-27, |J|<7lambda,
epsilon lambda<2000^-10/(60pi), and
lambda0/lambda=1+11epsilon lambda/6. These control BOTH upper
coefficients from the exact K and 2R4+48lambda0; no gradient is set zero.

Divide scalar and TT energies by the SAME factor
N_halves 4pi beta_T B^3 lambda/(8piG). With r_R=2.946e-9,
Q_plus=Q+tau z^3, Ls=Q s''+(Q'-2Q/z)s', use

    A_up=int_0^1 z^-4[6r_R(z^2 Ls/beta)^2
                       +6.0001 Q_plus z^2 s'^2/beta+24.001s^2]dz,
    L_low=.99 F_l int_0^c z^-4[
                      1.5 Q(z xi'-2xi)^2/beta+xi^2]dz.     (6)

Here s=p_ext and xi=p_xi. Q_plus=(1-z)^2 times the original positive
bracket, so extending the scalar NONNEGATIVE MAJORANT from z_h to 1
only enlarges it. The actual Q remains inside Ls and the TT integral.
This does not replace the finite-temperature operator by a double root.
The unknown interior TT energy is positive and is simply not needed
for this lower bound: l(xi_actual)>=lambda_tail^2 L_low times the factor.

For a reproducible continuum enclosure take N=80 and R=.4. On |z|<=R,
the exact coefficient triangle inequality gives Q_lower>0 and
|V0|<=Vmax<11.511. The integral equation for (4) has norm at most
alpha=Vmax R^2/14<.132, so |f|<=Fmax=1/(1-alpha)<1.152.
Writing V0=sum V_k z^k, its exact series is

    f_0=1, f_1=0,
    f_n=sum_(k=0)^(n-2) V_k f_(n-2-k)/[n(n+5)].             (7)

Cauchy's coefficient bound supplies the ENTIRE omitted analytic tail.
For j=0,1,2 let m_j=(4,12,52)_j and

    T_j=Fmax (N+6)^(j+1) m_j (c/R)^(N+1)/4.

Then |partial^j(p_s-p_s,N)|<=T_j z^(4-j), and the p_xi tail has
the bound 2T_j z^(4-j), throughout (0,c). The N+6 factor covers
the xi coefficient n+5 as well as the differentiated power n+4.
The moments are exactly sum_(k>=0)(k+1)^(j+1)/2^k.

Polynomial products in (6) are integrated exactly on (0,c). Minkowski
adds/subtracts the full energy tails using Q<.17, |Q'|<.3; the two
tail norm bounds are below 3e-22 in these normalized units. Enclose
s_c,s_c' by those same tails. On (c,d), 512 closed interval cells in
t cover EVERY point; outward interval integration encloses the complete
derivative-junction cost. On (d,1) integrate the exact Laurent polynomial

    (s_c c)^2 z^-6[6r_R(4Q-zQ')^2/beta^2
                           +6.0001 Q_plus/beta+24.001].    (8)

The frozen test targets are

    A_up<.000617, L_low>.00279,
    a(E_s xi)/l(xi)<=r=A_up/L_low<.221.                    (9)

Preview values are approximately .0006163, .002793465 and .22060;
only the outward frozen execution and ordinary inequalities certify (9).
This is a complete interval calculation of a specified extension, not
a finite spectral selection or a sampled quadrature.

## 4. Actual interior inverse and complete separated Schur entry

In coordinates (t,xi), the TT diagonal is

    d_E=l(xi)-a(E_s xi)+2 Re b(E_s xi,xi)
       >=delta l(xi), delta=1-r-2sqrt(.338 r)>0.             (10)

The mixed terms cancel in the real pairing with (-t,xi). The transformed
form is bounded on the closed product spaces just identified, and its
signed pairing is coercive with constant min(1,delta). Lax--Milgram
therefore applies there, for all dual sources. The triangular map
(t,xi)->(t+E_s xi,xi) and its transpose have squared norm at most

    C_R^2=1+r/2+sqrt(r+r^2/4).

Consequently the ACTUAL canonical interior inverse satisfies

    ||P_I||_(W_I*->W_I)<=C_R^2/delta<7.                    (11)

Both conversions are paid. The full inverse's bound one has not been
substituted for (11). No strong interior graph or bolt-value bound follows.

For every actual coupled input v, (1)-(2) now give

    E(P_I H_Iv)<=7sqrt(1.338) E(v).                         (12)

The resulting canonical field vanishes on (0,.2). Thus all local cross
terms with original outputs z<=.1 vanish by separated support, while
the nonlocal quantum return is exactly the whole entry estimated in
current-global-far-return-energy-bound.md. That argument includes all
angular and signed thermal modes, both physical paths and finite returns;
it bounds the original .05/.1 cutoff transition, not a new localization.
Contacts vanish only in THIS separated entry and remain in H_II and
the complementary blocks. Its constants times (12) give

    ||L_oI P_I H_Iv||H3             <3.712e-12 E(v),
    |quartic finite-part return|    <6.437e-25 E(v),
    |mass moment return|           <6.437e-26 E(v),
    |exact slow moment return|     <9.240e-22 E(v).           (13)

These are symmetric zero-centered enclosures of the entire interior
return, not computed signed centers. They apply simultaneously to the
mass lift, scalar quartic, exact correlated slow pair and fast inputs,
provided their ORIGINAL input maps are inserted. No phase is a fast input.
They are actual strong OUTER estimates obtained by separation, not an
invalid generic weak-to-strong embedding for the interior field.

For f=H h_lin with the established E(h_lin)<1.1e-10, the entire affine
interior term in the H3 row is <4.084e-22. For any identified input in
the earlier weighted three-jet norm N,

    E(v)<31[N(nu)+N(T)]

makes its H3 return <1.151e-10 N(v). Individual phase/feedback input
norms and subsequent finite-border amplification are still UNKNOWN.

## 5. Theorem, regulator and boundary hypotheses

The new inverse proof uses bounded forms on CLOSED Hilbert subspaces,
the explicit coercivity (10), and pays the congruence norms. It uses no
smooth-boundary inverse theorem or unproved positivity of a compression.
The Cauchy estimate extends only a one-sided, zero-free polynomial seed
ODE. Scalar form completion across an auxiliary junction is not smoothing
the actual reflected C2,1 mass cubic or changing its log domain.

The imported H is the already commonly restored fixed-seed form, with
all contacts and density corrections. We restrict that limiting form;
we do not interchange regulator removal with energy completion or compute
a determinant separately on W_I. On the separated output supports the
existing absolutely summable differentiated kernel agrees with that form
on the original core. Its energy-to-L2 bound extends the same distributional
identity continuously to the present W_I inputs, without requiring bounded
pointwise interior fields. No nonlinear action on arbitrary V fields is
being defined by this linear extension.

Re-inspected [Grubb, Theorem 2.3](https://web.math.ku.dk/~grubb/CPDE14.pdf):
identified C^tau S^d_(1,0) symbols, |s|<tau and N>dimension/2+1 are
required. Its retained seed use is output order zero; it does not supply
the missing strong quotient seminorms. [Gaarde--Grubb, Theorem 2.4 and
Proposition 2.9](https://arxiv.org/pdf/math/0703878) apply to the stated
smooth differential boundary realization and reflected logarithm, not
directly to an unidentified determinant Hessian. Their retained local
use still needs the finite kernel/contact construction. We do not apply
them to nonlinear rough backgrounds or to the artificial junction (5).

The separated proof's [McIntosh--Morris, Theorem 3.1](https://arxiv.org/pdf/1201.1818)
uses the original self-adjoint paired first-order system, its unitary
group, domain-preserving Lipschitz multiplication, bounded first
commutator and zero repeated scalar commutator. The restriction W_I is
on metric inputs and does not alter any of these Dirac hypotheses.
No new neighborhood result, independent boundary-phase C3 bound, or
separately regulated low-word action is asserted.

## 6. Full ledger, what is still missing and research decision

Retain J_m=int_0^z t^3/Q^2, n_m=-2beta QJ_m/z,
b_m=beta z^3/Q+beta Q'J_m. The reference f=2xz, alpha=0 cancels
BOTH local equations. The original chi_m (.1/.2) gives
f_cut=2xz chi_m-2F n_m Dchi_m and alpha_cut=n_m Dchi_m.
It retains zero reference bolt displacement, fixed temperature and
reference mass -xB per unit coefficient. It is not a global zero mode.
No derivative of the not-yet-quantified chi_m is assigned a value here.

The original conversion 2133<||d_chi||H3<2174 and LINEAR
|m|<1.627532e-15 debit <=3.538254568e-12. With distant shape
<1.24e-20, the near SHAPE allowance remains >6.4617454e-12.
The original scalar coefficient 20beta Fstar c4~3.333333294c4
is not erased by D-4. The slow pair must be evaluated as
B0(z^4-z^(4+e_*)) using its exact root and expm1, not separate B0 terms.

The proposed kappa_m<=5e-5 remains UNKNOWN. In the original Schur
block it requires the direct mass entry minus (13), followed by the
remaining original border/return maps. If their proved amplification
is A_o and direct/other debit is d_m, an explicit sufficient interior
allocation is

    d_m+A_o (3.712e-12) E(Phi_m)<5e-5.                    (14)

For d_m=0 only, A_o E(Phi_m)<1.3469e7 is a conservative allowance;
in N units A_o N(Phi_m)<4.344e5 suffices. Neither is an actual input
norm or a mass-column certificate. Phase observations in (13) must
also be multiplied by their actual border inverse entries.

At conditional kappa_m=5e-5 retain costs 2.14765e-12 affine,
.0024731 graph and .0195239 per independent energy increment.
The new factorization may allocate the affine interior H3 entry
4.084e-22 ONCE instead of the prior far allocation, leaving
5.3423499995916e-12 for other affine entries when A_o=1 and no
phase-to-graph amplification is present. Those are explicit conditions,
not facts about the remaining border. In general require

    Y_other+A_o 4.084e-22+Y_phase_return <=5.34235e-12,
    eta_other+.0195239 gamma_E
              +A_o 3.712e-12 J_E+eta_phase_return <.0012069. (15)

J_E is the actual energy norm of the quotient feedback INPUT lift;
gamma_E bounds the full solution energy increment by the original graph
data on a proved source-controlled solution family. They are different
unknown maps. An outer-data-only gamma_E on the whole domain is false.
With every other debit zero the sufficient ceilings are gamma_E<.0618165,
and separately A_o J_E<3.251e8. Their evident scale separation shows
that further optimization of (11)-(13) will not decide feasibility unless
the missing input/border maps are enormous. It does NOT bound those maps.

The first remaining continuum gap is quantitative control of the full
common-restored complementary quotient operator, including its near
coincidence and .1/.2 transition, original scalar finite-part/mass/slow
rows and matching reconstruction. The weak elimination (2) closes its
interior block and (13) its separated outer return. It does not identify
the quotient with an already bounded strong graph, control its interface
or bolt trace, or compute its signed retained centers. The previous
complete-spectrum Dirac tail bounds remain available for that calculation;
the strong error depends on actual U3,W4 input maps. There is no
certificate that all remaining entries fit (15).

The fast screen .028913+1e-11(1.0706e10)=.135973<.1599164
retains ONLY its two-fast admitted-family scope. The original conditions

    E_center+1e-8 alpha_embed<=1e-11,
    R_X alpha_embed*.028913+Z0_other<=.25,
    R_X alpha_embed^2*1.0706e10+Z2_other<=5e7,
    Y0+1e-8 Z0+5e-17 Z2<=1e-8, Z0+1e-8 Z2<1

remain unclosed; bare C5-to-H8log is not a bounded embedding. Compatibility,
physical matching, validity/noise and all-sector stability remain unproved.

RESEARCH DECISION: retain the construction and use this constrained
energy elimination with the existing exact Dirac transport method.
Stop refining the separated return. Sharpen the COMMON-RESTORED
COMPLEMENTARY QUOTIENT AND BORDER; its unknown strong maps and signed
entries, not a demonstrated excessive interior bound, decide feasibility.
No actual controlled solution or obstruction to construction is established.

ONE NEXT CALCULATION (PROPOSED, NOT EXECUTED): assemble and certify the
signed complementary quotient return for the affine source and ALL
mass/quartic/correlated-slow/fast columns, with original interface and
phase matching. Insert (13) as the entire interior error enclosure,
compute the input/border map norms and test (14)-(15) simultaneously.
This turn executes the interior inverse and separated Schur bound, not
that complete signed solve or its physical bolt/admission estimate.

Remaining to-do list: signed complete coupled quotient/border and physical
norm comparisons; nonlinear admission/contraction; exact compatibility;
physical matching; semiclassical validity/noise; all-sector stability;
charge/return, opening, transport, preparation, separate creation and
the full stable semiclassical wormhole objective.
