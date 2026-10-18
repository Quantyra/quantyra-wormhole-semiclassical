# The actual scalar phase and its physical source condition

2026-10-17 HST, E003/S139 supporting S074/S075. Ordinary argument
before focused verification. **PROGRESS; actual source compatibility open.**

Continue the [scalar outgoing inverse](stationary-scalar-outgoing-inverse.md),
[actual TT inverse](stationary-weighted-tt-inverse.md) and
[adjoint-range theorem](stationary-adjoint-phase.md). The source
conventions are those of the [centered seed](stationary-centered-covariant-seed.md),
[complete exterior source](stationary-centered-exterior-source.md) and
[finite-candidate residual](stationary-centered-complete-residual.md).
All results keep the original preparation interval and physical data.

The unique actual adjoint line has an explicit real construction
from the principal-value reference inverse and the bounded auxiliary
inverse. Its proper-radius phase is within 9e-12 of -pi/4 modulo pi.
In an outer logarithmic coordinate the center of that enclosure contains
an explicit, convergent seed-geometric quadrature. The actual outgoing
amplitude is exactly a nonzero known factor times the surviving signed
source moment. The retained unsigned source bounds do not decide whether
that moment is zero. The last point is a limitation of those bounds, not
a claim that the fixed retained source is undefined or freely variable.

## 1. Fixed operators, pairing and generalized-field domains

All geometry and operators below are at the actual analytic centered
seed, at its prescribed positive temperature, fixed flux and laboratory
clock, with the retained Dirac state/domain/subtraction. Put

    dmu=4pi beta r^2 dr,       T=-Box >=0,
    E=60pi/G+11lambda/3,       k^2=E-9lambda/4,
    Q=(T+4lambda)^(1/2),       c0=k(E+4lambda),
    U tau=sqrt(4pi beta r^2 sqrt(F)) tau.

U is unitary from physical scalar L2 to L2(d rho). The physical scalar
derivative becomes D_rho=partial_rho-(A_orb+2r'/r)/2. At the bolt its
domain has the regular sqrt(rho) Liouville branch and excludes the log
branch. The two scalar factors of the fourth-order operator use that
same regular realization, as in the outgoing theorem.

Write Fcal simply as V and Jcal as J in this note; this V is an
auxiliary vertex map, not a metric potential. Explicitly

    V=(sqrt(80pi^2) Bmix D^(-1/2), U_s*),
    J=diag(I,-J_s),              J*=J, ||J||<=1,
    P=(T+4lambda)(T-E)-V J V*.

Every operator is real: it commutes with complex conjugation. D is the
actual positive TT operator and its spectral inverse is the one proved
in the preceding TT theorem. No weighted bound for D^(-1/2) is assumed.

For w=B/r retain

    P_w={p=g+div f: w^-1g,w^-1f in L2(dmu)},
    V_w={regular-bolt tau: w tau,w grad tau in L2(dmu)}.

The quotient source norm and paired output norm make V:Haux->P_w and
V*:V_w->Haux bounded. The pairing with a weak source is

    <tau,g+div f>=integral conjugate(tau)g dmu
                       -integral grad(conjugate(tau)).f dmu.

Inner products are conjugate-linear in the first slot. Real fields give
the earlier real bilinear pairing. Thus each pairing below is defined
without claiming that an oscillatory generalized field belongs to L2.
The theorem gives the bounds

    R+=(T-E-i0)^-1:P_w->V_w,
    R4+=(R+-Q^-2)/(E+4lambda):P_w->V_w,
    A=V*Q^-2 V/(E+4lambda)>=0,
    L=V*R+ V/(E+4lambda),
    N=(I+A J)^-1,   ||N||<2,   ||L||<4e-12.              (1)

These include the derivative contact, actual half/full normalization,
physical volume, off-shell return and full TT block. In particular,
their constants are not bounds just for a local or high-frequency model.

Define the real principal-value inverse by the already established
operator-norm boundary values:

    R_PV=(R++R-)/2,
    R4_PV=(R_PV-Q^-2)/(E+4lambda),
    R-=(T-E+i0)^-1=conjugate(R+).                        (2)

Equations (2) define bounded maps P_w->V_w. Their pairings are symmetric.
No unweighted L2 principal-value inverse at the continuous-spectrum
energy is asserted. Each factor still has its regular bolt realization.

## 2. Unit-amplitude reference field and exact rank-one spectral density

There is a unique real regular solution phi_ref of (T-E)phi_ref=0 up to
scale. Choose its scale and sign so

    v_ref=U phi_ref=cos(k rho+delta_ref)+o(1),
    v_ref'=-k sin(k rho+delta_ref)+o(1).                 (3)

The sign makes delta_ref conventional modulo pi. Section 5 constructs
this solution directly and chooses the branch near -pi/4; it does not
assume the phase which is to be computed. Its nonzero amplitude follows
also from absence of a regular outgoing homogeneous solution.

For completeness the needed far-end scalar argument only uses the
integrable potential V_L-9lambda/4 on every tail [rho0,infinity).
The Volterra equation with kernel sin(k(s-rho))/k constructs the Jost
solution j_+(rho)=exp(i k rho)+o(1), j_+'=i k exp(i k rho)+o(1).
Iteration is controlled by the exponential of the L1 tail norm/k.
Continue this solution to the bolt through the ordinary equation.
Since the regular field is real, it has conjugate incoming/outgoing
coefficients; its nonzero amplitude permits normalization (3).

With theta=k rho+delta_ref, the Wronskian is

    W(v_ref,j_+)=i k exp(-i delta_ref).

The derivative jump -1 for -partial_rho^2 therefore gives

    U R+ U^-1(rho,sigma)
     =-v_ref(rho_<) j_+(rho_>)/W(v_ref,j_+),
    Im R+=(1/k)|phi_ref><phi_ref|.                       (4)

For example its far-output coefficient against a scalar source p is
i exp(i delta_ref)<phi_ref,p>/k. The rank-one identity follows either
by taking the imaginary part of the displayed Green kernel, or from
j_+ and its conjugate and their Wronskian. These are kernels with
respect to d rho after U, or with respect to dmu before U. The factor
4pi beta is inside U and has not been dropped.

The field phi_ref belongs to V_w. Away from the bolt both v_ref and
D_rho v_ref are bounded, while w is L2(d rho). At the bolt phi_ref is
a smooth even scalar and its physical derivative vanishes there. Thus

    g=V*phi_ref belongs to Haux,
    Im L=|g><g|/c0,
    ||g||^2/c0=||Im L||<=||L||<4e-12.                   (5)

The identity extends from compact scalar tests to P_w by the paired
weighted bounds. An equivalent proof differentiates the rank-one
kernel at its input/output vertices; the identity contact is real and
therefore contributes nothing to Im L.

Far-output coefficient formulas also extend to these weak sources. In
Liouville coordinates w^-1g_tilde,w^-1f_tilde are L2, hence g_tilde and
f_tilde are L1 by Cauchy--Schwarz. The Green kernel and its input
derivative are bounded. Split the input integral at a large fixed
radius, use the Jost asymptotic on its compact part and then its L1 tail.
Near the bolt regular phi_ref and its physical derivative are bounded.
The elliptic term Q^-2p is H1 at infinity and has vanishing Liouville
value there. Thus the outgoing VALUE amplitude exists in this argument;
pointwise derivative asymptotics for arbitrary rough f are unnecessary.
Smooth retained sources have the stronger asymptotics from the accepted
weighted elliptic bootstrap.

## 3. Exact real construction of the actual phase line

Put L_R=Re L=(L+L*)/2 and

    M_R=I+A J-L_R J,
    K_R=N L_R J,      rho_mix=8e-12.

Then M_R=(I+A J)(I-K_R), ||K_R||<rho_mix and

    ||M_R^-1||<=2/(1-rho_mix),
    M_R* J=J M_R,
    (J M_R^-1)*=J M_R^-1.                              (6)

The last identity uses only invertibility and selfadjointness of J; it
does not divide by J or assume that J has no zero eigenspaces. Define

    g_R=M_R^-1 g,
    a_mix=<g,J g_R>/c0 in R,
    Psi=phi_ref+R4_PV V J g_R.                          (7)

Since V*R4_PV V=L_R-A, applying V* to (7) gives
V*Psi=g+(I-M_R)g_R=g_R. Applying the fourth-order reference operator
then gives P Psi=0 weakly and with the regular bolt conditions.

The far-output coefficient of the real inverse (2) is
-sin(theta)<phi_ref,p>/[k(E+4lambda)]. Consequently

    U Psi=cos(theta)-a_mix sin(theta)+o(1),
    psi=Psi/sqrt(1+a_mix^2),
    delta_full=delta_ref+arctan(a_mix) modulo pi,
    U psi=cos(k rho+delta_full)+o(1),
    |a_mix|<8e-12/(1-8e-12)<8.1e-12.                   (8)

The ordinary global domain argument is the same one used for the
outgoing reconstruction, now with the smooth reference field added.
First psi is in V_w by (7). Weighted scalar ellipticity gives
psi in z^(1/2)H_b^3, with regular disk behavior at the bolt. The actual
mixed adjoint has order at most two and gains two outer powers, and
the sharp York transpose keeps its allowed z^4 tails. Therefore

    Bmix*psi in z^(5/2)H_b^1 + retained z^4 tails
              is contained in physical L2,
    c_psi=-D_crit^-1 Bmix*psi,
    chi=Rcal(psi,c_psi).                               (9)

This is the actual coupled null field with vanishing TT source
coefficient. The critical TT inverse, not a new extension, appears in
(9). Its expansion may contain the genuine mass cubic. The previously
proved adjoint asymptotic theorem upgrades the field to the allowed
oscillatory/mass/fast expansion. Since its scalar amplitude in (8) is
nonzero and the physical kernel is zero, chi spans the one actual
mass-orthogonal adjoint line.

Our normalization is UNIT LIOUVILLE COSINE amplitude. If the earlier
adjoint convention instead makes its disk leading scalar coefficient
have length one, multiply chi by the corresponding positive geometric
constant. The phase and zero/nonzero moment condition do not change;
the amplitude identity below uses precisely the normalization (8).

There is a direct incoming/outgoing cross-check. Since

    M=I+A J-L J=M_R-i |g><g|J/c0,
    M^-1g=g_R/(1-i a_mix),

the unit incident construction

    tau_inc=phi_ref+R4+ V J M^-1g=Psi/(1-i a_mix)

has incoming coefficient exp(-i delta_ref)/2 and outgoing coefficient
exp(i delta_ref)(1+i a_mix)/(2(1-i a_mix)). Thus

    S_full/S_ref=(1+i a_mix)/(1-i a_mix),
    exp(-i arctan(a_mix))tau_inc=psi.                  (10)

This proves both the phase sign and the absence of a spurious factor
two. The negative elliptic piece of R4+ remains in M_R throughout.

## 4. The exact moment is the outgoing amplitude

For real p in P_w set

    m0=<phi_ref,p>,       h_R=V*R4_PV p,
    m_num(p)=m0+<g_R,J h_R>.

Symmetry of R4_PV and J gives the exact convergent weak pairing

    <psi,p>=m_num(p)/sqrt(1+a_mix^2).                    (11)

The outgoing auxiliary source is
h_+=V*R4+p=h_R+i g m0/c0. If y=M^-1h_+, scalar rank-one elimination
gives

    m0+<g,Jy>=m_num(p)/(1-i a_mix).

Indeed substitute M=M_R-i |g><g|J/c0, then pair with g,J; the identity
(J M_R^-1)*=J M_R^-1 converts the cross term to <g_R,Jh_R>.
The outgoing Liouville coefficient of
tau=R4+(p+VJy) is therefore exactly

    A_out(p)=i exp(i delta_full)<psi,p>/[k(E+4lambda)].   (12)

For the complete coupled source, let (f_tau,f_c)=Rcal* f and

    f_eff=f_tau-Bmix D_crit^-1 f_c,
    p=80pi^2 f_eff.

For the fast retained sources, f_c is physical L2 and f_eff is P_w.
For example f_tau in Y_(7/2) embeds into w L2 near infinity, and
Bmix D_crit^-1 f_c=(80pi^2)^(-1/2)V_T D^(-1/2)f_c
has the already bounded vertex form. The precise source pairing is

    <chi,f>=<psi,f_eff>,
    A_out= i exp(i delta_full) 80pi^2
                   <chi,f>/[k(E+4lambda)].             (13)

Thus the unique fast-range condition is exactly A_out=0, equivalently
m_num(80pi^2 f_eff)=0. Necessity follows from (12). Sufficiency also
follows from the already established closed-range/cokernel theorem:
the unique compatible fast solution is outgoing and must equal this
outgoing construction. No extra growing elliptic or TT branch appears.

The construction supplies a convergent certificate representation:

    M_R^-1=sum_(n=0)^infinity K_R^n N,
    ||M_R^-1-sum_(n=0)^m K_R^n N||
                  <=2 rho_mix^(m+1)/(1-rho_mix).        (14)

The resulting error in a_mix is at most
8e-12 rho_mix^(m+1)/(1-rho_mix); the error in m_num is at most
2||g||||h_R|| rho_mix^(m+1)/(1-rho_mix). Formula (14) certifies the
series truncation if its actual signed terms are evaluated. It does
not certify evaluation errors in g,h_R,N or the source. In particular
the potentially large positive TT block inside N is not discarded.

## 5. Geometric reference phase and its proper-radius enclosure

The accepted Liouville theorem gives

    U T U^-1=H_B0+9lambda/4+U_geo,
    H_B0=-partial_rho^2-1/(4rho^2),
    U_geo bounded,    L_geo=||U_geo||1<6000sqrt(lambda),
    eta_B=L_geo/k<2e-14.                               (15)

H_B0 has the regular critical Friedrichs realization. Its normalized
regular field is phi_B=sqrt(pi k rho/2) J0(k rho). Nicholson's estimate
gives |phi_B|<=1 and its asymptotic is cos(k rho-pi/4).
Use on L2(d rho) the bounded multiplication vertex

    V_geo=|U_geo|^(1/2),      J_geo=-sign(U_geo),
    L_B=V_geo R_B0+ V_geo,    g_B=V_geo phi_B,
    M_B,R=I-Re(L_B)J_geo.

The same bounded Bessel kernel used in the reference proof yields
||L_B||<=eta_B and ||g_B||^2<=L_geo. This is a special case of the
rank-one argument of sections 2-3, with no A block and with c0 replaced
by k. Define

    a_geo=<g_B,J_geo M_B,R^-1g_B>/k in R,
    v_geo=phi_B+Re(R_B0+)V_geo J_geo M_B,R^-1g_B.

Then (H_B0+U_geo-k^2)v_geo=0, it has the regular bolt domain, and

    v_geo=cos(k rho-pi/4)-a_geo sin(k rho-pi/4)+o(1),
    v_ref=v_geo/sqrt(1+a_geo^2),
    delta_ref=-pi/4+arctan(a_geo),
    |a_geo|<=eta_B/(1-eta_B)<2.1e-14.                  (16)

Here the bounded-potential equation and its regular endpoint ensure
the ordinary smooth even scalar at the bolt, not a new boundary
parameter. The L1 tail argument in section 2 gives the derivative
asymptotic (3). Equations (8) and (16) prove

    delta_full=-pi/4+arctan(a_geo)+arctan(a_mix),
    |delta_full+pi/4|<9e-12 modulo pi.                  (17)

Both arctangent formulas are exact, including the full TT and return
terms. The signed perturbation in the Bessel construction is NEGATIVE
sign(U_geo), because the equation is H_B0-k^2+U_geo.

## 6. Actual logarithmic coordinate and geometric quadrature

A small proper-radius phase is not a small phase relative to an
arbitrarily anchored log coordinate. Define the exact seed constant

    ell_inf=lim_(R->infinity)[integral_(r_h)^R dr/sqrt(F)
                                    -log(R/B)/sqrt(lambda)]
     =rho(2B)-log(2)/sqrt(lambda)
         +integral_(2B)^infinity
                     [1/sqrt(F)-1/(sqrt(lambda)r)]dr.     (18)

The explicit rational seed polynomial and its lower bound
F>=lambda r^2/2 for r>=2B prove convergence. In fact |F-lambda r^2|
<=3lambda B^2 there, so the absolute remainder beyond R>=2B in (18)
is less than 3B^2/[sqrt(lambda)R^2]. The integral at the bolt is
regularized by r=r_h+s^2; its integrand 2s/sqrt(F(r_h+s^2)) has the
finite limit sqrt(2/kappa). These supply a concrete quadrature, using
the exact retained coefficients, rather than a formal integration
constant. More asymptotic seed terms can be subtracted if desired.

For u=B/r, write nu=k/sqrt(lambda). Then

    U psi=cos(nu log u-Phi_u)+o(1),
    Phi_u=k ell_inf+delta_full modulo pi.               (19)

If the fixed outer coordinate is z=c_u u[1+O(u^2)] with c_u>0, then

    Phi_z=k ell_inf+nu log c_u+delta_full modulo pi,
    dist_(R/pi Z)(Phi_z,
               k ell_inf+nu log c_u-pi/4)<9e-12.         (20)

For the explicit disk x=-log xi=kappa integral_r^infinity dr/F,
c_u=kappa/(lambda B). For z=c_z x, multiply this c_u by c_z.
The earlier global weighted note only requires z proportional to x;
one must retain that specified proportionality when quoting its phase.
The leading amplitude of the physical scalar itself, with z=c_u B/r,
is (4pi beta sqrt(lambda))^(-1/2)(c_u B)^(-3/2) z^(3/2).
It is positive and therefore does not change (20).

Equation (20) is a rigorous phase enclosure about an explicit geometric
function of the actual seed, NOT a decimal evaluation of Phi_z. The
seed radii are the exact homogeneous roots, d=B^2/A^2-1-6lambda B^2,
and the prescribed temperature is retained. The preparation records
state z_i in [1/4,7/25], with kappa=z_i nu_preparation; no new value of
z_i may be selected to tune (20). The formula holds for each prescribed
member of that range. Already nu>4e17, so a numerical phase requires
ell_inf with error small compared with 1/k, as well as comparable
control of nu log c_u. Coarse allowed parameter intervals do not give
a narrow numerical phase modulo pi.

## 7. What is actually specified about the preparation source

There are two distinct source questions, which must not be conflated.
At the analytic seed, the full fixed-volume negative gradient is the
inhomogeneity for its linear correction. At the finite prepared
candidate, its full residual is a different inhomogeneity. The present
phase belongs to the analytic-seed Hessian. To use a candidate residual
in (13), express it in the same fixed disk and seed-volume output,
then apply the same seed transpose Rcal*. The candidate's own moving
volume gradient or candidate-base adjoint must not be silently used.

The seed source is defined by the full common determinant, with the
exact action split Gamma_rest=Gamma_Dirac-Gamma_ref. The exterior-source
note defines the canonical massive orbit trace

    Theta(t)=B^2 tau_hat,mass(Bt),  t=r/B,  v=r_h/B,
    Bfun(t)=F(t)/t^2,
    I(t)=integral_v^t Bfun'(s)Theta(s)ds,
    J_trace(t)=integral_v^t F(s)Theta'(s)/s^2 ds,
    B^2 p2,mass=Theta/(2t^2)-J_trace/(2F),
    B^2 rho2,mass=-Theta/(2t^2)-J_trace/(2F).             (21)

HH bolt regularity fixes the integration constant. The local Weyl,
LLL, reference and anomaly contributions are retained explicitly.
The full angular pressure follows from the full anomaly trace.
The actual proved information on the nonlocal trace includes

    |Theta'(t)|<=1e15 q^2,  |Theta(2)|<=1e11 q^2,
    Theta(t)=L_trace t+R_trace(t),
    L_trace=-a_seed x/(5pi),
    |R_trace(t)|<=4e15 q^2(1+log t), t>=4.               (22)

These are anchor/gradient/tail bounds, not signed values of Theta on
the entire radial interval. The determinant definition determines
those values in principle. It is not evaluated by replacing R_trace
with zero or one of its upper bounds.

The four-dimensional conformal trace is a useful simplification:
the exact full conformal gradient at the seed is -E_tr/(8pi G), with
E_tr the local curvature/magnetic polynomial in the outgoing theorem.
Hence the NEGATIVE seed gradient has f_seed,tau=+E_tr/(8pi G).
But its effective scalar source is still

    f_seed,eff=E_tr/(8pi G)-Bmix D_crit^-1 f_seed,c.       (23)

The TT source f_seed,c contains the actual state-dependent stress
through (21) and the physical York transpose. Equations (22) do not
evaluate the second term of (23), and one cannot compute the coupled
moment by integrating only the anomaly polynomial.

There is a useful exact cancellation in the REFERENCE projection of
the local Newton source. Set epsilon=G/(60pi),
a_tr=1+11epsilon lambda/3, so E=a_tr/epsilon. Write J_Ric for the
trace-free Ricci tensor, distinct from the auxiliary sign operator J,
and C_W for the Weyl tensor. The exact trace identity is

    E_tr=(a_tr-epsilon T)Z
          +epsilon[(7/12)C_W^2+(11/6)J_Ric^2
                                    -20M_em-11Z^2/72].

Since (T-E)phi_ref=0, the first term pairs to zero. This integration
by parts is valid on the entire actual exterior: at the bolt both
phi_ref and Z are smooth even scalars and the radial Green coefficient
r^2F vanishes; at infinity Z=O(r^-5), Z_r=O(r^-6),
phi_ref=O(r^-3/2), phi_ref,r=O(r^-5/2), so the boundary product
r^2F(phi_ref Z_r-phi_ref,r Z) is O(r^-7/2) and vanishes.
The bulk integrals converge absolutely. Hence, exactly,

    <phi_ref,80pi^2 f_seed,tau>
     =(1/6)<phi_ref,(7/12)C_W^2+(11/6)J_Ric^2
                                    -20M_em-11Z^2/72>.  (23a)

Here 80pi^2 epsilon/(8pi G)=1/6, and M_em=q^2/(2r^4)
is the retained physical magnetic invariant. Thus the pointwise small
d-dependent linear trace contact is not a dominant term in this
reference projection: its whole (a_tr-epsilon T)Z combination cancels.
Equation (23a) concerns phi_ref and the negative seed gradient only.
It neither replaces phi_ref by the full adjoint psi nor removes
-Bmix D_crit^-1 f_seed,c or the second term of (27). It is not a
zero/nonzero evaluation of the complete moment or the candidate source.

The first-order prepared correction is also specified exactly in
terms of that same source. With F0 the seed, R=r_h and rho,p the full
sphere-integrated Gamma_rest pressures, it is

    a_lapse(r)=-G integral_r^infinity (rho+p)/(tF0)dt,
    D_hor=F0''(R)/2+kappa/R,
    s_R=[G rho(R)/R-kappa a_lapse(R)]/D_hor,
    f_metric(r)=-[2kappa R s_R+2G integral_R^r rho(t)dt]/r. (24)

The common-disk conformal and radius tangents u_geom,V_geom then define
the finite candidate family, writing s rather than lambda for its
bookkeeping parameter,

    S_s=S0+s Gamma_rest,
    DE0(g0)g1+J_rest(g0)=0,
    sigma_s=sigma0+s u_geom,
    r_s-B=(r0-B)exp(s V_geom), 0<=s<=1,
    R_s=E0(g_s)+sJ_rest(g_s)
       =[E0(g_s)-sDE0(g0)g1]
                    +s[J_rest(g_s)-J_rest(g0)].          (25)

The complete-residual theorem bounds both brackets with the fixed-disk
Jacobians and outer mass cancellation. It gives absolute profile and
constraint-image bounds, not their signed transform against chi.
The accepted adjoint note proves the full residual is O(z^4 log^N z),
smooth at the positive-temperature bolt, and lies in Y_(7/2).
Thus the actual moments are well defined. Neither this admissibility
nor the small equation defect proves their values or their signs.

## 8. Why the available absolute estimates cannot decide zero

The map f -> <chi,f> is a nonzero continuous real linear functional on
the admitted fast source space. Choose a smooth compact scalar bump
supported where psi is nonzero and take f_c=0. Its moment is nonzero;
changing its sign changes the moment's sign. A suitable linear
combination of two such bumps has exactly zero moment. Scaling makes
all these examples obey arbitrarily small absolute source norms and
the same regularity/tail class. The Ward witness in the actual adjoint
note completes compact reduced source pairs to conserved spherical
tensors with the permitted tails. Thus covariance supplies no missing
zero-moment identity. This is a logical statement about the information
content of unsigned estimates; these test sources are not substituted
for the actual retained Dirac source.

Even a positive smooth source need not have a comfortably sized
oscillatory moment. In the regular Bessel model the positive smooth
radial source with Liouville image sqrt(rho)exp(-b rho^2), b>0, has

    integral phi_B(rho)sqrt(rho)exp(-b rho^2)d rho
         =sqrt(pi k/2)/(2b) exp(-k^2/(4b)).              (26)

This follows from the two-dimensional Gaussian Fourier transform,
or the [DLMF Gaussian-Bessel integral](https://dlmf.nist.gov/10.22.E51) at
order zero. It is a model example, not an asserted asymptotic for the
actual preparation source. It shows why norm-small errors can be much
larger than a nonzero high-frequency moment. No positivity principle
turns (22) into a sign certificate at k/sqrt(lambda)>4e17.

Likewise (17) controls an asymptotic PHASE. It does not alone bound the
interior field difference psi-phi_ref by 9e-12 in a source-pairing norm.
The real principal-value construction includes the elliptic subtraction
and the possibly large low-frequency positive TT block. An absolute
bound on this difference, even if derived, would still have to be
compared with a certified signed leading moment before deciding sign.

## 9. Exact remaining calculation and the scope of a decision

The strongest justified source reduction is the signed real scalar

    m_num(80pi^2 f_eff)
     =<phi_ref,80pi^2 f_eff>
          +<M_R^-1g, J V*R4_PV(80pi^2 f_eff)>.          (27)

For the actual selected source, one must evaluate (27), retaining the
actual Theta/state contribution, fixed-density/York source conversion,
TT inverse response and both real resolvent pieces. A finite numerical
enclosure excluding zero would prove nonzero. Exact zero requires an
analytic identity or a separately justified exact compatibility
construction; an interval merely containing zero proves neither zero
nor nonzero. The series remainder (14) is available but is only one
part of the total error. All other signed evaluations and truncation
errors must be enclosed below the size relevant to the decision.

The explicit geometric phase formula (20) and the exact moment and
amplitude identities (11)-(13) are complete linear theorems at the
seed. The currently retained signed source data do not evaluate (27)
or bound it away from zero. This theorem therefore does not assert a
fast physical solution or an obstruction for the actual preparation.

A certified nonzero value would obstruct the particular fast LINEAR
seed equation with that fixed source. It would not by itself rule out
all nonlinear exteriors at nearby comparisons. A certified zero would
solve this one linear compatibility condition; quantitative compatible
fast inversion and nonlinear existence/control would still be needed.
An outgoing response with nonzero (12) is not a real fast branch or a
finite-mass error estimate for the order-reduced candidate.

## 10. Verification scope

The [focused checker](../checks/stationary_scalar_phase_moment.py),
[immutable receipt](../evidence/stationary-scalar-phase-moment.json) and
[raw log](../evidence/stationary-scalar-phase-moment.log) retain the
executed sources and finite exact checks. They audit rank-one/block
algebra, signs, normalization, phase budgets, source cancellation and
coordinate factors. Generalized-field domains, limiting absorption,
global adjoint identification and source admissibility remain ordinary
analytic arguments. No numerical value of the retained signed moment
or of the logarithmic phase is certified by those algebra checks.

The phase formula and source reduction do not close S139.1-S139.2 or
S074.2-S074.6. A nonzero linear moment would concern that fixed linear
source, not every nonlinear exterior; zero would still require nonlinear
control. The complete E003/S074-S080 objective remains active.

Remaining to-do list: evaluate the actual signed preparation moment or
prove its exact cancellation; obtain the compatible fast weighted
inverse bound and nonlinear control, or justify controlled spatial
reduction with observable errors; match exterior/state and full charge,
transfer S138, and retain driven stress, radius, controller, capture,
tail, spatial area and every S074-S080 gate.
