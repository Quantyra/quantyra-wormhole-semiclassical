# Complete scalar source for the accepted Green input

2026-10-28 HST, E003/S169. Ordinary argument and numerical contract before
execution. Preserve P32, the original norms, closed audit and accepted results.
The goal remains the full stable semiclassical wormhole construction.

This is the earlier scalar-source solve identified by
current-first-inverse-source-screen.md. Its quantitative contribution is
E_source <= K_scalar ||rho_s||E_s*, with K_scalar approximately
1.8608968038e-44 in ORIGINAL bolt units. It is not the independent signed
quantum solve or a complete bolt certificate.

## Actual space and complete equations

H_mu=L2((0,L),z^-4 dz), e_s=a/N, e_t=ell_t/N and N=|W_action| beta retain
their definitions. E_s is the positive scalar form completion of the
original zero-source mass/log core and regular disk bolt. The fixed paired
first-order matter domain, complementary adjoint and composed squares,
clock, temperature, spin and subtraction scale are unchanged. The actual
York inverse has Y(0)=Y(L)=0, equivalently the already accepted original
return conditions. No independent boundary phase is set to zero by an
estimate on this fixed-domain source map.

With the exact Green input t_hat of 0459e57, define

    xi=V^-1 t_hat,
    e_s(y,w)=<xi,bbar* w>H_mu/beta = g_s(t_hat)(w),
    rho_s(w)=g_s(t_hat)(w)-e_s(y_hat,w).

Existence of this source and its Riesz solution follows from the accepted
complete bound ||bbar* w||H_mu<=Kb ||w||E_s and bounded V^-1 on H_mu.
It does not assume that xi is an admissible C5 metric path. Both rows are
used below. The canonical nonlocal determinant has zero scalar/mixed row
by the retained fixed-domain Ward identity in
current-seed-full-hessian-realization.md, equations (11)-(12). The full
anomaly, moving-density, magnetic and gauge-return contacts do not vanish.
The common regulator and full angular/thermal sums remain in that identity;
this calculation is not a new proof of the independent quantum shape row.

Write d=D_R/A_r^2, epsbar=G/(60pi A_r^2), Cv=1+11epsbar beta/3,
Delta_bar=z^2[Q partial_s^2+(Q_s+2Q/z)partial_s], s=L-z, and
P=-Delta_bar+4beta. The EXACT scalar and mixed rows are

    Abar=6(Cv-d Delta_bar)P + Arem,
    bbar xi=4(Cv-d Delta_bar)(tau_bar xi)+brem xi.       (1)

Here tau_bar=(u+3beta+epsilon_seed)z^4
-epsilon_seed(5z^5-9z^6/2). The radial gradient coefficient of Arem is
-(2d+11epsbar/6)wbar+(22epsbar/3)Jbar_r, and its mass is 2wbar.
Use |wbar|<=w0=1.1e-27 beta, |Hess_bar wbar|<=w2=4.3e-25 beta^2,
|Jbar_r|<=tau_max=3.01beta. Equations (1) retain the off-shell vacuum
relation and are identities, not a zero-temperature substitution.

## Explicit trial and the first regularity implication

Use the SAME correlated exact parameters gamma,a,L,ell,nu,At as the
accepted Green input. Put Q2=2gamma s+a s^2, D2=(Q2 partial_s)',
t0=At[Q_ell(1+a s/gamma)-Q_nu(1+a s/gamma)]. The standard real Q
normalization is fixed. These profiles define a computational trial only.
All residuals below use the original finite Q and original endpoints.

For f=sum A_j Q_kj with sum A_j=0, set

    R_mu f=sum A_j/(mu-a k_j(k_j+1)) [Q_kj-Q_kmu],
    k_mu(k_mu+1)=mu/a.

Then (-D2+mu)R_mu f=f, with zero bolt flux and decay at infinity.
Define muY=4a/3, muS=4beta/L^2, tau0=tau_bar(L), c=2tau0/(3L^2),

    qY=R_muY t0,
    Y0=-2 Q2 qY_s/(3L^2),
    xi0=-2t0/3-muY qY/3,
    y0=c R_muS xi0.                                    (2)

Direct differentiation proves
-Y0_ss+(muY/Q2)Y0=-2t0_s/(3L^2),
xi0=-t0+L^2 Y0_s/2 and (-D2+muS)y0=c xi0.
The common factor Cv-d Delta_bar in (1) is retained through this scalar
solve; bounding its two occurrences separately would lose the cancellation.

Choose S=1e-12 and the C2 quintic chi_S, one through S and zero after 2S,
with |chi'|<=2/S, |chi''|<=15/S^2. Set tS=chi_S t0, Y_hat=chi_S Y0 and
y_hat=chi_S y0. The accepted input remains t_hat=chi_(1e-16)t0; its
replacement error is explicitly charged below. No physical cutoff moves.

The first domain issue is that these Q combinations individually have
logarithms. Each resolvent's coefficients sum to zero; y0 additionally
has sum A_j a k_j(k_j+1)=0. Thus y0 has the expansion const+c1 s+
O(s^2 log s), not a log source at the bolt. In regular disk radius
s proportional to rho^2, its nonsmooth term is rho^4 log rho. It lies
in the scalar second-derivative form domain and is approximable there
by regular radial core fields. Delta_bar y_hat is in H_mu. Integration
by parts therefore holds by form closure; no point-source flux is added.
Y0=O(s), Y0_s is bounded, and the York trial is in the Dirichlet energy
domain. All trials vanish before the outer endpoint, satisfying its
original homogeneous mass/log data. C2 cutoff transitions remain in
the weak fourth-order scalar equation; no third-derivative boundary term
is discarded. This proves E_s admissibility, NOT geometric C5 admission.

## Whole-profile bounds, not a selected matter-mode calculation

The positive integral for standard Q follows from
[DLMF 14.12.6](https://dlmf.nist.gov/14.12.E6) and
[14.3.10](https://dlmf.nist.gov/14.3.E10), with order zero, real k>0,
w>1. The new degrees kY and kS are positive but below one; the previous
degree>=1 software guard is NOT applied to them. For x=a s/gamma,

    Q_k(1+x)<=log(3/x)/2 (0<x<=1),
    Q_k(1+x)<=pi x^(-k-1) (x>=1),
    |Q2 partial_s Q_k|<=gamma,
    |Q2 partial_s Q_k|<=pi gamma(k+1)x^-k (x>=1).      (3)

For the last inequality differentiate the positive integral and use
partial_w log(w+sqrt(w^2-1)cosh u)<=w/(w^2-1).
The uniform flux bound follows from its monotone derivative and its
bolt limit gamma. For k>=0, Q_k<=Q_0 gives the first inequality.
The logarithmic endpoint expansion is
[DLMF 14.8.9](https://dlmf.nist.gov/14.8.E9), with that normalization.

For each profile f=sum A_j Q_kj use Minkowski to compute N(f)>=||f||L2,
M(f)>=||s f||L2(0,2S), H(f)>=||Q2 f_s||L2(0,2S), and point bounds
T(f), U(f) for |f|,|Q2 f_s| on [S,2S]. The code evaluates exactly

    N(Q_k)^2=(gamma/a) trigamma(k+1)/(2k+1),
    M(Q_k)^2<=(gamma/a)^3 [B2+pi^2 int_1^X x^(-2k) dx],
    H(Q_k)^2<=gamma^2(gamma/a)
                     [1+pi^2(k+1)^2 int_1^X x^(-2k) dx],
    B2=((log 3)^2/3+2 log 3/9+2/27)/4, X=2aS/gamma.

For k>1/2 the tail integral is bounded by 1/(2k-1), otherwise its
finite X expression is used. T,U use (3) at x=aS/gamma. These are
analytic full-support bounds, not quadrature inferred from sample points.
The finite cap moment for kS<1/2 is essential. Coefficient balls enclose
correlated exact definitions; none are replaced by independent midpoints.

## Actual York residual and return

Let zmin=L-2S, M3=1000, r=M3(2S)/(6a), so
|Q/Q2-1|<=r and Q_s<=2gamma+4aS+M3(2S)^2/2 on the cap.
The exact potential is PY=6/z^2+2(Q_ss+4Q_s/z)/(3Q).
There Q PY>=2(2a-2M3 S)/3>1; globally PY>=5/z^2, as accepted.
For RY=source(tS)-(-partial_s^2+PY)Y_hat,

 RY=chi[(2/(3L^2)-2/(3z^2))t0_s
        +(muY/Q2-PY)Y0-8t0/(3z^3)]
       -2chi't0/(3z^2)+2chi'Y0_s+chi''Y0.             (4)

Hence eY=||sqrt(Q)RY||L2 bounds the actual global York energy error.
This follows by dual testing on the entire original interval, using
Q PY>=1 only on the residual support. It retains the global boundary
return; the exact solution is not truncated at 2S. The induced H_mu
error in xi is <=(1/2+1/sqrt(5))eY.

For a reproducible bound of (4), use the following elementary majorants.
Write Qmax=(1+r)[4gamma S+4aS^2],
Dg=16/(3zmin), and
Ds=6(1+r)(2gamma+2aS)/zmin^2+2M3/3
   +8(2a+M3 S)/(3zmin)+muY M3/(6a).
Then |Q PY-muY Q/Q2|<=Dg gamma+Ds s. The three core terms are at most

 4sqrt(1+r)H(t0)/(3L zmin^2 sqrt(a)),
 [Dg gamma*2N(t0)/(3L^2 sqrt(muY))+Ds N(Y0)/sqrt(a)]/sqrt(1-r),
 8sqrt(1+r)[sqrt(2gamma)*sqrt(N(t0)M(t0))+sqrt(a)M(t0)]/(3zmin^3).

Use N(Y0)<=2H(qY)/(3L^2),
M(Y0_s)<=2[muY M(qY)+M(t0)]/(3L^2),
T(Y0_s)<=2[muY T(qY)+T(t0)]/(3L^2),
T(Y0)<=2U(qY)/(3L^2). The transition term is bounded by
sqrt(Qmax S)[4T(t0)/(3zmin^2 S)+4T(Y0_s)/S+15T(Y0)/S^2].
Finally, reconstructing xi rather than changing its basis costs

 e_rec <=N(Y0)/zmin+L M(Y0_s)/zmin^2+T(Y0)/sqrt(S).

This retains -zY and (z^2-L^2)Y_s/2 as well as chi'Y.
Set e_xi=(1/2+1/sqrt(5))eY+e_rec.

## Complete scalar residual and all error debits

Put rP=(2/3)tau_bar V^-1 tS-P y_hat. In the ORIGINAL H_mu norm,
use |tau_bar'|<=3 (verified from its coefficient sum), tau_max=3.01beta,
CE=2L+L^2 M3/(6a), CB=5L^2 M3/(6a)+2L(1+r). Then

 ||rP|| <=(2tau_max/3)e_xi+2M(xi0)/zmin^2
       +zmin^-2{CE[muS M(y0)+c M(xi0)]+CB H(y0)}
       +zmin^-2 sqrt(S){L^2 Qmax[15T(y0)/S^2+4U(y0)/(S Q2(S))]
                            +(L^2 Qsmax+2L Qmax)2T(y0)/S}.       (5)

This is the complete core/cutoff residual after the cancellation in (1).
For S0=1/sqrt(23.99), S1=sqrt(beta/5.99), S2=sqrt(beta/(6d)),
SH=sqrt(S2^2+13beta S1^2), its scalar-dual cost is
6(Cv S0+d S2)||rP||/beta. It requires only rP in H_mu, not derivatives
of rP. The weak identity moves Delta onto the scalar test function.

The remaining scalar form costs at most

 E_A=[K_A S1 G_y+2w0 S0 N(y0)/zmin^2]/beta,
 K_A=(2d+11epsbar/6)w0+(22epsbar/3)tau_max,
 G_y=sqrt(1+r)[sqrt(c N(y0)N(xi0))
                        +2sqrt(Q2(2S)/S)T(y0)]/zmin.

The full mixed remainder is bounded by (Kb_rem/beta)KV N(t0)/zmin^2.
Kb_rem includes ALL five remainder types from the accepted equation (5)
of current-complete-bolt-energy-bound.md: the small w contact in the
main coefficient, Hess w, curvature anomaly, magnetic term and T0 return.
For the curvature term we retain that proof's explicit tensor sums
instead of rounding successively to 1e12,1e14,1e19,1e20:

 Rjets=3*1e6[4400+2*4*2sqrt(2)*20+4*2], Cjets=100 Rjets,
 Kcurv=sqrt(22){100 Cjets(SH+2S1+S0)+1e6 S0},
 Kb_rem=4tau_max(11epsbar w0/36)S0+sqrt(22)d w2 S0
       +epsbar Kcurv+[20/(60pi)]4e-27 beta S0
       +(9e-28 beta)(5/sqrt(beta))S1.                  (6)

No new tensor estimate is asserted: this is the displayed bound of the
accepted proof evaluated before its power-of-ten rounding. The code
also reports the OLD 1e20(SH+S1+S0) version in the SAME complete budget.
This refinement is justified only by its direct contribution here.

The input-change charge is Ks ||tS-t_hat||H_mu, with the accepted Ks and

 ||tS-t_hat||H_mu <=At gamma^2/[sqrt(27)a^2 SP^(3/2) zmin^2], SP=1e-16.

The complete scalar residual upper report is the sum of this charge,
(5)'s scalar-dual cost, E_A and (6)'s mixed-remainder cost. All intervals
and arithmetic errors are enclosed. The analytic trial is defined by
(2), not by numerical interpolation. Later evaluation/quadrature of
b(y_hat,.) and the independent quantum solve must carry their OWN errors;
those unexecuted applications are not assigned zero error here.

Execute this bound with 512-bit outward intervals and frozen inputs.
Compare its original-output debit to the WHOLE 1.667e-9 allowance and
add it to the known partial subtotal. Neither comparison creates a free
margin or establishes the retained signed center, later inverse errors,
other physical/graph/source outputs, mass column or boundary phases.

Remaining to-do list: complete signed quantum residual and all application
errors; every original coupled/admission budget; nonlinear existence and
contraction; exact compatibility; matching/opening; semiclassical validity/
noise; all-sector stability; the full stable semiclassical wormhole objective.
