# Source-adapted first solve: full source error and the fractional flow gate

2026-10-28 HST, E003/S169. Ordinary argument and contract before execution.
The previous turn was progress: 42679b1 rejected every coefficient in the
old one-dimensional correction range on the actual source. Preserve that
result, d86adb2, P32, the closed audit and the full stable-wormhole objective.
This attempts the complete first reference solve. It is not a completed
inverse certificate if either finite-domain return or flow norm is unknown.

## Exact target and correction of the compression shorthand

Use e_t=ell_t/N on the original E_t and g_hat=b(y_hat,.)/N, with the actual
scalar trial of 5967b97. The target is an explicit x_hat satisfying

    ||g_hat-e_0(x_hat,.)||E_t* <9.955e32, e_0=ell_0/N.       (1)

A0 is the E_t Riesz representative of this COMPLETE form. The prior
checkpoint's shorthand A0=J* A(iX)J is incomplete: that compression
describes its geodesic nonlocal contribution; the original local affine
action, LLL, a0/a2 contacts and York pullback must also be present, followed
by the Riesz identification. No result there used the shortened equality
numerically. Neither the full A0 nor its compressed part is inverted by
inverting a scalar flow multiplier and compressing afterward.

Keep the actual paired first-order mass/log realization and both squares,
fixed clock/temperature/scale, full angular and thermal sums, common
regulator, coincidence contacts, original cutoffs and regular disk bolt.
The geodesic reference uses the actual C2,1 even double and all its returns.
Its .0022 difference from the actual quantum form is NOT set to zero;
the already frozen later ordered applications must account for it.

## Construct a source-adapted trial without replacing the physical domain

Set x=(L-z)/gamma. The auxiliary comparison H=-partial_x[(2x+a*x^2)
partial_x] is the nonnegative radial Friedrichs operator on [0,infinity),
with regular bolt flux and decay at infinity. It is used ONLY to define
analytic trial profiles. Write R_mu=(H+mu)^-1, muY=4a/3,
muS=4beta/L^2, lambda_l=a*ell*(ell+1), lambda_n=1/(d L^2),
c=2tau0/(3L^2), Cv=1+11 epsbar beta/3. These are the same correlated
parameters and scalar y0 already accepted. Its exact multiplier is

 y0=-(2c At/3)(lambda_n-lambda_l)(H+2a)
                  R_muS R_muY R_lambda_l R_lambda_n delta_b. (2)

All scalar profiles, including the large-degree one, are retained.
The bolt distribution has -lim q2 (Q_k)'=1; no new physical source is
introduced. These are auxiliary geometric profiles, not selected matter modes.

Freeze z=L and the non-derivative background phi=-log L only to define a
classical computational reference. Set V(p)=-3/2+a/(p+2a). Directly
differentiating the exact Einstein, magnetic and R2 generator gives the
following NORMALIZED energy symbol:

 K(p)=-[HE(p)+HM(p)+HR(p)]/beta,
 HE=[3p/2+2pV-a+(1+2V)^2/C]/L^2
                            +12lambda0 B^2(1+V)^2/(C L^4),
 HM=-4Z V^2/C, Z=pi q_ang^2/(g_R^2 B^2),
 HR=alpha_R{16V^2-12C p+C^2[16a^2 V^2+16a V D+2D^2+12a p]},
 D=p(2V+3), alpha_R=8pi c_R/(C B^2).                   (3)

The gamma terms, finite Q, original York return, Weyl and canonical contacts
are NOT asserted to vanish. Their full difference belongs to E_return in
the complete budget below. Polynomial division gives

 K(p)=A p+B+(r1 p+r0)/(p+2a)^2.                       (4)

The code retains the complete coefficients and tests A,B>0. Define
muT=B/A, Kx=4tau0/(beta L^4 A), and the analytic first profile

 xi0=Kx (Cv+d L^2 H)R_muT y0,
 v0=V(H)xi0=Kx c At(lambda_n-lambda_l)(Cv+d L^2 H)
                            R_muT R_muS R_lambda_l R_lambda_n delta_b.

Before norms, the identities V(H)=-(3/2)(H+muY)/(H+2a) and
lambda_n=1/(d L^2) cancel both York factors and the large auxiliary mass:

 v0=At K0 [1+(Cv-1)/(d L^2(H+lambda_n))]
                                  R_muT R_muS R_lambda_l delta_b,
 K0=Kx c (1-d L^2 lambda_l).                           (5)

The second term is retained, not assigned zero. No high-degree or tail
parameter is optimized. This cancellation is specific to the complete
source composition; it is not freedom to choose quantum stress.

To define an ACTUAL paired-domain field, put r=R_muY v0,
chi=-(2/3)partial_x r and

 J=z^5 chi(x)/(gamma L^5),
 I(z)=int_0^z (Q'/t^4)'J(t)dt,
 x_hat=(QJ)'-z^4 I,
 Ux_hat=x_hat/2+QJ/z, Vx_hat=-x_hat-(QJ)'/2.             (6)

Here Q is the exact finite seed on [0,L], with its original root and slope.
Both York equations and J(0)=0, (Q^2 J'/z^4)(L)=0 hold identically, as for
the accepted global trial. This introduces no cutoff or wall. The field
starts at z^4 at the paired outer seam. Its actual bolt regularity is
checked below; (6) is not silently identified with xi0 on the model.
The source y_hat retains its ORIGINAL cutoff at S=1e-12.

## Complete actual-source error, in the original dual norm

Write T(z)=tau_bar/z^4=u+3beta+eps-5eps*z+(9eps/2)z^2,
T0=tau0/L^4. The exact mixed source density in dz is

    j_hat=4T(z)[Cv y_hat-d Delta_bar y_hat]/beta+j_rem.

This uses the full mixed Ward identity and all its retained contacts.
The canonical determinant has no mixed row in this fixed chart only.
Compare with j_ref=4T0[Cv y0-d L^2 Delta2 y0]/beta on the actual interval,
where Delta2=(Q2 partial_s)' and s=L-z. This is a SOURCE comparison, not
a transfer of the quantum response to the quadratic geometry.

The actual TT lower form gives ||xi||H_mu<=Cxi||xi||Et with Cxi=1.01:
ell_t>=(1-epsilon).99 E_eff, E_eff/N>=F_l||xi||H_mu^2 and
F_l>1-1.179e-8. Thus also ||xi||L2(dz)<=Cxi||xi||Et.

Compute exact whole-model norms Ny0, Ng0, Nd0 from the full Gram matrix

 int_0^infinity Q_k Q_j ds=(gamma/a)
       [digamma(k+1)-digamma(j+1)]/[k(k+1)-j(j+1)],
 diagonal=(gamma/a)trigamma(k+1)/(2k+1).

For y0=sum c_i Q_ki, Ny0^2=c^T G c,
Ng0^2=-c^T G diag(lambda_i)c and
Nd0^2=(diag(lambda_i)c)^T G diag(lambda_i)c.
The exact cancellations sum c_i=sum lambda_i c_i=0 justify integration
by parts and remove the bolt flux; independent interval endpoints are
not substituted for these algebraic identities. These are continuous
full-profile integrals, not a finite quantum-mode approximation.

Reuse M3=1000, S, zmin=L-2S, r=M3*2S/(6a), Qmax, Qsmax and the accepted
whole-profile M(y), H(y), T(y), U(y). Let Mdy=sum |c_i|lambda_i M(Q_ki).
For f=y0 or Delta2 y0, its complete tail beyond S satisfies

 Tail(f)<=pi sqrt(S) sum |f_i| (aS/gamma)^(-ki-1)/sqrt(2ki+1).

On the cutoff support |Q/Q2-1|<=r. The original norms obey

 Ny<=Ny0/zmin^2,
 Ng<=sqrt(1+r)/zmin [Ng0+2sqrt(Q2(2S)/S)T(y)],
 Nd<=(1+r)Nd0+[5M3/(6a)+2(1+r)/zmin]H(y)+Cut,
 Cut=sqrt(S){4(1+r)U(y)/S
                 +[15Qmax/S^2+2(Qsmax+2Qmax/zmin)/S]T(y)}.

Nd bounds ||Delta_bar y_hat||H_mu; its z^2 cancels the H_mu weight.
The accepted complete scalar row consequently gives the ACTUAL norm bound

 Es(y_hat)^2 <= {6[d Nd^2+(Cv+4d beta)Ng^2+4beta Cv Ny^2]
                         +KA Ng^2+2w0 Ny^2}/beta,
 KA=(2d+11epsbar/6)w0+(22epsbar/3)(3.01beta).

This evaluates the existing trial norm; it does not change or optimize
the scalar solution. Its explicit role is the complete mixed-source debit
Erem=Cxi Kb_rem Es(y_hat)/beta, using the already accepted full Kb_rem.
Compare with the earlier coarse norm Ks||t_hat||+rho_s in the SAME debit.

Put CE=2L+L^2 M3/(6a), CB=5L^2 M3/(6a)+2L(1+r). Then

 Ndiff<=CE Mdy+CB H(y)+Cut+L^2 Tail(Delta2 y0)

bounds ||Delta_bar y_hat-L^2 Delta2 y0||L2 on the entire actual interval.
The source cutoff and the region outside it are both charged. With
|T(z)-T0|<=14|eps|, the COMPLETE source error is

 Esource=Erem+(4Cxi/beta){T0[Cv Tail(y0)+d Ndiff]
                                    +14|eps|(Cv Ny+d Nd)}.           (7)

This proves ||g_hat-j_ref||Et*<=Esource, including the original cutoff,
all complete mixed contacts and source rounding. It supplies no independent
quantum column or geometric-admission estimate.

## Actual retained geodesic response: a strong fractional estimate

Let P=iX on the ACTUAL complete Liouville flow. Define
N_9/2(v)=(int_SM | |P|^(9/2) lift(v)|^2/(2pi))^(1/2), with dz and
the original physical-half convention. The flow is the C2,1 paired one;
its strong unitary group supplies a self-adjoint P. Do not replace it by
a Laplacian multiplier. For g(s)=1-s^2/2+s^4/24-cos s,

 0<=g(s)<=min(s^6/720,s^4/24)<=c_g |s|^(9/2),
 c_g=720^(-1/4)24^(-3/4).

The positive Schwinger integral, [DLMF 10.32.10](https://dlmf.nist.gov/10.32.E10),
gives, for mu>2|nu|,
J(mu,nu)=int_0^infinity t^(mu-1)K_nu(t)^2 dt
 =2^(mu-3)Gamma(mu/2)^2 Gamma(mu/2+nu)Gamma(mu/2-nu)/Gamma(mu).
Use mu=13/2, nu=0,1, within that range. Thus the EXACT untruncated
positive geodesic multiplier has operator bound

 ||A_M(P)f||2<=c_* M^(-5/2)|| |P|^(9/2) f||2,
 c_*=c_g[J(13/2,1)-J(13/2,0)]/(2pi^2).

For q_ang=2000^4, monotone sum-plus-integral gives the COMPLETE angular sum

 S_5/2=sum d_n M_n^(-5/2)
 <=C^(5/4)[2(q_ang+2)/(q_ang+1)^(5/4)+8/(q_ang+1)^(1/4)].

Consequently the whole retained nonlocal geodesic column satisfies

 ||Qgeo,nl[x_hat,.]||Et* <= Kgeo N_9/2(Vx_hat),
 Kgeo=4pi c_* S_5/2 CLt/(B^2 beta).                    (8)

Dominated spectral integration proves convergence of the complete angular
sum and flow integral on Dom(|P|^(9/2)), and extension against ALL Et tests
uses the established ||V eta||2<=CLt||eta||Et. The exact ultraviolet
logarithm is bounded, not removed. All reference flow returns remain.
The common a0/a2, LLL and Weyl terms remain in the separate LOCAL form;
none is erased by (8). This is a theorem for the reference nonlocal part,
not a new bound on the actual-minus-reference operator or later applications.

## Regularity of the actual trial and the first unresolved comparison

The model boundary spectral measure can be derived without sampled modes.
The regularized Green diagonal is -digamma(k+1), mu=a k(k+1); its
resolvent differences are the accepted Q Gram identity. Stieltjes inversion,
with the digamma reflection identity [DLMF 5.5.4](https://dlmf.nist.gov/5.5.E4),
gives density (1/2)tanh(pi sqrt(p/a-1/4)) on p>=a/4, no lower atoms.
The hypotheses follow from the regular-flux/decaying Friedrichs realization:
apply its spectral theorem first to the L2 vector R_mu delta_b and recover
the boundary measure by multiplication by (p+mu)^2. Analytic continuation
of the resolvent difference has no poles below a/4; its cut yields the
stated density. In particular it is bounded by 1/2 on the whole spectrum.

For m=min(1,muT,muS,lambda_l)>0, (5) therefore proves

 ||(1+H)^(9/4)(v0/At)||L2(dx) <= K0 Cv/m^(5/2)=Nmodel. (9)

Indeed integrate (1+p)^(9/2)/(p+m)^6 with density <=1/2 and use
1+p<=(p+m)/m. Its upper integral is m^-5. Also
||xi0||L2(ds)<=At sqrt(gamma) Nmodel, since
|(2/3)(p+2a)/(p+muY)|<=1. The frozen coefficient remainder in (4) has
an original Et* density debit

 Emodel=Cxi At sqrt(gamma) Nmodel[|r1|/(8a)+|r0|/(4a^2)]. (10)

This is only the explicit model algebraic remainder in the full residual;
it is not the unknown finite-domain/return difference.

The positive K0 in (5) gives v0=-At K0 x^2 log(x)/32 plus smoother or
analytic terms. Formula (6) preserves that coefficient: the extra York
resolvent r starts its nonanalytic term at x^3 log(x), and its derivative
is multiplied by a coefficient vanishing at the bolt. Since proper bolt
radius rho=sqrt(2x)+O(x^(3/2)), the ACTUAL Vx_hat contains

    -(At K0/64) rho^4 log(rho).                       (11)

The coefficient is nonzero. Its fifth radial derivative is asymptotic to
-3 At K0/(8rho), so its squared fifth-derivative norm diverges logarithmically
in the regular two-dimensional disk measure. An H5 or C5-source quantum
theorem cannot be used on this trial. This is an AUXILIARY first-solve
field; no original geometric admission requirement is weakened.

The actual trial does belong to the paired energy form domain and to
Dom(|P|^(9/2)). Locally, rho^4 log(rho) lies in H^(4+s) for every s<1:
partition into dyadic annuli, whose H4/H5 interpolation costs are bounded
by a constant times (1+j)2^(-j(1-s)); their sum converges. A fixed smooth
bolt chart converts this to the angularly averaged flow norm. On the
outer seam the even extension starts at order four; its fourth derivative
matches and its fifth weak derivative is locally L2. The C2,1 metric's
possible higher distributional coefficients multiply vanishing lower
source jets. The remaining compact regions are smooth. This establishes
membership and the applicability of (8), NOT its numerical source norm.

The first unresolved QUANTITATIVE implication is a bound on the actual
N_9/2(Vx_hat), or a proved constant Cflow in

 N_9/2(Vx_hat)<=Cflow At sqrt(gamma) Nmodel.            (12)

The actual global York reconstruction (6), finite Q, rough paired seam
and regular-bolt return must all enter it. The quadratic spectral bound
(9) alone does not prove (12). A smooth surrogate or qualitative smoothing
would be insufficient. Reject an H5-based route because of (11), and
test the fractional route using its explicit allowance in the full budget.

## Complete budget and decision contract

The exact residual decomposition yields the sufficient first-stage test

 Esource+Emodel+Kgeo N_9/2(Vx_hat)+Ereturn+Enum <9.955e32. (13)

Ereturn includes the ENTIRE original local action/contacts and both York
return differences from the density used in (3)-(5). Enum includes any
later numerical application beyond the exact analytic profile definition.
Both remain UNKNOWN. The actual-minus-reference quantum applications and
later filters also remain unexecuted and are not included in this first
reference stage. Compute the attained source debit and model remainder,
Kgeo, Nmodel and nonzero bolt coefficient. Report the ceiling for Cflow
only with Ereturn=Enum=0, explicitly CONDITIONAL and not an available
physical margin. Report the stronger-H5 failure as an estimator failure,
not a construction obstruction. Do not call (13) certified if either
unknown remains. No further degree, tail or component optimization is
selected without an explicit contribution to this complete inequality.

Remaining to-do list: actual complete first reference residual and later
quantum applications; G_t and every original coupled/admission budget;
nonlinear existence, compatibility, matching, validity/noise and stability.
