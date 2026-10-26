# Angular-resummed actual quantum operator: strong coupled solve still unresolved

2026-10-25 HST, E003/S169. Ordinary argument before frozen execution.
The full stable-wormhole objective stays active; the feasibility audit is closed.

This is an attempt to assemble the missing signed operator, not a completed
coupled feasibility milestone. It gives an exact all-angular representation
of the massive canonical action on the ACTUAL paired orbit operator, and an
interval implementation of its angular resolvent/bubble weights. The complete
thermal/spatial operator application and original strong quotient solve are
NOT executed. No kappa_m, physical correction or nonlinear gate is certified.

## 1. Actual operator and the first unsupported assembly step

Let L_u be the nonnegative geometric orbit Dirac square with the original
complementary composed domains, paired mass/log phase and regular disk bolt.
Use its physical-half heat trace K_u(t), the original thermal period and
all thermal labels. Keep the fixed-Hilbert-space metric vertices L_i,L_ij
of current-canonical-heat-restoration-assembly.md, including the second
vertex. They are interpreted through the established domain/form argument;
the square's operator domain is not frozen by fiat.

Set q=2000^4, M_n^2=n(n+q)/C, d_n=2(q+2n), n>=1,
Theta(t)=sum d_n exp(-t M_n^2), ell=log(mu0 B),
mu_orbit=exp(ell)/sqrt(C). The exact massive canonical zeta function is

 Z(s,u)=mu_orbit^(2s)/Gamma(s) int_0^infinity
                                  t^(s-1) Theta(t) K_u(t) dt,
 Gamma_mass(u)=Z'(0,u)/2.                               (1)

Equation (1) is the SAME common regulator, not separate finite parts of
two divergent sums. The orbit spectrum below always includes all thermal
frequencies, the paired boundary and finite-domain returns. No product
geometry, smooth even double or finite mode selection replaces L_u.

The first unsupported passage in the proposed strong solve is from the
admitted-family WEAK metric Hessian to a bounded operator into the ORIGINAL
two-row strong source space, including transition equations and border.
The weak C2 theorem does not supply this map or its norm. A usable
all-angular operator application is needed to attempt it; summing
modewise positive bounds does not produce it. The angular resolvent
can be summed exactly, but inserting its finite part into an unregulated
orbit trace would still be unsupported. Sections 2-3 separate those issues.

## 2. Complete angular resolvent and two-variable bubble

For a scalar spectral argument lambda>=0, define

 R_s(lambda)=mu_orbit^(2s) sum_n d_n (lambda+M_n^2)^(-1-s).

Its pole at s=0 has residue 2C, independent of lambda. Write
a+b=q, ab=C lambda, choosing either ordering of the two roots. Then

 R(lambda):=FP_(s=0) R_s(lambda)
          =2C[2ell-psi(1+a)-psi(1+b)].                    (2)

Proof: at lambda=0 the retained angular zeta identity gives
R(0)=C[4 EulerGamma-2H_q+4ell]. Subtract that value BEFORE summing.
The difference is absolutely convergent, and partial fractions give
2C sum[1/(n+a)+1/(n+b)-1/n-1/(n+q)]. The convergent digamma identity
[DLMF 5.7.6](https://dlmf.nist.gov/5.7.E6) proves (2). Differentiation
is uniformly absolutely convergent on compact nonnegative lambda intervals.
For lambda beyond q^2/(4C), the roots are conjugate with real part q/2;
the symmetric sum is real. No pole crosses the physical spectrum.

The complete angular bubble weight is

 B(lambda,mu)=sum_n d_n/[(lambda+M_n^2)(mu+M_n^2)]
            =[R(lambda)-R(mu)]/(mu-lambda),
 B(lambda,lambda)=-R'(lambda).                          (3)

It is positive, symmetric and decreasing in each nonnegative argument.
At coincident arguments and at the root coalescence the values are
continuous; neither apparent divided difference is a physical singularity.
The same convergent series proves the WHOLE-spectrum bound

 0<B(lambda,mu)<=B(0,0)
   = (2C^2/q)[pi^2/6-psi'(q+1)].                       (4)

The telescoping identity behind (4) is
(q+2n)/[n^2(n+q)^2]=(1/q)[1/n^2-1/(n+q)^2].
The trigamma series is [DLMF 5.15.1](https://dlmf.nist.gov/5.15.E1).
This is an exact complete angular sum, not a finite-N estimate.

The implementation retains interval parameters. Away from coalescence it
uses the stable small root a=2C lambda/[q+sqrt(q^2-4C lambda)] and b=q-a.
Near coalescence, put r=1+q/2, d=q^2-4C lambda. The convergent formulas are

 R=4C[ell-psi(r)+sum_(k>=1)(d/4)^k zeta(2k+1,r)],
 B(lambda,lambda)=4C^2 sum_(k>=1)k(d/4)^(k-1)zeta(2k+1,r). (5)

For t=|d|/(4r^2)<1 and truncation after k=N, the omitted absolute
errors in (5) are at most

 4C t^(N+1)/(1-t)[1/r+1/(2(N+1))],
 (4C^2/r^2)t^N{[(N+1)-Nt]/[r(1-t)^2]+1/[2(1-t)]}.       (6)

These follow from the decreasing-integral bound
zeta(2k+1,r)<=r^(-2k-1)+r^(-2k)/(2k). They bound the entire remainder,
including interval balls crossing d=0. No spectral branch is rounded away.

## 3. Common regulator preserved: do not stop at the scalar finite part

Equations (2)-(4) alone do NOT authorize an ordinary orbit trace of
R(L)L_ij or the bubble with differential metric vertices. R(lambda)
grows logarithmically, and the vertices have order two. The scalar bound
(4) is not a strong quantum-response or trace bound. A separate angular
finite part followed by an arbitrary orbit cutoff can lose finite terms.

Here is an exact common-regulator regrouping that avoids that step. The
complete angular heat expansion, with the massive LLL-excluded convention,
is

 Theta(t)=2C/t+theta0+theta1 t+O_(q,C)(t^2),
 theta0=-q-1/3, theta1=(q^2/6-1/30)/C.                  (7)

For completeness apply [Euler--Maclaurin, DLMF 2.10.1](https://dlmf.nist.gov/2.10.E1) to
f(x)=2(q+2x)exp[-t x(x+q)/C], x>=0. Its integral is 2C/t;
the endpoint is -q; -f'(0)/12 contributes -1/3+q^2 t/(6C),
and f'''(0)/720 contributes -t/(30C)+O(t^2).
The next endpoint term is O(t^2). The remainder with six derivatives
is O(t^2): after scaling x by sqrt(t), its derivative L1 norm is
O(t^2) for fixed finite q,C. Gaussian decay controls infinity. This
establishes (7); it does not assign a small numerical remainder constant.

Use the fixed auxiliary shift 1, independent of u, and put

 A_u=L_u+1, b0=2C, b1=2C+theta0, b2=C+theta0+theta1,
 W_theta(t)=exp(t)Theta(t)-b0/t-b1-b2 t.

Since M_1^2>324, this shift introduces no infrared issue. For initially
Re s>2, absolute convergence and ordinary Fubini give EXACTLY

 Z(s,u)=mu_orbit^(2s){
   b0 zeta_Au(s-1)/(s-1)+b1 zeta_Au(s)+b2 s zeta_Au(s+1)
   +[1/Gamma(s)] int_0^infinity t^(s-1) W_theta(t)Tr exp(-t A_u)dt}. (8)

The original admitted-family paired heat argument gives K_u(t)=O(t^-1)
and the same bound for its first two metric variations. W_theta=O(t^2)
at zero, so the last integral and those variations converge there at
s=0. At infinity the shift gives exponential decay, including the
subtracted polynomials. Thus the remainder is holomorphic near s=0.
The established weak heat expansion supplies the meromorphic continuation
of the three zeta_A terms: K=a0/t+a1+a2 t+O(t^(3/2)), with two metric
variations, suffices. Multiplication by exp(-t) preserves that expansion.
The zeros of 1/Gamma(s) cancel its ordinary heat poles at s=0,-1;
zeta_A is regular there and has only a simple pole at s=1. No t log t
coefficient is silently removed. The actual seam's later fractional/log
terms remain in the remainder. No smooth-seam heat theorem is substituted.

With ell_A=ell-(log C)/2, (8) yields the complete massive action

 2 Gamma_mass =
 -b0[zeta_A'(-1)+(1+2ell_A)zeta_A(-1)]
 +b1[zeta_A'(0)+2ell_A zeta_A(0)]
 +b2[FP zeta_A(1)+2ell_A Res zeta_A(1)]
 +int_0^infinity W_theta(t)Tr exp(-t A)dt/t.             (9)

In particular the pole/residue term in the third line is retained.
Both metric derivatives of every term in (9), including the tadpole
and boundary contributions, are required. For compact smooth tangents
this follows by the common positive-regulator identity; the established
uniform admitted-family weak bounds pass it to the actual mass/log phase.
This gives a C2 identity on those ADMITTED families, not a strong
operator estimate or an embedding of an unknown nonlinear X ball.

Specifically, the hypotheses used here are the original fixed first-order
paired domain, positive ellipticity, bounded W1,infinity coefficients,
the cubic plus fixed-finite-log weighted C5 phase, and the bounded C5
interior/regular-disk atlas. The proper-coordinate map is C2 into
W3,infinity, not spatially analytic. The weak collar proof includes the
distributional completed coefficient and the entire thermal image sum;
the finite-atlas proof retains its differentiated gluing defect. Those
arguments provide a qualitative family-dependent constant. They do not
bound it in the requested strong norm. No third geometric neighborhood
variation or physical-X embedding is imported from them.

An explicit positive-time recipe retains both Hessian vertices. On the
original common closed form, write A_i=D_i A and A_ij=D_ij A. Then

 D_ij Tr exp(-t A)
   =t^2 int_0^1 Tr[exp(-(1-a)t A) A_i exp(-a t A) A_j] da
                         -t Tr[exp(-t A) A_ij].          (10)

Products in (10) mean their closed-form Duhamel limits, including the
actual composed boundary realization, not a newly imposed common square
domain. At positive regulator the identity follows by differentiating
the heat evolution and cyclic trace; the stated uniform weak bounds
justify its integrated C2 limit. Apply (10) to the last term of (9) and
differentiate ALL three zeta terms in the SAME identity. The second term
is the coincidence vertex; dropping it is not a permissible bubble
estimate. The finite disk heat kernel contains the entire thermal sum
and finite-domain returns. No numerical value for (10) is assigned here.

Equation (9) replaces the ENTIRE previous canonical massive allocation.
Do not add its old a0/a2 restorations a second time. Einstein, Maxwell,
the specified R2 action, four-dimensional Weyl conversion and the separate
LLL action are still included once in the full coupled operator. The
physical factor W=-N pi B^2/(G gamma), its sign and beta in energy/dual
norms remain unchanged. The auxiliary shift changes no physical mass,
state, boundary condition, temperature or renormalization scale.

## 4. Executed scope and remaining quantitative obstacle

The ordinary identities above evaluate the COMPLETE angular resolvent
and bubble for every nonnegative orbit spectral argument, and preserve
the common regulator for the full action. The supporting implementation
encloses their scalar values, coalescent limits and (4), using the actual
parameters and fixed ell. The external spectral arguments in its checks
are function evaluations, not selected physical thermal or angular modes.

NO full operator norm follows by inserting (4) between the differential
metric vertices: those vertices are unbounded, and the separate orbit
traces are not trace class. The first remaining unsupported numerical
step is evaluation and strong control of the combined metric derivatives
of (9) on the full near/transition quotient, with its original border.
This requires the actual orbit heat/resolvent kernels or conforming spectral
functional calculus, not an ordinary trace of (2) alone. Neither the
uncomputed zeta_A terms nor the last integral is assigned zero. Their
separate coefficients are NOT physical response bounds.

The former exact transport error remains available in its proved energy
and outer scopes. No further refinement is made. This new representation
does not shrink that error by renormalization, and no error is charged twice.
The previous energy-only reconstruction failure remains: the annular
equations cannot be replaced by the outer graph and three phase values.

The reference mass lift is J_m=int_0^z t^3/Q^2,
n_m=-2 beta Q J_m/z, b_m=beta z^3/Q+beta Q'J_m. Its uncut areal
variation f=2xz, alpha=0 gives (r f)'=alpha'=0, so both reference
equations cancel. The ORIGINAL chi_m is one through .1 and zero by .2:
f_cut=2xz chi_m-2F n_m Dchi_m, alpha_cut=n_m Dchi_m. These transition
terms stay in both equations. At the regular bolt the cut lift is zero,
so its reference bolt shift is zero; integrating the full transition
gives reference mass -xB per unit coefficient. It is not a global
fixed-temperature zero mode or an adjustable physical quantum source.

The ORIGINAL conversion g_old=g_lifted+m d_chi has
2133<||d_chi||H3<2174. Together with the known LINEAR
|m|<1.627532e-15, its debit <=3.538254568e-12 is preserved, leaving
near-shape allowance >6.4617454e-12. No derivative constant of the
separate mass cutoff is invented. The undifferentiated scalar quartic
and B0(z^4-z^(4+e_*))+Delta z^(4+e_*), with the exact unrounded e_*,
remain coupled to the fast remainder. No actual post-elimination kappa_m
has been computed.

At CONDITIONAL kappa_m<=5e-5, the residual requirements still are

 Y_other+4.084e-22 A_o+Y_phase_return<=5.34235e-12,
 eta_other+.0195239 gamma_E+3.712e-12 A_o J_E
                                      +eta_phase_return<.0012069.

These unknown entries have not been shown to fit. The conditional fast
screen .135973<.1599164 keeps its two-fast/admitted-family scope; it is
not control of independent mass/quartic/slow phases or original physical X.
The angular identities give zero SYMBOLIC angular truncation error; their
interval evaluations enclose only scalar function roundoff. They assign
NO error or center to the unevaluated orbit application, so neither budget
is reduced. This is the quantitative limitation of the supporting work.
The original nonlinear requirements remain, with unknown comparisons,

 E_center+1e-8 alpha_embed <= 1e-11,
 R_X alpha_embed*.028913+Z0_other <= .25,
 R_X alpha_embed^2*1.0706e10+Z2_other <= 5e7,
 Y0+1e-8 Z0+5e-17 Z2 <= 1e-8, Z0+1e-8 Z2 < 1.

All nonlinear admission/contraction, compatibility, matching,
semiclassical validity/noise and all-sector stability gates remain open.

RESEARCH DECISION: retain the construction. Use the exact angular-resummed
GEOMETRIC operator, with (8)-(9)'s common regulator, to assemble the signed
quantum part of the full strong solve. This removes the need to approximate
the angular tower by a finite physical selection; it does not certify the
coupled return. Do not use the digamma finite part as an unregulated orbit
trace, or optimize the already controlled cold/scalar/transport components.

ONE NEXT CALCULATION, PROPOSED NOT EXECUTED: apply the combined metric
Hessian of (9), alongside both local original equations, to a conforming
full quotient parametrix and certify its strong residual and original
row maps. Retain every thermal contribution, interface/phase equation
and finite-domain return, then test all coupled columns and the unchanged
mass/residual budgets together. The complete signed solve is still missing;
this turn establishes its angular operator representation, not a solution.

Remaining to-do list: complete signed strong quotient/border operator
application and residual; original physical-X/neighborhood comparisons;
nonlinear admission/contraction; exact compatibility; physical matching;
validity/noise and all-sector stability; charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
