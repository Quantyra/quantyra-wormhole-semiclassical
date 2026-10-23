# Common canonical restoration and the remaining geometric heat Hessian

2026-10-23, E003/S169. Ordinary argument before frozen execution.
Previous goal turn: progress. The full stable-wormhole goal is incomplete;
the feasibility audit stays closed and its physical thresholds are fixed.

This assembles the remaining canonical local restoration and exact LLL
weak rows on the actual finite seed. Keep the sufficient known enclosure
of the vacuum coefficient; do not commission a separate spectral-constant
optimization. The remaining MASSIVE quantum term is specified as the
commonly subtracted GEOMETRIC heat Hessian. It is not replaced by the
mass-vertex action alone or declared evaluated by the local calculation.

## 1. Common regulator, finite scale and physical half

Use the same affine canonical shape u, gamma0, Q, B, C, paired domain,
temperature and W=-N pi B^2/(G gamma) as
[the local action assembly](current-finite-local-action-assembly.md).
The canonical orbit is (B^2/C) exp(2u)gamma0. In dimensionless gamma0
units let L_u be the actual paired geometric orbit Dirac square, and
K_u(t) its physical-half rank-two heat trace, including all thermal modes.
Let

    Theta(t)=sum_(n>=1) d_n exp[-t n(n+q)/C],
    zeta_h(s)=sum d_n[n(n+q)]^-s,
    Lscale=log(mu0 B), mu0=1/b0.

The factor N for identical exteriors multiplies the action; it does not
change the queried common radius. On the smooth core the common Mellin
regulator is applied to Theta(t)K_u(t), BEFORE taking its finite part.
The positive massive squares avoid a zero-mode infrared issue here.
The canonical LLL is the separately declared, kernel-free massless
Polyakov allocation on the same paired HH domain.

For the existing heat representative, write

    K_u(t)=[A0(u)/t+A1(u)+t A2(u)]/(4pi)+R_u(t),
    A0(u)=2 integral exp(2u),
    A1(u)=-(1/6) integral [R-2Delta u], R=-Q'',
    A2(u)=-(1/120) integral exp(-2u)(R-2Delta u)^2.    (1)

The integrals include the canonical time measure. A2 in (1) is the
INTEGRATED paired completion: the bulk 2Delta R term and its seam
distribution are combined, as in the retained common action. A bare
half-disk bulk Laplacian must not be erased. The outer phase has
u=c3 z^3+O(z^4 log^p z), so u(0)=u'(0)=0; the geodesic boundary
representative and the divergence have no extra limiting Hessian on
these tangents. At the bolt use smooth disk coordinates, not an inner
wall. This is the existing one-sided/paired argument, not a smooth
reflection theorem for the cubic seam.

The local Mellin factors, apart from their heat coefficients, are

    C^-1 exp(2s Lscale) zeta_h(s-1)/(s-1),
    exp(2s Lscale) zeta_h(s),
    C exp(2s Lscale) s zeta_h(s+1).                  (2)

Use the fermionic effective action +zeta'(0)/2, including its scale.
Since zeta_h(-1)=1/30-q^2/6 and zeta_h(s+1)=2/s+Fang+O(s), (2) gives

    J0=zeta_h'(-1)+(1+2Lscale)zeta_h(-1),
    F_eff=Fang+4Lscale=4EulerGamma-2H_q+4Lscale,
    Gamma_a0= -J0/(4pi C) integral exp(2u),
    Gamma_a2= -C F_eff/(960pi)
                         integral exp(-2u)(R-2Delta u)^2. (3)

The a1 allocation is a fixed multiple of A1. Its first and second
variations vanish on the stated physical phase; no value of zeta_h'(0)
is needed to compute these rows. In particular the fourth finite action
coefficient has F_eff, without a spurious extra +2. The source's
mass-weighted Green formula has a different Mellin factor; its +2
cancels the local anomaly and is not an additional action coefficient.

These are the same allocations already retained in the action. They
are being assembled, not added as new finite couplings or counted twice.

### Fixed-scale correction discovered during execution

The first scientific run rejected the older claim |log(mu0 B)|<1e-21.
At the fixed physical scale mu0=1/b0 and actual recentered B, the correct
enclosure is Lscale=2.9452431301493e-9 with uncertainty <3.29e-23; in
particular 0<Lscale<3e-9. The c_R-dependent radius shift must be retained.
This agrees with the already executed bolt-contact and vacuum-budget
calculations; no renormalization scale is changed here.

This supersedes that narrow assertion in stationary-seed-quantum-form.md
section 6 and stationary-curved-quantum-tt.md sections 1 and 10. Their
particular fourth-order margins survive with the corrected hypothesis:
Fang<-1/6 and 4Lscale<12e-9 still imply -F_eff>3/25, so the flat
C p^4/(1000pi) lower margin holds. Direct actual-parameter evaluation
also verifies x=g B^2>4e9, (log x)/2-Lscale>9 and kappa4=-F_eff/4>13.9,
as required for the stated outer and interior reserves; |F_eff|<100
preserves the local curved allowance. These are the precise uses of the
incorrect scale claim in those two proofs. Their other domain/continuum
arguments are not certified by these scalar checks. Preserve the failed
run and old sources; rerun the complete current assembly after correction.

## 2. Implement source and Hessian together

The known bound |zeta_h'(-1)|<604q^2 gives the sufficient interval

    c0=2G J0/(pi C B^2),
    |c0|<=U0=2[604q^2+(1+2|Lscale|)(q^2/6-1/30)]/(pi C B^2). (4)

The signed value is not evaluated, but (4) is a KNOWN upper enclosure,
not an unknown estimate. Preserve its shared scalar identity across
all entries: the vacuum matrix is c0 times one exactly integrated Gram
matrix. It is not a collection of independent arbitrary coefficients.
This does not change the old vacuum-contact feasibility decision.

Let vi, vj be canonical shape directions after the frozen York return,
delta_i=Delta vi, beta2=GC F_eff/(480pi B^2), and integrate against dz.
In the SAME I/W normalization, the local gradient and Hessian pieces are

    Ga0[i]=(c0/2) integral vi,
    Ha0[i,j]=c0 integral vi vj,
    Ga2[i]=beta2 integral[-2R^2 vi-4R delta_i],
    Ha2[i,j]=beta2 integral[4R^2 vi vj
                 +8R(vi delta_j+vj delta_i)+8delta_i delta_j]. (5)

The exact canonical LLL cocycle gives

    GLLL[i]=Gq/(12pi B^2) integral R vi,
    HLLL[i,j]=Gq/(6pi B^2) integral Q vi' vj'.        (6)

The sign in (6) includes the negative W: the physical LLL Hessian is
-q/(12pi) integral |grad v|^2. The kernel-free paired domain eliminates
an additional zero-mode normalization term. No massive angular sector
is discarded. The four-dimensional Weyl term was already included in
the previous local action and is NOT inserted again here.

Evaluate (5)-(6) for the actual three TT fields generated by J=z^4,
J=z^5 and xi=z^4,J=0. Each is smooth at the bolt and has the retained
outer mass/fast behavior. Integrate the full domain 0<=z<=z_h using
polynomial primitives. Store exact Gram coefficients and separate a2,
LLL and bounded a0 contributions. Add them to the previous local weak
matrix to form the KNOWN part of H. Its massive remainder stays open.
These local gradients are not the full seed source or Newton defect r1.

Check the a2 bilinear form independently against its strong row

    (GC F_eff/(60pi B^2))[
      Delta^2 v+R Delta v+Delta(Rv)+(R^2/2)v].       (7)

The endpoint terms in this comparison vanish for the current fields.
Record every integral and interval width. The retained local-assembly
width target is 1e-14 in these units; it is not a physical error tolerance.

## 3. The exact remaining term and its boundary hypotheses

At a finite common regulator, the massive determinant equals the three
Mellin terms and the corresponding regulated integral of R_u. Define
the time-reduced trace remainder r_u=gamma R_u/(2pi). Formally after
the justified regulator removal, its normalized weak Hessian is

    Hrem[i,j]=-(G/B^2) integral_0^infinity
                                  Theta(t) D_iD_j r_u(t)|_0 dt/t. (8)

Equation (8) identifies what must be evaluated. It is not a computation
of that integral. A valid numerical implementation needs a uniform
bound on the DIFFERENTIATED paired remainder. A bound on the undifferentiated
heat kernel or on a single stress component does not supply it.

For precision about the regulator-level vertices, put the geometric
Dirac expression on the fixed spinor Hilbert space:

    Dtilde_u=exp(-u/2) D exp(-u/2), L_u=-Dtilde_u^2,
    D_i=-{vi,D}/2,
    D_ij=(vi vj D+D vi vj+vi D vj+vj D vi)/4,
    L_i=-(D_i D+D D_i),
    L_ij=-(D_ij D+D D_ij+D_i D_j+D_j D_i).           (9)

Scalar multiplications preserve the actual paired domain. In particular,
the fixed-source outer traces do not introduce a new moving projector.
The heat Hessian before local subtraction is the Duhamel expression

    D_iD_j K(t)=t^2 integral_0^1 Tr[
       exp(-(1-a)tL)L_i exp(-atL)L_j] da
                                      -t Tr[exp(-tL)L_ij].       (10)

The last term is required. Equations (9)-(10) retain both the geometric
operator and its density through the fixed-Hilbert-space transformation.
The earlier fixed-metric mass action -Tr log(D+M exp u) is only one part
of a different representation. It also needs its two-dimensional Weyl
conversion and matching local subtraction. It must not be identified
with (10) without that argument, even if the finite-mode ODE checks pass.

If, for example, the reduced differentiated remainder obeys
|D_iD_j r_u(t)|<=Aij t^(1+alpha), alpha>0, at 0<t<=t0, the actual
counting bound Theta(t)<=4C/t yields the explicit omitted heat-cap bound

    |Hrem,cap|<=4GC Aij t0^alpha/(alpha B^2).        (11)

For a desired normalized entry error delta_H, the explicit condition is
Aij t0^alpha/alpha <= delta_H B^2/(4GC). Record the actual coefficient
B^2/(4GC) as well as the illustrative delta_H=1e-14 assembly threshold.
This is not the unweighted mass-action endpoint tolerance and does not
become the physical residual tolerance without the dual norm conversion.

The C2,1 seam and weight-four logarithms must be accounted for in Aij
and alpha. They are UNKNOWN here; the smooth closed-manifold expansion
does not assign them. For the singular gravitational adjoint, use a
controlled form-domain approximation and bound its limiting remainder;
do not differentiate a nonexistent smooth metric family at that field.
The prior full-tower numerical bolt-cut certificate remains usable for
the mass-action representation but does not supply (11).

There is a precise test for how accurate an OPERATOR remainder must be.
For an implemented H0 and its unknown remainder R=H-H0, let

    ||R||_(V to V*)<=eta_H,
    d0=H0 e_app+r1, db0=L_b-H0 psi_app,
    Q0=L_b e_app-<psi_app,d0>,
    E=||e_app||_V, P=||psi_app||_V,
    D0=||d0||_V*, Db0=||db0||_V*.

The actual two-residual certificate is bounded by

    |Q0|+eta_H E P+(Db0+eta_H P)(D0+eta_H E).       (12)

With M=1.667e-9-|Q0|-Db0 D0>0, a=EP and
b=EP+Db0 E+D0 P, it suffices (for a>0) that

    eta_H < 2M/[b+sqrt(b^2+4aM)].                  (13)

For a=0,b>0 use eta_H<M/b. This includes the error in the signed pairing
as well as both residuals; an entry-wise numerical tolerance cannot
replace eta_H. The actual E,P,D0,Db0,Q0 are still UNKNOWN, so (13)
cannot honestly be reported as an executed numerical physical test.
If (11) holds in operator norm with constant Aop, its threshold is
Aop t0^alpha/alpha < [B^2/(4GC)] times the right side of (13).
The existing trial-entry integrals do not establish such an Aop.

The Mellin convention is the one reviewed in
[Vassilevich, equations 2.23 and 2.33](https://arxiv.org/html/hep-th/0306138v3#S2.SS2).
The same source's section 9.1 explains the conformal allocation and its
zero-mode qualification; its smooth heat coefficients are not a theorem
for this nonsmooth paired seam. The physical LLL allocation and actual
boundary cancellation are the retained repository arguments cited above.

## 4. Decision and execution boundary

Continue the physical weak primal-dual method using the geometric,
fixed-domain Duhamel representation, with its second vertex retained.
Prioritize its differentiated seam/form-domain remainder before further
source optimization or extension of a raw mass-vertex estimate chain.
This execution closes the
specified common local/LLL assembly using the already sufficient vacuum
enclosure. It does not close the massive remainder, its full-domain
operator norm, metric admission or exact nonlinear compatibility.

The single next calculation is a quantitative differentiated heat-remainder
bound for the actual paired seam, extended in the form norm to the
logarithmic bolt adjoint. It must make (8) evaluable in the norm required
by (12)-(13); a smooth-core calculation alone does not qualify. Its success
would permit combining the complete rows on the physical pair to test
|Q_app|+||d_b||_V*||d||_V*<1.667e-9, exterior return <6.67e-10,
Y0<=5e-9 and every admission/nonlinear requirement. That full physical
feasibility test remains unexecuted; a local matrix is not a substitute.

Remaining to-do list: massive geometric heat/curved/boundary remainder
and physical forward/adjoint residual certificate; X/admission, nonlinear
closure and exact compatibility; matching, causal realization, semiclassical
validity/noise, all-sector stability, charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
