# Complete fixed-input quantum comparison: kernel return test

2026-10-29 HST, E003/S169. Ordinary argument and numerical contract frozen
before execution. The preceding goal turn was NO PROGRESS: it restated the
plan. Preserve f384d0f, P32, the closed audit, original norms and the full
stable semiclassical wormhole objective. This tests a COMPLETE upper
estimator for the actual-minus-reference action on the accepted analytic
trial. A zero-centered norm enclosure is not a computed nonzero signed
center, and an excessive enclosure is not a construction obstruction.

## Operator, input, output and the domain implication

The input is the SAME finite-Q source-adapted field x_hat from
current-source-adapted-fractional-return.md, including both York equations,
the global Volterra integral, original paired mass/log phase, outer z^4
behavior and regular bolt. Write u=V x_hat. The output is the functional

    D[x_hat,eta]=(ell_actual-ell_0)[x_hat,eta]/N,
    ||D[x_hat,.]||Et* = sup_(||eta||Et<=1) |D[x_hat,eta]|.

Here Et^2=ell_t/N and N=N_halves*pi*B^2*beta/gamma. All orbital densities
below use one physical rank-two half, so time integration and division by
N give the SAME factor 2/(B^2 beta). No radius factor or physical output
comparison is removed. The final C5 geometric admission map is different
and remains unproved.

The trial has a nonzero rho^4 log(rho) term. Thus the spatially analytic
source theorem in current-localized-quantum-source-test.md does NOT apply.
There is nevertheless no missing definition of D on this fixed trial:
sections 2-4 of current-geodesic-two-sided-remainder.md give a symmetric,
commonly renormalized bounded form on the actual Et closure, including the
principal logarithm and all contacts; the accepted trial belongs to Et.
This is an energy-dual extension, not a claim that all energy tests are
smooth metric paths or that metric variations are controlled on a nonlinear
neighborhood.

One can also check the proposed fixed-operator starting point directly.
On the actual paired first-order domain, Dsharp is anti-self-adjoint and
||(Dsharp+M)^-1||<=1/M. The present u is bounded W2,infinity, has the
original vanishing seam data and is regular at the bolt. Multiplication
by exp(tu/2) preserves the first-order matching and the composed-square
matching. The conformal covariance identity then uses Dsharp+M exp(tu)
on a fixed domain. If exp(|t| ||u||infinity)-1<1/4, its resolvent Neumann
series converges uniformly in M. This finite-mass statement alone does
NOT justify the common angular finite part, a differentiated trace, or
interchanging its regulator limit. The existing bounded-form argument
supplies the definition here; a sharper matched trace estimate is still
needed for a useful evaluation. No smooth surrogate or selected spectrum
is substituted.

## Replace the global point-evaluation losses by an explicit kernel

Use precisely the actual-to-cap identities and buffers of sections 3-5
of stationary-curved-quantum-tt.md. The cap potential v satisfies
|v|<=.026 and its covariant kernel is

    B=e^(-sigma(x)/2)(Dflat+M+Mv)^-1 e^(-sigma(y)/2),
    .974<=exp(sigma)<=1.026.

For the real flat kernel use k(r)=(sqrt(2)/pi) exp(-Mr/sqrt(2))/r.
After weight exp(Mr/4), its L1 norm is L/M, where
L=2sqrt(2)/(1/sqrt(2)-1/4). The elementary triangle inequality gives

    sup_r r [k_weighted^{*n}](r) <= n*(sqrt(2)/pi)*(L/M)^(n-1).

Summing ALL Born insertions therefore gives kernel constant
Kreal=(sqrt(2)/pi)/(1-.026L)^2. This is an absolute convergent kernel
series, not a truncated Born expansion. Covariant endpoint densities and
the physical measure contribute at most 1.026/.974 to its squared row.
If the commutator slope is s and its proper separation from the endpoint
is d, Euclidean separation is at least d/1.026. Consequently

    ||ev_x B[D,chi]||2^2 <= 8 s^2 E1(Md/3).             (1)

The numerical check verifies .026L<1 and
2pi*(1.026/.974)*Kreal^2<8. The weaker exponent in (1) keeps the full
distance/density conversion. The estimate is uniform in both endpoints
of the retained local pair. Reverse rows use the same adjoint/grading.
The cap is a comparison component; A=(Dactual+M)^-1 still contains the
complete paired boundary conditions, thermal sum and regular-bolt returns.

For the local pair take s=260,d=1/132, so c_local=1/396. The exact identity
E=A-B=-B[D,chi]A[D,chi]B, together with ||A||,||B||<=1/M, gives

    |E(x,y)|<=R_M^2/M, ||E(x,.)||2<=260 R_M/M^2,
    local bubble Schur <=40R_M^2+2*260^2 R_M^2/M^2.

The factors include physical rank two and both cross terms, using the
accepted cap L1 norm 10/M. For the far pair use s=512,d=1/263;
the one-sided identity gives bubble Schur <=2R_far,M^2. This retains
the original partition into local and far pairs, without adding a wall.

The COMPLETE angular counting function satisfies N(M)<=4CM^2 and the
gap is M>18. Integration by parts, not a mode truncation, gives

    sum d E1(cM) <=4C/c^2,
    sum d E1(cM)/M^2 <=4C[Lc^2+(1+2/e)Lc+3/e],
    Lc=log(1/(18c)), 18c<1.                              (2)

For the second inequality split the positive integral at cM=1 and use
E1(t)<=log(1/t)+1/e below one and <=exp(-t)/t above it. The source and
test are never restricted to finite quantum modes. Also include the
entire long-flow reference kernel: its row is at most
2pi*18C*256^4. Doubling the sum of local, far and reference rows gives
the row-sum-subtracted global form coefficient G0 C. Its mixed bound is

    (2/(B^2 beta))*G0*C*||u||2*||Veta||2.                 (3)

This is the whole global gluing entry, not its raw bubble without the
diagonal row correction.

For the constant contact use the holomorphic identity on |z-M|<=M/4.
The flat kernel is bounded by (5/3) times the real kernel at mass .75M.
Weight exp(Mr/8). Set

    Lcpx=2sqrt(2)*(5/3)/(.75/sqrt(2)-1/8),
    Kcpx=(5sqrt(2)/(3pi))/(1-1.25*.026*Lcpx)^2.

The guards 1.25*.026*Lcpx<1 and
2pi*(1.026/.974)*Kcpx^2<32 give

    R_cpx,M^2<=32*260^2 E1(M/792).

Since ||A(z)||<=4/(3M), rank two gives
|tr(A-B)(x,x)|<=(8/3)R_cpx,M^2/M. Cauchy applied to this HOLomorphic
identity bounds -M^2 F'(M)-MF(M) by (40/3)R_cpx,M^2. Equation (2)
then gives a complete constant-contact gluing bound Pglue*C. This
improves only its gluing contribution. The actual cap tadpole, constant
bubble, finite common restoration and seam term are NOT removed.

## Entire remaining estimator, using established source comparisons

Reuse the scalar derivative constants S0,S1,S2 and connection D0 from
f384d0f, with Z=At sqrt(gamma) Nmodel. Since the lower-order input norms
are bounded by the accepted order-9/2 input norm, the proved reconstruction
formulas imply

    U0=S0 Z >= ||u||2,
    U1=sqrt(10) S1 Z >= ||grad_actual u||2,
    U2=10(S2+D0*S1)Z >= ||Hess_actual u||2.               (4)

These intentionally coarse but FINITE bounds retain both equations and
the global return; qualitative smoothing is not used as a constant.

The full test ball has ||Veta||L2(z^-4 dz)<=CV<2.916127 from 1b04eae.
The retained energy lower bound and original gradient return give

    ||grad Veta||2 <= CG,
    CG^2=1e5 beta/[(1-epsilon)*.99*(1-1.179e-8)],
    ||Delta Veta||2<=CD,

where CD is the accepted positive-time energy comparison. Thus the
averaged flow second derivative is <=sqrt(3CD^2/8+CG^2/4). Use U1 and
U2 also as conservative first/second flow bounds for the actual input.
The finite second contact is bounded separately by 1e13 C U1 CG.

For the complete singular gain-two remainder retain K=1e13,R0=1/400,
eps_split=1e-12, b=sqrt(eps_split/(4K)), alpha=eps_split/3 and
beta_gain=K[1/4+2log(R0/b)]. The accepted finite-difference proof, followed
by symmetric polarization in its positive comparison norm, gives

    local <=(2/(B^2 beta))*2pi*C
      *sqrt(alpha U2^2+beta_gain U1^2)
      *sqrt(alpha(3CD^2/8+CG^2/4)+beta_gain CG^2).          (5)

The 2pi uses the original unaveraged angular measure. Keeping the second
contact separately is conservative; no previously charged contact is
assigned zero. eps_split is a Taylor splitting choice, NOT a regulator,
source normalization, polynomial degree, or changed admission norm.

After replacing only its gluing part, section 5 of the curved proof
retains the global cap constant-contact bound
|Theta_cap|<=1e18 q^2(1+1/Y). Its cold bound fits this reserve. The
dimensionless orbit density is Theta/C. At the actual seed Q<=beta+u_mass
<1/5 (the cubic/quartic beta polynomial decreases on [0,1]), hence
Y>=sqrt(5)z. With the weighted test comparison this gives the WHOLE cap
contact contribution

    cap0 <=(2/(B^2 beta))*(1e18 q^2/C)
                          *(1+1/sqrt(5))*U0*CV.          (6)

The seam distribution vanishes on the original u(0)=u'(0)=0 phase, not
on a substituted domain. Its finite boundary representative is retained.
The contact gluing debit is (2/(B^2 beta))*Pglue*C*U0*CV.

Retain the explicit known curvature polynomial with |F_eff|<100:

    curvature <=(2/(B^2 beta))*100C
                *[U1*CG/(36pi)+816U0*CV/(960pi)].        (7)

Keep the shared vacuum uncertainty using its accepted relative E0 bound
<1e-21, divided by (1-epsilon)*.99*(1-1.179e-8), times the actual Et trial
upper. Add epsilon/(1-epsilon) times that upper for the retained energy
comparison. This deliberately overcounts known small local differences;
it does not choose the quantum stress to cancel a residual.

Add (3),(5),(6),(7), the finite second contact, contact gluing, vacuum and
energy comparison. This is a complete absolute comparison for D on the
explicit trial, with the exact UV logarithm canceled against the SAME
reference before norms. The common regulator and finite restoration are
those of the accepted continuum argument. It is not a new regulator
theorem for independent phases or neighboring metrics.

## Frozen decision test and error scope

Compute every debit with outward Arb arithmetic and compare their sum
plus the accepted first-reference residual with the sufficient actual
first-residual screen 1e33. Also report each debit against the more
conservative remaining residual allocation 4.308597e32. That allocation
is CONDITIONAL, not a free physical margin. The previously recorded
reference result is not changed.

Report the old and new global-return and contact-gluing coefficients,
the complete estimator's overrun, and the subtotal after removing the
largest term hypothetically. If that subtotal still overruns, optimizing
only the largest term is not a selected next step. The halved-account
thresholds for cap0 and the remaining subtotal are sensitivity numbers,
NOT established constants or sufficient budgets for the full wormhole.

All parameter and bound-arithmetic errors are enclosed for this exact
analytic trial. A nonzero signed center and numerical quantum application
have NOT been computed. No implementation errors are assigned zero.
P32's original defect, mass column, graph and physical outputs, nonlinear
admission/contraction, exact compatibility, physical matching, validity/
noise and all-sector stability remain separate unresolved requirements.

The specified alternative for a failed estimator is the FULL matched
unvaried heat/resolvent action on this same input, treating cap contact,
gain-two kernel and their subtraction together. Its first unsupported
implication is a quantitative energy-dual short-time remainder after
the full common angular sum. Finite-mass analyticity and small unvaried
heat defects do not by themselves prove that differentiated trace bound.
The present complete budget will determine whether this method change,
rather than further transport or isolated contact work, is justified.

Remaining to-do list: execute this full estimator and compare methods;
complete signed response and applications; original coupled/admission
budgets; nonlinear existence, compatibility, matching, validity/noise,
all-sector stability and the full stable semiclassical wormhole objective.
