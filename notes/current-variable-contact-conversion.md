# Variable-direction common-regulator conversion at the actual paired seed

2026-10-24, E003/S169. Ordinary argument before acceptance execution.
Full stable-wormhole objective ACTIVE and INCOMPLETE; feasibility audit CLOSED.

The proposed actual strong remainder is not yet bounded. Inspecting its
inputs exposes an earlier representation gap: the useful radial mass
kernels are not the geometric Hessian without their common-regulator
conversion. The constant test established that a massless 2D anomaly alone
does not supply that conversion. Here we compute its COMPLETE first two
jets at the actual seed, for the admitted stationary mass/log directions.
This is a seed response identity, not a nonlinear strong-neighborhood theorem.

## 1. Closed operators and an exact cutoff identity

Use the original paired skew-adjoint D, complementary first-order domains,
grading, regular bolt, period, angular masses M_n and multiplicities d_n.
The domains are those in current-rough-paired-variation.md; the square is
always defined by composition. For a real admitted b, scalar multiplication
by exp(tb) preserves both first-order domains. Set

 A_r(b)=exp(-r b/2) [D+M exp(b)] exp(-r b/2), 0<=r<=1,
 P_r=A_r^* A_r, Pbar_r=A_r A_r^*.

These are closed operators on the actual realization. Accretivity of
D+M exp(b), also for its adjoint, gives an inverse and a gap uniform in
n after scaling by M. The endpoint P_1=-D_b^2+M^2 is the geometric square;
P_0 is the variable-mass square. The grading interchanges P_r and Pbar_r,
and commutes with every scalar insertion. No determinant-product identity
or independently frozen square domain is assumed.

At a positive common proper-time cutoff e,

 d_r A_r=-{b,A_r}/2,
 d_r P_r=-{b,P_r}/2-A_r^* b A_r,
 d_r Gamma_e=1/2 Tr[b(exp(-e P_r)+exp(-e Pbar_r))].        (1)

To obtain (1), differentiate +1/2 integral_e^infinity Tr exp(-tP_r)dt/t,
use cyclicity at positive time and A_r exp(-tP_r)=exp(-tPbar_r)A_r, then
integrate P exp(-tP) in t. The massive gap removes the upper endpoint.
The two traces in (1) agree by grading. Consequently, including the
ENTIRE angular tower and the physical-half, time-reduced trace,

 J_e(b)=Gamma_g,e(b)-Gamma_m,e(b)
       =sum_n d_n integral_0^1 Tr_bar[b exp(-e P_r(b))] dr. (2)

All thermal modes, original boundary data and bolt are in these traces.
For each e>0 the angular sums converge. The trace normalizations here are
the original ones: the final normalized action multiplies J by -2G/B^2.
Equation (2) is not an auxiliary physical geometry or a finite-mode model.

Let K(t)=exp(-tL), L=-D^2, and define

 T_e(u)=Tr_bar[u K(e)],
 I_e(u,v)=integral_0^e Tr_bar[u K(e-s)v K(s)] ds.

At b=0 the vertex is r L_v+2(1-r)M^2 v+(1-r)M[v,D].
The last term is grading-odd and has zero scalar-smeared heat trace.
Integrating r in (2), before taking any limit, gives exactly

 J'_e[u]=Theta(e) T_e(u),
 J''_e[u,v]=Theta(e)/2 [D_v T_e(u)+D_u T_e(v)]
                          +2 Theta'(e) I_e(u,v).         (3)

Here D_v T is the GEOMETRIC variation on the fixed spinor Hilbert space.
The second term in (3) is essential; dropping it loses finite contacts.
Symmetry I_e(u,v)=I_e(v,u) follows by cyclicity and s -> e-s.

## 2. Smeared coefficients on the original rough paired collar

Use the actual proper distance y, g=dy^2+a(y)^2 dtau^2, f=1/a, and
the exact reflected C2,1 coefficient from current-curved-paired-heat-bridge.md.
Both directions have their original cubic plus weight-four logarithmic
phase, bounded weighted C5 norms, W3,infinity even extensions and regular
disk-bolt fields. They need not be analytic or smoothly even at the seam.

Write Delta=div grad, R=-Q'', and use the time-reduced physical measure.
The smeared geometric coefficients are

 T_e(u)=(4pi)^-1 integral u[2/e-R/6
                         +e(-R^2/120-Delta R/60)]+o(e). (4)

The divergence is the paired distributional completion. Its seam delta
is retained before pairing, not erased using a half-disk bulk formula.
For the second factor in (3) the needed identity is

 I_e(u,v)=(4pi)^-1 { integral 2uv
   +e integral[-Ruv/6-grad u.grad v/3]
   +e^2 integral[(-R^2/120-Delta R/60)uv
                  +R grad u.grad v/60+(Delta u)(Delta v)/30] }+o(e^2). (5)

Here is a direct rough-collar derivation and remainder argument for (5).
The Brownian formula already used for the actual squared Dirac operator
acquires the insertion integral_0^1 v(y+sqrt(e) xi_alpha)dalpha:

 I_e(u,v)=1/(2pi) integral dy u(y) E[
   (integral_0^1 v(y+sqrt(e)xi_alpha)dalpha)
   A_h^(-1/2) exp(e B_h^2/(4A_h))], h=sqrt(e),            (6)

for the zero thermal image, with the complete image sum restored as in
the retained exact formula. A_h=int f(y+h xi)^2 and B_h=int f'(y+h xi).
The cyclic heat split is precisely insertion at bridge time alpha.
This does not replace the curved f by a constant.

Expand the explicit expectation to h^4. The seven required bridge moments
are (1/3,1/6,2/5,1/5,1/5,11/90,1/12), for mu2,mu1^2,mu4,mu3 mu1,
mu2^2,mu2 mu1^2,mu1^4. Wick's rule with covariance 2(min(s,t)-st),
integrated over every ordered time simplex, derives them. In particular
the coefficient of h^4 acting on v, before the y integral, is

 [-f^4 v''''+2f^3 f'v'''+3f^3 f''v''+2f^3 f'''v'
  +f^3 f''''v-5f^2(f')^2v''-10f^2 f'f''v'
  -7f^2 f'f'''v-4f^2(f'')^2v+9f(f')^3v'
  +23f(f')^2 f''v-12(f')^4v]/(-60 f^5).             (7)

Twice (7), with R=2f''/f-4(f'/f)^2, equals
f^-1[(-R^2/120-Delta R/60)v
             -(R Delta v+grad R.grad v)/60+Delta^2 v/30].
Weak integration by parts gives (5). This also fixes the R-gradient
coefficient, which cannot be inferred from the flat multiplier alone.

For the remainder, use the integrated six-derivative argument of the
retained heat bridge with the two additional u,v factors. Total spatial
derivative count remains six. At most one factor can have order >3;
move its excess derivatives by integration by parts onto the other
factors and cutoff. Every resulting factor has order <=3. Finite sums,
bounded W3 coefficients and the finite sixth bridge-supremum moment give
O(h^6)=O(e^3) in (6), uniformly on bounded admitted directions. Odd terms
vanish after bridge reflection. The same argument with a smeared u gives
the differentiated remainder o(e) in (4); the moving proper coordinate
and insertion are differentiated using the original weighted C5 phase.
The extra fourth derivative of u is multiplied by the coordinate shift
O(y^4), so it is bounded just as in section 5 of the retained bridge.

The construction is first justified for physical-side approximants that
keep the cubic exactly; weak integration by parts and translation
continuity pass to the original coefficient and directions. Complete
nonzero thermal images retain their exponential bounds. On regular C5
charts, including the bolt, the finite Gaussian construction with the
scalar insertion gives O(e^(5/2)) for the remainder in I_e. Separated
patch commutators retain their exponential bound. Both orders are o(e^2).
There is no claimed strong-operator or nonlinear-ball uniformity here.

These arguments are needed in addition to the smooth heat coefficient
identities in [Vassilevich, sections 4.1 and 7.1](https://arxiv.org/html/hep-th/0306138v3#S4.SS1).
That reference alone does not prove (4)-(5) on this paired seam.

## 3. Remove the common regulator and compute all conversion contacts

The full angular coefficients in the retained convention are

 Theta(e)=2C/e+H0+H1 e+O(e^2),
 H0=-q-1/3, H1=(q^2/6-1/30)/C.

Use the differentiated angular expansion as well. In (3) the geometric
remainder is multiplied by O(e^-1), and the insertion remainder by
O(e^-2); the preceding estimates make both errors tend to zero uniformly
on bounded admitted directions. Thus this is common regulator removal
for the first two jets at the actual seed. There is no angular/thermal
truncation. The finite-cutoff identity also retains every large-time term:
their difference has already become the endpoint trace (2).

Let J denote this common finite part. Its first derivative is

 J'[u]=(4pi)^-1 integral u[2H1-H0 R/6-C R^2/60-C Delta R/30]. (8)

Vary (4), using the fixed-Hilbert-space insertion:

 delta A0[u]=4 integral uv, delta A1[u]=integral u Delta v/3,
 delta A2[u]=1/60 integral u[vR^2+2R Delta v+2Delta(Rv+Delta v)].

Substitute these and (5) in (3). The COMPLETE result is

 J''[u,v]=(4pi)^-1 integral {
   [8H1+C R^2/15+2C Delta R/15]uv
   -[H0/3+C R/5]grad u.grad v-C(Delta u)(Delta v)/15 }. (9)

All boundary integrations used here are legitimate on the actual phase:
u,v=O(y^3), their first derivatives O(y^2), and second derivatives O(y).
Every fourth-order boundary pairing vanishes. In particular the seam
delta in Delta R pairs with uv(0)=0. It was retained to establish that
fact. At the regular bolt there is no puncture boundary term. No extra
Robin condition or mass-dependent projector is introduced.

The small-cutoff difference has only inverse powers and a finite constant,
with no logarithmic divergence. The common zeta/Mellin prescription and
cutoff finite part therefore give the same difference: their universal
Euler-constant/scale correction multiplies the equal logarithmic coefficient.
No finite coupling is adjusted. For u=v=1, (9) recovers exactly the previous
constant test J''=(A0 H1-2C A2)/pi after the paired divergence completion.
Constants are only a consistency probe; the physical directions still
have vanishing outer traces.

After multiplication by -2G/B^2 the local strong conversion, on the
physical side and in the original dz measure, is

 Hconv v=fq { Delta^2 v/30-(R Delta v+Q R'v')/10
             -H0 Delta v/(6C)
             -[4H1/C+R^2/30+Delta R/15]v },
 fq=GC/(pi B^2), Delta v=(Qv')'.                       (10)

This identifies the seed mass/geometric conversion completely. It does
not prove a full variable-background conversion theorem in the strong
graph, nor a Z2 estimate or the admission of the whole nonlinear X ball.

## 4. Finite fourth contact: preserve the existing normal, do not add twice

The flat fourth conversion in (9) is -C k^4/(60pi). In the normalized
shape row its normal is

 Bconv(s)=K(s-4)s(s-1)(s-2)(s-3)/2,
 K=GC beta^2/(15pi B^2).                              (11)

This is nonzero, but does NOT mean the retained geometric normal is wrong.
Check it against the actual geometric heat vertex including its second
vertex. With a=alpha(1-alpha),

 F(x)=1/2+1/2 integral_0^1(1+ax)exp(-ax)dalpha,
 Rg(x)=F(x)-1+x^2/120,
 Rm(t)=log(1+t)-t+t^2/2,

 integral_0^infinity exp(-M^2 t) Rg(t k^2)dt/t^2
        =M^2/2 integral_0^1 Rm(a k^2/M^2)dalpha.       (12)

Differentiate each side in a k^2 before the alpha integral: both derivatives
are q/M^2-q/(M^2+q), with the common factor 1/2; both values vanish at q=0.
This proves (12), including its finite part. Positivity and the small-q
cubic bound dominate the angular sum by a constant times
k^6 sum d_n/M_n^4, which is finite. Restore the SAME a2 finite coefficient.
The result is exactly the existing Kflat, with -31/30 and hence the full
-6/5 after the Weyl term. The finite normal and slow-root definition do
not change. This is an equality of comparison multipliers justified by
the geometric heat integral; it is not an identification of the actual
curved determinant with a flat domain.

In a mass representation the normal must be Nmass=Ngeom-Hconv,0, so that
Hgeom-Ngeom=(Hmass-Nmass)+(Hconv-Hconv,0). Adding (11) to the already
geometric Ngeom would count it twice. The local a0/a2/LLL/Weyl blocks
already allocated must also be subtracted exactly once in either ledger.

## 5. Quantified conversion remainder and the remaining actual requirement

Use the actual degree-six Q, correlated z_h, sigma=4.005, a=.01, the
original output cutoff and the full normal graph from the preceding notes.
Apply (10) to z^s, multiply by z^4, then apply the original (D-4) shape
row. Its constant-in-z polynomial is exactly (11). All remaining terms
are sum_(p>=2) z^p c_p(s). They have s-degree <=5. No log or finite term
is discarded. The scalar row is zero in the fixed canonical variables.

Use the already proved whole-frequency constants C_b,j for
||partial_x^j B_sigma^-1||_(H3->H3), j=0,...,5. For c_p(sigma+t)=sum c_pj t^j,

 eta_conv <= fq sum_p a^p (p^3+12p^2+51p+100)
                                    sum_j |c_pj| C_b,j.       (13)

The same cutoff multiplication bound and unchanged physical normalization
are used as in current-local-quantum-contact-graph.md. Freeze acceptance
eta_conv<8.36e-5 and .18-eta_conv>.1799164. This is a strong bound for
the explicitly computed LOCAL conversion difference only.

For a converted-mass implementation of the canonical remainder, subtract
the already allocated geometric local terms in the same ledger. The sum
of ALL remaining mass-kernel terms must then have graph norm <.1799164
to fit the existing .18 allocation by triangle inequality. This is an
UNKNOWN estimate, not an excessive upper bound or a measured error.
The selected extra global allowance .00368, joint affine <=7.49e-12,
known-source <=1e-14 and independent Mellin Efar<=1e-10 remain unknown.
No conversion allowance is added to the previous combined .18368 budget.
Identical-collar domain-transfer conversion terms cancel locally, and
separated-input conversion terms vanish by locality; the preceding
mass-kernel transfer bounds therefore remain usable with this accounting.

The old action-form kernel envelope cannot be differentiated into (13)
for the nonlocal part. The first unsupported step is now the actual
near-diagonal mass-kernel strong estimate with its fixed contacts and
common subtraction. All-mode Riccati/Dirac identities are available, but
their required simultaneous/output derivatives and weighted cancellations
have not been bounded. An existing finite weak-form bound is not such a
bound. No result here says the actual remainder exceeds .1799164.

## 6. Research decision and execution boundary

Retain the construction and matched principal subtraction. Use the
conversion (9)-(10) to make the radial mass-kernel representation exact
at the seed; use the geometric heat identity (12) to anchor its normal.
This sharpens the selected method without changing theory, domain,
source, finite couplings or normalization. The contact accounting is
tractable: its executed target consumes less than .047 percent of .18.
The remaining nonlocal and nonlinear conditions have no certified size.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): derive and evaluate the
complete coincident, moment-subtracted radial mass-kernel response in the
coupled H3 graph on the admitted fast inputs, with the now explicit local
conversion and the already allocated blocks counted once. Test the
remaining .1799164 allocation, including the symmetric reflected kernel,
diagonal jump/seagull contacts, every angular/thermal tail and original
global coupling. Pointwise envelopes or finite mode samples will not
certify that norm. Do not optimize the controlled cold/scalar components.

The sufficient projection enclosure E_rest<2.123e-42<1.287e-40 is preserved;
its zero enclosure center is not a measured signed projection, and this
turn does not reevaluate it. No actual controlled correction is computed.
At rho=1e-8, Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1 remain open;
the allocations Y0<=5e-9, Z0<=.25, Z2<=5e7 are still unproved. The full
X-ball admission, mass/slow matching, exact compatibility, physical and
causal matching, semiclassical validity/noise and all-sector stability
remain unresolved. This is neither existence nor nonexistence evidence.

Remaining to-do list: actual strong near-diagonal response and global data;
mass/slow matching and complete nonlinear admission/contraction; exact
compatibility, physical matching/causal realization, validity/noise,
all-sector stability, charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
