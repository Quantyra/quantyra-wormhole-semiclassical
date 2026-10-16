# A global range condition for the stationary retained trace

2026-10-16, S139/S074. Ordinary argument independently reviewed before
focused verification.

The accepted centered candidate has a small complete equation defect.
That result does not yet give a stationary solution error. This note
isolates an additional global condition that any sufficiently decaying
exact solution of the selected full free-Dirac equations must obey.
The trace equation has a radial oscillatory channel at a proper
Planck momentum. Regularity at the nonextremal bolt and the accepted
finite-mass AdS asymptotics do not give an unrestricted scalar inverse.

The result is a necessary condition on the actual nonlinear metric,
and an exact identity for the candidate's projection on that condition.
It is not nonexistence of a coupled solution, a physical instability,
or a failure of the established equation-defect estimates. It adds
neither a finite coupling nor a spatial cutoff.

## Fixed action, scale and dimensions

Use the complete physical action at bookkeeping parameter one, with
the fixed magnetic flux, laboratory clock, prepared temperature,
paired Dirac domain and zero additional finite curvature coefficients
of [the centered seed](stationary-centered-covariant-seed.md).
In particular the common representative has

```
16pi^2 A4=c C4^2-a E4+b_an Box R4-(2/3)F4^2,
c=1/20, a=11/360, b_an=1/30,
R4-4Lambda+8piG A4=0.                              (1)
```

These are the bulk trace and normalization in
[the solver obstruction](stationary-exterior-solver-obstruction.md).
The independent Weyl-invariant determinant has no compact-variation
four-dimensional trace. Its nonlocal traceless stress is still present
in the other equations. A boundary variation is not used to change
the coefficient in (1).

Write Lambda=-3lambda0 and let lambda_e be the selected positive
outer-vacuum root

```
lambda_e+11G lambda_e^2/(360pi)=lambda0.
R_*=-12lambda_e,
epsilon=G/(60pi),
a_inf=1-Ga R_*/(6pi)=1+11G lambda_e/(180pi).         (2)
```

Here lambda_e has dimension length^-2, epsilon has dimension
length^2, and a_inf is dimensionless. The original completion's
bookkeeping parameter is not lambda_e; equations (1)-(2) concern the
full selected action. No new vacuum, reference radius or renormalization
scale is chosen.

## Exact nonlinear scalar identity

Put

```
Z=R4-R_*,
J_ab=Ric_ab-(R4/4)g_ab,
J^2=J_ab J^ab.
```

The algebraic identity E4=C4^2-2J^2+R4^2/6 holds for any
four-metric. Subtract (1) at its constant-curvature vacuum. This
gives the exact identity

```
P_g Z=S[g],
P_g=a_inf+epsilon Box_g,
S[g]=G/(2pi){(2/3)F4^2-(7/360)C4^2
                    -(11/180)J^2+(11/2160)Z^2}.    (3)
```

Every term in (3) has dimension length^-2. The right side contains
curvature but no derivative of curvature. No sign or positivity is
assumed for J^2 or for S[g]. This identity includes both physical
orbit pressures and the radius equation through the full trace.

In particular it does not set the centered seed curvature constant.
On that seed

```
B^2(R40-R_*)=4d/t^5-6d/t^6, t=r/B.
```

Linearizing (3) at any fixed seed gives

```
P_0 delta Z=delta S-epsilon(delta Box)Z0.           (4)
```

The last term is retained when d is nonzero. Replacing the seed by
classical RN before this variation would incorrectly remove it.
For a static scalar f in areal gauge,

```
Box_g f=F f''+[F'+F alpha'+2F/r]f',
(delta Box)Z0=delta F Z0''
 +[delta F'+delta F alpha'+F delta alpha'
                              +2delta F/r]Z0'.    (5)
```

Thus the extra term is explicit; it is not an unknown state response.

## Class of stationary geometries

Consider one exterior with

```
ds^2=-exp(2alpha)F dU^2+dr^2/F+r^2dOmega^2,
r_h<r<infinity, F>0,
F(r_h)=0, F'_h>0,
kappa=exp(alpha_h)F'_h/2>0,
alpha(infinity)=0.
```

Use the same positive prepared temperature as in the accepted
construction; no extremal limit is taken. Assume ordinary regularity
at the Euclidean bolt. For the asymptotic conclusions assume the
accepted finite-mass falloffs

```
F=lambda_e r^2+1-2m/r+O(r^-2[1+log(r/B)]),
alpha=O(r^-4[1+log(r/B)]),
Z=O(r^-4[1+log(r/B)]).                             (6)
```

The O terms have the differentiated bounds required to apply the
operators below: two derivatives of the metric remainders for its
Liouville potential, and two derivatives of Z for P_g Z. The actual
candidate's real-C5/open-collar conormal estimates provide these.
An asserted value asymptotic alone would not suffice. The metric has
at least the ordinary local regularity needed for the second-order
scalar equation; near infinity its stated logarithms are allowed.

For the magnetic bundle F4^2=q^2/(2r^4). Equations (3),(6) give

```
S[g]=Gq^2/(6pi r^4)+O(r^-6[1+log(r/B)]^2).         (7)
```

Indeed the AdS Weyl tensor has no constant part, its leading mass
term is order r^-3, J_ab is order r^-4 up to the retained logarithms,
and Z^2 is order r^-8. The exact source need not equal its leading
magnetic term on the finite exterior.

As one check of the scale and sign, its required leading curvature
coefficient is

```
Z=Gq^2/[6pi(a_inf+4epsilon lambda_e)] r^-4
                         +O(r^-6[1+log(r/B)]^2).  (8)
```

Here Box(r^-4)=4lambda_e r^-4+O(r^-6 log r).
To justify the remainder without differentiating an unspecified O
term, subtract the displayed multiple of r^-4 from the scalar
equation. Its forcing is O(r^-6 log^2 r). The fast inverse
constructed below gives precisely that remainder; the assumed
fast falloff rules out an added homogeneous solution.
Equation (8) is an asymptotic consequence, not construction of the
global metric or a bound on a solution's distance from the candidate.

## The actual regular homogeneous scalar

Define

```
w=exp(alpha)r^2, p=wF,
L_g=-w^-1 partial_r(p partial_r),
P_g=a_inf-epsilon L_g.
```

Let phi_g be the real homogeneous solution selected at the actual
nonextremal horizon by

```
P_g phi_g=0,
phi_g(r_h)=1,
phi'_g(r_h)=-a_inf/(epsilon F'_h).                 (9)
```

This is not an additional quantum state. It is a test solution of
the scalar differential expression determined by the metric.
Its regular existence follows directly from the Volterra equation

```
phi(r)=1-(a_inf/epsilon) integral_rh^r dt/p(t)
                                    integral_rh^t w(s)phi(s)ds.
```

Since p(t)=w_h F'_h(t-r_h)+O((t-r_h)^2), the apparent division
is removable. On a sufficiently short interval it is a contraction,
and the regular solution then continues uniquely along F>0. The
size of that starting interval need not be uniform as kappa tends
to zero; that limit is not being asserted. In proper bolt radius
rho, phi=1-a_inf rho^2/(4epsilon)+O(rho^4), so the large areal
derivative in (9) is compatible with smooth bolt regularity.

## Liouville reduction and the asymptotic channel

Let rho be physical proper radial distance, dr/d rho=sqrt(F), and
put

```
mu=w sqrt(F), psi=mu^(1/2) f,
V_L=.5(log mu)''+.25[(log mu)']^2.
```

Primes in this display mean rho derivatives. Direct conjugation
gives

```
mu^(1/2)L_g mu^(-1/2)=-partial_rho^2+V_L,
V_L=9lambda_e/4+O(r^-2[1+log(r/B)]),
rho=lambda_e^(-1/2)log(r/B)+constant+O(r^-2 log r).
```

The limiting proper momentum is real:

```
k_*^2=a_inf/epsilon-9lambda_e/4
     =60pi/G+17lambda_e/12>0,
omega=k_*/sqrt(lambda_e),
omega^2=60pi/(G lambda_e)+17/12.                   (10)
```

The perturbation of V_L from its limit is integrable in rho.
Variation of constants therefore yields, for any homogeneous
solution, two asymptotic coefficients in the basis

```
(r/B)^(-3/2) cos[omega log(r/B)],
(r/B)^(-3/2) sin[omega log(r/B)].                   (11)
```

The asymptotics hold with one differentiated remainder. For a
nonzero homogeneous solution the two coefficients are not both
zero: in the rotating constant-coefficient frame the evolution
matrix is integrable in rho, and its limiting fundamental matrix
and inverse are finite and nonsingular. In particular the pair
(A_phi,B_phi) of the regular phi in (9) is nonzero.

The forced equation with source (7) has the same two coefficients
and a fast particular remainder. After multiplying the source by
mu^(1/2), it decays as r^-5/2 and is integrable in rho. This
justifies the asymptotic coefficient and Wronskian limits below.
The accepted falloff (6) sets both oscillatory coefficients to zero.
These forbidden tails are curvature tails; there is no assertion
that each independent scalar homogeneous solution satisfies the
other gravitational or state equations.

At k=2000, the retained vacuum root has
.999 k^-10<G lambda_e<k^-10. Elementary pi bounds in (10)
give 4e17<omega<4.5e17, uniformly in the original preparation
range. This scale is a proper Planck momentum, not a new
macroscopic instability estimate or a justification for trusting
the retained theory at that momentum.

## What a scalar Hilbert-space inverse cannot supply

Take L_g with its nonnegative quadratic form
integral p|f'|^2dr in L2(w dr), regular at the bolt. Compactly
supported functions away from both ends belong to this domain,
independently of an admissible local bolt realization. At the
infinite end its Liouville form is a one-dimensional Schrodinger
operator with limiting potential 9lambda_e/4. Hence its essential
spectrum contains [9lambda_e/4,infinity). The inclusion needed
here has an elementary proof. In the rho coordinate take

```
psi_T(rho)=T^(-1/2) chi[(rho-R_T)/T] exp(i k_* rho),
||chi||_2=1, chi smooth and compactly supported in (0,1),
R_T -> infinity faster than T.
```

These normalized functions escape every compact set. Their
P_g residual is bounded by

```
epsilon{2k_*||chi'||_2/T+||chi''||_2/T^2
         +sup_support|V_L-9lambda_e/4|},            (12)
```

which tends to zero. Unitary transformation gives the same
statement in the physical weighted radial space. Thus zero lies
in the essential spectrum of this scalar P_g, and no everywhere
bounded L2 inverse exists on that space.

The scalar bolt endpoint is not limit-point. In proper radius
mu is a positive constant times rho to leading order and
V_L=-1/(4rho^2)+O(1). Both sqrt(rho) and
sqrt(rho)log(rho) are square integrable in d rho. Regularity
or the Friedrichs quadratic-form condition excludes the logarithm:
its scalar derivative contributes integral d rho/rho to the
energy. Mere L2 integrability would not select the physical
regular endpoint. None of these local endpoint choices removes
the Weyl sequence supported arbitrarily far out.

This does not rule out another weighted range problem, an
augmented inverse using coupled horizon data, or a nonlocal
selection justified by the physical variational problem. It does
not by itself rule out an inverse of the full coupled equations.
In particular (12) is not substituted for a counterexample in
the accepted mass/lapse observable norms.

## Exact nonlinear range condition and outward amplitudes

For any Z and phi with the indicated regularity,

```
w[phi P_g Z-Z P_g phi]
           =epsilon partial_r[p(phi Z'-phi' Z)].
```

At the bolt, regularity makes the right boundary concomitant
zero. For an exact solution obeying (6), its infinity limit is
also zero: p=O(r^4), phi=O(r^-3/2), phi'=O(r^-5/2),
Z=O(r^-4 log r) and Z'=O(r^-5 log r). Consequently every
such full solution must satisfy

```
C[g]:=integral_rh^infinity exp(alpha)r^2 phi_g S[g]dr=0. (13)
```

The integral converges by (7),(11). It has dimension length.
It is a nonlinear condition on the complete metric and magnetic
source, not a condition only on an independently prescribed trace.
It includes all the selected local trace terms and the same
matched outer vacuum. The HH matter condition alone is not a
proof of (13) on a trial metric.

For a fixed metric and fixed source S[g], consider the auxiliary
inhomogeneous scalar equation P_g Z_aux=S[g], REGULAR at
the bolt. Write its
oscillatory coefficients in the dimensionless radial basis (11)
as A_Z,B_Z. The regular homogeneous phi has coefficients
A_phi,B_phi in the same basis. The exact boundary identity gives

```
C[g]=epsilon lambda_e B^3 omega
                             (A_phi B_Z-B_phi A_Z). (14)
```

Here A_phi,B_phi are dimensionless and A_Z,B_Z have dimension
length^-2; the B^3 factor is necessary with basis (11).
With a dimensionful r^-3/2 basis it would instead be absorbed
in the amplitudes. The sign in (14) follows from
phi Z'-phi' Z=omega B^3 r^-4(A_phi B_Z-B_phi A_Z)+o(r^-4).

Changing the auxiliary regular horizon value shifts its two
coefficients only along (A_phi,B_phi). It cannot remove the
component measured by (14). If C[g]=0, one horizon value can
cancel the remaining parallel component in this FIXED linear
scalar problem. These statements are not a parameter-count
theorem for the full coupled problem: changing its horizon
geometry changes P_g, phi_g and S[g] simultaneously.

The condition is nonvacuous. For fixed P_g choose a compact
source supported in a sufficiently short interval where phi_g
has one sign. Its integral against w phi_g can be nonzero.
There is then no auxiliary scalar solution both regular at the
bolt and having (6)'s fast falloff. No such arbitrary source is
being substituted for the actual Dirac or magnetic source.

## Complete fixed-metric scalar range characterization

The fixed-coefficient scalar problem has an explicit fast inverse
on its compatible subspace. This assertion is separate from a full
gravitational inverse. Let f(r) be a prescribed real source with
regular bolt value and

```
f=O(r^-gamma[1+log(r/B)]^j), gamma>3/2,
```

with the ordinary differentiated bounds needed for any requested
solution derivatives. Let psi be any real independent homogeneous
solution normalized by

```
p(phi psi'-phi' psi)=1.
```

It exists globally on the open exterior. Zeros of phi do not
prevent continuation: define psi by independent Cauchy data and
the constant Wronskian, rather than using a reduction-of-order
formula across such a zero. Near the horizon
psi=log(r-r_h)/p'_h+O(1), where p'_h=w_h F'_h>0.
At infinity psi has the same oscillatory decay as phi.

Define

```
G_fast f(r)= [phi(r) integral_r^infinity psi(s)w(s)f(s)ds
            -psi(r) integral_r^infinity phi(s)w(s)f(s)ds]
                                                     /epsilon. (15)
```

Both integrals converge since their integrands are
O(r^(1/2-gamma)log^j r). Direct differentiation shows that
the endpoint terms cancel in the first derivative. In the second
derivative their Wronskian contributes w f/epsilon, so
P_g G_fast f=f with the sign displayed. The upper integrals are
O(r^(3/2-gamma)log^j r); multiplication by phi or psi proves
G_fast f=O(r^-gamma log^j r), with the corresponding first
derivative. This construction has no outgoing homogeneous tail.

Its horizon logarithm has coefficient

```
-[integral_rh^infinity w phi f dr]/(epsilon p'_h).   (16)
```

The other upper integral is finite at the horizon because log
is locally integrable. Thus the fast solution is regular iff
its compatibility moment is zero. When it is zero, the apparent
y log y terms cancel in (15); equivalently the regular
Volterra equation applies with its finite horizon value and
u'_h=[f_h-a_inf u_h]/(epsilon F'_h).

Any two fast solutions differ by a homogeneous solution whose
two oscillatory coefficients vanish, hence are identical. This
proves existence and uniqueness of a regular fast solution of
the fixed-metric scalar problem EXACTLY on

```
{ f : integral_rh^infinity w phi f dr=0 }.          (17)
```

No uniform numerical operator bound as h tends to zero, and no
uniform inverse constant for changes of metric, is asserted.
In (3) f=S[g] is itself a nonlinear function of the unknown
metric. Therefore (15)-(17) do not independently solve that
metric, the traceless equations or its physical state problem.

## Exact projection of the accepted candidate defect

For the original finite candidate at physical bookkeeping
parameter one put

```
E_trace[g_app]=P_gapp Z_app-S[g_app].
```

This is the geometric trace equation defect in the normalization
(1), with dimension length^-2. Use the actual candidate's phi,
not a homogeneous scalar on the unchanged classical RN metric.
The candidate is regular at its own bolt and has the strong
outer falloff, so the same boundary identity gives exactly

```
C[g_app]=-integral_rh^infinity
                 exp(alpha_app)r^2 phi_gapp E_trace[g_app]dr. (18)
```

No solution-error inverse was assumed in this formula. The
accepted [complete residual](stationary-centered-complete-residual.md)
controls stated mass, lapse and horizon norms. It neither asserts
that the oscillatory moment in (18) vanishes exactly nor supplies
a quantitative inverse for its complement. Even a very small
nonzero moment is not the exact condition C=0. Without a new
argument it cannot be converted into a nearby exact finite-mass
geometry or into an error bar on the renormalized charge.

## The branch decision still required

A spatially reduced equation can be defined by replacing the
higher derivatives using lower-order equations and consistently
truncating in the prescribed action expansion. The accepted
candidate and its complete bookkeeping-parameter O(ell^2) defect are useful inputs
to such a construction. Here ell multiplies Gamma_rest in the accepted
action expansion and has physical value one; this order statement is
not reinterpreted as a power bound in the Planck length. Formal truncation does not establish
that a solution of the full equations has been selected, that
the spectral moment (13) is satisfied, or that the full coupled
inverse is bounded in the required observable norms.

[Flanagan--Wald, section IV.4 and equations (115)-(122)](https://arxiv.org/html/gr-qc/9602052)
describe reduction of order as a modification accurate to a
specified retained order. Their limitations distinguish full
nonlocal equations and spatial derivatives from the time-order
reduction in their worked setting. Their equation (114) also
shows why smooth dependence on a small parameter by itself
need not select a unique exact prescription. Those results are
not an RN-AdS stationary branch theorem.

[Anderson--Molina-Paris--Mottola, sections II-IV](https://arxiv.org/html/gr-qc/0209075)
likewise distinguish low-energy semiclassical predictions from
additional Planck-scale solutions. Their flat-space linear
response analysis does not supply the boundary-value inverse
needed here. This note's global scalar condition uses the
specific retained coefficient and accepted AdS falloffs.

The next complete solver step must therefore specify and justify
the original action's gravitational boundary/branch prescription,
then establish a coupled weighted range/inverse result with (13)
included. If a reduced spatial model is used instead, it must
state its covariant constraint-compatible definition and justify
the discarded-order and branch-selection error for the actual
observables. An arbitrary bandlimit, changing b_an, or declaring
all fast modes absent does not supply that argument. The range
condition alone is a route correction, not a physical null or
completion of S139.1-S139.2.

## Focused verification contract

Before execution, the focused checker is restricted to exact algebra
and rational inequalities: anomaly decomposition and vacuum subtraction,
the seed curvature and delta-Box term, static scalar and Liouville
operators, critical horizon solutions, Green/Wronskian signs, the
dimensioned B^3 amplitude factor, candidate projection and the benchmark
frequency bracket. Every symbolic residual must vanish exactly and
every displayed rational inequality must be true; no numerical
tolerance is used. No shooting through the Planck-frequency channel
is planned or implied. Endpoint, spectrum, range and regularity
arguments above remain reviewed ordinary proofs. The checker cannot
establish a coupled inverse or actual zero moment. Execution outcomes
belong to its separate receipt and log.

Remaining to-do list: derive the physical gravitational boundary/branch
prescription and the coupled weighted inverse or a controlled
spatial reduction; complete exterior/state and charge matching,
then transfer the driven-source and remaining physical gates.
