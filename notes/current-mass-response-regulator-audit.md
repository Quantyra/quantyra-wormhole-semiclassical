# Complete third mass response: regulator inventory and a cancellation test

2026-10-24, E003/S169. Ordinary argument before frozen acceptance.
The stable-wormhole objective remains ACTIVE and INCOMPLETE; audit CLOSED.

The previous seed conversion is retained. The missing result is still a
strong, uniform bound on the nonlocal metric response. This note executes
the preceding regulator test: derive the complete positive-cutoff third
derivative, identify its common ultraviolet contacts, and test the signed
combination using the entire angular tower. It does NOT execute the full
strong feasibility calculation or certify an actual correction. This
distinction matters more than the number of identities checked below.

## 1. The actual operator, rather than a formally similar expression

Here D denotes the SKEW-ADJOINT representative Dsharp=J D0 on the original
paired first-order domain. The unpaired covariant D0 itself has the
complementary adjoint domain described in current-rough-paired-variation.md.
These statements are not interchangeable. Keep both identical copies,
their projector, the regular disk bolt and the thermal spin structure.
All scalar fields are identical on the copies. Write F=exp(b),

 A_M(b)=D+M F,  P_M(b)=A_M(b)^* A_M(b)
       =L+M^2 F^2+M[F,D],  L=-D^2.                    (1)

The mass/log phase has bounded F,F^-1 and dF. Multiplication by either
preserves Dom D, by the product rule and the unchanged scalar boundary
trace. For each M, the composed square in (1) consequently has Dom D^2:
D psi+MF psi is in Dom D exactly when D psi is. This is a deduction from
composition for the MASS representation, not a new Robin condition and
not a claim that the geometric square's domain is fixed. The all-angular
direct sum is treated with its closed form; M-dependent vertices are not
called bounded on the entire unweighted direct sum.

For real b with ||b||infinity<=B, accretivity of A_M and its adjoint gives
P_M>=M^2 exp(-2B). A bounded admitted W3,infinity mass family supplies
the weak common-regulator argument of current-third-metric-conversion.md.
Its existing conclusion is an integrated weak third derivative. The
five-dilation-jet W3,infinity method remains obstructed on the fast graph;
no such method is reintroduced here. The original X-ball admission and
the uniform STRONG graph limit do not follow from the weak argument.

For directions u_i, the exact multiplication vertices are

 P_I=2^|I| M^2 exp(2b) product_(i in I) u_i
                        +M[exp(b) product_(i in I) u_i,D]. (2)

In a Clifford convention [f,D]=-c(df), the last sign is minus. In
particular P_u=2M^2u-Mc(du), P_uv=4M^2uv-Mc(d(uv)), and
P_uvw=8M^2uvw-Mc(d(uvw)) at the seed. Dropping any gradient vertex
would change the operator. The skew representative's Clifford matrices,
not a separately selected unpaired matrix representation, are understood.

## 2. Exact common-cutoff third derivative

Let K_s=exp(-sP_M(b)), and use the original physical-half, time-reduced
trace. For t>0 put

 T3_M(t;u,v,w) = -t Tr(P_uvw K_t)
  +t^2 integral_0^1 sum_pairs Tr(P_uv K_(t a) P_w K_(t(1-a))) da
  -t^3 integral_(a,c>=0,a+c<=1) {
     Tr(P_u K_(t a) P_v K_(t c) P_w K_(t(1-a-c)))
    +Tr(P_u K_(t a) P_w K_(t c) P_v K_(t(1-a-c)))} da dc. (3)

The pairs are (uv,w),(uw,v),(vw,u). The triangle has area 1/2;
both cyclic orientations are present. It is not three independent
resolvents with arbitrarily inserted heat factors.

To derive (3), differentiate the Duhamel series three times. Before
cyclicity the triple term has six orderings and four heat intervals.
Combine the first and last intervals under the trace. The corresponding
length weight, summed over the three cyclic rotations of one ordering,
is a+c+(1-a-c)=1. This leaves two orientations on the triangle. The
same cyclic grouping reduces the pair term to the displayed interval.
The term with P_uvw comes from differentiating the operator itself.
For commuting vertices the coefficients become -tP3+3t^2P2P1-t^3P1^3,
but commutativity is not used in the derivation.

For every positive common cutoff e, the ACTUAL mass action satisfies

 D^3 Gamma_m,e[u,v,w]
       =1/2 integral_e^infinity sum_(n>=1) d_n T3_Mn(t;u,v,w) dt/t. (4)

All thermal modes, reflections and finite-domain propagation are inside
K. At fixed positive cutoff the mass gap and heat smoothing justify
the sums and derivatives; vertices only add finite mass powers. The
positive-time formulas on the admitted rough family follow by the
fixed-domain form limit. There is no thermal truncation or auxiliary
wall. At t=0, only the COMBINED expression can be passed to the common
finite part. Individual terms of (3) are not separately renormalized.

## 3. Common contacts, including a surviving cubic derivative term

Write the integrated heat trace as a0(b)t^-2+a1(b)t^-1+a2(b)+O(t).
These are four-dimensional coefficients (a2 here is often named a4).
With R=-Q'', Delta=(Q partial_z)', H0=-q-1/3 and
H1=(q^2/6-1/30)/C, their integrated representatives are

 a0(b)=C/pi integral exp(-2b) dz,
 a1(b)=integral {H0/(2pi)
       -C exp(-2b)[R+2Q(b')^2]/(12pi)} dz,
 a2(b)=1/(4pi) integral {2H1 exp(2b)-H0 R/6
             -C exp(-2b)[R-2Delta b]^2/60} dz.        (5)

Boundary completions are understood before variation. Here is how they
enter, rather than assuming a smooth boundaryless surrogate. The angular
isometry identifies the mass metric g2+exp(-2b)gS. Its volume gives a0.
Its scalar curvature is R+2exp(2b)/C+4Delta b-6|db|^2; the integral
of exp(-2b)Delta b is 2 integral exp(-2b)|db|^2 plus its endpoint flux.
This gives a1. For a2 use the integrated conformally invariant heat
coefficient and the product-frame orbit metric exp(2b)g2, whose curvature
is exp(-2b)(R-2Delta b). The zero-angular subtraction is q times the
FIXED orbit heat trace and has no b derivatives in this mass chart.
The retained anomaly convention fixes the curvature-divergence completion;
no finite fourth contact, scale or coupling is changed.

For the THIRD variations used here, every outer flux is at least cubic
in admitted directions and has total derivative order at most three.
Since each direction is O(z^3) with the fixed mass/log remainder, these
fluxes vanish at least as z^6 times fixed log powers. Retain the paired
seam distributions first; their direction traces vanish as well. At
the bolt the regular disk has no boundary. This justifies the third
variation of (5) on the actual domain. It does not discard a possibly
nonzero first-order boundary flux or claim a generic boundary theorem.
The earlier weak ordered-heat argument, not the smooth coefficient
algebra alone, supplies the extension to the mass/log family.

Let A_j=D^3 a_j(0)[u,v,w]. In particular

 A0=-8C/pi integral uvw,
 A1=2C/(3pi) integral {Ruvw+u grad v.grad w
                         +v grad u.grad w+w grad u.grad v}. (6)

The cubic Taylor density of a2 is

 1/(4pi){(8H1/3)b^3+(C R^2/45)b^3
                   +(2C/15)R b^2 Delta b+(2C/15)b(Delta b)^2}. (7)

Its polarization gives A2. Unlike the finite geometric-to-mass
CONVERSION, the mass action's logarithmic heat contact has a surviving
b(Delta b)^2 term. The conversion's cancellation cannot be transferred
to this different quantity. This coefficient determines part of the
matched principal subtraction; it is not an extra physical interaction.

For a fixed proper-time split T>0 the ultraviolet part of the third
derivative must contain the convergent combination

 1/2 integral_0^T [sum d_n T3_Mn(t)-A0/t^2-A1/t-A2] dt/t. (8)

The exact large-time integral, the fixed finite restorations and the
certified conversion are added to (8). Equivalently the divergent part
of the raw cutoff expression (4) is

 A0/(4e^2)+A1/(2e)-A2 log(e)/2.                     (9)

Subtract (9) with the SAME regulator. A matched strong principal
subtraction is still needed: (8), which is a weak functional identity,
does not evade the strong cutoff obstruction already established in
current-geometric-cutoff-graph-test.md.

An executed ACTUAL admissible-direction test is u=v=w=z^5 on
0<z<z_h, with the correlated finite-seed Q and z_h. It preserves the
outer traces and is regular at the disk bolt; a sufficiently small
multiple is in any specified positive directional neighborhood.
All polynomial integrals in (6)-(7) were enclosed, giving

 A0=-7.817640236496421...e9,
 A1=-6.569445576871031...e8,
 A2= 6.912334259590839...e13.                         (10)

The archived intervals are authoritative. These numbers are heat
contacts in the ORIGINAL time-reduced action units. They are neither
the bolt projection nor its error, and must not be compared directly
with 1.287e-40 or Z2. They demonstrate nonzero contacts on an actually
admissible direction, rather than using a nonadmitted constant field
as a physical counterexample.

## 4. Execute the complete-angular cancellation diagnostic

For the auxiliary constant direction the mass heat trace factorizes
EXACTLY as Tr_bar exp(-tL) Theta(t exp(2b)). This is a coefficient
diagnostic in the larger bounded-multiplier family, not a substitution
for the physical mass/log phase. The exact orbit trace retains its
boundary, bolt and all thermal frequencies; it cancels from the ratio
below and is never replaced by selected orbit modes. The angular factor
of its third derivative is

 8t Theta'(t)+24t^2 Theta''(t)+8t^3 Theta'''(t),
 Theta(t)=sum_(n>=1) 2(q+2n)exp[-t n(n+q)/C].         (11)

The three terms are respectively the single contact, all three pair
contacts, and both triangle orientations. Relative to 2C/t, their
ultraviolet limits are -8,+48,-48. The combined limit is -8: the
triangle alone overcounts by SIX, and the sum of absolute values loses
a factor THIRTEEN. The actual local a0 term in (6) shows the same leading
coefficient for admissible directions. This is not an accidental
cancellation in a finite selection of modes.

For numerical evaluation of (11), put a=qt/C, c=t/C, r=exp(-a).
Expand only exp(-c n^2) to degree J=12. Every resulting infinite sum is
computed exactly by S_m(r)=sum_(n>=1)n^m r^n, with
S_0=r/(1-r), S_(m+1)=r partial_r S_m. For the kth positive mass moment,
the omitted Taylor tail is bounded by

 c^(J+1)/(J+1)! sum_(n>=1)2(q+2n)
                [n(n+q)/C]^k n^(2J+2) r^n.         (12)

This is a convergent complete sum, evaluated by the same rational
functions. No angular cutoff is imposed. The bound is a fixed-t
certificate; it is not used as a uniform t->0 estimate. Equation (5)
and the prior weak remainder supply the latter within their stated scope.

At the actual C and q, interval evaluation gives the combined quantity
in (11), divided by 2C/t, in

 t=1e-6: [-7.9999999292644740,-7.9999999292644738],
 t=1e-8: [-7.999999999992926448,-7.999999999992926446]. (13)

At t=1e-6 the triangle/complete ratio is 6.000000053051645...,
and the absolute-sum/absolute-complete ratio is 13.000000106103289....
The archived errors include the complete tail (12) and interval
uncertainty in the actual parameters. These ratios are a regulator
inventory test, NOT the nonlocal strong norm or the physical projection.

## 5. What this decides, and the numerical target that remains

The closed first-order realization and weak combined derivative are
available on bounded admitted families. The first unsupported continuum
step is still the uniform strong estimate for (3) after its matched
principal/local subtraction, including finite-domain returns. Neither
the seed conversion nor a trace expansion proves it. A theorem about
smooth closed-manifold spectral Hessians, or scalar Schrödinger trace
regularity, does not supply it without these hypotheses. For example,
the scopes stated by [Okikiolu--Wang](https://arxiv.org/abs/math/0211101v2)
and [Smith](https://arxiv.org/abs/1809.05614v1) are not this paired
mass/log operator estimate; only their abstracts were consulted here.
No external theorem from those papers is asserted as a substitute.

Classify the evidence precisely:

* The separate raw heat-word integral is ultraviolet divergent. Equations
  (6),(9)-(13) expose a failure of that ESTIMATOR, not divergence of the
  common-subtracted action. Only the combined expression is admissible.
* The complete strong remainder constant is UNKNOWN. There is no finite
  certified full bound here that can honestly be labeled excessive.
* No obstruction to the actual quantum operator, controlled solution or
  wormhole existence has been established. The earlier strong raw-cutoff
  and uniform W3,infinity proof-method obstructions remain unchanged.

In the coupled fast graph let C3 bound the complete nonlocal bilinear
response into H3, uniformly on the admitted neighborhood. If eta_seed
is its seed linear feedback and r_graph the neighborhood radius, the
conditional nonlocal screen is

 eta_seed+C3 r_graph < .1599164.                    (14)

At r_graph=1e-11 this requires
C3 < 1.599164e10-1e11 eta_seed. Both constants are unknown. Even the
best-case eta_seed=0 ceiling 1.599164e10 is a CONDITIONAL target, not
an achieved bound. A seed allowance .08 would leave C3<7.99164e9.
This illustrates sensitivity without assigning the unknown seed norm.

For the original nonlinear X norm, let alpha bound the X-to-input-graph
map and beta bound the ORIGINAL preconditioned output return H3-to-X.
With the other sectors' second-derivative bound Z2_other, certification
requires

 beta alpha^2 (C3+C3_conversion+other graph terms)+Z2_other<=5e7. (15)

The conversion value <9.841e-33 is ONLY its seed outer fast value; it
cannot be inserted as a uniform global C3_conversion without proof.
Alpha, beta, global/slow contributions and Z2_other remain unknown.
Setting these comparisons to one, or changing normalization, would not
constitute an improvement. The thirteenfold raw cancellation removes an
incorrect decomposition; it gives no justified proximity to (14)-(15).

The seed contact targets M_A<=3.18e27 and M_V<=2.66e30 remain unproved.
The preserved-ball requirements at rho=1e-8 remain
Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1, with unproved allocations
Y0<=5e-9, Z0<=.25, Z2<=5e7. Preserve the earlier sufficient projection
enclosure E_rest<2.123e-42<1.287e-40: zero was its enclosure center,
not a measured signed value. No new full signed projection, actual global
correction or nonlinear solution is executed here. Exact compatibility,
physical matching/causal realization, semiclassical validity/noise and
all-sector stability remain open. The closed audit is not reopened.

RESEARCH DECISION: retain the construction, and use the complete
same-regulator Duhamel combination (3),(8) for the direct Sobolev method.
Reject triangle-only and separately absolutized ultraviolet estimates.
Do not optimize cold/scalar terms, the seed tube or the controlled
conversion contact. The cancellation test is executed; the decisive
strong feasibility bound is still NOT executed.

SINGLE NEXT CALCULATION (PROPOSED): bound the first metric variation of
the common-regulated, moment/principal-subtracted mass Hessian from the
coupled fast graph to H3, uniformly on the actual admitted mass/log
family, using (3) with all contacts combined BEFORE norms. Keep the
highest spatial derivative in L2 and include the entire angular/thermal
sums and finite-domain return. A finite explicit C3 and the original
norm comparison must test (14)-(15). Until then neither continuing to
existence nor switching construction is numerically justified; the
selected change is to the estimator, not the field theory or domain.

Remaining to-do list: complete strong nonlocal response and seed quantum
norms, actual global correction, original X/graph comparison and nonlinear
admission/contraction, mass/slow matching, exact compatibility, physical
matching/causal realization, semiclassical validity/noise, all-sector
stability, charge/return, opening, transport, preparation, separate
creation and the full stable-wormhole objective.
