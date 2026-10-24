# A common-regulator test of the metric-to-mass conversion

2026-10-24, E003/S169. Ordinary argument frozen before acceptance.
Full stable-wormhole objective ACTIVE and INCOMPLETE; feasibility audit CLOSED.

## 1. Decision being tested, and its limits

The previous checkpoint leaves the canonical geometric heat Hessian open.
The radial mass-vertex resolvent is not that Hessian until its change of
measure, diagonal contacts and common subtraction are matched. Test the
proposed simplification of adding only the two-dimensional massless Weyl
cocycle. This note derives an exact counterexample to that simplification
using a scalar rescaling of the ORIGINAL paired seed operator. It does not
replace the physical domain by a smooth double, a flat torus or finitely
many eigenfunctions. It does not evaluate the missing strong response.

Constant b is an auxiliary consistency direction, not an admissible
fixed-source correction: the actual phase requires b(0)=b'(0)=0. The
constant multiplication nevertheless preserves the original first-order
projector and its complementary adjoint, bolt regularity and thermal spin
structure. The test can disprove a purported operator conversion on this
larger multiplier family. It cannot certify the target on the smaller
physical phase. The local volume coefficient it detects also has to be
matched for admitted directions supported inside the physical domain;
their derivative contacts require a separate calculation.

The seed-domain and bounded-admitted-family weak C2 arguments remain those
of current-rough-paired-variation.md and current-curved-paired-heat-bridge.md.
Neither supplies a uniform nonlinear solution ball or the remaining H3
graph estimate. The published smooth-surface calculation of
[Namuduri and Bilal, section 5.1](https://arxiv.org/html/2308.05802v2#S5.SS1)
also keeps massive local terms alongside the conformal anomaly. Its scope
does not establish our rough paired domain or our common angular limit.
The calculation below uses the repository's actual paired heat trace.

## 2. Exact all-spectrum identity at a common cutoff

Let L be the nonnegative canonical two-dimensional square with its
original composition domain. For n>=1 put M_n^2=n(n+q)/C and d_n=2(q+2n).
In the constant real direction b, the geometric and mass-chart squares are

 P_g,n(b)=exp(-2b)L+M_n^2,
 P_m,n(b)=L+exp(2b)M_n^2,
 P_g,n(b)=exp(-2b)P_m,n(b).                         (1)

This is an equality of closed operators with the original domain. It
requires no factorization of a regularized determinant. All massive
sectors have a gap; no zero-mode projector is silently omitted. Let

 Kbar(t)=gamma/(2pi) Tr_phys-half exp(-tL),
 Theta(t)=sum_(n>=1) d_n exp(-tM_n^2),
 Hbar_m(t;b)=Kbar(t)Theta(exp(2b)t).

Kbar includes the entire antiperiodic thermal spectrum and the actual
bolt. The paired completion includes the cubic seam. In time-reduced
units the parity-even fermion action is +1/2 integral Hbar dt/t. Hence
at the SAME common proper-time cutoff epsilon,

 Gbar_g,epsilon(b)-Gbar_m,epsilon(b)
       =1/2 integral_(epsilon exp(-2b))^epsilon Hbar_m(t;b) dt/t. (2)

The sign also follows from Gamma=-logdet(P)/2. The change of variable
in (2) is exact before any limit or angular summation interchange.
Positive massive heat traces justify the combined sum for epsilon>0.
No level-dependent finite subtraction is allowed.

Use the already established paired, integrated coefficients

 Kbar(t)=(A0/t+A1+A2 t)/(4pi)+O(t^2),
 A0=2 z_h,
 A1=-integral R/6=(Q'(z_h)-Q'(0))/6=-gamma/3,
 A2=-integral R^2/120,
 R=-Q'',   Q'(0)=0, Q'(z_h)=-2gamma.                (3)

In particular A2 is the integrated PAIRED completion, not the bare
half-disk bulk term with its curvature-divergence boundary erased. The
prior weak rough-seam argument supplies this expansion without a C4
reflected metric. The constant direction merely rescales its heat time.
Thermal images and seam remainders remain in Kbar.

The complete massive angular heat coefficients in the same convention are

 Theta(t)=2C/t+H0+H1 t+O(t^2),
 H0=-q-1/3,   H1=(q^2/6-1/30)/C.                  (4)

Multiplying (3) and (4), write Hbar_m=a_-2/t^2+a_-1/t+a_0+O(t).
The three coefficients are

 a_-2=A0*2C exp(-2b)/(4pi),
 a_-1=[A0 H0+A1*2C exp(-2b)]/(4pi),
 a_0=[A0 H1 exp(2b)+A1 H0+2C A2 exp(-2b)]/(4pi). (5)

Substitution in (2) gives

 a_-2(exp(4b)-1)/(4 epsilon^2)
 +a_-1(exp(2b)-1)/(2 epsilon)+b a_0+O(epsilon).    (6)

Thus the common finite part of the conversion is exactly

 Jbar(b)=b[A0 H1 exp(2b)+A1 H0+2C A2 exp(-2b)]/(4pi). (7)

Equivalently, where the common zeta trace initially converges,
zeta_g(s;b)=exp(2bs)zeta_m(s;b); its derivative at zero gives b zeta_m(0;b).
Equation (6), rather than an assumed determinant-product law, justifies
that identity's finite part on this actual rough paired family. Its
O(epsilon) tends to zero uniformly for b in a bounded constant interval,
as do its first two b derivatives: use the same positive-time analytic
heat trace and the admitted-family derivative majorant. This proves
common regulator removal for this test, not for arbitrary variable b
in the strong norm needed by the correction equation.

## 3. Executed quantity and acceptance conditions

The retained normalized action multiplies Gbar by -2G/B^2, with G=1.
Do not change that normalization. Differentiate (7) twice at zero:

 Jbar''(0)=(A0 H1-2C A2)/pi,
 j_conversion=-2G(A0 H1-2C A2)/(pi B^2)
             =-4G[z_h H1-C A2]/(pi B^2).           (8)

Evaluate (8) using the correlated actual parameters, including the
latest finite geometry, the true z_h=1/(1+h), and

 Q=(1-z)^2[beta(1+2z+3z^2)+u z^2+eps z^4]-tau z^3. (9)

The R^2 integral is an exact polynomial primitive on [0,z_h], enclosed
by Arb. Keep the two contributions, their signs and their certified
width. Acceptance checks the sign j_conversion<0 and agreement of
the cutoff and zeta derivations, not a guessed target magnitude.
The independent verifier computes integral (Q'')^2 by an interval
Bernstein product on [0,z_h] instead of the power-basis primitive.

The first run's additional arithmetic-width guard of 1e-110 failed:
inherited parameter enclosures give a radius about 1.425e-50. That run
and a diagnostic replay are retained. The replay changes only this
nonphysical width guard to 1e-48, retaining more than nineteen relative
decimal places. The exact identity and required nonzero sign are unchanged;
no response/projection target is relaxed and no parameter interval narrowed.

The purely two-dimensional massless cocycle is linear under a constant
rescaling (Euler/boundary completion fixed); its Hessian in this test
is zero. Equation (8) is nonzero: H1>0 and A2<0 on the actual seed.
Consequently that cocycle alone cannot convert the commonly regulated
mass action. The extra terms are fixed by the regulator, not selectable
forces or finite counterterms. Existing notes already reserved these
contacts; no previously certified block is being changed or enlarged.

This is a complete signed scalar consistency test with a numerical
enclosure. It is NOT the requested full signed source projection, an H3
operator bound, or an error bar for the actual correction. Its numerical
size cannot be compared directly with .18368 or 1.287e-40: those are
different functionals/norms. No new requirement for a controlled solution
is declared closed from that comparison.

## 4. First unsupported step and research decision

For a variable multiplier b the first-order factorization
e^(-b/2)(D+M exp(b))e^(-b/2)=Dtilde_b+M holds on the fixed paired domain.
It does not make the squares scalar multiples, as in (1). Their regulated
determinants need a variable-background cocycle. A full four-dimensional
Weyl path through a variable spherical radius is a possible derivation,
but the existing canonical product heat proof does not establish its
uniform rough-domain differentiated limit. This is the first unsupported
step in the proposed mass-chart route. It is UNKNOWN. The simplified
two-dimensional-only conversion is demonstrably false; this is an
obstruction to that shortcut, not to the construction or wormhole existence.

RESEARCH DECISION: retain the construction and use the original geometric
Duhamel response directly for the remaining estimate. Do not extend a
mass-vertex estimate chain until the common variable conversion is proved.
The geometric representation already has the correct closed operator,
composition domain, second vertex and weak admitted-family regulator
argument. The new test supplies a reason to choose that representation;
it does not supply its missing quantitative strong bound.

Write H_g,epsilon for that geometric Hessian, with the full common a0,
a1, a2 subtraction/restoration. Subtract the full normal and all pieces
allocated in current-local-quantum-contact-graph.md exactly once. Denote
the remaining coupled feedback, after its displayed row and variable
transforms and original global matching, by R_epsilon. The decisive
requirement is an operator-norm limit R and

 ||R||_(H3 direct-sum H3 -> H3 direct-sum H3)<=.18368. (10)

This includes other unresolved global feedback. A useful explicit split
would be canonical heat <=.18 and every remaining global term <=.00368,
but BOTH are proposed allocations, not measured bounds. Uniform weak C2
action or a finite matrix cannot be differentiated to establish (10).
The prior generic weak remainder upper bound is excessive for that use;
the actual strong remainder is UNKNOWN, not known to exceed (10).
The controlled logarithmic principal part removes the previously proved
H6-only high-frequency obstruction. There is therefore no current proof
of an obstruction to the full graph method, but also no numerical basis
yet for asserting that the remaining constants fit. Sharpening already
passing cold/scalar or explicit contacts cannot settle this question.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): apply the original
geometric Duhamel formula, including its second vertex, to the actual
cubic collar; subtract the already allocated normal/logarithmic principal
terms BEFORE taking three output derivatives and common regulator removal.
Prove a uniform operator-norm Cauchy estimate and bound the resulting
canonical remainder by .18, retaining angular/thermal tails, contact
distributions and original global coupling. This is the calculation that
can determine whether the remaining graph method fits its allocation.

The prior sufficient projection screen E_rest<2.123e-42<1.287e-40 is
unchanged and not recomputed here; its zero enclosure center is not a
measured signed projection. Joint affine <=7.49e-12 and known-source
<=1e-14 remain UNKNOWN. Together with (10) they would imply joint
E<9.946377e-12, still conditional. Independent Mellin Efar<=1e-10 and
actual mass/slow matching coefficients remain unknown.

At rho=1e-8, preserved-ball and contraction still require

 Y0+1e-8 Z0+5e-17 Z2<=1e-8,   Z0+1e-8 Z2<1.

Y0<=5e-9, Z0<=.25, Z2<=5e7 are unproved allocations. No uniform nonlinear
X-ball admission, actual controlled solution, exact compatibility,
physical matching/causal realization, semiclassical validity/noise or
all-sector stability is established. The closed audit stays closed.

Remaining to-do list: direct geometric strong remainder and common limit;
global feedback and joint sources, mass/slow moments and independent
Mellin far tail; actual controlled correction and nonlinear admission;
exact compatibility, physical matching/causal realization, validity/noise,
all-sector stability, charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
