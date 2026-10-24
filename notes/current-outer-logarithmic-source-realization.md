# The complete outer source and the logarithmic weight of its correction

2026-10-24, E003/S169. Ordinary continuum argument before frozen execution.
The stable semiclassical wormhole objective remains active and incomplete;
the feasibility audit stays closed. The previous turn reconstructed h_M,
not the full first correction. This addresses the boundary prerequisite
for reconstructing the remaining prescribed source and its strong return.

## 1. The first unsupported inference

A C2 action on admitted weighted mass/log families supplies weak metric
variations. It does not identify their pointwise endpoint expansion or a
strong source norm. The earlier fixed-seed estimate only bounded the
massive trace after its cubic pole by O(log(1/z)); that estimate did not
determine its logarithmic coefficient. It cannot be used as if that
coefficient were zero. We determine it on the ACTUAL paired realization.

The result will also test a possible single-logarithm source weight. The
repository has not established such a strong inverse; we do not claim
that an accepted inverse theorem used that weight. In particular the
already rejected second reference-constraint iteration is not reinstated.

## 2. An all-angular estimate which removes the ambiguity

In the existing C-normalized proper coordinate Y, define

    E(n)=n(n+q)/C, d_n=2(q+2n),
    A(t)=sum_(n>=1) d_n E(n) exp[-t E(n)].

Put f(E)=E exp(-tE) and g(x)=E'(x) f(E(x)), so A=2C sum g(n)
and integral_0^infinity g(x)dx=1/t^2. The trapezoid identity on each
unit cell, followed by absolute integration of its Peano kernel, gives

    |sum g(n)-integral g| <= (1/8) integral |g''(x)| dx,

because g(0)=g(infinity)=0. This is a bound on the ENTIRE infinite sum.
Since E''=2/C,

    g''=3E'E'' f'+(E')^3 f'', (E')^2=(q^2+4CE)/C^2,
    integral |f'|dE<=2/t, integral |f''|dE<=3,
    integral E|f''|dE<=4/t.

It follows, for every t>0, that

    |A(t)-2C/t^2| <= 7/t+3q^2/(4C).                 (1)

The weaker previous rectangle bound had a t^-3/2 term. That loss produced
an O(log Y) allowance after the cubic subtraction; (1) removes that loss.
No finite physical mode selection or smooth boundary replacement occurs.

## 3. The actual paired collar: bounded remainder after the cubic pole

The exact paired spin-trace/even-double identity and physical-half rank
are those of stationary-massive-boundary-layer.md. They apply to the
closed paired square; the double retains its |Y|^3 jet and is not smooth.
For the actual analytic one-sided seed Q, write the time radius as

    a(Y)=1+a2 Y^2+c Y^3+a4 Y^4+O(Y^5),
    c=Q3 sqrt(beta)/2, Q3=-4beta-2/C-tau.

Its even realization and a smooth auxiliary continuation have the same
quartic jet. Their first difference is c(|Y|^3-Y^3). Use the exact curved
Brownian-bridge formula in current-curved-paired-heat-bridge.md, including
the observation density 1/a. Taylor expansion by parabolic weight shows:
the sole difference of weight three is the already computed cubic K3;
every weight-four term is equal on the two realizations. Those terms are
the single quartic insertion and the quadratic/quadratic products, also
including the connection exponential. The first remaining difference
has weight at least five. The density correction to K3 has weight five.

Taylor's integral remainder uses bounded one-sided fifth derivatives on
a fixed collar. On paths crossing the seam, Gaussian bridge moments give

    |partial_Y^j(k_actual-k_smooth-K3)|
      <= C_j t^((3-j)/2) exp[-Y^2/(c0 t)]
                 +C_j exp[-delta^2/(c0 t)], j=0,1,2,3.        (2)

For j derivatives one can differentiate the Gaussian bridge density and
integrate the bounded coefficient remainder against its Hermite factors;
this gives t^-j/2 and a weaker Gaussian. No smooth fifth derivative of
the reflected metric is used. Products from the nonlinear denominators
have the same weight bound because a is positive. Changes outside the
collar and all cutoff commutators give the displayed exponential term.
These constants are finite for the fixed seed; their full numerical
values and a global strong source error are NOT evaluated here.

The complete antiperiodic thermal images have the Gaussian factors of
the exact Poisson formula and are exponentially small on this collar.
The true regular bolt affects the smooth positive-time/localization part;
it is not replaced by a boundary at the collar edge. The angular bound
A(t)<=24C/t^2, already established on the required short-time interval,
makes (2) integrable for j=0. For the cubic/angular discrepancy, use (1)
and the exact positive cubic profile bound

    |K3(t,Y)|<=3|c| sqrt(t)/(32 sqrt(pi)).

For a short-time cutoff t0 this part is explicitly bounded by

    C integral_0^t0 |A-2C/t^2| |K3| dt
       <= |c|[84C sqrt(t0)+3q^2 t0^(3/2)]/(64 sqrt(pi)). (3)

The derivative profile is O(t^((1-j)/2)) with the same Gaussian.
Thus its first derivative and (2) cost at most 1+|log Y|; higher
derivatives j<=3 cost at most Y^(1-j)(1+|log Y|). Long time, the
smooth subtracted collar kernel and the common finite restorations are
regular. Dominated convergence therefore yields the stronger result

    Theta_mass(Y)=c C^2/(5pi Y)+R0+O(Y(1+|log Y|)),             (4)

with the corresponding differentiated remainder bounds through order
three. R0 is a finite, UNEVALUATED global constant. In particular there
is no leading log Y term in the massive trace. The claim is about this
analytic SEED, not arbitrary mass/log metrics throughout an X ball.

The common regulator is retained before these subtractions. Its completed
a2 is -integral R^2/(480pi): the seam distribution in the curvature
divergence cancels only with the actual paired boundary completion, as
proved in current-curved-paired-heat-bridge.md section 4. Its variation
has no extra cubic-mass boundary datum, since b(0)=b'(0)=0. The a0,
topological a1, LLL, product-gradient and fixed finite restorations remain
unchanged. No separate finite part is assigned to angular modes. The
auxiliary smooth continuation in (2) is a parametrix with its full defect,
not the physical operator or domain.

## 4. Execute the complete first-source coefficients

Let Q4=3beta+1/C+epsilon, x=lambda B^2=C beta and

    k=Gq^2/(12pi C B^2),
    c_nm=4G C beta Q4/(5pi B^2),
    c_log=4k-c_nm.                                      (5)

In the original normalized canonical rows, the massive pole in (4) is
+G C beta Q3/(5pi B^2 z). Direct Euler variation of the FULL retained
Weyl action gives -G C beta Q3/(5pi B^2 z); they cancel BEFORE bounds.
The same variation gives its leading logarithm c_log log z. All remaining
local sources and completed restorations have no leading logarithm.
Equation (4) proves that the canonical massive source does not add one.

For an explicit check, with phi=-log z, L=-Q'', M=(Qphi')', V=Qphi'^2,
write A=phi L^2/60+LM/30-M^2/10+LV/30-7MV/90-7V^2/180. Its shape row is

    G/(2pi C B^2){-2C^2 A
       -2C^2 Delta[(phi L+M+V)/30]+(C/9)M
                              +2(1/15-q^2/3)phi}.              (6)

The coefficient of log z at z=0 is
[C^2(8/C^2-96beta Q4)/60-2/15+2q^2/3]G/(2pi C B^2),
which is (5). The pure 2/15 contacts cancel exactly. In the scalar row
the vacuum relation cancels the z^-4 coefficient. Weyl curvature and
trace-free Ricci squared start at z^6 and z^8; the only constant is the
magnetic anomaly. Hence the COMPLETE prescribed seed source satisfies

    f_nu=-2k+O(z(1+|log z|)),
    f_b=c_log log z+C_b+O(z(1+|log z|)).              (7)

C_b remains unknown. The fixed-temperature reconstruction and Ward
condition from the previous reference-source theorem are retained.

The reference I0's leading Euler matrix on (nu,b)=z^s(A,B) is

    beta [[24+18s-6s^2, 24+2s-2s^2],
          [12+10s-2s^2, 12]].                      (8)

The actual first correction satisfies I0''h1=-f, with its original
constraint and boundary data. Its cubic mass coefficient remains global
and unknown, but (7)-(8) determine its quartic logarithms uniquely:

    nu1=a3 z^3+a_L z^4 log z+O(z^4),
    b1=-2a3 z^3+b_L z^4 log z+O(z^4),
    a_L=(6(-2k)+7c_log)/(40beta),
    b_L=-(2(-2k)+3c_log)/(8beta).                   (9)

The O(z^4) includes an unknown constant coefficient. After subtracting
it, (4) and the reference Green equations give O(z^5 log z), with the
weighted derivatives required here. These are COMPLETE h1 coefficients,
not coefficients obtained by replacing its quantum source with J_M.
The magnetic-only values are a_L,M=2k/(5beta), b_L,M=-k/beta;
their differences give the executed outer logarithms of h_rest.

## 5. Metric variation: why the strong source needs a squared logarithm

Retain the whole actual paired kernel. Its universal leading high-angular
collar response follows directly by linearizing the bridge formula:

    delta k(t;Y)=1/(2pi t) integral_0^1
       [exp(t u(1-u)partial_Y^2)eta(Y)-eta(Y)]du.       (10)

Here eta is the time-radius variation in proper distance. For the quartic
logarithm it equals b to leading order; the change of proper coordinate
adds terms of order Y^6 log Y. Combining (10) with 2C/t^2 and the SAME
local subtractions gives the nonlocal principal part

    delta Theta=-C^2/(30pi) partial_Y^4 log|D_Y| eta
                         +a finite local fourth derivative.   (11)

The sign and factor independently reproduce the cubic pole:
partial_Y^4|Y|^3=12delta_0 and log|D|delta_0=-1/(2|Y|) away from zero.
The angular discrepancy (1), curved coefficients, cutoff terms and thermal
images have lower weight and cannot change the leading squared logarithm.
Their source remainder is at most O(1+|log z|) at the order at issue.

For a cutoff equal to one near zero, the Fourier identities (away from
zero frequency) are F[log|Y|]=-pi/|p| and
F[(log|Y|)^2]=2pi(log|p|+EulerGamma)/|p|. Consequently

    partial_Y^4 log|D_Y|(Y^4 log|Y|)
                       =-12(log|Y|)^2+O(1+|log|Y||).           (12)

The cutoff changes only lower terms. The multiplier convention agrees
with [Laptev--Weth, section 1](https://www.ma.imperial.ac.uk/~alaptev/Papers/laptev-weth-2020-revised.pdf),
where log|D| is one half of log(-Delta). Their Dirichlet realization is
NOT imported as our paired domain. Equations (10)-(12) are the local
calculation on the actual paired kernel, not a replacement by a flat
operator on the full geometry.

Put K=G C beta^2/(15pi B^2). In normalized shape units, the massive
contribution to b_L z^4 log z is +12K b_L(log z)^2. The Weyl contact
is equally essential: its principal row is
K partial_z^2[-log z partial_z^2 b], giving
-24K b_L(log z)^2. Therefore the COMPLETE actual defect
r_new=J''h1/F_lambda has

    r_new,b = A_log2 (log z)^2+O(1+|log z|),
    A_log2=-12K b_L/F_lambda.                       (13)

All finite fourth-order contacts contribute only a single logarithm.
The remaining local, scalar and lower-order geometric terms do likewise.
The unknown cubic mass and the O(z^4) coefficient cannot cancel (13).
The finite-scale logarithm changes only its lower terms. The weak pairing
is finite: all admitted tests vanish at least cubically. Thus (13) is
consistent with the C2 action and the already passed projection screen.

## 6. Executed decision test and limits

The first frozen execution stopped at a structural equality comparison
between two differently factored symbolic matrices. Preserve that run.
The replay simplifies their difference entry by entry before testing the
identity; neither the matrix nor the physical targets is changed.

Freeze numerical tests b_L<0 and A_log2>0 using all actual parameter
intervals, and quantify the magnetic and nonmagnetic coefficients. A
nonzero (13) proves that the canonical shape residual is NOT bounded by
any finite constant times 1+|log z| on the entire actual exterior. It
requires at least a squared-logarithm source weight. This is an actual
obstruction to that particular proposed strong-source estimate, not an
excessive upper bound, and not an obstruction to a full-H solution.

Likewise I0'' maps the retained single-logarithm metric phase to at most
a single-logarithm canonical source. A second reference iteration on
r_new cannot return to that phase. The full quantum Hessian has the
additional term (13), so its inverse may still map this source to the
admitted phase; no impossibility or numerical bound for that inverse is
inferred. Merely enlarging a norm or changing normalization proves nothing.

Decision: retain the construction, but make the boundary part of the
full coupled inverse the next calculation. It must keep the paired
logarithmic kernel AND Weyl contact together and establish a quantitative
return from the forced squared-log source to the admitted mass/log metric
phase, with a remainder. This is a specified sharper operator method.
The calculated boundary asymptotic is executed; that inverse and the
global remaining-source reconstruction are still unexecuted. Constants
in (2), R0/C_b, the global strong error and the actual coupled correction
are unknown. The original excessive projection allowance has not become
a physical obstruction or a complete-source value.

Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1 remain unclosed;
Y0<=5e-9, Z0<=.25, Z2<=5e7 are still sufficient unproved allocations.
The closed audit, quantitative admission margins, exact compatibility,
physical matching/causal realization, semiclassical validity/noise and
all-sector stability retain their unresolved status. No controlled
solution or wormhole nonexistence result is established.

Remaining to-do list: complete coupled outer inverse with the required
logarithmic source weight, global prescribed-source reconstruction and
controlled correction; nonlinear admission and exact compatibility;
physical matching and causal realization; semiclassical validity/noise
and all-sector stability; charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
