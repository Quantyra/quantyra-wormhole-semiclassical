# Complete mixed-frequency loop in the retained zero-winding comparison

2026-10-24, E003/S169. Ordinary argument before frozen acceptance.
Full stable-wormhole objective ACTIVE/INCOMPLETE; feasibility audit CLOSED.

This executes the mixed-frequency calculation needed after the coefficient-
norm obstruction. The stationary flat comparison's massive third loop,
with its pair contacts, both triangle orientations and fixed fourth
restoration, has a finite all-frequency bilinear Sobolev bound. Its
apparent soft denominators cancel exactly. This is a NEW nonlinear
component, not optimization of the already controlled linear angular term.

It is not the complete curved third response, physical weighted graph
return, or signed bolt projection. In particular the original sharp
near-input split is not an H8 input to a whole-line estimate. Its local
annulus returns and separated mixed tails still need bounds. The actual
paired seam, bolt, state, lower contacts and curved difference stay in
the exact allocation; none is replaced by a smooth physical surrogate.

## 1. Scope and exact common-regulator inventory

Use the SAME translation-invariant zero-winding comparison whose Hessian
is K in outer-flat-tower-normal-remainder.md. Its stationary momenta are
p1+p2+p3=0, masses M_n^2=n(n+q)/C, d_n=2(q+2n), n>=1. The LLL and
explicit local action remain in their existing separate allocations.
Write a(alpha)=alpha(1-alpha), and

 R3(t)=log(1+t)-t+t^2/2,
 R4(t)=R3(t)-t^3/3,
 A(p)=sum_n d_n M_n^2/(2pi) int_0^1 R3(a p^2/M_n^2) d alpha,
 B(p)=sum_n d_n (2M_n^4/pi) int_0^1 R4(a p^2/M_n^2) d alpha. (1)

A and B are even smooth real functions. Both sums converge locally
with all needed frequency derivatives: their large-M remainders are
O(M^-4) at fixed p, and the entire angular count is N(M)<=4CM^2.
These are complete sums, not finite angular approximations.

For one massive Dirac field let S_M=(D+M)^-1 in the flat comparison.
Differentiating -Tr log(D+M exp b) gives the COMPLETE third variation

 -M Tr(S_M u v w)
 +M^2 Tr(S_M uv S_M w+S_M uw S_M v+S_M vw S_M u)
 -M^3 Tr(S_M u S_M v S_M w+S_M u S_M w S_M v).       (2)

Equation (2) is used with one common regulator, consistently with the
earlier positive-time Duhamel inventory. The first term and all three
pair terms are present. At fixed M the momentum-independent divergence
is local. No individual angular vertex is assigned its own finite part.

Let T_ge4 denote the third comparison kernel after subtracting its
degree-zero and degree-two momentum Taylor contacts, RETAINING the
prescribed fourth contact. The exact full allocation is

 T_actual = T_ge4 + T_low + T_curved/images/return + T_conversion/local,

with the original thermal and relative split understood. T_low contains
the degree-zero/two tadpole and local completions removed in defining
T_ge4. They are not set to zero or certified by the present calculation.
The known seed bound on a vacuum coefficient is not a computed uniform
third response. This note only bounds the defined T_ge4 comparison block.

## 2. Derive the loop before estimating it

For nonzero p_i the scalar bubble and triangle identities give

 L_M(p)=int_0^1 log(1+a p^2/M^2) d alpha,

 T_M(p1,p2,p3)-T_M(0,0,0)
   =M^2/(2pi){sum_i [1+4M^2/(p_j p_k)]L_M(p_i)-2}.   (3)

Here j,k are the two indices other than i. To check the matrix and
normalization content directly, the spinor bubble trace is
2[M^2-l.(l+p)]/(D_l D_(l+p)). Its difference from p=0 integrates to
L_M(p)/(2pi). The triangle trace is
2M[M^2-l1.l2-l2.l3-l3.l1]/(D1 D2 D3). Use
l_i.l_j=(D_i+D_j-2M^2-p_ij^2)/2 to reduce its numerator.
Label vertices D2,D3,D1 by collinear shifts (0,p1,p1+p2), so the
edge opposite D_i has transfer p_i. Partial fractions give

 1/(D1 D2 D3) = -sum_i [1/(p_j p_k)]/(D_j D_k),

where each remaining bubble has transfer p_i. Finally the elementary
bubble integrals obey

 (4M^2+p^2) int_0^1 [M^2+a p^2]^-1 d alpha=2L_M(p)+4.

Insertion into BOTH triangle orientations and all pair terms gives (3).
The singular-looking formula is only an intermediate identity; its
continuous extension is established below. The common degree-two term
of (3) cancels. Its degree-four term is
sum_(i<j) p_i^2 p_j^2/(60pi M^2).

Subtract this common quartic Taylor term before summing the angular
tower. The exact remainder is

 T_sub = sum_i A(p_i)+sum_i B(p_i)/(p_j p_k).         (4)

The fourth polynomial has not been discarded: section 3 supplies its
fixed common-regulator restoration. The Taylor removal in (4) can be
identified at the original common positive cutoff. Six external-momentum
derivatives of the combined Duhamel words give an integrable massive
bound O(M^-4), uniformly as the cutoff decreases, on compact momentum
sets. This follows by combining denominators before Taylor's integral
remainder; their positive massive denominators bound the remainder by
a polynomial in the fixed momentum radius times M^-4. Angular counting
makes that envelope summable. It proves common finite-part identification
of (4) as a distribution and on smooth test inputs. It does NOT prove
operator-norm convergence of a raw cutoff after a fixed unregulated
normal subtraction. The previous cutoff obstruction remains in force.

## 3. Restore the finite fourth contact in the same scheme

The retained quadratic comparison has

 K(p)=A(p)+c4 p^4,
 c4=C kappa4/(30pi),
 kappa4=H_q/2-gamma_E-log(mu0 B).                     (5)

The cubic fourth contact is not obtained by an arbitrary nonlinear
extension of K. A constant mass shift b=c scales the regulated flat
Hessian as exp(2c) H_(exp(2c)epsilon,0)(exp(-c)p).
Its a2 Hessian contact is -C exp(-2c)p^4/(30pi), so removal of the
SAME cutoff adds +c C exp(-2c)p^4/(30pi) to the finite part. Hence the
soft third identity is

 T_ge4(0,p,-p)=2K(p)-pK'(p)+C p^4/(30pi).            (6)

The last term is required by the regulator scale. Omitting it would
misidentify the finite fourth contact despite having the right leading
logarithm. At stationary cubic fourth order the only nonzero symmetric
integrated local density, modulo total derivatives, is b(b'')^2.
For example (b')^2 b'' is a total derivative. Thus the soft identity
fixes the entire stationary quartic polynomial:

 T_ge4=T_sub+d4 sum_(i<j) p_i^2 p_j^2,
 d4=C(1-2kappa4)/(30pi).                            (7)

Boundary representatives are not inferred from this flat polynomial.
They remain in the actual paired/local conversion allocation. The
compact comparison calculation has no new physical boundary condition.
The additive common-regulator term in (6)-(7) is retained, not fitted.

## 4. Remove every soft denominator, including opposing frequencies

Put f(p)=pB(p), an odd smooth function. For real a,b and c=-a-b,

 [f(a)+f(b)+f(c)]/(abc)
   =1/2 int_[0,1]^3 f'''(a(t+r-1)+b(s+r-1)) dt ds dr. (8)

This identity follows by two integrations of f'', then pairing (t,s)
with (1-t,1-s) and integrating f''' across the remaining centered
interval. It therefore holds continuously when any of a,b,c is zero.
The argument of f''' in (8) is a convex combination of vertices
0,+/-a,+/-b,+/-(a+b); its absolute value is at most
P=max(|a|,|b|,|c|). No 1/a, 1/b or 1/(a+b) loss remains.

In particular (4) defines a single continuous kernel at equal,
opposing, soft and zero frequencies. This is an all-frequency argument,
not an inference from evaluating a few modes. The useful exact identity

 (pB(p))'=-p^3 A'(p)                               (9)

also follows by integration by parts in alpha in (1). Equivalently,
the beta moments I_n=int a^n=(n!)^2/(2n+1)! obey
I_(n+1)/I_n=(n+1)/[2(2n+3)]. Power series prove (9) near zero and the
defining real integrals extend it to every real p. At a zero frequency,
(4) is 2A(p)-pA'(p), consistent with (6)-(7).
The beta identity is [DLMF 5.12.1](https://dlmf.nist.gov/5.12.E1);
its positive-parameter hypotheses hold for every moment used here.

## 5. A whole-frequency bound with the complete angular tail

For R3, R3(t)<=t^3/[2(1+t)] for t>=0; integrate
t^3 int_0^1 u^2/(1+tu) du and use 1+tu>=u(1+t).
For f_M(p)=p(2M^4/pi)int R4(a p^2/M^2)d alpha,

 |f_M'''(p)|=(2M^4/(pi p^2)) int
       t^4(126+180t+70t^2)/(1+t)^3 d alpha
 <=252M^4/(pi p^2) int t^4/(1+t) d alpha.            (10)

Both expressions extend at p=0. Set m=18<M1. Integration by parts
against N(M)<=4CM^2 gives, for A0>=0,

 sum_n d_n/[M_n^2(M_n^2+A0)]
 <=4C{log(1+A0/m^2)/A0+1/(m^2+A0)}
 <=8C log(1+A0/m^2)/A0.                            (11)

At A0=0 read the continuous value. This counts the whole tower.
Using I2=1/30 and I3=1/140 in (10)-(11) gives

 A(p)<=C p^4/(15pi) log(1+p^2/(4m^2)),
 sup_(|r|<=P)|f'''(r)|<=72C P^4/(5pi)
                                      log(1+P^2/(4m^2)).

Equations (4),(8) imply

 |T_sub|<=37C P^4/(5pi) log(1+P^2/1296).             (12)

The actual fixed parameters give 0<kappa4<16. Since
sum_(i<j)p_i^2 p_j^2<=3P^4, (7),(12) and
log(1+P^2/1296)<=2log(2+P) prove the convenient guard

 |T_ge4(p,q,-p-q)|<=7C P^4 log(2+P).                (13)

It includes the finite fourth contact. Lower Taylor contacts are still
in T_low, not covered by (13). This is a finite useful upper bound,
not yet an actual curved source norm.

With unitary Fourier convention, let T_ge4[u,v] be the bilinear
response kernel (7). Define H8,log by multiplier
<p>^8 log(2+|p|). Young's convolution inequality now proves

 ||partial_z T_ge4[u,v]||H3
 <=1280 C {||u||H8,log ||v||H1+||v||H8,log ||u||H1}. (14)

For completeness put F(t)=(1+t^2)^4 log(2+t), t>=0. It is increasing
and convex: differentiating twice gives a lower bound
(8log 2-1)(1+t^2)^3>0. Also F(2t)<=512F(t). Since |p+q|<=P,
the output derivative and H3 weight times (13) are at most 7C F(P).
Use P<=|p|+|q|, then convexity at the midpoint of 2|p| and 2|q|:
F(P)<=256[F(|p|)+F(|q|)]. This gives the factor 1792C without
an unaccounted loss from the two input allocations.
Fourier L1<=sqrt(pi) H1 and the convolution factor
1/sqrt(2pi) give 1792/sqrt(2)<1280. This proves (14) for all inputs
in the stated Sobolev spaces, not just wave packets. Keeping the highest
logarithmic derivative in L2 is sufficient for this comparison block.
The constant-beta replacement p->sqrt(beta)p only decreases the guard.

The physical paired comparison uses the same even input returns as the
quadratic comparison. All direct/reflected contributions are supplied
by those returns. Their norm factors and the actual weights must be
included when (14) is transferred to the physical half; they are not
silently set to one here.

## 6. Executed signed diagnostics with certified entire-spectrum error

These are evaluations of the comparison kernel, NOT selected physical
states or the actual remaining bolt projection. They check (4),(7),(8)
with all angular levels retained. Define h_r=sum d_n/M_n^(2r), r>=2.
Partial fractions give the exact entire-spectrum identity

 h_r=2C^r sum_(j=2)^r a_(r,j)
                      [zeta(j)+(-1)^(j+1)zeta(j,q+1)],
 a_(r,j)=(j-1)/(r-1) (-1)^(r-j)
                binom(2r-j-2,r-j) q^(-(2r-j-1)).     (15)

For (15), differentiate [n(n+q)]^(1-r) and expand its finite principal
part at n=0. The substitution n->-n-q fixes the other pole by
antisymmetry. Summing the resulting powers gives the two zeta functions;
the sum starts at j=2 and is absolutely convergent. There is no angular
cutoff. For P<2M1 the convergent expansion is

 T_sub=(1/(2pi)) sum_(n=3)^infinity (-1)^(n+1) I_n h_(n-1)/n
   {sum_i p_i^(2n)-[2n/(2n+3)] sum_i p_i^(2n+3)/(p1p2p3)}. (16)

Newton's power-sum recurrence with e1=0 makes the braces a polynomial
in e2,e3. It is evaluated that way at zero and nearly zero frequencies.
For truncation at J, r=P^2/(4M1^2)<1, a rigorous full tail is

 5 h2 P^6/(280pi) r^(J-2)
 { (J+2)^2/(1-r)+2(J+2)r/(1-r)^2+r(1+r)/(1-r)^3 }.   (17)

To prove it use h_(n-1)<=h2/M1^(2n-6), I_n<=I3/4^(n-3), and (8)
for the polynomial p^(2n+3). After division by n the absolute brace
is at most 5(n+1)^2 P^(2n). Sum that geometric second moment exactly.

At J=12 the interval executions give
  (p1,p2,p3)=(1,1,-2): T_sub=-19517129.872887967...,
      series tail <2.640e-15;
  (1,-1,0): T_sub=-751548.4028779413..., tail <3.924e-23;
  (1,-1+2^-20,-2^-20): T_sub=-751546.2530863132...,
      tail <3.924e-23.
The archived balls also include parameter and arithmetic uncertainty.
The fourth contact is large in these unnormalized action units and is
explicitly included in each T_ge4 result. Its omission is not justified
by the smooth behavior at a zero external frequency.

## 7. Quantitative consequence, remaining domain work and decision

The original normalization 2G/B^2 turns (14)'s prefactor into

 D_flat=2560GC/B^2
       =1.500000000030537656815...e-29.              (18)

Let L_transfer denote a CERTIFIED bound for applying this tame estimate
after the two original fast returns, actual weights, local coordinate
changes, even extensions and required near/far localization, with all
separated mixed tails accounted for. It is NOT presently known.
The sharp indicator 1_(z<.02) has a jump, so its whole-line H8 norm is
infinite even for a smooth input with nonzero value at .02. Directly
plugging that restricted field into (14) is INVALID. Use smooth analysis
cutoffs strictly inside the split and retain the separated remainder,
as in the quadratic angular argument. This changes no physical domain.

If C3_add bounds the other complete graph contributions, the unchanged
remaining sufficient requirement would be

 eta_zero+1e-11[D_flat L_transfer+C3_add]<.159916395999. (19)

Even in the best case eta_zero=C3_add=0, this requires
L_transfer<1.066109306638295679911...e39. An equal seed/variation split
would leave L_transfer<5.33054653319147839956...e38 if C3_add=0.
These are thresholds for missing estimates, NOT achieved values. There
is no complete finite actual upper bound to classify as excessive.
The Fourier bound removes a potential mixed-frequency singularity and
is consistent with the full graph's derivative count. It does not
prove that the weighted/curved transfer meets either threshold.

The first unsupported continuum step is now explicit: the actual
paired, weighted bilinear return, its separated mixed tails and curved
coefficient defect, with the complete low contacts and a matched strong
regulator limit uniform on admitted background GRAPH/mass-log fields.
A-only/C5 control remains insufficient. The closed first-order domain
and weak action theorem do not supply this strong result. No smooth
replacement domain, arbitrary stress or finite mode model is adopted.

The relative and thermal bounds stay parked. The contact targets
M_A<=3.18e27, M_V<=2.66e30, remaining global A<=.00094319 and original
X/graph comparisons are still unproved. At rho=1e-8 the exact nonlinear
requirements remain Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1,
with unproved allocations Y0<=5e-9, Z0<=.25, Z2<=5e7. No complete
|P_hat_rest|+E_rest<1.287e-40 screen, actual correction or controlled
solution is computed. That screen is sufficient, not necessary for
existence. Compatibility, matching, validity/noise and stability stay open.

RESEARCH DECISION: retain the construction and the direct full-graph
method. Use the contact-complete divided difference (4),(7),(8), rather
than separately bounding triangle denominators. The executed whole-
frequency estimate is finite and has the required logarithmic order;
the actual weighted and curved estimate remains undecided.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): execute the physical
weighted bilinear return of this third kernel using annular analysis
cutoffs, both original normal inverses and all separated mixed tails,
together with the retained curved-defect decomposition. Prove its
matched strong regulator limit on the stated background graph family,
retain T_low and finite boundary completions, and test (19). Do not
replace this transfer by an identity of Sobolev norms or a norm rescaling.

Remaining to-do list: that actual graph transfer and curved/low-contact
response, global admission/correction and original X comparison,
nonlinear closure, exact compatibility, physical matching/causal
realization, semiclassical validity/noise, all-sector stability,
charge/return, opening, transport, preparation, separate creation and
the full stable-wormhole objective.
