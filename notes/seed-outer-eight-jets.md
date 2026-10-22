# Quantitative outer source jets through order eight

2026-10-21, E003/S167. Ordinary argument before supporting verification.
Keep the centered seed, fixed physical scale, HH state, paired domain,
full angular tower and common four-dimensional subtraction. Put

    D=Y partial_Y, L(Y)=1+log(1/Y),
    S(Y)=c C^2/(5pi Y), c=-a sqrt(x)/C^(3/2),
    R(Y)=Theta_seed(Y)-S(Y).

For the original proper coordinate Y, the result is

    |D^j R|/(j!100^j q^2 L(Y))<1e44,
    |Y^j partial_Y^j R|/(j!100^j q^2 L(Y))<1e44,
    0<=j<=8, 0<Y<=.1.                            (1)

The much smaller original outer C3 allowance 1e20 remains available;
it is not replaced by this high-order reserve. The singular S is
retained in the physical stress reconstruction with its correlated
Weyl term. Removing it in (1) is not a new finite subtraction.

This complements [the middle eighth-order estimate](seed-middle-eight-jets.md).
The cold weighted counterpart, ten metric derivatives and the actual
finite-metric H6 Schur source are not supplied. No cancellation,
nonlinear geometry or physical stability is established.

## 1. Smooth auxiliary caps and a true eight-jet parameter family

Use the two cap constructions of
[the outer source proof](stationary-centered-conormal-source.md) and
[the fixed-cap heat calculation](stationary-reference-nonbolt-stress.md).
The even cap reflects at Y=0; the smooth cap continues the analytic
seed negatively and flattens on [-1,-1/2]. In particular its taper
starts at -1/2, not at -1. Both positive caps agree with the actual
seed through a seam at least 24 beyond the observation collar.

For the new high-derivative auxiliary family, taper w=f'/f with
the degree-21 C10 step from the middle estimate, of width 1/2.
The real taper has |eta'|<7.6 and |eta''|<160 in Y units.
Use the retained seed guards (|w|,|w'|,|w''|)<(1.01,2.021,13)
on the positive side and (1,3,9) on the negative side. Direct
Leibniz and f'''/f=w''+3ww'+w^3 give the global real allowances

    |f'|/f<2, |f''|/f<20, |f'''|/f<400.           (2)

The cap f is C11. Taper values remain between zero and one, so
the old ratio and circumference proof is unchanged:

    R_even<1e14, R_smooth<4e14, d_min>1e17.        (3)

These are auxiliary metrics, not other physical HH exteriors.
For the positive seam, calB/C>=(1-u)^2/8 gives its proper coordinate
below sqrt(8) log(1/zc)+1/2<94, zc=exp(-24)/1001. Thus every
nonconstant cap coefficient lies in -1<=Y<=100.

At each real taper point use its analytic polynomial branch. On a
complex spatial disk of radius .01 the step argument changes by
at most .02. With c21=21!/(10!)^2, the same beta estimate gives
|eta'|<=c21(1/4+2(.02)^2)^10<4 and |eta|<1.08. The analytic seed
has |w|<4 on its .05 strip. Hence the local capped branch has
|w|<4.32 and |f/f_real|<exp(.0432)<2. Cauchy gives

    |f^(m)(Y)|/f(Y)<=2 m!100^m, 0<=m<=11.        (4)

At taper endpoints these branches match the true jets through
eleven. Outside the taper use the analytic seed or constant branch.
There is no global spatial analyticity claim for a cutoff.

Choose the complex scale disk |lambda|<=rho=1e-5. On the negative
cap define the finite, parameter-holomorphic family

    P_lambda f=sum_(j=0)^8 lambda^j D^j f/j!.     (5)

Blend it to f(exp(lambda)Y) on [-1/2,-2/5], using the same C10
step at width .1. On the positive side use f(exp(lambda)Y) to
Y=4 and blend to the fixed positive cap by Y=4.5. Reflection
defines the even family's negative side. At lambda=0 the smooth
and even caps coincide on Y>=0 with their respective actual jets.

The finite negative family has the TRUE dilation jets through
eight, including spatial derivatives through three. Only f through
eleven is required. On its nonconstant support |Y|<=1, differentiate
D^j f=sum_k S(j,k)Y^k f^(k), and use (2) for j=0 and (4) for
the remaining terms. Exact rational sums give

    |P_lambda/f-1|<.003,
    |P_lambda'|/f<2.5, |P_lambda''|/f<150,
    |P_lambda'''|/f<50000.                        (6)

On the negative blend the analytic and finite families agree to
order eight in lambda. On |lambda|<=.01 their analytic spatial
jets through three have relative bound 1e6. The Taylor tail on
the smaller disk is below 1e6(.001)^9/(1-.001). All blend products
through three cost below 1e8, so their discrepancy is below 1e-6.

For the positive blend, |exp(lambda)Y-Y|<.00005 and the analytic
seed jets have relative bounds 2m!40^m on the .025 strip. Expand
in this displacement and include exp(a lambda) for a spatial
derivatives. With (2) as the real baseline, the difference bounds
through three are less than (.01, .4, 50, 7000) times f. The
width-.5 blend has derivatives below (7.6,160,45000); its product
terms still give the same global guards used by the outer heat proof:

    |f_lambda'/f_lambda|<5,
    |f_lambda''/f_lambda|<1e3,
    |f_lambda'''/f_lambda|<1e5, |arg f_lambda|<.05. (7)

The strict ratio and circumference margins in (3) absorb the small
parameter changes. On the analytic observation/anchor patches the
old spatial-strip test is strengthened by the smaller scale disk:
2.4(exp(rho)-1)+exp(rho)/40<.05.

## 2. Scale the angular weight and finite parts together

Write Atilde(tau)=sum d_n M_n^2 exp(-tau M_n^2). The exact cap heat
identity, including the covariant density and thermal circle, is

    k_f(exp(2lambda)tau,exp(lambda)Y;beta)
       =exp(-2lambda) k_(f_dil)(tau,Y;beta exp(-lambda)),
    f_dil(Y)=f(exp(lambda)Y).                     (8)

For the subtracted heat kernel Ksub the same identity holds:
k1_dil=exp(2lambda)k1(exp(lambda)Y) and
k2_dil=exp(4lambda)k2(exp(lambda)Y). Changing the proper-time
integration variable therefore gives the convergent expression

    integral_0^infinity Atilde(exp(2lambda)tau)
                         Ksub_(f_dil)(tau,Y;beta exp(-lambda))dtau.
                                                        (9)

The measure factor and heat factor cancel. Do NOT leave Atilde
unscaled in (9). Keep the common a0 restoration constant and
multiply the local k2 restoration by exp(-4lambda); its angular
finite part and fixed log(mu0 B) remain unchanged. Differentiating
(9) once reproduces the old 2 tau Atilde' Mellin term and the
restoration's -4k2 term. Thus (9) includes all their higher analogues.

The cubic boundary coefficient is c_lambda=exp(3lambda)c. Its
continuum angular profile scales by exp(-4lambda). The complete
principal term is consequently exp(-lambda)S(Y)=S(exp(lambda)Y),
as required. This cancellation is imposed before bounding the
remainder. No boundary coefficient or physical scale is reset.

Use the family (5) in place of f_dil in (9). Their negative jets
agree through eight; their positive defect is kept in section 4.
This supplies a holomorphic auxiliary finite remainder F(lambda,Y).
It is not an assertion that the full physical exterior has a global
complex dilation family.

For |lambda|<=rho, Re exp(2lambda)>.9999. Absolute values of the
scaled angular integrals are bounded by replacing tau with a tau,
a=.9999, in their exponential factors. Even after the largest retained
heat growth exp(100tau), their effective rescaling is above .9998
because M^2>=324. The short, long and winding budgets cost at most
.9998^-6<1.002.
For example the smooth M^-4 integral costs a^-3 and the angular
continuum/error heat weights cost a^-2 or a^-3/2. The reflected
principal profile is treated exactly as above, before the remainder.
The finite local and winding estimates also keep their strict margins.

With (7), the old complete auxiliary table is therefore still below
4e19 q^2 L(Y), allowing a factor two over its prior 2e19 budget.
This includes the common finite parts and thermal winding terms,
with the scaled circle in (8); it does not set windings to zero.
Cauchy gives

    |partial_lambda^j F(0,Y)|/(j!100^j q^2 L)
       <=A_j=4e19*1000^j, 0<=j<=8.               (10)

The distant actual/cap and positive-dilation defects are added next.

## 3. Actual/cap defect: one-sided observation jets through eight

Use the original two-commutator identity for the complementary
positive square. The actual inverse norm is M^-2. The cap and actual
operators agree through the inner collar, at distance D>=24 from
the observation. No actual off-diagonal kernel bound is introduced.

The local cap-row bootstrap extends to total observation order eight.
Here is a conservative budget for that extension. Tangential
derivatives commute with the constant paired boundary projectors.
On nested half-balls with widths 1/[20M(j+2)], j<=10, the cutoffs
cost at most 480M. The coefficient is within the old .07 energy
margin and its analytic half-side jets are bounded by 2r!40^r.
After normalizing the tangential coordinate at the observation,
the first-order equation is

    partial_Y v=-gamma_Y gamma_t f partial_t v-M gamma_Y v.

Tangential energy, this equation and the differentiated coefficient
products have a reserve j!(1e5 M)^j through order eight. For instance
the normalized coefficient convolution is bounded by
2 sum_(r=1)^10 10!/(10-r)! (40/(1e5*18))^r<.001;
the cutoff term is below .005. The point evaluation step uses two
additional tangential energy orders, already included in j<=10.
All normal derivatives are then obtained on the smooth half-side.
No eighth derivative of the C2 even reflection is assumed.

The undifferentiated physical point row is below 3e4 exp(-MD/2).
Incoming graph and fixed commutator factors fit below 1e8 M times
this exponential: the cap curvature is bounded by the real w,w'
budgets above, weighted inverse is below 2/M and the graph constant
ten is sufficient at M>=18. The much larger prior allowance
R*=1e12 C0, C0<=1e23, therefore covers all these endpoint factors.
With the original inverse powers of M retained, the two-commutator
identity and trace admit the permissive per-mode bound

    |partial_Y^m G_mode|<=1e74(m+1)m!(1e5)^m
                         M^(m+2) exp(-9.6M), m<=8.            (11)

It allows two excess mass powers and bounds both paired traces.
It follows by the m+1 Leibniz distributions; rank, commutators,
graph and density factors fit the additional 1e4 over (R*)^2.
No derivative is placed on a distant cutoff by the observation jet.

For 2<=p<=10, the actual benchmark M_n>=18sqrt(n), d_n<=6qn and
monotonicity of M^p exp(-9.6M) on M>=18 give the full moment

    sum d_n M_n^p exp(-9.6M_n)<=q H_p,
    H_p=6*18^p*1e-74*[1+2/(172.8-p-3)].           (12)

The integral comparison uses u^(p+3)<=exp[(p+3)(u-1)], u>=1.
The retained rational exponential certificate gives exp(-172.8)<1e-74.
Thus (12) is an infinite-tail bound, not an angular cutoff.

Since D^j=sum_(m=0)^j S(j,m)Y^m partial_Y^m, define

    G_j=(.00352)*1e74/(j!100^j)
          *sum_(m=0)^j S(j,m)(m+1)m!(1e4)^m H_(m+2).         (13)

For j=0 take S(0,0)=1. Equation (13) uses Y<=.1 and C/q<.00352.
It bounds the normalized actual/cap correction by G_j q^2 L,
since L>=1. The separation 24 is used only for this defect.

## 4. Positive dilation defect: retain its distance 3.9

The analytic family agrees with the true dilation only through Y=4.
Its defect for an observation Y<=.1 has distance at least 3.9.
It is not assigned the distance 24 from section 3.

The true cap dilation has real derivatives through eight. Local
analytic branches bound these jets on |lambda|<=rho because the
nonconstant cap is inside |Y|<=100 and its displacement is below
.0011<.01. The fixed blend is parameter analytic. After transport
of the scaled circle to a fixed circle, include the corresponding
coefficient factor and the scaled mass. The resulting first-order
vertices have graph bounds 2j!100000^j. The weighted graph inverse
costs at most ten. Every term in the difference of the ordered
inverse jets has a vertex supported in Y>=4.

For each j, telescoping the ordered products costs at most
4j!100000^j 21^(j-2)(20j+1), below j!1e7^j for 1<=j<=8.
The two endpoint rows, incoming graph factor, trace, paired rank
and explicit scalar prefactors fit the conservative mode bound

    1e16 j!1e7^j M^2 exp(-3.9M).

The original full angular bound sum d_n M_n^2 exp(-3.9M_n)
<2.1e-27q therefore gives normalized reserves

    Q_0=0, Q_j=1e16*(1e5)^j*(.00352)*(2.1e-27), j>=1.       (14)

The zero entry follows because both cap families agree at lambda=0.
The scaled angular masses and circle are included here; they are not
discarded when differentiating the defect. A common variation outside
the defect does not remove its required remote vertex.

## 5. Complete eighth-order aggregate and ordinary derivatives

The exact real scaling identity, negative jet agreement and the two
separate defect comparisons prove the Euler-derivative bounds

    E_j=A_j+G_j+Q_j, j=0,...,8.                  (15)

All three terms use the same physical subtraction and principal S.
For ordinary derivatives use the signed Stirling coefficients of
the first kind, s(j,k), in

    Y^j partial_Y^j=D(D-1)...(D-j+1),
    O_j=sum_(k=0)^j |s(j,k)| k!100^k E_k/(j!100^j).        (16)

Exact rational aggregation gives max(E_j,O_j)<1e44. This proves
both parts of (1), including the full ordinary eighth derivative.
It is stronger than merely asserting finite conormal regularity,
and remains an absolute estimate rather than an actual source value.

Freeze this proof, checker and review before execution. Verify the
new taper and finite-jet bounds, the exact scaled subtraction and
boundary profile, all inverse insertion and angular-tail factors,
and both derivative aggregations. Independent exponential heat
fixtures must include the angular scaling and finite k2 factor.
Continuum heat and row steps remain ordinary arguments, not a
machine proof. No Lean run is selected.

Remaining to-do list: cold weighted seed jets through eight; the
quantitative ten-derivative metric and finite-metric quantum/Schur
source bounds; positive transfer constants or a sharper actual
transfer; cancellation, signed compatibility, nonlinear existence;
physical validity, state/boundary, stability, matching, charge/return,
opening, transport, preparation, separate creation and audit.
