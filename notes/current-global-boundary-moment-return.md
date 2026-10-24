# The original global boundary response: extraction, conditioning and targets

2026-10-24, E003/S169. Ordinary argument before frozen execution.
The full stable-wormhole objective is active and incomplete. The audit stays
closed. This defines the missing global response on the actual fixed-seed
solution and quantifies what must be bounded. Its unknown global terms are
not evaluated or replaced by a model. No nonlinear solution is established.

## 1. An earlier requirement than assigning a matching coefficient

The preceding normal test did not prove an expansion separating z^(4+e_*)
from the actual remainder. A weight-four regularity conclusion alone does
not justify reading off that coefficient. First cross the one slow root in
the actual fixed-seed equation, retaining its original row and all moments.

Use the actual paired Hessian realization, canonical rows, fixed normal
family and coefficients of current-coupled-outer-return-test.md. Let u be
the existing weak solution H u=-f, with its regular bolt and mass tail.
This is a qualitative existence input for a LINEAR response, not the desired
controlled nonlinear solution. Set A=z^4 H, F=z^4 f, A=A0+R, where

    A0(s)=D0 C(3-s)^T diag(J(s),N(s)) C(s),
    D0=-120pi^2 K, B(s)=4D0(s-4)N(s).

All symbols, contacts and the physical measure retain their previous values.
Let m e_m be its actual correlated mass cubic with a fixed smooth cutoff,
w=u-m e_m, and v=chi w, for a smooth cutoff chi equal to one near zero.
The exact identity on the whole normal half-line is

    A0 v = chi[-F-A(m e_m)-R w]+[A0,chi]w = V.       (1)

The physical part is extended by the stated localization only; the
commutator includes its full separated tail. It is not a new boundary
condition at the cutoff. Apply T(V)=(D-4)V_b-(D-6)V_nu/3 and define G=T(V).
Then

    B(D) v_b=G,
    G=g0 z^4+O_b(z^5 log^N z),
    g0=-(c_log+2f_nu,0/3)=B(4)B0.                  (2)

The remote z^4 kernel moment cancels in T near zero, as already proved on
the actual paired realization. It survives in the localized original (1)
through its transition terms. None of A(m e_m), Rw or [A0,chi]w is omitted.
All angular and thermal modes, local contacts and bolt effects remain in
these terms. In particular G is not just the known seed source.

The retained separated kernel has degree -5, its difference degree -3,
and its first moment Taylor remainder gains a power after T. The actual
one-sided analytic seed and completed boundary contacts therefore give (2).
At the opposite end of the auxiliary normal half-line the commutator tail
is O(z^-1 log^N z). These estimates are for the actual localization, not
a smooth reflected physical geometry. Their numerical constants are unknown.

## 2. Finiteness of the needed derivatives at this seed only

For the already proved source expansion, the exact paired Brownian bridge
argument can be differentiated once more than its previously displayed
three derivatives. Four observation derivatives cost t^-2 and a Gaussian
Hermite factor. The weight-five coefficient remainder then costs t^-1/2;
the all-angular bound and the Gaussian cutoff give O(Y^-3(1+|log Y|))
after integration. The cubic/angular discrepancy has the same derivative
order. Thus the fourth conormal derivative of the pole-subtracted source
is still O(Y(1+|log Y|)). Positive-time and separated terms are regular.
The fixed analytic one-sided seed supplies the required bounded jets; no
smooth fifth derivative of its reflected cubic is invoked.

Commuting four conormal derivatives in the retained fixed-seed elliptic
graph argument gives eight disk derivatives (with the shape logarithmic
multiplier). The differentiated, moment-subtracted R and commutator in
(1) consequently have three conormal source derivatives. This extends the
finite fixed-seed bootstrap by one order. It is NOT a numerical estimate
for those derivatives, nor an assertion on arbitrary nonlinear C5 families.
The closed critical inverse supplies u; this argument cannot construct u
by solving only the differentiated equation.

## 3. Cross exactly the slow root, not a sampled strip

Put sigma=4+delta, delta=1/200, e=s-4 and L=A_T/d0. In
4<=Re s<=4.01, H(s)=1/e+h(e). For |Im s|<=2, the digamma series gives
|psi(s-1)|<2, while the pole-subtracted cotangent series has modulus <2;
hence |h|<18. Also |q(s)+2|<21. Outside |e|<.01 this implies
|(q+2)H|<3000, far below L>3e34.

On |e|<=.01 the same absolutely convergent series give |h|<17,
|h'|<4 and |h''|<4. The pole-cleared equation is

    R0(e)=(6+5e+e^2)(1+e h(e))-L e.

It satisfies |R0-(6-Le)|<108|e|, |R0'+L|<120 and |R0''|<300.
Rouche on |e|=.01 gives exactly one zero inside that disk. The previously
enclosed real e_* is that zero.

For t=|Im s|>=2, the decreasing-series comparison for Im psi(x+it),
3<=x<=3.01, gives 0<pi/2-Im psi<arctan(x/t). Also
|cot(pi s/2)+i|<=2/(exp(pi t)-1) for t>0. Hence |t Im H|<5.
Re H>5, and Re H<15 for 2<=t<=sqrt(7); these follow from
Re psi(x+it)<=psi(x)+t^2 sum_(n>=0)(n+x)^-3 and the same cotangent bound.
It follows that Re[(q+2)H]<100 for all t>=2. No zero lies there.
Complex conjugation treats negative t. Thus the strip from weight four
to sigma contains precisely the simple root 4+e_*, with no other zeros.
This is an entire-strip argument for the retained family; no finite
physical spectrum or contour samples substitute for it.

The identities used here are the digamma and cotangent series in
[DLMF 5.7.6](https://dlmf.nist.gov/5.7.E6) and
[DLMF 4.22.3](https://dlmf.nist.gov/4.22.E3). They provide no paired-domain
or determinant theorem; those hypotheses are addressed in (1)-(2).

## 4. An exact original-equation formula without catastrophic subtraction

Define Mellin transforms using integral z^-s (...) dz/z. Its continuation
through weight four is

    Ghat(s)=-g0/(s-4)+M_G(s),
    M_G(s)=integral_0^1 z^-s[G-g0 z^4] dz/z
                    +integral_1^infinity z^-s G dz/z.             (3)

(2) and the separated tail make M_G holomorphic near 4+e_* and on the
needed strip. This fixes the finite part; it is not a change of the physical
regulator or a freedom in the action. On the line Re s=sigma, Ghat is the
ordinary transform of G-g0 z^4 on the WHOLE half-line. The artificial
subtraction is only an exact Mellin identity and its entire tail is retained.

Move the Mellin inversion contour, initially with a Gaussian frequency
regularization and then by the weighted Plancherel bounds below. The sign
of a right-hand residue is negative in this convention. This gives, near
zero on the original physical side,

    b-b_mass=B0 z^4+C_* z^(4+e_*)+r,
    C_*=-Ghat(4+e_*)/B'(4+e_*),
    Delta=B0+C_*=B0[1+B(4)/(e_* B'(4+e_*))]
                              -M_G(4+e_*)/B'(4+e_*).             (4)

Here r belongs to weight sigma with the derivatives quantified below.
The coefficient and remainder are thus defined for the actual fixed-seed
solution. All unknown global dependence is in the actual G, not a freely
chosen forcing. Formula (4) does not evaluate that dependence.

Avoid forming B0+C_* from separate rounded numbers. Since B(4+e_*)=0,

    B(4)+e_* B'(4+e_*)=integral_0^e_* t B''(4+t)dt.

With B=4D0 d0 q R0, the displayed derivative bounds imply
sup_(0<=t<=e_*) |B''(4+t)|/B'(4+e_*)<4 and B'(4+e_*)>0. Therefore

    |Delta| <= 2|B0|e_*+|M_G(4+e_*)|/B'(4+e_*).    (5)

The numerical response factor is of Einstein size, not 1/K. This is an
actual cancellation in the original-equation extraction formula. It does
not assume that the global moment is small.

## 5. Quantitative normal remainder, including all frequencies

On s=sigma+it, for |t|<=2 use |e|>=delta and the bounds above; for |t|>=2
use the real-part bound. Since |e q|>=delta(1+t^2)^(3/2),

    |B(s)| >= c3 (1+t^2)^(3/2),
    c3=4|D0| A_T delta(1-5000/L).                    (6)

For |t|>=4 the stronger estimate is
Re[(q+2)H]<=61-5t^2, so L-Re[(q+2)H]>=L/2+4t^2.
Using |eq|>=|t|^3, (1+t^2)^2<=(17/16)^2 t^4 and the arithmetic-geometric
mean inequality gives a second whole-line bound

    |B(s)| >= c4 (1+t^2)^2,
    c4=min(c3/sqrt(17),
            8|D0|d0(16/17)^2 sqrt(2L)).             (7)

For g_c(x)=exp(-sigma x)[G(exp x)-g0 exp(4x)], Plancherel yields

    ||z^-sigma r||_H6 <= C3 ||g_c||_H3,
    ||z^-sigma r||_H6 <= C4 ||g_c||_H2,
    C3=1/c3, C4=1/c4.                              (8)

Use the smaller APPLICABLE bound; H3 data must actually be bounded before
using C3. Merely renaming a data norm is not an improvement. Fixed-seed
finiteness was addressed in section 2, its numerical value was not.

For 0<z<=a=.1, Fourier Cauchy--Schwarz gives
||partial_x^j w||_infinity<=||w||_H6/sqrt(2), j<=5. Expanding
the falling polynomials p_j(sigma+partial_x) and dividing by the existing
j!100^j shows

    A5,a(r) <= a^delta/[sqrt(2)(1+log(1/a))]
                            min(C3||g_c||_H3,C4||g_c||_H2).      (9)

This is the previously declared canonical shape phase norm, not the full
physical admission norm. Scalar, mass, bolt, interior and X-ball conditions
remain additional requirements.

## 6. Execute the known source contribution and expose the missing estimates

For a source split ONLY choose x=log z, x0=log a-1, chi0=1 below x0,
chi0=1-S(x-x0) for x0<x<log a and zero above, where
S(u)=35u^4-84u^5+70u^6-20u^7. This split changes no physical operator.
Write G=g0 z^4 chi0+G_R. Define the two unknown actual global errors

    E3=||exp(-sigma x) G_R(exp x)||_H3(R),
    Efar=||z^(-7/2) G_R||_L2((1,infinity),dz/z).     (10)

For the known term g_known=g0 exp(-delta x)(chi0-1), all derivatives and
the entire infinite tail can be bounded without quadrature truncation:

    ||g_known||_Hm^2 <= g0^2 a^(-2delta) {
       exp(2delta) sum_(j=0)^m binom(m,j)
           integral_0^1 |(partial_u-delta)^j S|^2 du
                       +(1+delta^2)^m/(2delta)}, m=2,3.          (11)

The integrals are exact rational polynomial integrals. This is a known
piece of the actual source decomposition; it does NOT set G_R to zero.
Its regular moment obeys

    |M_known(4+e_*)| <= |g0|(1+log(1/a))
                                  exp(e_*(1+log(1/a))).          (12)

Cauchy--Schwarz, with the two different end weights, then gives

    |M_G(4+e_*)| <= |M_known|
          +E3/sqrt(2(delta-e_*))+Efar/sqrt(2(1/2+e_*)).          (13)

Combine (5), (9), (11), (13) with the previous paired-profile bound:

    A5,a(b-b_mass) <= |B0|e_*+
       [2|B0|e_*+(|M_known|+c_o E3+c_f Efar)/B']/(1+log 10)
                   +a^delta C3(E3+||g_known||_H3)
                                             /[sqrt(2)(1+log 10)]. (14)

Freeze a sensitivity test of (14) at E3<=1e-11 and Efar<=1e-10, targeting
the existing diagnostic shape cap 1e-8. Also execute the H2 alternative
and quantify the single-extra-derivative sensitivity. These are explicit
targets for UNKNOWN actual errors, not demonstrated upper bounds. No
known upper bound on E3 or Efar is claimed to meet or fail them.

## 7. Decision and remaining scope

The executed work establishes a fixed-seed global coefficient formula and
an entire-frequency normal remainder certificate. It evaluates the known
endpoint-source part and the conditional thresholds. The actual G_R
contains unknown original-row, mass, kernel and commutator contributions;
its norm and moment have NOT been evaluated. Thus this is not an executed
global solution or a passed physical feasibility test.

Research decision: keep the construction and use this combined moment and
strong-data method, rather than subtract large fitted asymptotic amplitudes.
The single next calculation is to bound the ACTUAL G_R in (10), including
its original-row global terms, against the two numerical targets. The
one-extra-derivative method is quantitatively motivated by (8), not a
normalization change. Its nonlinear-neighborhood analogue remains unproved.
Do not optimize already passing cold/scalar components to evade this test.

The quantum/contact screen P_hat_rest=0, E_rest<2.123e-42<1.287e-40 retains
its sufficient-only scope. At rho=1e-8, Y0+1e-8 Z0+5e-17 Z2<=1e-8 and
Z0+1e-8 Z2<1 remain open, as do the unproved allocations Y0<=5e-9,
Z0<=.25, Z2<=5e7. C2 on bounded admitted stationary families is not
X-ball admission. Exact compatibility, physical matching/causal realization,
semiclassical validity/noise and all-sector stability stay unresolved.

Remaining to-do list: bound the actual original-row global remainder and
moment, then control the full correction and nonlinear admission; exact
compatibility; physical matching and causal realization; semiclassical
validity/noise and all-sector stability; charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
