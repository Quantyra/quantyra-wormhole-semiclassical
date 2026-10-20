# A quantitative source bound for the complete physical mass column

2026-10-20 HST, E003/S166 supporting S139/S074. Ordinary argument
and verification contract before execution. **Mass-column image
controlled; annulus matching and the full inverse remain open.**

The physical mass cubic is outside the fast input space. Its
normal response cancels exactly, however, and the actual-minus-normal
terms have enough gain to have a finite source norm. This gives a
bound for the FULL mass-column image, including output outside
the collar and all massive global/contact terms.

## 1. Chosen column and source norms

Use the action, domain, parameters, seed and half-density similarity
of [outer localization](outer-fast-localization-and-tails.md).
Fix epsilon<=1e-160, r=epsilon/32 and L=epsilon/8=4r.
With its beta(6,6) radial cutoff c_L, define in (theta,b) variables

    U_L(Y)=c_L(Y)Y^3, v_m=(-1/12,1), ehat_m=U_L v_m,
    e_m=M_epsilon^-1 ehat_m.                              (1)

Here e_m is in the original (theta,b) variables; conversion to
(S,T) is S=theta-b/6, T=b/2. Its support lies in Y<=2L=epsilon/4,
so it vanishes near the bolt. Its leading b coefficient is one
and its leading theta coefficient is -1/12. M_epsilon=1+O(Y^2)
leaves the cubic unchanged; the change from another cutoff
representative is a fast field. No new mass datum is added.

Let Lhat be the full fixed-canonical-area Hessian in these
auxiliary variables. Thus

    Ahat=C u^4/(4pi B^4) Lhat,
    lambda^-2 Ahat=C u^4/(4pi xi^2) Lhat,
    Q(ehat_m)=||C^-1 Lhat ehat_m||_0,
    ||f||_0^2=A0^-1 integral_gbar abs(f)^2,
    A0=beta sqrt(C lambda), xi=lambda B^2.

The canonical norm here is the two-row Hilbert norm on the
WHOLE physical half. The principal result is

    Q(ehat_m) <1e86 L^-1/2.                              (2)

For the original physical output (f_S,f_T), put
Fdisk(f)=||u^-2 f||_L2(dmu4) and gamma0=4pi beta B^3.
Then (2) implies the complete physical source bound

    Fdisk(A e_m) <1e88 lambda^2 sqrt(gamma0) L^-1/2.       (3)

In (3), A e_m means the same physical direction after the
displayed (theta,b)-to-(S,T) conversion. Both full source
equations, the off-shell density contact and global tails are
included. This is a bound for a specified mass representative,
not an inverse norm or a nonlinear contraction.

## 2. Exact normal cancellation and the cutoff source

Write N0 for the complete dimensionless normal operator and
q_s=s(s-3), K0=60pi/g+11/3. The matrix from the
[complete local assembly](outer-complete-disk-normal-comparison.md)
satisfies exactly

    N0(3) v_m=0.                                       (4)

The Mellin logarithm is regular at s=3 and P4(3)=q_s(3)=0.
The remaining entries cancel on (-1/12,1). The accepted Mellin
formula therefore gives N0[Y^3 v_m]=0 as an even homogeneous
distribution on Y>0. This is not an assertion that Y^3 belongs
to X_(7/2). It includes the cancellation of the massive 1/Y
term with its Weyl contact before estimation.

We first prove

    ||N0 ehat_m||_Z <1e80 L^-1/2,                       (5)

where Z is the weighted two-row L2 norm on the full Y half-line.
The normal operator is dilation invariant, so it is enough to
take L=1; its output norm then scales by L^(3-7/2)=L^-1/2.

For Y<1/2 subtract the full homogeneous power using (4).
The remaining input is supported in Z>=1. The normal separated
kernel is -Y^4[abs(Y-Z)^-5+(Y+Z)^-5]/(40pi^2).
Since Z-Y>=Z/2,

    abs(N0 ehat_m(Y))<=33Y^4 integral_1^infinity Z^-2 dZ
                    =33Y^4.

Its weighted norm on this interval is below 24. For Y>4,
the input is supported in Z<=2. The same kernel gives
abs(N0 ehat_m(Y))<=132/Y, whose weighted norm is below one.
Local terms vanish in both separated calculations.

On 1/2<=Y<=4 use Xi=1-c_(1/8). The field Xi U_1 v_m is
supported in 1/8<=Y<=2 and is in the fast mixed space. In
log coordinates its weighted profile is
Y^-1/2(1-c_(1/8))c_1 v_m. The two cutoff factors have log
derivatives through five below 1e14, the profile factor below
three, the interval has length log16<3 and abs(v_m)<2.
Leibniz gives every derivative through five below 1e32 per
component. Thus its H5 pair norm is below 1e34 and its X
norm below 1e35, using H4log<=2H5. The accepted forward
normal norm below 1e38 gives a source norm below 1e73.

The remaining input (1-Xi)U_1 is supported in Y<=1/4.
On the observation interval the separated kernel bounds its
output by 33/(1024Y), with weighted norm below one hundred.
Combining the three regions proves (5). No fast norm of the
uncut homogeneous mass has been used.

## 3. The massive difference on a cubic input

The fast comparison cannot be applied directly to (1).
Instead use its exact contact/Taylor allocations on this input.
The beta derivative bounds give, on the support of U_L,

    abs(partial_Y^k U_L)<=1e6 Y^(3-k), k<=3,
    abs(D^k U_L)<=1e18 Y^3, k<=4, D=Y partial_Y.

Put q_L=m^-1 U_L, m=sqrt(j). The preceding analytic inverse
disk gives abs(partial_Y^a m^-1)<1e5 for a<=3 on this support.
Product rule therefore gives abs(partial_Y^k q_L)<=1e15Y^(3-k).
With j<2, abs(h')<=2Y, abs(h'')<=2, the following canonical
norm bounds hold:

    ||Y^-2 q_L'||_0 <4e15 sqrt(L),
    ||Hess q_L||_0 <1e17 L^(3/2),
    ||nabla^3 q_L||_0 <1e18 sqrt(L),
    ||q_L||_0 <10 L^(7/2).                              (6)

For the actual-to-geodesic difference use the exact allocation
from [the strong curved remainder](outer-strong-curved-quantum-remainder.md):

    E q=p0 q+A1.grad q+A2:Hess q+T3 q+B0 q,
    abs(A1_Y)<1e60 C/Y^2, abs(A2)<1e30 C,
    ||T3 q||_0<=2e6 C||nabla^3 q||_0,
    ||B0||<=8e24 C,
    abs(p0)/C<=1.6e25(1+1/Y).

At outer output the physical weight converts Y^4 to L2(dY).
The output m, density and measure conversions cost less than
four after division by C. Equations (6), including the direct
p0 estimate on q_L, give the following safe source reserves:

| Actual-to-geodesic contribution | Weighted norm bound |
| --- | --- |
| First contact | 1e77 sqrt(L) |
| Second contacts | 1e48 L^(3/2) |
| Taylor remainder | 1e25 sqrt(L) |
| Constant contact | 1e30 L^(5/2) |
| Bounded remainder | 1e30 L^(7/2) |

Their sum is below 1e80 sqrt(L). For the conjugated geodesic-to-flat
difference, use its coefficients abs(c_k)<1e40 C/Y, k<=2,
Taylor bound 1e7 C||U_L'''||_2 and bounded tail 1e20 C.
The same explicit cubic derivatives give a bound below
1e50 sqrt(L) after the physical density.

Finally the exact flat-tower correction has the GLOBAL Fourier
bound abs(K-K_infty)<=1000C(1+p^2)^(3/2). The even U_L is
C2,1 and in H3, with H3 norm below 1e8 sqrt(L). Thus this
last contribution is below 1e12 sqrt(L). Enlarge the combined
massive-difference reserve to 1e81 sqrt(L).

These estimates extend to the mass representative by cutting
it off smoothly for Y<a and then sending a to zero. The H3,
Y^-2 first-derivative and Y^-1 second-derivative errors are
O(sqrt(a)); the bounded and constant terms converge as well.
Use these difference allocations for the limit. This does not
assert convergence of the truncated normal or full images in Z.
The full mass action uses the accepted
[mass-extended realization](stationary-global-weighted-map.md)
and (4). The separate massive and Weyl terms are not asserted
to converge in the weighted source norm before cancellation.

## 4. Remaining local error and the full outer source

The complete local coefficient ledger gives a sum below
6e45 eta(1+abs(O)), eta=1e7Y^2, O=-log(sqrt(xi/C)Y).
Every local tangent derivative is D^k with k<=4. The derivative
bound above and Y^2(1+abs(O))<4Y give

    ||(L_local-L_local,0)ehat_m||_Z <1e74 sqrt(L).       (7)

Here the local operators are in the dimensionless normalization
inside rho in the preceding full assembly. On Y<epsilon the
exact full equation is

    lambda^-2 Ahat ehat_m
       =massive_difference +rho(L_local-L_local,0)ehat_m
                               +rho N0 ehat_m.

Use rho<3/2, (5)-(7) and L<1. This proves

    ||1_(Y<epsilon)lambda^-2 Ahat ehat_m||_Z
                                          <1e83 L^-1/2. (8)

All actual massive global and constant terms are already in the
first term and are not added again. The pair's theta component
and every finite contact are retained in the local ledger.

## 5. The global image and physical forcing norm

On the collar, rho>1/2 and j<2. The exact density relation and
weight identity imply

    ||1_(Y<epsilon)C^-1 Lhat ehat_m||_0
         <=100 ||1_(Y<epsilon)lambda^-2 Ahat ehat_m||_Z.

For Y>=epsilon the input is supported in Y<=epsilon/4, so
only the actual nonlocal massive entry can contribute. The
separated actual rows from the localization proof, with gap at
least epsilon/2, are below 1e29 C epsilon^-4. Input/output
half-density factors and the b normalization are bounded.
Using ||U_L||_0<10L^(7/2) gives

    ||1_(Y>=epsilon) C^-1 Lhat ehat_m||_0
                       <1e32 epsilon^-4 L^(7/2)
                       <1e32 L^-1/2.                    (9)

This is the full canonical output norm, including the bolt,
global returns and the actual thermal/domain kernel. Adding
(8)-(9) proves (2), with reserve below 1e86.

For completeness the physical source conversion is exact before
the harmless matrix factors. Since dmu4=(4pi B^4/C)u^-4 dvol_gbar,

    ||u^-2 [C u^4/(4pi B^4)] Lhat ehat_m||_L2(dmu4)
       =C^(3/2)sqrt(A0)/(2sqrt(pi)B^2) Q(ehat_m)
       =lambda^2 sqrt(gamma0) k^-7/2/(4pi) Q(ehat_m),
    k=sqrt(xi/C)>.4.

Return to (f_S,f_T) using the inverse dual matrix, whose norm
is below three, and M_epsilon^-1, pointwise below two. Since
k^-7/2<25, the total conversion factor is below 100. This
proves (3) without changing the original physical forcing norm.

## 6. Applying localization to a mass-extended field

The preceding mass coefficient is m_Y=-3(xi/C)^(3/2)m_c,
with abs(m_Y)<abs(m_c)/4. Thus vhat_fast=vhat-m_Y ehat_m is
in the fast space; its source is the full source minus the
actual mass image just bounded. The accepted numerical bound
abs(m_c)<=5.40144e14 g lambda^-2 gamma0^-1/2 Fdisk gives

    Fdisk(A v_fast)
       <=[1+2e102 g L^-1/2] Fdisk(A v).                  (10)

The cutoffs in the chosen annulus also give
||sigma ehat_m||_X<1e35 L^-1/2 and V0(ehat_m)<10L^(7/2).
The outer factor c_(8r)=c_(2L) is one on the entire support of
U_L. Thus sigma U_L=(1-c_(L/16))c_L Y^3 has only two effective
cutoff factors. After scaling by L its weighted support is
[1/16,2], of log length below four, and Y^-1/2 is at most four.
Product rule through five gives X norm below 1e35 by the same
H5 argument as section 2. Combining (8) and these profile norms with
the preceding fast localization theorem yields

    ||chi(vhat-m_Y ehat_m)||_X
       <=2e20 ||chi fhat||_Z +1e-40 ||sigma vhat||_X
          +1e2000 r^-4 V0(vhat)+1e2005 abs(m_Y)L^-1/2.    (11)

This is now a collar estimate for the physical mass-extended
field. The large mass/value factors are retained; they do not
constitute a useful nonlinear contraction. The transition and
intervening annulus, higher reconstruction and full residual/
inverse composition are still required.

## 7. Frozen verification and remaining work

Freeze this proof and checker before execution. Pin the complete
normal assembly, fast localization, actual massive allocations,
flat multiplier and mass-coefficient receipts. Verify the exact
normal mass null vector, cutoff derivative and norm budgets,
all infinite-end integrals, local/contact and global-tail sums,
physical measure conversion and mass-extended localization.
Preserve every run and failed source snapshot. The finite checks
support the ordinary continuation, distribution and operator
arguments; they do not formally certify a global solution.

Remaining to-do list: transition and intervening annulus/bolt
matching; higher York reconstruction and full inverse/residual
composition; nonlinear family and signed compatibility;
exterior/state, charge, return, opening, stability, transport,
preparation and separate creation.
