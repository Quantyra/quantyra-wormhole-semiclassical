# Outer graph feedback before estimating the unknown global remainder

2026-10-24, E003/S169. Ordinary argument before execution.
The full stable-wormhole objective remains active and incomplete; the audit
stays closed. This addresses an operator prerequisite of the selected global
remainder calculation. It does not evaluate the unknown correction or its
full remainder, and it does not replace the paired domain or physical source.

## 1. The circular estimate to avoid

In current-global-boundary-moment-return.md, G_R includes Rw and the original
cutoff commutator acting on the UNKNOWN full linear response. Finiteness of
the fixed-seed source's derivatives does not bound those terms numerically.
The displayed normal estimate H3 -> H6 is a valid consequence for given
data. Using H6 alone to estimate the feedback into H3 would lose derivatives.

Retain x=log z, sigma=4.005 and the full normal multiplier B(s). After the
mass, quartic particular and slow mode have been extracted as in that note,
write w=z^-sigma r and B_sigma=B(sigma+partial_x). All finite-dimensional
coefficients and original-row moments remain part of the coupled problem.
They are not chosen to make a residual small.

## 2. The actual principal feedback and an obstruction to the weaker norm

The full all-angular static shape response has the normalized local log
principal symbol

    -[G C Q(z)^2/(15pi B^2)] xi_z^4 log|xi_z|.

This follows by freezing the actual metric dz^2/Q+Q d tau^2 at an interior
point in the previously identified paired determinant symbol: the proper
radial frequency is sqrt(Q) xi_z and the canonical area density is unchanged.
The fixed Weyl contact cancels the explicit outer log z before the normal
comparison. Finite fourth contacts, log Q and density/coordinate derivative
terms remain in the nonlogarithmic or lower-order remainder. They are not
discarded. The angular sum producing log frequency is the entire tower.

Multiplication by z^4 and the original differentiated row add one D=z d_z.
Consequently its log principal difference from the normal B is

    -K rho(z) (i xi_x)^5 log|xi_x|,
    rho(z)=(Q(z)/beta)^2-1, K=G C beta^2/(15pi B^2).       (1)

The conformal and mixed rows have no quantum log principal symbol, by the
retained trace/Ward identity. Their local fourth-order contacts cannot
cancel (1). The fixed affine chart and full actual Hessian are unchanged.

Choose a smooth compact interval in .005<z<.0075 where rho is nonzero, and
real compact oscillatory directions n^-6 psi(x) cos(nx), with nu=0. These
are allowed identical-mouth interior metric directions; they vanish near
the actual boundary and bolt. They test a map, not the prescribed source.
Their H6 norm stays bounded, while the H3 norm of the feedback in (1) grows
as a nonzero constant times K n^2 log n. The localized log symbol gives
this leading asymptotic; the order-five nonlog terms are smaller by log n,
and lower-order terms by powers of n. There is therefore NO finite bound
from H6 to H3 for the full differentiated actual-minus-normal operator on
that interval. Making the collar coefficient tiny does not fix this order
loss. This is not an obstruction to the actual solution or to a graph norm.
No finite physical mode approximation or adjustable forcing is introduced.

## 3. Retain the full graph order with its physical quantum coefficient

Use log(2+|xi|) for the positive logarithmic Sobolev weight. Its difference
from log|xi| at high frequency is lower order and remains in the residual.
The full B, including all contacts, is never replaced by its principal part.
The preceding entire-frequency proof already gives

    ||w||_H6 <= C3 ||B_sigma w||_H3, C3<1801.

We now retain the complementary graph control

    K ||w||_(H8,log) <= Cq ||B_sigma w||_H3, Cq<2.328.    (2)

Here H8,log has multiplier <xi>^8 log(2+|xi|). K is the fixed physical
coefficient, not a rescaling of the source or a new adjustable parameter.

An elementary whole-frequency proof supplies the constant. For s=sigma+it,
t>=4 and X=Re(s-1), compare the digamma sum with its integral. The total
variation of y/(y^2+t^2) on y>=X is <=1/t, so

    Re psi(X+it)>=log(sqrt(X^2+t^2))-1/t.

Together with the retained kappa_M>7.8 and cotangent exponential bound,
this gives Re H(s)>=log(2+t)+7. The previous |t Im H|<5 then yields

    L-Re[(q+2)H]>=L/2+(t^2/2)log(2+t), t>=4.

Since B=-K(s-4)q[(q+2)H-L],

    K <t>^5 log(2+t)/|B(s)| <= 2(17/16)^(5/2), t>=4.

For |t|<=4 the earlier c3<t>^3 lower bound gives the bound
K C3 17 log 6, which is much smaller. Complex conjugation handles negative
t. The maximum of these two expressions is Cq. This is not a sampled
frequency certificate. The digamma and cotangent identities are those of
[DLMF 5.7.6](https://dlmf.nist.gov/5.7.E6) and
[DLMF 4.22.3](https://dlmf.nist.gov/4.22.E3), already reviewed for B.

In particular, for the concrete Fourier operator
P=K(partial_x)^5 log(2+|D_x|),

    ||P B_sigma^-1||_(H3 -> H3) <= Cq.                 (3)

This estimate retains exactly the derivatives missing from the H6-only
argument. It uses the original operator and physical coefficients; it does
not declare a norm-one inverse by changing units. The complete graph domain
and mass/slow asymptotic components are still required.

## 4. Execute the actual outer log-principal feedback bound

Only the OUTPUT coefficient is localized here. Let a=.01 and choose
chi_a=1 for x<=log a-1, chi_a=1-S(x-log a+1) on the next unit interval,
and zero for x>=log a, with

    S(u)=35u^4-84u^5+70u^6-20u^7.

This cutoff is W3,infinity, sufficient for multiplication on H3. It is
not an input cutoff in the H8,log graph and imposes no physical wall.
Set alpha(x)=chi_a(x)rho(exp x). Define a_j=||partial_x^j alpha||_infinity.
The exact integer Sobolev product rule gives

    ||alpha v||_H3 <= (a0+3a1+3a2+a3)||v||_H3.       (4)

For example the terms with one derivative on alpha have squared derivative
weights 3,12,9, bounded by 9 times the H2 weights; the two-derivative terms
are bounded by 9 times H1. Thus (4) has no untracked product constant.

Use the full actual finite-temperature polynomial

    Q/beta=1+(u/beta)z^2+(-4-2u/beta-tau/beta)z^3
        +(3+u/beta+epsilon/beta)z^4
                              -2epsilon z^5/beta+epsilon z^6/beta.

Square it to obtain rho=sum rho_p z^p. For j=0,...,3 set
r_j=sum |rho_p| p^j a^p. These bound all derivatives of rho on the entire
collar, including z=0 by continuity. Certify |chi'|<=3, |chi''|<=10,
|chi'''|<=60 by exact polynomial Bernstein subdivision on [0,1]. Then

    a0<=r0,
    a1<=r1+3r0,
    a2<=r2+6r1+10r0,
    a3<=r3+9r2+30r1+60r0.                            (5)

The executed target is

    ||alpha P B_sigma^-1||_(H3 -> H3)
         <= eta_pr=Cq(a0+3a1+3a2+a3)<.008.            (6)

This is a WHOLE-FREQUENCY operator bound for the identified log-principal
component of the actual outer difference, with full coefficient derivatives.
It is not the norm of the complete difference. Off-diagonal kernel terms,
nonlogarithmic fifth contacts, scalar mixing, coefficient derivatives in
the operator, commutators, mass/slow returns and regular-bolt data remain
outside (6). The existing .0022 action-form estimate is in a different
norm and cannot be added to (6) as if it controlled these graph terms.

## 5. What would close the circular estimate, and what has not been measured

Let k3 be the previously certified norm of the known endpoint source,
k3<1.613e-27. A complete coupled remainder argument would have to produce
an inequality of the form

    E3 <= S3+(eta_pr+eta_other)(E3+k3).               (7)

Here eta_other must bound ALL remaining feedback in the same graph and
actual domain. S3 must bound the remaining affine/global contributions,
including finite-dimensional mass/scalar/slow components through their
original equations. Neither is known. Formula (7) is a requirement, not
an established closed inequality or a choice of a new quantum force.

For a concrete allocation, reserve eta_pr<=.01, demand eta_other<=.24,
and use the strict affine target S3<=7.49e-12. Then (7), IF established,
would imply E3<1e-11. The independent Efar<=1e-10 requirement stays open.
Substitute the actual bound (6) as well and report both results. These are
thresholds for unknown estimates, not known bounds declared to pass.

## 6. Research decision and scope

Executed: reject the H6-only feedback bound by an actual principal-order
argument; retain a quantified normal graph estimate; and bound the entire
log-principal outer feedback component using the actual Q polynomial.
This closes one operator requirement for the global response. It does not
compute actual E3, Efar, S3, eta_other, a full correction, or a preserved
nonlinear ball. The preceding conditional targets and moment formula remain
valid; they have not become actual upper estimates.

Research decision: retain the construction and use the full graph order
when estimating its remainder. The single next calculation is the complete
remaining coupled graph feedback and affine term in (7), including the
original boundary/global contributions, against eta_other<=.24 and
S3<=7.49e-12. It is proposed, not executed. A separate-coefficient or weak-
form bound cannot substitute for this calculation. No cold/scalar source
optimization is selected by this principal-order test.

The earlier sufficient quantum/contact screen remains passed, with
E_rest<2.123e-42<1.287e-40 and unmeasured sign. At rho=1e-8 the unclosed
conditions remain Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1, with
unproved allocations Y0<=5e-9, Z0<=.25, Z2<=5e7. Fixed-seed graph control
does not prove uniform nonlinear X-ball admission. Exact compatibility,
physical matching/causal realization, semiclassical validity/noise and
all-sector stability remain unresolved. No controlled wormhole, nonexistence
or instability result follows from this operator calculation.

Remaining to-do list: complete the actual coupled graph feedback and global
affine bound, actual E3/Efar and controlled correction; nonlinear admission
and exact compatibility; physical matching and causal realization;
semiclassical validity/noise and all-sector stability; charge/return,
opening, transport, preparation, separate creation and the full objective.
