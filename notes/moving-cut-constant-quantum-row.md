# The full constant quantum row at a moving cutoff

2026-10-21, E003/S168. Ordinary argument before focused verification.
The complete canonical constant-response row, including its two-vertex
part, has a finite translated-logarithmic weak cutoff limit on squared
auxiliary energy-class shapes. Its boundary contribution cannot supply
the persistent oscillation required by the coupled trace calculation.
This does not evaluate the remaining two-point kernel or its derivative
contacts. The full pure-shape cutoff limit and physical branch exclusion
remain unproved.

Keep the fixed charged Euclidean seed, prescribed period, paired domain,
original exponential mass chart, full angular tower and common subtraction.
All traces have the retained per-exterior paired normalization. Reuse
[the cubic-test weak lemma](energy-domain-and-chart-contact.md) and
[global locality](moving-cut-global-quantum-locality.md), without changing
the physical domain or adopting an action-averaging prescription.

## 1. The exact local mass representation

On an isothermal collar at a cut, write the canonical metric as
exp(2sigma_d(x))(dx^2+dtau^2), x>=0, sigma_d(0)=0. The time coordinate
is rescaled together with x; the prescribed physical period is retained.
In the flat half-density representation put

    Dflat=N partial_x+T partial_tau,
    N=-tau_1 tensor 1, T=-tau_2 tensor sigma_2,
    chi=tau_3 tensor sigma_1, Pi_-=(1-chi)/2,
    m_d(x)=M exp(sigma_d(x)).                          (1)

The accepted conformal map preserves Pi_- on the boundary. For the
nonlocal mass determinant use Dflat+m_d, with the fixed local measure
and restoration terms still allocated as in the original functional.

Extend m_d and the scalar shape evenly. Reflection x->-x followed by
chi commutes with Dflat+m_d, since {N,chi}=0 and [T,chi]=0. The
parallel involution C=1 tensor sigma_2 commutes with N,T,m_d and
anticommutes with chi. Hence C interchanges the two reflection sectors
while commuting with the operator and every even scalar insertion.
Their determinant variations are equal. Equivalently, the half-domain
nonlocal functional and its scalar variations are half the evenly
varied double functional, with the same paired normalization.

This identity holds for a positive-cut Lipschitz even mass coefficient;
it uses a bounded scalar mass perturbation of Dflat, not a smooth-double
metric expansion. For small smooth finite-cut variation parameter the
mass stays positive, or in a fixed right sector for a complex parameter,
so the resolvents exist. Symmetry-preserving regulation gives the same
identity before the common subtraction. Local restoration is retained
separately, and the stationary phase remains the previously fixed one.

In particular the double representation keeps BOTH terms

    Tr(S m v S m v)-Tr(S m v^2), S=(Dflat+m)^-1.       (2)

For time-dependent v, use its even extension at the cut. This is an
exact representation of the nonlocal response, not an identification
with the constant-mass normal model. The seed mass has a linear |x|
jet at positive cut; that is still present in S.

## 2. A heat weight for the entire constant row

Return to the actual curved fixed-seed realization A=Dsharp, with
A^*=-A, R_M=(A+M)^-1 and P=M^2-A^2. A preserved grading anticommutes
with A and commutes with scalar smearings. Odd local spin traces vanish.
For the uniform exponential orbit-Weyl direction, the exact local row is

    M^2 tr(R_M^2)-M tr R_M
                  =2M^4 tr P^-2-2M^2 tr P^-1.         (3)

This is the Hessian operator applied to the constant function, not just
the second mass vertex. It can be paired with an arbitrary scalar U;
U need not commute with A. The retained fixed-domain conformal map is
what justifies the mass insertions in (3).

Let k_d(t;x,x) be the actual orbit heat diagonal of -A^2, with boundary
contacts understood in the common scheme. Define

    Aang(t)=sum_(n>=1) d_n M_n^2 exp(-t M_n^2),
    B(t)=2[t sum d_n M_n^4 exp(-t M_n^2)-Aang(t)]
        =-2 partial_t[t Aang(t)].                     (4)

Then the nonlocal part of the full constant row is FP integral B(t)k_d dt.
Its fixed local measure/restoration row must be added. Equation (4)
evaluates the mass derivative exactly; no differentiated remainder
estimate is invoked.

The retained sphere expansion Aang=2C/t^2+O(1), with its differentiated
smooth heat expansion, gives

    B(t)=4C/t^2+O(1).                                 (5)

Only the leading coefficient equals twice that of Aang. In general
B is not 2 Aang. Counting N(E)<=4CE^2 gives

    sum d M^2 exp(-tM^2)<=12C/t^2,
    sum d M^4 exp(-tM^2)<=40C/t^3,
    |B(t)|<=104C/t^2<120C/t^2,
    |B(t)|<=736C t^-2 exp(-M_1^2 t/2)
                                      <800C t^-2 exp(-M_1^2 t/2). (6)

The last line extracts half the mass gap and applies the preceding
bounds at t/2. It is uniform in the moving cut because the angular
tower is unchanged. The benchmark M_1>18 recovers exp(-162t).

## 3. The leading nongeodesic and cubic boundary profiles

Use inward proper distance y from the cut, normalized lapse

    a(d+y)/a(d)=1+c1(d)y+c2(d)y^2+c3(d)y^3+...,
    c1(d)=a'(d)/a(d)=O(d), c3(d)->c3(0).              (7)

The half-density square has the retained form
-partial_y^2-a^-2 partial_tau^2-gamma_y gamma_tau(a^-1)' partial_tau,
with complementary Dirichlet/Neumann conditions. Its balanced image
spin trace has no weight-zero profile. At linear order in an odd jet
c_j y^j, j=1 or 3, comparison with the smooth continuation gives

    K_j(t,y)=c_j/(2pi t^2) integral_0^t
          D_j(y,sqrt(2s(t-s)/t)) ds,
    D_j(y,sigma)=E[|y+X|^j-(y+X)^j], X~N(0,sigma^2). (8)

The connection insertion is traceless at this order. The outside
density is identical on the observed physical half. This is the same
rank-two calculation as the retained cubic boundary layer, now also
keeping the positive-cut linear jet.

Since D_j>=0, Tonelli gives an independent coefficient calculation:

    integral_0^infinity t^-3 D_j(y,sqrt(vt)) dt
       =6v^2 integral_y^infinity (s-y)^j/s^5 ds
       =6v^2 y^(j-4) Beta(j+1,4-j),
    integral_0^1 [2z(1-z)]^2 dz=2/15.                 (9)

Multiplying (8) by the leading B(t)=4C/t^2 therefore yields

    Theta_linear(y)=2C c1(d)/(15pi) y^-3,
    Theta_cubic(y)=2C c3(d)/(5pi) y^-1.               (10)

These are positive-distance coefficients of the COMPLETE constant row,
with the Euclidean Hessian sign in (3). Their distributional extensions
are the fixed-scale finite parts and contacts of the original scheme.
They are not independently renormalized boundary actions. At positive
cut, other nonlinear jet combinations of weights two and three remain;
(10) is not asserted to be the complete finite-cut row.

For g2=B(u)dtau^2+b^2 du^2/B(u), the jet conversion provides a check:

    c1=B'/(2b sqrt(B)), c2=B''/(4b^2),
    c3=B'B''/(24b^3 sqrt(B))+sqrt(B)B'''/(12b^3).      (11)

At u=0 with B'=0 and B'''=-12a_m, this gives the retained
c3=-a_m sqrt(B(0))/b^3. This conversion applies to the actual smooth
seed polynomial as well as its RN reference; no reference replaces it.

## 4. The complete fixed-seed row has the cubic-test structure

The earlier fixed-seed heat decomposition is a statement about k_d
before multiplication by its angular mass weight. Replacing Aang by B
preserves its proof: both weights are bounded by a constant times
t^-2 at small t and have the retained massive gap at large t. The
new leading coefficient is (5), not an assumed equality of the weights.
Thus the low boundary weights n=1,2,3 give fixed-scale finite parts
FP[y^-p log^j y], p<=3. Their coefficients and scaled cut derivatives
are bounded and converge. All higher-weight remainders remain uniformly
L1 after the same bulk subtraction; off-collar kernels converge.
Global completion has the independent ordinary limit already proved.

The finite local restoration and anomaly row has contacts through
smearing order two. The possible third derivative in
[Vassilevich (5.33)](https://arxiv.org/html/hep-th/0306138v3)
is proportional to tr chi and vanishes in the balanced pair. Uniform
constant-Weyl differentiation of the retained local coefficients does
not increase the smearing-derivative order. This keeps the full local
allocation, including bulk restoration and boundary charges.

Let Theta_d denote this COMPLETE canonical constant row. For every
continuous U(Y)=O(Y^3), initially smooth for Y>0, pair Theta_d with
U_d(y)=U(d+y). Put delta=exp(-Z0) and
psi_delta(d)=phi(log(delta/d))/d for a smooth compact test phi.
Transposing each finite part onto psi_delta gives, on Y=O(delta),

    O(delta^(4-p)(1+|log delta|)^m), p<=3,             (12)

after multiplication by U. A contact of order q<=2 costs
O(delta^(3-q)(1+|log delta|)^m). Both vanish. Away from this shrinking
collar, dominated convergence gives the limiting ordinary pairing;
the uniformly L1 remainder is handled by first removing a small collar.
This is precisely the previously proved weak cubic-test lemma, now
applied to the full constant row rather than the first-variation contact.

Every auxiliary energy-class shape satisfies v=Y^(3/2)V(-log Y),
V in H2, so U=v^2 meets this hypothesis. Polarization covers v1 v2.
Stationarity reduces a fixed temporal factor to its mean square.
No pointwise limit of (v^2)'' at the cut is assumed. At the limiting
geodesic boundary the derivative contacts vanish; the remaining value
contact pairs with U(0)=0. Finite interior and state contributions are
retained. Hence

    lim_(Z0->infinity) integral phi(Z-Z0)
       <Theta_(exp(-Z)),v(exp(-Z)+.)^2> dZ
                  =<Theta_0,v^2> integral phi, finite. (13)

In particular the complete constant row cannot supply a persistent
nonconstant cutoff sinusoid. It is not set to zero, and (13) does not
assert a pointwise cutoff limit or a physical averaging prescription.

## 5. What remains in the actual shape Hessian

Writing the massive bubble in row-subtracted difference form moves
its full constant row into Theta_d. Equation (13) now controls that
entire zeroth-order allocation, not just the chart contact. The
row-subtracted near-diagonal/reflected two-point kernel, its fixed
derivative contacts, and the original Weyl/LLL/classical shape partners
must still be compared with the normal representative. No independent
finite part for a bare bubble or tadpole is introduced by this split.

The exact even-mass representation (1)-(2) supplies a concrete local
resolvent problem for that next step. Its linear |x| background jet
cannot be discarded. No general fast-branch exclusion, physical mode,
nonlinear solution or stability follows from the present row theorem.

Verification checks paired reflection symmetries, the mass-row heat
identity with noncommuting scalar probes, the two independent odd-jet
Gaussian integrals, proper-coordinate jet conversion, angular-weight
enclosures and moving-cut finite-part pairings. Numerical fixtures are
diagnostics, not actual charged-kernel or shape-solution evaluations.
The first frozen run failed its direct j=3 heat-tail quadrature at
50 digits. The t^-3/2 tail is now integrated using the exact coordinate
u=1/sqrt(t) for t>=1; the tolerance is unchanged. The failed source and
raw run are preserved. This changes numerical integration, not (9).

Remaining to-do list: row-subtracted local curved quantum kernel and
derivative contacts with all local shape partners; physical boundary phase
and source compatibility; continuous causal constraints, exact nonlinear
existence and all-sector stability; exterior/state, charge/return, opening,
transport, preparation, separate creation and full audit.
