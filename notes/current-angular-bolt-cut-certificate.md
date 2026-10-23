# Full-tower bolt cut with its actual finite-cut adjoint terms

2026-10-23, E003/S169. Ordinary argument before frozen execution.
Previous goal turn: progress. The full stable-wormhole goal is incomplete;
the feasibility audit stays closed and its physical targets are unchanged.

The existing fixed-cut tip allowance is not angularly summable. Resolve
that endpoint obstacle by a frozen mass-dependent numerical cut and a
uniform bound on the actual auxiliary action adjoint. Include its finite-
cut boundary terms, not only the omitted tip integral. This concerns the
canonical mass action difference and its thermal restoration. The common
vacuum/local restoration and the global angular remainder are still required.

## 1. Same geometry, action and boundary realization

Use the current finite seed, original scale and paired domain of
[the endpoint contract](current-endpoint-thermal-restoration.md).
For each massive angular level,

    M_n^2=n(n+q)/C, d_n=2(q+2n), q=2000^4, n>=1,
    omega=j+1/2, j>=0, a=sqrt(Q)/gamma,
    K=sqrt(M^2+omega^2/a^2),
    c=-omega M a'/[2(omega^2+M^2 a^2)],
    X'=-2KX-c(1+X^2), |X|<=1,
    g(X)=X^2/(1+X^2), A=2K+2cX,
    ell'=A ell-4K g'(X), ell(L)=0.                    (1)

Prime means canonical proper rho differentiation from the bolt to the
physical paired seam. ell is the AUXILIARY action-integration adjoint,
not the gravitational bolt adjoint. Its seam value is already derived
from the varied paired characteristic; it is not a new wall condition.

For a direction define the inner norm

    U_i=sup_(0<=rho<=delta/M) max(|u_i|,|rho u_i'|), delta=1e-3.

First prove the estimates for U_i<=1. Gradient and mixed Hessian errors
then scale by U_i and U_i U_j. Smooth even fields satisfy the endpoint
contract; the same integral equations extend it to the C1 fields with
the rho^2 log rho or rho^2 log log(1/rho) adjoint behavior. Their displayed
inner norm is finite. No general closed-boundary C5 hypothesis is added.
The norm of any proposed gravitational trial must still be checked.

## 2. A global bound for the actual auxiliary adjoint

The actual seed has Q>0 on (0,z_h), Q(z_h)=0, Q'(z_h)=-2gamma and
Q'(0)=0. Write beta,u0,eps,tau for its polynomial parameters, with
6beta+u0+eps=1. Directly,

    2-Q''=(1-z){12beta(1+3z)+12u0 z
                  +eps(2+2z-10z^2+30z^3)}+6tau z.    (2)

The current box has 12beta-44|eps|>1.9, u0>0, tau>=0 and z_h<1.
Hence Q''<=2 throughout the ACTUAL interval. Put E=Q'^2/4-Q.
On each Q'<0 component E'=Q'(Q''/2-1)>=0; on Q'>0 it is <=0.
At a critical point E=-Q<=0, at the seam E=-Q(0)<0, and at the
bolt E=gamma^2. The polynomial has finitely many such components.
Thus E<=gamma^2 everywhere, and exactly

    |a'|<=sqrt(1+a^2)<=1+a.                          (3)

This uses the correlated horizon equation, not subtraction of tiny
interval values near the root. Set v=Ma/omega. Since both
v/(1+v^2)^(3/2) and v^2/(1+v^2)^(3/2) have maximum 2/(3sqrt(3)),

    |c|/K<=1/(3sqrt(3)) (1/omega+1/M)<.4,
    A>=1.2 K, |g'|<=9/(8sqrt(3))<.65.               (4)

Here M>=18 and omega>=1/2 are essential. Backward variation of
constants from ell(L)=0 now gives |ell|<=4(.65)/1.2<3. This bound
holds along the entire physical geometry and for every thermal and
massive angular mode. No unknown outer coefficient is assigned a value.

## 3. Sharper approach to the regular bolt

The same factored chi geometry as the endpoint contract, with the actual
h enclosure, gives .999<a/rho<1.001 and |a'|<2 on the retained chi
cover [0,.001]. It covers rho0=delta/M for EVERY M>=18. On rho<=rho0,

    omega/(1.001rho)<=K<=1.002omega/rho,
    |c|<=M/omega,
    |X|<=2M rho/[omega(1+omega)],
    |Y_i|<=10M rho/[omega(1+omega)],
    |Z_ij|<=14M rho/[omega(1+omega)].                (5)

The last three are the already proved regular integral-equation barriers;
they retain both initial mass sensitivities. Inserting (5) into A gives

    A>=omega/rho [2/1.001-4delta^2/(omega^3(1+omega))]
      >=1.99omega/rho>=.99/rho,
    |4K g'|<=16.032 M/(1+omega)<11M.                (6)

Let k=99/100 and t=rho/rho0. Use |ell(rho0)|<3 in (1) and integrate
the comparison with A>=k/rho backward. Since 1-k=.01,

    |ell(rho)|<=3t^k+11M rho0 t^k/(1-k)
               <=4.1t^k.                          (7)

This bound includes the lowest thermal channel's logarithm without
declaring its unknown homogeneous coefficient zero.

The exact finite-cut identities include +ell(epsilon)Y_i(epsilon)
and +ell(epsilon)Z_ij(epsilon). The telescoping sum
sum_(j>=0) 1/[omega(1+omega)]=2 and (5),(7) therefore give, per angular
state and with BOTH thermal signs already in the action convention,

    E_boundary,H<=114.8 delta^(-k) (M epsilon)^(1+k). (8)

The gradient constant is 82, so (8) covers it as well for unit norms.
This is an actual bound on the finite-cut boundary term. The earlier
tip integral alone did not bound it. A code may evaluate these products
more sharply; if it omits them, it must include (8) as an error interval.

## 4. Sum the endpoint errors over every angular mode

Keep the numerical cuts FIXED when taking direction derivatives. Choose

    m0=18, epsilon_n=eta (m0/M_n)^3, eta=2e-17.       (9)

They are computational cuts, not new boundaries. Since M epsilon_n
<=eta m0<delta, all local hypotheses hold uniformly. The already proved
omitted difference and local thermal-restoration integrals are

    E_diff,n<=1550(6-pi^2/2)(M epsilon_n)^2,
    E_th,n<=4(M epsilon_n)^2[log(2/(M epsilon_n))+7/2]. (10)

The first sum includes every j exactly, and the second is the complete
nonzero-winding restoration. No zero-temperature or mass-independent
metric term is erased. The LLL has no mass vertex; its independent Weyl
and physical allocations are still present in the full operator.

For N(E)=sum_(M_n<=E)d_n the exact angular counting gives N(E)<=4CE^2
and N(E)=0 below m0. Indeed its value through n=m is 2m(q+m+1),
while CE^2>=m(m+q). Stieltjes integration then gives, for r>2,

    S_r=sum d_n/M_n^r<=4Cr/(r-2) m0^(2-r).          (11)

Let L0=log(2/(eta m0))+7/2 and r=2(1+k)=3.98. The decreasing
log-weighted version of (11) gives the following COMPLETE tower bounds:

    E_diff<=12400(6-pi^2/2) C eta^2 m0^4,
    E_th<=16 C eta^2 m0^4(2L0+1),
    E_boundary<=114.8 delta^(-k) (eta m0^3)^(1+k)
                         *4Cr/(r-2)*m0^(2-r).       (12)

For the thermal term, integrate -d[E^-4(L0+2log(E/m0))]
against 4CE^2. The resulting coefficient is 4C m0^-2(2L0+1)
before multiplication by 4eta^2 m0^6. Thus its logarithmic angular
tail is included, not bounded by a nonexistent uniform log maximum.

At any finite angular regulator the cap and finite-cut identity are
exact. The positive summable majorants (12) control their omitted pieces
uniformly. This permits the angular/cut limit for THESE errors. It does
not authorize removing the common regulator from the remaining raw mass
action or replacing the original a0/a1/a2 and two-dimensional Weyl terms.
Those cancellations must still occur in the same four-dimensional scheme.

## 5. Quantitative execution and research consequence

Freeze note, checker and review before running. Require

    N U_i U_j (E_diff+E_th+E_boundary)<1e-10

for N<=2 identical exteriors and U_i,U_j<=2, so test the factor 8.
Use the SAME unweighted action units as the previous endpoint execution;
no factor gamma is used to improve the tolerance. This is an endpoint
assembly precision, not the <1.667e-9 metric-error certificate.
Enclose the maximal permitted eta for factor 8 at this target, not only
the selected cut. Verify the actual global polynomial directions from
the finite local assembly have the stated cap norms. Save all pieces,
the mass thresholds, geometry enclosures and failed runs, if any.

The method changes resolution as mass increases and uses the actual
adjoint damping. It neither changes normalization nor optimizes another
cold source-value envelope. It closes the fixed-cap summability obstacle
for the specified endpoint contributions if the executed bound passes.

Decision: continue the physical weak primal-dual route with this adaptive
endpoint contract. The next decisive calculation remains the complete
canonical shape contribution, including common local/curved/boundary
restoration, on the physical forward/bolt-adjoint pair. Use the exact
local rows already assembled and test |Q_app|+||d_b||_V*||d||_V*<1.667e-9,
exterior return <6.67e-10 and Y0<=5e-9 with all X/admission, nonlinear
and exact compatibility requirements. That full test is not executed
by this endpoint certificate. Feasibility is still undecided.

Remaining to-do list: full canonical shape/curved/boundary operator and
physical forward/adjoint residual certificate; X/admission, nonlinear
closure and compatibility; matching, causal realization, semiclassical
validity/noise, all-sector stability, charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
