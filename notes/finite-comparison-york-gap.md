# Uniform global York inverse on the finite comparison path

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary proof and
verification contract before execution. **The actual finite comparison
has a quantitative global York inverse and a controlled physical-L2
TT projector. This is not positivity of its physical TT Hessian.**

The [global TT theorem](stationary-global-tt-reduction.md), section 2,
already constructs the York map at each finite comparison and proves
a positive gap there without a uniform numerical value. The present
result supplies a uniform gap and comparison estimates using the actual
first metric jets. It does not require a second-derivative comparison
of the singular vector potential, or any new metric regularity.

Keep the original path, fixed state/domain, temperature, clock and
flux. Use the notation and bounds of the
[finite-comparison scalar proof](finite-comparison-scalar-resolvent.md):

    sigma_ell=sigma0+ell u, r_ell=B+(r0-B)exp(ell V),
    S=log[1+((t-1)/t)(exp(ell V)-1)], w=ell u+S,
    g_ell=exp(2w)(d rho^2+a0^2 dtheta^2)+r0^2 exp(2S)dOmega^2,
    d rho_ell=exp(w)d rho, y=sqrt(lambda)rho,
    |w|,|S|<epsilon=3e-8,
    |w_y|,|S_y|<1e-3.                              (1)

These bounds hold on the whole physical exterior, including the bolt,
uniformly for 0<=ell<=1. The fields are even at the bolt and retain
the original outer conormal data and correlated mass cubic.

## 1. Physical Hilbert spaces and the exact transformed operator

Let H_V,ell be radial vector components Y in the physical orthonormal
frame, with scalar L2 measure dmu_ell. Let H_T,ell be spherical diagonal
tracefree relative tensors t=diag(t0,t1,t2,t2), with point norm
|t|^2=t0^2+t1^2+2t2^2 and the same physical measure. In common seed
coordinates,

    dmu_ell=exp(2w+2S)dmu0,
    U_V Y=exp(w+S)Y, U_T t=exp(w+S)t.               (2)

The relative orthonormal components in (2) are identified before
multiplication. Both U maps are unitary onto the seed physical Hilbert
spaces. U_V identifies the COMPONENT of a vector in its own physical
proper frame; identifying its coordinate coefficient without that
frame change would give a different map.

Let L_c,ell be the tracefree Lie derivative of the metric and let
P_ell=div_ell L_c,ell on radial vectors, with its regular bolt and
physical L2 realization. The adjoint identity is

    L_c,ell*=-2 div_ell, -P_ell=L_c,ell* L_c,ell/2. (3)

Define the dimensionless first-order map between the seed spaces

    A_ell=lambda^(-1/2) U_T L_c,ell U_V^-1.

Put A0_orb=partial_y log a0 and R0=partial_y log r0. For the seed,

    A_0 v=diag(3v_y/2-(A0_orb+2R0)v/2,
               -v_y/2+(3A0_orb/2-R0)v,
               -v_y/2+(R0-A0_orb/2)v,
               -v_y/2+(R0-A0_orb/2)v).

Direct substitution of d/d(sqrt(lambda)rho_ell)=exp(-w)d/dy and
the two frame logarithmic derivatives gives the EXACT identity

    A_ell v=exp(-w)[A_0 v+E v],
    E=diag(-2w_y-5S_y/2, 2w_y-S_y/2, 3S_y/2, 3S_y/2),
    |E|^2=8w_y^2+8w_y S_y+11S_y^2.                (4)

It includes the derivative of the unitary density in (2). In
particular there is no residual w/y or S/y error at the bolt.
As a separate check, for a COMMON coordinate vector with seed proper
component v/sqrt(lambda), its candidate proper component is
exp(w)v/sqrt(lambda). Before the density transformation its tracefree
Lie derivative is `A_0 v+(w_y-S_y)v diag(1,1,-1,-1)`.
Here w_y-S_y=ell u_y: the finite radius completion cancels exactly.

## 2. Uniform form gap and the regular domain

The explicit seed projection potential in the global TT theorem,
sections 2-3, gives

    ||A_0 v||^2/2
      >=(3/2)||v_y||^2+(149/100)||v||^2.            (5)

All these norms use dmu0. Set delta=1/250=.004. From (1),(4),(5),

    ||(A_ell-A_0)v||
      <=[exp(epsilon)-1
           +exp(epsilon)sqrt(27)*1e-3/sqrt(2.98)]||A_0 v||
      <delta ||A_0 v||.                            (6)

The bound is strictly uniform in ell. Thus

    (1-delta)||A_0 v|| <= ||A_ell v||
                         <=(1+delta)||A_0 v||,
    ||A_ell v||^2/2 >1.47||v||^2,                 (7)

for nonzero v. The second inequality follows from
`1.49(1-.004)^2>1.47`.

Start on smooth radial vectors compact away from infinity, with the
regular vector behavior Y=O(rho) at the bolt. The vector graph core
can also be approximated by punctured-bolt smooth vectors: the
O(rho) branch has vanishing cutoff cost. The excluded rho^-1 branch
has infinite form norm. Equations (2),(4) preserve the regular even/odd
parities; (6) proves equality of the two closed graph domains, not
only an estimate on formal expressions. The multiplier exp(-w) and
its inverse are bounded, and E is bounded. Thus the conjugated
Friedrichs operator is exactly

    Q_ell=U_V(-P_ell)U_V^-1=(lambda/2) A_ell* A_ell.

Consequently the actual finite-comparison York operator has

    -P_ell >=1.47 lambda,
    ||(-P_ell)^-1||_(H_V,ell -> H_V,ell)<.681/lambda. (8)

The seed and finite vector realizations have the same selected bolt
condition and no free homogeneous L2 vector. No division by A-R or
photon-sphere coefficient occurs. The actual P_ell potential may
contain second derivatives, but its form has been compared through
the exact first-order factorization, including the physical measure.

For later nonlinear neighborhoods, the same argument has an explicit
sufficient criterion. If |w|<=eps, |w_y|<=a, |S_y|<=b, set

    d=exp(eps)-1+exp(eps)sqrt(8a^2+8ab+11b^2)/sqrt(2.98).

Whenever d<1, the gap is at least `1.49(1-d)^2 lambda` on the same
regular domain. This is a conditional C1 neighborhood criterion;
membership of an unknown nonlinear correction has not been proved.

## 3. Uniform TT projection and comparison with the seed

The closed ranges Ran A_ell and Ran A_0 are closed by (5)-(7).
Let R_ell denote orthogonal projection onto Ran A_ell in H_T,0,
and Pi_TT,ell=I-R_ell. This is the unitarily identified actual
orthogonal projection onto ker div_ell. It is a kinematic tensor
projection, distinct from the fixed source-compatible Pi_cut of S167.

For x=A_0 v, distance to Ran A_ell is at most delta||x||. For
x=A_ell v, distance to Ran A_0 is at most
delta/(1-delta)||x||. Therefore

    ||(I-R_ell)R_0||<=delta,
    ||(I-R_0)R_ell||<=delta/(1-delta).

To obtain a full operator bound, use
`R_ell-R_0=R_ell(I-R_0)-(I-R_ell)R_0`.
The outputs are orthogonal, and the two input components from R_0
and I-R_0 are orthogonal. Apply the two directional estimates
(and the adjoint of the second one) to these components. It follows
that

    ||Pi_TT,ell-Pi_TT,0||=||R_ell-R_0||
                         <=delta/(1-delta)<.0041.  (9)

Each TT projector itself has norm at most one in its actual physical
tensor L2 space. Equation (9) compares them AFTER the unitaries (2).
It is not a bound for a projection in a Holder, H4-log or conormal norm.

For any t in H_T,0, the weak York extraction is

    v_ell=lambda^(-1/2)(A_ell* A_ell)^-1 A_ell* t,
    sqrt(lambda) A_ell v_ell=R_ell t.              (10)

The expression in (10) is the bounded left pseudoinverse, defined
weakly for every L2 tensor; it does not assume div t is L2. Returning
Y_ell=U_V^-1 v_ell solves `P_ell Y_ell=div_ell(U_T^-1 t)`.
Its quantitative value and graph bounds are

    sqrt(lambda)||v_ell||<.584||t||,
    sqrt(lambda)||partial_y v_ell||<.58||t||,
    lambda[3||partial_y v_ell||^2+2.98||v_ell||^2]
                                      <=(1-delta)^-2||t||^2. (11)

Strict inequalities in (11) concern nonzero t; extraction of zero is
zero. The derivative is the seed proper-coordinate derivative of the
UNITARILY transformed vector component, not an unidentified full
tensor-gradient norm.

For the same identified tensor t, the seed extraction v_0 satisfies

    sqrt(lambda)||A_0(v_ell-v_0)||<.0081||t||,
    sqrt(lambda)||v_ell-v_0||<.0047||t||.           (12)

Indeed A_ell applied to the difference is lambda^-1/2 times the
projector difference minus (A_ell-A_0)v_0. Equations (6),(7),(9)
give the exact coefficient
`delta(2-delta)/(1-delta)^2<.0081`, followed by (5).

The vector inverse itself has a useful form comparison. Define the
bounded positive operator B_ell by the form
`Q_0^-1/2 Q_ell Q_0^-1/2`; its spectrum is within
[(1-delta)^2,(1+delta)^2]. Its inverse differs from identity by less
than .0082. The form resolvent identity therefore yields

    ||Q_ell^-1-Q_0^-1||<.0056/lambda.             (13)

This compares the two York inverses in the fixed seed vector space.
It is not a perturbation estimate for the physical Hessian inverse.

## 4. Consequence and limits

The formerly qualitative finite-comparison York gap now has a uniform
number, an L2 tensor projector comparison, and controlled value/graph
extraction. These remove a low-order kinematic obstacle to transferring
the TT analysis. They include the bolt, entire cold region and outer
end. They do not truncate angular quantum channels or change any
physical source, state or preparation parameter.

The [curved quantum TT bound](stationary-curved-quantum-tt.md) remains
a SEED theorem. Smallness of (9) in L2 cannot transfer its positive
fourth-order logarithmic form. The actual finite-comparison TT
Hessian, mixed vertices, density and off-shell return contacts still
need quantitative estimates in their own form/domain norms. The
candidate's sharp logarithmic York map is also not proved by (11).

Pi_cut remains fixed at the seed, with all preceding source/correction
bounds unchanged. No norm of its cutoff direction, cold correction
bound, candidate coupled range, nonlinear family preservation,
projected solution or actual signed compatibility value is supplied.
The scalar reference resolvent and this vector York inverse are
components of the prospective coupled construction, not a full solution.

## Frozen verification contract

Freeze the note and checker before execution. Reconcile the current
finite-comparison, global TT and curved TT evidence/input identities.
Verify the physical measure, relative component transformation,
first-order factorization, exact E polynomial and all rational bounds.
Check projection and extraction identities using independent rectangular
matrix fixtures, including a change of physical inner product. Those
fixtures are not evaluations of the physical quantum Hessian.
The graph closures, Hilbert-space range estimates and weak inverse
construction remain ordinary proofs. Retain every terminal run and
failed source if any. No Lean/cloud campaign is selected.

Remaining to-do list: finite-comparison TT/mixed/contact forms and sharp
logarithmic mapping, candidate coupled range; fixed projection direction
and cold/higher inverse/source bounds; nonlinear preservation and
projected full solve; signed compatibility; exterior/state, charge and
return; opening, stability, transport, preparation, separate creation
and full-result audit.
