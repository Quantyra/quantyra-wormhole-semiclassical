# Weighted propagation on the paired domain: square-domain correction and limit

2026-10-28 HST, E003/S169. Ordinary argument before numerical execution.
The preceding trial/regularity calculation is accepted supporting progress,
not a complete quantum certificate. Preserve P32, the closed audit, original
physical norms and the full stable semiclassical wormhole objective.
This attempts the proposed localized propagation. The complete differentiated
heat remainder remains the target; an undifferentiated operator estimate
must not be substituted for its diagonal or common-regulator limit.

## Spaces and the first unsupported domain step

Let Hspin be the original L2 spinor space of the finite orbit geometry,
with the physical paired copies, thermal spin structure and regular bolt.
Use the accepted self-adjoint A=iJD0, L=A^2 on its composed domain. This
is NOT the metric energy space E_t. The original paired square has
complementary Dirichlet/Neumann data at the geodesic outer seam in the
fixed half-density frame. Keep both original field equations and their
York reconstruction in the source; no gravitational inverse is applied here.

For the accepted localized trial introduce the computational weight

    w=gamma/(gamma+Lbolt-z), W=w^3, f=log W.

W and W^-1 are bounded on this fixed finite geometry and identical on
both copies. They preserve the FIRST-ORDER domain by the weak product
rule and the scalar boundary projector, as in the accepted rough-domain
argument. They do not in general preserve Dom(A^2). At the seam, with
inward proper coordinate y,

    (partial_y log W)(0)=3 sqrt(beta)/(gamma+Lbolt)>0.  (1)

A smooth core field in the Neumann sector with nonzero boundary trace
and zero normal derivative is admissible for L, but W times that field
has normal derivative W' times its nonzero trace. A cutoff away from
the bolt retains the original regular-bolt condition. This is a domain
test, not a finite-mode approximation of the quantum trace.

Thus the proposed shortcut W^-1 L W on the UNCHANGED square domain is
undefined as that similarity. Its actual domain is W^-1 Dom L, with
Neumann sector condition partial_y u+(partial_y log W)u=0. The numerical
size of (1) is reported; the tiny value W(0) does not make this relative
Robin term disappear. This changes no physical boundary condition: it
is the exact computational similarity of the original operator.

## Correct form construction and complete propagation bound

Use the fixed form domain Dom A, which the scalar multipliers preserve.
For C_f=[A,f], ||C_f||<=alpha=||df||infinity, the exact form of
L_f=W^-1 L W is

    q_f(u,v)=<A W u,A W^-1 v>
            =<(A+C_f)u,(A-C_f)v>.                    (2)

In particular Re q_f(u,u)>=||Au||^2-alpha^2||u||^2 and
|Im q_f(u,u)|<=2alpha||Au|| ||u||. This includes the transformed boundary
condition through the form domain; no integration-by-parts boundary term
is discarded. The similarity is defined already because W is bounded and
invertible. The following energy calculation gives a bound independent
of its large condition number.

For zeta=x+iy with x>0, completion of the square gives

    Re[zeta q_f(u,u)]>=-alpha^2 |zeta|^2/x ||u||^2.

Apply the differential inequality to W^-1 exp(-s zeta L)W, first on
its generator domain and then by density. Gronwall yields

 ||W^-1 exp(-zeta L)W||Hspin->Hspin
                  <=exp(alpha^2 |zeta|^2/Re zeta).    (3)

This is an ordinary proof, not an import of a scalar positivity theorem
to the paired matrix operator. The original common angular mass squared
m_n=n(n+qang)/C commutes with the similarity. For every n and all thermal
and radial modes, (3) gains exp(-m_n Re zeta). No sector is omitted.

On |zeta-t|=t/2 one has |zeta|^2/Re zeta<=3t/2 and Re zeta>=t/2.
Cauchy's derivative formula therefore gives, for D_t=t partial_t,

 ||W^-1 exp[-t(L+m_n)] W||<=exp[(alpha^2-m_n)t],
 ||W^-1 D_t exp[-t(L+m_n)] W||
                     <=2exp(3alpha^2 t/2-m_n t/2).   (4)

The adjoint is the inverse-weight conjugation, so the identical bounds
hold with W and W^-1 exchanged. For products whose positive time lengths
sum to t, the undifferentiated factors in (4) multiply to exp(alpha^2 t),
not an independent penalty for each subdivision. These are Hspin operator
norms in both directions, NOT pointwise, Schur or diagonal-density bounds.
No norm-to-diagonal implication is asserted.

## Actual geometric constants, without a new physical normalization

Write delta=1-Lbolt=h/(1+h). The accepted exact finite geometry obeys

 Q=(1-z)^2 Pseed(z)-tau z^3,
 Pseed=beta+2beta z+(3beta+u)z^2+epsilon_seed z^4.

Pseed'>=2beta-4|epsilon_seed|>0 and Pseed(1)=1, so
0<=Q<=(1-z)^2 on the original interval. Consequently

 |df|=3sqrt(Q)/(gamma+Lbolt-z)
     <=3 max(1,delta/gamma)=alpha_upper.              (5)

This is finite at the bolt in its regular chart. The weight is an analysis
multiplier, not a small admitted metric perturbation; no small-b theorem
is applied to its large logarithmic value. All parameters retain their
accepted correlated definitions and outward enclosures.

The exact original-volume norm is

 ||W||L2(dz)^2=gamma(1-w0^5)/5,
 w0=gamma/(gamma+Lbolt).                              (6)

Compute (1), (3)-(6) at t0=.001. Also display the crude similarity bound
||W||infinity ||W^-1||infinity=1/W(0) in the SAME Hspin norm. Its large
value is not compared directly with a physical-output ceiling. Correct
form control, not a normalization change, avoids that particular estimate.

For sensitivity only, a COMPLETE pointwise quantum remainder bound
|r_v(t,z)|<=D_W W(z)t^3 would give D_v<=D_W ||W||2. Thus its conditional
D_W ceiling follows from the preceding actual-trial D_v ceiling divided
by (6), with ALL OTHER unknown errors zero. Neither D_W nor such a complete
kernel estimate is attained here. It includes more than propagation (4).
In particular the Hspin constants cannot be inserted in place of D_W.

## What the full attempt did not establish

The differentiated parametrix defect is an operator kernel, not a scalar
source density. Its first correction is P*F and the quantum observation
restricts it to coincident endpoints after metric/time differentiation.
Equations (3)-(4) neither bound the diagonal trace of that kernel nor
its complete subtracted t^(11/4)/t^3 remainder. The necessary additional
estimate must retain localization through the varied transport amplitudes,
the coincidence vertex, the paired boundary part and both weighted Schur
directions before restriction to the diagonal. The corresponding weighted
kernel/trace constants remain UNKNOWN. Summing the bare angular operator
bounds or taking norms before common subtraction does not establish it.
The exact ultraviolet logarithm and all finite returns remain obligations.

The source review rechecks [Smith, theorem 1.1](https://arxiv.org/pdf/1809.05614)
and its trace-class argument in section 2. It treats scalar potential
perturbations of a fixed Laplacian on a complete smooth manifold. It does
not furnish the paired metric-Hessian remainder or its numerical constant.
Its explicit distinction between operator and trace norms supports keeping
that gap open. This assessment does not adopt a new theorem for our domain.

Decision: reject the frozen-square-domain shortcut and use the corrected
form similarity if continuing weighted propagation. Keep the construction.
The modest bounds (4) alone supply no new coupled margin. The specified
alternative is a complete localized resolvent column evaluator on the
original paired first-order domain; it avoids this heat-diagonal propagation
step but still needs the commonly restored, all-angular trace remainder.
The existing quartic evaluator does not cover this localized trial. Neither
method currently has an attained complete column error.

The ONE next calculation remains the full weighted differentiated kernel
remainder for the explicit trial, including its diagonal observation,
boundary and numerical errors in E_t*. This is still PROPOSED. Domain
correction and whole-spectrum Hspin propagation are EXECUTED supporting
results; the requested complete quantum-feasibility milestone is UNRESOLVED.
No physical budget or existence/stability requirement is promoted.

Remaining to-do list: complete signed quantum column/first residual and
later applications; G_t and every coupled/source/graph/physical/boundary/mass
budget; nonlinear admission/contraction and existence; exact compatibility;
physical matching; semiclassical validity/noise; all-sector stability;
the full stable semiclassical wormhole objective.
