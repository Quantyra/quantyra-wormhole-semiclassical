# Coupled asymptotic obstruction and the compatibility derivative

2026-10-17, S139/S074. Ordinary argument before focused verification.

The [fast stationary phase](stationary-fast-phase-variation.md) has a
differentiable action. This does not imply that its equations have a
bounded inverse. Here the earlier scalar obstruction is extended to
the FULL tensor linearization on a specified unweighted Sobolev scale.
The scalar compatibility condition is also differentiated without
treating a scalar horizon amplitude as an independent metric datum.

The result rejects unweighted higher-order inversion as a shortcut.
It does not exclude the required weighted boundary problem, prove
nonexistence of an exterior, or establish a physical instability.
The oscillations involved have proper Planck momentum; no reliability
of the retained effective theory at that momentum is inferred.

## Equations, hypotheses and scope

Keep the full selected massless Dirac functional, common subtraction,
fixed magnetic q, physical G, mu0 and lambda=lambda_e. Write

```
E_ab=G_ab+Lambda0 g_ab-8pi G(T_Maxwell,ab+T_Dirac,ab),
Z=R4+12lambda, epsilon=G/(60pi),
a=1+11G lambda/(180pi), P_g=a+epsilon Box_g,
E_trace=P_g Z-S[g]=-g^ab E_ab,
S[g]=G/(2pi)[(2/3)F_em^2-(7/360)C4^2
                       -(11/180)J^2+(11/2160)Z^2].  (1)
```

These are the [exact retained trace conventions](stationary-trace-range-condition.md).
The bookkeeping completion parameter is not lambda. The fixed-flux
magnetic invariant in areal gauge is `F_em^2=q^2/(2r^4)`.

For the no-estimate theorem, suppose an exact stationary solution of
E_ab=0 exists with the selected HH prescription and paired domain,
a regular nonextremal bolt, fixed prepared temperature, and

```
F=lambda r^2+1-2m/r+O(r^-2 log^p r),
alpha=O(r^-4 log^p r),                             (2)
```

with five differentiated weighted bounds. Assume the solution is smooth
on the OPEN exterior. That extra local smoothness ensures each compact
test below is a legitimate smooth core variation of the determinant;
no sixth ASYMPTOTIC derivative bound is used. This note does not prove
such a solution or its regularity. The exact-solution assumption makes
the diffeomorphism variation of E vanish and supplies the required
outer stress limit from the equations themselves.

Use static spherical tensors and positive Euclidean component norms on
the regular Wick-rotated section. Divide out the constant time and
sphere volumes. The scalar measure is `w dr=exp(alpha)r^2 dr`.
Hj uses covariant derivatives through j in this positive metric.
There is no extra radial weight penalizing escape to infinity.
Let

```
C_g(h)_b=nabla^a h_ab-(1/2)nabla_b tr_g h            (3)
```

be de Donder gauge. The assertion is that no finite constant K can make

```
||h||_H4 <= K ( ||DE_g[h]||_L2 + ||C_g(h)||_H3 )    (4)
```

hold on all compact smooth static spherical variations. Adding a norm
of h on any FIXED compact interior region to the right does not rescue
(4). We do not assume the complete nonlocal logarithmic response is a
bounded operator H4 -> L2 on its entire domain. Failure of the necessary
estimate (4) on this core is sufficient for the stated conclusion.
No claim of a Fredholm index for the unknown weighted problem is made.

## Why the complete conformal response is local

Let sigma be compactly supported in the exterior and vary
`g -> exp(2sigma)g`, with the magnetic potential fixed. The massless
Dirac conformal transformation transports its Green function, state and
paired domain. Because sigma is static and vanishes near both endpoints,
the time period, regular bolt and boundary projector are unchanged.
The conformally transported HH family is the prescribed family here;
there is no independent variation of its state.

The exact finite Weyl identity is
`Gamma_D[exp(2sigma)g]=Gamma_D[g]+W[sigma;g]`.
Differentiate once in the base inverse metric. If
`T^W_ab=-2/sqrt(g) delta W/delta g^ab`, it gives

```
T_D,ab[exp(2sigma)g]=exp(-2sigma)(T_D,ab[g]+T^W_ab),
delta_sigma T_D,ab=-2sigma T_D,ab[g]+A_ab[sigma;g],
A_ab=-2/sqrt(g) delta_g^ab integral sqrt(g) sigma A4[g]. (5)
```

The last expression is a LOCAL operator of order at most four in
sigma, with coefficients made from curvature and its derivatives.
All boundary variations vanish for these compact tests. The independent
Weyl-invariant functional in Gamma_D contributes its base-stress
rescaling in (5); it is not discarded as a tensor response in general.
Maxwell stress scales in the same way without an anomaly term.
Thus the otherwise nonlocal tensor response is local on this special
conformal channel. The stress at a hypothetical exact solution tends
to the selected vacuum stress by E=0 and (2), so the coefficients of
this local channel tend to their AdS values.

For comparison, [Anderson--Molina-Paris--Mottola, sections III-IV,
equations (28)-(29) and (53)-(63)](https://arxiv.org/html/0907.0823v2)
separate local conformal response, Weyl-dependent nonlocal response
and independent state variations. Their conformally flat cosmological
calculation is not an RN/HH boundary-value theorem. Equation (5) here
uses the exact Weyl identity in the retained prescription.

## Full AdS tensor factorization

On the selected on-shell maximally symmetric vacuum, Ric=-3lambda g,
C4=J=Z=F_em=0. For h_ab=2sigma g_ab,

```
delta R4=-6(Box-4lambda)sigma,
D_ab f=nabla_a nabla_b f-g_ab Box f+3lambda g_ab f,
nabla^a D_ab f=0,
tr D f=-3(Box-4lambda)f.                           (6)
```

The divergence identity uses Ric=-3lambda g and the commutator on a
gradient. By (5) the conformal tensor response is local, of order at
most four. Maximal symmetry permits only the Hessian and metric
operators applied to scalar polynomials in Box. Conservation reduces
this general local scalar tensor to `D_ab(A+B Box)sigma`.
There is no independent nonlocal state projector under the stated
prescription. On shell the trace of DE is minus the variation of
E_trace. Since delta S=0 and Z=0 at this vacuum, (1),(6) give
`tr DE=6(Box-4lambda)(a+epsilon Box)sigma`.
Matching the two polynomial coefficients therefore proves

```
DE_AdS[2sigma g]=-2D_ab(a+epsilon Box)sigma.         (7)
```

This is a full tensor identity, not just its trace. In particular,
P sigma=0 solves every vacuum tensor component and
`delta R4=6(a/epsilon+4lambda)sigma` is nonzero for sigma nonzero.
The branch is not removed by the tensor constraints or a pure gauge.

In spherical AdS, restoring areal gauge with xi^r=-r sigma gives

```
delta F=2rF sigma'-rF' sigma,
delta alpha=sigma-r sigma'.                        (8)
```

For `sigma~r^s`,
`s=-3/2 +/- i omega`, `omega^2=a/(epsilon lambda)-9/4`,
the leading metric terms are `2lambda(s-1)r^(s+2)` and
`(1-s)r^s`. Their relative metric amplitudes tend to zero, but the
`r^(1/2)` oscillation of delta F is excluded by (2)'s finite-mass
falloff. A nonzero fixed monopole and a nonextremal round bolt are
not globally maximally symmetric AdS. Thus (7)-(8) describe the
vacuum bulk/outer limiting operator, not a constructed global RN mode.

## Escaping witnesses for the full gauge-fixed equations

The local response (5) converges coefficient by coefficient to (7)
on the end of any exact solution satisfying the hypotheses. Curvature
and its required derivatives converge by (2); the base stress does so
by E=0. Consequently the difference, acting on compact conformal
tests, is a local operator through order four whose coefficients
tend to zero in a physical orthonormal frame. The nonlocal shape
response does not need to be estimated for these tests.

Use a FIXED coordinate and amplitude, not a differentiated moving
proper coordinate:

```
x=lambda^(-1/2)log(r/B), mu_infinity=sqrt(lambda)r^3,
k_*^2=a/epsilon-9lambda/4=60pi/G+17lambda/12>0,
sigma_T=mu_infinity^(-1/2) T^(-1/2)
          chi((x-R_T)/T) exp(i k_* x),              (9)
```

where chi is smooth, supported in (0,1), with integral |chi|^2=1,
and R_T/T -> infinity. Fixed length units are understood in the norm
normalization. The actual measure is
`w dr=exp(alpha)mu_infinity dx`, so `||sigma_T||_L2 -> 1`.
Every support is compact and eventually leaves any fixed interior set.

The scalar limiting operator is
`P_infinity=a+epsilon(partial_x^2+3sqrt(lambda)partial_x)`.
Conjugation by mu_infinity^(-1/2) gives
`epsilon(partial_x^2+k_*^2)`. Hence exactly

```
P_infinity sigma_T=epsilon mu_infinity^(-1/2)T^(-1/2)e^(ik_*x)
            [2ik_* chi'/T+chi''/T^2].              (10)
```

The actual coefficients differ by decaying conormal terms. Their
derivatives through the orders below also decay. It follows that

```
||P_g sigma_T||_H4 -> 0,
||sigma_T||_H6 <= C,                               (11)
```

with C depending on the FIXED physical k_* and geometry, not T.
There is no uniform bound as G tends to zero being asserted.
To check the regularity count, Box on a scalar in the fixed x chart
has a principal coefficient involving F and a first-order coefficient
involving F',alpha'. Four derivatives of P therefore use at most
five metric derivatives. Covariant derivatives of sigma through six
likewise use at most five metric derivatives. The smooth compact
core is justified by open-exterior smoothness, not a new sixth-order
asymptotic estimate.

Equations (5),(7),(11) imply
`||DE_g[2sigma_T g]||_L2 -> 0`. Include the gauge correction

```
c=-2epsilon/(a+6epsilon lambda),
xi_T=c grad sigma_T,
h_T=2sigma_T g+L_xi_T g.                           (12)
```

It is compact, static and spherical. For a radial xi the fixed
magnetic potential has L_xi A=0. Covariance at the exact solution
therefore gives `DE_g[L_xi_T g]=L_xi_T E=0`.
The identities

```
C_g(2sigma g)=-2 grad sigma,
C_g(L_xi g)=(Box+Ric)xi,
C_g(h_T)=-2/(a+6epsilon lambda) grad(P_g sigma_T)
                +2c(Ric+3lambda g) grad sigma_T    (13)
```

show `||C_g(h_T)||_H3 -> 0`; derivatives of the last curvature
coefficient through three require precisely the five metric jets.
Meanwhile

```
tr h_T=[8+4a/(a+6epsilon lambda)]sigma_T
                    -4/(a+6epsilon lambda)P_g sigma_T. (14)
```

In the positive tensor norm `||tr h||<=2||h||`. Thus h_T has L2
norm bounded below by a positive constant and H4 norm bounded above.
The full equation and gauge residuals tend to zero. Normalizing h_T
in H4 proves failure of (4), including any fixed compact-interior
correction. Complex notation in (9) is only convenient: an estimate
for all real tensors would imply the same estimate for their real
and imaginary parts together, so the complex witness also contradicts
it. No complex physical metric is assumed.

Each witness satisfies the endpoint conditions because it is compact.
Its norm in a space with additional FAST radial weights need not be
bounded. The commutator with such a weight can also change the residual
estimate. Therefore this proof does not settle the selected weighted
inverse, its kernel/cokernel or its horizon matching. It specifically
closes the loophole that the full tensor equations might restore an
ordinary unweighted estimate absent in the scalar trace.

## Compatibility is a coupled projection, not an extra free amplitude

For any regular fast comparison define phi_g by P_g phi_g=0,
phi_g(r_h)=1 and regularity, and write
`C[g]=integral_(r_h)^infinity w phi_g S[g] dr`.
The accepted exact identity is
`C[g]=-integral w phi_g E_trace[g] dr`.
Differentiate a smooth fixed-source family in areal coordinates,
allowing its horizon to move. The full formula is

```
DC[h]=-integral w phi DE_trace[h] dr
      -integral delta(w phi) E_trace dr
      +(w phi E_trace)_h delta r_h,
DE_trace[h]=P_g delta Z+epsilon(delta Box_g)Z-delta S. (15)
```

For a fast tangent delta Z, the Green concomitant in the integral of
phi P_g delta Z vanishes at both ends. At a GENUINE solution only,
E_trace=0, and therefore

```
DC[h]=integral w phi[delta S-epsilon(delta Box_g)Z] dr,
delta S=G/(2pi)[-7 delta(C4^2)/360-11 delta(J^2)/180
                                      +11Z delta R4/1080]. (16)
```

Here delta(F_em^2)=0 at fixed q and fixed areal r; it is not an
adjustable forcing. The variation of each curvature square includes
its metric contractions, not merely an unraised tensor variation.
The delta-Box term is retained at the nonconstant centered curvature.
At the maximally symmetric vacuum all terms in (16) vanish. Thus the
leading vacuum calculation supplies no nonzero compatibility Jacobian.
That statement is about its vacuum limiting problem; it is not a
claim that C has zero derivative on the actual fixed-q exterior.

At the ACCEPTED CANDIDATE the defect is nonzero. The second integral
and moving-endpoint term of (15) cannot be omitted. In particular,
delta phi solves `P_g delta phi=-epsilon(delta Box_g)phi` with
the varied regular endpoint normalization. The existing mass/lapse
defect bounds give no uniform bound for this oscillatory test-scalar
response. Substituting a candidate into (16) as though it were an
exact solution would not justify a shooting Jacobian.

Finally the spherical Ward identity is

```
E_r' +(alpha'+F'/(2F))(E_r-E_t)+(2/r)(E_r-E_theta)=0. (17)
```

At an exact solution its linearization has the same coefficient form.
Solving the two independent orbit equations already forces the angular
equation. The scalar moment is a necessary range condition in a
trace-plus-constraint formulation, not a third independent equation
to append after solving all tensor equations. Counting only the
auxiliary scalar's horizon amplitude and its two forbidden outer tails
therefore proves no full coupled nonexistence theorem.

## Research decision and source comparison

Do not commission a scalar-only shooting calculation or assume a new
horizon datum has a nonzero compatibility derivative. The next full
solver must construct the weighted coupled boundary map, identify its
independent regular horizon data after constraints and fixed temperature,
and prove its range/transversality on the actual geometry. Equations
(15)-(17) state the required projection and prevent double counting.
A controlled spatial reduction remains an alternative only with its
covariant definition and observable-error/branch prescription.

[Lu--Perkins--Pope--Stelle, sections 6-7](https://arxiv.org/html/1508.00010v2)
provide a useful comparison for local curvature-squared gravity. Their
trace no-hair integral uses a positive scalar mass squared and
asymptotically flat boundary data; after that reduction an additional
horizon parameter appears in the traceless problem. Those assumptions
do not transfer to the present sign, AdS end and nonlocal Dirac theory.
Neither their parameter count nor their no-hair result supplies our
missing coupled map.

## Focused verification contract

The complete ordinary obstruction and projection argument precedes
execution. A focused checker will test the exact AdS radial tensor
divergence/trace, conformal areal-gauge curvature, de Donder correction,
limiting conjugation/cutoff identity, compatibility variation signs and
Ward algebra. All symbolic residuals must be exactly zero; rational
inequalities must hold. No numerical tolerance or Planck-frequency
shooting is selected. Conformal-domain transport, asymptotic coefficient
convergence and the Sobolev sequence estimates remain ordinary proofs.
Hash exact UTF8/LF inputs before/after, retain every execution and failed
input snapshot, and record the actual Python -B command and versions.
No predecessor checker or Lean campaign is rerun.

Remaining to-do list: verify this coupled unweighted obstruction;
construct the genuinely weighted coupled range map at fixed preparation
or justify controlled spatial reduction; select the physical boundary
branch and complete exterior/state matching, S138 transfer and S074-S080.
