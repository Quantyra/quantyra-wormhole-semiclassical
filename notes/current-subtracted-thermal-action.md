# Exact action subtraction and a complete discrete thermal sum

2026-10-23, E003/S169. Ordinary argument before execution.
The full stable semiclassical wormhole objective remains incomplete.

The closed audit selects the weak physical H e=-r1 calculation. Its
current executable action is regulated: the raw thermal Hessian diverges.
This implements subtraction INSIDE that action, followed by a complete
discrete thermal sum and a certified infinite tail, for each of the first
two angular levels. This is reusable weak-action assembly, not another
isolated frequency sample, a full angular sum, or the physical correction.

## Exact subtraction with the paired domain retained

Use the actual finite seed, gamma, M_n and compact C5 directions of
[the coincident action](current-coincident-dirac-action-assembly.md).
Let f=sqrt(Q), W=p/f, p=gamma(j+1/2), m=M exp(b), and prime denote
canonical proper differentiation, dz/d rho=-f. The regular solution is
x'=m x+W y, y'=W x-m y, with U=y/x in [0,1]. Write

    K=sqrt(m^2+W^2), theta=atan(W/m), alpha=atan(U),
    X=tan(alpha-theta/2), c=theta'/2.

Since alpha'=K sin(theta-2alpha), the exact transformed equation is

    X'=-2K X-c(1+X^2),
    c=-p(mf)'/[2(p^2+(mf)^2)].                         (1)

The physical paired action for BOTH frequency signs is -2 log sqrt(x^2+y^2)
at the seam, up to a mass-direction-independent inner constant. Its radial
derivative is -2K cos(2(alpha-theta/2)). Subtract the local flat-frequency
action -2 integral K d rho BEFORE differentiating or interval estimation:

    A_diff=4 integral K X^2/(1+X^2) d rho.              (2)

This exact difference has the same paired-endpoint variation as the raw
characteristic. There is no second independently prescribed outer datum.
All integrations extend from z=15/16 through the seam z=0, including the
propagated sensitivities AFTER the directions vanish. The physical inner
solution is enclosed by X in [-1,1]. The directions vanish before z=7/8,
so initial directional sensitivities are zero. No wall replaces the bolt.

For a direction u let Y=D_u X, Z=D_u^2 X, K_i=D_u^i K, c_i=D_u^i c.
Differentiating (1) gives

    a=2K+2cX,
    Y'=-aY-2K_1 X-c_1(1+X^2),
    Z'=-aZ-4K_1Y-2K_2X-c_2(1+X^2)-4c_1XY-2cY^2.    (3)

At the background put D=p^2+M^2 f^2 and l=M^2 f^2/D. Then

    K_1=K l u, K_2=K(2l-l^2)u^2,
    c_1=-pM/(2D) [f u'+f' u(1-2l)],
    c_2=-pM/(2D) [(2-4l)f u u'+(1-8l+8l^2)f'u^2]. (4)

With g=X^2/(1+X^2), g'=2X/(1+X^2)^2 and
g''=2(1-3X^2)/(1+X^2)^3, the first and second variations of (2) are

    4 integral [K_1 g+K g'Y],
    4 integral [K_2 g+2K_1 g'Y+K(g''Y^2+g'Z)].       (5)

Evaluate u1=[4w(1-w)]^6, u2=w u1, w=(z-.5)/.375, and u1+u2
on [.5,.875], zero elsewhere. Polarization gives H12=(Hsum-H11-H22)/2.
These functions serve the weak-action implementation; they are neither
the actual h1 nor the unknown correction. The method accepts other
direction/derivative enclosures with the corresponding tail proof.

The next replay removes the second sensitivity Z from numerical integration
by an exact weak identity, avoiding its extra interval loss. Define an
AUXILIARY ACTION-INTEGRATION adjoint ell with ell(seam)=0 and

    ell'=a ell-4K g'.

It obeys a positively damped equation when integrated backward. This is
not the gravitational bolt adjoint psi_b and has no gravitational delta
source. Since Z(inner)=0, integration by parts gives
integral 4K g' Z=integral ell[Z'+a Z]; both endpoint terms vanish.
Consequently, with K_i,c_i for each of the two directions and mixed K_ij,c_ij,

    grad_i=integral [(4g-2ell X)K_i-ell c_i(1+X^2)],
    Hess_ij=integral [(4g-2ell X)K_ij-ell c_ij(1+X^2)
       +(4g'-2ell)(K_i Y_j+K_j Y_i)
       -2ell X(c_i Y_j+c_j Y_i)+(4K g''-2ell c)Y_i Y_j]. (5a)

Here K_ij=K(2l-l^2)u_i u_j and
c_ij=-pM/(2D)[(1-2l)f(u_i u_j'+u_j u_i')
 +(1-8l+8l^2)f'u_i u_j]. The forward whole-cell X,Y enclosures and
backward whole-cell ell enclosures give a validated integral. Direct
mixed entries avoid the unnecessary width from numerical polarization.
The analytic tail still follows from the equivalent original formula (5).
The action-adjoint integration starts at the first nonzero direction
cell: before it both first and second sensitivities vanish identically,
so that prefix contributes exactly zero to every variation. Forward X
still propagates from the original inner cut through the whole prefix.
The square X^2 is explicitly intersected with its nonnegative range
before forming 1+X^2; a wide symmetric interval must not invent a pole.

## Continuous coefficients, signed enclosures

For v'=-a(rho)v+F(rho), a in [amin,amax]>0 and cell length in [Lmin,Lmax],
variation of constants encloses the endpoint by E I+Phi F, where

    E=[exp(-amax Lmax),exp(-amin Lmin)],
    Phi=[(1-exp(-amax Lmin))/amax,(1-exp(-amin Lmax))/amin].

For the whole cell use E=[exp(-amax Lmax),1] and
Phi=[0,(1-exp(-amin Lmax))/amin]. All products are SIGNED interval
products. These bounds allow variable a,F; they do not freeze coefficients.
The first numerical attempt used these independent E/Phi bounds and missed
the accuracy target. The replay uses a sharper signed scalar comparison.
Let q=F/a in [qlo,qhi]. For the lower barrier starting at Ilo, approach
qlo with amin if Ilo<qlo and amax if Ilo>qlo; for the upper barrier
starting at Ihi, approach qhi with amax if Ihi<qhi and amin otherwise.
These choices bound a(t)(q(t)-v) on either side of the equilibrium and
never cross it. Their endpoint bounds are the lower/upper endpoints of
qlo+(Ilo-qlo)E and qhi+(Ihi-qhi)E, respectively, using the same E interval.
Each barrier is monotone, so its initial/end hull encloses the whole cell.
This works for signed v and retains the damping/forcing cancellation.
First enclose X with damping 2K+cX using the invariant [-1,1], refine
with the resulting whole-cell enclosure, then bound Y and Z in triangular
order from (3). The inherited proper mesh encloses f,f',L,u,u' everywhere.
The known invariant follows directly from U and its equilibrium value in
[0,1]. Every damping interval must be verified positive during execution.

The exact integer j range of each frequency block is retained. Evaluate
the ODE on the whole p interval and multiply by gamma times the integer
count. This encloses the DISCRETE sum; it is not a zero-temperature
quadrature. The published quantity is gamma times the sum, to avoid a
large common period factor. Recovering an unweighted action divides by
the SAME gamma; this changes no physical tolerance.

## Uniform high-frequency tail, n=1,2 only

After the corrected action-adjoint run still missed .02 (its H11 width
was about .03725), the same mesh uses exact polynomial Bernstein ranges
for Q,Q' and the directions and their derivatives. In monomial order,
Q has coefficients [beta,0,u,-4beta-2u-tau,3beta+u+eps,-2eps,eps].
Its degree-six Bernstein coefficient k is sum_(i<=k) a_i binom(k,i)/binom(6,i).
The degree-12 u1 has the sole nonzero coefficient 4096/binom(12,6)
at index 6; the degree-13 u2 has 4096/binom(13,7) at index 7.
Derivative coefficients are n(b_(k+1)-b_k). Restrict to each cell by
the exact de Casteljau subdivision. The Bernstein basis is nonnegative
and sums to one, so the minimum/maximum of the restricted coefficients
encloses the entire polynomial. Arb encloses the coefficient arithmetic
and endpoint uncertainty. This removes interval dependency losses; it
does not alter the geometry, directions, normalization or accuracy target.
Proper lengths are enclosed using the resulting sharper Q ranges.
This brought the n=1 H11 width to .0205847753, still above the original
.02 target. Keep the same argument and coarse (256,64) run and refine
the final level to (640,160); no acceptance target is relaxed.

The current full parameter box has 18<M_n<26 for these levels. A
whole-cell geometry cover verifies .05<f<.45 and |f'|<.15 on
0<=z<=15/16, and Q>.01 on 0<=z<=7/8. The directions and their sum
have |u|<=2, |u'|<=100. On p>=100, elementary differentiation gives

    p/.45<=K<=21p, |c|<=2/p,
    |K_1|<=609/p, |K_2|<=2434/p,
    |c_1|<=600/p, |c_2|<=2500/p.                    (6)

For example, the leading c1 coefficient is
26(100*.45+2*.15)/2=588.9. The remaining coefficient is bounded by
26^3*2*.15*.45^2/100^2. For c2 the full coefficient bound is

    26(2*2*100*.45+2^2*.15)/2
    +2*26^3*.45^2*2*(100*.45+2*.15)/100^2
    +2*26^3*.15*.45^2*2^2/100^2
    +4*26^5*.15*.45^4*2^2/100^4 <2500.

The gap from the inner cut to direction support has proper length
>=log(2)>.25, using Q<=(1-z)^2. Before support Y=Z=0. The X damping
is >=4p, so at and after support

    |X|<=exp(-p)+.5/p^2<=1/p^2,
    |Y|<=151/p^2, |Z|<=636/p^2.                    (7)

For (7), Y forcing is <=601/p and Z forcing <=2540/p using (6),
and damping is >=4p. The same inequalities preserve the stated boxes.
Since |g|<=X^2, |g'|<=2|X| and |g''|<=2, the absolute second-variation
integrand is bounded by

    4[370270/p^5+984354/p^3] <4e6/p^3.

The path after support begins has length <9, by Q>.01. Thus a diagonal
directional second variation is <=3.6e7/p^3, a polarized entry <=5.4e7/p^3.
First variations are smaller. Use 6e7/p^3 for every retained gradient or
matrix entry. If p_J=gamma(J+1/2) is the first omitted frequency,

    gamma sum_(j>=J) p_j^-3 <= 1/(2p_J^2)+gamma/p_J^3. (8)

This supplies a symmetric signed error interval for the COMPLETE omitted
thermal spectrum. The cutoff near 2^20 is chosen to give tail <3e-5.
No claim extends the M<26 constants to the infinite angular tower.

## Execution contract and physical restoration

The initial frozen levels were (64 radial cells per factor two, 8 frequency
blocks per factor two) and (128,16). They missed the accuracy target. With
the sharper signed comparison, (128,32) and (256,64) also missed the
same accuracy target. Use the exact action-adjoint identity (5a), direct
mixed entries and replay at (256,64) and (512,128). Require overlapping narrower final
enclosures, thermal tails <3e-5, and each final matrix interval's absolute
width <.02 in the stated gamma-weighted action units. This is a development
accuracy for the full thermal primitive, NOT an admission tolerance.
If refinement is needed, preserve every attempted run and snapshot and
keep these targets. Compare one retained p=1 case with the independent
raw endpoint-characteristic Hessian MINUS its local flat action Hessian;
this checks the action and contacts, not the continuous-coefficient proof.

Concentrate frequency blocks where they matter: B equal blocks on [0,1],
B per factor two on [1,32], max(4,B/4) on [32,128], and four per factor
two above 128 through 2^20. Each block still encloses all its exact
discrete modes. The same absolute-width target includes any additional
uncertainty from the coarser high-frequency grouping.

Equation (2) is an exact add-and-subtract rearrangement. It does NOT set
the local flat action, physical vacuum or Weyl allocation to zero. The
flat thermal restoration, fixed-scale common a0/a1/a2 terms and angular
tail must still be assembled in the original scheme. Applying an arbitrary
modewise finite renormalization would change the problem. In particular,
these output matrices are NOT the complete renormalized Hessian.

The compact-direction restriction makes the inner normalization independent
of b. Bolt-reaching directions require their own regular-tip normalization
and restored flat thermal contribution; this execution does not silently
apply the compact argument to them. The previous physical j=0 bolt
calculation remains a useful boundary input, not an all-frequency bolt sum.
The weak solve must still retain both gravitational adjoint point sources,
the mass/log phase, canonical-to-York dual map and exact compatibility.

If (1)-(8) and the numerical enclosures pass, the thermal ultraviolet
divergence of this primitive has been removed constructively. Next use
this subtraction within complete angular/restoration and physical-basis
assembly for the selected weak solve. A finite matrix or this one full
thermal sum does not decide controlled-solution feasibility.

The initial startup failure converted an Arb integer without unique_fmpz;
the correction uses its uniquely enclosed integer. The first numerical
failure retains its source snapshot, traceback and coarse spectrum, but
its final failed interval was not appended before the assertion. Subsequent
runs append every spectrum and its bins BEFORE acceptance assertions.
No omitted failed interval is fabricated or described as retained.
The first action-adjoint replay produced nonfinite intervals in the
unnecessary zero-direction prefix. It was deliberately terminated after
that authoritative result. Its unchanged scientific source snapshot,
process identity, observed nonfinite output, diagnostic source and terminal
exit are reconstructed in the receipt with that reconstruction explicitly
marked. Its unavailable per-bin data and checks are not fabricated. The
next replay removes that exactly zero contribution, enforces nonnegative
squares, fails immediately on a nonfinite bin and saves partial bins on
ordinary exceptions. Targets and the original physical domain are unchanged.

Remaining to-do list: complete fixed-scheme restoration, angular tail and
physical bolt/outer weak-action assembly; solve the coupled physical defect
with certified dual error and X/admission; nonlinear closure and compatibility,
matching, causal realization, validity/noise, all-sector stability, charge/return,
opening, transport, preparation, separate creation and the stable-wormhole goal.
