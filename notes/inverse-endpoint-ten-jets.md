# Quantitative inverse maps and fixed endpoint metric chords

2026-10-21, E003/S167. Ordinary argument before supporting verification.
Keep the original finite comparison, state, subtraction, mass datum and
physical endpoints. The [outer reconstruction](outer-metric-ten-jets.md)
supplies its order-ten logarithms and forward order-eleven maps. Here
their inverse maps and fixed endpoint pullbacks are quantified. All
derivatives use the seed product/orbit metric in A units. There is no
derivative in the original path parameter ell in a moving proper gauge.

Write X=F(Y), G=F^{-1}, where F is either the exact proper map Phi_ell
or the explicit degree-23-cutoff auxiliary map Psi_ell from that proof.
Set L(Y)=1+log(1/Y) and

    N_(m,J)(f)=max_(j<=J) sup |f^(j)|/[j!100^j Y^(m-j)L(Y)].

On the corresponding image of 0<Y<=.1 the inverse satisfies

    N_(5,11)[G(X)-X+ell u3 X^4/4] < ell 1e16.       (1)

The weight and derivatives in (1) are in X. This is conormal C11 on
the open collar, with a C4 boundary extension. A Y^5 logY term is
allowed; ordinary closed-boundary C11 is not asserted. Psi is the
identity for Y>=.1 and has positive derivative, so its inverse is a
global real identification, identity toward the cold region.

For the exact proper endpoint, let v be the logarithm of the inverse
orbit-warping ratio and w the radius logarithm at fixed X. Then

    |nabla^j(v-v3 X^3)|, |nabla^j(w-w3 X^3)|
        <ell 1e17 j!100^j X^(4-j)L(X), j<=10,
    v3=-ell u3, w3=ell(nu3+u3/4), v3=4w3.         (2)

The derivatives in (2) are full ordered seed covariant tensors; the
cubic is subtracted as a scalar first. The same estimate holds for
the tangential and radius logarithms of the global Psi pullback.
Its additional radial logarithm has weighted order-ten reserve 1e19.
It vanishes on the proper outer subcollar. On the transition all three
logarithms have ordinary physical C10 reserve 1e28 and C5 reserve 1e-4.
Farther inward the original middle and cold bounds remain applicable.
These are coefficient inputs, not a finite-metric quantum H6 estimate.

## 1. A shared forward-map remainder and separate small low jets

Put c=ell u3/4, d=F-Y. The forward proof gives, for Phi,

    d=cY^4+R, N_(5,11)(R)<ell 1e14,
    |u3|,|nu3|<1e-11, u3+2nu3=0.

For Psi, R=eta R_Phi+(eta-1)cY^4. The weighted coefficient sum of
eta through eleven is below 2.084. The second term is supported on
Y>=.05. Including all polynomial/product coefficients, its N_(5,11)
norm is below ell 1e-9. Consequently both maps satisfy

    N_(5,11)(R)<ell H, H=3e14.                     (3)

This does not replace small low jets by H. For e=d/Y the unchanged
Phi low estimate, and explicit multiplication by eta, give

    max_(j<=6) |Y^j e^(j)|/(j!100^j)<ell delta,
    delta=3e-10, |F'-1|<2.1e-8.                    (4)

The product bound for Psi is 2.084e-10<delta; it is proved for the
new cutoff, not transferred from the old C6 cutoff. Thus |X/Y-1|<delta.
The ratios Y^m L(Y)/(X^m L(X)), m=4,5, are below 1.001. The same
holds after reversing X and Y. Real monotonicity alone supplies the
inverse; no holomorphic extension of the perturbed metric is used.

## 2. Inversion with one large high-derivative insertion

At a real Y set r=Y/100, a=F'(Y), and form the finite Taylor jet

    [F(Y+rt)-F(Y)]/(ra)=t(1+Q(t)),
    Q=sum_(k=2)^11 F^(k)(Y) r^(k-1)t^(k-1)/(k!a).

All algebra below is truncated after degree eleven. Decompose Q into
Q_L of degrees one through five and Q_H of degrees six through ten.
The coefficient l1 norms obey

    ||Q_L||<5*101*delta/(1-2.1e-8)<epsilon=2e-7,
    ||Q_H||<501 ell H Y^4L.                       (5)

The quartic mass term has no derivative of order seven. Lagrange's
finite inverse formula is

    b_n=(1/n) coeff_(t^(n-1)) (1+Q)^(-n), n=1,...,11.

Since Q_H^2 starts at degree twelve, it contributes nothing here.
The Q_H contribution is exactly

    b_n^H=-[t^(n-1)]Q_H(1+Q_L)^(-n-1).            (6)

Let B(t)=[G(X+Xt/100)-Y]/r. Its coefficients are b_n times
[X/(Ya)]^n. Write B=B_0+B_H, with B_0 obtained by setting Q_H=0.
Then B_H starts at degree seven, and

    ||B_0||<1.00001,
    ||B_H||<4000 ell H Y^4L.                      (7)

Indeed sum_(n=2)^11 ((1-epsilon)^(-n)-1)/n bounds the nonlinear
part of B_0 before its scale factor, and the high sum costs at most
5(1-epsilon)^(-12)||Q_H||. The scale factor to power eleven is
below 1.000001. These are finite-jet identities for a real C11 map.
No bound on a twelfth derivative is assumed.

For any scalar f with a jet through N<=11, split its Taylor polynomial
in r at degree five. If its low coefficient sum is at most epsilon_f,
composition with G obeys the useful tame estimate

    ||J_X(f o G)|| <=1.001 ||J_Y f||
                      +6 epsilon_f ||B_H||.       (8)

Two B_H factors vanish, and inserting B_H into a power of degree at
least six has degree at least twelve. Thus only the derivative of
the degree-five low polynomial multiplies B_H. Formula (8) does not
multiply two large high-derivative bounds.

To retain the inverse's quartic explicitly, use the exact identity

    G-X+cX^4 = T o G,
    T(Y)=-R(Y)+c[(Y+d(Y))^4-Y^4].                 (9)

Finite products, (3),(4), and |c|<=ell 2.5e-12 give
N_(5,11)(T)<1.01 ell H and low coefficient sum <2e-9 ell Y.
For the polynomial correction in (9), terms with R retain its Y^5L
factor; the term with two mass coefficients begins at Y^7. A coarse
coefficient multiplier 1000 on this correction is already sufficient:
1000|c|Y^3<2.5e-12 and its pure mass contribution is below 1e-20
in N_(5,11). Equation (8), with twelve coefficients and the weight
ratio, gives a reserve below 20H<1e16. This proves (1), including
ell=0. Differentiating a parameter-dependent inverse instead would
introduce an unnecessary twelfth derivative; it is not done here.

## 3. Exact pullbacks and the radius pole

Let n(Y)=log N0(Y), where N0 is the seed orbit warping, and let
k(Y)=log[z(Y)/Y]. Their constant parts cancel in all differences.
For X=F(Y), define three scalar fields before composition with G:

    P(Y)=ell u(Y)-log F'(Y),
    A(Y)=ell u(Y)+n(Y)-n(F(Y)),
    W(Y)=w_ell(Y)+log[F(Y)/Y]+k(F(Y))-k(Y).        (10)

The pulled-back radial metric logarithm is P o G, the tangential
one A o G, and the radius one W o G. The inverse orbit-warping
logarithm is v=-A o G. For F=Phi, P=0 exactly. Equation (10) follows
by pulling back e^(2ell u)[dY^2+N0(Y)^2dtheta^2] and r_ell(Y).
In particular z(F)/z(Y) has been split into F/Y and analytic factors;
there is no Taylor series at the singular function log X at zero.

The analytic seed guards imply |n'|,|k'|<|Y| on |Y|<=.15. For k,
write p=z_Y: |p'|=|p n'|<.21|Y| and
pY-z=integral_0^Y t p'(t)dt; the lower bound |z/Y|>1/3 then
gives |k'|<.21|Y|. Local analytic derivative coefficients therefore
have majorant 1.5Y*.02^j on radius-Y/2 disks. Along the real short
segment from Y to F(Y) the same majorant can be enlarged to 2Y*.021^j.

Apply the mean-value identity

    n(Y+d)-n(Y)=d integral_0^1 n'(Y+td)dt

and its k version before estimating. The low part of d is small;
its high part starts at degree seven. Finite analytic composition
has at most one high insertion. Even a reserve 1000H in weight six
fits these two differences, hence 10H in weight four since Y^2<=.01.
Their order is at least Y^5, so neither supplies a cubic.

For log(F/Y)=log(1+e), divide (3) by Y first. Multiplication by
1/Y has normalized coefficient sum 1/(1-.01). Remove cY^3.
The low norm of e through degree five is below 6delta; the high
ideal has degree six and square zero through degree ten. The log
remainder has N_(4,10) reserve 3H, including its quadratic mass term.
The original finite w_ell remainder has reserve 3e13; it is the
original log-gap completion, not a newly chosen radius interpolation.

Subtract cubics in the final X coordinate by forming

    T_v=-A+ell u3 F^3,
    T_w=W-ell(nu3+u3/4)F^3.                       (11)

The additional F^3-Y^3 starts at Y^6. Each expression in (11) has
N_(4,10) reserve 15H and low coefficient sum below ell 1e-8.
Using (8) and the weight ratio gives ordinary reserve below 200H.
The outer full-tensor conversion from the preceding proof costs
less than 1.002, so (2) follows with reserve 1e17.

For the radial logarithm, the normalized coefficient of F'-1 at
order j is 100(j+1)d_(j+1)/Y. The five high coefficients j=6,...,10
have sum below 4500 ell H Y^4L. Its low sum through five is below
21*101*delta<1e-6. The logarithm has at most one high insertion.
For Phi and the proper subcollar P=0 exactly. On the remaining
transition Y>=.05, its low sum is below 2e-6 and converts to the
Y^4L weight with cost below one. This retains the low mass cancellation
before taking a boundary-weighted bound. Including ell u now gives
||J_Y P||<6000 ell H Y^4L and low sum <2e-6.
Composition (8) and tensor conversion fit 1e19. On Psi's transition
X>.049 and X<=.1, and max X^(4-j)L(X)<3e8 through ten; hence its
unweighted physical reserve is 1e28. For Y>=.1 the map is identity.

Finally, u3+2nu3=0 yields v3=4w3 exactly. The common-coordinate
relation and the proper-coordinate relation are different equations;
neither coefficient pair is silently substituted for the other.

## 4. Verify the new transition's low bounds separately

The new cutoff does not inherit the old transition's estimates by
name. Here is a direct low-order ledger on .05<=Y<=.1. In unweighted
Taylor units 1/100 the absolute u coefficients through five are below

    U=(3.32e-11,3.32e-10,3.32e-9,3.32e-8,4.01e-7,8.01e-6).

These follow from the retained 1e-7 outer remainder and its mass
cubic, using L(.1)<3.31 and L(.05)<4. The value of D=Phi-Y is below
8e-13 ell by integrating exp(ell u)-1. Bound its coefficients by

    D_0=8e-13,
    D_j=[t^(j-1)]exp(U(t))/(100j), 1<=j<=6,

where for j=1 subtract 1 from the exponential's constant coefficient.
Use exp(U_0)<1+1.001U_0 and the finite Taylor expansion in the
positive-degree polynomial. Coefficients of eta have majorant
E_0=1, E_j=1.25*4^j. Convolve E and D, then differentiate once:

    d_j=(E*D)_j, h_j=100(j+1)d_(j+1), j<=5.       (12)

The checker evaluates (12) rationally. The coefficient sum for
log F' is bounded by sum(h)/(1-sum(h))<2e-5. The u sum is below
9e-6. The analytic warping difference has coefficient sum <1e-7;
the radius pole and analytic k difference together have sum <2e-7.
For example sum_(j<=5)d_j<4e-9 and division by Y has coefficient
sum at most 1/[.05(1-.2)]=25. The analytic derivative multipliers
above cost less than one in these unweighted units. The original
finite w_ell low sum is below 2e-5, by the retained u,nu,V low
remainders and the exact exp/log-gap formula.

Lagrange inversion at radius 1/100 now has Q coefficient sum below
100 sum_(j=2)^6 d_j/(1-2.1e-8)<1e-5. Its composition cost through
five is below 1.01. The seed connection conversion costs below 1.1:
on the transition's radius-.025 disks |beta|<.1, so use the ordered
tensor recurrence with Gamma_i=.3 i!40^i. Thus every transition
metric logarithm fits ell 1e-4 j!100^j through five, with slack.

For independent very low geometric checks, the sharper real bounds

    |D|<8e-13, |D'|<3.4e-11,
    |D''|<3.4e-8, |D'''|<6.7e-5

give |Psi'-1|<1.2e-10, |Psi''|<3.7e-7 and |Psi'''|<6e-4.
Here use |eta'|<100, |eta''|<4e5, |eta'''|<4.8e8.
The extra gauge logarithms have C0<1e-8, C1<1e-6, C2<1e-3,
including the inverse scale and analytic seed/radius shifts. These
checks do not replace the smaller original-field low contracts.
The exact Phi endpoint itself is unchanged, so its earlier 1e-7
proper C5 remainder contract continues to apply.

## Verification scope and next step

The supporting checker verifies exact inverse/cubic/pullback identities,
rational tame-jet and transition ledgers, ordered tensor amplification,
and independent inverse-series/composition coefficients in logarithmic
test families at 60 and 90 digits, with forty guard digits for contour
cancellation. Those examples are not samples of
the actual charged source. The continuum argument is an ordinary proof,
not machine-checked formal verification. Failed runs retain snapshots.
The first frozen run found that the last real transition sum is
.000502, above the initially stated .0005 cap. The corrected cap is
.0006; the same C2<.001 requirement still passes without changing
any input, map, or physical endpoint. The failed sources are retained.

All original regional C10 coefficients and the endpoint C11 inverse
maps/C10 pullbacks are now quantitative. Next bound differentiated
finite-metric quantum response, tensor return and the full H6 Schur
source, retaining the epsilon Box R term. No sufficient transfer target,
exact signed compatibility equation, nonlinear solution, admissibility
or physical stability gate is passed here. Every S167 acceptance task
remains unchecked; the full goal stays active.

Remaining to-do list: finite-metric quantum/Schur estimates; positive
or sharper actual transfer; signed compatibility, nonlinear existence,
physical state/boundary and stability; matching, charge/return,
opening, transport, preparation, separate creation and audit.
