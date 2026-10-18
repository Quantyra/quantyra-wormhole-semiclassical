# Seed basepoint, off-shell contacts and the flat quantum form

2026-10-17 HST, E003/S139 supporting S074/S075. Ordinary argument before
focused verification. **PROGRESS; the curved quantum TT bound is open.**

Use the analytic centered seed g0 as the next linear basepoint for the
FULL retained action. Its sharp mixed York projection is already proved.
This avoids requiring new regularity of the finite candidate before
studying the linear operator. It does not replace the action by S0, or
establish a nonlinear solution near the seed.

The full gradient's gauge and disk-chart contacts at g0 cost less than
1e-6 of the positive classical TT comparison energy. The canonical LLL
Hessian costs less than 1e-25. The common local a2 curvature allocation
has negative part below 1e8(G lambda)E0. The exact flat massive curvature multiplier
is positive after the existing common four-dimensional subtraction.
That last statement concerns a flat model, not the actual curved TT form.

## 1. Basepoint and conventions

Retain the prepared period, boundary clock, magnetic flux, paired Dirac
domain, common scale mu0 and full action allocation of the
[centered seed](stationary-centered-covariant-seed.md) and
[complete seed source](stationary-centered-exterior-source.md):

```
I=I0+Gamma_rest, I0=S_classical+Gamma_ref,
Gamma_rest=W4+Gamma_hat,LLL+Gamma_hat,mass-Gamma_ref,
I0'(g0)=0, I'(g0)=Gamma_rest'(g0).
```

The [actual weighted operator](stationary-global-weighted-map.md) and
[adjoint phase theorem](stationary-adjoint-phase.md) permit this fixed
admissible comparison. The seed is smooth and conormal, so the sharp
mixed projection of [the global TT theorem](stationary-global-tt-reduction.md)
applies here. A later transfer to a finite candidate still needs the
appropriate coefficient and nonlinear estimates. No candidate C6 or
higher-source-derivative campaign is selected by this decision.

Write lambda=lambda_e, t=r/B, x=lambda B^2 and

```
F=1-2m0/r+q0/r^2+lambda r^2
                    +d(B^2/r^2-2B^3/r^3+B^4/r^4),
m0=B(1+2x+delta/2), q0=B^2(1+3x),
x>=1000, 0<delta<=1, |d|<=1.
```

The actual |d|<1e-12 and thermal guards remain available. In proper
radius rho, A=a'/a, R=r'/r, a=sqrt(F)/kappa and m=a r^2. A prime
means d/d rho below unless marked by an r subscript. The harmless
constant in m gives the same integrations by parts as r^2 dr. Use
norms with measure r^2 dr; physical four-volume is 4pi beta r^2 dr.

The accepted TT lift is h=diag(-b-2c,b,c,c) in the
(Euclidean time, proper radius, sphere, sphere) orthonormal frame.
In particular b is the radial component. It obeys

```
(r^2 F b)_r=-(r^2 F_r-2rF)c,  ||b||<=2||c||,
b'=-2A(b+c)-2R(b-c), ||h||^2<=22||c||^2.
C0=||c||^2, D=||c'+2Rc||^2, E=D+lambda C0,
E0=beta/(2G)[3D/2+lambda C0] <= I0''[h,h].       (1)
```

For sections 2-5 initially c is smooth and compact away from both ends.
Its induced b has the permitted r^-4 tail. The established logarithmic
bolt and outer H1 closure extends those estimates. Section 7 instead
uses smooth bolt-regular approximants for its H2 form. In particular,
the separate radial z^4 datum is not discarded.

## 2. First- and second-derivative estimates for the gauge return

The seed polynomial and the previous convexity certificate give

```
lambda <= F_rr <=21lambda,
R^2=F/r^2<=5lambda, 0<=AR=F_r/(2r)<=11lambda.     (2)
```

For the lower bound use B^2F_rr>=52x/27-1/8-2delta-50|d|>x.
For the upper bound, drop the negative m0 term and use
2x+6(1+3x)+50|d|<21x. The other two bounds follow by coefficient
sums and F_r>0. Set d1=b+c, e1=b-c; both norms are at most 3||c||.
The TT identity is d1'+2(A+R)d1=c'+4Rc. Multiplying by A d1 and
integrating, using A'=F_rr/2-A^2, gives

```
integral [2A^2+AR-F_rr/4]d1^2
                         =integral A d1(c'+4Rc),
||A d1||^2 <= (1/2)||c'||^2+88lambda C0.          (3)
```

Indeed X||c'+4Rc||<=X^2+||c'+4Rc||^2/4, and the last norm
is at most 2||c'||^2+160lambda C0. The remaining contribution
is (21/4)*9lambda C0. No singular bound on A itself is used.
The full tensor derivative, including connection terms, is

```
|nabla h|^2=(-b'-2c')^2+b'^2+2c'^2
                         +8A^2(b+c)^2+4R^2(b-c)^2.
||nabla h||^2<=26||c'||^2+4076lambda C0
              <=5200 E, since ||c'||^2<=2D+40lambda C0. (4)
```

Return to the fixed disk slice with Y=eta partial_rho,

```
eta=-a integral_rho^infinity (b+c)/a d rho,
eta'-A eta=d1, h_disk=h-L_Yg.
||d1||^2=integral [eta'^2+(F_rr/2+A^2+2AR)eta^2]. (5)
```

The integral is convergent, eta=O(rho) at the bolt and it has the
selected decaying outer behavior. Consequently

```
||eta||^2<=18 C0/lambda<25 C0/lambda,
||nabla Y||<=15sqrt(C0), ||L_Yg||<=30sqrt(C0).     (6)
```

For example, ||eta'||^2 and ||A eta||^2 are at most 9C0;
2||R eta||^2<=180C0. The following exact second-derivative
identities exhibit all bolt pole cancellations:

```
|nabla^2Y|^2=eta''^2
       +sum_(k=A,R,R)[(k'eta+k eta')^2+2(k eta'-k^2eta)^2],
eta''=(F_rr/2)eta-A d1-2R d1+c'+4Rc,
A'eta+A eta'=(F_rr/2)eta+A d1,
A eta'-A^2eta=A d1,
R'eta+R eta'=(2AR-R^2)eta+R d1,
R eta'-R^2eta=(AR-R^2)eta+R d1.                   (7)
```

Using the looser ||eta||^2<=20C0/lambda, the five displayed squared
contributions to (7), including their multiplicities, are bounded by

| Term | coefficient of ||c'||^2 | coefficient of lambda C0 |
| --- | ---: | ---: |
| eta''^2 | 15/2 | 12765 |
| (A'eta+A eta')^2 | 1 | 4586 |
| 2(A eta'-A^2eta)^2 | 1 | 176 |
| two radial-sphere first squares | 0 | 58500 |
| four radial-sphere second squares | 0 | 41320 |

Thus ||nabla^2Y||^2<=19D+117727lambda C0<160000E.
It follows that ||nabla L_Yg||<=800sqrt(E). For k1=2h-L_Yg,

```
||k1||<=40sqrt(C0), ||nabla k1||<=1000sqrt(E),
||h_disk||<=35sqrt(C0).                            (8)
```

## 3. The full gradient, including the outer and angular components

Define the physical action gradient E_g by I'(v)=integral E_g:v dvol.
It equals plus or minus one half the full remainder stress, depending
on the Euclidean stress convention; the sign has no role in the bound.
Put U=q^2/(4pi B^4). The complete seed source supplies this value ledger:

| Allocation | bound for each physical orbit component |
| --- | ---: |
| canonical massive, source equation (10) | 1.00005e15 U |
| Weyl, including its vacuum | 1.1e4 U |
| canonical LLL | U |
| full reference subtraction | 5.4e4 U |

Here sphere-integrated stress is divided by 4pi r^2. For the Weyl
row use source equations (20)-(23): V/t^2<q^2, |L|/t^3<q^2,
|RW|/t^2<=2q^2 and |B^2 Delta rho2,W|/t^2<=1e4q^2.
The elementary inequality (1+log t)/t^j<=1 for t>=1,j>=1
controls every outer logarithm. The pressure adds less than 4U:
the explicit E_F formula gives |F Ebar_F|/(2pi t^2)
<=10xq(1+log t)/(pi t^3)<4q^2. The LLL row follows from
|Theta|<=qC, |Theta'|<=qC/4, F increasing and C<q.
The complete reference profiles (35) of the
[horizon source theorem](stationary-horizon-source-defect.md) obey
|pA_ref|,|pS_ref|<3341q^2/B^4<5.4e4 U.

The angular component is reconstructed with the FULL four-dimensional
trace anomaly, including its retained box R term:

```
pS_rest=A4/2-(p2,Dirac-rho2,Dirac)/(8pi r^2)-pS_ref. (9)
```

There is an elementary local anomaly reserve |A4|<1e6 U. To see it,
the seed also has F_r<=13lambda r and every orthonormal sectional
curvature has modulus at most 11lambda. Hence B^4|Rm|^2<=2904x^2,
B^4|Ric|^2<=11616x^2 and B^4 R4^2<=46464x^2. Moreover

```
B^2 R4=-12x+4d t^-5-6d t^-6,
B^4 box R4=t^-2 partial_t[t^2 F partial_t(B^2R4)].  (10)
```

Coefficient bounds F<=5xt^2,F_t<=13xt give |B^4 box R4|<1e5x^2.
The anomaly's fixed coefficients and magnetic term then give
|A4|<4e5q^2/(16pi^2 B^4)<4e4 U, well within the reserve.
No parameter or spatial derivative of a source remainder is inferred.
Every component in (9) is therefore less than 2e15 U in modulus.
The four-component Frobenius norm and the factor one half yield

```
||E_g||_infinity <=1e16 q^2/(4pi B^4).              (11)
```

This deliberately leaves a factor-five reserve. It retains both constant
vacua in a value estimate; it does not assign divergent separate mass
integrals or omit the reference action's finite terms.

## 4. Actual off-shell gauge and exponential-chart contacts

In an affine tensor chart, diffeomorphism invariance implies
I''(h,L_Yg)=-I'(L_Yh). Consequently the exact difference is

```
I''(h-L_Yg,h-L_Yg)-I''(h,h)=I'[L_Y(2h-L_Yg)].      (12)
```

This identity uses the full nonzero gradient at g0. From (6)-(8),
|L_Yk1|<=|Y||nabla k1|+2|nabla Y||k1| and Cauchy give

```
|I'(L_Yk1)| <=4pi beta ||E_g||_infinity *6200 E/lambda
          <2e5 (G||E_g||_infinity/lambda) E0.        (13)
```

The exponential disk chart adds I'(q1^2), where q1=g^-1 h_disk
is the relative tangent and the square is its tensor endomorphism
square. Since |q1^2|<=|q1|^2, (8) gives an additional bound
4e4 (G||E_g||_infinity/lambda) E0. These constants include the
physical 4pi beta measure: the exact upper coefficients used are
49600pi and 9800pi. A common reserve is therefore

```
|all these contacts| <=3e5 (G||E_g||_infinity/lambda) E0
                     <1e-6 E0,                    (14)
```

because the accepted benchmark has Gq^2/(B^2 x)<4e-27 and pi>3.
No inverse horizon displacement, derivative of E_g, sixth metric jet,
or assumption of full-action stationarity enters. On the compact core
the end terms vanish; (14) extends by the accepted TT form closure.
It controls these gauge/chart contacts, not Gamma_rest'' itself.

## 5. The canonical LLL lower form on the disk slice

The [declared canonical LLL allocation](stationary-lll-quantum-exterior.md)
has central charge q and the fixed HH disk domain. Its Polyakov cocycle
for a disk conformal tangent T gives exactly

```
Gamma_hat,LLL''=-q/(12pi) integral_disk |nabla T|^2.
```

The general normalization follows from
[Bilal, de Lacroix and Erbin, equations (1.3)-(1.4)](https://arxiv.org/pdf/2109.03637).
Their compact-surface massive expansion does not supply a global curved
massive estimate here. The disk domain and absence of an extra zero-mode
normalization term are those of the already specified canonical LLL.
Boundary cocycle terms are affine on compact tangents; the form closes
with the same bolt and decaying-end approximation.

For the return (5), T=-(b+3c)/2-(A-R)eta. Using the exact TT constraint,

```
T'=-(3/2)c'+2Rb-E1 eta,
E1=A'-R'+A(A-R)=F_rr/2-F_r/r+F/r^2,
|E1|<39lambda.                                    (15)
```

The coefficient of F/r^2 in this identity is one. Three-square
Cauchy, ||b||^2<=4C0 and ||eta||^2<=18C0/lambda give
||T'||^2<=7||c'||^2+83000lambda C0<=1e5 E.
Two-dimensional conformal invariance gives the disk Dirichlet measure
beta dr, so integral_disk |nabla T|^2=beta integral |T'|^2 dr
<=beta B^-2||T'||^2. Therefore

```
Gamma_hat,LLL'' >=-c_LLL E0,
c_LLL=1e5 Gq/(6pi B^2)<1e-25,                     (16)
```

using q/C<400 and GC/B^2<8.5e-33. W4 and the remaining canonical
functional retain their declared allocations. Equation (16) is a disk
chart Hessian; (14) is the separate full-action chart/gauge comparison.

## 6. Exact flat massive curvature multiplier in the common scheme

This section concerns the translation-invariant flat orbit model, with
the zero-momentum vacuum term kept in its local allocation. It does not
assert positivity of the raw massive determinant. For an anti-selfadjoint
massless Dirac operator A and S=(A+M)^-1, its mass-conformal second
variation before local anomaly terms is

```
M^2 Re Tr(SuSu)-M Re Tr(Su^2).
```

In an A eigenbasis with eigenvalues i lambda_j, its symmetrized
coefficient is -M^2(lambda_j+lambda_k)^2 /
[2(M^2+lambda_j^2)(M^2+lambda_k^2)]<=0. Dropping the one-vertex
contact or treating this raw first-order bubble as positive is invalid.

In the flat rank-two calculation the Feynman numerator is
2[M^2-k dot(k+p)]. With Delta=M^2+alpha(1-alpha)p^2, translation
gives 2(Delta-l^2). Subtracting the p=0 bubble then yields

```
Pi_M(p)-Pi_M(0)=M^2/(2pi) integral_0^1
                   log[1+alpha(1-alpha)p^2/M^2]d alpha.
```

The 2D anomaly contributes -p^2/(12pi), cancelling the first
coefficient. The next coefficient is -p^4/(120pi M^2). Removing
that coefficient leaves the nonnegative expression

```
R(u)=log(1+u)-u+u^2/2=integral_0^u t^2/(1+t)dt>=0,
0<=R(u)<=u^3/3 for u>=0.                           (17)
```

In ell=B/sqrt(C) units, M_n^2=n(n+q)/C and d_n=2(q+2n).
Each angular mode uses this rank-two orbit trace; d_n counts both
angular signs, with no extra spin factor. For fixed p the sum of
d_n M_n^2 R(alpha(1-alpha)p^2/M_n^2) converges absolutely, since
sum d_n/M_n^4=2C^2 H_q^(2)/q<4C^2/q. Thus the convergent positive
remainder can be summed before the common restoration.

The sole divergent p^4 sum is restored by the SAME common heat
prescription as the source notes, not by modewise renormalization.
Indeed zeta_h(s)=sum d_n[n(n+q)]^-s differs from 4zeta(2s-1)
by an analytic function near s=1 with value -2H_q. Therefore
Res zeta_h(1)=2 and Fang=FP zeta_h(1)=4gamma_E-2H_q.
For the two-spinor heat coefficient a2=-(R2^2+2box R2)/120,
the common term is

```
zeta_4,a2(s)=s B^(2s+2) zeta_h(1+s)/(4pi) integral a2,
Gamma_a2=B^2[Fang+4log(mu0 B)]/(8pi) integral a2.    (18)
```

This follows by differentiating mu0^(2s) zeta_4,a2(s)/2 at zero.
There is no additional finite +2 in the action coefficient. Since
linear R2=-2Delta u, (18) restores exactly the p^4 multiplier
-C[Fang+4log(mu0 B)]/(120pi) in ell units. The flat curvature
multiplier is consequently

```
K_flat(p)=-C[Fang+4log(mu0 B)]p^4/(120pi)
 +sum_(n>=1) d_n M_n^2/(2pi) integral_0^1
          R(alpha(1-alpha)p^2/M_n^2)d alpha >=0.     (19)
```

For q>=4, H_q>=25/12 and gamma_E<1 give Fang<-1/6.
The actual |log(mu0 B)|<1e-21 yields the explicit margin
K_flat(p)>=C p^4/(1000pi). This is an exact all-frequency flat
model statement. The zero-momentum local allocation, the Weyl action,
reference potential and curved/paired-boundary effects have not been
bounded below by (19).

## 7. The restored a2 allocation also has a curved lower bound

The explicit local term (18) is alpha2 integral_hat R2^2, where
alpha2=-B^2[Fang+4log(mu0 B)]/(960pi)>0. The integrated box R2
term is a boundary term. Start here with smooth bolt-regular T and
compact outer support. The bolt is a smooth interior point, so it
contributes no boundary term. With hat g varying as exp(2sT)hat g,
the exact Hessian is

```
Gamma_a2''=alpha2 integral_hat [8(Delta T)^2
       -16R2|nabla T|^2+(4R2^2+8Delta R2)T^2].     (20)
```

This follows directly from R2(s)=exp(-2sT)(R2-2sDelta T).
The seed's compactified orbit polynomial calB satisfies
calB<=C/5, |calB'|<=5C, |calB'''|<=17C, |calB''''|<=13C,
with derivatives in u=B/r. For ell=B/sqrt(C),

```
|ell^2 R2|<=2,
ell^4 Delta R2=-(calB calB''''+calB' calB''')/C^2,
|ell^4 Delta R2|<100.
```

Keep the positive 8alpha2||Delta T||^2. Since |Fang+4log(mu0B)|<100,
the rest of (20) is bounded below by
-2C integral_hat |nabla T|^2-30C/ell^2 integral_hat T^2.
From (5)-(6), ||T||<=||b+3c||/2+||A eta||+||R eta||
<15sqrt(C0), so the looser ||T||^2<=500C0 suffices.
The hat-volume integral is beta integral (B^2/r^2)T^2 dr
<=500beta C0/B^2. Combining this with section 5 and C/x<6.1 gives

```
Gamma_a2'' >=8alpha2||Delta T||_hat^2-1e8(G lambda)E0. (21)
```

The coefficient before enlargement is at most
4e5(C/x)+3e4(C/x)^2<4e6. This is an actual curved estimate for
the declared common local a2 allocation, not an estimate of the
nonlocal remainder or W4. Extend on the H2 form domain using smooth
bolt-regular approximants and outer cutoffs. The logarithmic cutoff
that removes the bolt in H1 is not an H2 approximation for a nonzero
bolt value; it is not used for this curvature-squared term.

## 8. Exact remaining obstacle and verification scope

A lower bound for the actual curved remainder is still needed. The
positive flat multiplier does not by itself give a positive global
curved operator or control their difference. In particular, local
normal-coordinate Green bounds must be joined with the complete finite
contact left after (20), a quantitative comparison to a positive global
model, and the W4/reference allocations. W4 contains the growing
background log(r/B); its fourth-order coefficients must be combined
with the common canonical operator, not bounded separately by a false
global constant. The accepted outer normal estimate controls tails
qualitatively but does not yet supply the needed interior lower form.
Existing C5-to-C0 source/response value bounds cannot replace that form
estimate. No numerical value is assigned to an unproved remainder.

The focused checker may verify the exact differential identities,
connection-square algebra, seed/anomaly formulas, rational contact
constants, LLL normalization, Feynman parameter coefficients, curved a2
variation and flat
common-heat algebra. Global closure, source-ledger transfer, Feynman
integration and angular convergence remain ordinary arguments. It
does not compute the curved determinant, physical kernel, adjoint
phase or candidate moments. Freeze reviewed sources before execution;
retain every run and failed input bytes with raw UTF8/LF/noBOM hashes.

Remaining to-do list: bound the actual full curved quantum TT remainder;
determine the coupled physical kernel and adjoint phase, evaluate source
moments and quantify the compatible inverse; transfer regularity and
close nonlinear physical logarithmic phase or justify controlled
spatial reduction; match exterior/state and charge, transfer S138, and
complete driven stress, radius, controller, capture, tail, spatial area
and S074-S080.
