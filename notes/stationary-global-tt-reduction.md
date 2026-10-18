# Global TT projection and a positive classical comparison form

2026-10-17 HST, E003/S139 supporting S074/S075. Ordinary argument before
focused verification. **PROGRESS; the full quantum inverse remains open.**

There is a global transverse-traceless (TT) projection on the actual
fixed comparison exterior, and a quantitative positive TT Hessian for
the exact centered potential seed. These results replace a failed
shortcut: holding the interior principal variable theta fixed gives a
shape diagonal with a real critical-weight normal zero. That diagonal
has no bounded inverse at physical L2 weight.

The positive classical form does not include the retained Dirac
remainder. A one-sided bound for that remainder, including the actual
off-shell contacts, is still needed before reducing the coupled problem
to its scalar phase. The [actual adjoint phase](stationary-adjoint-phase.md),
kernel and source moments remain unevaluated. All prepared data and the
subtraction scheme remain those of the preceding notes.

## 1. The principal shape coordinate is not a global TT coordinate

Use disk tangents (S,T) and theta=S+T/3. Write lambda=lambda_e and
g=G lambda. The [accepted normal family](stationary-disk-and-normal-inverse.md)
is lambda^2 H(s), where

```
H(s)=C(3-s)^T diag(J(s),N_T(s)) C(s),
C(s)=[[1,(s+3)/(3(s+1))],[0,2]], q=s(s-3),
J=3[1+11g/(180pi)+g q/(60pi)](q-4)/(4pi g),
N_T=q(q+2)H_M/(480pi^2)-[1/(48pi g)+23/(8640pi^2)]q,
H_M=F_M(s)+kappa_M,
F_M(s)=psi(s-1)+(pi/2)cot(pi s/2),
39/5<kappa_M<54/5, 0<g<2000^-10.                   (1)
```

On s=3/2+i t set R=t^2+9/4. The fixed-theta direction is
v=(-1/3,1), and C(s)v=(2/[3(s+1)],2). Since |s+1|^2=R+4,
its exact diagonal normal coefficient is

```
D_theta(t)=v^T H(s)v
 = (R-4)/(12pi g)+(35R-44)/(2160pi^2)
                     +R(R-2)H_M/(120pi^2).         (2)
```

The equality uses the reflected congruence; v is real, and C(3-s)
is the conjugate of C(s) on this line. Digamma reflection gives

```
H_M=Re psi(1/2+i t)-(pi/2)sech(pi t)+kappa_M.        (3)
```

It is real and strictly increasing for t>0. Indeed the derivative of
Re psi is the positive convergent series
sum_(n>=0) 2t(n+1/2)/[(n+1/2)^2+t^2]^2; the remaining derivative is
also positive. Since psi(1/2)=-gamma_E-2log 2>-3 and pi/2<2,
H_M>14/5. Hence D_theta is strictly increasing for t>0. At t=0,

```
D_theta(0)=-7/(48pi g)+139/(8640pi^2)+3H_M(0)/(640pi^2)<0.
D_theta(sqrt(7)/2)=2/(45pi^2)+H_M(sqrt(7)/2)/(15pi^2)>0. (4)
```

For the first sign, H_M(0)<54/5<11 makes the loop term less than
1/(10pi^2), already dominated when g<=1. Thus exactly one positive
zero t_* exists. Its negative partner is the only other real zero.
For 0<=t<=sqrt(7)/2<3/2, split off n=0 in the digamma series and
bound its decreasing positive tail by its integral from 1/2. This gives
Re psi(1/2+i t)<2+log(1+2t)<4, hence H_M<74/5. The positive loop
remainder in (2) is then less than 2/pi^2, so

```
0 < 7/4-t_*^2 < 24g/pi < 8g < 8*2000^-10.          (5)
```

This is an order-one log-radius frequency, distinct from the scalar
Planck frequency omega of the full normal system. At t_* the full
matrix is nonsingular: N_T>0 and J<0, since R<4 and g<1.

The diagonal failure is actual, not just a formal symbol observation.
Conjugate physical L2 by z^(-3/2), put x=-log z, and choose normalized
smooth packets L^(-1/2)chi((x-L^2)/L)exp(i t_* x). Multiply back by
z^(3/2). Their H^(4,log) domain norms stay finite and bounded away
from zero, and they converge weakly to zero. The normal diagonal
residual tends to zero as O(1/L), by smoothness of its multiplier near
t_*. The proved [actual-minus-normal tail estimate](stationary-global-weighted-map.md)
and its compact separated remainder make the actual diagonal residual
tend to zero as well. The compact-core term vanishes. A Fredholm
estimate would contradict this sequence. Thus the actual fixed-theta
diagonal is not Fredholm on the critical-weight shape domain. Its
off-diagonal residual need not vanish, so this is not a full coupled
kernel, a physical instability or a nonexistence theorem.

## 2. Global kinematic York projection at the actual comparison

Work at any one of the accepted fixed positive-temperature comparison
metrics, including the finite candidate, in proper radial coordinate:

```
ds^2=d rho^2+a(rho)^2 d varphi^2+r(rho)^2 dOmega^2,
A=a'/a, R=r'/r, w=A+2R, m=a r^2.                    (6)
```

Primes in this section mean proper-radius derivatives. The angular
period is fixed. At the bolt a=rho+O(rho^3) and r=r_h+O(rho^2), with
the harmless period normalization absorbed in a. At infinity A,R tend
to sqrt(lambda), with the accepted two-power differentiated remainder.

The tracefree part of a disk tangent is t=diag(T,T,-T,-T) in the
orthonormal (rho,varphi,sphere,sphere) frame. Its divergence is
T'+4RT. For a radial vector Y partial_rho let L_cY be the tracefree
Lie derivative. Its three distinct components and divergence are

```
p0=3Y'/2-wY/2,
p1=-Y'/2+(3A/2-R)Y,
p2=-Y'/2+(R-A/2)Y,
P Y=div(L_cY)=3(Y''+wY')/2-VY,
V=w'/2+2A^2+4R^2.                                  (7)
```

Solve P Y=T'+4RT on the smooth-bolt, decaying outer domain. Then

```
h_TT=t-L_cY,
tau=S+T/2-(Y'+wY)/4,
h=2tau g+L_Y g+h_TT.                               (8)
```

This is a genuine global kinematic decomposition, without dividing by
A-R, which can vanish in the exterior. Its existence follows directly
from the radial operator. With measure m d rho,

```
<Y,-PY> = (1/2)||L_cY||^2
         = integral m[3|Y'|^2/2+V|Y|^2]d rho.      (9)
```

Use its Friedrichs form domain. The bolt roots are rho and rho^-1;
the second has infinite form energy. After conjugation by sqrt(m) the
inverse-square coefficient, divided by the leading 3/2, is 3/4. Thus
the smooth vector bolt has the unique finite-energy realization. At
infinity the conjugated potential tends to 75lambda/8>0. Smooth
compact-core Rellich compactness and this positive limiting potential
put the essential spectrum strictly above zero. The same argument
applies with the finite C5 coefficient regularity already established.

If (9) vanishes, then p0-p1=2(Y'-AY)=0 and p0-p2=2(Y'-RY)=0.
A nonzero solution would require a proportional to r, contrary to
a=0 and r_h>0 at the bolt. Equivalently Y vanishes near the bolt and
its first-order equation extends that zero. The kernel is therefore
zero. A positive spectral gap exists for each fixed comparison, giving
a bounded physical-L2 inverse for -P. This argument does not give a
uniform numerical gap over the full candidate family.

The outer homogeneous roots of P are z^-1 and z^4. A bolt-regular
homogeneous solution and an outer-decaying homogeneous solution have
nonzero Wronskian, since a common solution would contradict (9).
Their Green kernel gives the weighted inverse for -1<delta<4; its
log-coordinate tails have positive rates delta+1 and 4-delta. Finite
coefficient corrections and the regular bolt do not change this
conclusion. It applies in particular to delta=7/2 and the mass cubic.

With the actual finite C5 coefficient contract, ordinary weighted
regularity gives Y in H5 and h_TT in H4 for T in H4. It does not by
itself prove the sharp logarithmic upgrade: the third derivative of
V contains the fifth metric derivative, whose logarithmic multiplier
bound is not supplied by boundedness alone. On the explicit smooth
conormal seed of section 3, the usual weighted multiplier estimates do
give T in H^(4,log) => Y in H^(5,log). There the leading part of (7)
gives Y'=2T/3 plus an operator of order minus one on T, so tau=theta
plus a lower-order term and h_TT retains the shape logarithm. Extending
this sharp mixed projection map to the finite candidate requires an
additional coefficient estimate. At exact AdS,

```
Y=-2T/[3sqrt(lambda)(s+1)],
tau=S+(s+3)T/[3(s+1)].                             (10)
```

This recovers (1). At s=4 the global homogeneous coefficient must be
retained; the cancellation in the particular normal solution does not
remove a global datum.

## 3. Quantitative projection bound for the centered seed

Now restrict to the exact [centered potential seed](stationary-centered-covariant-seed.md).
Write t=r/B, x=lambda B^2, v=r_h/B=1+h, and

```
F=1-2m0/r+q0/r^2+lambda r^2
                  +d(B^2/r^2-2B^3/r^3+B^4/r^4),
m0=B(1+2x+delta/2), q0=B^2(1+3x),
x>=1000, 0<delta<=1, |d|<=1.                       (11)
```

These deliberately coarse guards follow from the retained sharper
bounds: |d|<1e-12, h<1e-30, delta/C_star<=1.001h^2 and
C_star<8k^3 at k=2000. No preparation interval is narrowed.
Here a=sqrt(F)/kappa and proper differentiation is sqrt(F)d/dr.
Equation (7) becomes

```
V=F''/4+F'/(2r)+3F/r^2+3(F')^2/(8F),
F''/4+F'/(2r)=3lambda/2+q0/(2r^4)
             +(d/B^2)[1/(2t^4)-3/t^5+3/t^6].       (12)
```

Primes from now on denote areal-radius derivatives. Since F>0 outside
the bolt and the last bracket has modulus at most 13/2,
V>=(149/100)lambda. Hence ||(-P)^-1||_(L2->L2)<=100/(149lambda)
for this seed. This is an inverse for the projection operator, not
the physical Hessian.

## 4. A positive classical TT Hessian on the exact seed

Use relative orthonormal metric components h_TT=diag(A0,B0,C0,C0),
so A0=-B0-2C0. To avoid confusing C0 with the benchmark constant,
write c=C0 and b=B0 below. Transversality is exactly

```
(r^2 F b)'=-(r^2 F'-2rF)c,
b(r)=-[r^2F(r)]^-1 integral_(r_h)^r (y^2F'(y)-2yF(y))c(y)dy. (13)
```

This is the unique bolt-regular TT field for a sphere component c.
At the bolt b=A0=-c, and simple-zero division gives smooth even
components in proper radius. No photon-sphere denominator is used.

First take c smooth and compact away from the bolt and infinity.
The resulting b is zero near the bolt and has an allowed r^-4 tail.
The exact seed is stationary for S0=S_classical+Gamma_ref. Thus its
Hessian is invariant under an admissible infinitesimal radial
diffeomorphism. In areal gauge, the tangent is

```
delta alpha=-(3c+rc')/2,
delta F=-Fb+F(c+rc')-rF'c/2.
I0=beta/(2G) integral exp(alpha)[(rF)'-V0(r^2)/2]dr
                                                + boundary terms,
I0''=beta/G integral delta alpha (r delta F)'dr.    (14)
```

The action is the Euclidean reduced action; the boundary terms have
zero variation on this core. Formula (14) follows by integrating the
orbit curvature and using (rF)'=V0/2. It is not an inference from
thermal monotonicity. For comparison, at fixed conformal factor the
unintegrated reduced action is affine in f, so its raw f-f Hessian
vanishes for compact affine f variations.

Here is an explicit integration certificate for (14). Put

```
L=3rF+r^2F'/2, K=-F+rF'-r^2F''/2,
M=-2rF+r^2F'/2+r^3F''/2,
(r delta F)'=r^2F c''+L c'+K c+Fb,
V_F=6F-7rF'/2+r^2F''/2-r^3F'''/4,
J=-3r^2F cc'-r^3F(c')^2/2-rFbc+Mc^2/2,
-(3c+rc')(r delta F)'
    =3r^2F(c')^2/2+V_F c^2-4Fbc+J'.               (15)
```

The last identity uses (13). J vanishes at both ends on the core.
Completing c'+2c/r gives

```
I0''=beta/(2G) integral [3r^2F(c'+2c/r)^2/2
                                  +W_F c^2-4Fbc]dr,
W_F=3F-rF'/2+r^2F''/2-r^3F'''/4
 =3-12m0/r+13q0/r^2+3lambda r^2
                       +d(13/t^2-51/t^3+45/t^4).  (16)
```

The cross term has a useful lower bound. Let
H(r)=integral_r^infinity c(y)/y^2 dy and P0=r^4F'-2r^3F. Fubini
and integration by parts in (13) give the exact identity

```
-4 integral Fbc dr
 =2 integral P0' H^2 dr+2r_h^4F'(r_h)H(r_h)^2,
P0'>=-6r^2.                                       (17)
```

The second term is positive and is retained. The bound on P0' follows
from P0'=-6r^2+12m0 r-4q0+d[-4B^2+6B^4/r^2] and (11).
Also integration of (r^3H^2)' and Cauchy--Schwarz gives
integral r^2H^2 dr <=(4/9)integral c^2 dr; the lower endpoint term
has a favorable sign. Therefore (17) is at least
-(16/3)integral c^2 dr.

The potential obeys W_F-16/3>=lambda r^2. Indeed its x coefficient
is 3t^2-24/t+39/t^2>2t^2: the maximum of 24t-t^4 is
18*6^(1/3)<36<39. Its remaining RN part is
3-12/t+13/t^2>=3/13. The thermal term is at least -6delta/t and
the d term at least -109|d|/t^2. The coarse guards (11) then leave
more than x t^2 after subtracting 16/3. Thus

```
I0'' >= E0[c]
 := beta/(2G) integral [3r^2F(c'+2c/r)^2/2
                                    +lambda r^2 c^2]dr. (18)
```

This controls the full TT norm as well. Set U=r^2F and
K0=U'/r^2=F'+2F/r. One has F'>0 and K0'>0 throughout the exterior.
For explicit guards,

```
B^2F'' >= 52x/27-1/8-2delta-50|d|>0,
B^2K0' >= 4x-2-18|d|>0.                           (19)
```

The first bound splits the x and constant terms, whose minima for
t>=1 are 52/27 at t=3 and -1/8 at t=2. F'(r_h)>0 is already proved.
Since |r^2F'-2rF|<=U', the kernel of (13) is dominated by
U(r)^-1 U'(y)1_(y<=r). On L2(r^2dr), use the Schur function
U^-1/2. The first Schur integral is exactly 2. The second is at most
2 because K0 is increasing: r^2<=U'(r)/K0(y) for r>=y.
Consequently ||b||<=2||c|| and ||A0||<=4||c||, so

```
integral r^2(A0^2+b^2+2c^2)dr <=22 integral r^2c^2dr,
I0'' >= beta lambda/(44G)
              integral r^2 |h_TT|^2dr.            (20)
```

The physical four-volume norm includes 4pi beta, hence its gap is
lambda/(176pi G). This is the restricted S0 TT quadratic form,
not a claim that the unprojected Einstein-Maxwell operator preserves TT.

Extend from the core to smooth finite-energy TT fields as follows.
At the bolt a logarithmic cutoff of c in rho has gradient cost tending
to zero because r^2F(c')^2 dr is a constant times rho(c_rho)^2 d rho.
The region where it differs has vanishing volume. Equation (13), its
bounded integral operator and the regular Taylor expansion extend the
TT components in H1; in the transition the additional derivative is
controlled by the same logarithmic cutoff cost. At infinity ordinary
cutoffs converge for the stated finite-energy decay. The full S0
Hessian is a continuous first-derivative form on this smooth geometry.
Thus (18)-(20) pass to the closure, including smooth bolt values and
the physical fast/mass tails. No singular areal-gauge shift at a moving
bolt, or omitted entropy endpoint term, is used in this extension.

## 5. What the positive form does and does not resolve

For the seed, (8) and (13) supply a genuine global TT space and a
positive comparison form. A TT field can be returned to the fixed disk
slice without a singular photon-sphere gauge: in proper radius set
eta=-a integral_rho^infinity (b+c)/a d rho. Then
eta'-A eta=b+c, eta=O(rho) at the bolt, and subtracting L_eta g makes
the two orbit components equal. Its normal disk coefficients are
(S,T)=((s+3)/(2(s+1)),-3/2)c for the tangential TT branch, reproducing
tau=0 and b_normal=-3c in (1). The separate radial homogeneous z^4
term must also be retained; its coefficient is fixed by (13), and can
be nonzero even for compactly supported c. The displayed branch formula
does not remove this s=4 datum. The gauge integral converges for the
selected positive weights.

This kinematic decomposition does not annihilate gauge terms in the
actual off-shell full Hessian. Only S0 is stationary at its seed.
For the full action, the pullback must retain gradient, chart and Lie
contacts. In an affine tensor chart, diffeomorphism invariance gives
I''(h,L_Yg)=-I'(L_Yh); it does not set the left side to zero when
I'(g) is nonzero. The corresponding chain-rule terms in disk variables
are part of the required actual response.

Nor is the retained loop a bounded small perturbation in the classical
energy norm: compact frequency-p TT packets have a positive leading
loop/classical ratio of order G p^2 log p, unbounded as p grows. The
needed estimate is one-sided. At the seed a sufficient bound would be

```
Gamma_rest''[h_TT,h_TT] >= -eta0 E0[c], eta0<1,      (21)
```

with the actual pullback contacts and candidate transfer also controlled
before applying it there. An explicit lower bound of the form
-beta C_ren[lambda integral r^2F(c'+2c/r)^2
+lambda^2 integral r^2c^2], with 2G lambda C_ren<1, would suffice.
The positive fourth-order logarithmic term can remain on the left.
The normal symbol supplies control sufficiently far out; no accepted
estimate yet bounds the negative finite part and contact form throughout
the compact/bolt region. Existing C5-to-C0 response value bounds are
not that quadratic-form estimate.

[Monteiro and Santos, sections II and IV and Appendix B](https://arxiv.org/pdf/0812.1767)
show why thermodynamic and TT arguments must be distinguished in
Einstein-Maxwell theory. Their asymptotically flat magnetic RN result
uses a reduced canonical sector and cannot supply (21) for this AdS
potential and retained nonlocal action. In particular their Appendix B
warns that the unprojected Einstein-Maxwell Hessian need not preserve TT.
The present statement is explicitly a restricted quadratic form.
[Louko and Winters-Hilt, section IV](https://arxiv.org/pdf/gr-qc/9602003)
reduce RN-AdS to mass and charge after solving its classical constraints;
that reduction does not eliminate the retained quantum response here.

The next unresolved scientific step is the actual one-sided quantum
form bound and the sharp mixed projection coefficient estimate, or
another controlled route to the coupled inverse. The
physical kernel, actual adjoint angle and candidate source moments
remain open. No complete exterior, Hamiltonian, physical stability,
instability or nonexistence follows from the positive classical form.

## 6. Focused verification contract

Complete ordinary review precedes execution. A finite checker may
verify (2), endpoint signs and rational guards, York components and
normal reconstruction, the seed potentials, TT identity (15), square
completion, Hardy/cross identities and the gap constants. The escaping
packets, actual projection mapping, Schur estimates, bolt closure and
quantum limitation remain ordinary proofs. No actual determinant
evaluation, phase shooting, candidate-moment computation, numerical
root fit or Lean campaign is selected. Freeze reviewed source bytes;
retain all runs, failed input bytes and raw UTF8/LF/noBOM hashes.

Remaining to-do list: bound the actual quantum TT form and off-shell
contacts and justify the candidate's sharp mixed projection map;
determine the physical kernel and adjoint phase, evaluate
candidate moments and quantify the compatible inverse; close nonlinear
physical logarithmic phase or justify controlled spatial reduction;
match exterior/state and charge, transfer S138, and complete driven
stress, radius, controller, capture, tail, spatial area and S074-S080.
