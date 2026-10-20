# Quantitative weighted York return and physical source conversion

2026-10-20 HST, E003/S166. Ordinary proof and contract before execution.
**S166 prerequisite: its full fast TT inverse and mass estimate remain open.**

Keep the analytic seed and the physical inverse of
[S165](tt-form-and-scalar-fourth-graph.md). The previous scalar bounds
used Fsrc=||W fbar_tau||+||W fbar_c||, with
(fbar_tau,fbar_c)=lambda^-2 Rcal* f and W=(r/B)^2.
This note supplies the previously unquantified order-zero source map.
All scalar/vector norms below use dmu=4pi beta r^2 dr. Put

    omega_+=(r/B)^2, omega_-=(B/r)^2,
    R_c c=(S_c,T_c),
    Fdisk=||omega_+ (f_S,f_T)||_L2(dmu;R^2).

For the actual global York return, including its fixed radial tail,

    ||omega_+ R_c c|| <=64||omega_+ c||,
    ||omega_- R_c c|| <=123||omega_- c||,
    Fsrc <=124 lambda^-2 Fdisk.                   (1)

Thus the actual compatible scalar obeys, in the ORIGINAL disk source,

    ||W tau||         <=5.456e9 g lambda^-2 Fdisk,
    ||W gradbar tau|| <=4.96e10 g^(3/4) lambda^-2 Fdisk,
    ||W Tbar tau||    <=3.968e11 sqrt(g) lambda^-2 Fdisk,
    ||W Tbar^2 tau||  <=7.688e13 g^(-1/2) lambda^-2 Fdisk, (2)

where g=G lambda and Tbar=(-Box_4)/lambda. These retain the actual
signed range condition. They do not evaluate the preparation source.
No fast fourth/log norm of the mass-extended shape or numerical mass
coefficient follows from (1). In particular a physical c with a
nonzero mass cubic is not itself in the fast input space of (1).

## 1. Two global polynomial inequalities

Use the [seed reconstruction](stationary-weighted-tt-inverse.md):

    U=r^2F, M=r^2F_r-2rF,
    b=-U^-1 integral_(r_h)^r M c dr,
    eta=-sqrt(F(r)) integral_r^infinity (b+c)/F ds,
    S_c=c/2-R eta, T_c=-(b+3c)/2-(A-R)eta,
    A=F_r/(2sqrt(F)), R=sqrt(F)/r.

The established guards are F>0,F_r>0,F_rr>=lambda,
R^2<=5lambda and F_r<=13lambda r on r>=r_h>B.
Let t=r/B and Uhat=U/B^2. Exactly,

    Uhat=x(t^4-4t+3)+(t-1)^2-delta t
                                      +d(1-2/t+1/t^2).

For x>=1000, 0<delta<=1, |d|<=1, direct differentiation gives

    100t Uhat_t-399Uhat
      =(t-1)[x(t^3+t^2+t+1197)+399-199t
                              +d(599-399t)/t^2]+299delta t,
    100t Uhat_tt-299Uhat_t
      =x(4t^3+1196)-398t+598+299delta
                                      +d(1198/t^3-998/t^2).

Both are positive for t>=1: in the first use
|599-399t|/t^2<=998 and 1000t^3-199t>0;
in the second use the bound 2196 for its d term and
4000t^3-398t>0. The remaining constants are positive.
On the actual exterior U,U_r>0, hence

    r U_r/U >=399/100, r U_rr/U_r >=299/100.        (3)

No extrapolation through the negative-F interval below r_h is used.

## 2. The complete radial TT integral in both weights

The kernel for b is dominated by
U(r)^-1 U_r(s)1_(s<=r), because |M|<=U_r.
For the plus weight the radial measure is m_+=r^6/B^4.
Set nu=m_+/U_r. Equation (3) gives

    d log(nu)/d log(U)
       =(6/r-U_rr/U_r)/(U_r/U)
       <=301/399<4/5.

If the numerator is negative the inequality is immediate; otherwise
use both positive bounds in (3). Thus
nu(r)<=nu(s)[U(r)/U(s)]^(4/5) for r>=s.
[Weighted Schur's test, Q1](https://www.math.ucla.edu/~tao/247a.1.06f/notes2.pdf)
with test function U^(-9/10) now has first integral
10, and second integral at most
1/(9/10-4/5)=10. The endpoint U(r_h)=0 is retained.
Consequently

    ||omega_+ b||<=10||omega_+ c||.               (4)

For the minus weight m_-=B^4/r^2, the already proved monotonicity
of K0=U_r/r^2 gives increasing U_r/m_-=r^4 K0/B^4.
Schur's test with U^-1/2 has both integrals at most two, exactly
as in the [unweighted TT argument](stationary-global-tt-reduction.md).
Therefore

    ||omega_- b||<=2||omega_- c||.                (5)

These are estimates of the actual bolt-to-r integral. Its radial
r^-4 tail is present. No independent homogeneous constant is added.

## 3. Weighted outer-fixed return through the bolt

Write d1=b+c. In proper radius, eta'-A eta=d1. For a power weight
omega=(r/B)^p, integration by parts gives

    ||omega d1||^2
      =integral dmu omega^2[
         eta'^2+(F_rr/2+A^2+2(1+p)AR)eta^2].      (6)

The regular bolt endpoint term vanishes. Initially the sphere input
is smooth and compact at infinity, and the induced b tail is kept;
its return makes the outer term vanish for both p=2 and p=-2.
For p=2 the potential is positive, so

    ||omega_+ A eta||<=||omega_+ d1||,
    ||omega_+ sqrt(lambda)eta||<=sqrt(2)||omega_+ d1||,
    ||omega_+ R eta||<3.2||omega_+ d1||.          (7)

For p=-2 the potential in (6) need not be positive. Instead split
d1 at r=2B and use the exact integral representation of eta.
This is an estimate of the same return, not a change of boundary data.

For input supported below 2B, the unweighted p=0 identity controls
sqrt(lambda)eta by sqrt(2)d1 and A eta by d1.
Since 1/4<=omega_-<=1 in this region, its weighted operator bounds
are respectively 4sqrt(2) and 4. Its output vanishes above 2B.

For r>=2B, the seed supplies

    F>=lambda r^2/2, sqrt(F)<1.001sqrt(lambda)r,
    F_r<3lambda r.

For the upper bounds use F<=lambda r^2+1 and, on t>=2,
F_r/(lambda r)<=2+(4+3/1000)/8+6/(1000*32)
+4/(1000*64)<3. The charge term is nonpositive and is dropped.

The tail Hardy operator r integral_r^infinity d1(s)/s^2 ds has norm
at most two in L2(r^-2 dr). Indeed set v=d1/r and
H=integral_r^infinity v(s)/s ds; integration of (r H^2)' gives
||H||^2<=2||H||||v||, with a favorable lower-end term.
Thus an input above 2B has output ABOVE 2B bounded by

    ||omega_- sqrt(lambda)eta||<5||omega_- d1||,
    ||omega_- A eta||<=6||omega_- d1||.

For its output BELOW 2B put J=integral_(2B)^infinity d1/F ds.
Cauchy--Schwarz gives
|J|<=sqrt(2)/(lambda B^(5/2))||omega_- d1||_L2(r^2dr).
For these inner integrals use this radial norm explicitly; restoring
4pi beta in both input and output cancels in the operator constants.
Using F<=5lambda r^2, F_r<=13lambda r and r_h>B,
the squared inner-output constants are at most 10 and 169/2,
respectively. Round their square roots to 4 and 10.

Combining the disjoint outputs gives outer-input constants 7 and 12.
Combining the two disjoint input pieces by Cauchy--Schwarz then gives

    ||omega_- sqrt(lambda)eta||<=9||omega_- d1||,
    ||omega_- A eta||<13||omega_- d1||,
    ||omega_- R eta||<21||omega_- d1||.           (8)

The bounds first hold on smooth regular-bolt sphere inputs. Density
extends the two integral operators to the weighted L2 spaces. Their
distributional limits preserve the original integral choices and
the smooth bolt realization whenever the input has that regularity.
Neither a uniform pointwise bound on the singular A nor a new bolt
value is introduced.

## 4. Numerical return and physical transpose

By (4),(7), ||omega_+ d1||<=11||omega_+ c||, and hence

    ||omega_+ S_c||<=35.7||omega_+ c||,
    ||omega_+ T_c||<=52.7||omega_+ c||.

Their Euclidean combination is below 64. By (5),(8),
||omega_- d1||<=3||omega_- c|| and

    ||omega_- S_c||<=63.5||omega_- c||,
    ||omega_- T_c||<=104.5||omega_- c||.

Their Euclidean combination is below 123. The corresponding full
Rcal(tau,c)=(tau+S_c,T_c) operator bounds are 65 and 123.
All these statements are in weighted physical L2; no numerical
higher-derivative reconstruction bound is asserted.

The transpose is the one in the actual fixed physical pairing.
Since omega_+ omega_-=1, duality gives

    ||omega_+ R_c* f||<=123||omega_+ f||.

Explicitly, test against v and move the weights:
<omega_+ R_c* f,v>=<omega_+ f,omega_- R_c(omega_+ v)>.
The scalar source f_tau=f_S needs no transformation. Adding these
two bounds and the stipulated factor lambda^-2 proves (1).
Substitution into S164 and S165 proves (2).

For S163's original dual output variables,
f_S=f_theta and f_T=f_theta/3+2f_b.
The matrix has Frobenius norm sqrt(46/9)<3, so alternatively

    Fsrc<=372lambda^-2||omega_+(f_theta,f_b)||.    (9)

With z=B/r the physical weight conversion is exact:
||W f||_dmu=sqrt(4pi beta B^3)||z^-7/2 f||_L2(dz/z).
The common factor is not silently dropped or used to shrink a
Newton product. Compact patch and higher-derivative conversions
remain to be specified when joining the complete inverse.

As a lower-order consequence for the full disk metric,
the old unweighted return norm is at most sqrt(600)<25.
S165 gives ||c||<=5.1sqrt(g)||Dbar^(1/2)c||
<=17850g Fsrc. Therefore

    ||(S,T)||_dmu<=4.5e7 g Fsrc
                      <=5.58e9 g lambda^-2 Fdisk. (10)

The physical mass column is included in (10); this is an unweighted
value bound, not a bound for its coefficient or fast remainder.

## Frozen verification contract

Freeze the proof and checker before execution. Revalidate the exact
receipts and current input hashes of the global TT reduction, weighted
TT inverse, outgoing scalar inverse, S164 and S165.

Verify the two global polynomial certificates and all coefficient
budgets with exact rational arithmetic. Verify the weighted return
integration-by-parts identity with symbolic differentiation of its
density, including the endpoint term. Check Schur's two elementary
integrals for exponents 9/10 and 4/5, and the inverse-weight argument.
Check the inner/outer return constants, all matrix norms and every
conversion coefficient in (1),(2),(9),(10).

As an exact operator diagnostic, set U=r^4 and m=r^6 on (0,infinity).
For c=r^-7/2 cos(k log r) or its sine partner the dominated Hardy
operator is the two-dimensional amplitude matrix
4[(1/2)I+kJ]^-1, J=[[0,1],[-1,0]], of norm
4/sqrt(1/4+k^2)<=8.
Check this by direct differentiation at k=0,1,10,100 and verify the
sharp value 8 at zero. These scale-invariant test profiles are not
normalizable actual sources; they verify the homogeneous kernel and
the direction of the weight, not the whole exterior theorem.

The continuum Schur bound, density and pairing arguments remain
ordinary proofs. This test does not certify the missing full TT graph
estimate, mass coefficient, nonlinear closure or stability. Preserve
every scientific run and exact failed sources if any.

Remaining to-do list: finish S166's fast fourth/log TT inverse, mass
coefficient and full matched-norm residual composition; nonlinear family
and signed compatibility; exterior/state, charge and return; opening,
stability, transport, preparation, creation and full-result audit.
