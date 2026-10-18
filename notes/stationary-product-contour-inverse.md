# Exact product bubble and the complete frequency-contour inverse

2026-10-17 HST, E003/S139 supporting S074/S075. Ordinary proof before
focused verification. **PROGRESS; actual exterior and nonlinear control open.**

The complete stationary reduced product Hessian has an explicit inverse
at every frequency on both cold growth contours. The ordinary complete-H2
spectral operator also has an inverse, with normalized L2 bound 16000.
The exact massive pole residues and the restored finite quadratic
coefficient are retained; this is not a large-mass replacement.

The [weighted spatial response](stationary-product-weighted-response.md)
reconstructs this product inverse and gives a bolt-regular particular
response with explicit finite-neck weights. The actual variable-radius
seam and fixed mass/temperature return remain separate.

Continue the [reduced cold response](stationary-reduced-cold-response.md)
and [covariant reduced inverse](stationary-covariant-reduced-inverse.md),
with their fixed theory, state, angular multiplicities, common heat
subtraction, prepared temperature and actual reference coefficients.

## 1. Inputs and the normalized coupled matrix

Use the notation of the reduced cold-response theorem:

    lambda=lambda_e, x=lambda B^2>7.9e9,
    C=B^2/A^2=1+6x+d, |d|<1e-12,
    epsilon=G/(60pi), eta_G=GC/B^2<8.5e-33,
    q=2000^4, M_n^2=n(n+q)/C, d_n=2(q+2n), M_1>18.

The accepted angular guards are C/q<.00352, q/C<400 and

    sum_n d_n/M_n^3<2C/3,
    F_eff=4gamma_E-2H_q+4log(mu0 B),
    |F_eff|<100, |log(mu0 B)|<1e-21.                     (1)

The fourth moment follows directly, without an additional estimate:

    sum_n d_n/M_n^4
     =(2C^2/q) sum_n [1/n^2-1/(n+q)^2]
     <4C^2/q<.01408C.                                  (2)

The product metric is
g4=e^(2S)[e^(2T) A^2 g_H2+B^2 g_S2]. Its physical scalar pairing is
4pi B^2 A^2 dvol_H2. Put p=-Delta_H2+2 and gamma0=(8pi G A^2)^-1.
The EXACT stationary reduced product Hessian of the reduced cold-response theorem is gamma0 M(p),
where

    M(p)=[[-a(p),-2b p],[-2b p,q_cold(p)]],
    c=6+24epsilon lambda-22epsilon/(3B^2)-28epsilon d/B^2,
    b=1-8epsilon lambda-epsilon(11/3+4d)/B^2,
    V_U=4lambda-2d/B^2,
    P_S=24lambda-4d/B^2+88epsilon lambda^2+6epsilon V_U^2,
    a(p)=c(p-2)+A^2 P_S=c(p-r0),
    r0=2-A^2 P_S/c.                                    (3)

In particular, retain the actual d and the full stationary-product P_S.
Setting epsilon=0 with fixed d is not the fixed-reference I0 replacement.
The parameters in (1) imply the deliberately weak, convenient guards

    5.9<c<6.1, .999<b<1.001, 1.3<r0<1.4.              (4)

For completeness, 24x/C differs from 4 by less than 1e-9, while all
epsilon contributions to c,b and A^2 P_S are less than 1e-29. These
follow by writing G lambda=eta_G x/C and G/B^2=eta_G/C in (3).
Thus A^2 P_S lies between 3.999999 and 4.000001, which is ample for (4).

The reduced cold-response theorem proves the product Ward factor and
the low-frequency bound

    q_cold(p)=(2G/B^2) Qhat_red(p), Qhat_red(0)=0,
    sup_(|p|<=.01) |Qhat_red(p)/p|<1e24C.               (5)

Only the resulting bound on Qhat_red'(0) is needed here; the small-p
bound is not extrapolated to a high-frequency contour.

## 2. Primary bubble, Green normalization and physical pairing

The primary source is Carmi, Di Pietro and Komatsu, A Study of Quantum
Field Theories in AdS at Finite Coupling, arXiv:1810.04185v3:
https://arxiv.org/html/1810.04185
Equations (6.1), (6.4), (6.7), (6.9), (6.25), (D.10), and Appendix B,
especially (B.3), (B.5), (B.8), (B.11) and (B.31), were directly read.
Its bubble is the spin trace of two inverse first-order Dirac kernels.
The per-flavor effective action is minus a complex Dirac trace log,
without a Majorana factor. The d=1 representation has two complex spin
components. Its standard positive-mass boundary condition has
Delta=M+1/2. We use only its subtracted bubble, not a Gross-Neveu gap
equation, coupling, or a separately chosen Pauli-Villars finite part.

Here is the normalization bridge. The orbit operator Dsharp retained in
the accepted construction is a rank-two anti-self-adjoint geometric
Dirac operator. Its rotated Clifford multiplication is parallel and is
unitarily equivalent to the ordinary complex rank-two Clifford action.
On simply connected complete H2 the spin structure is unique. The
unitary identification preserves the spin connection and the delta
source. Consequently Dsharp+M and the primary gamma.grad+M have the
same inverse kernel up to that unitary identification. There is no
choice of a second boundary condition: M>18 and the decaying inverse
are retained. The kernel product and its spin trace are invariant under
the unitary conjugation. This also agrees with the accepted local
Dirac singular coefficient 1/(2pi d), fixing the Green normalization.

There is one complex rank-two orbit kernel per angular mode. The d_n
already count both nonzero sphere-Dirac signs on one physical exterior.
There is no extra factor of two and no extra half in the primary bubble.
If one starts from the folded rank-four presentation, division by two
returns precisely this convention. Dividing its sphere-integrated orbit
Hessian by 4pi B^2 A^2 and then by gamma0 gives the factor 2G/B^2 in (5).

At H2 radius one the primary normalized spherical transform of an
invariant kernel K(d) is

    Ktilde(zeta)=2pi integral_0^infinity sinh d K(d)
                       P_(-1/2+i zeta)(cosh d) dd.       (6)

The primary harmonic eigenvalue is -Delta=zeta^2+1/4. Hence

    p=zeta^2+9/4.

At p=0 choose zeta=-3i/2, so the spherical function in (6) is
P_1(cosh d)=cosh d. The other sign gives P_-2=P_1. The convention
therefore agrees exactly with the Ward subtraction in the reduced cold-response theorem.
Convolution multipliers multiply without an additional factor at unit
radius; the physical A and B factors were accounted for above.

At fixed two-dimensional M, the bubble has scaling degree two at the
diagonal. Its extension ambiguity at that degree is a delta term, a
constant transform multiplier. The primary equations (6.25) and
(6.29) leave precisely this constant free. Taking the difference
between p and zero removes it. No derivative-delta ambiguity is
introduced at fixed M. The higher four-dimensional common angular
heat contact is retained explicitly in section 4 below.

## 3. Exact positive poles and their two needed moments

Define W(x)=Gamma(x+1/2)/Gamma(x). Substitution of d=1 and
Delta=M+1/2 into the primary (6.25), using the exact OPE coefficient
(D.10), and cancellation by the gamma duplication formula gives

    tau_j=(2M+2j)(2M+2j+3),
    A_j=(2/pi)(2M+2j+3/2)
          W(j+1)W(2M+j+1)/W(M+j+1)^2 >0.              (7)

One way to audit the cancellation is to remove the pole numerator
2(2M+2j+3/2) first. The remaining primary gamma coefficient reduces
to pi^-1 W(j+1)W(2M+j+1)/W(M+j+1)^2. This fixes the rank-sensitive
factor 2/pi in (7). At large j, A_j~4j/pi, consistently with the
rank-two flat bubble logarithmic coefficient 1/(2pi).

The primary pole term is -A_j/(tau_j+p). Consequently

    BF_M(p)-BF_M(0)=p sum_(j>=0) A_j/[tau_j(tau_j+p)],
    Pi_M(p)=M^2 p sum_j A_j/[tau_j(tau_j+p)]-p/(12pi).   (8)

These sums converge absolutely and locally uniformly off their poles:
A_j=O(j) and tau_j grows quadratically. The subtraction at p=0 is
below every first two-fermion pole. In particular tau_j>1404 for M>18.
Equation (8) is the invariant H2 bubble and the already retained
ordinary orbit anomaly, not a new choice of subtraction.

The exact elementary decomposition is

    Pi_M(p)=ell_M p-c_M p^2+p^3 T_M(p),
    ell_M=M^2 sum_j A_j/tau_j^2-1/(12pi),
    c_M=M^2 sum_j A_j/tau_j^3,
    T_M(p)=M^2 sum_j A_j/[tau_j^3(tau_j+p)].             (9)

In particular c_M is not replaced by its leading common heat
coefficient 1/(120pi M^2). We need an explicit bound on the difference.

Log convexity of Gamma gives, for x>=1,

    sqrt(x-1/2)<=x/sqrt(x+1/2)<=W(x)<=sqrt(x).

Put t=M+j. Applying these inequalities to (7) gives

    (4/pi)(t+3/4)/(t+1) sqrt((t+1/2)^2-M^2) <= A_j
     <=(4/pi)(t+3/4)/(t+1/2) sqrt((t+1)^2-M^2),
    A_j<2(t+1).                                       (10)

Let f_M(s)=sqrt(s^2-M^2)/s^6, s>=M. Its integral and maximum are

    I_M=integral_M^infinity f_M(s)ds=2/(15M^4),
    max f_M=125/(216sqrt(5)M^5).

It increases once and then decreases. Comparing on the unit intervals
[M+j,M+j+1], for every shift theta in [0,1],

    |sum_j f_M(M+j+theta)-I_M|
         <=TV(f_M)=2max f_M <(3.9/M) I_M.              (11)

For the upper bound in (10) use s=t+1. Then
tau=4s^2[1-1/(2s)-1/(2s^2)], and for M>=18 the bracket is at least
1-.53/M. The reciprocal cube is at most 1+1.8/M. The other ratio
(t+3/4)/(t+1/2) is at most 1+.25/M. Their product is at most
1+2.1/M. The reciprocal-cube bound follows, for example, by expanding
(1+1.8z)(1-.53z)^3-1 for 0<=z<=1/18; its coefficient after division
by z is greater than .09 on this interval.

For the lower bound use s=t+1/2. Now
tau=4(s-1/2)(s+1)<=4s^2(1+.5/M). Its reciprocal cube is at least
1-1.5/M, and (t+3/4)/(t+1)>=1-.25/M. Their product is at least
1-1.75/M. Combining these estimates with (11) and
M^2 I_M/(16pi)=1/(120pi M^2) proves

    1-5.65/M <=120pi M^2 c_M<=1+6.5/M,
    |c_M-1/(120pi M^2)|<=7/(120pi M^3).               (12)

The second moment needed for the high-frequency proof is much easier.
By (10), tau>4t^2, and the integral test,

    M^2 sum_j A_j/tau_j^4
     <(M^2/128) sum_j (t+1)/t^8
     <=(M^2/128)(19/18)[M^-7+1/(6M^6)]
     <=38/(20736M^4)<1/(500M^4).

Thus, using (2),

    sum_(n,j) d_n M_n^2 A_(n,j)/tau_(n,j)^4<3e-5C.     (13)

## 4. The exact restored quadratic coefficient is positive

The common angular finite part and reduced local counterterm are
exactly those of the reduced cold-response theorem:

    Qhat_red(p)=-q p/(12pi)
      -C[F_eff/(120pi)+1/(180pi)]p^2
      +sum_n d_n[Pi_Mn(p)+p^2/(120pi M_n^2)].            (14)

Equation (12) and the accepted summable slope remainder allow (14)
to be regrouped, with locally uniform convergence, as

    Qhat_red(p)=A_s p+B_s p^2+p^3 S(p),
    A_s=-q/(12pi)+sum_n d_n ell_Mn,
    B_s=-C F_eff/(120pi)-C/(180pi)
                  -sum_n d_n[c_Mn-1/(120pi M_n^2)],
    S(p)=sum_(n,j) d_n M_n^2 A_(n,j)/[tau_(n,j)^3(tau_(n,j)+p)].
                                                                 (15)

The sum defining S converges normally away from the negative pole ray
by (13). The summability of A_s is the fixed local common-subtracted
VALUE estimate already proved in the reduced cold-response theorem; it is not inferred from a
quadratic-form estimate. In particular |A_s|<1e24C by (5).

By (1) and (12), the absolute finite correction to B_s is below
7C/(180pi). Moreover H_q>log(q+1)>28 and gamma_E<1, so F_eff<-51.
For the elementary logarithm bound, e<11/4 and (11/4)^7<2000 imply
log(2000)>7. This is much weaker than its actual value. Consequently

    B_s > C[51/(120pi)-1/(180pi)-7/(180pi)]
         =137C/(360pi)>C/10,
    B_s < C[100/(120pi)+7/(180pi)]<C/3.                (16)

The upper estimate drops the negative reduction term. The inequalities
use only 3<pi<22/7. They include the exact finite curvature correction;
positivity has not been deduced from the leading a2 coefficient alone.

## 5. Positive Stieltjes structure of the physical denominator

Set eta=2G/B^2 and define

    alpha=eta A_s, |alpha|<1.7e-8,
    beta=eta B_s, eta_G/5<beta<2eta_G/3,
    dmu=sum_(n,j) eta d_n M_n^2 A_(n,j)/tau_(n,j)^3
                            times the point mass at tau_(n,j),
    H(p)=beta+p integral dmu(tau)/(tau+p),
    g(p)=q_cold(p)/p=alpha+p H(p).                      (17)

The positive measure is supported in tau>1404 and satisfies

    integral dmu(tau)/tau<6e-5 eta_G.                   (18)

Its total mass need not be finite. The weighted integrability in (18)
is exactly what (17) requires. H is analytic off the negative pole ray.
For real r>=0 it is positive and increasing. For p=x+iy with x>0,
each p/(tau+p) has argument between zero and arg p; hence H has that
same sector property. Direct calculation also gives

    Re H(p)>=H(x),
    H(r)<= (r/x) H(x) if r>=x>0.                       (19)

For the first inequality subtract x/(tau+x) from
Re[p/(tau+p)]; the difference is
tau*y^2/[(tau+x)((tau+x)^2+y^2)]>=0. The second inequality follows
term by term and includes the positive constant beta.

The exact coupled determinant is

    D_r(p)=4b^2 p+a(p)g(p),
    det M(p)=-p D_r(p),
    D_r(p)/p=4b^2+a(p)H(p)+a(p)alpha/p.                (20)

All signs in (20) come from the original two by two coupled matrix.
No scalar-alone radius equation or commuting Schur approximation is
being substituted for it.

## 6. Every frequency on both growth contours

Choose w=1+.01 or 1-.01 and xi real. On either contour put

    nu_growth=w+i xi,
    p=2-nu_growth(nu_growth+1)=xi^2-delta-i kappa xi,
    delta=w(w+1)-2 in {.0301,-.0299},
    kappa=2w+1 in {3.02,2.98}, r=|p|.                  (21)

Since kappa^2-2delta>0, r has its minimum at xi=0. Thus

    r>=.0299>1/34.                                     (22)

First suppose |xi|<=4. Then r<21 and Re p>=-.0301. From (17)-(18),

    |H(p)| <= beta+[21/(1-21/1404)] integral dmu/tau
             <eta_G,
    |a(p)|<6.1(21+1.4)<137,
    |a(p)alpha/p|
        <=6.1(1+1.4/.0299)(1.7e-8)<5e-6.

It follows that |D_r(p)/p|>3.99. The same bound (18) gives
H(r)<eta_G in this range. Hence in particular

    |D_r(p)|>=r[3+.75r H(r)]                           (23)

for the entire compact-frequency part. The comparison here is between
3.99 and 3+21eta_G, not an unsupported use of a complex sign.

For |xi|>=4 write p=x+iy. Now x>15.9699 and (4),(21) give

    |y|/x<.76, |y|/(x-r0)<.83,
    5r<|a(p)|<7r.

Both p and p-r0 lie in the same open right half-plane. Their argument
sum has cosine at least

    (1-.76*.83)/sqrt((1+.76^2)(1+.83^2))>1/5.

The sector property of H therefore implies

    Re[a(p)H(p)]> |a(p)||H(p)|/5 >r|H(p)|.

Using (20), 4b^2>3.992 and |a alpha/p|<7|alpha|,

    |D_r(p)/p|>=Re[D_r(p)/p]>3+r|H(p)|.

Also r/x<sqrt(1+.76^2)<1.26. Equation (19) yields
|H(p)|>=H(x)>H(r)/1.26>.75H(r), proving (23) at every remaining
frequency. Thus the complete two-contour result is

    |D_r(p)|>=r[3+.75r H(r)]
             >=.75r[1+r H(r)]>0.                       (24)

The two exact roots identified near p=0 in the reduced cold-response theorem (the Ward root
and the physical coupled root, counting multiplicity) lie between
these contours. Equation (24) excludes any additional intersection
of either contour with a determinant zero at arbitrarily high frequency.
It does not assert that either root is an admissible gauge mode at an
actual exterior seam.

## 7. Quantitative inverse and its precise contour domain

Let W_r=1+rH(r). On the two contours one also has

    |H(p)|<=2H(r).                                     (25)

For |xi|>=4 this follows termwise from
|p/(tau+p)|<=r/(tau+x) and (tau+r)/(tau+x)<=r/x<1.26.
In the compact-frequency part use instead
(tau+r)/|tau+p|< (tau+21)/(tau-.0301)<1.02.
The beta term is preserved in both comparisons.

The inverse matrix entries are exactly

    (M^-1)_SS=-g(p)/D_r(p),
    (M^-1)_ST=(M^-1)_TS=-2b/D_r(p),
    (M^-1)_TT=a(p)/[p D_r(p)].                         (26)

Equations (4),(17),(24),(25) imply the useful separate bounds

    |(M^-1)_SS|<=3/r,
    |(M^-1)_ST|<=3/(r W_r),
    |(M^-1)_TT|<=12(1+r)/(r^2 W_r).                   (27)

The last coefficient uses 6.1(r+1.4)/.75<12(1+r).
By (22), the bounds on the entries are respectively 102,102,14280.
The maximum absolute row and column sums are below 14382. Therefore

    ||M(p)^-1||_(C^2->C^2)<16000                      (28)

uniformly on both contours, with all high-frequency weights in (27)
retained. In physical normalization, multiply the inverse by
gamma0^-1=8pi G A^2.

Here is an exact functional statement which does not presume a
shifted-transform theorem. On either line let L2_C be two-component
L2(dxi), or L2 with any fixed positive measurable density in xi.
Define the maximal multiplication domain

    D_C={U in L2_C : M(p(xi))U(xi) is in L2_C}.

Then M:D_C->L2_C is bijective, and its inverse obeys (27)-(28).
The graph norm satisfies ||U||+||MU||<=16001||F|| for MU=F.
This is a complete full-frequency contour-multiplier inverse; no
frequency interval and no unestimated large-frequency remainder is
omitted. It does not identify D_C with a weighted spatial bolt domain.

## 8. An actual complete-product L2 inverse

There is also a spatial result requiring no shifted contour. The
ordinary spectral resolution of the self-adjoint -Delta_H2 has
p=9/4+xi^2, xi real. On that spectrum p is positive, p>r0 and

    a(p)/p=c(1-r0/p)>5.9(1-1.4/(9/4))>2,
    |a(p)/p|<6.1.

The real positive H(p) and (20) therefore give

    D_r(p)/p>3+2pH(p)>3+.75pH(p).

Bounds (24)-(28) hold on this spectrum as well (with substantially
more room). Apply the ordinary H2 spectral theorem componentwise.
The real symmetric matrix M(p) defines the self-adjoint maximal
spectral operator with domain {U in L2(H2;C^2):M(P2)U in L2}.
Its inverse is the bounded everywhere-defined operator M(P2)^-1.
This remains true on the rotationally invariant subspace, with its
regular H2 origin, and with the fixed physical constant pairing.
No bolt extension parameter is introduced: this is the complete H2
Laplacian and its usual regular-origin domain. The inverse obeys

    ||H_red^-1||_(L2->L2)<16000(8pi G A^2).            (29)

The maximal domain description is intentional: the exact function
H retains the positive nonlocal pole sum and the exact common heat
contact. It has not been replaced by a constant fourth-order model.

## 9. What this closes, and the specific remaining transfer

The full product frequency-contour calculation no longer requires an
unknown bubble bound. Equations (7)-(8) prove the rank-correct positive
pole representation; (12),(15)-(16) control the exact restored p^2
coefficient; and (24),(27) give explicit all-frequency inverses for
the coupled radius and orbit shape. The low physical root and Ward
factor are retained. Radius, shape and physical pairing were never
decoupled in this argument.

The companion [weighted response proof](stationary-product-weighted-response.md)
now reconstructs the actual complete-H2 inverse through its resolvent,
retaining the spherical density, regular bolt normalization and the
whole crossed pole cluster. For sources supported in a cold ball of
radius at most 84, it gives explicit cosh-weighted C0, L2 and radial
C2 bounds below 1.3e8 in the normalized equations. The source cutoff
jets and complete-H2 exterior tail are part of that statement.

The remaining transfer must bound the actual variable-radius Hessian,
including its quantum nonlocal kernel, relative to this product
particular inverse. It must impose the fixed physical mass/temperature
return constraints at both matching seams. The actual finite-throat
seed is not globally the stationary product. The local VALUE-response
theorem does not by itself estimate that operator difference in the
required inverse norm.

No actual exterior solution, original signed-source compatibility,
nonlinear branch or physical error control follows from these product
estimates. All E003/S074-S080 physical gates remain open. No numerical
bubble or contour solver was used to derive the ordinary inequalities.

The shared [focused checker](../checks/stationary_product_weighted_inverse.py)
audits finite identities and rational reserves for the contour and
weighted-response proofs. Its [JSON receipt](../evidence/stationary-product-weighted-inverse.json)
and [raw log](../evidence/stationary-product-weighted-inverse.log) record
the actual execution outcome and exact input hashes. Analytic
continuation, infinite pole sums, spectral domains, contour deformation
and spatial reconstruction remain ordinary proofs.

Remaining to-do list: use the companion weighted reconstruction in the
actual coupled quantum seam transfer and fixed mass/temperature return;
prove the nonlinear source-to-observable and physical error estimates.
