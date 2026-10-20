# Complete retained source response on the comparison path

2026-10-20, E003/S160 supporting S139/S074. Ordinary argument and
predeclared verification before acceptance.

This completes the retained Dirac-minus-reference density, null and
angular SOURCE response on the specified comparison path. It combines
S155, S156, S158 and S159 with the remaining Weyl cold/middle stress
and the independent local anomaly/radius allocation. It does not
complete the E0 operator, a general-family source theorem, full
higher-order inverse, coupled solution or physical stability proof.

## 1. Independent radius variation and the exact trace allocation

Use the retained spherical action and chain rule in
[the anomaly source derivation](spherical-anomaly-sources.md).
Its underlying Weyl transgression is
[Mazur--Mottola](https://arxiv.org/html/hep-th/0106151v1), inspected
2026-10-20. That paper does not supply the bounds below.
No finite coefficient, scale, flux, state or boundary ensemble changes.

Use dimensionless product-frame derivatives in A units and set

    e=1/C, qbar=q^2/C^2<160000, eta_G=GC/B^2<8.5e-33,
    phi=log(r/B), aA=11/360, cA=1/20, bA=1/30,
    dA=bA-2aA/3,
    R=R_product, U=Box phi, V=(grad phi)^2,
    Q=bA R-6bA U-6dA V,
    Z=bA R+(2bA-4aA)e-6dA(U+V),
    E=cA(R+2e)^2/3-4aA e R-qbar/3
                            +Box Q-2 div(Z grad phi). (1)

E is the independent Weyl-field Euler derivative. It is not
replaced by the orbit trace. Let M_t,M_r be the mixed local
Euler components in S158. Their normalized sources are

    W_rho=eta_G M_t/pi, W_N=eta_G(M_t-M_r)/pi,
    S_W=-eta_G(M_t+M_r)/pi,
    A_W=eta_G E/(4pi)-S_W/2.                     (2)

For any physical source use

    rho=(2G/(B^2 C))r^2 rho2,
    N=(2G/(B^2 C))r^2(rho2+p2),
    A=(8piG/(B^2 C))r^4 pS.

The canonical massive/LLL contribution has A_can=-S_can/2.
Thus the full Dirac-minus-reference source satisfies exactly

    A_rest=rho_rest-N_rest/2+T,
    T=eta_G E/(4pi)
                  -K[2pv t^4+DP g(t)],
    K=8piG B^2/C, g(t)=-2/t+3/t^2, t=r/B.        (3)

Indeed the chain rule gives
J_X=(2G/X)exp(-2phi)(s_W-tau_hat,total), and pS=J_X/(16piG).
Also S_total=N_total-2rho_total and
t^4(pA_ref+pS_ref)=2pv t^4+DP g(t). This derivation keeps
the independent radius source; it does not differentiate an
uncontrolled spatial stress profile or infer the angular equation
from solved orbit constraints.

At the exact product, phi=0, R=-2, and all its covariant
derivatives vanish. The Weyl ORBIT stress is zero, but

    E_*=(1+e^2)/15+e/9-qbar/3,
    eta_G E_* /(4pi)=K(p0A+p0S)=K(2pv+DP).       (4)

The angular product value is generally nonzero. Consequently
the cold local allocation is exactly

    T=eta_G(E-E_*)/(4pi)
           -K[2pv(t^4-1)+DP(g(t)-1)].            (5)

It vanishes at the product before a norm is taken.

At the outer boundary use the actual proper coordinate and
correlated jets from S158,

    w=eY^2/2+cY^3+O_4(Y^4 L),
    phi=-log(dY)-eY^2/6-cY^3/4+O_4(Y^4 L).

Substitution in (1) gives the negative powers

    E=-11/(15Y^4)+22e/(45Y^2)+11c/(15Y)+O(L).

They are exactly those of
E_vac=-11d^4 exp(4phi)/15. The fixed vacuum identity is
eta_G E_vac/(4pi)=2K pv t^4. Thus

    T=eta_G(E-E_vac)/(4pi)-K DP g(t)              (6)

is integrable in the outer squared-log profile norm.
There is no separate mass cubic in T. For an independent phi
cubic f3 the 1/Y coefficient of E itself is still 11c/15;
subtracting its actual-radius vacuum leaves (11c+44f3)/(15Y).
The correlation f3=-c/4 is essential.

## 2. Finite conformal jet algebra and coefficient reserves

At a fixed real observation use conformal Cartesian coordinates,
g_product=exp(2sigma)(dx^2+dy^2). Treat H=exp(-2sigma) as
a formal coefficient with D_i H=-2H sigma_i. Then

    R=-2H(sigma_xx+sigma_yy),
    U=H(phi_xx+phi_yy), V=H(phi_x^2+phi_y^2).

Insert these in the retained action L,P,Q,Z and

    M_ij^mixed=-delta_ij L/2+R delta_ij P/2
       +delta_ij Box P-H Hess_ij P
       +H Z phi_i phi_j-H phi_(i Q_j)
       +delta_ij[H grad Q.grad phi+QU]/2.        (7)

The Hessian uses the full conformal connection. In particular
Hess_xx P=P_xx-sigma_x P_x+sigma_y P_y and
Hess_yy P=P_yy+sigma_x P_x-sigma_y P_y.
At a radial observation y is the angular/time direction.
Formula (1) is expanded in the same variables.

The resulting polynomials have total derivative order at most four,
at most one undifferentiated phi and H power at most two.
Separate the algebraic magnetic monomials M_t=qbar phi/6
and E=-qbar/3. Count absolute coefficients times
max(1, total field-jet degree+2 degree_H); the 2 degree_H
accounts for varying H. For M_t, M_t-M_r and E the count
is below 1e4. After one Cartesian derivative each analogous
count is below 1e5, with total derivative order at most five.
The checker expands the full two-variable formulas to verify
these finite counts; it does not assume a proper-coordinate
formula is regular at the bolt.

Background differentiated jets are bounded by j!100^j. The
undifferentiated phi bound is one in the cold region and four
in the middle. Products at total derivative order four or five
are bounded by 4!100^4 or 5!100^5. Along the cold interpolations
below |H|^k<2 for k<=2. Therefore a weighted field difference
of 100 a_x j!100^j gives regular polynomial differences

    order four: 2*1e4*100*(4!100^4)a_x,
    order five: 2*1e5*100*(5!100^5)a_x.          (8)

The magnetic terms are counted at derivative order zero/one;
the constant in E cancels from E-E_*.

## 3. Actual cold and middle charts

All unknown spatial arguments remain real. Allow |lambda|<=33
on the original path sigma_lambda=sigma0+lambda u,
r_lambda=B+(r0-B)exp(lambda V). Only finite jets are used.

For cold points s=t0-1<=.001, put a_x=s[1+log(.01/s)]<=.00331.
Use a fixed H2 normal conformal chart and compare with the exact
H2 x S2 product in that chart. The seed cold functions are analytic.
Here are explicit margins sufficient for their weighted jets.
In the transvection-centered bidisk |z|,|w|<=1/64,
the fractional change of cosh chi is below .017. The retained
inverse-s bound of two gives |s_complex-s0|<.04s0, within the
quarter s-disk. The exact cold q0,D formulas give
|sigma_seed-sigma_H2|<20a_x and |phi_seed|<3s0 there.
Changing to Cartesian complex coordinates halves that radius.
Cauchy therefore bounds their normalized j-jets, j<=5, by
20(1.28)^j a_x and 3(1.28)^j a_x.

The same seed normalizer and its inverse have a finite
composition cost below four through degree five. For example
use inverse linear coefficient 1.1, higher coefficients bounded
by 2 Cat_(j-1)(.3)^(j-1)/(.9)^(2j-1), and divide the
composition coefficients by 100^j. The once-made endpoint
transfer of S155/S159, with cold weighted low jets, then gives

    u: normalized jets <=4e-6 a_x,
    nu=tau V: normalized jets <=4e-5 a_x, tau=s/t0.

This is a finite identification of the actual tangent, not
parameter differentiation of a moving normalizer.

For completeness finite radius exponentiation preserves this
ideal also at complex lambda. In the larger transvection
bidisk |z|,|w|<=.1, the same inverse-s bound gives relative
s change <.22. The analytic inverse-tau normalized jet sum
through five is below 2/s, using the Cartesian radius .05
and the geometric sum of (20/100)^j. Thus lambda V has
normalized high jets at most A=33(8e-5)L_c<=.214.
Its scalar value remains below 4e-7. The finite exponential
coefficients through five are bounded by A(1+A)^4.
The tau product and logarithm in
log[1+tau(exp(lambda V)-1)] consequently cost less than
.02 a_x. The scalar logarithm denominator stays above .999999.
Together with the seed bounds and affine orbit path this proves

    |partial^alpha(delta sigma,phi)|
              <=100 a_x j!100^j, j<=5,           (9)

relative to H2 and phi=0. The sharper values give |delta sigma|<.05,
|phi|<.002. Background differentiated jets on the straight
interpolation to the product obey j!100^j, since
100a_x<.331 and the H2 normal metric coefficients are smaller
than the remaining margin.

On S152's middle interval, 1.001<=t0<=40, S155's fixed
normal-chart orbit bounds apply. The seed radius log has
complex modulus <6 on its 1/32 normal bidisk. One can check
this directly from z_Y^2=calB(z)/C: on the real strip
|calB'|/C<.3, |calB''|/C<=2, and in the .02 complex strip
|z_Y|<.5. The normal-coordinate/Jacobi recurrence with
complex |R|<3 gives |Y-Y0|<.035. Since z0>=1/40, the
relative z change is <.75; hence |log z|<log40+log4<6.
The low seed radial gradient is <20.1. For j>=2 the
normalized seed bound is 6(.32)^j.

Repeat S155's inverse-quadratic scalar transfer with radius
endpoint bound 1e-4 and seed modulus 6. Through degree four
the normalized nu jets at orders two, three and four are bounded
by .000342, .000796 and .001854. The nonlinear radius path needs
an additional finite composition; it is not affine in phi.

Put tau0=(t0-1)/t0>=1/1001. On the seed normal bidisk of radius
1/64, both tau/tau0 and its inverse have modulus below two.
Indeed tau=1-z and the factorization
calB/C=(1-z)^2[x(1+2z+3z^2)+z^2+d z^4]/C-delta z^3/C
gives |partial_Y log tau|<1.2 on its continuation: the bracket
is below 1.1 on the .02 complex z-strip, while the thermal
fraction is below 1e-50 for tau0>=1/1001. The normal-coordinate
displacement there is below .018, so |log(tau/tau0)|<.022.
Cauchy in the Cartesian radius 1/128 bounds each normalized
positive-degree jet of tau and 1/tau by
2 tau0(1.28)^j and 2 tau0^-1(1.28)^j, respectively.

The following positive finite series retain all cross terms
through degree four; their indeterminate measures 100 times
the sum of Cartesian increments:

    n(X)=1.2e-8+6.3e-8 X+.000342 X^2
                                  +.000796 X^3+.001854 X^4,
    J(X)=1+2 sum_(j=1)^4 (1.28 X)^j,
    v_j=33 [X^j](n J), j>=1.

The nonconstant jets of lambda V are bounded by v_j/tau0;
its constant value is below 4e-7. Set tau_min=1/1001 and
q(X)=sum_(j=1)^4 v_j X^j/tau_min. For positive degrees the
finite coefficients of

    H(X)=4e-7 J(X)+tau_min J(X) exp(4e-7)[exp(q(X))-1]

majorize tau[exp(lambda V)-1] uniformly over tau0>=tau_min.
Each term containing k differentiated V factors scales as
tau0^(1-k), so its maximum occurs at tau_min. The separate
constant term uses tau0<=1. Finally expand
sum_(k=1)^4 [(H-H(0))/(1-4e-7)]^k/k to bound the logarithm.
The resulting radius-change coefficients at orders two through
four are below .012, .085 and .42. Adding 6(.32)^j stays
below one at each order. At first order the sharper real
gradient identity retains the bound 33(6.3e-6), up to the
factor exp(4e-7)/(1-4e-7) and an additional 1e-12 for
[exp(lambda V)-1-lambda V exp(lambda V)]grad tau.
The seed gradient is below 20.1.
This proves the middle differentiated-jet guard for the actual
exponential radius path. The affine orbit jets also stay below one.
The point values obey |phi|<4
and |sigma_lambda|<4e-7 in the fixed seed normal chart.
Thus the middle coefficient bounds in section 2 apply on the
entire complex parameter disk.

## 4. Local weighted response bounds

Equations (8)-(9), the factors eta_G/pi and eta_G/(4pi),
and qbar<160000 give the following joint complex bounds:

    cold: max(|W_rho|,|W_N|,eta_G|E-E_*|/(4pi))
                                                   <1e-15 a_x,
    middle: max(|W_rho|,|W_N|,eta_G|E|/(4pi))<1e-16.

The Weyl orbit source vanishes on phi=0 for every metric,
so its cold comparison has no product constant. Rotational
regularity gives W_N(0)=0 for every complex parameter.
Use the order-five bound along the radial segment, where
a_x<=2h L_h and seed proper length<1.1chi. It gives

    |W_N|<1e-11 h L_h chi
                         <=1e-11 h L_h sqrt(chi). (10)

The fixed-chart polynomial argument establishes holomorphy
in these weighted spaces, including the bolt modulus.
Cauchy circles of radius 32 give first bounds .032B and
second bounds .002B from each complex value bound B.

For the local T in (3), the reference coefficients obey
K q^2/B^4<1.71e-26, |pv|B^4/q^2<1, |DP|B^4/q^2<3.
Use (5) in the cold region, (3) in the middle and (6) outside.
In the outer estimate use S158's proper remainder bound
100 j!100^j Y^(4-j)L through order four. The polynomial E
has no undifferentiated phi. Its variation and the explicit
vacuum exponential difference lose no more than L^2.
On the truncated reference jets, the remaining nonnegative
power coefficient majorant is below 1e8(1+qbar); the
vacuum exponential remainder is bounded by its convergent
Taylor estimate. Normalization puts that reference remainder
below 1e-15 L^2. The proper-distance and radius responses
remain in these formulas; they are not frozen.

Sufficient T bounds are:

| Region/weight | Complex value | First response | Second response |
| --- | ---: | ---: | ---: |
| cold / a_x | 1e-14 | 3.2e-16 | 2e-17 |
| middle / 1 | 1e-15 | 3.2e-17 | 2e-18 |
| outer / L_e^2 | 1e-12 | 3.2e-14 | 2e-15 |

All homogeneous and vacuum cancellations precede these norms.

## 5. Complete source assembly and its limits

Apply S152 to the new Weyl cold/middle profiles, with
square-root bolt column 486, d_g=3e-7 and the full product
rules for changing geometry. Its second-response contributions
to mass primitive, absolute mass density /BC and lapse L1 are

    (8e-17, 1.1e-19, 5e-12).                    (11)

Adding S155, S156, S158 and S159 separately, with their
distinct regional/logarithmic columns, gives for the COMPLETE
retained Dirac-minus-reference source:

| Constraint density output | Second response |
| --- | ---: |
| zero-datum mass primitive | 7.21e-8 |
| absolute mass density /BC | 7.41e-11 |
| absolute lapse density L1 | 3.61e-5 |

Taylor allowances are half these values times lambda^2.
The cold massive allowance still dominates. These are source-map
bounds only, not an E0 residual, inverse estimate or solution error.

For the full radius source, (3) and the now complete orbit
profiles give the following angular bounds:

| Region and respective weights | Value | First | Second |
| --- | ---: | ---: | ---: |
| cold, a_x throughout | 1.7e-6 | 2.6e-10 | 4.7e-9 |
| middle, unweighted | 2.6e-9 | 2.2e-12 | 4.7e-14 |
| outer, L_e^2 / L_e^4 / L_e^6 | 4.7e-7 | 4.9e-11 | 3.2e-12 |

The outer source is C2 in the slightly weaker L_e^7 space,
as in S157. No angular source spatial derivative is inferred
from these parameter bounds. The radius equation remains a
required equation of any eventual solution.

At the actual moving bolt the force ratio
8piG pS/(2kappa/r_lambda)=A_rest/[2gamma t_lambda^3],
gamma=kappa B/C>.99h. The known t_lambda derivative
factors are retained. With L_h<81, first/second angular-force
responses are below 1.1e-8 and 2e-7; the corresponding
orbit-pressure responses are below 7.1e-9 and 1.3e-7.
These coarse derivative bounds do not replace the much sharper
fixed-seed forcing value, set the thermal gap, or solve its return
condition.

This theorem concerns one specified affine common-chart path.
It supplies no preservation theorem for iterates of an unknown
solver. The E0 operator/complete residual assembly and a usable
full higher-order inverse remain necessary next steps.

## 6. Predeclared verification

Revalidate S155-S159 and S152 receipts and exact current inputs.
Expand the two-variable conformal action/stress/radius polynomials;
check their derivative order, H powers and coefficient reserves,
both metric Euler variations, scalar Euler variation and Ward
identity. Independently reduce to proper radial jets. Check product,
angular/reference and outer Laurent cancellations, finite chart
majorants, all regional/constraint budgets and moving-bolt factors.
Preserve every run and any failed source snapshots.

Use 70 digits. For radial proper-chart diagnostics take

    w0=.001 y^2/2+.02 y^3,
    phi0=-log(.41y)-.001 y^2/6-.02 y^3/4,
    u=.001 y^3+.002 y^4(1-log y),
    w_lambda=w0+lambda u,
    phi_lambda=phi0-.0005lambda y^3
                             +.003lambda y^4(1-log y).

The proper derivative is exp(-lambda u)partial_y.
At y=.002,.02,.06 and lambda=0,.4,1 compare the direct
radial Euler expressions and the conformal-Cartesian formula
after a one-dimensional conformal-coordinate conversion.
Compare eta-free Weyl density, null and E-minus-vacuum second
responses from symbolic differentiation with 64/128-node Cauchy
contours of radius .01 and real steps .001,.0005. Use
absolute tolerance 1e-6(1+|target|), refinement ratio .35,
or floor 1e-35 when both errors are below it. Use one-sided
second-order stencils at the parameter endpoints. Formula
comparisons have tolerance 1e-30.

Also test a regular Cartesian bolt fixture
sigma=.02(x^2+y^2)+lambda*.003(x^2+y^2),
phi=.01+.03(x^2+y^2)+lambda*.002(x^2+y^2).
Require isotropic stress at the origin, a finite nonzero angular
product value and the local Ward identity. Check a manufactured
complete trace/radius combination and require that omitting E_*
from the cold cancellation or omitting the vacuum correlation
leaves a detectable error.

These fixtures test local algebra and parameter calculus, not
continuum determinants, a physical state, a nonlinear branch or
stability. Uniform spatial/weighted arguments remain ordinary proofs.

Remaining to-do list: E0 and complete residual assembly; general
admissible family and usable full inverse; coupled exterior/state,
charge, return and signed-source solve; opening, stability, transport,
preparation, separate creation and full-result audit.
