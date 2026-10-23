# Current reference regularity and the closed feasibility audit

2026-10-23, E003/S169. Ordinary argument frozen before execution.
The full stable semiclassical wormhole objective remains incomplete.
This calculation decides one open regularity requirement of the selected
actual-direction method. It does not construct the unknown correction.
The source-value estimate chain remains closed.

## 1. Scope and the quantitative test

Keep the current c_R=1e22 action, original mu0=1/b0, full Dirac tower,
temperature, product enclosure and finite paired matter domain. The
direction is the already constructed original-G reference h1, with
orbit component u. Work only on J=[.49,.88] in z=B/r. The observations
[.5,.875] and their .01 proper collars are inside J. No conclusion
here supplies the bolt, outer seam, unknown correction or nonlinear ball.

The preceding [diagonal audit](current-seed-diagonal-restoration-contract.md)
gave H5<2e8 and the principal kernel coefficient
kappa5=4GC/(5pi B^2). At delta=1e-5, a principal differentiated integration
error <=1e-28 needs

    [D_rho^5 u]_(1/2) <= 6e-28/(kappa5 sqrt(delta)),

approximately 1.2716273e8. This is a numerical component target, not a
solution tolerance. We derive a CURRENT H6 bound. The elementary inequality
min(2H5,H6|t|)<=sqrt(2H5 H6)|t|^(1/2) then gives the required modulus.
Alternatively the Lipschitz remainder is bounded by kappa5 H6 delta/12.
The executed receipt decides whether these bounds meet the frozen target.
An insufficient bound would not be a lower bound on the actual modulus.

Also substitute the newly obtained H5 into the PREVIOUS all-mode value
majorant, rather than keeping its old 2e8 direction allowance:

    E_value(delta)<=H5/60 sum_(j=2)^5 kappa_j delta^(6-j)/(6-j).

Test both delta=1e-5 and the already permitted .01 collar against 1e-28,
simultaneously with the principal Lipschitz derivative error. This is
a budget substitution using the same proved kernel majorant, not an
extension of the source-estimate chain. If .01 passes, the old 1e-5
resolution diagnostic is no longer necessary for these components.
The 32-attenuation mode scale still is not a certified truncation tolerance.
Remaining local moments and the full differentiated remainder are unchanged.

## 2. Fourth derivative of the current fixed-cap source

Extend [the fixed-cap C3 proof](fixed-cap-seed-source-three.md), using
exactly the same real cap and fixed gluing commutators. The degree-nine
cutoff is C4. Its ordered Cartesian derivative bounds through four are
315,7e5,3e9,4e12; the analytic seed bounds are
.002,.067,3,100,1000. Thus for phi=chi sigma,

    T0=.002, T1=.697, T2=1445.21, T3=6143635,
    T4=1000+4*315*100+6*7e5*3+4*3e9*.067+4e12*.002.

For a real unit direction e, use the GLOBAL L-infinity polynomial
p_e(w)=sum_(j=0)^4 w^j partial_e^j phi/j! on |w|=rho=1/500.
This is a holomorphic parameter polynomial, not a complex translation
of a nonanalytic cutoff. Its four real derivatives agree with translations
at zero; the C4 cap gives the required fourth-order uniform Taylor jet.
It has global magnitude <.0205 and center magnitude <7e-6. Therefore
|exp(p-p(0))-1|<.021<.06; the centered mass has real part >.99999M
and the global multiplier has real part >.979M. The same centered
Born inverse converges. For r>=1/80, .021<2r. On the inner disk the gradient
change is bounded by

    3rho+100rho^2/2+1000rho^3/6
       +(.002*5!*(96sqrt(2))^5)rho^4/24 <.02.

The original .067 gradient therefore still gives the local linear
potential bound 2r. Inner divided jets k=2,...,5 change by at most

    .002(.96)^k sum_(j=1)^4 binom(k+j,j)(96sqrt(2)rho)^j <.011.

These are exactly the local remainder inputs in the retained plane Born
proof. Its vertices are bounded multiplication potentials; a small
global gradient is not an assumption. The global cap needs four finite
translation derivatives, while the inner fifth spatial jets come from
the uncapped analytic seed. There is no claim that the cap is C5.
Uniform subtracted Born and angular convergence justify differentiation;
the fixed vacuum value is retained separately. Thus the same uniform
nonconstant source bound 2.1e12 q^2 holds on this parameter circle.

Cauchy and real polarization give, now through j=4,

    L_j=2^ceil(j/2) j^j(2.1e12)/(j!100^j rho^j).
    c0=10, c_j=180(j+3)[600(j+3)]^j,
    G_j=8*256^2*88^2*(31/2500000)
          *sum_k binom(j,k)c_k c_(j-k)/(j!100^j).

The fixed-cap off-diagonal row recurrence is checked also at j=4.
It uses the same real seed local coefficients and the actual paired
inverse <=1/M; it does not replace that domain by the plane cap.
For j>=1 the coordinate scalar jets are bounded by
(G_j+L_j+1e12) j!100^j q^2; use the retained 1e14 q^2 point bound at j=0.
The 1e12 reserve includes the explicit LLL curvature jets. At a normalized
center sigma=sigma_x=0, proper radial differentiation gives

    Df=f_x, D^2f=f_xx, D^3f=f_xxx-sigma_xx f_x,
    D^4f=f_xxxx-4sigma_xx f_xx-sigma_xxx f_x.

Use the real seed bounds 3 and 100 for the last two coefficients.
Multiplication by 2 eta_G(q/C)^2, eta_G<8.5e-33, q/C<400,
gives proper source jets S_can through four in their ORIGINAL units.
The raw fourth derivative allowance, not a renamed norm, is used below.

## 3. Complete local terms and exact reconstruction

The physical source is S=2G r^2(p-rho)/(B^2 C), not just S_can.
The R2/reference terms are combined algebraically before estimating,
as in [the complete cold defect](current-seed-cold-action-defect.md).
Its exact polynomial coefficient allowance is <1.16e-22 on |z|<=1;
the vacuum term is <2e-35 z^-4. The trace polynomial T=r^2 T0/(2C)
has complex magnitude <1e-24 on |z|<=1, where T0=8piG Arest.
The checker rebuilds these polynomials and their coefficient majorants.

For completeness the Weyl source is rebuilt from the compact-variation
Lagrangian in [the radial source derivation](stationary-weyl-radial-sources.md).
Put f=F/C=Q/z^2, R=-Q'', Tn=-z f_z, a=1/C, m=q/C and
sigma=-log(z)+log(B/b0). Its normalized bracket is

    H=sigma[R^2/60-a R/18+a^2/15-m^2/3]
       +R Tn/30-Tn^2/10+11Tn f/90-11f^2/180-a f/18,
    Hr=sigma(R/30-a/18)+Tn/30,
    Ht=R/30-Tn/5+11f/90,
    E=H-R Hr-(z^2 f_z+2zf)Hr_z-2z^2f Hr_zz
       +zf Ht_z+(zf_z+f)Ht,
    Ef=[sigma t^2 f_tttt+(4sigma+2)t f_ttt+2f_tt-2f_t/t]/30,
    S_W=(eta_G/pi)[E-f Ef/z^2], t=1/z, partial_t=-z^2 partial_z.

This includes the magnetic potential and the original scale. It is a
bulk interior source identity, not a boundary-charge prescription.
For each z0 in J the complex disk |z-z0|<=.04 has |z|>=.45, |z|<=.92,
|log z|<1 and the retained scale displacement <1e-8. Laurent/log
coefficient sums therefore bound S_W on that disk; no unknown metric
regularity is assumed. The checker verifies the generous total local
allowance M_S=1e-20 and trace allowance M_T=1e-24. Cauchy gives
coordinate derivatives <=j! M/.04^j through four. Exact chain rules
with z'=-sqrt(Q) convert these to proper derivatives and add them to
S_can. This supplies the complete S and T needed by the reference.

Write U=u', nu=xi/r=(1-z)V, N=nu', fcal=delta F/C and P=Gp/C.
Primes in the following system mean canonical proper derivatives.
Define H0=-Q_z/(2sqrt(Q)), L=sqrt(Q)/z and
rpoly=-12beta+epsilon(4z^5-6z^6)=B^2 R4/C. The exact equations are

    u'=U,
    U'=-H0 U+Q_zz u+3fcal-1.5z^2 S+T-.5z Q_zzz nu,
    nu'=N,
    N'=-(H0+2L)N+2(2Q/z^2-Q_z/z)u-fcal+.5z^2 S
        +[z Q_zzz/6+rpoly_z/(6z)-rpoly/(3z^2)]nu,
    fcal'=-L fcal-2L P+L z^2 S,
    P'=-2(H0+L)P+.5H0 z^2 S+L T.

These follow from the invariant curvature/radius equations of
[the reference transfer](stationary-completed-metric-response-norm.md),
the linear mass constraint and the FULL Ward identity. In particular
delta Rhat in canonical units is -2T-6fcal+3z^2S+nu z Q_zzz;
the trace forcing cancels in the radius equation. The nonconstant R4
Lie term is retained. The lapse/pressure Ward terms are not dropped.

The current, already proved low bounds on J are

    |u|<=1.2e-8, |U|<=7.7e-6, |nu|<=1.2e-8,
    |N|<=3e-7, |fcal|<=8e-10, |P|<=1e-10.

For N, use N=L[(t-1)V_t+zV], L<sqrt(.2)/.49, and
|(t-1)V_t|<3e-7. For fcal, use |fbar/F|<7.5e-10,
F/C<=.2/.49^2 and the existing tiny horizon shift. For P, integrate
Qp<1.3e-10 from the bolt with p=-rho there:
|P|< (1/5.9)*(.65e-10)(.49^-2-1)+negligible horizon term <1e-10.
These bounds refer to the original-G reference, not an iteration of
the full quantum Hessian.

On each interval cell, let w=(u,U,nu,N,fcal,P) and w'=Aw+bS+dT.
Use outward interval Taylor arithmetic for all coefficient derivatives
through four. With W_k componentwise absolute derivative bounds,

    W_(k+1)<=sum_(j=0)^k binom(k,j)
      [|A^(j)|W_(k-j)+|b^(j)|S_(k-j)+|d^(j)|T_(k-j)].

Apply this for k=0,...,4. The second component of W_5 bounds u^(6).
Every cell starts from the uniform low bounds above; this is a local
derivative estimate, not propagation of invented initial data. Two
meshes, 64 and 128 cells, enclose ALL of J. No smooth extension through
the cubic outer seam is used. Symbolic identities and independent
series identities support the implementation; they do not prove the
Born, Ward or domain statements in this ordinary argument.

## 4. Decision boundary

If this test passes, the current reference's middle-band principal
derivative modulus is closed; the enormous old conditional high-jet
allowance is unnecessary for THIS purpose. The full local moments,
variable-coefficient differentiated remainder, bolt/outer weights,
actual inverse/admission and physical compatibility remain open.
No full differentiated operator norm follows from this one direction.

Close the feasibility audit with the full inequalities and failed versus
unknown constants preserved in the planning report. Select the physical
finite-domain defect/adjoint-bolt solve with fixed-scheme local moments
and a matched ultraviolet parametrix. Do not commission another generic
source estimate or separated sample. Its first calculation must establish
the same domain for the physical boundary work, range compatibility and
the differentiated quantum action; otherwise an inverse for a surrogate
domain cannot decide attainability. The static matter pairing does not
supply retarded gravitational boundary data or all-sector stability.

Remaining to-do list: evaluate the full prescribed local moments and
differentiated defect/adjoint-bolt pairing on the physical domain; prove
admission, nonlinear closure and exact compatibility; physical matching,
causal realization, validity, all-sector stability, charge/return, opening,
transport, preparation, separate creation and the full stable-wormhole goal.
