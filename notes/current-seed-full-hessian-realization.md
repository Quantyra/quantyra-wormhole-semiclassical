# Full Hessian realization on the current finite seed

2026-10-22, E003/S169. Ordinary argument before supporting verification.
The full stable semiclassical wormhole objective remains incomplete.

The current product-centered FINITE seed has a closed coupled Hessian
realization and a unique linear response on its physical mass/fast class.
Its inverse is bounded by one in the actual positive-diagonal action-form
norm. The mass-extended weighted graph inverse exists, but its numerical
constant and the physical C5 return are still unknown. Neither statement
certifies Y0, the cold error, a nonlinear solution or Lorentzian stability.

This is a direct finite-seed argument. We do not assume a small norm
difference from the RN Hessian or substitute the product L2 inverse.

## 1. Fixed seed, exact local split and relevant norms

Keep the action, original scale, flux, paired HH domain and prepared
positive temperature of [the fixed-clock test](fixed-clock-seed-bolt-feasibility.md).
Write x=lambda B^2, C=B^2/A^2=1+6x+d, t=r/B and

    F=1+x t^2-2a/t+(1+3x+d)/t^2-2d/t^3+d/t^4,
    R4=-12lambda+w, w=d(4t^-5-6t^-6)/B^2,
    F_l=1-384piG c_R lambda, epsilon=G/(60pi),
    D_R=96piG c_R-epsilon, C_v=1+11epsilon lambda/3.

The current guards are x>8e9, abs(d)<8.6e-19, abs(d)/C<1.71e-29,
F_l>1-1.179e-8, D_R lambda<2.946e-9, G lambda<2000^-10,
abs(log(mu0 B))<1e-8, C/q<.0035 and M_n>18 sqrt(n).

Let Z=piG q^2/(g_R^2 F_l), the physical effective charge coefficient,
and V_E=2+6lambda X-2Z/X, X=r^2. The covariant seed potential is

    V0=2+6lambda X-2B^2(1+3x+d)/X
                         +8d B^3/X^(3/2)-6d B^4/X^2,
    Delta V=V0-V_E=-2Delta Q/X+8d B^3/X^(3/2)-6d B^4/X^2,
    Delta Q=B^2(1+3x+d)-Z.                            (1)

Let I_ref be the Euclidean seed-potential action with prefactor F_l/G;
it is exactly stationary at this seed. The original action splits as

    I_full=I_ref-I_pot+I_sq+Gamma_q,
    I_pot=I_ref-I_E,eff,
    I_sq=-c_R integral w(g)^2 dvol_g,
    Gamma_q=Gamma_D minus the fixed vacuum-volume functional. (2)

The identity is the retained completion of the curvature square, including
its boundary action. The reference potential is only added and subtracted.
No finite coefficient, magnetic flux, physical norm or renormalization
scale changes. The seed is not stationary for I_full.

Use a genuine TT field k=diag(-b-2xi,b,xi,xi), with radial component b:

    (r^2 F b)'=-(r^2 F'-2rF)xi.

The regular solution retains its outer tail. In physical four-volume set

    Ebar=int [F(xi'+2xi/r)^2+lambda xi^2]dvol,
    E_eff=F_l/(8piG) int [(3/2)F(xi'+2xi/r)^2+lambda xi^2]dvol.

The exact rational-family proof in
[the global TT reduction](stationary-global-tt-reduction.md) gives
I_ref''>=E_eff, ||k||^2<=22||xi||^2, ||nabla k||^2<=5200 Ebar.
Its York projection potential is >=1.49lambda. These statements already
cover abs(d)<=1, x>=1000 and the present thermal polynomial; they do not
need a translation divided by the tiny preparation gap. The same
polynomial gives lambda<=F''<=21lambda, F/r^2<=5lambda,
0<=F'/(2r)<=11lambda. Thus the disk return and its off-shell contact
bound from stationary-seed-quantum-form.md apply with the current metric:

    contact cost <=3e5 (G||E_g||infinity/lambda) E_eff/F_l. (3)

## 2. Bound the previously missing local seed contributions

Let e=B^2/B_RN^2-1, abs(e)<5.5e-23, where the CURRENT c_R RN reference
obeys Z=B_RN^2+3lambda B_RN^4. Exactly,

    Delta Q/(lambda B^4)
       =[1-1/(1+e)]/x+3[1-1/(1+e)^2]+d/x.

Hence abs(Delta Q)/(lambda B^4)<3.4e-22. In the volume-preserving
tensor exponential chart X->X exp(tau xi), I_pot has density proportional
to Delta V/X. Its three powers of X have exponents -2,-5/2,-3.
Their second derivatives therefore give the explicit bound

    abs(I_pot'')/E_eff
       <=4 abs(Delta Q)/(lambda B^4)+52 abs(d)/x<1.37e-21. (4)

Writing U=Delta V/X, direct variation of the spherical volume and X gives
G||E_pot||/lambda <=F_l(2||U||+||X U'||)/(16pi lambda)<5.5e-23.
Its additional gauge/chart cost in (3) is <1e-15 E_eff. These are
subtracted contributions in (2), estimated in absolute value.

The scalar-curvature deviation has explicit covariant guards

    w0=||w||/lambda<1.1e-27,
    w1=||nabla w||/lambda^(3/2)<1.4e-26,
    w2=||Hess w||/lambda^2<4.3e-25.                  (5)

For (5), the absolute radial coefficient sums are 10,56,372 for
w, r w', r^2 w''. Use F/r^2<=5lambda and F'/(2r)<=11lambda.
The radial/time/sphere Hessian entries are Fw''+F'w'/2, F'w'/2,
Fw'/r. Their Frobenius bound is <4000 abs(d)/x. The gradient bound
uses sqrt(5)*56 abs(d)/x. No singular bolt connection is bounded alone.

The current traceless Ricci has eigenvalues (-tau,-tau,tau,tau), where

    tau=(1+3x+d)/(B^2 t^4)-d(5t^-5-(9/2)t^-6)/B^2,
    abs(tau)<3.01lambda, delta R4[k]=-4tau xi.

In the tensor exponential chart the volume is fixed, so the exact square
Hessian is

    I_sq''=-2c_R int [(delta R4)^2+w R4'']dvol.       (6)

The second term must be retained: it vanished on RN but not here.
Here is a bound that does not assume fourth derivatives of k. Let
K=C'(0), L=C''(0) for the connection difference of g exp(tau k) and g.
Then ||K||<=3||nabla k||/2, ||L||<=6|k||nabla k| and C^c_ac=0.
The exact scalar second variation is

    R4''=(k^2):Ric-2k^ab nabla_c K^c_ab
                         +g^ab nabla_c L^c_ab-2g^ab K^c_bd K^d_ac.

After one integration by parts against w, Cauchy--Schwarz gives

    abs(int w R4'')
      <=int [12|w||nabla k|^2+15|nabla w||k||nabla k|
                                         +|w||Ric||k|^2].          (7)

All constants follow from the two connection norms and |g|=2; using
|Ric|<100lambda is conservative. Inserting the TT bounds into (7)
shows that the wR4'' term in (6) costs <1e-30 E_eff. The first term
costs at most

    256pi c_R G lambda (3.01)^2/F_l<7.2e-8 E_eff.     (8)

The square gradient is also nonzero. Its explicit Euler tensor obeys
G||E_sq||/lambda <=c_R G lambda(200w0+w0^2+10w2).
Equations (3),(5) put its gauge/chart contribution below 1e-27 E_eff.
Thus neither zero background curvature deviation nor full stationarity
has been smuggled into the current-seed Hessian.

## 3. The actual quantum TT form on this seed

Apply the ordinary curved-form argument, not only a source value bound,
of stationary-curved-quantum-tt.md and recentered-quantum-tt-form.md.
The latter's d=0 specialization is replaced by the original rational
family guards in stationary-centered-exterior-source.md. In particular:

* The canonical polynomial has derivatives bounded by (5,2,17,13)C
  through order four, and fifth/sixth derivatives <=960|d|,720|d|.
  Normal curvature jets retain (2,17,56,115). The same cold complex cap,
  non-bolt separation, warp/taper bounds and reflected cubic apply.
* The actual paired first-order inverse is <=1/M_n. Both partner domains,
  heat contacts and the LLL HH allocation remain. A cap is not substituted
  for the actual operator. The fixed-scale common restoration has
  abs(Fang)<100 and 13.9<kappa4<17.5, including log(mu0 B).
* The normal-cap, two-buffer gluing, tadpole plus bubble, constant and
  second contacts have exactly the hypotheses of the retained proof.
  The geodesic comparison's positive sixth remainder, common fourth
  contact and the full outer logarithmic compensation are kept.
* The TT return constants use the geometric inequalities of section 1,
  which hold for this d family. The same physical conversions are
  C J1*/E0<1e-26, C J0*/E0<2e-28 and
  Gq^2/(lambda B^4)<4e-27, with E0=E_eff/F_l.
  The outer cutoff measure gain is retained; it is not replaced by a
  global L2 bound. The vacuum and actual physical stress gradient are
  bounded as in the seed source proof. I_pot and I_sq are separate rows,
  not double-counted in the quantum contact budget.

The resulting negative quantum budget in E0 is

    .001 + .002 + .0006 + 1e-12 + 1e-12 + 1e-25 + 1e-6
       <.003601000003.

The last term includes all quantum off-shell gauge/chart contacts; no
gradient is declared zero. Division by the fixed F_l gives <.005.
Together with (4),(8) and their contacts, the full returned TT form l
therefore satisfies

    l>.99 E_eff.                                    (9)

This is an inequality for the ACTUAL combined Hessian on a smooth
bolt-regular TT core and its retained stronger completion. It is not
yet an inference of closedness in its own form norm.

For the mixed estimate keep the charge potential instead of discarding
it. The classical potential in the global TT identity is

    V=3lambda-12m/r^3+13B^2(1+3x)/r^4-7/(3r^2)
                   +(d/B^2)(13t^-4-51t^-5+45t^-6).

The exact degree-eight Bernstein certificate in recentered-mixed-form.md
has minimum .09. Its remaining loss is <33/x+2/x^2. Replacing its charge
by tau above and adding the d potential costs <3e-26 lambda. Hence
V>=1.99 tau^2/lambda+.08lambda still holds. Subtracting the quantum,
potential and curvature-square budgets proves l>=X^2, where

    X^2=F_l/(8piG) int [ .99(3/2)F(xi'+2xi/r)^2
                          +1.98 tau^2 xi^2/lambda+.07lambda xi^2]. (10)

## 4. Complete scalar and mixed rows on the current metric

The exact nonconstant-curvature conformal density row is already derived
in product-centered-scalar-cold-response.md:

    A s=6D_R Delta^2 s-div(K grad s)+(2R4+48lambda0)s,
    K=[6-(2D_R+11epsilon/6)R4]g+(22epsilon/3)J.

Using (5), |J|<7lambda and the positive D_R gives

    a(s,s)=(8piG)^-1 int s A s
      >=S^2=(8piG)^-1 int [6D_R(Delta s)^2
                               +5.99|grad s|^2+23.99lambda s^2]. (11)

The ACTION scalar block is -a. The canonical determinant has no
conformal or mixed derivative in quotient coordinates; the local anomaly
and off-shell density/return contacts remain in the mixed block.
For a genuine TT field its exact trace variation is

    (C_v-11epsilon w/36-D_R Delta)delta R4[k]
     -D_R delta(Delta)[k]w
     +epsilon[(7/12)delta C4^2+(11/6)delta J^2-20delta I]
     -div(T0 eta),                                   (12)

where delta R4=-4tau xi, delta(Delta)[k]w=-k:Hess w and eta is the
outer-fixed TT disk return. The last identity uses div k=tr k=0;
the moving-density contact is not omitted. The complete current T0
obeys ||T0||/lambda<9e-28, not zero.

Integrate the main Laplacian onto s. From (10)-(11), its relative square
constant is at most

    Kmain^2 <=[16/(1.98 F_l)] [C_v^2/23.99+D_R lambda/6]<.337.

The current rational curvature retains |Rm|<60lambda, |Ric|<13lambda,
|J|<7lambda, |C4|<35lambda, |grad J|<100lambda^(3/2) and
|grad C4|<1e4lambda^(3/2). The d contributions are bounded by (5)
and their explicit projector derivatives, far inside those margins.
The retained integrated curvature-variation identities therefore give
the same <3e-28 curvature contribution. The magnetic term is <3e-29;
using tau in place of Z/r^4 changes its bound by <2e-22 relatively.
The return term is <9e-28 sqrt(258/(5.99 F_l))<6e-27.
The new -D_R delta(Delta)w term is bounded by
D_R lambda w2 sqrt(315/(23.99 F_l))<5e-33. The w delta R4 term is
smaller. Thus the remaining constants sum to <1e-25, giving

    |b(s,xi)|^2<.338 S^2 X^2<=.338 a(s,s)l(xi,xi).    (13)

All quantities retain the physical 8piG and volume factors. This bound
does not improve by renaming a norm. The old .336 RN value is not claimed
on the new seed; (13) is sufficient and still below 3/8.

## 5. Actual graph realization, not a substituted comparison form

Use z proportional to 1/r, critical physical weight 3/2, and

    Z_delta=z^delta H_b^4, W_delta=z^delta H_b^(4,log),
    Y_delta=z^delta L_b^2, dmu=gamma(z)z^-3 dz/z.

The scalar fourth coefficient is nonzero. The shape has the full
logarithmic fourth principal part. Mixed terms have at most two
derivatives after the exact quotient/Ward contacts; the extra local
w terms above do not change either principal order. At the regular
bolt use rotational disk charts, with no puncture boundary condition.

The current canonical coefficients and boundary operators are analytic
one-sided, with the same leading constant/quadratic jets, mass cubic
and fixed AdS scale. The actual-minus-normal finite-part kernel has the
same two-power gain as recentered-reference-inverse.md: comparable
annuli cost O(h^2 log^N h) from the stated graph norms to L2; separated
remainders have coordinate degree -3 in place of -5. The full common
contacts and state-dependent separated terms are retained. Their
constants are finite at the fixed positive temperature, not numerically
estimated. Full-order errors on a fixed annulus are not called compact.

This follows from the same finite direct/reflected expansion, now using
the polynomial and domain guards in section 3. The normalized-symbol
Sobolev bound uses [Grubb, Theorem 2.3](https://web.math.ku.dk/~grubb/CPDE14.pdf)
at output order zero with positive spatial Holder regularity and two
frequency derivatives in dimension one. The smooth differential-boundary
resolvent estimates and singular logarithmic term of
[Gaarde--Grubb, Theorem 2.4 and Proposition 2.9](https://arxiv.org/pdf/math/0703878)
are retained. These inputs do not alone constitute a determinant-Hessian
theorem; the finite kernel/contact construction above is essential.

The outer normal families are IDENTICAL to those in the recentered
reference inverse: J_R(s) and N_T(s) depend on lambda, F_l, D_R, C_v,
mu0, and q, all fixed. The new decaying curvature/charge terms have zero
normal family. On -1/2<=Re s<=7/2 the scalar has no roots; the only TT
roots are the simple 0 and 3. The removable points 1,2 are not roots.
The log symbol is elliptic at large frequency; at the critical line
there is no normal zero. The reflected York congruence is unchanged.

The actual compact-interior parametrix and normal inverse, followed by
Neumann inversion on a sufficiently short end, give left/right graph
parametrices and maximal-adjoint regularity. The remainder cutoff
commutators are compact; no numerical cutoff is claimed. Symmetry and
this regularity identify the TT self-adjoint domain W_(3/2), and the
positive scalar domain Z_(3/2). A compact sphere component's TT radial
tail can be removed by a distant bump of amplitude O(L^-2); its critical
graph norm is O(L^-1/2). The leading TT moment kernel is still -2r+O(1),
so the existing genuine-TT core argument applies. Smooth bolt values
remain in the core; a punctured-bolt fourth-order approximation is not
used. Consequently the actual positive forms are closed in their own
spectral form norms, and (9)-(13) extend to them.

Let V be their direct sum with ||(s,xi)||_V^2=a(s,s)+l(xi,xi).
The actual reduced, density-correct Hessian is H=[-A,B;B*,D_T].
For J_sign(s,xi)=(-s,xi), the mixed terms cancel in the real pairing:

    Re <J_sign u,H u>=||u||_V^2,
    ||H^-1 f||_V<=||f||_V*.                           (14)

Boundedness follows from (13), so Lax--Milgram supplies the inverse.
For L2 data graph regularity identifies it with the actual maximal
operator solution. Equation (14) is an ACTION-FORM inverse estimate,
not an admission-norm estimate or a Lorentzian stability claim.

## 6. Physical tail, constraints and what is still quantitatively open

The same finite-strip improvement from weight 3/2 to 7/2 crosses only
the TT root 3. It gives the bounded bijection

    H: Z_(7/2) x [W_(7/2)+span{cutoff z^3}]->Y_(7/2)^2. (15)

Its numerical inverse norm is uncomputed. For the actual seed residual,
the current fixed-seed source expansion has weight-four finite log
polynomials. The endpoint argument of recentered-physical-endpoint.md
applies: use the original equation and its actual remote moment, then
apply the row combination (D_z-4)f2-(D_z-6)f1/3 for regularity only.
The normal multiplier (s-4)N_T(s) has the same nonzero endpoint bound.
The resulting solution has a correlated mass cubic plus six conormal
disk jets of O(z^4 log^p z), hence five physical areal jets. No original
weight-four equation or global moment is discarded. The full linear
response therefore lies in the selected physical mass/fast phase.

The [current boundary check](current-seed-boundary-phase-and-trace-defect.md)
and recentered-mass-boundary.md identify zero fixed-source boundary work
on that phase, including its linearization. Thus (15) acquires no
additional mass or scalar boundary equation. This uses the retained C1
action/cut limit on the declared phase, not arbitrary finite-wall data.
The current/Ward identity holds on this finite seed. Off shell, a Newton
step still has the connection-times-seed-residual remainder; it is not
an exact solution of all nonlinear constraints. No old projection is
imported and no constraint is discarded.

For the previously constructed reference correction h1 and the actual
full linear response h_lin, the useful identity is now well defined:

    r1=F_eq(seed)+H h1,
    h_lin-h1=-H^-1 r1, ||h_lin-h1||_V<=||r1||_V*.      (16)

The leading outer trace coefficient of r1 was evaluated previously;
its full norm, shape part, cold return and higher weighted derivatives
are not bounded here. A signed bolt return can equivalently be evaluated
by the adjoint finite-domain response to the physical bolt functional.
Its norm and the needed full pairing remain to be computed. Do not use
the unit number in (14) as K_c or as a bound on the physical C5 return.

Research decision: retain this construction. The selected full finite-seed
linear problem now has an actual realization and physical-phase solution;
its quantitative correction is the remaining decisive test. Next evaluate
the signed full-defect/adjoint-bolt pairing and weighted remainder, keeping
the mass/clock/scalar phase fixed. Require cold difference <1.667e-9 and
total Y0<=5e-9. The rho=1e-8 ball still needs Z0<=.25, Z2<=5e7 and all
admission margins. Existence of a finite graph constant does not meet
those thresholds. The feasibility audit remains closed; the generic
unsigned source chain stays parked.

Remaining to-do list: numerical full-defect return and weighted error;
preserved nonlinear solution, physical matching, semiclassical validity,
all-sector stability and the remaining stable-wormhole requirements.
