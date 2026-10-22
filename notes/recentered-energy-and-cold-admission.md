# Actual linear energy bound and the cold admission obstruction

2026-10-22, E003/S169. Ordinary argument before supporting verification.
The full stable semiclassical wormhole objective remains incomplete.

The ACTUAL coupled reference response to its complete Dirac residual
has normalized energy less than 1.44e-10. This uses the actual reference
inverse, not the normal Volterra comparison or an Einstein-only response.
Nevertheless a generic conversion of this energy allowance to the desired
physical radius fails: any bound controlling even the bolt value has
constant greater than 2.6e15. There is no finite energy-to-C5 embedding.
These are failures of a proposed certificate, not a lower bound on the
actual correction or a nonexistence result. The next estimate must use
the actual cold source and its response, not an arbitrary energy field.

## 1. Fixed reference, normalization and available requirements

Keep the same c_R=10^22 recentered RN reference, fixed q=k^4, k=2000,
g_R=k^-4, mu0, original positive temperature and actual paired Dirac
domain. No parameter, subtraction or reference radius is changed here.
Put lambda B^2=x, C=1+6x, r_+=B(1+h),

    f=1-2M/r+Z/r^2+lambda r^2,
    F_lambda=1-384pi c_R G lambda,
    D_R=96piG c_R-G/(60pi), epsilon=G/(60pi),
    C_v=1+11epsilon lambda/3.

The [new-reference guards](recentered-quantum-tt-form.md) give

    x>8e9, C/x<6.1, 0<h<3*2000^-10<1e-30,
    log(1/h)<80, F_lambda>1-1.179e-8,
    Gq^2/(lambda B^4)<4e-27, D_R lambda<3e-9,
    C/q<.0035, q/C<400, M_n>18 sqrt(n),
    abs(log(mu0 B))<1e-8.                               (1)

Remove the common factor 4pi beta_T from radial forms and pairings,
without changing their relative normalization. Let a_rad,l_rad be the
actual positive scalar-trace and TT forms, and define

    E(u)^2=[8piG/(lambda B^3)]
                    [a_rad(s,s)+l_rad(xi,xi)].          (2)

This displays the physical energy unit; it does not alter a norm product.
The actual saddle inverse has norm <=1 from V* to V in these forms.
The actual response H u=-F0 is already known to have the physical fast
tail and regular bolt. Its size in the complete admission norm is open.

## 2. Full source bound on this reference

The local EH+c_R R^2+Maxwell tensor of this RN reference equals the
constant outer Dirac vacuum tensor. Thus its complete action gradient
is represented by the full Dirac stress minus that vacuum. This is not
the old product-centered Gamma_rest split. No product stress is silently
subtracted, and no first-order constraint solution is used as H^-1 F0.

The retained full canonical source proof applies at this d=0 reference
after checking (1) and its analytic polynomial geometry. In particular
its trace derivative and anchor are

    abs(Theta'(t))<=K q^2, K=1e15,
    abs(Theta(2))<=K0 q^2, K0=1e11, t=r/B,
    abs(B^2 rho2,mass),abs(B^2 p2,mass)
                              <=(K+K0/2)q^2.           (3)

Here rho2 and p2 are sphere-integrated stresses, not four-dimensional
components. The transfer retains the actual massive resolvents, cold
and non-bolt caps, angular sums and fixed-scale finite parts. The new
scale bound is 1e-8, not the old 1e-21. In the anchor formula

    Fang_eff=4gamma_E-2H_q+4log(mu0 B),
    a0 bracket=zeta_h'(-1)+[1+2log(mu0 B)]zeta_h(-1),
    zeta_h(-1)=1/30-q^2/6,

this change is far inside the retained <100 Fang guard and the strict
<5e10 q^2 anchor subtotal before its 1e11 allowance. The source derivative
has the retained <3e14 q^2 subtotal before its 1e15 allowance. All geometry
and domain guards were checked on this new reference; the bare old source
number is not being assigned to an arbitrary corrected geometry.

The [exact Weyl RN formulas](stationary-weyl-radial-sources.md) give after
vacuum removal

    B^2 rho2,W,ren=L/t+R_W,
    abs(L)<q^2, abs(R_W)<=2q^2(1+log t)/t^2.

Since t>=1 and (1+log t)/t^2<=1, its value is <=3q^2. The companion
pressure obeys p2,W,ren=-rho2,W,ren-f E_F/(2pi), where

    B^2 E_F=2[6Q_c log t-8Q_c+3a_m t]/(15t^4),
    a_m=M/B<q, Q_c=Z/B^2<q, f<5x t^2, x<q.

Using log(t)/t^2<=1/(2e)<1/2 bounds abs(B^2 f E_F/(2pi)) by
5q^2; hence abs(B^2 p2,W,ren)<=8q^2. The actual LLL trace and
derivative bounds give both component values below q^2, using C/q<.0035.
Thus a convenient complete allowance is

    abs(rho),abs(p)<=C_s U(r),
    C_s=1.01e15, U(r)=q^2/(4pi B^2 r^2),               (4)

where rho,p now denote FOUR-dimensional vacuum-subtracted components.
Only spatial source estimates at the stated reference are used.
The full [common-subtraction source proof](stationary-centered-exterior-source.md)
is the comparator; the mixed-boundary heat coefficients are those of
[Vassilevich, section 5.3](https://arxiv.org/html/hep-th/0306138v3#S5.SS3).
A heat coefficient alone would not prove (3) or the state-dependent stress.

The exact local trace gives

    tau=-(28M^2/r^6-56MZ/r^7+(106/3)Z^2/r^8-10q^2/r^4)/(480pi^2).

Signs do not matter for the estimate. Since a_m,Q_c<q, its absolute
coefficient sum is 388/3. Therefore abs(tau)<=U(r), because
388/(360pi)<1 and t>=1. The angular pressure satisfies
p_o=(tau+rho-p)/2, so abs(p_o)<=(C_s+1/2)U. Consequently

    ||T_ren||_L2,rad <= C_T q^2/(4pi B^(5/2)),
    C_T=2.021e15,
    2C_s^2+2(C_s+1/2)^2<C_T^2.                         (5)

We used integral_(r_+)^infinity dr/r^2=1/r_+<=1/B.
The vacuum and the entire matter angular tower remain counted once.

## 3. A cold-uniform bound for the genuine TT lift

For the genuine tensor k=diag(-b-2xi,b,xi,xi), regularity fixes

    b(r)=-W(r)^-1 integral_(r_+)^r A(t)xi(t)dt,
    W=r^2 f, A=r^2 f'-2rf.                             (6)

There is no free singular bolt solution. f'>0 on the exterior: the
numerator lambda r^4+Mr-Z is increasing and positive at r_+.
Thus W'>0, |A|<=W', and W runs from zero to infinity. Change variable
to w=W(r). The radial L2 measure is

    r^2 dr=omega(w)dw,
    omega=r^2/W'=1/(4lambda r+2/r-2M/r^2).

The denominator derivative is 4lambda-2/r^2+4M/r^3>0 by (1), hence
omega is nonincreasing. For the Hardy average P v(w)=integral_0^1 v(tw)dt,
Minkowski and omega(y/t)<=omega(y) give

    ||P v||_L2(omega dw)
      <=integral_0^1 t^-1/2 dt ||v||=2||v||.

Applying this to (6) proves, with NO inverse-temperature factor,

    ||b||_L2,rad<=2||xi||_L2,rad,
    ||k||_L2,rad^2<=22||xi||_L2,rad^2.                 (7)

The returned TT tangent differs by a radial diffeomorphism. At the
reference, the COMPLETE residual is conserved by the verified current/
Ward identity. Its first pairing with that diffeomorphism is a boundary
term, zero at the regular bolt and fixed-source fast end. Thus the
source pairing can be evaluated on k itself. This does not remove the
return/density contacts from the Hessian; those remain in a,l,H.

## 4. Bound the actual coupled energy correction

The TT source functional is the half stress pairing, so (5)-(7) and
the actual lower form l_rad>=(.07 F_lambda lambda/(8piG))||xi||^2 give

    ||F_TT||_(l*) <= (sqrt(22)/2)
            sqrt(8piG/(.07F_lambda lambda)) ||T_ren||_L2,rad.

In the normalization (2), its source allowance is therefore

    E_TT,source <= sqrt(22) C_T/sqrt(.07F_lambda)
                                  *Gq^2/(lambda B^4). (8)

The [exact scalar source integral](recentered-r2-scalar-boundary.md)
gives E_scalar,source<1.82e-29 in the SAME normalization. The dual of
the direct-sum form norm combines these by their squares. The actual
coupled saddle inverse bound 1, with every mixed term retained, yields

    E(H^-1 F0)
      <=sqrt(E_scalar,source^2+E_TT,source^2)<1.44e-10. (9)

This is an actual-response ENERGY upper bound. It is not a lower bound,
a physical C5 estimate, a mass-coefficient bound or a nonlinear radius.
The exact constants in (1),(5),(8) are substituted in the checker; no
normal model is used to stand in for H. The scalar/TT boundary and
current hypotheses have been proved in the preceding milestones.

## 5. Execute the cold admission test with an actual reference field

For d=B h put s(r)=exp[-(r-r_+)/d] and xi=0. This is a smooth radial
disk field at the bolt, where r-r_+ is a smooth even function of proper
distance. It has s(r_+)=1, preserves the thermal cone as a conformal
tangent and decays faster than the physical fast tail at infinity.
Smooth cutoffs approximate it in the scalar form. It is a test tangent,
not the actual inverse correction or a finite-amplitude solution.

Put w=(r-r_+)/d. The exact scalar form is

    a_rad=(8piG)^-1 integral r^2 {
       6D_R(Delta s)^2+K_r f(s')^2+24lambda C_v s^2}dr,
    K_r=6+24D_R lambda+22epsilon lambda-(22epsilon/3)Z/r^4
                                                         <6.000001.

At the actual reference f''_t<=21x, f'_t(1+h)=2C h K(h)<=2C h;
the temperature factor K(h)<=1 follows directly from its decreasing
positive-power formula. Taylor's integral bound then gives

    f/(B^2 h^2 lambda)<=P(w)=12.2w+10.5w^2,
    abs(Delta s)/(lambda exp(-w))
                           <=Q(w)=12.3+33.3w+10.6w^2.             (10)

For the second inequality use
Delta s=exp(-w)[f/d^2-(f'+2f/r)/d] and h/(r/B)<=1e-30.
Also (r/B)^2=[1+h(1+w)]^2<=1+3h(1+w)^2 and 24C_v<24.000001.
All moments are elementary: integral_0^infinity w^j exp(-2w)dw=j!/2^(j+1).
Inserting (10) into the ACTUAL form, not just its lower comparison, gives

    E(s,0)^2 <= h integral_0^infinity exp(-2w)
       [1+3e-30(1+w)^2]
       [18e-9 Q(w)^2+6.000001 P(w)+24.000001]dw
                                                        <47h.    (11)

This is a bound for the finite reference, not only a near-horizon limit.
Direct numerical evaluation of the exact rational RN form below checks
it independently; the proof is the displayed rational moment inequality.
The limiting coefficient is approximately 39, consistently below 47.

Any linear embedding estimate controlling even the physical bolt
radius change, |s(r_+)|<=K_c E(s,xi), therefore necessarily obeys

    K_c>1/sqrt(47h)>2.6e15.                            (12)

The second number uses only h<3*2000^-10. In contrast, to turn the
available bound (9) into the design Y0<=5e-9 by this generic route would
require K_c<=5e-9/(1.44e-10)=625/18, about 34.72. The resulting upper
certificate cannot meet that target. This combines a known norm LOWER
bound with a known response UPPER allowance; it does not lower-bound
the response's actual bolt value. Its signed source could couple very
weakly to these test fields. Computing that coupling is the next task.

For the full physical C5 norm the generic energy embedding is worse:
on a fixed compact annulus scalar oscillations have a(s,s)^1/2 of order
frequency^2 and fifth metric derivatives of order frequency^5. Thus
their ratio is unbounded even at fixed temperature. This is not merely
an unknown finite embedding constant. A graph/source-specific estimate
is indispensable. The previous mass-trace discontinuity likewise means
that (9) alone cannot bound the asymptotic mass coefficient.

## 6. Research decision

Keep this construction and the new actual energy bound. Do not attempt
to certify a physical neighborhood by multiplying it by a generic
cold/derivative embedding. The decisive test (11)-(12) has been executed
and rejects that proposed certificate without rejecting the solution.

Next isolate the actual product-centered cold source and bound its
response in the regular disk, retaining its signed/correlated terms and
the exact fixed-period condition. The full homogeneous c_R product
exists, but its previous 9.43e-8 allowance relative to the old classical
product is not the actual mismatch from this recentered RN reference.
Compute that mismatch and use it in a specified reference change if
needed; do not divide an uncentered radius shift by the tiny h or assume
that a translated bolt removes every cold coefficient. The outer
combined-tail method and the actual global mass/kernel moments remain
part of the same admission problem. Neither a normalization change nor
another bound on arbitrary high derivatives solves this obstacle.

The nonlinear targets rho=1e-8, Y0<=5e-9, Z0<=.25, Z2<=5e7 and every
positive admission margin remain uncertified. Physical compatibility,
same-source matching, omission/noise control and all-sector retarded
stability remain independent gates. S169.2 is still open. Freeze this
argument, checker and source review before running; preserve all runs
and source snapshots. The focused check should include exact Hardy
geometry and moment arithmetic, direct stable RN quadrature at 100/150
digits for the cold field, and independent integral TT-lift fixtures.
These diagnostics do not replace the quantum heat estimates or prove
a physical nonlinear ball. No Lean invocation is selected.

Remaining to-do list: actual product-centered cold response and complete
admission norm, including outer cancellation and mass/kernel moments;
finite-state transfer, nonlinear existence, compatibility, matching,
semiclassical validity, all-sector stability, charge/return, opening,
transport, preparation, separate creation and full-result audit.
