# The quantum TT lower form on the recentered reference

2026-10-22, E003/S169. Ordinary argument before frozen verification.
Previous goal turn: progress, technical c9ea661 and planning 0259bf2.
The full stable semiclassical wormhole objective remains incomplete.

At the SAME c_R=10^22 recentered RN reference, the complete quantum
remainder and off-shell chart/gauge contacts satisfy

    Gamma_rest,TT'' >= -.005 E_eff.                       (1)

The tangent class is the smooth bolt-regular TT core and its explicit
stronger-form completion described below. Thus the preceding <=1/4
negative-form target is met on this reference and class. Combining
the [local square budget](recentered-canonical-tt-audit.md) gives

    I_full,TT'' >= (1-.005-7.07e-8) E_eff > .994 E_eff.    (2)

This is a lower bound for a specified TT pullback of the actual Hessian,
not the full coupled physical boundary realization, an inverse in an
admission norm, finite-geometry transfer, or stability of a wormhole.
We do not deduce closedness of the actual Hessian in its own form norm
merely from continuity in the stronger norm. That issue and the physical
coupled domain remain separate. The new result is obtained by checking
the hypotheses of the retained continuum estimates at this reference,
not by assigning the old seed's .005 allowance to a changed metric.

## 1. Geometry and physical scale, with numerical transfer guards

Keep k=2000, q=k^4, g_R=k^-4, eta=1/k, c_R=10^22,
mu0=1/b_class and the ORIGINAL prepared kappa=z_i nu0,
z_i in [1/4,7/25]. Define

    g0=k^-10, g+11g^2/(360pi)=g0, lambda=g/G,
    F_lambda=1-384pi c_R g,
    C0=sqrt(1+12pi k^6), x0=(C0-1)/6,
    b_class^2=x0/lambda0, lambda0=g0/G,
    nu0=2q eta G C0/b_class^3,
    C=sqrt(1+12pi g k^16/F_lambda),
    x=(C-1)/6, B^2=x/lambda, ell=B/sqrt(C).             (3)

The exact reference is f=1-2M/r+Z/r^2+lambda r^2,
Z=B^2(1+3x). Its horizon is r_+=B(1+h),
M=B(1+2x+delta/2), delta=h^2[1+x((1+h)^2+2(1+h)+3)]/(1+h).
The canonical sphere B is FIXED during variations. Re-expressing the
determinant at this sphere retains the fixed mu0 in every finite part.

Using (1-1e-34)g0<g<g0, 3.14159<pi<3.14160 and (3) gives

| Required input | Bound on this new reference |
| --- | --- |
| F_lambda | >1-1.179e-8 |
| C | 49119000000<C<49200000000 |
| x | >8e9 |
| C/q; q/C | <.0035; <400 |
| M_n^2=n(n+q)/C | M_n>18 sqrt(n) |
| C/x | <6.1 |
| GC/B^2 | <8.5e-33 |
| Gq^2/(lambda B^4) | <4e-27 |
| abs(log(mu0 B)) | <1e-8, NOT <1e-21 |
| -log(mu0/sqrt(lambda)) | >9 |

For the scale guard, C/C0 lies between 1-1e-33 and 1+1e-8.
Then B^2/b_class^2=(C-1)g0/[(C0-1)g] differs from one by <1.01e-8.
The logarithm bound follows from its positive interval and the elementary
log inequalities. No renormalization scale has been reset.

The new gamma=kappa B/C has ratio in (.99,1.01) to the classical
gamma0=kappa b_class/C0. The exact d=0 temperature equation is

    gamma=h K(h),
    K(h)=[(1+3x)((1+h)^-2+(1+h)^-3)
                         +3x(1+(1+h)^-1)]/(2C),
    K(0)=1, abs(K')<3 for 0<=h<=.001.                 (4)

The original positive-temperature bound eta/(8k^9)<=h_class<=eta/k^9
therefore gives, without a zero-temperature limit,

    eta/(24k^9)<h<3eta/k^9<1e-30, log(1/h)<80,
    delta/C<1.001h^2.                                (5)

For the logarithm, e>2.7 and 2.7^80>24k^10 suffice. This checks the
actual small-h hypothesis used by the old cold cap; the looser previous
h<6e-30 bound alone would not have checked that hypothesis.

## 2. Direct canonical geometry and the actual matter realization

Put u=B/r, a0=1+2x+delta/2, Q0=1+3x. In normalized proper units,

    calB(u)=x+u^2-2a0u^3+Q0u^4,
    du/dY=sqrt(calB/C), Rbar=-calB''/C.                (6)

This is precisely the d=0 subfamily of the explicit polynomial guards
in [the centered source proof](stationary-centered-exterior-source.md),
not a small perturbation divided by the tiny horizon displacement.
Here a0/C<.334, Q0/C<.501, calB<=x+1<=C/5. On the physical interval,
the convex quadratic calB'' has upper endpoint value <=2C and minimum
2-3a0^2/Q0>-2C. Coefficient sums give

    (abs(calB'),abs(calB''),abs(calB'''),abs(calB''''))/C
        <=(5,2,17,13), calB^(5)=calB^(6)=0,
    (abs(Rbar),abs(Rbar_Y),abs(Rbar_YY),abs(Rbar_YYY))
        <=(2,17,56,115).                              (7)

For example Rbar_YY=-(calB calB''''+calB' calB'''/2)/C^2;
the third derivative adds no fifth-polynomial term here. The same
extremal factorization and (5), with z=1-u and zc=e^-16/1001, give
Cz^2/8<=calB<=Cz^2 on z>=zc/2 and the non-bolt proper separation
>=16ell. The non-bolt logarithmic warp guards are 1.01,2.021,13;
the tapered relative jets are 2,8,100. The even cap and smooth cap
remain distinct, with ratios <=1e10 and <=4e10 respectively.

The cold polynomial has a=4x/C<=2/3 and c=x/C<=1/6 in equations
(3)-(12) of [the cold proof](stationary-reference-cold-bolt.md).
All its RN factor identities hold exactly with (3)-(5). In particular
the small H2 cap has abs(sigma_cap)<.05, first/second derivative
bounds .4/1.5, and the same normal-chart complex disk. No degree-six
perturbation or cold inverse-gap translation is needed. The normal
caps used by [the curved form proof](stationary-curved-quantum-tt.md)
therefore obey abs(v)<=8rho^2, sup abs(v)<=.026,
abs(v-T2v)<=1e4 rho^3, with the same usable radius 1/128.

The [paired phase/domain argument](stationary-paired-dirac-phase.md)
is algebraic in the identical-copy frame and monopole connection; it
does not require the old Einstein field equation. On this reference
the boundary is still geodesic in the canonical metric, the projector
and thermal spin structure are unchanged, and the bolt is smooth by
(4). Its global intertwiner gives Dsharp^dagger=-Dsharp,
-Dsharp^2=D2^dagger D2 and

    ||(Dsharp+M_n)^-1||<=1/M_n.                       (8)

Thus the MASSIVE actual-domain bound follows directly, not from
continuity of an old spectral gap. It does not assert a massless
interacting-state gap or construct the coupled physical state. The LLL
uses its retained regular HH Polyakov allocation on this new disk.
The even image cancellation and physical-half normalization remain
exact on the even static conformal paths used for the disk TT Hessian.
Both adjoint partner domains and the full angular tower are retained.

## 3. What transfers in the quantum proof, and why

The result used here is the chain of ordinary estimates in sections
2-10 of stationary-curved-quantum-tt.md, with its inputs replaced by
(3)-(8). The argument needs more than a source value bound. Its
separated kernel, constant contact and second contact are separately
fixed in the common four-dimensional subtraction. The following ledger
states the hypotheses and allocations that are carried through.

| Piece | Applicable ordinary estimate on (6) |
| --- | --- |
| Angular count and summability | N(E)<=4CE^2; sum d_n/M_n^4<4C^2/q; d_n=2(q+2n) |
| Normal cap | covariance factors, quadratic tag and convergent Born series with ratio <.162; local kernel difference <=1e8 C/distance^4 |
| Actual gluing | two buffered endpoints, rows Rloc<=7.8e7 exp(-M/270), Rfar<=9.8304e8 exp(-M/530), inverse (8) |
| Gluing quadratic form | row-sum-subtracted difference costs <=6e24 C J0, not the raw 3e24 bubble counted alongside its contact |
| Constant probe | full tadpole plus bubble 2M^4 Tr(P^-2)-2M^2 Tr(P^-1), with P=M^2-Dsharp^2 |
| Heat/contact weights | abs(B(t))<=120C/t^2 and <=800C t^-2 exp(-162t); smooth/reflected remainder reserves 1e15/1e16 |
| Actual p0 contact | abs(Theta_p0)<=1e20 q^2(1+1/Y), seam delta retained until paired with T(0)=0 |
| Second contact | compact quadratic fixes the finite second moment; local Born, full geodesic tail and anchor change fit 1e13 C |
| Radial gain-two error | <=1e9 C/abs(Y-Y')^3; with IMS and scalar comparison use K=1e13 |

The geometric identity A-B=-B[D,chi]A[D,chi]B applies only where
both endpoint cutoffs equal one. The point bound and L2 row bound are
used separately before angular summation; squaring the large point
bound over the whole manifold would fail. Reverse rows use the same
paired realization. The constant probe includes the exact actual row
sum, so no contact is silently dropped. Local quadratic probes fix the
remaining second contact; a separated-kernel envelope alone would not.

The changed scale enters the COMMON restorations explicitly:

    F_eff=4gamma_E-2H_q+4log(mu0B),
    kappa4=-F_eff/4, 13.9<kappa4<17.5,
    zeta_h(-1)=1/30-q^2/6, abs(zeta_h'(-1))<1000q^2,
    A0=2[zeta_h'(-1)+(1+2log(mu0B))zeta_h(-1)]/(4pi B^2),
    A2=B^2 F_eff k2/(4pi),
    constant-contact restoration=-2A0+2A2.            (9)

The old abs(log(mu0B))<1e-21 is replaced by <1e-8. The anchor's
strict subtotal <5e10 q^2 remains below 1e11 q^2: the new a0 bracket
changes by <(1e-8/3)q^2 and abs(F_eff)<100 still holds. The cold and
non-bolt local common-polynomial reserves also use abs(F_eff)<100,
not the 1e-21 value separately. All other restoration changes are
already present in (9). There is no modewise renormalization or extra
finite fourth contact. The same arguments consequently give

    abs(Theta'(t))<=1e15 q^2,
    abs(Theta(2))<=1e11 q^2.                           (10)

These are source values/first derivatives on this reference, not
parameter-Hessian estimates. The TT Hessian instead uses the kernel
and contact proof above. Signed differences of divergent bare terms
are never estimated separately. The common finite parts are taken
before removing the angular regulator; the positive sixth remainder
is absolutely summable, while the unique fourth divergence is (9).

## 4. The unbounded outer logarithm is controlled, not truncated

On the complete even double the geodesic flow preserves Liouville
measure. Its C1,1 vector field (metric C2,1 at the cubic seam) defines
a strongly continuous unitary group with self-adjoint generator P=iX.
For B_M(r)=M^4[K0(Mr)^2-K1(Mr)^2]/(2pi^2)<0 and
g6(z)=1-z^2/2+z^4/24-cos z>=0, retain the positive spectral form

    P_M[T]=int_0^infinity r dr [-B_M(r)]<T,g6(rP)T>.

The Bessel integral representations in
[DLMF 10.32.9-10](https://dlmf.nist.gov/10.32) give the positive-real
sign and moments used in the retained proof. The base lift and
closed spectral square roots define the model form; an H3 core
controls the tower via its M^-4 bound. The common fourth restoration is
C kappa4 J2/(30pi)-C kappa4 int Rbar|grad T|^2/(90pi).

The all-frequency flat inequality remains

    Kflat(p)>=C p^4/(30pi)[log(p/(mu0 ell))-38/15].     (11)

It follows by differentiating the convergent angular remainder and
using [DLMF 5.9.13](https://dlmf.nist.gov/5.9.E13) for the shifted
digamma functions with positive real part. It is not a high-frequency
asymptotic or a numerical fit. The local Weyl Hessian contributes
-C int(omega+19/60)(Delta T)^2/(30pi), omega=log(r/B), plus its
retained positive scalar square and the bounded lower terms.

For the outer comparison, extend the actual log warp evenly to a
complete line as in the retained proof. The guards (7) imply
H=-partial_Y^2+V=(partial_Y-h'/2)^*(partial_Y-h'/2)>=0,
abs(V)<=22 and Lip(V)<=2000. Since H+24>=-partial_Y^2, operator
monotonicity of log, regularized first, gives for w=H v

    <w,log sqrt(H)w> >= <w,log|D_Y|w>-12<v,Hv>.

[Beckner, Corollary 2, equation (10)](https://arxiv.org/pdf/math/0612486)
at n=1,p=2 gives the remaining full-line logarithmic inequality with
D_H=log 2+psi(1/4)=-gamma_E-2log 2-pi/2>-4. Smooth approximation
and the range w=Hv justify removing the positive regularization; no
unproved gap or auxiliary Dirichlet wall is used.

The scalar comparison's heat weight is <=.8C pi^-1 t^-3 exp(-648t).
[Simon's Theorem 1.1](https://arxiv.org/html/math-ph/9907022) applies
to this continuous bounded potential. Its -Delta/2 convention must
be rescaled to the variance-two bridge for H=-partial^2+V. It then
bounds the heat-kernel difference by
(exp(24t)-1)(4pi t)^-1/2 exp(-s^2/(4t)). The retained integrated
gain-two bound, second contact and BOTH global escape terms follow.
The auxiliary scalar operator is not identified with the actual Dirac
operator; their actual-domain differences are included in the ledger.

At partition scale a=1/100, the outer positive fourth margin is

    C/(30pi)[9-38/15-19/60-4-(log 2)/2] > .0175 C,

and the interior margin is >7.5 C/(30pi). The common-scale identity
-log(mu0 ell)+.5log(x/C)=-log(mu0/sqrt(lambda))>9 is unchanged
by the choice of canonical sphere. Thus the spatial outer logarithm
is canceled in a form inequality, not bounded by a constant. The
gain-two/IMS fourth loss is <.00224 C per outer field, leaving at
least C/(100pi) times the localized J2 sum. No artificial frequency
cutoff is used to make this comparison finite.

## 5. Actual lower-order loss and off-shell contacts

Use E0=E_eff/F_lambda, with the SAME G as the original physical norm:

    E0=beta_T/(2G) int [(3/2)r^2f(xi'+2xi/r)^2
                                             +lambda r^2xi^2]dr.

The genuine TT reconstruction of the preceding note has
||b||<=2||xi||, ||S||^2<=100||xi||^2,
||T||^2<=500||xi||^2, ||S'||^2<=5e4 E,
||T'||^2<=1e5 E, E=||xi'+2Rphysical xi||^2+lambda||xi||^2.
These follow from the same RN convexity, f/r^2<=5lambda,
f'/(2r)<=11lambda and lambda<=f''<=21lambda, all checked by
x>1000, delta<1, d=0. No singular separate bound on the bolt's
A=a'/a is used. The new parameter ratios in section 1 give

    C J0/E0<5.2e-29, C J1/E0<1.7e-27,
    C J0,outer/E0<3.4e-37,
    C J1*/E0<1e-26, C J0*/E0<2e-28.                  (12)

The starred norms include localized fields and the auxiliary half-density.
The small outer measure u^4<=16a^4/25 is retained when differentiating
the partition; replacing it by a global norm loses this margin.
The full recomputed quantum budget in E0 is

| Contribution | Upper negative allowance divided by E0 |
| --- | ---: |
| All finite second contacts, curvature, Hardy and IMS gradient pieces | 1e23 * 1e-26 = .001 |
| Bounded gluing/zero pieces other than actual p0 | 1e25 * 2e-28 = .002 |
| Actual constant contact, including tadpole | 1500e20 * 4e-27 = .0006 |
| Weighted Weyl IMS | <1e-12 |
| Lower Weyl and canonical LLL | <1e-12 |
| New vacuum-volume subtraction | <1e-25 |
| Full off-shell gauge/exponential-chart contacts | <1e-6 |

For the last row, the reconstruction (10), the full local trace anomaly
and angular Ward identity bound the physical gradient by
||E_g||_infinity<=1e16 q^2/(4pi B^4), as in section 3 of
[the seed contact proof](stationary-seed-quantum-form.md). The old
spherical reference potential is NOT subtracted here. The new constant
vacuum functional is smaller than U=q^2/(4piB^4), so the same factor-five
reserve for the full gradient applies. The c_R square has zero gradient
at R=R0, and the effective Einstein-Maxwell comparison is stationary.

The exact affine contact identity remains
I''(h-L_Yg,h-L_Yg)-I''(h,h)=I'[L_Y(2h-L_Yg)]. The return estimates
give a total gauge/chart cost <=3e5(G||E_g||_infinity/lambda)E0,
which is <1e-6 E0 by (3). Full stationarity has not been assumed.
The local square's Hessian is unchanged by this return because
delta R[L_Yg]=L_Y R0=0, and its gradient supplies no chart contact.

For the new vacuum term let abs(Vvac)=3(lambda0-lambda)/(8piG).
In disk variables the volume exponent is 4S+2T. Since
||2S+T||^2<1800||xi||^2, its second variation costs at most
21600(lambda0-lambda)/lambda E0<1e-25 E0. This replaces the old
reference-potential estimate; it is not counted in addition to it.

The sum of the displayed rows is <.003601000003 E0. Therefore

    eta_Q <=.003601000003/F_lambda <.003602<.005.       (13)

This proves (1). All coefficients are in fixed physical norms; the
division by F_lambda is displayed rather than hidden by rescaling.
The local square loss <7.07e-8 E_eff then proves (2). The selected
1/4 target has more than .245 E_eff of remaining allowance.

## 6. Domain scope, remaining thresholds and research decision

Start with genuine TT fields smooth at the bolt and compact at the
outer end, reconstructed with their nonlocal radial b tail retained.
Use the fixed return to (S,T) and the stronger norm

    N^2=E_eff+sum d_n P_Mn[T]
        +C[J2[T]+J2[S+T/3]+int omega(Delta T)^2].      (14)

The positive spectral and weak-derivative components define a complete
space embedded in E_eff. The ordinary estimates above make the Hessian
continuous there, so (1)-(2) pass to its core completion. Smooth bolt
values are not removed by a punctured-bolt H1 cutoff in this fourth-order
argument. For the selected mass/fast tails, the TT-preserving correction
to a far cutoff has coefficient O(L^-2); the retained proof gives
||X^2 delta T||^2=O(L^-1), ||X^3 delta T||^2=O(L).
The full spectral tower is bounded by p^4(1+log_+p), so interpolation
gives O(L^(-1+2epsilon)) for epsilon<1/2. The largest absolute Weyl
term is O(log L/L). Both outer traces T(0)=T'(0)=0 and the limiting
boundary flux are thereby retained. This verifies the domain transfer
of the LOWER INEQUALITY for these tails, rather than assuming traces
alone control logarithmic boundary terms.

It does not prove that (14) is the complete physical metric variational
domain, that the actual form is closed in its own form norm, or that
the coupled constrained operator has the desired Fredholm realization.
The independent current/radial constraints and full boundary action
must still be included. No inverse bound for L in the nonlinear
admission norm is inferred from the positive TT core inequality.

DECISION: the negative quantum/contact TT target passes on the new
reference and stated class; continue to the mixed scalar/shape block
and physical coupled realization. The next numerical requirements are
||E||<=1/4 and ||L^-1 Cmix|| ||Bmix||<=3/8 in the fixed spaces.
The canonical determinant's scalar/mixed derivative is zero in quotient
coordinates, so the remaining mixing is local anomaly/reconstruction
plus the necessary off-shell contacts. These constants are still unknown.
Use the positive TT form obtained here and the scalar energy comparison
to calculate them, rather than restarting the failed uniform source chain.

The preserved-ball targets remain rho=1e-8, Y0<=5e-9, Z0<=.25,
Z2<=5e7 with positive margins for EVERY admission norm. The scalar
particular <1.82e-29 does not bound the complete correction. Reference
kernel control is not a transfer to the unknown finite geometry or a
construction of its interacting state. Compatibility, physical matching,
omitted-operator/noise bounds and all-sector causal stability remain
independent acceptance requirements. The objective remains an actual
controlled stable wormhole, and the current reference is not one.

## 7. Frozen supporting verification

Freeze this note, checker and source review before the scientific run.
Use exact rational parameter/thermal/budget bounds, symbolic RN curvature
and contact identities, and 80/110-digit evaluations of the actual
preparation endpoints and midpoint. Independently compare finite angular
resolvent sums with the digamma formula and an analytic positive tail,
including a complex-root regime; use 32 and 128 modes with fixed tail
enclosures rather than fitting a cutoff. Check the changed common-scale
restoration, its exact local variation and the aggregate E_eff budget.
Keep every run and source snapshot. These computations are diagnostics
for the ordinary transfer argument, not samples of the curved quantum
state, a numerical full-kernel calculation or machine verification of
geodesic flow, finite parts, domains or nonlinear existence. No Lean run
is selected.

Remaining to-do list: mixed scalar/shape estimate and full physical
coupled realization/inverse; finite state, residual, admission and
nonlinear solution; exact compatibility, matching, validity, all-sector
stability, charge/return, opening, transport, preparation, separate
creation and full-result audit.
