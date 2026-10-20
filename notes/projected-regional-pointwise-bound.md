# Regional pointwise values of the actual projected correction

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary proof and
verification contract before execution. **Small middle-region values;
cold and higher-jet neighborhood preservation remain unproved.**

The [first derivative result](projected-first-derivative-bound.md) and
[separate scalar/TT value bounds](projected-low-correction-and-mass.md)
give a useful pointwise estimate for the actual correction away from
the bolt. Keeping the two contributions separate improves the estimate
over applying a single embedding to their combined norms.

Retain the actual analytic seed, full Hessian A_m, physical mass column,
source f_ell, compatible projection Pi0, paired domain, fixed clock,
temperature, flux and common subtraction. Write

    h_ell=-A_m^-1 Pi0 f_ell=Rcal(tau,c)=(S,T),
    (S,T)=(tau,0)+(S_c,T_c), 0<=ell<=1,
    t=r/B, s=t-1, x=lambda B^2, g=G lambda,
    dmu=gamma0 t^2 dt, gamma0=4pi beta B^3,
    gradbar=lambda^-1/2 grad.

Here h_ell is the linear correction to the full comparison residual,
not a certified nonlinear solution error. The thermal gap h in the
seed bolt t_h=1+h is a different quantity from h_ell.
Strict estimates below concern 0<ell<=1; at ell=0 the residual and
correction vanish.

For every t>=1+100h the actual continuous radial representatives obey

    |tau(t)|<5.2e-7 ell^2/[t sqrt(s)],
    |(S_c,T_c)(t)|<4.95e-6 ell^2/[t sqrt(s)],
    |(S,T)(t)|<5.5e-6 ell^2/[t sqrt(s)].             (1)

Consequently, on the entire middle and outer region t>=1.001,

    |(S,T)(t)|<1.75e-4 ell^2.                       (2)

The value allowance is below .001 already for s>=3.1e-5. On t>=40
it also obeys `|(S,T)|<5.6e-6 ell^2 t^-3/2`. That last upper bound
does not replace the sharper physical mass and fast-tail asymptotics.
The mass column and outer-fixed return are included in every bound.
No arbitrary additive bolt value or discarded homogeneous mode is used.

## 1. Radial point evaluation from the physical integral norms

For a real or complex radial vector v put

    A_v^2=integral_(t_h)^infinity t^2 |v|^2 dt,
    B_v^2=integral_(t_h)^infinity t^2(F0/x)|v_t|^2 dt.

These are exactly `||v||^2/gamma0` and
`||gradbar v||^2/gamma0`, with no change of norm units. The seed has
F0>0 and F0_t>0 for t>t_h, as retained in
[the weighted York proof](weighted-york-source-conversion.md).
Thus w(t)=t^2 sqrt(F0(t)/x) is increasing. On any interval t>t_h,
radial H1 regularity gives an absolutely continuous representative.
The scalar derivative of |v|^2 costs at most 2|v||v_t|. Cauchy--Schwarz
with the TWO physical weights, followed by monotonicity of w, gives

    |v(t)|^2
       <=2 integral_t^infinity |v||v_t| du
       <=2 A_v B_v/[t^2 sqrt(F0(t)/x)].              (3)

For completeness, A_v finite supplies a sequence t_j tending to
infinity with v(t_j) tending to zero. The displayed integrable bound
on the derivative of |v|^2 permits the fundamental-theorem limit
along that sequence. No arbitrary boundary value at infinity is
assumed. The estimate also applies by local H1 approximation, and
uses the pair norm before its derivative is bounded. It is not a
pointwise estimate for derivatives of v.

The seed polynomial is

    F0=s^2 f(t)-delta/t,
    f(t)=x+2x/t+(3x+1)/t^2+d/t^4,
    delta=(1+h)h^2 f(1+h).

The retained guards x>=1000, |d|<=1 and 0<h<1e-30 suffice here.
For t>=1, f(t)/x>=1: the possible negative d term is covered by
the positive 1/t^2 term. Also f(1+h)/x<=6.002 and 1+h<1.001,
so delta/x<7h^2. If s>=100h,

    F0/x>=s^2-7h^2>=.9993s^2>.99s^2.               (4)

Inserting (4) in (3) gives the uniform regional estimate

    |v(t)|<sqrt(2 A_v B_v)/[.99^(1/4)t sqrt(s)].     (5)

This estimate exposes its cold dependence. It does not extend a
uniform middle-region allowance to t_h, and gives no cold source
weight or small value at the bolt. The earlier radial H1
non-embedding example remains valid.

## 2. Compose the separate actual scalar and TT estimates

For a compatible source, the prior value/derivative proofs give

    ||tau||<43 g Fsrc,
    ||gradbar tau||<4.1e6 g Fsrc,
    ||c||<1664 g Fsrc,
    ||(S_c,T_c)||<=25||c||<41600 g Fsrc,
    ||gradbar(S_c,T_c)||<391680 g Fsrc.              (6)

The scalar gradient comes from the separated outgoing/elliptic
estimate. The TT derivative uses the exact fixed-disk return before
estimating its coefficients. Both are for the actual coupled solution.
The TT source is `fbar_c-Bbar* tau`; it has not been replaced by the
uncoupled input or a mass-free remainder.

For the full comparison residual and its exact Pi0,

    g Fsrc/sqrt(gamma0)<124(2.2e-13)ell^2
                         =2.728e-11 ell^2.          (7)

Only Pi0's established norm-one property in the physical weighted L2
forcing space is used. No pointwise or derivative bound for Pi0 is
inferred. Let d0=2.728e-11. The two numerators in (5) are bounded by

    d0 sqrt(2*43*4.1e6)/.99^(1/4)<5.2e-7,
    d0 sqrt(2*41600*391680)/.99^(1/4)<4.95e-6.

Their sum is below 5.5e-6, proving (1). Since t sqrt(t-1) is
increasing, evaluation at t=1.001 gives (2). The .001 threshold
follows from `(5.5e-6)^2<(.001)^2(3.1e-5)`. For t>=40 use
s>=(39/40)t to obtain the stated outer upper bound.

For comparison, applying (5) directly to the previously rounded
full-pair bounds A=1.2e-6, B=1.3e-4 would give a numerator above
1.7e-5. The smaller result follows from the actual inverse split;
it is not a new norm normalization or an extra source assumption.

## 3. Consequence and remaining neighborhood requirements

The actual correction's middle-region disk-coordinate value is now
quantitatively small. In S167.1 this removes a possible failure of
the .001 VALUE allowance on that region. It does not establish the
entire low coefficient contract: pointwise first jets, cap gradients,
the fixed-chart fourth Holder norm, and all cold requirements remain.
The [fractional canonical theorem](fractional-interior-dirac-response.md)
needs those separate hypotheses before it can be applied to a solver
ball. Small values alone do not meet its P45 background condition.

These estimates do not promote the seed inverse to a candidate
Hessian inverse. They do not control nonlinear remainders, a preserved
global family, the full projected equation or its actual signed moment.
The cold loss in (5) points to the need for source-specific cold and
higher-regularity response estimates, rather than further use of a
uniform H1-to-bolt-value assertion. No conclusion about nonexistence
or instability follows from that limitation.

## Frozen verification contract

Freeze this note and the checker before execution. Reconcile the
fractional response, first derivative and value/mass receipts, with
their exact input identities. Verify the seed factorization and
thermal polynomial bounds, physical norm conversion, weighted
Cauchy identity, split constants, middle/outer thresholds and the
gain over the combined estimate. Use an independent two-component
power-law fixture for the integral point-evaluation inequality. It
is a mathematical norm test, not a physical solution. Preserve every
terminal run and failed source if any. Continuum regularity and the
fundamental-theorem/domain limits remain ordinary arguments. No
Lean/cloud campaign is selected; all S167 tasks stay open.

Remaining to-do list: cold and higher spatial control of the actual
correction; nonlinear family preservation/remainders and projected
full solve; actual signed compatibility; exterior/state, charge and
return; opening, stability, transport, preparation, separate creation
and full-result audit.
