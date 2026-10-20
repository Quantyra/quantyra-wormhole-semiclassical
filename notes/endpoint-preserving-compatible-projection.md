# A compatible projection that preserves the source at both endpoints

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary proof and
verification contract before execution. **A fixed compatible solver
projection with norm below 1.006; no nonlinear solution or cold metric
bound is certified.**

The original orthogonal projection remains valid, together with all
previous correction bounds. This note constructs an alternative whose
subtracted source is smooth and compactly supported in the open
exterior. It leaves the source unchanged on neighborhoods of the bolt
and infinity. Its norm is below 1.006 in the SAME physical forcing
space. The actual projected linear residual is then confined to that
compact region, with its signed coefficient retained.

## 1. Actual adjoint, physical Hilbert space and fixed quantiles

Use the analytic seed, full retained Hessian A_m, mass column, state,
paired domain, clock, temperature, flux and subtraction of
[the matched inverse theorem](matched-residual-inverse-composition.md).
Let H have the physical norm

    ||f||_H=Fdisk(f)=||W(f_S,f_T)||_L2(dmu),
    W=(r/B)^2, dmu=4pi beta r^2 dr.

The actual continuous range functional is

    L0(f)=integral chi0 dot f dmu,
    r0=W^-2 chi0, n=r0/||r0||_H,
    eta(f)=<n,f>_H=L0(f)/||r0||_H,
    Ran A_m=ker eta, Pi0=I-n eta.                    (1)

Here chi0 is the actual coupled adjoint, including its TT component,
not a reference second-order scalar wave. The
[phase construction](stationary-scalar-phase-moment.md) identifies
its nonzero oscillatory outer coefficient and original domain. The
retained local mixed elliptic regularity gives its smoothness on each
compact subset of the open exterior. These facts make n smooth there
and nonzero in H. Its global norm is finite by the retained range
functional, not by assuming the generalized adjoint is physical L2.

Write t=r/B and t_h=1+h. Define the exact cumulative H energy

    C(t)=integral_(r_h)^(Bt) W^2 |n|^2 dmu.

It is continuous, nondecreasing, starts at zero and tends to one.
Define a and b as its first crossing points for 1/200 and 199/200,
respectively. Absolute continuity and the endpoint limits imply

    t_h<a<b<infinity, C(a)=1/200, C(b)=199/200.        (2)

This definition also works if C has flat intervals. It uses the actual
seed adjoint and makes no arbitrary selection of a preparation or
phase parameter. No numerical values for a,b are asserted.

Let a_-=(t_h+a)/2 and b_+=2b-t_h. Choose the following fixed smooth
cutoff, so its dependence is fully specified:

    H(v)=0 for v<=0, H(v)=exp(-1/v) for v>0,
    S(v)=H(v)/[H(v)+H(1-v)],
    theta(t)=S((t-a_-)/(a-a_-))*S((b_+-t)/(b_+-b)).

The denominator never vanishes. The exponential has every right
derivative zero at zero, so S and theta are smooth. We have
0<=theta<=1, theta=1 on [a,b] and theta=0 outside (a_-,b_+).
This cutoff is fixed at the seed; it is not recomputed for a new
forcing or nonlinear candidate.

Put

    v=theta n, alpha=<n,v>_H, b2=||v||_H^2,
    p_cut=v/alpha, Pi_cut f=f-p_cut eta(f).           (3)

The probability mass in [a,b] is .99, hence

    .99<=alpha<=1, alpha^2<=b2<=alpha.                (4)

The lower b2 inequality is Cauchy--Schwarz; the upper one is
theta^2<=theta. Thus (3) is defined without estimating a possibly
tiny signed moment of a hand-chosen bump. Changing the normalization
or sign of chi0 leaves Pi_cut unchanged.

## 2. Exact range, endpoint preservation and norm

Since eta(p_cut)=1,

    Pi_cut^2=Pi_cut, Ran Pi_cut=ker eta=Ran A_m,
    eta(Pi_cut f)=0, Pi_cut f=f when eta(f)=0.        (5)

The source direction p_cut is smooth and compactly supported in
[a_-,b_+]. Therefore Pi_cut f equals f on neighborhoods of BOTH
endpoints, including every local derivative or source asymptotic
coefficient that f has there. For nonsmooth H sources this equality
is understood almost everywhere, or distributionally. The statement
does not assert that the cutoff avoids the entire named cold region
t<=1.001, nor that its transition has a quantified cold-weighted norm.

Here is the sharp Hilbert norm calculation. Decompose
p_cut=n+w with w perpendicular to n. For f=z n+y, y perpendicular
to n, Pi_cut f=y-z w. Consequently

    ||Pi_cut||=||p_cut||=sqrt(b2)/alpha
                           <=1/sqrt(alpha)<=sqrt(100/99)<1.006,
    ||Pi_cut-Pi0||=||w||=sqrt(b2/alpha^2-1)
                           <=1/sqrt(99)<.101.       (6)

The operator norms follow from the row map (z,y)->y-z w; its square
on n-perp is I+w tensor w. This proof does not assume that Pi_cut
is orthogonal. Its slightly larger norm is charged explicitly.

For comparison, Pi0 need not preserve even the outer zero germ of a
compact source. Take f=v in (3), which has eta(f)=alpha>0. Outside
its compact support, Pi0 f=-alpha n. The actual adjoint has a nonzero
oscillatory outer coefficient, so n does not vanish on any whole
outer tail. Thus Pi0 f is not compactly supported. This is a property
of the actual projection, not a claim that its established weighted
L2 mapping or physical fast inverse is wrong. We make no assertion
about the value of chi0 precisely at the bolt.

For a finer source Banach space Y continuously embedded in H and
containing compact smooth sources, p_cut belongs to Y whenever its
interior regularity and endpoint conditions allow such sources.
Then Pi_cut is bounded with

    ||Pi_cut||_(Y->Y)<=1+||p_cut||_Y ||eta||_(Y*).    (7)

The number 1.006 is an H norm bound ONLY. Equation (7) is not a
numerical Holder, logarithmic, cold-weighted or conormal bound;
those require the norms of this fixed p_cut. Endpoint equality alone
does not bound the cutoff's interior derivatives.

## 3. Transfer to the actual comparison correction

For the same full retained comparison residual f_ell define

    h_cut,ell=-A_m^-1 Pi_cut f_ell,
    h_0,ell=-A_m^-1 Pi0 f_ell.

The existing inverse is applied on precisely the same closed range.
The new source remains in its established forcing space. Its norm
is at most 1.006 Fdisk(f_ell); no differentiated source estimate is
needed. All bounds below retain the actual mass column and fixed
return. Strict inequalities concern 0<ell<=1; both corrections vanish
at ell=0. Composing (6) with the previous generic inverse estimates
and the actual residual bound gives

    ||h_cut,ell||/sqrt(gamma0)<1.21e-6 ell^2,
    H1(h_cut,ell)/sqrt(gamma0)<1.32e-4 ell^2,
    |h_cut,ell(t)|<5.6e-6 ell^2/[t sqrt(t-1)],
                                   t>=1+100h,
    |h_cut,ell(t)|<1.8e-4 ell^2, t>=1.001.          (8)

Here gamma0=4pi beta B^3, with beta the fixed physical period.
The value and H1 norms in (8) are the unweighted disk-pair norms;
the gradient uses the retained gradbar=lambda^-1/2 grad. The last
two lines use [the regional pointwise proof](projected-regional-pointwise-bound.md).
Since p_cut-n is in Ran A_m, the exact difference is

    h_cut,ell-h_0,ell=A_m^-1(p_cut-n)eta(f_ell).

The same estimates with the .101 bound in (6) yield

    |h_cut,ell(t)-h_0,ell(t)|<5.6e-7 ell^2/[t sqrt(t-1)],
                                      t>=1+100h,
    |h_cut,ell(t)-h_0,ell(t)|<1.8e-5 ell^2, t>=1.001. (9)

These are estimates for two different well-defined linear
corrections. The previous statements about h_0,ell are unchanged.
Neither correction is identified with a nonlinear solution error.

The new linear residual is exactly

    f_ell+A_m h_cut,ell=p_cut eta(f_ell).            (10)

Thus the linearized full equations hold on neighborhoods of the
bolt and infinity; the unevaluated residual is supported in the
fixed compact region. The inverse is global and nonlocal: (10) does
NOT say that the metric correction is supported there or unchanged
near the endpoints. Its cold values and higher jets still need proof.

## 4. The complete physical equation is unchanged

For the actual nonlinear fixed-density gradient mathfrakF(g),
using this fixed seed projection gives the exact equivalence

    mathfrakF(g)=0
       iff Pi_cut mathfrakF(g)=0 AND eta(mathfrakF(g))=0. (11)

After the projected nonlinear equation, the remaining residual would
be p_cut eta(mathfrakF(g)). Localizing it does not remove its signed
coefficient, evaluate it or turn it into a physical external source.
The action, state, clock and endpoint data are unchanged. The new
projection is a solver choice; the intersection (11) has the original
physical solutions. The projected branch alone can depend on that
choice, and has not yet been constructed for either projection.

The result provides a concrete alternative for S167 that preserves
source endpoint germs with a controlled physical forcing norm. To use
it in the full nonlinear family, still control p_cut in the chosen
interior source norm and prove the actual inverse's cold and higher
spatial bounds, outer/radius/local maps and family preservation. No
S167 checkbox, signed compatibility or physical completion follows.

## Frozen verification contract

Freeze the note and checker before execution. Reconcile the current
regional correction, matched inverse and actual phase receipts and
their input hashes. Verify the projection/range algebra, sharp norm
and difference formulas, quantile mass and every transferred constant.
Use an exact three-dimensional matrix with a cutoff direction to
check endpoint preservation and the failure of the orthogonal version
to preserve that component. An independent continuous probability
example checks the two quantiles and mass accounting. These fixtures
do not calculate the actual adjoint, quantile radii or physical moment.
The adjoint domain, smoothness, quantile existence and operator passage
remain ordinary proofs. Retain all terminal runs and failed sources.
No Lean/cloud campaign is selected.

Remaining to-do list: cold and higher spatial inverse/source norms,
including the fixed projection direction; nonlinear preservation and
projected full solve; actual signed compatibility; exterior/state,
charge and return; opening, stability, transport, preparation, separate
creation and full-result audit.
