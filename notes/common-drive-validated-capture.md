# Common laboratory drive: finite capture certificate contract

2026-10-15, E003/S136.1-S136.2. This is a bounded certificate attempt for
the stipulated unit-AdS2 free matter/minimal-JT system of
[S135](laboratory-drive-clock-feedback.md). It does not supply the missing
four-dimensional remainder. The receipt determines whether the finite
inequalities below pass; a failed enclosure is not a physical counterexample.

## Frozen question and exact reduction

Keep eta=1/2000, tau=1/4, the quintic S and central reference w of S135,
w(0)=53/200, w'(0)=0. The single laboratory waveform is unchanged.
Cover z(0) in [1/4,7/25], v(0) in [-1/2000,1/2000], with zero prior
source history and every measurable additive residual |r|<=1/1000000
through actual source completion. The targets, fixed before computation,
are z>=3/200 while driving, actual completion by global time 7,
e_complete<=-3/200 and positive final forcing. The source completes at
s=F(t)=pi+tau, not necessarily at global time pi+tau.

Set a(s)=Theta0(s)/(4*pi*eta), f(s)=rho0(s)/(4*eta), C=4*pi*eta,
and U(t)=integral_0^t rho(v)^2 dv/C. Nonnegative rho and exact pi-return
periodization imply

```
max(U,4*eta*a^2) <= I=Q/C <= ceil(t/pi)*U.
```

Indeed Q is the squared L2 norm of the pi-periodization of the truncated
nonnegative rate rho*1_[0,t]. Dropping cross terms gives its lower bound
by integral rho^2. At most ceil(t/pi) translates overlap almost everywhere;
Cauchy gives the upper bound. Also integral_(t-pi)^t D=Theta0, so another
Cauchy inequality gives Q>=Theta0^2/pi. Before the first return I=U exactly.
The periodization proof uses the correlated zero prior rate history;
it does not erase returns or introduce independent reservoirs.

While z,w>0, use the strictly increasing phase s=F(t). Define
L=z*U and J=z*I, and scale l=L/eta. Every physical solution is contained
in the nonsingular differential inclusion

```
z_s = z*v/w,
v_s = (-z^2+a*z/2-eta*lambda*l/2+z*r)/w,
l_s = v*l/w + (4/pi)*f^2*w,
t_s = z/w,
1 <= lambda <= 3.
```

This uses the bootstrap t<7<3*pi and nonnegative L. At L=0 the
product lambda*L is interpreted as zero. Dropping the stronger Cauchy lower
bound enlarges the inclusion and is legitimate. A separate upper preparation
subrectangle may use lambda=1 if its whole-tube completion bound is below pi;
then a first-return contradiction proves that specialization without
assuming the result. Residuals and lambda may be arbitrary bounded
measurable controls; no derivative of either is assumed.

At completion a=1, J=eta*lambda*l and
e=v^2+z^2-z+J, g_final=1-J/z. An upper bound e<0 with z>0
already implies g_final>0. If the residual then vanishes, exact source
completion makes Q constant and the ideal oscillator has the autonomous
all-future positive-root bound. Persistent residuals need the separate
[tail argument](common-drive-area-and-tail.md); none is inferred here.

## Exact reference before numerical certification

Write I_S(u)=integral_0^u S(v/tau)dv and R(u)=integral_0^u S(v/tau)^2dv.
On the three local phases, with alpha=181/462,

```
A: 0<=u<=tau:       f=S(u/tau), a=I_S(u)/pi, Rdiff=R(u),
B: 0<=u<=pi-tau:    f=1, a=(u+tau/2)/pi, Rdiff=u+tau*alpha,
C: 0<=u<=tau:       f=1-S(u/tau),
                     a=(pi+u-tau/2-I_S(u))/pi,
                     Rdiff=pi+u-tau+tau*alpha-R(u).
g0=a-(4*eta/pi)*Rdiff, w''+w=g0/2.
```

These are polynomial forcing pieces of degree at most 11. Taylor
coefficients obey w_(n+2)=(g_n/2-w_n)/((n+2)*(n+1)); beyond the forcing
degree the homogeneous tail is bounded by the geometric majorant of its
factorial denominators. The long affine-forcing middle piece can instead
use its exact sine/cosine solution. Pi is enclosed by rational Machin
alternating sums; elementary trigonometric enclosures use Taylor remainder
bounds. The checker must prove w>0 on every whole phase segment before
dividing by it; a sampled positive reference is insufficient.

## Finite outward-tube proof and bounded stopping rule

The complete preparation cover consists of 30 z-cells of width 1/1000
and four v-cells of width 1/4000, including their exact rational endpoints.
There are 120 closed cells; overlaps at endpoints are harmless.
No preparation cell, source or acceptance threshold may be fitted after
execution. The bounded development attempt uses order-eight center
polynomials and phase steps at most 1/40, with one 1/80 fallback for a
failed estimate. There are at most 12 radius iterations per segment.
The first QR-reboxed development attempt failed by excessive wrapping in
the lowest cell at s=2.125: its center was about .0491 while its z tube
extended below .005. That failed bound is retained. The sole fallback
carries the fundamental matrix basis through interfaces instead of
repeatedly rotating and reboxing it; it uses the finer predeclared step.
This changes the enclosure representation, not the source or problem.

A floating center polynomial c(h) and matrix polynomial P(h) are merely
exact dyadic approximants. For x=(z,v,l)=c+P*xi, direct coefficient algebra
encloses, throughout 0<=h<=H,

```
xi_s=(P^-1/w) * [f_num(c)-w*c_s +(A_num(c)*P-w*P_s)*xi
                    +quadratic(P*xi)+control(c+P*xi)].
```

The quadratic numerator is (delta_z*delta_v,-delta_z^2,
delta_v*delta_l). The only uncertain control numerator is the second
component -eta*(lambda-lambda_nom)*(c_l+delta_l)/2
+r*(c_z+delta_z). Polynomial subtraction precedes interval evaluation,
so cancellation is accounted for without assuming a small floating
integration error. All divisions and inverses require outward enclosures
excluding zero. A QR factor may be any dyadic matrix: its transpose is
never presumed to be its exact inverse. Interval inversion and the whole
interface containment include matrix and center rounding.

For a moving coordinate box |xi_i|<=R0_i+h*d_i, the sufficient strict
face inequality |xi_s,i|<d_i proves invariance by first exit, including
measurable controls via the integral version. All operations enclose the
entire segment, all preparations in the cell, and the full residual
interval. A finite verifier must replay these inequalities and every
interface containment independently from its recorded approximants.
Whole-tube z positivity and t<7 close the coordinate and memory bootstrap.

If any face, inverse, reference positivity, interface, time or endpoint
inequality fails, retain its cell, phase, attempted radii and violated
bound. Such a result leaves full-rectangle capture uncertified, even if
ordinary trajectories meet all targets. A subset certificate must be
labelled as such and cannot replace the original rectangle question.

## Retained bounded outcome

The full rectangle is **not certified**. Both predeclared enclosure
representations failed in the first cell
z(0) in [.25,.251], v(0) in [-.0005,-.00025]. No physical counterexample
or target violation was found. The inclusion and source were not retuned,
and the unvisited 119 cells cannot be credited as covered.

The first development attempt used step 1/40 and repeated QR reboxing.
At phase s=2.125 its proposed z tube was [.0049541,.0947269], crossing
the fixed driving threshold. Its center was approximately .0491.
This attempt predates the independent review's homogeneous-tail,
uncertain-duration and rational-threshold tightenings. Its 85-segment
prefix is retained as a development diagnostic, **not** accepted rigorous
evidence. The first receipt hash, execution source hashes and exact
failed inequality remain in the latest receipt's prior-attempt record.

The corrected fallback used step 1/80 and a carried fundamental basis.
Its first 165 whole segments reached phase s=2.0625; the last whole z
enclosure was [.0377128,.0558483], with prefix global-time upper bound
2.153011. The next segment [2.0625,2.075] failed the moving-box estimate.
The first proposed coordinate face rates
[.126778,.259721,.680731] were below required bounds
[.210118,.431498,1.122606]. Twelve permitted radius iterations did not
close; subsequent radii diverged. These are coordinate enclosure rates,
not observed physical derivatives. Refinement stopped at the declared
limit. No source-completion energy bound was reached.

The replay mode reads the corrected fallback's recorded dyadic center
and basis coefficients and rate boxes without calling their generator.
It reconstructs exact reference enclosures, checks every whole-segment
face, the complete initial cell, all interfaces, driving positivity and
time bootstraps. A local-prefix replay result does not certify completion,
all 120 cells or any of the final targets. Original generation hashes
remain distinct from current replay/reporting source hashes. The raw log
retains commands and terminal statuses; status 2 is the intentional
incomplete-enclosure outcome, not a passing certificate.
The exact corrected fallback generation source is preserved verbatim in
the log appendix; its recovered bytes match execution SHA256
`94dd17eda8040c7fb6cd2b87a0f7a4e6b1ea4baf21dda293242e1d9c5f93ff5f`.
The current verifier also checks reference positivity across the entire
source interval: its whole-segment lower bound is greater than .0339342.
This verifies the reference clock, not every actual clock or capture.
A replay whose linked physical context changed during execution was
rejected by the start/end hash check (terminal status 1); its raw failure
is retained. Final evidence requires replay after those inputs are frozen.

The useful result is the exact positive-prefix finite differential
inclusion and a precise limitation of this bounded interval method.
The numerical feasibility of the full rectangle remains separate from
a theorem. A future proof needs a sharper enclosure/comparison argument
or a justified memory constraint; scanning more preparations does not
repair this failed estimate. Physical remainder and tail conditions are
derived separately in [the physical map](physical-clock-remainder.md)
and [area/tail companion](common-drive-area-and-tail.md).

Remaining to-do list: prove or refute the unchanged full-rectangle
capture targets with a sharper argument; supply the physical clock
remainder and matching assumptions. No full capture certificate is claimed.
