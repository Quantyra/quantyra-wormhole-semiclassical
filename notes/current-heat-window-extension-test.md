# Retained-window test of the existing propagation estimator

2026-10-27, E003/S169. The preceding finite-window calculation established
complete omitted tails but no signed value. The intervening email review
supplied no continuum constant. The full-domain residual goal remains
active/incomplete. This calculation rejects a particular extension of the
existing estimator; it does not evaluate the retained window or complete
the requested signed milestone.

## Actual map and first unsupported implication

Keep E_t, its original dual, psi=z^4, v=-z^4, both original rows, and the
actual paired mass/log realization exactly as in
current-quartic-column-energy-dual.md. Its fixed-column regulator result
is retained. The signed target is Z_window of
current-quartic-positive-time-window.md, on .001<=t<=.05 with the full
thermal kernel. The finite angular part retains its certified infinite
tail. No contact, boundary, cutoff or bolt return is removed.

The existing parametrix P contains an exact paired capped kernel, not
numerically evaluated data. Its proven derivative and propagation constants
apply only at t<=.001. Thus substituting t=.05 into that proof does NOT
produce a certified kernel. The following is a necessary-precondition
test of its positive algebraic majorant recipe, even if every missing
larger-time geometric estimate could be supplied. No operator bound is
asserted outside the accepted time interval.

## Ordinary calculation: an unavoidable term in this recipe

In current_disk_heat_propagation.py, W0=(1-b)^(-1/2) has constant
coefficient one. All subsequent majorant coefficients are nonnegative.
The j=0, derivative-order-zero, n=0 term in cutoff_norm therefore gives

    c_cut(t) >= 1.14e8 t^(-6) exp[-1/(16t)].

This is a LOWER bound on the reported positive UPPER estimator, not on
the true cutoff defect. The exponential is the exact Gamma(1,1/(16t))
term for the existing cutoff at distance .5. Each of the four defect
jets receives this same allowance. Hence bbar>=6 c_cut and the geometric
Volterra ratio used by the code necessarily has

    x(t)=t^4 bbar/140 >= x_floor(t)
         =(6*1.14e8/140)t^(-2)exp[-1/(16t)].                 (1)

The required x<1 fails already above the first root near .00226603.
On [.002,.003], log(x_floor)'=-2/t+1/(16t^2)>0, so interval bisection
certifies that root. The floor at .05 exceeds 5.59912e8. Merely changing
the cutoff distance, keeping this recipe, would need distance >2.06849
for x<1 at .05, outside the existing .625 support and .75 analytic ball.
This is not a proof that a different cutoff argument cannot work.

Stopping the retained evaluation at the root does not rescue the existing
tail estimator. Substitution into the already proved all-positive-time
tail formula gives a one-dimensional coefficient-weighted error report
about 1.001215e-23: over 1.0012e7 times the sufficient 1e-30 application
screen, even with all other errors zero. This is an excessive upper
report; the actual tail may be much smaller. The coefficient is only
the established exact one-dimensional envelope, not a multicolumn map.

## Execute one sharper alternative against the same loss

The geometric summation itself is avoidable. If z=t^4 bbar, retain the
exact convolution coefficients instead:

    Phi(z)=sum_(m>=0) 6^(m+1) z^m/(4m+3)!
          =3(sinh(y)-sin(y))/y^3, y=(6z)^(1/4).             (2)

Indeed, convolving t^3 with another t^3 multiplies by B(4,4), and the
successive beta factors telescope to the displayed factorial. The four
positive jet sum uses precisely bbar of the accepted argument. Thus the
higher correction is bounded by sum(E1_jets)t^3[Phi(z)-1]. This retains
the entire series and is finite for every z>=0.

But its estimator is still excessive. At t=.05, using ONLY the floor
z=140*x_floor gives log10(Phi)>351.075. The first convolution's base
term has positive coefficient at least c_cut/(128*pi), since the
parametrix diagonal coefficient is at least 1/(4*pi). The implemented
mixed-jet bound sums all first-convolution jets before multiplying by
the higher-series factor. Including rank two and the selected factor
two, its reported uniform pointwise-to-L2 error guard therefore exceeds

    4 t^3 c_cut/(128*pi) [Phi(140*x_floor)-1].              (3)

The log10 of (3) exceeds 360.489. Relative to the previous CONDITIONAL
uniform-kernel error guard .002, the excess exceeds 10^363.188. That
guard was one sufficient allocation in the full-window error formula,
not an attained constant or remaining physical margin. No intermediate
kernel constant is compared directly to a physical-output ceiling.

Equations (1)-(3) test the specified majorant extension, not a different
signed estimator or the actual residual. Since all other terms in this
recipe are positive, omitting them makes this rejection conservative.

## Decision, consequence, and remaining calculation

Reject one-step extension of these positive cutoff/propagation majorants
to .05. Exact factorial summation fixes their geometric convergence test
but not their quantitative usefulness. This is a verified estimator
failure, not a construction obstruction. No new signed value, implemented
F_n, numerical application error, or full residual bound was obtained.
No coupled or physical budget closes; all unknown comparisons stay unknown.

Retain the actual paired geometric method and its accepted t<=.001 proof.
The specified alternative for further execution is signed short-step
semigroup propagation: .05 can be reached in 50 steps of .001, inside
the existing time range. This observation is NOT a bound on the composed
metric/time jets or numerical quadratures. It also does not supply the
missing capped-kernel center. No cost or accuracy is certified for it.

The single next calculation is a signed, boundary-compatible short-step
kernel with joint base/metric/time/mixed error enclosures, composed to
produce the complete Z_window L2 enclosure. This is still proposed.
Its errors must enter the original complete residual through the actual
coefficient maps, not through conditional one-column ceilings. Preserve
P32, the original physical norms, closed audit and accepted results.

Remaining to-do list: actual signed application and complete residual;
source/solve/coefficient/application errors, G_t and every coupled budget;
nonlinear admission/contraction, exact compatibility, physical matching,
semiclassical validity/noise and all-sector stability; the full stable
semiclassical wormhole objective.
