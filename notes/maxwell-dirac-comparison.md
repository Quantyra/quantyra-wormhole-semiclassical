# Dirac propagation after source-free Maxwell completion

2026-10-11, S053. The [coupled Maxwell field](coupled-maxwell-matching.md)
has a small potential correction, not merely small integrated field
energy. On the prescribed geometry, with an explicit optical comparison
condition, its Dirac Hamiltonian differs from the preceding
reference Hamiltonian by bounded Hermitian multiplication of norm
`epsilon <= C q p`. At `p=q^(-3/2)` this tends to zero. The comparison
controls full finite-duration propagation, including ground/higher-mode
mixing, and positive-resolution resolvents. It also controls the
compressed higher-mode load and its first two bounded-source derivatives.

The metric is still prescribed. This does not derive its Einstein
backreaction, a renormalized vacuum response, or a physical damping or
shortcut theorem. The optical condition is an additional hypothesis;
it does not follow from the previous Maxwell energy theorem. A direct
upper cap on N/R would exclude an important RN neck regime. The main
comparison below uses radial attenuation and covers the exact RN and
AdS throat pieces without that cap; their physical junction remains open.

## 1. Geometry, gauges and two elementary transformations

Use the geometry, fixed flux, orientation, no-applied-field condition and
notation h,v,f,Lambda_E,Lambda_T of the Maxwell note. Lengths are in units
of separation d=1. Keep `0<p<=1/10`. Let a_E be the exterior correction
to the exact point dipole, with boundary k=h-v; let a_T be the throat
correction to the isolated monopole. Both multiply dphi. Fix the same
flat holonomy and gauge patches as the reference. No electric component
or extra gauge jump is introduced.

In cylindrical coordinates, putting `a_E=r^2 U` gives the exact identity

    (partial_rr-r^-1 partial_r+partial_zz)(r^2 U)
       = r^2 (U_rr+3 r^-1 U_r+U_zz).                       (1)

Thus U is an SO(4)-invariant harmonic function on the exterior of two
radius-p balls in R^5, with centers on the z axis. This is a mathematical
change of variables, not an extra physical dimension. Smooth regular
azimuthal one-forms give smooth U at r=0. Conversely the decaying
five-dimensional Dirichlet solution produces the unique finite-energy
Maxwell perturbation with those traces. It decays as U=O(distance^-3).

On T put `a_T=sin(theta)^2 b`. Its Maxwell equation becomes

    b_xx+w(x)^2 [b_thetatheta+3 cot(theta)b_theta-2b]=0,
    w=N/R.                                                  (2)

The bracket is `Delta_(S^4) b-2b` for zonal b. Hence the axis is a
regular pole, not an extra boundary. Since w>0, the maximum principle
on the finite cylinder [0,L] x S^4 gives

    ||b||_infinity <= max_ports ||h/sin(theta)^2||_infinity. (3)

These assertions can first be applied to smooth data; smoothness of the
coupled solution follows from the preceding transmission theorem. The
constants in (3) do not require an upper bound on the middle profile w.

## 2. A pointwise coercive boundary estimate

Here is the step that the integrated energy bound alone did not provide.
For arbitrary smooth port data h write B=h/sin(theta)^2, and let
`M=max_ports B>0`. For its exterior extension define `F=p^2 U`, so F has
boundary data B in the lifted problem. The positive harmonic function

    W=(p/rho_1)^3+(p/rho_2)^3

is at least one on either sphere and decays at infinity. Comparison gives
`F<=M W`. On the outer sphere rho_i=2p of either disjoint annulus,

    F<=M(1/8+gamma),    gamma=(p/(1-2p))^3 <= 1/512.        (4)

Compare on that annulus with the radial harmonic function H having value
M at p and M(1/8+gamma) at 2p:

    H(rho)=M[8gamma/7+(1-8gamma/7)(p/rho)^3].               (5)

At a point where the inner boundary B reaches M, F=H and F<=H on the
annulus. The derivative into the exterior therefore satisfies F_rho<=H'.
The Maxwell outward flux is not the five-dimensional normal derivative
alone: differentiating a=rho^2 sin(theta)^2 U gives

    (Lambda_E h)/sin(theta)^2 = -2B/p-F_rho
       >= (1-24gamma/7) M/p >= (445/448) M/p.              (6)

On the throat, (2)-(3) imply that at the same boundary maximum its
outward derivative is nonnegative. Thus
`(Lambda_T h)/sin(theta)^2>=0` there. Reflections at the other mouth
preserve sin(theta)^2 and this argument. Applying (6) also to -h proves
the pointwise inverse bound

    ||h/sin(theta)^2||_infinity
      <= (448/445) p ||(Lambda_E+Lambda_T)h/sin(theta)^2||_infinity.
                                                               (7)

There is no multipole truncation or middle-profile derivative in (7).
At a pole the same derivative comparison is performed on the smooth
lifted sphere. Only ordinary harmonic comparison is used.

## 3. The actual forcing and global potential bound

The coupled equation is `(Lambda_E+Lambda_T)h=Lambda_E v-f`.
The rationalized remote-monopole formula implies, uniformly at each port,

    ||v/sin(theta)^2||_(C^k(S^4)) <= C_k q p^2,
    ||f/sin(theta)^2||_infinity <= C q p.                  (8)

For completeness, the lifted exterior Dirichlet data for v are
`v/(p^2 sin(theta)^2)=q/[2S(S+1-p cos(theta))]`, with the reflected
expression at the second mouth. Their fixed C^k norms are O(q), since
S>=1-p. The harmonic maximum principle bounds the lifted solution by
Cq everywhere. An elementary annular barrier supplies the needed normal
derivative estimate. Extend its smooth boundary value G constantly in
rho on p<=rho<=2p. Then Z=U_v-G is zero at p, is O(q) at 2p, and
`|Delta_5 Z|=|rho^-2 Delta_(S^4) G|<=Cq/p^2`. For sufficiently large
fixed A, `psi=Aq[1-exp(-8(rho-p)/p)]` is nonnegative, dominates |Z| at
2p, and satisfies `-Delta_5 psi>=Cq/p^2` throughout the annulus.
Indeed its radial Laplacian is proportional to `-8+4p/rho<=-4`.
Comparison gives |Z|<=psi and hence `|partial_rho U_v|<=C q/p` at p.
Consequently `||(Lambda_E v)/sin(theta)^2||_infinity<=C q p`.
This proves the boundary estimate for these smooth data; it does not
presume a bounded DtN operator on arbitrary L-infinity data.

Equations (7)-(8) give

    ||h/sin(theta)^2||_infinity <= C q p^2,
    ||a_T/sin(theta)^2||_infinity <= C q p^2.              (9)

For the exterior correction k=h-v, its lifted boundary data are O(q).
Comparison with W now yields the global estimate

    |a_E|/r = r |U_E| <= C q r sum_i (p/rho_i)^3
              <= C q sum_i p^3/rho_i^2 <= C q p.         (10)

The second inequality uses r<=rho_i, since both centers lie on the axis.
Equation (10) is the physical orthonormal one-form norm in the exterior.
In the throat, the corresponding optical norm is

    N |a_T dphi|_g = w |a_T|/sin(theta) <= C q p^2 w.      (11)

One sufficient condition would be the global optical cap

    sup_(0<=x<=L) w(x) <= C_w/p.                        (12a)

Here C_w is independent of q,p,L. This suffices for the earlier numerical
test profiles but is too restrictive for the selected physical neck.
For extremal RN, w peaks at `1/(4 r_e)` in physical units, whereas the
port lies far outside r_e. At the leading selected scales,
`p max w=4pi delta q^(3/2)` diverges. We therefore do not assume (12a)
for that neck.

A more useful alternative, independent of any upper cap, is

    y=1/w=R/N>0,       y''(x)<=2/y(x).                   (12b)

At the two Euclidean seams y=p. The function y/p is a positive
supersolution of (2), since `(partial_xx-2w^2)(y/p)<=0`. Apply the
maximum principle to b-M y/p and -b-M y/p, with
`M=max_ports |h/sin(theta)^2|`. Together with (3), this proves

    |b(x,theta)|<=M min(1,y(x)/p),
    w |a_T|/sin(theta)<=M/p<=C q p.                     (12c)

Thus a large w is compensated by attenuation of the nonspherical field.
The optical condition is checked directly on the relevant analytic
pieces. For an extremal RN neck, with areal radius r>r_e,

    f=(1-r_e/r)^2, dx/dr=1/f, y=r^2/(r-r_e),
    y'=1-2r_e/r, y''=2r_e(r-r_e)^2/r^4,
    y y''=2r_e(r-r_e)/r^2<=1/2<2.                       (12d)

For the standard AdS throat `w=delta sec(delta x)`,
`y y''=-cos(delta x)^2<=0`. On a Euclidean collar y''=0. These formulas
are independent symbolic checks, not a solution for a joined Einstein
metric. A physical interpolation must satisfy (12b), or supply another
positive supersolution with endpoint values at least one and bounded
product with w. Arbitrary positive middle profiles can violate both
(12a) and (12b); the old energy theorem alone does not exclude them.

For a unit-charge Dirac field, Clifford multiplication has norm equal to
the one-form norm. The reference full H_0 can use either the sharp
piecewise dipole/monopole connection or the smooth cutoff connection of
[the loaded realization](loaded-operator-realization.md). The latter
differs from the former by O(qp) in its two collars. The corrected
Hamiltonian on the same metric and bundle is therefore

    H_1=H_0+V,        V=V*,       ||V||<=epsilon=C q p.     (13)

Under (12a) the constant may depend on C_w; under (12b) no such cap is
used. Apart from the stated comparison condition and fixed collar
cutoffs, it is independent of flux, throat length and middle-profile
derivatives. Spinor
rescaling is identical for both operators. The Dirac coupling is to A
in the established flux normalization, so no additional g or q is
multiplied into (13). Restore physical units by dividing H and epsilon
by d, and multiplying time by d.

## 4. Full propagation and bounded-source response

The full reference H_0 is self-adjoint by the transmission realization.
Bounded Hermitian perturbation leaves its domain unchanged and makes
H_lambda=H_0+lambda V self-adjoint for 0<=lambda<=1. This follows also
directly by factoring H_lambda-z at |Im z|>epsilon and using the
surjectivity criterion for a symmetric operator at both imaginary signs.
Its resolvent Q_lambda(z) satisfies ||Q_lambda||<=1/nu for Im z=nu>0.
The resolvent identity gives

    ||Q_1-Q_0|| <= epsilon/nu^2.                          (14)

Duhamel's identity, first on the common domain and then by density, gives

    U_1(t)-U_0(t)=-i integral_0^t U_1(t-s) V U_0(s) ds,
    ||U_1(t)-U_0(t)|| <= min(2, epsilon |t|).              (15)

These are full-space bounds. They include every angular mode and all
ground/higher-mode couplings introduced by the nonspherical Maxwell field.
No angular gap or commutation with the old ground projector is assumed.
For any normalized initial wave packet and any common detector
`0<=D<=I`, the probability difference is at most `2 epsilon |t|`.
This is a comparison at the same exterior time coordinate and detector;
it does not assert a travel-time advantage or positive transmission.

Let W_1,...,W_j be fixed bounded Hermitian sources of norms M_1,...,M_j,
added to both operators, and differentiate in their real amplitudes at
zero. They may mix angular modes. Inserting V once and the j sources in
all orders, and integrating lambda from zero to one, proves

    ||partial_1...partial_j (Q_1-Q_0)||
       <= (j+1)! epsilon product_i M_i / nu^(j+2),         (16)
    ||partial_1...partial_j (U_1-U_0)||
       <= epsilon |t|^(j+1) product_i M_i.                (17)

For (17), (j+1)! ordered terms cancel the time-simplex factorial. The
same bound holds for bounded time-dependent sources with uniform norms.
These are operator/source-amplitude derivatives. A metric or length
derivative can be unbounded and is not included. Nor do (14)-(17)
control the discontinuous negative-energy projector of a gapless vacuum,
a continuum determinant, an ultraviolet trace, or renormalized stress.

## 5. What the loaded return comparison does and does not mean

On `L2(E) direct_sum P_H L2(T)`, compress the throat part of V with the
fixed reference P_H and keep the old ground port boundary maps. The
result is a bounded Hermitian perturbation of the reference dissipative
operator. Its domain, maximal accretivity and current identity persist;
the proof uses surjectivity first at large nu and then the unchanged
operator/adjoint current estimate. Hence its Poisson and outgoing
resolvent maps still have norms at most `(2nu)^(-1/2)`.

For the compressed loaded returns T_0,T_1, the exact identity is
`T_1-T_0=-Gamma_out Q_1 V P_0`. The same insertion argument gives

    ||partial_1...partial_j (T_1-T_0)||
       <= (j+1)! epsilon product_i M_i / (2 nu^(j+1)).     (18)

The sources here act on the load, for example exterior multiplication
sources as in the earlier return theorem. This perturbed load retains
all higher-mode mixing, but removes ground throat modes by definition.
It is not the exact elimination of the corrected full throat, which has
distributed ground/higher-mode coupling. Equations (14)-(17) are the
separate full-operator comparison that controls those couplings over
finite times or at positive resolution. Replacing that full problem by
the compressed return alone would omit physical terms.

## 6. Scale consequence and decision

Keep `p=q^(-3/2)`, `g^2=q^-2`, `delta=(log q)^-2`, `nu=delta^2`
and `a=c_a q^(-1/4)`. Then epsilon=O(q^(-1/2)). For j=0,1,2:

| Quantity | Maxwell-completion comparison bound |
| --- | --- |
| Full resolvent/source derivative | O(q^(-1/2) (log q)^(4j+8)) |
| Compressed return/source derivative | O(q^(-1/2) (log q)^(4j+4)) |
| Full evolution/source derivative, abs(t)<=C/delta | O(q^(-1/2) (log q)^(2j+2)) |

With fixed bounded source norms, each is smaller than every fixed power
of delta asymptotically. The compressed-return correction is also
subleading to the preceding q^(-3/8) times logarithmic losses. Thus
source-free Maxwell completion does not spoil that reference comparison
within (12a) or (12b), and the full dynamics are close on the throat/long-path time
scale. The unspecified finite-q threshold and earlier large sampled
reference-trial errors remain unresolved. No numerical fit replaces them.

Continue the physical long-return question. The next physical obstacle is
the metric and local stress, followed by the coupled quantum state and
renormalized full response. A small one-particle propagator perturbation
does not establish interacting damping, finite-disturbance survival, or
the sign of a physical vacuum-energy derivative.

## 7. Verification and established inputs

Run `checks/maxwell_dirac_comparison.py --output
evidence/maxwell-dirac-comparison.json`. Exact algebra checks both lifts,
the annular barrier, its constant, the RN/AdS optical inequalities and
the scale powers. A scalar neck family with arbitrarily large p max(w)
checks attenuation independently against an exact solution. Coupled exterior
multipoles and continuous throat boundary-value solutions test the
potential and maximum-principle bounds with resolution refinement.
Independent finite-dimensional unitary/resolvent calculations test
Duhamel and source insertions, including off-diagonal channel mixing.
Those finite calculations supplement the proof; they are not a
three-dimensional Dirac simulation, validated numerics or a vacuum test.

The five-dimensional stream-function substitution is established. See
Erhardt, Wahlen and Weber, [2202.01754v2](https://arxiv.org/abs/2202.01754v2),
Section 3.1, for its use to remove the axis singularity. Their fluid
bifurcation results are not imported here. Bounded self-adjoint
perturbations and the second resolvent identity are reviewed in
[Teschl, second edition](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe2.pdf),
Theorem 6.4 and Lemma 6.5. The source receipt records the precise reading
scope and remaining dependencies. The comparison proof supplies the
geometry-specific potential estimate; neither the substitution nor
Duhamel's formula is claimed as a new physical result.

Remaining to-do list: derive and control the Einstein metric and local
stress, then the state-dependent renormalized response including actual
ground/higher-mode feedback, exterior determinant and discarded tail;
finish source/proof/novelty review and the consequential endpoint;
independently verify it, certify suitable mathematics exclusively on GCP,
and prepare the reproducible manuscript.
