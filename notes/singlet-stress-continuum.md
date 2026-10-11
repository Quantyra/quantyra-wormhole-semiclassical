# Leading absorption with the charged singlet retained

2026-10-10, S052/S053. **In the specified Gaussian return model, the
leading smoothed throat-stress absorption survives the spherical photon.**
The estimate is uniform in the bounded singlet mass profile, including
profiles whose size grows with flux. It controls the quadratic stress
observable, not merely a boundary return matrix. It does not establish a
physical wormhole decay law, the renormalized real response, a prepared
four-dimensional state or novelty of the complete research target.

This closes the continuum absorption question left by
[the coupled singlet calculation](schwinger-return-feedback.md). Its
finite-channel counterexample remains valid: pointwise replacement of
`q` by `q-1` is generally incorrect. The theorem instead takes a joint
large-channel, shrinking-band limit, with an explicit error estimate.

## 1. Statement and unchanged model

Use the same `q` scalar current channels on a common throat `[0,T]` and
exterior return intervals of lengths `D_j`, with channelwise matching.
Only the common singlet `v=(1,...,1)/sqrt(q)` has a throat potential
`m(s)^2`, where `m` is real and bounded for each finite system. The
exterior is free. The clock family remains
`H(phi)=exp(phi) V+H_ext`, where `V` is the complete throat energy,
including the singlet mass. All positive-frequency oscillators start in
their vacuum; zero/charge sectors are held fixed as in the preceding note.

Measure frequencies in units of the reference inverse throat scale, so
`T=pi`. Let `u_j=(j+1/2)/q`, `j=0,...,q-1`, and initially take

`D_j=kappa delta/u_j`,
`kappa=1+asinh(sqrt(3))/(2sqrt(3))`.

Let `w>=0` be even and `C_c^4(R)`, supported in `[-1,1]` with
vanishing derivatives at the support boundary, with `int w=1`, and put
`w_h(Omega)=h^(-1) w((Omega-1)/h)`, `h=delta<1/4`. Define

`R_m(q,delta)=int_0^infinity w_delta(Omega) rho_m(Omega) dOmega`.

Here `rho_m` is exactly the vacuum stress spectral measure from equation
(5) of the preceding note, with the finite singlet mass included. There
is a constant depending only on `w,T,kappa` such that, for sufficiently
small `delta` and `q delta>=1`, uniformly over all such `m`,

`|R_m/(q delta)-C_*| <= C_w [1/(q delta^5)`
` + delta log(2+q delta^2) + log(2+q delta)/(q delta^2)`
` + delta^2 + 1/(q delta)]`,                              (1)

`C_*=kappa [8log(2)-13pi^2/24]/(3pi^2)>0`.

Consequently `R_m/(q delta)->C_*` if `delta=q^(-a)` with
`0<a<1/5`. The same limit holds with the actual point-dipole lengths
`D_j=delta f(u_j)`; in that case add the uniform geometric `o(1)` from
[the proved dipole Jacobian limit](collective-clock-response.md).
No rate for that additional term is asserted here.

The frequency window is of the order of the previous model's proposed
damping width. It is a measurement/smearing scale, not a damping term
inserted into the Hamiltonian. Every finite graph remains conservative.

## 2. A first-order representation of the stress pairs

Let `H_s` be the direct sum of scalar `L2` spaces on all throat and
exterior edges. Its form domain `V_0` is the direct sum of the `H1`
circle domains with the original channelwise mouth continuity. Introduce

`d_m f=(f', m v* f_throat)`

from `V_0` into `H_g = H_s direct_sum L2(0,T)`. It is closed and
densely defined, and `h_m=d_m* d_m` is the scalar graph Hamiltonian.
The extra component is energy bookkeeping, not an additional physical
field. On `H_s direct_sum H_g` define the self-adjoint operator

`D_m=[[0,d_m*],[d_m,0]]`,
`Dom D_m=Dom d_m direct_sum Dom d_m*`.

For completeness, the adjoint of this block operator has the same
domain: testing separately against its two domain components requires
membership in `Dom d_m**=Dom d_m` and `Dom d_m*`. This proves
self-adjointness without a bounded-mass-size assumption uniform in `q`.
Its nonzero eigenvalues are `+/-omega_a`, with normalized eigenvectors

`e_a^+ = (u_a, d_m u_a/omega_a)/sqrt(2)`,
`e_a^- = (u_a,-d_m u_a/omega_a)/sqrt(2)`.

`X` multiplies the scalar and derivative edge components by the throat
indicator and is the identity on the auxiliary mass component. Thus
`X` is an orthogonal projection independent of `m`, and

`<e_a^-,X e_b^+> = (A_ab-B_ab/(omega_a omega_b))/2`,

where `A,B` are the throat mass and energy matrices in the previous
note. The ordered-pair weight is therefore exactly

`|B_ab-omega_a omega_b A_ab|^2/(8 omega_a omega_b)`
` = (omega_a omega_b/2) |<e_a^-,X e_b^+>|^2`.              (2)

This is the key conversion of an unbounded stress vertex to a bounded
projection with explicit frequency weights. The large auxiliary kernel
of `D_m` is harmless: all formulas below vanish near zero in both
spectral variables, then remove that cutoff by a separate estimate.

## 3. Cutting only the singlet at the two mouths

Restrict `d_m` by the two scalar constraints
`L f=(v*f_throat(0),v*f_throat(T))=0`, and call it `d_m,K`.
It is still closed and densely defined. Construct the self-adjoint
`D_m,K` from it by the same block formula.

The restriction `S_m=D_m|{L f=0}` is a common closed, densely defined
symmetric restriction of `D_m` and `D_m,K`. The map `L` is continuous
and onto `C^2` in the graph norm of `D_m`: the derivative controls
endpoint traces and test functions near the two mouths prescribe them
independently. For nonreal `z`, `D_m-z` is a graph-norm isomorphism
onto the Hilbert space. Hence `Ran(S_m-z)` is closed of codimension two,
and `dim ker(S_m*-z)=2` (use the conjugate parameter in the range).
The difference of the two resolvents has range in this defect space.
Self-adjointness then gives the explicit bounds

`rank[(D_m,K-z)^(-1)-(D_m-z)^(-1)]<=2`,
`||(D_m,K-z)^(-1)-(D_m-z)^(-1)||_1<=4/|Im z|`.            (3)

No edge-length lower bound uniform in `q`, nor bound on the size of `m`,
enters (3). Every individual finite graph has positive edge lengths.
This is the usual finite-boundary-extension argument, supplied here
to identify exactly which two conditions are changed.

After the cut, rotate the common throat into its singlet and neutral
components. The scalar operator splits into the massless reflecting
graph from the previous note and an isolated Dirichlet singlet throat
interval carrying potential `m^2`. The energy representation and `X`
respect this splitting. On the isolated interval `X=I`, so its
negative-to-positive matrix elements vanish. Cross terms with the
reflecting graph also vanish by orthogonality. Therefore the **stress
pair response of the cut system is independent of `m`**, even though
its total spectrum and vacuum energy are not. This avoids taking an
uncontrolled infinite-mass limit.

## 4. A trace bound for a two-frequency observable

For `F` compactly supported in the negative-positive quadrant and away
from both axes, define

`I_D(F)=sum_(lambda,mu) F(lambda,mu)`
`                    * |<e_lambda,X e_mu>|^2`.

Only finitely many nonzero eigenvalues participate. For the stress
calculation `lambda<0<mu`. Apply the Cauchy-Green functional calculus
in both variables. One can use the explicit almost-analytic extension

`F_tilde=(1+i y partial_x)(1+i v partial_t) F(x,t) chi(y)chi(v)`,

with even `chi=1` on `[-1,1]`, zero outside `[-2,2]`, monotone on
each transition. For one variable the coefficients multiplying
`f,f',f''` in its antiholomorphic derivative are

`b_0=i chi'/2`, `b_1=-y chi'/2`, `b_2=i y chi/2`.

Their integrals `int |b_j(y)|/|y| dy` are bounded by `1,1,2`.
By (3), expanding the difference of
`R_D(z) X R_D(w) X` in its two factors gives trace norm at most
`8/(|Im z| |Im w|)`. Integrating yields

`|I_(D_m)(F)-I_(D_m,K)(F)|`
` <= (32/pi^2) sum_(i,j=0)^2 ||partial_x^i partial_t^j F||_L1`.

The same bound applies at `m=0`. The cut responses coincide for the
stress functions under consideration, so

`|I_(D_m)(F)-I_(D_0)(F)|`
` <= (64/pi^2) sum_(i,j=0)^2 ||partial_x^i partial_t^j F||_L1`. (4)

The resolvent products individually need not be trace class because of
the auxiliary zero kernel. We integrate their **difference** in trace
norm, where (3) applies, and only then identify the finite spectral
traces. We do not interchange a divergent individual trace and integral.
The two-variable formula follows either by two scalar Cauchy-Green
identities and spectral calculus, or by smooth separated-function
approximation with the displayed derivative norm. `C_c^4` functions
follow by approximation from smooth ones in these norms.

Choose a fixed `C^4` cutoff `b(r)=0` for `r<=1/2`, `b(r)=1` for `r>=1`,
with `0<=b<=1`, and insert

`F_(epsilon,h)(x,t)=(-x)t b(-x/epsilon)b(t/epsilon)`
`                  * w_h(t-x)/2`.                        (5)

For `epsilon=h`, all derivatives through order `(2,2)` have supremum
at most `C_w h^(-5)`. More precisely the order `(i,j)` costs at most
`h^(-1-i-j)`. Their support is a strip of area `O(h)` in a fixed
negative-positive rectangle, since `t-x` lies in `[1-h,1+h]`.
Thus (4) is at most `C_w h^(-4)`. This support factor is essential:
the bound is for the quadratic stress statistic, with constants
independent of `q`, the lengths and the mass profile.

## 5. Removing the low-frequency cutoff

Let `n_m(epsilon)` count positive scalar frequencies at most `epsilon`.
The positive potential implies `h_m>=h_0` as quadratic forms.
The free kernel has dimension `q`; a nonzero singlet potential removes
at most one constant mode. Min-max therefore gives

`n_m(epsilon)<=n_0(epsilon)+1`.                            (6)

The massless case is immediate. By (2), Bessel's inequality for the
positive eigenvectors, and `||X||=1`, the part omitted from a positive
output band below `Omega_max` is at most

`(Omega_max ||w||_infinity epsilon/h) n_m(epsilon)`.        (7)

This includes both choices of the small pair frequency; the factor
`1/2` in (2) compensates that union bound. It also handles arbitrarily
long returns with many soft modes, rather than assuming a uniform gap.

For either inverse-tail or point-dipole lengths `D_j<=C delta/u_j`.
Take `epsilon T<=pi`. A free circle contributes below `epsilon` only
if `u_j<=C delta epsilon/pi`. Summing its two-polarization count
`2 floor[epsilon(T+D_j)/(2pi)]` over those indices, and using the
half-integer harmonic sum, proves

`n_0(epsilon)<=C q delta epsilon [1+log(2+q delta epsilon)]`. (8)

If the index set is empty the bound is immediate; otherwise its
cardinality is at most a constant times `q delta epsilon`, so no
spurious `q epsilon` term appears. Combining (4)-(8), setting
`epsilon=h=delta`, and restoring the `pi` in `rho` gives

`|R_m-R_0|/(q delta)`
` <= C_w [1/(q delta^5)+delta log(2+q delta^2)+1/(q delta)]`. (9)

This comparison is uniform in the singlet mass, including the finite
Schwinger profile on any regular truncated throat.

## 6. Discrete channels, geometry and the limiting coefficient

For one free circle let `L=T+D`, `Omega_n=2pi n/L`, `n>=2`.
The independently reproduced stress weight can be written

`W_n=c_n Omega_n^2 sin^2(Omega_n T/2)`,
`c_n=(n^2-1)/(6pi^2 n^3)`.

For fixed `n`, `W_n w_h(Omega_n)` has variation at most `C_w/(hn)`
as a function of `u`: change variable to the monotone `Omega_n(u)`.
The same holds for the actual dipole once the support lies in its
strictly decreasing tail. Outside that tail the weight is zero for
small `delta`, since `f(u)<=sqrt(5)/u` excludes low-frequency roots.

The midpoint sum minus `q` times the integral is bounded by total
variation. The largest discrete harmonic in the band is `O(1+q delta)`
because `D_max<=C q delta`. Summing the variation bounds gives
`O(log(2+q delta)/h)`. The continuum harmonics above that threshold
contribute only `O(1)`, using the summable `q delta/n^2` bound.
Thus the discrete-to-continuum error, divided by `q delta`, is

`O(log(2+q delta)/(q delta h)+1/(q delta))`.                (10)

For the inverse tail the continuum density is exactly, on the band,

`rho_cont(Omega)/(q delta)`
` = (kappa/3) Omega^2 sin^2(Omega T/2)`
`   * sum_(n>=2) (1-1/n^2)/(2pi n-T Omega)^2`.             (11)

For `T=pi`, its value at `Omega=1` is `C_*`. The series is smooth
locally, so averaging against the normalized even kernel changes it
by `O(h^2)`. Equations (9)-(11) prove (1). The actual dipole Jacobian
converges uniformly over this band and all its resonant roots, as
proved in the earlier clock note; it supplies the additional geometric
`o(1)`. No unverified channel derivative is assumed in this transfer.

## 7. Consequence, verification and predecessors

In the earlier fixed-coupling scale comparison, `delta` is proportional
to `q^(nu-1)`. The condition `q delta^5->infinity` is `nu>4/5`,
consistent with the other model scale tests. At `nu=15/16`, the
boundary comparison contributes relative power `q^(-11/16)`, the
infrared term is `O(q^(-1/16) log q)`, and the discrete-channel term
is `O(q^(-7/8) log q)`. These are asymptotic model bounds, not small
errors certified at a physical finite parameter point. Running
couplings and higher angular photon exchange remain outside this model.

The result preserves the leading **positive-band absorption coefficient**
under the included interaction. It does not establish a retarded pole,
an exponential decay interval, the renormalized contact/real terms, or
nonlinear survival. Those distinctions prevent promoting (1) directly
to the proposed physical damping rate. The exterior field, actual
mouth matching, total metric vertex and quantum state still need control.

The finite-boundary and functional-calculus methods are established.
[Post, sections 2 and 6](https://arxiv.org/html/0712.1501v1) describes
first-order graph domains and boundary resolvent formulas; the present
two-constraint proof does not import its equilateral spectral reduction.
[Nonnenmacher's spectral-theory notes, sections 6.2.1-6.2.3](https://www.imo.universite-paris-saclay.fr/~stephane.nonnenmacher/enseign/Spectral%20Theory%202023%20chap1-8.pdf)
give the Cauchy-Green construction and order-one extension argument.
The needed portions and proofs were inspected. The new application is
supporting research toward a physical result; final novelty acceptance
and the full citation/proof review remain open.

Run `python checks/singlet_stress_continuum.py`. The checks independently
test the energy-representation identity, the singlet cut, the two-variable
Cauchy formula, cutoff derivative scaling, soft-mode counts and the
finite-channel continuum coefficient. Source/dependency hashes and all
numerical scope limits are retained in its receipt. The analytic proofs
above establish uniformity; finite numerical cases do not prove it.
In particular, the large-channel numerical sums check the **free**
discrete-to-continuum step; their `q delta^5` values are below one and
do not demonstrate a small interacting error bound. The mass-uniform
interaction conclusion rests on (2)-(9). The energy-representation
checks use a finite-element contraction for `X`, converging to the
continuum projection; they do not numerically certify extension domains.
No Lean is invoked. Later formal certification remains GCP-only.

Remaining to-do list: control the renormalized real response, physical
vertex, endpoint matching, remaining gauge interactions and finite-time
state; finish primary/fallback novelty selection and the consequential
physical result; independently verify, GCP-certify suitable endpoints
and prepare the reproducible manuscript.
