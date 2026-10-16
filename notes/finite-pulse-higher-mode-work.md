# Complete finite-pulse higher-mode work on the static thermal RN reference

2026-10-15, E003/S138, supporting S074. This ordinary proof was
independently reviewed before focused verification. It bounds the complete free fermion work,
including all pulse orders and the specified thermofield correlations,
on the two identical static classical thermal RN exteriors. It does not
give the local stress during the pulse, a corrected exterior, or an
interacting prepared state.

## Frozen operator, endpoint observable and state

Use the exact radial domain and canonical half-densities of
[the coupled-domain note](coupled-fermion-domain.md) and the actual RN
family of [S137](rn-higher-mode-boundary-spectral-control.md). The primary
geometric and paired-boundary inputs are Bintanja--Freivogel--Rolph,
[v2, sections 2.1.1--2.1.2 and appendix A.4](https://arxiv.org/html/2308.00871v2).
Their static model and boundary condition are the imported inputs, not
a claim that their paper proves the finite-pulse estimate below.

Write x=nu U for normalized laboratory time, h_n=H_n/nu, and

```
k=2000, q=k^4, N=2q, eta=1/2000,
1/4<=z_i<=7/25, kappa0=z_i nu,
lambda_n=sqrt(n(n+q)), M_n=lambda_n/sqrt(C),
d_n=2(q+2n), n>=1, delta=nu d<10^-32.
```

The entire nonzero angular tower is included, with the same conventions
for signs, folded exteriors and angular multiplicity as S137. There is
no positive spectral gap at the thermal horizon.

The anti-Hermitian fiber matrix G=i tau_2 tensor sigma_1 has norm one
and commutes with the radial differential expression on the identical
complete half-lines. Its constant rotation V(theta)=exp(theta G/2)
maps the zero-angle boundary domain to the actual angle-theta domain.
It need not preserve the zero-angle domain. With the spinor written as
psi=V(theta)phi, the exact equation on that fixed domain is

```
i partial_x phi=[h_n+rho(x)K]phi,
K=-iG/2=K^dagger, ||K||=1/2, rho=partial_x theta_b.    (1)
```

Bounded perturbation evolution makes (1) well defined. No bounded
operator commutator [h_n,K] is assumed. Its boundary form is retained
below. Let U_n denote the complete interaction-picture pulse propagator,
with the initial phase anchored at source onset x=0. The full propagator
to x_f is exp(-ih_n x_f)U_n, whose final free factor cancels in the
energy expectation. If an earlier initial time is selected, its free
factor must also transport the TFD covariance; it is not discarded
while holding its relative phases fixed. Such transport preserves (19).

The initial covariance is the specified free thermofield purification
of the uncoupled identical RN Hartle--Hawking states at
beta0=2pi/kappa0, with the fixed bundle/charge-conjugation pairing.
In the self-dual particle/hole one-particle representation write

```
C0=P_-+X, P_-=1_(h_n<0), sigma=beta0 nu/2=pi/z_i>11.
```

On each thermofield energy pair with absolute energy E, its difference
from the zero-temperature polarization has the canonical form, up to
fixed unitary pairing and phases,

```
X(E)=[[n_beta(E), exp(i phi) sqrt(n_beta(E)(1-n_beta(E)))],
      [exp(-i phi) sqrt(n_beta(E)(1-n_beta(E))), -n_beta(E)]],
n_beta(E)=1/(1+exp(2sigma E)), X(E)^2=n_beta(E) I.      (2)
```

Consequently X preserves the absolute-energy fibers and
`||X(E)||<=exp(-sigma E)`. This retains the coherent opposite-exterior
correlation. Replacing it by independent reservoirs is unnecessary.
The covariance normalization is the same as S137; any optional
self-dual factor one-half only reduces the conservative bounds here.

The endpoint matter work in units nu is the relative expectation of
the final physical Hamiltonian. Since
`H_theta_f=V(theta_f) H_0 V(theta_f)^dagger`, it becomes in frame (1)

```
W_n/nu = W_vac,n + Tr[U_n^dagger [h_n,U_n] X].         (3)
```

The first term is the positive vacuum excitation energy defined below.
The second trace will be proved absolutely convergent. Formula (3)
does not take a trace of either infinite thermal reservoir energy.
Constant-angle endpoint unitary equivalence is used only on the
specified identical complete RN exteriors. It is not asserted for a
joined evolving interior.

## Static trace-ideal input, with the full angular sum

The companion [shifted-moment proof](rn-shifted-boundary-moments.md)
extends S137's boundary heat comparison to all positive shifts. State
the input explicitly. Put L_n=M_n/(3sqrt(2)delta), b0=24/25,
epsilon_n=exp(-M_n/sqrt(6)), and k_s=sqrt(L_n^2+b0^2 s^2).
For the canonical laboratory boundary trace Gamma and s>=1,

```
mu_1(s)=Tr Gamma(h_n^2+s^2)^-1 Gamma^dagger <=4/k_s,
mu_r(s)<=4^r k_s^(1-2r)
 +2c_r epsilon_n L_n^-2 (b0 s)^(3-2r), r=2,3,4,
c_r=[8/(3sqrt(pi))] Gamma(r-3/2)/Gamma(r).             (4)
```

These are boundary moments, not a measure-wise comparison or a
real-axis resolvent estimate. Positive heat-kernel domination transfers
the scalar moments. The explicit gapless tail is essential.

For clarity the exact trace-ideal consequences used here are

```
a(h)=(1+h^2)^(-7/4), f(h)=h a(h),
||[f(h_n),K]||_1 <= C_n, ||[a(h_n),K]||_1 <= C_n,
C_n=10000 L_n^(-5/2)+100 epsilon_n L_n^-2
                       +1000 sqrt(epsilon_n)L_n^(-3/2),
(1/q) sum_(n>=1) d_n C_n <10^-45.                    (5)
```

The common envelope bounds both commutators separately. It does not
assert an inequality between their actual norms. One can see why (4)
is sufficient without treating the domain commutator as bounded.
The exact weak boundary form gives
`(E-F)K_EF=(-i/2)(Gamma^dagger C_boundary Gamma)_EF`,
where `||C_boundary||=1`. For f_s(h)=h/(h^2+s^2)^2 its divided difference
has numerator

```
s^4-2s^2 EF-E^3 F-E^2 F^2-E F^3
```

over `(E^2+s^2)^2(F^2+s^2)^2`. Writing
`m_j=Tr Gamma h^(2j)(h^2+s^2)^-4 Gamma^dagger`, the trace norm is at most

```
(1/2)[s^4 m_0+2s^2 m_1+2sqrt(m_3 m_1)+m_2].          (6)
```

For a_s(h)=(h^2+s^2)^-2 the analogous bound is
`sqrt(m_3 m_0)+sqrt(m_2 m_1)+2s^2 sqrt(m_1 m_0)`.
Each is a finite sum of boundary Hilbert--Schmidt outer products;
`m_j<=mu_(4-j)(s)` follows by spectral positivity. Finally use

```
(1+h^2)^(-7/4)
 = [Gamma(7/4)Gamma(1/4)]^-1
     integral_0^infinity t^(-3/4)(h^2+1+t)^-2 dt.     (7)
```

Keeping the local mass power in mixed terms before integrating proves
(5). The coefficient 1000 in its last term is a deliberately enlarged
common envelope. The companion supplies the scalar continuum estimate,
integral constants and ordinary infinite angular sum. In particular no
finite angular truncation, hard gap, or inverse zero-energy moment is
used in (5).

## The actual laboratory source has four finite variation norms

Retain the single waveform of
[the fixed-laboratory-source note](laboratory-drive-clock-feedback.md):

```
tau=1/4, T=pi+tau,
S(v)=10v^3-15v^4+6v^5 on [0,1], extended by0 and1,
d0(s)=4eta S(s/tau),
p(s)=rho0(s)=d0(s)-d0(s-pi),
Theta0(s)=integral_0^s p,
Q0(s)=integral_(s-pi)^s d0(v)^2 dv,
g0=(Theta0-Q0)/(4pi eta),
w''+w=g0/2, w(0)=53/200, w'(0)=0,
ds/dx=w(s), theta_b(x)=Theta0(s(x)), rho(x)=p(s)w(s). (8)
```

The fixed reference w is positive on [0,T], as established for this
same source in [the reference enclosure](common-drive-validated-capture.md).
That reference positivity does not certify the varied preparation box.
The present pulse is specified in x=nu U independently of z_i; it is
not retimed separately for each preparation.

Here Theta0 is the moving-window integral of d0, so
`0<=Q0<=4eta Theta0` and `0<=g0<=1`. Also
`Q0'=d0(s)^2-d0(s-pi)^2`. Since
`||p||inf<=.002`, `||p'||inf<=.015`,
`||d0'||inf<=.015`, direct differentiation gives
`||g0'||inf<1`, `||g0''||inf<3`. Duhamel's formula on T<3.4 gives

```
|w|,|w'|<2, |w''|<3, |w'''|<3, |w''''|<5.           (9)
```

The compact extension of p is C2, and D_s^4 p is a finite signed
measure. Elementary quintic bounds, including the four jumps of its
third derivative, give

```
F_j=||D_s^j p||_TV,
F_0=4pi eta<.007, F_1=.004, F_2=.06,
F_3<4, F_4<123.                                     (10)
```

For example TV(S')=15/4 gives F_2=.06. For F_3 use
`|S'''|<=60` on each ramp. For F_4 use `|S''''|<=360`
and both endpoint third-derivative jumps on each ramp, giving
`2(4eta/tau^3)(360+120)=122.88<123`.

Let `a_j=||D_x^j rho||_TV`, where a_0 is the L1 norm and derivatives
through order three are integrable functions. Change variables
`dx=ds/w`. For j>=1 the leading outer factor w in
`(w partial_s)^j(pw)` cancels this Jacobian. In particular the fourth
integrand divided by w is exactly

```
p'''' w^4 +10p''' w^3w'
+p''(25w^2 w'^2+10w^3 w'')
+p'(15w w'^3+30w^2 w'w''+5w^3 w''')
+p(w'^4+11w w'^2 w''+4w^2 w''^2+7w^2 w'w'''+w^3 w'''').
```

The last expression is interpreted as a measure in its first term.
Equations (9)-(10), or the lower-order versions of this expansion,
give the simple sufficient bounds

```
a_0=theta_f=4pi eta<.007,
a_1<1, a_2<1, a_3<100, a_4<3000.                    (11)
```

For example the respective four derivative bounds are below .022,
.36, 36 and 2700. No inverse lower bound for w is needed for these
variation estimates, because of the exact Jacobian cancellation.
The nonzero fourth measure is retained, rather than calling the pulse
infinitely smooth.

## Tagged Dyson series: complete pulse orders in trace norm

This step applies to either b=f(h_n) or b=a(h_n). Put
`c_b=||[b,K]||_1`, and let `K_I(x)=exp(ih_n x)K exp(-ih_n x)`.
The n-th Dyson term on the full ordered simplex is

```
(-i)^n integral_(x_1>...>x_n) product_j rho(x_j)
                                      K_I(x_1)...K_I(x_n) dx.
```

Commuting b through the product inserts exactly one trace-class
`[b,K_I]`, in n possible positions. The other n-1 factors have norm
at most 1/2. Thus the trace norm of the commutator term is at most
`n(1/2)^(n-1)c_b a_0^n/n!`. This proves absolute trace-norm
convergence and gives `||[b,U_n]||_1<=B_0 c_b`.

To obtain the endpoint energy derivatives translate the whole pulse.
Conjugation `exp(ih_n v)U_n exp(-ih_n v)` is the scattering operator
with waveform rho(x-v). Hence its m-th derivative is, up to powers of
i and signs, `ad_(h_n)^m U_n`. The same relation holds for [b,U_n].
For smooth compact pulses integrate by parts using
`D=sum_j partial_(x_j)`. This vector is tangent to every ordering face
`x_i=x_j`, so it creates no internal simplex-face boundary term.
Compact support removes the full-pulse endpoints. Equivalently one
can differentiate the translated waveform directly. The multinomial
sum of the absolute derivative products is symmetric in all x_j.
Its integral on one ordered simplex is therefore 1/n! of its full
product integral; this is where the factorial is retained.

Put `A(t)=sum_(j=0)^4 a_j t^j/j!`. The resulting majorant is

```
sum_(n>=1) n(1/2)^(n-1) A(t)^n/n! = A(t) exp(A(t)/2).
B_m=partial_t^m[A(t)exp(A(t)/2)] at t=0,
||ad_(h_n)^m[b,U_n]||_1 <=B_m c_b, m=0,1,4.          (12)
```

Explicitly

```
B_0=a_0 exp(a_0/2),
B_1=exp(a_0/2)(1+a_0/2)a_1,
B_4=exp(a_0/2)[(1+a_0/2)a_4
 +(1+a_0/4)(4a_1a_3+3a_2^2)
 +(3/4+a_0/8)6a_1^2a_2+(1/2+a_0/16)a_1^4],
B_0<.008, B_1<1.01, B_4<4000.                       (13)
```

The same translation argument without the tagged commutator gives a
bounded first commutator ad_h U_n, because a_1 is finite. In particular
the complete endpoint U_n preserves D(h_n), even though an individual
K need not do so. This justifies the energy-core algebra in (3) and
(18); no intermediate-time graph norm is estimated.

For the actual C2 pulse, convolve rho with a positive smooth compact
mollifier. Its variation norms do not increase. All bounds above hold
uniformly; the bounded-perturbation propagators converge in operator
norm by their L1 source convergence. Trace-ideal lower
semicontinuity against finite-rank tests gives (12), including the
bounded trace-class realization of the fourth commutator. Thus no
fourth L1 derivative or instantaneous graph preservation by K is
silently assumed. No finite-order Dyson truncation is made.

## Vacuum excitation: positive Sylvester inverses and dyadic sums

Write A=h_n restricted to positive energy and B=-h_n restricted to
negative energy. Let T=P_+U_nP_-. Its weighted pair energy is

```
E(T)=||A^(1/2)T||_2^2+||T B^(1/2)||_2^2.
```

The cross block of `[f(h_n),U_n]` is `f(A)T+Tf(B)`.
The fourth h commutator on a cross block is the fourth power of the
positive Sylvester map `L(T)=AT+TB`. If A>=a>0 and B>=0,

```
L^-r(Y)=1/Gamma(r) integral_0^infinity
                   t^(r-1) exp(-tA)Y exp(-tB)dt,
||L^-r||_(S1->S1)<=a^-r.                            (14)
```

The same formula applies to the f-Sylvester map whenever one
spectral band has positive minimum f. These maps commute as spectral
maps. This proves the trace-norm estimates below; bounded scalar
kernel ratios alone would not prove a trace-ideal multiplier bound.

First restrict both energies to [0,1]. There `f(E)>=2^(-7/4)E`.
Trace positivity gives

```
E(T_low)<=2^(7/4) Re Tr[T_low^dagger
                              (f(A)T_low+T_low f(B))]
          <=2^(7/4) B_0 C_n,                        (15)
```

since `||T_low||<=1`. The identity can first be evaluated with positive
spectral cutoffs away from zero and then extended by monotone
convergence. This uses no inverse zero-energy moment.

For high energies use the two overlapping positive covers, i>=0:

```
row_i: e in [2^i,2^(i+1)], 0<=f<=2^(i+1),
column_i: f in [2^i,2^(i+1)], 0<=e<=2^(i+1).
```

Their overlap is harmless because the energy integrand is positive.
On the high band `f(e)>=5^(-7/4)2^(-5i/2)`, while e+f>=2^i.
Apply (14) four times to the block of
`ad_(h_n)^4[f(h_n),U_n]`, and once to the f-Sylvester map. It gives

```
||T_i||_1<=5^(7/4)2^(-3i/2) B_4 C_n.
```

Since e+f<=2^(i+2) and `||T_i||<=1`,

```
E(T_i)<=4*5^(7/4)2^(-i/2) B_4 C_n.
```

Summing both covers gives a factor
`8*5^(7/4)/(1-2^(-1/2))<700`. Combining with (15),

```
E(T)<4B_0 C_n+700B_4 C_n.                            (16)
```

The reverse crossing block P_-U_nP_+ satisfies the identical bound.
The vacuum energy is the positive sum

```
W_vac,n=||h_+^(1/2)P_+U_nP_-||_2^2
              +||(-h_-)^(1/2)P_-U_nP_+||_2^2
        <=2000(B_0+B_4) C_n.                         (17)
```

This definition and bound are continuum statements. They prove finite
relative vacuum energy even without first asserting an unweighted
finite particle number in the gapless infrared. Neither an
intermediate-time energy norm nor an all-frequency small shift operator
has been used.

## Thermal redistribution and coherent thermofield work

The following algebra is the needed relative-trace subtraction. On
energy-compact cores, because f(h)=h a(h),

```
[h,U_n] a(h)
 =[f(h),U_n]-ad_h[a(h),U_n]-[a(h),U_n]h.             (18)
```

For X in (2), the elementary inequality `log(1+E^2)<=E` gives

```
||a(h)^-1 X|| <=1,
||h a(h)^-1 X|| <=1/[e(sigma-7/4)]<1.                (19)
```

Only the left weighted norms are needed: put Y=a(h)^-1X, multiply
(18) by Y on the right, and interpret hY as its bounded extension.
Absolute-energy preservation of the canonical TFD pairing ensures
(19) even when X does not commute with h itself.

Equations (5), (12) and (18) now extend `[h,U_n]X` to a trace-class
operator. Taking its trace after multiplication by U_n^dagger gives

```
|W_X,n| <=(2B_0+B_1) C_n.                            (20)
```

Thus both thermal redistribution and the coherent cross covariance are
included to all orders. In particular no estimate of the form
`Tr |h| exp(-sigma|h|)` on the infinite horizon is used. Such a thermal
volume trace would not be justified by the boundary moment.

There is a bounded self-adjoint endpoint energy-difference observable
even before pairing it with X. For the full Schrodinger-picture
propagator S(x) of h+rho(x)K, the common-domain energy identity and
the zero endpoint values of rho give

```
Delta h=U_n^dagger h U_n-h
       =integral S(x)^dagger rho'(x)K S(x) dx,
||Delta h||<=a_1/2.                                 (20a)
```

This follows by differentiating S^dagger(h+rho K)S on the common
energy domain and integrating; the bounded right side extends the
identity to the Hilbert space. The norm-C1 bounded perturbation
rho K suffices. Then Delta h X=U_n^dagger[h,U_n]X is trace class by
(18), so (20) is a relative-work pairing before any cycling of
infinite reservoir energy traces. Its trace is real: bounded
self-adjoint Delta h and X have trace-class product, and the trace
agrees with that of its reversed product.

One can also recover the thermal term as a smooth spectral energy
cutoff limit. Put b(h)=(1+h^2)^(-9/4). Since b/a and h b/a are bounded,
(18) proves

```
||[h,U_n]b||_1 <=(2B_0+B_1)C_n.
```

The thermal bound gives bounded b^-1 X and h b^-1 X. Choose a real
even smooth compactly supported chi with chi=1 near zero, and let
chi_R(h)=chi(h/R). The exact decomposition is

```
[h chi_R(h),U_n]X
 =chi_R(h)[h,U_n]X+[chi_R(h),U_n]hX.                 (20b)
```

For precision, define the divided-difference map T_chi_R through the
Fourier Duhamel formula for chi_R. That formula, rather than a bare
pointwise kernel bound, proves its S1 multiplier norm is at most
`c_chi/R`, with `c_chi=integral |t chi_hat(t)|dt` in the chosen Fourier
normalization. Right spectral multipliers commute with this map, so

```
[chi_R(h),U_n]hX
 =T_chi_R([h,U_n]b) h b^-1 X.
```

Its trace norm tends to zero at rate O(1/R), while the first term of
(20b) tends to [h,U_n]X in trace norm. Thus
`Tr U_n^dagger[h chi_R(h),U_n]X` tends to the thermal trace in (3).
This establishes an actual smooth energy-cutoff relative limit.

The definition in (3) uses the positive energy-polarization relative
trace for its vacuum term. Increasing positive spectral cutoffs in
the crossing quadratic forms converge to the finite value in (17).
The thermal term is the ordinary trace of the trace-class operator
in (18), hence is independent of any strongly converging uniformly
bounded cutoff inserted there. This is a relative matter-energy
prescription, not an assignment to either individual infinite
reservoir energy. It does not assert convergence of an arbitrary
unmatched subtraction of two spatially regulated reservoir energies
or identify their local stress without the additional Ward matching.

## Normalized result and physical boundary

Sum (17) and (20) with the full angular multiplicity. Equations
(5) and (13) imply

```
|W_higher,U|/(nu N eta)
 <= [2000(B_0+B_4)+2B_0+B_1]/(2eta)
                          *[(1/q)sum d_n C_n]
 <10^-34.                                          (21)
```

For example the displayed upper bounds make the coefficient smaller
than 8.1*10^9, so (5) actually leaves additional slack below (21).
The target 10^-8 is therefore met within this frozen free static
reference problem, uniformly over the retained initial thermal
z_i interval. No source, preparation rectangle or eta was retuned.
The small result is a boundary spectral/trace-ideal conclusion and
does not remove the gapless horizon spectrum or high-angular
real-frequency trapping.

The work is the endpoint matter energy in (3). Constant endpoint
frames have no rho K term, and their physical Hamiltonians are exactly
related by V(theta_f). Accordingly this calculation does not count an
endpoint frame term as particle production.

The intrinsic local source-kinetic contact in the coupled-domain note
has a separate exact endpoint conclusion on this fixed static metric.
For any finite matched coefficient c and smooth coefficient function F,

```
I_ct=c integral F(theta) rho^2 dx,
H_ct=c F(theta)rho^2,
-rho (delta I_ct/delta theta)=partial_x H_ct.
```

Both endpoints have rho=0, so this contact's complete-pulse source
work is Delta H_ct=0, independently of c. The same algebraic endpoint
cancellation holds for a regulated divergent coefficient before its
renormalization. Its instantaneous stress and work density need not
vanish. The cancellation does not cover arbitrary metric/extrinsic
terms, unequal collars, mixed interacting contacts, or a physical
controller with extra degrees of freedom. Those responses still need
their matched action and endpoint ledger; they are not bounded by
(21).

This theorem does not control pointwise driven stress, the radius
source, interactions, state-preparation error, source/controller
contacts, unequal clocks or exteriors, or metric/domain transfer to
the quantum-corrected evolving solution. It does not certify common
laboratory-drive capture or the finite-cut beta remainder. In
particular an energy bound is not substituted for the local stress
bound required by S074/S076. The imported geometric model and standard
trace-ideal methods carry no novelty claim.

## Focused verification contract, frozen before execution

The following bounded contract was frozen before checker execution.
Verification outcomes are recorded separately in
[the receipt](../evidence/finite-pulse-higher-mode-work.json) and
[the log](../evidence/finite-pulse-higher-mode-work.log). There is no Lean,
angular truncation, new source choice or capture trajectory campaign.

Check the rational divided differences and their continuous diagonal
limits exactly; the source chain rule through four derivatives; the
fourth derivative of A exp(A/2), including its factorials; the contact
Euler/endpoint identity; and the conservative moment, beta-integral,
angular, source-jet, B_j, dyadic-series and normalized-target inequalities.
Algebraic residuals must vanish. Exact rational comparisons or explicitly
bounded elementary series supply the uniform inequalities; floating
samples do not replace the continuum or angular proofs.

For the scalar comparison use d=1 and the two synthetic scalar parameter
triples (M,delta,s)=(18,.05,1),(32,.2,3), with nu=delta,
kappa0=(7/25)nu, mu=M/(3sqrt(2)), a=sqrt(nu^2 s^2-kappa0^2).
These test the comparison operator, not an asserted RN parameter family.
For r=2,3,4 compare the exact boundary moment from
`(-1)^(r-1) partial_z^(r-1) g(sqrt(z))/Gamma(r)` at z=a^2
against an independently integrated continuum density (3) of the
companion note, continued above the barrier by hyperbolic-to-trigonometric
continuation. Integrate in p through P=max(100,4mu). For p>=P the density
is at most 8/(3pi), so the omitted positive moment tail is bounded by

```
8/[3pi(2r-1)] P^(1-2r).
```

Use 60 decimal digits and refine to 90 digits. Scaled absolute/relative
comparison tolerances are respectively 2e-8 and 2e-9, with the displayed
analytic omitted-tail allowance explicitly added. Verify the scalar
moment upper bounds themselves. Preserve the raw differences and tail
allowances, rather than fitting the tail to agreement.

Use one deterministic noncommuting six-dimensional diagnostic:

```
h=diag(-3.7,-1.3,-.15,.15,1.3,3.7),
F_jk=exp(2pi i j k/6)/sqrt(6), j,k=0,...,5,
K=F diag(1/2,1/2,1/2,-1/2,-1/2,-1/2) F^dagger.
```

Evolve this finite system with the actual fixed source (8), not a
substitute pulse. Phase coordinates are allowed: evolve
`i S_s=(h/w+p K)S`, `x_s=1/w`, and remove the final free factor to obtain
U. Obtain w independently from its forced oscillator/reference formula.
Use canonical thermofield pairs of opposite equal absolute energy,
z_i=1/4,7/25, and phases 0,.37 in (2). Check unitarity; the covariance
and weighted norms; the exact thermal identity (18); the vacuum crossing
and finite total-work traces; the independently integrated energy-work
identity (20a); and the tagged endpoint bounds using independently
computed finite-matrix norms ||[f(h),K]||_1 and ||[a(h),K]||_1. The
synthetic matrix is not required to satisfy the tiny RN constant C_n;
that constant belongs exclusively to the proved radial comparison.
These finite checks test identities and the noncommuting pulse
implementation; they are not an RN discretization or a numerical proof
of the continuum theorem.

Use DOP853 with rtol=2e-10, atol=2e-12 and refine to rtol=2e-12,
atol=2e-14. Algebra/work/unitarity discrepancies must be below 2e-7
and 2e-8 on refinement, respectively. No tolerance or parameter is
retuned after a failure. Preserve any failed run's raw snapshot and
terminal result. The final receipt must record input SHA-256 hashes
with an explicit UTF-8/LF normalization policy, source precision,
exact command/cwd and exit status; the log must be UTF-8/LF. Any later
text-only provenance refresh must be identified separately from the
execution hashes. The ordinary theorem and the supplementary verification
outcomes retain their separate scopes.

Remaining to-do list: transfer the state/operator/source to the
corrected evolving exterior; control local stress, radius and controller
terms and the signed finite-cut response; complete the S074-S079 gates.
