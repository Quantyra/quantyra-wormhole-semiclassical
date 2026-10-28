# Short-step composition in the required diagonal norm

2026-10-27, E003/S169. Previous turn: progress, rejecting the one-step
positive estimator. This argument bounds the analytic bias of a complete
short-step parametrix composition. Its signed kernel is NOT numerically
evaluated here. The full T_t-(I+T_t)F_n goal remains active/incomplete.

## 1. Actual operator, domain and output

Use the original fixed seed spinor Hilbert space H, including both paired
copies, the fixed thermal period/spin structure, outer projector and regular
disk bolt. With the accepted intertwiner, A=i Dsharp is self-adjoint and
A^2=L is the original paired positive square. This does not make the
unpaired Euclidean Dirac self-adjoint. Put V=Dom(A), with its graph norm.
The original mass/log realization and complementary adjoint remain intact.

For the fixed smooth direction v=-z^4, the accepted actual-geometry bounds
are |v|<=1 and |grad v|<=1/2. Scalar multiplication E_s=exp(-s v/2)
preserves V and the complementary domain, by the distributional product
rule and commutation with the boundary projector. The regular bolt remains
a disk. Define A_s=E_s A E_s and L_s through its closed form on V, rather
than freezing the square's composed operator domain. For real s this is
the original metric family. Only this fixed-direction family is used;
no nonlinear neighborhood or independent phase variation is admitted.

Let K_s(t)=exp(-t L_s). Let P_s(h) be the COMPLETE boundary-compatible
parametrix in current-disk-heat-propagation.md, including its exact capped
kernel, cutoff transitions and regular-bolt transport. It maps into the
actual boundary realization; its full remainder was already proved for
0<h<=.001. P is not P32 and no polynomial degree is changed.

For fixed integer N and t in [.001,.05], define K_N(t)=P_0(t/N)^N.
All products integrate over the actual full domain. The target kernel
jets are (1,partial_s,D_t,partial_s D_t), at s=0, D_t=t partial_t.
For the heat-Hessian density, take -2 D_t partial_s of the diagonal spin
trace and retain the completed a0/a2 contacts. The same contacts cancel
in the kernel approximation error. The exact UV logarithm, common
regulator, full thermal spectrum, angular weight and local restorations
remain as in current-quartic-positive-time-window.md.

## 2. A quantitative bound for exact evolution between errors

This step is needed before replacing repeated coarse Schur bounds by
exact evolution. Let C_v=[A,v], ||C_v||<=g=1/2. Then

    A_s=exp(-s v)(A-s C_v/2).

Use the analytic form q_s(u,w)=<A_conj(s) u,A_s w> on V. For |s|<=r,
with d=exp(2r)-1, its quadratic difference satisfies

    |q_s(u,u)-||Au||^2|
      <=d||Au||^2+r exp(2r)g||Au||||u||
                         +r^2 exp(2r)g^2||u||^2/4.          (1)

Take r=.1, sigma=.25, theta=.1, and complex time zeta with
|zeta/t-1|<=sigma. Young's inequality gives

    Re[(zeta/t)q_s(u,u)] >= c||Au||^2-omega||u||^2,
    c=1-sigma-(1+sigma)d-theta >.37324,
    a=(1+sigma)r exp(2r)g,
    omega=a^2/(4theta)+(1+sigma)r^2 exp(2r)g^2/4 <.015523.   (2)

The form bounds give closedness and sectoriality on the common dense V;
the shifted real part is equivalent to the graph norm. Coefficients are
holomorphic in s. These verify the hypotheses of
[Vogt--Voigt, Theorems 0.1 and 3.1](https://arxiv.org/html/1707.00923),
not a smooth-boundary heat theorem. Lax--Milgram gives the shifted
m-accretive realization. Its energy inequality implies
||K_s(zeta)||_(H->H)<=exp(omega t). The argument also holds on slightly
larger parameter disks, by the strict c>0 margin. Cauchy's inequalities
therefore apply to the four jets on the displayed disks.

For ANY operator norm define the weighted jet norm

    J(B)=||B||+r||partial_s B||+sigma||D_t B||
                                      +r sigma||partial_s D_t B||.

Leibniz gives J(BC)<=J(B)J(C). For every exact evolution interval
0<t<=T=.05, the preceding bounds give

    J_(2->2)(K(t)) <= C=4 exp(omega T) <4.003106.            (3)

This bound is for the actual infinite-dimensional operator. It does NOT
assert a complex-parameter pointwise heat-kernel estimate.

## 3. Recover the diagonal norm, rather than infer it from (3)

From the accepted complete P and F tables, reconstruct each P jet's
pointwise constant A_j/h and both Schur constants B_j. For E=P-K the
accepted pointwise constants are e_j h^3. Its Schur constants are
d_j h^4: the first convolution uses B_P B_F/4, with the same six mixed
Leibniz terms as the accepted pointwise argument. Add the entire higher
Volterra allowance, bounded by the positive jet sum times x/(1-x),
x=.001^4 bbar/140. For the Schur convolution this is conservative because
its first extra beta factor is 1/280, smaller than the used 1/140.

For weights w=(1,r,sigma,r sigma), kernel Cauchy--Schwarz gives

    e(h)=J_(2->2)(E)<=e0 h^4, e0=sum w_j d_j,
    d(h)=max(J_(2->infty)(E),J_(1->2)(E))<=d0 h^(7/2),
    d0=sum w_j sqrt(e_j d_j),
    B(h)=B0 h^(-1/2)+d0 h^(7/2),
    B0=sum w_j sqrt(A_j B_j).                              (4)

B(h) bounds both smoothing norms of K(h), by P-E. For an exact run
K(kh), k>=1, split off K(h) and apply (3) to the remaining run, including
the differentiated factors; its smoothing norm is <=B(h)C. At a base
time t, the split fractions are fixed when taking D_t. The identity
operator is used for an empty run, with jet norm one.

Expand (K(h)+E(h))^N-K(Nh) into all words containing errors. Group the
exact K factors between errors BEFORE taking norms. With z=1+C e(h),
the resulting kernel sup-norm jet bound is

    U_N(h)=B(h)^2 C(z^N-1)
           +2 d(h)B(h)C z^(N-1)+d(h)^2 C z^(N-2).         (5)

The three terms bound words with no end error, one end error, and two
end errors; binomial summation overcounts some words, conservatively.
For example a word with m errors and two nonempty end runs is bounded
by B(h)^2 C^(m+1)e(h)^m. An end error uses its smoothing norm d(h).
Thus (5) is a kernel estimate, not an inference from operator norm.
The continuous positive-time kernel representatives supplied by the
parametrix and its convolutions permit restriction to the diagonal.

Rank two and the selected factor two give the ORIGINAL density bias

    rho_N <=4 U_N(T/N)/(r sigma).                          (6)

This holds uniformly on the whole time window: expanding z^N-1 shows
every power of h in (5) is nonnegative after B(h)'s h^(-1/2) factors.
The physical z interval has length <1, so the same guard bounds its L2
norm. No new derivative requirement on energy test fields is introduced.

## 4. Numerical consequence and the still missing application

Interval arithmetic gives e0<131.451, d0<37.115, B0<1.426561. At N=50,
the density-bias report is about .034363; at N=200 it is .00214449.
At N=500 it is <.000343016, below the earlier CONDITIONAL .002 allocation.
Unlike the rejected one-step extension, every local bound here is used
inside its proved interval. This certifies an analytic approximation bias
for the complete composed map, not a numerical quantum evaluation.

Charge it to the original column as

    d_bias=alpha CLt p W60 rho_500, p=1/(B^2 beta),          (7)

using the accepted coefficient envelope, angular/time weight, CLt and
tail debits. Alpha is only the exact one-dimensional Ritz envelope.
The complete multicolumn coefficient/source/solve errors remain unknown.

For an explicit numerical accuracy test, suppose each numerical step
has weighted jet errors <=eta_inf in kernel sup norm and <=eta_row in
BOTH Schur norms. These include its kernel, internal quadrature,
interpolation and roundoff errors; they are not established here.
The N-fold numerical product may use different step approximations.
Let h_min=.001/N, h_max=.05/N and

    Z=1+C e(h_max), C_P=C Z^N,
    B_max=B0/sqrt(h_min)+d0 h_max^(7/2),
    D_P=C[B_max Z^N+d(h_max) Z^(N-1)],
    eta_d=sqrt(eta_inf eta_row), Z_num=1+C_P eta_row.

Applying the same word argument with exact P runs bounds additional
numerical kernel-jet error by

    U_num=(D_P^2/C_P)(Z_num^N-1)
           +2 eta_d D_P Z_num^(N-1)+eta_d^2 C_P Z_num^(N-2).

The additional density error is 4 U_num/(r sigma). Execute this test for
eta_inf=1e-11, eta_row=1e-16, N=500, then include the separate weighted
time/spatial integration error 9e-32 through the existing full-window
formula. Passing this CONDITIONAL test is not evidence that any of these
numerical errors has been attained.

Attempted implementation: the available code still evaluates only the
transport MAJORANTS and capped-kernel path bounds. No routine returns
the signed P_s(h) and its four jets with the required errors. The exact
capped Brownian expectation is a defined continuum object, not numerical
data; substituting its majorant would be a different kernel. Consequently
no signed center, actual eta bounds, numerical product or F_n was obtained.
The complete quantum-application certificate remains unresolved.

Decision: retain short-step composition; its analytic bias is quantitatively
affordable in the original norm. The first missing numerical implication
is now explicit: enclose that actual paired kernel and its four jets in
both pointwise and Schur norms, rather than only sampled scalar pairings.
The single next calculation is that signed step evaluation and its
propagation into Z_window. No further step-count optimization is selected.

No coupled admission, nonlinear existence, exact compatibility, physical
matching, validity/noise or stability claim follows. P32, original physical
norms, the closed audit and all accepted results remain unchanged.

Remaining to-do list: signed kernel/application and complete residual;
all source/solve/coefficient errors, G_t and coupled/admission budgets;
nonlinear existence, compatibility, matching, semiclassical validity/noise,
all-sector stability and the full stable semiclassical wormhole objective.
