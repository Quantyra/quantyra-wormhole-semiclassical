# Nonlinear current influence with unequal returns

2026-10-11, E003/S053. This is a bound on the free-fermion influence of
prescribed electric and rotational sources in the existing LLL circle
model. It controls all terms beyond quadratic, with the fermion cutoff
removed, including zero temperature. Its length dependence is an average,
not the longest return. The actual fluctuating gauge measure, stress and
physical junction remain unresolved. No new physical endpoint is accepted.

## Model and result

Use the dimensionless lengths of [the driven-sea note](driven-fermion-sea.md).
There are q>=2 complex Dirac channels with antiperiodic spatial circles
of lengths L_i=T+D_i. Sources act only on their common interval I of
length T, and may mix channels there. Euclidean time is either a circle
of length beta with antiperiodic fermions, or the real line (beta=infinity).
At zero temperature take a bounded source with compact time support;
boundedness and a finite L3 norm suffice for the estimates. The source
is extended by zero outside I. This is the specified inner-source model,
not the elimination of a physically gauged exterior.

Let gamma_tau,gamma_s be Hermitian Pauli matrices with anticommutator
2 delta_mu,nu. Set D0=gamma_tau partial_tau+gamma_s partial_s. The inverse
on channel i at frequencies omega,p is -i(gamma_tau omega+gamma_s p)
divided by omega^2+p^2. Spatial momenta are p=(2r+1)pi/L_i. At finite
temperature omega=(2n+1)pi/beta; at zero temperature it is continuous.
For fixed q, all L_i are finite, so the spatial gap makes D0 invertible.
No lower gap uniform in q is assumed.

For spin j=(q-1)/2 use

```text
T_q=q(q^2-1)/12, E_0=I_q/sqrt(q), E_a=J_a/sqrt(T_q),
Tr(E_A E_B)=delta_AB, sum_(A=0,x,y,z) E_A^2=4 I_q/q.
```

There are eight real source components phi_mu,A. Define
C_mu=sum_A phi_mu,A E_A, V=i sum_mu gamma_mu C_mu, and X=D0^-1 V.
Here |phi| is the Euclidean eight-component norm, and integral means
integral over time and I. Write

```text
Lbar_beta=(1/q) sum_i min(beta,L_i),
r(phi)=4 sqrt(2) Lbar_beta/(pi sqrt(q)) integral |phi|^3.       (1)
```

At beta=infinity replace Lbar_beta by Lbar=(1/q)sum_i L_i. Then

```text
||X||_S3^3 <= r(phi).                                        (2)
```

If r<1, the analytic logarithm based at phi=0 satisfies

```text
log det_3(I+X)=sum_(n>=3) (-1)^(n+1) Tr(X^n)/n,
|log det_3(I+X)| <= r/[3(1-r^(1/3))].                        (3)
```

S3 denotes the Schatten ideal whose cubed norm sums cubed singular
values. The usual fermion effective action has the opposite sign to
this logarithm. Neither inequality diagonalizes the source with the
return nor restricts the source to equal-frequency vertices.

## Coincident kernel and interpolation proof

The diagonal of |D0|^-3 is scalar in Dirac spin and, on channel i, is

```text
c_i=(1/(beta L_i)) sum_(n,r) (omega_n^2+p_ri^2)^(-3/2).
```

Poisson summation for the spatial antiperiodic heat kernel gives

```text
(1/L)sum_r exp(-t p_r^2)
 =(4pi t)^(-1/2)[1+2sum_(k>=1)(-1)^k exp(-k^2 L^2/(4t))]
 <=(4pi t)^(-1/2).
```

The last sum is alternating with decreasing positive terms and starts
negative. Integrate against (2/sqrt(pi))t^(1/2) exp(-t omega^2) dt.
The spatial momentum sum divided by L is at most 1/(pi omega^2).
Since sum_n omega_n^-2=beta^2/4, c_i<=beta/(4pi). Interchanging the
two circles proves

```text
c_i <= min(beta,L_i)/(4pi).                                  (4)
```

At zero temperature, integral d omega/(2pi) of
(omega^2+p^2)^(-3/2) equals 1/(pi p^2), so c_i=L_i/(4pi)
exactly. These bounds include every spatial and temporal mode.

For completeness, apply the standard Schatten interpolation theorem
to the polar decomposition V=U|V| and the analytic family

```text
F(z)=|D0|^(-3z/2) U |V|^(3z/2), 0<=Re z<=1.
```

On the first boundary its operator norm is at most one. On the second
boundary its Hilbert-Schmidt norm squared is
Tr(|D0|^-3 U|V|^3 U^dagger). Imaginary powers are unitary on their
supports, so this bound does not depend on Im z. Interpolation at
z=2/3 gives

```text
||D0^-1 V||_S3^3
 <= integral Tr_spin,color[(diag(c_i) tensor I_2) |V^dagger|^3]. (5)
```

The polar unitary of D0 can be removed without changing singular values.
One rigorous regularization of this argument restricts |D0| to bounded
spectral windows and |V| to its spectral interval [epsilon,M], setting
it to zero off that interval. The endpoint trace is bounded by the
unregularized positive trace. Interpolate first, remove the windows
strongly, and use lower semicontinuity of Schatten norms. The resulting
bounded product is the stated X. This also handles kernels of V and
real Euclidean time. All mixing occurs in the common interval; no
identification of the unequal circles outside it is required.

To bound (5), put W=diag(c_i), cbar=Tr W/q. The rotational spectral
radius is r_q=sqrt(3(q-1)/(q+1)), hence

```text
||C_mu||op <= (|phi_mu,0|+r_q |phi_mu,rot|)/sqrt(q)
            <=2 |phi_mu|/sqrt(q),
||V||op <=2 sqrt(2)|phi|/sqrt(q).
```

For real coefficients the weighted Gram form of the four E_A is
positive. Its trace is sum_A Tr(W E_A^2)=4cbar. Its largest eigenvalue
is at most this trace, giving Tr(W C_mu^2)<=4cbar |phi_mu|^2.
Tracing over Dirac spin cancels cross terms with different mu. Thus
Tr[(W tensor I_2) V V^dagger]<=8cbar |phi|^2. The positive-operator
inequality |V^dagger|^3<=||V||op V V^dagger then gives

```text
||X||_S3^3 <=16 sqrt(2)cbar/sqrt(q) integral |phi|^3.
```

Together with (4), this proves (1)-(2). Alternatively W<=beta/(4pi)
and the exact unweighted identity Tr(V V^dagger)=2|phi|^2 give
the smaller bound sqrt(2) beta/(pi sqrt(q)) integral |phi|^3.
One may take the minimum of these two bounds at finite temperature.

## Summing the nonlinear terms and removing the regulator

Schatten Holder and the operator norm imply, for n>=3,
|Tr X^n|<=||X||_S3^3 ||X||op^(n-3). Also ||X||op<=||X||_S3.
Summing the resulting geometric series with 1/n<=1/3 proves (3).
For any increasing finite-rank orthogonal projections P_N tending
strongly to identity, P_N X P_N tends to X in S3. This follows by
approximating X by finite-rank operators in S3 and using the ideal
norm contraction for P_N. The series converges uniformly when r<1;
telescoping each power and Schatten Holder show its termwise trace
convergence. Thus the finite determinant after subtraction of its
linear and quadratic traces tends to (3), independently of the
fermion regulator.

The subtraction defines det_3; it does NOT set the physical one- or
two-point function to zero. Those terms still require the correct
renormalization, the gauged exterior determinant and induced quadratic
vertices from [the mouth interaction](gauged-mouth-interaction.md).
The remainder alone is not asserted to be gauge invariant. This
argument includes no nonlinear geometry vertex or higher Landau mode.

## Long returns and the retained time scale

For the existing point-dipole midpoint model u_i=(i+1/2)/q and
D_i=delta f(u_i), the established bound is f(u)<=kappa/u, where
kappa=1+asinh(sqrt(3))/(2sqrt(3)). With H_n=sum_(k=1)^n 1/k,

```text
Lbar <= T+kappa delta (2 H_(2q)-H_q).                        (6)
```

At delta=(log q)^-2 this is T+O(1/log q), even though L_max can be
of order q delta. For a source of bounded normalized amplitude on
a time interval tau=O(delta^-1), (1) is O(log^2 q/sqrt(q)) and
tends to zero. For prescribed rotational amplitude a=Omega ell,
the previous scaling gives a=log^3 q/(2 sqrt(pi) q^(3/2)). Since
phi_rot=sqrt(T_q) a, a fixed bounded profile has phi=O(log^3 q)
and the conservative nonlinear bound becomes O(log^11 q/sqrt(q)).
This too tends to zero. It is a Euclidean source estimate; it does
not itself analytically continue the driven real-time state or bound
stress derivatives. These asymptotics supply no useful finite-flux
physical threshold.

## What the estimate says about a candidate gauge measure

Suppose a regulated, centered Gaussian eight-component source has
pointwise covariance trace at most sigma_q^2 on volume tau T. This
is an assumption to test against the actual quadratic gauge action,
including its constraints, holonomies and exterior load. For a real
Gaussian vector, E|phi|^4=(Tr C)^2+2Tr C^2<=3 sigma_q^4. Hence

```text
E r <= rho_q := 4 sqrt(2)/(pi sqrt(q)) Lbar tau T
                 * 3^(3/4) sigma_q^3.                       (7)
```

If rho_q->0, Markov's inequality gives
P(r>sqrt(rho_q))<=sqrt(rho_q). On the complementary event (3)
is bounded by sqrt(rho_q)/[3(1-rho_q^(1/6))], for rho_q<1.
For example, IF sigma_q^2=O(log Lambda) for a smooth photon regulator
Lambda polynomial in q, and tau=O(log^2 q), then
rho_q=O(log^(7/2) q/sqrt(q))->0. No such covariance bound has been
derived here for the actual coupled mouth/rotation action. No photon
continuum cutoff has been removed by this pointwise argument.

Even a proven (7) gives convergence in probability, not control of the
determinant-reweighted measure. To get convergence of its normalization
one still needs, for example, integral over the bad event of
|det_3(I+X)| dmu_G ->0. On the good event the determinant is uniformly
near one; on the bad event small probability alone is insufficient.
Response and stress require their own derivative estimates. Restricting
the Gaussian to the good event does not construct a physical state.
The remaining decisive question is the actual gauge covariance and
weighted tail, followed by renormalized stress and smooth matching.

There is a useful explicit comparator for (7). Replace the actual
return geometry by a flat infinite Euclidean plane at R=1, freeze the
coefficients of [the rotational action](rotational-current-response.md),
fix transverse gauge and omit global zero modes. Integrate the free
infinite-line fermions only to quadratic order. For one rotation axis,
the transverse polarization of the normalized current is 1/pi. For
nonzero two-momentum p the bare w,V quadratic matrix is

```text
B(p)=p^2 diag(A_w,A_V)+D [[1,-1],[-1,1]],
k=(1,q)/(q+1), phi=sqrt(T_q) k dot (w,V),
C0(p)=T_q k^T B(p)^-1 k,
C(p)=C0(p)/(1+C0(p)/pi).
```

This follows either by a rank-one Gaussian update or direct inversion
of B+(T_q/pi)kk^T. Since B>=p^2 diag(A_w,A_V),

```text
C(p)<=e_rot^2/(p^2+e_rot^2/pi),
e_rot^2=T_q(k_w^2/A_w+k_V^2/A_V)
       =3 T_q (q^2+2)/(K(q+1)^2).                            (8)
```

The transverse projector has trace one. A circular momentum cutoff
|p|<=Lambda therefore gives, for all three rotational axes, the bound

```text
sigma_rot^2 <=3 e_rot^2/(4pi) log(1+pi Lambda^2/e_rot^2).     (9)
```

With K=pi q^4 in the retained g_c^2=q^-2 family, e_rot^2=O(1/q).
For Lambda polynomial in q, (9) is O(log q/q). For rotational fields
alone and tau=O(log^2 q), its contribution to (7) is consequently
O(log^(7/2) q/q^2). These are derived coefficients for this Gaussian
flat reference, not an assumed O(log Lambda) amplitude. The same
calculation does not supply the singlet covariance or the actual
finite-return polarization. It omits the exterior, variable geometry,
global sectors and every nonlinear gauge vertex. Even its cutoff limit
at fixed q has logarithmically divergent point variance. Neither (9)
nor restriction to a finite time window resolves the weighted-tail
condition. Thus this comparator passes the specified local size test;
it cannot replace the physical gauge-state calculation.

## Source comparison and verification

The methods are established: [Britz et al., exact v3](https://arxiv.org/html/2007.12834v3)
define the general regularized determinant in (1.3), including nonnormal
operators; this is the det_3 used above. Their product-anomaly theorem
is not needed. [Hartmann and Lesch, v3, section 4.1](https://arxiv.org/html/2106.02444v3)
give the self-adjoint version and derivative identities. That restricted
theorem is not imported for the generally nonnormal X here. The trace
series and cutoff limit above supply the argument in our conventions.
[Briet et al., section 3](https://doi.org/10.4171/JST/368) state the
standard Schatten interpolation framework. We use it with the explicit
circle kernel rather than import their graphene resolvent estimates.
These are scoped readings, not a completed literature or novelty exit.

The [checker](../checks/nonlinear_current_determinant.py) compares matrix
norms and traces, independent heat-kernel representations, spatial
Fourier integrals, nonnormal finite Dirac determinants and direct
Gaussian Schur inversion of (8). The
[receipt](../evidence/nonlinear-current-determinant.json) records exact
inputs and numerical limits. These tests support the ordinary proof;
they do not prove a coupled gauge state, wormhole stability or shortcut.
The [source record](../evidence/nonlinear-current-source-review.json)
records retrieval and reading scope. No Lean/Lake run is claimed.

Remaining to-do list: derive and control the coupled gauge measure,
including weighted tails and full exterior response; control local
renormalized stress and omitted modes; complete smooth matching, the
physical endpoint and novelty; independently verify, GCP-certify suitable
mathematics and prepare the reproducible manuscript.
