# Higher-mode boundary spectral control on the thermal RN reference

2026-10-15, E003/S137, supporting S074. A shifted spectral moment of the
actual classical near-extremal RN-AdS radial operator bounds its entire
nonzero angular tower at the frozen laboratory scale. This is a free
static-reference result. It does not solve the finite drive, establish
the physical stress remainder, or replace a quantum-corrected exterior
by the classical reference.

## Reference, physical source and comparison selected before checking

Use one massless four-dimensional Dirac field on two identical purely
magnetic RN-AdS exteriors, zero electric potential/chemical potential,
the fixed spin/bundle identification and actual paired radial projector of
[the coupled-domain note](coupled-fermion-domain.md). The source is
theta=2 atan(h_b), with the same single laboratory waveform as
[S135](laboratory-drive-clock-feedback.md). The primary inputs are
[Bintanja--Freivogel--Rolph v2](https://arxiv.org/html/2308.00871v2),
equations (2.3)-(2.12), (2.23)-(2.27) and appendix A.4. Their LLL
truncation is not used to remove the nonzero angular modes here.

Write U for physical exterior Killing time and x_*<=0 for its tortoise
coordinate, zero at conformal infinity. With extremal sphere radius b,

```
X=b^2/L4^2, C=1+6X, r_m^2=b^2(1+3X),
f_ext=1+r^2/L4^2-2b(1+2X)/r+r_m^2/r^2,
f=f_ext-D/r, D=2G DeltaM>=0, dx_*/dr=1/f.
```

The outer horizon is r_+=b(1+s_h); kappa0=f'(r_+)/2. Use the retained
classical thermal family, k>=2000 integer,

```
q=k^4, N=2q, eta=1/2000, X(1+3X)=pi*k^6,
b^2/G=X*k^10, phi_r=b^3/(G C), nu=N*eta/phi_r,
kappa0=z_i*nu, 1/4<=z_i<=7/25.
```

The [RN finite-cut comparison](rn-ads-finite-cut-matching.md) gives
s_h<=8*kappa0*b/C<=112*eta/(25*k^9)<1/1000. Define

```
d=b/sqrt(C X), delta_lab=nu*d=2eta*sqrt(C/X)/(k^6 X),
X>=k^3, C/q<7/k, delta_lab<=2eta*sqrt(7)/k^9.
lambda_n=sqrt(n(n+q)), M_n=lambda_n/sqrt(C),
d_n=2(q+2n), n>=1, M_1>=18.
```

The last mass inequality follows from X<sqrt(pi/3)k^3 and
C<(123/20)k^3, giving M_1^2>20k/123>18^2. This M is not a global
spectral gap: the thermal horizon has gapless scattering states.
The zero angular modes are excluded and remain the retained light sector.
The normalized laboratory profile and eta stay fixed throughout this
auxiliary geometric family; nu supplies the U-time unit at each k.
At each fixed k, all preparations use one U-time waveform because nu
is independent of z_i. At the physical benchmark k=2000, eta also
equals 1/k; that equality is not used to retune the source for other k.

The normalized laboratory time is x_lab=nu U=eta u, and the frozen
waveform is theta_b(x_lab)=Theta0(s(x_lab)), ds/dx_lab=w(s).
Its total angle is theta_f=4pi*eta=pi/500, monotonically increasing.
No initial thermal tangent is substituted for the evolving global clock.

## Actual boundary domain and the square

After the canonical radial half-density rescaling, a positive angular
channel has H_n=-i alpha partial_(x_*)+beta W_n,
W_n=lambda_n sqrt(f)/r, alpha=sigma_2, beta=sigma_1.
The sign-reversed angular channel has the same estimates. In the folded
two-exterior space, A=tau_3 tensor alpha and the potential is
1 tensor beta W_n. At zero source angle the allowed boundary projector
is Pi0=(1+tau_3 tensor beta)/2, of rank two in this radial/flavor space.

This projector is the physical radial gamma projector, not a guessed
condition on a scalar component. Covariantly, gamma^r anticommutes with
the radial velocity gamma^0 gamma^r and commutes with each angular
Hamiltonian matrix gamma^0 gamma^A: the latter contains two Clifford
anticommutations. Thus

```
{chi0,A}=0, [chi0,1 tensor beta W_n]=0,
(1-Pi0)psi(0)=0, Pi0 partial_(x_*)psi(0)=0             (1)
```

for the domain of H_n^2. The second condition follows by imposing the
first also on H_n psi; the angular potential's forbidden projection
vanishes. There is no Robin W_n term. A genuine four-dimensional Dirac
mass would have a different Clifford relation and is outside this proof.
Infinity is a regular endpoint in the canonical conformal problem;
the horizon is limit point, without an added reflecting boundary.
This is the static Hamiltonian square. The driven Euclidean operator
D_E^dagger D_E has the additional theta'-dependent mixed Robin term of
the coupled-domain note; the present result does not remove that term.

The square is -partial_(x_*)^2+W_n^2+Y W_n', where Y is a constant
Hermitian involution. Constant nonzero paired angles are related by the
collar rotation of the coupled-domain note. On two complete identical
RN half-lines that rotation extends through the horizon ends and gives
the same boundary trace estimate. This extension is not asserted across
the joined interior of an evolving wormhole.

## The finite-temperature potential bound

Put h=f'/2-f/r, so W_n'=W_n h. For u=b/r, A0=1+3X and
H_X=X+2Xu+A0 u^2, the exact extremal identities are

```
b h_ext=u(1-u)(2A0 u-1),
H_X-Cu^2=X(1-u)(1+3u)>=0,
|h_ext|<=C/(4b), h_ext^2<=C f_ext/r^2.
```

Indeed |2A0 u-1|<=C and u(1-u)<=1/4 for 0<=u<=1.
Since h=h_ext+3D/(2r^2), on r>=r_+>=b,

```
h^2-C f/r^2 <= 7 C D/(4b^3)+9D^2/(4b^4).            (2)
```

The horizon equation gives D=b*s_h^2 J(s_h)/(1+s_h),
J(s)=C+4Xs+Xs^2; its temperature gives
kappa0*b>=C*s_h/(1+s_h)^3. For s_h<=1/1000,
D<=1.001 Cb s_h^2, (1+s_h)^6<1.007, hence
C D/b^3<1.009 kappa0^2 and D^2/b^4<1.01e-6 kappa0^2.
Equation (2) is therefore <(4kappa0/3)^2. Consequently

```
|W_n'| <= W_n^2/M_n+(4kappa0/3)W_n,
W_n^2-|W_n'| >= W_n^2/2-kappa0^2.                    (3)
```

For the second step use M_n>=18 and
(4kappa0/3)W<=(4/9)W^2+kappa0^2. Shift by nu^2. Since
kappa0<=7nu/25, the square's matrix potential obeys the pointwise bound

```
W_n^2+Y W_n'+nu^2 >= W_n^2/2+(24nu/25)^2.            (4)
```

Thus the negative tail near the nonextremal horizon is absorbed at the
actual laboratory spectral scale, uniformly over all retained modes.

There is a useful fixed collar, without a new physical wall. Let
a=4kappa0/3 and y=-x_*>=0. Integrating W'<=W^2/M+aW backwards from
W(0)=M/d gives

```
1/W(-y) <= exp(a y)*d/M+[exp(a y)-1]/(a M).
```

Here a*d<1/1000. On 0<=y<=d, the power-series inequalities
exp(v)<=1/(1-v) and (exp(v)-1)/v<=1/(1-v) for v<=1/1000
give W(-y)>=M/(3d). The a=0 limit is continuous. Equations (1)-(4)
compare the shifted square to the scalar Neumann half-line operator

```
B_n=-partial_y^2+a0^2+mu_n^2 1_[0,d],
a0=24nu/25, mu_n=M_n/(3sqrt(2)d).                    (5)
```

## Matrix heat domination and exact boundary moment

For the free mixed D/N problem, its matrix kernel is
K_N Pi0+K_D(1-Pi0). The positive scalar kernels satisfy 0<=K_D<=K_N,
so its matrix operator norm is at most K_N. A matrix potential bounded
below by a scalar v has ||exp(-t V)||<=exp(-t v). The Trotter product,
kernel composition and these two inequalities dominate the norm of the
matrix heat kernel by the scalar Neumann heat kernel with potential v.
The potentials are continuous bounded matrices for each fixed n, so the
kernel limit follows from the bounded-potential expansion; the step in
(5) is allowed by scalar positive-potential approximation. No operator
monotonicity of the exponential or inverse square is assumed.

At the boundary the matrix kernel is supported on rank-two Pi0, and is
positive as a spectral diagonal. Thus its trace is at most twice the
scalar kernel. Integrating t exp(-t Q_n), Q_n=H_n^2+nu^2, proves

```
tr_boundary [Q_n^-2](0,0) <= 2 [B_n^-2](0,0).         (6)
```

This heat argument is needed: a form inequality alone would not imply
the asserted order of inverse squares.

Let k_n=sqrt(mu_n^2+a0^2) and D0=k_n sinh(k_n d)+a0 cosh(k_n d).
The exact scalar Green function with one endpoint at the boundary is

```
G(0,y)=[cosh(k_n(d-y))+(a0/k_n)sinh(k_n(d-y))]/D0, y<=d,
G(0,y)=exp[-a0(y-d)]/D0,                            y>=d.
```

The resolvent identity gives `R_n=[B_n^-2](0,0)=integral_0^infinity G^2dy`.
Since a0<=k_n and D0>=k_n exp(k_n d)(1-exp(-2k_n d))/2,

```
R_n <= 2/[1-exp(-2k_n d)]^2
             *[1/k_n^3+exp(-2k_n d)/(a0 k_n^2)]
    <= 4[1/mu_n^3+exp(-2mu_n d)/(a0 mu_n^2)].        (7)
```

For the last inequality M_n>=18 makes the prefactor smaller than four.
The tail term is retained; it accounts for the gapless thermal horizon.

## Complete angular trace at the physical laboratory scale

Define the positive boundary counting measure, including both energy
signs and both exterior allowed traces,

```
N_n(E)=nu^-1 tr[Gamma 1_(|H_n|<=nu E) Gamma^dagger], E>=0.
```

Gamma is the canonical boundary trace. The nu^-1 comes from the physical
temporal half-density psi_lab=psi_U/sqrt(nu). Spectral positivity,
(6)-(7) and (E^2+1)^2/(H_n^2/nu^2+1)^2>=1 on that band give

```
N_n(E) <= alpha_n(1+E^2)^2,
alpha_n=612 delta_lab^3/M_n^3
       +150 delta_lab^2 exp[-sqrt(2)M_n/3]/M_n^2.      (8)
```

The coefficients include the rank-two folded boundary. The exact sphere
counting function is D(Y)=sum_(M_n^2<=Y)d_n<=K Y,
K=2C(1+1/q), since the first m levels have 2m(q+m+1) states.
Integration against this positive measure, with c=sqrt(2)/3, yields

```
sum d_n/M_n^3 <= 3K/M_1,
sum d_n exp(-cM_n)/M_n^2
                <=K exp(-cM_1)[1+2/(cM_1)].          (9)
```

Use K/q<15/k, M_1>=18, exp(-cM_1)<1/4000 and
1+2/(cM_1)<5/4. The entire tower, without truncation, satisfies

```
A_tower=(1/q)sum d_n alpha_n
 <=1530 delta_lab^3/k+(45/64)delta_lab^2/k
 <230000*eta^3/k^28+20*eta^2/k^19 <1e-68, k>=2000.
(1/q)sum d_n N_n(E) <=1e-68 (1+E^2)^2.              (10)
```

This is an integrated spectral bound, valid at every frequency. It is
not a pointwise density bound or a hard excitation gap.

## A specified TFD source term and what remains unbounded

Take the free thermofield purification of the uncoupled identical
classical RN Hartle-Hawking states at beta0=2pi/kappa0, with the selected
bundle/spin map. Its cross-covariance singular values are
1/[2cosh(beta0 H/2)]. Relative time evolution adds phases but preserves
the bound. This explicitly specifies the reference state; an interacting
prepared state with gauge constraints has not been identified with it.

Let C_boundary=i A G be the Hermitian norm-one matrix in the exact
boundary work P=(theta_dot/2) integral v^dagger C_boundary v. It is
off-diagonal in the two exterior flavors. A product HH state has zero
unperturbed cross expectation. In the specified TFD, trace Cauchy gives

```
|<integral v^dagger C_boundary v>_0|/q
 <= integral exp(-s E) d[(sum d_n N_n(E))/q]
 <= A_tower[1+4/s^2+24/s^4] <1.05e-68,
s=beta0*nu/2=pi/z_i>=25pi/7.                         (11)
```

The cross expectation is finite before any coincident same-exterior
vacuum trace is assigned. The angle is monotone, so the unperturbed
coherent first-order Killing work satisfies

```
|W_coherent,U^(1)|/(nu N eta)
 <=[theta_f/(4eta)]1.05e-68 <4e-68.                 (12)
```

This uses the actual paired source, not a scalar electric pulse. It is
only the term evaluated in the unperturbed state. It does not bound the
second-order excitation energy or the finite drive's complete work.

On the identical full RN half-lines, freezing the domain by V(theta)
produces the bounded one-particle frame perturbation -i theta'_lab G/2.
Its integrated operator norm is theta_f/2=pi/1000. That small norm does
not bound the ultraviolet-summed energy-weighted Bogoliubov block.
Extending (12) to the actual finite drive requires a nonstationary
two-point/stress estimate with the TFD correlations and controlled
higher orders. Replacing the TFD by independent KMS reservoirs changes
the preparation. Source-kinetic contacts and their finite matched
coefficients are separate controller energy, not bounded by (10).

The [inherited higher-Landau result](higher-landau-boundary.md) has a
much sharper conditional pointwise exponential suppression on an
extremal exterior, but requires an unproved real-frequency inner load.
The present weaker bound removes that condition for this static thermal
RN reference and includes the full spectral tail. The global-AdS bulk
current theorem in [curved massive response](curved-massive-response.md)
has a different source, state and geometry and cannot be substituted.

An initially tempting extremal zero-energy moment must be rejected for
this thermal task. At extremality W>=M/(d-x_*), and a Neumann inverse-square
comparison has finite boundary H^-4 moment. At nonzero temperature
W decays exponentially toward the horizon. At omega=0, ingoing radial
transport gives r_boundary=i exp(-2S_n),
S_n=integral_(-infinity)^0 W_n dx_*<infinity, hence the spectral density
is sech(2S_n)>0. Its zero-energy inverse moment diverges even when this
number is tiny. The positive shift and explicit tail in (7) are essential.

Finally, the [S136 physical remainder map](physical-clock-remainder.md)
is exact about its stated classical reference. The separate
[S137 quantum centering](quantum-reference-clock-centering.md) supplies
the reference correction. This classical RN theorem needs a quantitative
operator/state/source transfer to that corrected exterior, then to the
evolving geometry and physical boundary matching.
It supplies no bound on S_T, beta, local contacts, bulk massive
polarization or the spatial area error by itself.

## Verification contract fixed before the checker

Verify the Clifford/domain identities and (2)-(4) symbolically; check
the rational candidate inequalities uniformly, with pi enclosed by
22/7 and elementary rational bounds. Verify the scalar Green formula,
its derivative matching and integral by independent ODE/spectral
representations before using its upper bound. Exact algebraic residuals
must vanish. Independent scalar comparisons use absolute/relative
tolerance 2e-8, refined to 2e-9; radial mass/potential comparisons
permit 2e-10 relative error. Numerical samples only supplement the
ordinary uniform proof, and do not replace any angular tail or state
assumption. Preserve any failed comparison, source hashes, precision,
command and terminal status. No historical checker or Lean campaign
is selected.
Before the first run select scalar cases (M,nu*d)=(18,.05),(32,.2),
(64,.001), plus the exact k=2000 candidate. Use spectral precision
80 digits and 110 digits on refinement, continuum integration through
P=max(500,4mu), and the analytic free-tail/difference remainder beyond P.
Backward ODE integration uses relative tolerances 2e-11 and 2e-12.
Geometry spot checks use k=2000,4000 at z_i=1/4,7/25 and synthetic
s_h=1/2000,1/1000 at X=1 and X=10^10. The uniform claims remain the
exact rational/symbolic proof, not these finite geometry cases.

The first frozen checker run passed: exact Clifford/domain/RN/Green
identities, 23 uniform rational inequalities, eight independent metric
cases, and all four scalar cases at both precisions. The largest scaled
ODE disagreement was below 1.6e-15. The largest raw continuum discrepancy
was below 4.7e-13 and inside its independently bounded omitted spectral
tail. The uniform full-tower comparison is the exact rational bound
in (10), not a numerical angular truncation. No checker failure occurred.
The receipt and raw log preserve the execution hashes and a separate
text-only refresh for this verification report, scope pointer and formatting.

Remaining to-do list: bound finite-drive energy/stress and local controller terms;
transfer to the quantum-corrected evolving exterior and common physical
state, then complete the S074-S079 gates.
