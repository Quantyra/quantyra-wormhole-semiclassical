# Causal shape-pole diagnostic and the order of work

2026-10-22, E003/S169. Ordinary argument before frozen verification.
The full stable semiclassical wormhole objective remains incomplete.

The retained frozen-coefficient quantum shape model has a pair of
upper-half-frequency-plane poles under the causal continuation. Their
proper frequencies have 0.80<abs(Re omega)sqrt(G)<0.82 and
0.012<Im omega sqrt(G)<0.014. This is an EXECUTED MODEL DIAGNOSTIC,
not a pole theorem for the actual finite geometry or a physical
Planck-scale prediction. It changes the immediate priority: resolve
the causal realization/validity question before optimizing the remaining
stationary shape source. No action coefficient or boundary condition
is changed, and no stable solution is claimed.

## 1. Retain the exact benchmark, finite contact and model scope

Use the explicitly labeled frozen-frequency model, equations (4)-(6) of
physical-metric-block-and-preconditioner.md. Per physical Einstein volume,
the Euclidean shape entry is

    F_* k^2/(12piG)+k^4[log(k/mu_E)-6/5]/(120pi^2),
    F_*=1-384pi(c_R-1/(5760pi^2))g, g=G lambda,
    mu_E=mu0 exp(-u_b), |u_b|<5.4e-37.

Keep c_R=1e22, q=2000^4, g_R=2000^-4, mu0=1/b_classical,
the complete Dirac tower and the original finite contact. The scalar
curvature mode has positive mass; its introduction did not cancel the
independent shape logarithm. At frozen flat principal order a genuine
TT variation has delta R=0, so a local R^2 square supplies no new
fourth-order shape coefficient. The factor F_* above remains.

Set k0=2000 (a benchmark integer, not frequency), and

    g0=k0^-10, g+11g^2/(360pi)=g0,
    x0=(sqrt(1+12pi k0^6)-1)/6,
    z=G mu_E^2=(g0/x0)exp(-2u_b),
    a=20pi F_*, zeta=z exp(12/5).

The normalized shape factor is exactly

    d(t)=1+(t/a)Log(t/zeta), t=G k_E^2.                 (1)

The term exp(12/5) records the existing finite contact; it is not a
new scale choice. For t>0 its minimum is 1-zeta/(a e)>1-7.70e-45.
This agrees with the earlier Euclidean test. Magnetic degeneracy q
is not inserted again as a number of four-dimensional matter species:
the logarithmic coefficient already came from the whole angular tower.

Equation (1) retains the leading flat symbol and its stated finite
contact. Curvature, state, finite boundary, mixing and lower-order
terms were explicitly NOT computed by that model. They are not zero
here. In particular a local symbol is not an actual global retarded
characteristic determinant.

## 2. Causal sheet and exact growing-pole construction

Use the convention exp(-i omega time), Omega=omega sqrt(G), and
pbar^2=G|p|^2. Causal continuation uses

    t=pbar^2-Omega^2, Im Omega>0,
    D(Omega)=1-Omega^2/a Log(-Omega^2/zeta)              (2)

at pbar=0. Log is the principal logarithm in t, continued from positive
Euclidean t. For positive real omega its boundary value is
log(omega^2 G/zeta)-i pi. Choosing the opposite sign would change
the response prescription. A growing inverse need not have an ordinary
real-axis Fourier transform: its causal inverse uses a Laplace contour
above the poles. One cannot delete upper poles while retaining the
same causal initial-value problem.

For a branch-independent computation put A=a/zeta. Define theta in
(pi/2,pi) as the unique root of

    f(theta)=-theta cot(theta)+log(theta/sin(theta))=log A. (3)

On this interval

    f'=theta csc(theta)^2+1/theta-2cot(theta)>0,

and f ranges from log(pi/2) to infinity. The benchmark log A is
well above log(pi/2), so the root exists uniquely. Let

    U=-theta cot(theta)>0, W=U-i theta,
    Omega_+=sqrt(a/W), Re Omega_+>0, Im Omega_+>0,
    Omega_-=-conj(Omega_+).                             (4)

Indeed W exp(W)=-A and, on the CAUSAL sheet,
Log[-a/(zeta W)]=W. Thus D(Omega_+)=0 exactly and the same holds
for Omega_-. The Lambert W_{-1}(-A) expression is only a secondary
numerical comparator. The monotone real equation fixes the correct
sheet without relying on a library's branch convention.

For interval evaluation let M=sqrt(U^2+theta^2). Then

    Re Omega_+=sqrt[a(M+U)/(2M^2)],
    Im Omega_+=theta sqrt[a/(2M^2(M+U))].               (5)

These positive real formulas enclose both parts without ambiguous
complex square roots. The whole |u_b| interval is included.
The frozen run uses interval endpoint signs in (3), followed by (5),
and a separate high-precision Lambert/residual comparison.
It certifies the bounds stated at the start. In the model the amplitude
therefore grows as exp(Im Omega time/sqrt(G)), with e-fold time

    71 < time_e/sqrt(G) <84.                            (6)

These are normalized proper-frequency statements at the frozen point,
not an asymptotic-clock eigenfrequency of the finite wormhole.

## 3. A numerical transfer test, with the hypotheses exposed

The root is simple. At it L0=Log(-Omega_+^2/zeta)=W and

    D'=-2Omega(L+1)/a,
    D'(Omega_+)=-2(1+1/W)/Omega_+,
    D''=-2(L+3)/a.

Use the certified coarse bounds .80<|Omega_+|<.83,
90<|W|<105 and a>62. On the complex disk
|Omega-Omega_+|<=r0=.005, the imaginary part remains positive.
Its logarithm is analytic and

    |L-L0|<=2r0/(.80-r0), |D'(Omega_+)|>2.38,
    sup |D''|<3.5.

Taylor's integral remainder on the boundary gives the explicit margin

    |D|>2.38(.005)-(3.5/2)(.005)^2
        =.01185625>.0118.                              (7)

Comparison with D'(Omega_+)(Omega-Omega_+) also proves exactly one
model zero inside this disk. The conjugate disk has the same margin.

Suppose, and only suppose, the actual physical causal problem admits
an analytic reduction to a characteristic factor

    D_phys(Omega)=D(Omega)+Sigma(Omega)

on this disk, after division by nonvanishing factors and inversion of
the other constrained blocks. Then the sufficient numerical condition

    sup_boundary |Sigma|<.01                           (8)

forces an upper-half-plane zero by Rouche's theorem. The condition must
include curvature, all quantum state/contact terms, physical initial
constraints and the finite boundary realization. A zero of that factor
must reconstruct a non-gauge physical mode. An arbitrary local symbol
error, a stationary Euclidean form bound, or small background curvature
does NOT establish these hypotheses or (8).

Neither this analytic physical reduction nor a bound on Sigma has been
proved. Sigma is an UNKNOWN remainder, not a known upper bound that
fails .01. A failure to bound it would not establish stability either.
If (8) is proved for the actual realization, it would rule out an
assertion of unrestricted causal stability for the retained equations.
The required physical stability theorem is instead tied to an established
semiclassical validity band; no new demand for physical Planck-scale
stability is imposed. The present result by itself proves neither a
finite-geometry obstruction nor failure inside that physical band.

## 4. Validity is a separate quantitative requirement

The pole lies near the already computed Euclidean crossover
k sqrt(G)=.8089653..., far above the scalar mass
m sqrt(G) about 5.76e-13. Semiclassical validity has not been established
at that pole. This neither removes the mathematical pole from (1) nor
turns it into a reliable prediction of a physical runaway.

For scale comparison only, the SAME model admits a whole complex-disk
bound on the principal sheet:

    |t|<=1e-6 implies |d(t)-1|<1.4e-6.                 (9)

To see it, |arg t|<pi and, writing s=|t|,
s|log(s/zeta)| is bounded by zeta/e below zeta and by
1e-6 log(1e-6/zeta) above it. Add pi*1e-6 and divide by a.
Thus there is no normalized shape-factor zero in this model band.
The massless k^2 factor is kept; (9) does not discard the graviton.
This calculation is not a selected cutoff, an EFT error bound, a
finite-geometry low-frequency stability theorem or a change of equations.

The distinction is established prior art. Anderson, Molina-Paris and
Mottola, [Appendix A (A22)-(A23) and sections IV-V](https://arxiv.org/pdf/gr-qc/0209075),
derive the retarded continuation and distinguish infrared stability in
their flat scalar-field example from Planck-scale semiclassical failure.
Hu, Roura and Verdaguer, [Appendix E.3 and Figure 2](https://arxiv.org/pdf/gr-qc/0402029),
explain the different causal growing and final-condition non-growing
inverses. The latter prescription can be acausal. Neither paper proves
the present finite charged Dirac-domain spectrum. The logarithmic
Lambert-pole calculation also has the comparator
[Calmet--Kuntz--Mohapatra, (1)-(4)](https://link.springer.com/article/10.1140/epjc/s10052-016-4265-8).
We do not import its field-count convention, discard a pole, or change
mu0 as suggested in a different model discussion.

An EFT or order-reduction route would require an explicit action/domain
decision and controlled omission, state and noise errors. It cannot be
silently substituted for the selected exact retained equations. In
particular the deliberately retained light curvature scalar cannot be
removed by blanket order reduction merely to obtain an easier stable
model. Changing renormalization scale alone, with compensating running
couplings, leaves the physical logarithmic combination unchanged.

## 5. Research decision and unchanged solution requirements

Prioritize the sharper CAUSAL REALIZATION/POLE-TRANSFER calculation (8)
before further stationary shape-source optimization. The concrete reason
is the executed upper-half-plane pole test of the very shape model used
to motivate retaining the full higher-order inverse. It exposes an
unresolved earlier stability/validity issue, not a theorem that this
construction or semiclassical wormholes are impossible.

The single next calculation is the physical constrained retarded
characteristic reduction and its complex-disk remainder in (8), with
the actual paired boundary/state prescription. This checks the causal
equations being used; it does not demand trust in their Planck-scale
physical predictions. A physical validity prescription must separately
control the effect of excluded scales within its stated band. The model pole test and
its transfer tolerance are executed; the actual reduction and transfer
are still UNEXECUTED. Do not infer physical instability from the model,
or proceed as if Euclidean positivity settled the causal question.

The construction, all fixed coefficients, stationary inverse and trace
defect result remain preserved. The full shape-defect/adjoint-bolt
pairing remains uncomputed. If the causal/validity gate permits continued
use of the construction, its cold difference <1.667e-9, total Y0<=5e-9,
Z0<=.25, Z2<=5e7 at rho=1e-8, all admission margins, matching and
nonlinear requirements still apply. No S169 checkbox is completed and
the full stable-wormhole objective is not narrowed to stationary existence
or model-band stability.

Remaining to-do list: physical causal realization and pole-transfer test;
full shape-defect/adjoint-bolt pairing and weighted error; nonlinear
solution, matching, validity, all-sector stability and the remaining
opening, transport, preparation and separate creation requirements.
