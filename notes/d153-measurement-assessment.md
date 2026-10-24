# D153 and measurements of the Dirac wormhole source

2026-10-23 Hawaii / 2026-10-24 UTC. Quantyra Space E003/S171.

**Decision: adapt the measurement discipline; park a direct optical test of
the current Dirac construction.** D153 supplies useful historical pointers
to electromagnetic correlation measurements. It supplies neither the
detector-to-Dirac-source relation nor a mechanism for preparing our paired
state. The bounded assessment is complete. The selected S169 construction
calculation remains the program priority.

## What the report establishes

D153 is *Quantum Tomography of Negative Energy States in the Vacuum*,
DIA-08-1102-007, dated 11 January 2011. Its 51 PDF pages, including 47
numbered pages, were read as extracted text from the exact private-S3
Release 06 bundle. Figures were not digitized or independently inspected;
no numerical result here is read from a figure. The original PDF and full
extracted text remain in S3. Exact keys, versions and SHA-256 hashes are
in [source custody](../evidence/d153-source-custody.json).

The important distinctions are already partly present in D153:

- Printed pp. 10-12, especially Eq. (9), distinguish positive mode energy
  from oscillating local energy density in a periodic-box squeezed state.
- Printed pp. 28-32 derive a quadrature measurement selected by the local
  oscillator, including a vacuum-loss model. The oscillator selects a
  spatial-temporal mode; it does not observe all field modes.
- Printed pp. 36-39, Eqs. (14)-(17), connect detector response to electric
  field two-point functions and a frequency/polarization-selected spectrum.
- Printed p. 43 expressly says the Hansen pulsed detector was not designed
  to measure individual pulse energy density. The proposed energy-density
  instrument is a development recommendation, not a demonstrated device.

Hansen et al.'s [original experiment](https://arxiv.org/abs/quant-ph/0104084)
demonstrates quadrature measurements and coherent-state tomography.
Marecki's [detector derivation](https://arxiv.org/abs/quant-ph/0703076)
and [Casimir calculation](https://arxiv.org/abs/0711.1541) make the measured
field, smearing, frequency selection and background explicit. The latter's
Eq. (21) relates output variance to a weighted electric spectral density;
it is not a general reconstruction of every stress-tensor component.

The broad language on energy conditions and alleged quantum-inequality
violations in D153 pp. 4-5 and 42-44 cannot serve as a premise of our
construction. In particular, a stationary cavity's reduced electric
spectrum relative to a no-plate reference does not, by itself, contradict
a free-space time-averaged stress-energy inequality. Boundary conditions,
field algebra, reference state and observable must match the theorem.
Marecki explicitly notes the background-dependent commutator in the
Casimir paper, p. 6. We have not established that every inequality claimed
in this literature is correct; we have established that D153's proposed
comparison does not supply a matched theorem violation.

## A useful but restricted optical witness

Use hbar=c=1. For one travelling electromagnetic mode in a periodic
normalization volume V, with all other modes in their vacuum, define

    X_theta=(a exp(-i theta)+a_dagger exp(i theta))/sqrt(2),
    n=<a_dagger a>, m=<aa>, v_theta=<X_theta^2>.

The raw second moments obey

    v_theta=n+1/2+Re(m exp(-2i theta)),
    n=(v_0+v_pi/2-1)/2,
    Re(m)=(v_0-v_pi/2)/2,
    Im(m)=v_pi/4-(v_0+v_pi/2)/2.                     (1)

Thus three phase settings determine this quadratic target. Full Wigner
tomography is unnecessary. This is standard moment algebra, not a new
tomography method. If experimental variances replace raw second moments,
include the measured means or justify the zero-mean assumption.

Choose the field phase convention so the vacuum-normal-ordered density is

    rho(t,x)=K[n+Re(m exp(2i(k.x-omega t)))], K=omega/V>0.

Here the electric and magnetic contributions are equal for this one
travelling mode. That equality is an assumption of this example, not a
property of arbitrary cavity modes or fields. For a normalized positive
time weight w(t), at x=0,

    rho_bar=K[n+Re(m W)], W=integral w(t)exp(-2i omega t)dt.

For a Gaussian of standard deviation tau centered at zero, W is real:

    rho_bar/K=n+exp(-2 omega^2 tau^2)Re(m).           (2)

Pure squeezed vacuum has n=sinh(r)^2 and
m=-exp(i phi)sinh(r)cosh(r). At r=0.7, phi=0 the checked values are

| Quantity | Value |
| --- | ---: |
| Minimum quadrature second moment; vacuum is 0.5 | 0.123298481971 |
| Instantaneous rho/K at the squeezed phase | -0.376701518029 |
| Gaussian-sampled rho/K at omega*tau=1 | +0.446589641163 |
| Positive total excitation energy divided by omega | +0.575449232697 |

Squeezing alone therefore does not fix the sign of a specified sampled
energy. A squeezed Gaussian state also has an everywhere nonnegative
Wigner function. Conversely, the one-photon state has Wigner negativity
but n=1,m=0 and positive density in this model. These are different tests.
The example neither violates a QEI nor bounds a curved Dirac source.

## Detector uncertainty and an honest negative certificate

For the same calibrated mode, scalar vacuum loss gives

    v_observed=eta*v+(1-eta)/2, n_observed=eta*n,
    m_observed=eta*m, F_observed=eta*F,
    F=n+Re(m W), 0<eta<=1.                          (3)

Pure loss alone preserves this witness's sign; uncertain efficiency
affects magnitude. Do not turn it into a spurious sign obstruction.
Unknown additive electronic noise, reference drift, gain, phase error,
unequal phase efficiencies and mode mismatch need their own model.

For independent Gaussian preparations at each phase, the usual chi-square
variance interval can bound v_theta under the zero-mean model. Allocate
error probabilities alpha/3 to the three phase intervals; if calibration
is also learned, include its failure probability in the total budget.
For W=A+iB, propagate the joint intervals through

    F=((1+A+B)/2)v_0+((1-A+B)/2)v_pi/2-B*v_pi/4-1/2.

Use interval extrema and the allowed calibration/sampler parameters.
This is established Gaussian inference; no new coverage advantage is
claimed or tested. It does not apply unchanged to correlated continuous
wave samples or arbitrary non-Gaussian preparations with unbounded tails.

Even perfect measurement of F is only a mode contribution. For a full
electromagnetic target E=K*F+R, where R includes unobserved field content,
a negative certificate requires an upper bound E_upper<0, including R.
No bound on R is supplied by measuring one optical mode. Converting this
to absolute units also requires the physical mode normalization K.
The Casimir problem additionally requires magnetic correlations, field
components, material/boundary effects and a consistent subtraction.

## Why this does not measure our Dirac source

The retained S169 calculation uses the charged Dirac field, full angular
and thermal spectrum, the specified paired boundary/state, and common
renormalization/restoration. Its finite radial expression is

    D=[[0,partial_rho+W],[partial_rho-W,0]],
    m=M_n exp(b), Gamma_M=-Tr log(D+m),
    M_n^2=n(n+q)/C, d_n=2(q+2n).

Gamma_M is only the regulated mass-vertex part. The Weyl allocation and
common four-dimensional restoration are essential; its raw matrix sign
does not determine a physical energy sign. See the actual
[coincident action](current-coincident-dirac-action-assembly.md),
[subtracted thermal action](current-subtracted-thermal-action.md), and
[physical dual conversion](current-physical-york-dual-assembly.md).

A candidate future observable must specify a Lorentzian state, an observer
or null vector, a smooth sampling weight and an energy/tensor component:

    E_psi[w,u]=integral w(t)<T_mu_nu^psi>_ren u^mu u^nu dt,
    N_psi[f,k]=integral f(lambda)<T_mu_nu^psi>_ren k^mu k^nu dlambda.

Neither equals the Euclidean weak action correction or its Hessian. A
Lorentzian realization and continuation/response argument are additional
requirements. Point-splitting the Dirac source schematically requires
D_mu_nu^psi(S_psi-H_psi)+local conserved terms, with a spinor two-point
function and the model's fixed finite counterterms. Electric correlations
do not determine that spinor function. A bosonic CCR covariance is not a
fermionic CAR covariance, and optical squeezing is not preparation of our
magnetically supported paired Dirac state.

The immediate S169 obstacle also concerns the metric derivative of the
source: the restored quantum Hessian and its boundary/contact terms.
Knowing a mean stress at one background would not determine that
derivative or its retarded realization. D153 supplies no replacement for
the current physical operator assembly or its controlled remainder.

There is a simple identifiability obstruction in a precisely limited
observation class. On a fixed background with independent field factors,
let the detector POVM be M_y tensor I_psi. For states
rho_EM tensor sigma_psi and rho_EM tensor sigma'_psi,

    p(y)=Tr(rho_EM M_y)

is identical, while Tr(sigma_psi A_psi) and
Tr(sigma'_psi A_psi) may differ. Therefore no estimator based on these
optical data universally recovers an independent Dirac target. This is a
standard tensor-product argument, not a new impossibility theorem about
interacting QED or self-consistent gravity. With a theory-fixed state the
target may instead be calculated; an interacting probe could add data.
Both require a specified physical connection absent from D153 and the
current observation model. For example, in flat space squeezed photons
can coexist with an independent free Dirac vacuum whose normal-ordered
stress is zero. Optical squeezing therefore cannot universally certify
a negative Dirac source. This example is not a competing self-consistent
state of the S169 geometry.

## Energy inequalities and subsequent measurements

[Fewster and Mistry](https://arxiv.org/abs/gr-qc/0307098), Eq. (3.2), give
the established flat-space free-Dirac comparison, with their conventions,

    integral g(t)^2<T_00>_normal dt >=
        -1/(12 pi^2) integral |g''(t)|^2 dt.

Their mass-dependent Eq. (2.29) is tighter. This timelike result is not
a null-energy bound or an unchanged bound for our charged, curved,
boundary-coupled problem. [Fewster and Verch](https://arxiv.org/abs/math-ph/0105027),
Theorem 4.1, establishes a reference-Hadamard-state difference QWEI on
globally hyperbolic Dirac backgrounds under its hypotheses. Applying an
appropriate result to our geometry still requires checking the domain,
state, external field, boundary conditions and sampler.

[Maclay and Davis](https://arxiv.org/abs/1806.01269) later claim an optical
QI conflict, but their section V explicitly requires identifying a
local-oscillator phase scan with a physical time sampler and discusses
the carrier/sideband discrepancy. Their claim is not accepted here as a
failure of Dirac or electromagnetic QEIs. Its alleged violation would
need a common observable, filter and sampling protocol first.

Later experiments do advance measurement access beyond a 2011 survey:
[Riek et al.](https://arxiv.org/abs/1611.06773) resolve subcycle electric
noise below a bare-vacuum reference, while
[Benea-Chelmus et al.](https://arxiv.org/abs/1809.01785) measure spatial and
temporal electric correlations within an explicitly modeled electro-optic
bandwidth. Neither paper measures our Dirac stress or a gravitational
response to negative energy. The October 2026
[Fewster-Ford preprint](https://arxiv.org/abs/2610.14927) improves explicit
stationary-worldline scalar/electromagnetic QEI bounds; its hypotheses
also differ from S169. The scoped search found useful established
comparators, not a cleared novelty claim for our proposed measurement.

## Verification and disposition

The [prespecified protocol](d153-measurement-protocol.md) preceded the
[checker](../checks/d153_measurement_assessment.py). Direct 512-level
Fock-vector quadrature operations, second-moment reconstruction, independent
Gaussian time integration and scalar-loss propagation passed all 108
scoped parameter cases. The maximum identity residual is 6.32e-15 against
the fixed 1e-10 threshold. The [receipt](../evidence/d153-measurement-assessment.json)
pins the inputs and controls. These are floating-point development checks,
not certified tails, experimental observations or formal acceptance.

The [source review](../evidence/d153-source-review.json) distinguishes full
report text, selected primary derivations and later measurement/QEI updates.
Source PDFs and full extracted texts are private S3 objects; authored
notes, code and hash/version receipts are local research artifacts.

**Use:** cite the primary detector work and retain its mode, bandwidth,
reference and sampler requirements in future observable design.
**Adapt:** if an independently justified electromagnetic diagnostic is
later selected, begin with a specified quadratic witness and bounds on
unobserved contributions, rather than full tomography by default.
**Park:** equipment design or an optical validation claim for the current
Dirac wormhole. Reopen only with a specified detector coupling or inverse
map, the relevant renormalized target and a quantitative uncertainty bound.
No such follow-up is launched by this assessment.

Remaining to-do list: none for S171's bounded assessment. The S169
physical correction, validity and stability obligations remain open.
