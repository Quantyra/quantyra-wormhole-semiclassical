# D153 measurement assessment protocol

2026-10-23 Hawaii / 2026-10-24 UTC. Quantyra Space E003/S171.
Ordinary argument and numerical scope fixed before checker execution.

The question is whether D153's electromagnetic homodyne observables can
constrain the renormalized Dirac source in the retained S169 model. Read
the report, the detector derivation, relevant energy-inequality proofs,
and selected later measurement work. A use, adapt or park decision is
the endpoint; this does not select an experiment or change the action.

## Ordinary argument

For a zero-mean single electromagnetic mode set hbar=c=1,
X_theta=(a exp(-i theta)+a_dagger exp(i theta))/sqrt(2),
n=<a_dagger a>, m=<aa>. Then

    v_theta=<X_theta^2>=n+1/2+Re(m exp(-2i theta)),
    n=(v_0+v_pi/2-1)/2,
    Re(m)=(v_0-v_pi/2)/2,
    Im(m)=v_pi/4-n-1/2.

These moments do not require a Gaussian state, though an exact chi-square
variance interval would require independent Gaussian observations. A
zero mean must be established or first moments included explicitly.
For one travelling plane-wave mode in a periodic normalization volume V,
all other modes vacuum, choose the field phase so its vacuum-normal-ordered
local energy density is

    rho(t,x)=(omega/V)[n+Re(m exp(2i(k.x-omega t)))].

This is a restricted optical model, not a curved-spacetime renormalization
prescription. A positive normalized Gaussian time sampler of width tau,
centered at zero at x=0, gives

    rho_bar/(omega/V)=n+exp(-2(omega tau)^2) Re(m).

For a pure squeezed vacuum n=sinh(r)^2 and
m=-exp(i phi)sinh(r)cosh(r). Local negative intervals can coexist with a
positive long-sampled energy and positive total excitation energy.
Vacuum-referenced scalar pure loss sends n,m to eta*n,eta*m: the witness
is multiplied by eta, so pure loss alone cannot reverse its sign. Unknown
additive calibration offsets, mode mismatch and missing modes require
separate bounds. Do not label all inefficiency as a sign ambiguity.

For a fixed-background uncoupled product experiment with a photon-only
POVM M_y tensor I_psi, two states rho_EM tensor sigma_psi and
rho_EM tensor sigma'_psi have identical optical statistics even when a
Dirac target has different expectations. This rules out universal recovery
of that target from those data in this observation class. It does not
exclude a specified interacting detector, a fixed state's theoretical
calculation, or a gravitationally mediated inverse problem. Neither D153
nor our existing detector model supplies those additional connections.

The S169 field operator, full angular/thermal tower, paired state,
Hadamard subtraction, local action terms and boundary contributions remain
the physical inputs. An electromagnetic covariance cannot replace them.

## Prespecified numerical checks

Use direct photon-number state vectors (512 levels) for squeezed states
r in {0.1,0.7,1.5}, phase in {0,0.37,1.1}; independently compute the
quadrature operator norms at 0, pi/4, pi/2. Compare the reconstructed
moments with the analytic values. Integrate the local state-vector energy
against Gaussian samplers with omega*tau in {0,0.1,1,3}, using composite
Simpson integration with 12,000 intervals on [-10tau,10tau]. Check losses
eta in {0.4,0.8,1}. All identity residuals must be below 1e-10. Report
the Fock tail/norm diagnostic, numeric residuals, parameters and hashes.

Retain explicit controls: r=0.7 squeezed covariance with negative
instantaneous witness but positive omega*tau=1 witness; the one-photon
state with negative Wigner function at the origin but positive energy;
and two Dirac-proxy states with exactly identical photon POVM probabilities
but distinct target values. The last is a finite-dimensional illustration
of the ordinary tensor-product proof, not a physical Dirac simulation.

No experimental data, fitting, confidence-coverage claim, new theorem
certificate, physical Dirac source measurement, or wormhole validation is
produced by this check. If the field/observable connection fails, park the
direct application and retain only a justified diagnostics reference.
