"""Exact algebra supporting the stationary scalar range characterization.

No coupled inverse or actual zero moment is computed. Spectrum, endpoint
and range estimates remain the reviewed ordinary arguments in the note.
"""
from __future__ import annotations

import base64
import contextlib
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import traceback

import sympy as S

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "evidence/stationary-trace-range-condition.json"
LOG = ROOT / "evidence/stationary-trace-range-condition.log"
INPUTS = [
    "checks/stationary_trace_range_condition.py",
    "notes/stationary-trace-range-condition.md",
    "notes/stationary-exterior-solver-obstruction.md",
    "notes/spherical-anomaly-action.md",
    "notes/spherical-anomaly-sources.md",
    "notes/stationary-centered-covariant-seed.md",
    "notes/stationary-centered-complete-residual.md",
    "notes/stationary-completed-metric-response-norm.md",
    "notes/stationary-centered-first-order-geometry.md",
    "notes/uniform-dirac-reference.md",
]
CONTRACT = {
    "acceptance": "all displayed symbolic residuals exactly zero and rational inequalities true",
    "numerical_tolerances": "none",
    "hash_policy": "SHA-256 of exact raw input bytes; strict UTF-8, LF, no BOM",
    "scope": "exact identities supporting the fixed-metric scalar range condition and fast inverse",
    "ordinary_arguments": ["regular-bolt existence", "asymptotic solution coefficients",
                           "essential-spectrum Weyl sequence", "compatible-source range and uniqueness"],
    "excluded": ["formal certification of the ordinary analytic proof", "coupled metric inverse",
                 "actual zero compatibility moment", "existence or error bar for a full stationary solution",
                 "boundary Hamiltonian", "physical state matching or gate closure"],
    "failure_policy": "append every run/log and retain failed-run input bytes",
}
ROWS = []
R = S.Rational


def check(label, condition, detail=""):
    row = {"label": label, "passed": bool(condition), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row, sort_keys=True), flush=True)


def zero(label, expression):
    value = S.factor(S.together(expression))
    check(label, value == 0, value)


def trace_identities():
    G, le, Z, boxZ, K, Ric2, Fem2 = S.symbols("G le Z boxZ K Ric2 Fem2")
    a, b, c = R(11, 360), R(1, 30), R(1, 20)
    Rstar = -12*le
    eps = G/(60*S.pi)
    ainf = 1-G*a*Rstar/(6*S.pi)
    lam0 = le+11*G*le**2/(360*S.pi)
    r4 = Rstar+Z
    C2 = K-2*Ric2+r4**2/3
    E4 = K-4*Ric2+r4**2
    J2 = Ric2-r4**2/4
    zero("invariant_Euler_decomposition", E4-(C2-2*J2+r4**2/6))
    zero("vacuum_root_trace", Rstar+12*lam0-G*a*Rstar**2/(12*S.pi))
    original = r4+12*lam0+G/(2*S.pi)*(c*C2-a*E4+b*boxZ-R(2, 3)*Fem2)
    source = G/(2*S.pi)*(R(2, 3)*Fem2-R(7, 360)*C2-R(11, 180)*J2+R(11, 2160)*Z**2)
    zero("full_nonlinear_trace_reduction", original-(ainf*Z+eps*boxZ-source))
    zero("retained_scalar_epsilon", G*b/(2*S.pi)-eps)
    zero("retained_scalar_a_inf", ainf-(1+11*G*le/(180*S.pi)))
    zero("proper_mass_squared", ainf/eps-(60*S.pi/G+R(11, 3)*le))
    zero("trace_source_divided_by_epsilon", source/eps
         -(20*Fem2-R(7, 12)*C2-R(11, 6)*J2+R(11, 72)*Z**2))
    q, r = S.symbols("q r", positive=True)
    zero("magnetic_leading_source", (G/(2*S.pi)*R(2, 3)*Fem2).subs(Fem2, q**2/(2*r**4))
         -G*q**2/(6*S.pi*r**4))
    scale = S.symbols("scale", positive=True)
    cc, jj, ff, zz = S.symbols("cc jj ff zz")
    abstract_source = G/(2*S.pi)*(R(2, 3)*ff-R(7, 360)*cc-R(11, 180)*jj+R(11, 2160)*zz**2)
    transformed = abstract_source.subs({G: scale**2*G, ff: ff/scale**4,
                                       cc: cc/scale**4, jj: jj/scale**4,
                                       zz: zz/scale**2}, simultaneous=True)
    zero("trace_source_length_dimension", transformed-abstract_source/scale**2)


def static_and_seed_identities():
    r, ell = S.symbols("r ell", positive=True)
    F, al, Z0, dF, dal, dZ = [S.Function(n)(r) for n in ("F", "al", "Z0", "dF", "dal", "dZ")]

    def box(fr, ar, value):
        weight = S.exp(ar)*r**2
        return S.diff(weight*fr*S.diff(value, r), r)/weight

    zero("static_density_operator", box(F, al, Z0)
         -(F*S.diff(Z0, r, 2)+(S.diff(F, r)+F*S.diff(al, r)+2*F/r)*S.diff(Z0, r)))
    varied = S.diff(box(F+ell*dF, al+ell*dal, Z0+ell*dZ), ell).subs(ell, 0)
    expected_delta = dF*S.diff(Z0, r, 2)+(S.diff(dF, r)+dF*S.diff(al, r)
                       +F*S.diff(dal, r)+2*dF/r)*S.diff(Z0, r)
    zero("delta_Box_including_seed_term", varied-box(F, al, dZ)-expected_delta)

    t, x, m, d, B = S.symbols("t x m d B", positive=True)
    seedF = 1+x*t**2-2*m/t+(1+3*x+d)/t**2-2*d/t**3+d/t**4
    scalar = -S.diff(seedF, t, 2)-4*S.diff(seedF, t)/t+2*(1-seedF)/t**2
    zero("centered_seed_curvature", scalar+12*x-(4*d/t**5-6*d/t**6))
    zr = 4*d*B**3/r**5-6*d*B**4/r**6
    witness = S.diff(zr, r, 2)+2*S.diff(zr, r)/r
    zero("nonconstant_seed_delta_Box_witness", witness-20*d*B**3*(4*r-9*B)/r**8)
    check("seed_delta_Box_not_identically_zero", S.factor(witness) != 0)

    eps, ainf, fsource, zh, fhp, wh, p1 = S.symbols("eps ainf fsource zh fhp wh p1", nonzero=True)
    zhprime = (fsource-ainf*zh)/(eps*fhp)
    zero("regular_horizon_inhomogeneous_relation", ainf*zh+eps*fhp*zhprime-fsource)
    zero("regular_horizon_homogeneous_derivative", zhprime.subs({fsource: 0, zh: 1})+ainf/(eps*fhp))
    y, rho = S.symbols("y rho", positive=True)
    horizon_phi = 1-ainf*y/(eps*fhp)
    zero("regular_horizon_proper_expansion", horizon_phi.subs(y, fhp*rho**2/4)
         -(1-ainf*rho**2/(4*eps)))
    zero("horizon_log_Wronskian_normalization", p1*y*S.diff(S.log(y)/p1, y)-1)


def liouville_and_channel():
    rho = S.symbols("rho", positive=True)
    mu = S.Function("mu")(rho)
    psi = S.Function("psi")(rho)
    scalar = psi/S.sqrt(mu)
    op = -S.diff(mu*S.diff(scalar, rho), rho)/mu
    potential = S.diff(S.log(mu), rho, 2)/2+S.diff(S.log(mu), rho)**2/4
    zero("Liouville_conjugation", S.sqrt(mu)*op+S.diff(psi, rho, 2)-potential*psi)
    horizonV = potential.subs(mu, rho).doit()
    zero("critical_bolt_potential", horizonV+1/(4*rho**2))
    for name, value in (("regular", S.sqrt(rho)), ("logarithmic", S.sqrt(rho)*S.log(rho))):
        zero("critical_bolt_"+name+"_solution", -S.diff(value, rho, 2)-value/(4*rho**2))
    check("both_bolt_branches_are_L2", S.integrate(rho, (rho, 0, 1)) == R(1, 2)
          and S.integrate(rho*S.log(rho)**2, (rho, 0, 1)) == R(1, 4))
    zero("log_branch_form_energy_density", rho*S.diff(S.log(rho), rho)**2-1/rho)
    check("log_branch_form_energy_diverges", S.limit(-S.log(rho), rho, 0, dir="+") == S.oo)

    r, le, G, B = S.symbols("r le G B", positive=True)
    Dproper = lambda value: S.sqrt(le)*r*S.diff(value, r)
    mu_outer = S.sqrt(le)*r**3
    Vouter = Dproper(Dproper(S.log(mu_outer)))/2+Dproper(S.log(mu_outer))**2/4
    zero("outer_Liouville_potential", Vouter-9*le/4)
    eps = G/(60*S.pi)
    ainf = 1+11*G*le/(180*S.pi)
    k2 = ainf/eps-Vouter
    zero("proper_channel_momentum", k2-(60*S.pi/G+R(17, 12)*le))
    omega2 = 60*S.pi/(G*le)+R(17, 12)
    zero("logarithmic_channel_frequency", k2/le-omega2)
    omega = S.symbols("omega", real=True)
    for sign, name in ((1, "plus"), (-1, "minus")):
        delta = R(3, 2)+sign*S.I*omega
        zero("outer_indicial_"+name,
             S.expand(ainf+eps*le*delta*(delta-3)).subs(omega**2, omega2))
    fr = r**-4
    box_leading = le*r**2*S.diff(fr, r, 2)+4*le*r*S.diff(fr, r)
    zero("fast_particular_curvature_multiplier", box_leading-4*le*r**-4)
    q = S.symbols("q")
    particular = G*q**2/(6*S.pi*(ainf+4*eps*le))*r**-4
    zero("fast_particular_curvature_source", (ainf+4*eps*le)*particular-G*q**2/(6*S.pi*r**4))

    T, shift, k, DV = S.symbols("T shift k DV", positive=True)
    chi = S.Function("chi")
    env = chi((rho-shift)/T)
    wave = S.exp(S.I*k*rho)*env/S.sqrt(T)
    eps1, Vinf = S.symbols("eps1 Vinf")
    Pwave = eps1*(S.diff(wave, rho, 2)+(k**2-DV)*wave)
    expected = eps1*S.exp(S.I*k*rho)/S.sqrt(T)*(S.diff(env, rho, 2)
                        +2*S.I*k*S.diff(env, rho)-DV*env)
    zero("Weyl_cutoff_residual", Pwave-expected)
    check("Weyl_cutoff_envelope_orders", S.diff(env, rho).has(1/T)
          and S.diff(env, rho, 2).has(T**-2))


def green_identities():
    r = S.symbols("r", positive=True)
    eps, ainf = S.symbols("eps ainf", nonzero=True)
    p, w, phi, psi, z, f, A, Q = [S.Function(n)(r) for n in ("p", "w", "phi", "psi", "z", "f", "A", "Q")]
    P = lambda value: ainf*value+eps*S.diff(p*S.diff(value, r), r)/w
    W = phi*S.diff(psi, r)-S.diff(phi, r)*psi
    zero("Green_Lagrange_identity", w*(phi*P(z)-z*P(phi))
         -eps*S.diff(p*(phi*S.diff(z, r)-S.diff(phi, r)*z), r))
    homogeneous = {
        S.diff(phi, r, 2): -S.diff(p, r)*S.diff(phi, r)/p-ainf*w*phi/(eps*p),
        S.diff(psi, r, 2): -S.diff(p, r)*S.diff(psi, r)/p-ainf*w*psi/(eps*p),
    }
    zero("homogeneous_Wronskian_conservation", S.diff(p*W, r).subs(homogeneous))
    upper = {
        S.diff(A, r, 2): -S.diff(psi*w*f, r),
        S.diff(Q, r, 2): -S.diff(phi*w*f, r),
        S.diff(A, r): -psi*w*f,
        S.diff(Q, r): -phi*w*f,
    }
    fast = (phi*A-psi*Q)/eps
    zero("fast_Green_first_derivative_endpoint_cancellation", S.diff(fast, r).subs(upper)
         -(S.diff(phi, r)*A-S.diff(psi, r)*Q)/eps)
    fast_residual = (P(fast)-f).subs(upper).subs(homogeneous)
    zero("fast_Green_source_sign", fast_residual-(p*W-1)*f)

    y, p1, moment, constA = S.symbols("y p1 moment constA", positive=True)
    horizon_expression = (constA-moment*S.log(y)/p1)/eps
    zero("fast_Green_horizon_log_coefficient", S.diff(horizon_expression, y)*y
         +moment/(eps*p1))
    trace_defect = P(z)-f
    zero("candidate_defect_projection_integrand", w*phi*trace_defect+w*phi*f
         -eps*S.diff(p*(phi*S.diff(z, r)-S.diff(phi, r)*z), r)-w*z*P(phi))

    B, le, omega = S.symbols("B le omega", positive=True)
    ap, bp, az, bz = S.symbols("ap bp az bz")
    angle = omega*S.log(r/B)
    weight = (r/B)**(-R(3, 2))
    ph = weight*(ap*S.cos(angle)+bp*S.sin(angle))
    zz = weight*(az*S.cos(angle)+bz*S.sin(angle))
    amp = S.trigsimp(eps*le*r**4*(ph*S.diff(zz, r)-S.diff(ph, r)*zz))
    zero("dimensioned_outward_amplitude_normalization", amp-eps*le*B**3*omega*(ap*bz-bp*az))
    scale = S.symbols("scale", positive=True)
    prefactor = eps*le*B**3*omega*(ap*bz-bp*az)
    scaled = prefactor.subs({eps: scale**2*eps, le: le/scale**2,
                            B: scale*B, az: az/scale**2, bz: bz/scale**2}, simultaneous=True)
    zero("compatibility_moment_length_dimension", scaled-scale*prefactor)


def frequency_bracket():
    k = 2000
    z0 = R(1, k**10)
    pi_low, pi_up = R(31, 10), R(22, 7)
    zl = R(999, 1000)*z0
    check("vacuum_root_positive_lower_bracket", zl+11*zl**2/(360*pi_low) < z0)
    check("vacuum_root_upper_bracket", z0+11*z0**2/(360*pi_up) > z0)
    low_frequency_squared = 60*pi_low/z0+R(17, 12)
    upper_frequency_squared = 60*pi_up/zl+R(17, 12)
    check("benchmark_frequency_lower", low_frequency_squared > (4*10**17)**2)
    check("benchmark_frequency_upper", upper_frequency_squared < (45*10**16)**2)
    check("channel_above_continuum_threshold", 60*pi_low/z0+R(11, 3) > R(9, 4))


def hashes(data):
    return {name: hashlib.sha256(raw).hexdigest() for name, raw in data.items()}


def main():
    before = {name: (ROOT/name).read_bytes() for name in INPUTS}
    existing = json.loads(RECEIPT.read_text(encoding="utf-8")) if RECEIPT.exists() else {"runs": []}
    run = {"started_utc": datetime.now(timezone.utc).isoformat(), "contract": CONTRACT,
           "command": sys.orig_argv, "dont_write_bytecode": sys.dont_write_bytecode, "cwd": str(Path.cwd()),
           "python": platform.python_version(), "sympy": S.__version__, "execution_inputs_sha256": hashes(before)}
    with LOG.open("a", encoding="utf-8", newline="\n") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        print("RUN "+json.dumps(run, sort_keys=True), flush=True)
        try:
            check("actual_bytecode_disabled", sys.dont_write_bytecode and "-B" in sys.orig_argv)
            check("inputs_strict_UTF8_LF", all(raw.decode("utf-8") is not None and b"\r" not in raw
                  and not raw.startswith(b"\xef\xbb\xbf") for raw in before.values()))
            trace_identities()
            static_and_seed_identities()
            liouville_and_channel()
            green_identities()
            frequency_bracket()
        except Exception:
            traceback.print_exc()
            check("uncaught_exception", False)
        run["final_inputs_sha256"] = hashes({name: (ROOT/name).read_bytes() for name in INPUTS})
        print("FINAL_INPUTS "+json.dumps(run["final_inputs_sha256"], sort_keys=True), flush=True)
        check("frozen_inputs_unchanged", run["execution_inputs_sha256"] == run["final_inputs_sha256"])
        run["checks"] = ROWS
        run["failed_checks"] = [row for row in ROWS if not row["passed"]]
        code = int(bool(run["failed_checks"]))
        run["status"] = "failed" if code else "passed"
        run["terminal_exit_status"] = code
        run["finished_utc"] = datetime.now(timezone.utc).isoformat()
        if code:
            run["failed_inputs_base64"] = {name: base64.b64encode(raw).decode("ascii") for name, raw in before.items()}
            for name, raw in before.items():
                print(json.dumps({"path": name, "sha256": hashlib.sha256(raw).hexdigest(),
                                  "base64": base64.b64encode(raw).decode("ascii")}, sort_keys=True))
        print("TERMINAL "+json.dumps({"status": run["status"], "checks": len(ROWS),
              "failures": len(run["failed_checks"]), "terminal_exit_status": code,
              "finished_utc": run["finished_utc"]}, sort_keys=True), flush=True)
    existing["runs"].append(run)
    RECEIPT.write_text(json.dumps(existing, indent=2, sort_keys=True)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": run["status"], "checks": len(ROWS), "failures": len(run["failed_checks"]),
                      "receipt": str(RECEIPT)}, sort_keys=True))
    return code


if __name__ == "__main__":
    sys.exit(main())
