"""Check the limited support calculation; never label it a stable wormhole.

Run with Python, SymPy and mpmath. No network, Lean, or source edits.
The JSON binds this checker, its ordinary argument, and its source review.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import platform

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / "notes/static-formation-support.md"
REVIEW = ROOT / "evidence/static-formation-source-review.json"
OUT = ROOT / "evidence/static-formation-support.json"


def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def main():
    r, s, eta, mass = sp.symbols("r s eta M", positive=True)
    u = sp.Function("u")(r)
    potential = 1 - mass/r - u**2/(4*r**2) - mass**2/u**2
    # Independently differentiate the conservation law m_s'=-8*pi*r*p.
    sigma = sp.Function("Sigma")(r)
    pressure = sp.Function("p")(r)
    ms_prime = -8*sp.pi*r*pressure
    ms_second = sp.diff(ms_prime, r).subs(
        sp.diff(pressure, r), -2*eta*(sigma+pressure)/r)
    sigma0 = (1-s)/(4*sp.pi)
    p0 = (1-s)**2/(16*sp.pi*s)
    stat = {r: 1, sigma: sigma0, pressure: p0}
    u1 = sp.simplify(ms_prime.subs(stat, simultaneous=True))
    u2 = sp.simplify(ms_second.subs(stat, simultaneous=True))
    values = {u: 1-s, sp.diff(u, r): u1,
              sp.diff(u, r, 2): u2, r: 1, mass: (1-s**2)/2}
    derivatives = [sp.factor(sp.diff(potential, r, n).subs(
        values, simultaneous=True)) for n in range(3)]
    require(derivatives[:2] == [0, 0], "Static potential/force nonzero")
    curvature = (12*eta*s**3 + 4*eta*s**2 + 3*s**3-s**2-s-1)/(2*s**2)
    require(sp.simplify(derivatives[2]-curvature) == 0, "Curvature mismatch")
    critical = (1-s)*(3*s**2+2*s+1)/(4*s**2*(3*s+1))
    require(sp.simplify(curvature.subs(eta, critical)) == 0, "Wrong radial root")

    # Independent equilibrium-sequence criterion, 2604.05980v1 (33).
    beta = (1-s)/(4*s)
    gamma_crit = sp.Rational(3, 2)+4*beta+4*beta**2
    require(sp.simplify(critical-gamma_crit*beta/(1+beta)) == 0,
            "Sequence criterion disagrees after thermodynamic conversion")
    causal_poly = 15*s**3+3*s**2-s-1
    require(sp.simplify((1-critical)*4*s**2*(3*s+1)-causal_poly) == 0,
            "Causal-boundary polynomial mismatch")
    require(causal_poly.subs(s, sp.Rational(39529, 100000)) < 0 and
            causal_poly.subs(s, sp.Rational(39530, 100000)) > 0,
            "Root bracket failed")
    positive_roots = [x for x in sp.nroots(causal_poly, maxsteps=100)
                      if abs(sp.im(x)) < 1e-12 and sp.re(x) > 0]
    require(len(positive_roots) == 1, "Unexpected positive roots")
    s_star = float(sp.re(positive_roots[0]))

    sample = {s: sp.Rational(2, 5), eta: sp.Rational(99, 100)}
    require(sp.simplify(curvature.subs(sample)) == sp.Rational(81, 1000),
            "Exact subluminal example failed")
    require(critical.subs(s, sample[s]) == sp.Rational(171, 176), "eta root")
    require(beta.subs(s, sample[s]) == sp.Rational(3, 8), "pressure ratio")
    require(sp.simplify(beta.subs(s, sp.sqrt(2)-1)) == sp.sqrt(2)/4,
            "Confinement pressure margin")
    # The derivative is strictly negative, so smaller s increases beta.
    require(sp.simplify(sp.diff(beta, s)+1/(4*s*s)) == 0, "Monotonicity")

    mp.mp.dps = 70
    sn = mp.mpf(2)/5
    cn = (1-sn**2)/2
    sig0 = (1-sn)/(4*mp.pi)
    pres0 = (1-sn)**2/(16*mp.pi*sn)
    numerical = []
    for en in [mp.mpf(99)/100, mp.mpf(9)/10]:
        # Exact integration of Sigma'=-2(Sigma+p)/r for the affine EOS.
        def v(x):
            sig = ((sig0+pres0)*x**(-2*(1+en)) + en*sig0-pres0)/(1+en)
            ms = 4*mp.pi*x*x*sig
            return 1-cn/x-ms*ms/(4*x*x)-cn*cn/(ms*ms)
        exact = (12*en*sn**3+4*en*sn**2+3*sn**3-sn**2-sn-1)/(2*sn**2)
        rows = []
        for h in [mp.mpf(1)/1000, mp.mpf(1)/2000, mp.mpf(1)/4000]:
            fd = (v(1+h)-2*v(1)+v(1-h))/(h*h)
            rows.append({"h": float(h), "curvature": float(fd),
                         "absolute_error": float(abs(fd-exact))})
        require(rows[-1]["absolute_error"] < 2e-5, "Potential difference error")
        require(all(rows[i+1]["absolute_error"] < rows[i]["absolute_error"]/3
                    for i in range(2)), "No second-order refinement")
        numerical.append({"eta": float(en), "expected": float(exact), "rows": rows})

    def endpoint(x):
        if x <= 1:
            return sn*(1-x)
        return mp.sqrt(1-2*cn/x)-(sn+x-1)
    h = mp.mpf("1e-8")
    left = (endpoint(1)-endpoint(1-h))/h
    right = (endpoint(1+h)-endpoint(1))/h
    require(abs(left+sn) < mp.mpf("1e-50"), "Interior slope")
    require(abs(right-(cn/sn-1)) < mp.mpf("1e-6"), "Exterior slope")
    require(endpoint(1-h)>0 and endpoint(1+h)>0 and endpoint(1)==0,
            "Cusp energy not locally minimal")

    # Conditional finite-thickness bound: supplied coefficients, not physics fits.
    b0, aa, bb, dd = mp.mpf(3)/8, mp.mpf(2), mp.mpf(3), mp.mpf(4)
    threshold = (2*b0/(aa+mp.sqrt(aa*aa+4*b0*(bb+dd))))**2
    thickness = []
    for epsilon in [mp.mpf("1e-2"), mp.mpf("1e-4"), mp.mpf("1e-6")]:
        angular = int(mp.ceil(1/mp.sqrt(epsilon)))
        kk = mp.sqrt(angular*(angular+1))
        upper = -b0*kk**2+aa*kk+bb+dd*epsilon**2*kk**4
        require(upper < 0, "Illustrative discrete-mode sign failed")
        thickness.append({"epsilon": float(epsilon), "ell": angular,
                          "K": float(kk), "k_h": float(epsilon*kk),
                          "W_squared_upper_bound": float(upper)})
    require(all(thickness[i+1]["k_h"] < thickness[i]["k_h"] for i in range(2)),
            "Chosen modes do not move below the thickness cutoff")
    discrepancy = sp.factor(curvature-2*s*(2*eta+1))
    require(discrepancy != 0, "Printed-source discrepancy unexpectedly vanished")
    result = {
        "status": "checks_passed_candidate_rejected",
        "date": "2026-10-12",
        "claim": "Radial overlap exists; positive-pressure fluid support has known warping instability",
        "stable_wormhole": False, "formation_rate_computed": False,
        "full_mode_solver_reproduced": False, "novelty_accepted": False,
        "python": platform.python_version(), "sympy": sp.__version__,
        "mpmath": mp.__version__,
        "normalized_sha256": {str(p.relative_to(ROOT)): digest(p)
                              for p in [Path(__file__).resolve(), NOTE, REVIEW]},
        "symbolic": {"potential_derivatives": list(map(str, derivatives)),
                     "eta_critical": str(critical),
                     "equilibrium_sequence_comparison": "exact equality",
                     "printed_1708_04588_eq40_difference": str(discrepancy)},
        "interval": {"C_lower_exclusive": float(sp.sqrt(2)-1),
                     "C_upper_exclusive": (1-s_star*s_star)/2,
                     "s_star": s_star, "DEC_s_min_inclusive": 0.2},
        "example": {"s": "2/5", "C": "21/50", "eta": "99/100",
                    "eta_crit": "171/176", "R2_V_second": "81/1000",
                    "pressure_over_density": "3/8",
                    "endpoint_slopes": [float(left), float(right)],
                    "eikonal_omega2_over_k2": "-3/8"},
        "independent_potential_checks": numerical,
        "conditional_thickness_check": {
            "assumptions_derived_for_physical_material": False,
            "illustrative_coefficients": {"beta0": 0.375, "A": 2, "B": 3, "D": 4},
            "continuum_sufficient_epsilon_threshold": float(threshold),
            "discrete_modes": thickness},
        "limitations": ["Probe string", "Cusp is not a smooth orbit",
                        "Radial test is not nonspherical stability",
                        "No matched black-hole pair or quantum state",
                        "Thickness implication requires its stated mode and error bounds"]}
    OUT.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ["status", "interval", "example"]}, indent=2))


if __name__ == "__main__":
    main()
