"""S146: outer mixed scalar response, cubic allocation and uniform budgets.

Exact algebra and deterministic numerical diagnostics do not certify the
continuum domains, heat remainder, or physical metric-transfer hypotheses.
"""
from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import json
import math
import sys
import traceback

import mpmath as mp
import numpy as np
import sympy as s
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[1]
STEM = "outer-finite-background-response"
INPUTS = [
    "checks/outer_finite_background_response.py",
    "notes/outer-finite-background-response.md",
    "notes/stationary-centered-canonical-response.md",
    "notes/stationary-centered-exterior-source.md",
    "notes/stationary-massive-boundary-layer.md",
    "notes/general-interior-dirac-response.md",
    "notes/retained-benchmark-parameter-recovery.md",
]
ROWS = []
NUMERICS = []


def check(name, passed, detail=""):
    assert not any(row["name"] == name for row in ROWS)
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not row["passed"]:
        raise AssertionError(name)


def symbolic_checks():
    F, G, Fu, Fw, Fuw, Gu, Gw, Guw, t, u0, w0 = s.symbols(
        "F G Fu Fw Fuw Gu Gw Guw t u0 w0", nonzero=True)
    a, b = s.symbols("a b")
    FF = F+a*Fu+b*Fw+a*b*Fuw
    GG = G+a*Gu+b*Gw+a*b*Guw
    logJ = a*u0+b*w0-s.log(FF)/2+t*GG**2/(4*FF)
    Lu = u0-Fu/(2*F)+t*(2*G*Gu/F-G**2*Fu/F**2)/4
    Lw = w0-Fw/(2*F)+t*(2*G*Gw/F-G**2*Fw/F**2)/4
    Luw = -(Fuw/F-Fu*Fw/F**2)/2+t*(
        2*(Gu*Gw+G*Guw)/F-2*G*(Gu*Fw+Gw*Fu)/F**2
        -G**2*Fuw/F**2+2*G**2*Fu*Fw/F**3)/4
    at0 = {a: 0, b: 0}
    check("exact_first_log_response",
          s.simplify(s.diff(logJ, a).subs(at0)-Lu) == 0)
    check("exact_mixed_log_response",
          s.simplify(s.diff(logJ, a, b).subs(at0)-Luw) == 0)
    JJ = s.exp(a*u0+b*w0)*FF**(-s.Rational(1, 2))*s.exp(t*GG**2/(4*FF))
    check("exact_mixed_integrand_product",
          s.simplify(s.diff(JJ, a, b).subs(at0)/JJ.subs(at0)-Luw-Lu*Lw) == 0)
    constant = {Fu: 2*F, Fw: 2*F, Fuw: 4*F,
                Gu: G, Gw: G, Guw: G, u0: 1, w0: 1}
    check("constant_warp_scalar_invariance",
          s.simplify(Lu.subs(constant)) == 0 and
          s.simplify(Lw.subs(constant)) == 0 and
          s.simplify(Luw.subs(constant)) == 0)
    check("omitted_observation_factor_is_detectable",
          s.simplify((Luw+(Lu-u0)*(Lw-w0)).subs(constant)) == 1)
    yy, p, g, a2, a3, a4, v4 = s.symbols("Y p gamma a2 a3 a4 v4")
    inverse_f = (1+a2*yy**2+a3*yy**3+a4*yy**4)**-1
    physical_a = inverse_f*s.exp(-p*g*yy**3-v4*yy**4)
    cubic = s.expand(s.series(physical_a, yy, 0, 4).removeO()).coeff(yy, 3)
    check("finite_mass_cubic_is_affine", s.simplify(cubic+a3+p*g) == 0)
    check("mixed_singular_coefficient_vanishes",
          s.diff(cubic, p, 2) == 0 and s.diff(cubic, p, v4) == 0)
    check("geodesic_boundary_retained", s.diff(physical_a, yy).subs(yy, 0) == 0)
    Y = s.symbols("Y", positive=True)
    profile = Y**4*(1-s.log(Y))
    for j in range(6):
        normalized = s.expand(s.diff(profile, Y, j)/Y**(4-j))
        z = s.symbols("z")
        coeff = s.Poly(normalized.subs(s.log(Y), z), z)
        # |a+b log Y| <= (|a|+|b|)(1+|log Y|).
        reserve = sum(abs(c) for c in coeff.all_coeffs())
        # At j=0 the exact expression already is L, without a factor two.
        good = j == 0 or reserve <= math.factorial(j)*100**j
        check("weighted_real_c5_example_"+str(j), good, str(normalized))
    check("weighted_example_is_not_closed_c5",
          s.simplify(s.diff(profile, Y, 5)+24/Y) == 0)
    z = s.symbols("z", nonnegative=True)
    check("exact_cubic_bridge_moment",
          s.integrate((2*z*(1-z))**2, (z, 0, 1)) == s.Rational(2, 15))
    def erfc_moment(n):
        return 2**(s.Rational(n+1, 2))*s.gamma(s.Rational(n+2, 2))/(
            (n+1)*s.sqrt(s.pi))
    def gaussian_moment(n):
        return 2**s.Rational(n-1, 2)*s.gamma(s.Rational(n+1, 2))
    radial = -erfc_moment(3)-3*erfc_moment(1)+s.sqrt(2/s.pi)*(
        gaussian_moment(2)+2*gaussian_moment(0))
    check("exact_cubic_radial_moment", s.simplify(radial) == s.Rational(3, 4))
    check("exact_full_cubic_coefficient",
          s.simplify(2*radial*s.Rational(2, 15)/s.pi) == 1/(5*s.pi))


def rational_and_integral_checks():
    mp.mp.dps = 70
    delta = mp.mpf(".02")
    check("complex_cap_cone",
          mp.exp(-2*delta)*mp.cos(2*delta) > mp.mpf(".95") and
          mp.exp(2*delta) < mp.mpf("1.05") and
          (2+delta)*mp.exp(delta) < mp.mpf("2.1"))
    check("complex_exponent_rate", Q(21, 10)**2/(4*Q(95, 100)) < 2)
    const = Q(100)+Q(210)/(2*Q(95, 100))
    slope = (2*Q(21, 10)*310/Q(95, 100)
             + Q(21, 10)**2*210/Q(95, 100)**2)/4
    check("global_first_log_bound", const < 1000 and slope < 1000,
          {"constant": str(const), "t_coefficient": str(slope)})
    check("global_kernel_prefactor",
          mp.mpf("1.1")*1000/(2*mp.pi) < 10000)
    check("long_time_angular_budget", 10**4*10**10*96*5100 < 10**20)
    rho = Q(1, 10000)
    # log(10)<2.303 is also verified, not silently rounded down.
    check("logarithm_guard", mp.log(10) < mp.mpf("2.303"))
    Lmax = Q(3303, 1000)
    check("weighted_jet_cone",
          15*rho*Q(1, 10)**4*Lmax < Q(5, 10**7) and
          3750*rho*Q(1, 10)**3*Lmax < Q(13, 10000))
    check("analytic_product_reserve",
          mp.exp(mp.mpf(".002002"))*mp.mpf(".95")**-10*2 < 10)
    moment5 = 1+mp.quad(lambda x: 10*(1+x)**4*mp.exp(-x*x), [0, 1, 3, mp.inf])
    moment10 = 1+mp.quad(
        lambda x: 2*mp.exp(-x*x+x)*(10*(1+x)**9+(1+x)**10),
        [0, 1, 3, 8, mp.inf])
    check("bridge_fifth_moment", moment5 < 3000, mp.nstr(moment5, 24))
    check("bridge_large_short_moment", moment10 < 10**9, mp.nstr(moment10, 24))
    check("safe_leading_remainder",
          3000*10**15*3000/(2*mp.pi) < 10**22)
    weighted_budget = sum([
        Q(96, 10)*10**23, Q(3072, 1000)*10**17, 10**18,
        Q(1536, 1000)*10**18, 10**20, 2*10**17])
    mass_budget = Q(24)*5*10**23*Q(1, 5)+10**20+10**19+10**18
    check("weighted_cap_budget", weighted_budget < 10**25, str(weighted_budget))
    check("analytic_mass_budget", mass_budget < 10**25, str(mass_budget))
    # Sum+integral bound on the complete tower, with no sampled tail.
    aa = Q(864, 5)
    def monomial_tail(k):
        return sum(Q(math.factorial(k), math.factorial(k-j))/aa**(j+1)
                   for j in range(k+1))
    polynomial = 18**5*(6+4*monomial_tail(6)+8*monomial_tail(8))
    exp_lower = sum(aa**j/math.factorial(j) for j in range(301))
    exact_tower_upper = polynomial/exp_lower
    check("complete_separated_tower_bound",
          exact_tower_upper < Q(12, 10**68),
          "sum plus monotone integral, rational exponential Taylor lower bound")
    check("actual_domain_first_budget", 10**64*Q(12, 10**68)*400 == Q(48, 100))
    check("actual_domain_second_budget", Q(48, 100)*10**8 < 10**8)
    radius = Q(1, 40000)
    max_norm = Q(1, 10**7)+2*radius
    check("parameter_bidisk_inside_cap_tube", max_norm < rho)
    mixed_cap = 10**25*max_norm/radius**2
    check("weighted_cauchy_budget", mixed_cap == Q(8016, 1000)*10**29)
    check("combined_mixed_budget", mixed_cap+10**8 < 10**30)
    check("scaled_scalar_budget", 2*Q(85, 10**34)*10**30 == Q(17, 1000))
    check("actual_global_graph_ratio", 10*10**4*(Q(1, 10**7)+Q(1, 10**8)) < Q(2, 100))
    for Y in [mp.mpf(".05"), mp.mpf(".01"), mp.mpf(".000001")]:
        # Integrate after t=Y^2*z; removes a misleading tiny integration scale.
        safe = mp.quad(lambda z: z**(-mp.mpf(".5")), [0, 1])
        unsafe = mp.quad(lambda z: z**-3*mp.exp(-1/(8*z)), [0, 1])
        fullunsafe = mp.quad(lambda z: z**-3*mp.exp(-1/(8*z)), [0, 1, mp.inf])
        check("scaled_time_integrals_"+str(Y),
              abs(safe-2) < mp.mpf("1e-30") and unsafe < 64 and
              abs(fullunsafe-64) < mp.mpf("1e-40"))
    NUMERICS.append({"bridge_moment5_upper": mp.nstr(moment5, 24),
                     "bridge_moment10_upper": mp.nstr(moment10, 24),
                     "weighted_cap_budget": str(weighted_budget),
                     "mass_budget": str(mass_budget),
                     "mixed_cap_budget": str(mixed_cap)})


def cubic_integral_diagnostic():
    mp.mp.dps = 70
    def integrand(u):
        return -(u**3+3*u)*mp.erfc(u/mp.sqrt(2))+mp.sqrt(2/mp.pi)*(
            u*u+2)*mp.exp(-u*u/2)
    integral = mp.quad(integrand, [0, 1, 4, 10, mp.inf])
    check("independent_gaussian_cubic_integral",
          abs(integral-mp.mpf(".75")) < mp.mpf("1e-50"), mp.nstr(integral, 55))
    for Y in [mp.mpf(".01"), mp.mpf(".2")]:
        for variance in [mp.mpf(".2"), mp.mpf(".5")]:
            result = 2*variance**2/Y*integral
            check("cubic_scale_"+str(Y)+"_"+str(variance),
                  abs(result-3*variance**2/(2*Y)) < mp.mpf("1e-48"))


def finite_path_checks():
    mp.mp.dps = 70
    # Deterministic continuous paths with bridge endpoints. These are
    # integrand tests, not Monte Carlo evidence for the heat expectation.
    cases = [
        ("positive", ".2", ".03", [".35", "-.12", ".08"], False),
        ("crossing", ".07", ".08", ["-.8", ".25", ".1"], False),
        ("complex", ".12", ".04", ["-.5", ".2", "-.08"], True),
    ]
    all_product_gaps = []
    for name, ytext, ttext, amplitudes, iscomplex in cases:
        Y, t = mp.mpf(ytext), mp.mpf(ttext)
        amplitudes = [mp.mpf(v) for v in amplitudes]
        def fields(x):
            ax = abs(x)
            if ax == 0:
                v, vp = mp.mpf(0), mp.mpf(0)
            else:
                v = x**4*(1-mp.log(ax))
                vp = x**3*(3-4*mp.log(ax))
            cubic = ax**3
            cubicp = 3*x*ax
            seed = -mp.mpf(".4")*x*x+mp.mpf(".2")*cubic
            seedp = -mp.mpf(".8")*x+mp.mpf(".2")*cubicp
            bg = mp.mpf(".002")*cubic+mp.mpf(".001")*v
            bgp = mp.mpf(".002")*cubicp+mp.mpf(".001")*vp
            if iscomplex:
                bg += mp.j*mp.mpf(".00003")*v
                bgp += mp.j*mp.mpf(".00003")*vp
            u, up = cubic+mp.mpf(".7")*v, cubicp+mp.mpf(".7")*vp
            w, wp = -mp.mpf(".3")*cubic+mp.mpf("1.2")*v, (
                -mp.mpf(".3")*cubicp+mp.mpf("1.2")*vp)
            return seed+bg, seedp+bgp, u, up, w, wp
        targets = []
        for resolution in [128, 256]:
            nodes, weights = leggauss(resolution)
            data = []
            for node, weight in zip(nodes, weights):
                z = (mp.mpf(float(node))+1)/2
                weight = mp.mpf(float(weight))/2
                path = sum(v*mp.sin((j+1)*mp.pi*z) for j, v in enumerate(amplitudes))
                x = Y+mp.sqrt(t)*path
                data.append((weight, fields(x)))
            fy = fields(Y)
            def scalar(sa=mp.mpf(0), sb=mp.mpf(0)):
                FF = mp.mpf(0)
                GG = mp.mpf(0)
                for weight, (logf, logfp, u, up, w, wp) in data:
                    f = mp.exp(logf+sa*u+sb*w)
                    FF += weight*f*f
                    GG += weight*f*(logfp+sa*up+sb*wp)
                return mp.exp(fy[0]+sa*fy[2]+sb*fy[4])/mp.sqrt(FF)*mp.exp(t*GG*GG/(4*FF))
            sums = [mp.mpf(0) for _ in range(8)]
            for weight, (logf, logfp, u, up, w, wp) in data:
                f = mp.exp(logf)
                fp = f*logfp
                vals = (f*f, fp, 2*f*f*u, 2*f*f*w, 4*f*f*u*w,
                        fp*u+f*up, fp*w+f*wp, fp*u*w+f*(up*w+u*wp))
                sums = [old+weight*value for old, value in zip(sums, vals)]
            F, G, Fu, Fw, Fuw, Gu, Gw, Guw = sums
            Lu = fy[2]-Fu/(2*F)+t*(2*G*Gu/F-G*G*Fu/F**2)/4
            Lw = fy[4]-Fw/(2*F)+t*(2*G*Gw/F-G*G*Fw/F**2)/4
            Luw = -(Fuw/F-Fu*Fw/F**2)/2+t*(
                2*(Gu*Gw+G*Guw)/F-2*G*(Gu*Fw+Gw*Fu)/F**2
                -G*G*Fuw/F**2+2*G*G*Fu*Fw/F**3)/4
            J = scalar()
            target = J*(Luw+Lu*Lw)
            targets.append(target)
            errors = []
            for steptext in [".001", ".0005"]:
                h = mp.mpf(steptext)
                observed = (scalar(h,h)-scalar(h,-h)-scalar(-h,h)+scalar(-h,-h))/(4*h*h)
                errors.append(abs(observed-target))
            check(name+"_mixed_difference_"+str(resolution),
                  max(errors) < mp.mpf("1e-8")*(1+abs(target)),
                  {"target": mp.nstr(target, 24), "errors": [mp.nstr(e, 10) for e in errors]})
            check(name+"_step_refinement_"+str(resolution),
                  errors[1] < mp.mpf(".35")*errors[0] or max(errors) < mp.mpf("1e-40"))
            product_gap = abs(J*Lu*Lw)
            all_product_gaps.append(product_gap)
            NUMERICS.append({"path": name, "nodes": resolution,
                             "mixed_response": mp.nstr(target, 35),
                             "errors": [mp.nstr(e, 16) for e in errors],
                             "omitted_product_error": mp.nstr(product_gap, 16)})
        check(name+"_quadrature_refinement",
              abs(targets[1]-targets[0]) < mp.mpf("2e-5")*(1+abs(targets[1])),
              mp.nstr(abs(targets[1]-targets[0]), 16))
    check("omitted_mixed_product_is_detected", max(all_product_gaps) > mp.mpf("1e-10"))


def hashes():
    return {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    before = hashes()
    error = None
    try:
        symbolic_checks()
        rational_and_integral_checks()
        cubic_integral_diagnostic()
        finite_path_checks()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr)
    after = hashes()
    ok = error is None and before == after and all(row["passed"] for row in ROWS)
    record = {
        "story": "S146", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if ok else "failed", "python": sys.version,
        "numpy": np.__version__, "sympy": s.__version__, "mpmath": mp.__version__,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "numerics": NUMERICS,
        "exception": error,
        "scope": "outer canonical scalar mixed response under the fixed coefficient/domain contract; diagnostics are not continuum proof",
        "unresolved": ["physical metric-to-coefficient transfer", "global exterior inverse norm",
                       "Weyl/LLL/radius mappings", "coupled exterior/state and physical gates"],
    }
    if error:
        record["failed_checker_source_utf8"] = Path(__file__).read_text(encoding="utf-8")
        record["failed_note_source_utf8"] = (ROOT/INPUTS[1]).read_text(encoding="utf-8")
    output = ROOT/"evidence"/f"{STEM}.json"
    history = json.loads(output.read_text(encoding="utf-8")) if output.exists() else {"runs": []}
    history["runs"].append(record)
    output.write_text(json.dumps(history, indent=2)+"\n", encoding="utf-8", newline="\n")
    with (ROOT/"evidence"/f"{STEM}.log").open("a", encoding="utf-8", newline="\n") as log:
        log.write(json.dumps({"created_utc": record["created_utc"], "status": record["status"],
                              "inputs": before})+"\n")
        for row in ROWS:
            log.write(json.dumps(row)+"\n")
        if error:
            log.write(error+"\n")
    print(json.dumps({"status": record["status"], "checks": len(ROWS),
                      "failed": sum(not row["passed"] for row in ROWS)}))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
