"""S154: divided high metric jets in the actual-domain scalar theorem.

Numerical fixtures check finite-jet transfer, not a Dirac determinant.
"""
from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json
import math
import sys
import traceback

import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = "divided-jet-dirac-response"
INPUTS = [
    "checks/divided_jet_dirac_response.py",
    "notes/divided-jet-dirac-response.md",
    "notes/general-interior-dirac-response.md",
    "notes/coincident-normal-dirac-response.md",
    "notes/cold-response-from-real-endpoints.md",
    "notes/stationary-completed-metric-response-norm.md",
    "evidence/general-interior-dirac-response.json",
    "evidence/cold-response-from-real-endpoints.json",
]
PROTOCOL = {
    "digits": 60, "nodes_per_panel": [32, 64],
    "frequencies": [400, 1000, 4000],
    "segment_endpoints": ["0.0001", "0.005", "0.015625"],
    "max_panel_length": "1/frequency",
    "nonanalytic_join": "0.005",
    "join_panels": "square-distance substitution",
    "normalized_identity_and_refinement_tolerance": "1e-20",
    "remainder_normalization": "2e9*r^5",
    "fixture_is_physical_Dirac_source": False,
}
ROWS, NUMERICS = [], []


def check(name, passed, detail=""):
    assert name not in {r["name"] for r in ROWS}
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not passed:
        raise AssertionError(name)


def predecessor_checks():
    for stem, count in [("general-interior-dirac-response", 61),
                        ("cold-response-from-real-endpoints", 152)]:
        last = json.loads((ROOT/"evidence"/(stem+".json")).read_text(encoding="utf-8"))["runs"][-1]
        check(stem+"_receipt", last["status"] == "passed"
              and len(last["checks"]) == count
              and all(r["passed"] for r in last["checks"]))
        check(stem+"_input_identity", last["unchanged_inputs"]
              and all(sha256((ROOT/p).read_bytes()).hexdigest() == h
                      for p, h in last["input_sha256_after"].items()))


def exp_coefficients(z, order=5):
    out = [Q(1)]
    for n in range(1, order+1):
        out.append(sum(k*z[k]*out[n-k] for k in range(1, n+1))/n)
    return out


def exact_checks():
    x = s.symbols("x")
    z = [Q(0), Q(16, 1000)] + [Q(11, 1000)*150**j for j in range(2, 6)]
    ec = exp_coefficients(z)
    independent = s.series(s.exp(sum(s.Rational(z[j].numerator, z[j].denominator)*x**j
                                    for j in range(1, 6))), x, 0, 6).removeO()
    check("finite_exponential_recursion", all(
        independent.coeff(x, j) == s.Rational(ec[j].numerator, ec[j].denominator)
        for j in range(6)))
    check("real_ray_high_jet_factor", 2*100**2 < 150**2)
    seed = [Q(1026, 1000)] + [Q(26, 1000)*64**j for j in range(1, 6)]
    inner = [Q(1001, 1000)*sum(seed[k]*ec[j-k] for k in range(j+1)) for j in range(6)]
    reserves = {2: 400, 3: 50000, 4: 7000000, 5: 2000000000}
    displayed = {2: 361, 3: 45369, 4: 6276299, 5: 910342924}
    for j in range(2, 6):
        check(f"inner_real_remainder_{j}", inner[j] < displayed[j] < reserves[j], str(inner[j]))
    center_seed = [Q(1), Q(0), Q(8), Q(26, 1000)*64**3, Q(26, 1000)*64**4]
    pc = [None, Q(2, 100), Q(300), Q(50000), Q(7000000)]
    center = [sum(center_seed[k]*ec[j-k] for k in range(j+1)) for j in range(5)]
    for j in range(1, 5):
        check(f"center_coefficient_{j}", center[j] < pc[j], str(center[j]))
    for degree in range(1, 5):
        outside = Q(202, 100)*64**degree + sum(pc[j]*64**(degree+1-j)
                                                               for j in range(2, degree+1))
        check(f"global_taylor_division_{degree}", outside < reserves[degree+1], str(outside))
    bb = [None, Q(2), Q(400), Q(50000), Q(7000000), Q(2000000000)]

    @lru_cache(None)
    def product(n, d):
        if n == 1:
            return bb[d+1]
        if d == n-1:
            return Q(2)**n
        return bb[d-n+2]*2**(n-1) + sum(pc[a]*product(n-1, d-a)
                                                     for a in range(1, d-n+2))
    expected = [Q(2000000000), Q(49140000), Q(446820), Q(2020202, 625)]
    for n, target in enumerate(expected, 1):
        check(f"complete_product_remainder_{n}", product(n, 4) == target)
    rm = Q(102, 100)
    ratio = 5*rm*Q(6, 100)
    check("born_tail_ratio", ratio < 1)
    budget = 20*rm**2*product(1, 4)
    budget += sum(28*math.factorial(5)*10**5*5**(n-2)*rm**(n+1)*product(n, 4)
                  for n in range(2, 5))
    tail = 28*math.factorial(5)*10**5*2**5*5**3*rm**6/(1-ratio)
    budget += tail
    check("whole_born_remainder", budget < Q(184, 100)*10**16, str(budget))
    check("scalar_area_factor", Q(103, 100)*budget < Q(190, 100)*10**16 < 10**17)
    # Independently sum the first tail terms and enclose the exact remainder.
    partial = sum(28*math.factorial(5)*10**5*2**5*5**3*rm**6*ratio**j for j in range(20))
    check("born_tail_finite_plus_remainder", partial+tail*ratio**20 == tail)
    A = 2*(10**8+220)
    T = 2*(10**8+66000)
    B4 = 4*(10**8+26400000)
    p = Q(11, 1000)
    qp = B4+3*A*A+8*p*T+8*p*p*A
    check("complete_contact_polynomial", qp < 2*10**17, str(qp))
    check("contact_scalar_factor", Q(8*100, 960*3)*Q(105, 100)*qp < 10**17)
    eta = Q(85, 10**34)
    K = 2*eta*(Q(2, 3)*10**17 + 7*10**11*24*256**3 + 10**17)
    check("unchanged_complete_complex_bound", K < Q(5, 10**12), str(K))
    check("first_response_constant", 200*K < Q(1, 10**9))
    check("mixed_response_constant", 40000*K < Q(2, 10**7))
    comparison = Q(2, 10**7)*(math.factorial(5)*100**5)**2
    check("raw_norm_rescaling_is_weaker", comparison == Q(288, 100)*10**17)
    # A C5 function need not have a bounded sixth derivative.
    y = s.symbols("y", positive=True)
    frac = y**s.Rational(11, 2)
    check("fractional_fifth_continuous", s.limit(s.diff(frac, y, 5), y, 0, dir="+") == 0)
    check("fractional_sixth_unbounded", s.limit(s.diff(frac, y, 6), y, 0, dir="+") == s.oo)
    NUMERICS.append({"loop_majorant_exact": str(budget), "source_K_exact": str(K),
                     "raw_rescaling_constant": str(comparison),
                     "global_family_admissibility_proved": False})


def numerical_checks():
    mp.mp.dps = PROTOCOL["digits"]
    check("inner_exponential_value", mp.exp(mp.mpf(".016")/64) < mp.mpf("1.001"))
    check("contact_exponential_value", mp.exp(mp.mpf(".044")) < mp.mpf("1.05"))
    nodes = {}
    for n in PROTOCOL["nodes_per_panel"]:
        xs, ws = mp.gauss_quadrature(n, "legendre")
        nodes[n] = (list(xs), list(ws))
    c = mp.mpf(".005")
    radius = mp.mpf(1)/32
    power = mp.mpf("5.5")
    falling = [mp.mpf(1)]
    for j in range(1, 6):
        falling.append(falling[-1]*(power-j+1))

    for frequency in PROTOCOL["frequencies"]:
        n = mp.mpf(frequency)
        amp = min(mp.mpf(".002")/n, mp.mpf("6e9")/n**5)

        def zjet(x, order):
            shift = x-c
            fractional = falling[order]*abs(shift)**(power-order)
            if order % 2:
                fractional *= mp.sign(shift)
            affine = mp.mpf("1e-6") + mp.mpf(".001")*x if order == 0 else (
                mp.mpf(".001") if order == 1 else mp.mpf(0))
            return affine + amp*n**order*mp.sin(n*x+order*mp.pi/2) + fractional

        def sigmajet(x, order):
            if order == 0:
                return -mp.log1p(-x*x/4)
            return mp.factorial(order-1)*((2-x)**(-order)+(-1)**order*(2+x)**(-order))

        def hjet_coefficients(x, order):
            g = [sigmajet(x, 0)+zjet(x, 0)-zjet(mp.mpf(0), 0)]
            g += [(sigmajet(x, j)+zjet(x, j))/mp.factorial(j) for j in range(1, order+1)]
            e = [mp.exp(g[0])]
            for j in range(1, order+1):
                e.append(sum(k*g[k]*e[j-k] for k in range(1, j+1))/j)
            e[0] -= 1
            return e

        # Analytic bounds on the entire interval, including the cap transition.
        span = radius+c
        low0 = mp.mpf("1e-6")+mp.mpf(".001")*radius+amp+span**power
        low1 = mp.mpf(".001")+amp*n+falling[1]*span**(power-1)
        cap1 = low1+120*low0
        jets = [(amp*n**j+falling[j]*span**(power-j))/(mp.factorial(j)*100**j)
                for j in range(2, 6)]
        P5 = max([low0, low1, cap1]+jets)
        check(f"fixture_full_coefficient_cone_{frequency}", P5 < mp.mpf(".011"), mp.nstr(P5, 25))
        check(f"fixture_background_after_scaling_{frequency}", P5/10 < mp.mpf(".001"))
        check(f"fixture_raw_fifth_not_small_{frequency}", abs(zjet(mp.mpf(0), 5))/10 > 1000)
        center = hjet_coefficients(mp.mpf(0), 4)

        for endpoint in PROTOCOL["segment_endpoints"]:
            r = mp.mpf(endpoint)
            direct = (hjet_coefficients(r, 0)[0]-sum(center[j]*r**j for j in range(5)))/(mp.mpf("2e9")*r**5)
            panels = int(mp.ceil(n*r))
            cuts = [mp.mpf(k)/panels for k in range(panels+1)]
            if 0 < c/r < 1:
                cuts.append(c/r)
            cuts = sorted(set(cuts))
            observed = []
            for count in PROTOCOL["nodes_per_panel"]:
                xs, ws = nodes[count]
                total = mp.mpf(0)
                for left, right in zip(cuts, cuts[1:]):
                    mid, half = (left+right)/2, (right-left)/2
                    for point, weight in zip(xs, ws):
                        if left == c/r or right == c/r:
                            v = (point+1)/2
                            t = left+(right-left)*v*v if left == c/r else right-(right-left)*v*v
                            measure = (right-left)*v*weight
                        else:
                            t = mid+half*point
                            measure = half*weight
                        fifth = mp.factorial(5)*hjet_coefficients(r*t, 5)[5]
                        total += measure*(1-t)**4*fifth/mp.factorial(4)
                value = total/mp.mpf("2e9")
                observed.append(value)
                check(f"fifth_integral_identity_{frequency}_{endpoint}_{count}",
                      abs(value-direct) < mp.mpf("1e-20"),
                      {"target": mp.nstr(direct, 25), "error": mp.nstr(abs(value-direct), 10)})
            check(f"fifth_integral_node_refinement_{frequency}_{endpoint}",
                  abs(observed[0]-observed[1]) < mp.mpf("1e-20"))
            check(f"fifth_integral_global_bound_{frequency}_{endpoint}", abs(direct) < 1)
            NUMERICS.append({"frequency": frequency, "endpoint": endpoint,
                             "normalized_remainder": mp.nstr(direct, 25), "panels": len(cuts)-1})

        def capped_h(x):
            r = abs(x)
            if r >= radius:
                return mp.expm1(-zjet(mp.mpf(0), 0))
            if r <= radius/2:
                cutoff = mp.mpf(1)
            else:
                t = (r-radius/2)/(radius/2)
                cutoff = 1-10*t**3+15*t**4-6*t**5
            return mp.expm1(cutoff*(sigmajet(x, 0)+zjet(x, 0))-zjet(mp.mpf(0), 0))

        radii = [mp.mpf("1e-5"), mp.mpf(".001"), c, radius/2, 3*radius/4,
                 radius, 2*radius, mp.mpf(".25"), mp.mpf(1)]
        for degree, reserve in [(1, 400), (2, 50000), (3, 7000000), (4, 2000000000)]:
            ratios = []
            for r in radii:
                for sign in [-1, 1]:
                    x = sign*r
                    error = capped_h(x)-sum(center[j]*x**j for j in range(degree+1))
                    ratios.append(abs(error)/(reserve*r**(degree+1)))
            check(f"capped_global_taylor_{frequency}_{degree}", max(ratios) < 1,
                  mp.nstr(max(ratios), 25))


def hashes():
    return {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    before = hashes()
    error = None
    try:
        predecessor_checks()
        exact_checks()
        numerical_checks()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    ok = error is None and before == after and all(r["passed"] for r in ROWS)
    record = {
        "story": "S154", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if ok else "failed", "python": sys.version,
        "mpmath": mp.__version__, "sympy": s.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "numerics": NUMERICS,
        "exception": error,
        "scope": "actual-domain interior scalar response on divided high-jet coefficient ball",
        "unresolved": ["full common-chart family admissibility", "outer/other source sectors",
                       "full higher-order inverse and coupled solve", "physical completion requirements"],
    }
    if error:
        record["failed_checker_source_utf8"] = Path(__file__).read_text(encoding="utf-8")
        record["failed_note_source_utf8"] = (ROOT/INPUTS[1]).read_text(encoding="utf-8")
    out = ROOT/"evidence"/(STEM+".json")
    history = json.loads(out.read_text(encoding="utf-8")) if out.exists() else {"runs": []}
    history["runs"].append(record)
    out.write_text(json.dumps(history, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": record["status"], "checks": len(ROWS),
                      "unchanged_inputs": before == after}), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
