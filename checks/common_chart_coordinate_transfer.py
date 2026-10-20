"""S151: weighted common-chart/proper-coordinate transfer.

Manufactured smooth metrics test coordinate identities, not a physical state.
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
from sympy.functions.combinatorial.numbers import stirling

ROOT = Path(__file__).resolve().parents[1]
STEM = "common-chart-coordinate-transfer"
INPUTS = [
    "checks/common_chart_coordinate_transfer.py",
    "notes/common-chart-coordinate-transfer.md",
    "notes/averaged-ward-matching.md",
    "evidence/averaged-ward-matching.json",
    "notes/outer-finite-background-response.md",
    "notes/outer-ward-and-constraint-transfer.md",
    "notes/stationary-centered-canonical-response.md",
    "notes/stationary-completed-metric-response-norm.md",
]
ROWS, NUMERICS = [], []
PROTOCOL = {
    "nodes": [32, 64], "digits": 55, "steps": [".001", ".0005"],
    "observations": [".025", ".05", ".09"],
    "mixed_tolerance": "2e-6", "quadrature_tolerance": "1e-11",
    "refinement_ratio": ".35", "refinement_floor": "1e-30",
}


def check(name, passed, detail=""):
    assert not any(r["name"] == name for r in ROWS)
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not passed:
        raise AssertionError(name)


def predecessor_checks():
    expected = {
        "notes/averaged-ward-matching.md":
            "8ba87d20fca5ea3c51fac7e5cd0bb9f575c682321ba4b854e746dbb1595e091f",
        "checks/averaged_ward_matching.py":
            "9ab874fb43104d01f89a222606455899361562a00f941627ef3e42688317242d",
        "evidence/averaged-ward-matching.json":
            "582466967a1a4c70514f9eb1b0a50ee960c2c6684c5b96ef0caa9653d7ea2341",
        "evidence/averaged-ward-matching.log":
            "b9b3779ca485a0091de01785c468885dc5c31f1f51af4dc1f56870f6db8ff3e8",
    }
    check("s150_committed_artifact_hashes", all(
        sha256((ROOT / p).read_bytes()).hexdigest() == h for p, h in expected.items()))
    receipt = json.loads((ROOT / INPUTS[3]).read_text(encoding="utf-8"))
    last = receipt["runs"][-1]
    check("s150_final_receipt", last["status"] == "passed"
          and len(last["checks"]) == 61 and all(r["passed"] for r in last["checks"])
          and last["unchanged_inputs"])
    check("s150_all_eight_inputs_unchanged", len(last["input_sha256_after"]) == 8
          and all(sha256((ROOT / p).read_bytes()).hexdigest() == h
                  for p, h in last["input_sha256_after"].items()))


def exact_checks():
    st, tt, y = s.symbols("st tt y")
    eb, bp, h, k, Ih, Ik, Ihk = s.symbols("eb bp h k Ih Ik Ihk", nonzero=True)
    yh, yk = -Ih/eb, -Ik/eb
    yhk = -Ihk/eb-h*yk-k*yh-bp*yh*yk
    shift = st*yh+tt*yk+st*tt*yhk
    F = eb*shift+eb*bp*shift**2/2 + st*(Ih+eb*h*shift) + tt*(Ik+eb*k*shift)+st*tt*Ihk
    zero = {st: 0, tt: 0}
    check("inverse_first_implicit", s.simplify(s.diff(F, st).subs(zero)) == 0)
    check("inverse_mixed_implicit", s.simplify(s.diff(F, st, tt).subs(zero)) == 0)
    Db, Ah, Ak, Ahk = s.symbols("Db Ah Ak Ahk")
    B = -Ahk+h*Ak+k*Ah-Db*Ah*Ak
    check("inverse_euler_mixed", s.simplify(yhk.subs({
        Ih: eb*y*Ah, Ik: eb*y*Ak, Ihk: eb*y*Ahk, bp: Db/y})-y*B) == 0)
    gp, gpp, hp, kp = s.symbols("gp gpp hp kp")
    beta = -gp*shift-gpp*shift**2/2-st*(h+hp*shift)-tt*(k+kp*shift)
    beta1 = -h-gp*yh
    beta2 = -hp*yk-kp*yh-gpp*yh*yk-gp*yhk
    check("beta_first_chain", s.simplify(s.diff(beta, st).subs(zero)-beta1) == 0)
    check("beta_mixed_chain", s.simplify(s.diff(beta, st, tt).subs(zero)-beta2) == 0)
    Dg, D2g, Dh, Dk = s.symbols("Dg D2g Dh Dk")
    euler2 = Dh*Ak+Dk*Ah-(D2g-Dg)*Ah*Ak-Dg*B
    check("beta_mixed_euler", s.simplify(beta2.subs({
        Ih: eb*y*Ah, Ik: eb*y*Ak, Ihk: eb*y*Ahk, bp: Db/y,
        gp: Dg/y, gpp: (D2g-Dg)/y**2, hp: Dh/y, kp: Dk/y})-euler2) == 0)
    nh, nk, nhp, nkp = s.symbols("nh nk nhp nkp")
    rho = gp*shift+gpp*shift**2/2+st*(nh+nhp*shift)+tt*(nk+nkp*shift)
    rho2 = nhp*yk+nkp*yh+gpp*yh*yk+gp*yhk
    check("radius_log_first_chain", s.simplify(s.diff(rho, st).subs(zero)-nh-gp*yh) == 0)
    check("radius_log_mixed_chain", s.simplify(s.diff(rho, st, tt).subs(zero)-rho2) == 0)
    check("inverse_acceleration_omission_detected",
          s.simplify(beta2-(-hp*yk-kp*yh-gpp*yh*yk)+gp*yhk) == 0
          and s.simplify(gp*yhk) != 0)
    Y, m, n3, d, c0, a2 = s.symbols("Y m n3 d c0 a2")
    yy = Y-m*Y**4/4
    aa = (1+m*yy**3)*(1+a2*yy**2+c0*yy**3)
    uu = (d*yy+d*a2*yy**3/3+d*c0*yy**4/4)*(1+n3*yy**3)
    aser = s.series(aa, Y, 0, 4).removeO().expand()
    user = s.series(uu, Y, 0, 5).removeO().expand()
    check("inverse_cubic_series", s.series(yy+m*yy**4/4-Y, Y, 0, 5).removeO() == 0)
    check("warp_cubic_series", aser.coeff(Y, 3) == c0+m)
    check("radius_quartic_series", s.simplify(user.coeff(Y, 4)-d*(c0/4+n3-m/4)) == 0)
    al = s.series(s.log(d*aser/s.diff(user, Y)), Y, 0, 4).removeO().expand()
    check("lapse_quadratic_zero", al.coeff(Y, 2) == 0)
    check("lapse_cubic_correlation", al.coeff(Y, 3) == 2*m-4*n3)
    check("correlated_radius_cubic_removes_lapse", al.coeff(Y, 3).subs(n3, m/2) == 0)
    check("independent_radius_cubic_leaves_lapse", al.coeff(Y, 3).subs(n3, m/4) == m)
    BC, x = s.symbols("BC x", positive=True)
    check("physical_mass_parameter", s.simplify((-BC*m/d)/(BC*d**2)+m/d**3) == 0)
    rh, rk, rhk, l, Dh0, Dk0, Dhk0 = s.symbols("rh rk rhk l Dh0 Dk0 Dhk0")
    u = s.exp(st*rh+tt*rk+st*tt*rhk)
    up = u*(l+st*Dh0+tt*Dk0+st*tt*Dhk0)/l
    target = rhk+rh*rk+(Dhk0+rh*Dk0+rk*Dh0)/l
    check("relative_radius_mixed", s.diff(u, st, tt).subs(zero) == rhk+rh*rk)
    check("relative_radius_gradient_mixed", s.simplify(s.diff(up, st, tt).subs(zero)-target) == 0)
    freq, center, primitive = s.symbols("freq center primitive", positive=True)
    osc = s.cos(freq*(Y-center))/freq**5
    bilinear = s.diff(osc, Y)*(Y-center+primitive)
    growth = s.diff(bilinear, Y, 5).subs(Y, center)
    check("c5_bilinear_unbounded_exact", s.simplify(growth+freq*primitive) == 0)
    check("c5_bilinear_growth_limit", s.limit(-growth, freq, s.oo) == s.oo)
    check("ordinary_euler_conversion_reserve", all(math.factorial(j) <= 100**j for j in range(6)))


def budget_checks():
    mp.mp.dps = PROTOCOL["digits"]
    r, eps = Q(11, 100), Q(1, 1000)
    check("mass_euler_norm", sum(Q(3**i, math.factorial(i)) for i in range(8)) < 21)
    check("log_weight_embedding", mp.mpf(".11")*(1-mp.log(mp.mpf(".11"))) < mp.mpf(".36"))
    check("exponential_jet_margin", mp.exp(mp.mpf(str(float(r**3*eps)))) < mp.mpf("1.001"))
    check("averaging_integral_weights", Q(1, 4) < Q(251, 1000)/Q(1001, 1000)
          and Q(1001, 1000)**2/7 < Q(144, 1000))
    logcap = Q(1001, 1000)/4/(1-Q(251, 1000)*eps*r**3)
    check("log_coordinate_margin", logcap < Q(252, 1000))
    check("coordinate_multiplier_margin",
          mp.exp(mp.mpf(str(float(Q(1252, 1000)*eps*r**3)))) < mp.mpf("1.001"))
    check("inverse_first_margin", Q(1001, 1000)**2/4 < Q(251, 1000))
    bcap = Q(144, 1000)+2*Q(251, 1000)+7*eps*Q(251, 1000)**2*r**3
    check("inverse_mixed_margin", bcap < Q(65, 100), str(bcap))
    check("target_collar_inside_image", mp.mpf(".11")*mp.exp(-mp.mpf(".001")*mp.mpf(".11")**3) > mp.mpf(".1"))
    check("weight_ratio_margin", mp.exp(mp.mpf("7")*mp.mpf(".252")*mp.mpf(".001")*mp.mpf(".11")**3) < mp.mpf("1.01"))
    check("composition_five_jet_margin", Q(101, 100)*Q(2001, 1000)**5 < 33)
    check("cubic_exponential_jet_margin", mp.mpf(".756")*mp.exp(mp.mpf(".756")*mp.mpf(".001")*mp.mpf(".11")**3) < mp.mpf(".757"))
    check("cubic_reexpression_margin", 21*Q(757, 1000) < 16 and Q(16, 21) < Q(762, 1000))
    # Affine in M>=1: compare both the constant and M coefficients.
    first0 = 33*(1+Q(762, 1000)*eps+7*eps*Q(251, 1000))
    first1 = 33*Q(251, 1000)
    check("proper_first_response_margin", first0 < 50 and first1 < 50)
    second0 = 33*(Q(3514, 1000)+49*eps*Q(251, 1000)**2+7*eps*Q(65, 100))
    second1 = 33*(Q(251, 1000)**2+Q(65, 100))
    check("proper_mixed_response_margin", second0 < 150 and second1 < 150)
    check("physical_mass_norm_margin", 21*Q(2, 5)**3 > 1)
    check("radius_precomposition_reserves", 3*(Q(2, 1000)+Q(2, 1000)) < 1
          and 3*(Q(12, 10**6)+Q(2, 10**6)+Q(4, 10**6)) < 1)
    check("radius_relative_and_scaling_reserves", 5+9 == 14 and 14 < 5**2)
    check("radius_seed_first_derivative_guards", Q(42, 100)**2+Q(5, 4) < Q(6, 5)**2
          and Q(2, 5)**2/2 > Q(28, 100)**2)
    check("radius_seed_logarithmic_coefficients", Q(6, 5)/Q(28, 100) < 5
          and r*Q(5, 2)/Q(28, 100) < 1
          and r*r*250/Q(28, 100) < 11)
    check("radius_seed_explicit_constant", 5+1+25 == 31
          and 5+3+11+75+15+250 == 359
          and 5+31+Q(359, 2) < 500 and 26+359+31 < 500)
    R = {1: Q(1800)}
    for l in range(2, 8):
        R[l] = 9*math.factorial(l-1)*200**(l-1)*r**(l-2)
    E = {j: sum(int(stirling(j, l, kind=2))*R[l] for l in range(1, j+1)) for j in range(1, 8)}
    M = max(Q(1), sum(E[j+1]/math.factorial(j) for j in range(7)),
            sum((E[j+2]+E[j+1])/math.factorial(j) for j in range(6)))
    V = {l: Q(10*math.factorial(4)*100**4, math.factorial(4-l)) for l in range(5)}
    V[5] = 10*math.factorial(5)*100**5*r
    P = V[0]+sum(sum(int(stirling(j, l, kind=2))*V[l] for l in range(1, j+1))/math.factorial(j)
                 for j in range(1, 6))
    C1, C2 = 51+50*M+P, 150*(1+M)
    check("explicit_seed_constants_finite", M > 1 and P > 0 and C1 > 0 and C2 > 0,
          {"M": str(M), "P": str(P), "C1": str(C1), "C2": str(C2),
           "sufficient_local_background": str(Q(1, 10**7)/C1)})
    check("scaled_coordinate_contract", C1**2 > C2 and C1 > 2500)
    check("seed_complex_first_derivative_margin", Q(25, 2)*Q(61, 10) < 81)
    check("unbounded_fixture_growth", [Q(n, 20) for n in [100, 1000, 10000]] == [5, 50, 500])
    NUMERICS.append({"seed_constants": {"M": str(M), "P": str(P), "C1": str(C1), "C2": str(C2)},
                     "M_R": 500,
                     "existing_C5_endpoint_in_C7_ball_established": False})


def numerical_checks():
    mp.mp.dps = PROTOCOL["digits"]
    Ysym = s.symbols("y", positive=True)
    polynomials = {
        "b": Ysym**3/5+Ysym**4/10, "h": 3*Ysym**3+2*Ysym**4,
        "k": -2*Ysym**3+Ysym**5,
        "n": Ysym**3/10+Ysym**4/5, "nh": 3*Ysym**3/2+2*Ysym**4/5,
        "nk": -Ysym**3+Ysym**5/2,
        "L": Ysym**2/5-Ysym**3/10,
        "R": s.log(s.Rational(41, 100)*Ysym)+Ysym**2/15-Ysym**3/40,
    }
    fun = {name: [s.lambdify(Ysym, s.diff(expr, Ysym, j), "mpmath") for j in range(3)]
           for name, expr in polynomials.items()}
    all_targets = []
    for count in PROTOCOL["nodes"]:
        # Nodes and weights must share the working precision. Binary64
        # quadrature data left a floor in the tiny inverse mixed response.
        nodes, weights = mp.gauss_quadrature(count, "legendre")

        def quad(f, y):
            return y/2*sum(w*f(y*(1+z)/2) for z, w in zip(nodes, weights))

        def values(Y, st, tt):
            field = lambda z: fun["b"][0](z)+st*fun["h"][0](z)+tt*fun["k"][0](z)
            root = mp.findroot(lambda z: quad(lambda w: mp.exp(field(w)), z)-Y, (Y*.99, Y*1.01))
            beta = fun["L"][0](Y)-fun["L"][0](root)-field(root)
            radiuslog = fun["n"][0](root)+st*fun["nh"][0](root)+tt*fun["nk"][0](root)+fun["R"][0](root)
            rho = radiuslog-fun["R"][0](Y)
            return [root, beta, rho, mp.exp(radiuslog)]

        targets_for_nodes = []
        for Ytext in PROTOCOL["observations"]:
            Y = mp.mpf(Ytext)
            base = values(Y, mp.mpf(0), mp.mpf(0))
            y = base[0]
            b, h, k = [fun[name][0](y) for name in ["b", "h", "k"]]
            Ih = quad(lambda z: mp.exp(fun["b"][0](z))*fun["h"][0](z), y)
            Ik = quad(lambda z: mp.exp(fun["b"][0](z))*fun["k"][0](z), y)
            Ihk = quad(lambda z: mp.exp(fun["b"][0](z))*fun["h"][0](z)*fun["k"][0](z), y)
            yh, yk = -mp.exp(-b)*Ih, -mp.exp(-b)*Ik
            yhk = -mp.exp(-b)*Ihk-h*yk-k*yh-fun["b"][1](y)*yh*yk
            gp = fun["b"][1](y)+fun["L"][1](y)
            gpp = fun["b"][2](y)+fun["L"][2](y)
            betah = -h-gp*yh
            betahk = -fun["h"][1](y)*yk-fun["k"][1](y)*yh-gpp*yh*yk-gp*yhk
            rp = fun["n"][1](y)+fun["R"][1](y)
            rpp = fun["n"][2](y)+fun["R"][2](y)
            rhoh = fun["nh"][0](y)+rp*yh
            rhok = fun["nk"][0](y)+rp*yk
            rhohk = fun["nh"][1](y)*yk+fun["nk"][1](y)*yh+rpp*yh*yk+rp*yhk
            targets = [yhk, betahk, rhohk, base[3]*(rhohk+rhoh*rhok)]
            first_targets = [yh, betah, rhoh, base[3]*rhoh]
            errs = [[] for _ in targets]
            first_errs = [[] for _ in targets]
            for steptext in PROTOCOL["steps"]:
                step = mp.mpf(steptext)
                four = [values(Y, st, tt) for st, tt in [(step, step), (step, -step), (-step, step), (-step, -step)]]
                two = [values(Y, st, mp.mpf(0)) for st in [step, -step]]
                for j in range(4):
                    mixed = (four[0][j]-four[1][j]-four[2][j]+four[3][j])/(4*step**2)
                    first = (two[0][j]-two[1][j])/(2*step)
                    errs[j].append(abs(mixed-targets[j]))
                    first_errs[j].append(abs(first-first_targets[j]))
            tag = f"{count}_{Ytext}"
            for j, name in enumerate(["inverse", "warp", "radius_log", "radius"]):
                check(f"{name}_mixed_difference_{tag}", max(errs[j]) < mp.mpf(PROTOCOL["mixed_tolerance"])*(1+abs(targets[j])),
                      {"target": mp.nstr(targets[j], 25), "errors": [mp.nstr(e, 16) for e in errs[j]]})
                check(f"{name}_mixed_refinement_{tag}", errs[j][1] < mp.mpf(".35")*errs[j][0]
                      or max(errs[j]) < mp.mpf("1e-30"))
            check(f"all_first_differences_{tag}", all(
                max(es) < mp.mpf(PROTOCOL["mixed_tolerance"])*(1+abs(target))
                for es, target in zip(first_errs, first_targets)))
            adaptive_residual = mp.quad(lambda z: mp.exp(fun["b"][0](z)), [0, y])-Y
            check(f"independent_adaptive_inverse_{tag}", abs(adaptive_residual) < mp.mpf("1e-11"))
            targets_for_nodes.extend(targets)
            NUMERICS.append({"nodes": count, "Y": Ytext, "mixed_targets": [mp.nstr(v, 30) for v in targets],
                             "manufactured_fixture_is_physical_solution": False})
        all_targets.append(targets_for_nodes)
    check("quadrature_refinement_all_coordinate_maps", all(
        abs(a-b) < mp.mpf(PROTOCOL["quadrature_tolerance"])*(1+abs(b))
        for a, b in zip(*all_targets)))


def hashes():
    return {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    before = hashes()
    error = None
    try:
        predecessor_checks()
        exact_checks()
        budget_checks()
        numerical_checks()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    ok = error is None and before == after and all(r["passed"] for r in ROWS)
    record = {
        "story": "S151", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if ok else "failed", "python": sys.version,
        "numpy": np.__version__, "sympy": s.__version__, "mpmath": mp.__version__,
        "protocol": PROTOCOL, "input_sha256_before": before,
        "input_sha256_after": after, "unchanged_inputs": before == after,
        "checks": ROWS, "numerics": NUMERICS, "exception": error,
        "scope": "local weighted C7-to-C5 coordinate transfer and radius/mass correlation",
        "unresolved": ["global common-family/source contract", "cold response",
                       "actual C7 input and global inverse", "other source sectors",
                       "exterior/state solve and all physical completion conditions"],
    }
    if error:
        record["failed_checker_source_utf8"] = Path(__file__).read_text(encoding="utf-8")
        record["failed_note_source_utf8"] = (ROOT/INPUTS[1]).read_text(encoding="utf-8")
    out = ROOT/"evidence"/f"{STEM}.json"
    history = json.loads(out.read_text(encoding="utf-8")) if out.exists() else {"runs": []}
    history["runs"].append(record)
    out.write_text(json.dumps(history, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": record["status"], "checks": len(ROWS), "unchanged_inputs": before == after}), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
