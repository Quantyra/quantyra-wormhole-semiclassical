"""S163 normal inverse scale and compatible-global lower-bound diagnostics."""
from __future__ import annotations
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import mpmath as mp
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "weighted-normal-inverse-scale"
OWN = ["checks/weighted_normal_inverse_scale.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("stationary-disk-and-normal-inverse", 61),
    ("stationary-global-weighted-map", 50),
    ("stationary-weighted-tt-inverse", 52),
    ("stationary-scalar-outgoing-inverse", 107),
    ("general-background-newton-operator", 174),
]
PROTOCOL = {
    "digits": 100, "identity_tolerance": "1e-70",
    "G_lambda_multipliers": [".999", "1"], "kappa_M": ["7.8", "10.8"],
    "frequencies": ["0", ".25", "1", "2", "3", "10", "1e4", "omega/2",
                    "omega-4", "omega-1", "omega", "omega+1", "omega+4",
                    "2omega", "10omega", "1000omega"],
    "packet_omega": ["10", "100", "1e4", "benchmark_midpoint"],
    "packet_widths": [4, 8, 16, 32], "packet_relative_error": ".02",
    "packet_refinement": ".35", "roundoff_floor": "1e-70",
    "fixtures_are_global_exteriors": False,
}
ROWS, NUMERICS = [], []
R = S.Rational


def check(name, ok, detail=""):
    assert name not in {r["name"] for r in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expr):
    if isinstance(expr, S.MatrixBase):
        check(name, all(S.factor(S.together(x)) == 0 for x in expr))
    else:
        check(name, S.factor(S.together(expr)) == 0)


def latest(stem):
    d = json.loads((ROOT/"evidence"/(stem+".json")).read_text(encoding="utf-8"))
    return d.get("runs", [d])[-1]


def inputs(r, after=False):
    return r.get("input_sha256_after" if after else "input_sha256_before",
                 r.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def paths():
    out = set(OWN)
    for stem, _ in PREDECESSORS:
        out.update([f"notes/{stem}.md", f"evidence/{stem}.json"])
        out.update(inputs(latest(stem)))
    return sorted(out)


def predecessors():
    for stem, count in PREDECESSORS:
        r = latest(stem)
        check(stem+"_receipt", r["status"] == "passed" and len(r["checks"]) == count
              and all(c["passed"] for c in r["checks"]) and inputs(r) == inputs(r, True))
        for p, h in inputs(r).items():
            check(stem+"_input_"+p, sha256((ROOT/p).read_bytes()).hexdigest() == h)


def exact():
    y, a, om = S.symbols("y a omega", positive=True)
    den = (y-a)**2+16*y
    maxi = (a*a+2*a-15)/(16*(a-4))
    ym = a+8+48/(a-7)
    zero("scalar_exact_maximizer", S.diff((1+y)**2/den, y).subs(y, ym))
    zero("scalar_exact_maximum", ((1+y)**2/den).subs(y, ym)-maxi)
    zero("scalar_maximum_completed_square",
         maxi*den-(1+y)**2-(a-7)**2/(16*(a-4))*(y-ym)**2)
    zero("scalar_maximum_omega_form", maxi.subs(a, om**2+4)-(om**2+10+9/om**2)/16)
    # Positive coefficients after substitution omega>=10.
    h = S.symbols("h", nonnegative=True)
    upper_gap = S.cancel(om**2*((om+1)**2-(om**2+10+9/om**2)))
    check("scalar_upper_simplification_positive", all(c > 0 for c in S.Poly(upper_gap.subs(om, 10+h), h).all_coeffs()))
    zero("q_modulus_gap", (y+R(49, 4))*(y+R(1, 4))-(1+y)**2-R(21, 2)*y-R(33, 16))
    zero("q_minus_four_modulus_gap", (y+R(81, 4))*(y+R(1, 4))-(1+y)**2-R(37, 2)*y-R(65, 16))
    zero("scalar_denominator_AMGM", (y+R(41, 4))**2-(y+R(1, 4))*(y+R(81, 4))-100)
    check("scalar_lower_reserve", R(41, 4) < 11)
    disk_gap = (om-11/om)/(1+1/(18*om**2))-(om-12/om)
    check("disk_lower_reserve", S.factor(disk_gap) == (17*om**2+12)/(om*(18*om**2+1)))
    zero("TT_series_denominator", (1+4*y*y)-(2+y)**2-(3*y*y-4*y-3))
    check("TT_log_comparison_at_two", 3*2**2-4*2-3 > 0)
    check("TT_low_interval_budget", 800-400 > 5*2)
    check("TT_mid_interval_budget", 800-8 > 10*2)
    positive = 5*y*y-R(75, 4)-R(19, 4)*y
    check("TT_high_interval_polynomial", all(c > 0 for c in S.Poly(positive.subs(y, 3+h), h).all_coeffs()))
    zero("TT_imaginary_recurrence_bound",
         2-y*(1/(y+R(1, 4))+1/(y+R(9, 4)))
         -(R(1, 4)/(y+R(1, 4))+R(9, 4)/(y+R(9, 4))))
    check("rank_one_factor_upper", R(1066, 27**2) < R(16, 9))
    check("shape_to_scalar_mixing_bound", R(2, 27)/R(2, 3) == R(1, 9))
    s, j, n = S.symbols("s J N", nonzero=True)
    C = S.Matrix([[1, (s+3)/(3*(s+1))], [0, 2]])
    disk = S.Matrix([[1, -R(1, 6)], [0, R(1, 2)]])
    D = S.Matrix([[1, 1/(3*(s+1))], [0, 1]])
    Dr = D.subs(s, 3-s)
    zero("theta_b_normal_transformation", C*disk-D)
    inv = (Dr.T*S.diag(j, n)*D).inv()
    cs, cr = D[0, 1], Dr[0, 1]
    target = S.diag(1/j, 0)+S.Matrix([-cs, 1])*S.Matrix([[-cr, 1]])/n
    zero("full_reflected_inverse", inv-target)
    q = s*(s-3)
    K = om**2+R(9, 4)
    jj = (q-4)*(q+K)/(80*S.pi**2)
    zero("pure_conformal_output_polynomial", jj/(3*(4-s))+(s+1)*(q+K)/(240*S.pi**2))
    alpha, beta, MM, LL, ee = S.symbols("alpha beta M L defect", positive=True)
    zero("Newton_constant_units_invariant", (alpha*MM/beta)**2*(beta*LL/alpha**2)*(beta*ee)-MM**2*LL*ee)
    lo, hi, k10 = R("3.14159"), R("3.14160"), S.Integer(2000)**10
    olo, ohi = R("4.39e17"), R("4.40e17")
    check("benchmark_omega_lower", olo**2 < 60*lo*k10+R(17, 12))
    check("benchmark_omega_upper", 60*hi*k10/R(".999")+R(17, 12) < ohi**2)
    check("benchmark_global_inverse_floor", 20*lo**2*(olo-12/olo) > R("8.6e19"))
    check("benchmark_normal_inverse_ceiling", 20*hi**2*(ohi+1)+640*hi**2 < R("8.7e19"))


def frequencies(om):
    fixed = [mp.mpf(x) for x in PROTOCOL["frequencies"][:7]]
    return fixed+[om/2, om-4, om-1, om, om+1, om+4, 2*om, 10*om, 1000*om]


def symbols(t, gl, kap):
    s = mp.mpf("3.5")+1j*t
    q = s*(s-3)
    K = 60*mp.pi/gl+mp.mpf(11)/3
    J = (q-4)*(q+K)/(80*mp.pi**2)
    H = mp.digamma(s-1)+mp.pi/2*mp.cot(mp.pi*s/2)+kap
    L = 10*mp.pi/gl+mp.mpf(23)/18
    NT = q*((q+2)*H-L)/(480*mp.pi**2)
    return J, NT, H, 1/(3*(s+1)), 1/(3*(4-s))


def matrix_norm(M):
    trace = sum(abs(M[i, j])**2 for i in range(2) for j in range(2))
    determinant_sq = abs(mp.det(M))**2
    return mp.sqrt((trace+mp.sqrt(trace**2-4*determinant_sq))/2)


def multipliers():
    for gmult in PROTOCOL["G_lambda_multipliers"]:
        gl = mp.mpf(gmult)/mp.mpf(2000)**10
        om = mp.sqrt(60*mp.pi/gl+mp.mpf(17)/12)
        cup = 20*mp.pi**2*mp.sqrt(om**2+10+9/om**2)
        clo = 20*mp.pi**2*(om-12/om)
        for kap in PROTOCOL["kappa_M"]:
            for label, t in zip(PROTOCOL["frequencies"], frequencies(om)):
                name = f"multiplier_{gmult}_{kap}_{label}"
                J, NT, H, cs, cr = symbols(t, gl, mp.mpf(kap))
                w, log = (1+t*t)**2, mp.log(2+abs(t))
                sj = 1/(t*t+mp.mpf(".25"))+1/(t*t+mp.mpf("2.25"))
                err = abs(H.imag+t*sj)
                check(name+"_imaginary_identity", err < mp.mpf("1e-70")*(1+abs(H)), mp.nstr(err, 20))
                rc, rt = w/J, w*log/NT
                check(name+"_scalar_upper", abs(rc) <= cup)
                check(name+"_shape_upper", abs(rt) <= 480*mp.pi**2)
                mat = mp.matrix([[rc+cs*cr*rt/log, -cs*rt/log], [-cr*rt, rt]])
                norm = matrix_norm(mat)
                check(name+"_coupled_upper", norm <= cup+640*mp.pi**2)
                if label == "omega":
                    witness = w/(abs(J)*mp.sqrt(1+abs(cr)**2))
                    check(name+"_physical_floor", witness > clo > mp.mpf("8.6e19"))
                    NUMERICS.append({"G_lambda_multiplier": gmult, "kappa_M": kap,
                                     "omega": mp.nstr(om, 35), "normal_norm_at_omega": mp.nstr(norm, 35),
                                     "physical_witness": mp.nstr(witness, 35),
                                     "normal_upper": mp.nstr(cup+640*mp.pi**2, 35)})


def gaussian_packets():
    t, om, width = S.symbols("t omega width", real=True, positive=True)
    s = R(7, 2)+S.I*t
    q = s*(s-3)
    j = S.expand((q-4)*(q+om**2+R(9, 4)))  # 80pi^2 J
    jb = S.expand(-(s+1)*(q+om**2+R(9, 4))/3)
    in_poly = (1+t*t)**4
    out_poly = S.expand(j*S.conjugate(j)+jb*S.conjugate(jb))
    check("packet_output_degree", S.degree(out_poly, t) == 8 and not out_poly.has(S.I))
    average = lambda p: sum(S.diff(p, t, 2*k).subs(t, om)/(S.factorial(k)*(4*width**2)**k)
                            for k in range(5))
    inp, out = average(in_poly), average(out_poly)
    point_in, point_out = in_poly.subs(t, om), out_poly.subs(t, om)
    zero("packet_zero_width_limit", S.limit(inp/out, width, S.oo)-point_in/point_out)
    norm = S.lambdify((om, width), 80*S.pi**2*S.sqrt(inp/out), "mpmath", cse=True)
    point = S.lambdify(om, 80*S.pi**2*S.sqrt(point_in/point_out), "mpmath", cse=True)
    for spec in PROTOCOL["packet_omega"]:
        omega = (mp.sqrt(60*mp.pi*mp.mpf(2000)**10/mp.mpf(".9995")+mp.mpf(17)/12)
                 if spec == "benchmark_midpoint" else mp.mpf(spec))
        target = point(omega)
        errors = []
        for H in PROTOCOL["packet_widths"]:
            got = norm(omega, H)
            err = abs(got/target-1)
            check(f"packet_{spec}_{H}_point_limit", err < mp.mpf(".02"), mp.nstr(err, 20))
            if errors:
                check(f"packet_{spec}_{H}_refinement",
                      err < mp.mpf(".35")*errors[-1] or max(err, errors[-1]) < mp.mpf("1e-70"))
            errors.append(err)
        NUMERICS.append({"packet_omega": spec, "point_ratio": mp.nstr(target, 35),
                         "relative_errors": [mp.nstr(x, 25) for x in errors]})


def main():
    # Keep original stdout UTF8/LF on Windows, before any scientific output.
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    mp.mp.dps = PROTOCOL["digits"]
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before, error = hashes(), None
    print(json.dumps({"event": "start", "story": "S163", "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        exact()
        multipliers()
        gaussian_packets()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(r["passed"] for r in ROWS)
    run = {
        "story": "S163", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "mpmath": mp.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "numerics": NUMERICS,
        "exception": error,
        "scope": "quantitative normal inverse and actual compatible-inverse lower bound",
        "unresolved": ["global inverse upper bound", "nonlinear map", "signed compatibility",
                       "coupled solution", "physical completion"],
    }
    if error:
        run["failed_checker_source_utf8"] = (ROOT/OWN[0]).read_text(encoding="utf-8")
        run["failed_note_source_utf8"] = (ROOT/OWN[1]).read_text(encoding="utf-8")
    dest = ROOT/"evidence"/(STEM+".json")
    data = json.loads(dest.read_text(encoding="utf-8")) if dest.exists() else {"runs": []}
    data["runs"].append(run)
    dest.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": run["status"], "checks": len(ROWS),
                      "unchanged_inputs": before == after}), flush=True)
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
