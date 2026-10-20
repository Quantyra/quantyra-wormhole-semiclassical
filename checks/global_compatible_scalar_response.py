"""S164: global scalar-component bounds, not the full fourth/log inverse."""
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
STEM = "global-compatible-scalar-response"
OWN = ["checks/global_compatible_scalar_response.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("stationary-weighted-tt-inverse", 52),
    ("stationary-scalar-outgoing-inverse", 107),
    ("stationary-scalar-phase-moment", 80),
    ("general-background-newton-operator", 174),
    ("weighted-normal-inverse-scale", 453),
]
PROTOCOL = {
    "digits": 80, "kernel_k": ["1", "10", "1e4"],
    "kernel_pairs": [["1e-8", "2e-8"], [".1", ".5"], [".5", "2"],
                     ["2", "3"], ["10", "40"], ["100", "100"]],
    "identity_tolerance": "1e-45",
    "source_k": ["1", "3", "10"], "source_rho": [".01", ".3", "1", "2"],
    "quadrature_endpoint": "12", "reconstruction_tolerance": "1e-40",
    "incompatible_rho": ["1e-8", "1e-10"], "bolt_slope_tolerance": "1e-10",
    "fixtures_are_physical_preparation_sources": False,
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
              and all(x["passed"] for x in r["checks"]) and inputs(r) == inputs(r, True))
        for p, h in inputs(r).items():
            check(stem+"_input_"+p, sha256((ROOT/p).read_bytes()).hexdigest() == h)


def exact():
    x = S.symbols("x", positive=True)
    v, w = S.Function("v")(x), S.Function("w")(x)
    zero("weighted_scalar_form_identity", S.diff(v, x)*S.diff(w*w*v, x)
         -S.diff(w*v, x)**2+S.diff(w, x)**2*v*v)
    zero("weighted_weak_source_identity", S.diff(w*w*v, x)
         -w*S.diff(w*v, x)-S.diff(w, x)*(w*v))
    t, E = S.symbols("T E")
    r4 = (1/(t-E)-1/(t+4))/(E+4)
    zero("fourth_order_resolvent_identity", r4-1/((t+4)*(t-E)))
    zero("second_graph_resolvent_identity", t*r4-(E/(t-E)+4/(t+4))/(E+4))
    k, sig = S.symbols("k sigma", positive=True)
    flat_terminal = -S.sin(k*(sig-x))/k
    zero("terminal_flat_homogeneous", -S.diff(flat_terminal, x, 2)-k*k*flat_terminal)
    zero("terminal_source_jump_sign", S.diff(flat_terminal, x).subs(x, sig)-1)
    f = S.exp(-x*x)
    aa, bb = -k*k*S.sqrt(x)*f, 2*x**R(3, 2)*f
    source = S.sqrt(x)*(4-4*x*x-k*k)*f
    zero("manufactured_weak_source", aa+S.diff(bb, x)+bb/(2*x)-source)
    vv = S.sqrt(x)*f
    zero("manufactured_Bessel_source", -S.diff(vv, x, 2)-vv/(4*x*x)-k*k*vv-source)
    check("elliptic_weight_coercivity", R(1, 16)-R(1, 1000) > R(1, 20))
    check("source_gradient_weight", 4+R(1, 1000) < 9)
    check("stronger_TT_vertex_weight", 3 >= 2)
    pi_hi = R(22, 7)
    check("stronger_scalar_vertex_g", 201**2+40*pi_hi*R("1.3e-22") < 1000**2)
    check("stronger_scalar_vertex_f", R(".001")**2+5 < 3**2)
    eta, kmin = R("2e-14"), R("4e17")
    check("terminal_kernel_value_budget", 2/(1-2*eta) < R("2.001"))
    check("terminal_kernel_input_derivative_budget",
          3/(1-2*eta)+8/(kmin*(1-2*eta)) < 4)
    check("outgoing_kernel_budget", 1/(1-eta) < R("1.001"))
    check("Bessel_large_input_derivative", 4*R(3, 2) < 9)
    # sqrt(pi/2+1)<1.61 and sqrt(pi/2)/2<.63, hence sum<3.
    check("Bessel_small_input_derivative", pi_hi/2+1 < R("1.61")**2
          and pi_hi/8 < R(".63")**2 and R("1.61")+R(".63") < 3)
    check("weighted_Helmholtz_value_budget", R("2.001")*173 < 347)
    check("weighted_Helmholtz_vector_budget", 4*173 < 700)
    check("fourth_reference_value_budget", 80+347/kmin < 81 and 700+320 < 1100)
    check("fourth_reference_graph_budget", 320/kmin < 1 and 1280/kmin**2 < 1)

    alpha = S.Integer(2000)**-5  # upper bound for sqrt(g)
    check("auxiliary_inverse_budget", 2/(1-R("8e-12")) < R("2.001"))
    yA = R("2.001")*(250+(R(173173, 13)+R(3117, 2))*alpha)/180
    yB = R("2.001")*(519500+1557*14+6*alpha)/180
    check("auxiliary_scalar_budget", yA < 3, str(yA))
    check("auxiliary_vector_budget", yB < 6100, str(yB))
    check("total_force_vector_budget", 1+18300*alpha < R("1.001"))
    check("full_value_scalar_budget", 81*3001+1100*9*alpha < 244000)
    check("full_value_vector_budget", 81*6100000+1100*R("1.001") < 495000000)
    check("full_value_after_E", R(244000, 180) < 1400 and R(495000000, 180) < 2800000)
    check("full_graph_scalar_budget", R(348*3001, 13)+701*9 < 100000)
    check("full_graph_vector_budget", R(348*6100000, 13)*alpha+701*R("1.001") < 1000)
    factor_sq = 80*pi_hi**2*(8*pi_hi/R(".99"))
    check("TT_source_scalar_normalization", factor_sq*201**2 < 30000**2)
    check("TT_source_vector_normalization", factor_sq*R(".001")**2 < R(".2")**2)
    check("trace_source_normalization", 80*pi_hi**2 < 800)
    check("physical_value_scalar", 1400*800 < R("4.4e7"))
    check("physical_value_shape", 1400*30000+2800000*R(".2")*alpha < R("4.4e7"))
    check("physical_graph_scalar", 100000*800 < R("3.2e9"))
    check("physical_graph_shape", 100000*30000+1000*R(".2") < R("3.2e9"))
    check("fourth_root_g_budget", alpha < R("6e-9")**2)
    grad_reserve = R("4e8")-2*5*R("4.4e7")*R("6e-9")
    check("physical_gradient_budget", R("4.4e7")*R("3.2e9") < grad_reserve**2)
    check("benchmark_value_constant", R("4.4e7")*alpha**2 < R("4.4e-26"))
    check("benchmark_graph_constant", R("3.2e9")*alpha < R("1.1e-7"))
    check("benchmark_gradient_constant", R("4e8")**4*alpha**6 < R("8e-17")**4)


def kernel(k, r, s):
    return mp.pi/2*mp.sqrt(r*s)*(mp.bessely(0, k*r)*mp.besselj(0, k*s)
                                -mp.besselj(0, k*r)*mp.bessely(0, k*s))


def kernel_input(k, r, s):
    return -mp.pi*k/2*mp.sqrt(r*s)*(mp.bessely(0, k*r)*mp.besselj(1, k*s)
                                   -mp.besselj(0, k*r)*mp.bessely(1, k*s))


def kernel_diagnostics():
    for ks in PROTOCOL["kernel_k"]:
        k = mp.mpf(ks)
        for a, b in PROTOCOL["kernel_pairs"]:
            r, s = mp.mpf(a)/k, mp.mpf(b)/k
            name = f"kernel_{ks}_{a}_{b}"
            val, der = kernel(k, r, s), kernel_input(k, r, s)
            direct = mp.diff(lambda y: kernel(k, r, y), s)-val/(2*s)
            check(name+"_value", abs(k*val) <= 2, mp.nstr(k*val, 20))
            check(name+"_input_derivative", abs(der) < 3, mp.nstr(der, 20))
            check(name+"_derivative_identity", abs(der-direct) < mp.mpf("1e-45")*(1+abs(der)))
        s = mp.mpf(".7")/k
        check("kernel_"+ks+"_jump", abs(mp.diff(lambda r: kernel(k, r, s), s)-1) < mp.mpf("1e-45"))
        check("kernel_"+ks+"_diagonal", kernel(k, s, s) == 0)
    tail = 1000*mp.exp(-144)
    check("Gaussian_tail_allowance", tail < mp.mpf("1e-59"), mp.nstr(tail, 20))


def integrate(fun, lower):
    cuts = [lower]+[mp.mpf(s) for s in [".1", ".5", "1", "2", "4", "8", "12"] if mp.mpf(s) > lower]
    return mp.quad(fun, cuts)


def source_diagnostics():
    for ks in PROTOCOL["source_k"]:
        k = mp.mpf(ks)
        moment = integrate(lambda s: mp.sqrt(mp.pi*k/2)*s*mp.besselj(0, k*s)
                            *(4-4*s*s-k*k)*mp.exp(-s*s), mp.mpf(0))
        check("manufactured_"+ks+"_zero_moment", abs(moment) < mp.mpf("1e-40"), mp.nstr(moment, 20))
        for rs in PROTOCOL["source_rho"]:
            r = mp.mpf(rs)
            target = mp.sqrt(r)*mp.exp(-r*r)
            direct = integrate(lambda s: kernel(k, r, s)*mp.sqrt(s)*(4-4*s*s-k*k)*mp.exp(-s*s), r)
            weak = integrate(lambda s: (kernel(k, r, s)*(-k*k*mp.sqrt(s))
                                      -kernel_input(k, r, s)*(2*s**mp.mpf("1.5")))*mp.exp(-s*s), r)
            err1, err2 = abs(direct-target), abs(weak-target)
            check(f"manufactured_{ks}_{rs}_direct", err1 < mp.mpf("1e-40")*(1+abs(target)), mp.nstr(err1, 20))
            check(f"manufactured_{ks}_{rs}_weak", err2 < mp.mpf("1e-40")*(1+abs(target)), mp.nstr(err2, 20))
            NUMERICS.append({"k": ks, "rho": rs, "value": mp.nstr(target, 30),
                             "direct_error": mp.nstr(err1, 20), "weak_error": mp.nstr(err2, 20)})
    # Nonzero moment: fast terminal response has a forbidden logarithmic bolt term.
    radii = [mp.mpf(x) for x in PROTOCOL["incompatible_rho"]]
    values = [integrate(lambda s: kernel(1, r, s)*mp.sqrt(s)*mp.exp(-s*s), r)/mp.sqrt(r) for r in radii]
    slope = (values[0]-values[1])/mp.log(radii[0]/radii[1])
    target = mp.exp(-mp.mpf(".25"))/2
    check("incompatible_source_logarithmic_bolt", abs(slope/target-1) < mp.mpf("1e-10"), mp.nstr(slope, 30))
    check("incompatible_source_moment_nonzero", target > mp.mpf(".38"))
    NUMERICS.append({"incompatible_bolt_slope": mp.nstr(slope, 35),
                     "predicted": mp.nstr(target, 35)})


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    mp.mp.dps = PROTOCOL["digits"]
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before, error = hashes(), None
    print(json.dumps({"event": "start", "story": "S164", "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        exact()
        kernel_diagnostics()
        source_diagnostics()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(x["passed"] for x in ROWS)
    run = {
        "story": "S164", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "mpmath": mp.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "numerics": NUMERICS,
        "exception": error,
        "scope": "global fast value/gradient/second-graph bounds for the actual coupled York scalar",
        "unresolved": ["whole-metric fourth/log inverse", "mass response and norm conversions",
                       "nonlinear family", "actual signed source", "physical completion"],
    }
    if error:
        run["failed_checker_source_utf8"] = (ROOT/OWN[0]).read_text(encoding="utf-8")
        run["failed_note_source_utf8"] = (ROOT/OWN[1]).read_text(encoding="utf-8")
    dest = ROOT/"evidence"/(STEM+".json")
    data = json.loads(dest.read_text(encoding="utf-8")) if dest.exists() else {"runs": []}
    data["runs"].append(run)
    dest.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": run["status"], "checks": len(ROWS), "unchanged_inputs": before == after}), flush=True)
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
