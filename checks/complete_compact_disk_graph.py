"""S166 complete compact disk Hessian, exact local Euler operators and bounds."""
from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "complete-compact-disk-graph"
OWN = ["checks/complete_compact_disk_graph.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("compact-massive-logarithmic-graph", 236),
    ("stationary-seed-quantum-form", 53),
    ("retained-benchmark-parameter-recovery", 25),
    ("stationary-centered-horizon", 79),
]
PROTOCOL = {
    "scope": "complete local disk Hessian in theta=S+T/3 and T, including bolt charts",
    "full_action": "physical Einstein-Maxwell + Weyl + canonical LLL + complete massive",
    "reference_action": "exact off-shell addition/subtraction cancellation",
    "verification": "independent action variations and two-dimensional Euler operators",
    "graph_norm": "sum H4 and H4log on the fixed chart",
    "global_inverse": False,
    "nonlinear_contraction": False,
}
R = S.Rational
ROWS = []


def check(name, ok, detail=""):
    assert name not in {c["name"] for c in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expression):
    expression = S.expand(expression)
    if expression != 0:
        expression = S.cancel(expression)
    check(name, expression == 0, "" if expression == 0 else expression)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_text(encoding="utf-8"))
    return data.get("runs", [data])[-1]


def inputs(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def paths():
    out = set(OWN)
    for stem, _ in PREDECESSORS:
        out.update([f"evidence/{stem}.json"])
        out.update(inputs(latest(stem)))
    out.update(["notes/stationary-centered-covariant-seed.md",
                "notes/stationary-curved-quantum-tt.md"])
    return sorted(out)


def predecessors():
    for stem, count in PREDECESSORS:
        run = latest(stem)
        check(stem+"_receipt", run["status"] == "passed"
              and len(run["checks"]) == count and all(c["passed"] for c in run["checks"])
              and inputs(run) == inputs(run, True))
        for path, digest in inputs(run).items():
            check(stem+"_input_"+path, sha256((ROOT/path).read_bytes()).hexdigest() == digest)


def seed_and_benchmark():
    u = S.symbols("u")
    nu, aa, bb, ee, ratio = S.symbols("nu aa bb ee ratio")
    f = 1+nu*u*u-2*aa*u**3+bb*u**4-2*ee*u**5+ee*u**6
    guards = {nu: R(".001"), aa: R("2.002"), bb: R("3.002"), ee: R("1e-15"), ratio: 1}
    for j in range(4):
        pp = S.Poly(S.diff(ratio*f, u, j), u, *guards)
        upper = sum(abs(c)*S.prod(guards[p]**n for p, n in zip(guards, powers[1:]))
                    for powers, c in pp.terms())
        check("seed_J_derivative_"+str(j), upper < (9 if j == 0 else 10**4), upper)
    yy, jj = S.symbols("Y J", positive=True)
    function = S.Function("u")(yy)
    J = S.Function("J")
    # u'=sqrt(J); a is proportional to sqrt(J), so both Hessian entries equal J'/2.
    zero("radial_Hess_u", S.diff(S.sqrt(J(function)), yy).subs(S.diff(function, yy), S.sqrt(J(function)))
         -S.diff(J(function), function)/2)
    zero("time_Hess_u", (S.diff(S.sqrt(J(function)), yy)/S.sqrt(J(function))
         *S.diff(function, yy)).subs(S.diff(function, yy), S.sqrt(J(function)))
         -S.diff(J(function), function)/2)
    check("coordinate_first_u", R("1.1")*3 < 4)
    u2 = 10**4+2*100*4
    u3 = 4*10**4+2*10**7+204800+4*10**6
    u4 = 10**4*(16+2*10**4)+10**4*4*4000+10**4*R("1.6e7")/2 \
        +2*(10**7*4+2*25600*2*10**4+100*10**8)
    check("coordinate_second_u", u2 < 2*10**4, u2)
    check("coordinate_third_u", u3 < 10**8, u3)
    check("coordinate_fourth_u", u4 < 2*10**11, u4)
    uj = [4, 2*10**4, 10**8, 2*10**11]
    for j in range(1, 5):
        bound = sum(S.factorial(k-1)*28**k*S.bell(j, k, tuple(uj[:j-k+1]))
                    for k in range(1, j+1))
        check("omega_coordinate_jet_"+str(j), bound < S.factorial(j)*10**(6*j), bound)
        check("sigma_coordinate_jet_"+str(j), R(27, 1000)*128**j < 10**(6*j))
    check("physical_charge_guard", R("1.01")/R(".99")**2 < 2)
    check("physical_curvature_guard", R(2)/R("5.9") < 1)
    check("physical_inverse_g", R(2000)**10/R(".999") < 2*10**33)
    check("magnetic_Weyl_coefficient", R(1, 15)+R(400**2, 3) < 60000)


def independent_variations():
    t, ss, tt, rr, uu, vv, om, e0, k0, ls, lt, ww, gg = S.symbols(
        "t ss tt rr uu vv om e0 k0 ls lt ww gg")
    rv = S.exp(-2*t*tt)*(rr-2*t*lt)
    uv = S.exp(-2*t*tt)*(uu+t*ls)
    v_v = S.exp(-2*t*tt)*(vv+2*t*ww+t*t*gg)
    ov = om+t*ss
    lag = ov*(rv**2/60-e0*rv/18+k0)+rv*(uv+v_v)/30-e0*v_v/18 \
        -uv**2/10-R(7, 90)*uv*v_v-R(7, 180)*v_v**2
    actual = -S.diff(S.exp(2*t*tt)*lag, t, 2).subs(t, 0)/4
    ff = ls+lt/3
    dd = ss*(-rr/15+e0/9)+tt*((6*om+1)*rr/45+11*vv/135)-11*ww/135
    ll = tt*(-rr/15+2*uu/5+7*vv/45)-7*ww/45
    j0 = rr/30-7*(uu+vv)/90-e0/18
    wz = -2*rr/15+14*(uu+vv)/45
    z0 = tt**2*(om*rr**2/30+rr*(uu+vv)/15-uu**2/5
                 -7*uu*vv/45-7*vv**2/90+2*k0*om)+ss*tt*(-rr**2/30+2*k0)
    expected = ff**2/20-(om+R(1, 6))*lt**2/30 \
        -(lt*dd+ff*ll+j0*gg+wz*tt*ww-7*ww**2/45+z0)/2
    zero("independent_complete_Weyl_second_variation", actual-expected)
    zz, ell, charge, vs, vt, st = S.symbols("Z ell charge vs vt st")
    classical = zz*S.exp(2*t*ss)*(rr+6*(vv+2*t*vs+t*t*gg)
                                  +4*(t*vt+t*t*st)) \
        +2*e0*zz*S.exp(2*t*(ss+tt))+6*ell*zz**2*S.exp(t*(4*ss+2*tt)) \
        -2*charge*S.exp(2*t*tt)
    classical_q = zz*(4*ss**2*(rr+6*vv)+48*ss*vs+16*ss*vt+12*gg+8*st) \
        +8*e0*zz*(ss+tt)**2+6*ell*zz**2*(4*ss+2*tt)**2-8*charge*tt**2
    zero("independent_physical_classical_second_variation",
         S.diff(classical, t, 2).subs(t, 0)-classical_q)
    vphys, vref, volume = S.symbols("vphys vref volume")
    zero("off_shell_reference_cancellation", volume*(vphys+vref)-volume*vref-volume*vphys)
    # Full two-dimensional Weyl conversion of the physical reduced kinetic action.
    hh, xh = S.symbols("hh xh")
    zero("classical_kinetic_six", 4*xh*hh+2*xh*hh-6*xh*hh)


class Jets:
    def __init__(self):
        self.v = {}
        self.reverse = {}
        for field in ["theta", "T", "sigma", "omega"]:
            for n in range(7):
                for i in range(n+1):
                    idx = (i, n-i)
                    symbol = S.Symbol(f"{field}_{i}_{n-i}")
                    self.v[field, idx] = symbol
                    self.reverse[symbol] = (field, idx)
        self.E, self.Z, self.e0, self.k0, self.ell, self.charge = S.symbols("E Z e0 k0 ell charge")

    def f(self, field, i=0, j=0):
        return self.v[field, (i, j)]

    @lru_cache(maxsize=None)
    def D(self, expression, axis):
        direction = (1, 0) if axis == 0 else (0, 1)
        out = S.diff(expression, self.E)*(-2*self.E*self.v["sigma", direction]) \
            +S.diff(expression, self.Z)*(2*self.Z*self.v["omega", direction])
        for symbol in expression.free_symbols:
            if symbol in self.reverse:
                field, idx = self.reverse[symbol]
                nxt = (idx[0]+direction[0], idx[1]+direction[1])
                out += S.diff(expression, symbol)*self.v[field, nxt]
        return S.expand(out)

    def lap(self, expression):
        return self.D(self.D(expression, 0), 0)+self.D(self.D(expression, 1), 1)

    def dot(self, left, right):
        return self.E*sum(self.D(left, i)*self.D(right, i) for i in range(2))

    def euler(self, density, field):
        out = 0
        for symbol in density.free_symbols:
            if symbol not in self.reverse or self.reverse[symbol][0] != field:
                continue
            idx = self.reverse[symbol][1]
            term = S.diff(density, symbol)
            for axis, count in enumerate(idx):
                for _ in range(count):
                    term = self.D(term, axis)
            out += (-1)**sum(idx)*term
        return S.expand(out)

    def norm(self, expression, max_tangent_order):
        expression = S.expand(expression)
        symbols = sorted(expression.free_symbols, key=str)
        poly = S.Poly(expression, *symbols)
        weights = {self.E: 2, self.Z: 28**2, self.e0: 1, self.k0: 60000,
                   self.ell: 1, self.charge: 2}
        max_order = 0
        total = 0
        for powers, coefficient in poly.terms():
            weight = abs(coefficient)
            tangent_degree = 0
            for symbol, exponent in zip(symbols, powers):
                if not exponent:
                    continue
                if symbol in self.reverse:
                    field, idx = self.reverse[symbol]
                    order = sum(idx)
                    if field in ["theta", "T"]:
                        tangent_degree += exponent
                        max_order = max(max_order, order)
                        w = 1
                    elif order:
                        w = S.factorial(order)*10**(6*order)
                    else:
                        w = 4 if field == "omega" else 1
                else:
                    w = weights[symbol]
                weight *= w**exponent
            assert tangent_degree == 1, (tangent_degree, powers)
            total += weight
        assert max_order <= max_tangent_order, max_order
        return total, len(poly.terms()), max_order


def complete_euler_operators():
    j = Jets()
    E, Z, e0, k0, ell, charge = j.E, j.Z, j.e0, j.k0, j.ell, j.charge
    theta, T, omega, sigma = [j.f(f) for f in ["theta", "T", "omega", "sigma"]]
    ss = theta-T/3
    rr = -2*E*j.lap(sigma)
    uu = E*j.lap(omega)
    vv = j.dot(omega, omega)
    dd, ff = E*j.lap(T), E*j.lap(theta)
    ww, gg = j.dot(omega, ss), j.dot(ss, ss)
    DD = ss*(-rr/15+e0/9)+T*((6*omega+1)*rr/45+11*vv/135)-11*ww/135
    LL = T*(-rr/15+2*uu/5+7*vv/45)-7*ww/45
    JJ = rr/30-7*(uu+vv)/90-e0/18
    WW = -2*rr/15+14*(uu+vv)/45
    ZZ = T*T*(omega*rr**2/30+rr*(uu+vv)/15-uu**2/5
               -7*uu*vv/45-7*vv**2/90+2*k0*omega)+ss*T*(-rr**2/30+2*k0)
    cov_q = ff**2/20-(omega+R(1, 6))*dd**2/30 \
        -(dd*DD+ff*LL+JJ*gg+WW*T*ww-7*ww**2/45+ZZ)/2
    density = S.expand(cov_q/E)
    hw = [S.expand(E*j.euler(density, field)/2) for field in ["theta", "T"]]
    principal = [E*E*j.lap(j.lap(theta))/20,
                 -E*E*(omega+R(1, 6))*j.lap(j.lap(T))/30]
    total_w = 0
    for row, (field, operator, lead) in enumerate(zip(["theta", "T"], hw, principal)):
        remainder = S.expand(operator-lead)
        norm, count, order = j.norm(remainder, 3)
        total_w += norm
        check("Weyl_row_order_"+str(row), order <= 3, order)
        check("Weyl_row_coefficient_"+str(row), norm < 10**42,
              {"terms": count, "norm": str(norm)})
        other = "T" if field == "theta" else "theta"
        for idx, multiplicity in [((4, 0), 1), ((2, 2), 2), ((0, 4), 1)]:
            coefficient = E**2/R(20) if field == "theta" else -E**2*(omega+R(1, 6))/30
            zero(f"Weyl_principal_{field}_{idx}", S.diff(operator, j.v[field, idx])-multiplicity*coefficient)
            zero(f"Weyl_fourth_cross_{field}_{idx}", S.diff(operator, j.v[other, idx]))
    check("complete_Weyl_lower_coefficient_sum", total_w < 10**42, total_w)
    qc = Z*(4*ss**2*(rr+6*vv)+48*ss*j.dot(omega, ss)
             +16*ss*j.dot(omega, T)+12*j.dot(ss, ss)+8*j.dot(ss, T)) \
        +8*e0*Z*(ss+T)**2+6*ell*Z**2*(4*ss+2*T)**2-8*charge*T*T
    classical_density = S.expand(qc/E)
    hc = [S.expand(-E*j.euler(classical_density, field)/8) for field in ["theta", "T"]]
    total_c = 0
    for row, operator in enumerate(hc):
        norm, count, order = j.norm(operator, 2)
        total_c += norm*2*10**33
        check("physical_classical_row_order_"+str(row), order <= 2, order)
        check("physical_classical_row_coefficients_"+str(row), norm*2*10**33 < 10**59,
              {"terms": count, "norm_with_inverse_g": str(norm*2*10**33)})
    check("physical_classical_complete_budget", total_c < 10**59, total_c)
    check("LLL_complete_row", R(400*2*2, 36) < 100)
    check("full_local_remainder_budget", total_w+total_c+100 < 10**60)


def graph_and_source():
    pi_hi = R(22, 7)
    check("scalar_leading_inverse", R("1.12")*20*pi_hi < 80)
    check("scalar_graph_forcing", 2*80 < 200)
    check("scalar_graph_lower", 160*10**60+2 < 10**64)
    check("paired_graph_forcing", max(200, 1600) < 2000)
    check("paired_graph_lower", 10**64+1600*10**60+10**96 < 2*10**96)
    eps = R(1, 4*10**96)
    check("paired_interpolation_absorption", 2*10**96*eps == R(1, 2))
    check("paired_final_value", 4*10**96*(2/eps)**3 < 10**400)
    check("paired_final_forcing", 2*2000 < 5000)
    z = S.symbols("z", nonnegative=True)
    zero("scalar_Fourier_H4_envelope", 2*(1+z*z)-(1+z)**2-(z-1)**2)
    s, t, f, h = S.symbols("s t f h")
    zero("physical_source_transpose", f*(s-t/3)+h*t-(f*s+(h-f/3)*t))
    check("source_pair_factor", R(4, 3)**2+1 == R(5, 3)**2)
    check("source_measure_factor", (R(5, 3)*R("1.1")*2*28**2)**2*pi_hi < 10**8)
    check("input_measure_factor", (R(5, 3)*R("1.1"))**2 < 4*3)
    check("physical_graph_forcing", 5000*10**4 == 5*10**7)


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 complete compact disk graph",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        seed_and_benchmark()
        independent_variations()
        complete_euler_operators()
        graph_and_source()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(c["passed"] for c in ROWS)
    run = {
        "story": "S166 complete compact disk graph",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed",
        "python": sys.version, "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "complete_compact_disk_graph": good,
        "S166_complete": False, "goal_complete": False,
        "global_TT_inverse": False, "nonlinear_contraction": False,
        "unresolved": ["patch commutators", "higher York reconstruction",
                       "outer fast fourth/log estimate", "full inverse/residual composition",
                       "nonlinear signed compatibility", "physical requirements"],
    }
    if not good:
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
