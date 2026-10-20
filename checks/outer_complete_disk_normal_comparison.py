"""S166 complete localized disk-to-normal comparison, exact coefficient ledger."""
from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "outer-complete-disk-normal-comparison"
OWN = ["checks/outer_complete_disk_normal_comparison.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("outer-complete-massive-model-comparison", 277),
    ("complete-compact-disk-graph", 217),
    ("weighted-normal-inverse-scale", 453),
]
PROTOCOL = {
    "operator": "complete physical affine disk Hessian including off-shell density contact",
    "input": "theta=S+T/3, b=2T; Hilbert H4/H4log at weight 7/2",
    "output": "f_theta=f_S, f_b=(f_T-f_S/3)/2; Hilbert L2 at weight 7/2",
    "auxiliary_variables": "same quantified half-density equation similarity; norms and inverses <1.01",
    "normal_identification": "all four matrix entries, reflected 3-s transpose and finite contacts",
    "localized_bound": "1e67 epsilon, epsilon<=1e-15, for supported fast input",
    "normal_composition": "<1e-30 at epsilon=1e-120",
    "global_inverse": False,
    "nonlinear_contraction": False,
}
R = S.Rational
ROWS = []


def check(name, ok, detail=""):
    assert name not in {r["name"] for r in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expression):
    expression = S.cancel(S.expand(expression))
    check(name, expression == 0, "" if expression == 0 else expression)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_bytes())
    return data.get("runs", [data])[-1]


def inputs(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def paths():
    out = set(OWN)
    for stem, _ in PREDECESSORS:
        out.add("evidence/"+stem+".json")
        out.update(inputs(latest(stem)))
    return sorted(out)


def predecessors():
    expected = {}
    for stem, count in PREDECESSORS:
        run = latest(stem)
        check(stem+"_receipt", run["status"] == "passed" and len(run["checks"]) == count
              and all(r["passed"] for r in run["checks"]) and inputs(run) == inputs(run, True))
        conflicts = [p for p, h in inputs(run).items() if p in expected and expected[p] != h]
        check(stem+"_shared_inputs", not conflicts, conflicts)
        expected.update(inputs(run))
    for p, h in sorted(expected.items()):
        check("input_"+p, sha256((ROOT/p).read_bytes()).hexdigest() == h)


def radial_and_density_identities():
    y = S.symbols("Y", positive=True)
    h, f = S.Function("h")(y), S.Function("f")(y)
    m = S.exp(h/2)
    dd = lambda a: y*S.diff(a, y)
    hh = dd(h)
    potential = (dd(hh)-hh)/2+hh**2/4
    zero("radial_first_conjugation", y*m*S.diff(f/m, y)-(dd(f)-hh*f/2))
    lap = lambda a: S.diff(a, y, 2)+S.diff(h, y)*S.diff(a, y)
    zero("radial_second_conjugation", y*y*m*lap(f/m)-(dd(dd(f))-dd(f)-potential*f))
    zero("rescaled_curvature", y*y*(-2*S.diff(S.exp(h), y, 2)/S.exp(h))
         +2*(dd(hh)-hh+hh*hh))
    g = S.Function("g")(y)
    zero("weighted_adjoint_D", -y**4*S.diff(y**-3*g, y)-(3*g-dd(g)))
    for order in range(3):
        adj = g
        for _ in range(order):
            adj = 3*adj-dd(adj)
        ordinary = y**4*(-1)**order*S.diff(y**(order-4)*g, y, order)
        # (Y partial_Y)^2 is Y^2 partial_Y^2+Y partial_Y, not Y^2 partial_Y^2.
        if order == 2:
            ordinary -= y**4*S.diff(y**-3*g, y)
        zero("weighted_Euler_order_"+str(order), adj-ordinary)
    alpha, w, c, xi, u = S.symbols("alpha w C xi u", positive=True)
    zero("physical_output_density", (c*u**4/(4*S.pi*xi**2)).subs(u, S.sqrt(xi/c)*y*w)
         -w**4*y**4/(4*S.pi*c))
    b, theta, fs, ft = S.symbols("b theta fS fT")
    zero("physical_b_dual", fs*(theta-b/6)+ft*b/2-(fs*theta+(ft-fs/3)*b/2))
    ell_o, bb = S.symbols("ell_o B", positive=True)
    zero("outer_log_scale", (bb/S.sqrt(c))/S.sqrt(xi/c)-bb/S.sqrt(xi))
    zero("finite_log_constant", -R(31, 30)-R(1, 6)+R(6, 5))


class Jet:
    def __init__(self):
        self.v = {}
        self.reverse = {}
        for field in ["theta", "b", "h", "v"]:
            for order in range(7):
                symbol = S.Symbol(field+str(order))
                self.v[field, order] = symbol
                self.reverse[symbol] = field, order
        self.O, self.A, self.B, self.e, self.z, self.nu, self.chi, self.ell = S.symbols(
            "O Aminus2 Bminus4 e z nu chi ell")

    def f(self, name, order=0):
        return self.v[name, order]

    @lru_cache(maxsize=None)
    def D(self, expression):
        expression = S.sympify(expression)
        out = -S.diff(expression, self.O)
        out += -2*self.A*self.f("v", 1)*S.diff(expression, self.A)
        out += -4*self.B*self.f("v", 1)*S.diff(expression, self.B)
        out += 2*self.e*S.diff(expression, self.e)+2*self.z*S.diff(expression, self.z)
        for symbol in expression.free_symbols:
            if symbol in self.reverse:
                field, order = self.reverse[symbol]
                out += S.diff(expression, symbol)*self.f(field, order+1)
        return S.expand(out)

    def euler(self, q, field):
        out = S.Integer(0)
        for order in range(3):
            term = S.diff(q, self.f(field, order))
            for _ in range(order):
                term = 3*term-self.D(term)
            out += term/2
        return S.expand(out)

    def normal(self, expression):
        sub = {self.A: 1, self.B: 1, self.e: 0, self.z: 0}
        sub.update({self.f(f, k): 0 for f in ["h", "v"] for k in range(7)})
        return S.expand(expression.subs(sub))

    def forms(self):
        theta, b = self.f("theta"), self.f("b")
        ss, tt = theta-b/6, b/2
        hh = self.f("h", 1)
        a1 = lambda f: self.D(f)-hh*f/2
        vh = (self.f("h", 2)-hh)/2+hh*hh/4
        b2 = lambda f: self.D(self.D(f))-self.D(f)-vh*f
        rr = -2*(self.f("h", 2)-hh+hh*hh)
        o1, o2 = -1-self.f("v", 1), -self.f("v", 2)
        uu, vv = o2+(hh-1)*o1, o1*o1
        kk = self.e**2/15-self.nu**2*self.z**2/3
        om = self.O-self.f("v")
        d, f, st, sp = b2(tt), b2(theta), a1(tt), a1(ss)
        ww = o1*sp
        d0 = ss*(-rr/15+self.e/9)+tt*((6*om+1)*rr/45+11*vv/135)-11*ww/135
        l0 = tt*(-rr/15+2*uu/5+7*vv/45)-7*ww/45
        j0 = rr/30-7*(uu+vv)/90-self.e/18
        w0 = -2*rr/15+14*(uu+vv)/45
        z0 = tt**2*(om*rr**2/30+rr*(uu+vv)/15-uu**2/5
                    -7*uu*vv/45-7*vv**2/90+2*kk*om)+ss*tt*(-rr**2/30+2*kk)
        qw = f*f/20-(om+R(1, 6))*d*d/30 \
            -(d*d0+f*l0+j0*sp*sp+w0*tt*ww-7*ww*ww/45+z0)/2
        qc = -(self.A*(4*ss*ss*(rr+6*vv)+48*ss*o1*sp+16*ss*o1*st
                       +12*sp*sp+8*sp*st+8*self.e*(ss+tt)**2)
               +6*self.ell*self.B*(4*ss+2*tt)**2-8*self.chi*self.z**2*tt**2)/4
        ql = -self.nu*self.z*st*st/12
        return [S.expand(q) for q in [qw, qc, ql]]

    def ledger(self, expression):
        da, db = S.symbols("deviation_A deviation_B")
        expression = S.expand(expression.subs({self.A: 1+da, self.B: 1+db}))
        if expression == 0:
            return 0, 0, 0, 0, 0
        small = {da, db, self.e, self.z}
        small.update(self.f(f, k) for f in ["h", "v"] for k in range(7))
        symbols = sorted(expression.free_symbols, key=str)
        weights = {self.nu: 400, self.chi: 2, self.ell: 2}
        total, max_tangent, max_background, max_O = 0, 0, 0, 0
        terms = S.Poly(expression, *symbols).terms()
        for powers, coeff in terms:
            if coeff == 0:
                continue
            weight, tangent_degree, vanishing = abs(coeff), 0, 0
            for symbol, power in zip(symbols, powers):
                if not power:
                    continue
                if symbol in small:
                    vanishing += power
                if symbol == self.O:
                    max_O = max(max_O, power)
                if symbol in self.reverse:
                    field, order = self.reverse[symbol]
                    if field in ["theta", "b"]:
                        tangent_degree += power
                        max_tangent = max(max_tangent, order)
                    else:
                        max_background = max(max_background, order)
                weight *= weights.get(symbol, 1)**power
            assert tangent_degree == 1 and vanishing >= 1, (powers, coeff)
            total += weight
        return total, len(terms), max_tangent, max_background, max_O


def complete_local_operators():
    j = Jet()
    forms = j.forms()
    operators = [[j.euler(q, field) for field in ["theta", "b"]] for q in forms]
    normal = [[j.normal(op) for op in pair] for pair in operators]
    for label, q, pair in zip(["Weyl", "classical", "LLL"], forms, normal):
        for row, field in enumerate(["theta", "b"]):
            zero(label+"_normal_commutes_Euler_"+str(row), pair[row]-j.euler(j.normal(q), field))
    s, g, F = S.symbols("s g F", nonzero=True)
    qs = s*(s-3)
    pp = s*(s-1)*(s-2)*(s-3)
    K0 = 60*S.pi/g+R(11, 3)
    local = [S.expand(normal[0][row]/(4*S.pi**2)
                      +normal[1][row].subs(j.ell, 1+11*g/(360*S.pi))/(4*S.pi*g))
             for row in range(2)]
    actual = S.zeros(2)
    for row in range(2):
        for col, field in enumerate(["theta", "b"]):
            actual[row, col] = sum(S.diff(local[row], j.f(field, n))*s**n for n in range(5))
    massive = (pp*(F+j.O+R(1, 6))-s*(s-1)*(2*s-5))/(480*S.pi**2)
    actual[1, 1] += massive
    direct = S.Matrix([
        [(qs-4)*(qs+K0)/(80*S.pi**2), (s-4)*(qs+K0)/(240*S.pi**2)],
        [-(s+1)*(qs+K0)/(240*S.pi**2), pp*F/(480*S.pi**2)
         -(1/(48*S.pi*g)+R(23, 8640)/S.pi**2)*qs-(qs+K0)/(720*S.pi**2)]])
    J = (qs-4)*(qs+K0)/(80*S.pi**2)
    NT = pp*F/(480*S.pi**2)-(1/(48*S.pi*g)+R(23, 8640)/S.pi**2)*qs
    mat = lambda x: S.Matrix([[1, 1/(3*(x+1))], [0, 1]])
    reflected = mat(3-s).T*S.diag(J, NT)*mat(s)
    for row in range(2):
        for col in range(2):
            zero(f"complete_normal_entry_{row}_{col}", actual[row, col]-direct[row, col])
            zero(f"reflected_normal_entry_{row}_{col}", reflected[row, col]-direct[row, col])
            zero(f"logarithm_cancellation_{row}_{col}", S.diff(actual[row, col], j.O))
    y = S.symbols("Y", positive=True)
    zero("massive_outer_two_derivatives", y**(4-s)*S.diff(
        s*(s-1)*y**(s-2)*(F-S.log(y)), y, 2)
         -(pp*(F-S.log(y))-s*(s-1)*(2*s-5)))
    for label, pair, lim in zip(["Weyl", "classical", "LLL"], operators, normal):
        total = 0
        for row in range(2):
            delta = S.expand(pair[row]-lim[row])
            budget, count, tangent, background, logdegree = j.ledger(delta)
            total += budget
            check(label+"_vanishing_full_ledger_"+str(row),
                  tangent <= 4 and background <= 4 and logdegree <= 1,
                  {"monomials": count, "coefficient_sum": str(budget),
                   "tangent_order": tangent, "background_order": background, "log_degree": logdegree})
        check(label+"_complete_error_budget", total < 10**12, total)
    check("LLL_normal_zero", normal[2] == [0, 0])


def analytic_and_operator_budgets():
    t = S.symbols("t", positive=True)
    series = t/(1-t)
    for _ in range(5):
        series = S.cancel(t*S.diff(series, t))
    check("analytic_positive_series", S.cancel((series-t)/t**2).subs(t, R(1, 2)) == 4326)
    check("analytic_Cauchy_budget", 2*4326*20**2 < 10**7)
    check("analytic_small_ratio", 20*R("1e-15") < R(1, 2))
    check("small_error_eta", 10**7*R("1e-15")**2 < 1)
    check("k_lower", R(".4")**2 < 1/R("6.1"))
    check("k_upper", 1/R("5.9") < R(".5")**2)
    check("inverse_g_guard", R(2000)**10/R(".999") < 2*10**33)
    check("corrected_lambda_ratio", 1+11/(R(360)*3*2000**10) < 2)
    check("combined_local_coefficient", 3*10**12*2*10**33 <= 6*10**45)
    for order in range(5):
        check("weighted_derivative_norm_"+str(order), 2*R(9, 2)**order < 10**3)
    check("logarithm_envelope", 2+1 < 4)
    check("local_operator_error", 6*10**45*10**7*10**3*4 < 10**58)
    yy = S.symbols("y", nonnegative=True)
    # Coefficient positivity proves the three quadratic inequalities on all y>=0.
    sq = [(R(7, 4)-yy)**2+16*yy,
          (-R(9, 4)-yy)**2+16*yy,
          (R(15, 4)-yy)**2+16*yy]
    for k, qsq in enumerate(sq):
        diff = S.Poly(36*(1+yy)**2-qsq, yy)
        check("normal_q_uniform_"+str(k), all(c > 0 for c in diff.all_coeffs()), diff.as_expr())
    check("normal_K0_guard", 60*R(22, 7)*2*10**33+R(11, 3) < 4*10**35)
    check("normal_K0_plus_q", 4*10**35+6 < 5*10**35)
    check("normal_log_budget", 4+R(54, 5)+R(4, 3)+1 < 29*R(2, 3))
    # Use direct matrix entries; all positive reciprocal pi denominators are dropped.
    entries = [6*5*10**35, 2*2*5*10**35, 5*5*10**35,
               36*30+2*6*(2*10**33+1)+2*5*10**35]
    check("normal_forward_bound", sum(entries) < 10**38, sum(entries))
    check("density_error_budget", 10**7*10**38 == 10**45)
    check("complete_localized_comparison", 10**65+R("1.5e58")+10**45 < 10**67)
    check("normal_inverse_composition", R("8.7e19")*10**67*R("1e-120") < R("1e-30"))
    check("collar_contract", R("1e-120") < R("1e-15"))


def main():
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S166 complete outer disk comparison",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        radial_and_density_identities()
        complete_local_operators()
        analytic_and_operator_budgets()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(c["passed"] for c in ROWS)
    run = {
        "story": "S166 complete outer disk comparison", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version, "sympy": S.__version__,
        "protocol": PROTOCOL, "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "complete_localized_disk_normal_comparison": good,
        "S166_complete": False, "goal_complete": False,
        "global_TT_inverse": False, "nonlinear_contraction": False,
        "unresolved": ["mass/fast localization and nonlocal cutoff terms", "annulus and bolt matching",
                       "higher York maps", "full inverse/residual composition", "nonlinear signed compatibility",
                       "exterior/state, charge, return, opening, stability, transport, preparation and creation"],
    }
    if not good:
        run["failed_checker_source_utf8"] = (ROOT/OWN[0]).read_text(encoding="utf-8")
        run["failed_note_source_utf8"] = (ROOT/OWN[1]).read_text(encoding="utf-8")
    dest = ROOT/"evidence"/(STEM+".json")
    data = json.loads(dest.read_bytes()) if dest.exists() else {"runs": []}
    data["runs"].append(run)
    dest.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": run["status"], "checks": len(ROWS),
                      "unchanged_inputs": before == after}), flush=True)
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
