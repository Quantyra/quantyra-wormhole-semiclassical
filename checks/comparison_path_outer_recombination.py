"""S158: local tensor, correlated cancellation and fixed-chart response."""
from __future__ import annotations
from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = "comparison-path-outer-recombination"
INPUTS = [
    "checks/comparison_path_outer_recombination.py",
    "notes/comparison-path-outer-recombination.md",
    "notes/comparison-path-outer-regular-response.md",
    "evidence/comparison-path-outer-regular-response.json",
    "notes/fixed-chart-global-constraint-transfer.md",
    "evidence/fixed-chart-global-constraint-transfer.json",
    "notes/stationary-completed-metric-response-norm.md",
    "notes/stationary-centered-first-order-geometry.md",
    "notes/stationary-centered-local-response.md",
    "notes/spherical-anomaly-sources.md",
    "notes/stationary-weyl-radial-sources.md",
]
PROTOCOL = {
    "digits": 70, "y": [".002", ".02", ".06"], "parameters": ["0", ".4", "1"],
    "gauss_nodes": [32, 64], "contour_nodes": [64, 128], "contour_radius": ".01",
    "recombination_tolerance": "1e-35", "quadrature_tolerance": "1e-30",
    "real_steps": [".001", ".0005"], "second_tolerance": "1e-6",
    "refinement_ratio": ".35", "refinement_floor": "1e-35",
    "endpoint_differences": "second-order one-sided within [0,1]",
    "normalization": "eta_G/pi", "fixture_is_a_physical_source": False,
    "distance_gauss_panels": ["0", ".5", "1"],
}
ROWS, NUMERICS = [], []
RAT = s.Rational


def check(name, passed, detail=""):
    assert name not in {row["name"] for row in ROWS}
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not passed:
        raise AssertionError(name)


def predecessors():
    for stem, count in [("comparison-path-outer-regular-response", 202),
                        ("fixed-chart-global-constraint-transfer", 99)]:
        run = json.loads((ROOT/"evidence"/(stem+".json")).read_text())["runs"][-1]
        check(stem+"_receipt", run["status"] == "passed" and len(run["checks"]) == count
              and all(x["passed"] for x in run["checks"])
              and len({x["name"] for x in run["checks"]}) == count)
        check(stem+"_current_inputs", run["unchanged_inputs"]
              and run["input_sha256_before"] == run["input_sha256_after"]
              and all(sha256((ROOT/p).read_bytes()).hexdigest() == h
                      for p, h in run["input_sha256_after"].items()))


def tensor_polynomials():
    w = s.symbols("w0:5")
    p = s.symbols("p0:5")
    e, qm = s.symbols("e qm")
    def D(f):
        return s.expand(sum(s.diff(f, xs[j])*xs[j+1]
                            for xs in (w, p) for j in range(4)))
    rr, uu, vv, pp = s.symbols("rr uu vv pp")
    aa, cc, bb = RAT(11, 360), RAT(1, 20), RAT(1, 30)
    dd = bb-2*aa/3
    lag = (pp*(cc*(rr+2*e)**2/3-4*aa*e*rr-qm/3)
           +bb*rr*uu-3*bb*uu**2+(bb*rr+(2*bb-4*aa)*e)*vv
           -6*dd*uu*vv-3*dd*vv**2)
    P, Q, Z = [s.diff(lag, z) for z in (rr, uu, vv)]
    R, U, V = -2*(w[2]+w[1]**2), p[2]+w[1]*p[1], p[1]**2
    sub = {rr: R, uu: U, vv: V, pp: p[0]}
    lag, P, Q, Z = [s.expand(f.subs(sub)) for f in (lag, P, Q, Z)]
    mt = s.expand(-lag/2+R*P/2+D(D(P))+(D(Q)*p[1]+Q*U)/2)
    my = s.expand(-lag/2+R*P/2+w[1]*D(P)+Z*V-D(Q)*p[1]/2+Q*U/2)
    # Direct one-dimensional Euler variation of a L under w and radial lapse.
    ew = lag+(D(4*w[1]*P-Q*p[1])+w[1]*(4*w[1]*P-Q*p[1])) \
        -2*(D(D(P))+2*w[1]*D(P)+(w[2]+w[1]**2)*P)
    en = lag-2*R*P-2*Q*U-2*Z*V-D(2*w[1]*P-p[1]*Q)-w[1]*(2*w[1]*P-p[1]*Q)
    check("independent_orbit_Euler_variation", s.expand(ew+2*mt) == 0)
    check("independent_radial_lapse_variation", s.expand(en+2*my) == 0)
    variables = (*w[1:], *p, e, qm)
    for label, expr in [("density", mt), ("null", mt-my)]:
        poly = s.Poly(expr, *variables)
        order = [1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 0]
        check(label+"_derivative_order_four", all(sum(a*b for a, b in zip(ex, order)) <= 4
                                                for ex, co in poly.terms()))
        check(label+"_at_most_one_phi", all(ex[4] <= 1 for ex, co in poly.terms()))
        norm = sum(abs(co)*max(1, sum(ex[:9])) for ex, co in poly.terms())
        check(label+"_polynomial_lipschitz_count", norm < 1000, str(norm))
        NUMERICS.append({"polynomial": label, "monomials": len(poly.terms()), "coefficient_degree_sum": str(norm)})
    return w, p, e, qm, mt, s.expand(mt-my)


def singular_and_seed_checks(data):
    w, p, e, qm, mt, mn = data
    y, c, f, d = s.symbols("y c f d", positive=True)
    basew = e*y**2/2+c*y**3
    basep = -s.log(d*y)-e*y**2/6+f*y**3
    subs = {x[j]: s.diff(base, y, j) for x, base in [(w, basew), (p, basep)] for j in range(5)}
    exprs = [s.expand(expr.subs(subs)) for expr in (mt, mn)]
    expected = [{-4: RAT(11, 120), -3: 0, -2: -11*e/180, -1: (12*c+11*f)/30},
                {-4: 0, -3: 0, -2: 0, -1: (23*c+44*f)/30}]
    for label, expr, values in zip(["density", "null"], exprs, expected):
        for k, val in values.items():
            check(f"{label}_Laurent_{k}", s.simplify(expr.coeff(y, k)-val) == 0)
    check("correlated_density_mass", s.simplify(expected[0][-1].subs(f, -c/4)-37*c/120) == 0)
    check("correlated_null_mass", s.simplify(expected[1][-1].subs(f, -c/4)-2*c/5) == 0)
    mismatch = s.symbols("mismatch")
    check("wrong_radius_cubic_detected",
          s.expand(expected[1][-1].subs(f, -c/4+mismatch)-2*c/5) == 22*mismatch/15)
    vac = RAT(11, 120)/y**4*s.exp(-2*e*y*y/3-c*y**3)
    principal_vac = s.series(vac, y, 0, 1).removeO()
    check("vacuum_quartic", principal_vac.coeff(y, -4) == RAT(11, 120))
    check("vacuum_quadratic", principal_vac.coeff(y, -2) == -11*e/180)
    check("vacuum_plus_mass_cubic", s.simplify(principal_vac.coeff(y, -1)+2*c/5-37*c/120) == 0)
    # Polynomial reference remainder: remove every negative power and count
    # coefficients with |c|<=2, e<=1, |Y|<=1, |log(dY)|<=2L.
    for label, expr in zip(["density", "null"], exprs):
        expr = s.expand(expr.subs(f, -c/4))
        regular = s.expand(expr-sum(expr.coeff(y, k)*y**k for k in range(-4, 0)))
        logd, logy = s.symbols("logd logy")
        regular = s.expand_log(regular, force=True).subs({s.log(d): logd, s.log(y): logy})
        norm = sum(abs(co)*2**ex[1]*2**(ex[4]+ex[5])
                   for ex, co in s.Poly(regular, y, c, e, qm, logd, logy).terms())
        check(label+"_reference_remainder_majorant", norm < 10**8, str(norm))
    A, q0, E5 = s.symbols("A q0 E5")
    zs = d*y
    for k in range(2, 7):
        ak = s.symbols("a"+str(k))
        ztrial = zs+ak*y**k
        eq = s.expand(s.diff(ztrial, y)**2-d*d-e*ztrial**2+2*A*ztrial**3-q0*ztrial**4-E5*ztrial**5)
        zs += s.solve(eq.coeff(y, k-1), ak)[0]*y**k
    ws = s.series(s.log(s.diff(zs, y)/d), y, 0, 6).removeO().expand()
    ps = -s.series(s.log(zs/(d*y)), y, 0, 6).removeO().expand()
    check("seed_cubic_correlation", s.simplify(ps.coeff(y, 3)+ws.coeff(y, 3)/4) == 0)
    for label, expr in [("w", ws), ("phi", ps)]:
        for k in [4, 5]:
            coef = s.Poly(expr.coeff(y, k), d, e, A, q0, E5)
            bound = sum(abs(co)*RAT(".42")**ex[0]*RAT(1, 6000)**ex[1]*
                        RAT(".334")**ex[2]*RAT(".501")**ex[3]*RAT("1e-20")**ex[4]
                        for ex, co in coef.terms())
            check(f"seed_{label}_coefficient_{k}", bound < 1, str(bound))
    # Inverse jets follow by differentiating F(G)=identity.
    F1, F2, F3, F4 = s.symbols("F1:5")
    fs = [F1, F2, F3, F4]
    inv = [1/F1, -F2/F1**3, 3*F2**2/F1**5-F3/F1**4,
           -15*F2**3/F1**7+10*F2*F3/F1**6-F4/F1**5]
    for j in range(3):
        derivative = sum(s.diff(inv[j], fs[k])*fs[k+1] for k in range(3))/F1
        check(f"inverse_jet_{j+2}", s.simplify(derivative-inv[j+1]) == 0)
    check("fourth_composition_coefficient_reserve", 1+6+3+4*(3+1)+(15+10+1) == 52)
    Y, Y1, Y2, cc, mm = s.symbols("Y Y1 Y2 cc mm")
    deriv = lambda h: s.diff(h, Y)*Y1+s.diff(h, Y1)*Y2+s.diff(h, cc)*mm
    check("changing_distance_second_formula",
          s.simplify(deriv(deriv(cc/Y))-(-2*mm*Y1/Y**2-cc*Y2/Y**2+2*cc*Y1**2/Y**3)) == 0)


def budget_checks():
    y = mp.mpf(".081")
    L = 1-mp.log(y)
    def uj(j):
        cubic = mp.mpf("1e-11")*y**3*(mp.factorial(3)/mp.factorial(3-j) if j <= 3 else 0)
        return cubic+mp.mpf("2e-8")*mp.factorial(j)*100**j*y**4*L
    check("complex_scalar_argument", 33*mp.mpf("1.2e-8") < mp.mpf("4e-7"))
    expbound = mp.exp(33*mp.mpf("1.2e-8"))
    guards = [expbound-1, expbound*33*uj(1),
              expbound*(33*uj(2)+(33*uj(1))**2),
              expbound*(33*uj(3)+3*33**2*uj(1)*uj(2)+(33*uj(1))**3)]
    for j, value, bound in zip(range(1, 5), guards, ["1e-6", "2e-8", "3e-6", ".001"]):
        check(f"proper_map_scaled_guard_{j}", value < mp.mpf(bound), mp.nstr(value, 30))
    # Ordinary analytic seed tail; closed generating function, no sampled tail.
    x = s.symbols("x")
    for j in range(5):
        total = x**j/(1-x)**(j+1)
        tail = (total-sum(s.binomial(k, j)*x**k for k in range(6)))/x**4
        evaluate = s.lambdify(x, tail, "mpmath")
        val = (sum(mp.binomial(k, j)*y**(k-4) for k in [4, 5])
               +mp.mpf(".00101")/mp.mpf(".15")**4*evaluate(y/mp.mpf(".15")))/(100**j*L)
        check(f"analytic_seed_remainder_{j}", val < mp.mpf(".7"), mp.nstr(val, 30))
    # Analytic multiplier derivative sum <1.021 in the retained seed proof.
    # Product of five normalized jets costs <=5 at order four. For nu=tau V,
    # V-nu = z nu/tau is a retained analytic-multiplier product.
    check("finite_log_radius_linear_remainder", 33*mp.mpf("2e-8") < mp.mpf("1e-6"))
    quadratic = (mp.mpf("1e-11")**2*y*y/L
                 +mp.mpf("4e-19")*y**3+mp.mpf("4e-16")*y**4*L)
    check("finite_log_radius_quadratic_reserve",
          10**6*33**2*quadratic < mp.mpf("1e-7"), mp.nstr(10**6*33**2*quadratic, 30))
    check("quadratic_jet_count_reserve", 2*5**6*RAT("1.05")**8 < 10**6)
    check("finite_radius_remainder_reserve", mp.mpf("1e-6")+mp.mpf("1e-7") < mp.mpf("1e-5"))
    check("distance_log_remainder_reserve",
          33*mp.mpf("2e-8")*mp.mpf(".3")+mp.mpf("1e-7") < mp.mpf("1e-5"))
    check("inverse_composition_reserve", 52*mp.mpf("1.01") < 100)
    check("actual_reference_value_difference", y**4*L < mp.mpf(".001"))
    # Direct jet-polynomial bound, stronger than the old 1e20 reserve.
    # The qmag term is separated; its order-zero difference gains y^4.
    eta = mp.mpf("8.5e-33")
    polyvariation = 1000*100*24*100**4*2
    local = eta/mp.pi*(polyvariation+160000*100*y**4)
    check("direct_complex_tensor_variation", local < mp.mpf("4e-10"), mp.nstr(local, 30))
    vacuum = eta/mp.pi*mp.mpf(11)/120*4*100*mp.exp(mp.mpf(".01"))
    mass = eta/mp.pi*mp.mpf(".8")*100*y**3*mp.exp(mp.mpf(".01"))
    check("explicit_vacuum_and_mass_remainders", local+vacuum+mass < mp.mpf("4e-10"))
    check("reference_normalized_majorant", eta/mp.pi*10**8*(1+160000) < mp.mpf("1e-15"))
    check("complex_ward_ratio", mp.exp(4*mp.mpf(".00056")+4*33*mp.mpf("1.2e-8")) < mp.mpf("1.01"))
    check("complex_ward_gradient", 1+2*33*uj(1)/y**2 < 3)
    check("cubic_pressure_and_J_budget", 10*eta*(1+100*y**3) < mp.mpf("1e-28"))
    check("joint_complex_recombination", mp.mpf("4e-10")+mp.mpf("1e-15")+mp.mpf("1e-28") < mp.mpf("1e-9"))
    check("first_Cauchy_budget", RAT("1e-9")/32 < RAT("3.2e-11"))
    check("second_Cauchy_budget", 2*RAT("1e-9")/32**2 < RAT("2e-12"))
    dg = RAT("3e-7")
    n0, n1, n2 = RAT("1e-9"), RAT("3.2e-11"), RAT("2e-12")
    m2 = RAT("1.003")*n2+RAT("4.02")*dg*n1+RAT("5.54")*dg**2*n0
    a2 = RAT("1.006")*n2+RAT("10.08")*dg*n1+RAT("26.5")*dg**2*n0
    outs = [RAT(".00016")*m2/2, RAT(".835")*m2/2, RAT(".00003")*a2/2]
    bounds = [RAT("1.7e-16"), RAT("8.4e-13"), RAT("3.1e-17")]
    old = [RAT("1.4e-14"), RAT("7.3e-11"), RAT("9e-16")]
    full = [RAT("1.5e-14"), RAT("7.4e-11"), RAT("1e-15")]
    for tag, val, b, prior, joint in zip(["mass", "density", "lapse"], outs, bounds, old, full):
        check("pair_second_"+tag, val < b, str(val))
        check("combined_second_"+tag, val+prior < joint, str(val+prior))
    NUMERICS.append({"pair_second_exact": [str(v) for v in outs],
                     "goal_complete": False, "all_source_sectors_complete": False})


class Jet:
    """Three Taylor coefficients, with independent product/reciprocal algebra."""
    def __init__(self, a, b=0, c=0):
        self.v = (a, b, c)
    @staticmethod
    def lift(x):
        return x if isinstance(x, Jet) else Jet(x)
    def __add__(self, other):
        o = self.lift(other)
        return Jet(*(a+b for a, b in zip(self.v, o.v)))
    __radd__ = __add__
    def __neg__(self):
        return Jet(*(-a for a in self.v))
    def __sub__(self, other):
        return self+-self.lift(other)
    def __rsub__(self, other):
        return self.lift(other)+-self
    def __mul__(self, other):
        o = self.lift(other)
        return Jet(*(sum(self.v[k]*o.v[j-k] for k in range(j+1)) for j in range(3)))
    __rmul__ = __mul__
    def __truediv__(self, other):
        a, b, c = self.lift(other).v
        return self*Jet(1/a, -b/a**2, b*b/a**3-c/a**2)
    def __rtruediv__(self, other):
        return self.lift(other)/self
    def __pow__(self, n):
        if n < 0:
            return 1/self**(-n)
        result = Jet(1)
        for _ in range(n):
            result *= self
        return result
    def exp(self):
        a, b, c = self.v
        ea = mp.exp(a)
        return Jet(ea, ea*b, ea*(c+b*b/2))


def ex(x):
    return x.exp() if isinstance(x, Jet) else mp.exp(x)


def fixtures(data):
    w, p, e, qm, mt, mn = data
    y, lam = s.symbols("y lam", positive=True)
    ee, dd, c0, mm, vv, nn = map(RAT, [".001", ".41", ".02", ".001", ".002", ".003"])
    direction = mm*y**3+vv*y**4*(1-s.log(y))
    wf = ee*y*y/2+c0*y**3+lam*direction
    pf = -s.log(dd*y)-ee*y*y/6-c0*y**3/4-lam*mm*y**3/2+lam*nn*y**4*(1-s.log(y))
    D = lambda f: s.exp(-lam*direction)*s.diff(f, y)
    jets = []
    for f in (wf, pf):
        local = [f]
        for _ in range(4):
            local.append(s.expand(D(local[-1])))
        jets.extend(local)
    jetfun = s.lambdify((y, lam), jets, "mpmath", cse=True)
    jetder = s.lambdify((y, lam), [[s.diff(f, lam, j)/s.factorial(j) for j in range(3)] for f in jets],
                       "mpmath", cse=True)
    tensor = s.lambdify((*w, *p), [mt.subs({e: ee, qm: 4}), mn.subs({e: ee, qm: 4})], "mpmath", cse=True)
    nodes = {n: mp.gauss_quadrature(n, "legendre") for n in PROTOCOL["gauss_nodes"]}
    top = mp.mpf(".08")
    def u(v):
        return mp.mpf(".001")*v**3+mp.mpf(".002")*v**4*(1-mp.log(v))
    def up(v):
        return mp.mpf(".003")*v*v+mp.mpf(".002")*v**3*(3-4*mp.log(v))
    def wbase(v):
        return mp.mpf(".001")*v*v/2+mp.mpf(".02")*v**3
    def wbprime(v):
        return mp.mpf(".001")*v+mp.mpf(".06")*v*v
    @lru_cache(maxsize=None)
    def moments(v):
        L = 1-mp.log(v)
        result = []
        for k in range(11):
            ans = mp.mpf(0)
            for j in range(k+1):
                power = 3*k+j
                integral = v**(power+1)*sum(mp.binomial(j, h)*L**(j-h)*mp.factorial(h)/
                                           (power+1)**(h+1) for h in range(j+1))
                ans += mp.binomial(k, j)*mp.mpf(".001")**(k-j)*mp.mpf(".002")**j*integral
            result.append(ans/mp.factorial(k))
        return result
    def distance(v, q):
        result = 0
        for coef in reversed(moments(v)):
            result = result*q+coef
        return result
    check("distance_series_uniform_tail",
          top*mp.exp(2*u(top))*(2*u(top))**11/mp.factorial(11) < mp.mpf("1e-65"))
    def gl(fun, lo, hi, count):
        xs, ws = nodes[count]
        return sum(weight*(hi-lo)/2*fun((hi+lo)/2+(hi-lo)*z/2) for z, weight in zip(xs, ws))
    def pressure_kernel(v, base, q):
        H = ex(2*(wbase(v)-wbase(base))+2*q*(u(v)-u(base)))
        T = mp.mpf(".4")*(mp.mpf(".02")+q*mp.mpf(".001"))/distance(v, q)
        return H*(wbprime(v)+q*up(v))*T
    def pressure(base, q, count):
        # Logarithmic coordinate resolves the near-boundary integration range.
        ell = mp.log(top/base)
        return -gl(lambda z: pressure_kernel(base*mp.exp(z), base, q)*base*mp.exp(z), 0, ell, count)
    def evaluate(base, q, count=64):
        if isinstance(q, Jet):
            local = [Jet(*vals) for vals in jetder(base, q.v[0])]
        else:
            local = jetfun(base, q)
        rho, null = tensor(*local)
        phi = local[5]
        Y = distance(base, q)
        c = mp.mpf(".02")+q*mp.mpf(".001")
        vac = mp.mpf(11)/120*mp.mpf(".41")**4*ex(4*phi)
        mass = mp.mpf(".4")*c*mp.mpf(".41")*ex(phi)
        T = mp.mpf(".4")*c/Y
        P = pressure(base, q, count)
        urho, unull, J = rho-vac-mass, null-mass, mass-T
        recombined = [P+J+urho, 2*P+J+unull]
        direct = [P-T+rho-vac, 2*P-T+null]
        return recombined, direct, T
    def finite_second(fun, q, h):
        if q == 0:
            return (2*fun(q)-5*fun(q+h)+4*fun(q+2*h)-fun(q+3*h))/h**2
        if q == 1:
            return (2*fun(q)-5*fun(q-h)+4*fun(q-2*h)-fun(q-3*h))/h**2
        return (fun(q+h)-2*fun(q)+fun(q-h))/h**2
    for ylabel in PROTOCOL["y"]:
        base = mp.mpf(ylabel)
        for qlabel in PROTOCOL["parameters"]:
            q = mp.mpf(qlabel)
            tag = f"fixture_{ylabel}_{qlabel}"
            Yad = mp.quad(lambda v: mp.exp(q*u(v)) if v else 1, [0, base/2, base])
            check(tag+"_distance_adaptive", abs(distance(base, q)-Yad) < mp.mpf("1e-30"))
            for count in PROTOCOL["gauss_nodes"]:
                integrand = lambda z: 16*base*z**15*mp.exp(q*u(base*z**16))
                yg = gl(integrand, 0, mp.mpf(".5"), count)+gl(integrand, mp.mpf(".5"), 1, count)
                check(tag+f"_distance_gauss_{count}", abs(yg-Yad) < mp.mpf("1e-30"), mp.nstr(abs(yg-Yad), 30))
                if ylabel == ".06" and qlabel == "1" and count == 32:
                    NUMERICS.append({"first_run_single_panel_error":
                                     mp.nstr(abs(gl(integrand, 0, 1, count)-Yad), 35)})
            ell = mp.log(top/base)
            Pad = -mp.quad(lambda z: pressure_kernel(base*mp.exp(z), base, q)*base*mp.exp(z),
                          [0, ell/2, ell])
            for count in PROTOCOL["gauss_nodes"]:
                check(tag+f"_pressure_gauss_{count}", abs(pressure(base, q, count)-Pad) < mp.mpf("1e-30"))
            val, direct, _ = evaluate(base, q)
            targets = evaluate(base, Jet(q, 1))[0]
            pressure_target = 2*pressure(base, Jet(q, 1), 64).v[2]
            pressure_ad = -mp.quad(lambda z: 2*pressure_kernel(base*mp.exp(z), base, Jet(q, 1)).v[2]*base*mp.exp(z),
                                  [0, ell/2, ell])
            check(tag+"_pressure_second_adaptive", abs(pressure_target-pressure_ad) < mp.mpf("1e-30"))
            for j, label in enumerate(["density", "null"]):
                check(tag+"_"+label+"_exact_recombination", abs(val[j]-direct[j]) < mp.mpf("1e-35"))
            radius = mp.mpf(".01")
            for count in PROTOCOL["contour_nodes"]:
                sums = [mp.mpc(0), mp.mpc(0)]
                for k in range(count):
                    phase = mp.exp(2j*mp.pi*k/count)
                    vals = evaluate(base, q+radius*phase)[0]
                    for j in range(2):
                        sums[j] += vals[j]/phase**2
                for j, label in enumerate(["density", "null"]):
                    target = 2*targets[j].v[2]
                    err = abs(2*sums[j]/count/radius**2-target)
                    check(tag+f"_{label}_contour_{count}", err < mp.mpf("1e-6")*(1+abs(target)), mp.nstr(err, 20))
            for j, label in enumerate(["density", "null"]):
                target = 2*targets[j].v[2]
                errors = []
                for step in PROTOCOL["real_steps"]:
                    approx = finite_second(lambda z: evaluate(base, z)[0][j], q, mp.mpf(step))
                    err = abs(approx-target)
                    check(tag+f"_{label}_difference_{step}", err < mp.mpf("1e-6")*(1+abs(target)), mp.nstr(err, 20))
                    errors.append(err)
                check(tag+f"_{label}_refinement", errors[1] < mp.mpf(".35")*errors[0]
                      or max(errors) < mp.mpf("1e-35"))
                NUMERICS.append({"y": ylabel, "lambda": qlabel, "component": label,
                                 "normalized_second": mp.nstr(target, 35)})
            tj = evaluate(base, Jet(q, 1))[2]
            check(tag+"_freezing_distance_omission_detected", abs(2*tj.v[2]) > mp.mpf("1e-15"))


def hashes():
    return {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    mp.mp.dps = 70
    before, error = hashes(), None
    print(json.dumps({"event": "start", "story": "S158", "input_count": len(before)}), flush=True)
    try:
        predecessors()
        data = tensor_polynomials()
        singular_and_seed_checks(data)
        budget_checks()
        fixtures(data)
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    ok = error is None and before == after and all(r["passed"] for r in ROWS)
    record = {
        "story": "S158", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if ok else "failed", "python": sys.version,
        "mpmath": mp.__version__, "sympy": s.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "numerics": NUMERICS,
        "exception": error, "scope": "actual-path massive canonical plus Weyl minus fixed vacuum outer response",
        "unresolved": ["other source sectors and radius equation", "general family and full inverse",
                       "coupled solution", "physical completion requirements"],
    }
    if error:
        record["failed_checker_source_utf8"] = Path(__file__).read_text(encoding="utf-8")
        record["failed_note_source_utf8"] = (ROOT/INPUTS[1]).read_text(encoding="utf-8")
    out = ROOT/"evidence"/(STEM+".json")
    history = json.loads(out.read_text()) if out.exists() else {"runs": []}
    history["runs"].append(record)
    out.write_text(json.dumps(history, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": record["status"], "checks": len(ROWS), "unchanged_inputs": before == after}), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
