"""S160: complete retained Dirac/reference source on the comparison path."""
from __future__ import annotations
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import importlib.util
import sys
import traceback
import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
# The pinned embedded Python does not put the script directory on sys.path.
_spec = importlib.util.spec_from_file_location("s158_taylor_arithmetic", ROOT/"checks/comparison_path_outer_recombination.py")
_module = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_module)
Jet = _module.Jet
STEM = "comparison-path-complete-source-response"
INPUTS = [
    "checks/comparison_path_complete_source_response.py",
    "notes/comparison-path-complete-source-response.md",
    "checks/comparison_path_outer_recombination.py",
    *[f"{folder}/{stem}.{ext}" for stem in [
        "comparison-path-interior-response", "comparison-path-cold-response",
        "comparison-path-outer-regular-response", "comparison-path-outer-recombination",
        "comparison-path-lll-reference-response", "fixed-chart-global-constraint-transfer"]
      for folder, ext in [("notes", "md"), ("evidence", "json")]],
    "notes/stationary-completed-metric-response-norm.md",
    "notes/stationary-centered-local-response.md",
    "notes/spherical-anomaly-sources.md",
    "notes/stationary-centered-covariant-seed.md",
    "notes/stationary-centered-exterior-source.md",
    "notes/stationary-reference-cold-bolt.md",
    "notes/stationary-horizon-source-defect.md",
]
PROTOCOL = {
    "digits": 70, "y": [".002", ".02", ".06"], "parameters": ["0", ".4", "1"],
    "formula_tolerance": "1e-30", "contour_nodes": [64, 128], "contour_radius": ".01",
    "real_steps": [".001", ".0005"], "second_tolerance": "1e-6",
    "refinement_ratio": ".35", "refinement_floor": "1e-35",
    "endpoint_differences": "second-order one-sided within [0,1]",
    "target": "symbolic jet parameter derivatives and independent Taylor polynomial arithmetic",
    "fixtures_are_physical_states": False,
}
ROWS, NUMERICS = [], []
Rat = s.Rational


def check(name, passed, detail=""):
    assert name not in {r["name"] for r in ROWS}
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not passed:
        raise AssertionError(name)


def predecessors():
    for stem, count in [
        ("comparison-path-interior-response", 126), ("comparison-path-cold-response", 242),
        ("comparison-path-outer-regular-response", 202), ("comparison-path-outer-recombination", 242),
        ("comparison-path-lll-reference-response", 376), ("fixed-chart-global-constraint-transfer", 99)]:
        run = json.loads((ROOT/"evidence"/(stem+".json")).read_text())["runs"][-1]
        check(stem+"_receipt", run["status"] == "passed" and len(run["checks"]) == count
              and len({r["name"] for r in run["checks"]}) == count and all(r["passed"] for r in run["checks"]))
        check(stem+"_current_inputs", run["unchanged_inputs"]
              and run["input_sha256_before"] == run["input_sha256_after"]
              and all(sha256((ROOT/p).read_bytes()).hexdigest() == h
                      for p, h in run["input_sha256_after"].items()))


def lagrangian(R, U, V, phi, e, qm):
    aa, cc, bb = Rat(11, 360), Rat(1, 20), Rat(1, 30)
    dd = bb-2*aa/3
    L = (phi*(cc*(R+2*e)**2/3-4*aa*e*R-qm/3)+bb*R*U-3*bb*U**2
         +(bb*R+(2*bb-4*aa)*e)*V-6*dd*U*V-3*dd*V**2)
    P = phi*(2*cc*(R+2*e)/3-4*aa*e)+bb*(U+V)
    Q = bb*R-6*bb*U-6*dd*V
    Z = bb*R+(2*bb-4*aa)*e-6*dd*(U+V)
    Lphi = cc*(R+2*e)**2/3-4*aa*e*R-qm/3
    return [s.expand(f) for f in (L, P, Q, Z, Lphi)]


def cartesian_polynomials():
    H, e, qm = s.symbols("H e qm")
    indices = [(j-k, k) for j in range(6) for k in range(j+1)]
    sig = {a: s.Symbol(f"s{a[0]}{a[1]}") for a in indices}
    phi = {a: s.Symbol(f"p{a[0]}{a[1]}") for a in indices}
    directions = []
    for axis in range(2):
        mapping = {H: -2*H*sig[(1, 0) if axis == 0 else (0, 1)]}
        for jets in (sig, phi):
            for (i, j), symbol in jets.items():
                if i+j < 5:
                    mapping[symbol] = jets[(i+1, j) if axis == 0 else (i, j+1)]
        directions.append(mapping)
    def D(f, axis):
        return s.expand(sum(s.diff(f, v)*directions[axis][v]
                            for v in f.free_symbols if v in directions[axis]))
    R = -2*H*(sig[2, 0]+sig[0, 2])
    U = H*(phi[2, 0]+phi[0, 2])
    V = H*(phi[1, 0]**2+phi[0, 1]**2)
    L, P, Q, Z, Lphi = lagrangian(R, U, V, phi[0, 0], e, qm)
    dp, dq = [[D(f, i) for i in range(2)] for f in [P, Q]]
    boxp = H*(D(dp[0], 0)+D(dp[1], 1))
    hxx = D(dp[0], 0)-sig[1, 0]*dp[0]+sig[0, 1]*dp[1]
    hyy = D(dp[1], 1)+sig[1, 0]*dp[0]-sig[0, 1]*dp[1]
    hxy = D(dp[0], 1)-sig[0, 1]*dp[0]-sig[1, 0]*dp[1]
    common = -L/2+R*P/2+boxp+(H*(dq[0]*phi[1, 0]+dq[1]*phi[0, 1])+Q*U)/2
    mxx = s.expand(common-H*hxx+H*Z*phi[1, 0]**2-H*phi[1, 0]*dq[0])
    myy = s.expand(common-H*hyy+H*Z*phi[0, 1]**2-H*phi[0, 1]*dq[1])
    mxy = s.expand(-H*hxy+H*Z*phi[1, 0]*phi[0, 1]-H*(phi[1, 0]*dq[1]+phi[0, 1]*dq[0])/2)
    E = s.expand(Lphi+H*(D(dq[0], 0)+D(dq[1], 1))
                 -2*H*(D(Z, 0)*phi[1, 0]+D(Z, 1)*phi[0, 1]+Z*(phi[2, 0]+phi[0, 2])))
    divx = D(mxx, 0)+D(mxy, 1)+sig[1, 0]*(mxx-myy)+2*sig[0, 1]*mxy
    divy = D(mxy, 0)+D(myy, 1)+sig[0, 1]*(myy-mxx)+2*sig[1, 0]*mxy
    check("two_dimensional_Ward_x", s.expand(divx+E*phi[1, 0]/2) == 0)
    check("two_dimensional_Ward_y", s.expand(divy+E*phi[0, 1]/2) == 0)
    check("independent_conformal_metric_variation",
          s.expand(2*L-2*R*P-2*Q*U-2*Z*V-2*boxp+2*(mxx+myy)) == 0)
    variables = [H, e, qm, *sig.values(), *phi.values()]
    orders = [0, 0, 0, *[sum(a) for a in indices], *[sum(a) for a in indices]]
    phizero = variables.index(phi[0, 0])
    for label, expr in [("density", myy), ("null", myy-mxx), ("radius", E)]:
        for degree in [4, 5]:
            cur = s.expand(expr if degree == 4 else D(expr, 0))
            poly = s.Poly(cur, *variables)
            terms = poly.terms()
            check(f"{label}_total_order_{degree}", all(sum(a*b for a, b in zip(ex, orders)) <= degree for ex, co in terms))
            check(f"{label}_H_power_{degree}", all(ex[0] <= 2 for ex, co in terms))
            check(f"{label}_phi_degree_{degree}", all(ex[phizero] <= 1 for ex, co in terms))
            count = sum(abs(co)*max(1, sum(ex[3:])+2*ex[0]) for ex, co in terms if ex[2] == 0)
            check(f"{label}_coefficient_reserve_{degree}", count < (10**4 if degree == 4 else 10**5), str(count))
            NUMERICS.append({"polynomial": label, "order": degree, "monomials": len(terms), "weighted_count": str(count)})
    check("density_magnetic_monomial", s.diff(myy, qm) == phi[0, 0]/6)
    check("null_magnetic_cancellation", s.diff(myy-mxx, qm) == 0)
    check("radius_magnetic_constant", s.diff(E, qm) == -Rat(1, 3))
    return H, e, qm, sig, phi, myy, s.expand(myy-mxx), E, mxy


def radial_polynomials():
    w, p = s.symbols("w0:5"), s.symbols("p0:5")
    e, qm = s.symbols("e qm")
    def D(f):
        return s.expand(sum(s.diff(f, a[j])*a[j+1] for a in (w, p) for j in range(4)))
    R = -2*(w[2]+w[1]**2)
    U, V = p[2]+w[1]*p[1], p[1]**2
    L, P, Q, Z, Lphi = lagrangian(R, U, V, p[0], e, qm)
    mt = s.expand(-L/2+R*P/2+D(D(P))+(D(Q)*p[1]+Q*U)/2)
    mr = s.expand(-L/2+R*P/2+w[1]*D(P)+Z*V-D(Q)*p[1]/2+Q*U/2)
    E = s.expand(Lphi+D(D(Q))+w[1]*D(Q)-2*(D(Z*p[1])+w[1]*Z*p[1]))
    adj1 = lambda f: D(f)+w[1]*f
    adj2 = lambda f: D(D(f))+2*w[1]*D(f)+(w[2]+w[1]**2)*f
    check("independent_radial_orbit_Euler", s.expand(L-adj1(s.diff(L, w[1]))+adj2(s.diff(L, w[2]))+2*mt) == 0)
    en = L-2*R*P-2*Q*U-2*Z*V-adj1(2*w[1]*P-p[1]*Q)
    check("independent_radial_lapse_Euler", s.expand(en+2*mr) == 0)
    ep = s.diff(L, p[0])-adj1(s.diff(L, p[1]))+adj2(s.diff(L, p[2]))
    check("independent_radius_Euler", s.expand(ep-E) == 0)
    y, c, f, d = s.symbols("y c f d", positive=True)
    wb, pb = e*y**2/2+c*y**3, -s.log(d*y)-e*y*y/6+f*y**3
    sub = {a[j]: s.diff(expr, y, j) for a, expr in [(w, wb), (p, pb)] for j in range(5)}
    outer = s.expand(E.subs(sub))
    expected = {-4: -Rat(11, 15), -3: 0, -2: 22*e/45, -1: 11*c/15}
    for k, val in expected.items():
        check(f"radius_Laurent_{k}", s.simplify(outer.coeff(y, k)-val) == 0)
    vacuum = -Rat(11, 15)/y**4*s.exp(-2*e*y*y/3+4*f*y**3)
    vp = s.series(vacuum, y, 0, 1).removeO().expand()
    for k in [-4, -3, -2, -1]:
        check(f"correlated_vacuum_cancellation_{k}", s.simplify((outer-vp).coeff(y, k).subs(f, -c/4)) == 0)
    mismatch = s.symbols("mismatch")
    check("uncorrelated_vacuum_pole_detected",
          s.simplify((outer-vp).coeff(y, -1).subs(f, -c/4+mismatch)-44*mismatch/15) == 0)
    regular = s.expand(outer.subs(f, -c/4)-sum(v*y**k for k, v in expected.items()))
    poly = s.Poly(regular, y, c, e, qm)
    count = sum(abs(co)*2**ex[1] for ex, co in poly.terms())
    check("radius_reference_regular_majorant", count < 10**8, str(count))
    Eproduct = lagrangian(-2, 0, 0, 0, e, qm)[4]
    check("nonzero_angular_product_value", s.expand(Eproduct-((1+e**2)/15+e/9-qm/3)) == 0)
    return w, p, e, qm, mt, s.expand(mt-mr), E


def angular_identities():
    E, mt, mr, kE, K, pa0, pv, dp, t, Pc, Sc = s.symbols("E mt mr kE K pa0 pv dp t Pc Sc")
    pa = pa0/t**4+pv*(1-t**-4)+dp*(-t**-4+4*t**-5-3*t**-6)
    ps = -pa0/t**4+pv*(1+t**-4)+dp*(t**-4-6*t**-5+6*t**-6)
    rho = Pc-Sc+4*kE*mt+K*t**4*pa
    null = 2*Pc-Sc+4*kE*(mt-mr)
    angular = -Sc/2+kE*E+2*kE*(mt+mr)-K*t**4*ps
    T = kE*E-K*(2*pv*t**4+dp*(-2/t+3/t**2))
    check("independent_angular_chain_rule", s.simplify(angular-rho+null/2-T) == 0)
    estar = s.symbols("estar")
    cold = kE*(E-estar)-K*(2*pv*(t**4-1)+dp*(-2/t+3/t**2-1))
    check("cold_trace_centering", s.expand((T-cold).subs(estar, K*(2*pv+dp)/kE)) == 0)
    check("missing_angular_product_detected", (kE*estar).subs({kE: 1, estar: 2}) != 0)
    ev = s.symbols("ev")
    check("outer_trace_vacuum_allocation", s.expand((kE*(E-ev)-K*dp*(-2/t+3/t**2)-T).subs(ev, 2*K*pv*t**4/kE)) == 0)


def budgets():
    X = s.symbols("X")
    H = Rat("1.1")*X+sum(2*s.catalan(j-1)*Rat(".3")**(j-1)/Rat(".9")**(2*j-1)*X**j for j in range(2, 6))
    comp = s.series(sum((100*H)**k for k in range(6)), X, 0, 6).removeO().expand()
    for j in range(1, 6):
        check(f"cold_chart_composition_{j}", comp.coeff(X, j)/100**j < 4)
        check(f"cold_seed_jet_{j}", 20*Rat("1.28")**j+33*Rat("4e-6") < 100)
    for radius, bound in [(Rat(1, 64), Rat(".04")), (Rat(".1"), Rat(".22"))]:
        check("complex_cold_s_displacement_"+str(radius), 2*(radius+radius**2/2)/(1-radius**2/4) < bound)
    check("inverse_tau_jet_sum", Rat("1.3")*sum(Rat(".2")**j for j in range(6)) < 2)
    A = 33*Rat("8e-5")*81
    check("complex_gap_direction_jets", A < Rat(".214"))
    finite_radius = 2*33*Rat("8e-5")*(1+A)**4*Rat("1.001")**4*Rat("1.00001")
    check("finite_cold_radius_ideal", finite_radius < Rat(".02"), str(finite_radius))
    check("cold_background_jet_guard", 100*Rat(".00331")+Rat(".01") < 1)
    check("cold_exponential_guard", mp.exp(mp.mpf(".2")) < 2)
    # S155 inverse-quadratic majorant for the radius scalar, no log Jacobian.
    amin = Rat(".999999")
    HH = Rat("1.000001")*X+sum(2*s.catalan(j-1)*Rat("4e-6")**(j-1)/amin**(2*j-1)*X**j for j in range(2, 5))
    DD = Rat("5e-9")+HH-X
    transfer = (Rat("1e-4")*sum((100*HH)**j for j in range(5))
                +6*((1-160*(X+DD))**-2-(1-160*X)**-2))
    series = s.series(transfer, X, 0, 5).removeO().expand()
    nu_limits = {2: Rat(".000342"), 3: Rat(".000796"), 4: Rat(".001854")}
    for j in range(2, 5):
        val = series.coeff(X, j)/100**j
        check(f"middle_radius_fixed_jet_{j}", val < nu_limits[j], str(val))
    n = Rat("1.2e-8")+Rat("6.3e-8")*X+sum(v*X**j for j,v in nu_limits.items())
    J = 1+2*sum((Rat("1.28")*X)**j for j in range(1,5))
    nj = s.expand(n*J)
    q = 33*1001*sum(nj.coeff(X,j)*X**j for j in range(1,5))
    # exp(4e-7)<1.000001, and all truncated coefficients are positive.
    eq = s.series(sum(q**k/s.factorial(k) for k in range(1,5)),X,0,5).removeO()
    hseries = s.series(Rat("4e-7")*J+J*Rat("1.000001")*eq/1001,X,0,5).removeO()
    zseries = (hseries-hseries.subs(X,0))/(1-Rat("4e-7"))
    logmajorant = s.series(sum(zseries**k/k for k in range(1,5)),X,0,5).removeO().expand()
    for j,lim in [(2,".012"),(3,".085"),(4,".42")]:
        v = logmajorant.coeff(X,j)
        check(f"middle_nonlinear_radius_{j}", v < Rat(lim), str(v))
        check(f"middle_radius_background_{j}", 6*Rat(".32")**j+v < 1)
    check("middle_tau_complex_ratio", mp.exp(mp.mpf(".022")) < 2)
    check("middle_tau_complex_displacement", Rat("1.2")*Rat(".018") < Rat(".022"))
    check("middle_tau_bracket_strip", 1+Rat("1.4")*Rat(".02")+Rat("1.1")*Rat(".02")**2 < Rat("1.1"))
    check("middle_low_radius_gradient", Rat("20.1")+33*Rat("6.3e-6")*Rat("1.000001")/(1-Rat("4e-7"))+Rat("1e-12") < 100)
    check("middle_complex_seed_radius_log", mp.log(40)+mp.log(4) < 6)
    check("middle_relative_complex_z", Rat(".5")*Rat(".035")*40 < Rat(".75"))
    check("middle_complex_strip_derivative", Rat(".3")+2*Rat(".02")+17*Rat(".02")**2/2+13*Rat(".02")**3/6 < Rat(".4"))
    check("middle_complex_curvature_guard", 2+17*Rat(".02")+13*Rat(".02")**2/2 < 3)
    check("middle_normal_coordinate_margin", Rat(1, 32)/(1-3*Rat(1, 32)**2) < Rat(".035"))
    check("middle_complex_z_speed", Rat(".2")+Rat(".4")*Rat(".02") < Rat(".5")**2)
    check("middle_value_radius_log", mp.log(40)+mp.mpf("4e-7") < 4)
    eta, qbar = Rat("8.5e-33"), 160000
    b4 = 2*10**4*100*24*100**4
    b5 = 2*10**5*100*120*100**5
    cold = eta/3*(b4+qbar*100/6)
    check("cold_Weyl_and_radius_complex", cold < Rat("1e-15"), str(cold))
    bolt = Rat("2.2")*eta/3*(b5+qbar*100*100/6)
    check("bolt_Weyl_modulus", bolt < Rat("1e-11"), str(bolt))
    mid = eta/3*(2*10**4*4*24*100**4+qbar*4/6)
    check("middle_Weyl_and_radius_complex", mid < Rat("1e-16"), str(mid))
    ref = Rat("1.71e-26")
    check("cold_T_budget", Rat("1e-15")+ref*(2*Rat("4.1")+3*Rat("8.1")) < Rat("1e-14"))
    check("middle_T_budget", Rat("1e-16")+ref*(2*41**4+3*(2*Rat("1.001")+3*Rat("1.001")**2)) < Rat("1e-15"))
    check("outer_reference_E_budget", eta/12*10**8*(1+qbar) < Rat("1e-15"))
    check("outer_T_budget", eta/12*b4+Rat("1e-15")+ref*3*(2*Rat(".026")+3*Rat(".026")**2) < Rat("1e-12"))
    check("outer_vacuum_exponential_remainder",
          mp.mpf(11)/30*(mp.mpf(2)/3+2*mp.mpf(".081"))**2*mp.exp(mp.mpf(".006")) < 1)
    dg = Rat("3e-7")
    B = list(map(Rat, ["1e-15", "1e-11", "1e-16"]))
    N1, N2 = [[v*factor for v in B] for factor in [Rat(".032"), Rat(".002")]]
    for i in range(3):
        check(f"Weyl_Cauchy_first_{i}", B[i]/32 < N1[i])
        check(f"Weyl_Cauchy_second_{i}", 2*B[i]/32**2 < N2[i])
    mass = [Rat("1.003")*c+Rat("4.02")*dg*b+Rat("5.54")*dg**2*a for a, b, c in zip(B, N1, N2)]
    lapse = [Rat("1.006")*c+Rat("10.08")*dg*b+Rat("26.5")*dg**2*a for a, b, c in zip(B, N1, N2)]
    outs = [(46*mass[0]+279*mass[2])/2, (Rat(".000002")*mass[0]+mass[2])/2,
            (13122*lapse[0]+486*lapse[1]+1300*lapse[2])/2]
    for name, val, limit in zip(["mass", "density", "lapse"], outs, ["8e-17", "1.1e-19", "5e-12"]):
        check("new_Weyl_second_"+name, val < Rat(limit), str(val))
    full = [Rat("7.2e-8")+Rat("4.4e-12")+Rat("1.5e-14")+Rat("2.9e-19")+outs[0],
            Rat("3.2e-15")+Rat("1.6e-14")+Rat("7.4e-11")+Rat("1.1e-21")+outs[1],
            Rat("3.6e-5")+Rat("2.1e-11")+Rat("1e-15")+Rat("1.4e-18")+outs[2]]
    for name, val, limit in zip(["mass", "density", "lapse"], full, ["7.21e-8", "7.41e-11", "3.61e-5"]):
        check("complete_source_second_"+name, val < Rat(limit), str(val))
    pairs = [[Rat(x) for x in row] for row in [
        ["1.11e-6", "1.71e-10", "3.11e-9"], ["1.71e-9", "1.41e-12", "3.11e-14"],
        ["3.12e-7", "3.21e-11", "2.1e-12"]]]
    trace = [[Rat(x) for x in row] for row in [
        ["1e-14", "3.2e-16", "2e-17"], ["1e-15", "3.2e-17", "2e-18"], ["1e-12", "3.2e-14", "2e-15"]]]
    angular = [["1.7e-6", "2.6e-10", "4.7e-9"], ["2.6e-9", "2.2e-12", "4.7e-14"],
               ["4.7e-7", "4.9e-11", "3.2e-12"]]
    component_rows = [
        [["1.1e-6", "1.7e-10", "3.1e-9"], ["1e-15", "3.2e-17", "2e-18"], ["1e-23", "3.2e-25", "2e-26"]],
        [["1.7e-9", "1.4e-12", "3.1e-14"], ["1e-16", "3.2e-18", "2e-19"], ["1e-18", "3.2e-20", "2e-21"]],
        [["3.1e-7", "7e-14", "9.2e-14"], ["1e-9", "3.2e-11", "2e-12"], ["1e-21", "3.2e-23", "2e-24"]],
    ]
    for i in range(3):
        for j in range(3):
            check(f"complete_orbit_region_{i}_order_{j}", sum(Rat(row[j]) for row in component_rows[i]) < pairs[i][j])
            check(f"angular_source_region_{i}_order_{j}", Rat(3, 2)*pairs[i][j]+trace[i][j] < Rat(angular[i][j]))
    hmax = Rat("1e-30")
    J1, J2 = 4*Rat("1e-8")*hmax, 20*Rat("1e-16")*hmax
    for label, vals, limits in [
        ("angular", list(map(Rat, angular[0])), ["1.1e-8", "2e-7"]),
        ("orbit", pairs[0], ["7.1e-9", "1.3e-7"])]:
        first = Rat(".506")*81*(Rat("1.000001")*vals[1]+J1*vals[0])
        second = Rat(".506")*81*(Rat("1.000001")*vals[2]+2*J1*vals[1]+J2*vals[0])
        check(label+"_moving_bolt_first", first < Rat(limits[0]))
        check(label+"_moving_bolt_second", second < Rat(limits[1]))
    NUMERICS.append({"new_Weyl_second_exact": [str(v) for v in outs],
                     "complete_source_second_exact": [str(v) for v in full],
                     "complete_E0_residual": False, "goal_complete": False})


def fixtures(cart, radial):
    H, ec, qc, sig, phi, mtc, mnc, Ec, mxy = cart
    w, p, er, qr, mtr, mnr, Er = radial
    y, lam = s.symbols("y lam", positive=True)
    w0 = Rat(".001")*y*y/2+Rat(".02")*y**3
    p0 = -s.log(Rat(".41")*y)-Rat(".001")*y*y/6-Rat(".02")*y**3/4
    u = Rat(".001")*y**3+Rat(".002")*y**4*(1-s.log(y))
    wf, pf = w0+lam*u, p0-Rat(".0005")*lam*y**3+Rat(".003")*lam*y**4*(1-s.log(y))
    Dx, Dr = lambda f: s.exp(w0)*s.diff(f, y), lambda f: s.exp(-lam*u)*s.diff(f, y)
    def jets(f, D):
        ans = [f]
        for _ in range(4):
            ans.append(s.expand(D(ans[-1])))
        return ans
    cjets, rjets = [jets(f, Dx) for f in (wf, pf)], [jets(f, Dr) for f in (wf, pf)]
    cv = [H, *sig.values(), *phi.values()]
    cexpr = [s.exp(-2*wf)]
    for family in cjets:
        cexpr += [family[i] if j == 0 and i <= 4 else 0 for i, j in sig]
    cfun = s.lambdify((y, lam), cexpr, "mpmath", cse=True)
    cder = s.lambdify((y, lam), [[s.diff(f, lam, j)/s.factorial(j) for j in range(3)]
                                for f in cexpr], "mpmath", cse=True)
    cpoly = s.lambdify(cv, [f.subs({ec: Rat(".001"), qc: 4}) for f in (mtc, mnc, Ec)], "mpmath", cse=True)
    rfun = s.lambdify((y, lam), rjets[0]+rjets[1], "mpmath", cse=True)
    rpoly = s.lambdify((*w, *p), [f.subs({er: Rat(".001"), qr: 4}) for f in (mtr, mnr, Er)], "mpmath", cse=True)
    phifun = s.lambdify((y, lam), pf, "mpmath")
    phider = s.lambdify((y, lam), s.diff(pf, lam), "mpmath")
    def value(yy, q):
        ans = cpoly(*cfun(yy, q))
        ans[2] += mp.mpf(11)/15*mp.mpf(".41")**4*mp.exp(4*phifun(yy, q))
        return ans
    def target(yy, q):
        ans = cpoly(*(Jet(*v) for v in cder(yy, q)))
        ph = Jet(phifun(yy, q), phider(yy, q))
        ans[2] += mp.mpf(11)/15*mp.mpf(".41")**4*(4*ph).exp()
        return [2*v.v[2] for v in ans]
    def fd(fun, q, step):
        if q == 0:
            return (2*fun(q)-5*fun(q+step)+4*fun(q+2*step)-fun(q+3*step))/step**2
        if q == 1:
            return (2*fun(q)-5*fun(q-step)+4*fun(q-2*step)-fun(q-3*step))/step**2
        return (fun(q+step)-2*fun(q)+fun(q-step))/step**2
    for ys in PROTOCOL["y"]:
        yy = mp.mpf(ys)
        for qs in PROTOCOL["parameters"]:
            q = mp.mpf(qs)
            tag = f"fixture_{ys}_{qs}"
            va, vb = cpoly(*cfun(yy, q)), rpoly(*rfun(yy, q))
            for j, label in enumerate(["density", "null", "radius"]):
                check(tag+"_"+label+"_radial_vs_Cartesian", abs(va[j]-vb[j]) < mp.mpf("1e-30"))
            targ = target(yy, q)
            radius = mp.mpf(".01")
            for n in PROTOCOL["contour_nodes"]:
                sums = [mp.mpc(0)]*3
                for k in range(n):
                    phase = mp.exp(2j*mp.pi*k/n)
                    vals = value(yy, q+radius*phase)
                    sums = [a+b/phase**2 for a, b in zip(sums, vals)]
                for j, label in enumerate(["density", "null", "radius_regular"]):
                    err = abs(2*sums[j]/n/radius**2-targ[j])
                    check(tag+f"_{label}_contour_{n}", err < mp.mpf("1e-6")*(1+abs(targ[j])), mp.nstr(err, 25))
            for j, label in enumerate(["density", "null", "radius_regular"]):
                errs = []
                for step in PROTOCOL["real_steps"]:
                    err = abs(fd(lambda z: value(yy, z)[j], q, mp.mpf(step))-targ[j])
                    errs.append(err)
                    check(tag+f"_{label}_difference_{step}", err < mp.mpf("1e-6")*(1+abs(targ[j])), mp.nstr(err, 25))
                check(tag+f"_{label}_refinement", errs[1] < mp.mpf(".35")*errs[0] or max(errs) < mp.mpf("1e-35"))
            NUMERICS.append({"y": ys, "lambda": qs, "second_responses": [mp.nstr(v, 35) for v in targ]})
    x, z = s.symbols("x z")
    r2 = x*x+z*z
    sf, pf = Rat(".02")*r2+lam*Rat(".003")*r2, Rat(".01")+Rat(".03")*r2+lam*Rat(".002")*r2
    sub = {H: 1, ec: Rat(".001"), qc: 4}
    for family, expr in [(sig, sf), (phi, pf)]:
        sub.update({var: s.diff(expr, x, i, z, j).subs({x: 0, z: 0}) for (i, j), var in family.items()})
    check("Cartesian_bolt_null_zero", s.expand(mnc.subs(sub)) == 0)
    check("Cartesian_bolt_offdiagonal_zero", s.expand(mxy.subs(sub)) == 0)
    product_sub = {H: 1, ec: Rat(".001"), qc: 4}
    sh = -s.log(1-r2/4)
    for family, expr in [(sig, sh), (phi, s.Integer(0))]:
        product_sub.update({var: s.diff(expr, x, i, z, j).subs({x: 0, z: 0}) for (i, j), var in family.items()})
    check("Cartesian_product_orbit_zero", mtc.subs(product_sub) == 0 and mnc.subs(product_sub) == 0)
    check("Cartesian_product_angular_nonzero", Ec.subs(product_sub) != 0)


def hashes():
    return {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    mp.mp.dps = 70
    before, error = hashes(), None
    print(json.dumps({"event": "start", "story": "S160", "input_count": len(before)}), flush=True)
    try:
        predecessors()
        cart = cartesian_polynomials()
        radial = radial_polynomials()
        angular_identities()
        budgets()
        fixtures(cart, radial)
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    ok = error is None and before == after and all(r["passed"] for r in ROWS)
    record = {
        "story": "S160", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if ok else "failed", "python": sys.version,
        "mpmath": mp.__version__, "sympy": s.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after, "unchanged_inputs": before == after,
        "checks": ROWS, "numerics": NUMERICS, "exception": error,
        "scope": "complete retained Dirac-minus-reference source response on the specified comparison path",
        "unresolved": ["E0 and complete residual assembly", "general family and full inverse",
                       "coupled solution and physical completion"],
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
