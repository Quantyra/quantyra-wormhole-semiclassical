"""S142: exact identities and independent quadratures for the flat vertex.

The numerical tests are diagnostics. Uniform bounds and the infinite angular
sum are justified by the ordinary proof, not by the sampled quadratures.
"""
from __future__ import annotations

from datetime import datetime, timezone
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback

import numpy as np
import sympy as s
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[1]
STEM = "flat-variable-dirac-three-point"
INPUTS = [
    "checks/flat_variable_dirac_three_point.py",
    "notes/flat-variable-dirac-three-point.md",
    "notes/stationary-reduced-parameter-form-bounds.md",
    "notes/stationary-product-contour-inverse.md",
    "notes/stationary-centered-canonical-response.md",
    "notes/retained-benchmark-parameter-recovery.md",
]
ROWS = []
NUMERICS = []


def check(name, passed, detail=""):
    if any(row["name"] == name for row in ROWS):
        raise ValueError("duplicate label")
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not row["passed"]:
        raise AssertionError(name)


def simplex_integral(poly, x, y, z):
    result = 0
    for (a, b, c), coeff in s.Poly(s.expand(poly), x, y, z).terms():
        result += coeff*s.factorial(a)*s.factorial(b)*s.factorial(c) / s.factorial(a+b+c+2)
    return s.expand(result)


def exact_checks():
    x, y, z, Q1, Q2, Q3 = s.symbols("x y z Q1 Q2 Q3")
    D = x*y*Q1+y*z*Q2+z*x*Q3
    J = (s.Rational(1, 2)-z)*Q1+(s.Rational(1, 2)-x)*Q2+(s.Rational(1, 2)-y)*Q3
    integ = lambda p: simplex_integral(p, x, y, z)
    check("simplex_area", integ(1) == s.Rational(1, 2))
    order2 = (Q1+Q2+Q3)/6+2*integ(J-4*D)
    order4 = -(Q1**2+Q2**2+Q3**2)/60+2*integ(6*D**2-2*J*D)
    order6 = (Q1**3+Q2**3+Q3**3)/420+2*integ(-8*D**3+3*J*D**2)
    cross = Q1*Q2+Q1*Q3+Q2*Q3
    sixth = -(2*sum(a*a*b for a, b in
        [(Q1,Q2),(Q1,Q3),(Q2,Q1),(Q2,Q3),(Q3,Q1),(Q3,Q2)])
        +5*Q1*Q2*Q3)/420
    check("two_derivative_seagull_cancellation", s.expand(order2) == 0)
    check("fourth_order_kernel", s.expand(order4-cross/30) == 0)
    check("sixth_order_kernel", s.expand(order6-sixth) == 0)
    a, b, c, alpha = s.symbols("a b c alpha")
    action = 4*alpha*s.exp(-2*(a+b+c))*(Q1*a+Q2*b+Q3*c)**2
    third = s.diff(action, a, b, c).subs({a:0,b:0,c:0})
    check("direct_curvature_action_third_variation",
          s.simplify(third+16*alpha*cross) == 0)
    M = s.symbols("M", positive=True)
    check("single_mode_contact_matching",
          s.simplify(third.subs(alpha, -1/(960*s.pi*M*M))
                     -cross/(60*s.pi*M*M)) == 0)
    d, j = s.symbols("d j")
    remainder = (2+j)/(1+d)**2-2-(j-4*d)-(6*d*d-2*j*d)
    exact_remainder = d*d*(j*(3+2*d)-2*d*(4+3*d))/(1+d)**2
    check("exact_triangle_remainder",
          s.factor(remainder-exact_remainder) == 0)
    check("sixth_order_global_coefficient",
          s.Rational(43,108)+s.Rational(1,280) == s.Rational(3037,7560))

    # Clifford trace is checked before either numerical representation uses it.
    p, q, r, t, u, v = s.symbols("p q r t u v", real=True)
    g1, g2 = s.Matrix([[0,1],[1,0]]), s.Matrix([[0,-s.I],[s.I,0]])
    numerator = lambda a, b: M*s.eye(2)-s.I*(a*g1+b*g2)
    trace3 = s.trace(numerator(p,q)*numerator(r,t)*numerator(u,v))
    check("rank_two_triangle_trace",
          s.expand(trace3-2*M*(M*M-p*r-q*t-p*u-q*v-r*u-t*v)) == 0)
    r1, r2, r3 = s.Matrix([0,0]), s.Matrix([p,q]), s.Matrix([r,t])
    center = x*r1+y*r2+(1-x-y)*r3
    aa = [rr-center for rr in (r1,r2,r3)]
    pair = sum((aa[i].dot(aa[k]) for i,k in ((0,1),(0,2),(1,2))),s.Integer(0))
    substitution = {Q1:r2.dot(r2), Q2:(r2-r3).dot(r2-r3),
                    Q3:r3.dot(r3), z:1-x-y}
    check("shifted_triangle_numerator",
          s.expand(pair+3*D.subs(substitution)-J.subs(substitution)) == 0)
    nn, qq, CC = s.symbols("n q C", positive=True)
    angular = 2*(qq+2*nn)*CC**2/(nn**2*(nn+qq)**2)
    telescope = 2*CC**2/qq*(1/nn**2-1/(nn+qq)**2)
    check("angular_fourth_moment_identity", s.factor(angular-telescope) == 0)
    for charge, cutoff in ((1,3),(3,7),(37,11)):
        partial = sum(s.Rational(1,n*n)-s.Rational(1,(n+charge)**2)
                      for n in range(cutoff+1,cutoff+100+1))
        edge = sum(s.Rational(1,n*n) for n in range(cutoff+1,cutoff+charge+1))
        endpoint = sum(s.Rational(1,n*n) for n in range(cutoff+101,cutoff+100+charge+1))
        check(f"finite_tail_telescope_q{charge}", partial == edge-endpoint)
    # Rational upper bounds use pi>3, log(2)>2/3, M1>18.
    a6, a4 = s.Rational(3037,7560*3), s.Rational(103,60*3)
    check("summed_logarithmic_constants",
          8*a4+16*a6 < s.Rational(68,10)
          and 16*a4 < s.Rational(92,10))
    check("normalized_low_momentum_bound",
          2*a6*s.Rational(1408,100000)*s.Rational(85,10**34)
          < s.Rational(33,10**36))
    check("fourth_log_energy_coefficient",
          120*3*s.Rational(36,100) < 130
          and 130*s.Rational(85,10**34) < s.Rational(111,10**32))
    check("fourier_L1_embedding",
          s.Rational(1,8) < s.Rational(36,100)**2)


@lru_cache(maxsize=None)
def unit_gauss(n):
    x, w = leggauss(n)
    return (x+1)/2, w/2


def invariants(p1, p2):
    momenta = np.array([p1, p2, -np.asarray(p1)-np.asarray(p2)], dtype=float)
    return np.einsum("ij,ij->i", momenta, momenta), momenta


def parameter_vertex(p1, p2, n=96):
    """M=1; general M follows by scaling momenta and multiplying by M^2."""
    Q, _ = invariants(p1, p2)
    a, w = unit_gauss(n)
    bubble = sum(np.dot(w, np.log1p(a*(1-a)*qi)) for qi in Q)
    x, t = a[:,None], a[None,:]
    y, z = (1-x)*t, (1-x)*(1-t)
    weights = w[:,None]*w[None,:]*(1-x)
    d = x*y*Q[0]+y*z*Q[1]+z*x*Q[2]
    j = (.5-z)*Q[0]+(.5-x)*Q[1]+(.5-y)*Q[2]
    # Algebraic difference avoids subtracting two near-constant integrals.
    triangle = np.sum(weights*(j-4*d-2*d*d)/(1+d)**2)
    return float((bubble+2*triangle)/(2*np.pi))


def parameter_remainder(p1, p2, n=96):
    Q, _ = invariants(p1, p2)
    a, w = unit_gauss(n)
    # Positive integral for log(1+x)-x+x^2/2 avoids cancellation at small x.
    bubble = 0.
    for qi in Q:
        xx = a*(1-a)*qi
        f = xx**3*np.sum(w[None,:]*a[None,:]**2/
                        (1+xx[:,None]*a[None,:]),axis=1)
        bubble += np.dot(w,f)
    x, t = a[:,None], a[None,:]
    y, z = (1-x)*t, (1-x)*(1-t)
    weights = w[:,None]*w[None,:]*(1-x)
    d = x*y*Q[0]+y*z*Q[1]+z*x*Q[2]
    j = (.5-z)*Q[0]+(.5-x)*Q[1]+(.5-y)*Q[2]
    rem = d*d*(j*(3+2*d)-2*d*(4+3*d))/(1+d)**2
    return float((bubble+2*np.sum(weights*rem))/(2*np.pi))


def direct_loop_vertex(p1, p2, nr, nang):
    """Rational spin traces, integrated on R2 via r=tan(pi*t/2)."""
    _, momenta = invariants(p1,p2)
    t,w = unit_gauss(nr)
    radius = np.tan(np.pi*t/2)[:,None]
    theta = 2*np.pi*np.arange(nang)[None,:]/nang
    px, py = radius*np.cos(theta), radius*np.sin(theta)
    p_squared = px*px+py*py
    base_den = 1+p_squared
    b0 = 2*(1-p_squared)/base_den**2
    c0 = 2*(1-3*p_squared)/base_den**3
    value = np.zeros_like(px)
    for shift in momenta:
        ax, ay = px+shift[0], py+shift[1]
        value += 2*(1-px*ax-py*ay)/(base_den*(1+ax*ax+ay*ay))-b0
    ax, ay = px+momenta[0,0], py+momenta[0,1]
    bx, by = px-momenta[2,0], py-momenta[2,1]
    triangle = 2*(1-px*ax-py*ay-px*bx-py*by-ax*bx-ay*by)/(
        base_den*(1+ax*ax+ay*ay)*(1+bx*bx+by*by))
    value -= 2*(triangle-c0)
    jac = radius*(np.pi/2)*(1+radius**2)
    return float(np.dot(w,np.mean(value*jac,axis=1))/(2*np.pi))


def numerical_checks():
    cases = [
        ((.7,.2),(-.1,.8)), ((2.,0.),(0.,1.)),
        ((3.,0.),(-1.25,.7)), ((1.,0.),(-1.,0.)),
        ((.1,.3),(.2,-.1)), ((4.,1.),(-2.,3.)),
    ]
    for i,(p1,p2) in enumerate(cases):
        par = parameter_vertex(p1,p2)
        par2 = parameter_vertex(p1,p2,144)
        direct = [direct_loop_vertex(p1,p2,n,2*n) for n in (96,192,384)]
        tol = 3e-8+3e-7*abs(par)
        Q,momenta = invariants(p1,p2)
        rem = parameter_remainder(p1,p2)
        local = (Q[0]*Q[1]+Q[0]*Q[2]+Q[1]*Q[2])/(60*np.pi)
        NUMERICS.append({"case":i,"momenta":momenta.tolist(),"parameter":par,
                         "parameter_refined":par2,"direct_resolutions":[96,192,384],
                         "direct":direct,"tolerance":tol,"remainder":rem,
                         "fourth_contact":local})
        check(f"direct_loop_case_{i}",
              abs(par-direct[-1])<tol and abs(direct[-1]-direct[-2])<tol,
              NUMERICS[-1])
        check(f"parameter_refinement_and_contact_case_{i}",
              abs(par-par2)<1e-10+1e-9*abs(par)
              and abs(par-rem-local)<1e-10+1e-9*abs(par))
        swapped = parameter_vertex(momenta[1],momenta[2])
        check(f"permutation_symmetry_case_{i}",abs(swapped-par)<1e-10+1e-9*abs(par))
    aa,ww=unit_gauss(144)
    for momentum in (.2,1.,4.):
        xx=aa*(1-aa)*momentum**2
        independent=float(np.dot(ww,np.log1p(xx)-xx/(1+xx))/np.pi)
        vertex=parameter_vertex((momentum,0),(-momentum,0),144)
        check(f"constant_insertion_limit_{momentum}",
              abs(vertex-independent)<1e-10+1e-9*abs(vertex))
    a6,a4=3037/(7560*np.pi),103/(60*np.pi)
    for scale in (.01,.2,1.,4.,20.):
        p1,p2=np.array([.7,.2])*scale,np.array([-.1,.8])*scale
        Q,_=invariants(p1,p2)
        K=float(np.sqrt(max(Q)))
        rem=abs(parameter_remainder(p1,p2,144))
        check(f"uniform_bound_diagnostic_scale_{scale}",
              rem<=a6*K**6 and (K<1 or rem<=a4*K**4),
              {"K":K,"absolute_remainder":rem,"sixth_bound":a6*K**6,
               "high_momentum_fourth_bound":a4*K**4})


def hashes():
    return {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    initial=hashes()
    error=None
    try:
        exact_checks()
        numerical_checks()
    except Exception:
        error=traceback.format_exc()
        print(error,file=sys.stderr)
    final=hashes()
    ok=error is None and initial==final and all(r["passed"] for r in ROWS)
    record={
        "story":"S142","created_utc":datetime.now(timezone.utc).isoformat(),
        "status":"passed" if ok else "failed","python":sys.version,
        "numpy":np.__version__,"sympy":s.__version__,
        "input_sha256_before":initial,"input_sha256_after":final,
        "unchanged_inputs":initial==final,"checks":ROWS,"numerics":NUMERICS,
        "exception":error,
        "scope":"flat basepoint third variation, common contact identity and "
                "all-mode bounds; quadratures are non-rigorous diagnostics",
        "unresolved":["curved and boundary completion","finite-background response",
                      "exterior solution and signed source","physical stability"],
    }
    if error:
        record["failed_checker_source_utf8"]=Path(__file__).read_text(encoding="utf-8")
    out=ROOT/"evidence"/f"{STEM}.json"
    history=json.loads(out.read_text(encoding="utf-8")) if out.exists() else {"runs":[]}
    history["runs"].append(record)
    out.write_text(json.dumps(history,indent=2)+"\n",encoding="utf-8",newline="\n")
    with (ROOT/"evidence"/f"{STEM}.log").open("a",encoding="utf-8",newline="\n") as log:
        log.write(json.dumps({"created_utc":record["created_utc"],"status":record["status"],
                              "inputs":initial})+"\n")
        for row in ROWS:
            log.write(json.dumps(row)+"\n")
        if error:
            log.write(error+"\n")
    print(json.dumps({"status":record["status"],"checks":len(ROWS),
                      "failed":sum(not r["passed"] for r in ROWS)}))
    return 0 if ok else 1


if __name__=="__main__":
    sys.exit(main())
