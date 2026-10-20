"""S143: exact constants and independent finite-matrix determinant diagnostics.

The continuum resolvent/domain/trace-class argument is an ordinary proof.
Finite matrices test the variation algebra, not the exterior PDE theorem.
"""
from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import math
import sys
import traceback

import mpmath as mp
import numpy as np
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = "curved-remote-dirac-response"
INPUTS = [
    "checks/curved_remote_dirac_response.py",
    "notes/curved-remote-dirac-response.md",
    "notes/flat-variable-dirac-three-point.md",
    "notes/retained-benchmark-parameter-recovery.md",
    "notes/stationary-curved-quantum-tt.md",
    "notes/stationary-centered-canonical-response.md",
    "notes/stationary-paired-dirac-phase.md",
    "notes/stationary-massive-boundary-layer.md",
    "notes/stationary-reduced-cold-response.md",
]
ROWS = []
NUMERICS = []


def check(name, passed, detail=""):
    if any(row["name"] == name for row in ROWS):
        raise ValueError(f"duplicate label: {name}")
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not row["passed"]:
        raise AssertionError(name)


def exact_checks():
    # For 0<x<1, exp(x)<1/(1-x); exp(-x)>1-x.
    check("weighted_inverse_reserve", F(1, 1)/(F(99,100)-F(1,2)) < 3)
    check("conformal_mass_value", 1/(1-F(6,100)) < F(107,100))
    check("flat_mass_gradient", F(107,100)*(10+F(1,100)/(1-F(5,100))) < 11)
    check("squared_local_potential", F(107,100)**2+F(11,18) < 3)
    # sqrt(2)<10/7 suffices for the two-variable Cauchy estimate.
    cauchy = F(20,7)*F(26,1000)/((F(1,32)-F(1,180))*(1-F(26,1000)))
    check("analytic_chart_gradient", cauchy < 3)
    check("analytic_chart_value", F(26,1000)/(1-F(26,1000)) < F(5,100))
    check("evaluation_disk_fits_chart", F(1,180) < F(1,32))
    check("caccioppoli_reserve", 4*40**2+6 < 81**2)

    t = s.symbols("t", real=True)
    H = 1-10*t**3+15*t**4-6*t**5
    check("quintic_endpoints", [H.subs(t,0),H.subs(t,1)] == [1,0])
    check("quintic_first_derivative", s.expand(s.diff(H,t)+30*t**2*(1-t)**2) == 0)
    # t(1-t)<=1/4 gives |H'|<=15/8. H''' zeros give |H''|<=10/sqrt(3).
    critical = [(3-s.sqrt(3))/6, (3+s.sqrt(3))/6]
    check("quintic_second_extrema",
          all(s.simplify(s.diff(H,t,3).subs(t,z)) == 0 for z in critical)
          and all(s.simplify(s.diff(H,t,2).subs(t,z)**2-s.Rational(100,3)) == 0
                  for z in critical))
    check("quintic_second_reserve", F(100,3) < 6**2)
    check("cutoff_gradient_reserves", F(15,8)*20 < 40 and F(15,8)*40 < 80)
    # Radial Laplacian in dimension two: H'' (2/R)^2 + H' (2/R)/r,
    # and r>=R/2; use |H''|<6 and |H'|<2, R^-2=400 M^2.
    check("cutoff_laplacian_reserve", (4*6+4*2)*400 < 16000)
    check("localized_graph_reserve", 4+2*80*81+16000 < 30000)
    # pi>3 and exp(.05)<20/19: even replacing sqrt(pi) by 1 suffices.
    check("point_evaluation_reserve", F(30000,2)*F(20,19) < 10**6)
    check("disk_proper_radius", F(20,19)/10 < F(1,9))
    check("row_exponential_reserve", 3*10**6*F(18,17) < 32*10**5)
    check("rank_two_hilbert_schmidt_reserve", 2*(32*10**5)**2 < 21*10**12)
    check("finite_background_vertex_reserve",
          1/(1-F(2,100))+2/(1-F(4,100)) < F(32,10))
    check("single_mode_constant", F(32,10)*21*10**12 == F(672,10)*10**12)

    N, q, C, E, d, a = s.symbols("N q C E d a", positive=True)
    n = s.symbols("n", integer=True, positive=True)
    count = s.summation(2*(q+2*n), (n,1,N))
    check("angular_count_exact", s.expand(count-2*N*(q+N+1)) == 0)
    gap = s.expand(4*N*(N+q)-count)
    check("angular_count_comparison", s.factor(gap-2*N*(N+q-1)) == 0)
    # N,q>=1 makes the displayed nonnegative gap explicit.
    check("angular_count_positive", s.expand(gap.subs({N:N+1,q:q+1}))
          == 2*N**2+2*N*q+4*N+2*q+2)

    P = a**4+2*a**3/d+6*a*a/d**2+12*a/d**3+12/d**4
    # Differentiation at the lower endpoint independently certifies the tail integral.
    check("stieltjes_tail_antiderivative",
          s.simplify(s.diff(s.exp(-a*d)*P,a)
                     +(d*a**4-2*a**3)*s.exp(-a*d)) == 0)
    check("tail_monotonicity_range", 24*18 > 2)
    check("normalized_tower_constant", F(672,10)*10**12*8 == F(5376,10)*10**12)
    p24 = F(18**4)+F(2*18**3,24)+F(6*18**2,24**2)+F(12*18,24**3)+F(12,24**4)
    e_lower = sum((F(1,math.factorial(j)) for j in range(11)), F(0))
    coefficient_upper = F(5376,10)*10**12*F(85,10**34)*p24/e_lower**432
    check("exact_exponential_tower_enclosure", coefficient_upper < F(2,10**200),
          {"rational_upper_scientific": f"{float(coefficient_upper):.16e}",
           "method": "(sum_{j=0}^10 1/j!)^-432, all comparisons rational"})
    check("physical_cold_separation", F(95,100)*(60-1) > 56)


def tr(A):
    return sum(A[i,i] for i in range(A.rows))


def diag(values):
    return mp.diag(list(values)*2)


def finite_case(M, scale):
    k = mp.matrix([[1,2,0,0],[-1,1,3,0],[0,-2,1,2],[0,0,-1,2]])*scale
    D = mp.zeros(8)
    for i in range(4):
        for j in range(4):
            D[i,j+4] = k[i,j]
            D[j+4,i] = -k[i,j]
    J = mp.diag([1]*4+[-1]*4)
    base = [mp.mpf(x) for x in (".003","-.004",".002","-.002")]
    u = diag([mp.mpf(x) for x in (".7","-.4","0","0")])
    w = diag([mp.mpf(x) for x in ("-.2",".6","0","0")])
    v = diag([mp.mpf(x) for x in ("0","0",".6","-.3")])
    px = mp.diag([1,1,0,0]*2)
    py = mp.eye(8)-px

    def mass(aa=0,bb=0,cc=0):
        return mp.diag([M*mp.exp(base[i%4]+aa*u[i,i]+bb*v[i,i]+cc*w[i,i])
                        for i in range(8)])

    def action(aa,bb,cc):
        return -mp.log(mp.det(D+mass(aa,bb,cc)))

    def resolvent(bb=0):
        m = mass(bb=bb)
        return (D+m)**-1, m

    def vertex(bb=0):
        R,m = resolvent(bb)
        seagull = tr(R*m*u*w*R*m*v)
        triangle1 = tr(R*m*u*R*m*v*R*m*w)
        triangle2 = tr(R*m*u*R*m*w*R*m*v)
        return seagull-triangle1-triangle2

    def hessian(bb=0):
        R,m = resolvent(bb)
        return -tr(R*m*u*w)+tr(R*m*u*R*m*w)

    prefix=f"matrix_M{M}_K{scale}"
    R,m = resolvent()
    check(prefix+"_operator_algebra", mp.norm(D.T+D) == 0 and mp.norm(J*D*J+D) == 0)
    check(prefix+"_disjoint_products", mp.norm(u*v) == 0 and mp.norm(v*w) == 0
          and mp.norm(u*w) > 0)
    check(prefix+"_adjoint_resolvent", mp.norm(R.T-J*R*J) < mp.mpf("1e-60"))
    cross = px*R*py
    reverse = py*R*px
    check(prefix+"_cross_HS_symmetry",
          abs(mp.norm(cross)**2-mp.norm(reverse)**2) < mp.mpf("1e-60"))
    opnorm = np.linalg.norm(np.array(R.tolist(), dtype=float), ord=2)
    check(prefix+"_resolvent_norm", opnorm <= float(mp.exp(mp.mpf(".01"))/M))
    target = vertex()
    errors=[]
    estimates=[]
    for hstr in (".0001",".00005"):
        h=mp.mpf(hstr)
        difference=mp.mpf(0)
        for a in (-1,1):
            for b in (-1,1):
                for c in (-1,1):
                    difference += a*b*c*action(a*h,b*h,c*h)
        difference /= 8*h**3
        errors.append(abs(difference-target))
        estimates.append(difference)
    check(prefix+"_independent_logdet_difference", max(errors) < mp.mpf("1e-7"),
          {"trace":mp.nstr(target,35),"central_differences":[mp.nstr(e,35) for e in estimates],
           "errors":[mp.nstr(e,12) for e in errors]})
    check(prefix+"_difference_refinement",
          errors[1] < mp.mpf(".35")*errors[0] or max(errors) < mp.mpf("1e-45"))
    # Stay in the rho=.01 value ball along this nonconstant background path.
    endpoint=mp.mpf(".01")
    check(prefix+"_path_value_ball",
          max(abs(base[i]+endpoint*v[i,i]) for i in range(4)) < mp.mpf(".01"))
    integrated=mp.quad(vertex,[0,endpoint])
    change=hessian(endpoint)-hessian()
    check(prefix+"_finite_hessian_fundamental_theorem",
          abs(integrated-change) < mp.mpf("1e-30"),
          {"quadrature":mp.nstr(integrated,35),"hessian_change":mp.nstr(change,35),
           "error":mp.nstr(abs(integrated-change),12)})
    # Actual algebraic HS inequality, independently using the matrix blocks.
    prefactor=mp.exp(mp.mpf(".02"))+2*mp.exp(mp.mpf(".04"))
    rhs=prefactor*M*M*mp.norm(cross)**2*max(abs(u[i,i]) for i in range(8))*\
        max(abs(v[i,i]) for i in range(8))*max(abs(w[i,i]) for i in range(8))
    check(prefix+"_trace_HS_inequality", abs(target) <= rhs)
    NUMERICS.append({"M":M,"scale":scale,"vertex":mp.nstr(target,35),
                     "logdet_errors":[mp.nstr(e,12) for e in errors],
                     "hessian_change":mp.nstr(change,35),
                     "quadrature_error":mp.nstr(abs(integrated-change),12)})


def hashes():
    return {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    initial=hashes()
    error=None
    try:
        exact_checks()
        mp.mp.dps=70
        finite_case(18,6)
        finite_case(23,13)
    except Exception:
        error=traceback.format_exc()
        print(error,file=sys.stderr)
    final=hashes()
    ok=error is None and initial==final and all(r["passed"] for r in ROWS)
    record={
        "story":"S143","created_utc":datetime.now(timezone.utc).isoformat(),
        "status":"passed" if ok else "failed","python":sys.version,
        "numpy":np.__version__,"sympy":s.__version__,"mpmath":mp.__version__,
        "input_sha256_before":initial,"input_sha256_after":final,
        "unchanged_inputs":initial==final,"checks":ROWS,"numerics":NUMERICS,
        "exception":error,
        "scope":"ordinary continuum proof of separated massive response, exact "
                "constant checks and non-certifying finite-matrix diagnostics",
        "unresolved":["coincident curved finite-background response",
                      "general metric and LLL response","coupled exterior and signed source",
                      "physical stability, transport and preparation"],
    }
    if error:
        record["failed_checker_source_utf8"]=Path(__file__).read_text(encoding="utf-8")
        record["failed_note_source_utf8"]=(ROOT/INPUTS[1]).read_text(encoding="utf-8")
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
