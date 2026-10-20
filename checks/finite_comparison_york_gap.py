"""S167: quantitative global York map on the actual finite comparison.

Exact metric algebra and rational budgets; continuum claims remain ordinary proofs.
"""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "finite-comparison-york-gap"
OWN = ["checks/finite_comparison_york_gap.py", "notes/"+STEM+".md"]
PREDECESSORS = [("finite-comparison-scalar-resolvent", 297),
                ("stationary-global-tt-reduction", 44),
                ("stationary-curved-quantum-tt", 157)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "uniform actual finite-comparison York inverse and physical-L2 TT projection",
    "relative_graph_error": ".004", "gap_over_lambda": "1.47",
    "York_inverse_times_lambda": ".681", "TT_projector_difference": ".0041",
    "metric_derivatives_used": 1, "physical_density_unitaries_retained": True,
    "fixture_scope": "exact metric and finite rectangular weighted projection examples",
    "finite_comparison_TT_Hessian_positive": False, "sharp_logarithmic_York_map": False,
    "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
    "nonlinear_family_preserved": False, "signed_compatibility_evaluated": False,
    "formal_verification": False,
}


def check(name, ok, detail=""):
    assert name not in {row["name"] for row in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expr):
    expr = S.cancel(expr)
    check(name, expr == 0, expr)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_bytes())
    return data.get("runs", [data])[-1]


def inputs(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def paths():
    found = set(OWN)
    for stem, _ in PREDECESSORS:
        found.add("evidence/"+stem+".json")
        found.update(inputs(latest(stem)))
    return sorted(found)


def predecessors():
    expected = {}
    for stem, count in PREDECESSORS:
        run = latest(stem)
        check(stem+"_receipt", run["status"] == "passed" and len(run["checks"]) == count
              and all(row["passed"] for row in run["checks"]) and inputs(run) == inputs(run, True))
        conflicts = [p for p, h in inputs(run).items() if p in expected and expected[p] != h]
        check(stem+"_shared_inputs", not conflicts, conflicts)
        expected.update(inputs(run))
    for p, h in sorted(expected.items()):
        check("input_"+p, sha256((ROOT/p).read_bytes()).hexdigest() == h)
    for p in ["finite-comparison-scalar-resolvent", "stationary-global-tt-reduction",
              "stationary-curved-quantum-tt", "stationary-weighted-tt-inverse",
              "stationary-completed-metric-response-norm", "endpoint-preserving-compatible-projection"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


def zero_matrix(name, mat):
    check(name, all(S.simplify(x)==0 for x in mat), mat)


def physical_operator_algebra():
    y=S.Symbol('y',positive=True)
    a,r,w,s,v=[S.Function(n)(y) for n in ['a','r','w','s','v']]
    alpha=S.diff(a,y)/a
    rr=S.diff(r,y)/r
    dz=lambda f:S.exp(-w)*S.diff(f,y)
    aa,rad=S.exp(w)*a,S.exp(s)*r
    ae,re=dz(aa)/aa,dz(rad)/rad
    beta=S.exp(w+s)
    def lc(f,deriv,orb,radial):
        return S.Matrix([R(3,2)*deriv(f)-(orb+2*radial)*f/2,
                         -deriv(f)/2+(R(3,2)*orb-radial)*f,
                         -deriv(f)/2+(radial-orb/2)*f])
    dy=lambda f:S.diff(f,y)
    seed=lc(v,dy,alpha,rr)
    weight=S.diag(1,1,2)
    ww,ss=S.diff(w,y),S.diff(s,y)
    ee=S.Matrix([-2*ww-R(5,2)*ss,2*ww-ss/2,R(3,2)*ss])
    actual=beta*lc(v/beta,dz,ae,re)
    zero('physical_measure_ratio',aa*rad**2*S.exp(w)/(a*r*r)-beta**2)
    zero_matrix('actual_unitary_York_operator',actual-S.exp(-w)*(seed+ee*v))
    zero('exact_error_polynomial',(ee.T*weight*ee)[0]-8*ww**2-8*ww*ss-11*ss**2)
    zero('error_tracefree',ee[0]+ee[1]+2*ee[2])
    zero_matrix('common_coordinate_vector_cancellation',lc(S.exp(w)*v,dz,ae,re)
                -seed-(ww-ss)*v*S.Matrix([1,1,-1]))
    ell=S.Symbol('ell')
    uu=S.Function('u')(y)
    zero('radius_completion_cancels',S.diff(ell*uu+s,y)-S.diff(s,y)-ell*S.diff(uu,y))
    # Independent divergence and formal-adjoint factorization in physical radial measure.
    j=alpha+2*rr
    bb=S.Matrix([R(3,2),-R(1,2),-R(1,2)])
    cc=S.Matrix([-j/2,R(3,2)*alpha-rr,rr-alpha/2])
    divergence=dy(seed[0])+alpha*(seed[0]-seed[1])+2*rr*(seed[0]-seed[2])
    expected=R(3,2)*(S.diff(v,y,2)+j*S.diff(v,y))-(S.diff(j,y)/2+2*alpha**2+4*rr**2)*v
    zero('York_second_order_divergence',divergence-expected)
    adjoint_lie=-(bb.T*weight*seed.diff(y))[0]-j*(bb.T*weight*seed)[0]+(cc.T*weight*seed)[0]
    zero('factorization_sign_and_half',adjoint_lie/2+divergence)
    t1,t2=S.Function('t1')(y),S.Function('t2')(y)
    tt=S.Matrix([-t1-2*t2,t1,t2])
    divt=dy(tt[0])+alpha*(tt[0]-tt[1])+2*rr*(tt[0]-tt[2])
    adjt=-(bb.T*weight*tt.diff(y))[0]-j*(bb.T*weight*tt)[0]+(cc.T*weight*tt)[0]
    zero('tensor_adjoint_minus_twice_divergence',adjt+2*divt)


def bolt_parity_fixture():
    y=S.Symbol('y',positive=True)
    c,d=S.symbols('c d',real=True)
    # Independent regular warp fixture, not the physical finite comparison.
    a=S.sinh(y); r=2+y*y
    w=c*y*y; s=d*y*y
    beta=S.exp(w+s)
    dz=lambda f:S.exp(-w)*S.diff(f,y)
    aa,rr=S.exp(w)*a,S.exp(s)*r
    orb,rad=dz(aa)/aa,dz(rr)/rr
    Y=y/beta
    raw=S.Matrix([R(3,2)*dz(Y)-(orb+2*rad)*Y/2,
                  -dz(Y)/2+(R(3,2)*orb-rad)*Y,
                  -dz(Y)/2+(rad-orb/2)*Y])
    seed=S.Matrix([R(3,2)-(S.coth(y)+4*y/(2+y*y))*y/2,
                   -R(1,2)+(R(3,2)*S.coth(y)-2*y/(2+y*y))*y,
                   -R(1,2)+(2*y/(2+y*y)-S.coth(y)/2)*y])
    for i,expected in enumerate([1,1,-1]):
        check('fixture_regular_bolt_component_'+str(i),S.limit(beta*raw[i],y,0)==expected)
    expected=S.Matrix([-4*c-5*d,4*c-d,3*d])
    for i in range(3):
        # Cancel the hyperbolic/exponential terms algebraically before the
        # singular quotient. The first retained run's direct limit returned
        # -oo; this identity independently establishes the removable zero.
        difference=S.simplify((beta*raw[i]-S.exp(-w)*seed[i]).rewrite(S.exp))
        zero('fixture_exact_bolt_difference_'+str(i),difference-S.exp(-w)*expected[i]*y**2)
        zero('fixture_bolt_error_order2_'+str(i),S.limit(difference/y**2,y,0)-expected[i])
    zero('fixture_vector_component_odd',Y.subs(y,-y)+Y)
    zero('fixture_density_even',beta.subs(y,-y)-beta)


def exact_budgets():
    eps=R('3e-8'); delta=R('.004'); expbound=1/(1-eps)
    check('E_first_jet_square',8+8+11==27)
    check('sqrt27_upper',27<R('5.2')**2)
    check('seed_gap_square_root_lower',R('1.7')**2<R('2.98'))
    graph=expbound-1+expbound*R('5.2')*R('1e-3')/R('1.7')
    check('relative_graph_error',graph<delta,graph)
    check('strict_graph_gap',0<delta<1)
    check('finite_York_gap',R('1.49')*(1-delta)**2>R('1.47'))
    check('finite_York_inverse',1/R('1.47')<R('.681'))
    angle=delta/(1-delta)
    check('TT_projector_difference',angle<R('.0041'),angle)
    check('York_value_extraction',1/(R('2.98')*(1-delta)**2)<R('.584')**2)
    check('York_derivative_extraction',1/(3*(1-delta)**2)<R('.58')**2)
    diff=delta*(2-delta)/(1-delta)**2
    check('extraction_graph_difference',diff<R('.0081'),diff)
    check('extraction_value_difference',R('.0081')**2/R('2.98')<R('.0047')**2)
    invrel=(2*delta+delta*delta)/(1-delta)**2
    check('inverse_relative_form_difference',invrel<R('.0082'),invrel)
    check('inverse_absolute_difference',R('.0082')/R('1.49')<R('.0056'))
    for d in [R(0),R('.001'),delta,R('.1'),R('.5')]:
        check('inverse_spectral_lower_'+str(d),1/(1-d)**2-1 <= (2*d+d*d)/(1-d)**2)
        check('inverse_spectral_upper_'+str(d),1-1/(1+d)**2 <= (2*d+d*d)/(1-d)**2)


def rectangular_physical_fixture():
    aa=S.Matrix([[2,0],[0,2],[0,0],[0,0]])
    perturb=R('.001')*S.Matrix([[1,2],[-1,1],[2,-1],[1,3]])
    ab=aa+perturb
    identity=S.eye(4)
    def left_inverse(mat):
        return (mat.T*mat).inv()*mat.T
    aplus,bplus=left_inverse(aa),left_inverse(ab)
    pa,pb=aa*aplus,ab*bplus
    def frobenius_squared(mat):
        return S.trace(mat.T*mat)
    check('fixture_relative_graph_error',frobenius_squared(perturb)/4<R('.004')**2)
    for name,p in [('seed',pa),('finite',pb)]:
        zero_matrix('fixture_'+name+'_idempotent',p*p-p)
        zero_matrix('fixture_'+name+'_orthogonal',p-p.T)
    zero_matrix('fixture_projection_difference_split',pb-pa-pb*(identity-pa)+(identity-pb)*pa)
    zero_matrix('fixture_difference_outputs_orthogonal',((identity-pb)*pa).T*pb*(identity-pa))
    check('fixture_TT_projector_difference',frobenius_squared(pb-pa)<R('.0041')**2)
    zero_matrix('fixture_seed_extraction_identity',aplus*aa-S.eye(2))
    zero_matrix('fixture_finite_extraction_identity',bplus*ab-S.eye(2))
    zero_matrix('fixture_reconstructed_range',ab*bplus-pb)
    zero_matrix('fixture_extraction_difference_identity',ab*(bplus-aplus)-(pb-pa)+perturb*aplus)
    check('fixture_extraction_graph_bound',frobenius_squared(aa*(bplus-aplus))<R('.0081')**2)
    qa,qb=aa.T*aa/2,ab.T*ab/2
    shifted=qb-R('1.47')*S.eye(2)
    check('fixture_gap_Sylvester',shifted[0,0]>0 and shifted.det()>0)
    check('fixture_inverse_difference',frobenius_squared(qb.inv()-qa.inv())<R('.0056')**2)
    # Non-Euclidean physical pairings are part of the fixture, not optional.
    uv=S.diag(2,3); ut=S.diag(1,2,3,4)
    gv,gt=uv*uv,ut*ut
    physical=ut.inv()*ab*uv
    adj=gv.inv()*physical.T*gt
    pp=physical*(adj*physical).inv()*adj
    zero_matrix('fixture_physical_adjoint_pairing',gv*adj-physical.T*gt)
    zero_matrix('fixture_physical_weighted_selfadjoint',gt*pp-pp.T*gt)
    check('fixture_physical_projector_not_Euclidean_symmetric',pp!=pp.T)
    zero_matrix('fixture_unitary_tensor_projection',ut*pp*ut.inv()-pb)
    zero_matrix('fixture_unitary_vector_factorization',uv*(adj*physical/2)*uv.inv()-qb)
    zero_matrix('fixture_orthogonal_TT_range',physical.T*gt*(S.eye(4)-pp))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S167 finite comparison quantitative York map",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        physical_operator_algebra()
        bolt_parity_fixture()
        exact_budgets()
        rectangular_physical_fixture()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S167 finite comparison quantitative York map",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "finite_comparison_uniform_York_gap": good,
        "physical_L2_TT_projector_comparison": good,
        "weak_York_extraction_bounds": good, "regular_vector_domain_retained": good,
        "finite_comparison_TT_Hessian_positive": False, "sharp_logarithmic_York_map": False,
        "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
        "signed_compatibility_evaluated": False, "nonlinear_family_preserved": False,
        "useful_Newton_certificate": False, "S167_complete": False, "goal_complete": False,
        "unresolved": ["finite-comparison TT/mixed/contact forms", "sharp logarithmic mapping and candidate range",
                       "projection direction and cold/higher bounds", "nonlinear preservation and solve",
                       "signed compatibility", "physical requirements"],
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
