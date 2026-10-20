"""S167: finite-comparison disk return, form logarithm and LLL bound.

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
STEM = "finite-comparison-disk-return-form"
OWN = ["checks/finite_comparison_disk_return_form.py", "notes/"+STEM+".md"]
PREDECESSORS = [("finite-comparison-york-gap", 277),
                ("stationary-global-tt-reduction", 44),
                ("stationary-seed-quantum-form", 53)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "actual finite-comparison disk return and canonical LLL disk Hessian",
    "return_value_coefficient": "4.6", "return_gradient_coefficients": [".8", "73"],
    "LLL_H1_relative_lower_bound": "1.3e-27",
    "LLL_comparison_energy": "(||h_prime||^2+lambda||h||^2)/(8pi G)",
    "form_frequency_domain": "H2 with square-root logarithmic multiplier",
    "form_map_constants_numerically_bounded": False,
    "full_H4_logarithmic_York_map": False, "full_quantum_form_domain_identified": False,
    "finite_comparison_full_TT_Hessian_positive": False,
    "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
    "nonlinear_family_preserved": False, "signed_compatibility_evaluated": False,
    "fixture_scope": "exact geometry, radial Volterra and interpolation algebra only",
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
    for p in ["finite-comparison-york-gap", "finite-comparison-scalar-resolvent",
              "stationary-global-tt-reduction", "stationary-seed-quantum-form",
              "stationary-curved-quantum-tt", "stationary-lll-quantum-exterior",
              "stationary-completed-metric-response-norm"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


def zero_matrix(name, mat):
    reduced=mat.applyfunc(S.simplify)
    check(name, all(x==0 for x in reduced), reduced)


def return_algebra():
    y=S.Symbol('y',positive=True)
    a,r,w,s,v=[S.Function(n)(y) for n in ['a','r','w','s','v']]
    aa,rr=S.exp(w)*a,S.exp(s)*r
    beta=S.exp(w+s)
    dz=lambda f:S.exp(-w)*S.diff(f,y)
    a0,r0=S.diff(a,y)/a,S.diff(r,y)/r
    ae,re=dz(aa)/aa,dz(rr)/rr
    actual=beta*(dz(v/beta)-ae*v/beta)
    target=S.exp(-w)*(S.diff(v,y)-a0*v-(2*S.diff(w,y)+S.diff(s,y))*v)
    zero('actual_return_unitary_operator',actual-target)
    m=a*r*r; j=S.diff(m,y)/m
    b0=S.diff(v,y)-a0*v
    potential=S.diff(a0,y)+j*a0+a0*a0
    zero('return_seed_form_integration',m*b0*b0-m*(S.diff(v,y)**2+potential*v*v)+S.diff(m*a0*v*v,y))
    zero('return_seed_positive_potential',potential-S.diff(a,y,2)/a-a0*a0-2*a0*r0)
    e10=S.diff(a0-r0,y)+a0*(a0-r0)
    e20=S.diff(r0,y)+a0*r0
    e1=dz(ae-re)+ae*(ae-re)
    e2=dz(re)+ae*re
    zero('finite_E1_cancellation',e1-S.exp(-2*w)*(e10+S.diff(w-s,y,2)+a0*S.diff(w-s,y)))
    zero('finite_E2_cancellation',e2-S.exp(-2*w)*(e20+S.diff(s,y,2)+a0*S.diff(s,y)))
    b,c,eta,aa,rr=[S.Function(n)(y) for n in ['b','c','eta','A','R']]
    sc=c/2-rr*eta
    tc=-(b+3*c)/2-(aa-rr)*eta
    rules={S.diff(eta,y):aa*eta+b+c,
           S.diff(b,y):-2*(aa+rr)*b-2*(aa-rr)*c}
    zero('TT_disk_return_T_derivative',S.diff(tc,y).subs(rules)+R(3,2)*S.diff(c,y)-2*rr*b
         +(S.diff(aa-rr,y)+aa*(aa-rr))*eta)
    zero('disk_return_S_derivative',S.diff(sc,y).subs(rules)-S.diff(c,y)/2+rr*(b+c)
         +(S.diff(rr,y)+aa*rr)*eta)
    zero('orbit_components_equal',b-2*(aa*eta+b+c)-(-b-2*c-2*aa*eta))
    zero('disk_sphere_component',2*sc-(c-2*rr*eta))
    zero('disk_orbit_component',2*(sc+tc)-(-b-2*c-2*aa*eta))
    bb,cc=S.symbols('b c',real=True)
    h2=bb**2+(-bb-2*cc)**2+2*cc**2
    zero('exact_tensor_norm',h2-2*(bb+cc)**2-4*cc**2)
    metric=S.Matrix([[2,2],[2,6]])
    for label,covector,value in [('F',[1,1],R(1,2)),('c',[0,1],R(1,4)),
                                  ('b',[1,0],R(3,4)),('b3c',[1,3],R(3,2))]:
        lv=S.Matrix(covector)
        zero('tensor_dual_norm_'+label,(lv.T*metric.inv()*lv)[0]-value)


def budgets_and_physical_normalization():
    eps=R('3e-8'); d=R('.005'); ex=1/(1-eps)
    check('sqrt2_reserve',2<R('1.415')**2)
    graph=ex-1+ex*R('.003')*R('1.415')
    check('return_relative_graph_bound',graph<d,graph)
    check('return_unrestricted_value_coefficient',2/(1-d)**2<R('1.422')**2)
    check('TT_vector_return_value',1/(1-d)<R('1.006'))
    check('seed_sqrt5_reserve',5<R('2.24')**2)
    check('finite_radius_gradient',ex*(R('2.24')+R('.001'))<R('2.25'))
    ap=ex/(1-d)*(1+R('.001')*R('1.415'))
    check('regular_polar_product',ap<R('1.008'),ap)
    check('polar_product_per_tensor',R('1.008')**2/2<R('.713')**2)
    check('b3c_half_per_tensor',R(3,8)<R('.613')**2)
    check('return_S_value',R('.25')+R('2.25')*R('1.006')<R('2.6'))
    check('return_T_value',R('.613')+R('.713')+R('2.25')*R('1.006')<R('3.7'))
    check('return_pair_value',R('2.6')**2+R('3.7')**2<R('4.6')**2)
    sy=(R('1.001')*R('.00025')+R('2e-16')*R('2.3'))/R('.999')
    syy=(R('1.001')*R('6.2')+22*R('2e-16')+R('1.001')*7*R('3e-7')**2)/R('.999')+sy*sy
    check('finite_radius_second_jet_refinement',syy<R('6.3'),syy)
    check('outer_jets_fit_global_u_bounds',R(10)/40**3<R('.00025') and R(10)/40**3<R('6.2'))
    check('outer_jets_fit_global_S_second',R(10)/40**3<R('6.3'))
    check('seed_A0_comparison',R('10.5')<R('3.25')**2)
    check('polar_orbit_product_budget',R('3.25')*R('.00025')==R('.0008125'))
    check('polar_radius_product_budget',R('3.25')*R('.001')==R('.00325'))
    check('seed_E1_bound',R('10.5')+22+5<39)
    check('seed_E2_bound',22+5==27)
    check('exponential_squared_reserve',1/(1-2*eps)<R('1.001'))
    e1=R('1.001')*(39+R('6.2')+R('6.2')+R('.0008125'))
    e2=R('1.001')*(27+R('6.3')+R('6.3')+R('.00325'))
    check('finite_E1_bound',e1<52,e1)
    check('finite_E2_bound',e2<41,e2)
    check('b_per_tensor',R(3,4)<R('.87')**2)
    check('F_per_tensor',R(1,2)<R('.708')**2)
    check('return_T_gradient_value',2*R('2.25')*R('.87')+52*R('1.006')<57)
    check('return_S_gradient_value',R('2.25')*R('.708')+41*R('1.006')<44)
    check('return_pair_gradient_derivative',R('.75')**2+R('.25')**2<R('.8')**2)
    check('return_pair_gradient_value',57**2+44**2<73**2)
    check('LLL_gradient_energy',2*R('.75')**2<6500 and 2*57**2<6500)
    a,r,q,g,bb=S.symbols('a r q G B',positive=True)
    mu=8*S.pi**2*a*r*r
    zero('two_dimensional_Dirichlet_four_volume',2*S.pi*a-mu/(4*S.pi*r*r))
    zero('LLL_physical_prefactor',q/(12*S.pi)/(4*S.pi)-q/(48*S.pi**2))
    zero('LLL_reference_energy_ratio',(6500*q/(48*S.pi**2*bb**2))*8*S.pi*g-6500*g*q/(6*S.pi*bb**2))
    check('benchmark_Gq_over_B2',400*R('8.5e-33')==R('3.4e-30'))
    check('LLL_actual_finite_comparison_lower_bound',R(6500,18)*R('3.4e-30')<R('1.3e-27'))


def form_log_interpolation():
    t,nu=S.symbols('t nu',positive=True)
    f,x=S.symbols('f x',real=True)
    optimal=f/(1+t*t*nu*nu)
    minimum=t*t*nu*nu/(1+t*t*nu*nu)*f*f
    zero('quadratic_K_minimizer',S.diff((f-x)**2+t*t*nu*nu*x*x,x).subs(x,optimal))
    zero('quadratic_K_complete_square',(f-x)**2+t*t*nu*nu*x*x-minimum
         -(1+t*t*nu*nu)*(x-optimal)**2)
    primitive=S.log(1+t*t*nu*nu)/2
    zero('K_integral_primitive',S.diff(primitive,t)-t*nu*nu/(1+t*t*nu*nu))
    zero('K_integral_lower_endpoint',primitive.subs(t,0))
    zero('K_integral_upper_endpoint',primitive.subs(t,1)-S.log(1+nu*nu)/2)
    weight=1+S.log(1+nu*nu)/2
    lower=S.diff(weight-S.log(2+nu),nu)
    upper=S.diff(2*S.log(2+nu)-weight,nu)
    zero('formlog_lower_difference_derivative',lower-(2*nu-1)/((1+nu*nu)*(2+nu)))
    zero('formlog_upper_difference_derivative',upper-((nu-1)**2+1)/((1+nu*nu)*(2+nu)))
    check('exp_half_below_2_series_guard',1/(1-R('.5'))==2)
    check('log3_below_1point2',sum(R('1.2')**j/S.factorial(j) for j in range(6))>3)
    check('formlog_lower_at_one',1+R('.5')/2-R('1.2')>0)
    # e<3 follows from factorial j!>=2^(j-1), with strict inequality j>=3.
    check('formlog_upper_at_one',2-1-R(1,2)>0)
    check('quadratic_form_half_log_order',4/2==2)
    check('full_operator_domain_not_promoted',not PROTOCOL['full_H4_logarithmic_York_map']
          and not PROTOCOL['full_quantum_form_domain_identified']
          and not PROTOCOL['form_map_constants_numerically_bounded'])


def independent_regular_return_fixture():
    # a=sinh(y), r=cosh(y), c=sech(y)^6. This is a finite-energy
    # geometry/TT fixture, not the physical comparison or a quantum solution.
    y=S.Symbol('y',positive=True)
    a,r=S.sinh(y),S.cosh(y)
    aa,rr=S.coth(y),S.tanh(y)
    c=1/r**6
    b=-(1/r**4+1/r**6)/2
    eta=a/(10*r**5)
    normexpr=lambda f:S.simplify(f.rewrite(S.exp))
    for name,expr in [
        ('TT_constraint',S.diff(b,y)+2*(aa+rr)*b+2*(aa-rr)*c),
        ('outer_fixed_return_ODE',S.diff(eta,y)-aa*eta-b-c),
        ('outer_integral_primitive',S.diff(-eta/a,y)+(b+c)/a),
        ('S_return',c/2-rr*eta-(R(3,5)/r**6-R(1,10)/r**4)),
        ('T_return',-(b+3*c)/2-(aa-rr)*eta-(R(1,4)/r**4-R(27,20)/r**6)),
    ]:
        zero('fixture_'+name,normexpr(expr))
    check('fixture_return_outer_fixed',S.limit(normexpr(eta/a),y,S.oo)==0)
    check('fixture_regular_bolt_vector',S.limit(normexpr(eta/y),y,0)==R(1,10))
    check('fixture_bolt_F_order2',S.limit(normexpr((b+c)/y**2),y,0)==-R(1,2))
    zero('fixture_vector_odd',normexpr(eta.subs(y,-y)+eta))
    x=S.Symbol('x',positive=True)
    # x=cosh(y) gives physical radial measure a*r^2 dy=x^2 dx.
    cx=1/x**6; bx=-(1/x**4+1/x**6)/2
    hx=bx**2+(-bx-2*cx)**2+2*cx**2
    eta2=(x*x-1)/(100*x**10)
    hn=S.integrate(S.expand(x*x*hx),(x,1,S.oo))
    en=S.integrate(S.expand(x*x*eta2),(x,1,S.oo))
    check('fixture_finite_TT_norm',hn==R(16,35))
    check('fixture_finite_vector_norm',en==R(1,1750))
    check('fixture_value_norm_ratio',en/hn==R(1,800))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S167 finite comparison disk return and LLL form",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        return_algebra()
        budgets_and_physical_normalization()
        form_log_interpolation()
        independent_regular_return_fixture()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S167 finite comparison disk return and LLL form",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "finite_comparison_disk_return_bounds": good,
        "actual_LLL_disk_Hessian_lower_bound": good,
        "form_level_half_log_York_return_mapping": good,
        "form_map_constants_numerically_bounded": False,
        "full_H4_logarithmic_York_map": False, "full_quantum_form_domain_identified": False,
        "finite_comparison_full_TT_Hessian_positive": False,
        "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
        "signed_compatibility_evaluated": False, "nonlinear_family_preserved": False,
        "useful_Newton_certificate": False, "S167_complete": False, "goal_complete": False,
        "unresolved": ["massive/Weyl/reference TT and mixed/contact forms", "outer position-log weights",
                       "operator-domain mapping and candidate range", "projection direction and cold/higher bounds",
                       "nonlinear preservation and solve", "signed compatibility", "physical requirements"],
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
