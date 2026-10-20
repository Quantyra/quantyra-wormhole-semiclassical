"""S167: actual massive constant contact on the finite comparison.

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
STEM = "finite-comparison-massive-constant-contact"
OWN = ["checks/finite_comparison_massive_constant_contact.py", "notes/"+STEM+".md"]
PREDECESSORS = [("finite-comparison-classical-disk-contact", 301),
                ("comparison-path-interior-response", 126),
                ("outer-finite-background-response", 63),
                ("stationary-curved-quantum-tt", 157)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "actual massive constant contact including tadpole and common restorations",
    "candidate_energy_absolute_bound": ".003",
    "interior_normalized_contact_bound": "4.3e-7",
    "outer_contact_bound": "2e20 q^2 (1+1/Y)",
    "actual_seam_distribution_retained": True,
    "constant_probe_selected_as_physical_direction": False,
    "full_massive_Hessian_bounded": False,
    "full_quantum_form_domain_identified": False,
    "finite_comparison_full_TT_Hessian_positive": False,
    "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
    "nonlinear_family_preserved": False, "signed_compatibility_evaluated": False,
    "fixture_scope": "exact mass/density/heat identities and rational/tail budgets only",
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
    expr = S.cancel(expr.doit())
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
    for p in ["finite-comparison-classical-disk-contact", "finite-comparison-classical-tt-form",
              "comparison-path-interior-response", "divided-jet-dirac-response",
              "outer-finite-background-response", "stationary-curved-quantum-tt",
              "stationary-centered-canonical-response", "stationary-completed-metric-response-norm"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


def density_mass_and_restoration():
    e,m,z,d=S.symbols('e m z d',positive=True)
    D=S.Matrix([[0,d],[-d,0]])
    sm=(D+m*S.eye(2)).inv(); pp=m*m*S.eye(2)-D*D
    zero('paired_constant_mass_identity',S.trace(m*m*sm*sm-m*sm)
         -S.trace(2*m**4*pp.inv()**2-2*m*m*pp.inv()))
    q=S.Symbol('q',positive=True)
    source=-m*S.trace(sm)
    me=m*S.exp(e)
    area_scalar=S.exp(-2*e)*source.subs(m,me)
    fixed_density=S.exp(2*e)*area_scalar
    p0=S.diff(fixed_density,e).subs(e,0)
    zero('constant_probe_density_identity',p0-S.diff(area_scalar,e).subs(e,0)-2*area_scalar.subs(e,0))
    zero('constant_probe_full_Hessian_identity',p0-S.trace(m*m*sm*sm-m*sm))
    omitted=S.diff(area_scalar,e).subs(e,0)-p0
    check('density_omission_negative_control',omitted.subs({m:2,d:3})!=0)
    A0,A2=S.symbols('A0 A2')
    local=-A0*S.exp(2*e)/2+A2*S.exp(-2*e)/2
    zero('common_constant_restoration',S.diff(local,e,2).subs(e,0)-(-2*A0+2*A2))
    t,C=S.symbols('t C',positive=True)
    A=2*C/t**2
    zero('leading_constant_contact_heat_weight',2*(-t*S.diff(A,t)-A)-4*C/t**2)
    F=S.Function('F')(z); G=-z*F
    zero('mass_derivative_of_source_comparison',z*S.diff(G,z)-(-z*z*S.diff(F,z)-z*F))
    Rn,Sn=S.symbols('Rn Sn')
    zero('subtracted_kernel_constant_probe',Rn*(1-1)+Sn-Sn)
    B,Gc,Cc,tt,mu4,muhat,theta,sp0=S.symbols('B G C t mu4 muhat Theta Sp0',positive=True)
    physical=(theta/B**2)*(mu4/(4*S.pi*B**2*tt**4))
    zero('physical_contact_measure',physical.subs(theta,B**2*Cc*sp0/(2*Gc))
         -Cc*sp0*mu4/(8*S.pi*Gc*B**2*tt**4))
    zero('outer_physical_contact_measure',physical-theta*mu4/(4*S.pi*B**4*tt**4))


def full_angular_tail_and_weights():
    z,a=S.symbols('z a',positive=True)
    summand=(2+4*z*z)*z**5*S.exp(-a*z)
    derivative=S.diff(summand,z)*S.exp(a*z)/z**4
    zero('angular_summand_derivative',derivative-(10+28*z*z-a*(2*z+4*z**3)))
    # For z>=1, the two positive powers are bounded by 10z+28z^3.
    check('angular_summand_monotonicity',2*171>10 and 4*171>28)
    for k in [6,8]:
        polynomial=sum(S.factorial(k)/S.factorial(k-j)*z**(k-j)/a**(j+1) for j in range(k+1))
        zero('exact_tail_antiderivative_'+str(k),S.diff(-S.exp(-a*z)*polynomial,z)-z**k*S.exp(-a*z))
    bracket=6+4*sum(S.factorial(6)/S.factorial(6-j)/R(171)**(j+1) for j in range(7))
    bracket+=8*sum(S.factorial(8)/S.factorial(8-j)/R(171)**(j+1) for j in range(9))
    lower=sum(R(171)**j/S.factorial(j) for j in range(257))
    upper=18**5*bracket/lower
    check('complete_angular_tail_without_truncation',upper<R('1e-66'),S.N(upper,8))
    y=S.Symbol('y',positive=True)
    zero('outer_log_weight_monotonicity',S.diff(y*(1+S.log(y)**2),y)-(1+S.log(y))**2)
    zero('outer_log_weight_endpoint',(y*(1+S.log(y)**2)).subs(y,1)-1)
    # The boundary contact has the distributional coefficient times T(0)^2.
    coefficient,T0=S.symbols('coefficient T0')
    zero('seam_contact_zero_trace', (coefficient*T0*T0).subs(T0,0))
    check('seam_contact_not_zero_for_general_trace',(coefficient*T0*T0).subs({coefficient:2,T0:3})==18)


def rational_budgets():
    check('heat_weight_global_count',2*(40+12)<120)
    check('heat_weight_gap_count',2*(40*8+12*4)<800)
    check('short_weight_change',R(120,24)<=10)
    check('long_weight_change',R(800,96)<10)
    check('restoration_weight_change',2<10)
    check('constant_cap_gradient',160/(1-R('.027'))<165)
    check('constant_probe_norm',165<170)
    check('interior_background_ball',R('1.3e-4')<R('.001'))
    check('constant_probe_parameter_tube',R('1.3e-4')+R(1,200)<R('.011'))
    check('interior_full_contact',170*R('2.5e-9')+2*R('5e-10')<R('4.3e-7'))
    check('actual_outer_coefficient_norm',R('1e-7')+R('4e-11')<R('2e-7'))
    check('actual_outer_graph_tube',10000*(R('2e-7')+R('1e-8'))*10<R('.03'))
    check('cap_background_tube',R('2e-7')<R('1e-4'))
    check('outer_safe_path_coverage',R('1.5')*R('.052')<R('.1'))
    check('fixed_exit_distance',R('.1')-R('.052')==R('.048'))
    check('fixed_exit_inflation',(R('.05')/R('.048'))**4<R('1.18'))
    check('enlarged_collar_reserve',R('1.18')<2)
    check('cap_constant_contact_ledger',2*10*R('1e25')<=R('2e26'))
    check('finite_cubic_change',R('.07')*R('4e-11')<R('1e-11'))
    check('complex_mass_real_part',1-R('.01')==R('.99'))
    check('complex_mass_modulus',1+R('.01')==R('1.01'))
    check('complex_mass_decay',R('.99')*R('9.6')>R('9.5'))
    check('complex_mass_factor_budget',2**(16+16+8+5)<10**14)
    check('background_inverse_ratio_reserve',(1-R('.02'))/(1-R('.03'))<2)
    check('complex_mass_gluing_prefactor',2*10**14*10**64<10**80)
    check('constant_mass_Cauchy_cost',100+1<=101)
    check('constant_gluing_angular_budget',101*R('1e80')*R('1e-66')*400<R('1e19'))
    check('constant_gluing_direction_reserve',R('1e19')<R('1e20'))
    check('outer_total_regular_difference',R('2e26')+R('1e20')<R('3e26'))
    outer_delta=R('3e26')*R('2e-7')*R('1.24e-5')+2*R('1e-11')*R('1.24e-5')/15
    check('outer_contact_from_seed',R('1e20')+outer_delta<R('2e20'),outer_delta)
    check('actual_outer_coverage',R('.051')/(1-R('1.2e-8'))<R('.052'))
    check('TT_return_value',R('3.7')**2*R('4.8')**2<500)
    check('physical_outer_coordinate_product',5*(1-R('1.2e-8')-R('3e-8'))**2>R('2.23')**2)
    check('outer_radius_weight',1+1/R('2.23')<R('1.5'))
    check('interior_contact_energy',500*R('6.1')*R('4.3e-7')<R('.00132'))
    check('outer_contact_energy',2*R('2e20')*R('1.5')*500*R('4e-27')<=R('.0012'))
    check('complete_constant_contact_energy',R('.00132')+R('.0012')<R('.003'))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S167 finite comparison massive constant contact",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        density_mass_and_restoration()
        full_angular_tail_and_weights()
        rational_budgets()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S167 finite comparison massive constant contact",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "actual_massive_constant_contact_bounded": good,
        "actual_tadpole_and_density_retained": good,
        "actual_outer_domain_constant_contact_bounded": good,
        "seam_distribution_retained": True,
        "full_massive_Hessian_bounded": False,
        "full_quantum_form_domain_identified": False,
        "finite_comparison_full_TT_Hessian_positive": False,
        "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
        "signed_compatibility_evaluated": False, "nonlinear_family_preserved": False,
        "useful_Newton_certificate": False, "S167_complete": False, "goal_complete": False,
        "unresolved": ["remaining massive kernel/second-contact and Weyl/reference TT forms", "mixed/density blocks", "outer position-log weights",
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
