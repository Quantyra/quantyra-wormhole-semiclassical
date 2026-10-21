"""S167: full massive form compared with the actual geodesic model.

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
STEM = "finite-comparison-massive-geodesic-form"
OWN = ["checks/finite_comparison_massive_geodesic_form.py", "notes/"+STEM+".md"]
PREDECESSORS = [("finite-comparison-massive-constant-contact", 286),
                ("stationary-curved-quantum-tt", 157),
                ("finite-comparison-classical-disk-contact", 301)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "entire actual massive disk Hessian versus the full geodesic spectral model",
    "comparison_bound": ".001 C J2 + .005 E_ell",
    "actual_massive_second_contact_reserve": "3e13 C",
    "global_row_subtracted_remainder": "3e25 C J0",
    "full_angular_tower_and_common_fourth_restoration": True,
    "bounded_tail_only_after_polynomial_allocation": True,
    "full_quantum_form_domain_identified": False,
    "finite_comparison_full_TT_Hessian_positive": False,
    "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
    "nonlinear_family_preserved": False, "signed_compatibility_evaluated": False,
    "fixture_scope": "exact finite algebra, moments and conservative rational budgets",
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
    for p in ["finite-comparison-massive-constant-contact", "finite-comparison-classical-disk-contact", "finite-comparison-classical-tt-form",
              "comparison-path-interior-response", "divided-jet-dirac-response",
              "outer-finite-background-response", "stationary-curved-quantum-tt",
              "stationary-centered-canonical-response", "stationary-completed-metric-response-norm"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


def geometry_and_kernel_algebra():
    ell,u,r0,dr,lap=S.symbols('ell u R0 dR lap_u')
    finite=S.exp(-2*ell*u)*(r0-2*ell*lap)
    zero('exact_finite_conformal_curvature',finite.subs(lap,-(dr+2*u*r0)/2)
         -S.exp(-2*ell*u)*(r0+ell*(dr+2*u*r0)))
    theta,a,b,d=S.symbols('theta a b d',real=True)
    contraction=a*S.cos(theta)**2+2*b*S.cos(theta)*S.sin(theta)+d*S.sin(theta)**2
    zero('angular_Hessian_identity',S.integrate(S.expand_trig(contraction**2),(theta,0,2*S.pi))
         -S.pi*((a+d)**2+2*(a*a+2*b*b+d*d))/4)
    j2,rg,cf=S.symbols('J2 curvature_gradient C_Feff')
    angular=S.pi*(3*j2-rg)/4
    zero('weak_Bochner_fourth_restoration',-cf*angular/(90*S.pi**2)
         -((-cf/4)*j2/(30*S.pi)-(-cf/4)*rg/(90*S.pi)))
    # Covariant inverse kernel: operator factors and source measure both matter.
    sx,sy,K=S.symbols('sx sy K',positive=True)
    op_inverse=sx**R(-1,2)*K*sy**R(3,2)
    zero('covariant_kernel_source_measure',op_inverse/sy**2-sx**R(-1,2)*K*sy**R(-1,2))
    check('positive_half_density_negative_control',
          (op_inverse/sy**2-sx**R(1,2)*K*sy**R(1,2)).subs({sx:4,sy:9,K:1})!=0)
    # Noncommuting resolvent fixture with [D,chi] commuting with chi.
    D=S.diag(2,3); chi=S.Matrix([[1,1],[0,1]]); B=(D+5*S.eye(2)).inv()
    Q=D*chi-chi*D
    zero('commutator_multiplier_fixture',(Q*chi-chi*Q).norm()**2)
    zero('two_sided_resolvent_fixture',
         (chi*B*chi-(B*chi**2+B*Q*chi*B-B*Q*B*Q*B)).norm()**2)
    zero('commutator_symmetric_identity_fixture',
         (B*Q*chi*B-(chi**2*B-B*chi**2)/2).norm()**2)
    x,y,k=S.symbols('x y k')
    raw=2*k*x*y; row=k*(x*x+y*y)
    zero('row_sum_subtracted_form',raw-row+k*(x-y)**2)
    check('raw_row_double_count_negative_control',(raw-(raw-row)).subs({k:2,x:3,y:4})!=0)


def spectral_moments_and_contacts():
    def kk(mu,nu):
        mu,nu=S.Integer(mu),S.Integer(nu)
        return S.sqrt(S.pi)*S.gamma(mu/2)*S.gamma(mu/2+nu)*S.gamma(mu/2-nu)/(4*S.gamma((mu+1)/2))
    zero('Bessel_second_moment',(kk(4,0)-kk(4,1))/(2*S.pi**2)+1/(6*S.pi**2))
    zero('Bessel_fourth_moment',(kk(6,0)-kk(6,1))/(2*S.pi**2)+4/(15*S.pi**2))
    zero('Bessel_fifth_absolute_moment',(kk(7,1)-kk(7,0))/(2*S.pi**2)-R(225,4096))
    z=S.symbols('z',real=True)
    g=1-z*z/2+z**4/24-S.cos(z)
    zero('positive_spectral_remainder_fourth_derivative',S.diff(g,z,4)-(1-S.cos(z)))
    for j in range(4): zero('positive_remainder_initial_jet_'+str(j),S.diff(g,z,j).subs(z,0))
    zero('finite_mode_fourth_cancellation',g-z**4/24-(1-z*z/2-S.cos(z)))
    zero('post_cancellation_bounded_tail_plus_second',g-z**4/24+z*z/2-(1-S.cos(z)))
    check('full_positive_tail_not_L2_bounded',S.limit(g/z**4,z,S.oo)==R(1,24))
    A,dq,cf=S.symbols('A Delta_q C_Feff')
    zero('a2_second_contact_difference',cf*A*dq/(20*S.pi)-cf*A*dq/(45*S.pi)-cf*A*dq/(36*S.pi))
    # At the normal center R=-2A and q=dq/2 times |z|^2.
    zero('curvature_gradient_operator_at_center',-cf*(-2*A)*dq/(72*S.pi)-cf*A*dq/(36*S.pi))
    r,aa,bb,cc=S.symbols('r a b c')
    square=(aa*r+bb*r**3+cc*r**4)**2
    zero('geodesic_quadratic_Taylor_remainder',S.expand(square-aa**2*r**2-2*aa*bb*r**4)
         -(2*aa*cc*r**5+bb**2*r**6+2*bb*cc*r**7+cc**2*r**8))
    x=S.symbols('x')
    tagged=(3-2*x)/(1-x)**2
    zero('tagged_Neumann_generating_function',tagged-(3/(1-x)+x/(1-x)**2))
    check('tagged_Neumann_finite_inflation',
          R(9,8)**3*tagged.subs(x,R('.135'))/tagged.subs(x,R('.13'))<2)


def normal_Born_and_global_budgets():
    check('finite_curvature_bound',(2+R('5.2e-8')+4*R('1.2e-8'))/(1-2*R('1.2e-8'))<R('2.001'))
    check('kappa4_upper',R(33,2)+R('1e-21')<R('17.5'))
    check('middle_directional_third',(600**2)*8<1700**2)
    check('cold_directional_third',(R('1e-6')*R('.00331')*6*100**3)**2*8<2000**2)
    delta_exp=1/(1-R('1e-6'))
    delta3=delta_exp*(2000+3*R('1e-4')*R('.05')+R('1e-4')**3)
    delta2=delta_exp*(R('.05')+R('1e-4')**2)
    product3=R('1.026')*delta3+3*R('1.664')*delta2
    product3+=3*R('212.992')*delta_exp*R('1e-4')+R('40894.464')*(delta_exp-1)
    check('finite_exponential_cubic_remainder',product3/6<400)
    check('outer_annular_cubic_remainder',2*R('.1')*64<=R('12.8'))
    check('finite_cubic_reserve',10000+400+R('12.8')<12000)
    check('finite_cap_sup',R('.026')*delta_exp+(delta_exp-1)<R('.027'))
    ratio=R('.027')*R('6.2')
    check('Born_ratio',ratio<R('.168'))
    born=9*R('.37')*R('6.2')/(1-ratio)/(1-R('.027'))
    born+=9*R('.37')/(1-R('.027'))**2
    check('covariant_endpoint_Born_budget',born<30)
    check('covariant_kernel_reserve',30<50)
    check('bubble_cross_coefficient',4*R('.48')*50<10000)
    check('bubble_square_coefficient',2*50**2<10000)
    check('complete_local_angular_envelope',10000*(256+R(2304,324))<10**7)
    check('local_eval_buffer',R('.19')/(2*18)<R(1,140))
    check('far_eval_buffer',R('.19')/(4*18)<R(1,280))
    check('local_exponential_buffer',1/(1-R('.19')/4)<R('1.05'))
    check('far_exponential_buffer',1/(1-R('.19')/8)<R('1.025'))
    check('local_row_reserve',60000*8*300*R('1.05')<R('1.6e8'))
    check('far_row_reserve',60000*64*550*R('1.025')<R('2.2e9'))
    check('local_cross_square_L2_budget',4*11+2*300**2/R(18)**2<1000)
    local=1000*R('1.6e8')**2*8/(R(2,280)**2)
    far=2*R('2.2e9')**2*8/(R(2,570)**2)
    check('local_full_angular_Schur',local<R('4.1e24'))
    check('far_full_angular_Schur',far<R('6.3e24'))
    check('combined_full_angular_Schur',R('4.1e24')+R('6.3e24')<R('1.1e25'))
    check('bounded_geodesic_tail',2*2*R(22,7)*18*128**4<R('2e12'))
    check('total_row_subtracted_remainder',2*R('1.1e25')+R('2e12')<R('3e25'))


def second_contact_and_form_budgets():
    check('one_Born_fifth',(64**3+12000+9**3)*1000<R('4e8'))
    check('two_Born_fifth',4*(12000+9**3)*12288000<R('7e11'))
    check('higher_Born_ratio',5*R('.027')==R('.135'))
    check('ultraviolet_second_contact',(R('4e8')+R('7e11'))*R(2,3)
          +R('2e13')*4*R('.00352')<R('1e13'))
    check('geodesic_second_derivative',16*R('1.02')**2<17)
    check('geodesic_third_derivative',7000*R('1.02')**3+2*16*17*R('1.02')<8000)
    fourth=6000000*R('1.02')**4+5*7000*17*R('1.02')**2+2*16*8000*R('1.02')+2*16*17**2
    check('geodesic_fourth_derivative',fourth<8000000)
    remainder=2*R(8000000,24)+R(1,128)*(5+R(8000000,24*128))**2
    check('geodesic_quadratic_remainder',remainder<800000)
    check('local_flow_second_contact',2*R(22,7)*800000*R(1,2)*R(2,3)<4000000)
    check('long_flow_quadratic_contact',2*R(22,7)*(R(18*128**4,1024)+36*128**2+10*10000)<R('4e9'))
    # log(12.5)<3 since a positive partial sum for exp(3) exceeds it.
    check('anchor_log_bound',sum(R(3)**j/S.factorial(j) for j in range(8))>R('12.5'))
    check('second_anchor_change',2*R(22,7)*R('4e8')*3<R('2e11'))
    check('all_unknown_second_contacts',R('1e13')+R('4e6')+R('4e9')+R('2e11')<R('3e13'))
    K=R('1e14'); eps=R('.001'); radius=R(1,400); bsq=eps/(4*K)
    check('form_split_inside_radius',bsq<radius**2)
    # Exact Young allocation: K*b*X*Z <= K X^2/4+K*b^2 Z^2.
    check('Young_fourth_loss',K*bsq*(1+R(1,8))<eps/3)
    logarg=radius**2/bsq
    check('form_log_argument',logarg==R('2.5e12'))
    check('form_log_upper',sum(R(30)**j/S.factorial(j) for j in range(71))>logarg)
    check('form_gradient_log_reserve',K*(R(1,4)+30)<R('1e16'))
    check('angular_fourth_loss',eps*R(22,7)/4<R('.001'))
    check('curvature_second_contact',70*R('2.001')/(72*3)<1)
    check('total_gradient_reserve',R('1e16')+R('3e13')+1+eps*R('2.001')*R(22,7)/12<R('2e16'))


def physical_energy_conversion():
    C,B,G,tt,lam,mu4,ec=S.symbols('C B G t lambda mu4 ec',positive=True)
    mu_bar=C*mu4/(4*S.pi*B**4*tt**4)
    zero('physical_J0_measure',C*mu_bar*8*S.pi*G/(lam*mu4)-2*G*C*C/(lam*B**4*tt**4))
    grad_scale=B*B*tt*tt/C
    zero('physical_J1_measure',C*mu_bar*grad_scale*8*S.pi*G/mu4-2*G*C/(B**2*tt**2))
    check('return_value_L2',R('3.7')**2*R('4.8')**2<500)
    derivative=R('1.5')*R('5.32')+R('4.5')*R('2.002')+R('39.001')*R('4.83')
    check('return_derivative_L2',derivative<210)
    check('return_derivative_squared',210**2<50000)
    check('physical_J0_energy',1000*R('6.1')*R('8.5e-33')<R('5.2e-29'))
    check('physical_J1_energy_majorant',100000*R('8.5e-33')==R('8.5e-28'))
    check('complete_massive_lower_order_cost',R('.003')+R('2e16')*R('8.5e-28')+R('3e25')*R('5.2e-29')<R('.005'))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S167 finite comparison massive geodesic form",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        geometry_and_kernel_algebra()
        spectral_moments_and_contacts()
        normal_Born_and_global_budgets()
        second_contact_and_form_budgets()
        physical_energy_conversion()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S167 finite comparison massive geodesic form",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "actual_massive_geodesic_form_comparison_bounded": good,
        "actual_massive_second_contact_bounded": good,
        "actual_global_gluing_rows_bounded": good,
        "full_angular_tower_and_common_fourth_restoration": True,
        "full_Weyl_outer_log_domain_identified": False,
        "full_quantum_form_domain_identified": False,
        "finite_comparison_full_TT_Hessian_positive": False,
        "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
        "signed_compatibility_evaluated": False, "nonlinear_family_preserved": False,
        "useful_Newton_certificate": False, "S167_complete": False, "goal_complete": False,
        "unresolved": ["Weyl/reference TT forms and outer logarithmic cancellation", "mixed/density blocks", "outer position-log weights",
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
