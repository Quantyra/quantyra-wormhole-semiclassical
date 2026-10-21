"""S167: physical TT contacts and fixed-density block at finite comparisons.

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
STEM = "finite-comparison-physical-tt-contacts"
OWN = ["checks/finite_comparison_physical_tt_contacts.py", "notes/"+STEM+".md"]
PREDECESSORS = [("finite-comparison-quantum-disk-tt", 378),
                ("comparison-path-complete-source-response", 305),
                ("general-background-newton-operator", 174),
                ("stationary-weighted-tt-inverse", 52)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "physical exponential/affine tensor TT forms and exact fixed-density disk block",
    "physical_tensor_TT_lower_bound": ".97 E_ell",
    "fixed_density_disk_TT_lower_bound": ".90 E_ell",
    "returned_density_contact_energy_dual_bound": ".025",
    "moving_volume_derivative_selfadjoint": False,
    "core": "smooth bolt-regular compact-outer TT and inherited explicit stronger completion",
    "independent_endpoint_and_action_parameters": True,
    "full_quantum_operator_domain_identified": False,
    "physical_mass_fast_domain_membership_proved": False,
    "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
    "nonlinear_family_preserved": False, "signed_compatibility_evaluated": False,
    "fixture_scope": "exact source normalization, contact/return/density algebra and rational budgets",
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
    for p in ["finite-comparison-quantum-disk-tt", "finite-comparison-classical-disk-contact",
              "finite-comparison-classical-tt-form", "comparison-path-complete-source-response",
              "comparison-path-interior-response", "comparison-path-cold-response",
              "comparison-path-outer-regular-response", "comparison-path-outer-recombination",
              "comparison-path-lll-reference-response", "general-background-newton-operator",
              "stationary-weighted-tt-inverse"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


def source_normalization():
    G,B,C,t,lam,rho,N,A=S.symbols('G B C t lambda rho N A',positive=True)
    radius=B*t
    rho2=rho*B**2*C/(2*G*radius**2)
    p2=(N-rho)*B**2*C/(2*G*radius**2)
    pS=A*B**2*C/(8*S.pi*G*radius**4)
    stress=S.Matrix([-rho2/(4*S.pi*radius**2),p2/(4*S.pi*radius**2),pS,pS])
    normalized=8*S.pi*G*stress/(2*lam)
    desired=C/(2*lam*B**2*t**4)*S.Matrix([-rho,N-rho,A,A])
    for j in range(4): zero('action_gradient_component_'+str(j),normalized[j]-desired[j])
    rr,nn=S.symbols('rho_bound N_bound',real=True)
    norm2=rr**2+(nn-rr)**2
    check('orbit_box_Frobenius_majorant',max(norm2.subs({rr:r,nn:n}) for r in [-1,1] for n in [-1,1])==5)
    check('missing_half_stress_negative_control',
          (2*normalized[0]-desired[0]).subs({G:1,B:2,C:3,t:4,lam:5,rho:6})!=0)
    # Calibrate the half-stress factor directly against the fixed reference action.
    X,z,h1,h2,hs=S.symbols('X z h1 h2 hs',positive=True)
    dv=S.Function('DeltaV')
    action=-S.exp(z*(h1+h2)/2)*dv(X*S.exp(z*hs))/(4*G)
    per_volume=S.diff(action,z).subs(z,0)/(4*S.pi*X)
    pA=dv(X)/(16*S.pi*G*X); ps=S.diff(dv(X),X)/(16*S.pi*G)
    zero('reference_calibrates_half_stress',per_volume+(pA*(h1+h2)+2*ps*hs)/2)


def source_values_and_physical_bound():
    check('cold_complete_orbit',R('1.1e-6')+R('1e-15')+R('1e-23')<R('1.2e-6'))
    check('middle_complete_orbit',R('1.7e-9')+R('1e-16')+R('1e-18')<R('1.8e-9'))
    check('outer_complete_orbit',R('3.1e-7')+R('1e-9')+R('1e-21')<R('3.2e-7'))
    cold2=(R('6.1')/2*R('.00331'))**2*(5*R('1.2e-6')**2+2*R('1.7e-6')**2)
    check('cold_physical_gradient',cold2<R('4e-8')**2)
    middle2=(R('6.1')/2)**2*(5*R('1.8e-9')**2+2*R('2.6e-9')**2)
    check('middle_physical_gradient',middle2<R('2e-8')**2)
    t=S.symbols('t',positive=True); L=1+S.log(t)
    zero('outer_weight_derivative',S.diff(L**2/t**4,t)-L*(2-4*L)/t**5)
    check('outer_radius_conversion',1/(1-4*R('3e-8'))<R('1.001'))
    check('outer_log_at_40',sum(R(4)**j/S.factorial(j) for j in range(13))>40)
    outer2=(R('6.1')/2*R('1.001')*25/40**4)**2*(5*R('3.2e-7')**2+2*R('4.7e-7')**2)
    check('outer_physical_gradient',outer2<R('3e-11')**2)
    check('uniform_remainder_gradient',R('2e-8')<R('4e-8') and R('3e-11')<R('4e-8'))
    check('uniform_full_gradient',R('1.8e-5')+R('4e-8')<R('1.81e-5'))
    # Exact cancellation of the cubic allocation, without bounding its pieces.
    P,Q,Pcub,T,W,Vvac,Jmass=S.symbols('P Q Pcub T W Vvac Jmass')
    U=W-Vvac-Jmass; J=Jmass-T
    zero('outer_density_recombination',(P-Q)+(Pcub+J+U)-(P+Pcub-Q-T+W-Vvac))
    zero('outer_null_recombination',(2*P-Q)+(2*Pcub+J+U)-(2*P+2*Pcub-Q-T+W-Vvac))
    check('missing_cubic_recombination_negative_control',(Jmass-T).subs({Jmass:3,T:2})!=0)


def gauge_and_chart_algebra():
    h,s,eta,hp,sp=S.symbols('h s eta hp sp')
    lie=2*eta*(hp-sp)+2*s*(2*h-2*s)
    Z=lie+(h-2*s)**2-h*h
    zero('combined_gauge_chart_acceleration',Z-2*eta*(hp-sp))
    check('dropping_chart_acceleration_negative_control',
          (lie-2*eta*(hp-sp)).subs({h:3,s:1,eta:2,hp:5,sp:4})!=0)
    H,grad,Liek,k2,h2=S.symbols('H grad Liek k2 h2')
    disk=H+grad*Liek+grad*k2; exp=H+grad*h2
    zero('disk_minus_tensor_exponential',disk-exp-grad*(Liek+k2-h2))
    zero('affine_from_exponential',exp-grad*h2-H)
    check('actual_combined_Z_L1',2*R('4.83')*(75+215)<3000)
    check('actual_h_square_L1',R('4.8')**2<24)
    check('remainder_gauge_chart_cost',3000*R('4e-8')==R('.00012'))
    check('full_affine_chart_cost',24*R('1.81e-5')==R('.0004344'))
    expmargin=R('.99')-(R('.01')+R('8e-24')+R('.00012'))
    check('physical_exponential_margin',expmargin>R('.9798'))
    check('physical_affine_margin',expmargin-R('.0004344')>R('.97'))
    check('uniform_action_parameter_margin',R('.99')>expmargin>R('.97'))


def exact_return_and_density():
    b,c,A,Rr,eta=S.symbols('b c A R eta')
    Sc=c/2-Rr*eta; Tc=-(b+3*c)/2-(A-Rr)*eta
    zero('orbit_return_cancellation',Sc+Tc+(b+2*c)/2+A*eta)
    hh=S.Matrix([b,-b-2*c,c,c]); ss=S.Matrix([A*eta+b+c,A*eta,Rr*eta,Rr*eta])
    kk=hh-2*ss
    disk=S.Matrix([2*(Sc+Tc),2*(Sc+Tc),2*Sc,2*Sc])
    for j in range(4): zero('exact_disk_tensor_component_'+str(j),kk[j]-disk[j])
    zero('TT_trace',sum(hh))
    zero('TT_point_norm',(hh.T*hh)[0]-2*(b+c)**2-4*c*c)
    st,tt=S.symbols('S T',real=True)
    Bmat=S.Matrix([[2,2],[2,2],[2,0],[2,0]])
    uv=S.Matrix([st,tt]); kval=Bmat*uv
    Eg=S.Matrix(S.symbols('Er Et Es1 Es2'))
    F=Bmat.T*Eg
    zero('exact_disk_gradient_pairing',(F.T*uv)[0]-(Eg.T*kval)[0])
    zero('disk_volume_direction',sum(kval)/2-(4*st+2*tt))
    zero('returned_tensor_norm',(kval.T*kval)[0]-8*((st+tt)**2+st**2))
    check('sharp_orbit_return_L2',(R('2.002')+2)/2+R('1.008')*R('3.002')<R('5.1'))
    check('returned_full_tensor_L2',8*(R('5.1')**2+12**2)<1400)
    check('returned_density_L2',2*(R('5.1')+12)<=R('34.2'))
    check('density_energy_dual_coefficient',R('34.2')**2*1400<1300**2)
    check('density_energy_dual_bound',1300*R('1.81e-5')<R('.025'))
    check('returned_moving_volume_quadratic_margin',R('.90')-R('.025')>R('.87'))
    # Fixed-density chain rule and non-symmetry are independent of a field equation.
    x,y=S.symbols('x y'); zz=S.Matrix([x,y]); ll=S.Matrix([4,2])
    aa=S.Matrix([1,3]); HH=S.diag(3,5)
    I=(aa.T*zz)[0]+(zz.T*HH*zz)[0]/2
    Jvol=S.exp((ll.T*zz)[0]); gradI=S.Matrix([S.diff(I,v) for v in zz])
    moving=gradI/Jvol
    DF=moving.jacobian(zz).subs({x:0,y:0})
    fixed=(Jvol*moving).jacobian(zz).subs({x:0,y:0})
    for i in range(2):
        for j in range(2):
            zero('fixed_density_chain_'+str(i)+str(j),fixed[i,j]-HH[i,j])
            zero('moving_density_contact_'+str(i)+str(j),DF[i,j]+aa[i]*ll[j]-HH[i,j])
    check('moving_gradient_not_selfadjoint',DF!=DF.T)
    zero('moving_gradient_adjoint',(DF.T-(HH-ll*aa.T)).norm()**2)
    check('omitted_density_not_Hessian',DF!=fixed)
    J,mu,a1,a2,d1,d2=S.symbols('J mu a1 a2 d1 d2')
    zero('fixed_output_dual_pairing',mu*(d1*J*a1+d2*J*a2)-(J*mu)*(d1*a1+d2*a2))
    check('fixed_density_disk_bound_not_spent',R('.90')>R('.87'))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S167 finite comparison physical TT contacts",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        source_normalization()
        source_values_and_physical_bound()
        gauge_and_chart_algebra()
        exact_return_and_density()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S167 finite comparison physical TT contacts",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "physical_exponential_TT_form_positive": good,
        "physical_affine_TT_form_positive": good,
        "actual_remainder_gradient_uniformly_bounded": good,
        "physical_tensor_chart_contacts_bounded": good,
        "returned_density_contact_bounded": good,
        "fixed_density_disk_TT_form_identified": good,
        "moving_volume_derivative_selfadjoint": False,
        "full_quantum_operator_domain_identified": False,
        "physical_mass_fast_domain_membership_proved": False,
        "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
        "signed_compatibility_evaluated": False, "nonlinear_family_preserved": False,
        "useful_Newton_certificate": False, "S167_complete": False, "goal_complete": False,
        "unresolved": ["mixed blocks",
                       "full logarithmic operator domains and physical-direction membership",
                       "candidate coupled range and adjoint", "projection direction and cold/higher bounds",
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
