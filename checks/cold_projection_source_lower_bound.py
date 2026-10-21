"""S167: actual cold projection-source lower bounds and necessary signed-coefficient budget.

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
STEM = "cold-projection-source-lower-bound"
OWN = ["checks/cold_projection_source_lower_bound.py", "notes/"+STEM+".md"]
PREDECESSORS = [("quantile-cutoff-localization-and-jets", 337)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "actual fixed projection source lower bounds; no metric lower bound or nonexistence claim",
    "source_normalization": "P=sqrt(gamma0) p_cut; y=sqrt(lambda)rho0; J=[1,2]",
    "derivative_lower_bound": "sup_J |partial_y^j P_S|>.006 kbar^j/sqrt(h), j=0,...,5",
    "source_norm_lower_bounds": "cold value >3e43; proper C5 >1e90; cold proper C5 >2e119",
    "signed_coefficient_evaluated": False,
    "nonexistence_claimed": False,
    "physical_oscillation_sampling": False,
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
    for p in ["quantile-cutoff-localization-and-jets", "on-shell-finite-range-response",
              "on-shell-adjoint-error-and-trace-gradient", "cold-reference-overlap-and-adjoint-budget",
              "endpoint-preserving-compatible-projection", "stationary-scalar-phase-moment",
              "stationary-centered-complete-residual"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)



def normalization_and_geometry():
    lam,gamma,t,G,h,N,alpha=S.symbols('lambda gamma t G h N alpha',positive=True)
    v=S.symbols('v',real=True)
    U=lam**R(1,4)*S.sqrt(gamma)*t*G**R(1,4)
    sigma=N/lam**R(1,4);chi=v/U
    P=S.sqrt(gamma)*t**-4*chi/(alpha*sigma)
    zero('physical_source_normalization',S.sqrt(h)*P-t**-5*(G/h**2)**R(-1,4)*v/(alpha*N))
    B,beta,F=S.symbols('B beta F',positive=True)
    gamma_phys=4*S.pi*beta*B**3;x=lam*B*B
    zero('physical_Liouville_volume',U.subs({gamma:gamma_phys,G:F/x})**2-4*S.pi*beta*(B*t)**2*S.sqrt(F))
    check('normalization_lower',(R('2.9')/R('.99'))**2<R('8.7'))
    check('thermal_h_upper',1/R(2000)**10<R('1e-30'))
    z2=4*R('6.02');terms=[z2**n/S.factorial(2*n) for n in range(21)]
    ratio=z2/(42*41)
    upper=sum(terms)+terms[-1]*ratio/(1-ratio)
    check('cold_annulus_upper',0<ratio<1 and upper<70,upper)
    check('cold_annulus_lower',1+R('5.98')/2+R('5.98')**2/24>5)
    hh=R('1e-30')
    check('cold_annulus_bootstrap',70*hh<R('.001'))
    check('cutoff_plateau_left',R('.2')<1)
    check('cutoff_plateau_right',1+70*hh<R('1.01'))
    check('amplitude_upper',R('5.98')*(25-1)*R('.3')**4>1)
    check('amplitude_lower',(1+70*hh)**20*R('6.02')*(70**2-1)<14**4)
    check('local_weight_conversion',(1+70*hh)**2<R('1.000001'))
    check('local_adjoint_error',R('1.000001')*R('.007')<R('.007001'))
    check('normalized_source_L2',R('.3')/R('2.9')*(R('1.000001')+R('.007001'))<R('.105'))


def window_and_duality():
    z=S.symbols('z',real=True);k=S.symbols('k',positive=True);phase=S.symbols('phase',real=True)
    Z=12012*z**6*(1-z)**6
    zero('window_unit_integral',S.integrate(Z,(z,0,1))-1)
    check('window_maximum',Z.subs(z,R(1,2))<3)
    zero('window_monotonicity_factor',S.diff(Z,z)-72072*z**5*(1-z)**5*(1-2*z))
    for j in range(6):
        check('window_endpoint_jets_'+str(j),all(S.diff(Z,z,j).subs(z,q)==0 for q in [0,1]))
    check('window_L2_enclosure',S.integrate(Z**2,(z,0,1))<3<R('1.733')**2)
    check('window_variation',2*Z.subs(z,R(1,2))<6)
    check('weighted_window_variation',R('.3')*6+3*R('.3')<3)
    lower=(R(1,28)-R(3,4)/(4*10**17)-R('.3e-9')-R('.3')*R('.007001')*R('1.733'))/R('4.007')
    check('actual_projection_source_detected_moment',lower>R('.007'),lower)
    coeff=sum(abs(c) for c in S.Poly(Z,z).all_coeffs())
    check('window_coefficient_sum',coeff==12012*2**6<10**6)
    for j in range(1,6):
        dsum=sum(abs(c) for c in S.Poly(S.diff(Z,z,j),z).all_coeffs())
        check('window_derivative_coefficient_'+str(j),dsum<10**6*12**j)
        phi=Z*S.cos(k*z+phase-j*S.pi/2)
        terms=sum(S.binomial(j,i)*S.diff(Z,z,i)*k**(j-i)*S.cos(k*z+phase-i*S.pi/2) for i in range(j+1))
        zero('phase_Leibniz_identity_'+str(j),S.expand_trig(S.diff(phi,z,j)-terms))
        rem=10**6*((1+R(12,4*10**17))**j-1)
        check('window_derivative_remainder_'+str(j),rem<R('2e-10'),rem)
    check('duality_lower_bound_reserve',R('.007')-R('.105')*R('2e-10')>R('.006'))
    # Independent polynomial field with actual trigonometric test;
    # its integration-by-parts equality is computed, not assumed.
    field=1+z+2*z*z+z**6
    for j in [1,3,5]:
        phi=Z*S.cos(3*z-j*S.pi/2)
        left=S.integrate(S.diff(field,z,j)*phi,(z,0,1))
        right=(-1)**j*S.integrate(field*S.diff(phi,z,j),(z,0,1))
        zero('independent_integration_by_parts_'+str(j),left-right)
    phi=Z*S.cos(3*z-S.pi/2)
    omitted=S.diff(phi,z)-3*Z*S.cos(3*z)
    check('window_derivative_omission_detected',omitted.subs(z,R(1,4)).simplify()!=0)


def cold_and_coefficient_budgets():
    check('cold_weight_bound',70*81==5670)
    check('cold_norm_prefactor',R('.006')/5670>R('1e-6'))
    check('source_value_lower_bound',R('.006')*2000**5>R('1.9e14'))
    check('cold_source_value_lower_bound',R('1e-6')*2000**15>R('3e43'))
    high=R('.006')*2000**5*(4*10**17)**5/(S.factorial(5)*100**5)
    check('unweighted_source_C5_lower_bound',high>R('1e90'),high)
    cold=R('1e-6')*2000**15*(4*10**17)**5/(S.factorial(5)*100**5)
    check('cold_source_C5_lower_bound',cold>R('2e119'),cold)
    check('current_coefficient_absolute_allowance',R('21.3')*R('2.2e-13')/R('3e-17')**2<R('5.3e21'))
    check('cold_unit_budget_allowance_ratio',R('5.3e21')*R('3e43')==R('1.59e65'))
    lam,gamma,sigma,moment,p,source=S.symbols('lambda gamma sigma moment p source',positive=True)
    compatible=source-p*moment/sigma
    coeff=moment/(lam**2*S.sqrt(gamma)); P=S.sqrt(gamma)*p
    zero('signed_coefficient_scaling',lam**-2*(source-compatible)-coeff/sigma*P)
    # The algebra is homogeneous and does not require a positive moment.
    signed=S.symbols('signed',real=True)
    zero('signed_scaling_negative_control',(lam**-2*(source-compatible)-coeff/sigma*P).subs(moment,signed))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 cold projection-source lower bound','input_count':len(ps)}),flush=True)
    try:
        predecessors();normalization_and_geometry();window_and_duality();cold_and_coefficient_budgets()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 cold projection-source lower bound','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'actual_cold_source_lower_bounds_proved':good,'necessary_signed_coefficient_budget_proved':good,
         'actual_signed_projection_coefficient_evaluated':False,'metric_lower_bound_proved':False,
         'nonexistence_proved':False,'cold_metric_control':False,'nonlinear_family_preserved':False,
         'projected_nonlinear_solve':False,'signed_physical_compatibility_evaluated':False,
         'S167_complete':False,'goal_complete':False}
    if not good:
        run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text()
        run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text()
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
