"""S167: cold reference overlap and conditional full-adjoint error budget.

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
STEM = "cold-reference-overlap-and-adjoint-budget"
OWN = ["checks/cold_reference_overlap_and_adjoint_budget.py", "notes/"+STEM+".md"]
PREDECESSORS = [("finite-comparison-adjoint-and-fixed-projection", 299),
                ("comparison-path-complete-equation-residual", 329),
                ("projected-regional-pointwise-bound", 243),
                ("finite-comparison-scalar-resolvent", 297)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "actual scalar reference overlap positive; full-adjoint transfer conditional",
    "dual_norm": "sqrt(lambda) integral t^-4 |field|^2 dmu0",
    "cold_interval": "h <= s <= 5e-15",
    "cold_phase_difference": ".13",
    "reference_overlap_lower": "1.5",
    "reference_norm_upper": "4",
    "sufficient_full_adjoint_errors": "e0 <= .1 and eell <= .1; NOT established",
    "conditional_overlap_ratio_lower": ".023",
    "conditional_functional_ratio_upper": "43",
    "fixed_projection_transversality": False,
    "fixture_scope": "thermal polynomials, exact integral budgets, Bessel/cosine and finite Hilbert fixtures",
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
    for p in ["finite-comparison-adjoint-and-fixed-projection", "stationary-completed-metric-response-norm",
              "comparison-path-complete-equation-residual", "projected-regional-pointwise-bound",
              "finite-comparison-scalar-resolvent", "endpoint-preserving-compatible-projection",
              "stationary-scalar-phase-moment"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)












def exp_bounds(x,n=300):
    term=S.Integer(1);lo=term
    for j in range(1,n+1):
        term=term*x/j;lo+=term
    assert x<S.Integer(n+2)
    return lo,lo+(term*x/(n+1))/(1-x/(n+2))


def thermal_geometry():
    t,h,x,d=S.symbols('t h x d',positive=True)
    f=x+2*x/t+(3*x+1)/t**2+d/t**4
    P=t*(t-1)**2*f
    pred=(t-1)*((t+1)*(1+3*x*(t*t+1))/t**2-d*(t-3)/t**4)
    zero('exact_thermal_derivative',S.diff(P,t)-pred)
    vmax=R('1.001');amin=6/vmax**2-R('.001')
    amax=((vmax+1)*(R('.001')+3*(vmax*vmax+1))+R('.002'))/2
    check('cold_lower_coefficient',amin/vmax>R('5.98'))
    check('cold_upper_coefficient',amax<R('6.02'))
    zero('outer_f_comparison',1+2/t+3/t**2-6/t**2-(t-1)*(t+3)/t**2)
    zero('outer_thermal_ratio_monotone',S.diff(t/(t-1)**2,t)+(t+1)/(t-1)**3)
    check('outer_thermal_margin',7*R('1e-30')**2*vmax/R('.001')**2<R('.01'))
    s,hh=S.symbols('s h',positive=True)
    lc=1+S.log(R('.01')/s);root=S.sqrt(s*s-hh*hh)
    zero('cold_phase_integration_by_parts',S.diff(lc*root,s)+root/s-s*lc/root)
    zero('coherent_budget_monotone',S.diff(s*(lc+1),s)-lc)
    H=S.log(t/(t-1))-1/t-1/(2*t*t)
    zero('weighted_outer_primitive',S.diff(H,t)+1/(t**3*(t-1)))
    # Three enlarged rational thermal geometries, evaluated independently.
    for i,(xx,dd,hh0) in enumerate([(1000,-1,R('1e-8')),(10000,0,R('1e-12')),(1000,1,R('1e-33'))]):
        ff=lambda tt:xx+2*xx/tt+(3*xx+1)/tt**2+dd/tt**4
        vv=1+hh0;delta=vv*hh0**2*ff(vv)
        for j,ss in enumerate([2*hh0,R('.0001'),R('.001')]):
            tt=1+ss;value=(ss**2*ff(tt)-delta/tt)/xx
            check('thermal_fixture_'+str(i)+'_'+str(j),R('5.98')*(ss*ss-hh0*hh0)<value<R('6.02')*(ss*ss-hh0*hh0))


def cold_and_error_budgets():
    sc=R('5e-15');hmax=R(1,2000**10);hmin=hmax/8
    check('benchmark_gap_order',0<hmin<hmax<R('1e-30')<sc<R('.001'))
    check('cold_w_completion',R('4e-6')+R('1.202e-8')<R('4.1e-6'))
    check('log_coherent_ratio_lower',exp_bounds(R(43))[1]<sc/hmax)
    check('log_cold_length_upper',exp_bounds(R(72))[0]>2*R('.001')/hmin)
    check('log_outer_integral_upper',exp_bounds(R('6.91'))[0]>1001)
    check('log_phase_budget_upper',exp_bounds(R('28.33'))[0]>2*10**12)
    check('sqrt_upper_cold',R('2.454')**2>R('6.02'))
    check('sqrt_lower_cold',R('2.445')**2<R('5.98'))
    check('sqrt_lower_outer',R('2.447')**2<R('5.99'))
    check('coherent_weight_lower',1/(1+sc)**4>R('.99999'))
    check('coherent_length_lower',R('.99999')*43/R('2.454')>R('17.5'))
    outer=(R('6.91')-1/R('1.001')-1/(2*R('1.001')**2))/R('2.447')
    check('outer_weighted_length',outer<R('2.22'),outer)
    check('full_weighted_length',72/R('2.445')+R('2.22')<R('31.7'))
    phase=5*10**17*R('1.001')*R('4.1e-6')*sc*R('30.33')/R('2.445')
    check('cold_phase_coherence',phase<R('.13'),phase)
    eta=R(72000,4*10**17)
    check('reference_Bessel_error',eta/(1-eta)+(eta/(1-eta))**2<R('2.1e-13'))
    cosine=R(1,4)/(4*10**17*R('1e-8')*(1-R('3e-8')))+R('2.1e-13')
    check('reference_cosine_error_away_bolt',cosine<R('1e-9'))
    check('volume_amplitude_error',R('9e-8')/(1-R('9e-8'))<R('1e-7'))
    check('oscillatory_variation_budget',4/R('1.99')+R('1.001')*R('.001')*R('31.7')<10)
    err=(R('3e-9')+R('1.1e-7'))*R('31.7')+R('3e-8')+R(10,4*10**17)
    check('quadratic_error_budget',err<R('1e-5'),err)
    check('reference_norm_upper',R('31.7')/2+R('1e-5')<16)
    overlap=((2-R('.13')**2/2)*R('17.5')-R('31.7'))/2-R('1e-5')
    check('actual_reference_overlap_margin',overlap>R('1.5'),overlap)
    q=R('.1');full=R('1.5')-(R('.4')+R('1.1')*q)**2-8*q-q*q-R('1e-5')
    check('conditional_full_overlap_margin',full>R('.4'),full)
    check('conditional_seed_overlap_ratio',R('.4')/R('4.1')**2>R('.023'))
    check('conditional_functional_ratio',R('4.1')**2/R('.4')<43)
    check('conditional_projected_inverse_factor',1+R('1.006')*43<R('44.3'))
    check('normalization_error_budget',2*R('.0247')**2<R('.0013'))
    zero('remaining_adjoint_budget',R('.0013')+R('.0987')-R('.1'))


def independent_fixtures():
    import mpmath as mp
    mp.mp.dps=70
    for i,x in enumerate([mp.mpf(1),mp.mpf(3),mp.mpf(20),mp.mpf(1000)]):
        v=mp.sqrt(mp.pi*x/2)*mp.besselj(0,x)
        check('Bessel_bound_fixture_'+str(i),abs(v)<=1)
        check('Bessel_Volterra_fixture_'+str(i),abs(v-mp.cos(x-mp.pi/4))<=1/(4*x))
    k=mp.mpf(17);delta=mp.mpf(1)/10000;length=mp.mpf(3)
    fun=lambda y:mp.cos(k*y-mp.pi/4)*mp.cos(k*(1+delta)*y-mp.pi/4)
    direct=mp.quad(fun,[length*j/24 for j in range(25)])
    exact=mp.sin(k*delta*length)/(2*k*delta)+(mp.sin(k*(2+delta)*length-mp.pi/2)+1)/(2*k*(2+delta))
    check('nonconstant_phase_integral',abs(direct-exact)<mp.mpf('1e-60'),mp.nstr(abs(direct-exact),8))
    coherent=length*(1-(k*delta*length)**2/2)/2-1/k
    check('independent_coherent_overlap',direct>coherent>0)
    f0=S.Matrix([3,0,0]);fl=S.Matrix([R(3,5),R(29,10),0])
    c0=f0+S.Matrix([0,R(3,50),R(2,25)])
    cl=fl+S.Matrix([R(3,50),0,-R(2,25)])
    theta=S.diag(1,R(1,2),0);I=S.eye(3)
    norm2=lambda v:(v.T*v)[0]
    check('fixture_full_error_hypotheses',norm2(c0-f0)==R('.01') and norm2(cl-fl)==R('.01'))
    check('fixture_reference_hypotheses',norm2(f0)<16 and norm2(fl)<16 and (fl.T*f0)[0]>R('1.5'))
    check('fixture_quantile_energy',(c0.T*(I-theta)*c0)[0]<=R('.01')*norm2(c0))
    check('fixture_full_cut_overlap',(cl.T*theta*c0)[0]>R('.4'))
    # Normalize with lambda=16, retaining the physical density conversion.
    alpha=(c0.T*theta*c0)[0]/norm2(c0);nu=S.sqrt(norm2(c0))/2
    sig0=(c0.T*theta*c0)[0]/(4*alpha*nu)
    sig=(cl.T*theta*c0)[0]/(4*alpha*nu)
    zero('physical_seed_normalization',sig0-nu)
    check('fixture_normalized_overlap',sig/nu>R('.023'))
    bad=fl-S.Matrix([2,0,0])
    check('missing_error_hypothesis_detected',norm2(bad-fl)>R('.01') and (bad.T*theta*c0)[0]<0)


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 cold reference overlap and adjoint budget','input_count':len(ps)}),flush=True)
    try:
        predecessors();thermal_geometry();cold_and_error_budgets();independent_fixtures()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 cold reference overlap and adjoint budget','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'cold_reference_coherence_proved':good,'actual_scalar_reference_overlap_positive':good,
         'conditional_full_adjoint_budget_proved':good,'full_adjoint_error_bounds_proved':False,
         'fixed_projection_transversality':False,'useful_full_weighted_inverse_constant':False,
         'signed_physical_compatibility_evaluated':False,'nonlinear_family_preserved':False,
         'S167_complete':False,'goal_complete':False}
    if not good:
        run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text()
        run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text()
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
