"""S167: fixed quantile cutoff localization and proper-radial jets.

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
STEM = "quantile-cutoff-localization-and-jets"
OWN = ["checks/quantile_cutoff_localization_and_jets.py", "notes/"+STEM+".md"]
PREDECESSORS = [("on-shell-finite-range-response", 340)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "original seed quantile cutoff only; full source and metric higher norms open",
    "location": ".06<y(a)<.2; .01h<a-t_h<.13h; 1.01<b<2; y(a_-)>.04; b_+<3",
    "jets": "|partial_y^j theta|<j! 10000^j, 1<=j<=5; seed dimensionless proper radius",
    "changed_projection": False,
    "physical_quantile_sampling": False,
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
    for p in ["on-shell-finite-range-response", "on-shell-adjoint-error-and-trace-gradient",
              "cold-reference-overlap-and-adjoint-budget", "endpoint-preserving-compatible-projection"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)



def cosh_bound(name,z2,bound):
    # cosh(sqrt(z2)), positive terms; after n=12 the ratio decreases.
    terms=[z2**n/S.factorial(2*n) for n in range(13)]
    ratio=z2/(26*25)
    upper=sum(terms)+terms[-1]*ratio/(1-ratio)
    check(name,0<ratio<1 and upper<bound,upper)


def quantiles_and_geometry():
    e=R('1e-5'); h=R('1e-30')
    check('reference_total_lower',R('2.958')**2<R('8.75')-e)
    check('full_total_lower',(R('2.958')-R('.007'))**2>R('8.7'))
    check('full_total_upper',(4+R('.007'))**2<R('16.2'))
    check('quantile_lower_mass',R('8.7')/200==R('.0435'))
    check('quantile_upper_mass',R('16.2')/200==R('.081'))
    check('initial_reference_upper',R('.06')/2+e<R('.174')**2)
    check('initial_full_upper',(R('.174')+R('.007'))**2<R('.0435'))
    check('initial_weight_lower',(1+R('1.13')*h)**4*R('.99999')<1)
    check('initial_reference_lower',R('.2')*R('.99999')/2-e>R('.316')**2)
    check('initial_full_lower',(R('.316')-R('.007'))**2>R('.081'))
    check('sqrt_upper',R('2.451')**2>R('6.002'))
    check('outer_interval_length_lower',R('2.3')/(R('1.1')**4*R('2.451'))>R('.6'))
    # log(10)>2.3 from exp(2.3)<10, with positive series tail.
    terms=[R('2.3')**n/S.factorial(n) for n in range(21)]
    ratio=R('2.3')/21
    check('log_ten_lower',sum(terms)+terms[-1]*ratio/(1-ratio)<10)
    check('tail_reference_lower',R('.3')-e>R('.54')**2)
    check('tail_full_lower',(R('.54')-R('.007'))**2>R('.081'))
    check('sqrt_lower',R('2.44')**2<R('5.99'))
    check('tail_reference_upper',1/(24*R('2.44'))+e<R('.0172')<R('.132')**2)
    check('tail_full_upper',(R('.132')+R('.007'))**2<R('.0435'))
    check('thermal_h_guard',1/R(2000)**10<h)
    cosh_bound('upper_first_quantile',R('.2')**2*R('6.02'),R('1.13'))
    cosh_bound('lower_transition_radius',R('.04')**2*R('6.02'),R('1.005'))
    check('first_quantile_thermal_width',R('5.98')*R('.06')**2/2>R('.01'))
    check('cold_bootstrap',R('1.13')*h<R('.001'))
    check('transition_order',1+R('1.13')*h<R('1.01'))
    check('outer_gap_ratio',R('.01')/(R('.01')-h)<R('1.000001'))
    t,ee,jj,dd=S.symbols('t e j delta',positive=True)
    original=(t-1)**2*(1+2/t+(3+ee)/t**2+jj/t**4)-dd/t
    expanded=t*t+ee-(4+2*ee+dd)/t+(3+ee+jj)/t**2-2*jj/t**3+jj/t**4
    zero('exact_seed_expansion',original-expanded)
    check('seed_first_derivative',R('12.004')+7*h<13)
    check('seed_fbar_derivative_negative',R('.004')<2)
    for n,bound in [(2,29),(3,100),(4,460)]:
        total=(2 if n==2 else 0)
        for power,coef in [(1,R('4.003')),(2,R('3.002')),(3,R('.002')),(4,R('.001'))]:
            total+=coef*S.rf(power,n)
        check('seed_derivative_'+str(n),total<bound,total)
    check('seed_inverse_t_coefficient',4+2*R('.001')+7*h*h<R('4.003'))


def coordinate_jets():
    t=S.symbols('t',positive=True); f=S.Function('f')(t)
    jets=[S.sqrt(f)]
    for _ in range(4):jets.append(S.diff(jets[-1],t)*S.sqrt(f))
    expected=[S.sqrt(f),S.diff(f,t)/2,S.diff(f,t,2)*S.sqrt(f)/2,
              (S.diff(f,t,3)*f+S.diff(f,t,2)*S.diff(f,t)/2)/2,
              (S.diff(f,t,4)*f**R(3,2)+3*S.diff(f,t,3)*S.sqrt(f)*S.diff(f,t)/2
               +S.diff(f,t,2)**2*S.sqrt(f)/2)/2]
    for n,(got,want) in enumerate(zip(jets,expected),1):zero('proper_chain_rule_'+str(n),got-want)
    h=R('1e-30')
    check('rise_first_jet',4*R('6.02')*R('2.01')/R('.01')<72**2)
    check('rise_second_jet',13*R('1.01')/R('.01')<1500)
    check('rise_t_first',R('6.02')*R('.13')*R('2.13')<4)
    check('rise_t_second',13*R('1.13')/2<8)
    check('rise_third_jet',R(29,2)*72<=1044)
    check('rise_fourth_jet',R(100,2)*(2*h)*72+R(29,2)*1500<22000)
    check('rise_fifth_jet',R(460,2)*(2*h)**2*72+R(300,2)*(2*h)*1500+R(29,2)*1044<16000)
    check('fall_first_jet',R('6.002')*(2*R('1.000001'))**2<25)
    check('fall_second_jet',13*R('1.000001')<14)
    check('fall_third_jet',R(29,2)*5<73)
    check('fall_fourth_jet',(100*25+29*14)/R(2)<1500)
    check('fall_fifth_jet',(460*125+300*5*14+29*73)/R(2)<41000)
    common=[100,1000,400,1300,500]
    for n,(a,b,c) in enumerate(zip([72,1500,1044,22000,16000],[5,14,73,1500,41000],common),1):
        check('common_divided_jet_'+str(n),max(a,b)<S.factorial(n)*c)
    check('rise_drift_denominator',R('5.98')*R('.005')*R('2.005')>R('.244')**2)
    check('rise_drift',13*R('1.13')/(2*R('.244'))+4*h<40)
    check('fall_drift',13*3/(2*R('2.44'))+2*R('2.451')<40)


def step_and_composition():
    # Independent direct differentiation checks the H recurrence.
    u,v,X=S.symbols('u v X',positive=True)
    P=S.Integer(1); A=[S.Integer(1)]; B=[S.Integer(1)]
    check('e_lower',sum(R(1,S.factorial(n)) for n in range(6))>R('2.7'))
    terms=[R(2)**n/S.factorial(n) for n in range(16)]
    ratio=R(2,16)
    check('step_denominator_lower',sum(terms)+terms[-1]*ratio/(1-ratio)<8)
    for n in range(1,6):
        P=S.expand(u*u*(P-S.diff(P,u)))
        direct=S.diff(S.exp(-1/v),v,n)/S.exp(-1/v)
        zero('independent_H_derivative_'+str(n),direct-P.subs(u,1/v))
        a=sum(abs(coef)*(R(m[0])/R('2.7'))**m[0] for m,coef in S.Poly(P,u).terms())/S.factorial(n)
        A.append(a)
        b=8*(a+sum(2*A[i]*B[n-i] for i in range(1,n+1)))
        B.append(b)
        check('step_divided_jet_'+str(n),b<50**n,b)
    Q=100*X+1000*X**2+400*X**3+1300*X**4+500*X**5
    envelope=S.Poly(sum(50**n*Q**n for n in range(1,6)),X)
    for n in range(1,6):
        coef=envelope.nth(n)
        check('cutoff_divided_jet_'+str(n),coef<10000**n,coef)
    multiplier=sum(100**i for i in range(6))
    check('weighted_multiplier_exact',multiplier==10101010101)
    check('weighted_multiplier_reserve',multiplier<R('1.011e10'))
    check('commutator_value_reserve',2*10000**2+40*10000<R('2.1e8'))
    check('commutator_gradient_coefficient',2*10000==20000)
    y=S.symbols('y',real=True); theta=y**3+2*y+1; field=S.exp(-y*y); drift=2+y*y
    def T(z):return -S.diff(z,y,2)-drift*S.diff(z,y)
    got=T(theta*field)-theta*T(field)
    want=-2*S.diff(theta,y)*S.diff(field,y)-(S.diff(theta,y,2)+drift*S.diff(theta,y))*field
    zero('independent_commutator',got-want)
    omitted=-(S.diff(theta,y,2)+drift*S.diff(theta,y))*field
    check('commutator_gradient_omission_detected',S.simplify((got-omitted).subs(y,1))!=0)


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths(); hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes(); error=None
    print(json.dumps({'event':'start','story':'S167 quantile cutoff localization and jets','input_count':len(ps)}),flush=True)
    try:
        predecessors();quantiles_and_geometry();coordinate_jets();step_and_composition()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 quantile cutoff localization and jets','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'original_quantile_cutoff_localized':good,'cutoff_proper_radial_C5_bounded':good,
         'cutoff_commutator_bounded':good,'full_projection_source_higher_norm_bounded':False,
         'cold_pointwise_control':False,'nonlinear_family_preserved':False,'projected_nonlinear_solve':False,
         'signed_physical_compatibility_evaluated':False,'S167_complete':False,'goal_complete':False}
    if not good:
        run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text()
        run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text()
    dest=ROOT/'evidence'/(STEM+'.json'); data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
