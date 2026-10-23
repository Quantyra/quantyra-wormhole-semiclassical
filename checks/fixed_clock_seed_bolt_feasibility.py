"""Frozen current-seed first reference response; full metric domain remains open."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback
import sympy as sp
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='fixed-clock-seed-bolt-feasibility';PREVIOUS='signed-cold-bubble-moment'
OWN=['notes/'+STEM+'.md','checks/fixed_clock_seed_bolt_feasibility.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def algebra():
    t,x,a,d=sp.symbols('t x a d',positive=True)
    F=1+x*t*t-2*a/t+(1+3*x+d)/t**2-2*d/t**3+d/t**4
    R=-sp.diff(F,t,2)-4*sp.diff(F,t)/t+2*(1-F)/t**2
    check('exact_seed_scalar_curvature',sp.factor(R-(-12*x+d*(4/t**5-6/t**6)))==0)
    primitive=t*sp.diff(R,t)-(R+12*x)
    f=-24/t**5+42/t**6
    check('R2_fixed_outer_clock_primitive',sp.factor(primitive-d*f)==0)
    check('R2_lapse_constraint_derivative',sp.factor(sp.diff(primitive,t)-t*sp.diff(R,t,2))==0)
    check('R2_outer_clock_zero',sp.limit(primitive,t,sp.oo)==0)
    check('lapse_extremum_derivative',sp.factor(sp.diff(f,t)-(120*t-252)/t**7)==0)
    critical=sp.Rational(21,10)
    check('lapse_maximum_at_one',f.subs(t,1)==18)
    check('lapse_negative_minimum',sp.factor(f.subs(t,critical)+4/critical**5)==0)
    check('lapse_variation_integral',sp.integrate(-sp.diff(f,t),(t,1,critical))+sp.integrate(sp.diff(f,t),(t,critical,sp.oo))==18+8/critical**5)
    # Eliminate the first constraint at the bolt and the fixed-temperature row.
    G,r,k,rho,aa,Fpp=sp.symbols('G R kappa rho a1 Fpp',positive=True)
    D=Fpp/2+k/r;sr=(G*rho/r-k*aa)/D
    fR=-2*k*sr;fpR=(-2*G*rho-fR)/r
    check('fixed_temperature_bolt_functional',sp.factor(fpR+Fpp*sr+2*k*aa)==0)
    FIXTURES.append(dict(fixture='exact_curvature_lapse_and_bolt_constraints',full_Hessian_solved=False))

def rational():
    k=Q(2000);beta=Q('1.18e-8');eps=Q('1.71e-29')
    wc_d=beta*eps/4
    check('R2_coefficient',wc_d<Q('5.045e-38'))
    check('R2_lapse_supremum',18*wc_d<Q('9.081e-37'))
    tv=18+8/Q('2.1')**5
    check('R2_total_variation_constant',tv<Q('18.2'))
    check('R2_total_variation',tv*wc_d<Q('9.182e-37'))
    check('finite_scale_Fang_guard',4+2*33+4*Q('1e-8')<71)
    check('finite_scale_a0_guard',(1000+(1+2*Q('1e-8'))/6)/(2*Q('3.14'))<167)
    scale=(3+1/Q('8e9'))/Q('3.14')/(1-Q('5.5e-23'))**2
    check('current_scale_ratio',scale<Q('.956'))
    massive=80*Q('1.001')/Q('.99')*Q('1e15')/k**8
    check('massive_lapse_bound',massive<Q('3.2e-10'))
    check('complete_Dirac_lapse',massive+Q('3e-28')+2/k**10<Q('3.3e-10'))
    hmax=k**-10
    # h/gamma differs by at most 12gamma; gamma<2h suffices.
    check('prepared_temperature_prefactor',1/((1-24*hmax)*(1-7*hmax))<Q('1.001'))
    radius=Q('1.001')*(Q('2e-12')+Q('3.3e-10')+Q('9.182e-37'))
    check('first_reference_bolt_value',radius<Q('3.33e-10'))
    check('remaining_cold_value_allowance',Q('2e-9')-radius>Q('1.667e-9'))
    check('whole_ball_regional_budget',Q('1.632e-11')*Q('61.2745098039')+Q('1.36e-9')*Q('1.4705882352')+Q('2e-9')<=Q('5e-9'))
    rr=Q('1e-8');y=Q('5e-9');z=Q('.25');zz=Q('5e7')
    check('ball_invariance_at_target',y+z*rr+zz*rr**2/2<=rr)
    check('strict_contraction_at_target',z+zz*rr<Q(1))
    check('interior_admission_threshold',Q('1.3e-4')+rr*87000<=Q('1e-3'))
    check('cold_admission_threshold',Q('4.9e-7')+rr<=Q('5e-7'))
    check('outer_admission_threshold',Q('1.4e-7')+rr<=Q('1.5e-7'))
    check('block_feedback_sufficient_target',(Q('3e-9')+Q('7.5e-10'))/(1-Q('.25'))<=Q('5e-9'))
    FIXTURES.append(dict(fixture='exact_rational_source_and_admission_substitution',
        first_reference_bolt_bound='3.33e-10',full_cold_value_target='2e-9',unproved_difference_target='1.667e-9',
        unknown_operator_constants_estimated=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    try:
        check('previous_signed_moment_passed',old['status']=='passed' and old['massive_signed_slope_exactly_zero'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        algebra();rational()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,
        unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        current_seed_quantum_lapse_transfer_checked=error is None,
        current_seed_R2_lapse_exactly_integrated=error is None,
        seed_reference_bolt_value_bounded=error is None,seed_reference_bolt_value_bound='3.33e-10',
        finite_seed_signed_radius_functional_computed=False,
        actual_finite_boundary_return_bounded=False,physical_coupled_inverse_proved=False,nonlinear_solution_constructed=False,
        semiclassical_validity_proved=False,physical_stability_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
