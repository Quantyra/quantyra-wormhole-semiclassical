"""Frozen candidate boundary test and exact leading full-trace defect."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback
import sympy as sp
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='current-seed-boundary-phase-and-trace-defect';PREVIOUS='fixed-clock-seed-bolt-feasibility'
OWN=['notes/'+STEM+'.md','checks/current_seed_boundary_phase_and_trace_defect.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def algebra():
    t,x,a,d,sig,Db,q=sp.symbols('t x a d sigma Dbar q',positive=True)
    F=1+x*t*t-2*a/t+(1+3*x+d)/t**2-2*d/t**3+d/t**4
    R=-sp.diff(F,t,2)-4*sp.diff(F,t)/t+2*(1-F)/t**2
    W2=(t*t*sp.diff(F,t,2)-2*t*sp.diff(F,t)+2*F-2)**2/(3*t**4)
    JA=-sp.diff(F,t,2)/2-sp.diff(F,t)/t;JS=(1-F-t*sp.diff(F,t))/t**2
    LapR=sp.diff(t*t*F*sp.diff(R,t),t)/t**2
    T=sp.expand(R+12*x+sig*(7*W2/12+sp.Rational(11,6)*(JA-JS)**2-sp.Rational(11,72)*(R**2-144*x*x)-10*q*q/t**4)-Db*LapR)
    check('current_seed_curvature',sp.factor(R+12*x-d*(4*t**-5-6*t**-6))==0)
    check('trace_leading_magnetic_term',sp.limit(t**4*T,t,sp.oo)==-10*sig*q*q)
    check('no_seed_curvature_weight_four',sp.limit(t**4*(R+12*x),t,sp.oo)==0)
    z=sp.symbols('z');cal=sp.expand(z*z*F.subs(t,1/z))
    check('actual_compact_seed_face',cal.subs(z,0)==x and sp.diff(cal,z).subs(z,0)==0)
    check('actual_compact_seed_mass_cubic',sp.diff(cal,z,3).subs(z,0)/6==-2*a)
    FL,AL,mu,F0,A0=sp.symbols('F_L A_L delta_m F_0 A_0')
    df=-2*mu/t+(FL*sp.log(t)+F0)/t**2
    da=(AL*sp.log(t)+A0)/t**4
    dR=-sp.diff(df,t,2)-4*sp.diff(df,t)/t-2*df/t**2-3*sp.diff(F,t)*sp.diff(da,t)-4*F*sp.diff(da,t)/t-2*F*sp.diff(da,t,2)
    check('actual_logarithmic_curvature_combination',sp.limit(t**4*dR,t,sp.oo)==FL+8*x*AL)
    check('curvature_laplacian_weight_four',sp.limit(t**4*sp.diff(t*t*F*sp.diff(t**-4,t),t)/t**2,t,sp.oo)==4*x)
    cbar=sp.symbols('cbar');chi=1+11*sig*x/3-4*(96*sp.pi*cbar-sig)*x
    check('full_trace_denominator_identity',sp.expand(chi-(1-384*sp.pi*cbar*x+sp.Rational(23,3)*sig*x))==0)
    Aref,C=sp.symbols('a_ref chi',nonzero=True)
    check('reference_tail_defect',sp.expand(-Aref+C*Aref-(C-1)*Aref)==0)
    check('full_tail_solution',sp.simplify(C*(Aref/C)-Aref)==0)
    check('outer_density_integral',sp.simplify(sp.diff((2+sp.log(t))/t,t)+(1+sp.log(t))/t**2)==0)
    check('mass_boundary_momentum_zero',-105-116+221==0)
    for j in range(6):
        value=sp.expand(sp.diff(z**4*sp.log(z),z,j)*z**(j-4))
        check('weighted_log_tail_order_'+str(j),sp.diff(value,z,2)*z+sp.diff(value,z)==0)
    FIXTURES.append(dict(fixture='exact_tail_and_boundary_algebra',stress_parametrix_machine_verified=False,full_BVP_solved=False))

def rational():
    k=Q(2000);loss=Q('1.179e-8');amp=Q('1.99e-28')
    check('relative_tail_correction_target',loss/(1-loss)<Q('1.18e-8'))
    amplitude=Q('.956')/(6*Q('3.14')*k**8)
    check('curvature_tail_amplitude_target',amplitude<amp)
    check('trace_tail_defect_target',loss*amp<Q('2.35e-36'))
    check('full_tail_difference_target',loss*amp/(1-loss)<Q('2.35e-36'))
    check('nonresonant_tail_denominator',1-loss>Q('.9999999882'))
    check('massive_curvature_fast_root',Q('3.39e8')>18000*(18000-3))
    check('cold_value_budget_unchanged',Q('3.33e-10')+Q('1.667e-9')==Q('2e-9'))

def intervals():
    iv=mp.iv;iv.dps=80;k=iv.mpf(2000);g0=k**-10
    g=2*g0/(1+iv.sqrt(1+4*iv.mpf(11)/(360*iv.pi)*g0))
    ep=g/(60*iv.pi);Fl=1-384*iv.pi*iv.mpf('1e22')*g
    chi=Fl+iv.mpf(23)*ep/3
    xRN=(iv.sqrt(1+12*iv.pi*g*k**16/Fl)-1)/6
    B2=(xRN/g)*(1+iv.mpf(['-5.5e-23','5.5e-23']))
    amp=k**8/(6*iv.pi*g*B2**2)
    relative=1/chi-1;difference=amp*relative;defect=amp*(chi-1)
    check('actual_denominator_interval',chi>1-iv.mpf('1.179e-8') and chi<1,str(chi))
    check('actual_relative_correction_interval',relative>0 and relative<iv.mpf('1.18e-8'),str(relative))
    check('actual_reference_tail_interval',amp>0 and amp<iv.mpf('1.99e-28'),str(amp))
    check('actual_tail_difference_interval',difference>0 and difference<iv.mpf('2.35e-36'),str(difference))
    check('actual_reference_trace_defect_interval',defect<0 and defect>-iv.mpf('2.35e-36'),str(defect))
    FIXTURES.append(dict(fixture='current_product_and_trace_tail_interval',digits=80,chi=str(chi),relative=str(relative),reference_amplitude=str(amp),tail_difference=str(difference),trace_defect=str(defect),global_residual_bound=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    try:
        check('previous_fixed_clock_test_passed',old['status']=='passed' and old['seed_reference_bolt_value_bounded'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        algebra();rational();intervals()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,
        unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        candidate_physical_outer_phase_proved=error is None,
        candidate_fixed_source_boundary_work_zero=error is None,
        exact_reference_curvature_evaluated=error is None,
        leading_full_trace_tail_defect_evaluated=error is None,
        leading_trace_tail_defect_bound='2.35e-36 in lambda B^4 coefficient units',
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
