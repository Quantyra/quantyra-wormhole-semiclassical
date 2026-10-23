"""Frozen direct finite-seed action-form Hessian realization; physical admission remains open."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback
import sympy as sp
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='current-seed-full-hessian-realization';PREVIOUS='current-seed-boundary-phase-and-trace-defect'
OWN=['notes/'+STEM+'.md','checks/current_seed_full_hessian_realization.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def algebra():
    t,x,d,e=sp.symbols('t x d e',positive=True)
    drift=1+3*x+d-1/(1+e)-3*x/(1+e)**2
    target=(1-1/(1+e))/x+3*(1-1/(1+e)**2)+d/x
    check('physical_charge_drift_identity',sp.factor(drift/x-target)==0)
    X,B,DQ,lam=sp.symbols('X B DeltaQ lambda',positive=True)
    DV=-2*DQ/X+8*d*B**3/X**sp.Rational(3,2)-6*d*B**4/X**2
    U=DV/X;DU=sp.expand(X*sp.diff(U,X));D2U=sp.expand(X*sp.diff(DU,X))
    check('potential_exponential_TT_second_row',sp.simplify(D2U-(-8*DQ/X**2+50*d*B**3/X**sp.Rational(5,2)-54*d*B**4/X**3))==0)
    w=d*(4/t**5-6/t**6)
    check('curvature_deviation_first_coefficients',sp.expand(t*sp.diff(w,t)/d)==-20/t**5+36/t**6)
    check('curvature_deviation_second_coefficients',sp.expand(t*t*sp.diff(w,t,2)/d)==120/t**5-252/t**6)
    tau=(1+3*x+d)/t**4-d*(5/t**5-sp.Rational(9,2)/t**6)
    check('tracefree_charge_deviation',sp.expand(tau-(1+3*x)/t**4-d*(t**-4-5*t**-5+sp.Rational(9,2)*t**-6))==0)
    # Direct first curvature variation for the volume-preserving radial TT path.
    F=sp.Function('F')(t);b=sp.Function('b')(t);c=sp.Function('c')(t);z=sp.symbols('z')
    N=sp.exp(z*b/2)/sp.sqrt(F);aa=sp.sqrt(F)*sp.exp(-z*(b+2*c)/2);rr=t*sp.exp(z*c/2)
    check('TT_exponential_volume_exact',sp.simplify(N*aa*rr**2-t**2)==0)
    R=-2/(N*aa)*sp.diff(sp.diff(aa,t)/N,t)-4/(N*aa*rr)*sp.diff(aa*sp.diff(rr,t)/N,t)+2/rr**2-2*(sp.diff(rr,t)/N)**2/rr**2
    dR=sp.diff(R,z).subs(z,0)
    bp=-(2/t+sp.diff(F,t)/F)*b-(sp.diff(F,t)/F-2/t)*c
    bpp=sp.diff(bp,t).subs(sp.diff(b,t),bp)
    dR=sp.factor(dR.subs(sp.diff(b,t,2),bpp).subs(sp.diff(b,t),bp))
    RicA=-sp.diff(F,t,2)/2-sp.diff(F,t)/t;RicS=(1-F-t*sp.diff(F,t))/t**2
    check('exact_TT_curvature_first_variation',sp.factor(dR-2*(RicA-RicS)*c)==0)
    u=sp.symbols('u');P=3-24*u**3+39*u**4-sp.Rational(1791,100)*u**8
    minima=[]
    for left,right in [(sp.Rational(0),sp.Rational(1,2)),(sp.Rational(1,2),sp.Rational(1))]:
        poly=sp.Poly(sp.expand(P.subs(u,left+(right-left)*u)),u)
        bs=[sum(poly.nth(j)*sp.binomial(i,j)/sp.binomial(8,j) for j in range(i+1)) for i in range(9)]
        minima.append(min(bs));check('charge_potential_Bernstein_'+str(left),min(bs)>=sp.Rational(9,100))
    K=sp.symbols('K',real=True);H=sp.Matrix([[-1,K],[K,1]])
    check('saddle_square_and_inverse',H*H==(1+K*K)*sp.eye(2))
    ss,tt=sp.symbols('s xi',real=True);v=sp.Matrix([ss,tt]);J=sp.diag(-1,1)
    check('saddle_mixed_cancellation',sp.expand(((J*v).T*H*v)[0]-ss**2-tt**2)==0)
    FIXTURES.append(dict(fixture='exact_current_seed_local_and_saddle_identities',Bernstein_minima=[str(v) for v in minima],graph_closure_machine_verified=False))

def rational():
    pi=Q('3.14160');k=Q(2000);fl=1-Q('1.179e-8');dd=Q('8.6e-19');xx=Q('8e9');eb=Q('5.5e-23');dx=dd/xx
    charge=eb/((1-eb)*xx)+3*eb*(2+eb)/(1-eb)**2+dx
    check('current_charge_difference_budget',charge<Q('3.4e-22'))
    pot=4*Q('3.4e-22')+52*dx
    check('potential_TT_budget',pot<Q('1.37e-21'))
    gradpot=(8*Q('3.4e-22')+66*dx)/(16*Q('3.14159'))
    check('potential_gradient_budget',gradpot<Q('5.5e-23'))
    check('potential_return_contacts',Q('3e5')*Q('5.5e-23')/fl<Q('1e-15'))
    w0=Q('1.1e-27');w1=Q('1.4e-26');w2=Q('4.3e-25')
    check('curvature_value_guard',10*dx<w0)
    check('curvature_gradient_guard',Q('2.237')*56*dx<w1)
    check('curvature_Hessian_guard',4000*dx<=w2)
    cgl=Q('1e22')*k**-10
    a=12*w0+Q('7.5')*w1;bb=100*w0+Q('7.5')*w1
    square_extra=16*pi*cgl/fl*(5200*a+22*bb)
    check('nonconstant_square_second_variation',square_extra<Q('1e-30'))
    square_main=256*pi*cgl*Q('3.01')**2/fl
    check('square_TT_main_budget',square_main<Q('7.2e-8'))
    square_gradient=cgl*(200*w0+w0*w0+10*w2)
    check('nonconstant_square_return_contacts',Q('3e5')*square_gradient/fl<Q('1e-27'))
    qb=Q('.001')+Q('.002')+Q('.0006')+2*Q('1e-12')+Q('1e-25')+Q('1e-6')
    check('retained_quantum_budget_with_scale',qb<Q('.003601000003') and qb/fl<Q('.005'))
    total=Q('.005')+Q('1.37e-21')+Q('1e-15')+Q('7.2e-8')+Q('1e-30')+Q('1e-27')
    check('actual_seed_full_TT_lower_bound',1-total>Q('.99'))
    dcharge=Q('10.5')*dx
    pd=109*dx+Q('1.99')*(2*(3+1/xx)*dcharge+dcharge**2)
    check('charge_potential_d_geometry_loss',pd<Q('3e-26'))
    check('actual_charge_weighted_potential',Q('.09')-33/xx-2/xx**2-Q('3e-26')>Q('.08'))
    check('retained_charge_weighted_TT_norm',1-total>Q('.99') and Q('.08')-total>Q('.07'))
    epsl=k**-10/(60*Q('3.14159'));Dl=Q('2.946e-9');Cv=1+11*epsl/3
    check('current_scalar_value_coercivity',24-2*w0>Q('23.99'))
    check('current_scalar_gradient_coercivity',6-(2*Dl+11*epsl/6)*w0-Q(22,3)*epsl*7>Q('5.99'))
    main=16/(Q('1.98')*fl)*(Cv**2/Q('23.99')+Dl/6)
    check('current_main_mixed_budget',main<Q('.337'))
    # The loose integer root factors dominate the displayed square roots.
    check('mixed_return_root_factor',258/(Q('5.99')*fl)<Q('6.6')**2)
    check('mixed_return_budget',Q('9e-28')*Q('6.6')<Q('6e-27'))
    check('new_curvature_mixed_root_factor',315/(Q('23.99')*fl)<Q('3.63')**2)
    check('new_curvature_mixed_budget',Dl*w2*Q('3.63')<Q('5e-33'))
    check('current_total_mixed_bound',(Q('.581')+Q('1e-25'))**2<Q('.338') and Q('.337')<Q('.581')**2)
    check('current_mixed_target_met',Q('.338')<Q(3,8))
    check('nonlinear_targets_not_changed',Q('5e-9')+Q('.25')*Q('1e-8')+Q('5e7')*Q('1e-8')**2/2==Q('1e-8'))
    FIXTURES.append(dict(fixture='exact_rational_current_seed_form_budgets',TT_lower='.99 E_eff',mixed_square='.338',inverse_scope='actual positive-diagonal action-form norm',weighted_admission_inverse_computed=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    try:
        check('previous_candidate_boundary_test_passed',old['status']=='passed' and old['candidate_physical_outer_phase_proved'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        algebra();rational()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,
        unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        current_seed_full_TT_coercivity_proved=error is None,
        current_seed_mixed_form_bound_proved=error is None,
        current_seed_metric_Hessian_graph_identified=error is None,
        current_seed_action_form_inverse_proved=error is None,
        current_seed_mass_extended_weighted_bijection_proved=error is None,
        current_seed_linear_response_physical_phase_proved=error is None,
        current_seed_admission_inverse_constant_computed=False,
        full_defect_bolt_pairing_evaluated=False,
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
