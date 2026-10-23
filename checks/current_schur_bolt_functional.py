"""Actual scalar point returns and Schur source thresholds; no global pair."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,traceback
import sympy as sp
import mpmath as mp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_uniform_bolt_certificate_obstruction as previous
local=previous.local;york=previous.york;prior=previous.prior;endpoint=previous.endpoint
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1]
STEM='current-schur-bolt-functional';PREVIOUS='current-uniform-bolt-certificate-obstruction'
OWN=['notes/'+STEM+'.md','checks/current_schur_bolt_functional.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];encoded=prior.encoded;Poly=local.Poly;zp=local.monomial
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)
prior.check=check;york.check=check

def algebra():
    z,h,eps,lam,rc=sp.symbols('z h epsilon lambda rc',positive=True)
    lam0=lam+sp.Rational(11,6)*eps*lam**2;Cv=1+sp.Rational(11,3)*eps*lam
    check('complete_constant_scalar_column',sp.expand((2*(-12*lam+rc)+48*lam0)/lam-(24*Cv+2*rc/lam))==0)
    zh=1/(1+h)
    check('actual_bolt_weight',(1-(1-h)*zh).simplify()==2*h/(1+h))
    q,p,qp,qpp,T1,T2,T3,j=sp.symbols('q p qp qpp T1 T2 T3 j',real=True)
    expr=T3*p**3+3*T2*p*qp/2+T1*qpp*p/2
    U=3*qp*T2/2+qpp*T1/2
    check('third_geodesic_derivative',sp.expand(expr-p*(T3*p*p+U))==0)
    th=sp.symbols('th',real=True)
    for k,value in [(1,1),(2,sp.Rational(3,4)),(3,sp.Rational(5,8))]:
        check('exact_directional_moment_'+str(k),sp.integrate(sp.cos(th)**(2*k),(th,0,2*sp.pi))==value*sp.pi)
    N,G,B,ga,s4,Z=sp.symbols('N G B gamma s4 Z',positive=True)
    moment=(2*sp.pi*N/ga)*sp.pi*Z*s4/(525*sp.pi**2)
    action=N*sp.pi*B**2/(G*ga)
    check('physical_sixth_moment_normalization',sp.simplify(moment/action-2*G*s4*Z/(525*sp.pi*B**2))==0)
    check('Schur_norm_upper',1+sp.Rational(338,1000)+sp.Rational(22,10000)<sp.Rational(1341,1000))

def execute(old):
    p,h,zh,env=york.actual_parameters();NUMERICS.extend(york.NUMERICS)
    Q=Poly({(power,0):env[sym] for power,sym in zip(york.powers,york.qs)})
    C,B2,lam=p['C'],p['B2'],p['g'];x=lam*B2;du=p['eps']*C/x
    ee=1/(60*arb.pi());Cv=1+arb(11)/3*ee*lam
    alpha=24*Cv+2*du*(4*zp(5)-6*zp(6));ab=alpha.value(zh)
    m=1-arb('1.1e-8');K=1/(6*arb('.63')*m*(1-arb('1e-26')))
    wb=(2*h/(1+h))**(arb(2)/5);ell=(1+h)/h;W=arb.pi()*B2/p['gamma']
    check('positive_actual_scalar_constant_column',ab>23 and ab<25)
    check('positive_full_scalar_barrier_constant',K>0 and K<arb('.265'))
    restoration=json.loads((ROOT/'evidence/current-canonical-heat-restoration-assembly.json').read_bytes())['runs'][-1]
    check('known_local_restoration_is_accepted',restoration['status']=='passed' and restoration['canonical_local_and_LLL_assembled'])
    diagonals={r['i']:endpoint.raw.decode(r['value']) for r in restoration['numerics'] if r['kind']=='assembled_known_TT_entry' and r['i']==r['j']}
    oldscalar=next(r for r in old['numerics'] if r['kind']=='exact_scalar_Schur_budget')
    qs=endpoint.raw.decode(oldscalar['scalar_bolt_upper']).upper()
    eta=arb('.0022');MS=arb('1.341');target=arb('1.667e-9');half=(target-qs)/2
    s4=4*C*C/arb(p['q'])
    NUMERICS.append(dict(kind='actual_Schur_point_parameters',h=encoded(h),zh=encoded(zh),C=encoded(C),B2=encoded(B2),
        lambda_value=encoded(lam),x=encoded(x),du=encoded(du),C_v=encoded(Cv),alpha_bolt=encoded(ab),
        barrier_K=encoded(K),bolt_weight=encoded(wb),ell_b=encoded(ell),physical_action_abs=encoded(W),
        full_angular_fourth_inverse_mass_sum_upper=encoded(s4),eta=encoded(eta),Schur_norm_upper=encoded(MS),
        scalar_bolt_upper=encoded(qs),post_scalar_half_budget=encoded(half)))
    for i,tr in enumerate(local.physical_trials(Q)[2:]):
        parts=local.trace_row_variation(p,Q,tr)
        g=-(C/x)*zp(4)*sum(parts.values(),Poly(0))
        check('mixed_trace_source_polynomial_'+str(i),all(k>=0 and l==0 for k,l in g.c))
        gb=g.value(zh);c=gb/ab;res=g-c*alpha
        check('correlated_exact_bolt_zero_enclosed_'+str(i),res.value(zh).contains(0))
        Cs=sum((k*v.abs_upper() for (k,l),v in res.c.items() if k>=1),arb(0))
        err=K*Cs*wb
        check('actual_scalar_point_return_accuracy_'+str(i),err<arb('1e-7'),encoded(err))
        xb=tr['xi'].value(zh);center=xb/2+c
        functional=center+arb(0,err.abs_upper());low=functional.abs_lower()
        check('effective_point_functional_is_not_canceled_'+str(i),low>arb('.99')*xb.abs_upper(),encoded(functional))
        T=tr['u'];T1=T.derivative();T2=T1.derivative();T3=T2.derivative()
        U=arb(3)/2*Q.derivative()*T2+Q.derivative().derivative()*T1/2
        density=Q*(arb(5)/8*Q**2*T3**2+arb(3)/2*Q*T3*U+U**2)
        Z=density.integral(zh)
        check('finite_actual_third_flow_energy_'+str(i),Z>0)
        positive_geo=2*s4*Z/(525*arb.pi()*B2)
        check('full_angular_smooth_trial_term_small_'+str(i),positive_geo<arb('1e-28'))
        known=diagonals[i];lupper=(-known.lower()+positive_geo.upper())/(1-eta)
        check('actual_TT_diagonal_upper_positive_'+str(i),lupper>0)
        P_lower=ell*low/(MS*(W*lupper).sqrt())
        threshold=target*MS*MS*lupper/(eta*(1+eta)*ell*low)
        halfthreshold=half*MS*MS*lupper/(eta*(1+eta)*ell*low)
        check('strict_reduced_source_threshold_'+str(i),threshold>0 and threshold<arb('1e-38'))
        check('conditional_half_budget_is_stricter_'+str(i),halfthreshold>0 and halfthreshold<threshold/2)
        NUMERICS.append(dict(kind='actual_scalar_mixed_point_return',i=i,name=tr['name'],
            g_coefficients=g.serialized(),g_at_bolt=encoded(gb),constant_response=encoded(c),
            response_remainder_coefficients=res.serialized(),weighted_source_upper=encoded(Cs),point_error_upper=encoded(err),
            xi_bolt=encoded(xb),effective_point_center=encoded(center),effective_point_coefficient=encoded(functional),
            effective_point_absolute_lower=encoded(low),actual_full_scalar_inverse_used=True,
            auxiliary_constant_outer_limits_cancel=True,scalar_Galerkin_inverse_used=False))
        NUMERICS.append(dict(kind='actual_Schur_trial_diagonal_and_source_threshold',i=i,name=tr['name'],
            third_flow_density_coefficients=density.serialized(),third_flow_integral=encoded(Z),third_flow_moments=density.moment_records(zh),
            full_angular_positive_form_upper=encoded(positive_geo),known_normalized_TT_diagonal=encoded(known),
            actual_l_diagonal_upper_over_action=encoded(lupper),Schur_adjoint_norm_lower=encoded(P_lower),
            necessary_reduced_source_entire_budget=encoded(threshold),necessary_reduced_source_half_budget=encoded(halfthreshold),
            reduced_source_projection_computed=False,adjoint_upper_norm_computed=False,full_physical_certificate_executed=False))
        print(json.dumps(dict(event='trial_complete',name=tr['name'],point_return=functional.str(12),source_threshold=threshold.str(12))),flush=True)

def main():
    ctx.prec=256;mp.mp.dps=55;review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|
        {'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};deps=endpoint.dependencies();error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',
        protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_all_trial_obstruction_passed',old['status']=='passed' and old['all_trial_pair_obstruction'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        algebra();execute(old)
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths};depafter=endpoint.dependencies()
    status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        unchanged_dependencies=deps==depafter,checks=ROWS,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),
        actual_scalar_mixed_point_returns_bounded=status=='passed',effective_bolt_functional_not_canceled_on_trials=status=='passed',
        actual_TT_diagonal_upper_bounds=status=='passed',necessary_reduced_source_thresholds_evaluated=status=='passed',
        signed_reduced_source_projection_computed=False,Schur_forward_adjoint_solve_executed=False,
        actual_physical_correction_computed=False,full_solution_feasibility_test_executed=False,
        metric_admission_proved=False,nonlinear_solution_constructed=False,physical_stability_proved=False,
        semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),
                    *ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
    if status!='passed':raise SystemExit(1)

if __name__=='__main__':main()
