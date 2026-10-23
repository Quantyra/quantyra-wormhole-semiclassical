"""Budget and regularity diagnostics; does not solve the quantum BVP."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback
import mpmath as mp
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
STEM='cold-defect-simultaneous-feasibility-audit'
PREVIOUS='current-seed-cold-action-defect'
OWN=['notes/'+STEM+'.md','checks/cold_defect_simultaneous_feasibility_audit.py','evidence/'+STEM+'-review.json']
ROWS=[]; FIXTURES=[]

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def budgets():
    rho=Q('1e-8'); y=Q('5e-9'); z0=Q('.25'); z2=Q('5e7')
    check('self_map_allocation',y+z0*rho+z2*rho*rho/2==rho)
    check('strict_contraction_allocation',z0+z2*rho==Q('.75'))
    check('seed_C3_allocations',Q('1.632e-11')*Q('61.2745098039')<=Q('1e-9') and Q('1.36e-9')*Q('1.4705882352')<=Q('2e-9'))
    gain=Q('578034.6820809248')
    check('cold_linear_envelope_allocation',Q('1.73e-15')*gain<=Q('1e-9'))
    check('cold_total_allocation',Q('3.33e-10')+Q('1e-9')+Q('6.67e-10')==Q('2e-9'))
    check('scalar_only_trace',Q(4,3)*Q('1e18')*Q('1.1e-31')<Q('1.466666667e-13'))
    check('coupled_trace_threshold',Q(4,3)*Q('1e18')*(Q('1.1e-31')+Q('1.25014e-27'))==Q('1.667e-9'))
    check('admission_thresholds',Q('1.3e-4')+rho*87000==Q('.001') and Q('4.9e-7')+rho==Q('5e-7') and Q('1.4e-7')+rho==Q('1.5e-7') and Q('.025')*rho==Q('2.5e-10'))
    check('local_sector_optimization_negligible',Q('1.06e-18')/Q('1.72e-15')<Q('.00062'))
    oldlim=Q('2.38393e-25')
    check('historical_interior_certificate_failure',Q('4.141944868545865e118')/oldlim>Q('1.7374e143'))
    check('historical_global_certificate_failure',Q('1.069626977892751e201')/oldlim>Q('4.4868e225'))
    # e^x <= 1/(1-x), 0<x<1, supplies a rational upper enclosure.
    density=(Q('2e10')+Q('.028')+4*Q('2.725e-12'))/(1-Q('2.45e-8'))
    M=Q('2.01e10'); ed=Q('.025')
    check('fixed_density_second_derivative_reserve',density<M,float(density))
    check('nonlinear_gain_threshold',z2/(M*ed*ed)>Q('3.980099502487562') and z2/(M*ed*ed)<Q('3.980099502487563'),float(z2/(M*ed*ed)))
    check('nonlinear_product_threshold',z2/M>Q('.002487562189054726') and z2/M<Q('.002487562189054727'),float(z2/M))
    detail=[]
    for label,g,ed_upper,m_lower,m_upper in [('classical_scale_only',Q('469.12'),Q('.002303'),Q('1.70532e8'),Q('170532061')),('linear_ceiling_only',gain,Q('6.561e-5'),Q('138400'),Q('138400.0000001'))]:
        allowance=M*g*ed*ed; need_ed2=z2/(M*g); need_m=z2/(g*ed*ed)
        check(label+'_nonlinear_allowance_exceeds_budget',allowance>z2,dict(allowance=float(allowance),ratio=float(allowance/z2)))
        check(label+'_embedding_threshold_scale',need_ed2<ed_upper**2,dict(exact_threshold_squared=str(need_ed2),rounded_column_is_not_sufficient_constant=True))
        check(label+'_source_second_derivative_threshold',m_lower<=need_m<m_upper,float(need_m))
        detail.append(dict(label=label,gain_is_proved=False,allowance=str(allowance),ratio=str(allowance/z2),required_ED_squared=str(need_ed2),required_M2=str(need_m)))
    check('classical_scale_exact_product',M*Q('469.12')*ed*ed==Q('5.89332e9'))
    FIXTURES.append(dict(fixture='exact_simultaneous_budget_arithmetic',diagnostics=detail,actual_inverse_evaluated=False))

def regularity():
    mp.mp.dps=70; diag=[]
    for n in [100,10000,1000000]:
        p=1+n**4; ell=mp.log(2+n)
        scalar=mp.mpf(n)**5/p; shape=scalar/ell
        # Invert the coupled principal model [-p,n^2;n^2,p log(2+n)].
        det=-(mp.mpf(p)**2*ell+n**4)
        inv00=p*ell/det; inv01=-mp.mpf(n)**2/det; inv11=-mp.mpf(p)/det
        err=max(abs(-p*inv00+n*n*inv01-1),abs(-p*inv01+n*n*inv11),abs(n*n*inv01+p*ell*inv11-1))
        check('coupled_principal_model_inverse_'+str(n),err<mp.mpf('1e-55'),str(err))
        check('inverse_fifth_derivative_growth_'+str(n),scalar>mp.mpf(n)/2 and shape>mp.mpf(n)/(2*ell) and abs(inv00)*n**5>mp.mpf(n)/2 and abs(inv11)*n**5>mp.mpf(n)/(2*ell))
        diag.append(dict(n=n,scalar_C5=str(scalar),shape_C5=str(shape),coupled_scalar_C5=str(abs(inv00)*n**5),coupled_shape_C5=str(abs(inv11)*n**5)))
    # Exact chi=n^-4 witness: chi^(1/4)/chi^(1/2)=n.
    check('quarter_power_not_sqrt_envelope',all((Q(1,n)**4)**1==Q(1,n**4) and Q(1,n)/Q(1,n*n)==n for n in [10,100,1000]))
    FIXTURES.append(dict(fixture='principal_order_and_bolt_exponent_diagnostics',diagnostics=diag,actual_full_kernel_sampled=False,actual_C5_correction_bounded=False))

def identities():
    # Nontrivial affine test of the fixed-density second derivative formula.
    # S(x,y)=c+a*x+b*y+d*x*y, bfield=x+y; coefficient xy in e^(2x+2y)S.
    c,a,b,d=map(Q,['.0002','.007','.009','4.3'])
    x,y=sp.symbols('x y')
    full=sp.exp(2*x+2*y)*(c+a*x+b*y+d*x*y)
    coefficient=sp.diff(full,x,y).subs({x:0,y:0})
    independent=4*c+2*a+2*b+d
    check('density_contact_polynomial_coefficient',coefficient==independent and coefficient>d)
    # Exact nonlinear residual is distinct from the linear defect.
    f0,H,h,M=Q(3),Q(5),Q(2,7),Q(11)
    linear=f0+H*h; nonlinear=linear+M*h*h/2
    check('recentered_nonlinear_remainder_nonzero',nonlinear-linear==Q(22,49))
    for a,b in [(Q(2),Q(3)),(Q(1,1000),Q(10**7)),(Q(23,7),Q(31,11))]:
        M,L,e=Q(7),Q(11),Q(13)
        check('normalization_invariant_'+str(a)+'_'+str(b),(a*M/b)**2*(b*L/a**2)*(b*e)==M*M*L*e)
    FIXTURES.append(dict(fixture='density_recentring_and_normalization_identities',full_nonlinear_map_estimated=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}; error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    try:
        check('previous_complete_cold_test_passed',old['status']=='passed' and old['complete_cold_fixed_density_action_defect_bounded'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        budgets(); regularity(); identities()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),simultaneous_budget_test_executed=error is None,inverse_order_regularity_test_executed=error is None,value_envelope_to_raw_C5_rejected=error is None,quarter_power_nonlinear_target_distinguished=error is None,global_full_Hessian_defect_bounded=False,full_defect_bolt_pairing_evaluated=False,actual_physical_correction_computed=False,raw_C5_admission_certified=False,composite_nonlinear_response_bounded=False,nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json'); receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run); out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
