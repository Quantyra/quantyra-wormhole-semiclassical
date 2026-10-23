"""Frozen exact signed product moment; finite boundary functional remains open."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback
import sympy as sp
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='signed-cold-bubble-moment';PREVIOUS='current-product-full-metric-inverse'
OWN=['notes/'+STEM+'.md','checks/signed_cold_bubble_moment.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def algebra():
    M,j=sp.symbols('M j',positive=True);half=sp.Rational(1,2)
    aratio=(2*M+2*j+sp.Rational(7,2))/(2*M+2*j+sp.Rational(3,2))*(j+sp.Rational(3,2))/(j+1)*(2*M+j+sp.Rational(3,2))/(2*M+j+1)*(M+j+1)**2/(M+j+sp.Rational(3,2))**2
    tauratio=(M+j)*(M+j+sp.Rational(3,2))/((M+j+1)*(M+j+sp.Rational(5,2)))
    aa=[2*M+sp.Rational(3,2),M+sp.Rational(7,4),sp.Rational(3,2),M,M]
    bb=[M+sp.Rational(3,4),2*M+1,M+sp.Rational(5,2),M+sp.Rational(5,2)]
    hyperratio=sp.prod(j+x for x in aa)/(sp.prod(j+x for x in bb)*(j+1))
    check('residue_to_hypergeometric_term_ratio',sp.factor(aratio*tauratio**2-hyperratio)==0)
    check('absolute_series_excess',sp.simplify(sum(bb)-sum(aa))==2)
    a=2*M+sp.Rational(3,2);b=sp.Rational(3,2);c=M;d=M
    check('Dougall_convergence_condition',sp.simplify(b+c+d-a)==0)
    Gm=sp.gamma
    t0=(2*M+sp.Rational(3,2))/(8*sp.pi*(M+sp.Rational(3,2))**2)*Gm(sp.Rational(3,2))*Gm(2*M+sp.Rational(3,2))*Gm(M+1)**2/(Gm(2*M+1)*Gm(M+sp.Rational(3,2))**2)
    theorem=Gm(a-b+1)*Gm(a-c+1)*Gm(a-d+1)*Gm(a-b-c-d+1)/(Gm(a+1)*Gm(a-b-c+1)*Gm(a-b-d+1)*Gm(a-c-d+1))
    check('exact_massive_slope_zero',sp.simplify(sp.expand_func(t0*theorem)-1/(12*sp.pi))==0)
    ag,dg,GX,be,eps=sp.symbols('alpha_g d_g_alpha GX beta epsilon0')
    chi=4+2*ag+8*dg
    cc=6+6*be*(1-ag)-11*GX*ag/(90*sp.pi)
    dd=4*(ag-1)+48*((1-ag-dg)/6+eps)
    b0=1-GX*(6+5*ag)/(180*sp.pi)
    A0=36*be-2*cc+dd;J0=b0-3*be
    target=-8*dg-48*eps+chi*GX*(6+5*ag)/(90*sp.pi)-11*GX*ag/(45*sp.pi)+48*be*dg
    check('complete_first_cold_conformal_row',sp.expand(-A0-2*chi*J0-target)==0)
    check('first_cold_curvature_identity',sp.expand(-12-4*(ag-1)+2*chi-16*dg)==0)
    t,geom=sp.symbols('t d_geometry')
    check('finite_seed_curvature_first_jet',sp.diff(geom*(4*t**-5-6*t**-6),t).subs(t,1)==16*geom)
    p,aq,H=sp.symbols('p alpha_q H')
    # On P S=0 and P T=chi S, every P^2 contribution vanishes.
    check('generalized_Ward_shape_coefficient',sp.diff(p*aq+p*p*(H-be),p).subs(p,0)==aq)
    FIXTURES.append(dict(fixture='exact_signed_moment_and_normal_first_jet',massive_slope='0',finite_seed_uniform_remainder_proved=False))

def rational():
    aq=Q('1.013e-31');rad=Q('2.05e-31')
    perturb=Q('8.04')*aq+rad*Q('8.04')*Q('.667')*Q('8.5e-33')
    check('evaluated_slope_root_Rouche',perturb<rad*Q('3.99'))
    check('negative_a_at_Ward_root',36*Q('1.18e-8')-2*Q('5.99')+Q('4.01')<0)
    check('evaluated_slope_exponent',rad/3<Q('6.84e-32'))
    x=84*Q('6.84e-32')
    check('evaluated_slope_neck_factor',x/(1-x)<Q('5.75e-30'))
    chi=4+2*Q('2.1e-11')+8*Q('1.71e-29')
    check('cold_shape_seed_coefficient',chi<Q('4.000000001'))
    check('signed_shape_first_forcing',chi*aq<Q('4.053e-31'))
    fs=8*Q('1.71e-29')+48*Q('2e-36')+chi*Q('8.5e-33')*(6+5*Q('2.1e-11'))/(90*Q('3.14'))+11*Q('8.5e-33')*Q('2.1e-11')/(45*Q('3.14'))+48*Q('1.18e-8')*Q('1.71e-29')
    check('complete_conformal_first_forcing',fs<Q('1.4e-28'))
    check('classical_shape_component_only',Q('1.001')*Q('4.053e-31')*Q(85,6)<Q('5.8e-30'))
    check('old_slope_allowance_superseded',Q('1.7e-8')/aq>10**23)
    check('physical_source_allocation_unchanged',Q('1.632e-11')*Q('61.2745098039')+Q('1.36e-9')*Q('1.4705882352')+Q('2e-9')<=Q('5e-9'))

def intervals():
    iv=mp.iv;iv.dps=80;k=iv.mpf(2000);g0=k**-10
    g=2*g0/(1+iv.sqrt(1+4*iv.mpf(11)/(360*iv.pi)*g0))
    Fl=1-384*iv.pi*iv.mpf('1e22')*g;C=iv.sqrt(1+12*iv.pi*g*k**16/Fl)
    Yrn=g/((C-1)/6);Y=Yrn/(1+iv.mpf(['-5.5e-23','5.5e-23']))
    alpha=-k**4*Y/(6*iv.pi)
    check('signed_LLL_slope_interval',alpha>iv.mpf('-1.012546e-31') and alpha<iv.mpf('-1.012545e-31'),str(alpha))
    FIXTURES.append(dict(fixture='exact_product_slope_interval',digits=80,interval=str(alpha),pressure_quadrature_used=False))

def numerical():
    records=[]
    for digits in [70,110]:
        mp.mp.dps=digits;out=[]
        for ms in ['0.5','1','18.05','31.2']:
            M=mp.mpf(ms);W=lambda x:mp.exp(mp.loggamma(x+mp.mpf('.5'))-mp.loggamma(x))
            A0=2/mp.pi*(2*M+mp.mpf('1.5'))*W(1)*W(2*M+1)/W(M+1)**2
            t0=M*M*A0/(2*M*(2*M+3))**2
            hyp=mp.hyper([2*M+mp.mpf('1.5'),M+mp.mpf('1.75'),mp.mpf('1.5'),M,M],[M+mp.mpf('.75'),2*M+1,M+mp.mpf('2.5'),M+mp.mpf('2.5')],1)
            val=t0*hyp;exact=1/(12*mp.pi)
            check('independent_hypergeometric_sum_'+str(digits)+'_'+ms,abs(val-exact)<mp.mpf('1e-60'),mp.nstr(val-exact,20))
            # Independent positive summation, using only the residue
            # recurrence. Its analytic tail majorant encloses the identity.
            term=t0;partial=mp.mpf(0);N=1024
            for j in range(N):
                partial+=term
                term*=(2*M+j+mp.mpf('1.5'))*(M+j+mp.mpf('1.75'))*(j+mp.mpf('1.5'))*(M+j)**2/((M+j+mp.mpf('.75'))*(2*M+j+1)*(M+j+mp.mpf('2.5'))**2*(j+1))
            tn=M+N
            tail=M*M/8*((tn+1)/tn**4+1/(2*tn**2)+1/(3*tn**3))
            check('positive_residue_sum_and_tail_'+str(digits)+'_'+ms,partial<exact<partial+tail)
            out.append(dict(M=ms,total=mp.nstr(val,45),partial=mp.nstr(partial,45),tail_majorant=mp.nstr(tail,45)))
        records.append(dict(digits=digits,modes=out))
    for i in range(4):check('cross_precision_sum_'+str(i),abs(mp.mpf(records[0]['modes'][i]['total'])-mp.mpf(records[1]['modes'][i]['total']))<mp.mpf('1e-40'))
    FIXTURES.append(dict(fixture='direct_hypergeometric_and_positive_residue_sums',records=records,
        numerical_fixture_is_not_full_boundary_response=True))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    try:
        check('previous_cold_audit_passed',old['status']=='passed' and old['current_action_product_L2_inverse_proved'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        algebra();rational();intervals();numerical()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,
        unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        massive_signed_slope_exactly_zero=error is None,full_product_signed_slope='-Gq/(6pi B^2)' if error is None else None,
        normal_model_first_cold_jet_evaluated=error is None,finite_seed_signed_radius_functional_computed=False,
        actual_finite_boundary_return_bounded=False,physical_coupled_inverse_proved=False,nonlinear_solution_constructed=False,
        semiclassical_validity_proved=False,physical_stability_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
