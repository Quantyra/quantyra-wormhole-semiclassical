"""Frozen fixed-cap source budget; no full metric inverse is computed."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
from math import factorial, comb
import json, traceback
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='fixed-cap-seed-source-three';PREVIOUS='auxiliary-physical-return-feasibility'
OWN=['notes/'+STEM+'.md','checks/fixed_cap_seed_source_three.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def rational():
    ee=Q('1.005');rr=Q(1,30);uz=Q('.284');sz=Q('.034')
    check('Goursat_uz',Q('.224')+ee*Q('7.08')/4*rr<uz)
    check('Goursat_sz',ee*rr<sz)
    uzz=rr*ee/4*(8*uz+2*Q('7.08')*sz)
    check('Goursat_uzz',uzz<Q('.03'))
    szzw=ee/8*(17*uz+16*sz)
    check('Goursat_szzw',szzw<Q('.68'))
    check('Goursat_szz',Q('.68')*rr<Q('.023'))
    szzzw=ee/8*(13*uz**2+17*Q('.03')+4*17*uz*sz+16*Q('.023')+32*sz**2)
    szzww=ee/8*(13*uz**2+17*Q('1.78')+4*17*uz*sz+16*ee+32*sz**2)
    check('Goursat_szzzw',szzzw<Q('.34'))
    check('Goursat_szzz',Q('.34')*rr<Q('.012'))
    check('Goursat_szzww',szzww<Q('6.1'))
    check('pure_fourth_Cauchy',480*Q('.012')<Q('6.1'))
    check('Cartesian_third_guard',8*3*Q('.68')<100)
    check('Cartesian_fourth_guard',16*4*Q('6.1')<1000)
    check('actual_normal_p_guard',Q(1,5)<(2*Q('.224'))**2)
    check('actual_small_d_guard',Q('8.6e-19')/(1+6*Q('8e9')-Q('8.6e-19'))<Q('1e-24'))
    # Coefficient sums of the normalized degree-six canonical polynomial.
    Cmin=1+6*Q('8e9')-Q('8.6e-19');din=Q('1e-24');zz=Q('1.02')
    b3=(2/Cmin+Q(2,3)+Q('1e-59'))
    b4=Q(1,2)+1/(2*Cmin)+din
    bp=2/Cmin*zz+3*b3*zz**2+4*b4*zz**3+10*din*zz**4+6*din*zz**5
    bppp=6*b3+24*b4*zz+120*din*zz**2+120*din*zz**3
    bpppp=24*b4+240*din*zz+360*din*zz**2
    check('actual_polynomial_first',bp<Q('7.08'))
    check('actual_polynomial_third',bppp<17)
    check('actual_polynomial_fourth',bpppp<13)
    check('complex_second_from_real',Q('2.001')+Q('.02')*17<8)
    eta=Q('.99');d=Q('.19')
    moments=[]
    for j,upper in [(1,Q('1.6')),(2,Q('5.3'))]:
        a=(2*Q(factorial(2*j))/eta**(2*j+1)+Q(3,4)*factorial(2*j-1)/eta**(2*j)+Q(9,64)*factorial(2*j-2)/eta**(2*j-1))/2
        check('weighted_Bessel_moment_'+str(j),a<upper**2)
        moments.append(str(a))
    check('remote_free_row',(2+Q(3,4)/d+Q(9,64)/d**2)/(2*eta)<9)
    check('matched_mass_lower',18*(1-Q('.002'))>Q('17.9'))
    check('weighted_decay',2-1/(1-Q('.002'))>eta)
    check('weighted_resolvent',1/(1-Q('.002'))/(1-Q('.002')-Q(1,2))<Q('2.1'))
    check('quadratic_global_h',Q('100')/128**2>Q('.0041'))
    check('remote_md',Q('17.9')/88>Q('.19'))
    row=3+Q('2.1')*(Q('.07')/Q('17.9')*Q('1.6')+100/Q('17.9')**2*Q('5.3'))
    check('sharp_row',row<7)
    check('endpoint_row_reserve',Q('1.01')*row<10)
    check('row_derivative_coefficient',10*Q('1.01')*(9+4*Q(75,64)+2*Q(75,64)**2)<180)
    for j in range(1,4):
        bb=600*(j+3);ll=512*(j+3)
        for k in range(j+2):
            val=Q(11,10)*(Q(ll,bb)+Q(21,10000)*sum(Q(factorial(k),factorial(k-b))*Q(128,bb)**b for b in range(1,k+1)))
            check('coercive_recurrence_'+str((j,k)),val<1)
    rho=Q(1,500)
    T=[Q('.002'),Q('.697'),Q('1445.21'),Q(6143635)]
    global_phi=sum(T[j]*rho**j/factorial(j) for j in range(4))
    center_phi=3*rho**2/2+100*rho**3/6
    check('translated_global_phi',global_phi<Q('.0145'))
    check('translated_center_phi',center_phi<Q('7e-6'))
    check('global_mass_potential',1/(1-global_phi-center_phi)-1<Q('.015'))
    check('outer_linear_h',Q('.015')*80<2)
    check('inner_gradient_change',3*rho+100*rho**2/2+1000*rho**3/6<Q('.0063'))
    inner=[]
    for k in range(2,6):
        val=Q('.002')*Q('.96')**k*sum(comb(k+j,j)*(96*Q('1.415')*rho)**j for j in range(1,4))
        check('inner_divided_jet_'+str(k),val<Q('.011'))
        inner.append(str(val))
    local=Q(5,3)*Q('1e17')*Q(31,2500000)
    check('local_Born_full_tower',local<Q('2.1e12'))
    cs=[Q(10)]+[Q(180*(j+3)*(600*(j+3))**j) for j in range(1,4)]
    pref=8*256**2*88**2*Q(31,2500000)
    Gs=[pref*sum(comb(j,k)*cs[k]*cs[j-k] for k in range(j+1))/(factorial(j)*100**j) for j in range(4)]
    Ls=[Q('2.1e12')]+[Q(2**((j+1)//2)*j**j)*Q('2.1e12')/(factorial(j)*100**j*rho**j) for j in range(1,4)]
    check('local_third_exact',Ls[3]==Q('4.725e15'))
    complete=[Q('1e14'),Gs[1]+Ls[1]+Q('1e12'),Gs[2]+Ls[2]+Q('1e12'),Gs[3]+Ls[3]+Q('1e12')+(Gs[1]+Ls[1]+Q('1e12'))/200]
    check('complete_fixed_cap_C3',max(complete)<Q('6e15'))
    norm=2*Q('8.5e-33')*400**2*Q('6e15')
    check('physical_canonical_source_normalization',norm==Q('1.632e-11'))
    check('old_source_term_diagnostic',Q('.06')*Q('6.8e-28')*Q('6e15')==Q('2.448e-13'))
    check('old_curvature_term_diagnostic',Q('6.6')*Q('6.8e-28')*Q('6e15')==Q('2.6928e-11'))
    FIXTURES.append(dict(fixture='rational_seed_resolvent_translation_and_source_budget',moments=moments,
        global_translation_phi=str(global_phi),inner_divided_jets=inner,rows=list(map(str,cs)),gluing_jets=list(map(str,Gs)),
        local_jets=list(map(str,Ls)),complete_jets=list(map(str,complete)),normalized_source=str(norm),
        full_budget_return_threshold=str(Q('5e-9')/norm),half_budget_return_threshold=str(Q('2.5e-9')/norm)))

def numerical():
    records=[]
    for digits in [40,60]:
        mp.mp.dps=digits
        for t in ['.01','.19','1','18','100']:
            tt=mp.mpf(t);base=mp.sqrt(mp.pi/(2*tt))*mp.exp(-tt)
            check('Bessel_point_'+str((digits,t)),mp.besselk(0,tt)<=base and mp.besselk(1,tt)<=base*(1+3/(8*tt)))
        # Independent scalar Bessel integral for the weighted free row
        # and its moments. These verify the envelopes, not the curved inverse.
        ratio=1/mp.mpf('.998')
        vals=[]
        for j in [1,2]:
            def f(t):return t**(2*j+1)*mp.exp(ratio*t)*(mp.besselk(0,t)**2+mp.besselk(1,t)**2)/mp.pi
            val=mp.quad(f,[0,1,5,20,mp.inf])
            check('Bessel_integrated_moment_'+str((digits,j)),val<mp.mpf('1.6' if j==1 else '5.3')**2)
            vals.append(mp.nstr(val,48))
        records.append(dict(digits=digits,weighted_moment_squares=vals))
    for j in range(2):
        check('cross_precision_moment_'+str(j),abs(mp.mpf(records[0]['weighted_moment_squares'][j])/mp.mpf(records[1]['weighted_moment_squares'][j])-1)<mp.mpf('1e-35'))
    FIXTURES.append(dict(fixture='independent_Bessel_integral_envelopes',records=records,actual_charged_kernel_computed=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_return_audit',old['status']=='passed' and old['N2_to_raw_C5_embedding_rejected'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        rational();numerical()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,
        unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        fixed_cap_middle_canonical_source_C3_bound='6e15 q^2 in j!100^j units' if error is None else None,
        normalized_middle_canonical_source_bound='1.632e-11' if error is None else None,
        withdrawn_moving_cap_combination_reinstated=False,full_physical_source_C3_proved=False,
        Einstein_frame_metric_response_bounded=False,actual_coupled_cold_response_bounded=False,
        full_transformed_boundary_domain_quantified=False,admission_norm_inverse_quantified=False,
        nonlinear_solution_constructed=False,semiclassical_validity_proved=False,physical_stability_proved=False,
        component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
