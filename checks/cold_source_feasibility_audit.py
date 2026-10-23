"""Frozen cold-row calculation and feasibility budgets, not a metric solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
from math import factorial, comb
import json, traceback
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='cold-source-feasibility-audit';PREVIOUS='fixed-cap-seed-source-three'
OWN=['notes/'+STEM+'.md','checks/cold_source_feasibility_audit.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def exp_bounds(x,n=90):
    """Exact positive series and geometric tail; x is nonnegative."""
    assert x>=0 and x<n+2
    lo=sum(x**j/factorial(j) for j in range(n+1))
    return lo,lo+x**(n+1)/factorial(n+1)/(1-x/(n+2))

def convert(H,gamma,frame=Q(1)):
    B={(0,k):factorial(k)*100**k*H[k] for k in range(4)}
    for r in range(3):
        for k in range(3-r):
            B[r+1,k]=B[r,k+1]+r*sum(comb(k,a)*gamma(a)*B[r,k-a] for a in range(k+1))
    return [frame*B[j,0]/(factorial(j)*100**j) for j in range(4)]

def glue(k=Q(15000),row=Q(400)):
    D=[Q(360000)]+[Q(360000)*241**(j-1) for j in range(1,4)]
    H=[Q('4505.6')*(row/400)**2*(D[0]+2400)]
    for j in range(1,4):
        T=sum((j-b+1)*150**(j-b)*k**b*D[b] for b in range(1,j+1))+1200*(j+1)*150**j
        H.append(Q('4505.6')*(row/400)**2*T/100**j)
    return convert(H,lambda a:20*factorial(a)*300**a,Q(6,5))

def domain():
    H=[]
    for j in range(4):
        L=6*18**j*Q('6e-8')*(1+2/(Q('16.776')-j-3))
        H.append(Q('1e19')*(j+1)*Q('.00352')*Q(1000,3)*L)
    return convert(H,lambda a:0 if a==0 else 2*factorial(a)*25**a)

def rational():
    check('quadratic_mass_cancellation',Q('2.31')*exp_bounds(Q('.075'))[1]/8<Q('.32'))
    check('complex_cap_potential',exp_bounds(Q('.025'))[1]-1<Q('.026'))
    check('weighted_cap_inverse',1/(1-Q('.026')-Q('.5'))<Q('2.11'))
    check('free_remote_exponential',1/exp_bounds(Q('2.25'))[0]<Q('.106'))
    d=Q('2.25');free2=(2+Q(3,4)/d+Q(9,64)/d**2)*Q('.106')/2
    check('free_remote_row',free2<Q('.354')**2)
    moment2=(2*factorial(4)+Q(3,4)*factorial(3)+Q(9,64)*factorial(2))/2
    check('quadratic_free_moment',moment2<Q('5.3')**2)
    check('frozen_distance_weight',exp_bounds(2*Q('1.06')*Q(5,32))[1]<Q('1.4'))
    check('normal_endpoint_factor',exp_bounds(Q('.025')/2)[1]<Q('1.1'))
    check('decay_reserve',1/exp_bounds(Q(18,32))[0]<Q('.57'))
    row=33*Q('1.4')*Q('1.1')*Q('.57')*(Q('.354')+Q('2.11')*Q('.32')*Q('5.3')/18**2)
    check('uniform_commutator_weighted_row_12',row<12,str(row))
    # Independent explicit third tensor recurrence; retains gamma_0.
    testH=list(map(Q,[2,3,5,7]))
    got=convert(testH,lambda a:20*factorial(a)*300**a,Q(6,5))
    explicit=Q(6,5)*(testH[3]+Q(20,100)*testH[2]+Q(6000+800,60000)*testH[1])
    check('third_full_connection_recurrence',got[3]==explicit)
    D=domain();G=glue();N=glue(row=Q(12));total=[N[j]+D[j]+Q('2e8') for j in range(4)]
    for j in range(4):
        check('quadratic_row_savings_'+str(j),N[j]==G[j]*Q(9,10000))
        check('complete_cold_jet_'+str(j),total[j]<Q('5e17'),str(total[j]))
    check('maximum_is_third',max(total)==total[3])
    check('tight_total_guard',max(total)<Q('4.603e17'))
    check('old_high_order_ledger_is_weaker',max(G[j]+D[j]+Q('2e8') for j in range(4))>Q('1e20'))
    lll=convert([Q(12)*Q(3,2)**j for j in range(4)],lambda a:20*factorial(a)*300**a,Q(6,5))
    check('LLL_curvature_tensor_reserve',max(lll)<Q('1e4'))
    check('LLL_trace_reserve',Q('.00352')*Q('1e4')<36<Q('1e8'))
    factor=2*Q('8.5e-33')*400**2
    check('unchanged_physical_normalization',factor==Q('2.72e-27'))
    middle=factor*Q('6e15');cold=factor*Q('5e17')
    check('middle_source',middle==Q('1.632e-11'))
    check('cold_source',cold==Q('1.36e-9'))
    check('real_improvement_same_units',Q('1e20')/Q('5e17')==200)
    check('allocated_middle_operator',Q('61.2745098039')*middle<=Q('1e-9'))
    check('allocated_cold_operator',Q('1.4705882352')*cold<=Q('2e-9'))
    check('conditional_unit_inverse_remainder',Q('5e-9')-middle-cold==Q('3.62368e-9'))
    target=Q('2.5e-9')/factor
    def budget(k,row):return max(glue(k,row)[j]+D[j]+Q('2e8') for j in range(4))
    check('row_threshold_enclosure',budget(Q(15000),Q('18.34149'))<target<budget(Q(15000),Q('18.34150')))
    check('derivative_threshold_enclosure',budget(Q('1921.26475'),Q(400))<target<budget(Q('1921.26476'),Q(400)))
    rho=Q('1e-8');Y=Q('5e-9');Z0=Q(1,4);Z2=Q('5e7')
    check('neighborhood_preservation',Y+Z0*rho+Z2*rho**2/2==rho)
    check('contraction',Z0+Z2*rho==Q(3,4)<1)
    check('a_posteriori_residual',Q('1.25e-9')/(1-Q(3,4))==Y)
    check('interior_conditional_embedding',(Q('.001')-Q('1.3e-4'))/rho==87000)
    check('cold_admission',(Q('5e-7')-rho)==Q('4.9e-7'))
    check('outer_admission',(Q('1e-7')-rho)==Q('9e-8'))
    check('old_outer_fails_certificate',Q('1.01e-7')>Q('1e-7'))
    check('block_allocation',(Q('3e-9')+Q('7.5e-10'))/(1-Q('.25'))==Y)
    check('LLL_domain_target',Q(1,2)/rho==Q('5e7'))
    g=Q(1,2000**10);CI=Q('9.1591500000000000015e103');CO=Q('2.365283519093731664e186')
    cs=Q('1421.40132786776365');cb=Q('4.52219350981899516e14');mhi=Q('2.577600186401461e12')
    start=Q('3.123143830205975e28');bc=Q('1e57')*g/800
    low_target=bc/start;source_target=bc*g/(Q('.97')*2*mhi)
    check('old_low_certificate_rejected',low_target<Q('3.90858e-8')<Q('.97'))
    check('old_source_target',source_target<Q('2.38393e-25'))
    check('latest_interior_loss',(cs+cb*CI)/source_target>Q('1.7374e143'))
    check('global_loss',(cs+cb*(CI+CO))/source_target>Q('4.4868e225'))
    check('scalar_value_tiny_at_N2_one',Q('3e20')*(Q('4.71875e-37')+Q('3.125e-36'))<Q('1.1e-15'))
    check('generic_energy_embedding_rejected',Q('5e-9')/Q('1.44e-10')<35<Q('2.6e15'))
    FIXTURES.append(dict(fixture='exact_row_source_and_solution_budgets',row_coefficient=str(row),
        old_gluing=list(map(str,G)),actual_domain=list(map(str,D)),new_gluing=list(map(str,N)),
        complete_cold=list(map(str,total)),normalized_cold=str(cold),normalized_middle=str(middle),
        actual_domain_normalized=str(max(D)*factor),old_low_target_upper=str(low_target),
        old_high_source_target_upper=str(source_target),full_metric_inverse_computed=False))

def numerical():
    records=[]
    for digits in [70,110]:
        mp.mp.dps=digits
        # Independent elementary evaluation of the new row bound; prior
        # free Bessel moments are retained, not expensively recomputed.
        d=mp.mpf(18)/8
        free=mp.sqrt((2+mp.mpf(3)/(4*d)+mp.mpf(9)/(64*d*d))*mp.exp(-d)/2)
        row=33*mp.exp(2*mp.mpf('1.06')*5/32)*mp.exp(mp.mpf('.025')/2)*mp.exp(-mp.mpf(18)/32)*(free+mp.mpf('2.11')*mp.mpf('.32')*mp.mpf('5.3')/18**2)
        check('independent_row_'+str(digits),row<12)
        h=mp.mpf(1)/(8*2000**10)
        Wr=2**(mp.mpf(2)/5)*(1+1/h)**(mp.mpf(3)/5)
        M=(mp.mpf('5e-9')/Wr-mp.mpf('4.71875e-37'))/(mp.mpf('.5')+mp.mpf('3.125e-36'))
        check('physical_radius_threshold_'+str(digits),mp.mpf('3.40058804e-29')<M<mp.mpf('3.40058805e-29'))
        records.append(dict(digits=digits,independent_row=mp.nstr(row,65),radius_value_threshold=mp.nstr(M,65)))
    for key in ['independent_row','radius_value_threshold']:
        check('cross_precision_'+key,abs(mp.mpf(records[0][key])/mp.mpf(records[1][key])-1)<mp.mpf('1e-60'))
    FIXTURES.append(dict(fixture='independent_elementary_row_and_return',records=records,
        actual_charged_kernel_computed=False,actual_metric_correction_computed=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_middle_passed',old['status']=='passed' and old['normalized_middle_canonical_source_bound']=='1.632e-11')
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        rational();numerical()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,
        unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        cold_canonical_seed_C3_bound='5e17 q^2 a_s in j!100^j units' if error is None else None,
        normalized_cold_canonical_source_bound='1.36e-9 a_s' if error is None else None,
        full_physical_source_C3_proved=False,Einstein_frame_metric_response_bounded=False,
        actual_coupled_cold_response_bounded=False,full_transformed_boundary_domain_quantified=False,
        admission_norm_inverse_quantified=False,nonlinear_solution_constructed=False,
        semiclassical_validity_proved=False,physical_stability_proved=False,
        component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
