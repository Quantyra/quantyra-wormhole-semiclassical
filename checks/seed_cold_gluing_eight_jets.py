"""Supporting algebra, budgets and first-order locality controls for cold gluing."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from functools import lru_cache
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='seed-cold-gluing-eight-jets';PREVIOUS='seed-cold-local-eight-jets'
OWN=['checks/seed_cold_gluing_eight_jets.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational;ROWS=[];FIXTURES=[]
PROTOCOL={'exact_arithmetic':'symbolic and rational','numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','gluing_order':8,'jet_radius':'1e-4',
          'normalizer_radius':'5e-3','row_radius':'1e-2',
          'ordinary_continuum_proof_separate':True,'actual_charged_source_evaluated':False,
          'matching_moving_local_cap':True,'actual_domain_error_included':False}


def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)


def zero(name,expr):
    v=S.simplify(expr);check(name,v==0,v)


def previous():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(previous()['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})


def provenance():
    p=previous()
    check('prior_local_source_passed',p['status']=='passed' and len(p['checks'])==72 and all(r['passed'] for r in p['checks']))
    check('prior_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('all_prior_inputs_unchanged',not bad,bad)
    check('retained_sources_pinned',set(review()['retained_inputs'])<=set(paths()))
    check('auxiliary_scope_only',review()['cold_gluing_C8_bounded'] and review()['auxiliary_cold_source_C8_bounded']
          and not review()['complete_actual_cold_source_C8_bounded'] and not review()['actual_domain_C8_bounded']
          and not review()['goal_complete'])


def exp_upper(x):
    n=16
    return sum(x**j/S.factorial(j) for j in range(n))+x**n/S.factorial(n)/(1-x/(n+1))


def normalizer_and_coefficient_budgets():
    a=R(331,100000);small=R(1,10000);large=R(5,1000)
    excess=(1-small/large)**-2-1
    check('two_parameter_tail',excess<R(41233,10**6),excess)
    check('axis_height_patch',R(1,8)+R(9,8)*(exp_upper(R(1,100))-1)<R(137,1000))
    check('normalizer_exp_guard',exp_upper(R(6,100))<R(1062,1000))
    check('normalizer_weighted_displacement',R(1062,1000)*36*R(23,100)<9)
    check('normalizer_large_image',R(1062,1000)*R(23,100)<R(245,1000))
    check('normalizer_small_image',R(105,100)*R(7,32)+9*a*excess<R(232,1000))
    margin=R(18,1000)
    first=R(1062,1000)*(36+2*R(2,100)/margin*9)
    second=2*((12/margin+R(4,100)/margin**2*9)*R(1062,1000)**2
              +R(2,100)/margin*2*R(1062,1000)*60)
    check('first_normalizer_difference',first<60,first)
    check('second_normalizer_difference',second<4400,second)
    check('map_difference_first',60+30*9==330)
    check('map_difference_second',4400+60*60+1500*9==21500)
    check('map_uniform_first',R(439,1000)+330*a*excess<R(485,1000))
    check('map_uniform_second',22+21500*a*excess<25)
    check('normalizer_uniform_first',R(5,100)+60*a*excess<R(59,1000))
    check('normalizer_uniform_second',R(25,10)+4400*a*excess<R(311,100))
    check('inverse_map_norm',1/(1-R(485,1000))<2)
    check('weighted_normalizer_first',36*R(1062,1000)<40)
    check('weighted_normalizer_second',2*12/margin*R(1062,1000)**2<1600)
    check('weighted_map_first',40+30*9<=310)
    check('weighted_map_second',1600+60*40+1500*9==17500)
    check('weighted_inverse_map',2*310<=620)
    z=R(1,4);den=1-z*z/4
    height=(1+z*z/4)/den-1+z/den
    boost=(1+height)*(exp_upper(R(1,100))-1)
    check('full_spatial_height_patch',height+boost<R(30,100))
    check('half_disk_sD_max',(R(406,100)*R(2303,1000)+9)/1000<R(184,10000))
    Hprime=(R(36,1000)+2*R(184,10000)+R(184,10000)**2)/(1-R(184,10000))**2+R(1,10**24)
    check('height_inverse_derivative',Hprime<R(81,1000),Hprime)
    check('half_disk_inverse_image',(1-R(81,1000))/2>R(459,1000))
    check('scalar_complex_amplitude',R(5,1000)+R(184,10000)/(1-R(184,10000))+R(1,10**24)<R(25,1000))
    check('scalar_weighted_constant',5+9/(1-R(184,10000))<20)
    check('scalar_weighted_log',R(406,100)/(1-R(184,10000))<20)
    check('scalar_uniform_gradient',R(149,100)*(2*R(25,1000)/margin+R(13,100))<R(65,10))
    check('scalar_weighted_value',20+R(13,100)*9<22)
    grad=R(149,100)*40/margin+310*R(13,100)+R(149,100)*2*9
    check('scalar_weighted_gradient',grad<3500,grad)
    check('spin_uniform_generator',(R(311,100)/R(941,1000)+30*R(62,1000))/2<3)
    check('spin_weighted_generator',(1600/R(941,1000)+30*36)/2<1500)
    check('spin_norm',exp_upper(R(31,1000))<R(1032,1000))
    check('principal_scalar_frame_factor',exp_upper(R(7,100))*R(1032,1000)**2<R(115,100))
    check('principal_uniform',R(23,10)*(R(97,100)+R(8,100)+R(128,1000))<3)
    check('principal_weighted',R(23,10)*(620+44+80)<2000)
    div=R(115,100)*(4*4*25+R(1415,1000)*2*R(65,10)+2*3*2*R(1415,1000))+R(3,10)
    check('complex_Piola_uniform_divergence',div<520,div)
    check('complex_Piola_weighted_divergence',R(115,100)*(4*4*17500+11000+8500+1000)<400000)
    check('complex_gradient_factor_retained',4>2*R(1415,1000))
    FIXTURES.append({'fixture':'normalizer_and_complete_coefficient_budgets','parameter_tail':str(excess),
                     'uniform_divergence':str(div),'uniform_divergence_reserve':520,
                     'weighted_divergence_reserve':400000,'complex_Piola_factor':4,
                     'fixed_C2_coordinate_extension':True,'full_density_and_spin_included':True})


def independent_transport_algebra():
    x,y,e,tau=S.symbols('x y e tau',real=True);coords=[x,y]
    g=S.Matrix([x+e*x*x,y+e*x*y]);jac=g.jacobian(coords);det=jac.det();K=jac.inv()
    gammas=[S.Matrix([[0,S.I],[S.I,0]]),S.Matrix([[0,1],[-1,0]])]
    spin=gammas[0]*gammas[1]
    for i in range(2):
        zero('nonlinear_Piola_'+str(i),sum(S.diff(det*K[a,i],coords[a]) for a in range(2)))
        check('spin_anticommutator_'+str(i),spin*gammas[i]+gammas[i]*spin==S.zeros(2))
    principal=[sum((gammas[i]*K[a,i] for i in range(2)),S.zeros(2)) for a in range(2)]
    direct=-sum((principal[a]*S.diff(det,coords[a])/det/2 for a in range(2)),S.zeros(2))
    divergence=sum((principal[a].diff(coords[a])/2 for a in range(2)),S.zeros(2))
    check('density_transport_equals_half_divergence',(direct-divergence).applyfunc(S.simplify)==S.zeros(2))
    origin=direct.subs({x:0,y:0})
    check('omitted_density_control',origin==-R(3,2)*e*gammas[0] and origin!=S.zeros(2))
    Q=S.cos(tau*x/2)*S.eye(2)+S.sin(tau*x/2)*spin
    Qi=S.cos(tau*x/2)*S.eye(2)-S.sin(tau*x/2)*spin
    rotated=Q*gammas[0]*Qi
    direct_spin=(Q*gammas[0]*Qi.diff(x)).subs(x,0)
    div_spin=(rotated.diff(x)/2).subs(x,0)
    check('independent_spin_transport',direct_spin==div_spin==-tau*gammas[0]*spin/2)
    check('omitted_spin_control',div_spin!=S.zeros(2))
    FIXTURES.append({'fixture':'independent_Piola_density_spin_transport',
                     'map':['x+e*x^2','y+e*x*y'],'density_term_at_origin':'-3e gamma1/2',
                     'spin_term_at_origin':'-tau gamma1 gamma1gamma2/2',
                     'omitted_density_nonzero':True,'omitted_spin_nonzero':True})


def graph_rows_and_ordered_inverse():
    d=R(6,5);G=12;u=20;w=25000
    check('frozen_weight_gradient',2*R(1439,1000)*exp_upper(R(4,100))<3)
    check('physical_cap_weighted_inverse',18/(18-2*exp_upper(R(3,100)))<d)
    check('homogeneous_weighted_inverse',18/(18-3)<=d)
    physical=exp_upper(R(3,100))*(1+exp_upper(R(3,100))*d+10*d/18)
    check('physical_gradient_reserve',physical<4,physical)
    graph=R(1439,1000)*R(42,10)+(44+3)*d/18+3*d/18+d
    check('transported_graph_inverse',graph<G,graph)
    vertex=3*R(1021,1000)+3*R(2,10)/18+R(520,36)+R(9,18)
    weighted=2000*R(1021,1000)+2000*R(2,10)/18+R(400000,36)+R(6000,18)
    check('conjugated_vertex',vertex<u,vertex)
    check('double_weighted_vertex',weighted<w,weighted)
    check('full_tensor_Cauchy',15000**2>2*10000**2 and 150**2>2*100**2)
    check('normal_row_parameter',R(5155,10000)*2*R(1003,1000)*R(1021,1000)<R(11,10))
    check('cold_row_difference',R(33,10)/R(89,10000)+R(1,10**20)<1000)
    check('row_weight_conversion',exp_upper(2*R(6,100)*R(5,32))<R(102,100))
    @lru_cache(None)
    def compositions(n):
        if n==0:return ((),)
        return tuple((first,)+tail for first in range(1,n+1) for tail in compositions(n-first))
    entries=[]
    for j in range(1,9):
        terms=compositions(j)
        check('ordered_compositions_'+str(j),len(terms)==2**(j-1))
        enumerated=sum(d*G**len(c)*w*u**(len(c)-1) for c in terms)
        formula=d*G*w*(1+G*u)**(j-1)
        zero('inverse_composition_budget_'+str(j),enumerated-formula)
        entries.append({'order':j,'composition_count':len(terms),'D':str(formula)})
    FIXTURES.append({'fixture':'weighted_graph_rows_and_ordered_inverse','graph_reserve':G,
                     'uniform_vertex_reserve':u,'weighted_vertex_reserve':w,
                     'inverse_entries':entries,'one_cold_vertex_per_product':True,
                     'rows_varied_independently':True,'spatial_cutoff_high_derivatives_used':False})


def full_gluing_ledger():
    d=R(6,5);eps=R(1,100);k=15000
    D=[S.Integer(360000)]+[S.Integer(360000)*241**(j-1) for j in range(1,9)]
    T=[D[0]+2000*d]
    for j in range(1,9):
        by_counts=sum((j-b+1)*eps**(j-b)*D[b] for b in range(1,j+1))+1000*d*(j+1)*eps**j
        explicit=sum(eps**(a+c)*D[b] if b else 1000*d*eps**j
                     for a in range(j+1) for b in range(j-a+1) for c in [j-a-b])
        zero('independent_row_Leibniz_'+str(j),by_counts-explicit)
        T.append(by_counts)
    factor=2*400**2*4*R(352,100000)
    check('rank_mass_angular_factor',factor==R(22528,5))
    H=[factor*150**j*T[j] for j in range(9)]
    @lru_cache(None)
    def B(r,m):
        if r==0:return S.factorial(m)*100**m*H[m]
        return B(r-1,m+1)+(r-1)*sum(S.binomial(m,a)*20*S.factorial(a)*300**a*B(r-1,m-a) for a in range(m+1))
    entries=[]
    for j in range(9):
        J=R(6,5)*B(j,0)/(S.factorial(j)*100**j)
        check('physical_gluing_'+str(j),J<10**44,S.N(J,20))
        check('complete_auxiliary_'+str(j),J+10**8<10**44,S.N(J+10**8,20))
        entries.append({'order':j,'coordinate':str(H[j]),'physical_exact':str(J),'physical_decimal':str(S.N(J,22))})
    last=R(entries[-1]['physical_exact'])
    check('eighth_refined_bound',last<R(237,100)*10**43,S.N(last,22))
    FIXTURES.append({'fixture':'full_row_angular_covariant_ledger','entries':entries,
                     'Gamma_zero_order_included':True,'higher_tensors_symmetrized':False,
                     'full_angular_tower_retained':True,'complete_auxiliary_reserve':'1e44',
                     'actual_domain_comparison_included':False})


def independent_circle_locality():
    r0,ra=S.symbols('r0 ra')
    blocks=[r0*(1+ra)/2,r0*(1+ra)/2,-ra,-r0*r0]
    diagonal=(ra-r0)/((1-ra)*(1-r0))
    zero('circle_exact_double_commutator',sum(blocks)/((1-r0)**2*(1-ra))+diagonal)
    check('single_commutator_omission_nonzero',S.simplify(blocks[0])!=0)
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            mass=mp.mpf(2);alpha=mp.mpf('.7');x=mp.mpf(1)
            free=mp.exp(-4*mass);actual=free*mp.exp(-1j*alpha)
            def kernel(z,w,phase):return mp.exp(-mass*(z-w))*(1 if z>=w else phase)/(1-phase)
            intervals=[(mp.mpf('.5'),mp.mpf('.75'),4),(mp.mpf('1.25'),mp.mpf('1.5'),-4)]
            total=mp.mpc(0)
            for iz,(lo,hi,slope) in enumerate(intervals):
                for iw,(wl,wh,wslope) in enumerate(intervals):
                    def inner(z):
                        f=lambda w:kernel(x,z,free)*slope*kernel(z,w,actual)*wslope*kernel(w,x,free)
                        return mp.quad(f,[wl,z,wh] if iz==iw else [wl,wh])
                    total+=mp.quad(inner,[lo,hi])
            difference=1/(1-actual)-1/(1-free)
            error=abs(total+difference)
            check('independent_circle_quadrature_'+str(dps),error<mp.mpf('1e-42'),mp.nstr(error,12))
            check('circle_nonzero_comparison_'+str(dps),abs(difference)>mp.mpf('1e-5'))
            records.append({'dps':dps,'M':'2','alpha':'.7','error':mp.nstr(error,18),
                            'diagonal_difference_real':mp.nstr(difference.real,18),
                            'diagonal_difference_imag':mp.nstr(difference.imag,18)})
    FIXTURES.append({'fixture':'independent_first_order_circle_locality','records':records,
                     'first_order_local_operator':True,'actual_charged_source_evaluated':False,
                     'generic_graph_matrix_used':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 cold gluing eight jets','input_count':len(ps)}),flush=True)
    try:
        provenance();normalizer_and_coefficient_budgets();independent_transport_algebra()
        graph_rows_and_ordered_inverse();full_gluing_ledger();independent_circle_locality()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 cold gluing eight jets','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'cold_gluing_C8_bounded':True,'auxiliary_cold_source_C8_bounded':True,'auxiliary_cold_C8_reserve':'1e44',
         'full_physical_covariant_tensors_bounded':True,'full_angular_tower_retained':True,
         'density_and_spin_transport_retained':True,'matching_moving_local_cap':True,
         'complete_actual_cold_source_C8_bounded':False,'actual_domain_C8_bounded':False,
         'actual_metric_tenth_derivative_bounded':False,'actual_finite_metric_source_H6_bounded':False,
         'actual_transfer_certificate_passed':False,'original_cold_C3_contract_changed':False,
         'endpoint_C5_contract_changed':False,'actual_signed_compatibility_zero_proved':False,
         'nonlinear_solution_constructed':False,'physical_stability_proved':False,
         'physical_source_projection_changed':False,'physical_boundary_condition_changed':False,
         'ordinary_proof_machine_verified':False,'formal_verification':False,
         'S167_complete':False,'S168_complete':False,'goal_complete':False}
    if not good:
        for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');saved=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    saved['runs'].append(run);dest.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
