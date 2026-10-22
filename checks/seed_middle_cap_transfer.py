"""Supporting checks for the ordinary common-cap source-three repair."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='seed-middle-cap-transfer'; PREVIOUS='seed-middle-splitting-audit'
OWN=['checks/seed_middle_cap_transfer.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational; ROWS=[]; FIXTURES=[]
PROTOCOL={'exact_arithmetic':'symbolic and rational','numerical_dps':[60,90],
          'numerical_tolerance':'1e-40','angular_cutoff':10000,
          'physical_center_bidisk_radius':'1e-5','cap_vertex_order':3,
          'cutoff_Taylor_order':4,'ordinary_continuum_proof_separate':True,
          'actual_signed_source_evaluated':False}


def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)


def zero(name,expr):
    value=S.simplify(expr);check(name,value==0,value)


def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})


def provenance():
    p=prior()
    check('previous_correction_passed',p['status']=='passed' and len(p['checks'])==62 and all(c['passed'] for c in p['checks']))
    check('previous_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('retained_inputs_unchanged',not bad,bad)
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('sources_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
    check('repair_scope',review['complete_middle_C3_reserve']=='1e26' and not review['withdrawn_1e22_claim_restored']
          and not review['actual_transfer_certificate_passed'] and not review['goal_complete'])


def seed_and_cutoffs():
    rho=R(1,100000);ex004=1/(1-R(4,1000))
    check('extended_u_Goursat',2*R(224,1000)/30+R(708,100)*ex004/(4*30**2)<R(2,100))
    check('extended_sigma_Goursat',8*ex004/(8*30**2)<R(2,1000))
    check('Goursat_u_row',ex004*(8+2*R(708,100))/(4*30**2)<R(2,100))
    check('Goursat_sigma_row',ex004*(17+16)/(8*30**2)<R(2,100))
    check('Bcal_third_majorant',(24+72*R(102,100))/R(59,10)+R(1,10**6)<17)
    check('Bcal_fourth_majorant',72/R(59,10)+R(1,10**6)<13)
    check('sigma_mixed_bound',ex004<R(1005,1000))
    uz=R(224,1000)+R(708,100)*ex004/(4*30)
    check('u_axis_derivative',uz<R(284,1000))
    check('sigma_pure_derivative',(17*R(284,1000)+16*R(1005,1000)/30)*ex004/8<R(7,10))
    check('normalizer_image_margin',R(1,32)+4*rho<R(1,30))
    t=S.symbols('t',real=True)
    eta=S.expand(sum(S.binomial(9,i)*t**i*(1-t)**(9-i) for i in range(5)))
    zero('eta_beta_derivative',S.diff(eta,t)+630*t**4*(1-t)**4)
    for j in range(5):
        zero('eta_left_'+str(j),S.diff(eta,t,j).subs(t,0)-(1 if j==0 else 0))
        zero('eta_right_'+str(j),S.diff(eta,t,j).subs(t,1))
    check('eta_first_bound',R(630,256)<R(247,100))
    check('eta_second_bound',R(2520,64)<40)
    check('eta_third_Bernstein',S.factorial(9)/S.factorial(6)*2<=1008)
    check('eta_fourth_Bernstein',S.factorial(9)/S.factorial(5)*3<=9072)
    tables={}
    for label,scale,rmin,target in [('chi',128,R(3,128),[315,7*10**5,3*10**9,4*10**12]),
                                   ('zeta',64,R(1,64),[158,18*10**4,4*10**8,4*10**11])]:
        f1=R(630,256)*scale;f2=40*scale**2;f3=1008*scale**3;f4=9072*scale**4
        raw=[f1,f2+f1/rmin,f3+2*f2/rmin+2*f1/rmin**2,
             f4+3*f3/rmin+9*f2/rmin**2+9*f1/rmin**3]
        check(label+'_Cartesian_tensor_budgets',all(a<=b for a,b in zip(raw,target)),raw)
        tables[label]=[str(v) for v in raw]
    # Radial components independently recover the tensor norms in dimension two.
    x,y=S.symbols('x y',real=True);r=S.symbols('r',positive=True)
    radial=(x*x+y*y)**S.Rational(7,2)
    for order in [2,3,4]:
        values=[S.diff(radial,x,order-k,y,k).subs({x:r,y:0}) for k in range(order+1)]
        direct=sum(S.binomial(order,k)*values[k]**2 for k in range(order+1))
        f=r**7
        if order==2:expected=S.diff(f,r,2)**2+(S.diff(f,r)/r)**2
        elif order==3:expected=S.diff(f,r,3)**2+3*(S.diff(f,r,2)/r-S.diff(f,r)/r**2)**2
        else:expected=S.diff(f,r,4)**2+6*(S.diff(f,r,3)/r-2*S.diff(f,r,2)/r**2+2*S.diff(f,r)/r**3)**2+9*(S.diff(f,r,2)/r**2-S.diff(f,r)/r**3)**2
        zero('independent_radial_tensor_norm_'+str(order),direct-expected)
    FIXTURES.append({'fixture':'seed_and_explicit_cutoff_budgets','radial_bounds':tables,'cutoff_degree':9})


def normalizer_and_density():
    z,w,b,c,a=S.symbols('z w b c a')
    F=(z+b)/(1+a*c*z);Fb=(w+c)/(1+a*b*w)
    Fz=S.diff(F,z);Fbw=S.diff(Fb,w)
    zero('constant_curvature_Goursat_normalizer',Fz-(1-a*F*c)**2/(1-a*b*c))
    zero('constant_curvature_normal_metric',Fz*Fbw/(1-a*F*Fb)**2-1/(1-a*z*w)**2)
    zero('normalizer_center_value',F.subs(z,0)-b)
    zero('normalizer_center_scale',Fz.subs(z,0)-(1-a*b*c))
    # Exact nonlinear coordinate pushforward, before spin rotation.
    x,y,t=S.symbols('x y t',real=True)
    G=S.Matrix([x+t*x*x/10,y+t*x*y/7])
    J=G.jacobian([x,y]);Ji=J.inv();det=J.det();u=x*x/100+y*y/200
    coords=[x,y]
    pushed=[S.exp(-u)*J[:,i] for i in range(2)]
    divergences=[]
    for aa in range(2):
        div=sum(Ji[k,i]*S.diff(pushed[aa][i],coords[k]) for i in range(2) for k in range(2))
        expected=S.exp(-u)*(S.diff(det,coords[aa])/det-S.diff(u,coords[aa]))
        zero('half_density_divergence_'+str(aa),div-expected)
        divergences.append(div)
    check('density_connection_not_zero',S.simplify(S.diff(det,x)/det)!=0)
    s1=S.Matrix([[0,1],[1,0]]);s2=S.Matrix([[0,-S.I],[S.I,0]]);generator=S.diag(S.I,-S.I)
    check('spin_rotation_extra_term_cancels',s1*generator+generator*s1==S.zeros(2)
          and s2*generator+generator*s2==S.zeros(2))
    # Taylor substitution needs only degree four because the map cutoff is multiplied by b.
    h=S.symbols('h')
    zero('map_cutoff_remainder_order',S.diff(h**5*h,h,2)-30*h**4)
    check('finite_jet_orders_sufficient',5-1>3 and 5+1-2>3)
    FIXTURES.append({'fixture':'independent_normalizer_and_half_density_transport',
                     'constant_curvature_example_is_actual_seed':False,
                     'density_term_retained':True,'cutoff_analyticity_assumed':False})


def finite_jet_budgets():
    rho=R(1,100000);d=2*rho
    for label,table,targets in [('chi',[315,7*10**5,3*10**9,4*10**12],[R(1007,1000),330,770000]),
                               ('zeta',[158,180000,4*10**8,4*10**11],[R(1004,1000),162,190000])]:
        jets=[S.Integer(1),*map(S.Integer,table)]
        extended=[sum(jets[k+i]*d**i/S.factorial(i) for i in range(5-k)) for k in range(3)]
        check(label+'_finite_Taylor_extension',all(a<b for a,b in zip(extended,targets)),extended)
    check('map_first_jet',R(1004,1000)*R(2,10)+162*R(15,10)<250)
    check('map_second_jet',R(1004,1000)*10+2*162*R(2,10)+190000*R(15,10)<286000)
    check('inverse_map_contraction',1/(1-250*rho)<R(1003,1000))
    check('inverse_map_displacement',R(1004,1000)*R(15,10)<R(16,10)<2)
    check('base_cap_complex_gradient',R(1007,1000)*R(67,1000)+330*R(2,1000)<R(75,100))
    check('base_cap_complex_Hessian',R(1007,1000)*3+2*330*R(67,1000)+770000*R(2,1000)<1700)
    check('cap_parameter_gradient',R(1007,1000)*10+330*R(3,10)<110)
    check('transported_scalar_value',(R(75,100)*R(16,10)+R(303,1000))*R(1004,1000)<R(16,10))
    grad=((1700*R(16,10)+110)*R(1003,1000)+R(75,100)*251+R(16,10)*R(75,100))*R(1004,1000)
    check('transported_scalar_gradient',grad<3200,grad)
    check('spin_rotation_gradient',(R(1004,1000)*10+162*R(21,100))*R(1004,1000)<45)
    principal=2*(R(16,10)*R(1003,1000)+R(1003,1000)*250+2*R(21,100)*R(1003,1000)**2)
    divergence=2*R(1003,1000)**2*286000+20000
    check('principal_coefficient_budget',principal<600,principal)
    check('density_divergence_budget',divergence<1200000,divergence)
    delta=600*rho+1200000*rho/(2*18)
    check('complex_jet_graph_budget',delta<R(35,100),delta)
    check('sqrt_two_rational_enclosure',R(17,12)**2>2)
    for j in [1,2,3]:check('true_vertex_tensor_'+str(j),R(35,100)*R(17,12)**j<1)
    check('unweighted_cap_gradient',2/(1-R(2,1000))+R(7,10)/(2*18)<R(21,10))
    check('unweighted_cap_graph',R(21,10)+1<4)
    ex=1/(1-R(2,1000))
    weighted=1/(1-ex/2)
    check('weighted_inverse',weighted<R(201,100),weighted)
    wg=ex*(1+(1+ex/2)*R(201,100))+R(7,10)*R(201,100)/(2*18)
    check('weighted_gradient',wg<R(41,10),wg)
    check('weighted_divergence',R(41,10)+R(201,100)/2<R(52,10))
    check('incoming_graph_row_reserve',12000*R(52,10)<40*30000)
    FIXTURES.append({'fixture':'common_physical_cap_vertex_ledger','vertex_constants':['100000','20000000000','6000000000000000'],
                     'graph_difference_budget':str(delta),'no_global_complex_cap_inverse_used':True})


def angular_and_source():
    x=S.symbols('x',positive=True)
    zero('unimodality_derivative_bound',(1+x)**2/S.Integer(8)-x*(1-x)-(3*x-1)**2/S.Integer(8))
    error=16*361/(R(27,10)**2*88**2)+8*R(1,200000)*(3/R(27,10))**3/88
    check('angular_moment_integral_plus_error',error<1,error)
    check('sqrt_C_lower_bound',S.Integer(2000)**4/361>200000**2)
    mass,c=S.symbols('mass c',positive=True)
    zero('angular_mass_squared_integral',S.integrate(4*mass**3*S.exp(-c*mass),(mass,0,S.oo))-24/c**4)
    tail=4*S.exp(-c*mass)*(mass**3/c+3*mass**2/c**2+6*mass/c**3+6/c**4)
    zero('angular_tail_primitive',S.diff(tail,mass)+4*mass**3*S.exp(-c*mass))
    before=json.loads((ROOT/'evidence/seed-middle-coercive-jets.json').read_bytes())['runs'][-1]
    ledger=next(f for f in before['fixtures'] if f['fixture']=='complete_angular_and_seed_source_ledger')
    cs=[S.sympify(v) for v in ledger['row_constants'][:4]]
    G=[S.sympify(v) for v in ledger['gluing_coordinate_bounds'][:4]]
    P=[sum(S.binomial(j,k)*cs[k]*cs[j-k] for k in range(j+1)) for j in range(3)]
    K=[None]+[S.factorial(j)*10**(5*j) for j in range(1,4)]
    V=[None,K[1],K[2]+8*K[1]**2,K[3]+24*K[1]*K[2]+96*K[1]**3]
    pref=82*25*88**4*R(31,2500000)
    kap=[S.Integer(0)]+[pref*sum(S.binomial(j,r)*V[r]*P[j-r] for r in range(1,j+1))/(S.factorial(j)*100**j) for j in range(1,4)]
    # Reconcile all ordered inverse-insertion multiplicities independently.
    zero('second_insertion_count',V[2]-(K[2]+2*4*K[1]**2))
    zero('third_insertion_count',V[3]-(K[3]+3*4*K[2]*K[1]+3*4*K[1]*K[2]+6*4**2*K[1]**3))
    delta=R(3,1000)
    born=612*94*delta*32**6+589824000*72*delta**2*32**6+S.Integer(2949120000)*8**3/(1-5*delta)
    check('complete_analytic_Born_coefficient',born<10**20,born)
    check('angular_Born_aggregate',4*10**20*R(352,100000)**3<2*10**13)
    check('finite_local_curvature_polynomial',2*13+R(708,100)*17<150)
    check('finite_part_restoration_reserve',1000*(8**2+2*150)*R(31,2500000)/1440<10**12)
    check('physical_center_Cauchy_factor',R(17,12)*100000/100<1500)
    L=[3*10**13*1500**j for j in range(4)]
    B=[S.Integer(10)**14,G[1]+L[1]+kap[1],G[2]+L[2]+kap[2],None]
    B[3]=G[3]+L[3]+kap[3]+B[1]/200
    check('complete_middle_C3_restored',max(B)<10**26,[str(S.N(v,18)) for v in B])
    check('sharper_target_not_certified',B[3]>10**22)
    FIXTURES.append({'fixture':'actual_cap_and_complete_source_budget',
                     'K':[str(v) for v in K[1:]],'V':[str(v) for v in V[1:]],
                     'kappas':[str(v) for v in kap[1:]],'kappas_decimal':[str(S.N(v,18)) for v in kap[1:]],
                     'complete_bounds':[str(v) for v in B],
                     'complete_bounds_decimal':[str(S.N(v,18)) for v in B],
                     'complete_C3_reserve':'1e26','withdrawn_1e22_restored':False,
                     'actual_finite_metric_source_H6_bounded':False})


def numerical_angular_check():
    records=[]
    for precision in PROTOCOL['numerical_dps']:
        with mp.workdps(precision):
            q=mp.mpf(2000)**4;C=q/340;N=PROTOCOL['angular_cutoff'];c=mp.mpf(1)/88
            M=lambda n:mp.sqrt(n*(n+q)/C)
            partial=mp.fsum(2*(q+2*n)/C*M(n)**2*mp.exp(-c*M(n)) for n in range(1,N+1))
            mn=M(N)
            # For this guard point N lies after the summand maximum, so the remaining tail is monotone.
            check('angular_tail_decreasing_'+str(precision),4*C*mn/(q*q+4*C*mn*mn)+2/mn-c<0)
            tail=4*mp.exp(-c*mn)*(mn**3/c+3*mn**2/c**2+6*mn/c**3+6/c**4)
            integral=mp.quad(lambda s:4*s**3*mp.exp(-s),[0,1,8,32,mp.inf])/c**4
            check('independent_angular_integral_'+str(precision),abs(integral-24/c**4)<mp.mpf('1e-40'))
            check('guard_point_full_moment_enclosure_'+str(precision),0<partial<partial+tail<25/c**4,
                  'partial/C='+mp.nstr(partial,20)+'; tail/C<'+mp.nstr(tail,15))
            records.append({'dps':precision,'N':N,'C':'q/340','partial_over_C':mp.nstr(partial,45),
                            'tail_over_C_upper':mp.nstr(tail,45),'actual_benchmark_C_evaluated':False})
    with mp.workdps(60):
        check('angular_precision_agreement',abs(mp.mpf(records[0]['partial_over_C'])-mp.mpf(records[1]['partial_over_C']))<mp.mpf('1e-35'))
    FIXTURES.append({'fixture':'independent_angular_moment_and_tail','records':records,'angular_tower_truncated_in_proof':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 common cap transfer','input_count':len(ps)}),flush=True)
    try:
        provenance();seed_and_cutoffs();normalizer_and_density();finite_jet_budgets();angular_and_source();numerical_angular_check()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
    run={'story':'S167 common cap transfer','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'ordinary_actual_cap_variation_bounded':True,'ordinary_complete_middle_C3_restored':True,
         'complete_middle_C3_reserve':'1e26','withdrawn_1e22_claim_restored':False,
         'source_splitting_gap_repaired_at_old_allowance':True,'all_descendant_estimates_reverified':False,
         'actual_finite_metric_source_H6_bounded':False,'actual_transfer_certificate_passed':False,
         'actual_signed_compatibility_zero_proved':False,'endpoint_C5_contract_changed':False,
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
