"""Supporting verification for complete middle seed jets through eight."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from itertools import product
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='seed-middle-eight-jets';PREVIOUS='repaired-source-metric-transfer'
OWN=['checks/seed_middle_eight_jets.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational;ROWS=[];FIXTURES=[]
PROTOCOL={'exact_arithmetic':'symbolic and rational','numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','source_order':8,'cutoff_degree':21,
          'cutoff_regularity':10,'physical_center_radius':'1e-5',
          'ordinary_continuum_proof_separate':True,'actual_signed_source_evaluated':False}


def check(name,ok,detail=''):
    assert name not in {row['name'] for row in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)


def zero(name,expr):
    value=S.simplify(expr);check(name,value==0,value)


def previous():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(previous()['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})


def provenance():
    p=previous()
    check('prior_transfer_review_passed',p['status']=='passed' and len(p['checks'])==53 and all(row['passed'] for row in p['checks']))
    check('prior_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('all_prior_inputs_unchanged',not bad,bad)
    check('retained_sources_pinned',set(review()['retained_inputs'])<=set(paths()))
    check('regional_scope',review()['complete_middle_source_order']==8 and review()['complete_middle_source_reserve']=='1e44'
          and not review()['actual_finite_metric_source_H6_bounded'] and not review()['goal_complete'])


def cutoff_branches():
    t=S.symbols('t');eta=S.expand(sum(S.binomial(21,i)*t**i*(1-t)**(21-i) for i in range(11)))
    c=S.factorial(21)/S.factorial(10)**2
    check('beta_normalization',c==3879876)
    zero('step_derivative',S.diff(eta,t)+c*t**10*(1-t)**10)
    for j in range(11):
        zero('step_left_'+str(j),S.diff(eta,t,j).subs(t,0)-(1 if j==0 else 0))
        zero('step_right_'+str(j),S.diff(eta,t,j).subs(t,1))
    rho=R(1,100000);d=2*rho;rmin=R(1,64);lower=R(9987,10000)
    qr=2*d/rmin+(d/rmin)**2
    check('complex_radial_square_root',lower**2<1-qr)
    check('complex_radial_displacement',(2+d/rmin)/(1+lower)<R(1002,1000))
    check('complex_radial_gradient',(1+d/rmin)/lower<R(1003,1000))
    check('complex_radial_Hessian',R(17,12)/lower+(1+d/rmin)**2/lower**3<R(5,2))
    delta=R(3,1000);B=R(1,4)+2*delta**2
    check('complex_step_argument',128*R(1002,1000)*d<delta)
    check('complex_step_first',c*B**10<R(38,10))
    check('second_step_root_enclosure',B/19<R(1,8)**2)
    second=20*c*(18*B/19)**9/R(8)+40*c*delta*B**9
    check('complex_step_second',second<40,second)
    budgets=[]
    for name,scale,rmin,targets in [('chi',128,R(3,128),[R(101,100),490,720000]),
                                  ('zeta',64,R(1,64),[R(1005,1000),245,210000])]:
        values=[1+R(38,10)*scale*R(1002,1000)*d,
                R(38,10)*scale*R(1003,1000),
                40*scale**2*R(1003,1000)**2+R(38,10)*scale*R(5,2)/rmin]
        check(name+'_analytic_branch_budgets',all(a<b for a,b in zip(values,targets)),values)
        budgets.append({'cutoff':name,'bounds':[str(v) for v in values]})
    check('common_operator_eight_derivatives',10-2>=8)
    FIXTURES.append({'fixture':'C10_cutoff_and_local_branches','budgets':budgets,
                     'global_cutoff_analyticity_assumed':False,'endpoint_jets_verified_through':10})


def vertex_budgets():
    rho=R(1,100000);fac=R(1004,1000)
    check('map_first',R(1005,1000)*R(2,10)+245*R(15,10)<370)
    check('map_second',R(1005,1000)*10+2*245*R(2,10)+210000*R(15,10)<316000)
    check('map_inverse',1/(1-370*rho)<fac)
    check('map_inverse_displacement',R(1005,1000)*R(15,10)<R(16,10)<2)
    check('base_gradient',R(101,100)*R(67,1000)+490*R(2,1000)<R(11,10))
    check('base_Hessian',R(101,100)*3+2*490*R(67,1000)+720000*R(2,1000)<1600)
    check('parameter_scalar_gradient',R(101,100)*10+490*R(3,10)<160)
    check('transported_scalar_value',(R(11,10)*R(16,10)+R(303,1000))*fac<R(21,10))
    gradient=((1600*R(16,10)+160)*fac+R(11,10)*372+R(21,10)*R(11,10))*fac
    check('transported_scalar_gradient',gradient<3200,gradient)
    check('spin_gradient',(R(1005,1000)*10+245*R(21,100))*fac<65)
    principal=2*(R(21,10)*fac**2+370*fac**3+2*R(21,100)*fac**3)
    divergence=2*fac**5*316000+20000
    check('full_principal_budget',principal<800,principal)
    check('full_divergence_budget',divergence<1200000,divergence)
    A=R(7,20);k=S.Integer(150000)
    check('branch_graph_difference',800*rho+1200000*rho/(2*18)<A)
    check('tensor_Cauchy_radius',R(17,12)*100000<k and R(17,12)**2>2)
    ex=1/(1-R(2,1000))
    check('base_unweighted_gradient',2*ex+R(11,10)/(2*18)<R(21,10))
    check('base_graph_inverse',R(21,10)+1<4)
    check('base_weighted_inverse',1/(1-ex/2)<R(201,100))
    weighted_grad=ex*(1+(1+ex/2)*R(201,100))+R(11,10)*R(201,100)/(2*18)
    check('base_weighted_gradient',weighted_grad<R(41,10),weighted_grad)
    check('weighted_incoming_divergence',R(41,10)+R(201,100)/2<R(52,10))
    check('incoming_row_reserve',12000*R(52,10)<40*30000)
    FIXTURES.append({'fixture':'true_eight_order_cap_vertex_budget','A':str(A),'k':str(k),
                     'half_density_divergence_retained':True,'global_complex_inverse_assumed':False})


def compositions(n):
    if n==0:yield ();return
    for a in range(1,n+1):
        for rest in compositions(n-a):yield (a,)+rest


def ordered_inverses():
    b=S.symbols('b');mat=S.Matrix([[2+b,b*b],[b+b**3,3-b+2*b*b]])
    inv=mat.inv();S0=mat.subs(b,0).inv()
    vertices={i:mat.diff(b,i).subs(b,0) for i in range(1,9)}
    check('independent_vertices_noncommute',vertices[1]*vertices[2]!=vertices[2]*vertices[1])
    series=inv.applyfunc(lambda v:S.series(v,b,0,9).removeO().expand())
    A=R(7,20);k=S.Integer(150000);counts=[]
    for r in range(1,9):
        terms=list(compositions(r));out=S.zeros(2);bound=S.Integer(0)
        for parts in terms:
            coeff=S.factorial(r)/S.prod(S.factorial(i) for i in parts)
            term=S0
            for i in parts:term=term*vertices[i]*S0
            out+=(-1)**len(parts)*coeff*term
            bound+=coeff*4**(len(parts)-1)*S.prod(A*S.factorial(i)*k**i for i in parts)
        direct=series.applyfunc(lambda v:S.factorial(r)*v.coeff(b,r))
        check('noncommuting_inverse_order_'+str(r),out==direct)
        zero('ordered_insertion_majorant_'+str(r),bound-A*S.factorial(r)*k**r*(1+4*A)**(r-1))
        check('ordered_composition_count_'+str(r),len(terms)==2**(r-1))
        counts.append(len(terms))
    FIXTURES.append({'fixture':'independent_noncommuting_inverse_jets','orders':list(range(1,9)),
                     'composition_counts':counts,'matrix_is_actual_charged_operator':False})


def complete_ledger():
    old=json.loads((ROOT/'evidence/seed-middle-coercive-jets.json').read_bytes())['runs'][-1]
    prior=next(f for f in old['fixtures'] if f['fixture']=='complete_angular_and_seed_source_ledger')
    cs=[S.sympify(v) for v in prior['row_constants'][:9]]
    check('interior_rows_cover_eight',len(cs)==9)
    P=[sum(S.binomial(j,a)*cs[a]*cs[j-a] for a in range(j+1)) for j in range(9)]
    A=R(7,20);k=S.Integer(150000)
    V=[S.Integer(0)]+[A*S.factorial(j)*k**j*(1+4*A)**(j-1) for j in range(1,9)]
    pref=82*25*88**4*R(31,2500000)
    kap=[S.Integer(0)]+[pref*sum(S.binomial(j,r)*V[r]*P[j-r] for r in range(1,j+1))/(S.factorial(j)*100**j) for j in range(1,9)]
    G=[8*256**2*88**2*R(31,2500000)*P[j]/(S.factorial(j)*100**j) for j in range(9)]
    check('fixed_gluing_ledger_reconciles',G==[S.sympify(v) for v in prior['gluing_coordinate_bounds'][:9]])
    L=[S.Integer(0)]+[3*10**13*1500**j for j in range(1,9)]
    H=[S.Integer(10)**14]+[G[j]+L[j]+kap[j] for j in range(1,9)]
    gamma=lambda a:R(67,500)*S.factorial(a)*100**a
    B={(0,j):S.factorial(j)*100**j*H[j] for j in range(9)}
    for r in range(8):
        for j in range(8-r):
            B[r+1,j]=B[r,j+1]+r*sum(S.binomial(j,a)*gamma(a)*B[r,j-a] for a in range(1,j+1))
    J=[S.factor(B[j,0]/(S.factorial(j)*100**j)) for j in range(9)]
    check('Christoffel_Cauchy_factor',64*R(17,12)<100)
    check('middle_source_eight_covariant_reserve',max(J)<10**44,[str(S.N(v,18)) for v in J])
    check('old_source_three_allowance_preserved',max(J[:4])<10**26)
    check('withdrawn_sharp_claim_not_restored',J[3]>10**22)
    check('covariant_conversion_keeps_first_two',J[:3]==H[:3])
    zero('covariant_third_conversion',J[3]-H[3]-gamma(1)*H[1]*100/(S.factorial(3)*100**3))
    FIXTURES.append({'fixture':'complete_middle_eight_source_ledger','V':[str(v) for v in V[1:]],
                     'cap_corrections':[str(v) for v in kap],
                     'coordinate_bounds':[str(v) for v in H],
                     'covariant_bounds':[str(v) for v in J],
                     'covariant_bounds_decimal':[str(S.N(v,18)) for v in J],
                     'complete_middle_reserve':'1e44','full_angular_tower_retained':True,
                     'finite_metric_H6_source_evaluated':False})


def independent_connection():
    x,y=S.symbols('x y',real=True);xs=[x,y];at0={x:0,y:0}
    sigma=(x*x+y*y)/100+x*y*y/1000
    f=x+2*y+x*x*y+x*y**3+x**4/7-y**4/5
    G={(a,i,j):(1 if a==i else 0)*S.diff(sigma,xs[j])+(1 if a==j else 0)*S.diff(sigma,xs[i])
                  -(1 if i==j else 0)*S.diff(sigma,xs[a]) for a,i,j in product(range(2),repeat=3)}
    zero('conformal_connection_tensor_norm',sum(v*v for v in G.values())-4*sum(S.diff(sigma,z)**2 for z in xs))
    tensors={1:{(i,):S.diff(f,xs[i]) for i in range(2)}}
    for r in range(1,4):
        nxt={}
        for i in range(2):
            for inds,val in tensors[r].items():
                correction=sum(G[a,i,inds[p]]*tensors[r][inds[:p]+(a,)+inds[p+1:]] for p in range(r) for a in range(2))
                nxt[(i,)+inds]=S.expand(S.diff(val,xs[i])-correction)
        tensors[r+1]=nxt
    errors3=[];errors4=[]
    for i,j,k in product(range(2),repeat=3):
        formula=S.diff(f,xs[i],xs[j],xs[k])-sum(S.diff(G[a,j,k],xs[i])*S.diff(f,xs[a]) for a in range(2))
        errors3.append(S.simplify((tensors[3][i,j,k]-formula).subs(at0)))
    check('independent_third_covariant_formula',all(e==0 for e in errors3))
    for i,j,k,l in product(range(2),repeat=4):
        formula=S.diff(f,xs[i],xs[j],xs[k],xs[l])
        for a in range(2):
            formula-=S.diff(G[a,k,l],xs[i],xs[j])*S.diff(f,xs[a])
            formula-=S.diff(G[a,k,l],xs[j])*S.diff(f,xs[i],xs[a])
            formula-=S.diff(G[a,k,l],xs[i])*S.diff(f,xs[j],xs[a])
            formula-=S.diff(G[a,j,k],xs[i])*S.diff(f,xs[a],xs[l])
            formula-=S.diff(G[a,j,l],xs[i])*S.diff(f,xs[k],xs[a])
        errors4.append(S.simplify((tensors[4][i,j,k,l]-formula).subs(at0)))
    check('independent_fourth_covariant_formula',all(e==0 for e in errors4))
    difference=S.simplify((tensors[3][0,1,1]-tensors[3][1,0,1]).subs(at0))
    check('higher_covariant_tensor_not_symmetrized',difference!=0,difference)
    FIXTURES.append({'fixture':'independent_two_dimensional_connection','third_components':8,'fourth_components':16,
                     'third_permutation_difference':str(difference),'actual_seed_evaluated':False})


def independent_beta_integrals():
    records=[]
    for dps in PROTOCOL['numerical_dps']:
        with mp.workdps(dps):
            c=mp.mpf(3879876)
            for real,imag in [('0','0.003'),('0.5','0.003'),('1','-0.003')]:
                t=mp.mpc(real,imag)
                poly=mp.fsum(mp.binomial(21,i)*t**i*(1-t)**(21-i) for i in range(11))
                integral=1-c*t*mp.quad(lambda s:(t*s)**10*(1-t*s)**10,[0,1])
                error=abs(poly-integral)
                check('independent_complex_beta_'+str(dps)+'_'+real,error<mp.mpf('1e-42'),mp.nstr(error,8))
                records.append({'dps':dps,'real':real,'imag':imag,'absolute_error':mp.nstr(error,15)})
    FIXTURES.append({'fixture':'independent_complex_beta_integrals','records':records,
                     'numerical_samples_replace_uniform_bound':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 middle eight jets','input_count':len(ps)}),flush=True)
    try:
        provenance();cutoff_branches();vertex_budgets();ordered_inverses();complete_ledger();independent_connection();independent_beta_integrals()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 middle eight jets','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'ordinary_actual_middle_seed_C8_bounded':True,'complete_middle_C8_reserve':'1e44',
         'cold_or_outer_C8_bounded':False,'actual_metric_tenth_derivative_bounded':False,
         'actual_finite_metric_source_H6_bounded':False,'actual_transfer_certificate_passed':False,
         'withdrawn_1e22_claim_restored':False,'endpoint_C5_contract_changed':False,
         'actual_signed_compatibility_zero_proved':False,'nonlinear_solution_constructed':False,
         'physical_stability_proved':False,'physical_source_projection_changed':False,
         'physical_boundary_condition_changed':False,'ordinary_proof_machine_verified':False,
         'formal_verification':False,'S167_complete':False,'S168_complete':False,'goal_complete':False}
    if not good:
        for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');saved=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    saved['runs'].append(run);dest.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
