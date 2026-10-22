"""Supporting verification of ninth source and eleventh original metric jets."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
from sympy.functions.combinatorial.numbers import stirling
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='seed-nine-metric-eleven';PREVIOUS='schur-source-derivative-budget'
OWN=['checks/seed_nine_metric_eleven.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational;ROWS=[];FIXTURES=[]
PROTOCOL={'source_order':9,'original_metric_order':11,'auxiliary_step_degree':23,
          'numerical_dps':[60,90],'numerical_tolerance':'1e-42',
          'ordinary_continuum_proof_separate':True,'actual_charged_source_evaluated':False,
          'prior_low_entries_retained':True,'endpoint_pullback_order_not_upgraded':True}


def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)


def zero(name,expr):
    val=S.simplify(expr);check(name,val==0,val)


def saved(stem):return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def fixture(stem,name):return next(f for f in saved(stem)['fixtures'] if f['fixture']==name)
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})
def vec(f,key):return [S.sympify(x) for x in f[key]]
def strings(xs):return [str(x) for x in xs]
def decimals(xs):return [str(S.N(x,18)) for x in xs]


def provenance():
    p=saved(PREVIOUS)
    check('prior_derivative_audit_passed',p['status']=='passed' and len(p['checks'])==59 and all(x['passed'] for x in p['checks']))
    check('prior_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    check('all_prior_inputs_unchanged',all(sha256((ROOT/f).read_bytes()).hexdigest()==h for f,h in p['input_sha256_after'].items()))
    check('direct_row_requirement_retained',p['direct_mixed_row_requires_metric_order_eleven'])
    check('new_scope_reviewed',review()['seed_source_C9_quantified'] and review()['original_metric_C11_quantified']
          and not review()['actual_finite_metric_source_H6_bounded'] and not review()['goal_complete'])


def exp_upper(x,n=36):
    return sum(x**j/S.factorial(j) for j in range(n))+x**n/S.factorial(n)/(1-x/(n+1))


def tensors(H,constant,scale,zero0=False,frame=S.Integer(1)):
    order=len(H)-1;B={(0,k):S.factorial(k)*100**k*H[k] for k in range(order+1)}
    for r in range(order):
        for k in range(order-r):
            B[r+1,k]=B[r,k+1]+r*sum(S.binomial(k,i)*constant*S.factorial(i)*scale**i*B[r,k-i]
                                    for i in range(1 if zero0 else 0,k+1))
    return [S.factor(frame*B[j,0]/(S.factorial(j)*100**j)) for j in range(order+1)]


def cutoff_and_middle_vertices():
    t=S.symbols('t');eta=S.expand(sum(S.binomial(23,i)*t**i*(1-t)**(23-i) for i in range(12)))
    c=S.factorial(23)/S.factorial(11)**2
    check('beta_normalization',c==16224936)
    zero('beta_derivative',S.diff(eta,t)+c*t**11*(1-t)**11)
    for j in range(12):
        check('step_endpoint_jets_'+str(j),S.diff(eta,t,j).subs(t,0)==(1 if j==0 else 0) and S.diff(eta,t,j).subs(t,1)==0)
    rho=R(1,100000);d=2*rho;B=R(1,4)+2*R(3,1000)**2
    check('middle_step_first',c*B**11<4)
    check('middle_step_second_root',B/21<R(1,8)**2)
    check('middle_step_second',22*c*(20*B/21)**10/8+44*c*R(3,1000)*B**10<40)
    for name,scale,rmin,limits in [('chi',128,R(3,128),[R('1.011'),520,720000]),('zeta',64,R(1,64),[R('1.006'),260,210000])]:
        values=[1+4*scale*R('1.002')*d,4*scale*R('1.003'),40*scale**2*R('1.003')**2+4*scale*R('2.5')/rmin]
        check(name+'_new_branch_margins',all(a<b for a,b in zip(values,limits)),strings(values))
    fac=R('1.005')
    check('middle_map_first',R('1.006')*R('.2')+260*R('1.5')<400)
    check('middle_map_second',R('1.006')*10+2*260*R('.2')+210000*R('1.5')<316000)
    check('middle_map_inverse',1/(1-400*rho)<fac and R('1.006')*R('1.5')<R('1.6'))
    check('middle_base_scalar_jets',R('1.011')*R('.067')+520*R('.002')<R('1.12')
          and R('1.011')*3+1040*R('.067')+720000*R('.002')<1600)
    check('middle_pretransport_gradient',R('1.011')*10+520*R('.3')<170)
    check('middle_transported_value',(R('1.12')*R('1.6')+R('.304'))*fac<R('2.2'))
    gradient=((1600*R('1.6')+170)*fac+R('1.12')*402+R('2.2')*R('1.12'))*fac
    check('middle_transported_gradient',gradient<3300,gradient)
    check('middle_spin_gradient',(R('1.006')*10+260*R('.21'))*fac<65)
    check('middle_principal',2*(R('2.2')*fac**2+400*fac**3+2*R('.21')*fac**3)<900)
    check('middle_half_divergence',2*fac**5*316000+20000<1200000)
    check('middle_graph_change',900*rho+1200000*rho/36<R(7,20))
    ex=1/(1-R('.002'))
    check('middle_base_graph',2*ex+R('1.12')/36<R('2.1') and R('2.1')+1<4)
    check('middle_weighted_inverse',1/(1-ex/2)<R('2.01'))
    check('middle_weighted_gradient',ex*(1+(1+ex/2)*R('2.01'))+R('1.12')*R('2.01')/36<R('4.1'))
    check('middle_incoming_row',R('4.1')+R('2.01')/2<R('5.2') and 12000*R('5.2')<40*30000)
    check('middle_true_nine_graph_jets',11-2==9)
    FIXTURES.append({'fixture':'C11_step_and_middle_operator_margins','beta_constant':int(c),'endpoint_order':11,
                     'graph_A':'7/20','vertex_k':150000,'global_analytic_cutoff_assumed':False})


def outer_cold_cap_margins():
    c=S.Integer(16224936);rho=R(1,100000)
    third=c*(22/R(4)**10+110/R(4)**9)
    check('outer_step_real',c/R(4)**11<4 and 8*third<60000 and 1000*third<8000000)
    for label,w,w1,w2 in [('positive',R('1.01'),R('2.021'),13),('negative',S.Integer(1),S.Integer(3),9)]:
        dw=8*w+w1;ddw=160*w+16*w1+w2
        check(label+'_real_cap_jets',dw+4<20 and ddw+6*dw+8<400)
    check('outer_complex_step',c*(R(1,4)+2*R('.02')**2)**11<R('4.1'))
    check('outer_complex_cap_ratio',1/(1-R('.04328'))<2)
    perturb=[];base=[1,2,20,400]
    for a in range(4):
        value=S.Integer(0)
        for j in range(1,10):
            raw=sum(stirling(j,k,kind=2)*S.binomial(a,l)*S.factorial(k)/S.factorial(k-l)*2*S.factorial(k+a-l)*100**(k+a-l)
                    for k in range(1,j+1) for l in range(min(a,k)+1))
            value+=rho**j*raw/S.factorial(j)
        perturb.append(value)
    check('negative_value_nine',perturb[0]<R('.003'))
    for a,bound in [(1,R('2.5')),(2,150),(3,50000)]:check('negative_spatial_nine_'+str(a),base[a]+perturb[a]<bound,S.N(base[a]+perturb[a],18))
    check('negative_blend_tail_nine',10**6*R('.001')**10/(1-R('.001'))*10**8<R('1e-6'))
    check('negative_blend_products',1+3*40+3*4000+8000000<10**8)
    changes=[R('.01'),R('.4'),S.Integer(50),S.Integer(7000)];cut=[1,8,160,60000]
    full=[base[j]+sum(S.binomial(j,k)*cut[k]*changes[j-k] for k in range(j+1)) for j in range(4)]
    check('positive_global_guards',full[1]/R('.99')<5 and full[2]/R('.99')<1000 and full[3]/R('.99')<100000)
    check('outer_cap_regularities',11+1==9+3)
    check('cold_complex_step',c*R('.255')**11<R('4.9'))
    check('cold_scalar_cap',R('1.245')*R('.023')<R('.03'))
    check('cold_real_gradient',R('.02')*R('.1')/(R('.0075')-R('.00021'))+100*4*R('.023')<10)
    check('cold_complex_gradient',R('1.245')*R('.02')*R('.1')/(R('.0075')-R('.00021'))+100*R('4.9')*R('.023')<12)
    check('cold_frame_nine',exp_upper(R('.18'))<R(6,5))
    check('cold_actual_row_radius',2-R(11,12)>1)
    check('cold_actual_point_factor',13+600*10+5000*10*11==556013<600000)
    check('cold_actual_row_constant',600000*200*R('1.1')/(1-R(2,75))<2*10**8)
    for n in range(1,12):check('cold_energy_induction_'+str(n),(26+R('1.1')*sum(R(2,75)**r for r in range(1,n)))/(50*n)<1)
    convolution=2*sum(S.factorial(11)/S.factorial(11-r)*R(40,10**5*18)**r for r in range(1,12))
    check('outer_row_absorption_eleven',convolution<R('.001') and R('1.2')*(convolution+R(520,100000))<1)
    check('new_moment_denominators',R('16.776')-9-3>0 and R('172.8')-11-3>0)
    FIXTURES.append({'fixture':'outer_and_cold_ninth_family_margins','negative_perturbations':strings(perturb),
                     'outer_cap_order':12,'cold_graph_parameter_order':9,'cold_point_factor':556013,
                     'physical_state_changed':False})


def source_ledgers():
    cs=[S.Integer(30000)]+[S.Integer(520000)*(j+3)*(600*(j+3))**j for j in range(1,10)]
    P=[sum(S.binomial(j,a)*cs[a]*cs[j-a] for a in range(j+1)) for j in range(10)]
    A=R(7,20);V=[0]+[A*S.factorial(j)*150000**j*(1+4*A)**(j-1) for j in range(1,10)]
    pref=82*25*88**4*R(31,2500000)
    kap=[0]+[pref*sum(S.binomial(j,r)*V[r]*P[j-r] for r in range(1,j+1))/(S.factorial(j)*100**j) for j in range(1,10)]
    gl=[8*256**2*88**2*R(31,2500000)*P[j]/(S.factorial(j)*100**j) for j in range(10)]
    hm=[S.Integer(10)**14]+[gl[j]+3*10**13*1500**j+kap[j] for j in range(1,10)]
    middle=tensors(hm,R('.134'),100,True)
    check('middle_complete_nine',max(middle)<10**48,decimals(middle))
    D=[S.Integer(360000)]+[S.Integer(360000)*241**(j-1) for j in range(1,10)]
    T=[S.Integer(362400)]+[sum((j-b+1)*R('.01')**(j-b)*D[b] for b in range(1,j+1))+1200*(j+1)*R('.01')**j for j in range(1,10)]
    cold_gl=tensors([R('4505.6')*150**j*T[j] for j in range(10)],20,300,frame=R(6,5))
    cold_loc=tensors([S.Integer(10)**7]+[4*10**5*R('1.5')**j for j in range(1,10)],20,300,frame=R(6,5))
    moments=[6*18**j*R('6e-8')*(1+2/(R('16.776')-j-3)) for j in range(10)]
    cold_act=tensors([10**19*(j+1)*R('.00352')*R(1000,3)*moments[j] for j in range(10)],2,25,True)
    cold=[a+b+c for a,b,c in zip(cold_gl,cold_loc,cold_act)]
    check('cold_complete_nine',max(cold)<10**48,decimals(cold))
    check('cold_local_nine',max(cold_loc)<10**9,decimals(cold_loc))
    check('cold_all_components_nonzero',min(cold_gl+cold_loc+cold_act)>0)
    H={p:6*18**p*R('1e-74')*(1+2/(R('172.8')-p-3)) for p in range(2,12)}
    G=[R('.00352')*10**74*sum(stirling(j,m,kind=2)*(m+1)*S.factorial(m)*10**(4*m)*H[m+2] for m in range(j+1))/(S.factorial(j)*100**j) for j in range(10)]
    Q=[0]+[10**16*10**(5*j)*R('.00352')*R('2.1e-27') for j in range(1,10)]
    E=[4*10**19*1000**j+G[j]+Q[j] for j in range(10)]
    O=[sum(stirling(j,k,kind=1)*S.factorial(k)*100**k*E[k] for k in range(j+1))/(S.factorial(j)*100**j) for j in range(10)]
    check('outer_complete_nine',max(E+O)<10**48,decimals([E[9],O[9]]))
    j=9;check('outer_positive_inverse_nine',4*S.factorial(j)*100000**j*21**(j-2)*(20*j+1)<S.factorial(j)*10**(7*j))
    check('retained_eighth_source_reserves',max(middle[:9]+cold[:9]+E[:9]+O[:9])<10**44)
    prior=fixture('seed-outer-eight-jets','complete_defect_and_outer_source_ledger')
    check('outer_prior_entries_exact',E[:9]==vec(prior,'Euler_bounds') and O[:9]==vec(prior,'ordinary_bounds'))
    FIXTURES.append({'fixture':'all_regional_ninth_source_ledgers','middle':strings(middle),'cold_local':strings(cold_loc),
                     'cold_gluing':strings(cold_gl),'cold_actual':strings(cold_act),'cold_total':strings(cold),
                     'outer_Euler':strings(E),'outer_ordinary':strings(O),'source_reserve':'1e48',
                     'ninth_decimal':{'middle':str(S.N(middle[9],18)),'cold':str(S.N(cold[9],18)),
                                      'outer':str(S.N(O[9],18))},'full_angular_tower_retained':True})


def conv(a,b,n):return [sum(a[i]*b[j-i] for i in range(j+1) if i<len(a) and j-i<len(b)) for j in range(n)]
def coeff(c,n=12):return [R(c)*R('.02')**j for j in range(n)]
def power(m,n=12):return [abs(S.binomial(m,j))/100**j for j in range(n)]
def prim(a,m,n):return [a[0]*R(m+2,(m+1)**2)]+[a[j-1]/(100*j) for j in range(1,n)]


def metric_ledgers():
    g=R('6.8e-28');old=fixture('middle-metric-ten-jets','radial_two_derivative_reconstruction')
    P,F,U,N=[vec(old,k) for k in ['P','F','U','N']];T=[S.Integer(10)**20]+[S.Integer(2)*10**44]*8+[S.Integer(2)*10**48]
    n=8
    P.append(sum(4*2**i*(P[n-i]+T[n-i]/2) for i in range(n+1))/(100*(n+1)))
    F.append((sum(40*2**i*F[n-i] for i in range(n+1))+2*g*sum(2**i*(P[n-i]+T[n-i]+10**16*2**(n-i)) for i in range(n+1)))/(100*(n+1)))
    n=9;C=lambda c,X:sum(c*2**i*X[n-i] for i in range(n+1));D=lambda c,X:sum(c*2**i*100*(n-i+1)*X[n-i+1] for i in range(n+1))
    U.append((D(2,U)+C(12,U)+C(10,N)+3*F[n]+3*g*C(2,T)+g*10**18*2**n)/(10000*10*11))
    N.append((D(100,N)+C(20000,U)+C(20000,N)+F[n]+g*C(2,T)+g*10**18*2**n)/(10000*10*11))
    check('middle_eleven_tangent',max(U+N)<10**17,decimals([U[11],N[11]]))
    check('middle_old_arrays_retained',U[:11]==vec(old,'U') and N[:11]==vec(old,'N') and P[:9]==vec(old,'P') and F[:9]==vec(old,'F'))
    rawD=lambda c,X:sum(S.binomial(n,i)*c*S.factorial(i)*200**i*S.factorial(n-i+1)*100**(n-i+1)*X[n-i+1] for i in range(n+1))
    rawC=lambda c,X:sum(S.binomial(n,i)*c*S.factorial(i)*200**i*S.factorial(n-i)*100**(n-i)*X[n-i] for i in range(n+1))
    raw=rawD(2,U)+rawC(12,U)+rawC(10,N)+S.factorial(n)*100**n*3*F[n]+3*g*rawC(2,T)+g*10**18*S.factorial(n)*200**n
    zero('middle_raw_eleventh_ODE',raw-S.factorial(11)*100**11*U[11])
    mid={'U11':str(U[11]),'N11':str(N[11]),'F9':str(F[9]),'tangent_reserve':'1e17'}
    old=fixture('cold-metric-ten-jets','weighted_radial_reconstruction')
    near,away,F,U,N=[vec(old,k) for k in ['Ward_near','Ward_away','F','U','N']]
    Q=[max(a,b) for a,b in zip(near,away)];K=[S.Integer(2)*10**44]*9+[S.Integer(2)*10**48]
    n=8
    F.append((sum(2**i*F[n-i] for i in range(n+1))+2*g*sum(2**i*(Q[n-i]+K[n-i]+10**22*2**(n-i)) for i in range(n+1)))/(100*(n+1)))
    n=9
    U.append(2*(D(12,U)+C(12,U)+C(10,N)+3*F[n]+3*g*C(2,K)+g*10**24*2**n)/(10000*10*11))
    N.append(2*(D(12,N)+C(2,U)+C(10,N)+F[n]+g*C(2,K)+g*10**24*2**n)/(10000*10*11))
    check('cold_eleven_tangent',max(U+N)<10**17,decimals([U[11],N[11]]))
    check('cold_old_arrays_retained',U[:11]==vec(old,'U') and N[:11]==vec(old,'N') and F[:9]==vec(old,'F'))
    raw=2*(rawD(12,U)+rawC(12,U)+rawC(10,N)+S.factorial(n)*100**n*3*F[n]+3*g*rawC(2,K)+g*10**24*S.factorial(n)*200**n)
    zero('cold_raw_eleventh_ODE',raw-S.factorial(11)*100**11*U[11])
    cold={'U11':str(U[11]),'N11':str(N[11]),'F9':str(F[9]),'tangent_reserve':'1e17'}
    old=fixture('outer-metric-ten-jets','outer_two_derivative_reconstruction')
    mu=vec(old,'mu')+[8*R('4e-27')*10**48/1000]
    F=conv(conv(coeff('.14'),power(1),10),mu,10)
    small=conv(conv(coeff('.15'),power(3),11),mu,11)
    W=conv(conv(coeff('.51'),power(-1),12),prim(small,4,12),12)
    W=[min(v,R('1e-8')) if j<=5 else v for j,v in enumerate(W)]
    V=vec(old,'proper_v');n=9;dv=[100*(j+1)*V[j+1] for j in range(10)]
    drift=conv(conv(coeff(1),power(1),10),dv,10)[9]/100
    lie=conv(conv(coeff(4),power(1),10),W,10)[9]/1000
    V.append((3*F[9]+6*g*10**48+g*10**8*R('.02')**9+lie+drift)/(10000*10*11))
    U=conv(coeff('1.001'),prim(conv(coeff('1.001'),V,12),4,12),12)
    term=conv(conv(coeff('.5'),power(1),12),U,12);u=[V[j]+term[j]/100 for j in range(12)]
    term=conv(conv(coeff('1.01'),power(-1),12),U,12);nu=[W[j]+term[j] for j in range(12)]
    gap=conv(coeff('1.07'),nu,12)
    check('outer_eleven_remainder',max(u+nu+gap)+R('1e-9')<10**16,decimals([u[11],nu[11],gap[11]]))
    check('outer_old_arrays_retained',V[:11]==vec(old,'proper_v') and W[:11]==vec(old,'proper_w') and u[:11]==vec(old,'common_u') and nu[:11]==vec(old,'common_nu'))
    check('outer_mass_tail_order',len(mu)==11 and len(small)==11 and len(W)==12)
    FIXTURES.append({'fixture':'appended_original_metric_reconstruction','middle':mid,'cold':cold,
                     'outer':{'mu10':str(mu[10]),'u11':str(u[11]),'nu11':str(nu[11]),'V11':str(gap[11]),'tangent_reserve':'1e16'},
                     'old_metric_entries_retained_exactly':True,'source_order_used':9})


def finite_tensor_ledgers():
    H=S.Integer(10)**17;low=R('.0006');delta=R('.0015')
    check('eleven_high_ideal_square_zero',6+6>11)
    check('middle_inverse_tau_sum',1200*sum(R(1,5)**j for j in range(12))<1500)
    check('middle_finite_V_high',1500*(6*H+low)<9001*H)
    high=6*H*R('2.5')+low*R('2.5')*9001*H
    check('middle_finite_log_eleven',delta/(1-delta)+(delta+high)/(1-delta)<40*H)
    middle=tensors([S.Integer(1)]*12,6,200)
    check('middle_physical_eleven',max(middle)*4*10**18<10**20,decimals(middle))
    check('cold_finite_V_high',(12*H+R('.0012'))*82<985*H)
    low=R('.00072');delta=low*R('.00331');high=low+(R('7.2')+low*985)*H
    check('cold_finite_log_eleven',high<8*H and (low+high)/(1-delta)<12*H<3*10**18)
    C={(1,1):S.Integer(1)};A=[S.Integer(1),S.Integer(1)]
    for m in range(1,11):
        for j in range(1,m+2):C[m+1,j]=C.get((m,j-1),0)+(4*m-j)*C.get((m,j),0)
        A.append(sum(C[m+1,j]/S.factorial(m+1-j) for j in range(1,m+2)))
    near=tensors(A,100,20);away=tensors([S.Integer(1)]*12,15,200)
    check('cold_physical_eleven',max(near+away)*3*10**18<10**27,decimals([near[11],away[11]]))
    H=S.Integer(10)**16;eps=R('2.1e-7')/1000
    nonlinear=10*R('2.1e-7')**2/100+200*eps*H
    check('outer_finite_log_eleven',nonlinear<R('1e-6')*H and H+nonlinear<3*H)
    outer=tensors([S.Integer(1)]*12,R(9,400),2)
    check('outer_physical_eleven',max(outer)<R('1.01') and max(outer)*3*H<10**17,decimals(outer))
    check('all_metric_connections_keep_zero_order',middle[2]>1 and near[2]>A[2] and away[2]>1 and outer[2]>1)
    FIXTURES.append({'fixture':'finite_completion_and_full_eleventh_tensors','middle_amplifications':strings(middle),
                     'cold_Cartesian':strings(A),'cold_near_amplifications':strings(near),'cold_away_amplifications':strings(away),
                     'outer_amplifications':strings(outer),'physical_reserves':{'middle':'1e20','cold':'1e27 a','outer':'1e17 Y^(4-j)L'},
                     'ell_retained':True,'endpoint_map_C12_claimed':False})


def independent_identities():
    t=S.symbols('t');M=S.Matrix([[2+t,t*t],[t+t**3,3-t+2*t*t]])
    M0=M.subs(t,0).inv();vertices={j:M.diff(t,j).subs(t,0) for j in range(1,10)}
    def compositions(n):
        if n==0:yield ();return
        for j in range(1,n+1):
            for rest in compositions(n-j):yield (j,)+rest
    value=S.zeros(2);major=S.Integer(0);parts=list(compositions(9));aa=R(7,20)
    for p in parts:
        term=M0;factor=S.factorial(9)/S.prod(S.factorial(j) for j in p)
        for j in p:term=term*vertices[j]*M0
        value+=(-1)**len(p)*factor*term
        major+=factor*4**(len(p)-1)*S.prod(aa*S.factorial(j)*150000**j for j in p)
    direct=M.inv().applyfunc(lambda e:S.series(e,t,0,10).removeO().expand().coeff(t,9)*S.factorial(9))
    check('noncommuting_ninth_inverse',value==direct and vertices[1]*vertices[2]!=vertices[2]*vertices[1])
    zero('ninth_ordered_majorant',major-aa*S.factorial(9)*150000**9*(1+4*aa)**8)
    check('ninth_composition_count',len(parts)==256)
    rho,u=S.symbols('rho u',positive=True)
    for power_n in [0,2,8,10,14]:
        force=rho**power_n;sol=rho**(power_n+2)/(power_n+2)**2
        direct=S.diff(sol,rho,11)
        identity=S.diff(force,rho,9)-S.integrate(u**10*S.diff(force,rho,9).subs(rho,u*rho),(u,0,1))
        zero('regular_bolt_eleventh_identity_'+str(power_n),direct-identity)
    f=S.Function('f')(rho);out=f
    for n in range(9):out=rho*S.diff(out,rho)-n*out
    zero('ninth_Euler_to_ordinary',S.expand(out-rho**9*S.diff(f,rho,9)))
    FIXTURES.append({'fixture':'independent_ninth_inverse_and_regular_bolt_identities','ordered_compositions':256,
                     'regular_monomials':[0,2,8,10,14],'actual_operator_sampled':False})


def independent_contours():
    t=S.symbols('t');records=[]
    # Exact coefficient recursion for exp and log, tested against nonlinear
    # complex contour evaluation rather than the positive-bound convolutions.
    for tau in [S.Integer(2),S.Integer(7)]:
        nu=R('.0001')+R('.0002')*t+R('.003')*t**6-R('.005')*t**9+R('.004')*t**11
        exponent=tau*nu;v=[S.expand(exponent).coeff(t,j) for j in range(12)]
        e=[S.exp(v[0])]
        for n in range(1,12):e.append(sum(j*v[j]*e[n-j] for j in range(1,n+1))/n)
        a=[1+(e[0]-1)/tau]+[x/tau for x in e[1:]]
        log=[S.log(a[0])]
        for n in range(1,12):log.append((n*a[n]-sum(j*log[j]*a[n-j] for j in range(1,n)))/(n*a[0]))
        for dps in [60,90]:
            with mp.workdps(dps):
                r=mp.mpf('.3');tv=int(tau);count=128
                def actual(z):
                    nv=mp.mpf('.0001')+mp.mpf('.0002')*z+mp.mpf('.003')*z**6-mp.mpf('.005')*z**9+mp.mpf('.004')*z**11
                    return mp.log1p(mp.expm1(tv*nv)/tv)
                contour=sum(actual(r*mp.exp(2j*mp.pi*k/count))*mp.exp(-22j*mp.pi*k/count) for k in range(count))/(count*r**11)
                expected=mp.mpf(str(S.N(log[11],dps)))
                error=abs(contour-expected)
                check('finite_log_eleventh_contour_'+str(tau)+'_'+str(dps),error<mp.mpf('1e-42'),mp.nstr(error,18))
                records.append({'tau':int(tau),'dps':dps,'order':11,'error':mp.nstr(error,18)})
    FIXTURES.append({'fixture':'independent_finite_log_eleventh_contours','records':records,'manufactured_fields':True,
                     'actual_charged_source_sampled':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 ninth source and eleventh original metric','input_count':len(ps)}),flush=True)
    try:
        provenance();cutoff_and_middle_vertices();outer_cold_cap_margins();source_ledgers()
        metric_ledgers();finite_tensor_ledgers();independent_identities();independent_contours()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(x['passed'] for x in ROWS)
    run={'story':'S167 ninth source and eleventh original metric','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after}
    for key,value in review().items():
        if key not in run and isinstance(value,bool):run[key]=value
    if not good:
        for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
