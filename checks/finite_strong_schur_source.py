"""Finite-comparison tensor regularity and actual strong source budgets."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
STEM='finite-strong-schur-source'
PREVIOUS='cold-projection-value-upper-bound'
OWN=['checks/finite_strong_schur_source.py','notes/'+STEM+'.md','evidence/finite-strong-schur-source-review.json']
R=S.Rational
ROWS=[]
FIXTURES=[]
PROTOCOL={'arithmetic':'exact symbolic and rational', 'general_warp':True,
          'strong_vertex_upper':'2e4/sqrt(g)', 'actual_source_norm_upper':'2e28 ell^2',
          'actual_signed_moment_evaluated':False, 'actual_transfer_certificate_passed':False}


def check(name, ok, detail=''):
    assert name not in {c['name'] for c in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)


def zero(name, expression):
    value=S.factor(S.expand(S.simplify(expression)))
    check(name,value==0,value)


def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})


def provenance():
    p=prior()
    check('predecessor_passed',p['status']=='passed' and len(p['checks'])==64 and all(c['passed'] for c in p['checks']))
    check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('predecessor_hashes',not bad,bad)
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('retained_inputs_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
    check('scope_retained',review['actual_weighted_starting_source_bounded'] and not review['actual_transfer_certificate_passed'] and not review['goal_complete'])


def general_return():
    A,rr,b,c,eta,c1,c2,K,K1,r1,r2=S.symbols('A R b c eta c1 c2 K K1 R1 R2')
    rules={A:K-A*A,rr:r1,r1:r2,b:-2*A*(b+c)-2*rr*(b-c),c:c1,c1:c2,eta:A*eta+b+c,K:K1}
    diff=lambda p:S.expand(sum(S.diff(p,v)*dv for v,dv in rules.items()))
    E1=K-r1-A*rr
    Tc=-(b+3*c)/2-(A-rr)*eta
    first=-3*c1/2+2*rr*b-E1*eta
    second=-3*(c2+A*c1)/2+(2*r1-2*A*rr-4*rr**2-E1)*b+(4*rr**2-4*A*rr-E1)*c-(diff(E1)+2*A*E1)*eta
    zero('general_return_first',diff(Tc)-first)
    zero('general_return_second',diff(diff(Tc))+A*diff(Tc)-second)
    rough=diff(rules[b])+(A+2*rr)*rules[b]-4*A*A*(b+c)-4*rr**2*(b-c)
    expected=-2*K*(b+c)-2*(r1+2*rr**2-A*rr)*(b-c)-2*(A-rr)*c1
    zero('general_TT_rough_substitution',rough-expected)
    seed_only=-2*K*(b+c)-2*rr**2*(b-c)-2*(A-rr)*c1
    zero('seed_substitution_defect',expected-seed_only+2*(r1+rr**2-A*rr)*(b-c))
    check('seed_only_identity_rejected',S.expand(expected-seed_only).subs({A:2,rr:1,r1:3,b:4,c:1,c1:2,K:5})!=0)
    rho=S.symbols('rho',positive=True)
    f,theta,w=[S.Function(n)(rho) for n in ['f','theta','w']]
    lap=lambda v:S.diff(v,rho,2)+w*S.diff(v,rho)
    chi=[S.cos(theta),S.sin(theta)]
    zero('finite_partition_identity',sum(v*lap(v*f) for v in chi)+S.diff(theta,rho)**2*f-lap(f))
    C,B,t,measure,L=S.symbols('C B t measure L',positive=True)
    ell2=B**2/C
    zero('finite_orbit_measure_conversion',C*(measure/(4*S.pi*B**2*t**4*ell2))*(ell2*t**2*L)**2-measure*L**2/(4*S.pi))
    FIXTURES.append({'fixture':'general_return_and_measure','seed_Rprime_identity_used':False})


def coordinate_rough():
    rho=S.symbols('rho',positive=True)
    time,angle,azimuth=S.symbols('time angle azimuth',real=True)
    coords=[time,rho,angle,azimuth]
    aa,rr,bb,cc=[S.Function(n)(rho) for n in ['a','r','b','c']]
    metric=[aa**2,S.Integer(1),rr**2,rr**2*S.sin(angle)**2]
    inv=[1/q for q in metric]
    frame=[-bb-2*cc,bb,cc,cc]
    H=[metric[i]*frame[i] for i in range(4)]
    Gamma={}
    for k in range(4):
        for i in range(4):
            for j in range(4):
                expr=((S.diff(metric[j],coords[i]) if k==j else 0)+(S.diff(metric[i],coords[j]) if k==i else 0)-(S.diff(metric[i],coords[k]) if i==j else 0))*inv[k]/2
                Gamma[k,i,j]=S.simplify(expr)
    nH={}
    for j in range(4):
        for a in range(4):
            for d in range(4):
                nH[j,a,d]=S.simplify((S.diff(H[a],coords[j]) if a==d else 0)-Gamma[d,j,a]*H[d]-Gamma[a,j,d]*H[a])
    Avec,Rvec=S.diff(aa,rho)/aa,S.diff(rr,rho)/rr
    scalar=lambda p:S.diff(p,rho,2)+(Avec+2*Rvec)*S.diff(p,rho)
    rb=scalar(bb)-4*Avec**2*(bb+cc)-4*Rvec**2*(bb-cc)
    rc=scalar(cc)+2*Rvec**2*(bb-cc)
    targets=[-rb-2*rc,rb,rc,rc]
    for a in range(4):
        for d in range(4):
            value=0
            for i in range(4):
                v=S.diff(nH[i,a,d],coords[i])
                v-=sum(Gamma[k,i,i]*nH[k,a,d]+Gamma[k,i,a]*nH[i,k,d]+Gamma[k,i,d]*nH[i,a,k] for k in range(4))
                value+=inv[i]*v
            zero(f'coordinate_rough_{a}_{d}',value-(targets[a]*metric[a] if a==d else 0))
    zero('coordinate_trace',sum(targets))
    FIXTURES.append({'fixture':'independent_four_dimensional_rough_laplacian','actual_geometry_evaluated':False})


def curvature_ledgers():
    m1,m2,m3,m4=S.symbols('M1 M2 M3 M4')
    G0=3*m1/2;G1=3*m2/2+6*m1**2
    G2=3*m3/2+18*m1*m2+48*m1**3
    G3=3*m4/2+24*m1*m3+18*m2**2+288*m1**2*m2+576*m1**4
    Q0=2*G1+8*G0**2;Q1=2*G2+16*G0*G1;Q2=2*G3+16*G0*G2+16*G1**2
    tensor=128*(Q2+(8*m1+36*G0)*Q1+(4*m2+16*G1+144*m1*G0+320*G0**2)*Q0)
    target=384*m4+34560*m1*m3+19968*m2**2+1241088*m1**2*m2+6236160*m1**4
    zero('tensor_Hessian_component_polynomial',tensor-target)
    check('tensor_Hessian_coefficient_reserve',max(S.Poly(target,m1,m2,m3,m4).coeffs())<10**7)
    p3=1117000;p4=442510000
    check('unweighted_Riemann_gradient_core',R('2.5')**3*10**6*p3<10**16)
    check('unweighted_Riemann_gradient_outer',1000*R('2.5')**3*10**6*p3/40**2+R(1000,40**3)<10**16)
    check('weighted_J_value_core',7+R('.001')*40**2<500)
    check('weighted_J_value_outer',R(12000,40)+R(7,40**2)<500)
    check('weighted_Weyl_value_core',40+R('.001')*40**2<500)
    check('weighted_Weyl_value_outer',R(6040,40)<500)
    check('weighted_curvature_gradient_core',40**3*R('1e15')**2<R('3e17')**2)
    check('weighted_curvature_gradient_outer',1000*R('2.5')**3*10**6*p3+R(1000,40)<R('3e17'))
    for region,factor in [('core',1600),('outer',1000)]:
        check('weighted_scalar_Hessian_'+region,factor*R('2.5')**4*10**6*p4+4<R('1e20'))
        check('weighted_tensor_Hessian_'+region,factor*R('2.5')**4*10**7*p4+10**6<R('1e25'))
    M=[R(1),R(30),R(3000),R(10**6),R(4*10**8),R(10**16)]
    N=[R(1)]
    for j in range(1,5):N.append(4*sum(S.binomial(j,k)*M[k]*N[j-k] for k in range(1,j+1)))
    G=[6*sum(S.binomial(j,k)*N[k]*M[j+1-k] for k in range(j+1)) for j in range(5)]
    Q=[2*G[j+1]+8*sum(S.binomial(j,k)*G[k]*G[j-k] for k in range(j+1)) for j in range(4)]
    SS=[64*sum(S.binomial(j,k)*N[k]*Q[j-k] for k in range(j+1)) for j in range(4)]
    B5=32*(N[1]*SS[2]+SS[3]+4*G[1]*SS[1]+4*G[0]*SS[2]+4*N[1]*G[0]*SS[1])
    check('fifth_jet_component_budget',B5<R('2.5e20'),B5)
    check('weighted_fifth_core',1600*R('2.5')**5*B5+10**5<R('1e27'))
    check('weighted_fifth_outer',20000*R('2.5')**5*B5+10**5<R('1e27'))
    FIXTURES.append({'fixture':'exact_curvature_component_ledgers','tensor_Hessian_polynomial':str(target),'B5':str(B5)})


def outer_trace_budget():
    t=R(40);q=R('1.001')
    U,V,n,w=R('8e-9')/t**2,R('2e-9')/t**2,R('1e-8')/t**2,R('1.1e-8')/t**2
    uy,ny,wy=R('3.2e-7')/t,R('1.5e-8')/t**2,R('3.3e-7')/t
    vy,vyy=R('4.3e-8')/t**2,R('2.1e-6')/t**2
    d,dy=R('2.01e-9')/t**3,R('4.6e-8')/t**3
    B,By=R('2.1e-18')/t**4,R('8.7e-17')/t**4
    sy,ty=R('1.6e-8')/t**2,q/t
    lnu,lnuy,ltau,ltauy=R('6.3e-7'),R(500),22/t,146/t
    rv,ry,dr,dry=R(13),R('2e-20')/t**5,R('1e-20')/t**4,R('3e-20')/t**4
    er=R('2.01')*(n*n+R('.502')*V*V/t)
    ery=2*q*(q*(d*ny+B*ty)+2*n*wy)
    vals=[rv*er,2*q*w*dr,6*q*d*lnu,6*q*q*B*ltau,6*q*q*vy**2,4*q*q*U**2/(1000*t**2)]
    ders=[ry*er+rv*ery,2*q*(wy*dr+w*dry),
          6*q*((2*wy*d+dy)*lnu+d*lnuy),
          6*q*((2*wy*q*B+q*(By+B*sy))*ltau+q*B*ltauy),
          6*q*(2*wy*q*vy**2+dy*vy**2+q*ty*vy**2+2*q*vy*vyy),
          2*q/(1000*t**2)*(R('4.004')*U*uy+2*(wy+q)*R('2.002')*U**2)]
    val,der=t*t*sum(vals),t*t*sum(ders)
    # Every term after multiplying by t^2 is a positive nonincreasing power of t.
    check('outer_curvature_remainder_value',val<R('1e-12'),val)
    check('outer_curvature_remainder_gradient',der<R('2e-7'),der)
    check('core_weighted_Z_value',40*R('3.1e-12')<R('1.3e-10'))
    check('core_weighted_Z_gradient',40*R('4e-5')<R('.0017'))
    check('outer_weighted_Z_value',val+R('1e-20')<R('1.3e-10'))
    check('outer_weighted_Z_gradient',R('1.000001')*(der+R('1e-20'))<R('.0017'))
    e=R('4e-17')**2/180
    tracevalue=(1+11*e/3)*R('1.3e-10')+e*(R('2e20')+11*R('1.3e-10')**2/72+R(7,12)*500**2+R(11,6)*500**2)+R('1e-26')
    tracegrad=(1+11*e/3)*R('.0017')+e*(R('1e27')+R('1e20'))+R('1e-24')
    check('weighted_full_trace_value',tracevalue<R('2e-10'),tracevalue)
    check('weighted_full_trace_gradient',tracegrad<R('.002'),tracegrad)
    check('differentiated_curvature_square_reserve',R(7,6)*50*R('3e17')+R(11,3)*8*R('3e17')+R(11,36)*R('1.3e-10')*R('.0017')<R('1e20'))
    FIXTURES.append({'fixture':'exact_outer_remainder_and_contact','weighted_N_upper':str(val),'weighted_gradient_N_upper':str(der)})


def final_budgets():
    sg=R('4e-17');g=sg**2
    check('positive_form_localization_loss',R(500)**4/(100*3)*R('3.4e-37')<R('1e-28'))
    check('physical_orbit_margin',800*R(22,7)**2<89**2)
    check('finite_c_value',8*R(22,7)/R('.90')<R('5.3')**2)
    check('finite_tensor_value',25**2/R('.90')<27**2)
    check('finite_tensor_gradient',420**2/R('.90')<443**2)
    check('finite_Aeta',R('.713')*27<20)
    check('finite_divergence',160**2/R('.90')<170**2)
    check('finite_b_value',R('2.002')*R('5.3')<11)
    check('finite_radial_curvature_derivative',3*R('1e16')+2*R('2.25')*35<R('4e16'))
    qc=R(2,3)*(89+(R('166.4')*11+R('120.4')*R('5.3')+104*20)*sg+R('4e16')*27*sg)+R('4.5')*50*sg
    check('finite_scalar_laplacian',qc<100,qc)
    rough=7*100+(26*50+1000*R('5.3'))*sg
    check('finite_tensor_rough',rough<701,rough)
    hess2=701**2+(15+4*55)*443**2*g+4*R('1e16')*27*443*g
    check('finite_tensor_Hessian',hess2<800**2,hess2)
    check('local_second_derivative_coefficient',2*500*(1+R(11,6))+R(7,12)*6*500<5000)
    vcoef=R('1.2')*(501*27+sg/180*5000*800+g/180*(R('6e17')*443+R('4e25')*27)+R('2e-10')*170+R('.002')*27)
    check('strong_weighted_vertex',vcoef<20000,vcoef)
    check('TT_source_factor',80*R(22,7)**2*8*R(22,7)/R('.90')<150**2)
    check('scalar_source_factor',80*R(22,7)**2<800)
    zero('Schur_TT_source_factor',20000*150-3*10**6)
    source=R('3e6')*27*R('2.2e-13')/R('3e-17')**2
    check('actual_starting_source_bound',source<R('2e28'),source)
    check('starting_norm_not_terminal_certificate',R('2e28')/25>R('1e-49'))
    FIXTURES.append({'fixture':'exact_form_and_strong_source_budgets','vertex_coefficient':str(vcoef),'source_upper':str(source),
                     'actual_signed_moment_evaluated':False,'actual_transfer_certificate_passed':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 finite strong Schur source','input_count':len(ps)}),flush=True)
    try:provenance();general_return();coordinate_rough();curvature_ledgers();outer_trace_budget();final_budgets()
    except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
    run={'story':'S167 finite strong Schur source','created_utc':datetime.now(timezone.utc).isoformat(),'status':'passed' if good else 'failed',
         'exception':error,'python':sys.version,'sympy':S.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'finite_TT_second_form_bound_proved':good,'finite_strong_weighted_vertex_proved':good,
         'actual_weighted_starting_source_bounded':good,'actual_signed_moment_evaluated':False,'actual_transfer_certificate_passed':False,
         'actual_weighted_graph_iterates_certified':False,'actual_signed_compatibility_zero_proved':False,
         'cold_metric_response_proved':False,'nonlinear_neighborhood_preserved':False,'nonlinear_solution_constructed':False,
         'physical_source_projection_changed':False,'physical_stability_proved':False,'formal_verification':False,
         'S167_complete':False,'S168_complete':False,'goal_complete':False}
    if not good:
        for kind,p in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');saved=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    saved['runs'].append(run);dest.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
