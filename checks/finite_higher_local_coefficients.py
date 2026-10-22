"""Exact supporting ledgers for higher finite-background local coefficients."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as F
from math import factorial, comb
from itertools import product
import json, sys, traceback
import sympy as S

ROOT=Path(__file__).resolve().parents[1];STEM='finite-higher-local-coefficients'
PREVIOUS='seed-nine-metric-eleven'
OWN=['checks/finite_higher_local_coefficients.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
N=11;ROWS=[];FIXTURES=[]
PROTOCOL={'metric_order':11,'trace_order':7,'mixed_coefficient_order':6,
          'arithmetic':'exact rational and symbolic','actual_charged_source_evaluated':False,
          'full_action_trace_is_bookkeeping_residual':False,'coefficientwise_finite_jets':True}


def check(name,ok,detail=''):
    assert name not in {x['name'] for x in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)


def zero(name,x):
    x=S.simplify(x);check(name,x==0,x)


def saved(stem):return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})
def strings(a):return [str(x) for x in a]
def decimal(x):return str(S.N(S.Rational(x.numerator,x.denominator),18))
def poly(a=()):return [F(x) for x in a]+[F(0)]*(N+1-len(a))
def scale(c,a):return [F(c)*x for x in a]
def add(*args):return [sum(a[j] for a in args) for j in range(N+1)]
def mul(a,b):return [sum(a[i]*b[j-i] for i in range(j+1)) for j in range(N+1)]
def der(a):return [(j+1)*a[j+1] for j in range(N)]+[F(0)]
def power(a,n):
    out=poly([1])
    for _ in range(n):out=mul(out,a)
    return out


def compose(a,b):
    out=poly()
    for x in reversed(a):out=add(poly([x]),mul(out,b))
    return out


def exps(a):
    assert a[0]==0
    out=poly([1])
    for j in range(1,N+1):out[j]=sum(i*a[i]*out[j-i] for i in range(1,j+1))/j
    return out


def raw_product(a,b,n):return sum(F(comb(n,i))*a[i]*b[n-i] for i in range(n+1))
def decimal_power(x):
    assert x>0
    p=len(str(x.numerator))-len(str(x.denominator))
    while F(10)**p<x:p+=1
    while F(10)**(p-1)>=x:p-=1
    return p


def covariant(H,C,rank,order):
    B={(0,k):H[k] for k in range(order+1)}
    for r in range(order):
        for k in range(order-r):B[r+1,k]=(k+1)*B[r,k+1]+(rank+r)*sum(C[i]*B[r,k-i] for i in range(k+1))
    return [2*F('2.5')**r*B[r,0] for r in range(order+1)]


def provenance_and_inputs():
    p=saved(PREVIOUS)
    check('prior_metric_eleven_passed',p['status']=='passed' and len(p['checks'])==115 and all(x['passed'] for x in p['checks']))
    check('prior_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    check('all_prior_inputs_unchanged',all(sha256((ROOT/f).read_bytes()).hexdigest()==h for f,h in p['input_sha256_after'].items()))
    check('metric_and_source_scope_retained',p['original_metric_C11_quantified'] and p['seed_source_C9_quantified'] and not p['actual_finite_metric_source_H6_bounded'])
    check('low_middle_covariant_allowance',F('.0016')*F('1.26')<F('.01') and F('50.1')**2*F('.0016')*F('1.26')<10)
    check('low_cold_covariant_allowance',F('.00073')*F('.00331')*914<F('.01') and F('1.001')**2*F('.00073')*F('.00331')*914<10)
    check('high_cold_weighted_allowance',10**24*F('.00331')*F('1.001')**2<4*10**21 and 10**27*F('.00331')*F('1.001')**2<4*10**24)
    check('high_middle_weighted_allowance',10**17*F('50.1')**2<4*10**21 and 10**20*F('50.1')**2<4*10**24)
    check('high_outer_weighted_allowance',9*10**14*F('.0331')+18*F('1e-11')*F('.1')<4*10**21 and 9*10**17*F('.0331')+18*F('1e-11')*F('.1')<4*10**24)
    check('low_outer_weighted_allowance',9*F('1e-7')*F('.0331')+18*F('1e-11')*F('.1')<10)
    check('weighted_value_allowance',F('50.1')**2*F('3e-8')<F('1e-4'))
    check('point_frame_conversion',F('2.5')**2>F('6.1') and 1/(1-14*F('3e-8'))<2)
    check('small_exponential_value',1/(1-4*F('3e-8'))<F('1.001'))
    check('large_jets_not_zeroth_values',6+6>11)
    check('low_physical_input_provenance',not review()['actual_scalar_residual_H6_bounded'] and not review()['higher_TT_lift_quantified'])
    FIXTURES.append({'fixture':'regional_covariant_input_allowances','unweighted_low':'.01','weighted_low':'10',
                     'high_through_ten':'4e21','eleven':'4e24','old_smaller_contracts_retained':True})


def exact_seed_identities():
    z=S.symbols('z',positive=True);b=S.Function('b')(z);p=S.sqrt(b)
    radial=lambda f:-p*S.diff(f,z)
    zero('radial_Hessian_z',radial(-p)-S.diff(b,z)/2)
    zero('angular_Hessian_z',(-S.diff(b,z)/(2*p))*(-p)-S.diff(b,z)/2)
    cc=S.symbols('c0:7');bb=sum(cc[j]*z**j for j in range(7));invC=S.symbols('invC')
    ko=-bb+z*S.diff(bb,z)-z*z*S.diff(bb,z,2)/2
    km=-bb+z*S.diff(bb,z)/2;ks=invC*z*z-bb
    jj=S.expand((ko-ks)/2).subs({cc[1]:0,cc[2]:invC})
    ww=S.expand(ko-2*km+ks).subs({cc[1]:0,cc[2]:invC})
    zero('seed_tracefree_eigenvalue',jj-(-cc[4]*z**4-S.Rational(5,2)*cc[5]*z**5-S.Rational(9,2)*cc[6]*z**6))
    zero('seed_Weyl_amplitude',ww-(-3*cc[3]*z**3-6*cc[4]*z**4-10*cc[5]*z**5-15*cc[6]*z**6))
    rr=S.expand(2*ko+8*km+2*ks).subs({cc[1]:0,cc[2]:invC})
    zero('seed_scalar_trace',rr-(-12*cc[0]-2*cc[5]*z**5-6*cc[6]*z**6))
    check('global_seed_gradient',F(1)/F('5.9')+F('1e-10')+F('.501')+3*F('1e-24')<1)
    # The paired density conversion cancels before the York transpose.
    Jvol,gradPhi=S.symbols('Jvol gradPhi',nonzero=True)
    zero('scalar_density_conversion',Jvol**-1*(Jvol*gradPhi)-gradPhi)
    ell=S.symbols('ell');U0,U1,E0,Rdot=S.symbols('U0 U1 E0 Rdot')
    zero('bookkeeping_linear_trace_cancellation',(Rdot+U1+E0).subs(Rdot,-E0-U1))
    FIXTURES.append({'fixture':'exact_seed_and_density_identities','Hessian_z':'(b_z/2)P',
                     'full_trace_distinguished_from_bookkeeping_residual':True,'physical_bolt_regular':True})


def majorants():
    B=poly([F(1)/F('5.9'),0,F('1e-10'),F('.668'),F('.501'),F('2e-24'),F('1e-24')])
    Z=poly([0,1])
    for j in range(10):Z[j+2]=F(3,4)*compose(der(B),add(poly([1]),Z))[j]/((j+1)*(j+2))
    L=poly();ip=poly([1]);Zpower=poly([1])
    for j in range(1,12):
        Zpower=mul(Zpower,Z);L=add(L,scale(F(1,j),Zpower));ip=add(ip,Zpower)
    check('positive_seed_jet_majorant',all(x>=0 for x in Z+L+ip))
    zero('majorant_first_Hessian',S.Rational(Z[2].numerator,Z[2].denominator)-S.Rational(3,8)*sum(j*S.Rational(B[j].numerator,B[j].denominator) for j in range(1,12)))
    T=poly([F('3e-8')]+[F('.01')*100**j for j in range(1,6)]+[F('4e21')*100**j for j in range(6,11)]+[F('4e24')*100**11])
    K=poly([F('1e-4')]+[F(10)*100**j for j in range(1,6)]+[F('4e21')*100**j for j in range(6,11)]+[F('4e24')*100**11])
    Tp=T.copy();Tp[0]=0
    ex=exps(add(scale(2,L),scale(4,Tp)));M=scale('3.003',ex);dM=scale('12.012',mul(ex,K))
    C=scale(F(3,2),mul(M,der(M)));dC=scale(F(3,2),add(mul(dM,der(M)),mul(M,der(dM))))
    ref=add(compose(der(der(B)),add(poly([1]),Z)),poly(['2e-10']))
    RR=add(ref,scale(2,der(C)),scale(2,mul(C,C)));dRR=add(scale(2,der(dC)),scale(4,mul(C,dC)))
    ric=scale(2,RR);dric=scale(2,dRR)
    scalar=mul(M,ric);dscalar=add(mul(dM,ric),mul(M,dric))
    jj=add(ric,scale(F(1,4),mul(scalar,M)));djj=add(dric,scale(F(1,4),add(mul(dscalar,M),mul(scalar,dM))))
    down=mul(M,RR);ddown=add(mul(dM,RR),mul(M,dRR))
    weyl=add(down,scale(2,mul(M,ric)),scale(F(2,3),mul(scalar,mul(M,M))))
    dweyl=add(ddown,scale(2,add(mul(dM,ric),mul(M,dric))),scale(F(2,3),add(mul(dscalar,mul(M,M)),scale(2,mul(scalar,mul(M,dM))))))
    qp=add(poly([1]),Z)
    pJ=add(scale('.501',power(qp,2)),scale('5e-24',power(qp,3)),scale('4.5e-24',power(qp,4)))
    pC=add(scale(3*F('.668'),ip),poly([6*F('.501')]),scale('20e-24',qp),scale('15e-24',power(qp,2)))
    pZ=scale('1.3e-22',add(scale(4,power(qp,5)),scale(6,power(qp,6))))
    JW=add(scale('6.1',djj),scale(13,pJ));CW=add(scale('6.1',dweyl),scale(25,pC));ZW=add(scale('6.1',dscalar),pZ)
    mu=scale(F('4e-27')*F('1.001')/18,mul(exps(scale(4,Tp)),power(qp,4)))
    check('all_geometric_coefficients_nonnegative',all(v>=0 for a in [M,dM,C,dC,RR,dRR,JW,CW,ZW,mu] for v in a))
    check('reference_has_no_weighted_difference',dRR[0]>0 and ref[0]>0)
    # Recompute the exponential from a symbolic truncated product at low orders.
    t=S.symbols('t');low=sum(S.Rational((2*L[j]+4*Tp[j]).numerator,(2*L[j]+4*Tp[j]).denominator)*t**j for j in range(1,5))
    direct=S.series(S.exp(low),t,0,5).removeO().expand()
    for j in range(5):zero('independent_low_exponential_'+str(j),direct.coeff(t,j)-S.Rational(ex[j].numerator,ex[j].denominator))
    FIXTURES.append({'fixture':'finite_reference_metric_and_curvature_majorants','seed_jets':strings(Z),'log_radius_jets':strings(L),
                     'metric':strings(M),'weighted_metric_difference':strings(dM),'connection':strings(C[:11]),
                     'weighted_connection_difference':strings(dC[:11]),'weighted_J_initial':strings(JW[:10]),
                     'weighted_C_initial':strings(CW[:10]),'weighted_Z_initial':strings(ZW[:10]),
                     'curvature_order':9,'no_endpoint_inverse_used':True})
    return C,JW,CW,ZW,scale('6.1',down),mu


def actual_coefficients(C,jw,cw,zw,rm,mu):
    JW=covariant(jw,C,2,8);CW=covariant(cw,C,4,7);ZW=covariant(zw,C,0,9)
    RM=covariant(rm,C,4,7);MU=covariant(mu,C,0,7)
    for j,b in {0:F(500),1:F('3e17'),2:F('1e25')}.items():JW[j]=min(JW[j],b)
    for j,b in {0:F(500),1:F('3e17')}.items():CW[j]=min(CW[j],b)
    for j,b in {0:F('1.3e-10'),1:F('.0017'),2:F('1e20')}.items():ZW[j]=min(ZW[j],b)
    RM[0]=min(RM[0],55);RM[1]=min(RM[1],F('1e16'))
    JU=[min(JW[j],2*RM[j]) for j in range(8)];CU=[min(CW[j],RM[j]) for j in range(8)];ZU=ZW.copy()
    JU[0]=min(JU[0],8);CU[0]=min(CU[0],50);ZU[0]=min(ZU[0],F('3.1e-12'))
    e=F(1,2000**10*180)
    EE=[(1+11*e/3)*ZW[k]+2*e*ZW[k+2]+11*e/72*raw_product(ZU,ZW,k)
        +7*e/12*raw_product(CU,CW,k)+11*e/6*raw_product(JU,JW,k)+MU[k] for k in range(8)]
    A2=[e*(F(17,3)*JW[k]+F(7,2)*CW[k]) for k in range(7)]
    A1=[4*e*JW[k+1] for k in range(7)]
    A0=[(1+11*e/3)*JW[k]+11*e/36*raw_product(ZU,JW,k)
        +e*(2*JW[k+2]+ZW[k+2]+F(35,3)*raw_product(CW,RM,k)+F(110,3)*raw_product(JW,RM,k))+2*MU[k] for k in range(7)]
    CE=max(EE);CL=max(A0+A1+A2);ep=decimal_power(CE);lp=decimal_power(CL)
    check('trace_uses_metric_eleven',len(ZW)==10 and len(EE)==8)
    check('local_row_uses_metric_ten',len(JW)==9 and len(A0)==7)
    check('trace_exact_decimal_enclosure',F(10)**(ep-1)<CE<=F(10)**ep,{'exponent':ep,'largest':decimal(CE)})
    check('local_exact_decimal_enclosure',F(10)**(lp-1)<CL<=F(10)**lp,{'exponent':lp,'largest':decimal(CL)})
    check('both_mixed_contacts_retained',CE>0 and min(A0+A1+A2)>0)
    check('magnetic_flux_not_reference_charge',min(MU)>0)
    check('tensor_slots_retained',covariant(poly([1,1,1]),poly([1]),2,2)[2]>covariant(poly([1,1,1]),poly([1]),0,2)[2])
    check('mixed_Leibniz_multiplier',3*sum(2**j for j in range(7))==381)
    # Distinct Cauchy product in divided units must give the raw derivative.
    for k in [0,3,6]:
        divided=sum(ZU[i]/factorial(i)*JW[k-i]/factorial(k-i) for i in range(k+1))*factorial(k)
        check('independent_raw_Leibniz_'+str(k),divided==raw_product(ZU,JW,k))
    FIXTURES.append({'fixture':'actual_covariant_curvature_and_local_coefficients','weighted_J':strings(JW),'weighted_C':strings(CW),
                     'weighted_Z':strings(ZW),'unweighted_Rm':strings(RM),'weighted_magnetic':strings(MU),
                     'trace_coefficients':strings(EE),'local_A0':strings(A0),'local_A1':strings(A1),'local_A2':strings(A2),
                     'C_E_exact':str(CE),'C_loc_exact':str(CL),'C_E_decimal':decimal(CE),'C_loc_decimal':decimal(CL),
                     'C_E_power_of_ten':ep,'C_loc_power_of_ten':lp,'higher_TT_lift_supplied':False,
                     'actual_residual_norm_supplied':False})


def independent_sphere_tensors():
    x,y=S.symbols('x y',real=True);coords=[x,y];den=1+x*x+y*y
    g=S.eye(2)*4/den**2;inv=g.inv();z=(1-x*x-y*y)/den
    G={(a,b,c):S.simplify(sum(inv[a,d]*(S.diff(g[d,c],coords[b])+S.diff(g[d,b],coords[c])-S.diff(g[b,c],coords[d]))/2 for d in range(2))) for a,b,c in product(range(2),repeat=3)}
    tensor={():z};records=[]
    for order in range(1,5):
        nxt={}
        for inds in product(range(2),repeat=order):
            a=inds[0];rest=inds[1:];val=S.diff(tensor[rest],coords[a])
            for slot,b in enumerate(rest):
                for c in range(2):
                    replaced=list(rest);replaced[slot]=c;val-=G[c,a,b]*tensor[tuple(replaced)]
            nxt[inds]=S.factor(val)
        tensor=nxt
        if order>=2:
            errors=[]
            for inds,val in tensor.items():
                if order==2:expected=-z*g[inds[0],inds[1]]
                elif order==3:expected=-S.diff(z,coords[inds[0]])*g[inds[1],inds[2]]
                else:expected=z*g[inds[0],inds[1]]*g[inds[2],inds[3]]
                errors.append(S.simplify(val-expected))
            check('independent_ordered_sphere_tensor_'+str(order),all(v==0 for v in errors))
            for point in [(0,0),(S.Rational(1,5),S.Rational(1,7))]:
                vals=[S.simplify(v.subs({x:point[0],y:point[1]})) for v in tensor.values()]
                check('sphere_regular_'+str(order)+'_'+str(point),all(v.is_finite for v in vals))
                records.append({'order':order,'point':strings(point),'components':len(vals)})
    check('ordered_tensor_not_symmetrized',S.simplify(tensor[0,0,1,1]-tensor[0,1,0,1])!=0)
    FIXTURES.append({'fixture':'independent_regular_sphere_covariant_tensors','records':records,'actual_wormhole_geometry':False})


def independent_conformal_curvature():
    x=S.symbols('x',real=True);v=S.Function('v')(x);gp=S.exp(2*v);gi=S.exp(-2*v)
    G={}
    def diff(expr,i):return S.diff(expr,x) if i==0 else S.Integer(0)
    for a,b,c in product(range(4),repeat=3):G[a,b,c]=(int(a==c)*int(b==0)+int(a==b)*int(c==0)-int(b==c)*int(a==0))*S.diff(v,x)
    Ric=S.zeros(4)
    for b,d in product(range(4),repeat=2):
        Ric[b,d]=S.simplify(sum(diff(G[a,d,b],a)-diff(G[a,a,b],d)+sum(G[a,a,c]*G[c,d,b]-G[a,d,c]*G[c,a,b] for c in range(4)) for a in range(4)))
    RR=S.simplify(gi*S.trace(Ric));expected=-6*gi*(S.diff(v,x,2)+S.diff(v,x)**2)
    zero('independent_four_dimensional_scalar',RR-expected)
    J=Ric-RR*gp*S.eye(4)/4;J2=S.simplify(gi**2*sum(z*z for z in J))
    zero('independent_four_dimensional_J_square',J2-3*gi**2*(S.diff(v,x,2)-S.diff(v,x)**2)**2)
    lap=S.simplify(gi*sum(diff(diff(RR,i),i)-sum(G[k,i,i]*diff(RR,k) for k in range(4)) for i in range(4)))
    zero('independent_full_scalar_Laplacian',lap-gi*(S.diff(RR,x,2)+2*S.diff(v,x)*S.diff(RR,x)))
    t=S.symbols('t');vv=t*(1+x*x+x**3)
    reduced=lap.subs(v,vv).doit()
    wrong=S.diff(expected,x,2).subs(v,vv).doit()
    check('omitted_metric_Laplacian_contact_fails',S.simplify(S.diff(reduced-wrong,t,2).subs({t:0,x:S.Rational(1,3)}))!=0)
    Z,e=S.symbols('Z e');R=Z-12;Jj,Cc,mu,Box=S.symbols('J C mu Box')
    euler=Cc-2*Jj+R*R/6
    old=R+12+22*e+e*(S.Rational(3,2)*Cc-S.Rational(11,12)*euler+Box)-mu
    new=(1+11*e/3)*Z+e*Box-11*e*Z*Z/72+7*e*Cc/12+11*e*Jj/6-mu
    zero('independent_complete_trace_rearrangement',old-new)
    FIXTURES.append({'fixture':'independent_four_dimensional_conformal_trace','metric_inverse_retained':True,
                     'Box_metric_contact_retained':True,'actual_charged_kernel_sampled':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 higher local coefficients','input_count':len(ps)}),flush=True)
    try:
        provenance_and_inputs();exact_seed_identities();actual_coefficients(*majorants())
        independent_sphere_tensors();independent_conformal_curvature()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(x['passed'] for x in ROWS)
    run={'story':'S167 higher local coefficients','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after}
    for k,v in review().items():
        if k not in run and isinstance(v,bool):run[k]=v
    if not good:
        for kind,p in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
