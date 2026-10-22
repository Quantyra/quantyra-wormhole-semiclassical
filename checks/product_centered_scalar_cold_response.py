"""Frozen direct scalar-row barrier; no coupled or nonlinear solution."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as Q
import json,traceback
import sympy as S
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
STEM='product-centered-scalar-cold-response';PREVIOUS='recentered-product-and-cold-residual'
OWN=['notes/'+STEM+'.md','checks/product_centered_scalar_cold_response.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def symbolic_data():
    t,x,d,delta=S.symbols('t x d delta',positive=True)
    fext=(t-1)**2/t**2*(1+x*(t*t+2*t+3)+d/t**2)
    f=S.expand(fext-delta/t)
    R=S.expand(-S.diff(f,t,2)-4*S.diff(f,t)/t-2*(f-1)/t**2)
    rfun=4/t**5-6/t**6
    check('independent_seed_scalar_curvature',S.simplify(R+12*x-d*rfun)==0)
    RA=-S.diff(f,t,2)/2-S.diff(f,t)/t
    RS=(1-f-t*S.diff(f,t))/t**2
    Weyl=(t*t*S.diff(f,t,2)-2*t*S.diff(f,t)+2*f-2)**2/(3*t**4)
    J2=(RA-RS)**2
    P=S.expand(S.Rational(7,12)*Weyl+S.Rational(11,6)*J2-S.Rational(11,72)*(R*R-144*x*x))
    Ps=2*(6*x+d)**2+S.Rational(22,3)*(6*x+d)+S.Rational(22,3)+22*x*x
    check('exact_product_curvature_P',S.simplify(P.subs({t:1,delta:0})-Ps)==0)
    z=S.symbols('z',nonnegative=True)
    remainder=S.expand(P+Ps*rfun/2)
    lapR=S.expand(S.diff(t*t*f*S.diff(R,t),t)/t**2)
    check('extremal_scalar_source_zero',S.simplify(remainder.subs({t:1,delta:0}))==0)
    check('extremal_laplacian_zero',S.simplify(lapR.subs({t:1,delta:0}))==0)
    sig,q=S.symbols('sigma q',positive=True)
    source=d*rfun+sig*(P-10*q*q/t**4)
    centered=-10*sig*q*q*(t-1)*(t+3)/t**6+sig*remainder
    check('magnetic_product_cancellation',S.simplify(source-centered-(d+5*sig*q*q-sig*Ps/2)*rfun)==0)
    ledgers=[]
    for name,expr,degree,target in [('curvature',remainder,2,100000),('laplacian',lapR/d,1,1000000)]:
        pp=S.Poly(S.expand(expr).subs(t,1/z),z,x,d,delta)
        ext=0;thermal=0
        for (nz,nx,nd,nm),c in pp.terms():
            assert nx+nd+nm<=degree and nz>=0
            if nm==0:ext+=abs(c)*abs(nz-4)
            else:thermal+=abs(c)*7**nm
        check('global_source_coefficient_bound_'+name,ext+thermal<target,dict(ext=str(ext),thermal=str(thermal)))
        ledgers.append(dict(channel=name,extremal_sum=str(ext),thermal_sum=str(thermal),allowance=target))
    # Contracted Bianchi relation and the radial conformal Hessian are
    # checked without assuming constant background curvature.
    Jr=S.simplify(RA-R/4)
    expectedJ=-(1+3*x+d)/t**4+d*(5/t**5-S.Rational(9,2)/t**6)
    check('nonconstant_seed_traceless_Ricci',S.simplify(Jr-expectedJ)==0)
    eps,D,lam,l0=S.symbols('epsilon D lambda lambda0')
    ss=S.Function('s')(t);RR=S.Function('R')(t);JJ=S.Function('Jr')(t);ff=S.Function('f')(t)
    lap=lambda u:ff*S.diff(u,t,2)+(S.diff(ff,t)+2*ff/t)*S.diff(u,t)
    div=lambda k:lap(ss)*k+ff*S.diff(k,t)*S.diff(ss,t)
    # div(J grad s)=J:Hess s+(grad R.grad s)/4.
    JH=div(JJ)-ff*S.diff(RR,t)*S.diff(ss,t)/4
    raw=6*D*lap(lap(ss))-6*lap(ss)+(2*RR+48*l0)*ss+S.Rational(11,6)*eps*RR*lap(ss)-S.Rational(22,3)*eps*JH+2*D*div(RR)
    K=6-(2*D+S.Rational(11,6)*eps)*RR+S.Rational(22,3)*eps*JJ
    target=6*D*lap(lap(ss))-div(K)+(2*RR+48*l0)*ss
    check('complete_nonconstant_conformal_Hessian',S.simplify(raw-target)==0)
    Rc,Jc,tc=S.symbols('Rc Jc tc',cls=S.Function)
    Cv=1+S.Rational(11,3)*eps*lam
    sub={RR:-12*lam+Rc(t),JJ:-tc(t)+Jc(t),l0:lam+S.Rational(11,6)*eps*lam*lam}
    kc=-(2*D+S.Rational(11,6)*eps)*Rc(t)+S.Rational(22,3)*eps*Jc(t)
    factored=6*(-D*lap(-lap(ss)+4*lam*ss)+Cv*(-lap(ss)+4*lam*ss))+S.Rational(22,3)*eps*div(tc(t))-div(kc)+2*Rc(t)*ss
    check('exact_factorization_with_lower_terms',S.simplify(target.subs(sub).doit()-factored)==0)
    FIXTURES.append(dict(fixture='independent_invariants_trace_Hessian_and_source_coefficients',ledgers=ledgers))
    return t,x,d,delta,fext,f,R,P,Ps,lapR

def barrier_certificate(data):
    t,x,d,delta,fext,f,*_=data
    h,y,X=S.symbols('h y X',nonnegative=True);V=1+h;a=1-h;p=S.Rational(2,5)
    dv=S.factor((t*fext).subs(t,V))
    lapw=p*a*(f*(p-1)*a+t*(t-a)*S.diff(f,t))/(t*t*(t-a)**2)
    num=S.cancel((S.Rational(337,100)*x-lapw.subs(delta,dv))*t**6*(t-a)**2*V**3)
    check('barrier_positive_denominator_cleared',S.denom(num)==1)
    poly=S.Poly(S.expand(num.subs({t:1+h+y,x:S.Integer(8*10**9)+X})),y,X,h,d)
    grouped={}
    for (ny,nX,nh,nd),c in poly.terms():
        assert nd<=1
        grouped.setdefault((ny,nX),{}).setdefault(nh,[S.Integer(0),S.Integer(0)])[nd]+=c
    hm=S.Rational('1e-30');dm=S.Rational('1e-12');records=[]
    for key,terms in sorted(grouped.items()):
        n0=min(terms);low=S.Integer(0)
        for n,(c0,cd) in terms.items():
            lo=c0-abs(cd)*dm
            if n==n0:low+=lo
            elif lo<0:low+=lo*hm**(n-n0)
        check('global_barrier_coefficient_'+str(key),low>=0,str(low))
        records.append(dict(y_power=key[0],X_power=key[1],factored_h_power=n0,lower_bound_exact=str(low)))
    check('barrier_margin_L4',4-Q('3.37')==Q('.63'))
    check('thermal_mass_delta_guard',(1+Q('1e-30'))**2+2*(1+Q('1e-30'))+3+Q(1,8*10**9)+Q('1e-12')/(8*10**9)<7)
    FIXTURES.append(dict(fixture='global_barrier_rational_coefficient_certificate',parameter_box={'x_min':'8e9','h_max':'1e-30','abs_d_max':'1e-12'},coefficients=records))

def budgets():
    pi=Q('3.14');xlo=Q('8e9');dmax=Q('8.6e-19');Dlam=Q('3e-9');epslam=Q(1,2000**10)/(60*pi)
    kap=Q('.63');H=1+4/kap;m=1-Q('1.1e-8');qratio=(Q('.0035')/6)**2
    check('product_Pstar_bound',(2*(6+Q('1e-12')/xlo)**2+Q(22,3)*(6/xlo+Q('1e-12')/xlo**2)+Q(22,3)/xlo**2+22)<100)
    check('product_d_to_magnetic_ratio',5+50*qratio<Q('5.0001'))
    coef=40+100000*qratio+1000000*Dlam*Q('5.0001')
    check('full_trace_weighted_source_coefficient',coef<41,coef)
    check('full_trace_weighted_source',41/(60*pi)*Q('4e-27')<Q('9e-28'))
    check('LD_barrier_margin',Q('3.37')*Dlam<Q('1.1e-8'))
    R=10*dmax/xlo;R1=56*dmax/xlo;J=10*dmax/xlo;J1=52*dmax/xlo
    k=(2*Dlam+Q(11,6)*epslam)*R+Q(22,3)*epslam*J
    k1=(2*Dlam+Q(11,6)*epslam)*R1+Q(22,3)*epslam*J1
    tmax=3+1/xlo+Q('1e-12')/xlo
    defect=(Q(22,3)*epslam*Q(7,3)*H*tmax+H*(k+k1/3)+2*R/kap)/(6*m)
    check('actual_scalar_weighted_transfer_defect',defect<Q('1e-26'),defect)
    U=Q('9e-28')/(6*m*(1-Q('1e-26')))
    check('direct_scalar_value',U/kap<Q('2.4e-28'))
    check('cold_disk_gradient',Q('1.06')*Q('1.003')/Q('.99')*H/Q('5.99')*U<Q('2.4e-28'))
    hmin=Q(1,8*2000**10);boost5=4/(64**5*hmin**3)
    check('cold_admission_two_components',(Q('2.4e-28')/Q('1.2e-9'))**5*boost5<1)
    general=1/(6*m*(1-Q('1e-26'))*kap)
    check('general_scalar_cold_response_factor',(general/(Q(4,3)*10**18))**5*boost5<1)
    check('mixed_source_threshold',(Q(4,3)*10**18)*(Q('9e-28')+Q('2.8e-27'))<Q('5e-9'))
    check('cold_weight_ratio_decreasing',Q(2,5)-1+Q(1,3)<0)
    check('Lh_min_guard',Q(1,2000**10)<Q('1e-30'))
    FIXTURES.append(dict(fixture='rational_actual_scalar_inverse_and_cold_budgets',weighted_defect_exact=str(defect),
        scalar_value_upper='2.4e-28 times w',cold_01_upper='1.2e-9',mixed_source_target='2.8e-27 times lambda in weighted value norm',coupled_response_bounded=False))

def numerical_fixtures(data):
    t,xx,dd,dt,fext,f,R,P,Ps,lapR=data
    pfun=S.lambdify((t,xx,dd,dt),P,'mpmath');lapfun=S.lambdify((t,xx,dd,dt),lapR,'mpmath')
    records=[]
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    eta_shift=mp.mpf(old['fixtures'][3]['records'][0]['sphere_squared_first_shift'])
    for digits in [90,140]:
        mp.mp.dps=digits;k=mp.mpf(2000);q=k**4;g0=k**-10
        g=2*g0/(1+mp.sqrt(1+44*g0/(360*mp.pi)));F=1-384*mp.pi*mp.mpf('1e22')*g
        C0=mp.sqrt(1+12*mp.pi*k**6);x0=(C0-1)*g/g0/6
        Crn=mp.sqrt(1+12*mp.pi*g*k**16/F);x=(Crn-1)/6*(1+eta_shift)
        sig=g/(60*mp.pi*x);Dbar=(96*mp.pi*mp.mpf('1e22')*g-g/(60*mp.pi))/x
        aa=sig;bb=sig*(12*x+mp.mpf(11)/3)-1;cc=sig*(47*x*x+22*x+mp.mpf(11)/3-5*q*q)
        d=2*cc/(-bb+mp.sqrt(bb*bb-4*aa*cc));C=1+6*x+d
        Pstar=2*(6*x+d)**2+mp.mpf(22)/3*(6*x+d)+mp.mpf(22)/3+22*x*x
        check('exact_trace_root_'+str(digits),abs((d+5*sig*q*q-sig*Pstar/2)/d)<mp.mpf('1e-80'))
        for zi in ['.25','.265','.28']:
            kap=mp.mpf(zi)*2*k**3*g*C0/x0**mp.mpf('1.5');gamma=kap*mp.sqrt(x)/C
            def K(h):
                v=1+h
                return ((1+3*x)*(v**-2+v**-3)+3*x*(1+v**-1)-d*(v**-4-3*v**-5))/(2*C)
            h=gamma
            for _ in range(5):h-=(h*K(h)-gamma)/(K(h)+h*mp.diff(K,h))
            V=1+h;delta=h*h/V*(1+x*(V*V+2*V+3)+d/V**2)
            wbolt=(2*h/V)**mp.mpf('.4');Lh=1+mp.log(mp.mpf('.01')/h)
            cold=mp.mpf('2.4e-28')*wbolt/(h*Lh)
            check('benchmark_cold_bound_'+str((digits,zi)),cold<mp.mpf('1.2e-9'))
            source_rows=[]
            for offset in [h,10*h,mp.mpf('1e-10'),mp.mpf('.001'),mp.mpf('.1'),mp.mpf('1'),mp.mpf('9')]:
                tt=1+offset;w=((offset+h)/tt)**mp.mpf('.4')
                rfun=4/tt**5-6/tt**6
                raw=d*rfun+sig*(pfun(tt,x,d,delta)-10*q*q/tt**4)-Dbar*lapfun(tt,x,d,delta)
                centered=-10*sig*q*q*(tt-1)*(tt+3)/tt**6+sig*(pfun(tt,x,d,delta)+Pstar*rfun/2)-Dbar*lapfun(tt,x,d,delta)
                check('independent_centered_source_'+str((digits,zi,str(offset))),abs(raw-centered)<mp.mpf('1e-80'))
                check('sample_full_trace_bound_'+str((digits,zi,str(offset))),abs(centered)/(x*w)<mp.mpf('9e-28'))
                source_rows.append(dict(offset=mp.nstr(offset,25),source_over_lambda_w=mp.nstr(centered/(x*w),45)))
            records.append(dict(digits=digits,prepared_zi=zi,h=mp.nstr(h,60),d=mp.nstr(d,60),
                scalar_bolt_upper=mp.nstr(mp.mpf('2.4e-28')*wbolt,60),cold_01_upper=mp.nstr(cold,60),source_rows=source_rows))
    for i in range(3):
        for key in ['h','d','scalar_bolt_upper','cold_01_upper']:
            check('cross_precision_'+str(i)+'_'+key,abs(mp.mpf(records[i][key])/mp.mpf(records[i+3][key])-1)<mp.mpf('1e-50'))
    FIXTURES.append(dict(fixture='benchmark_trace_and_cold_bounds',records=records,
        numerical_geometry='sphere radius from preceding illustrative pressure; d solves exact trace at that radius; analytic theorem covers entire guarded box',actual_coupled_solution_sampled=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_product_and_horizon_result',old['status']=='passed' and old['recentered_product_enclosure_proved'] and old['product_centered_horizon_value_target_met'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        data=symbolic_data();barrier_certificate(data);budgets();numerical_fixtures(data)
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,
        unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        product_centered_scalar_row_weighted_inverse_proved=error is None,
        product_centered_scalar_row_cold_01_bound='1.2e-9' if error is None else None,
        actual_coupled_cold_response_bounded=False,admission_norm_inverse_quantified=False,
        inverse_transferred_to_product_centered_seed=False,nonlinear_solution_constructed=False,
        semiclassical_validity_proved=False,physical_stability_proved=False,physical_coupled_inverse_proved=False,
        component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
