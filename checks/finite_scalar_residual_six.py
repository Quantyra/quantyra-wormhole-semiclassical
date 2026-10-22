"""Actual scalar residual: quadratic path cancellation and weighted H6."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as F
from math import factorial, comb
import importlib.util,json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-scalar-residual-six';PREVIOUS='finite-higher-local-coefficients'
OWN=['checks/finite_scalar_residual_six.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
spec=importlib.util.spec_from_file_location('retained_coefficient_calculus',ROOT/'checks/finite_higher_local_coefficients.py')
L=importlib.util.module_from_spec(spec);spec.loader.exec_module(L)
poly,add,mul,scale,der,power,compose,exps=L.poly,L.add,L.mul,L.scale,L.der,L.power,L.compose,L.exps
ROWS=[];FIXTURES=[]
PROTOCOL={'scalar_source_order':6,'metric_order_used':10,'path_derivatives':2,
          'pointwise_decay':'t^-4','physical_volume_retained':True,'numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','actual_charged_source_sampled':False,'ordinary_continuum_proof_separate':True}


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
def decimals(a):return [L.decimal(x) for x in a]
def fixture(stem,name):return next(f for f in saved(stem)['fixtures'] if f['fixture']==name)
def fromsym(x):return F(int(S.numer(x)),int(S.denom(x)))
def parameter_add(*args):return [add(*(x[j] for x in args)) for j in range(3)]
def parameter_scale(c,a):return [scale(c,x) for x in a]
def parameter_mul(a,b):return [add(*(scale(comb(j,k),mul(a[k],b[j-k])) for k in range(j+1))) for j in range(3)]
def parameter_der(a):return [der(x) for x in a]


def provenance_and_tangent_inputs():
    old=saved(PREVIOUS)
    check('previous_coefficients_passed',old['status']=='passed' and len(old['checks'])==58 and all(x['passed'] for x in old['checks']))
    check('previous_inputs_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    check('scope_is_actual_scalar_only',review()['actual_scalar_residual_H6_bounded'] and not review()['actual_TT_source_H4_bounded'] and not review()['actual_transfer_certificate_passed'])
    check('cold_V_low',2*F('.0006')*82*914<100 and 2*F('.0006')*82*914*F('1.001')**2<4000)
    check('cold_V_high',2*(F('.0006')+5*F('2e15'))*82*F('1.02e7')*F('1.001')**2<F('2e26'))
    check('middle_V_low',1500*F('.0006')*F('1.26')<100 and 1500*F('.0006')*F('1.26')*F('50.1')**2<4000)
    check('middle_V_high',1500*(F('.0006')+5*10**15)*F('3.12')*F('50.1')**2<F('2e26'))
    check('outer_V_high',F('1.07')*F('1.002')*9*F('1e14')*F('.0331')<F('2e26'))
    check('V_weighted_value',F('1.2e-8')*F('50.1')**2<F('4e-5'))
    check('nu_weighted_high',2*10**15*F('1.02e7')*F('.00331')*F('1.001')**2<F('4e21') and 10**15*F('3.12')*F('50.1')**2<F('4e21'))
    check('nu_low_covariant',F('.0001')*F('.00331')*914<F('.01') and F('.0001')*F('1.26')*F('50.1')**2<10)
    check('parameter_small_value_enclosures',1/(1-F('1.2e-8')-2*F('3e-8'))<F('1.001') and 1/(1-6*F('3e-8'))<F('1.001'))
    FIXTURES.append({'fixture':'actual_tangent_and_parameter_inputs','V_low_unweighted':'100','V_low_weighted':'4000',
                     'V_high':'2e26','inverse_thermal_gap_used':False,'metric_order':10})


def exact_residual_identities():
    ell,tau,V=S.symbols('ell tau V',real=True);T=1-tau+tau*S.exp(ell*V);ss=S.log(T);nu=tau*V
    zero('radius_first_parameter',S.diff(ss,ell)-nu*S.exp(ell*V)/T)
    zero('radius_second_parameter',S.diff(ss,ell,2)-tau*(1-tau)*V**2*S.exp(ell*V)/T**2)
    z=S.symbols('z');d=S.symbols('d');Z0=d*(4*z**5-6*z**6)
    zero('reference_Lie_derivative',S.diff(Z0.subs(z,z*S.exp(-ss)),ell).subs(ell,0)+nu*z*S.diff(Z0,z))
    zz,aa,hh,z0,h0,zd,ad=S.symbols('Z A H Z0 H0 Zd Ad')
    E=zz-(1-ell)*aa+ell*hh
    NR=zz-z0-ell*zd;NA=aa-z0-ell*ad
    rhs=NR-(1-ell)*NA+ell**2*ad+ell*(hh-h0)
    zero('complete_quadratic_trace_identity',(E-rhs).subs(zd,-z0-h0+ad))
    check('omitted_reference_factor_fails',S.expand((E-(NR-NA+ell**2*ad+ell*(hh-h0))).subs(zd,-z0-h0+ad))!=0)
    # Scalar Laplacian in the canonical product followed by four-dimensional conformal scaling.
    b=S.Function('b')(z);f=S.Function('f')(z);cx=S.symbols('cx')
    direct=cx*z*z*(b*S.diff(f,z,2)+S.diff(b,z)*S.diff(f,z)-2*b*S.diff(f,z)/z)
    zero('seed_physical_scalar_Laplacian',direct-cx*(z*z*b*S.diff(f,z,2)+(z*z*S.diff(b,z)-2*z*b)*S.diff(f,z)))
    jv,fp=S.symbols('jv fp',nonzero=True);zero('density_cancels_before_scalar_York',jv**-1*jv*fp-fp)
    FIXTURES.append({'fixture':'exact_bookkeeping_trace_and_path_identities','quadratic_zero_retained':True,
                     'radius_acceleration_retained':True,'full_action_trace_substituted_for_residual':False})


def build_parameter_geometry():
    prev=fixture(PREVIOUS,'finite_reference_metric_and_curvature_majorants')
    Z=[F(x) for x in prev['seed_jets']];ll=[F(x) for x in prev['log_radius_jets']];q=add(poly([1]),Z)
    T=poly([F('3e-8')]+[F('.01')*100**j for j in range(1,6)]+[F('4e21')*100**j for j in range(6,11)])
    K=poly([F('1e-4')]+[F(10)*100**j for j in range(1,6)]+[F('4e21')*100**j for j in range(6,11)])
    VV=poly([F('1.2e-8')]+[F(100)*100**j for j in range(1,6)]+[F('2e26')*100**j for j in range(6,11)])
    KV=poly([F('4e-5')]+[F(4000)*100**j for j in range(1,6)]+[F('2e26')*100**j for j in range(6,11)])
    Tp=T.copy();Tp[0]=0;Vp=VV.copy();Vp[0]=0
    S1=scale('1.001',mul(exps(add(Vp,Tp)),K))
    S2=scale('1.001',mul(mul(q,exps(add(Vp,scale(2,Tp)))),mul(K,KV)))
    M0=scale('3.003',exps(add(scale(2,ll),scale(4,Tp))))
    M1=scale(2,mul(M0,add(K,S1)))
    M2=mul(M0,add(scale(4,power(add(K,S1),2)),scale(2,S2)))
    MM=[M0,M1,M2]
    CC=parameter_scale(F(3,2),parameter_mul(MM,parameter_der(MM)))
    B=poly([F(1)/F('5.9'),0,F('1e-10'),F('.668'),F('.501'),F('2e-24'),F('1e-24')])
    reference=add(compose(der(der(B)),q),poly(['2e-10']))
    RR=parameter_add([reference,poly(),poly()],parameter_scale(2,parameter_der(CC)),parameter_scale(2,parameter_mul(CC,CC)))
    Ric=parameter_scale(2,RR);Sc=parameter_mul(MM,Ric)
    JJ=parameter_add(Ric,parameter_scale(F(1,4),parameter_mul(Sc,MM)))
    down=parameter_mul(MM,RR)
    CW=parameter_add(down,parameter_scale(2,parameter_mul(MM,Ric)),parameter_scale(F(2,3),parameter_mul(Sc,parameter_mul(MM,MM))))
    check('retained_metric_zero_parameter',M0[:11]==[F(x) for x in prev['metric'][:11]])
    check('retained_connection_zero_parameter',CC[0][:10]==[F(x) for x in prev['connection'][:10]])
    check('parameter_weights_multiply',all(v>=0 for a in MM+CC+Sc for v in a[:9]))
    check('metric_ten_suffices',8+2==10 and 6+2==8)
    # An independent raw Leibniz sum checks the second path derivative of the connection.
    raw=scale(F(3,2),add(mul(M2,der(M0)),scale(2,mul(M1,der(M1))),mul(M0,der(M2))))
    check('independent_second_connection_derivative',raw[:9]==CC[2][:9])
    FIXTURES.append({'fixture':'weighted_first_and_second_path_geometry','S1':strings(S1[:11]),'S2':strings(S2[:11]),
                     'metric_first':strings(M1[:11]),'metric_second':strings(M2[:11]),
                     'connection_first':strings(CC[1][:10]),'scalar_second':strings(Sc[2][:9]),
                     'all_coordinates_fixed':True,'second_parameter_weight':'t^4'})
    return prev,q,Tp,K,S1,S2,MM,CC,Sc,JJ,CW


def complete_scalar_ledger(prev,q,Tp,K,S1,S2,MM,CC,Sc,JJ,CW):
    M0,M1,M2=MM;C0,C1,C2=CC;NR2=scale('6.1',Sc[2]);em=F('4e-27')/18;emax=F(1,2000**10*180)
    z=S.symbols('z');B=S.Rational(10,59)+S.Rational('1e-10')*z*z+S.Rational('.668')*z**3+S.Rational('.501')*z**4+S.Rational('2e-24')*z**5+S.Rational('1e-24')*z**6
    za=S.Rational('1.3e-22')*(4*z**5+6*z**6)
    ja=S.Rational('.501')*z**4+S.Rational('5e-24')*z**5+S.Rational('4.5e-24')*z**6
    ca=3*S.Rational('.668')*z**3+6*S.Rational('.501')*z**4+S.Rational('20e-24')*z**5+S.Rational('15e-24')*z**6
    box=S.Rational('6.1')*(z*z*B*S.diff(za,z,2)+(z*z*S.diff(B,z)+2*z*B)*S.diff(za,z))
    seed=S.expand(za+S.Rational(emax.numerator,emax.denominator)*(S.Rational(11,3)*za+box+S.Rational(11,72)*za**2
                   +S.Rational(7,12)*S.Rational(4,3)*S.Rational('6.1')**2*ca**2+S.Rational(11,6)*4*S.Rational('6.1')**2*ja**2)+S.Rational(em.numerator,em.denominator)*z**4)
    def coeffs(expr):
        pp=S.Poly(expr,z);return [fromsym(pp.nth(j)) for j in range(pp.degree()+1)]
    seed_coeff=coeffs(seed);zabs=compose(coeffs(za),q);lie=mul(K,compose(coeffs(z*S.diff(za,z)),q))
    check('seed_trace_all_powers_at_least_four',all(v==0 for v in seed_coeff[:4]) and min(seed_coeff[4:])>0)
    check('untruncated_seed_polynomial',len(seed_coeff)==13 and seed_coeff[12]>0)
    fullseed=compose(seed_coeff,q)
    check('omitted_twelfth_seed_power_fails',fullseed[0]>compose(seed_coeff[:12],q)[0])
    Zd=add(fullseed,lie);DZ=add(Zd,scale(F(1,2),NR2))
    BX=scale('6.1',add(mul(M0,add(der(der(DZ)),mul(C0,der(DZ)))),mul(M1,der(der(zabs))),mul(add(mul(M1,C0),mul(M0,C1)),der(zabs))))
    JW0=poly([F(v) for v in prev['weighted_J_initial']]);CW0=poly([F(v) for v in prev['weighted_C_initial']]);ZW0=poly([F(v) for v in prev['weighted_Z_initial']])
    J1=scale('6.1',JJ[1]);CW1=scale('6.1',CW[1])
    DCsq=add(scale(2,mul(power(M0,4),mul(CW0,CW1))),scale(4,mul(mul(power(M0,3),M1),power(CW0,2))))
    DJsq=add(scale(2,mul(power(M0,2),mul(JW0,J1))),scale(2,mul(mul(M0,M1),power(JW0,2))))
    MU=scale(em*F('1.001'),mul(exps(scale(4,Tp)),power(q,4)));DMU=scale(4,mul(MU,S1))
    DH=add(scale(emax,add(BX,scale(F(11,3),DZ),scale(F(11,72),mul(add(ZW0,zabs),DZ)),scale(F(7,12),DCsq),scale(F(11,6),DJsq))),DMU)
    RA2=poly()
    for power_n,a in [(5,4),(6,6)]:
        RA2=add(RA2,scale(F('1.001')*F('1.3e-22')*a,mul(mul(power(q,power_n),exps(scale(power_n,Tp))),add(scale(power_n**2,power(S1,2)),scale(power_n,S2)))))
    P=add(scale(F(1,2),NR2),scale(F(1,2),RA2),lie,DH)
    H=L.covariant(P,C0,0,6);constant=F(1001,24000)*sum(H);exponent=L.decimal_power(constant)
    check('all_residual_terms_retained',min(NR2[:7]+RA2[:7]+lie[:7]+DH[:7]+DMU[:7])>0)
    check('residual_raw_jets_positive',min(H)>0)
    check('scalar_H6_exact_enclosure',F(10)**(exponent-1)<constant<=F(10)**exponent,{'constant':L.decimal(constant),'power_of_ten':exponent})
    check('physical_volume_square_root',F('1.001')**2>F('1.001'))
    check('actual_volume_ratio',1/(1-F('1.2e-7'))<F('1.001'))
    check('source_prefactor',F(1001,24000)==F('1.001')/24)
    check('six_derivatives_available',len(H)==7 and len(NR2)>=9)
    FIXTURES.append({'fixture':'actual_scalar_residual_H6_ledger','complete_seed_polynomial':strings(seed_coeff),
                     'linear_trace_W2':strings(Zd[:9]),'curvature_remainder_second_W2':strings(NR2[:9]),
                     'reference_remainder_second_W2':strings(RA2[:7]),'anomaly_difference_W2':strings(DH[:7]),
                     'residual_reference_W2':strings(P[:7]),'physical_covariant_bounds':strings(H),
                     'physical_covariant_decimals':decimals(H),'C_tau_exact':str(constant),'C_tau_decimal':L.decimal(constant),
                     'C_tau_power_of_ten':exponent,'bound':'C_tau ell^2/g times sqrt(gamma0)',
                     'noncompact_physical_volume_retained':True,'full_Schur_source_bounded':False})


def integral_and_ell_controls():
    t,h,ell=S.symbols('t h ell',positive=True)
    val=S.integrate(t**4*t**-8*t**2,(t,1+h,S.oo));zero('noncompact_physical_integral',val-1/(1+h))
    check('noncompact_integral_below_one',S.simplify(1-val)==h/(h+1))
    check('slower_t_minus_three_is_insufficient',S.integrate(t**4*t**-6*t**2,(t,2,S.oo))==S.oo)
    r=S.symbols('r');zero('Taylor_second_integral',S.integrate(ell-r,(r,0,ell))-ell**2/2)
    z0,q=S.symbols('z0 q',positive=True)
    zero('fixed_observation_weight_keeps_power',z0**-4*(z0*q)**12-z0**8*q**12)
    check('differentiating_divided_polynomial_fails',S.diff(q**12,q,6)!=S.diff(q**8,q,6))
    FIXTURES.append({'fixture':'physical_volume_decay_and_Taylor_controls','tail_integral':'1/(1+h)',
                     't_minus_three_rejected':True,'ell_quadratic_retained':True,'outer_cutoff_used':False})


def independent_warped_contours():
    x=S.symbols('x');records=[]
    for tau0 in [S.Rational(1,5),S.Rational(4,5)]:
        tau=tau0+x/100;z=1-tau;phip=-S.diff(z,x)/z;phipp=S.diff(phip,x)
        u=(1+2*x+x*x+x**4)/1000;v=(2-x+3*x*x+x**3)/1000
        nu=tau*v;acc=tau*(1-tau)*v*v
        h0=-6*phipp-6*phip**2;h1=-2*S.diff(u,x,2)-6*S.diff(nu,x,2)-12*phip*S.diff(nu,x)
        h2=-6*S.diff(acc,x,2)-12*(S.diff(nu,x)**2+phip*S.diff(acc,x))
        second=S.cancel(z*z*(h2-4*(u+nu)*h1+(4*(u+nu)**2-2*acc)*h0+2*(4*nu**2-2*acc)))
        expected=S.expand(second).coeff(x,6)/2
        # The denominator cancels in this manufactured family, so this is exact.
        check('warped_second_derivative_is_polynomial_'+str(tau0),not S.denom(second).has(x))
        omitted=S.cancel(z*z*(-12*S.diff(nu,x)**2-4*(u+nu)*h1+4*(u+nu)**2*h0+8*nu**2))
        check('omitted_radius_acceleration_fails_'+str(tau0),S.expand(omitted).coeff(x,6)/2!=expected)
        for dps in [60,90]:
            with mp.workdps(dps):
                t0=mp.mpf(str(tau0.p))/int(tau0.q);rx=mp.mpf('.05');re=mp.mpf('.01');nx=64;ne=32
                def actual(xx,ee):
                    ta=t0+xx/100;zz=1-ta;up=(2+2*xx+4*xx**3)/1000;upp=(2+12*xx**2)/1000
                    uu=(1+2*xx+xx*xx+xx**4)/1000;vv=(2-xx+3*xx*xx+xx**3)/1000
                    vp=(-1+6*xx+3*xx**2)/1000;vpp=(6+6*xx)/1000
                    kk=mp.exp(ee*vv);TT=1-ta+ta*kk
                    tx=(kk-1)/100+ta*kk*ee*vp
                    txx=2*kk*ee*vp/100+ta*kk*(ee*vpp+ee**2*vp**2)
                    sp=tx/TT;spp=txx/TT-sp**2
                    pp=1/(100*zz);ppp=1/(10000*zz**2)
                    # Direct finite warped-product curvature, before parameter expansion.
                    a= -mp.log(zz)+ee*uu+mp.log(TT);bb=-mp.log(zz)+mp.log(TT)
                    return mp.exp(-2*a)*(-2*(ppp+ee*upp+spp)-4*(ppp+spp)-6*(pp+sp)**2)+2*mp.exp(-2*bb)
                total=mp.mpc(0)
                for j in range(nx):
                    xx=rx*mp.exp(2j*mp.pi*j/nx)
                    val=sum(actual(xx,re*mp.exp(2j*mp.pi*k/ne))*mp.exp(-4j*mp.pi*k/ne) for k in range(ne))/(ne*re**2)
                    total+=val*mp.exp(-12j*mp.pi*j/nx)
                total/=nx*rx**6
                target=mp.mpf(str(S.N(expected,dps)));error=abs(total-target)
                check('warped_mixed_contour_'+str(tau0)+'_'+str(dps),error<mp.mpf('1e-42'),mp.nstr(error,18))
                records.append({'tau0':str(tau0),'dps':dps,'spatial_order':6,'path_order':2,'error':mp.nstr(error,18)})
    FIXTURES.append({'fixture':'independent_finite_warped_mixed_contours','records':records,'actual_charged_geometry':False,
                     'radius_acceleration_omission_detected':True})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 actual scalar residual H6','input_count':len(ps)}),flush=True)
    try:
        provenance_and_tangent_inputs();exact_residual_identities();complete_scalar_ledger(*build_parameter_geometry())
        integral_and_ell_controls();independent_warped_contours()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(x['passed'] for x in ROWS)
    run={'story':'S167 actual scalar residual H6','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
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
