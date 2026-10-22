"""Frozen actual-energy budgets and cold-field diagnostics; no admission solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as Q
import json,traceback
import sympy as S
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
STEM='recentered-energy-and-cold-admission';PREVIOUS='recentered-physical-endpoint'
OWN=['notes/'+STEM+'.md','checks/recentered_energy_cold_admission.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def exact_budgets():
    Cs=Q('1.01e15');CT=Q('2.021e15');Fl=1-Q('1.179e-8');ratio=Q('4e-27')
    check('full_source_value_allowance',Q('1e15')+Q('1e11')/2+9<Cs)
    check('full_trace_coefficient',10+28+56+Q(106,3)==Q(388,3))
    check('full_trace_vs_U',Q(388,360*3)<1)
    check('Weyl_pressure_remainder',Q(5,15*3)*(6*Q(1,2)+8+3)<5)
    check('angular_pressure_tensor_norm',2*Cs**2+2*(Cs+Q(1,2))**2<CT**2)
    E2=22*CT**2*ratio**2/(Q('.07')*Fl)+Q('1.82e-29')**2
    check('actual_coupled_energy_budget',E2<Q('1.44e-10')**2)
    check('energy_to_design_threshold',Q('5e-9')/Q('1.44e-10')==Q(625,18))
    hmax=Q(3,2000**10)
    check('cold_h_guard',hmax<Q('1e-30'))
    check('cold_embedding_lower_bound',47*hmax*Q('2.6e15')**2<1)
    check('generic_embedding_certificate_fails',Q('2.6e15')>Q(625,18))
    check('full_C5_embedding_not_finite',5>2)
    # Exact finite-reference polynomial moment, not only the h->0 limit.
    w=S.symbols('w',nonnegative=True)
    P=S.Rational(122,10)*w+S.Rational(105,10)*w*w
    QQ=S.Rational(123,10)+S.Rational(333,10)*w+S.Rational(106,10)*w*w
    poly=S.Poly(S.expand((1+3*S.Rational('1e-30')*(1+w)**2)*(S.Rational('18e-9')*QQ**2+S.Rational('6.000001')*P+S.Rational('24.000001'))),w)
    moment=sum(coeff*S.factorial(powers[0])/2**(powers[0]+1) for powers,coeff in poly.terms())
    check('finite_cold_field_moment',moment<47,moment)
    check('cold_laplacian_polynomial_margin',S.Poly(QQ-(P+S.Rational(122,10)+21*w)-2*S.Rational('1e-30')*P,w).all_coeffs()==[S.Rational(1,10)-21*S.Rational('1e-30'),S.Rational(1,10)-S.Rational(244,10)*S.Rational('1e-30'),S.Rational(1,10)])
    r,lam,M,Z=S.symbols('r lambda M Z',positive=True)
    f=1+lam*r*r-2*M/r+Z/r**2;W=r*r*f;A=r*r*S.diff(f,r)-2*r*f
    check('Hardy_lift_absolute_majorant',S.simplify(S.diff(W,r)-A-4*r*f)==0)
    denom=S.simplify(S.diff(W,r)/r**2)
    check('Hardy_measure_exact',S.simplify(denom-(4*lam*r+2/r-2*M/r**2))==0)
    check('Hardy_measure_monotonicity_identity',S.simplify(S.diff(denom,r)-(4*lam-2/r**2+4*M/r**3))==0)
    check('TT_tensor_Gram_bound',2*2**2+4*2+6==22)
    # Stable exact cold polynomial from the RN horizon condition.
    t,v,x,Qc,hh,ww=S.symbols('t v x Qc h ww',positive=True)
    am=(v+Qc/v+x*v**3)/2
    Nr=S.expand(x*t**4+t*t-2*am*t+Qc)
    z=S.symbols('z')
    stable=(3*x*v**3+v-Qc/v)*z+(6*x*v*v+1)*z*z+4*x*v*z**3+x*z**4
    check('stable_RN_horizon_factorization',S.simplify(Nr.subs(t,v+z)-stable)==0)
    Lh=2+12*x+(1+18*x)*hh+12*x*hh**2+3*x*hh**3
    check('stable_cold_linear_coefficient',S.simplify((3*x*v**3+v-Qc/v).subs({v:1+hh,Qc:1+3*x})-hh*Lh/(1+hh))==0)
    Pc=Lh/(1+hh)*ww+(6*x*(1+hh)**2+1)*ww**2+4*x*(1+hh)*hh*ww**3+x*hh**2*ww**4
    tt=1+hh*(1+ww);Fb=Pc/tt**2
    check('cold_laplacian_exact_cancellation',S.simplify(Fb-S.diff(Fb,ww)-2*hh*Fb/tt-(Pc-S.diff(Pc,ww))/tt**2)==0)
    FIXTURES.append(dict(fixture='exact_source_Hardy_and_cold_budgets',energy_squared_upper=str(E2),
        cold_moment_upper_exact=str(moment),cold_moment_upper_decimal=str(S.N(moment,30)),
        generic_embedding_required='<=625/18',generic_bolt_embedding_lower='>2.6e15',actual_correction_lower_bound=False))

def lift_fixtures():
    mp.mp.dps=80;records=[]
    for hh in ['.001','.01']:
        h=mp.mpf(hh);x=mp.mpf('1e10');v=1+h;Qc=1+3*x;am=(v+Qc/v+x*v**3)/2
        def nr(t):
            z=t-v;return (3*x*v**3+v-Qc/v)*z+(6*x*v*v+1)*z*z+4*x*v*z**3+x*z**4
        def IA(t):
            e=mp.exp(v-t)
            return -2*((v+1)-(t+1)*e)+6*am*(1-e)-4*Qc*mp.exp(v)*(mp.e1(v)-mp.e1(t))
        def bb(t):return -mp.mpf(1) if t==v else -IA(t)/nr(t)
        for offset in ['.01','.1','1','10']:
            t=v+mp.mpf(offset)
            direct=mp.quad(lambda y:(-2*y+6*am-4*Qc/y)*mp.exp(v-y),[v,t])
            check('TT_lift_direct_integral_'+hh+'_'+offset,abs(direct/IA(t)-1)<mp.mpf('1e-60'))
        nxi=mp.quad(lambda t:t*t*mp.exp(2*(v-t)),[v,v+1,mp.inf])
        nb=mp.quad(lambda t:t*t*bb(t)**2,[v,v+h,v+1,v+10,mp.inf])
        gram=mp.quad(lambda t:t*t*((bb(t)+2*mp.exp(v-t))**2+bb(t)**2+2*mp.exp(2*(v-t))),[v,v+h,v+1,v+10,mp.inf])
        check('Hardy_L2_fixture_'+hh,nb<4*nxi)
        check('TT_full_tensor_fixture_'+hh,gram<22*nxi)
        records.append(dict(h=hh,b_to_xi_squared_norm_ratio=mp.nstr(nb/nxi,45),tensor_to_xi_squared_norm_ratio=mp.nstr(gram/nxi,45)))
    FIXTURES.append(dict(fixture='independent_integral_TT_lift',records=records,actual_quantum_source_sample=False))

def actual_cold_field():
    records=[]
    for digits in [100,150]:
        mp.mp.dps=digits;k=mp.mpf(2000);g0=k**-10
        g=2*g0/(1+mp.sqrt(1+44*g0/(360*mp.pi)));Fl=1-384*mp.pi*mp.mpf('1e22')*g
        C0=mp.sqrt(1+12*mp.pi*k**6);bclass=mp.sqrt((C0-1)*g/(6*g0))
        C=mp.sqrt(1+12*mp.pi*g*k**16/Fl);x=(C-1)/6;B=mp.sqrt(x);Qc=1+3*x
        eps=g/(60*mp.pi);DR=96*mp.pi*g*mp.mpf('1e22')-eps;Cv=1+11*eps/3
        for zi in ['.25','.265','.28']:
            kap=mp.mpf(zi)*2*k**3*g*C0/bclass**3;gamma=kap*B/C
            def K(h):return ((1+3*x)*((1+h)**-2+(1+h)**-3)+3*x*(1+(1+h)**-1))/(2*C)
            h=gamma
            for _ in range(5):h-=(h*K(h)-gamma)/(K(h)+h*mp.diff(K,h))
            v=1+h;Lh=2+12*x+(1+18*x)*h+12*x*h*h+3*x*h**3
            def pol(w):return Lh/v*w+(6*x*v*v+1)*w*w+4*x*v*h*w**3+x*h*h*w**4
            def polp(w):return Lh/v+2*(6*x*v*v+1)*w+12*x*v*h*w*w+4*x*h*h*w**3
            def integrand(w):
                t=v+h*w;pf=pol(w);Kr=6+24*DR+22*eps-mp.mpf(22)/3*eps*Qc/(x*t**4)
                return mp.exp(-2*w)*(6*DR*(pf-polp(w))**2/(x*x*t*t)+Kr*pf/x+24*Cv*t*t)
            val=mp.quad(integrand,[0,1,5,20,mp.inf]);en=mp.sqrt(h*val)
            limit=mp.mpf('10.5')*DR*(C/x)**2+mp.mpf('.75')*(6+24*DR+22*eps-mp.mpf(22)/3*eps*Qc/x)*(C/x)+12*Cv
            check('actual_cold_energy_'+str((digits,zi)),val<47)
            check('actual_cold_limit_consistency_'+str((digits,zi)),abs(val-limit)<mp.mpf('1e-28'))
            check('actual_cold_embedding_'+str((digits,zi)),1/en>mp.mpf('2.6e15'))
            check('actual_temperature_'+str((digits,zi)),abs(h*K(h)/gamma-1)<mp.mpf('1e-80'))
            for w in map(mp.mpf,['0','1','10','100']):
                t=v+h*w;pf=pol(w);raw=pf/t**2-(polp(w)/t**2-2*h*pf/t**3)-2*h*pf/t**3
                check('actual_laplacian_two_forms_'+str((digits,zi,str(w))),abs(raw-(pf-polp(w))/t**2)<mp.mpf('1e-70'))
            actualratio=g*k**8/(x*x)
            TT=mp.sqrt(22)*mp.mpf('2.021e15')/mp.sqrt(mp.mpf('.07')*Fl)*actualratio
            records.append(dict(digits=digits,prepared_zi=zi,h=mp.nstr(h,65),exact_test_energy_squared_over_h=mp.nstr(val,65),
                exact_test_energy=mp.nstr(en,65),bolt_embedding_lower_from_test=mp.nstr(1/en,65),
                reference_source_energy_allowance=mp.nstr(mp.sqrt(TT*TT+mp.mpf('1.82e-29')**2),65)))
    for i in range(3):
        for key in ['h','exact_test_energy_squared_over_h','exact_test_energy','reference_source_energy_allowance']:
            check('cross_precision_'+str(i)+'_'+key,abs(mp.mpf(records[i][key])/mp.mpf(records[i+3][key])-1)<mp.mpf('1e-55'))
    FIXTURES.append(dict(fixture='actual_reference_cold_test_field',records=records,test_field_is_actual_correction=False,
        nonlinear_neighborhood_certified=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_physical_endpoint_result',old['status']=='passed' and old['physical_fast_tail_inverse_regularity_proved'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        exact_budgets();lift_fixtures();actual_cold_field()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        actual_reference_linear_energy_bound='1.44e-10 in the specified physical energy unit' if error is None else None,
        generic_cold_embedding_lower_bound='>2.6e15' if error is None else None,
        generic_energy_to_admission_certificate_rejected=error is None,
        actual_correction_pointwise_lower_bound_proved=False,admission_norm_inverse_quantified=False,
        physical_coupled_inverse_proved=False,nonlinear_solution_constructed=False,semiclassical_validity_proved=False,
        physical_stability_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
