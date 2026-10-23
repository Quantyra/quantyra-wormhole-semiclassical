"""Frozen auxiliary equivalence and conditional scalar/metric bound.

No metric inverse or nonlinear boundary-domain transfer is certified.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as Q
import json,traceback
import sympy as S
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1]
STEM='anomaly-complete-auxiliary-response';PREVIOUS='product-centered-scalar-cold-response'
OWN=['notes/'+STEM+'.md','checks/anomaly_complete_auxiliary_response.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def algebra():
    G,c,lam,l0,R,phi,C,H,Lphi,LC=S.symbols('G c lambda lambda0 R phi C H lap_phi lap_C')
    eps=G/(60*S.pi);dc=1/(5760*S.pi**2);ce=c-dc;D=96*S.pi*G*ce
    check('exact_effective_R2_coefficient',S.simplify(D-(96*S.pi*G*c-eps))==0)
    check('local_trace_derivative_cancellation',S.simplify(-12*dc+1/(480*S.pi**2))==0)
    check('total_action_repartition',S.expand(ce*R**2+dc*R**2-c*R**2)==0)
    check('total_boundary_repartition',S.simplify(4*ce*R+4*dc*R-4*c*R)==0)
    Phi=1+32*S.pi*G*ce*phi
    minus_trace=Phi*R+12*l0-32*S.pi*G*ce*phi**2-D*Lphi+eps*(H-S.Rational(11,72)*R**2)
    scalar=phi+12*l0-D*Lphi+eps*(H-S.Rational(11,72)*phi**2)
    bracket=Phi-S.Rational(11,36)*eps*phi-S.Rational(11,72)*eps*C
    check('exact_algebraic_constraint_recovery',S.simplify((minus_trace-scalar).subs(R,phi+C)-C*bracket)==0)
    oldPhi=1+32*S.pi*G*c*phi
    oldtrace=oldPhi*R+12*l0-32*S.pi*G*c*phi**2-96*S.pi*G*c*Lphi+eps*(H-S.Rational(11,72)*R**2+Lphi+LC)
    oldbracket=oldPhi-S.Rational(11,36)*eps*phi-S.Rational(11,72)*eps*C
    check('unrepartitioned_constraint_has_derivative',S.simplify((oldtrace-scalar).subs(R,phi+C)-C*oldbracket-eps*LC)==0)
    density,ddensity,dC=S.symbols('density delta_density delta_C')
    check('fixed_density_constraint_contact',S.simplify(S.diff(density*C*bracket,C).subs(C,0)-density*(Phi-S.Rational(11,36)*eps*phi))==0)
    u=S.symbols('u');Fs=1-384*S.pi*G*ce*lam;F0=1-384*S.pi*G*ce*l0
    K=1/(1024*S.pi**2*G**2*ce);pv=-3*(l0-lam)/(8*S.pi*G)
    check('vacuum_potential_identity',S.simplify(pv/K-(F0-Fs))==0)
    U=K*(Fs*Fs-2*Fs*S.exp(-2*u)+F0*S.exp(-4*u))-pv*S.exp(-4*u)
    expected=-3*Fs*lam/(8*S.pi*G)+K*Fs*(1-S.exp(-2*u))**2
    check('matched_frame_potential',S.simplify(U-expected)==0)
    check('local_mixed_potential_zero',S.simplify(S.diff(expected,u).subs(u,0))==0)
    phiu=(Fs*S.exp(2*u)-1)/(32*S.pi*G*ce)
    # Transformation of the scalar Laplacian: phi''=2phi' cancels
    # the gradient term from Delta_g=e^(2u)(Delta_E-2du.grad_E).
    check('auxiliary_laplacian_gradient_cancellation',S.simplify(S.diff(phiu,u,2)-2*S.diff(phiu,u))==0)
    check('scalar_row_rescale',S.simplify(16*S.pi*G*ce/(D*Fs)-1/(6*Fs))==0)
    Pot=S.exp(-4*u)*(phiu+12*l0-S.Rational(11,72)*eps*phiu**2)/(6*Fs)
    vac={l0:lam+S.Rational(11,6)*eps*lam**2}
    check('full_scalar_vacuum_zero',S.simplify(Pot.subs(vac).subs(u,0))==0)
    check('full_scalar_vacuum_mass',S.simplify(S.diff(Pot.subs(vac),u).subs(u,0)-(1+S.Rational(11,3)*eps*lam)/D)==0)
    check('frame_F_difference',S.simplify(Fs-(1-384*S.pi*G*c*lam)-G*lam/(15*S.pi))==0)
    # Independent arbitrary diagonal curvature/Hessian/gradient tensor
    # trace subtraction in four dimensions.
    r0,r1,r2,r3,k0,k1,k2,k3,z=S.symbols('r0 r1 r2 r3 k0 k1 k2 k3 z')
    Ric=S.diag(r0,r1,r2,r3);HH=S.diag(k0,k1,k2,k3);I=S.eye(4);gg=S.diag(z*z,0,0,0)
    transformed=Ric+2*HH+S.trace(HH)*I+2*gg-2*z*z*I
    Jg=transformed-S.trace(transformed)*I/4
    expectedJ=Ric-S.trace(Ric)*I/4+2*(HH-S.trace(HH)*I/4)+2*(gg-z*z*I/4)
    check('conformal_traceless_Ricci_identity',S.simplify(Jg-expectedJ)==S.zeros(4))
    FIXTURES.append(dict(fixture='action_trace_constraint_and_frame_algebra',derivative_constraint_removed_by_exact_repartition=True,
        fixed_density_contacts_retained=True,total_theory_unchanged=True))

def rational_budgets():
    pi=Q('3.14');piup=Q(22,7);g=Q(1,2000**10);c=Q('1e22');F=1-Q('1.179e-8');x=Q('8e9');d=Q('8.6e-19')
    epslam=g/(60*pi)
    check('local_split_relative_size',1/(5760*pi*pi*c)<Q('1.8e-27'))
    check('frame_coefficient_shift',g/(15*pi)<Q('2.1e-35'))
    constraint=1-32*piup*c*g*13-Q(11,36)*epslam*13-Q(11,72)*epslam
    check('constraint_bracket_positive',constraint>1-Q('2e-8'))
    check('constraint_recovery_inverse',1/(1-Q('2e-8'))<Q('1.00000003'))
    zmax=32*piup*c*g/F*10*d/x
    ub=zmax/(2*(1-zmax));rdu=16*piup*c*g/F*56*d/x/(1-zmax)
    r2du=16*piup*c*g/F*372*d/x/(1-zmax)+2*rdu*rdu
    check('seed_auxiliary_u_value',ub<Q('5.4e-37'))
    check('seed_auxiliary_gradient',7*rdu*rdu/(1-zmax)**2<Q('1e-35')**2)
    hrr=7*r2du+7*rdu+7*rdu*rdu;hoo=7*rdu+7*rdu*rdu
    check('seed_auxiliary_Hessian',(hrr*hrr+3*hoo*hoo)/(1-zmax)**4<Q('1e-33')**2)
    check('seed_Riemann_guard',4*(15**2+4*Q('6.1')**2+Q('8.1')**2)<50**2)
    check('barrier_Hessian_guard',Q('8.4')**2+3*Q('2.8')**2<10**2)
    principal=Q(11,9)*epslam/F*20;drift=Q(22,9)*epslam/F*20*Q('1e-35')
    check('scalar_principal_margin',principal<Q('1e-32'))
    check('scalar_drift_margin',drift<Q('1e-66'))
    check('mass_lower_margin',1/Q('3e-9')>Q('3.3e8'))
    check('scalar_barrier_margin',Q('3.3e8')-Q('3.38')-20*Q('1e-32')-2*Q('1e-66')>Q('3.2e8'))
    # Variation allowances include lower-index tensor contractions.
    dC2=2*100*5000+4*100**2;dZ2=2*20*6000+4*20**2
    check('Weyl_square_variation',dC2==1040000)
    check('traceless_square_variation',dZ2==241600)
    curv=epslam/(6*F)*(Q(7,12)*dC2+Q(11,6)*dZ2)
    check('curvature_mixed_source',curv<Q('2e-30'))
    mag=40/(360*pi*F)*Q('4e-27')
    check('magnetic_mixed_source',mag<Q('1.5e-28'))
    check('fixed_density_source_allowance',Q('9e-28')/(6*F*(1-zmax)**2)<Q('1.51e-28'))
    check('complete_mixed_row_budget',Q('2e-33')+Q('2e-30')+Q('1.5e-28')+Q('1.51e-28')<Q('1e-27'))
    hmin=Q(1,8*2000**10);boost5=4/(64**5*hmin**3)
    mixed_response=(Q('1.51e-28')+Q('1e-27'))/Q('3.2e8')
    source_response=Q('1.51e-28')/Q('3.2e8')
    check('conditional_auxiliary_cold_value',(mixed_response/Q('2e-17'))**5*boost5<1)
    check('source_only_auxiliary_cold_value',(source_response/Q('2.2e-18'))**5*boost5<1)
    check('Jordan_metric_return_HS_factor',4*2**2==4**2)
    FIXTURES.append(dict(fixture='rational_constraint_frame_and_response_budgets',constraint_bracket_lower=str(constraint),
        mixed_curvature_coefficient=str(curv),mixed_magnetic_coefficient=str(mag),cold_value_conditional='2e-17 if N2(h_E)<=1',
        metric_input_N2_certified=False,full_metric_response_certified=False))

def geometric_identities():
    t,x,d,delta=S.symbols('t x d delta',positive=True)
    fext=(t-1)**2/t**2*(1+x*(t*t+2*t+3)+d/t**2)
    fp_expected=2*x*(t-1)*(1+1/t+1/t**2+3/t**3)+2*(t-1)/t**3-2*d*(t-1)*(t-2)/t**5
    check('global_seed_fprime_factorization',S.simplify(S.diff(fext,t)-fp_expected)==0)
    check('global_f_upper_bracket',1+6*Q('8e9')+Q('1e-12')<7*Q('8e9'))
    check('global_fprime_bound',12+Q(2,8*10**9)+6*Q('1e-12')/(8*10**9)<13)
    check('thermal_fprime_reserve',7*Q('1e-30')<1)
    p=S.Rational(2,5);a=S.symbols('a',positive=True)
    logprime=p*a/(t*(t-a))
    second=S.diff(logprime,t)+logprime**2
    check('weight_second_derivative',S.simplify(second-p*a*((p+1)*a-2*t)/(t*t*(t-a)**2))==0)
    # Physical curvature and gradient remain explicit rather than being
    # absorbed into a redefinition of an input norm.
    check('weight_gradient_guard',7*Q(2,5)**2<4)
    FIXTURES.append(dict(fixture='global_seed_and_weight_geometry',actual_thermal_mass_retained=True,
        actual_quantum_metric_kernel_sampled=False))

def benchmark():
    records=[]
    old=json.loads((ROOT/'evidence/recentered-product-and-cold-residual.json').read_bytes())['runs'][-1]
    sphere_shift=old['fixtures'][3]['records'][0]['sphere_squared_first_shift']
    for digits in [90,140]:
        mp.mp.dps=digits;k=mp.mpf(2000);q=k**4;c=mp.mpf('1e22');dc=1/(5760*mp.pi**2);ce=c-dc
        g0=k**-10;g=2*g0/(1+mp.sqrt(1+44*g0/(360*mp.pi)))
        F=1-384*mp.pi*c*g;Fs=1-384*mp.pi*ce*g
        C0=mp.sqrt(1+12*mp.pi*k**6);x0=(C0-1)*g/g0/6
        Crn=mp.sqrt(1+12*mp.pi*g*k**16/F);x=(Crn-1)/6*(1+mp.mpf(sphere_shift))
        sig=g/(60*mp.pi*x);eps=g/(60*mp.pi);D=96*mp.pi*ce*g
        aa=sig;bb=sig*(12*x+mp.mpf(11)/3)-1;cc=sig*(47*x*x+22*x+mp.mpf(11)/3-5*q*q)
        d=2*cc/(-bb+mp.sqrt(bb*bb-4*aa*cc));C=1+6*x+d
        check('benchmark_frame_split_'+str(digits),abs((Fs-F)/(g/(15*mp.pi))-1)<mp.mpf('1e-50'))
        l0=1+mp.mpf(11)/6*eps
        def ph(u):return (Fs*mp.exp(2*u)-1)/(32*mp.pi*ce*g)
        def pot(u):return mp.exp(-4*u)*(ph(u)+12*l0-mp.mpf(11)/72*eps*ph(u)**2)/(6*Fs)
        def ub(t):return mp.log1p(32*mp.pi*ce*g/Fs*d/x*(4/t**5-6/t**6))/2
        check('benchmark_scalar_mass_'+str(digits),abs(mp.diff(pot,0)/((1+mp.mpf(11)/3*eps)/D)-1)<mp.mpf('1e-65'))
        for zi in ['.25','.265','.28']:
            kap=mp.mpf(zi)*2*k**3*g*C0/x0**mp.mpf('1.5');gamma=kap*mp.sqrt(x)/C
            def K(h):
                v=1+h
                return ((1+3*x)*(v**-2+v**-3)+3*x*(1+v**-1)-d*(v**-4-3*v**-5))/(2*C)
            h=gamma
            for _ in range(5):h-=(h*K(h)-gamma)/(K(h)+h*mp.diff(K,h))
            V=1+h;delta=h*h/V*(1+x*(V*V+2*V+3)+d/V**2)
            def fun(t):return (t-1)**2/t**2*(1+x*(t*t+2*t+3)+d/t**2)-delta/t
            samples=[]
            for tt in [V,1+mp.mpf('1e-8'),mp.mpf('1.001'),mp.mpf(2),mp.mpf(10)]:
                u=ub(tt);up=mp.diff(ub,tt);upp=mp.diff(ub,tt,2)
                ff=mp.mpf(0) if tt==V else fun(tt);fp=mp.diff(fun,tt)
                grad2=ff*up*up/x
                Hess=[(ff*upp+fp*up/2)/x-grad2,fp*up/(2*x)+grad2,ff*up/(x*tt)+grad2,ff*up/(x*tt)+grad2]
                normH=mp.exp(-2*u)*mp.sqrt(sum(z*z for z in Hess))
                check('actual_seed_u_guards_'+str((digits,zi,str(tt))),abs(u)<mp.mpf('5.4e-37') and mp.exp(-u)*mp.sqrt(grad2)<mp.mpf('1e-35') and normH<mp.mpf('1e-33'))
                check('actual_seed_mass_guard_'+str((digits,zi,str(tt))),mp.diff(pot,u)>mp.mpf('3.3e8'))
                samples.append(dict(t=mp.nstr(tt,50),u=mp.nstr(u,50),normalized_gradient=mp.nstr(mp.exp(-u)*mp.sqrt(grad2),50),normalized_Hessian=mp.nstr(normH,50)))
            boost=(2*h/V)**mp.mpf('.4')/(h*(1+mp.log(mp.mpf('.01')/h)))
            cold=(mp.mpf('1.51e-28')+mp.mpf('1e-27'))/mp.mpf('3.2e8')*boost
            check('actual_conditional_cold_value_'+str((digits,zi)),cold<mp.mpf('2e-17'))
            records.append(dict(digits=digits,prepared_zi=zi,relative_coefficient_split=mp.nstr(dc/c,60),F_difference=mp.nstr(Fs-F,60),
                scalar_mass_over_lambda=mp.nstr((1+mp.mpf(11)/3*eps)/D,60),h=mp.nstr(h,60),
                conditional_cold_value_upper=mp.nstr(cold,60),samples=samples))
    for i in range(3):
        for key in ['relative_coefficient_split','F_difference','scalar_mass_over_lambda','h','conditional_cold_value_upper']:
            check('cross_precision_'+str(i)+'_'+key,abs(mp.mpf(records[i][key])/mp.mpf(records[i+3][key])-1)<mp.mpf('1e-50'))
    FIXTURES.append(dict(fixture='benchmark_auxiliary_and_conditional_bounds',records=records,actual_metric_response_solved=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_direct_scalar_result',old['status']=='passed' and old['product_centered_scalar_row_weighted_inverse_proved'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        algebra();rational_budgets();geometric_identities();benchmark()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,
        unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        anomaly_complete_auxiliary_row_equivalence_proved=error is None,auxiliary_curvature_constraint_recovery='algebraic; bracket >1-2e-8' if error is None else None,
        auxiliary_scalar_metric_mixed_bound='1e-27 N2 in normalized weighted source norm' if error is None else None,
        auxiliary_cold_value_conditional_bound='2e-17 if N2(h_E)<=1' if error is None else None,
        Einstein_frame_metric_response_bounded=False,actual_coupled_cold_response_bounded=False,
        full_transformed_boundary_domain_quantified=False,admission_norm_inverse_quantified=False,
        nonlinear_solution_constructed=False,semiclassical_validity_proved=False,physical_stability_proved=False,
        component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
