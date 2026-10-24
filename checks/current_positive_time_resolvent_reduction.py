"""Complete product contact and conditional actual-domain resolvent targets.

This is NOT an execution of the complete signed finite-geometry projection.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import factorial
import json, sys, traceback
sys.path.insert(0, str(Path(__file__).resolve().parent))
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb, ctx
import sympy as sp
import current_disk_heat_propagation as previous
gy=previous.gy; bh=previous.bh; encoded=previous.encoded
root=Path(__file__).resolve().parents[1]
STEM='current-positive-time-resolvent-reduction'; PREVIOUS=previous.STEM
OWN=['notes/'+STEM+'.md','checks/current_positive_time_resolvent_reduction.py',
     'evidence/'+STEM+'-review.json']
ROWS=[]; NUMERICS=[]; N=2048

def check(name, ok, detail=''):
    assert name not in {r['name'] for r in ROWS}, name
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok: raise AssertionError(name+': '+str(detail))

gy.check=check; bh.prior.scalar.prev.prior.check=check

def digest(p): return sha256((root/p).read_bytes()).hexdigest()

def hull(lo,hi):
    lo=arb(lo.lower()); hi=arb(hi.upper())
    return (lo+hi)/2+arb(0,(hi-lo)/2)

def algebra():
    s,x,m=sp.symbols('s x m',positive=True)
    A=sp.Matrix([[2+s,1],[1,3+2*s]])
    P=sp.Matrix([[sp.Rational(1,4)+s/10,sp.Rational(1,7)],
                 [sp.Rational(1,9),sp.Rational(1,5)-s/20]])
    R=A.inv(); F=sp.eye(2)-A*P; Q=P+P.T-P.T*A*P
    simplify=lambda X:X.applyfunc(sp.factor)
    check('noncommuting_two_sided_identity',simplify(R-Q-F.T*R*F)==sp.zeros(2))
    dF=F.diff(s)
    dE=dF.T*R*F+F.T*R*dF-F.T*R*A.diff(s)*R*F
    check('complete_metric_error_derivative',simplify((R-Q).diff(s)-dE)==sp.zeros(2))
    P1=P*(sp.eye(2)+F); F1=sp.eye(2)-A*P1
    check('residual_squares_without_commutation',simplify(F1-F*F)==sp.zeros(2))
    check('both_metric_residual_terms',simplify(F1.diff(s)-dF*F-F*dF)==sp.zeros(2))
    for J in [7,8]:
        series=sum((-x)**j/m**(j+1) for j in range(J+1))
        rem=(-x)**(J+1)/(m**(J+1)*(m+x))
        check('global_Stieltjes_identity_'+str(J),sp.factor(1/(m+x)-series-rem)==0)
    t,lam,a,b=sp.symbols('t lam a b',positive=True)
    S=a*t**2+b*t**3
    full=sp.integrate(sp.exp(-lam*t)*(-2*sp.diff(S,t)),(t,0,sp.oo))
    transformed=-2*lam*sp.integrate(sp.exp(-lam*t)*S,(t,0,sp.oo))
    check('integration_by_parts_full_sign',sp.simplify(full-transformed)==0)
    t0=sp.symbols('t0',positive=True)
    split=sp.integrate(sp.exp(-lam*t)*(-2*sp.diff(S,t)),(t,t0,sp.oo))
    transformed_split=2*sp.exp(-lam*t0)*S.subs(t,t0)-2*lam*sp.integrate(sp.exp(-lam*t)*S,(t,t0,sp.oo))
    check('positive_time_boundary_term',sp.simplify(split-transformed_split)==0)
    z=sp.symbols('z',positive=True); h=sp.Function('h')(z); v=sp.Function('v')(z); Qz=sp.Function('Q')(z)
    check('smeared_a1_contact_flux',sp.expand(h*sp.diff(Qz*sp.diff(v,z),z)+Qz*sp.diff(h,z)*sp.diff(v,z)-sp.diff(h*Qz*sp.diff(v,z),z))==0)
    q,p=sp.symbols('q p',positive=True)
    numerator=2*x*(x+q)-p*(q+2*x)**2
    check('all_mass_moment_monotonicity',sp.expand(numerator+((4*p-2)*x*(x+q)+p*q*q))==0)
    NUMERICS.append(dict(kind='ordinary_argument_scope',
        algebra_examples_are_not_physical_modes=True,selected_domain_argument='notes/current-positive-time-resolvent-reduction.md sections 1, 2 and 4',
        actual_operator_columns_bounded=False,rough_family_C2_action_proved=False))

def moment(p,C,q,n=N):
    partial=sum((2*(q+2*k)*(C/(arb(k)*(k+q)))**p for k in range(1,n+1)),arb(0))
    tail=lambda a:2*C**p/((p-1)*(arb(a)*(a+q))**(p-1))
    lo=partial+tail(n+1); hi=partial+tail(n)
    return hull(lo,hi),dict(power=p,terms=n,partial=encoded(partial),
                           tail_lower=encoded(tail(n+1)),tail_upper=encoded(tail(n)),complete=encoded(hull(lo,hi)))

def heat_fraction(p,C,q,S,n=N):
    t0=arb(1)/1000
    gamma=lambda x:(-x).exp()*sum((x**j/factorial(j) for j in range(p)),arb(0))
    partial=arb(0)
    for k in range(1,n+1):
        lam=arb(k)*(k+q)/C
        partial+=2*(q+2*k)*gamma(lam*t0)/lam**p
    tail=2*C**p/((p-1)*(arb(n)*(n+q))**(p-1))*gamma(arb(n+1)*(n+1+q)/C*t0)
    return hull(partial,partial+tail)/S

def actual_operator_certificate(M1,K,kappa):
    old=json.loads((root/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    row=lambda kind:next(x for x in old['numerics'] if x['kind']==kind)
    uncut=row('full_uncut_transport_norms')['tables']
    cap=row('conditional_complete_cap_kernel_table')['table']
    extra=row('complete_commutator_and_thermal_tails')
    cut=row('complete_actual_domain_propagation_budget')['cutoff_jet_norm']
    de=previous.decode; thermal=de(extra['full_thermal_image_guard']); t0=arb(1)/1000
    actual={}
    for jet in ['base','metric']:
        actual[jet]={}
        for norm in ['inf','row']:
            f=de(uncut['F'][jet][norm])+de(cut[norm])+de(extra['commutator'][norm])+thermal
            pi=de(uncut['P'][jet][norm])+de(cut[norm])*t0**3
            pc=de(cap[jet][norm])+(thermal if norm=='inf' else 0)
            pp=arb(max((2*pc).upper(),pi.upper())) if norm=='inf' else pc+pi
            actual[jet][norm]={'P':pp,'F':f}
    guards={'base':{'inf':{'P':'.19','F':'2.3'},'row':{'P':'3','F':'34'}},
            'metric':{'inf':{'P':'.8','F':'27'},'row':{'P':'11','F':'450'}}}
    for jet in guards:
        for norm in guards[jet]:
            for kernel,guard in guards[jet][norm].items():
                check('complete_kernel_guard_'+jet+'_'+norm+'_'+kernel,actual[jet][norm][kernel]<arb(guard),actual[jet][norm][kernel])
    x=34*t0**4/4; dx=450*t0**4/4
    e0=3*x/(1-x); e1=11*x/(1-x)+3*dx/(1-x)**2
    r0=(-M1*M1*t0).exp()*(1+e0)+x
    s0=kappa*(3+arb.pi()/2)/arb(1).exp()+e1+dx
    check('actual_resolvent_operator_residual',r0<arb('.723'),r0)
    check('actual_resolvent_metric_operator_residual',s0<4,s0)
    check('mass_profile_endpoint_maximum_in_range',M1*M1*t0<2)
    tilt=4*arb(-2).exp()/(M1*M1*t0)**2
    a0=(2*arb('.19')*3/t0).sqrt()*tilt+(2*arb('2.3')*34).sqrt()*arb('3.5').gamma()/M1**7
    b0=(2*arb('.8')*11/t0).sqrt()*tilt+(2*27*450)**arb('.5')*arb('3.5').gamma()/M1**7
    check('actual_resolvent_complete_column_profile',a0<173,a0)
    check('actual_resolvent_metric_complete_column_profile',b0<680,b0)
    records=[]
    for m in [16,32,64]:
        a=173*arb('.723')**(m-1)
        b=680*arb('.723')**(m-1)+(m-1)*4*173*arb('.723')**(m-2)
        error=K*(2*a*b+kappa*a*a)
        passes=error<arb('5e-43')
        check('actual_residual_budget_classification_m'+str(m),bool(passes)==(m==64),error)
        records.append(dict(m=m,column_profile=encoded(a),metric_column_profile=encoded(b),
                            complete_error_allowance=encoded(error),passes_reserve=bool(passes)))
    NUMERICS.append(dict(kind='actual_complete_resolvent_certificate',
        reconstructed_kernel_tables={j:{n:{k:encoded(v) for k,v in a.items()} for n,a in d.items()} for j,d in actual.items()},
        guards=guards,operator_heat_difference=encoded(e0),metric_operator_heat_difference=encoded(e1),
        actual_initial_operator_residual=encoded(r0),actual_initial_metric_operator_residual=encoded(s0),
        actual_initial_column_profile=encoded(a0),actual_initial_metric_column_profile=encoded(b0),
        rounded_actual_guards={'a':'173','b':'680','r':'.723','s':'4'},
        residual_polynomial_tests=records,selected_m=64,selected_residual_squarings=6,
        actual_operator_columns_bounded=True,positive_time_extrapolation_used=False,
        finite_physical_mode_surrogate_used=False,projection_trace_evaluated=False))
    print(json.dumps(dict(event='actual_complete_resolvent_error',a0=str(a0),b0=str(b0),r0=str(r0),s0=str(s0),
        error_by_m={str(r['m']):str(previous.decode(r['complete_error_allowance'])) for r in records},
        actual_domain_error_certified=True,projection_trace_evaluated=False)),flush=True)

def execute():
    p,h,zh,env=gy.actual_parameters(); NUMERICS.extend(gy.NUMERICS)
    B2=p['B2']; C=p['C']; q=arb(p['q']); Fl=1-384*arb.pi()*arb('1e22')*p['g']
    lam1=(q+1)/C; M1=lam1.sqrt(); kappa=2+arb('.25')/M1
    Q=[p['beta'],arb(0),p['u'],-4*p['beta']-2*p['u']-p['tau'],3*p['beta']+p['u']+p['eps'],-2*p['eps'],p['eps']]
    # |grad(-z^4)|^2=16 z^6 Q. Full Bernstein covers preserve actual tau.
    grad2=[arb(0)]*6+[16*x for x in Q]
    from math import comb
    max_grad2=arb(0)
    for cell in range(128):
        lo=arb(cell)/128; width=arb(1)/128; degree=len(grad2)-1
        shifted=[sum((grad2[j]*comb(j,k)*lo**(j-k)*width**k for j in range(k,degree+1)),arb(0)) for k in range(degree+1)]
        bern=[sum((shifted[k]*arb(comb(i,k))/comb(degree,k) for k in range(i+1)),arb(0)) for i in range(degree+1)]
        max_grad2=arb(max(max_grad2.upper(),*(x.upper() for x in bern)))
    check('actual_v_metric_gradient_guard',max_grad2<arb('.25'),max_grad2)
    check('actual_positive_mass_gap',M1>18)
    moments={}; moment_rows=[]
    for power in range(2,11):
        moments[power],row=moment(power,C,q); moment_rows.append(row)
        check('complete_positive_mass_moment_'+str(power),moments[power]>0)
    def bose(power):
        m=(power+1)//2; b=abs(sp.bernoulli(2*m))/(4*m)
        return arb(int(sp.numer(b)))/int(sp.denom(b))
    lower=sum(((-1)**j*bose(5+2*j)*moments[j+2] for j in range(8)),arb(0))
    upper=sum(((-1)**j*bose(5+2*j)*moments[j+2] for j in range(9)),arb(0))
    cp=hull(lower,upper)*8/(arb.pi()*B2*Fl)
    contact=cp*arb('6.704e-12')
    check('complete_product_contact_positive',cp>0)
    check('complete_product_contact_below_reserved_precision',contact<arb('1e-43'),contact)
    fractions={power:heat_fraction(power,C,q,moments[power]) for power in [2,3]}
    for power,v in fractions.items():check('positive_time_weight_fraction_'+str(power),v>0 and v<1)
    K=2*arb('1.2e-8')*lam1**4*moments[4]/(B2*Fl)
    threshold=arb('5e-43')/K
    a=r=arb('.01'); b=s=arb('.2')
    a1=r*a; b1=s*a+r*b
    before=K*(2*a*b+kappa*a*a); after=K*(2*a1*b1+kappa*a1*a1)
    check('conditional_uncorrected_guard_exceeds_full_screen',before>arb('1.287e-40'),before)
    check('conditional_corrected_guard_meets_residual_reserve',after<arb('5e-43'),after)
    check('unchanged_sufficient_total_error_allocation',arb('1.5483364e-43')+arb('5e-43')+arb('2.5e-43')<arb('1e-42'))
    maximum_a={str(k):(threshold/(2*k+kappa)).sqrt() for k in [1,5,20]}
    NUMERICS.append(dict(kind='complete_mass_moments',B2=encoded(B2),C=encoded(C),q=str(p['q']),F_lambda=encoded(Fl),
        lambda1=encoded(lam1),M1=encoded(M1),rows=moment_rows,physical_modes_discarded=0))
    NUMERICS.append(dict(kind='complete_product_contact',coefficient=encoded(cp),
        coefficient_positive=True,projection_formula='-C_prod integral h v dz',
        actual_pairing_signed_value_evaluated=False,pairing_upper='6.704e-12',complete_projection_allowance=encoded(contact),
        product_spectral_tail='Global alternating Stieltjes inequalities at orders 7 and 8, no spectral cutoff',
        angular_tail='Monotone integral lower and upper enclosures after N=2048',actual_disk_replaced=False))
    NUMERICS.append(dict(kind='conditional_resolvent_targets',actual_parameters=True,
        max_v_gradient_squared=encoded(max_grad2),kappa1=encoded(kappa),K_res=encoded(K),bracket_threshold=encoded(threshold),
        maximum_a_at_b_over_a={k:encoded(v) for k,v in maximum_a.items()},
        design_guards={'a':'.01','b':'.2','r':'.01','s':'.2'},
        uncorrected_conditional_allowance=encoded(before),corrected_a=encoded(a1),corrected_b=encoded(b1),
        corrected_conditional_allowance=encoded(after),actual_column_and_operator_guards_verified=False,
        residual_error_budget='5e-43',other_numerical_tail_budget='2.5e-43',old_heat_cap_upper='1.5483364e-43',
        full_projection_target='1.287e-40',full_projection_evaluated=False))
    NUMERICS.append(dict(kind='positive_time_mass_weights',split='.001',
        beyond_split_fractions={str(k):encoded(v) for k,v in fractions.items()},
        interpretation='Complete positive mass-weight fractions only; no bound on the unknown actual heat response.'))
    actual_operator_certificate(M1,K,kappa)
    print(json.dumps(dict(event='contact_and_resolvent_targets',M1=str(M1),product_coefficient=str(cp),
        complete_contact_allowance=str(contact),K_res=str(K),threshold=str(threshold),before=str(before),after=str(after),
        positive_time_fractions={str(k):str(v) for k,v in fractions.items()},complete_signed_projection=False)),flush=True)

def main():
    ctx.prec=512
    old=json.loads((root/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    review=json.loads((root/OWN[2]).read_bytes())
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|
                 {'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log','checks/current_disk_heat_propagation.py'})
    before={p:digest(p) for p in paths}; deps=bh.prior.scalar.prev.endpoint.dependencies(); error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,
        invocation=sys.argv,dependencies=deps,frozen_sources_utf8={p:(root/p).read_text(encoding='utf-8') for p in OWN},protocol=review['protocol'])
    out=root/'evidence'/(STEM+'.json'); receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};receipt['runs'].append(run)
    def save():out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    save()
    try:
        check('previous_scoped_execution_passed',old['status']=='passed')
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('interval_library_pinned',bh.prior.scalar.prev.prior.flint.__version__=='0.9.0')
        algebra(); execute()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}; depafter=bh.prior.scalar.prev.endpoint.dependencies()
    status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,
        checks=ROWS,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),
        complete_product_contact_bounded=status=='passed',actual_resolvent_column_bounds_proved=status=='passed',
        continuum_argument_scope='Ordinary selected-domain resolvent and quadratic-form argument; checks only verify algebra and numerical inequalities.',
        complete_signed_projection_evaluated=False,full_projection_milestone_complete=False,
        rough_family_C2_action_proved=False,controlled_solution_certified=False,goal_complete=False)
    save()
    with (root/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:
            log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
    if status!='passed':raise SystemExit(1)

if __name__=='__main__':main()
