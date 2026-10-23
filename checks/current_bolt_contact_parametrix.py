"""Contact-complete local bolt parametrix; not the physical adjoint solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import sys, json, traceback
import sympy as sp
import mpmath as mp
sys.path.insert(0, str(Path(__file__).resolve().parent))
import current_physical_york_dual_assembly as york
principal=york.principal; endpoint=york.endpoint; prior=york.prior
from flint import arb, ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-bolt-contact-parametrix'
PREVIOUS='current-physical-york-dual-assembly'
OWN=['notes/'+STEM+'.md','checks/current_bolt_contact_parametrix.py',
     'evidence/'+STEM+'-review.json']
ROWS=[]; NUMERICS=[]; FIXTURES=[]; encoded=prior.encoded

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

prior.check=check; york.check=check

def algebra():
    C,N,p,K,E,a,q=sp.symbols('C N p K E a q',positive=True)
    check('TT_hankel_action_conversion',sp.simplify(sp.Rational(9,2)*sp.pi*N*C/(30*sp.pi)-3*N*C/20)==0)
    check('exact_Weyl_contact_not_inequality_reserve',sp.Rational(31,30)+sp.Rational(1,6)==sp.Rational(6,5))
    check('Einstein_to_log_ratio',sp.simplify((3*sp.pi*N*E/2)/(3*N*C/20)-10*sp.pi*E/C)==0)
    symbol=K*p**4*(sp.log(p)+a)+E*p*p
    check('both_TT_point_source_factors',sp.simplify(symbol*(q/2)/symbol-q/2)==0)
    r=sp.symbols('r',positive=True)
    lap=lambda v:sp.diff(v,r,2)+sp.diff(v,r)/r
    check('scalar_radial_delta_flux',sp.simplify(r*sp.diff(lap(r*r*sp.log(r)),r)-4)==0)
    y,b=sp.symbols('y b',positive=True)
    check('short_interval_log_moment',sp.integrate(y*sp.log(b/y),(y,0,b))==b*b/4)
    A,B,P=sp.symbols('A B P',positive=True)
    tail=A*A/(2*P**2)+A*B/(2*P**4)+B*B/(6*P**6)
    check('residual_tail_primitive',sp.simplify(-sp.diff(tail,P)-P*(A/P**2+B/P**4)**2)==0)
    middle=A*A*P**2/2+2*A*B*sp.log(P)-B*B/(2*P**2)
    check('finite_frequency_model_norm_primitive',sp.simplify(sp.diff(middle,P)-P*(A+B/P**2)**2)==0)
    # A finite contact mismatch cannot be discarded as an L2 error.
    check('contact_mismatch_L2_divergence',sp.limit(p/(sp.log(p)**2),p,sp.oo)==sp.oo)
    check('complex_digamma_sector_reserve',arb(2).sqrt()/6<arb(1)/4)

def parameters():
    p,h,zh,env=york.actual_parameters()
    NUMERICS.extend(york.NUMERICS)
    k=arb(2000); C0=(1+12*arb.pi()*k**6).sqrt()
    b02=k**10*(C0-1)/6
    L=(p['B2']/b02).log()/2
    a=p['C'].log()/2-L-(1+h).log()-arb(6)/5
    check('positive_exact_bolt_contact',a>10 and a<12,encoded(a))
    R2=p['B2']*(1+h)**2;ell=(1+h)/h
    A=10*(3*p['q']+5)/p['C']
    B=30*(3+arb(3*p['q']+5)/p['q']**2)*arb(p['q'])**4/p['C']**2
    threshold=2*arb(p['q'])/p['C'].sqrt()
    NUMERICS.append(dict(kind='exact_bolt_contact_and_angular_tail',a_b=encoded(a),
        L=encoded(L),h=encoded(h),ell_b=encoded(ell),R2=encoded(R2),
        A=encoded(A),B_tail=encoded(B),valid_above=encoded(threshold),
        exact_contact_definition='log(C)/2-log(B/b0)-log(1+h)-6/5'))
    return p,h,a,R2,ell,A,B,threshold

def asymptotic_checks(p):
    q=p['q'];t0=arb(q)**2
    for factor in [1,2,100,100000000]:
        t=t0*factor
        val=-principal.F(t,q)-2*t.log();bound=arb(3*q+5)/t
        NUMERICS.append(dict(kind='complex_digamma_tail_fixture',factor=factor,t=encoded(t),
            actual_difference=encoded(val),analytic_absolute_bound=encoded(bound)))
        check('complete_digamma_asymptotic_'+str(factor),val.abs_upper()<bound.lower())

def crossover(a,D,N):
    lo=arb(30);hi=arb(50)
    # Separate monotone extreme-parameter roots, including the whole seed box.
    for _ in range(170):
        mid=(lo+hi)/2;f=(2*mid).exp()*(mid+a)-D
        if f<0:lo=mid
        elif f>0:hi=mid
        else:break
    check('crossover_lower_sign_'+str(N),(2*lo).exp()*(lo+a)-D<0)
    check('crossover_upper_sign_'+str(N),(2*hi).exp()*(hi+a)-D>0)
    bracket=prior.hull(lo,hi)
    check('crossover_resolved_'+str(N),hi-lo<arb('1e-18'))
    return bracket,bracket.exp()

def value_at_zero(a,D,ell,KE,N,cells):
    step=arb(80)/cells;lower=arb(0);upper=arb(0);records=[]
    def f(s):return 1/(1+(2*s).exp()*(s+a)/D)
    left=f(arb(0))
    for i in range(cells):
        right=f((i+1)*step)
        lower+=step*right;upper+=step*left
        records.append(dict(index=i,left=encoded(left),right=encoded(right),
            lower_contribution=encoded(step*right),upper_contribution=encoded(step*left)))
        left=right
    tail=D*(-arb(160)).exp()/(2*(80+a))
    prefactor=ell/(2*KE)
    value=prefactor*prior.hull(lower.lower(),(upper+tail).upper())
    width=(value.upper()-value.lower())/value.lower()
    NUMERICS.append(dict(kind='local_TT_parametrix_bolt_value',exteriors=N,cells=cells,
        lower_sum=encoded(lower),upper_sum=encoded(upper),tail_upper=encoded(tail),
        prefactor=encoded(prefactor),value=encoded(value),relative_width=encoded(width),
        cell_records=records,actual_adjoint_value=False))
    check('finite_local_parametrix_'+str(N)+'_'+str(cells),value.is_finite() and value>0)
    if cells==4096:check('local_parametrix_width_'+str(N),width<arb('.001'),encoded(width))
    return value

def mpnum(x):return mp.mpf(x.mid().str(70,radius=False))

def independent_value(a,D,ell,KE):
    out=[]
    for digits in [70,100]:
        with mp.workdps(digits):
            aa,dd,ee,ke=map(mpnum,[a,D,ell,KE])
            val=ee/(2*ke)*mp.quad(lambda s:1/(1+mp.exp(2*s)*(s+aa)/dd),[0,20,30,35,40,50,80,mp.inf])
            out.append(mp.nstr(val,digits-10))
    return out

def numerical(p,h,a,R2,ell,A,B,threshold):
    values=[]
    for N in [1,2]:
        KT=3*N*p['C']/20;KE=3*arb.pi()*N*p['Fl']*R2/2;D=KE/KT
        DR=96*arb.pi()*arb('1e22')-1/(60*arb.pi());KS=6*arb.pi()*N*DR*p['C']
        cs=-ell/(4*KS)
        s_cross,p_cross=crossover(a,D,N)
        NUMERICS.append(dict(kind='physical_local_parametrix_coefficients',exteriors=N,
            KT=encoded(KT),KE=encoded(KE),KS=encoded(KS),D=encoded(D),scalar_singularity=encoded(cs),
            log_crossover=encoded(s_cross),crossover=encoded(p_cross),physical_pole=False))
        coarse=value_at_zero(a,D,ell,KE,N,1024)
        fine=value_at_zero(a,D,ell,KE,N,4096)
        check('monotone_grids_overlap_'+str(N),coarse.overlaps(fine))
        diagnostic=independent_value(a,D,ell,KE)
        check('independent_bolt_quadrature_'+str(N),all(fine.contains(arb(v)) for v in diagnostic))
        NUMERICS.append(dict(kind='independent_bolt_value_diagnostic',exteriors=N,precision_digits=[70,100],values=diagnostic,
            interval_certificate=False))
        # Exact full angular principal model; no finite angular cutoff.
        # Scale logs cancel in this difference before interval evaluation.
        delta0=-principal.F(arb(0),p['q'])-2*p['C'].log()+arb(62)/15
        M0=abs(delta0)+p['C']/p['q']+4*threshold.log()
        low_relative=M0*threshold**2/(4*D)
        high_relative=(A+B/threshold**2)/(4*D)
        relative_upper=max(low_relative.upper(),high_relative.upper())
        check('full_frozen_multiplier_relative_error_'+str(N),relative_upper<arb('1e-12'))
        split=arb('1e16')
        low_square=M0**2*(threshold**6-1)/(96*D**2)
        middle_square=(A*A*(split**2-threshold**2)/2+2*A*B*(split/threshold).log()
            +B*B*(1/threshold**2-1/split**2)/2)/(16*D**2)
        tail_square=(A*A/(2*split**2)+A*B/(2*split**4)+B*B/(6*split**6))/(16*(split.log()+a)**2)
        norm=(low_square+middle_square+tail_square).sqrt()
        check('complete_frozen_model_L2_residual_'+str(N),norm>0 and norm.is_finite())
        NUMERICS.append(dict(kind='complete_all_angular_frozen_model_mismatch',exteriors=N,
            M0=encoded(M0),relative_bound=encoded(relative_upper),
            low_squared_norm=encoded(low_square),middle_squared_norm=encoded(middle_square),tail_squared_norm=encoded(tail_square),
            unit_source_radial_L2_norm_upper=encoded(norm),physical_radial_L2_norm_upper=encoded(ell*norm/2),
            physical_operator_Z0=False,physical_Vdual_norm=False,curved_remainder_included=False))
        values.append(fine)
    check('identical_exterior_scaling',values[0].overlaps(2*values[1]))
    for exponent in [12,16,20]:
        P=arb(10)**exponent
        check('tail_threshold_'+str(exponent),P>threshold)
        square=(A*A/(2*P**2)+A*B/(2*P**4)+B*B/(6*P**6))/(16*(P.log()+a)**2)
        norm=square.sqrt();physical=ell*norm/2
        NUMERICS.append(dict(kind='all_angular_model_point_residual_tail',cutoff=encoded(P),
            unit_source_radial_L2_norm_upper=encoded(norm),physical_radial_L2_norm_upper=encoded(physical),
            physical_Vdual_norm=False,curved_remainder_included=False))
        check('finite_complete_angular_residual_tail_'+str(exponent),square>0 and square.is_finite())

def main():
    ctx.prec=256
    review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|
        {'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};deps=endpoint.dependencies();error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='no_progress',
        protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_conversion_passed',old['status']=='passed' and old['frozen_physical_York_transpose_executed'])
        check('previous_scientific_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        algebra();data=parameters();asymptotic_checks(data[0]);numerical(*data)
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths};depafter=endpoint.dependencies()
    status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,
        finished_utc=datetime.now(timezone.utc).isoformat(),exact_bolt_fourth_contact_identified=status=='passed',
        contact_complete_local_parametrix_executed=status=='passed',all_angular_principal_point_residual_tail_enclosed=status=='passed',
        complete_all_angular_frozen_model_mismatch_enclosed=status=='passed',
        actual_physical_correction_computed=False,full_solution_feasibility_test_executed=False,
        full_adjoint_residual_norm_bounded=False,full_forward_residual_norm_bounded=False,
        complete_renormalized_operator=False,full_angular_curved_remainder_certified=False,
        metric_admission_proved=False,nonlinear_solution_constructed=False,physical_stability_proved=False,
        semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),
                    *ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
    if status!='passed':raise SystemExit(1)

if __name__=='__main__':main()
