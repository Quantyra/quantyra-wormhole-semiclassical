"""Heat-cutoff small-form feasibility test, not a physical solution test."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,traceback
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_canonical_heat_restoration_assembly as rest
york=rest.york;prior=rest.prior;endpoint=rest.endpoint;principal=york.principal
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-heat-cutoff-form-obstruction';PREVIOUS='current-canonical-heat-restoration-assembly'
OWN=['notes/'+STEM+'.md','checks/current_heat_cutoff_form_obstruction.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];FIXTURES=[];encoded=prior.encoded
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)
prior.check=check;york.check=check;principal.check=check

def algebra():
    s,x,t,T,lam=sp.symbols('s x t T lambda',positive=True)
    rr=1-sp.exp(-s)-s+s*s/2
    check('positive_third_derivative',sp.diff(rr,s,3)==sp.exp(-s))
    check('triple_zero_at_origin',all(sp.diff(rr,s,k).subs(s,0)==0 for k in range(3)))
    R=sp.log(1+x)-x+x*x/2
    check('Laplace_representation_derivative',sp.simplify(sp.diff(R,x)-(1/(1+x)-1+x))==0)
    check('heat_long_moment',sp.simplify(sp.integrate(t*lam*sp.exp(-lam*t),(t,T,sp.oo))-sp.exp(-lam*T)*(T+1/lam))==0)
    a=sp.symbols('a',real=True);aa=a*(1-a)
    check('fourth_alpha_moment',sp.integrate(aa**2,(a,0,1))==sp.Rational(1,30))
    check('sixth_alpha_moment',sp.integrate(aa**3,(a,0,1))==sp.Rational(1,140))
    check('fixed_frequency_cap_coefficient',sp.Rational(12,12*140)==sp.Rational(1,140))
    ce,eta=sp.symbols('ce eta',positive=True)
    logratio=sp.Rational(38,15)+ce/(4*(1-eta))
    check('all_eta_threshold_identity',sp.simplify(4*logratio-sp.Rational(152,15)-ce/(1-eta))==0)
    z,f,xi,b=sp.symbols('r f xi b',positive=True)
    Fr=sp.Function('F')(z);Xi=sp.Function('xi')(z)
    integral=sp.Function('I')(z)
    bp=-integral/(z*z*Fr)
    w=z*z*sp.diff(Fr,z)-2*z*Fr
    check('exact_compact_TT_moment_constraint',sp.simplify(sp.diff(z*z*Fr*bp,z).subs(sp.diff(integral,z),w*Xi)+w*Xi)==0)
    eig=[-2*xi,0,xi,xi];proj=[eig[0],eig[2],eig[3]]
    check('nonzero_TT_spin_two_symbol',sp.expand(3*sum(v*v for v in proj)-sum(proj)**2)==18*xi*xi)

def expint_series(x,N=96):
    term=-x;series=term
    for k in range(2,N+1):
        term*=(-x)*(k-1)/(k*k);series+=term
    nxt=term*(-x)*N/((N+1)*(N+1))
    center=-arb.const_euler()-x.log()-series
    return center+prior.hull(min(arb(0),(-nxt).lower()),max(arb(0),(-nxt).upper())),center,nxt

def execute():
    p,h,zh,env=york.actual_parameters();NUMERICS.extend(york.NUMERICS)
    C=p['C'];k=arb(2000);C0=(1+12*arb.pi()*k**6).sqrt();b02=k**10*(C0-1)/6
    L=(p['B2']/b02).log()/2;kappa0=-principal.F(arb(0),p['q'])-4*L
    check('actual_scale_and_positive_principal_contact',L>0 and L<arb('3e-9') and kappa0>50)
    check('actual_mass_gap',arb(p['q']+1)/C>324)
    grid=principal.quadrature_grid()
    NUMERICS.append(dict(kind='heat_cutoff_parameters',C=encoded(C),B2=encoded(p['B2']),q=p['q'],
        L=encoded(L),kappa0=encoded(kappa0),angular_mass_square_gap='324',
        euler_gamma=encoded(arb.const_euler()),complete_angular_principal_model=True))
    for label,T in [('one_minus2',arb('1e-2')),('one_minus8',arb('1e-8')),('one_minus16',arb('1e-16'))]:
        x=324*T;E1=x.expint(1);independent,center,nextterm=expint_series(x)
        check('positive_E1_'+label,E1>0)
        check('independent_E1_series_'+label,E1.overlaps(independent))
        Jbar=4*((x+2)*(-x).exp()+E1);ceiling=kappa0+Jbar
        threshold=(L+arb(38)/15+ceiling/3).exp()/C.sqrt();witness=2*threshold
        full=principal.kappa(witness,p,L,grid)
        fraction=1-ceiling/full
        analytic_kappa_lower=4*((witness*C.sqrt()).log()-L)-arb(152)/15
        analytic_fraction=1-ceiling/analytic_kappa_lower
        check('strict_quarter_failure_from_analytic_bound_'+label,analytic_fraction>arb('.25'))
        check('strict_quarter_failure_all_angular_'+label,fraction>arb('.25'),encoded(fraction))
        check('retained_evaluator_width_'+label,full.upper()-full.lower()<arb('.002'))
        NUMERICS.append(dict(kind='heat_cutoff_uniform_failure_witness',label=label,T=encoded(T),x=encoded(x),
            exp_minus_x=encoded((-x).exp()),E1=encoded(E1),E1_series=encoded(independent),E1_series_center=encoded(center),
            E1_series_next_term=encoded(nextterm),E1_series_terms=96,Jbar=encoded(Jbar),retained_coefficient_ceiling=encoded(ceiling),
            quarter_threshold_p=encoded(threshold),witness_p=encoded(witness),full_kappa=encoded(full),
            omitted_fraction_lower_bound=encoded(fraction),analytic_kappa_lower=encoded(analytic_kappa_lower),
            analytic_omitted_fraction_lower=encoded(analytic_fraction),
            finite_witness_is_full_physical_Hessian=False,uniform_all_T_proof_in_note=True))
    NUMERICS.extend(principal.NUMERICS)
    for P in [arb(20),arb('1e4'),arb('1e16')]:
        Tallowed=7*kappa0/(24*P*P)
        check('fixed_band_threshold_identity_'+P.str(8),(6*P*P*Tallowed/(7*kappa0)).contains(arb('.25')))
        NUMERICS.append(dict(kind='fixed_frequency_band_comparator',P=encoded(P),sufficient_T_quarter=encoded(Tallowed),
            target='.25',physical_UV_cutoff_introduced=False,applies_to_complete_form_space=False))
    FIXTURES.append(dict(kind='decisive_method_feasibility_result',
        uniform_heat_cap_relative_supremum='1 for every T>0 in the complete principal model',
        actual_compact_core_norm_lower_bound='at least 1 if bounded; follows from local symbol and fixed-T H2 estimate',
        full_unremoved_principal_Aop_power_bound_possible=False,
        fixed_direction_heat_convergence_ruled_out=False,actual_minus_principal_remainder_bounded=False,
        physical_nonexistence_proved=False,physical_instability_proved=False,
        full_physical_residual_test_executed=False))

def main():
    ctx.prec=256;review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|
        {'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};deps=endpoint.dependencies();error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',
        protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_local_restoration_passed',old['status']=='passed' and old['canonical_local_and_LLL_assembled'])
        check('previous_scientific_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        algebra();execute()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths};depafter=endpoint.dependencies()
    status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,
        finished_utc=datetime.now(timezone.utc).isoformat(),uniform_heat_cutoff_obstruction_executed=status=='passed',
        full_unremoved_principal_Aop_bound_ruled_out=status=='passed',
        actual_minus_principal_paired_remainder_bounded=False,
        actual_physical_correction_computed=False,full_solution_feasibility_test_executed=False,
        full_adjoint_residual_norm_bounded=False,full_forward_residual_norm_bounded=False,
        complete_renormalized_operator=False,metric_admission_proved=False,
        nonlinear_solution_constructed=False,physical_stability_proved=False,
        semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),
                    *ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
    if status!='passed':raise SystemExit(1)

if __name__=='__main__':main()
