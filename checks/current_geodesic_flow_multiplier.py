"""Full flow multiplier with exact contact moments; not a global metric solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,traceback
import sympy as sp
import mpmath as mp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_geodesic_two_sided_remainder as previous
york=previous.york;prior=previous.prior;endpoint=previous.endpoint;principal=previous.principal
from flint import arb,acb,ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-geodesic-flow-multiplier';PREVIOUS='current-geodesic-two-sided-remainder'
OWN=['notes/'+STEM+'.md','checks/current_geodesic_flow_multiplier.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];encoded=prior.encoded

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)
prior.check=check;york.check=check

def algebra():
    m=sp.symbols('m',positive=True)
    a=(22+28*m-18*m*m)/45;b=(-22+28*m-6*m*m)/45
    check('elliptic_primitive_E_coefficient',sp.simplify(sp.diff(a,m)+a/(2*m)+b/(2*m*(1-m))-(1-m))==0)
    check('elliptic_primitive_K_coefficient',sp.simplify(sp.diff(b,m)-(a+b)/(2*m))==0)
    check('primitive_endpoint_zero',a.subs(m,0)+b.subs(m,0)==0)
    check('primitive_endpoint_one',a.subs(m,1)==sp.Rational(32,45) and b.subs(m,1)==0)
    s=sp.symbols('s',real=True)
    W=sp.pi*sp.gamma(s+2)*sp.gamma(s+3)/(2*sp.gamma(s+sp.Rational(7,2))*sp.gamma(s+sp.Rational(5,2)))
    check('Mellin_mass',sp.simplify(W.subs(s,0))==sp.Rational(32,45))
    mu=4*sp.log(2)-sp.Rational(97,30)
    check('Mellin_log_moment',sp.simplify(sp.diff(W,s).subs(s,0)/W.subs(s,0)-mu)==0)
    contact=2*mu-4*sp.log(2)
    check('angular_contact_recovery',sp.simplify(contact+4*(sp.Rational(7,12)-sp.log(2))+sp.Rational(62,15))==0)
    check('fourth_normalization',sp.Rational(3,4)/90==sp.Rational(1,120))
    check('sixth_single_mass_angular_moment',sp.Rational(1,525)*sp.Rational(5,8)==sp.Rational(1,840))
    check('absolute_form_fourth_conversion',sp.Rational(2,120)*100==sp.Rational(5,3))

def E(m):return acb(m).elliptic_e().real
def K(m):return acb(m).elliptic_k().real
def U(v):
    if v==0:return arb(0)
    if v==1:return arb(32)/45
    m=1-v
    return arb(32)/45-(22+28*m-18*m*m)/45*E(m)-(-22+28*m-6*m*m)/45*K(m)

def flow_grid():
    delta=arb(2).log()/256
    logs=[-j*delta for j in range(5120,-1,-1)]
    nodes=[s.exp() for s in logs];nodes[-1]=arb(1)
    cumulative=[U(v)*45/32 for v in nodes]
    weights=[cumulative[0]]+[b-a for a,b in zip(cumulative,cumulative[1:])]
    mids=[logs[0]]+[(a+b)/2 for a,b in zip(logs,logs[1:])]
    check('all_5121_weights_positive',all(w>0 for w in weights))
    check('exact_weight_mass',sum(weights,arb(0)).contains(arb(1)))
    mu=4*arb(2).log()-arb(97)/30
    mean=sum((w*s for w,s in zip(weights,mids)),arb(0));sigma=mu-mean
    corrected=[s+sigma for s in mids]
    check('all_shifted_nodes_inside_unit_interval',all(s<0 for s in corrected))
    check('exact_corrected_log_moment',sum((w*s for w,s in zip(weights,corrected)),arb(0)).overlaps(mu))
    bulk=delta*delta*(arb(1)/2+arb(5)/6*(5*delta/2).exp())
    tip=45*arb.pi()/64*nodes[0]**2
    eps=bulk+tip+4*sigma.abs_upper()
    check('uniform_coefficient_error_target',eps<arb('2.4e-5'),encoded(eps))
    NUMERICS.append(dict(kind='contact_corrected_flow_grid',delta=encoded(delta),mu=encoded(mu),
        raw_log_mean=encoded(mean),sigma=encoded(sigma),bulk_error=encoded(bulk),tip_error=encoded(tip),
        coefficient_error=encoded(eps),weight_mass=encoded(sum(weights,arb(0))),
        weights=[encoded(x) for x in weights],raw_log_nodes=[encoded(x) for x in mids],
        exact_formula_not_midpoint_operator=True))
    return weights,corrected,eps

def chi_num(w,p,L,grid):
    weights,logs,eps=grid
    if w==0:return -principal.F(arb(0),p['q'])-4*L
    args=[p['C']*w*w*s.exp()/4 for s in logs]
    vals=[-principal.F(t,p['q'])-4*L for t in args]
    return sum((a*b for a,b in zip(weights,vals)),arb(0))

def independent(w,C,q,L):
    def ff(v):
        t=C*w*w*v/4;disc=mp.sqrt(q*q-4*t)
        return 2*mp.re(mp.digamma(1+(q+disc)/2)+mp.digamma(1+2*t/(q+disc)))-4*L
    return mp.mpf(45)/32*mp.quad(lambda v:v*mp.ellipe(1-v)*ff(v),[0,mp.mpf('.000001'),mp.mpf('.01'),mp.mpf('.2'),1])

def execute(old):
    p,h,zh,env=york.actual_parameters();NUMERICS.extend(york.NUMERICS)
    C,B2,lam,Fl=p['C'],p['B2'],p['g'],p['Fl']
    k=arb(2000);C0=(1+12*arb.pi()*k**6).sqrt();b02=k**10*(C0-1)/6
    L=(B2/b02).log()/2;Feff=principal.F(arb(0),p['q'])+4*L
    check('retained_finite_seed_scale',L>0 and L<arb('3e-9') and -Feff>50)
    grid=flow_grid();eps=grid[2]
    j0=1000*C*C/(lam*B2*B2);j1=arb('2e5')*C/B2;jo=16*arb('.01')**4/25*j0
    fourth=5*eps/3;lower=eps*(500**4*jo/(60*arb.pi())+j1/(180*arb.pi()))/(arb('.99')*Fl)
    check('joint_fourth_margin_controls_quadrature',fourth>lower)
    check('quadrature_form_budget',fourth<arb('4e-5'),encoded(fourth))
    # Rounded UP from the preceding certified interval, not its midpoint.
    eta_previous=arb('.0021040238');combined=eta_previous+fourth
    check('combined_actual_form_budget_unchanged',combined<arb('.0022'),encoded(combined))
    contact=4*arb(2).log()-arb(7)/3
    check('rescaling_does_not_repair_contact',contact>arb('.4392') and contact<arb('.4393'))
    NUMERICS.append(dict(kind='actual_flow_quadrature_form_budget',C=encoded(C),B2=encoded(B2),
        lambda_value=encoded(lam),F_lambda=encoded(Fl),Lscale=encoded(L),F_eff=encoded(Feff),
        j1=encoded(j1),jo=encoded(jo),eta_quad=encoded(fourth),lower_order_ratio=encoded(lower),
        previous_eta_upper=encoded(eta_previous),combined_eta_upper=encoded(combined),
        unchanged_eta_star='.0022',wrong_normalized_flat_contact_difference=encoded(contact),
        physical_pair_computed=False,base_lift_and_actual_geodesics_still_required=True))
    mp.mp.dps=55
    fixtures=[('one',arb(1)),('twenty',arb(20)),('ten_thousand',arb('1e4')),
        ('mass_root_transition',arb(p['q'])/C.sqrt()),('one_sixteen',arb('1e16')),('one_thirtytwo',arb('1e32'))]
    for label,w in fixtures:
        val=chi_num(w,p,L,grid);cover=val+arb(0,eps.abs_upper())
        diag=independent(principal.mpnum(w),principal.mpnum(C),p['q'],principal.mpnum(L))
        check('independent_exact_integral_'+label,cover.contains(arb(mp.nstr(diag,52))))
        check('positive_flow_coefficient_'+label,cover>50)
        check('uniform_width_'+label,cover.upper()-cover.lower()<arb('4.8e-5'))
        NUMERICS.append(dict(kind='flow_frequency_fixture',label=label,w=encoded(w),chi_numerical=encoded(val),
            chi_exact_enclosure=encoded(cover),independent_mpmath=mp.nstr(diag,52),diagnostic_is_rigorous=False))
        print(json.dumps(dict(event='fixture_complete',label=label,chi=cover.str(12))),flush=True)
    check('exact_zero_contact',chi_num(arb(0),p,L,grid).overlaps(-Feff))

def main():
    ctx.prec=256;mp.mp.dps=55;review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|
        {'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};deps=endpoint.dependencies();error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',
        protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_comparison_passed',old['status']=='passed' and old['two_sided_geodesic_remainder_bounded'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        algebra();execute(old)
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths};depafter=endpoint.dependencies()
    status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        unchanged_dependencies=deps==depafter,checks=ROWS,numerics=NUMERICS,
        finished_utc=datetime.now(timezone.utc).isoformat(),
        complete_flow_multiplier_implemented=status=='passed',exact_contact_moments_preserved=status=='passed',
        actual_form_quadrature_error_bounded=status=='passed',physical_eta_star_unchanged=status=='passed',
        comparison_model_forward_adjoint_solve_executed=False,geodesic_flow_matrix_evaluated=False,
        actual_physical_correction_computed=False,full_solution_feasibility_test_executed=False,
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
