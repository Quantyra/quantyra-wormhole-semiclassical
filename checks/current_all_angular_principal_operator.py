"""Complete angular principal model, not the physical metric operator or solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import sys, json, traceback
import sympy as sp
import mpmath as mp
sys.path.insert(0, str(Path(__file__).resolve().parent))
import current_endpoint_thermal_restoration as endpoint
prior = endpoint.prior
from flint import arb, acb, ctx

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-all-angular-principal-operator'
PREVIOUS = 'current-endpoint-thermal-restoration'
OWN = ['notes/'+STEM+'.md', 'checks/current_all_angular_principal_operator.py',
       'evidence/'+STEM+'-review.json']
ROWS = []; NUMERICS = []; FIXTURES = []
encoded = prior.encoded

def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()

def check(name, ok, detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name, passed=bool(ok), detail=str(detail)))
    if not ok: raise AssertionError(name)

prior.check = check

def algebra():
    y = sp.symbols('y', positive=True)
    r = sp.sqrt(1-y)
    W = -sp.Rational(2,3)*r**3+sp.Rational(2,5)*r**5
    J = -2*(r**3/3-r**5/5)*sp.log(y)-sp.Rational(4,25)*r**5+sp.Rational(8,45)*r**3+sp.Rational(8,15)*r-sp.Rational(8,15)*sp.atanh(r)
    check('weight_primitive', sp.simplify(sp.diff(W,y)-y*r)==0)
    check('log_weight_primitive', sp.simplify(sp.diff(J,y)-y*r*sp.log(y))==0)
    alpha=sp.symbols('alpha'); a=alpha*(1-alpha)
    check('fourth_contact_normalization', sp.integrate(a*a,(alpha,0,1))==sp.Rational(1,30))
    check('sixth_remainder_normalization', sp.integrate(a**3,(alpha,0,1))==sp.Rational(1,140))
    n,q=sp.symbols('n q',positive=True)
    check('telescoping_fourth_tail',sp.factor(2*(q+2*n)/(n**2*(n+q)**2)-2/q*(1/n**2-1/(n+q)**2))==0)
    check('finite_cutoff_fourth_coefficient',sp.factor(2*(q+2*n)/(n*(n+q))-2/n-2/(n+q))==0)
    check('counting_reserve',sp.expand(4*n*(n+q)-2*n*(q+n+1)-2*n*(n+q-1))==0)
    check('single_integral_normalization',sp.integrate(y*sp.sqrt(1-y),(y,0,1))==sp.Rational(4,15))
    meanlog=2*sp.log(2)-sp.Rational(31,15)
    check('principal_high_frequency_constant',sp.simplify(-2*sp.log(4)+2*meanlog+sp.Rational(62,15))==0)

def F(t,q):
    """Real enclosure for the whole angular finite part at real t>=0."""
    assert t>=0
    disc=acb(arb(q*q)-4*t).sqrt()
    nuplus=(q+disc)/2
    numinus=2*t/(q+disc)
    val=-2*((1+nuplus).digamma()+(1+numinus).digamma())
    assert val.imag.contains(0)
    assert val.real.is_finite()
    return val.real

def W(y):
    if y==1: return arb(0)
    r=(1-y).sqrt()
    return -2*r**3/3+2*r**5/5

def J(y):
    if y==1: return arb(0)
    r=(1-y).sqrt()
    return -2*(r**3/3-r**5/5)*y.log()-4*r**5/25+8*r**3/45+8*r/15-8*r.atanh()/15

def quadrature_grid():
    delta=arb(2).log()/16
    nodes=[]
    for k in range(320,-1,-1):
        s=-k*delta
        y=arb(1) if k==0 else s.exp()
        nodes.append((y,s))
    weights=[]; good=True
    for (a,sa),(b,sb) in zip(nodes,nodes[1:]):
        ww=W(b)-W(a); jj=J(b)-J(a)
        left=(sb*ww-jj)/delta; right=(jj-sa*ww)/delta
        good=good and left>0 and right>0
        weights.append((left,right))
    tip=W(nodes[0][0])+arb(4)/15
    mass=sum((x+y for x,y in weights),arb(0))
    check('positive_exact_log_interpolation_weights',good)
    check('quadrature_weight_mass', (mass+tip).contains(arb(4)/15))
    ebulk=delta**2/2*arb(15)/4*mass
    etip=arb(15)/4*nodes[0][0]**2
    NUMERICS.append(dict(kind='frequency_uniform_quadrature_contract',delta=encoded(delta),
        tip_weight=encoded(tip),bulk_error=encoded(ebulk),tip_error=encoded(etip),
        nodes=[dict(y=encoded(y),log_y=encoded(s)) for y,s in nodes],
        weights=[[encoded(a),encoded(b)] for a,b in weights]))
    return nodes,weights,tip,ebulk,etip

def kappa(p,par,L,grid,save=True):
    if p==0: return -F(arb(0),par['q'])-4*L
    nodes,weights,tip,ebulk,etip=grid
    values=[-F(par['C']*p*p*y/4,par['q'])-4*L for y,s in nodes]
    total=tip*values[0]
    for i,(left,right) in enumerate(weights):
        total+=left*values[i]+right*values[i+1]
    center=arb(15)/4*total-etip/2
    result=center+arb(0,(ebulk+etip/2).abs_upper())
    if save:
        NUMERICS.append(dict(kind='complete_angular_principal_coefficient',p=encoded(p),
            kappa=encoded(result),quadrature_center=encoded(center),
            node_values=[encoded(x) for x in values],physical_operator=False))
    return result

def mpnum(x): return mp.mpf(x.mid().str(65,radius=False))

def independent(p,C,q,L):
    def fm(t):
        disc=mp.sqrt(q*q-4*t)
        plus=(q+disc)/2; minus=2*t/(q+disc)
        return -2*mp.re(mp.digamma(1+plus)+mp.digamma(1+minus))
    return mp.mpf(15)/4*mp.quad(lambda y:y*mp.sqrt(1-y)*(-fm(C*p*p*y/4)-4*L),[0,mp.mpf('.01'),mp.mpf('.1'),mp.mpf('.5'),1])

def low_modes(p,par,L,N):
    C=par['C'];q=par['q'];p2=p*p
    value=-F(arb(0),q)-4*L
    for n in range(1,N+1):
        m2=arb(n*(n+q))/C;d=2*(q+2*n);x=p2/m2
        b=2*(((4+x)/x).sqrt()*(x/(4+x)).sqrt().atanh()-1)
        value+=60*d*m2/(C*p2*p2)*(b-x/6+x*x/60)
    s4=2*C*C/q*(arb(2).zeta(N+1)-arb(2).zeta(N+q+1))
    tail=p2*s4/(7*C)
    return value,tail

def execution():
    par=prior.parameters();q=par['q'];C=par['C']
    k=arb(2000);C0=(1+12*arb.pi()*k**6).sqrt();b02=k**10*(C0-1)/6
    L=(par['B2']/b02).log()/2
    F0=F(arb(0),q);Feff=F0+4*L
    check('fixed_finite_part_at_zero',F0.overlaps(2*arb.const_euler()-2*arb(q+1).digamma()))
    check('positive_principal_minimum',-Feff>50)
    atcoalescence=F(arb(q*q)/4,q)
    check('root_coalescence_is_removable',atcoalescence.overlaps(-4*arb(1+q//2).digamma()))
    around=F(arb(q*q)/4+arb(0,'1e-30'),q)
    check('root_coalescence_interval_finite',around.is_finite() and around.overlaps(atcoalescence))
    grid=quadrature_grid();eps=grid[3]+grid[4]/2
    rel=eps/(-Feff-eps)
    check('uniform_principal_preconditioner_relative_target',rel<arb('2e-5'),encoded(rel))
    NUMERICS.append(dict(kind='current_parameters_and_uniform_bound',C=encoded(C),q=q,
        log_mu0B=encoded(L),F_eff=encoded(Feff),coefficient_error=encoded(eps),
        principal_relative_error=encoded(rel),full_physical_Z0_certified=False))
    fixtures=[('one',arb(1)),('twenty',arb(20)),('ten_thousand',arb('1e4')),
              ('coalescence_scale',q/C.sqrt()),('one_eight',arb('1e8')),
              ('one_sixteen',arb('1e16')),('one_thirtytwo',arb('1e32'))]
    for name,p in fixtures:
        val=kappa(p,par,L,grid)
        check('coefficient_width_'+name,val.upper()-val.lower()<arb('.002'))
        diag=independent(mpnum(p),mpnum(C),q,mpnum(L))
        check('independent_real_integral_'+name,val.contains(arb(mp.nstr(diag,62))))
        NUMERICS.append(dict(kind='independent_real_integral_diagnostic',fixture=name,
            p=encoded(p),mpmath_value=mp.nstr(diag,62),precision_digits=60,rigorous=False))
        if name in ['one','twenty']:
            low,tail=low_modes(p,par,L,128)
            cover=low+arb(0,tail.abs_upper()/2)+tail/2
            check('finite_modes_plus_exact_angular_tail_'+name,cover.overlaps(val))
            check('independent_integral_in_lowmode_tail_'+name,cover.contains(arb(mp.nstr(diag,62))))
            NUMERICS.append(dict(kind='finite_modes_with_exact_fourth_tail',p=encoded(p),
                N=128,low_coefficient=encoded(low),positive_tail_upper=encoded(tail),cover=encoded(cover)))
        print(json.dumps(dict(event='fixture_complete',name=name,kappa=val.str(12))),flush=True)
    for N in [2,10**6,10**16]:
        HN=arb(N+1).digamma()+arb.const_euler()
        HNq=arb(N+q+1).digamma()+arb.const_euler()
        Hq=arb(q+1).digamma()+arb.const_euler()
        ceiling=-Feff+2*(HN+HNq-Hq)
        p_quarter=(L+arb(38)/15+ceiling/3).exp()/C.sqrt()
        p=2*p_quarter
        full=kappa(p,par,L,grid)
        omitted=1-ceiling/full
        check('fixed_cutoff_quarter_failure_'+str(N),omitted>arb('.25'),encoded(omitted))
        check('witness_coefficient_width_'+str(N),full.upper()-full.lower()<arb('.002'))
        NUMERICS.append(dict(kind='fixed_angular_cutoff_obstruction_witness',N=N,
            cutoff_coefficient_ceiling=encoded(ceiling),quarter_threshold_p=encoded(p_quarter),
            witness_p=encoded(p),full_coefficient=encoded(full),
            certified_omitted_fraction_lower=encoded(omitted),
            scope='principal model only; not full lower-order physical Hessian'))
    FIXTURES.append(dict(fixture='whole_angular_principal_model',uniform_frequency_proof=True,
        normalized_coefficient_width_target='.002',relative_model_target='2e-5',full_curved_operator=False))
    FIXTURES.append(dict(fixture='fixed_angular_cutoff_uniform_form_obstruction',
        ordinary_proof_for_every_finite_N=True,finite_witnesses=3,physical_nonexistence=False))

def main():
    ctx.prec=256;mp.mp.dps=60
    review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|
        {'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};deps=endpoint.dependencies();error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
        previous_goal_turn='paused_at_user_request',protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('accepted_endpoint_result_retained',old['status']=='passed' and old['complete_local_thermal_restoration_implemented'])
        check('previous_scientific_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        algebra();execution()
    except Exception: error=traceback.format_exc()
    after={p:digest(p) for p in paths};depafter=endpoint.dependencies()
    status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,
        finished_utc=datetime.now(timezone.utc).isoformat(),all_angular_principal_model_enclosed=status=='passed',
        fixed_angular_cutoff_uniform_form_obstruction_executed=status=='passed',
        complete_renormalized_operator=False,full_angular_curved_remainder_certified=False,
        actual_physical_correction_computed=False,physical_restorations_complete=False,
        full_solution_feasibility_test_executed=False,metric_admission_proved=False,
        nonlinear_solution_constructed=False,physical_stability_proved=False,
        semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json')
    receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run)
    out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),
                    *ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
    if status!='passed': raise SystemExit(1)

if __name__=='__main__': main()
