"""Supporting arithmetic for exact all-angular quantum operator assembly.

No retained signed center, physical correction or mass-column certificate.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb,factorial
import sys,json,traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb,ctx
import sympy as S
ROOT=Path(__file__).resolve().parents[1];STEM='current-angular-resummed-quantum-operator'
OWN=['notes/'+STEM+'.md','checks/current_angular_resummed_quantum_operator.py','evidence/'+STEM+'-review.json',
 'operators/angular_resummed_response.py']
INPUTS=OWN+[
 'notes/current-canonical-heat-restoration-assembly.md',
 'notes/current-curved-paired-heat-bridge.md',
 'notes/current-rough-paired-variation.md',
 'notes/current-third-metric-conversion.md',
 'notes/current-strong-quotient-reconstruction-test.md',
 'notes/current-global-mass-conversion.md',
 'notes/current-coupled-global-return-attempt.md',
 'evidence/current-adaptive-dirac-trace-remainder.json']
sys.path.insert(0,str(ROOT))
from operators.angular_resummed_response import AngularResponse
from flint import acb
CHECKS=[]
def check(name,ok):
    assert name not in {x['name'] for x in CHECKS}
    CHECKS.append({'name':name,'passed':bool(ok)})
    if not ok:raise AssertionError(name)
def enc(x):
    m,r,e=x.mid_rad_10exp(55)
    return {'ball':str(x),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
    d=x['decimal_enclosure'];return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])
def prior(stem):
    r=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
    assert r['status']=='passed';return r['result']
def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    for p,h in review['scientific_sha256'].items():
        check('frozen_'+p,sha256((ROOT/p).read_bytes()).hexdigest()==h)
    pp=prior('current-adaptive-dirac-trace-remainder')['actual_parameters']
    C,B2=[dec(pp[k]) for k in ['C','B2']]
    k=arb(2000);q=k**4
    b02=k**10*((1+12*arb.pi()*k**6).sqrt()-1)/6
    ell=(B2/b02).log()/2
    ar=AngularResponse(C,q,ell)
    check('fixed_physical_scale_preserved',ell>arb('2.94524e-9') and ell<arb('3e-9'))
    R0,B0=ar.resolvent_and_diagonal(0)
    Hq=(q+1).digamma()+arb.const_euler()
    expected=C*(4*arb.const_euler()-2*Hq+4*ell)
    check('common_angular_finite_part_at_zero',R0.overlaps(expected))
    upper=ar.global_bubble_upper()
    check('entire_angular_zero_argument_bubble',B0.overlaps(upper))
    check('whole_spectrum_bubble_upper',upper>arb('4.9e8') and upper<arb('5e8'))
    crit=q*q/(4*C);rr=1+q/2
    Rc,Bc=ar.resolvent_and_diagonal(crit)
    check('coalescent_resolvent_limit',Rc.overlaps(4*C*(ell-rr.digamma())))
    check('coalescent_bubble_limit',Bc.overlaps(-2*C*C*acb(rr).polygamma(2).real))
    args=[('zero',arb(0)),('one',arb(1)),('324',arb(324)),('1e6',arb('1e6')),
          ('half_coalescence',crit/2),('coalescence',crit),('twice_coalescence',2*crit),('1e24',arb('1e24'))]
    vals=[]
    for tag,lam in args:
        r,b=ar.resolvent_and_diagonal(lam)
        check('positive_complete_bubble_'+tag,b>0)
        check('whole_spectrum_guard_'+tag,b<arb('5e8'))
        vals.append({'argument_label':tag,'argument':enc(lam),'finite_part_resolvent':enc(r),'diagonal_bubble':enc(b)})
    mixed=[]
    for a,b in [(arb(0),arb(1)),(arb(324),arb('1e6')),(crit/2,crit*2)]:
        v=ar.bubble(a,b);rev=ar.bubble(b,a)
        check('symmetric_bubble_'+str(len(mixed)),v.overlaps(rev) and v>0 and v<arb('5e8'))
        mixed.append({'lambda':enc(a),'mu':enc(b),'complete_angular_bubble':enc(v)})
    ball=arb(1,'1e-20');vb=ar.bubble(ball,arb(1))
    check('overlapping_argument_enclosure',vb.contains(ar.resolvent_and_diagonal(1)[1]))
    n,qs,cs,lam=S.symbols('n q C lambda',positive=True)
    check('exact_angular_telescoping',S.simplify((qs+2*n)/(n*n*(n+qs)**2)-(1/n**2-1/(n+qs)**2)/qs)==0)
    x,t=S.symbols('x t',positive=True)
    f=2*(qs+2*x)*S.exp(-t*x*(x+qs)/cs)
    em=-f.subs(x,0)/2-S.diff(f,x).subs(x,0)/12+S.diff(f,x,3).subs(x,0)/720
    poly=S.series(em,t,0,2).removeO()
    check('angular_heat_coefficients',S.simplify(poly-(-qs-S.Rational(1,3)+(qs**2/S.Integer(6)-S.Rational(1,30))*t/cs))==0)
    check('next_endpoint_starts_at_t_squared',S.series(S.diff(f,x,5).subs(x,0),t,0,2).removeO()==0)
    ss=S.symbols('s');z0,z1=S.symbols('z0 z1');la=S.symbols('ellA')
    # Taylor coefficients of the three original common Mellin factors.
    check('shifted_Mellin_first_factor',S.diff(S.exp(2*la*ss)*(z0+ss*z1)/(ss-1),ss).subs(ss,0)==-2*la*z0-z0-z1)
    check('shifted_Mellin_second_factor',S.diff(S.exp(2*la*ss)*(z0+ss*z1),ss).subs(ss,0)==2*la*z0+z1)
    res,fp=S.symbols('res fp')
    check('shifted_Mellin_pole_contact',S.diff(S.exp(2*la*ss)*(res+ss*fp),ss).subs(ss,0)==2*la*res+fp)
    check('conditional_fast_screen_unchanged',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    coeff=ar.common_heat_regrouping()
    return {'precision_bits':ctx.prec,'actual_parameters':{'C':enc(C),'B2':enc(B2),'q':enc(q),'log_mu_B':enc(ell)},
     'entire_angular_resolvent_finite_part_at_zero':enc(R0),'entire_angular_bubble_uniform_upper_expression':enc(upper),
     'common_heat_regrouping_coefficients':{k:enc(v) for k,v in coeff.items()},
     'spectral_function_evaluations':vals,'two_argument_evaluations':mixed,
     'scope':{'complete_angular_identities_executed':True,'common_regulator_action_regrouping_proved':True,
      'actual_paired_domain_preserved':True,'all_thermal_modes_retained_in_actual_orbit_operator':True,
      'orbit_operator_application_executed':False,'strong_quantum_response_bound_certified':False,
      'complete_signed_coupled_solve_executed':False,'actual_kappa_m_computed':False,
      'physical_X_comparisons':False,'nonlinear_admission':False,'nonlinear_contraction':False,
      'exact_compatibility':False,'physical_matching':False,'semiclassical_validity_noise':False,
      'all_sector_stability':False,'actual_controlled_solution':False,'construction_obstruction':False,
      'full_goal_complete':False}}

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
      'frozen_sources_utf8':{f:(ROOT/f).read_text(encoding='utf-8') for f in OWN},
      'runtime':{'python':sys.version,'executable_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
                 'flint_version':flint.__version__,'sympy_version':S.__version__}}
    try:
        run['result']=execute();run['input_sha256_after']={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    message=json.dumps({'mode':mode,'status':run['status'],'supporting_checks':len(CHECKS),
       'angular_bubble_upper':run.get('result',{}).get('entire_angular_bubble_uniform_upper_expression'),
       'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
