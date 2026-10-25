"""Supporting algebra/arithmetic for the actual global Dirac series.

No regulated trace, gravitational solve or physical graph error is computed.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import sys,json,traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb,ctx
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
STEM='current-global-dirac-transport-series'
OWN=['notes/'+STEM+'.md','checks/current_global_dirac_transport_series.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/current-coincident-dirac-action-assembly.md',
 'notes/stationary-massive-boundary-layer.md',
 'notes/stationary-curved-quantum-tt.md',
 'notes/current-rough-paired-variation.md',
 'notes/current-global-far-return-energy-bound.md',
 'notes/current-short-travel-strong-domain.md',
 'notes/current-strong-return-cauchy-test.md',
 'notes/current-global-mass-conversion.md',
 'evidence/current-global-mass-conversion.json',
 'notes/current-complete-weak-mass-residual.md',
 'evidence/current-complete-weak-mass-residual.json',
 'notes/current-actual-far-graph-energy-bound.md',
 'evidence/current-actual-far-graph-energy-bound.json',
 'notes/current-coupled-global-return-attempt.md']
CHECKS=[]
def check(name,ok):
    assert name not in {x['name'] for x in CHECKS}
    CHECKS.append({'name':name,'passed':bool(ok)})
    if not ok:raise AssertionError(name)
def enc(x):
    m,r,e=x.mid_rad_10exp(50)
    return {'ball':str(x),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
    d=x['decimal_enclosure']
    return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])
def prior(stem):
    r=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
    assert r['status']=='passed'
    return r['result']
def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('ordinary_argument_frozen',sha256((ROOT/OWN[0]).read_bytes()).hexdigest()==review['scientific_sha256'][OWN[0]])
    par=prior('current-global-mass-conversion')['parameters']
    beta,u,eps,tau=[dec(par[x]) for x in ['beta','u','eps','tau']]
    C=1/u;qflux=arb(2000)**4
    q2guard=12*beta+2*u+82*abs(eps)
    check('actual_seed_positive_parameters',beta>0 and u>0 and tau>0 and eps<0)
    check('actual_parameter_relation',(6*beta+u+eps-1).contains(0))
    check('whole_interval_second_derivative_upper',q2guard<arb('2.1'))
    check('all_angular_masses_above_eighteen',(1+qflux)/C>18**2)
    z,be,uu,ep,ta=S.symbols('z beta u eps tau')
    Q=be+uu*z**2+(-4*be-2*uu-ta)*z**3+(3*be+uu+ep)*z**4-2*ep*z**5+ep*z**6
    check('full_polynomial_second_derivative_identity',S.expand(S.diff(Q,z,2)-(be*(-24*z+36*z*z)+2*uu*(1-6*z+6*z*z)+ep*(12*z*z-40*z**3+30*z**4)-6*ta*z))==0)
    check('coefficient_bounds_on_unit_interval',S.expand(12-(-24*z+36*z*z)-12*(1-z)*(3*z+1))==0 and S.expand(1-(1-6*z+6*z*z)-6*z*(1-z))==0)
    # Frame calculation uses original first-order matrix, not a new square.
    sx=S.Matrix([[0,1],[1,0]]);sz=S.diag(1,-1);J=S.Matrix([[0,-1],[1,0]])
    C0=S.Matrix([[-1,-1],[1,-1]])/S.sqrt(2)
    check('original_constant_frame_orthogonal',C0.T*C0==S.eye(2))
    check('original_mass_and_frequency_transfer',C0.T*sx*C0==-sz and C0.T*sx*(-J)*C0==-sx)
    th=S.symbols('theta',real=True)
    O=S.Matrix([[S.cos(th/2),-S.sin(th/2)],[S.sin(th/2),S.cos(th/2)]])
    simp=lambda a:a.applyfunc(S.trigsimp)
    check('rotation_diagonalizes_actual_transfer',simp(O.T*(S.cos(th)*sz+S.sin(th)*sx)*O)==sz)
    check('moving_frame_term_retained',simp(O.T*O.diff(th))==J/2)
    B=simp(O.T*sz*O)
    check('fixed_frame_mass_vertex_is_isometric',simp(B.T*B)==S.eye(2))
    m,p,v,qp,gam,omega,t=S.symbols('m p v qp gamma omega t',positive=True)
    coupling=m*p*S.sqrt(v)*qp/(4*(m*m*v+p*p)**S.Rational(3,2))
    check('coupling_ratio_change_of_variable',S.simplify(coupling.subs(v,(p*t/m)**2)-qp*t/(4*p*(1+t*t)**S.Rational(3,2)))==0)
    rhs=t*t*(omega**-2+S.Rational(21,20)*t*t/m**2)/(4*(1+t*t)**3)
    from_q=((4*gam**2+S.Rational(21,5)*(p*t/m)**2)/(16*p*p)*t*t/(1+t*t)**3).subs(p,gam*omega)
    check('correlated_bolt_derivative_bound_conversion',S.simplify(rhs-from_q)==0)
    rho=S.symbols('rho',positive=True)
    alphapolar=-m*omega/(2*(m*m*rho*rho+omega*omega))
    kpolar=S.sqrt(m*m+omega*omega/rho**2)
    check('rotation_regular_at_each_fixed_bolt_mode',S.limit(alphapolar/kpolar,rho,0,dir='+')==0 and S.limit(rho*kpolar,rho,0,dir='+')==omega)
    a=arb(7)/2160;b=8-2*a
    ymax=8/(b+(b*b+16*a).sqrt())
    upper=(ymax*(4+a*ymax)/(4*(1+ymax)**3)).sqrt()
    yy,aa=S.symbols('y a',positive=True)
    f=yy*(4+aa*yy)/(4*(1+yy)**3)
    check('entire_frequency_maximum_derivative',S.simplify(S.diff(f,yy)*4*(1+yy)**4-(-aa*yy**2+(2*aa-8)*yy+4))==0)
    check('positive_unique_maximum_formula',(a*ymax*ymax+b*ymax-4).contains(0) and b>0 and ymax>0)
    check('maximum_location_enclosed',ymax>arb('.50030394253') and ymax<arb('.50030394254'))
    check('uniform_actual_coupling_contraction',upper>arb('.38497815535') and upper<arb('.38497815536') and upper<arb('.385'))
    qreal=arb('.385')+arb('.001').expm1()
    qcomplex=arb('.385')+arb('.011').expm1()
    check('admitted_supremum_family_ratio',qreal<arb('.38601'))
    check('complex_parameter_ratio',qcomplex<arb('.3961'))
    check('auxiliary_inverse_upper',1/(1-arb('.38601'))<arb('1.628692'))
    tails={}
    for n in (16,32,64,128):
        tails[str(n)]={'inverse':enc(arb('.38601')**n/(1-arb('.38601')))}
        for j in (1,2,3):tails[str(n)]['variation_'+str(j)]=enc(arb(100*j)**j*arb('.3961')**n/(1-arb('.3961')))
    for key,bound in {'inverse':'1.980e-53','variation_1':'5.472e-50','variation_2':'2.189e-47','variation_3':'1.478e-44'}.items():
        check('N128_'+key+'_guard',dec(tails['128'][key])<arb(bound))
    n,L,q=S.symbols('n L q',integer=True,positive=True)
    check('raw_angular_tail_estimator_diverges',S.simplify(S.summation(2*(q+2*n),(n,1,L))-(2*q*L+2*L*(L+1)))==0)
    check('summable_mass_power_three',S.Rational(12,18)==S.Rational(2,3))
    factor=arb('1.478e-44')*2*C/3
    budgets={'near_shape':'6.4617454e-12','other_affine':'5.34234999994958e-12','best_case_feedback':'.0012069'}
    ceilings={name:enc(arb(budget)/factor) for name,budget in budgets.items()}
    oldmass=prior('current-complete-weak-mass-residual')
    far=prior('current-actual-far-graph-energy-bound')
    check('retained_linear_mass_and_far',dec(oldmass['full_coupled_linear_mass_upper'])<arb('1.627532e-15') and dec(far['actual_linear_far_graph_upper'])<arb('1.24e-20'))
    check('original_mass_conversion_reserve',arb('1e-11')-2174*arb('1.627532e-15')-arb('1.24e-20')>arb('6.4617454e-12'))
    check('conditional_fast_scope_arithmetic',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    return {'precision_bits':ctx.prec,'actual_Q_second_derivative_upper':enc(q2guard),'actual_C':enc(C),
      'all_mode_maximum_location':enc(ymax),'coupling_upper_expression':enc(upper),
      'family_ratio_upper_expression':enc(qreal),'complex_ratio_upper_expression':enc(qcomplex),
      'family_inverse_bound':enc(1/(1-arb('.38601'))),'inverse_variation_tail_coefficients':tails,
      'conditional_summable_tail_factor':enc(factor),'conditional_original_output_B_ceilings':ceilings,
      'scope':{'actual_all_mode_dirac_transport_contraction':True,'regular_bolts_and_paired_transmission_retained':True,
        'actual_paired_mass_field_inverse_variations':True,'complete_renormalized_hessian_tail':False,
        'raw_uniform_trace_estimator_failed':True,'actual_kappa_m_computed':False,
        'full_coupled_signed_center_computed':False,'complete_original_graph_bound':False,
        'nonlinear_admission':False,'exact_compatibility':False,'physical_matching':False,
        'validity_noise':False,'all_sector_stability':False,'actual_controlled_solution':False,
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
        run['result']=execute()
        run['input_sha256_after']={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after'])
        run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),
       'coupling':run.get('result',{}).get('coupling_upper_expression'),
       'conditional_B_ceilings':run.get('result',{}).get('conditional_original_output_B_ceilings'),
       'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
