"""Supporting arithmetic for the actual original outer graph remainder.

No signed center, interior inverse, or post-elimination mass column.
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
ROOT=Path(__file__).resolve().parents[1];STEM='current-original-graph-dirac-tail'
OWN=['notes/'+STEM+'.md','checks/current_original_graph_dirac_tail.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/current-adaptive-dirac-trace-remainder.md','evidence/current-adaptive-dirac-trace-remainder.json',
 'notes/current-global-dirac-transport-series.md','evidence/current-global-dirac-transport-series.json',
 'notes/current-global-mass-conversion.md','evidence/current-global-mass-conversion.json',
 'notes/current-global-far-return-energy-bound.md','notes/current-coupled-global-return-attempt.md',
 'notes/current-mass-response-regulator-audit.md','notes/current-third-metric-conversion.md',
 'notes/current-complete-strong-response-attempt.md','notes/current-short-travel-strong-domain.md',
 'notes/current-rough-paired-variation.md','notes/current-coincident-dirac-action-assembly.md']
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
    check('ordinary_argument_frozen',sha256((ROOT/OWN[0]).read_bytes()).hexdigest()==review['scientific_sha256'][OWN[0]])
    p=prior('current-adaptive-dirac-trace-remainder')['actual_parameters']
    C,B2,gamma=[dec(p[x]) for x in ['C','B2','gamma']]
    par=prior('current-global-mass-conversion')['parameters']
    beta,u,eps,tau0=[dec(par[x]) for x in ['beta','u','eps','tau']]
    co={0:beta,2:u,3:-4*beta-2*u-tau0,4:3*beta+u+eps,5:-2*eps,6:eps}
    radius=arb('.101')
    remainder=sum((abs(v)*radius**j for j,v in co.items() if j),arb(0))
    qp=sum((j*abs(v)*radius**(j-1) for j,v in co.items() if j),arb(0))
    check('one_sided_complex_seed_Q_zero_free',beta-remainder>arb('.16'))
    check('one_sided_complex_seed_Q_upper',beta+remainder<arb('.18'))
    check('one_sided_complex_seed_Qprime',qp<arb('.1'))
    check('all_modes_have_M_above_eighteen',(arb(2000)**4+1)/C>18**2)
    check('retained_complete_L2_coupling',dec(prior('current-global-dirac-transport-series')['coupling_upper_expression'])<arb('.385'))
    e=arb('.011');q=arb('.401')
    kernel_ratio=arb('.385')+arb(2).sqrt()*(e.exp()-1)
    check('whole_interval_Linfinity_kernel_contraction',kernel_ratio<q)
    ej=[e.exp()-1,e.exp()*e,e.exp()*(e+e**2),e.exp()*(e+3*e**2+e**3)]
    check('complete_exponential_three_jets',all(x<arb('.02') for x in ej))
    check('analytic_seed_coefficient_guards',arb('.18').sqrt()/arb('.16')<8 and 1/arb('.16')<8 and 2/arb('.16').sqrt()<10 and arb('.1')/(8*arb('.16'))<1)
    # Exact products: the deliberately loose 32 d j! 2000^j covers
    # seed Cauchy jets and every derivative of the complex mass vertex.
    for r in range(4):
        V_over_d=arb(factorial(r)*1000**r)/18
        V_over_d+=sum((arb(comb(r,j))*ej[j]*10*factorial(r-j)*1000**(r-j) for j in range(r+1)),arb(0))
        check('actual_potential_jet_guard_'+str(r),V_over_d<32*factorial(r)*2000**r)
    for r in range(1,5):
        lhs=64*sum(2000**i*4096**(r-1-i) for i in range(r))
        check('kernel_derivative_induction_'+str(r),lhs<=r*4096**r)
    check('diagonal_product_jet_guard',8192+arb(1000)/18<10000)
    check('new_rule_preserves_previous_weak_rule',8/(1/q).log()>4/(1/arb('.3961')).log() and 6/(1/q).log()>2/(1/arb('.3961')).log())
    M,pf,lam=S.symbols('M pf lam',positive=True)
    expr=(8*S.log(M/18)+6*S.log(1+pf))/lam
    check('complete_spectral_power_identity',S.simplify(S.exp(-lam*expr)-(18/M)**8/(1+pf)**6)==0)
    check('complete_derivative_moment_reduction',S.simplify(M*(M*(1+pf))**4*(18/M)**8/(1+pf)**6-18**8/M**3/(1+pf)**2)==0)
    n,qf=S.symbols('n qf',positive=True)
    count=2*qf*n+2*n*(n+1)
    check('exact_angular_count_identity',S.expand(4*n*(n+qf)-count-2*n*(n+qf-1))==0)
    x=S.symbols('x',positive=True)
    check('entire_thermal_integral',S.integrate((1+x)**-2,(x,0,S.oo))==1)
    tau=480*arb(10000)**4*18**8*q**124/(1-q)*4*C*(1+gamma)/(3*arb.pi()*B2)
    check('whole_spectrum_C4_source_remainder',tau<arb('1.354e-53'))
    a,b,rho=arb('.1'),arb('.05'),arb('-.005')
    K=[arb(1)]+[arb(factorial(21)//factorial(21-j))*2**(j-1)/arb(2).log()**j for j in range(1,5)]
    bern=[1]*11+[0]*11
    for j in range(1,5):
        bern=[bern[k+1]-bern[k] for k in range(len(bern)-1)]
        check('original_beta11_cutoff_jet_'+str(j),max(map(abs,bern))<=2**(j-1))
    def rec(initial):
        rows=[initial]
        for r in range(3):
            out={}
            for j,v in rows[-1].items():
                out[j]=out.get(j,arb(0))+(rho+j)*v
                out[j+1]=out.get(j+1,arb(0))+v
            rows.append(out)
        return rows
    cc,dd=rec({0:arb(1)}),rec({1:arb(1)})
    def weight(t,tr):
        if tr:return ((a**(2*t)-b**(2*t))/(2*t)).sqrt() if not t.is_zero() else (a/b).log().sqrt()
        return a**t/(2*t).sqrt()
    L=[]
    for ell in range(4):
        val=arb(0)
        for j in range(ell+1):
            c=sum((abs(v)*weight(rho+r,True) for r,v in cc[ell-j].items()),arb(0))*K[j+1]
            d=sum((abs(v)*weight(rho+r,j>=1) for r,v in dd[ell-j].items()),arb(0))*K[j]
            val+=comb(ell,j)*(c+d)
        L.append(val)
    Cchi=sum((arb(comb(3,j)).sqrt()*L[j] for j in range(4)),arb(0))
    check('whole_original_cutoff_graph_constant',Cchi<arb('4368755'))
    check('unbounded_fourth_input_jet_has_unit_coefficient',dd[3][4].is_exact() and dd[3][4]==1)
    graph=116*Cchi*tau
    check('whole_original_H3_graph_J3_coefficient',graph<arb('6.861e-45'))
    quartic=116*tau;mass=arb('.1')*quartic
    emax=arb('1e-30')
    slow_multiplier=b**(-emax)+a**(1-emax)*(1+1/(1-emax))
    check('exact_slow_root_whole_interval_guard',slow_multiplier<arb('1.201'))
    slow=slow_multiplier*quartic
    check('original_quartic_tail_guard',quartic<arb('1.571e-51'))
    check('original_mass_moment_tail_guard',mass<arb('1.571e-52'))
    check('original_slow_moment_tail_guard',slow<arb('1.896e-51'))
    z=S.symbols('z',positive=True);Q=S.Function('Q')(z);J=S.Function('J')(z)
    chi=S.Function('chi')(z);f=S.Function('f')(z);betaS,CS=S.symbols('beta C',positive=True)
    D=lambda v:z*S.diff(v,z)
    nm=-2*betaS*Q*J/z;bm=betaS*z**3/Q+betaS*S.diff(Q,z)*J;F=CS*Q/z**2
    lift=lambda v:S.simplify(v.subs(S.diff(J,z),z**3/Q**2))
    areal=lambda nn,bb:2*F*(-D(nn)-bb)+D(F)*nn
    lapse=lambda nn,bb:nn+2*bb+D(nn)
    check('both_reference_local_rows_cancel',lift(areal(nm,bm)-2*CS*betaS*z)==0 and lift(lapse(nm,bm))==0 and S.diff(2*CS*betaS*z/z,z)==0)
    check('full_reference_mass_transition_retained',lift(areal(chi*nm,chi*bm)-chi*2*CS*betaS*z+2*F*nm*D(chi))==0 and lift(lapse(chi*nm,chi*bm)-nm*D(chi))==0)
    check('original_graph_transition_identity',S.simplify((D(chi*z**4*f)-4*chi*z**4*f)/z**4-(D(chi)*f+chi*z*S.diff(f,z)))==0)
    check('original_mass_moment_integrand',S.simplify(z*(S.diff(chi,z)*f+chi*S.diff(f,z))-(S.diff(z*chi*f,z)-chi*f))==0)
    # Noncommuting Schur identity uses only P-Ptilde=-P deltaD Ptilde.
    B,Cx,PP,Pt,dtB,dtC,dtD=S.symbols('B C P Pt dB dC dD',commutative=False)
    raw=B*PP*Cx-(B-dtB)*Pt*(Cx-dtC)
    expected=dtB*PP*Cx+(B-dtB)*PP*dtC-(B-dtB)*PP*dtD*Pt*(Cx-dtC)
    check('exact_Schur_error_identity',S.expand(raw-expected-(B-dtB)*(PP-Pt+PP*dtD*Pt)*(Cx-dtC))==0)
    budgets={'near_shape':arb('6.4617454e-12'),'other_affine':arb('5.34234999994958e-12'),
             'feedback_best_case':arb('.0012069'),'mass_column_screen':arb('5e-5')}
    ceilings={k:{'A_times_U3_half_split':enc(v/(2*arb('6.861e-45'))),
                 'A_times_W4_half_split':enc(v/(2*arb('1.354e-53')))} for k,v in budgets.items()}
    check('conditional_fast_screen_retained',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    return {'precision_bits':ctx.prec,'complex_collar_Q_remainder':enc(remainder),'complex_collar_Qprime_upper':enc(qp),
      'global_kernel_ratio_upper_expression':enc(kernel_ratio),
      'spectral_order_rule':'128+ceil((8log(M/18)+6log(1+abs(p)))/log(1/.401))',
      'C4_source_tail_coefficient':enc(tau),'original_cutoff_constant':enc(Cchi),
      'original_H3_tail_U3_coefficient':enc(graph),'original_H3_tail_W4_coefficient':enc(tau),
      'original_quartic_tail_U3_coefficient':enc(quartic),'original_mass_moment_tail_U3_coefficient':enc(mass),
      'original_slow_moment_tail_U3_coefficient':enc(slow),
      'conditional_product_ceilings_not_actual_constants':ceilings,
      'scope':{'actual_paired_mass_log_domain':True,'complete_angular_thermal_sums':True,
       'actual_bolt_and_finite_returns_in_kernel':True,'original_outer_H3_tail_certified':True,
       'original_outer_phase_observation_tails_certified':True,'exp_mass_tadpole_retained':True,
       'common_regulator_and_full_conversion':'previous weak trace identity plus absolutely summable differentiated remainder; conversion and all local restorations remain in retained center',
       'actual_U3_W4_input_maps_computed':False,'signed_center_computed':False,
       'interior_strong_inverse_computed':False,'complete_Schur_return_certified':False,
       'actual_kappa_m_computed':False,'nonlinear_admission':False,'nonlinear_contraction':False,
       'exact_compatibility':False,'physical_matching':False,'semiclassical_validity_noise':False,
       'all_sector_stability':False,'construction_obstruction':False,'actual_controlled_solution':False,'full_goal_complete':False}}

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
    message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),
       'graph_U3':run.get('result',{}).get('original_H3_tail_U3_coefficient'),
       'phase_U3':run.get('result',{}).get('original_slow_moment_tail_U3_coefficient'),
       'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
