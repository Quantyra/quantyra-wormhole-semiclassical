"""Complete separated return guards; not a physical coupled solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
import json, sys, traceback
import sympy as sp
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb, ctx
ctx.prec=512
ROOT=Path(__file__).resolve().parents[1]
STEM='current-global-far-return-energy-bound'
OWN=['notes/'+STEM+'.md','checks/current_global_far_return_energy_bound.py',
     'evidence/'+STEM+'-review.json']
sys.path.insert(0,str(ROOT/'checks'))
saved_argv=sys.argv[:]
sys.argv=[sys.argv[0]]
import current_physical_york_dual_assembly as york
sys.argv=saved_argv
CHECKS=[]

def check(name,value,detail=None):
    CHECKS.append(dict(name=name,passed=bool(value)))
    if not value: raise AssertionError(name)

york.check=check
york.prior.check=check

def enc(v):
    m,r,e=v.mid_rad_10exp(60)
    return dict(ball=v.str(40),decimal_enclosure=dict(mid=str(m),rad=str(r),exponent=str(e)))

def execute():
    p,h,zh,env=york.actual_parameters()
    coeff={k:env[v] for k,v in zip(york.powers,york.qs)}
    C,B2,gamma,Fl=p['C'],p['B2'],p['gamma'],p['Fl']
    beta,u,eps,tau=coeff[0],coeff[2],coeff[6],p['tau']
    a,b,c=map(arb,('.1','.05','.2'))
    def bernstein(co):
        n=max(co)
        return [sum((co.get(i,arb(0))*arb(comb(k,i))/comb(n,i) for i in range(k+1)),arb(0)) for k in range(n+1)]
    q2={k-2:k*(k-1)*v for k,v in coeff.items() if k>=2}
    check('global_Q_upper_on_entire_seed',all(v<arb('.17') for v in bernstein(coeff)))
    check('global_Q_second_upper',all(v<arb('2.1') for v in bernstein(q2)))
    check('Qprime_strictly_negative_on_c_to_bolt',12*beta-2*u/c-4*abs(eps)>0 and tau>0)
    qc=sum((v*c**k for k,v in coeff.items()),arb(0))
    check('interface_f_lower',qc>arb('.16'))
    M1=((p['q']+1)/C).sqrt()
    check('all_angular_mass_lower',M1>18)
    complex_qlo=beta-sum((abs(v)*arb('.101')**k for k,v in coeff.items() if k),arb(0))
    check('complex_one_sided_Q_lower',complex_qlo>arb('.16'))
    check('complex_output_disks_stay_inside_verified_region',a+1/(400*18)<arb('.101'))
    check('weighted_bolt_damping_guard',arb(3).sqrt()-arb('1.05').sqrt()/36>arb('1.7'))
    check('endpoint_L2_sensitivity_at_most_M',20/(17*arb('.4'))<18)
    check('collar_rate_3M_plus_8p',9*arb('.17')/4+64*arb('.17')**2/4<1)
    check('separation_exceeds_short_travel_support',(c-a)/arb('.17').sqrt()>arb(1)/64)
    # Exact minimization of sqrt(4M^2+t^2)-t/2.
    M,t=sp.symbols('M t',positive=True)
    damp=sp.sqrt(4*M*M+t*t)-t/2
    check('exact_weighted_damping_minimum',sp.simplify(damp.subs(t,2*M/sp.sqrt(3))-sp.sqrt(3)*M)==0 and sp.simplify(sp.diff(damp,t).subs(t,2*M/sp.sqrt(3)))==0 and sp.simplify(sp.diff(damp,t,2)-4*M*M/(4*M*M+t*t)**sp.Rational(3,2))==0)
    check('complex_flow_and_sensitivity_guards',arb('5.5')<20 and arb('13.8125')<20 and arb('13.75')<20 and arb(20)/400<arb('.1'))
    sourceguard=2*arb('2.5')*2*(2*arb('1.1')*(arb(20)/400).exp())/arb('.79')**2
    check('full_original_density_two_arm_complex_guard',sourceguard<64)
    U,V=sp.symbols('U V')
    ratio=(1-U*V)/(1+U*V)
    check('exact_coincident_spin_sensitivity',sp.simplify(sp.diff(ratio,U)+2*V/(1+U*V)**2)==0 and sp.simplify(sp.diff(ratio,V)+2*U/(1+U*V)**2)==0)
    n,q=sp.symbols('n q',positive=True)
    N=2*q*n+2*n*(n+1)
    check('exact_full_angular_count_envelope',sp.simplify(4*n*(n+q)-N-2*n*(n+q-1))==0)
    # G=1 is the unchanged length unit. Restore physical factors symbolically.
    G,B,lam,Cs,gam,NH,F=sp.symbols('G B lam C gam NH F',positive=True)
    period=2*sp.pi*B/(Cs*gam)
    j0=1000*G*Cs**2/(lam*B**4)
    L2perE=j0/Cs*(NH*4*sp.pi*period)/(NH*2*sp.pi/gam)/(sp.Rational(99,100)*F)*(lam*B**3/(8*sp.pi*G))
    check('exact_original_energy_to_L2_conversion',sp.simplify(L2perE-500/(sp.Rational(99,100)*F))==0)
    physical_mode=-2*M/sp.Symbol('f',positive=True)
    normalized=physical_mode/(-sp.pi*B**2/(G*gam))
    check('original_mode_source_sign_and_measure',sp.simplify(normalized-2*G*gam*M/(sp.pi*B**2*sp.Symbol('f',positive=True)))==0)
    ct=(arb(500)/(arb('.99')*Fl)).sqrt()
    check('energy_to_L2_guard',ct<arb('22.474'))
    S=[]
    for r in range(5):
        Z=sum((arb(comb(r,l)*factorial(4+r-l)*factorial(l))/(arb('.3')**(4+r-l)*arb('.8')**(l+1)) for l in range(r+1)),arb(0))
        S.append(512*C/(arb.pi()*B2)*(1+gamma*arb('.8'))*factorial(r)*400**r*Z)
    # Independent symbolic factorial integral identities behind every term.
    x,k=sp.symbols('x k',positive=True)
    for power in range(9):
        check('complete_positive_moment_integral_'+str(power),sp.integrate(x**power*sp.exp(-k*x),(x,0,sp.oo))==sp.factorial(power)/k**(power+1))
    z,rr=sp.symbols('z rr',positive=True)
    chi=sp.Function('chi')(z);R=sp.Function('R')(z)
    D=lambda f:z*sp.diff(f,z)
    Q=sp.Function('Q')(z);J=sp.Function('J')(z);bet=sp.symbols('beta',positive=True)
    nm=-2*bet*Q*J/z;bm=bet*z**3/Q+bet*sp.diff(Q,z)*J;Fref=Cs*Q/z**2
    lift=lambda expr:sp.simplify(expr.subs(sp.diff(J,z),z**3/Q**2))
    ff=lambda nn,bb:2*Fref*(-D(nn)-bb)+D(Fref)*nn
    al=lambda nn,bb:nn+2*bb+D(nn)
    check('exact_reference_mass_both_local_constraints',lift(ff(nm,bm)-2*Cs*bet*z)==0 and lift(al(nm,bm))==0 and sp.diff(2*Cs*bet*z/z,z)==0)
    check('exact_reference_cutoff_transition_retained',lift(ff(chi*nm,chi*bm)-chi*2*Cs*bet*z+2*Fref*nm*D(chi))==0 and lift(al(chi*nm,chi*bm)-nm*D(chi))==0)
    check('original_transition_term_retained',sp.simplify((D(chi*z**4*R)-4*chi*z**4*R)/z**4-(D(chi)*R+chi*z*sp.diff(R,z)))==0)
    check('mass_moment_exact_integrand',sp.simplify(z*(sp.diff(chi,z)*R+chi*sp.diff(R,z))-(sp.diff(z*chi*R,z)-chi*R))==0)
    t=sp.symbols('t')
    cutoff=sum(sp.binomial(21,k)*t**k*(1-t)**(21-k) for k in range(11))
    check('actual_beta11_cutoff_endpoints',cutoff.subs(t,0)==1 and cutoff.subs(t,1)==0)
    check('actual_beta11_cutoff_monotonicity',sp.expand(sp.diff(cutoff,t)+sp.factorial(21)/sp.factorial(10)**2*t**10*(1-t)**10)==0)
    check('actual_cutoff_first_ten_endpoint_jets',all(sp.diff(cutoff,t,k).subs(t,0)==0 and sp.diff(cutoff,t,k).subs(t,1)==0 for k in range(1,11)))
    bc=[1]*11+[0]*11
    for j in range(1,5):
        bc=[bc[k+1]-bc[k] for k in range(len(bc)-1)]
        check('Bernstein_cutoff_derivative_guard_'+str(j),max(map(abs,bc))<=2**(j-1))
    K=[arb(1)]+[arb(factorial(21)//factorial(21-j))*2**(j-1)/arb(2).log()**j for j in range(1,5)]
    rho=sp.Rational(-1,200)
    def recurrence(initial):
        rows=[initial]
        for m in range(3):
            prev=rows[-1]
            rows.append({r:sp.expand((rho+r)*prev.get(r,0)+prev.get(r-1,0)) for r in range(max(prev)+2)})
        return rows
    cc,dd=recurrence({0:sp.Integer(1)}),recurrence({1:sp.Integer(1)})
    for name,rows,expr in [('c',cc,R),('d',dd,z*sp.diff(R,z))]:
        for m,row in enumerate(rows):
            expected=sum(v*z**r*sp.diff(R,z,r) for r,v in row.items())
            check('exact_graph_jet_'+name+str(m),sp.expand(expr-expected)==0)
            expr=sp.expand(D(expr)+rho*expr)
    rhoA=arb('-.005')
    Wtr=lambda t:((a**(2*t)-b**(2*t))/(2*t)).sqrt()
    Wout=lambda t:a**t/(2*t).sqrt()
    aa=lambda x:abs(arb(int(sp.numer(x)))/int(sp.denom(x)))
    L=[]
    for ell in range(4):
        accum=arb(0)
        for j in range(ell+1):
            first=sum((aa(v)*S[r]*Wtr(rhoA+r) for r,v in cc[ell-j].items()),arb(0))*K[j+1]
            second=sum((aa(v)*S[r]*(Wout(rhoA+r) if j==0 else Wtr(rhoA+r)) for r,v in dd[ell-j].items() if r>=1),arb(0))*K[j]
            accum+=comb(ell,j)*(first+second)
        L.append(ct*accum)
    far=sum((arb(comb(3,l)).sqrt()*L[l] for l in range(4)),arb(0))
    e=arb('1e-30')
    slow=ct*(b**(-e)*S[0]+a**(1-e)*(1+1/(1-e))*S[1])
    vals={'graph':far,'quartic':ct*S[0],'mass_moment':a*ct*S[0],'slow_moment':slow}
    guards={'graph':'4.583e-13','quartic':'7.949e-26','mass_moment':'7.949e-27','slow_moment':'1.141e-22'}
    affine={'graph':'5.042e-23','quartic':'8.744e-36','mass_moment':'8.744e-37','slow_moment':'1.256e-32'}
    for name,val in vals.items():
        check('complete_original_'+name+'_guard',val<arb(guards[name]))
        check('affine_'+name+'_guard',arb(guards[name])*arb('1.1e-10')<arb(affine[name]))
    remainder=arb('5.34235e-12')-arb(affine['graph'])
    gamma_only_far=arb('.0012069')/arb(guards['graph'])
    gamma_with_mass=arb('.0012069')/(arb('.0195239')+arb(guards['graph']))
    check('remaining_affine_budget_positive',remainder>arb('5.34234999994e-12'))
    check('conservative_mass_plus_far_energy_comparison_target',arb('.06181654')<gamma_with_mass<arb('.06181655'))
    check('conditional_fast_screen_scope_unchanged',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    # Positive distance^-5 estimator cannot reach coincidence in L1.
    epsilon,d=sp.symbols('epsilon d',positive=True)
    primitive=sp.integrate(d**-5,d)
    check('absolute_near_estimator_diverges',sp.limit(primitive.subs(d,1)-primitive.subs(d,epsilon),epsilon,0,dir='+')==sp.oo)
    return dict(outcome='unresolved_complete_coupled_estimate_with_certified_actual_far_entry',
      original_energy_to_L2=enc(ct),source_derivative_upper_expressions=[enc(v) for v in S],
      graph_derivative_upper_expressions=[enc(v) for v in L],
      original_output_upper_expressions={k:enc(v) for k,v in vals.items()},
      upward_original_output_guards=guards,upward_affine_guards=affine,
      conditional_post_mass_and_far_affine_remaining=enc(remainder),
      conditional_energy_comparison_allowance_far_alone=enc(gamma_only_far),
      conditional_energy_comparison_allowance_mass_and_far_only=enc(gamma_with_mass),
      actual_energy_comparison=None,actual_mass=None,kappa_m=None,M_graph_4=None,
      actual_physical_operator_applications=0,coupled_signed_solve_executed=False,
      certified_operator_scope='all-mode fixed-seed input z>=.2 to original output z<=.1; regular bolt and both paired paths',
      near_absolute_estimator='nonintegrable; common subtraction required, not a physical obstruction',
      full_coupled_feasibility_decided=False,nonlinear_requirement_closed=False,
      full_objective_preserved=True,closed_audit_preserved=True)

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    out=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    review=json.loads((ROOT/OWN[2]).read_bytes())
    paths=set(OWN+review['retained_inputs'])
    for module in list(sys.modules.values()):
        f=getattr(module,'__file__',None)
        if f:
            try:rel=Path(f).resolve().relative_to(ROOT).as_posix()
            except ValueError:continue
            if rel.endswith('.py'):paths.add(rel)
    digest=lambda p:sha256(p.read_bytes()).hexdigest()
    before={p:digest(ROOT/p) for p in sorted(paths)}
    fr=Path(flint.__file__).resolve().parent
    runtime=dict(python=sys.version,python_executable_sha256=digest(Path(sys.executable)),flint_version=flint.__version__,
      flint_artifacts={p.relative_to(fr).as_posix():digest(p) for p in sorted(fr.rglob('*')) if p.is_file() and p.suffix in ['.py','.pyd','.dll']})
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),mode=mode,input_sha256_before=before,precision_bits=512,runtime=runtime,
      frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN})
    try:
        check('ordinary_argument_frozen',before[OWN[0]]==review['ordinary_argument_sha256'])
        run['result']=execute()
        run['input_sha256_after']={p:digest(ROOT/p) for p in before}
        check('frozen_inputs_unchanged',before==run['input_sha256_after'])
        run['status']='passed'
    except Exception:run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    ledger['runs'].append(run);out.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write((json.dumps(run)+'\n').encode())
    print(json.dumps(dict(mode=mode,status=run['status'],checks=len(CHECKS),result=run.get('result'),traceback=run.get('traceback')),indent=2))
    if run['status']!='passed':raise SystemExit(1)
