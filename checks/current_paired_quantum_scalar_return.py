"""Support the ordinary all-mode scalar response proof; no coupled solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb
import json, sys, traceback
import sympy as sp
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb, ctx
ctx.prec=512
ROOT=Path(__file__).resolve().parents[1]
STEM='current-paired-quantum-scalar-return'
OWN=['notes/'+STEM+'.md','checks/current_paired_quantum_scalar_return.py',
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
    def bernstein(co):
        n=max(co)
        return [sum((co.get(i,arb(0))*arb(comb(k,i))/comb(n,i) for i in range(k+1)),arb(0)) for k in range(n+1)]
    q2={k-2:k*(k-1)*v for k,v in coeff.items() if k>=2}
    c=arb('.2')
    beta,u,eps,tau=coeff[0],coeff[2],coeff[6],p['tau']
    check('actual_seed_global_Q_upper',all(v<arb('.17') for v in bernstein(coeff)))
    check('actual_seed_global_Qsecond_upper',all(v<arb('2.1') for v in bernstein(q2)))
    check('Qprime_negative_on_c_to_bolt',12*beta-2*u/c-4*abs(eps)>0 and tau>0)
    nearq=beta-sum((abs(v)*c**k for k,v in coeff.items() if k),arb(0))
    nearqp=sum((k*abs(v)*c**(k-1) for k,v in coeff.items() if k),arb(0))
    check('whole_near_Q_lower',nearq>arb('.16'))
    check('whole_near_Qprime_guard',nearqp<arb('.3'))
    check('whole_near_log_f_guard',arb('.3')/(2*arb('.4'))<arb('1.05').sqrt())
    U,V,W,M=sp.symbols('U V W M',positive=True)
    t=(1-U*V)/(1+U*V)
    a=(U+V)/(1+U*V)
    b=(U-V)/(1+U*V)
    Up=W*(1-U**2)-2*M*U
    Vp=-W*(1-V**2)+2*M*V
    deriv=lambda expr:sp.diff(expr,U)*Up+sp.diff(expr,V)*Vp
    check('exact_trace_equation',sp.simplify(deriv(t)-2*W*b)==0)
    check('exact_a_equation',sp.simplify(deriv(a)+2*M*b)==0)
    check('exact_b_equation',sp.simplify(deriv(b)-2*W*t+2*M*a)==0)
    check('exact_quadratic_constraint',sp.simplify(t*t+a*a-b*b-1)==0)
    A,B,T,y,yp,wp,v=sp.symbols('a b t y yp wp v')
    da=(B*yp/(2*W)-T*y)/A
    ypp=wp/W*yp+2*W*(2*W*y-2*M*da-2*M*v*A)
    expected=(wp/W-2*M*B/A)*yp+(4*W**2+4*M*W*T/A)*y-4*M*W*A*v
    check('exact_linearized_elimination',sp.simplify(ypp-expected)==0)
    check('exact_divergence_form_rhs',sp.simplify((-(ypp-(wp/W-2*M*B/A)*yp)+(4*W**2+4*M*W*T/A)*y)/(W*A)-4*M*v)==0)
    ell=sp.symbols('ell',positive=True)
    block=sp.Matrix([[-M*T,-M],[-M,ell/4]])
    check('exact_auxiliary_Schur_complement',sp.simplify(block[0,0]-block[0,1]*block[1,0]/block[1,1]+M*T+4*M**2/ell)==0)
    check('actual_paired_seam',sp.simplify(b.subs(V,U))==0)
    rho,omega=sp.symbols('rho omega',positive=True)
    leading=lambda f:-sp.diff(rho/omega*sp.diff(f,rho),rho)+4*omega/rho*f
    check('regular_bolt_indicial_roots',sp.simplify(leading(rho**(2*omega)))==0 and sp.simplify(leading(rho**(-2*omega)))==0)
    check('actual_bolt_nonresonant_coefficient',sp.simplify(leading(rho)-(4*omega**2-1)/omega)==0)
    check('lowest_thermal_resonant_log_coefficient',sp.simplify(leading(rho*sp.log(rho)).subs(omega,sp.Rational(1,2))+4)==0)
    check('log_branch_finite_form_integral',sp.integrate(rho*sp.log(rho)**2,(rho,0,1))==sp.Rational(1,4))
    check('singular_lowest_branch_not_L2',sp.integrate(rho**-2,(rho,0,1))==sp.oo)
    d,sigma,r=sp.symbols('d sigma r',positive=True)
    lower_derivative=-(2*W+sigma)*W/d
    F=W*(1-(W/d)**2)-2*M*W/d
    check('global_upper_Riccati_cone',sp.simplify((lower_derivative-F).subs(M,d+sigma/2)-W*(W/d-1)**2)==0)
    vlo=4*d*r/(2*r/(1+r*r))+4*M*((1-r*r)/(1+r*r))/(2*r/(1+r*r))**2
    check('global_scalar_potential_lower_expression',sp.simplify(vlo-(M/r**2+2*d+(2*d-M)*r**2))==0)
    check('global_scalar_potential_gap',sp.factor(vlo-4*d-(1-r*r)*(M/r**2-(2*d-M)))==0)
    aDlo=4*d*d*r*r+4*M*d*r*((1-r*r)/(1+r*r))/(2*r/(1+r*r))
    check('nondivergence_potential_lower_expression',sp.simplify(aDlo-(2*M*d+(4-2*M/d)*d*d*r*r))==0)
    x,Vnorm=sp.symbols('x Vnorm',positive=True)
    check('derivative_energy_completion_of_square',sp.expand(M*M/d*Vnorm**2-(4*M*Vnorm*x-4*d*x*x)-4*d*(x-M*Vnorm/(2*d))**2)==0)
    M1=((p['q']+1)/p['C']).sqrt()
    sg=arb('1.05').sqrt()
    ratio=M1/(M1-sg/2)
    inv=1/(4*(M1-sg/2))
    r18=1-sg/36
    check('all_angular_mass_lower',M1>18)
    check('uniform_scalar_inverse_guard',inv<arb('.014256564'))
    check('uniform_response_guard',ratio<arb('1.029218'))
    check('nondivergence_low_W_mass_guard',2*r18>arb('1.94'))
    check('nondivergence_low_W_frequency_guard',4-2/r18>arb('1.94'))
    check('nondivergence_high_W_guard',4*r18*r18/(1+r18*r18)>arb('1.94'))
    G,Bphysical,lam,C,gam,NH,Fl=sp.symbols('G B lam C gam NH Fl',positive=True)
    period=2*sp.pi*Bphysical/(C*gam)
    j0=1000*G*C*C/(lam*Bphysical**4)
    L2perE=j0/C*(NH*4*sp.pi*period)/(NH*2*sp.pi/gam)/(sp.Rational(99,100)*Fl)*(lam*Bphysical**3/(8*sp.pi*G))
    check('original_energy_to_canonical_L2',sp.simplify(L2perE-500/(sp.Rational(99,100)*Fl))==0)
    ct=(500/(arb('.99')*p['Fl'])).sqrt()
    nearin=arb('2.5').sqrt()*ct
    ret=ratio*nearin
    outer=arb('.17')**arb('.25')*ret
    check('original_energy_conversion_guard',ct<arb('22.474'))
    check('near_input_proper_L2_guard',nearin<arb('35.534'))
    check('global_auxiliary_return_guard',arb('1.029218')*arb('35.534')<arb('36.573'))
    check('outer_canonical_L2_guard',arb('.17')**arb('.25')*arb('36.573')<arb('23.485'))
    check('affine_auxiliary_per_mode_guard',arb('36.573')*arb('1.1e-10')<arb('4.024e-9'))
    dn,f=sp.symbols('dn f',positive=True)
    orig=-2*dn*M/f*(T*v+y)
    normalized=orig/(-sp.pi*Bphysical**2/(G*gam))
    check('original_source_sign_and_normalization',sp.simplify(normalized-2*G*gam*dn*M/(sp.pi*Bphysical**2*f)*(T*v+y))==0)
    check('outer_raw_mode_upper_coefficient',2*(arb('22.474')+arb('23.485'))/arb('.4')<arb('229.796') and 2*(sp.Rational('22.474')+sp.Rational('23.485'))/sp.Rational('.4')==sp.Rational('229.795'))
    J=sp.symbols('J',integer=True,nonnegative=True)
    check('complete_raw_thermal_norm_estimator_diverges',sp.limit((J+1)*M,J,sp.oo)==sp.oo)
    v0=sp.symbols('v0')
    logdensity=normalized.subs({T:-M*rho*sp.log(rho),y:-M*v0*rho*sp.log(rho),v:v0,f:gam*rho})
    check('actual_lowest_bolt_source_log_with_tadpole',sp.simplify(logdensity+4*G*dn*M*M*v0*sp.log(rho)/(sp.pi*Bphysical**2))==0)
    check('lowest_bolt_unregulated_source_not_bounded',sp.limit(-sp.log(rho),rho,0,dir='+')==sp.oo)
    z=sp.symbols('z',positive=True)
    Q=sp.Function('Q')(z);JJ=sp.Function('J')(z);chi=sp.Function('chi')(z);bet=sp.symbols('beta',positive=True)
    D=lambda e:z*sp.diff(e,z)
    nm=-2*bet*Q*JJ/z;bm=bet*z**3/Q+bet*sp.diff(Q,z)*JJ;Fref=C*Q/z**2
    lift=lambda e:sp.simplify(e.subs(sp.diff(JJ,z),z**3/Q**2))
    ff=lambda nn,bb:2*Fref*(-D(nn)-bb)+D(Fref)*nn
    al=lambda nn,bb:nn+2*bb+D(nn)
    check('exact_reference_mass_both_constraints',lift(ff(nm,bm)-2*C*bet*z)==0 and lift(al(nm,bm))==0 and sp.diff(2*C*bet*z/z,z)==0)
    check('actual_mass_lift_cutoff_transition',lift(ff(chi*nm,chi*bm)-2*C*bet*z*chi+2*Fref*nm*D(chi))==0 and lift(al(chi*nm,chi*bm)-nm*D(chi))==0)
    km=arb('5e-5')
    check('conditional_mass_affine_debit',km*arb('4.2953e-8')<arb('2.14765000000001e-12'))
    check('conditional_mass_graph_debit',km*arb('49.462')<arb('.00247310000000001'))
    check('conditional_mass_energy_debit',km*arb('390.478')<arb('.01952390000000001'))
    remainder=arb('5.34235e-12')-arb('5.042e-23')
    gammaE=arb('.0012069')/(arb('.0195239')+arb('4.583e-13'))
    check('conditional_remaining_affine_budget',remainder>arb('5.34234999994957e-12'))
    check('conditional_remaining_energy_comparison',arb('.06181654279993')<gammaE<arb('.06181654279994'))
    check('conditional_fast_screen_unchanged',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    return dict(outcome='unresolved_complete_coupled_return_with_exact_uniform_auxiliary_response',
      actual_M1=enc(M1),uniform_inverse_upper_expression=enc(inv),uniform_response_upper_expression=enc(ratio),
      original_energy_to_canonical_L2=enc(ct),near_input_upper_expression=enc(nearin),global_auxiliary_upper_expression=enc(ret),outer_auxiliary_upper_expression=enc(outer),
      upward_guards={'inverse':'.014256564','relative_response':'1.029218','near_input_per_E':'35.534','global_auxiliary_per_E':'36.573','outer_auxiliary_per_E':'23.485'},
      conditional_remaining_affine=enc(remainder),conditional_energy_comparison_ceiling=enc(gammaE),
      all_mode_scalar_domain_and_inverse=True,finite_mode_selection=False,common_subtracted_strong_sum_certified=False,
      raw_norm_sum_estimator='infinite; complete common-subtraction and derivative estimate still required',
      bolt_raw_source_estimator='unregulated omega=1/2 source has nonzero logarithm for admissible v0!=0; full renormalized response not declared singular',
      actual_physical_operator_applications=0,coupled_signed_solve_executed=False,actual_mass=None,kappa_m=None,
      other_coupled_columns_certified=False,actual_energy_feedback_comparison=None,M_graph_4=None,M_X_8=None,
      full_coupled_feasibility_decided=False,nonlinear_requirement_closed=False,physical_norm_target_changed=False,
      construction_obstruction=False,full_objective_preserved=True,closed_audit_preserved=True)

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
