"""Actual reference adjoint and complete-source envelope; no full quantum solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb,ctx
import flint
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
STEM='current-reference-clock-adjoint-trial'
OWN=['notes/'+STEM+'.md','checks/current_reference_clock_adjoint_trial.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+[
    'notes/current-coupled-mass-adjoint.md',
    'notes/current-coupled-global-return-attempt.md',
    'notes/current-seed-full-hessian-realization.md',
    'notes/current-global-return-energy-attempt.md',
    'notes/current-reference-rescaled-source-test.md',
    'notes/fixed-clock-seed-bolt-feasibility.md',
    'notes/stationary-centered-exterior-source.md',
    'notes/current-controlled-solution-audit-closeout.md',
    'evidence/current-global-mass-conversion.json',
    'evidence/current-physical-york-dual-assembly.json',
    'evidence/current-coupled-global-return-attempt.json',
    'evidence/current-coupled-mass-adjoint.json',
]
CHECKS=[]

def check(name,value):
    CHECKS.append({'name':name,'passed':bool(value)})
    if not value:raise AssertionError(name)

def enc(v):
    m,r,e=v.mid_rad_10exp(50)
    return {'ball':v.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}

def dec(v):
    d=v['decimal_enclosure']
    return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])

def prior(stem):
    r=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
    check('accepted_input_'+stem,r['status']=='passed')
    return r

def execute():
    ctx.prec=512
    z=sp.symbols('z',positive=True); be,C=sp.symbols('beta C',positive=True)
    Q=sp.Function('Q')(z); I=sp.Function('I')(z); J=sp.Function('J')(z)
    FF=C*Q/z**2
    nc=-Q*I/z;bc=1+sp.diff(Q,z)*I/2
    nm=-2*be*Q*J/z;bm=be*z**3/Q+be*sp.diff(Q,z)*J
    def red(ex):
        ex=ex.subs(sp.diff(I,z),1/Q).subs(sp.diff(J,z),z**3/Q**2)
        return sp.cancel(ex)
    f=lambda n,b:2*FF*(-z*sp.diff(n,z)-b)+z*sp.diff(FF,z)*n
    al=lambda n,b:n+2*b+z*sp.diff(n,z)
    check('exact_clock_reference_mass_row',red(f(nc,bc))==0)
    check('exact_clock_reference_lapse',red(al(nc,bc)-1)==0)
    check('exact_mass_reference_mass_row',red(f(nm,bm)-2*C*be*z)==0)
    check('exact_mass_reference_lapse',red(al(nm,bm))==0)

    zh,ga,q2,dd=sp.symbols('zh gamma q2 delta',positive=True)
    DA=zh**2*q2/2+3*ga*zh
    mt=ga**2/(be*zh*DA); cc=zh**2/(2*DA)
    c_log=3*zh**2/(4*ga**2)+q2*zh**3/(8*ga**3)
    check('exact_bolt_log_cancellation',sp.cancel(-1/(2*ga)+2*be*mt*c_log)==0)
    check('coordinate_temperature_matching',sp.cancel(1-mt*DA*be*zh/ga**2)==0)
    e2,e3,e4,e5=sp.symbols('e2 e3 e4 e5')
    U=2*ga+q2*dd/2+e2*dd**2+e3*dd**3+e4*dd**4+e5*dd**5
    numerator=dd*U+2*be*mt*(zh-dd)**3-cc*U**2
    check('regular_primitive_double_zero_value',sp.cancel(numerator.subs(dd,0))==0)
    check('regular_primitive_double_zero_derivative',sp.cancel(sp.diff(numerator,dd).subs(dd,0))==0)
    Qprime=-(U+dd*sp.diff(U,dd))
    bnum=cc*Qprime*U/2+be*mt*(zh-dd)**3
    check('regular_shape_numerator_zero',sp.cancel(bnum.subs(dd,0))==0)
    check('reference_bolt_radius_matching',sp.cancel(-2*ga*cc/zh+ga*zh/DA)==0)
    # Equality of the derivative of the pole-removed primitive, on delta>0.
    V=numerator/dd**2
    check('full_primitive_derivative_identity',
          sp.cancel(cc/dd**2+V/U**2-(dd*U+2*be*mt*(zh-dd)**3)/(dd*U)**2)==0)

    u,a,mt_s=sp.symbols('u a m_T')
    Qsmall=be+u*z**2-a*z**3
    Ismall=z/be-u*z**3/(3*be**2)+a*z**4/(4*be**2)
    ns=sp.series(-Qsmall*Ismall/z-mt_s*z**3/2,z,0,4).removeO().expand()
    bs=sp.series(1+sp.diff(Qsmall,z)*Ismall/2+mt_s*z**3,z,0,4).removeO().expand()
    check('auxiliary_weight_two_and_clock_root',sp.simplify(ns+1+2*u*z**2/(3*be)-(3*a/(4*be)-mt_s/2)*z**3)==0)
    check('canonical_cubic_is_not_m_T',sp.simplify(bs-1-u*z**2/be-(mt_s-3*a/(2*be))*z**3)==0)
    Fl,m1=sp.symbols('F_l m1')
    check('current_reference_prefactor_cancels_in_pairing',sp.cancel(2*Fl*(m1/Fl)-2*m1)==0)
    t=sp.symbols('t',positive=True)
    check('local_potential_mass_integral',sp.integrate(t**-2,(t,1,sp.oo))==1)
    check('local_square_mass_integral',sp.integrate(t**-3,(t,1,sp.oo))==sp.Rational(1,2))

    mass=prior('current-global-mass-conversion')['result']
    yp=prior('current-physical-york-dual-assembly')
    glob=prior('current-coupled-global-return-attempt')['result']
    adj=prior('current-coupled-mass-adjoint')['result']
    pp=mass['parameters']; bea=dec(pp['beta']);ua=dec(pp['u']);ta=dec(pp['tau']);ag=dec(pp['g'])
    ha=dec(glob['h']);gamma=dec(glob['gamma']);DAa=dec(glob['D_A_squared']);zha=1/(1+ha)
    oldep=next(v for v in yp['numerics'] if v['kind']=='actual_seed_correlated_endpoint')
    oldW=next(v for v in yp['numerics'] if v['kind']=='physical_action_and_adjoint_normalization' and v['identical_exteriors']==1)
    # W=-pi (B^2/G)/gamma for one half. This extracts the same old unit,
    # retaining correlations as enclosing intervals, not as midpoint equalities.
    B2G=-dec(oldW['W'])*dec(oldep['gamma'])/arb.pi()
    check('same_positive_temperature_enclosures',ha.overlaps(dec(oldep['h'])) and gamma.overlaps(dec(oldep['gamma'])))
    check('current_sphere_scale',B2G>arb('8.38e42') and B2G<arb('8.39e42'))
    mTa=gamma**2/(bea*zha*DAa)
    phase=mTa-3*(4*bea+2*ua+ta)/(2*bea)
    nbolt=-gamma*zha/DAa
    check('reference_areal_clock_mass_scale',mTa>arb('1e-66') and mTa<arb('2e-66'))
    check('auxiliary_canonical_mass_cubic_order_one',phase>arb('-6.000001') and phase<arb('-5.999999'))
    check('bolt_radius_change_sign',nbolt<0)
    check('reference_mass_parameter_matches_temperature_derivative',
          (mTa*dec(glob['reference_mass_relative_temperature_derivative'])).overlaps(arb(-1)))

    Lq=arb(2000)**8/(ag*B2G**2)
    ep=arb('5.5e-23')
    es=arb('1e22')*ag*(200*arb('1.1e-27')+arb('1.1e-27')**2+10*arb('4.3e-25'))
    local=8*arb.pi()*ep+4*arb.pi()*es
    bolt=(gamma/bea)*(1+ha)*ha*arb('3.33e-10')
    quantum=arb('8.001e15')*Lq
    Iabs=2*(quantum+local+bolt)
    target=arb('3.3e-15')
    Cmass=(target-2*(local+bolt))/(2*Lq)
    ratio=Iabs/target
    check('complete_source_envelope_upper',Iabs<arb('7e-11'))
    check('specified_unsigned_estimator_exceeds_budget',ratio>arb('1e4'))
    check('required_source_moment_coefficient_positive',Cmass>0)
    check('cold_value_refinement_cannot_resolve_this_estimator',2*bolt<arb('1e-70'))
    check('new_local_terms_do_not_dominate_flux_budget',2*local<arb('3e-21'))
    check('actual_mass_target_not_redefined',dec(adj['conditional_mass_upper'])<arb('2.29e-15'))
    return {'precision_bits':ctx.prec,'areal_clock_mass_parameter':enc(mTa),
            'auxiliary_canonical_cubic':enc(phase),'reference_bolt_conformal_value':enc(nbolt),
            'B_squared_over_G':enc(B2G),'dimensionless_complete_source_mass_factor':enc(Lq),
            'quantum_absolute_mass_allowance':enc(quantum),'local_absolute_mass_allowance':enc(local),
            'bolt_absolute_mass_allowance':enc(bolt),'complete_flux_absolute_envelope':enc(Iabs),
            'specified_flux_envelope_to_target_ratio':enc(ratio),
            'sufficient_replacement_absolute_quantum_mass_coefficient':enc(Cmass),
            'required_improvement_of_absolute_quantum_coefficient':enc(arb('8.001e15')/Cmass),
            'scope':{'regular_reference_adjoint_trial':True,'reference_rows_exactly_zero':True,
              'full_operator_residual':'-J_rest_second_variation Psi_ref; norm unknown',
              'complete_source_unsigned_envelope_evaluated':True,'narrow_signed_pairing_computed':False,
              'actual_adjoint_residual_norm_computed':False,'actual_mass_certified':False,
              'actual_kappa_m_computed':False,'complete_physical_coupled_return':'unresolved',
              'construction_obstruction':False,'nonlinear_admission':False,'full_objective_complete':False}}

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,
      'input_sha256_before':hashes,'frozen_sources_utf8':{f:(ROOT/f).read_text(encoding='utf-8') for f in OWN},
      'runtime':{'python':sys.version,'executable_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
        'flint_version':flint.__version__,'sympy_version':sp.__version__}}
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
    summary={'mode':mode,'status':run['status'],'checks':len(CHECKS),'result':run.get('result'),'traceback':run.get('traceback')}
    message=json.dumps(summary,indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
