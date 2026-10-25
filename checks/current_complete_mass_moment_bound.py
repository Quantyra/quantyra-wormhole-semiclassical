"""Complete seed-source moment bound, not a Hessian or nonlinear solution."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys,traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb,ctx
import flint
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
STEM='current-complete-mass-moment-bound'
OWN=['notes/'+STEM+'.md','checks/current_complete_mass_moment_bound.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/current-reference-clock-adjoint-trial.md',
 'notes/current-coupled-mass-adjoint.md',
 'notes/stationary-reference-nonbolt-stress.md',
 'notes/stationary-reference-finite-mass.md',
 'notes/stationary-reference-cold-bolt.md',
 'notes/stationary-massive-boundary-layer.md',
 'notes/stationary-centered-exterior-source.md',
 'notes/fixed-clock-seed-bolt-feasibility.md',
 'notes/current-global-return-energy-attempt.md',
 'notes/current-controlled-solution-audit-closeout.md',
 'evidence/current-global-mass-conversion.json',
 'evidence/current-reference-clock-adjoint-trial.json',
 'evidence/current-coupled-mass-adjoint.json']
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
    return r['result']

def execute():
    ctx.prec=512
    z=sp.symbols('z',positive=True)
    Q=sp.Function('Q')(z);I=sp.Function('I')(z)
    check('mass_weight_primitive',sp.simplify(sp.diff(Q*I,z).subs(sp.diff(I,z),1/Q)-1-sp.diff(Q,z)*I)==0)
    C=sp.symbols('C',positive=True)
    # Bfun=CQ, d/dt=-z^2 d/dz, H(t)=I_Q/C.
    check('original_mass_kernel_measure',sp.cancel((z**2+sp.Abs(-C*z**2*sp.diff(Q,z))*I/(2*C))/z**2-1-sp.Abs(sp.diff(Q,z))*I/2)==0)
    s=sp.symbols('s',positive=True)
    gaussian=sp.integrate(sp.exp(-s*z**2),(z,0,sp.oo))
    check('integrated_E1_weight_identity',sp.simplify(sp.integrate(gaussian/s,(s,1,sp.oo))-sp.sqrt(sp.pi))==0)
    check('bridge_moment_substitution',256*48+128*12+16*4==13888)

    pi=arb.pi();exp=lambda x:arb(x).exp()
    local=arb('4e11')*arb(13888).sqrt()/(2*pi)
    exit_ratio=arb(10000)*exp('-8.75')
    far_ratio=arb('1e10')*64/(9*exp(2)*arb('.1')**4)
    cubic_ratio=2*arb(2).sqrt()*exp('-.5')/arb('.1')
    check('local_boundary_value_remainder',local<arb('8e12'))
    check('exit_value_ratio',exit_ratio<arb('1.585'))
    check('far_crossing_value_ratio',far_ratio<arb('9.624e13'))
    check('far_cubic_value_ratio',cubic_ratio<18)
    Cv=arb('2e14')
    check('whole_reflected_value_constant_local',local+arb('2e10')*exit_ratio+1<Cv)
    check('whole_reflected_value_constant_far',far_ratio+cubic_ratio+1<Cv)
    check('global_FK_value_factor',(arb('5e10')*exp('.01')/(2*pi))<arb('1e10'))
    # Both local exit cap values and the cubic tail fit this prefactor.
    check('exit_FK_and_cubic_factor',arb('5e10')*exp('.01')/pi+1<arb('2e10'))
    check('local_density_remainder_constant',3*arb('1.3e11')+24*arb('.4')<arb('4e11'))
    check('all_fixed_restorations_fit_200',167+3+80*arb('.0035')**2+arb('1.6')*arb('.0035')+1<200)

    mass=prior('current-global-mass-conversion')
    trial=prior('current-reference-clock-adjoint-trial')
    adj=prior('current-coupled-mass-adjoint')
    pp=mass['parameters'];be=dec(pp['beta']);u=dec(pp['u']);eps=dec(pp['eps']);tau=dec(pp['tau'])
    q=arb(2000)**4;Ca=1/u;rc=Ca/q
    check('current_beta_lower',be>arb(1)/7)
    check('current_angular_ratio',rc<arb('.0035'))
    check('retained_positive_thermal_mass',tau>0 and tau<arb('1e-65'))
    check('increasing_reference_factor',2*be-4*abs(eps)>0)
    check('Qprime_negative_beyond_one_hundredth',-12*be*arb('.01')+2*u+4*abs(eps)<0)
    check('positive_Qprime_near_zero',2*u+4*abs(eps)<arb('5e-11'))
    check('actual_outer_Q_lower',arb('.99')**2*be-tau>arb('.16'))
    dc=arb(1)/1001;zc=1-dc
    relative_thermal=tau/(be*dc**2)
    check('nonbolt_thermal_relative_error',relative_thermal<arb('1e-58'))
    check('retained_antiderivative_factor',1/(1-relative_thermal)<arb('1.001'))
    qI_ratio=arb('1.001')*(2+3*tau/dc)/be
    check('negative_Qprime_times_primitive',qI_ratio<15)
    positive_weight=arb('5e-11')/arb('.16')*arb('.01')**3/3
    check('complete_mass_weight_integral',arb('1.5')+positive_weight<arb('1.501'))
    cold_weight=arb('1.5')*dc+arb('.5')*arb('1.001')*zc*dc/be
    check('cold_weight_without_inverse_temperature',cold_weight<arb('.006'))
    check('linear_trace_model_integral',15*zc+2*arb('.006')/zc<16)
    check('linear_trace_and_Weyl_model_cost',arb(16)/30<1)

    B0=arb('8e15')*rc**3+arb('1.4e13')*rc+200
    A0=24*Cv*rc**2+16*Ca**arb('1.5')/q**2
    check('complete_trace_constant',B0<arb('5e10'))
    E1_anchor_bound=(-125*zc**2).exp()/(125*zc**2)
    check('cold_anchor_reflection_cost',A0*E1_anchor_bound<1)
    check('existing_cold_derivative_sufficient',arb('5e10')+1+arb('.001')*(arb('3e14')+arb(1)/15)<arb('3.51e11'))
    kernel_weight=9*pi.sqrt()/(5*arb(5).sqrt())
    constant_mass=arb('1.501')*B0
    reflected_mass=kernel_weight*A0
    cold_mass=arb('2.106e9')
    Cmass=constant_mass+reflected_mass+cold_mass+10010
    Lq=dec(trial['dimensionless_complete_source_mass_factor'])
    locals=dec(trial['local_absolute_mass_allowance'])
    bolt=dec(trial['bolt_absolute_mass_allowance'])
    flux=2*(Cmass*Lq+locals+bolt)
    old=dec(trial['complete_flux_absolute_envelope'])
    Km=arb('1.9999999764')
    ceiling=(Km*arb('2.29e-15')-flux)/arb('1.1e-10')
    conditional_mass=(flux+arb('1.1e-10')*arb('1e-5'))/Km
    check('complete_quantum_mass_coefficient',Cmass<arb('1.7e11'))
    check('complete_pairing_budget_closed',flux<arb('1.3e-15') and flux<arb('3.3e-15'))
    check('new_bound_actually_improves_old_estimator',flux<old)
    check('sufficient_full_residual_ceiling',ceiling>arb('2.98e-5'))
    check('conditional_actual_mass_fits_unchanged_target',conditional_mass<arb('2.29e-15'))
    check('full_mass_flux_still_original',dec(adj['complete_mass_flux'])>Km)
    return {'precision_bits':ctx.prec,'current_C_over_q':enc(rc),
      'new_boundary_heat_value_constant':enc(Cv),'local_bridge_value_constant':enc(local),
      'far_crossing_value_constant':enc(far_ratio),'complete_trace_B0':enc(B0),'complete_trace_A0':enc(A0),
      'integrated_log_boundary_weight':enc(kernel_weight),'constant_trace_mass_allowance':enc(constant_mass),
      'reflected_remainder_mass_allowance':enc(reflected_mass),'cold_mass_allowance':enc(cold_mass),
      'complete_quantum_mass_coefficient':enc(Cmass),'complete_reference_adjoint_pairing_upper':enc(flux),
      'old_envelope_to_new_envelope_ratio':enc(old/flux),
      'sufficient_full_adjoint_residual_ceiling':enc(ceiling),
      'conditional_actual_mass_if_residual_le_1e_minus5':enc(conditional_mass),
      'scope':{'complete_source_pairing_allowance_closed':True,'narrow_signed_quadrature_executed':False,
        'pairing_sign_determined':False,'full_adjoint_residual_norm_computed':False,
        'actual_mass_certified':False,'actual_kappa_m_computed':False,'complete_coupled_return':'unresolved',
        'physical_norm_comparisons_closed':False,'nonlinear_admission':False,'all_sector_stability':False,
        'construction_obstruction':False,'full_objective_complete':False}}

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
