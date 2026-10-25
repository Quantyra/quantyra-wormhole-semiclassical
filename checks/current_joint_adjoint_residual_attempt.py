"""Support the joint residual attempt; no full residual or coupled pass."""
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
STEM='current-joint-adjoint-residual-attempt'
OWN=['notes/'+STEM+'.md','checks/current_joint_adjoint_residual_attempt.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/current-reference-clock-adjoint-trial.md',
 'notes/current-coupled-mass-adjoint.md',
 'notes/current-complete-mass-moment-bound.md',
 'notes/coincident-normal-dirac-response.md',
 'notes/general-interior-dirac-response.md',
 'notes/strong-interior-quantum-remainder.md',
 'notes/stationary-seed-quantum-form.md',
 'notes/current-seed-full-hessian-realization.md',
 'notes/current-global-return-energy-attempt.md',
 'notes/stationary-centered-exterior-source.md',
 'notes/stationary-centered-canonical-response.md',
 'notes/transported-clock-dirac-response.md',
 'notes/current-controlled-solution-audit-closeout.md',
 'evidence/current-global-mass-conversion.json',
 'evidence/current-reference-clock-adjoint-trial.json',
 'evidence/current-complete-mass-moment-bound.json']
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
    xi=sp.sqrt(Q)*I;bc=1+sp.diff(Q,z)*I/2
    dxi=sp.sqrt(Q)*sp.diff(xi,z).subs(sp.diff(I,z),1/Q)
    check('actual_clock_proper_radial_identity',sp.simplify(dxi-bc)==0)
    check('actual_clock_tangential_identity',sp.simplify(sp.diff(Q,z)/(2*sp.sqrt(Q))*xi+1-bc)==0)
    ga,zh,be,a=sp.symbols('gamma zh beta a',positive=True)
    DA=a*zh**2+3*ga*zh;mt=ga**2/(be*zh*DA)
    LJ=3*zh**2/(4*ga**2)+a*zh**3/(4*ga**3)
    check('joint_mass_log_exactly_minus_half',sp.cancel(mt*be*(-2*ga)*LJ+sp.Rational(1,2))==0)
    check('individual_mass_poles_cancel',sp.cancel(be*zh**3/(2*ga)-2*ga*be*zh**3/(4*ga**2))==0)
    dd=sp.symbols('delta',positive=True)
    check('separate_clock_energy_pole',sp.simplify(2*ga*dd*(sp.diff(sp.log(dd)/2,dd))**2-ga/(2*dd))==0)
    check('nonzero_temperature_log_energy_divergence',sp.limit(-sp.log(dd),dd,0,dir='+')==sp.oo)
    S,T=sp.symbols('S T',real=True)
    check('disk_shape_pointwise_lower',sp.expand(8*(S+T)**2+8*S**2-4*T**2-4*(2*S+T)**2)==0)
    check('disk_scalar_pointwise_lower',sp.expand(8*(S+T)**2+8*S**2-8*S**2)==sp.expand(8*(S+T)**2))
    # Read current original-unit seed intervals, never replace by midpoints.
    mass=prior('current-global-mass-conversion')
    trial=prior('current-reference-clock-adjoint-trial')
    moment=prior('current-complete-mass-moment-bound')
    p=mass['parameters'];beta=dec(p['beta']);C=1/dec(p['u'])
    B2G=dec(trial['B_squared_over_G']);etaG=C/B2G
    Lq=dec(trial['dimensionless_complete_source_mass_factor'])
    Flmin=1-arb('1.179e-8')
    CT=arb('17.5')/(arb('.99')*Flmin).sqrt()
    Cnu=1/arb('23.99').sqrt()+35/(8*arb('.99')*Flmin).sqrt()
    check('original_shape_output_dual_comparison',CT<arb('17.6'))
    check('original_scalar_output_dual_comparison',Cnu<arb('12.65'))
    check('retained_actual_commutator_separation',arb('-.05').exp()>arb('.95'))
    delta=arb('.95')/128
    K1=arb(2)/3*arb('2e19')+arb('1.4e14')*24/delta**3+arb('2e19')
    Kold=arb('6e22')
    check('complete_interior_first_variation_coefficient',K1<arb('8.252e21'))
    check('same_units_coefficient_improvement',K1<Kold and Kold/K1>arb('7.27'))
    cout=arb('35.2')*etaG/beta
    cint=cout*K1/arb(5).sqrt()
    dden=16*arb('1e15')*Lq*arb('17.6')/arb(3).sqrt()
    check('finite_area_density_must_be_budgeted',dden>0 and dden<arb('6.065e-10'))
    target=arb('3.27169e-5');fullceiling=dec(moment['sufficient_full_adjoint_residual_ceiling'])
    check('unchanged_mass_residual_target',target<fullceiling)
    umax=target/cint
    check('all_budget_uniform_U_allowance',umax>7164 and umax<arb('7164.01'))
    out_slice=(target-cint*1000-dden*10-arb('1e-5'))/cout
    check('hypothetical_outer_allowance_positive',out_slice>arb('1.46618e25'))
    example=cint*1000+cout*arb('1e25')+dden*10+arb('1e-5')
    check('hypothetical_joint_ledger_fits',example<arb('2.69480e-5') and example<target)
    source=dec(moment['complete_reference_adjoint_pairing_upper'])
    conditional_mass=(source+arb('1.1e-10')*example)/arb('1.9999999764')
    check('only_conditional_mass_fits_conversion_split',conditional_mass<arb('2.29e-15'))
    check('unchanged_conditional_mass_column_affine_cost',sp.Rational('5e-5')*sp.Rational('4.2953e-8')==sp.Rational('2.14765e-12'))
    check('unchanged_conditional_mass_column_graph_cost',sp.Rational('5e-5')*sp.Rational('49.462')==sp.Rational('.0024731'))
    check('unchanged_two_fast_slot_screen',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    return {'precision_bits':ctx.prec,'actual_eta_G':enc(etaG),
      'shape_output_comparison_pre_round':enc(CT),'scalar_output_comparison_pre_round':enc(Cnu),
      'retained_commutator_separation':enc(delta),'interior_all_mode_coefficient':enc(K1),
      'old_to_new_coefficient_ratio':enc(Kold/K1),'uniform_U_residual_coefficient':enc(cint),
      'outer_weighted_norm_residual_coefficient':enc(cout),'density_B0_residual_coefficient':enc(dden),
      'all_other_costs_zero_uniform_U_allowance':enc(umax),
      'hypothetical_slice_inputs':{'U0':'1000','B0':'10','R_local_and_boundary':'1e-5','B_out':'1e25','ALL_UNPROVED':True},
      'hypothetical_slice_outer_allowance':enc(out_slice),'hypothetical_slice_complete_residual_upper':enc(example),
      'hypothetical_slice_mass_upper':enc(conditional_mass),
      'actual_unknowns':['U0','B0','B_out','R_local_and_boundary','full_adjoint_residual','actual_mass','kappa_m','other_coupled_columns'],
      'scope':{'original_fixed_seed_source_to_energy_comparison':True,
        'separate_global_clock_mass_estimator_admissible':False,'joint_trial_domain_rejected':False,
        'full_adjoint_residual_norm_computed':False,'actual_mass_certified':False,'actual_kappa_m_computed':False,
        'complete_coupled_return':'unresolved','new_admission_requirement_closed':False,
        'physical_neighborhood_comparison':False,'nonlinear_admission':False,'all_sector_stability':False,
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
    message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),'result':run.get('result'),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
