"""Original source-to-point conversion; complete remaining quantum map uncomputed."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import sympy as sp
import flint
from flint import arb, ctx
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from operators.signed_local_core import SignedLocalCore
STEM='current-weighted-bolt-source-return'
OWN=['notes/'+STEM+'.md','checks/current_weighted_bolt_source_return.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+['operators/signed_local_core.py',
 'notes/current-global-mixed-schur-preconditioner.md','checks/current_global_mixed_schur_preconditioner.py',
 'evidence/current-global-mixed-schur-preconditioner.json',
 'notes/current-joint-adjoint-residual-attempt.md','evidence/current-joint-adjoint-residual-attempt.json',
 'notes/current-localized-bolt-adjoint-sensitivity.md','evidence/current-localized-bolt-adjoint-sensitivity.json',
 'notes/current-seed-full-hessian-realization.md',
 'notes/current-schur-polynomial-return-test.md',
 'notes/current-adjoint-uniform-estimator-test.md','evidence/current-adjoint-uniform-estimator-test.json']
CHECKS=[]
def sh(p):return sha256(p.read_bytes()).hexdigest()
def check(name,ok):
    CHECKS.append({'name':name,'passed':bool(ok)})
    print(name+(': passed' if ok else ': FAILED'),flush=True)
    if not ok:raise AssertionError(name)
def enc(x):
    m,r,e=x.mid_rad_10exp(50)
    return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(s):return arb(s.strip('[]'))
def prior(stem):
    r=json.loads((ROOT/'evidence'/f'{stem}.json').read_bytes())['runs'][-1]
    check('accepted_'+stem,r['status']=='passed')
    return r
def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('frozen_argument_and_checker',all(sh(ROOT/p)==h for p,h in review['scientific_sha256'].items()))
    mixed=prior('current-global-mixed-schur-preconditioner')
    check('retained_complete_mixed_remainder_guard',dec(mixed['result']['numerics']['complete_mixed_preconditioned_remainder_upper'])<arb('1.556e-23'))
    check('accepted_mixed_argument_unchanged',all(sh(ROOT/p)==mixed['input_sha256_before'][p] for p in
      ['notes/current-global-mixed-schur-preconditioner.md','checks/current_global_mixed_schur_preconditioner.py']))
    output=prior('current-joint-adjoint-residual-attempt')
    check('accepted_original_L2_comparison_unchanged',sh(ROOT/'notes/current-joint-adjoint-residual-attempt.md')==output['input_sha256_before']['notes/current-joint-adjoint-residual-attempt.md'])
    src=prior('current-localized-bolt-adjoint-sensitivity')
    prev=prior('current-adjoint-uniform-estimator-test')
    check('previous_uniform_failure_preserved',prev['result']['scope']['specified_complete_positive_estimator_rejected_for_every_trial'])
    pars=src['result']['parameters']
    h,gamma,beta,C,B2=[dec(pars[k]) for k in ['h','gamma','beta','C','B2']]
    g=C*beta/B2;Fl=1-384*arb.pi()*arb('1e22')*g
    d=(96*arb.pi()*arb('1e22')-1/(60*arb.pi()))*C/B2
    s=arb(7)/5;gap=12*beta-s*s;mD=1-d*s*s
    check('actual_parameter_and_diffusion_guards',0<h<arb('1e-30') and gamma>0 and beta<arb('.17') and arb('1e-8')<d<arb('2e-8'))
    check('positive_reference_factor',0<Fl<1 and Fl>1-arb('1.179e-8'))
    check('mixed_gap_and_L2_bounds',gap>arb('.03999999995') and 1/gap<26 and 1/gap.sqrt()+1/gap<32)
    Vmax=2/C+12*beta*(1+11*g/(360*arb.pi()))
    zeroth=Vmax+54*beta+10+s*s+2*s
    check('adjoint_zeroth_coefficient',Vmax<arb('2.001') and zeroth<26)
    lapP=1+2*(s+1)*32+26*26
    check('adjoint_laplacian_and_Hessian',lapP<850 and (arb(850)**2+13*beta*32**2).sqrt()<900)
    check('massive_energy_bounds',1/mD<arb('1.01') and 1/(d*mD).sqrt()<11000)
    lapD=arb('2.01')/d+2*s*11000+(10+s*s)*arb('1.01')
    check('massive_laplacian_and_Hessian',lapD<arb('2.02e8') and (arb('2.02e8')**2+13*beta*11000**2).sqrt()<arb('2.03e8'))
    R=arb('.01');pi=arb.pi()
    aR=2/((3*pi).sqrt()*R);bR=2/(3*pi).sqrt();cR=R/(2*pi.sqrt())
    KP=(2*pi).sqrt()*(26*aR+32*bR+900*cR)
    KD=(2*pi).sqrt()*(arb('1.01')*aR+11000*bR+arb('2.03e8')*cR)
    check('global_P_adjoint_point_factor',KP<5000)
    check('global_LD_point_factor',KD<arb('2e6'))
    minf=1-d*(10+s*s);rm=arb('1.556e-23')
    check('positive_weighted_maximum_principle_barrier',minf>arb('.99999976'))
    KC=(arb(5000)/minf+arb('2e6')*26*rm/(1-rm))/2
    check('complete_mixed_point_factor',KC<3000)
    # Independent exact reference T/T Hessian: every derivative entry is retained.
    core=SignedLocalCore();z=core.z;u=1/core.C;be=core.beta;ep=core.eps
    Vref=u/z**2+3*be/z**4-(u+3*be+ep)+4*ep*z-3*ep*z**2
    target=4*core.Fl*Vref
    weak=[[sp.expand(sp.diff(core.reference_density,core.fields[1][i],core.fields[1][j]).subs(core.zero))
          for j in range(3)] for i in range(3)]
    check('complete_reference_TT_row_is_multiplication',all(sp.simplify(weak[i][j]-(target if i==j==0 else 0))==0 for i in range(3) for j in range(3)))
    fact=(1-z)*(u*z**2*(1+z)+3*be*(1+z+z**2+z**3)-ep*z**4*(1-3*z))
    check('exact_whole_interval_reference_zero',sp.expand(z**4*Vref-fact)==0)
    check('reference_original_row_bound',4*(2+4*arb('1e-28'))<arb('8.1'))
    check('original_T_output_comparison',arb(35)/(2*(arb('.99')*Fl).sqrt())<arb('17.6'))
    # For x>=0, (1+x)^(7/5)>=1+x>=x, so no sampled supremum is used.
    zh=1/(1+h);yh=h/(1+h)
    # Preserve the common h: direct interval subtraction loses the O(h) strict gap.
    hh=sp.symbols('h',positive=True)
    identity=sp.Rational(5,2)-(1+1/(1+hh))/(sp.Rational(4,5)+hh)
    factored=hh*(sp.Rational(7,2)+sp.Rational(5,2)*hh)/((sp.Rational(4,5)+hh)*(1+hh))
    check('unchanged_bolt_weight_conversion',sp.cancel(identity-factored)==0 and h*(arb('3.5')+arb('2.5')*h)/((arb('.8')+h)*(1+h))>0)
    q=arb(169)/500;f=q/(1+q);x0=1+2/q;c0,c1=arb(1),x0
    for j in range(1,29):c0,c1=c1,2*x0*c1-c0
    h2=(1+q)**2/c1;h3=(1+q)**3/c1;pref=arb('1.1e-10')*h3
    ellb=(1+h)/h
    Kout=ellb*(arb('2.5')*h)**s*3000/gamma.sqrt()
    Mref=arb('8.1')*arb('17.6')*f**3*h**(-arb(2)/5)
    Cref=Kout*Mref;refdebit=pref*Cref;budget=arb('1.667e-9')
    ceiling=budget/(pref*Kout)-Mref
    check('complete_reference_tail_original_bolt_test',refdebit<arb('1e-18'))
    check('remaining_full_source_allowance_positive',ceiling>0)
    check('full_budget_identity',(pref*Kout*(Mref+ceiling)).overlaps(budget))
    return {'precision_bits':512,'physical_second_order_parameters':{'h':enc(h),'gamma':enc(gamma),'d':enc(d),'p0':enc(gap),'m_inf':enc(minf)},
      'actual_output_factors':{'P_adjoint_point_numerator':enc(KP),'LD_point_numerator':enc(KD),
       'complete_mixed_point_numerator_upper_expression':enc(KC),'complete_mixed_point_numerator_rounded_upper':3000,
       'original_physical_bolt_factor_from_weighted_L2':enc(Kout)},
      'reference_three_residual_map':{'M_ref_upper':enc(Mref),'C_b3_ref_upper':enc(Cref),
       'original_affine_bolt_debit_upper':enc(refdebit),'fraction_of_full_bolt_margin':enc(refdebit/budget),
       'per_independent_input_tail_coefficient':enc(Cref*h2)},
      'complete_remaining_map':{'definition':'||W (Dbar-Dbar_ref) V F^3||_(E to L2(z^-4 dz))',
       'actual_norm':None,'conditional_ceiling_full_bolt_margin':enc(ceiling),
       'physical_affine_debit_per_unit_norm':enc(pref*Kout),'full_bolt_margin_only':enc(budget)},
      'h2':enc(h2),'h3':enc(h3),'affine_prefactor':enc(pref),
      'retained_other_conditional_ceilings':prev['result']['conditional_original_space_ceilings'],
      'scope':{'actual_weighted_mixed_output_factor_certified':True,'reference_row_complete_F3_contribution_certified':True,
       'remaining_complete_quantum_source_map_computed':False,'joint_transferred_residual_executed':False,
       'full_Cb3_certified':False,'original_strong_delta_certified':False,'graph_physical_maps_certified':False,
       'actual_kappa_m_computed':False,'signed_coupled_columns_computed':False,
       'new_controlled_solution_requirement_closed':False,'actual_controlled_solution':False,'construction_obstruction':False,
       'nonlinear_admission_contraction':False,'exact_compatibility':False,'physical_matching':False,
       'semiclassical_validity_noise':False,'all_sector_stability':False,'full_goal_complete':False}}
if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={p:sh(ROOT/p) for p in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
      'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
      'runtime':{'python':sys.version,'python_sha256':sh(Path(sys.executable)),'flint_version':flint.__version__,'sympy_version':sp.__version__}}
    try:
        run['result']=execute();run['input_sha256_after']={p:sh(ROOT/p) for p in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    summary=json.dumps({'status':run['status'],'mode':mode,'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as out:out.write(summary.encode())
    print(summary,flush=True)
    if run['status']!='passed':raise SystemExit(1)
