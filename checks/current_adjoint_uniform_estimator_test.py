"""Support the all-trial estimator proof; never evaluate a quantum inverse."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
import sympy as sp
import flint
from flint import arb, ctx
ROOT=Path(__file__).resolve().parents[1]
STEM='current-adjoint-uniform-estimator-test'
OWN=['notes/'+STEM+'.md','checks/current_adjoint_uniform_estimator_test.py',
     'evidence/'+STEM+'-review.json']
SOURCE='current-localized-bolt-adjoint-sensitivity'
INPUTS=OWN+['notes/'+SOURCE+'.md','checks/current_localized_bolt_adjoint_sensitivity.py',
 'evidence/'+SOURCE+'.json','notes/current-coupled-adjoint-residual-audit.md',
 'evidence/current-coupled-adjoint-residual-audit.json',
 'notes/current-direct-triple-return-attempt.md',
 'notes/current-schur-polynomial-return-test.md',
 'notes/current-geodesic-two-sided-remainder.md',
 'notes/current-seed-full-hessian-realization.md']
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
def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('frozen_argument_and_checker',all(sh(ROOT/p)==h for p,h in review['scientific_sha256'].items()))
    src=json.loads((ROOT/('evidence/'+SOURCE+'.json')).read_bytes())['runs'][-1]
    check('accepted_actual_domain_source_test',src['status']=='passed')
    check('source_ordinary_proof_and_checker_unchanged',all(sh(ROOT/p)==src['input_sha256_before'][p]
      for p in ['notes/'+SOURCE+'.md','checks/current_localized_bolt_adjoint_sensitivity.py']))
    prior=json.loads((ROOT/'evidence/current-coupled-adjoint-residual-audit.json').read_bytes())['runs'][-1]
    check('accepted_complete_residual_propagation',prior['status']=='passed' and prior['result']['complete_triple_residual_inverse_upper']==1)
    # Noncommutative block algebra. Adjoint signs/normalization are in the ordinary proof.
    A,Ai,C,Cs,D,U,V,Vs=sp.symbols('A Ai C Cs D U V Vs',commutative=False)
    b=-A*U-C*V
    expr=sp.expand(-Ai*b-U-Ai*C*V).subs(Ai*A,sp.Integer(1))
    check('effective_source_minus_sign',sp.expand(expr)==0)
    # Equation-index recursion only; t represents the full self-adjoint T.
    t=sp.symbols('t',nonnegative=True)
    chain=sp.Matrix([[1+t,-t,0],[0,1+t,-t],[0,0,1+t]])
    row=t*chain.inv()[0,:]
    check('exact_joint_residual_transfer',all(sp.cancel(row[j]-(t/(1+t))**(j+1))==0 for j in range(3)))
    d,V0,eta,M=sp.symbols('d V0 eta M',positive=True)
    check('all_trial_minimum_identity',sp.expand(d+eta*V0-eta/M*(d+M*V0)-(1-eta/M)*d)==0)
    r=src['result'];pars=r['parameters']
    h,beta=dec(pars['h']),dec(pars['beta'])
    b_low=dec(r['effective_point_lower']);Lc=dec(r['cold_energy_upper_over_action'])
    ell_b=(1+h)/h;Jstar=ell_b*b_low*(beta/Lc).sqrt()
    check('actual_test_positive_and_finite',h>0 and beta>0 and b_low>0 and Lc>0 and Jstar.is_finite())
    normchecks=[]
    for entry in r['norms']:
        converted=dec(entry['effective_functional_dual_lower'])*(dec(entry['W_abs'])*beta).sqrt()
        normchecks.append({'folded_halves':entry['N'],'converted_E_dual_lower':enc(converted)})
        check('unchanged_original_units_N'+str(entry['N']),converted.overlaps(Jstar))
    q=arb(169)/500;f=q/(1+q);eta=arb('0.0022');M=1+q+eta
    x0=1+2/q;c0,c1=arb(1),x0
    for j in range(1,29):c0,c1=c1,2*x0*c1-c0
    D29=c1;h2=(1+q)**2/D29;h3=(1+q)**3/D29
    pref=arb('1.1e-10')*h3;budget=arb('1.667e-9');star=budget/pref
    d3floor=eta*Jstar/M
    sharpfloor=f**3*d3floor;oldfloor=q*d3floor
    physical=pref*sharpfloor;oldphysical=pref*oldfloor
    ratio=physical/budget;oldratio=oldphysical/budget
    check('stronger_uniform_estimator_fails_complete_bolt_margin',physical>budget)
    check('older_positive_estimator_also_fails',oldphysical>physical)
    threshold=star*(1+q)/(f**3*Jstar-star)
    oldthreshold=star*(1+q)/(q*Jstar-star)
    check('necessary_threshold_denominators_positive',f**3*Jstar>star and q*Jstar>star)
    check('threshold_substitution', (pref*f**3*threshold*Jstar/(1+q+threshold)).overlaps(budget))
    hypo=arb('1e-10');hypophysical=pref*f**3*hypo*Jstar/(1+q+hypo)
    check('hypothetical_1e_minus10_still_fails_estimator',hypophysical>budget)
    conditional={'Cdef3_times_IY':enc(arb('0.5')/h3),
      'Cg3_full_graph_margin':enc(arb('5.34235e-12')/pref),
      'CX3_full_physical_margin':enc(arb('5e-9')/pref),
      'Cb3_full_bolt_margin':enc(star),
      'Cmass3_times_Im_zero_center_and_error':enc(arb('5e-5')/h2),
      'Cg3_times_Ifb_zero_other_debits':enc(arb('0.0012069')/h2)}
    result={'precision_bits':512,'q':enc(q),'f':enc(f),'eta_available_upper':enc(eta),
      'L0_operator_upper':enc(M),'ell_b_original_bolt_factor':enc(ell_b),
      'accepted_source_test':{'b_low':enc(b_low),'Lc':enc(Lc),'Jstar_E_dual_lower':enc(Jstar),
                              'independent_normalization_checks':normchecks},
      'D29':enc(D29),'h2':enc(h2),'h3':enc(h3),'affine_prefactor':enc(pref),
      'all_trial_uniform_estimator':{'D3_report_floor':enc(d3floor),
       'improved_complete_Cb3_report_floor':enc(sharpfloor),
       'improved_physical_bolt_report_floor':enc(physical),
       'improved_ratio_to_full_bolt_margin':enc(ratio),
       'older_q_positive_physical_report_floor':enc(oldphysical),
       'older_ratio_to_full_bolt_margin':enc(oldratio),
       'known_power_gain_q_over_f3':enc(q/f**3),
       'necessary_new_uniform_eta_ceiling':enc(threshold),
       'necessary_old_uniform_eta_ceiling':enc(oldthreshold),
       'required_uniform_eta_improvement':enc(eta/threshold),
       'hypothetical_eta_1e_minus10_physical_report_floor':enc(hypophysical)},
      'conditional_original_space_ceilings':conditional,
      'complete_entries':{k:'unknown signed center, certified error, input norm and original output norm'
        for k in ['mass','scalar_quartic','correlated_slow_pair','fast_remainder','infinite_matching_and_boundary_phase']},
      'scope':{'specified_complete_positive_estimator_rejected_for_every_trial':True,
       'actual_Cb3_lower_bound_proved':False,'actual_quantum_error_lower_bound_proved':False,
       'actual_trial_constructed':False,'joint_quantum_residual_evaluated':False,
       'direct_composite_norms_enclosed':False,'original_strong_delta_certified':False,
       'actual_kappa_m_computed':False,'signed_coupled_columns_computed':False,
       'new_controlled_solution_requirement_closed':False,'actual_controlled_solution':False,
       'construction_obstruction':False,'nonlinear_admission_contraction':False,
       'exact_compatibility':False,'physical_matching':False,'semiclassical_validity_noise':False,
       'all_sector_stability':False,'full_goal_complete':False}}
    return result
if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={p:sh(ROOT/p) for p in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,
      'input_sha256_before':hashes,'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
      'runtime':{'python':sys.version,'python_sha256':sh(Path(sys.executable)),
                 'flint_version':flint.__version__,'sympy_version':sp.__version__}}
    try:
        run['result']=execute()
        run['input_sha256_after']={p:sh(ROOT/p) for p in INPUTS}
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
