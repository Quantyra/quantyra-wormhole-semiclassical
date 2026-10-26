"""Exact two-row assembly checks; not a quantum or coupled-solve certificate."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb
import json,sys,traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import sympy as sp
import flint
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1]
STEM='current-signed-strong-local-core'
OWN=['notes/'+STEM+'.md','operators/signed_local_core.py',
     'checks/current_signed_strong_local_core.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+['notes/current-finite-local-action-assembly.md',
 'notes/current-coupled-global-return-attempt.md',
 'notes/current-angular-resummed-quantum-operator.md',
 'notes/current-curved-paired-heat-bridge.md',
 'notes/current-rough-paired-variation.md',
 'notes/current-original-cutoff-coupled-estimator.md',
 'notes/current-fast-graph-domain-embedding.md']
sys.path.insert(0,str(ROOT))
from operators.signed_local_core import SignedLocalCore
CHECKS=[]
def check(name,ok):
    CHECKS.append({'name':name,'passed':bool(ok)})
    print(name+(': passed' if ok else ': FAILED'),flush=True)
    if not ok:raise AssertionError(name)
def zero(e):return sp.cancel(sp.expand(e))==0
def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[3]).read_bytes())
    for p,h in review['scientific_sha256'].items():
        check('frozen_'+p,sha256((ROOT/p).read_bytes()).hexdigest()==h)
    core=SignedLocalCore();z=core.z
    independent=core.independent_trace_variation()
    coefficients={}
    for name in core.densities:
        built=core.build(name);H=built['strong'];weak=built['weak']
        check(name+'_independent_complete_trace_row',all(
            zero(H[0][b][k]-independent[name][b][k]) for b in range(2) for k in range(5)))
        # Exact formal adjoint in ORIGINAL dz, including variable coefficients.
        check(name+'_two_row_formal_adjoint',all(zero(H[a][b][j]-sum(
            (-1)**j*comb(k,j)*core._adjoint_power(H[b][a][k],k-j)
            for k in range(j,5))) for a in range(2) for b in range(2) for j in range(5)))
        check(name+'_weak_contact_symmetry',all(zero(weak[a][i][b][j]-weak[b][j][a][i])
            for a in range(2) for b in range(2) for i in range(3) for j in range(3)))
        A=core.original_A_coefficients(name)
        check(name+'_original_z4_conversion',all(zero(A[a][b][j]-z**4*H[a][b][j])
            for a in range(2) for b in range(2) for j in range(5)))
        coefficients[name]={'arbitrary_seed_jet_coefficients':20,
            'strong_coefficient_sha256':sha256(str(H).encode()).hexdigest()}
    reference=core.build('reference');qsub=core.seed_substitution()
    check('both_reference_seed_equations',all(zero(v.subs(qsub)) for v in reference['source']))
    # J is an exact primitive, not a sampled mass profile.
    J=sp.Symbol('J')
    def mass_D(e):return sp.cancel(core.D(e)+z**4/core.Q[0]**2*sp.diff(e,J))
    directions=[[-2*core.beta*core.Q[0]*J/z],
                [core.beta*z**3/core.Q[0]+core.beta*core.Q[1]*J]]
    for row in directions:
        for _ in range(4):row.append(mass_D(row[-1]))
    applied=core.apply_symbolic('reference',directions)
    for a in range(2):
        rational=sp.cancel(applied[a])
        # Preserve the exact finite-seed polynomial identities before division.
        numerator,denominator=sp.fraction(rational)
        check('uncut_reference_mass_row_'+str(a),zero(numerator.subs(qsub)))
    # For arbitrary cutoff jets, all orders 1..4 survive as an exact commutator.
    chi=sp.symbols('chi0:5')
    cutjets=[[sum(comb(k,j)*chi[j]*row[k-j] for j in range(k+1))
              for k in range(5)] for row in directions]
    cut=core.apply_symbolic('reference',cutjets)
    commutator=[sum(reference['strong'][a][b][k]*sum(
        comb(k,j)*chi[j]*directions[b][k-j] for j in range(1,k+1))
        for b in range(2) for k in range(5)) for a in range(2)]
    check('both_arbitrary_cutoff_commutators',all(zero(cut[a]-chi[0]*applied[a]-commutator[a]) for a in range(2)))
    check('cutoff_transition_not_zero',all(not zero(v) for v in commutator))
    # Areal return, with F=C Q/z^2 and x=C beta; r derivative is -D.
    n,b=directions[0][0],directions[1][0];F=core.C*core.Q[0]/z**2
    alpha=n+2*b+directions[0][1]
    f=2*F*(-directions[0][1]-b)+core.D(F)*n
    check('uncut_areal_mass_return',zero(alpha) and zero(f-2*core.C*core.beta*z))
    alphac=chi[0]*alpha+chi[1]*n
    fc=chi[0]*f-2*F*n*chi[1]
    check('original_cutoff_areal_rows',zero(alphac-n*chi[1]) and zero(fc-(2*core.C*core.beta*z*chi[0]-2*F*n*chi[1])))
    check('cut_lift_regular_bolt_vanishes',all(v.subs({s:0 for s in chi})==0 for row in cutjets for v in row))
    # Boundary first variations are retained, not asserted zero.
    check('local_source_boundary_rows_retained',all(
        len(core.build(name)['boundary_source'])==2 for name in core.densities))
    check('weyl_boundary_representative_retained',core.weyl_endpoint!=0)
    # These are CONDITIONAL bookkeeping consequences, never actual columns.
    allowance=arb('5.34235e-12');reserve=arb('.0012069')
    check('conditional_fast_screen',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    check('original_mass_conversion_debit',sp.Rational('1.627532e-15')*2174==sp.Rational('3.538254568e-12'))
    ceilings={'A_o_if_other_debits_zero':str(allowance/arb('4.084e-22')),
     'A_o_J_E_if_other_debits_zero':str(reserve/arb('3.712e-12')),
     'gamma_E_if_other_debits_zero':str(reserve/arb('.0195239'))}
    return {'precision_bits':ctx.prec,'generic_local_operators':coefficients,
     'exact_reference_mass_rows':2,'original_cutoff_jets_retained':[0,1,2,3,4],
     'conditional_ceilings_not_actual_constants':ceilings,
     'scope':{'exact_signed_local_operator_executed':True,
      'both_original_local_rows':True,'finite_seed_and_original_dz_retained':True,
      'strong_quantum_operator_application_executed':False,
      'complete_coupled_return_executed':False,'actual_kappa_m_computed':False,
      'residual_budgets_certified':False,'physical_X_comparisons':False,
      'nonlinear_admission_contraction':False,'exact_compatibility':False,
      'physical_matching':False,'semiclassical_validity_noise':False,
      'all_sector_stability':False,'actual_controlled_solution':False,
      'construction_obstruction':False,'full_goal_complete':False}}
if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
     'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
     'runtime':{'python':sys.version,'executable_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
                'flint_version':flint.__version__,'sympy_version':sp.__version__}}
    try:
        run['result']=execute();run['input_sha256_after']={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
