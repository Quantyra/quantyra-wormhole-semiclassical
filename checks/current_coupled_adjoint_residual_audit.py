"""Equation-index algebra and residual budgets; no quantum inverse executed."""
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
STEM='current-coupled-adjoint-residual-audit'
OWN=['notes/'+STEM+'.md','checks/current_coupled_adjoint_residual_audit.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+['notes/current-direct-triple-return-attempt.md',
 'evidence/current-direct-triple-return-attempt.json',
 'notes/current-geodesic-two-sided-remainder.md',
 'notes/stationary-curved-quantum-tt.md',
 'notes/current-seed-full-hessian-realization.md']
CHECKS=[]
def check(name,ok):
    CHECKS.append({'name':name,'passed':bool(ok)})
    print(name+(': passed' if ok else ': FAILED'),flush=True)
    if not ok: raise AssertionError(name)
def enc(x):
    m,r,e=x.mid_rad_10exp(50)
    return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('frozen_argument_and_checker',all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in review['scientific_sha256'].items()))
    prev=json.loads((ROOT/INPUTS[4]).read_bytes())['runs'][-1]
    check('accepted_composition_input',prev['status']=='passed' and not prev['result']['scope']['direct_composite_norms_enclosed'])
    # These indices are three equations, never physical modes.
    H=sp.Matrix([[1,-sp.Rational(1,2),0],[-sp.Rational(1,2),1,-sp.Rational(1,2)],[0,-sp.Rational(1,2),1]])
    check('equation_index_positive_spectrum',set(H.eigenvals())=={sp.Integer(1),1-sp.sqrt(2)/2,1+sp.sqrt(2)/2})
    t=sp.symbols('t',nonnegative=True)
    A=sp.Matrix([[1+t,0,0],[-t,1+t,0],[0,-t,1+t]])
    check('full_operator_coercivity_identity',(A+A.T)/2==sp.eye(3)+t*H)
    p=A.T.inv()*sp.Matrix([0,0,1])
    check('adjoint_norm_representation',sp.cancel(t*p[0]-(t/(1+t))**3)==0)
    q=arb(169)/500; f=q/(1+q)
    x0=1+2/q;c0,c1=arb(1),x0
    for j in range(1,29):c0,c1=c1,2*x0*c1-c0
    h2=(1+q)**2/c1;h3=(1+q)**3/c1
    pref=arb('1.1e-10')*h3;star=arb('1.667e-9')/pref
    a=star/2;tau=star/4;sigma=star/(4*q)
    check('original_physical_allocation',(pref*(a+tau+q*sigma)).contains(arb('1.667e-9')))
    check('residual_debit_half_bolt_budget',(pref*(tau+q*sigma)).contains(arb('1.667e-9')/2))
    check('equation_index_coercivity_positive',1-1/arb(2).sqrt()>0)
    result={'precision_bits':512,'complete_triple_residual_inverse_upper':1,
      'coercivity_T_margin':enc(1-1/arb(2).sqrt()),'q':enc(q),'f':enc(f),
      'h2':enc(h2),'h3':enc(h3),'affine_prefactor':enc(pref),
      'conditional_residual_targets':{'bolt_Cstar_full_margin':enc(star),
       'trial_energy_half_margin':enc(a),'riesz_residual_quarter_margin':enc(tau),
       'triple_adjoint_residual_quarter_margin':enc(sigma)},
      'scope':{'complete_triple_residual_stability_proved':True,
       'actual_trial_constructed':False,'quantum_residual_evaluated':False,
       'direct_composite_norms_enclosed':False,'original_strong_delta_certified':False,
       'actual_kappa_m_computed':False,'signed_coupled_columns_computed':False,
       'cold_contact_improvement_proved':False,'actual_controlled_solution':False,
       'construction_obstruction':False,'nonlinear_admission_contraction':False,
       'exact_compatibility':False,'physical_matching':False,'semiclassical_validity_noise':False,
       'all_sector_stability':False,'full_goal_complete':False}}
    return result
if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,
      'input_sha256_before':hashes,'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
      'runtime':{'python':sys.version,'python_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
                 'flint_version':flint.__version__,'sympy_version':sp.__version__}}
    try:
        run['result']=execute()
        run['input_sha256_after']={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
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
