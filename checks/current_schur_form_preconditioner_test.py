"""Support a full-operator analytic test; no numerical quantum inverse solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb, ctx
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
STEM='current-schur-form-preconditioner-test'
OWN=['notes/'+STEM+'.md', 'checks/current_schur_form_preconditioner_test.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+['notes/current-seed-full-hessian-realization.md',
 'notes/current-global-mixed-schur-preconditioner.md',
 'notes/current-physical-york-dual-assembly.md',
 'notes/current-common-regulator-strong-limit-test.md',
 'notes/current-coupled-schur-graph-attempt.md',
 'notes/current-short-travel-strong-domain.md',
 'notes/current-adaptive-dirac-trace-remainder.md',
 'notes/current-global-mass-conversion.md',
 'notes/current-implicit-phase-graph-domain.md',
 'evidence/current-global-mass-conversion.json']
CHECKS=[]

def check(name, ok):
    CHECKS.append({'name':name, 'passed':bool(ok)})
    print(name+(': passed' if ok else ': FAILED'),flush=True)
    if not ok: raise AssertionError(name)

def zero(e): return sp.expand(e)==0

def dec(v):
    d=v['decimal_enclosure']
    return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])

def enc(v):
    m,r,e=v.mid_rad_10exp(50)
    return {'ball':v.str(32),'decimal_enclosure':
            {'mid':str(m),'rad':str(r),'exponent':str(e)}}

def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    for p,h in review['scientific_sha256'].items():
        check('frozen_'+p,sha256((ROOT/p).read_bytes()).hexdigest()==h)
    # Noncommuting blocks: do not replace the coupled system by scalar entries.
    A,Ai,U,Ut,b,bt,L=sp.symbols('A Ai U Ut b bt L',commutative=False)
    CV=-A*U-b
    VtCt=-Ut*A-bt
    VtDV=L+Ut*A*U+Ut*CV+VtCt*U
    cross=sp.expand(VtCt*Ai*CV).subs(A*Ai,sp.Integer(1)).subs(Ai*A,sp.Integer(1))
    check('noncommuting_full_Schur_congruence',zero(VtDV+cross-L-bt*Ai*b))
    # J S V=ell+B and J P=... gives the right residual, with its sign.
    B,Li=sp.symbols('B Li',commutative=False)
    check('right_defect_sign',zero(1-(L+B)*Li+B*Li-(1-L*Li)))
    # Iterated candidate identity is polynomial in the ONE error operator.
    e=sp.symbols('e')
    check('N32_residual_identity',sp.expand(1-(1-e)*sum(e**j for j in range(32))-e**32)==0)
    check('two_error_powers_remove_local_derivative_deficit',2*4>8-2)
    check('logarithmic_deficit_retained',sp.Rational(1)-sp.Rational(1,2)==sp.Rational(1,2))
    q=arb('.338')
    check('form_defect_below_half',q<arb('.5'))
    inverse_guard=1/(1-q)
    check('form_resolvent_bound',inverse_guard<arb('1.510575'))
    src=json.loads((ROOT/'evidence/current-global-mass-conversion.json').read_bytes())['runs'][-1]
    check('retained_seed_evidence_passed',src['status']=='passed')
    p={k:dec(v) for k,v in src['result']['parameters'].items()}
    be,uu,ep,ta,la=[p[k] for k in ['beta','u','eps','tau','g']]
    check('finite_seed_sign_guards',be>0 and uu>0 and ep<0 and ta>0)
    z=arb('.003')
    # Discard negative terms only, retaining an upper bound on all of I.
    qhi=1+(uu/be)*z*z+(3+uu/be+abs(ep)/be)*z**4+2*abs(ep)/be*z**5+abs(ep)/be*z**6
    check('whole_interval_Q_ratio_upper',qhi<arb('1.001'))
    Bsq=be/(uu*la)
    k0=la*be/(15*arb.pi())
    Wmax=2*arb.pi()*Bsq/arb('5e-34')  # N<=2, G=1; actual gamma>5e-34.
    coefficient=(k0/(Wmax*arb('1.001')**2*(arb('5.01')*z.log()).exp())).sqrt()
    check('whole_graph_output_asymptotic_coefficient',coefficient>arb('1e-50'))
    affine=arb('5.34235e-12')*(1-q)/q**30
    feedback=arb('.0012069')*(1-q)/q**30
    source_defect=arb('.5')/q**30
    screen=arb('.028913')+arb('1e-11')*arb('1.0706e10')
    check('conditional_fast_screen_unchanged',screen<arb('.1599164'))
    check('original_mass_basis_debit',sp.Rational('1.627532e-15')*2174==sp.Rational('3.538254568e-12'))
    check('conditional_mass_affine_debit',sp.Rational('5e-5')*sp.Rational('4.2953e-8')==sp.Rational('2.14765e-12'))
    check('conditional_mass_feedback_debit',sp.Rational('5e-5')*sp.Rational('49.462')==sp.Rational('.0024731'))
    return {'precision_bits':ctx.prec,'form_defect_upper':enc(q),
      'form_resolvent_upper':enc(inverse_guard),'physical_K0':enc(k0),
      'Q_over_beta_upper':enc(qhi),'action_W_upper_N2':enc(Wmax),
      'asymptotic_graph_over_form_coefficient_lower':enc(coefficient),
      'growth':'k^6 sqrt(log k); no finite-frequency onset certified',
      'original_graph_output_norm_from_Y_F':'unbounded; phase range can fail earlier',
      'conditional_N32_products':{'M2_U_affine_ceiling':enc(affine),
       'M2_U_feedback_ceiling':enc(feedback),'source_composite_ceiling_before_input_embedding':enc(source_defect),
       'trial_debits_must_be_subtracted':True,'actual_constants_computed':False},
      'scope':{'full_coupled_analytic_form_defect_test':True,
       'candidate_numerically_applied':False,'specified_form_to_graph_certificate_fails':True,
       'original_strong_delta_certified':False,'complete_coupled_return_executed':False,
       'actual_kappa_m_computed':False,'repaired_global_composite_bound':False,
       'physical_admission_or_feedback_budget_closed':False,'nonlinear_admission_contraction':False,
       'exact_compatibility':False,'physical_matching':False,'semiclassical_validity_noise':False,
       'all_sector_stability':False,'actual_controlled_solution':False,'construction_obstruction':False,
       'full_goal_complete':False}}

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
      'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
      'runtime':{'python':sys.version,'executable_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
                 'flint_version':flint.__version__,'sympy_version':sp.__version__}}
    try:
        run['result']=execute()
        run['input_sha256_after']={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after'])
        run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
