"""Original-norm composition of a specified estimator; not a coupled solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import factorial
import json, sys, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb, ctx
import sympy as sp
ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-original-admission-map-budget-test'
OWN = ['notes/'+STEM+'.md', 'checks/current_original_admission_map_budget_test.py',
       'evidence/'+STEM+'-review.json']
INPUTS = OWN + ['notes/'+p+'.md' for p in [
    'current-native-output-structured-source-test', 'current-interior-translation-jet-bound',
    'current-direct-triple-return-attempt', 'current-schur-residual-phase-domain',
    'current-symmetric-strong-response-attempt', 'divided-jet-dirac-response',
    'comparison-path-interior-response', 'finite-cold-admission-and-response',
    'finite-outer-admission-and-response', 'stationary-completed-metric-response-norm',
    'current-controlled-solution-audit-closeout']] + [
    'evidence/current-native-output-structured-source-test.json',
    'evidence/current-localized-bolt-adjoint-sensitivity.json']
CHECKS = []
def sha(p): return sha256(p.read_bytes()).hexdigest()
def check(name, ok):
    CHECKS.append({'name': name, 'passed': bool(ok)})
    print(name+(': passed' if ok else ': FAILED'), flush=True)
    if not ok: raise AssertionError(name)
def enc(x):
    m, r, e = x.mid_rad_10exp(50)
    return {'ball': x.str(35), 'decimal_enclosure': {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}
def dec(x):
    if isinstance(x, dict):
        d = x['decimal_enclosure']
        return arb(d['mid']+' +/- '+d['rad'])*arb(10)**int(d['exponent'])
    return arb(x.strip('[]'))
def log10(x): return x.log()/arb(10).log()
def execute():
    ctx.prec = 512
    review = json.loads((ROOT/OWN[2]).read_bytes())
    check('ordinary_argument_frozen', all(sha(ROOT/p)==h for p,h in review['scientific_sha256'].items()))
    prior = json.loads((ROOT/INPUTS[-2]).read_bytes())['runs'][-1]
    pars = json.loads((ROOT/INPUTS[-1]).read_bytes())['runs'][-1]
    check('native_and_parameter_receipts_accepted', prior['status']==pars['status']=='passed')
    h = dec(pars['result']['parameters']['h']); ellb = (1+h)/h
    q = arb('.338'); f = q/(1+q); xx = -1-2/q
    ch0, ch1 = arb(1), xx
    for _ in range(2,30): ch0,ch1=ch1,2*xx*ch1-ch0
    d29 = abs(ch1); h2=(1+q)**2/d29; h3=(1+q)**3/d29
    paff=arb('1.1e-10')*h3
    check('same_original_affine_prefactor', paff.overlaps(dec(prior['result']['budget']['P32_affine_prefactor'])))
    coeffs=[(arb(factorial(j)*factorial(5-j))/(4*arb.pi()*factorial(6))).sqrt() for j in range(6)]
    c0=coeffs[0]
    check('unscaled_gradient_embedding', coeffs[1]<c0)
    check('divided_high_jet_embedding', all(coeffs[j]/(factorial(j)*100**j)<c0 for j in range(2,6)))
    check('exact_beta_integrals', all(sp.factorial(j)*sp.factorial(5-j)/sp.factorial(6)==sp.integrate(sp.Symbol('t')**j*(1-sp.Symbol('t'))**(5-j),(sp.Symbol('t'),0,1)) for j in range(6)))
    # No unit conversion is omitted: this is a term of the specified positive estimator.
    K=arb(4)*arb(10)**4003; L=(2*K)**(arb(4)/3)
    M3=K*L**(arb(25)/4); M10=K*L**(arb(53)/4)
    check('native_absorption_identity', M3.overlaps(arb(2)**(arb(25)/3)*K**(arb(28)/3)))
    check('stronger_auxiliary_absorption_identity', M10.overlaps(arb(2)**(arb(53)/3)*K**(arb(56)/3)))
    check('half_absorption', (K/L**(arb(3)/4)).overlaps(arb('.5')))
    lowpref=2*c0*arb('1e19')*f**3*paff
    ap5=lowpref*M3; abolt=ellb*ap5
    Msym_same=8*K**4
    rows={}
    for name,margin,conversion in [('P5_whole',arb('.001'),arb(1)),
                                   ('P5_reference_margin',arb('.00087'),arb(1)),
                                   ('bolt_whole_shared_pair',arb('1.667e-9'),ellb)]:
        target=margin/(lowpref*conversion)
        bmax=(target/arb(2)**(arb(25)/3))**(arb(3)/28)/4000
        bsym=(target/8)**(arb(1)/4)/4000
        ratio=lowpref*conversion*M3/margin
        ratiosym=lowpref*conversion*Msym_same/margin
        check(name+'_generic_estimator_fails',ratio>1)
        check(name+'_order_gain_alone_fails',ratiosym>1)
        check(name+'_generic_threshold_identity',(arb(2)**(arb(25)/3)*(4000*bmax)**(arb(28)/3)).overlaps(target))
        check(name+'_symmetric_threshold_identity',(8*(4000*bsym)**4).overlaps(target))
        rows[name]={k:enc(v) for k,v in {
            'margin':margin,'original_output_conversion':conversion,'M_target_low_term_only':target,
            'log10_generic_excess':log10(ratio),'log10_order_gain_only_excess':log10(ratiosym),
            'B_H6_25_to_H3_ceiling':bmax,'Bsym_H5_25_to_H3_ceiling':bsym,
            'log10_required_generic_coefficient_improvement':4000-log10(bmax),
            'log10_required_symmetric_coefficient_improvement':4000-log10(bsym)}.items()}
    check('generic_failure_remains_enormous_in_native_norm',log10(ap5/arb('.001'))>37000)
    # The cutoff cubic has jump 12 in its third derivative, hence 12 delta in the fourth.
    y=sp.Symbol('y',positive=True); eps,rad=sp.symbols('eps rad',positive=True)
    jump=sp.diff(y**3,y,3)-sp.diff(-y**3,y,3)
    check('actual_mass_seam_delta_coefficient',jump==12)
    check('actual_log_squared_one_sided_derivative',sp.simplify(sp.diff(y**4*sp.log(y)**2,y,4)-(24*sp.log(y)**2+100*sp.log(y)+70))==0)
    moment=sp.integrate(1/y,(y,eps,rad))
    check('restored_second_moment_positive_envelope',sp.simplify(moment-sp.log(rad/eps))==0)
    check('unrestored_envelope_diverges',sp.limit(moment,eps,0,dir='+')==sp.oo)
    check('cubic_fourier_sobolev_endpoint',sp.solve_univariate_inequality(2*sp.Symbol('s',real=True)-8<-1,sp.Symbol('s',real=True))==(sp.Symbol('s',real=True)<sp.Rational(7,2)))
    defect=arb('.5')/h3; graph=arb('5.34235e-12')/paff
    check('form_defect_retained',q**3/d29<arb('7.351e-35'))
    check('quantum_screen_limited_value',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    check('exact_structured_source_c4_only',prior['result']['complete_structured_source']['quartic_source_finitepart']==0)
    return {'precision_bits':512,'original_constants':{k:enc(v) for k,v in {
        'D29':d29,'h2':h2,'h3':h3,'p_aff':paff,'f':f,'bolt_radius_conversion':ellb,
        'P5_local_H7_embedding':c0,'log10_M3':log10(M3),'log10_M10':log10(M10),
        'log10_affine_P5_generic_debit':log10(ap5),'log10_affine_bolt_generic_debit':log10(abolt),
        'conditional_original_defect_product_ceiling_no_error':defect,
        'conditional_graph_ceiling_whole_margin':graph}.items()},
        'original_output_budget_tests':rows,
        'domain_diagnostics':{'mass_seam_fourth_derivative_delta':12,'mass_global_Hs_strict_upper_endpoint':'7/2',
            'positive_second_moment':'log(r/epsilon)','complete_restored_moment_bound':None},
        'unknown_complete_maps':{k:None for k in ['Cdef3','I_Y','Cgraph3','original_X_completion','Cphysical3',
            'D_ell_embedding','outer_proper_identification','weighted_radius_output','actual_kappa_m',
            'signed_coupled_centers','certified_application_errors','nonlinear_neighborhood_comparison']},
        'scope':{'original_admission_ledger_recovered':True,'specified_generic_estimator_fails':True,
            'stronger_areal_estimator_required_for_admission':False,'symmetric_complete_gain_proved':False,
            'structured_source_c4_zero':True,'full_structured_source_norm_computed':False,
            'actual_polynomial_correction_applied':False,'complete_coupled_bound_certified':False,
            'original_defect_certified':False,'new_physical_gate_closed':False,'construction_obstruction':False,
            'degree_change_adopted':False}}
if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    out=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={p:sha(ROOT/p) for p in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
         'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
         'runtime':{'python':sys.version,'python_sha256':sha(Path(sys.executable)),
                    'flint_version':flint.__version__,'sympy_version':sp.__version__}}
    try:
        run['result']=execute();run['input_sha256_after']={p:sha(ROOT/p) for p in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};ledger['runs'].append(run)
    out.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    summary=json.dumps({'status':run['status'],'mode':mode,'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as stream:stream.write(summary.encode())
    print(summary,flush=True);sys.exit(0 if run['status']=='passed' else 1)
