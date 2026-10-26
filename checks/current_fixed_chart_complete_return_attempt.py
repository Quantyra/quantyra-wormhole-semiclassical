"""Vertex algebra and estimator sensitivity only; no quantum inverse is applied."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path') + 1])
import sympy as sp
from flint import arb, ctx
ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-fixed-chart-complete-return-attempt'
OWN = ['notes/' + STEM + '.md', 'checks/current_fixed_chart_complete_return_attempt.py',
       'evidence/' + STEM + '-review.json']
INPUTS = OWN + ['notes/' + n + '.md' for n in [
    'current-original-admission-map-budget-test', 'current-native-output-structured-source-test',
    'current-interior-translation-jet-bound', 'current-curved-paired-heat-bridge',
    'current-rough-paired-variation', 'current-common-regulator-conversion-test',
    'current-symmetric-strong-response-attempt', 'current-schur-residual-phase-domain',
    'current-controlled-solution-audit-closeout']] + [
    'evidence/current-original-admission-map-budget-test.json',
    'evidence/current-coupled-mass-adjoint-runtime.json']
CHECKS = []
def sha(p): return sha256(p.read_bytes()).hexdigest()
def check(name, condition):
    ok = bool(condition)
    CHECKS.append({'name': name, 'passed': ok})
    print(name + (': passed' if ok else ': FAILED'), flush=True)
    if not ok: raise AssertionError(name)
def enc(x):
    m, r, e = x.mid_rad_10exp(50)
    return {'ball': x.str(35), 'decimal_enclosure': {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}
def dec(x):
    d = x['decimal_enclosure']
    return arb(d['mid'] + ' +/- ' + d['rad']) * arb(10)**int(d['exponent'])
def execute():
    ctx.prec = 512
    review = json.loads((ROOT / OWN[2]).read_bytes())
    check('frozen_ordinary_sources', all(sha(ROOT / p) == h for p, h in review['scientific_sha256'].items()))
    rt = json.loads((ROOT / INPUTS[-1]).read_bytes())
    check('pinned_python', sha(Path(sys.executable)) == rt['python_sha256'])
    check('pinned_packages', all(sha(Path(pkg['root']) / p) == h
          for pkg in rt['packages'].values() for p, h in pkg['files'].items()))
    # Differential expressions acting on an arbitrary field, not a finite mode model.
    Y, s, t = sp.symbols('Y s t', real=True)
    a, b, u, v, w, zeta, psi = [sp.Function(n)(Y) for n in ['a','b','u','v','w','zeta','psi']]
    ar, br = a + s*u + t*w, b + s*v + t*zeta
    radial = sp.exp(-ar)*(sp.diff(psi,Y) - sp.diff(ar,Y)*psi/2)
    temporal = sp.exp(-br)  # f0 partial_tau is independent of the parameters.
    at0 = lambda expr: expr.subs({s: 0, t: 0})
    first = sp.exp(-a)*(-u*(sp.diff(psi,Y)-sp.diff(a,Y)*psi/2)-sp.diff(u,Y)*psi/2)
    second = sp.exp(-a)*(u*w*(sp.diff(psi,Y)-sp.diff(a,Y)*psi/2)+(u*sp.diff(w,Y)+w*sp.diff(u,Y))*psi/2)
    check('radial_first_vertex', sp.simplify(at0(sp.diff(radial,s))-first) == 0)
    check('radial_second_vertex', sp.simplify(at0(sp.diff(radial,s,t))-second) == 0)
    check('thermal_vertices', sp.simplify(at0(sp.diff(temporal,s))+v*sp.exp(-b)) == 0
          and sp.simplify(at0(sp.diff(temporal,s,t))-v*zeta*sp.exp(-b)) == 0)
    conf = sp.exp(-u/2)*sp.diff(sp.exp(-u/2)*psi,Y)
    check('conformal_endpoint_half_densities', sp.simplify(conf-sp.exp(-u)*(sp.diff(psi,Y)-sp.diff(u,Y)*psi/2)) == 0)
    q = arb('.338'); f = q/(1+q); x = -1-2/q
    t0, t1 = arb(1), x
    for j in range(2,30): t0, t1 = t1, 2*x*t1-t0
    D29 = abs(t1); h2 = (1+q)**2/D29; h3 = (1+q)**3/D29
    paff = arb('1.1e-10')*h3; c0 = 1/(24*arb.pi()).sqrt()
    log10 = lambda a: a.log()/arb(10).log()
    old = json.loads((ROOT / INPUTS[-2]).read_bytes())['runs'][-1]
    check('accepted_original_ledger', old['status'] == 'passed')
    check('unchanged_polynomial_and_output_units', all(v.overlaps(dec(old['result']['original_constants'][k]))
          for k,v in [('D29',D29),('h2',h2),('h3',h3),('p_aff',paff),('f',f),('P5_local_H7_embedding',c0)]))
    low_prefactor = 2*c0*paff*arb('1e19')*f**3
    logK = log10(arb(4000)) + 300
    logM = arb(25)/3*log10(arb(2)) + arb(28)/3*logK
    logMsym = log10(arb(8)) + 4*logK
    results = {}
    for name, margin in [('whole',arb('.001')),('reference',arb('.00087'))]:
        target = margin/low_prefactor
        B = (target/arb(2)**(arb(25)/3))**(arb(3)/28)/4000
        Bsym = (target/8)**(arb(1)/4)/4000
        excess = log10(low_prefactor/margin)+logM
        symexcess = log10(low_prefactor/margin)+logMsym
        check(name+'_jet_cost_removed_estimator_still_excessive', excess > 2800 and symexcess > 1190)
        results[name] = {k:enc(v) for k,v in {
            'log10_excess_if_H3_coefficient_granted_1e300':excess,
            'log10_excess_if_complete_order_gain_also_granted':symexcess,
            'required_B_13over4_low_term_only':B,
            'required_B_9over4_low_term_only':Bsym,
            'log10_coefficient_gain_from_1e300':300-log10(B),
            'log10_coefficient_gain_from_1e300_with_order_gain':300-log10(Bsym)}.items()}
    check('full_absorption_condition', arb(4000)/8000 == arb('.5'))
    screen = arb('.028913')+arb('1e-11')*arb('1.0706e10')
    check('limited_quantum_fast_screen_only', screen < arb('.1599164'))
    unknown = {k:None for k in ['strong_common_regulator_limit_constant','Cdef3','I_Y','Cg3',
        'Cphysical3','cold_embedding','global_radius_embedding','signed_coupled_centers',
        'independent_phase_comparisons','actual_kappa_m','application_errors','neighborhood_comparison']}
    return {'precision_bits':512,'new_computation':'fixed-chart differential vertices and estimator sensitivity',
        'complete_quantum_map_numerically_evaluated':False,
        'sensitivity_hypotheses':['grant H3 coefficient 1e300 without any jet-transfer cost',
            'grant principal inverse factor 4000/C','set all other positive debits to zero'],
        'sensitivity_results':results,
        'conditional_budget_values':{k:enc(v) for k,v in {
            'h2':h2,'h3':h3,'p_aff':paff,'original_defect_ceiling_if_e_def_zero':arb('.5')/h3,
            'graph_ceiling_if_whole_margin_free':arb('5.34235e-12')/paff,
            'two_fast_screen_in_its_proved_scope':screen,
            'B_sufficient_for_full_absorption_of_low_term':arb(1)/8000}.items()},
        'unknown_complete_maps':unknown,
        'scope':{'actual_domain_preserved':True,'full_sums_and_contacts_retained_in_formula':True,
            'strong_common_regulator_limit_proved':False,'new_consequential_coupled_bound':False,
            'complete_coupled_bound_certified':False,'original_defect_certified':False,
            'structured_map_disproved':False,'new_physical_gate_closed':False,
            'construction_obstruction':False,'degree_change_adopted':False,'full_goal_complete':False}}
if __name__ == '__main__':
    mode = 'acceptance' if '--acceptance' in sys.argv else 'development'
    out = ROOT/'evidence'/(STEM+('' if mode == 'acceptance' else '-development')+'.json')
    hashes = {p:sha(ROOT/p) for p in INPUTS}
    run = {'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
           'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN}}
    try:
        run['result'] = execute()
        run['input_sha256_after'] = {p:sha(ROOT/p) for p in INPUTS}
        check('inputs_unchanged', hashes == run['input_sha256_after'])
        run['status'] = 'passed'
    except Exception:
        run['status'] = 'failed'; run['traceback'] = traceback.format_exc()
    run['checks'] = CHECKS; run['finished_utc'] = datetime.now(timezone.utc).isoformat()
    ledger = json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    ledger['runs'].append(run); out.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    summary = json.dumps({'status':run['status'],'mode':mode,'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as stream: stream.write(summary.encode())
    print(summary,flush=True); sys.exit(0 if run['status']=='passed' else 1)
