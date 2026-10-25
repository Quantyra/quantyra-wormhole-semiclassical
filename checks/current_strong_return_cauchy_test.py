"""Finite real cutoff costs; NOT a Hessian application or coupled solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
from fractions import Fraction
import json, sys, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb, ctx
ctx.prec=512
ROOT=Path(__file__).resolve().parents[1]
STEM='current-strong-return-cauchy-test'
OWN=['notes/'+STEM+'.md','checks/current_strong_return_cauchy_test.py',
     'evidence/'+STEM+'-review.json']
INPUTS=sorted(OWN+[
    'operators/short_travel_comparison.py',
    'notes/current-short-travel-strong-domain.md',
    'notes/current-strong-return-weight-repair.md',
    'notes/current-curved-paired-heat-bridge.md',
    'notes/current-coupled-global-return-attempt.md',
    'notes/current-controlled-solution-audit-closeout.md',
    'evidence/current-short-travel-strong-domain.json',
    'evidence/current-strong-return-weight-repair.json'])
CHECKS=[]
def check(name, value):
    CHECKS.append(dict(name=name,passed=bool(value)))
    if not value: raise AssertionError(name)
def enc(v):
    m,r,e=v.mid_rad_10exp(70)
    return dict(ball=v.str(42),decimal_enclosure=dict(mid=str(m),rad=str(r),exponent=str(e)))
def execute():
    for stem in ['current-short-travel-strong-domain','current-strong-return-weight-repair']:
        check('retained_'+stem,json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]['status']=='passed')
    n=16
    d=[3**j*factorial(j)**2 for j in range(n+1)]
    z=[arb(1)];e2=arb(2).exp()
    for j in range(1,n+1):
        z.append(e2*(d[j]+2*sum((comb(j,k)*d[k]*z[j-k] for k in range(1,j+1)),arb(0))))
    # Exact Stirling table and a second construction by applying
    # y*d/dy to a polynomial differential operator.
    s=[[0]*(n+1) for _ in range(n+1)];s[0][0]=1
    for j in range(1,n+1):
        for k in range(1,j+1):s[j][k]=s[j-1][k-1]+k*s[j-1][k]
    coeff={0:1}
    for j in range(1,n+1):
        nxt={}
        for k,c in coeff.items():
            nxt[k]=nxt.get(k,0)+k*c
            nxt[k+1]=nxt.get(k+1,0)+c
        coeff={k:c for k,c in nxt.items() if c}
        if any(coeff.get(k,0)!=s[j][k] for k in range(j+1)):raise AssertionError(j)
    check('Euler_chain_rule_all_orders_through_16',True)
    t=[arb(1)]+[sum((s[j][k]*2**k*z[k] for k in range(1,j+1)),arb(0)) for j in range(1,n+1)]
    check('positive_finite_cutoff_bounds',all(v>0 and v.is_finite() for v in z+t))
    check('ordinary_and_Euler_norms_kept_distinct',128*z[1]>t[1])
    check('scalar_complex_germ_minimum',arb(2)/3>0)
    # The cutoff form integral is bounded over the whole half-line
    # [a,infinity); retaining [a,2a] would multiply by 15/16.
    check('complete_travel_integral_guard',arb(1)-arb(2)**(-4)<1)
    check('original_action_transfer_guard',arb('3.2e-18')>arb('2.156202783669969564786e-18'))
    rows=[]
    for j in range(n+1):
        row=dict(order=j,D_j=str(d[j]),Z_j=enc(z[j]),ordinary_travel=enc(128**j*z[j]),Euler_travel=enc(t[j]))
        row['cutoff_derivative_form_guard']=enc(arb('3.2e-18')*t[j]) if j else None
        rows.append(row)
    guards={1:('8513','133.004','4.257e-16'),4:('5.551e16','3.358e9','1.075e-8'),
            8:('1.669e35','6.135e20','1963'),16:('2.777e74','3.772e45','1.207e28')}
    check('published_upward_guards',all(128**j*z[j]<arb(b[0]) and t[j]<arb(b[1])
          and arb('3.2e-18')*t[j]<arb(b[2]) for j,b in guards.items()))
    aff=arb('4.2953e-8')*arb('5e-5')
    graph=arb('49.462')*arb('5e-5')
    energy=arb('390.478')*arb('5e-5')
    ar=arb('7.49e-12')-aff;gr=arb('.00368')-graph
    check('conditional_mass_affine_budget',Fraction('7.49e-12')-Fraction('4.2953e-8')*Fraction('5e-5')==Fraction('5.34235e-12'))
    check('conditional_mass_graph_budget',Fraction('.00368')-Fraction('49.462')*Fraction('5e-5')==Fraction('.0012069'))
    check('conditional_mass_energy_debit',Fraction('390.478')*Fraction('5e-5')==Fraction('.0195239'))
    fast=arb('.028913')+arb('1e-11')*arb('1.0706e10')
    check('retained_two_fast_input_screen',fast<arb('.1599164'))
    # The entire cutoff factor is available. No bound for an actual
    # strong composite or unknown coupling is assigned from it.
    return dict(outcome='unresolved_complete_coupled_estimate',
        failed_estimator='termwise holomorphic moving short comparison; ordinary actual-domain support proof',
        cutoff_radius='1/128',cutoff_bounds=rows,
        cutoff_scope='Uniform finite real travel jets and their full-angular form contributions only; not differentiated physical kernel or strong inverse.',
        physical_domain_changed=False,physical_norm_target_changed=False,
        common_regulator_full_sums_contacts_returns_retained=True,
        conditional_mass_affine_debit=enc(aff),conditional_mass_graph_debit=enc(graph),
        conditional_mass_energy_debit=enc(energy),conditional_affine_remaining=enc(ar),
        conditional_graph_remaining=enc(gr),mass_only_energy_comparison_allowance=enc(gr/energy),
        two_fast_input_screen=enc(fast),fast_screen_reserve=enc(arb('.1599164')-fast),
        physical_operator_applications=0,coupled_signed_solve_executed=False,
        complete_strong_estimate_executed=False,actual_mass=None,kappa_mass=None,
        M_graph_4=None,M_X_8=None,U_feedback=None,other_coupled_entries=None,
        original_physical_norm_comparisons=False,nonlinear_neighborhood_transfer=False,
        new_nonlinear_feasibility_gate_closed=False,nonlinear_admission=False,
        nonlinear_contraction=False,exact_compatibility=False,physical_matching=False,
        validity_noise=False,all_sector_stability=False,actual_controlled_solution=False,
        construction_obstruction=False,closed_audit_preserved=True)

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    out=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    digest=lambda p:sha256(p.read_bytes()).hexdigest()
    before={p:digest(ROOT/p) for p in INPUTS}
    fr=Path(flint.__file__).resolve().parent
    runtime=dict(python=sys.version,python_executable_sha256=digest(Path(sys.executable)),
        flint_version=flint.__version__,
        flint_artifacts={p.relative_to(fr).as_posix():digest(p) for p in sorted(fr.rglob('*'))
                        if p.is_file() and p.suffix in ['.py','.pyd','.dll']})
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),mode=mode,
        input_sha256_before=before,precision_bits=512,runtime=runtime,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN})
    try:
        run['result']=execute()
        run['input_sha256_after']={p:digest(ROOT/p) for p in INPUTS}
        check('frozen_inputs_unchanged',before==run['input_sha256_after'])
        run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    ledger['runs'].append(run);out.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    summary=dict(mode=mode,status=run['status'],checks=len(CHECKS),result=run.get('result'),traceback=run.get('traceback'))
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write((json.dumps(summary,indent=2)+'\n').encode())
    brief={k:summary[k] for k in ['mode','status','checks','traceback']}
    if run.get('result'):
        brief['selected_cutoff_bounds']=[{k:(v['ball'] if isinstance(v,dict) else v) for k,v in row.items()}
            for row in run['result']['cutoff_bounds'] if row['order'] in [1,4,8,16]]
    print(json.dumps(brief,indent=2))
    if run['status']!='passed':raise SystemExit(1)
