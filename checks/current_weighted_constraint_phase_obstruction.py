"""Supporting arithmetic for the constrained weighted global-return attempt; no complete Hessian bound."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys
import traceback

if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path') + 1])
from flint import arb, ctx
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-weighted-constraint-phase-obstruction'
OWN = ['notes/current-weighted-constraint-phase-obstruction.md', 'checks/current_weighted_constraint_phase_obstruction.py', 'evidence/current-weighted-constraint-phase-obstruction-review.json']
INPUTS = ['notes/current-weighted-constraint-phase-obstruction.md', 'checks/current_weighted_constraint_phase_obstruction.py', 'evidence/current-weighted-constraint-phase-obstruction-review.json', 'evidence/current-weighted-constraint-phase-obstruction-exploration.json', 'notes/current-restored-weight-commutator.md', 'evidence/current-restored-weight-commutator.json', 'notes/current-global-boundary-moment-return.md', 'notes/current-physical-york-dual-assembly.md', 'notes/current-lifted-boundary-graph-attempt.md', 'notes/current-coupled-global-return-attempt.md', 'notes/current-controlled-solution-audit-closeout.md', 'evidence/current-bolt-adjoint-domain.json', 'checks/current_all_angular_principal_operator.py', 'checks/current_coincident_dirac_action_assembly.py', 'checks/current_endpoint_thermal_restoration.py', 'checks/current_physical_york_dual_assembly.py', 'checks/current_seed_radial_green_enclosure.py', 'checks/current_seed_separated_quantum_response.py', 'checks/current_subtracted_thermal_action.py']
CHECKS = []


def check(name, condition, detail=""):
    CHECKS.append({'name': name, 'passed': bool(condition)})
    if not condition:
        raise AssertionError(name)


def encoded(value):
    m, r, e = value.mid_rad_10exp(60)
    return {'ball': value.str(45), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def decode(value):
    p = value['decimal_enclosure']
    return arb(int(p['mid']), int(p['rad'])) * arb(10) ** int(p['exponent'])


from math import comb, factorial
sys.path.insert(0, str(ROOT/'checks'))
# A legacy parameter helper parses argv on import. Keep this wrapper's
# acceptance mode local while importing it, then restore the exact argv.
original_argv = sys.argv[:]
sys.argv = [v for v in sys.argv if v != '--acceptance']
import current_physical_york_dual_assembly as york
sys.argv = original_argv
york.check = check
york.prior.check = check

def execute():
    ctx.prec=512;a=arb;enc=encoded
    p,h,zh,env=york.actual_parameters();be=p['beta'];sig=a('4.005');s=a('1.4')
    K=p['g']*be/(15*a.pi());B4=-24*K
    check('actual_seed_endpoint_guards',h>0 and p['tau']>0 and p['eps']<0 and h<a('1e-30'))
    check('negative_Qprime_on_bump_interval',12*be*a('.25')**2*a('.75')-2*abs(p['u'])-4*abs(p['eps'])>a('.09'))
    check('entire_finite_interval_Qprime_majorant',12*be+2*abs(p['u'])+4*abs(p['eps'])+3*p['tau']/(1-zh)<4)
    check('thermal_over_actual_bolt_gap',3*p['tau']/(1-zh)<a('1e-30'))
    gap=((h+a('.8'))/(h+a('.6')))**s-((h+a('.8'))/(h+a('.7')))**s
    check('actual_domain_phase_gap',gap>a('.29038048278'))
    check('full_normal_quartic_nonzero',B4<0 and abs(B4)>a('8.2893199525e-35'))
    leadg=4*p['g']*gap/(5*a.pi())
    check('induced_fast_source_nonzero',leadg>a('7.22117018e-35'))
    check('fixed_fast_phase_nonintegrable',2*(4-sig)<0)
    z=sp.symbols('z',positive=True);Q=sp.Function('Q')(z);A=sp.Function('A')(z);C=sp.symbols('C')
    S=Q*A/z;T=-(z**4/Q+sp.diff(Q,z)*A)/2;F=C*Q/z**2
    subs={sp.diff(A,z):z**4/Q**2}
    ff=(2*F*(-z*sp.diff(S,z)-T)+z*sp.diff(F,z)*S).subs(subs)
    aa=(S+2*T+z*sp.diff(S,z)).subs(subs)
    check('exact_quartic_local_areal_f',sp.simplify(ff+C*z*z)==0)
    check('exact_quartic_local_areal_alpha',sp.simplify(aa)==0)
    chi=sp.Function('chi')(z)
    fcut=(2*F*(-z*sp.diff(chi*S,z)-chi*T)+z*sp.diff(F,z)*chi*S).subs(subs)
    acut=(chi*S+2*chi*T+z*sp.diff(chi*S,z)).subs(subs)
    check('retained_cutoff_f_transition',sp.simplify(fcut+chi*C*z*z+2*F*S*z*sp.diff(chi,z))==0)
    check('retained_cutoff_alpha_transition',sp.simplify(acut-S*z*sp.diff(chi,z))==0)
    check('original_scalar_quartic_coordinate',sp.Rational(1,5)-sp.Rational(7,15)/2==-sp.Rational(1,30))
    # Exact weighted constraint difference, before any coefficient estimate.
    W=sp.Function('W')(z);br=sp.Function('br')(z);xi=sp.Function('xi')(z)
    constraint=lambda b,x:sp.diff(Q*b/z**4,z)+sp.diff(Q,z)*x/z**4
    check('weighted_TT_constraint_identity',sp.simplify(constraint(W*br,W*xi)-W*constraint(br,xi)-sp.diff(W,z)*Q*br/z**4)==0)
    # Whole analytic profile with a certified infinite coefficient tail.
    q=[a(1),a(0),p['u']/be,(-4*be-2*p['u']-p['tau'])/be,(3*be+p['u']+p['eps'])/be,-2*p['eps']/be,p['eps']/be]
    N=40;inv=[a(1)]
    for n in range(1,N+1):inv.append(-sum((q[j]*inv[n-j] for j in range(1,min(6,n)+1)),a(0)))
    inv2=[sum((inv[j]*inv[n-j] for j in range(n+1)),a(0)) for n in range(N+1)]
    bc={n+4:-inv[n]/(2*be) for n in range(N-3)}
    for j in range(1,7):
        for n in range(N+1):
            power=n+j+4
            if power<=N:bc[power]=bc.get(power,a(0))-j*q[j]*inv2[n]/(2*be*(n+5))
    bc[4]+=1/(2*be)
    check('profile_first_coefficients',bc[6].overlaps(3*p['u']/(10*be**2)) and bc[7].overlaps(q[3]/(5*be)))
    Rc=a('.2');core=a('.04')
    dq=sum((abs(q[j])*Rc**j for j in range(1,7)),a(0));dq1=sum((j*abs(q[j])*Rc**j for j in range(1,7)),a(0))
    guard=(dq/(Rc**2*(1-dq))+dq1/(5*Rc**2*(1-dq)**2))/(2*be)
    check('whole_complex_disk_profile_guard',dq<1 and guard<5 and guard<a('4.728011'))
    KL=2*be*(p['Fl']+23*p['g']/(180*a.pi()))/3
    def raw(n):
        k=a(n)-sig
        return core**k*(1+k*k)**a('1.5')/(2*k).sqrt()
    lead=abs(KL*bc[7]*84)*raw(7)
    others=sum((abs(KL*c*(n-4)*n*(n-3))*raw(n) for n,c in bc.items() if n>=6 and n!=7),a(0))
    n=N+1;ratio=core/Rc*(a(n+1)/n)**6
    tail=abs(KL)*5*Rc**6*core**(-sig)*a(n)**6*(core/Rc)**n/(1-ratio)
    check('entire_series_tail_ratio',ratio<1)
    wexp=6-sig;caps=[]
    for j in range(10):
        caps.append(5*a('1e9')*(wexp/2).exp()*a('.1')**wexp/(2*wexp).sqrt()*sum((comb(j,l)*(21/a(2).log())**l*factorial(j-l)*2**(j-l) for l in range(j+1)),a(0)))
    H9=sum((comb(9,j)*caps[j]**2 for j in range(10)),a(0)).sqrt();quantum=6000*K*H9;cost=lead-others-tail-quantum
    check('classical_leading_H3_lower',lead>a('.03747819294'))
    check('other_computed_terms_upper',others<a('.00222158030'))
    check('entire_analytic_tail_upper',tail<a('1.918e-18'))
    check('entire_H9_upper',H9<a('3.165345e21'))
    check('complete_normal_quantum_tail_upper',quantum<a('6.560e-11'))
    check('original_graph_conversion_lower',cost>a('.03525661258'))
    ceiling=2*a('1e-11')/cost
    check('necessary_phase_two_graph_targets',ceiling<a('5.672695e-10'))
    def primitive(y):return a('.9')**(2*s)*y**(3-2*s)/(3-2*s)-2*a('.9')**s*y**(3-s)/(3-s)+y**3/3
    total=a(0);cells=1024
    for j in range(cells):
        lo=a('.9')*j/cells;hi=a('.9')*(j+1)/cells
        piece=primitive(hi)-primitive(lo)
        # Positivity of the exact integral is analytic; Arb can enclose tiny zero-endpoint cells.
        total+=piece/(1-hi)**4
    phase=(16*total).sqrt();weakN=a('1.1e-10')/(a('.99')*p['Fl']).sqrt();allowance=phase*weakN
    check('integrable_entire_endpoint_majorant',3-2*s>0)
    check('entire_phase_functional_upper',phase<a('11.613'))
    check('original_unweighted_weak_norm_upper',weakN<a('1.105542e-10'))
    check('weak_actual_phase_allowance',allowance<a('1.283779e-9'))
    check('positive_estimator_allowance_exceeds_target',allowance/ceiling>a('2.26307') and allowance/ceiling<a('2.26309'))
    check('conditional_fast_screen_unchanged',a('.028913')+a('1e-11')*a('1.0706e10')<a('.1599164'))
    return dict(counterexample_phase_gap_lower=enc(gap),B4=enc(B4),counterexample_graph_leading_coefficient_lower=enc(leadg),phase_functional_upper=enc(phase),weak_seed_N_upper=enc(weakN),phase_allowance_from_weak_seed_bound=enc(allowance),profile_complex_guard=enc(guard),classical_leading_raw_H3=enc(lead),other_series_upper=enc(others),analytic_tail_upper=enc(tail),global_profile_H9_upper=enc(H9),whole_normal_quantum_upper=enc(quantum),original_graph_cost_lower_per_phase=enc(cost),necessary_phase_if_two_graph_radii=enc(ceiling),allowance_to_two_graph_ceiling=enc(allowance/ceiling),unbordered_fixed_phase_estimator_fails=True,explicit_boundary_alternative_quantified=True,actual_phase_moment_computed=False,full_coupled_resolvent_certified=False,post_elimination_mass_column_certified=False,actual_mass_computed=False,original_physical_strong_comparisons_closed=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,construction_obstruction=False,actual_controlled_solution=False,full_objective_complete=False,closed_audit_preserved=True,next_calculation_executed=False)

if __name__ == '__main__':
    mode = 'acceptance' if '--acceptance' in sys.argv else 'development'
    output = ROOT/'evidence'/(STEM+('' if mode == 'acceptance' else '-development')+'.json')
    hashes = {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
    import flint
    flint_root = Path(flint.__file__).resolve().parent
    runtime = {'python': sys.version, 'python_executable_sha256': sha256(Path(sys.executable).read_bytes()).hexdigest(),
               'flint_version': flint.__version__, 'sympy_version': sp.__version__,
               'flint_file_sha256': {str(p.relative_to(flint_root)): sha256(p.read_bytes()).hexdigest()
                                    for p in sorted(flint_root.rglob('*'))
                                    if p.is_file() and p.suffix in ['.py', '.pyd', '.dll']}}
    run = {'started_utc': datetime.now(timezone.utc).isoformat(), 'mode': mode,
           'input_sha256_before': hashes, 'precision_bits':512, 'runtime':runtime,
           'frozen_sources_utf8': {p: (ROOT/p).read_text(encoding='utf-8') for p in OWN}}
    try:
        run['result'] = execute()
        run['input_sha256_after'] = {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}
        check('all_frozen_inputs_unchanged', hashes == run['input_sha256_after'])
        run['status'] = 'passed'
    except Exception:
        run['status'] = 'failed'
        run['traceback'] = traceback.format_exc()
    run['checks'] = CHECKS
    run['finished_utc'] = datetime.now(timezone.utc).isoformat()
    ledger = json.loads(output.read_bytes()) if output.exists() else {'runs': []}
    ledger['runs'].append(run)
    output.write_bytes((json.dumps(ledger, indent=2)+'\n').encode())
    log = ROOT/'evidence'/(STEM+'.log')
    text = json.dumps({'mode': mode, 'status': run['status'], 'checks': len(CHECKS),
                       'result': run.get('result'), 'traceback': run.get('traceback')}, indent=2)+'\n'
    with log.open('ab') as f:
        f.write(text.encode())
    print(text)
    if run['status'] != 'passed':
        raise SystemExit(1)
