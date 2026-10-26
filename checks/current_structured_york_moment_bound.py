"""Complete weighted mixed/York return; no quantum inverse is applied."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path') + 1])
import sympy as sp
from flint import arb, ctx
ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-structured-york-moment-bound'
OWN = ['notes/' + STEM + '.md', 'checks/current_structured_york_moment_bound.py',
       'evidence/' + STEM + '-review.json']
INPUTS = OWN + ['notes/' + n + '.md' for n in [
    'current-schur-residual-phase-domain', 'current-complete-bolt-energy-bound',
    'current-coupled-input-energy-bound', 'current-original-admission-map-budget-test',
    'current-native-output-structured-source-test', 'current-fixed-chart-complete-return-attempt',
    'current-controlled-solution-audit-closeout']] + ['evidence/' + n + '.json' for n in [
    'current-schur-residual-phase-domain', 'current-complete-bolt-energy-bound',
    'current-global-mass-conversion', 'current-localized-bolt-adjoint-sensitivity',
    'current-complete-source-square-factorization', 'current-native-output-structured-source-test',
    'current-coupled-mass-adjoint-runtime']]
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
    if not isinstance(x, dict): return arb(x.strip('[]'))
    d = x['decimal_enclosure']
    return arb(d['mid'] + ' +/- ' + d['rad']) * arb(10)**int(d['exponent'])

def prior(stem):
    r = json.loads((ROOT / 'evidence' / (stem + '.json')).read_bytes())['runs'][-1]
    check('accepted_' + stem, r['status'] == 'passed')
    return r['result']

def evaluate(expr, env):
    if expr in env: return env[expr]
    if expr == sp.pi: return arb.pi()
    if expr.is_Rational: return arb(int(expr.p)) / int(expr.q)
    if expr.is_Add: return sum((evaluate(x, env) for x in expr.args), arb(0))
    if expr.is_Mul:
        ans = arb(1)
        for x in expr.args: ans *= evaluate(x, env)
        return ans
    if expr.is_Pow:
        b, e = expr.args
        return evaluate(b, env)**(int(e) if e.is_Integer else evaluate(e, env))
    raise ValueError('Unspecified coefficient: ' + str(expr))

def execute():
    ctx.prec = 512
    rev = json.loads((ROOT / OWN[2]).read_bytes())
    check('frozen_ordinary_sources', all(sha(ROOT/p) == h for p,h in rev['scientific_sha256'].items()))
    rt = json.loads((ROOT / INPUTS[-1]).read_bytes())
    check('pinned_python', sha(Path(sys.executable)) == rt['python_sha256'])
    check('pinned_packages', all(sha(Path(pkg['root'])/p) == h
          for pkg in rt['packages'].values() for p,h in pkg['files'].items()))
    mixed = prior('current-schur-residual-phase-domain')
    complete = prior('current-complete-bolt-energy-bound')
    pars = prior('current-global-mass-conversion')['parameters']
    bp = prior('current-localized-bolt-adjoint-sensitivity')['parameters']
    square = prior('current-complete-source-square-factorization')
    native = prior('current-native-output-structured-source-test')
    beta, u, eps, tau, g = [dec(pars[k]) for k in ['beta','u','eps','tau','g']]
    C, B2 = [dec(bp[k]) for k in ['C','B2']]
    check('same_physical_parameters', (1/C).overlaps(u) and (C*beta/B2).overlaps(g))
    d = dec(square['parameters']['d'])
    Kb = dec(complete['constants']['Kb'])
    z = sp.Symbol('z')
    parse_symbols = {name:sp.Symbol(name) for name in ['beta','C','B2','eps','tau','cR','qflux','gR']}
    coeff = {k:sp.sympify(v, locals={**parse_symbols,'z':z})
             for k,v in mixed['complete_mixed_coefficients'].items()}
    all_symbols = set().union(*(v.free_symbols for v in coeff.values()))
    values = {'beta':beta, 'C':C, 'B2':B2, 'eps':eps, 'tau':tau,
              'cR':arb('1e22'), 'qflux':arb(2000)**4, 'gR':arb(2000)**-4}
    check('all_coefficient_parameters_specified', all(str(s) in values or s == z for s in all_symbols))
    env = {s:values[str(s)] for s in all_symbols if s != z}
    def poly(expr): return {n[0]:evaluate(c,env) for n,c in sp.Poly(sp.cancel(expr),z).terms()}
    polys = {k:poly(v) for k,v in coeff.items() if v != 0}
    check('complete_polynomial_outer_gain', all(min(p) >= 3 for p in polys.values())
          and coeff['xi3'] == 0 and coeff['xi4'] == 0)
    D = lambda x:z*sp.diff(x,z)
    c0,c1,c2 = [coeff['xi'+str(j)] for j in range(3)]
    ps = [c0+3*c1-D(c1)+9*c2-6*D(c2)+D(D(c2)), -c1-6*c2+2*D(c2), c2]
    s = sp.Function('s')(z)
    adj = c0*s+(3*c1*s-D(c1*s))+9*c2*s-6*D(c2*s)+D(D(c2*s))
    check('whole_differential_adjoint', sp.expand(adj-ps[0]*s-ps[1]*D(s)-ps[2]*D(D(s))) == 0)
    a = arb('.1'); b = 2*a
    a0 = 4*beta+2*u+tau; b4 = 3*beta+u+eps
    lowerQ = lambda x:beta-a0*x**3-2*abs(eps)*x**5-abs(eps)*x**6
    q0,q02 = lowerQ(a),lowerQ(b)
    check('actual_seed_outer_denominators', beta > arb('.16') and b4 > 0 and q0 > arb('.16') and q02 > arb('.16'))
    Q1 = 2*u+3*a0*a+4*abs(b4)*a*a+10*abs(eps)*a**3+6*abs(eps)*a**4
    major = lambda p, shift:sum((abs(v)*a**(n-shift) for n,v in p.items()),arb(0))
    H = [major(poly(p),1) for p in ps]
    Cb = major(polys['b_radial'],3)
    Cmu = sum((abs(v)/arb(2*n-3).sqrt() for n,v in polys['J_Y'].items()),arb(0))
    S0 = 1/arb('23.99').sqrt(); S1 = (beta/arb('5.99')).sqrt(); S2 = (beta/(6*d)).sqrt()
    Sdd = (S2+(3*arb('.17')+a*a*Q1)*S1/q0.sqrt())/q0
    Kdiff = H[0]*S0+H[1]*S1/q0.sqrt()+H[2]*Sdd
    KRb = Q1*Cb*S0*a**4/(q0*arb(88).sqrt())
    KRJ = Cmu*S0*(a**(arb(3)/2)/(q0*arb(3).sqrt())
                  +Q1*a**(arb(7)/2)/(5*q0*q0*arb(7).sqrt()))
    Kminus = Kb/a+Kdiff+KRb+KRJ
    # Whole-interval Green/energy identities, with arbitrary functions.
    Q = sp.Function('Q')(z); y = sp.Function('y')(z); T = sp.Function('T')(z)
    V0 = sp.Rational(2,3)*(sp.diff(Q,z,2)-4*sp.diff(Q,z)/z)/Q
    Ly = lambda h:-sp.diff(h,z,2)+(6/z**2+V0)*h
    Ay = lambda h:sp.Rational(2,3)*(sp.diff(h,z)/z**2-4*h/z**3)
    Aystar = lambda h:-sp.Rational(2,3)*(sp.diff(h,z)/z**2+2*h/z**3)
    check('weighted_York_forcing_cancellation', sp.simplify(
        Ay(T)*y/z**2+sp.Rational(2,3)*T*sp.diff(y,z)/z**4
        -sp.diff(sp.Rational(2,3)*T*y/z**4,z)) == 0)
    check('weighted_York_energy_identity', sp.simplify(
        -sp.diff(y,z,2)*y/z**2-(sp.diff(y,z)**2/z**2-3*y*y/z**4)
        -sp.diff(-y*sp.diff(y,z)/z**2-y*y/z**3,z)) == 0)
    psi = sp.Function('psi')(z)
    check('Green_endpoint_identity', sp.simplify(psi*Ly(y)-y*Ly(psi)
          -sp.diff(y*sp.diff(psi,z)-psi*sp.diff(y,z),z)) == 0)
    check('York_adjoint_boundary_term', sp.simplify(psi*Ay(T)-Aystar(psi)*T
          -sp.diff(sp.Rational(2,3)*psi*T/z**2,z)) == 0)
    v = sp.Function('v')(z)
    check('correct_singular_dual_weight', sp.simplify(Aystar(v/z**2)
          +sp.Rational(2,3)*sp.diff(v,z)/z**4) == 0)
    barpsi = z**-2*Q**sp.Rational(2,3)
    check('whole_York_positive_supersolution', sp.simplify(Ly(barpsi)
          -sp.Rational(2,9)*(sp.diff(Q,z)/Q)**2*barpsi) == 0)
    check('v_energy_identity', sp.simplify(
        (-sp.diff(v,z,2)+4*sp.diff(v,z)/z)*z*z*v
        -(z*z*sp.diff(v,z)**2-3*v*v)
        -sp.diff(-z*z*v*sp.diff(v,z)+3*z*v*v,z)) == 0)
    # Exact seed polynomial controls the potential only on [0,2a].
    v0 = 4*u/q02
    v1 = arb(2)/3*(6*a0+4*abs(b4)*b+6*abs(eps)*b**3)/q02
    M = (arb('.17')/beta)**(arb(2)/3); K = 3+arb('1e-9')
    Dv = M*(K.sqrt()*(1+1/arb(7).sqrt())/a**2
            +v0*a.sqrt()/3+v1*a**(arb(3)/2)/(2*arb(3).sqrt()))
    KVm = 1+1/arb(3).sqrt()
    quartic = 2*Cmu*S0/(3*beta)+Dv*Kminus/3
    q = arb('.338'); f = q/(1+q); transfer = q.sqrt()*f**2
    CS2 = KVm*Kminus*transfer; C4S2 = quartic*transfer
    Cv = 1+11*g/(180*arb.pi()); rD = (96*arb.pi()*arb('1e22')-1/(60*arb.pi()))*g
    cstar4 = 20*beta*(Cv-4*rD)
    check('same_nonresonant_quartic_denominator', cstar4 > 3 and cstar4 < 4 and
          cstar4.overlaps(dec(native['complete_structured_source']['mixed_quartic_normal_value'])))
    x = -1-2/q; t0,t1 = arb(1),x
    for j in range(2,30): t0,t1 = t1,2*x*t1-t0
    D29 = abs(t1); h2=(1+q)**2/D29; h3=(1+q)**3/D29; paff=arb('1.1e-10')*h3
    Cnu4 = C4S2/cstar4
    check('all_complete_low_source_bounds_finite', all(v.is_finite() and v > 0
          for v in [Cmu,Kminus,Dv,quartic,CS2,C4S2,Cnu4]))
    check('same_P32_affine_units', paff.overlaps(dec(square['complete_budget']['p_aff'])))
    constants = {k:enc(v) for k,v in {'q0':q0,'q02':q02,'Q1':Q1,'Cmu':Cmu,'Cb':Cb,
        'H0':H[0],'H1':H[1],'H2':H[2],'S0':S0,'S1':S1,'S2':S2,'Sdd':Sdd,'Kb':Kb,
        'Kdiff':Kdiff,'KRb':KRb,'KRJ':KRJ,'Kminus':Kminus,'KVm':KVm,'v0':v0,'v1':v1,
        'M':M,'Dv':Dv,'quartic_per_scalar_energy':quartic,'F2_scalar_energy_transfer':transfer,
        'complete_bar_S2_weighted_bound':CS2,'complete_bar_S2_quartic_bound':C4S2,
        'cstar4':cstar4,'classical_branch_quartic_bound':Cnu4,'rho_a':8*arb.pi()*B2/C,
        'p_aff':paff,'affine_bar_S2_weighted_bound':paff*CS2,
        'affine_classical_branch_quartic_bound':paff*Cnu4}.items()}
    return {'precision_bits':512,'constants':constants,
        'complete_source_norm':'||rho_a (V*)^-1 Bop F^2 u / z||_L2(z^-4 dz), E(u)<=1',
        'unknown_required_maps':{k:None for k in ['structured_outer_H3','structured_regular_atlas_H3',
            'Cdef3','I_Y','Cg3','original_physical_output','cold_embedding','radius_atlas',
            'signed_coupled_entries','independent_phase_columns','actual_kappa_m','application_errors',
            'nonlinear_neighborhood']},
        'conditional_ledger':{k:enc(v) for k,v in {'h2':h2,'h3':h3,
            'zero_error_defect_product_ceiling':arb('.5')/h3,
            'whole_graph_margin_ceiling':arb('5.34235e-12')/paff,
            'two_fast_screen_only':arb('.028913')+arb('1e-11')*arb('1.0706e10')}.items()},
        'scope':{'complete_weighted_low_source_bound_proved':True,
            'both_global_moments_retained_and_bounded':True,
            'inadmissible_form_duality_not_used':True,'first_source_regular_derivatives_bounded':False,
            'complete_coupled_bound_certified':False,'actual_polynomial_application_executed':False,
            'original_defect_certified':False,'new_physical_gate_closed':False,
            'actual_controlled_solution':False,'construction_obstruction':False,'full_goal_complete':False}}

if __name__ == '__main__':
    mode = 'acceptance' if '--acceptance' in sys.argv else 'development'
    output = ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes = {p:sha(ROOT/p) for p in INPUTS}
    run = {'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,
        'input_sha256_before':hashes,'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
        'runtime':{'python':sys.version,'executable_sha256':sha(Path(sys.executable))}}
    try:
        run['result'] = execute()
        run['input_sha256_after'] = {p:sha(ROOT/p) for p in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after'])
        run['status'] = 'passed'
    except Exception:
        run['status'] = 'failed'; run['traceback'] = traceback.format_exc()
    run['checks'] = CHECKS; run['finished_utc'] = datetime.now(timezone.utc).isoformat()
    ledger = json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run); output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    message = json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),
                         'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as stream: stream.write(message.encode())
    print(message,flush=True)
    if run['status'] != 'passed': raise SystemExit(1)
