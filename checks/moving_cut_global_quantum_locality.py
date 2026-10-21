"""Global moving-cut locality diagnostics; local curved UV comparison remains open."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
STEM = 'moving-cut-global-quantum-locality'
PREVIOUS = 'energy-domain-and-chart-contact'
OWN = ['checks/moving_cut_global_quantum_locality.py', 'notes/'+STEM+'.md',
       'evidence/moving-cut-global-source-review.json']
ROWS = []; FIXTURES = []
PROTOCOL = {'precision_digits': [50, 80], 'masses': ['.75', '1.5', '3'],
            'matrix_dimension': 10, 'cut_powers': [4, 8, 16],
            'contour_nodes': [32, 64], 'mass_circle_fraction': '.25',
            'angular_cases': [['1', '1/9'], ['8', '1'], ['16000000000000', '320000000000/7']],
            'angular_partial_modes': 256, 'polynomial_degrees': [0, 4, 12],
            'exponential_rates': ['.5', '2'], 'scaled_identity_tolerance': '1e-40',
            'precision_tolerance': '1e-40', 'ordinary_continuum_proof': True,
            'matrices_are_actual_kernel_discretizations': False,
            'finite_matrix_double_commutator_retained': True}

def check(name, ok, detail=''):
    assert name not in {r['name'] for r in ROWS}
    row = {'name': name, 'passed': bool(ok), 'detail': str(detail)}
    ROWS.append(row); print(json.dumps(row), flush=True)
    if not ok: raise AssertionError(name)

def err(a, b): return abs(a-b)/(1+abs(b))
def out(x): return mp.nstr(x, mp.mp.dps)
def tr(A): return sum(A[i, i] for i in range(A.rows))
def fnorm(A): return mp.sqrt(sum(abs(x)**2 for x in A))
def prior(): return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths(): return sorted(set(prior()['input_sha256_after']) | set(OWN) | {'evidence/'+PREVIOUS+'.json'})
def zero(name, value):
    if isinstance(value, S.MatrixBase): value = sum(S.simplify(abs(x)**2) for x in value)
    value = S.simplify(value); check(name, value == 0, value)

def provenance():
    p = prior()
    check('predecessor_passed', p['status'] == 'passed' and len(p['checks']) == 164 and all(c['passed'] for c in p['checks']))
    check('predecessor_immutable', p['input_sha256_before'] == p['input_sha256_after'])
    bad = [f for f, h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest() != h]
    check('predecessor_hashes', not bad, bad)
    review = json.loads((ROOT/OWN[2]).read_bytes())
    check('review_inputs_pinned', all(f in p['input_sha256_after'] for f in review['retained_inputs']))
    check('scope', not review['actual_local_curved_cut_comparison_proved'] and not review['physical_fast_branch_excluded'])

def symbolic():
    x, d = S.symbols('x d', real=True)
    f = x+d*x*(1-x); fp = S.diff(f, x); u = S.Function('u')(x)
    transformed = S.sqrt(fp)*S.diff(u/S.sqrt(fp), x)/fp
    zero('unitary_pullback_half_density', transformed-u.diff(x)/fp-S.diff(1/fp, x)*u/2)
    cutoff = S.Function('cutoff')(x)
    commutator_u = S.diff(cutoff*u, x)-cutoff*S.diff(u, x)
    zero('differential_first_order_commutator', commutator_u-S.diff(cutoff, x)*u)
    zero('differential_double_commutator', S.diff(cutoff*cutoff*u,x)-2*cutoff*S.diff(cutoff*u,x)+cutoff**2*S.diff(u,x))
    s1 = S.Matrix([[0, 1], [1, 0]]); s2 = S.Matrix([[0, -S.I], [S.I, 0]]); s3 = S.diag(1, -1)
    chi = S.kronecker_product(s3, s1)
    normal = -S.kronecker_product(s1, S.eye(2))
    tangent = -S.kronecker_product(s2, s2)
    zero('paired_normal_boundary_green', normal*chi+chi*normal)
    zero('paired_tangential_boundary', tangent*chi-chi*tangent)
    zero('paired_elliptic_boundary_symbol', normal*tangent*chi+chi*normal*tangent)
    zero('paired_symbol_clifford', normal*tangent+tangent*normal)
    n, q, C = S.symbols('n q C', positive=True)
    count = 2*n*q+2*n*(n+1)
    zero('exact_counting_reserve', 4*n*(n+q)-count-2*n*(n+q-1))
    b, E = S.symbols('b E', positive=True)
    zero('exponential_count_integral', S.integrate(4*C*b*E**2*S.exp(-b*E), (E, 0, S.oo))-8*C/b**2)
    c = S.symbols('c', positive=True)
    zero('tail_exponent_split', c/2+c/4+c/4-c)
    zero('tail_count_factor', 8/(c/4)**2-128/c**2)

def operators(d):
    N = PROTOCOL['matrix_dimension']; A = mp.matrix(N); B = mp.matrix(N); Ad = mp.matrix(N)
    for i in range(N-1):
        base = mp.mpf('1.2')+mp.mpf(i+1)/10
        slope = mp.mpf((-1)**i)*(i+1)/50
        a = base+d*slope; b = a
        if i >= 7: b += mp.mpf('.4')*(i-6)
        A[i, i+1] = a; A[i+1, i] = -a
        B[i, i+1] = b; B[i+1, i] = -b
        Ad[i, i+1] = slope; Ad[i+1, i] = -slope
    chi = mp.diag([1, 1, 1, 1, mp.mpf('.8'), mp.mpf('.3'), 0, 0, 0, 0])
    xi = mp.diag([1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
    tt = []; td = []
    for i in range(N):
        y = mp.mpf(i+1)/8+d
        tt.append(y*y*mp.exp(-y) if i < 3 else 0)
        td.append((2*y-y*y)*mp.exp(-y) if i < 3 else 0)
    return A, B, Ad, chi, xi, mp.diag(tt), mp.diag(td)

def response(M, d):
    A, B, Ad, chi, xi, T, Td = operators(d); I = mp.eye(A.rows)
    R = (A+M*I)**-1; G = (B+M*I)**-1
    F = R-G; Q = A*chi-chi*A; E = xi*F*xi
    raw = M*M*(tr(T*R*T*R)-tr(T*G*T*G))-M*tr(T*T*F)
    Rp = -R*Ad*R; Gp = -G*Ad*G
    derivative = 2*M*M*(tr(Td*R*T*R)+tr(T*Rp*T*R)-tr(Td*G*T*G)-tr(T*Gp*T*G))
    derivative -= M*(2*tr(T*Td*F)+tr(T*T*(Rp-Gp)))
    return locals()

def matrix_fixtures():
    for j, mass in enumerate(PROTOCOL['masses']):
        levels = []
        for digits in PROTOCOL['precision_digits']:
            with mp.workdps(digits):
                M = mp.mpf(mass); d = mp.mpf(1)/16; v = response(M, d)
                A, B, Ad, chi, xi, T, R, G, E, Q = [v[k] for k in ['A','B','Ad','chi','xi','T','R','G','E','Q']]
                suffix = f'{j}_{digits}'; tol = mp.mpf(PROTOCOL['scaled_identity_tolerance'])
                check('common_patch_'+suffix, fnorm((A-B)*chi)+fnorm(chi*(A-B)) == 0)
                check('skew_adjoint_'+suffix, fnorm(A+A.T)+fnorm(B+B.T) == 0)
                grading = mp.diag([(-1)**i for i in range(A.rows)])
                check('grading_'+suffix, fnorm(grading*A+A*grading) == 0)
                # A finite matrix has [Q,chi]!=0; retain its exact defect.
                defect = xi*G*(Q*chi-chi*Q)*G*xi/2
                glued = -xi*G*Q*R*Q*G*xi+defect
                ge = fnorm(E-glued)/(1+fnorm(E))
                check('two_sided_gluing_'+suffix, ge < tol, out(ge))
                check('finite_commutator_defect_detected_'+suffix, fnorm(defect) > mp.mpf('1e-30'), out(fnorm(defect)))
                B0 = xi*G*xi
                expanded = M*M*tr(T*B0*T*E+T*E*T*B0+T*E*T*E)-M*tr(T*T*E)
                be = err(expanded, v['raw'])
                check('complete_bubble_expansion_'+suffix, be < tol, out(be))
                kernel = mp.matrix(A.rows)
                for i in range(A.rows):
                    for k in range(A.rows): kernel[i,k] = M*M*(R[i,k]*R[k,i]-G[i,k]*G[k,i])
                theta = [sum(kernel[i,k] for k in range(A.rows))-M*(R[i,i]-G[i,i]) for i in range(A.rows)]
                rowform = -sum(kernel[i,k]*(T[i,i]-T[k,k])**2 for i in range(A.rows) for k in range(A.rows))/2
                rowterm = sum(theta[i]*T[i,i]**2 for i in range(A.rows))
                re = err(rowform+rowterm, v['raw'])
                check('row_sum_and_contact_'+suffix, re < tol, out(re))
                check('omitted_diagonal_detected_'+suffix, abs(rowterm) > mp.mpf('1e-30'), out(abs(rowterm)))
                Fprime = -(R*R-G*G)
                ce = max(err(-M*M*Fprime[i,i]-M*(R[i,i]-G[i,i]), theta[i]) for i in range(A.rows))
                check('constant_probe_mass_derivative_'+suffix, ce < tol, out(ce))
                # Independent spectral diagonal and contour derivative of the actual-minus-cap resolvent.
                ea, Ua = mp.eighe(mp.j*A); eb, Ub = mp.eighe(mp.j*B); at = 1
                wa = [abs(Ua[at,k])**2 for k in range(A.rows)]; wb = [abs(Ub[at,k])**2 for k in range(A.rows)]
                def spectral(z):
                    return sum(wa[k]/(z-mp.j*ea[k])-wb[k]/(z-mp.j*eb[k]) for k in range(A.rows))
                check('spectral_diagonal_'+suffix, err(spectral(M), R[at,at]-G[at,at]) < tol)
                radius = M/4; contours = []
                for nodes in PROTOCOL['contour_nodes']:
                    fd = sum(spectral(M+radius*mp.exp(2*mp.pi*mp.j*k/nodes))*mp.exp(-2*mp.pi*mp.j*k/nodes) for k in range(nodes))/(nodes*radius)
                    er = err(fd, Fprime[at,at]); threshold = mp.mpf('1e-16' if nodes == 32 else '1e-34')
                    check('Cauchy_mass_derivative_'+suffix+'_'+str(nodes), er < threshold, out(er))
                    contours.append({'nodes': nodes, 'scaled_error': out(er)})
                envelope = sum(wa[k]/(abs(M-mp.j*ea[k])-radius)+wb[k]/(abs(M-mp.j*eb[k])-radius) for k in range(A.rows))
                ratio = abs(theta[at])/(5*M*envelope)
                check('holomorphic_probe_bound_'+suffix, ratio <= 1, out(ratio))
                v0 = response(M, mp.mpf(0)); cutvalues = []
                for power in PROTOCOL['cut_powers']:
                    cut = mp.mpf(2)**(-power); vd = response(M, cut)
                    diff = fnorm(vd['R']-v0['R']); bound = fnorm(vd['A']-v0['A'])/M**2
                    check('moving_resolvent_bound_'+suffix+'_'+str(power), diff <= bound, out(diff/bound))
                    remainder = abs(vd['raw']-v0['raw']-cut*v0['derivative'])
                    cutvalues.append({'power': power, 'response': out(vd['raw']), 'linear_remainder': out(remainder), 'resolvent_bound_ratio': out(diff/bound)})
                rem0 = mp.mpf(cutvalues[0]['linear_remainder']); rem1 = mp.mpf(cutvalues[-1]['linear_remainder'])
                check('moving_quadratic_remainder_'+suffix, rem1 < rem0/mp.mpf(1000000), out(rem1/rem0))
                levels.append({'precision': digits, 'response': out(v['raw']), 'gluing_error': out(ge),
                               'bubble_error': out(be), 'row_error': out(re), 'probe_error': out(ce),
                               'missing_diagonal': out(rowterm), 'finite_commutator_defect': out(fnorm(defect)), 'Cauchy_bound_ratio': out(ratio),
                               'contours': contours, 'cuts': cutvalues})
        with mp.workdps(80):
            pe = err(mp.mpf(levels[0]['response']), mp.mpf(levels[1]['response']))
            check('matrix_precision_'+str(j), pe < mp.mpf(PROTOCOL['precision_tolerance']), out(pe))
        FIXTURES.append({'fixture': 'moving_global_resolvents_and_diagonal', 'mass': mass, 'levels': levels, 'actual_kernel': False})

def maxpoly(p, a): return mp.mpf(1) if p <= a else (p/a)**p*mp.exp(-p+a)

def angular_fixtures():
    for j, (qs, cs) in enumerate(PROTOCOL['angular_cases']):
        levels = []
        for digits in PROTOCOL['precision_digits']:
            with mp.workdps(digits):
                q = mp.mpf(qs); C = mp.mpf(str(S.Rational(cs).p))/int(S.Rational(cs).q)
                values = []
                for n in [1, 2, 16, 256]:
                    count = 2*n*q+2*n*(n+1); bound = 4*n*(n+q)
                    check(f'count_{j}_{digits}_{n}', count <= bound, out(count/bound))
                for p in PROTOCOL['polynomial_degrees']:
                    for ctext in PROTOCOL['exponential_rates']:
                        c = mp.mpf(ctext); ap = maxpoly(p, c/2)
                        total = mp.mpf(0); tail = mp.mpf(0); L = mp.sqrt(16*(16+q)/C)
                        for n in range(1, PROTOCOL['angular_partial_modes']+1):
                            M = mp.sqrt(n*(n+q)/C); dn = 2*(q+2*n)
                            term = dn*(1+M)**p*mp.exp(-c*M); total += term
                            if M > L: tail += term
                        bound = 32*C*ap/c**2; tailbound = 128*C*ap/c**2*mp.exp(-c*L/4)
                        suffix = f'{j}_{digits}_{p}_{ctext}'
                        check('angular_envelope_'+suffix, total <= bound, out(total/bound))
                        check('angular_tail_'+suffix, tail <= tailbound, out(tail/tailbound))
                        if p > c/2:
                            xmax = p/(c/2)-1
                            check('polynomial_max_'+suffix, err((1+xmax)**p*mp.exp(-c*xmax/2), ap) < mp.mpf('1e-40'))
                        values.append({'p': p, 'c': ctext, 'partial_sum': out(total), 'bound_ratio': out(total/bound), 'tail_ratio': out(tail/tailbound)})
                levels.append({'precision': digits, 'values': values})
        with mp.workdps(80):
            pe = max(err(mp.mpf(a['partial_sum']), mp.mpf(b['partial_sum'])) for a,b in zip(levels[0]['values'], levels[1]['values']))
            check('angular_precision_'+str(j), pe < mp.mpf(PROTOCOL['precision_tolerance']), out(pe))
        FIXTURES.append({'fixture': 'full_tower_counting_and_tail_envelope', 'q': qs, 'C': cs, 'levels': levels})

def main():
    sys.stdout.reconfigure(encoding='utf-8', newline='\n'); sys.stderr.reconfigure(encoding='utf-8', newline='\n')
    ps = paths(); hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes(); error = None
    print(json.dumps({'event': 'start', 'story': 'S168 moving-cut global quantum locality', 'input_count': len(ps)}), flush=True)
    try: provenance(); symbolic(); matrix_fixtures(); angular_fixtures()
    except Exception: error = traceback.format_exc(); print(error, file=sys.stderr, flush=True)
    after = hashes(); good = error is None and before == after and all(c['passed'] for c in ROWS)
    run = {'story': 'S168 moving-cut global quantum locality', 'created_utc': datetime.now(timezone.utc).isoformat(),
           'status': 'passed' if good else 'failed', 'exception': error, 'python': sys.version,
           'sympy': S.__version__, 'mpmath': mp.__version__, 'protocol': PROTOCOL, 'checks': ROWS, 'fixtures': FIXTURES,
           'input_sha256_before': before, 'input_sha256_after': after, 'unchanged_inputs': before == after,
           'global_massive_cutoff_limit_derived': good, 'global_diagonal_terms_retained': True,
           'actual_local_curved_cut_comparison_proved': False, 'actual_nonlocal_finite_cut_comparison_proved': False,
           'pure_shape_finite_cut_limit_proved': False, 'physical_fast_branch_excluded': False,
           'actual_charged_kernel_evaluated': False, 'quantitative_uniform_temperature_bound_proved': False,
           'physical_action_averaging_adopted': False, 'bare_tadpole_independently_renormalized': False,
           'physical_exterior_state_constructed': False, 'continuous_causal_stability_proved': False,
           'finite_counterterm_changed': False, 'new_boundary_prescription_adopted': False,
           'order_reduction_adopted': False, 'formal_verification': False,
           'S167_complete': False, 'S168_complete': False, 'S077_complete': False, 'goal_complete': False}
    if not good:
        for kind, path in zip(['checker', 'note', 'review'], OWN): run['failed_'+kind+'_source_utf8'] = (ROOT/path).read_text(encoding='utf-8')
    dest = ROOT/'evidence'/(STEM+'.json'); data = json.loads(dest.read_bytes()) if dest.exists() else {'runs': []}
    data['runs'].append(run); dest.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8', newline='\n')
    print(json.dumps({'status': run['status'], 'checks': len(ROWS), 'unchanged_inputs': before == after}), flush=True)
    return 0 if good else 1

if __name__ == '__main__': raise SystemExit(main())
