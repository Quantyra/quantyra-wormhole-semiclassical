"""Full-source finite-output strong-residual test; not an inverse implementation."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction
from math import comb
import argparse
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-strong-residual-rank-test'
RUNTIME = ROOT / 'evidence/current-coupled-mass-adjoint-runtime.json'
runtime = json.loads(RUNTIME.read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx

INPUTS = [
    'notes/'+STEM+'.md', 'checks/current_strong_residual_rank_test.py',
    'notes/current-structured-source-graph-test.md',
    'notes/current-finite-local-action-assembly.md',
    'notes/current-canonical-heat-restoration-assembly.md',
    'notes/current-seed-full-hessian-realization.md',
    'notes/current-quartic-ritz-rank-test.md',
    'notes/current-geodesic-two-sided-remainder.md',
    'evidence/current-structured-source-graph-test.json',
    'evidence/current-geodesic-two-sided-remainder.json',
    'evidence/current-canonical-heat-restoration-assembly.json',
    'evidence/current-spatial-inverse-output-test.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]
BASELINE = Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def dec(v):
    d = v['decimal_enclosure']
    return arb(int(d['mid']), int(d['rad']))*arb(10)**int(d['exponent'])


def enc(v):
    m, r, e = v.mid_rad_10exp(60)
    return {'ball': v.str(40), 'decimal_enclosure':
            {'mid': str(m), 'rad': str(r), 'exponent': str(e)}}


def receipt(stem, status='passed'):
    obj = json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())
    run = obj.get('runs', [obj])[-1]
    assert run['status'] == status
    return run


def pick(run, kind):
    return next(x for x in run['numerics'] if x['kind'] == kind)


def add(*ps):
    out = {}
    for p in ps:
        for k, v in p.items():
            out[k] = out.get(k, arb(0))+v
    return {k: v for k, v in out.items() if not v.is_zero()}


def scale(p, x):
    return {k: v*x for k, v in p.items()}


def mul(p, q):
    out = {}
    for i, a in p.items():
        for j, b in q.items():
            out[i+j] = out.get(i+j, arb(0))+a*b
    return out


def diff(p):
    return {k-1: k*v for k, v in p.items() if k}


def shift(p, n):
    return {k+n: v for k, v in p.items()}


def value(p, z):
    return sum((v*z**k for k, v in p.items()), arb(0))


def fadd(*cs):
    return [add(*(c[i] for c in cs if i < len(c))) for i in range(max(map(len, cs)))]


def fscale(c, x):
    return [scale(p, x) for p in c]


def fmul(c, p):
    return [mul(q, p) for q in c]


def fdiff(c):
    return fadd([diff(p) for p in c], [{}]+c)


class Jet:
    """Value, two directional derivatives, mixed derivative: outward Arb."""
    def __init__(self, v=0, i=0, j=0, ij=0):
        self.v, self.i, self.j, self.ij = map(arb, [v, i, j, ij])

    @staticmethod
    def lift(x):
        return x if isinstance(x, Jet) else Jet(x)

    def __add__(self, x):
        x = Jet.lift(x)
        return Jet(self.v+x.v, self.i+x.i, self.j+x.j, self.ij+x.ij)
    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.v, -self.i, -self.j, -self.ij)

    def __sub__(self, x):
        return self+-Jet.lift(x)

    def __rsub__(self, x):
        return Jet.lift(x)+-self

    def __mul__(self, x):
        x = Jet.lift(x)
        return Jet(self.v*x.v, self.i*x.v+self.v*x.i,
                   self.j*x.v+self.v*x.j,
                   self.ij*x.v+self.i*x.j+self.j*x.i+self.v*x.ij)
    __rmul__ = __mul__

    def __truediv__(self, x):
        return self*(1/arb(x))

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        out = Jet(1)
        for _ in range(n):
            out = out*self
        return out

    def exp(self):
        ev = self.v.exp()
        return Jet(ev, ev*self.i, ev*self.j, ev*(self.ij+self.i*self.j))


def local_hessian(i, j, z, Q, Q1, Q2, C, B2, q):
    v = [Jet(0, int(i == k), int(j == k)) for k in range(6)]
    n, n1, n2, t, t1, t2 = v
    phi, phi1, phi2 = -z.log()+n, -1/z+n1, 1/z**2+n2
    L = -Q2-2*(Q*t2+Q1*t1)
    M = Q*phi2+Q1*phi1
    P = Q*phi1**2
    D = L-6*M-6*P
    en2, en4, et2, etm2 = (2*n).exp(), (4*n).exp(), (2*t).exp(), (-2*t).exp()
    LE = en2*D/(2*z**2)+en2*et2/(C*z**2)+3*arb(2000)**-10*B2/C*en4*et2/z**4
    LEM = -arb.pi()*q**2*arb(2000)**8/(C*B2)*et2
    LR = 8*arb.pi()*arb('1e22')/(C*B2)*(4*et2+4*C*D+C*C*etm2*D**2)
    bracket = phi*L**2/60+L*M/30-M**2/10+L*P/30-7*M*P/90-7*P**2/180
    LW = (C*C*etm2*bracket-C*phi*L/18+(arb(1)/15-q*q/3)*phi*et2-C*P/18)/(2*arb.pi()*C*B2)
    return dict(Einstein=LE.ij, magnetic=LEM.ij, R2=LR.ij, Weyl=LW.ij)


def execute():
    ctx.prec = 512
    assert digest(Path(sys.executable)) == runtime['python_sha256']
    for pkg in runtime['packages'].values():
        for name, expected in pkg['files'].items():
            assert digest(Path(pkg['root'])/name) == expected
    geo = receipt('current-geodesic-two-sided-remainder')
    ep = pick(geo, 'actual_seed_correlated_endpoint')
    gp = pick(geo, 'actual_finite_seed_two_sided_parameters')
    rc = pick(receipt('current-canonical-heat-restoration-assembly'), 'common_restoration_coefficients')
    sc = receipt('current-structured-source-graph-test',
                 'source_graph_bound_checked_complete_application_unresolved')['constants']
    sp = receipt('current-spatial-inverse-output-test',
                 'error_estimator_test_passed_signed_action_unresolved')['result']['constants']
    Q = {k: dec(v) for k, v in zip([0, 2, 3, 4, 5, 6], ep['Q_coefficients'])}
    C, B2, lam = dec(gp['C']), dec(gp['B2']), dec(gp['lambda_value'])
    beta, q = Q[0], dec(gp['q'])
    epsilon = dec(sp['epsilon'])
    d, Cv, Krest = dec(sc['d']), dec(sc['C_v']), dec(sc['mixed_remainders_total'])
    assert beta.overlaps(lam*B2/C) and arb('.99') < dec(ep['z_h'])
    assert 0 < epsilon < arb('1.227e-74')
    ds, x = Q[6]*C, lam*B2
    assert abs(ds) < arb('8.6e-19')
    m = diff(shift(diff(Q), -4))
    J = [scale(diff(m), 2), m]
    QJ = fmul(J, Q)
    Ap = fdiff(QJ)
    xi = fadd(Ap, [scale(shift(mul(m, m), 4), -1)])
    nu = fadd(fscale(xi, arb(1)/2), [shift(p, -1) for p in QJ])
    T = fadd(fscale(xi, -1), fscale(Ap, -arb(1)/2))
    # The original York differential equation as a Laurent/jet identity.
    qqjp = fmul(fdiff(J), mul(Q, Q))
    rhs = fmul(fadd([shift(scale(p, 4), -1) for p in xi], fscale(fdiff(xi), -1)), Q)
    york = fadd(fscale(fdiff(qqjp), -1), [shift(scale(p, 4), -1) for p in qqjp], fscale(rhs, -1))
    assert all(v.contains(0) for p in york for v in p.values())
    moment = fadd(fmul(J, m), fscale(fdiff([mul(m, m)]), -1))
    assert all(v.contains(0) for p in moment for v in p.values())

    # Exact polynomial integration for each bump-derivative norm.
    p = {2*i: Fraction((-1)**i*comb(8, i)) for i in range(9)}
    f, exact = [], []
    for k in range(6):
        moment2 = sum((a*b*Fraction(2, i+j+1) for i, a in p.items()
                       for j, b in p.items() if (i+j) % 2 == 0), Fraction(0))
        assert moment2 > 0
        exact.append(str(moment2))
        f.append((arb(moment2.numerator)/moment2.denominator).sqrt())
        p = {i-1: i*a for i, a in p.items() if i}

    a, b, cells = arb('.4'), arb('.6'), 128
    # The rational endpoints .4+j/(640) are represented by interval arithmetic.
    intervals = []
    for j in range(cells):
        lo, hi = a+(b-a)*j/cells, a+(b-a)*(j+1)/cells
        intervals.append(lo.union(hi))
    def sup(p):
        return max(value(p, z).abs_upper() for z in intervals)
    def inf(p):
        return min(value(p, z).lower() for z in intervals)
    bounds = {}
    for name, c in [('xi', xi), ('nu', nu), ('T', T)]:
        bounds[name] = []
        for ci in c:
            row, cur = [], ci
            for j in range(4):
                row.append(sup(cur))
                cur = diff(cur)
            bounds[name].append(row)
    Mdim = 1000001
    h = (b-a)/Mdim
    w = 2/h
    a2low = inf(xi[2])
    D0 = a2low*w*w*f[2]-bounds['xi'][1][0]*w*f[1]-bounds['xi'][0][0]*f[0]
    assert a2low > 0 and D0 > 0
    dz = {}
    for name, cs in bounds.items():
        dz[name] = []
        for k in range(4):
            nk = sum((comb(k, j)*cs[i][j]*w**(i+k-j)*f[i+k-j]
                      for i in range(3) for j in range(k+1)), arb(0))
            dz[name].append(b*b*nk/D0)
    RH1, RH2 = dz['xi'][1]/a**2, dz['xi'][2]/a**2
    Qmax, Qmin, Q1max, Q2max = sup(Q), inf(Q), sup(diff(Q)), sup(diff(diff(Q)))
    drift = add(diff(Q), scale(shift(Q, -1), -2))
    L1 = Qmax.sqrt()*b/beta.sqrt()*RH1
    L2 = b*b/beta*(Qmax*RH2+sup(drift)*RH1)
    glow = Qmin*a**4*arb.pi()**2/(beta*b*b*h*h)
    tau = {4: (1+3*x+ds)/x, 5: -5*ds/x, 6: arb(9)/2*ds/x}
    taulow = inf(tau)
    Dt = scale(shift(add(mul(Q, diff(diff(tau))), mul(drift, diff(tau))), 2), 1/beta)
    t2 = sup(Dt)
    # Check the scalar coefficient guards used in the ordinary upper bound.
    Rbar = {0: arb(-12), 5: 4*ds/x, 6: -6*ds/x}
    ee = lam/(60*arb.pi())
    Kupper = 6+(2*d+11*ee/6)*sup(Rbar)+(22*ee/3)*7
    Vupper = sup(add(scale(Rbar, 2), {0: 48*arb(2000)**-10/lam}))
    assert Kupper < arb('6.01') and Vupper < arb('24.01')
    Aup = 6*d*L2**2+arb('6.01')*L1**2+arb('24.01')
    blow = 4*(Cv*taulow+d*taulow*glow-d*t2/2)-Krest*Aup.sqrt()
    assert Qmin > 0 and taulow > 0 and blow > 0
    Blow = blow**2/Aup

    # Complete affine local Hessians on an interval cover; no point sampling.
    labels = ['Einstein', 'magnetic', 'R2', 'Weyl']
    matrices = {name: [[arb(0) for _ in range(6)] for _ in range(6)] for name in labels}
    for z in intervals:
        q0, q1, q2 = value(Q,z), value(diff(Q),z), value(diff(diff(Q)),z)
        for i in range(6):
            for j in range(i,6):
                co = local_hessian(i,j,z,q0,q1,q2,C,B2,q)
                for name in labels:
                    bound = co[name].abs_upper()
                    matrices[name][i][j] = max(matrices[name][i][j],bound)
                    matrices[name][j][i] = matrices[name][i][j]
    jets = dz['nu'][:3]+dz['T'][:3]
    local = {name: sum((mat[i][j]*jets[i]*jets[j] for i in range(6) for j in range(6)),arb(0))/beta
             for name,mat in matrices.items()}
    lt0, lt1, lt2, lt3 = dz['T']
    ld = Qmax*lt2+Q1max*lt1
    ea2 = dec(rc['beta2']).abs_upper()*(4*Q2max**2*lt0**2+16*Q2max*lt0*ld+8*ld**2)/beta
    elll = dec(rc['LLL_Hessian_coefficient']).abs_upper()*Qmax*lt1**2/beta
    evac = dec(rc['vacuum_coefficient_upper'])*lt0**2/beta
    A, D = Qmax*lt3, arb('1.5')*Q1max*lt2+Q2max*lt1/2
    s2 = 4*C*C/q
    integraltest = 2*C*C*(q+2)/(q+1)**2+2*C*C/(q+1)
    assert integraltest < s2
    egeo = 2*s2/(525*arb.pi()*B2*beta)*Qmax*(arb(5)/8*A*A+arb('1.5')*A*D+D*D)
    Eup = (1+epsilon)/(1-arb('.0022'))*(sum(local.values(),arb(0))+ea2+elll+evac+egeo)
    stronglower, weaklower = Blow/Eup.sqrt(), Blow/Eup
    strongscreen, weakscreen = arb('4.964e-32'), arb('4.989e-32')
    multiplier = dec(sc['strong_residual_to_COMPLETE_residual_multiplier'])
    weak_application = dec(sc['energy_dual_residual_to_application_multiplier'])*weakscreen
    weak_debit = multiplier*dec(sc['c_H']).sqrt()*weakscreen
    constants = dict(M=arb(Mdim), h=h, D0=D0, a2_lower=a2low, Q_lower=Qmin,
                     tau_lower=taulow, scalar_K_upper=Kupper, scalar_V_upper=Vupper,
                     scalar_L1_upper=L1, scalar_L2_upper=L2, gradient_squared_lower=glow,
                     scalar_energy_upper=Aup, mixed_main_lower=blow+Krest*Aup.sqrt(),
                     mixed_remainder_debit=Krest*Aup.sqrt(), complete_mixed_lower=blow,
                     complete_B_lower=Blow, canonical_a2_upper=ea2,
                     canonical_LLL_upper=elll, shared_vacuum_upper=evac,
                     full_angular_sum_upper=s2, complete_nonlocal_energy_upper=egeo,
                     complete_Et_energy_upper=Eup,
                     strong_residual_lower_bound=stronglower,
                     energy_dual_residual_lower_bound=weaklower,
                     strong_ratio_to_sufficient_screen=stronglower/strongscreen,
                     weak_ratio_to_sufficient_screen=weaklower/weakscreen,
                     strong_estimator_complete_source_debit_floor=multiplier*stronglower,
                     strong_debit_floor_to_ENTIRE_residual_screen=multiplier*stronglower/arb('1e-30'),
                     conditional_complete_source_allocation=multiplier*strongscreen,
                     conditional_weak_application_error_NOT_attained=weak_application,
                     conditional_weak_source_debit_NOT_attained=weak_debit)
    constants.update({'local_'+k+'_energy_upper':v for k,v in local.items()})
    return dict(constants={k:enc(v) for k,v in constants.items()},
                exact_bump_derivative_squared_norms=exact,
                ordinary_jet_norm_comparisons={k:[enc(x) for x in v] for k,v in dz.items()},
                local_coefficient_absolute_matrices={k:[[enc(x) for x in row] for row in mat] for k,mat in matrices.items()},
                interval_cover_cells=cells, maximum_fixed_output_dimension=Mdim-1,
                original_York_identity_checked=True, compact_York_moment_checked=True,
                strong_finite_output_screen_rejected=bool(stronglower > strongscreen),
                weak_test_rejects_screen=bool(weaklower > weakscreen),
                actual_S_application_error=None, actual_K_application_error=None,
                actual_complete_numerical_residual=None,
                certifies_inverse_or_wormhole_nonexistence=False)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output', default='evidence/'+STEM+'.json')
    args=parser.parse_args()
    output=ROOT/args.output
    if output.exists():
        raise RuntimeError('Receipts are immutable; choose a fresh output path.')
    protected=json.loads(BASELINE.read_bytes())
    verify=lambda: all(digest(Path(p))==h for p,h in protected.items())
    assert verify()
    before={p:digest(ROOT/p) for p in INPUTS}
    data={'started_utc':datetime.now(timezone.utc).isoformat(),
          'previous_goal_turn':'No progress: requested email recheck supplied no new bound.',
          'input_sha256_before':before,
          'ordinary_argument_utf8':(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
          'arithmetic_source_utf8':Path(__file__).read_text(encoding='utf-8')}
    try:
        data['result']=execute()
        data['status']='finite_output_estimator_test_executed_full_application_unresolved'
    except Exception as exc:
        data['status']='failed'
        data['error']=repr(exc)
        raise
    finally:
        data['input_sha256_after']={p:digest(ROOT/p) for p in INPUTS}
        data['inputs_unchanged']=before==data['input_sha256_after']
        data['protected_paths_unchanged']=verify()
        data['protected_path_count']=len(protected)
        data['finished_utc']=datetime.now(timezone.utc).isoformat()
        output.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    assert data['inputs_unchanged'] and data['protected_paths_unchanged']
    for k,v in data['result']['constants'].items():
        print(k+': '+v['ball'])
    print('Receipt: '+str(output))


if __name__=='__main__':
    main()
