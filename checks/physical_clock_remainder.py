"""S136.3: exact spherical-to-clock identity; no physical bound certificate."""
from pathlib import Path
import ast
import datetime
import hashlib
import json
import platform
import subprocess
import sys
import traceback

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'evidence/physical-clock-remainder.json'
LOG = ROOT / 'evidence/physical-clock-remainder.log'
DPS = 60
TOL = {'exact_numerical': '1e-35'}
SOURCES = [
    'notes/physical-clock-remainder.md',
    'checks/physical_clock_remainder.py',
    'notes/spherical-anomaly-sources.md',
    'notes/exact-spherical-em-reduction.md',
    'notes/sourced-spherical-mass-charge-balance.md',
    'notes/jt-prescribed-source-constraint.md',
    'notes/rn-ads-finite-cut-matching.md',
    'notes/retimed-returned-quantum-stress.md',
]


def write(path, text):
    with path.open('w', encoding='utf-8', newline='\n') as handle:
        handle.write(text)


def record(text):
    with LOG.open('a', encoding='utf-8', newline='\n') as handle:
        handle.write(text + '\n')
    print(text, flush=True)


def sha(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def close(error):
    assert abs(error) < mp.mpf(TOL['exact_numerical']), str(error)


def symbolic_checks():
    b, G, xi, phi = sp.symbols('b G xi phi', positive=True)
    X = b*b + 2*G*phi
    rm2 = b*b*(1+3*xi)
    a2 = b*b/(1+6*xi)
    potential = 2+6*xi*X/(b*b)-2*rm2/X
    remainder = potential/(4*G)-phi/a2
    assert sp.factor(remainder+2*G*rm2*phi*phi/(b**4*X)) == 0
    # Independent elimination of normal and curvature from the physical
    # moving-cut expression: variables are arbitrary local jets.
    sign, ell, vel, acc, wt, ws, pt, ps, p, cot = sp.symbols(
        'sign ell vel acc wt ws pt ps p cot', nonzero=True)
    sign = sp.Integer(1)  # s^2=1; check the other orientation explicitly.
    results = []
    for sign in [-1, 1]:
        norm = ell*sp.sqrt(1-vel*vel)
        nphi = sign*(vel*pt+ps)/norm
        nomega = sign*(vel*wt+ws)/norm
        curvature = sign*(acc/(1-vel*vel)+vel*wt+ws-cot)/norm
        rhs = norm*(nphi-(curvature-nomega)*p)+sign*(acc*p/(1-vel*vel)-vel*pt)
        assert sp.simplify(rhs-sign*(ps+cot*p)) == 0
        results.append(sign)
    return {'potential_identity': 'exact zero', 'moving_cut_orientations': results}


class Fields:
    """Synthetic smooth metric; stated bounds hold globally in time."""
    def __init__(self):
        self.b = mp.mpf('1.4')
        self.L4 = mp.mpf('2.1')
        self.G = mp.mpf('.02')
        self.B = mp.mpf('.2')
        self.eta = mp.mpf(1)/2000
        self.N = self.B/self.eta
        self.xi = self.b**2/self.L4**2
        self.C = 1+6*self.xi
        self.a = self.b/mp.sqrt(self.C)
        self.rm2 = self.b**2*(1+3*self.xi)
        self.mu = 2*self.G*self.B/self.b**2
        self.vartheta = (1+3*self.xi)/self.C

    def psi(self, t, s):
        return (mp.mpf('.4')+mp.mpf('.05')*mp.cos(mp.mpf('1.2')*t))/mp.sin(s)+mp.mpf('.02')*mp.sin(s)**2*mp.cos(mp.mpf('.6')*t)

    def omega(self, t, s):
        return self.mu*mp.mpf('.08')*(mp.cos(mp.mpf('.8')*t)/mp.sin(s)+mp.mpf('.2')*mp.sin(s))

    def phi(self, t, s):
        return self.B*self.psi(t, s)

    def X(self, t, s):
        return self.b**2+2*self.G*self.phi(t, s)

    def lapse(self, t, s):
        return self.a*mp.exp(self.omega(t, s))/mp.sin(s)

    def metric(self, i, t, s, th, ph):
        if i == 0:
            return -self.lapse(t, s)**2
        if i == 1:
            return self.lapse(t, s)**2
        return self.X(t, s)*(1 if i == 2 else mp.sin(th)**2)

    def jets(self, fn, t, s):
        return [fn(t, s), mp.diff(fn, (t, s), (1, 0)),
                mp.diff(fn, (t, s), (0, 1)), mp.diff(fn, (t, s), (2, 0)),
                mp.diff(fn, (t, s), (1, 1)), mp.diff(fn, (t, s), (0, 2))]

    def orbit(self, t, s):
        p, pt, ps, ptt, pts, pss = self.jets(self.phi, t, s)
        w, wt, ws, _, _, _ = self.jets(self.omega, t, s)
        x = self.X(t, s)
        u = -2*self.G*self.rm2*p*p/(self.b**4*x)
        q = (mp.exp(2*w)-1)*p+self.a**2*mp.exp(2*w)*u
        common = -self.G*(pt*pt+ps*ps)/(2*x)-wt*pt-ws*ps
        geom = [[common-q/mp.sin(s)**2, -self.G*pt*ps/x-wt*ps-ws*pt],
                [-self.G*pt*ps/x-wt*ps-ws*pt, common+q/mp.sin(s)**2]]
        linear = [[-pss-mp.cot(s)*ps+p/mp.sin(s)**2, -pts-mp.cot(s)*pt],
                  [-pts-mp.cot(s)*pt, -ptt-mp.cot(s)*ps-p/mp.sin(s)**2]]
        physical = [[linear[i][j]-geom[i][j] for j in range(2)] for i in range(2)]
        return linear, geom, physical


def four_dimensional_einstein(fields, point):
    """General diagonal metric Ricci tensor from derivative jets, no reduction."""
    diagonal = [fields.metric(i, *point) for i in range(4)]
    first = [[mp.mpf(0) for k in range(4)] for i in range(4)]
    second = [[[mp.mpf(0) for l in range(4)] for k in range(4)] for i in range(4)]
    for i in range(4):
        fn = lambda *coords, i=i: fields.metric(i, *coords)
        for k in range(3):
            orders = [0]*4
            orders[k] = 1
            first[i][k] = mp.diff(fn, point, tuple(orders))
            for l in range(3):
                orders2 = orders.copy()
                orders2[l] += 1
                second[i][k][l] = mp.diff(fn, point, tuple(orders2))
    gamma = [[[mp.mpf(0) for n in range(4)] for m in range(4)] for r in range(4)]
    dgamma = [[[[mp.mpf(0) for n in range(4)] for m in range(4)] for r in range(4)] for k in range(4)]
    for r in range(4):
        for m in range(4):
            for n in range(4):
                num = (first[r][m] if r == n else 0)+(first[r][n] if r == m else 0)-(first[m][r] if m == n else 0)
                gamma[r][m][n] = num/(2*diagonal[r])
                for k in range(4):
                    dnum = (second[r][k][m] if r == n else 0)+(second[r][k][n] if r == m else 0)-(second[m][k][r] if m == n else 0)
                    dgamma[k][r][m][n] = dnum/(2*diagonal[r])-first[r][k]*gamma[r][m][n]/diagonal[r]
    ricci = [[mp.mpf(0) for n in range(4)] for m in range(4)]
    for m in range(4):
        for n in range(4):
            ricci[m][n] = sum(dgamma[r][r][m][n]-dgamma[n][r][m][r] for r in range(4))
            ricci[m][n] += sum(gamma[r][r][l]*gamma[l][m][n]-gamma[r][n][l]*gamma[l][m][r] for r in range(4) for l in range(4))
    scalar = sum(ricci[i][i]/diagonal[i] for i in range(4))
    einstein = [[ricci[m][n]-(diagonal[m]*scalar/2 if m == n else 0) for n in range(4)] for m in range(4)]
    return diagonal, gamma, einstein


def full_tensor_checks(fields):
    rows = []
    for t, s in [('0', '.3'), ('.7', '.8'), ('2', '1.5'), ('7', str(mp.pi-mp.mpf('.4')))]:
        t, s = mp.mpf(t), mp.mpf(s)
        th = mp.mpf('.8')
        diagonal, _, einstein = four_dimensional_einstein(fields, (t, s, th, mp.mpf('.4')))
        x = fields.X(t, s)
        # Choose gc^2=4*pi*G, so P^2=r_m^2. The 4D magnetic tensor is
        # independently reconstructed from F_theta_phi=P*sin(theta).
        gc2 = 4*mp.pi*fields.G
        f2 = 2*fields.rm2/x**2
        mag = [-diagonal[0]*f2/(4*gc2), -diagonal[1]*f2/(4*gc2),
               fields.rm2/(2*gc2*x), fields.rm2*mp.sin(th)**2/(2*gc2*x)]
        remainder = [[einstein[i][j]+((-3/fields.L4**2)*diagonal[i]-8*mp.pi*fields.G*mag[i] if i == j else 0) for j in range(4)] for i in range(4)]
        linear, geom, physical = fields.orbit(t, s)
        errors = [abs(x*remainder[i][j]/(2*fields.G)-physical[i][j]) for i in range(2) for j in range(2)]
        for error in errors:
            close(error)
        J = 2*remainder[2][2]/x
        psi, pt, ps, ptt, _, pss = fields.jets(fields.psi, t, s)
        w, _, _, wtt, _, wss = fields.jets(fields.omega, t, s)
        boxp = mp.sin(s)**2*(-ptt+pss)
        boxw = mp.sin(s)**2*(-wtt+wss)
        normp = mp.sin(s)**2*(-pt*pt+ps*ps)
        angular = -2-2*boxw-fields.mu*boxp/(1+fields.mu*psi)+fields.mu**2*normp/(2*(1+fields.mu*psi)**2)
        angular += mp.exp(2*w)*(6*fields.xi/fields.C+2*fields.vartheta/(1+fields.mu*psi)**2+fields.a**2*J)
        close(angular)
        close(remainder[3][3]-mp.sin(th)**2*remainder[2][2])
        rows.append(dict(t=str(t),sigma=str(s),orbit_component_errors=list(map(str,errors)),angular_error=str(abs(angular)),J_X=str(J)))
    return rows


def moving_cut_checks(fields):
    rows = []
    for sign in [-1, 1]:
        cut = (lambda t: mp.mpf('.3')+mp.mpf('.015')*mp.sin(mp.mpf('.7')*t)) if sign == -1 else (lambda t: mp.pi-mp.mpf('.3')-mp.mpf('.01')*mp.cos(mp.mpf('.5')*t))
        for rawt in ['0', '.7', '2', '7']:
            t = mp.mpf(rawt)
            s = cut(t)
            v, acc = mp.diff(cut, t), mp.diff(cut, t, 2)
            ell = fields.lapse(t, s)
            _, wt, ws, _, _, _ = fields.jets(fields.omega, t, s)
            # Build the 2D connection from metric derivatives, then contract
            # the covariant acceleration with the spacelike unit normal.
            diag = [-ell*ell, ell*ell]
            dlog = [wt, ws-mp.cot(s)]
            dg = [[2*diag[i]*dlog[k] for k in range(2)] for i in range(2)]
            gamma = [[[( (dg[r][m] if r == n else 0)+(dg[r][n] if r == m else 0)-(dg[m][r] if m == n else 0) )/(2*diag[r]) for n in range(2)] for m in range(2)] for r in range(2)]
            tangent = [1, v]
            accel = [sum(gamma[r][m][n]*tangent[m]*tangent[n] for m in range(2) for n in range(2))+(acc if r == 1 else 0) for r in range(2)]
            normal = [sign*v/(ell*mp.sqrt(1-v*v)), sign/(ell*mp.sqrt(1-v*v))]
            curvature = sum(diag[r]*normal[r]*accel[r] for r in range(2))/(ell**2*(1-v*v))
            displayed = sign*(acc/(1-v*v)+v*wt+ws-mp.cot(s))/(ell*mp.sqrt(1-v*v))
            close(curvature-displayed)
            p, pt, ps, _, _, _ = fields.jets(fields.phi, t, s)
            nphi = normal[0]*pt+normal[1]*ps
            nw = normal[0]*wt+normal[1]*ws
            geometric = ell*mp.sqrt(1-v*v)*(nphi-(curvature-nw)*p)+sign*(acc*p/(1-v*v)-v*pt)
            q = sign*(ps+mp.cot(s)*p)
            close(geometric-q)
            radius = lambda tt, ss: mp.sqrt(fields.X(tt, ss))
            nradius = normal[0]*mp.diff(radius,(t,s),(1,0))+normal[1]*mp.diff(radius,(t,s),(0,1))
            close(nphi-radius(t,s)*nradius/fields.G)
            rows.append(dict(sign=sign,t=rawt,velocity=str(v),curvature_error=str(abs(curvature-displayed)),charge_error=str(abs(geometric-q)),full_four_boundary_trace=str(curvature+2*nradius/radius(t,s))))
    return rows


def finite_cut_and_weighted_checks(fields):
    rows = []
    eps = mp.mpf('.25')
    M0, Mt, M1 = map(mp.mpf, ['.5', '.08', '.5'])
    W0, Wt, W1 = map(mp.mpf, ['.1', '.07', '.1'])
    mx = 1-fields.mu*M0/mp.sin(eps)
    om = fields.mu*W0/mp.sin(eps)
    A = Mt**2/(4*mx)+Wt*Mt
    D = M1**2/(4*mx)+W1*M1+2*mp.exp(2*om)*W0*M0+fields.vartheta*mp.exp(2*om)*M0**2/mx
    i2 = 2*mp.cot(eps)
    i4 = 2*mp.cot(eps)+mp.mpf(2)/3*mp.cot(eps)**3
    bound = fields.mu*(A*i2+D*i4)
    # Global analytic jet bounds, using only |sin|,|cos|<=1:
    # .45+.02, .06+.012, .45+.04; .08*1.2,
    # .08*.8, .08*1.2. These dominate all six weighted jets.
    analytic = list(map(mp.mpf,['.47','.072','.49','.096','.064','.096']))
    constants = [M0,Mt,M1,W0,Wt,W1]
    assert all(a < b for a,b in zip(analytic,constants)) and mx > 0
    c = mp.mpf('.04')
    pchiral = lambda v:-c*(1+mp.mpf('.1')*mp.cos(2*v))/2
    tzero = lambda t,s:pchiral(t-s)+pchiral(t+s-mp.pi)
    E0 = -c*mp.pi
    gfree = -E0/fields.B
    z = lambda t:mp.mpf('.3')+mp.mpf('.02')*mp.sin(mp.mpf('.8')*t)
    for rawt in ['0','.7','2','7']:
        t = mp.mpf(rawt)
        left = mp.mpf('.3')+mp.mpf('.015')*mp.sin(mp.mpf('.7')*t)
        right = mp.pi-mp.mpf('.3')-mp.mpf('.01')*mp.cos(mp.mpf('.5')*t)
        def q(s,sign):
            return sign*(mp.diff(fields.phi,(t,s),(0,1))+mp.cot(s)*fields.phi(t,s))
        qsum = q(left,-1)+q(right,1)
        Emiss = mp.quad(lambda s:tzero(t,s),[0,left])+mp.quad(lambda s:tzero(t,s),[right,mp.pi])
        beta = (qsum-2*fields.B*(z(t)+mp.diff(z,t,2))-Emiss)/fields.B
        corr = mp.quad(lambda s:(fields.orbit(t,s)[2][0][0]-tzero(t,s)+fields.orbit(t,s)[1][0][0])/fields.B,[left,(left+right)/2,right])
        actual = mp.diff(z,t,2)+z(t)-gfree/2
        reconstructed = -(corr+beta)/2
        close(actual-reconstructed)
        wrong = -(corr+beta+Emiss/fields.B)/2
        close(wrong-actual+Emiss/(2*fields.B))
        assert abs(wrong-actual) > mp.mpf('1e-4'), 'collar omission did not test a nonzero error'
        # Split the positive absolute integral at pi/2. It is used only for
        # an inequality with ample slack, not an equality certificate.
        measured = mp.quad(lambda s:abs(fields.orbit(t,s)[1][0][0])/fields.B,[eps,mp.pi/2,mp.pi-eps])
        assert measured < bound
        for j in range(13):
            s = eps+(mp.pi-2*eps)*j/12
            psi,pt,ps,_,_,_ = fields.jets(fields.psi,t,s)
            w,wt,ws,_,_,_ = fields.jets(fields.omega,t,s)
            sampled = [abs(mp.sin(s)*psi),abs(mp.sin(s)*pt),abs(mp.sin(s)**2*ps),
                       abs(mp.sin(s)*w/fields.mu),abs(mp.sin(s)*wt/fields.mu),abs(mp.sin(s)**2*ws/fields.mu)]
            assert all(a <= b for a,b in zip(sampled,constants))
        rows.append(dict(t=rawt,finite_identity_error=str(abs(actual-reconstructed)),free_collar_energy=str(Emiss),
                         omission_error=str(wrong-actual),beta=str(beta),synthetic_residual=str(actual),
                         measured_S_G=str(measured),weighted_upper_bound=str(bound),slack=str(bound-measured)))
    close(mp.quad(lambda s:1/mp.sin(s)**2,[eps,mp.pi/2,mp.pi-eps])-i2)
    close(mp.quad(lambda s:1/mp.sin(s)**4,[eps,mp.pi/2,mp.pi-eps])-i4)
    return dict(analytic_global_jet_bounds=list(map(str,analytic)),accepted_constants=list(map(str,constants)),
                mX=str(mx),Omega=str(om),A=str(A),D=str(D),I2=str(i2),I4=str(i4),cases=rows,
                beta_is_not_independently_bounded=True,synthetic_reference_not_physical_experiment=True)


def candidate_checks():
    k = mp.mpf(2000)
    # Stable positive quadratic root of xi(1+3xi)=pi*k^6.
    xi = 2*mp.pi*k**6/(1+mp.sqrt(1+12*mp.pi*k**6))
    close((xi*(1+3*xi)/(mp.pi*k**6))-1)
    assert xi >= k**3
    mu = 4/(xi*k**7)
    upper = 4/k**10
    assert mu <= upper
    e = mp.mpf('1e-5')
    mx = 1-upper/mp.sin(mp.atan(e))
    E2 = mp.exp(2*upper/mp.sin(mp.atan(e)))
    assert mx >= mp.mpf('.99') and E2 <= mp.mpf('1.01')
    A = 1/(4*mx)+1
    D = 1/(4*mx)+1+2*E2+E2/mx
    assert A < mp.mpf('1.3') and D < mp.mpf('4.4')
    bound = upper*(A*2/e+D*(2/e+mp.mpf(2)/(3*e**3)))
    assert bound < mp.mpf('1.2e-17')
    K = mp.mpf('1e5')
    mxK = 1-upper*K/mp.sin(mp.atan(e))
    E2K = mp.exp(2*upper*K/mp.sin(mp.atan(e)))
    assert mxK >= mp.mpf('.99') and E2K <= mp.mpf('1.01')
    AK = K*K/(4*mxK)+K*K
    DK = K*K/(4*mxK)+K*K+2*E2K*K*K+E2K*K*K/mxK
    assert AK < mp.mpf('1.3')*K*K and DK < mp.mpf('4.4')*K*K
    boundK = upper*(AK*2/e+DK*(2/e+mp.mpf(2)/(3*e**3)))
    assert boundK < mp.mpf('1.2e-7')
    rows=[]
    for zi in [mp.mpf(1)/4,mp.mpf(7)/25]:
        lam = mu*zi/2
        olde = lam*k**4
        oldratio = mu/olde**3
        close(oldratio/(xi**2*k**2/(2*zi**3))-1)
        assert olde <= mp.mpf('8.75e-21') and oldratio > mp.mpf('5.83e27')
        delta = lam/e
        assert delta <= mp.mpf('5.46875e-29') and delta > lam
        close(lam/delta-e)
        rows.append(dict(z_initial=str(zi),lambda0=str(lam),old_e=str(olde),old_mu_over_e3=str(oldratio),new_delta=str(delta)))
    return dict(k=2000,eta='1/2000',xi=str(xi),mu=str(mu),mu_upper_bound=str(upper),new_e=str(e),
                new_mu_over_e3_upper_bound=str(upper/e**3),conditional_S_G_bound=str(bound),
                conditional_K_window=str(K),conditional_K_window_S_G_bound=str(boundK),
                thermal_static_mass_envelope=str(256*e**2),conditional_jet_hypotheses_established=False,cases=rows)


def main():
    if sys.argv[1:] == ['--contact-transfer-only']:
        contact_transfer_main()
        return
    assert not sys.argv[1:], 'unrecognized arguments'
    mp.mp.dps = DPS
    hashes = {name:sha(ROOT/name) for name in SOURCES}
    result = dict(story='S136.3',started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  dps=DPS,tolerances=TOL,input_sha256=hashes,
                  hash_policy='SHA-256 after CRLF-to-LF normalization; no other byte changes',
                  execution=dict(executable=sys.executable,argv=sys.orig_argv,cwd=str(Path.cwd()),command=subprocess.list2cmdline(sys.orig_argv)),
                  runtime=dict(python=sys.version,platform=platform.platform(),sympy=sp.__version__,mpmath=mp.__version__))
    if OUTPUT.exists():
        previous=json.loads(OUTPUT.read_text(encoding='utf-8'))
        if previous.get('status')=='failed':
            record('Retained previous failed receipt: '+json.dumps(previous))
    record('S136.3 physical remainder check: '+result['started_utc'])
    record('Execution: '+json.dumps(result['execution']))
    record('Frozen acceptance: '+json.dumps(dict(dps=DPS,tolerances=TOL)))
    try:
        result['symbolic']=symbolic_checks()
        record('Exact potential and moving-cut charge: PASS')
        fields=Fields()
        result['four_dimensional_tensor']=full_tensor_checks(fields)
        record('Independent four-dimensional orbit and angular Einstein equations: PASS')
        result['moving_cuts']=moving_cut_checks(fields)
        record('Moving-cut acceleration curvature and radius-normal charge: PASS')
        result['finite_cut_and_weighted_norm']=finite_cut_and_weighted_checks(fields)
        record('Finite-cut collar identity and weighted geometric bound: PASS')
        result['candidate_conditioning']=candidate_checks()
        record('Frozen candidate: old and conditioned internal seams: PASS')
        assert hashes == {name:sha(ROOT/name) for name in SOURCES},'input changed during run'
        result['status']='pass_exact_physical_remainder_identity_and_conditional_bound'
        result['claims']=dict(exact_spherical_identity=True,physical_stress_bound_established=False,
                              physical_boundary_defect_bound_established=False,weighted_jets_established_for_candidate=False,
                              matched_dynamical_clock_established=False,physical_residual_certified=False,
                              robustness_rectangle_certified=False,source_changed=False,physical_opening_achieved=False,
                              s074_gate_closed=False,lean_invoked=False)
        result['summary']=dict(four_dimensional_points=len(result['four_dimensional_tensor']),
                               moving_cut_cases=len(result['moving_cuts']),finite_cut_cases=4,
                               maximum_orbit_tensor_error=str(max(mp.mpf(error) for row in result['four_dimensional_tensor'] for error in row['orbit_component_errors'])),
                               maximum_angular_error=str(max(mp.mpf(row['angular_error']) for row in result['four_dimensional_tensor'])),
                               maximum_finite_identity_error=str(max(mp.mpf(row['finite_identity_error']) for row in result['finite_cut_and_weighted_norm']['cases'])))
        write(OUTPUT,json.dumps(result,indent=2)+'\n')
        record('Summary: '+json.dumps(result['summary']))
        record('Normalized input hashes: '+json.dumps(hashes))
        record('Terminal status: PASS, exit 0')
    except Exception:
        result['status']='failed'
        result['failure_traceback']=traceback.format_exc()
        write(OUTPUT,json.dumps(result,indent=2)+'\n')
        record(result['failure_traceback'].rstrip())
        record('Failure input hashes: '+json.dumps(hashes))
        record('Terminal status: FAIL, exit 1; retained in append-only log')
        raise


def contact_transfer_checks():
    c, left, right, sigma = sp.symbols('c left right sigma', real=True)
    qc = c*(-sp.cot(left)+sp.cot(right))
    anti = -c*sp.cot(sigma)
    assert sp.trigsimp(sp.diff(anti,sigma)-c/sp.sin(sigma)**2) == 0
    integral = anti.subs(sigma,right)-anti.subs(sigma,left)
    assert sp.simplify(integral+qc) == 0
    B, draw, beta = sp.symbols('B draw beta', nonzero=True)
    dshift = draw-integral
    bshift = beta-qc/B
    assert sp.simplify(dshift/B+bshift-draw/B-beta) == 0
    rows=[]
    for rawleft, rawright in [('.2','2.8'),('.31','2.72'),('.00001','3.14158')]:
        l,r=mp.mpf(rawleft),mp.mpf(rawright)
        cc=mp.mpf('.003')
        bb=mp.mpf('.2')
        q=cc*(-mp.cot(l)+mp.cot(r))
        integral_exact=cc*(mp.cot(l)-mp.cot(r))
        close(integral_exact+q)
        rawdelta=mp.mpf('.013')
        rawbeta=mp.mpf('-.07')
        shifted_delta=rawdelta-integral_exact
        shifted_beta=rawbeta-q/bb
        rawr=-(rawdelta/bb+rawbeta)/2
        shiftedr=-(shifted_delta/bb+shifted_beta)/2
        close(rawr-shiftedr)
        # Check the local tensor subtraction independently of the integral.
        # L0_tt c=+c csc^2, L0_ss c=-c csc^2, L0_ts c=0.
        mid=(l+r)/2
        p=lambda s:cc
        local=-mp.diff(p,mid,2)-mp.cot(mid)*mp.diff(p,mid)+p(mid)/mp.sin(mid)**2
        close(local-cc/mp.sin(mid)**2)
        rows.append(dict(left=rawleft,right=rawright,c=str(cc),q_c=str(q),
                         contact_integral=str(integral_exact),raw_residual=str(rawr),
                         shifted_residual=str(shiftedr),invariance_error=str(abs(rawr-shiftedr))))
    return dict(symbolic_contact_integral='exact zero',symbolic_pair_invariance='exact zero',cases=rows,
                raw_area_kept_in_geometric_remainder=True,constant_is_fixed_not_fitted=True)


def contact_transfer_main():
    """Validate only the additive contact change; retain prior full run."""
    mp.mp.dps=DPS
    previous=json.loads(OUTPUT.read_text(encoding='utf-8'))
    if previous.get('retained_full_validation'):
        if previous['status']=='failed':
            record('Retained prior failed targeted receipt: '+json.dumps(previous))
        previous=previous['retained_full_validation']
    assert previous['status']=='pass_exact_physical_remainder_identity_and_conditional_bound'
    hashes={name:sha(ROOT/name) for name in SOURCES}
    result=dict(story='S136.3',started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                dps=DPS,tolerances=TOL,input_sha256=hashes,
                hash_policy='SHA-256 after CRLF-to-LF normalization; no other byte changes',
                execution=dict(executable=sys.executable,argv=sys.orig_argv,cwd=str(Path.cwd()),command=subprocess.list2cmdline(sys.orig_argv)),
                runtime=dict(python=sys.version,platform=platform.platform(),sympy=sp.__version__,mpmath=mp.__version__),
                retained_full_validation=previous,
                validation_scope='Prior full tensor/cut/quadrature run retained; this run verifies new contact transfer and unchanged checker bodies, without repeating expensive checks.')
    record('S136.3 targeted contact-transfer check: '+result['started_utc'])
    record('Execution: '+json.dumps(result['execution']))
    record('Frozen acceptance unchanged: '+json.dumps(dict(dps=DPS,tolerances=TOL)))
    try:
        prefix='Prior full-run checker source snapshot: '
        snapshots=[json.loads(line[len(prefix):]) for line in LOG.read_text(encoding='utf-8').splitlines() if line.startswith(prefix)]
        original=next(source for source in snapshots if hashlib.sha256(source.encode('utf-8')).hexdigest()==previous['input_sha256']['checks/physical_clock_remainder.py'])
        current=(ROOT/'checks/physical_clock_remainder.py').read_text(encoding='utf-8')
        def definitions(source):
            return {node.name:ast.dump(node,include_attributes=False) for node in ast.parse(source).body
                    if isinstance(node,(ast.FunctionDef,ast.ClassDef)) and node.name!='main'}
        olddefs,newdefs=definitions(original),definitions(current)
        assert all(newdefs[name]==body for name,body in olddefs.items()),'original non-main checker body changed'
        oldmain=next(node for node in ast.parse(original).body if isinstance(node,ast.FunctionDef) and node.name=='main')
        newmain=next(node for node in ast.parse(current).body if isinstance(node,ast.FunctionDef) and node.name=='main')
        assert [ast.dump(node,include_attributes=False) for node in oldmain.body]==[ast.dump(node,include_attributes=False) for node in newmain.body[2:]], 'full-run main body changed beyond mode dispatch'
        def assignments(source):
            return [ast.dump(node,include_attributes=False) for node in ast.parse(source).body if isinstance(node,ast.Assign)]
        assert assignments(original)==assignments(current),'original checker constants changed'
        unchanged={name:hashlib.sha256(body.encode('utf-8')).hexdigest() for name,body in olddefs.items()}
        for name in SOURCES:
            if name not in ['notes/physical-clock-remainder.md','checks/physical_clock_remainder.py']:
                assert hashes[name]==previous['input_sha256'][name], 'prior full-run dependency changed: '+name
        result['retained_run_integrity']=dict(original_checker_sha256=previous['input_sha256']['checks/physical_clock_remainder.py'],
                                              original_source_preserved_in_log=True,unchanged_definition_AST_sha256=unchanged,
                                              full_run_main_body_unchanged_except_mode_dispatch=True,
                                              original_global_assignments_unchanged=True,all_other_full_run_source_hashes_unchanged=True)
        result['contact_transfer']=contact_transfer_checks()
        assert hashes=={name:sha(ROOT/name) for name in SOURCES},'input changed during targeted run'
        result['status']='pass_contact_transfer_with_retained_full_validation'
        result['claims']=previous['claims']
        result['summary']=dict(retained_full_run=previous['summary'],new_contact_cases=3,
                               maximum_contact_invariance_error=str(max(mp.mpf(row['invariance_error']) for row in result['contact_transfer']['cases'])),
                               physical_residual_certified=False)
        write(OUTPUT,json.dumps(result,indent=2)+'\n')
        record('Contact transfer, retained source hashes and unchanged original AST bodies: PASS')
        record('Summary: '+json.dumps(result['summary']))
        record('Current normalized input hashes: '+json.dumps(hashes))
        record('Terminal status: PASS, exit 0; full earlier run retained, no expensive rerun claimed')
    except Exception:
        result['status']='failed'
        result['failure_traceback']=traceback.format_exc()
        write(OUTPUT,json.dumps(result,indent=2)+'\n')
        record(result['failure_traceback'].rstrip())
        record('Failure input hashes: '+json.dumps(hashes))
        record('Terminal status: FAIL, exit 1; retained in append-only log')
        raise


if __name__ == '__main__':
    main()
