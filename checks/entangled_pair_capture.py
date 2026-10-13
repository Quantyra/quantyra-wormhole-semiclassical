"""S072: constrained JT reduced capture under a monotone smooth ramp.

No microscopic switching radiation or four-dimensional transfer included.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import scipy
from scipy.integrate import solve_ivp
import sympy as sp

z, w, a = sp.symbols('z w a', real=True)
energy = w*w+z*z-a*z
adot = sp.symbols('adot', real=True)
derivative = sp.diff(energy,z)*z*w+sp.diff(energy,w)*(-z*z+a*z/2)+sp.diff(energy,a)*adot
assert sp.simplify(derivative+adot*z) == 0
# Rational bounds used in the ordinary proof, including the largest endpoint.
assert sp.Rational(9,128)*sp.Rational(64,63) == sp.Rational(1,14)
assert max(x*x-sp.Rational(13,14)*x for x in (sp.Rational(1,4),sp.Rational(3,4))) == -sp.Rational(15,112)


def ramp(t, duration, kind):
    x = np.clip(t/duration, 0, 1)
    if kind == 'linear':
        return x, 1/duration if t < duration else 0
    if kind == 'cubic':
        return x*x*(3-2*x), 6*x*(1-x)/duration
    return (1-np.cos(np.pi*x))/2, np.pi*np.sin(np.pi*x)/(2*duration)


rows = []
for kind in ('linear', 'cubic', 'cosine'):
    for duration in (.05, .25, .5):
        for z0 in (.25, .5, .75):
            def rhs(t, v):
                zz, ww, work = v
                aa, da = ramp(t, duration, kind)
                return [zz*ww, -zz*zz+aa*zz/2, -da*zz]
            sol = solve_ivp(rhs, (0, duration), [z0, 0, 0], method='DOP853',
                            rtol=2e-11, atol=2e-13, max_step=duration/40)
            assert sol.success
            zf, wf, work = sol.y[:, -1]
            ef = wf*wf+zf*zf-zf
            balance = abs(ef-(z0*z0+work))
            assert balance < 3e-10
            assert ef <= -15/112+1e-10
            lower, upper = .5-np.sqrt(.25+ef), .5+np.sqrt(.25+ef)
            def autonomous(t, v):
                zz, ww = v
                return [zz*ww, -zz*zz+zz/2]
            follow = solve_ivp(autonomous, (0, 150), [zf, wf], method='DOP853',
                               rtol=2e-11, atol=2e-13, max_step=.1)
            assert follow.success
            zz, ww = follow.y
            drift = float(np.max(np.abs(ww*ww+zz*zz-zz-ef)))
            assert drift < 3e-9
            assert min(zz) >= lower-1e-9 and max(zz) <= upper+1e-9
            rows.append({'ramp': kind, 'duration_eta_times_u': duration, 'initial_y_over_eta': z0,
                         'final_energy_over_N_eta_squared': float(ef),
                         'work_balance_error': float(balance), 'autonomous_energy_drift': drift,
                         'turning_interval_y_over_eta': [float(lower), float(upper)]})
out = {'date': '2026-10-12', 'status': 'reduced_capture_theorem_verified',
       'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'numpy': np.__version__, 'scipy': scipy.__version__, 'sympy': sp.__version__,
       'theorem': 'For 1/4<=y0/eta_f<=3/4, phi_dot(0)=0, and any absolutely continuous monotone ramp 0 to eta_f over eta_f*Delta_u<=1/2, e_final<=-15/112 in the stated reduced model.',
       'energy_identity': 'd(w^2+z^2-a*z)/ds = -a_prime*z',
       'uniform_y_over_eta_lower_bound': float(.5-np.sqrt(13/112)),
       'numeric_cases': rows,
       'not_included': ['microscopic quench particle production', 'nonzero matter SL(2) charges',
                        'quantum state overlap errors', 'four-dimensional boundary matching',
                        'initial entangled-pair preparation', 'nonspherical stability']}
(Path(__file__).resolve().parents[1]/'evidence/entangled-pair-capture.json').write_text(
    json.dumps(out, indent=2)+'\n', encoding='utf-8', newline='\n')
print('PASS: reduced work identity, rational capture margin and 27 switching/post-switch dynamics checks')
