"""S071: fixed-background free boundary coupling and regulated entropy.

The gauge field is frozen. This is not the interacting MMP determinant.
"""
from pathlib import Path
import hashlib
import json
import math
import mpmath as mp
import sympy as s

h, C, sign = s.symbols('h C sign', real=True)
sx = s.Matrix([[0, 1], [1, 0]])
U = (s.eye(2)-s.I*h*sx)*(s.eye(2)+s.I*h*sx).inv()
assert s.simplify(U.conjugate().T*U-s.eye(2)) == s.zeros(2)
assert s.simplify(s.Abs(U[0, 1])**2-4*h*h/(1+h*h)**2) == 0
M = -s.I*s.eye(2)/C+s.Matrix([[0, s.I*h], [-s.I*h, 0]])
assert s.simplify(M.det()/M.subs(h, 0).det()-(1+C*C*h*h)) == 0
a = s.symbols('a', positive=True)
F = 1/(2*s.sinh(a))
entropy = s.simplify(F+a*s.diff(F, a))
assert s.limit(entropy, a, 0) == 0
assert s.limit(entropy/a, a, 0) == -s.Rational(1, 6)
mp.mp.dps = 45
rows = []
for beta_cutoff in (10, 30, 100, 300):
    aa = mp.pi/beta_cutoff
    last = math.ceil(35/float(aa))
    qsum = mp.fsum(mp.exp(-(2*n+1)*aa) for n in range(last))
    ssum = mp.fsum(mp.exp(-(2*n+1)*aa)*(1-(2*n+1)*aa) for n in range(last))
    for hh in (mp.mpf('.1'), mp.mpf('1'), mp.mpf('10')):
        A = mp.log1p(hh*hh)
        target = A/(2*mp.sinh(aa))
        starget = target*(1-aa/mp.tanh(aa))
        error = abs(A*qsum-target)
        serror = abs(A*ssum-starget)
        assert max(error, serror) < mp.mpf('1e-25')
        rows.append({'h': float(hh), 'beta_times_cutoff': beta_cutoff,
                     'positive_modes': last, 'logZ_absolute_error': float(error),
                     'entropy_absolute_error': float(serror),
                     'relative_entropy': float(starget)})
out = {'date': '2026-10-12', 'status': 'free_contact_only_verified',
       'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'sympy': s.__version__, 'mpmath': mp.__version__,
       'boundary_dimension': '1/2', 'contact_transmission': '4 h^2/(1+h^2)^2',
       'regulated_logZ_per_real_pair': 'log(1+lambda^2)/(2 sinh(pi/(beta Lambda)))',
       'relative_thermal_entropy_limit': 0,
       'checks': rows,
       'limitations': ['frozen gauge and geometry', 'no boundary zero mode or extra state sector',
                       'fixed antiperiodic thermal spin structure', 'no claim for arbitrary UV boundary constants',
                       'not the MMP interacting infrared or four-dimensional formation rate']}
(Path(__file__).resolve().parents[1]/'evidence/formation-boundary-flow.json').write_text(
    json.dumps(out, indent=2)+'\n', encoding='utf-8', newline='\n')
print('PASS: unitary contact, determinant identity, 12 regulated thermal sums and entropy limit')
