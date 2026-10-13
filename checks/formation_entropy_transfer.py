"""S071: conditional entropy/throat hierarchy compatibility identity."""
from pathlib import Path
import hashlib
import json
import sympy as s

G, q, g, Nf, ds = s.symbols('G q g Nf delta_s', positive=True)
re2 = s.pi*q*q*G/(g*g)
S0 = s.pi*re2/G
c = q*Nf
length_ratio = s.simplify(16*re2/(G*q*Nf))  # ell/r_e
fraction = s.simplify(2*c*ds/S0)
assert s.simplify(fraction-(32*ds/s.pi)/length_ratio) == 0
data = {'date': '2026-10-12', 'status': 'conditional_scaling_identity_verified',
        'sympy': s.__version__, 'code_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'conventions': 'Nf unit-charge four-dimensional Dirac fields; c=Nf*q; two Majoranas per 2D Dirac',
        'S0_per_mouth': str(s.simplify(S0)), 'ell_over_re': str(length_ratio),
        'hypothetical_entropy_fraction': str(fraction),
        'identity': 'DeltaS/S0 = (32 delta_s/pi) * r_e/ell',
        'assumption': 'delta_s is a uniformly bounded per-Majorana disk entropy change, independent of the hierarchy',
        'conclusion': 'r_e/ell -> 0 implies DeltaS/S0 -> 0 under the stated assumption',
        'not_computed': ['actual MMP boundary entropy change', 'interacting or higher-Landau-level effects',
                         'formation rate', 'universal topology-change bound', 'nonlinear wormhole stability']}
(Path(__file__).resolve().parents[1]/'evidence/formation-entropy-transfer.json').write_text(
    json.dumps(data, indent=2)+'\n', encoding='utf-8', newline='\n')
print('PASS: exact entropy-fraction / throat-hierarchy identity')
