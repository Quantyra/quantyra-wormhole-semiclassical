"""Exact algebra supporting a literature/domain comparison, not a gravity solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path') + 1])
import sympy as s

root = Path(__file__).resolve().parents[1]
t, u, rate, strength, cutoff, w = s.symbols('t u gamma Gamma Omega omega', positive=True)
variance = strength * (t - (1 - s.exp(-rate*t))/rate)
integrated = s.integrate(strength*rate*(t-u)*s.exp(-rate*u), (u, 0, t))
assert s.simplify(integrated - variance) == 0
assert s.simplify(s.diff(variance, t) - strength*(1-s.exp(-rate*t))) == 0
spectral = strength*rate**2/(rate**2+w**2)
derivative_cost = strength*rate**2/s.pi*(cutoff-rate*s.atan(cutoff/rate))
assert s.simplify(s.diff(derivative_cost, cutoff) - cutoff**2*spectral.subs(w, cutoff)/s.pi) == 0
assert s.simplify(s.limit(derivative_cost/cutoff, cutoff, s.oo) - strength*rate**2/s.pi) == 0
sigma2 = s.symbols('sigma2', positive=True)
assert s.simplify(s.limit(variance.subs(strength, 2*sigma2/rate), rate, 0) - sigma2*t**2) == 0
lower = s.Matrix([[0, 1], [0, 0]])
collective = s.kronecker_product(lower, s.eye(2)) + s.kronecker_product(s.eye(2), lower)
singlet = s.Matrix([0, 1, -1, 0])/s.sqrt(2)
assert collective*singlet == s.zeros(4, 1)
assert collective.T*singlet == s.zeros(4, 1)
paths = ['notes/current-noise-entanglement-lead-review.md',
         'checks/current_noise_entanglement_transfer.py', 'notes/coupled-fermion-domain.md']
record = {
    'executed_utc': datetime.now(timezone.utc).isoformat(),
    'paper': 'https://arxiv.org/html/2610.31947v1',
    'source_sha256': {p: sha256((root/p).read_bytes()).hexdigest() for p in paths},
    'symbolic_results': {
        'integrated_phase_variance': str(variance),
        'variance_derivative': str(s.simplify(s.diff(variance, t))),
        'OU_phase_derivative_cost_with_physical_bandwidth': str(derivative_cost),
        'linear_divergence_coefficient': str(strength*rate**2/s.pi),
        'fixed_variance_long_memory_limit': str(sigma2*t**2),
        'collective_raising_and_lowering_annihilate_singlet': True,
    },
    'scope': {
        'ordinary_argument_supported_by_exact_algebra': True,
        'full_paper_independently_verified': False,
        'reduced_qubit_entanglement_revivals_disproved': False,
        'direct_ideal_OU_boundary_transfer_meets_existing_hypotheses': False,
        'smooth_noise_construction_ruled_out': False,
        'wormhole_stress_or_stability_bound_obtained': False,
        'physical_budget_improvement': 0,
    },
}
out = root/'evidence/current-noise-entanglement-transfer.json'
assert not out.exists()
out.write_bytes((json.dumps(record, indent=2)+'\n').encode())
print(json.dumps(record['symbolic_results'], indent=2))
