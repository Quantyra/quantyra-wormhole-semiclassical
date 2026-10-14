"""Complete mixed dipole/l=3 vacuum pair kernel in the stated cubic representative.

Analytic tails and a closed quadratic form, not a full physical stability claim.
"""
from pathlib import Path
import argparse
import hashlib
import json
import platform

import mpmath as mp
import numpy as np
import sympy as s


def algebra():
    x = s.symbols('x', nonnegative=True)
    Dlower = 3+x+s.Rational(20, 13)*x*x
    gap = s.expand(Dlower*(Dlower-1)-(2+12*x))
    positive = s.Rational(5, 2)+s.Rational(41, 452)+s.Rational(113, 13)*(x-s.Rational(91, 226))**2+s.Rational(40, 13)*x**3+s.Rational(400, 169)*x**4
    assert s.expand(gap-positive) == 0
    assert s.expand(9+16*x-(3+2*x)**2-4*x*(1-x)) == 0
    assert s.expand((1+6*x)**2-(1+12*x+12*x*x)) == 24*x*x
    assert s.expand((1+12*x)**2-(1+24*x+24*x*x)) == 120*x*x
    assert s.expand((1+24*x+24*x*x)-(1+2*x)**2-20*x*(1+x)) == 0
    assert s.Rational(17, 24)**2 > s.Rational(1, 2)
    assert s.Rational(17, 24)+s.Rational(19, 32) == s.Rational(125, 96)
    assert s.Rational(125, 96) < s.Rational(4, 3)
    return {'gap_lower': '5/2+41/452+(113/13)(x-91/226)^2+(40/13)x^3+(400/169)x^4',
            'beta_upper': '125/96 < 4/3', 'beta_prime_upper': '8/15',
            'tail_after_N_terms': '19/[32*pi*(N-1)^2], N>=2',
            'scope': 'x=1/C in [0,1]; both axial branches; z between zero and the larger polar mass'}


def rho(a, b, j):
    j = mp.mpf(j)
    return (2*(a+b)+4*j-1)/(4*mp.pi)*mp.gamma(j+mp.mpf('.5'))*mp.gamma(j+a)*mp.gamma(j+b)*mp.gamma(j+a+b-mp.mpf('.5'))/(mp.gamma(j+1)*mp.gamma(j+a+mp.mpf('.5'))*mp.gamma(j+b+mp.mpf('.5'))*mp.gamma(j+a+b))


def weights(a, b, N):
    out = [rho(a, b, 0)]
    for j in range(N-1):
        ratio = (2*(a+b)+4*j+3)/(2*(a+b)+4*j-1)
        ratio *= (j+mp.mpf('.5'))*(j+a)*(j+b)*(j+a+b-mp.mpf('.5'))/((j+1)*(j+a+mp.mpf('.5'))*(j+b+mp.mpf('.5'))*(j+a+b))
        out.append(out[-1]*ratio)
    return out


def central_weight(a, n):
    return mp.gamma(n+mp.mpf('.5'))*mp.gamma(n+a)/(2*mp.pi*mp.gamma(n+1)*mp.gamma(n+a+mp.mpf('.5')))


def representation_checks():
    rows, kl = [], []
    for aa, bb in [('2', '1'), ('2.3', '1.6'), ('3', '3'), ('3', '5')]:
        a, b = mp.mpf(aa), mp.mpf(bb)
        for j in range(9):
            N = 2*j
            primary = [mp.mpf(1)]
            for p in range(N):
                primary.append(-primary[-1]*mp.sqrt((N-p)*(N-p-1+2*b)/((p+1)*(p+2*a))))
            norm = mp.sqrt(mp.fsum(v*v for v in primary))
            primary = [v/norm for v in primary]
            for p in range(N):
                err = primary[p+1]*mp.sqrt((p+1)*(p+2*a))+primary[p]*mp.sqrt((N-p)*(N-p-1+2*b))
                assert abs(err) < mp.mpf('1e-65')
            # Distinct fields: no identical-field Wick factor or sqrt(2).
            amp = (-1)**j*mp.fsum(primary[p]*mp.sqrt(central_weight(a, p//2)*central_weight(b, (N-p)//2)) for p in range(0, N+1, 2))
            actual = amp*amp/central_weight(a+b+N, 0)
            relative = abs(actual/rho(a, b, j)-1)
            assert relative < mp.mpf('1e-65')
            rows.append({'a': aa, 'b': bb, 'j': j, 'relative_weight_error': mp.nstr(relative, 8)})
        def W(h, t):
            return central_weight(h, 0)*mp.exp(-h*t)*mp.hyp2f1(mp.mpf('.5'), h, h+mp.mpf('.5'), mp.exp(-2*t))
        ws = weights(a, b, 100)
        for t in map(mp.mpf, ['.45', '.8', '1.2']):
            left = W(a, t)*W(b, t)
            right = mp.fsum(ws[j]*W(a+b+2*j, t) for j in range(100))
            error = abs(right/left-1)
            assert error < mp.mpf('1e-35')
            kl.append({'a': aa, 'b': bb, 'Euclidean_time': str(t), 'terms': 100, 'relative_error': mp.nstr(error, 8)})
    return {'independent_distinct_field_primary_weights': rows, 'Euclidean_product_checks': kl}


def physical(C, gamma, N):
    C, gamma = mp.mpf(C), mp.mpf(gamma)
    rp, ra = mp.sqrt(C*C+12*C+12), mp.sqrt(C*C+24*C+24)
    masses = [24/(C*(C+6+rp)), (C+6+rp)/C]
    axial = [120/(C*(C+12+ra)), (C+12+ra)/C]
    cp, sp = mp.sqrt((rp+C)/(2*rp)), mp.sqrt(12*(C+1))/mp.sqrt(2*rp*(rp+C))
    ca, sa = mp.sqrt((ra+C+2)/(2*ra)), mp.sqrt(20*(C+1))/mp.sqrt(2*ra*(ra+C+2))
    Up, Ua = mp.matrix([[-sp, cp], [cp, sp]]), mp.matrix([[-sa, ca], [ca, sa]])
    B = mp.matrix([[-2*mp.sqrt(42)*(2*C*C+22*C+27)/(35*C), -4*(C+2)*mp.sqrt(210*(C+1))/(35*C)], [-12*(C+8)*mp.sqrt(14*(C+1))/(35*C), -6*mp.sqrt(70)*(2*C+1)/(35*C)]])
    gc = mp.matrix(2, 2)
    for p in range(2):
        for branch in range(2):
            gc[p, branch] = mp.sqrt(8*mp.pi*gamma/C)*(Up[:, p].T*B*Ua[:, branch])[0]/mp.sqrt(4*mp.pi*C*masses[p]/5)
    # Compress the equivalent magnetic channels at fixed external m.
    gc *= mp.sqrt(mp.mpf(7)/3)
    q = gamma*(53+830/C+944/C**2)/15
    assert abs(mp.fsum(gc[p, branch]**2/masses[p] for p in range(2) for branch in range(2))/q-1) < mp.mpf('1e-60')
    delta = (1+mp.sqrt(9+16/C))/2
    hs = [(1+mp.sqrt(1+4*nu))/2 for nu in axial]
    baths = []
    for h in hs:
        w = weights(delta, h, N)
        lambdas = [(delta+h+2*j)*(delta+h+2*j-1) for j in range(N)]
        baths.append((w, lambdas))
    return masses, gc, q, delta, hs, baths


def norm2(matrix):
    return (matrix[0, 0]+matrix[1, 1]+mp.sqrt((matrix[0, 0]-matrix[1, 1])**2+4*matrix[0, 1]**2))/2


def spectral_checks():
    rows = []
    k = mp.mpf(2000)
    candidateC = mp.sqrt(1+12*mp.pi*k**6)
    N = 256
    tail = 19/(32*mp.pi*(N-1)**2)
    candidate = None
    for C in [mp.mpf(1), mp.mpf(4), mp.mpf(100), mp.mpf('1e12'), candidateC]:
        gamma = 6*C/(k**10*(C-1)) if C == candidateC else mp.mpf('.001')
        masses, g, q, delta, hs, baths = physical(C, gamma, N)
        epsilon = 4*q/3
        for branch, (h, (w, ls)) in enumerate(zip(hs, baths)):
            assert 2 <= delta <= 3 and 1 <= h <= 5
            assert w[0] <= 15/(2*mp.sqrt(6*mp.pi))
            for j in range(N):
                assert ls[j]-masses[1] >= mp.mpf('2.5')+10*j+4*j*j
                if j:
                    assert w[j] <= 19/(4*mp.pi*j)
                if j in [0, 1, 8, 100, 255]:
                    assert abs(w[j]/rho(delta, h, j)-1) < mp.mpf('1e-65')
        for z in [mp.mpf(0), masses[0], masses[1]]:
            sigma = mp.matrix(2)
            sigma_upper = mp.matrix(2)
            branch_rows = []
            for branch, (w, ls) in enumerate(baths):
                val = mp.fsum(wj/(lj-z) for wj, lj in zip(w, ls))
                deriv = mp.fsum(wj/(lj-z)**2 for wj, lj in zip(w, ls))
                assert val+tail < mp.mpf(4)/3
                assert deriv+2*tail/5 < mp.mpf(8)/15
                gg = g[:, branch]*g[:, branch].T
                sigma += val*gg
                sigma_upper += (val+tail)*gg
                branch_rows.append({'branch': branch, 'beta_lower': str(val), 'beta_upper': str(val+tail), 'beta_prime_upper': str(deriv+2*tail/5)})
            scale = mp.diag([1/mp.sqrt(m) for m in masses])
            relative_upper = norm2(scale*sigma_upper*scale)
            assert relative_upper < epsilon
            rows.append({'C': str(C), 'z': str(z), 'terms_per_branch': N, 'branches': branch_rows, 'relative_response_upper': str(relative_upper), 'uniform_epsilon': str(epsilon)})
        if C == candidateC:
            oldepsilon = 3*gamma/(4*mp.pi)*(1+25/C+28/C**2)
            candidate = {'C': str(C), 'G_over_a2': str(gamma), 'Qall': str(q), 'mixed_epsilon': str(epsilon), 'mixed_light_absolute_shift_bound': str(epsilon*masses[0]), 'old_dipole_epsilon': str(oldepsilon), 'two_selected_channels_epsilon': str(oldepsilon+epsilon), 'two_selected_channels_light_absolute_shift_bound': str((oldepsilon+epsilon)*masses[0])}
    return {'enclosed_cases': rows, 'tail': str(tail), 'candidate_k2000': candidate}


def bisect(fn, left, right):
    fl, fr = fn(left), fn(right)
    assert fl*fr < 0, (left, right, fl, fr)
    for _ in range(210):
        mid = (left+right)/2
        fm = fn(mid)
        if fl*fm > 0:
            left, fl = mid, fm
        else:
            right = mid
    return (left+right)/2


def finite_mass_checks():
    rows = []
    for C in [1, 4, 100]:
        for gamma in ['.001', '.003']:
            for N in [16, 48]:
                masses, g, q, delta, hs, baths = physical(C, gamma, N)
                epsilon = 4*q/3
                assert epsilon < 1
                bath_masses = [v for w, ls in baths for v in ls]
                diag = np.array([float(x) for x in masses+bath_masses])
                A = np.diag(diag)
                for branch, (w, ls) in enumerate(baths):
                    block = np.outer(np.array(list(g[:, branch]), dtype=float).reshape(2), np.sqrt(np.array(w, dtype=float)))
                    start = 2+branch*N
                    A[:2, start:start+N] = block
                    A[start:start+N, :2] = block.T
                eig, vectors = np.linalg.eigh(A)
                assert eig[0] > 0
                def selfenergy(z, prime=False):
                    out = mp.matrix(2)
                    for branch, (w, ls) in enumerate(baths):
                        beta = mp.fsum(wj/(lj-z)**(2 if prime else 1) for wj, lj in zip(w, ls))
                        out += beta*g[:, branch]*g[:, branch].T
                    return out
                def determinant(z):
                    return mp.det(mp.diag([m-z for m in masses])-selfenergy(z))
                roots = [bisect(determinant, (1-epsilon)*m, m) for m in masses]
                error = max(abs(eig[p]-float(roots[p])) for p in range(2))
                assert error < 2e-9, (C, gamma, N, error)
                residue_errors = []
                for p in range(2):
                    v = vectors[:2, p]
                    prime = np.array(selfenergy(mp.mpf(float(eig[p])), True).tolist(), dtype=float)
                    residue = float(v@(np.eye(2)+prime)@v)
                    residue_errors.append(abs(residue-1))
                    assert abs(residue-1) < 2e-10
                scaled = A/np.sqrt(diag[:, None]*diag[None, :])
                assert np.linalg.eigvalsh(scaled)[0] >= float(1-mp.sqrt(epsilon))-2e-12
                rows.append({'C': C, 'G_over_a2': gamma, 'terms_per_branch': N, 'mass_root_discrepancy': error, 'residue_normalization_error': max(residue_errors), 'epsilon': str(epsilon)})
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    mp.mp.dps = 80
    print('Checking analytic gap/tails and independent distinct-field weights', flush=True)
    result = {'algebra': algebra(), 'representation': representation_checks()}
    print('Enclosing both complete mixed spectral branches', flush=True)
    result['spectral'] = spectral_checks()
    print('Comparing matrix poles and residues with the Schur response', flush=True)
    result['finite_mass'] = finite_mass_checks()
    result.update({'status': 'PASS', 'scope': 'specified nondifferentiated cubic representative and fixed global Friedrichs vacuum; local action, physical state and full stability remain open', 'versions': {'python': platform.python_version(), 'sympy': s.__version__, 'mpmath': mp.__version__, 'numpy': np.__version__}, 'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print('PASS: complete mixed spectral bound, tail, positive comparison and pole intervals', flush=True)


if __name__ == '__main__':
    main()
