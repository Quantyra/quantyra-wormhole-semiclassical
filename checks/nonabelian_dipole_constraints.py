"""Independent classical checks of the angular projection, not quantum capture."""
from pathlib import Path
import argparse
import hashlib
import json
import math

import numpy as np
from numpy.polynomial.legendre import leggauss
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
RNG = np.random.default_rng(20261013)
CROSS = np.cross


def norm2(x):
    return float(np.vdot(x, x).real)


def error(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b)))) / max(
        1.0, float(np.max(np.abs(a))), float(np.max(np.abs(b))))


def sphere(nu=12, nphi=24):
    u, w = leggauss(nu)
    for z, weight in zip(u, w):
        st = math.sqrt(1-z*z)
        for phi in 2*math.pi*np.arange(nphi)/nphi:
            n = np.array([st*math.cos(phi), st*math.sin(phi), z])
            et = np.array([z*math.cos(phi), z*math.sin(phi), -st])
            ep = np.array([-math.sin(phi), math.cos(phi), 0.0])
            yield n, et, ep, st, weight*2*math.pi/nphi


def curvatures(a, h, da, dh, p):
    # da[i,j]=partial_i A_j; indices 0=t,1=x.
    fa = da[0, 1]-da[1, 0]-CROSS(a[0], a[1])
    d_h = dh[0, 1]-dh[1, 0]-CROSS(a[0], h[1])+CROSS(a[1], h[0])
    return fa, p*fa+d_h


def coordinate_maxwell(a, h, da, dh, p, f0, r, f, coupling):
    """Construct metric and d[-P cos(theta)dphi+n.(PA+H)+a0] directly."""
    total = 0.0
    c = p*a+h
    dc = p*(da[0, 1]-da[1, 0])+dh[0, 1]-dh[1, 0]
    det_error = 0.0
    for n, et, ep, st, weight in sphere():
        k = np.array([[np.dot(CROSS(aa, n), et),
                       np.dot(CROSS(aa, n), ep)/st] for aa in a]).T
        frame = np.eye(4)
        frame[2:, :2] = -k
        metric = frame.T @ np.diag([-f, f, r*r, r*r*st*st]) @ frame
        inv = np.linalg.inv(metric)
        det_error = max(det_error, error(np.linalg.det(metric), -f*f*r**4*st*st))
        field = np.zeros((4, 4))
        field[0, 1] = f0+np.dot(n, dc)
        for mu in range(2):
            field[mu, 2] = -np.dot(et, c[mu])
            field[mu, 3] = -st*np.dot(ep, c[mu])
        field[2, 3] = p*st
        field = field-field.T
        raised = inv @ field @ inv.T
        contraction = np.sum(field*raised)
        total += weight*(-f*r*r/(4*coupling**2))*contraction
    return total, det_error


def rodrigues(axis, theta):
    x, y, z = axis
    k = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    return np.eye(3)+math.sin(theta)*k+(1-math.cos(theta))*(k@k)


def rotated_jets(a, h, da, dh):
    axis = RNG.normal(size=3)
    axis /= np.linalg.norm(axis)
    rot = rodrigues(axis, RNG.normal())
    grad = RNG.normal(size=2)
    hessian = RNG.normal(size=(2, 2))
    hessian = (hessian+hessian.T)/2
    ap = np.array([rot@a[j]+grad[j]*axis for j in range(2)])
    hp = np.array([rot@h[j] for j in range(2)])
    dap, dhp = np.empty((2, 2, 3)), np.empty((2, 2, 3))
    for i in range(2):
        for j in range(2):
            dap[i, j] = rot@da[i, j]+grad[i]*CROSS(axis, rot@a[j])+hessian[i, j]*axis
            dhp[i, j] = rot@dh[i, j]+grad[i]*CROSS(axis, rot@h[j])
    return ap, hp, dap, dhp, rot


def symbolic_moment_map():
    a = sp.Matrix(sp.symbols('a0:3'))
    h = sp.Matrix(sp.symbols('h0:3'))
    e = sp.Matrix(sp.symbols('e0:3'))
    v = sp.Matrix(sp.symbols('v0:3'))
    charge = a.cross(e)+h.cross(v)
    coords, momenta = list(a)+list(h), list(e)+list(v)
    def bracket(x, y):
        return sp.expand(sum(sp.diff(x, q)*sp.diff(y, p)-sp.diff(x, p)*sp.diff(y, q)
                             for q, p in zip(coords, momenta)))
    for i, j, k in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        assert sp.expand(bracket(charge[i], charge[j])-charge[k]) == 0
    eps = sp.Matrix(sp.symbols('eps0:3'))
    generator = eps.dot(charge)
    for field in [a, h, e, v]:
        variation = sp.Matrix([bracket(component, generator) for component in field])
        assert (variation-eps.cross(field)).applyfunc(sp.expand) == sp.zeros(3, 1)
    return {'algebra_relations': 3, 'canonical_variations': 12, 'exact': True}


def spin_matrices(q):
    j = (q-1)/2
    m = np.arange(q)-j
    plus = np.zeros((q, q), dtype=complex)
    for k in range(q-1):
        plus[k+1, k] = math.sqrt((j-m[k])*(j+m[k]+1))
    return [(plus+plus.T)/2, (plus-plus.T)/(2j), np.diag(m)], m


def lll_projection(q):
    t, m = spin_matrices(q)
    j = (q-1)/2
    projected = np.zeros((3, q, q), dtype=complex)
    overlap = np.zeros((q, q), dtype=complex)
    for n, et, ep, st, weight in sphere(32, 64):
        phi = math.atan2(n[1], n[0])
        amps = np.array([math.sqrt(q*math.comb(q-1, k)/(4*math.pi))*
                         ((1-n[2])/2)**((j-mm)/2)*
                         ((1+n[2])/2)**((j+mm)/2)*np.exp(1j*mm*phi)
                         for k, mm in enumerate(m)])
        outer = np.outer(amps.conj(), amps)
        overlap += weight*outer
        projected += weight*n[:, None, None]*outer
    err = max(error(projected, 2*np.array(t)/(q+1)), error(overlap, np.eye(q)))
    comm = max(error(t[i]@t[jj]-t[jj]@t[i], 1j*t[k])
               for i, jj, k in [(0, 1, 2), (1, 2, 0), (2, 0, 1)])
    assert max(err, comm) < 3e-13
    return {'q': q, 'projection_error': err, 'spin_algebra_error': comm}


def run():
    rows = []
    for case in range(16):
        a, h = RNG.normal(size=(2, 2, 3))
        da, dh = RNG.normal(size=(2, 2, 2, 3))
        p, r, f, coupling, gn = RNG.uniform(.4, 2.0, size=5)
        f0 = RNG.normal()
        kg, ke, mu = r**4/(6*gn), 4*math.pi*r*r/(3*coupling**2), 8*math.pi/(3*coupling**2)
        fa, g = curvatures(a, h, da, dh, p)
        direct, determinant_error = coordinate_maxwell(a, h, da, dh, p, f0, r, f, coupling)
        reduced = (4*math.pi*r*r/coupling**2)*f0*f0/(2*f)+ke*norm2(g)/(2*f)
        reduced += mu*(norm2(h[0])-norm2(h[1]))/2 - 2*math.pi*f*p*p/(coupling**2*r*r)
        action_error = error(direct, reduced)
        ap, hp, dap, dhp, rot = rotated_jets(a, h, da, dh)
        fp, gp = curvatures(ap, hp, dap, dhp, p)
        covariance_error = max(error(fp, rot@fa), error(gp, rot@g))
        transformed, _ = coordinate_maxwell(ap, hp, dap, dhp, p, f0, r, f, coupling)
        covariance_error = max(covariance_error, error(transformed, direct))
        # General local jet of n, including base derivatives, tests the form identity.
        n = RNG.normal(size=3); n /= np.linalg.norm(n)
        dn = RNG.normal(size=(2, 3)); dn -= (dn@n)[:, None]*n
        dn_cov = dn-np.array([CROSS(aa, n) for aa in a])
        omega_cov = np.dot(n, CROSS(dn_cov[0], dn_cov[1]))+np.dot(n, fa)
        omega_local = np.dot(n, CROSS(dn[0], dn[1]))+np.dot(dn[0], a[1])-np.dot(dn[1], a[0])
        omega_local += np.dot(n, da[0, 1]-da[1, 0])
        form_error = error(omega_cov, omega_local)
        v = p*a+h; dv = p*da+dh
        fv, _ = curvatures(v/p, np.zeros_like(h), dv/p, np.zeros_like(dh), p)
        false_connection_error = error(p*fv, g-CROSS(h[0], h[1])/p)
        # Direct Legendre transform versus the integrated form plus its divergence.
        eh, ea = ke*g/f, kg*fa/f+p*ke*g/f
        deh, dea, rho = RNG.normal(size=(3, 3))
        alpha = .37
        lb = kg*norm2(fa)/(2*f)+ke*norm2(g)/(2*f)+mu*(norm2(h[0])-norm2(h[1]))/2
        h_local = np.dot(ea, da[0, 1])+np.dot(eh, dh[0, 1])-lb+np.dot(rho, a[0]+alpha*h[0])
        d_ea, d_eh = dea-CROSS(a[1], ea), deh-CROSS(a[1], eh)
        constraint = -d_ea+CROSS(h[1], eh)+rho
        b = -d_eh+alpha*rho
        h_integrated = f*norm2(ea-p*eh)/(2*kg)+f*norm2(eh)/(2*ke)+mu*norm2(h[1])/2
        h_integrated += np.dot(a[0], constraint)+np.dot(h[0], b)-mu*norm2(h[0])/2
        divergence = np.dot(da[1, 0], ea)+np.dot(a[0], dea)+np.dot(dh[1, 0], eh)+np.dot(h[0], deh)
        legendre_error = error(h_local, h_integrated+divergence)
        square_error = error(np.dot(h[0], b)-mu*norm2(h[0])/2,
                             norm2(b)/(2*mu)-mu*norm2(h[0]-b/mu)/2)
        dipole_spin = np.linalg.norm(CROSS(h[1], eh))
        spin_bound = r/math.sqrt(2*f)*(mu*norm2(h[1])/2+f*norm2(eh)/(2*ke))
        assert dipole_spin <= spin_bound*(1+1e-13)
        errors = [action_error, determinant_error, covariance_error, form_error,
                  false_connection_error, legendre_error, square_error]
        assert max(errors) < 2e-11, (case, errors)
        rows.append(dict(case=case, action_error=action_error, determinant_error=determinant_error,
                         covariance_error=covariance_error, form_error=form_error,
                         connection_identity_error=false_connection_error,
                         legendre_error=legendre_error, auxiliary_square_error=square_error))
    # Exact Cartan change of variables, not a repetition of the old radial solver.
    p, alpha, dw, dv, w, v = sp.symbols('P alpha dw dv w v')
    assert sp.expand(p*dw+(dv-p*dw)-dv) == 0
    q = sp.symbols('q', positive=True)
    assert sp.simplify((w+2*(v-q*w/2)/(q+1))-(w+2*v)/(q+1)) == 0
    quadrupoles = []
    for case in range(5):
        u = RNG.normal(size=3)
        mean, second = 0.0, 0.0
        for n, et, ep, st, weight in sphere():
            quad = np.dot(n, u)**2-norm2(u)/3
            mean += weight*quad/(4*math.pi)
            second += weight*quad*quad/(4*math.pi)
        ratio = math.sqrt(second)/(norm2(u)/3)
        assert abs(mean)/max(1, norm2(u)) < 2e-14
        assert error(second, 4*norm2(u)**2/45) < 2e-13
        assert abs(ratio-2/math.sqrt(5)) < 2e-13
        quadrupoles.append({'case': case, 'mean': mean, 'second': second, 'ratio': ratio})
    return {
        'status': 'PASS', 'scope': 'classical projected action and constraints; no quantum-state or stress error bound',
        'seed': 20261013, 'floating_tolerance': 2e-11,
        'coordinate_sphere_quadrature': {'gauss_nodes': 12, 'azimuth_nodes': 24},
        'noncommuting_cases': rows, 'moment_map': symbolic_moment_map(),
        'lll_projection': [lll_projection(q) for q in [2, 3, 4, 7, 10]],
        'cartan_symbolic_identities': 2, 'quadrupole_cases': quadrupoles,
        'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'limitations': ['No numerical four-dimensional Einstein evolution.',
                        'No regulator, quantum Hamiltonian construction or interacting correlator.',
                        'No bound relative to capture/area margins or wormhole creation claim.']}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT/'evidence/nonabelian-dipole-constraints.json')
    args = parser.parse_args()
    result = run()
    args.output.write_bytes((json.dumps(result, indent=2)+'\n').encode('utf-8'))
    print(json.dumps({'status': result['status'], 'noncommuting_cases': len(result['noncommuting_cases']),
                      'lll_projections': len(result['lll_projection']), 'quadrupole_cases': len(result['quadrupole_cases']),
                      'maximum_action_error': max(x['action_error'] for x in result['noncommuting_cases']),
                      'maximum_covariance_error': max(x['covariance_error'] for x in result['noncommuting_cases']),
                      'output': str(args.output)}, indent=2))
