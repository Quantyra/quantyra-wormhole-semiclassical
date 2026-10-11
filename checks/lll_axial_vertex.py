"""Derive the axial LLL Dirac vertex from a tetrad and its spin connection.

Warped spherical metric, arbitrary t,z,theta dependent axial shift. Derivatives
of the conformal factor and radius are retained before projection. Mode integrals are
checked separately in lll_projection.py.
"""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def verify():
    t, z, th, ph, eps = sp.symbols('t z theta phi epsilon', real=True)
    coords = [t, z, th, ph]
    vt, vz = [sp.Function(n)(t, z, th) for n in ['v_t', 'v_z']]
    N, R = [sp.Function(n)(t, z) for n in ['N', 'R']]
    co = sp.diag(N, N, R, R*sp.sin(th))
    co[3, 0], co[3, 1] = eps*R*sp.sin(th)*vt, eps*R*sp.sin(th)*vz
    flat = sp.diag(-1, 1, 1, 1)
    metric = co.T*flat*co
    frame = co.inv()
    inverse = frame*flat*frame.T
    g0, gi0 = metric.subs(eps, 0), inverse.subs(eps, 0)
    dg, dgi = metric.diff(eps).subs(eps, 0), inverse.diff(eps).subs(eps, 0)
    connection, variation = {}, {}
    for a in range(4):
        for b in range(4):
            for c in range(4):
                value = delta = 0
                for d in range(4):
                    jet = sp.diff(g0[d, c], coords[b])+sp.diff(g0[b, d], coords[c])-sp.diff(g0[b, c], coords[d])
                    djet = sp.diff(dg[d, c], coords[b])+sp.diff(dg[b, d], coords[c])-sp.diff(dg[b, c], coords[d])
                    value += gi0[a, d]*jet/2
                    delta += (dgi[a, d]*jet+gi0[a, d]*djet)/2
                connection[a, b, c] = sp.simplify(value)
                variation[a, b, c] = sp.simplify(delta)
    sx = sp.Matrix([[0, 1], [1, 0]])
    sy = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    sz = sp.diag(1, -1)
    gamma = [sp.I*sp.kronecker_product(sx, sp.eye(2)),
             sp.kronecker_product(sy, sp.eye(2)),
             sp.kronecker_product(sz, sx), sp.kronecker_product(sz, sy)]
    for a in range(4):
        for b in range(4):
            assert gamma[a]*gamma[b]+gamma[b]*gamma[a] == 2*flat[a, b]*sp.eye(4)
    co0, dco = co.subs(eps, 0), co.diff(eps).subs(eps, 0)
    fr0, dfr = frame.subs(eps, 0), frame.diff(eps).subs(eps, 0)
    spin, dspin = [], []
    for mu in range(4):
        omega, domega = sp.zeros(4), sp.zeros(4)
        for a in range(4):
            for b in range(4):
                val = delta = 0
                for nu in range(4):
                    jet = sp.diff(fr0[nu, b], coords[mu])+sum(connection[nu, mu, la]*fr0[la, b] for la in range(4))
                    djet = sp.diff(dfr[nu, b], coords[mu])+sum(variation[nu, mu, la]*fr0[la, b]+connection[nu, mu, la]*dfr[la, b] for la in range(4))
                    val += co0[a, nu]*jet
                    delta += dco[a, nu]*jet+co0[a, nu]*djet
                omega += flat[a, a]*val*gamma[a]*gamma[b]/4
                domega += flat[a, a]*delta*gamma[a]*gamma[b]/4
        spin.append(omega.applyfunc(sp.simplify))
        dspin.append(domega.applyfunc(sp.simplify))
    q, m = sp.symbols('q m', real=True)
    perturbation = sp.zeros(4)
    for mu in range(4):
        gg = sum((gamma[a]*fr0[mu, a] for a in range(4)), sp.zeros(4))
        dgg = sum((gamma[a]*dfr[mu, a] for a in range(4)), sp.zeros(4))
        perturbation += gg*dspin[mu]+dgg*spin[mu]
        if mu == 3:
            perturbation += dgg*sp.I*(m-q*sp.cos(th)/2)
    # Full warped-background vertex, including the chirality-changing curl.
    angular_shifts = [v*(m-q*sp.cos(th)/2)*sp.eye(4)
                      +sp.sin(th)*sp.diff(v, th)*sp.kronecker_product(sp.eye(2), sz)/4
                      for v in [vt, vz]]
    algebraic = -sp.I*(gamma[0]*angular_shifts[0]+gamma[1]*angular_shifts[1])/N
    curl = R*sp.sin(th)*(sp.diff(vz, t)-sp.diff(vt, z))/(4*N**2)*gamma[0]*gamma[1]*gamma[3]
    assert (perturbation-algebraic-curl).applyfunc(sp.simplify) == sp.zeros(4)
    # The static spatial source has no longitudinal curl. Check both
    # transverse chiralities, as needed for the massive-sector bound.
    vv = sp.Function('v')(z, th)
    static = perturbation.subs({vt: 0, vz: vv, N: 1, R: 1}).doit().applyfunc(sp.simplify)
    angular = (vv*(m-q*sp.cos(th)/2)*sp.eye(4)
               +sp.sin(th)*sp.diff(vv, th)*sp.kronecker_product(sp.eye(2), sz)/4)
    assert (static+sp.I*gamma[1]*angular).applyfunc(sp.simplify) == sp.zeros(4)
    # q>0 zero modes have negative transverse chirality (indices 1,3).
    reduced = perturbation.extract([1, 3], [1, 3]).applyfunc(sp.simplify)
    shifts = [v*(m-q*sp.cos(th)/2)-sp.sin(th)*sp.diff(v, th)/4 for v in [vt, vz]]
    expected = -sp.I*(sp.I*sx*shifts[0]+sy*shifts[1])/N
    residual = (reduced-expected).applyfunc(sp.simplify)
    assert residual == sp.zeros(2)
    return {'Clifford_relations': 16, 'projected_Dirac_matrix_entries': 4,
            'static_spatial_full_Dirac_matrix_entries': 16,
            'full_warped_Dirac_vertex_matrix_entries': 16,
            'curl_vertex': 'R sin(theta) (partial_t v_z-partial_z v_t) gamma0 gamma1 gamma3/(4 N^2)',
            'all_residuals': 'exact_zero',
            'effective_gauge_shift_before_mode_integral': 'v_a (m-q cos(theta)/2) - sin(theta) partial_theta v_a /4',
            'axial_shift_convention': 'v_a=-h_a P_l\'(cos(theta))/R^2',
            't_and_z_derivatives_retained_before_projection': True,
            'background': 'N(t,z)^2*(-dt^2+dz^2)+R(t,z)^2*dOmega^2; A_phi=q cos(theta)/2',
            'not_computed': 'Higher Landau levels, second-order contact terms, or quantum vacuum stress'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    data = {'status': 'passed_independent_axial_Dirac_vertex',
            'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'sympy_version': sp.__version__, 'checks': verify(), 'lean_certified': False}
    args.output.write_text(json.dumps(data, indent=2)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
