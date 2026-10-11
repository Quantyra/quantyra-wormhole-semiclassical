"""Regulated Schur/Dirac-port checks; no physical continuum-limit claim.

Run from any directory. Writes a deterministic receipt only after all checks
pass. A failed assertion exits without replacing a previous receipt.
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
ALG_TOL = 2e-11
FD_TOL = 3e-6
ODE_TOL = 2e-9
BOUND_TOL = 3e-11
FD_STEPS = (0.002, 0.001, 0.0005)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.diag([1., -1.]).astype(complex)
I2 = np.eye(2, dtype=complex)


def err(a, b):
    return float(np.max(np.abs(a - b)))


def opnorm(a):
    return float(np.linalg.norm(a, 2))


def cdata(z):
    return [float(np.real(z)), float(np.imag(z))]


def ring(n, mass, theta, va, vb):
    """Wilson r=1 ring, four jets [value,a,b,ab] in Peierls phases."""
    jets = [np.zeros((2*n, 2*n), dtype=complex) for _ in range(4)]
    for i in range(n):
        s = slice(2*i, 2*i+2)
        # Positive Hermitian part guarantees a nonsingular finite regulator.
        jets[0][s, s] = (mass + 1 + 0.03*i/n)*I2 + 0.07j*SZ
        j = (i + 1) % n
        t = slice(2*j, 2*j+2)
        f = (-I2+SX)/2 * np.exp(1j*theta[i])
        b = (-I2-SX)/2 * np.exp(-1j*theta[i])
        for out, cf, cb in zip(jets, [1, 1j*va[i], 1j*vb[i], -va[i]*vb[i]],
                              [1, -1j*va[i], -1j*vb[i], -va[i]*vb[i]]):
            out[s, t] += cf*f
            out[t, s] += cb*b
    return jets


def hessian(j):
    r = np.linalg.inv(j[0])
    return np.trace(r @ j[2] @ r @ j[1] - r @ j[3])


def schur(jets, core):
    ext = sorted(set(range(jets[0].shape[0])) - set(core))
    a, b, c, e = ([j[np.ix_(rows, cols)] for j in jets]
                   for rows, cols in [(core, core), (core, ext),
                                      (ext, core), (ext, ext)])
    r0 = np.linalg.inv(e[0])
    ra = -r0 @ e[1] @ r0
    rb = -r0 @ e[2] @ r0
    rab = r0 @ e[1] @ r0 @ e[2] @ r0 + r0 @ e[2] @ r0 @ e[1] @ r0 - r0 @ e[3] @ r0
    s0 = a[0] - b[0] @ r0 @ c[0]
    sa = a[1] - b[1] @ r0 @ c[0] - b[0] @ ra @ c[0] - b[0] @ r0 @ c[1]
    sb = a[2] - b[2] @ r0 @ c[0] - b[0] @ rb @ c[0] - b[0] @ r0 @ c[2]
    sab = a[3] - (b[3] @ r0 @ c[0] + b[1] @ rb @ c[0] + b[1] @ r0 @ c[2]
                  + b[2] @ ra @ c[0] + b[0] @ rab @ c[0] + b[0] @ ra @ c[2]
                  + b[2] @ r0 @ c[1] + b[0] @ rb @ c[1] + b[0] @ r0 @ c[3])
    s = [s0, sa, sb, sab]
    return s, e


def action(d):
    sign, logabs = np.linalg.slogdet(d)
    return -logabs - np.log(sign)


def fd_mixed(fun, h):
    return (fun(h, h) - fun(h, -h) - fun(-h, h) + fun(-h, -h))/(4*h*h)


def regulated_checks():
    rows = []
    gauge_rows = []
    for n, mass in [(6, 0.35), (8, 0.7), (10, 1.1)]:
        theta = 0.4*np.arange(1, n+1)/(n*(n+1)/2)
        chi = np.sin(2*np.pi*np.arange(n)/n) + 0.31*np.cos(4*np.pi*np.arange(n)/n)
        gauge = chi - np.roll(chi, -1)
        sources = [np.ones(n)/n, np.eye(n)[1], np.eye(n)[n//2-1],
                   np.cos(0.8*np.arange(n)), gauge]
        partitions = [list(range(4)), list(range(2*(n//2))),
                      [k for i in range(0, n, 2) for k in (2*i, 2*i+1)]]
        for ia, ib in [(0, 0), (1, 2), (2, 3), (3, 3), (4, 4), (4, 0)]:
            va, vb = sources[ia], sources[ib]
            jets = ring(n, mass, theta, va, vb)
            full = hessian(jets)
            def value(t, s):
                return action(ring(n, mass, theta+t*va+s*vb, va, vb)[0])
            fd_errors = [float(abs(fd_mixed(value, h)-full)) for h in FD_STEPS]
            assert max(fd_errors) < FD_TOL, (n, ia, ib, fd_errors)
            for ic, core in enumerate(partitions):
                sj, ej = schur(jets, core)
                split = hessian(sj)+hessian(ej)
                determinant = np.linalg.det(jets[0])
                determinant_error = abs(determinant-np.linalg.det(sj[0])*np.linalg.det(ej[0]))/abs(determinant)
                assert abs(split-full) < ALG_TOL
                assert determinant_error < ALG_TOL
                rows.append(dict(n=n, mass=mass, sources=[ia, ib], partition=ic,
                                 hessian=cdata(full), split_error=float(abs(split-full)),
                                 determinant_relative_error=float(determinant_error),
                                 fd_errors=fd_errors))
                if ia == ib == 4:
                    si = np.linalg.inv(sj[0])
                    bubble = np.trace(si @ sj[2] @ si @ sj[1])
                    assert abs(full) < ALG_TOL and abs(hessian(sj)) < ALG_TOL and abs(hessian(ej)) < ALG_TOL
                    gauge_rows.append(dict(n=n, partition=ic, full=float(abs(full)),
                                           exterior=float(abs(hessian(ej))),
                                           schur=float(abs(hessian(sj))),
                                           omitted_quadratic_defect=float(abs(bubble))))
        holonomy = hessian(ring(n, mass, theta, sources[0], sources[0]))
        assert abs(holonomy) > 1e-5, (n, holonomy)
    assert max(x['omitted_quadratic_defect'] for x in gauge_rows) > 0.01

    # Bare-linear exterior sources: the contact term is entirely induced.
    n = 6
    zero = np.zeros(n)
    d = ring(n, 0.35, zero, zero, zero)[0]
    core = list(range(6))
    v = np.zeros_like(d)
    w = np.zeros_like(d)
    v[6:, 6:] = np.kron(np.diag([0.4, -0.2, 0.7]), SZ)
    w[6:, 6:] = np.kron(np.diag([-0.1, 0.6, 0.3]), SX)
    linear = []
    for a, b in [(v, v), (v, w), (w, w)]:
        jets = [d, a, b, np.zeros_like(d)]
        sj, ej = schur(jets, core)
        si = np.linalg.inv(sj[0])
        full = hessian(jets)
        split = hessian(sj) + hessian(ej)
        no_contact = np.trace(si @ sj[2] @ si @ sj[1]) + hessian(ej)
        fd_errors = [float(abs(fd_mixed(lambda t, s: action(d+t*a+s*b), h)-full)) for h in FD_STEPS]
        assert abs(full-split) < ALG_TOL and max(fd_errors) < FD_TOL
        linear.append(dict(full=cdata(full), split_error=float(abs(full-split)),
                           induced_vertex_norm=opnorm(sj[3]),
                           omitted_contact_defect=float(abs(full-no_contact)),
                           fd_errors=fd_errors))
    assert max(x['omitted_contact_defect'] for x in linear) > 1e-4

    # Same zero-source operator, different exterior source response.
    fixed_s0 = schur([d, v, v, np.zeros_like(d)], core)[0][0]
    null_s0 = schur([d, np.zeros_like(d), np.zeros_like(d), np.zeros_like(d)], core)[0][0]
    assert err(fixed_s0, null_s0) == 0
    response_difference = float(abs(hessian([d, v, v, np.zeros_like(d)])))
    assert response_difference > 1e-4

    # Even the complete gauged Schur kernel does not encode a spectator determinant.
    spectator = ring(6, 0.35, np.zeros(6), np.ones(6)/6, np.ones(6)/6)
    spectator_response = hessian(spectator)
    from scipy.linalg import block_diag
    augmented = [block_diag(x, y) for x, y in zip([d, v, v, np.zeros_like(d)], spectator)]
    old_s = schur([d, v, v, np.zeros_like(d)], core)[0]
    new_s = schur(augmented, core)[0]
    same_s_error = max(err(x, y) for x, y in zip(old_s, new_s))
    difference = hessian(augmented) - hessian([d, v, v, np.zeros_like(d)])
    assert same_s_error < ALG_TOL and abs(difference-spectator_response) < ALG_TOL
    assert abs(difference) > 1e-4
    return dict(cases=rows, pure_gauge=gauge_rows, bare_linear=linear,
                same_zero_source_different_response=response_difference,
                spectator=dict(schur_jet_error=same_s_error, response=cdata(difference)))


def jetmul(x, y):
    return np.array([x[0]*y[0], x[1]*y[0]+x[0]*y[1],
                     x[2]*y[0]+x[0]*y[2],
                     x[3]*y[0]+x[1]*y[2]+x[2]*y[1]+x[0]*y[3]])


def jetinv(x):
    return np.array([1/x[0], -x[1]/x[0]**2, -x[2]/x[0]**2,
                     2*x[1]*x[2]/x[0]**3-x[3]/x[0]**2])


def scatter_jets(j):
    a, b, c, d = (j[:, i, k] for i, k in [(0, 0), (0, 1), (1, 0), (1, 1)])
    di = jetinv(d)
    out = np.zeros((4, 2, 2), dtype=complex)
    out[:, 0, 0] = -jetmul(c, di)
    out[:, 0, 1] = di
    out[:, 1, 0] = a - jetmul(b, jetmul(c, di))
    out[:, 1, 1] = jetmul(b, di)
    return out


def scatter(m):
    a, b, c, d = m[0, 0], m[0, 1], m[1, 0], m[1, 1]
    return np.array([[-c/d, 1/d], [a-b*c/d, b/d]])


def transfer_jets(z, mass, segments, method='expm'):
    """Each segment (length,Va,Vb) has H=-i SZ d/dx+mass SX+t Va+s Vb."""
    state = np.zeros((8, 2), dtype=complex)
    state[:2] = I2
    for length, va, vb in segments:
        g = 1j*SZ @ (z*I2-mass*SX)
        ga = -1j*SZ @ va
        gb = -1j*SZ @ vb
        gen = np.zeros((8, 8), dtype=complex)
        for k in range(4):
            gen[2*k:2*k+2, 2*k:2*k+2] = g
        gen[2:4, :2] = ga
        gen[4:6, :2] = gb
        gen[6:8, 2:4] = gb
        gen[6:8, 4:6] = ga
        if method == 'expm':
            state = expm(length*gen) @ state
        else:
            sol = solve_ivp(lambda _, y: (gen @ y.reshape(8, 2)).ravel(), (0, length),
                            state.ravel(), method='DOP853', rtol=2e-12, atol=2e-13)
            assert sol.success
            state = sol.y[:, -1].reshape(8, 2)
    return state.reshape(4, 2, 2)


def direct_port(z, mass, segments, t, s):
    m = I2.copy()
    for length, va, vb in segments:
        m = expm(length*(1j*SZ @ (z*I2-mass*SX-t*va-s*vb))) @ m
    return scatter(m)


def continuum_checks():
    rows = []
    configs = [(0.12, 0.35, 0.2, 0.7), (0.2, -0.7, 0.8, 1.9),
               (0.65, 0.1, 0.4, 1.2), (0.9, 0.85, 0.6, 2.1),
               (0.05, 0.3, 0.5, 2.5), (1.3, -0.4, 0.3, 0.6)]
    for eta, omega, mass, length in configs:
        z = omega+1j*eta
        segments = [(0.4*length, 0.3*SX+0.2*SZ, 0.4*SY-0.1*I2),
                    (0.6*length, -0.2*SY+0.1*I2, 0.25*SX+0.2*SZ)]
        tj = transfer_jets(z, mass, segments)
        sj = scatter_jets(tj)
        oj = scatter_jets(transfer_jets(z, mass, segments, 'ode'))
        ode_error = err(sj, oj)
        assert ode_error < ODE_TOL
        va = max(opnorm(a) for _, a, _ in segments)
        vb = max(opnorm(b) for _, _, b in segments)
        ratios = [opnorm(sj[0]), opnorm(sj[1])/(va/(2*eta)),
                  opnorm(sj[2])/(vb/(2*eta)), opnorm(sj[3])/(va*vb/eta**2)]
        assert max(ratios) < 1+BOUND_TOL
        fd_errors = []
        for h in FD_STEPS:
            da = (direct_port(z, mass, segments, h, 0)-direct_port(z, mass, segments, -h, 0))/(2*h)
            db = (direct_port(z, mass, segments, 0, h)-direct_port(z, mass, segments, 0, -h))/(2*h)
            dab = fd_mixed(lambda t, s: direct_port(z, mass, segments, t, s), h)
            fd_errors.append(max(err(da, sj[1]), err(db, sj[2]), err(dab, sj[3])))
        assert max(fd_errors) < FD_TOL
        rows.append(dict(eta=eta, omega=omega, mass=mass, length=length,
                         norm_bound_ratios=ratios, ode_error=ode_error, fd_errors=fd_errors))

    gauge_rows = []
    for chi0, chim, chi1 in [(0., 0.7, 0.), (-0.3, 0.4, 0.8), (0.2, -0.5, 0.2)]:
        z, mass, length = 0.4+0.3j, 0.6, 1.7
        segments = [(0.4*length, -(chim-chi0)/(0.4*length)*SZ, -(chim-chi0)/(0.4*length)*SZ),
                    (0.6*length, -(chi1-chim)/(0.6*length)*SZ, -(chi1-chim)/(0.6*length)*SZ)]
        sj = scatter_jets(transfer_jets(z, mass, segments))
        q = np.diag([chi0, chi1])
        d1 = 1j*(q @ sj[0] - sj[0] @ q)
        d2 = -q @ q @ sj[0] + 2*q @ sj[0] @ q - sj[0] @ q @ q
        gauge_error = max(err(sj[1], d1), err(sj[3], d2))
        assert gauge_error < ALG_TOL
        gauge_rows.append(dict(phases=[chi0, chim, chi1], error=gauge_error,
                               first_norm=opnorm(sj[1]), second_norm=opnorm(sj[3])))

    chiral = []
    for eta in [0.02, 0.2, 1.]:
        for length in [0.2/eta, 1/eta, 2/eta, 5/eta]:
            # H=-i d/dx-a on the right-moving scalar interval.
            value = np.exp((0.4j-eta)*length)
            first = 1j*length*value
            second = -length**2*value
            ratios = [abs(first)*2*eta, abs(second)*eta**2]
            assert max(ratios) < 1+BOUND_TOL
            chiral.append(dict(eta=eta, length=length, first_bound_ratio=float(ratios[0]),
                               second_bound_ratio=float(ratios[1])))
    return dict(noncommuting_ports=rows, gauge_ports=gauge_rows, chiral_intervals=chiral)


def main():
    regulated = regulated_checks()
    continuum = continuum_checks()
    inputs = ['checks/gauged_mouth_interaction.py', 'notes/gauged-mouth-interaction.md',
              'evidence/gauged-mouth-source-review.json']
    report = dict(date='2026-10-11', status='passed',
                  scope='finite regulated exact identities and conditional continuum bounded-source ports; no physical matching or renormalized stress acceptance',
                  environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__),
                  tolerances=dict(algebra=ALG_TOL, finite_differences=FD_TOL, ode=ODE_TOL, norm_bounds=BOUND_TOL),
                  finite_difference_steps=FD_STEPS,
                  input_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
                  regulated=regulated, continuum=continuum,
                  lean_run=False, physical_continuum_limit_proved=False,
                  physical_damping_proved=False, novelty_accepted=False)
    path = ROOT/'evidence/gauged-mouth-interaction.json'
    path.write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status='passed', regulated_cases=len(regulated['cases']),
                          gauge_cases=len(regulated['pure_gauge']), bare_linear_cases=len(regulated['bare_linear']),
                          port_cases=len(continuum['noncommuting_ports']), gauge_ports=len(continuum['gauge_ports']),
                          chiral_cases=len(continuum['chiral_intervals']),
                          max_split_error=max(x['split_error'] for x in regulated['cases']),
                          max_regulated_fd_error=max(max(x['fd_errors']) for x in regulated['cases']),
                          max_port_fd_error=max(max(x['fd_errors']) for x in continuum['noncommuting_ports']),
                          max_port_ode_error=max(x['ode_error'] for x in continuum['noncommuting_ports'])), indent=2))


if __name__ == '__main__':
    main()
