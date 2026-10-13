"""S072: joint compact thermal transfer and prescribed matrix junction.

This does not identify a compact heat kernel with the full 4D state.
The finite spectral sums check an analytic theorem, not a large-q truncation.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import math
import numpy as np
import scipy
import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def allowed(q, charge, two_l):
    return (two_l - (q - 1) * charge) % 2 == 0


def fusion(two_l, two_lp, two_j):
    return (abs(two_l - two_lp) <= two_j <= two_l + two_lp
            and (two_l + two_lp + two_j) % 2 == 0)


def spectral(q, x, a, tail=90):
    two_j = q - 1
    charge_cut = math.ceil(math.sqrt(tail / x)) + 3
    spin_cut = math.ceil(2 * math.sqrt(tail / a)) + two_j + 3
    charges = np.arange(-charge_cut, charge_cut + 1)
    spins = np.arange(spin_cut + 1)
    cas = spins * (spins + 2) / 4
    dims = spins + 1
    energy = x * charges[:, None] ** 2 + a * cas[None, :]
    mask = allowed(q, charges[:, None], spins[None, :])
    z = np.sum(mask * dims[None, :] ** 2 * np.exp(-energy))
    pair_mask = ((np.abs(spins[:, None] - spins[None, :]) <= two_j)
                 & (spins[:, None] + spins[None, :] >= two_j)
                 & ((spins[:, None] + spins[None, :] + two_j) % 2 == 0))
    weights, gaps, sum_rule = [], [], 0.0
    for charge in charges:
        e0 = x * charge ** 2 + a * cas
        e1 = x * (charge + 1) ** 2 + a * cas
        active = (pair_mask & allowed(q, charge, spins)[:, None]
                  & allowed(q, charge + 1, spins)[None, :])
        multiplicity = dims[:, None] * dims[None, :] * active / (q * z)
        gap = e1[None, :] - e0[:, None]
        weight = multiplicity * np.exp(-(e0[:, None] + e1[None, :]) / 2)
        sum_rule += np.sum(multiplicity * gap
                           * (np.exp(-e0)[:, None] - np.exp(-e1)[None, :]))
        weights.append(weight[active])
        gaps.append(gap[active])
    weights, gaps = np.concatenate(weights), np.concatenate(gaps)
    lam = x + a * two_j * (two_j + 2) / 4
    c0 = np.sum(weights)
    moment2 = np.sum(weights * gaps ** 2)
    assert abs(sum_rule - 2 * lam) < 3e-12 * max(1, lam)
    assert 0 <= c0 <= 1 + 2e-13
    assert c0 >= 1 - lam / 4 - 2e-13
    assert moment2 <= 2 * lam + 3e-12
    times = []
    for tau in (0.0, .04, .08, .5, 2.0):
        corr = np.sum(weights * np.exp(1j * tau * gaps))
        lower = 1 - lam / 4 - lam * tau ** 2
        assert corr.real >= lower - 2e-12
        assert abs(corr.imag) < 2e-13
        times.append(dict(tau=tau, correlation=float(corr.real), lower=float(lower)))
    return dict(q=q, x=x, a=a, lambda_transfer=lam, C0=float(c0),
                second_moment=float(moment2), sum_rule_error=float(abs(sum_rule-2*lam)),
                charge_cut=charge_cut, twice_spin_cut=spin_cut, times=times), (charges, spins, energy, mask, z)


def heat_kernel_integral(q, data):
    # Independent character-space quadrature. These grids exactly integrate
    # the finite Fourier/character polynomials apart from floating error.
    charges, spins, energy, mask, z = data
    nphi = 4 * max(abs(charges)) + 7
    ntheta = 2 * max(spins) + q + 5
    phi = 2 * np.pi * np.arange(nphi) / nphi
    theta = np.pi * np.arange(1, ntheta + 1) / (ntheta + 1)
    chi = np.sin((spins[:, None] + 1) * theta) / np.sin(theta)
    basis = (spins[:, None] + 1) * chi
    kernel = np.exp(1j * phi[:, None] * charges) @ (mask * np.exp(-energy / 2)) @ basis
    haar = 2 * np.sin(theta) ** 2 / (ntheta + 1)
    norm = np.mean(np.abs(kernel) ** 2 @ haar)
    transfer = np.sin(q * theta) / (q * np.sin(theta))
    c0 = np.mean(np.exp(1j * phi) * (np.abs(kernel) ** 2 @ (haar * transfer))) / z
    assert abs(norm / z - 1) < 2e-12
    assert abs(c0.imag) < 2e-13
    return float(c0.real), float(abs(norm / z - 1))


def lattice_checks():
    count = 0
    for q in range(2, 13):
        for charge in range(-3, 4):
            for tl in range(15):
                for tlp in range(15):
                    if allowed(q, charge, tl) and fusion(tl, tlp, q-1):
                        assert allowed(q, charge+1, tlp)
                        count += 1
    # Fixed individual charge projection annihilates a charge-one insertion.
    shift = np.diag(np.ones(10), -1)
    for index in range(11):
        projector = np.zeros((11, 11)); projector[index, index] = 1
        assert np.linalg.norm(projector @ shift @ projector) == 0
    return dict(allowed_transfer_triples=count, fixed_charge_projection="zero")


def spin_matrices(q):
    j = mp.mpf(q-1)/2
    plus = mp.matrix(q)
    for i in range(q-1):
        m = -j+i
        plus[i+1, i] = mp.sqrt((j-m)*(j+m+1))
    return (plus+plus.T)/2, mp.diag([-j+i for i in range(q)])


def junction_checks():
    mp.mp.dps = 65
    rows = []
    for q in (2, 3, 5):
        jx, jz = spin_matrices(q)
        eye = mp.eye(q)
        zero = mp.matrix(q)
        def block(aa, bb, cc, dd):
            return mp.matrix([[([aa, bb, cc, dd][2*(i//q)+(j//q)])[i%q, j%q]
                              for j in range(2*q)] for i in range(2*q)])
        J = block(zero, -eye, eye, zero)
        def unitary(t):
            return mp.exp(1j*(mp.mpf('.2')*t+mp.mpf('.11')*t*t)) * mp.expm(1j*t*jx) * mp.expm(1j*mp.mpf('.3')*t*t*jz)
        def angle(t):
            return mp.mpf('.07')+mp.mpf('.13')*t+mp.mpf('.09')*t*t
        def rotation(t):
            u=unitary(t); c=mp.cos(angle(t)); ss=mp.sin(angle(t))
            return block(c*eye, -ss*u, ss*u.H, c*eye)
        def trace(mat):
            return sum(mat[i,i] for i in range(mat.rows))
        for tstr in ('.1', '.4', '.9'):
            t=mp.mpf(tstr)
            u=unitary(t)
            ud=mp.matrix([[mp.diff(lambda v:unitary(v)[i,j],t) for j in range(q)] for i in range(q)])
            cd=mp.re(trace(ud))/q
            coh=mp.re(trace(u))/q
            velocity=mp.re(trace(ud*ud.H))/q
            td=mp.diff(angle,t); th=angle(t)
            expected=(td*td+mp.sin(th)**2*velocity-td*coh-mp.sin(th)*mp.cos(th)*cd)/(4*mp.pi)
            approximations=[]
            for ds in ('1e-5','1e-7'):
                d=mp.mpf(ds)
                def traced(v):
                    gin=mp.eye(2*q)/(4*mp.pi*1j*mp.sin(v/2))+1j*J/(4*mp.pi*mp.cos(v/2))
                    return trace(rotation(t+v/2)*gin*rotation(t-v/2).H-gin)
                value=-1j*mp.diff(traced,d)/(2*q)
                approximations.append(value)
            err=abs(approximations[-1]-expected)
            assert err < mp.mpf('2e-13')
            assert err < abs(approximations[0]-expected)/5000
            rows.append(dict(q=q, t=float(t), flux_per_real=float(expected),
                             point_split_error=str(err)))
    return rows


def scale_checks():
    X = sp.symbols('X', positive=True)
    assert sp.expand((1+6*X)*(X+sp.Rational(5,6))-2*X*(1+3*X)) == 4*X+sp.Rational(5,6)
    assert sp.expand((1+6*X)-6*X) == 1
    dep=ROOT/'checks/angular_boundary_response.py'
    spec=importlib.util.spec_from_file_location('angular_boundary_response',dep)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    rows=[]
    for k in (2000,4000,8000):
        old=module.scale_family(k)
        lam=old['x_e_bare_upper']+old['chi_J_static_upper']
        for delta in (.25,.5):
            # T Delta_U=delta z0/(2 pi); pair evolution has tau=2 T U.
            tau=delta*.5/np.pi
            deficit=lam*(.25+tau*tau)
            assert deficit < 2.14e-5
            rows.append(dict(k=k,delta=delta,lambda_upper=lam,
                             maximum_pair_tau=tau,initial_deficit=lam/4,
                             free_evolution_deficit=deficit,
                             frozen_link_energy_upper=-3/32+11*lam/64))
    return rows,hashlib.sha256(dep.read_bytes()).hexdigest()


def main():
    lattice=lattice_checks()
    rows=[]
    for q in (2,3,4,5,8):
        for x,a in ((.12,.08),(.7,.35),(2.3,1.1)):
            row,data=spectral(q,x,a)
            integral,norm_error=heat_kernel_integral(q,data)
            row['character_integral_error']=abs(integral-row['C0'])
            row['character_norm_error']=norm_error
            assert row['character_integral_error']<3e-12
            rows.append(row)
    refined,_=spectral(4,.12,.08,tail=125)
    baseline=next(r for r in rows if r['q']==4 and r['x']==.12)
    refinement=abs(refined['C0']-baseline['C0'])
    assert refinement<2e-13
    junction=junction_checks()
    scales,dependency=scale_checks()
    receipt=dict(date='2026-10-13',scope='joint compact heat-kernel state and externally prescribed matrix junction; physical quantum matching remains conditional',
                 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                 dependencies={'checks/angular_boundary_response.py':dependency},
                 numpy=np.__version__,scipy=scipy.__version__,mpmath=mp.__version__,
                 uniform_scale_proof='positive polynomial differences: 1 and 4 X+5/6',
                 lattice=lattice,quantum_spectral_cases=rows,refinement_error=refinement,
                 prescribed_junction_cases=junction,single_species_comparisons=scales,
                 limitations=['The compact thermal Hamiltonian and operator factorization are not derived from the 4D junction',
                              'Free compact evolution omits the fermion-transfer backreaction',
                              'The matrix flux treats a prescribed classical group history; it is not a quantum stress substitution',
                              'No preparation, evolving-geometry or sufficient stability result'],passed=True)
    (ROOT/'evidence/quantum-junction-state.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: joint charge-spin lattice, 15 spectral/character comparisons, 9 matrix point-splits, 6 scale comparisons')


if __name__=='__main__':
    main()
