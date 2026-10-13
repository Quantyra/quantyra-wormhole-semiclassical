"""CAR Ward and KZ checks of S072 angular transfer matching.

Run with Python, numpy, scipy, sympy and mpmath. No Lean invocation.
The finite Fermi sea is used only for states through energy two, where
the retained +/-1/2,+/-3/2 modes are sufficient exactly.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import platform
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
import sympy as sp
from scipy.integrate import solve_ivp
from sympy.physics.wigner import clebsch_gordan, wigner_6j

ROOT = Path(__file__).resolve().parents[1]


# Reuse the established S053 CAR implementation without running its main.
_spec = importlib.util.spec_from_file_location("prior_angular_coset", ROOT/'checks/angular_current_coset.py')
_prior = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_prior)


def spin(q):
    return _prior.generators(q)[0]


def ward_case(q):
    t = spin(q)
    sj = q*(q*q-1)/12
    vac = {(1 << (2*q))-1: 1.0+0j}
    triple = _prior.current(vac, t[2], -1, q)
    triple = _prior.current(triple, t[1], 0, q)
    triple = _prior.current(triple, t[0], 1, q)
    value = _prior.inner(vac, triple)
    assert abs(value-1j*sj) < 1e-11
    return {"q":q, "three_current_imaginary_part":float(value.imag),
            "expected_S_J":sj, "three_current_error":float(abs(value-1j*sj))}


def recoupled_pair_matrices(q):
    j = sp.Rational(q-1, 2)
    rec = np.array([[float(sp.sqrt((2*l+1)*(2*r+1))*wigner_6j(j,j,l,j,j,r))
                     for r in range(q)] for l in range(q)])
    pair = np.diag([(l*(l+1)-2*float(j*(j+1)))/2 for l in range(q)])
    assert np.linalg.norm(rec@rec.T-np.eye(q)) < 1e-12
    return pair, rec@pair@rec.T


def direct_pair_check(q):
    """Independent tensor-product construction of the four-spin singlets."""
    j = sp.Rational(q-1, 2)
    ms = [j*-1+i for i in range(q)]
    basis = []
    for l in range(q):
        v = np.zeros((q,q,q,q), complex)
        for ia,a in enumerate(ms):
            for ib,b in enumerate(ms):
                for ic,c in enumerate(ms):
                    d = -a-b-c
                    if d not in ms or abs(a+b) > l:
                        continue
                    coeff = clebsch_gordan(j,j,l,a,b,a+b)*clebsch_gordan(j,j,l,c,d,c+d)
                    coeff *= (-1)**int(l-a-b)/sp.sqrt(2*l+1)
                    v[ia,ib,ic,ms.index(d)] = complex(coeff)
        basis.append(v)
    t = spin(q)
    direct = np.array([[sum(np.vdot(a, np.einsum('bj,ck,ajkd->abcd', g, g, b)) for g in t)
                        for b in basis] for a in basis])
    _, rec = recoupled_pair_matrices(q)
    # Clebsch-Gordan channel phases may differ from the Racah convention.
    parity = np.diag([(-1)**l for l in range(q)])
    error = min(np.linalg.norm(direct-rec), np.linalg.norm(direct-parity@rec@parity))
    assert error < 1e-12
    return {"q": q, "direct_tensor_recoupling_error": float(error)}


def kz_case(q):
    j = (q-1)/2
    level = q*(q*q-1)/6
    a, b = recoupled_pair_matrices(q)
    factor = 2/(level+2)
    def rhs(x, u):
        return (factor*(a/x+b/(x-1))@u.reshape(q,q)).reshape(-1)
    lo, hi = 0.2, 0.8
    sol = solve_ivp(rhs, (lo,hi), np.eye(q).reshape(-1), method='DOP853',
                    rtol=2e-12, atol=2e-14)
    assert sol.success
    u = sol.y[:,-1].reshape(q,q)
    refined = solve_ivp(rhs, (lo,hi), np.eye(q).reshape(-1), method='DOP853',
                        rtol=1e-13, atol=2e-15)
    assert refined.success
    refine_error = np.linalg.norm(u-refined.y[:,-1].reshape(q,q), ord=2)
    length = np.log(hi/lo)+np.log((1-lo)/(1-hi))
    bound = np.expm1(factor*j*(j+1)*length)
    error = np.linalg.norm(u-np.eye(q), ord=2)
    assert error < bound + 1e-10
    assert refine_error < 1e-10
    commutator = np.linalg.norm(a@b-b@a, ord=2)
    offdiag = np.linalg.norm(u-np.diag(np.diag(u)), ord=2)
    assert commutator > 0.1 and offdiag > 1e-3
    return {"q": q, "x_interval": [lo,hi], "log_pair_length": float(length),
            "transport_distance_from_identity": float(error),
            "analytic_bound": float(bound), "refinement_error": float(refine_error),
            "pair_commutator_norm": float(commutator), "offdiagonal_transport_norm": float(offdiag),
            "scope": "KZ transport on ambient invariant tensors; no physical sewing constants inferred"}


def exact_and_candidate():
    table = []
    for q in [1,2,3,4,5,6,10,100,2000**4]:
        level = Fraction(q*(q*q-1),6)
        angular_c = 3*level/(level+2)
        angular_h = Fraction(q*q-1,4)/(level+2)
        neutral_h = Fraction(1,2)-Fraction(1,2*q)-angular_h
        cc = q-1-angular_c
        assert cc >= 0 and neutral_h >= 0
        assert 1+angular_c < 4
        assert Fraction(1,q)+2*angular_h <= Fraction(4,q)
        if q in [1,2,3]: assert cc == 0 and neutral_h == 0
        if q == 4: assert cc == Fraction(1,2) and neutral_h == Fraction(1,16)
        table.append({"q": q,"level": str(level),"c_coset": str(cc),
                      "fermion_h_charge":str(Fraction(1,2*q)),
                      "fermion_h_angular":str(angular_h),"fermion_h_coset":str(neutral_h)})
    mp.mp.dps=80
    k=mp.mpf(2000);q=k**4;level=q*(q*q-1)/6
    b=(q*q-1)/(2*(level+2))
    deficit=1/q+b
    four=mp.expm1(b*mp.log(16))
    # Illustrative fixed-regulator norm counting, NOT the physical h match.
    r=q/k
    pairs=(2*r)*(2*r-1)/2
    exponent=b*pairs
    return {"exact_table":table,"candidate_k":2000,"candidate_q":int(q),
            "charge_plus_angular_transfer_weight":mp.nstr(deficit,30),
            "four_insertion_transport_bound_x_point2_to_point8":mp.nstr(four,30),
            "illustrative_r_equal_q_over_k":mp.nstr(r,30),
            "unit_pair_log_length_exponent_at_that_order":mp.nstr(exponent,30),
            "limitation":"The last order is a norm-counting diagnostic, not a distribution of physical events or a lower bound on error."}


def main():
    result={"date":"2026-10-13","status":"PASS",
            "scope":"Exact current algebra/stress split and separated-insertion KZ transport; not the interacting wormhole stress",
            "python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,"sympy":sp.__version__,
            "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "ward_cases":[ward_case(q) for q in [2,3,4,5,6,8]],
            "inherited_dependencies_sha256":{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in ["checks/angular_current_coset.py","checks/gluing_support_compatibility.py","evidence/angular-current-coset.json"]},
            "direct_recoupling":[direct_pair_check(q) for q in [2,3,4]],
            "kz_cases":[kz_case(q) for q in [3,4,5,8,12]],
            "exact_and_candidate":exact_and_candidate()}
    target=ROOT/'evidence/angular-transfer-matching.json'
    target.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: 6 CAR Ward cases, 3 direct recouplings, 5 KZ transports and exact candidate bounds.')
    print('Maximum Ward coefficient error:',max(c['three_current_error'] for c in result['ward_cases']))
    print('Maximum KZ refinement error:',max(c['refinement_error'] for c in result['kz_cases']))
    print('Candidate transfer weight:',result['exact_and_candidate']['charge_plus_angular_transfer_weight'])
    print('Candidate four-insertion transport bound:',result['exact_and_candidate']['four_insertion_transport_bound_x_point2_to_point8'])


if __name__=='__main__':
    main()
