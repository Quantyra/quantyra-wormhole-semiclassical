"""Check the conditional mouth-operator energy identity; no Lean invocation."""

import argparse
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import quad


def symbolic_checks():
    r, mass, e, ell, screen = sp.symbols("r M0 e L s", positive=True)
    k = ell - 2
    f = 1 - 2 * mass / r + e**2 / r**2
    prefactor = f / r**3
    raw = sp.Matrix([
        [ell * r - 6 * mass + 4 * e**2 / r, -2 * e, 0],
        [-2 * e * k, ell * r + 4 * e**2 / r, ell * screen * r],
        [0, r, screen * r],
    ])
    weight = sp.diag(k, 1, ell * screen)
    deform = sp.diag(f / r, 0, 0)
    potential = prefactor * raw
    adjusted = potential + f * deform.diff(r) - deform**2
    gram_rows = sp.Matrix([[k * r, -2 * e, 0], [0, 1, screen]])
    gram = prefactor * gram_rows.T * sp.diag(1 / r, ell * r) * gram_rows
    identities = {
        "weighted_symmetry": weight * potential - potential.T * weight,
        "weighted_gram_factorization": weight * adjusted - gram,
        "source_determinant": sp.Matrix([raw.det() -
            8 * screen * e**2 / r * (r**2 - 3 * mass * r + 2 * e**2)]),
    }
    result = {}
    for name, expression in identities.items():
        reduced = expression.applyfunc(sp.factor)
        assert reduced == sp.zeros(*reduced.shape), (name, reduced)
        result[name] = "exact_zero"
    return result


def quadrature_case(mass, screen, multipole, inner_nonzero):
    # These functions are evaluated directly, without symbolic lambdification.
    e, a, b = 1.0, 2.2 if mass > 1 else 1.2, 8.0
    ell = multipole * (multipole + 1)
    k = ell - 2
    weights = np.array([k, 1.0, ell * screen])

    def fields(r):
        t = (r - a) / (b - a)
        envelope = (1 - t)**2 * (1 if inner_nonzero else t)
        envelope_prime = (-2 * (1 - t) if inner_nonzero
                          else (1 - t) * (1 - 3 * t)) / (b - a)
        h = np.array([1 + t, np.sin(2 * t), 0.3 - t * t])
        hp = np.array([1, 2 * np.cos(2 * t), -2 * t]) / (b - a)
        return envelope * h, envelope_prime * h + envelope * hp

    def values(r):
        f = 1 - 2 * mass / r + e * e / (r * r)
        assert f > 0
        z, zr = fields(r)
        raw = np.array([
            [ell*r - 6*mass + 4*e*e/r, -2*e, 0],
            [-2*e*k, ell*r + 4*e*e/r, ell*screen*r],
            [0, r, screen*r],
        ])
        original = f * np.dot(weights * zr, zr) + np.dot(weights*z, raw @ z) / r**3
        dz = zr.copy()
        dz[0] += z[0] / r
        squares = f * np.dot(weights * dz, dz)
        squares += (k*r*z[0] - 2*e*z[1])**2 / r**4
        squares += ell * (z[1] + screen*z[2])**2 / r**2
        return original, squares

    original, error_a = quad(lambda r: values(r)[0], a, b, epsabs=1e-10, epsrel=1e-12)
    squares, error_b = quad(lambda r: values(r)[1], a, b, epsabs=1e-10, epsrel=1e-12)
    boundary = k * (1 - 2*mass/a + e*e/a**2) / a * fields(a)[0][0]**2
    residual = abs(original - squares - boundary)
    tolerance = 2e-10 * max(1.0, abs(original))
    assert residual < tolerance, (residual, tolerance)
    assert original >= boundary - tolerance
    return {
        "M0": mass, "screening_squared": screen, "l": multipole,
        "r_interval": [a, b], "nonzero_inner_value": inner_nonzero,
        "original_energy": original, "squares_energy": squares,
        "inner_boundary_term": boundary, "absolute_residual": residual,
        "quadrature_error_estimates": [error_a, error_b],
    }


def negative_potential_example():
    r, mass, e, ell, screen = 1.4, 1.0, 1.0, 6.0, 100.0
    k = ell - 2
    matrix = np.array([
        [ell*r - 6*mass + 4*e*e/r, -2*e*np.sqrt(k), 0],
        [-2*e*np.sqrt(k), ell*r + 4*e*e/r, np.sqrt(ell*screen)*r],
        [0, np.sqrt(ell*screen)*r, screen*r],
    ])
    eigenvalues = np.linalg.eigvalsh(matrix) * (1 - 2*mass/r + e*e/r**2) / r**3
    assert eigenvalues[0] < 0
    return {"r": r, "M0": mass, "e": e, "L": ell,
            "screening_squared": screen, "eigenvalues": eigenvalues.tolist()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = {
        "status": "passed_conditional_operator_checks_not_physical_certification",
        "source": "https://arxiv.org/html/2610.09847v1",
        "source_equations": ["5.20", "5.21", "5.22"],
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "python": platform.python_version(),
        "versions": {p: importlib.metadata.version(p) for p in
                     ["numpy", "scipy", "sympy", "mpmath"]},
        "symbolic": symbolic_checks(),
        "quadrature": [quadrature_case(mass, screen, multipole, nonzero)
                       for mass in [0.8, 1.0, 1.2]
                       for screen in [0.01, 3.0, 100.0]
                       for multipole in [2, 5]
                       for nonzero in [False, True]],
        "negative_pointwise_potential": negative_potential_example(),
        "lean_certified": False,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": output["status"], "symbolic_identities": len(output["symbolic"]),
                      "quadrature_cases": len(output["quadrature"]),
                      "max_relative_residual": max(c["absolute_residual"] /
                           max(1, abs(c["original_energy"])) for c in output["quadrature"])}))


if __name__ == "__main__":
    main()
