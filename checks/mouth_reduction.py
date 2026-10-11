"""Reproduce the published radial-to-master reduction, including constraints.

This verifies algebra of the stated radial equations. A separate geometric
derivation is needed to verify the radial equations themselves.
"""

import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp


def radial_equations(r, mass, e, ell, screen, p, h0, h1, w0, w1, electric, charge):
    delta = r*r - 2*mass*r + e*e
    dr = lambda value, order=1: sp.diff(value, r, order)
    return {
        "Einstein_t_phi_A1": (ell*r*r - 4*mass*r + 6*e*e)*h0
            - r*r*delta*dr(h0, 2) + 2*p*r*delta*h1
            + p*r*r*delta*dr(h1) + 2*r*r*e*e*w0,
        "Einstein_r_phi_A2": 2*p*r**5*h0 - p*r**6*dr(h0)
            + (r**6*p*p + delta*((ell-2)*r*r+4*e*e))*h1
            + 2*r*r*e*e*delta*w1,
        "Einstein_angular_A3": p*r**5*h0 - 2*delta*(mass*r-e*e)*h1
            - r*delta**2*dr(h1),
        "Maxwell_t_A13": 2*ell*h0 + ell*r*r*w0 - 2*r*delta*dr(w0)
            - r*r*delta*dr(w0,2) + 2*p*r*delta*w1
            + p*r*r*delta*dr(w1) - screen*delta*dr(charge),
        "Maxwell_r_A14": -2*ell*delta*h1 + p*r**6*dr(w0)
            - (r**6*p*p+ell*r*r*delta)*w1 + screen*p*r**4*charge,
        "Maxwell_angular_A6": 2*p*r**5*h0
            + 4*delta*(r*r-3*mass*r+2*e*e)*h1
            - 2*r*delta**2*dr(h1) + p*r**7*w0
            - 2*r*r*delta*(mass*r-e*e)*w1 - r**3*delta**2*dr(w1),
        "anomaly_A15": -r**3*delta*electric - p*p*r**5*charge
            + 2*delta*(mass*r-e*e)*dr(charge) + r*delta**2*dr(charge,2),
        "electric_field_definition": electric - p*w1 + dr(w0),
    }


def verify():
    r, mass, e, ell, screen, p = sp.symbols("r M0 e L s p", nonzero=True)
    f = 1 - 2*mass/r + e*e/r**2
    g, d, xi = [sp.Function(name)(r) for name in ["g", "d", "xi"]]
    fields = sp.Matrix([g, d, xi])
    potential = f/r**3 * sp.Matrix([
        [ell*r - 6*mass + 4*e*e/r, -2*e, 0],
        [-2*e*(ell-2), ell*r + 4*e*e/r, ell*screen*r],
        [0, r, screen*r],
    ])
    second = (p*p*fields + potential*fields - f*sp.diff(f,r)*fields.diff(r))/f**2
    rule2 = {sp.diff(z,r,2): sp.factor(value) for z,value in zip(fields,second)}
    rule3 = {sp.diff(z,r,3): sp.diff(value,r).xreplace(rule2)
             for z,value in zip(fields,second)}
    h1 = r/f*g
    h0 = f/p*sp.diff(f*h1,r)
    electric = ell/(p*e*r*r)*(d+screen*xi)
    charge = ell/(p*e)*xi
    displacement = ell/(p*e)*d
    w1 = -2*h1/r**2 - p/(ell*f)*displacement
    w0 = -2*h0/r**2 - f/ell*sp.diff(displacement,r)
    residuals = radial_equations(r,mass,e,ell,screen,p,h0,h1,w0,w1,electric,charge)
    result = {}
    for name, expression in residuals.items():
        reduced = sp.factor(expression.xreplace(rule3).xreplace(rule2))
        assert reduced == 0, (name, reduced)
        result[name] = "exact_zero"
        print(name, "exact_zero", flush=True)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    data = {
        "status": "passed_published_radial_equation_reduction",
        "source": "https://arxiv.org/abs/2610.09847v1",
        "source_equations": ["3.5", "5.11", "5.20", "5.21", "A.1-A.15"],
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "sympy_version": sp.__version__,
        "restrictions": ["r != 0", "e != 0", "L != 0", "p != 0", "f != 0"],
        "residuals": verify(),
        "geometric_derivation_verified_by_this_script": False,
        "quantum_effective_action_validity_verified": False,
        "lean_certified": False,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")


if __name__ == "__main__":
    main()
