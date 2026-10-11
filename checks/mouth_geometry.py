"""Derive the mouth radial equations from linearized tensors, independently.

Coordinates (t,r,u,phi), u=cos(theta), avoid trigonometric simplification.
Choose q=4 so the paper's rescaled gauge amplitudes equal the perturbation
amplitudes; kappa/g_c^2=e^2/2 then enforces its magnetic-charge convention.
The first-order quantum stress is zero in this mouth approximation, while
the prescribed Schwinger current is retained. This checks the reduction's
classical geometry, not the validity of the quantum truncation.
"""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import sympy as sp

spec = importlib.util.spec_from_file_location(
    "mouth_reduction", Path(__file__).with_name("mouth_reduction.py"))
reduction = importlib.util.module_from_spec(spec)
spec.loader.exec_module(reduction)
radial_equations = reduction.radial_equations


def verify():
    r, u, mass, e, ell, screen, p = sp.symbols("r u M0 e L s p", nonzero=True)
    funcs = [sp.Function(n)(r) for n in ["h0", "h1", "w0", "w1", "X"]]
    h0, h1, w0, w1, charge = funcs
    y = sp.Function("Y")(u)
    yp = sp.diff(y,u)
    y2 = (2*u*yp-ell*y)/(1-u*u)
    angular = {sp.diff(y,u,3): sp.diff(y2,u).subs(sp.diff(y,u,2),y2),
               sp.diff(y,u,2): y2}

    def simplify(expression):
        return sp.factor(expression.xreplace(angular).subs(sp.diff(y,u,2),y2))

    def derivative(expression, index, perturbation=False):
        if index == 0:
            return p*expression if perturbation else sp.S.Zero
        if index == 1:
            return sp.diff(expression,r)
        if index == 2:
            return sp.diff(expression,u)
        return sp.S.Zero

    f = 1-2*mass/r+e*e/r**2
    metric = sp.diag(-f, 1/f, r*r/(1-u*u), r*r*(1-u*u))
    inverse = sp.diag(-1/f, f, (1-u*u)/r**2, 1/(r*r*(1-u*u)))
    h = sp.zeros(4)
    h[0,3] = h[3,0] = -h0*(1-u*u)*yp
    h[1,3] = h[3,1] = -h1*(1-u*u)*yp
    dinverse = -inverse*h*inverse
    # The first-order volume variation vanishes, so the Maxwell density
    # below needs no delta(sqrt(-g)) term.
    assert sp.factor(sp.trace(inverse*h)) == 0
    gamma, dgamma = {}, {}
    for a in range(4):
        for b in range(4):
            for c in range(b,4):
                base = linear = sp.S.Zero
                for d in range(4):
                    dg = derivative(metric[d,c],b)+derivative(metric[d,b],c)-derivative(metric[b,c],d)
                    dh = derivative(h[d,c],b,True)+derivative(h[d,b],c,True)-derivative(h[b,c],d,True)
                    base += inverse[a,d]*dg/2
                    linear += (inverse[a,d]*dh+dinverse[a,d]*dg)/2
                gamma[a,b,c] = gamma[a,c,b] = sp.factor(base)
                dgamma[a,b,c] = dgamma[a,c,b] = simplify(linear)
    print("linear_connection_derived",flush=True)
    ricci, dricci = sp.zeros(4), sp.zeros(4)
    for a in range(4):
        for b in range(a,4):
            base = linear = sp.S.Zero
            for c in range(4):
                base += derivative(gamma[c,a,b],c)-derivative(gamma[c,a,c],b)
                linear += derivative(dgamma[c,a,b],c,True)-derivative(dgamma[c,a,c],b,True)
                for d in range(4):
                    base += gamma[c,d,c]*gamma[d,a,b]-gamma[c,d,b]*gamma[d,a,c]
                    linear += dgamma[c,d,c]*gamma[d,a,b]+gamma[c,d,c]*dgamma[d,a,b]
                    linear -= dgamma[c,d,b]*gamma[d,a,c]+gamma[c,d,b]*dgamma[d,a,c]
            ricci[a,b] = ricci[b,a] = sp.factor(base)
            dricci[a,b] = dricci[b,a] = simplify(linear)
    scalar = sp.factor(sp.trace(inverse*ricci))
    dscalar = simplify(sp.trace(inverse*dricci)+sp.trace(dinverse*ricci))
    assert scalar == dscalar == 0
    print("linear_curvature_derived",flush=True)
    field = sp.zeros(4)
    field[2,3], field[3,2] = sp.Integer(2), -sp.Integer(2)
    dfield = sp.zeros(4)
    electric = p*w1-sp.diff(w0,r)
    dfield[0,1] = electric*y
    dfield[0,2] = -w0*yp
    dfield[1,2] = -w1*yp
    dfield = dfield-dfield.T
    field_up = inverse*field*inverse
    background_maxwell = sp.Matrix([
        sp.factor(sum(derivative(r*r*field_up[a,b],a) for a in range(4)))
        for b in range(4)])
    assert background_maxwell == sp.zeros(4,1)
    dfield_up = dinverse*field*inverse+inverse*dfield*inverse+inverse*field*dinverse
    invariant = sum(field[a,b]*field_up[a,b] for a in range(4) for b in range(4))
    dinvariant = sum(dfield[a,b]*field_up[a,b]+field[a,b]*dfield_up[a,b]
                     for a in range(4) for b in range(4))
    stress = field*inverse*field.T-metric*invariant/4
    dstress = dfield*inverse*field.T+field*dinverse*field.T+field*inverse*dfield.T
    dstress -= h*invariant/4+metric*dinvariant/4
    background_residual = ricci-e*e/2*stress
    assert background_residual.applyfunc(sp.factor) == sp.zeros(4)
    einstein = (dricci-e*e/2*dstress).applyfunc(simplify)
    current_density = sp.Matrix([screen*sp.diff(charge,r)*y,-screen*p*charge*y,0,0])
    maxwell = sp.Matrix([simplify(sum(derivative(r*r*dfield_up[a,b],a,True)
                                    for a in range(4))-current_density[b])
                         for b in range(4)])
    assert simplify(derivative(current_density[0],0,True)+derivative(current_density[1],1,True)) == 0
    forbidden = [(a,b) for a in range(4) for b in range(a,4) if (a,b) not in [(0,3),(1,3),(2,3)]]
    assert all(einstein[a,b] == 0 for a,b in forbidden)
    assert maxwell[3] == 0
    paper = radial_equations(r,mass,e,ell,screen,p,h0,h1,w0,w1,electric,charge)
    matches = {}
    pairs = [
        ("Einstein_t_phi_A1",einstein[0,3],h0),
        ("Einstein_r_phi_A2",einstein[1,3],w1),
        ("Einstein_angular_A3",einstein[2,3],h0),
        ("Maxwell_t_A13",maxwell[0],w0),
        ("Maxwell_r_A14",maxwell[1],charge),
        ("Maxwell_angular_A6",maxwell[2],h0),
    ]
    for name,actual,probe in pairs:
        zero_derivatives = {sp.diff(z,r,n):sp.S.Zero for z in funcs for n in [1,2,3]}
        constants = {z:sp.Integer(z==probe) for z in funcs}
        ratio = simplify(actual.xreplace(zero_derivatives).xreplace(constants)/
                         paper[name].xreplace(zero_derivatives).xreplace(constants))
        assert ratio != 0
        residual = simplify(actual-ratio*paper[name])
        assert residual == 0,(name,residual)
        matches[name] = {"residual":"exact_zero","normalization_factor":str(ratio)}
        print(name,"exact_zero",flush=True)
    return {"radial_matches":matches,"background_Einstein_Maxwell":"exact_zero",
            "linear_volume_variation":"exact_zero",
            "scalar_curvature_and_linear_variation":"exact_zero",
            "forbidden_Einstein_components":len(forbidden),
            "forbidden_Maxwell_phi":"exact_zero","linear_charge_conservation":"exact_zero"}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    result={"status":"passed_independent_linear_geometry_derivation",
            "source":"https://arxiv.org/abs/2610.09847v1",
            "coordinates":"t,r,u=cos(theta),phi; time dependence exp(p t)",
            "normalization":"q=4; kappa/g_c^2=e^2/2; gauge amplitudes w",
            "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "comparison_module_sha256":hashlib.sha256(Path(__file__).with_name('mouth_reduction.py').read_bytes()).hexdigest(),
            "sympy_version":sp.__version__,"checks":verify(),
            "assumed_quantum_closure":"zero first-order quantum stress in the mouth; line-density Schwinger current",
            "quantum_truncation_validity_established":False,"lean_certified":False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")


if __name__=="__main__":
    main()
