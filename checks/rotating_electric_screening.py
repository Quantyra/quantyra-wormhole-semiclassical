"""Independent finite-box PDE checks of the prescribed MMP screening equation.

No Lean, physical boundary closure, or infinite-domain numerical certificate.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import scipy
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import spsolve
import sympy as sp


def exact_checks():
    h, c, r, g, q, om, b, n = sp.symbols('h c r g q om b n', positive=True)
    x, rho, b0 = sp.symbols('x rho b0', positive=True)
    mono = c / (c + 2*h*h)
    barrier = sp.cosh(sp.sqrt(b0)*x/h) / sp.cosh(sp.sqrt(b0)*rho/h)
    checks = {
        'monopole_spherical_harmonic': sp.simplify((2*h*h+c)*mono-c),
        'barrier': sp.simplify(-h*h*sp.diff(barrier,x,2)+b*barrier-(b-b0)*barrier),
        'potential_normalization': sp.simplify((g*g*q*b/(2*sp.pi**2))*(-4*sp.pi**2*om*n/g**2)+2*om*q*b*n),
        'h_squared': sp.simplify((2*sp.pi**2/(g*g*q))*(g*g*q/(2*sp.pi**2))-1),
    }
    # Poisson identity for the unscreened potential, separately from the solver.
    xx, yy, zz = sp.symbols('xx yy zz', real=True)
    # Clear each positive monopole denominator before simplifying. A single
    # simplify() of two differently factored square roots can leave zero unevaluated.
    for sign in [-1,1]:
        radius2 = sp.expand(xx*xx+yy*yy+(zz-sp.Rational(sign,2))**2)
        phi = radius2**(-sp.Rational(1,2))
        a0 = -yy*phi/2
        by = -sp.diff(phi,yy)/2
        residual = sum(sp.diff(a0,v,2) for v in (xx,yy,zz))-2*by
        checks['unscreened_poisson_source_'+str(sign)] = sp.factor(sp.expand(residual*radius2**sp.Rational(5,2)))
    nu = sp.Rational(7,8)
    exponents = {
        'port_outside_extremal_radius': nu-sp.Rational(3,4),
        'electric_error': -nu/2,
        'electric_relative_to_joint_error': -nu/2+sp.Rational(3,8),
        'orbital_adiabaticity': 1-3*nu/2,
    }
    assert all(v == 0 for v in checks.values()), checks
    assert exponents['port_outside_extremal_radius'] > 0
    assert all(v < 0 for k,v in exponents.items() if k != 'port_outside_extremal_radius')
    return {k: str(v) for k,v in checks.items()}, {k:str(v) for k,v in exponents.items()}


def asymptotic_functions():
    r, z = sp.symbols('r z', positive=True)
    plus = r*r+(z-sp.Rational(1,2))**2
    minus = r*r+(z+sp.Rational(1,2))**2
    br = r*(plus**(-sp.Rational(3,2))-minus**(-sp.Rational(3,2)))/2
    bz = ((z-sp.Rational(1,2))*plus**(-sp.Rational(3,2))-(z+sp.Rational(1,2))*minus**(-sp.Rational(3,2)))/2
    b = sp.sqrt(br*br+bz*bz)
    n = br/b
    lap_n = sp.diff(n,r,2)+sp.diff(n,r)/r+sp.diff(n,z,2)-n/r**2
    return sp.lambdify((r,z),(n,lap_n/b),modules='numpy',cse=True)


POINTS = [(.25,.25),(.5,.25),(.75,.25),(.5,.75)]
HVALUES = [.2,.1,.05,.025]
MESHES = [32,64,128]


def solve_box(radius, cells_per_unit, h, functions):
    step=1/cells_per_unit
    nr=int(round(radius*cells_per_unit))
    r=np.arange(1,nr)*step
    z=np.arange(-nr+1,nr)*step
    rr,zz=np.meshgrid(r,z,indexing='ij')
    plus=rr*rr+(zz-.5)**2
    minus=rr*rr+(zz+.5)**2
    br=rr*(plus**(-1.5)-minus**(-1.5))/2
    bz=((zz-.5)*plus**(-1.5)-(zz+.5)*minus**(-1.5))/2
    b=np.hypot(br,bz)
    # L = partial_r^2 + (1/r) partial_r + partial_z^2 - 1/r^2.
    lr=diags([np.full(len(r)-1,1/step**2)-1/(2*step*r[1:]),
              -2/step**2-1/r**2,
              np.full(len(r)-1,1/step**2)+1/(2*step*r[:-1])],[-1,0,1],format='csr')
    lz=diags([np.full(len(z)-1,1/step**2),np.full(len(z),-2/step**2),
              np.full(len(z)-1,1/step**2)],[-1,0,1],format='csr')
    lap=kron(lr,eye(len(z)),format='csr')+kron(eye(len(r)),lz,format='csr')
    op=-h*h*lap+diags(b.ravel())
    rhs=br.ravel()
    start=time.monotonic()
    u=spsolve(op.tocsc(),rhs).reshape(b.shape)
    duration=time.monotonic()-start
    residual=float(np.max(np.abs(op@u.ravel()-rhs))/max(1,np.max(np.abs(rhs))))
    bound=float(np.max(np.abs(u)))
    assert residual < 1e-9, residual
    assert bound <= 1+1e-10, bound
    samples=[]
    for rp,zp in POINTS:
        i=int(round(rp/step))-1
        j=int(round(zp/step))+nr-1
        assert abs(r[i]-rp)<1e-12 and abs(z[j]-zp)<1e-12
        n,correction=map(float,functions(rp,zp))
        val=float(u[i,j])
        samples.append(dict(r=rp,z=zp,value=val,algebraic=n,second_order=n+h*h*correction,
                            error=abs(val-n),second_error=abs(val-n-h*h*correction)))
    return dict(radius=radius,cells_per_unit=cells_per_unit,h=h,unknowns=int(u.size),
                solve_seconds=duration,linear_relative_residual=residual,max_abs=bound,samples=samples)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    exact,exponents=exact_checks()
    functions=asymptotic_functions()
    cases=[]
    for h in HVALUES:
        for mesh in MESHES:
            cases.append(solve_box(2,mesh,h,functions))
        cases.append(solve_box(3,64,h,functions))
    refinement=[]
    for h in HVALUES:
        row=[next(c for c in cases if c['h']==h and c['radius']==2 and c['cells_per_unit']==m) for m in MESHES]
        outer=next(c for c in cases if c['h']==h and c['radius']==3)
        e0=max(abs(a['value']-b['value']) for a,b in zip(row[0]['samples'],row[1]['samples']))
        e1=max(abs(a['value']-b['value']) for a,b in zip(row[1]['samples'],row[2]['samples']))
        assert e1 < .7*e0+1e-8, (h,e0,e1)
        refinement.append(dict(h=h,coarse_change=e0,fine_change=e1,
            outer_box_change_same_mesh=max(abs(a['value']-b['value']) for a,b in zip(row[1]['samples'],outer['samples'])),
            max_error_over_h2=max(s['error'] for s in row[2]['samples'])/h**2,
            max_second_error_over_h4=max(s['second_error'] for s in row[2]['samples'])/h**4))
    # Eventual asymptotic improvement, not a fit or an infinite-domain claim.
    finest=[next(c for c in cases if c['h']==h and c['radius']==2 and c['cells_per_unit']==128) for h in HVALUES]
    leading=[max(s['error'] for s in c['samples']) for c in finest]
    second=[max(s['second_error'] for s in c['samples']) for c in finest]
    assert all(b < a for a,b in zip(leading,leading[1:])), leading
    assert second[-1] < leading[-1], (leading,second)
    result=dict(status='passed_prescribed_screening_checks',
        versions={'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__},
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        exact_identities=exact,exact_exponents=exponents,cases=cases,refinement=refinement,
        global_weak_solution_proved_for_prescribed_equation=True,
        finite_box_numerics_certify_infinite_domain=False,physical_electric_state_closed=False,
        rotating_dirac_return_proved=False,complete_energy_derived=False,
        novelty_established=False,lean_invocations=0)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ['status','exact_identities','exact_exponents','refinement']}))


if __name__=='__main__':
    main()
