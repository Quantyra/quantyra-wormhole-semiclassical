"""Reconstruct corrected throat fields and verify every radial constraint.

Connects the independently tensor-checked radial equations to the published
three-channel master system at orders alpha^0 and alpha^1. Requires p != 0.
"""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import sympy as sp


def module(name):
    spec=importlib.util.spec_from_file_location(name,Path(__file__).with_name(name+'.py'))
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result


geometry=module('throat_geometry')
energy=module('throat_energy')


def verify():
    r,ell,screen,p,alpha=sp.symbols('rho L s p alpha',real=True)
    f=1+r*r
    gamma=-4*r*r/3-r*(r*r+3)*sp.atan(r)+sp.log(f)
    phi=r*sp.atan(r)+sp.Rational(4,3)
    sf=f+alpha*gamma
    def linear(v):
        return sp.factor(v.subs(alpha,0))+alpha*sp.factor(sp.diff(v,alpha).subs(alpha,0))
    def solve(eq,target):
        coeff=sp.diff(eq,target)
        assert sp.diff(coeff,target)==0
        return linear(-eq.subs(target,0)/coeff)

    h0,h1,w0,w1,charge,electric=[sp.Function(n)(r) for n in ['h0','h1','w0','w1','X','E']]
    r0,r1=geometry.published_radial(r,ell,screen,p,h0,h1,w0,w1,charge)
    raw=[a+alpha*b for a,b in zip(r0,r1)]
    h0_from_h1=solve(raw[2],h0)
    # Use E=p*w1-w0' before solving the Maxwell constraints. This removes
    # derivatives of the gauge potentials from the algebraic reconstruction.
    electric_rules={sp.diff(w0,r,2):p*sp.diff(w1,r)-sp.diff(electric,r),
                    sp.diff(w0,r):p*w1-electric}
    w1_from_fields=solve(sp.expand(raw[4].xreplace(electric_rules)),w1)
    w0_from_fields=solve(sp.expand(raw[3].xreplace(electric_rules)),w0)
    print('triangular_constraints_solved',flush=True)

    g,d,xi=[sp.Function(n)(r) for n in ['g','d','xi']]
    fields=sp.Matrix([g,d,xi])
    n,p1,q1,_,_=energy.matrices(r,ell,screen,p*p)
    q=-p*p*sp.eye(3)-f*n+alpha*q1
    second=(-(sf*sp.diff(sf,r)*sp.eye(3)+alpha*p1*sf)*fields.diff(r)-q*fields)/sf**2
    second=second.applyfunc(linear)
    rule2={sp.diff(z,r,2):v for z,v in zip(fields,second)}
    rule3={sp.diff(z,r,3):linear(sp.diff(v,r).xreplace(rule2)) for z,v in zip(fields,second)}

    h1_recon=linear((1+alpha*phi/2)/sf*g)
    e_recon=ell/p*(1-alpha*phi)*(d+screen*xi)
    x_recon=ell/p*xi
    def replace_function(v,old,new):
        return v.xreplace({sp.diff(old,r,j):sp.diff(new,r,j) for j in range(4)})
    h0_recon=linear(replace_function(h0_from_h1,h1,h1_recon))
    def reconstruct(v):
        for old,new in [(h0,h0_recon),(h1,h1_recon),(electric,e_recon),(charge,x_recon)]:
            v=replace_function(v,old,new)
        return linear(v)
    w1_recon=reconstruct(w1_from_fields)
    w0_recon=reconstruct(w0_from_fields)
    checks={}
    def check(name,v):
        v=linear(linear(v).xreplace(rule3).xreplace(rule2))
        for order in [0,1]:
            residual=sp.factor(v.coeff(alpha,order))
            assert residual==0,(name,order,residual)
        checks[name]='both_orders_exact_zero'
        print(name,'both_orders_exact_zero',flush=True)
    rr0,rr1=geometry.published_radial(r,ell,screen,p,h0_recon,h1_recon,w0_recon,w1_recon,x_recon)
    names=['Einstein_tau_phi','Einstein_rho_phi','Einstein_angular','Maxwell_tau','Maxwell_rho','Maxwell_angular','charge_anomaly']
    for name,v0,v1 in zip(names,rr0,rr1):check(name,v0+alpha*v1)
    check('electric_definition',e_recon-p*w1_recon+sp.diff(w0_recon,r))
    check('gravitational_master_definition',sf*h1_recon/(1+alpha*phi/2)-g)
    check('displacement_master_definition',p/ell*((1+alpha*phi)*e_recon-screen*x_recon)-d)
    return checks


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    data={'status':'passed_corrected_throat_master_reconstruction',
          'source':'https://arxiv.org/abs/2610.09847v1',
          'source_equations':['3.5','7.3-7.5','B.12-B.14','C.1-C.8','D.1-D.8'],
          'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
          'dependency_sha256':{n:hashlib.sha256(Path(__file__).with_name(n+'.py').read_bytes()).hexdigest() for n in ['throat_geometry','throat_energy']},
          'sympy_version':sp.__version__,'restrictions':['p != 0','L != 0','formal expansion through alpha^1'],
          'residuals':verify(),'quantum_closure_derived_from_Dirac_theory':False,
          'global_matching_verified':False,'lean_certified':False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
