"""Actual two-sphere multipoles coupled to a spherical Maxwell throat.

Numerical truncations supplement the all-mode variational proof. No
Einstein equation, quantum state, or Dirac-return error is computed here.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_bvp
from scipy.linalg import solve_banded
from scipy.special import eval_legendre, roots_legendre
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]


def text_hash(path):
    return hashlib.sha256(path.read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest()


def angular(ell,c):
    # Y_l=(1-c^2) P_l', with its exact weighted norm.
    norm=np.sqrt(2*ell*(ell+1)/(2*ell+1))
    y=ell*(eval_legendre(ell-1,c)-c*eval_legendre(ell,c))/norm
    yp=-ell*(ell+1)*eval_legendre(ell,c)/norm
    return y,yp


def exterior(p,cutoff,order):
    c,weights=roots_legendre(order);weights=weights/(1-c*c)
    basis=np.array([angular(l,c)[0] for l in range(1,cutoff+1)])
    trace=np.zeros((2,order,2*cutoff));flux=trace.copy()
    v=np.zeros((2,order));f=v.copy()
    for side in range(2):
        r=p*np.sqrt(1-c*c);z=side-.5+p*c
        for source in range(2):
            zi=z-(source-.5);rho=np.hypot(r,zi);ci=zi/rho;si=r/rho
            for l in range(1,cutoff+1):
                y,yp=angular(l,ci);column=source*cutoff+l-1
                value=(p/rho)**l*y
                dr=(p/rho)**l/rho*(-l*y*si-yp*ci*si)
                dz=(p/rho)**l/rho*(-l*y*ci+yp*(1-ci*ci))
                trace[side,:,column]=value
                flux[side,:,column]=-(dr*np.sqrt(1-c*c)+dz*c)
        sign=1 if side==0 else -1
        dist=np.sqrt(1+p*p-2*sign*p*c)
        sr=(p-sign*c)/dist
        v[side]=p*p*(1-c*c)/(2*dist*(dist+1-sign*p*c))
        f[side]=-v[side]*(2/p-sr/dist-(sr-sign*c)/(dist+1-sign*p*c))
    project=lambda a:np.einsum('ln,n,snk->slk',basis,weights,a).reshape(2*cutoff,-1)
    T=project(trace);D=project(flux)
    LE=np.linalg.solve(T.T,D.T).T
    vp=np.einsum('ln,n,sn->sl',basis,weights,v).ravel()
    fp=np.einsum('ln,n,sn->sl',basis,weights,f).ravel()
    assert np.linalg.norm(LE-LE.T)<1e-7/p
    assert np.linalg.eigvalsh((LE+LE.T)/2)[0]*p>.43
    return dict(LE=LE,T=T,D=D,v=vp,f=fp,trace=trace,flux=flux,
                v_values=v,f_values=f,basis=basis,weights=weights)


def profile(t,length,core):
    distance=np.minimum(t,length-t)
    y=(distance-.25)/.25
    blend=np.zeros_like(y);blend[y>=1]=1
    mask=(y>0)&(y<1)
    left=np.exp(-1/y[mask]);right=np.exp(-1/(1-y[mask]))
    blend[mask]=left/(left+right)
    return (1-blend)/(1-np.minimum(distance,.5))+blend*core


def throat_fem(ell,length,core,level):
    edge=np.linspace(0,.5,32*level+1)
    mid=np.linspace(.5,length-.5,max(32,int(length*2))*level+1)
    grid=np.unique(np.r_[edge,mid,length-edge])
    h=np.diff(grid);nodes,weights=roots_legendre(8)
    local=(nodes+1)/2;shape=np.array([1-local,local])
    pts=grid[:-1,None]+h[:,None]*local
    potential=ell*(ell+1)*profile(pts,length,core)**2
    mass=np.einsum('en,n,in,jn,e->eij',potential,weights,shape,shape,h/2)
    diag=np.zeros(len(grid));off=-1/h+mass[:,0,1]
    diag[:-1]+=1/h+mass[:,0,0];diag[1:]+=1/h+mass[:,1,1]
    band=np.zeros((3,len(grid)-2));band[1]=diag[1:-1]
    band[0,1:]=off[1:-1];band[2,:-1]=off[1:-1]
    rhs=np.zeros((len(grid)-2,2));rhs[0,0]=-off[0];rhs[-1,1]=-off[-1]
    interior=solve_banded((1,1),band,rhs)
    C=np.diag([diag[0],diag[-1]])
    C[0]+=off[0]*interior[0];C[1]+=off[-1]*interior[-1]
    assert np.max(abs(C-C.T))<1e-9
    assert min(np.linalg.eigvalsh(C))>0
    return C


def throat_bvp(ell,length,core):
    def ode(t,y):
        mu=ell*(ell+1)*profile(t,length,core)**2
        return np.array([y[2],y[3],mu*y[0],mu*y[1]])
    def bc(a,b):return np.array([a[0]-1,a[1],b[0],b[1]-1])
    edge=np.linspace(0,.5,80)
    mesh=np.unique(np.r_[edge,np.linspace(.5,length-.5,320),length-edge])
    guess=np.array([1-mesh/length,mesh/length,-np.ones_like(mesh)/length,np.ones_like(mesh)/length])
    sol=solve_bvp(ode,bc,mesh,guess,tol=1e-10,max_nodes=18000)
    assert sol.success,(ell,length,core,sol.message)
    return np.array([[-sol.y[2,0],-sol.y[3,0]],[sol.y[2,-1],sol.y[3,-1]]])


def solve_case(p,cutoff,profile_index,throats,order):
    e=exterior(p,cutoff,order);LT=np.zeros_like(e['LE'])
    for l in range(1,cutoff+1):
        orientation=np.diag([1.,(-1.)**(l+1)])
        block=orientation@throats[profile_index,l]@orientation/p
        ids=[l-1,cutoff+l-1];LT[np.ix_(ids,ids)]=block
    rhs=e['LE']@e['v']-e['f']
    h=np.linalg.solve(e['LE']+LT,rhs);k=h-e['v'];coef=np.linalg.solve(e['T'],k)
    field_cost=float(k@e['LE']@k+h@LT@h)
    energy_shift=float(2*e['f']@k+field_cost)
    stationarity=float(np.linalg.norm(e['f']+e['LE']@k+LT@h)/(p+1e-300))
    boundary=e['v_values']+np.einsum('snk,k->sn',e['trace'],coef)
    expected=(e['basis'].T@h.reshape(2,cutoff).T).T
    jump=boundary-expected
    current=e['f_values']+np.einsum('snk,k->sn',e['flux'],coef)
    current+=(e['basis'].T@(LT@h).reshape(2,cutoff).T).T
    weighted=lambda a:float(np.sqrt(np.sum(a*a*e['weights'])))
    return dict(p=p,cutoff=cutoff,profile=profile_index,
        boundary_scaled=float(np.linalg.norm(h)/p**2),
        field_cost_scaled=field_cost/p**3,energy_shift_scaled=energy_shift/p**3,
        boundary_jump_scaled=weighted(jump)/p**2,current_jump_scaled=weighted(current)/p,
        stationarity_scaled=stationarity,first_four_scaled=(h.reshape(2,cutoff)[:,:4]/p**2).tolist())


def verify():
    c=sp.symbols('c',real=True)
    identities=[]
    for ell in range(1,7):
        Y=(1-c*c)*sp.diff(sp.legendre(ell,c),c)
        assert sp.expand(-(1-c*c)*sp.diff(Y,c,2)-ell*(ell+1)*Y)==0
        assert sp.integrate(sp.cancel(Y*Y/(1-c*c)),(c,-1,1))==sp.Rational(2*ell*(ell+1),2*ell+1)
        identities.append(ell)
    radius,height=sp.symbols('radius height',positive=True)
    distance=sp.sqrt(radius**2+height**2)
    for ell in range(1,5):
        potential=distance**(-ell)*(1-c*c)*sp.diff(sp.legendre(ell,c),c)
        potential=potential.subs(c,height/distance)
        laplace=sp.diff(potential,radius,2)-sp.diff(potential,radius)/radius+sp.diff(potential,height,2)
        assert sp.simplify(laplace)==0
    lapse,areal,theta=sp.symbols('lapse areal theta',positive=True)
    measure=lapse**2*areal**2*sp.sin(theta)
    assert sp.simplify(measure/lapse**2/(areal**2*sp.sin(theta)**2)-1/sp.sin(theta))==0
    assert sp.simplify(measure/areal**2/(areal**2*sp.sin(theta)**2)-lapse**2/(areal**2*sp.sin(theta)))==0
    x=sp.symbols('x',real=True);rr=sp.Function('r')(x);phi=sp.Function('phi')(x);mu=sp.symbols('mu')
    # Independent transformation of the source's electrostatic equation.
    z=rr*phi
    assert sp.simplify(rr*(sp.diff(z,x,2)-(sp.diff(rr,x,2)/rr+mu/rr**2)*z)
                       -(sp.diff(rr**2*sp.diff(phi,x),x)-mu*phi))==0
    profiles=[(2.,.2),(6.,1.),(20.,.05),(1000.,.005)]
    throat={};checks=[]
    for index,(length,core) in enumerate(profiles):
        for ell in range(1,9):
            coarse=throat_fem(ell,length,core,4)
            fine=throat_fem(ell,length,core,8)
            change=float(np.max(abs(coarse-fine)))
            assert change<.002,(index,ell,change)
            reference=throat_bvp(ell,length,core)
            independent=float(np.max(abs(reference-fine)))
            assert independent<.0003,(index,ell,independent)
            throat[index,ell]=reference
            checks.append(dict(profile=index,length_over_p=length,core_scaled_w=core,
                ell=ell,refinement=change,independent_bvp_error=independent,
                smallest_eigenvalue=float(min(np.linalg.eigvalsh(fine)))))
    cases=[]
    for p in [.1,.04,.015]:
        for index in range(len(profiles)):
            rows=[solve_case(p,n,index,throat,160) for n in [4,6,8]]
            change=max(abs(np.array(rows[-1]['first_four_scaled'])-np.array(rows[-2]['first_four_scaled'])).ravel())
            assert change<2e-8,(p,index,change)
            fine=rows[-1]
            independent=solve_case(p,8,index,throat,240)
            quad_change=abs(fine['field_cost_scaled']-independent['field_cost_scaled'])
            assert quad_change<1e-9
            assert fine['boundary_jump_scaled']<1e-6
            assert fine['current_jump_scaled']<2e-5
            assert fine['stationarity_scaled']<1e-12
            assert fine['field_cost_scaled']<10
            fine['multipole_refinement']=float(change);fine['angular_quadrature_refinement']=quad_change
            fine['truncation_current_jump_scaled']=[v['current_jump_scaled'] for v in rows]
            assert rows[-1]['current_jump_scaled']<rows[0]['current_jump_scaled']/100+1e-11
            cases.append(fine)
    scales={}
    for alpha in [Fraction(1),Fraction(3,2)]:
        scales[str(alpha)]={
            'field_energy_q_power':str(2-3*alpha),
            'physical_energy_q_power_with_g_squared_q_minus_two':str(4-3*alpha),
            'relative_to_q_over_ell_q_power':str(3-3*alpha),
            'relative_log_power':2,
            'physical_port_over_re_q_power':str(3-alpha)}
    assert Fraction(3)-3*Fraction(3,2)==-Fraction(3,2)
    inputs=['notes/coupled-maxwell-matching.md','notes/loaded-operator-realization.md','notes/mouth-mode-matching.md']
    return dict(status='passed_coupled_maxwell_matching_checks',exact_angular_modes=identities,
        exact_exterior_multipoles=4,exact_static_metric_weights=2,
        source_electrostatic_transform_checked=True,throat_maps_used='independent continuous BVP, cross-checked against finite elements',
        throat_dtn_cases=checks,coupled_exterior_cases=cases,scale_exponents=scales,
        source_free_maxwell_on_prescribed_geometry_proved=True,
        field_energy_error_order='q^2 p^3',physical_energy_error_order='q^2 p^3/(g^2 d)',
        selected_port='p=q^-3/2',metric_einstein_equations_solved=False,
        corrected_dirac_return_controlled=False,local_stress_backreaction_controlled=False,
        energy_length_derivatives_controlled=False,novelty_established=False,
        numerical_truncations_are_supplementary=True,lean_invocations=0,
        input_sha256={p:text_hash(ROOT/p) for p in inputs},
        versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True)
    args=parser.parse_args();result=verify();result['script_sha256']=text_hash(Path(__file__))
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],throat_cases=len(result['throat_dtn_cases']),
        exterior_cases=len(result['coupled_exterior_cases']))))
