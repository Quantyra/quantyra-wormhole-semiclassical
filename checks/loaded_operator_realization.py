"""Supplementary checks for the loaded realization and collar smoothing.

Finite coupled-channel BVPs test flux, feedback and forcing algebra. They
do not prove infinite-dimensional existence or simulate the dipole PDE.
"""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_bvp, simpson
from scipy.linalg import expm
from scipy.special import roots_legendre
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]


def text_hash(path):
    return hashlib.sha256(path.read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest()


def real_matrix(a):
    return np.block([[a.real,-a.imag],[a.imag,a.real]])


def load_matrix(m,s,length):
    k=np.sqrt(s*s+m*m)
    denominator=k+s*np.tanh(k*length)
    reflection=m*np.tanh(k*length)/denominator
    transmission=k/(np.cosh(k*length)*denominator)
    return np.array([[-reflection,transmission],[transmission,reflection]])


def case(modes,q,s,forced=False):
    high=modes-1;dim=2*modes+2*high;le=.7;lt=.4
    J=np.diag(np.r_[np.ones(modes),-np.ones(modes)])
    angular=np.diag(np.arange(modes,dtype=float)*.3)
    angular+=.12/(1+abs(np.arange(modes)[:,None]-np.arange(modes)[None,:]))
    potential=np.block([[np.zeros_like(angular),1j*angular],[-1j*angular,np.zeros_like(angular)]])
    shift=.08*np.cos(np.arange(modes)[:,None]+np.arange(modes)[None,:])
    potential+=np.block([[shift,np.zeros_like(shift)],[np.zeros_like(shift),shift]])
    AE=-J@(s*np.eye(2*modes)+1j*potential)
    masses=2*np.sqrt(np.arange(1,modes)*(np.arange(1,modes)+q))
    AT=np.block([[-s*np.eye(high),np.diag(masses)],[np.diag(masses),s*np.eye(high)]])
    matrix=np.zeros((dim,dim),complex);matrix[:2*modes,:2*modes]=le*AE
    matrix[2*modes:,2*modes:]=lt*AT
    phi=np.exp(.37j*np.arange(1,modes))
    ground=np.zeros(2,complex) if forced else np.array([1.,.31+.27j])
    force=np.zeros(dim,complex)
    if forced:
        force[:2*modes]=np.exp(.3j*np.arange(2*modes))/np.sqrt(2*modes)
    forcing=force.copy();forcing[:2*modes]=le*J@force[:2*modes]
    A0=np.zeros((dim,dim),complex);A1=A0.copy();target=np.zeros(dim,complex)
    A0[0,0]=1;A1[1,modes]=1
    target[:2]=ground
    # E(right)_H=T(left)_H; E(left)_H=phase*T(right)_H.
    eh=np.r_[np.arange(1,modes),np.arange(modes+1,2*modes)]
    th=np.arange(2*modes,dim)
    for j,(e,t) in enumerate(zip(eh,th)):
        row=2+j;A1[row,e]=1;A0[row,t]=-1
        row=2+2*high+j;A0[row,e]=1;A1[row,t]=-phi[j%high]
    rm=real_matrix(matrix);r0=real_matrix(A0);r1=real_matrix(A1)
    rf=np.r_[forcing.real,forcing.imag];rt=np.r_[target.real,target.imag]
    def ode(x,y):return rm@y+rf[:,None]*np.sin(np.pi*x)
    def jac(x,y):return np.repeat(rm[:,:,None],len(x),axis=2)
    def bc(a,b):return r0@a+r1@b-rt
    solutions=[];node_counts=[]
    for size,tol in [(60,2e-7),(100,2e-9)]:
        mesh=np.linspace(0,1,size);guess=np.zeros((2*dim,size))
        sol=solve_bvp(ode,bc,mesh,guess,tol=tol,max_nodes=16000,
                      fun_jac=jac,bc_jac=lambda a,b:(r0,r1))
        assert sol.success,(modes,q,s,forced,sol.message,sol.x.size)
        solutions.append(sol);node_counts.append(int(sol.x.size))
    grid=np.linspace(0,1,2001)
    def complex_values(sol,x):
        values=sol.sol(x);return values[:dim]+1j*values[dim:]
    y=complex_values(solutions[-1],grid)
    refinement=float(np.max(abs(y-complex_values(solutions[0],grid))))
    assert refinement<2e-7,(modes,q,s,refinement)
    bulk=le*simpson(np.sum(abs(y[:2*modes])**2,axis=0),x=grid)
    bulk+=lt*simpson(np.sum(abs(y[2*modes:])**2,axis=0),x=grid)
    outgoing=np.array([y[modes,0],y[0,-1]])
    work=le*simpson(np.real(np.sum(y[:2*modes].conj()*force[:2*modes,None],axis=0))*np.sin(np.pi*grid),x=grid)
    balance=abs(np.vdot(outgoing,outgoing).real+2*s.real*bulk-np.vdot(ground,ground).real-2*work)
    assert balance<2e-7,(modes,q,s,balance)
    trace_bound=float(np.sqrt(bulk))
    feedback_error=None
    if forced:
        forcing_norm=np.sqrt(le/2)*np.linalg.norm(force[:2*modes])
        assert np.sqrt(bulk)<=forcing_norm/s.real+1e-9
    else:
        assert bulk<=np.vdot(ground,ground).real/(2*s.real)+1e-9
        transfer=expm(le*AE)
        aa,bb=transfer[:modes,:modes],transfer[:modes,modes:]
        cc,dd=transfer[modes:,:modes],transfer[modes:,modes:]
        di=np.linalg.inv(dd)
        S=np.block([[-di@cc,di],[aa-bb@di@cc,bb@di]])
        gidx=np.array([0,modes]);hidx=np.setdiff1d(np.arange(2*modes),gidx)
        feedback=np.zeros((2*high,2*high),complex)
        swap=np.array([[0.,1.],[1.,0.]])
        for j,mass in enumerate(masses):
            phase=np.diag([phi[j],1.])
            block=phase@swap@load_matrix(mass,s,lt)@swap@phase.conj().T
            ids=[j,j+high];feedback[np.ix_(ids,ids)]=block
        high_in=np.linalg.solve(np.eye(2*high)-feedback@S[np.ix_(hidx,hidx)],
                                feedback@S[np.ix_(hidx,gidx)]@ground)
        incoming=np.zeros(2*modes,complex);incoming[gidx]=ground;incoming[hidx]=high_in
        predicted=(S@incoming)[gidx]
        feedback_error=float(np.max(abs(predicted-outgoing)))
        assert feedback_error<3e-8,(modes,q,s,feedback_error)
    return dict(modes=modes,q=q,real_s=float(s.real),imag_s=float(s.imag),forced=forced,
        refinement=refinement,current_error=float(balance),bulk_norm=trace_bound,
        feedback_ground_error=feedback_error,node_counts=node_counts,
        max_high_trace=float(np.max(abs(y[np.r_[eh,th]][:,[0,-1]]))))


def verify():
    s,m,x=sp.symbols('s m x',positive=True)
    k=sp.sqrt(s*s+m*m);r=m*sp.tanh(k*x)/(k+s*sp.tanh(k*x))
    t=k/(sp.cosh(k*x)*(k+s*sp.tanh(k*x)))
    defect=2*s*sp.tanh(k*x)/(k+s*sp.tanh(k*x))
    assert sp.simplify((1-r*r-t*t-defect).rewrite(sp.exp))==0
    rho,theta,B0=sp.symbols('rho theta B0',positive=True)
    chi=sp.Function('chi')(rho)
    ap=B0*rho*rho*sp.sin(theta)**2*chi/2
    br=sp.diff(ap,theta)/(rho*rho*sp.sin(theta))
    bt=-sp.diff(ap,rho)/(rho*sp.sin(theta))
    curl=(sp.diff(rho*bt,rho)-sp.diff(br,theta))/rho
    assert sp.simplify(curl+B0*(rho*sp.diff(chi,rho,2)+4*sp.diff(chi,rho))*sp.sin(theta)/2)==0
    ca,cb=sp.symbols('ca cb')
    vacuum=ca+cb/rho**3
    assert sp.simplify(rho*sp.diff(vacuum,rho,2)+4*sp.diff(vacuum,rho))==0
    rows=[case(n,q,s) for n in [2,4,7] for q in [1,40]
          for s in [.05+0j,.1+.7j,.4+2j]]
    forcing=[case(n,40,.1+.7j,True) for n in [2,4,7]]
    limits=[]
    for mass in [10.,100.,1000.,10000.]:
        # Stable exact imaginary-frequency defect; transmission may underflow.
        eta=.1;k=np.hypot(eta,mass);h=np.tanh(k*.4)
        loss=2*eta*h/(k+eta*h)
        limits.append(dict(mass=mass,defect=float(loss),scaled=float(loss*mass/(2*eta))))
    assert abs(limits[-1]['scaled']-1)<2e-5
    assert all(a['defect']>b['defect'] for a,b in zip(limits,limits[1:]))
    collar=[]
    nodes,weights=roots_legendre(160)
    for q in [10**4,10**6,10**8]:
        for alpha in [.75,1.]:
            p=q**(-alpha)
            rho=p*(1.5+.5*nodes)
            # Bound the norm on any unit angular ground profile plus its
            # known lower corrector; there is no angular cutoff in this bound.
            upper=q*rho/(4*(1-rho)**2)
            hnorm=27/np.sqrt(160)*rho*rho
            integral=p*.5*np.sum(weights*upper**2*(1+hnorm**2))
            envelope=(q*p**1.5)**2
            assert integral<envelope
            for radius in [p,1.5*p,2*p]:
                theta=np.linspace(.001,np.pi-.001,257)
                c=np.cos(theta);dist=np.sqrt(1+radius*radius-2*radius*c)
                potential=q*radius*np.sin(theta)/(2*dist*(dist+1-radius*c))
                assert max(abs(potential))<=q*radius/(4*(1-radius)**2)*(1+1e-12)
            collar.append(dict(q=q,port_exponent=alpha,port=p,
                squared_norm_upper=float(integral),order_envelope=float(envelope),
                ratio=float(integral/envelope)))
    exact_exponents={}
    for a in [Fraction(3,4),Fraction(1,1)]:
        costs=[1-4*a,Fraction(3,2)-3*a,-4*a,Fraction(1,2)-3*a,1-5*a]
        exact_exponents[str(a)]=dict(smoothing=str(1-Fraction(3,2)*a),
            boundary_cost_terms=[str(v) for v in costs],
            boundary_loss_sqrt=str(max(costs)/2),
            bulk=str(-Fraction(3,8)),physical_port_over_re_q_power=str(3-a))
        assert max(costs)==Fraction(3,2)-3*a
    assert 1-Fraction(3,2)==-Fraction(1,2)
    inputs=['notes/loaded-operator-realization.md','notes/throat-load-feedback.md',
            'notes/joint-mouth-matching.md','notes/dipole-source-matching.md',
            'notes/logarithmic-long-return.md','notes/mouth-mode-matching.md']
    return dict(status='passed_loaded_realization_component_checks',coupled_network_cases=rows,
        forced_resolvent_cases=forcing,high_mode_defect_cases=limits,
        collar_smoothing_cases=collar,exact_port_exponents=exact_exponents,
        numerical_checks_are_supplementary=True,numerics_are_three_dimensional_dipole=False,
        ordinary_prescribed_loaded_operator_realized=True,
        smooth_prescribed_gauge_completion_controlled=True,
        exact_spherical_flat_vacuum_interpolation_obstructed=True,
        einstein_maxwell_interface_derived=False,full_quantum_state_and_stress_derived=False,
        novelty_established=False,finite_flux_accuracy_established=False,lean_invocations=0,
        input_sha256={p:text_hash(ROOT/p) for p in inputs},
        versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True)
    args=parser.parse_args();result=verify();result['script_sha256']=text_hash(Path(__file__))
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],networks=len(result['coupled_network_cases']),
        forcing=len(result['forced_resolvent_cases']),collars=len(result['collar_smoothing_cases']))))
