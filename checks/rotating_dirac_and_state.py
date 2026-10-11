"""Verify local rotating Weyl congruence and the free-state UV diagnostic."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import expit
import sympy as sp

PAULI=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]],complex)
ID=np.eye(2,dtype=complex)

def sig(v):return np.einsum('j,jab->ab',v,PAULI)

def exact_checks():
    # Exact matrices, independent of the floating numerical representation.
    pauli=[sp.Matrix([[0,1],[1,0]]),sp.Matrix([[0,-sp.I],[sp.I,0]]),sp.diag(1,-1)]
    ident=sp.eye(2)
    t=sp.symbols('t',positive=True)
    v=2*t/(1+t*t);gam=(1+t*t)/(1-t*t)
    s=(ident+t*pauli[2])/sp.sqrt(1-t*t)
    checks={}
    def add(name,m):
        vals=[sp.simplify(e) for e in m]
        assert all(e==0 for e in vals),(name,vals)
        checks[name]=True
    add('boost_density',s*s-gam*(ident+v*pauli[2]))
    checks['unit_determinant']=sp.simplify(s.det()-1)==0
    for j in range(3):
        expected=pauli[j] if j<2 else pauli[2]/gam
        add('boost_principal_'+str(j),s*(pauli[j]-(v if j==2 else 0)*ident)*s-expected)
    bx,by,bz=sp.symbols('bx by bz',real=True)
    b=[bx,by,bz];sb=sum((b[j]*pauli[j] for j in range(3)),sp.zeros(2))
    for j in range(3):add('linear_cancellation_'+str(j),(pauli[j]*sb+sb*pauli[j])/2-b[j]*ident)
    om=sp.symbols('om',real=True)
    # b=(-om*y, om*x, 0); independently differentiate each boost coefficient.
    db=sp.Matrix([[0,-om,0],[om,0,0],[0,0,0]])
    derivative=-sp.I*sum((pauli[j]*pauli[k]*db[k,j]/2 for j in range(3) for k in range(3)),sp.zeros(2))
    add('linear_spin_connection',derivative-om*pauli[2])
    at,ax,ay,az=sp.symbols('at ax ay az',real=True)
    checks['coordinate_gauge_conversion']=sp.expand(sum(b[j]*a for j,a in enumerate([ax,ay,az]))-(at+bx*ax+by*ay+bz*az)+at)==0
    mu,p,lam,speed=sp.symbols('mu p lam speed',positive=True)
    phase=2*sp.pi/(2*sp.pi)**3*sp.integrate(p**3,(p,0,lam))*sp.integrate(1,(mu,1/speed,1))
    checks['uv_cap_coefficient']=sp.simplify(phase-lam**4*(1-1/speed)/(16*sp.pi**2))==0
    nu=sp.Rational(15,16)
    exponents={'static_bulk':-sp.Rational(3,8),'covariant_gradient':sp.Rational(7,8),
        'uncancelled_rotation':sp.Rational(7,8)-nu/2,
        'boosted_rotation_bulk':sp.Rational(7,8)-nu,
        'electric':-nu/2,'port_to_extremal_radius':nu-sp.Rational(3,4),
        'orbital_adiabaticity':1-3*nu/2,
        'transit_to_buffer':(1-3*nu/2)/2,
        'buffer_to_light_cylinder':(1-3*nu/2)/2,
        'separation_to_buffer':nu-(sp.Rational(1,2)+3*nu/4)}
    assert all(checks.values())
    assert exponents['uncancelled_rotation']>0 and exponents['boosted_rotation_bulk']<0
    assert all(exponents[k]<0 for k in ['transit_to_buffer','buffer_to_light_cylinder','separation_to_buffer'])
    return checks,{k:str(v) for k,v in exponents.items()}


def boost(x,omega):
    b=np.cross(omega,x)
    gam=1/np.sqrt(1-np.dot(b,b))
    aa=np.sqrt((gam+1)/2);cc=gam/(2*aa)
    s=aa*ID+cc*sig(b)
    tt=np.eye(3)-gam/(gam+1)*np.outer(b,b)
    ds=[];dt=[]
    for j in range(3):
        db=np.cross(omega,np.eye(3)[j])
        dg=gam**3*np.dot(b,db)
        da=dg/(4*aa);dc=dg/(2*aa)-gam*da/(2*aa**2)
        ds.append(da*ID+dc*sig(b)+cc*sig(db))
        dt.append(-dg/(gam+1)**2*np.outer(b,b)-gam/(gam+1)*(np.outer(db,b)+np.outer(b,db)))
    return b,gam,s,tt,np.array(ds),np.array(dt)


def phi(x,k):
    amplitude=np.exp(-.3*np.dot(x,x)+1j*np.dot(k,x))
    val=np.array([1+.2*x[0]+.1j*x[1],.3-.15*x[2]+.2j*x[0]])
    dval=np.array([[.2,.2j],[.1j,0],[0,-.15]],complex)
    grad=amplitude*(dval+(-.6*x+1j*k)[:,None]*val)
    return amplitude*val,grad


def direct_case(x,omega,k,q,xi,normal):
    b,gam,s,tt,ds,dt=boost(x,omega)
    val,grad=phi(x,k)
    field=q*np.array([.3,-.2,.7])
    gauge=np.cross(field,x)/2
    at=.13*np.linalg.norm(omega)*np.cos(x[0])
    z=-sig(omega)/2-at*ID
    mm=PAULI-b[:,None,None]*ID
    c=-1j*sum(s@mm[j]@ds[j] for j in range(3))+s@z@s
    piv=-1j*grad-gauge[:,None]*val
    reduced=xi*s@s@val+1j*(sum(sig(tt[:,j])@piv[j] for j in range(3))+c@val)
    def psi(xx):return boost(xx,omega)[2]@phi(xx,k)[0]
    raw_errors=[]
    for step in [1e-4,5e-5]:
        gradpsi=np.array([(psi(x+step*np.eye(3)[j])-psi(x-step*np.eye(3)[j]))/(2*step) for j in range(3)])
        ps=s@val
        original=xi*ps+1j*(sum(mm[j]@(-1j*gradpsi[j]-gauge[j]*ps) for j in range(3))+z@ps)
        original=s@original
        raw_errors.append(float(np.linalg.norm(original-reduced)/(1+np.linalg.norm(reduced))))
    assert raw_errors[1]<5e-7,raw_errors
    assert raw_errors[1]<.4*raw_errors[0]+1e-10,raw_errors
    principal=max(float(np.linalg.norm(s@mm[j]@s-sig(tt[:,j]))) for j in range(3))
    div_t=np.array([sum(dt[j,i,j] for j in range(3)) for i in range(3)])
    green=float(np.linalg.norm(c-c.conj().T+1j*sig(div_t)))
    ps=s@val
    flux=float(abs(np.vdot(ps,(sig(normal)-np.dot(b,normal)*ID)@ps)-np.vdot(val,sig(tt@normal)@val)))
    assert max(principal,green,flux)<2e-11,(principal,green,flux)
    normal_size=float(np.linalg.norm(tt@normal))
    assert 1/gam-1e-12<=normal_size<=1+1e-12
    return dict(speed=float(np.linalg.norm(b)),q=q,xi=xi,
                direct_refinement_errors=raw_errors,principal_error=principal,
                green_identity_error=green,flux_error=flux,transformed_normal_size=normal_size)


def state_case(speed,cutoff):
    # One particle polarization, beta=1. Compare analytic angular integration
    # with direct adaptive angular quadrature; neither uses the cap asymptote.
    def angular(t):
        aa=cutoff*t
        if aa==0:return 1.0
        return float((np.logaddexp(0,aa*(speed-1))-np.logaddexp(0,-aa*(1+speed)))/(aa*speed))
    analytic=quad(lambda t:t**3*angular(t),0,1,epsabs=1e-11,epsrel=1e-10)[0]/(4*np.pi**2)
    def nested(t):
        return t**3*quad(lambda mu:expit(-cutoff*t*(1-speed*mu)),-1,1,
                         points=[1/speed],epsabs=1e-11,epsrel=1e-10)[0]
    direct=quad(nested,0,1,epsabs=1e-11,epsrel=1e-10)[0]/(4*np.pi**2)
    error=abs(analytic-direct)
    assert error<1e-9,error
    limit=(1-1/speed)/(16*np.pi**2)
    return dict(speed=speed,momentum_cutoff=cutoff,energy_over_cutoff4=analytic,
                direct_quadrature=direct,quadrature_difference=error,limit_coefficient=limit,
                distance_to_limit=abs(analytic-limit))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    exact,exponents=exact_checks()
    rng=np.random.default_rng(20261010)
    differential=[]
    for speed in [.02,.2,.7]:
        for orientation in range(2):
            x=rng.normal(size=3);direction=rng.normal(size=3);direction/=np.linalg.norm(direction)
            omega=direction*speed/np.linalg.norm(np.cross(direction,x))
            k=rng.normal(size=3)*2
            normal=rng.normal(size=3);normal/=np.linalg.norm(normal)
            for q in [1,1000]:
                for xi in [.05,4]:differential.append(direct_case(x,omega,k,q,xi,normal))
    state=[state_case(v,cutoff) for v in [1.1,1.5,2] for cutoff in [8,32,128]]
    for speed in [1.1,1.5,2]:
        errors=[s['distance_to_limit'] for s in state if s['speed']==speed]
        assert all(b<a for a,b in zip(errors,errors[1:])),errors
    result=dict(status='passed_rotating_local_operator_and_state_diagnostic',
        versions={'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sp.__version__},
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        exact_checks=exact,exact_exponents=exponents,differential_cases=differential,state_cases=state,
        local_boost_identity_verified=True,state_obstruction_is_known_result=True,
        complete_rotating_return_proved=False,physical_hadamard_state_constructed=False,
        full_energy_proved=False,novelty_established=False,lean_invocations=0)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],exact_checks=len(exact),exponents=exponents,
        differential_cases=len(differential),state_cases=len(state),
        max_differential_error=max(c['direct_refinement_errors'][1] for c in differential),
        max_green_error=max(c['green_identity_error'] for c in differential),
        max_quadrature_difference=max(c['quadrature_difference'] for c in state))))

if __name__=='__main__':main()
