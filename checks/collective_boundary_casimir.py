"""Supplementary checks of the ideal collective-boundary Casimir response.

Run with the pinned requirements.txt environment. No Lean or physical-state claim.
The commensurate-graph check uses an independently assembled unitary bond map.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from pathlib import Path
import hashlib
import json
import platform
import time
import itertools
import numpy as np
import scipy
from scipy.integrate import quad_vec
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def frame(q, rank, rng):
    if rank == 0:
        return np.zeros((q, 0))
    return np.linalg.qr(rng.normal(size=(q, rank)))[0]


def kernel(k, T, D, W):
    """log F, its two T derivatives, and log F minus the exterior limit."""
    r = W.shape[1]
    x = np.exp(-k*T)
    y = np.exp(-k*D)
    z = x*y
    den = -np.expm1(-k*(T+D))
    a = -np.expm1(-2*k*T)
    out = np.array([r*np.log(a), 2*r*k*x*x/a,
                    -4*r*k*k*x*x/a**2, r*np.log(a)])
    for sign in [1, -1]:
        num = 1+y if sign == 1 else -np.expm1(-k*D)
        c = num/den
        C = W.T@(c[:, None]*W)
        Cp = W.T@((-k*c*z/den)[:, None]*W)
        Cpp = W.T@((k*k*c*z*(1+z)/den**2)[:, None]*W)
        _, logdet = np.linalg.slogdet(C)
        _, logext = np.linalg.slogdet(W.T@(num[:, None]*W))
        first = np.linalg.solve(C, Cp)
        out += [logdet, np.trace(first),
                np.trace(np.linalg.solve(C, Cpp)-first@first), logdet-logext]
    return out


def energies(T, D, W, tol=2e-10):
    # Break at all relevant inverse lengths so a short exterior is resolved.
    cuts = np.unique(np.r_[0., 1/T, 1/np.asarray(D), 45/min(T, min(D))])
    result = np.zeros(4)
    error = 0.
    for a, b in zip(cuts[:-1], cuts[1:]):
        val, err = quad_vec(lambda k: kernel(k, T, D, W), a, b,
                            epsabs=tol, epsrel=tol, limit=250)
        result += val
        error += err
    return result/(2*np.pi), error/(2*np.pi)


def symbolic_checks():
    x,y,T,D,r,phi,q,dc=sp.symbols('x y T D r phi q dc', positive=True)
    F=(1-x*x)*(1+y)/(1-x*y)*(1-y)/(1-x*y)
    E=r*(-sp.pi/(24*T)-sp.pi/(24*D)+sp.pi/(6*(T+D)))
    K=sp.log((1-x*x)*(1-y*y)/(1-x*y)**2)
    op=lambda f:sp.simplify(-x*sp.diff(f,x))
    checks={
      'equal_return_determinant':sp.simplify(F-(1-x*x)*(1-y*y)/(1-x*y)**2),
      'first_equal_kernel':sp.simplify(op(K)-(2*x*x/(1-x*x)-2*x*y/(1-x*y))),
      'second_equal_kernel':sp.simplify(op(op(K))-(-4*x*x/(1-x*x)**2+2*x*y/(1-x*y)**2)),
      'equal_energy_first':sp.simplify(sp.diff(E,T)-r*(sp.pi/(24*T*T)-sp.pi/(6*(T+D)**2))),
      'clock_chain':sp.simplify(sp.diff(E.subs(T,T*sp.exp(-phi)),phi,2).subs(phi,0)-(T*sp.diff(E,T)+T*T*sp.diff(E,T,2))),
      'deep_integral':sp.simplify(2/(4*sp.pi**2*q)*q*q-q/(2*sp.pi**2)),
      'counterterm_first':sp.simplify(sp.diff(dc*sp.exp(phi)*q/(2*sp.pi**2),phi).subs(phi,0)-dc*q/(2*sp.pi**2)),
      'counterterm_second':sp.simplify(sp.diff(dc*sp.exp(phi)*q/(2*sp.pi**2),phi,2).subs(phi,0)-dc*q/(2*sp.pi**2))}
    assert all(v == 0 for v in checks.values())
    return {k:str(v) for k,v in checks.items()}


def compressed_vs_scattering():
    rng=np.random.default_rng(115903)
    rows=[]
    for q,r in [(1,1),(3,1),(5,2),(8,3),(8,8)]:
        W=frame(q,r,rng)
        P=W@W.T;Q=np.eye(q)-P
        S=np.block([[-P,Q],[Q,-P]])
        assert np.linalg.norm(S.T@S-np.eye(2*q))<2e-14
        D=np.geomspace(.13,7.1,q)
        for T,k in itertools.product([.7,2.9],[.02,.3,2.,7.]):
            L=np.diag(np.r_[np.full(q,np.exp(-k*T)),np.exp(-k*D)])
            sign,ld=np.linalg.slogdet(np.eye(2*q)-(S@L)@(S@L))
            direct=ld-2*np.sum(np.log(-np.expm1(-k*(T+D))))
            expected=kernel(k,T,D,W)[0]
            error=abs(direct-expected)
            assert sign>0 and error<3e-11
            rows.append(dict(q=q,rank=r,T=T,k=k,log_determinant_error=error))
    return rows


def equal_returns():
    rng=np.random.default_rng(4257);rows=[]
    for q,r in [(1,1),(7,2),(9,4)]:
        W=frame(q,r,rng)
        for T,D in [(2.3,.07),(.8,4.7),(np.pi,np.pi),(3.1,.6)]:
            val,err=energies(T,np.full(q,D),W)
            expected=np.array([r*(-np.pi/(24*T)-np.pi/(24*D)+np.pi/(6*(T+D))),
              r*(np.pi/(24*T*T)-np.pi/(6*(T+D)**2)),
              r*(-np.pi/(12*T**3)+np.pi/(3*(T+D)**3)),
              r*(-np.pi/(24*T)+np.pi/(6*(T+D)))])
            error=np.max(abs(val-expected))
            assert error<3e-9
            rows.append(dict(q=q,rank=r,T=T,D=D,max_error=error,quadrature_error=err))
    return rows


def bond_map(Tsteps,Dsteps,W):
    """Equilateral subdivision; no boundary covariance/determinant formula."""
    q=len(Dsteps);r=W.shape[1];P=W@W.T;Q=np.eye(q)-P
    mouth=np.block([[-P,Q],[Q,-P]])
    incident={0:[None]*(2*q),1:[None]*(2*q)}
    reverse=[];vertex=2
    for edge,count in enumerate([Tsteps]*q+list(Dsteps)):
        nodes=[0,*range(vertex,vertex+count-1),1];vertex+=count-1
        for a,b in zip(nodes[:-1],nodes[1:]):
            ia=len(reverse);ib=ia+1;reverse.extend([ib,ia])
            if a==0:incident[0][edge]=ia
            else:incident.setdefault(a,[]).append(ia)
            if b==1:incident[1][edge]=ib
            else:incident.setdefault(b,[]).append(ib)
    U=np.zeros((len(reverse),len(reverse)))
    for vertex,ids in incident.items():
        S=mouth if vertex in [0,1] else np.array([[0.,1.],[1.,0.]])
        U[np.ix_(ids,[reverse[i] for i in ids])]=S
    assert np.linalg.norm(U.T@U-np.eye(len(U)))<3e-14
    return U


def bond_energy(Tsteps,Dsteps,W,ell):
    U=bond_map(Tsteps,Dsteps,W)
    lam=np.linalg.eigvals(U)
    assert max(abs(abs(lam)-1))<2e-14
    # Each phase is one ladder of positive frequencies. Zero contributes zero.
    alpha=(-np.angle(lam)/(2*np.pi))%1
    return -np.pi/(2*ell)*np.sum(alpha*alpha-alpha+1/6),len(U)


def independent_spectra():
    rng=np.random.default_rng(618201);rows=[]
    for ds in [[1,3],[1,2,4],[2,3,5]]:
        q=len(ds)
        for r in [1,q]:
          W=frame(q,r,rng)
          for steps,ell in [(2,.7),(3,1.1)]:
            ew,n=bond_energy(steps,ds,W,ell)
            e0,_=bond_energy(steps,ds,np.zeros((q,0)),ell)
            free=-sum(np.pi/(6*(steps+np.array(ds))*ell))
            T=steps*ell;D=np.array(ds)*ell
            val,err=energies(T,D,W)
            error=abs((ew-e0)-val[0])
            assert abs(e0-free)<2e-12 and error<2e-10
            rows.append(dict(q=q,rank=r,throat_steps=steps,exterior_steps=ds,
              step_length=ell,bond_dimension=n,independent_energy=ew-e0,
              integral_energy=val[0],energy_error=error,free_error=abs(e0-free)))
    return rows


def bounds_and_derivatives():
    rng=np.random.default_rng(70419);rows=[]
    for q,r in [(2,1),(5,2),(19,3),(71,4)]:
      W=frame(q,r,rng)
      for dlo,dhi in [(1e-5,.002),(.2,7.),(100.,1e5)]:
        D=np.geomspace(dlo,dhi,q);T=2.7
        # Derivatives and exterior-subtracted energy require only k=O(1/T).
        vals,err=quad_vec(lambda k:kernel(k,T,D,W)[[1,2,3]],0.,45/T,
                          epsabs=3e-10,epsrel=3e-10,limit=300)
        vals/=2*np.pi;err/=2*np.pi
        b1=np.pi*r/(6*T*T);b2=2*np.pi*r/(3*T**3);be=np.pi*r/(6*T)
        assert abs(vals[0])<=b1+2e-9 and abs(vals[1])<=b2+2e-9
        assert abs(vals[2])<=be+2e-9
        max_ratio=0.
        for k in np.geomspace(1e-5/T,30/T,60):
            z=kernel(k,T,D,W);x=np.exp(-k*T);den=-np.expm1(-k*T)
            ratio=max(abs(z[1])/(2*r*k*x/den),abs(z[2])/(4*r*k*k*x/den**2))
            max_ratio=max(max_ratio,ratio)
        assert max_ratio<=1+2e-8
        rows.append(dict(q=q,rank=r,D_min=dlo,D_max=dhi,T=T,
           energy_first=vals[0],energy_second=vals[1],energy_without_exterior=vals[2],
           first_bound=b1,second_bound=b2,energy_bound=be,
           maximum_kernel_ratio=max_ratio,quadrature_error=err))
    # Differentiate only the integrated exterior-subtracted energy numerically.
    fd=[]
    for q,r in [(2,1),(5,2),(11,3)]:
      W=frame(q,r,rng);D=np.geomspace(.15,5.,q);T=2.3
      center=energies(T,D,W)[0]
      errors=[]
      for h in [0.02,0.01,0.005,0.0025]:
        lo=energies(T-h,D,W)[0][3];hi=energies(T+h,D,W)[0][3]
        e1=abs((hi-lo)/(2*h)-center[1])
        e2=abs((hi-2*center[3]+lo)/(h*h)-center[2])
        errors.append([e1,e2])
      assert max(errors[-1])<2e-7
      assert max(errors[0])/max(errors[-1])>10
      fd.append(dict(q=q,rank=r,steps=[.02,.01,.005,.0025],errors=errors))
    return rows,fd


def scales():
    rows=[]
    for q in [1e4,1e8,1e16,1e32]:
        delta=np.log(q)**-2;T=np.pi
        rows.append(dict(q=q,delta=delta,
          ideal_rank_one_clock_curvature_bound=5*np.pi/(6*T),
          ideal_bound_over_q_delta=(5*np.pi/(6*T))/(q*delta),
          local_unit_coefficient_curvature=q/(2*np.pi**2),
          local_curvature_over_q_delta=1/(2*np.pi**2*delta)))
    assert rows[-1]['ideal_bound_over_q_delta']<rows[0]['ideal_bound_over_q_delta']
    assert rows[-1]['local_curvature_over_q_delta']>rows[0]['local_curvature_over_q_delta']
    return rows


def main():
    start=time.time()
    data={}
    for name,check in [('symbolic',symbolic_checks),('scattering',compressed_vs_scattering),
                       ('equal_returns',equal_returns),('independent_spectra',independent_spectra)]:
        data[name]=check()
        print(json.dumps(dict(stage=name,cases=len(data[name]),status='passed')),flush=True)
    data['uniform_bounds'],data['finite_difference']=bounds_and_derivatives()
    data['scales']=scales()
    sha=lambda p:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    inputs=['notes/collective-boundary-casimir.md','notes/gauss-mouth-boundary.md',
            'notes/singlet-static-curvature.md','notes/rotational-current-response.md',
            'evidence/collective-boundary-casimir-source-review.json','requirements.txt']
    data.update(status='passed_ideal_boundary_response_local_matching_unresolved',date='2026-10-11',
       script_sha256=sha(Path(__file__)),input_sha256={p:sha(ROOT/p) for p in inputs},
       versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__),
       elapsed_seconds=time.time()-start,finite_smooth_mouth_stress_controlled=False,
       physical_renormalization_matched=False,angular_state_constructed=False,
       physical_damping_proved=False,novelty_established=False,lean_invocations=0)
    (ROOT/'evidence/collective-boundary-casimir.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status=data['status'],elapsed_seconds=data['elapsed_seconds'],
       counts={k:len(data[k]) for k in ['symbolic','scattering','equal_returns','independent_spectra','uniform_bounds','finite_difference','scales']},
       max_scattering_error=max(x['log_determinant_error'] for x in data['scattering']),
       max_spectral_error=max(x['energy_error'] for x in data['independent_spectra']))))


if __name__=='__main__':
    main()
