"""Supplementary checks for the projected electric mouth boundary limit.
No compact-sector/state/stress or physical stability claim. No Lean invocation.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import hashlib, json, platform, time
import numpy as np
import scipy
from scipy.linalg import cho_factor, cho_solve, cholesky, eigvalsh
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]

def constants(a,eps,D):
    T=np.pi-2*eps;h=T/4
    beta=a*h/(2*eps*(eps+h));gamma=a*np.log1p(h/eps)
    C=2*(1/T+np.mean(1/np.asarray(D)))+(T+np.mean(D))/2
    c=np.pi**2*a
    bound=min(1.,max(c/(1+c),gamma/(1-gamma)+C/(beta+C))) if gamma<1 else None
    return dict(T=float(T),beta=float(beta),gamma=float(gamma),C=float(C),bound=bound)

def symbolic():
    R,lam,g,n,b,a=sp.symbols('R lam g n b a',positive=True)
    x=sp.symbols('x',positive=True)
    F=lam*R+1/(1+lam*R)-1/b
    E=g*n*n*R/(4*sp.pi)+g*n*n*(1/(1+lam*R)-1/b)/(4*sp.pi*lam)
    tests={
      'paired_Maxwell_energy':sp.simplify(E-n*n*g*F/(4*sp.pi*lam)),
      'voltage_energy_identity':sp.simplify(sp.diff(E,n)*n/2-E),
      'sphere_optical_mass_integral':sp.simplify(sp.integrate(a/(sp.cos(x)**2),(x,-sp.atan(R),sp.atan(R)))-2*a*R),
      'electric_inertia_scaling':sp.simplify(2*sp.pi/(16*sp.pi*x**3)/(x**-2)-1/(8*x)),
      'angular_inertia_scaling':sp.simplify(sp.diff(16*n*n/x,n)-n/(x/32)),
      'overlap_cancellation':sp.simplify(lam*R+1/(1+lam*R)-1-(lam*R)**2/(1+lam*R))}
    assert all(z==0 for z in tests.values())
    return dict(count=len(tests),identities={k:str(v) for k,v in tests.items()})

def capacitance():
    out=[]
    for lam,R,b,g,n in [(1e-3,17,30,.02,2),(1e-5,200,100,.03,-3),(2e-4,130,15,.001,1)]:
        F=lam*R+1/(1+lam*R)-1/b
        # Energy density and independently integrated electric voltage.
        ei=quad(lambda r:g*n*n/(8*np.pi),-R,R)[0]
        eo=2*quad(lambda r:g*n*n/(8*np.pi*lam*r*r),1+lam*R,b,epsabs=1e-10)[0]
        voltage=quad(lambda r:g*n/(4*np.pi),-R,R)[0]
        voltage+=2*quad(lambda r:g*n/(4*np.pi*lam*r*r),1+lam*R,b,epsabs=1e-10)[0]
        expected=g*n*n*F/(4*np.pi*lam)
        err=max(abs(ei+eo-expected),abs(n*voltage/2-expected))/expected
        assert err<1e-12
        out.append(dict(lam=lam,R=R,b=b,g_squared=g,n=n,energy=expected,inertia=2*np.pi*lam/(g*F),relative_error=err))
    return out

def matrix_identity():
    rng=np.random.default_rng(80317);out=[]
    for dim in [4,9,15]:
      Z=rng.normal(size=(dim,dim))+1j*rng.normal(size=(dim,dim))
      H=Z.conj().T@Z/dim;R=np.linalg.inv(H+np.eye(dim))
      B=rng.normal(size=(2,dim))+1j*rng.normal(size=(2,dim))
      G=B@R@B.conj().T;C=1/eigvalsh(G)[0]
      RD=R-R@B.conj().T@np.linalg.solve(G,B@R)
      for beta in [.003,.3,30,3000]:
        direct=np.linalg.inv(H+np.eye(dim)+beta*B.conj().T@B)
        formula=R-R@B.conj().T@np.linalg.solve(G+np.eye(2)/beta,B@R)
        err=np.linalg.norm(direct-formula,2)
        ev=eigvalsh((direct-RD+direct.conj().T-RD.conj().T)/2)
        bound=C/(beta+C)
        assert err<3e-11 and ev[0]>-3e-11 and ev[-1]<=bound+3e-11
        out.append(dict(dim=dim,beta=beta,identity_error=err,observed_norm=ev[-1],bound=bound))
    return out

def continuum_lift():
    out=[]
    for q in [1,2,7,31,100]:
      for width in [2,8]:
        D=np.geomspace(10.**(-width/2),10.**(width/2),q);T=2.8;L=T+D
        g0=np.mean((1+np.exp(-L))/(2*(-np.expm1(-L))))
        g1=np.mean((np.exp(-T)+np.exp(-D))/(2*(-np.expm1(-L))))
        C=2*(1/T+np.mean(1/D))+(T+np.mean(D))/2
        smallest=g0-g1
        assert C*smallest>=1-1e-12
        out.append(dict(q=q,D_min=D.min(),D_max=D.max(),C=C,G_min=smallest,margin=C*smallest))
    return out

def assemble(q,D,a,eps,n):
    """Conforming linear elements; geometric mouth mesh and Gauss integration."""
    T=np.pi-2*eps
    left=np.geomspace(eps,eps+T/2,n+1)-eps
    nodes=np.r_[left,T-left[-2::-1]];ni=len(nodes)-1
    ne=max(8,n//2);nc=ni+ne;size=nc*q
    K=np.zeros((size,size));M=K.copy();V=K.copy()
    gx,gw=np.polynomial.legendre.leggauss(16)
    def element(ids,length):
        K[np.ix_(ids,ids)]+=np.array([[1,-1],[-1,1]])/length
        M[np.ix_(ids,ids)]+=length*np.array([[2,1],[1,2]])/6
    for j in range(q):
      off=j*nc
      for k in range(ni):element([off+k,off+k+1],nodes[k+1]-nodes[k])
      ext=[off+ni,*range(off+ni+1,off+nc),off]
      for k in range(ne):element(ext[k:k+2],D[j]/ne)
    for k in range(ni):
      x0,x1=nodes[k:k+2];length=x1-x0;t=(gx+1)/2;x=x0+length*t
      weight=a/np.sin(eps+np.minimum(x,T-x))**2*gw*length/2
      N=np.vstack((1-t,t));local=(N*weight)@N.T/q
      for j in range(q):
       for l in range(q):V[np.ix_([j*nc+k,j*nc+k+1],[l*nc+k,l*nc+k+1])]+=local
    B=np.zeros((2,size))
    B[0,np.arange(q)*nc]=q**-.5;B[1,np.arange(q)*nc+ni]=q**-.5
    return K,M,V,B

def fem_case(q,D,eps,n):
    a=eps**(1/3);K,M,V,B=assemble(q,D,a,eps,n)
    factor0=cho_factor(K+M);factor=cho_factor(K+M+V)
    S=cholesky(M,lower=True)
    R0=S.T@cho_solve(factor0,S)
    U=cho_solve(factor0,B.T);G=B@U
    RD=R0-(S.T@U)@np.linalg.solve(G,U.T@S)
    Ra=S.T@cho_solve(factor,S)
    diff=(Ra-RD+Ra.T-RD.T)/2
    lo=eigvalsh(diff,subset_by_index=[0,0])[0]
    hi=eigvalsh(diff,subset_by_index=[len(diff)-1,len(diff)-1])[0]
    c=constants(a,eps,D);bound=c['bound'];norm=max(abs(lo),abs(hi))
    assert norm<=bound+1e-7
    # Both form inequalities are tested independently on arbitrary vectors.
    rng=np.random.default_rng(934+q+n)
    min_lower=1e300;max_hardy=0.
    for _ in range(12):
      f=rng.normal(size=len(M))
      lower=f@(V+c['gamma']*K-c['beta']*(B.T@B))@f
      min_lower=min(min_lower,lower)
      f=f-B.T@np.linalg.solve(B@B.T,B@f)
      ratio=(f@V@f)/(f@K@f)
      max_hardy=max(max_hardy,ratio)
    assert min_lower>=-1e-7 and max_hardy<=np.pi**2*a+1e-7
    top=eigvalsh((Ra+Ra.T)/2,subset_by_index=[len(Ra)-q-3,len(Ra)-1])[::-1]
    zero_error=max(abs(top[:q-1]-1)) if q>1 else 0.
    assert zero_error<4e-7
    return dict(q=q,D=list(D),epsilon=eps,a=a,mesh_half_cells=n,dofs=len(M),
        constants=c,norm=norm,eigenvalue_interval=[lo,hi],lower_form_min=min_lower,
        Hardy_ratio_max=max_hardy,neutral_zero_error=zero_error,
        first_positive_eigenvalue=1/top[q-1]-1)

def finite_elements():
    out=[]
    for D in [[.35,.35],[.17,.8,2.4],[.13,.4,.7,1.2,3.]]:
      eps_values=[1e-2,1e-4,1e-6] if len(D)<5 else [1e-4]
      for eps in eps_values:
        meshes=[32,64,128] if len(D)==2 else [16,32,64]
        runs=[fem_case(len(D),D,eps,n) for n in meshes]
        difference=abs(runs[-1]['norm']-runs[-2]['norm'])
        assert difference<.004*runs[-1]['norm']+8e-6,(D,eps,difference)
        out.append(dict(q=len(D),D=D,epsilon=eps,refinements=runs,successive_norm_difference=difference))
        print('FEM',len(D),eps,runs[-1]['norm'],difference,flush=True)
    return out

def independent_ODE(fem):
    out=[]
    for case in fem:
      if case['q']!=2:continue
      eps=case['epsilon'];a=eps**(1/3);T=np.pi-2*eps;D=case['D'][0]
      roots=[]
      for tol in [2e-9,2e-11]:
        def residual(k):
          sol=solve_ivp(lambda s,y:[y[1],(a/np.cos(s)**2-k*k)*y[0]],
                        [0,T/2],[1.,0.],method='DOP853',rtol=tol,atol=tol*.01)
          assert sol.success
          y,dy=sol.y[:,-1]
          return dy*np.cos(k*D/2)-k*y*np.sin(k*D/2)
        roots.append(brentq(residual,.5,1.7,xtol=2e-12))
      errors=[abs(r['first_positive_eigenvalue']-roots[-1]**2) for r in case['refinements']]
      assert abs(roots[0]-roots[1])<1e-8
      assert errors[-1]<.002 and errors[0]/errors[-1]>10
      out.append(dict(epsilon=eps,ODE_frequencies=roots,FEM_eigenvalue_errors=errors,
                      error_improvement=errors[0]/errors[-1]))
    for q in [2,3]:
      rows=[x for x in fem if x['q']==q]
      assert rows[-1]['refinements'][-1]['norm']<.08*rows[0]['refinements'][-1]['norm']
    return out

def scales():
    out=[];kap=1+np.arcsinh(np.sqrt(3))/(2*np.sqrt(3))
    for q in [1e4,1e8,1e16,1e32]:
      eps=np.arctan(q**-2);a=1/(4*np.pi**2*q);T=np.pi-2*eps;h=T/4
      delta=np.log(q)**-2;Dbar=kap*delta*(np.log(q)+2+2*np.log(2))
      C=2*(1/T+1/delta)+(T+Dbar)/2
      beta=a*h/(2*eps*(eps+h));gam=a*np.log1p(h/eps);c=np.pi**2*a
      E=min(1.,max(c/(1+c),gam/(1-gam)+C/(C+beta)))
      z=(1+1j*delta)**2
      complex_bound=E*(1+abs(1+z)/(2*delta))**2
      out.append(dict(q=q,epsilon=eps,a=a,integrated_potential=2*a*q*q,
         overlap_parameter=1/(16*np.pi*q),delta=delta,C_upper=C,beta=beta,gamma=gam,
         resolvent_bound=E,frequency_one_bound=complex_bound,
         rigid_singlet_rayleigh=2*a*q*q/(T+delta),charged_trial_bound=(np.pi/T)**2+2*a,
         bare_electric_inertia=1/(8*q),bare_relative_angular_inertia=q/32))
    assert out[-1]['resolvent_bound']<out[0]['resolvent_bound']
    return out

def main():
    start=time.time()
    data=dict(symbolic=symbolic(),capacitance=capacitance(),matrix_identity=matrix_identity(),
              continuum_lift=continuum_lift(),finite_elements=finite_elements(),scales=scales())
    data['independent_ODE']=independent_ODE(data['finite_elements'])
    inputs=['notes/gauss-mouth-boundary.md','notes/singlet-static-curvature.md',
            'notes/rotational-current-response.md','notes/stationary-electric-dipole.md']
    sha=lambda p:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    data.update(status='passed_projected_electric_boundary_limit_physical_state_pending',
         date='2026-10-11',script_sha256=sha(Path(__file__)),input_sha256={p:sha(ROOT/p) for p in inputs},
         versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__),
         elapsed_seconds=time.time()-start,compact_sectors_resolved=False,physical_exterior_photon_derived=False,
         renormalized_stress_controlled=False,physical_damping_proved=False,novelty_established=False,lean_invocations=0)
    (ROOT/'evidence/gauss-mouth-boundary.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status=data['status'],elapsed_seconds=data['elapsed_seconds'],
       max_matrix_error=max(x['identity_error'] for x in data['matrix_identity']),
       max_fem_refinement=max(x['successive_norm_difference'] for x in data['finite_elements']))))

if __name__=='__main__':main()
