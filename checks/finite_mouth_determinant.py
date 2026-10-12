"""Finite-cutoff Gaussian determinant checks; no Lean or physical state claim."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import hashlib,json,time,platform
from functools import lru_cache
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from numpy.polynomial.legendre import leggauss
import mpmath as mp
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
mp.mp.dps=55
H=.125

def symbolic():
    x,e,a,h,k,L=sp.symbols('x e a h k L',positive=True)
    b=sp.sqrt(1-e*e*x*x);T=2*sp.acos(e*x)/x
    D=lambda v:sp.simplify(x*sp.diff(v,x))
    checks={'T_first':D(T)+T+2*e/b,
      'T_second':D(D(T))-T-2*e/b+2*e**3*x*x/b**3,
      'potential_split':2*a*b/e-2*a*(1/e-1/(e+h))-2*a*(1/(e+h)-e*x*x/(1+b))}
    plus=k*sp.tanh(k*L/2);minus=k/sp.tanh(k*L/2)
    checks.update(K_plus_first=sp.diff(plus,L)/plus-k/sp.sinh(k*L),
      K_minus_first=sp.diff(minus,L)/minus+k/sp.sinh(k*L),
      K_plus_second=sp.diff(plus,L,2)/plus+k*k/(2*sp.cosh(k*L/2)**2),
      K_minus_second=sp.diff(minus,L,2)/minus-k*k/(2*sp.sinh(k*L/2)**2))
    checks={name:sp.simplify(sp.trigsimp(v.rewrite(sp.exp))) for name,v in checks.items()}
    assert all(v==0 for v in checks.values()),checks
    integrals=[mp.quad(lambda t:t/mp.sinh(t),[0,1,mp.inf]),
      mp.quad(lambda t:t*t/(2*mp.sinh(t/2)**2)+2*t*t/mp.sinh(t)**2,[0,1,mp.inf])]
    assert abs(integrals[0]-mp.pi**2/4)<mp.mpf('1e-45')
    assert abs(integrals[1]-mp.pi**2)<mp.mpf('1e-45')
    return {name:str(v) for name,v in checks.items()},[float(v) for v in integrals]

def split_checks():
    rows=[]
    for eps in map(mp.mpf,['.02','.00001','.000000001']):
      for x in map(mp.mpf,['.5','1','2']):
        aa=mp.mpf('.01');hh=mp.mpf('.125');b=mp.sqrt(1-(eps*x)**2)
        half=mp.acos(eps*x)/x
        def w(d,ph=0):
            xx=x*mp.exp(ph);bb=mp.sqrt(1-(eps*xx)**2)
            p=eps*mp.cos(xx*d)+bb*mp.sin(xx*d)/xx
            return aa/p**2-aa/(eps+d)**2
        points=[eps*t for t in map(mp.mpf,['.01','.1','1','10']) if eps*t<hh]+[hh/2,hh]
        ratios=[]
        for d in points:
            value=w(d);assert value>0
            ratios.append([float(abs(mp.diff(lambda ph:w(d,ph),0,j))/value) for j in [1,2]])
            assert value<2*aa*x*x*d/(eps+d)
        eta=2*aa*(1/(eps+hh)-eps*x*x/(1+b))
        # Integrate the small difference directly, avoiding subtraction of large integrals.
        cells=[mp.mpf(0)]+sorted(set(p for p in [eps,10*eps,100*eps,hh] if p<=hh))
        direct=2*mp.quad(lambda d:w(d),cells)+2*mp.quad(lambda s:aa*x*x/mp.cos(x*s)**2,[0,half-hh])
        error=abs(direct-eta)
        assert error<mp.mpf('1e-38')
        assert max(max(r) for r in ratios)<8
        rows.append(dict(epsilon=float(eps),x=float(x),M=float(2*aa*(1/eps-1/(eps+hh))),eta=float(eta),
          integral_error=float(error),max_log_derivative_ratio=max(max(r) for r in ratios)))
    return rows

def free_scalar(k,length,order=0):
    u=k*length/2;t=np.tanh(u)
    if order==0:p,m=k*t,k/t
    elif order==1:p,m=k*k*(1-t*t)/2,-k*k*(1/t**2-1)/2
    else:p,m=-k**3*(1-t*t)*t/2,k**3*(1/t**2-1)/t/2
    return np.array([[(p+m)/2,(p-m)/2],[(p-m)/2,(p+m)/2]])

def free_matrix(k,length,q,order=0):return np.kron(free_scalar(k,length,order),np.eye(q))

@lru_cache(maxsize=4096)
def cell_scalar(k,a,eps):
    kk,aa,ee=map(mp.mpf,map(str,[k,a,eps]));nu=mp.sqrt(mp.mpf('.25')+aa)
    def basis(t):
        z=kk*t
        iv=mp.besseli(nu,z);kv=mp.besselk(nu,z)
        ip=(mp.besseli(nu-1,z)+mp.besseli(nu+1,z))/2
        kp=-(mp.besselk(nu-1,z)+mp.besselk(nu+1,z))/2
        return mp.matrix([[mp.sqrt(t)*iv,mp.sqrt(t)*kv],
          [iv/(2*mp.sqrt(t))+mp.sqrt(t)*kk*ip,kv/(2*mp.sqrt(t))+mp.sqrt(t)*kk*kp]])
    transfer=basis(ee+mp.mpf('.125'))*basis(ee)**-1
    A,B,C,D=transfer[0,0],transfer[0,1],transfer[1,0],transfer[1,1]
    assert abs(A*D-B*C-1)<mp.mpf('1e-40')
    return np.array([[float(A/B),float(-1/B)],[float(-1/B),float(D/B)]])

def return_matrix(k,returns):
    n=len(returns);out=np.zeros((2*n,2*n))
    for j,d in enumerate(returns):out[np.ix_([j,n+j],[j,n+j])]=free_scalar(k,float(d))
    return out

def outside(k,a,eps,returns):
    q=len(returns);proj=np.ones((q,q))/q
    free=free_scalar(k,H);charged=cell_scalar(float(k),float(a),float(eps))
    cell=np.kron(free,np.eye(q))+np.kron(charged-free,proj)
    mat=np.zeros((4*q,4*q))
    # Cell endpoint order is mouth then inner; the left cell is traversed backwards.
    left=list(range(q,2*q))+list(range(q));right=list(range(2*q,4*q))
    mat[np.ix_(left,left)]+=cell;mat[np.ix_(right,right)]+=cell
    mid=list(range(q,3*q));ends=list(range(q))+list(range(3*q,4*q))
    mat[np.ix_(mid,mid)]+=return_matrix(k,returns)
    N=mat[np.ix_(ends,ends)]-mat[np.ix_(ends,mid)]@np.linalg.solve(mat[np.ix_(mid,mid)],mat[np.ix_(mid,ends)])
    N0=return_matrix(k,np.asarray(returns)+2*H)
    return (N+N.T)/2,N0

def force_data(k,L,N,N0):
    q=len(N)//2;K=free_matrix(k,L,q);Kp=free_matrix(k,L,q,1);Kpp=free_matrix(k,L,q,2)
    R=np.linalg.inv(K+N);R0=np.linalg.inv(K+N0);delta=R0-R
    ld=np.linalg.slogdet(K+N)[1]-np.linalg.slogdet(K+N0)[1]
    first=-np.trace(delta@Kp)
    second=-np.trace(delta@Kpp)+np.trace(R0@Kp@R0@Kp-R@Kp@R@Kp)
    return np.array([ld,first,second])

def boundary_checks():
    rows=[]
    for q in [1,3,9]:
      returns=np.geomspace(.01,40,q) if q>1 else np.array([.7])
      for a,eps in [(.003,.02),(.02,.001),(.01,.00001)]:
       for k in [.002,.1,1.,6.]:
        N,N0=outside(k,a,eps,returns);ev=np.linalg.eigvalsh(N-N0)
        rank=int(np.sum(ev>1e-8));assert rank<=4
        assert ev[0]>-2e-9
        for L in [1.1,2.7]:
            vals=force_data(k,L,N,N0);r=min(4,2*q)
            c1=k/np.sinh(k*L);c2=k*k/(2*np.sinh(k*L/2)**2)
            ratios=[abs(vals[1])/(r*c1),abs(vals[2])/(r*(c2+2*c1*c1))]
            assert max(ratios)<1+1e-7,ratios
            rows.append(dict(q=q,a=a,epsilon=eps,k=k,L=L,rank=rank,minimum_increment_eigenvalue=float(ev[0]),bound_ratios=ratios))
    return rows

def precise_boundary_logdet(k,L,N,N0):
    # Near the zero oscillator, forming K in double precision loses small
    # symmetric eigenvalues before the second finite difference is taken.
    kk,ll=mp.mpf(float(k)),mp.mpf(float(L));q=len(N)//2
    plus=kk*mp.tanh(kk*ll/2);minus=kk/mp.tanh(kk*ll/2)
    K=mp.zeros(2*q)
    for j in range(q):
        K[j,j]=K[q+j,q+j]=(plus+minus)/2
        K[j,q+j]=K[q+j,j]=(plus-minus)/2
    return float(mp.log(mp.det(K+mp.matrix(N.tolist())))-mp.log(mp.det(K+mp.matrix(N0.tolist()))))

def force_integrals():
    rows=[]
    for a,eps,returns in [(.01,.001,np.array([.03,.7,25.])),(.005,.00001,np.array([.01,.03,.3,3.,30.]))]:
      L=2.7;hist=[]
      for n in [16,32]:
        xx,ww=leggauss(n);data=np.zeros(3);fds={h:np.zeros(2) for h in [.01,.005]}
        for lo,hi in zip([0,.05,.5,2,8],[.05,.5,2,8,24]):
          for k,w in zip((lo+hi)/2+(hi-lo)*xx/2,(hi-lo)*ww/2):
            N,N0=outside(float(k),a,eps,returns);val=force_data(k,L,N,N0)
            log_at=lambda length:precise_boundary_logdet(k,length,N,N0) if k<.05 else force_data(k,length,N,N0)[0]
            val[0]=log_at(L)
            inf=np.linalg.slogdet(N+k*np.eye(len(N)))[1]-np.linalg.slogdet(N0+k*np.eye(len(N0)))[1]
            # One lifted constant mode gives F(L)-F(infinity)=-log(k)+O(1).
            # Integrate that known endpoint term exactly on the first panel.
            energy=val[0]-inf+(np.log(k) if lo==0 else 0)
            data+=w*np.array([energy,val[1],val[2]])/(2*np.pi)
            for h in fds:
                plus=log_at(L+h);minus=log_at(L-h)
                fds[h]+=w*np.array([(plus-minus)/(2*h),(plus-2*val[0]+minus)/h**2])/(2*np.pi)
        data[0]+=.05*(1-np.log(.05))/(2*np.pi)
        hist.append(dict(nodes_per_interval=n,energy_and_derivatives=data.tolist(),finite_difference_errors={str(h):abs(v-data[1:]).tolist() for h,v in fds.items()}))
      err=max(abs(np.asarray(hist[0]['energy_and_derivatives'])-hist[1]['energy_and_derivatives']))
      assert err<2e-7,err
      assert max(hist[-1]['finite_difference_errors']['0.005'])<2e-7
      assert all(hist[-1]['finite_difference_errors']['0.01'][j]/hist[-1]['finite_difference_errors']['0.005'][j]>3 for j in [0,1])
      r=min(4,2*len(returns));bounds=[r*np.pi/(8*L*L),r*np.pi/(2*L**3)]
      assert all(abs(hist[-1]['energy_and_derivatives'][j+1])<bounds[j] for j in [0,1])
      rows.append(dict(a=a,epsilon=eps,returns=returns.tolist(),refinements=hist,quadrature_difference=float(err),derivative_bounds=bounds))
    return rows

def full_ode(k,a,eps,x,returns):
    T=2*np.arccos(eps*x)/x;half=T/2
    def rhs(s,y):
        w=k*k+a*x*x/np.cos(x*s)**2
        return [y[1],w*y[0],y[3],w*y[2]]
    sol=solve_ivp(rhs,[0,half],[1,0,0,1],method='DOP853',rtol=2e-12,atol=2e-13)
    assert sol.success
    e,ep,o,op=sol.y[:,-1];lp,lm=ep/e,op/o
    Nt=np.array([[(lp+lm)/2,(lp-lm)/2],[(lp-lm)/2,(lp+lm)/2]])
    q=len(returns);Nf=free_matrix(k,T,q);Ns=Nf+np.kron(Nt-free_scalar(k,T),np.ones((q,q))/q)
    ext=return_matrix(k,returns)
    return float(np.log(2*e*o/(np.sinh(k*T)/k))+np.linalg.slogdet(Ns+ext)[1]-np.linalg.slogdet(Nf+ext)[1])

def full_nystrom(k,a,eps,x,returns,n):
    T=2*np.arccos(eps*x)/x;half=T/2;edges=[0.,eps]
    while edges[-1]*4<half:edges.append(edges[-1]*4)
    edges.append(half);g,w=leggauss(n);points=[];weights=[]
    for lo,hi in zip(edges[:-1],edges[1:]):
        dd=(lo+hi)/2+(hi-lo)*g/2;ww=(hi-lo)*w/2
        points.extend(dd);weights.extend(ww);points.extend(T-dd);weights.extend(ww)
    s=np.asarray(points);weights=np.asarray(weights)
    dist=abs(s[:,None]-s[None,:]);G=np.zeros_like(dist)
    for ret in returns:
        length=T+ret
        G+=(np.exp(-k*dist)+np.exp(-k*(length-dist)))/(2*k*(-np.expm1(-k*length)))/len(returns)
    v=a*x*x/np.cos(x*(s-half))**2;fac=np.sqrt(weights*v)
    matrix=np.eye(len(s))+fac[:,None]*G*fac[None,:]
    chol=np.linalg.cholesky(matrix)
    return float(2*np.log(np.diag(chol)).sum()),len(s)

def full_checks():
    rows=[]
    for eps in [.02,.002]:
      a=.01;x=.9;returns=np.array([.03,.7,12.])
      for k in [.1,1.,4.]:
        exact=full_ode(k,a,eps,x,returns);ref=[]
        for n in [12,24,48]:
            val,nodes=full_nystrom(k,a,eps,x,returns,n)
            ref.append(dict(order=n,nodes=nodes,logdet=val,error=abs(val-exact)))
        assert ref[-1]['error']<2e-5,(eps,k,ref)
        assert ref[-1]['error']<.09*ref[0]['error'],(eps,k,ref)
        rows.append(dict(epsilon=eps,k=k,integrated_potential=2*a*np.sqrt(1-(eps*x)**2)/eps,ode_logdet=exact,refinements=ref))
    return rows

def main():
    start=time.time();symbolic_result,integrals=symbolic()
    data=dict(symbolic=symbolic_result,universal_envelope_integrals=integrals,split=split_checks())
    print('Symbolic and positive small-remainder checks passed',flush=True)
    data['boundary']=boundary_checks();print('Finite barrier boundary ranks and force inequalities passed',flush=True)
    data['force_integrals']=force_integrals();print('Integrated force and derivative refinements passed',flush=True)
    data['full_determinants']=full_checks();print('Independent full-profile ODE and Fredholm determinants passed',flush=True)
    sha=lambda p:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    inputs=['notes/finite-mouth-determinant.md','notes/proper-mouth-local-stress.md','notes/singlet-static-curvature.md',
      'notes/gauss-mouth-boundary.md','evidence/finite-mouth-determinant-source-review.json','requirements.txt']
    data.update(status='passed_finite_product_gaussian_static_bound_physical_response_pending',date='2026-10-11',
      script_sha256=sha(Path(__file__)),input_sha256={p:sha(ROOT/p) for p in inputs},
      versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__,mpmath=mp.__version__,mpmath_digits=mp.mp.dps),
      elapsed_seconds=time.time()-start,finite_product_gaussian_static_premise_controlled=True,
      corrected_background_nonlocal_response_controlled=False,physical_moving_clock_vertex_derived=False,
      physical_damping_proved=False,novelty_established=False,lean_invocations=0)
    (ROOT/'evidence/finite-mouth-determinant.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status=data['status'],elapsed_seconds=data['elapsed_seconds'],counts={k:len(data[k]) for k in ['symbolic','split','boundary','force_integrals','full_determinants']},
       max_full_determinant_error=max(r['refinements'][-1]['error'] for r in data['full_determinants']))))

if __name__=='__main__':main()
