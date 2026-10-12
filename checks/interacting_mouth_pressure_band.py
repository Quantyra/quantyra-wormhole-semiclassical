"""Finite-band pressure proof diagnostics; no ultraviolet or physical claim.

Run with the pinned SpaceProofs Python. All matrices below are supplementary
conforming finite elements; analytic estimates, not fits, supply uniformity.
"""
from __future__ import annotations
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import hashlib,json,time,platform
import numpy as np
import scipy
from scipy.linalg import eigh,solve,cho_factor,cho_solve,null_space
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import eigsh
from scipy.integrate import solve_ivp,quad
from scipy.optimize import brentq
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]

def require(value,message):
    if not value:raise AssertionError(message)

def relative(x,y):
    return float(np.linalg.norm(x-y)/max(1.,np.linalg.norm(y)))

def exact_checks():
    u,z,w,k,d,a,b,eps=sy.symbols('u z w k d a b eps',positive=True)
    C,B=sy.symbols('C B',real=True)
    out={
      'ordered_pair_normalization':sy.simplify((C-w*k*B)**2/(8*w*k)-w*k*(B-C/(w*k))**2/8),
      'mouth_potential_capacity_scale':sy.simplify((eps/b)*(a/eps**2)*(eps/a)-1/b),
      'derivative_pressure_scale':sy.simplify((eps/b)**2*(a/eps)**2-a*a/(b*b)),
      'resolvent_factorization':sy.expand((w-k)*(w+k)-(w*w-k*k)),
      'neutral_value_identity':sy.simplify(sy.integrate(2*(1-z)*sy.cos(k*z),(z,0,1))
           -sy.integrate((1-z)**2*(-k*sy.sin(k*z)),(z,0,1))-1),
      'neutral_derivative_identity':sy.simplify(sy.integrate(2*(1-z)*k*sy.cos(k*z),(z,0,1))
           +k*k*sy.integrate((1-z)**2*sy.sin(k*z),(z,0,1))-k),
      'stieltjes_soft_weight':sy.simplify(u**2/u+sy.integrate(z**2/z**2,(z,0,u))-2*u),
    }
    require(all(x==0 for x in out.values()),out)
    return {k:str(v) for k,v in out.items()}

def assemble(D,a,eps,n):
    """Exact geometric cut; 16-point positive Gauss potential quadrature."""
    q=len(D);A=np.arccos(eps);T=2*A
    left=np.geomspace(eps,eps+A,n+1)-eps
    nodes=np.r_[left,T-left[-2::-1]];nt=len(nodes)
    ns=[max(8,int(np.ceil(n*d/T))) for d in D]
    sizes=[nt+m-1 for m in ns];offsets=np.r_[0,np.cumsum(sizes)]
    size=offsets[-1];K=np.zeros((size,size));M=K.copy();KT=K.copy();MT=K.copy();V=K.copy()
    gx,gw=np.polynomial.legendre.leggauss(16)
    shapes=np.vstack(((1-gx)/2,(1+gx)/2))
    def element(ids,h,throat):
        ke=np.array([[1.,-1.],[-1.,1.]])/h
        me=h*np.array([[2.,1.],[1.,2.]])/6
        K[np.ix_(ids,ids)]+=ke;M[np.ix_(ids,ids)]+=me
        if throat:KT[np.ix_(ids,ids)]+=ke;MT[np.ix_(ids,ids)]+=me
    for j,d in enumerate(D):
        off=offsets[j]
        for i in range(nt-1):element([off+i,off+i+1],nodes[i+1]-nodes[i],True)
        ext=[off+nt-1,*range(off+nt,offsets[j+1]),off]
        for i in range(ns[j]):element(ext[i:i+2],d/ns[j],False)
    for i in range(nt-1):
        lo,hi=nodes[i:i+2];zs=(lo+hi)/2+(hi-lo)*gx/2-A
        local=(shapes*(gw*(hi-lo)*a/(2*np.cos(zs)**2)))@shapes.T/q
        for j in range(q):
            for l in range(q):
                V[np.ix_([offsets[j]+i,offsets[j]+i+1],
                         [offsets[l]+i,offsets[l]+i+1])]+=local
    trace=np.zeros((2*q,size))
    for j in range(q):trace[j,offsets[j]]=1;trace[q+j,offsets[j]+nt-1]=1
    charged=np.kron(np.eye(2),np.ones((1,q))/np.sqrt(q))
    H=K+V
    # Exact neutral constant modes and a positive rank-(q-1) stabilization.
    constants=np.zeros((size,q))
    for j in range(q):constants[offsets[j]:offsets[j+1],j]=1
    Z=constants@null_space(np.ones((1,q)))
    if q>1:
        ze,zu=eigh(Z.T@M@Z);Z=Z@zu@np.diag(ze**-.5)
    stabilized=H+(M@Z)@(M@Z).T
    chtrace=charged@trace
    capacity=chtrace@cho_solve(cho_factor(stabilized),chtrace.T)
    # Local positive-energy DtN, extracted from the same conforming form.
    throat_ids=np.arange(nt)
    scalar=q*V[np.ix_(throat_ids,throat_ids)]+KT[np.ix_(throat_ids,throat_ids)]+MT[np.ix_(throat_ids,throat_ids)]
    ends=np.array([0,nt-1]);inside=np.arange(1,nt-1)
    poisson=np.zeros((nt,2));poisson[ends]=np.eye(2)
    poisson[inside]=-solve(scalar[np.ix_(inside,inside)],scalar[np.ix_(inside,ends)],assume_a='pos')
    dtn=scalar[np.ix_(ends,np.arange(nt))]@poisson
    return dict(K=K,M=M,H=H,KT=KT,MT=MT,V=V,trace=trace,charged=charged,
                capacity=capacity,dtn=dtn,nodes=nodes,A=A,zero=Z,poisson=poisson,
                scalar=scalar,nt=nt,offsets=offsets,ends=ends)

def cutoff(E,R):
    t=2*E/R-1;out=np.zeros_like(t);out[t<=0]=1
    use=(t>0)&(t<1);u=t[use]
    p=np.exp(-1/u);s=np.exp(-1/(1-u));out[use]=s/(p+s)
    return out

def fem_case(D,a,eps,n,R=4.):
    q=len(D);g=assemble(D,a,eps,n);H,M=g['H'],g['M']
    # Shift-invert resolves low frequencies without an absolute tolerance
    # scaled by the very large highest eigenvalue of the graded mesh.
    count=min(len(H)-2,int(np.ceil(q+R*(2*q*g['A']+sum(D))/np.pi))+4)
    eig,vec=eigsh(csc_matrix(H),M=csc_matrix(M),k=count,sigma=-.25,
                 which='LM',tol=2e-12,v0=np.linspace(.7,1.3,len(H)))
    order=np.argsort(eig);eig=eig[order];vec=vec[:,order]
    keep=eig<=R*R;eig=eig[keep];vec=vec[:,keep]
    zeros=eig[:q-1];eig=eig[q-1:];vec=vec[:,q-1:]
    require(len(eig)>0 and min(eig)>0,('positive spectrum',q,a,eps,n,eig[:3]))
    zero_error=max(abs(zeros),default=0.);require(zero_error<2e-5,('zero modes',zero_error))
    frequencies=np.sqrt(eig);trace=g['trace']@vec;U=g['charged']@trace
    # Weak endpoint flux is accurate under refinement even across strong layers.
    derivatives=g['trace']@((g['KT']+g['V'])@vec-(g['MT']@vec)*eig)
    p=g['charged']@derivatives
    cap=eigh(g['capacity'],eigvals_only=True);Mbar=a/eps
    require(max(cap)<=8/Mbar*(1+1e-7),('capacity',cap,8/Mbar))
    SU=(U/frequencies)@U.T;Spi=(U*frequencies)@U.T;Sp=(p/frequencies)@p.T
    require(max(eigh(SU,eigvals_only=True))<=8*R/Mbar*(1+1e-7),'charged soft trace bound')
    require(max(eigh(Spi,eigvals_only=True))<=8*R**3/Mbar*(1+1e-7),'charged momentum bound')
    om=frequencies[:,None]*frequencies[None,:];E=frequencies[:,None]+frequencies[None,:]
    norm=np.sqrt(8*om);r=eps/np.sqrt(1-eps*eps);Vb=a/eps**2
    B=vec.T@g['MT']@vec;C=vec.T@(g['KT']+g['V'])@vec
    J=(C-om*B)/norm
    P=r*(derivatives.T@derivatives-om*(trace.T@trace)-Vb*(U.T@U))/norm
    Pch=r*(p.T@p-om*(U.T@U)-Vb*(U.T@U))/norm
    neutral=P-Pch
    weighted=cutoff(E,R)/E
    AJ=float(np.sum(J*J*weighted));AP=float(np.sum(P*P*weighted))
    w=.95+.06j;kernel=2*w*w/(w*w-E*E)
    difference=np.sum(kernel*weighted*((J+P)**2-J**2))
    kmax=2*abs(w)**2/(w.imag*w.real)
    cs_bound=kmax*(AP+2*np.sqrt(AJ*AP))
    require(abs(difference)<=cs_bound*(1+1e-10),'weighted cross bound')
    whole=(vec.T@H@vec-om*(vec.T@M@vec))/norm
    whole_error=float(np.max(abs(whole)));require(whole_error<2e-6,('whole Hamiltonian pairs',whole_error))
    pair_rows=[]
    for cutoff_R in [.25,.5,1.,2.,R]:
        mask=E<=cutoff_R
        NP=float(np.sum(P[mask]**2));NJ=float(np.sum(J[mask]**2))
        Psplit=float(2*np.sum(Pch[mask]**2+neutral[mask]**2))
        require(NP<=Psplit+1e-12,'charged-neutral triangle')
        Lbar=2*g['A']+np.mean(D)
        count_bound=q*(1+Lbar*cutoff_R/np.pi)
        require(sum(frequencies<=cutoff_R)<=count_bound+1e-9,'min-max count')
        pair_rows.append(dict(R=cutoff_R,pressure_weight=NP,canonical_weight=NJ,
            pressure_bound_shape_ratio=NP/(cutoff_R**2*(1+q*eps**2*(1+cutoff_R**2)**4)),
            canonical_bound_ratio=NJ/(q*cutoff_R**2*(1+Lbar*cutoff_R)),triangle_bound=Psplit))
    # Direct finite element Green identity with the discrete Poisson columns.
    nt=g['nt'];uch=sum(vec[off:off+nt] for off in g['offsets'][:-1])/np.sqrt(q)
    end=g['ends'];inside=np.arange(1,nt-1)
    scalar_mass=g['MT'][:nt,:nt]
    # The test Gamma has endpoint data, and residual zero in the interior.
    from_green=g['dtn']@U-(g['poisson'].T@scalar_mass@uch)*(eig+1)
    green_error=relative(p,from_green);require(green_error<2e-7,('Green identity',green_error))
    row=dict(q=q,D=list(D),a=a,epsilon=eps,Mbar=Mbar,half_cells=n,dofs=len(M),
      positive_modes=len(eig),neutral_zero_error=float(zero_error),first_frequency=float(frequencies[0]),
      charged_capacity_eigenvalues=cap.tolist(),capacity_bound=8/Mbar,
      charged_soft_trace_max=float(max(eigh(SU,eigvals_only=True))),
      charged_derivative_scale_ratio=float(max(eigh(Sp,eigvals_only=True))/(R*(Mbar+(1+R*R)**2))),
      dtn=g['dtn'].tolist(),whole_H_pair_max=whole_error,green_identity_error=green_error,
      filtered_pressure_A=AP,filtered_canonical_A=AJ,
      filtered_response_difference=[float(difference.real),float(difference.imag)],
      cross_bound=float(cs_bound),pair_bands=pair_rows)
    return row

def fundamental(a,eps,k,parity,rtol=2e-11):
    A=np.arccos(eps)
    k2=complex(k)**2;require(k2.imag==0,'real spectral parameter');k2=k2.real
    y0=[1.,0.] if parity==0 else [0.,1.]
    result=solve_ivp(lambda z,y:[y[1],(a/np.cos(z)**2-k2)*y[0]],
       [0,A],y0,method='DOP853',rtol=rtol,atol=rtol*.01,dense_output=True)
    require(result.success,result.message)
    return result

def ode_dtn(a,eps):
    ends=[fundamental(a,eps,1j,p).y[:,-1] for p in [0,1]]
    even,odd=[float((y[1]/y[0]).real) for y in ends]
    return np.array([[even+odd,even-odd],[even-odd,even+odd]])/2

def first_even_frequency(a,eps,D):
    def residual(k):
        y= fundamental(a,eps,k,0).y[:,-1]
        return y[1]*np.cos(k*D/2)-k*y[0]*np.sin(k*D/2)
    grid=np.linspace(.01,1.9,30);values=[residual(k) for k in grid]
    bracket=next((grid[i],grid[i+1]) for i in range(len(grid)-1) if values[i]*values[i+1]<0)
    return brentq(residual,*bracket,xtol=2e-12)

def finite_elements():
    out=[]
    cases=[([.8],.05,1e-3),([.8],.05,1e-4),([.13,2.4],.03,2e-4),([.1,1.2,12.],.07,1e-3)]
    for D,a,eps in cases:
        refinements=[]
        meshes=[64,128,256] if eps==1e-4 else [32,64,128]
        for n in meshes:
            row=fem_case(D,a,eps,n);refinements.append(row)
            print('pressure FEM',len(D),eps,n,row['filtered_pressure_A'],flush=True)
        finest=refinements[-1];prev=refinements[-2]
        AP_change=abs(finest['filtered_pressure_A']-prev['filtered_pressure_A'])/finest['filtered_pressure_A']
        response_change=relative(np.array(finest['filtered_response_difference']),np.array(prev['filtered_response_difference']))
        require(AP_change<.025,('pressure refinement',D,eps,AP_change))
        require(response_change<.005,('response refinement',D,eps,response_change))
        independent=ode_dtn(a,eps)
        dtn_errors=[relative(np.array(r['dtn']),independent) for r in refinements]
        require(dtn_errors[-1]<3e-4,('ODE DtN',D,eps,dtn_errors))
        require(dtn_errors[-1]<.4*dtn_errors[-2],('DtN convergence',dtn_errors))
        extra={}
        if len(D)==1:
            exact=first_even_frequency(a,eps,D[0]);errors=[abs(r['first_frequency']-exact) for r in refinements]
            require(errors[-1]<3e-4 and errors[-1]<.4*errors[-2],('ODE eigenfrequency',errors))
            extra=dict(independent_frequency=exact,frequency_errors=errors)
        out.append(dict(q=len(D),D=D,a=a,epsilon=eps,refinements=refinements,
            finest_pressure_A_relative_change=AP_change,finest_response_normalized_change=response_change,
            independent_dtn=independent.tolist(),dtn_errors=dtn_errors,**extra))
    return out

def scaling():
    out=[]
    for logq in [32.,48.,64.,96.,128.]:
        q=np.exp(logq);delta=logq**-2;eps=q**-2;a=1/(4*np.pi**2*q)
        pressure=1/(np.sqrt(q)*delta**2);canonical=1/(q*delta**5)
        require(eps<=a and a<=1,'retained strong-barrier range')
        out.append(dict(log_q=logq,epsilon=eps,a=a,delta=delta,
             pressure_scale=pressure,canonical_scale=canonical,total_bound_shape=pressure+canonical))
    require(all(out[i+1]['total_bound_shape']<out[i]['total_bound_shape'] for i in range(len(out)-1)),'asymptotic decay')
    return out

def main():
    started=time.perf_counter();exact=exact_checks();fem=finite_elements();rates=scaling()
    sha=lambda p:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    files=['notes/interacting-mouth-pressure-band.md','notes/moving-mouth-clock-vertex.md','notes/singlet-complex-response.md']
    data=dict(date='2026-10-11',status='passed',script_sha256=sha(Path(__file__)),
        input_sha256={p:sha(ROOT/p) for p in files},python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        symbolic=exact,finite_elements=fem,scaling=rates,runtime_seconds=time.perf_counter()-started,
        fixed_band_interacting_pressure_controlled=True,renormalized_uv_response_controlled=False,
        interacting_local_contacts_complete=False,physical_damping_proved=False,novelty_established=False,
        lean_invocations=0,goal_complete=False)
    target=ROOT/'evidence/interacting-mouth-pressure-band.json';target.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status='passed',symbolic=len(exact),mesh_families=len(fem),
        mesh_runs=sum(len(x['refinements']) for x in fem),scalings=len(rates),runtime=data['runtime_seconds'])),flush=True)

if __name__=='__main__':main()
