"""Quadratic angular-current dipole action; not the full semiclassical response.
Run: python checks/rotational_current_response.py
Requires sympy, numpy, scipy. No Lean, cloud or external action is invoked.
"""
from pathlib import Path
import hashlib,json,platform,os,gc
os.environ['OPENBLAS_NUM_THREADS']='1'
import sympy as s
import numpy as np
import scipy
from scipy.sparse import diags,kron,eye,csc_matrix
from scipy.sparse.linalg import eigsh
ROOT=Path(__file__).resolve().parents[1]

def tensors(kind):
 t,r,c=s.symbols('t r c',real=True)
 V,W,Vr,Wr=[s.Function(n)(t,r) for n in ['V','W','Vr','Wr']]
 f=(r-1)**2/r**2 if kind=='RN' else 1+r*r
 b=r if kind=='RN' else s.S.One
 g=s.diag(-f,1/f,b*b/(1-c*c),b*b*(1-c*c)); inv=g.inv()
 h=s.zeros(4); h[0,3]=h[3,0]=-b*b*W*(1-c*c)
 h[1,3]=h[3,1]=-b*b*Wr*(1-c*c)
 di=-inv*h*inv
 def d(x,k): return s.diff(x,t) if k==0 else s.diff(x,r) if k==1 else s.diff(x,c) if k==2 else s.S.Zero
 G,D={},{}
 for i in range(4):
  for j in range(4):
   for k in range(4):
    G[i,j,k]=s.factor(sum(inv[i,l]*(d(g[l,k],j)+d(g[l,j],k)-d(g[j,k],l))/2 for l in range(4)))
    D[i,j,k]=s.factor(sum((di[i,l]*(d(g[l,k],j)+d(g[l,j],k)-d(g[j,k],l))+inv[i,l]*(d(h[l,k],j)+d(h[l,j],k)-d(h[j,k],l)))/2 for l in range(4)))
 checks={}
 Ric,dRic=s.zeros(4),s.zeros(4)
 for i in range(4):
  for j in range(i,4):
   base=delta=s.S.Zero
   for k in range(4):
    base+=d(G[k,i,j],k)-d(G[k,i,k],j)
    delta+=d(D[k,i,j],k)-d(D[k,i,k],j)
    for l in range(4):
     base+=G[k,i,j]*G[l,k,l]-G[l,i,k]*G[k,j,l]
     delta+=D[k,i,j]*G[l,k,l]+G[k,i,j]*D[l,k,l]-D[l,i,k]*G[k,j,l]-G[l,i,k]*D[k,j,l]
   Ric[i,j]=Ric[j,i]=s.factor(base); dRic[i,j]=dRic[j,i]=s.factor(delta)
 scalar=s.trace(inv*Ric); ds=s.factor(s.trace(di*Ric+inv*dRic))
 F,dF=s.zeros(4),s.zeros(4); F[2,3]=-1;F[3,2]=1
 dF[1,0]=(s.diff(V,r)-s.diff(Vr,t))*c;dF[0,1]=-dF[1,0];dF[2,0]=V;dF[0,2]=-V
 dF[2,1]=Vr;dF[1,2]=-Vr
 FU=inv*F*inv; dFU=di*F*inv+inv*dF*inv+inv*F*di
 square=sum(F[i,j]*FU[i,j] for i in range(4) for j in range(4))
 dsquare=sum(dF[i,j]*FU[i,j]+F[i,j]*dFU[i,j] for i in range(4) for j in range(4))
 stress=2*(F*inv*F.T-g*square/4)
 dstress=2*(dF*inv*F.T+F*di*F.T+F*inv*dF.T-h*square/4-g*dsquare/4)
 E=dRic-h*scalar/2-g*ds/2-dstress
 for i in range(4):
  for j in range(i,4):
   assert s.factor(Ric[i,j]-g[i,j]*scalar/2-stress[i,j])==0
   checks[f'background_E{i}{j}']='exact_zero'
 vol=s.trace(inv*h)/2
 raw={f'E{i}{j}':s.factor(E[i,j]) for i in range(4) for j in range(i,4)}
 for j in range(4):
  assert s.factor(sum(d(b*b*FU[i,j],i) for i in range(4)))==0
  checks[f'background_M{j}']='exact_zero'
  raw[f'M{j}']=s.factor(sum(d(b*b*(vol*FU[i,j]+dFU[i,j]),i) for i in range(4)))


 Ew=s.diff(Wr,t)-s.diff(W,r);Ev=s.diff(Vr,t)-s.diff(V,r)
 rw0=s.diff(b**4*Ew/6,r)-s.Rational(2,3)*(V-W)/f
 rw1=-s.diff(b**4*Ew/6,t)+s.Rational(2,3)*f*(Vr-Wr)
 rv0=s.diff(b**2*Ev/3,r)+s.Rational(2,3)*(V-W)/f
 rv1=-s.diff(b**2*Ev/3,t)-s.Rational(2,3)*f*(Vr-Wr)
 expected={'E03':-3*(1-c*c)*f*rw0/(b*b),'E13':3*(1-c*c)*rw1/(b*b*f),
           'M0':3*c*rv0,'M1':3*c*rv1,
           'M2':s.Rational(3,2)*(1-c*c)*(s.diff(rv0,t)+s.diff(rv1,r))}
 for key,expr in raw.items():
  assert s.factor(expr-expected.get(key,0))==0,(kind,key,s.factor(expr-expected.get(key,0)))
  checks['dynamic_action_'+key]='exact_zero'
 print(kind,'dynamic_tensor_checks',len(checks),flush=True)
 return checks

def action_normalization():
    r,c,a,v,R=s.symbols('r c a v R',real=True,nonzero=True)
    g=s.diag(-1,1,R**2/(1-c*c),R**2*(1-c*c))
    g[0,0]+=R**2*(1-c*c)*a*a*r*r
    g[0,3]=g[3,0]=-R**2*(1-c*c)*a*r
    inv=g.inv().applyfunc(s.factor)
    def d(x,k):return s.diff(x,r) if k==1 else s.diff(x,c) if k==2 else s.S.Zero
    conn={}
    for i in range(4):
      for j in range(4):
       for k in range(4):
        conn[i,j,k]=s.factor(sum(inv[i,l]*(d(g[l,k],j)+d(g[l,j],k)-d(g[j,k],l))/2 for l in range(4)))
    ric=s.zeros(4)
    for i in range(4):
      for j in range(4):
       ric[i,j]=s.factor(sum(d(conn[k,i,j],k)-d(conn[k,i,k],j)+sum(conn[k,i,j]*conn[l,k,l]-conn[l,i,k]*conn[k,j,l] for l in range(4)) for k in range(4)))
    scalar=s.factor(s.trace(inv*ric))
    assert s.factor(scalar-2/R**2-R**2*(1-c*c)*a*a/2)==0
    F=s.zeros(4);F[2,3]=-1;F[3,2]=1
    F[1,0]=v*c;F[0,1]=-v*c;F[2,0]=v*r;F[0,2]=-v*r
    up=inv*F*inv
    square=s.factor(sum(F[i,j]*up[i,j] for i in range(4) for j in range(4)))
    assert s.factor(square-(2/R**4-2*v*v*c*c-2*(v-a)**2*r*r*(1-c*c)/R**2))==0
    assert s.factor(g.det()+R**4)==0
    lag=s.integrate(R**2*(scalar-2/R**2-square+2/R**4)/8,(c,-1,1))
    assert s.factor(lag-(R**4*a*a/12+R**2*v*v/6+(v-a)**2*r*r/3))==0
    print('action_normalization',4,flush=True)
    return dict(count=4,curvature=str(scalar),action_density=str(s.factor(lag)))


def angular_checks():
    u=s.symbols('u',real=True)
    rows=[]
    for q in range(1,13):
      j=s.Rational(q-1,2); moments=[]
      for k in range(q):
        m=-j+k
        # Opposite monopole orientation to finite-flux-response-audit.md.
        rho=(1+u)**int(j-m)*(1-u)**int(j+m)
        norm=s.integrate(rho,(u,-1,1))
        avg=s.integrate(u*rho,(u,-1,1))/norm
        assert avg==-2*m/(q+1)
        metric=s.simplify(-m-s.Rational(q,2)*avg)
        charge=s.Rational(q,2)*avg
        assert metric==-m/(q+1) and charge==-q*m/(q+1)
        assert metric+charge==-m
        moments.append(m)
      T=sum(m*m for m in moments)
      assert T==s.Rational(q*(q*q-1),12) and sum(moments)==0
      # Independent spin matrices, including the off-axis response normalization.
      Jz=s.diag(*moments); Jplus=s.zeros(q)
      for k,m in enumerate(moments[:-1]):Jplus[k+1,k]=s.sqrt((j-m)*(j+m+1))
      Jx=(Jplus+Jplus.T)/2; Jy=(Jplus-Jplus.T)/(2*s.I)
      assert (Jx*Jy-Jy*Jx-s.I*Jz).applyfunc(s.simplify)==s.zeros(q)
      for A in [Jx,Jy,Jz]:
       for B in [Jx,Jy,Jz]:assert s.simplify(s.trace(A*B)-(T if A==B else 0))==0
      rows.append(dict(q=q,modes=q,T=str(T)))
    print('angular_projection',sum(x['modes'] for x in rows),flush=True)
    return rows


def reduction_checks():
    checks={}
    def zero(name,expr):
      assert s.simplify(expr)==0,(name,s.simplify(expr))
      checks[name]='exact_zero'
    Aw,Av,D,S,ga,kw,kv=s.symbols('Aw Av D S gamma kw kv',positive=True)
    Ew,Ev,pt,pr,ch,dc,xc,ww,vv,pw,pv,pc,dp,Q=s.symbols('Ew Ev Pt Pr chi dc xc ww vv pw pv pc dp Q',real=True)
    lag=Aw*Ew**2/2+Av*Ev**2/2+D*pt**2/(2*S)-D*S*pr**2/2+dc**2/(2*S)-S*xc**2/2+ga*ch*(kw*Ew+kv*Ev)
    zero('momentum_w',s.diff(lag,Ew)-(Aw*Ew+ga*kw*ch))
    zero('momentum_V',s.diff(lag,Ev)-(Av*Ev+ga*kv*ch))
    H=pw*(Ew+ww)+pv*(Ev+vv)+pc*dc-lag
    H=H.subs({Ew:(pw-ga*kw*ch)/Aw,Ev:(pv-ga*kv*ch)/Av,dc:S*pc})
    # Integration by parts: the common multiplier imposes (pw+pv)'=0.
    H=s.expand(H-pw*ww-pv*vv-dp*pt)
    zero('eliminate_temporal_difference',s.diff(H,pt).subs(pt,-S*dp/D))
    reduced=H.subs({pt:-S*dp/D,pw:Q-pv})
    target=(Q-pv-ga*kw*ch)**2/(2*Aw)+(pv-ga*kv*ch)**2/(2*Av)+S*(pc**2+xc**2)/2+S*dp**2/(2*D)+D*S*pr**2/2
    zero('Legendre_reduced_energy',reduced-target)
    y=s.symbols('y',real=True)
    potential=target.subs({Q:0,pv:s.sqrt(D)*y,pc:0,xc:0,dp:0,pr:0})
    mat=s.hessian(potential,[ch,y])
    zero('positive_determinant',mat.det().subs(kv,1-kw)-D*ga**2/(Aw*Av))
    q,K,R,h=s.symbols('q K R h',positive=True)
    nq=s.symbols('nq',integer=True,positive=True);ik=s.symbols('ik',integer=True)
    zero('general_spin_trace',s.summation((ik-(nq-1)/2)**2,(ik,0,nq-1))-nq*(nq*nq-1)/12)
    explicit=s.Matrix([[3*h*(2/R**4+q*q/R**2)/(q+1)**2,s.sqrt(6*h)*(2/R**4-q/R**2)/(q+1)],
                       [s.sqrt(6*h)*(2/R**4-q/R**2)/(q+1),4/R**4+2/R**2]])
    reduced_mat=mat.subs({Aw:K*R**4/6,Av:K*R**2/3,D:2*K/3,ga:s.sqrt(h*K),kw:1/(q+1),kv:q/(q+1)})
    for i in range(2):
     for j in range(2):zero(f'variable_radius_matrix_{i}{j}',reduced_mat[i,j]-explicit[i,j])
    zero('variable_radius_determinant',explicit.det()-12*h/R**6)
    c=(q-2)/(3*(q+1));B=6+(2+9*c*c)*h
    zero('product_trace',s.trace(explicit.subs(R,1))-B)
    X,ss=s.symbols('x s0',real=True)
    poly=(X*s.eye(2)-explicit.subs(R,1)).det()
    zero('product_polynomial',poly-(X*X-B*X+12*h))
    hp=ss*(1-1/q**2)/3
    zero('finite_flux_trace',B.subs(h,hp)-(6+ss*(q-1)*(q*q+2)/(q*q*(q+1))))
    zero('finite_flux_determinant',(12*h).subs(h,hp)-4*ss*(1-1/q**2))
    zero('published_large_flux_polynomial',s.limit(poly.subs(h,hp),q,s.oo)-(X*X-(6+ss)*X+4*ss))
    # General reduced local energy identity, with nonconstant S and M.
    t,r=s.symbols('t r',real=True);Sf=s.Function('S')(r)
    z=s.Matrix([s.Function('z1')(t,r),s.Function('z2')(t,r)])
    M=s.Matrix([[s.Function('m11')(r),s.Function('m12')(r)],[s.Function('m12')(r),s.Function('m22')(r)]])
    energy=((z.diff(t).T*z.diff(t))[0]/Sf+Sf*(z.diff(r).T*z.diff(r))[0]+(z.T*M*z)[0])/2
    flux=(z.diff(t).T*Sf*z.diff(r))[0]
    evo=s.diff(Sf*z.diff(r),r)-M*z
    residual=s.diff(energy,t)-s.diff(flux,r)
    residual=residual.subs({s.diff(z[i],t,2):Sf*evo[i] for i in range(2)})
    zero('local_energy_flux',residual)
    # Positive Hamiltonian equations reconstruct both unreduced vector equations.
    W,V,Wr,Vr,chi=[s.Function(n)(t,r) for n in ['W','V','Wr','Vr','chi']]
    A=s.Function('Aw')(r);Bf=s.Function('Av')(r)
    EW=s.diff(Wr,t)-s.diff(W,r);EV=s.diff(Vr,t)-s.diff(V,r)
    Pw=A*EW+ga*kw*chi;PV=Bf*EV+ga*kv*chi
    rw0=s.diff(Pw,r)-D*(V-W)/Sf;rv0=s.diff(PV,r)+D*(V-W)/Sf
    rw1=-s.diff(Pw,t)+D*Sf*(Vr-Wr);rv1=-s.diff(PV,t)-D*Sf*(Vr-Wr)
    zero('Gauss_spatial',rw0+rv0-s.diff(Pw+PV,r))
    zero('Gauss_temporal',rw1+rv1+s.diff(Pw+PV,t))
    zero('common_rotation_Ward',s.diff(rw0+rv0,t)+s.diff(rw1+rv1,r))
    yf=s.Function('y')(t,r);pf=s.sqrt(D)*yf
    Phit=-Sf*s.diff(pf,r)/D;Phir=-s.diff(pf,t)/(D*Sf)
    ew=(Q-pf-ga*kw*chi)/A;ev=(pf-ga*kv*chi)/Bf
    for name,expr in {
      'w_time_constraint':s.diff(Q-pf,r)-D*Phit/Sf,
      'V_time_constraint':s.diff(pf,r)+D*Phit/Sf,
      'w_radial_equation':-s.diff(Q-pf,t)+D*Sf*Phir,
      'V_radial_equation':-s.diff(pf,t)-D*Sf*Phir}.items():zero(name,expr)
    Uy=(pf-ga*kv*chi)/Bf*s.sqrt(D)-(Q-pf-ga*kw*chi)/A*s.sqrt(D)
    Uch=-ga*(kw*ew+kv*ev)
    EY=s.diff(yf,t,2)/Sf-s.diff(Sf*s.diff(yf,r),r)+Uy
    zero('relative_curl_reconstruction',ev-ew-s.diff(Phir,t)+s.diff(Phit,r)-EY/s.sqrt(D))
    zero('scalar_anomaly_reconstruction',Uch+ga*(kw*ew+kv*ev))
    print('reduction_identities',len(checks),flush=True)
    return checks


def spectral_checks():
    rows=[]
    for q,h in [(2,.2),(3,.7),(10,4.)]:
      c=(q-2)/(3*(q+1));B=6+(2+9*c*c)*h
      hi=(B+np.sqrt(B*B-48*h))/2;lo=12*h/hi
      exact=np.sort(np.array([(n+(1+np.sqrt(1+4*m))/2)**2 for m in [lo,hi] for n in range(6)]))[:6]
      M=np.array([[(2+9*c*c)*h,-c*np.sqrt(54*h)],[-c*np.sqrt(54*h),6.]])
      errors=[];samples=[]
      for n in [255,511,1023,2047]:
        dx=np.pi/(n+1);x=-np.pi/2+dx*np.arange(1,n+1)
        lap=diags([-np.ones(n-1),2*np.ones(n),-np.ones(n-1)],[-1,0,1])/dx**2
        op=kron(lap,eye(2))+kron(diags(1/np.cos(x)**2),csc_matrix(M))
        vals=np.sort(eigsh(op,k=6,sigma=0,which='LM',tol=2e-11,return_eigenvectors=False))
        err=float(np.max(np.abs(vals-exact)/np.maximum(1,exact)))
        errors.append(err);samples.append(dict(points=n,eigenvalues=vals.tolist(),relative_error=err))
      assert all(b<a*.85 for a,b in zip(errors,errors[1:])),errors
      assert errors[-1]<2e-4,errors
      rows.append(dict(q=q,h=h,masses_squared=[lo,hi],exact_eigenvalues=exact.tolist(),grids=samples))
    print('coupled_spectrum',len(rows),'max_error',max(x['grids'][-1]['relative_error'] for x in rows),flush=True)
    return rows


def main():
    retained=json.loads((ROOT/'evidence/lll-axial-vertex.json').read_text())
    assert hashlib.sha256((ROOT/'checks/lll_axial_vertex.py').read_text().encode()).hexdigest()==retained['script_sha256']
    result=dict(status='passed_projected_rotational_current_not_full_semiclassical_response',
                action=action_normalization(),tensors={kind:tensors(kind) for kind in ['throat','RN']},
                angular=angular_checks(),reduction=reduction_checks())
    result['retained_Dirac_source_matches_receipt']=True
    s.core.cache.clear_cache();gc.collect()
    result['spectra']=spectral_checks()
    result['script_sha256']=hashlib.sha256(Path(__file__).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest()
    result['input_sha256']={name:hashlib.sha256((ROOT/name).read_text(encoding='utf-8').replace('\r\n','\n').encode()).hexdigest() for name in [
      'notes/rotational-current-response.md','notes/stationary-electric-dipole.md',
      'notes/finite-flux-response-audit.md','notes/casimir-cylinder-ward-completion.md','notes/weak-coupling-window.md',
      'checks/lll_axial_vertex.py','evidence/lll-axial-vertex.json']}
    result['versions']=dict(python=platform.python_version(),sympy=s.__version__,numpy=np.__version__,scipy=scipy.__version__)
    result['limitations']=['Free LLL projected quadratic current action only','Physical boundary spin/holonomy and state not selected',
      'Higher-level/contact and background Casimir-state Hessian not controlled','No full wormhole stability, transport or shortcut theorem',
      'Rotational screening is established prior work; consequential novelty not established']
    result['lean_invocations']=0
    (ROOT/'evidence/rotational-current-response.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')


if __name__=='__main__':main()
