"""Finite-frequency RN dipole reconstruction and prescribed classical matching.

No binary exterior, complete quantum-open background or state is constructed.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import scipy
from scipy.integrate import solve_ivp
import sympy as s

ROOT=Path(__file__).resolve().parents[1]

def text_hash(path):
    return hashlib.sha256(path.read_text(encoding='utf-8').encode()).hexdigest()

def geometry_equations():
    r,c,w=s.symbols('r c w',real=True)
    T,R,J,a=[s.Function(x)(r) for x in ['T','R','J','a']]
    jets=[s.diff(v,r,k) for v in [T,R,J,a] for k in range(3)]
    def linear_factor(expr):
        poly=s.Poly(s.expand(expr),*jets)
        out=0
        for powers,coeff in poly.terms():
            assert sum(powers)<=1
            out+=s.factor(coeff)*s.prod(v**k for v,k in zip(jets,powers))
        return out
    f=(r-1)**2/r**2
    g=s.diag(-f,1/f,r*r/(1-c*c),r*r*(1-c*c))
    inv=g.inv()
    h=s.zeros(4)
    h[0,0]=-f*T*c; h[1,1]=R*c/f; h[0,1]=h[1,0]=J*c
    di=-inv*h*inv
    def d(v,k,linear=False):
        if k==0: return -s.I*w*v if linear else s.S.Zero
        if k==1: return s.diff(v,r)
        if k==2: return s.diff(v,c)
        return s.S.Zero
    G,D={},{ }
    for i in range(4):
     for j in range(4):
      for k in range(4):
       G[i,j,k]=s.factor(sum(inv[i,l]*(d(g[l,k],j)+d(g[l,j],k)-d(g[j,k],l))/2 for l in range(4)))
       D[i,j,k]=linear_factor(sum((inv[i,l]*(d(h[l,k],j,True)+d(h[l,j],k,True)-d(h[j,k],l,True))+di[i,l]*(d(g[l,k],j)+d(g[l,j],k)-d(g[j,k],l)))/2 for l in range(4)))
    print('connections_derived',flush=True)
    Ric,dRic=s.zeros(4),s.zeros(4)
    for i in range(4):
     for j in range(i,4):
      b=v=s.S.Zero
      for k in range(4):
       b+=d(G[k,i,j],k)-d(G[k,i,k],j)
       v+=d(D[k,i,j],k,True)-d(D[k,i,k],j,True)
       for l in range(4):
        b+=G[k,i,j]*G[l,k,l]-G[l,i,k]*G[k,j,l]
        v+=D[k,i,j]*G[l,k,l]+G[k,i,j]*D[l,k,l]-D[l,i,k]*G[k,j,l]-G[l,i,k]*D[k,j,l]
      Ric[i,j]=Ric[j,i]=s.factor(b)
      dRic[i,j]=dRic[j,i]=linear_factor(v)
      print('Ricci',i,j,flush=True)
    scalar=s.trace(inv*Ric)
    ds=linear_factor(s.trace(di*Ric+inv*dRic))
    F,dF=s.zeros(4),s.zeros(4)
    F[2,3]=-1; F[3,2]=1
    for i in range(3): dF[i,3]=d(a*(1-c*c),i,True); dF[3,i]=-dF[i,3]
    FU=inv*F*inv
    dFU=di*F*inv+inv*dF*inv+inv*F*di
    square=sum(F[i,j]*FU[i,j] for i in range(4) for j in range(4))
    dsquare=sum(dF[i,j]*FU[i,j]+F[i,j]*dFU[i,j] for i in range(4) for j in range(4))
    dStress=2*(dF*inv*F.T+F*di*F.T+F*inv*dF.T-h*square/4-g*dsquare/4)
    E=dRic-h*scalar/2-g*ds/2-dStress
    out={}
    for i in range(4):
     for j in range(i,4):
      expr=linear_factor(E[i,j]); out[f'E{i}{j}']=expr
    vol=s.trace(inv*h)/2
    for j in range(4): out[f'M{j}']=linear_factor(sum(d(r*r*(vol*FU[i,j]+dFU[i,j]),i,True) for i in range(4)))
    out['scalar']=ds
    background={}
    stress=2*(F*inv*F.T-g*square/4)
    for i in range(4):
        for j in range(i,4):
            assert s.factor(Ric[i,j]-g[i,j]*scalar/2-stress[i,j])==0
            background[f'background_Einstein_Maxwell_{i}_{j}']='exact_zero'
    for j in range(4):
        assert s.factor(sum(d(r*r*FU[i,j],i) for i in range(4)))==0
        background[f'background_Maxwell_{j}']='exact_zero'
    return out,background,(r,c,w,T,R,J,a)

def physical_reconstruction(raw,symbols):
    r,c,w,T,R,J,a=symbols
    A,Z,V=s.symbols('A Z V')
    z=r-1; v=3*r-2; f=z*z/r**2
    W=3*z*z/(r*(v/r)**2)+r/v
    pot=f/r**3*(2-2/r+6*W)
    ap=V+(r/z**2+2/(r*v))*Z-r*A/z**2
    zp=V
    vp=-s.diff(f,r)/f*V+(pot-w*w)/f**2*Z
    def simp(expr):
        return sum(s.factor(s.diff(expr,state))*state for state in [A,Z,V])
    def d(expr):
        return simp(s.diff(expr,r)+s.diff(expr,A)*ap+s.diff(expr,Z)*zp+s.diff(expr,V)*vp)
    Rsol=2*v*(A-Z)/z**2
    Jsol=s.I*w*r*Rsol
    Tsol=(w*w*r**4-z)*Rsol/v+4*z*z/(r*v)*V+(4/v+8*z*z/(r*r*v*v))*Z
    # Before imposing the master equation, Maxwell gives exactly that equation.
    V2=s.symbols('V2')
    a2=s.diff(ap,r)+s.diff(ap,A)*ap+s.diff(ap,Z)*V+s.diff(ap,V)*V2
    independent=raw['M3'].xreplace({a:A,s.diff(a,r):ap,s.diff(a,r,2):a2,R:Rsol,T:Tsol})
    expected=f*V2+s.diff(f,r)*V+(w*w-pot)*Z/f
    for state in [A,Z,V,V2]:
        assert s.factor(s.diff(independent-expected,state))==0
    assert (independent-expected).subs({A:0,Z:0,V:0,V2:0})==0
    subs={}
    for field,val in [(a,A),(R,Rsol),(J,Jsol),(T,Tsol)]:
        for k in range(3):
            subs[s.diff(field,r,k)]=val
            if k<2: val=d(val)
        print('reconstructed',field,flush=True)
    checks={}
    for name,expr in raw.items():
        replaced=expr.xreplace(subs)
        assert replaced.subs({A:0,Z:0,V:0})==0
        result=simp(replaced)
        assert result==0,(name,str(result))
        checks['linear_'+name]='exact_zero'
    print('finite_frequency_reconstruction_passed',flush=True)
    return checks

def bound(lo,hi):
    return ((lo**-2-hi**-2)/6+2.5*(lo**-1-hi**-1)+6.2*np.log(hi/lo)
            +5.7*(hi-lo)+.9*(hi*hi-lo*lo))

def load_bound(lo,hi,nu):
    primitive=lambda z:-1/(5*z**5)-1.5/z**4-3/z**3
    return abs(nu)*lo**5/(1+3*lo)**2*(primitive(hi)-primitive(lo))

def transfer_checks():
    cases=[]
    for lo,hi in [(.003,1.),(.003,100.),(.03,1.),(.03,100.)]:
     for eta in [.003,.03,.15]:
      for nu in [0.,.02,-.02]:
       omega=np.sqrt(eta/bound(lo,hi))
       y=np.linspace(np.log(lo),np.log(hi),161)
       def rhs(y,w):
           z=np.exp(y)
           return [w[1],-(5+9*z)/(1+3*z)*w[1]-omega**2*(1+z)**4/z**2*w[0]]
       results=[]
       for method,tol in [('DOP853',2e-11),('Radau',2e-10)]:
           sol=solve_ivp(rhs,(y[0],y[-1]),[1,nu],method=method,t_eval=y,rtol=tol,atol=tol/100)
           assert sol.success,sol.message
           results.append(sol.y[0])
       error=max(abs(results[0]-results[1]))
       zeta=load_bound(lo,hi,nu)
       boundW=(zeta+eta)/(1-eta)
       actual=max(abs(results[0]-1))
       normed=results[0]/results[0][-1]
       normerror=max(abs(normed-1))
       normbound=2*boundW/(1-boundW)
       assert error<2e-8 and actual<=boundW and normerror<=normbound
       if nu==0:
           assert min(results[0])>=1-eta-1e-10 and max(results[0])<=1+1e-10
           assert normerror<=eta/(1-eta)
       cases.append(dict(z_inner=lo,z_outer=hi,omega=omega,eta=eta,nu=nu,zeta=zeta,
                         solver_discrepancy=float(error),ratio_error=float(actual),ratio_bound=float(boundW),
                         outer_normalized_error=float(normerror),outer_normalized_bound=float(normbound)))
    return cases

def throat(R,sigma):
    x=np.arctan(R); f=1+R*R; K=1-sigma*sigma+3*R*R
    sn=np.sin(sigma*x); cs=np.cos(sigma*x); sx=x*np.sinc(sigma*x/np.pi)
    E=K*cs-3*sigma*R*sn
    O=K*sx+3*R*cs
    Ep=(6*R-3*sigma*sigma*R/f)*cs-sigma*(K/f+3)*sn
    Op=6*R*sx+(K/f+3)*cs-3*sigma*R*sn/f
    even=E/(6*np.cos(sigma*np.pi/2)); odd=O/(3*np.pi*np.sinc(sigma/2))
    evenp=Ep/(6*np.cos(sigma*np.pi/2)); oddp=Op/(3*np.pi*np.sinc(sigma/2))
    return np.array([even+odd,even-odd]),np.array([evenp+oddp,evenp-oddp])

def matched_checks():
    cases=[]
    for lam in [.003,.001,.0003,.0001]:
        R=lam**(-1/3); lo=lam*R; hi=3.; z0=(1+lo)*lo**2/(1+3*lo)
        slope0=2-2*lo/((1+lo)*(1+3*lo)); fluxratio=(1+lo)**2*(1+1/R**2)
        static=None
        for sigma in [0.,.25,.5]:
            v,d=throat(R,sigma)
            vals=lam**2*v/z0
            dots=lam**2*(fluxratio*R*d-slope0*v)/z0
            initial=np.r_[vals,dots]
            omega=lam*sigma
            def rhs(y,s):
                z=np.exp(y)
                return np.r_[s[2:],-(5+9*z)/(1+3*z)*s[2:]-omega**2*(1+z)**4/z**2*s[:2]]
            matrices=[]
            for method,tol in [('DOP853',2e-11),('Radau',2e-10)]:
                sol=solve_ivp(rhs,(np.log(lo),np.log(hi)),initial,method=method,rtol=tol,atol=tol/100)
                assert sol.success
                row=sol.y[:2,-1]
                matrices.append(np.array([row,row[::-1]]))
            M=matrices[0]
            if sigma==0: static=M.copy()
            delta=float(np.linalg.norm(M-np.eye(2),np.inf))
            dyn=float(np.linalg.norm(M-static,np.inf))
            discrepancy=float(np.max(abs(matrices[0]-matrices[1])))
            eta=omega**2*bound(lo,hi)
            gap_squared=min(6.,2/(3*lam*lam)*min(lo*lo/(1+lo)**4,hi*hi/(1+hi)**4))
            assert discrepancy<1e-8 and delta<.5
            assert gap_squared>sigma*sigma
            cases.append(dict(lam=lam,R=R,z_inner=lo,z_outer=hi,sigma=sigma,eta=eta,
                              matrix=M.tolist(),matrix_error=delta,dynamic_matrix_error=dyn,
                              geometric_error_ratio=delta/(lo+R**-2+eta),
                              dynamic_error_ratio=(dyn/(sigma*sigma/R**2+eta) if sigma else None),
                              solver_discrepancy=discrepancy,
                              conservative_gap_squared=float(gap_squared),
                              minimum_singular_value=float(np.min(np.linalg.svd(M,compute_uv=False)))))
    return cases

def algebra_checks():
    z,w=s.symbols('z w',positive=True)
    r=1+z; v=1+3*z; f=z*z/(r*r)
    W=3*z*z/(r*(v/r)**2)+r/v
    potential=f/r**3*(2-2/r+6*W)
    Z0=r*z*z/v
    assert s.factor(f*s.diff(f*s.diff(Z0,z),z)-potential*Z0)==0
    numerator=9*r**4-18*r*r+16*r-4
    assert s.factor(potential-2*z*z*numerator/(r**6*v*v))==0
    assert s.expand(numerator-3*r**4)==6*z**4+24*z**3+18*z*z+4*z
    p=z**6/v**2; weight=z*z*r**4/v**2
    assert s.factor(f*Z0**2-p)==0
    assert s.factor(Z0**2/f-weight)==0
    assert s.factor(z*s.diff(p,z)/p-1-(5+9*z)/v)==0
    assert s.factor(z*z*weight/p-r**4/z**2)==0
    integral_upper=z**3/s.Integer(3)+z**4/s.Integer(2)+z**5/s.Integer(5)
    primitive=9*z*z/10+57*z/10+s.Rational(31,5)*s.log(z)-s.Rational(5,2)/z-1/(6*z*z)
    assert s.factor(s.diff(primitive,z)-integral_upper/p)==0
    J=-1/(5*z**5)-s.Rational(3,2)/z**4-3/z**3
    assert s.factor(s.diff(J,z)-1/p)==0
    I=s.symbols('I')
    amouth=z*z*(s.Rational(1,2)+I)
    Rmouth=2*(v*I-1)
    assert s.factor(amouth-z*z*Rmouth/(2*v)-s.Rational(3,2)*Z0)==0
    assert s.limit(Z0/(z*z),z,0)==1
    chi=s.symbols('chi')
    chip=(1-1/f)*chi/r
    agauge=chi/r
    Rgauge=chi*s.diff(f,z)/f-2*chip
    assert s.factor(agauge-z*z*Rgauge/(2*v))==0
    assert s.factor(s.diff(Z0,z)*z/Z0-(2-2*z/(r*v)))==0
    return dict(static_master='r(r-1)^2/(3r-2)',static_coefficient='3/2',
                residual_gauge_annihilated=True,Volterra_primitives_checked=True,
                static_wave_residual='exact_zero',factorization_checked=True,
                positive_potential_lower_bound_checked=True)

def scale_checks():
    cases=[]
    for q in [10.,100.,1000.,1000000.]:
        lam=1/(16*np.pi*q**3)
        sigma=np.log(q)**3/(2*np.sqrt(np.pi)*q**1.5)
        lo=lam**(2/3)
        hi=16*np.pi*q**1.5/np.log(q)**2-1
        eta=(sigma*lam)**2*bound(lo,hi)
        assert hi>lo>0 and 0<eta<1
        cases.append(dict(q=q,lam=lam,sigma=sigma,z_inner=lo,z_outer=hi,
                          eta=float(eta),throat_end_dynamic_error_scale=float(sigma*sigma*lam**(2/3))))
    assert all(a['eta']>b['eta'] for a,b in zip(cases,cases[1:]))
    return cases

def verify():
    raw,background,symbols=geometry_equations()
    reconstruction=physical_reconstruction(raw,symbols)
    return dict(status='passed_dynamic_RN_dipole_reconstruction_and_prescribed_matching',
                tensor=dict(checks={**background,**reconstruction},count=len(background)+len(reconstruction)),
                algebra=algebra_checks(),transfer_cases=transfer_checks(),
                matched_cases=matched_checks(),scale_cases=scale_checks(),
                physical_RN_master_and_reconstruction_verified=True,
                Maxwell_equation_reduces_to_master_before_substitution=True,
                prescribed_piecewise_classical_master_matching_verified=True,
                actual_binary_outer_drive_derived=False,complete_four_dimensional_junction_proved=False,
                boundary_frames_and_holonomy_fixed=False,semiclassical_source_response_derived=False,
                time_dependent_quantum_state_derived=False,physical_stability_or_transport_proved=False,
                novelty_established=False,lean_invocations=0,
                versions=dict(python=sys.version.split()[0],numpy=np.__version__,scipy=scipy.__version__,sympy=s.__version__))

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    result=verify()
    result['script_sha256']=text_hash(Path(__file__))
    result['input_sha256']={name:text_hash(ROOT/name) for name in [
        'notes/dynamic-rn-dipole.md','notes/dipole-throat-matching.md','notes/magnetic-mouth-force.md']}
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],tensor_checks=result['tensor']['count'],
                         transfer_cases=len(result['transfer_cases']),matched_cases=len(result['matched_cases']))))
