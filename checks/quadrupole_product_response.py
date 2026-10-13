"""Fixed-product polar Einstein-Maxwell response; not the driven quantum state."""
from pathlib import Path
import argparse
import hashlib
import json
import math

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.special import eval_gegenbauer, gammaln, roots_jacobi
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def symbolic_tensor_check():
    # AdS2 radius=1, sphere radius squared=C, kappa^2=2.
    r,c,w,C,e=s.symbols('r c w C e',real=True,nonzero=True)
    T,R,J,F,u,v,W=[s.Function(x)(r) for x in ['T','R','J','F','u','v','W']]
    fields=[T,R,J,F,u,v,W]
    jets=[s.diff(x,r,n) for x in fields for n in range(5)]
    B=1+r*r; Y=s.legendre(2,c); L=s.Integer(6)
    def back(expr):
        expr=s.cancel(expr)
        num,den=s.fraction(expr)
        rem=s.rem(s.Poly(num,e),s.Poly(e*e-(C+1)/(2*C),e)).as_expr()
        return s.factor(rem/den)
    def lin(expr):
        poly=s.Poly(s.expand(expr),*jets)
        out=0
        for powers,coeff in poly.terms():
            assert sum(powers)<=1
            out+=back(coeff)*s.prod(z**k for z,k in zip(jets,powers))
        return out
    def d(expr,k,pert=False):
        if k==0: return -s.I*w*expr if pert else s.S.Zero
        if k==1: return s.diff(expr,r)
        if k==2: return s.diff(expr,c)
        return s.S.Zero
    g=s.diag(-B,1/B,C/(1-c*c),C*(1-c*c)); gi=g.inv()
    h=s.zeros(4)
    h[0,0]=-B*T*Y; h[1,1]=R*Y/B; h[0,1]=h[1,0]=J*Y
    h[2,2]=2*C*F*Y/(1-c*c); h[3,3]=2*C*F*Y*(1-c*c)
    dgi=-gi*h*gi
    gam,dgam={},{}
    for i in range(4):
     for j in range(4):
      for k in range(4):
       gam[i,j,k]=s.factor(sum(gi[i,z]*(d(g[z,k],j)+d(g[z,j],k)-d(g[j,k],z))/2 for z in range(4)))
       dgam[i,j,k]=lin(sum((gi[i,z]*(d(h[z,k],j,True)+d(h[z,j],k,True)-d(h[j,k],z,True))
                           +dgi[i,z]*(d(g[z,k],j)+d(g[z,j],k)-d(g[j,k],z)))/2 for z in range(4)))
    print('tensor connections derived',flush=True)
    ric,dric=s.zeros(4),s.zeros(4)
    for i in range(4):
     for j in range(i,4):
      base=pert=s.S.Zero
      for k in range(4):
       base+=d(gam[k,i,j],k)-d(gam[k,i,k],j)
       pert+=d(dgam[k,i,j],k,True)-d(dgam[k,i,k],j,True)
       for z in range(4):
        base+=gam[k,i,j]*gam[z,k,z]-gam[z,i,k]*gam[k,j,z]
        pert+=dgam[k,i,j]*gam[z,k,z]+gam[k,i,j]*dgam[z,k,z]-dgam[z,i,k]*gam[k,j,z]-gam[z,i,k]*dgam[k,j,z]
      ric[i,j]=ric[j,i]=s.factor(base)
      dric[i,j]=dric[j,i]=lin(pert)
    scalar=s.trace(gi*ric); dscalar=lin(s.trace(dgi*ric+gi*dric))
    em,dem=s.zeros(4),s.zeros(4)
    em[0,1]=e; em[1,0]=-e
    potential=[u*Y,v*Y,0,0]
    for i in range(4):
     for j in range(4): dem[i,j]=d(potential[j],i,True)-d(potential[i],j,True)
    # Independent radial electric component permits a magnetic-current source
    # after dualizing the original magnetic background. Do not impose dF=0.
    dem[0,1]=W*Y; dem[1,0]=-W*Y
    magnetic=s.I*w*v+s.diff(u,r)+W
    raised=gi*em*gi
    draised=dgi*em*gi+gi*dem*gi+gi*em*dgi
    square=sum(em[i,j]*raised[i,j] for i in range(4) for j in range(4))
    dsquare=sum(dem[i,j]*raised[i,j]+em[i,j]*draised[i,j] for i in range(4) for j in range(4))
    stress=2*(em*gi*em.T-g*square/4)
    dstress=2*(dem*gi*em.T+em*dgi*em.T+em*gi*dem.T-h*square/4-g*dsquare/4)
    lam=(1/C-1)/2
    assert (ric-g*scalar/2+lam*g-stress).applyfunc(back)==s.zeros(4)
    residual=(dric-h*scalar/2-g*dscalar/2+lam*h-dstress).applyfunc(lin)
    vol=s.trace(gi*h)/2
    maxwell=[lin(sum(d(vol*raised[i,j]+draised[i,j],i,True) for i in range(4))) for j in range(4)]
    print('four-dimensional residuals derived',flush=True)
    # The KI current uses div on the second antisymmetric index: J_KI=-maxwell.
    At=lin((-C*W+e*C*(T+R-4*F)/2)/L)
    jt=lin(-C*maxwell[1]/(L*Y)); jr=lin(C*maxwell[0]/(L*Y))
    ST=-(T+R)
    trace=lin((-residual[0,0]/B+B*residual[1,1])/Y)
    wave=lambda x: w*w*x/B+s.diff(B*s.diff(x,r),r)
    f_source=lin(wave(F)-(L/C+2)*F-2*L*e*At/C-(trace/2-L*ST/(4*C)))
    a_source=lin(wave(At)-L*At/C-2*e*F-(e*ST/2-s.I*w*jt/B-s.diff(B*jr,r)+magnetic))
    assert f_source==0,('metric source identity',str(f_source))
    assert a_source==0,('Maxwell source identity',str(a_source))
    print('covariant source identities exact',flush=True)
    # Source-free reconstruction with arbitrary frequency and radial master data.
    ff,fp,aa,ap=s.symbols('ff fp aa ap')
    state=[ff,fp,aa,ap]
    rhs=[fp,((L/C+2-w*w/B)*ff+2*L*e*aa/C-2*r*fp)/B,
         ap,((L/C-w*w/B)*aa+2*e*ff-2*r*ap)/B]
    def drv(expr):
        return s.diff(expr,r)+sum(s.diff(expr,x)*dx for x,dx in zip(state,rhs))
    xx=4*C/L*((w*w/B-L/C-1)*ff+r*fp)-4*e*aa
    tt=xx+2*ff
    zz=-s.I*w*4*C/L*(-r*ff+B*fp)
    solutions=[tt,-tt,zz/B,ff,B*ap,-s.I*w*aa/B,-L*aa/C-2*e*ff]
    subs={}
    for field,value in zip(fields,solutions):
        for order in range(3):
            subs[s.diff(field,r,order)]=value
            if order<2: value=drv(value)
    count=0
    for name,expr in [(f'E{i}{j}',residual[i,j]) for i in range(4) for j in range(i,4)]+[(f'M{i}',z) for i,z in enumerate(maxwell)]:
        val=s.expand(expr.xreplace(subs))
        assert val.subs(dict.fromkeys(state,0))==0
        for x in state: assert back(s.diff(val,x))==0,(name,x)
        count+=1
    print('all reconstructed tensor components exact',flush=True)
    return {'background_equations':'exact_zero','source_identities':2,
            'source_free_tensor_components':count,'frequency':'arbitrary, including zero',
            'electric_and_magnetic_current_sources':True,
            'sphere_index':2,'radius_ratio':'symbolic C>=1','method':'direct linearized four-dimensional metric, Ricci, Maxwell and stress tensors'}


def masses(L,C):
    y=(C+1)/2
    R=math.sqrt(C*C+4*L*y)
    small=L*(L-2)/(C*(C+L+R))
    large=(C+L+R)/C
    return small,large


def exterior_limit_check():
    # Independent RN-AdS axial/polar potentials from Berti-Kokkotas (13)-(17).
    xx,nn=s.symbols('xx nn',nonnegative=True)
    mass=1+2*xx; q2=1+3*xx; C=1+6*xx
    delta=s.sqrt(9*mass**2+8*nn*q2)
    p1,p2=3*mass+delta,3*mass-delta
    varpi=nn+3*mass-2*q2
    W=(nn+mass-2*xx)/varpi
    U=(2*nn+3*mass)*W+varpi-nn-mass+2*xx
    assert s.simplify(W-1)==0
    for sign in [-1,1]:
        target=(C+2*(nn+1)+sign*delta)/C
        polar=(U+sign*delta*W)/C
        axial=(2*(nn+1)-3*mass+sign*delta+4*q2)/C
        assert s.simplify(polar-target)==0
        assert s.simplify(axial-target)==0
    assert s.simplify(delta**2-(C*C+4*2*(nn+1)*q2))==0
    return {'exact_extremal_horizon_limits':4,'mass_discriminant':'exact','no_global_isospectrality_claim':True}


def mode_checks():
    rows=[]
    for C in [1.,4.,61.,6001.]:
      for mu in masses(6,C):
        d=(1+math.sqrt(1+4*mu))/2
        x,weight=roots_jacobi(80,d-.5,d-.5)
        modes=[]
        for n in range(5):
            log_h=math.log(math.pi)+(1-2*d)*math.log(2)+gammaln(n+2*d)-gammaln(n+1)-math.log(n+d)-2*gammaln(d)
            norm=math.exp(-log_h/2)
            modes.append(norm*eval_gegenbauer(n,d,x))
        mat=np.array(modes)
        gram=(mat*weight)@mat.T
        gram_error=float(np.max(abs(gram-np.eye(5))))
        z=np.linspace(-1.31,1.31,101); sn=np.sin(z); cs=np.cos(z)
        residual=0.
        for n in range(5):
            poly=eval_gegenbauer(n,d,sn)
            dp=2*d*eval_gegenbauer(n-1,d+1,sn) if n else np.zeros_like(z)
            ddp=4*d*(d+1)*eval_gegenbauer(n-2,d+2,sn) if n>1 else np.zeros_like(z)
            uu=cs**d*poly
            udd=cs**d*((d*d*np.tan(z)**2-d/cs**2)*poly-(2*d+1)*sn*dp+cs**2*ddp)
            res=-udd+mu*uu/cs**2-(n+d)**2*uu
            residual=max(residual,float(np.max(abs(res)))/max(1,float(np.max(abs((n+d)**2*uu)))))
        assert max(gram_error,residual)<3e-12,(C,mu,gram_error,residual)
        rows.append({'C':C,'mu':mu,'Delta':d,'gram_error':gram_error,'mode_residual':residual})
    return rows


def retarded_pulses():
    rows=[]; T=math.pi+.25
    for C in [1.,61.,6001.]:
        vals=masses(6,C)
        frequencies=np.concatenate([np.arange(5)+(1+np.sqrt(1+4*v))/2 for v in vals])
        b=np.arange(1.,11.)**-1; b/=np.linalg.norm(b)
        def drive(t): return math.sin(math.pi*t/T)**2
        def rhs(t,y): return np.concatenate([y[10:],-frequencies**2*y[:10]+b*drive(t)])
        sol=solve_ivp(rhs,(0,T),np.zeros(20),method='DOP853',rtol=2e-12,atol=2e-14,dense_output=True)
        assert sol.success
        qref=np.array([bb*quad(lambda t:math.sin(om*(T-t))*drive(t)/om,0,T,epsabs=2e-13)[0] for om,bb in zip(frequencies,b)])
        pref=np.array([bb*quad(lambda t:math.cos(om*(T-t))*drive(t),0,T,epsabs=2e-13)[0] for om,bb in zip(frequencies,b)])
        discrepancy=float(max(np.max(abs(sol.y[:10,-1]-qref)),np.max(abs(sol.y[10:,-1]-pref))))
        energy=float(np.sum(pref**2+frequencies**2*qref**2)/2)
        bound=T*T/8
        work=quad(lambda t:float(np.dot(sol.sol(t)[10:],b))*drive(t),0,T,epsabs=2e-12)[0]
        assert discrepancy<2e-10 and abs(work-energy)<2e-10 and energy<=bound
        rows.append({'C':C,'duration':T,'Duhamel_ODE_error':discrepancy,'final_master_energy':energy,'source_impulse_energy_bound':bound,'work_balance_error':abs(work-energy)})
    return rows


def candidate_checks():
    mp.mp.dps=80; rows=[]
    for k0 in [2000,4000,10000]:
        k=mp.mpf(k0); C=mp.sqrt(1+12*mp.pi*k**6); L=mp.mpf(6); Y=(C+1)/2
        rad=mp.sqrt(C*C+4*L*Y)
        small=L*(L-2)/(C*(C+L+rad)); large=(C+L+rad)/C
        dsmall_minus_one=2*small/(1+mp.sqrt(1+4*small))
        pmetric=2*L*Y/(rad*(rad+C))
        assert abs(small-(C+L-rad)/C)<mp.mpf('1e-70')
        assert abs(small*large-L*(L-2)/C**2)<mp.mpf('1e-70')
        Lmax=C+1+mp.sqrt((C+1)**2+C*C)
        ell=int(mp.floor((mp.sqrt(1+4*Lmax)-1)/2))
        def mass_l(l):
            ll=mp.mpf(l*(l+1)); rr=mp.sqrt(C*C+4*ll*Y)
            return ll*(ll-2)/(C*(C+ll+rr))
        assert mass_l(ell)<=1 and mass_l(ell+1)>1
        count=(ell+1)**2-4
        rows.append({'k':k0,'C':str(C),'quadrupole_mu_minus':str(small),'quadrupole_mu_plus':str(large),
                     'Delta_minus_minus_one':str(dsmall_minus_one),'metric_source_light_projector_weight':str(pmetric),
                     'ell_max_mu_minus_le_one':ell,'angular_multiplicity_one_polar_branch':count,
                     'multiplicity_over_LLL_complex_channels':str(mp.mpf(count)/k**4),
                     'count_is_not_stress_estimate':True})
    return rows


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'evidence/quadrupole-product-response.json')
    args=parser.parse_args()
    out={'status':'PASS','scope':'classical linear sourced Einstein-Maxwell theory on fixed AdS2 x S2',
         'tensor_check':symbolic_tensor_check(),'exterior_check':exterior_limit_check(),
         'global_modes':mode_checks(),'retarded_pulses':retarded_pulses(),'candidate':candidate_checks(),
         'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'versions':{'sympy':s.__version__,'numpy':np.__version__,'scipy':scipy.__version__},
         'limitations':['Friedrichs global boundaries are a comparator, not physical horizon gluing.',
                        'Master energy normalization is mathematical; no physical stress bound is inferred.',
                        'Actual second-order dipole/fermion source correlations and renormalization remain.',
                        'No evolving wormhole, quantum capture, preparation, full stability or creation theorem.']}
    args.output.write_bytes((json.dumps(out,indent=2)+'\n').encode())
    print(json.dumps({'status':'PASS','tensor_components':out['tensor_check']['source_free_tensor_components'],
                      'sourced_identities':2,'global_mode_cases':len(out['global_modes']),'retarded_pulses':len(out['retarded_pulses']),
                      'candidate':out['candidate'][0]},indent=2),flush=True)
