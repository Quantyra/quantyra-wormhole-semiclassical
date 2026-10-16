"""S137 static RN higher-mode boundary moment; no finite-drive certificate."""
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_key]='1'
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time
import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import solve_ivp,quad
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
ETA=Q(1,2000)
SCALAR_CASES=[('synthetic18','18','.05'),('synthetic32','32','.2'),('synthetic64','64','.001'),('candidate2000',None,None)]
INPUTS=['notes/rn-higher-mode-boundary-spectral-control.md','checks/rn_higher_mode_boundary_spectral_control.py',
        'notes/coupled-fermion-domain.md','notes/dirac-spherical-source-reduction.md',
        'notes/higher-landau-boundary.md','notes/curved-massive-response.md',
        'notes/exterior-lll-delayed-boundary.md','notes/rn-ads-finite-cut-matching.md',
        'notes/laboratory-drive-clock-feedback.md','notes/physical-clock-remainder.md']


def hashes():return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def check(condition,name):
    if not condition:raise AssertionError(name)


def exact_checks():
    sx=S.Matrix([[0,1],[1,0]]);sy=S.Matrix([[0,-S.I],[S.I,0]]);sz=S.diag(1,-1);id2=S.eye(2)
    kron=S.kronecker_product;id4=S.eye(4)
    gamma0=S.I*kron(sx,id2);gammar=kron(sy,id2)
    for ga in (kron(sz,sx),kron(sz,sy)):
        check(gammar*gamma0*ga-gamma0*ga*gammar==S.zeros(4),'physical projector angular commutator')
    check(gammar*gamma0*gammar+gamma0*gammar*gammar==S.zeros(4),'physical projector radial anticommutator')
    A=kron(sz,sy);B=kron(id2,sx);chi=kron(sz,sx);Pi=(id4+chi)/2;G=S.I*kron(sy,sx);C=S.I*A*G
    for expression in (chi*A+A*chi,chi*B-B*chi,(id4-Pi)*B*Pi,(id4-Pi)*A-A*Pi,A*G-G*A,B*G-G*B,C*C-id4):
        check(expression==S.zeros(4),'folded Clifford/domain identity')
    check(S.trace(Pi)==2 and Pi.rank()==2,'folded allowed boundary rank')
    check(C.H==C and C[:2,:2]==S.zeros(2) and C[2:,2:]==S.zeros(2),'source work is Hermitian off-diagonal')
    X,u,D=S.symbols('X u D',positive=True);CC=1+6*X;AA=1+3*X
    fe=X/u**2+1-2*(1+2*X)*u+(1+3*X)*u**2
    he=-u**2*S.diff(fe,u)/2-u*fe
    check(S.factor(he-u*(1-u)*(2*AA*u-1))==0,'extremal radial h')
    HH=X+2*X*u+AA*u*u
    check(S.expand(HH-CC*u*u-X*(1-u)*(1+3*u))==0,'radial positive polynomial')
    fn=fe-D*u;hn=-u*u*S.diff(fn,u)/2-u*fn
    check(S.factor(hn-he-3*D*u*u/2)==0,'mass perturbation h')
    residual=S.expand((hn*hn-CC*fn*u*u)-(he*he-CC*fe*u*u)-3*D*he*u*u-S.Rational(9,4)*D*D*u**4-CC*D*u**3)
    check(residual==0,'finite temperature square expansion')
    W,kappa,M=S.symbols('W kappa M',positive=True)
    young=(1-1/M)*W*W-S.Rational(4,3)*kappa*W-(W*W/2-kappa*kappa)
    check(S.expand(young-((2*W/3-kappa)**2+(S.Rational(1,18)-1/M)*W*W))==0,'shifted Young square')
    s=S.symbols('s',positive=True);JJ=CC+4*X*s+X*s*s
    ff=S.factor(fe.subs(u,1/(1+s)))
    check(S.factor(ff-s*s*JJ/(1+s)**2)==0,'horizon mass factor')
    # Uniform rational family comparisons, retaining frozen eta.
    k0=2000;r=Q(1001,1000);assertions={
        'X_ge_k3':Q(31,10)>3+Q(1,k0**3),
        'C_lt_123k3_over20':36*Q(22,21)<(Q(123,20)-Q(1,k0**3))**2,
        'M1_gt18':Q(20*k0,123)>324,
        'C_over_q_lt7_over_k':Q(123,20)<7,
        'counting_K_over_q_lt15_over_k':14*(1+Q(1,k0**4))<15,
        'sixth_power_bound':r**6<Q(1007,1000),
        'CD_bound':r**7<Q(1009,1000),
        'D_square_bound':r**8<Q(101,100),
        'thermal_h_bound':Q(7,4)*Q(1009,1000)+Q(9,4)*Q(101,100)*Q(1,10**6)<Q(16,9),
        'thermal_horizon_small':Q(112,25)*ETA/k0**9<Q(1,1000),
        'collar_a_d_small':Q(4,3)*Q(7,25)*2*Q(8,3)*ETA/k0**9<Q(1,1000),
        'collar_lower_mass':Q(2000,999)<3,
        'sqrt2_lower':Q(141,100)**2<2,
        'sqrt2_upper':Q(17,12)**2>2,
        'sqrt7_upper':Q(8,3)**2>7,
        'exponential_tower_bound':sum(Q(423,50)**j/Q(math.factorial(j)) for j in range(51))>4000,
        'angular_tail_correction':1+Q(2,1)/Q(423,50)<Q(5,4),
        'moment_prefactor':1/(1-Q(1,4000))**2<2,
        'tower_target':Q(230000)*ETA**3/k0**28+Q(20)*ETA**2/k0**19<Q(1,10**68),
        'first_coefficient':85680*Q(8,3)<230000,
        'second_coefficient':Q(45,64)*28<20,
        'TFD_laplace':1+4/Q(775,70)**2+24/Q(775,70)**4<Q(21,20),
        'coherent_work':Q(21,20)*Q(22,7)<4,
    }
    for name,truth in assertions.items():check(truth,name)
    z,y,k,a,ell=S.symbols('z y k a ell',positive=True)
    den=k*S.sinh(k*ell)+a*S.cosh(k*ell)
    gi=(S.cosh(k*(ell-y))+a/k*S.sinh(k*(ell-y)))/den
    go=S.exp(-a*(y-ell))/den
    for value in (-S.diff(gi,y,2)+k*k*gi,-S.diff(go,y,2)+a*a*go,
                  gi.subs(y,ell)-go.subs(y,ell),S.diff(gi,y).subs(y,ell)-S.diff(go,y).subs(y,ell),S.diff(gi,y).subs(y,0)+1):
        check(S.simplify(value)==0,'step barrier Green identity')
    xx,ss=S.symbols('xx ss',positive=True)
    switch=10*xx**3-15*xx**4+6*xx**5
    check(S.factor(S.diff(switch,xx)-30*xx**2*(1-xx)**2)==0,'monotone frozen quintic')
    check(switch.subs(xx,0)==0 and switch.subs(xx,1)==1 and S.integrate(switch,(xx,0,1))==S.Rational(1,2),'frozen source endpoints/area')
    check(S.simplify(ss*S.integrate(S.exp(-ss*xx)*(1+xx**2)**2,(xx,0,S.oo))-(1+4/ss**2+24/ss**4))==0,'TFD positive-measure Laplace factor')
    print('EXACT Clifford/domain, RN, scalarGreen and',len(assertions),'rational inequalities PASS',flush=True)
    return {'symbolic_identities':'pass','rational_assertions':list(assertions),'boundary_allowed_rank':2,
            'tower_upper_rational':str(Q(230000)*ETA**3/k0**28+Q(20)*ETA**2/k0**19),
            'tower_upper_decimal':float(Q(230000)*ETA**3/k0**28+Q(20)*ETA**2/k0**19),
            'frozen_source_total_angle':'pi/500','integrated_frame_operator_norm':'pi/1000'}


def candidate(k):
    kk=mp.mpf(k);X=(-1+mp.sqrt(1+12*mp.pi*kk**6))/6;C=1+6*X
    delta=2*mp.mpf(1)/2000*mp.sqrt(C/X)/(kk**6*X)
    M=mp.sqrt((kk**4+1)/C)
    return X,C,delta,M


def parameters(label,M,delta):
    if M is None:
        _,_,delta,M=candidate(2000)
    return mp.mpf(M),mp.mpf(delta)


def closed_moment(M,delta):
    mu=M/(3*mp.sqrt(2));a=mp.mpf(24)*delta/25;k=mp.sqrt(mu*mu+a*a)
    den=k*mp.sinh(k)+a*mp.cosh(k)
    inner=(1-a*a/k**2)/2+(1+a*a/k**2)*mp.sinh(2*k)/(4*k)+a*mp.sinh(k)**2/k**2
    R=(inner+1/(2*a))/den**2
    bound=2/(1-mp.exp(-2*k))**2*(1/k**3+mp.exp(-2*k)/(a*k*k))
    simple=4*(1/mu**3+mp.exp(-2*mu)/(a*mu*mu))
    return R,bound,simple,mu,a,k


def ode_moment(M,delta,refined):
    R,bound,simple,mu,a,k=closed_moment(M,delta)
    af=float(a);kf=float(k)
    tol=2e-12 if refined else 2e-11
    sol=solve_ivp(lambda y,v:[v[1],kf*kf*v[0]],(1.,0.),[1.,-af],method='DOP853',rtol=tol,atol=tol*1e-3,dense_output=True,max_step=.01 if refined else .02)
    check(sol.success,'backward scalar ODE')
    den=-sol.y[1,-1]
    inner=quad(lambda y:(sol.sol(y)[0]/den)**2,0.,1.,epsabs=tol,epsrel=tol)[0]
    return mp.mpf(inner)+1/(2*a*mp.mpf(den)**2)


def continuum_moment(M,delta,dps):
    with mp.workdps(dps):
        R,bound,simple,mu,a,k=closed_moment(M,delta)
        P=max(mp.mpf(500),4*mu)
        check(P>=2*a and P>=2*mu,'spectral UV comparison range')
        def density(p):
            if p==mu:return 2/mp.pi
            if p<mu:
                g=mp.sqrt(mu*mu-p*p)
                den=mp.cosh(g)**2+(g/p)**2*mp.sinh(g)**2
            else:
                q=mp.sqrt(p*p-mu*mu)
                den=mp.cos(q)**2+(q/p)**2*mp.sin(q)**2
            return 2/(mp.pi*den)
        vtop=mp.log(mu/a)
        breaks=sorted(set([mp.mpf(-50),vtop]+[mp.mpf(v) for v in (-10,0,2,10,20,40,60) if -50<v<float(vtop)]))
        low_scaled=mp.quad(lambda v:a*density(a*mp.exp(v))*(a*mp.exp(v))/(a*a+(a*mp.exp(v))**2)**2,breaks)
        high_breaks=[mu];step=mp.pi/(4 if dps>80 else 2)
        while high_breaks[-1]+step<P:high_breaks.append(high_breaks[-1]+step)
        high_breaks.append(P)
        high=mp.quad(lambda p:density(p)/(a*a+p*p)**2,high_breaks)
        # Exact free UV tail plus rigorous bound on its scattering correction.
        free_tail=2/mp.pi*sum((-1)**j*(j+1)*a**(2*j)/((2*j+3)*P**(2*j+3)) for j in range(14))
        free_series_error=4/mp.pi*15*a**28/(31*P**31)
        uv_error=8*mu*mu/(15*mp.pi*P**5)+free_series_error
        pmin=a*mp.exp(-50)
        ir_error=8*pmin**3/(9*mp.pi*a**4*mu**2*mp.sinh(mp.sqrt(3)*mu/2)**2)
        value=low_scaled/a+high+free_tail
        return +value,+(uv_error+ir_error),{'P':mp.nstr(P,20),'precision_digits':dps,'positive_scattering_measure':'Neumann step barrier, full continuum','UV_difference_bound':mp.nstr(uv_error,20),'IR_bound':mp.nstr(ir_error,20)}


def scalar_checks():
    out=[]
    with mp.workdps(120):
        for label,Mraw,draw in SCALAR_CASES:
            M,delta=parameters(label,Mraw,draw);R,bound,simple,mu,a,k=closed_moment(M,delta)
            check(0<R<=bound<=simple,'analytic scalar bounds')
            errors=[];spectral=[]
            for refined,dps,tolerance in ((False,80,mp.mpf('2e-8')),(True,110,mp.mpf('2e-9'))):
                od=ode_moment(M,delta,refined)
                sd,tail,meta=continuum_moment(M,delta,dps)
                ode_error=abs(od-R)/(1+abs(R));spectral_error=max(mp.mpf(0),abs(sd-R)-tail)/(1+abs(R))
                check(ode_error<tolerance and spectral_error<tolerance,'independent scalar comparison '+label)
                errors.append({'precision_digits':dps,'ODE_error_scaled':mp.nstr(ode_error,20),'spectral_error_after_tail_scaled':mp.nstr(spectral_error,20),'raw_spectral_difference':mp.nstr(abs(sd-R),20)})
                spectral.append(meta)
            out.append({'case':label,'M':mp.nstr(M,35),'nu_d':mp.nstr(delta,35),'moment':mp.nstr(R,35),
                        'exact_upper':mp.nstr(bound,35),'simple_upper':mp.nstr(simple,35),'comparisons':errors,'spectral_metadata':spectral})
            print('SCALAR',label,'ODE/spectral/upper bounds PASS',flush=True)
    return out


def geometry_checks():
    rows=[]
    with mp.workdps(140):
        specs=[]
        for kk in (2000,4000):
            X,C,delta,M=candidate(kk)
            for zi in (mp.mpf(1)/4,mp.mpf(7)/25):
                nu=delta*mp.sqrt(C*X);kap=zi*nu
                temp=lambda s:(2*C*s+(1+18*X)*s*s+12*X*s**3+3*X*s**4)/(2*(1+s)**3)
                sh=mp.findroot(lambda s:temp(s)-kap,(kap/C,2*kap/C))
                specs.append(('physical',X,sh,nu,M,kk,zi))
        for X in (mp.mpf(1),mp.mpf(10)**10):
            for sh in (mp.mpf(1)/2000,mp.mpf(1)/1000):specs.append(('synthetic',X,sh,mp.mpf(1),mp.mpf(18),None,None))
        for kind,X,sh,nu,M,kk,zi in specs:
            C=1+6*X;D=sh*sh*(C+4*X*sh+X*sh*sh)/(1+sh)
            f=lambda r:X*r*r+1-(2*(1+2*X)+D)/r+(1+3*X)/(r*r)
            kap=mp.diff(f,1+sh)/2
            lam=M*mp.sqrt(C);W=lambda r:lam*mp.sqrt(f(r))/r
            max_chain=mp.mpf(0);min_margin=mp.inf
            for rr in [1+sh*mp.mpf(v) for v in ('1.01','2','10')]+[mp.mpf('1.01'),mp.mpf(2),mp.mpf(10),mp.mpf(10)**6]:
                if rr<=1+sh:continue
                val=W(rr);wp=f(rr)*mp.diff(W,rr);h=mp.diff(f,rr)/2-f(rr)/rr
                chain=abs(wp-val*h)/(1+abs(wp));max_chain=max(max_chain,chain)
                margin=val*val-abs(wp)-(val*val/2-kap*kap);min_margin=min(min_margin,margin/(1+val*val+kap*kap))
                check(chain<mp.mpf('2e-10') and margin>=-mp.mpf('2e-10')*(1+val*val+kap*kap),'RN geometry comparison')
            rows.append({'kind':kind,'k':kk,'z_i':None if zi is None else mp.nstr(zi),'X':mp.nstr(X,25),'s_h':mp.nstr(sh,25),'kappa0':mp.nstr(kap,25),'maximum_radial_chain_error_scaled':mp.nstr(max_chain,20),'minimum_shifted_margin_scaled':mp.nstr(min_margin,20)})
    print('GEOMETRY',len(rows),'independent direct metric cases PASS',flush=True)
    return rows


def main():
    started=time.time();frozen=hashes();path=ROOT/'evidence/rn-higher-mode-boundary-spectral-control.json'
    receipt={'scope':'classical static magnetic near-extremal RN free Dirac; first-order specified TFD work only',
             'eta':'1/2000','finite_drive_work_proved':False,'quantum_exterior_transfer_proved':False,
             'source_sha256':frozen,'precision':{'symbolic':'exact SymPy and Fraction','spectral_digits':[80,110],'geometry_digits':140},
             'predeclared_tolerances':{'scalar':[2e-8,2e-9],'geometry_relative':2e-10},'failures':[]}
    if path.exists():
        previous=json.loads(path.read_text(encoding='utf-8'))
        receipt['prior_attempts']=previous.get('prior_attempts',[])+[{k:v for k,v in previous.items() if k!='prior_attempts'}]
    try:
        receipt['exact']=exact_checks();receipt['geometry']=geometry_checks();receipt['scalar']=scalar_checks()
        check(hashes()==frozen,'sources changed during run')
        receipt['outcome']='PASS_STATIC_REFERENCE_COMPONENT_ONLY';status=0
    except Exception as exc:
        receipt['failures'].append({'type':type(exc).__name__,'message':str(exc)})
        receipt['outcome']='FAIL';status=1
        import traceback
        traceback.print_exc()
    receipt.update({'elapsed_seconds':time.time()-started,'terminal_exit_status':status,'python':sys.version,'platform':platform.platform(),
                    'versions':{'numpy':np.__version__,'scipy':scipy.__version__,'sympy':S.__version__,'mpmath':mp.__version__}})
    path.write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(receipt['outcome'],'elapsed',receipt['elapsed_seconds'],flush=True)
    return status


if __name__=='__main__':sys.exit(main())
