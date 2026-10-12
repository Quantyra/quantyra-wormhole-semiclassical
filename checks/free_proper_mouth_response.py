"""Checks the complete prescribed free matter kernel; no Lean/cloud calls."""
from pathlib import Path
import hashlib,json,platform,time,importlib.util
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import erf
import sympy as sp
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]


def module(name,file):
    spec=importlib.util.spec_from_file_location(name,ROOT/'checks'/file)
    out=importlib.util.module_from_spec(spec);spec.loader.exec_module(out)
    return out


prior=module('moving','moving_mouth_clock_vertex.py')
cov=module('covariant','covariant_mouth_inertia.py')


def digest(path):return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def pair(z):return [float(z.real),float(z.imag)]


def symbolic():
    results=[]
    def zero(name,value):
        assert sp.simplify(value)==0,(name,value)
        results.append(name)
    t,z=sp.symbols('t z',real=True)
    a,b,g=sp.Function('a')(z),sp.Function('b')(z),sp.Function('g')(z)
    phi=sp.Function('phi')(t);e=sp.symbols('e')
    h=sp.Matrix([[-2*a*phi,g*sp.diff(phi,t)],[g*sp.diff(phi,t),2*b*phi]])
    eta=sp.diag(-1,1);coords=[t,z]
    ga=[[[sum(eta[i,l]*(sp.diff(h[l,j],coords[k])+sp.diff(h[l,k],coords[j])-sp.diff(h[j,k],coords[l])) for l in range(2))/2 for k in range(2)] for j in range(2)] for i in range(2)]
    # Vary both Christoffel terms in covariant conservation of T^{mu nu}.
    source=[e*(sum(ga[mu][mu][nu] for mu in range(2))+sum(ga[nu][mu][mu] for mu in range(2))) for nu in range(2)]
    zero('time Ward source from connection variation',source[0]-e*(2*a+2*b-sp.diff(g,z))*sp.diff(phi,t))
    zero('space Ward source with varied trace',source[1]+2*e*sp.diff(a-b,z)*phi-e*(4*sp.diff(a,z)*phi+g*sp.diff(phi,t,2)))
    an,bn,gn,k,w=sp.symbols('an bn gn k w')
    F=an-bn+sp.I*k*gn
    U=-4*e*an+2*e*w*w*F/(w*w-k*k)
    P=sp.I*e*w*gn+2*e*k*w*F/(w*w-k*k)
    zero('retarded energy Ward solution',-sp.I*w*U+sp.I*k*P-sp.I*w*e*(2*an+2*bn-sp.I*k*gn))
    zero('retarded momentum Ward solution',-sp.I*w*P+sp.I*k*U-e*(-4*sp.I*k*an+w*w*gn))
    am,bm,gm=sp.symbols('am bm gm');Fm=am-bm-sp.I*k*gm
    dyn=-(am-bm)*(U+4*e*an)+sp.I*w*gm*P+e*w*w*gm*gn
    zero('state metric seagull cancels local inertia',dyn+2*e*w*w*Fm*F/(w*w-k*k))
    eps,p,pd=sp.symbols('epsilon p pd',positive=True)
    bcut=sp.sqrt(1-eps*eps);s1=-eps*p/bcut
    s2=-eps*p*p/(2*bcut**3)+eps**3*pd*pd/(2*bcut**3)
    zero('proper-cut first amplitude coefficient',-2*eps*bcut*s1-2*eps*eps*p)
    zero('proper-cut second amplitude coefficient',-2*eps*bcut*s2+(1-2*eps*eps)*s1*s1-eps*eps*(2*p*p-eps*eps*pd*pd/(bcut*bcut)))
    rho=sp.Function('rho')(z)
    Tt=e+(2*sp.diff(rho,z,2)-sp.diff(rho,z)**2)/(24*sp.pi)
    Tz=e-sp.diff(rho,z)**2/(24*sp.pi)
    zero('static covariant stress conservation',sp.diff(Tz,z)+sp.diff(rho,z)*(Tt-Tz))
    Z=sp.Function('Z')(z);A,theta,d=sp.symbols('A theta d')
    zero('trajectory integration by parts with seam',sp.diff(sp.diff(rho,z)**2*Z,z).subs(sp.diff(Z,z),theta-A*d)-(2*sp.diff(rho,z)*sp.diff(rho,z,2)*Z+sp.diff(rho,z)**2*theta-A*sp.diff(rho,z)**2*d))
    pexp,qexp=sp.symbols('pexp qexp')
    Q=sp.I*(1-pexp)*(1-qexp)/(2*(1-pexp*qexp))
    zero('return-uniform length derivative',sp.diff(Q,pexp)*sp.I*w*pexp-w*pexp*(1-qexp)**2/(2*(1-pexp*qexp)**2))
    return results


def ward_modes():
    rng=np.random.default_rng(9461);rows=[]
    for L in [3.4,7.,23.]:
        count=12
        def coeff():
            pos=(rng.normal(size=count)+1j*rng.normal(size=count))/(1+np.arange(1,count+1)**2)
            return np.r_[pos[::-1].conj(),rng.normal(),pos]
        a,b,g=coeff(),coeff(),coeff();k=2*np.pi*np.arange(-count,count+1)/L
        e=-np.pi/(6*L*L);f=a-b;F=f+1j*k*g
        def restoring(w):
            U=np.zeros_like(a);P=np.zeros_like(a)
            for j,kk in enumerate(k):
                if w==0 and kk==0:
                    U[j]=-2*e*(a[j]+b[j]);continue
                mat=np.array([[-1j*w,1j*kk],[1j*kk,-1j*w]])
                rhs=np.array([1j*w*e*(2*a[j]+2*b[j]-1j*kk*g[j]),e*(-4j*kk*a[j]+w*w*g[j])])
                U[j],P[j]=np.linalg.solve(mat,rhs)
            density=-a.conj()*(U+e*(a+b))+b.conj()*(U+2*e*f+e*(a+b))+1j*w*g.conj()*P
            seagull=2*e*(abs(b)**2-abs(a)**2)+e*w*w*abs(g)**2
            return -L*np.sum(density+seagull)
        baseline=restoring(0.)
        mask=k!=0
        for w in [1.1+.03j,.8+.2j,2.1+.09j]:
            direct=restoring(w)-baseline
            spectrum=2*e*L*np.sum(w*w*abs(F[mask])**2/(w*w-k[mask]**2))
            err=abs(direct-spectrum)/(1+abs(spectrum))
            assert err<3e-13
            rows.append(dict(L=L,w=pair(w),ward=pair(direct),state_spectrum=pair(spectrum),normalized_error=err))
    return rows


def profiles(z,T,D,h):
    L=T+D;A=T/2
    y=z-np.arange(-2,3)*L
    gp=np.exp(-(y+A)**2/(2*h*h))/(np.sqrt(2*np.pi)*h)
    gm=np.exp(-(y-A)**2/(2*h*h))/(np.sqrt(2*np.pi)*h)
    theta=(erf((y+A)/(np.sqrt(2)*h))-erf((y-A)/(np.sqrt(2)*h)))/2
    return theta.sum(),(gp+gm).sum(),(y*theta+h*h*(gp-gm)).sum()


def trajectory_case(epsilon,D,h):
    A=np.arccos(epsilon);T=2*A;L=T+D;nu=2*np.pi/L
    # Independent general-background check. This periodic cosine has the
    # required endpoint jets; it is not an Einstein-matched AdS throat.
    amp=-np.sqrt(1-epsilon**2)/(epsilon*nu*np.sin(nu*A))
    def rhojets(z):return -amp*nu*np.sin(nu*z),-amp*nu*nu*np.cos(nu*z)
    centers=[s*A+j*L for j in [-1,0,1] for s in [-1,1]]
    points=sorted(set([-L/2,L/2]+[float(c+m*h) for c in centers for m in [-8,0,8] if -L/2<c+m*h<L/2]))
    def integrate(fn):
        return sum(quad(fn,l,r,epsabs=2e-10,epsrel=3e-12,limit=150)[0] for l,r in zip(points,points[1:]))
    e0=-np.pi/(6*L*L)
    def rho_d(z):
        p,_=rhojets(z);_,d,_=profiles(z,T,D,h);return p*p*d
    edge=integrate(rho_d)
    expected=2*e0-edge/(24*np.pi)
    refinements=[]
    for step in [.004,.002,.001]:
        def jprofile(z,sgn):
            th,_,_=profiles(z,T,D,h);return np.exp(sgn*step*th/A)
        jp=integrate(lambda z:jprofile(z,1));jm=integrate(lambda z:jprofile(z,-1))
        state_fd=np.pi*(1/jp-1/jm)/(12*step)
        def odd(z):
            p,pp=rhojets(z);th,d,Z=profiles(z,T,D,h)
            ap=(pp*Z+p*(th-A*d))/A;t=th/A
            return -np.sinh(step*t)/step*(p*p+step*step*ap*ap)+2*p*ap*np.cosh(step*t)
        fd=state_fd+integrate(odd)/(24*np.pi)
        err=abs(fd-expected)/(1+abs(expected))
        refinements.append(dict(step=step,derivative=fd,normalized_error=err))
    assert all(b['normalized_error']<.3*a['normalized_error'] for a,b in zip(refinements,refinements[1:]))
    assert refinements[-1]['normalized_error']<3e-5
    analytic_edge=amp*amp*nu*nu*(1-np.cos(2*nu*A)*np.exp(-2*h*h*nu*nu))
    assert abs(edge-analytic_edge)/(1+abs(edge))<2e-11
    ds=epsilon**3/(2*(1-epsilon**2)**1.5)
    return dict(epsilon=epsilon,D=D,h=h,background='Periodic cosine with product-mouth jets; a check of the general stress identity',stress_derivative=expected,endpoint_integral=edge,analytic_endpoint_error=abs(edge-analytic_edge)/(1+abs(edge)),finite_trajectory_coefficient=-2*ds*expected,refinements=refinements)


def trajectory_coefficient(eps,D):
    b=np.sqrt(1-eps*eps);L=2*np.arccos(eps)+D
    return eps/(12*np.pi*b)+np.pi*eps**3/(3*L*L*b**3)


def complete(eps,D,w):
    T=2*np.arccos(eps);r=eps/np.sqrt(1-eps*eps)
    return prior.finite_part(T,D,r,w)+w*w*trajectory_coefficient(eps,D)


def zero_frequency():
    mp.mp.dps=60;rows=[]
    for eps,D in [('.03','.2'),('.005','8')]:
        ep,dd=mp.mpf(eps),mp.mpf(D);b=mp.sqrt(1-ep*ep)
        T=2*mp.acos(ep);L=T+dd;r=ep/b;nu=2*mp.pi/L
        coefficient=ep/(12*mp.pi*b)+mp.pi*ep**3/(3*L*L*b**3)
        values=[]
        for scale in ['.01','.005','.0025']:
            w=mp.mpf(scale)*(1+mp.j/5);U=w*L/2;Z=w*(L/2-T)
            q1=(L-2*T)/4
            q2=(q1-L*mp.sin(Z)/(4*mp.sin(U)))/(w*w)
            q3=L*(mp.cos(U)+mp.cos(Z))/(8*w*mp.sin(U))-1/(2*w*w)
            sharp=2/(3*mp.pi*L)*(-nu*nu*T*dd/8-(w*w-nu*nu)*L*mp.sin(w*T/2)*mp.sin(w*dd/2)/(4*w*mp.sin(U)))
            val=sharp+2*r*w*w/(3*mp.pi*L)*(-q1+(w*w-nu*nu)*q2)+r*r*w*w/(3*mp.pi*L)+2*r*r*w*w*(w*w-nu*nu)*q3/(3*mp.pi*L)+w*w*coefficient
            values.append(dict(scale=scale,value=pair(val),over_w2=pair(val/(w*w))))
        ratios=[complex(*v['over_w2']) for v in values]
        assert abs(ratios[2]-ratios[1])<.3*abs(ratios[1]-ratios[0])
        rows.append(dict(epsilon=eps,D=D,precision=60,refinements=values))
    return rows


def run():
    start=time.perf_counter();syms=symbolic();ward=ward_modes()
    print(json.dumps(dict(phase='state_Ward',status='passed')),flush=True)
    trajectory=[trajectory_case(e,D,h) for e in [.1,.03] for D in [.5,2.] for h in [.025,.008]]
    print(json.dumps(dict(phase='trajectory_energy_derivative',status='passed')),flush=True)
    cut=[]
    for eps in [.2,.03,.003]:
        b=np.sqrt(1-eps*eps);A=np.arccos(eps)
        for frequency in [.7,2.1]:
            p=np.cos(frequency*.37);pd=-frequency*np.sin(frequency*.37);pdd=-frequency**2*p
            refinements=[]
            for amp in [.04,.02,.01]:
                s=A-eps/b*amp*p+eps/(2*b**3)*amp**2*(-p*p+eps**2*pd*pd)
                sd=-eps/b*amp*pd+eps/b**3*amp**2*(-p*pd+eps**2*pd*pdd)
                residual=abs(np.cos(s)**2-eps**2*(np.exp(2*amp*p)-sd*sd))/eps**2
                refinements.append(dict(amplitude=amp,residual=residual))
            assert all(v['residual']<.14*u['residual'] for u,v in zip(refinements,refinements[1:]))
            cut.append(dict(epsilon=eps,frequency=frequency,refinements=refinements))
    limits=[]
    for eps in [.03,.005]:
        for D in [.2,2.,9.]:
            A=np.arccos(eps);T=2*A;L=T+D;nu=2*np.pi/L;r=eps/np.sqrt(1-eps**2)
            amp=-np.sqrt(1-eps**2)/(eps*nu*np.sin(nu*A));w=1.2+.1j
            target=complete(eps,D,w);rows=[]
            for h in [.008,.004,.002,.001]:
                edge=amp*amp*nu*nu*(1-np.cos(2*nu*A)*np.exp(-2*h*h*nu*nu))
                Ih=-np.pi/(3*L*L)-edge/(24*np.pi)
                ds=eps**3/(2*(1-eps*eps)**1.5)
                value=cov.stable_completed(T,D,r,h,w)-2*ds*Ih*w*w
                rows.append(dict(h=h,value=pair(value),error=abs(value-target)))
            assert all(v['error']<.7*u['error'] for u,v in zip(rows,rows[1:]))
            limits.append(dict(epsilon=eps,D=D,w=pair(w),sharp=pair(target),refinements=rows))
    uniform=[]
    for eps in [.03,.003,.0003]:
        for D in [.0001,.01,.3,3.,30.,300.]:
            for eta in [.02,.2]:
                w=1.1+1j*eta
                value=complete(eps,D,w);reference=prior.sharp(np.pi,D,w)
                expr=eps*(1+eta**-2)
                uniform.append(dict(epsilon=eps,D=D,eta=eta,error=abs(value-reference),bound_expression=expr,diagnostic_ratio=abs(value-reference)/expr))
    zeros=zero_frequency();scales=[]
    for q in [100,10000,1000000,100000000]:
        eps=q**-2;delta=np.log(q)**-2
        scales.append(dict(q=q,epsilon=eps,delta=delta,normalized_complete_bound_expression=eps*(1+delta**-2)/delta))
    deps=['notes/free-proper-mouth-response.md','notes/covariant-mouth-inertia.md','notes/moving-mouth-clock-vertex.md','notes/proper-clock-anomaly.md','checks/covariant_mouth_inertia.py','checks/moving_mouth_clock_vertex.py']
    data=dict(date='2026-10-11',status='passed',scope='Complete statically-subtracted free CFT matter response of the prescribed regularized proper-mouth family',symbolic=syms,state_ward=ward,trajectory=trajectory,cut_amplitude=cut,sharp_limits=limits,uniform=uniform,zero_frequency=zeros,scaling=scales,
              environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__,mpmath=mp.__version__),
              free_state_contacts_complete=True,quadratic_proper_trajectory_included=True,prescribed_free_matter_response_complete=True,physical_metric_matched=False,interacting_pressure_controlled=False,physical_damping_proved=False,novelty_established=False,lean_invocations=0,
              script_sha256=digest(Path(__file__)),input_sha256={p:digest(ROOT/p) for p in deps},elapsed_seconds=time.perf_counter()-start)
    (ROOT/'evidence/free-proper-mouth-response.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',seconds=data['elapsed_seconds'],symbolic=len(syms),state_ward=len(ward),trajectory=len(trajectory),cut=len(cut),limits=len(limits),uniform=len(uniform),zero_frequency=len(zeros),max_state_ward_error=max(r['normalized_error'] for r in ward),max_trajectory_fd_error=max(r['refinements'][-1]['normalized_error'] for r in trajectory),max_uniform_ratio=max(r['diagnostic_ratio'] for r in uniform))))


if __name__=='__main__':run()
