"""Independent ordinary checks of the anomaly-completed response part."""
from pathlib import Path
import hashlib,json,platform,time,importlib.util
import numpy as np
import scipy
from scipy.integrate import quad
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
prior_spec=importlib.util.spec_from_file_location('moving_mouth',ROOT/'checks/moving_mouth_clock_vertex.py')
prior=importlib.util.module_from_spec(prior_spec)
prior_spec.loader.exec_module(prior)
spectral,sharp,finite_part=prior.spectral,prior.sharp,prior.finite_part


def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def pair(z):
    return [float(z.real),float(z.imag)]


def symbols():
    checks=[]
    def zero(name,value):
        assert sp.simplify(value)==0,(name,value)
        checks.append(name)
    t,z=sp.symbols('t z',real=True)
    a,b,g=sp.Function('a')(z),sp.Function('b')(z),sp.Function('g')(z)
    p=sp.Function('p')(t)
    h=sp.Matrix([[-2*a*p,g*sp.diff(p,t)],[g*sp.diff(p,t),2*b*p]])
    eta=sp.diag(-1,1);coords=[t,z]
    gamma=[[[sum(eta[i,l]*(sp.diff(h[l,j],coords[k])+sp.diff(h[l,k],coords[j])-sp.diff(h[j,k],coords[l])) for l in range(2))/2 for k in range(2)] for j in range(2)] for i in range(2)]
    R=sum(eta[j,k]*sum(sp.diff(gamma[i][j][k],coords[i])-sp.diff(gamma[i][j][i],coords[k]) for i in range(2)) for j in range(2) for k in range(2))
    zero('curvature from linearized Christoffel symbols',R-2*((b-sp.diff(g,z))*sp.diff(p,t,2)-sp.diff(a,z,2)*p))
    k,w,an,am,cn,cm=sp.symbols('k w an am cn cm')
    raw=(k*k*an-w*w*cn)*(k*k*am-w*w*cm)/(w*w-k*k)+k*k*an*am
    zero('curvature spectral polynomial split',raw-w*w*k*k*(an-cn)*(am-cm)/(w*w-k*k)-w*w*cn*cm)
    c,rhog,gv=sp.symbols('C rhog gamma')
    zero('background Weyl kinetic square',c*c-2*rhog*gv*c+rhog*rhog*gv*gv-(c-rhog*gv)**2)
    zero('pure spatial coordinate variation curvature',R.subs({a:0,b:sp.diff(g,z)}).doit())
    A,kap,rr=sp.symbols('A kap r',nonzero=True,real=True)
    zero('proper throat local kinetic combination',kap*(1+z*sp.tan(z))-kap-sp.tan(z)*kap*z)
    theta=2*sp.sin(k*A)/k
    gam=sp.I*kap*sp.diff(theta,k)
    zero('smoothed coordinate derivative profile',sp.I*k*gam-kap*theta+2*kap*A*sp.cos(k*A))
    L,S=sp.symbols('L S',positive=True)
    zero('endpoint Parseval cancellation including zero mode',rr**2*w*w/(12*sp.pi)*(4/L+8*S/L)-2*rr**2*w*w*S/(3*sp.pi*L)-rr**2*w*w/(3*sp.pi*L))
    return checks


def nonlinear_curvature_checks():
    # Fully nonlinear metric curvature by Christoffel jets, independently
    # compared with the coframe identity used to integrate the cocycle.
    t,z,scale=sp.symbols('t z scale',real=True)
    N=sp.exp(scale*(sp.sin(t)/5+sp.cos(z)/10))
    B=sp.exp(scale*(sp.cos(t)/9+sp.sin(z)/8))
    u=scale*sp.sin(t+z)/7
    metric=sp.Matrix([[-N*N+B*B*u*u,B*B*u],[B*B*u,B*B]])
    coords=[t,z]
    components=list(metric)
    raw=components+[sp.diff(v,x) for x in coords for v in components]+[sp.diff(v,x,y) for x in coords for y in coords for v in components]
    jets=sp.lambdify((t,z,scale),raw,'numpy')
    omz=(sp.diff(B,t)-sp.diff(B*u,z))/N
    omt=sp.diff(N,z)/B+u*omz
    cartan=sp.lambdify((t,z,scale),2*(sp.diff(omz,t)-sp.diff(omt,z))/(N*B),'numpy')
    rows=[]
    for sc in [.01,.3,1.]:
        for tt,zz in [(.13,.31),(.77,-.4),(1.2,1.8)]:
            vals=np.array(jets(tt,zz,sc),float)
            met=vals[:4].reshape(2,2); inv=np.linalg.inv(met)
            d=vals[4:12].reshape(2,2,2); dd=vals[12:].reshape(2,2,2,2)
            dinv=np.array([-inv@d[x]@inv for x in range(2)])
            ga=np.zeros((2,2,2));dga=np.zeros((2,2,2,2))
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        for l in range(2):
                            comb=d[j,l,k]+d[k,l,j]-d[l,j,k]
                            ga[i,j,k]+=inv[i,l]*comb/2
                            for x in range(2):
                                dga[x,i,j,k]+=(dinv[x,i,l]*comb+inv[i,l]*(dd[x,j,l,k]+dd[x,k,l,j]-dd[x,l,j,k]))/2
            ric=np.zeros((2,2))
            for j in range(2):
                for k in range(2):
                    for i in range(2):
                        ric[j,k]+=dga[i,i,j,k]-dga[k,i,j,i]
                        for l in range(2):
                            ric[j,k]+=ga[i,i,l]*ga[l,j,k]-ga[i,k,l]*ga[l,j,i]
            value=np.sum(inv*ric)
            ref=cartan(tt,zz,sc)
            err=abs(value-ref)
            assert err<3e-13
            rows.append(dict(scale=sc,t=tt,z=zz,christoffel=value,coframe=ref,error=err))
    return rows


def independent_fourier():
    rng=np.random.default_rng(73629)
    rows=[]
    for L in [3.4,7.,27.]:
        n=np.arange(1,19); k=2*np.pi*n/L
        a=(rng.normal(size=18)+1j*rng.normal(size=18))/(1+n*n)
        b=(rng.normal(size=18)+1j*rng.normal(size=18))/(1+n*n)
        g=(rng.normal(size=18)+1j*rng.normal(size=18))/(1+n*n)
        a0,b0=.13,.23
        C=b-1j*k*g; F=a-C
        grid=np.arange(8192)*L/8192
        phase=np.exp(1j*np.outer(grid,k))
        realC=b0+2*np.real(phase@C)
        J=L*np.mean(realC**2)
        for w in [1.1+.04j,.8+.3j,2.3+.1j]:
            raw=(k*k*a-w*w*C)*(k*k*np.conj(a)-w*w*np.conj(C))/(w*w-k*k)+k*k*np.abs(a)**2
            curvature=L/(6*np.pi)*raw.sum()+L*w*w*b0*b0/(12*np.pi)
            spectral_value=L/(6*np.pi)*np.sum(w*w*k*k*np.abs(F)**2/(w*w-k*k))
            alternative=spectral_value+w*w*J/(12*np.pi)
            err=abs(curvature-alternative)/(1+abs(curvature))
            assert err<3e-13
            rows.append(dict(L=L,w=pair(w),curvature=pair(curvature),spectral_plus_realspace=pair(alternative),normalized_error=err))
    return rows


def endpoint_integral(T,D,h):
    L=T+D; A=T/2
    centers=np.array([sign*A+j*L for j in range(-2,3) for sign in [-1,1]])
    def profile(z):
        return np.exp(-(z-centers)**2/(2*h*h)).sum()/(np.sqrt(2*np.pi)*h)
    points=sorted(set([-L/2,L/2]+[float(c+m*h) for c in centers for m in [-10,-5,0,5,10] if -L/2<c+m*h<L/2]))
    spatial=sum(quad(lambda z:profile(z)**2,left,right,epsabs=1e-10,epsrel=2e-12)[0] for left,right in zip(points,points[1:]))
    images=np.arange(-4,5)*L
    poisson=(np.exp(-images**2/(4*h*h)).sum()+np.exp(-(T-images)**2/(4*h*h)).sum())/(np.sqrt(np.pi)*h)
    assert abs(spatial-poisson)/(1+poisson)<2e-10
    return spatial,poisson


def stable_completed(T,D,r,h,w):
    L=T+D; nu=2*np.pi/L
    n=np.arange(1,int(np.ceil(9/(h*nu)))+1,dtype=float); k=nu*n
    heat=np.exp(-h*h*k*k); s=np.sin(k*T/2); co=np.cos(k*T/2)
    base=2*w*w/(3*np.pi*L)*np.sum((k*k-nu*nu)*s*s*heat/(k*k*(w*w-k*k)))
    cross=2*r*w*w/(3*np.pi*L)*np.sum((-1/k+(w*w-nu*nu)/(k*(w*w-k*k)))*np.sin(k*T)*heat)
    end=2*r*r*w*w*(w*w-nu*nu)/(3*np.pi*L)*np.sum(co*co*heat/(w*w-k*k))
    return base+cross+end+r*r*w*w/(3*np.pi*L)


def run():
    start=time.perf_counter()
    sym=symbols();curvature=nonlinear_curvature_checks();fourier=independent_fourier()
    cancellations=[]
    for D in [.0001,.1,1.,12.]:
        for h in [.03,.003,.0003]:
            T=3.;r=.04;w=1.1+.2j
            spatial,poisson=endpoint_integral(T,D,h)
            isolated,meta=spectral(T,D,r,h,w)
            completed=isolated+r*r*w*w*spatial/(12*np.pi)
            stable=stable_completed(T,D,r,h,w)
            err=abs(completed-stable)/(1+abs(stable))
            assert err<2e-10
            cancellations.append(dict(T=T,D=D,h=h,r=r,w=pair(w),spectral=pair(isolated),local=pair(r*r*w*w*spatial/(12*np.pi)),completed=pair(completed),stable=pair(stable),normalized_error=err,spatial_poisson_relative_error=abs(spatial-poisson)/(1+poisson),**meta))
    limits=[]
    for D in [.1,1.,7.]:
        for w in [1.1+.2j,2.1+.07j]:
            r=.05;T=2*np.arccos(r/np.sqrt(1+r*r))
            limit=finite_part(T,D,r,w)
            rows=[]
            for h in [.008,.004,.002,.001,.0005]:
                value=stable_completed(T,D,r,h,w)
                rows.append(dict(h=h,value=pair(value),error=abs(value-limit)))
            assert all(b['error']<.7*a['error'] for a,b in zip(rows,rows[1:]))
            limits.append(dict(D=D,w=pair(w),limit=pair(limit),refinements=rows))
    uniform=[]
    for D in [.0001,.01,1.,10.,100.]:
        for h,r in [(.005,.003),(.02,.02),(.005,.05)]:
            for eta in [.03,.2]:
                T=3.;w=1.3+1j*eta
                value=stable_completed(T,D,r,h,w)
                err=abs(value-sharp(T,D,w))
                expression=h+h*h/eta+(r+r*r)*(1+1/eta)
                uniform.append(dict(D=D,h=h,r=r,eta=eta,difference=err,bound_expression=expression,diagnostic_ratio=err/expression))
    scales=[]
    for q in [100,10000,1000000,100000000]:
        epsilon=q**-2;delta=np.log(q)**-2;r=epsilon/np.sqrt(1-epsilon**2)
        scales.append(dict(q=q,epsilon=epsilon,delta=delta,sharp_normalized_bound_expression=(r+r*r)*(1+1/delta)/delta,isolated_over_sharp_r2_over_h=r*r/epsilon**4))
    deps=['notes/covariant-mouth-inertia.md','notes/moving-mouth-clock-vertex.md','checks/moving_mouth_clock_vertex.py','notes/proper-clock-anomaly.md','notes/collar-conservative-response.md']
    result=dict(date='2026-10-11',status='passed',scope='Covariant anomaly Hessian plus known free vacuum commutator, not the full physical metric response',symbolic=sym,nonlinear_curvature=curvature,fourier_hessian=fourier,cancellations=cancellations,sharp_limits=limits,uniform=uniform,scaling=scales,
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__),
                anomaly_contact_derived=True,sharp_completed_response_part_controlled=True,full_state_contacts_complete=False,physical_family_matched=False,interacting_pressure_controlled=False,physical_damping_proved=False,novelty_established=False,lean_invocations=0,
                elapsed_seconds=time.perf_counter()-start,script_sha256=digest(Path(__file__)),input_sha256={p:digest(ROOT/p) for p in deps})
    (ROOT/'evidence/covariant-mouth-inertia.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',seconds=result['elapsed_seconds'],symbolic=len(sym),nonlinear_curvature=len(curvature),fourier=len(fourier),cancellations=len(cancellations),sharp_limits=len(limits),uniform=len(uniform),max_curvature_error=max(x['error'] for x in curvature),max_fourier_error=max(x['normalized_error'] for x in fourier),max_cancellation_error=max(x['normalized_error'] for x in cancellations),max_uniform_diagnostic_ratio=max(x['diagnostic_ratio'] for x in uniform))))


if __name__=='__main__':run()
