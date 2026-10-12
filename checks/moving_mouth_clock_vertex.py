"""Supplementary ordinary checks; no Lean, cloud, or physical-stability claim."""
from pathlib import Path
import hashlib
import json
import platform
import time
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.special import erf
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def symbolic():
    checks = []

    def zero(name, expr):
        assert sp.simplify(sp.expand_trig(expr)) == 0, (name, expr)
        checks.append(name)

    J, x, u, p, y, d, a, z, kap, phi = sp.symbols('J x u p y d a z kap phi', real=True)
    lag = J / (2*x) * (d-u*y)**2 - x*y*y/(2*J)-a*x*J*sp.sec(J*z)**2/2
    vel = x*p/J+u*y
    ham = (p*d-lag).subs(d, vel)
    zero('Legendre transform including shift', ham-(x*(p*p+y*y)/(2*J)+a*x*J*sp.sec(J*z)**2/2+u*p*y))
    zero('linear kinetic coefficient', sp.diff(sp.exp(phi)/(1+kap*phi), phi).subs(phi, 0)-(1-kap))
    zero('linear potential coefficient', sp.diff(a*sp.exp(phi)*(1+kap*phi)*sp.sec((1+kap*phi)*z)**2, phi).subs(phi, 0)-a*sp.sec(z)**2*(1+kap+2*kap*z*sp.tan(z)))
    e = sp.symbols('epsilon', positive=True)
    zero('proper-cut first variation', sp.diff(sp.acos(e*sp.exp(phi)), phi).subs(phi, 0)+e/sp.sqrt(1-e*e))
    D = sp.symbols('D', positive=True)
    length = D+2*sp.acos(e*sp.exp(phi))/sp.exp(phi)
    zero('free vacuum clock derivative agrees with pressure', sp.diff(-sp.pi/(6*length),phi).subs(phi,0)+sp.pi*(2*sp.acos(e)+2*e/sp.sqrt(1-e*e))/(6*(D+2*sp.acos(e))**2))
    f, g, pot = sp.Function('f')(z), sp.Function('g')(z), sp.Function('V')(z)
    wa, wb = sp.symbols('wa wb', positive=True)
    fp, gp = sp.diff(f,z), sp.diff(g,z)
    B = z*(fp*gp-(pot+wa*wb)*f*g)
    db = sp.diff(B,z).subs({sp.diff(f,z,2):(pot-wa*wa)*f, sp.diff(g,z,2):(pot-wb*wb)*g})
    target = fp*gp-(pot+wa*wb)*f*g-z*sp.diff(pot,z)*f*g-(wa+wb)*z*(wa*f*gp+wb*g*fp)
    zero('local on-shell Ward identity', db-target)
    k, w, nu = sp.symbols('k w nu')
    zero('cross spectral decomposition', (k*k-nu*nu)/(k*(w*w-k*k))+1/k-(w*w-nu*nu)/(k*(w*w-k*k)))
    zero('endpoint spectral decomposition', (k*k-nu*nu)/(w*w-k*k)+1-(w*w-nu*nu)/(w*w-k*k))
    t, r = sp.symbols('t r', real=True)
    ray = 1-sp.exp(2*sp.I*t)-2*sp.I*t*r*(1+sp.exp(2*sp.I*t))
    zero('ray form factor', sp.expand_complex(ray+2*sp.I*sp.exp(sp.I*t)*(sp.sin(t)+2*t*r*sp.cos(t))))
    return checks


def ward_case(epsilon, a, parity, frequencies):
    A = np.arccos(epsilon)
    kap = -epsilon/(A*np.sqrt(1-epsilon**2))
    wa, wb = frequencies
    init = [1.,0.] if parity == 'even' else [0.,1.]

    def rhs(z, y):
        fa, pa, fb, pb = y[:4]
        V = a/np.cos(z)**2
        K = pa*pb-wa*wb*fa*fb
        mass = V*fa*fb
        variation = V*(1+2*z*np.tan(z))*fa*fb
        shift = z*(wa*fa*pb+wb*fb*pa)
        return [pa,(V-wa*wa)*fa,pb,(V-wb*wb)*fb,K,mass,variation,shift]

    refinements = []
    for tolerance in [2e-9, 2e-11]:
        sol = solve_ivp(rhs,(0,A),init+init+[0.]*4,rtol=tolerance,atol=tolerance/50,method='DOP853')
        assert sol.success
        fa,pa,fb,pb = sol.y[:4,-1]
        K,M,Q,I = 2*sol.y[4:,-1]
        B = 2*A*(pa*pb-(a/epsilon**2+wa*wb)*fa*fb)
        norm = 4*np.sqrt(wa*wb)
        vertex = ((1-kap)*K+M+kap*Q+kap*(wa+wb)*I)/norm
        pressure = (K+M-kap*B)/norm
        error = abs(vertex-pressure)/(1+abs(vertex)+abs(pressure))
        assert error < 3e-8, (epsilon,a,parity,error)
        refinements.append(dict(tolerance=tolerance,vertex=vertex,pressure=pressure,normalized_error=error,evaluations=sol.nfev))
    return dict(epsilon=epsilon,a=a,parity=parity,frequencies=frequencies,mouth_potential=a/epsilon**2,refinements=refinements)


def profile_coefficient(T, D, r, h, n):
    L=T+D
    k=2*np.pi*n/L
    # Spatial representation of the periodic convolution, independent of
    # its Fourier multiplier. Split at both edges to resolve narrow mouths.
    def profile(s):
        out=0.
        for j in range(-2,3):
            y=s-j*L
            out += (erf((y+T/2)/(np.sqrt(2)*h))-erf((y-T/2)/(np.sqrt(2)*h)))/2
            out += r/(np.sqrt(2*np.pi)*h)*(np.exp(-(y-T/2)**2/(2*h*h))+np.exp(-(y+T/2)**2/(2*h*h)))
        return out
    value=quad(lambda s: profile(s)*np.cos(k*s),-L/2,L/2,points=[-T/2,T/2],epsabs=1e-11,epsrel=1e-11,limit=400)[0]/L
    expected=2/L*(np.sin(k*T/2)/k+r*np.cos(k*T/2))*np.exp(-h*h*k*k/2)
    assert abs(value-expected)<2e-11
    return dict(T=T,D=D,r=r,h=h,n=n,error=abs(value-expected))


def complex_pair(z):
    return [float(z.real),float(z.imag)]


def sharp(T,D,w):
    L=T+D
    nu=2*np.pi/L
    return 2/(3*np.pi*L)*(-nu**2*T*D/8-(w*w-nu**2)*L*np.sin(w*T/2)*np.sin(w*D/2)/(4*w*np.sin(w*L/2)))


def finite_part(T,D,r,w):
    L=T+D
    nu=2*np.pi/L
    U=w*L/2
    Z=w*(L/2-T)
    Q1=(L-2*T)/4
    Q2=(Q1-L/4*np.sin(Z)/np.sin(U))/(w*w)
    Q3=L*(np.cos(U)+np.cos(Z))/(8*w*np.sin(U))-1/(2*w*w)
    return sharp(T,D,w)+2*r*w*w/(3*np.pi*L)*(-Q1+(w*w-nu*nu)*Q2)+r*r*w*w/(3*np.pi*L)+2*r*r*w*w*(w*w-nu*nu)*Q3/(3*np.pi*L)


def spectral(T,D,r,h,w,cut=8.):
    L=T+D
    nu=2*np.pi/L
    n=np.arange(1,int(np.ceil(cut/(h*nu)))+1,dtype=float)
    k=nu*n
    sine=np.sin(k*T/2)
    cosine=np.cos(k*T/2)
    heat=np.exp(-h*h*k*k)
    weights=2*(n*n-1)/(3*L*L*n)*(sine+r*k*cosine)**2*heat
    direct=np.sum(2*w*w*weights/(k*(w*w-k*k)))
    base=np.sum(2*w*w/(3*np.pi*L)*(k*k-nu*nu)*sine*sine*heat/(k*k*(w*w-k*k)))
    cross=2*r*w*w/(3*np.pi*L)*np.sum((-1/k+(w*w-nu*nu)/(k*(w*w-k*k)))*np.sin(k*T)*heat)
    S=np.sum(cosine*cosine*heat)
    endpoint=-w*w*2*r*r/(3*np.pi*L)*S+2*r*r*w*w*(w*w-nu*nu)/(3*np.pi*L)*np.sum(cosine*cosine*heat/(w*w-k*k))
    decomposition_error=abs(direct-base-cross-endpoint)
    assert decomposition_error<1e-10*(1+abs(direct))
    # Poisson summation includes images approaching when D is tiny.
    images=np.arange(-4,5,dtype=float)*L
    ps=L/(8*np.sqrt(np.pi)*h)*(np.exp(-images**2/(4*h*h)).sum()+np.exp(-(T-images)**2/(4*h*h)).sum())-.5
    poisson_error=abs(S-ps)/(1+abs(ps))
    assert poisson_error<2e-11
    return direct,dict(modes=len(n),decomposition_error=decomposition_error,poisson_normalized_error=poisson_error)


def run():
    start=time.perf_counter()
    syms=symbolic()
    ward=[ward_case(e,a,p,f) for e in [.03,.005,.0005] for a in [.003,.05]
          for p,f in [('even',(.7,1.3)),('odd',(1.1,2.2))]]
    profiles=[profile_coefficient(3.,D,.035,h,n) for D in [.04,2.,9.]
              for h in [.012,.08] for n in [2,7,19]]
    uv=[]
    for D in [.1,1.,7.]:
        for w in [1.1+.2j,2.1+.07j]:
            r=.05
            T=2*np.arccos(r/np.sqrt(1+r*r))
            F=finite_part(T,D,r,w)
            rows=[]
            for h in [.008,.004,.002,.001,.0005]:
                value,meta=spectral(T,D,r,h,w)
                leading=-r*r*w*w/(12*np.pi**1.5*h)
                err=abs(value-leading-F)
                rows.append(dict(h=h,value=complex_pair(value),finite_part_error=err,leading=complex_pair(leading),**meta))
            for previous,current in zip(rows,rows[1:]):
                assert current['finite_part_error'] < .7*previous['finite_part_error']
            slope=(complex(*rows[-1]['value'])-complex(*rows[-2]['value']))/(complex(*rows[-1]['leading'])-complex(*rows[-2]['leading']))
            assert abs(slope-1)<.003
            uv.append(dict(D=D,w=complex_pair(w),r=r,T=T,finite_part=complex_pair(F),slope_ratio=complex_pair(slope),refinements=rows))
    uniform=[]
    for D in [.0001,.01,.1,1.,10.,100.]:
        for h,r in [(.005,.003),(.02,.02),(.005,.05)]:
            for eta in [.03,.2]:
                T=3.
                w=1.3+1j*eta
                value,meta=spectral(T,D,r,h,w)
                refined,_=spectral(T,D,r,h,w,cut=10.)
                diff=abs(value-sharp(T,D,w))
                expression=h+h*h/eta+r*(1+1/eta)+r*r*(1/h+1+1/eta)
                tail=abs(value-refined)
                assert tail<1e-11*(1+abs(value))
                # A diagnostic ratio, not a fitted constant in the theorem.
                uniform.append(dict(D=D,h=h,r=r,eta=eta,error=diff,bound_expression=expression,ratio=diff/expression,tail_refinement_error=tail,**meta))
    scales=[]
    for q in [100,10000,1000000,100000000]:
        epsilon=q**-2
        delta=np.log(q)**-2
        r=epsilon/np.sqrt(1-epsilon**2)
        h=epsilon
        bound=(h+h*h/delta+r*(1+1/delta)+r*r*(1/h+1+1/delta))/delta
        scales.append(dict(q=q,epsilon=epsilon,delta=delta,normalized_sufficient_expression=bound,over_sharp_r2_over_h=r*r/epsilon**4))
    assert all(b['normalized_sufficient_expression']<a['normalized_sufficient_expression'] for a,b in zip(scales,scales[1:]))
    dependencies=['notes/moving-mouth-clock-vertex.md','notes/proper-clock-anomaly.md',
                  'notes/collar-conservative-response.md','notes/transparent-clock-flux.md',
                  'notes/finite-mouth-determinant.md','notes/singlet-complex-response.md']
    result=dict(date='2026-10-11',status='passed',scope='Product Gaussian linear vertex and specified smoothed free CFT spectral response',
                environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__),
                symbolic=syms,ward=ward,profiles=profiles,ultraviolet=uv,uniform=uniform,scaling=scales,
                product_gaussian_linear_vertex_derived=True,free_smoothing_bound_proved_in_note=True,
                physical_full_vertex_derived=False,interacting_pressure_response_controlled=False,
                physical_damping_proved=False,novelty_established=False,lean_invocations=0,
                script_sha256=digest(Path(__file__)),input_sha256={p:digest(ROOT/p) for p in dependencies},
                elapsed_seconds=time.perf_counter()-start)
    (ROOT/'evidence/moving-mouth-clock-vertex.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],seconds=result['elapsed_seconds'],symbolic=len(syms),ward=len(ward),profiles=len(profiles),uv_families=len(uv),uniform=len(uniform),max_ward_error=max(r['refinements'][-1]['normalized_error'] for r in ward),max_uv_slope_error=max(abs(complex(*r['slope_ratio'])-1) for r in uv),max_uniform_ratio=max(r['ratio'] for r in uniform))))


if __name__ == '__main__':
    run()
