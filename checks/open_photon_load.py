"""Independent checks of a photon core with nonnegative exterior loads.

No Lean or cloud invocation. Finite samples supplement the ordinary proof.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from pathlib import Path
import hashlib
import importlib.util
import json
import platform
import numpy as np
import scipy
from scipy.linalg import block_diag
from scipy.integrate import solve_ivp
from scipy.special import spherical_jn, spherical_yn
import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / 'notes/open-photon-load.md'
OUT = ROOT / 'evidence/open-photon-load.json'
DEP = ROOT / 'checks/propagating_photon_polarization.py'
spec = importlib.util.spec_from_file_location('previous_photon', DEP)
previous = importlib.util.module_from_spec(spec)
spec.loader.exec_module(previous)
ODE_TOLERANCE = .001
PROBE_TOLERANCE = .01


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def spectral_factor(w):
    v, eta = abs(w.real), w.imag
    assert eta > 0
    return 1. if v <= eta else abs(w)**2 / (2*v*eta)


def bound(w, xmin=1.):
    return 1 + max(1., abs(w)**2/(2*xmin*xmin))*spectral_factor(w)


def exact_checks():
    lam, v, eta = sy.symbols('lambda v eta', positive=True)
    f = lam**2 / ((lam-v*v+eta*eta)**2+4*v*v*eta*eta)
    lm = (v*v+eta*eta)**2/(v*v-eta*eta)
    assert sy.factor(sy.diff(f, lam).subs(lam, lm)) == 0
    assert sy.factor(f.subs(lam, lm)-(v*v+eta*eta)**2/(4*v*v*eta*eta)) == 0
    r, w = sy.symbols('r w', nonzero=True)
    u = sy.exp(sy.I*w*r)*(1+sy.I/(w*r))
    assert sy.simplify(-sy.diff(u, r, 2)+(2/r**2-w*w)*u) == 0
    zeta = sy.I*w-sy.I/(r*(w*r+sy.I))
    assert sy.simplify(sy.diff(u, r)/u-zeta) == 0
    a, b = sy.symbols('a b', positive=True)
    assert sy.expand(-b+1+2*a).subs({a:sy.Rational(1,16), b:2}) == -sy.Rational(7,8)
    assert sy.expand(-b+1+3*a).subs({a:sy.Rational(1,16), b:2}) == -sy.Rational(13,16)
    return ['spectral_maximum_stationarity', 'spectral_maximum_value',
            'outgoing_l1_wave_equation', 'outgoing_l1_log_derivative',
            'epsilon_power', 'absolute_correction_power']


def matrix_checks():
    """Full resolvent versus spectral decomposition and exterior elimination."""
    rng = np.random.default_rng(531011)
    rows = []
    for nc, ne in [(3, 5), (5, 11), (8, 19)]:
        G = rng.normal(size=(nc+1, nc))
        # Positive couplings include arbitrary angular/port mixing.
        springs = rng.normal(size=(ne+3, nc+ne))
        springs *= .13
        H = block_diag(G.T@G+2*np.eye(nc), np.diag(np.geomspace(1e-8, .4, ne)))
        H += springs.T@springs
        Hcore = block_diag(G.T@G+2*np.eye(nc), np.zeros((ne, ne)))
        assert np.linalg.eigvalsh(H-Hcore)[0] > -2e-13
        E, V = np.linalg.eigh(H)
        assert E[0] > 0 and E[0] < .1
        for w in [1+.2j, .7+.02j, -.7+.02j, .4j, 2.5+.1j]:
            K = H-w*w*np.eye(nc+ne)
            R = np.linalg.inv(K)
            Rs = (V/(E-w*w)[None, :])@V.T
            Rc = np.linalg.inv(K[:nc,:nc]-K[:nc,nc:]@np.linalg.solve(K[nc:,nc:],K[nc:,:nc]))
            err = max(np.linalg.norm(R-Rs,2),np.linalg.norm(R[:nc,:nc]-Rc,2))
            assert err < 2e-10
            B = np.hstack([G.T, w*np.eye(nc)])
            L = np.vstack([G, -w*np.eye(nc)])
            D = np.eye(2*nc+1)-L@Rc@B
            val = float(np.linalg.norm(D,2))
            assert val <= bound(w)*(1+2e-12)
            actual_r = float(np.max(E/abs(E-w*w)))
            assert actual_r <= spectral_factor(w)*(1+2e-12)
            rows.append(dict(core_size=nc, exterior_size=ne, w=[w.real,w.imag],
                lowest_eigenvalue=float(E[0]), reduction_error=float(err),
                weighted_norm=val, bound=bound(w), spectral_ratio=actual_r/spectral_factor(w)))
    return rows


def resonance_checks():
    H = np.array([[3.,-1.],[-1.,1.]])
    lam = 2-np.sqrt(2.)
    root = np.sqrt(lam)
    weight = lam/4
    limiting = root*weight/2
    rows=[]
    for eta in [1e-1,1e-2,1e-3,1e-4,1e-5]:
        w=root+1j*eta
        R=np.linalg.inv(H-w*w*np.eye(2))[0,0]
        response=abs(1+w*w*R)
        assert response <= bound(w)
        rows.append(dict(eta=eta, response=float(response),
                         eta_times_response=float(eta*response),bound=bound(w)))
    assert abs(rows[-1]['eta_times_response']/limiting-1) < 1e-6
    assert rows[-1]['response'] > 1000
    return dict(eigenvalue=float(lam),core_weight=float(weight),
                eta_times_response_limit=float(limiting),cases=rows)


def zeta(l,w,rb):
    z=w*rb
    hankel=spherical_jn(l,z)+1j*spherical_yn(l,z)
    hp=spherical_jn(l,z,True)+1j*spherical_yn(l,z,True)
    return 1/rb+w*hp/hankel


def core_matrices(n,a,l):
    """Piecewise linear scalar, orthonormal constant potential sources."""
    h=2*a/n
    nodes=np.linspace(-a,a,n+1)
    centers=(nodes[:-1]+nodes[1:])/2
    G=np.zeros((n,n+1)); J=np.zeros_like(G)
    for i in range(n):
        G[i,i:i+2]=[-1/np.sqrt(h),1/np.sqrt(h)]
        J[i,i:i+2]=np.sqrt(h)/2
    mass=np.zeros((n+1,n+1)); pot=np.zeros_like(mass)
    source=np.zeros(n+1)
    gx,gw=np.polynomial.legendre.leggauss(6)
    shape=np.array([(1-gx)/2,(1+gx)/2])
    for i,c in enumerate(centers):
        p=c+h*gx/2
        measure=h*gw/2
        mass[i:i+2,i:i+2]+=(shape*measure)@shape.T
        V=l*(l+1)/np.cos(p)**2
        pot[i:i+2,i:i+2]+=(shape*(measure*V))@shape.T
        source[i:i+2]+=shape@(measure*np.cos(np.pi*p/(2*a))/np.sqrt(a))
    return nodes, centers, h, G, J, mass, G.T@G+pot,source


def open_core(n,a,l,w):
    nodes,centers,h,G,J,mass,H,source=core_matrices(n,a,l)
    rb=np.cos(a)
    load=zeta(l,w,rb)
    K=(H-w*w*mass).astype(complex)
    K[0,0]-=load; K[-1,-1]-=load
    B=np.hstack([G.T,w*J.T]); L=np.vstack([G,-w*J])
    U=np.r_[np.ones(n),-1j*np.ones(n)]
    weighted=(np.eye(2*n)-L@np.linalg.solve(K,B))/U[:,None]/U[None,:]
    norm=float(np.linalg.norm(weighted,2))
    assert norm <= bound(w)*(1+2e-11)
    probe=source@np.linalg.solve(K,source)
    return weighted,dict(n=n,a=a,l=l,w=[w.real,w.imag],weighted_norm=norm,
        bound=bound(w),scalar_probe=[float(probe.real),float(probe.imag)]),probe


def ode_probe(a,l,w):
    """Continuum shooting, independent of core Galerkin/Schur matrices."""
    load=zeta(l,w,np.cos(a))
    def rhs(s,y):
        f=np.cos(np.pi*s/(2*a))/np.sqrt(a)
        pot=l*(l+1)/np.cos(s)**2-w*w
        u,du,p,dp,ih,ip=y
        return [du,pot*u,dp,pot*p-f,f*u,f*p]
    result=solve_ivp(rhs,[-a,a],np.array([1,-load,0,0,0,0],dtype=complex),
                     rtol=2e-12,atol=2e-13,method='DOP853')
    assert result.success
    u,du,p,dp,ih,ip=result.y[:,-1]
    amplitude=-(dp-load*p)/(du-load*u)
    residual=abs(amplitude*(du-load*u)+(dp-load*p))
    assert residual < 1e-9
    return amplitude*ih+ip,residual


def radiative_checks():
    rows=[]
    for a in [.7,1.3,1.5]:
        for l in [1,2,4]:
            for w in [1+.15j,.4j]:
                _,row,_=open_core(24,a,l,w)
                if w.real==0:
                    load=zeta(l,w,np.cos(a))
                    assert abs(load.imag)<1e-10 and load.real<0
                if l==1:
                    rb=np.cos(a)
                    formula=1j*w-1j/(rb*(w*rb+1j))
                    assert abs(formula-zeta(l,w,rb))<2e-12
                rows.append(row)
    refinement=[]
    a,l,w=1.3,1,1+.15j
    reference,residual=ode_probe(a,l,w)
    for n in [12,24,48,96,192]:
        _,row,value=open_core(n,a,l,w)
        row['continuum_ode_relative_error']=float(abs(value-reference)/abs(reference))
        refinement.append(row)
    assert refinement[-1]['continuum_ode_relative_error']<ODE_TOLERANCE
    assert all(refinement[i]['continuum_ode_relative_error'] < refinement[i-1]['continuum_ode_relative_error'] for i in range(1,len(refinement)))
    return rows,dict(reference=[reference.real,reference.imag],boundary_residual=float(residual),
                     tolerance=ODE_TOLERANCE,cases=refinement)


def coupled(q,n,cut,mixed):
    """Exact half-line loads, mixed canonical dual fields, all LLL vertices."""
    a,w,g2=1.3,1+.15j,1e-4
    E,degrees,sigmas,_,gram=previous.angular.angular_basis(q)
    ls=degrees+1; kappas=ls*(ls+1); taus=sigmas*kappas
    Ts=E*np.sqrt(taus)[:,None,None]
    count=q*q-1
    data={l:core_matrices(n,a,l) for l in range(1,q)}
    nodes,centers,h,G,J,_,_,source=data[1]
    K=block_diag(*[data[int(l)][6]-w*w*data[int(l)][5] for l in ls]).astype(complex)
    rng=np.random.default_rng(11053+q)
    rotation=np.linalg.qr(rng.normal(size=(count,count)))[0] if mixed else np.eye(count)
    Z=rotation@np.diag([zeta(int(l),w,np.cos(a)) for l in ls])@rotation.T
    for end in [0,n]:
        ids=np.arange(count)*(n+1)+end
        K[np.ix_(ids,ids)]-=Z
    B=block_diag(*[np.hstack([G.T,w*J.T])]*count)
    L=block_diag(*[np.vstack([G,-w*J])]*count)
    U=np.tile(np.r_[np.ones(n),-1j*np.ones(n)],count)
    size=2*n
    mhalf=np.sqrt(np.repeat(kappas,size)); phalf=np.sqrt(np.repeat(taus,size))
    weighted=(np.eye(count*size)-L@np.linalg.solve(K,B))/U[:,None]/U[None,:]
    assert np.linalg.norm(weighted,2)<=bound(w)*(1+2e-11)
    D=weighted/mhalf[:,None]/mhalf[None,:]
    lengths=2*a+np.array([.3,1.2,4.7][:q])
    Sig=np.zeros((count*size,count*size),dtype=complex)
    for i in range(q):
        for j in range(q):
            tij=Ts[:,i,j]
            pair=previous.pair_current(lengths[i],lengths[j],cut,w,centers,centers,h)
            Sig+=np.kron(tij[:,None]*tij.conj()[None,:],pair)
    contact=np.tile(np.r_[np.zeros(n),np.ones(n)],count)
    Sig+=np.diag(np.repeat(taus/np.pi,size)*contact)
    _,AT=previous.frame_constant(2*a,2*a)
    BJ=4*AT*(1+abs(w)/w.imag)+1/np.pi
    epsilon=g2*bound(w)*BJ*max(taus/kappas)
    assert epsilon<1
    Sweighted=Sig/mhalf[:,None]/mhalf[None,:]
    # Invert I-g^2 D Sigma, not D itself; the theorem does not need that inverse.
    dressed=np.linalg.solve(np.eye(count*size)-g2*weighted@Sweighted,weighted)
    delta=dressed-weighted
    residual=np.linalg.norm(dressed-weighted-g2*weighted@Sweighted@dressed,2)
    assert residual<2e-12
    actual=float(np.linalg.norm(delta,2)); predicted=bound(w)*epsilon/(1-epsilon)
    assert actual<=predicted
    Dg=dressed/mhalf[:,None]/mhalf[None,:]
    current=g2*Sig@Dg@Sig/phalf[:,None]/phalf[None,:]
    assert np.linalg.norm(current,2)<=BJ*epsilon/(1-epsilon)
    probe=np.zeros(count*size)
    probe[n:2*n]=np.sqrt(h)*np.cos(np.pi*centers/(2*a))*np.sinc(h/(4*a))/np.sqrt(a)
    value=probe@delta@probe
    mixing=float(np.linalg.norm(Z-np.diag(np.diag(Z)),2))
    if mixed and q==3: assert mixing>1e-4
    return dict(q=q,n=n,cutoff=cut,mixed_exterior=mixed,angular_gram_error=gram,
        exterior_angular_mixing_norm=mixing,epsilon=float(epsilon),
        weighted_correction_norm=actual,correction_bound=float(predicted),
        current_correction_norm=float(np.linalg.norm(current,2)),current_bound=float(BJ*epsilon/(1-epsilon)),
        resolvent_identity_error=float(residual),probe=[float(value.real),float(value.imag)])


def run():
    exact=exact_checks(); matrices=matrix_checks(); resonance=resonance_checks()
    radiative,ode=radiative_checks()
    angular=[coupled(q,8,64.,mixed) for q in [2,3] for mixed in [False,True]]
    uv=[coupled(2,48,cut,False) for cut in [16.,32.,64.,128.]]
    spatial=[coupled(2,n,128.,False) for n in [12,24,48,96]]
    def relative(rows):
        values=[complex(*r['probe']) for r in rows]
        return [float(abs(values[i]-values[i-1])/abs(values[i])) for i in range(1,len(values))]
    uv_errors,spatial_errors=relative(uv),relative(spatial)
    assert uv_errors[-1]<PROBE_TOLERANCE,uv_errors
    assert spatial_errors[-1]<PROBE_TOLERANCE,spatial_errors
    result=dict(status='passed_open_photon_load_checks',date='2026-10-11',
        scope='Projected quadratic photon response for nonnegative exterior extensions with local core energy domination',
        note_sha256=sha(NOTE),script_sha256=sha(Path(__file__)),
        dependencies={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in [DEP,previous.DEP]},
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sy.__version__),
        physical_coupled_load_derived=False,renormalized_clock_stress_controlled=False,
        physical_damping_proved=False,novelty_established=False,lean_invocations=0,
        exact_checks=exact,matrix_cases=matrices,resonance=resonance,
        radiative_cases=radiative,continuum_ode=ode,angular_cases=angular,
        ultraviolet_refinement=uv,spatial_refinement=spatial,
        ultraviolet_relative_changes=uv_errors,spatial_relative_changes=spatial_errors,
        quadratic_probe_tolerance=PROBE_TOLERANCE)
    OUT.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],exact_checks=len(exact),matrix_cases=len(matrices),
        radiative_cases=len(radiative),angular_cases=len(angular),
        final_ode_error=ode['cases'][-1]['continuum_ode_relative_error'],
        ultraviolet_relative_changes=uv_errors,spatial_relative_changes=spatial_errors)))


if __name__=='__main__':
    run()
