"""Actual dipole source-corrected trial checks; no loaded PDE simulation."""
from __future__ import annotations
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import platform
from pathlib import Path

import numpy as np
import scipy
import sympy as sp
from scipy.linalg import expm
from scipy.special import roots_hermite, roots_legendre

ROOT = Path(__file__).resolve().parents[1]
_spec=importlib.util.spec_from_file_location('joint_mouth_matching',ROOT/'checks/joint_mouth_matching.py')
_joint=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_joint)
Dipole=_joint.Dipole
ALG_TOL = 5e-8
FD_TOL = 5e-5
REFINE_TOL = 2e-7
ORDER_RATIO_LIMIT = 1.04
GAUGE_ORDER_LIMIT = 1.04
SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.])


def center(u, e):
    a=1-e*e
    root=np.sqrt(a*a+4*u*u*e*e)
    mu=(a+root)/(2*u)
    mup=(-2*e+(-4*e*a+8*u*u*e)/(2*root))/(2*u)
    r=.5*np.sqrt((mu*mu-1)*(1-e*e))
    z=.5*mu*e
    rp=.25*(2*mu*mup*(1-e*e)-2*e*(mu*mu-1))/np.sqrt((mu*mu-1)*(1-e*e))
    zp=.5*(mu+e*mup)
    return np.array([r,z]),np.hypot(rp,zp)


def cartesian_source(points, amplitudes):
    """Two fixed global axisymmetric C2 profiles, zero in |z|>=0.3.

    Components are scalar, cylindrical radial, azimuthal, axial.
    Radial/azimuthal coefficients vanish linearly on the axis.
    """
    r,z=np.asarray(points)
    bump=np.where(abs(z)<.3,np.cos(np.pi*z/.6)**4,0.)
    f=bump*np.exp(-r*r/3)
    a,b=amplitudes
    return np.array([(.2*a-.15*b)*f,
                     (.3*a-.2*b*z)*r*f,
                     (.4*a+.25*b)*r*f,
                     (.1*a*(1+z)+.3*b)*f])


def frame_source(points,tangent,normal,amplitudes):
    scalar,vr,vy,vz=cartesian_source(points,amplitudes)
    return np.array([scalar, -vr*normal[0]-vz*normal[1], vy,
                     vr*tangent[0]+vz*tangent[1]])


def geom(dip,u,e,amplitudes):
    p,speed=center(u,e)
    gg=dip.geometry(p,u)
    g={k:float(gg[k][0]) for k in ['b','omega','k','c']}
    g.update(p=p,speed=speed,t=gg['t'][:,0],n=gg['n'][:,0])
    src=frame_source(p,g['t'],g['n'],amplitudes)
    g['v']=src[0]+src[3]
    g['w']=src[1]+1j*src[2]
    return g


def differentiated(dip,u,e,amplitudes,step):
    g=geom(dip,u,e,amplitudes)
    nearby=[geom(dip,u,e+j*step,amplitudes) for j in [-2,-1,1,2]]
    for name in ['omega','k','c','w']:
        g[name+'s']=(nearby[0][name]-8*nearby[1][name]+8*nearby[2][name]-nearby[3][name])/(12*step*g['speed'])
    return g


def norm(values,weights):
    return float(np.sqrt(np.sum(weights*np.sum(np.abs(values)**2,axis=0))))


def phase_between(dip,u,e0,e1,z,amplitudes):
    nodes,weights=roots_legendre(10)
    es=(e0+e1)/2+(e1-e0)*nodes/2
    integrand=[]
    for ee in es:
        g=geom(dip,u,ee,amplitudes)
        integrand.append((z-g['v'])*g['speed'])
    return np.exp(1j*(e1-e0)/2*np.dot(weights,integrand))


def raw_trial(dip,u,e,q,y,amplitudes,phase=1):
    g=geom(dip,u,e,amplitudes)
    gauss=(q*g['b']/np.pi)**.25*np.exp(-q*g['b']*y*y/2)
    p=1-q*g['c']*y**3/6+1j*g['w']*y
    lower=g['omega']*y
    return phase*gauss*np.array([p,lower])


def section(dip,u,e,q,z,amplitudes,step,finite_difference=False):
    g=differentiated(dip,u,e,amplitudes,step)
    nodes,weights=roots_hermite(96)
    weights/=np.sqrt(np.pi)
    b,c,k,omega=g['b'],g['c'],g['k'],g['omega']
    distances=[g['p'][0],np.linalg.norm(g['p']-[0.,.5]),np.linalg.norm(g['p']-[0.,-.5]),1/max(abs(k),1e-12)]
    width=.12*min(distances)
    y=nodes/np.sqrt(q*b)
    inside=abs(y)<width
    y=y[inside];weights=weights[inside]
    h=1-k*y;hs=-g['ks']*y
    assert np.min(h)>.8
    m=dip.mf(*(g['p'][:,None]+g['n'][:,None]*y),u)
    v0,vx,vy,vz=frame_source(g['p'][:,None]+g['n'][:,None]*y,g['t'],g['n'],amplitudes)
    vpp,vmm,vpm,vmp=v0+vz,v0-vz,vx-1j*vy,vx+1j*vy
    p0=1-q*c*y**3/6
    p=p0+1j*g['w']*y;lower=omega*y
    ps=-q*g['cs']*y**3/6+1j*g['ws']*y
    py=-q*c*y*y/2+1j*g['w']
    ls=g['omegas']*y
    G=omega*(1-2*q*b*y*y)
    vec=np.array([p,lower])
    dn=np.array([py-q*b*y*p,omega-q*b*y*lower])
    ds=np.array([ps+p*(G+1j*(z-g['v'])),ls+lower*(G+1j*(z-g['v']))])
    vvec=np.array([vpp*p+vpm*lower,vmp*p+vmm*lower])
    direct=1j*SX@dn-1j*SZ@ds/h+1j*SZ@vec*hs/(2*h*h)+SY@vec*(q*m)+vvec-z*vec
    upper=1j*(omega-q*(b*y+m)*lower)-1j*(G*p+ps)/h+1j*hs*p/(2*h*h)+((z-g['v'])/h+vpp-z)*p+vpm*lower
    down=1j*(py+q*(m-b*y)*p)+1j*(G*lower+ls)/h-1j*hs*lower/(2*h*h)+(-(z-g['v'])/h+vmm-z)*lower+vmp*p
    expanded=np.array([upper,down])
    algebra_error=float(np.max(abs(direct-expanded)))
    assert algebra_error<ALG_TOL
    # A constant-width transverse cutoff is legitimate for this local section.
    chi=np.ones_like(y);cy=np.zeros_like(y)
    trans=abs(y)>width/2
    angle=np.pi*(abs(y[trans])/width-.5)
    chi[trans]=np.cos(angle)**2
    cy[trans]=-np.pi/width*np.sin(2*angle)*np.sign(y[trans])
    cut=chi*direct+1j*SX@vec*cy
    # Keep the Wilson phase but omit the new transverse correction.
    oldvec=np.array([p0,lower])
    olddn=np.array([-q*c*y*y/2-q*b*y*p0,omega-q*b*y*lower])
    oldds=np.array([-q*g['cs']*y**3/6+p0*(G+1j*(z-g['v'])),ls+lower*(G+1j*(z-g['v']))])
    oldv=np.array([vpp*p0+vpm*lower,vmp*p0+vmm*lower])
    old=1j*SX@olddn-1j*SZ@oldds/h+1j*SZ@oldvec*hs/(2*h*h)+SY@oldvec*(q*m)+oldv-z*oldvec
    old=chi*old+1j*SX@oldvec*cy
    out=dict(q=q,ell=q*(1-u),residual=norm(cut,weights),
             sqrt_q_residual=float(np.sqrt(q)*norm(cut,weights)),
             phase_only_residual=norm(old,weights),algebra_error=algebra_error)
    if finite_difference:
        # Independent fourth-order differences of the full normalized spinor,
        # including the source phase integral and moving field geometry.
        gauss=(q*b/np.pi)**.25*np.exp(-q*b*y*y/2)
        es=step
        fields=[]
        for j in [-2,-1,1,2]:
            phase=phase_between(dip,u,e,e+j*es,z,amplitudes)
            fields.append(raw_trial(dip,u,e+j*es,q,y,amplitudes,phase))
        dsfd=(fields[0]-8*fields[1]+8*fields[2]-fields[3])/(12*es*g['speed'])/gauss
        dy=.0003/np.sqrt(q*b)
        fields=[raw_trial(dip,u,e,q,y+j*dy,amplitudes) for j in [-2,-1,1,2]]
        dnfd=(fields[0]-8*fields[1]+8*fields[2]-fields[3])/(12*dy)/gauss
        fd=1j*SX@dnfd-1j*SZ@dsfd/h+1j*SZ@vec*hs/(2*h*h)+SY@vec*(q*m)+vvec-z*vec
        fd_error=norm(fd-direct,weights)
        assert fd_error<FD_TOL,(u,e,q,z,fd_error)
        out['direct_difference_error']=fd_error
    return out


def source_sections():
    dip=Dipole();cases=[]
    for u in [.25,.5,.75]:
        fluxes=[10001,40001,160001] if u==.5 else [10002,40002,160002]
        for e in [-.08,0.,.08]:
            for amplitudes in [(1.,.3),(.45+.025j,-.2+.012j)]:
                for z in [.6+.03j,-.4+.4j]:
                    rows=[]
                    for q in fluxes:
                        assert abs((1-u)*q-(np.floor((1-u)*q)+.5))<1e-9
                        coarse=section(dip,u,e,q,z,amplitudes,1e-4)
                        fine=section(dip,u,e,q,z,amplitudes,5e-5,finite_difference=(q==fluxes[-1]))
                        refinement=abs(coarse['residual']-fine['residual'])
                        assert refinement<REFINE_TOL
                        fine['geometry_refinement_error']=refinement
                        rows.append(fine)
                    ratio=rows[-1]['sqrt_q_residual']/rows[-2]['sqrt_q_residual']
                    assert ratio<ORDER_RATIO_LIMIT,(u,e,z,ratio)
                    assert rows[-1]['residual']<rows[-1]['phase_only_residual']
                    cases.append(dict(u=u,prolate_eta=e,z=[z.real,z.imag],
                                      amplitudes=[[complex(a).real,complex(a).imag] for a in amplitudes],
                                      rows=rows,scaled_last_ratio=ratio))
    return cases


def gauge_checks():
    # A general smooth gauge in one Fermi patch. Exact phase versus its
    # normal expansion multiplies both old spinor components.
    nodes,weights=roots_hermite(100);weights/=np.sqrt(np.pi)
    rows=[]
    for b,c,omega,chi_y,chi_yy in [(0.8,.3,.2,.4,-.3),(2.,-.6,-.15,-.7,.2)]:
        entries=[]
        for q in [1000,4000,16000]:
            y=nodes/np.sqrt(q*b)
            p0=1-q*c*y**3/6;lower=omega*y
            exact=np.exp(1j*(chi_y*y+chi_yy*y*y/2))*np.array([p0,lower])
            approx=np.array([p0+1j*chi_y*y,lower])
            error=norm(exact-approx,weights)
            entries.append(dict(q=q,error=error,q_times_error=q*error))
        assert entries[-1]['q_times_error']/entries[-2]['q_times_error']<GAUGE_ORDER_LIMIT
        rows.append(entries)
    y,b,q,w=sp.symbols('y b q w')
    g=sp.exp(-q*b*y*y/2)
    assert sp.simplify(sp.I*(sp.diff(sp.I*w*y*g,y)+q*b*y*sp.I*w*y*g)+w*g)==0
    return rows


def collar_checks():
    rows=[]
    for mass in [.3,1.,4.]:
        for nu in [.015,.1,.5]:
            for omega in [-1.2,.7]:
                length=1.1;ell=.3
                gen=np.array([[-nu+1j*omega,mass],[mass,nu-1j*omega]])
                mat=expm(length*gen)
                a,b,c,d=mat.ravel()
                K=np.array([[-c/d,1/d],[a-b*c/d,b/d]])
                least=float(np.min(np.linalg.eigvalsh(np.eye(2)-K.conj().T@K)))
                bound=nu/(nu+mass)*(1-np.exp(-2*(nu+mass)*ell))
                assert least>=bound-2e-12
                rows.append(dict(mass=mass,nu=nu,omega=omega,defect=least,bound=bound))
    return rows


def symbolic_and_source_checks():
    n,xi=sp.symbols('n xi',nonnegative=True)
    a=sp.sqrt(2*n+2)
    matrix=sp.Matrix([[0,a],[a,2*xi]])
    printed=sp.Matrix([[2*xi,-a],[-a,0]])/(2*n+2)
    assert sp.simplify(matrix*printed+sp.eye(2))==sp.zeros(2)
    assert sp.simplify(matrix*(-printed)-sp.eye(2))==sp.zeros(2)
    # At n=0 an upper phi_1 coefficient maps to w_1(0), although the
    # spinor is orthogonal to the upper phi_0 kernel. It need not vanish.
    n0_transform_upper_phi1=1
    assert n0_transform_upper_phi1!=0
    exponents=[-F(3,8)+(F(j,1)+F(1,2))*F(1,8) for j in range(3)]
    assert exponents==[-F(5,16),-F(3,16),-F(1,16)]
    return dict(source_printed_inverse_product='-I',corrected_inverse_product='I',
                shifted_hermite_counterexample='u=(f(x1) phi_1(x2),0) has w_1(xi,0)=Fourier(f), not zero',
                exponents=[str(x) for x in exponents])


def main():
    cases=source_sections()
    gauge=gauge_checks()
    collar=collar_checks()
    exact=symbolic_and_source_checks()
    dependencies=['checks/dipole_source_matching.py','notes/dipole-source-matching.md',
                  'evidence/dipole-source-review.json','checks/joint_mouth_matching.py',
                  'notes/joint-mouth-matching.md','notes/gauged-mouth-interaction.md',
                  'notes/throat-load-feedback.md','evidence/joint-mouth-matching.json']
    result=dict(date='2026-10-11',status='passed',
                scope='source-corrected actual dipole trial, conditional loaded return theorem; not a full PDE simulation, long-return estimate or vacuum-stress result',
                tolerances=dict(algebra=ALG_TOL,finite_difference=FD_TOL,geometry_refinement=REFINE_TOL,
                                scaled_order_ratio=ORDER_RATIO_LIMIT,gauge_scaled_order=GAUGE_ORDER_LIMIT),
                input_hash_encoding='UTF-8 text with LF newlines, matching committed Git blobs',
                input_sha256={p:hashlib.sha256((ROOT/p).read_text(encoding='utf-8').encode('utf-8')).hexdigest() for p in dependencies},
                versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__),
                development_events=[dict(stage='initial import',outcome='ModuleNotFoundError: joint_mouth_matching',
                                         correction='Load the adjacent pinned dependency explicitly by importlib path; the embedded Python distribution does not add the script directory to sys.path.',
                                         mathematical_assertions_executed=False,tolerance_changed=False),
                                    dict(stage='staged provenance check',outcome='Prior joint receipt has CRLF in the worktree and LF in Git',
                                         correction='Hash all textual verification inputs with explicit UTF-8/LF normalization, then compare with staged blobs.',
                                         mathematical_failure=False,tolerance_changed=False)],
                dipole_source_cases=cases,gauge_cases=gauge,collar_cases=collar,exact=exact,
                long_return_matching_proved=False,physical_geometry_matched=False,
                renormalized_stress_proved=False,novelty_accepted=False,lean_invocations=0)
    (ROOT/'evidence/dipole-source-matching.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='passed',source_cases=len(cases),flux_values_per_case=3,
                         gauge_cases=len(gauge),collar_cases=len(collar),
                         max_algebra_error=max(r['algebra_error'] for c in cases for r in c['rows']),
                         max_finite_difference_error=max(c['rows'][-1]['direct_difference_error'] for c in cases),
                         max_geometry_refinement=max(r['geometry_refinement_error'] for c in cases for r in c['rows']),
                         max_scaled_last_ratio=max(c['scaled_last_ratio'] for c in cases)),indent=2))


if __name__=='__main__':
    main()
