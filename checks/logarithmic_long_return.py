"""Supplementary global dipole checks at shrinking resonant flux labels.

No loaded PDE solver, stress calculation, or finite-q theorem is claimed.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
from scipy.special import roots_legendre
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('joint', ROOT/'checks/joint_mouth_matching.py')
joint = importlib.util.module_from_spec(spec)
spec.loader.exec_module(joint)


def text_hash(path):
    return hashlib.sha256(path.read_text(encoding='utf-8').replace('\r\n', '\n').encode()).hexdigest()


def stable_flux(r, z):
    """Exact rationalization; no quadrature or far-field truncation."""
    r, z = np.broadcast_arrays(np.asarray(r, dtype=np.longdouble),
                               np.asarray(z, dtype=np.longdouble))
    am, ap = z + np.longdouble('.5'), z - np.longdouble('.5')
    rm, rp = np.sqrt(r*r+am*am), np.sqrt(r*r+ap*ap)
    value = np.empty_like(r)
    same = abs(z) > .5
    value[same] = (z[same]*r[same]**2 /
                   (rm[same]*rp[same]*(am[same]*rp[same]+ap[same]*rm[same])))
    value[~same] = (am[~same]/rm[~same]-ap[~same]/rp[~same])/2
    return value


class StableDipole(joint.Dipole):
    def __init__(self):
        super().__init__()
        self.uf = stable_flux
        self.mf = lambda r,z,u: (stable_flux(r,z)-np.longdouble(u))/np.asarray(r, dtype=np.longdouble)

    def geometry(self, points, u):
        return super().geometry(np.asarray(points, dtype=np.longdouble), u)

    def on_sphere(self,rho,u):
        if rho == 1/(2*u):
            return np.arccos(.5/rho),np.array([np.sqrt(rho*rho-.25),0.])
        return super().on_sphere(rho,u)

    def section(self,rho,U,Q,Xi,step_fraction=1e-4,normal_order=256):
        # Resolve the compact support even when it is narrower than the
        # finite-q Gaussian. Split at both taper junctions.
        theta,p=self.on_sphere(rho,U);g=self.geometry(p,U)
        t=g['t'][:,0];n=g['n'][:,0];step=step_fraction*rho
        shifted=[self.geometry(p+j*step*t,U) for j in [-2,-1,1,2]]
        for name in ['omega','k','c']:
            g[name+'s']=(shifted[0][name]-8*shifted[1][name]+8*shifted[2][name]-shifted[3][name])/(12*step)
        b,c,a,k=[float(g[name][0]) for name in ['b','c','omega','k']]
        cs,ks,aas=[float(g[name][0]) for name in ['cs','ks','omegas']]
        def width(point,curvature):
            R,Z=point
            return .15/np.sqrt(1/R**2+1/(R*R+(Z+.5)**2)+1/(R*R+(Z-.5)**2)+1+curvature**2)
        tube=width(p,k)
        widths=[width(p+j*step*t,float(gg['k'][0])) for j,gg in zip([-2,-1,1,2],shifted)]
        tube_s=(widths[0]-8*widths[1]+8*widths[2]-widths[3])/(12*step)
        radius=min(tube,10/np.sqrt(Q*b))
        edges=sorted(set([-radius,0.,radius]+[v for v in [-tube/2,tube/2] if abs(v)<radius]))
        nodes,weights=roots_legendre(normal_order//4)
        ys=[];ws=[]
        for left,right in zip(edges[:-1],edges[1:]):
            ys.append((left+right)/2+(right-left)/2*nodes)
            ws.append((right-left)/2*weights)
        N=np.concatenate(ys);w=np.concatenate(ws)*np.sqrt(Q*b/np.pi)*np.exp(-Q*b*N*N)
        h=1-k*N;hs=-ks*N
        assert min(h)>.8
        mm=self.mf(*(p[:,None]+n[:,None]*N),U)
        P=1-Q*c*N**3/6;Z=a*N;G=a*(1-2*Q*b*N*N)
        delta=mm-b*N;rem=delta-c*N*N/2
        upper=G*(-k*N+Q*c*N**3/6)/h+Q*cs*N**3/(6*h)+Xi*(1/h-1)*P+hs*P/(2*h*h)-Q*a*N*delta
        lower=Q*rem-Q*Q*c*N**3*delta/6+aas*N/h+Z*G/h-Xi*Z*(1/h+1)-hs*Z/(2*h*h)
        sx=np.array([[0,1],[1,0]],complex);sy=np.array([[0,-1j],[1j,0]]);sz=np.diag([1.,-1.])
        vec=np.array([P,Z])
        dn=np.array([-Q*c*N*N/2-Q*b*N*P,a-Q*b*N*Z])
        ds=np.array([-Q*cs*N**3/6+P*(G-Xi),aas*N+Z*(G-Xi)])
        direct=1j*sx@dn-1j*sz@ds/h+1j*sz@vec*hs/(2*h*h)+sy@vec*(Q*mm)-1j*Xi*vec
        assert np.max(abs(direct-1j*np.array([upper,lower])))<1e-7
        chi=np.ones_like(N);chi_y=np.zeros_like(N);chi_s=np.zeros_like(N)
        transition=abs(N)>tube/2;angle=np.pi*(abs(N[transition])/tube-.5)
        chi[transition]=np.cos(angle)**2
        chi_y[transition]=-np.pi/tube*np.sin(2*angle)*np.sign(N[transition])
        chi_s[transition]=np.pi*abs(N[transition])*tube_s/tube**2*np.sin(2*angle)
        upper,lower=chi*upper+chi_y*Z-chi_s*P/h,chi*lower+chi_y*P+chi_s*Z/h
        norm=np.sqrt(np.sum(w*(upper*upper+lower*lower)))
        radial_cos=float(t@np.array([np.sin(theta),np.cos(theta)]))
        assert radial_cos>0
        return float(norm),radial_cos


def verify_algebra():
    x,e,u=sp.symbols('x e u', positive=True)
    flux=x*(1-e**2)/(x*x-e*e)
    esq=x*(1-u*x)/(x-u)
    assert sp.factor((flux-u).subs(e**2,esq)) == 0
    rsq=(x*x-1)*(1-e*e)/4
    assert sp.factor(rsq.subs(e**2,esq)-u*(x*x-1)**2/(4*(x-u))) == 0
    # Prolate metric: grad(U)^2/r^2 is the physical magnetic magnitude squared.
    b2=sp.factor(4*((x*x-1)*sp.diff(flux,x)**2+(1-e*e)*sp.diff(flux,e)**2)
                 /((x*x-e*e)*rsq))
    target=16*(x*x+e*e*(3*x*x-1)+e**4)/(x*x-e*e)**4
    assert sp.factor(b2-target)==0
    # Numerator >= x^2, denominator <= x^8: b >= 4/x^3 >= 4u^3.
    mp.mp.dps=70
    points=[('0.00001','-0.5001'),('.00003','.50001'),
            ('10000','-10000'),('100000','0.2'),('.5','0')]
    errors=[]
    for rr,zz in points:
        # Compare the same binary inputs; near a center, decimal-to-binary
        # input rounding otherwise dominates the algorithm error.
        r0,z0=float(rr),float(zz)
        r,z=mp.mpf(r0),mp.mpf(z0)
        exact=((z+mp.mpf('.5'))/mp.sqrt(r*r+(z+mp.mpf('.5'))**2)
              -(z-mp.mpf('.5'))/mp.sqrt(r*r+(z-mp.mpf('.5'))**2))/2
        actual=mp.mpf(str(stable_flux(np.longdouble(r0),np.longdouble(z0))))
        errors.append(float(abs((actual-exact)/exact)))
    # Windows longdouble may have only double precision; tolerance covers that ABI.
    assert max(errors)<2e-12, errors
    return dict(prolate_identities=3,stable_flux_relative_errors=errors,
                longdouble_bits=int(np.finfo(np.longdouble).nmant+1))


def verify():
    algebra=verify_algebra()
    model=StableDipole()
    kappa=1+np.arcsinh(np.sqrt(3))/(2*np.sqrt(3))
    bulk=[];overlaps=[];geometry=[]
    for exponent in [20,24,28]:
        q=2**exponent;delta=np.log(q)**-2;nu=delta*delta
        a=.4*q**(-.25)
        for harmonic in [2,3,5]:
            target=kappa*delta/(np.pi*(2*harmonic-1))
            u=(round(q*target-.5)+.5)/q
            assert q*(1-u)%1==.5
            assert abs(u-target)<=.5/q
            integrals=[]
            for order,step,normal_order in [(96,1e-4,256),(160,5e-5,512)]:
                nodes,weights=roots_legendre(order)
                lo,hi=np.log(a),np.log(1/(2*u))
                radii=np.exp((lo+hi)/2+(hi-lo)/2*nodes)
                values=[model.section(r,u,q,nu,step,normal_order) for r in radii]
                integrals.append(float((hi-lo)*np.sum(weights*radii*
                    np.array([v[0]**2/v[1] for v in values]))))
            relative=abs(integrals[0]-integrals[1])/max(integrals[1],1e-300)
            assert relative<.015,(q,harmonic,integrals,relative)
            bulk.append(dict(q=q,harmonic=harmonic,flux_label=u,delta=delta,nu=nu,
                residual_squared=integrals[1],relative_refinement=relative,
                scaled_residual=float(np.sqrt(integrals[1]*q*a))))
            for factor in [1.,2.]:
                coarse=model.overlap(factor*a,u,q,nu,1025)
                fine=model.overlap(factor*a,u,q,nu,2049)
                error=max(abs(coarse[k]-fine[k]) for k in coarse)
                assert error<2e-5,(q,harmonic,factor,error)
                overlaps.append(dict(q=q,harmonic=harmonic,radius=factor*a,
                                     refinement=error,**fine))
            for rho in [a,1.,1/(2*u)]:
                _,point=model.on_sphere(rho,u)
                geom=model.geometry(point,u)
                r,z=point;rm=np.hypot(r,z+.5);rp=np.hypot(r,z-.5)
                xi=rm+rp;eta=rm-rp
                b2=16*(xi*xi+eta*eta*(3*xi*xi-1)+eta**4)/(xi*xi-eta*eta)**4
                rel=abs(geom['b'][0]**2/b2-1)
                assert rel<2e-8,(q,harmonic,rho,rel)
                assert geom['b'][0]>=4*u**3*(1-1e-9)
                geometry.append(dict(q=q,harmonic=harmonic,rho=rho,
                    b_identity_relative_error=float(rel)))
    # Exact asymptotic bookkeeping: x=log(q), delta=x^-2.
    # e^(-alpha*x)*x^beta tends to zero for every alpha>0, finite beta.
    scales={
        'free_channel_count': {'q_power':1,'log_power':-6,'limit':'infinity'},
        'response_channel_count': {'q_power':1,'log_power':-10,'limit':'infinity'},
        'soft_pair_bound': {'q_power':0,'log_power':-1,'limit':'zero'},
        'open_photon_weighted_error': {'q_power':-1,'log_power':6,'limit':'zero'},
        'rotation_collective_error': {'q_power':-.5,'log_power':5,'limit':'zero'},
        'spherical_reference_relative_error': {'q_power':-.25,'log_power':3,'limit':'zero'},
    }
    for row in scales.values():
        growth=row['q_power'] if row['q_power'] else row['log_power']
        assert (growth>0)==(row['limit']=='infinity')
    inputs=['checks/joint_mouth_matching.py','checks/dipole_source_matching.py',
            'notes/joint-mouth-matching.md','notes/dipole-source-matching.md',
            'notes/logarithmic-long-return.md']
    return dict(status='passed_logarithmic_long_return_checks',algebra=algebra,
        bulk_cases=bulk,overlap_cases=overlaps,geometry_cases=geometry,scales=scales,
        rejected_scaling='delta=1/log(q): delta*log(2+q*delta^2) tends to 1',
        accepted_scaling='delta=1/log(q)^2; g^2=q^-2; nu=delta^2',
        theorem_scope='Conditional loaded flat dipole; finite retained long-return harmonics',
        numerics_are_supplementary=True,explicit_trim_exponent_computed=False,
        sampled_bulk_errors_are_small=False,
        numerical_limitation='Sampled far-field Gaussian is not confined inside the narrow tube; large bulk residuals are retained, not evidence for finite-q return accuracy.',
        finite_flux_threshold_certified=False,loaded_operator_existence_proved=False,
        full_vacuum_response_derived=False,physical_endpoint_achieved=False,
        novelty_established=False,lean_invocations=0,
        input_sha256={p:text_hash(ROOT/p) for p in inputs},
        versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__,mpmath=mp.__version__))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True)
    args=parser.parse_args();result=verify()
    result['script_sha256']=text_hash(Path(__file__))
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],bulk_cases=len(result['bulk_cases']),
        overlap_cases=len(result['overlap_cases']),geometry_cases=len(result['geometry_cases']),
        max_bulk_refinement=max(r['relative_refinement'] for r in result['bulk_cases']))))
