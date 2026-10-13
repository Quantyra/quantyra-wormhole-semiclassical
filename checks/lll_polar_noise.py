"""Ordered free-LLL polar stress noise, not a constrained interacting state."""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import mpmath as mp
import numpy as np
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'checks'))
from lll_projection import form_factor, matrix_element


def tensor_response(ell):
    # Independent direct 4D calculation: a=1, b^2=C, kappa^2=2.
    r, c, w, C, p, m = s.symbols('r c w C p m', real=True)
    B = 1+r*r
    Y = s.legendre(ell, c)
    L = s.Integer(ell*(ell+1))
    def clean(z): return s.factor(s.cancel(z))
    def d(z, k, pert=False):
        if k == 0: return -s.I*w*z if pert else s.S.Zero
        if k == 1:
            return s.diff(z, r) + (s.I*w/B*(p*s.diff(z,p)-m*s.diff(z,m)) if pert else 0)
        if k == 2: return s.diff(z,c)
        return s.S.Zero
    g=s.diag(-B,1/B,C/(1-c*c),C*(1-c*c)); gi=g.inv()
    t=s.zeros(4)
    t[0,0]=p+m; t[0,1]=t[1,0]=(-p+m)/B; t[1,1]=(p+m)/B**2
    h=4*t*Y/L
    assert clean(s.trace(gi*h)) == 0
    dgi=-gi*h*gi
    gam={}; dg={}
    for i in range(4):
      for j in range(4):
       for k in range(4):
        gam[i,j,k]=clean(sum(gi[i,z]*(d(g[z,k],j)+d(g[z,j],k)-d(g[j,k],z))/2 for z in range(4)))
        dg[i,j,k]=clean(sum((gi[i,z]*(d(h[z,k],j,True)+d(h[z,j],k,True)-d(h[j,k],z,True))
                        +dgi[i,z]*(d(g[z,k],j)+d(g[z,j],k)-d(g[j,k],z)))/2 for z in range(4)))
    ric=s.zeros(4); dric=s.zeros(4)
    for i in range(4):
      for j in range(i,4):
        z=dz=0
        for k in range(4):
            z+=d(gam[k,i,j],k)-d(gam[k,i,k],j)
            dz+=d(dg[k,i,j],k,True)-d(dg[k,i,k],j,True)
            for v in range(4):
                z+=gam[k,i,j]*gam[v,k,v]-gam[v,i,k]*gam[k,j,v]
                dz+=dg[k,i,j]*gam[v,k,v]+gam[k,i,j]*dg[v,k,v]-dg[v,i,k]*gam[k,j,v]-gam[v,i,k]*dg[k,j,v]
        ric[i,j]=ric[j,i]=clean(z); dric[i,j]=dric[j,i]=clean(dz)
    scalar=s.trace(gi*ric)
    ds=clean(s.trace(dgi*ric+gi*dric))
    em=s.zeros(4); em[0,1]=1;em[1,0]=-1
    raised=gi*em*gi; dr=dgi*em*gi+gi*em*dgi
    sq=sum(em[i,j]*raised[i,j] for i in range(4) for j in range(4))
    dsq=sum(em[i,j]*dr[i,j] for i in range(4) for j in range(4))
    # Include the metric variation of the electromagnetic stress.
    dstress=(C+1)/C*(em*dgi*em.T-h*sq/4-g*dsq/4)
    residual=dric-h*scalar/2-g*ds/2+(1/C-1)*h/2-dstress-2*t*Y/C
    for i in range(4):
      for j in range(i,4): assert clean(residual[i,j])==0, (ell,i,j)
    for j in range(4): assert clean(sum(d(dr[i,j],i,True) for i in range(4)))==0
    for b in range(2):
        div=sum(gi[a,a]*(d(t[a,b],a,True)-sum(gam[v,a,a]*t[v,b]+gam[v,a,b]*t[a,v] for v in range(4))) for a in range(4))
        assert clean(div)==0
    print('direct tensor response passed, l='+str(ell), flush=True)
    return {'ell':ell,'Einstein_components':10,'Maxwell_components':4,
            'stress_divergences':2,'frequency':'arbitrary including zero',
            'radius_ratio':'symbolic C','chi_and_Acal':0}


def oscillator_checks():
    # Enumerate the actual filled-sea particle-hole transitions, retaining
    # noncommuting flavour matrices. No independent bosonic rotors are added.
    rng=np.random.default_rng(72813); rows=[]
    for q in [2,3,5]:
        a=rng.normal(size=(q,q))+1j*rng.normal(size=(q,q)); a=(a+a.conj().T)/2
        b=rng.normal(size=(q,q))+1j*rng.normal(size=(q,q)); b=(b+b.conj().T)/2
        tr=np.trace(a@b); err=0.
        for n in range(1,13):
            jj=tt=tj=0j
            for h in range(n):
                r=-n+.5+h; weight=r+n/2
                for i in range(q):
                    for j in range(q):
                        z=a[i,j]*b[j,i]
                        jj+=z; tt+=weight*weight*z; tj+=weight*z
            err=max(err,abs(jj-n*tr),abs(tt-n*(n*n-1)*tr/12),abs(tj))
        assert err<2e-10
        rows.append({'q':q,'frequencies':12,'max_absolute_error':float(err),
                     'matrix_commutator_norm':float(np.linalg.norm(a@b-b@a))})
    z=s.symbols('z'); derivative=lambda v:s.factor(z*s.diff(v,z))
    first=derivative(1/(1-z)); cubic=derivative(derivative(first))-first
    assert s.factor(first-z/(1-z)**2)==0
    assert s.factor(cubic-6*z*z/(1-z)**4)==0
    return {'particle_hole_cases':rows,'two_ordered_generating_functions':'exact'}


def angular_checks():
    cases=0
    for q in range(3,15):
      for ell in range(2,q+1):
        entries=[matrix_element(q,i,ell) for i in range(q)]
        assert sum(entries)==0
        assert s.factor(sum(v*v for v in entries)-q*form_factor(q,ell)/(2*ell+1))==0
        cases+=1
      aq=sum(s.Rational(2*l+1,l*l*(l+1)**2)*form_factor(q,l) for l in range(2,q))
      remainder=(s.harmonic(q-1)+s.harmonic(q)-s.Rational(5,2))/q+s.Rational(1,q*q)
      assert 0<=s.Rational(1,4)-aq<=remainder
    l=s.symbols('l',positive=True)
    assert s.factor((2*l+1)/(l*l*(l+1)**2)-(1/l**2-1/(l+1)**2))==0
    # Local zero-mode identity, using D1 eta=-i D2 eta and its adjoint.
    f,g,df,dg=s.symbols('f g df dg')
    assert s.expand(f*(-s.I*dg)-(s.I*df)*g+s.I*(df*g+f*dg))==0
    return {'independent_LLL_matrix_variances':cases,'angular_sum_bounds':12,
            'axial_mixed_stress_identity':'exact first-order zero-mode identity'}


def smearing_and_candidate():
    mp.mp.dps=80; rows=[]
    for width in [mp.mpf('.125'),mp.mpf('.25'),mp.mpf('.5')]:
        circumference=mp.pi; decay=(2*mp.pi*width/circumference)**2
        N=128
        total=mp.fsum(n*(n*n-1)*mp.exp(-decay*n*n) for n in range(2,N+1))
        # For N beyond the maximum of x^3 exp(-decay*x^2), omitted
        # terms are bounded by its integral from N to infinity.
        assert N*N>3/(2*decay)
        tail=mp.exp(-decay*N*N)*(N*N/decay+1/decay**2)/2
        independent=mp.nsum(lambda n:n*(n*n-1)*mp.exp(-decay*n*n),[2,mp.inf])
        assert abs(independent-total)<mp.mpf('1e-70')
        rows.append({'gaussian_time_width':str(width),'circumference':'pi',
                     'stress_spectral_sum':str(total),'positive_tail_upper':str(tail),
                     'nsum_difference':str(abs(independent-total))})
    S=mp.mpf(rows[1]['stress_spectral_sum'])+mp.mpf(rows[1]['positive_tail_upper'])
    candidates=[]
    for k0 in [2000,4000,10000]:
        k=mp.mpf(k0); q=k**4; C=mp.sqrt(1+12*mp.pi*k**6)
        a2=k**10*(C-1)/(6*C) # Planck units, L4=k^5.
        variance=(8*mp.pi)**2*q*S/(24*a2*a2*mp.pi**4)
        rms=mp.sqrt(variance)
        asym=mp.sqrt(96*S)/mp.pi*C/(C-1)*k**-8
        assert abs(rms/asym-1)<mp.mpf('1e-70')
        threshold=C+1+mp.sqrt((C+1)**2+C*C)
        lmax=int(mp.floor((mp.sqrt(1+4*threshold)-1)/2))
        eps=mp.mpf(lmax*(lmax+1))/q
        quad=1-((q-1)*(q-2))/((q+1)*(q+2))
        assert eps<1 and quad<6/q
        candidates.append({'k':k0,'a_squared_Planck_units':str(a2),
                           'smeared_polar_fractional_gtautau_RMS_upper':str(rms),
                           'quadrupole_relative_variance_form_factor_loss':str(quad),
                           'soft_branch_l_max':lmax,'soft_branch_form_factor_loss_upper':str(eps),
                           'physical_capture_probability_bound':False})
    return {'gaussian_sums':rows,'candidate':candidates}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    out={'status':'PASS','scope':'free NS Dirac LLL vacuum on fixed round product; polar constrained metric only',
         'tensor_response':[tensor_response(l) for l in [2,3]],
         'ordered_kernels':oscillator_checks(),'angular':angular_checks(),
         'smearing':smearing_and_candidate(),
         'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'dependencies':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in
                         ['notes/lll-polar-noise.md','checks/lll_projection.py','notes/quadrupole-product-response.md',
                          'notes/angular-current-renormalization.md','notes/lll-state-difference-source.md']},
         'versions':{'sympy':s.__version__,'numpy':np.__version__,'mpmath':mp.__version__},
         'limitations':['Free comparison is not the interacting angular-Gauss state.',
                        'Zero direct polar master noise does not eliminate constrained metric or axial current noise.',
                        'Pointwise stress variance diverges; the stated metric estimate is time smeared.',
                        'Induced Hall/anomaly response, nonlinear dipole sources, massive levels and boundary matching remain.',
                        'No physical capture, preparation, stability, transport or creation bound.']}
    args.output.write_bytes((json.dumps(out,indent=2)+'\n').encode())
    print(json.dumps({'status':'PASS','candidate':out['smearing']['candidate'][0]},indent=2),flush=True)
