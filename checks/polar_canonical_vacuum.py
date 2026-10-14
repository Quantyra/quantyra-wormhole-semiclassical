"""Constrained quadratic polar action and free intrinsic vacuum; no formation claim."""
from pathlib import Path
import argparse
import hashlib
import json
import math

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad
from scipy.special import eval_gegenbauer, gammaln
import sympy as s

ROOT = Path(__file__).resolve().parents[1]


def zero(expr):
    assert s.simplify(s.expand(expr)) == 0, str(s.factor(expr))


def adm_tensor_checks():
    y, c = s.symbols('y c', real=True)
    b, N, k2, e = s.symbols('b N k2 e', positive=True)
    dc, dD, dA, beta, bp, alpha, ap = s.symbols('dc dD dA beta bp alpha ap')
    chi, D = s.Function('chi')(y), s.Function('D')(y)
    A = s.symbols('A')
    g = s.diag(1, b*b/(1-c*c), b*b*(1-c*c))
    gi = g.inv()
    def d(f, i):
        return s.diff(f, (y, c)[i]) if i < 2 else s.S.Zero
    gam = {}
    for i in range(3):
        for j in range(3):
            for k in range(3):
                gam[i,j,k] = s.factor(sum(gi[i,z]*(d(g[z,k],j)+d(g[z,j],k)-d(g[j,k],z))/2 for z in range(3)))
    ric = s.zeros(3)
    for i in range(3):
        for j in range(3):
            ric[i,j] = s.factor(sum(d(gam[k,i,j],k)-d(gam[k,i,k],j)+sum(gam[k,i,j]*gam[z,k,z]-gam[z,i,k]*gam[k,j,z] for z in range(3)) for k in range(3)))
    rows = []
    for ell in [2,3,4]:
        L = s.Integer(ell*(ell+1))
        Y = s.legendre(ell,c)
        norm = s.integrate(Y*Y,(c,-1,1))
        def avg(f):
            return s.factor(s.integrate(s.cancel(f),(c,-1,1))/norm)
        h = s.diag(D*Y,2*b*b*chi*Y/(1-c*c),2*b*b*chi*Y*(1-c*c))
        hd = s.diag(dD*Y,2*b*b*dc*Y/(1-c*c),2*b*b*dc*Y*(1-c*c))
        shift = s.Matrix([beta*Y,alpha*s.diff(Y,c),0])
        def ds(i,j):
            partial = shift[i].subs({beta:bp,alpha:ap}) if j == 0 else d(shift[i],j)
            return partial-sum(gam[k,i,j]*shift[k] for k in range(3))
        K = s.Matrix(3,3,lambda i,j: (hd[i,j]-ds(i,j)-ds(j,i))/(2*N))
        kinetic = avg(b*b*N/(2*k2)*(s.trace(gi*K*gi*K)-s.trace(gi*K)**2))
        expected = b*b/(2*k2*N)*(-2*dc*dc-2*L*alpha*dc/b**2-L*alpha**2/b**4-(dD-2*bp)*(2*dc+L*alpha/b**2)+L*(beta+ap)**2/(2*b**2))
        zero(kinetic-expected)
        V = s.Matrix([0,0,-(1-c*c)*s.diff(Y,c)])
        F = s.zeros(3); F[1,2]=-e*b*b; F[2,1]=e*b*b
        electric = dA*V+s.Matrix([sum((gi*shift)[j]*F[i,j] for j in range(3)) for i in range(3)])
        zero(avg(b*b/(2*N)*(electric.T*gi*electric)[0])-L*(dA-e*alpha)**2/(2*N))
        dgi = -gi*h*gi
        dg = {}
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    dg[i,j,k]=s.factor(sum((gi[i,z]*(d(h[z,k],j)+d(h[z,j],k)-d(h[j,k],z))+dgi[i,z]*(d(g[z,k],j)+d(g[z,j],k)-d(g[j,k],z)))/2 for z in range(3)))
        dric=s.zeros(3)
        for i in range(3):
            for j in range(3):
                dric[i,j]=s.factor(sum(d(dg[k,i,j],k)-d(dg[k,i,k],j)+sum(dg[k,i,j]*gam[z,k,z]+gam[k,i,j]*dg[z,k,z]-dg[z,i,k]*gam[k,j,z]-gam[z,i,k]*dg[k,j,z] for z in range(3)) for k in range(3)))
        dscalar=s.factor(s.trace(dgi*ric+gi*dric))
        zero(dscalar/Y-(-4*s.diff(chi,y,2)+L*D/b**2+(2*L-4)*chi/b**2))
        dF=s.zeros(3);dF[1,2]=A*L*Y;dF[2,1]=-A*L*Y
        dB2=s.S.Zero
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for z in range(3):
                        dB2+=(dF[i,j]*F[k,z]*gi[i,k]*gi[j,z]+F[i,j]*F[k,z]*dgi[i,k]*gi[j,z])
        zero(dB2/Y-(-2*e*L*A/b**2-4*e*e*chi))
        rows.append({'ell':ell,'gravity_kinetic':'exact','Maxwell_kinetic':'exact','linear_spatial_Ricci':'exact','linear_magnetic_energy':'exact'})
        print('ADM tensor ell',ell,'PASS',flush=True)
    return rows


def reduction_checks():
    y=s.symbols('y',real=True)
    b,k2,C,L,e=s.symbols('b k2 C L e',positive=True)
    N,v,w,D,beta,alpha=[s.Function(t)(y) for t in ['N','v','w','D','beta','alpha']]
    n=s.diff(N,y)/N
    kin=b*b/(2*k2*N)*(-2*v*v-2*L*alpha*v/b**2-L*alpha**2/b**4-(D-2*s.diff(beta,y))*(2*v+L*alpha/b**2)+L*(beta+s.diff(alpha,y))**2/(2*b**2))+L/(2*N)*(w-e*alpha)**2
    eb=s.diff(kin,beta)-s.diff(s.diff(kin,s.diff(beta,y)),y)
    ea=s.diff(kin,alpha)-s.diff(s.diff(kin,s.diff(alpha,y)),y)
    beta0=4*b*b/L*(s.diff(v,y)-n*v)
    D0=4*b*b/L*s.diff(v,y,2)-(2+4*C/L)*v-2*k2*e*w
    def gauge(expr):
        return s.simplify(expr.subs(alpha,0).doit().subs(beta,beta0).doit().subs(D,D0).doit().subs(s.diff(N,y,2),C*N/b**2))
    zero(gauge(eb));zero(gauge(ea))
    reduced=(b*b/k2*(1+2*C/L)*v*v+2*b*b*e*v*w+L*w*w/2)/N
    # Compare densities, retaining the exact integration-by-parts boundary term.
    boundary=4*b**4/(k2*L)*(v*s.diff(v,y)/N-3*s.diff(N,y)*v*v/(2*N*N))
    zero(gauge(kin)-reduced-s.diff(boundary,y).subs(s.diff(N,y,2),C*N/b**2))
    W=s.Matrix([[2*b*b/k2*(1+2*C/L),2*b*b*e],[2*b*b*e,L]])
    K=s.Matrix([[(L+2*C)/b**2,L*k2*e/b**2],[2*e,L/b**2]])
    back=lambda x:s.factor(s.expand(x).subs(e**2,(C+1)/(k2*b**2)))
    assert (W*K-K.T*W).applyfunc(back)==s.zeros(2)
    zero(back(W.det())-2*b*b*(L-2)/k2)
    S=s.diag(s.sqrt(2*b*b/(L*k2)),1)
    Wz=(S.inv()*W*S.inv()).subs(e,s.sqrt((C+1)/(k2*b*b))).applyfunc(s.simplify)
    CM=s.Matrix([[L+2*C,s.sqrt(2*L*(C+1))],[s.sqrt(2*L*(C+1)),L]])
    assert (Wz-CM).applyfunc(s.simplify)==s.zeros(2)
    zero(s.simplify(CM.inv()[0,0])-1/(L-2))
    source_v=s.Matrix([s.Rational(9,4)/b**2,-3*e/4])
    assert (W.subs(L,6)*source_v-s.Matrix([3/k2,0])).applyfunc(back)==s.zeros(2,1)
    return {'radial_and_angular_shift_constraints':'exact','kinetic_boundary_term':str(boundary),'Wz':'C M_L','det_W':'2 b^2 (L-2)/kappa^2','inverse_Wz_11':'1/(L-2)','quadrupole_Wv_over_G':['3/kappa^2','0']}


def spectral_checks():
    mp.mp.dps=65
    rows=[]
    for C in map(mp.mpf,['1','4','100','1e10']):
        for ell in [2,3,20,500]:
            L=mp.mpf(ell*(ell+1));R=mp.sqrt(C*C+2*L*(C+1))
            mum=L*(L-2)/(C*(C+L+R));mup=(C+L+R)/C
            A=mp.matrix([[C,mp.sqrt(2*L*(C+1))],[mp.sqrt(2*L*(C+1)),-C]])/R
            pp=(mp.eye(2)+A)/2;pm=(mp.eye(2)-A)/2
            CM=mp.matrix([[L+2*C,mp.sqrt(2*L*(C+1))],[mp.sqrt(2*L*(C+1)),L]])
            cov=pm/(C*mum)+pp/(C*mup)
            err=max(abs(v) for v in CM*cov-mp.eye(2))
            assert err<mp.mpf('1e-40')
            if C==1:
                assert abs(mum-ell*(ell-1))<mp.mpf('1e-50')
                assert abs(mup-(ell+1)*(ell+2))<mp.mpf('1e-50')
            rows.append({'C':str(C),'ell':ell,'CCR_matrix_error':str(err)})
    # Check the normalized radial mode independently by quadrature.
    radial=[]
    for h in [1.,1.7,3.,8.]:
        for n in range(6):
            k=2*n
            lognorm=math.lgamma(k+2*h)-math.lgamma(k+1)-math.log(k+h)-2*math.lgamma(h)+math.log(math.pi)+(1-2*h)*math.log(2)
            def mode(sig):
                return math.exp(-lognorm/2)*math.cos(sig)**h*eval_gegenbauer(k,h,math.sin(sig))
            norm,error=quad(lambda x:mode(x)**2,-math.pi/2,math.pi/2,epsabs=1e-11,epsrel=1e-11)
            assert abs(norm-1)<2e-10
            actual=mode(0)**2/(2*(2*n+h))
            expected=math.exp(gammaln(n+.5)+gammaln(n+h)-gammaln(n+1)-gammaln(n+h+.5))/(2*math.pi)
            assert abs(actual/expected-1)<2e-12
            radial.append({'h':h,'n':n,'norm_error':abs(norm-1),'central_weight_relative_error':abs(actual/expected-1)})
    return {'CCR':rows,'radial_normalization':radial}


def infinite_sum_checks():
    t=s.symbols('t',nonnegative=True)
    zero((t/2+s.Rational(17,8))**2-(1+4*t)-(4*t-15)**2/64)
    constant=s.Rational(6)*s.Rational(9,10)*s.Rational(128,127)*40
    assert constant==s.Rational(27648,127) and constant<225
    assert (1+mp.sqrt(mp.pi))/mp.pi<mp.mpf('.9')
    assert mp.exp(mp.mpf(1)/128)<mp.mpf(128)/127
    # Finite partial sums are diagnostics; the ordinary proof controls all tails.
    rows=[]
    for C in [1.,4.,100.,10000.]:
        alpha=1/(32*C)
        ell=np.arange(0,int(math.ceil(12/math.sqrt(alpha)))+1,dtype=float)
        heat=float(np.sum((2*ell+1)*np.exp(-alpha*ell*(ell+1))))
        assert heat<40*C
        rows.append({'C':C,'ell_max':int(ell[-1]),'heat_partial_sum':heat,'analytic_heat_bound':40*C})
    direct=[]
    for C in [1.,4.,100.]:
        total=0.
        lmax=int(math.ceil(80*math.sqrt(C)))
        n=np.arange(50,dtype=float)
        for ell in range(2,lmax+1):
            L=float(ell*(ell+1));R=math.sqrt(C*C+2*L*(C+1))
            mum=L*(L-2)/(C*(C+L+R));mup=(C+L+R)/C
            # Safe small branch projector avoids subtracting near-equal numbers.
            pm=L*(C+1)/(R*(R+C));pp=1-pm
            for mu,p in [(mum,pm),(mup,pp)]:
                h=(1+math.sqrt(1+4*mu))/2
                an=np.exp(gammaln(n+.5)+gammaln(n+h)-gammaln(n+1)-gammaln(n+h+.5))/(2*math.pi)
                V=float(np.sum(an*np.exp(-(2*n+h)**2/16)))
                total+=4*(2*ell+1)*L*p*V/(C*C*mu)
        # total is Var(h_sphere)/(G/a^2), retaining all multiplicities.
        assert 0<total<float(constant)
        direct.append({'C':C,'ell_max':lmax,'even_radial_terms':50,'dimensionless_variance_partial_sum':total,'proven_all_mode_upper_bound':str(constant)})
    k=mp.mpf(2000);C=mp.sqrt(1+12*mp.pi*k**6);a2=k**10*(C-1)/(6*C)
    mu=24/(C*(C+6+mp.sqrt(C*C+12*(C+1))))
    p=6*(C+1)/(mp.sqrt(C*C+12*(C+1))*(mp.sqrt(C*C+12*(C+1))+C))
    return {'constant':str(constant),'RMS_bound':'15 sqrt(G)/a','candidate_k2000':{'C':str(C),'a2_over_G':str(a2),'RMS_bound':str(15/mp.sqrt(a2)),'light_metric_projector':str(p),'light_kinetic_coefficient':str(C*mu),'projector_over_kinetic':str(p/(C*mu))},'heat_diagnostics':rows,'vacuum_diagnostics':direct,'tail_proof':'See note: positive radial gamma bound, mass inequality, unimodal sphere heat sum. No truncation is used to prove the all-mode result.'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    out={'story':'S072','scope':'Free intrinsic l>=2 polar Einstein-Maxwell vacuum on fixed global product, Friedrichs boundaries; central Gaussian time smear s=1/4.',
         'ADM':adm_tensor_checks(),'reduction':reduction_checks(),'spectral':spectral_checks(),'infinite_sums':infinite_sum_checks()}
    out['versions']={'sympy':s.__version__,'mpmath':mp.__version__,'numpy':np.__version__,'scipy':scipy.__version__}
    out['sources']={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in ['checks/polar_canonical_vacuum.py','notes/polar-canonical-vacuum.md','checks/quadrupole_product_response.py','notes/quadrupole-product-response.md']}
    out['status']='PASS';out['goal_complete']=False
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_bytes((json.dumps(out,indent=2)+'\n').encode())
    print('PASS: constrained kinetic action; canonical covariance; complete quadratic polar sphere bound',flush=True)


if __name__=='__main__':
    main()
