"""Quadratic bosonic dipole source and ordered vacuum kernel on a fixed product."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import numpy as np
import sympy as s
from scipy.special import eval_gegenbauer, gammaln, roots_legendre

ROOT=Path(__file__).resolve().parents[1]


def tensor_source_check():
    t,r,c,C,P=s.symbols('t r c C P',real=True)
    W,H,V=[s.Function(n)(t,r) for n in ['W','H','V']]
    B=1+r*r; S=1-c*c
    coords=[t,r,c,None]
    def d(x,k):return s.diff(x,coords[k]) if k<3 else s.S.Zero
    def clean(x):return s.factor(s.cancel(x))
    def back(x):
        x=s.cancel(x);n,den=s.fraction(x)
        return s.factor(s.rem(s.Poly(n,P),s.Poly(P*P-C*(C+1)/2,P)).as_expr()/den)
    g=[s.diag(-B,1/B,C/S,C*S),s.zeros(4),s.zeros(4)]
    g[1][0,3]=g[1][3,0]=-C*S*W;g[2][0,0]=C*S*W*W
    gi=[g[0].inv(),s.zeros(4),s.zeros(4)]
    gi[1][0,3]=gi[1][3,0]=-W/B;gi[2][3,3]=-W*W/B
    for n in range(3):
        target=sum((gi[k]*g[n-k] for k in range(n+1)),s.zeros(4))
        assert target.applyfunc(clean)==(s.eye(4) if n==0 else s.zeros(4))
    gam=[]
    for n in range(3):
        gmn={}
        for i in range(4):
          for j in range(4):
           for k in range(4):
            gmn[i,j,k]=clean(sum(gi[a][i,v]*(d(g[n-a][v,k],j)+d(g[n-a][v,j],k)-d(g[n-a][j,k],v))/2 for a in range(n+1) for v in range(4)))
        gam.append(gmn)
        print('connection order',n,flush=True)
    ric=[]
    for n in range(3):
        rr=s.zeros(4)
        for i in range(4):
          for j in range(i,4):
            z=sum(d(gam[n][k,i,j],k)-d(gam[n][k,i,k],j) for k in range(4))
            for a in range(n+1):
             for k in range(4):
              for v in range(4):z+=gam[a][k,i,j]*gam[n-a][v,k,v]-gam[a][v,i,k]*gam[n-a][k,j,v]
            rr[i,j]=rr[j,i]=clean(z)
        ric.append(rr)
        print('Ricci order',n,flush=True)
    scalar=[clean(sum(s.trace(gi[k]*ric[n-k]) for k in range(n+1))) for n in range(3)]
    em=[s.zeros(4) for _ in range(3)]
    em[0][2,3]=-P;em[0][3,2]=P
    G=d(V,0)-d(P*W+H,1)
    em[1][0,1]=c*G;em[1][1,0]=-c*G
    em[1][0,2]=-(P*W+H);em[1][2,0]=P*W+H
    em[1][1,2]=-V;em[1][2,1]=V
    raised=[]; square=[]; stress=[]; E=[]; J=[]
    for n in range(3):
        ff=sum((gi[a]*em[b]*gi[n-a-b] for a in range(n+1) for b in range(n-a+1)),s.zeros(4)).applyfunc(clean)
        raised.append(ff)
        square.append(clean(sum(em[a][i,j]*raised[n-a][i,j] for a in range(n+1) for i in range(4) for j in range(4))))
        tt=sum((em[a]*gi[b]*em[n-a-b].T for a in range(n+1) for b in range(n-a+1)),s.zeros(4))
        tt-=sum((g[a]*square[n-a]/4 for a in range(n+1)),s.zeros(4))
        stress.append(tt.applyfunc(clean))
        ee=ric[n]-sum((g[a]*scalar[n-a]/2 for a in range(n+1)),s.zeros(4))+(1/C-1)*g[n]/2-2*stress[n]
        E.append(ee.applyfunc(back))
        J.append([clean(sum(d(raised[n][i,j],i) for i in range(4))) for j in range(4)])
        print('Einstein-Maxwell order',n,flush=True)
    assert E[0]==s.zeros(4)
    horizontal_tt=E[2][0,0]+2*W*E[1][0,3]+W*W*E[0][3,3]
    strace=clean(horizontal_tt/B-B*E[2][1,1])
    h2=-H*H/B+B*V*V
    predtrace=-2*c*c*G*G-C*S*d(W,1)**2/2
    assert back(strace-predtrace)==0,('base trace',str(clean(strace-predtrace)))
    coef=clean(-E[2][2,2]+E[2][3,3]/S**2)
    predcoef=2*h2-C*C*d(W,1)**2/2
    assert back(coef-predcoef)==0,('sphere tracefree',str(clean(coef-predcoef)))
    trace_l2=s.expand(strace).coeff(c,2)*s.Rational(2,3)
    ST=s.Rational(2,3)*coef
    rF=clean(trace_l2/2-6*ST/(4*C))
    predF=-s.Rational(2,3)*G*G-2*h2/C+2*C*d(W,1)**2/3
    assert back(rF-predF)==0
    assert clean(J[2][3]-W*J[1][0]-c*d(W,1)*G)==0
    assert clean(J[1][0]-c*(d(G,1)+2*H/(C*B)))==0
    assert clean(J[1][1]+c*(d(G,0)+2*B*V/C))==0
    # The polar Maxwell equation has the same scalar master as the
    # dual-electric convention: Acal=original axial potential/g.
    aa,ch,hh=[s.Function(n)(t,r) for n in ['aa','ch','hh']]
    Y=s.legendre(2,c)
    polarF=s.zeros(4)
    potential=[0,0,0,-S*s.diff(Y,c)*aa]
    for i in range(4):
      for j in range(4):polarF[i,j]=d(potential[j],i)-d(potential[i],j)
    metric=s.zeros(4)
    metric[0,0]=-B*hh*Y/2;metric[1,1]=hh*Y/(2*B)
    metric[2,2]=2*C*ch*Y/S;metric[3,3]=2*C*ch*S*Y
    inv=-gi[0]*metric*gi[0]
    vol=(hh+4*ch)*Y/2
    variation=inv*em[0]*gi[0]+gi[0]*polarF*gi[0]+gi[0]*em[0]*inv+vol*raised[0]
    jphi=clean(sum(d(variation[i,3],i) for i in range(4)))
    box=lambda z:-s.diff(z,t,2)/B+s.diff(B*s.diff(z,r),r)
    expected=-s.diff(Y,c)/C*(box(aa)-6*aa/C-2*P*ch/C+P*hh/(2*C))
    assert clean(jphi-expected)==0,('linear polar Maxwell',str(clean(jphi-expected)))
    # Setting hh=-S_T adds the established E0*S_T/2 source. The
    # quadratic angular residual supplies a further C*W_r*G/3.
    out_polar={'direct_quadratic_Maxwell_source':'b^2 W_r G/(3g)',
               'linear_polar_Maxwell_identity':'exact, including base trace'}
    z=s.Function('z')(t,r)
    kg_rhs=B*(s.diff(B*s.diff(z,r),r)-(2+4/C)*z)
    def on_shell(expr):
        expr=expr.subs({H:-B*s.diff(z,r),V:-s.diff(z,t)/B},simultaneous=True).doit()
        rules={}
        for derivative in expr.atoms(s.Derivative):
            if derivative.expr==W and derivative.variables.count(r):
                nr=derivative.variables.count(r);nt=derivative.variables.count(t)
                rules[derivative]=4*P/C**2*s.diff(z,t,nt,r,nr-1)
        expr=expr.xreplace(rules)
        rules={}
        for derivative in expr.atoms(s.Derivative):
            if derivative.expr==z and derivative.variables.count(t)>=2:
                nr=derivative.variables.count(r);nt=derivative.variables.count(t)
                rules[derivative]=s.diff(kg_rhs,t,nt-2,r,nr)
        return back(expr.xreplace(rules).doit())
    for i in range(4):
      for j in range(i,4):assert on_shell(E[1][i,j])==0,('linear dipole Einstein',i,j)
    for i in range(4):assert on_shell(J[1][i])==0,('linear dipole Maxwell',i)
    out={'status':'PASS','quadratic_base_trace':str(predtrace),'quadratic_sphere_tracefree_coefficient':str(predcoef),
         'rF':str(predF),'horizontal_angular_Maxwell_residual':str(c*d(W,1)*G),
         'units':'a=1,b^2=C,kappa^2=2,g=1,P^2=C(C+1)/2',
         'method':'direct coordinate metric/inverse, second-order Ricci, Maxwell tensor and stress, horizontal projection'}
    out['polar_Maxwell']=out_polar
    out['on_shell_linear_dipole_tensor_components']=14
    return out


def reduction_checks():
    b2,Gn,ee,PP=s.symbols('b2 Gn ee PP',positive=True)
    kg=b2**2/(6*Gn);ke=4*s.pi*b2/(3*ee**2);mu=8*s.pi/(3*ee**2)
    ym=4*s.pi*Gn*PP**2/(ee**2*b2)
    mass=s.factor(mu*(1/ke+PP**2/kg))
    assert s.factor(mass-(2+4*ym)/b2)==0
    grad,phi2=s.symbols('grad phi2')
    RF=-8*s.pi*Gn*mu*phi2/(3*ee**2*ke**2)+8*s.pi*Gn*grad/(ee**2*b2*mu)+2*b2*PP**2*mu*phi2/(3*kg**2)
    target=3*Gn*grad/b2+4*Gn*(4*ym-1)*phi2/b2**2
    assert s.factor(RF-target)==0
    E0=PP/(ee*b2)
    ST=-2*Gn*grad-8*Gn*ym*phi2/b2
    RA=E0*ST/2+b2*PP*mu*phi2/(3*ee*kg*ke)
    assert s.factor(RA-(-Gn*E0*grad+4*Gn*E0*(1-ym)*phi2/b2))==0
    boxq,cc=s.symbols('boxq cc')
    polar=s.expand(s.Rational(3,2)*target.subs(grad,boxq/2-mass*phi2))
    expected=9*Gn*boxq/(4*b2)+3*Gn*(2*ym-5)*phi2/b2**2
    assert s.factor(polar-expected)==0
    polarA=s.expand(s.Rational(3,2)*RA.subs(grad,boxq/2-mass*phi2))
    assert s.factor(polarA-(-3*Gn*E0*boxq/4+9*Gn*E0*phi2/b2))==0
    # Complete source includes the anisotropic sphere constraint. Its
    # derivative piece is a local composite field contribution.
    Csym,E=s.symbols('Csym E',positive=True)
    v=s.Matrix([9*Gn/(4*b2),-3*Gn*E/4])
    K=s.Matrix([[(6+2*Csym)/b2,6*8*s.pi*Gn*E/b2],[2*E,6/b2]])
    Kv=(K*v).subs(E**2,(Csym+1)/(8*s.pi*Gn*b2)).applyfunc(s.factor)
    assert Kv==s.Matrix([9*Gn/b2**2,0])
    # Radial gauge: exact round-projection bosonic Hamiltonian after Gauss.
    f=s.symbols('f',positive=True)
    p=s.Matrix(s.symbols('p0:3'));x=s.Matrix(s.symbols('x0:3'));flux=s.Matrix(s.symbols('K0:3'))
    before=f*(flux-PP*s.sqrt(mu)*x).dot(flux-PP*s.sqrt(mu)*x)/(2*kg)+f*mu*x.dot(x)/(2*ke)+p.dot(p)/2
    after=p.dot(p)/2+f*mass*x.dot(x)/2-f*PP*s.sqrt(mu)*flux.dot(x)/kg+f*flux.dot(flux)/(2*kg)
    assert s.expand(before-after)==0
    rho=x.cross(p)
    def pb(a,b):return s.expand(sum(s.diff(a,x[i])*s.diff(b,p[i])-s.diff(a,p[i])*s.diff(b,x[i]) for i in range(3)))
    for i in range(3):
      for j in range(3): assert s.simplify(pb(rho[i],rho[j])-sum(s.LeviCivita(i,j,k)*rho[k] for k in range(3)))==0
    return {'Proca_mass':'mu*(1/kappae+P^2/kappag)=(2+4Y)/b^2',
            'canonical_rF_and_rA':'exact, including E0*S_T/2','composite_wave_reduction':'exact',
            'local_composite_redefinition':'K v=(9G/b^4,0)',
            'Gauss_reduced_cubic_quartic_Hamiltonian':'exact','rotation_Poisson_brackets':9}


def modes(C,N,x):
    mass=2+4/C;delta=(1+math.sqrt(1+4*mass))/2
    uu=[];du=[];ddu=[]
    for n in range(N):
        norm=math.exp(-.5*(math.log(math.pi)+(1-2*delta)*math.log(2)+gammaln(n+2*delta)-gammaln(n+1)-math.log(n+delta)-2*gammaln(delta)))
        co=math.cos(x);si=math.sin(x);poly=eval_gegenbauer(n,delta,si)
        dp=2*delta*eval_gegenbauer(n-1,delta+1,si) if n else 0.
        v=norm*co**delta*poly
        dv=norm*(-delta*si*co**(delta-1)*poly+co**(delta+1)*dp)
        uu.append(v);du.append(dv);ddu.append((mass/co**2-(n+delta)**2)*v)
    return delta,np.arange(N)+delta,np.array(uu),np.array(du),np.array(ddu)


def source_profiles(C,N,x):
    delta,om,u,du,ddu=modes(C,N,x)
    pair=np.outer(u,u);second=np.outer(ddu,u)+2*np.outer(du,du)+np.outer(u,ddu)
    freq=om[:,None]+om[None,:]
    # G_N=1,a=1,b^2=C for numerical normalization checks only.
    E0=math.sqrt((C+1)/(8*math.pi*C))
    F=9/(4*C)*math.cos(x)**2*(freq**2*pair+second)+3*(C-4)/C**2*pair
    A=-3*E0/4*math.cos(x)**2*(freq**2*pair+second)+9*E0/C*pair
    return delta,om,pair,np.array([F,A])


def pair_vectors(T,profiles,om):
    N=len(om);out={}
    for a in range(3):
      for b in range(3):
       for p in range(N):
        for q in range(N):
            i=a*N+p;j=b*N+q;key=tuple(sorted((i,j)))
            value=T[a,b]*profiles[:,p,q]/(2*math.sqrt(om[p]*om[q]))
            if i==j:value=value*math.sqrt(2)
            out[key]=out.get(key,np.zeros(len(profiles),complex))+value
    return out


def spectral_kernel_checks():
    Ts=[]
    Ts.append(np.diag([1.,-1.,0.])/math.sqrt(2))
    Ts.append(np.diag([1.,1.,-2.])/math.sqrt(6))
    for a,b in [(0,1),(0,2),(1,2)]:
        v=np.zeros((3,3));v[a,b]=v[b,a]=1/math.sqrt(2);Ts.append(v)
    Ts=[T*math.sqrt(8*math.pi/15) for T in Ts]
    assert max(abs(np.trace(a@b)-(8*math.pi/15 if i==j else 0)) for i,a in enumerate(Ts) for j,b in enumerate(Ts))<1e-14
    rows=[];N=6;time=.43
    for C in [1.,4.,1.e8]:
        delta,om,pair,profiles=source_profiles(C,N,-.37)
        _,_,pair2,profiles2=source_profiles(C,N,.21)
        vp=pair_vectors(Ts[0],profiles,om);vq=pair_vectors(Ts[0],profiles2,om)
        cross=pair_vectors(Ts[2],profiles2,om)
        actual=np.zeros((2,2),complex);orth=np.zeros((2,2),complex)
        for key,amp in vp.items():
            energy=sum(om[i%N] for i in key);phase=np.exp(-1j*time*energy)
            actual+=np.outer(amp.conj(),vq[key])*phase
            orth+=np.outer(amp.conj(),cross[key])*phase
        target=np.zeros((2,2),complex)
        for p in range(N):
          for q in range(N):
            target+=(16*math.pi/15)*np.outer(profiles[:,p,q],profiles2[:,p,q])*np.exp(-1j*time*(om[p]+om[q]))/(4*om[p]*om[q])
        err=float(np.max(abs(actual-target)))/(1+float(np.max(abs(target))))
        assert err<1e-12 and np.max(abs(orth))<1e-9
        # Gauss charge rho_0 = phi_1 pi_2 - phi_2 pi_1.
        rhoT=np.zeros((3,3));rhoT[1,2]=1;rhoT[2,1]=-1
        rho_profiles=np.array([1j*om[None,:]*pair])
        rho_profiles2=np.array([1j*om[None,:]*pair2])
        vr=pair_vectors(rhoT,rho_profiles,om);vr2=pair_vectors(rhoT,rho_profiles2,om)
        rho_actual=sum(vr[key][0].conjugate()*vr2[key][0]*np.exp(-1j*time*sum(om[i%N] for i in key)) for key in vr)
        rho_target=sum((om[p]-om[q])**2*pair[p,q]*pair2[p,q]*np.exp(-1j*time*(om[p]+om[q]))/(4*om[p]*om[q]) for p in range(N) for q in range(N))
        # pair_vectors supplies the canonical 1/(2 sqrt(omega_p omega_q)).
        assert abs(rho_actual-rho_target)<1e-10,('Gauss spectrum',rho_actual,rho_target)
        rows.append({'C':C,'Delta_dipole':delta,'modes':N,'source_matrix_relative_error':err,
                     'Gauss_kernel_error':float(abs(rho_actual-rho_target))})
    return {'STF_matrix_orthogonality':25,'independent_bosonic_pair_state_checks':rows}


def smeared_Gauss_flux():
    rows=[];sigma=.25;N=32
    z,weight=roots_legendre(128)
    for C in [1.,4.,1.e8]:
        for cut in [-.4,0.,.35,math.pi/2]:
            xx=-math.pi/2+(z+1)*(cut+math.pi/2)/2
            ww=weight*(cut+math.pi/2)/2
            us=np.array([modes(C,N,float(x))[2] for x in xx]).T
            delta,om,_,_,_=modes(C,N,0.)
            overlap=(us*ww)@us.T
            freqs=om[:,None]+om[None,:]
            kernel=(om[:,None]-om[None,:])**2/(4*om[:,None]*om[None,:])*overlap**2*np.exp(-sigma*sigma*freqs**2)
            value=float(np.sum(kernel))
            low=float(np.sum(kernel[:20,:20]))
            def tail(M):
                R=M-1
                assert R*R>3/(2*sigma*sigma)
                return math.exp(-sigma*sigma*R*R)*(R*R/sigma**2+1/sigma**4)/(4*delta*delta)
            assert 0<=value-low<=tail(20)+1e-14
            bound=sum((m+1)*m*m/(4*delta*delta)*math.exp(-sigma*sigma*(m+2*delta)**2) for m in range(1,128))+tail(128)
            assert 0<=value<=bound
            if cut==math.pi/2:assert value<1e-18
            rows.append({'C':C,'cut':cut,'Gaussian_width':sigma,'per_axis_flux_variance_N32':value,
                         'all_mode_remainder_upper':tail(32),'cut_uniform_upper':bound})
    return {'cases':rows,'tail_method':'Cauchy-Schwarz spatial overlaps and decreasing Gaussian polynomial integral',
            'no_global_gap_assigned_to_horizon':True}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    out={'status':'PASS','scope':'leading bosonic dipole quadratic polar source and Gaussian ordered kernel; fixed product',
         'direct_four_dimensional_source':tensor_source_check(),'canonical_reduction':reduction_checks(),
         'ordered_kernel':spectral_kernel_checks(),'smeared_Gauss_flux':smeared_Gauss_flux(),
         'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'dependencies':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in
                         ['notes/dipole-quadrupole-source.md','notes/nonabelian-dipole-constraints.md','notes/quadrupole-product-response.md','notes/lll-polar-noise.md']},
         'versions':{'sympy':s.__version__,'numpy':np.__version__},
         'limitations':['Round bosonic projection is not a consistent nonlinear 4D truncation.',
                        'Gaussian kernel is a leading expansion; cubic/quartic interactions and fermionic polarization remain.',
                        'Global Friedrichs vacuum is not the physical horizon or boundary-driven state.',
                        'Local response contacts and full renormalized stress remain.',
                        'No capture, creation, preparation, full stability or transport bound.']}
    args.output.write_bytes((json.dumps(out,indent=2)+'\n').encode())
    print(json.dumps({'status':'PASS','kernel':out['ordered_kernel']},indent=2),flush=True)
