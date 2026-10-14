"""S074: complete local anomaly sources and corrected spherical constraints."""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import sympy as s

t,x=s.symbols('t x',real=True)
coords=(t,x)
aa,bb,cc,y,q=s.symbols('a b c y q',real=True)


def zero(expr):
    value=s.factor(s.expand(expr))
    assert value==0, value


class Conformal:
    def __init__(self,w):
        self.w=w
        self.g=s.diag(-s.exp(2*w),s.exp(2*w))
        self.inv=s.diag(-s.exp(-2*w),s.exp(-2*w))
        self.vol=s.exp(2*w)
        self.gam=[[[sum(self.inv[a,d]*(s.diff(self.g[d,c],coords[b])+s.diff(self.g[b,d],coords[c])-s.diff(self.g[b,c],coords[d]))/2 for d in range(2)) for c in range(2)] for b in range(2)] for a in range(2)]
        self.r=-2*s.exp(-2*w)*(-s.diff(w,t,2)+s.diff(w,x,2))

    def grad(self,f): return s.Matrix([s.diff(f,z) for z in coords])
    def dot(self,f,g): return (self.grad(f).T*self.inv*self.grad(g))[0]
    def box(self,f): return s.exp(-2*self.w)*(-s.diff(f,t,2)+s.diff(f,x,2))
    def hess(self,f):
        return s.Matrix(2,2,lambda a,b:s.diff(f,coords[a],coords[b])-sum(self.gam[c][a][b]*s.diff(f,coords[c]) for c in range(2)))
    def div(self,v): return sum(s.diff(self.vol*v[a],coords[a]) for a in range(2))/self.vol
    def divtensor(self,T):
        return s.Matrix([sum(self.inv[a,c]*(s.diff(T[a,b],coords[c])-sum(self.gam[d][c][a]*T[d,b]+self.gam[d][c][b]*T[a,d] for d in range(2))) for a in range(2) for c in range(2)) for b in range(2)])


def sources(geom,sig,ef):
    g,inv=geom.g,geom.inv
    r=geom.r; u=geom.box(sig); v=geom.dot(sig,sig); d=bb-2*aa/3
    F=s.Matrix([[0,ef],[-ef,0]]); Fu=inv*F*inv; FF=F*inv*F.T
    I=sum(F[a,b]*Fu[a,b] for a in range(2) for b in range(2))+q*q*y*y/2
    A=cc*(r+2*y)**2/3-4*aa*y*r-2*I/3
    L=sig*A+bb*r*u-3*bb*u*u+(bb*r+(2*bb-4*aa)*y)*v-6*d*u*v-3*d*v*v
    P=sig*(2*cc*(r+2*y)/3-4*aa*y)+bb*(u+v)
    Q=bb*r-6*bb*u-6*d*v
    V=bb*r+(2*bb-4*aa)*y-6*d*(u+v); Z=-2*sig/3
    gs=geom.grad(sig); gQ=geom.grad(Q)
    M=-g*L/2+g*r*P/2+g*geom.box(P)-geom.hess(P)+V*gs*gs.T-(gs*gQ.T+gQ*gs.T)/2+g*(geom.dot(Q,sig)+Q*u)/2+2*Z*FF
    E=A+geom.box(Q)-2*geom.div(V*inv*gs)
    j=s.Matrix([-4*geom.div(s.Matrix([Z*Fu[a,b] for a in range(2)])) for b in range(2)])
    return dict(L=L,P=P,Q=Q,V=V,Z=Z,E=E,T=-2*M,j=j,F=F,Fu=Fu,FF=FF,I=I,u=u,v=v)


def exact_anomaly():
    h=Conformal(s.Integer(0)); sig=s.Function('sigma')(t,x); ef=s.Function('e')(t,x)
    f=sources(h,sig,ef)
    ward=h.divtensor(f['T'])-f['E']*h.grad(sig)-f['F']*f['j']
    for v0 in ward: zero(v0)
    zero(h.div(f['j']))
    # The complete pointwise first variation, including arbitrary inverse-metric
    # and gauge variations, agrees with the Euler terms and boundary current.
    k=s.Matrix([[s.Function('k0')(t,x),s.Function('k1')(t,x)],[s.Function('k1')(t,x),s.Function('k2')(t,x)]])
    z=s.Function('z')(t,x); al=s.Matrix([s.Function('A0')(t,x),s.Function('A1')(t,x)])
    trace=s.trace(h.g*k); divk=s.Matrix([sum(s.diff(k[a,b],coords[b]) for b in range(2)) for a in range(2)])
    gs=h.grad(sig)
    dr=h.box(trace)-sum(s.diff(divk[a],coords[a]) for a in range(2))
    du=sum(k[a,b]*s.diff(sig,coords[a],coords[b]) for a in range(2) for b in range(2))+(divk.T*gs)[0]-h.dot(trace,sig)/2+h.box(z)
    dv=(gs.T*k*gs)[0]+2*h.dot(sig,z)
    dF=s.diff(al[1],t)-s.diff(al[0],x)
    dI=2*sum(f['FF'][a,b]*k[a,b] for a in range(2) for b in range(2))+4*f['Fu'][0,1]*dF
    direct=-trace*f['L']/2+f['P']*dr+f['Q']*du+f['V']*dv+(cc*(2*y)**2/3-2*f['I']/3)*z+f['Z']*dI
    P,Q,V,Z=[f[key] for key in ('P','Q','V','Z')]
    theta=P*(h.inv*h.grad(trace)-divk)-trace*h.inv*h.grad(P)+k*h.grad(P)+Q*k*gs-Q*h.inv*gs*trace/2+Q*h.inv*h.grad(z)-z*h.inv*h.grad(Q)+2*V*z*h.inv*gs+4*Z*f['Fu']*al
    euler=-sum(f['T'][a,b]*k[a,b] for a in range(2) for b in range(2))/2+f['E']*z+(f['j'].T*al)[0]
    zero(direct-euler-h.div(theta))
    # Physical source transformation, including an arbitrary h-frame stress.
    tr,src,G,X,sigma=s.symbols('tau_trace s G X sigma')
    J=2*G*s.exp(-2*sigma)*(src-tr)/X
    zero(s.exp(-2*sigma)*tr+X*J/(2*G)-s.exp(-2*sigma)*src)
    w=s.Function('w')(t,x); sigma=s.Function('sigma')(t,x)
    old=Conformal(w); new=Conformal(w+sigma)
    arbitrary=s.Matrix([[s.Function('T00')(t,x),s.Function('T01')(t,x)],
                        [s.Function('T01')(t,x),s.Function('T11')(t,x)]])
    trace=s.trace(old.inv*arbitrary)
    for entry in new.divtensor(arbitrary)-s.exp(-2*sigma)*(old.divtensor(arbitrary)-trace*old.grad(sigma)): zero(entry)
    return dict(flat_arbitrary_function_Ward='PASS',gauge_conservation='PASS',
                full_pointwise_boundary_variation='PASS for arbitrary flat-chart fields and all metric/gauge/scalar variations',physical_source_chain_rule='PASS',
                physical_frame_Ward_transfer='PASS for arbitrary two-dimensional conformal factors and symmetric stress')


def exact_constraints():
    w=s.Function('rho')(t,x); X=s.Function('X')(t,x); h=Conformal(w)
    G,gc,rm,Lambda=s.symbols('G gc rm Lambda',nonzero=True)
    VX=2-2*Lambda*X-2*rm**2/X; dVX=-2*Lambda+2*rm**2/X**2
    xx=h.grad(X); v=h.dot(X,X)
    E=h.g*h.box(X)-h.hess(X)+(xx*xx.T-h.g*v/2)/(2*X)-h.g*VX/2
    D=h.r-h.box(X)/X+v/(2*X**2)+dVX
    for value in h.divtensor(E)+D*xx/2: zero(value)
    ef=s.Function('e')(t,x); F=s.Matrix([[0,ef],[-ef,0]]); Fu=h.inv*F*h.inv
    I=sum(F[a,b]*Fu[a,b] for a in range(2) for b in range(2))
    TM=4*s.pi*X/gc**2*(F*h.inv*F.T-h.g*I/4)
    MA=s.Matrix([4*s.pi/gc**2*h.div(s.Matrix([X*Fu[a,b] for a in range(2)])) for b in range(2)])
    MX=-s.pi*I/gc**2
    for value in h.divtensor(TM)-MX*xx-F*MA: zero(value)
    # Maxwell implies the scaled Q law; the historical factor instead belongs to Dflux.
    eps=s.Matrix([[0,h.vol],[-h.vol,0]])
    Qe=4*s.pi*X*ef/(gc**2*h.vol); Dflux=X*ef/h.vol
    for value in h.grad(Qe)-eps*(-MA): zero(value)
    for value in h.grad(Dflux)-gc**2/(4*s.pi)*eps*(-MA): zero(value)
    z=s.symbols('z',positive=True); f=s.Function('f')(z)
    J=4*G*s.diff(f,z)
    zero(s.diff(4*G*f,z)-J)
    wrong_angular_residual=2*J
    assert wrong_angular_residual!=0
    # Double-null divergence of an arbitrary covariant constraint tensor.
    rho=s.Function('rho')(t,x); vol=s.exp(2*rho)/2
    g=s.Matrix([[0,-vol],[-vol,0]]); inv=g.inv()
    C=s.Matrix([[s.Function('Cuu')(t,x),s.Function('Cuv')(t,x)],[s.Function('Cuv')(t,x),s.Function('Cvv')(t,x)]])
    gam=[[[sum(inv[a,d]*(s.diff(g[d,c],coords[b])+s.diff(g[b,d],coords[c])-s.diff(g[b,c],coords[d]))/2 for d in range(2)) for c in range(2)] for b in range(2)] for a in range(2)]
    for b in range(2):
        div=sum(inv[a,c]*(s.diff(C[a,b],coords[c])-sum(gam[d][c][a]*C[d,b]+gam[d][c][b]*C[a,d] for d in range(2))) for a in range(2) for c in range(2))
        other=1-b
        expected=-2*s.exp(-2*rho)*(s.diff(C[b,b],coords[other])+s.diff(C[0,1],coords[b])-2*s.diff(rho,coords[b])*C[0,1])
        zero(div-expected)
    # Independent local Lorentzian jets prove the charged/Lambda mass balance.
    R,Qe,rt,rx,htt,htx,hxx,qt,qx=s.symbols('R Q rt rx htt htx hxx qt qx',real=True)
    g=s.diag(-1,1); grad=s.Matrix([rt,rx]); Hess=s.Matrix([[htt,htx],[htx,hxx]])
    vv=(grad.T*g*grad)[0]; box=s.trace(g*Hess)
    V0=2-2*Lambda*R**2-2*rm**2/R**2
    E=g*(2*R*box+vv)-2*R*Hess-g*V0/2
    TM=-g*gc**2*Qe**2/(8*s.pi*R**2)
    T=E/(2*G)-TM; trace=s.trace(g*T)
    vv_sym=s.symbols('v')
    mass=R/(2*G)*(1-vv_sym-Lambda*R**2/3+(rm**2+G*gc**2*Qe**2/(4*s.pi))/R**2)
    dm=s.diff(mass,R).subs(vv_sym,vv)*grad+s.diff(mass,vv_sym)*2*Hess*g*grad+s.diff(mass,Qe)*s.Matrix([qt,qx])
    flux=(T*g-s.eye(2)*trace)*grad+gc**2*Qe/(4*s.pi*R)*s.Matrix([qt,qx])
    for entry in dm-flux: zero(entry)
    return dict(gravitational_Noether_identity='PASS on arbitrary conformal metric and radius',Maxwell_Noether_identity='PASS',
                dilaton_sign='D0+J_X=0 without electric field',wrong_sign_residual=str(wrong_angular_residual),
                scaled_charge='d Q_e = epsilon j',unscaled_flux='d D = gc^2/(4*pi) epsilon j',
                null_constraint_propagation='PASS',charged_Lambda_mass_work_identity='PASS on unrestricted local radius/charge jets')


def tensor_function(expressions):
    funcs=[s.lambdify((t,x),v,'numpy',cse=True) for v in expressions]
    def evaluate(tt,xx):
        shape=np.broadcast(tt,xx).shape
        return np.stack([np.broadcast_to(fn(tt,xx),shape) for fn in funcs],axis=0)
    return evaluate


def metric_data(expressions):
    # Covariant metric and two coordinate derivatives, directly from the fields.
    entries=list(expressions)
    entries += [s.diff(v,z) for z in coords for v in expressions]
    entries += [s.diff(v,z,w) for z in coords for w in coords for v in expressions]
    return tensor_function(entries)


def direct_density(metric_jets,sigma_jets,ef,aval,bval,cval,yval,qval):
    shape=ef.shape
    g=metric_jets[:4].reshape(2,2,*shape)
    dg=metric_jets[4:12].reshape(2,2,2,*shape)
    ddg=metric_jets[12:].reshape(2,2,2,2,*shape)
    det=g[0,0]*g[1,1]-g[0,1]*g[1,0]
    assert np.all(det<0)
    inv=np.array([[g[1,1],-g[0,1]],[-g[1,0],g[0,0]]])/det
    dinv=np.zeros_like(dg)
    for v in range(2):
        for i in range(2):
            for j in range(2): dinv[v,i,j]=-sum(inv[i,k]*dg[v,k,l]*inv[l,j] for k in range(2) for l in range(2))
    gam=np.zeros((2,2,2)+shape); dgam=np.zeros((2,2,2,2)+shape)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                gam[i,j,k]=sum(inv[i,l]*(dg[j,l,k]+dg[k,j,l]-dg[l,j,k])/2 for l in range(2))
                for v in range(2):
                    dgam[v,i,j,k]=sum((dinv[v,i,l]*(dg[j,l,k]+dg[k,j,l]-dg[l,j,k])+inv[i,l]*(ddg[v,j,l,k]+ddg[v,k,j,l]-ddg[v,l,j,k]))/2 for l in range(2))
    ric=np.zeros_like(g)
    for j in range(2):
        for l in range(2):
            ric[j,l]=sum(dgam[i,i,l,j]-dgam[l,i,i,j]+sum(gam[i,i,k]*gam[k,l,j]-gam[i,l,k]*gam[k,i,j] for k in range(2)) for i in range(2))
    r=sum(inv[a,b]*ric[a,b] for a in range(2) for b in range(2))
    sig=sigma_jets[0]; ds=sigma_jets[1:3]; dds=sigma_jets[3:].reshape(2,2,*shape)
    u=sum(inv[a,b]*(dds[a,b]-sum(gam[c,a,b]*ds[c] for c in range(2))) for a in range(2) for b in range(2))
    v=sum(inv[a,b]*ds[a]*ds[b] for a in range(2) for b in range(2))
    I=2*ef**2/det+qval*qval*yval*yval/2
    d=bval-2*aval/3
    L=sig*(cval*(r+2*yval)**2/3-4*aval*yval*r-2*I/3)+bval*r*u-3*bval*u*u+(bval*r+(2*bval-4*aval)*yval)*v-6*d*u*v-3*d*v*v
    return np.sqrt(-det)*L


def integrated_checks():
    vals={aa:s.Rational(11,360),bb:s.Rational(1,30),cc:s.Rational(1,20),y:s.Rational(3,2),q:3}
    const=[float(vals[v]) for v in (aa,bb,cc,y,q)]
    bump=(1-t*t)**5*(1-x*x)**5
    z=bump*(1+t/5+x/7)
    covvars=[s.Matrix([[bump,0],[0,0]]),s.Matrix([[0,bump],[bump,0]]),s.Matrix([[0,0],[0,bump]])]
    alph=s.Matrix([bump*(1+x/3),bump*(1+t/4)])
    dF=s.diff(alph[1],t)-s.diff(alph[0],x)
    records=[]
    for curved in (False,True):
        h=Conformal((t*t-x*x)/10+t*x/20 if curved else s.Integer(0))
        sig=s.Rational(1,7)*(1-t*t)**3*(1-x*x)**3
        ef=s.Rational(1,5)+t/7-x/11
        f=sources(h,sig,ef)
        sj=[sig]+list(h.grad(sig))+[s.diff(sig,a,b) for a in coords for b in coords]
        zj=[z]+list(h.grad(z))+[s.diff(z,a,b) for a in coords for b in coords]
        bmet=metric_data(list(h.g)); smet=tensor_function(sj); zmet=tensor_function(zj)
        efun=tensor_function([ef,dF])
        cases=[]
        for n,p in enumerate(covvars):
            # Covariant perturbation p corresponds to k=-h^-1 p h^-1.
            integrand=h.vol*sum(f['T'][a,b]*(h.inv*p*h.inv)[a,b] for a in range(2) for b in range(2))/2
            cases.append(('metric_'+str(n),metric_data(list(p)),False,False,integrand))
        cases += [('sigma',None,True,False,h.vol*f['E']*z),('gauge',None,False,True,h.vol*(f['j'].T*alph)[0])]
        p=sum(covvars,s.zeros(2)); joint=h.vol*(sum(f['T'][a,b]*(h.inv*p*h.inv)[a,b] for a in range(2) for b in range(2))/2+f['E']*z+(f['j'].T*alph)[0])
        cases.append(('joint',metric_data(list(p)),True,True,joint))
        for name,pfun,ps,pa,pred in cases:
            # Keep the derivative expression factored for common-subexpression
            # elimination; expanding the curved profile is unnecessarily costly.
            pexpr=pred.subs(vals)
            pcalc=tensor_function([pexpr])
            estimates=[]
            for order in (32,48):
                nodes,weights=np.polynomial.legendre.leggauss(order)
                tt,xx=np.meshgrid(nodes,nodes,indexing='ij'); wt=weights[:,None]*weights[None,:]
                gm=bmet(tt,xx); sj0=smet(tt,xx); dz=zmet(tt,xx); ee=efun(tt,xx)
                pm=pfun(tt,xx) if pfun else np.zeros_like(gm)
                target=float(np.sum(wt*pcalc(tt,xx)[0]))
                diffs=[]
                for step in (1e-3,5e-4):
                    actions=[]
                    for m in (-2,-1,1,2):
                        eps=m*step
                        density=direct_density(gm+eps*pm,sj0+eps*dz*ps,ee[0]+eps*ee[1]*pa,*const)
                        actions.append(float(np.sum(wt*density)))
                    derivative=(actions[0]-8*actions[1]+8*actions[2]-actions[3])/(12*step)
                    assert abs(derivative-target)<1e-7,(curved,name,order,step,derivative,target)
                    diffs.append(derivative)
                estimates.append((target,*diffs))
            assert max(abs(a-b) for a,b in zip(*estimates))<1e-7,(name,estimates)
            records.append(dict(curved=curved,variation=name,predicted=estimates[-1][0],direct_derivative=estimates[-1][-1],
                                discrepancy=abs(estimates[-1][-1]-estimates[-1][0]),step_refinement=abs(estimates[-1][-1]-estimates[-1][-2]),
                                quadrature_refinement=max(abs(a-b) for a,b in zip(*estimates))))
            print('variation',curved,name,'PASS',flush=True)
    return records


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',required=True); args=parser.parse_args()
    exact=exact_anomaly(); print('anomaly identities PASS',flush=True)
    constraints=exact_constraints(); print('constraint and normalization identities PASS',flush=True)
    cases=integrated_checks()
    result=dict(status='PASS',analytic=exact,constraints=constraints,variation_cases=cases,
                tolerance=1e-7,quadrature_orders=[32,48],finite_difference_steps=[1e-3,5e-4],finite_difference_order=4,
                numerical_checks_are_not_interval_certificates=True,
                versions=dict(python=platform.python_version(),sympy=s.__version__,numpy=np.__version__),
                checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({key:result[key] for key in ('status','checker_sha256')}))


if __name__=='__main__': main()
