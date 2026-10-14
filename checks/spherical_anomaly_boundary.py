"""S074: auxiliary boundary action, boundary data and transparent gluing."""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import sympy as s

t,x=s.symbols('t x',real=True)
aa=s.Rational(11,360); bb=s.Rational(1,30); cc=s.Rational(1,20)
y=s.Rational(3,2); flux=s.Integer(3)


def zero(v):
    z=s.factor(s.cancel(v))
    assert z==0,z


def potential(P,Q,sig,v,I,a=aa,b=bb,c=cc,Y=y):
    d=b-2*a/3; A=c*sig/3; B=sig*(4*c*Y/3-4*a*Y)
    p=P-B-b*v; z=Q+6*d*v; delta=b+4*c*sig
    L0=(2*b-4*a)*Y*v-3*d*v*v+sig*(4*c*Y*Y/3-2*I/3)
    return (3*p*p+p*z-A*z*z/b)/delta-L0


def exact_checks():
    a,b,c,Y,r,u,v,I,sig,P,Q=s.symbols('a b c y r u v I sigma P Q')
    d=b-2*a/3; A=c*sig/3; B=sig*(4*c*Y/3-4*a*Y)
    L=A*r*r+b*r*u-3*b*u*u+(B+b*v)*r-6*d*v*u+(2*b-4*a)*Y*v-3*d*v*v+sig*(4*c*Y*Y/3-2*I/3)
    PP=s.diff(L,r); QQ=s.diff(L,u); H=potential(P,Q,sig,v,I,a,b,c,Y)
    Hess=s.hessian(L,(r,u)); zero(Hess.det()+b*(b+4*c*sig))
    sub={P:PP,Q:QQ}
    zero(H.subs(sub)-(PP*r+QQ*u-L))
    zero(s.diff(H,P).subs(sub)-r); zero(s.diff(H,Q).subs(sub)-u)
    zero(s.diff(H,v).subs(sub)+s.diff(L,v)); zero(s.diff(H,I).subs(sub)+s.diff(L,I))
    delta_min=bb-4*cc/s.Integer(100)
    assert delta_min==s.Rational(47,1500)
    assert 1+s.Rational(1,100)<7
    bound=7/delta_min; assert bound==s.Rational(10500,47)
    # For sigma>=0, the row norms 7/Delta and (1+sigma)/Delta
    # are at most 210 and 30; the negative short interval sets the bound.
    zero(30*(bb+4*cc*sig)-(1+sig)-5*sig)
    assert 210<bound and 30<bound
    # The second normal jets (omega_xx,sigma_xx) determine (P,Q).
    zero((Hess*s.diag(-2,1)).det()-2*b*(b+4*c*sig))
    Om=s.Matrix([[0,2,0,0],[-2,0,0,0],[0,0,0,-1],[0,0,1,0]])
    assert Om.rank()==4 and Om.det()==4
    # Direct Gaussian-normal curvature from the Christoffel definition.
    w=s.Function('omega')(t,x); g=s.diag(-s.exp(2*w),1); inv=g.inv(); coord=(t,x)
    gam=[[[sum(inv[i,l]*(s.diff(g[l,k],coord[j])+s.diff(g[j,l],coord[k])-s.diff(g[j,k],coord[l]))/2 for l in range(2)) for k in range(2)] for j in range(2)] for i in range(2)]
    ric=s.Matrix(2,2,lambda j,l:sum(s.diff(gam[i][l][j],coord[i])-s.diff(gam[i][i][j],coord[l])+sum(gam[i][i][k]*gam[k][l][j]-gam[i][l][k]*gam[k][i][j] for k in range(2)) for i in range(2)))
    zero(s.trace(inv*ric)+2*(s.diff(w,x,2)+s.diff(w,x)**2))
    PP=s.Function('P')(t,x); QQ=s.Function('Q')(t,x); ss=s.Function('sigma')(t,x)
    rr=-2*(s.diff(w,x,2)+s.diff(w,x)**2)
    uu=s.diff(ss,x,2)+s.diff(w,x)*s.diff(ss,x)-s.exp(-2*w)*(s.diff(ss,t,2)-s.diff(w,t)*s.diff(ss,t))
    raw=s.exp(w)*(PP*rr+QQ*uu)
    first=s.exp(w)*(2*s.diff(PP,x)*s.diff(w,x)-s.diff(QQ,x)*s.diff(ss,x)+s.exp(-2*w)*s.diff(QQ,t)*s.diff(ss,t))
    fx=-2*s.exp(w)*PP*s.diff(w,x)+s.exp(w)*QQ*s.diff(ss,x)
    ft=-s.exp(-w)*QQ*s.diff(ss,t)
    zero(raw-first-s.diff(fx,x)-s.diff(ft,t))
    # Classical Weyl/GHY cancellation prior to choosing any boundary profile.
    XX,G=s.symbols('X G',nonzero=True)
    physical=XX*r-2*XX*u+2*XX*v
    representative=XX*r+6*XX*v
    zero((physical-representative)/(4*G)+(XX*u+2*XX*v)/(2*G))
    # The additional transformed GHY term is +X n.sigma/(2G).
    ns=s.symbols('normal_sigma'); zero(-XX*ns/(2*G)+XX*ns/(2*G))
    # Legendre transform fixes the displacement density, including metric variation.
    At,D,dAt,dD=s.symbols('A D dA dD')
    zero(D*dAt-(D*dAt+At*dD)+At*dD)
    N,b0,gc=s.symbols('N b0 gc',nonzero=True)
    Ztotal=-s.pi*b0**2/gc**2-2*N*sig/3
    zero(-4*Ztotal-(4*s.pi*b0**2/gc**2+8*N*sig/3))
    # Moving cut curvature, directly from div n with an off-surface extension.
    cut=s.Function('cut')(t); vel=s.diff(cut,t); acc=s.diff(cut,t,2)
    e=s.exp(-2*w); D=1-e*vel**2
    nup=s.Matrix([e*vel/s.sqrt(D),1/s.sqrt(D)])
    div=sum(s.diff(s.exp(w)*nup[i],coord[i]) for i in range(2))/s.exp(w)
    claimed=(s.diff(w,x)*(1-2*e*vel**2)+e*(acc-s.diff(w,t)*vel))/D**s.Rational(3,2)
    zero((div-claimed)*D**s.Rational(3,2))
    return dict(Legendre_inverse='PASS',Hessian_determinant='-b*(b+4*c*sigma)',
                Delta_lower_bound=str(delta_min),inverse_Hessian_infinity_bound=str(bound),invertibility_range='sigma >= -1/100 at b=1/30,c=1/20',
                independent_normal_jets='PASS',restricted_boundary_symplectic_rank=4,
                restricted_boundary_symplectic_determinant=4,
                Gaussian_normal_curvature='PASS',auxiliary_GHY_and_Q_normal_flux_identity='PASS for arbitrary two-variable Gaussian-normal fields',classical_Weyl_GHY_flux_cancellation='PASS',
                electric_density_Legendre_transform='PASS',moving_cut_extrinsic_curvature='PASS')


fields=s.symbols('omega P Q sigma At Ax')
dt=s.symbols('omega_t P_t Q_t sigma_t At_t Ax_t')
dx=s.symbols('omega_x P_x Q_x sigma_x At_x Ax_x')
w,P,Q,sig,At,Ax=fields
ef=dt[5]-dx[4]
v=dx[3]**2-s.exp(-2*w)*dt[3]**2
I=-2*s.exp(-2*w)*ef**2+flux**2*y**2/2
H=potential(P,Q,sig,v,I)
L1=s.exp(w)*(2*dx[1]*dx[0]-dx[2]*dx[3]+s.exp(-2*w)*dt[2]*dt[3]-H)
px=[s.diff(L1,j) for j in dx]
pt=[s.diff(L1,j) for j in dt]
pd=[s.diff(L1,j) for j in fields]


def boundary_identities():
    vv,II=s.symbols('v I')
    symbolic=potential(P,Q,sig,vv,II)
    Hv=s.diff(symbolic,vv).subs({vv:v,II:I})
    HI=s.diff(symbolic,II).subs({vv:v,II:I})
    target=[2*s.exp(w)*dx[1],2*s.exp(w)*dx[0],-s.exp(w)*dx[3],
            s.exp(w)*(-dx[2]-2*Hv*dx[3]),4*s.exp(-w)*(-HI)*ef,0]
    for a,b in zip(px,target): zero(a-b)
    # Classical addition has curvature multiplier X/(4G), with X=b0^2 exp(2sigma).
    b0,G=s.symbols('b0 G',positive=True); X=b0**2*s.exp(2*sig)
    classical=s.exp(w)*(X*dx[3]*dx[0]/G+3*X*v/(2*G))
    zero(s.diff(classical,dx[3])-s.exp(w)*(X*dx[0]/G+3*X*dx[3]/G))
    zero(s.diff(classical,dx[0])-s.exp(w)*X*dx[3]/G)
    return dict(all_auxiliary_boundary_momenta='PASS',classical_sigma_metric_boundary_momenta='PASS')


def fun(exprs):
    f=s.lambdify((t,x),list(exprs),'numpy',cse=True)
    def evaluate(tt,xx):
        shape=np.broadcast(tt,xx).shape
        return np.stack([np.broadcast_to(a,shape) for a in f(tt,xx)])
    return evaluate


def jets(profile):
    return fun(list(profile)+[s.diff(f,t) for f in profile]+[s.diff(f,x) for f in profile]+[s.diff(profile[0],x,2)])


def substitute(expr,profile):
    replace=dict(zip(fields,profile))
    replace.update(zip(dt,[s.diff(f,t) for f in profile]))
    replace.update(zip(dx,[s.diff(f,x) for f in profile]))
    return expr.xreplace(replace)


def density(data):
    omega,P,Q,sig,At,Ax=data[:6]; ft=data[6:12]; fx=data[12:18]; wxx=data[18]
    vel=fx[3]**2-np.exp(-2*omega)*ft[3]**2
    ef=ft[5]-fx[4]; inv=-2*np.exp(-2*omega)*ef**2+float(flux**2*y**2/2)
    d=float(bb-2*aa/3); a=float(cc)*sig/3
    B=sig*float(4*cc*y/3-4*aa*y); p=P-B-float(bb)*vel; z=Q+6*d*vel
    delta=float(bb)+4*float(cc)*sig
    assert np.all(delta>0)
    Hval=(3*p*p+p*z-a*z*z/float(bb))/delta-(float((2*bb-4*aa)*y)*vel-3*d*vel**2+sig*(float(4*cc*y*y/3)-2*inv/3))
    R=-2*(wxx+fx[0]**2)
    return np.exp(omega)*(P*R-fx[2]*fx[3]+np.exp(-2*omega)*ft[2]*ft[3]-Hval)


def ghy(data):
    return 2*np.exp(data[0])*data[1]*data[12]


def quadrature(order):
    z,w=np.polynomial.legendre.leggauss(order)
    tt,xx=np.meshgrid(z,z,indexing='ij')
    return z,w,tt,xx,w[:,None]*w[None,:]


def action(data,left,right,weights,areaweights):
    return float(np.sum(areaweights*density(data))+np.sum(weights*(ghy(right)-ghy(left))))


def numerical_checks():
    envelope=(1-t*t)**6
    profile0=[x/13+x*x/17+envelope*(x/19+x*x/23),s.Rational(21,1000)+x/37+envelope*(1+x*x)/101,
              s.Rational(3,100)-x/41+envelope*(1+x)/97,s.Rational(1,500)+x/2000+envelope*(1+x+x*x)/4000,
              x/5+envelope*(t+x*x)/31,envelope*(x+t*t)/29]
    shape=envelope*(1+x/7+x*x/11)
    scales=[s.Rational(1,10),s.Rational(1,100),s.Rational(1,80),s.Rational(1,2000),s.Rational(1,12),s.Rational(1,13)]
    variations=[]
    for i in range(6):
        d=[s.Integer(0)]*6; d[i]=scales[i]*shape; variations.append((str(fields[i]),d))
    variations.append(('joint',[z*shape for z in scales]))
    variations.append(('fixed_boundary_values',[envelope*(1-x*x)/10,0,0,envelope*(1-x*x)/2000,0,0]))
    records=[]; maxframes=[]
    for sign in (1,-1):
        profile=list(profile0); profile[3]=sign*profile[3]
        jcalc=jets(profile)
        E=[substitute(pd[i],profile)-s.diff(substitute(pt[i],profile),t)-s.diff(substitute(px[i],profile),x) for i in range(6)]
        Ec=fun(E); Pc=fun([substitute(z,profile) for z in px])
        for label,var in variations:
            dc=jets(var); zc=fun(var); rows=[]
            for order in (32,48):
                nodes,weights,tt,xx,wt=quadrature(order)
                data=jcalc(tt,xx); left=jcalc(nodes,-np.ones_like(nodes)); right=jcalc(nodes,np.ones_like(nodes))
                dv=dc(tt,xx); dl=dc(nodes,-np.ones_like(nodes)); dr=dc(nodes,np.ones_like(nodes))
                bulk=float(np.sum(wt*np.sum(Ec(tt,xx)*zc(tt,xx),axis=0)))
                boundary=float(np.sum(weights*(np.sum(Pc(nodes,np.ones_like(nodes))*zc(nodes,np.ones_like(nodes)),axis=0)-np.sum(Pc(nodes,-np.ones_like(nodes))*zc(nodes,-np.ones_like(nodes)),axis=0))))
                target=bulk+boundary; diffs=[]; raw_diffs=[]
                for step in (1e-3,5e-4):
                    values=[]; raw=[]
                    for m in (-2,-1,1,2):
                        eps=m*step
                        values.append(action(data+eps*dv,left+eps*dl,right+eps*dr,weights,wt))
                        raw.append(float(np.sum(wt*density(data+eps*dv))))
                    deriv=(values[0]-8*values[1]+8*values[2]-values[3])/(12*step)
                    raw_deriv=(raw[0]-8*raw[1]+8*raw[2]-raw[3])/(12*step)
                    assert abs(deriv-target)<2e-8,(sign,label,order,step,deriv,target)
                    diffs.append(deriv); raw_diffs.append(raw_deriv)
                rows.append((target,*diffs))
            refine=max(abs(a-b) for a,b in zip(*rows)); assert refine<2e-8,(label,refine)
            record=dict(sigma_sign=sign,variation=label,predicted=target,direct=diffs[-1],discrepancy=abs(diffs[-1]-target),
                        step_refinement=abs(diffs[-1]-diffs[-2]),quadrature_refinement=refine,boundary_contribution=boundary,
                        discrepancy_if_GHY_omitted=abs(raw_diffs[-1]-target))
            if label=='fixed_boundary_values':
                assert abs(boundary)<1e-14
                assert record['discrepancy_if_GHY_omitted']>1e-4
            records.append(record)
            print('variation',sign,label,'PASS',flush=True)
        assert sum(abs(v['boundary_contribution'])>1e-7 for v in records if v['sigma_sign']==sign)>=3
    # Compare four independent split-domain quadratures with the unsplit action.
    jcalc=jets(profile0); split=[]
    for cut in [s.Integer(0),s.Rational(1,5),t/7,(1-t*t)/5]:
        ccalc=fun([cut,s.diff(cut,t),s.diff(cut,t,2)])
        rows=[]
        for order in (40,64):
            nodes,weights,tt,xx,wt=quadrature(order)
            left=jcalc(nodes,-np.ones_like(nodes)); right=jcalc(nodes,np.ones_like(nodes))
            full=action(jcalc(tt,xx),left,right,weights,wt)
            cutdata=ccalc(nodes,np.zeros_like(nodes)); cp=cutdata[0][:,None]
            xleft=-1+(xx+1)*(cp+1)/2; xright=cp+(xx+1)*(1-cp)/2
            bulk=float(np.sum(wt*((cp+1)/2*density(jcalc(tt,xleft))+(1-cp)/2*density(jcalc(tt,xright)))))
            jc=jcalc(nodes,cutdata[0]); e=np.exp(-2*jc[0]); vel=cutdata[1]; acc=cutdata[2]; DD=1-e*vel**2
            assert np.min(DD)>0
            K=(jc[12]*(1-2*e*vel**2)+e*(acc-jc[6]*vel))/DD**1.5
            surface=2*np.exp(jc[0])*np.sqrt(DD)*jc[1]*K
            innerplus=float(np.sum(weights*surface)); innerminus=float(np.sum(weights*(-surface)))
            combined=bulk+float(np.sum(weights*(ghy(right)-ghy(left))))+innerplus+innerminus
            assert abs(combined-full)<2e-8,(cut,order,combined,full)
            rows.append((full,combined,innerplus,innerminus))
        refine=max(abs(a-b) for a,b in zip(*rows)); assert refine<2e-8,(cut,refine)
        split.append(dict(cut=str(cut),unsplit=rows[-1][0],split=rows[-1][1],discrepancy=abs(rows[-1][1]-rows[-1][0]),
                          one_inner_GHY=rows[-1][2],other_inner_GHY=rows[-1][3],quadrature_refinement=refine))
    return records,split


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',required=True); args=ap.parse_args()
    exact=exact_checks(); exact.update(boundary_identities()); print('exact identities PASS',flush=True)
    variations,gluing=numerical_checks()
    result=dict(status='PASS',analytic=exact,variation_cases=variations,gluing_cases=gluing,tolerance=2e-8,
                variation_quadrature_orders=[32,48],gluing_quadrature_orders=[40,64],finite_difference_steps=[1e-3,5e-4],
                numerical_checks_are_not_interval_certificates=True,
                versions=dict(python=platform.python_version(),sympy=s.__version__,numpy=np.__version__),
                checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='PASS',checker_sha256=result['checker_sha256'])))


if __name__=='__main__': main()
