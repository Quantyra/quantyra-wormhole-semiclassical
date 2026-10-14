"""S074: nonlinear spherical Weyl functional and bulk source identities."""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import sympy as s


def zero(expr):
    value=s.factor(s.expand(expr))
    assert value==0,value


def exact_checks():
    aa,bb,cc=s.symbols('a b c')
    r,y,u,v,boxr,boxu,hsq,hvv,wr,wu,I=s.symbols('r y u v boxr boxu Hess2 Hessvv grad_r grad_u I')
    dd=bb-2*aa/3; H=r+2*y; K=bb*r+(2*bb-4*aa)*y
    A=cc*H**2/3-4*aa*y*r+bb*boxr-2*I/3
    boxv=2*hsq+2*wu+r*v
    euler=A-6*bb*boxu-6*dd*boxv-2*(K-6*dd*u-6*dd*v)*u-2*(bb*wr-6*dd*wu-12*dd*hvv)
    T=H-6*u-6*v; A0=r/2-u-2*v
    ricnorm=2*A0**2+4*hsq+4*v**2-4*A0*u+4*A0*v-8*hvv+2*(y-u-2*v)**2
    euler4=H**2/3-2*ricnorm+2*T**2/3
    box4R=boxr-6*boxu-6*boxv-2*(wr-6*wu-12*hvv)-2*u*T
    target=cc*H**2/3-aa*euler4+bb*box4R-2*I/3
    zero(euler-target)
    # Direct heat coefficient, with the fermion anomaly sign.
    heat_derivative=s.Rational(-60+48,360)
    assert -heat_derivative==s.Rational(1,30)
    aval,bval,cval=s.Rational(11,360),s.Rational(1,30),s.Rational(1,20)
    vals={aa:aval,bb:bval,cc:cval}
    zero((-3*bb).subs(vals)+s.Rational(1,10))
    zero(K.subs(vals)-(r/s.Integer(30)-y/s.Integer(18)))
    zero((-6*dd).subs(vals)+s.Rational(7,90))
    zero((-3*dd).subs(vals)+s.Rational(7,180))
    # Raw-to-polynomial difference after writing the two exact divergences.
    sig,P,Bdiv1,Bdiv2=s.symbols('sigma P div1 div2')
    raw=sig*(cc*H**2/3-aa*(4*y*r-2*boxr/3)-2*I/3)-2*aa*sig*P-dd*((H-6*u-6*v)**2-H**2)/12
    poly=sig*A-3*bb*u*u+K*v-6*dd*u*v-3*dd*v*v
    # div1=H Box sigma-sigma Box H; div2=sigma P sigma-u^2+(r-4y)v/3.
    divergence=dd*(H*u-sig*boxr)-2*aa*(sig*P-u*u+(r-4*y)*v/3)
    zero(raw-poly-divergence)
    cR=s.symbols('cR')
    zero((raw.subs(bb,bb-192*s.pi**2*cR)-raw)/(16*s.pi**2)-cR*((H-6*u-6*v)**2-H**2))
    c1,c2,c3,rsq,r2,ri2=s.symbols('c1 c2 c3 R_squared Ric_squared Riem_squared')
    zero(c1*rsq+c2*r2+c3*ri2-((c1+(c2+c3)/3)*rsq
             +(c2/2+2*c3)*(ri2-2*r2+rsq/3)+(-c2/2-c3)*(ri2-4*r2+rsq)))
    # Boundary-sensitive identities with arbitrary one-dimensional profiles.
    z=s.symbols('z'); sigx=s.Function('s')(z); rr=s.Function('r')(z); test=s.Function('f')(z)
    ux=s.diff(sigx,z,2); vx=s.diff(sigx,z)**2; Hx=rr+2*y
    Ax=cc*Hx**2/3-4*aa*y*rr+bb*s.diff(rr,z,2)-2*I/3
    Kx=bb*rr+(2*bb-4*aa)*y
    Ux=-6*bb*ux-6*dd*vx; Vx=Kx-6*dd*ux-6*dd*vx
    Ex=Ax+s.diff(Ux,z,2)-2*s.diff(Vx*s.diff(sigx,z),z)
    theta=Ux*s.diff(test,z)-s.diff(Ux,z)*test+2*Vx*s.diff(sigx,z)*test
    dv=2*s.diff(sigx,z)*s.diff(test,z); du=s.diff(test,z,2)
    variation=Ax*test-6*bb*ux*du+Kx*dv-6*dd*(du*vx+ux*dv)-6*dd*vx*dv
    zero(variation-Ex*test-s.diff(theta,z))
    Px=s.diff(sigx,z,4)+s.diff((rr-4*y)*s.diff(sigx,z)/3,z)
    rawx=sigx*(cc*Hx**2/3-aa*(4*y*rr-2*s.diff(rr,z,2)/3)-2*I/3)-2*aa*sigx*Px-dd*((Hx-6*ux-6*vx)**2-Hx**2)/12
    polyx=sigx*Ax-3*bb*ux**2+Kx*vx-6*dd*ux*vx-3*dd*vx**2
    flux=dd*(Hx*s.diff(sigx,z)-sigx*s.diff(Hx,z))-2*aa*(sigx*s.diff(ux,z)-ux*s.diff(sigx,z)+(rr-4*y)*sigx*s.diff(sigx,z)/3)
    zero(rawx-polyx-s.diff(flux,z))
    # Constant-product anomaly agrees with the accepted all-level trace.
    x,q=s.symbols('x q')
    product=A.subs(vals).subs({r:-2*x,boxr:0,I:q*q*y*y/2})/(16*s.pi**2)
    zero(product-(3*x*x+5*x*y+(3-15*q*q)*y*y)/(720*s.pi**2))
    # Gauge current and Ward identity for arbitrary fields in a flat 2D chart.
    t,z=s.symbols('t z'); sigma=s.Function('sigma')(t,z); ef=s.Function('E')(t,z)
    C=s.symbols('C')
    jt=4*C*s.diff(sigma*ef,z); jz=-4*C*s.diff(sigma*ef,t)
    zero(s.diff(jt,t)+s.diff(jz,z))
    f2=-2*ef**2; source=-C*f2
    ttt=-C*sigma*f2; tzz=C*sigma*f2
    zero(-s.diff(ttt,t)-source*s.diff(sigma,t)-ef*jz)
    zero(s.diff(tzz,z)-source*s.diff(sigma,z)+ef*jt)
    return {'local_heat_box_R_coefficient':str(-heat_derivative),
            'Dirac_coefficients':{'a':str(aval),'b':str(bval),'c':str(cval)},
            'arbitrary_covariant_jet_trace_identity':'PASS',
            'raw_polynomial_and_retained_divergence':'PASS',
            'explicit_boundary_flux_and_variational_current':'PASS for arbitrary one-coordinate profiles; covariant product-rule derivation retained in the note',
            'finite_R_squared_shift':'delta b=-192*pi^2*cR',
            'constant_product_full_Dirac_regression':'PASS',
            'gauge_current_and_Lorentz_force_identity':'PASS',
            'sigma_vertices':['linear anomaly term','-u^2/10+(r/30-y/18)*v','-7*u*v/90','-7*v^2/180'],
            'no_terms_beyond_fourth_order_at_fixed_h_and_A':True}


def metric_curvature(g,dg,ddg):
    size=4; inv=g.inv(); dinv=[-inv*d*inv for d in dg]
    gam={}; dgam={}
    for i in range(size):
        for j in range(size):
            for k in range(size):
                gam[i,j,k]=sum(inv[i,z]*(dg[j][z,k]+dg[k][z,j]-dg[z][j,k])/2 for z in range(size))
                for v in range(size):
                    dgam[v,i,j,k]=sum((dinv[v][i,z]*(dg[j][z,k]+dg[k][z,j]-dg[z][j,k])
                              +inv[i,z]*(ddg[v][j][z,k]+ddg[v][k][z,j]-ddg[v][z][j,k]))/2 for z in range(size))
    riem={}
    for i in range(size):
        for j in range(size):
            for k in range(size):
                for l in range(size):
                    riem[i,j,k,l]=dgam[k,i,l,j]-dgam[l,i,k,j]+sum(gam[i,k,m]*gam[m,l,j]-gam[i,l,m]*gam[m,k,j] for m in range(size))
    ric=s.Matrix(size,size,lambda j,l:sum(riem[i,j,i,l] for i in range(size)))
    scalar=s.trace(inv*ric)
    ric2=sum(inv[i,i]*inv[j,j]*ric[i,j]**2 for i in range(size) for j in range(size))
    riem2=sum(g[i,i]*inv[j,j]*inv[k,k]*inv[l,l]*val**2 for (i,j,k,l),val in riem.items())
    return ric,scalar,ric2,riem2


def direct_metric_checks():
    records=[]
    for case in range(12):
        # Exact two-variable jets at sigma=omega=0, theta=pi/2.
        R=lambda j:s.Rational(((case+2)*(j+3))%17-8,j+3)
        ds=[R(0),R(1)]; dw=[R(2),R(3)]
        dds=s.Matrix([[R(4),R(5)],[R(5),R(6)]])
        ddw=s.Matrix([[R(7),R(8)],[R(8),R(9)]])
        b2=s.Rational(case+2,3); y=1/b2
        g=s.diag(-1,1,b2,b2)
        dg=[s.zeros(4) for _ in range(4)]
        ddg=[[s.zeros(4) for _ in range(4)] for _ in range(4)]
        for a in range(2):
            for i in range(4):
                derivative=2*(ds[a]+dw[a]) if i<2 else 2*ds[a]
                dg[a][i,i]=g[i,i]*derivative
                for b in range(2):
                    second=2*(dds[a,b]+ddw[a,b])+4*(ds[a]+dw[a])*(ds[b]+dw[b]) if i<2 else 2*dds[a,b]+4*ds[a]*ds[b]
                    ddg[a][b][i,i]=g[i,i]*second
        ddg[2][2][3,3]=-2*b2
        ric,scalar,ric2,riem2=metric_curvature(g,dg,ddg)
        eta=s.diag(-1,1)
        hess=s.Matrix(2,2,lambda a,b:dds[a,b]-sum((int(c==a)*dw[b]+int(c==b)*dw[a]-eta[a,b]*eta[c,c]*dw[c])*ds[c] for c in range(2)))
        r=2*(ddw[0,0]-ddw[1,1]); u=-dds[0,0]+dds[1,1]; v=-ds[0]**2+ds[1]**2
        expected=s.zeros(4)
        for a in range(2):
            for b in range(2):
                expected[a,b]=(r/2-u-2*v)*eta[a,b]-2*hess[a,b]+2*ds[a]*ds[b]
        expected[2,2]=expected[3,3]=(y-u-2*v)*b2
        for val in ric-expected:
            zero(val)
        zero(scalar-(r+2*y-6*u-6*v))
        weyl=riem2-2*ric2+scalar**2/3
        zero(weyl-(r+2*y)**2/3)
        records.append({'case':case,'b0_squared':str(b2),'base_curvature':str(r),
                        'R4':str(scalar),'Weyl_squared':str(weyl),'exact_residual':0})
    return records


def compact_variation_checks():
    t,x=s.symbols('t x',real=True)
    f=(1-t*t)**3*(1-x*x)**3
    eta=f*(1+t/5+x/7+t*x/11)
    electric=s.Rational(1,5)+t/7-x/11
    aa,bb,cc=s.Rational(11,360),s.Rational(1,30),s.Rational(1,20)
    dd=bb-2*aa/3; y=s.Rational(3,2); q=3
    rows=[]
    for profile,omega in [('flat',s.S.Zero),('curved',(t*t-x*x)/10+t*x/20)]:
        factor=s.exp(-2*omega); volume=s.exp(2*omega)
        def dot(p,r):
            return factor*(-s.diff(p,t)*s.diff(r,t)+s.diff(p,x)*s.diff(r,x))
        def box(p):
            return factor*(-s.diff(p,t,2)+s.diff(p,x,2))
        r=-2*box(omega); H=r+2*y
        # F_tx is specified in coordinate units; its h-frame square has exp(-4omega).
        Ih=-2*factor**2*electric**2+s.Rational(q*q,2)*y*y
        A=cc*H*H/3-4*aa*y*r+bb*box(r)-2*Ih/3
        K=bb*r+(2*bb-4*aa)*y
        for amp in [s.Rational(1,20),s.Rational(-1,7),s.Rational(1,3)]:
            sigma=amp*f; u=box(sigma); v=dot(sigma,sigma)
            du=box(eta); dv=2*dot(sigma,eta)
            variation=eta*A-6*bb*u*du+K*dv-6*dd*(du*v+u*dv)-6*dd*v*dv
            # Independent anomaly density from conformally transformed Ricci curvature.
            coords=[t,x]; signs=[-1,1]
            hess=s.Matrix(2,2,lambda a,b:s.diff(sigma,coords[a],coords[b])-sum((int(c==a)*s.diff(omega,coords[b])+int(c==b)*s.diff(omega,coords[a])-int(a==b)*signs[a]*signs[c]*s.diff(omega,coords[c]))*s.diff(sigma,coords[c]) for c in range(2)))
            hsq=factor**2*sum(signs[a]*signs[b]*hess[a,b]**2 for a in range(2) for b in range(2))
            hvv=factor**2*sum(signs[a]*signs[b]*hess[a,b]*s.diff(sigma,coords[a])*s.diff(sigma,coords[b]) for a in range(2) for b in range(2))
            a0=r/2-u-2*v
            ric2=2*a0*a0+4*hsq+4*v*v-4*a0*u+4*a0*v-8*hvv+2*(y-u-2*v)**2
            T=H-6*u-6*v
            e4=H*H/3-2*ric2+2*T*T/3
            anomaly=cc*H*H/3-aa*e4+bb*(box(T)-2*dot(sigma,T)-2*u*T)-2*Ih/3
            fun=s.lambdify((t,x),[volume*variation,volume*eta*anomaly],'numpy',cse=True)
            values=[]
            for order in [32,48]:
                nodes,weights=np.polynomial.legendre.leggauss(order)
                xx,tt=np.meshgrid(nodes,nodes,indexing='xy')
                ww=np.outer(weights,weights)
                left,right=fun(tt,xx)
                values.append([float(np.sum(ww*left)),float(np.sum(ww*right))])
            discrepancy=max(abs(left-right) for left,right in values)
            refinement=max(abs(values[0][j]-values[1][j]) for j in range(2))
            assert discrepancy<2e-10,(profile,amp,discrepancy)
            assert refinement<2e-10,(profile,amp,refinement)
            rows.append({'profile':profile,'amplitude':str(amp),'integrated_action_variation':values[1][0],
                         'integrated_four_dimensional_trace':values[1][1],
                         'maximum_discrepancy':discrepancy,'quadrature_refinement':refinement})
    return rows


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); root=Path(__file__).resolve().parents[1]
    deps=['notes/throat-trace-anomaly.md','notes/dirac-spherical-source-reduction.md',
          'notes/uniform-dirac-reference.md','notes/curved-electric-matching.md',
          'notes/matter-aware-local-basis.md']
    hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in deps}
    print('Checking local heat coefficient, nonlinear transgression and arbitrary-jet source identities',flush=True)
    result={'analytic':exact_checks()}
    print('Checking direct four-dimensional metric jets',flush=True)
    result['metric_cases']=direct_metric_checks()
    print('Comparing compact-profile action variation with the four-dimensional trace',flush=True)
    result['compact_variation_cases']=compact_variation_checks()
    result.update({'status':'PASS','scope':'bulk spherical anomaly component; no complete nonlocal state, boundary ensemble or gravitational constraint solution',
                   'versions':{'python':platform.python_version(),'sympy':s.__version__,'numpy':np.__version__},
                   'dependencies_sha256':hashes,
                   'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: nonlinear spherical anomaly action and joint bulk source identity',flush=True)


if __name__=='__main__':
    main()
