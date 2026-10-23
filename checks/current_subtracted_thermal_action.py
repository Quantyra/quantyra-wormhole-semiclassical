"""Action-level subtraction and complete thermal sums; not the full metric operator."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb
import json, sys, traceback
import sympy as sp
sys.path.insert(0, str(Path(__file__).resolve().parent))
import current_coincident_dirac_action_assembly as raw
prior=raw.prior; radial=raw.radial
from flint import arb, ctx

ROOT=Path(__file__).resolve().parents[1]; STEM='current-subtracted-thermal-action'
PREVIOUS='current-controlled-solution-audit-closeout'
OWN=['notes/'+STEM+'.md','checks/current_subtracted_thermal_action.py','evidence/'+STEM+'-review.json']
ROWS=[]; NUMERICS=[]; FIXTURES=[]; CURRENT={}
SCHEDULE=[(256,64),(640,160)]
encoded=prior.encoded; hull=prior.hull

def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name, ok, detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok: raise AssertionError(name)
prior.check=check

def intersect(x, lo, hi):
    a=max(x.lower(),arb(lo)); b=min(x.upper(),arb(hi)); assert a<=b
    return hull(a,b)

def exponential(a,L):
    al,ah=a.lower(),a.upper(); ll,lh=L.lower(),L.upper()
    assert al>0 and ll>0
    return hull((-ah*lh).exp().lower(),(-al*ll).exp().upper())

def linear_enclosures(a,F,I,L,E=None):
    E=exponential(a,L) if E is None else E
    q=F/a; qlo,qhi=q.lower(),q.upper(); ilo,ihi=I.lower(),I.upper()
    lo=(qlo+(ilo-qlo)*E).lower(); hi=(qhi+(ihi-qhi)*E).upper()
    end=hull(lo,hi); whole=hull(min(ilo,lo),max(ihi,hi))
    return end,whole

def bernstein_split(coeffs,t):
    row=list(coeffs); left=[row[0]]; right=[row[-1]]
    while len(row)>1:
        row=[(1-t)*a+t*b for a,b in zip(row,row[1:])]
        left.append(row[0]); right.append(row[-1])
    return left,list(reversed(right))

def bernstein_range(coeffs,lo,hi):
    assert lo>=0 and hi<=1 and lo<1 and hi>=lo
    _,right=bernstein_split(coeffs,lo)
    left,_=bernstein_split(right,(hi-lo)/(1-lo))
    return hull(min(x.lower() for x in left),max(x.upper() for x in left))

def derivative_coeffs(cs):
    return [(len(cs)-1)*(b-a) for a,b in zip(cs,cs[1:])]

def basis(z,f):
    if z.upper()<=arb('.5') or z.lower()>=arb('.875'):
        return [(arb(0),arb(0))]*3
    w=intersect((z-arb('.5'))/arb('.375'),0,1)
    cs=[arb(0) for _ in range(13)]; cs[6]=arb(4096)/comb(12,6)
    ds=[arb(0) for _ in range(14)]; ds[7]=arb(4096)/comb(13,7)
    lo=max(arb(0),w.lower()); hi=min(arb(1),w.upper())
    u=intersect(bernstein_range(cs,lo,hi),0,1); uw=bernstein_range(derivative_coeffs(cs),lo,hi)
    v=intersect(bernstein_range(ds,lo,hi),0,1); vw=bernstein_range(derivative_coeffs(ds),lo,hi)
    up=-f*uw/arb('.375'); vp=-f*vw/arb('.375')
    return [(u,up),(v,vp),(u+v,up+vp)]

def grid(par,N):
    nodes,_=radial.mesh(par,N); out=[]; cells=[]
    beta,u,e,t=par['beta'],par['u'],par['eps'],par['tau']
    mono=[beta,arb(0),u,-4*beta-2*u-t,3*beta+u+e,-2*e,e]
    cs=[sum((mono[i]*comb(k,i)/comb(6,i) for i in range(k+1)),arb(0)) for k in range(7)]
    ds=derivative_coeffs(cs)
    for a,b in zip(nodes,nodes[1:]):
        z=hull(b.lower(),a.upper()); lo=max(arb(0),z.lower()); hi=min(arb(1),z.upper())
        Q=bernstein_range(cs,lo,hi); Qp=bernstein_range(ds,lo,hi)
        f=Q.sqrt(); fp=-Qp/2; us=basis(z,f)
        L=(a-b)/f; cells.append((1/f,L))
        assert f>arb('.05') and f<arb('.45') and fp.abs_upper()<arb('.15')
        if a.upper()<=arb('.875'): assert Q>arb('.01')
        assert all(u.abs_upper()<=2 and up.abs_upper()<100 for u,up in us)
        out.append((f,fp,L,us))
    return nodes,cells,out

def coefficients(M,p,f,fp,u,up):
    D=p*p+(M*f)**2; lam=(M*f)**2/D; K=D.sqrt()/f
    c=-p*M*fp/(2*D)
    k1=K*lam*u; k2=K*(2*lam-lam*lam)*u*u
    c1=-p*M/(2*D)*(f*up+fp*u*(1-2*lam))
    c2=-p*M/(2*D)*((2-4*lam)*f*u*up+(1-8*lam+8*lam*lam)*fp*u*u)
    return K,c,k1,k2,c1,c2

def integrate(M,p,cells):
    X=arb(0,1); Ys=[arb(0),arb(0)]; forward=[]; active=False
    for f,fp,L,us in cells:
        D=p*p+(M*f)**2; K=D.sqrt()/f; c=-p*M*fp/(2*D)
        A=2*K+c*arb(0,1)
        _,whole=linear_enclosures(A,-c,X,L)
        Xw=intersect(whole,-1,1)
        A=2*K+c*Xw
        Xnext,refined=linear_enclosures(A,-c,X,L)
        Xw=hull(max(Xw.lower(),refined.lower()),min(Xw.upper(),refined.upper()))
        Xnext=intersect(Xnext,-1,1)
        active=active or any(not u.is_zero() or not up.is_zero() for u,up in us[:2])
        if not active:
            X=Xnext
            continue
        damping=2*K+2*c*Xw; xx=Xw**2; den=1+xx
        E=exponential(damping,L); lam=(M*f)**2/D; factor=-p*M/(2*D)
        ys=[]; ks=[]; cs=[]
        for j,(u,up) in enumerate(us[:2]):
            k1=K*lam*u
            c1=factor*(f*up+fp*u*(1-2*lam))
            FY=-2*k1*Xw-c1*(1+xx)
            Yn,Yw=linear_enclosures(damping,FY,Ys[j],L,E)
            Ys[j]=Yn; ys.append(Yw); ks.append(k1); cs.append(c1)
        forward.append((K,c,damping,E,lam,factor,Xw,ys,ks,cs,f,fp,L,us))
        X=Xnext
    adj=arb(0); gradients=[arb(0),arb(0)]; matrix=[arb(0) for _ in range(3)]
    for K,c,a,E,lam,factor,X,ys,ks,cs,f,fp,L,us in reversed(forward):
        xx0=X**2; xx=hull(max(arb(0),xx0.lower()),xx0.upper())
        den=1+xx; assert den>0
        g=xx/den; gp=2*X/den**2; gpp=2*(1-3*xx)/den**3
        adj,aw=linear_enclosures(a,4*K*gp,adj,L,E)
        C0=4*g-2*aw*X; C1=4*gp-2*aw; C2=4*K*gpp-2*aw*c
        for i in range(2): gradients[i]+=(C0*ks[i]-aw*cs[i]*(1+xx))*L
        for entry,(i,j) in enumerate([(0,0),(0,1),(1,1)]):
            ui,upi=us[i]; uj,upj=us[j]
            kij=K*(2*lam-lam*lam)*ui*uj
            cij=factor*((1-2*lam)*f*(ui*upj+uj*upi)+(1-8*lam+8*lam*lam)*fp*ui*uj)
            matrix[entry]+=(C0*kij-aw*cij*(1+xx)+C1*(ks[i]*ys[j]+ks[j]*ys[i])-2*aw*X*(cs[i]*ys[j]+cs[j]*ys[i])+C2*ys[i]*ys[j])*L
    return gradients,matrix

def blocks(par,B):
    # Exact integer half-open blocks cover every j from zero to J-1.
    targets=[arb(i)/B for i in range(1,B+1)]
    for k in range(20):
        count=B if k<5 else (max(4,B//4) if k<7 else 4)
        targets += [(arb(2).log()*(arb(k)+arb(i)/count)).exp() for i in range(1,count+1)]
    targets[-1]=arb(2)**20
    indices=[0]+[int((x/par['gamma_cl']).floor().unique_fmpz()) for x in targets]
    assert all(b>a for a,b in zip(indices,indices[1:]))
    out=[]
    for a,b in zip(indices,indices[1:]):
        lo=par['gamma']*(arb(a)+arb('.5')); hi=par['gamma']*(arb(b)-arb('.5'))
        out.append((a,b,hull(lo.lower(),hi.upper()),par['gamma']*(b-a)))
    return out,indices[-1]

def identities_and_tail():
    z=sp.symbols('z'); b,u,e,t=sp.symbols('b u e t')
    mono=[b,0,u,-4*b-2*u-t,3*b+u+e,-2*e,e]
    cs=[sum(mono[i]*sp.Rational(comb(k,i),comb(6,i)) for i in range(k+1)) for k in range(7)]
    poly=sum(cs[k]*comb(6,k)*z**k*(1-z)**(6-k) for k in range(7))
    check('Bernstein_geometry_identity',sp.expand(poly-((1-z)**2*(b*(1+2*z+3*z*z)+u*z*z+e*z**4)-t*z**3))==0)
    check('Bernstein_compact_direction_identity',sp.expand(4096*z**6*(1-z)**6-(4*z*(1-z))**6)==0)
    linear=bernstein_range([arb(0),arb(1)],arb('.25'),arb('.625'))
    NUMERICS.append(dict(kind='Bernstein_linear_fixture',enclosure=encoded(linear),exact_range=['.25','.625']))
    check('Bernstein_restriction_linear_fixture',linear.contains(arb('.25')) and linear.contains(arb('.625')) and linear.lower()>arb('.24999999') and linear.upper()<arb('.62500001'),encoded(linear))
    X,e,K,c=sp.symbols('X e K c',real=True)
    U=(X+e)/(1-X*e); m=K*(1-e*e)/(1+e*e); W=2*K*e/(1+e*e)
    transformed=sp.diff(U,X)*(-2*K*X-c*(1+X*X))+sp.diff(U,e)*c*(1+e*e)
    check('exact_rotated_Riccati',sp.cancel(transformed-(W*(1-U*U)-2*m*U))==0)
    growth=(m*(1-U*U)+2*W*U)/(1+U*U)
    check('exact_action_subtraction',sp.cancel(-2*growth+2*K-4*K*X*X/(1+X*X))==0)
    t,M,p,f,fp,u,up=sp.symbols('t M p f fp u up',positive=True)
    mt=M*sp.exp(t*u); Dt=p*p+mt*mt*f*f; kt=sp.sqrt(Dt)/f
    ct=-p*mt*(fp+t*f*up)/(2*Dt); D=p*p+M*M*f*f; lam=M*M*f*f/D; kk=sp.sqrt(D)/f
    for name,expr,want in [
        ('K1',sp.diff(kt,t).subs(t,0),kk*lam*u),
        ('K2',sp.diff(kt,t,2).subs(t,0),kk*(2*lam-lam*lam)*u*u),
        ('c1',sp.diff(ct,t).subs(t,0),-p*M/(2*D)*(f*up+fp*u*(1-2*lam))),
        ('c2',sp.diff(ct,t,2).subs(t,0),-p*M/(2*D)*((2-4*lam)*f*u*up+(1-8*lam+8*lam*lam)*fp*u*u))]:
        check('coefficient_'+name,sp.simplify(expr-want)==0)
    xx,Y,Z,k1,k2,c1,c2=sp.symbols('xx Y Z k1 k2 c1 c2')
    variables=[xx,K,c,Y,k1,c1]; shifts=[Y,k1,c1,Z,k2,c2]
    F=-2*K*xx-c*(1+xx*xx)
    delta=lambda ex:sp.expand(sum(sp.diff(ex,a)*b for a,b in zip(variables,shifts)))
    one=delta(F); two=delta(one)
    check('sensitivity_Y',sp.expand(one-(-(2*K+2*c*xx)*Y-2*k1*xx-c1*(1+xx*xx)))==0)
    check('sensitivity_Z',sp.expand(two-(-(2*K+2*c*xx)*Z-4*k1*Y-2*k2*xx-c2*(1+xx*xx)-4*c1*xx*Y-2*c*Y*Y))==0)
    ell=sp.symbols('ell'); gg=xx*xx/(1+xx*xx); gp=sp.diff(gg,xx); gpp=sp.diff(gg,xx,2)
    adjoint_integrand=4*k2*gg+8*k1*gp*Y+4*K*gpp*Y*Y+ell*(-2*k2*xx-c2*(1+xx*xx)-4*k1*Y-4*c1*xx*Y-2*c*Y*Y)
    grouped=(4*gg-2*ell*xx)*k2-ell*c2*(1+xx*xx)+(4*gp-2*ell)*2*k1*Y-4*ell*xx*c1*Y+(4*K*gpp-2*ell*c)*Y*Y
    check('exact_action_adjoint_grouping',sp.cancel(adjoint_integrand-grouped)==0)
    A=arb; m=A(26); f=A('.45'); fp=A('.15'); U=A(2); V=A(100); P=A(100)
    C1=m*(V*f+U*fp)/2+m**3*U*fp*f*f/P**2
    C2=m*(2*U*V*f+U*U*fp)/2+2*m**3*f*f*U*(V*f+U*fp)/P**2+2*m**3*fp*f*f*U*U/P**2+4*m**5*fp*f**4*U*U/P**4
    check('explicit_c1_c2_tail_constants',C1<600 and C2<2500)
    check('X_prefix_and_invariant',(-P).exp()+A('.5')/P**2<1/P**2 and 2/A('.45')-2/P**2>4)
    check('Y_tail_forcing',2*609/P**2+600*(1+P**-4)<601 and A(601)/4<151)
    zforcing=4*609*151/P**2+2*2434/P**2+2500*(1+P**-4)+4*600*151/P**4+4*151**2/P**4
    check('Z_tail_forcing',zforcing<2540 and A(2540)/4<636)
    check('Hessian_tail_coefficient',4*(A(370270)/P**2+984354)<A('4e6') and A('3.6e7')*A('1.5')<A('6e7'))
    check('gradient_tail_coefficient',9*4*(609/P**2+21*2*151)<A('6e7'))
    FIXTURES.append(dict(fixture='exact_action_and_derivative_identities_and_infinite_tail',full_operator=False))

def independent_case(par,mesh):
    nodes,cells,full=mesh; n=1; p=arb(1); M=((par['q']+1)/par['C']).sqrt()
    _,ours=integrate(M,p,full); independent=raw.discrete_characteristic_matrix(par,n,p,(nodes,cells))
    flats=[]
    for direction in range(3):
        total=arb(0)
        for k,(iq,L) in enumerate(cells):
            us=raw.directions(((nodes[k]+nodes[k+1])/2).mid()); u=(us[direction] if direction<2 else us[0]+us[1]).mid()
            mm=raw.Jet2(M.mid(),M.mid()*u,M.mid()*u*u); kk=(mm*mm+(p*iq).mid()**2).sqrt()
            total-=2*kk.dd*L.mid()
        flats.append(total)
    flat=[flats[0],(flats[2]-flats[0]-flats[1])/2,flats[1]]
    ref=[independent[0][0]-flat[0],independent[0][1]-flat[1],independent[1][1]-flat[2]]
    check('independent_raw_minus_flat_action',all(x.contains(y) for x,y in zip(ours,ref)))
    NUMERICS.append(dict(kind='independent_action_check',subtracted_continuum=[encoded(x) for x in ours],piecewise_raw_minus_flat=[encoded(x) for x in ref],piecewise_is_continuum_certificate=False))

def calculation():
    par=prior.parameters(); meshes={N:grid(par,N) for N,B in SCHEDULE}
    check('whole_path_and_direction_guards',True)
    independent_case(par,meshes[SCHEDULE[-1][0]])
    for n in [1,2]:
        M=(n*(n+par['q'])/par['C']).sqrt(); check('mass_tail_scope_'+str(n),M>18 and M<26)
        previous=None
        for N,B in SCHEDULE:
            bins,J=blocks(par,B); total=[arb(0) for _ in range(5)]; records=[]
            CURRENT.clear(); CURRENT.update(kind='partial_thermal_progress',angular_n=n,radial_N=N,bins=records)
            for i,(a,b,p,w) in enumerate(bins):
                grad,mat=integrate(M,p,meshes[N][2]); vals=grad+mat
                CURRENT.update(current_j_start=str(a),current_j_end_exclusive=str(b),current_values=[encoded(x) for x in vals])
                assert all(x.is_finite() for x in vals), 'nonfinite thermal action entry'
                total=[x+w*y for x,y in zip(total,vals)]
                records.append(dict(j_start=str(a),j_end_exclusive=str(b),p=encoded(p),weight=encoded(w),gradient=[encoded(v) for v in grad],hessian=[encoded(v) for v in mat]))
                if i%48==0: print(json.dumps(dict(n=n,radial_N=N,thermal_blocks_done=i,total_blocks=len(bins))),flush=True)
            pj=par['gamma']*(arb(J)+arb('.5')); tail=arb('6e7')*(1/(2*pj*pj)+par['gamma']/pj**3)
            total=[v+arb(0,tail.abs_upper()) for v in total]
            NUMERICS.append(dict(kind='full_subtracted_thermal_action',angular_n=n,radial_N=N,thermal_blocks_per_factor_two=B,gradient=[encoded(x) for x in total[:2]],hessian=[encoded(x) for x in total[2:]],infinite_tail_upper=encoded(tail),first_omitted_j=str(J),first_omitted_p=encoded(pj),bins=records,full_discrete_thermal_spectrum_included=True,full_angular_sum_included=False,physical_restorations_complete=False))
            check(f'infinite_thermal_tail_{n}_{N}',pj>100 and tail<arb('3e-5'))
            if previous is not None:
                check(f'refinement_{n}_{N}',all(x.overlaps(y) and x.rad()<y.rad() for x,y in zip(total,previous)))
                check(f'final_matrix_absolute_width_{n}',all(x.upper()-x.lower()<arb('.02') for x in total[2:]))
            previous=total
            print(json.dumps(dict(n=n,N=N,status='complete_thermal_spectrum_enclosed',matrix=[x.str(10) for x in total[2:]])),flush=True)
            CURRENT.clear()
    FIXTURES.append(dict(fixture='continuous_signed_flow_and_complete_discrete_thermal_sum',full_physical_correction=False))

def main():
    ctx.prec=256; review=json.loads((ROOT/OWN[2]).read_bytes()); old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}; deps=prior.dependency_identity(); error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_audit_closed',old['status']=='passed' and old['audit_closed'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        identities_and_tail(); calculation()
    except Exception:
        error=traceback.format_exc()
        if CURRENT: NUMERICS.append(dict(CURRENT))
    after={p:digest(p) for p in paths}; depafter=prior.dependency_identity(); status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),exact_action_subtraction_implemented=status=='passed',complete_thermal_difference_sum_enclosed=status=='passed',full_thermal_angular_sum_evaluated=False,complete_renormalized_operator=False,physical_restorations_complete=False,actual_physical_correction_computed=False,raw_C5_admission_certified=False,nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json'); result=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}; result['runs'].append(run)
    out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(result['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]: log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
    if status!='passed': raise SystemExit(1)

if __name__=='__main__': main()
