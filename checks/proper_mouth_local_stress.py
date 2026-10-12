"""Local proper-mouth comparison; does not compute the nonlocal determinant."""
from pathlib import Path
import hashlib,json,platform,time
import mpmath as mp
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
mp.mp.dps=65


def symbolic():
    x,e,a,L=sp.symbols('x e a L',positive=True)
    b=sp.sqrt(1-e*e*x*x);A=sp.acos(e*x)
    I=2*a*b/e;J=2*a*(b/e*(L-1)+x*A)
    D=lambda v:sp.simplify(x*sp.diff(v,x))
    Ip=-2*a*e*x*x/b
    Ipp=-2*a*e*(2*x*x/b+e*e*x**4/b**3)
    R,U,V=sp.symbols('R U V',real=True,positive=True)
    f=1+R*R;H=sp.log(f);AR=sp.atan(R)
    cg=1+R*AR-U
    s=(V+2*U)*R*R-R*R-R*(R*R+3)*AR+H
    ci=(sp.Rational(1,2)-U)*R+f*AR/2
    B=AR-V*R/2-H/(2*R)
    KR=sp.Function('K')(R)
    Q=L*B+V*(R-AR)/2-AR*H/2+KR
    checks={
      'proper_I_first':D(I)-Ip,
      'proper_I_second':D(Ip)-Ipp,
      'proper_J_first':D(J)-(2*a*x*A+L*Ip),
      'proper_J_second':D(2*a*x*A+L*Ip)-(2*a*(x*A-e*x*x/b)+L*Ipp),
      'sphere_primitive':sp.diff(ci,R)-cg,
      'moving_cut_cancellation':-s/(2*R)-ci-B,
      'logarithmic_integrand_cancellation':R*ci+s/2-(V*R*R/2-R*AR+H/2),
      'logarithmic_primitive':sp.diff(Q,R).subs(sp.diff(KR,R),H/f)-L*sp.diff(B,R)-(V*R*R/2-R*AR+H/2)/f,
      'B_clock_derivative':B-(R+1/R)*sp.diff(B,R)-(AR+V/(2*R)-H/R-H/(2*R**3)),
      'xR_first':R-(R+1/R)+1/R,
      'xR_second':-1/R-(R+1/R)/R**2+2/R+1/R**3}
    checks={k:sp.simplify(v) for k,v in checks.items()}
    assert all(v==0 for v in checks.values()),checks
    return {k:str(v) for k,v in checks.items()}


def geometry(x,eps):
    z=eps*x;b=mp.sqrt(1-z*z);R=b/z
    return b,R,mp.atan(R)


def pieces(R):
    # Geometric integration cells resolve the center and long overlap equally.
    p=[mp.mpf(0)];v=mp.mpf(1)
    while v<R:p.append(v);v*=10
    p.append(R)
    return p


def product(x,eps,lam,a=1):
    b,R,A=geometry(x,eps);Lb=mp.log(lam/eps)
    return 2*a*b/eps,2*a*(b/eps*(Lb-1)+x*A)


def product_integrals(x,eps,lam,a=1):
    _,R,_=geometry(x,eps)
    return (2*a*x*mp.quad(lambda r:1,pieces(R)),
            2*a*x*mp.quad(lambda r:mp.log(lam*x*mp.sqrt(1+r*r)),pieces(R)))


def product_derivatives(x,eps,lam,a=1):
    b,R,A=geometry(x,eps);Lb=mp.log(lam/eps)
    Ip=-2*a*eps*x*x/b
    Ipp=-2*a*eps*(2*x*x/b+eps*eps*x**4/b**3)
    Jp=2*a*x*A+Lb*Ip
    Jpp=2*a*(x*A-eps*x*x/b)+Lb*Ipp
    return ((Ip,Ipp),(Jp,Jpp))


def product_checks():
    rows=[]
    for eps in map(mp.mpf,['.2','.01','.000001']):
      for x in map(mp.mpf,['.8','1.1']):
        lam=mp.mpf('.00003');direct=product_integrals(x,eps,lam)
        closed=product(x,eps,lam)
        err=max(abs(direct[i]-closed[i]) for i in [0,1])
        assert err<mp.mpf('1e-48')
        rows.append(dict(epsilon=float(eps),x=float(x),absolute_integral_error=float(err)))
    fd=[]
    for eps in map(mp.mpf,['.1','.001','.000001']):
        x=mp.mpf('.9');lam=mp.mpf('.00003');center=product_integrals(x,eps,lam)
        expected=product_derivatives(x,eps,lam);errors=[]
        for h in map(mp.mpf,['.002','.001','.0005']):
            hi=product_integrals(x*mp.exp(h),eps,lam)
            lo=product_integrals(x*mp.exp(-h),eps,lam)
            err=[]
            for k in [0,1]:
                first=(hi[k]-lo[k])/(2*h)
                second=(hi[k]-2*center[k]+lo[k])/(h*h)
                err.extend([abs(first-expected[k][0]),abs(second-expected[k][1])])
            errors.append([float(t) for t in err])
        assert max(errors[-1])<2e-6
        assert max(errors[0])/max(errors[-1])>12
        fd.append(dict(epsilon=float(eps),steps=[.002,.001,.0005],errors=errors))
    return rows,fd


def fields(r,U,V):
    f=1+r*r;cg=1+r*mp.atan(r)-U
    s=(V+2*U)*r*r-r*r-r*(r*r+3)*mp.atan(r)+mp.log(f)
    return f,cg,s


def coefficients(x,eps,lam,V):
    _,R,A=geometry(x,eps);H=mp.log(1+R*R);Lb=mp.log(lam/eps)
    K=mp.quad(lambda t:-2*mp.log(mp.cos(t)),[0,A/2,A])
    B=A-V*R/2-H/(2*R)
    Q=Lb*B+V*(R-A)/2-A*H/2+K
    return 2*x*B,2*x*Q


def corrected_integrals(alpha,x,eps,lam,U,V):
    _,R0,_=geometry(x,eps)
    target=1/(eps*x)**2
    def lapse(r):
        f,cg,s=fields(r,U,V)
        return f+alpha*s
    R=mp.findroot(lambda r:lapse(r)-target,(R0*(1-mp.mpf('.01')),R0*(1+mp.mpf('.01'))),tol=mp.mpf('1e-55'))
    def sphere(r):return 1+alpha*fields(r,U,V)[1]
    assert R>0 and sphere(0)>0 and sphere(R)>0
    return (2*x*mp.quad(lambda r:1/sphere(r),pieces(R)),
            2*x*mp.quad(lambda r:mp.log(lam*x*mp.sqrt(lapse(r)))/sphere(r),pieces(R)),R)


def first_order_checks():
    rows=[]
    for eps in map(mp.mpf,['.15','.003','.000001']):
      for U,V in [(mp.mpf('-0.333333333'),mp.mpf('0.333333333')),(mp.mpf('.7'),mp.mpf('-.4'))]:
        x=mp.mpf('.9');lam=mp.mpf('.000003');Lb=mp.log(lam/eps)
        _,R,_=geometry(x,eps);R1=-fields(R,U,V)[2]/(2*R)
        # Independently integrate the uncombined large terms.
        Bi=R1-mp.quad(lambda r:fields(r,U,V)[1],pieces(R))
        Qi=R1*Lb+mp.quad(lambda r:fields(r,U,V)[2]/(2*(1+r*r))-
          fields(r,U,V)[1]*mp.log(lam*x*mp.sqrt(1+r*r)),pieces(R))
        closed=coefficients(x,eps,lam,V)
        err=max(abs(2*x*Bi-closed[0]),abs(2*x*Qi-closed[1]))
        assert err<mp.mpf('1e-42')
        # U drops out although each separate integral depends on U.
        rows.append(dict(epsilon=float(eps),U=float(U),V_p=float(V),absolute_coefficient_error=float(err),
             radial_cut_shift=float(R1),sphere_integral=float(R1-Bi),I_coefficient=float(closed[0]),J_coefficient=float(closed[1])))
    refined=[]
    for eps in map(mp.mpf,['.15','.03']):
      x=mp.mpf('.9');lam=mp.mpf('.00003');U=mp.mpf('-.3');V=mp.mpf('.4')
      _,R,_=geometry(x,eps);base=product(x,eps,lam);co=coefficients(x,eps,lam,V)
      errors=[]
      for scale in map(mp.mpf,['.001','.0005','.00025']):
        alpha=scale/R;actual=corrected_integrals(alpha,x,eps,lam,U,V)
        residual=[abs(actual[k]-base[k]-alpha*co[k]) for k in [0,1]]
        errors.append(dict(alpha=float(alpha),remainders=[float(v) for v in residual],
                           normalized=[float(v/(alpha*alpha*R**3)) for v in residual]))
      ratios=[errors[0]['remainders'][k]/errors[-1]['remainders'][k] for k in [0,1]]
      assert min(ratios)>14 and max(ratios)<18
      refined.append(dict(epsilon=float(eps),refinements=errors,quadratic_ratios=ratios))
    return rows,refined


def collar(eps,x,lam,order):
    Lb=mp.log(lam/eps)
    def integrand(u,phi):
        xx=x*mp.exp(phi);z=eps*xx;b=mp.sqrt(1-z*z)
        rr=mp.log(1-u)-mp.log(mp.cos(z*u)-b*mp.sin(z*u)/z)
        v=u+1;chi=1-10*v**3+15*v**4-6*v**5
        cr=chi*rr;ex=mp.exp(2*cr);den=(1-u)**2
        return mp.matrix([mp.expm1(2*cr)/den,
          (mp.expm1(2*cr)*(Lb-mp.log(1-u))+ex*cr)/den])
    vals=[]
    for k in [0,1]:
        vals.append(2/eps*mp.quad(lambda u:mp.diff(lambda ph:integrand(u,ph)[k],0,order),[-1,0]))
    return vals


def collar_checks():
    rows=[]
    for eps in map(mp.mpf,['.1','.001','.00001']):
      x=mp.mpf('1.1');lam=mp.mpf('.000003')
      for j in [0,1,2]:
        val=collar(eps,x,lam,j);norm=[abs(val[0])/(eps*x*x),abs(val[1])/(eps*x*x*(1+abs(mp.log(lam/eps))))]
        assert max(norm)<4
        rows.append(dict(epsilon=float(eps),derivative=j,I=float(val[0]),J=float(val[1]),normalized=[float(v) for v in norm]))
    return rows


def scales():
    rows=[]
    for q in map(mp.mpf,['1e4','1e8','1e16','1e32']):
        a=1/(4*mp.pi**2*q);alpha=1/(4*mp.pi**2*q**3);eps=q**-2;lam=1/(16*mp.pi*q**3)
        der=product_derivatives(mp.mpf(1),eps,lam,a)
        local=[der[0][j]+der[1][j]/(4*mp.pi) for j in [0,1]] # c=1 example
        delta=mp.log(q)**-2
        rows.append(dict(q=float(q),I= float(product(1,eps,lam,a)[0]),
           canonical_unit_coefficient_curvature=float(product(1,eps,lam,a)[0]),
           proper_unit_coefficient_curvature=float(der[0][1]),local_first=float(local[0]),
           local_second=float(local[1]),leading_a_over_4=float(a/4),
           local_second_over_q_delta=float(abs(local[1])/(q*delta)),
           formal_first_order_scale=float(a*alpha*(1+abs(mp.log(lam/eps))+abs(mp.log(eps))))))
    assert rows[-1]['local_second_over_q_delta']<rows[0]['local_second_over_q_delta']
    return rows


def main():
    start=time.time();data={'symbolic':symbolic()}
    print('Symbolic identities passed',flush=True)
    data['product_integrals'],data['clock_refinement']=product_checks()
    print('Product integrals and independent clock differences passed',flush=True)
    data['first_order'],data['affine_geometry_refinement']=first_order_checks()
    print('Corrected-throat coefficients and implicit-cutoff refinements passed',flush=True)
    data['collars']=collar_checks();data['scales']=scales()
    sha=lambda p:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
    inputs=['notes/proper-mouth-local-stress.md','notes/singlet-static-curvature.md',
      'notes/proper-clock-anomaly.md','notes/casimir-reference-matching.md',
      'notes/collective-boundary-casimir.md','evidence/proper-mouth-local-source-review.json','requirements.txt']
    data.update(status='passed_proper_clock_local_comparison_full_determinant_pending',date='2026-10-11',
      script_sha256=sha(Path(__file__)),input_sha256={p:sha(ROOT/p) for p in inputs},
      versions=dict(python=platform.python_version(),sympy=sp.__version__,mpmath=mp.__version__,decimal_precision=mp.mp.dps),
      elapsed_seconds=time.time()-start,nonlocal_finite_mouth_determinant_controlled=False,
      physical_clock_vertex_derived=False,higher_order_physical_geometry_controlled=False,
      physical_damping_proved=False,novelty_established=False,lean_invocations=0)
    (ROOT/'evidence/proper-mouth-local-stress.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(status=data['status'],elapsed_seconds=data['elapsed_seconds'],
       counts={k:len(data[k]) for k in ['symbolic','product_integrals','clock_refinement','first_order','affine_geometry_refinement','collars','scales']},
       max_product_error=max(x['absolute_integral_error'] for x in data['product_integrals']),
       max_first_order_error=max(x['absolute_coefficient_error'] for x in data['first_order']))))


if __name__=='__main__':main()
