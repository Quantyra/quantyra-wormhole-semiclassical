"""Full cold action-defect assembly, not the finite-domain inverse solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback, math
import sympy as sp
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='current-seed-cold-action-defect';PREVIOUS='current-seed-directional-cold-response'
OWN=['notes/'+STEM+'.md','checks/current_seed_cold_action_defect.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def identities():
    r,ell,G=sp.symbols('r ell G',positive=True);F=sp.Function('F')(r);a=sp.Function('a')(r)
    R=sp.Function('R')(r);O=sp.Function('O')(r);T=sp.Function('T')(r);f=sp.Function('f')(r)
    rho,pp=sp.symbols('rho p')
    # Direct covariant Hessian in the static warped geometry, before varying it.
    htt=(sp.diff(F,r)/2+F*sp.diff(a,r))*sp.diff(R,r)
    hrr=F*sp.diff(R,r,2)+sp.diff(F,r)*sp.diff(R,r)/2
    box=htt+hrr+2*F*sp.diff(R,r)/r
    Hsum=2*R*O-R**2+2*(2*box-htt-hrr)
    target=2*R*O-R**2+2*F*sp.diff(R,r,2)+2*(sp.diff(F,r)+F*sp.diff(a,r)+4*F/r)*sp.diff(R,r)
    check('static_R2_shape_tensor',sp.simplify(Hsum-target)==0)
    Rv=R-ell*T;Ov=O+ell*(2*G*(pp-rho)/r**2-T);Fv=F+ell*f;av=ell*a
    Hv=2*Rv*Ov-Rv**2+2*Fv*sp.diff(Rv,r,2)+2*(sp.diff(Fv,r)+Fv*sp.diff(av,r)+4*Fv/r)*sp.diff(Rv,r)
    expected=-2*O*T+4*G*R*(pp-rho)/r**2-2*F*sp.diff(T,r,2)-2*(sp.diff(F,r)+4*F/r)*sp.diff(T,r)+2*f*sp.diff(R,r,2)+2*(sp.diff(f,r)+F*sp.diff(a,r)+4*f/r)*sp.diff(R,r)
    check('complete_areal_shape_variation',sp.simplify(sp.diff(Hv,ell).subs(ell,0)-expected)==0)
    xi=sp.Function('xi')(r);SL=sp.Function('SL')(r);dSL=sp.Function('dSL')(r);u=sp.symbols('u')
    pulled=sp.exp(2*ell*u)*(SL.subs(r,r+ell*xi)+ell*dSL.subs(r,r+ell*xi))
    check('shape_density_and_pullback',sp.simplify(sp.diff(pulled,ell).subs(ell,0)-(2*u*SL+xi*sp.diff(SL,r)+dSL))==0)
    aa,xxp,nu=sp.symbols('a xi_prime nu');uu=(aa+xxp)/2-nu
    check('trace_density_same_chart',sp.expand(aa+xxp+2*nu-2*uu-4*nu)==0)
    x,y=sp.symbols('x y');I0=x*x/2+2*y*y;Gr=3*x-2*y+x*y+2*x*x
    grad=lambda I:sp.Matrix([sp.diff(I,x),sp.diff(I,y)])
    at0=lambda v:v.subs({x:0,y:0})
    h=-sp.hessian(I0,(x,y)).inv()*at0(grad(Gr))
    defect=at0(grad(I0+Gr))+sp.hessian(I0+Gr,(x,y))*h
    check('exact_reference_action_defect_identity',defect==sp.hessian(Gr,(x,y))*h)
    path={x:ell*h[0]+ell**2,y:ell*h[1]+2*ell**2}
    accel=sp.diff((I0+Gr).subs(path,simultaneous=True),ell,2).subs(ell,0)-(h.T*sp.hessian(I0+Gr,(x,y))*h)[0]
    check('nonlinear_path_acceleration_is_not_Hessian_row',accel!=0 and accel==(at0(grad(Gr)).T*sp.Matrix([2,4]))[0])
    ss,lap,Rg=sp.symbols('s lap Rg')
    conformalR=sp.exp(-2*ell*ss)*(Rg-2*ell*lap)
    check('LLL_exact_density_curvature_contact',sp.simplify(sp.diff(sp.exp(2*ell*ss)*conformalR,ell).subs(ell,0)+2*lap)==0)
    FIXTURES.append(dict(fixture='exact_tensor_density_action_identities',reference_iteration_convergence_claimed=False))

def polynomial_inputs():
    z,u,d,beta,c,e=sp.symbols('z u d beta c e');hm=Q(2000)**-10
    vals=[Q(1)/Q('8e9'),Q('8.6e-19'),7*hm**2,Q('4.91e-10'),Q('3.4e-22')]
    def cb(p):return sum(abs(Q(int(sp.numer(a)),int(sp.denom(a))))*math.prod(v**n for v,n in zip(vals,exps)) for exps,a in sp.Poly(sp.expand(p),u,d,beta,c,e).terms())
    def ac(p):return sum(cb(a) for a in sp.Poly(sp.expand(p),z).all_coeffs())
    dt=lambda p:sp.expand(-z*z*sp.diff(p,z))
    fn=z**-2+u-(4+2*u+beta)*z+(3+u+d*u)*z*z-2*d*u*z**3+d*u*z**4
    wr=4*z**5-6*z**6;rn=-12+d*u*wr;tn=(3+u+d*u)*z**4-d*u*(5*z**5-sp.Rational(9,2)*z**6)
    hr=sp.expand(-2*rn*tn+2*(dt(fn)/2+2*z*fn)*dt(rn));ht=sp.expand(hr+2*fn*dt(dt(rn)))
    J=(24*c*(3+u+d*u)+(1-24*c)*e)*z**4-4*d*u*z**5+3*d*u*z**6
    pdiff=sp.expand(J-c*(hr+ht)/2);L=sp.expand(pdiff/z**4)
    for j,limit in enumerate(['1.16e-22','4e-28','1e-27']):
        val=Q(2)/Q('5.9')*ac(L)
        check('combined_local_shape_jet_'+str(j),val<Q(limit),float(val));L=dt(L)
    vacuum=6/Q('5.9')*11*hm/(360*Q('3.14159'))*Q('1.001')**4
    check('local_reference_vacuum',vacuum<Q('2e-35'),float(vacuum))
    cn=sp.expand(-z**3*sp.diff(fn,z)+rn/6)
    ps=2*(6+d*u)**2+sp.Rational(22,3)*(6+d*u)*u+sp.Rational(22,3)*u*u+22
    bp=sp.expand(7*cn**2+sp.Rational(22,3)*tn**2-sp.Rational(11,72)*(rn**2-144)+ps*wr/2)
    aa=-10*z**4-20*z**5+30*z**6;lr=sp.expand(z**4*sp.diff(fn*sp.diff(wr,z),z))
    weights=[Q('1.001')/(60*Q('3.14159')*Q(2000)**8),hm/(60*Q('3.14159')),Q('2.946e-9')*vals[0]*vals[1]]
    polys=[aa,bp,lr]
    check('trace_degree_and_product_zero',all(sp.degree(p,z)<=12 and sp.simplify(p.subs({z:1,beta:0}))==0 for p in polys))
    first=sum(w*ac(z*sp.diff(p,z)) for w,p in zip(weights,polys))
    thermal=sum(w*cb(p.subs(z,1)) for w,p in zip(weights,polys))/(hm/8)
    check('cold_trace_linear_zero',first+thermal<Q('8e-27'),dict(first=float(first),thermal=float(thermal)))
    for j,limit in [(1,'7e-27'),(2,'1e-25'),(3,'1.3e-24')]:
        val=0
        for w,p in zip(weights,polys):
            for k in range(j):p=dt(p)
            val+=w*ac(p)
        check('trace_ordinary_derivative_'+str(j),val<Q(limit),float(val))
    on=rn/2-2*tn
    check('orbit_Ricci_inputs',ac(on)<13 and ac(dt(on))<30)
    FIXTURES.append(dict(fixture='current_exact_polynomial_inputs',old_separate_R2_reference_bound_used=False))

def numerical_budgets():
    s=Q('.001');t=1+s;du=Q('8.6e-19')/Q('8e9');chi=Q('4.91e-10');hm=Q(2000)**-10;pi=Q('3.14159')
    qr=Q('1.3e-10');rhob=2*Q('1e-12')*7*t;rhoder=qr*t
    check('complete_cold_seed_energy',max(rhob,rhoder)<Q('1.31e-10'))
    S=4/Q('5.9')*t*t*Q('1.31e-10')
    St=2/Q('5.9')*(4*t*Q('1.31e-10')*s+2*qr*t**3)
    check('complete_cold_seed_shape',S<Q('8.9e-11') and St<Q('9e-11'),dict(value=float(S),derivative=float(St)))
    f=Q('1e-8');ft=Q('4e-10');f2=Q('2.7e-10');D=Q('1.1e-9');D1=Q('5.5e-10')
    check('cold_f_value',Q('7.5e-10')*7+Q('3.33e-10')*14<f)
    check('cold_f_first',f*s+2*Q('1.31e-10')<ft)
    check('cold_D_value',ft+Q('2.62e-10')+4*f*s<D)
    check('cold_f_second',ft*s+f*s*s+2*qr+Q('2.62e-10')*s<f2)
    check('cold_D_first',f2+2*qr+Q('2.62e-10')*s+4*ft*s+4*f*s*s<D1)
    T0=Q('8e-27');T1=Q('7e-27');T2=Q('1e-25');T3=Q('1.3e-24')
    bracket=2*(13*T0+7*s*T2+(14+28*s)*T1+f*s*372*du+D*56*du)
    bracketd=2*(30*T0*s+13*T1+14*s*T2+7*s*s*T3+(Q('12.1')+56*s+28*s*s)*T1+(14*s+28*s*s)*T2+ft*s*372*du+f*s*s*2856*du+D1*56*du+D*s*372*du)
    check('R2_trace_curvature_bracket',bracket<Q('4.1e-25') and bracketd<Q('4e-25'),dict(value=float(bracket),derivative=float(bracketd)))
    contacts=2*Q('1e-6')*(Q('1.16e-22')+Q('2e-35'))+Q('1.2e-8')*(Q('4e-28')+8*Q('1e-35'))
    principal=2*chi*Q('12.001')*Q('8.9e-11')
    local=principal+contacts+chi/Q('5.9')*t**4*bracket
    check('complete_R2_reference_shape_value',contacts<Q('2.4e-28') and local<Q('1.06e-18'),float(local))
    # chi derivative / a_s; the final factor 2 uses a_s <= 2 h L_h.
    main_d=4*chi*(Q('12.001')*Q('9e-11')+56*du*Q('8.9e-11')*s)
    contact_d=2*Q('8e-6')*(Q('1.16e-22')+Q('2e-35'))+4*Q('1e-6')*(Q('4e-28')+8*Q('1e-35'))*s+Q('.005')*(Q('4e-28')+8*Q('1e-35'))+2*Q('1.2e-8')*s*(Q('1e-27')+24*Q('1e-35'))
    trace_d=chi/Q('5.9')*(8*t**3*s*bracket+2*t**4*bracketd)
    check('R2_reference_shape_bolt_oscillation',2*(main_d+contact_d+trace_d)<Q('5e-18'),float(2*(main_d+contact_d+trace_d)))
    eta=Q('8.5e-33'); direction=Q('4e-5')
    weyl4=eta/pi*2*Q('1e4')*direction*math.factorial(4)*100**4
    weyl5=eta/pi*2*Q('1e5')*direction*math.factorial(5)*100**5
    mag=eta/pi*160000*direction
    check('Weyl_direct_variation_polynomial',weyl4+mag<Q('6e-24') and weyl5+mag<Q('3e-20'))
    weyl=3*Q('6e-24')+2*Q('1e-6')*Q('3e-15')*Q('.00331')
    gradseed=3*eta/pi*2*Q('1e5')*100*math.factorial(5)*100**5
    check('Weyl_seed_gradient_for_contact',gradseed<Q('2e-13'))
    weylosc=2*Q('1.2')*3*Q('3e-20')+4*(2*Q('8e-6')*Q('3e-15')+2*Q('1e-6')*Q('1.2')*Q('2e-13'))*81*hm
    check('Weyl_fixed_density_value_and_bolt',weyl<Q('1e-22') and weylosc<Q('1e-18'),dict(value=float(weyl),osc=float(weylosc)))
    kl=eta*400/(12*pi)
    check('LLL_fixed_density_bounds',kl<Q('9.1e-32') and 2*kl*Q('.16')<Q('3e-32') and 2*kl*100*Q('1.2')*2<Q('1e-27'))
    ep=hm/(60*pi);dl=Q('2.946e-9');Af=490*ep+316*dl*du;AQ=535*ep+56*dl*du
    et=dl*(7*s*T2+(14+14*s)*T1)+100*ep*T0+Af*f*s+AQ*Q('2.62e-10')
    disk=et+(2*Q('1e-6')*Q('.00331')+Q('4.8e-8')*s)*T0+Q('1.2e-8')*T1
    check('complete_cold_trace_value',disk<Q('5e-34') and t**4/Q('5.9')*disk<Q('8.6e-35'),float(disk))
    lapd=14*s*T2+7*s*s*T3+(Q('12.1')+28*s+14*s*s)*T1+(14*s+14*s*s)*T2
    etdt=dl*lapd+ep*(1000*T0*s+100*T1)+Af*(12*f*s*s+ft*s)+AQ*(12*Q('2.62e-10')*s+2*qr*t)
    diskchi=2*etdt+(2*Q('8e-6')+4*Q('.0048'))*T0*s+(2*Q('1e-6')+Q('4.8e-8'))*2*s*T1+Q('.005')*T1+2*Q('1.2e-8')*s*T2
    check('complete_trace_chi_derivative',diskchi<Q('1e-28'),float(diskchi))
    check('complete_trace_bolt_oscillation',2*(t**4/Q('5.9')*diskchi+8*t**3/Q('5.9')*s*disk)<Q('1e-27'))
    total=Q('1.72e-15')+Q('1.06e-18')+Q('1e-22')+Q('3e-32')
    osc=Q('1.72e-14')+Q('5e-18')+Q('1e-18')+Q('1e-27')
    check('complete_two_row_cold_defect',total<Q('1.73e-15') and osc<Q('1.73e-14'),dict(shape=float(total),osc=float(osc)))
    check('new_return_allocation_not_old_C3_constant',Q('1.73e-15')*Q('578034.6820809248')<=Q('1e-9') and Q('1e-9')+Q('6.67e-10')==Q('1.667e-9'))

def manufactured():
    r,v=sp.symbols('r v');G=sp.Rational(7,10);c=sp.Rational(13,100)
    F=1+r*r-2/r+sp.Rational(3,10)/r**2+sp.Rational(1,500)/r**4
    f=sp.Rational(1,25)/r+sp.Rational(1,50)/r**2
    xi=r/(200*(1+r*r));a=sp.diff(xi,r)-f/F-xi*sp.diff(F,r)/F
    u=(a+sp.diff(xi,r))/2-xi/r
    def curvature(ff,aa):
        R2=-sp.diff(ff,r,2)-3*sp.diff(aa,r)*sp.diff(ff,r)-2*ff*(sp.diff(aa,r,2)+sp.diff(aa,r)**2)
        RR=R2+2*(1-ff)/r**2-4*(sp.diff(ff,r)+ff*sp.diff(aa,r))/r
        OO=R2-2*(sp.diff(ff,r)+ff*sp.diff(aa,r))/r
        HH=2*RR*OO-RR**2+2*ff*sp.diff(RR,r,2)+2*(sp.diff(ff,r)+ff*sp.diff(aa,r)+4*ff/r)*sp.diff(RR,r)
        return RR,OO,HH
    RR,OO,HH=curvature(F+v*f,v*a)
    T=-sp.diff(RR,v).subs(v,0);R0=RR.subs(v,0);O0=OO.subs(v,0)
    rho=-sp.diff(r*f,r)/(2*G);pp=r*F*sp.diff(a,r)/G-rho
    rhs=-2*O0*T+4*G*R0*(pp-rho)/r**2-2*F*sp.diff(T,r,2)-2*(sp.diff(F,r)+4*F/r)*sp.diff(T,r)+2*f*sp.diff(R0,r,2)+2*(sp.diff(f,r)+F*sp.diff(a,r)+4*f/r)*sp.diff(R0,r)
    ref=sp.Rational(1,50)*r*r+sp.Rational(1,100)/r**2-sp.Rational(3,1000)/r**3+sp.Rational(1,250)/r**4
    SL=-16*sp.pi*G*c*r**4*HH+4*G*r*r*ref
    base=SL.subs(v,0);target=2*u*base+xi*sp.diff(base,r)-16*sp.pi*G*c*r**4*rhs
    # Direct full-coordinate pullback includes area and Jacobian, rather than manually adding contacts.
    moved=r+v*xi
    full=sp.exp(v*a.subs(r,moved))*(1+v*sp.diff(xi,r))*(r/moved)**2*SL.subs(r,moved)
    mp.mp.dps=70
    fn=sp.lambdify((v,r),full,'mpmath');tar=sp.lambdify(r,target,'mpmath')
    direct=sp.lambdify(r,sp.diff(SL,v).subs(v,0),'mpmath');simplified=sp.lambdify(r,-16*sp.pi*G*c*r**4*rhs,'mpmath')
    diagnostics=[]
    for rr in [mp.mpf('1.2'),mp.mpf('2'),mp.mpf('4')]:
        expected=tar(rr);errs=[]
        for hh in [mp.mpf('1e-12'),mp.mpf('5e-13')]:
            val=(fn(hh,rr)-fn(-hh,rr))/(2*hh);errs.append(abs(val-expected)/(1+abs(expected)))
        check('manufactured_full_pullback_'+str(rr),abs(direct(rr)-simplified(rr))<mp.mpf('1e-45') and max(errs)<mp.mpf('1e-20') and (errs[1]<mp.mpf('.26')*errs[0] or max(errs)<mp.mpf('1e-55')),list(map(str,errs)))
        diagnostics.append(dict(r=str(rr),normalized_errors=list(map(str,errs))))
    contact=sp.lambdify(r,2*u*base+xi*sp.diff(base,r),'mpmath')
    check('manufactured_contacts_are_nonzero',abs(contact(mp.mpf('1.2')))>mp.mpf('1e-6'))
    FIXTURES.append(dict(fixture='independent_full_curvature_and_coordinate_pullback',diagnostics=diagnostics,actual_quantum_state_sampled=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    try:
        check('previous_directional_test_passed',old['status']=='passed' and old['massive_directional_cold_value_and_bolt_bound'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        identities();polynomial_inputs();numerical_budgets();manufactured()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),complete_cold_fixed_density_action_defect_bounded=error is None,cold_R2_reference_action_contacts_included=error is None,global_full_Hessian_defect_bounded=False,full_local_action_defect_pairing_evaluated=False,full_quantum_shape_defect_evaluated=False,full_defect_bolt_pairing_evaluated=False,actual_physical_correction_computed=False,raw_C5_admission_certified=False,nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
