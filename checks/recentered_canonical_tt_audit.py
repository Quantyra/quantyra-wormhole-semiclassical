"""Supporting diagnostics for the recentered canonical/TT audit; no quantum solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='recentered-canonical-tt-audit'; PREVIOUS='recentered-r2-scalar-boundary'
OWN=['notes/'+STEM+'.md','checks/recentered_canonical_tt_audit.py','evidence/'+STEM+'-review.json']
ROWS=[]; FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def exact_algebra():
    R,R0,La,G,c=S.symbols('R R0 La G c'); F=1+32*S.pi*G*c*R0
    lhs=(R-2*La)/(16*S.pi*G)+c*R**2
    rhs=F*(R-R0/2)/(16*S.pi*G)+c*(R-R0)**2+(R0-4*La)/(32*S.pi*G)
    check('exact_local_action_split',S.simplify(lhs-rhs)==0)
    e,dR=S.symbols('e dR')
    check('Euclidean_square_Hessian',S.diff(-c*(e*dR)**2,e,2)==-2*c*dR**2)
    b,xi,t=S.symbols('b xi t'); k=S.diag(-b-2*xi,b,xi,xi); J=S.diag(-t,-t,t,t)
    check('TT_scalar_contraction',S.trace(k)==0 and S.expand(S.trace(J*k))==4*t*xi)
    al,qq=S.symbols('grad_alpha grad_logF')
    check('Einstein_frame_kinetic_identity',S.expand(6*(al-qq/2)**2+6*qq*(al-qq/2)-6*al**2+S.Rational(3,2)*qq**2)==0)
    r,M,Z,lam=S.symbols('r M Z lam',positive=True); f=1-2*M/r+Z/r**2+lam*r*r
    check('RN_scalar_curvature',S.simplify(-S.diff(f,r,2)-4*S.diff(f,r)/r+2*(1-f)/r**2+12*lam)==0)
    check('York_potential_lower_part',S.simplify(S.diff(f,r,2)/4+S.diff(f,r)/(2*r)-3*lam/2-Z/(2*r**4))==0)
    W=3*f-r*S.diff(f,r)/2+r*r*S.diff(f,r,2)/2-r**3*S.diff(f,r,3)/4
    check('completed_TT_potential',S.simplify(W-(3-12*M/r+13*Z/r**2+3*lam*r*r))==0)
    P=r**4*S.diff(f,r)-2*r**3*f
    check('Hardy_cross_derivative',S.simplify(S.diff(P,r)-(-6*r*r+12*M*r-4*Z))==0)
    x0,x1,x2,b0=S.symbols('x0 x1 x2 b0')
    bp=-(S.diff(r*r*f,r)*b0+(r*r*S.diff(f,r)-2*r*f)*x0)/(r*r*f)
    df=-f*b0+f*(x0+r*x1)-r*S.diff(f,r)*x0/2
    def total(q):return S.diff(q,r)+S.diff(q,b0)*bp+S.diff(q,x0)*x1+S.diff(q,x1)*x2
    L=3*r*f+r*r*S.diff(f,r)/2; K=-f+r*S.diff(f,r)-r*r*S.diff(f,r,2)/2
    check('original_areal_TT_variation',S.simplify(total(r*df)-(r*r*f*x2+L*x1+K*x0+f*b0))==0)
    MM=-2*r*f+r*r*S.diff(f,r)/2+r**3*S.diff(f,r,2)/2
    JJ=-3*r*r*f*x0*x1-r**3*f*x1*x1/2-r*f*b0*x0+MM*x0*x0/2
    VF=6*f-7*r*S.diff(f,r)/2+r*r*S.diff(f,r,2)/2-r**3*S.diff(f,r,3)/4
    actual=-(3*x0+r*x1)*total(r*df)
    rebuilt=S.Rational(3,2)*r*r*f*x1*x1+VF*x0*x0-4*f*b0*x0+total(JJ)
    check('exact_classical_TT_integration_identity',S.simplify(actual-rebuilt)==0)
    FIXTURES.append(dict(fixture='local_action_and_TT_identities',square_Hessian='-32 c_R (Z/r^4)^2 xi^2 before angular/time integration',quantum_remainder_evaluated=False))

def boundary_algebra():
    e,t=S.symbols('e t'); a0,ap0,s0,sp0,sg,b2=S.symbols('a ap s sp sigma b2')
    u,u1,u2,j,j1,j2,h,h1=S.symbols('u u1 u2 j j1 j2 h h1')
    a=S.exp(-t*h)*(a0+t*j1); ss=S.exp(-t*h)*(s0+e*u1)
    ap=S.exp(-2*t*h)*(ap0+t*j2-t*h1*(a0+t*j1))
    sp=S.exp(-2*t*h)*(sp0+e*u2-t*h1*(s0+e*u1))
    sig=sg+e*u
    expr=S.exp(t*j)*(sig*(12*a*ap+S.Rational(44,7)*a**3+10*a*b2)-6*ap*ss-S.Rational(24,7)*a*a*ss-6*a*sp-S.Rational(96,7)*a*ss*ss-18*ss*sp-S.Rational(85,21)*ss**3)
    mixed=S.expand(S.diff(expr,e,t).subs({e:0,t:0}))
    Y=S.symbols('Y',positive=True);u4,v4,j3,h3=S.symbols('u4 v4 j3 h3')
    U=Y**4*(u4*S.log(Y)+v4); jj=j3*Y**3;hh=h3*Y**3
    lead=S.expand(mixed.subs({a0:0,ap0:0,s0:1/Y,sp0:1/Y**2,sg:-S.log(Y),u:U,u1:-S.diff(U,Y),u2:S.diff(U,Y,2),j:jj,j1:-S.diff(jj,Y),j2:S.diff(jj,Y,2),h:hh,h1:-S.diff(hh,Y)}))
    L=S.symbols('logY'); poly=S.Poly(lead.subs(S.log(Y),L),Y)
    expected=((5028*h3-452*j3)*(u4*L+v4)+(2769*h3-113*j3)*u4)/7
    check('full_mixed_boundary_leading_coefficient',S.simplify(poly.coeff_monomial(Y**4)-expected)==0)
    check('mixed_boundary_no_lower_weight',all(power[0]>=4 for power,coeff in poly.terms()))
    check('radial_metric_mixed_term_is_nonzero',expected.subs(j3,0)!=0)
    FIXTURES.append(dict(fixture='unexpanded_finite_cut_boundary_variation',leading=str(expected),radial_metric_varied=True,full_physical_domain_proved=False))
    args=(a0,ap0,s0,sp0,sg,b2,u,u1,u2,j,j1,j2,h,h1)
    return S.lambdify(args,mixed,'mpmath')

def rational_bounds():
    pl=Q('3.14159');ph=Q('3.14160');k=2000;g0=Q(1,k**10);gl=(1-Q('1e-34'))*g0;c=Q('1e22')
    flo=1-384*ph*c*g0;beta=Q('1.223e-10');Cmin=Q(49119000000)
    check('recentered_vacuum_root_interval',gl+11*gl*gl/(360*pl)<g0 and g0+11*g0*g0/(360*ph)>g0)
    check('positive_effective_Newton_factor',flo>1-Q('1.179e-8'))
    check('C_RN_lower_guard',1+12*pl*gl*k**16>Cmin**2)
    check('C_RN_upper_guard',1+12*ph*g0*k**16/flo<Q('1e11')**2)
    check('inverse_radius_guard',6/(Cmin-1)<beta)
    # C_new/C_class<1.01, B_new/B_old<2 follows with much room.
    Cclassmin=6*k**3
    check('new_to_classical_C_ratio',1/flo< Q('1.01')**2)
    check('new_to_old_B_ratio',Q('1.01')*Cclassmin/(Cclassmin-1)/(1-Q('1.4e-22'))<4)
    check('prescribed_temperature_gamma_transfer',2*(8*k**3)/Cmin<3)
    hh=Q('6e-30');xmax=Q('1e11')/6
    delta=hh**2*(1+xmax*((1+hh)**2+2*(1+hh)+3))
    check('thermal_delta_below_one',delta<1)
    check('classical_x_above_1000',(Cmin-1)/6>1000)
    eta=256*ph*c*g0*(3+beta)**2/flo
    check('actual_local_TT_relative_budget',eta<Q('7.07e-8'))
    check('quantum_target_positive_margin',1-Q('7.07e-8')-Q(1,4)==Q('.7499999293'))
    rho=Q('1e-8');Y0=Q('5e-9');Z0=Q(1,4);Z2=Q('5e7')
    check('full_radius_inequality',Y0+Z0*rho+Z2*rho*rho/2<=rho)
    check('full_contraction_inequality',Z0+Z2*rho==Q(3,4))
    check('isolated_scalar_embedding_target',Q('2.74725e20')*Q('1.82e-29')<Q('5e-9'))
    check('Schur_scalar_embedding_target',Q('1.03021e20')*Q(8,3)*Q('1.82e-29')<Q('5e-9'))
    check('mixed_schur_inverse_allocation',Q(4,3)*Q(3,8)==Q(1,2) and Q(4,3)/(1-Q(1,2))==Q(8,3))
    FIXTURES.append(dict(fixture='rational_physical_budgets',eta_R_upper=str(eta),F_lambda_lower=str(flo),thermal_delta_upper=str(delta),eta_Q='unknown; selected target <=1/4',nonlinear_ball_certified=False))

def boundary_numerical(mixed):
    records=[]
    for prec in [60,90]:
        mp.mp.dps=prec
        for Ys in ['.03','.07']:
            Y=mp.mpf(Ys); b2=mp.mpf('.4')
            lnN=lambda y:mp.mpf('.3')*y*y-mp.mpf('.2')*y**3
            sig=lambda y:-mp.log(y)
            uu=lambda y:y**4*(mp.mpf('.7')*mp.log(y)+mp.mpf('.2'))
            jj=lambda y:mp.mpf('.6')*y**3+mp.mpf('.1')*y**4*mp.log(y)
            hh=lambda y:-mp.mpf('.4')*y**3+mp.mpf('.2')*y**4*mp.log(y)
            out=lambda fun,n=1:(-1)**n*mp.diff(fun,Y,n)
            vals=[out(lnN),out(lnN,2),out(sig),out(sig,2),sig(Y),b2,uu(Y),out(uu),out(uu,2),jj(Y),out(jj),out(jj,2),hh(Y),out(hh)]
            N=mp.exp(lnN(Y)); expected=N*mixed(*vals)
            def direct(e,t):
                logN=lambda y:lnN(y)+t*jj(y)
                sigma=lambda y:sig(y)+e*uu(y)
                deriv=lambda fun,y:-mp.exp(-t*hh(y))*mp.diff(fun,y)
                aa=deriv(logN,Y);ss=deriv(sigma,Y)
                ap=deriv(lambda y:deriv(logN,y),Y);sp=deriv(lambda y:deriv(sigma,y),Y)
                gg=sigma(Y)*(12*aa*ap+mp.mpf(44)/7*aa**3+10*aa*b2)-6*ap*ss-mp.mpf(24)/7*aa*aa*ss-6*aa*sp-mp.mpf(96)/7*aa*ss*ss-18*ss*sp-mp.mpf(85)/21*ss**3
                return mp.exp(logN(Y))*gg
            errors=[]
            for step in ['1e-5','5e-6']:
                d=mp.mpf(step); fd=(direct(d,d)-direct(d,-d)-direct(-d,d)+direct(-d,-d))/(4*d*d)
                errors.append(abs((fd-expected)/expected))
            ratio=errors[1]/errors[0]
            check('boundary_fd_error_'+Ys+'_'+str(prec),errors[1]<mp.mpf('1e-8'))
            check('boundary_fd_refinement_'+Ys+'_'+str(prec),mp.mpf('.24')<ratio<mp.mpf('.26'))
            records.append(dict(digits=prec,Y=Ys,mixed=mp.nstr(expected,50),errors=[mp.nstr(x,30) for x in errors],refinement=mp.nstr(ratio,20)))
    FIXTURES.append(dict(fixture='independent_finite_cut_mixed_differences',records=records,full_domain_proved=False))

def radial_numerical():
    records=[]
    rr,aa,ee=S.symbols('r a e'); poly=S.Poly(S.expand((rr-aa)**3*(ee-rr)**3),rr)
    for prec in [60,90]:
        mp.mp.dps=prec
        for xs,hs in [('1000','.001'),('10000','.002')]:
            x=mp.mpf(xs);rp=1+mp.mpf(hs);Z=1+3*x;M=(rp+Z/rp+x*rp**3)/2
            start=rp+mp.mpf('.2');end=start+1
            coeff=[mp.mpf(str(z)) for z in S.Poly(poly.as_expr().subs({aa:S.Rational(str(start)),ee:S.Rational(str(end))}),rr).all_coeffs()][::-1]
            xi=lambda r:(r-start)**3*(end-r)**3
            f=lambda r:1-2*M/r+Z/r**2+x*r*r
            def primitive(r):
                return sum(v*(6*M*(r**(i+1)-start**(i+1))/(i+1)-2*(r**(i+2)-start**(i+2))/(i+2)-4*Z*(mp.log(r/start) if i==0 else (r**i-start**i)/i)) for i,v in enumerate(coeff))
            b=lambda r:-primitive(r)/(r*r*f(r))
            df=lambda r:-f(r)*b(r)+f(r)*(xi(r)+r*mp.diff(xi,r))-r*mp.diff(f,r)*xi(r)/2
            original=mp.quad(lambda r:-(3*xi(r)+r*mp.diff(xi,r))*mp.diff(lambda y:y*df(y),r),[start,end])
            square=lambda r:mp.mpf('1.5')*r*r*f(r)*(mp.diff(xi,r)+2*xi(r)/r)**2
            W=lambda r:3-12*M/r+13*Z/r**2+3*x*r*r
            rebuilt=mp.quad(lambda r:square(r)+W(r)*xi(r)**2-4*f(r)*b(r)*xi(r),[start,end])
            energy=mp.quad(lambda r:square(r)+x*r*r*xi(r)**2,[start,end])
            error=abs(original/rebuilt-1)
            check('independent_radial_identity_'+xs+'_'+str(prec),error<mp.mpf('1e-40'))
            check('classical_radial_energy_lower_bound_'+xs+'_'+str(prec),original>energy)
            check('nonzero_retained_TT_tail_'+xs+'_'+str(prec),abs(primitive(end))>0)
            records.append(dict(digits=prec,x=xs,h=hs,identity_error=mp.nstr(error,20),energy_ratio=mp.nstr(original/energy,50),tail_numerator=mp.nstr(primitive(end),50)))
    for xs in ['1000','10000']:
        pair=[r for r in records if r['x']==xs]
        check('radial_precision_'+xs,abs(mp.mpf(pair[0]['energy_ratio'])/mp.mpf(pair[1]['energy_ratio'])-1)<mp.mpf('1e-40'))
    FIXTURES.append(dict(fixture='independent_RN_TT_radial_quadrature',records=records,actual_quantum_state_evaluated=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_recentered_scalar_result',old['status']=='passed' and old['recentered_scalar_source_admissible'])
        check('all_previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        exact_algebra();mixed=boundary_algebra();rational_bounds();boundary_numerical(mixed);radial_numerical()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,
        unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        local_TT_relative_loss_bound='7.07e-8' if error is None else None,mixed_AdS_boundary_limit_zero=error is None,
        weighted_log_regularization_already_available=True,old_quantum_constants_transferred=False,
        quantum_eta_Q_known=False,physical_boundary_completed=False,physical_coupled_inverse_proved=False,
        coupled_correction_norm_proved=False,nonlinear_solution_constructed=False,semiclassical_validity_proved=False,
        physical_stability_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
