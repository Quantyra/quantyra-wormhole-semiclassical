"""Whole finite-seed trace defect; no full quantum shape or coupled bolt claim."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
STEM='current-seed-full-trace-defect'; PREVIOUS='current-seed-full-hessian-realization'
OWN=['notes/'+STEM+'.md','checks/current_seed_full_trace_defect.py','evidence/'+STEM+'-review.json']
ROWS=[]; FIXTURES=[]

def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name, ok, detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok: raise AssertionError(name)

def algebra():
    r,G,D,eps,T,rho,pp=sp.symbols('r G D epsilon T rho p', nonzero=True)
    F=sp.Function('F')(r); f=sp.Function('f')(r); a=sp.Function('a')(r)
    RR=sp.Function('R4')(r); tau,c=sp.symbols('tau c')
    fp=(-2*G*rho-f)/r; ap=G*(rho+pp)/(r*F)
    # Angular Ricci of the general static areal metric, differentiated first.
    v=sp.symbols('v'); Fv=F+v*f; av=v*a
    RicSv=(1-Fv-r*sp.diff(Fv,r)-r*Fv*sp.diff(av,r))/r**2
    dRicS=sp.diff(RicSv,v).subs(v,0).subs({sp.diff(f,r):fp,sp.diff(a,r):ap})
    q=G*(rho-pp)
    check('angular_Ricci_constraint_identity',sp.simplify(dRicS-q/r**2)==0)
    dRicT=-2*G*rho/r**2-T/2; dRicR=2*G*pp/r**2-T/2
    check('Ricci_trace_constraint',sp.simplify(dRicT+dRicR+2*dRicS+T)==0)
    dJ=sp.expand(2*tau*(-dRicT-dRicR+2*dRicS))
    check('full_tracefree_Ricci_variation',sp.simplify(dJ-(2*tau*T+8*tau*q/r**2))==0)
    dCS=-f/r**2-dRicS-T/6
    dC2=24*c*dCS
    check('full_Weyl_square_variation',sp.simplify(dC2+24*c*((f+q)/r**2+T/6))==0)
    lapv=sp.exp(-av)/r**2*sp.diff(sp.exp(av)*r**2*Fv*sp.diff(RR,r),r)
    dlap=sp.diff(lapv,v).subs(v,0).subs({sp.diff(f,r):fp,sp.diff(a,r):ap})
    expected=f*(sp.diff(RR,r,2)+sp.diff(RR,r)/r)-q*sp.diff(RR,r)/r
    check('nonconstant_curvature_Laplacian_variation',sp.simplify(dlap-expected)==0)
    lapT=sp.symbols('DeltaT')
    raw=D*lapT+eps*(sp.Rational(7,12)*dC2+sp.Rational(11,6)*dJ+sp.Rational(11,36)*RR*T)-D*dlap
    KT=sp.Rational(11,36)*RR-sp.Rational(7,3)*c+sp.Rational(11,3)*tau
    Af=-14*eps*c-D*(r*r*sp.diff(RR,r,2)+r*sp.diff(RR,r))
    AQ=eps*(sp.Rational(44,3)*tau-14*c)+D*r*sp.diff(RR,r)
    check('complete_areal_trace_defect_identity',sp.simplify(raw-(D*lapT+eps*KT*T+Af*f/r**2+AQ*q/r**2))==0)
    tt=sp.Function('T0')(r); xi=sp.Function('xi')(r)
    # Differentiate the pulled-back density, including its Jacobian.
    Rv=r+v*xi
    density=sp.exp(v*a)*Rv**2*(1+v*sp.diff(xi,r))*tt.subs(r,Rv)
    dd=sp.diff(density,v).subs(v,0)/r**2
    check('fixed_disk_density_pullback',sp.simplify(dd-((a+sp.diff(xi,r)+2*xi/r)*tt+xi*sp.diff(tt,r)))==0)
    FIXTURES.append(dict(fixture='exact_full_trace_and_density_identities',full_shape_variation_evaluated=False))

def polynomials_and_budgets():
    z,u,d,beta=sp.symbols('z u d beta')
    k=Q(2000); hm=k**-10; hn=hm/8; pi=Q('3.14160'); pi_lo=Q('3.14159')
    ub=Q(1)/Q('8e9'); db=Q('8.6e-19'); bb=7*hm*hm
    eta=hm/(60*pi_lo); sb=Q('1.001')/(60*pi_lo*k**8); dl=Q('2.946e-9')
    ff=z**-2+u-(4+2*u+beta)*z+(3+u+d*u)*z**2-2*d*u*z**3+d*u*z**4
    rr=4*z**5-6*z**6; rn=-12+d*u*rr
    tn=(3+u+d*u)*z**4-d*u*(5*z**5-sp.Rational(9,2)*z**6)
    cn=sp.expand(-z**3*sp.diff(ff,z)+rn/6)
    P=7*cn**2+sp.Rational(22,3)*tn**2-sp.Rational(11,72)*(rn**2-144)
    Ps=2*(6+d*u)**2+sp.Rational(22,3)*(6+d*u)*u+sp.Rational(22,3)*u*u+22
    Bp=sp.expand(P+Ps*rr/2); A=-10*z**4-20*z**5+30*z**6
    def L(p): return sp.expand(z**4*sp.diff(ff*sp.diff(p,z),z))
    def dt(p): return sp.expand(-z*z*sp.diff(p,z))
    def exactq(v): return Q(int(sp.numer(v)),int(sp.denom(v)))
    def cbound(p):
        poly=sp.Poly(sp.expand(p),u,d,beta)
        return sum(abs(exactq(c))*ub**a*db**b*bb**c0 for (a,b,c0),c in poly.terms())
    def abscoeff(p):
        return sum(cbound(c) for c in sp.Poly(sp.expand(p),z).all_coeffs())
    def derivbound(p):
        return sum(n[0]*cbound(c) for n,c in sp.Poly(sp.expand(p),z).terms())
    ledger={}
    def wbound(name,p):
        p=sp.expand(p); p0=sp.expand(p.subs(beta,0))
        check(name+'_extremal_zero',sp.simplify(p0.subs(z,1))==0)
        terms=sp.Poly(p0,z).terms(); deg=sp.degree(p0,z)
        choices=[sum(abs(n[0]-j)*cbound(c) for n,c in terms) for j in range(int(deg)+1)]
        best=min(choices); jj=choices.index(best)
        # Exact thermal parameter polynomial and conservative w >= h_min.
        therm=abscoeff(p-p0)/hn
        result=best+therm
        ledger[name]={'pivot':jj,'zero_temperature_bound':str(best),
                      'thermal_bound':str(therm),'total_decimal':float(result)}
        return result
    check('product_trace_polynomial_cancellation',sp.simplify(Bp.subs({z:1,beta:0}))==0)
    LR=L(rr)
    tbound=sb*wbound('A',A)+eta*wbound('Bp',Bp)+dl*db*ub*wbound('Lr',LR)
    d1=sb*derivbound(A)+eta*derivbound(Bp)+dl*db*ub*derivbound(LR)
    lbound=sb*wbound('LA',L(A))+eta*wbound('LBp',L(Bp))+dl*db*ub*wbound('LLr',L(LR))
    check('whole_exterior_trace_bound',tbound<Q('9e-28'),float(tbound))
    check('weighted_trace_transport_derivative',d1<Q('7e-27'),float(d1))
    check('whole_exterior_trace_Laplacian',lbound<Q('4e-25'),float(lbound))
    # H^a_a/(lambda^2) is expressed entirely in normalized geometry.
    hr=-2*rn*tn+2*(dt(ff)/2+2*z*ff)*dt(rn)
    ht=hr+2*ff*dt(dt(rn))
    qt=sp.expand(z*dt(ht/z**2)); qr=sp.expand(z*dt(hr/z**2))
    bht=abscoeff(qt); bhr=abscoeff(qr)
    check('R2_energy_derivative_polynomial',bht<Q(200),float(bht))
    check('R2_pressure_derivative_polynomial',bhr<Q(200),float(bhr))
    qquant=(31*Q('1e15')+Q('1e6'))*Q('1.001')/k**8
    qr2=200*8*pi*Q('1e22')*hm
    qref=Q('3.4e-22')+Q('1.179e-8')*Q('3.001')+12*db*ub
    check('complete_current_energy_derivatives',qquant+qr2+qref<Q('1e-7'),float(qquant+qr2+qref))
    srat=Q('3.33e-10'); lapse=Q('3.3e-10')+Q('9.182e-37'); q0=Q('1e-7')
    fbar=lapse+3*q0+80*srat*hm
    hh=lapse+fbar; vv=srat+10*hh
    check('current_moving_bolt_relative_F',fbar<Q('3.004e-7'),float(fbar))
    check('current_H_value',hh<Q('3.01e-7'),float(hh))
    check('current_disk_displacement',vv<Q('3.02e-6'),float(vv))
    innerH=3*q0/10+Q('.7')*lapse+Q('1e-25')
    check('inner_H_linear_zero',innerH<Q('1e-7'))
    innerXi=Q('3.01e-7')+Q('2.21')*(Q('1e-7')+Q('2.8')*Q('3.01e-7'))
    outerXi=261*Q('3.01e-7')
    check('current_disk_displacement_derivative',max(innerXi,outerXi)<Q('8e-5'))
    # h/w <= h_max^(3/5) = 2000^-6.
    fw=7*Q('3.004e-7')+22*srat*k**-6
    qw=2*q0+Q('2.8e-11')*k**-6
    check('weighted_F_correction',fw<Q('2.11e-6'))
    check('weighted_pressure_difference',qw<Q('2.01e-7'))
    density=lapse+Q('8e-5')+2*Q('3.02e-6')
    check('moving_density_contact_coefficient',density<Q('8.7e-5'))
    check('curvature_trace_coefficient',Q(11,36)*Q('12.001')+Q(7,3)*35+Q(11,3)*Q('3.01')<100)
    af=490*eta+316*dl*db*ub; aq=535*eta+56*dl*db*ub
    areal=dl*Q('4e-25')+100*eta*Q('9e-28')+af*Q('2.11e-6')+aq*Q('2.01e-7')
    disk=areal+Q('8.7e-5')*Q('9e-28')+Q('3.02e-6')*Q('7e-27')
    check('complete_disk_trace_defect_bound',disk<Q('1.1e-31'),float(disk))
    response=Q(4,3)*Q('1e18')*Q('1.1e-31')
    check('scalar_only_cold_response',response<Q('1.467e-13'),float(response))
    threshold=Q('1.667e-9')/(Q(4,3)*Q('1e18'))-Q('1.1e-31')
    check('remaining_shape_to_trace_threshold',threshold==Q('1.25014e-27'))
    FIXTURES.append(dict(fixture='exact_rational_whole_exterior_bounds',polynomial_ledger=ledger,
        current_energy_derivative_bound=float(qquant+qr2+qref),trace_defect_bound=str(Q('1.1e-31')),
        scalar_only_cold_bound=str(response),unknown_shape_to_trace_target=str(threshold),
        full_shape_defect_evaluated=False,full_coupled_bolt_error_bounded=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}; error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
        previous_goal_turn='progress',protocol=review['protocol'])
    try:
        check('previous_finite_seed_inverse_run_passed',old['status']=='passed' and old['current_seed_action_form_inverse_proved'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        algebra(); polynomials_and_budgets()
    except Exception: error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,
        input_sha256_after=after,unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,
        finished_utc=datetime.now(timezone.utc).isoformat(),
        full_local_trace_defect_identity_derived=error is None,
        full_disk_trace_defect_bounded=error is None,scalar_only_trace_defect_cold_response_bounded=error is None,
        full_quantum_shape_defect_evaluated=False,full_defect_bolt_pairing_evaluated=False,
        current_seed_admission_inverse_constant_computed=False,actual_physical_correction_computed=False,
        actual_finite_boundary_return_bounded=False,raw_C5_admission_certified=False,
        nonlinear_solution_constructed=False,semiclassical_validity_proved=False,
        physical_stability_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json'); receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,
                    dict(event='finish',status=run['status'],exception=error)]: log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed': raise SystemExit(1)

if __name__=='__main__': main()
