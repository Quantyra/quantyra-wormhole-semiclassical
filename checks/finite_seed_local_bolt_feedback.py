"""Signed local reference-constraint feedback; not the full Hessian return."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json,traceback
import sympy as sp
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='finite-seed-local-bolt-feedback';PREVIOUS='causal-boundary-feasibility-decision'
OWN=['notes/'+STEM+'.md','checks/finite_seed_local_bolt_feedback.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def algebra():
    r=sp.symbols('r',positive=True);F=sp.Function('F')(r);a=sp.Function('alpha')(r);R4=sp.Function('R4')(r)
    # Independent orbit connection and warped-product Ricci calculation.
    t=sp.symbols('t');xx=[t,r];g=sp.diag(-sp.exp(2*a)*F,1/F);gi=g.inv()
    C=[[[sp.simplify(sum(gi[i,l]*(sp.diff(g[l,k],xx[j])+sp.diff(g[l,j],xx[k])-sp.diff(g[j,k],xx[l])) for l in range(2))/2) for k in range(2)] for j in range(2)] for i in range(2)]
    Ric2=sp.Matrix(2,2,lambda i,j:sp.simplify(sum(sp.diff(C[k][i][j],xx[k])-sp.diff(C[k][i][k],xx[j])+sum(C[k][k][l]*C[l][i][j]-C[k][j][l]*C[l][i][k] for l in range(2)) for k in range(2))))
    Ric4=sp.Matrix(2,2,lambda i,j:Ric2[i,j]+2*C[1][i][j]/r)
    Ric=sp.simplify(gi*Ric4)
    Hess=sp.Matrix(2,2,lambda i,j:sp.diff(R4,xx[i],xx[j])-C[1][i][j]*sp.diff(R4,r))
    Hm=sp.simplify(gi*Hess)
    box=sp.simplify(sp.trace(Hm)+2*F*sp.diff(R4,r)/r)
    H=sp.simplify(2*R4*Ric-sp.eye(2)*R4**2/2+2*(sp.eye(2)*box-Hm))
    ht=2*R4*Ric[0,0]-R4**2/2+2*(F*sp.diff(R4,r,2)+(sp.diff(F,r)/2+2*F/r)*sp.diff(R4,r))
    hd=2*F*((2*R4/r+sp.diff(R4,r))*sp.diff(a,r)-sp.diff(R4,r,2))
    check('actual_static_Htt_identity',sp.simplify(H[0,0]-ht)==0)
    check('actual_static_Hrr_minus_Htt',sp.simplify(H[1,1]-H[0,0]-hd)==0)
    k,Db,R,s,sig,ar,rc,rcp,rcpp,T,Tp,I=sp.symbols('k D_b R s sigma a_R R4_R R4p R4pp T Tp I')
    ric=-Db-k/R;fpp=2*(Db-k/R)
    energy=-R*ric*T/k-2*rc*sig-2*R*s*rcpp-(R*fpp*s/(2*k)+R*ar+4*s)*rcp-R*Tp
    lapse=-R*Tp+T-(2*rc+R*rcp)*ar-I
    result=sp.expand((k/Db)*(energy-lapse)).subs(sig,Db*s/k+ar)
    target=R*T-2*rc*s-2*k*R*s*rcpp/Db-(R+3*k/Db)*s*rcp+k*I/Db
    check('exact_signed_bolt_cancellation',sp.simplify(result-target)==0)
    check('trace_derivative_cancels',sp.diff(result,Tp)==0)
    check('lapse_endpoint_derivative_contact_cancels',sp.simplify(sp.diff(result,ar))==0)
    lam,G,c=sp.symbols('lambda G c_R')
    check('constant_curvature_feedback',sp.simplify((16*sp.pi*G*c*target).subs({rc:-12*lam,rcp:0,rcpp:0,I:0})-(384*sp.pi*G*c*lam*s+16*sp.pi*G*c*R*T))==0)
    w=4/r**5-6/r**6
    check('nonconstant_integrand',sp.simplify(3*sp.diff(w,r)+r*sp.diff(w,r,2)-(60/r**6-144/r**7))==0)
    check('integral_coefficient_budget',Q(60,5)+Q(144,6)==36)
    # Derivative of a scalar row conversion independently checks the off-shell contact.
    z=sp.symbols('z');M=sp.Function('M')(z);E=sp.Function('E')(z);u=sp.symbols('u')
    check('off_shell_row_conversion',sp.expand(M*E+u*sp.diff(M*E,z)-M*(E+u*sp.diff(E,z))-u*sp.diff(M,z)*E)==0)
    FIXTURES.append(dict(fixture='curvature_and_signed_functional_identities',full_action_defect_pairing_evaluated=False))

def bounds():
    z,u,d,beta=sp.symbols('z u d beta');hm=Q(2000)**-10;hn=hm/8
    ub=Q(1)/Q('8e9');db=Q('8.6e-19');bb=7*hm**2
    eta=hm/(60*Q('3.14159'));S=Q('1.001')/(60*Q('3.14159')*Q(2000)**8);dl=Q('2.946e-9')
    ff=z**-2+u-(4+2*u+beta)*z+(3+u+d*u)*z**2-2*d*u*z**3+d*u*z**4
    rr=4*z**5-6*z**6;rn=-12+d*u*rr
    tn=(3+u+d*u)*z**4-d*u*(5*z**5-sp.Rational(9,2)*z**6)
    cn=sp.expand(-z**3*sp.diff(ff,z)+rn/6)
    Ps=2*(6+d*u)**2+sp.Rational(22,3)*(6+d*u)*u+sp.Rational(22,3)*u*u+22
    Bp=sp.expand(7*cn**2+sp.Rational(22,3)*tn**2-sp.Rational(11,72)*(rn**2-144)+Ps*rr/2)
    A=-10*z**4-20*z**5+30*z**6;Lr=sp.expand(z**4*sp.diff(ff*sp.diff(rr,z),z))
    def cb(p):
        return sum(abs(Q(int(sp.numer(c)),int(sp.denom(c))))*ub**i*db**j*bb**k for (i,j,k),c in sp.Poly(sp.expand(p),u,d,beta).terms())
    def ac(p):return sum(cb(c) for c in sp.Poly(sp.expand(p),z).all_coeffs())
    derivative=S*ac(z*sp.diff(A,z))+eta*ac(z*sp.diff(Bp,z))+dl*db*ub*ac(z*sp.diff(Lr,z))
    for name,p in [('A',A),('Bp',Bp),('Lr',Lr)]:
        check('exact_product_zero_'+name,sp.simplify(p.subs({z:1,beta:0}))==0)
    thermal=(S*cb(A.subs(z,1))+eta*cb(Bp.subs(z,1))+dl*db*ub*cb(Lr.subs(z,1)))/hn
    check('ordinary_log_radial_trace_derivative',derivative<Q('7e-27'),float(derivative))
    check('thermal_bolt_offset_over_gap',thermal<Q('1e-60'),float(thermal))
    check('actual_bolt_trace_over_gap',derivative+thermal<Q('8e-27'))
    check('nonconstant_s_coefficient',20+56*(1+3*Q('1.001')*hm)+744*Q('1.001')*hm<78)
    check('physical_curvature_coefficient',16*Q('3.14160')*Q('1e22')*hm<Q('4.91e-10'))
    leading=Q('1.179e-8')*Q('3.33e-10')
    remainder=Q('4.91e-10')*(Q('1.001')*Q('8e-27')+db*ub*(78*Q('3.33e-10')+36*Q('1.001')*(Q('3.3e-10')+Q('9.182e-37'))))
    check('local_signed_reference_bolt_bound',leading+remainder<Q('3.927e-18'),float(leading+remainder))
    check('comparison_margin_not_full_error',Q('1.667e-9')/Q('3.927e-18')>Q('4e8'))
    mp.mp.dps=65;rr0=mp.mpf('1.2');kk=mp.mpf('.13');DD=mp.mpf('.9');ss=mp.mpf('.04')
    aa=lambda x:mp.mpf('.03')/x**2
    curv=lambda x:-mp.mpf('3.6')+mp.mpf('.002')*(4/x**5-6/x**6)
    tr=lambda x:mp.mpf('.007')*(x**-4-x**-5)
    sigma=DD*ss/kk+aa(rr0)
    lapse=mp.quad(lambda x:x*mp.diff(tr,x,2)+(2*curv(x)+x*mp.diff(curv,x))*mp.diff(aa,x),[rr0,2,5,mp.inf])
    ric=-DD-kk/rr0;fpp=2*(DD-kk/rr0)
    energy=-rr0*ric*tr(rr0)/kk-2*curv(rr0)*sigma-2*rr0*ss*mp.diff(curv,rr0,2)-(rr0*fpp*ss/(2*kk)+rr0*aa(rr0)+4*ss)*mp.diff(curv,rr0)-rr0*mp.diff(tr,rr0)
    ii=mp.quad(lambda x:(3*mp.diff(curv,x)+x*mp.diff(curv,x,2))*aa(x),[rr0,2,5,mp.inf])
    closed=rr0*tr(rr0)-2*curv(rr0)*ss-2*kk*rr0*ss*mp.diff(curv,rr0,2)/DD-(rr0+3*kk/DD)*ss*mp.diff(curv,rr0)+kk*ii/DD
    check('independent_integrated_functional',abs(kk*(energy-lapse)/DD-closed)<mp.mpf('1e-55'))
    FIXTURES.append(dict(fixture='exact_benchmark_budget_and_independent_integral',leading=str(leading),remainder=str(remainder),trace_derivative=float(derivative),thermal_offset=float(thermal),manufactured_quadrature_error=mp.nstr(abs(kk*(energy-lapse)/DD-closed),10),actual_quantum_column_evaluated=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    try:
        check('previous_boundary_audit_passed',old['status']=='passed' and old['kinematic_boundary_implication_rejected'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        algebra();bounds()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),local_R2_reference_bolt_feedback_bounded=error is None,full_local_action_defect_pairing_evaluated=False,full_quantum_shape_defect_evaluated=False,full_defect_bolt_pairing_evaluated=False,actual_physical_correction_computed=False,raw_C5_admission_certified=False,nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
