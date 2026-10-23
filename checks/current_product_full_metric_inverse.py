"""Frozen full-product and cold-envelope audit; no finite-seed solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback
import sympy as sp
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='current-product-full-metric-inverse';PREVIOUS='physical-metric-block-and-preconditioner'
OWN=['notes/'+STEM+'.md','checks/current_product_full_metric_inverse.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def algebra():
    S,T,Ls,Lt,Gs,e=sp.symbols('S T Ls Lt Gs e')
    X,Y,k2,G,cR,l0,q,gR=sp.symbols('X Y k2 G cR l0 q gR',positive=True)
    Rh=-2*X;R=2*(Y-X);dc=1/(5760*sp.pi**2);cs=cR-dc
    raw=sp.exp(-2*e*T)*(Rh-2*e*Lt-6*e*Ls-6*e*e*Gs)**2+4*Y*(Rh-2*e*Lt-6*e*Ls-6*e*e*Gs)+4*Y**2*sp.exp(2*e*T)
    quadratic=sp.diff(raw,e,2).subs(e,0)/2
    quadratic=sp.expand(quadratic.subs({Ls:-k2*S,Lt:-k2*T,Gs:k2*S*S}))
    HR=sp.hessian(quadratic,(S,T))
    expected=sp.Matrix([[72*k2*k2-24*R*k2,24*k2*(k2+2*X)],[24*k2*(k2+2*X),8*k2*k2+32*X*k2+16*X*X+16*Y*Y]])
    check('R2_Hessian_from_reduced_density',all(sp.simplify(x)==0 for x in HR-expected))
    tau=(6*X*X+10*X*Y+(6-30*q*q)*Y*Y)/(1440*sp.pi**2)
    gs=11*Y*k2/(720*sp.pi**2)
    gt=(X*X+(5*q*q-1)*Y*Y)/(120*sp.pi**2)+(6*X+5*Y)*k2/(720*sp.pi**2)
    ehss=-3*k2/(4*sp.pi*G)-(R+24*l0)/(4*sp.pi*G)
    ehst=-(k2+2*Y+12*l0)/(4*sp.pi*G)
    ell,p=sp.symbols('ell p')
    beta=64*sp.pi*G*cs*X;cc=6+6*beta*(1-Y/X)-11*G*Y/(90*sp.pi)
    dd=2*(R+24*l0)/X;bb=1-G*(6*X+5*Y)/(180*sp.pi)
    gamma=X/(8*sp.pi*G)
    check('full_conformal_matrix_entry',sp.simplify(((ehss+gs-cs*HR[0,0])/gamma).subs(k2,X*ell)+9*beta*ell**2+cc*ell+dd)==0)
    mixed=(ehst+gt-cs*HR[0,1])/gamma
    check('full_mixed_entry_with_stationary_trace',sp.simplify((mixed+2*(k2/X+2)*(bb+sp.Rational(3,2)*beta*k2/X)).subs(l0,(X-Y-4*sp.pi*G*tau)/6))==0)
    pA,qr=sp.symbols('pA qr')
    eh_tt=4*pA-16*cR*(X*X-Y*Y)
    gd_tt=gamma*qr+8*dc*(k2+2*X)**2-4*pA
    check('full_shape_entry_with_stationary_orbit',sp.simplify((eh_tt-cR*HR[1,1]+gd_tt)/gamma-qr+beta*(k2/X+2)**2)==0)
    be,c,d,b,H,alpha=sp.symbols('beta c d b H alpha')
    a=9*be*ell**2+c*ell+d;j=b+sp.Rational(3,2)*be*ell
    B=4*b*b-be*d+be*ell*(12*b-c);D=B+a*(H+alpha/p)
    M=sp.Matrix([[-a,-2*p*j],[-2*p*j,p*alpha+p*p*H-be*p*p]])
    check('determinant_quartic_cancellation',sp.expand(4*j*j-be*a-B)==0)
    check('complete_product_determinant',sp.simplify(M.det()+p*p*D)==0)
    inv=sp.Matrix([[be-H-alpha/p,-2*j/p],[-2*j/p,a/p**2]])/D
    check('full_matrix_inverse_identity',all(sp.simplify(x)==0 for x in M*inv-sp.eye(2)))
    x,L0=sp.symbols('x L0',positive=True)
    J=(x**3*(L0-sp.log(x)+sp.Rational(1,3))-(L0+sp.Rational(1,3)))/3
    check('cold_forcing_integral',sp.simplify(sp.diff(J,x)-x*x*(L0-sp.log(x)))==0 and J.subs(x,1)==0)
    FIXTURES.append(dict(fixture='exact_current_action_matrix_and_cold_Green_integral',finite_seed_stationarity_assumed=False))

def rational():
    be=Q('1.18e-8');aq=Q('1.7e-8');plo=Q('2.25')
    check('B_constant_reserve',4*Q('.999')**2-be*Q('4.01')>Q('3.99'))
    check('B_linear_reserve',12*Q('.999')-Q('6.01')>Q('5.97'))
    check('a_over_p',Q('6.01')+Q('4.01')/plo<Q('7.8'))
    check('D_constant_reserve',Q('3.99')-aq*Q('7.8')>Q('3.98'))
    check('D_linear_reserve',Q('5.97')-9*aq>Q('5.96'))
    check('a_lower',Q('5.99')/4+Q('3.99')==Q('5.4875'))
    check('inverse_SS',1/Q('5.4875')+be/Q('3.98')+aq/(plo*Q('3.98'))<Q('.183'))
    check('inverse_ST_intercept',2*Q('1.001')/Q('3.98')<Q('.504'))
    check('inverse_ST_slope',3/Q('5.96')<Q('.504'))
    check('inverse_ST',Q('.504')/plo==Q('.224'))
    check('inverse_TT_zero_beta',Q('6.01')**2/(4*Q('3.98')*(2*Q('6.01')-Q('4.01')))<Q('.284'))
    check('inverse_TT_large_beta',9/(8*Q('5.96'))<Q('.189')<Q('.284'))
    check('inverse_operator_bound',max(Q('.183')+Q('.224'),Q('.224')+Q('.284'))==Q('.508')<Q('.51'))
    rad=Q('4e-8')
    aa=Q('8.03')+36*be+(Q('6.01')+36*be)*rad+9*be*rad*rad
    check('cold_disk_a',aa<Q('8.04'))
    Blow=4*Q('.999')**2-be*Q('4.01')-be*(2+rad)*Q('6.022')
    check('cold_disk_B',Blow>Q('3.99'))
    HH=Q(2,3)+rad/(1-rad/1404)*Q('6e-5')
    check('cold_disk_H',HH<Q('.667'))
    perturb=Q('8.04')*aq+rad*Q('8.04')*Q('.667')*Q('8.5e-33')
    check('cold_root_Rouche',perturb<rad*Q('3.99'))
    check('cold_exponent_root_guard',9-4*rad>Q('2.9999999')**2 and 2*rad/(3+Q('2.9999999'))<Q('1.34e-8'))
    lower=Q('.72')/6*(72*(69-1)-Q(69**2-1,2))
    check('envelope_rigorous_gain_lower',lower==Q('301.92'))
    check('envelope_response_lower',lower*Q('1.36e-9')==Q('4.106112e-7'))
    check('envelope_rejects_unit_target',lower>200*Q('1.4705882352'))
    check('source_budget_unchanged',Q('1.632e-11')*Q('61.2745098039')+Q('1.36e-9')*Q('1.4705882352')+Q('2e-9')<=Q('5e-9'))

def intervals():
    iv=mp.iv;iv.dps=70;k=iv.mpf(2000);g0=k**-10;cR=iv.mpf('1e22')
    g=2*g0/(1+iv.sqrt(1+4*iv.mpf(11)/(360*iv.pi)*g0))
    Fl=1-384*iv.pi*cR*g;Cr=iv.sqrt(1+12*iv.pi*g*k**16/Fl);xr=(Cr-1)/6
    Yr=g/xr;Xr=Cr*Yr
    X=Xr/(1+iv.mpf(['1.56e-29','1.71e-29']))
    Y=Yr/(1+iv.mpf(['-5.5e-23','5.5e-23']))
    be=64*iv.pi*(cR-1/(5760*iv.pi**2))*X
    c=6+6*be*(1-Y/X)-11*Y/(90*iv.pi)
    d=2*(2*(Y-X)+24*g0)/X;b=1-(6*X+5*Y)/(180*iv.pi)
    check('exact_product_beta_guard',be>0 and be<iv.mpf('1.18e-8'),str(be))
    check('exact_product_c_guard',c>iv.mpf('5.99') and c<iv.mpf('6.01'),str(c))
    check('exact_product_d_guard',d>iv.mpf('3.99') and d<iv.mpf('4.01'),str(d))
    check('exact_product_b_guard',b>iv.mpf('.999') and b<iv.mpf('1.001'),str(b))
    hmax=g0;neck=iv.mpf('.001')/hmax
    Ll=iv.log(neck+iv.sqrt(neck*neck-1));L0=1+iv.log(iv.mpf('.01')/hmax)
    check('envelope_neck_lower',Ll>69,str(Ll))
    check('envelope_log_lower',L0>72,str(L0))
    cosh_one=(iv.exp(iv.mpf(1))+iv.exp(iv.mpf(-1)))/2
    check('envelope_hyperbolic_lower',1-1/cosh_one**3>iv.mpf('.72'))
    check('cold_exponent_factor',iv.exp(iv.mpf(84)*iv.mpf('1.34e-8'))<iv.mpf('1.00000113'))
    FIXTURES.append(dict(fixture='interval_product_and_cold_guards',digits=70,full_exact_product_enclosure_used=True))

def numerical():
    records=[]
    for digits in [70,110]:
        mp.mp.dps=digits;be=mp.mpf('1.178e-8');c=mp.mpf('6.00000007');d=mp.mpf(4);b=mp.mpf(1);alpha=mp.mpf('-1.6e-8')
        worst=mp.mpf(0);maxerr=mp.mpf(0)
        for ell in map(mp.mpf,['.25','2','1e8','1e16']):
            for H in map(mp.mpf,['0','1e-30','1e-4','1e8']):
                p=ell+2;a=9*be*ell**2+c*ell+d;j=b+mp.mpf('1.5')*be*ell
                B=4*b*b-be*d+be*ell*(12*b-c);D=B+a*(H+alpha/p)
                M=mp.matrix([[-a,-2*p*j],[-2*p*j,p*alpha+p*p*H-be*p*p]])
                closed=mp.matrix([[be-H-alpha/p,-2*j/p],[-2*j/p,a/p**2]])/D
                direct=mp.inverse(M)
                err=max(abs(direct[i,j]-closed[i,j]) for i in range(2) for j in range(2))
                maxerr=max(maxerr,err)
                worst=max(worst,max(sum(abs(direct[i,j]) for j in range(2)) for i in range(2)))
        check('independent_matrix_inverse_'+str(digits),maxerr<mp.mpf('1e-45'),mp.nstr(maxerr,20))
        check('matrix_control_with_large_local_terms_'+str(digits),worst<mp.mpf('.51'))
        gains=[]
        for h in [mp.mpf(1)/(8*2000**10),mp.mpf(1)/2000**10]:
            L=mp.acosh(mp.mpf('.001')/h);L0=1+mp.log(mp.mpf('.01')/h)
            def integrand(u):
                if not u:return mp.mpf(0)
                x=mp.cosh(u)
                J=(x**3*(L0-mp.log(x)+mp.mpf(1)/3)-(L0+mp.mpf(1)/3))/3
                return J/(x*x*mp.sinh(u))/2
            gain=mp.quad(integrand,[0,1,10,L])
            check('model_gain_lower_'+str(digits)+'_'+str(len(gains)),gain>mp.mpf('301.92'))
            gains.append(dict(h=mp.nstr(h,40),gain=mp.nstr(gain,45),radius_at_retained_envelope=mp.nstr(gain*mp.mpf('1.36e-9'),45),source_for_two_e_minus_nine=mp.nstr(mp.mpf('2e-9')/gain,45)))
        records.append(dict(digits=digits,matrix_worst_row_sum=mp.nstr(worst,45),gains=gains))
    for i in range(2):check('cross_precision_cold_gain_'+str(i),abs(mp.mpf(records[0]['gains'][i]['gain'])/mp.mpf(records[1]['gains'][i]['gain'])-1)<mp.mpf('1e-40'))
    FIXTURES.append(dict(fixture='independent_matrix_inverse_and_classical_cold_envelope',records=records,
        quadrature_is_interval_certificate=False,actual_quantum_cold_response_computed=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    try:
        check('previous_metric_audit_passed',old['status']=='passed' and old['exact_physical_bulk_rows_derived'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        algebra();rational();intervals();numerical()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,
        unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        current_action_product_L2_inverse_proved=error is None,product_L2_inverse_bound='.51' if error is None else None,
        current_action_product_cold_root_enclosed=error is None,classical_envelope_gain_lower='301.92' if error is None else None,
        actual_quantum_cold_response_computed=False,finite_seed_signed_radius_functional_computed=False,
        physical_coupled_inverse_proved=False,nonlinear_solution_constructed=False,semiclassical_validity_proved=False,
        physical_stability_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
