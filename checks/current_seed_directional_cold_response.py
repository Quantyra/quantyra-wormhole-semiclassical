"""Current reference path and massive cold direction; not a full-Hessian solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback, math
import sympy as sp
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='current-seed-directional-cold-response'; PREVIOUS='finite-seed-local-bolt-feedback'
OWN=['notes/'+STEM+'.md','checks/current_seed_directional_cold_response.py','evidence/'+STEM+'-review.json']
ROWS=[]; FIXTURES=[]

def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok: raise AssertionError(name)

def exact_profiles():
    z,u,d,beta,c,e=sp.symbols('z u d beta c e')
    hm=Q(2000)**-10; bounds=[Q(1)/Q('8e9'),Q('8.6e-19'),7*hm**2,Q('4.91e-10'),Q('3.4e-22')]
    params=[u,d,beta,c,e]
    def cb(p):
        return sum(abs(Q(int(sp.numer(a)),int(sp.denom(a))))*math.prod(v**n for v,n in zip(bounds,exps)) for exps,a in sp.Poly(sp.expand(p),*params).terms())
    def ac(p): return sum(cb(a) for a in sp.Poly(sp.expand(p),z).all_coeffs())
    dt=lambda p:sp.expand(-z*z*sp.diff(p,z))
    f=z**-2+u-(4+2*u+beta)*z+(3+u+d*u)*z*z-2*d*u*z**3+d*u*z**4
    rn=-12+d*u*(4*z**5-6*z**6)
    tn=(3+u+d*u)*z**4-d*u*(5*z**5-sp.Rational(9,2)*z**6)
    hr=sp.expand(-2*rn*tn+2*(dt(f)/2+2*z*f)*dt(rn))
    ht=sp.expand(hr+2*f*dt(dt(rn)))
    # Independent seed Ricci formula, including scalar subtraction.
    rictt=-(dt(dt(f))/2+z*dt(f))
    check('seed_Ricci_tracefree_identity',sp.simplify(rictt-rn/4+tn)==0)
    ht_direct=2*rn*rictt-rn**2/2+2*(f*dt(dt(rn))+(dt(f)/2+2*z*f)*dt(rn))
    check('independent_Htt_formula',sp.simplify(ht-ht_direct)==0)
    qseed=3+u+d*u; Z=qseed-e; fl=1-24*c
    charge=24*c*qseed+(1-24*c)*e
    check('effective_charge_identity',sp.expand(qseed-fl*Z-charge)==0)
    J=charge*z**4-4*d*u*z**5+3*d*u*z**6
    profiles=[sp.expand((c*ht-J)/2),sp.expand((-c*hr+J)/2)]
    for label,p,sgn in zip(['rho','p'],profiles,[-1,1]):
        check('leading_cancellation_'+label,sp.expand(p.subs(d,0)-sgn*(1-24*c)*e*z**4/2)==0)
        coef=ac(p); der=ac(-z**3*sp.diff(p/z**2,z)); low=sp.expand(p/z**2)
        check('normalized_local_coefficient_'+label,coef<Q('1.701e-22'),float(coef))
        check('energy_derivative_coefficient_'+label,der<Q('3.401e-22'),float(der))
        powers=[m[0] for m,a in sp.Poly(low,z).terms()]
        check('local_decay_and_degree_'+label,min(powers)>=2 and max(powers)<=10)
        local=Q('1.052')*Q(2000)**8*ac(low)
        mass=Q('1.052')*Q(2000)**8*sum(cb(a)/Q(m[0]-1) for m,a in sp.Poly(low,z).terms())
        thermal=Q('1.052')*Q(2000)**8*cb((low-low.subs(beta,0)).subs(z,1))/(hm/8)
        check('local_q_density_and_mass_'+label,max(local,mass)<Q('4.59e4'),dict(coefficient=float(local),mass=float(mass)))
        expected=(108 if label=='rho' else 24)*c*d*u*beta
        check('exact_thermal_offset_'+label,sp.simplify((p-p.subs(beta,0)).subs(z,1)-expected)==0)
        check('local_thermal_over_gap_'+label,thermal<Q('1e-39'),float(thermal))
    check('Hq_conversion',Q('3.14160')/(Q('.99999998')*Q('2.99'))<Q('1.052'))
    for label,k in [('rho',31),('p',29)]:
        value=(k*Q('1e15')+Q('1e6'))*Q('1.001')/Q(2000)**8+Q('3.401e-22')
        check('complete_Q_'+label,value<Q('1.3e-10'),float(value))
    # Full trace high-jet input retains the previous exact polynomial.
    rr=4*z**5-6*z**6; cn=sp.expand(-z**3*sp.diff(f,z)+rn/6)
    ps=2*(6+d*u)**2+sp.Rational(22,3)*(6+d*u)*u+sp.Rational(22,3)*u*u+22
    bp=sp.expand(7*cn**2+sp.Rational(22,3)*tn**2-sp.Rational(11,72)*(rn**2-144)+ps*rr/2)
    aa=-10*z**4-20*z**5+30*z**6; lr=sp.expand(z**4*sp.diff(f*sp.diff(rr,z),z))
    ss=Q('1.001')/(60*Q('3.14159')*Q(2000)**8); eta=hm/(60*Q('3.14159'))
    trace=ss*ac(aa/z**2)+eta*ac(bp/z**2)+Q('2.946e-9')*bounds[0]*bounds[1]*ac(lr/z**2)
    check('local_full_trace_spatial_reserve',trace<Q('1e-24'),float(trace))
    check('trace_product_zeros',all(sp.simplify(p.subs({z:1,beta:0}))==0 for p in [aa,bp,lr]))
    FIXTURES.append(dict(fixture='exact_current_polynomials',full_coupled_defect_computed=False))

def tensor_and_integral():
    z,theta=sp.symbols('z theta',real=True); b=sp.Function('b')(z); L=sp.Function('L')(z)
    coords=[z,theta]; g=sp.diag(1/b,b); inv=g.inv()
    conn=[[[sp.simplify(sum(inv[i,k]*(sp.diff(g[k,j],coords[l])+sp.diff(g[k,l],coords[j])-sp.diff(g[j,l],coords[k])) for k in range(2))/2) for l in range(2)] for j in range(2)] for i in range(2)]
    hs=sp.Matrix(2,2,lambda i,j:sp.diff(L,coords[i],coords[j])-sum(conn[k][i][j]*sp.diff(L,coords[k]) for k in range(2)))
    frame=[-sp.sqrt(b),1/sp.sqrt(b)]
    cov3=lambda i,j,k:sp.simplify(frame[i]*frame[j]*frame[k]*(sp.diff(hs[j,k],coords[i])-sum(conn[l][i][j]*hs[l,k]+conn[l][i][k]*hs[j,l] for l in range(2))))
    p=sp.sqrt(b); bz=sp.diff(b,z); bzz=sp.diff(b,z,2); l1=sp.diff(L,z);l2=sp.diff(L,z,2);l3=sp.diff(L,z,3)
    check('radial_tensor_Hessian',sp.simplify(frame[0]**2*hs[0,0]-(b*l2+bz*l1/2))==0 and sp.simplify(frame[1]**2*hs[1,1]-bz*l1/2)==0)
    targets={(0,0,0):-p*(b*l3+sp.Rational(3,2)*bz*l2+bzz*l1/2),(0,1,1):-p*(bzz*l1+bz*l2)/2,(1,0,1):-p*bz*l2/2,(1,1,0):-p*bz*l2/2}
    check('ordered_third_tensor_no_bolt_division',all(sp.simplify(cov3(*ijk)-target)==0 for ijk,target in targets.items()))
    cap=Q('4.59e4'); s=Q('.001')
    global_jets=[cap,10*cap/100,(90+50)*cap/20000,(720+15*90+30)*cap/6000000]
    cold_jets=[10*cap+Q('1e-39'),10*cap/100,(90*s+30)*cap/20000,(720*s*s+810*s+30)*cap/6000000]
    check('local_C3_normal_reserve',2*max(global_jets)<Q('1e8'))
    check('local_cold_C3_reserve',2*max(cold_jets)<Q('1e8'))
    # Positive primitive for the total-variation estimate, not a value-to-TV inference.
    t=sp.symbols('t',positive=True); primitive=sp.log((t-1)/t)+1/t
    check('positive_lapse_integral_primitive',sp.simplify(sp.diff(primitive,t)-1/(t*t*(t-1)))==0)
    hm=Q(2000)**-10
    mass=80*Q('1.001')/Q('.99')*Q('1e15')/Q(2000)**8
    lll=80*Q('1.001')/Q('.99')*(Q('.0035')/4)/Q(2000)**8
    weyl=hm*((Q('2.01')+Q(19,8)*Q('3.01'))/(15*Q('3.14159'))+3*Q('8.6e-19')/(2*Q('3.14159')*Q('8e9')))
    check('massive_absolute_lapse',mass<Q('3.2e-10'),float(mass))
    check('LLL_absolute_lapse',lll<Q('3e-28'),float(lll))
    check('Weyl_absolute_lapse',weyl<2*hm,float(weyl))
    check('complete_absolute_lapse',mass+lll+weyl+Q('9.182e-37')<Q('3.3e-10'))
    check('low_geometry_transfer',Q('3.3e-10')+3*Q('1.3e-10')+80*Q('3.33e-10')*hm<Q('7.5e-10'))
    check('current_C3_inputs_fit',Q('5e17')<Q('1e20') and Q('6e15')<Q('1e26'))
    FIXTURES.append(dict(fixture='independent_tensor_connection_and_positive_lapse_primitive',generic_elliptic_C3_to_C5_asserted=False))

def endpoint_and_targets():
    base=Q('1e23')*2*Q('8.5e-33')*Q('1.000001')*Q('1e-6')
    check('finite_endpoint_to_seed_derivative',base<Q('1.71e-15') and 10*base<Q('1.71e-14'),float(base))
    seed=Q('1.0001e15')*Q('2.72e-27');contact=2*Q('1e-6')*Q('2.721e-12')
    check('cold_seed_value',seed<Q('2.721e-12'),float(seed))
    check('fixed_density_value_bound',Q('1.71e-15')+contact<Q('1.72e-15'))
    osc=2*Q('16e-6')*Q('2.721e-12')+Q('1.4e-12')*81*Q(2000)**-10
    check('fixed_density_bolt_contact',osc<Q('1e-16'),float(osc))
    check('fixed_density_bolt_bound',Q('1.71e-14')+osc<Q('1.72e-14'))
    mp.mp.dps=70
    # This is only a manufactured analytic family with the same limit/contact structure.
    ell=sp.symbols('ell');U=sp.Rational(2,100);S0=sp.Rational(3,10);k=sp.Rational(1,1000)
    S=S0+k*ell/(1+ell);density=sp.exp(2*ell*U)*S
    exact=sp.diff(density,ell).subs(ell,0)
    check('manufactured_density_contact_identity',exact==k+2*U*S0)
    fun=sp.lambdify(ell,density,'mpmath');target=mp.mpf(str(exact.evalf(72)))
    errors=[]
    for h in [mp.mpf('1e-12'),mp.mpf('5e-13')]:
        errors.append(abs((-3*fun(0)+4*fun(h)-fun(2*h))/(2*h)-target))
    check('one_sided_endpoint_derivative_fixture',max(errors)<mp.mpf('1e-24') and errors[1]<mp.mpf('.26')*errors[0],list(map(str,errors)))
    rho=Q('1e-8');Y=Q('5e-9');Z0=Q('.25');Z2=Q('5e7')
    check('preserved_neighborhood_targets',Y+Z0*rho+Z2*rho*rho/2==rho and Z0+Z2*rho==Q('.75'))
    check('regional_C3_inverse_allocations',Q('1.632e-11')*Q('61.2745098039')+Q('1.36e-9')*Q('1.4705882352')+Q('2e-9')<=Y)
    check('admission_thresholds',Q('1.3e-4')+rho*87000==Q('1e-3') and Q('4.9e-7')+rho==Q('5e-7') and Q('1.4e-7')+rho==Q('1.5e-7') and rho*Q('.025')==Q('2.5e-10'))
    check('coupled_trace_target',Q(4,3)*Q('1e18')*(Q('1.1e-31')+Q('1.25014e-27'))==Q('1.667e-9'))
    FIXTURES.append(dict(fixture='manufactured_endpoint_density_and_admission_arithmetic',actual_quantum_determinant_sampled=False,square_root_norm_convergence_asserted=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    try:
        check('previous_signed_local_test_passed',old['status']=='passed' and old['local_R2_reference_bolt_feedback_bounded'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        exact_profiles();tensor_and_integral();endpoint_and_targets()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),current_reference_path_C5_hypotheses_checked=error is None,massive_directional_cold_value_and_bolt_bound=error is None,directional_response_C3_bounded=False,full_local_action_defect_pairing_evaluated=False,full_quantum_shape_defect_evaluated=False,full_defect_bolt_pairing_evaluated=False,actual_physical_correction_computed=False,raw_C5_admission_certified=False,nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__': main()
