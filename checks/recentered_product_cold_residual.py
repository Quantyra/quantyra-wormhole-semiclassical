"""Frozen product enclosure and exact finite-bolt residual budgets.

The pressure quadrature illustrates the product; it is not an interval
certificate or a sample of the finite HH quantum source.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback, sys, importlib.util
import sympy as S
import mpmath as mp
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
STEM = 'recentered-product-and-cold-residual'
PREVIOUS = 'recentered-energy-and-cold-admission'
OWN = ['notes/'+STEM+'.md', 'checks/recentered_product_cold_residual.py', 'evidence/'+STEM+'-review.json']
ROWS = []; FIXTURES = []

def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()

def check(name, ok, detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name, passed=bool(ok), detail=str(detail)))
    if not ok: raise AssertionError(name)

def product_identities():
    a,e,b,F,eta,v,A,av,a2 = S.symbols('alpha e beta F eta v A av a2')
    chi=a*eta+v
    # Independent expansion of the original orbit equation divided by X0.
    raw=-a*eta+F*(1+a)*(2*eta+eta**2)/2+b*(2*chi+chi**2-a*a*(2*eta+eta**2))/2-e*((1+eta)**2*A-av)
    target=F*eta+b*(1+a*eta)*v+b*v*v/2+F*(1+a)*eta**2/2-e*((1+eta)**2*A-av)
    check('product_orbit_exact_recenter',S.expand((raw-target).subs(F,1-b*(1-a)))==0)
    X,Y,lam,c,G,Qm=S.symbols('X Y lambda c G Qm')
    X0=Y+6*lam;Fl=1-384*S.pi*G*c*lam
    baseline=-Y-3*lam+Fl*(Y+3*lam)+32*S.pi*G*c*(X0**2-Y**2)
    check('local_RN_vacuum_tensor_exact',S.expand(baseline)==0)
    pA,pS,l0=S.symbols('pA pS lambda0')
    EA=-Y-3*l0+Qm*Y**2+32*S.pi*G*c*(X**2-Y**2)-8*S.pi*G*pA
    ES=X-3*l0-Qm*Y**2-32*S.pi*G*c*(X**2-Y**2)-8*S.pi*G*pS
    check('full_product_trace_R2_cancels',S.expand(EA+ES-(X-Y-6*l0-8*S.pi*G*(pA+pS)))==0)
    C,x=S.symbols('C x',positive=True)
    Cs=C*(1+chi)/(1+eta);xs=x/(1+eta)
    check('independent_radius_seed_defect',S.simplify((Cs-1-6*xs-C*v/(1+eta)).subs({a:1/C,x:(C-1)/6}))==0)
    check('seed_relative_defect',S.simplify((C*v/(1+eta))/Cs-v/(1+chi))==0)
    FIXTURES.append(dict(fixture='independent_full_product_algebra',R2_trace_cancellation=True,physical_flux_unchanged=True))

def exact_budgets():
    al=Q('2.1e-11');be=Q('1.18e-8');Fl=1-Q('1.179e-8')
    elo=Q('1.5e-26');ehi=Q('1.6e-26');av=Q('4e-10');B=Q(3334);p=Q('1e-20')
    amin=Q('.0021');amax=Q('.00212');da=Q('2e-8');q=Q(2000**4)
    pilo=Q('3.14');pihi=Q(22,7);rm=Q('.0036')
    check('alpha_bound',Q(1,6*2000**3)<al)
    check('exact_e_identity_lower',4*Q(1,2000**8)*Fl>elo)
    check('exact_e_identity_upper',4*Q(1,2000**8)*(1+al)<ehi)
    check('beta_guard',Q('1.179e-8')/(1-al)<be)
    check('vacuum_pressure_guard',11*(Q('.0035')/6)**2/(960*pilo**2)<av)
    check('heat_trace_lower',(15-3/q**2-5*rm/q-3*rm**2)/(720*pihi**2)>amin)
    check('heat_trace_upper',15/(720*pilo**2)<amax)
    check('heat_trace_derivative',(5*rm/q+6*rm**2)/(720*pilo**2)<da)
    check('pressure_and_scale_guard',100*(1+32+Q('2e-8')+2*p)<B)
    ch=(1+al)*p
    veta=ehi*((1+p)**2*amax/2+2*av)
    etamap=ehi*((1+p)**2*B+av)/Fl+be*((1+al*p)*p+p*p/2)/Fl+(1+al)*p*p/2
    check('product_box_eta_preserved',etamap<p)
    check('product_box_v_preserved',veta<p)
    # Explicit Lipschitz row bounds, including coordinate denominators.
    dr_eta=al/(1-ch)+1/(1-p);dr_v=1/(1-ch)
    l_eta=1/(1-p)
    row_eta=ehi*(2*(1+p)*B+(1+p)**2*(B*(dr_eta+dr_v)+l_eta))/Fl+be*(al*p+1+al*p+p)/Fl+(1+al)*p
    row_v=ehi*(2*(1+p)*amax+(1+p)**2*da*(dr_eta+dr_v))/2
    check('contraction_row_eta',row_eta<Q('1.2e-8'))
    check('contraction_row_v',row_v<Q('1e-27'))
    vmin=elo*((1-p)**2*amin/2-2*av)
    check('fixed_trace_upper',veta<Q('1.7e-29'))
    check('fixed_trace_signed_lower',vmin>Q('1.57e-29'))
    ve=Q('1.7e-29');re=Q('5.4e-23')
    eta_fixed=(ehi*((1+p)**2*B+av)+be*((1+al*p)*ve+ve*ve/2))/Fl+(1+al)*p*p/2
    check('fixed_eta_upper',eta_fixed<re)
    check('squared_sphere_radius_bound',re/(1-re)<Q('5.5e-23'))
    ce=al*re+ve
    check('squared_orbit_radius_lower',(Q('1.57e-29')-al*re)/(1+ce)>Q('1.56e-29'))
    check('squared_orbit_radius_upper',ce/(1-ce)<Q('1.71e-29'))
    # Lipschitz remainders about exact product pressure, not quadrature.
    log_r=ce/(1-ce)+re/(1-re);log_y=re/(1-re)
    delta_a2=da*log_r
    dv=ehi*((2*re+re*re)*amax+(1+re)**2*delta_a2)/2
    check('exact_v_first_order_remainder',dv<Q('2e-51'))
    delta_A=B*log_r+log_y
    de=(ehi*((2*re+re*re)*B+(1+re)**2*delta_A)+be*(dv+al*re*ve+ve*ve/2))/Fl+(1+al)*re*re/2
    check('exact_eta_first_order_remainder',de<Q('2e-44'))
    check('sphere_radius_first_order_remainder',de+re*re/(1-re)<Q('3e-44'))
    check('orbit_radius_first_order_remainder',al*de+dv+ce*ce/(1-ce)<Q('3e-51'))
    check('seed_d_guard',Q('5e10')*ve/(1-re)<Q('8.6e-19'))
    check('seed_relative_d_guard',ve/(1-ce)<Q('1.71e-29'))
    check('trace_vs_cold_gap_lower',Q('1.57e-29')/(2*Q(3,2000**10))>2500)
    check('R2_local_coefficient_guard',16*pihi*Q('1e22')*Q(1,2000**10)*(1+Q(1,8*10**9)+Q('1e-12')/(8*10**9))<Q('1e-9'))
    check('vacuum_mismatch_relative',11*Q(1,2000**10)/(360*pilo)*Q(1,5)<Q('2e-36'))
    check('local_horizon_residual',68*Q('1.71e-29')+13*Q('2e-36')<Q('1.2e-27'))
    check('massive_local_subtotal',Q('1.2e16')*81*Q('.00352')**3<Q('4.24e10'))
    check('massive_gluing_subtotal',Q('4e12')*81*Q('.00352')<Q('1.141e12'))
    check('massive_local_restoration_subtotal',Q('1e12')*81*Q('.00352')**2<Q('1.004e9'))
    check('full_product_change_allowance',B*4+2+4<Q('2e4'))
    check('full_horizon_source_allowance',Q('1e14')+Q('2e4')<Q('1.001e14'))
    full=Q('1.2e-27')+Q('1.004')*Q('1.001e14')*ehi
    check('full_horizon_residual',full<Q('2e-12'))
    check('conditional_direct_response_target',Q('5e-9')/Q('2e-12')==2500)
    check('nonlinear_preserved_radius',Q('5e-9')+Q('.25')*Q('1e-8')+Q('5e7')*Q('1e-8')**2/2==Q('1e-8'))
    check('nonlinear_contraction',Q('.25')+Q('5e7')*Q('1e-8')==Q('.75'))
    FIXTURES.append(dict(fixture='rational_product_and_cold_budgets',eta_upper='5.4e-23',v_interval=['-1.7e-29','-1.57e-29'],
        full_horizon_residual_upper_exact=str(full),direct_response_factor_target='<=2500 conditional on a valid complete source inequality',admission_certified=False))

def horizon_identities():
    t,V,x,d,a,Qp,l0,wc=S.symbols('t V x d a Qp lambda0_b2 wc')
    Qc=1+3*x+d;C=1+6*x+d
    f=1+x*t*t-2*a/t+Qc/t**2-2*d/t**3+d/t**4
    av=(V+x*V**3+Qc/V-2*d/V**2+d/V**3)/2
    R=S.expand(-S.diff(f,t,2)-4*S.diff(f,t)/t-2*(f-1)/t**2)
    expected=-12*x+d*(4/t**5-6/t**6)
    check('full_seed_scalar_curvature',S.simplify(R-expected)==0)
    check('actual_bolt_root',S.simplify(f.subs(t,V).subs(a,av))==0)
    j=1+3*x*V**2-Qc/V**2+4*d/V**3-3*d/V**4
    w=-1+2*Qc/V**2-10*d/V**3+9*d/V**4
    check('exact_bolt_temperature',S.simplify((t*S.diff(f,t)).subs(t,V).subs(a,av)-j)==0)
    check('exact_bolt_orbit_curvature',S.simplify((t*t*S.diff(f,t,2)/2).subs(t,V).subs(a,av)-w)==0)
    rb=expected.subs(t,V);rp=S.diff(expected,t).subs(t,V)
    RA=-(w+j)/V**2;RS=(1-j)/V**2
    check('curvature_trace_at_bolt',S.simplify(2*(RA+RS)-rb)==0)
    HA=S.expand(2*rb*RA-rb*rb/2+j*rp/V)
    HS=S.expand(2*rb*RS-rb*rb/2+2*j*rp/V)
    LA=(j-1)/V**2-3*l0+Qp/V**4+wc*HA
    LS=(w+j)/V**2-3*l0-Qp/V**4+wc*HS
    ta=d*(4/V**5-3/V**6-1/V**4)-3*(l0-x)*(1-1/V**4)+wc*(HA-HA.subs(V,1)/V**4)
    ts=d*(-6/V**5+6/V**6)-3*(l0-x)*(1-1/V**4)+wc*(HS-HS.subs(V,1)/V**4)
    check('finite_local_orbit_cancellation',S.simplify(LA-LA.subs(V,1)/V**4-ta)==0)
    check('finite_local_sphere_cancellation',S.simplify(LS-LS.subs(V,1)/V**4-ts)==0)
    check('RN_R2_orbit_scaling',S.simplify((HA-24*x*Qc/V**4).subs(d,0))==0)
    check('RN_R2_sphere_scaling',S.simplify((HS+24*x*Qc/V**4).subs(d,0))==0)
    z=S.symbols('z');ledgers=[]
    for name,expr in [('orbit',HA-24*x*Qc/V**4),('sphere',HS+24*x*Qc/V**4)]:
        pp=S.Poly(S.expand(expr/d).subs(V,1/z),z,x,d)
        led={0:0,1:0,2:0}
        for (nz,nx,nd),coeff in pp.terms():
            assert (nx,nd) in [(0,0),(1,0),(0,1)]
            led[1 if nx else 2 if nd else 0]+=abs(coeff)*(nz+4)
        check('R2_polynomial_majorant_'+name,max(led.values())<100000,led)
        ledgers.append(dict(channel=name,coefficients=list(map(str,led.values()))))
    # Conservative endpoint inequalities use only x>=1000, |d|<1e-12.
    vm=Q('1.001');dm=Q('1e-12');xm=Q(1000)
    check('thermal_j_derivative_lower',6*xm+2*(1+3*xm-dm)/vm**3-24*dm>(1+6*xm+dm))
    check('thermal_denominator_guard',vm*vm<Q('1.003'))
    check('local_EH_difference_A',20+18+4<66)
    check('local_EH_difference_S',30+36==66)
    check('local_horizon_rounding_d',Q('1.003')*(66+Q('1e-4'))<68)
    check('local_horizon_rounding_vacuum',Q('1.003')*12<13)
    FIXTURES.append(dict(fixture='independent_R2_field_equations_and_polynomial_bounds',R2_ledgers=ledgers,actual_HH_source_computed_numerically=False))

def numerical_illustration():
    spec=importlib.util.spec_from_file_location('retained_product',ROOT/'checks/uniform_dirac_reference.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    vacuum=module.vacuum
    records=[]
    mp.mp.dps=90;k=mp.mpf(2000);q=k**4;g0=k**-10
    g=2*g0/(1+mp.sqrt(1+44*g0/(360*mp.pi)))
    F=1-384*mp.pi*mp.mpf('1e22')*g
    C0=mp.sqrt(1+12*mp.pi*k**6);x0=(C0-1)*g/g0/6
    C=mp.sqrt(1+12*mp.pi*g*k**16/F);x=(C-1)/6
    L0=-mp.log(q)+mp.log(x/x0)
    vals=[vacuum(float(q),float(C/q),float(L0),t0=cut) for cut in [.125,.2]]
    check('product_pressure_two_UV_splits',abs(vals[0]['scaled_pA']-vals[1]['scaled_pA'])<1e-9)
    A0=mp.mpf(vals[0]['scaled_pA'])
    for digits in [90,140]:
        mp.mp.dps=digits
        k=mp.mpf(2000);q=k**4;g0=k**-10;g=2*g0/(1+mp.sqrt(1+44*g0/(360*mp.pi)))
        F=1-384*mp.pi*mp.mpf('1e22')*g
        C0=mp.sqrt(1+12*mp.pi*k**6);x0=(C0-1)*g/g0/6
        C=mp.sqrt(1+12*mp.pi*g*k**16/F);x=(C-1)/6
        alpha=1/C;e=8*mp.pi*g*q*q/(C*x);beta=64*mp.pi*mp.mpf('1e22')*g*C/x
        av=-11*x*x/(960*mp.pi**2*q*q);r=C/q
        a2=(15-3/q**2-5*r/q-3*r*r)/(720*mp.pi**2)
        v1=e*(-a2/2-2*av);eta1=(e*(A0-av)-beta*v1)/F;chi1=alpha*eta1+v1
        check('numeric_e_two_form_identity_'+str(digits),abs(e/(4*k**-8*F*(1+alpha))-1)<mp.mpf('1e-80'))
        xs=x/(1+eta1);Cs=C*(1+chi1)/(1+eta1);ds=C*v1/(1+eta1)
        check('seed_metric_shape_guard_'+str(digits),abs(ds)<mp.mpf('8.6e-19'))
        check('seed_independent_curvature_guard_'+str(digits),abs(ds/Cs)<mp.mpf('1.71e-29'))
        b=mp.sqrt(xs);b0=mp.sqrt(x0)
        for zi in ['.25','.265','.28']:
            kap=mp.mpf(zi)*2*k**3*g*C0/b0**3;gamma=kap*b/Cs
            def K(h):
                V=1+h
                return ((1+3*xs)*(V**-2+V**-3)+3*xs*(1+V**-1)-ds*(V**-4-3*V**-5))/(2*Cs)
            h=gamma
            for _ in range(5):h-=(h*K(h)-gamma)/(K(h)+h*mp.diff(K,h))
            check('fixed_original_temperature_'+str((digits,zi)),abs(h*K(h)/gamma-1)<mp.mpf('1e-80'))
            check('actual_cold_guards_'+str((digits,zi)),1/(8*k**10)<h<1/k**10 and mp.log(1/h)<80)
            V=1+h;Qc=1+3*xs+ds
            j=1+3*xs*V**2-Qc/V**2+4*ds/V**3-3*ds/V**4
            # The exact polynomial expression is intentionally evaluated at
            # enough precision to resolve the thermal cancellation.
            check('independent_temperature_polynomial_'+str((digits,zi)),abs(j/(2*kap*b*V)-1)<mp.mpf('1e-45'))
            records.append(dict(digits=digits,prepared_zi=zi,e=mp.nstr(e,60),beta=mp.nstr(beta,60),
                pressure_A_illustrative=str(A0),sphere_squared_first_shift=mp.nstr(-eta1,60),
                orbit_squared_first_shift=mp.nstr(-chi1,60),product_v_first=mp.nstr(v1,60),
                seed_d_first=mp.nstr(ds,60),h=mp.nstr(h,60),trace_shift_over_thermal_gap=mp.nstr(abs(v1)/(2*h),60)))
    for i in range(3):
        for key in ['e','beta','sphere_squared_first_shift','orbit_squared_first_shift','h','trace_shift_over_thermal_gap']:
            check('cross_precision_'+str(i)+'_'+key,abs(mp.mpf(records[i][key])/mp.mpf(records[i+3][key])-1)<mp.mpf('1e-50'))
    FIXTURES.append(dict(fixture='benchmark_and_two_UV_splits',pressure_records=vals,records=records,
        pressure_error_is_interval_certified=False,full_quantum_source_sampled=False,finite_corrected_solution=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_energy_result',old['status']=='passed' and old['generic_energy_to_admission_certificate_rejected'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        product_identities();exact_budgets();horizon_identities();numerical_illustration()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,
        unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        recentered_product_enclosure_proved=error is None,product_centered_horizon_value_target_met=error is None,
        actual_correction_pointwise_lower_bound_proved=False,admission_norm_inverse_quantified=False,
        inverse_transferred_to_product_centered_seed=False,nonlinear_solution_constructed=False,
        semiclassical_validity_proved=False,physical_stability_proved=False,physical_coupled_inverse_proved=False,
        component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
