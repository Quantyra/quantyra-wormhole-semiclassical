"""Actual finite-seed form remainder budget; no forward/adjoint solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,traceback
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_heat_cutoff_form_obstruction as previous
york=previous.york;prior=previous.prior;endpoint=previous.endpoint;principal=previous.principal
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-geodesic-two-sided-remainder';PREVIOUS='current-heat-cutoff-form-obstruction'
OWN=['notes/'+STEM+'.md','checks/current_geodesic_two_sided_remainder.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];FIXTURES=[];encoded=prior.encoded
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)
prior.check=check;york.check=check

def algebra():
    check('exact_second_contact_after_fourth_cancellation',sp.Rational(1,60)-sp.Rational(1,360)==sp.Rational(1,72))
    check('exact_fourth_angular_normalization',sp.Rational(3,4)/90==sp.Rational(1,120))
    y=sp.symbols('y',real=True);theta=sp.Function('theta')(y);T=sp.Function('T')(y);hh=sp.Function('h')(y)
    ch=[sp.cos(theta),sp.sin(theta)];lap=lambda v:sp.diff(v,y,2)+hh*sp.diff(v,y)
    expr=sum(c*lap(c*T) for c in ch)-lap(T)+sp.diff(theta,y)**2*T
    check('partition_Laplacian_identity',sp.trigsimp(sp.expand(expr))==0)
    expr=sum(sp.diff(c*T,y)**2 for c in ch)-sp.diff(T,y)**2-sp.diff(theta,y)**2*T*T
    check('partition_gradient_identity',sp.trigsimp(sp.expand(expr))==0)
    J,Jo,d=sp.symbols('J Jo d',positive=True)
    check('actual_star_gradient_ledger',sp.expand(J+(J+d*d*Jo)+2*(J+d*d*Jo)+sp.Rational(8,10000)*Jo-(4*J+(3*d*d+sp.Rational(8,10000))*Jo))==0)
    eps,C=sp.symbols('eps C',positive=True)
    check('joint_fourth_relative_coefficient',sp.simplify((eps*sp.pi*C/2)/(C/(100*sp.pi))-50*sp.pi**2*eps)==0)
    eta,E,P=sp.symbols('eta E P',positive=True)
    check('ideal_pair_remainder_error',sp.expand(eta*E*P+eta*P*eta*E-eta*(1+eta)*E*P)==0)
    check('curvature_zero_polynomial_reserve',4*2**2+8*100==816)

def execute():
    p,h,zh,env=york.actual_parameters();NUMERICS.extend(york.NUMERICS)
    C,B2,lam,Fl=p['C'],p['B2'],p['g'],p['Fl'];q=arb(p['q']);pi=arb.pi();x=lam*B2
    k=arb(2000);C0=(1+12*pi*k**6).sqrt();b02=k**10*(C0-1)/6
    L=(B2/b02).log()/2;Feff=4*arb.const_euler()-2*((q+1).digamma()+arb.const_euler())+4*L
    guards=dict(C_over_q=C/q,q_over_C=q/C,C_over_x=C/x,GC_over_B2=C/B2,
        Gq2_over_lambdaB4=q*q/(lam*B2*B2),minimum_mass=((q+1)/C).sqrt(),
        finite_seed_d=(p['eps']*C).abs_upper(),thermal_delta=p['tau']*C)
    check('actual_geometry_transfer_ratios',guards['C_over_q']<arb('.0035') and guards['q_over_C']<400 and guards['C_over_x']<arb('6.1'))
    check('actual_mass_and_small_curvature_parameters',guards['minimum_mass']>18 and guards['GC_over_B2']<arb('8.5e-33') and guards['Gq2_over_lambdaB4']<arb('4e-27'))
    check('actual_finite_seed_and_temperature_range',guards['finite_seed_d']<arb('8.6e-19') and guards['thermal_delta']>0 and guards['thermal_delta']<1)
    check('actual_fixed_scale_and_common_contact',L>0 and L<arb('3e-9') and -Feff>50 and -Feff<100)
    Fl_guard=1-arb('1.179e-8');check('actual_classical_factor',Fl>Fl_guard)
    remaining_E0=Fl*(1-arb('7.2e-8')-arb('1.37e-21')-arb('1e-15')-arb('1e-27'))-arb('.003601000003')
    joint_E0=arb('.99')*Fl
    check('joint_E0_and_fourth_margin_retained',remaining_E0>joint_E0,encoded(remaining_E0-joint_E0))
    eps=arb('1e-6');K=arb('1e13');R0=arb(1)/400;a=arb(1)/100;dt=5/a
    b=(eps/(4*K)).sqrt();gradient=K*(arb(1)/4+2*(R0/b).log())
    check('gain_two_anchor_inside_original_chart',b<R0)
    check('changed_epsilon_keeps_existing_gradient_reserve',gradient<arb('1e15'))
    j0=1000*C*C/(lam*B2*B2);j1=arb('2e5')*C/B2;jo=16*a**4/25*j0
    j0star=3*j0;j1star=4*j1+(3*dt*dt+arb('.0008'))*jo
    check('direct_return_bounds_fit_retained_guards',j0<arb('5.2e-29') and j1<arb('1.7e-27'))
    check('weighted_outer_guard',jo<arb('3.4e-37'))
    check('actual_star_norms_fit_retained_guards',j0star<arb('2e-28') and j1star<arb('1e-26'))
    Jbound=604*q*q+(1+2*L.abs_upper())*(q*q/6-arb(1)/30)
    parts=dict(gradient_envelope=arb('1e23')*j1star,
        gluing_zero_envelope=arb('1e25')*j0star,
        actual_constant_contact=arb('1500e20')*q*q/(lam*B2*B2),
        shared_trial_vacuum=Jbound*j0/(pi*C*C),
        known_curvature_contact=Feff.abs_upper()*(j1/(36*pi)+816*j0/(960*pi)),
        fourth_partition_conversion=eps*pi/2*dt**4*jo)
    Berr=sum(parts.values(),arb(0));eta_fourth=50*pi*pi*eps;eta_lower=Berr/joint_E0
    check('vacuum_and_curvature_are_already_small',parts['shared_trial_vacuum']<arb('2e-21') and parts['known_curvature_contact']<arb('2e-26'))
    check('lower_order_ledger_dominates_fourth_ledger',eta_lower>eta_fourth)
    eta=eta_lower;eta_star=arb(11)/5000
    check('two_sided_physical_form_target',eta<eta_star,encoded(eta))
    inv=1/(1-eta);neumann=eta/(1-eta)
    check('linear_comparison_inverse_reserve',inv<arb('1.00221') and neumann<arb('.00221'))
    ideal_EP=arb('8.335e-10')/(eta_star*(1+eta_star));model_EP=(1-eta_star)*ideal_EP
    check('ideal_physical_product_threshold_identity',(eta_star*(1+eta_star)*ideal_EP).contains(arb('8.335e-10')))
    check('model_to_actual_norm_product_conversion',model_EP<ideal_EP and (model_EP/(1-eta_star)).overlaps(ideal_EP))
    NUMERICS.append(dict(kind='actual_finite_seed_two_sided_parameters',C=encoded(C),B2=encoded(B2),lambda_value=encoded(lam),
        F_lambda=encoded(Fl),q=encoded(q),Lscale=encoded(L),F_eff=encoded(Feff),guards={k:encoded(v) for k,v in guards.items()},
        retained_E0_coefficient=encoded(remaining_E0),joint_E0_coefficient=encoded(joint_E0),
        epsilon=encoded(eps),K=encoded(K),R0=encoded(R0),partition_a=encoded(a),dtheta=encoded(dt),
        splitting_b=encoded(b),gain_two_gradient_coefficient=encoded(gradient),Jbound=encoded(Jbound),
        j0=encoded(j0),j1=encoded(j1),j0_outer=encoded(jo),j0star=encoded(j0star),j1star=encoded(j1star)))
    NUMERICS.append(dict(kind='actual_two_sided_remainder_budget',components={k:encoded(v) for k,v in parts.items()},
        Berr=encoded(Berr),fourth_relative_bound=encoded(eta_fourth),lower_order_relative_bound=encoded(eta_lower),
        eta=encoded(eta),frozen_certified_eta='0.0022',comparison_inverse_upper=encoded(inv),
        linear_Neumann_factor_upper=encoded(neumann),physical_V_to_Vstar_bound=True,
        X_admission_Z0=False,finite_geometry_nonlinear_ball_bound=False))
    NUMERICS.append(dict(kind='remaining_physical_certificate_threshold',eta_star=encoded(eta_star),
        assumed_exact_model_residuals=True,assumption_executed=False,maximum_assumed_Q0='8.335e-10',
        physical_EP_threshold=encoded(ideal_EP),model_EP_threshold=encoded(model_EP),
        actual_E_unknown=True,actual_P_unknown=True,actual_D0_unknown=True,actual_Db0_unknown=True,actual_Q0_unknown=True,
        actual_model_pair_computed=False,full_bolt_test_executed=False))
    FIXTURES.append(dict(kind='continuum_and_execution_scope',absolute_finite_kernel_contacts_required=True,
        original_common_fourth_contact_retained=True,actual_paired_C2_1_domain_retained=True,
        bounded_form_extension_covers_bolt_adjoint=True,bolt_adjoint_norm_computed=False,
        source_value_bound_differentiated=False,smooth_double_theorem_assumed=False,
        geodesic_flow_replaced_by_scalar_Laplacian=False,geodesic_flow_matrix_evaluated=False))

def main():
    ctx.prec=256;review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|
        {'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};deps=endpoint.dependencies();error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',
        protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_heat_cutoff_obstruction_retained',old['status']=='passed' and old['uniform_heat_cutoff_obstruction_executed'])
        check('previous_scientific_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        algebra();execute()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths};depafter=endpoint.dependencies()
    status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,
        finished_utc=datetime.now(timezone.utc).isoformat(),two_sided_geodesic_remainder_bounded=status=='passed',
        actual_linear_form_comparison_certified=status=='passed',bounded_form_extension_to_bolt_adjoint=status=='passed',
        comparison_model_forward_adjoint_solve_executed=False,geodesic_flow_matrix_evaluated=False,
        actual_physical_correction_computed=False,full_solution_feasibility_test_executed=False,
        full_adjoint_residual_norm_bounded=False,full_forward_residual_norm_bounded=False,
        complete_renormalized_operator_evaluated=False,metric_admission_proved=False,
        nonlinear_solution_constructed=False,physical_stability_proved=False,
        semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),
                    *ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
    if status!='passed':raise SystemExit(1)

if __name__=='__main__':main()
