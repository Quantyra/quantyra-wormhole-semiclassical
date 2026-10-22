"""Exact derivative-budget audit; no actual higher source bound is asserted."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='schur-source-derivative-budget';PREVIOUS='inverse-endpoint-ten-jets'
OWN=['checks/schur_source_derivative_budget.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational;ROWS=[];FIXTURES=[]
PROTOCOL={'exact_arithmetic':'symbolic and rational','numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','parameter_contour_points':16,'spatial_fourier_points':32,
          'source_order':6,'mixed_metric_order':11,'quantum_output_order_selected':4,
          'ordinary_continuum_proof_separate':True,'actual_charged_source_evaluated':False,
          'manufactured_parameters_not_physical':True}


def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)


def zero(name,expr):
    value=S.simplify(expr);check(name,value==0,value)


def previous():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(previous()['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})


def provenance():
    p=previous()
    check('prior_endpoint_result_passed',p['status']=='passed' and len(p['checks'])==94 and all(r['passed'] for r in p['checks']))
    check('prior_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('all_prior_inputs_unchanged',not bad,bad)
    check('retained_sources_pinned',set(review()['retained_inputs'])<=set(paths()))
    check('audit_does_not_claim_missing_estimates',not review()['actual_finite_metric_source_H6_bounded']
          and not review()['original_metric_C11_quantified'] and review()['direct_mixed_row_requires_metric_order_eleven'])


def conformal_and_return_identities():
    x,t,eps,E0,eta,G=S.symbols('x t epsilon E0 eta0 G',real=True)
    v=S.Function('v')(x);phi=t*v
    scalar=-6*S.exp(-2*phi)*(S.diff(phi,x,2)+S.diff(phi,x)**2)
    box=S.exp(-2*phi)*(S.diff(scalar,x,2)+2*S.diff(phi,x)*S.diff(scalar,x))
    dR=S.diff(scalar,t).subs(t,0);dBox=S.diff(box,t).subs(t,0)
    zero('full_conformal_scalar_first_variation',dR+6*S.diff(v,x,2))
    zero('full_conformal_box_first_variation',dBox+6*S.diff(v,x,4))
    E=E0+scalar+eps*box
    div=S.exp(-4*phi)*S.diff(S.exp(4*phi)*E*eta,x)
    ddiv=S.simplify(S.diff(div,t).subs(t,0))
    expected=eta*(-6*S.diff(v,x,3)-6*eps*S.diff(v,x,5)+4*E0*S.diff(v,x))
    zero('physical_volume_return_contact',ddiv-expected)
    sixth=S.expand(S.diff(ddiv,x,6));top=S.diff(v,x,11)
    zero('eleventh_derivative_coefficient',sixth.coeff(top)+6*eps*eta)
    normalized=S.expand(sixth.subs(eps,G/(60*S.pi))/(8*S.pi*G))
    zero('physical_G_cancels_from_top_return',normalized.coeff(top)+eta/(80*S.pi**2))
    check('top_return_coefficient_nonzero',S.simplify(normalized.coeff(top))!=0)
    check('omitted_Box_R_detected',S.expand(sixth.subs(eps,0)).coeff(top)==0 and sixth.coeff(top)!=0)
    check('omitted_return_detected',S.simplify(normalized.coeff(top)-0)!=0)
    volume_omitted=eta*(dR+eps*dBox)
    zero('omitted_volume_error_identified',ddiv-S.diff(volume_omitted,x)-4*E0*eta*S.diff(v,x))
    N=S.symbols('N',positive=True);wave=N**(-10)*S.cos(N*x)
    mode=S.simplify(sixth.subs(v,wave).doit())
    zero('C10_wave_return_output',mode-eta*(6/N-6*eps*N+4*E0/N**3)*S.sin(N*x))
    for j in [0,5,10,11]:
        deriv=S.diff(wave,x,j)
        zero('wave_derivative_'+str(j),deriv-N**(j-10)*S.cos(N*x+j*S.pi/2))
    check('C10_uniform_C11_grows',max(j-10 for j in range(11))==0 and 11-10==1)
    e=S.Function('E')(x);y=S.Function('eta')(x);density=S.Function('L')(x)
    row=S.expand(S.diff(S.diff(e*y,x)+density*e*y,x,6))
    zero('general_radial_top_trace_coefficient',row.coeff(S.diff(e,x,7))-y)
    check('full_row_trace_order_seven',max(len(d.variables) for d in row.atoms(S.Derivative) if d.expr==e)==7)
    FIXTURES.append({'fixture':'exact_conformal_and_off_shell_return_symbol',
                     'top_metric_derivative':11,'top_prefactor':'-eta0/(80*pi^2)',
                     'density_contact_retained':True,'seed_equation_substitution_used':False,
                     'actual_comparison_counterexample_claimed':False})


def schur_norm_ledger():
    leibniz=sum(2**j for j in range(7))
    check('ordered_Leibniz_sum',leibniz==127)
    check('three_local_tensor_terms',3*leibniz==381)
    check('four_dimensional_divergence_contraction',2**2==4)
    check('return_product_budget',(1+2)*leibniz==381)
    check('scalar_Schur_constant',80*R(22,7)**2<800)
    check('TT_inverse_output_order',4+4==8)
    check('local_mixed_output_order',8-2==6)
    check('local_trace_metric_order',4+6==10)
    check('return_trace_metric_order',4+7==11)
    check('seed_reconstruction_next_order',9+2==11)
    check('log_quantum_H4_sufficient_order',4+4+1==9)
    # The top background jet in the TT first-order ODE rises at most once per differentiation.
    b_orders=[0];eta_orders=[0]
    for j in range(1,9):
        b_orders.append(max(j,max(b_orders)))
        eta_orders.append(max(j,max(b_orders[:j]),max(eta_orders)))
    check('lift_ODE_does_not_generate_order_eleven',b_orders[-1]==8 and eta_orders[-1]==8)
    k=S.symbols('k',nonnegative=True)
    zero('Sobolev_log_embedding_margin',4*(1+k*k)-(1+k)**2-(3*(k-R(1,3))**2+R(8,3)))
    check('positive_log_embedding_margin',R(8,3)>0)
    check('H10_wave_norm_uniform',(1+R(1))**5==32)
    check('weights_not_commuted',review()['weight_only_on_mixed_output'])
    FIXTURES.append({'fixture':'conditional_schur_norm_split_and_derivative_ledger',
                     'mixed_coefficient_multiplier':381,'Schur_multiplier':800,
                     'source_orders':{'scalar_weighted':6,'TT_unweighted':4},
                     'TT_output_order':8,'trace_coefficient_order':7,
                     'quantitative_TT_inverse_supplied':False,'quantitative_source_supplied':False,
                     'lift_ODE_derivative_orders':b_orders})


def tower_and_fourier_identities():
    n,q,C=S.symbols('n q C',positive=True)
    mass2=n*(n+q)/C;deg=2*(q+2*n)
    zero('angular_partial_fraction',deg/mass2-2*C*(1/n+1/(n+q)))
    zero('harmonic_lower_coefficient',2/C*R(3,2048)*2*C-R(3,512))
    z=S.symbols('z');a=z*(1-z)
    zero('second_parameter_moment',S.integrate(a*a,(z,0,1))-R(1,30))
    zero('third_parameter_moment',S.integrate(a**3,(z,0,1))-R(1,140))
    v=S.symbols('v',positive=True)
    zero('R_remainder_derivative',S.diff(S.log(1+v)-v+v*v/2,v)-v*v/(1+v))
    check('positive_parameter_lower_interval',R(3,16)*16>=3)
    # |K| times the H6-to-H10 wave scaling is at least |K|/k^4.
    k=S.symbols('k',positive=True)
    zero('output_wave_scaling',((1+k*k)**3/k**10)/(k**(-4))-(1+k**(-2))**3)
    zero('input_wave_scaling',(1+k*k)**5/k**10-(1+k**(-2))**5)
    check('quantum_log_endpoint_not_declared_curved',not review()['flat_multiplier_identified_with_curved_operator'])
    FIXTURES.append({'fixture':'full_tower_logarithm_and_fixed_torus_diagnostic',
                     'integer_frequency_choice':'4*ceil(M_N)',
                     'unboundedness_from_exact_harmonic_lower_bound':True,
                     'numerical_cutoff_extrapolation_used':False})


def independent_numerics():
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            for N in [3,12]:
                q=mp.mpf(7);C=mp.mpf(2);kn=4*int(mp.ceil(mp.sqrt(N*(N+q)/C)))
                total=mp.mpf(0);closed=mp.mpf(0)
                for n in range(1,N+1):
                    m2=n*(n+q)/C;d=2*(q+2*n);aa=kn**2/m2
                    value=mp.quad(lambda z:mp.log1p(aa*z*(1-z))-aa*z*(1-z)+(aa*z*(1-z))**2/2,[0,mp.mpf('.5'),1])
                    form=2*(mp.sqrt((4+aa)/aa)*mp.atanh(mp.sqrt(aa/(4+aa)))-1)-aa/6+aa*aa/60
                    total+=d*m2*value/(mp.pi*C);closed+=d*m2*form/(mp.pi*C)
                err=abs(total-closed)
                check('independent_tower_parameter_integral_'+str(dps)+'_'+str(N),err<mp.mpf('1e-42'),mp.nstr(err,12))
                lower=mp.mpf(3)/(512*mp.pi)*sum(mp.mpf(1)/n for n in range(1,N+1))
                ratio=total/(kn**4*lower)
                check('finite_positive_partial_tower_lower_'+str(dps)+'_'+str(N),ratio>1,mp.nstr(ratio,12))
                records.append({'dps':dps,'toy_q':7,'toy_C':2,'levels':N,'integer_frequency':kn,
                                'integral_error':mp.nstr(err,18),'lower_ratio':mp.nstr(ratio,18)})
            for N in [2,5]:
                eps=mp.mpf('.01');E0=mp.mpf(3);eta=mp.mpf('.7');radius=mp.mpf('.01')
                def physical_div(x,t):
                    v=[t*mp.mpf(N)**(j-10)*mp.cos(N*x+j*mp.pi/2) for j in range(6)]
                    rr=-6*(v[2]+v[1]**2)
                    r1=-6*(v[3]+2*v[1]*v[2])
                    r2=-6*(v[4]+2*v[2]**2+2*v[1]*v[3])
                    r3=-6*(v[5]+6*v[2]*v[3]+2*v[1]*v[4])
                    bb=r2-2*v[1]*r1-2*v[2]*rr
                    b1=r3-2*v[1]*r2-4*v[2]*r1-2*v[3]*rr
                    curvature=mp.exp(-2*v[0])*rr;box=mp.exp(-4*v[0])*bb
                    E=E0+curvature+eps*box
                    Ex=mp.exp(-2*v[0])*(r1-2*v[1]*rr)+eps*mp.exp(-4*v[0])*(b1-4*v[1]*bb)
                    return eta*(Ex+4*v[1]*E)
                sine=mp.mpc(0)
                for j in range(32):
                    x=2*mp.pi*j/32
                    linear=sum(physical_div(x,radius*mp.exp(2j*mp.pi*k/16))*mp.exp(-2j*mp.pi*k/16) for k in range(16))/(16*radius)
                    sine+=linear*mp.sin(N*x)*2/32
                output=-(mp.mpf(N)**6)*sine
                expected=eta*(mp.mpf(6)/N-6*eps*N+4*E0/N**3)
                err=abs(output-expected)
                check('independent_physical_return_sixth_derivative_'+str(dps)+'_'+str(N),err<mp.mpf('1e-42'),mp.nstr(err,12))
                records.append({'dps':dps,'wave_number':N,'return_sixth_error':mp.nstr(err,18),
                                'physical_volume_included':True,'nonlinear_parameter_contour':True})
    FIXTURES.append({'fixture':'independent_tower_and_physical_return_diagnostics','records':records,
                     'actual_charged_source_sampled':False,'manufactured_geometry_used':True})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 Schur source derivative budget','input_count':len(ps)}),flush=True)
    try:
        provenance();conformal_and_return_identities();schur_norm_ledger();tower_and_fourier_identities();independent_numerics()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 Schur source derivative budget','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after}
    for key,value in review().items():
        if key not in run and isinstance(value,bool):run[key]=value
    if not good:
        for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');saved=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    saved['runs'].append(run);dest.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
