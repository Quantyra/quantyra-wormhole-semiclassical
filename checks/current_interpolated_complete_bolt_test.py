"""Complete bolt estimator sensitivity, not a numerical inverse application."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction
import json, sys, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb, ctx
ROOT=Path(__file__).resolve().parents[1]
STEM='current-interpolated-complete-bolt-test'
OWN=['notes/'+STEM+'.md','checks/current_interpolated_complete_bolt_test.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+['notes/current-coupled-truncation-energy-error.md',
 'notes/current-complete-bolt-energy-bound.md','evidence/current-complete-bolt-energy-bound.json',
 'notes/current-coupled-application-residual-bound.md','evidence/current-coupled-application-residual-bound.json',
 'evidence/current-localized-bolt-adjoint-sensitivity.json',
 'notes/current-strong-quotient-reconstruction-test.md',
 'notes/current-implicit-phase-graph-domain.md']
CHECKS=[]
def sh(p):return sha256(p.read_bytes()).hexdigest()
def check(name,ok):
    CHECKS.append({'name':name,'passed':bool(ok)})
    print(name+(': passed' if ok else ': FAILED'),flush=True)
    if not ok:raise AssertionError(name)
def enc(x):
    m,r,e=x.mid_rad_10exp(50)
    return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
    if isinstance(x,dict):
        e=x['decimal_enclosure'];return arb(e['mid']+' +/- '+e['rad'])*arb(10)**int(e['exponent'])
    return arb(x.strip('[]'))
def prior(stem):
    r=json.loads((ROOT/'evidence'/f'{stem}.json').read_bytes())['runs'][-1]
    check('accepted_'+stem,r['status']=='passed')
    for p,h in r['input_sha256_before'].items():
        if p.startswith('notes/') and p in INPUTS:
            check('unchanged_'+p,sh(ROOT/p)==h)
    return r['result']
def aa(x):return arb(x.numerator)/arb(x.denominator)
def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('ordinary_argument_checker_frozen',all(sh(ROOT/p)==h for p,h in review['scientific_sha256'].items()))
    prev=prior('current-complete-bolt-energy-bound')
    appprev=prior('current-coupled-application-residual-bound')
    pars=prior('current-localized-bolt-adjoint-sensitivity')['parameters']
    h,gamma,beta,Cgeom,B2=[dec(pars[k]) for k in ['h','gamma','beta','C','B2']]
    Fl=1-384*arb.pi()*arb('1e22')*Cgeom*beta/B2
    CL=(500/(arb('.99')*Fl)).sqrt()
    CD=(100*arb.pi()*B2*beta/Cgeom+2*arb(500)**4*CL**2).sqrt()
    rstar=CL*gamma/(4*CD);Cnew=4*(CL*CD/gamma).sqrt()
    check('averaging_interval_inside_actual_half_domain',rstar>0 and rstar<arb('.99'))
    avg=CL/rstar.sqrt()+4*CD*rstar.sqrt()/gamma
    check('optimized_averaging_formula',avg.overlaps(Cnew))
    check('averaging_terms_balance',(CL/rstar.sqrt()).overlaps(4*CD*rstar.sqrt()/gamma))
    Cold=arb('7.477e50')
    check('strict_shape_comparison_improvement',Cnew<Cold)
    const=prev['constants'];Cs=dec(const['scalar_point_Cs']);Kxi=dec(const['Kxi'])
    J1=dec(const['J1']);ellb=(1+h)/h;q=arb(169)/500;f=q/(1+q)
    J0=ellb*(Kxi*Cnew/2+Cs*q.sqrt());Cb3=J0*f**3+J1*f**2
    Knew=ellb*(Cs**2+(Kxi*Cnew/2)**2).sqrt()
    budget=arb('1.667e-9');half=budget/2;energy=arb('1.1e-10')
    e1=arb('9.204e-79')*arb('.3961')**128
    e2=arb('9.204e-77')*arb('.3961')**128
    def app(C):
        eps=e2*C*C/beta;ef=e1*C/beta
        K=ellb*(Cs*Cs+(Kxi*C/2)**2).sqrt()
        return K*(ef+eps*energy)/(1-eps)
    eps=e2*Cnew*Cnew/beta;ef=e1*Cnew/beta
    check('complete_retained_form_inverse_perturbation_admissible',eps<1)
    Us=Knew*ef/(1-eps);Uh=Knew*eps*energy/(1-eps);U=app(Cnew)
    check('both_source_and_Hessian_errors_retained',U.overlaps(Us+Uh))
    check('old_bound_reproduces_excessive_estimator',app(Cold)>half)
    check('new_bound_closes_complete_transport_to_bolt_slice',U<half)
    # Monotone crossing, not fitting the estimate to a result.
    lo=Fraction(1);hi=Fraction(10**40)
    check('crossing_bracket_has_invertible_operator',e2*aa(hi)**2/beta<1)
    check('crossing_bracket_straddles_half_margin',app(aa(lo))<half and app(aa(hi))>half)
    for _ in range(80):
        mid=(lo+hi)/2;value=app(aa(mid))
        if value<half:lo=mid
        elif value>half:hi=mid
        else:raise ArithmeticError('Crossing enclosure requires a wider certified bracket')
    check('certified_crossing_endpoints',app(aa(lo))<half and app(aa(hi))>half)
    check('attained_gain_exceeds_required_gain',Cnew<aa(lo))
    # Every polynomial denominator is rational and exact before conversion.
    xr=1+Fraction(2)/Fraction(169,500);ds=[Fraction(1),xr]
    for k in range(2,98):ds.append(2*xr*ds[-1]-ds[-2])
    rows=[]
    for N in range(3,100):
        D=aa(ds[N-3]);h2=(1+q)**2/D;h3=(1+q)**3/D
        rows.append((N,h2,h3,energy*h3*Cb3,D/aa(ds[29])))
    full=next(r for r in rows if r[3]<budget)
    halfrow=next(r for r in rows if r[3]<half)
    for label,r,margin in [('full',full,budget),('half',halfrow,half)]:
        prevrow=next(v for v in rows if v[0]==r[0]-1)
        check('first_degree_'+label,prevrow[3]>margin and r[3]<margin)
    row32=next(r for r in rows if r[0]==32)
    check('P32_complete_bolt_upper_estimator_still_excessive',row32[3]>budget)
    check('scalar_trace_now_dominates_first_output_term',Cs*q.sqrt()>Kxi*Cnew/2)
    selected={}
    for N,h2,h3,tail,gain in rows:
        if N not in {32,full[0]-1,full[0],halfrow[0]-1,halfrow[0],76}:continue
        pref=energy*h3
        selected[str(N)]={k:enc(v) for k,v in {
          'h2':h2,'h3':h3,'bolt_tail_upper':tail,'denominator_gain_from32':gain,
          'conditional_defect_product_ceiling_application_zero':arb('.5')/h3,
          'conditional_graph_coefficient_ceiling_other_debits_zero':arb('5.34235e-12')/pref,
          'conditional_physical_coefficient_ceiling_other_debits_zero':arb('5e-9')/pref,
          'conditional_mass_product_ceiling_center_application_zero':arb('5e-5')/h2,
          'conditional_feedback_product_ceiling_other_debits_zero':arb('.0012069')/h2}.items()}
        selected[str(N)]['actual_coefficients_centers_and_application_errors']=None
    src=appprev['source_and_output_maps'];MRbeta=dec(src['physical_source_to_energy_MR_over_beta']);Krow=dec(src['Krow_global_row_recovery'])
    cf=Knew*MRbeta*Krow;cm=Knew*MRbeta/arb(5).sqrt()
    oldU=dec(appprev['actual_retained_operator_application_estimator']['whole_bolt_error_upper'])
    return {'precision_bits':512,'constants':{k:enc(v) for k,v in {
       'CL':CL,'CD':CD,'averaging_radius':rstar,'Cnew':Cnew,'point_bound_gain':Cold/Cnew,
       'required_point_gain_lower':Cold/aa(hi),'required_point_gain_upper':Cold/aa(lo),
       'point_half_margin_crossing_lower':aa(lo),'point_half_margin_crossing_upper':aa(hi),
       'J0new':J0,'J1_unchanged':J1,'scalar_point_Cs_unchanged':Cs,
       'complete_Cb3_upper':Cb3,'complete_unfiltered_Knew':Knew,
       'P32_bolt_upper':row32[3],'P32_remaining_required_gain':row32[3]/budget,
       'complete_Cb3_estimator_gain':dec(const['complete_Cb3_upper'])/Cb3}.items()},
      'complete_retained_inverse_transport_test':{k:enc(v) for k,v in {
       'epsilon_new':eps,'source_error_energy_dual':ef,'source_bolt_debit':Us,
       'Hessian_bolt_debit':Uh,'whole_bolt_debit':U,'fraction_half_bolt_margin':U/half,
       'slack_factor_half_margin':half/U,'improvement_old_complete_estimator':oldU/U}.items()},
      'conditional_residual_targets':{k:enc(v) for k,v in {
       'energy_dual_cap':half/Knew,'fast_residual_bolt_coefficient':cf,'finitepart_bolt_coefficient':cm,
       'global_fast_cap_moment_zero':half/cf,'finitepart_cap_fast_zero':half/cm}.items()},
      'degree_test':{'first_full_bolt_degree':full[0],'first_half_bolt_degree':halfrow[0],
       'rows':selected,'higher_degree_adopted':False},
      'scope':{'complete_two_row_bolt_map_bounded':True,'complete_transport_to_bolt_slice_passes':True,
       'P32_bolt_upper_estimator_fails':True,'numerical_inverse_applied':False,
       'polynomial_application_errors_certified':False,'actual_strong_residual_evaluated':False,
       'actual_quartic_finitepart_evaluated':False,'original_source_graph_physical_maps_closed':False,
       'original_defect_certified':False,'actual_mass_target_certified':False,
       'every_coupled_entry_certified':False,'complete_coupled_estimate_certified':False,
       'nonlinear_admission_or_contraction_closed':False,'physical_or_stability_gate_closed':False,
       'construction_obstruction':False,'normalization_changed':False,'quantum_fast_scope_unchanged':True}}
if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    path=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={p:sh(ROOT/p) for p in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
      'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
      'runtime':{'python':sys.version,'python_sha256':sh(Path(sys.executable)),'flint_version':flint.__version__}}
    try:
        run['result']=execute();run['input_sha256_after']={p:sh(ROOT/p) for p in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
    except Exception:run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(path.read_bytes()) if path.exists() else {'runs':[]};ledger['runs'].append(run)
    path.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    summary=json.dumps({'status':run['status'],'mode':mode,'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as stream:stream.write(summary.encode())
    print(summary,flush=True);sys.exit(0 if run['status']=='passed' else 1)
