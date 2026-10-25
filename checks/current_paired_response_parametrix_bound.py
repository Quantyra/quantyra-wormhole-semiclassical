"""Verify scalar comparison constants; no renormalized coupled application."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
import json, sys, traceback
import sympy as sp
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb,ctx
ctx.prec=512
ROOT=Path(__file__).resolve().parents[1]
STEM='current-paired-response-parametrix-bound'
OWN=['notes/'+STEM+'.md','checks/current_paired_response_parametrix_bound.py',
     'evidence/'+STEM+'-review.json']
sys.path.insert(0,str(ROOT/'checks'))
saved_argv=sys.argv[:];sys.argv=[sys.argv[0]]
import current_physical_york_dual_assembly as york
sys.argv=saved_argv
CHECKS=[]

def check(name,value,detail=None):
    CHECKS.append(dict(name=name,passed=bool(value)))
    if not value:raise AssertionError(name)

york.check=check;york.prior.check=check

def enc(v):
    m,r,e=v.mid_rad_10exp(60)
    return dict(ball=v.str(40),decimal_enclosure=dict(mid=str(m),rad=str(r),exponent=str(e)))

def execute():
    par,horizon,zh,env=york.actual_parameters()
    co={k:env[v] for k,v in zip(york.powers,york.qs)}
    def bern(co,c):
        n=max(co)
        return [sum((co.get(i,arb(0))*c**i*comb(k,i)/comb(n,i) for i in range(k+1)),arb(0)) for k in range(n+1)]
    def deriv(j):
        return {k-j:v*factorial(k)/factorial(k-j) for k,v in co.items() if k>=j}
    check('global_actual_Q_upper',all(v<arb('.17') for v in bern(co,arb(1))))
    check('original_collar_Q_lower',all(v>arb('.16') for v in bern(co,arb('.2'))))
    geometry={}
    for end,ql,qp,qpp in [('.4','.13','.2','.75'),('.6','.08','.3','.9')]:
        a=arb(end)
        check('whole_interval_Q_lower_'+end,all(v>arb(ql) for v in bern(co,a)))
        check('whole_interval_Qprime_'+end,all(abs(v)<arb(qp) for v in bern(deriv(1),a)))
        check('whole_interval_Qsecond_'+end,all(abs(v)<arb(qpp) for v in bern(deriv(2),a)))
        geometry[end]=[ql,qp,qpp]
    check('whole_interval_Qthird_near',all(abs(v)<arb('4.1') for v in bern(deriv(3),arb('.4'))))
    check('exact_geodesic_seam_Qprime',co.get(1,arb(0))==0)
    M1=((par['q']+1)/par['C']).sqrt()
    check('entire_angular_mass_gap',M1>18)
    x,z=sp.symbols('x z',positive=True);M,p=sp.symbols('M p',positive=True)
    W=sp.Function('W')(x);aa=sp.Function('a')(x);bb=sp.Function('b')(x);tt=sp.Function('t')(x)
    seed={sp.diff(tt,x):-2*W*bb,sp.diff(aa,x):2*M*bb,sp.diff(bb,x):2*M*aa-2*W*tt}
    k2=M*M+W*W
    bpp=sp.diff(seed[sp.diff(bb,x)],x).subs(seed)
    check('exact_actual_coefficient_equation',sp.simplify(-bpp+4*k2*bb-2*sp.diff(W,x)*tt)==0)
    d_equation=sp.diff(2*sp.diff(W,x)*tt,x).subs(seed)-sp.diff(4*k2,x)*bb
    check('exact_actual_coefficient_derivative_equation',sp.simplify(d_equation-2*sp.diff(W,x,2)*tt+12*W*sp.diff(W,x)*bb)==0)
    caux=W*aa+M*tt
    check('exact_caux_derivative_cancellation',sp.simplify(sp.diff(caux,x).subs(seed)-sp.diff(W,x)*aa)==0)
    check('exact_caux_constraint',sp.simplify(caux**2+seed[sp.diff(bb,x)]**2/4-k2*(aa**2+tt**2))==0)
    check('exact_a_reconstruction',sp.simplify(W*caux/k2+M*seed[sp.diff(bb,x)]/(2*k2)-aa)==0)
    Q=sp.Function('Q')(z);f=sp.sqrt(Q);Wz=p/f
    dx=lambda v:f*sp.diff(v,z)
    w2=sp.diff(Q,z)**2/(2*Q)-sp.diff(Q,z,2)/2
    w3=5*sp.diff(Q,z)*sp.diff(Q,z,2)/(4*f)-3*sp.diff(Q,z)**3/(4*Q**sp.Rational(3,2))-f*sp.diff(Q,z,3)/2
    check('exact_frequency_second_jet',sp.simplify(dx(dx(Wz))/Wz-w2)==0)
    check('exact_frequency_third_jet',sp.simplify(dx(dx(dx(Wz)))/Wz-w3)==0)
    check('near_Wprime_bound',arb('.2')/(2*arb('.13').sqrt())<arb('.3'))
    check('near_Wsecond_bound',arb('.2')**2/(2*arb('.13'))+arb('.75')/2<arb('.55'))
    check('buffer_Wprime_bound',arb('.3')/(2*arb('.08').sqrt())<arb('.6'))
    check('buffer_Wsecond_bound',arb('.3')**2/(2*arb('.08'))+arb('.9')/2<arb('1.1'))
    check('near_Wthird_bound',5*arb('.2')*arb('.75')/(4*arb('.13').sqrt())+3*arb('.2')**3/(4*arb('.13')**arb('1.5'))+arb('.17').sqrt()*arb('4.1')/2<2)
    sigma=arb('1.05').sqrt();d18=18-sigma/2;rd=1-sigma/36
    delta=arb('.2')/arb('.17').sqrt()
    check('actual_endpoint_derivative_guard',2+4/rd<arb('6.12'))
    check('spectral_envelopes_decrease_above_18',18*delta>4)
    tails={}
    for tag,qmin,c1,c2,B,D in [('near','.13','.3','.55','.201','1.7'),('buffer','.08','.6','1.1','.401','4')]:
        r=(arb('.17')/arb(qmin)).sqrt();c1=arb(c1);c2=arb(c2)
        check(tag+'_barrier_weight_second_jet',c2+3*c1*c1<3)
        bt=r**3*18**2/d18*(-36*delta).exp()
        dt=2*r**3*(arb('6.12')+8*c1*r*r/d18)*18**2*(-18*delta).exp()
        check(tag+'_actual_b_return_guard',2*c1/3+bt<arb(B))
        check(tag+'_actual_bprime_return_guard',(2*c2+8*c1*c1)/3+dt<arb(D))
        tails[tag]=dict(b_tail=enc(bt),bprime_tail=enc(dt))
    ec=arb('.201')**2+arb('1.7')**2/(4*18**2)
    check('caux_two_power_guard',ec<arb('.043'))
    check('a_two_power_guard',ec+arb('1.7')/2<arb('.9'))
    check('P_two_power_guard',arb('.9')/(1-arb('.9')/324)<arb('.904'))
    check('V_two_power_guard',(arb('.043')+2*arb('.9')+arb('.9')**2/324)/(1-arb('.9')/324)**2<arb('1.86'))
    check('V_reverse_two_power_guard',(arb('.043')+2*arb('.9')+arb('.9')**2/324)/(1-arb('.043')/324)<arb('1.85'))
    check('bulk_tprime_guard',2*sp.Rational('.401')==sp.Rational('.802'))
    check('bulk_tsecond_guard',2*sp.Rational('.3')*sp.Rational('.401')+2*4==sp.Rational('8.2406'))
    b0forcing=(arb(2)+2*arb('.55')*arb('.802')+arb('.3')*arb('8.2406')+2*(arb('.55')+arb('.3')*arb('.802'))*arb('.6')+arb('.3')*arb('2.1'))/2
    check('actual_b0_second_forcing',b0forcing<arb('3.5'))
    check('w4_barrier_second_jet',3*arb('.55')+14*arb('.3')**2<3)
    rr=(arb('.17')/arb('.13')).sqrt()
    etail=rr**5*(arb(18)**4/d18+arb('.15')*18**2)*(-36*delta).exp()
    check('actual_b_minus_b0_endpoint_guard',arb('1.2')+etail<arb('1.21'))
    loga=(2*arb('1.21')+arb('.3')*arb('1.7')/2)/(1-arb('.9')/324)
    logv=2*arb('2.69')+arb('.3')*(arb('.9')+arb('.043'))/(1-arb('.043')/324)
    check('relative_first_P_jet_gain_two',loga<arb('2.69'))
    check('relative_first_V_jet_gain_two',logv<arb('5.67'))
    b0=sp.diff(W,x)*tt/(2*k2)
    ra=2*M*bb/aa;ra0=M*M*sp.diff(W,x)/(W*k2)
    cancel=2*M*(bb-b0)/aa-M*(sp.diff(W,x)/W)*seed[sp.diff(bb,x)]/(2*k2*aa)
    check('exact_log_a_derivative_cancellation',sp.simplify(ra-ra0-cancel)==0)
    t0=M/sp.sqrt(k2)
    tdiffprime=sp.diff(tt-t0,x).subs(seed)
    check('exact_tadpole_derivative_cancellation',sp.simplify(tdiffprime+W*sp.diff(W,x)/k2*(tt-t0)+2*W*(bb-b0))==0)
    check('tadpole_value_gain_two',arb('.043')+arb('1.7')/2<arb('.9'))
    check('tadpole_first_gain_two',sp.Rational('.3')*sp.Rational('.9')+2*sp.Rational('1.21')==sp.Rational('2.69'))
    c0=arb(2).sqrt()/(4*rd)
    invdiff=c0*arb('1.85')/(1-arb('1.85')/324)
    response=4*invdiff+arb('.9')
    check('global_equivalent_forms_delta',arb('1.85')/324<1)
    check('global_inverse_difference_guard',invdiff<arb('.678'))
    check('complete_mode_response_difference_guard',response<arb('3.609'))
    K1=21/arb(2).log();K2=840/arb(2).log()**2
    chi1=arb('.17').sqrt()*K1/arb('.05')
    chi2=arb('.17')*(K2+K1)/arb('.05')**2+arb('.1')*K1/arb('.05')
    check('original_cutoff_first_proper_derivative',chi1<250)
    check('original_cutoff_second_proper_derivative',chi2<122000)
    check('base_conormal_field_first_derivative',1+arb('.1')*arb('.2')/(2*arb('.16'))<arb('1.1'))
    check('base_conormal_field_second_derivative',arb('.17').sqrt()/2*((arb('.2')+arb('.1')*arb('.75'))/arb('.16')+arb('.1')*arb('.2')**2/arb('.16')**2)<arb('.4'))
    check('conormal_h_first_guard',chi1*arb('.25')+arb('1.1')<arb('63.6'))
    check('conormal_h_second_guard',chi2*arb('.25')+2*chi1*arb('1.1')+arb('.4')<31100)
    hh=sp.Function('h')(x);PP=sp.Function('P')(x);VV=sp.Function('V')(x);uu=sp.Function('u')(x)
    X=lambda v:hh*sp.diff(v,x)
    LL=lambda v:-sp.diff(PP*sp.diff(v,x),x)+VV*v
    comm=-sp.diff((X(PP)-2*PP*sp.diff(hh,x))*sp.diff(uu,x),x)-PP*sp.diff(hh,x,2)*sp.diff(uu,x)+X(VV)*uu
    check('exact_conormal_commutator_including_first_order',sp.simplify(X(LL(uu))-LL(X(uu))-comm)==0)
    ar=2*arb('.201')/(1-arb('.9')/324)
    check('actual_logP_jet_bound',arb('.3')+ar<arb('.704'))
    check('actual_logV_jet_bound',arb('.3')+2*ar<arb('1.107'))
    check('actual_gradient_to_value_form_ratio',((1+arb('.9')/324)/(4*(1-arb('.043')/324))).sqrt()<arb('.502'))
    check('difference_Pprime_bound',(1+arb('.9')/324)*arb('2.69')+arb('.704')*arb('.9')<arb('3.332'))
    check('difference_Vprime_bound',(1+arb('1.85')/324)*arb('5.67')+arb('1.107')*arb('1.85')<arb('7.752'))
    check('first_order_difference_factor',sp.Rational('.9')*sp.Rational('.502')/18==sp.Rational('.0251'))
    delta_comm=arb('.25')*(arb('2.75')*250+arb('11.084'))+arb('1.8')*arb('63.6')+arb('.0251')*31100
    actual_comm=arb('.25')*arb('.704')+2*arb('63.6')+arb('.25')*arb('1.107')+31100*arb('.502')/18
    check('complete_difference_commutator_form_guard',delta_comm<1070)
    check('complete_actual_commutator_form_guard',actual_comm<996)
    invcomm=c0*(2*996*arb('1.85')+1070)/(1-arb('1.85')/324)**2
    hcomm=4*invcomm+arb('.25')*(250*arb('.9')+arb('2.69'))
    check('global_resolvent_commutator_guard',invcomm<1751)
    check('complete_mode_response_commutator_guard',hcomm<7059)
    A0,B0,A1,B1=map(arb,['2.752','.85','7004.311','53.794'])
    check('retain_mass_numerator_value',4*invdiff+arb('.043')<A0)
    check('retain_mass_numerator_commutator',4*invcomm+arb('.25')*(250*arb('.043')+arb('.0129'))<A1)
    check('retain_mass_numerator_contact_derivative',arb('.25')*(250*arb('.85')+arb('2.675'))<B1)
    check('exact_tadpole_mass_numerator_derivatives',sp.Rational('.3')*sp.Rational('.043')==sp.Rational('.0129') and sp.Rational('.3')*sp.Rational('.85')+2*sp.Rational('1.21')==sp.Rational('2.675'))
    s=sp.symbols('s',positive=True)
    check('complete_thermal_integral_three_halves',sp.integrate((1+s*s)**(-sp.Rational(3,2)),(s,0,sp.oo))==1)
    check('complete_thermal_integral_one',sp.integrate(1/(1+s*s),(s,0,sp.oo))==sp.pi/2)
    check('decreasing_thermal_envelopes',sp.simplify(sp.diff((1+s*s)**(-sp.Rational(3,2)),s)+3*s/(1+s*s)**sp.Rational(5,2))==0 and sp.simplify(sp.diff(1/(1+s*s),s)+2*s/(1+s*s)**2)==0)
    check('actual_positive_temperature_guard',0<par['gamma']<arb('1e-33'))
    thermal0=arb('.17').sqrt()*(A0+arb.pi()*B0/2)+arb('1e-33')*(A0+B0)/18
    thermal1=arb('.17').sqrt()*(A1+arb.pi()*B1/2)+arb('1e-33')*(A1+B1)/18
    check('complete_thermal_value_guard',thermal0<arb('1.686'))
    check('complete_thermal_commutator_guard',thermal1<2923)
    orig=2*arb('2.5').sqrt()*arb('3.609')*arb('35.534')
    check('original_energy_and_action_conversion',orig<406)
    origthermal=2*arb('2.5').sqrt()*arb('35.534')*arb('1.686')
    origthermalcomm=2*arb('2.5').sqrt()*arb('35.534')*(2923+arb('.0625')*arb('1.686'))
    check('original_complete_thermal_source_guard',origthermal<190)
    check('original_complete_thermal_commutator_density_guard',origthermalcomm<329000)
    n,N,q=sp.symbols('n N q',integer=True,positive=True)
    angular=sp.summation(2*(q+2*n),(n,1,N))
    check('complete_angular_upper_estimator_count',sp.expand(angular-2*q*N-2*N*(N+1))==0)
    check('remaining_absolute_angular_estimator_diverges',sp.limit(angular,N,sp.oo)==sp.oo)
    freq=sp.symbols('freq',positive=True)
    check('complete_thermal_estimator_harmonic_tail',sp.limit(freq/sp.sqrt(M*M+freq*freq/sp.Rational('.17')),freq,sp.oo)==sp.sqrt(sp.Rational('.17')))
    remain=arb('5.34235e-12')-arb('5.042e-23')
    ge=arb('.0012069')/(arb('.0195239')+arb('4.583e-13'))
    check('conditional_affine_allowance_unchanged',remain>arb('5.34234999994957e-12'))
    check('conditional_energy_feedback_ceiling_unchanged',arb('.06181654279993')<ge<arb('.06181654279994'))
    check('conditional_fast_screen_scope_unchanged',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    return dict(outcome='unresolved_complete_coupled_return_with_two_power_response_and_first_commutator',
      actual_M1=enc(M1),geometry_guards=geometry,actual_endpoint_tail_expressions=tails,b0_endpoint_tail=enc(etail),
      relative_log_P_derivative=enc(loga),relative_log_V_derivative=enc(logv),global_inverse_difference=enc(invdiff),
      complete_mode_response_difference=enc(response),original_cutoff_derivative_expressions=[enc(chi1),enc(chi2)],
      global_resolvent_commutator=enc(invcomm),complete_mode_response_commutator=enc(hcomm),original_density_per_energy_coefficient=enc(orig),
      complete_thermal_response=enc(thermal0),complete_thermal_commutator=enc(thermal1),
      original_thermal_source_coefficient=enc(origthermal),original_thermal_commutator_coefficient=enc(origthermalcomm),
      upward_guards={'response_difference':'3.609/K','first_conormal_commutator':'7059/K','original_near_source_difference':'406 G gamma d_n/(pi B^2 K) per E'},
      actual_paired_and_bolt_realization_retained=True,all_angular_and_thermal_labels=True,finite_mode_selection=False,
      complete_relative_thermal_sum_certified=True,common_regulator_angular_sum_executed=False,common_subtracted_strong_sum_certified=False,auxiliary_comparison_is_physical_quantum_action=False,
      raw_absolute_sum_estimator='coarse 1/K thermal failure repaired by retaining mass ratios; remaining angular upper bound grows as 2qN+2N(N+1), not a physical divergence theorem',
      actual_physical_operator_applications=0,coupled_signed_solve_executed=False,kappa_m=None,actual_mass=None,other_coupled_columns_certified=False,
      conditional_remaining_affine=enc(remain),conditional_energy_feedback_ceiling=enc(ge),actual_energy_feedback_comparison=None,
      physical_norm_comparisons_closed=False,nonlinear_requirement_closed=False,full_coupled_feasibility_decided=False,
      construction_obstruction=False,full_objective_preserved=True,closed_audit_preserved=True)

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    out=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    review=json.loads((ROOT/OWN[2]).read_bytes())
    paths=set(OWN+review['retained_inputs'])
    for module in list(sys.modules.values()):
        f=getattr(module,'__file__',None)
        if f:
            try:rel=Path(f).resolve().relative_to(ROOT).as_posix()
            except ValueError:continue
            if rel.endswith('.py'):paths.add(rel)
    digest=lambda p:sha256(p.read_bytes()).hexdigest()
    before={p:digest(ROOT/p) for p in sorted(paths)}
    fr=Path(flint.__file__).resolve().parent
    runtime=dict(python=sys.version,python_executable_sha256=digest(Path(sys.executable)),flint_version=flint.__version__,
      flint_artifacts={p.relative_to(fr).as_posix():digest(p) for p in sorted(fr.rglob('*')) if p.is_file() and p.suffix in ['.py','.pyd','.dll']})
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),mode=mode,input_sha256_before=before,precision_bits=512,runtime=runtime,
      frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN})
    try:
        check('ordinary_argument_frozen',before[OWN[0]]==review['ordinary_argument_sha256'])
        run['result']=execute()
        run['input_sha256_after']={p:digest(ROOT/p) for p in before}
        check('frozen_inputs_unchanged',before==run['input_sha256_after'])
        run['status']='passed'
    except Exception:run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    ledger['runs'].append(run);out.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write((json.dumps(run)+'\n').encode())
    print(json.dumps(dict(mode=mode,status=run['status'],checks=len(CHECKS),result=run.get('result'),traceback=run.get('traceback')),indent=2))
    if run['status']!='passed':raise SystemExit(1)
