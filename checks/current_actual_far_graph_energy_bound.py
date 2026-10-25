"""Original-cutoff distant return for the actual coupled linear solution.

Supporting exact identities and outward arithmetic for the ordinary proof.
No finite angular/thermal selection and no numerical solution is substituted.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import factorial, comb
import sys, json, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb, ctx
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
STEM='current-actual-far-graph-energy-bound'
OWN=['notes/'+STEM+'.md','checks/current_actual_far_graph_energy_bound.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/current-complete-weak-mass-residual.md',
 'notes/current-joint-adjoint-residual-attempt.md',
 'notes/current-separated-quantum-graph-bound.md',
 'notes/current-global-mass-conversion.md',
 'notes/current-global-return-energy-attempt.md',
 'notes/current-coupled-global-return-attempt.md',
 'notes/current-controlled-solution-audit-closeout.md',
 'evidence/current-global-mass-conversion.json',
 'evidence/current-reference-clock-adjoint-trial.json',
 'evidence/current-coupled-global-return-attempt.json',
 'evidence/current-complete-weak-mass-residual.json']
CHECKS=[]
def check(name, ok):
    assert name not in {x['name'] for x in CHECKS}
    CHECKS.append({'name':name,'passed':bool(ok)})
    if not ok:raise AssertionError(name)
def enc(v):
    mid,rad,exp=v.mid_rad_10exp(50)
    return {'ball':str(v),'decimal_enclosure':{'mid':str(mid),'rad':str(rad),'exponent':str(exp)}}
def dec(v):
    d=v['decimal_enclosure'];return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])
def prior(name):
    r=json.loads((ROOT/'evidence'/(name+'.json')).read_bytes())['runs'][-1]
    assert r['status']=='passed'
    return r['result']
def execute():
    ctx.prec=512
    par=prior('current-global-mass-conversion')['parameters']
    beta,u,eps,tau,g,Fl=[dec(par[k]) for k in ['beta','u','eps','tau','g','Fl']]
    tr=prior('current-reference-clock-adjoint-trial')
    ret=prior('current-coupled-global-return-attempt')
    old=prior('current-complete-weak-mass-residual')
    gamma=dec(ret['gamma']);zh=1/(1+dec(ret['h']));C=1/u
    etaG=C/dec(tr['B_squared_over_G']);pi=arb.pi()
    check('same_actual_positive_bolt',gamma>0 and gamma<arb('1e-30') and zh>arb('.99') and zh<1)
    check('same_original_mass_certificate',dec(old['full_coupled_linear_mass_upper'])<arb('1.627532e-15'))
    qt=[dec(v) for v in old['actual_bolt_polynomial_coefficients']]
    check('bolt_coefficients_cover_full_polynomial',qt[0]>arb('.999999') and qt[0]<arb('1.000001')
          and qt[1]>arb('-1.334') and qt[1]<arb('-1.333') and qt[2]>arb('.5') and qt[2]<arb('.501')
          and abs(qt[3])<arb('1e-27') and abs(qt[4])<arb('1e-27'))
    # P(delta)=a+b delta+c delta^2+... is decreasing on [0,1],
    # using its coefficient enclosure; endpoints then give these reserves.
    check('bolt_polynomial_P_decreasing',arb('-1.333')+2*arb('.501')+7*arb('1e-27')<0)
    check('whole_interval_P_lower',arb('.999999')-arb('1.334')+arb('.5')-2*arb('1e-27')>arb('.16'))
    check('whole_interval_P_upper',arb('1.000001')+2*arb('1e-27')<arb('1.001'))
    check('lower_Q_includes_thermal_and_quadratic',arb('.16')>arb(1)/8)
    # The old whole-interval Q<.17 bound is retained, not inferred from
    # the much looser delta upper polynomial at delta=1.
    q=[beta,arb(0),u,-4*beta-2*u-tau,3*beta+u+eps,-2*eps,eps]
    radius=arb('.11')
    complex_error=sum((abs(q[j]/beta)*radius**j for j in range(1,7)),arb(0))
    check('actual_observation_germ_polynomial_guard',complex_error<arb('.01'))
    check('observation_disk_stays_in_guard',arb('.1')+1/(1000*arb(18))<radius)
    check('inverse_sqrt_Q_guard',1/(beta*(1-complex_error)).sqrt()<4)
    check('inverse_Q_guard',1/(beta*(1-complex_error))<16)
    check('riccati_rhs_and_lipschitz',16*(1+arb('1.1')**2)<36 and 8*arb('1.1')<36 and 32*arb('1.1')<36)
    check('riccati_picard_radius',arb(36)/1000<arb('.1'))
    check('denominator_germ_reserve',2*arb('.036')+arb('.036')**2<arb('.21'))
    check('actual_kernel_cauchy_prefactor',4*arb('1.1')/arb('.79')**2*arb('.036').exp()<8 and 4*8*(4+16)<=640)
    Flmin=1-arb('1.179e-8')
    gradient_squared=arb('1e5')*beta/(arb('.99')*Flmin)
    check('actual_gradient_energy_conversion',Fl>Flmin and gradient_squared<130**2)
    anchor=arb('17.6')/arb('.1').sqrt()+130*(arb('.8')/arb('.69')**2).sqrt()
    check('anchor_bound',anchor<300)
    check('middle_value_bound',300+130*arb(8).sqrt()<700)
    check('middle_growth_rate',arb('3.001').sqrt()<2)
    check('middle_L1_coefficient',700/arb('.79').sqrt()<800)
    cold=300*arb(2).sqrt()+130*(2*arb(8).sqrt()+pi.sqrt())/arb(2).sqrt()
    check('cold_integrable_energy_majorant',cold<1110)
    gu=arb('1e-30')
    sg=(arb('.79')/(2*gu)).log()/arb('1.001')
    check('cold_weight_product_decreases_with_k',sg>1/arb(17))
    check('cold_weight_bound_monotone_in_gamma',18/arb('1.001')-arb('.5')>0)
    cold_weight=1110*17/gu.sqrt()*(2*gu/arb('.79'))**(18/arb('1.001'))
    check('entire_cold_weighted_integral',cold_weight<arb('1e-500'))
    check('whole_weighted_return_constant',arb(801)*18/17<900)
    a=arb('.05');b=arb('.1');c=arb('.2');d=arb(5).sqrt()*(c-b)
    cm=arb(2).sqrt()*d;cp=arb(10).sqrt()*d
    check('actual_separation_positive',d>0 and cm>0 and cp>0)
    Aj=[]
    for j in range(5):
        Aj.append(2560*etaG/pi*factorial(j)*1000**j*(1+gamma*cp)*sum(
            (comb(j,k)*factorial(j-k+3)*factorial(k)/(cm**(j-k+3)*cp**(k+1)) for k in range(j+1)),arb(0)))
    t=sp.symbols('t')
    cutoff=1-sp.Rational(factorial(21),factorial(10)**2)*sum(
        ((-1)**j*comb(10,j)*t**(11+j)/sp.Integer(11+j) for j in range(11)))
    coeff=sp.Poly(cutoff,t)
    bern=[sum(coeff.nth(j)*sp.Rational(comb(i,j),comb(21,j)) for j in range(i+1)) for i in range(22)]
    check('original_beta_cutoff_bernstein_identity',bern==[sp.Integer(1)]*11+[sp.Integer(0)]*11)
    check('original_cutoff_endpoints',cutoff.subs(t,0)==1 and cutoff.subs(t,1)==0)
    check('cutoff_endpoint_jets_retained',all(sp.diff(cutoff,t,j).subs(t,0)==0 and sp.diff(cutoff,t,j).subs(t,1)==0 for j in range(1,11)))
    hh=[arb(1)]+[arb(factorial(21)//factorial(21-j))*2**(j-1)/arb(2).log()**j for j in range(1,5)]
    def recurrence(alpha):
        rows=[[arb(1)]]
        for m in range(3):
            old=rows[-1]
            rows.append([(alpha+j)*(old[j] if j<len(old) else 0)+(old[j-1] if j else 0) for j in range(m+2)])
        return rows
    # Verify the derivative recurrence as a formal identity on arbitrary f jets.
    z,al=sp.symbols('z al');jets=sp.symbols('f0:5');expr=jets[0]
    symbolic=[[sp.Integer(1)]]
    for k in range(1,4):
        expr=sp.expand(al*expr+z*sp.diff(expr,z)+sum(z*sp.diff(expr,jets[j])*jets[j+1] for j in range(4)))
        prev=symbolic[-1]
        row=[(al+j)*(prev[j] if j<len(prev) else 0)+(prev[j-1] if j else 0) for j in range(k+1)]
        symbolic.append(row)
        check('Euler_derivative_identity_'+str(k),sp.expand(expr-sum(row[j]*z**j*jets[j] for j in range(k+1)))==0)
    r1=recurrence(arb('.995'));r0=recurrence(arb('-.005'))
    L1=[];L0=[]
    for k in range(4):
        one=arb(0);zero=arb(0)
        for l in range(k+1):
            for j in range(k-l+1):
                one+=comb(k,l)*hh[l]*abs(r1[k-l][j])*b**(arb('.995')+j)*Aj[j+1]/(2*(arb('.995')+j)).sqrt()
                power=arb('-.005')+j
                maxpow=a**power if j==0 else b**power
                zero+=comb(k,l)*hh[l+1]*abs(r0[k-l][j])*maxpow*Aj[j]*arb(2).log().sqrt()
        L1.append(one);L0.append(zero)
    nf=lambda vals:sum((comb(3,k)*vals[k]**2 for k in range(4)),arb(0)).sqrt()
    Cfar=900*nf([L1[k]+L0[k] for k in range(4)])
    actual=Cfar*arb('1.1e-10');row=900*arb('1.1e-10')*Aj[0]
    quartic=20*beta*(Fl+23*g/(180*pi))
    check('complete_actual_far_graph_per_energy',Cfar<arb('1.124e-10'))
    check('complete_actual_linear_far_graph_target',actual<arb('1.24e-20'))
    check('original_shape_finite_part_retained',row<arb('1.27e-34'))
    check('original_scalar_quartic_row_coefficient',quartic>arb('3.3333332'))
    check('far_row_only_algebraic_quartic_allowance',row/quartic<arb('3.82e-35'))
    remaining=arb('1e-11')-2174*arb('1.627532e-15')-arb('1.24e-20')
    check('remaining_near_allowance_is_positive',remaining>arb('6.4617454e-12'))
    check('far_cost_not_decision_dominant',actual/remaining<arb('2e-9'))
    aux=arb('1.02e-46')+arb('6.50e-39')*arb('1.627532e-15')
    check('auxiliary_Mellin_tail_after_actual_mass_bound',aux<arb('1.020000106e-46'))
    return {'precision_bits':ctx.prec,'complex_actual_polynomial_error':enc(complex_error),
        'gradient_energy_constant_squared':enc(gradient_squared),'anchor_upper':enc(anchor),
        'cold_unweighted_coefficient':enc(cold),'cold_weighted_reserve_test':enc(cold_weight),
        'whole_input_weighted_L1_constant':'900/k for k>=18',
        'actual_kernel_derivative_coefficients':[enc(v) for v in Aj],
        'original_cutoff_derivative_guards':[enc(v) for v in hh],
        'gradient_output_coefficient_per_energy':enc(900*nf(L1)),
        'original_row_transition_coefficient_per_energy':enc(900*nf(L0)),
        'complete_far_graph_coefficient_per_energy':enc(Cfar),
        'actual_linear_far_graph_upper':enc(actual),'original_shape_row_finite_part_upper':enc(row),
        'far_row_only_quartic_algebraic_allowance':enc(row/quartic),
        'remaining_near_graph_allowance':enc(remaining),
        'remaining_allowance_scope':'Sufficient shape graph only; original scalar graph remains independent and unbounded here.',
        'auxiliary_Mellin_tail_upper':enc(aux),
        'scope':{'actual_linear_far_graph_bounded':True,'original_shape_finite_part_retained':True,
          'normal_far_cancels_commutator_exactly':True,'complete_angular_and_thermal_sums':True,
          'common_regulator_and_actual_domain_retained':True,'actual_near_graph_bounded':False,
          'actual_kappa_m_computed':False,'independent_phase_columns_bounded':False,
          'complete_coupled_graph_certified':False,'nonlinear_admission':False,'nonlinear_contraction':False,
          'exact_compatibility':False,'physical_matching':False,'validity_noise':False,
          'all_sector_stability':False,'actual_controlled_wormhole':False,
          'construction_obstruction':False,'full_objective_complete':False}}

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,
      'input_sha256_before':hashes,'frozen_sources_utf8':{f:(ROOT/f).read_text(encoding='utf-8') for f in OWN},
      'runtime':{'python':sys.version,'executable_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
                 'flint_version':flint.__version__,'sympy_version':sp.__version__}}
    try:
        run['result']=execute()
        run['input_sha256_after']={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after'])
        run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),
        'result':run.get('result'),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
