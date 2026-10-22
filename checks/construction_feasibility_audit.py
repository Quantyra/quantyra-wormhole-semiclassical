"""Budget arithmetic and normal-model diagnostics, not continuum certification."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as F
from math import factorial
from itertools import product
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='construction-feasibility-audit'
PREVIOUS='finite-outer-shape-source'
OWN=['checks/construction_feasibility_audit.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
ROWS=[]; FIXTURES=[]
PROTOCOL={'digits':[60,90], 'identity_tolerance':'1e-45',
          'epsilon_c':['1','1e-4'], 'design_metric_radius':'1e-8',
          'gaussian_nodes':32,'packet_lengths':[4,8,16],
          'actual_quantum_state_sampled':False, 'continuum_proof_machine_verified':False}

def check(name,ok,detail=''):
    assert name not in {x['name'] for x in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row); print(json.dumps(row),flush=True)
    if not ok: raise AssertionError(name)

def saved(stem): return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def fixture(stem,name): return next(x for x in saved(stem)['fixtures'] if x['fixture']==name)
def review(): return json.loads((ROOT/OWN[2]).read_bytes())
def paths(): return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})
def mf(v): return mp.mpf(v.numerator)/v.denominator if isinstance(v,F) else mp.mpf(v)
def dec(v):
    with mp.workdps(40): return mp.nstr(mf(v),25)
def pair(v): return {'exact':str(v),'decimal':dec(v)}

def inputs_and_budgets():
    for stem in [PREVIOUS,'finite-interior-shape-source','finite-refreshed-local-mixed','finite-higher-york-source','finite-radius-operator-source']:
        old=saved(stem)
        check(stem+'_retained_success',old['status']=='passed' and all(r['passed'] for r in old['checks']))
        check(stem+'_unchanged_inputs',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    outer=fixture(PREVIOUS,'global_physical_shape_and_TT_source')
    low=fixture('finite-refreshed-local-mixed','low_TT_inverse_composition')
    york=fixture('finite-higher-york-source','density_conversion_and_complete_shape_budget')
    interior=fixture('finite-interior-shape-source','complete_interior_shape_budget')
    scalar=fixture('finite-radius-operator-source','complete_refined_scalar_and_anomaly')
    ci=F(interior['C_I_exact']); co=F(outer['C_O']); cb=F(outer['C_b']); cf=F(outer['C_F4'])
    cs=F(york['C_sigma_exact']); cy=F(york['C_b4_exact'])
    ct=F(scalar['C_tau_exact']); start=F(low['Cstart_exact']); hi=F(low['mhi_exact'])
    start_tight=start-F('1e20')+ct
    cft=cs+cy*cb; cfi=cs+cy*ci
    g0=F(1,2000**10); gl=(1-F('1e-34'))*g0; gu=g0
    check('physical_root_enclosure',F(11,1080)*g0<F('1e-34'))
    check('retained_benchmark_guards',F('9e-34')<gl<gu<F('1.6e-33') and F('3e-17')**2<gl and gu<F('4e-17')**2)
    check('complete_interior_and_outer_reconcile',ci==F(outer['C_I']) and cb==ci+co and cf==F('1e4')+F('1e15')*cb)
    check('exact_York_recombination',cs<F('1e4') and cy<F('1e15') and cft<cf and cft>F('1e201'))
    check('latest_scalar_in_start',0<ct<F('1e20') and start_tight>F('3e28') and start_tight<=start)
    check('interior_and_high_mixed_sizes',cfi>F('4e118') and hi>F('2e12'))
    check('known_low_composition_identity',start==F('1e20')+F('6e-12')*F(low['Klow_bar_exact']))
    lower_r=F('.97'); lower_d=F(2); records=[]
    for eps in map(F,PROTOCOL['epsilon_c']):
        for label,g in [('lower_g',gl),('upper_g',gu)]:
            bc=F('1e57')*g*eps/800
            # Unknown upper constants remain symbolic; no claimed measured norm.
            thresholds={
                'Bc':bc, 'R_if_high_term_zero':bc/start_tight,
                'each_Rs_if_other_two_at_point99_and_high_term_zero':bc/(F('.99')**2*start_tight),
                'Cstart_if_high_term_zero_and_R_lower':bc/lower_r,
                'signed_high_budget_at_R_lower':bc/lower_r-start_tight,
                'signed_F4_threshold_with_current_low_and_CD4_lower':(bc/lower_r-start_tight)/(hi*lower_d),
                'signed_CD4_threshold_with_current_low_and_tight_source':(bc/lower_r-start_tight)*g/(hi*cft),
                'signed_Cb_threshold_with_current_low_and_norms_lower':(g*(bc/lower_r-start_tight)/(hi*lower_d)-cs)/cy,
                'R_CD4_if_low_term_zero_published':bc*g/(hi*cf),
                'R_CD4_if_low_term_zero_tight':bc*g/(hi*cft),
                'R_CD4_if_low_term_zero_interior_only':bc*g/(hi*cfi),
                'C_F_if_low_term_zero_and_norms_lower':bc*g/(lower_r*hi*lower_d),
                'F4_if_low_term_zero_and_norms_lower':bc/(lower_r*hi*lower_d),
                'scalar_only_R_using_exact_C_tau':bc/ct,
                'old_rounded_equal_split_R':F('6.25e24')*g*eps,
                'old_rounded_equal_split_R_CD4_F4':F('6.25e40')*g*eps,
                'old_rounded_equal_split_R_CD4':F('6.25e-162')*g*g*eps,
                'beta_target':F('1e-49')*eps,
                'exact_transfer_conservative_target':F('1e57')*eps,
            }
            for name in ['published','tight','interior_only']:
                cap={'published':cf,'tight':cft,'interior_only':cfi}[name]
                thresholds['allowance_over_target_'+name]=lower_r*(start_tight+hi*lower_d*cap/g)/bc
            check(f'negative_remaining_high_budget_{eps}_{label}',thresholds['signed_high_budget_at_R_lower']<0)
            check(f'known_low_term_cannot_certify_{eps}_{label}',thresholds['R_if_high_term_zero']<lower_r)
            check(f'global_allowance_miss_{eps}_{label}',thresholds['allowance_over_target_tight']>F('1e225')/eps)
            check(f'interior_allowance_miss_{eps}_{label}',thresholds['allowance_over_target_interior_only']>F('1e142')/eps)
            records.append({'epsilon_c':str(eps),'g_endpoint':label,'g':pair(g),'thresholds':{k:pair(v) for k,v in thresholds.items()}})
    # Most generous broad g bound still proves these claims without decimals.
    broad=F('1.6e-33'); bc_broad=F('1e57')*broad/800
    check('robust_low_allowance_required_epsilon',lower_r*start_tight/bc_broad>F('1e7'))
    check('robust_R_CD4_threshold',bc_broad*broad/(hi*cft)<F('1.6e-225'))
    hypothetical=cs+cy*(ci+co/F('1e160'))
    check('outer_only_160_orders_cannot_remove_interior',hypothetical>=cfi)
    FIXTURES.append({'fixture':'complete_allowances_and_thresholds',
        'constants':{k:pair(v) for k,v in {'g_lower':gl,'g_upper':gu,'C_I':ci,'C_O':co,'C_b':cb,'C_F_published':cf,'C_F_tight':cft,'C_F_interior_only':cfi,'C_sigma':cs,'C_b4':cy,'C_tau':ct,'Cstart':start,'Cstart_tight':start_tight,'mhi':hi,'York_recombination_factor':cf/cft,'outer_over_interior':co/ci}.items()},
        'records':records,'counterfactual_outer_160_order_reduction':pair(hypothetical),
        'counterfactual_is_proved_source_bound':False,'actual_correction_lower_bound':False,
        'missing_upper_constants_assigned_values':False})
    return gl,gu

def exact_inequalities():
    r=F('1e-8'); y=r/2; z0=F(1,4); z2=1/(2*r)
    check('Newton_ball_allocation',y+z0*r+z2*r*r/2<=r and z0+z2*r==F(3,4))
    check('Newton_design_targets',y==F('5e-9') and z2==F('5e7'))
    check('design_interior_margin',F('.001')-F('1.3e-4')==F('8.7e-4') and r<F('8.7e-4'))
    check('design_cold_and_outer_thresholds',F('5e-7')-r==F('4.9e-7') and F('1e-7')-r==F('9e-8'))
    check('low_response_allowance_ratios',F('1.1e-7')/y==22 and F('2.2e-6')/y==440)
    a,b,M,L,e=S.symbols('a b M L e',positive=True)
    check('normalization_invariant_Newton_product',S.simplify((a*M/b)**2*(b*L/a**2)*(b*e)-M*M*L*e)==0)
    E=S.symbols('E',positive=True); r0,r1,r2,q3=S.symbols('r0 r1 r2 q3')
    check('three_transfer_error_weights',S.expand(r0+(r1+(r2+q3/E)/E)/E-r0-r1/E-r2/E**2-q3/E**3)==0)
    emin=F('180')/F('1.6e-33')
    error=(F('2e-49')+F('2e-14')/emin+F('2e21')/emin**2+F('2e56')/emin**3)/25
    check('approximate_transfer_residual_allocation',error<F('1e-49'))
    check('retained_exact_transfer_allocation',F('1e57')/(25*emin**3)<F('.281e-49'))
    transfer=(1+F('1e-32'))**-3
    check('positive_transfer_lower_product',transfer>F('.97'))
    check('normal_log_bound_exp40',sum(F(40)**j/factorial(j) for j in range(101))>1+2/F('3e-17'))
    check('normal_A_bound',F(1,144)+F(23,8640*9)*F('1.6e-33')<F('.007'))
    upper=9*(F('.007')*F('1.001')+F(54,480*9)*F('1.001')**2)
    check('normal_d_over_nu4_upper',upper<F(1,5),str(upper))
    check('physical_H4_coefficient_upper',sum(F(108)**j/F(100)**4 for j in range(5))<2)
    check('TT_inverse_lower_reserve',F(1,2)/upper>F(5,2)>2)
    check('compatibility_parameter_half_width',(F(7,25)-F(1,4))/2==F('.015'))
    FIXTURES.append({'fixture':'exact_certificate_and_lower_bound_arithmetic','design_radius':str(r),
                     'design_Y0':str(y),'design_Z0':str(z0),'design_Z2':str(z2),
                     'normal_d_upper_over_nu4':str(upper),'transfer_lower_product':str(transfer),
                     'actual_metric_ball_admitted':False,'beta_source_budget_is_metric_budget':False})

def tensor_polynomials():
    s,v=S.symbols('s v',real=True); tensors=[{():S.Integer(1)}]; polys=[]
    for j in range(5):
        tensor=tensors[j]
        norm=S.expand(sum(p.subs(s,S.Rational(3,2)+S.I*v)*p.subs(s,S.Rational(3,2)-S.I*v) for p in tensor.values()))
        polys.append(norm)
        check(f'horospherical_radial_component_{j}',S.expand(tensor[(0,)*j]-(-s)**j)==0)
        check(f'horospherical_norm_real_{j}',S.im(norm)==0)
        if j==4: break
        new={}
        for axis in range(4):
            for ids in product(range(4),repeat=j):
                if axis==0: value=-s*tensor[ids]
                else:
                    value=0
                    for slot,index in enumerate(ids):
                        if index==0:
                            repl=ids[:slot]+(axis,)+ids[slot+1:]; value-=tensor[repl]
                        elif index==axis:
                            repl=ids[:slot]+(0,)+ids[slot+1:]; value+=tensor[repl]
                new[(axis,)+ids]=S.expand(value)
        tensors.append(new)
    check('independent_physical_Hessian_norm',S.expand(polys[2]-((v*v+S.Rational(9,4))**2+3*(v*v+S.Rational(9,4))))==0)
    check('connection_slots_are_not_radial_only',S.expand(polys[4]-(v*v+S.Rational(9,4))**4)!=0)
    FIXTURES.append({'fixture':'full_physical_normal_tensor_polynomials','dimension':4,
                     'squared_tensor_amplitudes':[str(p) for p in polys],
                     'tensor_slot_counts':[len(t) for t in tensors],
                     'radial_derivative_proxy_used':False})
    return [S.lambdify(v,p,'mpmath') for p in polys]

def numerical_normal(gl,gu,polys):
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            for gval in [gl,gu]:
                g=mf(gval); nu=1/mp.sqrt(g); kappa=mp.mpf('10.8')
                s=mp.mpf('1.5')+1j*nu
                full=mp.digamma(s-1)+mp.pi*mp.cot(mp.pi*s/2)/2+kappa
                real=mp.re(mp.digamma(mp.mpf('.5')+1j*nu))-mp.pi/(2*mp.cosh(mp.pi*nu))+kappa
                check(f'digamma_critical_identity_{dps}_{gval}',abs(full-real)<mp.mpf('1e-45'))
                aa=1/(48*mp.pi*g)+23/(8640*mp.pi**2); rr=nu*nu+mp.mpf('2.25')
                dd=9*(aa*rr+rr*(rr-2)*real/(480*mp.pi**2))
                low=abs(s)**8/(dd*sum(mp.sqrt(p(nu)) for p in polys))
                check(f'normal_inverse_lower_diagnostic_{dps}_{gval}',dd>0 and real<54 and dd/nu**4<mp.mpf('.2') and low>mp.mpf('2.5'))
                fast=mp.sqrt(60*mp.pi+17*g/12)
                check(f'fast_frequency_above_Planck_{dps}_{gval}',13<fast<14)
                records.append({'dps':dps,'g':str(gval),'normal_test_nu':mp.nstr(nu,25),
                                'd_over_nu4':mp.nstr(dd/nu**4,25),
                                'physical_ratio_lower':mp.nstr(low,25),'fast_k_sqrtG':mp.nstr(fast,25),
                                'exact_Nq3_over_epsilon':mp.nstr(25*(60*mp.pi/g+mp.mpf(11)/3)**3*mp.mpf('1e-49'),25)})
            for vv in [mp.mpf(0),mp.mpf(1),mp.mpf(100)]:
                amp=[mp.sqrt(p(vv)) for p in polys]
                check(f'frame_generator_norm_bound_{dps}_{vv}',all(amp[j+1]<=mp.sqrt(vv*vv+mp.mpf('2.25')+3*j*j)*amp[j]+mp.mpf('1e-45') for j in range(4)))
            # Independently sum the convergent psi difference and bound the tail.
            vv=mp.mpf('1.7'); nmax=100
            term=lambda x:vv**2/(x*(x*x+vv*vv))
            finite=sum(term(n+mp.mpf('.5')) for n in range(nmax))
            exact=mp.re(mp.digamma(mp.mpf('.5')+1j*vv))-mp.digamma(mp.mpf('.5'))
            integ=mp.log(1+vv**2/(nmax+mp.mpf('.5'))**2)/2
            check(f'positive_digamma_series_integral_tail_{dps}',finite+integ<exact<finite+integ+term(nmax+mp.mpf('.5')))
    FIXTURES.append({'fixture':'independent_normal_symbol_and_validity','records':records,'physical_instability_inferred':False})

def normal_packet_model(gl,polys):
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            g=mf(gl); nodes,weights=mp.gauss_quadrature(PROTOCOL['gaussian_nodes'],'hermite')
            def symbol(v):
                rr=v*v+mp.mpf('2.25')
                hm=mp.re(mp.digamma(mp.mpf('.5')+1j*v))-mp.pi/(2*mp.cosh(mp.pi*v))+9
                return 9*((1/(48*mp.pi*g)+23/(8640*mp.pi**2))*rr+rr*(rr-2)*hm/(480*mp.pi**2))
            for center in [mp.mpf(10),mp.mpf(100)]:
                point=(center*center+mp.mpf('2.25'))**4/(symbol(center)*sum(mp.sqrt(p(center)) for p in polys))
                errors=[]
                for length in PROTOCOL['packet_lengths']:
                    frequencies=[center+z/length for z in nodes]
                    def integral(fun): return sum(w*fun(v) for w,v in zip(weights,frequencies))/mp.sqrt(mp.pi)
                    numerator=mp.sqrt(integral(lambda v:(v*v+mp.mpf('2.25'))**8))
                    denominator=sum(mp.sqrt(integral(lambda v,p=p:p(v)*symbol(v)**2)) for p in polys)
                    errors.append(abs(numerator/denominator/point-1))
                check(f'normal_packet_concentration_{dps}_{center}',errors[0]<mp.mpf('.02') and all(errors[j+1]<mp.mpf('.35')*errors[j] or errors[j+1]<mp.mpf('1e-45') for j in range(2)))
                records.append({'dps':dps,'center':str(center),'relative_errors':[mp.nstr(x,20) for x in errors]})
    FIXTURES.append({'fixture':'normal_Fourier_packet_diagnostic','records':records,
                     'actual_operator_packet_rate_claimed':False,'compact_packet_continuum_proof_machine_verified':False})

def noncommuting_positive_model():
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            W=mp.matrix([[2,mp.mpf('.2')],[mp.mpf('.2'),1]])
            Q=mp.matrix([[2,mp.mpf('.3')],[mp.mpf('.3'),1]])
            u=mp.matrix([1,mp.mpf('.7')])
            for s in [0,2,4]:
                ratios=[]
                for n in [10,100,1000]:
                    T=mp.diag([n*n,(n+1)**2]); D=T*(T+4*mp.eye(2)); K=T+4*mp.eye(2)+Q
                    q=K*u; recovered=mp.lu_solve(K,q)
                    norm=lambda x,k:sum(mp.norm(W*mp.diag([mp.mpf(n)**j,mp.mpf(n+1)**j])*x) for j in range(k+1))
                    ratio=norm(D*u,s)/norm(q,s+2); ratios.append(ratio)
                    check(f'actual_model_K_inverse_identity_{dps}_{s}_{n}',mp.norm(recovered-u)<mp.mpf('1e-45') and mp.norm(T*Q-Q*T)>0)
                check(f'positive_transfer_model_limit_{dps}_{s}',abs(ratios[-1]-1)<mp.mpf('.002') and abs(ratios[-1]-1)<abs(ratios[0]-1))
                records.append({'dps':dps,'order':s,'ratios':[mp.nstr(v,20) for v in ratios]})
    FIXTURES.append({'fixture':'noncommuting_positive_transfer_model','records':records,
                     'Q_positive_and_noncommuting':True,'actual_quantum_operator_sampled':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths(); hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes(); error=None
    print(json.dumps({'event':'start','story':'S167 construction feasibility audit','input_count':len(ps)}),flush=True)
    try:
        gl,gu=inputs_and_budgets();exact_inequalities();polys=tensor_polynomials()
        numerical_normal(gl,gu,polys);normal_packet_model(gl,polys);noncommuting_positive_model()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
    run={'story':'S167 construction feasibility audit','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
         'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after}
    for k,v in review().items():
        if isinstance(v,bool) and k not in run:run[k]=v
    for kind,p in zip(['checker','note','review'],OWN):run['frozen_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
