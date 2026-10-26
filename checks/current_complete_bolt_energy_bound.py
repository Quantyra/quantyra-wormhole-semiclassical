"""Complete two-row bolt bound. Does not apply the quantum inverse."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import sys, json, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import sympy as sp
import flint
from flint import arb, ctx
ROOT=Path(__file__).resolve().parents[1]
STEM='current-complete-bolt-energy-bound'
OWN=['notes/'+STEM+'.md','checks/current_complete_bolt_energy_bound.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+['notes/current-direct-triple-return-attempt.md','notes/current-seed-full-hessian-realization.md',
 'notes/recentered-mixed-form.md','notes/stationary-seed-quantum-form.md',
 'notes/current-global-mixed-schur-preconditioner.md','notes/current-schur-residual-phase-domain.md',
 'notes/current-coupled-input-energy-bound.md','notes/current-coupled-truncation-energy-error.md',
 'notes/current-schur-polynomial-return-test.md','notes/current-weighted-bolt-source-return.md',
 'evidence/current-weighted-bolt-source-return.json','evidence/current-localized-bolt-adjoint-sensitivity.json',
 'evidence/current-complete-source-square-factorization.json']
CHECKS=[]
def sh(p):return sha256(p.read_bytes()).hexdigest()
def check(n,ok):
    CHECKS.append({'name':n,'passed':bool(ok)});print(n+(': passed' if ok else ': FAILED'),flush=True)
    if not ok:raise AssertionError(n)
def enc(x):
    m,r,e=x.mid_rad_10exp(50)
    return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
    if isinstance(x,dict):
        e=x['decimal_enclosure'];return arb(e['mid']+' +/- '+e['rad'])*arb(10)**int(e['exponent'])
    return arb(x.strip('[]'))
def prior(stem):
    r=json.loads((ROOT/'evidence'/f'{stem}.json').read_bytes())['runs'][-1]
    check('accepted_'+stem,r['status']=='passed');return r
def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('ordinary_argument_checker_frozen',all(sh(ROOT/p)==h for p,h in review['scientific_sha256'].items()))
    src=prior('current-localized-bolt-adjoint-sensitivity')
    out=prior('current-weighted-bolt-source-return')
    square=prior('current-complete-source-square-factorization')
    pars=src['result']['parameters']
    h,gamma,beta,C,B2=[dec(pars[k]) for k in ['h','gamma','beta','C','B2']]
    # Actual bolt relation tau=h^2 K(1+h)/(1+h)^3; K<2 on this seed.
    check('finite_seed_parameter_envelopes',C*beta>arb('8e9') and 2*h*h/beta<arb('1e-50') and beta<arb('.167'))
    # Entire finite seed sectional curvatures, including d and temperature.
    t,u,e,v=sp.symbols('t u e v',positive=True)
    lapse=t*t+u-(4+2*u+v)/t+(3+u+e)/t**2-2*e/t**3+e/t**4
    curvs=[-sp.diff(lapse,t,2)/2,-sp.diff(lapse,t)/(2*t),(u-lapse)/t**2]
    majorants=[]
    for curvature in curvs:
        for j in range(3):
            expr=sp.expand(t**j*sp.diff(curvature,t,j));total=sp.Rational(0)
            for term in sp.Add.make_args(expr):
                powers=term.as_powers_dict();check('sectional_nonpositive_radial_power',powers.get(t,0)<=0)
                cterm=term/t**powers.get(t,0)
                coeff,monomial=cterm.as_coeff_Mul()
                total+=abs(coeff)*monomial.subs({u:sp.Rational(1,8000000000),e:sp.Rational(1,10**27),v:sp.Rational(1,10**50)})
            majorants.append(str(total));check('whole_sectional_two_jet_coefficient_sum',total<10000)
    # Quantitative contraction counts in the ordinary tensor proof.
    wedge2=4*(2*arb(2).sqrt()*100+2*20**2)
    check('projector_second_product_guard',wedge2<4400)
    check('full_curvature_two_jet_guard',3*arb('1e6')*(4400+2*4*2*arb(2).sqrt()*20+4*2)<arb('1e12'))
    check('integrated_curvature_and_TT_guard',arb('1e19')*arb(22).sqrt()<arb('1e20'))
    d=dec(square['result']['parameters']['d']);g=C*beta/B2;epsbar=C/(60*arb.pi()*B2)
    q=arb(169)/500;f=q/(1+q);ellb=(1+h)/h
    alpha=(1+arb(21).sqrt())/2;Kxi=1+4/(3*(alpha-2))
    check('poisson_profile_integrable',alpha>2 and Kxi<3)
    S0=1/arb('23.99').sqrt();S1=(beta/arb('5.99')).sqrt();S2=(beta/(6*d)).sqrt()
    SH=(S2*S2+13*beta*S1*S1).sqrt()
    Cs=S0/arb('.99').sqrt()+4*(S2+2*arb('.17').sqrt()*S1)/gamma
    w0=arb('1.1e-27')*beta;w2=arb('4.3e-25')*beta**2;Cv=1+11*g/(180*arb.pi())
    pieces=[4*arb('3.01')*beta*((Cv+11*epsbar*w0/36)*S0+d*S2),
        arb(22).sqrt()*d*w2*S0,epsbar*arb('1e20')*(SH+S1+S0),
        20/(60*arb.pi())*arb('4e-27')*beta*S0,arb('9e-28')*beta*5/beta.sqrt()*S1]
    Kb=sum(pieces);check('complete_mixed_adjoint_below_frozen_guard',Kb<1)
    c=2*arb(17).sqrt()/9;KV=1+c/arb('4.5').sqrt()+c/2
    Kout=ellb*(arb('2.5')*h)**(arb(7)/5)*3000/gamma.sqrt()
    check('same_original_weighted_point_map',Kout.overlaps(dec(out['result']['actual_output_factors']['original_physical_bolt_factor_from_weighted_L2'])))
    J0=ellb*(Kxi*arb('7.477e50')/2+Cs*q.sqrt())
    J1=Kout*h**(-arb(7)/5)*KV*q.sqrt()
    U3=J0*f**3+J1*f**2;U1=J0*f+J1
    # Exact polynomial denominators, no repeated interval polynomial evaluation.
    xr=1+sp.Rational(2)/sp.Rational(169,500)
    ds=[sp.Integer(1),xr]
    for m in range(2,130):ds.append(sp.cancel(2*xr*ds[-1]-ds[-2]))
    def aa(r):return arb(str(sp.numer(r)))/arb(str(sp.denom(r)))
    budget=arb('1.667e-9');paff=arb('1.1e-10')*(1+q)**3/aa(ds[29]);debit=paff*U3
    check('same_P32_affine_prefactor',paff.overlaps(dec(square['result']['complete_budget']['p_aff'])))
    check('specified_P32_upper_estimator_exceeds_full_allowance',debit>budget)
    rows=[]
    for N in range(3,100):
        D=aa(ds[N-3]);h2=(1+q)**2/D;h3=(1+q)**3/D;tail=arb('1.1e-10')*h3*U3
        rows.append((N,h2,h3,tail,D/aa(ds[29])))
    full=next(r for r in rows if r[3]<budget);half=next(r for r in rows if r[3]<budget/2)
    for label,r,margin in [('full',full,budget),('half',half,budget/2)]:
        previous=next(x for x in rows if x[0]==r[0]-1)
        check('minimal_degree_'+label,previous[3]>margin and r[3]<margin)
    ceilings={}
    for r in [next(x for x in rows if x[0]==32),full,half]:
        N,h2,h3,tail,ratio=r;pref=arb('1.1e-10')*h3
        ceilings[str(N)]={'h2':enc(h2),'h3':enc(h3),'bolt_tail_upper':enc(tail),'denominator_gain':enc(ratio),
          'conditional_defect_Cdef3_IY_ceiling_application_zero':enc(arb('.5')/h3),
          'conditional_graph_Cg3_ceiling_other_errors_zero':enc(arb('5.34235e-12')/pref),
          'conditional_physical_CX3_ceiling_other_errors_zero':enc(arb('5e-9')/pref),
          'conditional_mass_Cmass3_Im_ceiling_center_application_zero':enc(arb('5e-5')/h2),
          'conditional_feedback_Cg3_Ifb_ceiling_other_errors_zero':enc(arb('.0012069')/h2),
          'actual_Cdef3_IY':None,'actual_Cg3':None,'actual_CX3':None,'actual_mass_column':None,
          'actual_feedback':None,'certified_application_error':None,'actual_remaining_margins':None}
    eps_target=(budget/2)/U1;epsH=arb('1.020e-25');epsin=arb('1.1e-10')*epsH/(1-epsH)
    appdebit=U1*epsin
    check('favorable_F_filtered_application_upper_estimator_fails',appdebit>budget/2)
    check('scalar_point_term_subdominant',Cs*q.sqrt()<Kxi*arb('7.477e50')/2)
    check('full_mixed_B_return_subdominant',J1*f*f<J0*f**3)
    return {'precision_bits':512,'finite_seed_sectional_majorants':majorants,
      'constants':{k:enc(v) for k,v in {'Kxi':Kxi,'scalar_point_Cs':Cs,'scalar_Hessian_SH':SH,'Kb':Kb,'KV':KV,
          'J0':J0,'J1':J1,'complete_Cb3_upper':U3,'complete_Cb1_upper':U1,'P32_affine_bolt_upper':debit,
          'P32_estimator_improvement_required':debit/budget,'first_term_fraction':J0*f**3/U3}.items()},
      'mixed_adjoint_terms':[enc(x) for x in pieces],
      'polynomial_test':{'first_degree_full_bolt_margin':full[0],'first_degree_half_bolt_margin':half[0],
          'conditional_ceilings_are_not_coefficients':True,'rows':ceilings,'higher_degree_recommended':False},
      'application_test':{'assumption':'error passes through F; unfiltered final errors still lack a numerical output comparison',
          'half_bolt_margin_E_error_ceiling':enc(eps_target),'existing_energy_perturbation_upper':enc(epsH),
          'affine_inverse_energy_error_expression':enc(epsin),'F_filtered_application_upper_expression':enc(appdebit),
          'application_estimator_improvement_required':enc(appdebit/(budget/2)),
          'actual_original_norm_application_error':None,'actual_quantum_inverse_applied':False},
      'scope':{'complete_three_residual_bolt_map_finite_bound':True,'P32_positive_upper_certificate_fails':True,
          'actual_P32_norm_lower_bound':False,'original_defect_certified':False,'graph_physical_maps_certified':False,
          'mass_column_certified':False,'complete_coupled_estimate_certified':False,'construction_obstruction':False,
          'nonlinear_or_physical_gate_closed':False,'quantum_fast_scope_unchanged':True}}
if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    path=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={p:sh(ROOT/p) for p in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
      'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
      'runtime':{'python':sys.version,'python_sha256':sh(Path(sys.executable)),'flint_version':flint.__version__,'sympy_version':sp.__version__}}
    try:
        run['result']=execute();run['input_sha256_after']={p:sh(ROOT/p) for p in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
    except Exception:run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(path.read_bytes()) if path.exists() else {'runs':[]}
    ledger['runs'].append(run);path.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    summary=json.dumps({'status':run['status'],'mode':mode,'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as stream:stream.write(summary.encode())
    print(summary,flush=True);sys.exit(0 if run['status']=='passed' else 1)
