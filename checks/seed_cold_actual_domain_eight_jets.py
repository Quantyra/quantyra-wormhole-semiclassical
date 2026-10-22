"""Supporting checks for the actual-domain cold C8 source completion."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
from functools import lru_cache
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='seed-cold-actual-domain-eight-jets';PREVIOUS='seed-cold-gluing-eight-jets'
OWN=['checks/seed_cold_actual_domain_eight_jets.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational;ROWS=[];FIXTURES=[]
PROTOCOL={'exact_arithmetic':'symbolic and rational','numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','source_order':8,'energy_order':10,
          'ordinary_continuum_proof_separate':True,'actual_charged_source_evaluated':False,
          'actual_inverse_differentiated':False,'boundary_or_state_changed':False,
          'Hermitian_principal_Clifford_convention':True,'prior_convention_correction_recorded':True}


def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)


def zero(name,expr):
    v=S.simplify(expr);check(name,v==0,v)


def previous():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(previous()['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})


def provenance():
    p=previous()
    check('prior_gluing_passed',p['status']=='passed' and len(p['checks'])==113 and all(r['passed'] for r in p['checks']))
    check('prior_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('all_prior_inputs_unchanged',not bad,bad)
    check('retained_sources_pinned',set(review()['retained_inputs'])<=set(paths()))
    check('prior_convention_error_explicit',
          'gamma matrices are anti-Hermitian' in (ROOT/'notes/seed-cold-gluing-eight-jets.md').read_text(encoding='utf-8')
          and review()['prior_principal_Clifford_convention_corrected'])
    check('complete_cold_not_nonlinear',review()['complete_actual_cold_source_C8_bounded']
          and review()['all_regional_seed_source_C8_bounded'] and not review()['actual_metric_tenth_derivative_bounded']
          and not review()['nonlinear_solution_constructed'] and not review()['goal_complete'])


def exp_upper(x,n=32):
    return sum(x**j/S.factorial(j) for j in range(n))+x**n/S.factorial(n)/(1-x/(n+1))


def corrected_physical_symbol():
    gam=[S.Matrix([[0,1],[1,0]]),S.Matrix([[0,-S.I],[S.I,0]])]
    J=gam[0]*gam[1];unit=S.eye(2);nil=S.zeros(2)
    for i in range(2):
        check('Hermitian_principal_'+str(i),gam[i].H==gam[i] and gam[i]*gam[i]==unit)
        check('spin_anticommutation_'+str(i),J*gam[i]+gam[i]*J==nil)
    check('positive_Clifford_pair',gam[0]*gam[1]+gam[1]*gam[0]==nil)
    check('unitary_spin_generator',J.H==-J and J*J==-unit)
    a,b,theta,lam=S.symbols('a b theta lam',real=True)
    symbol=S.I*(a*gam[0]+b*gam[1])
    check('anti_self_adjoint_Fourier_symbol',symbol.H==-symbol)
    check('Plancherel_gradient_identity',S.simplify(symbol.H*symbol-(a*a+b*b)*unit)==nil)
    wrong=S.I*(a*S.I*gam[0]+b*S.I*gam[1])
    check('wrong_antiHermitian_principal_detected',(wrong.H-wrong).applyfunc(S.simplify)==nil
          and (wrong.H+wrong).applyfunc(S.simplify)!=nil)
    Q=S.cos(theta/2)*unit+S.sin(theta/2)*J;Qi=Q.H
    K=S.Matrix([[S.cos(theta),S.sin(theta)],[-S.sin(theta),S.cos(theta)]])/lam
    for row in range(2):
        aligned=Q*sum((gam[i]*K[row,i] for i in range(2)),nil)*Qi
        check('physical_frame_alignment_'+str(row),(aligned-gam[row]/lam).applyfunc(S.trigsimp)==nil)
    x,y,e,tau=S.symbols('x y e tau',real=True);coords=[x,y]
    G=S.Matrix([x+e*x*x,y+e*x*y]);DG=G.jacobian(coords);det=DG.det();Ki=DG.inv()
    P=[sum((gam[i]*Ki[row,i] for i in range(2)),nil) for row in range(2)]
    direct=-sum((P[row]*S.diff(det,coords[row])/det/2 for row in range(2)),nil)
    halfdiv=sum((P[row].diff(coords[row])/2 for row in range(2)),nil)
    check('physical_density_transport',(direct-halfdiv).applyfunc(S.simplify)==nil)
    check('formal_anti_self_adjoint_half_density',all(v.H==v for v in P) and (halfdiv.H-halfdiv).applyfunc(S.simplify)==nil)
    check('nonzero_density_term',direct.subs({x:0,y:0})==-R(3,2)*e*gam[0])
    q=S.cos(tau*x/2)*unit+S.sin(tau*x/2)*J
    z=(q*gam[0]*q.H.diff(x)).subs(x,0)
    check('nonzero_spin_term',z==-tau*gam[0]*J/2 and z!=nil)
    sigma=S.Function('sigma')(x,y);f=S.Matrix([S.Function('f')(x,y),S.Function('g')(x,y)])
    D=lambda v:sum((gam[i]*v.diff(coords[i]) for i in range(2)),nil[:,0])
    physical=S.exp(-sigma)*(D(f)+sum((gam[i]*sigma.diff(coords[i])*f/2 for i in range(2)),nil[:,0]))
    check('normal_conformal_mass_potential',(S.exp(R(3,2)*sigma)*physical-D(S.exp(sigma/2)*f)).applyfunc(S.simplify)==S.zeros(2,1))
    FIXTURES.append({'fixture':'corrected_physical_symbol_and_transport',
                     'principal_matrices':'Hermitian positive Clifford','Fourier_symbol':'anti-Hermitian',
                     'rotation_alignment_verified':True,'density_and_spin_retained':True,
                     'prior_norm_ledger_changed':False,'prior_frozen_sources_preserved':True})


def row_derivative_budgets():
    check('centered_separation_guard',R(94,100)*R(1006,1000)*exp_upper(R(51,1000))<1)
    check('normal_disk_margin',R(2,18)+R(2,16)<R(1,4))
    check('full_tensor_Cauchy_factor',24**2>2*16**2)
    check('scaled_coefficient_guard',R(24,18)==R(4,3) and exp_upper(R(25,1000))<R(11,10))
    check('observation_distance_loss',2*exp_upper(R(55,1000))<R(22,10))
    check('spin_Hilbert_row_factor',10*R(1415,1000)*exp_upper(R(3,100))<16)
    check('row_initial_L2_budget',4*16*exp_upper(R(55,100)+R(125,10000))<200)
    U=[S.Integer(1)];ratios=[]
    for n in range(1,11):
        value=26*U[n-1]+R(11,10)*sum(S.binomial(n-1,r)*S.factorial(r)*R(4,3)**r*U[n-1-r] for r in range(1,n))
        U.append(value)
        check('nested_energy_'+str(n),value<S.factorial(n)*50**n,S.N(value/(S.factorial(n)*50**n),18))
        ratio=(26+R(11,10)*sum(R(2,75)**r for r in range(1,n)))/(50*n)
        check('factorial_energy_induction_'+str(n),ratio<1,ratio);ratios.append(str(ratio))
    check('ten_energy_orders_leave_point_disk',2-R(10,12)>1)
    x,y,t=S.symbols('x y t',real=True);eta=(1-x*x-y*y)**3
    lap=S.diff(eta,x,2)+S.diff(eta,y,2)
    zero('point_cutoff_laplacian',lap-(-12+48*(x*x+y*y)-36*(x*x+y*y)**2))
    check('point_cutoff_laplacian_bound',all(abs((-12+48*t-36*t*t).subs(t,k))<=12 for k in [0,R(2,3),1]))
    for j in range(9):
        budget=13+600*(j+1)+5000*(j+1)*(j+2)
        check('point_evaluation_'+str(j),budget<500000,budget)
    check('endpoint_product_and_row_reserve',500000*200*R(11,10)/(1-R(4,3)/50)<2*10**8)
    check('trace_commutator_prefactor',2*9**2*(2*10**8)**2<10**19)
    check('normal_connection_analytic_norm',2*R(1415,1000)*R(25,1000)/R(1,16)<2)
    check('normal_connection_tensor_factor',25**2>2*16**2)
    FIXTURES.append({'fixture':'interior_row_derivative_budgets','energy_induction_ratios':ratios,
                     'energy_order':10,'row_order':8,'row_prefactor':'2e8','row_derivative_scale':'100M',
                     'actual_inverse_differentiated':False,'distant_cutoff_differentiated':False,
                     'physical_spin_Hilbert_factor_included':True})


def weight_and_angular():
    tlo=R(11511,5000);c=R(233,250);A=18*c
    check('log10_lower_certificate',exp_upper(tlo)<10)
    check('cold_decay_exponent',R(47,100)*tlo-R(15,100)>c)
    check('cold_weight_denominator',1+tlo>3)
    check('row_separation_through_disk',R(94,100)*tlo-R(3,10)-R(22,10)/18>1)
    check('cold_weight_monotonicity',R(47,100)*18-1>0)
    check('angular_exponential_certificate',sum(A**j/S.factorial(j) for j in range(81))>R(10**8,6))
    n,q=S.symbols('n q',positive=True)
    zero('degeneracy_comparison_identity',6*q*n-2*(q+2*n)-(2*q*(n-1)+4*n*(q-1)))
    entries=[]
    for j in range(9):
        check('mass_moment_monotonicity_'+str(j),18*c>j)
        check('index_moment_monotonicity_'+str(j),A>j+2)
        check('angular_tail_integrable_'+str(j),A-j-3>0)
        L=6*18**j*R(6,10**8)*(1+2/(A-j-3))
        entries.append({'order':j,'moment_bound':str(L)})
    FIXTURES.append({'fixture':'cold_weight_and_infinite_angular_bounds','log10_rational_lower':str(tlo),
                     'exponential_rate':str(c),'A':str(A),'entries':entries,
                     'all_angular_modes_included':True,'bolt_radial_division_used':False})


def physical_and_regional_ledger():
    A=R(2097,125)
    H=[10**19*(j+1)*R(352,100000)*R(1000,3)*6*18**j*R(6,10**8)*(1+2/(A-j-3)) for j in range(9)]
    @lru_cache(None)
    def B(r,k):
        if r==0:return S.factorial(k)*100**k*H[k]
        return B(r-1,k+1)+(r-1)*sum(S.binomial(k,a)*2*S.factorial(a)*25**a*B(r-1,k-a) for a in range(1,k+1))
    gluing=previous()['fixtures'][3]['entries']
    local=json.loads((ROOT/'evidence/seed-cold-local-eight-jets.json').read_bytes())['runs'][-1]['fixtures'][3]['entries']
    entries=[]
    for j in range(9):
        J=B(j,0)/(S.factorial(j)*100**j)
        check('actual_domain_physical_'+str(j),J<10**24,S.N(J,22))
        # Both predecessor entries are already physical tensors with the same a_x weight.
        combined=J+R(gluing[j]['physical_exact'])+R(local[j]['exact'])
        check('complete_cold_source_'+str(j),combined<10**44,S.N(combined,22))
        entries.append({'order':j,'coordinate_exact':str(H[j]),'physical_exact':str(J),
                        'physical_decimal':str(S.N(J,22)),'complete_cold_exact':str(combined),
                        'complete_cold_decimal':str(S.N(combined,22))})
    for stem in ['seed-middle-eight-jets','seed-outer-eight-jets']:
        prior=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
        check('regional_input_passed_'+stem,prior['status']=='passed' and all(c['passed'] for c in prior['checks']))
    check('strict_combination_margin',10**24+R(237,100)*10**43+10**8<10**44)
    FIXTURES.append({'fixture':'full_physical_and_regional_ledger','entries':entries,
                     'actual_domain_reserve':'1e24','complete_actual_cold_reserve':'1e44',
                     'Gamma_zero_order_vanishes_in_normal_chart':True,'higher_tensors_symmetrized':False,
                     'predecessor_tensors_combined_after_physical_conversion':True,
                     'metric_tenth_derivative_bounded':False})


def independent_mass_solution():
    x=S.symbols('x',real=True);kap=R(1,5);polys=[S.Integer(1)]
    for j in range(10):polys.append(S.expand(S.diff(polys[-1],x)-(1+kap*x)*polys[-1]))
    for j in [1,4,8,10]:
        herm=(-S.sqrt(kap/2))**j*S.hermite(j,(1+kap*x)/S.sqrt(2*kap))
        zero('independent_Hermite_derivative_'+str(j),herm-polys[j])
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            coeffs=[[mp.mpf(str(c.p))/mp.mpf(str(c.q)) for c in S.Poly(p,x).all_coeffs()] for p in polys]
            norms=[]
            for n in range(11):
                radius=2-mp.mpf(n)/12
                f=lambda z:4*mp.sqrt(max(mp.mpf(0),radius*radius-z*z))*mp.exp(-2*z-z*z/5)*mp.polyval(coeffs[n],z)**2
                norms.append(mp.sqrt(mp.quad(f,[-radius,0,radius])))
            for n in range(1,11):
                rhs=26*norms[n-1]+mp.mpf('1.1')*sum(mp.binomial(n-1,r)*mp.factorial(r)*(mp.mpf(4)/3)**r*norms[n-1-r] for r in range(1,n))
                check('independent_energy_'+str(dps)+'_'+str(n),norms[n]<rhs,mp.nstr(norms[n]/rhs,12))
            radius=mp.mpf('.05')
            for j in [1,4,8,10]:
                f=lambda theta:mp.exp(-radius*mp.exp(1j*theta)-(radius*mp.exp(1j*theta))**2/10)*mp.exp(-1j*j*theta)
                contour=mp.quad(f,[0,mp.pi/2,mp.pi,3*mp.pi/2,2*mp.pi])/(2*mp.pi*radius**j)
                rational=polys[j].subs(x,0)
                target=mp.mpf(str(S.numer(rational)))/mp.mpf(str(S.denom(rational)))/mp.factorial(j)
                error=abs(contour-target)
                check('independent_mass_solution_contour_'+str(dps)+'_'+str(j),error<mp.mpf('1e-42'),mp.nstr(error,12))
                records.append({'dps':dps,'order':j,'contour_error':mp.nstr(error,18)})
            for j in [0,4,8]:
                exact=abs(mp.polyval(coeffs[j],0))*mp.sqrt(2)
                rhs=13*norms[j]+12*norms[j+1]+2*norms[j+2]
                check('independent_point_estimate_'+str(dps)+'_'+str(j),exact<rhs,mp.nstr(exact/rhs,12))
    FIXTURES.append({'fixture':'independent_nonconstant_mass_potential_solution','solution':'exp(-x-x^2/10)*(1,1)',
                     'mass_potential':'1+x/5','records':records,'actual_seed_sampled':False,
                     'independent_Hermite_and_contour_routes':True})


def independent_tail_quadrature():
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            A=mp.mpf('16.776')
            for j in [0,4,8]:
                p=j+3
                integral=mp.quad(lambda u:u**p*mp.exp(-A*u),[1,2,4,mp.inf])
                exact=mp.exp(-A)*sum(mp.factorial(p)/mp.factorial(p-k)/A**(k+1) for k in range(p+1))
                error=abs(integral-exact)
                check('independent_tail_integral_'+str(dps)+'_'+str(j),error<mp.mpf('1e-42'),mp.nstr(error,12))
                major=mp.exp(-A)/(A-p)
                check('infinite_tail_majorant_'+str(dps)+'_'+str(j),integral<major,mp.nstr(integral/major,12))
                records.append({'dps':dps,'order':j,'absolute_error':mp.nstr(error,18),'ratio_to_majorant':mp.nstr(integral/major,18)})
    FIXTURES.append({'fixture':'independent_infinite_angular_tail_quadrature','records':records,
                     'finite_angular_cutoff_used':False,'actual_charged_kernel_sampled':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 cold actual domain eight jets','input_count':len(ps)}),flush=True)
    try:
        provenance();corrected_physical_symbol();row_derivative_budgets();weight_and_angular()
        physical_and_regional_ledger();independent_mass_solution();independent_tail_quadrature()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 cold actual domain eight jets','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'actual_domain_C8_bounded':True,'actual_domain_C8_reserve':'1e24',
         'complete_actual_cold_source_C8_bounded':True,'complete_cold_C8_reserve':'1e44',
         'all_regional_seed_source_C8_bounded':True,'full_physical_covariant_tensors_bounded':True,
         'full_angular_tower_retained':True,'prior_principal_Clifford_convention_corrected':True,
         'prior_gluing_numerical_ledger_changed':False,'actual_inverse_differentiated':False,
         'actual_metric_tenth_derivative_bounded':False,'actual_finite_metric_source_H6_bounded':False,
         'actual_transfer_certificate_passed':False,'original_cold_C3_contract_changed':False,
         'endpoint_C5_contract_changed':False,'actual_signed_compatibility_zero_proved':False,
         'nonlinear_solution_constructed':False,'physical_stability_proved':False,
         'physical_source_projection_changed':False,'physical_boundary_condition_changed':False,
         'ordinary_proof_machine_verified':False,'formal_verification':False,
         'S167_complete':False,'S168_complete':False,'goal_complete':False}
    if not good:
        for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');saved=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    saved['runs'].append(run);dest.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
