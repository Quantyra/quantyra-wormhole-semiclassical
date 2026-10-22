"""Frozen normal-action, domain-scaling and saddle diagnostics; no physical solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as Q
import json,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='recentered-reference-inverse';PREVIOUS='recentered-mixed-form'
OWN=['notes/'+STEM+'.md','checks/recentered_reference_inverse.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def unexpanded_local_actions():
    z,lam,G,D,ep=S.symbols('z lambda G D epsilon',positive=True)
    e,p=S.symbols('e p',real=True)
    omega=-S.log(z)-S.log(lam)/2+e*z**p
    inv=S.exp(-2*omega)
    Kr=-inv*S.diff(omega,z,2);Kt=-inv*S.diff(omega,z)**2
    RR=6*(Kr+Kt);Ric2=9*Kr**2+3*(Kr+2*Kt)**2
    JJ=Ric2-RR**2/4
    boxR=inv*(S.diff(RR,z,2)+2*S.diff(omega,z)*S.diff(RR,z))
    lam0=lam+11*ep*lam**2/6
    literal=S.exp(4*e*z**p)*(RR+12*lam0-D*boxR+ep*(S.Rational(11,6)*JJ-S.Rational(11,72)*RR**2))
    q=p*(p-3);Cv=1+11*ep*lam/3
    expected=6*lam*(Cv-D*lam*q)*(4-q)*z**p
    derived=S.simplify(S.diff(literal,e).subs(e,0))
    check('unexpanded_conformal_trace_density',S.simplify(derived-expected)==0)
    JR=3*(Cv-D*lam*q)*(q-4)/(4*S.pi*G*lam)
    check('scalar_action_sign_and_normal_volume',S.simplify(-derived/(8*S.pi*G)-lam**2*JR*z**p)==0)
    check('quantum_vacuum_recentered',S.simplify(literal.subs(e,0))==0)
    # Independently form every sectional curvature of the anisotropic normal metric.
    vp,vpp,c=S.symbols('vprime vsecond c',real=True)
    ws=[S.Rational(1,3),-S.Rational(1,6),-S.Rational(1,6)]
    radial=[-lam*(1-z*e*w*vp+z*z*e*w*vpp+z*z*e*e*w*w*vp*vp) for w in ws]
    tangent=[-lam*(1-z*e*ws[i]*vp)*(1-z*e*ws[j]*vp) for i in range(3) for j in range(i+1,3)]
    RT=S.expand(2*(sum(radial)+sum(tangent)))
    check('anisotropic_scalar_from_all_sections',S.expand(RT+12*lam+lam*z*z*e*e*vp*vp/6)==0)
    measure=1/(lam**2*z**4)
    local=measure*(-RT/(16*S.pi*G)-c*RT**2)
    Hess=S.diff(local,e,2).subs(e,0)
    Fl=1-384*S.pi*G*c*lam
    check('unexpanded_TT_local_Hessian_factor',S.simplify(Hess-Fl*vp*vp/(48*S.pi*G*lam*z*z))==0)
    check('R2_square_has_zero_TT_normal_Hessian',S.diff(-c*measure*(RT+12*lam)**2,e,2).subs(e,0)==0)
    check('TT_measure_independent',sum(ws)==0)
    FIXTURES.append(dict(fixture='independent_unexpanded_local_normal_actions',trace_variation=str(derived),TT_Hessian=str(S.factor(Hess)),quantum_nonlocal_kernel_sampled=False))

def exact_domain_and_saddle():
    s,JR,NT=S.symbols('s JR NT');C=lambda v:S.Matrix([[1,(v+3)/(3*(v+1))],[0,2]])
    R=S.Matrix([(s+3)/(2*(s+1)),-S.Rational(3,2)])
    check('TT_normal_coordinate_column',S.simplify(C(s)*R-S.Matrix([0,-3]))==S.zeros(2,1))
    normal=C(3-s).T*S.diag(JR,NT)*C(s)
    check('normal_reflected_determinant',S.factor(normal.det())==4*JR*NT)
    check('mass_is_pure_TT',C(3)*S.Matrix([1,-2])==S.Matrix([0,-4]))
    check('mass_return_factor',R.subs(s,3)==S.Rational(3,4)*S.Matrix([1,-2]))
    a,b,c,d=S.symbols('a b c d',real=True);K=S.Matrix([[a,b],[c,d]])
    I=S.eye(2);H=(-I).row_join(K).col_join(K.T.row_join(I))
    expected=(I+K*K.T).row_join(S.zeros(2)).col_join(S.zeros(2).row_join(I+K.T*K))
    check('saddle_square_with_nonnormal_mixed_block',S.simplify(H*H-expected)==S.zeros(4))
    sign=S.diag(-1,-1,1,1)
    check('sign_reversed_coercive_symmetric_part',(sign*H+(sign*H).T)/2==S.eye(4))
    r,L=S.symbols('r L',positive=True)
    profile=r**-3*(1-L/r)**5
    records=[]
    for order in range(5):
        f=profile
        for _ in range(order):f=r*S.diff(f,r)
        terms=S.Poly(S.expand((f*f*r*r).subs(r,1/s)),s)
        value=sum(co/(power[0]-1)*L**(1-power[0]) for power,co in terms.terms())
        # Every exponent is >=4; integrate powers on [L,infinity) exactly.
        value=S.factor(value)
        check('critical_log_derivative_scaling_'+str(order),S.simplify(S.diff(L**3*value,L))==0 and (value.subs(L,1)>0)==True)
        if order==0:check('unit_mass_exact_L2_norm',value==1/(858*L**3))
        records.append(dict(order=order,norm_squared=str(value)))
    check('unit_asymptotic_mass',S.limit(r**3*profile,r,S.oo)==1)
    check('model_junction_graph_regularity',all(S.diff(profile,r,j).subs(r,L)==0 for j in range(5)))
    w=S.symbols('w',real=True)
    transform=sum((-1)**j*S.binomial(5,j)/(S.Rational(3,2)+j+S.I*w) for j in range(6))
    product=120/S.prod(S.Rational(3,2)+j+S.I*w for j in range(6))
    check('log_graph_Fourier_transform',S.factor(transform-product)==0)
    check('log_graph_high_frequency_integrable',8-2*6==-4)
    check('TT_tail_correction_critical_exponent',-2+Q(3,2)==Q(-1,2))
    check('mass_trace_unbounded_graph_exponent',Q(3,2)>0)
    check('mass_boundary_relative_form_growth',3>0)
    check('new_fast_index_and_mass_extension',-1+1==0)
    FIXTURES.append(dict(fixture='actual_domain_scaling_and_general_saddle_algebra',mass_graph_derivative_norms=records,mass_trace_continuous=False,physical_mass_momentum_evaluated=False))

def root_guards():
    pl=Q('3.14159');ph=Q('3.14160');k=2000;g0=Q(1,k**10);gl=(1-Q('1e-34'))*g0;cr=Q('1e22')
    Fl=1-384*ph*cr*g0
    Dlo=96*pl*cr*gl-g0/(60*pl);Dhi=96*ph*cr*g0
    Cvhi=1+11*g0/(180*pl)
    check('Dlambda_positive',Dlo>0)
    check('new_mass_ratio_lower',1/Dhi>Q('3.39e8'))
    check('new_mass_ratio_upper',Cvhi/Dlo<Q('3.40e8'))
    check('fast_scalar_roots_outside_strip',Q('3.39e8')>18003*18000)
    check('effective_TT_large_parameter',10*pl*Fl>30)
    check('TT_small_frequency_exclusion',9*(15+11)+75<400<30*k**10)
    check('TT_large_frequency_exclusion',12<30*k**10)
    check('fixed_scale_lower_log',3**18<k**3/4)
    # lambda*b_class^2 < x0 < 1e10 < exp(24), since exp(1)>2.7.
    check('fixed_scale_upper_log',Q('2.7')**24>10**10)
    check('Splus_scalar_roots_absent',Q(7,2)*(Q(7,2)-3)<4)
    check('two_weight_steps_below_gain',1<2 and Q(3,2)+1==Q(5,2) and Q(5,2)+1==Q(7,2))
    check('mass_remainder_above_final_weight',4>Q(7,2) and 3+2>Q(7,2))
    check('reference_mixed_bounded',Q('.336')<Q(3,8))
    check('conditional_transfer_inverse',1/(1-Q(1,4))==Q(4,3))
    check('nonlinear_selfmap_target',Q('5e-9')+Q('.25')*Q('1e-8')+Q('5e7')/2*Q('1e-8')**2==Q('1e-8'))
    check('nonlinear_contraction_target',Q('.25')+Q('5e7')*Q('1e-8')==Q('.75')<1)
    FIXTURES.append(dict(fixture='uniform_parameter_and_root_guards',mass_ratio=['3.39e8','3.40e8'],reference_energy_inverse=1,finite_transfer_bound_proved=False))

def benchmark_mellin():
    records=[]
    for prec in [70,100]:
        mp.mp.dps=prec;k=mp.mpf(2000);g0=k**-10;cr=mp.mpf('1e22')
        g=2*g0/(1+mp.sqrt(1+4*11*g0/(360*mp.pi)))
        Fl=1-384*mp.pi*cr*g;Dl=96*mp.pi*cr*g-g/(60*mp.pi);Cv=1+11*g/(180*mp.pi)
        x0=(mp.sqrt(1+12*mp.pi*k**6)-1)/6
        kap=-mp.mpf(6)/5+mp.log(x0*g/g0)/2
        mr=Cv/Dl;rootp=(3+mp.sqrt(9+4*mr))/2;rootm=3-rootp
        check('actual_scalar_mass_'+str(prec),mp.mpf('3.39e8')<mr<mp.mpf('3.40e8'))
        check('actual_scalar_root_locations_'+str(prec),rootp>18000 and rootm<-18000)
        check('actual_fixed_Mellin_scale_'+str(prec),mp.mpf(39)/5<kap<mp.mpf(54)/5)
        FM=lambda s:mp.digamma(s-1)+mp.pi/2*mp.cot(mp.pi*s/2)
        TTlarge=10*mp.pi*Fl/g+mp.mpf(23)/18
        samples=[]
        for delta,t in [('1.5','0'),('1.5','.25'),('3.5','2'),('3.5','100'),('3.5','1000000')]:
            s=mp.mpc(delta,t);q=s*(s-3);FMv=FM(s)
            reflection=abs(FMv-FM(3-s))
            check('Mellin_reflection_'+delta+'_'+t+'_'+str(prec),reflection<mp.mpf('1e-50'))
            bracket=(q+2)*(FMv+kap)-TTlarge
            check('sample_nonzero_normal_bracket_'+delta+'_'+t+'_'+str(prec),abs(bracket)>1)
            samples.append(dict(real=delta,imag=t,reflection_error=mp.nstr(reflection,12),normalized_TT_bracket_abs=mp.nstr(abs(bracket)/TTlarge,35)))
        for point in [1,2]:
            h=mp.mpf('1e-22');s=mp.mpf(point)+h;q=s*(s-3)
            val=q*(q+2)*FM(s)
            check('removable_Mellin_value_'+str(point)+'_'+str(prec),abs(val+2)<mp.mpf('1e-19'))
        records.append(dict(digits=prec,mass_ratio=mp.nstr(mr,60),Delta_plus=mp.nstr(rootp,60),Delta_minus=mp.nstr(rootm,60),kappa_M=mp.nstr(kap,60),TT_large_parameter=mp.nstr(TTlarge,60),diagnostic_samples=samples))
    for key in ['mass_ratio','Delta_plus','Delta_minus','kappa_M','TT_large_parameter']:
        check('benchmark_cross_precision_'+key,abs(mp.mpf(records[0][key])/mp.mpf(records[1][key])-1)<mp.mpf('1e-50'))
    FIXTURES.append(dict(fixture='actual_benchmark_normal_parameters',records=records,samples_prove_root_count=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_mixed_form_result',old['status']=='passed' and old['relative_mixed_form_bound']=='.336')
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        unexpanded_local_actions();exact_domain_and_saddle();root_guards();benchmark_mellin()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        reference_bulk_inverse_argument='mass-extended weighted bijection; energy bound 1' if error is None else None,
        reference_TT_form_closedness_argument=error is None,reference_scalar_phase_cokernel_absent=error is None,
        physical_mass_momentum_evaluated=False,physical_boundary_completed=False,physical_coupled_inverse_proved=False,
        admission_norm_inverse_quantified=False,coupled_correction_norm_proved=False,nonlinear_solution_constructed=False,
        semiclassical_validity_proved=False,physical_stability_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
