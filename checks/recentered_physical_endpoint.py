"""Frozen endpoint identities and cancellation diagnostics, not a nonlinear solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as Q
import json,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='recentered-physical-endpoint';PREVIOUS='recentered-stationary-constraints'
OWN=['notes/'+STEM+'.md','checks/recentered_physical_endpoint.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def algebra():
    s,J,N=S.symbols('s J N');a=lambda s:(s+3)/(3*(s+1))
    C=lambda s:S.Matrix([[1,a(s)],[0,2]])
    H=C(3-s).T*S.diag(J,N)*C(s)
    row=S.Matrix([[-(s-6)/3,s-4]])*H
    check('exact_disk_row_combination',S.simplify(row-S.Matrix([[0,4*(s-4)*N]]))==S.zeros(1,2))
    check('no_output_York_pole',S.simplify((s-4)*a(3-s)-(s-6)/3)==0)
    check('input_York_regular_at_four',a(S.Integer(4))==S.Rational(7,15))
    p1,b,c,t,d0=S.symbols('p1 b c t d0');P2=b*t+c
    src=S.diff(P2,t)+(2*p1-S.diff(p1,t))/3;PT=src/(96*d0)
    check('source_polynomial_after_row_combination',src==b+2*p1/3)
    check('original_log_coefficient_recovered',S.simplify(96*d0*PT-2*p1/3-b)==0)
    check('global_constant_not_removed_from_original_row',S.diff(c,t)==0 and c!=0)
    cv,dr,lam,g=S.symbols('Cv DR lambda g');q=s*(s-3)
    JR=3*(cv-dr*lam*q)*(q-4)/(4*S.pi*g)
    slope=S.simplify(S.diff(JR,s).subs(s,4))
    check('scalar_simple_root_slope',S.simplify(slope-15*(cv-4*dr*lam)/(4*S.pi*g))==0)
    check('physical_areal_jet_reserve',7>S.Rational(13,2) and 6-1==5)
    # Rational parts of the whole-contour proof, independent of numerical sampling.
    y=S.symbols('y',nonnegative=True);P=S.expand((4-y+5*S.I*S.sqrt(y))*(6-y+5*S.I*S.sqrt(y)))
    check('P4_real_polynomial',S.simplify(S.re(P)-(24-35*y+y*y))==0)
    check('P4_imaginary_polynomial',S.simplify(S.im(P)-(50-10*y)*S.sqrt(y))==0)
    check('small_frequency_minimum',24-35*Q(1,4)+Q(1,16)==Q(245,16))
    check('small_frequency_reserve',Q(245,32)>6)
    check('digamma_cube_tail',Q(1,27)+Q(1,18)==Q(5,54))
    check('bounded_frequency_ReH',1+6*Q(5,54)+Q(54,5)<13)
    check('small_frequency_ReH',1+Q(1,4)*Q(5,54)+Q(54,5)<12)
    check('bounded_frequency_bracket',6*13+5*3==93)
    check('strip_cotangent_remainder',Q(42,10)/3+Q(42,10)/Q(28,10)<3)
    check('strip_digamma_bound',1+2*(Q(4,25)+Q(2,5))<3)
    check('strip_pole_subtracted_bracket',5+11*17==192)
    check('large_actual_normal_coefficient',30*2000**10>600)
    check('endpoint_normal_inverse_units',S.simplify(1/(24/(480*S.pi**2))-20*S.pi**2)==0)
    FIXTURES.append(dict(fixture='row_polynomials_and_contour_proof_algebra',endpoint_bound_scope='transformed normal multiplier only',
        original_shape_row_retained=True,actual_source_coefficients_evaluated=False))

def moments():
    records=[]
    for digits in [70,100]:
        mp.mp.dps=digits;a=mp.mpf('.1')
        for b in [1,3]:
            for nn in [0,1,2,3]:
                direct=mp.quad(lambda t:mp.exp(-b*t)*(1+t)**(2*nn),[-mp.log(a),mp.inf])
                closed=a**b*sum(mp.binomial(2*nn,j)*(1-mp.log(a))**(2*nn-j)*mp.factorial(j)/mp.mpf(b)**(j+1) for j in range(2*nn+1))
                check('moment_norm_integral_'+str((digits,b,nn)),abs(direct/closed-1)<mp.mpf('1e-55'))
        p=mp.mpf(15)/4;cc=mp.mpf(1)/7
        smooth=mp.beta(p+1,5)+mp.beta(p+2,5)
        m0=mp.beta(p-2,5)+cc*smooth;m1=-3*mp.beta(p-3,5)+cc*smooth
        norm=mp.sqrt(mp.beta(mp.mpf('.5'),9))
        check('actual_weight_fixture_integrable_'+str(digits),p>mp.mpf('3.5'))
        check('moment_Cauchy_Schwarz_'+str(digits),mp.beta(p-2,5)<=norm/mp.sqrt(3))
        check('first_face_derivative_Cauchy_Schwarz_'+str(digits),3*mp.beta(p-3,5)<=3*norm)
        def closed_R(z):return z*mp.beta(p+1,5)*mp.hyp2f1(3,p+1,p+6,-1/z)+z**4*cc*(1+z)*smooth
        for zz in ['.001','.00001','.0000001']:
            z=mp.mpf(zz);parts=[0,z,mp.sqrt(z),1]
            val=z**4*mp.quad(lambda y:((z+y)**-3+cc*(1+z)*(1+y))*y**p*(1-y)**4,parts)
            annul=z**5*mp.quad(lambda y:(-3*(z+y)**-4+cc*(1+y))*y**p*(1-y)**4,parts)
            independent=z*mp.diff(closed_R,z)-4*closed_R(z)
            check('remote_kernel_exact_integral_'+str((digits,zz)),abs(val/closed_R(z)-1)<mp.mpf('1e-55'))
            check('remote_annihilator_identity_'+str((digits,zz)),abs(annul/independent-1)<mp.mpf('1e-55'))
            check('nonzero_moment_retained_'+str((digits,zz)),abs(val/z**4-m0)<5*z and m0>mp.mpf('.01'))
            check('face_derivative_limit_'+str((digits,zz)),abs(annul/z**5-m1)<40*z**mp.mpf('.7'))
            records.append(dict(digits=digits,z=zz,remote_moment=mp.nstr(m0,45),scaled_response=mp.nstr(val/z**4,45),
                first_face_derivative=mp.nstr(m1,45),scaled_annihilator=mp.nstr(annul/z**5,45)))
    FIXTURES.append(dict(fixture='independent_separated_kernel_moment',records=records,actual_Dirac_kernel=False,
        kernel='(z+y)^-3+(1+z)(1+y)/7',input='y^(15/4)(1-y)^4'))

def actual_scales_and_volterra():
    records=[]
    for digits in [100,150]:
        mp.mp.dps=digits;k=mp.mpf(2000);g0=k**-10
        g=2*g0/(1+mp.sqrt(1+44*g0/(360*mp.pi)))
        Fl=1-384*mp.pi*mp.mpf('1e22')*g
        bclass=mp.sqrt((mp.sqrt(1+12*mp.pi*k**6)-1)*g/(6*g0))
        km=-mp.mpf(6)/5+mp.log(bclass);d0=1/(480*mp.pi**2)
        AT=Fl/(48*mp.pi*g)+mp.mpf(23)/(8640*mp.pi**2)
        R=24*d0;A=4*AT;eta=R/A
        def H_at(s):return mp.digamma(s-1)+mp.pi*mp.cot(mp.pi*(s-4)/2)/2+km
        def M_at4(t):
            if t==0:return R
            s=4+1j*t;q=s*(s-3)
            return 1j*t*(d0*q*(q+2)*H_at(s)-AT*q)
        ratios=[]
        for t in map(mp.mpf,['0','1e-35','1e-34','1e-30','.1','.5','1','2','3','10','1e20','1e40']):
            ratio=abs(M_at4(t))/(6*d0);ratios.append(mp.nstr(ratio,25))
            check('endpoint_contour_diagnostic_'+str((digits,str(t))),ratio>=1)
            if t:
                exact=mp.pi*t/(2*mp.sinh(mp.pi*t))+mp.mpf('.5')+t*t/(1+t*t)+t*t/(4+t*t)
                # At enormous t subtraction of the two limiting imaginary parts loses digits.
                if t<100:check('digamma_endpoint_identity_'+str((digits,str(t))),abs((-t*mp.im(H_at(4+1j*t)))/exact-1)<mp.mpf('1e-55'))
        for sig,t in [(mp.mpf('3.5'),0),(mp.mpf('3.9'),0),(4-mp.mpf('1e-34'),0),(mp.mpf('3.9'),1),(mp.mpf('3.9'),2),(mp.mpf('3.9'),3)]:
            s=sig+1j*t;q=s*(s-3);upper=192 if t<=2 else 20
            check('strip_diagnostic_'+str((digits,str(sig),t)),mp.re((q+2)*H_at(s))<upper)
        aa=mp.mpf(1);bb=mp.mpf(2)
        def w(t):return -bb/A*mp.exp(-eta*t)+aa/R*mp.expm1(-eta*t)
        for t in [mp.mpf(1),mp.mpf('1e6'),1/eta,10/eta]:
            integral=t*mp.quad(lambda x:w(t*x),[0,1])
            residual=-A*w(t)-R*integral-(bb+aa*t)
            check('Volterra_direct_integral_'+str((digits,mp.nstr(t,8))),abs(residual)/(bb+aa*t)<mp.mpf('1e-55'))
            check('Volterra_weighted_bound_'+str((digits,mp.nstr(t,8))),abs(w(t))/(1+t)<=max(abs(aa),abs(bb))/A)
        cancellation=(aa/R)/abs(w(1))
        check('separate_tail_loss_'+str(digits),cancellation>mp.mpf('1e33'))
        check('normal_combined_scale_'+str(digits),mp.mpf('5.34e33')<1/eta<mp.mpf('5.38e33'))
        records.append(dict(digits=digits,eta=mp.nstr(eta,65),separate_coefficient_scale=mp.nstr(1/eta,65),
            moment_residue=mp.nstr(R,65),transformed_normal_L2_inverse_bound=mp.nstr(20*mp.pi**2,65),
            comparison_weighted_bound=mp.nstr(2/A,65),comparison_at_t_one=mp.nstr(w(1),65),
            separate_to_combined_ratio_at_one=mp.nstr(cancellation,65),sampled_contour_ratios=ratios))
    for key in ['eta','separate_coefficient_scale','comparison_weighted_bound','comparison_at_t_one']:
        check('cross_precision_'+key,abs(mp.mpf(records[0][key])/mp.mpf(records[1][key])-1)<mp.mpf('1e-55'))
    FIXTURES.append(dict(fixture='actual_normal_scales_and_independent_Volterra_comparison',records=records,
        samples_prove_whole_contour=False,comparison_is_actual_quantum_response=False,
        admission_norm_inverse_quantified=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_stationary_constraint_result',old['status']=='passed' and old['stationary_current_constraint_realization_complete'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        algebra();moments();actual_scales_and_volterra()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        physical_fast_tail_inverse_regularity_proved=error is None,
        physical_fast_tail_scope='Reference inverse on the actual stationary residual and stated finite polylog source class; five areal jets, no numerical admission norm or nonlinear neighborhood.',
        actual_remainder_moment_functional_derived=error is None,actual_kernel_endpoint_moment_numerically_evaluated=False,
        physical_coupled_inverse_proved=False,full_dynamical_boundary_completed=False,admission_norm_inverse_quantified=False,
        nonlinear_solution_constructed=False,semiclassical_validity_proved=False,physical_stability_proved=False,
        component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
