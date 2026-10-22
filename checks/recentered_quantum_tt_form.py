"""Frozen diagnostics for the actual-reference quantum TT transfer argument."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as Q
import json,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='recentered-quantum-tt-form';PREVIOUS='recentered-canonical-tt-audit'
OWN=['notes/'+STEM+'.md','checks/recentered_quantum_tt_form.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def rational_transfer():
    k=2000;q=k**4;g0=Q(1,k**10);gl=(1-Q('1e-34'))*g0;c=Q('1e22')
    pl=Q('3.14159');ph=Q('3.14160');fl=1-384*ph*c*g0;cl=Q(49119000000);ch=Q(49200000000)
    check('vacuum_root_bracket',gl+11*gl**2/(360*pl)<g0 and g0+11*g0**2/(360*ph)>g0)
    check('C_lower_interval',cl**2<1+12*pl*gl*k**16)
    check('C_upper_interval',1+12*ph*g0*k**16/fl<ch**2)
    check('positive_Newton_factor',fl>1-Q('1.179e-8'))
    check('massive_gap',Q(q+1,1)/ch>18**2)
    check('angular_ratio_guard',ch/q<Q('.0035') and q/cl<400)
    xl=(cl-1)/6;xh=(ch-1)/6
    check('x_lower_guard',xl>8*k**3/8 and xl>Q('8e9'))
    check('C_over_x_guard',6+1/xl<Q('6.1'))
    check('physical_GC_over_B2',g0*(6+1/xl)<Q('8.5e-33'))
    check('physical_charge_curvature_ratio',g0*q*q/xl**2<Q('4e-27'))
    check('canonical_C_ratio_upper',1/fl<(1+Q('1e-8'))**2)
    check('canonical_C_ratio_lower',gl/g0>(1-Q('1e-33'))**2)
    C0min=6*k**3;Cratio=1+Q('1e-8')
    B2hi=(1+(Cratio-1)*C0min/(C0min-1))*g0/gl
    B2lo=1-Q('1e-33')*C0min/(C0min-1)
    check('fixed_scale_radius_ratio',B2hi<1+Q('1.01e-8') and B2lo>1-Q('2e-33'))
    # |log B/b_class|<1e-8 using log(1+x)<=x and -log(1-x)<=x/(1-x).
    check('fixed_scale_log_guard',(B2hi-1)/2<Q('1e-8') and (1-B2lo)/(2*B2lo)<Q('1e-8'))
    check('thermal_ratio_upper',B2hi/(1-Q('1e-33'))**2<Q('1.01')**2)
    check('thermal_ratio_lower',B2lo/(1+Q('1e-8'))**2>Q('.99')**2)
    hlo=Q(1,24*k**10);hhi=Q(3,k**10)
    check('actual_small_h_guard',hhi<Q('1e-30'))
    check('cold_log_guard',Q(27,10)**80>24*k**10)
    check('thermal_delta_ratio',(1+Q(2,3)*hhi+Q(1,6)*hhi*hhi)/(1+hhi)<Q('1.001'))
    delta=xh*0+hhi*hhi*(1+xh*((1+hhi)**2+2*(1+hhi)+3))/(1+hhi)
    check('thermal_delta_less_than_one',delta<1)
    check('RN_coefficient_a',(1+2*xl+Q(1,2))/(1+6*xl)<Q('.334'))
    check('RN_coefficient_Q',(1+3*xl)/(1+6*xl)<Q('.501'))
    check('RN_second_derivative_lower',6*Q('.334')**2<2)
    check('canonical_curvature_derivative_ledger',Q(1,5)*17**2<17**2 and Q(13,5)+Q(5,2)*17<56 and Q(1,5)*(Q(3,2)*5*13+17)**2<115**2)
    check('outer_physical_scale_guard',(6*k**3-1)/6*(gl/g0)>3**18)
    check('kappa4_lower_guard',2*Q('7.5')-1-Q('1e-8')>Q('13.9'))
    check('kappa4_upper_guard',2*8+Q(1,2)+Q('1e-8')<Q('17.5'))
    J0=1000*g0*(6+1/xl)**2;J1=200000*g0*(6+1/xl)
    check('actual_return_J0',J0<Q('5.2e-29'))
    check('actual_return_J1',J1<Q('1.7e-27'))
    check('outer_measure_J0',Q('5.2e-29')*Q(16,25)*Q(1,100)**4<Q('3.4e-37'))
    contact=Q(3*10**5*10**16,1)*Q('4e-27')/(4*pl)
    check('full_gradient_contact_budget',contact<Q('1e-6'))
    vacuum=21600*11*g0/(360*pl)
    check('changed_vacuum_Hessian_budget',vacuum<Q('1e-25'))
    # -log(mu0/sqrt(lambda))>9 and log 2<.7 yield the outer margin.
    outer=(9-Q(38,15)-Q(19,60)-4-Q(7,20))/(30*ph)
    check('outer_fourth_margin',outer>Q('.0175'))
    check('fourth_gain_two_loss',Q('1e-3')*(ph/2+Q(2,3))<Q('.00224'))
    check('positive_fourth_reserve',outer-Q('.00224')>1/(100*pl))
    loss=Q('.001')+Q('.002')+Q('.0006')+Q('2e-12')+Q('1e-25')+Q('1e-6')
    check('complete_negative_budget',loss<Q('.003601000003') and loss/fl<Q('.003602'))
    check('quantum_target_passes',loss/fl<Q('.005')<Q(1,4))
    check('complete_TT_lower_form',1-Q('.005')-Q('7.07e-8')>Q('.994'))
    FIXTURES.append(dict(fixture='rational_reference_and_loss_transfer',C_interval=[str(cl),str(ch)],F_lower=str(fl),h_interval=[str(hlo),str(hhi)],log_scale_absolute_bound='1e-8',J0_upper=str(J0),J1_upper=str(J1),quantum_loss_upper=str(loss/fl),physical_coupled_inverse_proved=False))

def identities():
    u,x,delta=S.symbols('u x delta',positive=True);C=1+6*x;a=1+2*x+delta/2;Q0=1+3*x
    BB=x+u*u-2*a*u**3+Q0*u**4
    check('exact_RN_factorization',S.expand(BB-((1-u)**2*(x*(1+2*u+3*u*u)+u*u)-delta*u**3))==0)
    check('exact_geodesic_seam',S.diff(BB,u).subs(u,0)==0)
    check('exact_second_derivative_vertex',S.simplify(S.diff(BB,u,2).subs(u,a/(2*Q0))-(2-3*a*a/Q0))==0)
    check('second_derivative_endpoint',S.simplify(S.diff(BB,u,2).subs(u,1)-(2*C-6*delta))==0)
    DY=lambda v:S.sqrt(BB/C)*S.diff(v,u)
    RR=-S.diff(BB,u,2)/C
    expected2=-(BB*S.diff(BB,u,4)+S.diff(BB,u)*S.diff(BB,u,3)/2)/C**2
    check('independent_canonical_second_curvature_jet',S.simplify(DY(DY(RR))-expected2)==0)
    expected3=-S.sqrt(BB/C)*(S.Rational(3,2)*S.diff(BB,u)*S.diff(BB,u,4)+S.diff(BB,u,2)*S.diff(BB,u,3)/2)/C**2
    # Common-denominator factorization reduces the exact radical identity;
    # simplify alone leaves an equivalent unreduced expression on this build.
    check('independent_canonical_third_curvature_jet',S.factor(S.together(DY(DY(DY(RR)))-expected3))==0)
    s,t,e,V,A0,A2=S.symbols('S T e V A0 A2')
    check('new_vacuum_disk_second_variation',S.simplify(S.diff(V*S.exp(e*(4*s+2*t)),e,2).subs(e,0)-4*V*(2*s+t)**2)==0)
    act=-A0*S.exp(2*e)/2+A2*S.exp(-2*e)/2
    check('complete_common_constant_restoration',S.diff(act,e,2).subs(e,0)==-2*A0+2*A2)
    n,q,z=S.symbols('n q z',positive=True);P=n*(n+q)
    f=2*z*(2*n+q)/(P*(P+z))
    check('positive_angular_difference',S.simplify(2*(2*n+q)*(1/P-1/(P+z))-f)==0)
    check('exact_integrated_angular_tail',S.simplify(S.diff(2*S.log(1+z/P),n)+f)==0)
    logder=S.diff(S.log(f),n)
    check('angular_tail_monotonicity_certificate',S.simplify(logder-(2/(2*n+q)-(2*n+q)/P-(2*n+q)/(P+z)))==0)
    check('positive_term_for_monotonicity',S.expand((2*n+q)**2-2*P)==2*n*n+2*n*q+q*q)
    FIXTURES.append(dict(fixture='RN_curvature_and_common_contact_identities',tail_formula='2 log(1+z/[N(N+q)])',old_reference_potential_reused=False))

def actual_parameters():
    records=[]
    for prec in [80,110]:
        mp.mp.dps=prec;k=mp.mpf(2000);q=k**4;g0=k**-10
        g=2*g0/(1+mp.sqrt(1+4*11*g0/(360*mp.pi)));F=1-384*mp.pi*mp.mpf('1e22')*g
        C0=mp.sqrt(1+12*mp.pi*k**6);x0=(C0-1)/6;b0=mp.sqrt(x0*g/g0)
        C=mp.sqrt(1+12*mp.pi*g*k**16/F);x=(C-1)/6;B=mp.sqrt(x)
        nu=2*q/k*g*C0/b0**3
        for zs in ['.25','.265','.28']:
            zi=mp.mpf(zs);kap=zi*nu;gam=kap*B/C
            KK=lambda h:((1+3*x)*((1+h)**-2+(1+h)**-3)+3*x*(1+(1+h)**-1))/(2*C)
            hh=gam
            for _ in range(5):hh-= (hh*KK(hh)-gam)/(KK(hh)+hh*mp.diff(KK,hh))
            de=hh*hh*(1+x*((1+hh)**2+2*(1+hh)+3))/(1+hh)
            # Direct horizon derivative in a cancellation-free polynomial identity.
            v=1+hh;kap_direct=hh*((v+1)*(1+3*x*(v*v+1)))/(2*B*v**3)
            check('actual_temperature_equation_'+zs+'_'+str(prec),abs(hh*KK(hh)/gam-1)<mp.mpf('1e-65'))
            check('independent_horizon_temperature_'+zs+'_'+str(prec),abs(kap_direct/kap-1)<mp.mpf('1e-65'))
            check('actual_cold_guards_'+zs+'_'+str(prec),1/(24*k**10)<hh<3/k**10 and mp.log(1/hh)<80)
            logscale=mp.log(B/b0);k4=mp.harmonic(q)/2-mp.euler-logscale
            check('actual_scale_and_spectral_guards_'+zs+'_'+str(prec),abs(logscale)<mp.mpf('1e-8') and mp.mpf('13.9')<k4<mp.mpf('17.5') and mp.sqrt((q+1)/C)>18)
            loss=mp.mpf('.003601000003')/F
            records.append(dict(digits=prec,z_i=zs,h=mp.nstr(hh,70),delta=mp.nstr(de,70),temperature=mp.nstr(kap,70),C=mp.nstr(C,70),M1=mp.nstr(mp.sqrt((q+1)/C),70),log_mu0B=mp.nstr(logscale,70),kappa4=mp.nstr(k4,70),quantum_loss_upper=mp.nstr(loss,70)))
    for zs in ['.25','.265','.28']:
        pp=[r for r in records if r['z_i']==zs]
        for key in ['h','temperature','C','M1','log_mu0B','kappa4','quantum_loss_upper']:
            check('precision_'+zs+'_'+key,abs(mp.mpf(pp[0][key])/mp.mpf(pp[1][key])-1)<mp.mpf('1e-60'))
    FIXTURES.append(dict(fixture='actual_preparation_and_new_reference',records=records,curved_quantum_kernel_sampled=False))

def angular_enclosures():
    records=[]
    for prec in [80,110]:
        mp.mp.dps=prec;q=mp.mpf(2000)**4;F0=4*mp.euler-2*mp.harmonic(q)
        for label,z in [('q/5',q/5),('q^2/2',q*q/2),('4q^2',4*q*q)]:
            root=mp.sqrt(q*q-4*z);np=(q+root)/2;nm=(q-root)/2
            F=-2*(mp.digamma(1+np)+mp.digamma(1+nm));D=mp.re(F0-F)
            check('digamma_real_'+label+'_'+str(prec),abs(mp.im(F))<mp.mpf('1e-65'))
            widths=[];entries=[]
            for N in [32,128]:
                partial=sum(2*z*(2*n+q)/(n*(n+q)*(n*(n+q)+z)) for n in range(1,N+1))
                lo=partial+2*mp.log1p(z/((N+1)*(N+1+q)))
                hi=partial+2*mp.log1p(z/(N*(N+q)))
                check('full_angular_enclosure_'+label+'_'+str(N)+'_'+str(prec),lo<D<hi)
                widths.append(hi-lo);entries.append(dict(N=N,lower=mp.nstr(lo,65),upper=mp.nstr(hi,65)))
            check('angular_tail_refinement_'+label+'_'+str(prec),widths[1]<widths[0]/3)
            records.append(dict(digits=prec,shift=label,complex_roots=z>q*q/4,full_difference=mp.nstr(D,65),enclosures=entries))
    for label in ['q/5','q^2/2','4q^2']:
        pp=[r for r in records if r['shift']==label]
        check('angular_precision_'+label,abs(mp.mpf(pp[0]['full_difference'])/mp.mpf(pp[1]['full_difference'])-1)<mp.mpf('1e-60'))
    FIXTURES.append(dict(fixture='independent_finite_sum_and_full_angular_resolvent',records=records,curved_spectral_data=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_local_TT_result',old['status']=='passed' and old['local_TT_relative_loss_bound']=='7.07e-8')
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        rational_transfer();identities();actual_parameters();angular_enclosures()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        quantum_TT_reference_lower_allowance='.005 E_eff' if error is None else None,
        full_TT_reference_lower_form='>.994 E_eff' if error is None else None,
        scope='smooth bolt-regular TT core and specified stronger-form completion at recentered RN',
        actual_form_closedness_in_own_norm_proved=False,physical_boundary_completed=False,
        physical_coupled_inverse_proved=False,coupled_correction_norm_proved=False,nonlinear_solution_constructed=False,
        semiclassical_validity_proved=False,physical_stability_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
