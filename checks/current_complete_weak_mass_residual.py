"""Complete weak residual and linear mass; original graph and nonlinear goal open."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys,traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
from flint import arb,ctx
import flint
import sympy as sp
s=sp

ROOT=Path(__file__).resolve().parents[1]
STEM='current-complete-weak-mass-residual'
OWN=['notes/'+STEM+'.md','checks/current_complete_weak_mass_residual.py',
     'evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/current-reference-clock-adjoint-trial.md',
 'notes/current-coupled-mass-adjoint.md',
 'notes/current-complete-mass-moment-bound.md',
 'notes/current-joint-adjoint-residual-attempt.md',
 'notes/current-geodesic-two-sided-remainder.md',
 'notes/stationary-curved-quantum-tt.md',
 'notes/stationary-reference-cold-bolt.md',
 'notes/stationary-seed-quantum-form.md',
 'notes/current-seed-full-hessian-realization.md',
 'notes/current-global-return-energy-attempt.md',
 'notes/current-finite-local-action-assembly.md',
 'notes/current-full-domain-contact-assembly.md',
 'notes/current-controlled-solution-audit-closeout.md',
 'evidence/current-global-mass-conversion.json',
 'evidence/current-reference-clock-adjoint-trial.json',
 'evidence/current-coupled-global-return-attempt.json',
 'evidence/current-complete-mass-moment-bound.json']
CHECKS=[]

def check(name,value):
    CHECKS.append({'name':name,'passed':bool(value)})
    if not value:raise AssertionError(name)

def enc(v):
    m,r,e=v.mid_rad_10exp(50)
    return {'ball':v.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}

def dec(v):
    d=v['decimal_enclosure']
    return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])

def prior(stem):
    r=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
    check('accepted_input_'+stem,r['status']=='passed')
    return r['result']

def local_rows():
    z=s.symbols('z',positive=True)
    ph=s.symbols('ph0:5');b=s.symbols('b0:5')
    be,uu,aa,bb,ep,C=s.symbols('be uu aa bb ep C',real=True)
    al,kv,dd,ar,h,q=s.symbols('al kv dd ar h q',real=True)
    Q=be+uu*z*z-aa*z**3+bb*z**4-2*ep*z**5+ep*z**6
    def D(f):return s.diff(f,z)+sum(s.diff(f,ph[j])*ph[j+1]+s.diff(f,b[j])*b[j+1] for j in range(4))
    Lap=lambda p:Q*p[2]+s.diff(Q,z)*p[1]
    L=-s.diff(Q,z,2)-2*Lap(b);M=Lap(ph);N=Q*ph[1]**2;Di=L-6*M-6*N
    Wd=s.exp(-2*ph[0])*(2+C*s.exp(-2*b[0])*Di)+12*C*be
    sq=al*s.exp(4*ph[0]+2*b[0])*Wd**2
    pot=(kv*s.exp(4*ph[0])-dd*(4*s.exp(-ph[0])-3*s.exp(-2*ph[0]))+ar)*s.exp(2*b[0])
    we=h*(C*C*s.exp(-2*b[0])*(ph[0]*L**2/s.Integer(60)+L*M/30-M*M/10+L*N/30-7*M*N/90-7*N*N/180)-C*ph[0]*L/18+(s.Rational(1,15)-q*q/3)*ph[0]*s.exp(2*b[0])-C*N/18)
    seed={ph[0]:-s.log(z),**{ph[j]:s.diff(-s.log(z),z,j) for j in range(1,5)},**{v:0 for v in b}}
    out={}
    for name,expr in [('square',sq),('potential',pot),('weyl',we)]:
     for label,vars in [('nu',ph),('shape',b)]:
      grad=s.diff(expr,vars[0])-D(s.diff(expr,vars[1]))+D(D(s.diff(expr,vars[2])))
      f=s.expand(grad.subs(seed))
      out[name+'_'+label]=str(f)
    return out

def execute():
    ctx.prec=512
    p=prior('current-global-mass-conversion')['parameters'];beta,uu,eps,tau,g,Fl=[dec(p[k]) for k in ['beta','u','eps','tau','g','Fl']]
    tr=prior('current-reference-clock-adjoint-trial');ret=prior('current-coupled-global-return-attempt')
    ga=dec(ret['gamma']);h=dec(ret['h']);zh=1/(1+h);C=1/uu
    z=s.symbols('z',positive=True);d=s.symbols('d',real=True)
    be,us,aa,bb,ep,cs=s.symbols('be uu aa bb ep C',real=True)
    Q=be+us*z*z-aa*z**3+bb*z**4-2*ep*z**5+ep*z**6
    av=4*beta+2*uu+tau;bv=3*beta+uu+eps
    vals={be:beta,us:uu,aa:av,bb:bv,ep:eps,cs:C}
    def ev(ex,vs):
     if ex.is_Number:return arb(str(ex.p))/int(ex.q) if ex.is_Rational else arb(str(ex))
     if ex in vs:return vs[ex]
     if ex.is_Add:return sum((ev(x,vs) for x in ex.args),arb(0))
     if ex.is_Mul:
      o=arb(1)
      for x in ex.args:o*=ev(x,vs)
      return o
     if ex.is_Pow:return ev(ex.base,vs)**int(ex.exp)
     raise ValueError(ex)
    qt=[arb(0),2*ga]+[ev(s.diff(Q,z,k),{**vals,z:zh})*((-1)**k)/int(s.factorial(k)) for k in range(2,7)]
    check('positive_actual_temperature_and_bolt',ga>0 and ga<arb('1e-30') and zh>arb('.99') and zh<1)
    check('actual_quadratic_bolt_coefficient',qt[2]>arb('.999999') and qt[2]<arb('1.000001'))
    check('actual_cubic_bolt_coefficient',qt[3]>arb('-1.334') and qt[3]<arb('-1.333'))
    check('actual_quartic_bolt_coefficient',qt[4]>arb('.5') and qt[4]<arb('.501'))
    check('actual_degree_five_six_terms_retained',abs(qt[5])<arb('1e-27') and abs(qt[6])<arb('1e-27'))
    check('whole_real_P_lower',arb('.999999')-arb('1.334')+arb('.5')-arb('2e-27')>arb('.16'))
    check('whole_real_P_decreasing',-arb('1.333')+2*arb('.501')+arb('7e-27')<0)
    check('whole_real_P_upper',qt[2]+arb('2e-27')<arb('1.001'))
    check('retained_Q_real_global_upper',beta+uu+abs(eps)<arb('.17'))
    check('outer_tenth_Q_lower',beta-av*arb('.1')**3-2*abs(eps)>arb('.16'))
    check('exact_beta_parameter_relation_interval_consistency',(6*beta+uu+eps-1).contains(0))
    check('middle_relative_Q_comparison',12*(uu+abs(eps))+6*tau/arb('1e-16')<arb('1e-8'))
    check('middle_Qprime_comparison',(arb(22)/3)*uu+(arb(28)/3)*abs(eps)+3*tau/arb('1e-8')<arb('1e-9'))
    check('middle_clock_comparison_error',(arb('3e-9')+arb('6e-8'))/(1-arb('1e-8'))<arb('1e-7'))
    def istar(v):
     return -4*(1-v).log()/3+2*(3*v*v+2*v+1).log()/3+7*arb(2).sqrt()/6*((3*v+1)/arb(2).sqrt()).atan()-7*arb(2).sqrt()/6*(1/arb(2).sqrt()).atan()+1/(1-v)-1
    Ip=1/(1-z)-1-s.Rational(4,3)*s.log(1-z)+s.Rational(2,3)*s.log(1+2*z+3*z*z)+7*s.sqrt(2)/6*(s.atan((3*z+1)/s.sqrt(2))-s.atan(1/s.sqrt(2)))
    Qp=(1-z)**2*(1+2*z+3*z*z)/6
    check('exact_comparison_primitive',s.simplify(s.diff(Ip,z)-1/Qp)==0 and s.simplify(Ip.subs(z,0))==0)
    Hprime=s.diff(s.Rational(3,2)/(z*z*(1-z))-Ip,z)
    check('whole_interval_critical_polynomial',s.simplify(Hprime*z**3*(1-z)**2*(1+2*z+3*z*z)-s.Rational(3,2)*(5*z**3-z-2))==0)
    Avec=s.Function('A')(z);mp=s.symbols('mp')
    xiz=Q*Avec;bjoint=1+s.diff(Q,z)*Avec/2+mp*z**3/Q
    check('local_density_Ward_factor',s.simplify((s.diff(xiz,z)+1-2*bjoint).subs(s.diff(Avec,z),(Q+2*mp*z**3)/Q**2))==0)
    # unique critical point in [.827,.829] for 5z^3-z-2; use a wider interval only if proved below.
    rr=arb('.827', '.002')
    Hmin=arb('1.5')/(rr*rr*(1-rr))-istar(rr)
    check('unique_H_critical_root_bracket',5*arb('.825')**3-arb('.825')-2<0 and 5*arb('.829')**3-arb('.829')-2>0)
    check('whole_interval_clock_comparator',Hmin>3)
    zc=1-arb('1e-8');dc=zh-zc
    ac=qt[2];DA=ac*zh**2+3*ga*zh;massp=ga**2/(zh*DA);cc=zh**2/(2*DA)
    Uc=sum((qt[k]*dc**(k-1) for k in range(1,7)),arb(0))
    V0=ac+6*massp*zh-cc*(ac*ac+4*ga*qt[3]);K=V0/ac
    Im=istar(zc)/(1-arb('1e-8'));Il=istar(zc)/(1+arb('1e-8'))
    Ibound=arb.union(Il,Im)
    Rc=Ibound-cc/dc-K/Uc
    Jcor=128*massp/(3*dc**3)
    check('middle_correlated_mass_not_dropped',100*massp/dc**2<arb('1e-48'))
    check('cold_anchor_with_integral_comparison',abs(Rc)+Jcor<arb('1e6'))
    check('cold_pole_removed_anchor',abs(Ibound-cc/dc)+Jcor<arb('1e9'))
    gs,zs=s.symbols('gs zs',positive=True);qq=s.symbols('q2:7',real=True)
    Us=2*gs+sum(qq[k-2]*d**(k-1) for k in range(2,7))
    DAs=qq[0]*zs*zs+3*gs*zs;ps=gs*gs/(zs*DAs);ccs=zs*zs/(2*DAs)
    Vn=s.Poly(s.expand(d*Us+2*ps*(zs-d)**3-ccs*Us*Us),d)
    assert s.cancel(Vn.nth(0))==s.cancel(Vn.nth(1))==0
    Vs=sum(s.cancel(Vn.nth(k))*d**(k-2) for k in range(2,Vn.degree()+1))
    Wn=s.Poly(s.expand(-ccs*(Us+d*s.diff(Us,d))*Us/2+ps*(zs-d)**3),d)
    assert s.cancel(Wn.nth(0))==0
    Ws=sum(s.cancel(Wn.nth(k))*d**(k-1) for k in range(1,Wn.degree()+1))
    Ks=s.cancel(Vs.subs(d,0)/qq[0])
    Ns=s.Poly(s.expand(Vs-Ks*s.diff(Us,d)),d)
    assert s.cancel(Ns.nth(0))==0
    N1=sum(s.cancel(Ns.nth(k))*d**(k-1) for k in range(1,Ns.degree()+1))
    Ds=s.Poly(s.expand(Us-Ks*(Us+d*s.diff(Us,d))/2+Ws),d)
    vs={gs:ga,zs:zh,**{qq[k-2]:qt[k] for k in range(2,7)}}
    def polynorm(ex,radius):return sum((abs(ev(s.cancel(co),vs))*radius**degree for (degree,),co in s.Poly(ex,d).terms()),arb(0))
    vn=polynorm(Vs,arb('2e-8'));n1=polynorm(N1,arb('2e-8'))
    wn=sum((abs(ev(s.cancel(Ws.coeff(d,k)/(gs if k==0 else 1)),vs))*arb('2e-8')**max(k-1,0) for k in range(s.degree(Ws,d)+1)),arb(0))
    dn=sum((abs(ev(s.cancel(Ds.nth(k)/(gs**2 if k==0 else gs if k==1 else 1)),vs))*arb('2e-8')**max(k-2,0) for k in range(Ds.degree()+1)),arb(0))
    check('cold_V_uniform_polynomial',vn<100)
    check('cold_N1_uniform_polynomial',n1<arb('1e4'))
    check('cold_W_uniform_polynomial',wn<100)
    check('cold_Dnum_factor_and_bound',dn<5000)
    cold_b=arb('1e4')*(ga+arb('1e-8'))+arb('1.5')*((ga+arb('1e-8'))*arb('1e6')+arb('4e4')*(2*ga+arb('1e-8'))/arb(1).exp())
    check('entire_regular_bolt_bounded',cold_b<arb('.02'))
    check('near_outer_n_bound',arb('.17')/arb('.16')+arb('1e-40')<100)
    check('middle_n_bound',9*arb('1.001')/arb('.1')+18*ga/arb('.1')<100)
    check('cold_n_bound',3*(ga+arb('1e-8'))*(1+arb('1e-8')*arb('1e9')+400)/arb('.99')<100)
    Rchart=arb(1)/32
    check('complex_Q_second_derivative_cone',2*arb('1.000001')+6*arb('1.334')*arb('1.02')+12*arb('.501')*arb('1.02')**2+arb('1e-20')<20)
    check('complex_Q_third_derivative_cone',6*arb('1.334')+24*arb('.501')*arb('1.02')+arb('1e-20')<30)
    check('Goursat_relative_displacement',arb(2).sqrt()*Rchart+6*arb('.052').exp()*Rchart**2/4<arb('.047'))
    check('Goursat_absolute_displacement',arb('.17').sqrt()*Rchart+arb('.7')*arb('.052').exp()*Rchart**2/4<arb('.014'))
    check('Goursat_sigma_cone',20*arb('.052').exp()*Rchart**2/8<arb('.003'))
    check('Goursat_Picard_rows',32*arb('.052').exp()*Rchart**2/4<arb('.02') and 71*arb('.052').exp()*Rchart**2/8<arb('.02'))
    check('complex_middle_P_real_lower',arb('.16')-arb('2.36')*arb('.02')>arb('.11'))
    check('complex_middle_Q_lower',arb('.874')*(arb('.874')*arb('.11')-arb('.126')*arb('.05'))>arb(1)/32)
    check('complex_middle_A_increment',arb('.126')*32+arb('1e-30')<5)
    check('complex_middle_b_bound',1+6*14/2+arb('1e-30')<50)
    check('complex_cold_U_lower',2-arb('1.01')/8>arb('.8') and arb('.98')-arb('1.01')/8>arb('.8'))
    check('complex_cold_A_increment',100/(8*arb('.8')**2)<20)
    check('complex_cold_W_quotient',100*arb('1.125')/arb('.8')<142)
    check('complex_cold_b_bound',1+arb('1.5')*(arb('1e9')*(ga+arb('1e-8'))+420)+142<800)
    jets={k:arb(2)**k*int(s.factorial(k))*32**k*arb('1e6') for k in range(1,5)}
    K1,K2,K3,K4=map(arb,['1e8','1e11','1e14','1e16'])
    check('joint_geodesic_first_jet',jets[1]<K1)
    check('joint_geodesic_second_jet',jets[2]<K2)
    check('joint_geodesic_third_jet',jets[3]+15*jets[1]<K3)
    check('joint_geodesic_fourth_jet',jets[4]+60*jets[2]+arb('4e6')*jets[1]<K4)

    # Local source rows, with exact constant-curvature and vacuum cancellations.
    exprs=local_rows()
    al,kv,dd,ar,hs,q=s.symbols('al kv dd ar h q',real=True)
    loc={}
    for key,st in exprs.items():loc[key]=s.sympify(st,locals={str(v):v for v in [z,be,us,aa,bb,ep,cs,al,kv,dd,ar,hs,q]})
    hv=g/(2*arb.pi()*C*C*beta);alph=8*arb.pi()*arb('1e22')*g/(C*C*beta)
    lv={be:beta,cs:C,aa:av,bb:bv,ep:eps,al:alph,dd:Fl*eps,ar:arb(0,'6e-23'),hs:hv,q:arb(2000)**4}
    N0=arb(100);B0=arb('1.000001')
    def norm(poly):
     out=arb(0);ll=s.symbols('ll')
     pp=s.Poly(s.expand(poly.subs(s.log(z),ll)*z**4),z,ll)
     for (j,k),coeff in pp.terms():
      power=j-4
      assert power>=-2,(power,k,coeff)
      moment=arb(int(s.factorial(2*k)))/(2*power+5)**(2*k+1)
      out+=abs(ev(coeff,lv))*moment.sqrt()
     return out
    local_bound={}
    for label in ['nu','shape']:
     f=s.expand(sum(loc[name+'_'+label] for name in ['square','potential','weyl']).subs({us:1/cs,kv:hs*11*cs*cs*be*be/60}))
     nn=norm(f);nd=norm(z*s.diff(f,z));nm=norm(s.diff(f,aa))
     r=(2*B0*nn+N0*nd+2*massp*nm)/beta
     local_bound[label]={'source':nn,'radial':nd,'mass_parameter':nm,'residual':r}

    check('area_coefficient_uncertainty_retained',Fl<1 and arb('3.4e-22')*beta<arb('6e-23'))
    local_total=arb('12.65')*local_bound['nu']['residual']+arb('17.6')*local_bound['shape']['residual']
    check('complete_both_local_rows',local_total<arb('1e-18'))
    pi=arb.pi();etaG=C/dec(tr['B_squared_over_G']);alpha=2*etaG/beta
    Lq=dec(tr['dimensionless_complete_source_mass_factor']);Flmin=1-arb('1.179e-8')
    T1=(pi*arb('1e5')*beta/(arb('.99')*Flmin)).sqrt()
    j0=1000*etaG/beta;jo=arb('16e-8')/25*j0;j1=arb('2e5')*etaG
    T2sq=150*pi*pi+((3*pi/2)*500**4*jo+(pi/2)*j1)/(arb('.99')*Flmin)
    check('physical_fourth_joint_reserve',T2sq<39**2)
    interp=(1/(525*pi*pi))**arb('.125')*(1/(90*pi*pi))**arb('.875')*36*arb(18)**arb('-.25')
    check('complete_angular_fractional_multiplier',interp<arb('.02'))
    check('complete_principal_multiplier_with_common_fourth',interp+100/(90*pi*pi)<2)
    h0=arb(1)/32
    check('single_seam_crossing_short_time',h0<pi/arb(2).sqrt())
    check('seam_normalized_flux_reserve',8*beta.sqrt()<32)
    check('fractional_difference_coefficient_less_than_one',8/(3*pi.sqrt())>1)
    Bfrac=((8*pi/3)*K4*K4*h0**arb('1.5')+512*K3*K3*h0.sqrt()+16*pi*K3*K3/h0.sqrt()).sqrt()
    epsmix=arb('1e-24')
    anchor=arb(1)/400;small=(epsmix/(4*arb('1e13'))).sqrt()
    check('full_geodesic_tail_gradient_reserve',2*pi*36*128**2<arb('1e7'))
    gradient_raw=arb('1e13')*(arb('.25')+2*(anchor/small).log())+arb('1e13')+arb('1e7')
    check('complete_gain_two_gradient_reserve',gradient_raw<arb('1e15'))
    def gain_product(e):
     return alpha*((2*pi)*(e*K2*K2/3+arb('1e15')*K1*K1)).sqrt()*(e*39**2/(3*alpha)+arb('1e15')*T1*T1).sqrt()
    old_product=gain_product(arb('1e-6'))
    check('old_epsilon_valid_product_does_not_certify_target',old_product>arb('3.27169e-5'))
    new_product=gain_product(epsmix)
    check('corrected_combined_norm_product_fits',new_product<arb('4.054e-6'))
    costs={
     'principal':2*alpha*((2*pi).sqrt()*K1+Bfrac)*T1,
     'gain_two_complete_mixed_product':new_product,
     'actual_row_subtracted_gluing':arb('17.6')*alpha*arb('6e24')*B0,
     'actual_complete_constant_contact':arb('35.2')*Lq*arb('1e20')*B0*(arb(1)/5+1/(2*arb(5).sqrt())+arb(1)/15).sqrt(),
     'local_both_rows_reserve':arb('1e-18'),'LLL_reserve':arb('1e-17')}
    LLL=alpha*(arb(2000)**4/C)*K1*T1/(12*pi*pi.sqrt())
    check('full_fixed_domain_LLL',LLL<costs['LLL_reserve'])
    residual=sum(costs.values(),arb(0))
    check('COMPLETE_RESIDUAL_CERTIFIED',residual<arb('2.1e-5') and residual<arb('3.27169e-5'))
    moment=prior('current-complete-mass-moment-bound')
    source=dec(moment['complete_reference_adjoint_pairing_upper'])
    mass=(source+arb('1.1e-10')*residual)/arb('1.9999999764')
    remaining=arb('1e-11')-2174*mass
    check('FULL_COUPLED_LINEAR_MASS_TARGET_CLOSED',mass<arb('2.29e-15'))
    check('sufficient_remaining_lifted_allowance',remaining>arb('6.46e-12'))
    return {'precision_bits':ctx.prec,'joint_supremum':enc(B0),'joint_n_supremum':enc(N0),
     'actual_bolt_polynomial_coefficients':[enc(v) for v in qt[2:]],
     'whole_interval_comparator_H_lower':enc(Hmin),'pole_removed_cold_anchor':enc(Rc),
     'retained_mass_integral_error':enc(Jcor),'cold_supremum_bound':enc(cold_b),
     'cold_polynomial_guards':{'V':enc(vn),'N1':enc(n1),'W_scaled':enc(wn),'Dnum_scaled':enc(dn)},
     'geodesic_jets':{str(k):enc(v) for k,v in enumerate([K1,K2,K3,K4],1)},
     'fractional_trial_bound':enc(Bfrac),'complete_angular_interpolation_coefficient':enc(interp),
     'physical_gradient_test_constant':enc(T1),'physical_fourth_test_constant_squared':enc(T2sq),
     'mixed_estimator_review':{'epsilon':enc(epsmix),'gradient_coefficient_upper':enc(gradient_raw),
       'old_epsilon_1e_minus_6_corrected_product':enc(old_product),
       'separate_derivative_product_draft_valid':False,'combined_norm_product_used':True},
     'local_row_norms':{k:{j:enc(v) for j,v in a.items()} for k,a in local_bound.items()},
     'complete_local_residual_upper':enc(local_total),'LLL_residual_upper':enc(LLL),
     'complete_residual_costs':{k:enc(v) for k,v in costs.items()},
     'complete_adjoint_residual_upper':enc(residual),'full_coupled_linear_mass_upper':enc(mass),
     'sufficient_remaining_lifted_graph_allowance':enc(remaining),
     'scope':{'complete_weak_adjoint_residual_certified':True,'complete_source_pairing_retained':True,
     'full_coupled_LINEAR_mass_target_closed':True,'actual_kappa_m_computed':False,
     'original_coupled_graph_certified':False,'independent_phase_columns_bounded':False,
     'nonlinear_mass_certified':False,'nonlinear_admission':False,'exact_compatibility':False,
     'physical_matching':False,'validity_noise':False,'all_sector_stability':False,
     'actual_controlled_wormhole':False,'construction_obstruction':False,'full_objective_complete':False}}

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
    message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),'result':run.get('result'),'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
