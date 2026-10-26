"""Supporting arithmetic for constrained canonical interior elimination.

No signed complementary center, physical bolt bound, or mass-column certificate.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb,factorial
import sys,json,traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb,ctx
import sympy as S
ROOT=Path(__file__).resolve().parents[1];STEM='current-canonical-interior-schur-bound'
OWN=['notes/'+STEM+'.md','checks/current_canonical_interior_schur_bound.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+[
 'notes/current-seed-full-hessian-realization.md',
 'notes/current-physical-york-dual-assembly.md',
 'notes/current-coupled-input-energy-bound.md',
 'notes/current-coupled-schur-graph-attempt.md',
 'notes/current-global-mass-conversion.md','evidence/current-global-mass-conversion.json',
 'notes/current-global-far-return-energy-bound.md','evidence/current-global-far-return-energy-bound.json',
 'notes/current-coupled-global-return-attempt.md',
 'notes/current-original-graph-dirac-tail.md','evidence/current-original-graph-dirac-tail.json',
 'notes/current-mass-response-regulator-audit.md',
 'notes/current-rough-paired-variation.md','notes/current-third-metric-conversion.md',
 'notes/charged-conformal-response-and-constraints.md']
CHECKS=[]
def check(name,ok):
    assert name not in {x['name'] for x in CHECKS}
    CHECKS.append({'name':name,'passed':bool(ok)})
    if not ok:raise AssertionError(name)
def enc(x):
    m,r,e=x.mid_rad_10exp(55)
    return {'ball':str(x),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
    d=x['decimal_enclosure'];return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])
def prior(stem):
    r=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
    assert r['status']=='passed';return r['result']
def execute():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('ordinary_argument_frozen',sha256((ROOT/OWN[0]).read_bytes()).hexdigest()==review['scientific_sha256'][OWN[0]])
    p=prior('current-global-mass-conversion')['parameters']
    beta,u,eps,tau,Fl=[dec(p[k]) for k in ['beta','u','eps','tau','Fl']]
    Q={0:beta,2:u,3:-4*beta-2*u-tau,4:3*beta+u+eps,5:-2*eps,6:eps}
    Qplus={**Q,3:Q[3]+tau}
    def add(*args):
     out={}
     for d in args:
      for j,v in d.items():out[j]=out.get(j,arb(0))+v
     return out
    def scale(d,v):return {j:x*v for j,x in d.items()}
    def shift(d,n):return {j+n:x for j,x in d.items()}
    def mul(a,b):
     out={}
     for j,x in a.items():
      for k,y in b.items():out[j+k]=out.get(j+k,arb(0))+x*y
     return out
    def der(a):return {j-1:j*x for j,x in a.items() if j}
    def ev(a,z):return sum((x*z**j for j,x in a.items()),arb(0))
    def integral(a,lo,hi):
     total=arb(0)
     for j,x in a.items():
      total+=x*((hi/lo).log() if j==-1 else (hi**(j+1)-lo**(j+1))/(j+1))
     return total
    N=80;R=arb('.4');c=arb('.2');delta=arb('.0002');end=c+delta
    A={j-2:j*(j-5)*x for j,x in Q.items() if j>=2}
    Qlower=beta-sum((abs(x)*R**j for j,x in Q.items() if j),arb(0))
    Vmax=arb(2)/3*sum((abs(x)*R**j for j,x in A.items()),arb(0))/Qlower
    alpha=Vmax*R*R/14;Fmax=1/(1-alpha)
    invQ=[1/beta]
    for j in range(1,N):invQ.append(-sum((v*invQ[j-k] for k,v in Q.items() if k and k<=j),arb(0))/beta)
    V=[arb(2)/3*sum((v*invQ[j-k] for k,v in A.items() if k<=j),arb(0)) for j in range(N-1)]
    f=[arb(1),arb(0)]
    for j in range(2,N+1):f.append(sum((V[k]*f[j-2-k] for k in range(j-1)),arb(0))/(j*(j+5)))
    s={j+4:(j+1)*v/4 for j,v in enumerate(f)}
    xi={j+4:-(j+5)*v/2 for j,v in enumerate(f)}
    sp=der(s);spp=der(sp);xip=der(xi);Qp=der(Q)
    rr=arb('2.946e-9');Kup=arb('6.0001');Vup=arb('24.001')
    LS=add(mul(Q,spp),mul(add(Qp,scale(shift(Q,-1),-2)),sp))
    Da=scale(shift(LS,2),1/beta)
    as_int=shift(add(scale(mul(Da,Da),6*rr),scale(shift(mul(Qplus,mul(sp,sp)),2),Kup/beta),scale(mul(s,s),Vup)),-4)
    dxi=add(shift(xip,1),scale(xi,-2))
    lx_int=scale(shift(add(scale(mul(Q,mul(dxi,dxi)),arb('1.5')/beta),mul(xi,xi)),-4),arb('.99')*Fl)
    Apoly=integral(as_int,arb(0),c);Lpoly=integral(lx_int,arb(0),c)
    mom=[4,12,52]
    Ts=[Fmax/4*(N+6)**(j+1)*mom[j]*(c/R)**(N+1) for j in range(3)]
    Tx=[2*x for x in Ts]
    qm=arb('.17');qpm=arb('.3')
    laperr=(qm*Ts[2]+2*qm*Ts[1]+qpm*c*Ts[1])/beta
    Aerr=c**5/5*(6*rr*laperr**2+Kup*qm/beta*Ts[1]**2+Vup*Ts[0]**2)
    Lerr=arb('.99')*Fl*c**5/5*(arb('1.5')*qm/beta*(Tx[1]+2*Tx[0])**2+Tx[0]**2)
    Aout=(Apoly.sqrt()+Aerr.sqrt())**2;Llo=(Lpoly.sqrt()-Lerr.sqrt())**2
    sc=ev(s,c)+arb(0,(Ts[0]*c**4).abs_upper());spc=ev(sp,c)+arb(0,(Ts[1]*c**3).abs_upper())
    amp=spc+sc/c
    def junction(cells):
     out=arb(0)
     for i in range(cells):
      tt=arb(2*i+1,1)/(2*cells);z=c+delta*tt
      P=tt*(1-tt)**2;P1=1-4*tt+3*tt**2;P2=-4+6*tt
      sv=sc*c/z+amp*delta*P
      sv1=-sc*c/z**2+amp*P1
      sv2=2*sc*c/z**3+amp/delta*P2
      qz=ev(Q,z);qpz=ev(Qp,z);qplus=ev(Qplus,z)
      lap=z*z/beta*(qz*sv2+(qpz-2*qz/z)*sv1)
      term=z**-4*(6*rr*lap**2+Kup*qplus/beta*z*z*sv1**2+Vup*sv**2)
      out+=term*delta/cells
     return out
    Ah=junction(512)
    qpcomb={j:(4-j)*v for j,v in Q.items()}
    far_integrand=shift(add(scale(mul(qpcomb,qpcomb),6*rr/beta**2),scale(Qplus,Kup/beta),{0:Vup}),-6)
    Afar=(sc*c)**2*integral(far_integrand,end,arb(1))
    Aup=Aout+Ah+Afar;ratio=Aup/Llo
    margin=1-ratio-2*(arb('.338')*ratio).sqrt()
    conversion=1+ratio/2+(ratio+ratio*ratio/4).sqrt()
    inverse=conversion/margin
    normH=arb('1.338').sqrt()

    check('actual_seed_parameter_family',beta>arb('.1666') and beta<arb('.1667') and tau>0 and tau<arb('1e-65') and abs(eps)<arb('2e-29') and Fl>arb('.99999998'))
    check('complex_seed_Q_zero_free',Qlower>0)
    check('analytic_V_majorant',Vmax<arb('11.511'))
    check('Volterra_contraction',alpha<arb('.132'))
    check('analytic_f_bound',Fmax<arb('1.152'))
    check('Qplus_nonnegative_factor',beta-abs(eps)>0 and u>0)
    # On [0,.2], the negative cubic dominates the positive quartic and
    # tiny eps powers; beta+u*c^2 is a convenient real Q upper bound.
    check('outer_Q_upper',beta+u*c*c+4*abs(eps)<arb('.17'))
    check('outer_Qprime_upper',sum((abs(v)*c**j for j,v in Qp.items()),arb(0))<arb('.3'))
    epl=arb(2000)**-10/(60*arb.pi())
    check('scalar_gradient_upper',6+(2*rr+11*epl/6)*(12+arb('1e-27'))+22*epl*7/3<Kup)
    check('scalar_potential_upper',24+2*arb('1.1e-27')+88*epl<Vup)
    check('scalar_tail_energy_bound',Aerr.sqrt()<arb('3e-22'))
    check('TT_tail_energy_bound',Lerr.sqrt()<arb('3e-22'))
    check('TT_reverse_triangle_positive',Lpoly.sqrt()>Lerr.sqrt())
    check('whole_scalar_extension',Aup<arb('.000617'))
    check('outer_TT_lower',Llo>arb('.00279'))
    check('extension_ratio',ratio<arb('.221'))
    rg=arb('.221');dg=1-rg-2*(arb('.338')*rg).sqrt()
    cg=1+rg/2+(rg+rg*rg/4).sqrt()
    check('signed_test_coercive',dg>0)
    check('interior_inverse_original_energy',cg/dg<7 and inverse<7)
    graph=7*normH*arb('4.583e-13')
    quartic=7*normH*arb('7.949e-26')
    mass=7*normH*arb('7.949e-27')
    slow=7*normH*arb('1.141e-22')
    check('entire_outer_Schur_H3',graph<arb('3.712e-12'))
    check('entire_quartic_Schur',quartic<arb('6.437e-25'))
    check('entire_mass_moment_Schur',mass<arb('6.437e-26'))
    check('entire_slow_moment_Schur',slow<arb('9.240e-22'))
    check('entire_affine_interior_H3',graph*arb('1.1e-10')<arb('4.084e-22'))
    check('weighted_input_conversion',graph*31<arb('1.151e-10'))
    check('conditional_mass_energy_ceiling',arb('3.712e-12')*arb('1.3469e7')<arb('5e-5'))
    check('conditional_mass_N_ceiling',arb('1.151e-10')*arb('4.344e5')<arb('5e-5'))
    check('conditional_feedback_input_ceiling',arb('3.712e-12')*arb('3.251e8')<arb('.0012069'))
    check('conditional_energy_feedback_ceiling',arb('.0195239')*arb('.0618165')<arb('.0012069'))
    check('conditional_affine_residual_exact',S.Rational('5.34235e-12')-S.Rational('4.084e-22')==S.Rational('5.3423499995916e-12'))
    check('conditional_fast_scope_screen',arb('.028913')+arb('1e-11')*arb('1.0706e10')<arb('.1599164'))
    z,t,vv=S.symbols('z t vv');ff=S.Function('f')(z)
    yy=z**3*ff
    check('exterior_ODE_identity',S.expand((-S.diff(yy,z,2)+(6/z**2+vv)*yy)/(-z**3)-(S.diff(ff,z,2)+6*S.diff(ff,z)/z-vv*ff))==0)
    check('exact_scalar_tail',S.expand(z*z*S.diff(yy,z)/4-z*yy/2-(z**4*ff+z**5*S.diff(ff,z))/4)==0)
    check('exact_TT_tail',S.expand(-z*yy-z*z*S.diff(yy,z)/2+(5*z**4*ff+z**5*S.diff(ff,z))/2)==0)
    pp=t*(1-t)**2
    check('junction_value_slope_matching',pp.subs(t,0)==pp.subs(t,1)==0 and S.diff(pp,t).subs(t,0)==1 and S.diff(pp,t).subs(t,1)==0)
    k=S.symbols('k',integer=True,nonnegative=True)
    check('complete_series_tail_moments',all(S.summation((k+1)**(j+1)/2**k,(k,0,S.oo))==mom[j] for j in range(3)))
    qfun=S.Function('Q')(z);sf=S.Symbol('kappa')/z
    lap=S.expand(qfun*S.diff(sf,z,2)+(S.diff(qfun,z)-2*qfun/z)*S.diff(sf,z))
    check('exact_far_scalar_Laplacian',S.simplify(lap-S.Symbol('kappa')*(4*qfun-z*S.diff(qfun,z))/z**3)==0)
    check('indefinite_compression_counterexample',(S.Matrix([[1,1]])*S.diag(-1,1)*S.Matrix([1,1]))[0]==0)
    return {'precision_bits':ctx.prec,'parameters':p,'polynomial_degree':6,'series_order':N,'junction_cells':512,
     'analytic_Vmax':enc(Vmax),'Volterra_ratio':enc(alpha),'analytic_f_bound':enc(Fmax),
     'outer_scalar_energy_upper':enc(Aout),'junction_scalar_energy':enc(Ah),'far_scalar_energy':enc(Afar),
     'scalar_extension_upper':enc(Aup),'outer_TT_energy_lower':enc(Llo),'extension_ratio':enc(ratio),
     'signed_test_margin':enc(margin),'interior_inverse_upper_expression':enc(inverse),
     'rounded_ratio_inverse_expression':enc(cg/dg),'H_form_norm_upper_expression':enc(normH),
     'original_outer_Schur_H3_per_energy':enc(graph),'quartic_Schur_per_energy':enc(quartic),
     'mass_moment_Schur_per_energy':enc(mass),'exact_slow_moment_Schur_per_energy':enc(slow),
     'affine_interior_H3_upper':enc(graph*arb('1.1e-10')),'input_N_H3_coefficient':enc(graph*31),
     'scalar_tail_energy':enc(Aerr.sqrt()),'TT_tail_energy':enc(Lerr.sqrt()),
     'conditional_ceilings_not_actual_maps':{
      'mass_A_times_E':enc(arb('5e-5')/arb('3.712e-12')),
      'mass_A_times_N':enc(arb('5e-5')/arb('1.151e-10')),
      'feedback_A_times_JE':enc(arb('.0012069')/arb('3.712e-12')),
      'solution_energy_gamma_best_case':enc(arb('.0012069')/arb('.0195239'))},
     'scope':{'actual_paired_mass_log_domain':True,'common_regulator_full_angular_thermal_contacts_returns':True,
      'canonical_interior_energy_inverse_certified':True,'original_separated_outer_Schur_bound_certified':True,
      'original_outer_phase_observations_certified':True,'actual_input_map_norms_computed':False,
      'signed_complementary_centers_computed':False,'physical_bolt_or_strong_interior_inverse_computed':False,
      'complete_coupled_bound_certified':False,'actual_kappa_m_computed':False,
      'nonlinear_admission':False,'nonlinear_contraction':False,'exact_compatibility':False,
      'physical_matching':False,'semiclassical_validity_noise':False,'all_sector_stability':False,
      'actual_controlled_solution':False,'construction_obstruction':False,'full_goal_complete':False}}

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
      'frozen_sources_utf8':{f:(ROOT/f).read_text(encoding='utf-8') for f in OWN},
      'runtime':{'python':sys.version,'executable_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
                 'flint_version':flint.__version__,'sympy_version':S.__version__}}
    try:
        run['result']=execute();run['input_sha256_after']={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),
       'interior_inverse':run.get('result',{}).get('interior_inverse_upper_expression'),
       'outer_Schur':run.get('result',{}).get('original_outer_Schur_H3_per_energy'),
       'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
