"""Full source enclosure and fractional-reference coefficient; inverse still open."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb,ctx,dec,enc
from flint import acb
from current_localized_first_residual_range_test import pa,ps,pm,pp,pv,mag

ROOT=Path(__file__).resolve().parents[1]
STEM='current-source-adapted-geodesic-first-attempt'
BASELINE=Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
INPUTS=['notes/'+STEM+'.md','checks/current_source_adapted_geodesic_first_attempt.py',
 'checks/current_quartic_ritz_rank_test.py','checks/current_localized_first_residual_range_test.py',
 'notes/current-first-quantum-solve-attempt.md','notes/current-finite-local-action-assembly.md',
 'notes/current-geodesic-two-sided-remainder.md','notes/current-seed-full-hessian-realization.md',
 'notes/current-scalar-green-source-residual.md','notes/current-quartic-ritz-rank-test.md',
 'evidence/current-global-green-trial-enclosure.json','evidence/current-global-mass-conversion.json',
 'evidence/current-scalar-green-source-residual.json','evidence/current-first-inverse-source-screen.json',
 'evidence/current-first-quantum-solve-attempt.json','evidence/current-positive-time-quantum-operator.json',
 'evidence/current-localized-first-residual-range-test-verification.json',
 'evidence/current-coupled-mass-adjoint-runtime.json']
def digest(p):return sha256(p.read_bytes()).hexdigest()
def result(name):
 d=json.loads((ROOT/'evidence'/(name+'.json')).read_bytes());r=d.get('runs',[d])[-1]
 assert r['status']!='failed';return r['result']
def constants(name):return {k:dec(v) for k,v in result(name)['constants'].items()}

def execute():
 ctx.prec=512
 p={k:dec(v) for k,v in result('current-global-green-trial-enclosure')['trial_parameters'].items()}
 L,a,ga,beta,d,C,B2=[p[k] for k in ['L','a','gamma','beta','d','C','B2']]
 At=dec(result('current-global-green-trial-enclosure')['constants']['A_t'])
 scalar=result('current-scalar-green-source-residual');sc=constants('current-scalar-green-source-residual')
 fs=constants('current-first-inverse-source-screen');plus=constants('current-positive-time-quantum-operator')
 yprof=scalar['profiles']['y0'];ks=[dec(v['degree']) for v in yprof]
 coeff=[dec(v['coefficient']) for v in yprof];lam=[a*k*(k+1) for k in ks]
 assert sum(coeff,arb(0)).contains(0) and sum((c*l for c,l in zip(coeff,lam)),arb(0)).contains(0)
 assert arb('1.01')**2*(1-arb('1e-70'))*arb('.99')*(1-arb('1.179e-8'))>1
 Cxi=arb('1.01');S=arb('1e-12');zm=L-2*S;M3=arb(1000)
 rr=M3*2*S/(6*a);Q2S=2*ga*S+a*S*S;Q22S=4*ga*S+4*a*S*S
 Qmax=(1+rr)*Q22S;Qsmax=2*ga+4*a*S+M3*(2*S)**2/2
 epsbar=C/(60*arb.pi()*B2);Cv=1+11*epsbar*beta/3
 w0=arb('1.1e-27')*beta;taumax=arb('3.01')*beta
 KA=(2*d+11*epsbar/6)*w0+(22*epsbar/3)*taumax
 gram=[]
 for i,k in enumerate(ks):
  row=[]
  for j,h in enumerate(ks):
   if i==j:g=acb(k+1).polygamma(1).real/(2*k+1)
   else:g=((k+1).digamma()-(h+1).digamma())/(k*(k+1)-h*(h+1))
   row.append(ga/a*g)
  gram.append(row)
 def quad(x,y):return sum((x[i]*gram[i][j]*y[j] for i in range(4) for j in range(4)),arb(0))
 dcoeff=[c*l for c,l in zip(coeff,lam)]
 n2=quad(coeff,coeff);g2=-quad(coeff,dcoeff);dd2=quad(dcoeff,dcoeff)
 assert n2>0 and g2>0 and dd2>0
 Ny0,Ng0,Nd0=n2.sqrt(),g2.sqrt(),dd2.sqrt()
 old={k:dec(v) for k,v in scalar['profile_bounds']['y0'].items()}
 cut=S.sqrt()*(4*(1+rr)*old['U']/S+(15*Qmax/S**2+2*(Qsmax+2*Qmax/zm)/S)*old['T'])
 Ny=Ny0/zm**2
 Ng=(1+rr).sqrt()/zm*(Ng0+2*(Q22S/S).sqrt()*old['T'])
 Nd=(1+rr)*Nd0+(5*M3/(6*a)+2*(1+rr)/zm)*old['H']+cut
 Es=( (6*(d*Nd**2+(Cv+4*d*beta)*Ng**2+4*beta*Cv*Ny**2)+KA*Ng**2+2*w0*Ny**2)/beta ).sqrt()
 oldEs=fs['scalar_source_uniform_upper']*fs['actual_trial_Hmu_norm_upper']+sc['full_scalar_residual_Es_dual_upper']
 Erem=Cxi*sc['complete_mixed_remainder_constant']*Es/beta
 oldErem=Cxi*sc['complete_mixed_remainder_constant']*oldEs/beta
 xS=a*S/ga;X=2*xS
 B2moment=(arb(3).log()**2/3+2*arb(3).log()/9+arb(2)/27)/4
 tail_y=arb(0);tail_dy=arb(0);Mdy=arb(0)
 for k,c,l in zip(ks,coeff,lam):
  tail=arb.pi()*S.sqrt()*abs(c)*xS**(-k-1)/(2*k+1).sqrt()
  tail_y+=tail;tail_dy+=l*tail
  tint=1/(2*k-1) if k>arb('.5') else (X**(1-2*k)-1)/(1-2*k)
  Mdy+=abs(c)*l*((ga/a)**3*(B2moment+arb.pi()**2*tint)).sqrt()
 CE=2*L+L*L*M3/(6*a);CB=5*L*L*M3/(6*a)+2*L*(1+rr)
 ndiff=CE*Mdy+CB*old['H']+cut+L*L*tail_dy
 T0=sc['tau0']/L**4
 Efinite=4*Cxi/beta*(T0*(Cv*tail_y+d*ndiff)+14*abs(p['epsilon_seed'])*(Cv*Ny+d*Nd))
 Esource=Erem+Efinite
 # Exact polynomial numerator of the classical comparison, retaining R2.
 den=[2*a,arb(1)];den2=pm(den,den);v=[-2*a,-arb(3)/2]
 P=[arb(0),arb(1)];D=pm(P,pa(ps(v,2),ps(den,3)))
 E=pa(ps(pm(pa(ps(P,arb(3)/2),[-a]),den2),1/L**2),
      ps(pm(pm(P,v),den),2/L**2),
      ps(pp(pa(den,ps(v,2)),2),1/(C*L**2)),
      ps(pp(pa(den,v),2),12*arb(2000)**-10*B2/(C*L**4)))
 ZoverC=arb.pi()*arb(2000)**16/(C*B2)
 EM=ps(pm(v,v),-4*ZoverC)
 alphaR=8*arb.pi()*arb('1e22')/(C*B2)
 R2=ps(pa(ps(pm(v,v),16),ps(pm(P,den2),-12*C),
    ps(pa(ps(pm(v,v),16*a*a),ps(pm(v,D),16*a),ps(pm(D,D),2),ps(pm(P,den2),12*a)),C*C)),alphaR)
 num=ps(pa(E,EM,R2),-1/beta)
 # The D numerator is exactly 2a*p; remove no uncertain numerical terms.
 assert D[0].is_zero() and D[1].overlaps(2*a) and all(t.contains(0) for t in D[2:])
 assert all(t.contains(0) for t in num[4:])
 A=num[3];B=num[2]-4*a*A
 rem=pa(num,ps(pm([B,A],den2),-1))
 assert all(t.contains(0) for t in rem[2:])
 assert A>0 and B>0
 muT=B/A;muY=4*a/3;muS=4*beta/L**2
 lambda_l=a*p['ell']*(p['ell']+1);lambda_n=1/(d*L**2)
 Kx=4*sc['tau0']/(beta*L**4*A);c=2*sc['tau0']/(3*L**2)
 K0=Kx*c*(1-d*L**2*lambda_l)
 assert K0>0 and Cv>1 and lambda_n>lambda_l
 m=arb(min(arb(1).lower(),muT.lower(),muS.lower(),lambda_l.lower()))
 Nmodel=K0*Cv/m**(arb(5)/2)
 scale=At*ga.sqrt()*Nmodel
 Emodel=Cxi*scale*(abs(rem[1])/(8*a)+abs(rem[0])/(4*a*a))
 # Complete untruncated fractional geodesic coefficient.
 order=arb(9)/2;mu=order+2
 def moment(nu):
  return 2**(mu-3)*(mu/2).gamma()**2*(mu/2+nu).gamma()*(mu/2-nu).gamma()/mu.gamma()
 cg=arb(720)**(-arb(1)/4)*arb(24)**(-arb(3)/4)
 cstar=cg*(moment(1)-moment(0))/(2*arb.pi()**2)
 qa=arb(2000)**4
 Sfrac=C**(arb(5)/4)*(2*(qa+2)/(qa+1)**(arb(5)/4)+8/(qa+1)**(arb(1)/4))
 Kgeo=4*arb.pi()*cstar*Sfrac*plus['CL_t']/(B2*beta)
 coef=Kgeo*scale;target=arb('9.955e32')
 residual_after_known=target-Esource-Emodel
 constants_out=dict(
  actual_scalar_model_L2=Ny0,actual_scalar_model_gradient_L2=Ng0,actual_scalar_model_Delta_L2=Nd0,
  actual_cutoff_scalar_Es_upper=Es,old_uniform_scalar_Es_upper=oldEs,
  actual_full_mixed_remainder_debit=Erem,old_full_mixed_remainder_debit=oldErem,
  scalar_cutoff_and_finite_geometry_source_debit=Efinite,complete_actual_source_debit=Esource,
  source_debit_over_first_stage_screen=Esource/target,
  old_mixed_debit_over_first_stage_screen=oldErem/target,
  norm_gain=oldEs/Es,cutoff_Delta_cost=cut,whole_source_Delta_difference=ndiff,
  classical_A=A,classical_B=B,muT=muT,muY=muY,muS=muS,
  classical_r0=rem[0],classical_r1=rem[1],complete_model_density_remainder=Emodel,
  K0=K0,Cv=Cv,large_mass_lambda=lambda_n,
  model_fractional_profile_upper=Nmodel,original_model_fractional_scale=scale,
  actual_bolt_rho4_log_coefficient=-At*K0/64,
  fractional_geodesic_kernel_coefficient=cstar,complete_angular_fractional_upper=Sfrac,
  actual_reference_fractional_output_coefficient=Kgeo,
  full_first_residual_cost_per_Cflow=coef,
  known_first_stage_subtotal=Esource+Emodel,
  remainder_after_known_NOT_available_margin=residual_after_known,
  conditional_Cflow_ceiling_RETURN_AND_NUM_ZERO=residual_after_known/coef,
  first_residual_to_bolt_multiplier=fs['first_residual_to_original_bolt_multiplier'],
  known_source_bolt_estimator_debit=fs['first_residual_to_original_bolt_multiplier']*Esource)
 assert all(v.is_finite() for v in constants_out.values())
 status='source_adapted_full_source_enclosed_actual_reference_inverse_unresolved'
 if Esource+Emodel>=target:status='source_adapted_known_estimator_excessive_actual_inverse_unresolved'
 return status,dict(precision_bits=ctx.prec,constants={k:enc(v) for k,v in constants_out.items()},
  exact_scalar_Gram=[[enc(v) for v in row] for row in gram],
  classical_symbol_numerator=[enc(v) for v in num],
  source_profiles=[{'degree':enc(k),'coefficient':enc(c)} for k,c in zip(ks,coeff)],
  norm_conventions=dict(original_energy='E_t^2=ell_t/N',source_dual='E_t*',
   original_measure='dz=gamma dx retained',flow_norm='(integral_SM |P|^(9/2) lift(v) squared / (2pi))^1/2',
   model_norm='L2(dx) norm of (1+H)^(9/4)(v0/At); NOT the actual flow norm'),
  actual_flow_norm=None,actual_flow_comparison_constant=None,complete_finite_return_error=None,
  actual_numerical_application_error=None,complete_first_reference_residual=None,
  actual_complete_quantum_correction=None,signed_Gt=None,
  scope=dict(actual_source_error_certified=True,actual_seed_fractional_operator_coefficient=True,
   actual_H5_source_estimator_fails=True,actual_fractional_norm_quantified=False,
   reference_inverse_certified=False,actual_minus_reference_application_certified=False,
   complete_original_admission=False,nonlinear_existence=False,compatibility=False,
   physical_matching=False,semiclassical_validity_noise=False,all_sector_stability=False,
   construction_obstruction=False,goal_complete=False))

if __name__=='__main__':
 path=ROOT/'evidence'/(STEM+'.json');assert not path.exists(),'Preserve immutable evidence.'
 before={p:digest(ROOT/p) for p in INPUTS};protected=json.loads(BASELINE.read_bytes())
 assert all(digest(Path(p))==h for p,h in protected.items())
 rt=json.loads((ROOT/INPUTS[-1]).read_bytes());assert digest(Path(sys.executable))==rt['python_sha256']
 pkg=rt['packages']['flint'];assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
 rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
  arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),python_sha256=rt['python_sha256'])
 try:
  rec['status'],rec['result']=execute()
  after={p:digest(ROOT/p) for p in INPUTS};assert after==before
  assert all(digest(Path(p))==h for p,h in protected.items())
  rec.update(input_sha256_after=after,protected_paths_unchanged=len(protected))
 except Exception:
  import traceback
  rec.update(status='failed',traceback=traceback.format_exc())
 rec['finished_utc']=datetime.now(timezone.utc).isoformat()
 path.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
 print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
  constants={k:v['ball'] for k,v in rec.get('result',{}).get('constants',{}).items()}),indent=2))
 if rec['status']=='failed':raise SystemExit(1)
