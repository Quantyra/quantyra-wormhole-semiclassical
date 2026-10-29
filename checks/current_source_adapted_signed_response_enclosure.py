"""Complete actual-trial response with spatially resolved source comparisons."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb
import json,sys,argparse,traceback
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_source_adapted_fractional_return as fr
lr=fr.lr;arb,ctx,dec,enc,J=fr.arb,fr.ctx,fr.dec,fr.enc,fr.J
ROOT=fr.ROOT;STEM='current-source-adapted-signed-response-enclosure'
INPUTS=list(dict.fromkeys(['notes/'+STEM+'.md','checks/current_source_adapted_signed_response_enclosure.py',
 'notes/current-source-adapted-complete-response-test.md',
 'evidence/current-source-adapted-complete-response-test.json',
 'evidence/current-source-adapted-complete-response-test-verification.json',
 'evidence/current-source-adapted-fractional-return.json',
 'evidence/current-source-adapted-fractional-return-verification.json',
 'notes/stationary-curved-quantum-tt.md','notes/stationary-centered-exterior-source.md',
 'notes/stationary-reference-cold-bolt.md','notes/current-seed-full-hessian-realization.md',
 *fr.INPUTS]))
def digest(p):return sha256(p.read_bytes()).hexdigest()
def up(v):return arb(abs(v).upper())
def minimum(*v):return arb(min(up(x).upper() for x in v))
def compute():
 ctx.prec=512;p,c,bs=fr.load();L,a,ga=[p[k] for k in ['L','a','gamma']]
 assert c['muY'].overlaps(4*a/3) # exact identity used in factored coefficient B
 sc=arb('.0004');assert L>sc and 1/(L-sc)<arb('1.001')
 profiles={n:lr.profile_components(p,c,n) for n in ['r','f']}
 scale=c['original_model_fractional_scale'];At=p['At']
 m=arb(min(arb(1).lower(),c['muT'].lower(),c['muS'].lower(),(a*p['ell']*(p['ell']+1)).lower()))
 facf=At*ga.sqrt()*c['K0']*c['Cv']
 Nf=up(facf/(10*m**5).sqrt());Gf=up(facf/(40*m**4).sqrt())
 mp={k:dec(v) for k,v in lr.result('current-global-mass-conversion')['parameters'].items()}
 beta=p['beta'];um,eps,tau=[mp[k] for k in ['u','eps','tau']]
 Q=[beta,arb(0),um,-4*beta-2*um-tau,3*beta+um+eps,-2*eps,eps]
 Dmax=sum((j*(j-1)*abs(5-j)*up(Q[j]) for j in range(2,7)),arb(0))
 edges=[arb(0),ga]
 while edges[-1]<sc:
  n=arb((2*edges[-1]).mid())
  if n>=sc:break
  edges.append(n)
 edges.append(sc)
 while edges[-1]+arb('.002')<L:edges.append(edges[-1]+arb('.002'))
 edges.append(L)
 queue=list(zip(edges,edges[1:]));cells=[];err2=arb(0);grad2=arb(0);tailerr2=arb(0)
 while queue:
  s0,s1=queue.pop(0);x0=s0/ga;x1=s1/ga
  x=J([arb((x0+x1)/2,(x1-x0)/2),arb(1)]+[arb(0)]*5)
  z=L-ga*x;q2=2*x+a*x*x
  dR=sum((bs[j]*ga**(j-2)*x**(j-1)/(2+a*x) for j in range(3,7)),0)
  R=1+dR;q=q2*R
  if not (R.c[0]>0 and R.c[0]<2):
   assert s1-s0>ga/1024
   mid=(s0+s1)/2;queue[0:0]=[(s0,mid),(mid,s1)];continue
  # Exact paired two-row cancellation, including all finite coefficients.
  A=-(ga*x/L)*sum(((z/L)**j for j in range(5)),0)+(z/L)**5*dR
  B=z**4/L**5*(2*z*(2*q2.d()*R.d()+q2*R.d().d())/3+8*ga*q.d()/3)
  CC=z**4/L**5*(5*ga*R-z*R.d());E=-2*z**4/(3*L**5)
  Dx=z*q.d().d().d()+3*ga*q.d().d()
  qr=up(q2.c[0]).sqrt()
  nr=minimum(lr.profile_bounds(p,c,profiles['r'],s0,s1),scale)
  nf=minimum(lr.profile_bounds(p,c,profiles['f'],s0,s1),scale)
  hr=minimum(lr.profile_bounds(p,c,profiles['r'],s0,s1,'flux'),qr*scale)
  wr=up(Dmax*(s1-s0).sqrt()*lr.profile_bounds(p,c,profiles['r'],s0,s1,'tail_integral'))
  flux=[]
  for n in ['r','f']:
   b=scale
   if s0>0:
    qlo=2*x0+a*x0*x0
    b=minimum(scale,lr.profile_bounds(p,c,profiles[n],s0,s1,'flux')/qlo.sqrt())
   flux.append(b)
  value=sum((up(co.c[0])*bb for co,bb in zip([B,CC,A,E],[nr,hr,nf,wr])),arb(0))
  gc=[qr*up((B.d()+(4*a/3)*CC+E*Dx).c[0]),up((B+q2*CC.d()).c[0]),
      qr*up((A.d()-CC).c[0]),up(A.c[0]),qr*up(E.d().c[0])]
  gradient=sum((co*bb for co,bb in zip(gc,[nr,flux[0],nf,flux[1],wr])),arb(0))
  assert value.is_finite() and gradient.is_finite()
  err2+=value*value;grad2+=gradient*gradient
  outside=not (s1<sc) # Include BOTH boundary-touching cells; Arb equality is not order.
  if outside:tailerr2+=value*value
  cells.append(dict(s0=enc(s0),s1=enc(s1),outside_cold=bool(outside),
    reconstruction_error_L2=enc(value),reconstruction_error_gradient=enc(gradient)))
  if len(cells)%150==0:print(json.dumps(dict(cells=len(cells),queued=len(queue))),flush=True)
 err=up(err2.sqrt());gerr=up(grad2.sqrt());U0=up(Nf+err);U1=up(arb(2).sqrt()*(Gf+gerr))
 tail=up(lr.profile_bounds(p,c,profiles['f'],sc,L)+tailerr2.sqrt())
 old=lr.result('current-source-adapted-complete-response-test')
 oc={k:dec(v) for k,v in old['constants'].items()};od={k:dec(v) for k,v in old['debits'].items()}
 U2,CV,CG,CD=[oc[k] for k in ['source_Hessian_upper','test_weighted_L2_constant','test_gradient_constant','test_Laplacian_constant']]
 C,B2,beta=[p[k] for k in ['C','B2','beta']];factor=2/(B2*beta);qang=arb(2000)**4;pi=arb.pi()
 K=arb('1e10');sp=arb('1e-24');b=(sp/(4*K)).sqrt();bg=K*(arb(1)/4+2*((arb(1)/400)/b).log())
 assert b<arb(1)/400 and K>4*arb('1e8')
 debits=dict(
  local_gain_two=factor*2*pi*C*(sp*U2**2/3+bg*U1**2).sqrt()*(sp*(3*CD**2/8+CG**2/4)/3+bg*CG**2).sqrt(),
  full_local_annular_kernel=factor*2*pi*arb('4e8')*C*(arb(400)/128).log()*U1*CG,
  finite_second_contact=factor*arb('1e13')*C*U1*CG,
  global_row_subtracted_gluing=factor*oc['complete_global_gluing_coefficient']*C*U0*CV,
  actual_constant_contact_cold_cap=factor*arb('2.1e8')*qang*qang/C*U0*CV,
  actual_constant_contact_OUTSIDE_cold=factor*arb('1e18')*qang*qang/C*(1+1/arb(5).sqrt())*tail*CV,
  actual_constant_contact_gluing=factor*oc['constant_contact_gluing_coefficient']*C*U0*CV,
  known_curvature_polynomial=factor*100*C*(U1*CG/(36*pi)+816*U0*CV/(960*pi)),
  shared_vacuum_uncertainty=od['shared_vacuum_uncertainty'],retained_energy_comparison=od['retained_energy_comparison'])
 quantum=sum(debits.values(),arb(0))
 ref=lr.constants('current-source-adapted-fractional-return')['complete_exact_analytic_first_reference_residual']
 first=lr.constants('current-first-quantum-solve-attempt');stage=first['conditional_common_stage_output_error_ceiling']
 actual=quantum+ref;account=arb('4.308597e32')
 Kfirst=c['first_residual_to_bolt_multiplier']
 vals=dict(model_L2_upper=Nf,model_gradient_upper=Gf,reconstruction_L2_error=err,
  reconstruction_gradient_error=gerr,actual_source_L2_upper=U0,actual_source_gradient_upper=U1,
  complete_outside_cold_source_L2_upper=tail,retained_actual_source_Hessian_upper=U2,
  source_L2_report_improvement=oc['source_L2_upper']/U0,source_gradient_report_improvement=oc['source_gradient_upper']/U1,
  complete_actual_minus_reference_response_upper=quantum,response_over_remaining_conditional_account=quantum/account,
  complete_actual_first_residual_upper=actual,actual_residual_over_screen=actual/arb('1e33'),
  actual_inverse_energy_error_upper=actual,conditional_stage_output_allowance=stage,
  remaining_UNEXECUTED_caller_Et_account=stage-actual,
  actual_first_energy_error_bolt_coefficient_product_NOT_complete_output=Kfirst*actual,
  complete_response_report_improvement=oc['complete_actual_minus_reference_upper']/quantum)
 assert all(v.is_finite() for v in [*debits.values(),*vals.values()])
 passed=actual<arb('1e33') and quantum<account and actual<stage
 return ('complete_exact_analytic_actual_first_residual_passes_caller_and_later_filters_open' if passed
  else 'complete_actual_fixed_trial_estimator_excessive'),dict(constants={k:enc(v) for k,v in vals.items()},
  debits={k:enc(v) for k,v in debits.items()},cell_count=len(cells),cells=cells,
  actual_signed_nonzero_center=None,numerical_caller_error=None,later_filter_application_errors=None,
  proposed_t_5_over_4_trace_theorem_proved=False,signed_Gt=None,
  scope=dict(actual_fixed_trial_response_enclosed=True,exact_analytic_actual_first_residual_passes=bool(passed),
   full_common_regulator_angular_thermal_contacts_returns_retained=True,actual_H5_source_assumed=False,
   numerical_first_inverse_applied=False,full_P32_or_Gt_certified=False,original_admission=False,
   nonlinear_existence=False,compatibility=False,physical_matching=False,validity_noise=False,
   all_sector_stability=False,construction_obstruction=False,goal_complete=False))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--receipt',default=STEM);args=ap.parse_args()
 out=ROOT/'evidence'/(args.receipt+'.json');assert not out.exists(),'Preserve immutable receipts'
 before={p:digest(ROOT/p) for p in INPUTS};protected=json.loads(lr.BASELINE.read_bytes())
 assert all(digest(Path(p))==h for p,h in protected.items())
 rt=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes());pkg=rt['packages']['flint']
 assert digest(Path(sys.executable))==rt['python_sha256']
 assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
 rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),
  python_sha256=rt['python_sha256'])
 try:
  rec['status'],rec['result']=compute();after={p:digest(ROOT/p) for p in INPUTS};assert after==before
  assert all(digest(Path(p))==h for p,h in protected.items())
  rec.update(input_sha256_after=after,protected_paths_unchanged=len(protected),pinned_flint_files_unchanged=len(pkg['files']))
 except Exception:rec.update(status='failed',traceback=traceback.format_exc())
 rec['finished_utc']=datetime.now(timezone.utc).isoformat()
 out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
 print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
  constants={k:v['ball'] for k,v in rec.get('result',{}).get('constants',{}).items()},
  debits={k:v['ball'] for k,v in rec.get('result',{}).get('debits',{}).items()}),indent=2),flush=True)
 if rec['status']=='failed':raise SystemExit(1)
