"""Whole correlated reference trial in the actual retained energy norm."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction
import json,sys,argparse,traceback,ast,inspect
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_coupled_last_filter_test as last
num=last.num;old=last.old;lr=last.lr;fr=last.fr
arb,ctx,dec,enc,acb=last.arb,last.ctx,last.dec,last.enc,last.acb
up=last.up;cs=last.cs;ROOT=last.ROOT
qball,qenc,middle,cmiddle,cball,cenc=num.qball,num.qenc,num.middle,num.cmiddle,num.cball,num.cenc
polyadd,polymul,polyeval,polypoles=num.polyadd,num.polymul,num.polyeval,num.polypoles
STEM='current-correlated-signed-trial-test'
INPUTS=list(dict.fromkeys(['notes/'+STEM+'.md','checks/current_correlated_signed_trial_test.py',
 'evidence/current-coupled-last-filter-test.json','evidence/current-coupled-last-filter-test-verification.json',
 'checks/current_coupled_last_filter_review.py','notes/current-signed-adjoint-trial-attempt.md',
 'evidence/current-signed-adjoint-trial-attempt.json','evidence/current-correlated-signed-trial-test.json',
 'notes/current-signed-observation-regularization-test.md',
 'notes/current-spatial-inverse-output-test.md','notes/current-filtered-mixed-green-residual.md',
 *last.INPUTS]))
digest=lambda p:sha256(p.read_bytes()).hexdigest()
STAGES={}

def rational_trial(p,c):
 prior=lr.result('current-coupled-last-filter-test')['nominal_profile_parameters']
 Q=lambda q:Fraction(int(q['numerator']),int(q['denominator']))
 v={k:Q(x) for k,x in prior.items()}
 alpha=cs('current-filtered-mixed-green-residual')['alpha_original_radius_cost']
 v.update({k:middle(x) for k,x in dict(alpha=alpha,gamma=p['gamma'],A=c['classical_A'],
  muY=c['muY'],scalar_c=c['scalar_c'],twice_a=2*p['a']).items()})
 lin=[v['Cv'],v['delta']];quad=[v['d0'],v['d1'],Fraction(1)]
 n1=[v['K0']*x for x in polymul(polypoles([v['Cv']*v['lambda_n'],v['twice_a']]),quad)]
 cj=-3*v['alpha']*v['kappa']/(4*v['gamma']*v['A']*v['At'])
 n2=[cj*x for x in polymul(lin,polypoles([v['muY'],v['muS']+2*v['scalar_c'],v['lambda_l'],v['lambda_n']]))]
 joint=polyadd(n1,n2)
 N=[v['kappa']**2*x for x in polymul(polymul(lin,lin),joint)]
 D=polypoles([v[k] for k in ['muT','muS','lambda_l','lambda_n','twice_a']])
 for _ in range(3):D=polymul(D,quad)
 assert len(N)==8 and len(D)==12 and all(x>0 for x in D)
 return v,(N,D),dict(common_signed_bracket=[qenc(x) for x in joint],
  numerator=[qenc(x) for x in N],denominator=[qenc(x) for x in D],
  j_reference_coefficient=qenc(cj),first_bracket=[qenc(x) for x in n1],second_bracket=[qenc(x) for x in n2])


def contour_representation(p,nominal,nd,count=256):
 Nn,Dn=nd;B=[Fraction(1,4),Fraction(1,2),Fraction(1),Fraction(3)]
 Bpoly=polypoles(B);Gnum=polymul(Bpoly,Nn);ginf=Gnum[-1]
 ln=nominal['lambda_n'];other=polymul(polypoles([nominal[k] for k in ['muT','muS','lambda_l','twice_a']]),[nominal['d0'],nominal['d1'],1])
 other=polymul(other,polymul([nominal['d0'],nominal['d1'],1],[nominal['d0'],nominal['d1'],1]))
 Rn=polyeval(Gnum,-ln)/polyeval(other,-ln)
 field=dict(amplitude=qenc(nominal['At']),base_poles=[qenc(v) for v in B],
  constant=qenc(ginf),large_pole=qenc(ln),large_pole_residue=qenc(Rn),upper_half_nodes=[])
 Gn=[qball(v) for v in Gnum];Dnball=[qball(v) for v in Dn]
 circles=[(Fraction(2),Fraction(1,4)),(Fraction(4,3),Fraction(1,5)),(Fraction(2,3),Fraction(1,10))]
 qcost=rcost=arb(0);details=[]
 h=qball(nominal['d1'])/2;e=qball(nominal['d0'])-h*h
 for center,radius in circles:
  C,R=qball(center),qball(radius)
  inside=[]
  for key in ['muT','muS','lambda_l','lambda_n','twice_a']:
   dist=abs(qball(nominal[key])-C)
   assert dist<R/2 or dist>3*R/2
   if dist<R/2:inside.append(key)
  dh=abs(C-h)
  assert dh+abs(e).sqrt()<R/2 or dh-abs(e).sqrt()>3*R/2
  bounds=[]
  for rr in [R/2,3*R/2]:
   mindeg=arb(min(arb(1).lower(),(C-rr).lower()));assert mindeg>0
   den=arb(1)
   for key in ['muT','muS','lambda_l','lambda_n','twice_a']:
    dist=abs(qball(nominal[key])-C);gap=rr-dist if dist<rr else dist-rr
    assert gap>0;den*=gap
   gap=rr-dh if dh<rr else dh-rr
   dg=gap*gap-abs(e);assert dg>0;den*=dg**3
   numerator=sum((up(v)*(C+rr)**j for j,v in enumerate(Gn)),arb(0))
   bounds.append(up(rr*numerator/den*256/(3*mindeg)))
  quad=bounds[1]/((arb(3)/2)**count-1)+bounds[0]/(arb(2)**count-1)
  qcost+=quad;round_circle=arb(0)
  for j in range(count//2):
   theta=2*arb.pi()*(arb(j)+arb('.5'))/count
   pole=acb(C)+R*acb(theta.cos(),theta.sin())
   weight=-(pole-C)*polyeval(Gn,-pole)/polyeval(Dnball,-pole)/count
   pn,wn=cmiddle(pole),cmiddle(weight);pb,wb=cball(pn),cball(wn)
   dmu=up(abs(pole-pb));dw=up(abs(weight-wb))
   mm=arb(min(arb(1).lower(),pole.real.lower(),pb.real.lower()));assert mm>0
   term=256*dw/(3*mm)+256*up(abs(wb))*dmu/(arb(13).sqrt()*mm**2)
   round_circle+=2*term
   field['upper_half_nodes'].append(dict(pole=cenc(pn),weight=cenc(wn),circle_center=qenc(center)))
  rcost+=round_circle
  details.append(dict(center=qenc(center),radius=qenc(radius),enclosed_linear_poles=inside,
   annular_norm_bounds=[enc(v) for v in bounds],quadrature_X_before_scale=enc(quad),rounding_X_before_scale=enc(round_circle)))
 scale=qball(nominal['At'])*p['gamma'].sqrt()
 return field,up(scale*qcost),up(scale*rcost),details

def reuse_reference_trial():
 raw=json.loads((ROOT/'evidence/current-correlated-signed-trial-test.json').read_bytes())
 assert raw['status']=='specified_complete_joint_reference_certificate_rejected'
 changed={'notes/'+STEM+'.md','checks/current_correlated_signed_trial_test.py'}
 assert all(digest(ROOT/p)==h for p,h in raw['input_sha256_after'].items() if p not in changed)
 tree=ast.parse(raw['arithmetic_source_utf8'])
 for fun in [rational_trial,contour_representation]:
  node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==fun.__name__)
  assert ast.get_source_segment(raw['arithmetic_source_utf8'],node).strip()==inspect.getsource(fun).strip()
 return raw['result']

def averaged_reference_error(p,v,CE):
 lin=[v['Cv'],v['delta']];quad=[v['d0'],v['d1'],Fraction(1)]
 P=polymul(polymul(polymul(lin,lin),lin),polypoles([v['muY'],v['muS']+2*v['scalar_c']]))
 coeff=-3*v['alpha']*v['kappa']**3/(4*v['gamma']*v['A']*v['At'])
 P=[coeff*x for x in P]
 D=polypoles([v[k] for k in ['twice_a','muS','muT']])
 for _ in range(3):D=polymul(D,quad)
 from math import comb
 assert len(P)==6 and len(D)==10 and all(x>0 for x in D)
 g=min(x/Fraction(comb(9,j)) for j,x in enumerate(D));M=sum(abs(x) for x in P)
 XH=qball(v['At'])*(8*p['gamma']/15).sqrt()*qball(M/g)
 error=up(CE*arb('1e-110')*XH/(2*p['gamma']))
 return error,dict(numerator=[qenc(x) for x in P],denominator=[qenc(x) for x in D],
  positive_Bernstein_guard=qenc(g),numerator_coefficient_sum=qenc(M),
  complete_H_profile_X_upper=enc(XH),actual_paired_energy_error=enc(error),
  prescribed_averaging_radius=enc(arb('1e-110')))

def execute():
 ctx.prec=768;p,c,bs=fr.load();lr.model_rows(p,c);c=old.setup(p,c)
 nominal,nd,polydata=rational_trial(p,c)
 STAGES['polynomial_trial']=polydata
 CE=cs('current-filter-numerical-representation-attempt3')['original_X_to_Et_comparison']
 reused=reuse_reference_trial()
 field=reused['numerical_field'];contours=reused['contour_bounds']
 oldc={k:dec(v) for k,v in reused['constants'].items()}
 Xquad,Xround=oldc['contour_X_error'],oldc['node_weight_X_error']
 # Reuse the original outward error for the identical field. Reapplying
 # up() to serialized intermediate upper endpoints can enlarge it slightly
 # and need not overlap its previously stored, much narrower endpoint ball.
 assert CE.overlaps(oldc['original_X_to_Et_comparison'])
 erepr=oldc['complete_representation_energy_error']
 center={k:dec(v) for k,v in reused['point_center'].items()}
 queries=reused['whole_cell_queries'];lower=oldc['joint_trial_energy_lower']
 avgerr,average=averaged_reference_error(p,nominal,CE)
 STAGES['averaged_reference_comparison']=average
 lower=arb((lower-avgerr).lower())
 print(json.dumps(dict(averaged_reference_error=str(avgerr),whole_trial_lower=str(lower))),flush=True)
 # Retain the EXACT full quartic contribution through a rigorous norm upper.
 oldtrial=lr.result('current-signed-adjoint-trial-attempt')['known_complete_certificate_debits']
 vcap=dec(oldtrial['psi_retained_energy_upper'])
 alpha=cs('current-filtered-mixed-green-residual')['alpha_original_radius_cost']
 eps=cs('current-filtered-mixed-green-residual')['epsilon'];qt=arb('.338')/(1-eps)
 f=arb('.338')/arb('1.338')
 quartic=up(alpha*f*f*qt*vcap)
 full_lower=arb((lower-quartic).lower())
 pref=cs('current-filtered-mixed-green-residual')['pref_t'];bolt=arb('1.667e-9')
 previous=cs('current-coupled-last-filter-test')
 known=previous['known_cumulative_error_subtotal_NOT_complete_bolt']
 physical=full_lower*pref
 gap=full_lower-bolt/pref
 # The entire known subtotal is a conservative overcharge for the first
 # contribution's existing error. It is never treated as unused margin.
 decisive_gap=gap-known/pref
 passed=full_lower>bolt/pref
 values=dict(original_X_to_Et_comparison=CE,contour_X_error=Xquad,node_weight_X_error=Xround,
  complete_representation_energy_error=erepr,averaged_reference_energy_error=avgerr,joint_trial_energy_lower=lower,
  quartic_complete_energy_upper=quartic,complete_trial_energy_lower=full_lower,
  original_bolt_trial_report_lower=physical,whole_bolt_allowance=bolt,
  trial_report_over_whole_allowance=physical/bolt,minimum_center_change_in_energy=gap,
  known_error_subtotal_NOT_margin=known,conditional_signed_allowance_other_errors_zero=(bolt-known)/pref,
  known_error_subtotal_divided_by_pref=known/pref,
  decisive_remaining_j_comparison_threshold=decisive_gap,
  proposed_half_gap_j_error_contract=decisive_gap/2)
 return ('specified_complete_joint_reference_certificate_rejected' if passed else
         'joint_reference_trial_lower_test_inconclusive'),dict(constants={k:enc(v) for k,v in values.items()},
  exact_nominal_parameters={k:qenc(v) for k,v in nominal.items()},polynomial_trial=polydata,
  numerical_field=field,contour_bounds=contours,point_center={k:enc(v) for k,v in center.items()},
  whole_cell_queries=queries,averaged_reference_comparison=average,
  reused_point_reference_receipt="evidence/current-correlated-signed-trial-test.json",full_actual_j_response=None,actual_Zr=None,
  scope=dict(actual_paired_domain_trial=True,complete_quantum_form_used_in_lower_comparison=True,
   all_spectrum_representation_error=True,whole_cell_not_pointwise_lower=True,
   full_quartic_term_retained=True,prescribed_averaged_reference_error_included=True,
   specified_reference_certificate_rejected=bool(passed),
   actual_j_Riesz_or_filtered_response_evaluated=False,actual_signed_Zr_bounded=False,
   actual_bolt_screen_rejected=False,original_admission=False,nonlinear_existence=False,
   exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,
   construction_obstruction=False,goal_complete=False))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--receipt',default=STEM);args=ap.parse_args()
 out=ROOT/'evidence'/(args.receipt+'.json');assert not out.exists(),'Preserve all receipts'
 before={p:digest(ROOT/p) for p in INPUTS};protected=json.loads(lr.BASELINE.read_bytes())
 assert all(digest(Path(p))==h for p,h in protected.items())
 rt=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes());pkg=rt['packages']['flint']
 assert digest(Path(sys.executable))==rt['python_sha256']
 assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
 receipt=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  ordinary_argument_utf8=(ROOT/'notes'/(STEM+'.md')).read_text(encoding='utf-8'),
  arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),python_sha256=rt['python_sha256'])
 try:
  receipt['status'],receipt['result']=execute();after={p:digest(ROOT/p) for p in INPUTS};assert after==before
  assert all(digest(Path(p))==h for p,h in protected.items())
  receipt.update(input_sha256_after=after,protected_paths_unchanged=len(protected),pinned_flint_files_unchanged=len(pkg['files']))
 except Exception:receipt.update(status='failed',traceback=traceback.format_exc(),partial_results=STAGES)
 receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
 out.write_bytes((json.dumps(receipt,indent=2)+'\n').encode('utf-8'))
 print(json.dumps(dict(status=receipt['status'],error=receipt.get('traceback'),
  constants={k:v['ball'] for k,v in receipt.get('result',{}).get('constants',{}).items()}),indent=2),flush=True)
 if receipt['status']=='failed':raise SystemExit(1)
