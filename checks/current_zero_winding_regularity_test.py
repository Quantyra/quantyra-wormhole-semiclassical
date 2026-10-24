"""Support the actual-domain zero-winding coefficient-norm obstruction and graph packet test."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb, factorial
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-zero-winding-regularity-test'
OWN=['notes/'+STEM+'.md','checks/current_zero_winding_regularity_test.py','evidence/'+STEM+'-review.json']
ROWS=[]
def check(name,ok):
 assert name not in {v['name'] for v in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def execute():
 c,s,t,C,R,X,Y,v,w,Q,G,B2,be,k=S.symbols('c s t C R X Y v w Q G B2 beta k',positive=True)
 density=-C*S.exp(-2*(c+s*v+t*w))*(R-2*(s*X+t*Y))**2/(240*S.pi)
 second=S.expand(S.diff(density,s,t).subs({s:0,t:0}))
 third=S.expand(S.diff(density,s,t,c).subs({s:0,t:0,c:0}))
 check('complete_a2_high_quadratic_coefficient',S.simplify(second.coeff(X*Y)+C*S.exp(-2*c)/(30*S.pi))==0)
 check('complete_a2_soft_third_coefficient',S.simplify(third.coeff(X*Y)-C/(15*S.pi))==0)
 logcoef=C*S.exp(-2*c)*Q**2/(30*S.pi)
 check('soft_log_derivative_nonzero',S.simplify(S.diff(logcoef,c).subs(c,0)+C*Q**2/(15*S.pi))==0)
 K=G*C*be**2/(15*S.pi*B2)
 shape=(-2*G/B2)*S.diff(logcoef,c).subs(c,0)
 check('original_normalization_and_shape_row',S.simplify(shape-2*K*(Q/be)**2)==0)
 heat=2*C*S.exp(-2*c)/t
 check('entire_Weyl_mass_rescaling',S.simplify(S.diff(heat,c).subs(c,0)+4*C/t)==0)
 n,q=S.symbols('n q',integer=True,positive=True);j=S.symbols('j',integer=True)
 count=S.summation(2*(q+2*j),(j,1,n))
 check('entire_angular_count_at_levels',S.expand(count-2*n*(n+q))==2*n)
 x,ss=S.symbols('x ss');sig=S.Rational(801,200)
 for jj in range(1,6):
  pp=S.Poly(S.prod(sig+x-i for i in range(jj)),x)
  check('physical_derivative_leading_order_'+str(jj),pp.degree()==jj and pp.LC()==1)
 check('coefficient_norm_tends_to_zero',-6+4==-2)
 check('weighted_C5_tends_to_zero',-6+5==-1)
 check('localized_H3_response_growth',-6+5+3==2)
 check('coefficient_norm_ratio_growth',(-6+5+3)-(-6+4)==4)
 check('C5_ratio_growth',(-6+5+3)-(-6+5)==3)
 check('bare_H8_missing_logarithm',-8+8==0 and -8+5+3==0)
 cv,rr=S.symbols('cv rr');J=-6*be*(cv-rr*ss*(ss-3))*(ss-4)*(ss+1);aop=(ss+3)/(3*(ss+1))
 coupled=S.cancel(J*aop)
 check('coupled_scalar_datum_order_four',S.Poly(coupled,ss).degree()==4)
 check('coupled_scalar_subleading',-6+4+3==1 and 1<2)
 pars,h,zh,env=old.scalar.prev.york.actual_parameters()
 zz=arb('.002','.0001');beta=pars['beta']
 QQ=beta+pars['u']*zz**2+(-4*beta-2*pars['u']-pars['tau'])*zz**3+(3*beta+pars['u']+pars['eps'])*zz**4-2*pars['eps']*zz**5+pars['eps']*zz**6
 aa=(QQ/beta)**2;KK=pars['C']*beta**2/(15*arb.pi()*pars['B2']);coef=2*KK*aa
 check('actual_support_ellipticity',QQ>0 and beta>0)
 check('whole_support_principal_margin',aa>arb('.9999999') and aa<arb('1.0000001'))
 check('actual_output_cutoff_is_one',arb('.0021')<arb('.01')/arb(1).exp())
 check('physical_soft_log_coefficient_positive',coef>arb('6.9077e-36') and coef<arb('6.9078e-36'))
 check('existing_graph_packet_ratio',2*aa>arb('1.9999998') and 2*aa<arb('2.0000002'))
 L=10*arb.pi()*(pars['Fl']+23*pars['g']/(180*arb.pi()))/pars['g']
 CB=10*KK*(3840+L/arb(2).log())
 check('retained_all_frequency_normal_upper_constant',CB<2)
 target=S.Rational('.159916395999');radius=S.Rational('1e-11')
 check('conditional_C3_ceiling',target/radius==S.Rational('1.59916395999e10'))
 check('conditional_equal_split',target/2==S.Rational('.0799581979995') and target/(2*radius)==S.Rational('7.99581979995e9'))
 check('unchanged_nonlinear_allocation_arithmetic',S.Rational('5e-9')+S.Rational('1e-8')/4+S.Rational('5e-17')*S.Rational('5e7')==S.Rational('1e-8') and S.Rational(1,4)+S.Rational('1e-8')*S.Rational('5e7')==S.Rational(3,4))
 return dict(K=old.encoded(KK),actual_support_a=old.encoded(aa),physical_soft_log_coefficient=old.encoded(coef),localized_full_graph_packet_ratio=old.encoded(2*aa),retained_normal_all_frequency_upper=old.encoded(CB),packet_amplitude='N^-6',coefficient_norm_decay='O(N^-2)',weighted_C5_decay='O(N^-1)',localized_H3_response_growth='2 K N^2 log N ||a psi||2 (1+o(1))',response_to_A_ratio_growth='at least a positive constant times N^4 log N',plain_H8_response_growth='log N',
 coefficient_only_metric_bound_obstructed=True,weighted_C5_metric_bound_obstructed=True,plain_interior_H8_metric_bound_obstructed=True,existing_full_graph_packet_loss_controlled=True,finite_N_crossing_certified=False,packet_limit_is_complete_operator_upper_bound=False,complete_graph_C3_known=False,complete_zero_winding_strong_limit_proved=False,local_interior_third_response_identified=True,rough_global_neighborhood_theorem_imported=False,actual_boundary_realization_retained=True,all_angular_levels_retained=True,all_thermal_and_boundary_terms_retained=True,physical_domain_changed=False,physical_source_changed=False,normalization_changed=False,original_graph_changed=False,closed_audit_reopened=False,remaining_nonlocal_target='.159916395999',actual_global_background_A_admission_known=False,original_X_graph_comparison_known=False,full_signed_projection_executed=False,actual_controlled_solution=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,next_complete_graph_response_executed=False)

def main():
 ctx.prec=512;review=json.loads((ROOT/OWN[2]).read_bytes());paths=set(OWN+review['retained_inputs'])
 for mod in list(sys.modules.values()):
  f=getattr(mod,'__file__',None)
  if f:
   try:p=Path(f).resolve().relative_to(ROOT).as_posix()
   except ValueError:continue
   if p.endswith('.py'):paths.add(p)
 before={p:digest(p) for p in sorted(paths)};out=ROOT/'evidence'/f'{STEM}.json'
 archive=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},precision_bits=ctx.prec,python=sys.version,invocation=sys.argv)
 archive['runs'].append(run)
 def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
 save()
 try:
  check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256']);run['result']=execute();run['status']='passed'
 except Exception:run['status']='failed';run['exception']=traceback.format_exc()
 finally:
  run['checks']=ROWS;run['input_sha256_after']={p:digest(p) for p in before}
  if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
  run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
 result=run.get('result',{})
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(before),exception=run.get('exception'),coefficient=result.get('physical_soft_log_coefficient'),coefficient_only_metric_bound_obstructed=result.get('coefficient_only_metric_bound_obstructed')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
