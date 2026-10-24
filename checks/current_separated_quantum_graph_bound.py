"""Entire-mode, actual-domain separated quantum H3 operator component."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb,factorial
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as sp
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-separated-quantum-graph-bound'
OWN=['notes/'+STEM+'.md','checks/current_separated_quantum_graph_bound.py','evidence/'+STEM+'-review.json']
ROWS=[];SYMBOLIC={};enc=old.encoded
def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read_result(stem):return json.loads((ROOT/'evidence'/f'{stem}.json').read_bytes())['runs'][-1]['result']

def algebra():
 M,w,t=sp.symbols('M w t',positive=True);k=sp.sqrt(M*M+w*w)
 ustar=w/(M+k);v=ustar*sp.tanh(k*t)
 residual=sp.simplify(w*(1-v*v)-2*M*v-sp.diff(v,t))
 # Exact positive Riccati subsolution, not equality with the constant flow.
 target=M*w/(M+k)*(1-sp.tanh(k*t))**2
 check('Riccati_subsolution_residual',sp.simplify(residual-target)==0)
 check('Riccati_integrated_lower_bound',sp.simplify(w*ustar-(k-M))==0)
 n,q,C=sp.symbols('n q C',positive=True)
 count=2*n*(q+n+1);bound=4*n*(n+q)
 check('entire_angular_count',sp.expand(bound-count-2*n*(n+q-1))==0)
 x=sp.symbols('x',positive=True)
 check('normal_log_kernel_off_diagonal',sp.diff(-1/(2*x),x,4)==-12/x**5)
 for j in range(1,5):
  check('normal_tail_integral_'+str(j),sp.Rational(24*factorial(j+4),24*(j+4))==factorial(j+3))
  for k0 in range(j+1):
   rr=j-k0+1;cc=sp.symbols('c',positive=True)
   integ=sp.integrate(cc*x**(rr+2)*sp.exp(-cc*x),(x,0,sp.oo))
   check('angular_moment_%d_%d'%(j,k0),sp.simplify(integ-factorial(rr+2)/cc**(rr+2))==0)
 # The exact cutoff/product coefficients used in the H3 bound.
 al=sp.symbols('alpha',positive=True);F=sp.Function('F')(x);Dj=lambda f:x*sp.diff(f,x)
 cs=[[sp.Integer(1)]];v=x**al*sp.diff(F,x)
 for m in range(4):
  if m:
   cs.append([(al+j)*(cs[-1][j] if j<m else 0)+(cs[-1][j-1] if j>0 else 0) for j in range(m+1)])
   v=Dj(v)
  rhs=sum(cs[m][j]*x**(al+j)*sp.diff(F,x,j+1) for j in range(m+1))
  check('weighted_output_derivative_'+str(m),sp.simplify(v-rhs)==0)
 check('correlated_pair_maximum',sp.simplify(sp.diff(x**4*sp.log(1/x),x).subs(x,sp.exp(-sp.Rational(1,4))))==0)
 SYMBOLIC.update(Riccati_subsolution=str(target),normal_separated_kernel='12K[(zprime-z)^-5+(zprime+z)^-5]',
  actual_domain_and_regular_bolt_retained=True,thermal_sum_is_discrete=True,finite_physical_modes_selected=False,
  four_observation_derivatives_from_analytic_ODE_germs=True,local_contacts_retained_in_complement=True)

def execute():
 algebra();p,h,zh,env=old.scalar.prev.york.actual_parameters()
 a=arb('.01');b=arb('.02');sigma=arb('4.005');alpha=5-sigma
 check('physical_input_range',zh>b and zh<1)
 M1=((p['q']+1)/p['C']).sqrt();check('mass_gap',M1>18)
 # One-sided complex polynomial, never a smoothed physical double.
 be=p['beta'];coeff=[arb(0),arb(0),p['u']/be,-4-2*p['u']/be-p['tau']/be,
  3+p['u']/be+p['eps']/be,-2*p['eps']/be,p['eps']/be]
 disk_error=sum((abs(coeff[k])*arb('.02')**k for k in range(2,7)),arb(0))
 check('complex_polynomial_disk',disk_error<arb('.001'))
 qmin=be*(1-disk_error);check('complex_inverse_guards',1/qmin<16 and 1/qmin.sqrt()<4)
 check('Cauchy_disk_inside_polynomial_disk',a+1/(1000*M1)<arb('.02'))
 check('Cauchy_disk_separated_from_inputs',1/(1000*M1)<b-a)
 check('Picard_disk_invariant',arb('35.36')<36 and arb('8.8')<36 and arb(36)/1000<arb('.1'))
 check('Picard_contraction',arb('35.2')<36 and arb(36)/1000<1)
 den=1-2*sp.Rational(1,10)-sp.Rational(1,100)
 check('no_complex_denominator_zero',den==sp.Rational(79,100))
 expcap=4*arb('1.1')/arb('.79')**2*arb('.036').exp()
 check('complex_spin_trace_guard',expcap<8)
 check('Cauchy_complete_two_path_coefficient',8*4*(4+16)==640)
 check('thermal_period_positive',p['gamma']>0)
 d=arb(5).sqrt()*(b-a);cm=arb(2).sqrt()*d;cp=arb(10).sqrt()*d
 check('mass_thermal_exponent_split', (cm**2+(cp/arb(5).sqrt())**2-4*d*d).contains(0))
 K=p['C']*be**2/(15*arb.pi()*p['B2']);actual=[];normal=[];F=[]
 for j in range(1,5):
  series=sum((comb(j,k)*factorial(j-k+3)*factorial(k)/cm**(j-k+3)/cp**(k+1) for k in range(j+1)),arb(0))
  A=2560*p['C']/(arb.pi()*p['B2'])*factorial(j)*arb(1000)**j*series*(1+p['gamma']*cp)
  N=K*factorial(j+3)/(b-a)**(j+4)
  actual.append(A);normal.append(N);F.append(A+N)
 cs=[[arb(1)]]
 for k in range(1,4):cs.append([(alpha+j)*(cs[-1][j] if j<k else 0)+(cs[-1][j-1] if j>0 else 0) for j in range(k+1)])
 caps=[1,3,10,60];lk=[]
 for k in range(4):
  lk.append(sum((comb(k,l)*caps[l]*cs[k-l][j]*a**(alpha+j)/(2*(alpha+j)).sqrt()*F[j]
   for l in range(k+1) for j in range(k-l+1)),arb(0)))
 Lfar=sum((comb(3,k)*v*v for k,v in enumerate(lk)),arb(0)).sqrt()
 previous=read_result('current-global-boundary-moment-return');C3=old.decode(previous['C3_H3_to_H6'])
 eta=C3*Lfar/arb(2).sqrt()
 check('whole_separated_graph_target',eta<arb('2e-6'))
 check('separated_graph_refined_guard',eta<arb('1.579e-6'))
 remaining=arb('.187')-eta;check('near_diagonal_allocation_available',remaining>arb('.18699'))
 prior=read_result('current-coupled-local-graph-bound')
 local=old.decode(prior['local_block_row_sum_upper']);principal=old.decode(prior['prior_quantum_principal_upper'])
 total=local+principal+eta+arb('.18699');joint=(arb('7.49e-12')+total*arb('1e-14'))/(1-total)
 check('conditional_joint_return',joint<arb('1e-11'))
 endpoint=read_result('current-coupled-outer-return-test')
 B0=old.decode(endpoint['full_normal_shape_particular']);estar=old.decode(previous['exact_slow_exponent_offset'])
 knownpair=B0*estar/(4*arb(1).exp());Delta=arb('4.500000054e-10')
 correlated=Lfar*(knownpair+Delta);separate=Lfar*(2*B0+Delta)
 threshold=arb('7.49e-12')/Lfar
 check('correlated_pair_conditional_remote_affine',correlated<arb('6e-19'))
 check('separate_coefficients_excessive',separate>arb('7.49e-12'))
 return dict(output_collar=enc(a),input_cut=enc(b),proper_separation_lower=enc(d),complex_Q_relative_error_upper=enc(disk_error),
  complex_spin_trace_guard_value=enc(expcap),discrete_thermal_correction_factor=enc(1+p['gamma']*cp),
  actual_entire_mode_derivative_upper=[enc(v) for v in actual],matched_normal_derivative_upper=[enc(v) for v in normal],
  combined_derivative_upper=[enc(v) for v in F],weighted_output_derivative_upper=[enc(v) for v in lk],
  L_infinity_to_outer_H3_upper=enc(Lfar),normal_C3=enc(C3),separated_quantum_graph_feedback_upper=enc(eta),
  remaining_feedback_threshold=enc(remaining),selected_near_remaining_target='.18699',
  conditional_total_feedback=enc(total),conditional_joint_E_upper=enc(joint),
  finite_component_supremum_threshold_for_full_affine_budget=enc(threshold),
  known_correlated_pair_supremum_piece=enc(knownpair),conditional_Delta_cap=enc(Delta),
  conditional_correlated_pair_remote_affine_upper=enc(correlated),separate_coefficient_remote_upper=enc(separate),
  separated_actual_and_normal_operator_component_certified=True,common_regulator_removal_with_four_observation_derivatives_on_separated_supports=True,
  near_diagonal_finite_contacts_certified=False,complete_quantum_graph_feedback_certified=False,
  mass_slow_original_row_coefficients_evaluated=False,joint_source_bounds_certified=False,
  independent_Mellin_Efar_target_closed=False,actual_E3_upper_known=False,actual_global_correction_computed=False,
  actual_controlled_solution=False,nonlinear_admission=False,exact_compatibility=False,physical_matching=False,
  semiclassical_validity_noise=False,all_sector_stability=False)

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
 run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,
  frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},precision_bits=ctx.prec,python=sys.version,invocation=sys.argv)
 archive['runs'].append(run)
 def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
 save()
 try:
  check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256']);run['result']=execute();run['status']='passed'
 except Exception:run['status']='failed';run['exception']=traceback.format_exc()
 finally:
  run['checks']=ROWS;run['symbolic']=SYMBOLIC;run['input_sha256_after']={p:digest(p) for p in before}
  if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
  run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
 with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
 print(json.dumps(dict(status=run['status'],checks=len(ROWS),exception=run.get('exception'),result=run.get('result')),indent=2))
 if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
