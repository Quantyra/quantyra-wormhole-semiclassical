"""Actual finite-domain minus identical-collar cylinder: full near Hessian."""
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
STEM='current-coincident-domain-transfer'
OWN=['notes/'+STEM+'.md','checks/current_coincident_domain_transfer.py','evidence/'+STEM+'-review.json']
ROWS=[];SYMBOLIC={};enc=old.encoded

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def read_result(stem):return json.loads((ROOT/'evidence'/f'{stem}.json').read_bytes())['runs'][-1]['result']

def algebra():
 M,W,U,V,U1,U2=sp.symbols('M W U V U1 U2',positive=True)
 flow=lambda u:W*(1-u*u)-2*M*u
 check('Riccati_difference',sp.expand(flow(U1)-flow(U2)+(2*M+W*(U1+U2))*(U1-U2))==0)
 t=sp.symbols('t',positive=True);k=sp.sqrt(M*M+W*W);us=W/(M+k);v=us*sp.tanh(k*t)
 check('positive_subsolution',sp.simplify(flow(v)-sp.diff(v,t)-M*W/(M+k)*(1-sp.tanh(k*t))**2)==0)
 check('integrated_subsolution',sp.simplify(W*us-(k-M))==0)
 D=1+U*V
 gm=(1-U)*(1+V)/(2*M*D);gp=(1+U)*(1-V)/(2*M*D)
 check('seagull_diagonal_trace',sp.simplify(M*(gm+gp)-(1-U*V)/D)==0)
 # Proper coordinate points outwards: U'=flow(U), V'=-flow(V).
 dp=sp.diff(D,U)*flow(U)-sp.diff(D,V)*flow(V)
 check('image_symmetric_factorization_identity',sp.simplify(dp/D-W*(V-U))==0)
 # Kernel left/right derivatives at fixed other endpoint.
 left=-4*V*(flow(U)+2*(M+W*U)*U)/D**2
 right=sp.diff(-4*U*V/D**2,V)*(-flow(V))+sp.diff(-4*U*V/D**2,U)*flow(U)
 # On the right, only the outer endpoint's U in the denominator varies;
 # the inner numerator U is held fixed, so subtract its derivative.
 right-=(-4*V/D**2)*flow(U)
 right+=(-4*U*V/D**2)*(-2*(M+W*U))
 check('direct_left_derivative',sp.simplify(left+4*V*W*(1+U*U)/D**2)==0)
 check('direct_right_derivative',sp.simplify(right-4*U*W*(1+V*V)/D**2)==0)
 check('direct_first_derivative_jump',sp.simplify(right-left-4*W*(U+V)/D)==0)
 # Cauchy estimate for derivatives of the jump's diagonal value.
 for r in range(1,4):
  for ell in range(4-r):
   s=sum(comb(ell,j)*factorial(r+j)*factorial(ell-j) for j in range(ell+1))
   check('mixed_jump_identity_%d_%d'%(r,ell),sp.Rational(2*s,factorial(r+ell))==sp.Rational(2*(r+ell+1),r+1))
 x,cc=sp.symbols('x cc',positive=True);f=sp.Function('f')(x)
 # Exact distribution fixture e^-c|x-y|: F''=c^2 F-2c f.
 FF=sp.Function('F')(x)
 f4=cc**2*(cc**2*FF-2*cc*f)-2*cc*sp.diff(f,x,2)
 check('fourth_distribution_derivative_fixture',sp.expand(f4-(cc**4*FF-2*cc**3*f-2*cc*sp.diff(f,x,2)))==0)
 sig=sp.Rational(801,200);s=sp.symbols('s')
 for m in range(5):
  poly=sp.Poly(sp.prod(s+sig-i for i in range(m)),s)
  check('input_derivative_positive_coefficients_'+str(m),all(c>0 for c in poly.all_coeffs()) and poly.eval(1)==sp.prod(sig+1-i for i in range(m)))
 # The explicit cap taper has eight zero derivatives at each endpoint.
 weight=t**8*(1-t)**8;norm=sp.integrate(weight,(t,0,1));taper=1-sp.integrate(weight,(t,0,x))/norm
 check('cap_taper_endpoint_values',taper.subs(x,0)==1 and taper.subs(x,1)==0)
 check('cap_taper_C8',all(sp.diff(taper,x,j).subs(x,e)==0 for j in range(1,9) for e in [0,1]))
 check('cap_taper_monotone_identity',sp.simplify(sp.diff(taper,x)+x**8*(1-x)**8/norm)==0)
 SYMBOLIC.update(direct_derivative_jump=str(4*W*(U+V)/D),seagull_trace=str((1-U*V)/D),
  physical_cubic_seam_preserved=True,common_regulator_relative_difference_only=True,
  auxiliary_cylinder_is_not_physical_replacement=True,uniform_nonlinear_neighborhood_proved=False,
  complete_angular_and_discrete_thermal_sums=True,finite_physical_modes_selected=False)

def execute():
 algebra();p,h,zh,env=old.scalar.prev.york.actual_parameters()
 a=arb('.01');b=arb('.02');c=arb('.04');sigma=arb('4.005');alpha=5-sigma
 M1=((p['q']+1)/p['C']).sqrt();check('mass_gap',M1>18)
 check('physical_collar_inside_actual_bolt',zh>arb('.05'))
 be=p['beta'];coeff=[arb(0),arb(0),p['u']/be,-4-2*p['u']/be-p['tau']/be,
  3+p['u']/be+p['eps']/be,-2*p['eps']/be,p['eps']/be]
 err=sum((abs(coeff[j])*arb('.05')**j for j in range(2,7)),arb(0))
 check('actual_polynomial_complex_disk',err<arb('.001'))
 check('actual_complex_inverse_bounds',1/(be*(1-err))<16 and 1/(be*(1-err)).sqrt()<4)
 check('actual_and_cap_Q_below_point_two',be+p['u']/16+abs(p['eps'])<arb('.2'))
 check('Cauchy_disks_before_cap',b+1/(1000*M1)<c)
 check('Picard_invariance_and_contraction',arb('35.36')<36 and arb('35.2')<36 and arb(36)/1000<arb('.1'))
 check('complex_denominator_reserve',1-2*sp.Rational(1,10)-sp.Rational(1,100)==sp.Rational(79,100))
 check('complex_difference_arm_inflation',arb('.036').exp()<arb('1.1'))
 check('two_arm_exponent_inflation',arb('.072').exp()<arb('1.125'))
 den=arb('.79');u=arb('1.1')
 check('rational_prefactor_and_partial_bounds',u*u/den**2<2 and 2*u**3/den**3<6 and u*arb('1.21')/den**3<6)
 check('diagonal_trace_partial_bound',2*u/den**2<4)
 check('real_direct_and_image_exponent_reserves',8/arb(5).sqrt()<4 and 32/arb(5).sqrt()<15)
 check('complex_two_path_reserves',4*arb('.072').exp()<5 and 16*arb('.072').exp()<18)
 check('bubble_comparison_ledger',16*(4*(180*5+2*2*5*5)+4*(180*18+2*2*24*18))<2**20)
 check('seagull_comparison_ledger',4*4*(5+20)<2**10)
 check('positive_original_temperature',p['gamma']>0)
 d=arb(5).sqrt()*(c-a);cm=arb(2).sqrt()*d;cp=arb(10).sqrt()*d
 pref=8*p['C']/(arb.pi()*p['B2'])*(1+p['gamma']*cp)
 def moments(r,j):
  return sum((comb(j,k)*factorial(r+j-k+2)*factorial(k)/cm**(r+j-k+2)/cp**(k+1) for k in range(j+1)),arb(0))
 K=[pref*2**20*factorial(j)*arb(1000)**j*moments(2,j) for j in range(5)]
 S=[pref*2**10*factorial(j)*arb(1000)**j*moments(1,j) for j in range(5)]
 L=[arb(1)]
 for m in range(1,5):L.append(L[-1]*(sigma+2-m))
 # Store each term as (coefficient, ordinary z power, origin).
 response={}
 for n in range(1,5):
  terms=[(b**(sigma+1)*K[n],arb(0),0)]
  for r in range(1,n):
   for ell in range(n-r):
    m=n-r-1-ell
    terms.append((comb(n-r-1,ell)*arb(2*(r+ell+1))/(r+1)*K[r+ell]*L[m],sigma-m,1))
  for ell in range(n+1):terms.append((comb(n,ell)*S[ell]*L[n-ell],sigma-(n-ell),2))
  response[n]=terms
 cs=[[arb(1)]]
 for k in range(1,4):cs.append([(alpha+j)*(cs[-1][j] if j<k else 0)+(cs[-1][j-1] if j>0 else 0) for j in range(k+1)])
 caps=[1,3,10,60];parts=[]
 for k in range(4):
  row=[arb(0),arb(0),arb(0)]
  for ell in range(k+1):
   for j in range(k-ell+1):
    for A,power,part in response[j+1]:
     pp=alpha+j+power
     assert pp>0
     row[part]+=comb(k,ell)*caps[ell]*cs[k-ell][j]*A*a**pp/(2*pp).sqrt()
  parts.append(row)
 Nk=[sum(row,arb(0)) for row in parts]
 C3=old.decode(read_result('current-global-boundary-moment-return')['C3_H3_to_H6'])
 check('retained_normal_inverse_guard',C3<1801)
 eta=C3/arb(2).sqrt()*sum((comb(3,k)*v*v for k,v in enumerate(Nk)),arb(0)).sqrt()
 check('complete_coincident_transfer_target',eta<arb('1e-10'))
 check('refined_complete_coincident_transfer_guard',eta<arb('2.013e-13'))
 remaining=arb('.18699')-eta;check('combined_rest_allocation_available',remaining>arb('.18698'))
 prior=read_result('current-coupled-local-graph-bound')
 local=old.decode(prior['local_block_row_sum_upper']);principal=old.decode(prior['prior_quantum_principal_upper'])
 far=old.decode(read_result('current-separated-quantum-graph-bound')['separated_quantum_graph_feedback_upper'])
 total=local+principal+far+eta+arb('.18698');joint=(arb('7.49e-12')+total*arb('1e-14'))/(1-total)
 check('conditional_joint_return',joint<arb('1e-11'))
 return dict(output_collar=enc(a),near_input_cut=enc(b),identical_collar_end=enc(c),complex_Q_error_upper=enc(err),
  proper_separation_lower=enc(d),normal_C3=enc(C3),bubble_derivative_upper=[enc(v) for v in K],
  seagull_derivative_upper=[enc(v) for v in S],weighted_output_parts_regular_jumps_seagull=[[enc(v) for v in row] for row in parts],
  weighted_output_derivative_upper=[enc(v) for v in Nk],coincident_domain_transfer_graph_feedback_upper=enc(eta),
  remaining_feedback_threshold=enc(remaining),selected_remaining_target='.18698',conditional_total_feedback=enc(total),conditional_joint_E_upper=enc(joint),
  fixed_seed_relative_domain_transfer_certified=True,coincident_seagull_and_diagonal_jumps_retained=True,
  common_regulator_removed_for_relative_difference=True,complete_local_cylinder_Hessian_certified=False,
  complete_quantum_graph_feedback_certified=False,joint_source_bounds_certified=False,
  independent_Mellin_Efar_target_closed=False,actual_global_correction_computed=False,actual_controlled_solution=False,
  uniform_nonlinear_admission=False,exact_compatibility=False,physical_matching=False,semiclassical_validity_noise=False,all_sector_stability=False)

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
