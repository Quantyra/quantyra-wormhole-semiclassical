"""R-squared candidate window; scoped vacuum/product calculation, not a wormhole."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as F
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='finite-r2-scalar-window'
PREVIOUS='finite-cold-local-gluing-and-admission'
OWN=['notes/'+STEM+'.md','checks/finite_r2_scalar_window.py','evidence/'+STEM+'-review.json']
ROWS=[]; FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
 assert name not in {x['name'] for x in ROWS}
 ROWS.append({'name':name,'passed':bool(ok),'detail':str(detail)})
 if not ok:raise AssertionError(name)
def algebra():
 t=S.symbols('t', real=True);s=S.Function('s')(t)
 d=lambda u,n=1:S.diff(u,t,n)
 L=36*(d(s,2)+d(s)**2)**2
 el=S.simplify(-d(S.diff(L,d(s)))+d(S.diff(L,d(s,2)),2))
 R=6*S.exp(-2*s)*(d(s,2)+d(s)**2)
 boxR=-S.exp(-2*s)*(d(R,2)+2*d(s)*d(R))
 check('nonlinear_conformal_R2_variation',S.simplify(el+12*S.exp(4*s)*boxR)==0)
 check('nonlinear_fourth_derivative_contact',S.simplify(el-72*(d(s,4)-6*d(s)**2*d(s,2)))==0)
 G,c,lam,z=S.symbols('G c lam z',positive=True)
 eps=G/(60*S.pi); Cv=1+S.Rational(11,3)*eps*lam;chi=96*S.pi*c-1/(60*S.pi)
 check('Lorentzian_trace_sign',S.simplify(eps-96*S.pi*G*c+G*chi)==0)
 check('anomaly_coefficient_shift',S.simplify(G/(2*S.pi)*(S.Rational(1,30)-192*S.pi**2*c)+G*chi)==0)
 check('threshold',S.simplify(chi.subs(c,1/(5760*S.pi**2)))==0)
 m2=Cv/(G*chi)
 check('original_negative_mass',S.simplify(m2.subs(c,0)+(60*S.pi/G+S.Rational(11,3)*lam))==0)
 psi=-6*(z-4*lam)
 trlift=6*(z-4*lam)*(Cv-G*chi*z)
 check('full_conformal_tensor_trace',S.simplify(trlift+(Cv-G*chi*z)*psi)==0)
 check('positive_mass_lift',S.simplify((Cv-G*chi*z).subs(z,m2))==0)
 x,y,p=S.symbols('x y p',positive=True)
 LR=4*c*(y-x)**2
 pa=S.expand(LR-x*S.diff(LR,x));ps=S.expand(LR-y*S.diff(LR,y))
 check('homogeneous_R2_pressures',S.simplify(pa-4*c*(y*y-x*x))==0 and S.simplify(pa+ps)==0)
 RR,phi=S.symbols('RR phi')
 aux=c*(2*phi*RR-phi*phi)
 check('auxiliary_equation',S.simplify(S.diff(aux,phi)-2*c*(RR-phi))==0)
 check('auxiliary_elimination',S.simplify(aux.subs(phi,RR)-c*RR**2)==0)
 slope=1/(16*S.pi*G)+2*c*phi
 check('local_boundary_coefficient',S.simplify(2*slope-(1/(8*S.pi*G)+4*c*phi))==0)
 # Direct logarithmic derivative of the actual global-AdS radial mode.
 D=S.symbols('D',positive=True);xx=S.symbols('xx',positive=True)
 logder=S.cot(xx)-(D-1)*S.tan(xx)
 res=-S.diff(logder,xx)-logder**2+(2+D*(D-3))/S.cos(xx)**2-D**2
 check('AdS_ground_radial_mode',S.trigsimp(S.expand_trig(res))==0)
 FIXTURES.append({'fixture':'exact_action_and_radial_identities','scalar_tower_truncated':False,'charged_response_evaluated':False})
def rational_window():
 plo=F('3.14159');phi=F('3.14160');k=2000
 g0=F(1,k**10);gl=(1-F('1e-34'))*g0
 # Verify the root enclosure by substitution for every pi in the bracket.
 check('outer_vacuum_root_lower',gl+11*gl*gl/(360*plo)<g0)
 check('outer_vacuum_root_upper',g0+11*g0*g0/(360*phi)>g0)
 Clo=F(49119000000);Chi=F(49120000000)
 check('C0_rational_interval',Clo*Clo<1+12*plo*k**6 and Chi*Chi>1+12*phi*k**6)
 Lvhi=48*phi*g0*Clo/(Clo-1)
 Lvlo=48*plo*g0*Chi/(Chi-1)
 def mass_bounds(c):
  chil=96*plo*c-1/(60*plo);chih=96*phi*c-1/(60*phi)
  cvl=1+11*gl/(180*phi);cvh=1+11*g0/(180*plo)
  return cvl/chih,cvh/chil,(1/g0+11/(180*phi))/chih,(1/gl+11/(180*plo))/chil
 low=mass_bounds(F('4e21'));high=mass_bounds(F('3e24'))
 check('whole_window_subplanck_target',low[1]<F('1e-24'))
 check('whole_window_curvature_separation',high[2]>F('1e6'))
 check('whole_window_product_existence',Lvhi*32*F('3e24')<F('1e-4'))
 c=F('1e22');mb=mass_bounds(c)
 check('selected_physical_mass_interval',F('3.31e-25')<mb[0] and mb[1]<F('3.32e-25'))
 check('selected_curvature_mass_interval',F('3.39e8')<mb[2] and mb[3]<F('3.40e8'))
 budget=Lvhi*32*c;Dref=F('7e-23')+budget;rad=F('9.43e-8')
 check('selected_product_budget',budget<F('4.713e-8'))
 check('selected_radius_certificate',2*Dref<rad)
 check('selected_product_F_positive',8*c*Lvhi/(1-rad)<F('1.20e-8'))
 check('selected_vacuum_F_positive',384*phi*c*g0<1)
 q=k**4;Mlo=F(q+1)*(1-rad)/(Chi*(1+rad));Mhi=F(q+1)*(1+rad)/(Clo*(1-rad))
 check('selected_lowest_mass_preserved',F('18.04')**2<Mlo and Mhi<F('18.06')**2)
 check('limit_point_scalar',mb[2]+2>F(3,4))
 # With s=Delta-1, the L2 radial exponents are s and 1-s.
 check('unique_KG_boundary_branch',mb[2]>F('1e6'))
 FIXTURES.append({'fixture':'rational_window','pi_bracket':[str(plo),str(phi)],'g_bracket':[str(gl),str(g0)],
 'C0_bracket':[str(Clo),str(Chi)],'Lambda_ref_bracket':[str(Lvlo),str(Lvhi)],
 'mass_intervals':[str(v) for v in mb],'product_budget_upper':str(budget),
 'radius_upper':str(2*Dref),'M1_squared_interval':[str(Mlo),str(Mhi)],
 'window':['4e21','3e24'],'selected_c':'1e22'})
def numeric():
 records=[]
 for digits in [60,90]:
  mp.mp.dps=digits;pi=mp.pi;k=mp.mpf(2000);g0=k**-10;a=11/(360*pi)
  g=2*g0/(1+mp.sqrt(1+4*a*g0));C0=mp.sqrt(1+12*pi*k**6);c=mp.mpf('1e22')
  Cv=1+11*g/(180*pi);chi=96*pi*c-1/(60*pi)
  mg=Cv/chi;ml=mg/g;Delta=mp.mpf('1.5')+mp.sqrt(mp.mpf('2.25')+ml)
  Lref=48*pi*g0*C0/(C0-1);budget=Lref*32*c
  mode_errors=[]
  for scale in ['.1','1','3']:
   x=mp.mpf(scale)/mp.sqrt(Delta)
   u=lambda xx:mp.sin(xx)*mp.cos(xx)**(Delta-1)
   val=u(x);res=-mp.diff(u,x,2)+(2+ml)*val/mp.cos(x)**2-Delta**2*val
   mode_errors.append(abs(res)/(Delta**2*abs(val)))
  check('actual_vacuum_mode_'+str(digits),max(mode_errors)<mp.mpf('1e-45'))
  # Compute product local pressures directly by varying the fixed-flux volume density.
  aa=mp.mpf('1.3');bb=mp.mpf('2.7')
  I=lambda ar,br:ar**2*br**2*c*(2/br**2-2/ar**2)**2
  pa=mp.diff(lambda ar:I(ar,bb),aa)/(2*aa*bb**2)
  expected=4*c*(bb**-4-aa**-4)
  err=abs((pa-expected)/expected)
  check('direct_product_pressure_'+str(digits),err<mp.mpf('1e-45'))
  records.append({'digits':digits,'g':mp.nstr(g,50),'c_star':mp.nstr(1/(5760*pi*pi),50),
   'mass_squared_G':mp.nstr(mg,50),'mass_squared_over_lambda':mp.nstr(ml,50),
   'mass_sqrtG':mp.nstr(mp.sqrt(mg),50),'Delta_plus':mp.nstr(Delta,50),
   'coefficient_min_mass':mp.nstr((mp.mpf('1e24')*Cv+1/(60*pi))/(96*pi),50),
   'coefficient_max_curvature':mp.nstr((Cv/(mp.mpf('1e6')*g)+1/(60*pi))/(96*pi),50),
   'coefficient_max_product':mp.nstr(mp.mpf('1e-4')/(32*Lref),50),
   'product_budget':mp.nstr(budget,50),'M1_classical':mp.nstr(mp.sqrt((k**4+1)/C0),50),
   'ground_mode_relative_errors':[mp.nstr(v,8) for v in mode_errors],'pressure_relative_error':mp.nstr(err,8)})
 for key in ['g','c_star','mass_squared_G','mass_squared_over_lambda','Delta_plus','product_budget']:
  check('precision_agreement_'+key,abs(mp.mpf(records[0][key])/mp.mpf(records[1][key])-1)<mp.mpf('1e-45'))
 FIXTURES.append({'fixture':'benchmark_and_actual_AdS_ground_mode','records':records,'scope':'exact vacuum comparison; no charged metric/state samples'})
def main():
 review=json.loads((ROOT/OWN[2]).read_bytes())
 old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
 paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})
 before={p:digest(p) for p in paths}
 run={'started_utc':datetime.now(timezone.utc).isoformat(),'input_sha256_before':before,
  'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
  'protocol':{'digits':[60,90],'numerical_relative_tolerance':'1e-45','continuum_machine_verified':False},
  'previous_goal_turn':'progress'}
 error=None
 try:
  check('preceding_milestone',old['status']=='passed' and old['general_cold_family_admitted'])
  check('all_previous_input_identities',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
  ref=json.loads((ROOT/'evidence/finite-local-reference-window.json').read_bytes())
  check('retained_full_reference_receipt',ref['status']=='PASS')
  check('reference_input_identities',all(digest(p)==h for p,h in ref['dependencies_sha256'].items()))
  check('reference_checker_identity',digest('checks/finite_local_reference_window.py')==ref['checker_sha256'])
  algebra();rational_window();numeric()
 except Exception:
  error=traceback.format_exc()
 after={p:digest(p) for p in paths}
 run.update(status='passed' if error is None and before==after else 'failed',exception=error,
  input_sha256_after=after,unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,
  finished_utc=datetime.now(timezone.utc).isoformat(),
  original_construction_preserved=True,selected_new_candidate_c='1e22',
  positive_vacuum_scalar=True if error is None else False,
  full_Dirac_homogeneous_theorem_applied=True if error is None else False,
  actual_charged_response_proved=False,actual_correction_norm_proved=False,
  nonlinear_solution_constructed=False,physical_boundary_completed=False,
  semiclassical_validity_proved=False,physical_stability_proved=False,
  goal_complete=False,component_complete=False,formal_verification=False)
 out=ROOT/'evidence'/(STEM+'.json')
 receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
 receipt['runs'].append(run)
 out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
 log=ROOT/'evidence'/(STEM+'.log')
 with log.open('a',encoding='utf-8',newline='\n') as stream:
  stream.write(json.dumps({'run':len(receipt['runs']),'started_utc':run['started_utc'],'input_sha256':before})+'\n')
  for row in ROWS:stream.write(json.dumps(row)+'\n')
  stream.write(json.dumps({'status':run['status'],'exception':error,'fixtures':FIXTURES})+'\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'inputs':len(paths),'fixtures':len(FIXTURES),'runs':len(receipt['runs']),'exception':error}))
 if run['status']!='passed':sys.exit(1)
if __name__=='__main__':main()
