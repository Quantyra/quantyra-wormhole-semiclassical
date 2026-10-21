"""Actual seed volume-norm obstruction: ordinary proof and independent fixtures."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='outer-volume-norm-singular-sequence'
OWN=['checks/outer_volume_norm_singular_sequence.py','notes/'+STEM+'.md','evidence/volume-singular-source-review.json']
PREVIOUS='annular-retarded-operator-limit'
ROWS=[];NUMERICS=[]
PROTOCOL={'scope':'No closed-range finite-kernel realization in equal spacetime-volume weights containing compact conformal tests',
 'mpmath_precision':55,'identity_tolerance':'1e-35','carrier_frequencies':[3,7,20],
 'window_lengths':[8,16,32,64],'residual_refinement_max':'0.55',
 'envelope':'(1-x^2)^6 on [-1,1], normalized in L2; H5 and approximable in H4 by smooth compact tests',
 'curved_fixture':{'Phi':'-log Y+Y^2/10','V':'-Y^2/20','q':2,'Q2':3,'k':3,
 'Y':['0.05','0.025','0.0125'],'decay_refinement_max':'0.27'},
 'fixtures_are_actual_exterior_spectrum':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
 print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,expr):
 value=S.simplify(S.sympify(expr).doit());check(name,value==0,value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(OWN)|set(prior()['input_sha256_after'])|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior()
 check('predecessor_passed',p['status']=='passed' and len(p['checks'])==90 and all(r['passed'] for r in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('all_predecessor_hashes',not bad,bad)
 required=['local-coupled-frequency-polynomial','charged-conformal-response-and-constraints','outer-coupled-scalar-separation','fast-scalar-boundary-solvability']
 check('required_local_columns_pinned',all('notes/'+f+'.md' in p['input_sha256_after'] for f in required))
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('primary_source_scope',len(review['primary_sources'])==1 and not review['physical_metric_boundary_selected']
  and not review['actual_charged_instability_proved'] and not review['order_reduction_adopted'])

def algebra():
 Y=S.symbols('Y',positive=True);r,z,k,G=S.symbols('r z k G')
 e=60*S.pi/G+S.Rational(11,3)
 D=lambda f:Y**2*S.diff(f,Y,2)-2*Y*S.diff(f,Y)
 Bstar=lambda f:Y**2*S.diff(f,Y,2)-4*Y*S.diff(f,Y)-6*f
 for name,op,target in [
  ('D',D,z*z-S.Rational(9,4)),
  ('A',lambda f:D(f)-4*f,z*z-S.Rational(25,4)),
  ('B',lambda f:Y**2*S.diff(f,Y,2)-12*f,z*z+2*z-S.Rational(45,4)),
  ('Bstar',Bstar,z*z-2*z-S.Rational(45,4))]:
  zero('log_radius_conjugation_'+name,op(Y**r)/Y**r-target.subs(z,r-S.Rational(3,2)))
 zero('weighted_B_transpose',Y**4*S.diff(Y**-2*Y**r,Y,2)-12*Y**r-Bstar(Y**r))
 zero('fast_wave_number',e-S.Rational(9,4)-(60*S.pi/G+S.Rational(17,12)))
 P=[(z*z-S.Rational(25,4))*(z*z+k*k)/(20*S.pi),
    (z*z-2*z-S.Rational(45,4))*(z*z+k*k)/(60*S.pi)]
 for i,p in enumerate(P):zero('both_rows_annihilate_carrier_'+str(i),p.subs(z,S.I*k))
 zero('volume_density_unitarity',Y**-4*Y**3*Y-1)
 zero('wave_density_after_conjugation',Y**-2*Y**3*Y-Y**2)
 w=S.Function('w')(Y)
 zero('physical_wave_rescaling',-(D(Y*w)+e*Y*w)/Y**3-(-S.diff(w,Y,2)-(e-2)*w/Y**2))
 L,R,x=S.symbols('L R x',positive=True)
 zero('support_end_exponent',2*(-R+L)+2*(R-L))
 n=S.symbols('n',positive=True,integer=True)
 # Successive log supports are eventually disjoint for R_n=n^3+10,L_n=n.
 check('log_supports_eventually_disjoint',S.expand(((n+1)**3-(n+1))-(n**3+n))==3*n*n+n)
 check('sequence_error_tends_to_zero',S.limit(1/n+S.exp(-2*(n**3+10-n)),n,S.oo)==0)
 check('finite_projection_of_escaping_ONB',S.limit(1/n,n,S.oo)==0)
 check('lapse_multiplier_unbounded',S.limit(Y**-2,Y,0,dir='+')==S.oo)
 check('phase_scale_above_180',60*S.pi>180)
 # A small unweighted residual does not bound its lapse-multiplied value.
 zero('lapse_weight_counter_scaling',Y*Y**-2-Y**-1)
 return P,z,k

def action_columns():
 """Differentiate the exact local density, independently of normal factorization."""
 P=S.symbols('P0:5');V=S.symbols('V0:5')
 Y,r,G,la,Q,q=S.symbols('Y r G la Q q',positive=True)
 def total(f):return S.expand(sum(S.diff(f,P[i])*P[i+1]+S.diff(f,V[i])*V[i+1] for i in range(4)))
 FC=S.exp(2*P[0])*(6*P[1]**2+4*P[1]*V[1])+S.exp(2*V[0])*(2*S.exp(2*P[0])+6*la*S.exp(4*P[0])-2*Q)
 FW=S.exp(-2*V[0])*(P[0]*V[2]**2/15-P[2]*V[2]/15-P[2]**2/10-V[2]*P[1]**2/15
   -7*P[2]*P[1]**2/90-7*P[1]**4/180)+P[0]*V[2]/9-P[1]**2/18+(S.Rational(1,15)-q*q/3)*P[0]*S.exp(2*V[0])
 density=-FC/(4*G)-FW/(4*S.pi)
 bg={P[i]:S.diff(-S.log(Y),Y,i) for i in range(5)}|{v:0 for v in V}|{la:1+11*G/(360*S.pi)}
 e=60*S.pi/G+S.Rational(11,3);fast=r*(r-3)+e
 targets=[fast*(r*(r-3)-4)/(20*S.pi),(r*r-5*r-6)*fast/(60*S.pi)]
 columns=[]
 for i,J in enumerate([P,V]):
  E=S.diff(density,J[0])-total(S.diff(density,J[1]))+total(total(S.diff(density,J[2])))
  coeff=[S.diff(E,P[j]) for j in range(5)]
  col=sum(coeff[j]*S.diff(Y**r,Y,j) for j in range(5))
  reference=S.simplify(col.subs(bg)*Y**(4-r))
  expected=targets[i]+Y**2*(-r*(r-1)/(36*S.pi)-2/G)
  if i==1:expected+=Y**4*(5*q*q-1)/(30*S.pi)
  zero('full_local_action_column_'+str(i),reference-expected)
  zero('normal_indicial_column_'+str(i),S.limit(reference,Y,0)-targets[i])
  # Unmatched cosmological coefficient is an actual order-one error, not a tail.
  delta=S.symbols('delta')
  wrong=S.simplify(col.subs(bg|{la:1+11*G/(360*S.pi)+delta})*Y**(4-r)-reference)
  zero('wrong_lambda_contact_detected_'+str(i),wrong+([24,12][i])*delta/G)
  columns.append(col)
 # Curved background: exact retained local action including charge and sphere.
 kval=3;gval=60*S.pi/(S.Rational(kval*kval)-S.Rational(17,12))
 curved={P[i]:S.diff(-S.log(Y)+Y**2/10,Y,i) for i in range(5)}
 curved|={V[i]:S.diff(-Y**2/20,Y,i) for i in range(5)}
 curved|={la:1+11*gval/(360*S.pi),G:gval,Q:3,q:2,r:S.Rational(3,2)+S.I*kval}
 functions=[]
 for i,col in enumerate(columns):
  expr=S.simplify((col*Y**(4-r)).subs(curved))
  zero('curved_column_constant_vanishes_'+str(i),S.limit(expr,Y,0))
  leading=S.simplify(S.limit(expr/Y**2,Y,0))
  check('curved_column_has_nonzero_second_order_remainder_'+str(i),leading!=0)
  functions.append(S.lambdify(Y,expr,'mpmath'))
 old=None
 for yy in PROTOCOL['curved_fixture']['Y']:
  y=mp.mpf(yy);vals=[f(y) for f in functions];norm=mp.sqrt(sum(abs(v)**2 for v in vals))
  if old is not None:check('curved_coefficient_refinement_'+yy,norm<mp.mpf('.27')*old,mp.nstr(norm/old,15))
  check('curved_coefficient_finite_'+yy,mp.isfinite(norm) and norm>0)
  NUMERICS.append({'fixture':'full_local_curved_column','Y':yy,'column_norm':mp.nstr(norm,22),'divided_by_Y_squared':mp.nstr(norm/y**2,22)})
  old=norm
 # Fixed frequency terms in the pure AdS normal column decay as Y^2,Y^4.
 sig=S.symbols('sig')
 D=lambda f:Y**2*(S.diff(f,Y,2)-sig*f)-2*Y*S.diff(f,Y)
 Bs=lambda f:Y**2*(S.diff(f,Y,2)-sig*f)-4*Y*S.diff(f,Y)-6*f
 pp=S.simplify((D(D(Y**r)-4*Y**r)+e*(D(Y**r)-4*Y**r))/(20*S.pi*Y**r))
 vp=S.simplify(Bs(D(Y**r)+e*Y**r)/(60*S.pi*Y**r))
 for i,col in enumerate([pp,vp]):
  zero('fixed_frequency_limit_'+str(i),S.limit(col,Y,0)-targets[i])
  zero('quartic_frequency_weight_'+str(i),S.expand(col).coeff(sig,2)-Y**4/([20,60][i]*S.pi))
 # Carrier does NOT solve the finite-frequency spatial equation exactly.
 check('carrier_is_not_an_exact_mode',S.simplify((D(Y**r)+e*Y**r).subs(r,S.Rational(3,2)+S.I*3).subs(G,gval))!=0)

def envelope_fixtures(P,z,k):
 x=S.symbols('x',real=True);chi=(1-x*x)**6
 mass=S.integrate(chi*chi,(x,-1,1));chi=chi/S.sqrt(mass)
 derivatives=[S.diff(chi,x,j) for j in range(5)]
 for j in range(6):zero('compact_envelope_endpoint_'+str(j),S.diff(chi,x,j).subs(x,1))
 zero('envelope_unit_norm',S.integrate(chi*chi,(x,-1,1))-1)
 norms=[mp.sqrt(mp.mpf(str(S.N(S.integrate(f*f,(x,-1,1)),60)))) for f in derivatives]
 for kval in PROTOCOL['carrier_frequencies']:
  old=None
  for length in PROTOCOL['window_lengths']:
   polys=[];bounds=[]
   for pp in P:
    pp=pp.subs(k,kval)
    envelope=S.expand(sum(S.diff(pp,z,j).subs(z,S.I*kval)*derivatives[j]/(S.factorial(j)*length**j) for j in range(1,5)))
    # Independently apply the differential polynomial to the modulated envelope.
    modulated=chi*S.exp(S.I*kval*length*x)
    direct=sum(c*S.diff(modulated,x,j)/length**j for (j,),c in S.Poly(pp,z).terms())
    zero('direct_column_envelope_'+str(kval)+'_'+str(length)+'_'+str(len(polys)),S.expand(direct/S.exp(S.I*kval*length*x))-envelope)
    polys.append(envelope)
    bounds.append(sum(abs(complex(S.N(S.diff(pp,z,j).subs(z,S.I*kval),18)))*norms[j]/(mp.factorial(j)*length**j) for j in range(1,5)))
   # Exact polynomial moment integration versus independent numerical quadrature.
   square=S.expand(sum(f*S.conjugate(f) for f in polys))
   exact=S.integrate(square,(x,-1,1))
   exact_norm=mp.sqrt(mp.mpf(str(S.N(exact,60))))
   funcs=[S.lambdify(x,f,'mpmath') for f in polys]
   quad_norm=mp.sqrt(mp.quad(lambda y:sum(abs(f(y))**2 for f in funcs),[-1,0,1]))
   error=abs(quad_norm-exact_norm)/(1+exact_norm)
   check('independent_residual_quadrature_'+str(kval)+'_'+str(length),error<mp.mpf('1e-35'),mp.nstr(error,12))
   check('explicit_residual_bound_'+str(kval)+'_'+str(length),exact_norm<=sum(bounds),mp.nstr(exact_norm/sum(bounds),15))
   if old is not None:check('residual_refinement_'+str(kval)+'_'+str(length),exact_norm<mp.mpf('.55')*old,mp.nstr(exact_norm/old,15))
   NUMERICS.append({'fixture':'compact_conformal_two_row_residual','k':kval,'L':length,
    'norm':mp.nstr(exact_norm,22),'quadrature_relative_error':mp.nstr(error,12),'L_times_norm':mp.nstr(length*exact_norm,22)})
   old=exact_norm

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 mp.mp.dps=PROTOCOL['mpmath_precision'];ps=paths()
 hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
 before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 volume-norm singular sequence','input_count':len(ps)}),flush=True)
 try:provenance();P,z,k=algebra();action_columns();envelope_fixtures(P,z,k)
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
 run={'story':'S168 volume-norm singular sequence','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'checks':ROWS,'exception':error,'numerics':NUMERICS,
  'actual_seed_volume_singular_sequence_derived':good,'both_metric_rows_controlled':good,
  'volume_norm_closed_range_obstruction_derived':good,'physical_lapse_weight_distinguished':good,
  'all_local_contacts_retained':good,'physical_metric_boundary_selected':False,
  'actual_metric_causal_inverse_proved':False,'actual_metric_mode_constructed':False,
  'quantitative_physical_collar_radius_proved':False,'actual_charged_stability_proved':False,
  'actual_wormhole_instability_proved':False,'retained_fast_factor_removed':False,
  'order_reduction_adopted':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text(encoding='utf-8')
  run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
