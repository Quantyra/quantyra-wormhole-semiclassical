"""Charged reconstruction trace identity and explicit model boundary tests."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import importlib.util,json,subprocess,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='coupled-fast-metric-boundary-traces'
OWN=['checks/coupled_fast_metric_boundary_traces.py','notes/'+STEM+'.md','evidence/fast-metric-traces-source-review.json']
PREVIOUS='coupled-thermal-boundary-index';ROWS=[];FIXTURES=[];R=S.Rational
PROTOCOL={'scope':'Fast metric traces and auxiliary boundary-domain indices; physical selection unresolved',
 'precision_digits':[50,80],'model_lambda':'2','model_cases':[['1.5','.7'],['3','1.2'],['7','2.8']],
 'model_columns':['k','gamma'],'observation_points':['.02','.01','.005'],
 'shape_exponent':'1.25','shape_amplitude_in_units_of_trace_coefficient':'0.05',
 'green_identity_relative_tolerance':'1e-35','precision_relative_tolerance':'1e-35',
 'fine_curvature_trace_error_limit':'.001','fine_full_metric_trace_error_limit':'.01',
 'full_metric_error_refinement_ratio':'.8','model_is_charged_seed':False,
 'actual_shape_kernel_evaluated':False,'formal_verification':False}

def load(name,path):
 spec=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
old=load('pinned_reconstruction','checks/physical_curvature_reconstruction.py')

def check(name,ok,detail=''):
 assert name not in {c['name'] for c in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)
def zero(name,value):
 value=S.simplify(value);check(name,value==0,value)
def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())

def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'}|set(review()['additional_committed_inputs']))

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==52 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 rev=review();extra=rev['additional_committed_inputs']
 check('additional_inputs_exactly_identified',set(extra)==set(rev['retained_inputs'])-set(p['input_sha256_after']))
 bad=[]
 for f,item in extra.items():
  committed=subprocess.check_output(['git','show',item['commit']+':'+f],cwd=ROOT)
  if committed!=(ROOT/f).read_bytes() or sha256(committed).hexdigest()!=item['sha256']:bad.append(f)
 check('additional_inputs_match_committed_sources',not bad,bad)
 check('review_inputs_pinned',all(f in paths() for f in rev['retained_inputs']))
 check('physical_scope',not rev['physical_boundary_selected'] and not rev['physical_mode_constructed'] and not rev['full_boundary_action_evaluated'])

def exact_checks():
 r,lam,m,Q,c3,c4,E,gamma=S.symbols('r lam m Q c3 c4 E gamma',positive=True)
 F=lam*r*r+1-2*m/r+Q/r**2+c3/r**3+c4/r**4
 RR,W,j,c,dr,er,HR,D,L,A1,A0,B1,B0=old.coeffs(F,r);Z=RR+12*lam
 h=c-(23*Z+6*r*S.diff(Z,r))/36;d0=E+4*lam
 Vtilde=F*(S.diff(F,r)/r-E+11*Z/36-c)
 zero('exact_charged_operator_difference',F*D-Vtilde-F*(d0+h))
 zero('gap_correction_outer_limit',S.limit(F*h,r,S.oo))
 zero('gap_correction_outer_leading',S.limit(r*F*h,r,S.oo)-R(14,3)*lam*m)
 mh=r*(lam*r*r+1+Q/r**2+c3/r**3+c4/r**4)/2
 zero('gap_correction_horizon_limit',S.factor(F*h).subs(m,mh))
 A,U,FF,Hh,Kq=S.symbols('A U FF Hh Kq')
 zero('resolvent_correction_sign',FF*U/6-(A+FF*(d0+Hh))*U/(6*d0)+(A+FF*Hh)*U/(6*d0))
 k=S.symbols('k',positive=True);a=k*k+R(1,4);alpha=R(1,2)+S.I*k
 zero('leading_reconstruction_denominator',6-alpha*(alpha-1)-(a+6))
 zero('trace_coefficient_matches_vacuum',6*lam*(a+6)-6*(lam*(a+2)+4*lam))
 y=S.symbols('y',positive=True);gg=S.Function('g')(y);b=S.symbols('b',positive=True)
 x=S.symbols('x',positive=True)
 z=(y**-2*S.Integral(x**3*S.Function('g')(x),(x,0,y))+y**3*S.Integral(x**-2*S.Function('g')(x),(x,y,b)))/5
 zero('local_Green_representation',-S.diff(z,y,2)+6*z/y**2-gg)
 zero('local_Green_derivative',S.diff(z,y)-(-2*y**-3*S.Integral(x**3*S.Function('g')(x),(x,0,y))+3*y*y*S.Integral(x**-2*S.Function('g')(x),(x,y,b)))/5)
 zero('head_Cauchy_weight',S.integrate(x**6,(x,0,y))-y**7/7)
 zero('tail_Cauchy_weight',S.integrate(x**-4,(x,y,b))-(y**-3-b**-3)/3)
 uc=y**alpha*S.exp(-y*y)
 zero('model_curvature_source',(-S.diff(uc,y,2)-a*uc/y**2+gamma**2*uc)-uc*(4+gamma**2-4*y*y+4*S.I*k))
 power=S.symbols('power',positive=True);w=y**power*S.exp(-y*y)
 Bw=-S.diff(w,y,2)-2*S.diff(w,y)/y+(gamma**2+12/y**2)*w
 zero('model_shape_source',Bw-w*((12-power*(power+1))/y**2+4*power+6+gamma**2-4*y*y))
 zero('model_shape_leading_coefficient',-(12-power*(power+1))/(3*(6-power*(power-1)))+(power+4)/(3*(power+2)))
 check('shape_control_H1',2*R(5,4)-2>-1)
 check('shape_control_fourth_form',2*R(5,4)-2>-1)
 zero('metric_second_derivative',S.diff(y**(R(3,2)+S.I*k),y,2)-(R(3,2)+S.I*k)*(R(1,2)+S.I*k)*y**(-R(1,2)+S.I*k))
 check('nonzero_metric_second_derivative_coefficient',S.expand_complex(abs((R(3,2)+S.I*k)*(R(1,2)+S.I*k))**2)>0)
 zero('diagonal_scalar_trace_Wronskian',S.sqrt(y)*S.cos(k*S.log(y))*S.diff(S.sqrt(y)*S.sin(k*S.log(y)),y)-S.diff(S.sqrt(y)*S.cos(k*S.log(y)),y)*S.sqrt(y)*S.sin(k*S.log(y))-k)

def green_value(source,y,gamma):
 aa=lambda x:mp.sqrt(x)*mp.besseli(mp.mpf('2.5'),gamma*x)
 bb=lambda x:mp.sqrt(x)*mp.besselk(mp.mpf('2.5'),gamma*x)
 head=mp.quad(lambda x:aa(x)*source(x),[0,y/16,y/4,y])
 tail=mp.quad(lambda x:bb(x)*source(x),[y,mp.mpf('.1'),1,4,mp.inf])
 return bb(y)*head+aa(y)*tail

def model_fixtures():
 for ix,row in enumerate(PROTOCOL['model_cases']):
  levels=[]
  for precision in PROTOCOL['precision_digits']:
   with mp.workdps(precision):
    lam=mp.mpf(PROTOCOL['model_lambda']);k,gamma=map(mp.mpf,row);a=k*k+mp.mpf('.25');d0=lam*(a+6);coef=1/(6*d0)
    power=mp.mpf(PROTOCOL['shape_exponent'])
    U=lambda y:mp.sqrt(y)*mp.exp(-y*y)*mp.cos(k*mp.log(y))
    hsrc=lambda y:mp.sqrt(y)*mp.exp(-y*y)*((4+gamma*gamma-4*y*y)*mp.cos(k*mp.log(y))-4*k*mp.sin(k*mp.log(y)))
    shape=lambda y:coef*mp.mpf(PROTOCOL['shape_amplitude_in_units_of_trace_coefficient'])*y**power*mp.exp(-y*y)
    Bshape=lambda y:shape(y)*((12-power*(power+1))/y**2+4*power+6+gamma*gamma-4*y*y)
    vals=[]
    for jy,ys in enumerate(PROTOCOL['observation_points']):
     y=mp.mpf(ys);direct=green_value(lambda s:U(s)/(6*lam*s*s),y,gamma)
     correction=-coef*green_value(hsrc,y,gamma);rhs=coef*U(y)+correction
     qshape=-green_value(Bshape,y,gamma)/3
     scale=coef*mp.sqrt(y);err=abs(direct-rhs)/scale
     tag=str(ix)+'_'+str(precision)+'_'+str(jy)
     check('Green_trace_identity_'+tag,err<mp.mpf(PROTOCOL['green_identity_relative_tolerance']),mp.nstr(err,15))
     curvature_error=abs(direct/scale-mp.cos(k*mp.log(y)))
     metric_error=abs((direct+qshape)/scale-mp.cos(k*mp.log(y)))
     check('nonzero_shape_tail_'+tag,abs(qshape/scale)>mp.mpf('1e-6'),mp.nstr(qshape/scale,15))
     vals.append({'y':ys,'q_curvature':mp.nstr(direct,precision),'q_correction':mp.nstr(correction,precision),
      'q_shape':mp.nstr(qshape,precision),'identity_error':mp.nstr(err,20),'curvature_trace_error':mp.nstr(curvature_error,20),
      'full_metric_trace_error':mp.nstr(metric_error,20),'shape_scaled_tail':mp.nstr(qshape/scale,20)})
    check('fine_curvature_trace_'+str(ix)+'_'+str(precision),mp.mpf(vals[-1]['curvature_trace_error'])<mp.mpf(PROTOCOL['fine_curvature_trace_error_limit']))
    check('fine_full_metric_trace_'+str(ix)+'_'+str(precision),mp.mpf(vals[-1]['full_metric_trace_error'])<mp.mpf(PROTOCOL['fine_full_metric_trace_error_limit']))
    check('full_metric_trace_refines_'+str(ix)+'_'+str(precision),all(mp.mpf(vals[j+1]['full_metric_trace_error'])<mp.mpf(PROTOCOL['full_metric_error_refinement_ratio'])*mp.mpf(vals[j]['full_metric_trace_error']) for j in range(2)))
    levels.append({'precision':precision,'values':vals})
  with mp.workdps(80):
   for jy,ys in enumerate(PROTOCOL['observation_points']):
    a0=levels[0]['values'][jy];a1=levels[1]['values'][jy]
    err=max(abs(mp.mpf(a0[key])-mp.mpf(a1[key]))/(1+abs(mp.mpf(a1[key]))) for key in ['q_curvature','q_correction','q_shape'])
    check('precision_refinement_'+str(ix)+'_'+str(jy),err<mp.mpf(PROTOCOL['precision_relative_tolerance']),mp.nstr(err,20))
  FIXTURES.append({'fixture':'inverse_square_reconstruction_trace','k':row[0],'gamma':row[1],'levels':levels,'actual_charged_kernel':False})

def domain_controls():
 cases=[S.Matrix([[1,0,0],[0,1,0],[0,0,1],[1,-2,3]]),S.Matrix([[1,0,0],[0,1,0],[1,1,0],[2,-1,0]])]
 for ix,A in enumerate(cases):
  zs=A.T.nullspace();columns=[zs[0],zs[1] if len(zs)>1 else S.Matrix([1,0,0,0])]
  values=[];M=A
  for d in range(3):
   if d:M=M.row_join(columns[d-1])
   rank=M.rank();ker=M.cols-rank;coker=M.rows-rank
   check('domain_index_'+str(ix)+'_'+str(d),ker-coker==d-1,[ker,coker])
   if d==2:check('maximal_domain_has_kernel_'+str(ix),ker>=1)
   if ix==1 and d==1:check('index_zero_not_invertibility',ker==coker==1)
   values.append({'trace_dimension':d,'kernel_dimension':ker,'cokernel_dimension':coker,'index':ker-coker})
  FIXTURES.append({'fixture':'infinite_shift_domain_extension','case':ix,'domains':values,'tail':'unchanged isometry from input coordinates >=3 to output coordinates >=4','charged_operator':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 coupled fast metric boundary traces','input_count':len(ps)}),flush=True)
 try:provenance();exact_checks();model_fixtures();domain_controls()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 coupled fast metric boundary traces','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'exact_maximal_domain_reconstruction_identity_derived':good,'fast_metric_trace_transfer_derived':good,
  'shape_cannot_cancel_leading_fast_trace_derived':good,'leading_boundary_metric_insufficient_for_trace_selection_derived':good,
  'auxiliary_domain_indices_minus_one_zero_plus_one_derived':good,'auxiliary_maximal_domain_kernel_exists_derived':good,
  'physical_boundary_domain_selected':False,'full_boundary_variational_action_evaluated':False,
  'physical_matter_domain_for_all_fast_traces_proved':False,'physical_metric_mode_constructed':False,
  'continuous_causal_inverse_proved':False,'actual_shape_kernel_evaluated':False,'actual_charged_stability_proved':False,
  'actual_wormhole_instability_proved':False,'physical_curvature_norm_changed':False,'finite_counterterm_changed':False,
  'order_reduction_adopted':False,'conformal_contour_rotation_adopted':False,
  'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in [('checker',OWN[0]),('note',OWN[1]),('review',OWN[2])]:run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
