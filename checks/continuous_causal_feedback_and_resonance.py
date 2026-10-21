"""Common-form continuity, logarithmic transport and actual-source algebra."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='continuous-causal-feedback-and-resonance'
PREVIOUS='global-continuous-causal-shape-inverse'
OWN=['checks/continuous_causal_feedback_and_resonance.py','notes/'+STEM+'.md','evidence/continuous-causal-feedback-source-review.json']
ROWS=[];FIXTURES=[]
PROTOCOL={'precision_digits':[50,80],'identity_tolerance':'1e-28','precision_tolerance':'1e-28',
 'spectral_anchors':['10','14','18','22'],'spectral_s':['12.3','19.1'],
 'spectral_measure':'t^2/(1+t) times [[1,1/(1+t)],[1/(1+t),1+1/(1+t)^2]] dt',
 'subtraction_point':'1','finite_difference_steps':['.002','.001'],
 'finite_difference_tolerance':'1e-6','refinement_ratio_ceiling':'.3',
 'transport_A':[['3','-1'],['-1','2']],'transport_B':[['1','0'],['0','80']],
 'transport_s':['100','180'],'transport_bound':'4/s',
 'coupled_s':['4.2','9.7'],'coupled_G':'1','coupled_epsilon':'.2',
 'coupled_K0':[['9','-1','.2'],['-1','12','-.5'],['.2','-.5','16']],
 'coupled_B_offset':[['1','.6','0'],['-.3','1','.2'],['.1','0','2']],
 'coupled_F':['2','6','15'],'canonical_A':[['2','-1','0'],['-1','3','-1'],['0','-1','2']],
 'canonical_B':['1','3','9'],'canonical_W':['.3','.5','.2'],
 'mixed_D0':[['.2','.1','0'],['0','.3','-.1'],['.1','0','.4']],
 'mixed_D1':['.02','.03','.01'],'mixed_D2':[['.4','-.1','0'],['.2','.1','.1'],['0','-.2','.3']],
 'curvature_C0':[['.2','.1','0'],['0','.1','.1'],['.1','0','.3']],
 'curvature_C1':[['1','.2','0'],['-.1','.7','.1'],['0','.3','.9']],
 'curvature_P':['.1','.2','.3'],'curvature_W':['.3','-.2','.1'],
 'curvature_base':['.5','.8','1.1'],'source_f':['.7','-.2','1'],'source_ell':['.1','.3','-.4'],
 'actual_charged_kernel_evaluated':False,'actual_numeric_threshold_certified':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)
def zero(name,value):
 value=S.simplify(S.sympify(value).doit());check(name,value==0,value)
def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})
def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==100 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes());check('review_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('scope',review['common_shape_form_domain_proved'] and review['actual_continuous_Schur_family_verified']
       and review['strong_feedback_locally_norm_Lipschitz_proved'] and review['fixed_phase_auxiliary_resonance_sequence_proved']
       and not review['physical_boundary_state_admissibility_proved'] and not review['charged_Lorentzian_instability_proved'] and not review['goal_complete'])

def exact_algebra():
 y,s=S.symbols('y s',real=True);a=S.Function('a')(y);v=S.Function('v')(y);B=a**-2
 A=lambda f:-S.diff(a*S.diff(f,y),y)/a;H=lambda f:A(f)+s*B*f
 P=lambda f:-2*S.diff(B,y)/B*S.diff(f,y)+(2*(S.diff(B,y)/B)**2+A(B)/B)*f
 zero('exact_multiplier_transport',H(B*v)-P(B*v)-B*H(v))
 zero('transport_first_order_coefficients',P(v)-4*S.diff(a,y)/a*S.diff(v,y)-(4*S.diff(a,y)**2+2*a*S.diff(a,y,2))*B*v)
 zero('bounded_inverse_multiplier_identity',H(v/B)-(H(v)-P(v))/B)
 x=S.symbols('x',nonnegative=True);z=x+1
 zero('dyadic_weight_endpoint_factor',9*z-2*(1+z)**2-(z-2)*(1-2*z))
 h,c=S.symbols('h c',positive=True)
 zero('central_derivative_remainder_constant',((s+h)**3-(s-h)**3)/(2*h)-3*s*s-6*h*h/6)
 zero('reconstruction_u_norm',S.Rational(4,5)+S.Rational(1,3)-S.Rational(17,15))
 F,r,Fr,HR=S.symbols('F r Fr HR',positive=True)
 ao=F/r;ay=-F*Fr/r+F**2/r**2;VK=F*(Fr/r-HR/6)
 zero('conformal_square_coefficient',VK+ay+ao**2-F*(-HR/6+2*F/r**2))
 zero('compensated_Weyl_contact',S.Rational(1,30)*(S.Rational(19,60)-S.Rational(1,6))-S.Rational(1,200))
 K0=S.Matrix([[4,1],[1,3]]);B0=S.Matrix([[5,2],[0,4]]);I=S.eye(2);R=(K0+s*I).inv();Q=-R*(B0+s*I)/3
 for i in range(2):
  for j in range(2):zero(f'exact_metric_reconstruction_derivative_{i}_{j}',S.diff(Q[i,j],s)+(R*(Q+I/3))[i,j])
 posF=S.diag(2,7);RR=R*posF/6
 for i in range(2):
  for j in range(2):zero(f'exact_weighted_reconstruction_derivative_{i}_{j}',S.diff(RR[i,j],s)+(R*RR)[i,j])
 check('smooth_bolt_transport_threshold',4/mp.sqrt(200)+mp.sqrt(2)*4/100<mp.mpf('.5'))

def parse(M):return mp.matrix([[mp.mpf(x) for x in row] for row in M])
def diag(xs):return mp.diag(list(map(mp.mpf,xs)))
def vec(xs):return mp.matrix(list(map(mp.mpf,xs)))
def eigmin(M):return min(mp.eigsy((M+M.T)/2,eigvals_only=True))
def opnorm(M):return mp.sqrt(max(mp.eigsy(M.T*M,eigvals_only=True)))
def matrixfun(M,fun):
 ev,U=mp.eigsy(M);return U*mp.diag([fun(x) for x in ev])*U.T
def matrixlog(M):return matrixfun(M,mp.log)
def serial(M,digits):return [[mp.nstr(M[i,j],digits) for j in range(M.cols)] for i in range(M.rows)]
def quadmatrix(fun,n,grid):return mp.matrix([[mp.quad(lambda t:fun(t)[i,j],grid) for j in range(n)] for i in range(n)])

def spectral_derivatives():
 for case,ss in enumerate(PROTOCOL['spectral_s']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    anchors=list(map(mp.mpf,PROTOCOL['spectral_anchors']));s=mp.mpf(ss);grid=[0,1,mp.inf]
    measure=lambda t:t*t/(1+t)*mp.matrix([[1,1/(1+t)],[1/(1+t),1+1/(1+t)**2]])
    poly=lambda x:mp.matrix([[1+x+x*x/3,mp.mpf('.2')+x/5],[mp.mpf('.2')+x/5,2-x+x*x/7]])
    polyp=lambda x:mp.matrix([[1+2*x/3,mp.mpf(1)/5],[mp.mpf(1)/5,-1+2*x/7]])
    value=lambda x:poly(x)+(x-1)**3*quadmatrix(lambda t:measure(t)/((t+x)*(t+1)**3),2,grid)
    dd=sum((value(x)/mp.fprod(x-y for y in anchors if y!=x) for x in anchors),mp.zeros(2))
    integral=quadmatrix(lambda t:measure(t)/mp.fprod(t+x for x in anchors),2,grid)
    err=mp.norm(dd-integral)/(1+mp.norm(integral));check(f'independent_subtracted_divided_difference_{case}_{digits}',err<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(err,digits))
    third=6*quadmatrix(lambda t:measure(t)/(t+s)**4,2,grid);ratio=mp.fprod(max(1,x/s) for x in anchors)
    gap=eigmin(6*ratio*integral-third);check(f'positive_measure_derivative_enclosure_{case}_{digits}',eigmin(third)>0 and gap>0,mp.nstr(gap,digits))
    errs=[];records=[]
    for hs in PROTOCOL['finite_difference_steps']:
     h=mp.mpf(hs);fd=(value(s+2*h)-2*value(s+h)+2*value(s-h)-value(s-2*h))/(2*h**3)
     err=mp.norm(fd-third)/(1+mp.norm(third));errs.append(err)
     check(f'spectral_third_finite_difference_{case}_{digits}_{hs}',err<mp.mpf(PROTOCOL['finite_difference_tolerance']),mp.nstr(err,digits))
     records.append({'step':hs,'scaled_error':mp.nstr(err,digits)})
    ratiofd=errs[1]/errs[0];check(f'spectral_refinement_{case}_{digits}',ratiofd<mp.mpf(PROTOCOL['refinement_ratio_ceiling']),mp.nstr(ratiofd,digits))
    h=mp.mpf('.1');first=polyp(s)+quadmatrix(lambda t:measure(t)/(t+1)**3*(3*(s-1)**2/(t+s)-(s-1)**3/(t+s)**2),2,grid)
    centered=(value(s+h)-value(s-h))/(2*h);bound=6*mp.fprod(max(1,x/(s-h)) for x in anchors)*opnorm(integral)*h*h/6
    rr=opnorm(centered-first)/bound;check(f'central_Taylor_operator_bound_{case}_{digits}',0<rr<1,mp.nstr(rr,digits))
    levels.append({'precision':digits,'third_derivative':serial(third,digits),'derivative_gap_min':mp.nstr(gap,digits),'finite_differences':records})
  with mp.workdps(80):
   aa=parse(levels[0]['third_derivative']);bb=parse(levels[1]['third_derivative']);err=mp.norm(aa-bb)/(1+mp.norm(bb))
   check('spectral_precision_'+str(case),err<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(err,80))
  FIXTURES.append({'fixture':'positive_subtracted_measure','s':ss,'levels':levels,'charged_spectral_measure':False})

def logarithmic_transport():
 for case,ss in enumerate(PROTOCOL['transport_s']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    A=parse(PROTOCOL['transport_A']);B=parse(PROTOCOL['transport_B']);s=mp.mpf(ss);H=A+s*B;T=B*H**-1;U=T/(4/s);I=mp.eye(2)
    check(f'noncommuting_transport_{case}_{digits}',mp.norm(H*B-B*H)>1)
    check(f'transport_ordinary_contraction_{case}_{digits}',opnorm(U)<1,mp.nstr(opnorm(U),digits))
    check(f'transport_graph_contraction_{case}_{digits}',opnorm(H*U*H**-1)<1,mp.nstr(opnorm(H*U*H**-1),digits))
    L=matrixlog(I+H*H);integral=quadmatrix(lambda t:I/(1+t)-(I+H*H+t*I)**-1,2,[0,1,mp.mpf('1e6'),mp.inf])
    err=mp.norm(L-integral)/(1+mp.norm(L));check(f'independent_log_transport_integral_{case}_{digits}',err<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(err,digits))
    jgap=eigmin(matrixlog(I+U.T*H*H*U)-U.T*L*U);gap=eigmin(L-U.T*L*U)
    check(f'contractive_log_Jensen_gap_{case}_{digits}',jgap>0,mp.nstr(jgap,digits))
    check(f'logarithmic_transport_gap_{case}_{digits}',gap>0,mp.nstr(gap,digits))
    levels.append({'precision':digits,'logarithm':serial(L,digits),'Jensen_gap_min':mp.nstr(jgap,digits),'transport_gap_min':mp.nstr(gap,digits)})
  with mp.workdps(80):
   aa=parse(levels[0]['logarithm']);bb=parse(levels[1]['logarithm']);err=mp.norm(aa-bb)/(1+mp.norm(bb))
   check('transport_precision_'+str(case),err<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(err,80))
  FIXTURES.append({'fixture':'noncommuting_logarithmic_transport','s':ss,'levels':levels,'charged_operator':False})

def coupled_continuity():
 for case,ss in enumerate(PROTOCOL['coupled_s']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    K0=parse(PROTOCOL['coupled_K0']);B0=K0+parse(PROTOCOL['coupled_B_offset']);F=diag(PROTOCOL['coupled_F']);I=mp.eye(3)
    Ac=parse(PROTOCOL['canonical_A']);Bc=diag(PROTOCOL['canonical_B']);Wc=diag(PROTOCOL['canonical_W'])
    D0=parse(PROTOCOL['mixed_D0']);D1=diag(PROTOCOL['mixed_D1']);D2=parse(PROTOCOL['mixed_D2'])
    C0=parse(PROTOCOL['curvature_C0']);C1=parse(PROTOCOL['curvature_C1']);P=diag(PROTOCOL['curvature_P']);W=diag(PROTOCOL['curvature_W'])
    base=diag(PROTOCOL['curvature_base']);eps=mp.mpf(PROTOCOL['coupled_epsilon']);G=mp.mpf(PROTOCOL['coupled_G']);s=mp.mpf(ss)
    def build(x):
     Kinv=(K0+x*I)**-1;Q=-Kinv*(B0+x*I)/3;R=Kinv*F/6;H=Ac+x*Bc;L=matrixlog(I+H)
     AA=H*L*H-H*Wc*H+Q.T*Q+(1+x)*I;D=D0+x*D1+D2*Q;b=-R.T*D/(2*G);Kd=AA**-1*b.T
     CC=C0+C1*R;V=W+CC-P*R/eps+D*Kd/eps
     return {'Kinv':Kinv,'Q':Q,'R':R,'H':H,'L':L,'A':AA,'D':D,'b':b,'Kd':Kd,'V':V,'L0':base+x*I+W+CC-P*R/eps}
    z=build(s);check(f'coupled_positive_shape_{case}_{digits}',eigmin(z['A'])>0,mp.nstr(eigmin(z['A']),digits))
    check(f'coupled_strong_feedback_nonnormal_{case}_{digits}',mp.norm(z['V']-z['V'].T)>mp.mpf('.01'))
    err=mp.norm(z['A']*z['Kd']-z['b'].T)/(1+mp.norm(z['b']));check(f'energy_adjoint_identity_{case}_{digits}',err<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(err,digits))
    ref=build(s+mp.mpf('.3'))['A'];Rhalf=matrixfun(ref,mp.sqrt);Rinv=Rhalf**-1;As=Rinv*z['A']*Rinv
    common=Rinv*As**-1*Rinv*z['b'].T;err=mp.norm(common-z['Kd'])/(1+mp.norm(z['Kd']))
    check(f'common_energy_representation_{case}_{digits}',err<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(err,digits))
    block=mp.matrix(6);f=vec(PROTOCOL['source_f']);ell=vec(PROTOCOL['source_ell'])
    for i in range(3):
     for j in range(3):
      block[i,j]=z['L0'][i,j];block[i,j+3]=-z['D'][i,j]/eps;block[i+3,j]=z['b'].T[i,j];block[i+3,j+3]=z['A'][i,j]
    full=mp.lu_solve(block,mp.matrix(list(f)+list(ell)));wl=z['A']**-1*ell
    UU=mp.lu_solve(z['L0']+z['D']*z['Kd']/eps,f+z['D']*wl/eps);ww=wl-z['Kd']*UU
    err=mp.norm(full-mp.matrix(list(UU)+list(ww)))/(1+mp.norm(full))
    check(f'independent_full_block_elimination_{case}_{digits}',err<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(err,digits))
    Qp=-z['Kinv']*(z['Q']+I/3);Rp=-z['Kinv']*z['R'];Dp=D1+D2*Qp;bp=-(Rp.T*z['D']+z['R'].T*Dp)/(2*G)
    Hp=Bc;Lp=quadmatrix(lambda t:(I+z['H']+t*I)**-1*Bc*(I+z['H']+t*I)**-1,3,[0,1,mp.inf])
    Ap=Hp*z['L']*z['H']+z['H']*Lp*z['H']+z['H']*z['L']*Hp-Hp*Wc*z['H']-z['H']*Wc*Hp+Qp.T*z['Q']+z['Q'].T*Qp+I
    Kdp=-z['A']**-1*Ap*z['Kd']+z['A']**-1*bp.T;Vp=C1*Rp-P*Rp/eps+(Dp*z['Kd']+z['D']*Kdp)/eps
    errors=[];records=[]
    for hs in PROTOCOL['finite_difference_steps']:
     h=mp.mpf(hs);fd=(build(s+h)['V']-build(s-h)['V'])/(2*h);err=mp.norm(fd-Vp)/(1+mp.norm(Vp));errors.append(err)
     check(f'coupled_feedback_derivative_{case}_{digits}_{hs}',err<mp.mpf(PROTOCOL['finite_difference_tolerance']),mp.nstr(err,digits))
     records.append({'step':hs,'scaled_error':mp.nstr(err,digits)})
    ratio=errors[1]/errors[0];check(f'coupled_feedback_refinement_{case}_{digits}',ratio<mp.mpf(PROTOCOL['refinement_ratio_ceiling']),mp.nstr(ratio,digits))
    wrong=z['b']*z['Kd'];check(f'strong_variational_feedback_distinguished_{case}_{digits}',mp.norm(z['D']*z['Kd']-wrong)>mp.mpf('1e-12'))
    levels.append({'precision':digits,'V_derivative':serial(Vp,digits),'finite_differences':records})
  with mp.workdps(80):
   aa=parse(levels[0]['V_derivative']);bb=parse(levels[1]['V_derivative']);err=mp.norm(aa-bb)/(1+mp.norm(bb))
   check('coupled_precision_'+str(case),err<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(err,80))
  FIXTURES.append({'fixture':'nonnormal_coupled_form_continuity','s':ss,'levels':levels,'charged_operator':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 continuous causal feedback and resonance','input_count':len(ps)}),flush=True)
 try:provenance();exact_algebra();spectral_derivatives();logarithmic_transport();coupled_continuity()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 continuous causal feedback and resonance','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'common_shape_form_domain_proved':good,'actual_continuous_Schur_family_verified':good,
  'strong_feedback_locally_norm_Lipschitz_proved':good,'fixed_phase_auxiliary_resonance_sequence_proved':good,
  'analytic_poles_or_simple_roots_proved':False,'actual_numeric_threshold_certified':False,'physical_boundary_state_admissibility_proved':False,
  'actual_charged_two_point_kernel_evaluated':False,'physical_fast_branch_excluded':False,'finite_counterterm_changed':False,
  'new_boundary_prescription_adopted':False,'physical_source_projection_adopted':False,'order_reduction_adopted':False,
  'physical_exterior_state_constructed':False,'continuous_causal_stability_proved':False,'charged_Lorentzian_instability_proved':False,
  'formal_verification':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]};data['runs'].append(run)
 dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n');print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
