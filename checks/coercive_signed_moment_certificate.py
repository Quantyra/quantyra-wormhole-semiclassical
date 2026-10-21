"""Full coupled moment transfer, positive solves and explicit error budgets."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='coercive-signed-moment-certificate'
PREVIOUS='seed-residual-constraint-obstruction'
OWN=['checks/coercive_signed_moment_certificate.py','notes/'+STEM+'.md',
     'evidence/coercive-signed-moment-source-review.json']
ROWS=[];FIXTURES=[]
PROTOCOL={'precision_digits':[50,80],'identity_tolerance':'1e-28','steps':3,
 'E_lower':'1.125e35','cold_coefficient_target':'1/(3e43)',
 'residual_norm_targets':['2e-43','2e-8','2e27'],'terminal_norm_target':'2e62',
 'exact_terminal_norm_target':'1e63','moment_norm_constant':'1/25',
 'trial_error_amplitudes':['.003','.006','.009'],
 'T_cases':[[['4','-.3','.2'],['-.3','11','-.5'],['.2','-.5','24']],
            [['9','.4','-.6'],['.4','31','.8'],['-.6','.8','67']]],
 'C_cases':[[['2','.5'],['.3','1'],['1.1','-.6']],
            [['8','1'],['-2','4'],['3','-1']]],
 'TT_D':[['4','.4'],['.4','6']],'local_vertex':['.02','.03','-.01'],
 'source':['1','.3','-.2','.4','-.7'],'volume_ratio':['1.1','.9','1.2','.8','1.3'],
 'lambda':'1.7','weight':['1','2','3'],'trial_error_vector':['.2','-.4','.7'],
 'actual_source_evaluated':False,'actual_certificate_passed':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,value):
 value=S.simplify(S.sympify(value).doit());check(name,value==0,value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==111 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes());check('review_inputs_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('scope',review['full_coupled_moment_certificate_proved'] and not review['actual_certificate_passed']
       and not review['actual_signed_coefficient_bounded_to_cold_target'] and not review['goal_complete'])

def certificate_bound(E_lower,residual_norms,terminal_norm):
 """Combine certified normalized norms; this does not certify the inputs."""
 E=S.Rational(E_lower);rr=[S.Rational(x) for x in residual_norms];last=S.Rational(terminal_norm)
 if E<=0 or last<0 or any(x<0 for x in rr):raise ValueError('Positive E and nonnegative norm bounds required')
 return (sum(v/E**j for j,v in enumerate(rr))+last/E**len(rr))/25

def exact_budgets():
 R=S.Rational;xi=R('30.1')*R('6e-4')+R('1e-8')
 gt=R('1.026')*xi+R('3e-7')*2520;gs=R('4e-7')*xi+R('1.000002')*2520
 err=R('1.5')*gt+R('3e-8')*gs+R('.0013')
 check('scalar_adjoint_error_budget',err<R('.031'),err)
 check('actual_volume_scalar_norm_budget',R('1.000001')*(4+err)<R('4.04'))
 # pi>3.1 and sqrt(8.7)>2.9 already prove the rounded constant.
 constant=R('4.04')/(80*R('3.1')**2*R('.048')*R('2.9'))
 check('moment_norm_constant',constant<R(1,25),constant)
 check('E_divisor_local_form_bound',R('2.6e-13')*R('1.000001')<R('3e-13'))
 E0=R('180')/R('4e-17')**2;zero('benchmark_E_lower',E0-R(PROTOCOL['E_lower']))
 target=1/R('3e43');bound=certificate_bound(PROTOCOL['E_lower'],PROTOCOL['residual_norm_targets'],PROTOCOL['terminal_norm_target'])
 check('inexact_three_step_cold_budget',bound<R('2.71e-44') and bound<target,str(bound))
 exact=certificate_bound(PROTOCOL['E_lower'],['0','0','0'],PROTOCOL['exact_terminal_norm_target'])
 check('exact_three_step_cold_budget',exact<R('2.81e-44') and exact<target,str(exact))
 failed=certificate_bound(PROTOCOL['E_lower'],['1e-40','0','0'],'0')
 check('loose_early_residual_rejected',failed>target,str(failed))
 check('old_absolute_bound_insufficient',R('5.3e21')>target)
 E=S.symbols('E',positive=True);q0,q1,q2,q3,r0,r1,r2,a0,a1,a2=S.symbols('q0 q1 q2 q3 r0 r1 r2 a0 a1 a2')
 expression=q0-q3/E**3-r0-r1/E-r2/E**2+(a0+a1/E+a2/E**2)/E
 zero('inexact_telescoping_three_steps',expression.subs({q0:q1/E+r0-a0/E}).subs({q1:q2/E+r1-a1/E}).subs({q2:q3/E+r2-a2/E}))
 n=S.symbols('n',positive=True);zero('uncoupled_filter_unbounded_symbol',n*(n+4)/(n+4)-n)
 zero('resonant_transfer_no_gain',(E**3)/E**3-1)
 lam,volume,ft,fc,B,D,psi=S.symbols('lambda volume f_tau f_c B D psi',nonzero=True)
 cpsi=-B*psi/D;p=80*S.pi**2*(ft-B*fc/D)
 zero('physical_full_source_normalization',lam**2*(psi*ft+cpsi*fc)-lam**2*psi*p/(80*S.pi**2))
 zero('fixed_density_before_transpose',volume*(ft/volume)-ft)
 for bad in [('0',['0'],'0'),('1',['-1'],'0'),('1',['0'],'-1')]:
  rejected=False
  try:certificate_bound(*bad)
  except ValueError:rejected=True
  check('invalid_certificate_input_'+str(bad),rejected)
 FIXTURES.append({'fixture':'exact_charged_budget','inexact_bound':str(bound),'exact_bound':str(exact),
  'cold_target':str(target),'actual_source_certificate':False})

def mat(rows):return mp.matrix([[mp.mpf(x) for x in row] for row in rows])
def vec(xs):return mp.matrix(list(map(mp.mpf,xs)))
def fmat(A,fn):
 vals,U=mp.eigsy(A);return U*mp.diag([fn(x) for x in vals])*U.T
def scaled_error(a,b):return mp.norm(a-b)/(1+mp.norm(b))

def noncommuting_blocks():
 for case in range(2):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    tol=mp.mpf(PROTOCOL['identity_tolerance']);T=mat(PROTOCOL['T_cases'][case]);H=T+4*mp.eye(3);D0=T*H
    C=mat(PROTOCOL['C_cases'][case]);Dt=mat(PROTOCOL['TT_D']);vs=vec(PROTOCOL['local_vertex']);signed=C*Dt**-1*C.T-vs*vs.T
    Hm=fmat(H,lambda x:1/mp.sqrt(x));vals,U=mp.eigsy(Hm*(D0-signed)*Hm);E=vals[2];psi=Hm*U[:,2]
    K=H+signed/E;P=D0-E*K;cpsi=-Dt**-1*C.T*psi;chi=mp.matrix(list(psi)+list(cpsi))
    check(f'positive_K_{case}_{digits}',min(mp.eigsy(K,eigvals_only=True))>0)
    check(f'noncommuting_filter_{case}_{digits}',mp.norm(D0*K-K*D0)>1,mp.nstr(mp.norm(D0*K-K*D0),digits))
    # Full tensor block retains the independent TT source and its response.
    block=mp.matrix(5);AA=D0-E*H+vs*vs.T
    for i in range(3):
     for j in range(3):block[i,j]=AA[i,j]
     for j in range(2):block[i,j+3]=C[i,j];block[j+3,i]=C[i,j]
    for i in range(2):
     for j in range(2):block[i+3,j+3]=Dt[i,j]
    err=mp.norm(block*chi)/(1+mp.norm(block)*mp.norm(chi));check(f'full_block_adjoint_{case}_{digits}',err<tol,mp.nstr(err,digits))
    Rmap=mp.eye(5);Rmap[0,3]=mp.mpf('.2');Rmap[1,4]=mp.mpf('-.3');Rmap[3,2]=mp.mpf('.1')
    Jvol=mp.diag(list(map(mp.mpf,PROTOCOL['volume_ratio'])));seed=vec(PROTOCOL['source']);lam=mp.mpf(PROTOCOL['lambda'])
    converted=Rmap.T*(Jvol**-1*seed)/lam**2;fs=converted[:3,0];fc=converted[3:,0];q0=fs-C*Dt**-1*fc
    moment=(psi.T*q0)[0];direct=(chi.T*converted)[0];density=((Rmap*chi).T*(Jvol**-1*seed))[0]/lam**2
    for name,value in [('full_Schur_source',direct),('fixed_density_source',density)]:
     err=abs(value-moment)/(1+abs(moment));check(f'{name}_{case}_{digits}',err<tol,mp.nstr(err,digits))
    omission=abs((psi.T*fs)[0]-moment);check(f'TT_source_omission_detected_{case}_{digits}',omission>mp.mpf('1e-5'),mp.nstr(omission,digits))
    wrongdensity=Rmap.T*seed/lam**2;wrong=(chi.T*wrongdensity)[0]
    check(f'density_omission_detected_{case}_{digits}',abs(wrong-moment)>mp.mpf('1e-5'))
    # Independent block solve for the positive Schur system E K.
    positiveblock=mp.matrix(5);upper=E*H-vs*vs.T
    for i in range(3):
     for j in range(3):positiveblock[i,j]=upper[i,j]
     for j in range(2):positiveblock[i,j+3]=C[i,j];positiveblock[j+3,i]=C[i,j]
    for i in range(2):
     for j in range(2):positiveblock[i+3,j+3]=-Dt[i,j]
    q=q0;residuals=[];steps=[];ww=mp.diag(list(map(mp.mpf,PROTOCOL['weight'])));errvec=vec(PROTOCOL['trial_error_vector'])
    for j in range(3):
     full=mp.lu_solve(positiveblock,mp.matrix(list(E*q)+[0,0]));sol=mp.lu_solve(K,q)
     err=scaled_error(full[:3,0],sol);check(f'independent_positive_block_solve_{case}_{digits}_{j}',err<tol,mp.nstr(err,digits))
     du=mp.mpf(PROTOCOL['trial_error_amplitudes'][j])*errvec;trial=sol+du;r=q-K*trial
     err=scaled_error(r,-K*du);check(f'actual_trial_residual_{case}_{digits}_{j}',err<tol,mp.nstr(err,digits))
     nextq=D0*trial
     defect=q-nextq/E-r+P*trial/E;err=mp.norm(defect)/(1+mp.norm(q));check(f'transfer_vector_identity_{case}_{digits}_{j}',err<tol,mp.nstr(err,digits))
     residuals.append(r);steps.append({'residual_norm':mp.nstr(mp.norm(ww*r),digits),'next_source_norm':mp.nstr(mp.norm(ww*nextq),digits)})
     q=nextq
    rhs=sum((psi.T*r)[0]/E**j for j,r in enumerate(residuals))+(psi.T*q)[0]/E**3
    err=abs(rhs-moment)/(1+abs(moment));check(f'inexact_full_moment_identity_{case}_{digits}',err<tol,mp.nstr(err,digits))
    bound=mp.norm(ww**-1*psi)*(sum(mp.norm(ww*r)/E**j for j,r in enumerate(residuals))+mp.norm(ww*q)/E**3)
    check(f'inexact_weighted_moment_bound_{case}_{digits}',bound>=abs(moment),mp.nstr(bound-abs(moment),digits))
    omitted=(psi.T*q)[0]/E**3;check(f'solve_error_omission_detected_{case}_{digits}',abs(omitted-moment)>mp.mpf('1e-6'),mp.nstr(abs(omitted-moment),digits))
    wrongq=q0
    for j in range(3):wrongq=mp.lu_solve(K,D0*wrongq)
    wrongmoment=(psi.T*wrongq)[0]/E**3
    check(f'noncommuting_order_error_detected_{case}_{digits}',abs(wrongmoment-moment)>mp.mpf('1e-6'),mp.nstr(abs(wrongmoment-moment),digits))
    resonant=K*psi;filtered=D0*mp.lu_solve(K,resonant)
    err=scaled_error(filtered,E*resonant);check(f'resonant_control_no_suppression_{case}_{digits}',err<tol,mp.nstr(err,digits))
    levels.append({'precision':digits,'E':mp.nstr(E,digits),'direct_moment':mp.nstr(moment,digits),
                   'transferred_moment':mp.nstr(rhs,digits),'weighted_bound':mp.nstr(bound,digits),'steps':steps})
  with mp.workdps(80):
   err=abs(mp.mpf(levels[0]['direct_moment'])-mp.mpf(levels[1]['direct_moment']))/(1+abs(mp.mpf(levels[1]['direct_moment'])))
   check('block_precision_'+str(case),err<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(err,80))
  FIXTURES.append({'fixture':'noncommuting_full_block_moment','case':case,'levels':levels,'charged_operator':False,'actual_source':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S167 coercive signed moment certificate','input_count':len(ps)}),flush=True)
 try:provenance();exact_budgets();noncommuting_blocks()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S167 coercive signed moment certificate','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'full_coupled_moment_certificate_proved':good,'positive_auxiliary_realization_proved':good,
  'inexact_three_step_cold_budget_proved':good,'actual_source_evaluated':False,'actual_certificate_passed':False,
  'actual_signed_coefficient_bounded_to_cold_target':False,'actual_signed_compatibility_zero_proved':False,
  'actual_weighted_graph_iterates_certified':False,'physical_fast_branch_excluded':False,
  'nonlinear_neighborhood_preserved':False,'nonlinear_solution_constructed':False,
  'finite_counterterm_changed':False,'new_boundary_prescription_adopted':False,'physical_source_projection_changed':False,
  'order_reduction_adopted':False,'physical_exterior_state_constructed':False,'physical_stability_proved':False,
  'formal_verification':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]};data['runs'].append(run)
 dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
