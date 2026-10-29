"""Complete localized kernel return map; actual varied source inputs remain unknown."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json, sys, traceback
sys.path.insert(0,str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb,ctx,dec,enc
ROOT=Path(__file__).resolve().parents[1]
STEM='current-localized-kernel-return-attempt'
INPUTS=['notes/'+STEM+'.md','checks/current_localized_kernel_return_attempt.py',
 'checks/current_quartic_ritz_rank_test.py','notes/current-disk-heat-propagation.md',
 'checks/current_disk_heat_propagation.py','evidence/current-disk-heat-propagation.json',
 'notes/current-localized-short-time-trial.md','notes/current-weighted-paired-propagation-test.md',
 'notes/current-localized-boundary-layer-test.md','notes/current-quartic-resolvent-dual-test.md',
 'evidence/current-global-green-trial-enclosure.json',
 'evidence/current-localized-short-time-trial.json',
 'evidence/current-weighted-paired-propagation-test.json',
 'evidence/current-localized-boundary-layer-test.json',
 'evidence/current-positive-time-quantum-operator.json',
 'evidence/current-first-quantum-solve-attempt.json',
 'evidence/current-first-inverse-source-screen.json',
 'evidence/current-coupled-mass-adjoint-runtime.json']
BASELINE=Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
KEYS=['A_delta','B_delta','A_deltaD','B_deltaD','a_delta','b_delta','a_deltaD','b_deltaD']

def digest(p):return sha256(p.read_bytes()).hexdigest()
def result(stem):
 d=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())
 r=d.get('runs',[d])[-1];assert r['status']!='failed'
 return r['result']
def pair(tab):return [dec(tab['inf']),dec(tab['row'])]
def uppermax(a,b):return arb(max(a.upper(),b.upper()))
def cv(p,f):return p[0]*f[1]/32+7*p[1]*f[0]/24

def execute():
 ctx.prec=512
 r=json.loads((ROOT/'evidence/current-disk-heat-propagation.json').read_bytes())['runs'][-1]
 assert r['status']=='passed' and r['complete_actual_domain_remainder_bounded']
 rows={n['kind']:n for n in r['numerics']}
 raw=rows['full_uncut_transport_norms']['tables']
 cap=rows['conditional_complete_cap_kernel_table']['table']
 cut=rows['complete_actual_domain_propagation_budget']['cutoff_jet_norm']
 tails=rows['complete_commutator_and_thermal_tails']
 assert all(dec(v)<arb('1e-25') for v in tails['commutator'].values())
 assert dec(tails['full_thermal_image_guard'])<arb('1e-100')
 t0=arb('.001');P={};F={}
 for key in ['base','time']:
  pi=pair(raw['P'][key]);ca=pair(cap[key]);fi=pair(raw['F'][key])
  for j,w in enumerate(['inf','row']):
   pi[j]+=dec(cut[w])*t0**3
   fi[j]+=dec(cut[w])+arb('1e-25')+arb('1e-100')
  ca[0]+=arb('1e-100')
  P[key]=[uppermax(2*ca[0],pi[0]),ca[1]+pi[1]]
  F[key]=fi
 p={k:dec(v) for k,v in result('current-global-green-trial-enclosure')['trial_parameters'].items()}
 L,gamma,beta,C,B2=[p[k] for k in ['L','gamma','beta','C','B2']]
 alpha=3*uppermax(arb(1),(1-L)/gamma);c=(4*alpha).exp()
 w0=gamma/(gamma+L);Wnorm=(gamma*(1-w0**5)/5).sqrt()
 PW={k:[c*x for x in v] for k,v in P.items()}
 FW={k:[c*x for x in v] for k,v in F.items()}
 x0=t0**4/140*FW['base'][1]
 xT=t0**4/140*(FW['base'][1]+FW['time'][1])
 assert 0<x0<arb('.001')
 E0=cv(PW['base'],FW['base'])
 ET=E0+cv(PW['time'],FW['base'])+cv(PW['base'],FW['time'])
 def response(v):
  Pd,PdT,Fd,FdT=v[:2],v[2:4],v[4:6],v[6:8]
  Ed=cv(Pd,FW['base'])+cv(PW['base'],Fd)
  EdT=Ed+cv(PdT,FW['base'])+cv(PW['time'],Fd)+cv(Pd,FW['time'])+cv(PW['base'],FdT)
  xd=t0**4/140*Fd[1];xdT=t0**4/140*(Fd[1]+FdT[1])
  S0=1/(1-x0);ST=xT/(1-x0)**2
  Sd=xd/(1-x0)**2;SdT=xdT/(1-x0)**2+2*xd*xT/(1-x0)**3
  return 4*(EdT*S0+Ed*ST+ET*Sd+E0*SdT)
 Lcoef=[]
 for i in range(8):Lcoef.append(response([arb(int(j==i)) for j in range(8)]))
 assert response([arb(1)]*8).overlaps(sum(Lcoef,arb(0)))
 old=result('current-localized-short-time-trial')['constants']
 layer=result('current-localized-boundary-layer-test')['constants']
 plus=result('current-positive-time-quantum-operator')['constants']
 CLt=dec(plus['CL_t']);kp=dec(plus['total_operator_difference_upper'])
 Elower=dec(old['trial_Et_lower']);tau6=dec(old['attained_leading_column_upper'])
 tau5=dec(layer['explicit_layer_relative_column_upper'])
 req=dec(old['conditional_relative_column_screen'])
 pref=4*C*CLt/(B2*beta)*Wnorm/Elower*t0*t0/2
 cc=[pref*x for x in Lcoef]
 Ghat=dec(result('current-first-quantum-solve-attempt')['constants']['proved_source_dual_upper'])
 Kfirst=dec(result('current-first-inverse-source-screen')['constants']['first_residual_to_original_bolt_multiplier'])
 physical=Kfirst*Ghat/(1-kp)
 constants=dict(alpha=alpha,support_weight_factor=c,trivial_kernel_weight_factor=1/w0**3,
  original_W_L2=Wnorm,Volterra_base_ratio=x0,Volterra_time_ratio=xT,
  complete_return_coefficient_sum=sum(Lcoef,arb(0)),
  relative_column_coefficient_sum=sum(cc,arb(0)),
  conditional_simultaneous_kernel_input_ceiling_OTHER_ZERO=(req-tau6-tau5)/sum(cc,arb(0)),
  conditional_relative_column_screen=req,attained_leading_charge=tau6,attained_layer_charge=tau5,
  retained_original_first_error_multiplier=physical)
 return dict(precision_bits=ctx.prec,constants={k:enc(v) for k,v in constants.items()},
  actual_metric_source_kernel_inputs=None,actual_cap_remainder=None,actual_boundary_remainder=None,
  actual_numerical_application_error=None,
  input_order=KEYS,complete_return_coefficients={k:enc(v) for k,v in zip(KEYS,Lcoef)},
  original_relative_column_charge_coefficients={k:enc(v) for k,v in zip(KEYS,cc)},
  original_first_error_charge_coefficients={k:enc(physical*v) for k,v in zip(KEYS,cc)},
  scope=dict(actual_paired_domain=True,complete_unvaried_thermal_kernels=True,
   complete_Volterra_return=True,full_angular_return_charge=True,
   defined_conditional_kernel_to_energy_dual_map=True,
   actual_varied_source_bound=False,complete_actual_quantum_remainder=False,
   complete_signed_first_residual=False,new_full_admission_margin=False,
   nonlinear_existence_matching_validity_stability=False,goal_complete=False))

if __name__=='__main__':
 out=ROOT/'evidence'/(STEM+'.json');assert not out.exists()
 before={p:digest(ROOT/p) for p in INPUTS}
 protected=json.loads(BASELINE.read_bytes());assert all(digest(Path(p))==h for p,h in protected.items())
 rt=json.loads((ROOT/INPUTS[-1]).read_bytes());pkg=rt['packages']['flint']
 assert digest(Path(sys.executable))==rt['python_sha256']
 assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
 rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
  arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),python_sha256=rt['python_sha256'])
 try:
  rec['result']=execute()
  after={p:digest(ROOT/p) for p in INPUTS};assert before==after
  assert all(digest(Path(p))==h for p,h in protected.items())
  rec.update(status='localized_kernel_return_map_bounded_actual_quantum_source_unresolved',
   input_sha256_after=after,protected_inputs_unchanged=len(protected))
 except Exception:rec.update(status='failed',traceback=traceback.format_exc())
 rec['finished_utc']=datetime.now(timezone.utc).isoformat()
 out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
 print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
  constants={k:v['ball'] for k,v in rec.get('result',{}).get('constants',{}).items()},
  relative_coefficients={k:v['ball'] for k,v in rec.get('result',{}).get('original_relative_column_charge_coefficients',{}).items()}),indent=2))
 if rec['status']=='failed':raise SystemExit(1)
