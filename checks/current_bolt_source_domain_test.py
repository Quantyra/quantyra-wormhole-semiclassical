"""Actual full-column domain counterexample and complete residual comparison."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import sympy as sp
import flint
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from operators.signed_local_core import SignedLocalCore
STEM='current-bolt-source-domain-test'
OWN=['notes/'+STEM+'.md','checks/current_bolt_source_domain_test.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+['operators/signed_local_core.py','notes/current-finite-local-action-assembly.md',
 'notes/current-coupled-application-residual-bound.md','evidence/current-coupled-application-residual-bound.json',
 'notes/current-interpolated-complete-bolt-test.md','evidence/current-interpolated-complete-bolt-test.json',
 'evidence/current-localized-bolt-adjoint-sensitivity.json','notes/current-strong-quotient-reconstruction-test.md']
CHECKS=[]
def sh(p):return sha256(p.read_bytes()).hexdigest()
def check(n,ok):
 CHECKS.append({'name':n,'passed':bool(ok)});print(n+(': passed' if ok else ': FAILED'),flush=True)
 if not ok:raise AssertionError(n)
def enc(x):
 m,r,e=x.mid_rad_10exp(50)
 return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(x):
 if isinstance(x,dict):
  e=x['decimal_enclosure'];return arb(e['mid']+' +/- '+e['rad'])*arb(10)**int(e['exponent'])
 return arb(x.strip('[]'))
def prior(stem):
 r=json.loads((ROOT/'evidence'/f'{stem}.json').read_bytes())['runs'][-1]
 check('accepted_'+stem,r['status']=='passed')
 for p,h in r['input_sha256_before'].items():
  if p in INPUTS and p.startswith(('notes/','operators/')):check('unchanged_'+p,sh(ROOT/p)==h)
 return r['result']
def execute():
 ctx.prec=512
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('ordinary_argument_and_checker_frozen',all(sh(ROOT/p)==h for p,h in review['scientific_sha256'].items()))
 latest=prior('current-interpolated-complete-bolt-test')
 app=prior('current-coupled-application-residual-bound')
 pars=prior('current-localized-bolt-adjoint-sensitivity')['parameters']
 # Compute complete local rows independently of the leading-density argument.
 core=SignedLocalCore();t,L,ga,p=sp.symbols('t L gamma p',positive=True)
 subst={core.z:L-t,core.Q[0]:2*ga*t,core.Q[1]:-2*ga}
 jets=sp.symbols('q2:9')
 subst.update({core.Q[j]:jets[j-2] for j in range(2,9)})
 # Keep higher seed jets symbolic; the leading indicial coefficient
 # must be independent of them. Assemble each exact original row first.
 pieces={};rows=[sp.Integer(0),sp.Integer(0)]
 for allocation in ['Einstein','Maxwell','R2','Weyl']:
  op=core.original_A_coefficients(allocation);arows=[]
  for a in range(2):
   lead=sp.Integer(0)
   for k in range(5):
    c=sp.factor(op[a][0][k].subs(subst,simultaneous=True))
    if k<=1:continue
    # Each leading coefficient has at least k-2 vanishing powers.
    ck=sp.simplify(sp.limit(c/t**(k-2),t,0,dir='+'))
    check('finite_bolt_coefficient_'+allocation+'_'+str(a)+'_'+str(k),not ck.has(sp.oo,-sp.oo,sp.zoo,sp.nan))
    lead+=ck*(-L)**k*sp.ff(p,k)
   lead=sp.factor(lead);arows.append(str(lead));rows[a]+=lead
  pieces[allocation]=arows
 d=(96*sp.pi*core.cR-1/(60*sp.pi))*core.C/core.B2
 expected=[6*d,2*d]
 for a in range(2):
  target=expected[a]*4*ga**2*L**4*p**2*(p-1)**2
  check('complete_original_indicial_row_'+str(a),sp.simplify(rows[a]-target)==0)
 check('complete_mixed_row_retained',sp.simplify(rows[1]-rows[0]/3)==0)
 # Explicit power and measure audit; no numerical cutoff is a continuum proof.
 alpha=sp.Rational(17,2);alpha_source=alpha-4;sigma=sp.Rational(801,200)
 check('regular_disk_H7_field_margin',alpha-7+1>0)
 check('regular_disk_H3_source_margin',alpha_source-3+1>0)
 check('global_radial_graph_square_not_integrable',2*(alpha_source-6)+1==-2)
 check('three_conormal_derivative_factor',alpha_source*(alpha_source-2)*(alpha_source-4)==sp.Rational(45,8))
 h,gamma,Cgeom,B2=[dec(pars[k]) for k in ['h','gamma','C','B2']]
 Lb=1/(1+h);dd=(96*arb.pi()*arb('1e22')-1/(60*arb.pi()))*Cgeom/B2
 aa=arb(17)/2;sig=arb(801)/200
 Anu=6*dd*Lb**4*aa**2*(aa-2)**2;AT=Anu/3
 Kdiv=arb(45)/8*Anu*Lb**(arb(5)/2-sig)*gamma**(-arb(5)/2)
 check('nonzero_full_column_indicial_coefficient',dd>0 and Anu>0 and Kdiv>0)
 c=app['source_and_output_maps'];MRbeta=dec(c['physical_source_to_energy_MR_over_beta']);Krow=dec(c['Krow_global_row_recovery'])
 Knew=dec(latest['constants']['complete_unfiltered_Knew'])
 co=Knew*MRbeta*Krow;cm=Knew*MRbeta/arb(5).sqrt();ci=Knew*MRbeta*400
 half=arb('1.667e-9')/2
 check('original_interior_weight_factor',sp.Rational(1,20)**(-2)==400)
 check('outer_residual_factor_unchanged',co.overlaps(dec(latest['conditional_residual_targets']['fast_residual_bolt_coefficient'])))
 check('actual_finitepart_factor_unchanged',cm.overlaps(dec(latest['conditional_residual_targets']['finitepart_bolt_coefficient'])))
 check('prior_complete_transport_slice_still_passes',dec(latest['complete_retained_inverse_transport_test']['whole_bolt_debit'])<half)
 return {'precision_bits':512,'complete_symbolic_indicial_rows':{'allocation_rows':pieces,
  'full_rows':[str(sp.factor(x)) for x in rows],'expected_d':str(d)},
  'actual_domain_growth_coefficients':{k:enc(v) for k,v in {'d':dd,'A_nu':Anu,'A_T':AT,
   'Kdiv_unit_amplitude':Kdiv,'liminf_epsilon_sqrt_Rfast_per_amplitude':Kdiv/arb(2).sqrt()}.items()},
  'complete_mixed_chart_residual_test':{k:enc(v) for k,v in {
   'outer_H3_pair_to_bolt':co,'actual_finitepart_to_bolt':cm,'whole_interior_L2_pair_to_bolt':ci,
   'outer_cap_other_debits_zero':half/co,'finitepart_cap_other_debits_zero':half/cm,
   'interior_cap_other_debits_zero':half/ci,'energy_dual_cap':half/Knew}.items()},
  'scope':{'global_radial_H3_comparison_from_regular_disk_H3_unbounded':True,
   'counterexample_uses_actual_full_coupled_operator':True,'both_original_rows_retained':True,
   'physical_operator_undefined':False,'original_IY_itself_proved_infinite':False,
   'finite_complete_mixed_chart_source_to_bolt_map':True,'accepted_global_sufficient_test_invalidated':False,
   'actual_retained_quantum_center_evaluated':False,'actual_trial_residual_evaluated':False,
   'original_defect_certified':False,'graph_physical_output_maps_certified':False,
   'actual_mass_target_certified':False,'every_coupled_entry_certified':False,
   'complete_coupled_estimate_certified':False,'higher_degree_adopted':False,
   'normalization_changed':False,'construction_obstruction':False,'new_physical_gate_closed':False}}
if __name__=='__main__':
 mode='acceptance' if '--acceptance' in sys.argv else 'development'
 out=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
 hashes={p:sh(ROOT/p) for p in INPUTS};run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,
  'input_sha256_before':hashes,'frozen_sources_utf8':{p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
  'runtime':{'python':sys.version,'python_sha256':sh(Path(sys.executable)),'flint_version':flint.__version__,'sympy_version':sp.__version__}}
 try:
  run['result']=execute();run['input_sha256_after']={p:sh(ROOT/p) for p in INPUTS}
  check('frozen_inputs_unchanged',hashes==run['input_sha256_after']);run['status']='passed'
 except Exception:run['status']='failed';run['traceback']=traceback.format_exc()
 run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
 ledger=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};ledger['runs'].append(run)
 out.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
 summary=json.dumps({'status':run['status'],'mode':mode,'checks':len(CHECKS),'traceback':run.get('traceback')},indent=2)+'\n'
 with (ROOT/'evidence'/(STEM+'.log')).open('ab') as stream:stream.write(summary.encode())
 print(summary,flush=True);sys.exit(0 if run['status']=='passed' else 1)
