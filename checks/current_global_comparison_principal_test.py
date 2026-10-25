"""Actual two-row defect test for a globally truncated comparison.

This is not a solution, a near quantum evaluation or a stability test.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,math,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
import sympy as S
import flint
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1]
STEM='current-global-comparison-principal-test'
OWN=['notes/'+STEM+'.md','checks/current_global_comparison_principal_test.py','evidence/'+STEM+'-review.json']
INPUTS=OWN+['notes/current-finite-local-action-assembly.md',
 'notes/current-coupled-local-graph-bound.md','notes/current-original-cutoff-coupled-estimator.md',
 'notes/current-short-travel-strong-domain.md','notes/current-strong-return-weight-repair.md',
 'notes/current-seed-full-hessian-realization.md','notes/current-global-mass-conversion.md',
 'evidence/current-global-mass-conversion.json','notes/current-complete-weak-mass-residual.md',
 'evidence/current-complete-weak-mass-residual.json','notes/current-actual-far-graph-energy-bound.md',
 'evidence/current-actual-far-graph-energy-bound.json','notes/current-coupled-schur-graph-attempt.md']
CHECKS=[]
def check(name,ok):
 assert name not in {v['name'] for v in CHECKS}
 CHECKS.append({'name':name,'passed':bool(ok)})
 if not ok:raise AssertionError(name)
def enc(v):
 m,r,e=v.mid_rad_10exp(50)
 return {'ball':str(v),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}
def dec(v):
 d=v['decimal_enclosure'];return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])
def prior(stem):
 r=json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
 assert r['status']=='passed'
 return r['result']
def execute():
 ctx.prec=512
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('ordinary_argument_frozen',sha256((ROOT/OWN[0]).read_bytes()).hexdigest()==review['scientific_sha256'][OWN[0]])
 par=prior('current-global-mass-conversion')['parameters']
 bv,uv,ev,tv,gv=[dec(par[k]) for k in ['beta','u','eps','tau','g']]
 z,s=S.symbols('z s');be,u,ep,ta,ell,cr,pp=S.symbols('beta u eps tau ell cr pp')
 Q=be+u*z**2+(-4*be-2*u-ta)*z**3+(3*be+u+ep)*z**4-2*ep*z**5+ep*z**6
 der=lambda v:S.diff(v,z)+s*v/z
 lap=lambda v:der(Q*der(v))
 D0=-S.diff(Q,z,2)+6*S.diff(Q,z)/z-12*Q/z**2
 d=-6*lap(S.Integer(1))+12*Q*der(S.Integer(1))/z
 FnR=-12*der(der(Q*d))-der((24*Q/z-12*S.diff(Q,z))*d-12*Q*(4*u+2*D0)*der(S.Integer(1)))
 FbR=-4*D0*d+4*der(S.diff(Q,z)*d)-4*der(der(Q*d))
 FnE=z**-2*(2*D0+d)+4*u*z**-2+48*ell*z**-4+der(der(-6*Q*z**-2))-der((12*Q/z-6*S.diff(Q,z))*z**-2-6*Q*z**-2*der(S.Integer(1)))
 FbE=4*u*z**-2+24*ell*z**-4-2*der(der(Q*z**-2))+2*der(S.diff(Q,z)*z**-2)
 AcN=S.Poly(S.expand(z**4*(FnE+cr*FnR)),z)
 AcB=S.Poly(S.expand(z**4*(FbE+cr*FbR)),z)
 dtilde=-S.diff(Q,z,2)*Q/15+7*S.diff(Q,z)*Q/(45*z)-7*Q**2/(15*z**2)-u*Q/9
 W=S.expand(z**4*pp*(-lap(lap(S.Integer(1)))/5-der(dtilde*der(S.Integer(1)))))
 actual=S.Poly(S.expand(AcN.as_expr()+W),s)
 mapping={be:bv,u:uv,ep:ev,ta:tv,ell:arb(2000)**-10*bv/gv,cr:8*arb.pi()*arb('1e22')*gv/bv,pp:gv/(2*arb.pi()*bv)}
 def val(e):
  if e in mapping:return mapping[e]
  if e.is_Rational:return arb(int(e.p))/int(e.q)
  if e.is_Add:return sum((val(v) for v in e.args),arb(0))
  if e.is_Mul:
   ans=arb(1)
   for v in e.args:ans*=val(v)
   return ans
  if e.is_Pow:return val(e.base)**int(e.exp)
  raise ValueError(e)
 r=(96*arb.pi()*arb('1e22')-1/(60*arb.pi()))*gv
 Cv=1+11*gv/(180*arb.pi())
 J=S.Poly(S.expand(-6*be*(S.Symbol('cv')-S.Symbol('rr')*s*(s-3))*(s-4)*(s+1)),s)
 mapping[S.Symbol('cv')]=Cv;mapping[S.Symbol('rr')]=r
 cn=S.Poly(S.expand(J.as_expr()+sum(AcN.nth(k)*z**k for k in (3,4))),s)
 cb=S.Poly(S.expand(S.cancel((s-6)*J.as_expr()/(3*(s-4)))+sum(AcB.nth(k)*z**k for k in (3,4))),s)
 check('both_original_rows_same_principal_zero',S.simplify(cb.nth(4)-cn.nth(4)/3)==0)
 check('complete_scalar_principal_identity',S.simplify(actual.nth(4)-(72*cr-pp/5)*Q**2)==0)
 def caps(row):
  return [[sum((abs(val(c))*arb('.63')**p[0]*p[0]**j for p,c in S.Poly(row.nth(k),z).terms()),arb(0)) for j in range(4)] for k in range(5)]
 cap_a,cap_n,cap_b=map(caps,(actual,cn,cb))
 t=S.Symbol('t');psi=(1-t*t)**10
 norms=[]
 for j in range(8):
  pol=S.Poly(S.diff(psi,t,j)**2,t)
  integ=sum(c*S.Rational(2,p[0]+1) for p,c in pol.terms() if p[0]%2==0)
  norms.append(val(integ).sqrt())
 norms=[v/norms[0] for v in norms]
 N=arb('1e8')
 dn=[N**(k-7)*sum((math.comb(k,j)*N**(-arb(j)/2)*norms[j] for j in range(k+1)),arb(0)) for k in range(8)]
 def upper(caps):
  L=[]
  for h in range(4):
   total=arb(0)
   for k in range(5):
    for j in range(h+1):
     guard=caps[k][j]
     if k==4 and j==0:guard=caps[4][1]/N.sqrt()
     total+=math.comb(h,j)*guard*dn[k+h-j]
   L.append(total)
  return sum((math.comb(3,h)*L[h]**2 for h in range(4)),arb(0)).sqrt()
 un,ub=upper(cap_n),upper(cap_b)
 actual_principal_lower=arb('7.36e-10')
 main_lower=actual_principal_lower*(2-dn[7])
 variation=cap_a[4][1]/N.sqrt()*dn[7]
 rest=sum((math.comb(3,j)*cap_a[k][j]*dn[k+3-j] for k in range(5) for j in range(4) if (k,j)!=(4,0)),arb(0))
 actual_lower=main_lower-variation-rest

 q=[arb(1),arb(0),uv/bv,-4-2*uv/bv-tv/bv,3+uv/bv+ev/bv,-2*ev/bv,ev/bv]
 da=1/(5760*arb.pi()**2*arb('1e22'))
 a3=2*q[3];a4=2*q[4]+q[2]**2
 v=lambda z:1-da+a3*z**3+a4*z**4
 qfull=lambda z:sum((q[k]*z**k for k in range(7)),arb(0))
 lo=arb('.61');hi=arb('.62')
 check('candidate_root_endpoint_signs',v(lo)>0 and v(hi)<0)
 interval=arb('.615','.005')
 check('candidate_root_unique_on_bracket',3*a3+4*a4*interval<0)
 for _ in range(150):
  mid=(lo+hi)/2
  if v(mid)>0:lo=mid
  elif v(mid)<0:hi=mid
  else:break
 check('outward_root_bracket_keeps_endpoint_signs',v(lo)>0 and v(hi)<0)
 zr=(lo+hi)/2+arb(0,((hi-lo)/2).abs_upper())
 check('candidate_root_final_interval',zr>arb('.61427243184514') and zr<arb('.61427243184516'))
 principal_scale=576*arb.pi()*arb('1e22')*gv*bv
 real_principal=principal_scale*(1-da)*qfull(zr)**2
 check('actual_principal_strictly_positive_at_artificial_zero',real_principal>arb('7.36e-10') and real_principal<arb('7.37e-10'))
 check('actual_normalized_principal_quarter', (1-da)*qfull(zr)**2>arb('.2499') and (1-da)*qfull(zr)**2<arb('.2501'))
 check('simple_zero_derivative',3*a3*zr**2+4*a4*zr**3>arb('-3.494') and 3*a3*zr**2+4*a4*zr**3<arb('-3.493'))
 check('collar_positive_but_global_negative',v(arb('.1'))>arb('.9925') and v(arb('.7'))<arb('-.3033'))
 check('packet_support_on_actual_interior',zr*(-1/N.sqrt()).exp()>arb('.6') and zr*(1/N.sqrt()).exp()<arb('.63'))
 # Independent strong-row check from the raw anomaly bilinear.
 M0=-S.diff(Q,z)/z+Q/z**2
 raw_d=-S.diff(Q,z,2)*Q/15-7*M0*Q/45-7*Q**2/(15*z**2)-u*Q/9
 vv=7*Q/(45*z)
 Wraw=z**4*pp*(-lap(lap(S.Integer(1)))/5+lap(vv*der(S.Integer(1)))-der(vv*lap(S.Integer(1)))-der(raw_d*der(S.Integer(1))))
 check('raw_and_reduced_anomaly_nu_column_agree',S.expand(Wraw-W)==0)
 principal_symbol=S.Poly(actual.nth(4),z)
 check('principal_matches_covariant_conformal_row',S.expand((72*cr-pp/5)*Q**2-actual.nth(4))==0)
 expected=72*cr*be**2-pp*be**2/5+72*cr*(2*be*(-4*be-2*u-ta)*z**3+(2*be*(3*be+u+ep)+u**2)*z**4)
 substitutions={S.Symbol('rr'):96*S.pi*S.Integer(10)**22*S.Symbol('gg')-S.Symbol('gg')/(60*S.pi),
                cr:8*S.pi*S.Integer(10)**22*S.Symbol('gg')/be,pp:S.Symbol('gg')/(2*S.pi*be)}
 check('candidate_principal_polynomial_identity',S.simplify((cn.nth(4)-expected).subs(substitutions))==0)
 check('all_actual_and_candidate_coefficients_polynomial',all(p[0]>=0 for row in (actual,cn,cb) for k in range(5) for p,c0 in S.Poly(row.nth(k),z).terms()))
 check('compact_packet_endpoint_jets',all(S.diff(psi,t,j).subs(t,1)==0 and S.diff(psi,t,j).subs(t,-1)==0 for j in range(10)))
 check('finite_derivative_norms',all(n>0 and n<arb('1e6') for n in norms))
 check('coefficient_jets_bounded_on_whole_interval',max(v for cc in (cap_a,cap_n,cap_b) for row in cc for v in row)<58)
 check('comparison_both_rows_upper',un+ub<arb('3.853e-12'))
 check('actual_scalar_row_lower',actual_lower>arb('7.27e-10'))
 ratio=(actual_lower-un-ub)/(un+ub)
 check('actual_defect_ratio_fails_candidate',ratio>187)
 # These are retained form statements, not a new strong inverse bound.
 eta=arb('.0022000000000000032')
 check('full_geometry_comparison_form_inverse_only',1/(1-eta)<arb('1.002205'))
 check('retained_actual_linear_mass',dec(prior('current-complete-weak-mass-residual')['full_coupled_linear_mass_upper'])<arb('1.627532e-15'))
 check('retained_actual_distant_graph',dec(prior('current-actual-far-graph-energy-bound')['actual_linear_far_graph_upper'])<arb('1.24e-20'))
 return {'precision_bits':ctx.prec,'root_enclosure':enc(zr),'root_derivative':enc(3*a3*zr**2+4*a4*zr**3),
   'actual_normalized_principal_at_root':enc((1-da)*qfull(zr)**2),'actual_principal_at_root':enc(real_principal),
   'candidate_at_original_cutoff':enc(v(arb('.1'))),'candidate_at_point_seven':enc(v(arb('.7'))),
   'packet_frequency':'100000000','packet_width_in_log_z':'0.0001',
   'comparison_scalar_H3_upper':enc(un),'comparison_shape_H3_upper':enc(ub),
   'comparison_joint_H3_upper':enc(un+ub),'actual_scalar_H3_lower':enc(actual_lower),
   'actual_candidate_defect_ratio_lower':enc(ratio),
   'psi_relative_derivative_L2':[enc(v) for v in norms],
   'actual_scalar_coefficient_jet_guards':[[enc(v) for v in row] for row in cap_a],
   'candidate_scalar_coefficient_jet_guards':[[enc(v) for v in row] for row in cap_n],
   'candidate_shape_coefficient_jet_guards':[[enc(v) for v in row] for row in cap_b],
   'symbolic_actual_scalar_row':str(actual.as_expr()),'symbolic_candidate_scalar_row':str(cn.as_expr()),
   'symbolic_candidate_shape_row':str(cb.as_expr()),
   'retained_full_geometry_form_inverse_upper':enc(1/(1-eta)),
   'scope':{'specified_global_comparison_defect_failure':True,'both_original_rows_in_probe':True,
     'complete_nonlocal_conformal_column_zero_by_action_identity':True,'full_anomaly_conformal_column_kept':True,
     'full_near_shape_quantum_estimate_executed':False,'physical_hessian_obstruction':False,
     'complete_coupled_return_certified':False,'actual_kappa_m_computed':False,'independent_phase_columns_bounded':False,
     'nonlinear_admission':False,'nonlinear_contraction':False,'exact_compatibility':False,'physical_matching':False,
     'validity_noise':False,'all_sector_stability':False,'actual_controlled_solution':False,'full_objective_complete':False}}

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    output=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    hashes={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'mode':mode,'input_sha256_before':hashes,
      'frozen_sources_utf8':{f:(ROOT/f).read_text(encoding='utf-8') for f in OWN},
      'runtime':{'python':sys.version,'executable_sha256':sha256(Path(sys.executable).read_bytes()).hexdigest(),
                 'flint_version':flint.__version__,'sympy_version':S.__version__}}
    try:
        run['result']=execute()
        run['input_sha256_after']={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in INPUTS}
        check('frozen_inputs_unchanged',hashes==run['input_sha256_after'])
        run['status']='passed'
    except Exception:
        run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    ledger['runs'].append(run);output.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    message=json.dumps({'mode':mode,'status':run['status'],'checks':len(CHECKS),
       'defect_ratio':run.get('result',{}).get('actual_candidate_defect_ratio_lower'),
       'candidate_root':run.get('result',{}).get('root_enclosure'),
       'traceback':run.get('traceback')},indent=2)+'\n'
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write(message.encode())
    print(message,flush=True)
    if run['status']!='passed':raise SystemExit(1)
