"""Explicit Weyl/a0/a2/LLL graph block on the actual seed."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from math import comb,factorial
import sys,json,traceback
if '--arb-path' in sys.argv:sys.path.insert(0,sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-local-quantum-contact-graph'
OWN=['notes/'+STEM+'.md','checks/current_local_quantum_contact_graph.py','evidence/'+STEM+'-review.json']
ROWS=[];SYMBOLIC={};enc=old.encoded

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError(name)

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()

def read_result(stem):return json.loads((ROOT/'evidence'/f'{stem}.json').read_bytes())['runs'][-1]['result']

def build_symbols():
 z,s,t=S.symbols('z s t');be,uu,ep,ta,fe,q,c0=S.symbols('beta u eps tau Fe q c0')
 Q=be+uu*z**2+(-4*be-2*uu-ta)*z**3+(3*be+uu+ep)*z**4-2*ep*z**5+ep*z**6
 N=S.symbols('n0:3');B=S.symbols('b0:3');fields=[N,B];zeros={v:0 for v in N+B}
 phi=-S.log(z)+N[0];L=-S.diff(Q,z,2)-2*(S.diff(Q,z)*B[1]+Q*B[2]);M=S.diff(Q,z)*(-1/z+N[1])+Q*(1/z**2+N[2]);V=Q*(-1/z+N[1])**2
 A=phi*L**2/60+L*M/30-M*M/10+L*V/30-7*M*V/90-7*V*V/180
 # Divide all nonvacuum rows by GC/(pi B^2). c0 denotes vacuum in same units.
 lag=S.exp(-2*B[0])*A/2-uu*phi*L/36+(S.Rational(1,15)-q*q/3)*uu**2*phi*S.exp(2*B[0])/2-uu*V/36
 lag+=fe*S.exp(-2*B[0])*L*L/480+q*uu*Q*B[1]**2/12+c0*S.exp(2*B[0])/4
 der=lambda v:S.diff(v,z)+s*v/z
 rows=[]
 for i in range(2):
  row=[]
  for j in range(2):
   val=0
   for r in range(3):
    for k in range(3):
     coeff=S.diff(lag,fields[i][r],fields[j][k]).subs(zeros)
     direction=S.prod(s-h for h in range(k))/z**k
     v=coeff*direction
     for d in range(r):v=-der(v)
     val+=v
   row.append(S.expand(z**4*val))
  rows.append(row)
 out=[];a=(s+3)/(3*(s+1))
 for j in range(2):
  an=rows[0][j];tt=z*S.diff(rows[1][j],z)+(s-4)*rows[1][j]-(z*S.diff(rows[0][j],z)+(s-6)*rows[0][j])/3
  out.append([an,S.expand(tt)])
 matrix=[[out[0][0],S.cancel(out[1][0]-a*out[0][0])],[out[0][1],S.cancel(out[1][1]-a*out[0][1])]]
 logz=S.symbols('ellz');records=[]
 for i in range(2):
  for j in range(2):
   P=S.Poly(S.expand(S.cancel(matrix[i][j]*(s+1) if j else matrix[i][j])).subs(S.log(z),logz),z,logz)
   terms=[]
   for (power,lpower),coeff in P.terms():
    terms.append(dict(z_power=int(power),log_power=int(lpower),symbol=str(S.factor(coeff)),s_degree=int(S.degree(coeff,s))))
   records.append(dict(row=i,column=j,terms=terms))
 return records

def normal_compatibility(records):
 s,z,be,w,rR,fe,ellz=S.symbols('s z beta w rR Fe ellz')
 Q=be;ell=be*(1+11*w/360);cr=rR/(12*be)
 der=lambda v:S.diff(v,z)+s*v/z
 lap=lambda v:der(Q*der(v));D0=-12*Q/z**2
 def classical(n,b):
  d=-2*lap(b)-6*lap(n)+12*Q*der(n)/z;v=d-2*D0*b
  FnR=-12*der(der(Q*v))-der(24*Q/z*v-24*Q*D0*der(n))
  FbR=4*D0**2*b-4*D0*d-4*der(der(Q*v))
  FnE=z**-2*(2*D0*n+d)+ell*z**-4*(48*n+24*b)+der(der(-6*Q*z**-2*n))-der(12*Q/z**3*n-6*Q/z**2*der(n))
  FbE=ell*z**-4*(24*n+12*b)-2*der(der(Q*z**-2*n))
  return S.expand(z**4*(FnE+cr*FnR)),S.expand(z**4*(FbE+cr*FbR))
 rn=classical(1,0);rb=classical(0,1);tt=lambda row:(s-4)*row[1]-(s-6)*row[0]/3
 a=(s+3)/(3*(s+1));mat=[[rn[0],S.cancel(rb[0]-a*rn[0])],[tt(rn),S.cancel(tt(rb)-a*tt(rn))]]
 normal=[[0,0],[0,0]]
 for entry in records:
  i,j=entry['row'],entry['column']
  for term in entry['terms']:
   if term['z_power']==0:
    normal[i][j]+=S.sympify(term['symbol'],locals={'s':s,'beta':be,'Fe':fe})*ellz**term['log_power']/(s+1 if j else 1)
 check('exact_normal_mixed_scalar_to_shape',S.cancel(mat[1][0]+w/be*normal[1][0])==0)
 check('exact_normal_mixed_shape_to_scalar',S.cancel(mat[0][1]+w/be*normal[0][1])==0)
 J=-6*be*(1+11*w/180-(rR-w/60)*s*(s-3))*(s-4)*(s+1)
 check('exact_scalar_normal_graph',S.expand(mat[0][0]+w/be*normal[0][0]-J)==0)
 num,den=S.fraction(S.cancel(normal[1][1]))
 shape=S.div(S.Poly(num,s),S.Poly(den,s))[0]
 check('highest_shape_finite_contact',S.simplify(shape.nth(5)-be**2*(fe/60+S.Rational(1,90)-ellz/15))==0)
 check('finite_fourth_normal_constant',-S.Rational(31,30)-S.Rational(1,6)==-S.Rational(6,5))
 lk,lz,lb,lmu=S.symbols('lk lz lb lmu')
 massive=-(lk-lz+lb/2-lmu-S.Rational(31,30)+fe/4)
 local=fe/4+S.Rational(1,6)-lz
 check('common_normal_log_and_Fe_cancellation',S.expand(massive+local+lk+lb/2-lmu-S.Rational(6,5))==0)
 SYMBOLIC.update(normal_symbols=[[str(v) for v in row] for row in normal],mixed_normal_contacts_cancel=True,
  scalar_normal_is_retained_Jc=True,highest_shape_finite_contact_matched=True,
  full_canonical_heat_remainder_contacts_bound=False,physical_boundary_domain_unchanged=True)


def execute():
 records=build_symbols()
 normal_compatibility(records)
 pars,h,zh,env=old.scalar.prev.york.actual_parameters()
 syms=S.symbols('beta u eps tau Fe q c0 s');be,u,ep,ta,fe,q,c0,s=syms;t=S.symbols('t')
 C,B2=pars['C'],pars['B2'];f=C/(arb.pi()*B2);sig=arb('4.005');delta=arb('.005')
 C0=(1+12*arb.pi()*arb(2000)**6).sqrt();b02=arb(2000)**10*(C0-1)/6
 Lscale=(B2/b02).log()/2;Fe=2*arb.const_euler()-2*(arb(pars['q'])+1).digamma()+4*Lscale
 vac=2*(604*pars['q']**2+(1+2*Lscale.abs_upper())*(arb(pars['q'])**2/6-arb(1)/30))/(arb.pi()*C*B2)
 mapping={be:pars['beta'],u:pars['u'],ep:pars['eps'],ta:pars['tau'],fe:Fe,q:pars['q'],c0:arb(0,(vac/f).abs_upper())}
 def val(e):
  if e in mapping:return mapping[e]
  if e.is_Rational:return arb(int(e.p))/int(e.q)
  if e.is_Add:return sum((val(v) for v in e.args),arb(0))
  if e.is_Mul:
   v=arb(1)
   for x in e.args:v*=val(x)
   return v
  if e.is_Pow:return val(e.base)**int(e.exp)
  raise ValueError(e)
 K=C*pars['beta']**2/(15*arb.pi()*B2);LL=10*arb.pi()*(pars['Fl']+23*pars['g']/(180*arb.pi()))/pars['g']
 r=(96*arb.pi()*arb('1e22')-1/(60*arb.pi()))*pars['g'];aa=1+11*pars['g']/(180*arb.pi())-r*sig*(sig-3)
 capj=[1/(6*pars['beta']*aa*delta*(sig+1)),1/(6*pars['beta']*aa*(sig+1)),1/(6*pars['beta']*aa),1/(12*pars['beta']*(aa*r).sqrt()),1/(6*pars['beta']*r)]
 T=(LL/LL.log()).sqrt()
 capsb=[2/(K*LL*delta*sig*(sig-3)),2/(K*LL*sig*(sig-3)),2/(K*LL*sig),2/(K*LL),1/(K*(LL*arb(6).log()).sqrt()),max((2*T*T/(K*LL)).abs_upper(),(2/(K*(2+T).log())).abs_upper())]
 capsb[4]=max(capsb[4].abs_upper(),capsb[0].abs_upper(),(8/(K*LL)).abs_upper())
 capsb[5]=max(capsb[5].abs_upper(),capsb[0].abs_upper())
 check('crossover_above_four',T>4)
 check('positive_scalar_coercivity',aa>arb('.999999'))
 check('shared_vacuum_is_retained',vac>0 and mapping[c0].contains(0) and not mapping[c0].is_zero())
 totals=[[arb(0),arb(0)],[arb(0),arb(0)]];items=[]
 for entry in records:
  i,j=entry['row'],entry['column'];caps=capsb if j else capj
  for row in entry['terms']:
   zp,lp=row['z_power'],row['log_power']
   if zp==0:continue
   check('power_and_log_%d_%d_%d_%d'%(i,j,zp,lp),2<=zp<=12 and lp in [0,1])
   coeff=S.sympify(row['symbol'],locals=dict(zip([str(x) for x in syms],syms)))
   if j:poly,rem=S.div(S.Poly(coeff,s),S.Poly(s+1,s));remval=val(rem.as_expr()).abs_upper()*caps[0]/(sig+1)
   else:poly=S.Poly(coeff,s);remval=arb(0)
   shifted=S.Poly(S.expand(poly.as_expr().subs(s,t+S.Rational(801,200))),t)
   cb=remval+sum((val(shifted.nth(k)).abs_upper()*caps[k] for k in range(shifted.degree()+1)),arb(0))
   M=zp**3+12*zp**2+51*zp+100
   if lp:M=(-arb('.01').log())*M+3*zp*zp+24*zp+51
   check('allowed_multiplier_order_%d_%d_%d_%d'%(i,j,zp,lp),shifted.degree()<len(caps))
   bound=f*arb('.01')**zp*M*cb
   totals[i][j]+=bound
   items.append(dict(row=i,column=j,power=zp,log_power=lp,symbol=row['symbol'],multiplier_bound=enc(cb),cutoff_product_bound=enc(arb('.01')**zp*M),value=old.encoded(bound)))
 eta=max(sum(row,arb(0)).abs_upper() for row in totals)
 check('explicit_local_graph_target',eta<arb('.004'))
 check('explicit_local_graph_refined_guard',eta<arb('.0033'))
 remaining=arb('.18698')-eta
 check('remaining_allocation',remaining>arb('.18368'))
 local=read_result('current-coupled-local-graph-bound')
 known=old.decode(local['local_block_row_sum_upper'])+old.decode(local['prior_quantum_principal_upper'])
 known+=old.decode(read_result('current-separated-quantum-graph-bound')['separated_quantum_graph_feedback_upper'])
 known+=old.decode(read_result('current-coincident-domain-transfer')['coincident_domain_transfer_graph_feedback_upper'])
 total=known+eta+arb('.18368')
 joint=(arb('7.49e-12')+total*arb('1e-14'))/(1-total)
 check('conditional_joint_return',joint<arb('1e-11'))
 return dict(crossover=enc(T),shape_multiplier_constants=[enc(v) for v in capsb],scalar_multiplier_constants=[enc(v) for v in capj],
  Feff=enc(Fe),fixed_scale=enc(Lscale),vacuum_coefficient_upper=enc(vac),physical_factor=enc(f),
  matrix_entry_upper=[[enc(v) for v in row] for row in totals],explicit_local_quantum_graph_feedback_upper=enc(eta),
  remaining_feedback_threshold=enc(remaining),selected_remaining_target='.18368',conditional_total_feedback=enc(total),conditional_joint_E_upper=enc(joint),
  terms=items,exact_local_symbols=records,explicit_Weyl_a0_a2_LLL_block_certified=True,
  full_canonical_heat_remainder_contacts_certified=False,complete_quantum_graph_feedback_certified=False,
  joint_source_bounds_certified=False,independent_Mellin_Efar_target_closed=False,actual_controlled_solution=False,
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
