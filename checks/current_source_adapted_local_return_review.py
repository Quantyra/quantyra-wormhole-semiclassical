"""Author ordinary-argument, independent action-row and provenance review."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_source_adapted_local_return as m
arb,ctx,dec,enc,J=m.arb,m.ctx,m.dec,m.enc,m.J
ROOT=m.ROOT;STEM=m.STEM

def digest(p):return sha256(p.read_bytes()).hexdigest()

class B:
 """Two independent variations; coefficients are spatial Taylor jets."""
 def __init__(self,v=0,i=0,j=0,ij=0):self.v,self.i,self.j,self.ij=[J.cast(x,10) for x in [v,i,j,ij]]
 @staticmethod
 def lift(v):return v if isinstance(v,B) else B(v)
 def __add__(self,o):
  o=B.lift(o);return B(self.v+o.v,self.i+o.i,self.j+o.j,self.ij+o.ij)
 __radd__=__add__
 def __neg__(self):return B(-self.v,-self.i,-self.j,-self.ij)
 def __sub__(self,o):return self+-B.lift(o)
 def __rsub__(self,o):return B.lift(o)+-self
 def __mul__(self,o):
  o=B.lift(o);return B(self.v*o.v,self.i*o.v+self.v*o.i,self.j*o.v+self.v*o.j,
   self.ij*o.v+self.i*o.j+self.j*o.i+self.v*o.ij)
 __rmul__=__mul__
 def __truediv__(self,o):return self*B.lift(o).v.inv()
 def __pow__(self,n):
  out=B(1)
  for _ in range(n):out=out*self
  return out

def independent_rows(p,c,bs,x,initial):
 L,a,ga,beta,C,B2=[p[k] for k in ['L','a','gamma','beta','C','B2']]
 z=L-ga*x;q2=2*x+a*x*x
 rq=1+sum((bs[j]*ga**(j-2)*x**(j-1)/(2+a*x) for j in range(3,7)),0)
 q=q2*rq;D=z*q.d().d().d()+3*ga*q.d().d()
 # Local formal jets satisfy the defining ODEs, with independent initial data.
 coeff=[[arb(v)] for v in initial]
 for n in range(6):
  r,hr,f,hf,g,hg,k,W=[J(v) for v in coeff]
  rhs=[hr/q2,c['muY']*r-f,hf/q2,-g,hg/q2,-k,J([arb(0)]*(n+1)),D*r]
  for i in range(8):coeff[i].append(rhs[i].c[n]/(n+1))
 r,hr,f,hf,g,hg,k,W=[J(v) for v in coeff]
 MM=z*q.d().d()+4*ga*q.d()
 nu=z**4/(3*L**5)*(-z*rq*f+(z*rq*c['muY']-MM)*r-(7*ga*rq-z*rq.d())*hr+W)
 u=z**4/L**5*(z*rq*f+(2*MM/3-z*rq*c['muY'])*r+(5*ga*rq-z*rq.d())*hr-2*W/3)
 def H(v):return q2*v.d()
 def Delta(v):return q*v.d().d()+q.d()*v.d()
 ni,ui=nu,u;dni,dui=Delta(nu),Delta(u);hni=H(nu)
 ld=[]
 for j in range(5):
  nj,uj,dnj,duj,hnj=[int(i==j) for i in range(5)]
  phi=B(-z.log(),ni,nj)
  LL=B(-q.d().d(),-2*dui,-2*duj)
  MM=B(ga*q.d()/z+ga**2*q/z**2,dni,dnj)
  NN=B(ga**2*q/z**2,2*ga*rq/z*hni,2*ga*rq/z*hnj,2*rq/q2*hni*hnj)
  DD=LL-6*MM-6*NN
  def expjet(i,j,k):return B(1,k*i,k*j,k*k*i*j)
  en2,en4,eu2,eum2=expjet(ni,nj,2),expjet(ni,nj,4),expjet(ui,uj,2),expjet(ui,uj,-2)
  le=en2*DD/(2*z**2)+en2*eu2/(C*z**2)+3*arb(2000)**-10*B2/C*en4*eu2/z**4
  lm=-arb.pi()*arb(2000)**16/(C*B2)*eu2
  lr=8*arb.pi()*arb('1e22')/(C*B2)*(4*eu2+4*C*DD+C*C*eum2*DD**2)
  wb=phi*LL**2/60+LL*MM/30-MM**2/10+LL*NN/30-7*MM*NN/90-7*NN**2/180
  lw=(C*C*eum2*wb-C*phi*LL/18+(arb(1)/15-arb(2000)**8/3)*phi*eu2-C*NN/18)/(2*arb.pi()*C*B2)
  ld.append(((le+lm+lr).ij,lw.ij))
 out=[]
 for allocation in [0,1]:
  d=[row[allocation] for row in ld]
  fn=d[0]+Delta(d[2])-H(d[4])-q2.d()*d[4]
  fu=d[1]+Delta(d[3])
  fn=-fn/beta;fu=-fu/beta
  if allocation==0:
   fn=fn-c['fn_r']*r-c['fn_f']*f-c['fn_g']*g
   fu=fu-c['fu_r']*r-c['fu_f']*f-c['fu_g']*g
  else:
   R=-q.d().d();du=Delta(u)
   U0=2*(604*arb(2000)**8+(1+2*arb('1e-8'))*(arb(2000)**8/6-arb(1)/30))/(arb.pi()*C*B2)
   a2=arb(0,100)*C/(480*arb.pi()*B2)
   p4=4*R**2*u+8*R*du+8*Delta(R*u)+8*Delta(du)
   fu=fu-arb(0,m.up(U0))/beta*u-a2/beta*p4+arb(2000)**4/(6*arb.pi()*B2*beta)*du
  out.extend([fn.c[0]*z.c[0]**2,fu.c[0]*z.c[0]**2])
 return out

def main():
 ctx.prec=512
 path=ROOT/'evidence'/(STEM+'-verification.json');assert not path.exists()
 source=ROOT/'evidence'/(STEM+'-attempt4.json');d=json.loads(source.read_bytes())
 assert d['status']=='complete_local_return_enclosed_fractional_entry_unknown'
 assert d['input_sha256_before']==d['input_sha256_after']
 for p,h in d['input_sha256_before'].items():assert digest(ROOT/p)==h,p
 assert d['ordinary_argument_utf8']==(ROOT/'notes'/(STEM+'.md')).read_text(encoding='utf-8')
 assert d['arithmetic_source_utf8']==(ROOT/'checks/current_source_adapted_local_return.py').read_text(encoding='utf-8')
 assert digest(Path(sys.executable))==d['python_sha256']
 baseline=json.loads(m.BASELINE.read_bytes());assert all(digest(Path(p))==h for p,h in baseline.items())
 r=d['result'];v={k:dec(t) for k,t in r['constants'].items()}
 assert v['complete_local_return']<arb('4.394648e29')
 assert v['complete_classical_two_row_difference']<arb('4.394062e29')
 assert v['complete_local_quantum_rows']<arb('5.732109e25')
 assert v['complete_adjoint_York_return']<arb('1.246910e24')
 assert v['local_return_over_remaining_account']<arb('.001018931')
 assert v['account_after_local_NOT_physical_margin']>arb('4.308605e32')
 assert v['local_original_bolt_estimator_debit']<arb('1.406657e-14')
 assert r['cell_count']==152 and len(r['cells'])==152
 assert dec(r['cells'][0]['s0']).is_zero()
 for previous,current in zip(r['cells'],r['cells'][1:]):assert previous['s1']==current['s0']
 for key in ['actual_fractional_flow_norm','actual_application_error','complete_first_reference_residual']:
  assert r[key] is None
 p={k:dec(t) for k,t in m.result('current-global-green-trial-enclosure')['trial_parameters'].items()}
 c=m.constants('current-source-adapted-geodesic-first-attempt');m.model_rows(p,c)
 mp={k:dec(t) for k,t in m.result('current-global-mass-conversion')['parameters'].items()}
 L,a,ga,beta=[p[k] for k in ['L','a','gamma','beta']]
 u,e,t=[mp[k] for k in ['u','eps','tau']]
 Q=[beta,arb(0),u,-4*beta-2*u-t,3*beta+u+e,-2*e,e]
 bs=[sum((Q[j]*m.comb(j,i)*L**(j-i)*(-1)**i for j in range(i,7)),arb(0)) for i in range(7)]
 bs[:3]=[arb(0),2*ga,a]
 assert dec(r['cells'][-1]['s1']).overlaps(L-arb('.01'))
 assert (1+4*arb('1.75')/arb('1.001')).sqrt()>arb('2.8')
 fixtures=[]
 for xp in [arb('.3'),arb(3),arb('1e20')]:
  x=J([xp,arb(1)]+[arb(0)]*8);rr=m.rows(p,c,bs,x)
  for i in range(8):
   initial=[int(k==i) for k in range(8)]
   alt=independent_rows(p,c,bs,x,initial)
   for j in range(4):assert rr[j][i].overlaps(alt[j]),(str(xp),i,j,str(rr[j][i]),str(alt[j]))
   fixtures.append(dict(x=enc(xp),basis_index=i,all_four_rows_overlap=True))
  print('Independent action rows checked at x='+xp.str(8),flush=True)
 failed=[]
 for i in [1,2,3]:
  f=ROOT/'evidence'/(STEM+'-attempt'+str(i)+'.json');a0=json.loads(f.read_bytes())
  assert a0['status']=='failed'
  failed.append(dict(path=str(f.relative_to(ROOT)),sha256=digest(f),reason=a0['traceback'].splitlines()[-1]))
 reviews=dict(
  domain='The actual analytic profile is reconstructed by both exact finite-Q York rows. Its rho^4 log term has fourth strong derivatives in L2; local integration-by-parts fluxes vanish at the regular bolt and on the original outer mass/log phase. No H5 source assertion or new boundary wall.',
  row_comparison='The transformed Dirichlet equation has Pxi>=4.5/z^2 and QPxi>=1 by the stated polynomial guards. Testing it gives CU,CV in the original z^-4 L2 norm, not C5 physical admission.',
  adjoint='The original-endpoint trial Y=z^2 Q2 rf_z/2 gives identity (4) after one exact adjoint transfer; the h-prime terms cancel using both rows. All three residual terms and the global Volterra return are charged.',
  profile='Near-coincident masses are enclosed on a single radius-.25 contour; positive Q integrals dominate complex degrees by degree .9, with the retained flux factor two. All other poles, including the large-mass Cv-1 term, remain. Spectral norm caps are auxiliary-profile bounds only.',
  full_domain='All 152 whole real cells and all 64 complex outer coefficient arcs are included. Weighted coefficients are holomorphic through zero after extracting the sole linear log; original log-square integration and one-sided boundary fields are unchanged.',
  local_quantum='Complete Weyl, LLL, shared bounded a0 and actual a2 restorations retained. The exact nonlocal ultraviolet logarithm remains a separate UNKNOWN numerical flow-norm entry, not a discarded contact or a finite mode approximation.',
  independent_rows='96 coefficient comparisons from separate two-direction differentiation of the accepted action check the strong-row implementation, including the second gradient contact. These fixtures support the ordinary generator/adjoint argument; they are not a continuum proof.',
  arithmetic='Integer constants are converted to Arb before division; broad outer cells are subdivided without dropping regions. Three failed development receipts are retained; the fourth is the completed frozen enclosure.',
  decision='The complete local return fits 0.101894% of the remaining first-stage account. Stop optimizing it; next quantify the actual fractional geodesic norm together with application errors. No first inverse or original admission certificate.',
  scope='Fixed seed and phase only. No independent boundary-phase, nonlinear neighborhood, original P32 defect/mass/graph/physical output, existence, compatibility, matching, validity/noise or stability result. Author review, not independent or Lean certification.')
 output=dict(status='complete_local_return_reviewed_fractional_response_unresolved',
  finished_utc=datetime.now(timezone.utc).isoformat(),receipt_sha256=digest(source),receipt_path=str(source.relative_to(ROOT)),
  review_source_sha256=digest(Path(__file__)),review_source_utf8=Path(__file__).read_text(encoding='utf-8'),
  checked_inputs=len(d['input_sha256_before']),protected_paths_unchanged=len(baseline),
  failed_development_receipts=failed,independent_action_fixtures=fixtures,review=reviews,
  constants={k:enc(t) for k,t in v.items()},
  complete_first_reference_inverse=False,complete_original_admission=False,goal_complete=False)
 path.write_bytes((json.dumps(output,indent=2)+'\n').encode('utf-8'))
 print(json.dumps(dict(status=output['status'],checked_inputs=output['checked_inputs'],
  independent_row_comparisons=len(fixtures)*4,local_return=v['complete_local_return'].str(12)),indent=2))

if __name__=='__main__':main()
