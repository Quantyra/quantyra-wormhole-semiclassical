"""Complete quartic energy-dual application test; not the full F_n residual."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import factorial
import json, sys, traceback

ACCEPTANCE='--acceptance' in sys.argv
if ACCEPTANCE:sys.argv.remove('--acceptance')

ROOT=Path(__file__).resolve().parents[1]
STEM='current-quartic-resolvent-dual-test'
runtime=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
sys.path.insert(0,str(Path(runtime['packages']['flint']['root']).parent))
sys.path.insert(0,str(ROOT/'checks'))
from flint import arb,ctx
import current_signed_resolvent_projection as old
import current_quartic_column_energy_dual as heat

INPUTS=['notes/'+STEM+'.md','checks/current_quartic_resolvent_dual_test.py',
 'notes/current-signed-resolvent-projection.md','checks/current_signed_resolvent_projection.py',
 'checks/current_signed_resolvent_certificate.py','evidence/current-signed-resolvent-projection.json',
 'notes/current-geometric-mass-conversion.md','notes/current-quartic-column-energy-dual.md',
 'checks/current_quartic_column_energy_dual.py','evidence/current-quartic-column-energy-dual-assembly.json',
 'notes/current-short-step-heat-composition.md','evidence/current-short-step-heat-composition.json',
 'evidence/current-coupled-mass-adjoint-runtime.json']

def digest(p):return sha256(p.read_bytes()).hexdigest()
def dec(x):
 d=x['decimal_enclosure'];return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])
def enc(x):
 m,r,e=x.mid_rad_10exp(60)
 return {'ball':x.str(35),'decimal_enclosure':{'mid':str(m),'rad':str(r),'exponent':str(e)}}

def extend_cover(p,N,joint,previous):
 """Rebound complete coefficients on the added genuine finite-geometry interval."""
 T,P,V,op=joint
 lo,hi=arb(2)**-87,arb(2)**-56
 geo=old.geometry_numerators(p,N+1)
 AJ,VJ,F=old.rational_geometry_jets(geo,lo,hi)
 assert F>arb('.99')
 def bound(poly,weight):
  total=arb(0)
  for factor,bs,neg in old.t_control_polynomials(poly):
   # Old positive polynomial controls cover [0,2^-12], including this interval.
   pr=old.positive_y_tree(bs)[-1]
   power=lambda k:arb(max((lo**k).upper(),(hi**k).upper()))
   value=pr*power(factor)
   for k,c in neg.items():value+=abs(c)*power(k)
   total=arb(max(total,value.upper()))
  return total/arb(F.lower())**weight
 PB=[arb(0)]+[arb(max(dec(previous['P'][j]).upper(),bound(P[j],j).upper())) for j in range(1,N+1)]
 VB=[arb(0)]+[arb(max(dec(previous['V'][j]).upper(),bound(V[j],j).upper())) for j in range(1,N+1)]
 RB=arb(max(dec(previous['D_P_N']).upper(),bound(op(P[N],N),N+1).upper()))
 RVB=arb(max(dec(previous['D_V_N']).upper(),bound(op(V[N],N),N+1).upper()))
 amax=arb(max(dec(previous['a_jets'][0]).upper(),abs(AJ[0]).upper()))
 vmax=arb(max(dec(previous['v_jets'][1]).upper(),abs(VJ[1]).upper()))
 assert amax<arb('1.01') and vmax<arb('.5') and F<arb('1.01')
 return PB,VB,RB,RVB,F,amax,vmax

def execute():
 ctx.prec=512
 assert digest(Path(sys.executable))==runtime['python_sha256']
 for pkg in runtime['packages'].values():
  for name,h in pkg['files'].items():assert digest(Path(pkg['root'])/name)==h
 prior=json.loads((ROOT/'evidence/current-signed-resolvent-projection.json').read_bytes())['runs'][-1]
 assert prior['status']=='passed'
 rows=prior['numerics'];par=next(x for x in rows if x['kind']=='retained_exact_product_parameter_correlations')
 p={k:dec(v) for k,v in par['correlated_parameters'].items()};p['q']=arb(2000)**4
 N=9
 prev=next(x for x in rows if x['kind']=='complete_joint_radial_coefficient_bounds' and x['N']==N)
 cert=next(x for x in rows if x['kind']=='whole_signed_projection_certificate' and x['N']==N)
 joint=old.joint_coefficients(p,N)
 PB,VB,R,RV,Fext,Aactual,Vactual=extend_cover(p,N,joint,prev)
 M1=((p['q']+1)/p['C']).sqrt();A=arb('1.01');cc=A/4;cv=arb('1.51')/4
 J=sum((PB[j]/M1**j for j in range(1,N+1)),arb(0))/2
 JV=sum((VB[j]/M1**j for j in range(1,N+1)),arb(0))/2
 for i in range(1,N+1):
  for j in range(1,N+1):
   if i+j>=N:
    R+=Aactual*PB[i]*PB[j]/(2*M1**(i+j-N))
    RV+=((Vactual+Aactual)*PB[i]*PB[j]+Aactual*(VB[i]*PB[j]+PB[i]*VB[j]))/(2*M1**(i+j-N))
 assert J<arb('.1') and JV<arb('.1')
 eta=cc*max(1+J,2*J)/(2*M1);b0=2*(1-eta);dN=b0-(N+2)*A/M1
 assert dN>0
 Ca=(2+(cv*(1+J)+5*cc)/M1)/dN
 Kerr=((9*Ca+arb('9.72'))*R+6*RV)/dN
 fm=arb('.409');yc=arb(2)**-80
 # All accepted bulk cover cells have F>.4; extension has F>.99.
 assert all(dec(c['F'])>arb('.4') for c in prev['cells'])
 fmin=arb('.4')*yc
 assert p['beta']+p['u']/16<fm**2 and p['eps']<0 and p['tau']>0
 # Actual near-root deficit is far below the observation/start cuts.
 assert p['gamma']<arb(2)**-100 and p['tau']/((arb(2)**-87)**2)<arb('1e-12')
 masses=[(arb(n)*(n+p['q'])/p['C']).sqrt() for n in range(1,81)]
 ds=[2*(p['q']+2*n) for n in range(1,81)]
 sm=lambda k:sum((d/M**k for d,M in zip(ds,masses)),arb(0))
 je=lambda e:arb.pi().sqrt()*(e-arb('.5')).gamma()/(2*e.gamma())
 profile=lambda power,e:fm*je(e)*sm(power-1)+p['gamma']*sm(power)/2
 ec=json.loads((ROOT/'evidence/current-quartic-column-energy-dual-assembly.json').read_bytes())['result']['constants']
 CL=dec(ec['CLt']);alpha=dec(ec['sharper_exact_one_dimensional_coefficient_envelope'])
 norm=CL/(arb.pi()*p['B2']*p['beta']);pref=1/(p['B2']*p['beta'])
 Erad=norm*Kerr*(je(arb(N+2)/2)*sm(N-1)+p['gamma']*sm(N)/(2*fmin))

 # Whole rational majorant with the newly extended coefficient bounds.
 lim=6
 def mul(x,y):return [sum((x[j]*y[k-j] for j in range(k+1)),arb(0)) for k in range(lim+1)]
 plus=lambda x,y:[a+b for a,b in zip(x,y)]
 scale=lambda x,b:[a*b for a in x]
 one=[arb(1)]+[arb(0)]*lim
 qa=mul(PB[:7],PB[:7]);qva=scale(mul(PB[:7],VB[:7]),2)
 for j in [1,3,5]:qa[j]=qva[j]=arb(0)
 sa=[2*PB[j] if j%2==0 else arb(0) for j in range(7)]
 sva=[2*VB[j] if j%2==0 else arb(0) for j in range(7)]
 na=plus(scale(qa,4),scale(sa,2))
 nva=plus(plus(scale(qva,4),scale(qa,8)),plus(scale(sva,2),scale(sa,2)))
 q=scale(qa,arb(1)/4);qv=scale(qva,arb(1)/4);qq=mul(q,q);qqq=mul(qq,q)
 inv=plus(plus(plus(one,q),qq),qqq)
 inv2=plus(plus(plus(one,scale(q,2)),scale(qq,3)),scale(qqq,4))
 hm=plus(mul(nva,inv),mul(mul(na,qv),inv2))
 qval=sum((PB[i]*PB[j]/M1**(i+j) for i in range(1,N+1) for j in range(1,N+1) if (i+j)%2==0),arb(0))
 qvval=sum((2*PB[i]*VB[j]/M1**(i+j) for i in range(1,N+1) for j in range(1,N+1) if (i+j)%2==0),arb(0))
 sval=2*sum((PB[j]/M1**j for j in range(2,N+1,2)),arb(0))
 svval=2*sum((VB[j]/M1**j for j in range(2,N+1,2)),arb(0))
 qn=qval/4;assert qn<arb('.01')
 hval=(4*qvval+8*qval+2*svval+2*sval)/(1-qn)+(4*qval+2*sval)*qvval/4/(1-qn)**2
 Ktail=M1**8*(hval-sum((hm[j]/M1**j for j in [2,4,6]),arb(0)))
 assert Ktail>0
 Eapp=norm*Ktail*(arb(16)/315*sm(6)+p['gamma']/fmin*sm(7))
 Etherm=norm*p['gamma']/fmin*(sum((145*hm[2*k]*sm(2*k-1) for k in [1,2,3]),arb(0))+3*sm(-1))

 # Both true directed regular-branch returns.
 buffer=7*arb(2).log()/A;half=buffer/2
 xe=cc/(2*M1-cc);xb=xe+(1+xe)*(-(2*M1-cc)*half).exp()
 ve=(arb('.008')+cv*(1+arb('.008')**2)/(2*M1))/(1-cc*arb('.008')/M1)
 vb=ve+(arb('2.5')+ve)*(-(2*M1-2*cc*arb('.008'))*half).exp()
 assert xb<arb('.008') and vb<arb('.02')
 amin=b0*M1*buffer;assert amin>13
 Ca_old=(1+(cv*(1+J)+5*cc)/(2*M1))/(1-eta)
 ti=(3*(arb('2.5')+JV)+(3*Ca_old*amin+arb('3.24'))*(1+J))*(-amin).exp()
 Einner=4*norm/fmin*ti*profile(arb(-1),arb(6))

 # Convert ALL actual seam terms directly to L2(dz).
 D0=2*b0-A/M1;assert D0>0
 distance_rows=[(3*(arb('.008')+J)/4,4,0),(3*(arb('.008')+J)*A/20,5,0),
  (3*cc,5,0),(6*cc,4,1),(12*cc,3,2),(18*cc,2,3),(18*cc,1,4),
  (arb(33),4,0),(arb(48),3,1),(arb(72),2,2),(arb(72),1,3),(arb(36),0,4)]
 boundary=arb(0);brows=[]
 for row in cert['seam_terms']:
  j=row['j'];SB,SD=dec(row['P']),dec(row['V'])
  def mom(k,m=0,extra=0):
   power=arb(j+m+k-extra)-arb('.5')
   exponent=arb(j+k-extra)/2+arb('.75')
   return arb(factorial(2*k)).sqrt()/D0**(arb(k)+arb('.5'))*b0**extra*profile(power,exponent)
  term=6*SD*mom(0)
  val=3/(1-eta)*(mom(4)+mom(5,extra=1))
  val+=sum((coef*mom(k,m) for coef,k,m in distance_rows),arb(0))
  term+=2*SB*fm**4*val
  term*=norm/p['beta']**arb('.25');boundary+=term
  brows.append({'j':j,'energy_dual_bound':enc(term)})

 # Original canonical polynomial center and pointwise contact guards.
 beta,u,ep,tau=[p[k] for k in ['beta','u','eps','tau']]
 Q=[beta,arb(0),u,-4*beta-2*u-tau,3*beta+u+ep,-2*ep,ep]
 RR=heat.scale(heat.derivative(heat.derivative(Q)),-1);v=[arb(0)]*4+[arb(-1)]
 D=lambda x:heat.derivative(heat.mul(Q,heat.derivative(x)))
 HH=lambda x:heat.padd(D(x),heat.mul(RR,x))
 O6=heat.padd(heat.scale(heat.padd(HH(D(HH(v))),heat.mul(D(RR),HH(v)),D(heat.mul(v,D(RR)))),-arb(1)/105),
  heat.scale(heat.padd(heat.scale(heat.mul(heat.mul(heat.mul(RR,RR),RR),v),16),
  heat.scale(heat.mul(heat.mul(RR,RR),D(v)),24),heat.scale(D(heat.mul(heat.mul(RR,RR),v)),24),
  heat.scale(D(heat.mul(RR,D(v))),24)),-arb(1)/1008))
 P4=heat.padd(heat.scale(heat.mul(heat.mul(RR,RR),v),4),heat.scale(heat.mul(RR,D(v)),8),
  heat.scale(D(heat.mul(RR,v)),8),heat.scale(D(D(v)),8))
 osup=sum((abs(x) for x in O6),arb(0));p4sup=sum((abs(x) for x in P4),arb(0))
 center_coeff=heat.scale(O6,-pref*sm(4)/(4*arb.pi()))
 center_mid=[arb(x.mid()) for x in center_coeff]
 # A numerical polynomial center, with ALL coefficient/parameter radii charged.
 center_round=CL*sum((arb(x.rad())/arb(2*j+1).sqrt() for j,x in enumerate(center_coeff)),arb(0))
 center_norm=CL*heat.l2(center_mid)
 coldcenter=CL*pref*sm(4)*osup/(4*arb.pi())*yc.sqrt()
 hs=json.loads((ROOT/'evidence/current-short-step-heat-composition.json').read_bytes())['result']['constants']
 B0=dec(hs['weighted_short_step_smoothing_constant']);d0=dec(hs['weighted_short_step_error_smoothing_constant'])
 omega=dec(hs['complex_semigroup_growth']);delta=arb('.00025');t0=arb('.001')
 Bh=B0/delta.sqrt()+d0*delta**arb('3.5')
 H0=16*Bh**2/arb('.025');a=(p['q']+1)/p['C']
 assert a>324 and a>omega
 coldpos=CL*pref*yc.sqrt()*4*p['C']/t0*(H0*(omega*t0).exp()/(t0*(a-omega))+
  2/(arb.pi()*a*t0**2)+p4sup/(480*arb.pi()*a))
 Cb=dec(ec['boundary_Linfinity_collar_constant'])+dec(ec['boundary_Linfinity_exit_constant'])+arb('1e-100')
 colduv=CL*pref*yc.sqrt()*4*p['C']*(osup*t0/(4*arb.pi())+arb(2)/3*Cb*t0**arb('1.5')+arb(953)/2*t0**2)
 # n>80: common heat cap and complete positive-time angular tail.
 x=(-a*t0).exp();n=80
 ageom=(2*p['q']+4*(n+1))/(1-x)+4*x/(1-x)**2
 Mheat=dec(ec['heat_diagonal_bound_M']);P4norm=dec(ec['P4_L2'])
 CH,C0,C2=64*Mheat*CL/3,2*CL/(3*arb.pi()),CL*P4norm/(480*arb.pi())
 tailpos=pref*ageom*x**(n+1)/(a*(n+1))*(CH/t0+C0/t0**2+C2)
 tailuv=dec(ec['canonical_cap_at_001'])
 errors={'bulk_forced_resolvent':Erad,'entire_rational_tail':Eapp,
  'complete_thermal_difference':Etherm,'paired_seam':boundary,'true_bolt_data_return':Einner,
  'cold_positive_time':coldpos,'cold_common_UV':colduv,'cold_center_extension':coldcenter,
  'higher_angular_common_UV':tailuv,'higher_angular_positive_time':tailpos,
  'numerical_center_coefficient_and_parameter_error':center_round}
 total=sum(errors.values(),arb(0));zero=center_norm+total
 transport=dec(ec['retained_transport_column_error'])
 debit=alpha*(zero+transport);target=arb('1e-30')
 return {'N':N,'angular_cut':80,'observation_y':'2^-80','true_inner_data_y':'2^-87',
  'constants':{k:enc(x) for k,x in {'dN':dN,'Ca':Ca,'Kerr':Kerr,'Ktail':Ktail,'fmin':fmin,
   'R':R,'Rv':RV,'J':J,'Jv':JV,'extended_F':Fext,'Hcold0':H0,'O6_sup':osup,'P4_sup':p4sup,
   'CLt':CL,'alpha_one_dimensional':alpha,'signed_center_energy_dual_norm':center_norm,
   'complete_error_energy_dual':total,'zero_centered_massive_column_bound':zero,
   'one_column_residual_debit':debit,'one_column_sufficient_target':target,
   'one_column_slack':target-debit,'improvement_over_old_column_report':dec(ec['canonical_complete_massive_column_bound'])/zero}.items()},
  'error_contributions':{k:enc(x) for k,x in errors.items()},'seam_contributions':brows,
  'canonical_action_Hessian_center_coefficients':[enc(x) for x in center_mid],
  'canonical_action_Hessian_exact_coefficient_enclosures':[enc(x) for x in center_coeff],
  'energy_operator_center_sign':'negative of canonical action-Hessian center',
  'one_column_screen_passes':bool(debit<target),'actual_F_n':None,'full_residual_bound':None,
  'scope':{'complete_quartic_massive_column_test':True,'finite_thermal_mode_substitution':False,
   'complete_coupled_or_physical_admission':False,'rough_neighborhood_theorem':False,
   'existence_matching_validity_noise_stability':False}}

if __name__=='__main__':
 suffix='-acceptance' if ACCEPTANCE else ''
 out=ROOT/'evidence'/(STEM+suffix+'.json');assert not out.exists(),'Preserve immutable attempts.'
 before={p:digest(ROOT/p) for p in INPUTS}
 receipt={'started_utc':datetime.now(timezone.utc).isoformat(),'input_sha256_before':before,
  'ordinary_argument_utf8':(ROOT/INPUTS[0]).read_text(encoding='utf-8'),'source_utf8':Path(__file__).read_text(encoding='utf-8')}
 try:
  receipt['result']=execute();receipt['status']='complete_column_arithmetic_executed_requires_argument_audit'
 except BaseException:
  receipt['status']='failed';receipt['exception']=traceback.format_exc();raise
 finally:
  receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
  receipt['input_sha256_after']={p:digest(ROOT/p) for p in INPUTS}
  out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
 assert before==receipt['input_sha256_after']
 print(json.dumps({'status':receipt['status'],'constants':{k:v['ball'] for k,v in receipt['result']['constants'].items()},
  'errors':{k:v['ball'] for k,v in receipt['result']['error_contributions'].items()}},indent=2))
