"""Complete last coupled filter on the implemented first-filter input."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction
from math import comb
import json,sys,argparse,traceback
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_filter_numerical_representation as num
old=num.old;lr=old.lr;fr=old.fr
arb,ctx,dec,enc,J=old.arb,old.ctx,old.dec,old.enc,old.J
acb=num.acb;pa,pm,ps,pv=old.pa,old.pm,old.ps,old.pv
up=old.up;minimum=old.minimum;cs=old.cs;ROOT=old.ROOT
STEM='current-coupled-last-filter-test'
INPUTS=list(dict.fromkeys(['notes/'+STEM+'.md','checks/current_coupled_last_filter_test.py',
 'evidence/current-filter-numerical-representation-attempt3.json',
 'evidence/current-filter-numerical-representation-attempt3-verification.json',
 'checks/current_filter_numerical_representation_review.py',*num.INPUTS]))
digest=lambda p:sha256(p.read_bytes()).hexdigest()
STAGES={}
setup=old.setup;point_tail=old.point_tail

def terms(p,c,name):
 return [(ps(pm(coeff,[c['Cv'],c['delta']]),c['kappa']),poles)
         for coeff,poles in old.terms(p,c,name)]


def components(p,c,name):
 c=setup(p,c);tt=terms(p,c,name)
 circles=[(arb(2),arb('.25'),arb('.9')),(arb(4)/3,arb('.2'),arb('.65'))]
 out=[];isolated=[]
 for center,radius,degree in circles:
  weight=arb(0)
  dh=abs(center-c['h'])
  gap=radius-dh if dh<radius else dh-radius
  dg=gap**2-abs(c['e']);assert gap>0 and dg>0
  # Bound complex square-root degree and differentiated Q integral.
  lower_degree=((1+4*(center-radius)/p['a']).sqrt()-1)/2
  kplus_upper=(1+(1+4*(center+radius)/p['a']).sqrt())/2
  assert lower_degree>degree and kplus_upper<2*(degree+1)
  for num,poles in tt:
   den=dg**2
   for mu in poles:
    dist=abs(mu-center)
    den*=radius-dist if dist<radius else dist-radius
    assert den>0
    if not any(abs(mu-C)<R for C,R,_ in circles):
     if not any(mu.overlaps(v) for v in isolated):isolated.append(mu)
   weight+=radius*sum((up(v)*(center+radius)**j for j,v in enumerate(num)),arb(0))/den
  out.append((degree,up(weight),arb(2)))
 for mu in isolated:
  residue=arb(0)
  for num,poles in tt:
   ix=[j for j,v in enumerate(poles) if mu.overlaps(v)]
   if not ix:continue
   assert len(ix)==1
   den=(mu**2-c['d1']*mu+c['d0'])**2
   for j,v in enumerate(poles):
    if j!=ix[0]:den*=v-mu
   assert not den.contains(0)
   residue+=pv(num,-mu)/den
  out.append((((1+4*mu/p['a']).sqrt()-1)/2,up(residue),arb(1)))
 return out

def scalar_spectral(p,c,count=512):
 ll=p['a']*p['ell']*(p['ell']+1);ln=c['large_mass_lambda']
 sums=[arb(0)]*4;cells=[]
 for i in range(count):
  t=arb(arb(2*i+1)/(2*count),arb(1)/(2*count));w=1-t
  Dt=t*t+c['d1']*t*w+c['d0']*w*w
  # Direct polynomial D_t may lose positivity in coarse intervals; its
  # Bernstein form supplies the rigorous lower bound, with no root division.
  guard=min(arb(1).lower(),(c['d1']/2).lower(),c['d0'].lower())
  assert guard>0
  Dlower=arb(guard)
  if Dt>Dlower:Dlower=Dt
  # Use the interval denominator only when wholly positive; otherwise a
  # universal Bernstein lower bound is safe (the upper is not needed).
  invD=1/Dt if Dt>0 else 1/Dlower
  coeff=(-2*c['scalar_c']*c['K0']/3)*(t+2*p['a']*w)
  coeff*=1+(c['Cv']-1)*ln*w/(ln*w+t)
  for mu in [c['muY'],c['muS'],ll]:coeff/=mu*w+t
  Fref=c['kappa']*(c['Cv']*w+c['delta']*t)*w*invD
  coeff*=invD*Fref
  costs=[up(t**j*w**(6-j)*coeff**2)/count for j in range(3)]
  xcoeff=(-2*c['K0']/3)*(t+2*p['a']*w)
  xcoeff*=1+(c['Cv']-1)*ln*w/(ln*w+t)
  for mu in [c['muY'],ll]:xcoeff/=mu*w+t
  xcoeff*=invD*Fref
  costs.append(up(w**4*xcoeff**2)/count)
  sums=[x+y for x,y in zip(sums,costs)]
  cells.append([enc(v) for v in costs])
 norms=[up(p['At']*(p['gamma']*v/2).sqrt()) for v in sums]
 return norms,cells

def reconstruct(p,c,bs,Fmax):
 L,a,ga=[p[k] for k in ['L','a','gamma']];sc=arb('.0004')
 old=cs('current-source-adapted-signed-response-enclosure-attempt3')
 scale=c['original_model_fractional_scale'];Dmax=cs('current-source-adapted-local-return-attempt4')['Dmax']
 prof={n:components(p,c,n) for n in ['f','r','diff','rdiff']}
 edges=[arb(0),ga]
 while edges[-1]<sc:
  n=arb((2*edges[-1]).mid())
  if n>=sc:break
  edges.append(n)
 edges.append(sc)
 while edges[-1]+arb('.002')<L:edges.append(edges[-1]+arb('.002'))
 edges.append(L);queue=list(zip(edges,edges[1:]));cells=[]
 eu2=eg2=et2=ex2=ew2=arb(0)
 while queue:
  s0,s1=queue.pop(0);x0=s0/ga;x1=s1/ga
  x=J([arb((x0+x1)/2,(x1-x0)/2),arb(1)]+[arb(0)]*5)
  z,q2,R,q=fr.geometry_jets(p,bs,x)
  if not (R.c[0]>0 and R.c[0]<2):
   assert s1-s0>ga/1024;mid=(s0+s1)/2;queue[0:0]=[(s0,mid),(mid,s1)];continue
  dR=R-1
  # Use factored R-1 rather than subtracting interval values near the bolt.
  dR=sum((bs[j]*ga**(j-2)*x**(j-1)/(2+a*x) for j in range(3,7)),0)
  A=-(ga*x/L)*sum(((z/L)**j for j in range(5)),0)+(z/L)**5*dR
  B=z**4/L**5*(2*z*(2*q2.d()*R.d()+q2*R.d().d())/3+8*ga*q.d()/3)
  C=z**4/L**5*(5*ga*R-z*R.d());E=-2*z**4/(3*L**5)
  Dx=z*q.d().d().d()+3*ga*q.d().d();qr=up(q2.c[0]).sqrt()
  def bounds(fn,rn,sp):
   nf=minimum(lr.profile_bounds(p,c,prof[fn],s0,s1),sp)
   nr=minimum(lr.profile_bounds(p,c,prof[rn],s0,s1),sp)
   hr=minimum(lr.profile_bounds(p,c,prof[rn],s0,s1,'flux'),qr*sp)
   wr=up(Dmax*(s1-s0).sqrt()*lr.profile_bounds(p,c,prof[rn],s0,s1,'tail_integral'))
   return nf,nr,hr,wr
  nf,nr,hr,wr=bounds('f','r',scale*Fmax**2)
  vu=sum((up(co.c[0])*v for co,v in zip([A,B,C,E],[nf,nr,hr,wr])),arb(0))
  fl=[]
  for name in ['r','f']:
   zf=scale*Fmax**2
   if s0>0:zf=minimum(zf,lr.profile_bounds(p,c,prof[name],s0,s1,'flux')/(2*x0+a*x0*x0).sqrt())
   fl.append(zf)
  gc=[qr*up((B.d()+c['muY']*C+E*Dx).c[0]),up((B+q2*C.d()).c[0]),
      qr*up((A.d()-C).c[0]),up(A.c[0]),qr*up(E.d().c[0])]
  vg=sum((b*v for b,v in zip(gc,[nr,fl[0],nf,fl[1],wr])),arb(0))
  nf,nr,hr,wr=bounds('diff','rdiff',scale*Fmax)
  vx=sum((up(co.c[0])*v for co,v in zip([-2*A/3,-B-4*a*A/9,-2*C/3,-E],[nf,nr,hr,wr])),arb(0))
  Aw=(-(ga*x/L)*sum(((z/L)**j for j in range(3)),0)+(z/L)**3*dR)/L**2
  Bw=z**2/L**5*(2*z*(2*q2.d()*R.d()+q2*R.d().d())/3+8*ga*q.d()/3)
  Cw=z**2/L**5*(5*ga*R-z*R.d());Ew=-2*z**2/(3*L**5)
  vw=sum((up(co.c[0])*v for co,v in zip([-2*Aw/3,-Bw-4*a*Aw/9,-2*Cw/3,-Ew],[nf,nr,hr,wr])),arb(0))
  eu2+=vu**2;eg2+=vg**2;ex2+=vx**2;ew2+=vw**2
  outside=not(s1<sc)
  if outside:et2+=vu**2
  cells.append(dict(s0=enc(s0),s1=enc(s1),outside=outside,
    u_error=enc(vu),gradient_error=enc(vg),xi_diff_error=enc(vx),xi_weighted_error=enc(vw)))
  if len(cells)%200==0:print(json.dumps(dict(reconstruction_cells=len(cells))),flush=True)
 result=dict(U0=Fmax**2*old['model_L2_upper']+eu2.sqrt(),
  U1=arb(2).sqrt()*(Fmax**2*old['model_gradient_upper']+eg2.sqrt()),
  U2=Fmax**2*old['retained_actual_source_Hessian_upper'],
  Utail=lr.profile_bounds(p,c,prof['f'],sc,L)+et2.sqrt(),
  exi=ex2.sqrt(),exiw=ew2.sqrt(),u_error=eu2.sqrt(),grad_error=eg2.sqrt())
 return {k:up(v) for k,v in result.items()},cells

def local_return(p,c,comparison,Fmax):
 geometry=lr.result('current-coupled-first-filter-test-attempt2')['local_return']
 const={k:dec(v) for k,v in geometry['constants'].items()}
 CU,CV=const['CU'],const['CV'];L,a,ga=[p[k] for k in ['L','a','gamma']]
 Dmax=cs('current-source-adapted-local-return-attempt4')['Dmax']
 profiles={n:components(p,c,n) for n in ['r','f','g','k','fn','fu','rf']}
 scale=c['original_model_fractional_scale']*Fmax**2
 def basis(s0,s1):
  ns={n:minimum(lr.profile_bounds(p,c,profiles[n],s0,s1),scale) for n in ['r','f','g','k']}
  qmax=2*s1/ga+a*(s1/ga)**2
  hs={n:minimum(lr.profile_bounds(p,c,profiles[n],s0,s1,'flux'),qmax.sqrt()*scale) for n in ['r','f','g']}
  W=Dmax*(s1-s0).sqrt()*lr.profile_bounds(p,c,profiles['r'],s0,s1,'tail_integral')
  return [ns['r'],hs['r'],ns['f'],hs['f'],ns['g'],hs['g'],ns['k'],W]
 sums=[arb(0)]*4;cells=[]
 for cell,rows in zip(geometry['cells'],comparison['whole_cell_coefficients']):
  s0,s1=dec(cell['s0']),dec(cell['s1']);bb=basis(s0,s1)
  costs=[sum((dec(co)*v for co,v in zip(row,bb)),arb(0)) for row in rows]
  sums=[s+v*v for s,v in zip(sums,costs)]
  cells.append(dict(s0=enc(s0),s1=enc(s1),costs=[enc(v) for v in costs]))
 bb=basis(L-arb('.01'),L);outer=[]
 for j,row in enumerate(comparison['outer_coefficients_after_factoring_z_squared']):
  val=sum((dec(co)*b for co,b in zip(row,bb)),arb(0));sums[j]+=val*val;outer.append(enc(val))
 mm={n:lr.moments(p,c,profiles[n]) for n in ['fn','fu','rf']}
 def wf(n):
  N,M,H=mm[n];return (2*ga).sqrt()*(N*M).sqrt()+a.sqrt()*M
 m=arb(min(arb(1).lower(),c['muT'].lower(),c['muS'].lower(),(a*p['ell']*(p['ell']+1)).lower()))
 poly=pa([c['fu_r']],pm([c['fu_f'],c['fu_g']],[c['muY'],arb(1)]))
 ratio=sum((up(v)/m**(2-i) for i,v in enumerate(poly)),arb(0))
 Crf=c['K0']*c['Cv']*ratio
 Grf=Fmax**2*p['At']*ga.sqrt()*Crf/(arb(40).sqrt()*m*m)
 Cxi=arb('1.01');cy=1+2/arb('4.5').sqrt();CN=lr.maximum(arb(8),arb(4000)/3+8*a+2000*L)
 Hrf=mm['rf'][2]
 Eadj=2*Cxi*Hrf+cy*Cxi*(arb(2).sqrt()*(wf('fn')+2*(2*a*wf('rf')+wf('fu')))
     +2*arb('.17').sqrt()*Hrf+CN/(2*arb('.1').sqrt())*(2*ga*Grf+2*Hrf/a.sqrt()))
 norms=[v.sqrt() for v in sums]
 classical=CU*norms[0]+CV*norms[1];quantum=CU*norms[2]+CV*norms[3]
 return classical+quantum+Eadj,dict(whole_cells=cells,outer_rows=outer,
  constants={k:enc(v) for k,v in dict(CU=CU,CV=CV,complete_classical_two_row_difference=classical,
   complete_local_quantum_rows=quantum,complete_adjoint_York_return=Eadj,
   complete_local_return=classical+quantum+Eadj).items()},
  profile_moments={n:[enc(v) for v in mm[n]] for n in mm})


def compute():
 ctx.prec=512;p,c,bs=fr.load();lr.model_rows(p,c);c=setup(p,c)
 L,a,ga,beta,d=[p[k] for k in ['L','a','gamma','beta','d']]
 assert abs(c['h']-arb(4)/3)<arb('.01') and abs(c['e'])<arb('.001')
 Fmax=arb('.251')
 Fguards=[Fmax*c['d0']-c['kappa']*c['Cv'],Fmax*c['d1']-c['kappa']*c['delta'],Fmax]
 assert all(v>0 for v in Fguards) and c['kappa']>0 and c['delta']>0
 profiles={n:components(p,c,n) for n in ['f','diff','xidiff','sig','Hsig']}
 (Ns,Gs,Ds,Nx),spectral_cells=scalar_spectral(p,c)
 rec,cells=reconstruct(p,c,bs,Fmax)
 Cenergy,comparison=num.reuse_output_comparison()
 Elocal,local=local_return(p,c,comparison,Fmax)
 S=arb('1e-12');zm=L-2*S;M3=arb(1000);rr=M3*2*S/(6*a)
 Q2S=2*ga*S+a*S*S;Q22S=4*ga*S+4*a*S*S
 Qmax=(1+rr)*Q22S;Qsmax=2*ga+4*a*S+M3*(2*S)**2/2
 CE=2*L+L*L*M3/(6*a);CB=5*L*L*M3/(6*a)+2*L*(1+rr)
 ms=lr.moments(p,c,profiles['sig']);md=lr.moments(p,c,profiles['Hsig']);mx=lr.moments(p,c,profiles['xidiff'])
 T,U,Tail=point_tail(p,profiles['sig'],S)
 TailD=point_tail(p,profiles['Hsig'],S)[2];TailX=point_tail(p,profiles['xidiff'],S)[2]
 cut=S.sqrt()*(4*(1+rr)*U/S+(15*Qmax/S**2+2*(Qsmax+2*Qmax/zm)/S)*T)
 Ny=Ns/zm**2;Ng=(1+rr).sqrt()/zm*(Gs+2*(Q22S/S).sqrt()*T)
 Nd=(1+rr)*Ds+(5*M3/(6*a)+2*(1+rr)/zm)*ms[2]+cut
 C,B2=p['C'],p['B2'];epsbar=C/(60*arb.pi()*B2);Cv=c['Cv']
 w0=arb('1.1e-27')*beta;taumax=arb('3.01')*beta
 KA=(2*d+11*epsbar/6)*w0+22*epsbar*taumax/3
 Es=((6*(d*Nd**2+(Cv+4*d*beta)*Ng**2+4*beta*Cv*Ny**2)+KA*Ng**2+2*w0*Ny**2)/beta).sqrt()
 S0=1/arb('23.99').sqrt();S1=(beta/arb('5.99')).sqrt();S2=(beta/(6*d)).sqrt()
 sc=cs('current-scalar-green-source-residual');Kb=sc['complete_mixed_remainder_constant']
 oldfirst=cs('current-source-adapted-signed-response-enclosure-attempt3')
 xiH=minimum(Nx,oldfirst['model_L2_upper'])/L**2+rec['exiw']
 cutP=S.sqrt()/zm**2*(L**2*Qmax*(15*T/S**2+4*U/(S*Q2S))+(L**2*Qsmax+2*L*Qmax)*2*T/S)
 rPparts=dict(actual_reconstruction=2*taumax*rec['exi']/(3*zm**2),
  entire_source_outside_scalar_cutoff=2*taumax*(TailX+rec['exi'])/3,
  tau_variation=2*mx[1]/zm**2,finite_geometry=(CE*md[1]+CB*ms[2])/zm**2,cutoff=cutP)
 rP=sum(rPparts.values(),arb(0))
 scalar_debits=dict(core_and_cutoff=6*(Cv*S0+d*S2)*rP/beta,
  scalar_remainder=(KA*S1*Ng+2*w0*S0*Ny)/beta,mixed_remainder=Kb*xiH/beta)
 Rs=sum(scalar_debits.values(),arb(0))
 ndiff=CE*md[1]+CB*ms[2]+cut+L**2*TailD
 source_debits=dict(mixed_remainder=arb('1.01')*Kb*Es/beta,
  cutoff_and_finite_geometry=4*arb('1.01')/beta*(sc['tau0']/L**4*(Cv*Tail+d*ndiff)
  +14*abs(p['epsilon_seed'])*(Cv*Ny+d*Nd)))
 Esource=sum(source_debits.values(),arb(0))
 oc=cs('current-source-adapted-complete-response-test');frac=cs('current-source-adapted-fractional-return')
 CV,CG,CD=[oc[k] for k in ['test_weighted_L2_constant','test_gradient_constant','test_Laplacian_constant']]
 U0,U1,U2,Utail=[rec[k] for k in ['U0','U1','U2','Utail']]
 factor=2/(B2*beta);qang=arb(2000)**4;pi=arb.pi();K=arb('1e10');sp=arb('1e-24')
 br=(sp/(4*K)).sqrt();bg=K*(arb(1)/4+2*((arb(1)/400)/br).log())
 quantum=dict(local_gain_two=factor*2*pi*C*(sp*U2**2/3+bg*U1**2).sqrt()*(sp*(3*CD**2/8+CG**2/4)/3+bg*CG**2).sqrt(),
  full_local_annular_kernel=factor*2*pi*arb('4e8')*C*(arb(400)/128).log()*U1*CG,
  finite_second_contact=factor*arb('1e13')*C*U1*CG,
  global_row_subtracted_gluing=factor*oc['complete_global_gluing_coefficient']*C*U0*CV,
  actual_constant_contact_cold_cap=factor*arb('2.1e8')*qang**2/C*U0*CV,
  actual_constant_contact_OUTSIDE_cold=factor*arb('1e18')*qang**2/C*(1+1/arb(5).sqrt())*Utail*CV,
  actual_constant_contact_gluing=factor*oc['constant_contact_gluing_coefficient']*C*U0*CV,
  known_curvature_polynomial=factor*100*C*(U1*CG/(36*pi)+816*U0*CV/(960*pi)))
 epsilon=cs('current-complete-polynomial-weak-transfer')['epsilon']
 ef=(1-epsilon)*arb('.99')*(1-arb('1.179e-8'))
 vac=arb('1e-21')/ef;ecoef=epsilon/(1-epsilon);ce=vac+ecoef;assert ce<1
 # Full f1 geodesic debit is named in its authoritative receipt.
 Egeo=Fmax**2*c['actual_reference_fractional_output_coefficient']*frac['actual_fractional_flow_norm_upper']
 Emodel=Fmax**2*c['complete_model_density_remainder']
 tt=dict(complete_source_comparison=Esource,classical_symbol_remainder=Emodel,
  complete_local_return=Elocal,complete_geodesic_remainder=Egeo)
 Rt_other=sum(tt.values(),arb(0))+sum(quantum.values(),arb(0))
 filtered=cs('current-filtered-mixed-green-residual');qt=arb('.338')/(1-filtered['epsilon']);ft=qt/(1+qt)
 Yenergy=(qt.sqrt()*Es+Rt_other)/(1-ce)
 quantum['shared_vacuum_uncertainty']=vac*Yenergy;quantum['retained_energy_comparison']=ecoef*Yenergy
 tt['complete_actual_minus_reference']=sum(quantum.values(),arb(0))
 Rt=sum(tt.values(),arb(0));Rfilter=Rt+qt.sqrt()*Rs
 Kfirst=c['first_residual_to_bolt_multiplier'];Klast=Kfirst/(ft*ft)
 screen=ft*ft*arb('1e33')
 previous=cs('current-filter-numerical-representation-attempt3')
 input_error=previous['complete_numerical_caller_Et_error'];inherited=ft*input_error
 vals=dict(kappa=c['kappa'],D_d1=c['d1'],D_d0=c['d0'],Fref_uniform_upper=Fmax,
  model_scalar_L2=Ns,model_scalar_gradient=Gs,model_scalar_H_norm=Ds,actual_cutoff_scalar_Es=Es,
  model_difference_xi_L2=Nx,actual_difference_xi_Hmu_upper=xiH,
  complete_scalar_residual=Rs,complete_TT_residual=Rt,complete_analytic_filter_error=Rfilter,
  actual_TT_trial_energy_upper=Yenergy,actual_filter_norm_upper=ft,
  last_filter_to_original_bolt_multiplier=Klast,last_filter_screen=screen,
  implemented_input_energy_error=input_error,inherited_input_filter_error=inherited,
  analytic_plus_implemented_input_error=Rfilter+inherited,
  analytic_report_over_screen=(Rfilter+inherited)/screen,
  remaining_UNEXECUTED_output_caller_allowance=screen-Rfilter-inherited,
  analytic_original_bolt_error=Klast*Rfilter,**rec)
 assert all(v.is_finite() for v in vals.values())
 data=dict(constants={k:enc(v) for k,v in vals.items()},
  scalar_residual_debits={k:enc(v) for k,v in scalar_debits.items()},
  scalar_second_order_debits={k:enc(v) for k,v in rPparts.items()},
  TT_source_comparison_debits={k:enc(v) for k,v in source_debits.items()},
  TT_residual_debits={k:enc(v) for k,v in tt.items()},quantum_debits={k:enc(v) for k,v in quantum.items()},
  spectral_cells=spectral_cells,reconstruction_cells=cells,local_return=local,
  filter_polynomial_positivity_guards=[enc(v) for v in Fguards],
  profiles={n:[dict(degree=enc(k),weight=enc(w),flux_factor=enc(ff)) for k,w,ff in v] for n,v in profiles.items()})
 STAGES['complete_analytic_last_filter']=data
 print(json.dumps(dict(analytic_error=str(Rfilter),input_debit=str(inherited),screen=str(screen),
  scalar_debits={k:str(v) for k,v in scalar_debits.items()},TT_debits={k:str(v) for k,v in tt.items()})),flush=True)
 return p,c,bs,Cenergy,comparison,vals,data


# Reuse exact rational and polynomial arithmetic, not another caller bound.
qball,qenc,middle,cmiddle,cball,cenc=num.qball,num.qenc,num.middle,num.cmiddle,num.cball,num.cenc
polyadd,polymul,polyeval,polypoles=num.polyadd,num.polymul,num.polyeval,num.polypoles


def profile_parameters(p,c):
 actual={k:c[k] for k in ['K0','kappa','Cv','delta','muT','muS','d1','d0']}
 actual.update(lambda_l=p['a']*p['ell']*(p['ell']+1),lambda_n=c['large_mass_lambda'],At=p['At'])
 nominal={k:middle(v) for k,v in actual.items()}
 def polys(v):
  N=polymul([v['Cv'],v['delta']],[v['Cv']*v['lambda_n'],1])
  N=[v['K0']*v['kappa']*x for x in N]
  D=polymul(polypoles([v[k] for k in ['muT','muS','lambda_l','lambda_n']]),[v['d0'],v['d1'],1])
  N=[v['kappa']*x for x in polymul(N,[v['Cv'],v['delta']])]
  D=polymul(D,[v['d0'],v['d1'],1])
  return N,D
 return actual,nominal,polys(actual),polys(nominal)

def parameter_error(p,actual,nominal,ad,nd):
 Na,Da=ad;Da=[arb(v) for v in Da];Nn,Dn=nd;amplitude=qball(nominal['At'])
 ratio=actual['At']/amplitude
 diff=polyadd([ratio*v for v in polymul(Na,[qball(v) for v in Dn])],
              [-v for v in polymul([qball(v) for v in Nn],Da)])
 assert len(diff)==12
 gard=lambda d:min((arb(v).lower()/comb(8,j) for j,v in enumerate(d)))
 Dnb=[qball(v) for v in Dn];ga=arb(gard(Da));gn=arb(gard(Dnb));assert ga>0 and gn>0
 total=arb(0);moments=[arb(0)]*3;cells=[];count=1024
 for i in range(count):
  tl=arb(i)/count;tu=arb(i+1)/count;wl=1-tu;wu=1-tl
  # Exact integration endpoints are rational; use outward positive guards.
  tl=arb(max(arb(0),tl.lower()));wl=arb(max(arb(0),wl.lower()))
  tu=up(tu);wu=up(wu)
  dl=sum((arb(v.lower())*tl**j*wl**(8-j) for j,v in enumerate(Da)),arb(0))
  nl=sum((arb(v.lower())*tl**j*wl**(8-j) for j,v in enumerate(Dnb)),arb(0))
  dl=arb(max(ga.lower(),dl.lower()));nl=arb(max(gn.lower(),nl.lower()))
  num=sum((up(v)*tu**j*wu**(11-j) for j,v in enumerate(diff)),arb(0))
  cost=up(wu**arb('3.5')*(num/(dl*nl))**2/count)
  for j in range(3):moments[j]+=up(tu**(2*j)*wu**(8-2*j)*(num/(dl*nl))**2/count)
  total+=cost;cells.append(enc(cost))
 X=up(amplitude*(p['gamma']*total/2).sqrt())
 return X,dict(cross_numerator=[enc(v) for v in diff],whole_spectral_cells=cells,
  denominator_Bernstein_guards=[enc(ga),enc(gn)],parameter_X_error=enc(X),
  model_H_power_norms=[enc(up(amplitude*(p['gamma']*v/2).sqrt())) for v in moments])

def direct_parameter_output(p,c,actual,nominal,pardata,comparison,X):
 delta=[dec(v) for v in pardata['cross_numerator']]
 nl={k:qball(v) for k,v in nominal.items()}
 circles=[(arb(2),arb('.25'),arb('.9')),(arb(4)/3,arb('.2'),arb('.65')),
          (arb(2)/3,arb('.1'),arb('.4')),(nl['lambda_n'],nl['lambda_n']/4,arb(6000))]
 profiles={}
 for name,power in [('f',0),('r',0),('g',1),('k',2)]:
  num=[arb(0)]*power+delta;parts=[]
  for center,radius,degree in circles:
   den=arb(1)
   for par in [actual,nl]:
    for key in ['muT','muS','lambda_l','lambda_n']:
     dist=abs(par[key]-center);gap=radius-dist if dist<radius else dist-radius
     assert gap>0;den*=gap
    h=par['d1']/2;e=par['d0']-h*h;dist=abs(center-h)
    gap=radius-dist if dist<radius else dist-radius
    dg=gap*gap-abs(e);assert gap>0 and dg>0;den*=dg**2
   if name=='r':
    dist=abs(c['muY']-center);gap=radius-dist if dist<radius else dist-radius
    assert gap>0;den*=gap
   kd=((1+4*(center-radius)/p['a']).sqrt()-1)/2
   kp=(1+(1+4*(center+radius)/p['a']).sqrt())/2
   assert kd>degree and kp<2*(degree+1)
   weight=radius*sum((up(v)*(center+radius)**j for j,v in enumerate(num)),arb(0))/den
   parts.append((degree,up(weight),arb(2)))
  profiles[name]=parts
 pp=dict(p);pp['At']=qball(nominal['At'])
 Nf,Ng,Nk=[old.minimum(dec(v),X) for v in pardata['model_H_power_norms']]
 bounds={'f':Nf,'r':up(Nf/c['muY']),'g':Ng,'k':Nk}
 Dmax=cs('current-source-adapted-local-return-attempt4')['Dmax']
 def basis(s0,s1):
  ns={name:old.minimum(lr.profile_bounds(pp,c,profiles[name],s0,s1),bounds[name]) for name in profiles}
  qmax=2*s1/p['gamma']+p['a']*(s1/p['gamma'])**2
  hs={name:old.minimum(lr.profile_bounds(pp,c,profiles[name],s0,s1,'flux'),qmax.sqrt()*X) for name in ['r','f','g']}
  W=Dmax*(s1-s0).sqrt()*lr.profile_bounds(pp,c,profiles['r'],s0,s1,'tail_integral')
  return [ns['r'],hs['r'],ns['f'],hs['f'],ns['g'],hs['g'],ns['k'],W]
 geometry=lr.result('current-coupled-first-filter-test-attempt2')['local_return']
 sums=[arb(0)]*4;cells=[]
 for cell,rr in zip(geometry['cells'],comparison['whole_cell_coefficients']):
  s0,s1=dec(cell['s0']),dec(cell['s1']);bb=basis(s0,s1)
  costs=[sum((dec(co)*b for co,b in zip(row,bb)),arb(0)) for row in rr]
  sums=[s+v*v for s,v in zip(sums,costs)]
  cells.append(dict(s0=enc(s0),s1=enc(s1),costs=[enc(v) for v in costs]))
 bb=basis(p['L']-arb('.01'),p['L']);outer=[]
 for j,row in enumerate(comparison['outer_coefficients_after_factoring_z_squared']):
  val=sum((dec(co)*b for co,b in zip(row,bb)),arb(0));sums[j]+=val*val;outer.append(enc(val))
 costs=[v.sqrt() for v in sums];model=[]
 for label in ['fn','fu']:
  model.append(p['L']**2*(up(c[label+'_r'])*bounds['r']+up(c[label+'_f'])*Nf+up(c[label+'_g'])*Ng))
 CU,CV=[dec(geometry['constants'][k]) for k in ['CU','CV']]
 local=CU*(costs[0]+costs[2]+model[0])+CV*(costs[1]+costs[3]+model[1])
 comp={k:dec(v) for k,v in comparison['constants'].items()}
 quantum=sum((dec(v) for v in comparison['quantum_entries'].values()),arb(0))*X
 geo=comp['complete_geodesic_dual']*X
 error=up((local+quantum+geo)/(1-comp['energy_feedback']))
 return error,dict(complete_local_dual=enc(local),complete_quantum_difference_except_feedback=enc(quantum),
  complete_geodesic_dual=enc(geo),complete_parameter_Et_error=enc(error),
  effective_Et_comparison_on_THIS_error=enc(error/X),classical_model_rows=[enc(v) for v in model],
  whole_cells=cells,outer_rows=outer,model_norms=[enc(v) for v in [Nf,Ng,Nk]],
  profiles={n:[dict(degree=enc(k),weight=enc(w),flux_factor=enc(ff)) for k,w,ff in v] for n,v in profiles.items()})

def contour_representation(p,nominal,nd,count=256):
 Nn,Dn=nd;B=[Fraction(1,4),Fraction(1,2),Fraction(1),Fraction(3)]
 Bpoly=polypoles(B);Gnum=polymul(Bpoly,Nn);ginf=Fraction(0)
 ln=nominal['lambda_n'];other=polymul(polypoles([nominal[k] for k in ['muT','muS','lambda_l']]),[nominal['d0'],nominal['d1'],1])
 other=polymul(other,[nominal['d0'],nominal['d1'],1])
 Rn=polyeval(Gnum,-ln)/polyeval(other,-ln)
 field=dict(amplitude=qenc(nominal['At']),base_poles=[qenc(v) for v in B],
  constant=qenc(ginf),large_pole=qenc(ln),large_pole_residue=qenc(Rn),upper_half_nodes=[])
 Gn=[qball(v) for v in Gnum];Dnball=[qball(v) for v in Dn]
 circles=[(Fraction(2),Fraction(1,4)),(Fraction(4,3),Fraction(1,5)),(Fraction(2,3),Fraction(1,10))]
 qcost=rcost=arb(0);details=[]
 h=qball(nominal['d1'])/2;e=qball(nominal['d0'])-h*h
 for center,radius in circles:
  C,R=qball(center),qball(radius)
  inside=[]
  for key in ['muT','muS','lambda_l','lambda_n']:
   dist=abs(qball(nominal[key])-C)
   assert dist<R/2 or dist>3*R/2
   if dist<R/2:inside.append(key)
  dh=abs(C-h)
  assert dh+abs(e).sqrt()<R/2 or dh-abs(e).sqrt()>3*R/2
  bounds=[]
  for rr in [R/2,3*R/2]:
   mindeg=arb(min(arb(1).lower(),(C-rr).lower()));assert mindeg>0
   den=arb(1)
   for key in ['muT','muS','lambda_l','lambda_n']:
    dist=abs(qball(nominal[key])-C);gap=rr-dist if dist<rr else dist-rr
    assert gap>0;den*=gap
   gap=rr-dh if dh<rr else dh-rr
   dg=gap*gap-abs(e);assert dg>0;den*=dg**2
   numerator=sum((up(v)*(C+rr)**j for j,v in enumerate(Gn)),arb(0))
   bounds.append(up(rr*numerator/den*256/(3*mindeg)))
  quad=bounds[1]/((arb(3)/2)**count-1)+bounds[0]/(arb(2)**count-1)
  qcost+=quad;round_circle=arb(0)
  for j in range(count//2):
   theta=2*arb.pi()*(arb(j)+arb('.5'))/count
   pole=acb(C)+R*acb(theta.cos(),theta.sin())
   weight=-(pole-C)*polyeval(Gn,-pole)/polyeval(Dnball,-pole)/count
   pn,wn=cmiddle(pole),cmiddle(weight);pb,wb=cball(pn),cball(wn)
   dmu=up(abs(pole-pb));dw=up(abs(weight-wb))
   mm=arb(min(arb(1).lower(),pole.real.lower(),pb.real.lower()));assert mm>0
   term=256*dw/(3*mm)+256*up(abs(wb))*dmu/(arb(13).sqrt()*mm**2)
   round_circle+=2*term
   field['upper_half_nodes'].append(dict(pole=cenc(pn),weight=cenc(wn),circle_center=qenc(center)))
  rcost+=round_circle
  details.append(dict(center=qenc(center),radius=qenc(radius),enclosed_linear_poles=inside,
   annular_norm_bounds=[enc(v) for v in bounds],quadrature_X_before_scale=enc(quad),rounding_X_before_scale=enc(round_circle)))
 scale=qball(nominal['At'])*p['gamma'].sqrt()
 return field,up(scale*qcost),up(scale*rcost),details

def execute():
 p,c,bs,CE,comparison,v,analytic=compute()
 scope=dict(actual_paired_domain=True,complete_two_equation_residual=True,
  full_common_regulator_angular_thermal_contacts_returns=True,numerical_last_filter=False,
  complete_last_filter_screen_passes=False,full_signed_response=False,full_bolt_certificate=False,
  original_P32_defect=False,actual_mass_column=False,original_graph_physical_boundary_admission=False,
  nonlinear_existence=False,compatibility=False,physical_matching=False,validity_noise=False,
  all_sector_stability=False,construction_obstruction=False,goal_complete=False)
 data=dict(analytic=analytic,constants={k:enc(x) for k,x in v.items()},numerical_output_error=None,
  numerical_field=None,full_signed_Gt=None,scope=scope)
 if not v['analytic_plus_implemented_input_error']<v['last_filter_screen']:
  return 'complete_last_filter_analytic_estimator_excessive_numerical_stage_not_executed',data
 ctx.prec=768
 # Reload all geometry at the higher evaluation precision without replacing
 # its physical interval, domain or normalization.
 p,c,bs=fr.load();lr.model_rows(p,c);c=setup(p,c)
 actual,nominal,ad,nd=profile_parameters(p,c)
 Xpar,pardata=parameter_error(p,actual,nominal,ad,nd)
 STAGES['parameter_error']=pardata
 Epar,direct=direct_parameter_output(p,c,actual,nominal,pardata,comparison,Xpar)
 STAGES['direct_parameter_output']=direct
 field,Xquad,Xround,contours=contour_representation(p,nominal,nd)
 STAGES['numerical_field']=field;STAGES['contour_bounds']=contours
 error=up(Epar+CE*(Xquad+Xround));full=v['analytic_plus_implemented_input_error']+error
 K=v['last_filter_to_original_bolt_multiplier']
 previous=cs('current-filter-numerical-representation-attempt3')
 known=previous['cumulative_known_error_subtotal_NOT_complete_output']+K*(v['complete_analytic_filter_error']+error)
 v.update(original_X_to_Et_comparison=CE,parameter_X_error=Xpar,parameter_Et_error=Epar,
  contour_X_error=Xquad,node_and_weight_X_error=Xround,numerical_output_error=error,
  generic_output_error=CE*(Xpar+Xquad+Xround),complete_implemented_last_filter_error=full,
  complete_report_over_screen=full/v['last_filter_screen'],remaining_last_filter_account=v['last_filter_screen']-full,
  complete_last_filter_bolt_error=K*full,incremental_global_bolt_error_charged=K*(v['complete_analytic_filter_error']+error),
  known_cumulative_error_subtotal_NOT_complete_bolt=known)
 print(json.dumps(dict(complete_last_filter_error=str(full),screen=str(v['last_filter_screen']),
  output_error=str(error),known_incomplete_subtotal=str(known))),flush=True)
 model=num.PairedFilterProfile(p,bs,field);queries=[]
 for label,s,bolt in [('bolt',arb(0),True),('thermal_scale',p['gamma'],False),
                      ('scalar_cutoff_scale',arb('1e-12'),False),('cold_split',arb('.0004'),False)]:
  values=model.at(s,bolt=bolt);assert all(x.is_finite() for x in values.values())
  queries.append(dict(label=label,s=enc(s),values={k:enc(x) for k,x in values.items()}))
  print(json.dumps(dict(query=label,xi=str(values['xi']),u=str(values['u']))),flush=True)
 queries.append(dict(label='original_outer_boundary',values={k:enc(x) for k,x in model.at_outer().items()}))
 passed=full<v['last_filter_screen'];scope.update(numerical_last_filter=True,complete_last_filter_screen_passes=bool(passed))
 data.update(constants={k:enc(x) for k,x in v.items()},numerical_output_error=enc(error),
  parameter_error=pardata,direct_parameter_output=direct,contour_bounds=contours,
  numerical_field=field,nominal_profile_parameters={k:qenc(x) for k,x in nominal.items()},
  point_and_endpoint_queries=queries,scope=scope)
 return ('complete_implemented_last_filter_passes_signed_output_open' if passed else
         'complete_implemented_last_filter_estimator_excessive'),data

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--receipt',default=STEM);args=ap.parse_args()
 out=ROOT/'evidence'/(args.receipt+'.json');assert not out.exists(),'Preserve immutable receipts'
 before={p:digest(ROOT/p) for p in INPUTS};protected=json.loads(lr.BASELINE.read_bytes())
 assert all(digest(Path(p))==h for p,h in protected.items())
 runtime=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes());pkg=runtime['packages']['flint']
 assert digest(Path(sys.executable))==runtime['python_sha256']
 assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
 receipt=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  ordinary_argument_utf8=(ROOT/'notes'/(STEM+'.md')).read_text(encoding='utf-8'),
  arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),python_sha256=runtime['python_sha256'])
 try:
  receipt['status'],receipt['result']=execute();after={p:digest(ROOT/p) for p in INPUTS};assert after==before
  assert all(digest(Path(p))==h for p,h in protected.items())
  receipt.update(input_sha256_after=after,protected_paths_unchanged=len(protected),pinned_flint_files_unchanged=len(pkg['files']))
 except Exception:receipt.update(status='failed',traceback=traceback.format_exc(),partial_results=STAGES)
 receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
 out.write_bytes((json.dumps(receipt,indent=2)+'\n').encode('utf-8'))
 print(json.dumps(dict(status=receipt['status'],error=receipt.get('traceback'),
  constants={k:x['ball'] for k,x in receipt.get('result',{}).get('constants',{}).items()}),indent=2),flush=True)
 if receipt['status']=='failed':raise SystemExit(1)
