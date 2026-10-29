"""Complete coupled first-filter residual on the accepted analytic input."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import argparse, json, sys, traceback
sys.path.insert(0, str(Path(__file__).resolve().parent))
import current_source_adapted_signed_response_enclosure as prior
import current_scalar_green_source_residual as scalar_prior
import current_source_adapted_complete_response_test as quantum_prior
fr=prior.fr;lr=prior.lr;arb,ctx,dec,enc,J=prior.arb,prior.ctx,prior.dec,prior.enc,prior.J
pa,pm,ps,pv=lr.pa,lr.pm,lr.ps,lr.pv
ROOT=prior.ROOT;STEM='current-coupled-first-filter-test'
INPUTS=list(dict.fromkeys(['notes/'+STEM+'.md','checks/current_coupled_first_filter_test.py',
 'evidence/current-source-adapted-signed-response-enclosure-attempt3.json',
 'evidence/current-source-adapted-signed-response-enclosure-verification.json',
 'notes/current-filtered-mixed-green-residual.md',
 'evidence/current-filtered-mixed-green-residual.json',
 'evidence/current-scalar-green-source-residual-verification.json',
 'evidence/current-coupled-first-filter-test.json',
 'evidence/current-coupled-first-filter-test-verification.json',
 *prior.INPUTS,*lr.INPUTS,*scalar_prior.INPUTS,*quantum_prior.INPUTS]))
digest=lambda p:sha256(p.read_bytes()).hexdigest()
up=prior.up;minimum=prior.minimum
cs=lr.constants

def setup(p,c):
 L,beta,d=p['L'],p['beta'],p['d']
 tau0=cs('current-scalar-green-source-residual')['tau0']
 c=dict(c);c['scalar_c']=2*tau0/(3*L**2);c['delta']=d*L**2
 c['kappa']=4*tau0*c['scalar_c']/(beta*L**4*c['classical_A'])
 c['d1']=c['muT']+c['muS']+c['kappa']*c['delta']
 c['d0']=c['muT']*c['muS']+c['kappa']*c['Cv']
 c['h']=c['d1']/2;c['e']=c['d0']-c['h']**2
 return c

def terms(p,c,name):
 a=p['a'];ll=a*p['ell']*(p['ell']+1);ln=c['large_mass_lambda']
 if name in ['diff','rdiff','xidiff','sig','Hsig']:
  poles=[ll];num=[c['K0']]
  if name=='rdiff':poles.append(c['muY'])
  if name not in ['diff','rdiff']:
   poles.append(c['muY']);num=ps(pm(num,[2*a,arb(1)]),-arb(2)/3)
  if name in ['sig','Hsig']:
   poles.append(c['muS']);num=ps(num,c['scalar_c'])
  if name=='Hsig':num=pm(num,[arb(0),arb(1)])
 else:
  poles=[c['muT'],c['muS'],ll]
  num=ps([c['Cv'],c['delta']],c['K0']*c['kappa'])
  if name in ['r','rf','fn','fu']:poles.append(c['muY'])
  if name in ['fn','fu','rf']:
   label='fu' if name=='rf' else name
   num=pm(num,pa([c[label+'_r']],pm([c[label+'_f'],c[label+'_g']],[c['muY'],arb(1)])))
  if name=='rf':poles.append(2*a)
  if name=='g':num=pm(num,[arb(0),arb(1)])
  if name=='k':num=pm(num,[arb(0),arb(0),arb(1)])
 return [(num,poles),(ps(num,(c['Cv']-1)*ln),poles+[ln])]

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
   den=dg
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
   den=mu**2-c['d1']*mu+c['d0']
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
  coeff*=invD
  costs=[up(t**j*w**(6-j)*coeff**2)/count for j in range(3)]
  xcoeff=(-2*c['K0']/3)*(t+2*p['a']*w)
  xcoeff*=1+(c['Cv']-1)*ln*w/(ln*w+t)
  for mu in [c['muY'],ll]:xcoeff/=mu*w+t
  xcoeff*=invD
  costs.append(up(w**4*xcoeff**2)/count)
  sums=[x+y for x,y in zip(sums,costs)]
  cells.append([enc(v) for v in costs])
 norms=[up(p['At']*(p['gamma']*v/2).sqrt()) for v in sums]
 return norms,cells

def point_tail(p,comps,S):
 X=p['a']*S/p['gamma'];assert X>1
 T=U=tail=arb(0)
 for k,w,ff in comps:
  T+=w*arb.pi()*X**(-k-1)
  U+=w*ff*arb.pi()*p['gamma']*(k+1)*X**(-k)
  tail+=w*arb.pi()*S.sqrt()*X**(-k-1)/(2*k+1).sqrt()
 return [up(p['At']*v) for v in [T,U,tail]]

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
  nf,nr,hr,wr=bounds('f','r',scale*Fmax)
  vu=sum((up(co.c[0])*v for co,v in zip([A,B,C,E],[nf,nr,hr,wr])),arb(0))
  fl=[]
  for name in ['r','f']:
   zf=scale*Fmax
   if s0>0:zf=minimum(zf,lr.profile_bounds(p,c,prof[name],s0,s1,'flux')/(2*x0+a*x0*x0).sqrt())
   fl.append(zf)
  gc=[qr*up((B.d()+c['muY']*C+E*Dx).c[0]),up((B+q2*C.d()).c[0]),
      qr*up((A.d()-C).c[0]),up(A.c[0]),qr*up(E.d().c[0])]
  vg=sum((b*v for b,v in zip(gc,[nr,fl[0],nf,fl[1],wr])),arb(0))
  nf,nr,hr,wr=bounds('diff','rdiff',scale)
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
 result=dict(U0=Fmax*old['model_L2_upper']+eu2.sqrt(),
  U1=arb(2).sqrt()*(Fmax*old['model_gradient_upper']+eg2.sqrt()),
  U2=Fmax*old['retained_actual_source_Hessian_upper'],
  Utail=lr.profile_bounds(p,c,prof['f'],sc,L)+et2.sqrt(),
  exi=ex2.sqrt(),exiw=ew2.sqrt(),u_error=eu2.sqrt(),grad_error=eg2.sqrt())
 return {k:up(v) for k,v in result.items()},cells

def compute():
 ctx.prec=512;p,c,bs=fr.load();lr.model_rows(p,c);c=setup(p,c)
 L,a,ga,beta,d=[p[k] for k in ['L','a','gamma','beta','d']]
 assert abs(c['h']-arb(4)/3)<arb('.01') and abs(c['e'])<arb('.001')
 Fmax=arb('.251')
 Fguards=[Fmax*c['d0']-c['kappa']*c['Cv'],Fmax*c['d1']-c['kappa']*c['delta'],Fmax]
 assert all(v>0 for v in Fguards) and c['kappa']>0 and c['delta']>0
 profiles={n:components(p,c,n) for n in ['f','diff','xidiff','sig','Hsig']}
 (Ns,Gs,Ds,Nx),spectral_cells=scalar_spectral(p,c)
 # The explicit trial and every non-xi-norm entry are unchanged. Reuse the
 # COMPLETE immutable calculation, not a conditional bound or margin.
 base=lr.result('current-coupled-first-filter-test')
 bc={k:dec(v) for k,v in base['constants'].items()}
 rec={k:bc[k] for k in ['U0','U1','U2','Utail','exi','exiw','u_error','grad_error']}
 cells=base['reconstruction_cells']
 for n,profile in profiles.items():
  assert len(profile)==len(base['profiles'][n])
  for (k,w,ff),old in zip(profile,base['profiles'][n]):
   assert k.overlaps(dec(old['degree'])) and w.overlaps(dec(old['weight'])) and ff.overlaps(dec(old['flux_factor']))
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
 # Both original rows, outer logarithms and adjoint return were executed
 # in the retained complete first receipt for these IDENTICAL profiles.
 local=base['local_return']
 lc={k:dec(v) for k,v in local['constants'].items()}
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
 Egeo=Fmax*c['actual_reference_fractional_output_coefficient']*frac['actual_fractional_flow_norm_upper']
 Emodel=Fmax*c['complete_model_density_remainder']
 tt=dict(complete_source_comparison=Esource,classical_symbol_remainder=Emodel,
  complete_local_return=lc['complete_local_return'],complete_geodesic_remainder=Egeo)
 Rt_other=sum(tt.values(),arb(0))+sum(quantum.values(),arb(0))
 filtered=cs('current-filtered-mixed-green-residual');qt=arb('.338')/(1-filtered['epsilon']);ft=qt/(1+qt)
 Yenergy=(qt.sqrt()*Es+Rt_other)/(1-ce)
 quantum['shared_vacuum_uncertainty']=vac*Yenergy;quantum['retained_energy_comparison']=ecoef*Yenergy
 tt['complete_actual_minus_reference']=sum(quantum.values(),arb(0))
 Rt=sum(tt.values(),arb(0));Rfilter=Rt+qt.sqrt()*Rs
 Kfirst=c['first_residual_to_bolt_multiplier'];Kfilter=Kfirst/ft
 screen=ft*arb('1e33');debit=Kfilter*Rfilter
 known=sc['known_partial_subtotal_not_free_margin']+Kfirst*oldfirst['complete_actual_first_residual_upper']+debit
 xi_threshold=(screen-Rt-qt.sqrt()*(scalar_debits['core_and_cutoff']+scalar_debits['scalar_remainder']))/(qt.sqrt()*Kb/beta)
 assert xi_threshold>0
 vals=dict(kappa=c['kappa'],D_d1=c['d1'],D_d0=c['d0'],D_center=c['h'],D_cluster_error=c['e'],Fref_uniform_upper=Fmax,
  model_scalar_L2=Ns,model_scalar_gradient=Gs,model_scalar_H_norm=Ds,actual_cutoff_scalar_Es=Es,
  model_difference_xi_L2=Nx,coarse_difference_xi_Hmu_upper=bc['actual_difference_xi_Hmu_upper'],
  conditional_difference_xi_Hmu_ceiling=xi_threshold,
  required_xi_comparison_improvement=bc['actual_difference_xi_Hmu_upper']/xi_threshold,
  attained_xi_comparison_improvement=bc['actual_difference_xi_Hmu_upper']/xiH,
  preceding_complete_filter_report=bc['complete_filter_energy_error'],
  complete_filter_report_improvement=bc['complete_filter_energy_error']/Rfilter,
  actual_difference_xi_Hmu_upper=xiH,complete_scalar_residual=Rs,complete_TT_residual=Rt,
  complete_filter_energy_error=Rfilter,actual_TT_trial_energy_upper=Yenergy,
  actual_filter_norm_upper=ft,filter_error_to_original_bolt_multiplier=Kfilter,
  analytic_filter_screen=screen,filter_report_over_screen=Rfilter/screen,
  filter_error_original_bolt_component=debit,whole_bolt_allowance=arb('1.667e-9'),
  known_cumulative_error_subtotal_NOT_complete_bolt=known,
  account_before_unknown_terms_NOT_margin=arb('1.667e-9')-known,
  remaining_UNEXECUTED_caller_Et_account=screen-Rfilter,
  inherited_first_solve_error_after_this_filter=ft*oldfirst['complete_actual_first_residual_upper'],
  physical_screen_for_this_filter=Kfirst*arb('1e33'),**rec)
 assert all(v.is_finite() for v in vals.values())
 # Strip historical first-source subtotal labels from reused row evaluator.
 local={k:local[k] for k in ['cell_count','cells','outer_row_norm_upper','model_row_coefficients','profile_components','profile_moments']}
 local['constants']={k:enc(lc[k]) for k in ['CU','CV','complete_classical_two_row_difference','complete_local_quantum_rows','complete_adjoint_York_return','complete_local_return']}
 result=dict(constants={k:enc(v) for k,v in vals.items()},
  scalar_residual_debits={k:enc(v) for k,v in scalar_debits.items()},
  scalar_second_order_debits={k:enc(v) for k,v in rPparts.items()},
  TT_source_comparison_debits={k:enc(v) for k,v in source_debits.items()},
  TT_residual_debits={k:enc(v) for k,v in tt.items()},quantum_debits={k:enc(v) for k,v in quantum.items()},
  spectral_cells=spectral_cells,reconstruction_cells=cells,local_return=local,
  filter_polynomial_positivity_guards=[enc(v) for v in Fguards],
  profiles={n:[dict(degree=enc(k),weight=enc(w),flux_factor=enc(ff)) for k,w,ff in v] for n,v in profiles.items()},
  reused_complete_trial_receipt='evidence/current-coupled-first-filter-test.json',
  numerical_caller_error=None,remaining_filter_error=None,complete_signed_Gt=None,
  scope=dict(complete_analytic_coupled_filter_residual=True,analytic_filter_screen_passes=bool(Rfilter<screen),
   full_common_regulator_angular_thermal_contacts_returns_retained=True,
   numerical_filter_application=False,full_bolt_certificate=False,original_P32_defect=False,
   actual_mass_column=False,original_graph_physical_admission=False,nonlinear_existence=False,
   compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,
   construction_obstruction=False,goal_complete=False))
 return ('complete_analytic_coupled_filter_passes_caller_and_signed_output_open' if Rfilter<screen else
         'complete_analytic_coupled_filter_estimator_excessive'),result

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--receipt',default=STEM);args=ap.parse_args()
 out=ROOT/'evidence'/(args.receipt+'.json');assert not out.exists(),'Preserve every receipt'
 before={p:digest(ROOT/p) for p in INPUTS};protected=json.loads(lr.BASELINE.read_bytes())
 assert all(digest(Path(p))==h for p,h in protected.items())
 runtime=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes());pkg=runtime['packages']['flint']
 assert digest(Path(sys.executable))==runtime['python_sha256']
 assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
 receipt=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  ordinary_argument_utf8=(ROOT/'notes'/(STEM+'.md')).read_text(encoding='utf-8'),
  arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),python_sha256=runtime['python_sha256'])
 try:
  receipt['status'],receipt['result']=compute();after={p:digest(ROOT/p) for p in INPUTS};assert after==before
  assert all(digest(Path(p))==h for p,h in protected.items())
  receipt.update(input_sha256_after=after,protected_paths_unchanged=len(protected),pinned_flint_files_unchanged=len(pkg['files']))
 except Exception:receipt.update(status='failed',traceback=traceback.format_exc())
 receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
 out.write_bytes((json.dumps(receipt,indent=2)+'\n').encode('utf-8'))
 print(json.dumps(dict(status=receipt['status'],error=receipt.get('traceback'),
  constants={k:v['ball'] for k,v in receipt.get('result',{}).get('constants',{}).items()},
  scalar_debits={k:v['ball'] for k,v in receipt.get('result',{}).get('scalar_residual_debits',{}).items()},
  TT_debits={k:v['ball'] for k,v in receipt.get('result',{}).get('TT_residual_debits',{}).items()}),indent=2),flush=True)
 if receipt['status']=='failed':raise SystemExit(1)
