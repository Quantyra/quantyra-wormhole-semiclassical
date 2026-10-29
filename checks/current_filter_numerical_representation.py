"""Complete caller error for a finite product-preserving filter representation."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction
from math import comb
import json,sys,argparse,traceback,ast,inspect
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_coupled_first_filter_test as old
arb,ctx,dec,enc,J=old.arb,old.ctx,old.dec,old.enc,old.J
from flint import acb
ROOT=old.ROOT;lr=old.lr;fr=old.fr;cs=old.cs;up=old.up
sys.path.insert(0,str(ROOT/'operators'))
from paired_filter_profile import PairedFilterProfile
STEM='current-filter-numerical-representation'
INPUTS=list(dict.fromkeys(['notes/'+STEM+'.md','checks/current_filter_numerical_representation.py',
 'operators/paired_filter_profile.py','evidence/current-coupled-first-filter-test-attempt2.json',
 'evidence/current-filter-numerical-representation-attempt2.json',
 'evidence/current-coupled-first-filter-test-attempt2-verification.json',*old.INPUTS]))
digest=lambda p:sha256(p.read_bytes()).hexdigest()
STAGES={}

def qball(q):return arb(q.numerator)/arb(q.denominator)
def qenc(q):return dict(numerator=str(q.numerator),denominator=str(q.denominator))
def middle(v):
 m,_,e=v.mid_rad_10exp(110);return Fraction(int(m))*Fraction(10)**int(e)
def cmiddle(v):return middle(v.real),middle(v.imag)
def cball(v):return acb(qball(v[0]),qball(v[1]))
def cenc(v):return dict(real=qenc(v[0]),imag=qenc(v[1]))
def polyadd(*vs):
 out=[0]*max(map(len,vs))
 for v in vs:
  for j,a in enumerate(v):out[j]+=a
 return out
def polymul(v,w):
 out=[0]*(len(v)+len(w)-1)
 for j,a in enumerate(v):
  for k,b in enumerate(w):out[j+k]+=a*b
 return out
def polyeval(v,x):
 out=0
 for a in reversed(v):out=out*x+a
 return out
def polypoles(poles):
 out=[1]
 for pole in poles:out=polymul(out,[pole,1])
 return out

def output_comparison(p,c,bs):
 L,a,ga=p['L'],p['a'],p['gamma'];lr.model_rows(p,c)
 base=lr.result('current-coupled-first-filter-test-attempt2')['local_return']
 maxima=[[arb(0)]*8 for _ in range(4)];rows=[]
 for i,cell in enumerate(base['cells']):
  s0,s1=dec(cell['s0']),dec(cell['s1']);x0,x1=s0/ga,s1/ga
  x=J([arb((x0+x1)/2,(x1-x0)/2),arb(1)]+[arb(0)]*8)
  rr=lr.rows(p,c,bs,x)
  for h in range(4):
   for j in range(8):maxima[h][j]=lr.maximum(maxima[h][j],up(rr[h][j]))
  rows.append([[enc(up(v)) for v in row] for row in rr])
  if i%40==0:print(json.dumps(dict(uniform_row_cells=i,total=len(base['cells']))),flush=True)
 outerA=[[arb(0)]*8 for _ in range(4)];outerB=[[arb(0)]*8 for _ in range(4)]
 for i in range(64):
  theta=arb((2*i+1)*arb.pi()/64,arb.pi()/64)
  z=acb(theta.cos(),theta.sin())*arb('.02')
  x=J([(L-z)/ga,acb(1)]+[acb(0)]*8)
  r0=lr.rows(p,c,bs,x,acb(0));r1=lr.rows(p,c,bs,x,acb(1))
  for h in range(4):
   for j in range(8):
    outerA[h][j]=lr.maximum(outerA[h][j],abs(r0[h][j]))
    outerB[h][j]=lr.maximum(outerB[h][j],abs(r1[h][j]-r0[h][j]))
 b=arb('.01');outer=[]
 for h in range(4):
  row=[]
  for j in range(8):
   v=(b/arb('.02'))**2*(outerA[h][j]+abs(b.log())*outerB[h][j])
   maxima[h][j]=lr.maximum(maxima[h][j],v);row.append(enc(v))
  outer.append(row)
 Dmax=cs('current-source-adapted-local-return-attempt4')['Dmax']
 flux=(2*ga*L+a*L*L).sqrt()/ga
 bb=[arb(1),flux,arb(1),flux,arb(1),flux,arb(1),L*Dmax]
 norms=[sum((v*w for v,w in zip(row,bb)),arb(0)) for row in maxima]
 for j,name in enumerate(['fn','fu']):
  norms[j]+=L**2*sum((up(c[name+'_'+k]) for k in ['r','f','g']),arb(0))
 CU,CV=[dec(base['constants'][k]) for k in ['CU','CV']]
 local=CU*(norms[0]+norms[2])+CV*(norms[1]+norms[3])
 shape=lr.result('current-source-adapted-fractional-return');fc=cs('current-source-adapted-fractional-return')
 S=[dec(v) for v in shape['derivative_accounting']['reconstructed_reference_derivative_constants']]
 D0=dec(shape['connection_tensor_derivative_bounds'][0])
 U0=S[0];U1=arb(10).sqrt()*S[1];U2=10*(S[2]+D0*S[1])
 oc=cs('current-source-adapted-complete-response-test');C,B2,beta=p['C'],p['B2'],p['beta']
 CG,CD=oc['test_gradient_constant'],oc['test_Laplacian_constant']
 factor=2/(B2*beta);qang=arb(2000)**4;pi=arb.pi();K=arb('1e10');sp=arb('1e-24')
 br=(sp/(4*K)).sqrt();bg=K*(arb(1)/4+2*((arb(1)/400)/br).log())
 qparts=dict(local_gain_two=factor*2*pi*C*(sp*U2**2/3+bg*U1**2).sqrt()*(sp*(3*CD**2/8+CG**2/4)/3+bg*CG**2).sqrt(),
  full_annulus=factor*2*pi*arb('4e8')*C*(arb(400)/128).log()*U1*CG,
  finite_second_contact=factor*arb('1e13')*C*U1*CG,
  global_row_return=factor*oc['complete_global_gluing_coefficient']*C*U0*CV,
  complete_cap_contact=factor*arb('1e18')*qang**2/C*(1+1/arb(5).sqrt())*U0*CV,
  constant_contact_return=factor*oc['constant_contact_gluing_coefficient']*C*U0*CV,
  curvature_polynomial=factor*100*C*(U1*CG/(36*pi)+816*U0*CV/(960*pi)))
 eps=cs('current-complete-polynomial-weak-transfer')['epsilon']
 ce=arb('1e-21')/((1-eps)*arb('.99')*(1-arb('1.179e-8')))+eps/(1-eps)
 assert ce<1
 geo=c['actual_reference_fractional_output_coefficient']*fc['actual_reconstruction_fractional_comparison']
 CE=up((local+geo+sum(qparts.values(),arb(0)))/(1-ce))
 return CE,dict(constants={k:enc(v) for k,v in dict(complete_local_dual=local,
  complete_geodesic_dual=geo,energy_feedback=ce,original_X_to_Et=CE).items()},
  quantum_entries={k:enc(v) for k,v in qparts.items()},
  whole_cell_coefficients=rows,outer_coefficients_after_factoring_z_squared=outer,
  global_coefficient_maxima=[[enc(v) for v in row] for row in maxima],
  basis_norm_comparisons=[enc(v) for v in bb],row_norms=[enc(v) for v in norms])

def profile_parameters(p,c):
 actual={k:c[k] for k in ['K0','kappa','Cv','delta','muT','muS','d1','d0']}
 actual.update(lambda_l=p['a']*p['ell']*(p['ell']+1),lambda_n=c['large_mass_lambda'],At=p['At'])
 nominal={k:middle(v) for k,v in actual.items()}
 def polys(v):
  N=polymul([v['Cv'],v['delta']],[v['Cv']*v['lambda_n'],1])
  N=[v['K0']*v['kappa']*x for x in N]
  D=polymul(polypoles([v[k] for k in ['muT','muS','lambda_l','lambda_n']]),[v['d0'],v['d1'],1])
  return N,D
 return actual,nominal,polys(actual),polys(nominal)

def parameter_error(p,actual,nominal,ad,nd):
 Na,Da=ad;Da=[arb(v) for v in Da];Nn,Dn=nd;amplitude=qball(nominal['At'])
 ratio=actual['At']/amplitude
 diff=polyadd([ratio*v for v in polymul(Na,[qball(v) for v in Dn])],
              [-v for v in polymul([qball(v) for v in Nn],Da)])
 assert len(diff)==9
 gard=lambda d:min((arb(v).lower()/comb(6,j) for j,v in enumerate(d)))
 Dnb=[qball(v) for v in Dn];ga=arb(gard(Da));gn=arb(gard(Dnb));assert ga>0 and gn>0
 total=arb(0);moments=[arb(0)]*3;cells=[];count=1024
 for i in range(count):
  tl=arb(i)/count;tu=arb(i+1)/count;wl=1-tu;wu=1-tl
  # Exact integration endpoints are rational; use outward positive guards.
  tl=arb(max(arb(0),tl.lower()));wl=arb(max(arb(0),wl.lower()))
  tu=up(tu);wu=up(wu)
  dl=sum((arb(v.lower())*tl**j*wl**(6-j) for j,v in enumerate(Da)),arb(0))
  nl=sum((arb(v.lower())*tl**j*wl**(6-j) for j,v in enumerate(Dnb)),arb(0))
  dl=arb(max(ga.lower(),dl.lower()));nl=arb(max(gn.lower(),nl.lower()))
  num=sum((up(v)*tu**j*wu**(8-j) for j,v in enumerate(diff)),arb(0))
  cost=up(wu**arb('1.5')*(num/(dl*nl))**2/count)
  for j in range(3):moments[j]+=up(tu**(2*j)*wu**(6-2*j)*(num/(dl*nl))**2/count)
  total+=cost;cells.append(enc(cost))
 X=up(amplitude*(p['gamma']*total/2).sqrt())
 return X,dict(cross_numerator=[enc(v) for v in diff],whole_spectral_cells=cells,
  denominator_Bernstein_guards=[enc(ga),enc(gn)],parameter_X_error=enc(X),
  model_H_power_norms=[enc(up(amplitude*(p['gamma']*v/2).sqrt())) for v in moments])

def reuse_output_comparison():
 raw=json.loads((ROOT/'evidence/current-filter-numerical-representation-attempt2.json').read_bytes())
 assert raw['status']=='failed' and 'output_comparison' in raw['partial_results']
 tree=ast.parse(raw['arithmetic_source_utf8'])
 previous=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='output_comparison')
 assert ast.get_source_segment(raw['arithmetic_source_utf8'],previous).strip()==inspect.getsource(output_comparison).strip()
 changed={'notes/'+STEM+'.md','checks/current_filter_numerical_representation.py','operators/paired_filter_profile.py'}
 assert all(digest(ROOT/p)==h for p,h in raw['input_sha256_before'].items() if p not in changed)
 result=raw['partial_results']['output_comparison']
 return dec(result['constants']['original_X_to_Et']),result

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
    dg=gap*gap-abs(e);assert gap>0 and dg>0;den*=dg
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
 Bpoly=polypoles(B);Gnum=polymul(Bpoly,Nn);ginf=Gnum[-1]
 ln=nominal['lambda_n'];other=polymul(polypoles([nominal[k] for k in ['muT','muS','lambda_l']]),[nominal['d0'],nominal['d1'],1])
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
   dg=gap*gap-abs(e);assert dg>0;den*=dg
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
 ctx.prec=768;p,c,bs=fr.load();lr.model_rows(p,c);c=old.setup(p,c)
 # Pin branch convention against elementary exact solutions.
 w=acb(3);q0=w.legendre_q(0,0,type=3);q1=w.legendre_q(1,0,type=3)
 assert q0.overlaps(((w+1)/(w-1)).log()/2) and q1.overlaps(w*q0-1)
 CE,comparison=reuse_output_comparison()
 STAGES['output_comparison']=comparison
 print(json.dumps(dict(complete_original_output_comparison=str(CE))),flush=True)
 actual,nominal,ad,nd=profile_parameters(p,c)
 Xpar,pardata=parameter_error(p,actual,nominal,ad,nd)
 STAGES['parameter_error']=pardata
 Epar,direct=direct_parameter_output(p,c,actual,nominal,pardata,comparison,Xpar)
 STAGES['direct_parameter_output']=direct
 field,Xquad,Xround,contours=contour_representation(p,nominal,nd)
 STAGES['numerical_field']=field;STAGES['contour_bounds']=contours
 generic=up(CE*(Xpar+Xquad+Xround));err=up(Epar+CE*(Xquad+Xround));target=arb('1.70e31')
 first=cs('current-coupled-first-filter-test-attempt2');K=first['filter_error_to_original_bolt_multiplier']
 print(json.dumps(dict(comparison=str(CE),parameter_error=str(Xpar),quadrature_error=str(Xquad),
  original_energy_caller_error=str(err))),flush=True)
 # Queries act on the returned finite-product field, not on an interpolant.
 model=PairedFilterProfile(p,bs,field);queries=[]
 for label,s,bolt in [('bolt',arb(0),True),('thermal_scale',p['gamma'],False),
                       ('scalar_cutoff_scale',arb('1e-12'),False),('cold_split',arb('.0004'),False)]:
  vals=model.at(s,bolt=bolt)
  assert all(v.is_finite() for v in vals.values())
  queries.append(dict(label=label,s=enc(s),values={k:enc(v) for k,v in vals.items()}))
  print(json.dumps(dict(evaluated=label,xi=str(vals['xi']),u=str(vals['u']))),flush=True)
 queries.append(dict(label='original_outer_boundary',values={k:enc(v) for k,v in model.at_outer().items()}))
 vals=dict(original_X_to_Et_comparison=CE,parameter_X_error=Xpar,contour_X_error=Xquad,
  direct_parameter_Et_error=Epar,generic_caller_Et_upper=generic,generic_report_over_target=generic/target,
  direct_parameter_effective_comparison=Epar/Xpar,generic_to_direct_report_gain=generic/err,
  coefficient_and_node_X_error=Xround,complete_numerical_caller_Et_error=err,caller_target=target,
  report_over_caller_target=err/target,caller_original_bolt_error=K*err,
  complete_analytic_plus_caller_filter_error=first['complete_filter_energy_error']+err,
  combined_filter_screen=first['analytic_filter_screen'],
  remaining_filter_account=first['analytic_filter_screen']-first['complete_filter_energy_error']-err,
  cumulative_known_error_subtotal_NOT_complete_output=first['known_cumulative_error_subtotal_NOT_complete_bolt']+K*err,
  conditional_required_X_error=target/CE,required_comparison_if_X_fixed=target/(Xpar+Xquad+Xround))
 passed=err<target and vals['complete_analytic_plus_caller_filter_error']<vals['combined_filter_screen']
 return ('complete_finite_product_filter_caller_passes' if passed else 'complete_finite_product_filter_caller_estimator_excessive'),dict(
  precision_bits=ctx.prec,constants={k:enc(v) for k,v in vals.items()},output_comparison=comparison,
  parameter_error=pardata,direct_parameter_output=direct,nominal_profile_parameters={k:qenc(v) for k,v in nominal.items()},
  contour_bounds=contours,numerical_field=field,point_and_endpoint_queries=queries,
  numerical_state='finite exact-rational resolvent-product field; no sampled interpolation',
  full_signed_Gt=None,remaining_filter_error=None,
  scope=dict(actual_paired_domain=True,full_quantum_regulator_contacts_and_returns=True,
   complete_numerical_representation_error=True,numerical_caller_target_passes=bool(passed),
   arbitrary_interpolant_certified=False,full_bolt_certificate=False,original_P32_defect=False,
   mass_column=False,original_physical_graph_admission=False,nonlinear_existence=False,
   compatibility=False,physical_matching=False,validity_noise=False,all_sector_stability=False,
   construction_obstruction=False,goal_complete=False))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--receipt',default=STEM);args=ap.parse_args()
 out=ROOT/'evidence'/(args.receipt+'.json');assert not out.exists(),'Preserve immutable evidence'
 before={p:digest(ROOT/p) for p in INPUTS};protected=json.loads(lr.BASELINE.read_bytes())
 assert all(digest(Path(p))==h for p,h in protected.items())
 rt=json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes());pkg=rt['packages']['flint']
 assert digest(Path(sys.executable))==rt['python_sha256']
 assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
 receipt=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
  ordinary_argument_utf8=(ROOT/'notes'/(STEM+'.md')).read_text(encoding='utf-8'),
  arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),
  numerical_operator_source_utf8=(ROOT/'operators/paired_filter_profile.py').read_text(encoding='utf-8'),
  python_sha256=rt['python_sha256'])
 try:
  receipt['status'],receipt['result']=execute();after={p:digest(ROOT/p) for p in INPUTS};assert after==before
  assert all(digest(Path(p))==h for p,h in protected.items())
  receipt.update(input_sha256_after=after,protected_paths_unchanged=len(protected),pinned_flint_files_unchanged=len(pkg['files']))
 except Exception:receipt.update(status='failed',traceback=traceback.format_exc(),partial_results=STAGES)
 receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
 out.write_bytes((json.dumps(receipt,indent=2)+'\n').encode('utf-8'))
 print(json.dumps(dict(status=receipt['status'],error=receipt.get('traceback'),
  constants={k:v['ball'] for k,v in receipt.get('result',{}).get('constants',{}).items()}),indent=2),flush=True)
 if receipt['status']=='failed':raise SystemExit(1)
