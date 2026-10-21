"""Literal first-order conformal/shape mixing; no physical branch exclusion."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1];STEM='first-order-mixed-boundary-response'
OWN=['checks/first_order_mixed_boundary_response.py','notes/'+STEM+'.md','evidence/first-order-mixed-source-review.json']
PREVIOUS='coupled-fast-quadratic-boundary-obstruction';ROWS=[];FIXTURES=[];R=S.Rational
PROTOCOL={'precision_digits':[50,80],'charged_radii':['80','160','320'],
 'charged_cases':[['2','.01','1.1','2','2'],['3','-.02','1.2','3','5'],['4','1e-14','1.05','5','9']],
 'charged_columns':['lambda','d','horizon','magnetic_flux','k'],'contour_nodes':[12,16],
 'contour_radius':'.05','contour_scaled_limit':'1e-30','precision_scaled_limit':'1e-35',
 'normal_betas':['2','2.25','2.75'],'normal_shape_frequencies':['3','4','6'],
 'normal_annuli':[0,1,2],'quadrature_scaled_limit':'1e-35','spike_indices':[2,4,8],
 'static_fixtures':True,'thermal_transfer':'ordinary derivative count using full boundary heat coefficient',
 'physical_benchmark_evaluated':False,'canonical_determinant_evaluated':False,
 'physical_action_averaging_adopted':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,value):
 value=S.factor(S.cancel(S.sympify(value).doit()));check(name,value==0,value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==201 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 rev=json.loads((ROOT/OWN[2]).read_bytes())
 check('retained_inputs_pinned',all(f in p['input_sha256_after'] for f in rev['retained_inputs']))
 check('scope',not rev['physical_fast_branch_excluded'] and not rev['pure_shape_finite_cut_limit_proved'] and not rev['new_boundary_prescription_adopted'] and not rev['physical_stability_proved'])

def charged_identities():
 r,lam,G,qm=S.symbols('r lam G qm',positive=True)
 fj=S.symbols('F:6',real=True);p=S.symbols('p:4');s=S.symbols('s:5');v=S.symbols('v:5');w=S.symbols('w:3')
 F,Fr,Frr=fj[:3];a,b=S.symbols('a b')
 def D(expr):
  return S.diff(expr,r)+sum(S.diff(expr,x[j])*x[j+1] for x in [fj,p,s,v] for j in range(len(x)-1))
 def Dw(expr):return D(expr)+sum(S.diff(expr,w[j])*w[j+1] for j in range(2))
 speed=r*S.sqrt(F)*S.exp(-v[0]);N=S.exp(v[0])*S.sqrt(F)/r
 aa=S.exp(-v[0])*((r*Fr/2-F)/S.sqrt(F)+r*S.sqrt(F)*v[1]);ax=speed*D(aa)
 ss=S.exp(-v[0])*S.sqrt(F)*(1+r*w[1]);sx=speed*Dw(ss)
 Rh=S.expand(-2*(ax+aa*aa));u=S.expand(sx+aa*ss);V=ss*ss;sig=S.log(r)+w[0]
 Rs=-Frr-4*Fr/r+2*(1-F)/r**2;HR=2*Rs+r*D(Rs);Cv=2*Frr+8*Fr/r+4*F/r**2
 Rbase=S.exp(-2*v[0])*(Rs-2/r**2-2*F*v[2]-2*Fr*v[1])+2/r**2
 lam0=lam+11*G*lam**2/(360*S.pi)
 LC=S.exp(2*v[0])*r*r*(S.exp(2*w[0])*(Rbase-6*S.exp(-2*v[0])*(F*w[2]+(Fr+2*F/r)*w[1]+F*w[1]**2))+6*lam0*S.exp(4*w[0]))/(4*G)
 LW=S.exp(2*v[0])/(4*S.pi*r*r)*(sig*((Rh+2)**2/60-R(11,90)*Rh-qm*qm/3)+Rh*u/30-u*u/10+(Rh/30-R(1,18))*V-R(7,90)*u*V-R(7,180)*V*V)
 GG=sig*(12*aa*ax+R(44,7)*aa**3+10*aa)-6*ax*ss-R(24,7)*aa*aa*ss-6*aa*sx-R(96,7)*aa*ss*ss-18*ss*sx-R(85,21)*ss**3
 GHY=S.exp(2*w[0])*r*r*(Fr/2+2*F/r+F*(3*w[1]+v[1]))/(2*G)
 CT=-S.sqrt(F)*(S.sqrt(lam)*r*r*S.exp(3*w[0]+v[0])+S.exp(w[0]+v[0])/(2*S.sqrt(lam)))/G
 c0=-221*lam**R(3,2)/(15120*S.pi**2);c1=29*S.sqrt(lam)/(20160*S.pi**2)
 loop=4*S.pi*S.sqrt(F)*(c0*r*r*S.exp(3*w[0]+v[0])+2*c1*S.exp(w[0]+v[0]))
 BD=GHY+CT+loop-N*GG/(360*S.pi)
 sub={**{w[j]:a*p[j]+b*s[j] for j in range(3)},**{v[j]:b*v[j] for j in range(3)}}
 cross=lambda f:S.expand(S.diff(f.subs(sub,simultaneous=True),a,b).subs({a:0,b:0}))
 lc=cross(LC+LW);bc=cross(BD)
 mom2=S.diff(lc,p[2]);mom1=S.diff(lc,p[1])-D(mom2)
 end=S.expand(bc+p[0]*mom1+p[1]*mom2)
 rep={s[2]:-(Fr+2*F/r)*s[1]/F-HR*s[0]/(6*F)-v[2]/3-Fr*v[1]/(3*F)+Cv*v[0]/(6*F)}
 rep[s[3]]=S.expand(D(rep[s[2]]).subs(rep))
 er=S.factor(end.subs(rep,simultaneous=True))
 for j in [s[2],s[3],v[3]]:zero('charged_high_jet_cancellation_'+str(j),S.diff(er,j))
 J=-F*p[0]*r*(r*Fr-2*F)/(60*S.pi)
 zero('charged_remaining_second_jet',S.diff(er,v[2])-J)
 check('nonvacuum_second_jet_not_identically_zero',J!=0)
 zero('exact_finite_annulus_identity',lc-p[0]*(S.diff(lc,p[0])-D(mom1))-D(p[0]*mom1+p[1]*mom2))
 rep[s[4]]=S.expand(D(rep[s[3]]).subs(rep,simultaneous=True))
 row=S.expand(S.diff(lc,p[0])-D(mom1)).subs(rep,simultaneous=True)
 eps=G/(60*S.pi);Wc=-Frr+2*Fr/r+2*(1-F)/r**2;j=-Frr/4+(F-1)/(2*r*r)
 dR=7*Wc/(3*r)+44*j/(3*r)+D(Rs)
 eR=44*j/(3*r*r)+D(D(Rs))+2*D(Rs)/r
 T=Rs+12*lam0+eps*(7*Wc*Wc/36+22*j*j/3-11*Rs*Rs/72-10*qm*qm/r**4+F*D(D(Rs))+(Fr+2*F/r)*D(Rs))
 A1=r*F*dR;A0=r*F*(dR*Fr+eR*F);B1=2*r*F*eR
 LL=Fr+r*Frr+r*HR/6;B0=2*eR*F*F+r*F*eR*Fr+F*dR*LL
 Uv=r*(-2*(F*v[2]+Fr*v[1])+Cv*v[0]);qw=r*s[0];ww=r*v[0]
 gw=A1*Uv/6+B1*F*D(qw)-B0*qw-2*A0*v[0]
 dw=eps*gw+F*((r*D(T)+4*T)*qw+2*T*ww)
 zero('independent_residual_trace_row',row-r*dw/(2*G*F))
 zero('magnetic_bulk_mixed_coefficient',S.diff(S.exp(2*b*v[0])*(-S.pi*G*qm**2)/(2*G*r*r),a,b))
 power=S.symbols('power',real=True)
 normalF={fj[j]:S.diff(lam*r*r,r,j) for j in range(6)}
 scaledp={p[j]:S.prod(-R(3,2)-i for i in range(j))*r**(-R(3,2)-j) for j in range(3)}
 for field,scale in [(s[0],R(3,2)),(v[0],R(3,2)),(s[1],R(5,2)),(v[1],R(5,2))]:
  value=S.simplify(S.limit(S.diff(er,field).subs(normalF,simultaneous=True).subs(scaledp,simultaneous=True)/r**scale,r,S.oo))
  check('low_jet_coefficient_growth_'+str(field),not value.has(S.oo,-S.oo,S.zoo,S.nan),value)
 m,Q,c3,c4=S.symbols('m Q c3 c4')
 charged=lam*r*r+1-2*m/r+Q/r**2+c3/r**3+c4/r**4
 zero('charged_seed_radial_mismatch',r*S.diff(charged,r)-2*charged-(-2+6*m/r-4*Q/r**2-5*c3/r**3-6*c4/r**4))
 zero('curvature_free_second_derivative_coefficient',S.diff(rep[s[2]],v[2])+R(1,3))
 # The temporal term in reconstruction changes only zeroth jets.
 gam=S.symbols('gam');temporal=gam**2*(s[0]+v[0]/3)/F**2
 zero('thermal_high_jet_coefficient_unchanged',S.diff(temporal,v[2]))
 args=[r,lam,G,qm,*fj[:3],*p[:3],*s[:3],*v[:3]]
 coeff=S.lambdify(args,[lc,bc],'mpmath')
 raw=S.lambdify([*args,a,b],[(LC+LW).subs(sub,simultaneous=True),BD.subs(sub,simultaneous=True)],'mpmath')
 return coeff,raw

def normal_identities():
 z,t=S.symbols('z t',real=True)
 p=S.symbols('p:5'); s=S.symbols('s:5'); v=S.symbols('v:5'); w=S.symbols('w:3')
 aa,bb=S.symbols('aa bb')
 def D(f):
  return S.diff(f,z)+sum(S.diff(f,a[j])*a[j+1] for a in [p,s,v] for j in range(4))
 u=1+w[1]+w[2];q=(1+w[1])**2;b=v[1]+v[2]
 LC=S.exp(2*w[0])*(-12-2*b-6*(w[2]+3*w[1]+w[1]**2))/4+R(3,2)*(1+11*t/6)*S.exp(4*w[0]+2*v[0])
 LW=15*t*S.exp(-2*v[0])*((z+w[0])*b*b/15-u*b/15-u*u/10-b*q/15-7*u*q/90-7*q*q/180)
 ax=v[2]+v[1]-v[1]**2;sx=(1-v[1])*(1+w[1])+w[2]
 GG=(z+w[0])*(12*v[1]*ax+R(44,7)*v[1]**3)-6*ax*(1+w[1])-R(24,7)*v[1]**2*(1+w[1])-6*v[1]*sx-R(96,7)*v[1]*(1+w[1])**2-18*(1+w[1])*sx-R(85,21)*(1+w[1])**3
 BD=S.exp(2*w[0])*(3+3*w[1]+v[1])/2-(1+R(221,63)*t)*S.exp(3*w[0]+v[0])-t*S.exp(-2*v[0])*GG/6
 sub={**{w[j]:aa*p[j]+bb*s[j] for j in range(3)},**{v[j]:bb*v[j] for j in range(3)}}
 cross=lambda f:S.expand(S.diff(f.subs(sub,simultaneous=True),aa,bb).subs({aa:0,bb:0}))
 lc=cross(LC+LW);bc=cross(BD)
 mom2=S.diff(lc,p[2]);mom1=S.diff(lc,p[1])-D(mom2)-3*mom2
 row=S.diff(lc,p[0])-D(mom1)-3*mom1
 end=S.expand(bc+p[0]*mom1+p[1]*mom2)
 constraint=3*s[2]+9*s[1]-12*s[0]+v[2]+v[1]-12*v[0]
 rep={s[2]:(12*s[0]-9*s[1]-v[2]-v[1]+12*v[0])/3}
 rep[s[3]]=S.expand(D(rep[s[2]]).subs(rep))
 rep[s[4]]=S.expand(D(rep[s[3]]).subs(rep))
 er=S.factor(end.subs(rep,simultaneous=True))
 zero('normal_mixed_row_zero',row.subs(rep,simultaneous=True))
 target=-(3+221*t/7)*p[0]*s[0]+((3+11*t)*p[0]+33*t*p[1]/7+3*t*p[2])*s[1]
 target+=(-(3+125*t/7)*p[0]-54*t*p[1]/7-6*t*p[2])*v[0]+((1+11*t/3)*p[0]-3*t*p[1]/7+t*p[2])*v[1]
 zero('normal_complete_endpoint',er-target)
 for field in [s[2],s[3],v[2],v[3]]:zero('normal_no_high_jet_'+str(field),S.diff(er,field))
 args=[z,t,*p[:3],*s[:3],*v[:3]]
 return S.lambdify(args,[lc,bc,er],'mpmath')


def endpoint_controls():
 x,d=S.symbols('x d',positive=True)
 zero('H2_value_bound_integral',S.integrate((d-x)**2,(x,0,d))-d**3/3)
 zero('logarithmic_L1_weight',S.integrate(x**3,(x,0,d))-d**4/4)
 zero('pure_conformal_bulk_cross_power',-4+R(3,2)+R(5,2))
 zero('pure_conformal_boundary_cross_power',-3+R(3,2)+R(5,2)-1)
 zero('thermal_first_jet_power',R(3,2)+R(1,2)-2)
 psi=x*x*(1-x*x)**6/2
 zero('spike_center_second_jet',S.diff(psi,x,2).subs(x,0)-1)
 for j in range(6):zero('spike_compact_join_'+str(j),S.diff(psi,x,j).subs(x,1))
 norm=[S.integrate(S.diff(psi,x,j)**2,(x,-1,1)) for j in range(3)]
 for j,n in enumerate(norm):check('spike_norm_positive_'+str(j),n>0,n)
 zero('spike_H2_energy_power',2*R(27,2)-3*8-3)
 zero('spike_weighted_trace_power',R(5,2)+R(27,2)-2*8)
 zero('spike_weighted_H1_derivative_power',2*(R(27,2)-8-1)+8-17)
 zero('spike_energy_geometric_sum',S.summation(S.Rational(1,8)**S.Symbol('n',integer=True),(S.Symbol('n',integer=True),2,S.oo))-R(1,56))
 return S.lambdify(x,S.diff(psi,x,2),'mpmath'),[mp.mpf(str(S.N(n,90))) for n in norm]

def fall(x,n):
 out=1
 for j in range(n):out*=x-j
 return out

def charged_fixtures(coeff,raw):
 for ix,row in enumerate(PROTOCOL['charged_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    lam,d,rh,qm,k=map(mp.mpf,row);t=1/(k*k-mp.mpf(17)/12);G=60*mp.pi*t/lam
    Q=1+3*lam+d;m=rh*(1+lam*rh**2+Q/rh**2-2*d/rh**3+d/rh**4)/2
    records=[]
    for ri,rstring in enumerate(PROTOCOL['charged_radii']):
     r=mp.mpf(rstring);Fj=[sum(c*fall(n,j)*r**(n-j) for c,n in [(lam,2),(1,0),(-2*m,-1),(Q,-2),(-2*d,-3),(d,-4)]) for j in range(4)]
     F,Fr,Frr,Frrr=Fj;p=[mp.re(fall(-mp.mpf('1.5')+1j*k,j)*r**(-mp.mpf('1.5')+1j*k-j)) for j in range(3)]
     nu=-mp.mpf('2.25')+2j*k;v=[mp.re(fall(nu,j)*r**(nu-j)) for j in range(3)]
     ss=[mp.re(-(nu-3)/(3*(nu-1))*fall(nu,j)*r**(nu-j)) for j in range(2)]
     Rseed=-Frr-4*Fr/r+2*(1-F)/r**2;Rprime=-Frrr-4*Frr/r+2*Fr/r**2-4*(1-F)/r**3
     HR=2*Rseed+r*Rprime;Cv=2*Frr+8*Fr/r+4*F/r**2
     ss.append(-(Fr+2*F/r)*ss[1]/F-HR*ss[0]/(6*F)-v[2]/3-Fr*v[1]/(3*F)+Cv*v[0]/(6*F))
     args=[r,lam,G,qm,*Fj[:3],*p,*ss,*v];exact=coeff(*args);contours=[]
     for n in PROTOCOL['contour_nodes']:
      radius=mp.mpf(PROTOCOL['contour_radius']);acc=[0j,0j]
      for j in range(n):
       aa=radius*mp.exp(2j*mp.pi*j/n)
       for h in range(n):
        bb=radius*mp.exp(2j*mp.pi*h/n);value=raw(*args,aa,bb)
        for ch in range(2):acc[ch]+=value[ch]/(aa*bb*n*n)
      errors=[abs(acc[ch]-exact[ch])/(1+abs(exact[ch])) for ch in range(2)]
      tag=str(ix)+'_'+str(digits)+'_'+str(ri)+'_'+str(n)
      check('charged_unexpanded_mixed_'+tag,max(errors)<mp.mpf(PROTOCOL['contour_scaled_limit']),mp.nstr(max(errors),24))
      contours.append({'nodes_per_amplitude':n,'scaled_errors':[mp.nstr(e,25) for e in errors]})
     records.append({'radius':rstring,'coefficients':[mp.nstr(x,digits) for x in exact],'contours':contours})
    levels.append({'precision':digits,'values':records})
  with mp.workdps(80):
   err=max(abs(mp.mpf(a)-mp.mpf(b))/(1+abs(mp.mpf(b))) for x,y in zip(levels[0]['values'],levels[1]['values']) for a,b in zip(x['coefficients'],y['coefficients']))
   check('charged_precision_'+str(ix),err<mp.mpf(PROTOCOL['precision_scaled_limit']),mp.nstr(err,25))
  FIXTURES.append({'fixture':'charged_original_mixed_Taylor','case':row,'levels':levels,'physical_solution':False})

def normal_fixtures(coeff):
 for ix,row in enumerate(PROTOCOL['charged_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    k=mp.mpf(row[-1]);t=1/(k*k-mp.mpf(17)/12);beta=mp.mpf(PROTOCOL['normal_betas'][ix]);om=mp.mpf(PROTOCOL['normal_shape_frequencies'][ix]);nu=-beta+1j*om
    def values(z):
     p=[mp.re((-mp.mpf('1.5')+1j*k)**j*mp.exp((-mp.mpf('1.5')+1j*k)*z)) for j in range(3)]
     v=[mp.re(nu**j*mp.exp(nu*z)) for j in range(3)];s=[mp.re(-(nu-3)/(3*(nu-1))*nu**j*mp.exp(nu*z)) for j in range(3)]
     return [mp.exp(3*z)*a for a in coeff(z,t,*p,*s,*v)]
    records=[]
    for n in PROTOCOL['normal_annuli']:
     left=mp.log(80)+n*mp.pi/k;right=left+mp.pi/(4*k)
     integral=mp.quad(lambda z:values(z)[0],[left,(left+right)/2,right])
     actual=integral+values(right)[1]-values(left)[1];expected=values(right)[2]-values(left)[2]
     err=abs(actual-expected)/(1+abs(expected));tag=str(ix)+'_'+str(digits)+'_'+str(n)
     check('normal_annular_identity_'+tag,err<mp.mpf(PROTOCOL['quadrature_scaled_limit']),mp.nstr(err,25))
     records.append({'annulus':n,'actual':mp.nstr(actual,digits),'expected':mp.nstr(expected,digits),'scaled_error':mp.nstr(err,25)})
    levels.append({'precision':digits,'values':records})
  with mp.workdps(80):
   err=max(abs(mp.mpf(a['actual'])-mp.mpf(b['actual']))/(1+abs(mp.mpf(b['actual']))) for a,b in zip(levels[0]['values'],levels[1]['values']))
   check('normal_precision_'+str(ix),err<mp.mpf(PROTOCOL['precision_scaled_limit']),mp.nstr(err,25))
  FIXTURES.append({'fixture':'normal_curvature_free_annular_action','beta':PROTOCOL['normal_betas'][ix],'levels':levels,'physical_solution':False})

def spike_fixtures(psi2,norms):
 levels=[]
 for digits in PROTOCOL['precision_digits']:
  with mp.workdps(digits):
   rows=[]
   for n in PROTOCOL['spike_indices']:
    d=mp.mpf(2)**(-n);ell=d**8;amp=d**mp.mpf('13.5')
    energy=amp**2/ell**3*mp.quad(lambda x:psi2(x)**2,[-1,0,1])
    logmass=mp.quad(lambda x:abs((d+ell*x)**mp.mpf('1.5')*amp/ell**2*psi2(x))*ell,[-1,-mp.mpf('.75'),-mp.mpf('.5'),-mp.mpf('.25'),0,mp.mpf('.25'),mp.mpf('.5'),mp.mpf('.75'),1])
    bound=(d+ell)**2*mp.sqrt(energy)/2
    center=d**mp.mpf('2.5')*amp/ell**2*psi2(0)
    check('spike_center_'+str(digits)+'_'+str(n),abs(center-1)<mp.mpf('1e-40'),mp.nstr(center,20))
    check('spike_log_L1_bound_'+str(digits)+'_'+str(n),logmass<bound,mp.nstr(logmass/bound,25))
    check('spike_energy_decay_'+str(digits)+'_'+str(n),energy<100*d**3,mp.nstr(energy,25))
    rows.append({'index':n,'energy':mp.nstr(energy,digits),'logarithmic_L1':mp.nstr(logmass,digits),'bound':mp.nstr(bound,digits),'bound_ratio':mp.nstr(logmass/bound,25),'weighted_center':mp.nstr(center,25)})
   levels.append({'precision':digits,'values':rows})
 FIXTURES.append({'fixture':'Sobolev_spikes_not_actual_shape_solution','levels':levels,'actual_shape_solution':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 first-order mixed boundary response','input_count':len(ps)}),flush=True)
 try:
  provenance();coeff,raw=charged_identities();normal=normal_identities();psi2,norms=endpoint_controls()
  charged_fixtures(coeff,raw);normal_fixtures(normal);spike_fixtures(psi2,norms)
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 first-order mixed boundary response','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'exact_charged_mixed_endpoint_derived':good,'independent_mixed_trace_row_verified':good,
  'mixed_logarithmic_L1_limit_derived':good,'pure_conformal_lift_persistent_term_preserved':good,
  'pure_shape_finite_cut_limit_proved':False,'full_second_order_trace_with_first_order_shape_proved':False,
  'ordinary_mixed_cutoff_convergence_proved':False,'physical_action_averaging_adopted':False,
  'physical_fast_branch_excluded':False,'full_first_order_shape_response_included':False,'new_boundary_prescription_adopted':False,
  'rough_finite_amplitude_action_differentiability_proved':False,'full_shape_equation_solved':False,
  'static_source_compatibility_proved':False,'regulator_independence_for_moving_cuts_proved':False,
  'continuous_causal_stability_proved':False,'physical_constraint_mode_constructed':False,'actual_wormhole_instability_proved':False,
  'finite_counterterm_changed':False,'physical_curvature_norm_changed':False,'order_reduction_adopted':False,
  'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in [('checker',OWN[0]),('note',OWN[1]),('review',OWN[2])]:run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
