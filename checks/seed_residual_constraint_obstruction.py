"""Off-shell density Ward identity and the seed's full-kernel obstruction."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='seed-residual-constraint-obstruction'
PREVIOUS='continuous-causal-feedback-and-resonance'
OWN=['checks/seed_residual_constraint_obstruction.py','notes/'+STEM+'.md',
     'evidence/seed-residual-constraint-source-review.json']
ROWS=[];FIXTURES=[]
PROTOCOL={'precision_digits':[50,80],'identity_tolerance':'1e-28',
 'finite_difference_steps':['.0002','.0001'],'finite_difference_tolerance':'1e-28',
 'finite_difference_contract':'The orbit divergence is quadratic in the path parameter; symmetric secants are exact apart from roundoff.',
 'ward_points':[['.2','.4','.9'],['.35','.8','1.1']],
 'ward_gamma':'1.7','oscillatory_k':['1.3','4.7'],'asymptotic_tau':['20','40'],
 'oscillatory_correction_power':'.2','oscillatory_correction_amplitude':'.1',
 'roots_per_fixture':4,'charged_kernel_evaluated':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)}
 ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,value):
 value=S.simplify(S.sympify(value).doit());check(name,value==0,value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==112 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('review_inputs_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('scope',review['off_shell_full_auxiliary_kernel_zero_proved'] and not review['physical_fast_branch_excluded']
       and not review['charged_Lorentzian_instability_proved'] and not review['goal_complete'])

def symbolic_ward():
 t,y,th,ph,z=S.symbols('t y theta phi z',real=True);co=[t,y,th,ph]
 F=S.Function('F')(y);r=S.Function('r')(y)
 p=S.Function('p')(t,y);v=S.Function('v')(t,y);b=p+v;w0=4*p+2*v
 metric=[-F*S.exp(2*z*b),F*S.exp(2*z*b),r*r*S.exp(2*z*p),r*r*S.sin(th)**2*S.exp(2*z*p)]
 def gamma(a,b,c):
  return S.simplify(((S.diff(metric[a],co[b]) if a==c else 0)
   +(S.diff(metric[a],co[c]) if a==b else 0)
   -(S.diff(metric[b],co[a]) if b==c else 0))/(2*metric[a]))
 full=[[[gamma(a,b,c) for c in range(4)] for b in range(4)] for a in range(4)]
 g=[[[full[a][b][c].subs(z,0) for c in range(4)] for b in range(4)] for a in range(4)]
 dg=[[[S.diff(full[a][b][c],z).subs(z,0) for c in range(4)] for b in range(4)] for a in range(4)]
 Et,Ey,Eo=[S.Function(k)(y) for k in ['E_t','E_y','E_o']]
 E=S.diag(Et,Ey,Eo,Eo);Z=Et+Ey+2*Eo;O=Et+Ey
 e=S.Function('e')(t,y);j=S.Function('j')(t,y)
 Sm=S.diag(-e,e,0,0);Sm[1,0]=j;Sm[0,1]=-j
 dE=Sm-w0*E
 def div(M,G):
  return [sum(S.diff(M[a,b],co[a])+sum(G[a][a][c]*M[c,b]-G[c][a][b]*M[a,c] for c in range(4)) for a in range(4)) for b in range(4)]
 divE=div(E,g);divS=div(Sm,g);direct=div(dE,g)
 for bidx in range(2):
  direct[bidx]+=sum(dg[a][a][c]*E[c,bidx]-dg[c][a][bidx]*E[a,c] for a in range(4) for c in range(4))
  force=Z*S.diff(p,co[bidx])+O*S.diff(v,co[bidx])
  zero('direct_off_shell_Ward_'+str(bidx),direct[bidx]-(divS[bidx]-force-w0*divE[bidx]))
  zero('connection_trace_'+str(bidx),sum(dg[a][a][bidx] for a in range(4))-S.diff(w0,co[bidx]))
  for aidx in range(4):zero(f'diagonal_connection_{aidx}_{bidx}',dg[aidx][aidx][bidx]-S.diff(b if aidx<2 else p,co[bidx]))
 W=F*r*r
 zero('time_constraint_density',W*divS[0]+S.diff(W*e,t)-S.diff(W*j,y))
 zero('radial_constraint_density',W*divS[1]+S.diff(W*j,t)-S.diff(W*e,y))
 zero('optical_trace_row',sum(dE[a,a] for a in range(4))+w0*Z)
 zero('optical_orbit_row',dE[0,0]+dE[1,1]+w0*O)
 for aidx in range(3):zero('full_chart_density_gradient_'+str(aidx),S.diff(S.exp(z*w0)*(E[aidx,aidx]+z*dE[aidx,aidx]),z).subs(z,0)-Sm[aidx,aidx])
 zero('residual_Wronskian_elimination',(Z*p+O*v)*S.diff(v,y)-(Z*S.diff(p,y)+O*S.diff(v,y))*v-Z*(p*S.diff(v,y)-v*S.diff(p,y)))
 check('density_contact_is_nonzero',S.simplify(w0*Z)!=0)
 FIXTURES.append({'fixture':'symbolic_four_dimensional_Ward','arbitrary_background':True,'physical_solution':False})

def trace_and_boundary_algebra():
 u,x,mu,d,e0,q=S.symbols('u x mu d e0 q',positive=True)
 F=x/u**2+1-2*mu*u+(1+3*x+d)*u**2-2*d*u**3+d*u**4
 Dr=lambda f:-u*u*S.diff(f,u)
 R=-Dr(Dr(F))-4*u*Dr(F)-2*u*u*(F-1)
 Wc=-Dr(Dr(F))+2*u*Dr(F)+2*u*u*(1-F)
 jr=-Dr(Dr(F))/4+(F-1)*u*u/2
 BoxR=F*Dr(Dr(R))+(Dr(F)+2*u*F)*Dr(R)
 trace=S.expand(R+12*(x+S.Rational(11,6)*e0*x*x)+e0*(S.Rational(7,36)*Wc**2+S.Rational(22,3)*jr**2-S.Rational(11,72)*R**2-10*q*q*u**4+BoxR))
 zero('seed_scalar_curvature',R-(-12*x+4*d*u**5-6*d*u**6))
 for n in range(4):zero('seed_trace_lower_coefficient_'+str(n),trace.coeff(u,n))
 zero('magnetic_trace_residual_coefficient',-trace.coeff(u,4)-10*e0*q*q)
 P=S.cancel((-trace/u**4-10*e0*q*q)/u)
 check('trace_remainder_is_polynomial',P.is_polynomial(u),S.degree(P,u))
 zero('exact_trace_polynomial_reconstruction',trace+u**4*(10*e0*q*q+u*P))
 c,N=S.symbols('c N',positive=True)
 check('nonzero_trace_collar_margin',S.simplify(c-c*N/(2*(1+N))-c/2)>0)
 y,k=S.symbols('y k',positive=True);L=1-S.log(y);pp=y**S.Rational(3,2)*S.sin(k*S.log(y));vv=pp/L;ZZ=S.Function('Z')(y);OO=-ZZ*L
 zero('time_constraint_logarithmic_cancellation',ZZ*pp+OO*vv)
 zero('radial_constraint_logarithmic_failure',ZZ*S.diff(pp,y)+OO*S.diff(vv,y)+ZZ*pp/(y*L))
 zero('fast_trace_derivative',S.diff(pp,y)-S.sqrt(y)*(S.Rational(3,2)*S.sin(k*S.log(y))+k*S.cos(k*S.log(y))))
 zero('simple_zero_slope',S.diff(pp,y).subs(y,1)-k)
 cL,cR=S.symbols('cL cR')
 zero('piecewise_ratio_derivative_jump',(cR-cL)*S.diff(pp,y).subs(y,1)-k*(cR-cL))
 zero('H2_boundary_value_constant',S.integrate((y-u)**2,(u,0,y))-y**3/3)
 zero('H2_boundary_derivative_constant',S.integrate(1,(u,0,y))-y)
 zero('on_shell_residual_source_disappears',(ZZ*pp+OO*vv).subs(ZZ,0))
 FIXTURES.append({'fixture':'exact_seed_trace_polynomial','trace_coefficient':'10 e0 q_m^2',
  'P':str(P),'nonzero_flux_required':True,'physical_solution':False})

def finite_metric_ward():
 for case,point in enumerate(PROTOCOL['ward_points']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    pt=list(map(mp.mpf,point))+[mp.mpf('.1')];gam=mp.mpf(PROTOCOL['ward_gamma'])
    def at(x):
     t,y,theta,phi=x;F=(1+y)**-2;r=2+1/(1+y);Fr=2*(r-2)
     Et=(r*Fr+F-1)/r**2-mp.mpf('.39');Eo=1+Fr/r-mp.mpf('.39')
     p=mp.exp(gam*t)*y*y*mp.exp(-y)/7;v=mp.exp(gam*t)*y**3*mp.exp(-2*y)/11
     ee=mp.exp(gam*t)*y/(1+y*y);jj=mp.exp(gam*t)*(1+y)/5
     return F,r,Et,Eo,p,v,ee,jj
    def metric(x,z):
     F,r,Et,Eo,p,v,ee,jj=at(x);eb=mp.exp(2*z*(p+v));ep=mp.exp(2*z*p)
     return [-F*eb,F*eb,r*r*ep,r*r*mp.sin(x[2])**2*ep]
    def tensor(x,z):
     F,r,Et,Eo,p,v,ee,jj=at(x);E=mp.diag([Et,Et,Eo,Eo]);Sm=mp.diag([-ee,ee,0,0]);Sm[1,0]=jj;Sm[0,1]=-jj
     return E+z*(Sm-(4*p+2*v)*E)
    def derivative(fn,x,i):
     def one(a):
      xx=list(x);xx[i]=a;return fn(xx)
     return mp.diff(one,x[i])
    def divergence(z):
     gm=metric(pt,z);dg=[[derivative(lambda xx:metric(xx,z)[j],pt,i) for j in range(4)] for i in range(4)]
     def G(a,b,c):return ((dg[b][a] if a==c else 0)+(dg[c][a] if a==b else 0)-(dg[a][b] if b==c else 0))/(2*gm[a])
     M=tensor(pt,z)
     return mp.matrix([sum(derivative(lambda xx:tensor(xx,z)[a,b],pt,a)+sum(G(a,a,c)*M[c,b]-G(c,a,b)*M[a,c] for c in range(4)) for a in range(4)) for b in range(2)])
    F,r,Et,Eo,p,v,ee,jj=at(pt);Z=2*Et+2*Eo;O=2*Et;W=F*r*r
    Ecal=lambda xx:at(xx)[0]*at(xx)[1]**2*at(xx)[6]
    Jcal=lambda xx:at(xx)[0]*at(xx)[1]**2*at(xx)[7]
    divS=mp.matrix([(-derivative(Ecal,pt,0)+derivative(Jcal,pt,1))/W,(-derivative(Jcal,pt,0)+derivative(Ecal,pt,1))/W])
    force=mp.matrix([Z*derivative(lambda xx:at(xx)[4],pt,b)+O*derivative(lambda xx:at(xx)[5],pt,b) for b in range(2)])
    pred=divS-force
    check(f'background_conservation_{case}_{digits}',mp.norm(divergence(0))<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(mp.norm(divergence(0)),digits))
    check(f'residual_contact_detected_{case}_{digits}',mp.norm(force)>mp.mpf('.001'),mp.nstr(mp.norm(force),digits))
    errors=[];records=[]
    for hs in PROTOCOL['finite_difference_steps']:
     h=mp.mpf(hs);fd=(divergence(h)-divergence(-h))/(2*h);err=mp.norm(fd-pred)/(1+mp.norm(pred));errors.append(err)
     check(f'independent_metric_Ward_FD_{case}_{digits}_{hs}',err<mp.mpf(PROTOCOL['finite_difference_tolerance']),mp.nstr(err,digits))
     records.append({'step':hs,'scaled_error':mp.nstr(err,digits)})
    # All contracted orbit Christoffels are affine in z and the tensor is
    # affine in z. A central secant removes the quadratic term exactly.
    check(f'metric_Ward_exact_secants_{case}_{digits}',max(errors)<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(max(errors),digits))
    levels.append({'precision':digits,'predicted_variation':[mp.nstr(a,digits) for a in pred],'finite_differences':records})
  with mp.workdps(80):
   aa=mp.matrix(list(map(mp.mpf,levels[0]['predicted_variation'])));bb=mp.matrix(list(map(mp.mpf,levels[1]['predicted_variation'])))
   err=mp.norm(aa-bb)/(1+mp.norm(bb));check('metric_Ward_precision_'+str(case),err<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(err,80))
  FIXTURES.append({'fixture':'independent_finite_metric_divergence','point':point,'levels':levels,'charged_operator':False})

def oscillatory_traces():
 for case,ks in enumerate(PROTOCOL['oscillatory_k']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    k=mp.mpf(ks);delta=mp.mpf('.2');amp=mp.mpf('.1')
    p=lambda y:y**mp.mpf('1.5')*(mp.sin(k*mp.log(y))+amp*y**delta*mp.cos(2*k*mp.log(y)))
    v=lambda y:y**mp.mpf('1.7')/(1+y)
    P=lambda tau:p(mp.exp(-tau))/mp.exp(-mp.mpf('1.5')*tau)
    errors=[]
    for ts in PROTOCOL['asymptotic_tau']:
     tau=mp.mpf(ts);y=mp.exp(-tau)
     ev=abs(P(tau)+mp.sin(k*tau));ed=abs(mp.diff(P,tau)+k*mp.cos(k*tau));errors.append(max(ev,ed))
     check(f'differentiated_trace_remainder_{case}_{digits}_{ts}',max(ev,ed)<=amp*(1+delta+2*k)*y**delta,mp.nstr(max(ev,ed),digits))
    ratio=errors[1]/errors[0];check(f'trace_remainder_decay_{case}_{digits}',ratio<mp.mpf('.1'),mp.nstr(ratio,digits))
    roots=[];start=int(mp.floor(20*k/mp.pi))+1
    for n in range(start,start+PROTOCOL['roots_per_fixture']):
     target=n*mp.pi/k;root=mp.findroot(P,(target-mp.mpf('.05')/k,target+mp.mpf('.05')/k));y=mp.exp(-root)
     residual=abs(P(root));slope=abs(mp.diff(P,root));wr=p(y)*mp.diff(v,y)-v(y)*mp.diff(p,y)
     scaled=abs(wr)/y**mp.mpf('2.2')
     check(f'simple_fast_zero_{case}_{digits}_{n}',residual<mp.mpf(PROTOCOL['identity_tolerance']) and slope>k/2 and abs(root-target)<mp.mpf('.05')/k,mp.nstr(slope,digits))
     check(f'nonproportional_shape_Ward_failure_{case}_{digits}_{n}',scaled>k/2,mp.nstr(scaled,digits))
     roots.append({'tau':mp.nstr(root,digits),'slope':mp.nstr(slope,digits),'scaled_Wronskian':mp.nstr(scaled,digits)})
    levels.append({'precision':digits,'roots':roots,'trace_error_ratio':mp.nstr(ratio,digits)})
  with mp.workdps(80):
   err=max(abs(mp.mpf(a['tau'])-mp.mpf(b['tau']))/(1+abs(mp.mpf(b['tau']))) for a,b in zip(levels[0]['roots'],levels[1]['roots']))
   check('oscillatory_precision_'+str(case),err<mp.mpf(PROTOCOL['identity_tolerance']),mp.nstr(err,80))
  FIXTURES.append({'fixture':'oscillatory_C1_trace_and_Ward_failure','k':ks,'levels':levels,'charged_mode':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 seed residual constraint obstruction','input_count':len(ps)}),flush=True)
 try:provenance();symbolic_ward();trace_and_boundary_algebra();finite_metric_ward();oscillatory_traces()
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 seed residual constraint obstruction','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'off_shell_density_Ward_identity_proved':good,'seed_outer_trace_nonzero_proved':good,
  'differentiated_full_metric_trace_proved':good,'off_shell_full_auxiliary_kernel_zero_proved':good,
  'reduced_resonances_full_seed_Hessian_lift_excluded':good,
  'physical_fast_branch_excluded':False,'full_tensor_inverse_proved':False,'general_gauge_coverage_proved':False,
  'actual_numeric_threshold_certified':False,'physical_boundary_state_admissibility_proved':False,
  'actual_charged_two_point_kernel_evaluated':False,'finite_counterterm_changed':False,
  'new_boundary_prescription_adopted':False,'physical_source_projection_adopted':False,'order_reduction_adopted':False,
  'physical_exterior_state_constructed':False,'continuous_causal_stability_proved':False,'charged_Lorentzian_instability_proved':False,
  'formal_verification':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]};data['runs'].append(run)
 dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
