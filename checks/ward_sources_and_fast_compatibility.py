"""Original-source adjoints and exact Ward completion diagnostics."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='ward-sources-and-fast-compatibility'
PREVIOUS='combined-shape-cutoff-obstruction'
OWN=['checks/ward_sources_and_fast_compatibility.py','notes/'+STEM+'.md',
     'evidence/ward-source-compatibility-source-review.json']
ROWS=[]; FIXTURES=[]
PROTOCOL={'precision_digits':[50,80], 'precision_tolerance':'1e-32',
 'thermal_harmonics':[2,4], 'kappa':'1', 'source_interval':['-2','-1'],
 'integration_interval':['-3','-.5'], 'rk4_steps':[64,128,256],
 'rk4_scaled_error_limit':'1e-3', 'rk4_refinement_ratio_limit':'.12',
 'model_geometry':'r=1/(1-exp(2x)), N^2=exp(2x)/(1-exp(2x))^2',
 'model_geometry_is_actual_charged_seed':False,
 'numeric_bump':'4096*(x+2)^6*(-1-x)^6 on [-2,-1], zero elsewhere; C5 fixture',
 'theorem_sources':'arbitrary smooth compact original pairs, completed by the exact Green formula',
 'actual_Dirac_state_source_evaluated':False, 'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)}
 ROWS.append(row); print(json.dumps(row),flush=True)
 if not ok: raise AssertionError(name)

def zero(name,value):
 value=S.simplify(S.sympify(value).doit())
 check(name,value==0,value)

def prior():
 return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]

def paths():
 return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior()
 check('predecessor_passed',p['status']=='passed' and len(p['checks'])==133 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('retained_review_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('scope',review['original_source_pullback_derived'] and review['Ward_conditions_alone_insufficient']
       and not review['actual_Dirac_source_incompatibility_proved'])

def geometric_identities():
 t,x,th,ph=S.symbols('t x theta phi',real=True); coords=[t,x,th,ph]
 gamma=S.symbols('gamma',positive=True)
 N=S.Function('N')(x); r=S.Function('r')(x); W=N*N*r*r
 E,J,tau,o=[S.Function(n)(x) for n in ['E','J','tau','o']]
 g=S.diag(-N*N,N*N,r*r,r*r*S.sin(th)**2)
 inv=g.inv()
 christoffel=[[[S.simplify(sum(inv[a,d]*(S.diff(g[d,c],coords[b])+S.diff(g[d,b],coords[c])-S.diff(g[b,c],coords[d])) for d in range(4))/2)
                for c in range(4)] for b in range(4)] for a in range(4)]
 T=S.exp(gamma*t)*S.Matrix([[o-E/W,-J/W,0,0],[J/W,E/W,0,0],
                            [0,0,(tau-o)/2,0],[0,0,0,(tau-o)/2]])
 cov=g*T
 check('covariant_stress_symmetric',cov==cov.T)
 F=gamma*W*o; G1=-W*((S.diff(N,x)/N-S.diff(r,x)/r)*o+S.diff(r,x)/r*tau)
 targets=[(F-gamma*E+S.diff(J,x))/W,(-gamma*J+S.diff(E,x)+G1)/W,0,0]
 for b in range(4):
  div=sum(S.diff(T[a,b],coords[a]) for a in range(4))
  div+=sum(christoffel[a][a][c]*T[c,b]-christoffel[c][a][b]*T[a,c] for a in range(4) for c in range(4))
  zero('four_dimensional_divergence_'+str(b),S.exp(-gamma*t)*div-targets[b])
 zero('original_trace',S.exp(-gamma*t)*S.trace(T)-tau)
 zero('original_orbit_trace',S.exp(-gamma*t)*(T[0,0]+T[1,1])-o)
 p,v=S.symbols('p v'); h=S.diag(2*(p+v),2*(p+v),2*p,2*p)
 zero('source_work_contraction',S.exp(-gamma*t)*S.trace(T*h)/2-(tau*p+o*v))
 eps,GN,ff,gg,rr,ww=S.symbols('eps GN f g r W',nonzero=True)
 zero('trace_source_normalization',4*S.pi*ww*eps*rr*ff/(8*S.pi*GN*ww)-eps*rr*ff/(2*GN))
 zero('orbit_source_normalization',4*S.pi*ww*rr*gg/(4*S.pi*ww)-rr*gg)
 cp,cm,fp,fm=S.symbols('cp cm fp fm',real=True)
 Ec=(cp+cm)/2; Jc=(cp-cm)/2
 # Transport derivatives cp'=gamma cp-fp, cm'=fm-gamma cm.
 dx=lambda f:S.diff(f,cp)*(gamma*cp-fp)+S.diff(f,cm)*(fm-gamma*cm)
 zero('characteristic_first_equation',gamma*Ec-dx(Jc)-(fp+fm)/2)
 zero('characteristic_second_equation',gamma*Jc-dx(Ec)-(fp-fm)/2)
 zero('energy_pair_norm',Ec*Ec+Jc*Jc-(cp*cp+cm*cm)/2)
 zero('boundary_flux',Ec*Jc-(cp*cp-cm*cm)/4)
 kap,V,amp,rad=S.symbols('kappa V amplitude radius',positive=True)
 for n in [1,2,3,5]:
  coefficient=amp*V**n/(2*rad**2)/(kap*V)**2
  zero('Kruskal_coefficient_'+str(n),coefficient-amp*V**(n-2)/(2*kap**2*rad**2))
  power=n-2
  check('Kruskal_regular_'+str(n),((power>=0)==(n>=2)),power)
 check('first_harmonic_pole',S.limit(V**(-1),V,0,dir='+')==S.oo)

def weighted_source_identities():
 I=S.I; eps=S.Rational(2,7); GN=S.Rational(3,5); alpha=eps/(2*GN)
 H=S.Matrix([[3,1,I],[1,4,0],[-I,0,5]])
 A=S.Matrix([[4,1+I],[1-I,3]])
 Q=S.Matrix([[1,I],[2,-1],[I,3]])
 D=S.Matrix([[2,1+I],[I,3],[1,2-I]])
 j=S.Matrix([[1,0],[0,1],[1,-1]])
 R=S.Matrix([[2,I,1],[0,3,1-I],[I,0,4]])
 U=S.Matrix([1+I,2-I,3]); w=S.Matrix([2+I,-1+2*I])
 f=S.Matrix([1-I,2+I,-3]); g=S.Matrix([2,1+3*I,4-I]); z=S.Matrix([2-I,1,3+I])
 inner=lambda a,b:(a.conjugate().T*H*b)[0]
 form=lambda a,b:(a.conjugate().T*A*b)[0]
 check('Hilbert_Gram_positive',all(H[:i,:i].det()>0 for i in range(1,4)))
 check('energy_Gram_positive',all(A[:i,:i].det()>0 for i in range(1,3)))
 check('shape_embedding_injective',j.rank()==2)
 beta=alpha*Q.conjugate().T*H*f+j.conjugate().T*H*g
 well=A.inv()*beta; nu=A.inv()*D.conjugate().T*H*z/eps
 mu=z+alpha*Q*nu
 zero('Riesz_shape_source',form(w,well)-(alpha*inner(Q*w,f)+inner(j*w,g)))
 zero('Riesz_mixed_adjoint',form(w,nu)-inner(D*w,z)/eps)
 zero('full_original_source_work',alpha*inner(R*U+Q*w,f)+inner(j*w,g)
      -(alpha*inner(R*U,f)+form(w,well)))
 zero('full_original_source_compatibility',inner(z,f+D*well/eps)-(inner(mu,f)+inner(j*nu,g)))
 K=-H.inv()*R.conjugate().T*H*D/(2*GN)
 Kdag=A.inv()*K.conjugate().T*H
 zero('cross_form_adjoint',form(w,Kdag*U)-inner(K*w,U))
 wtotal=well-Kdag*U
 zero('shape_elimination',form(w,wtotal)+inner(K*w,U)-form(w,well))
 missed=S.simplify(inner(z,D*A.inv()*(alpha*Q.conjugate().T*H*f)/eps))
 check('omitted_induced_shape_source_changes_overlap',missed!=0,missed)
 zero('tracefree_nonzero_adjoint_overlap',inner(z,D*A.inv()*j.conjugate().T*H*j*nu/eps)-inner(j*nu,j*nu))
 normnu=S.simplify(inner(j*nu,j*nu))
 check('tracefree_incompatible_control',normnu>0,normnu)
 z0=(D.conjugate().T*H).nullspace()[0]
 nu0=(A.inv()*D.conjugate().T*H*z0/eps).applyfunc(S.simplify)
 check('zero_shape_adjoint_control',nu0==S.zeros(2,1))
 zero('zero_shape_adjoint_traceless_compatibility',inner(z0,D*A.inv()*j.conjugate().T*H*g/eps))
 zero('zero_shape_adjoint_traceful_overlap',inner(z0,z0+D*A.inv()*alpha*Q.conjugate().T*H*z0/eps)-inner(z0,z0))
 check('zero_shape_adjoint_traceful_incompatibility',S.simplify(inner(z0,z0))>0)
 FIXTURES.append({'fixture':'complex_weighted_source_pullback',
  'omitted_Q_contribution':str(missed),'traceless_nonzero_nu_overlap':str(normnu),
  'actual_cokernel_profile_used':False,'finite_matrix_index_claimed':False})

def transport_fixtures():
 for harmonic in PROTOCOL['thermal_harmonics']:
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    gamma=mp.mpf(harmonic)
    def data(x):
     if x<=-2 or x>=-1:return mp.mpf(0),mp.mpf(0)
     bump=4096*(x+2)**6*(-1-x)**6
     return bump*(1+x/7),bump*(x+mp.mpf('1.4'))
    def geometry(x):
     u=mp.exp(2*x); rad=1/(1-u); N2=u*rad**2
     return N2*rad**2,1+2*u/(1-u),2*u/(1-u)
    def forcing(x):
     tau,o=data(x)
     if tau==0 and o==0:return mp.mpf(0),mp.mpf(0)
     W,nd,ld=geometry(x)
     return gamma*W*o,-W*((nd-ld)*o+ld*tau)
    def integral(fun,a,b):
     if a>=b:return mp.mpf(0)
     return mp.quad(fun,[a,(a+b)/2,b])
    Cminus0=integral(lambda s:mp.exp(gamma*s)*(forcing(s)[0]-forcing(s)[1]),mp.mpf(-2),mp.mpf(-1))
    def green(x):
     cm=integral(lambda s:mp.exp(-gamma*(x-s))*(forcing(s)[0]-forcing(s)[1]),mp.mpf(-2),min(x,mp.mpf(-1)))
     cp=mp.exp(gamma*x)*Cminus0+integral(lambda s:mp.exp(gamma*(x-s))*(forcing(s)[0]+forcing(s)[1]),max(x,mp.mpf(-2)),mp.mpf(-1))
     return (cp+cm)/2,(cp-cm)/2
    E0,J0=green(mp.mpf(0)); check(f'outer_zero_flux_{harmonic}_{digits}',J0==0)
    Eh,Jh=green(mp.mpf(-3)); check(f'horizon_characteristic_{harmonic}_{digits}',Eh==Jh)
    Eh2,Jh2=green(mp.mpf(-4))
    er=abs(Eh*mp.exp(3*gamma)-Eh2*mp.exp(4*gamma))/(1+abs(Eh*mp.exp(3*gamma)))
    check(f'horizon_amplitude_{harmonic}_{digits}',er<mp.mpf('1e-32'),mp.nstr(er,digits))
    # Integrate the radial covariant conservation equations directly.
    # This does not call the characteristic forcing routine.
    def rhs(x,y):
     tau,o=data(x); W,nd,ld=geometry(x)
     return mp.matrix([gamma*y[1]+W*((nd-ld)*o+ld*tau),gamma*y[0]-gamma*W*o])
    target=mp.matrix(green(mp.mpf('-.5'))); errors=[]
    for count in PROTOCOL['rk4_steps']:
     x=mp.mpf(-3); h=mp.mpf('2.5')/count; val=mp.matrix([Eh,Jh])
     for step in range(count):
      x=mp.mpf(-3)+step*h
      k1=rhs(x,val); k2=rhs(x+h/2,val+h*k1/2)
      k3=rhs(x+h/2,val+h*k2/2); k4=rhs(x+h,val+h*k3)
      val+=h*(k1+2*k2+2*k3+k4)/6
     error=mp.norm(val-target)/(1+mp.norm(target)); errors.append(error)
    ratios=[errors[i+1]/errors[i] for i in range(len(errors)-1)]
    check(f'Ward_ODE_refinement_{harmonic}_{digits}',max(ratios)<mp.mpf(PROTOCOL['rk4_refinement_ratio_limit']),mp.nstr(max(ratios),digits))
    check(f'Ward_ODE_final_error_{harmonic}_{digits}',errors[-1]<mp.mpf(PROTOCOL['rk4_scaled_error_limit']),mp.nstr(errors[-1],digits))
    levels.append({'precision':digits,'endpoint_E':mp.nstr(target[0],digits),'endpoint_J':mp.nstr(target[1],digits),
                   'outer_E':mp.nstr(E0,digits),'horizon_amplitude':mp.nstr(2*Eh*mp.exp(3*gamma),digits),
                   'rk4_scaled_errors':[mp.nstr(a,digits) for a in errors],
                   'rk4_refinement_ratios':[mp.nstr(a,digits) for a in ratios]})
  with mp.workdps(80):
   error=max(abs(mp.mpf(levels[0][k])-mp.mpf(levels[1][k]))/(1+abs(mp.mpf(levels[1][k])))
             for k in ['endpoint_E','endpoint_J','outer_E','horizon_amplitude'])
   check('transport_precision_'+str(harmonic),error<mp.mpf(PROTOCOL['precision_tolerance']),mp.nstr(error,80))
  FIXTURES.append({'fixture':'conserved_stress_Green_vs_covariant_ODE','harmonic':harmonic,'levels':levels,
                   'model_geometry_is_actual_charged_seed':False,'Dirac_state_realizability_claimed':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n'); sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths(); hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
 before=hashes(); error=None
 print(json.dumps({'event':'start','story':'S168 Ward sources and fast compatibility','input_count':len(ps)}),flush=True)
 try:
  provenance(); geometric_identities(); weighted_source_identities(); transport_fixtures()
 except Exception:
  error=traceback.format_exc(); print(error,file=sys.stderr,flush=True)
 after=hashes(); good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 Ward sources and fast compatibility','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
  'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
  'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'original_source_pullback_derived':good,'Ward_conditions_alone_insufficient':good,
  'actual_Dirac_source_incompatibility_proved':False,'actual_adjoint_profiles_evaluated':False,
  'actual_charged_two_point_kernel_evaluated':False,'first_harmonic_arbitrary_completion_claimed':False,
  'physical_fast_branch_excluded':False,'finite_counterterm_changed':False,'new_boundary_prescription_adopted':False,
  'order_reduction_adopted':False,'physical_exterior_state_constructed':False,'continuous_causal_stability_proved':False,
  'formal_verification':False,'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json'); data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run); dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
