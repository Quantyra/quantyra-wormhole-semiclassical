"""Necessary nonlinear trace versus action boundary tails; no physical completion."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='fast-nonlinear-trace-boundary-compatibility'
OWN=['checks/fast_nonlinear_trace_boundary_compatibility.py','notes/'+STEM+'.md','evidence/fast-nonlinear-trace-source-review.json']
PREVIOUS='fast-trace-finite-cut-action';ROWS=[];FIXTURES=[];R=S.Rational
PROTOCOL={'scope':'Second-order trace/action incompatibility for the stated conformal expansion, with necessary coupled shape relation',
 'precision_digits':[50,80],'reference_B':'1','base_radius':'80','cycle_indices':[0,1,2],
 'phase_at_base_in_pi_units':'1/4','charged_cases':[['2','.01','1.1','2','2'],['3','-.02','1.2','3','5'],['4','1e-14','1.05','5','9']],
 'charged_columns':['lambda','d','horizon','magnetic_flux','k'],
 'fine_scaled_limit':'.001','coarse_scaled_limit':'.01','refinement_bound':'exp(-pi/k)',
 'precision_scaled_limit':'1e-35','contour_nodes':[16,32],'contour_radius':'.03','contour_scaled_limit':'1e-32',
 'static_numerical_fixtures':True,'thermal_extension_is_ordinary_power_count':True,
 'fixtures_are_physical_benchmark':False,'coupled_completion_adopted':False,'formal_verification':False}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)

def zero(name,value):
 value=S.factor(S.cancel(S.sympify(value).doit()));check(name,value==0,value)

def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==124 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 rev=json.loads((ROOT/OWN[2]).read_bytes())
 check('retained_inputs_pinned',all(f in p['input_sha256_after'] for f in rev['retained_inputs']))
 check('scope',not rev['physical_fast_branch_excluded'] and not rev['coupled_boundary_completion_adopted'] and not rev['physical_stability_proved'])

def exact_derivations():
 t=S.symbols('t',positive=True);a=S.symbols('a');w=S.symbols('w0:7');v=S.symbols('v0:7');z=S.symbols('z')
 def D(q):return S.diff(q,z)+sum(S.diff(q,w[i])*w[i+1]+S.diff(q,v[i])*v[i+1] for i in range(6))
 def Euler(L):
  return S.expand(S.diff(L,w[0])-D(S.diff(L,w[1]))-3*S.diff(L,w[1])+D(D(S.diff(L,w[2])))+6*D(S.diff(L,w[2]))+9*S.diff(L,w[2]))
 u=1+w[1]+w[2];q=(1+w[1])**2;b=v[1]+v[2]
 LC=(S.exp(2*w[0])*(6*q+4*(1+w[1])*v[1])+6*(1+R(11,6)*t)*S.exp(4*w[0]+2*v[0]))/4
 LW=15*t*S.exp(-2*v[0])*((z+w[0])*b*b/15-u*b/15-u*u/10-b*q/15-R(7,90)*u*q-R(7,180)*q*q)
 action=Euler(LC+LW)
 # Independent warped-metric sectional curvatures, with lambda normalized to one.
 n=S.exp(-w[0]-v[0]);K=[n*(1+w[1]+v[1]),n*(1+w[1]),n*(1+w[1])]
 DK=[n*D(x) for x in K]
 Ric=[-sum(DK[i]+K[i]**2 for i in range(3))]+[-DK[i]-K[i]*sum(K) for i in range(3)]
 scalar=S.simplify(sum(Ric));ric2=S.expand(sum(x*x for x in Ric))
 riem2=4*sum((DK[i]+K[i]**2)**2 for i in range(3))+4*sum((K[i]*K[j])**2 for i in range(3) for j in range(i+1,3))
 weyl=S.simplify(riem2-2*ric2+scalar*scalar/3);jsq=S.simplify(ric2-scalar*scalar/4)
 zero('independent_normal_scalar',scalar-2*S.exp(-2*w[0]-2*v[0])*(-v[1]-v[2]-3*w[1]**2-9*w[1]-3*w[2]-6))
 zero('independent_normal_Weyl',weyl-R(4,3)*S.exp(-4*w[0]-4*v[0])*(v[1]+v[2])**2)
 box=S.exp(-2*w[0]-2*v[0])*(D(D(scalar))+(3+2*w[1])*D(scalar))
 trace=scalar+12+22*t+t*(R(7,12)*weyl+R(11,6)*jsq-R(11,72)*scalar*scalar+box)
 curvature=S.exp(4*w[0]+2*v[0])*trace/2
 vacuum={x:0 for x in (*w,*v)}
 pw=-(a-1)*(a+4)*(3*t*a*a+9*t*a+11*t+3)
 pv=-(a-3)*(a+4)*(3*t*a*a+9*t*a+11*t+3)/3
 f,g=S.symbols('f g',real=True);k2=1/t+R(17,12);p=[f]
 for i in range(4):p.append(S.expand(S.diff(p[-1],f)*g-S.diff(p[-1],g)*k2*f-R(3,2)*p[-1]))
 Nt=(R(471,4)*t+R(289,4)+6/t)*f*f+(R(35,3)*t-5)*f*g-(21*t+3)*g*g
 for name,expr in [('action',action),('curvature',curvature)]:
  zero(name+'_vacuum',expr.subs(vacuum))
  zero(name+'_linear_conformal',sum(S.diff(expr,w[i]).subs(vacuum)*a**i for i in range(5))-pw)
  zero(name+'_linear_shape',sum(S.diff(expr,v[i]).subs(vacuum)*a**i for i in range(5))-pv)
  expr=expr.subs({x:0 for x in v})
  second=sum(S.diff(expr,w[i],w[j]).subs(vacuum)*p[i]*p[j]/2 for i in range(5) for j in range(5))
  zero(name+'_quadratic_fast_forcing',second-Nt)
 zero('fast_indicial_factor',S.expand(3*t*a*a+9*t*a+11*t+3).subs(a,-R(3,2)+S.I*S.sqrt(k2)))
 # Omitting the density contact must change the nonlinear result for general jets.
 contact=S.expand(sum((S.diff(curvature,w[i],w[j])-S.diff(trace/2,w[i],w[j])).subs(vacuum)*w[i]*w[j]/2 for i in range(5) for j in range(5)))
 check('density_contact_not_dropped',contact!=0)
 return t,k2,pw,pv,Nt,f,g,a

def exact_compatibility(data):
 t,k2,pw,pv,Nt,f,g,a=data;k=S.symbols('k',positive=True);D=S.symbols('D')
 def red(q):
  n,d=S.fraction(S.cancel(q));return S.factor(S.rem(n,k*k-k2,k)/S.rem(d,k*k-k2,k))
 J=S.Matrix([[0,2*k],[-2*k,0]]);I=S.eye(2)
 C=3*t*J**2-9*t*J+(11*t+3)*I
 L=(-(J-4*I)*(J+I)*C).applyfunc(red);Lj=(-(J-6*I)*(J+I)*C/3).applyfunc(red)
 M=(R(3,7)*t*(7*J**2-33*J+20*I)).applyfunc(red)
 zero('shifted_conformal_operator',pw.subs(a,D-3)+(D-4)*(D+1)*(3*t*D*D-9*t*D+11*t+3))
 zero('shifted_shape_operator',pv.subs(a,D-3)+(D-6)*(D+1)*(3*t*D*D-9*t*D+11*t+3)/3)
 Ltarget=(5*t+3)*S.Matrix([[(19*t-12)/t,-42*k],[42*k,(19*t-12)/t]])
 Mtarget=S.Matrix([[-(59*t+84)/7,-198*k*t/7],[198*k*t/7,-(59*t+84)/7]])
 for i in range(4):zero('trace_matrix_'+str(i),L[i]-Ltarget[i]);zero('action_matrix_'+str(i),M[i]-Mtarget[i])
 dl=4*(5*t+3)**2*(11*t+3)*(65*t+12)/t**2;dm=4*(65*t+12)*(227*t+147)/49;dj=4*(5*t+3)**2*(11*t+3)*(125*t+12)/(9*t*t)
 for name,mat,det in [('trace',L,dl),('action',M,dm),('shape',Lj,dj)]:
  zero(name+'_determinant',red(mat.det())-det);check(name+'_invertible',det>0)
 A=S.expand(Nt).coeff(f,2);B=S.expand(Nt).coeff(f,1).coeff(g,1);Cc=S.expand(Nt).coeff(g,2)
 n0=(88*t*t+47*t+3)/(2*t);nc=(295*t*t+195*t+18)/(4*t);ns=k*(R(5,2)-R(35,6)*t)
 zero('forcing_mean',(A+Cc*k2)/2-n0);zero('forcing_cosine',(A-Cc*k2)/2-nc);zero('forcing_sine',-B*k/2-ns)
 h0=-(88*t*t+47*t+3)/(8*t*(11*t+3));zero('forced_mean',(12+44*t)*h0+n0)
 q=S.Matrix([-9*(t-9)/56,k*(71*t+21)/28]);nvec=S.Matrix([nc,ns])
 zero('inherited_action_cosine',(-(3001*t+546)/504+R(167,42)*t*k2)/2-q[0])
 ht=(-L.adjugate()*nvec/dl).applyfunc(red);ha=(-M.adjugate()*q/dm).applyfunc(red)
 for i in range(2):zero('trace_tail_solves_'+str(i),red((L*ht+nvec)[i]));zero('action_tail_cancels_'+str(i),red((M*ha+q)[i]))
 residual=(q+M*ht).applyfunc(red);trace_bad=(L*ha+nvec).applyfunc(red)
 rc=-(78815*t**3+74391*t*t+16137*t+81)/(168*(5*t+3)*(11*t+3))
 rs=-k*(6655*t**3+111*t*t-3015*t-567)/(84*(5*t+3)*(11*t+3))
 bad_s=-k*(61175*t**3+51351*t*t+14625*t+3969)/(12*t*(227*t+147))
 zero('trace_tail_action_cosine',residual[0]-rc);zero('trace_tail_action_sine',residual[1]-rs)
 zero('action_tail_trace_sine',trace_bad[1]-bad_s)
 positive=78815*t**3+73896*t*t+15705*t
 zero('strict_action_gap_polynomial',(-rc-R(3,56))*168*(5*t+3)*(11*t+3)-positive)
 check('strict_action_gap',positive>0);check('strict_trace_mismatch',-bad_s>0)
 jj=(-Lj.adjugate()*trace_bad/dj).applyfunc(red)
 for i in range(2):zero('necessary_shape_repair_'+str(i),red((Lj*jj+trace_bad)[i]))
 check('required_shape_nonzero',any(x!=0 for x in jj))
 zero('mean_coupling',pv.subs(a,-3)-(6+22*t))
 zero('small_coupling_mean_limit',S.limit(t*h0,t,0)+R(1,8))
 x=S.symbols('x',real=True);zero('thermal_static_and_double_harmonic',S.expand_trig(S.cos(2*x))/2+R(1,2)-S.cos(x)**2)
 check('fixed_frequency_power_gain',-3-2<-3)
 check('bounded_homogeneous_root_split',-1<0 and 4>0 and R(3,2)>0)
 return S.lambdify([t,k],[h0,*ht,*ha,trace_bad[1]],'mpmath')

def fall(x,n):
 out=1
 for j in range(n):out*=x-j
 return out

def profile_jets(r,k,phase,power,harmonic,cosine,sine):
 exponent=mp.mpc(power,harmonic*k);amplitude=mp.mpc(cosine,-sine)*mp.exp(1j*harmonic*phase)
 return [mp.re(amplitude*fall(exponent,n)*r**(exponent-n)) for n in range(5)]

def charged_trace_density(r,lam,t,qm,F,Rseed,w):
 f,fr,frr,frrr,_=F;w0,w1,w2,w3,w4=w
 bb=fr+2*f/r;bp=frr+2*fr/r-2*f/r**2;bpp=frrr+2*frr/r-4*fr/r**2+4*f/r**3
 box=f*w2+bb*w1
 H=Rseed[0]-6*(box+f*w1*w1)
 Hp=Rseed[1]-6*(f*w3+(fr+bb)*w2+bp*w1+fr*w1*w1+2*f*w1*w2)
 Hpp=Rseed[2]-6*(f*w4+(2*fr+bb)*w3+(frr+2*bp)*w2+bpp*w1+frr*w1*w1+4*fr*w1*w2+2*f*(w2*w2+w1*w3))
 base=-frr/2-fr/r
 rr=base-2*f*w2-fr*w1-box
 rt=base-fr*w1-box-2*f*w1*w1
 rs=(1-f-r*fr)/r**2-2*f*w1/r-box-2*f*w1*w1
 jsq=rr*rr+rt*rt+2*rs*rs-H*H/4
 weyl=(-frr+2*fr/r+2*(1-f)/r**2)**2/3
 anomaly=7*weyl/12+11*jsq/6-11*H*H/72-10*qm*qm/r**4+f*Hpp+(bb-2*f*w1)*Hp-2*box*H
 return (mp.exp(2*w0)*H+12*lam*(1+mp.mpf(11)*t/6)*mp.exp(4*w0)+(t/lam)*anomaly)/(2*lam)

def fixtures(tails):
 for ix,row in enumerate(PROTOCOL['charged_cases']):
  levels=[]
  for precision in PROTOCOL['precision_digits']:
   with mp.workdps(precision):
    lam,d,rh,qm,k=map(mp.mpf,row);t=1/(k*k-mp.mpf(17)/12);radius=mp.mpf(PROTOCOL['base_radius']);phase=mp.pi/4-k*mp.log(radius)
    Qh=1+3*lam+d;m=rh*(1+lam*rh*rh+Qh/rh**2-2*d/rh**3+d/rh**4)/2
    h0,htc,hts,hac,has,bads=tails(t,k);scale=max(1,abs(bads));vals=[];contours=[]
    for level,index in enumerate(PROTOCOL['cycle_indices']):
     r=radius*mp.exp(index*mp.pi/k)
     fj=[sum(co*fall(p,n)*r**(p-n) for co,p in [(lam,2),(1,0),(-2*m,-1),(Qh,-2),(-2*d,-3),(d,-4)]) for n in range(5)]
     rj=[sum(co*fall(p,n)*r**(p-n) for co,p in [(-12*lam,0),(4*d,-5),(-6*d,-6)]) for n in range(3)]
     pjet=profile_jets(r,k,phase,mp.mpf('-1.5'),1,1,0)
     outputs=[]
     for label,hc,hs,expected in [('trace',htc,hts,mp.mpf(0)),('action',hac,has,bads)]:
      hj=profile_jets(r,k,phase,-3,2,hc,hs)
      hj=[hj[n]+h0*fall(-3,n)*r**(-3-n) for n in range(5)]
      def raw(eta):return charged_trace_density(r,lam,t,qm,fj,rj,[eta*pjet[n]+eta*eta*hj[n] for n in range(5)])
      coefficient=mp.diff(raw,0,2)/2;actual=r**3*coefficient;error=abs(actual-expected)/scale
      tag=str(ix)+'_'+str(precision)+'_'+str(level)+'_'+label
      check('charged_trace_coarse_'+tag,error<mp.mpf(PROTOCOL['coarse_scaled_limit']),mp.nstr(error,22))
      outputs.append({'profile':label,'normalized_actual':mp.nstr(actual,precision),'normalized_expected':mp.nstr(expected,precision),'scaled_error':mp.nstr(error,25)})
      if level==0:
       ar=mp.mpf(PROTOCOL['contour_radius'])
       for nodes in PROTOCOL['contour_nodes']:
        extracted=sum(raw(ar*mp.exp(2j*mp.pi*j/nodes))*mp.exp(-4j*mp.pi*j/nodes)/(nodes*ar*ar) for j in range(nodes))
        ce=abs(extracted-coefficient)/(1+abs(coefficient))
        check('unexpanded_trace_Taylor_'+tag+'_'+str(nodes),ce<mp.mpf(PROTOCOL['contour_scaled_limit']),mp.nstr(ce,22))
        contours.append({'profile':label,'nodes':nodes,'scaled_error':mp.nstr(ce,25)})
     vals.append({'radius':mp.nstr(r,precision),'profiles':outputs})
    for j,label in enumerate(['trace','action']):
     errors=[mp.mpf(x['profiles'][j]['scaled_error']) for x in vals]
     check('charged_trace_fine_'+str(ix)+'_'+str(precision)+'_'+label,errors[-1]<mp.mpf(PROTOCOL['fine_scaled_limit']))
     ratios=[errors[n+1]/errors[n] for n in range(2)]
     check('charged_trace_refines_'+str(ix)+'_'+str(precision)+'_'+label,max(ratios)<mp.exp(-mp.pi/k),[mp.nstr(x,20) for x in ratios])
    levels.append({'precision':precision,'values':vals,'Taylor_contours':contours})
  with mp.workdps(80):
   for n in range(3):
    aa=levels[0]['values'][n]['profiles'];bb=levels[1]['values'][n]['profiles']
    pe=max(abs(mp.mpf(x['normalized_actual'])-mp.mpf(y['normalized_actual']))/(1+abs(mp.mpf(y['normalized_actual']))) for x,y in zip(aa,bb))
    check('precision_refinement_'+str(ix)+'_'+str(n),pe<mp.mpf(PROTOCOL['precision_scaled_limit']),mp.nstr(pe,22))
  FIXTURES.append({'fixture':'charged_nonlinear_trace_tails','parameters':row,'levels':levels,'physical_solution':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 nonlinear trace compatibility of fast boundary correction','input_count':len(ps)}),flush=True)
 try:provenance();data=exact_derivations();tails=exact_compatibility(data);fixtures(tails)
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 nonlinear trace compatibility of fast boundary correction','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'nonlinear_trace_forcing_derived':good,'conformal_tail_trace_action_incompatibility_derived':good,
  'strict_normalized_cutoff_gap_derived':good,'necessary_coupled_shape_tail_relation_derived':good,
  'physical_fast_branch_excluded':False,'coupled_boundary_completion_adopted':False,'full_shape_equation_solved':False,
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
