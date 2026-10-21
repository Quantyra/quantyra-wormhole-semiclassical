"""Energy-domain inclusion and complete chart contact; nonlocal cut comparison open."""
from datetime import datetime,timezone
from hashlib import sha256
from pathlib import Path
import json,sys,traceback
import sympy as S
import mpmath as mp
ROOT=Path(__file__).resolve().parents[1];STEM='energy-domain-and-chart-contact'
PREVIOUS='finite-cut-normal-shape-limit'
OWN=['checks/energy_domain_and_chart_contact.py','notes/'+STEM+'.md','evidence/energy-domain-chart-source-review.json']
ROWS=[];FIXTURES=[];R=S.Rational
PROTOCOL={'precision_digits':[50,80],'matrix_cases':[['1','.2'],['2','.35'],['4','.5']],
 'Gaussian_cases':[['.5','0','.2'],['1','2','.3'],['2','5','.1']],
 'finite_part_orders':[1,2,3],'cut_powers':[3,5,7],'spike_powers':[2,4,6],
 'contour_nodes':[24,48],'contour_radius':'.02','scaled_limit':'1e-30',
 'precision_scaled_limit':'1e-32','spike_relative_limit':'1e-18',
 'actual_charged_hessian_evaluated':False,'continuum_proofs_are_ordinary':True,
 'physical_action_averaging_adopted':False,'finite_regulator_matrices_are_diagnostics':True}

def check(name,ok,detail=''):
 assert name not in {r['name'] for r in ROWS}
 row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)
def zero(name,value):
 value=S.factor(S.cancel(S.sympify(value).doit()));check(name,value==0,value)
def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})
def err(a,b):return abs(a-b)/(1+abs(b))
def textnum(a):return mp.nstr(a,mp.mp.dps)
def dot(a,b):return sum(a[i]*b[i] for i in range(len(a)))
def tr(A):return sum(A[i,i] for i in range(A.rows))
def oplog(A):
 e,Q=mp.eigsy(A);return Q*mp.diag([mp.log(x) for x in e])*Q.T
def opnorm(A):return max(abs(x) for x in mp.eigsy(A,eigvals_only=True))

def provenance():
 p=prior();check('predecessor_passed',p['status']=='passed' and len(p['checks'])==179 and all(c['passed'] for c in p['checks']))
 check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
 bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
 check('predecessor_hashes',not bad,bad)
 review=json.loads((ROOT/OWN[2]).read_bytes())
 check('review_inputs_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
 check('scope',not review['actual_nonlocal_finite_cut_comparison_proved'] and not review['physical_fast_branch_excluded'] and not review['bare_tadpole_independently_renormalized'])

def symbolic():
 y=S.symbols('y',positive=True);h=S.Function('h')(y);u=S.Function('u')(y)
 T=S.exp(-h/2)*u
 raw=-S.exp(h/2)*(S.diff(T,y,2)+S.diff(h,y)*S.diff(T,y))
 zero('scalar_half_density',raw+S.diff(u,y,2)-(S.diff(h,y,2)/2+S.diff(h,y)**2/4)*u)
 m=S.Function('m')(y);v=S.Function('v')(y)
 zero('density_second_derivative',S.diff(m*v,y,2)-m*S.diff(v,y,2)-2*S.diff(m,y)*S.diff(v,y)-S.diff(m,y,2)*v)
 eta,a=S.symbols('eta a');f=S.Function('f')
 chain=S.diff(f(S.log(1+eta*a)),eta,2).subs(eta,0)
 zero('complete_chart_chain_rule',chain-a*a*S.diff(f(eta),eta,2).subs(eta,0)+a*a*S.diff(f(eta),eta).subs(eta,0))
 chi=S.kronecker_product(S.diag(1,-1),S.Matrix([[0,1],[1,0]]))
 zero('balanced_third_smearing_contact',S.trace(chi))
 x,Y,s=S.symbols('x Y s');b=630*(Y-1)**4*(2-Y)**4
 zero('averaging_test_mass',S.integrate(b,(Y,1,2))-1)
 poly=S.Poly(S.expand(b.subs(Y,Y-x)),x);coeff=[poly.nth(j) for j in range(9)];kernels={}
 for p in PROTOCOL['finite_part_orders']:
  residue=coeff[p-1]
  zero('finite_part_residue_'+str(p),residue-(-1)**(p-1)*S.diff(b,Y,p-1)/S.factorial(p-1))
  regular=S.cancel(sum(coeff[j]*(Y-1)**(j-p+1)/S.Rational(j-p+1) for j in range(9) if j!=p-1))
  regular=S.expand(regular.subs(Y,1+s));residue=S.expand(residue.subs(Y,1+s))
  check('finite_part_local_polynomial_'+str(p),regular.is_polynomial(s) and residue.is_polynomial(s))
  kernels[p]=(S.lambdify(s,regular,'mpmath'),S.lambdify(s,residue,'mpmath'))
 bump=x*x*(1-x*x)**6/2;jets=[S.diff(bump,x,j) for j in range(4)]
 zero('bump_value',jets[0].subs(x,0));zero('bump_slope',jets[1].subs(x,0));zero('bump_curvature',jets[2].subs(x,0)-1)
 norm=S.integrate(jets[2]**2,(x,-1,1))
 zero('bump_energy_exponent',2*13-3*8-2)
 zero('squared_shape_contact_exponent',2+13-2*8+1)
 return kernels,S.lambdify(x,jets,'mpmath'),S.lambdify((),norm,'mpmath')

def matrix_fixtures():
 for ix,(ns,amps) in enumerate(PROTOCOL['matrix_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    nu=mp.mpf(ns);amp=mp.mpf(amps);P=mp.matrix(4);W=mp.diag([amp,-amp/2,3*amp/4,-amp]);I=mp.eye(4)
    for i in range(4):
     P[i,i]=nu*nu+2
     if i<3:P[i,i+1]=P[i+1,i]=-1
    H=P+W;LH=oplog(H);LP=oplog(P);f=mp.matrix([1,mp.mpf('.3'),mp.mpf('-.2'),mp.mpf('.7')])
    m=nu*nu-amp;ratio=opnorm(LH-LP)*m/amp
    check('operator_log_bound_'+str(ix)+'_'+str(digits),ratio<=1,mp.nstr(ratio,25))
    integral=mp.quad(lambda t:dot(mp.lu_solve(H+t*I,f),W*mp.lu_solve(P+t*I,f)),[0,1,10,100,mp.inf])
    target=dot(f,(LH-LP)*f);e=err(integral,target)
    check('noncommuting_resolvent_'+str(ix)+'_'+str(digits),e<mp.mpf(PROTOCOL['scaled_limit']),mp.nstr(e,25))
    hu=H*f;pu=P*f;wu=W*f
    lhs=dot(hu,LH*hu)-dot(pu,LP*pu)
    rhs=dot(hu,(LH-LP)*hu)+2*dot(pu,LP*wu)+dot(wu,LP*wu)
    qe=err(lhs,rhs);check('log_energy_split_'+str(ix)+'_'+str(digits),qe<mp.mpf(PROTOCOL['scaled_limit']),mp.nstr(qe,25))
    M=mp.diag([mp.mpf('1.1'),mp.mpf('1.6'),mp.mpf('2.2'),mp.mpf('3')])
    D=mp.matrix([[0,amp,0,amp/2],[-amp,0,2*amp,0],[0,-2*amp,0,3*amp/2],[-amp/2,0,-3*amp/2,0]])
    vv=[mp.mpf('.2')*(ix+1),mp.mpf('-.4'),mp.mpf('.6'),mp.mpf('-.3')]
    V=mp.diag(vv);weights=[mp.mpf('.3'),mp.mpf('-.2'),mp.mpf('.15'),mp.mpf('.25')]
    inv=(D+M)**-1;vertex=M*V;bare=-tr(inv*M*V*V)
    restored=bare+sum(weights[i]*vv[i]**2 for i in range(4))
    local2=mp.mpf('.2')*sum((vv[i+1]-vv[i])**2 for i in range(3))
    exp2=tr(inv*vertex*inv*vertex)+bare+local2
    aff2=tr(inv*vertex*inv*vertex)-sum(weights[i]*vv[i]**2 for i in range(4))+local2
    def action(eta,affine):
     vals=[mp.log(1+eta*v) if affine else eta*v for v in vv]
     mass=mp.diag([M[i,i]*mp.exp(vals[i]) for i in range(4)])
     return -mp.log(mp.det(D+mass))+sum(weights[i]*vals[i] for i in range(4))+mp.mpf('.1')*sum((vals[i+1]-vals[i])**2 for i in range(3))
    contours=[]
    for n in PROTOCOL['contour_nodes']:
     radius=mp.mpf(PROTOCOL['contour_radius'])
     derivatives=[2*mp.fsum(action(radius*mp.exp(2j*mp.pi*j/n),aff)*mp.exp(-4j*mp.pi*j/n) for j in range(n))/(n*radius**2) for aff in [False,True]]
     ce=max(err(derivatives[0],exp2),err(derivatives[1],aff2),err(derivatives[0]-derivatives[1],restored))
     check('complete_chart_contour_'+str(ix)+'_'+str(digits)+'_'+str(n),ce<mp.mpf(PROTOCOL['scaled_limit']),mp.nstr(ce,25))
     contours.append({'nodes':n,'scaled_error':mp.nstr(ce,25)})
    check('chart_contact_not_zero_'+str(ix)+'_'+str(digits),abs(restored)>mp.mpf('.01'),mp.nstr(restored,25))
    levels.append({'precision':digits,'log_pairing':textnum(target),'bound_ratio':mp.nstr(ratio,25),'resolvent_scaled_error':mp.nstr(e,25),'energy_split_scaled_error':mp.nstr(qe,25),'restored_chart_contact':textnum(restored),'contours':contours})
  with mp.workdps(80):
   pe=max(err(mp.mpf(levels[0][k]),mp.mpf(levels[1][k])) for k in ['log_pairing','restored_chart_contact'])
   check('matrix_precision_'+str(ix),pe<mp.mpf(PROTOCOL['precision_scaled_limit']),mp.nstr(pe,25))
  FIXTURES.append({'fixture':'noncommuting_log_and_chart','case':ix,'levels':levels,'actual_quantum_response':False})

def Mellin_fixtures():
 for ix,(ss,ks,als) in enumerate(PROTOCOL['Gaussian_cases']):
  levels=[]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    s,k,al=map(mp.mpf,[ss,ks,als]);s2=s/mp.sqrt(1+2*s*s)
    ft=lambda xi,sw:sw/2*(mp.exp(-sw*sw*(xi-k)**2/2)+mp.exp(-sw*sw*(xi+k)**2/2))
    symbol=lambda xi:mp.re(mp.digamma(mp.mpf('.5')+1j*xi))-mp.digamma(mp.mpf('.5'))
    splits=sorted(set([mp.mpf(0),k,k+1/s,k+4/s,k+10/s]))+[mp.inf]
    energy=2*mp.quad(lambda xi:symbol(xi)*ft(xi,s)**2,splits)
    multiplied=2*mp.quad(lambda xi:symbol(xi)*(ft(xi,s)+al*ft(xi,s2))**2,splits)
    derivative2=2*mp.quad(lambda xi:xi*xi*ft(xi,s)**2,splits)
    norm=mp.sqrt(mp.pi)*s/2*(1+mp.exp(-k*k*s*s))
    def corr_gap(r):return norm*(-mp.expm1(-r*r/(4*s*s)))+mp.sqrt(mp.pi)*s*mp.exp(-r*r/(4*s*s))*mp.sin(k*r/2)**2
    real=mp.quad(lambda r:corr_gap(r)/(2*mp.sinh(r/2)),[0,1,2,4,8,16,mp.inf])
    e=err(real,energy);check('Mellin_kernel_energy_'+str(ix)+'_'+str(digits),e<mp.mpf(PROTOCOL['scaled_limit']),mp.nstr(e,25))
    Lip=al*mp.sqrt(2/mp.e);bound=2*(1+al)**2*energy+14*mp.zeta(3)*Lip**2*norm
    ratio=multiplied/bound;check('Mellin_multiplier_bound_'+str(ix)+'_'+str(digits),ratio<=1,mp.nstr(ratio,25))
    h1ratio=energy/(7*mp.zeta(3)*derivative2)
    check('Mellin_H1_inclusion_'+str(ix)+'_'+str(digits),h1ratio<=1,mp.nstr(h1ratio,25))
    levels.append({'precision':digits,'energy':textnum(energy),'multiplied_energy':textnum(multiplied),'kernel_scaled_error':mp.nstr(e,25),'multiplier_bound_ratio':mp.nstr(ratio,25),'H1_bound_ratio':mp.nstr(h1ratio,25)})
  with mp.workdps(80):
   pe=err(mp.mpf(levels[0]['energy']),mp.mpf(levels[1]['energy']))
   check('Mellin_precision_'+str(ix),pe<mp.mpf(PROTOCOL['precision_scaled_limit']),mp.nstr(pe,25))
  FIXTURES.append({'fixture':'Mellin_density_domain_transfer','case':[ss,ks,als],'levels':levels,'actual_shape_solution':False})
 with mp.workdps(80):
  mass=mp.quad(lambda r:r*r/(2*mp.sinh(r/2)),[0,1,4,16,mp.inf])
  check('Mellin_second_moment',err(mass,14*mp.zeta(3))<mp.mpf('1e-60'),mp.nstr(mass,25))

def gammafp(n):
 if n>0:return mp.factorial(n-1)
 m=-n;return (-1)**m/mp.factorial(m)*(mp.harmonic(m)-mp.euler)
def analytic_fp(p,d):return mp.exp(-d)*sum(mp.binomial(4,j)*d**(4-j)*gammafp(j-p+1) for j in range(5))
def direct_fp(p,d):
 bs=[mp.binomial(4,j)*d**(4-j) for j in range(5)]
 cs=[sum(bs[j]*(-1)**(l-j)/mp.factorial(l-j) for j in range(l+1)) for l in range(p)]
 def rem(x):
  return sum(bs[j]*(-1)**(p-j)/mp.factorial(p-j)*mp.hyp1f1(1,p-j+1,-x) for j in range(p))+sum(bs[j]*x**(j-p)*mp.exp(-x) for j in range(p,5))
 return mp.exp(-d)*(mp.quad(rem,[0,1])+mp.quad(lambda x:mp.exp(-x)*(d+x)**4/x**p,[1,mp.inf])+sum(cs[l]/(l-p+1) for l in range(p-1)))

def distribution_fixtures(kernels):
 for p in PROTOCOL['finite_part_orders']:
  levels=[];cf,res=kernels[p]
  for digits in PROTOCOL['precision_digits']:
   with mp.workdps(digits):
    de=err(direct_fp(p,mp.mpf(1)/8),analytic_fp(p,mp.mpf(1)/8))
    check('direct_finite_part_'+str(p)+'_'+str(digits),de<mp.mpf(PROTOCOL['scaled_limit']),mp.nstr(de,25))
    values=[]
    for n in PROTOCOL['cut_powers']:
     delta=mp.mpf(2)**(-n);U=lambda y:y**4*mp.exp(-y)
     direct=mp.quad(lambda a:630*(a-1)**4*(2-a)**4*analytic_fp(p,delta*a),[1,2])
     local=mp.quad(lambda s:U(delta*(1+s))*(cf(s)+(mp.log(s)+mp.log(delta))*res(s)),[0,mp.mpf('.5'),1])
     tail=mp.quad(lambda y:U(delta*y)/(y-1)**p*mp.hyp2f1(p,5,10,1/(y-1)),[2,3,1/delta,mp.inf])
     dual=delta**(1-p)*(local+tail)
     scalecontact=delta**(1-p)*mp.log(delta)*mp.quad(lambda s:U(delta*(1+s))*res(s),[0,1])
     e=err(direct,dual);check('transposed_cut_'+str(p)+'_'+str(digits)+'_'+str(n),e<mp.mpf(PROTOCOL['scaled_limit']),mp.nstr(e,25))
     check('scale_contact_required_'+str(p)+'_'+str(digits)+'_'+str(n),abs(scalecontact)>mp.mpf('1e-14'),mp.nstr(scalecontact,25))
     limit=mp.factorial(4-p);limitratio=abs(direct-limit)/delta
     check('cubic_envelope_diagnostic_'+str(p)+'_'+str(digits)+'_'+str(n),limitratio<10,mp.nstr(limitratio,25))
     values.append({'power':n,'average':textnum(direct),'scaled_duality_error':mp.nstr(e,25),'scale_contact':mp.nstr(scalecontact,25),'limit_error_over_cut':mp.nstr(limitratio,25)})
    levels.append({'precision':digits,'direct_fp_scaled_error':mp.nstr(de,25),'values':values})
  with mp.workdps(80):
   pe=max(err(mp.mpf(a['average']),mp.mpf(b['average'])) for a,b in zip(levels[0]['values'],levels[1]['values']))
   check('finite_part_precision_'+str(p),pe<mp.mpf(PROTOCOL['precision_scaled_limit']),mp.nstr(pe,25))
  FIXTURES.append({'fixture':'translated_cubic_test_finite_parts','order':p,'levels':levels,'actual_heat_kernel_evaluated':False})

def spike_fixtures(bump,normfun):
 levels=[]
 for digits in PROTOCOL['precision_digits']:
  with mp.workdps(digits):
   values=[]
   for n in PROTOCOL['spike_powers']:
    d=mp.mpf(2)**(-n);ell=d**8;amplitude=d**13
    def parts(x):
     y=d+ell*x;j=bump(x);b=amplitude*j[0];b1=amplitude*j[1]/ell;b2=amplitude*j[2]/ell**2
     du=2*y*y*b+b*b;du2=4*b+8*y*b1+2*y*y*b2+2*b1*b1+2*b*b2
     return y,du,du2,b2,amplitude*j[3]/ell**3
    direct=ell*mp.quad(lambda x:parts(x)[2]/parts(x)[0],[-1,0,1])
    dual=2*ell*mp.quad(lambda x:parts(x)[1]/parts(x)[0]**3,[-1,0,1])
    re=abs(direct-dual)/abs(dual)
    check('squared_shape_weak_contact_'+str(digits)+'_'+str(n),re<mp.mpf(PROTOCOL['spike_relative_limit']),mp.nstr(re,25))
    N=ell*mp.quad(lambda x:parts(x)[3]**2,[-1,0,1])
    N1=ell*mp.quad(lambda x:(parts(x)[3]/2+parts(x)[0]*parts(x)[4])**2,[-1,0,1])
    check('spike_H2_energy_'+str(digits)+'_'+str(n),err(N,d*d*normfun())<mp.mpf(PROTOCOL['scaled_limit']))
    peak=12*d*d+2/d;ratio=abs(dual)/d**20
    check('weak_contact_decay_'+str(digits)+'_'+str(n),ratio<10 and peak>2/d,mp.nstr(ratio,25))
    values.append({'power':n,'pointwise_second_jet':textnum(peak),'weak_contact':textnum(dual),'relative_duality_error':mp.nstr(re,25),'H2_cost':textnum(N),'log_energy_Jensen_quantity':textnum(N*(1+mp.log(1+N1/N))),'decay_scaled_ratio':mp.nstr(ratio,25)})
   levels.append({'precision':digits,'values':values})
 with mp.workdps(80):
  pe=max(abs(mp.mpf(a['weak_contact'])-mp.mpf(b['weak_contact']))/abs(mp.mpf(b['weak_contact'])) for a,b in zip(levels[0]['values'],levels[1]['values']))
  check('spike_precision',pe<mp.mpf('1e-30'),mp.nstr(pe,25))
 FIXTURES.append({'fixture':'squared_shape_concentrated_jets','levels':levels,'actual_shape_solution':False})

def main():
 sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
 print(json.dumps({'event':'start','story':'S168 energy domain and chart contact','input_count':len(ps)}),flush=True)
 try:
  provenance();kernels,bump,norm=symbolic();matrix_fixtures();Mellin_fixtures();distribution_fixtures(kernels);spike_fixtures(bump,norm)
 except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
 run={'story':'S168 energy domain and chart contact','created_utc':datetime.now(timezone.utc).isoformat(),
  'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
  'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
  'actual_energy_domain_transfer_derived':good,'complete_chart_contact_weak_limit_derived':good,
  'actual_nonlocal_finite_cut_comparison_proved':False,'pure_shape_finite_cut_limit_proved':False,
  'full_boundary_action_solved':False,'physical_action_averaging_adopted':False,'physical_fast_branch_excluded':False,
  'new_boundary_prescription_adopted':False,'bare_tadpole_independently_renormalized':False,'complete_actual_second_variation_evaluated':False,
  'rough_finite_amplitude_action_differentiability_proved':False,'full_shape_equation_solved':False,
  'static_source_compatibility_proved':False,'regulator_independence_for_moving_cuts_proved':False,
  'continuous_causal_stability_proved':False,'physical_constraint_mode_constructed':False,'actual_wormhole_instability_proved':False,
  'finite_counterterm_changed':False,'order_reduction_adopted':False,'formal_verification':False,
  'S167_complete':False,'S168_complete':False,'S077_complete':False,'goal_complete':False}
 if not good:
  for kind,path in [('checker',OWN[0]),('note',OWN[1]),('review',OWN[2])]:run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
 dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
 data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
 print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
 return 0 if good else 1
if __name__=='__main__':raise SystemExit(main())
