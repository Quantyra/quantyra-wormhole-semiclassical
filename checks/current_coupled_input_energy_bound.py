"""Universal complete seed input-energy bound; not a strong coupled solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb
import json, sys, traceback
import sympy as sp
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
import flint
from flint import arb, ctx
ctx.prec=512
ROOT=Path(__file__).resolve().parents[1]
STEM='current-coupled-input-energy-bound'
OWN=['notes/'+STEM+'.md','checks/current_coupled_input_energy_bound.py',
     'evidence/'+STEM+'-review.json']
sys.path.insert(0,str(ROOT/'checks'))
saved_argv=sys.argv[:]
sys.argv=[sys.argv[0]]
import current_physical_york_dual_assembly as york
sys.argv=saved_argv
CHECKS=[]
def check(name,value,detail=None):
    CHECKS.append(dict(name=name,passed=bool(value)))
    if not value: raise AssertionError(name)
york.check=check
york.prior.check=check
def enc(v):
    m,r,e=v.mid_rad_10exp(60)
    return dict(ball=v.str(40),decimal_enclosure=dict(mid=str(m),rad=str(r),exponent=str(e)))

class Jet:
    """Positive integral product majorants justified in the ordinary argument."""
    def __init__(self,v=0,d=0,dd=0): self.v,self.d,self.dd=map(arb,(v,d,dd))
    @staticmethod
    def lift(x):return x if isinstance(x,Jet) else Jet(x)
    def __add__(self,x):
        x=self.lift(x);return Jet(self.v+x.v,self.d+x.d,self.dd+x.dd)
    __radd__=__add__
    def __mul__(self,x):
        x=self.lift(x)
        return Jet(self.v*x.v,self.d*x.v+self.v*x.d,self.dd*x.v+2*self.d*x.d+self.v*x.dd)
    __rmul__=__mul__
    def __truediv__(self,x):return self*(1/arb(x))
    def __pow__(self,k):
        out=Jet(1)
        for _ in range(k):out=out*self
        return out

def execute():
    z=sp.symbols('z',positive=True)
    Q=sp.Function('Q')(z);X=sp.Function('X')(z)
    T=sp.Function('T')(z);nu=sp.Function('nu')(z)
    A=sp.diff(Q,z,2)-4*sp.diff(Q,z)/z
    x2=4*sp.diff(X,z)/z+sp.Rational(2,3)*A*X/Q-sp.Rational(2,3)*sp.diff(T,z)+sp.Rational(8,3)*T/z
    J=X/Q;xi=-T-sp.diff(X,z)/2
    residual=-sp.diff(Q**2/z**4*sp.diff(J,z),z)-Q/z**4*(4*xi/z-sp.diff(xi,z))
    check('exact_inverse_York_equation',sp.simplify(residual.subs(sp.diff(X,z,2),x2))==0)
    scalar=nu+T/2+sp.diff(X,z)/4-X/z
    lap=Q*sp.diff(scalar,z,2)+(sp.diff(Q,z)-2*Q/z)*sp.diff(scalar,z)
    lap=lap.subs(sp.diff(X,z,3),sp.diff(x2,z)).subs(sp.diff(X,z,2),x2)
    y=sp.Function('y')(z)
    lap=lap.subs({sp.diff(X,z):2*z*y+z*z*sp.diff(y,z),X:z*z*y})
    by=Q+z*z*sp.diff(Q,z,2)/6-2*z*sp.diff(Q,z)/3
    cy=z**3*sp.diff(Q,z,3)/6-2*z*z*sp.diff(Q,z,2)/3+5*z*sp.diff(Q,z)/3-2*Q
    expected=Q*sp.diff(nu,z,2)+(sp.diff(Q,z)-2*Q/z)*sp.diff(nu,z)+sp.diff(Q*sp.diff(T,z),z)/3+(2*sp.diff(Q,z)/(3*z)-2*Q/z**2)*T+by*sp.diff(y,z)+cy*y/z
    check('exact_laplacian_cancellation_no_inverse_Q',sp.simplify(lap-expected)==0)
    p,h,zh,env=york.actual_parameters()
    coeff={k:env[v] for k,v in zip(york.powers,york.qs)}
    beta,u,eps=coeff[0],coeff[2],coeff[6]
    aa,b4=-coeff[3],coeff[4]
    check('potential_coercivity_whole_interval',beta-3*u-5*abs(eps)>0 and 3*aa-5*b4/3>0)
    # Exact Bernstein range over [0,1], including the physical [0,z_h].
    def polynomial_derivative(c):return {k-1:k*v for k,v in c.items() if k}
    def bernstein(c):
        n=max(c)
        return [sum((c.get(i,arb(0))*arb(comb(k,i))/comb(n,i) for i in range(k+1)),arb(0)) for k in range(n+1)]
    q1=polynomial_derivative(coeff);q2=polynomial_derivative(q1);q3=polynomial_derivative(q2)
    ac={k-2:k*(k-5)*v for k,v in coeff.items() if k>=2}
    check('global_Q_guard',all(abs(v)<arb('.17') for v in bernstein(coeff)))
    check('global_Qprime_guard',12*beta*arb(4)/27+u/2+6*abs(eps)+3*p['tau']<arb('.3'))
    check('global_Qsecond_guard',all(abs(v)<arb('2.1') for v in bernstein(q2)))
    check('global_Qthird_guard',all(abs(v)<arb('8.1') for v in bernstein(q3)))
    check('global_A_guard',all(abs(v)<3 for v in bernstein(ac)))
    zc=arb('1e-9')
    qlo=beta-sum((abs(v)*zc**k for k,v in coeff.items() if k),arb(0))
    check('negative_potential_part_guard',qlo>arb('.16') and (6*u+6*abs(eps))/qlo<arb('1e-9') and zc*(6*aa-4*b4)-6*u-6*abs(eps)>0)
    c=2*arb(17).sqrt()/9
    s0=arb('.5')+c/4+c/(2*(arb(9)/2).sqrt())
    grad=(2*arb('.17').sqrt()/3+c*((3*(arb('1.5')+arb('2e-9')/(arb(9)/2))).sqrt()/6+(arb('.17')/(arb(9)/2)).sqrt()))
    lap0=arb('.54')+arb('.72')*c+arb('3.59')*c/(arb(9)/2).sqrt()
    check('weighted_value_gradient_laplacian_guards',s0<1 and grad<1 and lap0<3 and 3+2*arb('.17').sqrt()<4)
    check('complete_sixth_moment_Z_guard',arb(5)/8+arb('1.5')*arb('.3825').sqrt()+arb('.3825')<2)
    C,B2,q=p['C'],p['B2'],arb(p['q'])
    phi,L,MM,NN,DD=Jet(1,1),Jet(3,2),Jet(1,4),Jet(1,1,2),Jet(6,20,12)
    exponent=lambda k:Jet(1,abs(k),k*k)
    en2,en4,eu2,em2=map(exponent,(2,4,2,-2))
    la0=arb(2000)**-10
    parts={
      'Einstein':(en2*DD/2+en2*eu2/C+3*la0*B2/C*en4*eu2).dd,
      'magnetic':(arb.pi()*q*q*arb(2000)**8/(C*B2)*eu2).dd,
      'R2':(8*arb.pi()*arb('1e22')/(C*B2)*(4*eu2+4*C*DD+C*C*em2*DD**2)).dd,
      'Weyl':((C*C*em2*(phi*L**2/60+L*MM/30+MM**2/10+L*NN/30+7*MM*NN/90+7*NN**2/180)+C*phi*L/18+(arb(1)/15+q*q/3)*phi*eu2+C*NN/18)/(2*arb.pi()*C*B2)).dd,
      'vacuum_interval':2*(604*q*q+(1+arb('6e-9'))*(q*q/6-arb(1)/30))/(arb.pi()*C*B2),
      'a2':100*C/(480*arb.pi()*B2)*(4*arb('2.1')**2+16*arb('2.1')+8),
      'LLL':q*arb('.17')/(6*arb.pi()*B2),
      'complete_positive_geodesic':2*(4*C*C/q)*2/(525*arb.pi()*B2)}
    upper=sum(parts.values(),arb(0))
    conversion=(2*upper/(beta*(1-arb('.0022')))).sqrt()
    check('complete_comparison_diagonal_guard',upper<arb('78.000004'))
    check('original_energy_input_constant',conversion<31)
    B,G,gamma,NH,lam,Cs=sp.symbols('B G gamma NH lam C',positive=True)
    W=NH*sp.pi*B**2/(G*gamma);period=2*sp.pi*B/(Cs*gamma)
    check('original_action_to_radial_energy_units',sp.simplify(8*sp.pi*G/(lam*B**3)*W/(NH*4*sp.pi*period)-Cs/(lam*B**2))==0)
    eta=arb('.0022000000000000032');qs=eta/(1-eta)
    check('actual_Rs_input_conversion_upward_guard',31*eta<arb('.0682000000001'))
    check('conditional_graph_tail_product_allowance',arb('2.2085e157')*arb('.0682000000001')<arb('1.5062e156'))
    return dict(outcome='complete_strong_coupled_estimate_unresolved',
      normalized_component_upper_bounds={k:enc(v) for k,v in parts.items()},
      normalized_single_diagonal_upper=enc(upper),original_energy_conversion_upper_expression=enc(conversion),
      original_energy_conversion_guard='31',weighted_York_value=enc(s0),weighted_York_gradient=enc(grad),weighted_York_Laplacian=enc(lap0),
      complete_H_inverse_Rs_energy_per_N=enc(31*eta),complete_Hs_inverse_Rs_energy_per_N=enc(31*qs),
      conditional_Mgraph_times_Nfeedback_allowance='2.2085e157; only for source map R_s times the specified N input',
      norm_definition='N(nu)+N(T), weighted value/first/second and Q^(3/2) third derivatives; see ordinary argument',
      independent_phase_energy_scope=True,individual_cutoff_column_norms_evaluated=False,
      physical_operator_applications=0,coupled_signed_solve_executed=False,
      complete_strong_estimate_executed=False,actual_mass=None,kappa_mass=None,
      M_graph_4=None,M_X_8=None,U_original_matching_feedback=None,
      original_physical_X_comparisons=False,nonlinear_neighborhood_transfer=False,
      new_nonlinear_feasibility_gate_closed=False,nonlinear_admission=False,
      nonlinear_contraction=False,exact_compatibility=False,physical_matching=False,
      validity_noise=False,all_sector_stability=False,actual_controlled_solution=False,
      construction_obstruction=False,closed_audit_preserved=True,
      complete_angular_thermal_contacts_boundary_returns_retained=True)

if __name__=='__main__':
    mode='acceptance' if '--acceptance' in sys.argv else 'development'
    out=ROOT/'evidence'/(STEM+('' if mode=='acceptance' else '-development')+'.json')
    review=json.loads((ROOT/OWN[2]).read_bytes())
    paths=set(OWN+review['retained_inputs'])
    for module in list(sys.modules.values()):
        f=getattr(module,'__file__',None)
        if f:
            try:rel=Path(f).resolve().relative_to(ROOT).as_posix()
            except ValueError:continue
            if rel.endswith('.py'):paths.add(rel)
    digest=lambda p:sha256(p.read_bytes()).hexdigest()
    before={p:digest(ROOT/p) for p in sorted(paths)}
    fr=Path(flint.__file__).resolve().parent
    runtime=dict(python=sys.version,python_executable_sha256=digest(Path(sys.executable)),flint_version=flint.__version__,
      flint_artifacts={p.relative_to(fr).as_posix():digest(p) for p in sorted(fr.rglob('*')) if p.is_file() and p.suffix in ['.py','.pyd','.dll']})
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),mode=mode,input_sha256_before=before,precision_bits=512,runtime=runtime,
      frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN})
    try:
        check('ordinary_argument_frozen',before[OWN[0]]==review['ordinary_argument_sha256'])
        run['result']=execute()
        run['input_sha256_after']={p:digest(ROOT/p) for p in before}
        check('frozen_inputs_unchanged',before==run['input_sha256_after'])
        run['status']='passed'
    except Exception:run['status']='failed';run['traceback']=traceback.format_exc()
    run['checks']=CHECKS;run['finished_utc']=datetime.now(timezone.utc).isoformat()
    ledger=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    ledger['runs'].append(run);out.write_bytes((json.dumps(ledger,indent=2)+'\n').encode())
    with (ROOT/'evidence'/(STEM+'.log')).open('ab') as f:f.write((json.dumps(run)+'\n').encode())
    print(json.dumps(dict(mode=mode,status=run['status'],checks=len(CHECKS),result=run.get('result'),traceback=run.get('traceback')),indent=2))
    if run['status']!='passed':raise SystemExit(1)
