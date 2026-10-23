"""Endpoint-aware mass action and full local thermal restoration, not a metric solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import comb
import sys,json,traceback
import sympy as sp
import mpmath as mp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_subtracted_thermal_action as thermal
prior=thermal.prior; raw=thermal.raw
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]; STEM='current-endpoint-thermal-restoration'
PREVIOUS='current-subtracted-thermal-action'
OWN=['notes/'+STEM+'.md','checks/current_endpoint_thermal_restoration.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];FIXTURES=[];encoded=prior.encoded;hull=prior.hull

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)
prior.check=check

def dependencies():
    d=prior.dependency_identity()
    for name,module in [('sympy',sp),('mpmath',mp)]:
        path=Path(module.__file__).resolve().parent
        d[name]=dict(version=module.__version__,root=str(path),source_sha256={str(q.relative_to(path)):sha256(q.read_bytes()).hexdigest() for q in path.rglob('*.py')})
    d['mpmath_dps']=mp.mp.dps
    return d

def small_restoration(y,K=16):
    """Return S0,S1,a*T for any positive ball with upper edge <=1/8."""
    assert y>0 and y<=arb(1)/8
    c0=-(arb.const_euler()+(2*y).log())/2
    s0=c0;s1=arb(1)/2
    value=-arb(1)/12+y*y*(arb(1)/2-arb.const_euler()-(2*y).log())
    for k in range(1,K+1):
        ck=-arb((-1)**k*comb(2*k,k)*(2**(2*k+1)-1))/(2*4**k)*arb(2*k+1).zeta()
        term=ck*y**(2*k);s0+=term;s1-=2*k*term
        value+=2*term*y*y/(k+1)
    r=4*y*y; E0=arb('1.203')*r**(K+1)/(1-r)
    E1=arb('2.406')*r**(K+1)*((K+1)/(1-r)+r/(1-r)**2)
    ET=2*y*y*E0/(K+2)
    return s0+arb(0,E0.abs_upper()),s1+arb(0,E1.abs_upper()),value+arb(0,ET.abs_upper()),dict(method='convergent_half_integer_series',terms=K,tails=[encoded(x) for x in [E0,E1,ET]])

def bessel_restoration(y):
    """Alternating Bessel sums with next-term remainders, including all k."""
    assert y>=arb(1)/16
    z=2*arb.pi()*y;s0=arb(0);s1=arb(0);value=arb(0)
    for k in range(1,4097):
        arg=z*k;k0=arg.bessel_k(0);k1=arg.bessel_k(1);sgn=(-1)**(k+1)
        s0+=sgn*k0;s1+=sgn*arg*k1;value-=sgn*2*y*k1/(arb.pi()*k)
        nxt=z*(k+1);n0=nxt.bessel_k(0);n1=nxt*nxt.bessel_k(1)
        nt=2*y*nxt.bessel_k(1)/(arb.pi()*(k+1))
        if all(x.abs_upper()<arb('1e-22') for x in [n0,n1,nt]):
            return s0+arb(0,n0.abs_upper()),s1+arb(0,n1.abs_upper()),value+arb(0,nt.abs_upper()),dict(method='alternating_Bessel_series',terms=k,tails=[encoded(x) for x in [n0,n1,nt]])
    raise AssertionError('Bessel next-term target not met')

def restoration_dimensionless(y):
    """The dimensionless local action a*T and its two b derivatives."""
    if y<=arb(1)/8: s0,s1,value,meta=small_restoration(y)
    elif y>=arb(1)/16: s0,s1,value,meta=bessel_restoration(y)
    else:raise ValueError('Split a wide y interval at the analytic overlap before evaluation')
    return [value,4*y*y*s0,4*y*y*(2*s0-s1)],meta

def finite_cut_terms(ell,first,second):
    """Add these to grouped integrals AND the separately enclosed tip."""
    return [ell*y for y in first],[ell*z for z in second]

def initial_boxes(M,epsilon,omega):
    assert M>0 and epsilon>0 and M*epsilon<=arb('1e-3') and omega>=arb('.5')
    q=omega*(1+omega)
    return [arb(0,c*M*epsilon/q) for c in [2,10,14]]

def cap_budgets(M,epsilon):
    assert M>0 and epsilon>0 and M*epsilon<=arb('1e-3')
    z=M*epsilon
    diff=1550*(6-arb.pi()**2/2)*z*z
    restore=4*z*z*((2/z).log()+arb(7)/2)
    return diff,restore,diff+restore

def algebra():
    r,w,m=sp.symbols('r w m',positive=True)
    k=m/(2*w*(1+2*w));c=-m/(2*w)
    check('regular_rotated_coefficient',sp.simplify(k-(-2*w*k-c))==0)
    j=sp.symbols('j',integer=True,nonnegative=True);b=sp.symbols('b',real=True)
    logC=sp.log(2)/2+sp.log(sp.factorial(j))+j*(sp.log(2)-sp.log(m)-b)
    check('regular_normalization_first_variation',sp.diff(-2*logC,b)==2*j)
    check('regular_normalization_second_variation',sp.diff(-2*logC,b,2)==0)
    check('flat_polar_action_tip_coefficient',sp.simplify(2*(m*m/(4*w)-m*m*(w+1)/(1+2*w)**2)-m*m/(2*w*(1+2*w)**2))==0)
    C=sp.symbols('C');ell=r*(C-2*m*sp.log(r))
    check('lowest_channel_adjoint_log',sp.simplify(sp.diff(ell,r)-ell/r+2*m)==0)
    e,Y,Z=sp.symbols('e Y Z');ap,lp,yp,zp=sp.symbols('ap lp yp zp')
    check('finite_cut_gradient_boundary_sign',sp.expand((ap*e-lp)*Y-e*(yp+ap*Y)+(lp*Y+e*yp))==0)
    check('finite_cut_hessian_boundary_sign',sp.expand((ap*e-lp)*Z-e*(zp+ap*Z)+(lp*Z+e*zp))==0)
    check('thermal_cap_partial_fraction',sp.cancel(1/(w*(1+w)**2)-(1/w-1/(w+1)-1/(w+1)**2))==0)
    delta=arb('1e-3');w0=arb('.5');q0=w0*(1+w0)
    check('positive_sensitivity_damping',4*delta**2/(w0**3*(1+w0))<arb('.5'))
    check('first_sensitivity_barrier_forcing',4+8*delta**2/q0<5)
    check('second_sensitivity_barrier_forcing',6+96*delta**2/q0+360*delta**2/q0**2<7)
    check('all_frequency_gradient_Hessian_cap_coefficient',480+32*delta**2/w0**2<3100 and 3072+704*delta**2/w0**2<3100)
    check('exact_thermal_cap_sum_bound',6-arb.pi()**2/2>1 and 6-arb.pi()**2/2<arb('1.1'))
    check('series_zeta_majorant',arb(3).zeta()<arb('1.203'))
    z=sp.symbols('z',positive=True)
    check('small_mass_restoration_primitive',sp.simplify(z*sp.diff(z*z*(sp.Rational(1,2)-sp.EulerGamma-sp.log(2*z)),z)-4*z*z*(-sp.EulerGamma-sp.log(2*z))/2)==0)
    n,q,Cc=sp.symbols('n q Cc',positive=True)
    angular=sp.expand(2*(q+2*n)*n*(n+q)/Cc)
    check('fixed_cap_angular_allowance_is_not_summable',sp.expand(angular-4*n**3/Cc-6*q*n*n/Cc-2*q*q*n/Cc)==0)
    NUMERICS.append(dict(kind='exact_endpoint_and_tail_contract',normalization_first_variation='2 j u(0)',normalization_is_not_optional=True,lowest_channel_adjoint='rho*(C_ell-2*m0*log(rho)+o(1))',finite_cut_terms=['+ell(epsilon)*Y_i(epsilon)','+ell(epsilon)*Z_ij(epsilon)'],thermal_cap_sum=encoded(6-arb.pi()**2/2),angular_allowance_polynomial=str(angular)))
    FIXTURES.append(dict(fixture='regular_endpoint_boundary_terms_and_uniform_thermal_cap',continuum_proof_in_note=True))

def actual_geometry_and_cap():
    par=prior.parameters()
    old=json.loads((ROOT/'evidence/current-bolt-adjoint-domain.json').read_bytes())['runs'][-1]
    h=raw.decode(next(x for x in old['numerics'] if x['kind']=='actual_current_bolt_normalization')['h'])
    beta,u,e=par['beta'],par['u'],par['eps'];v=1+h;chi=hull(arb(0),arb('.001'))
    Kv=(v+1)*(u+3*beta*(v*v+1))/(2*v**3)-e*(v-3)/(2*v**5)
    yy=2*(chi/2).sinh()**2; yy=hull(max(arb(0),yy.lower()),yy.upper());t=1+h*(1+yy)
    AA=3*beta+u+e
    H=beta*(t+2*v)+AA/(t*v*v)-2*e*(2*t+v)/(t*t*v**3)+e*(3*t*t+2*t*v+v*v)/(t**3*v**4)
    Ht=beta-AA/(t*t*v*v)+4*e/(t*t*v**3)+4*e/(t**3*v*v)-3*e/(t*t*v**4)-4*e/(t**3*v**3)-3*e/(t**4*v*v)
    B=2*Kv*v+yy*H;g=((yy+2)/(t*B)).sqrt()
    apchi=(chi/2).cosh()*B.sqrt()/(arb(2).sqrt()*Kv*t*t.sqrt())
    Bprime=chi.sinh()*(H+yy*h*Ht);tprime=h*chi.sinh()
    apchi+=arb(2).sqrt()*(chi/2).sinh()*(Bprime/(2*Kv*t*t.sqrt()*B.sqrt())-3*B.sqrt()*tprime/(2*Kv*t*t*t.sqrt()))
    aproper=apchi/g
    sinch=hull(arb(1),(chi.abs_upper()/2).cosh().upper())
    achi=sinch*B.sqrt()/(arb(2).sqrt()*Kv*t*t.sqrt());arho=achi/g
    check('actual_whole_cap_geometry',g>arb('.5') and g<2 and arho>arb('.5') and arho<2 and aproper.abs_upper()<2)
    epsilon=arb('5e-9');M=arb(26);diff,restore,total=cap_budgets(M,epsilon)
    check('cap_inside_actual_geometry',epsilon<arb('.001')*g.lower())
    check('combined_tip_precision',total<arb('1e-10'))
    budgets=[]
    for n in [1,2]:
        mass=(n*(n+par['q'])/par['C']).sqrt();assert mass<26
        budgets.append(dict(angular_n=n,mass=encoded(mass),budgets=[encoded(x) for x in cap_budgets(mass,epsilon)]))
    NUMERICS.append(dict(kind='actual_seed_endpoint_cap',h=encoded(h),chi_cover=['0','.001'],proper_over_chi=encoded(g),a_over_rho=encoded(arho),a_proper_derivative=encoded(aproper),epsilon=encoded(epsilon),mass_upper='26',difference_cap=encoded(diff),thermal_restoration_cap=encoded(restore),combined_cap=encoded(total),per_angular_level=budgets,units='unweighted paired-frequency action per angular state; all thermal j',includes_finite_cut_adjoint_terms=False,full_angular_remainder_certified=False))

def independent_bolt_diagnostic():
    M=mp.mpf(20);rho=mp.mpf('5e-9')
    results=[]
    NUMERICS.append(dict(kind='independent_flat_polar_regular_solution',cases=results,not_a_finite_geometry_operator=True,not_an_interval_certificate=True))
    for j in [0,1,8,100]:
        omega=mp.mpf(j)+mp.mpf('.5')
        def action(b,normalized=True):
            mass=M*mp.exp(b);x=mass*rho
            logR=mp.log(rho)/2+mp.log(mp.besseli(j,x)**2+mp.besseli(j+1,x)**2)/2
            if normalized:logR+=mp.loggamma(j+1)+j*mp.log(2/mass)
            root=mp.sqrt(omega**2+x*x)
            phi=root-omega+omega*mp.log(2*omega*rho/(omega+root))
            return -2*logR+2*phi
        val=action(0);grad=mp.diff(action,0);hess=mp.diff(action,0,2)
        lost=grad-mp.diff(lambda b:action(b,False),0)
        cap=mp.mpf(1550)*(M*rho)**2/(omega*(1+omega)**2)
        results.append(dict(j=j,subtracted_action=mp.nstr(val,45),gradient=mp.nstr(grad,45),hessian=mp.nstr(hess,45),omitted_normalizer_gradient=mp.nstr(lost,45),comparison_cap=mp.nstr(cap,45)))
        check('independent_regular_bolt_'+str(j),abs(grad)<cap and abs(hess)<cap and abs(lost-2*j)<mp.mpf('1e-45'))
        leading=M*M*rho*rho/(2*omega*(1+2*omega)**2)
        check('regular_bolt_leading_coefficient_'+str(j),abs(val/leading-1)<mp.mpf('1e-10'))
    FIXTURES.append(dict(fixture='actual_finite_seed_cap_and_independent_regular_bolt',all_thermal_cap_included=True,physical_correction=False))

def independent_fermi(y):
    y=mp.mpf(y);end=max(mp.mpf(2),mp.log(160/(mp.pi*y)))
    def f(v):return 1/(mp.exp(2*mp.pi*y*mp.cosh(v))+1)
    s0=mp.quad(f,[0,end/4,end/2,3*end/4,end])
    s1=mp.quad(lambda v:2*mp.pi*y*mp.cosh(v)*f(v)*(1-f(v)),[0,end/4,end/2,3*end/4,end])
    value=-4*y*y*mp.quad(lambda v:mp.sinh(v)**2*f(v),[0,end/4,end/2,3*end/4,end])
    return [value,4*y*y*s0,4*y*y*(2*s0-s1)]

def restoration_execution():
    cases=['1e-40','1e-8','.01','.0625','.125','.25','.5','1','4','20']
    for y in cases:
        vals,meta=restoration_dimensionless(arb(y)); ref=independent_fermi(y)
        NUMERICS.append(dict(kind='full_local_thermal_restoration',y=y,dimensionless_action_gradient_hessian=[encoded(x) for x in vals],algorithm=meta,independent_Fermi_values=[mp.nstr(x,55) for x in ref],independent_quadrature_is_not_certificate=True,independent_comparison_absolute_tolerance='1e-45',all_thermal_frequencies=True,physical_vacuum_restoration=False))
        check('restoration_width_'+y,all(x.upper()-x.lower()<arb('1e-12') for x in vals))
        # mpmath is an independent diagnostic, not a rigorous quadrature.
        check('independent_Fermi_integral_'+y,all(x.overlaps(arb(mp.nstr(r,55),'1e-45')) for x,r in zip(vals,ref)))
    for y in ['.0625','.1','.125']:
        aa=small_restoration(arb(y));bb=bessel_restoration(arb(y))
        check('independent_series_overlap_'+y,all(x.overlaps(z) for x,z in zip(aa[:3],bb[:3])))
        NUMERICS.append(dict(kind='two_rigorous_restoration_series_overlap',y=y,small=[encoded(x) for x in aa[:3]],Bessel=[encoded(x) for x in bb[:3]]))
    FIXTURES.append(dict(fixture='complete_flat_thermal_restoration_and_independent_Fermi_integral',analytic_tails=True,common_angular_restoration=False))

def main():
    ctx.prec=256;mp.mp.dps=60
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};deps=dependencies();error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_thermal_replay_accepted',old['status']=='passed' and old['complete_thermal_difference_sum_enclosed'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        algebra();actual_geometry_and_cap();independent_bolt_diagnostic();restoration_execution()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths};depafter=dependencies();status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),canonical_mass_endpoint_contract_executed=status=='passed',all_thermal_tip_allowance_enclosed=status=='passed',complete_local_thermal_restoration_implemented=status=='passed',full_thermal_angular_sum_evaluated=False,complete_renormalized_operator=False,physical_restorations_complete=False,actual_physical_correction_computed=False,raw_C5_admission_certified=False,nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};receipt['runs'].append(run)
    out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
    if status!='passed':raise SystemExit(1)

if __name__=='__main__':main()
