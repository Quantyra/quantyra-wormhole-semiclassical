"""Frozen arithmetic and independent mechanism checks; no continuum certification."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as F
from math import factorial, comb
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='finite-outer-admission-and-response'
PREVIOUS='construction-feasibility-audit'
OWN=['checks/finite_outer_admission_and_response.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
ROWS=[]; FIXTURES=[]
PROTOCOL={'digits':[60,90], 'identity_tolerance':'1e-40',
          'slab_relative_tolerance':'1e-20', 'cauchy_nodes':[32,64],
          'slab_mass_frequency':[[18,2],[24,3]],
          'actual_quantum_state_sampled':False,'continuum_proof_machine_verified':False}

def check(name,ok,detail=''):
    assert name not in {x['name'] for x in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)

def saved(stem):return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})
def mf(x):return mp.mpf(x.numerator)/x.denominator if isinstance(x,F) else mp.mpf(x)

def budgets():
    old=saved(PREVIOUS)
    check('previous_feasibility_audit_preserved',old['status']=='passed' and old['current_uniform_allowance_route_parked'])
    check('previous_source_hashes_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    t=S.symbols('t',real=True);cut=1-10*t**3+15*t**4-6*t**5
    check('quintic_C2_joins',all(S.diff(cut,t,j).subs(t,k)==(1-k if j==0 else 0) for j in range(3) for k in [0,1]))
    check('quintic_first_derivative_max',S.factor(S.diff(cut,t))==-30*t*t*(t-1)**2 and F(30,16)<2)
    # Second derivative extrema occur at (3+-sqrt(3))/6, giving 10/sqrt(3).
    check('quintic_second_derivative_max',S.simplify(S.diff(cut,t,2).subs(t,(3-S.sqrt(3))/6))==-10/S.sqrt(3) and F(100,3)<36)
    check('three_cutoff_gradient_budgets',F(30,16)/F('.5')<4 and F(30,16)/F('.25')<8)
    z=F('.063');bp=2*z/1000+6*F('2.01')*z*z+4*F('3.01')*z**3
    check('seed_polynomial_derivative_majorant',bp<F('.052'))
    rem=F('.075')*((F('1.001')**3-1)+F('.0011'))/F('.9989')+F('.052')*F('.412')**4*F('.15')/(8*F('.9989')**2)
    check('fixed_mass_direction_remainder_majorant',rem<F('.001'))
    leib=[]
    for j in range(6):
        # d^j(Y^5 H)/(j!100^j Y^(4-j)L), H derivative/k! <=Hmax20^k.
        bound=F('.001')/F('.15')**2 * F('.1')/3 /100**j * sum(F(comb(5,k))*2**(j-k) for k in range(min(j,5)+1))
        leib.append(bound)
        check('analytic_psi_weighted_remainder_'+str(j),bound<F('.002'))
    check('mass_coefficient_guard',F('1e-11')/F('.4')**3<F('1.6e-10'))
    check('finite_b_value_guard',2*F('3e-8')+F('1.2e-8')==F('7.2e-8'))
    check('finite_b_gradient_guard',13*F('3e-8')+F('7.7e-6')*F('1.00000002')+2*F('1.3e-8')<F('8.2e-6'))
    capglob=F('8.2e-6')+5*F('7.2e-8')
    check('actual_matching_cap_global_C1',capglob<F('8.6e-6'))
    check('mass_removed_cap_global_C1',F('8.6e-6')+100*F('1.6e-10')<F('1e-5'))
    v=F('1e-7')+F('.002')*F('1.6e-10');candidate=v+F('1.6e-10')
    check('actual_outer_coefficient_admission',v<F('1.000004e-7') and candidate<F('1.01e-7'))
    check('global_identification_first_derivative',2*F('4e-8')+F('1.3e-8')<F('1e-7'))
    check('global_identification_second_derivative',6*F('4e-8')+4*F('1.3e-8')+F('7.8e-6')<F('8.2e-6'))
    radial=(F('7.7e-6')+F('8.2e-6')/(1-F('1e-7')))/(1-F('1e-7'))
    angular=(F('7.7e-6')+2*F('1e-7'))/(1-F('1e-7'))+13*F('4e-8')
    check('actual_identification_log_derivatives',radial<F('1.7e-5') and angular<F('9e-6'))
    check('actual_global_log_sum',F('1.7e-5')+F('9e-6')+F('1.2e-7')+F('1e-7')<F('1e-4'))
    bg=F('1.5e-7');r=F('1e-8');caprho=F(1,40000);small=F('1e-7')
    cap1=F('1e25')*(bg+caprho)/caprho
    cap2=F('1e25')*(bg+2*caprho)/caprho**2
    check('large_cap_complex_tube',bg+2*caprho==F('.00005015') and bg+2*caprho<F('1e-4'))
    check('cap_Cauchy_coefficients',cap1==F('1.006e25') and cap2==F('8.024e29'))
    delta=10000*(bg+2*small)
    weighted=10*delta*F(4,3)+delta/3
    check('small_actual_graph_tube',delta==F('.0035') and 10*delta==F('.035'))
    check('weighted_actual_Neumann_margin',weighted<14*delta==F('.049') and 14*delta<F('.05'))
    check('weighted_inverse_below_two_over_M',F(4,3)/(1-weighted)<2)
    check('geometric_separation_and_exponent',F('1.1')-F('.05')==F('1.05') and 2*F('1.05')/5>F('.4'))
    with mp.workdps(60):
        check('point_evaluation_factor_reserve',34000/(2*mp.sqrt(mp.pi))*2*mp.sqrt(2)*mp.exp(mp.mpf('.05'))<mp.mpf('3e4'))
    check('trace_rank_commutator_density',2*8*F('1.1')*F('3e5')**2<F('2e12'))
    moment=24/F('.4')**3
    glued=F('2e12')*moment; norm=2*F('8.5e-33')
    g1=norm*glued/small;g2=norm*glued/small**2
    check('complete_angular_moment_and_glue',moment==375 and glued==F('7.5e14'))
    check('normalized_glue_Cauchy',norm*glued==F('1.275e-17') and g1==F('1.275e-10') and g2==F('.001275'))
    first=norm*cap1+g1;second=norm*cap2+g2
    check('full_scalar_first_response_budget',first==F('1.711475e-7') and first+F('5e-26')<F('1.8e-7'))
    check('full_scalar_mixed_response_budget',second==F('.0149158') and second+F('4e-27')<F('.02'))
    check('LLL_first_allowance',F('1e-31')*(2*F('2e5')+4*F('2.001')*100)<F('5e-26'))
    check('LLL_mixed_allowance',F('1e-31')*4*100**2==F('4e-27'))
    check('positive_outer_background_margin',bg-F('1.01e-7')==F('4.9e-8') and bg-F('1.01e-7')-r==F('3.9e-8'))
    check('hypothetical_outer_correction_source_targets',F('1.8e-7')*r==F('1.8e-15') and F('.02')*r*r/2==F('1e-18'))
    FIXTURES.append({'fixture':'outer_admission_and_response_budgets','candidate_coefficient_upper':str(candidate),
                    'fixed_mass_remainder':str(rem),'weighted_psi_bounds':list(map(str,leib)),
                    'first_response_before_LLL':str(first),'mixed_response_before_LLL':str(second),
                    'admission_margin':'4.9e-8','correction_norm_assumed':'1e-8',
                    'correction_norm_proved':False,'full_metric_ball_proved':False})

def symbolic_identities():
    x=S.symbols('x');a=S.Function('a')(x);b=S.Function('b')(x);w=S.Function('w')(x)
    density=S.sqrt(a*b)
    transformed=S.simplify(density/a*(S.diff(w/density,x)+S.diff(b,x)/(2*b)*w/density))
    check('half_density_Dirac_radial_identity',S.simplify(transformed-(S.diff(w,x)/a-S.diff(a,x)*w/(2*a*a)))==0)
    l=S.Function('l')(x);u=S.Function('u')(x);v=S.Function('v')(x);s,t=S.symbols('s t')
    warp=S.exp(-l-s*u-t*v);curv=S.simplify(-2*S.diff(warp,x,2)/warp)
    check('proper_scalar_curvature_identity',S.simplify(curv-2*S.diff(l+s*u+t*v,x,2)+2*S.diff(l+s*u+t*v,x)**2)==0)
    check('LLL_mixed_curvature_identity',S.simplify(S.diff(curv,s,t)+4*S.diff(u,x)*S.diff(v,x))==0)
    E,alpha=S.symbols('E alpha',positive=True)
    check('Stieltjes_mass_moment_identity',S.integrate(4*alpha*E**3*S.exp(-alpha*E),(E,0,S.oo))==24/alpha**3)
    # Independent formal solution of z'=d sqrt(1+B2 z^2+B3 z^3+...) through Y^4.
    Y,d,B2,B3=S.symbols('Y d B2 B3');zz=d*Y+B2*d**3*Y**3/6+B3*d**4*Y**4/8
    ode=S.series(S.diff(zz,Y)-d*S.sqrt(1+B2*zz**2+B3*zz**3),Y,0,4).removeO()
    check('mass_seed_proper_series_equation',S.expand(ode)==0)
    psifirst=S.series(zz**3/(1+B2*zz**2+B3*zz**3),Y,0,5).removeO()
    check('mass_seed_cubic_and_absent_quartic',S.expand(psifirst)==d**3*Y**3)
    FIXTURES.append({'fixture':'independent_symbolic_geometry','half_density':str(transformed),
                    'curvature':str(curv),'mass_cubic':str(psifirst),'actual_boundary_realization_proved_by_checker':False})

def finite_proper_pullback():
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            eps=mp.mpf('1e-11');ell=mp.mpf('.7')
            # Real C5, non-C6 at Y=1; logarithmic boundary jets are only open-collar C5.
            u=lambda y:eps*(y**3+y**4*mp.log(y)+max(y-1,mp.mpf(0))**mp.mpf('5.5')) if y else mp.mpf(0)
            ell0=lambda y:mp.mpf('.01')*y*y+mp.mpf('.005')*y**3/(1+y)
            def phi(y):return mp.quad(lambda t:mp.exp(ell*u(t)),[0,y] if y<=1 else [0,1,y])
            for y in map(mp.mpf,['.02','.1','1.3','2']):
                X=phi(y)
                def pulled(xx):
                    yy=mp.findroot(lambda z:phi(z)-xx,(y-mp.mpf('.001'),y+mp.mpf('.001')),solver='secant')
                    return ell0(yy)-ell0(xx)-ell*u(yy)
                value=ell0(y)-ell0(X)-ell*u(y)
                derivative=(mp.diff(ell0,y)-ell*mp.diff(u,y))*mp.exp(-ell*u(y))-mp.diff(ell0,X)
                h=mp.mpf('1e-9')
                fd=(-pulled(X+2*h)+8*pulled(X+h)-8*pulled(X-h)+pulled(X-2*h))/(12*h)
                err=abs(fd-derivative)
                check(f'fixed_proper_pullback_{dps}_{y}',abs(pulled(X)-value)<mp.mpf('1e-40') and err<mp.mpf('1e-38'))
                records.append({'dps':dps,'Y':str(y),'X':mp.nstr(X,25),'derivative_error':mp.nstr(err,15)})
    FIXTURES.append({'fixture':'real_C5_non_C6_proper_pullback','records':records,'actual_charged_geometry_sampled':False})

def mass_direction():
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            x=mp.mpf('8e9');C=1+6*x;a=1+2*x;q=1+3*x
            B=lambda z,p=0:x+z*z-2*(a+x*p)*z**3+q*z**4
            proper=lambda z,p:mp.sqrt(C)*mp.quad(lambda t:1/mp.sqrt(B(t,p)),[0,z])
            for Y in map(mp.mpf,['.02','.08','.1']):
                z=mp.findroot(lambda z:proper(z,0)-Y,(Y*mp.mpf('.4'),Y*mp.mpf('.42')))
                exact=x*z**3/B(z)+x*mp.diff(B,z)/(2*mp.sqrt(B(z)))*mp.quad(lambda t:t**3/B(t)**mp.mpf('1.5'),[0,z])
                def fixed_logf(p):
                    zp=mp.findroot(lambda zz:proper(zz,p)-Y,(z*mp.mpf('.99'),z*mp.mpf('1.01')))
                    return mp.log(x/B(zp,p))/2
                h=mp.mpf('1e-8')
                fd=(-fixed_logf(2*h)+8*fixed_logf(h)-8*fixed_logf(-h)+fixed_logf(-2*h))/(12*h)
                err=abs(fd-exact)/abs(exact)
                check(f'fixed_proper_mass_variation_{dps}_{Y}',err<mp.mpf('1e-30'))
                records.append({'dps':dps,'Y':str(Y),'relative_error':mp.nstr(err,15),'psi':mp.nstr(exact,25)})
    FIXTURES.append({'fixture':'independent_mass_direction_variation','records':records,'benchmark_thermal_geometry_sampled':False})

def dirac_slab():
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            I=mp.eye(2);g1=mp.diag([1,-1]);g2=mp.matrix([[0,1],[1,0]])
            aa=mp.mpf('1.1');width=mp.mpf('.25')
            for mass,omega in PROTOCOL['slab_mass_frequency']:
                M=mp.mpf(mass);om=mp.mpf(omega);mu=mp.sqrt(M*M+om*om);k=om/(M+mu)
                H=lambda f:-g1*(M*I+1j*om*f*g2)
                def transfer(f,distance):
                    lam=mp.sqrt(M*M+(om*f)**2)
                    return mp.cosh(lam*distance)*I+mp.sinh(lam*distance)*H(f)/lam
                left=mp.matrix([-1j*k,1]);right=mp.matrix([1,-1j*k])
                def green(delta):
                    rr=transfer(1,-aa)*transfer(1+delta,-width)*right
                    mat=mp.matrix([[rr[0],-left[0]],[rr[1],-left[1]]]);cc=mat**-1*g1
                    return (rr*cc[0,:]+left*cc[1,:])/2
                tr=lambda mat:mat[0,0]+mat[1,1]
                scalar=lambda dd:tr(green(dd))
                expected=(M*I-1j*om*g2)/(2*mu)
                check(f'slab_free_Green_{dps}_{mass}',mp.norm(green(0)-expected)<mp.mpf('1e-40'))
                sp=(M*I-1j*om*g2+mu*g1)/(2*mu)
                sm=(M*I-1j*om*g2-mu*g1)/(2*mu)
                born=-tr(sm*(1j*om*g2)*sp)*mp.quad(lambda y:mp.exp(-2*mu*y),[aa,aa+width])
                exact=mp.diff(scalar,mp.mpf(0))
                check(f'slab_resolvent_first_variation_{dps}_{mass}',abs(exact-born)/abs(born)<mp.mpf('1e-20'))
                center=mp.mpf('.001');radius=mp.mpf('.01');second=mp.diff(scalar,center,2);errors=[]
                base=scalar(center)
                for count in PROTOCOL['cauchy_nodes']:
                    coefficient=sum((scalar(center+radius*mp.exp(2j*mp.pi*j/count))-base)*mp.exp(-4j*mp.pi*j/count) for j in range(count))/count
                    estimate=2*coefficient/radius**2;error=abs(estimate-second)/abs(second);errors.append(error)
                    check(f'slab_Cauchy_second_{dps}_{mass}_{count}',error<mp.mpf('1e-20'))
                # A local cap cannot simply inherit a distance-24 exponent.
                false_allowance=mp.mpf('1e10')*mp.exp(-48*mu)
                check(f'false_long_separation_rejected_{dps}_{mass}',abs(born)>mp.mpf('1e100')*false_allowance)
                records.append({'dps':dps,'mass':mass,'omega':omega,
                                'first_variation':mp.nstr(born,25),'second_variation':mp.nstr(second,25),
                                'Cauchy_relative_errors':[mp.nstr(e,15) for e in errors],
                                'false_distance24_allowance':mp.nstr(false_allowance,15)})
    FIXTURES.append({'fixture':'exact_nonlocal_Dirac_slab','records':records,
                    'actual_quantum_state_sampled':False,'original_distance24_theorem_refuted':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 finite outer admission and response','input_count':len(ps)}),flush=True)
    try:budgets();symbolic_identities();finite_proper_pullback();mass_direction();dirac_slab()
    except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
    run={'story':'S167 finite outer admission and response','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
         'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after}
    for k,v in review().items():
        if isinstance(v,bool) and k not in run:run[k]=v
    for kind,p in zip(['checker','note','review'],OWN):run['frozen_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
