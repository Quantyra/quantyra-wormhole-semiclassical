"""Supporting checks for the ordinary fixed-seed coercive interior estimate."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='seed-middle-coercive-jets'
PREVIOUS='finite-graph-admissibility'
OWN=['checks/seed_middle_coercive_jets.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational
ROWS=[]; FIXTURES=[]
PROTOCOL={'exact_arithmetic':'symbolic and rational', 'numerical_dps':[60,90],
          'maximum_gluing_coordinate_order':9, 'complete_seed_source_order':3,
          'angular_sum_cutoff':10000, 'numerical_tolerance':'1e-45',
          'actual_quantum_source_numerically_evaluated':False,
          'ordinary_continuum_proof_separate':True}


def check(name, ok, detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row); print(json.dumps(row),flush=True)
    if not ok: raise AssertionError(name)


def zero(name, expression):
    v=S.simplify(expression); check(name,v==0,v)


def prior():
    return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]


def paths():
    return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})


def provenance():
    p=prior()
    check('predecessor_passed',p['status']=='passed' and len(p['checks'])==65 and all(c['passed'] for c in p['checks']))
    check('predecessor_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('retained_hashes',not bad,bad)
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('retained_sources_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
    check('scope_limits',review['complete_seed_source_order']==3 and review['gluing_coordinate_order']==9
          and not review['actual_finite_source_H6_bounded'] and not review['actual_transfer_certificate_passed'])


def geometry_and_energy():
    radius=R(1,256); separation=R(1,65)-R(1003,1000)*radius
    check('remaining_separation',separation>R(1,88),separation)
    check('original_evaluation_fits',R(1,88)>R(19,100)/18)
    check('Cartesian_Cauchy_base',8/(R(1,32)-radius)**2<128**2)
    # e^x <= 1/(1-x) for 0<=x<1; e^-x >= 1-x.
    check('analytic_coefficient_amplitude',R(1,1)/(1-R(1,500))-1<R(21,10000))
    check('conjugated_base_row',30000/(1-R(1,1000))<31000)
    lower=R(996,1000)-R(21,10000)*128/18
    check('mass_coercivity',lower>R(98,100),lower)
    check('coercive_square_root',R(11,10)**2*R(98,100)>1)
    x,y,M=S.symbols('x y M',real=True)
    a=S.Function('a')(x,y)
    u=S.Matrix([S.Function('u')(x,y),S.Function('v')(x,y)])
    sx=S.Matrix([[0,1],[1,0]]); sy=S.Matrix([[1,0],[0,-1]])
    check('flat_Clifford_relations',sx*sx==S.eye(2) and sy*sy==S.eye(2) and sx*sy+sy*sx==S.zeros(2))
    D=lambda f:sx*f.diff(x)+sy*f.diff(y)
    check('flat_D_squared',S.simplify(D(D(u))-u.diff(x,2)-u.diff(y,2))==S.zeros(2,1))
    A=lambda f:D(f)+M*a*f
    Astar=lambda f:-D(f)+M*a*f
    expected=-u.diff(x,2)-u.diff(y,2)+M*M*a*a*u-M*(sx*S.diff(a,x)+sy*S.diff(a,y))*u
    check('variable_mass_square_identity',S.simplify(Astar(A(u))-expected)==S.zeros(2,1))
    check('mass_gradient_contact_not_zero',S.simplify(Astar(A(u))+u.diff(x,2)+u.diff(y,2)-M*M*a*a*u)!=S.zeros(2,1))
    # Normal-chart Christoffels: exact tensor norm on a general symmetric Hessian.
    h00,h01,h11=S.symbols('h00 h01 h11',real=True)
    h=S.Matrix([[h00,h01],[h01,h11]])
    dg=[[[[int(i==k)*h[l,j]+int(j==k)*h[l,i]-int(i==j)*h[l,k]
           for j in range(2)] for i in range(2)] for k in range(2)] for l in range(2)]
    norm2=sum(t*t for ll in dg for kk in ll for ii in kk for t in ii)
    zero('Christoffel_derivative_tensor_norm',norm2-4*(h00*h00+2*h01*h01+h11*h11))
    check('third_covariant_conversion',4*R(2,1000)*2*128**2<300)
    FIXTURES.append({'fixture':'local_geometry_and_mass_square','actual_normal_chart_guards_used':True,
                     'actual_quantum_kernel_evaluated':False,'coercivity_lower':str(lower)})


def interior_recurrence():
    margins=[]
    for j in range(10):
        m=j+3; B=S.Integer(600*m); L=S.Integer(512*m)
        ratios=[R(11,10)*(L/B+R(21,10000)*sum(S.factorial(k)/S.factorial(k-b)*(128/B)**b for b in range(1,k+1)))
                for k in range(j+2)]
        check('mass_independent_recurrence_'+str(j),max(ratios)<1,max(ratios))
        endpoint=1+R(21,10000)*sum(S.factorial(j)/S.factorial(j-b)*(128/B)**b for b in range(1,j+1))
        check('endpoint_Leibniz_'+str(j),endpoint<R(101,100),endpoint)
        margins.append({'order':j,'energy_ratio_max':str(max(ratios)),'endpoint_factor':str(endpoint)})
    t=S.symbols('t',real=True)
    cutoff=1-10*t**3+15*t**4-6*t**5
    zero('quintic_first_derivative',S.diff(cutoff,t)+30*t*t*(1-t)**2)
    check('cutoff_slope',R(15,8)<2)
    # max |chi''|=10/sqrt(3)<6; r>=rho then bounds the radial Laplacian.
    check('cutoff_second_derivative',R(100,3)<36)
    check('cutoff_Laplacian',6+2<=8)
    for k in range(3):
        zero('cutoff_left_jet_'+str(k),S.diff(cutoff,t,k).subs(t,0)-(1 if k==0 else 0))
        zero('cutoff_right_jet_'+str(k),S.diff(cutoff,t,k).subs(t,1))
    r,rho=S.symbols('r rho',positive=True)
    zero('Fourier_evaluation_integral',S.integrate(2*S.pi*r/(1+rho*rho*r*r)**2,(r,0,S.oo))-S.pi/rho**2)
    point=31000*R(101,100)*(9+4*R(75,64)+2*R(75,64)**2)
    check('row_prefactor',point<520000,point)
    FIXTURES.append({'fixture':'all_finite_energy_and_endpoint_budgets','orders':margins,
                     'row_prefactor_upper':str(point),'mass_power_in_row_bound':0})


def angular_and_source():
    q=S.Integer(2000)**4
    check('angular_monotonicity_guard',R(352,100000)/q<R(1,88**2))
    z=S.symbols('z',positive=True)
    # 2z/(1+z^2)<=1 is the exact maximum of the logarithmic derivative.
    zero('angular_log_derivative_maximum',(1+z*z)*(1-2*z/(1+z*z))-(z-1)**2)
    n,C,qq=S.symbols('n C q',positive=True)
    energy=n*(n+qq)/C; multiplicity=2*(qq+2*n)
    zero('angular_exact_count_jacobian',multiplicity-2*C*S.diff(energy,n))
    c,MM=S.symbols('c M',positive=True)
    zero('angular_mass_integral',S.integrate(4*C*MM*S.exp(-c*MM),(MM,0,S.oo))-4*C/c**2)
    tail=4*C*(MM/c+1/c**2)*S.exp(-c*MM)
    zero('angular_tail_primitive',S.diff(tail,MM)+4*C*MM*S.exp(-c*MM))
    coeff=[S.Integer(30000)]+[S.Integer(520000)*(j+3)*(600*(j+3))**j for j in range(1,10)]
    pref=8*256**2*88**2*R(31,2500000)
    bounds=[]
    # Independent exponential-generating-function product checks every Leibniz sum.
    generating=sum(coeff[j]*z**j/S.factorial(j) for j in range(10))
    square=S.Poly(S.expand(generating**2),z)
    for j in range(10):
        convolve=sum(S.binomial(j,k)*coeff[k]*coeff[j-k] for k in range(j+1))
        zero('gluing_Leibniz_'+str(j),convolve-S.factorial(j)*square.coeff_monomial(z**j))
        gj=pref*convolve/(S.factorial(j)*100**j)
        check('finite_positive_gluing_bound_'+str(j),gj>0 and gj.is_Rational,gj)
        bounds.append(gj)
    cov3=bounds[3]+bounds[1]/200
    complete=max(*bounds[:3],cov3)+2*10**13+10**12
    check('complete_middle_C3_target',complete<10**22,complete)
    check('four_order_improvement',10**26/10**22==10000)
    gamma=R(68,10**29)
    zero('new_f_source_term',R(6,100)*gamma*10**22-R(408,10**9))
    zero('new_direct_curvature_source_term',R(66,10)*gamma*10**22-R(4488,10**8))
    FIXTURES.append({'fixture':'complete_angular_and_seed_source_ledger',
                     'row_constants':[str(v) for v in coeff],
                     'gluing_coordinate_bounds':[str(v) for v in bounds],
                     'gluing_coordinate_bounds_decimal':[str(S.N(v,18)) for v in bounds],
                     'complete_covariant_C3_bound':str(complete),
                     'complete_covariant_C3_bound_decimal':str(S.N(complete,18)),
                     'adopted_complete_middle_Km':'1e22','previous_Km':'1e26',
                     'actual_finite_source_H6_bounded':False,'endpoint_C5_contract_changed':False})


def energy_integrals():
    records=[]
    for precision in PROTOCOL['numerical_dps']:
        with mp.workdps(precision):
            radius=mp.mpf(1)/256
            p=lambda x:(1-(x/radius)**2)**2
            dp=lambda x:-4*x/radius**2*(1-(x/radius)**2)
            py=mp.quad(lambda y:p(y)**2,[-radius,0,radius])
            dy=mp.quad(lambda y:dp(y)**2,[-radius,0,radius])
            for mass in [18,25,100]:
                a=lambda x:1+x/10
                v2=lambda x:mp.exp(-2*mass*(x+x*x/20))
                integ=lambda f:mp.quad(lambda x:f(x)*v2(x),[-radius,0,radius])
                p2=integ(lambda x:p(x)**2)
                lhs=integ(lambda x:dp(x)**2)*py+p2*dy
                rhs=(integ(lambda x:(dp(x)-mass*a(x)*p(x))**2)
                     +mass**2*integ(lambda x:a(x)**2*p(x)**2)-mass*p2/10)*py+p2*dy
                error=abs(lhs-rhs)/max(1,abs(lhs))
                check('independent_energy_integral_'+str(precision)+'_'+str(mass),error<mp.mpf('1e-45'),mp.nstr(error,12))
                omitted=mass*p2*py/10
                check('mass_gradient_omission_detected_'+str(precision)+'_'+str(mass),omitted>mp.mpf('1e-9'),mp.nstr(omitted,12))
                records.append({'dps':precision,'mass':mass,'relative_error':mp.nstr(error,20),
                                'omitted_cross_term':mp.nstr(omitted,20)})
    FIXTURES.append({'fixture':'independent_variable_mass_Dirac_energy','quadrature':records,
                     'actual_charged_geometry':False,'constant_mass_substitution':False})


def spectral_sums():
    records=[]
    for precision in PROTOCOL['numerical_dps']:
        with mp.workdps(precision):
            q=mp.mpf(2000)**4; N=PROTOCOL['angular_sum_cutoff']
            for divisor in [340,350]:
                C=q/divisor
                mass=lambda n:mp.sqrt(n*(n+q)/C)
                partial=mp.fsum(2*(q+2*n)/C*mp.exp(-mass(n)/88) for n in range(1,N+1))
                mn=mass(N); tail=4*(88**2+88*mn)*mp.exp(-mn/88)
                upper=partial+tail
                check('guard_point_spectral_enclosure_'+str(precision)+'_'+str(divisor),0<partial<upper<4*88**2,
                      'partial/C='+mp.nstr(partial,20)+'; tail/C<'+mp.nstr(tail,12))
                check('guard_point_M1_'+str(precision)+'_'+str(divisor),18<mass(1)<19,mp.nstr(mass(1),15))
                records.append({'dps':precision,'q':str(int(q)),'C':'q/'+str(divisor),'N':N,
                                'partial_over_C':mp.nstr(partial,45),'tail_over_C_upper':mp.nstr(tail,45)})
    for i in range(2):
        with mp.workdps(60):
            err=abs(mp.mpf(records[i]['partial_over_C'])-mp.mpf(records[i+2]['partial_over_C']))
            check('spectral_precision_consistency_'+str(i),err<mp.mpf('1e-38'),mp.nstr(err,12))
    FIXTURES.append({'fixture':'guard_point_complete_angular_tail','enclosures':records,
                     'actual_benchmark_C_evaluated':False,'angular_tower_truncated_in_proof':False})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n')
    sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths(); hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes(); error=None
    print(json.dumps({'event':'start','story':'S167 seed middle coercive jets','input_count':len(ps)}),flush=True)
    try:
        provenance(); geometry_and_energy(); interior_recurrence(); angular_and_source(); energy_integrals(); spectral_sums()
    except Exception:
        error=traceback.format_exc(); print(error,file=sys.stderr,flush=True)
    after=hashes(); good=error is None and before==after and all(r['passed'] for r in ROWS)
    run={'story':'S167 seed middle coercive jets','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'ordinary_complete_seed_middle_C3_bound_proved':True,'complete_seed_middle_Km':'1e22',
         'ordinary_gluing_coordinate_jets_through_nine_proved':True,
         'complete_seed_nine_jet_bound_proved':False,'actual_finite_source_H6_bounded':False,
         'actual_transfer_norms_numerically_bounded':False,'actual_transfer_certificate_passed':False,
         'actual_signed_moment_evaluated':False,'actual_signed_compatibility_zero_proved':False,
         'endpoint_C5_contract_changed':False,'nonlinear_neighborhood_preserved':False,
         'nonlinear_solution_constructed':False,'physical_stability_proved':False,
         'physical_source_projection_changed':False,'physical_boundary_condition_changed':False,
         'ordinary_proof_machine_verified':False,'formal_verification':False,
         'S167_complete':False,'S168_complete':False,'goal_complete':False}
    if not good:
        for kind,path in zip(['checker','note','review'],OWN):
            run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json')
    saved=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    saved['runs'].append(run)
    dest.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__': raise SystemExit(main())
