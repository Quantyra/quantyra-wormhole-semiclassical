"""S167: finite trace gradient, repaired radius pullback and fixed-projection transversality.

Exact metric algebra and rational budgets; continuum claims remain ordinary proofs.
"""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "finite-trace-gradient-and-projection-transversality"
OWN = ["checks/finite_trace_gradient_and_projection_transversality.py", "notes/"+STEM+".md"]
PREDECESSORS = [("on-shell-adjoint-error-and-trace-gradient", 329),
                ("finite-comparison-scalar-curvature-remainder", 314),
                ("finite-comparison-weighted-mixed-form", 335),
                ("cold-reference-overlap-and-adjoint-budget", 316)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "original fixed projection transverse; useful range-response constant and nonlinear solve open",
    "coordinate_correction": "Delta_hat=alpha^-2 D2; L4=alpha^2(Delta_hat+2 grad phi dot grad)",
    "preserved_numeric_inputs": "D2nu<2.3e-6, D2S<2.31e-6 and L4nu<2.2e-6",
    "trace_gradient_upper": "6e-4 (sharper intermediate <5.7e-4); 3e-4 NOT established",
    "full_errors": "e0<.007, e_ell<.13",
    "full_cutoff_overlap_lower": ".78",
    "normalized_projection_overlap_lower": ".048",
    "functional_ratio_upper": "21.3",
    "projected_to_range_inverse_factor_upper": "22.5",
    "fixture_scope": "coordinate and derivative identities, rational budgets, independent warped and projection fixtures",
    "formal_verification": False,
}


def check(name, ok, detail=""):
    assert name not in {row["name"] for row in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expr):
    expr = S.cancel(expr.doit())
    check(name, expr == 0, expr)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_bytes())
    return data.get("runs", [data])[-1]


def inputs(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def paths():
    found = set(OWN)
    for stem, _ in PREDECESSORS:
        found.add("evidence/"+stem+".json")
        found.update(inputs(latest(stem)))
    return sorted(found)


def predecessors():
    expected = {}
    for stem, count in PREDECESSORS:
        run = latest(stem)
        check(stem+"_receipt", run["status"] == "passed" and len(run["checks"]) == count
              and all(row["passed"] for row in run["checks"]) and inputs(run) == inputs(run, True))
        conflicts = [p for p, h in inputs(run).items() if p in expected and expected[p] != h]
        check(stem+"_shared_inputs", not conflicts, conflicts)
        expected.update(inputs(run))
    for p, h in sorted(expected.items()):
        check("input_"+p, sha256((ROOT/p).read_bytes()).hexdigest() == h)
    for p in ["on-shell-adjoint-error-and-trace-gradient", "cold-reference-overlap-and-adjoint-budget",
              "finite-comparison-scalar-curvature-remainder", "finite-comparison-weighted-mixed-form",
              "finite-comparison-adjoint-and-fixed-projection", "stationary-completed-metric-response-norm",
              "finite-comparison-classical-disk-contact", "finite-comparison-classical-tt-form",
              "stationary-centered-covariant-seed", "stationary-centered-constraint-residual"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)
















def coordinate_repair():
    v1,v2,aa,rr,alpha=S.symbols('v1 v2 A R alpha',nonzero=True)
    dzz=(v2+rr*v1)/alpha**2;Ahat=(aa-rr)/alpha
    zero('canonical_scalar_laplacian_chain_rule',dzz+Ahat*v1/alpha-(v2+aa*v1)/alpha**2)
    zero('correct_four_dimensional_laplacian',alpha**2*(dzz+Ahat*v1/alpha+2*rr*v1/alpha**2)-(v2+(aa+2*rr)*v1))
    zero('omitted_chain_term_is_nonzero',v2/alpha**2+Ahat*v1/alpha-(v2+aa*v1)/alpha**2+rr*v1/alpha**2)
    check('old_missing_chain_term_detected',S.cancel(-rr*v1/alpha**2)!=0)
    K=S.symbols('K');Rprime=aa*rr-rr**2
    zero('radius_force_coefficient_correction',Rprime+aa*rr+rr**2-2*aa*rr)
    x,d,delta,t=S.symbols('x d delta t',positive=True)
    F=x*t**2+1-(2+4*x+delta)/t+(1+3*x+d)/t**2-2*d/t**3+d/t**4
    ext=F+delta/t
    zero('extremal_gap_derivative',S.diff(ext,t)-2*x*(t-1)*(t**3+t**2+t+3)/t**3-2*(t-1)/t**3+2*d*(t-1)*(t-2)/t**5)
    check('seed_mass_charge_negative_remainder',-1-R('7.9e9')+4*R('1e-12')<0)
    check('sharper_seed_radial_square',1+1/R('7.9e9')<R('1.001'))
    force=R('6.1')*R('5.2e-8')/6+(R('6.1e-20')+R('1.2e-8')*R('.056'))/6+4*R('1.2e-8')*11+R('1.2e-8')*R('12.01')/3
    check('correct_radius_force_budget',force<R('6.3e-7'),force)
    repaired=R('6.3e-7')+2*R('1.001')*(R('3e-7')+R('1.2e-8'))
    check('repaired_D2nu_sharper_bound',repaired<R('1.26e-6'),repaired)
    check('previous_D2nu_input_preserved',repaired<R('2.3e-6'))
    check('previous_L4nu_input_preserved',R('6.3e-7')<R('2.2e-6'))
    Sbound=(R('1.001')*repaired+R('2e-16')*22+R('1.001')*7*R('3e-7')**2)/R('.999')+R('7.1e-7')**2
    check('previous_D2S_input_preserved',Sbound<R('2.31e-6'),Sbound)
    check('actual_nu_second_derivative',R('1.26e-6')+R('3.432e-6')<R('5.7e-6'))
    check('gap_derivative_coefficient',12+2/R('7.9e9')+4*R('1.3e-22')+R('7.1e-30')<R('12.1'))


def low_and_third_budgets():
    check('V_first_derivative',7*R('3e-7')**2<R('8e-7')**2)
    check('relative_tau_first',7<R('2.65')**2)
    check('relative_tau_second',R('6.05')+14<21)
    v2=R('5.7e-6')+2*R('2.65')*R('8e-7')+21*R('1.2e-8')
    check('weighted_V_second',v2<R('1.1e-5'),v2)
    check('L4tau_derivative',5*(21+44)**2<146**2)
    check('core_n_weight',40*R('2.4e-8')<R('1e-6'))
    check('outer_n_weight',R('1e-8')/4<R('1e-6'))
    square=40*(R('2.65')*R('8e-7')**2+2*R('8e-7')*R('1.1e-5'))
    check('weighted_gradient_square_derivative',square<R('8e-10'),square)
    vtt=R('2.5e-7')*R(4,3)+R('1.4e-8')*R(16,9)
    check('outer_V_tt',vtt<R('3.6e-7'))
    check('outer_V_y',5*R('1.9e-8')**2<R('4.3e-8')**2)
    check('outer_V_yy',5*R('3.6e-7')+11*R('1.9e-8')<R('2.1e-6'))
    check('canonical_nu_gradient',R('7e-7')**2/R('5.9')<R('3e-7')**2)
    check('canonical_nu_hessian',(R('5.7e-6')+3*R('2.24')*R('7e-7'))/R('5.9')<R('3e-6'))
    check('cold_normal_third_jet',R('1e-5')*R('.00331')*6*100**3<R('.199'))
    tail=sum(R(1,10)**j for j in range(1,6))
    check('normal_jet_triangular_inverse',1/(1-tail)<2)
    bracket=24+40*R('3e-7')+2*5*R('3e-7')+2*R('3e-6')+2*(2*R('3e-6')+2*R('3e-7'))
    check('canonical_L4nu_derivative',bracket<R('24.001'))
    check('core_L4nu_derivative',R('2.5')**3*bracket<376<500)
    Y=R('.075');L=R('3.6')
    # L(.075)<3.6: exp(2.6)>40/3, by a positive Taylor lower bound.
    check('outer_log_endpoint',sum(R('2.6')**j/S.factorial(j) for j in range(30))>R(40,3))
    n1=3*R('1e-11')*R('1.1')*Y**3+R('2e-8')*100*R('1.1')*Y**4*L
    n2=(6*R('1.1')**2+9)*R('1e-11')*Y**3+R('2e-8')*(2*100**2*R('1.1')**2+300)*Y**4*L
    n3=139*R('1e-11')*Y**3+R('2e-8')*(6*100**3*R('1.1')**3+3*2*100**2*R('1.1')*3+2400)*Y**4*L
    check('outer_nu_third',n3<R('2e-5'),n3)
    check('outer_L4nu_derivative',n3+(2+2*R('2.24'))*n2+(15+32)*n1<R('.001')<500)


def nonlinear_derivative_budgets():
    U=R('1.2e-8');vy=R('8e-7');wy=R('1.03e-6');sy=R('7.1e-7');a=R('1.001')
    check('radius_small_coefficient',a*a*U<R('1.21e-8'))
    check('radius_small_coefficient_derivative',a*a*(vy+U*(R('2.24')+sy))<R('8.5e-7'))
    check('radius_B_derivative',a*U*vy<R('1e-14'))
    linearS=R('1.21e-8')*R('7e-7')+a*R('2e-16')*R('2.24')
    first=R('12.01')*2*a*(linearS+2*R('1e-6')*wy)+R('1e-30')
    check('value_remainder_derivative',first<R('1e-9'),first)
    dot=2*a*wy*R('6.2e-20')+R('6.1e-8')*R('3e-20')
    check('linear_trace_remainder_derivative',dot<R('1e-23'))
    main=6*a*(R('1.21e-8')*500+R('8.5e-7')*R('2.2e-6')+2*wy*R('1.21e-8')*R('2.2e-6'))
    check('dominant_L4nu_derivative_budget',main<R('3.64e-5'),main)
    tau=6*a*a*((R('1e-14')+sy*R('2e-16'))*22+R('2e-16')*146+2*wy*R('2e-16')*22)
    check('L4tau_derivative_budget',tau<R('1e-10'))
    square=6*a*a*(R('8e-10')+(2*wy+vy+sy)*40*vy**2)
    check('gradient_square_derivative_budget',square<R('5e-9'),square)
    end=2*a*(5*U*R('3.2e-7')+4*(wy+R('2.24'))*U**2)/R('7.9e9')
    check('final_exponential_derivative_budget',end<R('1e-20'))
    total=R('1e-9')+R('1e-23')+R('3.64e-5')+R('1e-10')+R('5e-9')+R('1e-20')
    check('physical_scalar_gradient_bound',a*(total+R('5e-20'))<R('4e-5'),total)


def fifth_jet_counts():
    q=S.symbols('q');K,K1,K2,K3=S.symbols('K K1 K2 K3')
    J=q-K*q**3/6-K1*q**4/12+(K*K/120-K2/40)*q**5+(K*K1/120-K3/180)*q**6
    curvature=K+K1*q+K2*q*q/2+K3*q**3/6
    defect=S.expand(S.diff(J,q,2)+curvature*J)
    check('Jacobi_normal_metric_coefficients',all(defect.coeff(q,j)==0 for j in range(5)))
    zero('normal_metric_fifth_coefficient',S.expand((J/q)**2).coeff(q,5)-(2*K*K1/45-K3/90))
    check('seed_normal_metric_fifth',5**5*(R(2,45)*2*10**6+R(1,90)*10**6)<10**10)
    check('seed_log_radius_fifth_count',2**5*S.factorial(5)*10**4<10**8)
    L1,L2,L3,L4,L5=10,1000,300000,10**8,2*10**10
    bell5=L5+5*L1*L4+10*L2*L3+10*L1**2*L3+15*L1*L2**2+10*L1**3*L2+L1**5
    check('log_fifth_reserve',R('.01')*S.factorial(5)*100**5+2*10**8< L5)
    check('exponential_fifth_Bell_polynomial',bell5<3*10**10,bell5)
    bell=[1,L1,L2+L1**2,L3+3*L1*L2+L1**3,L4+4*L1*L3+3*L2**2+6*L1**2*L2+L1**4,bell5]
    normal=[1,0,10,10000,10**7,10**10]
    check('metric_fifth_component_reserve',100*sum(S.binomial(5,j)*bell[j]*normal[5-j] for j in range(6))<10**16)
    M=[1,30,3000,10**6,4*10**8,10**16];N=[S.Integer(1)]
    for j in range(1,5):N.append(4*sum(S.binomial(j,p)*M[p]*N[j-p] for p in range(1,j+1)))
    G=[6*sum(S.binomial(j,p)*N[p]*M[j+1-p] for p in range(j+1)) for j in range(5)]
    Q=[2*G[j+1]+8*sum(S.binomial(j,p)*G[p]*G[j-p] for p in range(j+1)) for j in range(4)]
    sc=[64*sum(S.binomial(j,p)*N[p]*Q[j-p] for p in range(j+1)) for j in range(4)]
    B=32*(N[1]*sc[2]+sc[3]+4*G[1]*sc[1]+4*G[0]*sc[2]+4*N[1]*G[0]*sc[1])
    check('gradient_laplacian_component_count',B<R('2.5e20'),B)
    check('compact_fifth_weight',40*(R('2.5')**5*B)**2<R('1.6e23')**2)
    check('outer_fifth_product_count',S.factorial(5)**2*2**5<10**6)
    check('outer_fifth_metric_difference',10**6*R('1.2e-6')*S.factorial(5)*100**5/40<R('.001')*M[5])
    check('outer_telescoping_product_count',20*2**19<2*10**7)
    check('outer_fifth_weight',(2*10**4*R('2.5')**5*B)**2/40**3<R('2e24')**2)
    check('global_anomaly_gradient_reserve',R('2e24')+R('9e4')*R('1.3e-22')<R('1e28'))


def complete_trace_and_projection():
    sg=R('4e-17');r200=R('14.143')
    check('sqrt200_reserve',r200**2>200)
    squares=R(7,6)*50*10**15+R(11,3)*8*10**15+R(11,36)*R('3.1e-12')*R('4e-5')
    check('curvature_square_gradient_reserve',squares<10**17)
    trace=R('4e-5')*r200+sg**2/180*(10**28+10**17)*r200+R('1e-24')
    check('complete_finite_trace_gradient',trace<R('5.7e-4')<R('6e-4'),trace)
    err=R('6.7')*(R('30.1')*R('6e-4')+R('1e-8'))+R('.0064')
    check('full_finite_adjoint_error',err<R('.13'),err)
    e0=R('.007');el=R('.13')
    overlap=R('1.5')-(R('.4')+R('1.1')*e0)**2-4*(e0+el)-e0*el-R('1e-5')
    check('actual_cutoff_overlap_lower',overlap>R('.7848')>R('.78'),overlap)
    check('normalized_overlap_ratio',R('.78')/R('4.007')**2>R('.048'))
    check('functional_ratio_bound',R('4.007')*R('4.13')/R('.78')<R('21.3'))
    check('projected_range_inverse_factor',1+R('1.006')*R('21.3')<R('22.5'))
    check('symmetric_budget_not_required',el>R('.1') and overlap>0)
    f0=S.Matrix([3,0,0]);fl=S.Matrix([R('.6'),R('2.9'),0])
    c0=f0+S.Matrix([0,R('.006'),0]);cl=fl+S.Matrix([0,0,R('.12')]);theta=S.diag(1,R(1,2),0)
    check('projection_fixture_error_bounds',(c0-f0).dot(c0-f0)<e0**2 and (cl-fl).dot(cl-fl)<el**2)
    check('projection_fixture_reference_bounds',f0.dot(fl)>R('1.5') and f0.dot(f0)<16 and fl.dot(fl)<16)
    check('projection_fixture_cutoff_energy',(c0.T*(S.eye(3)-theta)*c0)[0]<R('.01')*c0.dot(c0))
    check('projection_fixture_positive_overlap',(cl.T*theta*c0)[0]>R('.78'))
    p=theta*c0/(c0.T*theta*c0)[0];proj=S.eye(3)-p*c0.T
    B=S.Matrix([[-cl[1]/cl[0],-cl[2]/cl[0]],[1,0],[0,1]])
    check('projection_fixture_full_range',cl.T*B==S.zeros(1,2) and B.rank()==2)
    check('projection_fixture_fixed_isomorphism',(proj*B).rank()==2)
    y=S.Matrix([-c0[1]/c0[0],1,2]);effective=y-p*(cl.dot(y)/cl.dot(p))
    solution=effective[1:3,0]
    check('projection_fixture_actual_inverse',B*solution==effective and proj*B*solution==y)
    bad=S.Matrix([cl[1],-cl[0],0]);badproj=S.eye(3)-bad*c0.T/c0.dot(bad)
    check('zero_overlap_projection_failure',cl.dot(bad)==0 and (badproj*B).rank()==1)


def independent_warped_fixtures():
    import mpmath as mp
    mp.mp.dps=85
    a=lambda q:q*(1+q*q/7+q**4/101)
    radius=lambda q:2+q*q/5+q**4/103
    u=lambda q:mp.mpf('.02')*q*q/(1+q*q)
    V=lambda q:mp.mpf('.03')/(1+q*q)
    tau=lambda q:1-1/radius(q)
    nu=lambda q:tau(q)*V(q)
    A=lambda q:mp.diff(a,q)/a(q)
    Rr=lambda q:mp.diff(radius,q)/radius(q)
    L=lambda f,q:mp.diff(f,q,2)+(A(q)+2*Rr(q))*mp.diff(f,q)
    D2=lambda f,q:mp.diff(f,q,2)+A(q)*mp.diff(f,q)
    base=lambda q:-2*mp.diff(a,q,2)/a(q)-4*(mp.diff(Rr,q)+A(q)*Rr(q))-6*Rr(q)**2+2/radius(q)**2
    linear=lambda q:-2*(u(q)+nu(q))*base(q)-2*D2(u,q)-6*L(nu,q)+4*u(q)/radius(q)**2
    for i,ell in enumerate([mp.mpf('.2'),mp.mpf('.7'),mp.mpf(1)]):
        k=lambda q:mp.exp(ell*V(q));T=lambda q:1+tau(q)*(k(q)-1)
        ss=lambda q:mp.log(T(q));w=lambda q:ell*u(q)+ss(q)
        ae=lambda q:mp.exp(w(q))*a(q);re=lambda q:T(q)*radius(q)
        def actual(q):
            lapse=mp.exp(w(q));wp=mp.diff(w,q)
            ar=mp.diff(ae,q)/(lapse*ae(q));rr=mp.diff(re,q)/(lapse*re(q))
            kk=(mp.diff(ae,q,2)-wp*mp.diff(ae,q))/(lapse**2*ae(q))
            e2=(mp.diff(re,q,2)-wp*mp.diff(re,q))/(lapse**2*re(q))-rr**2+ar*rr
            return -2*kk-4*e2-6*rr**2+2/re(q)**2
        def remainder(q):
            exp=mp.exp(-2*w(q));bb=k(q)-1-ell*V(q)*k(q)
            return base(q)*(exp*(1+2*ell*(u(q)+nu(q)))-1)+ell*(exp-1)*linear(q)-6*exp*(ell*(k(q)/T(q)-1)*L(nu,q)+bb/T(q)*L(tau,q)+ell**2*k(q)/T(q)*tau(q)*mp.diff(V,q)**2)+2*exp*(mp.exp(2*ell*u(q))-1-2*ell*u(q))/radius(q)**2
        missing=lambda q:6*mp.exp(-2*w(q))*ell**2*k(q)/T(q)*tau(q)*mp.diff(V,q)**2
        for j,q in enumerate([mp.mpf(1)/5,mp.mpf(2)/3,mp.mpf(7)/4]):
            observed=mp.diff(actual,q)-mp.diff(base,q)-ell*mp.diff(linear,q)
            expected=mp.diff(remainder,q);error=abs(observed-expected)
            check('independent_warped_derivative_'+str(i)+'_'+str(j),error<mp.mpf('1e-65'),mp.nstr(error,10))
            check('gradient_square_omission_'+str(i)+'_'+str(j),abs(mp.diff(missing,q))>mp.mpf('1e-10'))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 finite trace gradient and fixed projection','input_count':len(ps)}),flush=True)
    try:
        predecessors();coordinate_repair();low_and_third_budgets();nonlinear_derivative_budgets();fifth_jet_counts();complete_trace_and_projection();independent_warped_fixtures()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 finite trace gradient and fixed projection','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'radius_coordinate_conversion_repaired':good,'downstream_numeric_inputs_preserved':good,
         'finite_trace_gradient_bounded':good,'finite_full_adjoint_error_bounded':good,
         'fixed_projection_transversality':good,'projected_inverse_isomorphism':good,
         'trace_gradient_3e_minus4_proved':False,'useful_full_weighted_inverse_constant':False,
         'signed_physical_compatibility_evaluated':False,'nonlinear_family_preserved':False,
         'projected_nonlinear_solve':False,'S167_complete':False,'goal_complete':False}
    if not good:
        run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text()
        run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text()
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
