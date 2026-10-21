"""S167: weighted mixed form and complete trace at finite comparisons.

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
STEM = "finite-comparison-weighted-mixed-form"
OWN = ["checks/finite_comparison_weighted_mixed_form.py", "notes/"+STEM+".md"]
PREDECESSORS = [("finite-comparison-scalar-curvature-remainder", 314),
                ("finite-comparison-scalar-resolvent", 297),
                ("finite-comparison-quantum-disk-tt", 378)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "actual curvature derivative, complete trace and weighted mixed energy-dual estimates",
    "trace_bound": "abs(Etr)<3.2e-12 lambda/t",
    "weight": "omega=t^(-1/2)",
    "mixed_f_coefficient": "4e6 sqrt(lambda)",
    "mixed_g_coefficient": "1000 sqrt(lambda/G)",
    "input_normalization": "sqrt(.90 E_ell[c])",
    "outgoing_contraction_conditional_on_TT_realization": ".025",
    "full_candidate_operator_realized": False,
    "fixture_scope": "exact finite metric jet and rational budgets, compact weighted divergence pairings",
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
    for p in ["finite-comparison-scalar-curvature-remainder", "finite-comparison-scalar-resolvent",
              "finite-comparison-quantum-disk-tt", "stationary-completed-metric-response-norm",
              "finite-comparison-classical-disk-contact", "stationary-scalar-outgoing-inverse",
              "stationary-weighted-tt-inverse", "stationary-global-weighted-map"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)






def coefficient_lemma():
    a,b,c,d=S.symbols('M1 M2 M3 M4',nonnegative=True)
    N1=a;N2=b+8*a*a;N3=c+24*a*b+96*a**3
    G0=R(3,2)*a
    G1=R(3,2)*b+6*a*a
    G2=R(3,2)*c+18*a*b+48*a**3
    G3=R(3,2)*d+24*a*c+18*b*b+288*a*a*b+576*a**4
    zero('connection_first_product_rule',G1-(6*N1*a+R(3,2)*b))
    zero('connection_second_product_rule',G2-(6*N2*a+12*N1*b+R(3,2)*c))
    zero('connection_third_product_rule',G3-(6*N3*a+18*N2*b+18*N1*c+R(3,2)*d))
    Q0=2*G1+8*G0**2; Q1=2*G2+16*G0*G1; Q2=2*G3+16*G0*G2+16*G1**2
    grad=32*(Q1+(4*a+16*G0)*Q0)
    Hess=4*(16*Q2+128*a*Q1+64*N2*Q0+4*G0*(16*Q1+64*a*Q0))
    zero('full_covariant_curvature_gradient',grad-(96*c+4992*a*b+34560*a**3))
    zero('full_scalar_Hessian_count',Hess-(192*d+8064*a*c+5376*b*b+165888*a*a*b+506880*a**4))
    p3=c+a*b+a**3;p4=d+a*c+b*b+a*a*b+a**4
    check('J_C_gradient_coefficient_majorant',all(k>=0 for k in S.Poly(10**6*p3-2*grad,a,b,c,d).coeffs()))
    check('scalar_Hessian_coefficient_majorant',all(k>=0 for k in S.Poly(10**6*p4-Hess,a,b,c,d).coeffs()))
    check('outer_product_telescope',16*2**15<10**6)
    vals={a:30,b:3000,c:10**6,d:4*10**8}
    zero('P3_value',p3.subs(vals)-1117000)
    zero('P4_value',p4.subs(vals)-442510000)
    check('core_gradient_physical_weight',(40*(R('2.5')**3*10**6*p3.subs(vals))**2)<R('1e15')**2)
    check('core_Hessian_physical_weight',40*R('2.5')**4*10**6*p4.subs(vals)<R('1e20'))
    # Use sqrt(40)<7 and t>=40 for the outer weight comparison.
    check('outer_gradient_physical_weight',R(1000)*R('2.5')**3*10**6*p3.subs(vals)*7/40**2+R(1000)*7/40**3<R('1e15'))
    check('outer_Hessian_physical_weight',R(1000)*R('2.5')**4*10**6*p4.subs(vals)/40+R(4)/40**4<R('1e20'))
    check('connection_omission_control',S.expand(Hess-4*(16*Q2+128*a*Q1+64*N2*Q0))!=0)


def geometric_jets():
    # Coefficient sums for f=F/C; thermal data have 0<delta<=1.
    xx=R('7.9e9');cx=R('5.9');dd=R('1e-12')
    co={2:1/cx,0:1/(cx*xx),-1:(4+3/xx)/cx,-2:(3+(1+dd)/xx)/cx,-3:2*dd/(cx*xx),-4:dd/(cx*xx)}
    D=lambda j:sum(abs(p)**j*v for p,v in co.items())
    check('seed_Da',D(2)/2<2)
    check('seed_Db',D(2)/2+D(1)<5)
    check('seed_DDa',D(3)/2<4)
    check('seed_DDb',D(3)/2+D(2)<8)
    check('seed_phi_third',2+5+2*3<=13)
    check('seed_phi_fourth',2*13+(4+8+2*5)+3*13<100)
    # Rhat=-D^2 f+3 Df-2f; constant and t^2 derivative terms vanish.
    RD=lambda j:sum(abs(-p*p+3*p-2)*abs(p)**j*v for p,v in co.items() if p!=0)
    check('seed_curvature_D',RD(1)<40)
    check('seed_curvature_DD',RD(2)<80)
    check('seed_Hess_curvature',R('.85')*80+(R('1.87')+R('2.72'))*40<252)
    check('seed_normal_phi3',13+10*2<40)
    check('seed_normal_phi4',100+100*(2*5+40)<6000)
    phi=[1,5,40,6000];L=[10,1000,300000,100000000]
    for j in range(1,5):
        perturb=R('.01')*S.factorial(j)*100**j
        check('log_multiplier_jet_'+str(j),2*(phi[j-1]+2*perturb)<L[j-1])
    l1,l2,l3,l4=L
    bells=[1,l1,l2+l1*l1,l3+3*l1*l2+l1**3,l4+4*l1*l3+3*l2*l2+6*l1*l1*l2+l1**4]
    zero('Bell4_value',S.Integer(bells[4])-115610000)
    base=[1,0,10,10000,10000000];Ms=[1,30,3000,10**6,4*10**8]
    for j in range(1,5):
        total=R('1.001')*sum(S.binomial(j,k)*bells[k]*base[j-k] for k in range(j+1))
        check('actual_metric_component_jet_'+str(j),total<Ms[j],total)
    # Outer cubic/remainder and finite normal products, before physical derivatives.
    check('outer_point_log_jets',6*R('1e-9')*27+R('1e-7')*81*R('.12')<R('1.2e-6'))
    for j in range(5):
        check('outer_metric_difference_jet_'+str(j),R('.0012')*S.factorial(j)*100**j/40<R('.001')*Ms[j])
    check('core_cover',1/(R('.412')*40)>R('.05'))
    check('physical_derivative_conversion',R('6.1')*R('1.001')**2<R('2.5')**2)


def value_and_trace():
    # Small exponentials leave the derivative and measure terms in each pullback.
    exp=R('1.001'); changed=R('2.01')*10
    K=exp*(10+2*10)+R('10.5')*changed
    E2=exp*(10+2*10)+22*changed
    AR=exp*(2*10+R('2.24')*10+100/R(40)**3)+11*changed
    Rr2=exp*(2*R('2.24')*10+100/R(40)**3)+5*changed
    check('outer_K_value',K<250);check('outer_E2_value',E2<500)
    check('outer_AR_value',AR<300);check('outer_Rr2_value',Rr2<150)
    check('outer_radial_sphere',500+300+150<=950)
    check('outer_sphere_sphere',Rr2+R('2.01')*10/R('7.9e9')/40**2<151)
    check('outer_Rm_value',2*(250+2*950+2*300+151)<6000)
    check('outer_J_value',2*(2*3*950)<12000)
    check('inner_weighted_J',7+R('.001')*40<8)
    check('outer_weighted_J',R(12000)/40**2+R(7)/40**3<8)
    check('inner_weighted_C',40+R('.001')*40<50)
    check('outer_weighted_C',R(6040)/40**2<50)
    check('Schouten_orthogonal',R(8)**2/4+R(13)**2/144<6**2)
    gg=R(1,2000**10);z=R('3.1e-12')
    trace=(1+11*gg/540)*z+gg/180*(R('2e20')+R(11,72)*z*z+R(7,12)*50**2+R(11,6)*8**2)+R('10.01')*R('4e-27')/180
    check('complete_physical_trace',trace<R('3.2e-12'),trace)
    check('anomaly_derivative_retained',gg*R('2e20')/180>0)


def mixed_and_outgoing():
    check('mixed_vector_derivative_coefficient',R(17,6)*8+R(7,2)*50<200)
    check('mixed_source_derivative_coefficient',R(7,2)+R(11,6)<6)
    weyl=R(7,12)*(2*50*55+8*50*6+4*50**2)
    ric=R(11,6)*(8*(2*55+2*15+R(13,2))+2*8**2)
    classical=R(11,36)*13*8
    check('curvature_algebra_reserve',weyl+ric+classical<R('1e6'))
    check('Hessian_plus_algebra',R('1e20')+R('1e6')<R('2e20'))
    cov=75**2+8*R('9.6')**2+4*R('2.25')**2*R('3.002')**2
    check('full_covariant_TT_gradient',cov<82**2,cov)
    check('h_units',R('4.8')**2*8*R(22,7)<25**2)
    check('Y_units',R('4.83')**2*8*R(22,7)<25**2)
    check('grad_h_units',82**2*8*R(22,7)<420**2)
    check('mixed_normalization',R(5)/(4*R('.90'))<R('1.2')**2)
    af=R('1.2')*((R('4e-17')/180)*(200*420+R('1e15')*25)+25*R('3.2e-12')/R('3e-17'))
    gg=R(1,2000**10)
    ag=R('1.2')*(8*25+(gg/180)*(R('6e15')*420+R('2e20')*25)+R('20.02')*25/180*R('4e-27'))
    check('mixed_f_energy_bound',af<R('3.3e6'),af)
    check('mixed_g_energy_bound',ag<241,ag)
    check('local_scalar_vector',R('3.1e-12')/36+R(11,9)*8<4**2)
    check('local_scalar_value',R('3.1e-12')*20*R(22,7)<R('2e-5')**2)
    check('combined_f_bound',R('3.3e6')**2+4**2<R('4e6')**2)
    check('combined_g_bound',R(241)**2+R('2e-5')**2<1000**2)
    eta=R('2e-13'); k=R('4e17')
    check('physical_single_derivative',(2+5/k)/(1-eta)<3)
    check('physical_double_derivative',2+(4*72000+20+25/k)/(k*(1-eta))<3)
    check('weight_integral',178+R('1.5')/40<200)
    check('weighted_scalar_kernel',R('1.001')*200<201)
    check('weighted_double_kernel',1+600*k<601*k)
    outgoing=R(201)*1000**2/(180*k)+R(1200,180)*1000*R('4e6')*R('4e-17')+601*R('4e6')**2/k
    check('conditional_outgoing_norm',outgoing<R('.02405'),outgoing)
    check('conditional_outgoing_contraction',R('1.000002')*R('.02405')<R('.025'))
    check('identity_distribution_not_discarded',1+600*k>600*k)


def compact_pairing_fixtures():
    s=S.symbols('s',real=True)
    m=s*(1+s*s)*(2+s*s)**2
    for power in [3,4,5]:
        tau=(1-s*s)**power
        Y=s*(1-s*s)**4; E=1+s*s
        G=s*(1-s*s)**3; f=(1-s*s)**4
        div=lambda v:S.diff(m*v,s)/m
        Ttrace=-div(G)+f; B=-Ttrace+div(E*Y)
        difference=S.cancel(tau*B*m+S.diff(tau,s)*(G+E*Y)*m+tau*f*m)
        zero('physical_divergence_pairing_'+str(power),S.integrate(difference,(s,0,1)))
        contact=S.integrate(S.cancel(tau*div(E*Y)*m),(s,0,1))
        check('off_shell_contact_omission_'+str(power),contact!=0,contact)
        zero('bolt_and_outer_flux_'+str(power),(m*tau*(G+E*Y)).subs(s,0)+(m*tau*(G+E*Y)).subs(s,1))
    # The antisymmetric local part is repaired only by the original density term.
    F=S.Matrix([2,3]);l=S.Matrix([4,2]);H=S.Matrix([[5,7],[7,11]])
    DF=H-F*l.T
    check('density_term_retained_for_symmetry',DF!=DF.T and DF+F*l.T==H)


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 weighted mixed form','input_count':len(ps)}),flush=True)
    try:
        predecessors();coefficient_lemma();geometric_jets();value_and_trace();mixed_and_outgoing();compact_pairing_fixtures()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 weighted mixed form','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'weighted_curvature_derivatives_bounded':good,'complete_trace_defect_bounded':good,
         'mixed_coefficient_estimate_complete':good,'conditional_outgoing_contraction_bounded':good,
         'full_candidate_operator_realized':False,'full_candidate_Hessian_inverse':False,
         'physical_mass_fast_domain_membership_proved':False,'signed_compatibility_evaluated':False,
         'nonlinear_family_preserved':False,'S167_complete':False,'goal_complete':False}
    if not good:
        run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text()
        run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text()
    dest=ROOT/'evidence'/(STEM+'.json')
    data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
