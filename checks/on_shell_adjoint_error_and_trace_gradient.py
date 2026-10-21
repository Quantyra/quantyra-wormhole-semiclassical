"""S167: on-shell full-adjoint error and local trace-gradient reduction.

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
STEM = "on-shell-adjoint-error-and-trace-gradient"
OWN = ["checks/on_shell_adjoint_error_and_trace_gradient.py", "notes/"+STEM+".md"]
PREDECESSORS = [("cold-reference-overlap-and-adjoint-budget", 316),
                ("finite-comparison-weighted-mixed-form", 335),
                ("finite-comparison-adjoint-and-fixed-projection", 299),
                ("stationary-scalar-outgoing-inverse", 107)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "seed full-adjoint bound; finite full error reduced to an unproved local threshold",
    "dual_norm": "sqrt(lambda) integral t^-4 |field|^2 dmu0",
    "trace_gradient": "d_tr = L2(dz) norm of partial_z(Etr/lambda), z=sqrt(lambda)rho_ell",
    "TT_vertex_bound": "xi < 30.1 d_tr + 1e-8",
    "full_error_bound": "e_ell < 6.7 xi + .0064",
    "seed_full_error_upper": ".007",
    "finite_sufficient_threshold": "d_tr <= 3e-4; NOT proved",
    "finite_trace_gradient_threshold_proved": False,
    "fixture_scope": "exact trace and integration identities, rational operator majorants, finite block systems",
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
    for p in ["cold-reference-overlap-and-adjoint-budget", "finite-comparison-weighted-mixed-form",
              "finite-comparison-adjoint-and-fixed-projection", "stationary-scalar-outgoing-inverse",
              "stationary-weighted-tt-inverse", "finite-comparison-disk-return-form",
              "finite-comparison-classical-disk-contact", "finite-comparison-tt-realization-and-range",
              "finite-comparison-scalar-curvature-remainder"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)














def on_shell_identities():
    eps,lam,Z,r,T_r,hess,dC,dJ,dM=S.symbols('eps lam Z r T_r hess dC dJ dM')
    a=1+11*eps*lam/3;R4=Z-12*lam
    actual=(1-11*eps*R4/36)*r-eps*T_r-eps*hess+eps*(R(7,12)*dC+R(11,6)*dJ-20*dM)
    deltaH=R(7,12)*dC+R(11,6)*dJ-20*dM-R(11,36)*Z*r
    zero('exact_on_shell_decomposition',actual-((a*r-eps*T_r)+eps*(-hess+deltaH)))
    zero('reference_eigenvalue_cancels_entire_trace_vertex',a-eps*(a/eps))
    # Fourier integration of smooth endpoint-vanishing functions; not a physical TT geometry.
    x=S.symbols('x',real=True);phi=S.cos(3*x);r0=S.sin(x)**4*(1+S.cos(3*x))
    integrate=lambda f:S.integrate(S.expand_trig(f).expand(),(x,0,2*S.pi))
    classical=integrate(phi*R(9,5)*r0)
    box=integrate(phi*R(1,5)*S.diff(r0,x,2))
    zero('independent_on_shell_integration',classical+box)
    zero('nonzero_classical_term_retained',classical-R(27,40)*S.pi)
    check('dropping_laplacian_detected',classical!=0)
    trace=2+S.cos(2*x);Y=S.sin(x)
    lhs=integrate(trace*Y*S.diff(phi,x))
    value=integrate(phi*trace*S.diff(Y,x))
    gradient=integrate(phi*S.diff(trace,x)*Y)
    zero('independent_return_integration',lhs+value+gradient)
    zero('nonzero_trace_gradient_term',gradient-S.pi)
    check('dropping_trace_gradient_detected',lhs+value!=0)
    check('fixture_endpoint_terms_zero',all(S.simplify(f.subs(x,p))==0 for f in [r0,S.diff(r0,x),Y] for p in [0,2*S.pi]))


def source_budgets():
    sg=R('4e-17');gamma=sg**2;root200=R('14.143')
    check('sqrt200_majorant',root200**2>200)
    check('pairing_factor',R('1.2')**2>R(5,4)/R('.90'))
    check('physical_h_constant',(R('4.8')**2*8*R(22,7))<25**2)
    check('physical_Y_constant',(R('4.83')**2*8*R(22,7))<25**2)
    check('physical_covariant_gradient_constant',82**2*8*R(22,7)<420**2)
    div=R('3.002')+2*R('.713')*R('4.8')+2*R('2.25')*R('4.83')
    check('divergence_return_constant',div**2*8*R(22,7)<160**2,div)
    check('reference_derivative_bessel_reserve',2+2*72000/(4*10**17)/(1-R('2e-13'))+5/(4*10**17)*R('1.000001')<3)
    check('trace_gradient_pairing',R('1.2')*25*R('1.000001')<R('30.1'))
    value=R('1.2')*160*R('3.2e-12')*R('1.000001')*root200
    grad=R('1.2')/180*200*420*3*root200*(14*sg)
    curvature=R('1.2')/180*R('6e15')*420*gamma*R('1.000001')*root200
    hess=R('1.2')/180*R('2e20')*25*gamma*R('1.000001')*root200
    magnetic=R('1.2')/180*R('20.02')*25*R('4e-27')*R('1.000001')*root200
    check('trace_value_return_budget',value<R('9e-9'),value)
    check('on_shell_derivative_budget',grad<R('2e-11'),grad)
    check('curvature_gradient_budget',curvature<R('4e-13'),curvature)
    check('hessian_algebra_budget',hess<R('8e-13'),hess)
    check('magnetic_vertex_budget',magnetic<R('2e-25'),magnetic)
    check('non_gradient_TT_vertex_budget',value+grad+curvature+hess+magnetic<R('1e-8'),value+grad+curvature+hess+magnetic)
    scalar=(4*3+R('2e-5')*R('1.000001')/13)*root200
    check('separate_local_reference_vertex',scalar<180,scalar)
    check('local_reference_scaled_budget',180*14==2520)


def block_majorants():
    b=R('2.6e-13');sg=R('4e-17');kb=4*10**17
    check('N_TT_diagonal',1+b/(4*(1-b))<R('1.000001'))
    check('N_off_diagonal',b/(4*(1-b)**2)<R('2.6e-7')**2)
    check('N_local_diagonal',1/(1-b)<R('1.000001'))
    f=[R('4e6'),S.Integer(4)];g=[S.Integer(1000),R('2e-5')]
    N=S.Matrix([[R('1.000001'),R('2.6e-7')],[R('2.6e-7'),R('1.000001')]])
    L=S.Matrix([[R('.02405'),R('2.405e-8')],[R('2.405e-8'),R('2.405e-14')]])
    for i in range(2):
        for j in range(2):
            val=R(201,180)*g[i]*g[j]/kb+R(600,180)*(f[i]*g[j]+g[i]*f[j])*sg+601*f[i]*f[j]/kb
            check('L_block_'+str(i)+str(j),val<L[i,j],val)
    K=N*L
    check('positive_majorant_series_converges',max(sum(K[i,j] for j in range(2)) for i in range(2))<R('.025'))
    inverse=(S.eye(2)-K).inv()*N
    bound=S.Matrix([[R('1.026'),R('3e-7')],[R('4e-7'),R('1.000002')]])
    for i in range(2):
        for j in range(2):check('real_inverse_block_'+str(i)+str(j),0<=inverse[i,j]<bound[i,j],inverse[i,j])
    # Independent positive finite block with a deliberately large TT block.
    A=S.Matrix([[10**6,R('1e-5')],[R('1e-5'),R('1e-14')]])
    J=S.diag(1,-1);LR=S.Matrix([[R('.01'),R('1e-9')],[R('1e-9'),R('1e-15')]])
    check('fixture_A_positive',A[0,0]>0 and A.det()>0 and A[1,1]<b)
    check('fixture_large_TT_block',A[0,0]>10**5)
    MR=S.eye(2)+A*J-LR*J;direct=MR.inv();n=(S.eye(2)+A*J).inv()
    check('fixture_reconciles_positive_elimination',direct==(S.eye(2)-n*LR*J).inv()*n)
    for i in range(2):
        for j in range(2):check('fixture_inverse_bound_'+str(i)+str(j),abs(direct[i,j])<bound[i,j])
    source=S.Matrix([1,10**18]);response=direct*source
    check('fixture_full_equation',MR*response==source)
    check('fixture_direct_TT_budget',abs(response[0])<bound[0,0]*source[0]+bound[0,1]*source[1])
    check('fixture_direct_local_budget',abs(response[1])<bound[1,0]*source[0]+bound[1,1]*source[1])
    bad=(S.eye(2)-LR*J).inv()*source
    check('dropping_positive_TT_elimination_detected',MR*bad!=source and abs(bad[0]-response[0])>10**6)


def reconstruction_and_seed():
    sg=R('4e-17');gamma=sg**2;u=R('1.000001')
    check('finite_disk_return_value',R('4.6')*R('4.8')<25)
    check('TT_return_constant',R('25.1')**2/(9*3)<25)
    check('scalar_weight_product',32*200==80**2)
    for i,(f,g,bound) in enumerate([(R('4e6'),S.Integer(1000),R('1.5')),(S.Integer(4),R('2e-5'),R('3e-8'))]):
        elliptic=u*(g/4+f*sg/4)/180
        outgoing=u*(80*R('1.001')*g/13+240*f)*sg/180
        check('scalar_reconstruction_channel_'+str(i),elliptic+outgoing<bound,elliptic+outgoing)
    xicoeff=R('6.5')*R('1.026')+R('3e-8')*R('4e-7')
    constant=2520*(R('6.5')*R('3e-7')+R('3e-8')*R('1.000002'))
    check('full_error_xi_coefficient',xicoeff<R('6.7'),xicoeff)
    check('local_channel_error_constant',constant<R('.0051'),constant)
    check('amplitude_normalization_cost',2*R('.0247')**2<R('.0013'))
    check('full_error_constant',constant+R('.0013')<R('.0064'))
    conditional=R('6.7')*(R('30.1')*R('3e-4')+R('1e-8'))+R('.0064')
    check('conditional_finite_full_error',conditional<R('.07'),conditional)
    t=S.symbols('t',positive=True);H=4/t**5-6/t**6
    for j,expect in enumerate([56,372,2856],1):
        derivative=S.diff(H,t,j).expand()
        total=sum(abs(term.as_coeff_Mul()[0]) for term in S.Add.make_args(derivative))
        zero('seed_scalar_derivative_sum_'+str(j),total-expect)
    check('seed_first_proper_derivative',R('2.237')**2>5 and 56*R('2.237')*R('1.3e-22')<R('1.7e-20'))
    check('seed_laplacian_first_coefficient',2*5+22==32)
    check('seed_laplacian_coefficient_derivative',2*22+2*5+21==75)
    box=22*372+5*2856+75*56+32*372
    zero('seed_differentiated_laplacian_sum',S.Integer(box)-38568)
    check('seed_differentiated_laplacian_reserve',box*R('2.237')<90000)
    # Full seed trace derivative, including magnetic and quadratic anomaly terms.
    first=(1+R(11,3)*gamma/180)*R('1.7e-20')
    lap=gamma/180*90000*R('1.3e-22')
    quadratic=gamma/180*(R(7,12)*2*40*1000+R(11,6)*2*7*100+R(11,36)*(10*R('1.3e-22'))*R('1.7e-20'))
    mag=R(40,180)*R('2.24')*R('4e-27')
    check('seed_complete_trace_gradient',first+lap+quadratic+mag<R('2e-20'),first+lap+quadratic+mag)
    check('seed_gradient_L2_bound',R('2e-20')**2*200<R('3e-19')**2)
    seed=R('6.7')*(R('30.1')*R('3e-19')+R('1e-8'))+R('.0064')
    check('unconditional_full_seed_adjoint_error',seed<R('.007'),seed)
    # A small value does not bound its derivative: a concrete exact counterexample.
    x=S.symbols('x',real=True);small=R('1e-12')*S.sin(10**20*x)
    zero('trace_value_does_not_bound_gradient',S.diff(small,x).subs(x,0)-10**8)


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 on-shell adjoint error and trace gradient','input_count':len(ps)}),flush=True)
    try:
        predecessors();on_shell_identities();source_budgets();block_majorants();reconstruction_and_seed()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 on-shell adjoint error and trace gradient','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'on_shell_TT_vertex_bounded':good,'channelwise_adjoint_error_bound':good,
         'seed_full_adjoint_error_bounded':good,'finite_error_reduced_to_trace_gradient':good,
         'finite_trace_gradient_threshold_proved':False,'finite_full_adjoint_budget_proved':False,
         'fixed_projection_transversality':False,'useful_full_weighted_inverse_constant':False,
         'signed_physical_compatibility_evaluated':False,'nonlinear_family_preserved':False,
         'S167_complete':False,'goal_complete':False}
    if not good:
        run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text()
        run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text()
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
