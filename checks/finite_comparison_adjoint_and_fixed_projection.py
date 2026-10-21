"""S167: finite adjoint phase and density-correct fixed projection overlap.

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
STEM = "finite-comparison-adjoint-and-fixed-projection"
OWN = ["checks/finite_comparison_adjoint_and_fixed_projection.py", "notes/"+STEM+".md"]
PREDECESSORS = [("finite-comparison-tt-realization-and-range", 319),
                ("stationary-scalar-phase-moment", 80),
                ("endpoint-preserving-compatible-projection", 262),
                ("finite-comparison-physical-tt-contacts", 304)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "finite actual adjoint phase and exact fixed-projection overlap criterion",
    "phase_coordinate": "each metric's own bolt-anchored proper radius",
    "phase_enclosure_radians_mod_pi": ".025",
    "fixed_source_conversion": "Jvol^-1 p_cut before actual York transpose",
    "projection_fixed_at_seed": True,
    "overlap_numerically_evaluated": False,
    "fixed_projection_transversality": False,
    "series_order_9_tail_coefficient": "less than 1e-16; other evaluation errors excluded",
    "fixture_scope": "finite matrix real/rank-one resolvent and projected range; rational phase and geometric budgets",
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
    for p in ["finite-comparison-tt-realization-and-range", "stationary-scalar-phase-moment",
              "endpoint-preserving-compatible-projection", "finite-comparison-physical-tt-contacts",
              "finite-comparison-scalar-resolvent", "finite-comparison-weighted-mixed-form",
              "stationary-global-weighted-map", "general-background-newton-operator"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)










def phase_and_geometric_budgets():
    beta=R('1.000002');q=R('.025');l0=R('.02405')
    check('real_auxiliary_contraction',beta*l0<q)
    check('real_auxiliary_inverse',beta/(1-q)<R('1.026'))
    a=beta*l0/(1-q)
    check('mixed_phase_coefficient',a<R('.0247'))
    eta=R(72000,4*10**17)
    check('geometric_scattering',eta<R('2e-13'))
    check('geometric_phase_coefficient',eta/(1-eta)<R('2.1e-13'))
    check('proper_radius_phase_enclosure',a+eta/(1-eta)<R('.025'))
    tail=beta*q**10/(1-q)
    check('ten_term_series_tail',tail<R('1e-16'),tail)
    check('tail_ratio_each_term',q==R(1,40))
    x=R('3.22')
    check('log25_upper_bound',sum(x**j/S.factorial(j) for j in range(21))>25)
    d=R('1.001')*(177*R('3e-8')+R('1.5')*R('3e-8')*x+R('5e-9'))
    check('proper_length_shift_upper',d<R('5.5e-6'),d)
    t=S.symbols('t',positive=True)
    zero('outer_shift_tail',S.integrate(15/t**4,(t,1000,S.oo))-R('5e-9'))
    gmin=R('.999')/2000**10
    check('frequency_upper',60*R(22,7)/gmin+R(17,12)<(5*10**17)**2)
    zero('phase_advance_allowance',R('5.5e-6')*(5*10**17)-R('2.75e12'))
    check('advance_allowance_not_small',R('2.75e12')>R('.025'))
    dens=R('1.2e-7')
    check('density_ratio_margin',dens/(1-dens)<R('1.3e-7'))


def real_resolvent_fixture():
    # Enlarged finite fixture; no physical kernel evaluation is implied.
    A=S.diag(2,R('1e-14'));J=S.diag(1,-1)
    LR=S.Matrix([[R(1,1000),R(1,3000)],[R(1,3000),R(1,2000)]])
    g=S.Matrix([R(1,1000),R(1,500)]);h=S.Matrix([R(2,7),-R(3,11)])
    c0=R(10);m0=R(5,13);I=S.eye(2)
    N=(I+A*J).inv();K=N*LR*J;MR=I+A*J-LR*J
    zero('real_factorization',(MR-(I+A*J)*(I-K)).norm()**2)
    zero('J_real_symmetry',(MR.T*J-J*MR).norm()**2)
    zero('J_inverse_symmetry',((J*MR.inv()).T-J*MR.inv()).norm()**2)
    gR=MR.inv()*g;a=(g.T*J*gR)[0]/c0
    L=LR+S.I*g*g.T/c0;M=I+A*J-L*J
    zero('rank_one_factorization',(M-(MR-S.I*g*g.T*J/c0)).norm()**2)
    zero('rank_one_inverse',(M.inv()*g-gR/(1-S.I*a)).norm()**2)
    hp=h+S.I*g*m0/c0;y=M.inv()*hp
    moment=m0+(gR.T*J*h)[0]
    zero('moment_outgoing_identity',m0+(g.T*J*y)[0]-moment/(1-S.I*a))
    scat=(1+S.I*a)/(1-S.I*a)
    zero('unitary_scattering',S.expand_complex(scat*S.conjugate(scat))-1)
    check('mixed_phase_is_real',S.im(a)==0)
    wrong=m0+(g.T*J*(I-LR*J).inv()*h)[0]
    check('elliptic_subtraction_matters',wrong!=moment)
    check('wrong_rank_one_normalization_detected',M.inv()*g!=gR/(1-S.I*a/2))
    beta=R('1.000002');q=R('.025')
    check('fixture_N_bound',max(abs(N[i,i]) for i in range(2))<beta)
    check('fixture_K_bound',sum(v*v for v in K)<q*q)
    approx=sum((K**j*N*g for j in range(10)),S.zeros(2,1))
    diff=gR-approx;tail=beta*q**10/(1-q)
    check('direct_inverse_series_vector',sum(v*v for v in diff)<tail**2*(g.T*g)[0])
    err=(diff.T*J*h)[0]
    check('direct_inverse_series_moment',err**2<tail**2*(g.T*g)[0]*(h.T*h)[0])
    check('actual_fixture_tail_nonzero',diff!=S.zeros(2,1))


def fixed_projection_algebra():
    a,b=S.symbols('a b',real=True)
    B=S.Matrix([[1,0],[0,1],[a,b]])
    chi=S.Matrix([-a,-b,1]);eta=S.Matrix([[0,0,1]])
    pc=S.Matrix([R(1,2),0,1]);Pi=S.eye(3)-pc*eta
    zero('projection_idempotent',(Pi*Pi-Pi).norm()**2)
    zero('seed_projection_normalization',(eta*pc)[0]-1)
    zero('physical_cokernel',(chi.T*B).norm()**2)
    sigma=(chi.T*pc)[0]
    PB=(Pi*B)[:2,:]
    zero('projected_determinant_is_overlap',PB.det()-sigma)
    y=S.Matrix([S.Symbol('y1',real=True),S.Symbol('y2',real=True),0])
    r=y-pc*(chi.T*y)[0]/sigma;u=r[:2,0]
    zero('range_correction',(chi.T*r)[0])
    zero('full_inverse_formula',(B*u-r).applyfunc(S.cancel).norm()**2)
    zero('projected_inverse_formula',(Pi*B*u-y).applyfunc(S.cancel).norm()**2)
    B1=B.subs({a:1,b:2});chi1=chi.subs({a:1,b:2})
    check('transverse_fixture',B1.rank()==2 and (Pi*B1).rank()==2 and (chi1.T*pc)[0]==R(1,2))
    Bz=B.subs({a:2,b:2});chiz=chi.subs({a:2,b:2})
    check('zero_overlap_fixture',Bz.rank()==2 and (Pi*Bz).rank()==1 and (chiz.T*pc)[0]==0)
    u0=S.Matrix([R(1,2),0])
    zero('projected_kernel_witness',(Bz*u0-pc).norm()**2)
    zero('projected_kernel_zero',(Pi*Bz*u0).norm()**2)
    # Physical endpoint volume is D; fixed output is multiplied by D.
    D=S.diag(2,3,5);A=D.inv()*Bz
    zero('density_fixed_derivative',(D*A-Bz).norm()**2)
    zero('moving_adjoint_equation',(A.T*D*chiz).norm()**2)
    zero('correct_fixed_overlap',(chiz.T*D*(D.inv()*pc))[0]-(chiz.T*pc)[0])
    check('omitted_density_gives_false_transversality',(chiz.T*D*pc)[0]==3)
    Rt=S.Matrix([[1,1,0],[0,1,1]])
    check('density_does_not_commute_with_transpose',Rt*D.inv()*pc!=S.diag(R(1,2),R(1,3))*Rt*pc)
    scale=S.symbols('scale',positive=True)
    zero('normalization_cannot_improve_ratio',scale**2*(chi1.T*chi1)[0]/(scale*(chi1.T*pc)[0])**2-(chi1.T*chi1)[0]/(chi1.T*pc)[0]**2)


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 finite adjoint and fixed projection','input_count':len(ps)}),flush=True)
    try:
        predecessors();phase_and_geometric_budgets();real_resolvent_fixture();fixed_projection_algebra()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 finite adjoint and fixed projection','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'actual_finite_adjoint_constructed':good,'proper_radius_phase_enclosed':good,
         'density_correct_overlap_formula':good,'projected_inverse_criterion':good,
         'series_tail_bounded':good,'overlap_evaluated':False,
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
