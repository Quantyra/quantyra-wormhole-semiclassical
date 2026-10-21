"""S168: canonical frequency core bound.

Exact metric algebra and rational budgets; continuum claims remain ordinary proofs.
"""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
STEM = "canonical-frequency-core-bound"
OWN = ["checks/canonical_frequency_core_bound.py", "notes/"+STEM+".md", "evidence/canonical-frequency-core-source-review.json"]
PREDECESSORS = [("local-coupled-frequency-polynomial", 326)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope":"Actual centered-seed canonical core comparison and auxiliary logarithmic inequality",
    "matrix_seed":16805,"matrix_frequencies":[0,.5,2,20,100],
    "matrix_absolute_tolerance":1e-10,
    "torus_grids":[64,128],"torus_frequencies":[1,3,10],
    "torus_relative_tolerance":1e-10,
    "fixture_h":".03cos(Y)+.02sin(2Y)","fixture_t":".3cos(2Y)+.1sin(3Y)",
    "fixtures_are_exterior_spectral_data":False,
    "actual_wormhole_instability_proved":False,"actual_causal_inverse_proved":False,
    "actual_uniform_logarithmic_comparison_proved":False,
    "order_reduction_adopted":False,"formal_verification":False,
}
NUMERICS=[]


def check(name, ok, detail=""):
    assert name not in {row["name"] for row in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expr):
    expr = S.cancel(S.sympify(expr).doit())
    check(name, expr == 0, expr)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_bytes())
    return data.get("runs", [data])[-1]


def inputs(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def paths():
    found = set(OWN) | {"notes/stationary-coupled-asymptotic-obstruction.md", "notes/spherical-anomaly-action.md", "notes/dirac-spherical-source-reduction.md", "notes/stationary-lll-quantum-exterior.md", "notes/coupled-fermion-domain.md", "notes/exact-spherical-em-reduction.md", "notes/stationary-shape-logarithmic-response.md", "notes/stationary-curved-quantum-tt.md"}
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
    for p in ["general-background-newton-operator", "stationary-scalar-phase-moment",
              "on-shell-adjoint-error-and-trace-gradient", "on-shell-finite-range-response",
              "cold-projection-source-lower-bound"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)







def budgets_and_identities():
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('primary_source_scope',review['story']=='S168' and len(review['primary_sources'])==1 and not review['actual_wormhole_instability_proved'])
    eps=R(1,1000);K=S.Integer(10)**13;radius=R(1,400);bsq=eps/(4*K)
    check('split_inside_local_radius',bsq<radius**2)
    check('log_ratio_below_28',radius**2/bsq<R(8,3)**28)
    pi_hi=R(22,7);pi_lo=S.Integer(3)
    check('two_dimensional_fourth_loss',eps*pi_hi/4<R(1,1000))
    grad=K*(R(1,4)+28)*pi_hi+eps*pi_hi/6+S.Integer(10)**13+68/(36*pi_lo)+400/(12*pi_lo)
    check('pre_seam_gradient_reserve',grad<S.Integer(10)**17,grad)
    check('gluing_zero_reserve',6*S.Integer(10)**24<S.Integer(10)**25)
    check('p0_density_conversion',S.Integer(10)**20*400**2<2*S.Integer(10)**25)
    # exp(.02)<=1/(1-.02); the geometric-series majorant is exact.
    check('collar_weight_ratio',1/(1-R(2,100))<R(5,4))
    zero('collar_Hardy_coefficient',R(1,10)*R(5,4)-R(1,8))
    check('two_sided_gradient_budget',S.Integer(10)**17+2*S.Integer(10)**25/8<3*S.Integer(10)**24)
    check('two_sided_zero_budget',S.Integer(10)**25+22*S.Integer(10)**25<3*S.Integer(10)**26)
    check('log2000_less_than_8',R(8,3)**8>2000)
    check('kappa_upper',R(33,2)+R(1,10**21)<17)
    check('Garding_fourth_reserve',(R(139,10)-3)/(30*pi_hi)>R(1,1000))
    check('Garding_gradient_reserve',3*S.Integer(10)**24+17/(45*pi_lo)<S.Integer(10)**25)
    check('Garding_zero_reserve',3*S.Integer(10)**26<S.Integer(10)**27)
    z=S.Integer(10)**28
    check('frequency_threshold',S.Integer(10)**25/z+S.Integer(10)**27/z**2<1/(20*pi_hi))
    check('outer_logarithmic_margin',9-R(38,15)-R(19,60)-4-R(7,20)>R(33,20))
    check('auxiliary_positive_shift',24>22)
    check('shift_log_cost',R(24,2)==12)
    # Arbitrary-jet differential identity, not the fixture below.
    y=S.symbols('y',real=True);nu=S.symbols('nu',real=True)
    h=S.Function('h')(y);t=S.Function('t')(y);a=S.exp(h);v=S.sqrt(a)*t
    L=-S.diff(a*S.diff(t,y),y)/a+nu**2*t/a**2
    V=S.diff(h,y,2)/2+S.diff(h,y)**2/4
    Hv=-S.diff(v,y,2)+(V+nu**2/a**2)*v
    zero('general_half_density_frequency_identity',S.simplify(Hv-S.sqrt(a)*L))
    zero('positive_factorization',S.diff(h,y,2)/2+S.diff(h,y)**2/4-V)
    # An explicit even collar shows why the old bounded-potential
    # scalar comparison cannot simply discard the frequency dependence.
    yy,nn=S.symbols('Y n',real=True)
    potential=nn**2*(S.exp(-2*yy**2)-1)
    zero('frequency_potential_second_jet',S.diff(potential,yy,2).subs(yy,0)+4*nn**2)
    check('frequency_potential_not_uniform',S.diff(potential,yy,2).subs(yy,0).coeff(nn,2)!=0)
    # Three-anchor signs; lower bounds alone do not control the middle node.
    ss=S.symbols('s');nodes=[S.Integer(1),S.Integer(2),S.Integer(3)]
    weights=[S.prod((ss-b)/(a-b) for b in nodes if b!=a) for a in nodes]
    zero('anchor_partition_unity',sum(weights)-1)
    check('middle_anchor_requires_upper_bound',weights[0].subs(ss,4)>0 and weights[1].subs(ss,4)<0 and weights[2].subs(ss,4)>0)


def matrix_log_fixture():
    rng=np.random.default_rng(PROTOCOL['matrix_seed']);U,_=np.linalg.qr(rng.normal(size=(5,5)))
    D2=np.diag([1.,3.,7.,15.,31.])
    V=np.diag([-.6,-.2,.1,.3,.4])+np.diag([.1,.12,.2,.15],1)+np.diag([.1,.12,.2,.15],-1)
    A=U@np.diag([.4,.7,1.1,1.7,2.2])@U.T
    def logm(M):
        e,Q=np.linalg.eigh(M);assert min(e)>0
        return (Q*np.log(e))@Q.T
    logD=logm(D2)
    check('fixture_noncommuting_potential',np.linalg.norm(D2@A-A@D2)>1e-3)
    check('unshifted_log_shortcut_fails',np.linalg.eigvalsh(logm(D2+V)-logD)[0]<-.01)
    for nu in PROTOCOL['matrix_frequencies']:
        H=D2+V+nu**2*A;w0=rng.normal(size=5);w=H@w0
        logH=logm(H);shift=logm(H+24*np.eye(5));inv=np.linalg.inv(H)
        m1=np.linalg.eigvalsh(shift-logD)[0]
        m2=np.linalg.eigvalsh(24*inv-(shift-logH))[0]
        m3=np.linalg.eigvalsh(logH-logD+24*inv)[0]
        bound=.5*w@(logH-logD)@w+12*w0@H@w0
        for label,value in [('shift_monotonicity',m1),('spectral_shift_cost',m2),('combined_operator_bound',m3),('range_vector_bound',bound)]:
            check(label+'_'+str(nu),value>=-1e-10,value)
        NUMERICS.append({'case':'noncommuting_log','nu':nu,'margins':[m1,m2,m3,bound]})


def two_dimensional_fixture():
    # Periodic orbit fixture: evaluate the full two-dimensional Laplacian,
    # then independently the radial half-density operator by FFT.
    for N in PROTOCOL['torus_grids']:
        y=2*np.pi*np.arange(N)/N;u=2*np.pi*np.arange(N)/N
        h=.03*np.cos(y)+.02*np.sin(2*y);hp=-.03*np.sin(y)+.04*np.cos(2*y)
        hpp=-.03*np.cos(y)-.08*np.sin(2*y);a=np.exp(h)
        t=.3*np.cos(2*y)+.1*np.sin(3*y);tp=-.6*np.sin(2*y)+.3*np.cos(3*y)
        tpp=-1.2*np.cos(2*y)-.9*np.sin(3*y)
        k=np.fft.fftfreq(N,1/N);v=np.sqrt(a)*t
        vpp=np.fft.ifft(-k*k*np.fft.fft(v)).real
        vp=np.fft.ifft(1j*k*np.fft.fft(v)).real
        V=hpp/2+hp*hp/4
        for nu in PROTOCOL['torus_frequencies']:
            co=np.cos(nu*u);si=np.sin(nu*u)
            T=t[:,None]*co;TY=tp[:,None]*co;TU=-nu*t[:,None]*si
            lap=(tpp+hp*tp)[:,None]*co-nu**2/a[:,None]**2*T
            J0=2*np.mean(a[:,None]*T*T)
            J1=2*np.mean(a[:,None]*(TY*TY+TU*TU/a[:,None]**2))
            J2=2*np.mean(a[:,None]*lap*lap)
            Hv=-vpp+(V+nu**2/a**2)*v
            targets=[np.mean(v*v),np.mean(vp*vp+(V+nu**2/a**2)*v*v),np.mean(Hv*Hv)]
            errors=np.abs(np.array([J0,J1,J2])-targets)
            for j,e in enumerate(errors):check('torus_norm_'+str(j)+'_'+str(nu)+'_'+str(N),e<1e-10*(1+abs(targets[j])),e)
            z=nu**2/np.max(a)**2
            check('torus_mode_J1_bound_'+str(nu)+'_'+str(N),J1<=J2/z+1e-10)
            check('torus_mode_J0_bound_'+str(nu)+'_'+str(N),J0<=J2/z**2+1e-10)
            NUMERICS.append({'case':'two_dimensional_norms','N':N,'nu':nu,'J':[J0,J1,J2],'errors':errors.tolist()})

def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S168 canonical frequency core bound','input_count':len(ps)}),flush=True)
    try:
        predecessors();budgets_and_identities();matrix_log_fixture();two_dimensional_fixture()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S168 canonical frequency core bound','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'seed_canonical_frequency_core_bound_derived':good,'auxiliary_uniform_frequency_log_bound_derived':good,
         'actual_uniform_logarithmic_comparison_proved':False,'nonlinear_exterior_frequency_bound_proved':False,
         'numerics':NUMERICS,'actual_positive_frequency_operator_values_evaluated':False,'actual_causal_inverse_proved':False,
         'actual_charged_retarded_transfer_proved':False,'actual_wormhole_instability_proved':False,
         'order_reduction_adopted':False,'S167_complete':False,'S168_complete':False,
         'S077_complete':False,'goal_complete':False}
    if not good:
        run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text(encoding="utf-8")
        run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text(encoding="utf-8")
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
