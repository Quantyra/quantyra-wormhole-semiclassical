"""S168: exact AdS conformal-vacuum scalar stability preflight.

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
STEM = "ads-fast-scalar-stability-preflight"
OWN = ["checks/ads_fast_scalar_stability_preflight.py", "notes/"+STEM+".md", "evidence/ads-fast-scalar-stability-source-review.json"]
PREDECESSORS = [("cold-projection-source-lower-bound", 321)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "exact AdS4 conformal-vacuum comparison; actual charged retarded transfer open",
    "metric_lift": "delta E_ab[2s g]=-2 epsilon Q_ab(Box+E)s under stated vacuum assumptions",
    "BF_parameter": "nu^2=-60pi/(G lambda)-17/12<0",
    "radial_operator": "A=-d_z^2-a/sin(z)^2; a=60pi/g+5/3>1.125e35",
    "compact_test": "smoothed Rayleigh<-a/(16d^2)<-7e33/d^2",
    "actual_wormhole_instability_proved": False,
    "order_reduction_adopted": False,
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
    for p in ["general-background-newton-operator", "stationary-scalar-phase-moment",
              "on-shell-adjoint-error-and-trace-gradient", "on-shell-finite-range-response",
              "cold-projection-source-lower-bound"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)



def sources_and_trace():
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('primary_source_scope',review['story']=='S168' and not review['actual_wormhole_instability_proved'] and not review['order_reduction_adopted'])
    check('primary_source_count',len(review['primary_sources'])==2 and review['direct_quotes']==[])
    for source in review['primary_sources']:
        check('primary_identity_'+source['id'],source['url'].startswith('https://arxiv.org/pdf/') and len(source['sha256'])==64 and source['bytes']>100000)
    lam,eps,G,D,Ric,psi=S.symbols('lambda epsilon G D R psi',nonzero=True)
    lam0=lam+R(11,6)*eps*lam**2
    trace=Ric+12*lam0-R(11,12)*eps*Ric**2/6
    zero('AdS_vacuum_trace_cancellation',trace.subs(Ric,-12*lam))
    zero('Euler_first_variation',S.diff(Ric**2/6,Ric).subs(Ric,-12*lam)+4*lam)
    E=1/eps+R(11,3)*lam
    variation=psi+eps*(R(11,3)*lam*psi+D*psi)
    zero('Lorentzian_trace_factor',variation-eps*(D+E)*psi)
    zero('physical_fast_mass',E.subs(eps,G/(60*S.pi))-(60*S.pi/G+R(11,3)*lam))
    conformal_R=-6*D+24*lam
    zero('conformal_curvature_factor',conformal_R+6*(D-4*lam))
    # Local scalar-to-tensor map: Hess A + metric B. Conservation fixes B.
    aa,bb=S.symbols('a b');A=aa*D+bb;B=-(D-3*lam)*A
    zero('tensor_divergence_constraint',(D-3*lam)*A+B)
    zero('tensor_Q_trace',D*A+4*B+3*(D-4*lam)*A)
    response_A=-2*eps*(D+E)
    zero('tensor_lift_trace',-3*(D-4*lam)*response_A+eps*(D+E)*conformal_R)
    zero('tensor_lift_classical_limit',S.limit(response_A,eps,0)+2)
    zero('branch_nonzero_curvature',conformal_R.subs(D,-E)-6*(E+4*lam))
    zero('branch_full_tensor_zero',response_A.subs(D,-E))
    check('wrong_mass_sign_detected',S.simplify(response_A.subs(D,E))!=0)


def radial_and_compact_test():
    x=S.symbols('x',positive=True);mu=S.symbols('mu',real=True);u=S.Function('u')(x)
    f=u/S.tan(x)
    # Spatial part of Box-m^2 in the displayed global AdS metric.
    spatial=S.diff(S.tan(x)**2*S.diff(f,x),x)/S.tan(x)**2-mu*f/S.cos(x)**2
    want=S.diff(u,x,2)-(mu+2)*u/S.cos(x)**2
    zero('global_AdS_radial_conjugation',S.trigsimp(S.expand_trig(S.tan(x)*spatial-want)))
    g=S.symbols('g',positive=True);Ebar=60*S.pi/g+R(11,3)
    nu2=R(9,4)-Ebar;a=Ebar-2
    zero('BF_parameter',nu2+60*S.pi/g+R(17,12))
    zero('inverse_square_coefficient',a+nu2-R(1,4))
    z,p=S.symbols('z p',positive=True)
    # u=z^(p-1), psi~z^p; the indicial equation is p(p-3)+Ebar=0.
    zero('scalar_indicial_equation',-(p-1)*(p-2)-a-(-p*(p-3)-Ebar))
    k=S.symbols('k',positive=True)
    for sign in [-1,1]:zero('static_fast_exponent_'+str(sign),(p*(p-3)+k*k+R(9,4)).subs(p,R(3,2)+sign*S.I*k))
    lower=180/R('4e-17')**2
    check('benchmark_a_lower',lower==R('1.125e35') and lower>200)
    check('smooth_witness_reserve',lower/16>R('7e33'))
    v=S.symbols('v',real=True);eta=(v-1)**2*(2-v)**2
    mass=S.integrate(eta**2,(v,1,2));grad=S.integrate(S.diff(eta,v)**2,(v,1,2))
    zero('profile_L2_mass',mass-R(1,630))
    zero('profile_gradient_mass',grad-R(2,105))
    zero('profile_gradient_ratio',grad/mass-12)
    check('profile_minimal_domain',all(S.diff(eta,v,j).subs(v,t)==0 for j in [0,1] for t in [1,2]))
    # Exact inverse-square integral, independently integrated.
    moment=S.integrate(S.expand(eta**2)/v**2,(v,1,2))
    check('profile_inverse_square_lower',bool(S.N(moment-mass/4,60)>0),S.N(moment/mass,25))
    check('benchmark_compact_margin',12-lower/4 < -lower/8)
    check('perturbed_tail_margin',24+1-lower/4 < -lower/8)
    d=S.symbols('d',positive=True)
    ud=d**R(-1,2)*eta.subs(v,z/d)
    zero('shell_L2_scaling',S.simplify(d*ud.subs(z,d*v)**2-eta**2))
    zero('shell_gradient_scaling',S.simplify(d*S.diff(ud,z).subs(z,d*v)**2-S.diff(eta,v)**2/d**2))


def independent_rayleigh_refinement():
    import mpmath as mp
    mp.mp.dps=65
    gg=mp.mpf(1999)/2000*mp.mpf(2000)**-10
    a=60*mp.pi/gg+mp.mpf(5)/3
    def eta(v):return (v-1)**2*(2-v)**2
    def etap(v):return 2*(v-1)*(2-v)*(3-2*v)
    for case,d in enumerate([mp.mpf('.1'),mp.mpf('.05'),mp.mpf('.025')]):
        mass=mp.quad(lambda v:eta(v)**2,[1,2])
        kinetic=mp.quad(lambda v:etap(v)**2,[1,2])/mass
        potential=mp.quad(lambda v:eta(v)**2/mp.sin(d*v)**2,[1,2])/mass
        exact=kinetic/a-d*d*potential
        check('direct_integral_negative_'+str(case),exact < -mp.mpf(1)/8,mp.nstr(exact,30))
        errors=[]
        for n in [32,64,128]:
            step=mp.mpf(1)/n
            values=[eta(1+i*step) for i in range(n+1)]
            norm=sum(v*v for v in values[1:-1])
            # Discrete Laplacian energy with zero boundary values, not analytic eta'.
            kinetic_fd=sum((values[i+1]-values[i])**2 for i in range(n))/(step*step*norm)
            pot_fd=sum(values[i]**2/mp.sin(d*(1+i*step))**2 for i in range(1,n))/norm
            got=kinetic_fd/a-d*d*pot_fd
            err=abs(got-exact);errors.append(err)
            check('finite_difference_negative_'+str(case)+'_'+str(n),got < -mp.mpf(1)/8,mp.nstr(got,30))
        check('rayleigh_refinement_'+str(case),errors[1]<errors[0]/8 and errors[2]<errors[1]/8 and errors[2]<mp.mpf('1e-10'),[mp.nstr(v,12) for v in errors])
        check('negative_potential_omission_'+str(case),kinetic/a>0)
        check('opposite_mass_sign_control_'+str(case),kinetic/a+d*d*potential>0)


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S168 AdS fast scalar stability preflight','input_count':len(ps)}),flush=True)
    try:
        predecessors();sources_and_trace();radial_and_compact_test();independent_rayleigh_refinement()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S168 AdS fast scalar stability preflight','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'vacuum_trace_branch_identified':good,'vacuum_conformal_metric_lift_proved':good,
         'vacuum_radial_operator_unbounded_below':good,'conditional_tail_transfer_criterion_proved':good,
         'actual_charged_retarded_transfer_proved':False,'actual_wormhole_instability_proved':False,
         'order_reduction_adopted':False,'S167_complete':False,'S168_complete':False,
         'S077_complete':False,'goal_complete':False}
    if not good:
        run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text()
        run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text()
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
