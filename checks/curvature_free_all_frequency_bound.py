"""S168: curvature-free all-frequency bound.

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
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
STEM = "curvature-free-all-frequency-bound"
OWN = ["checks/curvature_free_all_frequency_bound.py", "notes/"+STEM+".md", "evidence/curvature-free-frequency-source-review.json"]
PREDECESSORS = [("outer-coupled-scalar-separation", 351)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope":"Uniform curvature-free normal positivity and positive-axis Friedrichs inverse",
    "quadrature_grids":[32,64],"frequencies":[0,.5,2,8],
    "quadrature_identity_tolerance":1e-11,"quadrature_refinement_tolerance":1e-10,
    "profile":"(eta-1)^3(2-eta)^3 on [1,2], zero outside",
    "energy_equivalence_lower":"1/81","energy_equivalence_upper":"1/9",
    "normal_inverse_bound":"150pi/(19 kappa^4)",
    "retained_fast_factor_removed":False,
    "fixtures_are_exterior_spectral_data":False,
    "actual_wormhole_instability_proved":False,"actual_causal_inverse_proved":False,
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










def weighted_identities():
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('primary_source_reuse_scope',review['story']=='S168' and len(review['primary_sources'])==1 and not review['actual_wormhole_instability_proved'])
    y=S.symbols('y',positive=True);k=S.symbols('k',nonnegative=True)
    v=S.Function('v')(y);F=S.Function('F')(y)
    A=lambda f:y*y*S.diff(f,y,2)-2*y*S.diff(f,y)-4*f-k*k*y*y*f
    B=lambda f:y*y*S.diff(f,y,2)-12*f-k*k*y*y*f
    C=lambda f:2*y*S.diff(f,y)-8*f
    zero('B_equals_negative_S_plus_C',B(v)-A(v)-C(v))
    zero('C_weighted_pairing',y**-4*v*C(v)+5*y**-4*v*v-S.diff(y**-3*v*v,y))
    u0=-2*F/y
    zero('zero_frequency_resolvent_source',(-A(u0)).subs(k,0)-C(S.diff(F,y)))
    vv=S.diff(F,y)
    zero('J0_exact_with_tail',y**-4*C(vv)*u0-4*y**-4*vv*vv+4*S.diff(y**-4*F*vv,y))
    square=y**-2*(S.diff(v,y)-R(3,2)*v/y)**2
    zero('weighted_Hardy_square',square-y**-2*S.diff(v,y)**2+R(9,4)*y**-4*v*v+R(3,2)*S.diff(y**-3*v*v,y))
    zero('primitive_Hardy_integration',y**-4*F*F-R(2,3)*y**-3*F*S.diff(F,y)+S.diff(y**-3*F*F,y)/3)
    zero('primitive_tail_cost',4*R(4,9)-R(16,9))
    w=-S.diff(v,y,2)+k*k*v
    density=y**-4*B(v)**2-w*w-24*y**-2*(S.diff(v,y)**2+k*k*v*v)-72*y**-4*v*v
    zero('all_frequency_curvature_norm',density+24*S.diff(y**-2*v*S.diff(v,y)+y**-3*v*v,y))
    # Explicit norm algebra; the two lower estimates must both be used.
    G,X,P,J,E=S.symbols('G X P J E',real=True)
    lower1=G+X-4*P;lower2=G-R(7,9)*X
    zero('averaged_two_lower_bounds',(lower1+lower2)/2-(G-2*P+X/9))
    zero('Hardy_substituted_uniform_bound',((lower1+lower2)/2).subs(P,4*G/9)-(G+X)/9)
    zero('uniform_lower_constant',R(1,9)/9-R(1,81))
    zero('upper_constant_from_J0',(G+X-4*P+J).subs(J,4*P)/9-(G+X)/9)
    for label,subs in [('first',{G:R(9,4),P:1,X:0}),('second',{G:R(9,4),P:1,X:10})]:
        expr=lower1 if label=='first' else lower2
        check('single_lower_bound_insufficient_'+label,expr.subs(subs)<0)
    Qlog,W,g=S.symbols('Qlog W g',positive=True)
    direct=Qlog+2*(G+X)/(15*S.pi)-(6/g+11/(30*S.pi))*P
    removed=((W+24*(G+X)+72*P)/9-(60*S.pi/g+R(23,3))*E)/(20*S.pi)
    target=Qlog-W/(180*S.pi)+(3/g+23/(60*S.pi))*(E-2*P)
    zero('exact_nonzero_frequency_cancellation',direct-removed-target)
    zero('uniform_gradient_coefficient',(3/g+23/(60*S.pi))/81-1/(27*g)-23/(4860*S.pi))
    check('logarithmic_uniform_margin',9-R(6,5)-4==R(19,5))
    zero('inverse_constant',1/(R(19,150)/S.pi)-150*S.pi/19)
    # Elementary resolvent comparison for every positive spectral value.
    z=S.symbols('z',nonnegative=True)
    zero('resolvent_difference_positive',1-1/(1+z)-z/(1+z))
    zero('resolvent_difference_upper',z-z/(1+z)-z*z/(1+z))


def stable_green_identities():
    y,k=S.symbols('y k',positive=True)
    ul=15*y/k**3*((1+3/(k*y)**2)*S.sinh(k*y)-3/(k*y)*S.cosh(k*y))
    ur=S.exp(-k*y)*(k*k*y/3+k+1/y)
    A=lambda f:y*y*S.diff(f,y,2)-2*y*S.diff(f,y)-(4+k*k*y*y)*f
    simplify=lambda f:S.simplify(S.expand(f.rewrite(S.exp)))
    zero('regular_Green_solution',simplify(A(ul)))
    zero('decaying_Green_solution',simplify(A(ur)))
    W=simplify(ul*S.diff(ur,y)-S.diff(ul,y)*ur)
    zero('normalized_Green_Wronskian',W+5*y*y)
    zero('weighted_Green_delta_normalization',-y*y*(W/5)*y**-4-1)
    series=S.series(ul,y,0,7).removeO()
    zero('regular_endpoint_series',series-y**4-k*k*y**6/14)
    zero('decaying_zero_frequency_tail',S.limit(ur,k,0)-1/y)
    Bstar=lambda f:y*y*S.diff(f,y,2)+4*y*S.diff(f,y)-(10+k*k*y*y)*f
    zero('unweighted_endpoint_action_bounded',S.limit(Bstar(series/y**4),y,0)+10)
    zero('zero_frequency_outside_operator_cancels',Bstar(y**-5).subs(k,0))
    # Squared H_k norm implies the positive-axis L2 gap.
    v=S.Function('v')(y)
    zero('H_squared_form_identity',(-S.diff(v,y,2)+k*k*v)**2-S.diff(v,y,2)**2-2*k*k*S.diff(v,y)**2-k**4*v*v+2*k*k*S.diff(v*S.diff(v,y),y))


def green_quadrature_fixtures():
    def jets(x):
        t=x-1
        return t**3-3*t**4+3*t**5-t**6,3*t**2-12*t**3+15*t**4-6*t**5,6*t-36*t**2+60*t**3-30*t**4
    def left_right(k,y,x):
        if k==0:return y**4/(5*x)
        z=k*y
        ul=15*y/k**3*((1+3/z**2)*np.sinh(z)-3/z*np.cosh(z))
        ur=np.exp(-k*x)*(k*k*x/3+k+1/x)
        return ul*ur/5
    previous={}
    for N in PROTOCOL['quadrature_grids']:
        node,weight=np.polynomial.legendre.leggauss(N)
        xx=1+(node+1)/2;wx=weight/2
        yy=1+(xx[:,None]-1)*(node[None,:]+1)/2
        wy=(xx[:,None]-1)*weight[None,:]/2
        vx,dx,ddx=jets(xx);vy,dy,ddy=jets(yy)
        P=np.sum(wx*vx*vx/xx**4);Z=np.sum(wx*vx*vx/xx**2);G0=np.sum(wx*dx*dx/xx**2)
        for k in PROTOCOL['frequencies']:
            kernel=left_right(k,yy,xx[:,None])
            cx=2*xx*dx-8*vx;cy=2*yy*dy-8*vy
            bx=xx**2*(ddx-k*k*vx)-12*vx;by=yy**2*(ddy-k*k*vy)-12*vy
            pair=lambda sx,sy:2*np.sum(wx*sx/xx**4*np.sum(wy*kernel*sy/yy**4,axis=1))
            J=pair(cx,cy);E=pair(bx,by)/9;G=G0+k*k*Z;expected=(G+14*P+J)/9
            err=abs(E-expected);diff=E-2*P;ratio=diff/G
            check('independent_Green_energy_'+str(k)+'_'+str(N),err<1e-11*(1+abs(E)),err)
            check('positive_Green_resolvent_'+str(k)+'_'+str(N),J>=-1e-12 and J<=4*P+1e-12,J)
            check('Green_resolvent_difference_'+str(k)+'_'+str(N),4*P-J<=k*k*16*Z/9+1e-12,4*P-J)
            check('uniform_energy_lower_'+str(k)+'_'+str(N),ratio>=1/81-1e-9,ratio)
            check('uniform_energy_upper_'+str(k)+'_'+str(N),ratio<=1/9+1e-9,ratio)
            if k==0:check('zero_frequency_exact_energy_'+str(N),abs(diff-G0/9)<1e-11,abs(diff-G0/9))
            if k in previous:
                check('Green_grid_refinement_'+str(k),max(abs(E-previous[k][0]),abs(J-previous[k][1]))<1e-10)
            previous[k]=(E,J)
            NUMERICS.append({'case':'stable_Green_quadrature','N':N,'kappa':k,'E_direct':E,'E_from_C':expected,'J':J,'P':P,'G':G,'energy_ratio':ratio,'identity_error':err})

def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S168 curvature-free all-frequency bound','input_count':len(ps)}),flush=True)
    try:
        predecessors();weighted_identities();stable_green_identities();green_quadrature_fixtures()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S168 curvature-free all-frequency bound','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'normal_curvature_free_all_frequency_bound_derived':good,'positive_axis_normal_Friedrichs_inverse_derived':good,
         'actual_charged_constraints_reduced':False,'actual_complex_frequency_poles_excluded':False,
         'numerical_collar_threshold_evaluated':False,'retained_fast_factor_removed':False,
         'nonlinear_exterior_frequency_bound_proved':False,
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
