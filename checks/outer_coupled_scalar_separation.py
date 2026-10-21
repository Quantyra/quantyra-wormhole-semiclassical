"""S168: outer coupled scalar separation.

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
STEM = "outer-coupled-scalar-separation"
OWN = ["checks/outer_coupled_scalar_separation.py", "notes/"+STEM+".md", "evidence/outer-scalar-separation-source-review.json"]
PREDECESSORS = [("uniform-logarithmic-frequency-comparison", 338)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope":"Actual seed direct-shape annular limit and coupled AdS normal curvature-scalar separation",
    "mpmath_digits":40,"digamma_q":[3,20],"digamma_t":[10000,1000000,100000000],
    "digamma_error_bound":"2(q+1)/t","refinement_ratio":.02,
    "profile":"(eta-1)^3(2-eta)^3 on [1,2], zero outside",
    "profile_approximation_scope":"H3 compact approximation, not a compact exact eigenmode",
    "auxiliary_inverse_gap":"25/4","retained_fast_factor_inverted":False,
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









def local_action_columns():
    y=S.symbols('y',positive=True);g=S.symbols('g',positive=True);k=S.symbols('k',nonnegative=True)
    amp,co,si=S.symbols('amp co si');P0=S.symbols('P0')
    p=S.Function('p')(y);v=S.Function('v')(y)
    P=-S.log(y)+P0;P1=-1/y;P2=1/y**2
    DP=P2+amp*(S.diff(p,y,2)-k*k*p)*co
    DV=amp*(S.diff(v,y,2)-k*k*v)*co
    KP=P1**2+2*amp*P1*S.diff(p,y)*co+amp**2*(S.diff(p,y)**2*co**2+k*k*p*p*si**2)
    gradcross=(P1+amp*S.diff(p,y)*co)*amp*S.diff(v,y)*co+amp**2*k*k*p*v*si**2
    ratio=1+11*g/(360*S.pi)
    FC=(1+2*amp*p*co+2*amp**2*p*p*co**2)*(6*KP+4*gradcross)/y**2
    FC+=6*ratio/y**4*(1+amp*(4*p+2*v)*co+amp**2*(4*p+2*v)**2*co**2/2)
    FW=(1-2*amp*v*co+2*amp**2*v*v*co**2)*((P+amp*p*co)*DV**2/15-DP*DV/15-DP**2/10-DV*KP/15-R(7,90)*DP*KP-R(7,180)*KP**2)
    coef=lambda f,n:S.expand(f).coeff(amp,n).subs({co**2:1,si**2:1,co:1,si:1}).expand()
    Q=-coef(FC,2)/(2*g)-coef(FW,2)/(2*S.pi)
    e=60*S.pi/g+R(11,3)
    D=lambda f:y*y*(S.diff(f,y,2)-k*k*f)-2*y*S.diff(f,y)
    A=lambda f:D(f)-4*f
    B=lambda f:y*y*(S.diff(f,y,2)-k*k*f)-12*f
    EL=lambda f,u:S.diff(f,u)-S.diff(S.diff(f,S.diff(u,y)),y)+S.diff(S.diff(f,S.diff(u,y,2)),y,2)
    target=2/y**4*((D(A(p))+e*A(p))/(20*S.pi)+(D(B(v))+e*B(v))/(60*S.pi))
    zero('independent_local_action_conformal_and_mixed_columns',S.expand(EL(Q,p)-target))
    vv=Q.subs(p,0).doit()
    dv=S.diff(v,y,2)-k*k*v
    expected=-6*ratio/(g*y**4)*v*v-P/(30*S.pi)*dv*dv-2/(15*S.pi*y**2)*v*dv+13/(60*S.pi*y**4)*v*v
    zero('independent_local_action_direct_shape_block',S.expand(vv-expected))
    # The normal first variation is static. Its two local rows vanish
    # with the prescribed lambda0 correction; no residual is omitted.
    first=(-coef(FC,1)/(4*g)-coef(FW,1)/(4*S.pi)).subs(k,0)
    zero('normal_background_conformal_stationarity',S.expand(EL(first,p)))
    zero('normal_background_shape_stationarity',S.expand(EL(first,v)))
    # Scalar curvature from the independent conformal curvature formula.
    # lambda=1 here; the result is multiplied by lambda in the note.
    raw=(1-2*amp*p*co)*y*y*(1-2*amp*v*co)*(-2*DV-6*DP-6*KP)
    dR=coef(raw,1)
    zero('normal_curvature_map',S.expand(dR+2*(3*A(p)+B(v))))
    zero('quantum_weighted_potential',-6*(ratio-1)/g-R(11,60)/S.pi+R(11,30)/S.pi)


def ultraviolet_and_budgets():
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('primary_source_scope',review['story']=='S168' and len(review['primary_sources'])==2 and not review['actual_wormhole_instability_proved'])
    al=S.symbols('alpha',positive=True);a=al*(1-al)
    zero('Feynman_a2',S.integrate(a*a,(al,0,1))-R(1,30))
    zero('Feynman_a2_log',S.integrate(a*a*(S.log(al)+S.log(1-al)),(al,0,1))+R(47,900))
    zero('K_second_constant',1-R(47,30)+R(17,30))
    zero('K_integrated_constant',R(17,30)+R(3,2)-R(31,15))
    zero('momentum_log_constant',R(31,15)/2-R(31,30))
    z=S.symbols('z',positive=True);scale=S.symbols('scale',positive=True)
    zero('twice_integrated_UV_formula',S.diff(z*z*(S.log(z/scale)-R(31,15))/2,z,2)-(S.log(z/scale)-R(17,30)))
    check('benchmark_x_upper',(R(352,100000)*2000**4-1+R(1,10**12))/6<10**10)
    check('log_x_upper',(R(8,3))**47>10**20)
    check('Lbar_upper',R(47,4)+R(1,10**21)<12)
    check('curvature_free_log_margin',9-R(6,5)-4==R(19,5))
    zero('curvature_free_fourth_coefficient',R(19,5)/30-R(19,150))
    check('stable_gap_inverse',R(4,25)<R(1,4))
    eta=S.symbols('eta',positive=True);phi=(eta-1)**3*(2-eta)**3
    moments=[S.integrate(S.diff(phi,eta,j)**2,(eta,1,2)) for j in range(4)]
    expected=[R(1,12012),R(1,770),R(2,35),R(36,7)]
    for j in range(4):zero('profile_derivative_norm_'+str(j),moments[j]-expected[j])
    I0,I1,I2,I3=moments
    bound=((I3/2+R(35,3)*I2)/90+2*I1/45)/I0
    zero('direct_block_quantum_upper',bound-R(97396,225))
    check('quantum_upper_less_than_434',bound<434)
    gmax=(R(4,10**17))**2
    check('negative_direct_block_margin',8*gmax*bound<1)
    zero('logarithmic_quartic_shift',R(31,30)+R(1,6)-R(6,5))


def curvature_separation():
    eta=S.symbols('eta',positive=True);k=S.symbols('k',nonnegative=True)
    F=S.Function('F')(eta);v=S.diff(F,eta)
    D=lambda f:eta**2*S.diff(f,eta,2)-2*eta*S.diff(f,eta)
    A=lambda f:D(f)-4*f
    B=lambda f:eta**2*S.diff(f,eta,2)-12*f
    T=(v+2*F/eta)/3
    zero('explicit_zero_frequency_separation',S.expand(A(T)-B(v)/3))
    zero('auxiliary_Green_Wronskian',eta**4*S.diff(eta**-1,eta)-S.diff(eta**4,eta)*eta**-1+5*eta**2)
    tt=S.symbols('t',real=True);u=S.Function('u')(tt);p=S.exp(3*tt/2)*u
    conj=S.simplify(S.exp(-3*tt/2)*(S.diff(p,tt,2)-3*S.diff(p,tt)-(4+k*k*S.exp(2*tt))*p))
    zero('stable_auxiliary_gap',conj-(S.diff(u,tt,2)-(R(25,4)+k*k*S.exp(2*tt))*u))
    # Explicit polynomial plus its required noncompact stable tail.
    phi=(eta-1)**3*(2-eta)**3;xx=S.symbols('x',positive=True)
    prim=S.integrate(phi.subs(eta,xx),(xx,1,eta));total=S.integrate(phi,(eta,1,2))
    Tin=(phi+2*prim/eta)/3;Tout=2*total/(3*eta)
    Einner=S.integrate(S.diff(Tin,eta)**2/eta**2+4*Tin**2/eta**4,(eta,1,2))
    Etail=S.integrate(S.diff(Tout,eta)**2/eta**2+4*Tout**2/eta**4,(eta,2,S.oo))
    Ig=S.integrate(S.diff(phi,eta)**2/eta**2,(eta,1,2));Ip=S.integrate(phi**2/eta**4,(eta,1,2))
    IB=S.integrate(B(phi)**2/eta**4,(eta,1,2));I2=S.integrate(S.diff(phi,eta,2)**2,(eta,1,2))
    zero('stable_energy_including_tail',S.simplify(Einner+Etail-Ig/9-2*Ip))
    check('omitted_stable_tail_detected',Etail>0,Etail)
    zero('weighted_shape_curvature_norm',S.simplify(IB-I2-24*Ig-72*Ip))
    # Algebra for arbitrary norms, independent of the chosen profile.
    G,P,J,Qlog,g=S.symbols('G P J Qlog g',positive=True)
    direct=Qlog+2*G/(15*S.pi)-(6/g+11/(30*S.pi))*P
    scalar=((J+24*G+72*P)/9-(60*S.pi/g+R(23,3))*(G/9+2*P))/(20*S.pi)
    target=Qlog-J/(180*S.pi)+(1/(3*g)+23/(540*S.pi))*G
    zero('curvature_free_potential_cancellation',direct-scalar-target)
    # Abstract exact square completion uses A T=Bv/3, not a fast inverse.
    aa,ee,pp,bv,vv=S.symbols('a e p b v')
    orig=aa*(aa+ee+4)*pp**2+R(2,3)*(aa+ee+4)*pp*bv+vv
    shifted=aa*(aa+ee+4)*(pp+bv/(3*aa))**2+vv-(aa+ee+4)*bv*bv/(9*aa)
    zero('commuting_factor_square_completion',orig-shifted)


def digamma_fixtures():
    mp.mp.dps=40
    def f(t,q):
        root=mp.sqrt(q*q-4*t)
        return -2*(mp.digamma(1+(q+root)/2)+mp.digamma(1+(q-root)/2))
    for q in PROTOCOL['digamma_q']:
        old=None
        for t0 in PROTOCOL['digamma_t']:
            t=mp.mpf(t0);value=f(t,q);der=mp.diff(lambda x:f(x,q),t)
            errs=[abs(value+2*mp.log(t)),abs(t*der+2)]
            check('digamma_real_'+str(q)+'_'+str(t0),abs(mp.im(value))+abs(mp.im(der))<mp.mpf('1e-30'))
            for j,err in enumerate(errs):
                check('digamma_asymptotic_'+str(j)+'_'+str(q)+'_'+str(t0),err<2*(q+1)/t,mp.nstr(err,16))
                if old is not None:check('digamma_refinement_'+str(j)+'_'+str(q)+'_'+str(t0),err<mp.mpf('.02')*old[j])
            old=errs
            NUMERICS.append({'case':'digamma_UV','q_fixture':q,'t':t0,'errors':[str(x) for x in errs]})

def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S168 outer coupled scalar separation','input_count':len(ps)}),flush=True)
    try:
        predecessors();local_action_columns();ultraviolet_and_budgets();curvature_separation();digamma_fixtures()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S168 outer coupled scalar separation','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'actual_seed_direct_shape_annular_limit_derived':good,'normal_curvature_scalar_separation_derived':good,
         'zero_scaled_frequency_curvature_free_positivity_derived':good,'direct_shape_core_unbounded_below_derived':good,
         'nonzero_scaled_frequency_curvature_free_bound_proved':False,'actual_charged_constraints_reduced':False,
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
