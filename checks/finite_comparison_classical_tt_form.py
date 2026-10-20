"""S167: classical exponential TT form throughout the actual finite comparison.

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
STEM = "finite-comparison-classical-tt-form"
OWN = ["checks/finite_comparison_classical_tt_form.py", "notes/"+STEM+".md"]
PREDECESSORS = [("finite-comparison-disk-return-form", 291),
                ("stationary-global-tt-reduction", 44),
                ("stationary-seed-quantum-form", 53)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "fixed I0 classical plus reference exponential TT Hessian on actual finite path",
    "candidate_energy_lower_bound": ".99", "mutual_energy_ratio_bound": "1.002",
    "LLL_disk_candidate_energy_lower_bound": "8e-24",
    "actual_TT_lift_L2_difference_bound": "7.2e-6 times scalar L2 norm",
    "stationarity_at_finite_comparison_assumed": False,
    "classical_disk_contact_numerically_bounded": False,
    "finite_comparison_classical_disk_Hessian_positive": False,
    "finite_comparison_full_TT_Hessian_positive": False,
    "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
    "nonlinear_family_preserved": False, "signed_compatibility_evaluated": False,
    "fixture_scope": "exact action/geometry identities and rational budgets only",
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
    expr = S.cancel(expr)
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
    for p in ["finite-comparison-disk-return-form", "finite-comparison-york-gap", "finite-comparison-scalar-resolvent",
              "stationary-global-tt-reduction", "stationary-seed-quantum-form",
              "stationary-curved-quantum-tt", "stationary-lll-quantum-exterior",
              "stationary-completed-metric-response-norm"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


def canonical_gradient_and_lift():
    t,v,C,kB,x=S.symbols('t v C kB x',positive=True)
    F=S.Function('F')(t)
    N=S.sqrt(C)*S.sqrt(F)/(kB*t)
    dzdt=S.sqrt(C)/(t*S.sqrt(F))
    primitive=C/kB*(1/v-1/t)
    zero('canonical_warp_measure',N*dzdt-C/(kB*t*t))
    zero('regular_Volterra_primitive',S.diff(primitive,t)-N*dzdt)
    zero('regular_Volterra_bolt_value',primitive.subs(t,v))
    zero('canonical_Volterra_ratio',primitive/N-S.sqrt(C)*(t-v)/(v*S.sqrt(F)))
    zero('canonical_coordinate_jacobian',dzdt*S.sqrt(F/x)-S.sqrt(C/x)/t)
    Ay,Ry,J,uz,uzz,Fu=S.symbols('A R J uz uzz Fu')
    uyy=J**2*uzz-Ry*J*uz
    zero('conformal_Laplacian_cancellation',
         (uyy+Ay*J*uz).subs(uzz,Fu-(Ay-Ry)*uz/J)-J**2*Fu)
    V=S.Function('V')(t); ell=S.Symbol('ell',real=True)
    tau=(t-1)/t; nu=tau*V; kk=S.exp(ell*V)
    TT=1+tau*(kk-1); BB=kk-1-ell*V*kk
    zero('radius_T_derivative',S.diff(TT,t)-ell*kk*S.diff(nu,t)-BB*S.diff(tau,t))
    zero('linear_radius_derivative',t*S.diff(nu,t)-(t-1)*S.diff(V,t)-V/t)
    y=S.Symbol('y',positive=True)
    a,r,f=[S.Function(n)(y) for n in ['a','r','f']]
    A=S.diff(a,y)/a; RR=S.diff(r,y)/r; m=a*r*r
    zero('D0_adjoint',-S.diff(m,y)/m+2*(A+RR)-A)
    zero('D0_real_part',2*(A+RR)-S.diff(m,y)/(2*m)-(R(3,2)*A+RR))
    zero('D0_homogeneous_solution',S.diff(1/(a*a*r*r),y)+2*(A+RR)/(a*a*r*r))
    zero('D0_adjoint_homogeneous_solution',-S.diff(a,y)+A*a)
    d0=S.diff(f,y)+2*(A+RR)*f
    pot=4*(A+RR)**2-2*(S.diff(A+RR,y)+(A+2*RR)*(A+RR))
    zero('D0_energy_integration',m*(d0*d0-S.diff(f,y)**2-pot*f*f)
         -S.diff(2*m*(A+RR)*f*f,y))
    aa,rr,ff=S.symbols('A R Frr')
    expanded=4*(aa+rr)**2-2*(ff/2-aa*aa+aa*rr-rr*rr+(aa+2*rr)*(aa+rr))
    zero('D0_seed_energy_potential',expanded-(4*aa*aa+2*rr*rr-ff))
    b0,be,c,p,q=S.symbols('b0 bell c p q')
    dell=-2*(aa+rr+p)*be-2*(aa-rr+q)*c
    dseed=-2*(aa+rr)*b0-2*(aa-rr)*c
    zero('difference_of_actual_TT_lifts',dell-dseed+2*(aa+rr)*(be-b0)+2*p*be+2*q*c)


def exponential_action_and_bolt():
    e,b,c,by,cy,A,RR,P=S.symbols('e b c by cy A R P',real=True)
    K=S.exp(-e*b)*((RR+e*cy/2)**2
                       +2*(A-e*(by+2*cy)/2)*(RR+e*cy/2))
    qpre=R(3,2)*cy**2+by*cy+2*b*(A-RR)*cy-2*RR*b*by-b*b*(RR*RR+2*A*RR)-P*c*c
    zero('direct_exponential_kinetic_second_variation',-S.diff(K,e,2).subs(e,0)-(qpre+P*c*c))
    X=S.Symbol('X',positive=True); V=S.Function('V')
    Pe=(V(X)-X*S.diff(V(X),X)+X*X*S.diff(V(X),X,2))/(2*X)
    action_potential=S.exp(-e*c)*V(X*S.exp(e*c))/(2*X)
    zero('direct_exponential_potential_second_variation',S.diff(action_potential,e,2).subs(e,0).doit()-Pe*c*c)
    zero('tracefree_exponential_volume',S.exp(e*b/2)*S.exp(-e*(b+2*c)/2)*S.exp(e*c)-1)
    y=S.Symbol('y',positive=True)
    a,r,bf,cf,pf=[S.Function(n)(y) for n in ['a','r','b','c','P']]
    AA=S.diff(a,y)/a; RRR=S.diff(r,y)/r; m=a*r*r
    E1=S.diff(AA-RRR,y)+AA*(AA-RRR)
    pre=qpre.subs({b:bf,c:cf,by:-2*(AA+RRR)*bf-2*(AA-RRR)*cf,
                   cy:S.diff(cf,y),A:AA,RR:RRR,P:pf})
    post=R(3,2)*S.diff(cf,y)**2-4*RRR*bf*S.diff(cf,y)
    post+=(E1+2*AA*RRR-2*RRR**2-pf)*cf**2+(2*AA*RRR+3*RRR**2)*bf**2
    post+=4*(AA*RRR-RRR**2)*bf*cf
    zero('nonstationary_TT_integration_by_parts',m*(pre-post)+S.diff(m*(AA-RRR)*cf**2,y))
    # Relative tensor exponential acceleration is h^2, including disk k.
    h=S.diag(b,-b-2*c,c,c)
    accelerations=S.diag(*[S.diff(S.exp(e*h[i,i]),e,2).subs(e,0) for i in range(4)])
    check('exponential_tensor_acceleration',all(S.expand(z)==0 for z in accelerations-h*h))
    rh,cbolt,G=S.symbols('rh cbolt G',positive=True)
    entropy=-S.pi*rh**2*S.exp(e*cbolt)/G
    zero('bolt_entropy_second_variation_cancellation',S.diff(entropy,e,2).subs(e,0)+S.pi*rh**2*cbolt**2/G)
    zero('regular_bolt_equal_radial_time_exponents',(b-(-b-2*c)).subs(b,-c))
    # Algebraic composition of the retained diffeomorphism Hessian identity.
    Ha,Hk,gradLie,gradh2,gradk2=S.symbols('Ha Hk gradLie gradh2 gradk2')
    zero('exact_off_shell_disk_contact',((Hk+gradk2)-(Ha+gradh2)).subs(Hk,Ha+gradLie)
         -(gradLie+gradk2-gradh2))


def potential_seed_and_finite_geometry():
    r,B,lam,q0,d=S.symbols('r B lambda q0 d',positive=True)
    VV=2+6*lam*r*r-2*q0/r**2+2*d*(-B*B/r**2+4*B**3/r**3-3*B**4/r**4)
    dx=lambda f:S.diff(f,r)/(2*r)
    PP=(VV-r*r*dx(VV)+r**4*dx(dx(VV)))/(2*r*r)
    target=1/r**2-4*q0/r**4+d*(-4*B**2/r**4+25*B**3/r**5-27*B**4/r**6)
    zero('fixed_reference_potential_P',PP-target)
    zero('fixed_reference_potential_log_derivative',r*S.diff(PP,r)
         -(-2/r**2+16*q0/r**4+d*(16*B**2/r**4-125*B**3/r**5+162*B**4/r**6)))
    F=S.Function('F')(r); c=S.Function('c')(r); b=S.Function('b')(r)
    seedV=2*(F+r*S.diff(F,r))
    seedP=(seedV-r*r*dx(seedV)+r**4*dx(dx(seedV)))/(2*r*r)
    zero('seed_P_identity',seedP-(F/r**2-S.diff(F,r)/(2*r)+r*S.diff(F,r,3)/4))
    qexp=R(3,2)*F*S.diff(c,r)**2-4*F/r*b*S.diff(c,r)
    qexp+=(S.diff(F,r,2)/2-F/r**2-seedP)*c*c
    qexp+=(S.diff(F,r)/r+3*F/r**2)*b*b+(2*S.diff(F,r)/r-4*F/r**2)*b*c
    W=3*F-r*S.diff(F,r)/2+r*r*S.diff(F,r,2)/2-r**3*S.diff(F,r,3)/4
    qold=R(3,2)*r*r*F*(S.diff(c,r)+2*c/r)**2+W*c*c-4*F*b*c
    boundary=-r*F*(b+c)*(b+3*c)
    brule={S.diff(b,r):-(2/r+S.diff(F,r)/F)*b-(S.diff(F,r)/F-2/r)*c}
    zero('independent_seed_TT_form_identity',(r*r*qexp-qold-S.diff(boundary,r)).subs(brule))
    y=S.Symbol('y',positive=True)
    a,r,w,s=[S.Function(n)(y) for n in ['a','r','w','s']]
    A0=S.diff(a,y)/a; R0=S.diff(r,y)/r
    dz=lambda f:S.exp(-w)*S.diff(f,y)
    ae,re=S.exp(w)*a,S.exp(s)*r
    AA,RR=dz(ae)/ae,dz(re)/re
    e1=dz(AA-RR)+AA*(AA-RR)
    j0=S.diff(A0-R0,y)+A0*(A0-R0)
    zero('finite_E1_exact_cancellation',S.exp(2*w)*e1-j0-S.diff(w-s,y,2)-A0*S.diff(w-s,y))
    zero('finite_R_pullback',S.exp(w)*RR-R0-S.diff(s,y))
    zero('finite_AR_pullback',S.exp(2*w)*AA*RR-(A0+S.diff(w,y))*(R0+S.diff(s,y)))
    zero('finite_physical_measure',S.exp(w)*ae*re**2/(a*r**2)-S.exp(2*w+2*s))
    zero('finite_derivative_measure',S.exp(2*w+2*s)*S.exp(-2*w)-S.exp(2*s))


def rational_budgets():
    eps=R('3e-8'); ex2=1/(1-2*eps); ex4=1/(1-4*eps)
    u=R('1.2e-8'); fu=R('5e-8')
    check('curvature_source_budget',2*u+R('5.2e-8')/2<=fu)
    check('canonical_u_gradient',fu**2*R('6.1')/R('.95')<R('1.27e-7')**2)
    check('sqrt_point95_lower',R('.97')**2<R('.95'))
    check('physical_u_gradient',fu*R('6.1')/R('.97')<R('3.2e-7'))
    check('conformal_Laplacian_budget',R('6.1')*fu<=R('3.05e-7'))
    check('sqrt5_upper',5<R('2.24')**2)
    check('nu_gradient',R('2.24')*(R('3e-7')+u)<R('7e-7'))
    check('polar_nu_gradient',11*(R('3e-7')+u)<=R('3.432e-6'))
    sy=R('7.1e-7'); wy=R('1.03e-6'); pp=R('1.74e-6'); qq=R('3.2e-7')
    check('finite_S_gradient',(R('1.001')*R('7e-7')+R('2e-16')*R('2.24'))/R('.999')<sy)
    check('finite_polar_S_gradient',(R('1.001')*R('3.432e-6')+R('2e-16')*11)/R('.999')<R('3.45e-6'))
    check('finite_w_gradient',qq+sy<=wy)
    check('finite_p_gradient',qq+2*sy<=pp)
    check('D0_inverse_constant',R(8,9)<R('.943')**2)
    check('D0_energy_constant',1+21*R(8,9)<20)
    inv=R('.943'); bn=R('2.001'); db=R('7.2e-6')
    check('TT_lift_perturbation_invertible',2*inv*pp<1)
    check('candidate_TT_lift_norm',(2+2*inv*qq)/(1-2*inv*pp)<bn)
    f=2*pp*bn+2*qq
    check('TT_difference_forcing',f<R('7.604e-6'),f)
    check('TT_lift_difference_norm',inv*f<db)
    check('potential_value_bound',12+R(61,1000)<13)
    check('potential_log_derivative_bound',48+R(321,1000)<49)
    check('small_exponential_budget',ex2-1<R('6.1e-8'))
    check('radius_square_change',2*R('2.24')*sy+sy**2<R('3.2e-6'))
    check('AR_product_change',R('3.45e-6')+R('2.24')*wy+wy*sy<R('6e-6'))
    check('potential_pullback_change',49*eps+13*R('6.1e-8')<R('2.3e-6'))
    check('c_square_coefficient_change',R('3.05e-7')+2*R('6e-6')+2*R('3.2e-6')+R('2.3e-6')<R('2.2e-5'))
    check('b_square_coefficient_change',2*R('6e-6')+3*R('3.2e-6')<R('2.2e-5'))
    check('bc_coefficient_change',4*(R('6e-6')+R('3.2e-6'))<R('3.7e-5'))
    check('seed_j_bound',R('10.5')+2*11+5<39)
    check('seed_c_square_coefficient',39+2*11+2*5+13<=84)
    check('seed_b_square_coefficient',2*11+3*5<=37)
    check('seed_bc_coefficient',4*max(11,5)<=44)
    L=4*(sy*bn+R('2.24')*db)
    M=R('2.2e-5')*(1+bn**2)+R('3.7e-5')*bn+37*(bn+2)*db+44*db
    check('form_gradient_difference',L<R('7.03e-5'),L)
    check('form_value_difference',M<R('.00157'),M)
    check('difference_D_square_margin',R('7.03e-5')/2<R('.002')*R(3,2))
    check('difference_C_square_margin',R('.00157')+R('7.03e-5')*(R('4.48')+R('.5'))<R('.002'))
    check('seed_absolute_cross_coefficient',4*R('2.24')*2<=R('17.92'))
    check('seed_absolute_C_square_coefficient',84+37*4+44*2<=320)
    coeffD=3+R('17.92')/2
    coeffC=3*R('4.48')**2+R('17.92')*(R('.5')+R('4.48'))+320
    check('seed_absolute_D_square_bound',coeffD<R(3,2)*500,coeffD)
    check('seed_absolute_C_square_bound',coeffC<500,coeffC)
    check('candidate_absolute_density',500+R('.002')<501)
    check('total_classical_form_change',R('.002')+501*R('6.1e-8')<R('.003'))
    delta=2*sy
    check('energy_forward_derivative',ex2*R('1.001')<R('1.002'))
    check('energy_forward_value',ex2*(ex2+R('1501.5')*delta**2)<R('1.002'))
    check('energy_reverse_derivative',ex2*R('1.001')<R('1.002'))
    check('energy_reverse_value',ex4*(1+R('1501.5')*delta**2)<R('1.002'))
    check('positive_classical_candidate_energy',R('.997')/R('1.002')>R('.99'))
    check('seed_radial_gradient_sqrt',5200<R('72.12')**2)
    check('sqrt40_upper',40<R('6.325')**2)
    check('candidate_radial_gradient',R('72.12')+R('6.325')*R('7.604e-6')<73)
    check('candidate_tensor_value',bn**2+(bn+2)**2+2<23)
    check('candidate_radial_H1_energy',ex2*73**2+ex4*23<5400)
    check('H1_to_candidate_energy',5400*R('1.002')<5500)
    check('LLL_disk_candidate_energy',R('1.3e-27')*5500<R('8e-24'))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S167 finite comparison classical exponential TT form",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        canonical_gradient_and_lift()
        exponential_action_and_bolt()
        potential_seed_and_finite_geometry()
        rational_budgets()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S167 finite comparison classical exponential TT form",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "classical_exponential_TT_candidate_energy_positive": good,
        "actual_TT_lift_and_candidate_energy_comparison": good,
        "LLL_disk_lower_bound_on_candidate_energy": good,
        "exact_off_shell_disk_contact_identified": good,
        "classical_disk_contact_numerically_bounded": False,
        "finite_comparison_classical_disk_Hessian_positive": False,
        "finite_comparison_full_TT_Hessian_positive": False,
        "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
        "signed_compatibility_evaluated": False, "nonlinear_family_preserved": False,
        "useful_Newton_certificate": False, "S167_complete": False, "goal_complete": False,
        "unresolved": ["classical off-shell disk contact", "massive/Weyl/reference TT and mixed/density forms", "outer position-log weights",
                       "operator-domain mapping and candidate range", "projection direction and cold/higher bounds",
                       "nonlinear preservation and solve", "signed compatibility", "physical requirements"],
    }
    if not good:
        run["failed_checker_source_utf8"] = (ROOT/OWN[0]).read_text(encoding="utf-8")
        run["failed_note_source_utf8"] = (ROOT/OWN[1]).read_text(encoding="utf-8")
    dest = ROOT/"evidence"/(STEM+".json")
    data = json.loads(dest.read_bytes()) if dest.exists() else {"runs": []}
    data["runs"].append(run)
    dest.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": run["status"], "checks": len(ROWS),
                      "unchanged_inputs": before == after}), flush=True)
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
