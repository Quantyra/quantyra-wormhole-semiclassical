"""S167: off-shell contact and classical disk TT positivity at the finite comparison.

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
STEM = "finite-comparison-classical-disk-contact"
OWN = ["checks/finite_comparison_classical_disk_contact.py", "notes/"+STEM+".md"]
PREDECESSORS = [("finite-comparison-classical-tt-form", 305),
                ("finite-comparison-disk-return-form", 291),
                ("stationary-seed-quantum-form", 53)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "fixed I0 classical disk-chart TT form on actual finite path",
    "classical_disk_contact_energy_bound": ".06",
    "classical_disk_energy_lower_bound": ".93",
    "classical_plus_canonical_LLL_disk_lower_bound": ".92",
    "dimensionless_classical_gradient_bound": "1.8e-5 lambda",
    "stationarity_at_finite_comparison_assumed": False,
    "finite_comparison_full_TT_Hessian_positive": False,
    "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
    "nonlinear_family_preserved": False, "signed_compatibility_evaluated": False,
    "fixture_scope": "exact action/curvature/contact identities and rational budgets only",
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
    for p in ["finite-comparison-classical-tt-form", "stationary-centered-constraint-residual", "finite-comparison-disk-return-form", "finite-comparison-york-gap", "finite-comparison-scalar-resolvent",
              "stationary-global-tt-reduction", "stationary-seed-quantum-form",
              "stationary-curved-quantum-tt", "stationary-lll-quantum-exterior",
              "stationary-completed-metric-response-norm"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


def radius_and_curvature_identities():
    y=S.Symbol('y',positive=True)
    a,r,nu,w,s,t,V=[S.Function(n)(y) for n in ['a','r','nu','w','s','t','V']]
    A0=S.diff(a,y)/a; R0=S.diff(r,y)/r
    zero('canonical_radius_Laplacian_conversion',S.diff(nu,y,2)+(A0-R0)*S.diff(nu,y)
         +2*R0*S.diff(nu,y)-(S.diff(nu,y,2)+(A0+R0)*S.diff(nu,y)))
    drhat,dr4,phiL,uu,nn,rr4,nuL,tt,lam=S.symbols('drhat dr4 phiL u nu R40 nuL t lambda')
    deltaL=nuL-2*uu*phiL
    scalar_variation=(drhat-6*deltaL)/tt**2-2*nn*rr4
    radius_rhs=(drhat-tt**2*dr4)/6+2*uu*phiL-nn*tt**2*rr4/3
    zero('linearized_scalar_curvature_radius_equation',(radius_rhs-nuL).subs(dr4,scalar_variation))
    ell=S.Symbol('ell',real=True); tau=(t-1)/t; kk=S.exp(ell*V)
    TT=1+tau*(kk-1); BB=kk-1-ell*V*kk; nv=tau*V
    lhs=S.diff(TT,y,2)+A0*S.diff(TT,y)
    rhs=ell*kk*(S.diff(nv,y,2)+A0*S.diff(nv,y))+BB*(S.diff(tau,y,2)+A0*S.diff(tau,y))
    rhs+=ell**2*kk*tau*S.diff(V,y)**2
    zero('finite_radius_Laplacian_identity',lhs-rhs)
    zero('log_radius_Laplacian_identity',S.diff(S.log(TT),y,2)+A0*S.diff(S.log(TT),y)
         -lhs/TT+(S.diff(TT,y)/TT)**2)
    tt,xx,dd=S.symbols('t x d',positive=True)
    fext=(tt-1)**2/tt**2*(1+xx*(tt**2+2*tt+3)+dd/tt**2)
    zero('finite_radius_quadratic_factor',fext/(xx*tt*(tt-1))
         -(1-1/tt)*(1+2/tt+(3+1/xx)/tt**2+dd/(xx*tt**4)))
    rr4=-12+dd/xx*(4*tt**-5-6*tt**-6)
    zero('seed_four_curvature_log_derivative',tt*S.diff(rr4,tt)-dd/xx*(-20*tt**-5+36*tt**-6))
    dz=lambda f:S.exp(-w)*S.diff(f,y)
    aa,rr=S.exp(w)*a,S.exp(s)*r
    AA,RR=dz(aa)/aa,dz(rr)/rr
    zero('actual_time_warp_curvature',S.exp(2*w)*dz(dz(aa))/aa
         -S.diff(a,y,2)/a-S.diff(w,y,2)-A0*S.diff(w,y))
    zero('actual_radius_warp_curvature',S.exp(2*w)*(dz(RR)+AA*RR)
         -(S.diff(R0,y)+A0*R0+S.diff(s,y,2)+A0*S.diff(s,y)))
    zero('time_curvature_E1_E2',dz(AA)+AA**2-(dz(AA-RR)+AA*(AA-RR))-(dz(RR)+AA*RR))


def action_gradient_and_potential():
    a,a1,a2,r,r1,r2,N,N1,B,lam,q,d=S.symbols('a a1 a2 r r1 r2 N N1 B lambda q d',positive=True)
    VV=2+6*lam*r*r-2*q/r**2+2*d*(-B**2/r**2+4*B**3/r**3-3*B**4/r**4)
    VX=S.diff(VV,r)/(2*r)
    lag=a*r1*r1/N+2*a1*r*r1/N+a*N*VV/2
    jets={a:a1,a1:a2,r:r1,r1:r2,N:N1}
    D=lambda f:sum(S.diff(f,v)*vp for v,vp in jets.items())
    els={v:S.diff(lag,v)-D(S.diff(lag,v1)) for v,v1 in [(a,a1),(r,r1)]}
    en=(-S.diff(lag,N)/(2*a*r*r)).subs(N,1)
    et=(-els[a]/(2*r*r)).subs({N:1,N1:0})
    es=(-els[r]/(4*a*r)).subs({N:1,N1:0})
    A=a1/a; RR=r1/r; K=a2/a; E2=r2/r-r1*r1/r**2+A*RR
    zero('direct_action_radial_gradient',en-(A*RR+RR**2/2-VV/(4*r*r)))
    zero('direct_action_time_gradient',et-(E2-A*RR+R(3,2)*RR**2-VV/(4*r*r)))
    zero('direct_action_sphere_gradient',es-(K+E2+RR**2-VX/2)/2)
    F,F1,F2=S.symbols('F F1 F2')
    seedV=2*(F+r*F1); seedVX=2*F1/r+F2
    zero('radial_seed_gradient_zero',F1/(2*r)+F/(2*r*r)-seedV/(4*r*r))
    zero('time_seed_gradient_zero',(F1/r-F/r**2)-F1/(2*r)+R(3,2)*F/r**2-seedV/(4*r**2))
    zero('sphere_seed_gradient_zero',(F2/2+F1/r-F/r**2+F/r**2-seedVX/2)/2)
    U=VV/(4*lam*r*r); HH=VX/(2*lam)
    zero('potential_U_expansion',U-(R(3,2)+1/(2*lam*r**2)-q/(2*lam*r**4)
         +d/(2*lam)*(-B**2/r**4+4*B**3/r**5-3*B**4/r**6)))
    zero('potential_H_expansion',HH-(3+q/(lam*r**4)+d/lam*(B**2/r**4-6*B**3/r**5+6*B**4/r**6)))
    zero('potential_U_log_derivative',r*S.diff(U,r)-(-1/(lam*r**2)+2*q/(lam*r**4)
         +d/(2*lam)*(4*B**2/r**4-20*B**3/r**5+18*B**4/r**6)))
    zero('potential_H_log_derivative',r*S.diff(HH,r)-(-4*q/(lam*r**4)
         +d/lam*(-4*B**2/r**4+30*B**3/r**5-36*B**4/r**6)))


def combined_contact_and_polar_identities():
    y=S.Symbol('y',positive=True)
    a,r,eta,b,c=[S.Function(n)(y) for n in ['a','r','eta','b','c']]
    A=S.diff(a,y)/a; RR=S.diff(r,y)/r
    metric=[S.Integer(1),a*a,r*r,r*r]
    h=[b,-b-2*c,c,c]
    ss=[S.diff(eta,y),A*eta,RR*eta,RR*eta]
    for i in range(4):
        v=metric[i]*(2*h[i]-2*ss[i])
        lie=eta*S.diff(v,y)/metric[i]+(2*S.diff(eta,y)*v/metric[i] if i==0 else 0)
        z=lie+(h[i]-2*ss[i])**2-h[i]**2
        zero('combined_contact_component_'+str(i),z-2*eta*(S.diff(h[i],y)-S.diff(ss[i],y)))
    aa,rr,kk,ee,et,bv,cv,cp=S.symbols('A R K E2 eta b c cp')
    bp=-2*(aa+rr)*bv-2*(aa-rr)*cv
    etap=aa*et+bv+cv
    etap2=(kk-aa*aa)*et+aa*etap+bp+cp
    zero('regular_return_second_derivative',etap2-(kk*et-aa*(bv+cv)-2*rr*(bv-cv)+cp))
    zero('regular_time_polar_derivative',(kk-aa*aa)*et+aa*etap-(kk*et+aa*(bv+cv)))
    zero('regular_sphere_polar_derivative',(ee-aa*rr)*et+rr*etap-(ee*et+rr*(bv+cv)))
    d1,d1p,uu=S.symbols('d1 d1p U')
    expr=(2*aa*aa+aa*rr-kk/2)*d1*d1-aa*d1*uu
    derivative_half=(kk+2*aa*rr)*d1*d1/2+aa*d1*d1p
    zero('polar_TT_integration_identity',(expr+derivative_half).subs(uu,d1p+2*(aa+rr)*d1))
    zero('TT_tensor_point_norm',bv*bv+(-bv-2*cv)**2+2*cv*cv-2*(bv+cv)**2-4*cv*cv)
    mat=S.Matrix([[2,2],[2,6]])
    for label,vec,bound in [('difference',[1,-1],R(3,2)),('sum',[1,1],R(1,2))]:
        vv=S.Matrix(vec)
        zero('TT_dual_norm_'+label,(vv.T*mat.inv()*vv)[0]-bound)
    # Independent regular bolt jets check that the contact has zero bolt value.
    aa1,rr1,e1,e3,c0,c2,b2=S.symbols('a1 r1 eta1 eta3 c0 c2 b2')
    AA=1/y+aa1*y; RRR=rr1*y; ET=e1*y+e3*y**3
    cc=c0+c2*y*y; bb=-c0+b2*y*y
    hh=[bb,-bb-2*cc,cc,cc]
    sv=[S.diff(ET,y),AA*ET,RRR*ET,RRR*ET]
    for i in range(4):
        z=S.expand(2*ET*S.diff(hh[i]-sv[i],y))
        zero('regular_bolt_contact_zero_'+str(i),z.subs(y,0))
        check('regular_bolt_contact_even_'+str(i),S.expand(z.subs(y,-y)-z)==0)


def rational_budgets():
    eps=R('3e-8'); ex2=1/(1-2*eps); ex4=1/(1-4*eps)
    check('seed_radius_E20_bound',max(2*11,5)<=22)
    check('seed_four_curvature_bound',12+R(10,1000)<=R('12.01'))
    check('seed_four_curvature_derivative_bound',R(56,1000)<=R('.056'))
    rhs=R('6.1')*R('5.2e-8')/6+(R('6.1e-20')+R('1.2e-8')*R('.056'))/6
    rhs+=2*R('1.2e-8')*22+R('1.2e-8')*R('12.01')/3
    check('linear_radius_equation_force',rhs<R('6.3e-7'),rhs)
    check('linear_radius_orbit_Laplacian',rhs+R('2.24')*R('7e-7')<R('2.3e-6'))
    check('finite_radius_quadratic_seed_bound',1+2+3+R(2,1000)<7)
    check('tau_orbit_Laplacian_bound',2*max(11,5)<=22)
    finite_lap=(R('1.001')*R('2.3e-6')+R('2e-16')*22+R('1.001')*7*R('3e-7')**2)/R('.999')+R('7.1e-7')**2
    check('finite_radius_orbit_Laplacian',finite_lap<R('2.31e-6'),finite_lap)
    check('actual_time_curvature_bound',ex2*(R('10.5')+R('3.05e-7')+R('2.31e-6'))<R('10.501'))
    check('actual_radius_curvature_bound',ex2*(22+R('2.31e-6'))<R('22.001'))
    check('actual_time_warp_positive',R('.7')**2<R(1,2) and R('.7')>R('1.03e-6'))
    check('actual_radius_warp_positive',1-R('3e-7')>0)
    check('potential_U_value',3+R(5,1000)<4)
    check('potential_U_derivative',6+R(24,1000)<7)
    check('potential_H_value',6+R(14,1000)<7)
    check('potential_H_derivative',12+R(74,1000)<13)
    check('U_finite_change',7*eps+4*R('6.1e-8')<R('4.6e-7'))
    check('H_finite_change',13*eps+7*R('6.1e-8')<R('8.2e-7'))
    check('radial_gradient_bound',ex2*(R('6e-6')+R('3.2e-6')/2+R('4.6e-7'))<R('8.2e-6'))
    check('time_gradient_bound',ex2*(R('2.31e-6')+R('6e-6')+R(3,2)*R('3.2e-6')+R('4.6e-7'))<R('1.36e-5'))
    check('sphere_gradient_bound',ex2*(R('3.05e-7')+2*R('2.31e-6')+R('3.2e-6')+R('8.2e-7'))/2<R('4.5e-6'))
    check('four_component_classical_gradient',R('8.2e-6')**2+R('1.36e-5')**2+2*R('4.5e-6')**2<R('1.8e-5')**2)
    check('actual_measure_TT_lift',ex4*R('2.001')<R('2.002'))
    check('actual_tensor_value',R('2.002')**2+R('4.002')**2+2<23 and 23<R('4.8')**2)
    check('actual_tensor_radial_gradient',5500<75**2)
    check('actual_return_vector',R('1.006')*R('4.8')<R('4.83'))
    check('energy_D_coefficient',R(2,3)<R('.817')**2)
    check('polar_time_curvature_sqrt',R('10.501')/2<R('2.292')**2)
    check('polar_TT_product_bound',R('.817')/2+R('2.25')+R('2.292')*R('3.002')<R('9.6'))
    check('scalar_radial_gradient',R('.817')+R('4.5')<R('5.32'))
    check('return_curvature_vector',2*R('10.501')**2+2*R('22.001')**2<R('34.5')**2)
    check('sqrt2_reserve',2<R('1.415')**2)
    check('sqrt7_reserve',7<R('2.646')**2)
    sp=R('34.5')*R('4.83')+R('1.415')*R('9.6')+R('2.646')*R('2.25')*R('4.8')+R('5.32')
    check('return_diagonal_radial_derivative',sp<215,sp)
    check('combined_contact_L1_energy',2*R('4.83')*(75+215)<3000)
    check('classical_contact_energy_margin',R('1.8e-5')*3000<R('.06'))
    check('classical_disk_positive_margin',R('.99')-R('.06')>=R('.93'))
    check('classical_plus_LLL_disk_margin',R('.93')-R('8e-24')>R('.92'))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S167 finite comparison classical disk contact and positivity",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        radius_and_curvature_identities()
        action_gradient_and_potential()
        combined_contact_and_polar_identities()
        rational_budgets()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S167 finite comparison classical disk contact and positivity",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "classical_gradient_uniformly_bounded": good,
        "classical_disk_contact_numerically_bounded": good,
        "finite_comparison_classical_disk_Hessian_positive": good,
        "classical_plus_canonical_LLL_disk_Hessian_positive": good,
        "finite_comparison_full_TT_Hessian_positive": False,
        "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
        "signed_compatibility_evaluated": False, "nonlinear_family_preserved": False,
        "useful_Newton_certificate": False, "S167_complete": False, "goal_complete": False,
        "unresolved": ["remaining massive/Weyl/reference TT and mixed/density forms", "outer position-log weights",
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
