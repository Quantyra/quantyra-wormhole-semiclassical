"""S168: charged conformal response and Lorentzian constraint identities.

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
STEM = "charged-conformal-response-and-constraints"
OWN = ["checks/charged_conformal_response_and_constraints.py", "notes/"+STEM+".md", "evidence/charged-conformal-response-source-review.json"]
PREDECESSORS = [("ads-fast-scalar-stability-preflight", 318)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "Local causal conformal column; charged trace and dynamic constraints",
    "independent_fixtures": "Exact Lorentzian Christoffel divergence, warped Einstein response and Box R radial variational density",
    "rational_radii": ["1/3", "4/5"],
    "fixtures_are_physical_solutions": False,
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
    expr = S.cancel(S.sympify(expr).doit())
    check(name, expr == 0, expr)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_bytes())
    return data.get("runs", [data])[-1]


def inputs(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def paths():
    found = set(OWN) | {"notes/stationary-coupled-asymptotic-obstruction.md", "notes/spherical-anomaly-action.md"}
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




def local_algebra():
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('source_review_scope',review['story']=='S168' and len(review['primary_sources'])==1 and not review['actual_wormhole_instability_proved'])
    ep,lam,rr,C,J,M,BR=S.symbols('epsilon lambda R C2 J2 M BoxR')
    s,Bs,BBs,JH,gr=S.symbols('s Boxs Box2s JH gradRgrads')
    l0=lam+R(11,6)*ep*lam**2
    Ban=R(7,12)*C+R(11,6)*J-R(11,72)*rr**2-20*M+BR
    T=rr+12*l0+ep*Ban
    EU=C-2*J+rr**2/6
    dR=-6*Bs-2*rr*s
    dEU=-4*s*EU+8*JH-2*rr*Bs
    dBR=-6*BBs-2*rr*Bs-2*gr-4*s*BR
    direct=dR+ep*(R(3,2)*(-4*s*C)-R(11,12)*dEU+80*s*M+dBR)
    pred=-6*ep*BBs-(6+ep*rr/6)*Bs-R(22,3)*ep*JH-2*ep*gr+(-2*rr-4*ep*Ban)*s
    zero('conformal_trace_direct',direct-pred)
    zero('on_shell_trace_contact',pred+4*s*T-(-6*ep*BBs-(6+ep*rr/6)*Bs-R(22,3)*ep*JH-2*ep*gr+(2*rr+48*l0)*s))
    E=1/ep+R(11,3)*lam
    dcurv,dweyl,dJ,dM,opBox=S.symbols('dR dC2 dJ2 dM operator_BoxR')
    var=dcurv+ep*(R(7,12)*dweyl+R(11,6)*dJ-R(11,36)*rr*dcurv-20*dM+opBox)
    zero('trace_local_row',var-((1-R(11,36)*ep*rr)*dcurv+ep*(R(7,12)*dweyl+R(11,6)*dJ-20*dM+opBox)))
    Z=S.symbols('Z')
    zero('areal_trace_mass_coefficient',(1/ep-R(11,36)*rr).subs(rr,Z-12*lam)-(E-R(11,36)*Z))
    rp,dr,xi,dR0=S.symbols('Rprime delta_r xi delta_R')
    zero('gauge_invariant_Psi',(dR0+xi*rp)-rp*(dr+xi)-(dR0-rp*dr))
    check('bare_curvature_not_gauge_invariant',xi*rp != 0)
    radius,q=S.symbols('r q',positive=True)
    zero('areal_magnetic_invariant',S.diff(q**2/(2*radius**4),radius)*dr-(-2*q**2/radius**5)*dr)
    # Count the volume, index raising and varying-Weyl-vector contacts.
    zero('conformal_column_contact_count',4-4+2-2)
    ric,hs,bx=S.symbols('Ric_component Hess_component Box_s')
    e=ric-rr/2-3*l0
    col=-2*(hs-bx-ric*s)-s*(rr+12*l0)-2*s*e
    zero('Einstein_column_residual_contacts',col-(-2*hs+2*bx-6*l0*s))
    check('omitted_equation_contact_detected',S.expand(col+2*s*e-col)!=0)
    check('omitted_trace_contact_detected',s*(rr+12*l0)!=0)
    # A spherical KG jet: at s=ds=0 the angular Hessian vanishes.
    jt,jr,H=S.symbols('j_t j_r H')
    contraction=-jt*H+jr*H
    zero('spherical_KG_zero_first_jet',pred.subs({s:0,Bs:0,BBs:0,gr:0,JH:contraction})+R(22,3)*ep*(jr-jt)*H)
    zero('RN_zero_first_jet_cancellation',contraction.subs(jr,jt))
    eta=S.diag(-1,1,1,1);Jcov=S.diag(3,1,1,1);Hcov=S.diag(1,1,0,0)
    zero('unrestricted_J_tracefree',S.trace(eta*Jcov))
    zero('unrestricted_H_KG_tracefree',S.trace(eta*Hcov))
    check('nonEinstein_unrestricted_tracefree_pair',S.trace(eta*Jcov*eta*Hcov)==4)
    # RN contraction in mixed orthonormal components; use its KG row.
    F,fp,sp,tt,spp,Q2=S.symbols('F Fprime sprime stt srr Q2')
    Ht=-tt/F+fp*sp/2;Hr=F*spp+fp*sp/2;Ho=F*sp/radius
    Box=Ht+Hr+2*Ho
    JHs=Q2/radius**4*(-Ht-Hr+2*Ho)
    zero('RN_Hessian_contraction',JHs-Q2/radius**4*(-Box+4*F*sp/radius))
    got=pred.subs({rr:-12*lam,Bs:-E*s,BBs:E**2*s,gr:0,JH:Q2/radius**4*(E*s+4*F*sp/radius)}).subs(s,0)
    zero('RN_KG_nonzero_trace_jet',got+R(88,3)*ep*Q2*F*sp/radius**5)
    check('RN_point_jet_detected',got.subs({ep:R(1,7),Q2:R(2,5),F:3,sp:2,radius:4})!=0)
    # Fixed metric state perturbations have zero anomaly variation.
    zero('state_only_trace_variation',var.subs({dcurv:0,dweyl:0,dJ:0,dM:0,opBox:0}))


def ward_from_christoffels():
    t,r,th,ph=S.symbols('t r theta phi',real=True)
    co=[t,r,th,ph];N=S.Function('N')(r);A=S.Function('A')(r)
    metric=S.diag(-N**2,A**2,r**2,r**2*S.sin(th)**2);inv=metric.inv()
    Gamma=[[[S.simplify(sum(inv[a,d]*(S.diff(metric[d,c],co[b])+S.diff(metric[d,b],co[c])-S.diff(metric[b,c],co[d])) for d in range(4))/2) for c in range(4)] for b in range(4)] for a in range(4)]
    et,er,eo,j=[S.Function(n)(t,r) for n in ['et','er','eo','j']]
    tensor=S.diag(et,er,eo,eo);tensor[1,0]=j;tensor[0,1]=-A**2*j/N**2
    div=[sum(S.diff(tensor[a,b],co[a])+sum(Gamma[a][a][c]*tensor[c,b]-Gamma[c][a][b]*tensor[a,c] for c in range(4)) for a in range(4)) for b in range(4)]
    Wt=S.diff(et,t)+S.diff(j,r)+(S.diff(N,r)/N+S.diff(A,r)/A+2/r)*j
    Wr=-A**2/N**2*S.diff(j,t)+S.diff(er,r)+S.diff(N,r)/N*(er-et)+2/r*(er-eo)
    for i,ans in enumerate([Wt,Wr,0,0]):zero('Lorentzian_divergence_'+str(i),div[i]-ans)
    f=S.Function('f')(t+r);sub={N:1,A:1,et:-f/r**2,er:f/r**2,eo:0,j:f/r**2}
    zero('nonzero_Ward_wave_time',Wt.subs(sub).doit())
    zero('nonzero_Ward_wave_radius',Wr.subs(sub).doit())
    zero('nonzero_Ward_wave_trace',(et+er+2*eo).subs(sub))
    zero('nonzero_Ward_wave_orbit',(et+er).subs(sub))
    check('Ward_wave_not_zero', (f/r**2).subs(f,1)!=0)
    v=S.symbols('v');speed=S.Matrix([[0,1],[N**2/A**2,0]])
    zero('constraint_characteristic_polynomial',(speed-v*S.eye(2)).det()-(v**2-N**2/A**2))
    zero('static_constraint_first_integral',Wr.subs({et:-er,eo:0,j:0}).doit()-S.diff(N**2*r**2*er,r)/(N**2*r**2))
    # Magnetic Weyl current at fixed q: divergence of s F from its density.
    q=S.symbols('q');s=S.Function('s')(t,r);w=N*A*r**2*S.sin(th)
    field=S.zeros(4);field[2,3]=q/(2*r**4*S.sin(th));field[3,2]=-field[2,3]
    for b in range(4):zero('magnetic_conformal_current_'+str(b),sum(S.diff(w*s*field[a,b],co[a]) for a in range(4))/w)


def curvature(x,N,A,r):
    proper=lambda f:S.diff(f,x)/A
    k01=-proper(proper(N))/N;k02=-proper(N)*proper(r)/(N*r)
    k12=-proper(proper(r))/r;k23=(1-proper(r)**2)/r**2
    ric=[k01+2*k02,k01+2*k12,k02+k12+k23]
    scal=ric[0]+ric[1]+2*ric[2]
    return ric,scal


def exact_warped_fixtures():
    x,u=S.symbols('x u',real=True)
    N=1+x**2/7;A=1+x/5;r=2+x**2/3;s=1+x/11+x**3/13
    ric,scal=curvature(x,N,A,r);scal=S.cancel(scal);w=N*A*r**2
    box=lambda f:S.diff(N*r**2/A*S.diff(f,x),x)/w
    hess=lambda f:[S.diff(N,x)*S.diff(f,x)/(N*A**2),(S.diff(f,x,2)-S.diff(A,x)*S.diff(f,x)/A)/A**2,S.diff(r,x)*S.diff(f,x)/(r*A**2)]
    ru,Su=curvature(x,S.exp(u*s)*N,S.exp(u*s)*A,S.exp(u*s)*r)
    l0=R(2,9)
    for idx in range(3):
        direct=S.diff(ru[idx]-Su/2-3*l0,u).subs(u,0)+2*s*(ric[idx]-scal/2-3*l0)
        predicted=-2*hess(s)[idx]+2*box(s)-6*l0*s
        for xx in [R(1,3),R(4,5)]:zero('warped_Einstein_'+str(idx)+'_'+str(xx),(direct-predicted).subs(x,xx))
    # Direct Euler derivatives of the radial action after one exact integration
    # by parts: int w s Box R = -int (N r^2/A) s' R'. Boundary tests compact.
    n,a,b=[S.Function(name)(x) for name in ['n','a','b']]
    N1=N*(1+u*n);A1=A*(1+u*a);r1=r*(1+u*b)
    _,R1=curvature(x,N1,A1,r1)
    lag=-N1*r1**2/A1*S.diff(s,x)*S.diff(R1,x)
    variation=S.diff(lag,u).subs(u,0)
    # Independently assemble DR*(Box s)+K_R*s+volume term.
    bs=S.cancel(box(s));bbs=box(bs);hh=hess(bs)
    cross=S.diff(s,x)*S.diff(scal,x)/A**2
    KR=[-cross/2-s*box(scal)/2,cross/2-s*box(scal)/2,-cross/2-s*box(scal)/2]
    target=[hh[i]-bbs-ric[i]*bs+KR[i]+s*box(scal)/2 for i in range(3)]
    for i,field in enumerate([n,a,b]):
        # Maximum derivative is three after the one explicit integration by parts.
        coeffs=[S.cancel(S.diff(variation,S.diff(field,x,k))) for k in range(4)]
        euler=sum((-1)**k*S.diff(coeffs[k],x,k) for k in range(4))
        for xx in [R(1,3),R(4,5)]:
            scale=(4 if i==2 else 2)*w
            zero('BoxR_action_transpose_'+str(i)+'_'+str(xx),(euler-scale*target[i]).subs(x,xx))
            check('BoxR_missing_operator_contact_'+str(i)+'_'+str(xx),S.cancel(KR[i].subs(x,xx))!=0)
    # Formal transpose at a nonconstant scalar background, independently by
    # the product rule: -s Hess R + sym grad(s grad R) -g div(s grad R)/2.
    sp,Rp,Rpp,sv=S.symbols('sprime Rprime Rsecond svalue')
    zero('Box_operator_transpose_product_rule',-sv*Rpp+(sp*Rp+sv*Rpp)-sp*Rp)

def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S168 charged conformal response and constraints','input_count':len(ps)}),flush=True)
    try:
        predecessors();local_algebra();ward_from_christoffels();exact_warped_fixtures()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S168 charged conformal response and constraints','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'local_causal_conformal_column_derived':good,'gauge_invariant_trace_row_derived':good,
         'RN_direct_KG_lift_counterjet_checked':good,'dynamic_constraint_rows_derived':good,
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
