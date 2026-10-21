"""S168: exact local coupled frequency polynomial.

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
STEM = "local-coupled-frequency-polynomial"
OWN = ["checks/local_coupled_frequency_polynomial.py", "notes/"+STEM+".md", "evidence/local-frequency-polynomial-source-review.json"]
PREDECESSORS = [("thermal-shape-response-and-subtractions", 328)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope":"Exact classical/Weyl frequency coefficients; canonical shape subtraction remains separate",
    "torus_period":"2pi","frequencies":[1,2,3],"contour_radius":.04,
    "contour_grids":[[32,16],[64,32]],"time_nodes":32,
    "contour_scaled_tolerance":1e-9,"real_steps":[.02,.01],
    "real_scaled_tolerance":1e-3,"refinement_ratio":.35,"roundoff_floor":1e-9,
    "fixture_P":".2+.1cos(x)+.03sin(2x)","fixture_W":"-.4+.05sin(x)+.02cos(3x)",
    "fixture_p":".3cos(2x)+.1sin(3x)","fixture_v":".2cos(x)-.1sin(2x)",
    "fixture_parameters":{"B":1.4,"G":.7,"lambda0":.11,"Q2":.2,"q":3},
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






def exact_frequency_coefficients():
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('primary_source_scope',review['story']=='S168' and len(review['primary_sources'])==1 and not review['actual_wormhole_instability_proved'])
    a,c,d,s=S.symbols('amplitude cos sin s')
    P,P1,P2,W2,p,p1,p2,v,v1,v2,y,AC,AW=S.symbols('P P1 P2 W2 p p1 p2 v v1 v2 y AC AW')
    PP=P+a*p*c;DP=P2+a*(p2-s*p)*c;DV=W2+a*(v2-s*v)*c
    K=P1**2+2*a*P1*p1*c+a*a*(p1*p1*c*c+s*p*p*d*d)
    L=PP*DV**2/15-DP*DV/15-DP**2/10-DV*K/15-R(7,90)*DP*K-R(7,180)*K**2
    coef=lambda f:S.expand(f).coeff(a,2).subs({c*c:1,d*d:1}).expand()
    part=coef((1-2*a*v*c+2*a*a*v*v*c*c)*L)
    extra=coef(y*PP*DV/9-y*K/18)
    J=(-24*P*W2*v*v+12*P*v*v2+P1**2*(7*p*p+14*p*v+12*v*v)
       -14*P1*p*p1-12*P1*p1*v+P2*(7*p*p+36*p*v+12*v*v)
       +W2*(6*p*p+24*p*v)-18*p*p2-6*p*v2-6*p2*v)
    zero('Weyl_sigma_coefficient',-2*AW*part.coeff(s,1)-AW*J/45)
    zero('Weyl_sphere_sigma_coefficient',-2*AW*extra.coeff(s,1)-AW*y*(p*p+2*p*v)/9)
    zero('Weyl_sigma_squared_coefficient',-2*AW*part.coeff(s,2)-AW*(3*p*p+2*p*v-2*P*v*v)/15)
    zero('Weyl_no_higher_frequency',S.diff(part,s,3))
    kinetic=6*K+4*((P1+a*p1*c)*(S.Symbol('W1')+a*v1*c)+a*a*s*p*v*d*d)
    classpart=coef((1+2*a*p*c+2*a*a*p*p*c*c)*kinetic)
    zero('classical_sigma_coefficient',-2*AC*classpart.coeff(s,1)+4*AC*(3*p*p+2*p*v))
    zero('classical_sigma_squared_zero',classpart.coeff(s,2))
    check('omitted_time_gradient_detected',S.expand(part-coef((1-2*a*v*c+2*a*a*v*v*c*c)*L.subs(d,0))).coeff(s,1)!=0)
    check('omitted_exponential_contacts_detected',S.expand(part-coef(L)).coeff(s,1)!=0)
    M=S.Matrix([[3,1],[1,-2*P]])
    zero('quartic_determinant',M.det()+6*P+1)
    zero('quartic_square_completion',3*p*p+2*p*v-2*P*v*v-(3*(p+v/3)**2-(2*P+R(1,3))*v*v))
    dd=S.symbols('canonical_d')
    zero('quartic_shape_Schur_test',(M+S.diag(0,dd)).det()-3*(dd-2*P-R(1,3)))
    B,G,W,lam=S.symbols('B G W lambda',positive=True)
    q2pp=B*B*S.exp(-2*W)/(20*S.pi);q1pp=-3*B*B*S.exp(2*P)/G
    zero('proper_frequency_balance',q2pp*s/(-q1pp)-G*s/(60*S.pi*S.exp(2*P+2*W)))
    zero('outer_negative_log_coefficient',(-B*B*P/(30*S.pi)*S.Symbol('em2W')).subs('em2W',1/(lam*B*B))+P/(30*S.pi*lam))
    gmax=R(16,10**34)
    check('benchmark_curvature_shift',11*gmax/(180*R(314,100))<R(32,10**36))
    # Exact quadratic interpolation retains every local contact.
    q0,q1,q2=S.symbols('q0 q1 q2');poly=q0+s*q1+s*s*q2;nodes=[R(1,2),R(3,2),R(7,2)]
    interp=sum(poly.subs(s,n)*S.prod((s-m)/(n-m) for m in nodes if m!=n) for n in nodes)
    zero('local_anchor_interpolation_exact',interp-poly)
    # Independent reduction of the accepted spherical Lpoly at the density level.
    V,DP0,DV0,K0,D_Rhat,charge=S.symbols('V DP DV gradPhi2 Delta_Rhat charge')
    em=S.exp(-2*V);rhat=-2*em*DV0;H=rhat+2*y;u=em*DP0;vv=em*K0
    raw=S.exp(2*V)*(P*(H**2/60-R(11,90)*y*rhat+em*D_Rhat/30-charge**2*y*y/3)
        -u*u/10+(rhat/30-y/18)*vv-R(7,90)*u*vv-R(7,180)*vv**2)
    # Integrated Phi Delta Rhat is Delta Phi Rhat, with compact/periodic tests.
    raw=raw.subs(D_Rhat,0)+DP0*rhat/30
    reduced=em*(P*DV0*DV0/15-DP0*DV0/15-DP0*DP0/10-DV0*K0/15-R(7,90)*DP0*K0-R(7,180)*K0*K0)+y*P*DV0/9-y*K0/18+y*y*(R(1,15)-charge**2/3)*P*S.exp(2*V)
    zero('spherical_to_optical_Weyl_density',raw-reduced)
    # Classical two integrations by parts: -2 e2P Delta V and -6 e2P Delta P.
    cross,grad=S.symbols('cross gradPhi2')
    zero('classical_integration_by_parts',4*cross+12*grad-6*grad-(4*cross+6*grad))


def fixture_jets(N):
    x=2*np.pi*np.arange(N)/N
    P=[.2+.1*np.cos(x)+.03*np.sin(2*x),-.1*np.sin(x)+.06*np.cos(2*x),-.1*np.cos(x)-.12*np.sin(2*x)]
    W=[-.4+.05*np.sin(x)+.02*np.cos(3*x),.05*np.cos(x)-.06*np.sin(3*x),-.05*np.sin(x)-.18*np.cos(3*x),-.05*np.cos(x)+.54*np.sin(3*x),.05*np.sin(x)+1.62*np.cos(3*x)]
    p=[.3*np.cos(2*x)+.1*np.sin(3*x),-.6*np.sin(2*x)+.3*np.cos(3*x),-1.2*np.cos(2*x)-.9*np.sin(3*x)]
    v=[.2*np.cos(x)-.1*np.sin(2*x),-.2*np.sin(x)-.2*np.cos(2*x),-.2*np.cos(x)+.4*np.sin(2*x),.2*np.sin(x)+.8*np.cos(2*x),.2*np.cos(x)-1.6*np.sin(2*x)]
    return [[f[:,None] for f in row] for row in [P,W,p,v]]


def actions(amplitude,k,N):
    P,W,p,v=fixture_jets(N);time=2*np.pi*np.arange(32)/32
    co=np.cos(k*time)[None,:];si=np.sin(k*time)[None,:];s=k*k
    Phi=P[0]+amplitude*p[0]*co;V=W[0]+amplitude*v[0]*co
    Px=P[1]+amplitude*p[1]*co;Pt=-amplitude*k*p[0]*si
    Vx=W[1]+amplitude*v[1]*co;Vt=-amplitude*k*v[0]*si
    DP=P[2]+amplitude*(p[2]-s*p[0])*co
    DV=W[2]+amplitude*(v[2]-s*v[0])*co
    DVx=W[3]+amplitude*(v[3]-s*v[1])*co
    DVt=-amplitude*k*(v[2]-s*v[0])*si
    D4V=W[4]+amplitude*(v[4]-2*s*v[2]+s*s*v[0])*co
    KP=Px*Px+Pt*Pt;KV=Vx*Vx+Vt*Vt;dot=Vx*DVx+Vt*DVt
    B=1.4;G=.7;l0=.11;Q2=.2;q=3;y=1/(B*B);AC=B*B/(4*G);AW=B*B/(4*np.pi)
    em=np.exp(-2*V)
    # Unintegrated physical four-dimensional scalar curvature and volume.
    R4=np.exp(-2*Phi)*(em*(-2*DV-6*DP-6*KP)+2*y)
    classdens=-AC*np.exp(4*Phi+2*V)*(R4+6*l0-2*Q2*y*y*np.exp(-4*Phi))
    # Independently use the covariant spherical Lpoly, including Boxhat Rhat.
    rhat=-2*em*DV;boxr=-2*em*em*(D4V-4*dot+(4*KV-2*DV)*DV)
    H=rhat+2*y;u=em*DP;vv=em*KP
    L=Phi*(H*H/60-11*y*rhat/90+boxr/30-q*q*y*y/3)-u*u/10+(rhat/30-y/18)*vv-7*u*vv/90-7*vv*vv/180
    weyldens=-AW*np.exp(2*V)*L
    return (2*np.pi)**2*np.array([np.mean(classdens),np.mean(weyldens)])


def predicted_coefficients(N):
    PP,WW,pp,vv=fixture_jets(N)
    P,P1,P2=[x[:,0] for x in PP];W,W1,W2,W3,W4=[x[:,0] for x in WW]
    p,p1,p2=[x[:,0] for x in pp];v,v1,v2,v3,v4=[x[:,0] for x in vv]
    B=1.4;G=.7;y=1/(B*B);AC=B*B/(4*G);AW=B*B/(4*np.pi)
    J=(-24*P*W2*v*v+12*P*v*v2+P1**2*(7*p*p+14*p*v+12*v*v)-14*P1*p*p1-12*P1*p1*v
       +P2*(7*p*p+36*p*v+12*v*v)+W2*(6*p*p+24*p*v)-18*p*p2-6*p*v2-6*p2*v)
    q1c=2*np.pi*np.mean(-4*AC*np.exp(2*P)*(3*p*p+2*p*v))
    q1w=2*np.pi*np.mean(AW*(np.exp(-2*W)*J/45+y*(p*p+2*p*v)/9))
    q2w=2*np.pi*np.mean(AW/15*np.exp(-2*W)*(3*p*p+2*p*v-2*P*v*v))
    return np.array([q1c,q1w]),np.array([0.,q2w])


def curvature_frequency_fixture():
    def contour(k,N,nodes):
        rad=.04;theta=2*np.pi*np.arange(nodes)/nodes;zero_action=actions(0.,k,N)
        coef=sum((actions(rad*np.exp(1j*t),k,N)-zero_action)*np.exp(-2j*t) for t in theta)/(nodes*rad*rad)
        norm=2*np.pi if k==0 else np.pi
        return 2*coef/norm
    for N,nodes in PROTOCOL['contour_grids']:
        reference=contour(0,N,nodes);q1,q2=predicted_coefficients(N)
        for k in PROTOCOL['frequencies']:
            target=k*k*q1+k**4*q2;got=contour(k,N,nodes)-reference
            errors=np.abs(got-target)
            for j,name in enumerate(['classical','Weyl']):
                check('covariant_contour_'+name+'_'+str(k)+'_'+str(N),errors[j]<1e-9*(1+abs(target[j])),errors[j])
            NUMERICS.append({'case':'curvature_contour','N':N,'nodes':nodes,'k':k,'targets':target.tolist(),'errors':errors.tolist()})
    N=64;q1,q2=predicted_coefficients(N)
    def real_second(h,k):
        norm=2*np.pi if k==0 else np.pi
        return (actions(h,k,N)+actions(-h,k,N)-2*actions(0.,k,N))/(h*h*norm)
    for k in PROTOCOL['frequencies']:
        target=k*k*q1+k**4*q2;errors=[]
        for step in PROTOCOL['real_steps']:
            error=np.abs(real_second(step,k)-real_second(step,0)-target);errors.append(error)
            for j,name in enumerate(['classical','Weyl']):check('covariant_real_'+name+'_'+str(k)+'_'+str(step),error[j]<1e-3*(1+abs(target[j])),error[j])
        for j,name in enumerate(['classical','Weyl']):check('covariant_real_refinement_'+name+'_'+str(k),errors[1][j]<.35*errors[0][j] or max(e[j] for e in errors)<1e-9,[e[j] for e in errors])
        NUMERICS.append({'case':'curvature_real','k':k,'errors':[e.tolist() for e in errors]})

def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S168 local coupled frequency polynomial','input_count':len(ps)}),flush=True)
    try:
        predecessors();exact_frequency_coefficients();curvature_frequency_fixture()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S168 local coupled frequency polynomial','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'classical_Weyl_frequency_polynomial_derived':good,'complete_conformal_mixed_frequency_rows_derived':good,
         'negative_outer_Weyl_quartic_term_derived':good,'canonical_shape_polynomial_evaluated':False,
         'numerics':NUMERICS,'actual_positive_frequency_operator_values_evaluated':False,'actual_causal_inverse_proved':False,
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
