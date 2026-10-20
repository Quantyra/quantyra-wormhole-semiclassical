"""S166 prerequisite: weighted York conversion; full TT inverse is still open."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
STEM="weighted-york-source-conversion"
OWN=["checks/weighted_york_source_conversion.py","notes/"+STEM+".md"]
PREDECESSORS=[
    ("stationary-global-tt-reduction",44),
    ("stationary-weighted-tt-inverse",52),
    ("stationary-scalar-outgoing-inverse",107),
    ("global-compatible-scalar-response",344),
    ("tt-form-and-scalar-fourth-graph",238),
]
PROTOCOL={
    "arithmetic":"exact rational and symbolic",
    "seed_box":{"x_min":1000,"delta":[0,1],"abs_d_max":1},
    "weights":["(r/B)^2","(B/r)^2"],
    "Schur_exponents":["9/10","4/5","1/2"],
    "homogeneous_frequencies":[0,1,10,100],
    "diagnostic_profiles_are_physical_sources":False,
    "scope":"source conversion prerequisite; S166 inverse incomplete",
}
ROWS=[]
R=S.Rational


def check(name,ok,detail=""):
    assert name not in {x["name"] for x in ROWS}
    row={"name":name,"passed":bool(ok),"detail":str(detail)}
    ROWS.append(row)
    print(json.dumps(row),flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name,expr):
    z=S.factor(S.expand(S.simplify(expr)))
    check(name,z==0,"" if z==0 else z)


def latest(stem):
    d=json.loads((ROOT/"evidence"/(stem+".json")).read_text(encoding="utf-8"))
    return d.get("runs",[d])[-1]


def inputs(run,after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256",{}))


def paths():
    names=set(OWN)
    for stem,_ in PREDECESSORS:
        names.update([f"notes/{stem}.md",f"evidence/{stem}.json"])
        names.update(inputs(latest(stem)))
    return sorted(names)


def predecessors():
    for stem,n in PREDECESSORS:
        r=latest(stem)
        check(stem+"_receipt",r["status"]=="passed" and len(r["checks"])==n
              and all(c["passed"] for c in r["checks"]) and inputs(r)==inputs(r,True))
        for p,h in inputs(r).items():
            check(stem+"_input_"+p,sha256((ROOT/p).read_bytes()).hexdigest()==h)


def exact():
    t,x,d,delta=S.symbols("t x d delta",positive=True)
    U=x*(t**4-4*t+3)+(t-1)**2-delta*t+d*(1-2/t+1/t**2)
    first=(t-1)*(x*(t**3+t*t+t+1197)+399-199*t+d*(599-399*t)/t**2)+299*delta*t
    second=x*(4*t**3+1196)-398*t+598+299*delta+d*(1198/t**3-998/t**2)
    zero("seed_first_global_identity",100*t*S.diff(U,t)-399*U-first)
    zero("seed_second_global_identity",100*t*S.diff(U,t,2)-299*S.diff(U,t)-second)
    check("first_polynomial_lower_constants",1000-199>0 and 1000*1197+399-998>0)
    check("second_polynomial_lower_constants",4000-398>0 and 1000*1196+598-2196>0)
    # Explicit shifted polynomials certify the powers on the entire half-line.
    q=S.symbols("q",nonnegative=True)
    for name,p in [("first",1000*t**3-199*t),("second",4000*t**3-398*t)]:
        check(name+"_halfline_positive",all(v>=0 for v in S.Poly(p.subs(t,1+q),q).all_coeffs()))
    check("Schur_growth_exponent",R(301,399)<R(4,5))
    u,v=S.symbols("u v",positive=True)
    for label,p,gamma in [("fast",R(9,10),R(4,5)),("reciprocal",R(1,2),S.Integer(0))]:
        row=u**(p-1)*S.integrate(v**(-p),(v,0,u))
        col=u**(p-gamma)*S.integrate(v**(gamma-p-1),(v,u,S.oo))
        zero(label+"_Schur_row",row-1/(1-p))
        zero(label+"_Schur_column",col-1/(p-gamma))
        check(label+"_Schur_finiteness",0<p<1 and gamma<p)
    check("fast_Schur_constant",1/(1-R(9,10))/(R(9,10)-R(4,5))==100)
    check("reciprocal_Schur_constant",1/(1-R(1,2))/(R(1,2)-0)==4)

    # The total derivative includes the endpoint contribution.
    A,rr,F2,m,p,eta,ep=S.symbols("A R F2 m p eta eta_prime")
    dm=m*(A+2*(1+p)*rr)
    dA=F2/2-A*A
    boundary=dm*A*eta**2+m*dA*eta**2+2*m*A*eta*ep
    lhs=m*(ep-A*eta)**2
    rhs=m*(ep**2+(F2/2+A*A+2*(1+p)*A*rr)*eta**2)
    zero("weighted_return_identity_with_boundary",lhs-rhs+boundary)
    zero("fast_return_potential", (F2/2+A*A+2*(1+p)*A*rr).subs(p,2)
         -(F2/2+A*A+6*A*rr))
    zero("reciprocal_return_potential", (F2/2+A*A+2*(1+p)*A*rr).subs(p,-2)
         -(F2/2+A*A-2*A*rr))
    H,V,r=S.symbols("H V r")
    zero("tail_Hardy_endpoint_identity",H**2+2*r*H*(-V/r)-(H**2-2*H*V))
    K=S.Matrix(2,2,S.symbols("k0:4"))
    w0,w1=S.symbols("w0 w1",positive=True)
    W=S.diag(w0,w1)
    check("physical_transpose_reciprocal_weights",
          (W*K.T*W.inv()).T==W.inv()*K*W)
    M=S.Matrix([[1,0],[R(1,3),2]])
    check("dual_disk_source_matrix_norm",sum(a*a for a in M)==R(46,9)<9)
    B,z,beta=S.symbols("B z beta",positive=True)
    # r=B/z, |dr|=B dz/z^2: coefficient of dz on both sides.
    zero("physical_fast_measure_conversion",(B/z)**2*(B/z/B)**4*B/z**2-B**3*z**-8)


def budgets():
    check("fast_Reta",10<R("3.2")**2)
    check("outer_sqrtF",1+R(1,4000)<R("1.001")**2)
    check("outer_Fr",2+(4+R(3,1000))/8+R(6,32000)+R(4,64000)<3)
    check("outer_eta_Hardy",4*R("1.001")<5)
    check("inner_eta_square",2*5==10<4**2)
    check("inner_Aeta_square",R(13**2,2)<10**2)
    check("outer_input_eta",5**2+4**2<7**2)
    check("outer_input_Aeta",6**2+10**2<12**2)
    check("combined_eta",4**2*2+7**2==9**2)
    check("combined_Aeta",4**2+12**2<13**2)
    check("combined_Reta",5*9**2<21**2)
    check("fast_S",R(1,2)+R("3.2")*11==R("35.7"))
    check("fast_T",R(10,2)+R(3,2)+(1+R("3.2"))*11==R("52.7"))
    check("fast_return_pair",R("35.7")**2+R("52.7")**2<64**2)
    check("fast_full_York_map",1+64**2<65**2)
    check("reciprocal_S",R(1,2)+21*3==R("63.5"))
    check("reciprocal_T",R(2,2)+R(3,2)+(13+21)*3==R("104.5"))
    check("reciprocal_return_pair",R("63.5")**2+R("104.5")**2<123**2)
    check("reciprocal_full_York_map",1+R("63.5")**2+R("104.5")**2<123**2)
    check("physical_source_sum",1+123==124)
    check("dual_disk_source_sum",124*3==372)
    for name,a,b in [("value","4.4e7","5.456e9"),("gradient","4e8","4.96e10"),
                     ("second_graph","3.2e9","3.968e11"),("fourth_graph","6.2e11","7.688e13")]:
        check("source_conversion_"+name,R(a)*124==R(b))
    check("unweighted_return",600<25**2)
    check("full_shape_L2",R("5.1")*3500==17850)
    check("full_metric_L2",R("4.4e7")+25*17850<R("4.5e7"))
    check("full_metric_source_conversion",R("4.5e7")*124==R("5.58e9"))


def homogeneous_diagnostic():
    r=S.symbols("r",positive=True)
    J=S.Matrix([[0,1],[-1,0]])
    for k in PROTOCOL["homogeneous_frequencies"]:
        M=4*(S.eye(2)/2+k*J).inv()
        zero(f"homogeneous_{k}_matrix_norm",(M.T*M)[0,0]-16/(R(1,4)+k*k))
        check(f"homogeneous_{k}_isotropic",M.T*M==16/(R(1,4)+k*k)*S.eye(2))
        check(f"homogeneous_{k}_bound",16/(R(1,4)+k*k)<=64)
        basis=S.Matrix([[S.cos(k*S.log(r)),S.sin(k*S.log(r))]])
        for col in range(2):
            c=r**-R(7,2)*basis[col]
            out=r**-R(7,2)*(basis*M)[col]
            zero(f"homogeneous_{k}_{col}_integral_derivative",
                 S.diff(r**4*out,r)-4*r**3*c)
        if k==0:
            check("homogeneous_zero_sharp",M==8*S.eye(2))


def main():
    sys.stdout.reconfigure(encoding="utf-8",newline="\n")
    sys.stderr.reconfigure(encoding="utf-8",newline="\n")
    names=paths()
    hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in names}
    before,error=hashes(),None
    print(json.dumps({"event":"start","story":"S166.1","input_count":len(names)}),flush=True)
    try:
        predecessors(); exact(); budgets(); homogeneous_diagnostic()
    except Exception:
        error=traceback.format_exc()
        print(error,file=sys.stderr,flush=True)
    after=hashes()
    good=error is None and before==after and all(c["passed"] for c in ROWS)
    run={"story":"S166.1","created_utc":datetime.now(timezone.utc).isoformat(),
         "status":"passed" if good else "failed","python":sys.version,"sympy":S.__version__,
         "protocol":PROTOCOL,"input_sha256_before":before,"input_sha256_after":after,
         "unchanged_inputs":before==after,"checks":ROWS,"exception":error,
         "scope":"weighted York return and numerical conversion to original physical disk source",
         "S166_complete":False,"goal_complete":False,
         "unresolved":["fast fourth/log TT inverse","mass coefficient","matched residual/inverse composition",
                       "nonlinear family","signed preparation source","physical completion"]}
    if not good:
        run["failed_checker_source_utf8"]=(ROOT/OWN[0]).read_text(encoding="utf-8")
        run["failed_note_source_utf8"]=(ROOT/OWN[1]).read_text(encoding="utf-8")
    dest=ROOT/"evidence"/(STEM+".json")
    data=json.loads(dest.read_text(encoding="utf-8")) if dest.exists() else {"runs":[]}
    data["runs"].append(run)
    dest.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps({"status":run["status"],"checks":len(ROWS),"unchanged_inputs":before==after}),flush=True)
    return 0 if good else 1


if __name__=="__main__":
    raise SystemExit(main())
