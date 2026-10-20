"""S166 forcing prerequisite; does not certify the full fast TT inverse."""
from datetime import datetime, timezone
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
STEM="strong-effective-tt-source"
OWN=["checks/strong_effective_tt_source.py","notes/"+STEM+".md"]
PREDECESSORS=[
 ("stationary-weighted-tt-inverse",52),
 ("stationary-scalar-outgoing-inverse",107),
 ("global-compatible-scalar-response",344),
 ("tt-form-and-scalar-fourth-graph",238),
 ("weighted-york-source-conversion",293),
]
PROTOCOL={"arithmetic":"exact symbolic and rational",
 "coordinate_metric":"diag(a(rho)^2,1,r(rho)^2,r(rho)^2 sin(theta)^2)",
 "projector":"orbit two-plane; full covariant first and second derivatives",
 "g_interval":[".999*2000^-10","2000^-10"],
 "physical_source_fixtures":False,
 "scope":"strong effective fast TT source; full S166 inverse remains open"}
ROWS=[]
R=S.Rational


def check(name,ok,detail=""):
 assert name not in {v["name"] for v in ROWS}
 row={"name":name,"passed":bool(ok),"detail":str(detail)}
 ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)


def zero(name,expr):
 z=S.factor(S.expand(S.simplify(expr)))
 check(name,z==0,"" if z==0 else z)


def latest(stem):
 d=json.loads((ROOT/"evidence"/(stem+".json")).read_text(encoding="utf-8"))
 return d.get("runs",[d])[-1]


def inputs(r,after=False):
 return r.get("input_sha256_after" if after else "input_sha256_before",
              r.get("final_inputs_sha256" if after else "execution_inputs_sha256",{}))


def paths():
 out=set(OWN)
 for stem,_ in PREDECESSORS:
  out.update([f"notes/{stem}.md",f"evidence/{stem}.json"])
  out.update(inputs(latest(stem)))
 return sorted(out)


def predecessors():
 for stem,n in PREDECESSORS:
  r=latest(stem)
  check(stem+"_receipt",r["status"]=="passed" and len(r["checks"])==n
        and all(v["passed"] for v in r["checks"]) and inputs(r)==inputs(r,True))
  for p,h in inputs(r).items():
   check(stem+"_input_"+p,sha256((ROOT/p).read_bytes()).hexdigest()==h)


def coordinate_geometry():
 time,rho,theta,azimuth=S.symbols("time rho theta azimuth",real=True)
 coords=[time,rho,theta,azimuth]
 a,r=[S.Function(n)(rho) for n in ["a","r"]]
 metric=[a*a,S.Integer(1),r*r,r*r*S.sin(theta)**2]
 inv=[1/g for g in metric]
 Gamma={}
 for k,i,j in product(range(4),repeat=3):
  v=((S.diff(metric[j],coords[i]) if k==j else 0)
     +(S.diff(metric[i],coords[j]) if k==i else 0)
     -(S.diff(metric[i],coords[k]) if i==j else 0))*inv[k]/2
  Gamma[k,i,j]=S.simplify(v)
 def covariant(tensor,rank):
  out={}
  for i in range(4):
   for inds in product(range(4),repeat=rank):
    expr=S.diff(tensor.get(inds,S.Integer(0)),coords[i])
    for slot,j in enumerate(inds):
     for k in range(4):
      other=inds[:slot]+(k,)+inds[slot+1:]
      expr-=Gamma[k,i,j]*tensor.get(other,S.Integer(0))
    expr=S.simplify(expr)
    if expr!=0:out[(i,)+inds]=expr
  return out
 def norm2(tensor):
  return sum(S.prod(inv[i] for i in inds)*v*v for inds,v in tensor.items())
 P={(0,0):metric[0],(1,1):metric[1]}
 nP=covariant(P,2);nnP=covariant(nP,3)
 RR=S.diff(r,rho)/r;A=S.diff(a,rho)/a
 zero("projector_first_covariant_norm",norm2(nP)-4*RR**2)
 zero("projector_second_covariant_norm",
      norm2(nnP)-(4*S.diff(RR,rho)**2+4*(A*RR)**2+24*RR**4))
 # Scalar Hessian and product identities in the same coordinate realization.
 v=S.Function("v")(rho)
 H=covariant(covariant({():v},0),1)
 for i in range(4):
  expected=[A*S.diff(v,rho),S.diff(v,rho,2),RR*S.diff(v,rho),RR*S.diff(v,rho)][i]
  zero("scalar_Hessian_frame_"+str(i),H.get((i,i),0)/metric[i]-expected)
 ww=r*r
 Wrr=S.diff(ww,rho,2)/ww
 zero("W_radial_Hessian",Wrr-(2*S.diff(r,rho,2)/r+2*RR**2))
 zero("W_time_Hessian",A*S.diff(ww,rho)/ww-2*A*RR)
 zero("W_sphere_Hessian",RR*S.diff(ww,rho)/ww-2*RR**2)
 zero("W_laplacian",
      (S.diff(ww,rho,2)+(A+2*RR)*S.diff(ww,rho))/ww
      -(2*S.diff(r,rho,2)/r+2*A*RR+6*RR**2))
 f,w,m,j,h=[S.Function(n)(rho) for n in ["f","w","m","j","h"]]
 lap=lambda z:S.diff(m*S.diff(z,rho),rho)/m
 zero("scalar_Laplacian_product",lap(w*f)-w*lap(f)-2*S.diff(w,rho)*S.diff(f,rho)-f*lap(w))
 zero("scalar_Hessian_product",S.diff(w*f,rho,2)-w*S.diff(f,rho,2)
      -2*S.diff(w,rho)*S.diff(f,rho)-f*S.diff(w,rho,2))
 zero("coefficient_Box_transpose_boundary",
      m*(f*lap(j*h)-j*h*lap(f))
      -S.diff(m*(f*S.diff(j*h,rho)-j*h*S.diff(f,rho)),rho))
 K=S.Matrix(2,2,S.symbols("k0:4"));w0,w1=S.symbols("w0 w1",positive=True);W=S.diag(w0,w1)
 check("physical_reciprocal_transpose",(W*K.T*W.inv()).T==W.inv()*K*W)


def budgets():
 alpha=S.Integer(2000)**-5
 e=alpha**2/180;v=R("4e-27")/180
 check("W_gradient",4*5<R("4.5")**2)
 check("W_Hessian_radial",2*11+2*5==32)
 check("W_Laplacian",4*11+6*5==74)
 check("W_Hessian_norm",32**2+22**2+2*10**2<75**2)
 check("weighted_Hessian_gradient_coefficient",4*R("4.5")+4==22)
 check("weighted_Hessian_value_coefficient",74+4*R("4.5")+75==167)
 check("quarter_root_g",alpha<R("6e-9")**2)
 check("weighted_Hessian_final",R("3.2e9")+22*R("4e8")*R("6e-9")
       +167*R("4.4e7")*alpha<R("3.3e9"))
 check("projector_second_norm",4*16**2+4*11**2+24*5**2<300**2)
 check("W0_first_derivative",48**2*5<108**2)
 check("W0_second_derivative",24*(4*5)+24*300==7680)
 check("j_second_coefficient",30*R("3.6")==108)
 check("w_second_coefficient",30*R("5.1")==153)
 check("j_Hessian_budget",5*108+(10+11)*20+11*20+2*5*20<1500)
 check("w_Hessian_budget",5*153+(10+11)*40+11*40+2*5*40<3000)
 check("PminusS_first_derivative",16*5<9**2)
 check("J_second_budget",2*1500+2*30*9+R("3.5")*600<10000)
 check("C_second_budget",7*3000+2*43*108+R("5.1")*7680<100000)
 check("tensor_Laplacian_trace_factor",4==2**2)
 check("trace_adjoint_Hessian",14+R(7,12)*6*40+R(11,6)*14<200)
 check("trace_adjoint_gradient",R(7,12)*12*1000+R(11,6)*2*100<8000)
 alg=R(7*11*13,36)+4+R(7,12)*(6*100000+12320)+R(11,6)*(20000+7*(2*54+2*14)+R(13,2)*7+2*7**2)
 check("trace_adjoint_value",alg<500000,alg)
 check("TT_lift_reciprocal",2*2**2+4*2+6==22<25)
 check("return_reciprocal",9*(2+1)==27)
 # Strict lower bound on g gives the inverse fourth-root reserve.
 check("lower_g_inverse_quarter_root",R(".999")*alpha**2*S.Integer(200000000)**4>1)
 value=(5*((7+20*v+500000*e)*R("4.4e7")
           +R(200,180)*R("3.3e9")*alpha
           +R(8000,180)*R("4e8")*alpha*R("6e-9"))
        +27*R("6e-22")*R("4e8")*200000000)/24
 check("strong_mixed_transpose_final",value<R("6.5e7"),value)
 check("effective_TT_source_final",1+R("6.5e7")<R("6.6e7"))
 check("original_disk_TT_source",R("6.6e7")*124==R("8.184e9"))


def main():
 sys.stdout.reconfigure(encoding="utf-8",newline="\n")
 sys.stderr.reconfigure(encoding="utf-8",newline="\n")
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
 before,error=hashes(),None
 print(json.dumps({"event":"start","story":"S166 forcing prerequisite","input_count":len(ps)}),flush=True)
 try:
  predecessors();coordinate_geometry();budgets()
 except Exception:
  error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c["passed"] for c in ROWS)
 run={"story":"S166 forcing prerequisite","created_utc":datetime.now(timezone.utc).isoformat(),
      "status":"passed" if good else "failed","python":sys.version,"sympy":S.__version__,
      "protocol":PROTOCOL,"input_sha256_before":before,"input_sha256_after":after,
      "unchanged_inputs":before==after,"checks":ROWS,"exception":error,
      "S166_complete":False,"goal_complete":False,
      "scope":"strong effective TT source in fast weighted L2, with actual mixed transpose",
      "unresolved":["global fast fourth/log TT inverse","mass coefficient","full residual composition",
                    "nonlinear family","signed preparation condition","physical requirements"]}
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
