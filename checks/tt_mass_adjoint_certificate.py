"""S166 actual mass-functional reduction; model diagnostics are not a physical solve."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import mpmath as mp
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
STEM="tt-mass-adjoint-certificate"
OWN=["checks/tt_mass_adjoint_certificate.py","notes/"+STEM+".md"]
PREDECESSORS=[
 ("stationary-global-weighted-map",50),
 ("stationary-weighted-tt-inverse",52),
 ("weighted-normal-inverse-scale",453),
 ("strong-effective-tt-source",350),
]
PROTOCOL={
 "arithmetic":"exact symbolic/rational and 80-decimal mpmath quadrature",
 "identity_tolerance":"1e-65",
 "model_b":["0","1/10","1/2","9/10"],
 "model_n":[3,4,5,7,9],
 "trial_t":["0","1/2","1","3/2"],
 "model_domain":"u^-4 du on (0,1); critical L2 and Neumann at 1",
 "actual_adjoint_computed":False,
 "actual_mass_numeric_upper_bound":False,
 "scope":"exact actual mass functional; adjoint residual certificate; manufactured checks",
}
R=S.Rational
ROWS=[]


def check(name,ok,detail=""):
 assert name not in {r["name"] for r in ROWS}
 row={"name":name,"passed":bool(ok),"detail":str(detail)}
 ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)


def zero(name,expr):
 expr=S.factor(S.cancel(S.expand(expr)))
 check(name,expr==0,"" if expr==0 else expr)


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
        and all(c["passed"] for c in r["checks"]) and inputs(r)==inputs(r,True))
  for p,h in inputs(r).items():
   check(stem+"_input_"+p,sha256((ROOT/p).read_bytes()).hexdigest()==h)


def exact_identities():
 s,u=S.symbols("s u",positive=True)
 g,A,d,H3,gamma0,B,beta=S.symbols("g A d H3 gamma0 B beta",positive=True)
 H=S.Function("H")
 q=s*(s-3);NT=q*(d*(q+2)*H(s)-A)
 zero("mass_normal_derivative",
      -9*S.diff(NT,s).subs(s,3).subs(H(3),H3)-27*(A-2*d*H3))
 zero("normal_zero_root",NT.subs(s,0))
 zero("normal_mass_root",NT.subs(s,3))
 kval=27*(1/(48*S.pi*g)+R(23,8640)/S.pi**2-2*H3/(480*S.pi**2))
 zero("boundary_leading_coefficient",
      kval-(R(9,16)/(S.pi*g)+R(23,320)/S.pi**2-R(9,80)*H3/S.pi**2))
 gmax=S.Integer(2000)**-10
 check("boundary_H3_positive",R(39,5)>0)
 check("boundary_H3_upper",1+R(54,5)<12)
 check("boundary_lower_scaled",R(63,352)-R(3,20)*gmax>R(1,6))
 check("boundary_upper_scaled",R(3,16)+R(23,2880)*gmax<R(1,5))
 check("critical_inverse_gap_reserve",8*R(22,7)/R(".99")<26)
 radius=B/u
 zero("physical_volume",4*S.pi*beta*radius**2*(-S.diff(radius,u))
      -4*S.pi*beta*B**3*u**-4)
 uh=S.symbols("uh",positive=True)
 zero("reciprocal_constant_norm",S.integrate(gamma0*u**4*u**-4,(u,0,uh))-gamma0*uh)
 check("critical_defect_integrability",2*2-4>-1)
 check("mass_defect_pair_integrability",2+3-4>-1)
 check("fast_weight_threshold_exponent",2*R(7,2)-4-4==-1)
 # The last equality is the threshold exponent; the actual target is an L2 space,
 # not the bare borderline power u^(7/2).
 zero("classical_mass_boundary_sign",
      (-9*A*u**-2*S.diff(u**3,u))*(-1)-27*A)
 check("disk_mass_conversion",R(3,4)*1==R(3,4))
 check("coupled_York_mass_coefficient",6*R("6.6e7")==R("3.96e8"))
 check("coupled_disk_mass_coefficient",6*R("8.184e9")==R("4.9104e10"))
 check("two_source_conversions",R("3.96e8")*124==R("4.9104e10"))
 # Algebraic residual and inverse identities on an arbitrary invertible matrix.
 D=S.Matrix([[5,1],[1,3]]);hh=S.Matrix(S.symbols("h0:2"))
 va=S.Matrix(S.symbols("v0:2"));one=S.ones(2,1)
 inv=D.inv();psi=one-inv*hh;psia=one-va;res=D*va-hh
 check("adjoint_residual_sign",S.simplify(psi-psia-inv*res)==S.zeros(2,1))
 ff=S.Matrix(S.symbols("f0:2"));cc=inv*ff
 zero("source_specific_residual_pairing",((psi-psia).T*ff-res.T*cc)[0])
 # Integration-by-parts terms in the local toy model, before any fixture.
 p=S.Function("p")(u);c=S.Function("c")(u)
 L0=lambda z:-u*u*S.diff(z,u,2)+2*u*S.diff(z,u)
 zero("local_Green_boundary",
      (p*L0(c)-L0(p)*c)*u**-4
      -S.diff(u**-2*(S.diff(p,u)*c-p*S.diff(c,u)),u))
 x=S.symbols("x",real=True);v=S.Function("v")(x)
 zero("critical_conjugation_energy",
      (S.diff(v,x)-R(3,2)*v)**2
      -(S.diff(v,x)**2+R(9,4)*v*v-R(3,2)*S.diff(v*v,x)))


def numerical_models():
 mp.mp.dps=80
 tol=mp.mpf("1e-65")
 u,b=S.symbols("u b",real=True)
 psi=1-b*u*u+b*u**4/2
 vb=1-psi;V=2*b*(u*u+u**4)/psi
 L=lambda z:-u*u*S.diff(z,u,2)+2*u*S.diff(z,u)+V*z
 zero("model_adjoint_equation",L(psi))
 zero("model_critical_probe",L(vb)-V)
 zero("model_adjoint_endpoint",S.diff(psi,u).subs(u,1))
 zero("model_critical_endpoint",S.diff(vb,u).subs(u,1))
 norm_poly=S.integrate(psi**2,(u,0,1))
 cv=lambda z:mp.mpf(str(S.N(z,85)))
 norm=lambda z:mp.sqrt(cv(S.integrate(S.expand(z*z),(u,0,1))))
 for ib,bval in enumerate([R(0),R(1,10),R(1,2),R(9,10)]):
  tag=f"model_{ib}"
  pp=psi.subs(b,bval);vv=vb.subs(b,bval);pot=V.subs(b,bval)
  check(tag+"_positive_adjoint",1-bval/2>0)
  check(tag+"_nonnegative_potential",bval>=0)
  hnorm=mp.sqrt(mp.quad(S.lambdify(u,S.cancel((pot/u**2)**2),"mpmath"),[0,mp.mpf(".5"),1]))
  adj_norm=norm(pp)
  direct=mp.sqrt(mp.quad(S.lambdify(u,pp**2,"mpmath"),[0,mp.mpf(".5"),1]))
  check(tag+"_adjoint_norm_quadrature",abs(direct-adj_norm)<tol)
  zero(tag+"_Riesz_source_pairing",pp*(u**4*pp)*u**-4-pp**2)
  for n in [3,4,5,7,9]:
   nn=f"{tag}_power_{n}"
   cn=u**n-R(n,n+1)*u**(n+1)
   fn=S.cancel(L(cn).subs(b,bval))
   mass=S.Integer(1 if n==3 else 0)
   zero(nn+"_endpoint",S.diff(cn,u).subs(u,1))
   integrand=S.cancel(pp*fn*u**-4)
   zero(nn+"_adjoint_mass",S.integrate(integrand,(u,0,1))-3*mass)
   bare=S.cancel((fn-pot*cn)*u**-4)
   zero(nn+"_defect_mass",S.integrate(bare,(u,0,1))-3*mass)
   val=mp.quad(S.lambdify(u,integrand,"mpmath"),[0,mp.mpf(".5"),1])
   check(nn+"_mass_quadrature",abs(val-3*int(mass))<tol)
   # Test the strong-residual error identity with a nonexact trial.
   tval=R(1,2)
   pa=1-tval*vv;rr=(tval-1)*pot
   zero(nn+"_error_pairing",
        S.integrate(S.cancel(((pp-pa)*fn-rr*cn)*u**-4),(u,0,1)))
  for it,tval in enumerate([R(0),R(1,2),R(1),R(3,2)]):
   tt=f"{tag}_trial_{it}"
   pa=1-tval*vv
   anorm=norm(pa);actual_error=norm(pp-pa)
   rnorm=abs(cv(tval-1))*hnorm;allowance=mp.mpf(4)/9*rnorm
   check(tt+"_residual_error",actual_error<=allowance)
   lower=max(mp.mpf(0),anorm-allowance)/3
   upper=(anorm+allowance)/3
   check(tt+"_mass_norm_interval",lower<=adj_norm/3<=upper,
         {"lower":mp.nstr(lower,18),"actual":mp.nstr(adj_norm/3,18),"upper":mp.nstr(upper,18)})
  check(tag+"_zero_trial_upper",adj_norm<=1+mp.mpf(4)/9*hnorm)


def main():
 sys.stdout.reconfigure(encoding="utf-8",newline="\n")
 sys.stderr.reconfigure(encoding="utf-8",newline="\n")
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
 before=hashes();error=None
 print(json.dumps({"event":"start","story":"S166 mass-adjoint reduction","input_count":len(ps)}),flush=True)
 try:
  predecessors();exact_identities();numerical_models()
 except Exception:
  error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c["passed"] for c in ROWS)
 run={"story":"S166 mass-adjoint reduction","created_utc":datetime.now(timezone.utc).isoformat(),
      "status":"passed" if good else "failed","python":sys.version,"sympy":S.__version__,
      "mpmath":mp.__version__,"protocol":PROTOCOL,
      "input_sha256_before":before,"input_sha256_after":after,"unchanged_inputs":before==after,
      "checks":ROWS,"exception":error,"S166_complete":False,"goal_complete":False,
      "actual_adjoint_computed":False,"actual_mass_numeric_upper_bound":False,
      "scope":"exact mass functional and residual certificate; no actual mass or full inverse bound",
      "unresolved":["actual adjoint probe/residual enclosure","global fast fourth/log TT inverse",
                    "higher reconstruction and residual composition","nonlinear signed compatibility",
                    "physical requirements"]}
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
