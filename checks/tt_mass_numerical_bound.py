"""S166 actual mass upper bound; no full fast-remainder inverse or physical solve."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
STEM="tt-mass-numerical-bound"
OWN=["checks/tt_mass_numerical_bound.py","notes/"+STEM+".md"]
PREDECESSORS=[
 ("stationary-global-tt-reduction",44),
 ("stationary-curved-quantum-tt",157),
 ("tt-form-and-scalar-fourth-graph",238),
 ("strong-effective-tt-source",350),
 ("tt-mass-adjoint-certificate",388),
]
PROTOCOL={
 "arithmetic":"exact symbolic and rational; no fitted constants",
 "probe":"York projection of fixed disk (S,T)=(3/2,-3/2)",
 "mass_coefficient":"scalar sphere c at u^3; disk cubic is 3/4 of this",
 "actual_adjoint_solve":False,
 "quantum_response":"full previously bounded constant-contact row; no angular truncation",
 "independent_local_check":"mixed differentiation of complete integrated Weyl density",
 "scope":"numerical actual TT mass upper bound; global fast fourth/log remainder remains open",
}
R=S.Rational
ROWS=[]


def check(name,ok,detail=""):
 assert name not in {v["name"] for v in ROWS}
 row={"name":name,"passed":bool(ok),"detail":str(detail)}
 ROWS.append(row);print(json.dumps(row),flush=True)
 if not ok:raise AssertionError(name)


def zero(name,expr):
 expr=S.factor(S.cancel(S.expand(S.simplify(expr))))
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
  rr=latest(stem)
  check(stem+"_receipt",rr["status"]=="passed" and len(rr["checks"])==n
        and all(v["passed"] for v in rr["checks"]) and inputs(rr)==inputs(rr,True))
  for p,h in inputs(rr).items():
   check(stem+"_input_"+p,sha256((ROOT/p).read_bytes()).hexdigest()==h)


def projection_and_potential():
 r,rh,lam,B,x,d,qc=S.symbols("r rh lambda B x d qc",positive=True)
 F=S.Function("F")(r)
 A=S.diff(F,r)/(2*S.sqrt(F));RR=S.sqrt(F)/r;w=A+2*RR
 drho=lambda z:S.sqrt(F)*S.diff(z,r)
 div=lambda z:drho(z)+w*z
 pot=drho(w)/2+2*A*A+4*RR*RR
 P=lambda z:R(3,2)*(drho(drho(z))+w*drho(z))-pot*z
 Y0=(r**3-rh**3)/(r*r*S.sqrt(F))
 ss=S.diff(F,r,2)+2*S.diff(F,r)/r
 zero("comparison_vector_divergence",div(Y0)-3)
 zero("comparison_projection_residual",P(Y0)+ss*Y0)
 zero("comparison_sphere_component",
      R(3,2)+drho(Y0)/2-RR*Y0+A*Y0/2-(1+2*(rh/r)**3))
 zero("projection_energy_potential",
      pot-(S.diff(F,r,2)/4+S.diff(F,r)/(2*r)+3*RR**2+R(3,2)*A**2))
 Z=S.Function("Z")(r);ric=-ss/2
 zero("radial_Ricci_divergence_identity",P(Z)-R(3,2)*drho(div(Z))-2*ric*Z)
 rho,fh=S.symbols("rho fh",positive=True)
 # r-rh=fh*rho^2/4 at leading order; the coefficient is independent of fh.
 yy=(3*fh*rho**2/4)/S.sqrt(fh**2*rho**2/4)
 zero("regular_bolt_vector_leading",yy-R(3,2)*rho)
 X=S.symbols("X",positive=True)
 V0=2+6*lam*X-2*qc/X+2*d*(-B**2/X+4*B**3/X**R(3,2)-3*B**4/X**2)
 pp=3*(X*S.diff(V0,X)-V0)
 p_exp=-6+12*(qc+d*B**2)/X-60*d*B**3/X**R(3,2)+54*d*B**4/X**2
 zero("I0_scaling_potential",pp-p_exp)
 zero("I0_S_source",-pp-X*S.diff(pp,X)-(6-30*d*B**3/X**R(3,2)+54*d*B**4/X**2))
 # Independently polarize the reference density, keeping the second X variation.
 alpha,epsilon,sv,tv=S.symbols("alpha epsilon sv tv")
 deltaV=S.Function("DeltaV")
 ss0=R(3,2)*alpha+epsilon*sv
 tt0=-R(3,2)*alpha+epsilon*tv
 density=S.exp(2*(ss0+tt0))*deltaV(X*S.exp(2*ss0))
 mixed=S.diff(density,alpha,epsilon).subs({alpha:0,epsilon:0})
 expected=6*X*S.diff(deltaV(X),X)*(2*sv+tv)+6*X**2*S.diff(deltaV(X),X,2)*sv
 zero("reference_constant_direction_chain_rule",mixed-expected)
 pv=-R(11,960)/S.pi**2
 zero("reference_vacuum_S",-12*pv-R(11,80)/S.pi**2)
 zero("reference_vacuum_T",-6*pv-R(11,160)/S.pi**2)


def polynomial_and_weyl():
 u,nu,a,b,e,v,Q=S.symbols("u nu a b e v Q",positive=True)
 f=1+nu*u**2-2*a*u**3+b*u**4-2*e*u**5+e*u**6
 rb=-S.diff(f,u,2)
 ub=S.diff(-f/u,u);vb=f/u**2
 ss=6+2*b*u**4-12*e*u**5+12*e*u**6
 FF=f/u**2
 # r partial_r=-u partial_u; S0/lambda=r^-2[(rD)^2+rD]F/lambda.
 ed=lambda z:-u*S.diff(z,u)
 zero("seed_S0_polynomial",u**2*(ed(ed(FF))+ed(FF))-ss)
 NN=S.expand(6*f-ss*(1-v**3*u**3))
 pars=[nu,a,b,e,v,Q]
 guards={nu:R(".001"),a:R("2.002"),b:R("3.002"),e:R("1e-15"),v:R("1.001"),Q:S.Integer(2440)}
 def envelope(expr,extra=()):
  poly=S.Poly(S.expand(expr),u,*extra,*pars)
  return sum(abs(c)*S.prod(guards[p]**power for p,power in zip(pars,mon[-len(pars):]))
             for mon,c in poly.terms())
 check("seed_f_coefficient_sum",envelope(f)<9)
 check("seed_f_minus_one_sum",envelope(f-1)<R("7.01"))
 check("seed_Rbar_sum",envelope(rb)<61)
 check("seed_Ubar_sum",envelope(u**2*ub)<19)
 check("seed_Vbar_sum",envelope(u**2*vb)<9)
 check("projection_residual_coefficient_sum",envelope(NN)<62,envelope(NN))
 check("projection_residual_degree",all(2<=mon[0]<=9 for mon,_ in S.Poly(NN,u,*pars).terms()))
 check("outer_f_lower",1-2*guards[a]/8-3*guards[e]>R(".49"))
 check("seed_S0_bound",6+2*guards[b]+24*guards[e]<13)
 check("seed_S0_derivative",(8*guards[b]+132*guards[e])**2*5<55**2)
 # Independent mixed derivative of the complete local Weyl density after
 # integrating sigma Box R by parts, with compact variable test directions.
 aa,bb=S.symbols("aa bb")
 om,rr,uu,vv,yy,kk=S.symbols("om rr uu vv yy kk")
 sx,tx,ls,lt,ww=S.symbols("sx tx ls lt ww")
 sf=R(3,2)*aa+bb*sx;tf=-R(3,2)*aa+bb*tx
 sig=om+sf;rv=S.exp(-2*tf)*(rr-2*bb*lt)
 uv=S.exp(-2*tf)*(uu+bb*ls);v_v=S.exp(-2*tf)*(vv+2*bb*ww)
 lag=sig*(rv**2/60-yy*rv/18+kk)+rv*(uv+v_v)/30-yy*v_v/18-uv**2/10-R(7,90)*uv*v_v-R(7,180)*v_v**2
 mixed=S.expand(S.diff(S.exp(2*tf)*lag,aa,bb).subs({aa:0,bb:0}))
 om_u=-S.log(u);kb=nu**2/15-Q**2/3
 subs={om:om_u,rr:rb,uu:ub,vv:vb,yy:nu,kk:kb}
 coeff={z:S.diff(mixed,z).subs(subs) for z in [sx,tx,ls,lt,ww]}
 lap=lambda z:S.diff(f*S.diff(z,u),u)
 independent_S=coeff[sx]+lap(coeff[ls])-S.diff(coeff[ww]*f*(-1/u),u)
 independent_T=coeff[tx]+lap(coeff[lt])
 DP=R(3,2)*(-rb/15+nu/9)-R(3,2)*((6*om_u+1)*rb/45+11*vb/135)
 LP=-R(3,2)*(-rb/15+2*ub/5+7*vb/45)
 WB=-2*rb/15+14*(ub+vb)/45
 ZT=om_u*rb**2/30+rb*(ub+vb)/15-ub**2/5-7*ub*vb/45-7*vb**2/90+2*kb*om_u
 ZST=-rb**2/30+2*kb
 divW=-S.diff(f*WB/u,u)
 qS=lap(LP)+R(3,2)*divW-R(3,2)*ZST
 qT=lap(DP)+lap(LP)/3-3*ZT+R(3,2)*ZST
 zero("Weyl_independent_Euler_S",independent_S-qS)
 zero("Weyl_independent_Euler_T",independent_T-qT)
 logu=S.symbols("logu")
 residues=[]
 for name,expr,vac in [("S",-u**4*qS,R(11,5)),("T",-u**4*qT,R(11,10))]:
  poly=S.expand(expr).subs(S.log(u),logu)
  zero("Weyl_vacuum_"+name,poly.subs({u:0,logu:0})-vac)
  rem=S.expand(poly-vac);pol=S.Poly(rem,u,logu,*pars)
  check("Weyl_remainder_degree_"+name,all(mon[0]>=2 and mon[1]<=1 for mon,_ in pol.terms()))
  residues.append(envelope(rem,(logu,)))
 check("Weyl_full_remainder_envelope",sum(residues)<S.Integer(10)**8,sum(residues))
 zero("normal_physical_volume_ratio",(u**2)/(4*S.pi*(1/u)**2)-u**4/(4*S.pi))
 check("Weyl_L2_reserve",R(3,2)*10**8/(16*9)<10**8)


def budgets():
 gmax=S.Integer(2000)**-10
 check("Y0_value",18<R("4.3")**2)
 check("Y0_A_bound",R(3,2)*42==63)
 check("Y0_radial_derivative",3+63<=66)
 check("J0_R_value",36*5<R("13.5")**2)
 check("J0_global_value",(R("13.5")+13*R("4.3"))<70)
 check("J0_global_radial",6*16+55*R("4.3")+13*66<1200)
 check("J0_global_time",6*11+13*63==885)
 check("J0_global_sphere",5*70**2<157**2)
 check("J0_global_covariant",1200+885+2*157<2500)
 check("J0_outer_value",62/R(".7")<89)
 check("J0_outer_radial",9*62+62*43/(2*R(".49"))<4000)
 check("J0_outer_time",11*62/R(".49")<1500)
 check("J0_outer_sphere",5*89**2<200**2)
 check("J0_outer_covariant",4000+1500+2*200<10000)
 check("J0_inner_decay",4*70<300 and 4*2500==10000)
 check("Z_value",R(100,149)*300<202)
 check("Z_radial_and_time",R(200,447)*300**2<201**2)
 check("Z_sphere",R(100,447)*300**2<143**2)
 check("Z_covariant",R(600,447)*300**2<350**2)
 check("tau_value",R(201+201+2*143,4)<180)
 check("probe_cubic_norm",4*R("1.001")**3/3<R("1.2")**2)
 check("probe_weighted_value",1+R("1.2")+R(201,2)+143+R(201,2)<400)
 check("tau_gradient",R(300+28*202,6)<1000)
 check("Ricci_derivative",100+R(1,2)<200)
 check("tau_Hessian",R(10000,6)+R(200*202+14*350,3)<20000)
 # Exact endpoint integrals; k=2 gives the largest allowed norm.
 z=S.symbols("z",positive=True)
 zero("critical_u2_norm",S.integrate(z**4*z**-4,(z,0,1))-1)
 zero("critical_u2_log_norm",S.integrate(S.log(z)**2,(z,0,1))-2)
 check("sqrt_two_reserve",2<R(3,2)**2)
 check("classical_coefficient_sum",12*R(".001")+12*R("3.002")+198*R("1e-15")<37)
 check("classical_L2_reserve",R(37,24)<2)
 check("angular_ratio",2440**2<6000000)
 check("massive_contact_point",R(3,2)*10**20*2440**2*R(3,2)/12<10**27)
 check("reference_pressure",3301+13*3==3340)
 check("reference_derivative",4*3301+70*3<14000)
 check("reference_row",12*3340+3*14000+6*3340<110000)
 check("reference_L2_reserve",110000*2440**2<10**12)
 check("complete_disk_source",2+gmax*(10**27+10**12+10**8)<3)
 check("unweighted_return_transpose",100+500<25**2)
 check("unweighted_TT_lift",22<5**2)
 check("unweighted_vector_return",18<5**2)
 ee=gmax/180;vmag=R("4e-27")/180
 mixed=R(1,24)*(5*((7+20*vmag+500000*ee)*180+200*ee*20000+8000*ee*1000)
                  +5*R("6e-22")*1000)
 check("actual_mixed_transpose_budget",mixed<300,mixed)
 check("actual_TT_probe",25*3+300<400)
 check("actual_adjoint_norm",400+26*400<11000)
 check("actual_mass_norm",6*11000==66000)
 check("coupled_York_mass",66000*R("6.6e7")==R("4.356e12"))
 check("coupled_disk_mass",66000*R("8.184e9")==R("5.40144e14"))
 check("coupled_source_consistency",R("4.356e12")*124==R("5.40144e14"))


def main():
 sys.stdout.reconfigure(encoding="utf-8",newline="\n")
 sys.stderr.reconfigure(encoding="utf-8",newline="\n")
 ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
 before=hashes();error=None
 print(json.dumps({"event":"start","story":"S166 numerical TT mass bound","input_count":len(ps)}),flush=True)
 try:
  predecessors();projection_and_potential();polynomial_and_weyl();budgets()
 except Exception:
  error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
 after=hashes();good=error is None and before==after and all(c["passed"] for c in ROWS)
 run={"story":"S166 numerical TT mass bound","created_utc":datetime.now(timezone.utc).isoformat(),
      "status":"passed" if good else "failed","python":sys.version,"sympy":S.__version__,
      "protocol":PROTOCOL,"input_sha256_before":before,"input_sha256_after":after,
      "unchanged_inputs":before==after,"checks":ROWS,"exception":error,
      "S166_complete":False,"goal_complete":False,"actual_mass_upper_bound":good,
      "actual_adjoint_solve":False,
      "scope":"numerical actual TT mass functional bound through constant canonical probe",
      "unresolved":["global fast fourth/log TT remainder","higher reconstruction and residual composition",
                    "nonlinear signed compatibility","physical requirements"]}
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
