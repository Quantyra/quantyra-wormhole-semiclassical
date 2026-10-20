"""S152: global fixed-chart source-to-constraint transfer.

Manufactured profiles are diagnostics, not physical quantum states.
"""
from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import importlib.util
import json
import sys
import traceback

import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
STEM = "fixed-chart-global-constraint-transfer"
INPUTS = [
    "checks/fixed_chart_global_constraint_transfer.py",
    "notes/fixed-chart-global-constraint-transfer.md",
    "notes/stationary-centered-complete-residual.md",
    "notes/stationary-centered-covariant-seed.md",
    "notes/stationary-centered-constraint-residual.md",
    "notes/interior-pressure-and-bolt-primitive.md",
    "checks/interior_pressure_and_bolt_primitive.py",
    "notes/common-chart-coordinate-transfer.md",
    "evidence/common-chart-coordinate-transfer.json",
]
ROWS, NUMERICS = [], []
PROTOCOL = {
    "digits": 55, "nodes": [32, 64], "steps": [".001", ".0005"],
    "gaps": ["1e-5", "1e-30"], "endpoints": [".75", "2.5"],
    "mixed_tolerance": "2e-6", "identity_tolerance": "1e-11",
    "refinement_ratio": ".35", "refinement_floor": "1e-30",
}


def check(name, passed, detail=""):
    assert not any(r["name"] == name for r in ROWS)
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not passed:
        raise AssertionError(name)


def predecessor_checks():
    last = json.loads((ROOT/INPUTS[-1]).read_text(encoding="utf-8"))["runs"][-1]
    check("s151_final_receipt", last["status"] == "passed"
          and len(last["checks"]) == 111 and all(r["passed"] for r in last["checks"]))
    check("s151_current_inputs", last["unchanged_inputs"] and all(
        sha256((ROOT/p).read_bytes()).hexdigest() == h
        for p, h in last["input_sha256_after"].items()))


def exact_checks():
    b,V,z,tau = s.symbols("b V z tau", real=True)
    hb,hV,hz,kb,kV,kz = s.symbols("hb hV hz kb kV kz", real=True)
    ss,tt = s.symbols("ss tt", real=True)
    T = 1+tau*(s.exp(V)-1)
    a = tau*s.exp(V)/T
    JM = s.exp(V)*(1+z)/T**2
    JA = s.exp(2*b-V)/(T*(1+z))
    def dh(expr):
        return hb*s.diff(expr,b)+hV*s.diff(expr,V)+hz*s.diff(expr,z)
    def dk(expr):
        return kb*s.diff(expr,b)+kV*s.diff(expr,V)+kz*s.diff(expr,z)
    lm = (1-2*a)*hV+hz/(1+z)
    la = 2*hb-(1+a)*hV-hz/(1+z)
    lmk = -2*a*(1-a)*hV*kV-hz*kz/(1+z)**2
    lak = -a*(1-a)*hV*kV+hz*kz/(1+z)**2
    check("mass_multiplier_first", s.simplify(dh(JM)-JM*lm) == 0)
    check("lapse_multiplier_first", s.simplify(dh(JA)-JA*la) == 0)
    check("mass_log_mixed", s.simplify(dk(lm)-lmk) == 0)
    check("lapse_log_mixed", s.simplify(dk(la)-lak) == 0)
    check("mass_multiplier_mixed", s.simplify(dk(dh(JM))-JM*(lm*lm.xreplace({hb:kb,hV:kV,hz:kz})+lmk)) == 0)
    check("lapse_multiplier_mixed", s.simplify(dk(dh(JA))-JA*(la*la.xreplace({hb:kb,hV:kV,hz:kz})+lak)) == 0)
    C,r0,rp0,F0 = s.symbols("C r0 rp0 F0", positive=True)
    r = r0*T
    rp = rp0*s.exp(V)*(1+z)
    F = F0*s.exp(2*V-2*b)*(1+z)**2/T**2
    check("mass_density_exact", s.simplify(rp/r**2/(rp0/r0**2)-JM) == 0)
    check("lapse_density_exact", s.simplify(rp/(r**3*F)/(rp0/(r0**3*F0))-JA) == 0)
    H0,Hh,Hk,Hhk,J0,Jh,Jk,Jhk = s.symbols("H0 Hh Hk Hhk J0 Jh Jk Jhk")
    source = H0+ss*Hh+tt*Hk+ss*tt*Hhk
    geom = J0+ss*Jh+tt*Jk+ss*tt*Jhk
    check("full_source_geometry_product_rule", s.diff(source*geom,ss,tt).subs({ss:0,tt:0})
          == H0*Jhk+Hh*Jk+Hk*Jh+Hhk*J0)
    c0,ch,ck,chk,D,x = s.symbols("c0 ch ck chk D x")
    cut = c0+ss*ch+tt*ck+ss*tt*chk
    coeff = 1+ss+tt**2+ss*tt
    lower,upper = coeff*cut**3/3, coeff*(D**3-cut**3)/3
    zero={ss:0,tt:0}
    f = coeff*x**2
    boundary = (s.diff(f,ss).subs(zero).subs(x,c0)*ck
                +s.diff(f,tt).subs(zero).subs(x,c0)*ch
                +s.diff(f,x).subs(zero).subs(x,c0)*ch*ck
                +f.subs(zero).subs(x,c0)*chk)
    check("moving_lower_boundary_chain", s.simplify(s.diff(lower,ss,tt).subs(zero)-c0**3/3-boundary) == 0)
    check("moving_upper_boundary_chain", s.simplify(s.diff(upper,ss,tt).subs(zero)-(D**3-c0**3)/3+boundary) == 0)
    check("moving_interface_full_cancellation", s.simplify(s.diff(lower+upper,ss,tt).subs(zero)-D**3/3) == 0)
    eps = s.symbols("eps", positive=True)
    cuspquot = (s.Rational(2,3)*(2*eps)**s.Rational(3,2))/(4*eps**2)
    check("rough_split_second_difference", s.simplify(cuspquot-s.sqrt(2)/(3*s.sqrt(eps))) == 0)
    check("rough_lower_piece_not_c2", s.limit(cuspquot,eps,0,dir="+") == s.oo)
    q = s.symbols("q",positive=True)
    for power in [2,5]:
        antiderivative = -q**(1-power)*((1+s.log(q))**2/(power-1)
                        +2*(1+s.log(q))/(power-1)**2+s.Rational(2,(power-1)**3))
        check(f"outer_integral_primitive_{power}", s.simplify(s.diff(antiderivative,q)-(1+s.log(q))**2/q**power) == 0)
    sc = s.symbols("sc", positive=True)
    check("middle_primitive_maximum", s.simplify(s.diff((q-sc)/q**2,q).subs(q,2*sc)) == 0
          and s.simplify(((q-sc)/q**2).subs(q,2*sc)-1/(4*sc)) == 0)
    chi,h,B,G = s.symbols("chi h B G",positive=True)
    rt=B*(1+h*s.cosh(chi)); rpt=s.diff(rt,chi)
    Ft=C*(h*s.sinh(chi)/(1+h*s.cosh(chi)))**2
    phi=s.symbols("phi")
    normalized=G*C**2*h*phi*rpt/(rt**3*Ft)/(G*C/B**2)
    check("toy_gap_cancellation_exact", s.simplify(normalized-phi/((1+h*s.cosh(chi))*s.sinh(chi))) == 0)
    check("missing_modulus_bolt_singularity", s.limit(chi*normalized.subs(phi,1),chi,0,dir="+") == 1/(1+h))
    check("dini_power_integral", s.integrate(chi**(-s.Rational(1,2)),(chi,0,1)) == 2)
    check("missing_gap_factor_exact", s.simplify(normalized/h-phi/(h*(1+h*s.cosh(chi))*s.sinh(chi))) == 0)


def budget_checks():
    mp.mp.dps=PROTOCOL["digits"]
    eps=Q(1,1000); inv=1/(1-eps)
    m1,a1=1+inv,4+inv
    m2,a2=m1*m1+Q(1,2)+inv**2,a1*a1+Q(1,4)+inv**2
    jM,jA=Q(1003,1000),Q(1006,1000)
    check("mass_value_multiplier_bound", mp.exp(mp.mpf(".001"))/mp.mpf(".999") < mp.mpf("1.003"))
    check("lapse_value_multiplier_bound", mp.exp(mp.mpf(".004"))/mp.mpf(".999") < mp.mpf("1.006"))
    check("mass_first_multiplier_bound", jM*m1<Q(201,100))
    check("mass_mixed_multiplier_bound", jM*m2<Q(554,100))
    check("lapse_first_multiplier_bound", jA*a1<Q(504,100))
    check("lapse_mixed_multiplier_bound", jA*a2<Q(265,10))
    check("cold_mass_primitive_bound", Q(1001,1000)*81/Q(18,10)<46)
    check("middle_mass_primitive_bound", 1/(Q(36,10)*Q(1001,1000)*Q(1,1000))<279)
    L=mp.log(40)
    ie=(L**2+4*L+5)/40
    ia=( (1+L)**2/4+(1+L)/8+mp.mpf(1)/32 )/40**4
    check("outer_absolute_mass_integral", ie<mp.mpf(".835"),mp.nstr(ie,30))
    check("outer_mass_primitive_bound", mp.mpf("12.2")*ie/40**3<mp.mpf(".00016"))
    check("outer_lapse_integral", mp.mpf("12.2")*ia<mp.mpf(".00003"),mp.nstr(mp.mpf("12.2")*ia,30))
    coldmass=mp.mpf(".001")**2*(mp.mpf("1.5")+mp.log(10))/2
    check("cold_absolute_mass_integral", coldmass<mp.mpf(".000002"))
    check("cold_bolt_lapse_bound", 6*81==486)
    check("cold_remaining_lapse_bound", 2*81**2==13122)
    check("regional_common_envelopes", 46+279+Q(16,100000)<330
          and Q(2,10**6)+1+Q(835,1000)<2
          and 13122+486+1300+Q(3,100000)<15000)
    check("mixed_profile_cross_coefficients", 2*Q(201,100)==Q(402,100) and 2*Q(504,100)==Q(1008,100))
    check("thermal_tail_guard", 2*Q(1,10**60)*Q(1001,1000)/Q(1,10**6)<Q(1,10**53))
    eta=Q(85,10**34)
    NUMERICS.append({"regional_matrix":[["46","0","279",".00016"],
                                      [".000002","0","1",".835"],
                                      ["13122","486","1300",".00003"]],
                     "common_envelopes_times_eta":[str(a*eta) for a in [330,2,15000]],
                     "actual_quantum_source_jet_constants_proved":False,
                     "full_equation_inverse_proved":False})


def numerical_checks():
    mp.mp.dps=PROTOCOL["digits"]
    B,C,G=mp.mpf(2),mp.mpf(3),mp.mpf(".01")
    endpoints=list(map(mp.mpf,PROTOCOL["endpoints"]))
    D=endpoints[-1]
    all_results={}

    def fields(x):
        b=[mp.mpf(".0002")+mp.mpf(".00001")*x*x,
           mp.mpf(".03")+mp.mpf(".004")*x*x,-mp.mpf(".02")+mp.mpf(".003")*x*x,mp.mpf(0)]
        v=[mp.mpf(".0003")+mp.mpf(".00002")*x*x,
           mp.mpf(".02")+mp.mpf(".003")*x*x,-mp.mpf(".015")+mp.mpf(".002")*x*x,mp.mpf(0)]
        xcoth=x*mp.coth(x) if x else mp.mpf(1)
        z=[2*mp.mpf(a)*xcoth for a in [".00002",".003",".002","0"]]
        rho=[mp.mpf(".8")+x*x/10,mp.mpf(".2")+mp.mpf(".03")*x*x,
             -mp.mpf(".1")+mp.mpf(".02")*x*x,mp.mpf(".11")+mp.mpf(".01")*x*x]
        nn=[1+mp.mpf(".05")*x*x,mp.mpf(".1")+mp.mpf(".02")*x*x,
            mp.mpf(".07")-mp.mpf(".01")*x*x,mp.mpf(".03")+mp.mpf(".005")*x*x]
        return b,v,z,rho,nn

    def expand(v,ss,tt):
        return v[0]+ss*v[1]+tt*v[2]+ss*tt*v[3]

    for count in PROTOCOL["nodes"]:
        nodes,weights=mp.gauss_quadrature(count,"legendre")
        def quad(f,end):
            return sum(w*end/2*f(end*(1+x)/2) for x,w in zip(nodes,weights))
        for gaptext in PROTOCOL["gaps"]:
            gap=mp.mpf(gaptext)
            def seed(x):
                r0=B*(1+gap*mp.cosh(x))
                rp0=B*gap*mp.sinh(x)
                f0=C*(rp0/r0)**2
                return r0,rp0,f0
            def fixed_jet(x,freeze=False):
                b,v,z,rho,nn=fields(x)
                if freeze: b,v,z=[a[:1]+[0,0,0] for a in [b,v,z]]
                b,v,z=Jet(*b),Jet(*v),Jet(*z)
                r0,rp0,f0=seed(x)
                T=(1+gap*mp.cosh(x)*v.exp())/(1+gap*mp.cosh(x))
                JM=v.exp()*(1+z)/(T*T)
                JA=(2*b-v).exp()/(T*(1+z))
                HR=Jet(*rho)*(C*C*gap*mp.cosh(x))
                HN=Jet(*nn)*(C*C*gap*x*x/(1+x*x))
                return JM*HR*(G*rp0/r0**2),JA*HN*(G*rp0/(r0**3*f0))
            def direct_scalar(x,ss,tt):
                if not x: return mp.mpf(0),mp.mpf(0)
                b,v,z,rho,nn=[expand(a,ss,tt) for a in fields(x)]
                r=B*(1+gap*mp.cosh(x)*mp.exp(v))
                rp=B*gap*mp.sinh(x)*mp.exp(v)*(1+z)
                F=C*mp.exp(-2*b)*(rp/r)**2
                HR=C*C*gap*mp.cosh(x)*rho
                HN=C*C*gap*x*x/(1+x*x)*nn
                return G*HR*rp/r**2,G*HN*rp/(r**3*F)
            def normalizer(end):
                r0,_,f0=seed(end)
                return r0*f0
            def maps_scalar(ss,tt):
                mass=[quad(lambda x:direct_scalar(x,ss,tt)[0],end) for end in endpoints]
                lapse=quad(lambda x:direct_scalar(x,ss,tt)[1],D)
                # Scale the toy infinity mass by gap^2 to exercise its
                # response at ordinary size; this is only a test scale.
                return [mass[0]/normalizer(endpoints[0]),mass[1]/normalizer(D),
                        mass[1]/(B*C*gap**2),lapse]
            def maps_jet(freeze=False):
                mass=[quad(lambda x:fixed_jet(x,freeze)[0],end) for end in endpoints]
                lapse=quad(lambda x:fixed_jet(x,freeze)[1],D)
                return [mass[0]/normalizer(endpoints[0]),mass[1]/normalizer(D),
                        mass[1]/(B*C*gap**2),lapse]
            jets=maps_jet()
            base=maps_scalar(mp.mpf(0),mp.mpf(0))
            targets=[a.v[3] for a in jets]
            tag=f"{count}_{gaptext}"
            check(f"direct_physical_vs_fixed_values_{tag}", all(abs(a-b.v[0])<mp.mpf("1e-11")*(1+abs(a)) for a,b in zip(base,jets)))
            errors=[[] for _ in targets]
            firsterrors=[[] for _ in targets]
            for steptext in PROTOCOL["steps"]:
                step=mp.mpf(steptext)
                four=[maps_scalar(ss,tt) for ss,tt in [(step,step),(step,-step),(-step,step),(-step,-step)]]
                two=[maps_scalar(ss,mp.mpf(0)) for ss in [step,-step]]
                for j,target in enumerate(targets):
                    errors[j].append(abs((four[0][j]-four[1][j]-four[2][j]+four[3][j])/(4*step**2)-target))
                    firsterrors[j].append(abs((two[0][j]-two[1][j])/(2*step)-jets[j].v[1]))
            for name,errs,target in zip(["mass_inner","mass_outer","mass_total_scaled","lapse"],errors,targets):
                check(f"{name}_mixed_difference_{tag}", max(errs)<mp.mpf("2e-6")*(1+abs(target)),
                      {"target":mp.nstr(target,25),"errors":[mp.nstr(e,16) for e in errs]})
                check(f"{name}_mixed_refinement_{tag}", errs[1]<mp.mpf(".35")*errs[0] or max(errs)<mp.mpf("1e-30"))
            check(f"all_first_differences_{tag}", all(max(err)<mp.mpf("2e-6")*(1+abs(j.v[1])) for err,j in zip(firsterrors,jets)))
            adaptive=mp.quad(lambda x:direct_scalar(x,0,0)[1],[0,D/2,D])
            check(f"independent_adaptive_lapse_{tag}",abs(adaptive-base[-1])<mp.mpf("1e-11")*(1+abs(adaptive)))
            frozen=maps_jet(True)
            check(f"mass_jacobian_response_omission_detected_{tag}",abs(frozen[0].v[3]-targets[0])>mp.mpf("1e-10"))
            check(f"lapse_jacobian_response_omission_detected_{tag}",abs(frozen[-1].v[3]-targets[-1])>mp.mpf("1e-10"))
            all_results[(count,gaptext)]=targets
            NUMERICS.append({"nodes":count,"gap":gaptext,"mixed_targets":[mp.nstr(a,30) for a in targets],
                             "toy_infinity_mass_test_divisor":"B*C*h^2",
                             "manufactured_profiles_are_physical_sources":False})
    for gaptext in PROTOCOL["gaps"]:
        check(f"quadrature_refinement_{gaptext}",all(
            abs(a-b)<mp.mpf("1e-11")*(1+abs(b))
            for a,b in zip(all_results[(32,gaptext)],all_results[(64,gaptext)])))
    cold=[]
    for gaptext in PROTOCOL["gaps"]:
        gap=mp.mpf(gaptext)
        val=mp.quad(lambda x:x*x/((1+x*x)*(1+gap*mp.cosh(x))*mp.sinh(x)) if x else 0,[0,1])
        cold.append(val)
    check("cold_weight_uniform_small_gap",abs(cold[1]-cold[0])<mp.mpf("1e-5") and min(cold)>mp.mpf(".1"))
    check("dropping_gap_weight_loses_uniformity",
          (cold[1]/mp.mpf("1e-30"))/(cold[0]/mp.mpf("1e-5"))>mp.mpf("1e24"))
    quot=[mp.sqrt(2)/(3*mp.sqrt(mp.mpf(e))) for e in [".001",".0001",".00001"]]
    check("rough_split_divergence_refinement",all(b>3*a for a,b in zip(quot,quot[1:])))
    for steptext in PROTOCOL["steps"]:
        step=mp.mpf(steptext)
        integrals=[]
        for cut in [2*step,mp.mpf(0),mp.mpf(0),-2*step]:
            lower=mp.quad(mp.sqrt,[0,cut]) if cut>0 else mp.mpf(0)
            upper=mp.quad(mp.sqrt,[max(cut,mp.mpf(0)),mp.mpf(".1"),1])
            integrals.append((lower,upper))
        lower_difference=(integrals[0][0]-integrals[1][0]-integrals[2][0]+integrals[3][0])/(4*step**2)
        whole_difference=sum(sign*(lo+hi) for sign,(lo,hi) in zip([1,-1,-1,1],integrals))/(4*step**2)
        check(f"rough_quadrature_interface_cancellation_{steptext}",
              lower_difference>10 and abs(whole_difference)<mp.mpf("1e-11"),
              {"lower_mixed":mp.nstr(lower_difference,25),"whole_mixed":mp.nstr(whole_difference,25)})


def hashes():
    return {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    global Jet
    before=hashes()
    error=None
    try:
        spec=importlib.util.spec_from_file_location("s152_jet",ROOT/INPUTS[6])
        dependency=importlib.util.module_from_spec(spec); spec.loader.exec_module(dependency)
        Jet=dependency.Jet
        predecessor_checks()
        exact_checks()
        budget_checks()
        numerical_checks()
    except Exception:
        error=traceback.format_exc()
        print(error,file=sys.stderr,flush=True)
    after=hashes()
    ok=error is None and before==after and all(r["passed"] for r in ROWS)
    record={"story":"S152","created_utc":datetime.now(timezone.utc).isoformat(),
            "status":"passed" if ok else "failed","python":sys.version,
            "mpmath":mp.__version__,"sympy":s.__version__,"protocol":PROTOCOL,
            "input_sha256_before":before,"input_sha256_after":after,"unchanged_inputs":before==after,
            "checks":ROWS,"numerics":NUMERICS,"exception":error,
            "scope":"global fixed-chart source-to-constraint map on stated regional profile norms",
            "unresolved":["actual quantum source jets in the weighted space","full higher-order inverse",
                          "coupled exterior/state, charge and signed source","all physical completion requirements"]}
    if error:
        record["failed_checker_source_utf8"]=Path(__file__).read_text(encoding="utf-8")
        record["failed_note_source_utf8"]=(ROOT/INPUTS[1]).read_text(encoding="utf-8")
    out=ROOT/"evidence"/f"{STEM}.json"
    history=json.loads(out.read_text(encoding="utf-8")) if out.exists() else {"runs":[]}
    history["runs"].append(record)
    out.write_text(json.dumps(history,indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps({"status":record["status"],"checks":len(ROWS),"unchanged_inputs":before==after}),flush=True)
    return 0 if ok else 1


if __name__=="__main__":
    raise SystemExit(main())
