"""S153: real endpoint smallness -> cold canonical response.

Analytic manufactured families verify the transfer, not the Dirac determinant.
"""
from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback

import mpmath as mp
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
STEM="cold-response-from-real-endpoints"
INPUTS=[
    "checks/cold_response_from_real_endpoints.py",
    "notes/cold-response-from-real-endpoints.md",
    "notes/general-interior-dirac-response.md",
    "notes/stationary-centered-canonical-response.md",
    "notes/interior-pressure-and-bolt-primitive.md",
    "notes/fixed-chart-global-constraint-transfer.md",
    "evidence/fixed-chart-global-constraint-transfer.json",
    "notes/retained-benchmark-parameter-recovery.md",
]
ROWS,NUMERICS=[],[]
PROTOCOL={"digits":60,"nodes":[64,128],"steps":["1e-4","5e-5"],
          "thermal_exponents":[20,50,79],"bolt_exponents":[0,20,80],
          "difference_tolerance":"1e-7","contour_tolerance":"1e-20",
          "refinement_ratio":".35","refinement_floor":"1e-40",
          "diagnostics_normalized_by_amplitude":True,
          "difference_frequency_divisor":16}


def check(name,passed,detail=""):
    assert not any(r["name"]==name for r in ROWS)
    row={"name":name,"passed":bool(passed),"detail":str(detail)}
    ROWS.append(row); print(json.dumps(row),flush=True)
    if not passed: raise AssertionError(name)


def predecessor_checks():
    last=json.loads((ROOT/INPUTS[6]).read_text(encoding="utf-8"))["runs"][-1]
    check("s152_final_receipt",last["status"]=="passed" and len(last["checks"])==99
          and all(r["passed"] for r in last["checks"]))
    check("s152_current_inputs",last["unchanged_inputs"] and all(
        sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in last["input_sha256_after"].items()))


def exact_checks():
    x,y=s.symbols("x y",real=True)
    omega=2/s.pi*s.atan(2*y/(1-x*x-y*y))
    check("half_disk_harmonic_measure",s.simplify(s.diff(omega,x,2)+s.diff(omega,y,2))==0)
    check("harmonic_measure_diameter",omega.subs(y,0)==0)
    rr=s.symbols("rr",positive=True)
    imaginary_axis=4/s.pi*s.atan(rr)
    check("harmonic_measure_arc_limit",s.limit(imaginary_axis,rr,1,dir="-")==1)
    check("logarithmic_map_derivative",s.simplify(s.diff(s.log((1+x)/(1-x)),x)-2/(1-x*x))==0)
    n=s.symbols("n",positive=True)
    f=s.exp(-n)*s.cos(n*x)
    check("log_loss_second_derivative",s.simplify(s.diff(f,x,2).subs(x,0)/s.exp(-n)+n*n)==0)
    check("constant_times_real_bound_insufficient",s.limit(n*n,n,s.oo)==s.oo)
    t=s.symbols("t",nonnegative=True)
    for order in [1,2]:
        weight=s.exp(-t/4)*(90+t/2)**order
        ratio=s.factor(s.diff(weight,t)/weight)
        check(f"bolt_log_power_monotone_{order}",
              s.simplify(ratio-(s.Rational(order,2)/(90+t/2)-s.Rational(1,4)))==0
              and s.Rational(order,180)<s.Rational(1,4))
    A,B,C,uu,vv=s.symbols("A B C uu vv")
    quadratic=lambda a,b:A*a*a+2*B*a*b+C*b*b
    check("mixed_polarization",s.expand((quadratic(uu,vv)-quadratic(uu,-vv))/4-B*uu*vv)==0)
    q=s.symbols("q",positive=True)
    L=s.symbols("L",positive=True)
    # f_q'' / sqrt(q) diverges, whereas the quarter-weighted response vanishes.
    second_over_sqrt=(L-s.log(q)/2)**2
    check("square_root_response_space_not_inferred",s.limit(second_over_sqrt,q,0,dir="+")==s.oo)
    check("quarter_response_space_tail_vanishes",s.limit(q**s.Rational(1,4)*second_over_sqrt,q,0,dir="+")==0)
    check("quarter_bolt_integral",s.integrate(q**(-s.Rational(3,4)),(q,0,1))==4)
    check("log_bolt_integral_zeroth",s.integrate(s.exp(-t/2),(t,0,s.oo))==2)
    check("log_bolt_integral_second",s.integrate(s.exp(-t/2)*(L+t/2)**2,(t,0,s.oo))==2*L**2+4*L+4)
    z=s.symbols("z",real=True)
    alpha=s.Rational(1,2)-s.I*z/32
    analytic_trace=s.exp(s.I*L*z/16)*q**alpha
    analytic_pressure=analytic_trace/(alpha+2)
    check("analytic_ward_primitive",s.simplify(s.diff(q*q*analytic_pressure,q)-q*analytic_trace)==0)
    integrated=(2/(alpha+2)-1)*s.exp(s.I*L*z/16)/alpha
    check("analytic_integrated_null_source",s.simplify(
        (2*analytic_pressure-analytic_trace)/q
        -s.diff((2/(alpha+2)-1)*s.exp(s.I*L*z/16)*q**alpha/alpha,q))==0)
    NUMERICS.append({"toy_integrated_null_formula":str(integrated),
                     "manufactured_family_is_Dirac_source":False})


def budget_checks():
    mp.mp.dps=PROTOCOL["digits"]
    E=mp.mpf("1e-6"); av=mp.mpf("1.7e-15"); ao=mp.mpf("1.7e-14"); BC=mp.mpf("2e-11")
    check("harmonic_measure_linear_reserve",16/(3*mp.pi)<2)
    check("unit_cauchy_radius_stays_inside",mp.mpf(1)/(2*90)<mp.mpf(".5"))
    check("value_thermal_log_reserve",80+mp.log(BC/av)<90)
    check("oscillation_thermal_log_reserve",80+mp.log(BC/ao)<90)
    check("real_bounds_below_complex_bound",av*mp.mpf(".00331")<BC and ao*mp.mpf("1e-30")*81<BC)
    first=2*mp.e*av*90/(E/2)
    pure=8*mp.e*av*90**2/(E/2)**2
    check("scalar_first_budget",first<mp.mpf("2e-6"),mp.nstr(first,25))
    check("scalar_quadratic_budget",pure<1200,mp.nstr(pure,25))
    check("scalar_mixed_polarization_budget",2*pure<2500)
    check("scalar_oscillation_first_budget",10*first<mp.mpf("2e-5"))
    check("scalar_oscillation_mixed_budget",20*pure<25000)
    eta=Q(85,10**34)
    seed=2*eta*10**20*400**2
    check("seed_cold_value_budget",seed+Q(17,10**16)<Q(3,10**7))
    check("seed_bolt_oscillation_budget",seed*100*8*2+Q(17,10**15)<Q(5,10**4))
    ac=Q(331,10**5); ab=Q(2,10**28)
    S0,S1,S2=Q(3,10**7),Q(2,10**6),Q(2500)
    P0=Q(12,10)*S0
    P1=Q(12,10)*S1+700*ac*S0
    P2=Q(12,10)*S2+1400*ac*S1+5500*ac*ac*S0
    check("cold_pressure_first_budget",P1<Q(32,10**7),str(P1))
    check("cold_pressure_mixed_budget",P2<3001,str(P2))
    check("cold_profile_value_budget",2*P0+S0<=Q(11,10**7))
    check("cold_profile_first_budget",2*P1+S1<Q(9,10**6))
    check("cold_profile_mixed_budget",2*P2+S2<8600)
    O0,O1,O2=Q(5,10**4),Q(2,10**5),Q(25000)
    check("bolt_profile_value_budget",Q(34,10)*O0==Q(17,10**4))
    check("bolt_profile_first_budget",Q(34,10)*O1+1400*ab*O0<Q(7,10**5))
    check("bolt_profile_mixed_budget",Q(34,10)*O2+2800*ab*O1+11000*ab*ab*O0<86000)
    check("quarter_bolt_transfer_coefficient",3*81*4==972)
    m2=Q(1003,1000)*8600+Q(402,100)*Q(9,10**6)+Q(554,100)*Q(11,10**7)
    c2=Q(1006,1000)*8600+Q(1008,100)*Q(9,10**6)+Q(265,10)*Q(11,10**7)
    b2=Q(1006,1000)*86000+Q(1008,100)*Q(7,10**5)+Q(265,10)*Q(17,10**4)
    QM=46*m2/2; QI=Q(2,10**6)*m2/2; QA=(13122*c2+972*b2)/2
    check("cold_component_mass_mixed",QM<200000,str(QM))
    check("cold_component_infinity_mixed",QI<Q(87,10000),str(QI))
    check("cold_component_lapse_mixed",QA<10**8,str(QA))
    check("coarse_bound_is_not_a_contraction_claim",Q(10**8,2)*Q(1,10**6)**2==Q(5,10**5))
    check("conditional_small_step_taylor",Q(10**8,2)*Q(1,10**10)**2==Q(5,10**13))
    NUMERICS.append({"mixed_cold_constraint_bounds":{"mass":str(QM),"infinity":str(QI),"lapse":str(QA)},
                     "full_inverse_or_contraction_proved":False,
                     "constructed_endpoint_admissibility_proved":False})


def numerical_checks():
    mp.mp.dps=PROTOCOL["digits"]
    BC=mp.mpf("2e-11")
    for radius in [mp.mpf(".01"),mp.mpf(".1"),mp.mpf(".5")]:
        vals=[]
        for j in range(1,32):
            z=radius*mp.exp(mp.j*mp.pi*j/32)
            omega=2/mp.pi*mp.im(mp.log((1+z)/(1-z)))
            vals.append(omega)
        check(f"harmonic_measure_circle_{radius}",min(vals)>0 and max(vals)<2*radius)
    contours={}
    for thermal in PROTOCOL["thermal_exponents"]:
        h=mp.exp(-thermal); Lh=1+mp.log(mp.mpf(".01")/h)
        for bolt in PROTOCOL["bolt_exponents"]:
            chi=mp.exp(-bolt)
            amp=mp.mpf("1.7e-14")*h*Lh*mp.sqrt(chi)
            loss=mp.log(BC/amp)
            freq=loss/PROTOCOL["difference_frequency_divisor"]
            radius=1/(2*max(mp.mpf(1),loss))
            # All derivative diagnostics are divided by amp.
            fun=lambda z:mp.mpf(".6")*mp.cos(freq*z)+mp.mpf(".4")*mp.sin(freq*z)
            targets=[mp.mpf(".4")*freq,-mp.mpf(".6")*freq**2]
            tag=f"{thermal}_{bolt}"
            err1=[];err2=[]
            for steptext in PROTOCOL["steps"]:
                step=mp.mpf(steptext)
                err1.append(abs((fun(step)-fun(-step))/(2*step)-targets[0]))
                err2.append(abs((fun(step)-2*fun(0)+fun(-step))/step**2-targets[1]))
            for name,errs,target in zip(["first","second"],[err1,err2],targets):
                check(f"analytic_{name}_difference_{tag}",max(errs)<mp.mpf("1e-7")*(1+abs(target)),
                      {"normalized_target":mp.nstr(target,25),"errors":[mp.nstr(a,16) for a in errs]})
                check(f"analytic_{name}_step_refinement_{tag}",errs[1]<mp.mpf(".35")*errs[0]
                      or max(errs)<mp.mpf("1e-40"))
            for count in PROTOCOL["nodes"]:
                values=[radius*mp.exp(2*mp.pi*mp.j*j/count) for j in range(count)]
                observed=[mp.factorial(order)*sum(fun(z)/z**order for z in values)/count for order in [1,2]]
                check(f"cauchy_contour_derivatives_{count}_{tag}",all(
                    abs(v-target)<mp.mpf("1e-20")*(1+abs(target)) for v,target in zip(observed,targets)))
                check(f"small_circle_value_budget_{count}_{tag}",max(abs(fun(z)) for z in values)<mp.e)
                contours[(count,thermal,bolt)]=observed
            check(f"cauchy_node_refinement_{tag}",all(abs(a-b)<mp.mpf("1e-20")*(1+abs(b))
                  for a,b in zip(contours[(64,thermal,bolt)],contours[(128,thermal,bolt)])))
            # This check uses the full frequency, independent of the
            # finite-difference frequency reduction.
            full_second=amp*loss**2
            check(f"quarter_weight_dominates_full_log_response_{tag}",
                  full_second/(mp.mpf("1.7e-14")*h*Lh*chi**mp.mpf(".25"))<90**2)
            NUMERICS.append({"thermal_exponent":thermal,"bolt_exponent":bolt,
                             "loss":mp.nstr(loss,25),"normalized_second_target":mp.nstr(targets[1],25)})
    # Analytic necessity: uniform complex bound with arbitrary real smallness.
    growth=[mp.mpf(n)**2 for n in [20,50,100]]
    check("logarithmic_loss_counterexample_growth",growth==[400,2500,10000])
    log_weights=[(90+mp.mpf(t)/2)**2 for t in [0,20,80,400]]
    quarter_weights=[mp.exp(-mp.mpf(t)/4)*(90+mp.mpf(t)/2)**2 for t in [0,20,80,400]]
    check("square_root_weight_failure",all(b>a for a,b in zip(log_weights,log_weights[1:])))
    check("quarter_weight_tail_control",all(b<a for a,b in zip(quarter_weights,quarter_weights[1:])))
    # Exactly soluble Ward family, integrated in logarithmic distance.
    for Lam in [20,50,90]:
        lam=mp.mpf(Lam)
        def complex_integrated(z):
            alpha=mp.mpf(".5")-mp.j*z/32
            return (2/(alpha+2)-1)*mp.exp(mp.j*lam*z/16)/alpha
        exact=mp.re(mp.diff(complex_integrated,mp.mpf(0),2))
        def density_log(t):
            # n=lam/16+t/32; pressure denominator is 5/2-i z/32.
            n=lam/16+t/32
            c0=mp.mpf("2")/mp.mpf("2.5")-1
            c1=2*mp.j/(32*mp.mpf("2.5")**2)
            c2=-4/(32**2*mp.mpf("2.5")**3)
            second=mp.re(-n*n*c0+2*mp.j*n*c1+c2)
            return mp.exp(-t/2)*second
        observed=mp.quad(density_log,[0,4,20,mp.inf])
        absolute=mp.quad(lambda t:abs(density_log(t)),[0,4,20,mp.inf])
        bound=mp.quad(lambda t:mp.exp(-t/2)*(lam+t/2)**2,[0,4,20,mp.inf])
        check(f"ward_logarithmic_integral_{Lam}",abs(exact-observed)<mp.mpf("1e-20")*(1+abs(exact)),
              {"signed":mp.nstr(observed,25),"absolute":mp.nstr(absolute,25)})
        check(f"ward_absolute_bolt_integral_{Lam}",absolute<bound and mp.isfinite(absolute))
        errs=[]
        for steptext in PROTOCOL["steps"]:
            step=mp.mpf(steptext)
            fd=mp.re((complex_integrated(step)-2*complex_integrated(0)+complex_integrated(-step))/step**2)
            errs.append(abs(fd-exact))
        check(f"ward_second_difference_{Lam}",max(errs)<mp.mpf("1e-7")*(1+abs(exact)))
        check(f"ward_difference_refinement_{Lam}",errs[1]<mp.mpf(".35")*errs[0]
              or max(errs)<mp.mpf("1e-40"))


def hashes():
    return {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    before=hashes();error=None
    try:
        predecessor_checks();exact_checks();budget_checks();numerical_checks()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes()
    ok=error is None and before==after and all(r["passed"] for r in ROWS)
    record={"story":"S153","created_utc":datetime.now(timezone.utc).isoformat(),
            "status":"passed" if ok else "failed","python":sys.version,
            "mpmath":mp.__version__,"sympy":s.__version__,"protocol":PROTOCOL,
            "input_sha256_before":before,"input_sha256_after":after,"unchanged_inputs":before==after,
            "checks":ROWS,"numerics":NUMERICS,"exception":error,
            "scope":"canonical cold response on intersected real endpoint and fixed-domain geometry contracts",
            "unresolved":["admissible global solver family","outer recombined and other sector responses",
                          "full higher-order inverse and coupled exterior/state solve","all physical completion gates"]}
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
