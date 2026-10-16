"""Frozen exact checks for the centered seed and complete horizon source.

Heat, complex-domain and all-order remainder estimates are ordinary proofs.
No numerical sampling certifies them, and no corrected exterior is solved.
"""
from __future__ import annotations

import base64
import contextlib
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sys
import traceback

import sympy as S

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "evidence/stationary-centered-horizon.json"
LOG = ROOT / "evidence/stationary-centered-horizon.log"
INPUTS = [
    "checks/stationary_centered_horizon.py",
    "notes/stationary-centered-covariant-seed.md",
    "notes/stationary-horizon-source-defect.md",
    "notes/stationary-reference-cold-bolt.md",
    "notes/stationary-reference-nonbolt-stress.md",
    "notes/stationary-weyl-radial-sources.md",
    "notes/stationary-horizon-product-matching.md",
    "notes/uniform-dirac-reference.md",
    "notes/spherical-anomaly-action.md",
]
CONTRACT = {
    "acceptance": "all symbolic residuals exactly zero and all rational inequalities true",
    "numerical_tolerances": "none",
    "scope": "exact centered covariant seed and both-pressure source VALUE on two explicit HH references",
    "source_bound": "complete Gamma_rest horizon forcing / (2 kappa/R) < 5.1e-12",
    "ordinary_arguments": ["analytic cap construction and complex normal disks",
                           "uniform complete Born and angular remainder bounds",
                           "actual-domain localization", "HH regularity and state allocation"],
    "excluded": ["corrected exterior existence", "full metric response estimate",
                 "bounded solver inverse", "full boundary Hamiltonian", "physical gate closure"],
    "failure_policy": "append all runs/logs and preserve every failed run's input bytes",
}
ROWS = []
R = S.Rational


def check(label, condition, detail=""):
    row = {"label": label, "passed": bool(condition), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row, sort_keys=True), flush=True)


def zero(label, expr):
    value = S.factor(S.together(expr))
    check(label, value == 0, value)


def seed_identities():
    t, x, d, delta, B, A, G, Q, le, l0 = S.symbols(
        "t x d delta B A G Q lambda_e lambda0", positive=True)
    C = 1+6*x+d
    ext = (t-1)**2/t**2*(1+x*(t*t+2*t+3)+d/t**2)
    f = ext-delta/t
    P = t*ext
    zero("extremal_double_root_value", ext.subs(t,1))
    zero("extremal_double_root_slope", S.diff(ext,t).subs(t,1))
    zero("independent_orbit_curvature", S.diff(ext,t,2).subs(t,1)/2-C)
    bracket = (t+1)*(1+3*x*(t*t+1))/t**2-d*(t-3)/t**4
    zero("global_exterior_positive_factor", S.diff(P,t)-(t-1)*bracket)
    v = S.symbols("v",positive=True)
    kv = S.diff(P,t).subs(t,v)/(2*v)
    zero("thermal_root_clock", (S.diff(f,t)/2).subs(delta,P.subs(t,v)).subs(t,v)-kv)
    K = ((1+3*x)*(v**-2+v**-3)+3*x*(1+v**-1)
         -d*(v**-4-3*v**-5))/(2*C)
    zero("thermal_clock_factorization", kv-C*(v-1)*K)
    zero("thermal_factor_at_center", K.subs(v,1)-1)
    zero("global_temperature_derivative", S.diff(kv,v)
         -(3*x-v**-2+3*(1+3*x)*v**-4
           +d*(3*v**-4-16*v**-5+15*v**-6))/2)
    # Potential derived by differentiating rF, before imposing a root.
    V0 = 2*S.diff(t*ext,t)
    Vphys = 2+6*l0*B**2*t**2-2*Q/(B**2*t**2)
    qe = B**2*(1+3*x)
    Vdisplay = 2+6*x*t**2-2*qe/(B**2*t**2) \
        +2*d*(-t**-2+4*t**-3-3*t**-4)
    zero("exact_covariant_potential", V0-Vdisplay)
    zero("thermal_member_same_potential", 2*S.diff(t*f,t)-V0)
    zero("center_potential_value", V0.subs(t,1))
    zero("center_potential_derivative", (S.diff(V0,t)/(2*B**2*t)).subs(t,1)-2*C/B**2)
    DV = V0-Vphys
    pA = DV/(16*S.pi*G*B**2*t**2)
    pS = S.diff(DV,t)/(32*S.pi*G*B**2*t)
    zero("joint_reference_Ward", pS-pA-t*S.diff(pA,t)/2)
    zero("potential_radius_partner", S.diff(DV,t)/(2*B**2*t)
         -16*S.pi*G*pS)
    # Do not confuse the geometric mass with the full-action charge.
    r = B*t
    mgeo = r*(1+l0*r*r+Q/r**2-f)/2
    Av = 3*(l0-x/B**2)/(2*G)
    mu = mgeo-G*Av*r**3/3
    zero("finite_seed_mass_asymptote", S.limit(mu,t,S.oo)-B*(1+2*x+delta/2))
    zero("vacuum_allocation_identity", l0-2*G*Av/3-x/B**2)
    # Exact Gamma_ref profiles after both product equations are used.
    p0A,p0S,pv,DP = S.symbols("p0A p0S pv DP")
    substitutions = {Q:qe-3*(x/B**2-l0)*B**4+8*S.pi*G*B**4*p0A,
                     d:8*S.pi*G*B**2*(p0A+p0S)-6*(x/B**2-l0)*B**2}
    pa = pA.subs(substitutions, simultaneous=True)
    ps = pS.subs(substitutions, simultaneous=True)
    profiles = [p0A/t**4+pv*(1-t**-4)+DP*(-t**-4+4*t**-5-3*t**-6),
                -p0A/t**4+pv*(1+t**-4)+DP*(t**-4-6*t**-5+6*t**-6)]
    repl = {DP:p0A+p0S-2*pv,pv:3*(x/B**2-l0)/(8*S.pi*G)}
    profile_values = [p.subs(DP,p0A+p0S-2*pv).subs(pv,repl[pv]) for p in profiles]
    zero("complete_reference_orbit_profile", pa-profile_values[0])
    zero("complete_reference_sphere_profile", ps-profile_values[1])
    zero("reference_center_orbit", profiles[0].subs(t,1)-p0A)
    zero("reference_center_sphere", profiles[1].subs(t,1).subs(DP,p0A+p0S-2*pv)-p0S)
    zero("profile_Ward_identity", profiles[1]-profiles[0]-t*S.diff(profiles[0],t)/2)


def cap_identities():
    s,h,x,d,e = S.symbols("s h x d e",positive=True)
    C = 1+6*x+d
    a,c = 4*x/C,x/C
    g = (1+s)**-2-1
    f = 1+a*s+c*s*s+d/C*g
    ext = s*s/(1+s)**2*(1+x*((1+s)**2+2*(1+s)+3)+d/(1+s)**2)
    zero("normalized_seed_cold_polynomial", (1+s)**2*ext/C-s*s*f)
    fh = f.subs(s,h)
    dd = S.cancel((f-fh)/(s-h))
    q0 = f+h*h/(s+h)*(dd-fh/(1+h))
    zero("cold_polynomial_factorization", s*s*f-h*h*fh*(1+s)/(1+h)-(s*s-h*h)*q0)
    gr = (1+s)**-2-1
    ddg = -(2+s+h)/((1+s)**2*(1+h)**2)
    zero("rational_extra_divided_difference", (gr-gr.subs(s,h))-(s-h)*ddg)
    qe = gr+h*h/(s+h)*(ddg+h*(h+2)/(1+h)**3)
    zero("extra_cold_factor", qe-(gr+h*h/(s+h)*(ddg-gr.subs(s,h)/(1+h))))
    ar,cr = 4*x/(1+6*x),x/(1+6*x)
    qp = ar*(s+h*h/((s+h)*(1+h)))+cr*(s*s+h*h-h**4/((s+h)*(1+h)))
    zero("RN_polynomial_dependence", qp-(ar*s+cr*s*s+h*h/(s+h)*(ar+cr*(s+h)-(ar*h+cr*h*h)/(1+h))))
    rn = 1-h*h/((s+h)*(1+h))+qp
    zero("seed_perturbation_decomposition", q0-rn-d/C*(qe-qp))
    # Seed horizon derivatives independently computed from the thermal metric.
    t = S.symbols("t",positive=True)
    ext_t = (t-1)**2/t**2*(1+x*(t*t+2*t+3)+d/t**2)
    vt = 1+h
    thermal = ext_t-vt*ext_t.subs(t,vt)/t
    w = (t*t*S.diff(thermal,t,2)/2).subs(t,vt)
    l3 = (t**3*S.diff(thermal,t,3)).subs(t,vt)
    wd = -1+2*(1+3*x+d)/vt**2-10*d/vt**3+9*d/vt**4
    ld = 6+6*x*vt**2-18*(1+3*x+d)/vt**2+108*d/vt**3-114*d/vt**4
    zero("exact_seed_horizon_second_jet",w-wd)
    zero("exact_seed_horizon_third_jet",l3-ld)
    zero("independent_center_curvature",w.subs(h,0)-C)
    zero("thermal_q0_clock",(S.diff(thermal,t).subs(t,vt)/2)-C*h*q0.subs(s,h)/(1+h)**2)
    # Constant physical Weyl and point-density normalization.
    z,sigma,M = S.symbols("z sigma M")
    Fz = S.exp(-z*sigma)
    zero("normalizer_density_Jacobian",Fz**2-S.exp(-2*z*sigma))
    zero("normalizer_endpoint_derivative",S.diff(Fz,z)+sigma*Fz)


def curvature_and_weyl():
    # Direct Christoffels/Riemann of the alpha=0 metric, not the displayed
    # generic horizon anomaly or the RN constant-scalar simplification.
    U,r,th,ph = S.symbols("U r theta phi",real=True)
    coordinates = (U,r,th,ph)
    F = S.Function("F")(r)
    diagonal = (-F,1/F,r*r,r*r*S.sin(th)**2)
    metric, inverse = S.diag(*diagonal),S.diag(*(1/v for v in diagonal))
    gam = {}
    for k in range(4):
        for i in range(4):
            for j in range(4):
                gam[k,i,j] = S.simplify(sum(inverse[k,l]*(S.diff(metric[l,j],coordinates[i])
                    +S.diff(metric[l,i],coordinates[j])-S.diff(metric[i,j],coordinates[l]))
                    for l in range(4))/2)
    sections = {}
    for i in range(4):
        for j in range(i+1,4):
            curv = S.diff(gam[i,j,j],coordinates[i])-S.diff(gam[i,j,i],coordinates[j]) \
                +sum(gam[i,i,l]*gam[l,j,j]-gam[i,j,l]*gam[l,j,i] for l in range(4))
            sections[i,j] = S.simplify(curv/diagonal[j])
    scalar = S.simplify(2*sum(sections.values()))
    ricci2 = sum(sum(sections[min(i,j),max(i,j)] for j in range(4) if j!=i)**2 for i in range(4))
    riem2 = 4*sum(v*v for v in sections.values())
    c2 = S.simplify(riem2-2*ricci2+scalar**2/3)
    eu = S.simplify(riem2-4*ricci2+scalar**2)
    box = F*S.diff(scalar,r,2)+(S.diff(F,r)+2*F/r)*S.diff(scalar,r)
    j,w,l3,rr,q,sig = S.symbols("j w l3 R q sigma",real=True)
    horizon = {F:0,S.diff(F,r):j/r,S.diff(F,r,2):2*w/r**2,S.diff(F,r,3):l3/r**3}
    def bolt(expr):
        return S.factor(S.expand(expr).subs(horizon,simultaneous=True).subs(r,rr))
    zero("direct_generic_scalar_at_bolt",bolt(scalar)*rr**2-(-2*w-4*j+2))
    zero("direct_generic_Ricci_square",bolt(ricci2)*rr**4-2*(w+j)**2-2*(1-j)**2)
    zero("direct_generic_Weyl_square",bolt(c2)*rr**4-R(4,3)*(w-j-1)**2)
    zero("direct_generic_Euler_density",bolt(eu)*rr**4-(4*j*j-8*w))
    zero("direct_generic_derivative_anomaly",bolt(box)*rr**4-j*(-l3-8*w+2*j-4))
    anomaly = (R(1,20)*c2-R(11,360)*eu+box/30-q*q/(3*r**4))
    product = (w*w+1)/15+w/9-q*q/3
    zero("complete_generic_anomaly_difference",bolt(anomaly)*rr**4-product-j*(j-36*w-3*l3)/90)
    # Derive the lapse Euler derivative from the complete off-shell Weyl
    # representative, taking all alpha variations BEFORE imposing alpha=0.
    alpha = S.Function("alpha")(r)
    T = r*(S.diff(F,r)+F*S.diff(alpha,r))
    Rbar = -r*r*S.diff(F,r,2)-3*r*r*S.diff(alpha,r)*S.diff(F,r) \
        -2*r*r*F*(S.diff(alpha,r,2)+S.diff(alpha,r)**2)+2*T-2*F
    sigma = S.log(r)  # b is fixed; log(r/b) differs only by its constant.
    L = S.exp(alpha)/r**2*(sigma*(Rbar**2/60-Rbar/18+R(1,15)-q*q/3)
        +Rbar*T/30-T*T/10+11*T*F/90-11*F*F/180-F/18)
    a0 = {alpha:0,S.diff(alpha,r):0,S.diff(alpha,r,2):0}
    # Substitution after partial variation, before radial differentiation,
    # is valid for the evaluation on the identically-zero lapse family.
    EA = S.diff(L,alpha).subs(a0,simultaneous=True) \
        -S.diff(S.diff(L,S.diff(alpha,r)).subs(a0,simultaneous=True),r) \
        +S.diff(S.diff(L,S.diff(alpha,r,2)).subs(a0,simultaneous=True),r,2)
    EAh = bolt(EA)*rr**2
    EAh = EAh.subs(S.log(rr),sig)
    generic = sig*((1-(w-j)**2)/15-q*q/3+j*l3/30)+j*((w-j)/15+R(1,18))
    zero("full_offshell_lapse_Euler_at_bolt",EAh-generic)
    prod_weyl = sig*((1-w*w)/15-q*q/3)
    delta_weyl = j*(sig*(4*w-2*j+l3)/30+(w-j)/15+R(1,18))
    zero("matched_product_generic_Weyl_transfer",generic-prod_weyl-delta_weyl)
    x,y = S.symbols("x y")
    zero("RN_generic_anomaly_specialization",(j*(j-36*w-3*l3)/90).subs(
        {w:2*y-1,l3:6+6*x-18*y}).subs(y,1+3*x-j)-j*(19*j-72*x)/90)


def rational_margins():
    k = S.Integer(2000)
    ratio = R(22,6250)  # .00352, valid on both explicit references.
    eps = R(1,10**24)
    # q0 perturbation: h<=s<=.01, the denominators are >=1 or s+h.
    # These inequalities are elementary positive-term bounds, not samples.
    qe_bound = 2+R(204,100)/2
    qp_bound = R(2,3)*R(3,2)+R(1,6)*R(2,100)
    check("q0_r_value_envelope",qe_bound+qp_bound<R(42,10))
    qe_prime = 2+R(204,100)/4+R(303,100)*R(1,200)
    qp_prime = R(2,3)+R(1,6)*(R(2,100)+R(1,40000))
    check("q0_r_derivative_envelope",qe_prime+qp_prime<R(16,5)<4)
    Rprime = R(16,5)+R(42,10)*R(1,100)*R(921,1000)
    check("J_seed_difference_envelope",Rprime*R(103,100)/(1-R(42,100)*eps)<5)
    check("cap_global_C2_perturbation",10**10*eps<R(1,10**14)+R(1,10**30))
    check("cap_seed_local_radius_guard",R(19,1000)<R(25,1000))
    check("cap_inverse_weighted_similarity",R(18,1)/(18-R(211,100))<R(6,5))
    # exp(.051)-1 < .053 from its positive series and a geometric tail.
    zz = R(51,1000)
    exp_upper = 1+zz+zz*zz/(2*(1-zz/3))
    check("global_seed_exponential_guard",exp_upper<1+R(131,2500))
    check("global_seed_offdiagonal_row",49*R(131,2500)/(1-8*R(131,2500))<R(9,2))
    check("local_amplitude_Born_remainder",R(14,10)*10**11*10**4<3*10**15)
    check("global_inverse_variation",R(6,5)+R(36,25)==R(264,100))
    check("compact_transport_variation",R(6,5)*20*45+56<10**4)
    check("finite_gluing_derivative",2*(2*10**5*R(6,5)+10**4)*400**2+3*10**6<10**12)
    check("local_massive_remainder",R(12,10)*10**16*81*ratio**3<R(424,100)*10**10)
    check("global_gluing_remainder",4*10**12*81*ratio<R(1141,1000)*10**12)
    check("common_restoration_remainder",10**12*81*ratio**2<R(1004,1000)*10**9)
    massive = R(12,10)*10**16*81*ratio**3+4*10**12*81*ratio+10**12*81*ratio**2+1
    check("complete_massive_value_ledger",massive<10**14)
    check("actual_domain_tail_commutators",2*(10*9)**2<2*10**4)
    # Direct exponential tail requires h<=10^-30 and M>=18. The log
    # exponent is minimized at those endpoints after retaining exp(-M/4).
    check("actual_domain_tail_log_margin",(R(9,20)*18-1)*30-R(9,10)*18-R(18,4)>180)
    # Both-pressure budget; the factor1/(8pi) is retained as pi>31/10.
    orbit = 10**14/(8*R(31,10))+10**5/(8*R(31,10))+2+2*10**4
    check("complete_orbit_pressure_budget",orbit<5*10**12)
    # j<3Ch,w<2C,|l3|<20C bounds generic derivative anomaly.
    anom = 3*(3+72+60)*ratio**2/(90*16*R(31,10)**2)
    check("generic_anomaly_pressure_budget",anom<1)
    check("complete_sphere_pressure_budget",orbit+anom/2<10**14)
    check("reference_pressure_derivatives",4*3300+4+70*3<14000)
    check("product_pressure_derivatives",4*3400+4*3400+2<30000)
    check("integrated_product_minus_reference",14000+30000<10**5)
    check("thermal_radius_dilution",R(1001,1000)**4<R(1005,1000))
    check("product_reference_R_units",10**5*R(1005,1000)<2*10**6)
    check("complete_Gamma_rest_pressure",10**14+2*10**6<R(1001,1000)*10**14)
    check("complete_relative_gap_forcing",4*R(22,7)*R(1001,1000)*10**14/k**8
          *(1+R(1,10**21))<R(51,10**13))
    check("physical_horizon_target",R(51,10**13)<R(1,100))
    check("small_temperature_root_bracket",1-6*R(1,1000)>0)
    check("thermal_K_derivative_guard",(5+18*1+19*R(1,10**12))/(2*(1+6*1-R(1,10**12)))<3)
    # |log(mu0 B)|<1e-21 cannot exhaust the retained Fang<100 margin;
    # the exact scale dependence is separately checked symbolically below.
    z,L = S.symbols("z L")
    check("same_scale_small_log",4*R(1,10**21)<R(1,10**18))
    zero("a0_scale_restoration",S.diff((1+z)*S.exp(2*L*z),z).subs(z,0)-(1+2*L))
    Fp = S.symbols("Fang")
    zero("a4_scale_restoration",S.limit((Fp+2/z)*S.exp(2*L*z)-2/z,z,0)-(Fp+4*L))


def hashes(data):
    return {name:hashlib.sha256(raw).hexdigest() for name,raw in data.items()}


def main():
    before = {name:(ROOT/name).read_bytes() for name in INPUTS}
    existing = json.loads(RECEIPT.read_text(encoding="utf-8")) if RECEIPT.exists() else {"runs":[]}
    run = {"started_utc":datetime.now(timezone.utc).isoformat(),"contract":CONTRACT,
           "command":sys.orig_argv,"dont_write_bytecode":sys.dont_write_bytecode,"cwd":str(Path.cwd()),
           "python":platform.python_version(),"sympy":S.__version__,"execution_inputs_sha256":hashes(before)}
    with LOG.open("a",encoding="utf-8",newline="\n") as log, contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
        print("RUN "+json.dumps(run,sort_keys=True),flush=True)
        try:
            check("actual_bytecode_disabled",sys.dont_write_bytecode and "-B" in sys.orig_argv)
            check("inputs_strict_UTF8_LF",all(raw.decode("utf-8") is not None and b"\r" not in raw
                  and not raw.startswith(b"\xef\xbb\xbf") for raw in before.values()))
            seed_identities()
            cap_identities()
            curvature_and_weyl()
            rational_margins()
        except Exception:
            traceback.print_exc()
            check("uncaught_exception",False)
        run["final_inputs_sha256"] = hashes({name:(ROOT/name).read_bytes() for name in INPUTS})
        print("FINAL_INPUTS "+json.dumps(run["final_inputs_sha256"],sort_keys=True),flush=True)
        check("frozen_inputs_unchanged",run["execution_inputs_sha256"]==run["final_inputs_sha256"])
        run["checks"] = ROWS
        run["failed_checks"] = [row for row in ROWS if not row["passed"]]
        code = int(bool(run["failed_checks"]))
        run["status"] = "failed" if code else "passed"
        run["terminal_exit_status"] = code
        run["finished_utc"] = datetime.now(timezone.utc).isoformat()
        if code:
            run["failed_inputs_base64"] = {name:base64.b64encode(raw).decode("ascii") for name,raw in before.items()}
            for name,raw in before.items():
                print(json.dumps({"path":name,"sha256":hashlib.sha256(raw).hexdigest(),
                                  "base64":base64.b64encode(raw).decode("ascii")},sort_keys=True))
        print("TERMINAL "+json.dumps({"status":run["status"],"checks":len(ROWS),
              "failures":len(run["failed_checks"]),"terminal_exit_status":code,
              "finished_utc":run["finished_utc"]},sort_keys=True),flush=True)
    existing["runs"].append(run)
    RECEIPT.write_text(json.dumps(existing,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps({"status":run["status"],"checks":len(ROWS),"failures":len(run["failed_checks"]),
                      "receipt":str(RECEIPT)},sort_keys=True))
    return code


if __name__ == "__main__":
    sys.exit(main())
