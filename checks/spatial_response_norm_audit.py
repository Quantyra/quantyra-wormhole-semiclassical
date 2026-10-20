"""S149: exact flat regularity endpoint and sufficient Fourier response norm.

Finite diagnostics do not prove a curved mapping or a physical solution.
The infinite-frequency conclusion is the harmonic lower-bound argument.
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
import numpy as np
import sympy as s
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[1]
STEM = "spatial-response-norm-audit"
INPUTS = [
    "checks/spatial_response_norm_audit.py",
    "notes/spatial-response-norm-audit.md",
    "notes/flat-variable-dirac-three-point.md",
    "checks/flat_variable_dirac_three_point.py",
    "notes/stationary-curved-quantum-tt.md",
    "notes/general-interior-dirac-response.md",
    "notes/interior-pressure-and-bolt-primitive.md",
    "notes/stationary-centered-canonical-response.md",
    "notes/stationary-centered-complete-residual.md",
]
ROWS, NUMERICS = [], []
PROTOCOL = {"gauss_nodes":[48,96],"tower_levels":[512,1024],
            "quadrature_tolerance":"1e-10","finite_difference_steps":[".0001",".00005"],
            "finite_difference_tolerance":"2e-7","refinement_ratio":".35"}


def check(name, passed, detail=""):
    assert not any(row["name"] == name for row in ROWS)
    row = {"name":name,"passed":bool(passed),"detail":str(detail)}
    ROWS.append(row)
    print(json.dumps(row),flush=True)
    if not row["passed"]:
        raise AssertionError(name)


def exact_checks():
    v,t,z,n,q,C,M,k,cc = s.symbols("v t z n q C M k cc", positive=True)
    R = s.log(1+v)-v+v**2/2
    check("positive_log_remainder_derivative", s.simplify(s.diff(R,v)-v**2/(1+v)) == 0)
    check("quadratic_remainder_upper_identity", s.simplify(v**2/2-R-(v-s.log(1+v))) == 0)
    check("cubic_remainder_upper_derivative", s.simplify(s.diff(v**3/3-R,v)-v**3/(1+v)) == 0)
    a=z*(1-z)
    check("second_and_third_parameter_moments",
          s.integrate(a**2,(z,0,1))==s.Rational(1,30) and s.integrate(a**3,(z,0,1))==s.Rational(1,140))
    check("central_parameter_lower_interval", (s.Rational(1,4)*(1-s.Rational(1,4)))==s.Rational(3,16))
    check("large_argument_quadratic_lower", s.simplify(v**2/2-v-v**2/6-v*(v-3)/3) == 0)
    check("mode_lower_coefficient", Q(1,2)*Q(1,2)*Q(1,6)*Q(3,16)**2==Q(3,2048))
    Mn2=n*(n+q)/C
    dn=2*(q+2*n)
    check("inverse_square_partial_fraction", s.simplify(dn/Mn2-2*C*(1/n+1/(n+q))) == 0)
    check("inverse_fourth_partial_fraction", s.simplify(dn/Mn2**2-2*C**2/q*(1/n**2-1/(n+q)**2)) == 0)
    check("counting_function_bound", s.simplify(4*n*(n+q)-2*n*(q+n+1)-2*n*(n+q-1)) == 0)
    check("split_tower_constant", Q(4,120)+Q(8,840)==Q(3,70))
    check("normalized_harmonic_lower_coefficient", 2*Q(3,2048)*2==Q(3,512))
    check("contact_reserve", Q(100,60*3)==Q(5,9))
    check("high_frequency_linear_upper_budget", Q(5,9)+2*Q(3,70)/3 < Q(6,10) and Q(2,45)<Q(5,100))
    check("low_frequency_linear_upper_budget", Q(5,9)+Q(16,840*3) < Q(6,10))
    check("third_contact_normalization", Q(32*100*3,960*3)==Q(10,3))
    # The actual simpler estimate: 10/3 < 5 log 2, using log 2 > 2/3.
    check("third_contact_log_reserve", Q(10,3)==5*Q(2,3))
    check("mixed_scalar_multiplier_reserve", 20+5+2+2 < 30)
    check("output_gradient_frequency_budget", 30*2*16*2==1920 < 2000)
    x,ss,tt=s.symbols("x ss tt")
    f0,fs,ft,fst,h0,w0=s.symbols("f0 fs ft fst h0 w0")
    scalar=s.exp(-2*(ss*h0+tt*w0))*(f0+ss*fs+tt*ft+ss*tt*fst)
    check("scalar_density_mixed_conversion", s.simplify(s.diff(scalar,ss,tt).subs({ss:0,tt:0})-
          (fst-2*h0*ft-2*w0*fs+4*h0*w0*f0)) == 0)
    vac=s.symbols("vac")
    check("vacuum_scalar_response_cancellation",
          s.simplify((fst-2*h0*ft-2*w0*fs+4*h0*w0*f0).subs(
              {f0:vac,fs:2*h0*vac,ft:2*w0*vac,fst:4*h0*w0*vac})) == 0)
    Pmode=M**2/(2*s.pi)*R.subs(v,a*k**2/M**2)
    Vconst=M**2/s.pi*(s.log(1+a*k**2/M**2)-(a*k**2/M**2)/(1+a*k**2/M**2)-
                                (a*k**2/M**2)**2/2)
    check("constant_action_vertex_identity", s.simplify(Vconst-2*Pmode+k*s.diff(Pmode,k)) == 0)
    check("common_contact_constant_identity", s.simplify(-2*cc*k**4-2*cc*k**4+k*s.diff(cc*k**4,k)) == 0)
    # Every derivative of k^-5 cos(k x), through five, has amplitude k^(j-5).
    check("c5_sequence_derivative_amplitudes",
          all(s.simplify(s.diff(k**-5*s.exp(s.I*k*x),x,j)/(s.I**j*s.exp(s.I*k*x))-k**(j-5)) == 0
              for j in range(6)))
    check("source_gradient_sequence_cancellation", s.simplify(k*k**-5*k**4-1) == 0)
    check("source_value_sequence_decay", s.limit(s.log(2+k)/k,k,s.oo) == 0)
    check("harmonic_lower_bound_diverges", s.limit(s.log(n+1),n,s.oo) == s.oo)
    check("extra_power_single_wave_decay", s.limit(s.log(2+k)/s.sqrt(k),k,s.oo) == 0)


def Rvalue(v):
    return mp.log1p(v)-v+v*v/2


def closed_parameter(t):
    if t == 0:
        return mp.mpf(0)
    y=mp.sqrt(t/(t+4))
    return 2*mp.atanh(y)/y-2-t/6+t*t/60


def mode(k,M=1):
    M=mp.mpf(M)
    return M*M*closed_parameter(k*k/(M*M))/(2*mp.pi)


def diagnostic_checks():
    mp.mp.dps=60
    check("logarithm_budget_guards", mp.log(2)>mp.mpf(2)/3 and mp.log(18)>2)
    for index,(momentum,mass) in enumerate([(".2","1"),("1","1"),("4","1"),("19","3")]):
        k,M=mp.mpf(momentum),mp.mpf(mass)
        target=mode(k,M)
        adaptive=M*M/(2*mp.pi)*mp.quad(lambda z:Rvalue(z*(1-z)*k*k/(M*M)),[0,.5,1])
        check(f"closed_parameter_vs_adaptive_{index}", abs(target-adaptive)<mp.mpf("1e-45")*(1+abs(target)))
        for nodes in PROTOCOL["gauss_nodes"]:
            xx,ww=leggauss(nodes)
            observed=M*M/(2*mp.pi)*sum(mp.mpf(float(w))/2*Rvalue(
                ((mp.mpf(float(x))+1)/2)*((1-mp.mpf(float(x)))/2)*k*k/(M*M)) for x,w in zip(xx,ww))
            check(f"independent_gauss_parameter_{index}_{nodes}",
                  abs(observed-target)<mp.mpf(PROTOCOL["quadrature_tolerance"])*(1+abs(target)))
        lower=3*k**4/(2048*mp.pi*M*M)
        check(f"mode_envelopes_{index}", 0<=target<=min(k**4/(120*mp.pi*M*M),k**6/(840*mp.pi*M**4)) and
              (k<4*M or target>=lower))
        contact=mp.mpf(".13")
        K=lambda p:mode(p,M)+contact*p**4
        vertex=M*M/mp.pi*mp.quad(lambda z:mp.log1p(z*(1-z)*k*k/M**2)-
              (z*(1-z)*k*k/M**2)/(1+z*(1-z)*k*k/M**2),[0,.5,1])-k**4/(60*mp.pi*M*M)-2*contact*k**4
        mixed=vertex-2*K(k)
        target_mixed=-k*mp.diff(K,k)
        check(f"independent_constant_scalar_identity_{index}",
              abs(mixed-target_mixed)<mp.mpf("1e-45")*(1+abs(target_mixed)))
        errors=[]
        for step_text in PROTOCOL["finite_difference_steps"]:
            step=mp.mpf(step_text)
            observed=(K(k*mp.exp(-step))-K(k*mp.exp(step)))/(2*step)
            errors.append(abs(observed-mixed))
        check(f"constant_scale_finite_difference_{index}",
              max(errors)<mp.mpf(PROTOCOL["finite_difference_tolerance"])*(1+abs(mixed)),
              {"mixed":mp.nstr(mixed,28),"errors":[mp.nstr(e,18) for e in errors]})
        check(f"constant_scale_step_refinement_{index}",
              errors[1]<mp.mpf(PROTOCOL["refinement_ratio"])*errors[0] or max(errors)<mp.mpf("1e-35"))
        check(f"density_omission_detected_{index}", abs(vertex-mixed)>mp.mpf("1e-5"))
    # Synthetic parameters test arithmetic only, without replacing the benchmark.
    q,C,F=7,mp.mpf(2),mp.mpf(-11)
    def tower(k,N):
        return -F*k**4/(60*mp.pi)+2/C*sum(2*(q+2*n)*mode(k,mp.sqrt(n*(n+q)/C)) for n in range(1,N+1))
    def tail(k,N):
        inv4=2*C*C/q*sum(mp.mpf(1)/j**2 for j in range(N+1,N+q+1))
        return 2/C*k**6/(840*mp.pi)*inv4
    for ktext in ["1","7","25","80"]:
        k=mp.mpf(ktext)
        vals=[tower(k,N) for N in PROTOCOL["tower_levels"]]
        tails=[tail(k,N) for N in PROTOCOL["tower_levels"]]
        check(f"positive_tower_tail_{ktext}", 0<vals[1]-vals[0]<tails[0] and tails[1]<tails[0])
        check(f"nested_tower_enclosures_{ktext}", vals[0]+tails[0]>vals[1]+tails[1])
        # For this fixture M1=2, use the unreduced counting upper bound.
        M1=mp.sqrt((1+q)/C)
        if k>=M1:
            bound=k**4*(abs(F)/(60*mp.pi)+2/mp.pi*(mp.mpf(3)/70+mp.log(k/M1)/15))
        else:
            bound=abs(F)*k**4/(60*mp.pi)+16*k**6/(840*mp.pi*M1*M1)
        check(f"full_tower_enclosure_vs_upper_{ktext}", vals[1]+tails[1]<bound)
        NUMERICS.append({"fixture":"q=7,C=2,F=-11; eta factored out","k":ktext,
                         "partial_sums":[mp.nstr(x,30) for x in vals],
                         "rigorous_tail_upper":[mp.nstr(x,25) for x in tails]})
    seq=[]
    for power in [4,8,16,32,64]:
        N=2**power
        H=mp.digamma(N+1)+mp.euler
        lower=-F/(60*mp.pi)+3*H/(512*mp.pi)
        seq.append(lower)
        check(f"harmonic_lower_comparison_{power}", H>=mp.log(N+1))
    check("increasing_gradient_lower_sequence", all(b>a for a,b in zip(seq,seq[1:])))
    NUMERICS.append({"gradient_lower_eta_factored_out":[mp.nstr(x,30) for x in seq],
                     "infinite_growth_proof":"harmonic lower bound, not finite extrapolation"})
    for index,(p1,p2) in enumerate([((.7,.2),(-.1,.8)),((3.,.5),(-1.,2.)),((0.,0.),(2.,0.))]):
        p1,p2=np.array(p1),np.array(p2)
        mags=[np.linalg.norm(p1),np.linalg.norm(p2),np.linalg.norm(p1+p2)]
        squares=[v*v for v in mags]
        contact=float(F)/(30*np.pi)*(squares[0]*squares[1]+squares[0]*squares[2]+squares[1]*squares[2])
        vertices=[]
        for nodes in PROTOCOL["gauss_nodes"]:
            rem=0.
            for n in range(1,65):
                Mn=np.sqrt(n*(n+q)/float(C))
                rem+=2*(q+2*n)*Mn*Mn*parameter_remainder(p1/Mn,p2/Mn,nodes)
            vertices.append(contact+2/float(C)*rem)
        K1,K2=float(tower(mp.mpf(float(mags[0])),64)),float(tower(mp.mpf(float(mags[1])),64))
        Bmix=vertices[1]-2*K1-2*K2
        check(f"mode_pair_vertex_quadrature_{index}",
              abs(vertices[1]-vertices[0])<1e-10*(1+abs(vertices[1])))
        check(f"mode_pair_scalar_density_omission_{index}", abs(Bmix-vertices[1])>1e-5)
        Kmax=max(mags)
        inv4=2*float(C)**2/q*sum(1/j**2 for j in range(65,65+q))
        a6=3037/(7560*np.pi)
        Btail=2/float(C)*a6*Kmax**6*inv4+2*float(tail(mp.mpf(float(mags[0])),64))+2*float(tail(mp.mpf(float(mags[1])),64))
        W=lambda p:(1+p)**5*np.log(2+p)
        check(f"mode_pair_gradient_norm_with_tail_{index}",
              mags[2]*(abs(Bmix)+Btail)<2000*(W(mags[0])+W(mags[1])))
        NUMERICS.append({"mode_pair":index,"density_vertex":vertices[1],"scalar_mixed":Bmix,
                         "retained_tail_upper":Btail})


def hashes():
    return {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    global parameter_remainder
    before=hashes()
    error=None
    try:
        # The pinned Windows Python uses an isolated module path.
        spec=importlib.util.spec_from_file_location("s149_flat_vertex_input",
                    ROOT/"checks"/"flat_variable_dirac_three_point.py")
        dependency=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(dependency)
        parameter_remainder=dependency.parameter_remainder
        exact_checks()
        diagnostic_checks()
    except Exception:
        error=traceback.format_exc()
        print(error,file=sys.stderr)
    after=hashes()
    ok=error is None and before==after and all(row["passed"] for row in ROWS)
    record={"story":"S149","created_utc":datetime.now(timezone.utc).isoformat(),
            "status":"passed" if ok else "failed","python":sys.version,"numpy":np.__version__,
            "sympy":s.__version__,"mpmath":mp.__version__,"protocol":PROTOCOL,
            "input_sha256_before":before,"input_sha256_after":after,"unchanged_inputs":before==after,
            "checks":ROWS,"numerics":NUMERICS,"exception":error,
            "scope":"flat comparison spatial regularity; no curved finite-background mapping or wormhole solution",
            "unresolved":["actual curved spatial response and norm transfer","cold weighted mixed response",
                          "other source maps and global solve","all physical completion gates"]}
    if error:
        record["failed_checker_source_utf8"]=Path(__file__).read_text(encoding="utf-8")
        record["failed_note_source_utf8"]=(ROOT/INPUTS[1]).read_text(encoding="utf-8")
    output=ROOT/"evidence"/f"{STEM}.json"
    history=json.loads(output.read_text(encoding="utf-8")) if output.exists() else {"runs":[]}
    history["runs"].append(record)
    output.write_text(json.dumps(history,indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps({"status":record["status"],"checks":len(ROWS),"unchanged_inputs":before==after}),flush=True)
    return 0 if ok else 1


if __name__=="__main__":
    raise SystemExit(main())
