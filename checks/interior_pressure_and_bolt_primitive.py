"""S148: pressure response, moving proper cut and signed bolt primitive.

Manufactured radial sources test identities, not a quantum state. Tolerances,
steps and resolutions below are fixed before the first acceptance run.
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
import numpy as np
import sympy as s
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[1]
STEM = "interior-pressure-and-bolt-primitive"
INPUTS = [
    "checks/interior_pressure_and_bolt_primitive.py",
    "notes/interior-pressure-and-bolt-primitive.md",
    "notes/general-interior-dirac-response.md",
    "notes/outer-ward-and-constraint-transfer.md",
    "notes/stationary-centered-exterior-source.md",
    "notes/stationary-centered-complete-residual.md",
    "notes/stationary-centered-first-order-geometry.md",
    "notes/retained-benchmark-parameter-recovery.md",
]
ROWS, NUMERICS = [], []
FD_STEPS = [".001", ".0005"]
QUADRATURE_NODES = [32, 64]
FD_TOL = "2e-7"
QUAD_TOL = "1e-11"


def check(name, passed, detail=""):
    assert not any(row["name"] == name for row in ROWS)
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not row["passed"]:
        raise AssertionError(name)


def symbolic_checks():
    W, Wp, I, SS = s.symbols("W Wp I S", nonzero=True)
    P = I/(2*W)
    Pp = s.diff(P, I)*Wp*SS+s.diff(P, W)*Wp
    check("regular_pressure_ward_identity", s.simplify(Pp-Wp/W*(SS/2-P)) == 0)
    ss, tt, dh, dk, hp, kp, ell = s.symbols("ss tt dh dk hp kp ell")
    kernel = s.exp(2*(ss*dh+tt*dk))*(ell+2*ss*hp+2*tt*kp)/2
    zero = {ss: 0, tt: 0}
    check("first_ward_operator_kernel",
          s.simplify(s.diff(kernel, ss).subs(zero)-(dh*ell+hp)) == 0)
    check("mixed_ward_operator_kernel",
          s.simplify(s.diff(kernel, ss, tt).subs(zero)-
                     (2*dh*dk*ell+2*dh*kp+2*dk*hp)) == 0)
    # W(0)=0 makes F[constant]=constant/2 for every parameter.
    wc, const = s.symbols("wc const", nonzero=True)
    check("constant_trace_exact_cancellation", s.simplify(const*(wc-0)/(2*wc)-const/2) == 0)
    ch, ck, chk, E, Eh, Ek, Echi, Ih, Ik, Ihk = s.symbols(
        "ch ck chk E Eh Ek Echi Ih Ik Ihk", nonzero=True)
    cutfirst = -E*ch+Ih
    cutmixed = -E*chk-Eh*ck-Ek*ch-Echi*ch*ck+Ihk
    check("implicit_cut_first", s.simplify(cutfirst.subs(ch, Ih/E)) == 0)
    check("implicit_cut_mixed", s.simplify(cutmixed.subs(
        chk, Ihk/E-Eh*ck/E-Ek*ch/E-Echi*ch*ck/E)) == 0)
    # An arbitrary second jet avoids SymPy's inability to canonicalize
    # equivalent Subs(Derivative(...)) dummy variables in the first run.
    xx, c0 = s.symbols("xx c0")
    p0, ps, pt, px, pst, psx, ptx, pxx = s.symbols("p0 ps pt px pst psx ptx pxx")
    pf = p0+ps*ss+pt*tt+px*(xx-c0)+pst*ss*tt+psx*ss*(xx-c0)+ptx*tt*(xx-c0)+pxx*(xx-c0)**2/2
    cf = c0+ch*ss+ck*tt+chk*ss*tt
    composite = pf.subs(xx, cf)
    chain = (s.diff(pf, ss, tt)+s.diff(pf, ss, xx)*s.diff(cf, tt)+
             s.diff(pf, tt, xx)*s.diff(cf, ss)+
             s.diff(pf, xx, xx)*s.diff(cf, ss)*s.diff(cf, tt)+
             s.diff(pf, xx)*s.diff(cf, ss, tt)).subs(xx, cf)
    check("moving_cut_full_chain_rule",
          s.simplify((s.diff(composite, ss, tt)-chain).subs(zero)) == 0)
    B, C, r, rp, F, sig, PP = s.symbols("B C r rp F sig P", nonzero=True)
    source = B**2*C*rp*(2*PP-SS)/(2*r**3*F)
    HH = B**2*C*rp*W/(r**3*F*Wp)
    check("bolt_lapse_ward_sign", s.simplify(source+HH*Wp/W*(SS/2-PP)) == 0)
    check("bolt_coefficient_geometric_transfer", s.simplify(
        HH.subs(F, C*s.exp(-2*sig)*(rp/r)**2)-
        B**2*s.exp(2*sig)*W/(r*rp*Wp)) == 0)
    rh, r2, w2, chi = s.symbols("rh r2 w2 chi", positive=True)
    HHseries = B**2*s.exp(2*sig)*w2*chi**2/(
        (rh+r2*chi**2)*(2*r2*chi)*(2*w2*chi))
    H0 = B**2*s.exp(2*sig)/(4*rh*r2)
    check("regular_bolt_coefficient_limit", s.simplify(s.limit(HHseries, chi, 0)-H0) == 0)
    Fhprime = 4*C*s.exp(-2*sig)*r2/rh**2
    check("horizon_slope_coefficient", s.simplify(B**2*C/(Fhprime*rh**3)-H0) == 0)
    Hv, Hprime, Zv, Zprime = s.symbols("H Hp Z Zp")
    check("signed_primitive_integration_by_parts",
          s.simplify(-Hprime*Zv-Hv*Zprime+Hprime*Zv+Hv*Zprime) == 0)
    gap = s.symbols("gap", positive=True)
    fixtureH0 = 1/(2*gap*(1+gap))
    check("cold_coefficient_inverse_gap", s.limit(gap*fixtureH0, gap, 0) == s.Rational(1, 2))
    # Oscillatory counterexample: after integration by parts in (t+v)^2,
    # |P(t)| <= a(0)+integral |a'|=1/t^2, where a decreases to zero.
    t, v = s.symbols("t v", positive=True)
    amp = s.exp(-2*v)/(2*(t+v)**2)
    check("counterexample_pressure_ibp_amplitude",
          s.simplify(s.diff(amp, v)+s.exp(-2*v)*((t+v)**-2+(t+v)**-3)) == 0)
    check("counterexample_pressure_integrable_bound", s.simplify(2*amp.subs(v, 0)-1/t**2) == 0)
    n = s.symbols("n", positive=True)
    check("counterexample_absolute_trace_harmonic_lower_bound",
          s.integrate(s.sin(v), (v, 0, s.pi)) == 2 and
          s.simplify((s.Rational(1, 2)*2/((n+1)*s.pi))-1/((n+1)*s.pi)) == 0)


def quantitative_checks():
    mp.mp.dps = 60
    check("interior_radius_guard", 1/(1-mp.exp(-mp.mpf(".04")/mp.sqrt(8))) < 80)
    seed = 2*4*10**15*Q(68, 10**29)*80
    check("seed_normalized_trace_budget", seed == Q(4352, 10**13))
    K = 2*Q(85, 10**34)*(Q(2, 3)*10**17+7*10**11*24*256**3+10**17)
    check("complex_scalar_budget", K < Q(5, 10**12), str(K))
    check("first_cauchy_budget", 200*K < Q(1, 10**9))
    check("mixed_cauchy_budget", 40000*K < Q(2, 10**7))
    S0, S1, S2 = Q(5, 10**10), Q(1, 10**9), Q(2, 10**7)
    check("finite_background_trace_budget", seed+Q(1, 1000)*S1 < S0)
    check("ward_kernel_integral_guard", mp.exp(mp.mpf(".004"))*mp.mpf("1.01")*82 < 84)
    check("ward_total_variation_guard", mp.exp(mp.mpf(".004"))*(
        mp.mpf("1.02")+2*mp.mpf(".008")*mp.mpf("1.01")*82) < mp.mpf("2.4"))
    F0, F1, F2 = Q(12, 10), Q(700), Q(5500)
    check("first_operator_budget", 2*Q(24, 10)+8*84 < F1)
    check("mixed_operator_budget", 8*Q(24, 10)+64*84 < F2)
    P0, P1, P2 = Q(6, 10**10), Q(36, 10**8), Q(44, 10**7)
    check("pressure_value_budget", F0*S0 == P0)
    check("pressure_first_budget", F0*S1+F1*S0 < P1)
    check("pressure_mixed_budget", F0*S2+2*F1*S1+F2*S0 == Q(439, 10**8) < P2)
    emin = mp.mpf(".39")*mp.exp(-mp.mpf(".001"))
    cut1 = mp.mpf(".05")/emin
    check("cut_first_budget", cut1 < mp.mpf(".13"))
    check("cut_mixed_budget", cut1+2*mp.mpf(".13")+3*mp.mpf(".13")**2 < mp.mpf(".45"))
    check("proper_cut_inside_interior_cover", mp.mpf(".05")*mp.exp(-mp.mpf(".001")) > mp.mpf(".04"))
    seamY = mp.mpf(".05")*mp.exp(mp.mpf(".001"))
    check("seam_connection_guard", mp.mpf("7.31")*2*seamY+mp.coth(1)+mp.mpf(".008") < 3)
    check("seam_log_warp_first_guard", 2*mp.mpf("7.31")*2*seamY+mp.mpf(".016") < 2)
    check("seam_log_warp_second_guard", 2*mp.mpf("7.31")*(
        3*2*seamY+mp.mpf("7.31")*(1+(2*seamY)**2))+mp.mpf(".2") < 120)
    check("chi_second_chart_conversion_guard", Q(731, 100)**2+3*Q(731, 100) < 100)
    diff = S0/2+P0
    Pc, Phc, Pcc0 = Q(2, 10**9), Q(74, 10**8), Q(11, 10**8)
    check("pressure_spatial_first_budget", 2*diff < Pc)
    check("pressure_mixed_spatial_first_budget", 16*diff+2*(S1/2+P1) < Phc)
    check("pressure_spatial_second_conditional_budget", 120*diff+2*Pc < Pcc0)
    check("proper_cut_first_pressure_budget", P1+Pc*Q(13, 100) < Q(37, 10**8))
    base2 = P2+2*Phc*Q(13, 100)+Pcc0*Q(13, 100)**2+Pc*Q(45, 100)
    check("proper_cut_mixed_conditional_constant", base2 < Q(46, 10**7), str(base2))
    check("proper_cut_unresolved_spatial_coefficient", Q(13, 100)**2 < Q(17, 1000))
    check("bolt_pressure_subtraction_budgets",
          P0+S0/2 == Q(85, 10**11) and P1+S1/2 < Q(37, 10**8) and
          P2+S2/2 == Q(45, 10**7))
    check("mixed_primitive_product_rule_constants",
          2*Q(37, 10**8) == Q(74, 10**8))
    check("dini_response_operator_constants", 2*F0 == Q(24, 10) and 2*F1 == 1400 and 2*F2 == 11000)
    NUMERICS.append({"complex_scalar_budget": str(K), "moving_cut_mixed_constant": str(base2),
                     "spatial_trace_bound_J_supplied": False,
                     "cold_geometry_B0_B1_B2_supplied": False})


class Jet:
    """Value, ds, dt, dsdt, without factorial scaling."""
    def __init__(self, value, ds=0, dt=0, dst=0):
        self.v = [mp.mpf(v) for v in [value, ds, dt, dst]]
    @staticmethod
    def cast(x):
        return x if isinstance(x, Jet) else Jet(x)
    def __add__(self, other):
        return Jet(*(a+b for a,b in zip(self.v, self.cast(other).v)))
    __radd__ = __add__
    def __neg__(self):
        return Jet(*(-x for x in self.v))
    def __sub__(self, other):
        return self+-self.cast(other)
    def __rsub__(self, other):
        return self.cast(other)+-self
    def __mul__(self, other):
        a,b,c,d = self.v
        e,f,g,h = self.cast(other).v
        return Jet(a*e, b*e+a*f, c*e+a*g, d*e+b*g+c*f+a*h)
    __rmul__ = __mul__
    def reciprocal(self):
        a,b,c,d = self.v
        return Jet(1/a, -b/a**2, -c/a**2, 2*b*c/a**3-d/a**2)
    def __truediv__(self, other):
        return self*self.cast(other).reciprocal()
    def exp(self):
        a,b,c,d = self.v
        e = mp.exp(a)
        return Jet(e, e*b, e*c, e*(d+b*c))


def fixture(x):
    """Smooth even radial data, deliberately unrelated to a Dirac state."""
    b = mp.mpf(".003")*x*x
    bp = mp.mpf(".006")*x
    bpp = mp.mpf(".006")
    h, k = mp.mpf(".03")+x*x/100, -mp.mpf(".02")+mp.sin(x)**2/40
    hp, kp = x/50, mp.sin(x)*mp.cos(x)/20
    hpp, kpp = mp.mpf(".02"), mp.cos(2*x)/20
    trace = [1+x*x/5, mp.mpf(".3")+x*x/10, -mp.mpf(".2")+mp.mpf(".07")*x**4,
             mp.mpf(".11")+mp.mpf(".03")*x*x]
    return b,bp,bpp,h,k,hp,kp,hpp,kpp,trace


def manufactured_checks():
    mp.mp.dps = 55
    summaries = []
    for n in QUADRATURE_NODES:
        rawx, raww = leggauss(n)
        nodes = [(mp.mpf(float(x))+1)/2 for x in rawx]
        weights = [mp.mpf(float(w))/2 for w in raww]
        def quad(fun, lo, hi):
            return sum((hi-lo)*wt*fun(lo+(hi-lo)*xx) for xx,wt in zip(nodes,weights))
        def jet_fields(x, freeze=False):
            b,bp,bpp,h,k,hp,kp,hpp,kpp,T = fixture(x)
            if freeze:
                h=k=hp=kp=hpp=kpp=0
            return Jet(b,h,k), Jet(bp,hp,kp), Jet(bpp,hpp,kpp), Jet(*T)
        def scalar_fields(x, ss=0, tt=0):
            b,bp,bpp,h,k,hp,kp,hpp,kpp,T = fixture(x)
            return b+ss*h+tt*k, bp+ss*hp+tt*kp, T[0]+ss*T[1]+tt*T[2]+ss*tt*T[3]
        def pressure_jet(c, freeze=False):
            b,_,_,_ = jet_fields(c, freeze)
            denom = (2*b).exp()*mp.sinh(c)**2
            def integrand(x):
                b,bp,_,T = jet_fields(x, freeze)
                return (2*b).exp()*mp.sinh(x)*(mp.cosh(x)+bp*mp.sinh(x))*T/denom
            return quad(integrand, mp.mpf(0), c)
        def pressure_scalar(c, ss=0, tt=0):
            bc,_,_ = scalar_fields(c,ss,tt)
            denom = mp.exp(2*bc)*mp.sinh(c)**2
            def integrand(x):
                b,bp,T = scalar_fields(x,ss,tt)
                return mp.exp(2*b)*mp.sinh(x)*(mp.cosh(x)+bp*mp.sinh(x))*T/denom
            return quad(integrand, mp.mpf(0), c)
        def cut(ss=0,tt=0):
            # Independent adaptive quadrature and nonlinear root solve.
            return mp.findroot(lambda c: mp.quad(
                lambda x: mp.exp(scalar_fields(x,ss,tt)[0]), [c,2])-mp.mpf(".05"), (mp.mpf("1.94"),mp.mpf("1.97")))
        c = cut()
        b,bp,bpp,h,k,hp,kp,hpp,kpp,T = fixture(c)
        Ec = mp.exp(b)
        ch = quad(lambda x: mp.exp(fixture(x)[0])*fixture(x)[3],c,mp.mpf(2))/Ec
        ck = quad(lambda x: mp.exp(fixture(x)[0])*fixture(x)[4],c,mp.mpf(2))/Ec
        chk = quad(lambda x: mp.exp(fixture(x)[0])*fixture(x)[3]*fixture(x)[4],c,mp.mpf(2))/Ec-h*ck-k*ch-bp*ch*ck
        pj = pressure_jet(c)
        pv,ph,pk,phk = pj.v
        ell = 2*mp.coth(c)+2*bp
        ellp = -2/mp.sinh(c)**2+2*bpp
        pc = ell*(T[0]/2-pv)
        phc = 2*hp*(T[0]/2-pv)+ell*(T[1]/2-ph)
        pkc = 2*kp*(T[0]/2-pv)+ell*(T[2]/2-pk)
        pcc = ellp*(T[0]/2-pv)+ell*(c/5-pc)
        moving = phk+phc*ck+pkc*ch+pcc*ch*ck+pc*chk
        bc = scalar_fields(c)[0]
        independent = mp.quad(lambda x: mp.exp(2*(scalar_fields(x)[0]-bc))*mp.sinh(x)*(
            mp.cosh(x)+scalar_fields(x)[1]*mp.sinh(x))*scalar_fields(x)[2]/mp.sinh(c)**2, [0,c/2,c])
        check(f"independent_pressure_quadrature_{n}", abs(independent-pv) < mp.mpf(QUAD_TOL))
        check(f"nonzero_hh_pressure_retained_{n}", pv > mp.mpf(".5"))
        check(f"kernel_variation_omission_detected_{n}", abs(phk-pressure_jet(c,True).v[3]) > mp.mpf("1e-7"))
        check(f"cut_transport_omission_detected_{n}", abs(moving-phk) > mp.mpf("1e-7"))
        fixed_errors, moving_errors, cut_errors = [],[],[]
        for step_text in FD_STEPS:
            step = mp.mpf(step_text)
            fvals,mvals,cvals = [],[],[]
            for ss,tt in [(step,step),(step,-step),(-step,step),(-step,-step)]:
                cb = cut(ss,tt)
                fvals.append(pressure_scalar(c,ss,tt))
                mvals.append(pressure_scalar(cb,ss,tt))
                cvals.append(cb)
            mixed = lambda vals: (vals[0]-vals[1]-vals[2]+vals[3])/(4*step**2)
            fixed_errors.append(abs(mixed(fvals)-phk))
            moving_errors.append(abs(mixed(mvals)-moving))
            cut_errors.append(abs(mixed(cvals)-chk))
        for label, errors, target in [("fixed_pressure",fixed_errors,phk),
                                      ("moving_pressure",moving_errors,moving),
                                      ("proper_cut",cut_errors,chk)]:
            check(f"{label}_independent_mixed_difference_{n}",
                  max(errors) < mp.mpf(FD_TOL)*(1+abs(target)),
                  {"target":mp.nstr(target,25),"errors":[mp.nstr(e,15) for e in errors]})
            check(f"{label}_step_refinement_{n}",
                  errors[1] < mp.mpf(".35")*errors[0] or max(errors) < mp.mpf("1e-35"))
        # Constant traces test the integrated kernel and its responses.
        def constant_kernel(x):
            bx,bpx,_,_ = jet_fields(x)
            bcj,_,_,_ = jet_fields(c)
            return (2*(bx-bcj)).exp()*mp.sinh(x)*(mp.cosh(x)+bpx*mp.sinh(x))/mp.sinh(c)**2
        constj = quad(constant_kernel, mp.mpf(0),c)
        check(f"constant_trace_operator_annihilation_{n}",
              abs(constj.v[0]-mp.mpf(".5")) < mp.mpf(QUAD_TOL) and
              max(abs(v) for v in constj.v[1:]) < mp.mpf(QUAD_TOL))
        # Bolt fixture r=1+gap*cosh(chi), B=C=1; geometry response is in b.
        gap, end = mp.mpf(".2"), mp.mpf(".7")
        def bolt_coefficients(x):
            bj,bpj,bppj,Tj = jet_fields(x)
            rr, rp, rpp = 1+gap*mp.cosh(x), gap*mp.sinh(x), gap*mp.cosh(x)
            ej = 2*mp.coth(x)+2*bpj
            epj = -2/mp.sinh(x)**2+2*bppj
            Hj = (2*bj).exp()/(rr*rp*ej)
            Hpj = Hj*(2*bpj-rp/rr-rpp/rp-epj/ej)
            qj = (2*bj).exp()/(2*rr*rp)
            return Hj,Hpj,qj,Tj
        Pbolt = Jet(*fixture(mp.mpf(0))[-1])/2
        direct = quad(lambda x: bolt_coefficients(x)[2]*(2*pressure_jet(x)-bolt_coefficients(x)[3]),mp.mpf(0),end)
        primitive = -bolt_coefficients(end)[0]*(pressure_jet(end)-Pbolt)+quad(
            lambda x: bolt_coefficients(x)[1]*(pressure_jet(x)-Pbolt),mp.mpf(0),end)
        check(f"signed_bolt_primitive_all_parameter_components_{n}",
              max(abs(a-b) for a,b in zip(direct.v,primitive.v)) < mp.mpf(QUAD_TOL),
              [mp.nstr(v,25) for v in direct.v])
        def scalar_lapse(ss,tt):
            def integrand(x):
                bx,_,Tx = scalar_fields(x,ss,tt)
                rr,rp = 1+gap*mp.cosh(x),gap*mp.sinh(x)
                return mp.exp(2*bx)/(2*rr*rp)*(2*pressure_scalar(x,ss,tt)-Tx)
            return quad(integrand,mp.mpf(0),end)
        lapse_errors = []
        for step_text in FD_STEPS:
            step = mp.mpf(step_text)
            vals = [scalar_lapse(ss,tt) for ss,tt in [(step,step),(step,-step),(-step,step),(-step,-step)]]
            observed = (vals[0]-vals[1]-vals[2]+vals[3])/(4*step**2)
            lapse_errors.append(abs(observed-direct.v[3]))
        check(f"bolt_lapse_independent_mixed_difference_{n}",
              max(lapse_errors) < mp.mpf(FD_TOL)*(1+abs(direct.v[3])),
              [mp.nstr(v,15) for v in lapse_errors])
        check(f"bolt_lapse_step_refinement_{n}",
              lapse_errors[1] < mp.mpf(".35")*lapse_errors[0] or max(lapse_errors) < mp.mpf("1e-35"))
        summaries.append([phk,moving,chk,direct.v[3]])
        NUMERICS.append({"nodes":n,"fixed_pressure_mixed":mp.nstr(phk,35),
                         "moving_pressure_mixed":mp.nstr(moving,35),"proper_cut_mixed":mp.nstr(chk,35),
                         "bolt_lapse_mixed":mp.nstr(direct.v[3],35),
                         "cut_transport_contribution":mp.nstr(moving-phk,35)})
    check("all_mixed_responses_quadrature_refinement",
          max(abs(a-b) for a,b in zip(*summaries)) < mp.mpf(QUAD_TOL))


def hashes():
    return {p:sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS}


def main():
    before = hashes()
    error = None
    try:
        symbolic_checks()
        quantitative_checks()
        manufactured_checks()
    except Exception:
        error = traceback.format_exc()
        print(error,file=sys.stderr)
    after = hashes()
    ok = error is None and before == after and all(row["passed"] for row in ROWS)
    record = {
        "story":"S148","created_utc":datetime.now(timezone.utc).isoformat(),
        "status":"passed" if ok else "failed","python":sys.version,
        "numpy":np.__version__,"sympy":s.__version__,"mpmath":mp.__version__,
        "input_sha256_before":before,"input_sha256_after":after,"unchanged_inputs":before==after,
        "protocol":{"finite_difference_steps":FD_STEPS,"quadrature_nodes":QUADRATURE_NODES,
                    "finite_difference_tolerance":FD_TOL,"quadrature_tolerance":QUAD_TOL},
        "checks":ROWS,"numerics":NUMERICS,"exception":error,
        "scope":"fixed-cut canonical pressure response, conditional moving-cut Hessian and signed bolt primitive; manufactured sources are diagnostics",
        "unresolved":["spatial trace bound at moving proper cut","cold geometric/source error control",
                      "other local and radius source maps","global exterior/state solve and physical gates"],
    }
    if error:
        record["failed_checker_source_utf8"] = Path(__file__).read_text(encoding="utf-8")
        record["failed_note_source_utf8"] = (ROOT/INPUTS[1]).read_text(encoding="utf-8")
    output = ROOT/"evidence"/f"{STEM}.json"
    history = json.loads(output.read_text(encoding="utf-8")) if output.exists() else {"runs":[]}
    history["runs"].append(record)
    output.write_text(json.dumps(history,indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps({"status":record["status"],"checks":len(ROWS),"unchanged_inputs":before==after}),flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
