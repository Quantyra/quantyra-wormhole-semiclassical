"""S147: exact proper-gauge/constraint transfer and finite-geometry Ward maps.

Manufactured-source quadratures are diagnostics, not a quantum state or
continuum verification of the inherited Dirac response estimates.
"""
from __future__ import annotations

from datetime import datetime, timezone
from fractions import Fraction as Q
from hashlib import sha256
from pathlib import Path
import json
import math
import sys
import traceback

import mpmath as mp
import numpy as np
import sympy as s
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[1]
STEM = "outer-ward-and-constraint-transfer"
INPUTS = [
    "checks/outer_ward_and_constraint_transfer.py",
    "notes/outer-ward-and-constraint-transfer.md",
    "notes/outer-finite-background-response.md",
    "notes/stationary-centered-conormal-source.md",
    "notes/stationary-centered-covariant-seed.md",
    "notes/stationary-centered-first-order-geometry.md",
    "notes/stationary-centered-complete-residual.md",
    "notes/retained-benchmark-parameter-recovery.md",
]
ROWS = []
NUMERICS = []


def check(name, passed, detail=""):
    assert not any(row["name"] == name for row in ROWS)
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not row["passed"]:
        raise AssertionError(name)


def symbolic_checks():
    u, up, upp, a, ap, B, C, G, x, rho, p = s.symbols(
        "u up upp a ap B C G x rho p", nonzero=True)
    F0, F0r = s.symbols("F0 F0r")
    r = B/u
    rY = -B*up/u**2
    FF = C*up**2/u**2
    expalpha = s.sqrt(x/C)*a/up
    check("exact_radial_metric_factor", s.simplify(rY**2/FF-B**2/(C*u**2)) == 0)
    check("exact_time_metric_factor", s.simplify(expalpha**2*FF-x*a**2/u**2) == 0)
    mass = r*(F0-FF)/2
    massY = s.diff(mass, u)*up+s.diff(mass, up)*upp+s.diff(mass, F0)*F0r*rY
    mass_equation = 2*C*u*upp-3*C*up**2+u**2*(F0+r*F0r)-2*G*rho*u**2
    check("exact_mass_constraint",
          s.simplify(massY-G*rho*rY+B*up*mass_equation/(2*u**4)) == 0)
    lapseY = ap/a-upp/up
    lapse_equation = C*(upp-ap*up/a)-G*u*(rho+p)
    check("exact_lapse_constraint",
          s.simplify(lapseY-G*(rho+p)*rY/(r*FF)+lapse_equation/(C*up)) == 0)
    pressure_equation = 3*C*up**2-2*C*u*ap*up/a-u**2*(F0+r*F0r)-2*G*p*u**2
    check("pressure_constraint_combination",
          s.simplify(mass_equation-2*u*lapse_equation+pressure_equation) == 0)
    P, PY, TT, AA = s.symbols("P PY T Abar")
    rhoph = u**2*(P-TT)/B**2
    pph = u**2*P/B**2
    psph = u**2*(AA-TT)/(2*B**2)
    ward = (2*u*up*P+u**2*PY)/B**2+(rhoph+pph)*(ap/a-up/u)+2*psph*up/u
    check("full_radius_and_trace_ward",
          s.simplify(B**2*ward/u**2-(PY+2*ap*P/a-ap*TT/a+up*AA/u)) == 0)
    check("proper_mass_integrand",
          s.simplify(-G*rhoph*rY/(B*x)-G*up*(P-TT)/(B**2*x)) == 0)
    check("proper_lapse_integrand",
          s.simplify(G*(rhoph+pph)*rY/(r*FF)+G*u**3*(2*P-TT)/(C*B**2*up)) == 0)

    Y, d, a2, cc, e3, e4, m = s.symbols("Y d a2 c e3 e4 m", nonzero=True)
    uu = d*Y+e3*Y**3+e4*Y**4
    aa = 1+a2*Y**2+cc*Y**3
    al = s.series(s.log(d*aa/s.diff(uu, Y)), Y, 0, 4).removeO().expand()
    check("lapse_quadratic_and_cubic_jets",
          s.simplify(al.coeff(Y, 2)-(a2-3*e3/d)) == 0 and
          s.simplify(al.coeff(Y, 3)-(cc-4*e4/d)) == 0)
    residualF = s.series((C*s.diff(uu, Y)**2-C*d*d)/uu**2, Y, 0, 2).removeO().expand()
    check("geometric_mass_outer_coefficients",
          s.simplify(residualF.coeff(Y, 0)-6*C*e3/d) == 0 and
          s.simplify(residualF.coeff(Y, 1)-8*C*e4/d) == 0)
    correlation = {a2: 1/(2*C), e3: d/(6*C), e4: d*cc/4}
    check("correlated_lapse_jets_vanish", s.simplify(al.subs(correlation)) == 0)
    check("correlated_physical_mass",
          s.simplify((residualF-1+2*m*d*Y/B).subs(correlation).subs(m, -B*C*cc/d)) == 0)
    pp = s.symbols("p_mass")
    check("mass_parameter_normalization",
          s.simplify(-B*C*(-d**3*pp)/d-B*(C*d*d)*pp) == 0)
    check("correlated_radius_quartic", s.simplify(d*(-d**3*pp)/4+d**4*pp/4) == 0)

    hs, ht, hY, kY, w0, hp, kp, ss, tt = s.symbols("hs ht hY kY w0 hp kp ss tt")
    kernel = s.exp(2*(ss*(hY-hs)+tt*(kY-ht)))
    factor = kernel*(w0-ss*hp-tt*kp)
    mixed = s.diff(factor, ss, tt).subs({ss: 0, tt: 0})
    expected = 4*(hY-hs)*(kY-ht)*w0-2*(hY-hs)*kp-2*(kY-ht)*hp
    check("exact_mixed_ward_kernel", s.simplify(mixed-expected) == 0)
    v, vp, w, wp = s.symbols("v vp w wp")
    radial = (u+ss*v+tt*w)**3/(up+ss*vp+tt*wp)
    mixed_radial = s.diff(radial, ss, tt).subs({ss: 0, tt: 0})/(u**3/up)
    expected_radial = 6*v*w/u**2-3*(v*wp+w*vp)/(u*up)+2*vp*wp/up**2
    check("exact_mixed_radius_prefactor", s.simplify(mixed_radial-expected_radial) == 0)
    HH, Tconst, Pmatch = s.symbols("H Tconst Pmatch")
    const_pressure = HH*Pmatch+(1-HH)*Tconst/2
    check("constant_trace_state_memory_cancellation",
          s.simplify(const_pressure.subs(Pmatch, Tconst/2)-Tconst/2) == 0)
    check("omitted_memory_changes_constant_trace",
          s.simplify(const_pressure-const_pressure.subs(Pmatch, 0)) == HH*Pmatch)


def quantitative_checks():
    mp.mp.dps = 65
    y0 = Q(1, 20)
    # These rational log endpoints enclose log(20); no unchecked decimal
    # integral is used in the algebraic budgets below.
    loglo, loghi = Q(2995732273, 10**9), Q(2995732274, 10**9)
    check("logarithm_rational_enclosure", mp.mpf(loglo.numerator)/loglo.denominator < mp.log(20) <
          mp.mpf(loghi.numerator)/loghi.denominator)
    def upper_I(j):
        n = j+1
        return y0**n/Q(n)*(1+loghi**2+2*loghi/n+Q(2, n*n))
    I0, I1, I3 = [upper_I(j) for j in [0, 1, 3]]
    J3 = y0**4/4
    check("outer_integral_enclosures",
          I0 < Q(899, 1000) and I1 < Q(1684, 10**5) and I3 < Q(1813, 10**8))
    for j in [0, 1, 3]:
        # Y=y0*z^6 makes the logarithmic endpoint easy to integrate.
        numeric = mp.quad(lambda z: 6*mp.mpf(".05")**(j+1)*z**(6*j+5)*(
            1+(mp.log(mp.mpf(".05"))+6*mp.log(z))**2), [0, .25, 1])
        exact = mp.mpf(".05")**(j+1)/(j+1)*(1+mp.log(mp.mpf(".05"))**2
             -2*mp.log(mp.mpf(".05"))/(j+1)+mp.mpf(2)/(j+1)**2)
        check("independent_weight_integral_"+str(j), abs(numeric-exact) < mp.mpf("1e-55"),
              mp.nstr(numeric, 28))
    check("warp_ratio_guard", (Q(101, 100)/Q(99, 100))**2 < Q(105, 100))
    check("seed_warp_value_guards",
          mp.cos(mp.mpf(".05"))*mp.exp(-mp.mpf(".00001")) > mp.mpf(".99") and
          mp.cosh(mp.mpf(".05"))*mp.exp(mp.mpf(".00001")) < mp.mpf("1.01"))
    check("seed_warp_derivative_guard", mp.tan(mp.mpf(".05")) < mp.mpf(".1"))
    check("finite_background_warp_derivative",
          2+24*y0+Q(1, 10**5)*y0**2*(1+loghi) < 4)
    check("ward_value_bound", Q(105, 100)*4*I1 < Q(71, 1000))
    check("ward_first_bound", Q(105, 100)*(40*I1+100*I0) < 100)
    check("ward_second_bound", Q(105, 100)*(400*I1+2000*I0) < 2000)
    check("matching_kernel_bounds",
          Q(105, 100)*10 < 11 and Q(105, 100)*100 < 110)
    K = 2000
    cmin = 6*K**3
    ustarlo = 1-Q(4, 10**26)
    lamlo = 1-Q(1, K**10)
    eta_x_upper = (6*Q(cmin, cmin-1))**2/(K**10*lamlo*ustarlo**2)
    check("retained_physical_mass_prefactor", eta_x_upper < Q(5, 10**32))
    K0, K1, K2 = map(Q, [4*10**25, 2*10**25, 10**30])
    check("retained_scalar_value_budget",
          2*10**20*400**2+10**25*Q(1, 10**7)+1 < K0)
    check("first_source_cauchy_budget",
          10**25*(Q(1, 10**7)+Q(1, 40000))*40000+Q(48, 100) < K1)
    p0 = Q(71, 1000)*K0
    p1 = Q(71, 1000)*K1+100*K0
    p2 = Q(71, 1000)*K2+200*K1+2000*K0
    check("composite_pressure_constants",
          p0 == Q(284, 100)*10**24 and p1 == Q(400142, 100000)*10**27 and
          p2 == Q(155, 100)*10**29)
    dmax = Q(21, 50)
    mass = 2*dmax*Q(5, 10**32)*(y0*(p2+2*p1)+I0*(K2+2*K1))
    lapse = 16*dmax**2*Q(85, 10**34)*(
        J3*(2*p2+16*p1+28*p0)+I3*(K2+8*K1+14*K0))
    check("mass_second_response", mass < Q(39, 1000), float(mass))
    check("lapse_second_response", lapse < Q(46, 10**8), float(lapse))
    check("component_taylor_bounds",
          Q(39, 2000)*Q(1, 10**14) == Q(195, 10**18) and
          Q(23, 10**8)*Q(1, 10**14) == Q(23, 10**22))
    memory_m = [2*dmax*y0*t for t in [Q(105, 100), Q(105, 100)+11, 132]]
    memory_a = [32*dmax**2*J3*t for t in [
        Q(105, 100), 4*Q(105, 100)+11, 14*Q(105, 100)+88+110]]
    check("mass_memory_coefficient_bounds",
          all(a < b for a, b in zip(memory_m, [Q(5,100), Q(51,100), Q(56,10)])))
    check("lapse_memory_coefficient_bounds",
          all(a < b for a, b in zip(memory_a, [Q(1,10**5), Q(14,10**5), Q(2,1000)])))
    NUMERICS.append({"mass_second_upper": str(mass), "lapse_second_upper": str(lapse),
                     "pressure_constants": list(map(str, [p0,p1,p2])),
                     "mass_memory": list(map(str,memory_m)),
                     "lapse_memory": list(map(str,memory_a))})


class Jet:
    """Value, ds, dt, dsdt. No factorial-scaled coefficients."""
    def __init__(self, value, ds=0, dt=0, dst=0):
        self.v = [mp.mpf(x) for x in [value, ds, dt, dst]]
    @staticmethod
    def cast(other):
        return other if isinstance(other, Jet) else Jet(other)
    def __add__(self, other):
        other = self.cast(other)
        return Jet(*(a+b for a,b in zip(self.v, other.v)))
    __radd__ = __add__
    def __neg__(self):
        return Jet(*(-x for x in self.v))
    def __sub__(self, other):
        return self + -self.cast(other)
    def __rsub__(self, other):
        return self.cast(other) + -self
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


def manufactured_checks():
    mp.mp.dps = 55
    y0, d = mp.mpf(".05"), mp.mpf(".4")
    targets = []
    frozen_targets = []
    for n in [32, 64]:
        nodes, weights = leggauss(n)
        nodes = [(mp.mpf(float(x))+1)/2 for x in nodes]
        weights = [mp.mpf(float(w))/2 for w in weights]
        def field_data(Y):
            z = Y/y0
            loga0 = mp.mpf(".03")*z*z-mp.mpf(".02")*z**3
            w0 = (mp.mpf(".06")*z-mp.mpf(".06")*z*z)/y0
            h = z**3
            k = -mp.mpf(".3")*z**3+mp.mpf(".2")*z**4
            hp = 3*z*z/y0
            kp = (-mp.mpf(".9")*z*z+mp.mpf(".8")*z**3)/y0
            u0 = d*Y*(1+mp.mpf(".01")*z*z-mp.mpf(".005")*z**3)
            u0p = d*(1+mp.mpf(".03")*z*z-mp.mpf(".02")*z**3)
            du = d*Y*(-z**3/4+mp.mpf(".05")*z**4)
            dup = d*(-z**3+mp.mpf(".25")*z**4)
            dv = d*Y*(mp.mpf(".075")*z**3-mp.mpf(".02")*z**4)
            dvp = d*(mp.mpf(".3")*z**3-mp.mpf(".1")*z**4)
            H = 1+mp.log(Y)**2
            T = [H, H*mp.mpf(".2")*(1+z), H*mp.mpf(".3")*(1-z),
                 H*mp.mpf(".17")*(1+z*z)]
            return loga0,w0,h,k,hp,kp,u0,u0p,du,dup,dv,dvp,T
        # Store a fixed quadrature data set, used by separately written
        # jet and scalar evaluators. Y=y0*z^6 resolves the log endpoint.
        outer = []
        for z, wt in zip(nodes, weights):
            Y = y0*z**6
            ow = wt*6*y0*z**5
            inner = [(wt2*(y0-Y), field_data(Y+(y0-Y)*z2))
                     for z2,wt2 in zip(nodes,weights)]
            outer.append((ow,field_data(Y),inner))
        boundary = field_data(y0)
        def jets(data, freeze=False):
            l,w,h,k,hp,kp,u,up,du,dup,dv,dvp,T = data
            if freeze:
                h=k=hp=kp=mp.mpf(0)
            return Jet(l,-h,-k),Jet(w,-hp,-kp),Jet(u,du,dv),Jet(up,dup,dvp),Jet(*T)
        def jet_maps(freeze=False, no_memory=False):
            mb,ab = Jet(0),Jet(0)
            lb = jets(boundary,freeze)[0]
            Pm = Jet(0) if no_memory else Jet(".7",".2","-.15",".11")
            for wt,data,inner in outer:
                l,w,u,up,T = jets(data,freeze)
                P = (2*(lb-l)).exp()*Pm
                for iw,idata in inner:
                    il,iwlog,_,_,iT = jets(idata,freeze)
                    P -= iw*(2*(il-l)).exp()*iwlog*iT
                mb += wt*up*(P-T)
                ab -= wt*(u*u*u/up)*(2*P-T)
            return mb,ab
        def scalar_maps(ss,tt,no_memory=False):
            def values(data):
                l,w,h,k,hp,kp,u,up,du,dup,dv,dvp,T = data
                return l-ss*h-tt*k,w-ss*hp-tt*kp,u+ss*du+tt*dv,(
                    up+ss*dup+tt*dvp),T[0]+ss*T[1]+tt*T[2]+ss*tt*T[3]
            lb = values(boundary)[0]
            Pm = mp.mpf(0) if no_memory else mp.mpf(".7")+mp.mpf(".2")*ss-mp.mpf(".15")*tt+mp.mpf(".11")*ss*tt
            mass,lapse = mp.mpf(0),mp.mpf(0)
            for wt,data,inner in outer:
                l,w,u,up,T = values(data)
                P = mp.exp(2*(lb-l))*Pm
                for iw,idata in inner:
                    il,iwlog,_,_,iT = values(idata)
                    P -= iw*mp.exp(2*(il-l))*iwlog*iT
                mass += wt*up*(P-T)
                lapse -= wt*u**3/up*(2*P-T)
            return mass,lapse
        jmass,jlapse = jet_maps()
        analytic = [jmass.v[3],jlapse.v[3]]
        targets.append(analytic)
        frozen_targets.append([v.v[3] for v in jet_maps(freeze=True)])
        errors = [[],[]]
        for hs in [".001",".0005"]:
            step = mp.mpf(hs)
            pp = scalar_maps(step,step)
            pm = scalar_maps(step,-step)
            mpv = scalar_maps(-step,step)
            mm = scalar_maps(-step,-step)
            for j in range(2):
                observed = (pp[j]-pm[j]-mpv[j]+mm[j])/(4*step**2)
                errors[j].append(abs(observed-analytic[j]))
        for j,label in enumerate(["mass","lapse"]):
            check(label+"_independent_mixed_difference_"+str(n),
                  max(errors[j]) < mp.mpf("2e-6")*(1+abs(analytic[j])),
                  {"target":mp.nstr(analytic[j],28),"errors":[mp.nstr(e,12) for e in errors[j]]})
            check(label+"_step_refinement_"+str(n),
                  errors[j][1] < mp.mpf(".35")*errors[j][0] or max(errors[j]) < mp.mpf("1e-35"))
        base = scalar_maps(mp.mpf(0),mp.mpf(0))
        nomem = scalar_maps(mp.mpf(0),mp.mpf(0),no_memory=True)
        check("nonzero_matching_pressure_detected_"+str(n),
              abs(base[0]-nomem[0]) > mp.mpf(".001") and abs(base[1]-nomem[1]) > mp.mpf("1e-8"))
        NUMERICS.append({"quadrature_nodes":n,
                         "mixed_mass":mp.nstr(analytic[0],35),"mixed_lapse":mp.nstr(analytic[1],35),
                         "mass_errors":[mp.nstr(e,18) for e in errors[0]],
                         "lapse_errors":[mp.nstr(e,18) for e in errors[1]],
                         "matching_mass_effect":mp.nstr(base[0]-nomem[0],25),
                         "matching_lapse_effect":mp.nstr(base[1]-nomem[1],25)})
    for j,label in enumerate(["mass","lapse"]):
        check(label+"_quadrature_refinement",
              abs(targets[1][j]-targets[0][j]) < mp.mpf("2e-6")*(1+abs(targets[1][j])),
              mp.nstr(abs(targets[1][j]-targets[0][j]),18))
    check("omitted_ward_kernel_variation_detected",
          abs(targets[-1][0]-frozen_targets[-1][0]) > mp.mpf("1e-4") and
          abs(targets[-1][1]-frozen_targets[-1][1]) > mp.mpf("1e-8"))


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
        "story":"S147","created_utc":datetime.now(timezone.utc).isoformat(),
        "status":"passed" if ok else "failed","python":sys.version,
        "numpy":np.__version__,"sympy":s.__version__,"mpmath":mp.__version__,
        "input_sha256_before":before,"input_sha256_after":after,"unchanged_inputs":before==after,
        "checks":ROWS,"numerics":NUMERICS,"exception":error,
        "scope":"outer regular canonical Ward/constraint transfer with explicit matching memory; manufactured sources are diagnostics",
        "unresolved":["interior/state matching pressure and bolt response","other local and radius source maps",
                      "actual global inverse composition","coupled exterior/state and physical gates"],
    }
    if error:
        record["failed_checker_source_utf8"] = Path(__file__).read_text(encoding="utf-8")
        record["failed_note_source_utf8"] = (ROOT/INPUTS[1]).read_text(encoding="utf-8")
    output = ROOT/"evidence"/f"{STEM}.json"
    history = json.loads(output.read_text(encoding="utf-8")) if output.exists() else {"runs":[]}
    history["runs"].append(record)
    output.write_text(json.dumps(history,indent=2)+"\n",encoding="utf-8",newline="\n")
    with (ROOT/"evidence"/f"{STEM}.log").open("a",encoding="utf-8",newline="\n") as log:
        log.write(json.dumps({"created_utc":record["created_utc"],"status":record["status"],"inputs":before})+"\n")
        for row in ROWS:
            log.write(json.dumps(row)+"\n")
        if error:
            log.write(error+"\n")
    print(json.dumps({"status":record["status"],"checks":len(ROWS),
                      "failed":sum(not row["passed"] for row in ROWS)}))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
