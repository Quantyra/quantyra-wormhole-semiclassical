"""S141: constraint checks and outward rational benchmark enclosures.

Run with Python -B. SymPy performs direct horizon checks; mpmath independently
solves the defining polynomials. Numerical agreement is supplementary to the
rational enclosures. No exterior equation or stability claim is tested here.
"""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json
import platform
import sys
import traceback

import mpmath as mp
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
STEM = "retained-benchmark-parameter-recovery"
INPUTS = [
    "checks/retained_benchmark_parameter_recovery.py",
    "notes/retained-benchmark-parameter-recovery.md",
    "notes/uniform-dirac-reference.md",
    "notes/stationary-horizon-product-matching.md",
    "notes/stationary-centered-covariant-seed.md",
    "notes/stationary-reference-cold-bolt.md",
    "notes/stationary-scalar-phase-moment.md",
    "notes/stationary-disk-and-normal-inverse.md",
    "notes/classical-radius-clock-width.md",
    "checks/classical_radius_clock_width.py",
]
GRID = 10**180
ROWS: list[dict] = []
VALUES: dict[str, dict] = {}


def floor_grid(x: F) -> F:
    return F(x.numerator * GRID // x.denominator, GRID)


def ceil_grid(x: F) -> F:
    return -floor_grid(-x)


@dataclass(frozen=True)
class I:
    lo: F
    hi: F

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError("reversed enclosure")

    @staticmethod
    def point(x):
        if isinstance(x, I):
            return x
        if isinstance(x, float):
            raise TypeError("binary floats are excluded from interval inputs")
        return I(F(x), F(x))

    @staticmethod
    def enclose(lo, hi):
        return I(floor_grid(F(lo)), ceil_grid(F(hi)))

    def __add__(self, other):
        b = I.point(other)
        return I.enclose(self.lo + b.lo, self.hi + b.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -I.point(other)

    def __rsub__(self, other):
        return I.point(other) - self

    def __mul__(self, other):
        b = I.point(other)
        products = [a * c for a in (self.lo, self.hi)
                    for c in (b.lo, b.hi)]
        return I.enclose(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        b = I.point(other)
        if b.lo <= 0 <= b.hi:
            raise ZeroDivisionError("interval denominator includes zero")
        return self * I.enclose(1 / b.hi, 1 / b.lo)

    def __rtruediv__(self, other):
        return I.point(other) / self

    def __pow__(self, n: int):
        if n < 0:
            return 1 / (self ** -n)
        result = I.point(1)
        for _ in range(n):
            result = result * self
        return result

    def sqrt(self):
        if self.lo < 0:
            raise ValueError("negative square root")
        low = isqrt(self.lo.numerator * GRID**2 // self.lo.denominator)
        high = isqrt(self.hi.numerator * GRID**2 // self.hi.denominator)
        # The extra endpoint is safe even when hi is a perfect square.
        return I(F(low, GRID), F(high + 1, GRID))

    def contains(self, x: F) -> bool:
        return self.lo <= x <= self.hi


def atan_bounds(x: F, terms: int) -> I:
    """Alternating Taylor theorem, 0<x<1."""
    if not 0 < x < 1:
        raise ValueError("arctangent range")
    s = sum(((-1)**j * x**(2*j + 1) / (2*j + 1)
             for j in range(terms)), F(0))
    next_term = (-1)**terms * x**(2*terms + 1) / (2*terms + 1)
    return I.enclose(min(s, s + next_term), max(s, s + next_term))


def log_unit(x: F, terms: int = 150) -> I:
    """log(x) for 1<=x<=2, with the positive atanh-series tail."""
    if not 1 <= x <= 2:
        raise ValueError("log reduction range")
    y = I.point((x - 1) / (x + 1))
    y2 = y * y
    power, total = y, I.point(0)
    for j in range(terms):
        total += power / (2*j + 1)
        power *= y2
    tail = 2 * power / ((2*terms + 1) * (1 - y2))
    return I.enclose(2 * total.lo, 2 * total.hi + tail.hi)


def log_point(x: F) -> I:
    if x <= 0:
        raise ValueError("log positive argument required")
    if x < 1:
        return -log_point(1 / x)
    exponent = 0
    while x >= 2:
        x /= 2
        exponent += 1
    return log_unit(x) + exponent * log_unit(F(2))


def log_interval(x: I) -> I:
    return I(log_point(x.lo).lo, log_point(x.hi).hi)


def decimal_bound(x: F, rounding: str) -> str:
    with localcontext() as ctx:
        ctx.prec = 34
        ctx.rounding = rounding
        return str(Decimal(x.numerator) / Decimal(x.denominator))


def keep(name: str, x: I):
    VALUES[name] = {
        "lower": decimal_bound(x.lo, ROUND_FLOOR),
        "upper": decimal_bound(x.hi, ROUND_CEILING),
        "exact_rational_lower": str(x.lo),
        "exact_rational_upper": str(x.hi),
    }


def check(name: str, condition, detail=""):
    row = {"name": name, "passed": bool(condition), "detail": str(detail)}
    if any(r["name"] == name for r in ROWS):
        raise ValueError("duplicate check")
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not row["passed"]:
        raise AssertionError(name)


def mp_fraction(x) -> F:
    """Exact rational value of an mpmath binary mantissa/exponent."""
    sign, mantissa, exponent, _ = x._mpf_
    value = F((-1)**sign * mantissa)
    return value * 2**exponent if exponent >= 0 else value / 2**(-exponent)


def run():
    # Direct field-equation substitution precedes the parameter evaluation.
    r, b, lam, Q, m = sp.symbols("r b lambda0 Q m", positive=True)
    f = 1 + lam*r**2 - 2*m/r + Q/r**2
    root_data = {Q: b**2 + 3*lam*b**4, m: b*(1 + 2*lam*b**2)}
    fext = f.subs(root_data)
    check("classical_double_horizon",
          sp.simplify(fext.subs(r, b)) == 0
          and sp.simplify(sp.diff(fext, r).subs(r, b)) == 0)
    factor = (r-b)**2/r**2 * (1 + lam*(r*r + 2*b*r + 3*b*b))
    check("independent_exterior_factorization",
          sp.factor(fext - factor) == 0)
    check("orbit_radius_from_second_metric_derivative",
          sp.simplify(sp.diff(fext, r, 2).subs(r, b)/2
                      - (1+6*lam*b*b)/(b*b)) == 0)
    y = sp.symbols("y", positive=True)
    check("classical_radius_unique_positive_root",
          sp.diff(y + 3*lam*y*y, y).is_positive)
    A = sp.symbols("A", positive=True)
    vacuum_root = 2*lam/(1 + sp.sqrt(1+4*A*lam))
    check("vacuum_root_solves_defining_equation",
          sp.simplify(vacuum_root + A*vacuum_root**2-lam) == 0)
    check("vacuum_positive_root_unique",
          sp.diff(y+A*y*y, y).is_positive)

    # pi is enclosed by Machin's identity and alternating-series remainders.
    pi = 16*atan_bounds(F(1, 5), 160)-4*atan_bounds(F(1, 239), 50)
    check("pi_rational_enclosure",
          F(333, 106) < pi.lo < pi.hi < F(22, 7)
          and pi.hi-pi.lo < F(1, 10**170))
    K, eta = 2000, F(1, 2000)
    q = K**4
    lambda0 = I.point(F(1, K**10))
    C0 = (1 + 12*pi*K**6).sqrt()
    xi0 = (C0-1)/6
    b02 = K**10*xi0
    b0 = b02.sqrt()
    a0 = (b02/C0).sqrt()
    mu0 = 1/b0
    charge = pi*K**16
    check("derived_radius_satisfies_charge_constraint",
          (b02+3*lambda0*b02*b02-charge).contains(F(0)))
    check("retained_classical_charge_and_orbit_bounds",
          K**3 < xi0.lo < xi0.hi < 2*K**3
          and C0.hi < F(123, 20)*K**3)

    vacuum_A = 11/(360*pi)
    lambda_e = 2*lambda0/(1+(1+4*vacuum_A*lambda0).sqrt())
    delta_e = vacuum_A*lambda_e*lambda_e/lambda0
    check("enclosed_vacuum_residual_contains_zero",
          (lambda_e+vacuum_A*lambda_e*lambda_e-lambda0).contains(F(0)))
    check("physical_curvature_is_in_original_guard",
          lambda_e.lo > F(999, 1000)*lambda0.hi
          and lambda_e.hi < lambda0.lo)
    check("fractional_vacuum_shift_below_1e_minus_34",
          0 < delta_e.lo < delta_e.hi < F(1, 10**34))
    lambda_test = F(1999, 2000)*lambda0
    test_residual = (lambda_test+vacuum_A*lambda_test**2-lambda0)/lambda0
    check("old_spectral_point_rejected_by_vacuum_equation",
          F(-1, 2000) < test_residual.lo <= test_residual.hi
          < F(-1, 2000) + F(1, 10**34))

    x_classical = lambda_e*b02
    kappa_M = F(-6, 5) + log_interval(x_classical)/2
    check("normal_symbol_constant_determined_and_narrow",
          kappa_M.hi-kappa_M.lo < F(1, 10**20)
          and F(39, 5) < kappa_M.lo < kappa_M.hi < F(54, 5))
    nu = 2*q*eta*C0/(b0**3)
    kappa_low, kappa_high = nu/4, 7*nu/25
    beta_low, beta_high = 50*pi/(7*nu), 8*pi/nu
    check("full_prescribed_temperature_family_positive",
          0 < kappa_low.lo < kappa_low.hi < kappa_high.lo)
    check("full_prescribed_period_family_positive",
          0 < beta_low.lo < beta_low.hi < beta_high.lo)
    check("period_temperature_endpoint_pairing",
          (beta_low*kappa_high-2*pi).contains(F(0))
          and (beta_high*kappa_low-2*pi).contains(F(0)))
    E = 60*pi+F(11, 3)*lambda_e
    k_rad = (E-F(9, 4)*lambda_e).sqrt()
    c0 = k_rad*(E+4*lambda_e)
    check("scalar_phase_normalization_fixed",
          F(258792, 100) < c0.lo < c0.hi < F(258793, 100))

    dz, du = F(14, 10**23), F(4, 10**26)
    Cstar = C0 * I(1-dz, 1+dz)/I(1-du, 1+du)
    M1 = ((1+q)/Cstar).sqrt()
    check("proper_product_mass_gap",
          18 < M1.lo < M1.hi < 19)
    n = sp.symbols("n", integer=True, positive=True)
    check("all_level_mass_ratio_identity",
          sp.simplify(n*(n+q)/(1+q)-n
                      - n*(n-1)/(1+q)) == 0,
          "n(n-1)/(1+q)>=0 for every integer n>=1")
    xstar = x_classical*I(1-dz, 1+dz)
    check("loose_guard_witness_excluded_from_product_box",
          xstar.hi < 82*10**8 < 10**10)
    heat = 130000000/(60*pi*M1**3)
    check("scalar_heat_product_stays_large",
          100 < heat.lo < heat.hi < 124)
    for name, value in {
        "C0": C0, "xi0": xi0, "b0_planck": b0, "a0_planck": a0,
        "mu0_inverse_planck": mu0, "lambda0": lambda0,
        "lambda_e": lambda_e, "fractional_vacuum_shift": delta_e,
        "old_point_normalized_residual": test_residual,
        "lambda_e_b0_squared": x_classical, "kappa_M": kappa_M,
        "nu": nu, "kappa_lower_endpoint": kappa_low,
        "kappa_upper_endpoint": kappa_high,
        "beta_lower_endpoint": beta_low, "beta_upper_endpoint": beta_high,
        "scalar_c0": c0, "Cstar_product_box": Cstar,
        "M1_product_box": M1, "xstar_product_box": xstar,
        "formal_heat_product": heat,
    }.items():
        keep(name, value)

    # Independent method: solve polynomials, rather than evaluating their
    # closed radical expressions. Scaling resolves the 10^-35 vacuum shift.
    for precision in (200, 260):
        with mp.workdps(precision):
            kk = mp.mpf(K)
            xi_num = mp.findroot(lambda x: x+3*x*x-mp.pi*kk**6,
                                 (kk**3, 2*kk**3))
            lam0_num = kk**-10
            a_num = 11/(360*mp.pi)
            # lambda_e=lambda0*(1-A*lambda0*t);
            # the defining equation is t=(1-A*lambda0*t)^2.
            t = mp.findroot(lambda v: v-(1-a_num*lam0_num*v)**2,
                            (mp.mpf(".9"), mp.mpf("1.1")))
            lam_num = lam0_num*(1-a_num*lam0_num*t)
            b_num = mp.sqrt(kk**10*xi_num)
            nu_num = 2*kk**4*mp.mpf(eta.numerator)/eta.denominator * (
                1+6*xi_num)/b_num**3
            km_num = -mp.mpf(6)/5+mp.log(lam_num*b_num*b_num)/2
            comparisons = [
                pi.contains(mp_fraction(mp.pi)),
                xi0.contains(mp_fraction(xi_num)),
                b0.contains(mp_fraction(b_num)),
                lambda_e.contains(mp_fraction(lam_num)),
                delta_e.contains(mp_fraction(a_num*lam0_num*t)),
                nu.contains(mp_fraction(nu_num)),
                kappa_M.contains(mp_fraction(km_num)),
            ]
            check(f"independent_polynomial_solve_{precision}_digits",
                  all(comparisons), comparisons)

    # Audit input provenance separately; these checks are not mathematical
    # substitutes for the derivation or root computations above.
    source = (ROOT/"notes/uniform-dirac-reference.md").read_text("utf-8")
    check("preexisting_explicit_radius_definition",
          "b0^2=L4^2(C0-1)/6" in source
          and "L4/l_p=k^5" in source)


def hashes():
    return {path: sha256((ROOT/path).read_bytes()).hexdigest()
            for path in INPUTS}


def main():
    before = hashes()
    error = None
    try:
        run()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr)
    after = hashes()
    ok = error is None and all(r["passed"] for r in ROWS) and before == after
    receipt = {
        "story": "S141", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if ok else "failed",
        "python": sys.version, "platform": platform.platform(),
        "sympy": sp.__version__, "mpmath": mp.__version__,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after,
        "interval_method": "outward exact rational arithmetic at 10^-180; "
                           "Machin alternating bounds; integer square roots; "
                           "range-reduced atanh log with geometric tail",
        "checks": ROWS, "enclosures": VALUES, "exception": error,
        "scope": "retained parameter recovery; inherited product bounds "
                 "are inputs, not re-proved determinant estimates",
        "unresolved": ["actual quantum exterior and boundary completion",
                       "signed source and seam trace", "variable vertices",
                       "nonlinear remainder", "opening and area",
                       "stability, transport, preparation and creation"],
    }
    # Append runs; retain any failed execution rather than overwrite it.
    target = ROOT/"evidence"/f"{STEM}.json"
    history = json.loads(target.read_text("utf-8")) if target.exists() else {"runs": []}
    history["runs"].append(receipt)
    target.write_text(json.dumps(history, indent=2)+"\n",
                      encoding="utf-8", newline="\n")
    with (ROOT/"evidence"/f"{STEM}.log").open(
            "a", encoding="utf-8", newline="\n") as log:
        log.write(json.dumps({"created_utc": receipt["created_utc"],
                              "status": receipt["status"],
                              "input_sha256": before})+"\n")
        for row in ROWS:
            log.write(json.dumps(row)+"\n")
        if error:
            log.write(error+"\n")
    print(json.dumps({"status": receipt["status"], "checks": len(ROWS),
                      "failed": sum(not r["passed"] for r in ROWS)}))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
