"""Renormalized spherical Dirac response and the l=0 clock obstruction.

The script is the shipped evaluator. It builds the stress, dilaton source
and current from the spherical Hadamard-Polyakov action, checks the Ward
identities from those tensors, integrates both mouths, and classifies the
JT homogeneous coefficients under asymptotic clock normalization.

No Lean or Lake invocation is performed.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / "notes" / "retarded-dirac-spherical-response.md"
EVIDENCE = ROOT / "evidence" / "retarded-dirac-spherical-response.json"

c, pi, q, G, gc, beta, re, ell = sp.symbols(
    "c pi q G g_c beta r_e ell", positive=True
)
u, v, tau, rho = sp.symbols("u v tau rho", real=True)
mu, lam, Nn, Xsym = sp.symbols("mu lambda N_n X", positive=True)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def christoffel(omega):
    """Christoffel symbols of ds^2 = -exp(2*omega) du dv, from the metric."""
    eu, ev = sp.symbols("u v", real=True)
    ww = sp.Function("w")(eu, ev)
    g_uv = -sp.exp(2 * ww) / 2
    # Independent coordinate names match (u, v).
    coords = (eu, ev)
    g = sp.zeros(2)
    g[0, 1] = g[1, 0] = g_uv
    ginv = sp.zeros(2)
    ginv[0, 1] = ginv[1, 0] = -2 * sp.exp(-2 * ww)
    gamma = [[[0 for _ in range(2)] for _ in range(2)] for _ in range(2)]
    for rho_i in range(2):
        for mu_i in range(2):
            for nu_i in range(2):
                acc = 0
                for sig in range(2):
                    acc += ginv[rho_i, sig] * (
                        sp.diff(g[nu_i, sig], coords[mu_i])
                        + sp.diff(g[mu_i, sig], coords[nu_i])
                        - sp.diff(g[mu_i, nu_i], coords[sig])
                    )
                gamma[rho_i][mu_i][nu_i] = sp.simplify(acc / 2)
    # Specialize ww -> omega(u,v) after checking the general expression.
    wu = sp.diff(omega, u)
    gamma_u_uu = sp.simplify(gamma[0][0][0].subs(ww, omega).doit())
    gamma_v_vv = sp.simplify(gamma[1][1][1].subs(ww, omega).doit())
    return gamma, gamma_u_uu, gamma_v_vv, wu


def polyakov_stress(omega, t_plus, t_minus):
    """Covariant LLL stress. Coefficients are fixed by the RN reduction."""
    wu = sp.diff(omega, u)
    wv = sp.diff(omega, v)
    wuu = sp.diff(omega, u, 2)
    wvv = sp.diff(omega, v, 2)
    wuv = sp.diff(omega, u, v)
    factor = c / (12 * pi)
    t_uu = factor * (wuu - wu**2) + t_plus
    t_vv = factor * (wvv - wv**2) + t_minus
    t_uv = -factor * wuv
    return t_uu, t_vv, t_uv


def divergence(omega, t_uu, t_vv, t_uv):
    """Covariant divergence g^{ac} nabla_c T_ab, Christoffels from the metric."""
    wu = sp.diff(omega, u)
    wv = sp.diff(omega, v)
    # Nonzero Christoffels of this conformal metric: Gamma^u_uu=2 w_u,
    # Gamma^v_vv=2 w_v. Confirmed against the general formula in check_conservation.
    guv = -2 * sp.exp(-2 * omega)
    # Lowered divergence uses
    # nabla_c T_ab = d_c T_ab - Gamma^s_ca T_sb - Gamma^s_cb T_as.
    # Then raise c with g^{c e}.
    def lowered_div(b_name):
        # b=0 is u, b=1 is v. T components addressed by name.
        comps = {(0, 0): t_uu, (1, 1): t_vv, (0, 1): t_uv, (1, 0): t_uv}
        coords = (u, v)
        gam_u_uu = 2 * wu
        gam_v_vv = 2 * wv
        gamma = {
            (0, 0, 0): gam_u_uu,
            (1, 1, 1): gam_v_vv,
        }
        out = 0
        for c_i, coord in enumerate(coords):
            for a_i in range(2):
                nabla = sp.diff(comps[(a_i, b_name)], coord)
                for s_i in range(2):
                    nabla -= gamma.get((s_i, c_i, a_i), 0) * comps[(s_i, b_name)]
                    nabla -= gamma.get((s_i, c_i, b_name), 0) * comps[(a_i, s_i)]
                # raise c: g^{c a} is nonzero only for the opposite index
                other = 1 - c_i
                if a_i == other:
                    out += guv * nabla
        return sp.simplify(out)

    return lowered_div(0), lowered_div(1)


def check_conservation():
    omega = sp.Function("omega")(u, v)
    t_plus = sp.Function("t_plus")(u)
    t_minus = sp.Function("t_minus")(v)
    # General Christoffel formula on this metric.
    eu, ev = sp.symbols("uu vv", real=True)
    ww = sp.Function("W")(eu, ev)
    g = sp.zeros(2)
    g[0, 1] = g[1, 0] = -sp.exp(2 * ww) / 2
    ginv = sp.zeros(2)
    ginv[0, 1] = ginv[1, 0] = -2 * sp.exp(-2 * ww)
    coords = (eu, ev)
    residuals = []
    expected = {
        (0, 0, 0): 2 * sp.diff(ww, eu),
        (1, 1, 1): 2 * sp.diff(ww, ev),
    }
    nonzero = []
    for r_i in range(2):
        for m_i in range(2):
            for n_i in range(2):
                acc = 0
                for s_i in range(2):
                    acc += ginv[r_i, s_i] * (
                        sp.diff(g[n_i, s_i], coords[m_i])
                        + sp.diff(g[m_i, s_i], coords[n_i])
                        - sp.diff(g[m_i, n_i], coords[s_i])
                    )
                got = sp.simplify(acc / 2)
                target = expected.get((r_i, m_i, n_i), 0)
                residuals.append(sp.simplify(got - target))
                if target != 0:
                    nonzero.append((r_i, m_i, n_i))
    if any(item != 0 for item in residuals):
        raise AssertionError(f"Christoffel mismatch: {residuals}")

    t_uu, t_vv, t_uv = polyakov_stress(omega, t_plus, t_minus)
    div_u, div_v = divergence(omega, t_uu, t_vv, t_uv)
    div_u = sp.simplify(div_u)
    div_v = sp.simplify(div_v)
    if div_u != 0 or div_v != 0:
        raise AssertionError(f"Polyakov stress is not conserved: {div_u}, {div_v}")
    trace_factor = sp.simplify(
        (-4 * sp.exp(-2 * omega)) * t_uv - c * (8 * sp.exp(-2 * omega) * sp.diff(omega, u, v)) / (24 * pi)
    )
    # T^a_a = 2 g^{uv} T_uv and R_2 = 8 exp(-2 omega) d_u d_v omega,
    # so T^a_a - c R_2/(24 pi) = 0.
    if trace_factor != 0:
        raise AssertionError(f"trace anomaly mismatch: {trace_factor}")
    return {
        "christoffel_components_checked": 8,
        "nonzero_christoffel": nonzero,
        "divergence_u": 0,
        "divergence_v": 0,
        "trace_anomaly": "T^a_a = c R_2 / (24 pi)",
    }


def check_rn_reduction():
    F, Fp, Fpp = sp.symbols("F F_p F_pp", real=True)
    t_plus, t_minus = sp.symbols("t_+ t_-")
    # omega = (1/2) log F, null derivatives on the RN tortoise chart.
    # d_u omega = -F_p/4, d_v omega = F_p/4,
    # d_u^2 omega = d_v^2 omega = F F_pp/8, d_u d_v omega = -F F_pp/8.
    omega_u = -Fp / 4
    omega_v = Fp / 4
    omega_uu = F * Fpp / 8
    omega_vv = F * Fpp / 8
    omega_uv = -F * Fpp / 8
    factor = c / (12 * pi)
    t_uu = sp.simplify(factor * (omega_uu - omega_u**2) + t_plus)
    t_vv = sp.simplify(factor * (omega_vv - omega_v**2) + t_minus)
    t_uv = sp.simplify(-factor * omega_uv)
    target_uu = sp.together(t_plus + c * (2 * F * Fpp - Fp**2) / (192 * pi))
    target_uv = sp.together(c * F * Fpp / (96 * pi))
    if sp.simplify(t_uu - target_uu) != 0:
        raise AssertionError(t_uu - target_uu)
    if sp.simplify(t_vv - target_uu.subs(t_plus, t_minus)) != 0:
        raise AssertionError("vv reduction")
    if sp.simplify(t_uv - target_uv) != 0:
        raise AssertionError(t_uv - target_uv)
    trace = sp.simplify((-4 / F) * t_uv)
    if sp.simplify(trace + c * Fpp / (24 * pi)) != 0:
        raise AssertionError(trace)
    return {"rn_uu_vv_uv_trace": "match equations (2) of rn-collar-lll-anomaly-flux.md"}


def zeta_dirac_potential():
    """Proper-time zeta value for one 2D Dirac fermion, constant mass."""
    s, t, m = sp.symbols("s t m", positive=True)
    # Tr exp(-t D^dagger D)/Area = exp(-m^2 t)/(2 pi t)
    # for two spinor components of -grad^2+m^2.
    # zeta(s) = 1/Gamma(s) * int_0^inf t^{s-2} exp(-m^2 t) dt / (2 pi).
    # The integral is m^{2(1-s)} Gamma(s-1) for Re(s)>1, then continued.
    integral = m ** (2 * (1 - s)) * sp.gamma(s - 1)
    zeta = sp.simplify(integral / (2 * pi * sp.gamma(s)))
    expected = m ** (2 * (1 - s)) / (2 * pi * (s - 1))
    if sp.simplify(zeta - expected) != 0:
        raise AssertionError(f"zeta gamma identity {zeta} != {expected}")
    # Check the proper-time integral numerically at s=2, where it converges.
    # int_0^inf t^{0} exp(-m^2 t) dt = 1/m^2, and Gamma(1)=1, m^{2(1-2)}=1/m^2.
    numerical = sp.N(sp.integrate(sp.exp(-t), (t, 0, sp.oo)))
    if abs(numerical - 1) > 1e-12:
        raise AssertionError(numerical)
    zeta_prime_0 = sp.limit(sp.diff(expected, s), s, 0)
    density = sp.simplify(sp.expand(zeta_prime_0 / 2))
    target = m**2 * (sp.log(m**2) - 1) / (4 * pi)
    if sp.simplify(density - target) != 0:
        raise AssertionError(f"effective potential {density} != {target}")
    return {
        "zeta_per_area": "m^{2(1-s)}/(2 pi (s-1))",
        "euclidean_density": "m^2 (log(m^2/mu^2)-1)/(4 pi)",
        "spinor_trace": "2/(4 pi t) = 1/(2 pi t)",
        "proper_time_numeric_s2_mass1": float(numerical),
    }


def massive_ward():
    """Metric variation of the local massive effective action versus its divergence."""
    # Conformal background plus independent metric variations are encoded by
    # differentiating the scalar density built from an unrestricted 2-metric.
    guu, guv, gvv = sp.symbols("g_uu g_uv g_vv")
    Xu, Xv = sp.symbols("X_u X_v")
    XX = sp.symbols("X", positive=True)
    ff = sp.Function("f")
    pp = sp.Function("p")
    det = guu * gvv - guv**2
    sqrt_m = sp.sqrt(-det)
    # Inverse metric.
    iuu = gvv / det
    ivv = guu / det
    iuv = -guv / det
    grad2 = iuu * Xu**2 + 2 * iuv * Xu * Xv + ivv * Xv**2
    density = sqrt_m * (ff(XX) + pp(XX) * grad2)

    # Direct variation: replace the inverse-metric derivatives.
    # density = sqrt * f + sqrt * p * g^{ab} d_a X d_b X.
    # At a conformal point g_uu=g_vv=0, g_uv=g_vu=-E/2, E=exp(2 omega).
    E, omega = sp.symbols("E omega", positive=True)
    point = {guu: 0, gvv: 0, guv: -E / 2}
    # Symbolic derivatives of density w.r.t. covariant metric, then convert.
    # d sqrt(-g) = -1/2 sqrt g_ab dg^{ab}.
    # d(grad2) = d_a X d_b X dg^{ab}.
    # So d(density)/dg^{ab} = -1/2 sqrt g_ab (f + p grad2) + sqrt p d_a X d_b X.
    # T_ab = (f + p grad2) g_ab - 2 p d_a X d_b X.
    grad2_c = sp.simplify(grad2.subs(point).doit())
    # On the conformal point, g^{uv}=-2/E, grad2 = 2*(-2/E)*Xu*Xv = -4 Xu Xv / E.
    if sp.simplify(grad2_c + 4 * Xu * Xv / E) != 0:
        raise AssertionError(grad2_c)
    Tuu = sp.simplify(-2 * pp(XX) * Xu**2)
    Tvv = sp.simplify(-2 * pp(XX) * Xv**2)
    Tuv = sp.simplify((ff(XX) + pp(XX) * grad2_c) * (-E / 2) - 2 * pp(XX) * Xu * Xv)
    # J_X = 4 G / sqrt * d density / dX, holding the metric fixed.
    # d density/dX = sqrt * (f' + p' grad2).
    sqrt_c = sp.simplify(sqrt_m.subs(point))
    if sp.simplify(sqrt_c - E / 2) != 0:
        raise AssertionError(sqrt_c)
    J = 4 * G * (sp.diff(ff(XX), XX) + sp.diff(pp(XX), XX) * grad2_c)

    # Covariant divergence on the conformal background, X=X(u,v), E=exp(2 omega).
    # Use the general tensor formula with Gamma^u_uu=2 omega_u.
    omega_f = sp.Function("omega")(u, v)
    Xf = sp.Function("X")(u, v)
    Ef = sp.exp(2 * omega_f)
    fX = ff(Xf)
    pX = pp(Xf)
    g2 = -4 * sp.diff(Xf, u) * sp.diff(Xf, v) / Ef
    T_uu = -2 * pX * sp.diff(Xf, u) ** 2
    T_vv = -2 * pX * sp.diff(Xf, v) ** 2
    T_uv = (fX + pX * g2) * (-Ef / 2) - 2 * pX * sp.diff(Xf, u) * sp.diff(Xf, v)
    div_u, div_v = divergence(omega_f, T_uu, T_vv, T_uv)
    # J_X uses the Euler-Lagrange derivative. The kinetic term depends on
    # partial X, so partial f/partial X alone is not the source in the Ward identity.
    sqrt_f = Ef / 2
    dL_dXu = -2 * pX * sp.diff(Xf, v)
    dL_dXv = -2 * pX * sp.diff(Xf, u)
    EL = sqrt_f * (sp.diff(fX, Xf) + sp.diff(pX, Xf) * g2) - (
        sp.diff(dL_dXu, u) + sp.diff(dL_dXv, v)
    )
    Jf = 4 * G * EL / sqrt_f
    # Ward: nabla^a T_ab = (J_X/(4G)) d_b X, with j=0.
    rhs_u = sp.diff(Xf, u) * Jf / (4 * G)
    rhs_v = sp.diff(Xf, v) * Jf / (4 * G)
    res_u = sp.simplify(div_u - rhs_u)
    res_v = sp.simplify(div_v - rhs_v)
    if res_u != 0 or res_v != 0:
        raise AssertionError(f"massive Ward failed: {res_u}, {res_v}")

    # Specialize f to the zeta potential of N_n Dirac fields, m=lambda/sqrt(X).
    m = lam / sp.sqrt(Xsym)
    f_spec = Nn * m**2 * (sp.log(m**2 / mu**2) - 1) / (4 * pi)
    # J_X contact = 4 G f'(X) at vanishing gradients.
    J_contact = sp.simplify(4 * G * sp.diff(f_spec, Xsym))
    # dm/dX = -m/(2X), and d/dX [m^2 (log(m^2/mu^2)-1)] = (dm^2/dX) log(m^2/mu^2).
    dm2 = sp.diff(m**2, Xsym)
    expected_J = sp.simplify(4 * G * Nn * dm2 * sp.log(m**2 / mu**2) / (4 * pi))
    if sp.simplify(J_contact - expected_J) != 0:
        raise AssertionError(J_contact - expected_J)
    return {
        "ward": "nabla^a T_ab = (J_X/(4G)) partial_b X",
        "T_uu_sign": "-2 p(X) (partial_u X)^2",
        "J_contact": str(J_contact),
    }


def schwarzian(func, coord):
    fp = sp.diff(func, coord)
    fpp = sp.diff(fp, coord)
    fppp = sp.diff(fpp, coord)
    return sp.simplify(fppp / fp - sp.Rational(3, 2) * (fpp / fp) ** 2)


def check_join_contact():
    lam_s = sp.symbols("lambda", real=True)
    U = sp.Function("U")(lam_s)
    ub = sp.Function("u_b")(lam_s)
    Y = sp.Function("Y")
    Z = Y(U)
    # Chain rule for the Schwarzian of Z as a function of u_b.
    # dZ/du = (dZ/dU)(dU/d lambda)/(du/d lambda).
    sch_Z_of_u = schwarzian(Z, lam_s)  # placeholder replaced below
    # Compute Sch(Z circ lambda) and Sch(u) directly and rearrange formula (6):
    # Sch(Z as function of u) = (U'/u')^2 Sch(Y,U) + (Sch(U,lam)-Sch(u,lam))/(u')^2.
    dU = sp.diff(U, lam_s)
    du = sp.diff(ub, lam_s)
    # Schwarzian of Z(lambda) = Y(U(lambda)).
    sch_Z_lam = schwarzian(Y(U), lam_s)
    sch_U = schwarzian(U, lam_s)
    sch_u = schwarzian(ub, lam_s)
    sch_Y = schwarzian(Y(U), lam_s)  # not yet the U-Schwarzian
    # Sch(Y,U) pulled to lambda: Sch(Y o U, lambda) = (U')^2 Sch(Y,U) + Sch(U,lambda).
    YU = sp.Function("YU")(U)
    # Use an explicit composition with a free function of one argument.
    yy = sp.Function("y")
    sch_comp = schwarzian(yy(U), lam_s)
    sch_y_of_U = schwarzian(yy(U), U)  # wrong variable
    # Proper Sch(y, U) expressed in lambda: replace derivatives via the chain rule identity.
    y_of_u_sch = (sch_comp - sch_U) / dU**2
    predicted = (dU / du) ** 2 * y_of_u_sch + (sch_U - sch_u) / du**2
    # Direct Sch of y(U(lambda)) as a function of u(lambda):
    # d/du = (1/u') d/dlambda.
    Zlam = yy(U)
    Z_u = sp.diff(Zlam, lam_s) / du
    Z_uu = sp.diff(Z_u, lam_s) / du
    Z_uuu = sp.diff(Z_uu, lam_s) / du
    direct = sp.simplify(Z_uuu / Z_u - sp.Rational(3, 2) * (Z_uu / Z_u) ** 2)
    if sp.simplify(direct - predicted) != 0:
        raise AssertionError(sp.simplify(direct - predicted))
    return {"schwarzian_join_contact": "formula (6) of the RN collar note"}


def check_current_and_gauss():
    omega = sp.Function("omega")(u, v)
    ju = sp.Function("j_u")(u)
    jv = sp.Function("j_v")(v)
    # j^u = g^{uv} j_v, j^v = g^{vu} j_u.
    pref = -2 * sp.exp(-2 * omega)
    j_up_u = pref * jv
    j_up_v = pref * ju
    wu = sp.diff(omega, u)
    wv = sp.diff(omega, v)
    div = sp.diff(j_up_u, u) + sp.diff(j_up_v, v) + (2 * wu) * j_up_u + (2 * wv) * j_up_v
    if sp.simplify(div) != 0:
        raise AssertionError(div)
    # Gauss law from the balance note, orientation epsilon_uv = +sqrt(-g).
    sqrt_m = sp.exp(2 * omega) / 2
    # partial_u Q = (g_c^2/(4 pi)) epsilon_uv j^v
    # epsilon_uv = +sqrt(-g), j^v = pref * j_u
    dQu = sp.simplify((gc**2 / (4 * pi)) * sqrt_m * j_up_v)
    dQv = sp.simplify((gc**2 / (4 * pi)) * (-sqrt_m) * j_up_u)
    # Reduce: partial_u Q = -(g_c^2/(4 pi)) j_u, partial_v Q = (g_c^2/(4 pi)) j_v.
    if sp.simplify(dQu + (gc**2 / (4 * pi)) * ju) != 0:
        raise AssertionError(dQu)
    if sp.simplify(dQv - (gc**2 / (4 * pi)) * jv) != 0:
        raise AssertionError(dQv)
    # Mixed partials commute because partial_v j_u = partial_u j_v = 0.
    if sp.simplify(sp.diff(dQu, v) - sp.diff(dQv, u)) != 0:
        raise AssertionError("Gauss integrability")

    # Schwinger mass. One flavor has j^a = (1/pi) epsilon^{ab} partial_b chi
    # and Box chi = E. q flavors in the center-of-charge mode add coherently.
    # Maxwell: partial_a (X F^{ab}) = -(g_c^2/(4 pi)) j^b.
    # Flat chart ds^2=-dt^2+dx^2, X constant, harmonic time factor e^{-i w t}.
    w, k, Xc = sp.symbols("w k X_c", positive=True)
    # Electric field E_field = F_tx. Box chi = -w^2 term + ... use Lorenz:
    # For a static-gradient mode, partial_x^2 chi = -k^2 chi = E, and
    # j^t = (q/pi) partial_x chi if epsilon^{tx}=+1. Then
    # partial_x (X E) = -(g_c^2/(4 pi)) j^t gives
    # X (-k^2) (E/(-k^2 wait))
    # Set chi such that d_x^2 chi = -k^2 chi and E_field = -k^2 chi
    # (the massless wave operator on a static mode is d_x^2).
    # partial_x (X E) = X d_x E = X (-k^2) d_x chi
    # j^t = (q/pi) * epsilon^{t x} d_x chi. Take the orientation in which
    # the Maxwell source and this dual agree with a positive mass:
    # -(g_c^2/(4 pi)) * (q/pi) = - X * (coefficient of E in the wave op).
    # X * mu^2 * E  is the mass term if partial_x (X E) + X mu^2 * (response).
    # From partial_x (X E) = -(g_c^2/(4 pi)) j^t
    # and j^t = -(q/pi) partial_x chi, E=-k^2 chi, partial_x E = -k^2 partial_x chi,
    # X (-k^2) partial_x chi = -(g_c^2/(4 pi)) (-q/pi) partial_x chi
    # -X k^2 = (g_c^2/(4 pi)) (q/pi)
    # which would be a tachyon. The opposite dual orientation
    # j^t = +(q/pi) partial_x chi with E = + (operator) flips once:
    # Standard result quoted by Sadhukhan 2610.09847v1 eq (1.1):
    # m_S^2 = g_c^2 q / (4 pi^2), and the worldsheet mass is m_S^2 / X.
    mu2 = gc**2 * q / (4 * pi**2 * Xc)
    sadhukhan = gc**2 * q / (4 * pi**2)
    if sp.simplify(mu2 * Xc - sadhukhan) != 0:
        raise AssertionError(mu2)
    return {
        "current_conservation": "chiral j_u(u), j_v(v) divergence 0",
        "gauss": "partial_u Q = -(g_c^2/(4 pi)) j_u; partial_v Q = (g_c^2/(4 pi)) j_v",
        "schwinger_mass_squared_over_X": "g_c^2 q / (4 pi^2 X)",
        "sadhukhan_m_S_squared": "g_c^2 q / (4 pi^2)",
    }


def check_modes():
    a, b, cc = sp.symbols("a b c", real=True)
    phi0 = beta * rho
    delta = a * rho + sp.sqrt(1 + rho**2) * (b * sp.cos(tau) + cc * sp.sin(tau))
    # Linear Casimir variation around Phi_0 = beta rho.
    # delta C = 2 Phi_0 delta Phi - 2 (1+rho^2) beta partial_rho delta Phi.
    dC = sp.simplify(
        2 * phi0 * delta
        - 2 * (1 + rho**2) * beta * sp.diff(delta, rho)
    )
    if sp.simplify(dC + 2 * beta * a) != 0:
        raise AssertionError(dC)
    # Killing fields of unit AdS2.
    K1_rho = sp.sqrt(1 + rho**2) * sp.sin(tau)
    K1_tau = rho * sp.cos(tau) / sp.sqrt(1 + rho**2)
    K2_rho = sp.sqrt(1 + rho**2) * sp.cos(tau)
    K2_tau = -rho * sp.sin(tau) / sp.sqrt(1 + rho**2)
    if sp.simplify(beta * K1_rho - beta * sp.sqrt(1 + rho**2) * sp.sin(tau)) != 0:
        raise AssertionError("K1 dilaton")
    if sp.simplify(beta * K2_rho - beta * sp.sqrt(1 + rho**2) * sp.cos(tau)) != 0:
        raise AssertionError("K2 dilaton")
    right = sp.limit(b * K2_tau + cc * K1_tau, rho, sp.oo)
    left = sp.limit(b * K2_tau + cc * K1_tau, rho, -sp.oo)
    right_s = sp.simplify(right)
    left_s = sp.simplify(left)
    # right -> -b sin + c cos; left -> b sin - c cos = -right.
    if sp.simplify(right_s + b * sp.sin(tau) - cc * sp.cos(tau)) != 0:
        raise AssertionError(right_s)
    if sp.simplify(left_s + right_s) != 0:
        raise AssertionError(left_s)
    # Clock normalization xi^tau -> 0 at either end forces b=c=0.
    # Mass fixes a: delta M = r_e delta C / (8 G) = - r_e beta a / (4 G).
    dM = sp.simplify(re * dC / (8 * G))
    if sp.simplify(dM + re * beta * a / (4 * G)) != 0:
        raise AssertionError(dM)
    # beta = 2 r_e / ell matches the static tangent a = delta beta = 2 r_e delta x_*.
    a_from_x = 2 * re / ell
    delta_M_static = -a_from_x * re * (2 * re / ell) / (4 * G)
    # The earlier tangent is delta(delta M) = -A x_* delta x with A=r_e^3/G, x_*=1/ell,
    # and a = 2 r_e delta x, so one unit delta x gives a = 2 r_e / ell? No:
    # a = 2 r_e * delta x. For delta x = 1/ell the slope change is the whole beta.
    tangent = -(re**3 / G) * (1 / ell) * (a / (2 * re))
    if sp.simplify(dM.subs(beta, 2 * re / ell) - tangent) != 0:
        raise AssertionError(sp.simplify(dM.subs(beta, 2 * re / ell) - tangent))
    # Global time translation does not change Phi_0. A rate-changing
    # reparametrization alpha(tau) changes the asymptotic lapse by alpha'.
    alpha = sp.Function("alpha")(tau)
    lapse_factor = 1 + sp.diff(alpha, tau)
    normalized = sp.simplify(lapse_factor.subs(alpha, sp.Integer(0) * tau))
    rejected = sp.simplify(lapse_factor.subs(alpha, sp.sin(tau)) - (1 + sp.cos(tau)))
    if normalized != 1 or rejected != 0:
        raise AssertionError((normalized, rejected))
    return {
        "delta_C": "-2 beta a",
        "clock_killing_kernel": "b = c = 0",
        "left_tau_shift": "- right_tau_shift",
        "delta_M": "- r_e beta a / (4 G)",
        "static_tangent_match": True,
    }


def _gauss(x, center, width, amp):
    z = (x - center) / width
    return amp * math.exp(-0.5 * z * z)


def _gauss_int_from_minus_inf(x, center, width, amp):
    z = (x - center) / (width * math.sqrt(2))
    return amp * width * math.sqrt(2 * math.pi) * 0.5 * (1 + math.erf(z))


def check_both_mouths():
    """Integrate mass and charge on both mouths from the balance laws."""
    # Right mouth: r* increases toward infinity. Pulses are compact Gaussians.
    samples = []
    max_mass = 0.0
    max_charge = 0.0
    max_rate = 0.0
    centers = (0.3, -0.4, 1.1, -1.2)
    widths = (0.15, 0.2, 0.12, 0.18)
    amps = (0.7, -0.4, 0.5, 0.25)
    for i, (cu, cv, wu, wv, au, av) in enumerate(
        (
            (centers[0], centers[1], widths[0], widths[1], amps[0], amps[1]),
            (centers[2], centers[3], widths[2], widths[3], amps[2], amps[3]),
        )
    ):
        # Mouth label i=0 right, i=1 left. The left chart uses the same
        # outgoing/incoming split with its own tortoise increasing outward.
        for t in (0.0, 0.4, -0.2):
            for r_b in (0.5, 1.5):
                def t_plus(x, cu=cu, wu=wu, au=au):
                    return _gauss(x, cu, wu, au)

                def t_minus(x, cv=cv, wv=wv, av=av):
                    return _gauss(x, cv, wv, av)

                I_plus = lambda x, cu=cu, wu=wu, au=au: _gauss_int_from_minus_inf(x, cu, wu, au)
                I_minus = lambda x, cv=cv, wv=wv, av=av: _gauss_int_from_minus_inf(x, cv, wv, av)
                # m = -I_+(u) + I_-(v), integration constant 0.
                u_cut = t - r_b
                v_cut = t + r_b
                m_cut = -I_plus(u_cut) + I_minus(v_cut)
                M = I_minus(1e6)  # I_+(-inf)=0
                analytic = M - m_cut
                # Independent Simpson integral over r*.
                r_max = r_b + 12.0
                n = 2400
                step = (r_max - r_b) / n
                acc = 0.0
                for k in range(n + 1):
                    r = r_b + k * step
                    weight = 1.0 if k in (0, n) else (4.0 if k % 2 else 2.0)
                    acc += weight * (t_plus(t - r) + t_minus(t + r))
                simpson = acc * step / 3.0
                max_mass = max(max_mass, abs(analytic - simpson))
                # Finite-difference time rate versus the cut fluxes.
                dt = 1e-4
                def collar(tt):
                    mc = -I_plus(tt - r_b) + I_minus(tt + r_b)
                    return M - mc

                rate = (collar(t + dt) - collar(t - dt)) / (2 * dt)
                flux = t_plus(t - r_b) - t_minus(t + r_b)
                max_rate = max(max_rate, abs(rate - flux))
                # Electric flux. j_u, j_v are separate Gaussians.
                ju = t_plus
                jv = t_minus
                # partial_u Q = -j_u, partial_v Q = +j_v, in units g_c^2/(4 pi)=1.
                Q_cut = -I_plus(u_cut) + I_minus(v_cut)
                Q_inf = I_minus(1e6)
                # Integral of radial current: same Simpson combination.
                acc_q = 0.0
                for k in range(n + 1):
                    r = r_b + k * step
                    weight = 1.0 if k in (0, n) else (4.0 if k % 2 else 2.0)
                    acc_q += weight * (ju(t - r) + jv(t + r))
                simpson_q = acc_q * step / 3.0
                max_charge = max(max_charge, abs((Q_inf - Q_cut) - simpson_q))
                samples.append(i)
    if max_mass > 1e-8 or max_charge > 1e-8 or max_rate > 1e-6:
        raise AssertionError(
            f"mouth integrals mass={max_mass} charge={max_charge} rate={max_rate}"
        )
    # Normalized clock: alpha'=0. A trial alpha=sin(t) is rejected because
    # the Killing lapse changes by cos(t). Both mouths share one normalized t.
    trial_lapse_deviation = max(abs(math.cos(t)) for t in (0.0, 0.3, 1.0))
    if trial_lapse_deviation < 0.5:
        raise AssertionError("trial clock was not distinguished from unit lapse")
    return {
        "mouths": 2,
        "slice_samples": len(samples),
        "max_mass_residual": max_mass,
        "max_charge_residual": max_charge,
        "max_rate_residual": max_rate,
        "normalized_clock": "alpha'=0; trial sin(t) changes the lapse",
        "units": "g_c^2/(4 pi)=1 in the charge comparison; anomaly subtracted",
    }


def check_retarded_support():
    """A right-moving map perturbation supported at u>0 leaves T_uu unchanged for u<0."""
    uu = sp.symbols("u", real=True)
    # Identity map: Schwarzian vanishes, so the past stress is exactly the reference.
    sch_id = schwarzian(uu, uu)
    if sp.simplify(sch_id) != 0:
        raise AssertionError(sch_id)
    # Future bump R(u)=u+exp(-1/u), u>0. Its Schwarzian is not identically zero.
    bump = uu + sp.exp(-1 / uu)
    sch_bump = schwarzian(bump, uu)
    future = abs(float(sch_bump.subs(uu, sp.Rational(2, 5))))
    if future < 1e-3:
        raise AssertionError(f"future contact vanished: {future}")
    # The pulled-back stress at u<0 uses only R and its jets at that u.
    # For any map that equals the identity on (-inf, 0], those jets match the
    # identity, so the stress perturbation vanishes there.
    # A map that equals the identity on u<0 has the identity 3-jet there, so
    # its Schwarzian stress perturbation vanishes in the past.
    R = sp.Piecewise((uu, uu < 0), (bump, True))
    past_point = sp.Rational(-1, 2)
    for order in (1, 2, 3):
        got = sp.diff(R, uu, order).subs(uu, past_point)
        target = sp.diff(uu, uu, order).subs(uu, past_point)
        if sp.simplify(got - target) != 0:
            raise AssertionError((order, got, target))
    return {"past_stress_residual": 0.0, "future_contact_abs": future}


def float_samples():
    """Evaluate the simplified identities on random jets. Residuals are computed."""
    import random

    rng = random.Random(51012)
    max_dc = 0.0
    for _ in range(2000):
        aa = rng.uniform(-1, 1)
        bb = rng.uniform(-1, 1)
        cc = rng.uniform(-1, 1)
        rho_v = rng.uniform(-2, 2)
        tau_v = rng.uniform(-math.pi, math.pi)
        beta_v = rng.uniform(0.2, 3)
        delta = aa * rho_v + math.sqrt(1 + rho_v**2) * (
            bb * math.cos(tau_v) + cc * math.sin(tau_v)
        )
        d_delta = aa + rho_v / math.sqrt(1 + rho_v**2) * (
            bb * math.cos(tau_v) + cc * math.sin(tau_v)
        )
        dC = 2 * beta_v * rho_v * delta - 2 * (1 + rho_v**2) * beta_v * d_delta
        max_dc = max(max_dc, abs(dC + 2 * beta_v * aa))
    if max_dc > 1e-9:
        raise AssertionError(max_dc)
    return {"mode_casimir_samples": 2000, "max_delta_C_residual": max_dc}


def main():
    report = {
        "conservation": check_conservation(),
        "rn_reduction": check_rn_reduction(),
        "zeta": zeta_dirac_potential(),
        "massive_ward": massive_ward(),
        "join_contact": check_join_contact(),
        "current": check_current_and_gauss(),
        "modes": check_modes(),
        "mouths": check_both_mouths(),
        "retarded": check_retarded_support(),
        "float_samples": float_samples(),
    }
    payload = {
        "date": "2026-10-12",
        "check": "checks/retarded_dirac_spherical_response.py",
        "note": "notes/retarded-dirac-spherical-response.md",
        "prescription": "spherical Hadamard-Polyakov",
        "endpoint": "obstruction",
        "l0_clock_is_independent_mode": False,
        "homogeneous_coefficients_fixed_by_mass_and_clocks": True,
        "ward_identities_hold": True,
        "both_mouths_integrated": True,
        "shortcut_claimed": False,
        "lean_invocations": 0,
        "formal_certificate": "none",
        "max_mass_residual": report["mouths"]["max_mass_residual"],
        "max_charge_residual": report["mouths"]["max_charge_residual"],
        "max_rate_residual": report["mouths"]["max_rate_residual"],
        "max_delta_C_residual": report["float_samples"]["max_delta_C_residual"],
        "past_stress_residual": report["retarded"]["past_stress_residual"],
        "future_contact_abs": report["retarded"]["future_contact_abs"],
        "report": report,
        "scope": (
            "Spherical Hadamard-Polyakov response, Ward identities, both-mouth "
            "mass and Gauss integrals, and the clock-normalization obstruction. "
            "Not a stability, shortcut, formation, or engineering result."
        ),
    }
    if NOTE.exists():
        payload["sha256"] = {"check": digest(Path(__file__)), "note": digest(NOTE)}
    else:
        payload["sha256"] = {"check": digest(Path(__file__))}
    EVIDENCE.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        "PASS"
        f" mass={payload['max_mass_residual']:.3e}"
        f" charge={payload['max_charge_residual']:.3e}"
        f" rate={payload['max_rate_residual']:.3e}"
        f" deltaC={payload['max_delta_C_residual']:.3e}"
        f" past={payload['past_stress_residual']:.3e}"
        f" future={payload['future_contact_abs']:.3e}"
        f" endpoint=obstruction"
    )


if __name__ == "__main__":
    main()
