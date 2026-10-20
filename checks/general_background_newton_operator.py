"""S162 operator identities and diagnostics; no physical inverse certificate."""
from __future__ import annotations
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import mpmath as mp
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "general-background-newton-operator"
OWN = ["checks/general_background_newton_operator.py", "notes/"+STEM+".md"]
PREDECESSORS = [
    ("stationary-global-weighted-map", 50),
    ("stationary-adjoint-phase", 29),
    ("stationary-scalar-outgoing-inverse", 107),
    ("stationary-scalar-phase-moment", 80),
    ("stationary-shape-logarithmic-response", 45),
    ("stationary-disk-and-normal-inverse", 61),
    ("comparison-path-complete-equation-residual", 329),
]
PROTOCOL = {
    "digits": 80, "radii": ["1/5", "2/3", "7/4"],
    "contour_nodes": [64, 128], "contour_radius": ".005",
    "contour_tolerance": "1e-45", "real_steps": [".0001", ".00005"],
    "matrix_steps": [".0002", ".0001"], "heat_times": [".07", ".4"],
    "difference_tolerance": "1e-5", "refinement_ratio": ".35",
    "roundoff_floor": "1e-50", "fixtures_are_physical_states": False,
}
ROWS, NUMERICS = [], []
R = S.Rational


def check(name, passed, detail=""):
    assert name not in {r["name"] for r in ROWS}
    row = {"name": name, "passed": bool(passed), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not passed:
        raise AssertionError(name)


def zero(name, expr):
    if isinstance(expr, S.MatrixBase):
        check(name, all(S.factor(S.together(x)) == 0 for x in expr))
    else:
        check(name, S.factor(S.together(expr)) == 0)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_text(encoding="utf-8"))
    return data.get("runs", [data])[-1]


def inputs(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def names():
    out = set(OWN)
    for stem, _ in PREDECESSORS:
        out.update([f"notes/{stem}.md", f"evidence/{stem}.json"])
        out.update(inputs(latest(stem)))
    out.add("notes/exact-spherical-em-reduction.md")
    return sorted(out)


def predecessors():
    for stem, count in PREDECESSORS:
        r = latest(stem)
        check(stem+"_receipt", r["status"] == "passed" and len(r["checks"]) == count
              and all(x["passed"] for x in r["checks"]) and inputs(r) == inputs(r, True))
        for p, h in inputs(r).items():
            check(stem+"_input_"+p, sha256((ROOT/p).read_bytes()).hexdigest() == h)


def algebra():
    ep, lam, rr, C, Eu, M, BR = S.symbols("epsilon lambda R C2 Euler M BoxR")
    s, Bs, BBs, JHs, grs = S.symbols("s Boxs Box2s JHesss gradRgrads")
    lam0 = lam+R(11, 6)*ep*lam**2
    an = R(3, 2)*C-R(11, 12)*Eu-20*M+BR
    E = rr+12*lam0+ep*an
    dR = -6*Bs-2*rr*s
    dan = (R(3, 2)*(-4*s*C)
           -R(11, 12)*(-4*s*Eu+8*JHs-2*rr*Bs)
           +80*s*M-6*BBs-2*rr*Bs-2*grs-4*s*BR)
    dE = dR+ep*dan
    full = (6*ep*BBs+(6+ep*rr/6)*Bs+R(22, 3)*ep*JHs
            +2*ep*grs-(2*rr+48*lam0)*s)
    zero("off_shell_conformal_hessian", -(dE+4*s*E)-full)
    zero("contracted_Bianchi_divergence", ep/6+R(22, 3)*ep/4-2*ep)
    check("density_contact_not_identically_zero", S.expand(-dE-full) != 0)
    p, Z = S.symbols("T Z")
    lhs = (6*ep*p**2-(6+ep*rr/6)*p-(2*rr+48*lam0))/(6*ep)
    rhs = (p+4*lam)*(p-(1/ep+R(11, 3)*lam))-Z*p/36-Z/(3*ep)
    zero("general_factorization", lhs.subs(rr, Z-12*lam)-rhs)
    zero("tracefree_Ricci_coefficient", R(22, 3)/6-R(11, 9))

    ell, U, rUp = S.symbols("ell U rUp")
    El = rr+(1-ell)*U+12*ell*lam0+ell*ep*an
    dEl = dR+(1-ell)*rUp*s+ell*ep*dan
    Hl = (6*ell*ep*BBs+(6+ell*ep*rr/6)*Bs+R(22, 3)*ell*ep*JHs
          +2*ell*ep*grs+(-2*rr-(1-ell)*(rUp+4*U)-48*ell*lam0)*s)
    zero("bookkeeping_full_hessian", -(dEl+4*s*El)-Hl)
    zero("bookkeeping_seed_specialization", Hl.subs({ell: 0, rr: -U})
         -(6*Bs-(2*U+rUp)*s))

    x, y = S.symbols("Phi V")
    density = S.exp(4*x+2*y)
    action = density*(x**3+x*y+y**4)
    F = S.Matrix([S.diff(action, z)/density for z in (x, y)])
    DF = F.jacobian((x, y))
    A = DF+F*S.Matrix([[4, 2]])
    zero("fixed_density_Hessian_symmetry", A-A.T)
    zero("fixed_reference_gradient_derivative", (density*F).jacobian((x, y))-density*A)
    zero("moving_volume_adjoint_contact", DF.T-(A-S.Matrix([4, 2])*F.T))
    check("off_shell_DF_asymmetry_detected", DF != DF.T)
    zero("correlated_mass_density_cancellation", 4+2*(-2))
    er, et, ea = S.symbols("E_rho E_varphi E_angle")
    zero("independent_radius_row", (er+et+2*ea)-(er+et)-2*ea)
    # An injective index-minus-one map, not a model of the exterior.
    Am = S.Matrix([[1, 2], [0, 3], [0, 0]])
    L = S.Matrix([[0, 0, 1]])
    q = S.Matrix([2, -1, 1])
    Pi = S.eye(3)-q*L
    border = Am.row_join(q)
    f = S.Matrix(S.symbols("f0:3"))
    uv = border.inv()*f
    zero("bordered_eta_is_signed_condition", uv[2]-(L*f)[0])
    zero("bordered_projection_reconstructs", Am*uv[:2, 0]-Pi*f)
    zero("projection_idempotent", Pi*Pi-Pi)
    check("augmented_nonphysical_source_detected", uv[2].subs({f[2]: 1}) != 0)


def warped_fixture():
    x, u = S.symbols("rho u", positive=True)
    a = x*(1+x**2/7+x**4/101)
    rad = 2+x**2/5+x**4/103
    tau = R(1, 13)+x**2/17+x**4/109+x**6/601
    n, au, ru = S.exp(u*tau), S.exp(u*tau)*a, S.exp(u*tau)*rad
    proper = lambda f: S.diff(f, x)/n
    # Independent metric construction, including the varying radial lapse.
    k01 = S.simplify(-proper(proper(au))/au)
    k02 = S.simplify(-proper(au)*proper(ru)/(au*ru))
    k12 = S.simplify(-proper(proper(ru))/ru)
    k23 = S.simplify((1-proper(ru)**2)/ru**2)
    ric0, ric1, ric2 = k01+2*k02, k01+2*k12, k02+k12+k23
    scal = 2*k01+4*k02+4*k12+2*k23
    ric_sq = ric0**2+ric1**2+2*ric2**2
    rm_sq = 4*(k01**2+2*k02**2+2*k12**2+k23**2)
    weyl_sq, euler = rm_sq-2*ric_sq+scal**2/3, rm_sq-4*ric_sq+scal**2
    H = S.diff(a, x)/a+2*S.diff(rad, x)/rad
    box = lambda f: S.diff(f, x, 2)+H*S.diff(f, x)
    boxu = lambda f: (S.diff(f, x, 2)+(H+2*u*S.diff(tau, x))*S.diff(f, x))/n**2
    ep, lam, charge = R(3, 50), R(2, 11), R(3, 7)
    lam0 = lam+R(11, 6)*ep*lam**2
    trace = scal+12*lam0+ep*(R(3, 2)*weyl_sq-R(11, 12)*euler
                             -20*charge**2/(2*ru**4)+boxu(scal))
    dens = S.exp(4*u*tau)*trace
    r0 = scal.subs(u, 0)
    JH = ((ric1.subs(u, 0)-r0/4)*S.diff(tau, x, 2)
          +((ric0.subs(u, 0)-r0/4)*S.diff(a, x)/a
            +2*(ric2.subs(u, 0)-r0/4)*S.diff(rad, x)/rad)*S.diff(tau, x))
    pred = (6*ep*box(box(tau))+(6+ep*r0/6)*box(tau)
            +R(22, 3)*ep*JH+2*ep*S.diff(r0, x)*S.diff(tau, x)
            -(2*r0+48*lam0)*tau)
    fun = S.lambdify((x, u), dens, "mpmath", cse=True)
    answer = S.lambdify(x, pred, "mpmath", cse=True)
    omitted = S.lambdify(x, 4*tau*trace.subs(u, 0), "mpmath", cse=True)
    for xs in PROTOCOL["radii"]:
        xx = mp.mpf(xs)
        target = answer(xx)
        scale = 1+abs(target)
        contours = []
        for nodes in PROTOCOL["contour_nodes"]:
            radius = mp.mpf(PROTOCOL["contour_radius"])
            val = -sum(fun(xx, radius*mp.exp(2j*mp.pi*j/nodes))*mp.exp(-2j*mp.pi*j/nodes)
                       for j in range(nodes))/(nodes*radius)
            error = abs(val-target)
            check(f"warped_{xs}_contour_{nodes}", error < mp.mpf("1e-45")*scale, mp.nstr(error, 20))
            contours.append(mp.nstr(val, 30))
        errs = []
        for hs in PROTOCOL["real_steps"]:
            h = mp.mpf(hs)
            error = abs(-(fun(xx, h)-fun(xx, -h))/(2*h)-target)
            errs.append(error)
            check(f"warped_{xs}_difference_{hs}", error < mp.mpf("1e-5")*scale, mp.nstr(error, 20))
        check(f"warped_{xs}_refinement", errs[1] < mp.mpf(".35")*errs[0] or max(errs) < mp.mpf("1e-50"))
        check(f"warped_{xs}_omitted_density_detected", abs(omitted(xx)) > mp.mpf("1e-4"))
        NUMERICS.append({"fixture": xs, "target": mp.nstr(target, 30),
                         "contours": contours, "real_errors": [mp.nstr(e, 20) for e in errs]})


def matrix_fixture():
    D = S.Matrix([[2, R(1, 3), 0], [0, 3, R(2, 5)], [R(1, 7), 0, 4]])
    v, w = S.diag(R(1, 5), R(-1, 4), R(2, 7)), S.diag(R(-1, 3), R(1, 6), R(1, 8))
    d1, d2 = -(v*D+D*v)/2, -(w*D+D*w)/2
    d12 = (v*w*D+v*D*w+w*D*v+D*v*w)/4
    Q = D.T*D
    q1, q2 = d1.T*D+D.T*d1, d2.T*D+D.T*d2
    q12 = d12.T*D+D.T*d12+d1.T*d2+d2.T*d1
    # Entrywise exact differentiation of the conjugated matrix, independently.
    s, t = S.symbols("s t")
    family = S.diag(*[S.exp(-(s*v[i, i]+t*w[i, i])/2) for i in range(3)])
    ds = family*D*family
    qs = ds.T*ds
    zero("Dirac_first_vertex", ds.diff(s).subs({s: 0, t: 0})-d1)
    zero("Dirac_contact_vertex", ds.diff(s, t).subs({s: 0, t: 0})-d12)
    zero("square_first_vertex", qs.diff(s).subs({s: 0, t: 0})-q1)
    zero("square_contact_vertex", qs.diff(s, t).subs({s: 0, t: 0})-q12)
    bubble = S.trace(Q.inv()*q1*Q.inv()*q2)/2
    contact = -S.trace(Q.inv()*q12)/2
    zero("finite_logdet_contact_cancellation", bubble+contact)
    check("finite_logdet_omitted_contact_detected", bubble != 0)
    check("Dirac_fixture_nonnormal", D.T*D != D*D.T)
    check("multipliers_commute_but_D_does_not", v*w == w*v and v*D != D*v)

    conv = lambda a: mp.matrix([[mp.mpf(str(x)) for x in row] for row in a.tolist()])
    dd, qq, qa, qb, qab = map(conv, (D, Q, q1, q2, q12))
    vv, ww = [mp.mpf(str(v[i, i])) for i in range(3)], [mp.mpf(str(w[i, i])) for i in range(3)]
    eig, U = mp.eigh(qq)
    a, b, ab = U.T*qa*U, U.T*qb*U, U.T*qab*U
    tr = lambda m: sum(m[i, i] for i in range(3))
    for ts in PROTOCOL["heat_times"]:
        heat_t = mp.mpf(ts)
        ct = -heat_t*sum(ab[i, i]*mp.exp(-heat_t*eig[i]) for i in range(3))
        bt = mp.mpf(0)
        for i in range(3):
            for j in range(3):
                integ = (mp.exp(-heat_t*eig[i]) if i == j else
                         (mp.exp(-heat_t*eig[i])-mp.exp(-heat_t*eig[j]))
                         /(heat_t*(eig[j]-eig[i])))
                bt += heat_t**2*a[i, j]*b[j, i]*integ
        target = ct+bt

        def heat(s, t):
            mm = mp.diag([mp.exp(-(s*vv[i]+t*ww[i])/2) for i in range(3)])
            dv = mm*dd*mm
            return tr(mp.expm(-heat_t*dv.T*dv))

        errors = []
        for hs in PROTOCOL["matrix_steps"]:
            h = mp.mpf(hs)
            val = (heat(h, h)-heat(h, -h)-heat(-h, h)+heat(-h, -h))/(4*h*h)
            error = abs(val-target)
            errors.append(error)
            check(f"matrix_heat_{ts}_{hs}", error < mp.mpf("1e-5")*(1+abs(target)), mp.nstr(error, 20))
        check(f"matrix_heat_{ts}_refinement",
              errors[1] < mp.mpf(".35")*errors[0] or max(errors) < mp.mpf("1e-50"))
        check(f"matrix_heat_{ts}_omitted_contact_detected", abs(ct) > mp.mpf("1e-4"))
        NUMERICS.append({"heat_time": ts, "target": mp.nstr(target, 30),
                         "contact": mp.nstr(ct, 30), "real_errors": [mp.nstr(e, 20) for e in errors]})


def main():
    mp.mp.dps = PROTOCOL["digits"]
    paths = names()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    before, error = hashes(), None
    print(json.dumps({"event": "start", "story": "S162", "input_count": len(paths)}), flush=True)
    try:
        predecessors()
        algebra()
        warped_fixture()
        matrix_fixture()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(r["passed"] for r in ROWS)
    run = {
        "story": "S162", "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "mpmath": mp.__version__, "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "numerics": NUMERICS,
        "exception": error,
        "scope": "general-background retained Hessian and compatible range formulation",
        "unresolved": ["quantitative compatible inverse", "nonlinear family and map",
                       "signed compatibility", "coupled solution", "physical completion"],
    }
    if error:
        run["failed_checker_source_utf8"] = (ROOT/OWN[0]).read_text(encoding="utf-8")
        run["failed_note_source_utf8"] = (ROOT/OWN[1]).read_text(encoding="utf-8")
    dest = ROOT/"evidence"/(STEM+".json")
    data = json.loads(dest.read_text(encoding="utf-8")) if dest.exists() else {"runs": []}
    data["runs"].append(run)
    dest.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": run["status"], "checks": len(ROWS), "unchanged_inputs": before == after}), flush=True)
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
