"""Moving angular identification and its specified classical loaded response."""
import hashlib
import json
import platform
import runpy
from pathlib import Path

import numpy as np
import scipy
import sympy as sy
from scipy.integrate import quad, solve_bvp, solve_ivp
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
basis = runpy.run_path(str(ROOT/'checks/gluing_support_compatibility.py'))['basis']
ETA = np.diag([-1., 1., 1., 1.])


def rotation(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.]])


def boost(v):
    v = np.array(v)
    g = 1/np.sqrt(1-v@v)
    out = np.eye(4)
    out[0, 0] = g
    out[0, 1:] = out[1:, 0] = g*v
    out[1:, 1:] += g*g/(g+1)*np.outer(v, v)
    return out


def exact():
    t, om, omr = sy.symbols('t om omr', real=True)
    c, s = sy.cos(om*t), sy.sin(om*t)
    Rot = sy.Matrix([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    F0 = sy.diag(-1, 1, 1)
    F = Rot*F0*Rot.T
    Z = sy.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 0]])
    items = {
        'reflection velocity': sy.diff(F, t)*F.T-2*om*Z,
        'corotating reflection': Rot.T*F*Rot-F0,
        'axial pullback': F0.det()*F0.T*sy.Matrix([0, 0, -om])+sy.Matrix([0, 0, -om]),
        'relative proper rotation': sy.trigsimp(F*F0)-sy.Matrix([
            [sy.cos(2*om*t), -sy.sin(2*om*t), 0],
            [sy.sin(2*om*t), sy.cos(2*om*t), 0], [0, 0, 1]]),
    }
    for name, value in items.items():
        assert value.applyfunc(sy.trigsimp) == sy.zeros(*value.shape), name
    # A common shift w,V -> w-Omega,V-Omega leaves Phi and all curls unchanged.
    w, V, Omega = sy.symbols('w V Omega')
    assert sy.expand((V-Omega)-(w-Omega)-(V-w)) == 0
    return list(items)+['common Maxwell/metric frame shift']


def fermi_walker():
    rows = []
    for v in [0.02, 0.2, 0.6, 0.8]:
        g = 1/np.sqrt(1-v*v)
        for end in [1.7, 2*np.pi]:
            def velocity(t):
                return v*np.array([-np.sin(t), np.cos(t), 0.])
            def rhs(t, y):
                u = g*np.r_[1., velocity(t)]
                a = g*g*np.array([0., -v*np.cos(t), -v*np.sin(t), 0.])
                generator = (np.outer(u, a@ETA)-np.outer(a, u@ETA))/g
                return (generator@y.reshape(4, 3)).ravel()
            init = boost(velocity(0))[:, 1:]
            errors = []
            norms = []
            for tol in [2e-8, 2e-11]:
                sol = solve_ivp(rhs, [0, end], init.ravel(), method='DOP853', rtol=tol, atol=tol/10)
                assert sol.success
                frame = sol.y[:, -1].reshape(4, 3)
                rest = boost(-velocity(end))@frame
                expected = np.vstack((np.zeros(3), rotation(-(g-1)*end)))
                errors.append(float(np.max(abs(rest-expected))))
                norms.append(float(np.max(abs(frame.T@ETA@frame-np.eye(3)))))
            assert errors[-1] < 1e-9 and errors[-1] < errors[0]
            rows.append(dict(speed=v, orbital_angle=end, thomas_angle=-(g-1)*end,
                             errors=errors, orthonormality_errors=norms))
    return rows


def wilson_checks():
    sx = np.array([[0, 1], [1, 0]], complex)/2
    sy = np.array([[0, -1j], [1j, 0]], complex)/2
    sz = np.diag([1., -1.])/2
    rows = []
    for t in [0.1, 0.7, 1.3, 2.2]:
        def at(x):
            return x*x*sz+np.sin(t)*sx
        def rhs(x, y):
            W, dW, integral = y.reshape(3, 2, 2)
            Ax = t*sx+x*sy
            F = sx-2*x*sz+1j*(at(x)@Ax-Ax@at(x))
            return np.stack((-1j*Ax@W, -1j*sx@W-1j*Ax@dW, W.conj().T@F@W)).ravel()
        errors = []
        initial = np.stack((np.eye(2), np.zeros((2, 2)), np.zeros((2, 2)))).astype(complex)
        for tol in [1e-8, 1e-11]:
            sol = solve_ivp(rhs, [0, 1], initial.ravel(), method='DOP853', rtol=tol, atol=tol/10)
            assert sol.success
            W, dW, integ = sol.y[:, -1].reshape(3, 2, 2)
            residual = dW+1j*at(1)@W-1j*W@at(0)+1j*W@integ
            errors.append(float(np.linalg.norm(residual)))
        assert errors[-1] < 1e-10 and errors[-1] < errors[0]
        rows.append(dict(time=t, covariant_identity_errors=errors))
    return rows


def matching():
    rows = []
    # omega_hat=P Omega is set to one by linearity. e/omega_hat=chi_b.
    for lam in [.01, .003, .001, .0003, .0001]:
        R = lam**(-1/3)
        z0, rs = lam*R, 1+lam*R
        f = z0*z0/(rs*rs)
        H, g, gp = R*(1+R*R), 1+3*R*R, 6*R
        u0, u0p, u1, u1p = basis(rs)
        mat = np.array([[u1, -g], [f*u1p, -lam*(1+R*R)*gp]])
        bE, kE = np.linalg.solve(mat, [-u0, -f*u0p])
        bJ, kJ = np.linalg.solve(mat, [2*(1/rs-1), -2*f/rs**2])
        for rout, chi_b in [(3., .03), (21., .001), (101., .0001)]:
            def frame(e, J, b, k):
                def integrand(x):
                    z = np.exp(x)
                    r = 1+z
                    u = e*basis(r)[0]+b*basis(r)[2]-2*J/r
                    return z*(-6*J-2*u)/r**4
                return lam*(-2*k*H-2*J*R)+quad(integrand, np.log(z0), np.log(rout-1),
                                             epsabs=2e-12, epsrel=2e-12)[0]
            wE = frame(1, 0, bE, kE)
            wJ = frame(0, 1, bJ, kJ)
            e = chi_b
            J = (-1-e*wE)/wJ
            b, k = e*bE+J*bJ, e*kE+J*kJ
            assert abs(frame(e, J, b, k)+1) < 1e-11

            def rhs(x, y, param):
                z = np.exp(x)
                r = 1+z
                return np.vstack((y[1], (2+4/r**2)*y[0]-(2/r-1)*y[1]+12*param[0]/r**2,
                                  z*(-6*param[0]-2*y[0])/r**4, z*y[0]/r**2))

            def bc(left, right, param):
                j = param[0]
                kin = (left[0]+2*j)/g
                a0, a0p, a1, a1p = basis(rout)
                return np.array([left[1]-z0/f*lam*(1+R*R)*gp*kin,
                                 left[2]-lam*(-2*kin*H-2*j*R),
                                 left[3]-lam*(kin*H-2*j*R),
                                 right[1]/(rout-1)-e*a0p-(right[0]-e*a0+2*j/rout)*a1p/a1-2*j/rout**2,
                                 right[2]+1])
            refs = []
            for size, tol in [(120, 3e-7), (360, 3e-9)]:
                mesh = np.linspace(np.log(z0), np.log(rout-1), size)
                rr = 1+np.exp(mesh)
                guess = np.vstack((e*basis(rr)[0], e*(rr-1)*basis(rr)[1], np.zeros(size), np.zeros(size)))
                sol = solve_bvp(rhs, bc, mesh, guess, p=[0.], tol=tol, max_nodes=18000)
                assert sol.success, sol.message
                xs = np.linspace(mesh[0], mesh[-1], 901)
                rr = 1+np.exp(xs)
                vals = sol.sol(xs)
                expected = e*basis(rr)[0]+b*basis(rr)[2]-2*J/rr
                error = float(np.max(abs(vals[0]-expected)/np.maximum(1, abs(expected))))
                jerror = abs(sol.p[0]-J)
                assert max(error, jerror) < 5e-6
                # Phi=V-w must also agree with the flux derivative reconstruction.
                invariant = ((rr-1)/rr)**2*vals[1]/(rr-1)/2
                perr = float(np.max(abs(vals[3]-vals[2]-invariant)))
                assert perr < 2e-6
                refs.append(dict(nodes=int(sol.x.size), tolerance=tol, flux_error=error, charge_error=jerror,
                                 invariant_potential_error=perr,
                                 max_V_over_omega_hat=float(np.max(abs(vals[3]))),
                                 max_w_over_omega_hat=float(np.max(abs(vals[2])))))
            assert refs[-1]['charge_error'] < refs[0]['charge_error']
            w0 = 2/rout-(2*rout-1)/rout**4
            wkn = (2*rout-1)/rout**4
            rows.append(dict(lam=lam, R=R, r_out=rout, chi_b=chi_b, J_over_omega_hat=J,
                             asymptotic_ratio=(1-chi_b*(1-w0))/(1-wkn), refinements=refs))
    return rows


def spin_and_scales():
    rows = []
    for q in [2, 3, 6, 11, 22]:
        m = np.arange(q)-(q-1)/2
        Jz = np.diag(m)
        angles = [.01, .3, 1.2]
        for angle in angles:
            U = expm(-2j*angle*Jz)
            error = np.linalg.norm(U-np.eye(q), 2)
            assert error <= min(2, (q-1)*abs(angle))+1e-13
            rows.append(dict(q=q, orbital_angle=angle, norm=float(error)))
    scales = []
    for q in [1e3, 1e6, 1e9, 1e12, 1e18, 1e30]:
        delta = np.log(q)**-2
        lam = 1/(16*np.pi*q**3)
        chi = lam/delta
        sigma = 2*np.sqrt(lam/delta**3)
        K = np.pi*q**4
        omega_hat = lam*sigma
        e = chi*omega_hat
        assert np.isclose(K*omega_hat, q*sigma/16, rtol=1e-14, atol=0)
        scales.append(dict(q=q, delta=delta, Omega_ell=sigma, q_Omega_damping_time=q*sigma/delta,
                           classical_KJ_scale=K*omega_hat, electric_e=e, omega_hat=omega_hat,
                           worst_gap_proxy=q*q*sigma*delta))
    return dict(spin=rows, scales=scales)


def digest(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def main():
    result = dict(status='passed_moving_reflection_and_classical_load_not_quantum_state',
                  python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, sympy=sy.__version__,
                  exact=exact(), fermi_walker=fermi_walker(), wilson=wilson_checks(), matched=matching(),
                  spin_and_scales=spin_and_scales(), script_sha256=digest(Path(__file__)),
                  input_sha256={f: digest(ROOT/f) for f in ['notes/moving-reflection-load.md',
                    'checks/gluing_support_compatibility.py', 'notes/gluing-support-compatibility.md',
                    'notes/rotational-current-response.md', 'notes/rotating-dirac-and-state.md',
                    'notes/binary-mouth-drive.md', 'notes/weak-coupling-window.md']},
                  moving_physical_return_proved=False, physical_quantum_state_selected=False,
                  full_stress_bound=False, smooth_semiclassical_junction=False,
                  novelty_established=False, lean_invocations=0, goal_complete=False)
    (ROOT/'evidence/moving-reflection-load.json').write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(dict(status=result['status'], matching_cases=len(result['matched']),
                         max_refined_charge_error=max(x['refinements'][-1]['charge_error'] for x in result['matched']),
                         max_FW_error=max(x['errors'][-1] for x in result['fermi_walker']))))


if __name__ == '__main__':
    main()
