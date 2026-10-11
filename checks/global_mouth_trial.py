"""Check the global tube's two spherical traces; no PDE scattering simulation."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import roots_legendre
import sympy as sp


def verify():
    r, z, u = sp.symbols('r z u', real=True)
    rm = sp.sqrt(r*r+(z+sp.Rational(1, 2))**2)
    rp = sp.sqrt(r*r+(z-sp.Rational(1, 2))**2)
    flux = ((z+sp.Rational(1, 2))/rm-(z-sp.Rational(1, 2))/rp)/2
    mass = (flux-u)/r
    field = sp.Matrix([r/rm**3-r/rp**3,
                       (z+sp.Rational(1, 2))/rm**3-(z-sp.Rational(1, 2))/rp**3])/2
    bf = sp.lambdify((r, z), field, 'numpy', cse=True)
    jf = sp.lambdify((r, z), field.jacobian([r, z]), 'numpy', cse=True)
    hf = sp.lambdify((r, z, u), sp.hessian(mass, [r, z]), 'numpy', cse=True)
    uf = sp.lambdify((r, z), flux, 'numpy', cse=True)

    a = sp.symbols('a', real=True)
    rot = lambda x: sp.Matrix([[sp.cos(x/2), -sp.sin(x/2)],
                              [sp.sin(x/2), sp.cos(x/2)]])
    cc = sp.symbols('c', positive=True)
    beta, theta = sp.symbols('beta theta', real=True)
    etheta = sp.Matrix([sp.cos(theta), -sp.sin(theta)])
    exact = {
        'source_spin_rotation': sp.simplify(rot(a).T*rot(a)) == sp.eye(2),
        'receiver_desired_phase': sp.simplify(rot(-sp.pi-a)[1, 0]+sp.cos(a/2)) == 0,
        'receiver_wrong_spin': sp.simplify(rot(-sp.pi-a)[0, 0]+sp.sin(a/2)) == 0,
        'flux_difference': sp.simplify((1+cc)/(2*cc)-(1-cc)/(2*cc)) == 1,
        'normalized_wrong_spin': sp.simplify(((1-cc)/(2*cc))/((1+cc)/(2*cc))-(1-cc)/(1+cc)) == 0,
        'sphere_normal_coordinate_derivative': sp.trigsimp(sp.Matrix([-sp.cos(beta), sp.sin(beta)]).dot(etheta)+sp.cos(beta-theta)) == 0,
        'sphere_longitudinal_coordinate_derivative': sp.trigsimp(sp.Matrix([sp.sin(beta), sp.cos(beta)]).dot(etheta)-sp.sin(beta-theta)) == 0,
        'gaussian_overlap_expansion': sp.series(sp.sqrt(2*sp.sqrt(sp.cos(a))/(1+sp.cos(a))), a, 0, 6).removeO() == 1-a**4/64,
    }
    assert all(exact.values())
    # A finite geometry bound, followed by a large-q existence statement in the note.
    eps = Fraction(1, 200)
    k = eps**2/(1-eps)**2
    tau_bound = k/(2*(1-k))
    shape_bound = (k/(1-k))**2
    assert tau_bound == Fraction(1, 79200)
    assert shape_bound == Fraction(1, 1568160000)
    assert tau_bound+shape_bound < Fraction(1, 50000)
    assert shape_bound**2/2 < Fraction(1, 10**18)

    def geom(points, U):
        R, Z = points
        B = np.asarray(bf(R, Z), float).reshape(2, -1)
        J = np.asarray(jf(R, Z), float).reshape(2, 2, -1)
        H = np.asarray(hf(R, Z, U), float).reshape(2, 2, -1)
        b = np.linalg.norm(B, axis=0)
        t = B/b
        n = np.array([-t[1], t[0]])
        bs = np.einsum('in,ijn,jn->n', t, J, t)
        curvature = np.einsum('in,ijn,jn->n', n, J, t)/b
        c2 = np.einsum('in,ijn,jn->n', n, H, n)
        return b, t, n, bs/(4*b), curvature, c2

    cases = []
    geometry_cases = []
    for rho in [.02, .05, .1]:
        for U in [.25, .5, .75]:
            th0 = brentq(lambda th: float(uf(rho*np.sin(th), -.5+rho*np.cos(th)))-U,
                         .01, np.pi-.01, xtol=5e-15)
            p0 = np.array([rho*np.sin(th0), -.5+rho*np.cos(th0)])
            b0, t0, n0, *_ = geom(p0.reshape(2, 1), U)
            t0 = t0[:, 0]
            er = np.array([np.sin(th0), np.cos(th0)])
            cosine = float(t0@er)
            signed_angle = float(np.arctan2(t0[0], t0[1])-th0)
            tau = np.tan(signed_angle/2)
            overlap = np.sqrt(2*np.sqrt(cosine)/(1+cosine))
            A = float(rho*rho*b0[0]*cosine)
            normalization = np.sqrt(2*cosine/(1+cosine))
            assert cosine > .999

            def rhs(_, point):
                B = np.asarray(bf(*point), float).ravel()
                return B/np.linalg.norm(B)

            def trajectories(tolerance):
                opts = dict(rtol=tolerance, atol=tolerance*rho*.03, dense_output=True,
                            max_step=rho/30)
                forward = solve_ivp(rhs, [0, rho*.55], p0, **opts)
                backward = solve_ivp(rhs, [0, -rho*.55], p0, **opts)
                assert forward.success and backward.success

                def curve(ss):
                    ss = np.atleast_1d(ss)
                    points = np.empty((2, len(ss)))
                    pos = ss >= 0
                    if np.any(pos):
                        points[:, pos] = forward.sol(ss[pos])
                    if np.any(~pos):
                        points[:, ~pos] = backward.sol(ss[~pos])
                    return points
                return curve

            curve = trajectories(2e-13)
            fine_curve = trajectories(3e-14)
            ss = np.linspace(-rho*.5, rho*.5, 101)
            refinement = float(np.max(np.linalg.norm(curve(ss)-fine_curve(ss), axis=0))/rho)
            assert refinement < 3e-11
            assert np.max(abs(uf(*fine_curve(ss))-U)) < 5e-12
            reflected = fine_curve(ss)*np.array([1., -1.])[:, None]
            left_geometry = geom(fine_curve(ss), U)
            right_geometry = geom(reflected, U)
            symmetry_error = max(
                float(np.max(abs(right_geometry[idx]-sign*left_geometry[idx])))/max(1., float(np.max(abs(left_geometry[idx]))))
                for idx, sign in [(0, 1), (3, -1), (4, 1), (5, 1)])
            assert symmetry_error < 1e-12
            # A whole normal segment at centerline radius <=3rho/4 is inside the ball.
            inner = brentq(lambda s: np.linalg.norm(fine_curve([s])[:, 0]-np.array([0., -.5]))-.7*rho,
                           -.5*rho, 0, xtol=1e-15)
            pp = fine_curve([inner])
            nn = geom(pp, U)[2][:, 0]
            test_points = pp[:, 0, None]+nn[:, None]*np.linspace(-rho/8, rho/8, 31)
            containment = float(np.max(np.linalg.norm(test_points-np.array([0., -.5])[:, None], axis=0))/rho)
            assert containment < .825+1e-12
            geometry_cases.append(dict(rho=rho, flux_label=U, theta=th0, cosine=cosine,
                                       signed_half_angle_tangent=float(tau), gaussian_overlap=float(overlap),
                                       ode_refinement_relative=refinement, actual_field_reflection_error=symmetry_error,
                                       cutoff_section_max_radius_over_rho=containment))

            for Xi in [.05, 4.]:
                rows = []
                qs = [10001, 40001, 160001] if U == .5 else [10002, 40002, 160002]
                for Q in qs:
                    ell = (1-U)*Q
                    assert abs(ell-(np.floor(ell)+.5)) < 1e-10
                    outputs = []
                    for order in [480, 640]:
                        T, wt = roots_legendre(order)
                        T = 9*T
                        wt = 9*wt
                        theta = th0+T/np.sqrt(Q*A)
                        assert np.all((theta > 0) & (theta < np.pi))
                        points = np.array([rho*np.sin(theta), -.5+rho*np.cos(theta)])
                        s = np.einsum('i,in->n', t0, points-p0[:, None])
                        for _ in range(7):
                            central = fine_curve(s)
                            b, t, n, omega, curvature, c2 = geom(central, U)
                            displacement = points-central
                            y = np.sum(displacement*n, axis=0)
                            s += np.sum(displacement*t, axis=0)/(1-curvature*y)
                        central = fine_curve(s)
                        b, t, n, omega, curvature, c2 = geom(central, U)
                        y = np.sum((points-central)*n, axis=0)
                        coordinate_error = float(np.max(np.linalg.norm(points-central-n*y, axis=0))/rho)
                        assert coordinate_error < 5e-11
                        metric = 1-curvature*y
                        assert np.min(metric) > .99
                        width = rho*.05
                        cutoff = np.zeros_like(y)
                        cutoff[abs(y) <= width/2] = 1
                        transition = (abs(y) > width/2) & (abs(y) < width)
                        cutoff[transition] = np.cos(np.pi*(abs(y[transition])/width-.5))**2
                        g = (Q*b/np.pi)**.25*np.exp(-Q*b*y*y/2)
                        P = 1-Q*c2*y**3/6
                        V = omega*y
                        beta = np.arctan2(t[0], t[1])
                        angle = beta-theta
                        scalar = normalization*np.sqrt(rho/metric)*cutoff*g
                        # Divide the receiver trace by its central attenuation exp(-Xi L).
                        source_in = scalar*np.exp(-Xi*s)*(np.cos(angle/2)*P-np.sin(angle/2)*V)
                        source_out = scalar*np.exp(-Xi*s)*(np.sin(angle/2)*P+np.cos(angle/2)*V)
                        receiver_angle = -beta-(np.pi-theta)
                        receiver_in = scalar*np.exp(Xi*s)*(np.cos(receiver_angle/2)*P+np.sin(receiver_angle/2)*V)
                        receiver_out = -scalar*np.exp(Xi*s)*(np.sin(receiver_angle/2)*P-np.cos(receiver_angle/2)*V)
                        assert np.max(abs(receiver_out*np.exp(-Xi*s)-source_in*np.exp(Xi*s))) < 3e-11
                        assert np.max(abs(receiver_in*np.exp(-Xi*s)+source_out*np.exp(Xi*s))) < 3e-11
                        jac = 1/np.sqrt(Q*A)
                        norm = lambda f: float(np.sqrt(np.sum(wt*f*f)*jac))
                        # Independently evaluate and normalize the exact intrinsic angular mode.
                        def logprofile(th):
                            th = np.asarray(th, dtype=np.longdouble)
                            rr = np.longdouble(rho)
                            S = np.sqrt(1+rr*rr-2*rr*np.cos(th))
                            W = Q/2*np.log((S+1-rr)/((1-rr)*(S+1+rr)))
                            return ell*np.log(np.sin(th/2))+(Q-ell)*np.log(np.cos(th/2))-W
                        intrinsic = np.exp(np.asarray(logprofile(theta)-logprofile(th0), float))
                        intrinsic /= norm(intrinsic)
                        gaussian = (Q*A*cosine/np.pi)**.25*np.exp(-cosine*T*T/2)
                        intrinsic_gaussian = (Q*A/np.pi)**.25*np.exp(-T*T/2)
                        out_projection = float(np.sum(wt*intrinsic*receiver_out)*jac)
                        outputs.append(dict(
                            coordinate_relative_error=coordinate_error,
                            source_in_to_limit=norm(source_in-gaussian),
                            source_out_to_limit=norm(source_out-tau*gaussian),
                            receiver_in_to_limit=norm(receiver_in+tau*gaussian),
                            receiver_out_to_limit=norm(receiver_out-gaussian),
                            intrinsic_to_limit=norm(intrinsic-intrinsic_gaussian),
                            source_in_mismatch=norm(source_in-intrinsic),
                            receiver_in_norm=norm(receiver_in),
                            receiver_projection_relative=out_projection,
                            source_net_flux=float(np.sum(wt*(source_in**2-source_out**2))*jac)))
                    refinement = max(abs(outputs[0][key]-outputs[1][key]) for key in outputs[0])
                    assert refinement < 3e-7, (rho, U, Xi, Q, refinement)
                    out = outputs[-1]
                    out.update(q=Q, ell=ell, quadrature_refinement=refinement)
                    rows.append(out)
                for key in ['source_in_to_limit', 'source_out_to_limit', 'receiver_in_to_limit',
                            'receiver_out_to_limit', 'intrinsic_to_limit']:
                    assert rows[-1][key] < rows[-2][key]*.57, (rho, U, Xi, key, rows)
                cases.append(dict(rho=rho, flux_label=U, xi_d=Xi,
                                  limiting_wrong_spin_norm=float(abs(tau)),
                                  limiting_receiver_projection_relative=float(overlap), flux_refinements=rows))
    exceptional = []
    for Q in [3, 17, 101, 1001, 10002]:
        labels = (np.arange(Q)+.5)/Q
        for delta in [.01, .05, .2]:
            count = int(np.sum((labels < delta) | (labels > 1-delta)))
            assert count <= 2*delta*Q+2
            exceptional.append(dict(q=Q, trim=delta, excluded=count, upper_bound=2*delta*Q+2))
    return dict(status='passed_global_tube_spherical_trace_checks', exact_identities=exact,
                endpoint_geometry_cases=geometry_cases, two_mouth_trace_cases=cases,
                exceptional_count_cases=exceptional,
                asymptotic_budget=dict(rho_over_d_max=str(eps), wrong_spin_bound=str(tau_bound),
                                      source_shape_bound=str(shape_bound), finite_q_threshold_certified=False),
                global_flat_trial_constructed=True, ordinary_uniform_fixed_geometry_trace_limit=True,
                conditional_trimmed_large_flux_return=True, full_exterior_boundary_realization_accepted=False,
                full_physical_return_error_bounded=False, curved_throat_matching_derived=False,
                novelty_established=False, lean_certified=False, lean_invocations=0,
                numerical_checks_are_supplementary=True,
                versions=dict(numpy=np.__version__, scipy=scipy.__version__, sympy=sp.__version__))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    result = verify()
    result['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({key: len(result[key]) for key in
                      ['exact_identities', 'endpoint_geometry_cases', 'two_mouth_trace_cases', 'exceptional_count_cases']}))
