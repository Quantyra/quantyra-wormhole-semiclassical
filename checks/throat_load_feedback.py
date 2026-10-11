"""Independent finite-model checks for the ordinary throat-load estimate.

No Lean, physical MMP solve, or all-mode numerical truncation claim.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad, solve_bvp
from scipy.special import eval_jacobi, roots_jacobi
import sympy as sp


def constant_load(m, xi, length):
    k = np.hypot(m, xi)
    v = np.tanh(k * length)
    r = m * v / (k + xi * v)
    # Stable also when the transmission is below floating-point range.
    sech = 2 * np.exp(-k * length) / (1 + np.exp(-2 * k * length))
    t = k * sech / (k + xi * v)
    return np.array([[-r, t], [t, r]])


def bvp_load(lam, xi, profile, tol):
    """Solve two independent incoming columns; integrate the actual density."""
    def rhs(x, y):
        m = lam * profile(x)
        return np.array([-xi * y[0] + m * y[1], m * y[0] + xi * y[1]])

    outputs, solutions = [], []
    for incoming in np.eye(2):
        def bc(left, right):
            return np.array([left[0] - incoming[0], right[1] - incoming[1]])
        grid = np.linspace(0, 1, 81)
        guess = np.array([incoming[0] * np.exp(-xi * grid),
                          incoming[1] * np.exp(-xi * (1 - grid))])
        sol = solve_bvp(rhs, bc, grid, guess, tol=tol, max_nodes=12000)
        assert sol.success, sol.message
        solutions.append(sol)
        outputs.append([sol.y[1, 0], sol.y[0, -1]])
    K = np.array(outputs).T
    gram = np.array([[quad(lambda x: np.dot(a.sol(x), b.sol(x)), 0, 1,
                           epsabs=2e-11, epsrel=2e-11)[0]
                      for b in solutions] for a in solutions])
    residual = np.linalg.norm(np.eye(2) - K.T @ K - 2 * xi * gram, 2)
    assert residual < 2e-8
    return K, residual, max(len(sol.x) for sol in solutions)


def verify():
    xi, m = sp.symbols('xi m', positive=True)
    A = sp.Matrix([[-xi, m], [m, xi]])
    J = sp.diag(1, -1)
    assert A.T * J + J * A == -2 * xi * sp.eye(2)
    assert A * A == (xi**2 + m**2) * sp.eye(2)
    r, t = sp.symbols('r t', real=True)
    K = sp.Matrix([[-r, t], [t, r]])
    assert K.T * K == (r**2 + t**2) * sp.eye(2)
    k, C, S = sp.symbols('k C S', positive=True)
    # Exact defect identity, using k^2=m^2+xi^2 and C^2-S^2=1.
    expr = sp.expand((k*C+xi*S)**2-m**2*S**2-k**2-2*xi*S*(k*C+xi*S))
    assert sp.expand(expr.subs(m**2, k**2-xi**2)).factor() == k**2*(C**2-S**2-1)
    exact = ['radial_current', 'transfer_generator_square', 'constant_load_singular_values',
             'constant_load_defect']

    profiles = {
        'constant': (lambda x: np.ones_like(x), 1.),
        'varying_optical_mass': (lambda x: .2 + .8*np.cos(2*np.pi*x)**2, 1.),
        'asymmetric_optical_mass': (lambda x: .6 + .3*x + .1*np.sin(2*np.pi*x), 1.),
    }
    slab = []
    for name, (profile, ceiling) in profiles.items():
        for lam in [.7, 3., 12.]:
            for s in [.05, .7, 3.]:
                coarse, _, _ = bvp_load(lam, s, profile, 2e-7)
                fine, current, nodes = bvp_load(lam, s, profile, 2e-9)
                difference = np.linalg.norm(coarse-fine, 2)
                assert difference < 2e-8
                c = s + lam * ceiling
                # Two disjoint collars, each of optical length 1/4.
                delta = s/c * (-np.expm1(-c/2))
                defect = np.linalg.eigvalsh(np.eye(2)-fine.T@fine).min()
                adjoint_defect = np.linalg.eigvalsh(np.eye(2)-fine@fine.T).min()
                assert min(defect, adjoint_defect) >= delta - 2e-9
                exact_error = None
                if name == 'constant':
                    exact_error = np.linalg.norm(fine-constant_load(lam, s, 1), 2)
                    assert exact_error < 2e-9
                slab.append(dict(profile=name, lam=lam, xi=s, collar_bound=delta,
                                 smallest_defect=float(defect), current_residual=float(current),
                                 refinement=float(difference), constant_formula_error=exact_error,
                                 fine_mesh_nodes=nodes))

    rng = np.random.default_rng(20261010)
    completion, loaded = [], []
    for n in [1, 3, 7]:
        for repetition in range(12):
            left, _ = np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
            right, _ = np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
            singular = np.linspace(.15, .999, n)
            K = left @ np.diag(singular) @ right.conj().T
            b = rng.normal(size=n)+1j*rng.normal(size=n)
            y = rng.normal(size=n)+1j*rng.normal(size=n)
            D = np.eye(n)-K.conj().T@K
            E = np.eye(n)-K@K.conj().T
            z = np.linalg.solve(D, K.conj().T@b)
            lhs = np.vdot(y,y).real - np.vdot(K@y+b,K@y+b).real
            rhs = np.vdot(y-z,D@(y-z)).real - np.vdot(b,np.linalg.solve(E,b)).real
            err = abs(lhs-rhs)
            assert err < 2e-8
            # At the maximizing response this inverse-defect cost is attained.
            sharp = np.vdot(K@z+b,K@z+b).real - np.vdot(z,z).real
            cost = np.vdot(b,np.linalg.solve(E,b)).real
            assert abs(sharp-cost) < 2e-7 * max(1,cost)
            completion.append(dict(size=n, identity_residual=float(err), cost=float(cost)))

            # Independent direct closure of a passive exterior network.
            # Trial solves the exterior exactly but uses no returning high input.
            ng = 2
            raw = rng.normal(size=(ng+n,ng+n))+1j*rng.normal(size=(ng+n,ng+n))
            ext = raw / np.linalg.norm(raw,2) * .995
            gin = rng.normal(size=ng)+1j*rng.normal(size=ng)
            out_trial = ext @ np.r_[gin,np.zeros(n)]
            high = np.linalg.solve(np.eye(n)-ext[ng:,ng:]@K, out_trial[ng:])
            out_exact = ext @ np.r_[gin,K@high]
            residual = K@out_trial[ng:]
            budget = np.vdot(residual,np.linalg.solve(E,residual)).real
            error = np.linalg.norm(out_exact[:ng]-out_trial[:ng])**2
            assert error <= budget+2e-10
            loaded.append(dict(high_dimension=n, ground_error_squared=float(error),
                               inverse_defect_budget=float(budget)))

    hierarchy = []
    for mass in [1.,10.,100.,1000.,10000.]:
        s = .05
        K = constant_load(mass,s,1.)
        gap = 2*s*np.tanh(np.hypot(mass,s))/(np.hypot(mass,s)+s*np.tanh(np.hypot(mass,s)))
        assert abs(1-np.linalg.norm(K,2)**2-gap) < 2e-14
        hierarchy.append(dict(mass=mass, reflection=float(abs(K[0,0])),
                              transmission=float(K[1,0]), defect=float(gap),
                              defect_times_mass_over_2xi=float(gap*mass/(2*s))))
    assert hierarchy[-1]['reflection'] > .99999

    # The front is the throat length, not an external distance D>L.
    fronts=[]
    for s in [10.,20.,40.]:
        L,D,mass=.25,1.,3.
        transmission=constant_load(mass,s,L)[1,0]
        ratio=transmission*np.exp(s*D)
        assert ratio>1
        fronts.append(dict(xi=s,length=L,external_distance=D,
                           transmission_over_exterior_front=float(ratio)))

    # L2 convergence does not imply convergence in the required trace form.
    tails=[]
    for lam in [10.,100.,1000.,10000.]:
        s=.05;b=1.;ell=.25
        delta=s/(s+lam*b)*(-np.expm1(-2*(s+lam*b)*ell))
        error_norm_squared=lam**(-.5)
        tails.append(dict(lam=lam,L2_error_squared=error_norm_squared,
                          weighted_error_squared=error_norm_squared/delta))
    assert all(tails[i+1]['weighted_error_squared']>tails[i]['weighted_error_squared']
               for i in range(len(tails)-1))
    determinants=[]
    for mass in [.5,3.,8.]:
        for s in [.1,1.]:
            for length in [.2,.7,1.3]:
                raw=rng.normal(size=(4,4))+1j*rng.normal(size=(4,4))
                ext=.8*raw/np.linalg.norm(raw,2)
                G,B,C,H=ext[:2,:2],ext[:2,2:],ext[2:,:2],ext[2:,2:]
                def determinants_at(L):
                    K=constant_load(mass,s,L)
                    Q=-np.exp(-s*L)*np.eye(2)
                    throat=np.zeros((4,4));throat[:2,:2]=Q;throat[2:,2:]=K
                    R=G+B@K@np.linalg.solve(np.eye(2)-H@K,C)
                    full=np.linalg.det(np.eye(4)-ext@throat)
                    factors=np.linalg.det(np.eye(2)-H@K)*np.linalg.det(np.eye(2)-R@Q)
                    return full,factors,R,K,Q
                full,factors,R,K,Q=determinants_at(length)
                assert abs(full-factors)<2e-14
                kap=np.hypot(mass,s);v=np.tanh(kap*length)
                rp=mass*kap**2*(1-v*v)/(kap+s*v)**2
                tp=-kap*(kap*v+s)/(kap+s*v)*K[0,1]
                Kp=np.array([[-rp,tp],[tp,rp]])
                Qp=-s*Q
                Rp=B@np.linalg.solve(np.eye(2)-K@H,Kp)@np.linalg.solve(np.eye(2)-H@K,C)
                high=-np.trace(np.linalg.solve(np.eye(2)-H@K,H@Kp))
                low=-np.trace(np.linalg.solve(np.eye(2)-R@Q,Rp@Q+R@Qp))
                step=1e-5
                derivative=(np.log(determinants_at(length+step)[0])-np.log(determinants_at(length-step)[0]))/(2*step)
                assert abs(derivative-high-low)<2e-8
                frozen=-np.trace(np.linalg.solve(np.eye(2)-R@Q,R@Qp))
                determinants.append(dict(mass=mass,xi=s,length=length,
                    factorization_residual=float(abs(full-factors)),
                    derivative_residual=float(abs(derivative-high-low)),
                    missed_real_derivative_if_load_frozen=float((high+low-frozen).real)))
    basis=[]
    for q in [3,16,128]:
        for azimuth in [0,q//2,q-1]:
            for eps in [.01,.04,.1]:
                weighted=[]
                for order in [96,160]:
                    z,w=roots_jacobi(order,q-azimuth-1,azimuth)
                    x=(z+1)/2;w=w/w.sum();c=1-2*x
                    distance=np.sqrt(1+eps*eps-2*eps*c)
                    denominator=distance*(distance+1-eps*c)
                    W=q/2*np.log((distance+1-eps)/((1-eps)*(distance+1+eps)))
                    g=np.exp(-W);g/=np.sqrt(np.sum(w*g*g))
                    ground=np.sum(w*g)
                    high_norm_sq=max(0.,1-ground**2)
                    actual_derivative_sq=np.sum(w*g*g*(q*eps**2*np.sqrt(x*(1-x))/denominator)**2)
                    bound=q*eps**2/(4*(1-eps)**2)
                    gap=np.sqrt(q+1)
                    assert actual_derivative_sq<=bound**2+1e-15
                    assert high_norm_sq<=actual_derivative_sq/gap**2+2e-14
                    coefficients=[]
                    for j in range(1,29):
                        p=eval_jacobi(j,azimuth,q-azimuth-1,1-2*x)
                        p/=np.sqrt(np.sum(w*p*p))
                        coefficients.append(np.sum(w*g*p))
                    lam=np.sqrt(np.arange(1,29)*(np.arange(1,29)+q))
                    half_derivative_sq=np.sum(lam*np.array(coefficients)**2)
                    assert half_derivative_sq<=actual_derivative_sq/gap+2e-12
                    xi=.05;beta=1/eps;ell=eps/4
                    delta=xi/(xi+lam*beta)*(-np.expm1(-2*(xi+lam*beta)*ell))
                    cost=np.sum(np.array(coefficients)**2/delta)
                    d0=-np.expm1(-2*(xi+gap*beta)*ell)
                    upper=(bound**2/(q+1)+beta*bound**2/(xi*gap))/d0
                    assert cost<=upper+2e-10
                    weighted.append(float(cost))
                assert abs(weighted[0]-weighted[1])<2e-8*max(1.,weighted[1])
                basis.append(dict(q=q,azimuth=azimuth,epsilon=eps,
                    high_projection_norm_squared=float(high_norm_sq),
                    partial_28_mode_cost=weighted[1],all_mode_upper=float(upper),
                    quadrature_difference=abs(weighted[0]-weighted[1])))
    return dict(status='passed',exact_identities=exact,warped_load_cases=slab,
                nonnormal_completion_cases=completion,loaded_network_cases=loaded,
                large_mass_cases=hierarchy,front_counterexamples=fronts,trace_tail_cases=tails,
                determinant_cases=determinants,
                mouth_basis_cases=basis,
                absolute_window_target=float(.001*np.exp(-4)),
                ordinary_all_mode_load_defect_bound=True,
                loaded_residual_estimate_conditional_on_exterior_solution=True,
                physical_dipole_matching_proved=False,full_energy_derived=False,
                novelty_established=False,lean_invocations=0,
                numerical_checks_are_supplementary=True,
                versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__))


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    result=verify()
    result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:len(result[k]) for k in ['exact_identities','warped_load_cases',
          'nonnormal_completion_cases','loaded_network_cases','large_mass_cases',
          'front_counterexamples','trace_tail_cases','determinant_cases','mouth_basis_cases']}))
