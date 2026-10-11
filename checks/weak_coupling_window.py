"""Independent finite checks for the angular kernel and joint scale decision.

Run from any directory with Python, numpy, scipy and sympy. No Lean is run.
The ordinary proofs and physical scope are in notes/weak-coupling-window.md.
"""
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy import sparse
from scipy.sparse.linalg import cg
from scipy.special import sph_harm_y
import sympy as sy

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT/'notes/weak-coupling-window.md'
OUT = ROOT/'evidence/weak-coupling-window.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_sum(q):
    beta = Fraction(1)
    normalization = Fraction(1)
    total = Fraction(0)
    for degree in range(1, q):
        beta *= Fraction(q-degree, q+degree)
        normalization += (2*degree+1)*beta
        total += Fraction(2*degree+1, degree*(degree+1))*beta
    harmonic = sum((Fraction(1, k) for k in range(1, q+1)), Fraction(0))-1
    assert total == harmonic and normalization == q
    return dict(q=q, exact_sum=str(total), normalization=str(normalization))


def angular(q, refinement):
    # Direct integration of normalized LLL wavefunctions, independent of
    # a Clebsch-Gordan or matrix-form-factor implementation.
    ntheta, nphi = refinement*(2*q+8), refinement*(4*q+8)
    u, wu = np.polynomial.legendre.leggauss(ntheta)
    phi = 2*np.pi*np.arange(nphi)/nphi
    theta_grid, phi_grid = np.meshgrid(np.arccos(u), phi, indexing='ij')
    theta, phi = theta_grid.ravel(), phi_grid.ravel()
    weights = np.repeat(wu, nphi)*2*np.pi/nphi
    j = (q-1)/2
    eta = np.array([
        math.sqrt(q*math.comb(q-1,k)/(4*np.pi))
        * np.sin(theta/2)**k * np.cos(theta/2)**(q-1-k)
        * np.exp(1j*(j-k)*phi)
        for k in range(q)]).T
    basis_error = float(np.max(np.abs(eta.conj().T@(weights[:,None]*eta)-np.eye(q))))
    square_sum = np.zeros((q,q), complex)
    conjugate_kernel = np.zeros((q*q,q*q), complex)
    pair_kernel = np.zeros_like(conjugate_kernel)
    trace_error, hermitian_error, degree_error = 0., 0., 0.
    beta = 1.
    degree_one = []
    for degree in range(1,q):
        beta *= (q-degree)/(q+degree)
        real_harmonics = [sph_harm_y(degree,0,theta,phi).real]
        for m in range(1,degree+1):
            harmonic = sph_harm_y(degree,m,theta,phi)
            real_harmonics.extend([np.sqrt(2)*harmonic.real,np.sqrt(2)*harmonic.imag])
        degree_sum = np.zeros((q,q), complex)
        for harmonic in real_harmonics:
            vertex = eta.conj().T@((weights*harmonic)[:,None]*eta)
            hermitian_error = max(hermitian_error,float(np.max(np.abs(vertex-vertex.conj().T))))
            degree_sum += vertex@vertex
            trace_error = max(trace_error,abs(float(np.trace(vertex@vertex).real)-q*beta/(4*np.pi)))
            lam = degree*(degree+1)
            conjugate_kernel += np.kron(vertex,vertex.conj())/lam
            pair_kernel += np.kron(vertex,vertex)/lam
            if degree == 1:
                degree_one.append(vertex)
        degree_error = max(degree_error,float(np.max(np.abs(degree_sum-(2*degree+1)*beta/(4*np.pi)*np.eye(q)))))
        square_sum += degree_sum/(degree*(degree+1))
    exact = (sum(1/k for k in range(1,q+1))-1)/(4*np.pi)
    square_error = float(np.max(np.abs(square_sum-exact*np.eye(q))))
    norm = float(np.max(np.abs(np.linalg.eigvalsh(conjugate_kernel))))
    pair_norm = float(np.max(np.abs(np.linalg.eigvalsh(pair_kernel))))
    bell = np.eye(q).ravel()/np.sqrt(q)
    bell_error = float(np.linalg.norm(conjugate_kernel@bell-exact*bell))
    commutator = float(np.linalg.norm(degree_one[0]@degree_one[1]-degree_one[1]@degree_one[0],ord=2))
    assert max(basis_error,trace_error,hermitian_error,degree_error,square_error,bell_error,abs(norm-exact)) < 3e-12
    assert pair_norm <= exact+3e-12 and commutator > 1e-5
    return dict(q=q,refinement=refinement,quadrature=[ntheta,nphi],basis_error=basis_error,
        trace_error=trace_error,hermitian_error=hermitian_error,degree_sum_error=degree_error,
        square_sum_error=square_error,kernel_norm=norm,exact_norm=exact,
        unconjugated_norm=pair_norm,bell_error=bell_error,commutator_norm=commutator)


def poisson_checks():
    x,y,z = sy.symbols('x y z', real=True)
    r = sy.sqrt(x*x+y*y+z*z)
    lap = sum(sy.diff(y/r,v,2) for v in [x,y,z])
    assert sy.simplify(-lap-2*y/r**3) == 0

    def values(point):
        xp,yp,zp = point
        rp = math.sqrt(xp*xp+yp*yp+(zp-.5)**2)
        rm = math.sqrt(xp*xp+yp*yp+(zp+.5)**2)
        return -.5*(yp/rp-yp/rm), -(yp/rp**3-yp/rm**3)

    rows = []
    for point in [np.array([.3,.4,.7]),np.array([-.7,.2,1.1]),np.array([.8,.9,.3])]:
        exact, source = values(point)
        errors = []
        for h in [.02,.01,.005]:
            lap_fd = 0.
            for dim in range(3):
                shift = np.zeros(3); shift[dim]=h
                lap_fd += (values(point+shift)[0]-2*exact+values(point-shift)[0])/h**2
            errors.append(abs(-lap_fd-source))
        assert errors[-1] < errors[0]/12
        assert -.5 <= exact <= 0 and source < 0
        rows.append(dict(point=point.tolist(),steps=[.02,.01,.005],errors=errors))
    sharpness = []
    for epsilon in [.1,.01,.001,.0001]:
        value = abs(values(np.array([0.,epsilon,.5]))[0])
        assert value < .5
        sharpness.append(dict(distance=epsilon,absolute_potential=value))
    assert sharpness[-1]['absolute_potential']>.4999
    return dict(symbolic_identity=True,finite_difference_cases=rows,sharpness_cases=sharpness)


def box_comparison(n):
    # Domain [-2,2] x [0,2] x [0,2], zero Dirichlet data.
    # These are finite-box discrete comparison checks, not whole-space
    # convergence evidence. q*Omega=1 and gamma=g^2*q/(2*pi^2).
    axes = [np.linspace(-2,2,n+2)[1:-1],np.linspace(0,2,n+2)[1:-1],np.linspace(0,2,n+2)[1:-1]]
    steps = [4/(n+1),2/(n+1),2/(n+1)]
    ident = sparse.eye(n,format='csr')
    diff = [sparse.diags([-np.ones(n-1),2*np.ones(n),-np.ones(n-1)],[-1,0,1],format='csr')/h**2 for h in steps]
    lap = sparse.kron(sparse.kron(diff[0],ident),ident,format='csr')+sparse.kron(sparse.kron(ident,diff[1]),ident,format='csr')+sparse.kron(sparse.kron(ident,ident),diff[2],format='csr')
    x,y,z = [part.ravel() for part in np.meshgrid(*axes,indexing='ij')]
    rp = np.sqrt(x*x+y*y+(z-.5)**2)
    rm = np.sqrt(x*x+y*y+(z+.5)**2)
    bx = .5*x*(rp**-3-rm**-3)
    by = .5*y*(rp**-3-rm**-3)
    bz = .5*((z-.5)/rp**3-(z+.5)/rm**3)
    bmag = np.sqrt(bx*bx+by*by+bz*bz)
    source = -2*by
    previous = None
    rows = []
    for gamma in [0.,.02,1.,100.]:
        operator = lap+sparse.diags(gamma*bmag,format='csr')
        solution,info = cg(operator,source,rtol=1e-11,atol=0.,maxiter=3000)
        residual = float(np.linalg.norm(operator@solution-source)/np.linalg.norm(source))
        assert info == 0 and residual < 2e-10
        assert solution.max() < 1e-11 and solution.min() > -.5
        monotone = None if previous is None else float((solution-previous).min())
        if monotone is not None:
            assert monotone >= -1e-10
        rows.append(dict(n=n,gamma=gamma,relative_linear_residual=residual,
            minimum=float(solution.min()),maximum=float(solution.max()),minimum_change=monotone))
        previous = solution
    return rows


def scale_checks():
    a,b = sy.symbols('a b',real=True)
    re = (2+b)/2
    ell = 3*re-1
    distance = ell-a
    omega = (re-3*distance)/2
    td = ell+a
    # Construct from dimensional quantities, then compare to the note.
    powers = {
        'gaussian_rank':-1+5*a,
        'soft_pairs_log':-a,
        'quadrature_log':-1+2*a,
        'hypothetical_return':a-(1-a)/2,
        'nominal_collar':(re+ell)/2-distance,
        'orbital_change':omega+td,
        'electric_duhamel':1+omega+td,
        'angular_proxy_log':1-b+a,
        'singlet_mass_squared':1-b+a,
        'ell_over_re':ell-re,
        'distance_over_re':distance-re,
    }
    expected = [-1+5*a,-a,-1+2*a,-sy.Rational(1,2)+3*a/2,
        -(1+b)/2+a,-(1+b)/2+5*a/2,(1-b)/2+5*a/2,
        1-b+a,1-b+a,1+b,1+b-a]
    for actual,want in zip(powers.values(),expected):
        assert sy.simplify(actual-want) == 0
    cases = []
    for av,bv in [(sy.Rational(1,16),sy.Integer(2)),(sy.Rational(1,8),sy.Integer(2)),(sy.Rational(3,16),sy.Integer(3))]:
        assert 0 < av < sy.Rational(1,5) and bv > 1+5*av
        values = {key:sy.simplify(value.subs({a:av,b:bv})) for key,value in powers.items()}
        for key,value in values.items():
            if key not in ['ell_over_re','distance_over_re']:
                assert value < 0
        cases.append(dict(a=str(av),b=str(bv),powers={k:str(v) for k,v in values.items()}))
    assert cases[0]['powers']['electric_duhamel']=='-11/32'
    assert cases[0]['powers']['distance_over_re']=='47/16'
    return dict(symbolic_powers={key:str(sy.simplify(value)) for key,value in powers.items()},cases=cases)


def mass_comparison():
    # Flat constant-field formula only; C=1, C1=1.82, N=1.
    rows=[]
    for q in [10,100,1000,10000]:
        a,b = 1/16,2.
        alpha = q**(-b)/(4*np.pi)
        direct = .5*math.log(q/2)+math.log(alpha)/3-np.pi/(alpha*math.log(1.82/alpha))+a*math.log(q)
        grouped = (.5+a-b/3)*math.log(q)-4*np.pi**2*q**b/(b*math.log(q)+math.log(4*np.pi*1.82))-.5*math.log(2)-math.log(4*np.pi)/3
        relative_error=abs(direct-grouped)/max(1.,abs(direct))
        assert relative_error < 1e-14 and direct < 0
        rows.append(dict(q=q,log_mass_radius_over_delta=float(direct),relative_formula_error=float(relative_error),flavors=1))
    return rows


def main():
    started=time.perf_counter()
    report=dict(status='passed_weak_coupling_window_checks',
        script_sha256=sha(Path(__file__)),note_sha256=sha(NOTE),
        dependency_hashes={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in [
            ROOT/'notes/finite-flux-response-audit.md',ROOT/'notes/rotating-electric-screening.md',
            ROOT/'notes/singlet-stress-continuum.md']},
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,sympy=sy.__version__),
        exact_sum_cases=[exact_sum(q) for q in range(1,49)],
        angular_cases=[angular(q,refinement) for q in [2,3,5,8] for refinement in [1,2]],
        poisson=poisson_checks(),
        box_comparison_cases=[row for n in [10,16,24] for row in box_comparison(n)],
        scale_checks=scale_checks(),flat_mass_comparison=mass_comparison(),
        angular_norm_proved=True,local_screening_potential_bound_proved=True,
        joint_listed_scale_window=True,interacting_stress_error_bound_proved=False,
        physical_state_constructed=False,physical_damping_proved=False,
        curved_mass_gap_bound_proved=False,novelty_established=False,lean_invocations=0)
    report['elapsed_seconds']=time.perf_counter()-started
    OUT.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=report['status'],exact_sum_cases=len(report['exact_sum_cases']),
        angular_cases=len(report['angular_cases']),box_cases=len(report['box_comparison_cases']),
        max_angular_error=max(abs(row['kernel_norm']-row['exact_norm']) for row in report['angular_cases']),
        max_box_residual=max(row['relative_linear_residual'] for row in report['box_comparison_cases']),
        seconds=report['elapsed_seconds'])))


if __name__ == '__main__':
    main()
