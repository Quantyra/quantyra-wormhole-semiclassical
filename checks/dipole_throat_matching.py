"""Classical l=1 mouth gauge and driven AdS2 x S2 reconstruction checks.

This does not solve the moving binary or its quantum state.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import quad, solve_bvp, solve_ivp
from scipy.special import roots_laguerre
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]


def text_hash(path):
    return hashlib.sha256(path.read_text(encoding='utf-8').encode()).hexdigest()


def tensor_checks():
    r,c,sigma=sp.symbols('rho c sigma',real=True)
    a=sp.Function('a')(r)
    f,sphere=1+r*r,1-c*c
    def derivative(expr,i,linear=False):
        if i==0:
            return -sp.I*sigma*expr if linear else sp.S.Zero
        if i==1: return sp.diff(expr,r)
        if i==2: return sp.diff(expr,c)
        return sp.S.Zero
    metric=sp.diag(-f,1/f,1/sphere,sphere)
    inv=sp.diag(-1/f,f,sphere,1/sphere)
    h=sp.diag(-4*a*c*f,4*a*c/f,0,0)
    dinv=-inv*h*inv
    gamma,dgamma={},{}
    for i in range(4):
        for j in range(4):
            for k in range(4):
                base=derivative(metric[i,k],j)+derivative(metric[i,j],k)-derivative(metric[j,k],i)
                linear=derivative(h[i,k],j,True)+derivative(h[i,j],k,True)-derivative(h[j,k],i,True)
                gamma[i,j,k]=sp.factor(inv[i,i]*base/2)
                dgamma[i,j,k]=sp.factor((inv[i,i]*linear+dinv[i,i]*base)/2)
    ric,dric=sp.zeros(4),sp.zeros(4)
    for i in range(4):
        for j in range(i,4):
            base=linear=sp.S.Zero
            for k in range(4):
                base+=derivative(gamma[k,i,j],k)-derivative(gamma[k,i,k],j)
                linear+=derivative(dgamma[k,i,j],k,True)-derivative(dgamma[k,i,k],j,True)
                for n in range(4):
                    base+=gamma[k,i,j]*gamma[n,k,n]-gamma[n,i,k]*gamma[k,j,n]
                    linear+=(dgamma[k,i,j]*gamma[n,k,n]+gamma[k,i,j]*dgamma[n,k,n]
                             -dgamma[n,i,k]*gamma[k,j,n]-gamma[n,i,k]*dgamma[k,j,n])
            ric[i,j]=ric[j,i]=sp.factor(base)
            dric[i,j]=dric[j,i]=sp.factor(linear)
    scalar=sp.factor(sp.trace(inv*ric))
    dscalar=sp.factor(sp.trace(dinv*ric+inv*dric))
    einstein=ric-metric*scalar/2
    deinstein=dric-h*scalar/2-metric*dscalar/2
    F,dF=sp.zeros(4),sp.zeros(4)
    F[2,3],F[3,2]=-1,1
    for i in range(3):
        dF[i,3]=derivative(a*sphere,i,True)
        dF[3,i]=-dF[i,3]
    upper=inv*F*inv
    dupper=dinv*F*inv+inv*dF*inv+inv*F*dinv
    square=sum(F[i,j]*upper[i,j] for i in range(4) for j in range(4))
    dsquare=sum(dF[i,j]*upper[i,j]+F[i,j]*dupper[i,j] for i in range(4) for j in range(4))
    stress=2*(F*inv*F.T-metric*square/4)
    dstress=2*(dF*inv*F.T+F*dinv*F.T+F*inv*dF.T-h*square/4-metric*dsquare/4)
    a2=((6-sigma*sigma/f)*a-2*r*sp.diff(a,r))/f
    checks={}
    def zero(name,expr):
        out=sp.factor(expr.subs(sp.diff(a,r,2),a2))
        assert out==0,(name,str(out))
        checks[name]='exact_zero'
    for i in range(4):
        for j in range(i,4):
            zero(f'background_Einstein_Maxwell_{i}_{j}',einstein[i,j]-stress[i,j])
            zero(f'linear_Einstein_Maxwell_{i}_{j}',deinstein[i,j]-dstress[i,j])
    dvolume=sp.trace(inv*h)/2
    for j in range(4):
        zero(f'background_Maxwell_{j}',sum(derivative(upper[i,j],i) for i in range(4)))
        zero(f'linear_Maxwell_{j}',sum(derivative(dvolume*upper[i,j]+dupper[i,j],i,True) for i in range(4)))
    zero('linear_scalar_curvature',dscalar)
    print('four_dimensional_reconstruction_passed',flush=True)
    return dict(checks=checks,count=len(checks))


def algebra_checks():
    x,h,c,sigma=sp.symbols('x h c sigma',real=True)
    r=1+h
    U=sp.Function('U')(h)
    Up=r-(1+2*h)*U/(r*h*h)
    # Lie derivative of the static RN background, units P=1 and BP factored out.
    f=h*h/r**2
    base=sp.diag(-f,1/f,r*r/(1-c*c),r*r*(1-c*c))
    rel=[2*h*c,2*h*c,2*r*c,2*r*c]
    perturbation=sp.diag(*[base[i,i]*rel[i] for i in range(4)])
    xi=[0,U*c,(U/r-r)*(1-c*c),0]
    coords=[sp.Symbol('t'),h,c,sp.Symbol('phi')]
    lie=sp.Matrix(4,4,lambda i,j:sum(xi[k]*sp.diff(base[i,j],coords[k])
             +base[k,j]*sp.diff(xi[k],coords[i])+base[i,k]*sp.diff(xi[k],coords[j]) for k in range(4)))
    transformed=(perturbation-lie).subs(sp.diff(U,h),Up).applyfunc(sp.factor)
    assert transformed[2,2]==0 and transformed[3,3]==0 and transformed[1,2]==0
    I=sp.symbols('I')
    tt=sp.factor((transformed[0,0]/base[0,0]/c).subs(U,h*h*r*I))
    rr=sp.factor((transformed[1,1]/base[1,1]/c).subs(U,h*h*r*I))
    assert sp.factor(tt-2*h*(1-I))==0
    assert sp.factor(rr-2*((1+3*h)*I-1))==0
    anew=(r*r+1)/2+U/r-r
    assert sp.factor(anew.subs(U,h*h*r*I)-h*h*(sp.Rational(1,2)+I))==0
    # Integral expansion and a bounded third-order remainder give the overlap.
    truncated=1-3*h+12*h*h-60*h**3
    assert sp.limit((h*h*(sp.Rational(1,2)+truncated))/(h*h),h,0)==sp.Rational(3,2)
    assert sp.limit(tt.subs(I,truncated)/(h*h),h,0)==6
    assert sp.limit(rr.subs(I,truncated)/(h*h),h,0)==6
    Ip=(1-(1+3*h)*I)/(h*h)
    Uprime=sp.diff(h*h*r,h)*I+h*h*r*Ip
    assert sp.factor(Uprime-Up.subs(U,h*h*r*I))==0
    # Darboux solutions are verified in x, avoiding any imported master equation.
    seed=sp.exp(sp.I*sigma*x)
    step1=sp.diff(seed,x)+sp.tan(x)*seed
    step2=sp.diff(step1,x)+2*sp.tan(x)*step1
    expected=(1-sigma*sigma+3*sp.tan(x)**2+3*sp.I*sigma*sp.tan(x))*seed
    assert sp.trigsimp(step2-expected)==0
    assert sp.trigsimp(sp.diff(step2,x,2)+(sigma*sigma-6/sp.cos(x)**2)*step2)==0
    w=3*sp.tan(x)
    assert sp.trigsimp(w*w-sp.diff(w,x)+9-6/sp.cos(x)**2)==0
    ground=sp.cos(x)**3
    assert sp.trigsimp(-sp.diff(ground,x,2)+6*ground/sp.cos(x)**2-9*ground)==0
    # Independent predecessor: Burko (1995), equations (40), (41), (87).
    # This checks its extremal throat limit, not the full RN reconstruction.
    radius,P=sp.symbols('radius P',positive=True)
    delta=(radius-P)**2
    varpi=3*P-2*P**2/radius
    W=3*P*delta/(radius*varpi**2)+P/varpi
    V=delta/radius**5*(2*P-2*P**2/radius+6*P*W)
    assert sp.limit(V/(delta/radius**2),radius,P)==6/P**2
    exponent=sp.symbols('exponent')
    assert set(sp.solve(exponent*(exponent-1)+2*exponent-6,exponent))=={-3,2}
    rho=sp.symbols('rho',real=True)
    even=(1-sigma**2+3*rho**2)*sp.cos(sigma*sp.atan(rho))-3*sigma*rho*sp.sin(sigma*sp.atan(rho))
    odd=(1-sigma**2+3*rho**2)*sp.sin(sigma*sp.atan(rho))+3*sigma*rho*sp.cos(sigma*sp.atan(rho))
    assert even.subs(rho,0)==1-sigma**2
    assert sp.factor(sp.diff(odd,rho).subs(rho,0)-sigma*(4-sigma**2))==0
    assert sp.limit(sigma*(4-sigma**2)/(6*sp.sin(sp.pi*sigma/2)),sigma,0)==4/(3*sp.pi)
    q=sp.symbols('q',positive=True)
    lam=1/(16*sp.pi*q**3)
    separation=1/sp.log(q)**2
    core=lam**4/separation**2
    opening=1/(4*sp.pi**2*q**3)
    assert sp.simplify(core-sp.log(q)**4/(65536*sp.pi**4*q**12))==0
    assert sp.simplify(core/opening-sp.log(q)**4/(16384*sp.pi**2*q**9))==0
    return dict(gauge_conditions_checked=3,matching_a_coefficient='3/2',matching_metric_coefficient=6,
                exact_gap_squared=9,first_order_core_scale='O((BP)(P/ell)^2)',
                physical_scale='O(q^-12 (log q)^4)',burko_extremal_mass_squared='6/P^2',
                burko_static_indicial_exponents=[-3,2],central_transfer_and_scale_identities_checked=True)


def gauge_checks():
    def integral(h):
        return quad(lambda s:np.exp(-s)/(1+h*s)**3,0,np.inf,epsabs=2e-14,epsrel=2e-14)[0]
    cases=[]
    for h in [1e-5,1e-4,.001,.003,.01,.03,.1]:
        value=integral(h)
        approx=[]
        for order in [32,64]:
            nodes,weights=roots_laguerre(order)
            approx.append(float(np.dot(weights,(1+h*nodes)**-3)))
        assert abs(value-approx[-1])<1e-12
        assert 0<value<=1 and value>=1-3*h-1e-13
        assert abs(value-(1-3*h))<=12*h*h+1e-13
        assert abs(value-(1-3*h+12*h*h))<=60*h**3+1e-13
        cases.append(dict(h=h,I=value,matching_a_over_eps_h2=.5+value,
                          integral_discrepancy=abs(value-approx[-1]),
                          quadrature_refinement=abs(approx[1]-approx[0])))
    def ode(h,y):
        return [(1-(1+3*h)*y[0])/(h*h)]
    grid=np.geomspace(.001,.1,121)
    sol=solve_ivp(ode,(grid[0],grid[-1]),[integral(grid[0])],method='Radau',
                  t_eval=grid,rtol=1e-11,atol=1e-13)
    assert sol.success,sol.message
    error=max(abs(sol.y[0,i]-integral(h)) for i,h in enumerate(grid))
    assert error<1e-9
    return dict(integral_cases=cases,independent_ode_error=float(error))


def basis(r,sigma):
    angle=np.arctan(r)
    k=1-sigma*sigma+3*r*r
    even=k*np.cos(sigma*angle)-3*sigma*r*np.sin(sigma*angle)
    odd=k*angle*np.sinc(sigma*angle/np.pi)+3*r*np.cos(sigma*angle)
    return even,odd


def driven_checks():
    cases=[]
    core=np.linspace(-1,1,201)
    for sigma in [0.,.25,.5]:
        for minus,plus in [(1.,1.),(-.4,1.)]:
            even,odd=basis(core,sigma)
            target=(plus+minus)*even/(6*np.cos(sigma*np.pi/2))+(plus-minus)*odd/(3*np.pi*np.sinc(sigma/2))
            previous=None
            for end in [3.,10.,30.]:
                ymax=np.arcsinh(end)
                def ode(y,z):
                    return np.array([z[1],(6-sigma*sigma/np.cosh(y)**2)*z[0]-np.tanh(y)*z[1]])
                def bc(left,right):
                    return np.array([left[0]-minus*end*end,right[0]-plus*end*end])
                solved=[]
                for count,tol in [(120,1e-8),(240,1e-10)]:
                    mesh=np.linspace(-ymax,ymax,count)
                    guess=np.array([np.linspace(minus,plus,count)*end*end,np.zeros(count)])
                    sol=solve_bvp(ode,bc,mesh,guess,tol=tol,max_nodes=20000)
                    assert sol.success,sol.message
                    solved.append(sol.sol(np.arcsinh(core))[0])
                er,od=basis(end,sigma)
                exact_finite=(plus+minus)*end*end*even/(2*er)+(plus-minus)*end*end*odd/(2*od)
                solver_error=float(np.max(abs(solved[-1]-exact_finite)))
                finite_error=float(np.max(abs(solved[-1]-target)))
                refinement=float(np.max(abs(solved[-1]-solved[0])))
                assert solver_error<1e-8 and refinement<1e-7
                assert finite_error*end*end<3
                if previous is not None: assert finite_error<previous
                previous=finite_error
                cases.append(dict(sigma=sigma,D_minus=minus,D_plus=plus,overlap_radius=end,
                                  ode_exact_finite_error=solver_error,solver_refinement=refinement,
                                  core_infinite_boundary_error=finite_error,
                                  scaled_overlap_error=finite_error*end*end))
    return cases


def verify():
    return dict(status='passed_classical_dipole_throat_matching',tensor=tensor_checks(),
                algebra=algebra_checks(),gauge=gauge_checks(),driven_cases=driven_checks(),
                static_RN_to_AdS_leading_coefficient_derived=True,
                dynamic_RN_exterior_matching_proved=False,
                boundary_gauges_and_holonomy_fixed=False,quantum_open_background_corrected=False,
                time_dependent_quantum_state_derived=False,physical_transport_proved=False,
                novelty_established=False,lean_invocations=0,
                versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    result=verify()
    result['script_sha256']=text_hash(Path(__file__))
    result['input_sha256']={name:text_hash(ROOT/name) for name in [
        'notes/dipole-throat-matching.md','notes/magnetic-mouth-force.md',
        'notes/casimir-reference-matching.md']}
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=result['status'],tensor_checks=result['tensor']['count'],driven_cases=len(result['driven_cases']))))
