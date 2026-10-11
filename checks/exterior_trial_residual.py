"""Dipole interior trial residual and weighted-current checks; no full return bound."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np
import scipy
from scipy.special import roots_hermite
import sympy as sp

SX=np.array([[0,1],[1,0]],complex)
SY=np.array([[0,-1j],[1j,0]],complex)
SZ=np.diag([1.,-1.])


def verify():
    r,z,u=sp.symbols('r z u',real=True)
    rm=sp.sqrt(r*r+(z+sp.Rational(1,2))**2)
    rp=sp.sqrt(r*r+(z-sp.Rational(1,2))**2)
    flux=((z+sp.Rational(1,2))/rm-(z-sp.Rational(1,2))/rp)/2
    mass=(flux-u)/r
    B=sp.Matrix([r/rm**3-r/rp**3,(z+sp.Rational(1,2))/rm**3-(z-sp.Rational(1,2))/rp**3])/2
    bf=sp.lambdify((r,z),B,'numpy',cse=True)
    jf=sp.lambdify((r,z),B.jacobian([r,z]),'numpy',cse=True)
    hf=sp.lambdify((r,z,u),sp.hessian(mass,[r,z]),'numpy',cse=True)
    mf=sp.lambdify((r,z,u),mass,'numpy',cse=True)
    q,n,b,bs,c=sp.symbols('q n b bs c',positive=True)
    aa=bs/(4*b)
    exact={
        'width_source_cancellation':sp.simplify(aa*(1-2*q*b*n*n)-(bs/(4*b)-q*bs*n*n/2))==0,
        'quadratic_mass_source_cancellation':sp.diff(-q*c*n**3/6,n)+q*c*n*n/2==0,
        'feedback_column_small':sp.simplify((1-sp.sqrt(1-sp.Symbol('tau')**2))**2+sp.Symbol('tau')**2-2*(1-sp.sqrt(1-sp.Symbol('tau')**2)))==0,
    }
    assert all(exact.values())
    moment_powers=[(0,1),(0,2),(1,3),(1,4),(2,5),(2,6)]
    assert all(k>=2*p+1 for p,k in moment_powers)
    exp_upper=Fraction(65,24)+Fraction(1,100)
    budget=Fraction(121,16)*Fraction(103,1000000)+Fraction(14641,256)*Fraction(1,1000000)
    assert exp_upper<Fraction(11,4) and budget<Fraction(1,1000)
    assert Fraction(5,100000)**2+20*Fraction(2,100000)**2<Fraction(103,1000000)**2

    def basic(U,E):
        A=1-E*E;root=np.sqrt(A*A+4*U*U*E*E);mu=(A+root)/(2*U)
        mup=(-2*E+(-4*E*A+8*U*U*E)/(2*root))/(2*U)
        R=.5*np.sqrt((mu*mu-1)*(1-E*E));Z=.5*mu*E
        Rp=.25*(2*mu*mup*(1-E*E)-2*E*(mu*mu-1))/np.sqrt((mu*mu-1)*(1-E*E))
        Zp=.5*(mu+E*mup);speed=np.hypot(Rp,Zp)
        field=np.asarray(bf(R,Z),float).ravel();jac=np.asarray(jf(R,Z),float);bb=np.linalg.norm(field)
        tangent=field/bb;normal=np.array([-tangent[1],tangent[0]])
        assert np.linalg.norm(np.array([Rp,Zp])/speed-tangent)<2e-12
        bprime=float(tangent@jac@tangent)
        curvature=float(normal@jac@tangent/bb)
        cc=float(normal@np.asarray(hf(R,Z,U),float)@normal)
        assert abs(mf(R,Z,U))<2e-13
        return dict(r=R,z=Z,b=bb,bs=bprime,c=cc,a=bprime/(4*bb),k=curvature,speed=speed,normal=normal)

    def geometry(U,E,step):
        g=basic(U,E)
        values=[basic(U,E+i*step) for i in [-2,-1,1,2]]
        for name in ['a','c','k']:
            g[name+'s']=(values[0][name]-8*values[1][name]+8*values[2][name]-values[3][name])/(12*step*g['speed'])
        distances=[g['r'],np.hypot(g['r'],g['z']+.5),np.hypot(g['r'],g['z']-.5),1/max(abs(g['k']),1e-12)]
        g['width']=.18*min(distances)
        return g

    y,weights=roots_hermite(96);weights=weights/np.sqrt(np.pi)
    cases=[]
    for U in [.25,.5,.75]:
        for E in [-.6,-.3,0.,.3,.6]:
            geo=geometry(U,E,1e-4);fine=geometry(U,E,5e-5)
            derivative_error=max(abs(geo[key]-fine[key])/max(1,abs(fine[key])) for key in ['as','cs','ks'])
            assert derivative_error<2e-8
            g=fine;bb=g['b'];w=g['width'];a=g['a'];cc=g['c'];k=g['k']
            for Xi in [.05,4.]:
                rows=[]
                for Q in ([10001,40001,160001] if U==.5 else [10002,40002,160002]):
                    assert abs((1-U)*Q-(np.floor((1-U)*Q)+.5))<1e-10
                    nn=y/np.sqrt(Q*bb);inside=abs(nn)<w;N=nn[inside];W=weights[inside]
                    metric=1-k*N;hs=-g['ks']*N
                    assert np.min(metric)>.8
                    mm=mf(g['r']+N*g['normal'][0],g['z']+N*g['normal'][1],U)
                    gs=a*(1-2*Q*bb*N*N)
                    P=1-Q*cc*N**3/6;Z=a*N;vec=np.array([P,Z])
                    dn=np.array([-Q*cc*N*N/2-Q*bb*N*P,a-Q*bb*N*Z])
                    ds=np.array([-Q*g['cs']*N**3/6+P*(gs-Xi),g['as']*N+Z*(gs-Xi)])
                    actual=1j*SX@dn-1j*SZ@ds/metric+1j*SZ@vec*hs/(2*metric**2)+SY@vec*(Q*mm)-1j*Xi*vec
                    delta=mm-bb*N;rem=delta-cc*N*N/2
                    upper=1j*(gs*(-k*N+Q*cc*N**3/6)/metric+Q*g['cs']*N**3/(6*metric)+Xi*(1/metric-1)*P+hs*P/(2*metric**2)-Q*a*N*delta)
                    lower=1j*(Q*rem-Q*Q*cc*N**3*delta/6+g['as']*N/metric+Z*gs/metric-Xi*Z*(1/metric+1)-hs*Z/(2*metric**2))
                    algebra_error=np.max(abs(actual-np.array([upper,lower])))
                    assert algebra_error<3e-9
                    cut=np.ones_like(N);cp=np.zeros_like(N);transition=abs(N)>w/2
                    angle=np.pi*(abs(N[transition])/w-.5)
                    cut[transition]=np.cos(angle)**2
                    cp[transition]=-np.pi/w*np.sin(2*angle)*np.sign(N[transition])
                    Rcut=cut*actual+1j*SX@vec*cp
                    norm=float(np.sqrt(np.sum(W*np.sum(abs(Rcut)**2,axis=0))))
                    # Leading Gaussian alone in the same exact differential operator.
                    v0=np.array([np.ones_like(N),np.zeros_like(N)])
                    dn0=np.array([-Q*bb*N,np.zeros_like(N)])
                    ds0=np.array([gs-Xi,np.zeros_like(N)])
                    R0=1j*SX@dn0-1j*SZ@ds0/metric+1j*SZ@v0*hs/(2*metric**2)+SY@v0*(Q*mm)-1j*Xi*v0
                    R0=cut*R0+1j*SX@v0*cp
                    leading=float(np.sqrt(np.sum(W*np.sum(abs(R0)**2,axis=0))))
                    rows.append(dict(q=Q,ell=(1-U)*Q,normal_residual_norm=norm,sqrt_q_times_residual=float(np.sqrt(Q)*norm),leading_gaussian_residual_norm=leading,expanded_operator_residual=float(algebra_error)))
                ratio=rows[-1]['sqrt_q_times_residual']/rows[-2]['sqrt_q_times_residual']
                # This is an upper-order bound; cancellation can give faster decay.
                assert ratio<1.03,(U,E,Xi,ratio)
                assert rows[-1]['normal_residual_norm']<rows[-1]['leading_gaussian_residual_norm']
                cases.append(dict(flux_label=U,prolate_eta=E,xi_d=Xi,tube_half_width=w,b=bb,curvature=k,geometry_derivative_refinement_error=float(derivative_error),flux_refinements=rows,last_scaled_residual_ratio=float(ratio)))
    print('Dipole normal-section residual checks passed.',flush=True)

    feedback=[]
    for tau in [.03,.01,.001]:
        r=np.sqrt(1-tau*tau);S=np.array([[r,tau],[tau,-r]])
        assert np.linalg.norm(S.T@S-np.eye(2))<1e-14
        for Xi in [.05,1.,4.]:
            # A causal passive load with a small positive round-trip delay.
            delay=tau**4;load=-np.exp(-Xi*delay)
            loaded=r+tau*load*tau/(1+r*load)
            column_error=np.linalg.norm(S[:,0]-np.array([1.,0.]))
            assert abs(loaded-1)>1.9
            feedback.append(dict(coupling=tau,xi=Xi,load_delay=delay,column_error=float(column_error),loaded_response=float(loaded),error_to_uncoupled_response=float(abs(loaded-1))))
    weighted=[];rng=np.random.default_rng(260710)
    for Xi in [.05,1.,4.]:
        for alpha in [.25,.5,.75]:
            for _ in range(4):
                direction=rng.normal(size=3);direction=alpha*Xi*direction/np.linalg.norm(direction)
                M=Xi*np.eye(2)-direction[0]*SX-direction[1]*SY-direction[2]*SZ
                least=np.linalg.eigvalsh(M)[0]
                assert least>=(1-alpha)*Xi-2e-15
                weighted.append(dict(xi=Xi,weight_fraction=alpha,least_eigenvalue=float(least),coercivity_bound=float((1-alpha)*Xi)))
    return dict(status='passed_dipole_interior_trial_and_weighted_current_checks',exact_identities=exact,moment_power_conditions=[dict(q_power=p,normal_power=k,at_least_inverse_sqrt_q=True) for p,k in moment_powers],conditional_return_budget=dict(weight_fraction=.5,boundary_input_error=5e-5,dimensionless_weighted_residual=2e-5,receiver_trace_error=1e-6,proved_upper_rational=str(budget),target=.001,physical_hypotheses_verified=False),dipole_section_cases=cases,passive_feedback_cases=feedback,weighted_current_cases=weighted,interior_trial_residual_order='q^(-1/2) upper bound on a fixed regular trimmed tube, plus exponential cutoff error',global_matched_trial_constructed=False,full_exterior_return_error_bounded=False,global_boundary_proof_review_complete=False,novelty_established=False,lean_certified=False,lean_invocations=0,numerical_checks_are_supplementary=True,versions=dict(numpy=np.__version__,scipy=scipy.__version__,sympy=sp.__version__))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);args=p.parse_args()
    result=verify();result['script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:len(result[k]) for k in ['exact_identities','dipole_section_cases','passive_feedback_cases','weighted_current_cases']}))
