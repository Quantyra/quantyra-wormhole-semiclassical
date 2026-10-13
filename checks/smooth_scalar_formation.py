"""S069 stationary Einstein--complex-scalar backgrounds and endpoint diagnostic.

Units G=c=hbar=mu_scalar=1; potential/action from
arXiv:2504.17775v3 (1). This is not a wormhole evolution or stability proof.
"""
import argparse
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_bvp, solve_ivp, cumulative_trapezoid, simpson
from scipy.optimize import minimize_scalar, brentq
from scipy.special import expit


def equations(sigma):
    def fun(r, y, p):
        m, nu, f, h = y
        w = expit(p[0])
        g = 1 - 2*m/r
        t = w*w*np.exp(-2*nu)
        v = sigma*sigma*f*f*(1-2*f*f)**2
        vx = (1-2*f*f)*(1-6*f*f)
        rho = .5*(sigma*sigma*(t*f*f+g*h*h)+v)
        pr = .5*(sigma*sigma*(t*f*f+g*h*h)-v)
        mp = 4*np.pi*r*r*rho
        np_ = (m+4*np.pi*r**3*pr)/(r*r*g)
        gp = -2*mp/r+2*m/r**2
        hp = -(2/r+np_+gp/(2*g))*h+(vx-t)*f/g
        return np.array([mp, np_, h, hp])
    return fun


def boundary(sigma, fc, r0, rend):
    def bc(a, b, p):
        w = expit(p[0])
        tc = w*w*np.exp(-2*a[1])
        vc = sigma*sigma*fc*fc*(1-2*fc*fc)**2
        rhoc = .5*(sigma*sigma*tc*fc*fc+vc)
        fpp = ((1-2*fc*fc)*(1-6*fc*fc)-tc)*fc/3
        kap = np.sqrt(1-w*w)
        eps = b[0]*(1-2*w*w)/kap
        return np.array([a[0]-4*np.pi*rhoc*r0**3/3,
                         a[2]-fc-.5*fpp*r0*r0,
                         a[3]-fpp*r0,
                         b[1]-.5*np.log(1-2*b[0]/rend),
                         b[3]+(kap+(1+eps)/rend)*b[2]])
    return bc


def guess(r, sigma, central, mass, wall, nuc):
    # Smooth monotone seed; it is not used as evidence unless BVP converges.
    nuout = .5*np.log(1-2*mass/r[-1])
    nu = nuc+(nuout-nuc)*(1-np.exp(-(r/wall)**2))/(1-np.exp(-(r[-1]/wall)**2))
    w = np.sqrt((1-2*(central/sigma)**2)*(1-6*(central/sigma)**2))*np.exp(nuc)
    t = w*w*np.exp(-2*nu)
    f = np.sqrt((2+np.sqrt(1+3*t))/6)*np.sqrt(expit(-2*(r-wall)))
    h = np.gradient(f, r)
    # Bounded compactness for initial iterate.
    m = mass*r**3/(r**3+wall**3)
    return np.array([m, nu, f, h]), [np.log(w/(1-w))]


def solve_model(model, tol=1e-6, r0=1e-4, rend=90., mesh=1201, seed=None):
    sigma, central, mass, radius = (.08,.06,5.65,17.9) if model=='S08A06' else (.06,.044,10.43,31.7)
    r = np.linspace(r0, rend, mesh)
    if seed is None:
        y,p = guess(r,sigma,central,mass,radius-3,-1.2)
    else:
        y = seed.sol(r)
        p = seed.p
        outside = r > seed.x[-1]
        if np.any(outside):
            edge = seed.x[-1]
            mm, _, ff, _ = seed.y[:,-1]
            w = expit(p[0])
            kap = np.sqrt(1-w*w)
            eps = mm*(1-2*w*w)/kap
            y[0,outside] = mm
            y[1,outside] = .5*np.log(1-2*mm/r[outside])
            y[2,outside] = ff*np.exp(-kap*(r[outside]-edge))*(r[outside]/edge)**(-1-eps)
            y[3,outside] = -(kap+(1+eps)/r[outside])*y[2,outside]
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        sol = solve_bvp(equations(sigma),boundary(sigma,central/sigma,r0,rend),
                        r,y,p=p,tol=tol,max_nodes=50000)
    if not sol.success:
        raise RuntimeError(f'{model}: {sol.message}; nodes={sol.x.size}')
    return sol, sigma


def diagnose(sol, sigma):
    r0, rend = sol.x[[0,-1]]
    w = expit(sol.p[0])
    _, nu0, fc, _ = sol.y[:,0]
    pc = .5*sigma*sigma*fc*fc*(w*w*np.exp(-2*nu0)-(1-2*fc*fc)**2)
    def ff(r, F):
        y = sol.sol(r)
        return 1/np.sqrt(1-2*y[0]/r)-F*equations(sigma)(r,y,sol.p)[1]
    integral = solve_ivp(ff,(r0,rend),[r0-4*np.pi*pc*r0**3/3],
                         method='DOP853',rtol=2e-12,atol=1e-13,dense_output=True,max_step=.15)
    assert integral.success
    def dd(r):
        y = sol.sol(r)
        return integral.sol(r)[0]*np.sqrt(1-2*y[0]/r)*equations(sigma)(r,y,sol.p)[1]-1
    # Search every local peak, then refine the largest. No root inferred from
    # samples alone; spatial-resolution and alternative-integral checks follow.
    r = np.linspace(r0,rend,24001)
    y = sol.sol(r)
    m, nu, f, h = y
    A = np.exp(nu)
    g = 1-2*m/r
    dr = dd(r)
    indices = np.flatnonzero((dr[1:-1]>=dr[:-2]) & (dr[1:-1]>=dr[2:]))+1
    peaks = [minimize_scalar(lambda x: -float(dd(x)),bounds=(r[i-1],r[i+1]),method='bounded',options={'xatol':1e-11}) for i in indices]
    peak = max(peaks,key=lambda x: -x.fun)
    r99 = brentq(lambda x: sol.sol(x)[0]-.99*m[-1],r0,rend,xtol=1e-12)
    N = simpson(4*np.pi*r*r*w*sigma*sigma*f*f/(A*np.sqrt(g)),x=r)
    # Independent integration of the original lapse integral, with h halving.
    quad_errors = []
    for n in (4001,8001,16001):
        rq = np.linspace(r0,peak.x,n)
        yq = sol.sol(rq)
        iq = np.exp(nu0)*r0+simpson(np.exp(yq[1])/np.sqrt(1-2*yq[0]/rq),x=rq)
        yp = sol.sol(peak.x)
        dquad = iq*np.sqrt(1-2*yp[0]/peak.x)*equations(sigma)(peak.x,yp,sol.p)[1]/np.exp(yp[1])-1
        quad_errors.append(float(abs(dquad+peak.fun)))
    # Angular Einstein equation, not part of the first-order BVP constraints.
    # Midpoints avoid derivative ambiguities at cubic collocation knots.
    rc = (sol.x[1:]+sol.x[:-1])/2
    yc, dy, ddy = sol.sol(rc),sol.sol(rc,1),sol.sol(rc,2)
    mc,nuc,fc,hc = yc
    gc = 1-2*mc/rc
    gp = -2*dy[0]/rc+2*mc/rc**2
    pt = .5*sigma*sigma*(w*w*np.exp(-2*nuc)*fc*fc-gc*hc*hc-fc*fc*(1-2*fc*fc)**2)
    gt = gc*(ddy[1]+dy[1]**2+dy[1]/rc)+gp*(dy[1]+1/rc)/2
    angular_error = float(np.max(abs(gt-8*np.pi*pt)))
    ode_error = float(np.max(abs(sol.sol(rc,1)-equations(sigma)(rc,yc,sol.p))))
    # Tail envelope at R>=70: |psi|,|psi'| <= 1e-12 exp[-(r-R)/2].
    # This deliberately loose asymptotic envelope is checked at the matching
    # point and against domain/Robin variations, not interval-certified.
    amp_envelope = 1e-12
    assert max(abs(sigma*f[-1]),abs(sigma*h[-1])) < amp_envelope/100
    gR = g[-1]
    AR = A[-1]
    tail_rho_coefficient = .5*((w/AR)**2+1+1.0001)*amp_envelope**2
    tail_mass_bound = 4*np.pi*tail_rho_coefficient*(rend**2+2*rend+2)
    tail_pr_coefficient = .5*((w/AR)**2+1)*amp_envelope**2
    tail_D_upper = (m[-1]+tail_mass_bound)/(rend*np.sqrt(gR)*AR)+4*np.pi*rend**2*tail_pr_coefficient/(np.sqrt(gR)*AR)-1
    # Light rings are a distinct observable. Record them to guard against
    # confusing published ultracompactness with the endpoint criterion.
    lr = r*equations(sigma)(r,y,sol.p)[1]-1
    light_rings = [brentq(lambda z: z*equations(sigma)(z,sol.sol(z),sol.p)[1]-1,r[i],r[i+1])
                   for i in np.flatnonzero(lr[:-1]*lr[1:]<0)]
    result = dict(mass=float(m[-1]),r99=float(r99),Noether_N=float(N),omega=float(w),
                  central_lapse=float(A[0]),min_g=float(g.min()),D_max=float(-peak.fun),
                  D_peak_radius=float(peak.x),D_local_maxima=len(peaks),
                  grid_max_deficit=float(-peak.fun-dr.max()),
                  quadrature_D_differences=quad_errors,
                  angular_Einstein_residual=angular_error,ODE_pointwise_residual=ode_error,
                  collocation_RMS_residual=float(np.max(sol.rms_residuals)),nodes=int(sol.x.size),
                  boundary_residual=float(np.max(abs(boundary(sigma, (.06/.08 if sigma==.08 else .044/.06),r0,rend)(sol.y[:,0],sol.y[:,-1],sol.p)))),
                  r0=float(r0),rend=float(rend),scalar_tail=float(sigma*f[-1]),
                  conditional_tail_mass_bound=float(tail_mass_bound),conditional_tail_D_upper=float(tail_D_upper),
                  light_rings=light_rings)
    return result, dd


def study(output, profiles):
    evidence = {'date':'2026-10-12','status':'numerical_endpoint_rejection',
                'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
                'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                'numerical_D_error_allowance':1e-4,
                'error_status':'convergence-tested numerical allowance, not interval-certified',
                'models':{}}
    profiles.mkdir(parents=True,exist_ok=True)
    for model in ('S08A06','S06A044'):
        base,sigma = solve_model(model)
        fine,_ = solve_model(model,tol=1e-9,seed=base)
        runs = [('base',base),('fine',fine)]
        for name,settings in [('short_domain',dict(rend=70.,r0=1e-3,mesh=901)),
                              ('long_domain',dict(rend=110.,r0=1e-5,mesh=1801))]:
            sol,_ = solve_model(model,tol=1e-10,seed=fine,**settings)
            runs.append((name,sol))
        reports = {}
        for name,sol in runs:
            report,dd = diagnose(sol,sigma)
            reports[name] = report
            assert report['D_max']+1e-4 < 0
            assert report['conditional_tail_D_upper']+1e-4 < 0
            assert report['angular_Einstein_residual'] < 1e-6
            assert len(report['light_rings'])==2
            print(model,name,json.dumps({k:report[k] for k in ('mass','r99','Noether_N','D_max','angular_Einstein_residual')}),flush=True)
        ref = reports['long_domain']
        spreads = {key:max(abs(report[key]-ref[key]) for report in reports.values())
                   for key in ('mass','r99','Noether_N','D_max','D_peak_radius','central_lapse')}
        assert spreads['D_max'] < 1e-6
        table = dict(mass=5.65,r99=17.9,Noether_N=26.1) if model=='S08A06' else dict(mass=10.43,r99=31.7,Noether_N=65.9)
        for key,value in table.items():
            assert abs(ref[key]-value) < (.005 if key=='mass' else .05)
        rs = np.linspace(runs[-1][1].x[0],110.,2201)
        ys = runs[-1][1].sol(rs)
        _,dd = diagnose(runs[-1][1],sigma)
        path = profiles/(model+'.csv')
        with path.open('w',encoding='utf-8',newline='\n') as stream:
            np.savetxt(stream,np.column_stack((rs,ys[0],np.exp(ys[1]),sigma*ys[2],sigma*ys[3],dd(rs))),
                       delimiter=',',header='r,mass,lapse,scalar,scalar_derivative,D',comments='',fmt='%.15e')
        evidence['models'][model] = {'runs':reports,'refinement_spreads':spreads,'published_rounded_values':table,
                                    'profile_path':path.as_posix(),'profile_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    output.write_text(json.dumps(evidence,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: both selected backgrounds fail the zero-energy static endpoint test.',flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--study',action='store_true')
    ap.add_argument('--profiles',type=Path,default=Path('evidence/smooth-scalar-profiles'))
    ap.add_argument('--model', default='S08A06', choices=['S08A06','S06A044'])
    ap.add_argument('--wall', type=float)
    ap.add_argument('--nuc', type=float, default=-1.2)
    ap.add_argument('--tol', type=float, default=1e-6)
    ap.add_argument('--rend', type=float, default=90.)
    ap.add_argument('--r0', type=float, default=1e-4)
    ap.add_argument('--seed', type=Path)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    if args.study:
        study(args.out,args.profiles)
        return
    sigma, central, mass, radius = (.08,.06,5.65,17.9) if args.model=='S08A06' else (.06,.044,10.43,31.7)
    r = np.linspace(args.r0,args.rend,1201)
    if args.seed:
        old = np.load(args.seed)
        y = np.array([np.interp(r,old['r'],row) for row in old['y']])
        p = old['p']
    else:
        y,p = guess(r,sigma,central,mass,args.wall or radius-3,args.nuc)
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        sol = solve_bvp(equations(sigma),boundary(sigma,central/sigma,args.r0,args.rend),
                        r,y,p=p,tol=args.tol,max_nodes=30000)
    print(args.model,sol.status,sol.message,'nodes',sol.x.size,'p',sol.p,flush=True)
    if not sol.success:
        raise SystemExit(1)
    rr = np.linspace(args.r0,args.rend,18001)
    yy = sol.sol(rr)
    m,nu,f,h = yy
    g = 1-2*m/rr
    A = np.exp(nu)
    I = A[0]*args.r0 + cumulative_trapezoid(A/np.sqrt(g),rr,initial=0)
    D = I*np.sqrt(g)*equations(sigma)(rr,yy,sol.p)[1]/A-1
    N = np.trapezoid(4*np.pi*rr*rr*expit(sol.p[0])*sigma*sigma*f*f/(A*np.sqrt(g)),rr)
    print('M',m[-1],'r99',np.interp(.99*m[-1],m,rr),'N',N,
          'omega',expit(sol.p[0]),'lapse0',A[0],'min_g',g.min(),
          'Dmax',D.max(),'rDmax',rr[D.argmax()],
          'residual',np.max(sol.rms_residuals),'tail',f[-1],flush=True)
    np.savez(args.out,r=sol.x,y=sol.y,p=sol.p)


if __name__=='__main__':
    main()
