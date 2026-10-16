"""S133 exact cutoff action, geometric checks and uniform C2 estimates."""
import hashlib
import json
from pathlib import Path
import platform
import sys

import mpmath as mp
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT/'notes/cutoff-clock-contact-matching.md'
OUTPUT = ROOT/'evidence/cutoff-clock-contact-matching.json'
DPS = 50
FD_STEPS = ('2e-4','1e-4')
FIRST_TOL = '1e-8'
SECOND_TOL = '1e-6'
GEOMETRY_TOL = '1e-35'
BOUND_ALLOWANCE = '1e-30'


def zero(value):
    if s.simplify(value) != 0:
        raise AssertionError(str(value))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_geometry():
    time, depth, re, ep = s.symbols('T Z r_e epsilon', positive=True)
    x,p,y = s.symbols('x P y', positive=True)
    f2,p2,rate = s.symbols('f2 p2 theta_prime', real=True)
    metric = s.diag(-re**2/depth**2,re**2/depth**2)
    inverse = metric.inv()
    coordinates = (time,depth)
    gamma = [[[s.simplify(sum(inverse[a,d]*(s.diff(metric[d,c],coordinates[b])+s.diff(metric[d,b],coordinates[c])-s.diff(metric[b,c],coordinates[d]))/2 for d in range(2))) for c in range(2)] for b in range(2)] for a in range(2)]
    for a in range(2):
        for b in range(2):
            for c in range(2):
                expected = -1/depth if (a,b,c) in ((0,0,1),(0,1,0),(1,0,0),(1,1,1)) else 0
                zero(gamma[a][b][c]-expected)
    tangent = s.Matrix([x,ep*y])
    second = s.Matrix([f2,ep*p2])
    delta = x**2-ep**2*y**2
    lapse2 = re**2*delta/(ep**2*p**2)
    local_metric = metric.subs(depth,ep*p)
    zero((tangent.T*local_metric*tangent)[0]+lapse2)
    # Evaluate normal/curvature using the Christoffels derived above.
    normal = -ep*p/(re*s.sqrt(delta))*s.Matrix([ep*y,x])
    normal_cov = local_metric*normal
    zero((normal.T*local_metric*normal)[0]-1)
    zero((normal_cov.T*tangent)[0])
    acceleration = second+s.Matrix([sum(gamma[a][b][c].subs(depth,ep*p)*tangent[b]*tangent[c] for b in range(2) for c in range(2)) for a in range(2)])
    christoffel_K = s.simplify((normal_cov.T*acceleration)[0]/lapse2)
    numerator = x*delta+ep**2*p*(y*f2-x*p2)
    expected_K = numerator/(re*delta**s.Rational(3,2))
    zero(christoffel_K-expected_K)
    zero(expected_K.subs({y:0,f2:0,p2:0})-1/re)
    # Product sphere metric is independent of T,Z: angular K is zero.
    angular_metric = re**2
    zero(s.diff(angular_metric,time))
    zero(s.diff(angular_metric,depth))
    lapse = re*s.sqrt(delta)/(ep*p)
    zeta,c,d = s.symbols('zeta_s c_theta d_theta', real=True)
    intrinsic_inverse = -1/lapse**2
    lhat_upup = s.Rational(2,3)*expected_K*intrinsic_inverse
    sphere_volume = 4*s.pi*lapse*re**2
    source_density = -zeta*sphere_volume*intrinsic_inverse*rate**2/2
    curvature_density = c*sphere_volume*lhat_upup*rate**2
    scheme_density = -d*sphere_volume*expected_K*intrinsic_inverse*rate**2/3
    exact = source_density+curvature_density+scheme_density
    A = p*rate**2
    compact = ep*(2*s.pi*zeta*re*A/s.sqrt(delta)+(-8*s.pi*c/3+4*s.pi*d/3)*A*numerator/delta**2)
    zero(exact-compact)
    # P=F' sets p=x,y=f2,p2=f3; it is only asymptotic fixed lapse.
    f3 = s.symbols('f3',real=True)
    special_K = expected_K.subs({p:x,y:f2,p2:f3})
    closed_K = (1-ep**2*f3/x)/(re*(1-ep**2*(f2/x)**2)**s.Rational(3,2))
    # Sympy square-root branch simplification is valid for x>0, delta>0.
    zero(special_K-closed_K)
    curvature_series = s.series(closed_K,ep,0,4).removeO()
    schwarzian = f3/x-s.Rational(3,2)*(f2/x)**2
    zero(curvature_series-(1-ep**2*schwarzian)/re)
    special_lapse = re/ep*s.sqrt(1-ep**2*(f2/x)**2)
    weight_series = s.series(re**2*closed_K/special_lapse,ep,0,5).removeO()
    zero(weight_series-(ep-ep**3*(f3/x-2*(f2/x)**2)))
    # Fixed-lapse embedding equation through its first correction.
    u=s.symbols('u',real=True)
    F=s.Function('F')(u)
    p_approx=s.diff(F,u)-ep**2*s.diff(F,u,2)**2/(2*s.diff(F,u))
    residual=p_approx**2+ep**2*s.diff(p_approx,u)**2-s.diff(F,u)**2
    zero(s.series(residual,ep,0,4).removeO())
    return {'Christoffels':'derived from g_2; four nonzero entries equal -1/Z',
            'exact_K':str(expected_K),'exact_total_density':str(compact),
            'P_equals_Fprime_K':str(closed_K),
            'P_equals_Fprime_lapse':str(special_lapse),
            'curvature_series':str(curvature_series),
            'weight_series':str(weight_series),
            'exact_fixed_lapse_condition':'P^2+epsilon^2*Pprime^2=Fprime^2'}


def product_bounds(X,Y):
    return (X[0]*Y[0],X[1]*Y[0]+X[0]*Y[1],X[2]*Y[0]+2*X[1]*Y[1]+X[0]*Y[2])


def uniform_bounds(B,delta_lower,re,zeta,c,d):
    A=(B**5,5*B**5,22*B**5)
    M=(3*B**3,11*B**2,24*B)
    def power(alpha):
        return (delta_lower**(-alpha),4*alpha*B*delta_lower**(-alpha-1),16*alpha*(alpha+1)*B**2*delta_lower**(-alpha-2)+4*alpha*delta_lower**(-alpha-1))
    T=product_bounds(A,power(mp.mpf('0.5')))
    G=product_bounds(product_bounds(A,M),power(mp.mpf(2)))
    return tuple(2*mp.pi*abs(zeta)*re*T[j]+(8*mp.pi*abs(c)/3+4*mp.pi*abs(d)/3)*G[j] for j in range(3))


def smooth_families():
    u,a,b,ep=s.symbols('u a b epsilon',real=True)
    h=s.sin(3*u)/50+u**2/100
    j=s.cos(u)/50
    k=s.cos(2*u)/50-u/100
    ell=s.sin(u)/50+u**2/100
    entries=[
        ('sinusoidal',u+u**2/25+s.sin(u)/50,1+7*s.cos(2*u)/100+3*u**2/100),
        ('polynomial',u+3*u**3/100+u**2/50,1+u/25+u**2/50),
    ]
    for name,F0,P0 in entries:
        F=F0+a*h+b*k
        P=P0+a*j+b*ell
        x,y,f2,p2=s.diff(F,u),s.diff(P,u),s.diff(F,u,2),s.diff(P,u,2)
        # Smooth global-time source; its u-rate is varied with F.
        q=s.Rational(1,5)+s.cos(2*F)/25+3*F**2/100
        rate=q*x
        delta=x**2-ep**2*y**2
        numerator=x*delta+ep**2*P*(y*f2-x*p2)
        re,zeta,c,d=s.Rational(13,10),s.Rational(3,10),s.Rational(1,5),1/(60*s.pi**2)
        compact=ep*(2*s.pi*zeta*re*P*rate**2/s.sqrt(delta)+(-8*s.pi*c/3+4*s.pi*d/3)*P*rate**2*numerator/delta**2)
        # Differentiate the exact density, not its epsilon expansion.
        analytic=[expr.subs({a:0,b:0}) for expr in (compact,s.diff(compact,a),s.diff(compact,a,2),s.diff(compact,a,b))]
        analytic_fns=[s.lambdify((u,ep),expr,'mpmath',cse=True) for expr in analytic]
        jets_fn=s.lambdify((u,a,b),(F,P,x,y,f2,p2,q,rate),'mpmath',cse=True)
        variations=s.lambdify(u,(h,s.diff(h,u),s.diff(h,u,2),j,s.diff(j,u),s.diff(j,u,2),k,s.diff(k,u),s.diff(k,u,2),ell,s.diff(ell,u),s.diff(ell,u,2)),'mpmath',cse=True)
        # Bounds on an explicit |a|,|b|<=1 neighborhood, u in [0,1]:
        # |F'|>=1-.08-.05>.8, P>=.88>.8; all needed jets,
        # Theta',Theta'',Theta''' are bounded by B=2.
        yield name,jets_fn,analytic_fns,variations


def numerical_checks():
    mp.mp.dps=DPS
    first_tol,second_tol=mp.mpf(FIRST_TOL),mp.mpf(SECOND_TOL)
    geometry_tol,allowance=mp.mpf(GEOMETRY_TOL),mp.mpf(BOUND_ALLOWANCE)
    re,zeta,c,d=mp.mpf(13)/10,mp.mpf(3)/10,mp.mpf(1)/5,1/(60*mp.pi**2)
    B,v,p_lower=mp.mpf(2),mp.mpf(4)/5,mp.mpf(4)/5
    delta_lower=v*v/2
    C=uniform_bounds(B,delta_lower,re,zeta,c,d)
    results=[]
    largest_first=mp.mpf(0)
    largest_second=mp.mpf(0)
    for name,jets,analytic,variations in smooth_families():
        for ep in (mp.mpf(1)/5,mp.mpf(1)/10,mp.mpf(1)/20):
            max_geometry_error=mp.mpf(0)
            def direct_density(u,a,b):
                _,P,x,y,f2,p2,_,rate=jets(u,a,b)
                depth=ep*P
                delta=x*x-(ep*y)**2
                if P<=0 or x<=0 or delta<=0:
                    raise AssertionError('curve left timelike positive branch')
                lapse=re*mp.sqrt(delta)/depth
                # Geometric reconstruction via normal and Christoffel acceleration.
                normal_T_cov=re*ep*y/(depth*mp.sqrt(delta))
                normal_Z_cov=-re*x/(depth*mp.sqrt(delta))
                accel_T=f2-2*x*ep*y/depth
                accel_Z=ep*p2-(x*x+(ep*y)**2)/depth
                K=(normal_T_cov*accel_T+normal_Z_cov*accel_Z)/lapse**2
                return 2*mp.pi*zeta*re**2*rate**2/lapse+(-8*mp.pi*c/3+4*mp.pi*d/3)*re**2*K*rate**2/lapse
            def action(a,b):
                return mp.quad(lambda u:direct_density(u,a,b),[0,mp.mpf('0.5'),1])
            expected=[mp.quad(lambda u,fn=fn:fn(u,ep),[0,mp.mpf('0.5'),1]) for fn in analytic]
            center=action(0,0)
            if abs(center-expected[0])>geometry_tol:
                raise AssertionError('integrated geometric density disagreement')
            # These grid checks supplement the analytic neighborhood estimates.
            for i in range(21):
                u=mp.mpf(i)/20
                F,P,x,y,f2,p2,q,rate=jets(u,0,0)
                delta=x*x-ep**2*y*y
                if x<v or P<p_lower or delta<delta_lower:
                    raise AssertionError('declared geometric lower bound failed')
                if max(abs(P),abs(x),abs(y),abs(f2),abs(p2),abs(q))>B:
                    raise AssertionError('declared jet bound failed')
                # Theta''=-2 sin(2F)/25+3F/50; Theta'''=-4 cos(2F)/25+3/50.
                if max(abs(-2*mp.sin(2*F)/25+3*F/50),abs(-4*mp.cos(2*F)/25+mp.mpf(3)/50))>B:
                    raise AssertionError('declared source jet bound failed')
                point_error=abs(direct_density(u,0,0)-analytic[0](u,ep))
                max_geometry_error=max(max_geometry_error,point_error)
                if point_error>geometry_tol:
                    raise AssertionError('pointwise independent geometry disagreement')
                vals=variations(u)
                H=max(abs(x) for x in vals[:6])
                Jv=max(abs(x) for x in vals[6:])
                # Use whole-interval certified simple direction upper bounds 1.
                if H>1 or Jv>1:
                    raise AssertionError('direction bound failed')
                for index,Cindex in ((0,0),(1,1),(2,2),(3,2)):
                    if abs(analytic[index](u,ep))>ep*C[Cindex]+allowance:
                        raise AssertionError('pointwise uniform action/variation bound failed')
            step_records=[]
            for step in FD_STEPS:
                h=mp.mpf(step)
                plus,minus=action(h,0),action(-h,0)
                fd1=(plus-minus)/(2*h)
                fd2=(plus-2*center+minus)/h**2
                mixed=(action(h,h)-action(h,-h)-action(-h,h)+action(-h,-h))/(4*h*h)
                errors=[abs(fd1-expected[1]),abs(fd2-expected[2]),abs(mixed-expected[3])]
                largest_first=max(largest_first,errors[0])
                largest_second=max(largest_second,errors[1],errors[2])
                if errors[0]>first_tol or max(errors[1:])>second_tol:
                    raise AssertionError('finite difference derivative disagreement')
                step_records.append({'step':step,'first_error':str(errors[0]),'second_error':str(errors[1]),'mixed_error':str(errors[2])})
            if abs(expected[0])>ep*C[0]+allowance or abs(expected[1])>ep*C[1]+allowance or max(abs(expected[2]),abs(expected[3]))>ep*C[2]+allowance:
                raise AssertionError('integrated uniform bound failed')
            results.append({'family':name,'epsilon':str(ep),'action':str(expected[0]),'first_variation':str(expected[1]),'second_variation':str(expected[2]),'mixed_variation':str(expected[3]),'action_over_epsilon':str(expected[0]/ep),'geometry_max_error':str(max_geometry_error),'finite_differences':step_records})
    return {'synthetic_cases':results,'largest_first_difference_error':str(largest_first),'largest_second_difference_error':str(largest_second),'uniform_bound_example':{'B':str(B),'v':str(v),'P_lower':str(p_lower),'Delta_lower':str(delta_lower),'interval_length':'1','direction_norm_upper_bounds':['1','1'],'C0_C1_C2':list(map(str,C)),'meaning':'bounds are epsilon*interval_length*Cj times direction norms; no physical coefficient assignment or clock-solution error inferred'}}


def main():
    symbolic=exact_geometry()
    numeric=numerical_checks()
    sources=['notes/source-clock-radiation.md','notes/curved-fermion-boundary.md','notes/source-clock-joint.md','notes/boundary-vertex-matching.md','notes/cutoff-clock-contact-matching.md','checks/cutoff_clock_contact_matching.py']
    receipt={'status':'pass','milestone':'E003/S133','symbolic':symbolic,'numerical':numeric,'predeclared_tolerances':{'dps':DPS,'steps':list(FD_STEPS),'first_absolute':FIRST_TOL,'second_absolute':SECOND_TOL,'geometry_absolute':GEOMETRY_TOL,'bound_rounding_allowance':BOUND_ALLOWANCE},'claims':{'uniform_action_first_second_variations_O_epsilon':True,'zeta_s_throat_cut_also_suppressed':True,'fixed_lapse_gauge_solved':False,'exterior_boundary_matching_supplied':False,'finite_coefficients_physically_assigned':False,'full_Weyl_anomaly_computed':False,'S074_complete':False,'new_Lambda_interval':False,'area_budget_evaluated':False,'physical_opening_achieved':False},'inputs_sha256':{name:sha(ROOT/name) for name in sources},'runtime':{'python':sys.version,'platform':platform.platform(),'sympy':s.__version__,'mpmath':mp.__version__}}
    OUTPUT.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
    print('S133 exact cutoff contact matching: PASS')
    print('Independent Christoffel geometry, exact action and fixed-lapse distinction passed.')
    print('Six smooth finite-cut cases and two finite-difference scales passed.')
    print('Largest first-difference error: '+numeric['largest_first_difference_error'])
    print('Largest second-difference error: '+numeric['largest_second_difference_error'])
    print('Finite bounded throat-cut contacts and their C2 clock variations are O(epsilon).')
    print('Exterior matching and the full S074 physical gate remain open.')
    print('Receipt: '+str(OUTPUT))


if __name__=='__main__':
    main()
