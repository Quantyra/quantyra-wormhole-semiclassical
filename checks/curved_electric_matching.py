"""S074: common zeta/spectral subtraction and the massive constant-field response."""
import argparse
import hashlib
import json
import platform
from pathlib import Path

import mpmath as mp
import sympy as s


def zero(expr):
    assert s.simplify(expr)==0, expr


def analytic_checks():
    M=s.symbols('M',positive=True)
    A=(M*s.polygamma(1,M)-1-1/(2*M))/(2*s.pi)
    r2=-(s.polygamma(1,M)/M+s.polygamma(2,M)+1/(2*M**3))/(2*s.pi)
    zero(s.diff(A,M)/(2*M)+r2/2)
    n,q,b,a=s.symbols('n q b a',positive=True)
    zero(2*(q+2*n)/(n*(n+q))-2/n-2/(n+q))
    zero(2*(q+2*n)/(n**2*(n+q)**2)-2/q*(1/n**2-1/(n+q)**2))
    xx=s.symbols('x',nonnegative=True)
    first=s.factor(1-(1+xx)**-2)
    rem=s.factor((1+xx)**-2-(1-2*xx))
    top=s.factor(3*xx**2-rem)
    zero(first-xx*(xx+2)/(xx+1)**2)
    zero(rem-xx**2*(2*xx+3)/(xx+1)**2)
    zero(top-xx**3*(3*xx+4)/(xx+1)**2)
    constants=[4*s.zeta(2)/(2*s.pi)**2,
               8*s.zeta(4)/(2*s.pi)**4,12*s.zeta(6)/(2*s.pi)**6]
    for value,target in zip(constants,[s.Rational(1,6),s.Rational(1,180),s.Rational(1,5040)]):
        zero(value-target)
    laplace=[constants[j]*s.factorial(2*j+1) for j in range(3)]
    assert laplace==[s.Rational(1,6),s.Rational(1,30),s.Rational(1,42)]
    # Independent Lorentzian tensor contractions with both parallel fields.
    x,y,ee,bb,kap=s.symbols('x y E B kappa2',real=True)
    cs=s.symbols('c1:10')
    c1,c2,c3,c4,c5,c6,c7,c8,c9=cs
    signs=[-1,1,1,1]; eta=s.diag(*signs)
    F=s.zeros(4); F[0,1]=ee; F[1,0]=-ee; F[2,3]=bb; F[3,2]=-bb
    Fu=eta*F*eta; X=F*eta*F.T
    Ric=s.diag(x,-x,y,y)
    I=sum(F[i,j]*Fu[i,j] for i in range(4) for j in range(4))
    ricX=sum(signs[i]*signs[j]*Ric[i,j]*X[i,j] for i in range(4) for j in range(4))
    def R(i,j,k,l):
        curv=-x if all(v<2 for v in [i,j,k,l]) else y if all(v>=2 for v in [i,j,k,l]) else 0
        return curv*(eta[i,k]*eta[j,l]-eta[i,l]*eta[j,k])
    rFF=sum(R(i,j,k,l)*Fu[i,j]*Fu[k,l] for i in range(4) for j in range(4) for k in range(4) for l in range(4))
    I4=s.trace((F*Fu)**2)
    for actual,target in [(I,2*(bb**2-ee**2)),(I4,2*(bb**4+ee**4)),
                          (ricX,2*x*ee**2+2*y*bb**2),(rFF,4*x*ee**2+4*y*bb**2)]:
        zero(actual-target)
    L=kap*(c4*2*(y-x)*I+c5*ricX+c6*rFF)+kap**2*(c7*I**2+c8*I4)
    theta=s.diff(L,ee,2).subs(ee,0)
    K=2*c4+c5+2*c6
    zero(theta-(4*kap*(K*x-2*c4*y)-16*kap**2*c7*bb**2))
    cj=s.symbols('cJ')
    zero(theta.subs({c4:0,c5:0,c6:0,c7:-cj/8,c8:cj/4})-2*kap**2*cj*bb**2)
    # Exact box and family bounds, independent of numerical candidate samples.
    lo=s.Rational(99,100); alpha=s.Rational(1,4)
    assert 4/lo<5 and 8*alpha/lo<3 and 16*(1+alpha)/lo**2<21
    flat_upper=s.Rational(36,108)/s.Integer(2000)**8
    curve_upper=s.Rational(16,540)/s.Integer(2000)**9
    assert flat_upper+curve_upper<s.Rational(2,10**27)
    assert s.Rational(1,200)+s.Rational(2,10**27)<s.Rational(1,100)
    lam_upper=48*s.Rational(22,7)*(1+s.Rational(1,10**9))/s.Integer(2000)**10
    assert lam_upper*49*(2*10**24)<s.Rational(1,200)
    return {'mass_derivative_identity':'PASS','sphere_telescope_identities':'PASS',
            'positive_rational_remainders':[str(first),str(rem),str(top)],
            'Laplace_coefficients':list(map(str,laplace)),
            'local_tensor_and_J_squared_regression':'PASS',
            'uniform_loop_fraction_bound':'2e-27 for all k>=2000 on the one-percent reference box',
            'flat_fraction_rational_upper':str(flat_upper),'curvature_fraction_rational_upper':str(curve_upper),
            'conditional_local_range':'lambda*NE<=1/200 together with the preceding reference condition',
            'scope':'fixed-metric covariantly constant field; the massless LLL and physical boundary susceptibility are excluded'}


def cm_dimensionless(M):
    return (M*mp.polygamma(1,M)-1-1/(2*M))/mp.pi


def resolvent_closed(M,beta):
    nu=mp.sqrt(M*M+beta*beta)
    return -(mp.digamma(nu+beta)+mp.digamma(nu-beta)-2*mp.digamma(M)+(nu-M)/(M*M))/(2*mp.pi)


def resolvent_integral(M,beta):
    sine2=mp.sin(mp.pi*beta)**2
    def integrand(v):
        rho0=v*mp.coth(mp.pi*v)/mp.pi if v else 1/mp.pi**2
        delta_rho=-rho0*sine2/(mp.sinh(mp.pi*v)**2+sine2)
        left=v*v+M*M+beta*beta; right=v*v+M*M
        return delta_rho/left-rho0*beta*beta/(left*right)
    cuts=sorted(set([mp.mpf(0),beta,10*beta,mp.mpf(1),mp.mpf(5)]))+[mp.inf]
    return beta/(2*mp.pi*M*M)+mp.quad(integrand,cuts)


def F_laplace(t):
    if abs(t)<mp.mpf('.001'):
        return sum(mp.bernoulli(2*j)*t**(2*j-1)/mp.factorial(2*j-1) for j in range(1,13))
    return mp.coth(t/2)/2-t/(4*mp.sinh(t/2)**2)


def one_channel_checks():
    resolvents=[]
    for M in map(mp.mpf,['.3','1','4','18.05']):
        for beta in map(mp.mpf,['.07','.23']):
            direct=resolvent_integral(M,beta)
            closed=resolvent_closed(M,beta)
            error=abs(direct-closed)
            assert error<mp.mpf('1e-45'),(M,beta,error)
            resolvents.append({'M':str(M),'beta':str(beta),'direct_resolvent':str(direct),'absolute_error':str(error)})
    laplace=[]
    for M in map(mp.mpf,['.25','1','4','18.05','100','1000']):
        integral=mp.quad(lambda v:mp.exp(-v)*F_laplace(v/M)/M,[0,1,5,20,mp.inf])/mp.pi
        closed=cm_dimensionless(M)
        error=abs(integral-closed)
        assert error<mp.mpf('1e-45'),(M,error)
        flat=1/(6*mp.pi*M*M); bound=1/(30*mp.pi*M**4)
        next_bound=1/(42*mp.pi*M**6)
        assert 0<closed<=flat and 0<=flat-closed<=bound
        assert 0<=closed-flat+bound<=next_bound
        laplace.append({'M':str(M),'coefficient_over_a_squared':str(closed),'absolute_error':str(error),
                        'curvature_error_over_bound':str((flat-closed)/bound)})
    return {'direct_resolvent_cases':resolvents,'independent_Laplace_cases':laplace,
            'quadratures_are_not_interval_certificates':True}


def spectral_zeta_near_one(t,q,b):
    N=max(64,16*q)
    value=sum(2*(q+2*n)/(mp.mpf(n)*(n+q))**(1+t) for n in range(1,N+1))
    for j in range(41):
        coeff=(-1)**j*mp.rf(1+t,j)/mp.factorial(j)
        if j:
            coeff+=mp.mpf('.5')*(-1)**(j-1)*mp.rf(1+t,j-1)/mp.factorial(j-1)
        value+=4*coeff*mp.mpf(q)**j*mp.zeta(1+2*t+j,N+1)
    return b**(2+2*t)*value


def scheme_checks():
    rows=[]; b=mp.mpf('1.7'); mu=mp.mpf('.3')
    for q in [1,3,8,64]:
        exact=b*b/(3*mp.pi)*(2*mp.log(b*mu)+2*mp.euler-mp.harmonic(q))
        def response_zeta(t):
            return t*mu**(2*t)*spectral_zeta_near_one(t,q,b)/(6*mp.pi)
        h=mp.mpf('0.00001')
        d1=(response_zeta(h)-response_zeta(-h))/(2*h)
        d2=(response_zeta(2*h)-response_zeta(-2*h))/(4*h)
        derivative=(4*d1-d2)/3
        error=abs(derivative-exact)
        assert error<mp.mpf('1e-12'),(q,error)
        cutoff=[]
        for N in [100,1000,10000]:
            finite=b*b/(3*mp.pi)*(mp.harmonic(N)+mp.harmonic(N+q)-mp.harmonic(q)
                       -2*mp.log(mp.sqrt(N*(N+q))/(b*mu)))
            bound=b*b/(6*mp.pi)*(mp.mpf(1)/N+mp.mpf(1)/(N+q))
            assert 0<finite-exact<bound
            cutoff.append({'N':N,'positive_error':str(finite-exact),'bound':str(bound)})
        rows.append({'q':q,'b':str(b),'mu':str(mu),'common_coefficient':str(exact),
                     'independent_zeta_derivative_error':str(error),'spectral_cutoffs':cutoff})
    return {'cases':rows,'zeta_numerics':'paired finite difference and binomial/Hurwitz tail through j=40 with q/(N+1)<1/16; analytic finite part established separately',
            'tolerance':'1e-12','finite_F_squared_scheme_conversion':0,
            'inverse_g_squared_log_mu_derivative':'-1/(6*pi^2)'}


def tower_checks():
    rows=[]
    for q,C in [(1,mp.mpf('.2')),(3,mp.mpf(1)),(8,mp.mpf(10)),(64,mp.mpf(100))]:
        # a=1, b^2=C. Sum a finite tower independently of the harmonic bound.
        N=128
        partial=mp.mpf(0)
        for n in range(1,N+1):
            M=mp.sqrt(mp.mpf(n)*(n+q)/C)
            partial+=2*(q+2*n)*(cm_dimensionless(M)-1/(6*mp.pi*M*M))
        H2=mp.zeta(2)-mp.polygamma(1,q+1)
        bound=C*C*H2/(15*mp.pi*q)
        tail=C*C/(15*mp.pi*q)*(mp.polygamma(1,N+1)-mp.polygamma(1,N+q+1))
        assert -bound<=partial<0
        # Bound the leading-curvature remainder in the finite sum directly.
        exact_sum4=2*C*C/q*(H2-mp.polygamma(1,N+1)+mp.polygamma(1,N+q+1))
        rem=partial+exact_sum4/(30*mp.pi)
        assert 0<=rem<=C**3*H2/(21*mp.pi*q*(q+1))
        rows.append({'q':q,'C':str(C),'levels':N,'partial_curvature_coefficient':str(partial),
                     'absolute_angular_tail_bound':str(tail),'absolute_full_curvature_bound':str(bound)})
    return rows


def candidates():
    rows=[]
    for kval in [2000,4000,20000]:
        k=mp.mpf(kval); q=k**4; g2=k**-8
        C=mp.sqrt(1+12*mp.pi*k**6)
        H2=mp.zeta(2)-mp.polygamma(1,q+1)
        flat=g2*(2*mp.euler-mp.harmonic(q))/(12*mp.pi**2)
        bound=g2*C*H2/(60*mp.pi**2*q)
        remainder=bound*5*C/(7*(q+1))
        lam=48*mp.pi*C/(k**10*(C-1))
        assert abs(flat)+bound<mp.mpf('2e-27')
        assert remainder/bound<mp.mpf('.0022')
        rows.append({'k':kval,'M1':str(mp.sqrt((q+1)/C)),
                     'flat_loop_fraction':str(flat),
                     'curvature_fraction_lower':str(-bound),
                     'curvature_fraction_upper':str(-bound+remainder),
                     'leading_curvature_relative_remainder_bound':str(remainder/bound),
                     'NE_limit':str(1/(200*lam))})
    return rows


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(); mp.mp.dps=70
    root=Path(__file__).resolve().parents[1]
    deps=['notes/uniform-dirac-reference.md','evidence/uniform-dirac-reference.json',
          'notes/massive-bulk-polarization.md','evidence/massive-bulk-polarization.json',
          'notes/finite-local-reference-window.md','notes/finite-local-electromagnetic-input.md']
    hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in deps}
    print('Checking exact mass, tensor, remainder and uniform family identities',flush=True)
    result={'analytic':analytic_checks()}
    print('Comparing continuum/discrete resolvents and independent Laplace coefficients',flush=True)
    result['one_channel']=one_channel_checks()
    print('Comparing zeta finite part, spectral subtraction and complete-tower bounds',flush=True)
    result['scheme_matching']=scheme_checks(); result['tower_checks']=tower_checks()
    result['candidate_family']=candidates()
    result.update({'status':'PASS','dependencies_sha256':hashes,
                   'scope':'one-loop massive constant-field bulk electric coefficient in the reference scheme; no LLL local replacement or physical boundary/drive matching',
                   'versions':{'python':platform.python_version(),'sympy':s.__version__,'mpmath':mp.__version__},
                   'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('PASS: common electric subtraction, exact AdS coefficient and signed tower bound',flush=True)


if __name__=='__main__':
    main()
