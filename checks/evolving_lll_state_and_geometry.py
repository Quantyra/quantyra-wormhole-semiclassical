"""S137.2 exact free-LLL transport; no physical geometry certificate."""
from pathlib import Path
import datetime
import hashlib
import json
import platform
import subprocess
import sys
import traceback

import mpmath as mp
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'evidence/evolving-lll-state-and-geometry.json'
LOG=ROOT/'evidence/evolving-lll-state-and-geometry.log'
DPS=60
TOL={'exact_numerical':'1e-35','finest_point_split':'2e-10','missing_cocycle_discrepancy_min':'1e-7'}
SEPARATIONS=['1e-4','1e-5','1e-6']
SOURCES=['notes/evolving-lll-state-and-geometry.md','checks/evolving_lll_state_and_geometry.py',
         'notes/lll-state-difference-source.md','notes/retimed-returned-quantum-stress.md',
         'notes/dirac-spherical-source-reduction.md','notes/spherical-anomaly-action.md',
         'notes/spherical-anomaly-sources.md','notes/physical-clock-remainder.md',
         'notes/exterior-lll-delayed-boundary.md']


def write(path,value):
    with path.open('w',encoding='utf-8',newline='\n') as f:
        f.write(value)


def record(value):
    with LOG.open('a',encoding='utf-8',newline='\n') as f:
        f.write(value+'\n')
    print(value,flush=True)


def sha(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def close(error):
    assert abs(error)<mp.mpf(TOL['exact_numerical']),str(error)


def symbolic_point_split_and_covariance():
    d,h1,h2,h3,a1,N=sp.symbols('d h1 h2 h3 a1 N',positive=True)
    delta=h1*d+h3*d**3/24
    pref=h1+d*d*(h3/8-h2*h2/(8*h1))
    covariance=pref*(sp.cos(a1*d)/(4*sp.pi*sp.I*sp.sin(delta/2))-sp.I*sp.sin(a1*d)/(4*sp.pi*sp.cos(delta/2)))
    coefficient=sp.expand(sp.series(covariance,d,0,2).removeO()).coeff(d,1)
    measured=-sp.I*N*coefficient/2
    schwarz=h3/h1-sp.Rational(3,2)*(h2/h1)**2
    target=N*(-h1*h1/(96*sp.pi)-schwarz/(48*sp.pi)+(a1*a1-h1*a1)/(8*sp.pi))
    assert sp.simplify(measured-target)==0
    rp,rpp,fp,fpp,fppp=sp.symbols('rp rpp fp fpp fppp',nonzero=True)
    tilde_rp=fp*rp+fpp/(2*fp)
    tilde_rpp=fp*fp*rpp+fpp*rp+fppp/(2*fp)-fpp*fpp/(2*fp*fp)
    schf=fppp/fp-sp.Rational(3,2)*(fpp/fp)**2
    assert sp.simplify(-(tilde_rp**2-tilde_rpp)+fp**2*(rp**2-rpp)-schf/2)==0
    hp,hpp,hppp=sp.symbols('hp hpp hppp',nonzero=True)
    composed=[hp*fp,hpp*fp*fp+hp*fpp,hppp*fp**3+3*hpp*fp*fpp+hp*fppp]
    schcomp=composed[2]/composed[0]-sp.Rational(3,2)*(composed[1]/composed[0])**2
    assert sp.simplify(schcomp-fp*fp*(hppp/hp-sp.Rational(3,2)*(hpp/hp)**2)-schf)==0
    x=sp.symbols('x',real=True)
    R=sp.Function('R')(x)
    rho=sp.Function('rho')(x)
    psi=sp.Function('psi')(x)
    rescale=sp.exp(-rho/2)/R
    assert sp.simplify(sp.diff(rescale*psi,x)+(sp.diff(R,x)/R+sp.diff(rho,x)/2)*rescale*psi-rescale*sp.diff(psi,x))==0
    return dict(point_split_coefficient='exact zero',Weyl_cocycle='exact zero',Schwarzian_chain_rule='exact zero',
                canonical_radius_connection_cancellation='exact zero',absolute_four_dimensional_radius_source_inferred=False)


def symbolic_ward_and_system():
    p,m=sp.symbols('p m',real=True)
    N,G=sp.symbols('N G',positive=True)
    rho=sp.Function('rho')(p,m)
    X=sp.Function('X')(p,m)
    e=sp.exp(2*rho)
    metric=sp.Matrix([[0,-e/2],[-e/2,0]])
    inverse=metric.inv()
    coords=[p,m]
    gamma=[[[sp.simplify(sum(inverse[a,d]*(sp.diff(metric[d,b],coords[c])+sp.diff(metric[d,c],coords[b])-sp.diff(metric[b,c],coords[d])) for d in range(2))/2) for c in range(2)] for b in range(2)] for a in range(2)]
    rp,rm=sp.diff(rho,p),sp.diff(rho,m)
    rpm=sp.diff(rho,p,m)
    alpha=N/(24*sp.pi)
    tp,tm=sp.Function('tp')(p),sp.Function('tm')(m)
    T=sp.Matrix([[tp-alpha*(rp*rp-sp.diff(rho,p,2)),-alpha*rpm],[-alpha*rpm,tm-alpha*(rm*rm-sp.diff(rho,m,2))]])
    divergences=[]
    for b in range(2):
        div=sum(inverse[a,c]*(sp.diff(T[a,b],coords[c])-sum(gamma[d][c][a]*T[d,b]+gamma[d][c][b]*T[a,d] for d in range(2))) for a in range(2) for c in range(2))
        assert sp.simplify(div)==0
        divergences.append('exact zero')
    scalar=8*sp.exp(-2*rho)*rpm
    ricci=sp.Matrix([[sum(sp.diff(gamma[c][a][b],coords[c])-sp.diff(gamma[c][a][c],coords[b]) for c in range(2))
                     +sum(gamma[c][c][d]*gamma[d][a][b]-gamma[c][b][d]*gamma[d][a][c] for c in range(2) for d in range(2))
                     for b in range(2)] for a in range(2)])
    assert sp.simplify(sum(inverse[a,b]*ricci[a,b] for a in range(2) for b in range(2))-scalar)==0
    assert sp.simplify(sum(inverse[a,b]*T[a,b] for a in range(2) for b in range(2))-N*scalar/(48*sp.pi))==0
    state=sp.diag(tp,tm)
    assert sp.simplify(sum(inverse[a,b]*state[a,b] for a in range(2) for b in range(2)))==0
    V=sp.Function('V')(X)
    hessian=sp.Matrix([[sp.diff(X,coords[a],coords[b])-sum(gamma[c][a][b]*sp.diff(X,coords[c]) for c in range(2)) for b in range(2)] for a in range(2)])
    box=sum(inverse[a,b]*hessian[a,b] for a in range(2) for b in range(2))
    norm=sum(inverse[a,b]*sp.diff(X,coords[a])*sp.diff(X,coords[b]) for a in range(2) for b in range(2))
    E=sp.Matrix([[metric[a,b]*box-hessian[a,b]+(sp.diff(X,coords[a])*sp.diff(X,coords[b])-metric[a,b]*norm/2)/(2*X)-metric[a,b]*V/2 for b in range(2)] for a in range(2)])
    assert sp.simplify(E[0,1]-sp.diff(X,p,m)-e*V/4)==0
    assert sp.simplify(E[0,0]+sp.diff(X,p,2)-2*rp*sp.diff(X,p)-sp.diff(X,p)**2/(2*X))==0
    angular=scalar-box/X+norm/(2*X*X)+sp.diff(V,X)
    target_angular=2*rpm+sp.diff(X,p,m)/X-sp.diff(X,p)*sp.diff(X,m)/(2*X*X)+e*sp.diff(V,X)/4
    assert sp.simplify(e*angular/4-target_angular)==0
    xx,ee,v,vp,xp,xm,rest,J,kappa=sp.symbols('xx ee v vp xp xm rest J kappa',nonzero=True)
    r2,x2=sp.symbols('r2 x2')
    mixed=x2+ee*v/4+kappa*r2-2*G*rest
    angle=2*r2+x2/xx-xp*xm/(2*xx**2)+ee*vp/4+ee*J/4
    solution=sp.solve([mixed,angle],[r2,x2])
    expected=(ee*(v/xx-vp)/4+xp*xm/(2*xx**2)-2*G*rest/xx-ee*J/4)/(2-kappa/xx)
    assert sp.factor(solution[r2]-expected)==0
    assert sp.factor(solution[x2]+ee*v/4+kappa*expected-2*G*rest)==0
    # The finite state difference is traceless, so its physical-frame
    # radius chain rule has zero s and zero trace, without dropping W.
    state_trace=sp.simplify(sum(inverse[a,b]*state[a,b] for a in range(2) for b in range(2)))
    assert sp.simplify(2*G*(0-state_trace)/X)==0
    return dict(ward_divergences=divergences,trace_anomaly='exact zero',state_radius_chain_rule='exact zero',
                direct_metric_tensor_components='exact zero',angular_equation='exact zero',mixed_system_elimination='exact zero',
                absolute_J_LLL_zero_only_in_stated_2D_representative=True)


def smooth_step(x):
    if x<=0:return mp.mpf(0)
    if x>=1:return mp.mpf(1)
    return 10*x**3-15*x**4+6*x**5


def step_integral(x):
    if x<=0:return mp.mpf(0)
    if x>=1:return x-mp.mpf('.5')
    return mp.mpf('2.5')*x**4-3*x**5+x**6


def theta(t):
    eta=mp.mpf(1)/2000
    tau=mp.mpf(1)/4
    return 4*eta*tau*(step_integral(t/tau)-step_integral((t-mp.pi)/tau))


def rate(t):
    eta=mp.mpf(1)/2000
    tau=mp.mpf(1)/4
    return 4*eta*(smooth_step(t/tau)-smooth_step((t-mp.pi)/tau))


class RetardedExample:
    """A C3 optical bump creates an actual timelike paired-domain history.

    Its actual lab/global clock is the frozen reference clock, hence the
    fixed laboratory source at event t is exactly Theta0(t).
    """
    def hjets(self,v):
        end=mp.mpf(8)
        if v<=0 or v>=end:return [v,mp.mpf(1),mp.mpf(0),mp.mpf(0)]
        delta=mp.mpf('.1')
        q=v*(end-v)/end**2
        qp=(end-2*v)/end**2
        qpp=-2/end**2
        return [v+delta*q**4,1+4*delta*q**3*qp,
                delta*(12*q*q*qp*qp+4*q**3*qpp),
                delta*(24*q*qp**3+36*q*q*qp*qpp)]

    def inverse(self,w):
        if w<=0 or w>=8:return w
        return mp.findroot(lambda v:self.hjets(v)[0]-w,(w-mp.mpf('.001'),w+mp.mpf('.001')))

    def iterated(self,v,j):
        h,hp,hpp,hppp=self.hjets(v)
        w=v if j==0 else self.inverse(h-j*mp.pi)
        _,wp,wpp,wppp=self.hjets(w)
        d1=hp/wp
        d2=(hpp-wpp*d1*d1)/wp
        d3=(hppp-wppp*d1**3-3*wpp*d1*d2)/wp
        return [w,d1,d2,d3]

    def scattering(self,v,j):
        cur=self.iterated(v,j)
        past=self.iterated(v,j+1)
        return [(cur[0]+past[0]+mp.pi)/2]+[(cur[k]+past[k])/2 for k in range(1,4)]

    def phase(self,v):
        value=mp.mpf(0)
        derivative=mp.mpf(0)
        count=0
        for j in range(20):
            event=self.scattering(v,j)
            if event[0]<=0:break
            value+=theta(event[0])
            derivative+=rate(event[0])*event[1]
            count+=1
        else:raise AssertionError('causal source sum exceeded frozen twenty-return cap')
        return value,derivative,count

    def stress(self,v,N=1):
        _,hp,hpp,hppp=self.hjets(v)
        _,ap,_=self.phase(v)
        sch=hppp/hp-mp.mpf('1.5')*(hpp/hp)**2
        return N*(-hp*hp/(96*mp.pi)-sch/(48*mp.pi)+(ap*ap-hp*ap)/(8*mp.pi))

    def kernel(self,v1,v2):
        h1,hp1,_,_=self.hjets(v1)
        h2,hp2,_,_=self.hjets(v2)
        a1,_,_=self.phase(v1)
        a2,_,_=self.phase(v2)
        return mp.sqrt(hp1*hp2)*rotation(a1)*vacuum(h1-h2)*rotation(a2).T


def rotation(a):
    return mp.matrix([[mp.cos(a),-mp.sin(a)],[mp.sin(a),mp.cos(a)]])


def vacuum(d):
    J=mp.matrix([[0,-1],[1,0]])
    return mp.eye(2)/(4*mp.pi*1j*mp.sin(d/2))+1j*J/(4*mp.pi*mp.cos(d/2))


def schwarz(jets):
    return jets[3]/jets[1]-mp.mpf('1.5')*(jets[2]/jets[1])**2


def max_matrix_error(matrix):
    return max(abs(matrix[i,j]) for i in range(matrix.rows) for j in range(matrix.cols))


def causal_and_point_split(example):
    J=mp.matrix([[0,-1],[1,0]])
    close(max_matrix_error(vacuum(mp.mpf('.37')+mp.pi)-J*vacuum(mp.mpf('.37'))))
    rows=[]
    source_timing_errors=[]
    for rawv in ['.1','.4','1.3','3.5','5','6.9']:
        v=mp.mpf(rawv)
        h=example.hjets(v)
        b=example.iterated(v,1)
        hb=example.hjets(b[0])
        event=example.scattering(v,0)
        av,ap,count=example.phase(v)
        ab,apb,_=example.phase(b[0])
        errors=[abs(hb[0]-h[0]+mp.pi),abs(h[1]-hb[1]*b[1]),
                abs(schwarz(h)-b[1]**2*schwarz(hb)-schwarz(b)),
                abs(av-theta(event[0])-ab),abs(ap-rate(event[0])*event[1]-apb*b[1])]
        for error in errors:close(error)
        st=(b[1]-1)/(b[1]+1)
        stt=4*b[2]/(b[1]+1)**3
        sttt=8*b[3]/(b[1]+1)**4-24*b[2]**2/(b[1]+1)**5
        sch_boundary=(2*sttt/(1-st*st)+6*st*stt*stt/(1-st*st)**2)/(1-st)**2
        close(schwarz(b)-sch_boundary)
        assert abs(st)<1
        second=v+mp.mpf('.37')
        original=example.kernel(v,second)
        transported=mp.sqrt(b[1])*rotation(theta(event[0]))*J*example.kernel(b[0],second)
        kernel_error=max_matrix_error(original-transported)
        close(kernel_error)
        splits=[]
        for rawd in SEPARATIONS:
            d=mp.mpf(rawd)
            measured=-1j*(example.kernel(v+d/2,v-d/2)[0,0]-1/(2*mp.pi*1j*d))/(2*d)
            error=abs(measured-example.stress(v))
            splits.append(dict(separation=rawd,stress=str(measured),error=str(error)))
        assert mp.mpf(splits[-1]['error'])<mp.mpf(TOL['finest_point_split'])
        naive=sum(theta(v-j*mp.pi) for j in range(4))
        source_timing_errors.append(abs(av-naive))
        rows.append(dict(v=rawv,returns=count,endpoint_velocity=str(st),source_scattering_time=str(event[0]),
                         canonical_stress=str(example.stress(v)),causal_errors=list(map(str,errors)),
                         holonomy_kernel_error=str(kernel_error),boundary_Schwarzian_error=str(abs(schwarz(b)-sch_boundary)),point_splits=splits))
    return dict(cases=rows,maximum_actual_vs_naive_phase_difference=str(max(source_timing_errors)),
                actual_lab_clock_is_frozen_reference=True,physical_4D_geometry_solution_asserted=False)


def coordinate_covariance(example):
    N=mp.mpf(13)
    f=lambda x:x+mp.mpf('.02')*x**3
    g=lambda x:x+mp.mpf('.015')*x**3
    rho=lambda xp,xm:-mp.log(mp.sin((xp-xm)/2))+mp.mpf('.03')*mp.sin(mp.mpf('.4')*xp)*mp.cos(mp.mpf('.7')*xm)
    rows=[]
    for rawp,rawm in [('3','1.2'),('3.4','.8')]:
        p,m=mp.mpf(rawp),mp.mpf(rawm)
        xp,xm=f(p),g(m)
        tilde=lambda pp,mm:rho(f(pp),g(mm))+mp.log(mp.diff(f,pp)*mp.diff(g,mm))/2
        values=[]
        for orientation,fn,xi,x,shift in [(0,f,p,xp,mp.pi),(1,g,m,xm,mp.mpf(0))]:
            fj=[fn(xi),mp.diff(fn,xi),mp.diff(fn,xi,2),mp.diff(fn,xi,3)]
            oldh=example.hjets(x-shift)
            newh=[oldh[0],oldh[1]*fj[1],oldh[2]*fj[1]**2+oldh[1]*fj[2],
                  oldh[3]*fj[1]**3+3*oldh[2]*fj[1]*fj[2]+oldh[1]*fj[3]]
            ap=example.phase(x-shift)[1]
            apnew=ap*fj[1]
            oldt=example.stress(x-shift,N)
            newt=N*(-newh[1]**2/(96*mp.pi)-schwarz(newh)/(48*mp.pi)+(apnew*apnew-newh[1]*apnew)/(8*mp.pi))
            orders1=(1,0) if orientation==0 else (0,1)
            orders2=(2,0) if orientation==0 else (0,2)
            oldq=-N*(mp.diff(rho,(xp,xm),orders1)**2-mp.diff(rho,(xp,xm),orders2))/(24*mp.pi)
            newq=-N*(mp.diff(tilde,(p,m),orders1)**2-mp.diff(tilde,(p,m),orders2))/(24*mp.pi)
            canonical_error=abs(newt-fj[1]**2*oldt+N*schwarz(fj)/(48*mp.pi))
            Weyl_error=abs(newq-fj[1]**2*oldq-N*schwarz(fj)/(48*mp.pi))
            total_error=abs(newt+newq-fj[1]**2*(oldt+oldq))
            for error in [canonical_error,Weyl_error,total_error]:close(error)
            omitted=abs(newt+fj[1]**2*oldq-fj[1]**2*(oldt+oldq))
            assert omitted>mp.mpf(TOL['missing_cocycle_discrepancy_min'])
            values.append(dict(orientation=orientation,canonical_error=str(canonical_error),Weyl_error=str(Weyl_error),
                               total_tensor_error=str(total_error),missing_Weyl_cocycle_error=str(omitted)))
        oldmixed=-N*mp.diff(rho,(xp,xm),(1,1))/(24*mp.pi)
        newmixed=-N*mp.diff(tilde,(p,m),(1,1))/(24*mp.pi)
        close(newmixed-mp.diff(f,p)*mp.diff(g,m)*oldmixed)
        rows.append(dict(xi_plus=rawp,xi_minus=rawm,components=values,mixed_error=str(abs(newmixed-mp.diff(f,p)*mp.diff(g,m)*oldmixed))))
    return rows


def signed_Weyl_charge():
    # Symbolic local identity; the endpoint term is retained.
    s,t=sp.symbols('s t',real=True)
    w=sp.Function('w')(t,s)
    local=-(sp.diff(w,t)**2+sp.diff(w,s)**2-2*sp.cot(s)*sp.diff(w,s)-2*sp.diff(w,s,2))/(48*sp.pi)
    divergence=sp.diff(sp.diff(w,s)+sp.cot(s)*w,s)/(24*sp.pi)
    bulk=-(sp.diff(w,t)**2+sp.diff(w,s)**2-2*w/sp.sin(s)**2)/(48*sp.pi)
    assert sp.trigsimp(local-divergence-bulk)==0
    eta=mp.mpf(1)/2000
    mu=mp.mpf(4)/mp.mpf(2000)**10
    K=mp.mpf('1e5')
    e=mp.mpf('1e-5')
    eps=mp.atan(e)
    I2=2*mp.cot(eps)
    I4=2*mp.cot(eps)+mp.mpf(2)/3*mp.cot(eps)**3
    I3=mp.csc(eps)*mp.cot(eps)-mp.log(mp.tan(eps/2))
    close((mp.quad(lambda x:mp.csc(x)**3,[eps,mp.pi/2,mp.pi-eps])-I3)/I3)
    bound=2*mu*(K+K)*mp.csc(eps)**2/(24*mp.pi*eta)+(mu*mu*(K*K*I2+K*K*I4)+2*mu*K*I3)/(48*mp.pi*eta)
    assert bound<mp.mpf('5.2e-16') and bound/2<mp.mpf('2.6e-16')
    rows=[]
    for rawt in ['0','.7','2']:
        tt=mp.mpf(rawt)
        amp=mu*K
        w=lambda x:amp*(mp.mpf('.2')*mp.cos(mp.mpf('.7')*tt)*mp.csc(x)+mp.mpf('.1')*mp.sin(x))
        wt=lambda x:-amp*mp.mpf('.14')*mp.sin(mp.mpf('.7')*tt)*mp.csc(x)
        ws=lambda x:amp*(-mp.mpf('.2')*mp.cos(mp.mpf('.7')*tt)*mp.cos(x)*mp.csc(x)**2+mp.mpf('.1')*mp.cos(x))
        wss=lambda x:amp*(mp.mpf('.2')*mp.cos(mp.mpf('.7')*tt)*(mp.csc(x)+2*mp.cos(x)**2*mp.csc(x)**3)-mp.mpf('.1')*mp.sin(x))
        direct=mp.quad(lambda x:-(wt(x)**2+ws(x)**2-2*mp.cot(x)*ws(x)-2*wss(x))/(48*mp.pi),[eps,mp.pi/2,mp.pi-eps])
        endpoint=((ws(mp.pi-eps)+mp.cot(mp.pi-eps)*w(mp.pi-eps))-(ws(eps)+mp.cot(eps)*w(eps)))/(24*mp.pi)
        remainder=mp.quad(lambda x:-(wt(x)**2+ws(x)**2-2*mp.csc(x)**2*w(x))/(48*mp.pi),[eps,mp.pi/2,mp.pi-eps])
        close(direct-endpoint-remainder)
        normalized=abs(direct)/eta
        assert normalized<=bound
        rows.append(dict(t=rawt,signed_charge_per_N=str(direct),endpoint_per_N=str(endpoint),
                         identity_error=str(abs(direct-endpoint-remainder)),normalized_signed_charge=str(normalized)))
    kappa_ratio=mu/(24*mp.pi*eta*mp.mpf('.99'))
    assert kappa_ratio<mp.mpf('1.05e-31')
    return dict(symbolic_signed_identity='exact zero',mu_upper=str(mu),epsilon=str(eps),conditional_W_constant=str(K),
                signed_charge_bound=str(bound),clock_share_bound=str(bound/2),kappa_over_X_conditional_upper=str(kappa_ratio),
                absolute_spatial_stress_bound_asserted=False,candidate_weighted_jets_established=False,cases=rows)


def ray_error_bound_checks(example):
    # Direct finite comparison to the fixed unit-pi reference, with the
    # same actual clock/source, so no source protocol is changed.
    rows=[]
    for rawv in ['.4','1.3','3.5','5','6.9']:
        v=mp.mpf(rawv)
        h=example.hjets(v)
        ap=example.phase(v)[1]
        d0=sum(rate(v-j*mp.pi) for j in range(4))
        M,M0,H1,H0=abs(ap),abs(d0),abs(h[1]),mp.mpf(1)
        ea,eh,es=abs(ap-d0),abs(h[1]-1),abs(schwarz(h))
        bound=((M+M0+H0)*ea+M*eh)/(8*mp.pi)+(H1+H0)*eh/(96*mp.pi)+es/(48*mp.pi)
        reference=-1/(96*mp.pi)+(d0*d0-d0)/(8*mp.pi)
        measured=abs(example.stress(v)-reference)
        assert measured<=bound+mp.mpf(TOL['exact_numerical'])
        rows.append(dict(v=rawv,eps_a=str(ea),eps_h=str(eh),eps_S=str(es),chiral_error=str(measured),bound=str(bound)))
    close(mp.pi*(mp.mpf('2.4e-8')/(48*mp.pi))/(mp.mpf(1)/2000)-mp.mpf('1e-6'))
    return rows


def main():
    mp.mp.dps=DPS
    hashes={name:sha(ROOT/name) for name in SOURCES}
    result=dict(story='S137.2',started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                dps=DPS,tolerances=TOL,point_split_separations=SEPARATIONS,input_sha256=hashes,
                hash_policy='SHA-256 after CRLF-to-LF normalization; no other byte changes',
                execution=dict(executable=sys.executable,argv=sys.orig_argv,cwd=str(Path.cwd()),command=subprocess.list2cmdline(sys.orig_argv)),
                runtime=dict(python=sys.version,platform=platform.platform(),sympy=sp.__version__,mpmath=mp.__version__),
                primary_sources=[dict(url='https://arxiv.org/html/1807.04726v3',scope='appendices C and F: exact canonical reduction and 2D Weyl completion'),
                                 dict(url='https://arxiv.org/html/1804.00491v3',scope='appendix C.211: known constant AdS anomaly/dilaton shift')])
    if OUTPUT.exists():
        previous=json.loads(OUTPUT.read_text(encoding='utf-8'))
        if previous.get('status')=='failed':record('Retained prior failed receipt: '+json.dumps(previous))
    record('S137.2 evolving LLL state/geometry check: '+result['started_utc'])
    record('Execution: '+json.dumps(result['execution']))
    record('Frozen acceptance: '+json.dumps(dict(dps=DPS,tolerances=TOL,separations=SEPARATIONS)))
    try:
        result['symbolic_point_split_and_cocycles']=symbolic_point_split_and_covariance()
        record('Symbolic point splitting, coordinate cocycles and canonical radius cancellation: PASS')
        result['symbolic_Ward_and_bulk_system']=symbolic_ward_and_system()
        record('Direct metric tensor, Ward/trace and sourced bulk system: PASS')
        example=RetardedExample()
        result['retarded_source_and_point_splitting']=causal_and_point_split(example)
        record('Actual shaped retarded source, moving rays, holonomy and independent point splitting: PASS')
        result['coordinate_covariance']=coordinate_covariance(example)
        record('Complete stress covariance and missing-cocycle nonzero control: PASS')
        result['signed_Weyl_charge']=signed_Weyl_charge()
        record('Independent signed Weyl-charge quadratures and conditional candidate bound: PASS')
        result['ray_stress_bound']=ray_error_bound_checks(example)
        record('Finite ray-stress error bounds in the fixed comparison chart: PASS')
        assert hashes=={name:sha(ROOT/name) for name in SOURCES},'input changed during execution'
        result['status']='pass_exact_evolving_LLL_state_and_scoped_geometry_system'
        result['claims']=dict(finite_free_LLL_bulk_state_radius_source_zero=True,full_four_dimensional_radius_source_zero=False,
                              physical_characteristic_map_established=False,known_Weyl_signed_charge_bound_conditional=True,
                              absolute_spatial_Weyl_stress_bound_established=False,physical_boundary_defect_bounded=False,
                              complete_retarded_RN4_geometry_solved=False,physical_residual_certified=False,
                              original_preparation_rectangle_certified=False,source_protocol_changed=False,
                              full_4D_anomaly_added_twice=False,physical_opening_achieved=False,s074_gate_closed=False,lean_invoked=False)
        splits=[row['point_splits'][-1]['error'] for row in result['retarded_source_and_point_splitting']['cases']]
        result['summary']=dict(retarded_cases=6,independent_point_splits=18,coordinate_covariance_points=2,
                               signed_charge_cases=3,maximum_finest_point_split_error=str(max(map(mp.mpf,splits))),
                               maximum_tensor_covariance_error=str(max(mp.mpf(comp['total_tensor_error']) for row in result['coordinate_covariance'] for comp in row['components'])),
                               conditional_known_Weyl_clock_share=result['signed_Weyl_charge']['clock_share_bound'])
        write(OUTPUT,json.dumps(result,indent=2)+'\n')
        record('Summary: '+json.dumps(result['summary']))
        record('Normalized input hashes: '+json.dumps(hashes))
        record('Terminal status: PASS, exit 0')
    except Exception:
        result['status']='failed'
        result['failure_traceback']=traceback.format_exc()
        write(OUTPUT,json.dumps(result,indent=2)+'\n')
        record(result['failure_traceback'].rstrip())
        record('Failure input hashes: '+json.dumps(hashes))
        record('Terminal status: FAIL, exit 1; retained in append-only log')
        raise


if __name__=='__main__':
    main()
