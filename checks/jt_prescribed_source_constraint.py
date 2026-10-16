"""S135 minimal JT metric constraint and midpoint source work, no Lean."""
import datetime
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import traceback

import mpmath as mp
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'evidence/jt-prescribed-source-constraint.json'
LOG=ROOT/'evidence/jt-prescribed-source-constraint.log'
DPS=60
STEPS=['2e-5','1e-5']
TOL=dict(exact_numeric='1e-35',laboratory_energy_difference='5e-8',midpoint_covariance_limit='5e-10')
SOURCES=['notes/jt-prescribed-source-constraint.md','checks/jt_prescribed_source_constraint.py',
         'notes/prescribed-source-clock-variation.md','notes/switching-radiation-capture.md',
         'notes/returned-switch-capture.md','notes/driven-bulk-geometry.md',
         'notes/boundary-vertex-matching.md']


def sha(path):
    return hashlib.sha256(path.read_bytes().replace(b'\r\n',b'\n')).hexdigest()


def write(path,value):
    with path.open('w',encoding='utf-8',newline='\n') as stream:
        stream.write(value)


def record(value):
    with LOG.open('a',encoding='utf-8',newline='\n') as stream:
        stream.write(value+'\n')
    print(value,flush=True)


def zero(value):
    if s.simplify(s.trigsimp(value)) != 0:
        raise AssertionError(str(value))


def symbolic_geometry():
    t,x=s.symbols('t sigma',real=True)
    phi=s.Function('Phi')(t,x)
    g=s.diag(-s.sin(x)**-2,s.sin(x)**-2)
    inv=g.inv()
    coordinates=(t,x)
    christ=[[[s.simplify(sum(inv[a,d]*(s.diff(g[d,c],coordinates[b])+s.diff(g[d,b],coordinates[c])-s.diff(g[b,c],coordinates[d]))/2 for d in range(2))) for c in range(2)] for b in range(2)] for a in range(2)]
    hessian=s.Matrix(2,2,lambda a,b:s.diff(phi,coordinates[a],coordinates[b])-sum(christ[c][a][b]*s.diff(phi,coordinates[c]) for c in range(2)))
    box=sum(inv[a,b]*hessian[a,b] for a in range(2) for b in range(2))
    tensor=g*box-hessian-g*phi
    expected=s.Matrix([[-s.diff(phi,x,2)-s.cot(x)*s.diff(phi,x)+phi/s.sin(x)**2,
                        -s.diff(phi,t,x)-s.cot(x)*s.diff(phi,t)],
                       [-s.diff(phi,t,x)-s.cot(x)*s.diff(phi,t),
                        -s.diff(phi,t,2)-s.cot(x)*s.diff(phi,x)-phi/s.sin(x)**2]])
    for a in range(2):
        for b in range(2):
            zero(tensor[a,b]-expected[a,b])
    zero(tensor[0,0]+s.diff(s.diff(phi,x)+s.cot(x)*phi,x))
    trace=sum(inv[a,b]*tensor[a,b] for a in range(2) for b in range(2))
    zero(trace-(s.sin(x)**2*(-s.diff(phi,t,2)+s.diff(phi,x,2))-2*phi))

    BL,AL,CL,DL=[s.Function(name)(t) for name in ['B_L','A_L','C_L','D_L']]
    expansion=BL/x+AL+CL*x+DL*x*x
    trace_series=-s.diff(expansion,t,2)+s.diff(expansion,x,2)-2*expansion/s.sin(x)**2
    zero(s.limit(x*x*trace_series,x,0)+2*AL)
    zero(s.limit(x*trace_series.subs(AL,0),x,0)+s.diff(BL,t,2)+2*CL+2*BL/3)
    Cvalue=-s.diff(BL,t,2)/2-BL/3
    left=s.limit((s.diff(expansion,x)+s.cot(x)*expansion).subs({AL:0,CL:Cvalue}),x,0)
    zero(left+s.diff(BL,t,2)+BL)

    BR,AR,CR,DR=[s.Function(name)(t) for name in ['B_R','A_R','C_R','D_R']]
    right_exp=BR/(s.pi-x)+AR+CR*(s.pi-x)+DR*(s.pi-x)**2
    v=s.symbols('v',positive=True)
    right_trace=(-s.diff(right_exp,t,2)+s.diff(right_exp,x,2)-2*right_exp/s.sin(x)**2).subs(x,s.pi-v)
    zero(s.limit(v*v*right_trace,v,0)+2*AR)
    zero(s.limit(v*right_trace.subs(AR,0),v,0)+s.diff(BR,t,2)+2*CR+2*BR/3)
    right=s.limit((s.diff(right_exp,x)+s.cot(x)*right_exp).subs({AR:0,CR:-s.diff(BR,t,2)/2-BR/3}).subs(x,s.pi-v),v,0)
    zero(right-s.diff(BR,t,2)-BR)
    return dict(metric='csc(sigma)^2 diag(-1,1)',tensor_components=str(expected),
                left_boundary_limit=str(left),right_boundary_limit=str(right),
                energy_identity='E=-(B_Ltt+B_L+B_Rtt+B_R)',
                boundary_normalization='sigma_cut=epsilon*y+O(epsilon^3), Phi_cut=N/epsilon => B_L=B_R=N*y')


def symbolic_charges_and_work():
    u=s.symbols('u',real=True)
    N=s.symbols('N',positive=True)
    clock=s.Function('t')(u)
    y=s.diff(clock,u)
    lagrangian=N*(s.diff(clock,u,2)/y)**2-N*y*y
    p1=s.diff(lagrangian,s.diff(clock,u,2))
    p0=s.diff(lagrangian,y)-s.diff(p1,u)
    logy=s.log(y)
    expected_charge=-2*N*(y+s.diff(logy,u,2)/y)
    zero(p0-expected_charge)
    hamiltonian=p0*y+p1*s.diff(clock,u,2)-lagrangian
    schwarz=s.diff(clock,u,3)/y-s.Rational(3,2)*(s.diff(clock,u,2)/y)**2+y*y/2
    zero(hamiltonian+2*N*schwarz)
    velocity,pu,puu,E,Et=s.symbols('y phi_u phi_uu E E_t',real=True)
    H=N*(pu**2-2*puu-velocity**2)
    constraint=-velocity**2-E*velocity/(2*N)
    zero(H.subs(puu,constraint)-(N*(pu**2+velocity**2)+velocity*E))
    derivative=2*N*pu*puu+2*N*velocity**2*pu+velocity*pu*E+velocity**2*Et
    zero(derivative.subs(puu,constraint)-velocity**2*Et)

    gamma=s.symbols('gamma',real=True)
    J=s.Matrix([[0,-1],[1,0]])
    I=s.eye(2)
    R=(I-gamma*J/2).inv()*(I+gamma*J/2)
    for entry in R.T*R-I:
        zero(entry)
    for entry in -s.I*(R-I)+s.I*gamma*J*(I+R)/2:
        zero(entry)
    theta=s.symbols('theta',real=True)
    vac=N*(-theta/(4*s.pi)+theta**2/(4*s.pi**2))
    sourcegamma=2*s.tan(theta/2)
    static_X=N*s.cos(theta/2)**2*(1-2*theta/s.pi)/(4*s.pi)
    zero(s.diff(vac,theta)/s.diff(sourcegamma,theta)+static_X)
    r,Dret,y0=s.symbols('r Dret y',real=True)
    X=N*s.cos(theta/2)**2*(1-r-2*Dret)/(4*s.pi)
    gamma_lab_rate=y0*r/s.cos(theta/2)**2
    bulk_power=N*r*(r+2*Dret-1)/(4*s.pi)
    zero(-gamma_lab_rate*y0*X-y0*y0*bulk_power)
    O,Iflux,anom=s.symbols('O I anomaly',real=True)
    zero((y0*y0*O+anom)-(y0*y0*Iflux+anom)-y0*y0*(O-Iflux))
    return dict(global_gravity_charge=str(expected_charge),laboratory_hamiltonian=str(hamiltonian),
                constrained_energy='H_u=N(phi_u^2+y^2)+y E',laboratory_work='dH_u/du=y^2 E_t',
                Cayley_jump=str(R),static_Hellmann_Feynman='dE_vac/dgamma=-X_static',
                source_work='-gamma_b_prime*X_u=y^2 E_t',normal_order_flux_anomaly='cancels between both chiral directions')


def numerical_dilatons():
    tolerance=mp.mpf(TOL['exact_numeric'])
    rows=[]
    amplitude=mp.mpf('0.005')
    H=lambda v: amplitude*v**6 if v>0 else mp.mpf(0)
    Hp=lambda v: 6*amplitude*v**5 if v>0 else mp.mpf(0)
    Hpp=lambda v: 30*amplitude*v**4 if v>0 else mp.mpf(0)
    p=lambda v: -amplitude*(120*v**3+6*v**5) if v>0 else mp.mpf(0)
    hom=lambda t: mp.mpf('0.3')*mp.cos(t)+mp.mpf('0.1')*mp.sin(t)
    for name in ['homogeneous','static','retarded_polynomial']:
        def phi(t,x):
            value=hom(t)/mp.sin(x)
            if name=='static':
                value+=mp.mpf('0.04')*(1+(mp.pi/2-x)*mp.cot(x))
            if name=='retarded_polynomial':
                value+=Hp(t-x)+Hp(t+x-mp.pi)+mp.cot(x)*(H(t-x)-H(t+x-mp.pi))
            return value
        cases=[]
        for t in [mp.mpf('0.17'),mp.mpf('1.13'),mp.pi+mp.mpf('0.37')]:
            largest=mp.mpf(0)
            for x in [mp.mpf('0.23'),mp.mpf('0.91'),mp.mpf('1.31'),mp.pi-mp.mpf('0.31')]:
                dt=mp.diff(lambda tt:phi(tt,x),t)
                ds=mp.diff(lambda xx:phi(t,xx),x)
                dtt=mp.diff(lambda tt:phi(tt,x),t,2)
                dss=mp.diff(lambda xx:phi(t,xx),x,2)
                dts=mp.diff(lambda xx:mp.diff(lambda tt:phi(tt,xx),t),x)
                actual=[-dss-mp.cot(x)*ds+phi(t,x)/mp.sin(x)**2,
                        -dtt-mp.cot(x)*ds-phi(t,x)/mp.sin(x)**2,
                        -dts-mp.cot(x)*dt]
                if name=='static':
                    expected=[-mp.mpf('0.04'),-mp.mpf('0.04'),mp.mpf(0)]
                elif name=='retarded_polynomial':
                    expected=[p(t-x)+p(t+x-mp.pi),p(t-x)+p(t+x-mp.pi),p(t+x-mp.pi)-p(t-x)]
                else:
                    expected=[mp.mpf(0)]*3
                error=max(abs(a-b) for a,b in zip(actual,expected))
                error=max(error,abs(-dtt+dss-2*phi(t,x)/mp.sin(x)**2))
                assert error<tolerance,(name,t,x,error)
                largest=max(largest,error)
            if name=='retarded_polynomial':
                cuts=sorted(set([mp.mpf(0),mp.pi]+[v for v in [t,mp.pi-t] if 0<v<mp.pi]))
                energy=mp.quad(lambda x:p(t-x)+p(t+x-mp.pi),cuts)
                B=lambda tt: hom(tt)+H(tt)-H(tt-mp.pi)
                independent=-2*(Hpp(t)-Hpp(t-mp.pi)+H(t)-H(t-mp.pi))
            elif name=='static':
                energy=mp.quad(lambda x:-mp.mpf('0.04'),[0,mp.pi])
                B=lambda tt:hom(tt)+mp.pi*mp.mpf('0.04')/2
                independent=-mp.pi*mp.mpf('0.04')
            else:
                energy=mp.mpf(0)
                B=hom
                independent=mp.mpf(0)
            charge=mp.diff(B,t,2)+B(t)
            assert abs(energy-independent)<tolerance
            assert abs(charge+energy/2)<tolerance
            cases.append(dict(t=str(t),maximum_tensor_error=str(largest),energy=str(energy),
                              boundary_charge_error=str(abs(charge+energy/2))))
        rows.append(dict(family=name,cases=cases,physical_quantum_state_assigned=False))
    return rows


def numerical_laboratory_energy():
    N=mp.mpf(11)
    rows=[]
    for name in ['polynomial','sinusoidal']:
        clock=(lambda u:u+mp.mpf('0.05')*u**3) if name=='polynomial' else (lambda u:mp.mpf('1.1')*u+mp.mpf('0.07')*mp.sin(2*u))
        y=lambda u:mp.diff(clock,u)
        phi=lambda u:mp.log(y(u))
        E=lambda u:-2*N*(y(u)+mp.diff(phi,u,2)/y(u))
        def energy(u):
            q=lambda v:mp.tan(clock(v)/2)
            return -2*N*(mp.diff(q,u,3)/mp.diff(q,u)-mp.mpf('1.5')*(mp.diff(q,u,2)/mp.diff(q,u))**2)
        cases=[]
        for u in [mp.mpf('0.1'),mp.mpf('0.5'),mp.mpf('0.9')]:
            identity=energy(u)-(N*(mp.diff(phi,u)**2+y(u)**2)+y(u)*E(u))
            exact_power=y(u)*mp.diff(E,u)
            exact_error=abs(mp.diff(energy,u)-exact_power)
            assert abs(identity)<mp.mpf(TOL['exact_numeric'])
            assert exact_error<mp.mpf(TOL['exact_numeric'])
            differences=[]
            for raw_step in STEPS:
                step=mp.mpf(raw_step)
                error=abs((energy(u+step)-energy(u-step))/(2*step)-exact_power)
                assert error<mp.mpf(TOL['laboratory_energy_difference']),error
                differences.append(dict(step=raw_step,error=str(error)))
            cases.append(dict(u=str(u),y=str(y(u)),E=str(E(u)),laboratory_energy=str(energy(u)),
                              algebra_error=str(abs(identity)),exact_work_error=str(exact_error),differences=differences))
        rows.append(dict(clock_family=name,cases=cases,quantum_stress_realization_asserted=False))
    return rows


def midpoint_work():
    J=mp.matrix([[0,-1],[1,0]])
    I=mp.eye(2)
    R=lambda angle:mp.cos(angle)*I+mp.sin(angle)*J
    theta=lambda t:mp.mpf('0.1')+mp.mpf('0.07')*t+mp.mpf('0.03')*mp.sin(t)
    past=lambda t:mp.mpf('0.03')*t+mp.mpf('0.02')*mp.sin(mp.mpf('0.6')*t)
    rows=[]
    for N in [1,11]:
        for t in [mp.mpf('0.1'),mp.mpf('0.7'),mp.mpf('1.4')]:
            r=mp.diff(theta,t)
            Dret=mp.diff(past,t)
            expected=N*mp.cos(theta(t)/2)**2*(1-r-2*Dret)/(4*mp.pi)
            limit_checks=[]
            for raw_sep in STEPS:
                d=mp.mpf(raw_sep)
                t1,t2=t+d/2,t-d/2
                initial=I/(4*mp.pi*1j*mp.sin(d/2))+1j*J/(4*mp.pi*mp.cos(d/2))
                incoming=R(past(t1))*initial*R(past(t2)).T
                M1,M2=(I+R(theta(t1)))/2,(I+R(theta(t2)))/2
                midpoint=M1*incoming*M2.T
                measured=N*1j*midpoint[0,1]
                error=abs(measured-expected)
                assert error<mp.mpf(TOL['midpoint_covariance_limit']),error
                limit_checks.append(dict(separation=raw_sep,measured_X=str(measured),error=str(error)))
            # Independently compare the midpoint contact work with the
            # incoming/outgoing stress difference, including the return term.
            stress=lambda rate:N*(rate*rate-rate)/(8*mp.pi)
            global_power=2*(stress(Dret+r)-stress(Dret))
            y=mp.mpf('1.1')+mp.mpf('0.2')*mp.cos(t)
            gamma_b_prime=y*r/mp.cos(theta(t)/2)**2
            laboratory_source_power=-gamma_b_prime*y*expected
            error=abs(laboratory_source_power-y*y*global_power)
            assert error<mp.mpf(TOL['exact_numeric'])
            rows.append(dict(N=N,t=str(t),theta=str(theta(t)),source_rate=str(r),return_rate=str(Dret),
                             X_exact=str(expected),midpoint_limits=limit_checks,y=str(y),
                             laboratory_power=str(laboratory_source_power),work_error=str(error),
                             synthetic_local_covariance_not_full_return_solution=True))
    return rows


def main():
    mp.mp.dps=DPS
    hashes={name:sha(ROOT/name) for name in SOURCES}
    result=dict(story='S135.2',started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                dps=DPS,tolerances=TOL,difference_steps=STEPS,input_sha256=hashes,
                hash_policy='SHA-256 after CRLF-to-LF normalization; no other byte changes',
                execution=dict(executable=sys.executable,argv=sys.orig_argv,cwd=str(Path.cwd()),command=subprocess.list2cmdline(sys.orig_argv)),
                runtime=dict(python=sys.version,platform=platform.platform(),sympy=s.__version__,mpmath=mp.__version__),
                primary_source=dict(url='https://arxiv.org/html/1804.00491v3',scope='Section 2.2 action/boundary data; appendix B Schwarzian charges; appendix C.2-C.3 static free comparison',dynamic_source_constraint_imported=False))
    if OUTPUT.exists():
        previous=json.loads(OUTPUT.read_text(encoding='utf-8'))
        if previous.get('status')=='failed':
            record('Retained prior failed receipt: '+json.dumps(previous))
    record('S135.2 minimal JT constraint run: '+result['started_utc'])
    record('Execution: '+json.dumps(result['execution']))
    record('Predeclared checks: '+json.dumps(dict(dps=DPS,tolerances=TOL,steps=STEPS)))
    try:
        result['metric_and_boundary_series']=symbolic_geometry()
        record('Metric tensor and both boundary series: PASS')
        result['charges_and_work']=symbolic_charges_and_work()
        record('Schwarzian charges, total energy, Cayley jump and work: PASS')
        result['dilaton_checks']=numerical_dilatons()
        record('Homogeneous, static and retarded source dilatons: PASS')
        result['laboratory_energy_checks']=numerical_laboratory_energy()
        record('Independent laboratory-energy derivatives: PASS')
        result['midpoint_source_work']=midpoint_work()
        record('Midpoint point splitting, source work and returned stress: PASS')
        assert hashes=={name:sha(ROOT/name) for name in SOURCES},'input changed during execution'
        result['status']='pass_minimal_JT_constraint_and_prescribed_source_work'
        result['claims']=dict(minimal_JT_symmetric_constraint_derived=True,retimed_source_required=True,
                             independently_added_contact_energy=False,full_RN4_clock_derived=False,
                             physical_source_controller_implemented=False,interacting_state_derived=False,
                             s074_gate_closed=False,physical_opening_achieved=False,lean_invoked=False)
        result['summary']=dict(dilaton_families=3,dilaton_tensor_points=36,laboratory_clock_families=2,
                               midpoint_cases=6,midpoint_limits=12,
                               maximum_lab_difference_error=str(max(mp.mpf(d['error']) for row in result['laboratory_energy_checks'] for case in row['cases'] for d in case['differences'])),
                               maximum_midpoint_limit_error=str(max(mp.mpf(d['error']) for case in result['midpoint_source_work'] for d in case['midpoint_limits'])))
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
