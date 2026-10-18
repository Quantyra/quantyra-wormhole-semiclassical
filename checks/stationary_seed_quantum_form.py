"""Finite algebra checks for seed contacts, LLL and the flat quantum model.

No actual curved quantum lower bound, kernel or physical matching is computed.
Analytic source, angular convergence and closure arguments remain ordinary.
"""
from __future__ import annotations
import base64,contextlib,hashlib,json,platform,sys,traceback
from datetime import datetime,timezone
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/stationary-seed-quantum-form.json'
LOG=ROOT/'evidence/stationary-seed-quantum-form.log'
INPUTS=[
    'checks/stationary_seed_quantum_form.py',
    'notes/stationary-seed-quantum-form.md',
    'notes/stationary-global-tt-reduction.md',
    'notes/stationary-centered-exterior-source.md',
    'notes/stationary-horizon-source-defect.md',
    'notes/stationary-centered-covariant-seed.md',
    'notes/stationary-lll-quantum-exterior.md',
    'notes/stationary-reference-nonbolt-stress.md',
    'notes/stationary-global-weighted-map.md',
    'notes/stationary-adjoint-phase.md',
]
CONTRACT={
    'acceptance':'every symbolic residual exactly zero and rational inequality true',
    'numerical_tolerances':'none',
    'scope':'seed contact identities, LLL and flat/common-local curvature algebra',
    'ordinary_arguments':['source ledger','form closure','Feynman integration',
        'angular meromorphic continuation and convergence','basepoint decision'],
    'unresolved':['full curved quantum lower form','physical kernel','adjoint phase',
        'source moments','quantitative coupled inverse','nonlinear and physical matching'],
    'failure_policy':'append every run and retain failed input bytes',
    'hash_policy':'SHA256 of raw UTF8/LF/noBOM bytes, no normalization',
}
ROWS=[]
R=S.Rational


def check(label,condition,detail=''):
    row={'label':label,'passed':bool(condition),'detail':str(detail)}
    ROWS.append(row)
    print(json.dumps(row,sort_keys=True),flush=True)


def zero(label,expression):
    value=S.simplify(expression)
    check(label,value==0,value)


def gauge_and_contacts():
    aa,rr,f2,b,c,cp,eta=S.symbols('A R Frr b c cp eta',real=True)
    d=b+c; e=b-c
    ap=f2/2-aa**2; rp=aa*rr-rr**2; w=aa+2*rr
    bp=-2*aa*d-2*rr*e; ep=aa*eta+d
    zero('TT_d_constraint',bp+cp+2*(aa+rr)*d-cp-4*rr*c)
    zero('TT_Ad_integrated_coefficient',2*aa*(aa+rr)-(ap+w*aa)/2-(2*aa**2+aa*rr-f2/4))
    zero('Y_energy_potential',ap+w*aa+aa**2-(f2/2+aa**2+2*aa*rr))
    zero('Y_second_radial',ap*eta+aa*ep+bp+cp-(f2*eta/2-aa*d-2*rr*d+cp+4*rr*c))
    zero('Y_second_time_first',ap*eta+aa*ep-(f2*eta/2+aa*d))
    zero('Y_second_time_second',aa*ep-aa**2*eta-aa*d)
    zero('Y_second_sphere_first',rp*eta+rr*ep-((2*aa*rr-rr**2)*eta+rr*d))
    zero('Y_second_sphere_second',rr*ep-rr**2*eta-((aa*rr-rr**2)*eta+rr*d))
    tprime=-(bp+3*cp)/2-(ap-rp)*eta-(aa-rr)*ep
    ee=f2/2-2*aa*rr+rr**2
    zero('LLL_exact_T_derivative',tprime-(-3*cp/2+2*rr*b-ee*eta))
    zero('LLL_E_coefficient',ap-rp+aa*(aa-rr)-ee)
    check('seed_Frr_lower',R(52*1000,27)-R(1,8)-2-50>1000)
    check('seed_Frr_upper',2*1000+6*(1+3*1000)+50<21*1000)
    check('Ad_square_bound',40+R(21*9,4)<88)
    check('TT_gradient_ledger',32*R(1,2)+10==26 and 32*88+28*5*9==4076)
    check('TT_gradient_energy',26*2<=5200 and 26*40+4076<5200)
    check('Y_energy_L2',2*9==18<20<25)
    check('Y_first_derivative',9+9+2*5*18<15**2)
    check('Y_radial_square_ledger',5*(R(21**2,4)*20+88+4*5*9+16*5)==12765)
    check('Y_time_squares_ledger',2*(R(21**2,4)*20+88)==4586 and 2*88==176)
    check('Y_sphere_squares_ledger',2*2*(27**2*20+5*9)==58500 and 4*2*(16**2*20+5*9)==41320)
    check('Y_full_Hessian_energy',12765+4586+176+58500+41320+R(19,2)*40==117727<160000 and 19<160000)
    check('return_tensor_bounds',22<25 and 2*5+30<=40 and 5+30==35 and 2*73+800<1000 and 5200<73**2)
    check('contact_physical_volume_constants',49600*R(22,7)<200000 and 9800*R(22,7)<40000)
    check('complete_contact_budget',300000*R(10**16,1)*R(4,10**27)/12==R(1,10**6))
    check('LLL_gradient_ledger',3*R(9,4)<7 and 3*4*5*4+3*39**2*18<83000)
    check('LLL_energy_bound',14<10**5 and 7*40+83000<10**5)
    check('LLL_benchmark_budget',R(10**5*400,18)*R(85,10**34)<R(1,10**25))


def seed_source_and_flat():
    t,x,a,qc,d=S.symbols('t x a qc d',positive=True)
    ff=1+x*t**2-2*a/t+qc/t**2-2*d/t**3+d/t**4
    scalar=-S.diff(ff,t,2)-4*S.diff(ff,t)/t+2*(1-ff)/t**2
    zero('seed_scalar_curvature',scalar-(-12*x+4*d/t**5-6*d/t**6))
    lap=S.diff(t**2*ff*S.diff(scalar,t),t)/t**2
    zero('seed_scalar_laplacian',lap-((S.diff(ff,t)+2*ff/t)*S.diff(scalar,t)+ff*S.diff(scalar,t,2)))
    # For |d|<=1, |R_t|<=56/t^6 and |R_tt|<=372/t^7.
    check('boxR_coefficient_guard',(13+2*5)*56+5*372<10**5)
    check('source_orbit_sum',R(100005,100000)*10**15+11000+1+54000<2*10**15)
    check('source_angular_sum',R(100005,100000)*10**15+11000+1+54000+500000<2*10**15)
    check('anomaly_sectional_guards',24*11**2==2904 and 4*2904==11616 and 4*11616==46464)
    check('gradient_Frobenius_reserve',2*10**15<10**16)
    v=S.symbols('v',nonnegative=True)
    rem=S.log(1+v)-v+v**2/2
    zero('flat_remainder_positive_derivative',S.diff(rem,v)-v**2/(1+v))
    zero('flat_remainder_zero_origin',rem.subs(v,0))
    zero('flat_cubic_upper_derivative',S.diff(v**3/3-rem,v)-v**3/(1+v))
    al=S.symbols('alpha',real=True)
    zero('flat_p2_parameter_moment',S.integrate(al*(1-al),(al,0,1))-R(1,6))
    zero('flat_p4_parameter_moment',S.integrate(al**2*(1-al)**2,(al,0,1))-R(1,30))
    zero('flat_p6_parameter_moment',S.integrate(al**3*(1-al)**3,(al,0,1))-R(1,140))
    mm,j,k=S.symbols('M j k',real=True)
    sym=mm**2*(mm**2-j*k)/((mm**2+j**2)*(mm**2+k**2))-mm**2*(1/(mm**2+j**2)+1/(mm**2+k**2))/2
    zero('raw_Dirac_symmetrized_sign',sym+mm**2*(j+k)**2/(2*(mm**2+j**2)*(mm**2+k**2)))
    nn,qq=S.symbols('n q',positive=True)
    zero('angular_FP_telescoping_term',2*(qq+2*nn)/(nn*(nn+qq))-4/nn-2*(1/(nn+qq)-1/nn))
    zero('angular_M4_telescoping_term',2*(qq+2*nn)/(nn**2*(nn+qq)**2)-2*(1/nn**2-1/(nn+qq)**2)/qq)
    ss,fg,ll,zz=S.symbols('s Fang logscale next')
    heat=S.exp(2*ss*ll)*(2+fg*ss+zz*ss**2)
    zero('common_heat_finite_coefficient',S.diff(heat,ss).subs(ss,0)-(fg+4*ll))
    check('flat_positive_margin',R(25,12)==sum(R(1,i) for i in range(1,5)) and (R(1,6)-R(4,10**21))/120>R(1,1000))
    rv,tt,dt=S.symbols('R T DeltaT',real=True)
    curv=S.exp(-2*ss*tt)*(rv-2*ss*dt)**2
    zero('curved_a2_second_variation',S.diff(curv,ss,2).subs(ss,0)-(4*rv**2*tt**2+16*rv*tt*dt+8*dt**2))
    check('curved_a2_geometry',R(13,5)+5*17<100)
    check('curved_a2_gradient_constant',R(32*100,960*3)<2)
    check('curved_a2_potential_constant',R((4*2**2+8*100)*100,960*3)<30)
    check('curved_a2_energy_constant',400000*R(61,10)+30000*R(61,10)**2<4000000<10**8)


def hashes(data):
    return {name:hashlib.sha256(raw).hexdigest() for name,raw in data.items()}


def main():
    before={name:(ROOT/name).read_bytes() for name in INPUTS}
    existing=json.loads(RECEIPT.read_text(encoding='utf-8')) if RECEIPT.exists() else {'runs':[]}
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'contract':CONTRACT,
         'command':sys.orig_argv,'dont_write_bytecode':sys.dont_write_bytecode,'cwd':str(Path.cwd()),
         'python':platform.python_version(),'sympy':S.__version__,'execution_inputs_sha256':hashes(before)}
    with LOG.open('a',encoding='utf-8',newline='\n') as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
        print('RUN '+json.dumps(run,sort_keys=True),flush=True)
        try:
            check('actual_bytecode_disabled',sys.dont_write_bytecode and '-B' in sys.orig_argv)
            check('inputs_strict_UTF8_LF',all(raw.decode('utf-8') is not None and b'\r' not in raw
                  and not raw.startswith(b'\xef\xbb\xbf') for raw in before.values()))
            gauge_and_contacts()
            seed_source_and_flat()
        except Exception:
            traceback.print_exc()
            check('uncaught_exception',False)
        run['final_inputs_sha256']=hashes({name:(ROOT/name).read_bytes() for name in INPUTS})
        print('FINAL_INPUTS '+json.dumps(run['final_inputs_sha256'],sort_keys=True),flush=True)
        check('frozen_inputs_unchanged',run['execution_inputs_sha256']==run['final_inputs_sha256'])
        run['checks']=ROWS
        run['failed_checks']=[row for row in ROWS if not row['passed']]
        code=int(bool(run['failed_checks']))
        run['status']='failed' if code else 'passed'
        run['terminal_exit_status']=code
        run['finished_utc']=datetime.now(timezone.utc).isoformat()
        if code:
            run['failed_inputs_base64']={name:base64.b64encode(raw).decode('ascii') for name,raw in before.items()}
            for name,raw in before.items():
                print(json.dumps({'path':name,'sha256':hashlib.sha256(raw).hexdigest(),
                                  'base64':base64.b64encode(raw).decode('ascii')},sort_keys=True))
        print('TERMINAL '+json.dumps({'status':run['status'],'checks':len(ROWS),
              'failures':len(run['failed_checks']),'terminal_exit_status':code,
              'finished_utc':run['finished_utc']},sort_keys=True),flush=True)
    existing['runs'].append(run)
    RECEIPT.write_text(json.dumps(existing,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'failures':len(run['failed_checks']),
                      'receipt':str(RECEIPT)},sort_keys=True))
    return code


if __name__=='__main__':
    sys.exit(main())
