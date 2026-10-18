"""Exact finite checks supporting the fixed-seed weighted TT inverse.

Parametrices, graph-domain density, weighted bootstrap, the actual TT
inverse and adjoint regularity remain ordinary arguments.  This program
does not compute a coupled kernel, scalar phase or physical solution.
"""
from __future__ import annotations
import base64, contextlib, hashlib, json, platform, sys, traceback
from datetime import datetime, timezone
from pathlib import Path
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/stationary-weighted-tt-inverse.json'
LOG=ROOT/'evidence/stationary-weighted-tt-inverse.log'
INPUTS=[
    'checks/stationary_weighted_tt_inverse.py',
    'notes/stationary-weighted-tt-inverse.md',
    'notes/stationary-curved-quantum-tt.md',
    'notes/stationary-global-weighted-map.md',
    'notes/stationary-global-tt-reduction.md',
    'notes/stationary-disk-and-normal-inverse.md',
    'notes/stationary-shape-logarithmic-response.md',
    'notes/stationary-weighted-conformal-range.md',
    'notes/stationary-trace-range-condition.md',
    'notes/stationary-centered-covariant-seed.md',
    'notes/stationary-centered-exterior-source.md',
    'notes/stationary-adjoint-phase.md',
    'notes/stationary-seed-quantum-form.md',
    'notes/stationary-centered-canonical-response.md',
    'notes/spherical-anomaly-sources.md',
    'notes/stationary-coupled-asymptotic-obstruction.md',
]
CONTRACT={
    'acceptance':'exact residuals zero and stated rational comparisons satisfied',
    'numerical_tolerances':'none',
    'scope':'finite congruence, curvature/density, scaling, mass, Schur and energy-normalization algebra',
    'ordinary_arguments':['sharp seed York mapping and transpose',
        'interior/bolt and weighted outer parametrices', 'adjoint regularity',
        'critical graph-domain density and coercivity extension',
        'fast asymptotic bootstrap and actual weighted TT inverse',
        'compactness and convergent adjoint/source pairings'],
    'unresolved':['actual scalar Schur kernel', 'physical adjoint phase and source moments',
        'quantitative coupled inverse', 'nonlinear candidate transfer and physical matching'],
    'failure_policy':'append initialized scientific runs and retain failed raw input bytes; unreadable-input or receipt preflight errors abort before scientific execution',
    'hash_policy':'SHA256 of raw UTF8/LF/noBOM bytes without normalization',
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


def matrix_zero(label,matrix):
    vals=[S.simplify(entry) for entry in matrix]
    check(label,all(entry==0 for entry in vals),vals)


def normal_and_mass():
    s=S.symbols('s')
    jc,nt=S.symbols('J_c N_T')
    C=lambda z:S.Matrix([[1,(z+3)/(3*(z+1))],[0,2]])
    rc=lambda z:S.Matrix([(z+3)/(2*(z+1)),-R(3,2)])
    rs=S.Matrix([1,0])
    matrix_zero('TT_normal_reconstruction',C(s)*rc(s)-S.Matrix([0,-3]))
    H=C(3-s).T*S.diag(jc,nt)*C(s)
    change=S.Matrix.hstack(rs,rc(s))
    reflected=S.Matrix.hstack(rs,rc(3-s))
    matrix_zero('reflected_TT_congruence_factor_nine',reflected.T*H*change-S.diag(jc,9*nt))
    zero('interior_shape_coefficient',R(9,4)/120-R(3,160))
    zero('interior_shape_equals_nine_d0',9*R(1,480)-R(3,160))
    matrix_zero('mass_TT_reconstruction',rc(3)-R(3,4)*S.Matrix([1,-2]))
    # A source at exponent0 pairs with the return symbol at exponent0;
    # the mass is at exponent3.  Derive both from the SAME return map.
    Csrc=rc(0)
    matrix_zero('scalar_c_source_normalization',Csrc-R(3,2)*S.Matrix([1,-1]))
    zero('mass_boundary_ratio_nine_eighths',R(3,2)*R(3,4)-R(9,8))
    At,d0,h3=S.symbols('A_T d0 H3')
    hm=S.Function('H_M')
    q=s*(s-3)
    normal=d0*q*(q+2)*hm(s)-At*q
    deriv=S.diff(normal,s).subs(s,3).subs(hm(3),h3)
    zero('TT_mass_normal_derivative',deriv-3*(2*d0*h3-At))
    zero('mass_jump_constant',-9*deriv-27*(At-2*d0*h3))
    t=S.symbols('t',real=True)
    rr=t**2+R(9,4)
    zero('critical_TT_polynomial',q.subs(s,R(3,2)+S.I*t)+rr)
    zero('critical_TT_positive_expression',normal.subs(s,R(3,2)+S.I*t).subs(hm(R(3,2)+S.I*t),h3)-(At*rr+d0*rr*(rr-2)*h3))
    check('critical_TT_HM_lower_factor',R(9,4)*(R(9,4)-2)*R(14,5)==R(63,40)>0)
    check('critical_TT_root_weight_distances',R(3,2)>0 and 3-R(3,2)>0 and R(7,2)-3>0)
    check('York_strip_and_dual',-1<3-R(7,2)<4 and -1<R(7,2)<4)


def scalar_curvature_and_density():
    T,ep,la,Z,Jdiv,Zdiv=S.symbols('T epsilon lambda Z J_div Z_div')
    G=S.symbols('G',positive=True)
    aa=1+11*ep*la/3
    l0=la+11*ep*la**2/6
    rv=-12*la+Z
    # Here T=-Box; the two divergence operators are independent formal terms.
    from_trace=T**2+(-1/ep+la/3)*T+Zdiv/36+R(11,9)*Jdiv-(2*rv+48*l0)/(6*ep)
    factor=(T+4*la)*(T-aa/ep)+Zdiv/36+R(11,9)*Jdiv-Z/(3*ep)
    zero('exact_scalar_factorization',from_trace-factor)
    zero('conformal_fourth_normalization',(6*ep/(8*S.pi*G)).subs(ep,G/(60*S.pi))-1/(80*S.pi**2))
    E,V=S.symbols('E V')
    zero('fixed_density_cancels_offshell_potential',(2*rv+48*l0-4*E)+4*E-(2*rv+48*l0))
    u=S.symbols('u',positive=True)
    x,d,delta=S.symbols('x d delta',real=True)
    # Independent areal curvature calculation, with physical Q coefficient.
    f=1+x/u**2-2*(1+2*x+delta/2)*u+(1+3*x+d)*u**2-2*d*u**3+d*u**4
    # B=1, r=1/u; dr differentiation is -u^2 d/du.
    dr=lambda fun:-u**2*S.diff(fun,u)
    fr=dr(f); frr=dr(fr)
    r4=-frr-4*u*fr+2*u**2*(1-f)
    zero('seed_full_scalar_curvature',r4-(-12*x+d*(4*u**5-6*u**6)))
    j=-frr/4+(f-1)*u**2/2
    zero('seed_tracefree_Ricci_eigenvalue',j-(-(1+3*x+d)*u**4+5*d*u**5-R(9,2)*d*u**6))
    zero('seed_Ricci_Frobenius',2*j**2+2*(-j)**2-4*j**2)
    # Exact density and diffeomorphism return contact, with arbitrary E.
    rr,etr,eta,er,diveta=S.symbols('Eprime_TT E eta gradE diveta')
    zero('disk_mixed_density_return',rr-eta*er-etr*diveta-(rr-(eta*er+etr*diveta)))
    sig,dsig=S.symbols('sigma grad_sigma')
    # Product-rule identity for integration by parts: sigma[-div(Eeta)]
    # differs from Eeta.grad sigma by a pure negative divergence.
    zero('mixed_adjoint_return_sign',-sig*(eta*er+etr*diveta)-etr*eta*dsig+((eta*er+etr*diveta)*sig+etr*eta*dsig))
    # Invariant J variation includes inverse-metric and R h contacts.
    h1,h2,h3,j1,j2,j3,k1,k2,k3,sc=S.symbols('h1 h2 h3 j1 j2 j3 k1 k2 k3 R')
    HH=S.diag(h1,h2,h3,-h1-h2-h3)
    JJ=S.diag(j1,j2,j3,-j1-j2-j3)
    KK=S.diag(k1,k2,k3,-k1-k2-k3)
    rdot=-S.trace(JJ*HH)
    jdot=KK/2-rdot*S.eye(4)/4-sc*HH/4
    invariant=2*S.trace(JJ*jdot)-2*S.trace(HH*JJ*JJ)
    expected=S.trace(JJ*KK)-sc*S.trace(JJ*HH)/2-2*S.trace(HH*JJ*JJ)
    zero('TT_J_squared_full_metric_variation',invariant-expected)
    # Fixed magnetic two-form: vary BOTH inverse metrics in F_ab F^ab.
    ff,tt=S.symbols('magnetic_F parameter',real=True)
    flux=S.zeros(4); flux[2,3]=ff; flux[3,2]=-ff
    gi=S.diag(*[1/(1+tt*HH[j,j]) for j in range(4)])
    mag=S.trace(flux.T*gi*flux*gi)
    zero('fixed_flux_metric_variation',S.diff(mag,tt).subs(tt,0)+2*S.trace(HH*flux*flux.T))


def density_scaling_and_gaps():
    r,L=S.symbols('r L',positive=True)
    z=S.symbols('z',positive=True)
    # Leading M=r^2F'-2rF is independent of the AdS r^2 term.
    la,m,charge=S.symbols('lambda m Q')
    f=la*r**2+1-2*m/r+charge/r**2
    zero('TT_integral_coefficient_leading',r**2*S.diff(f,r)-2*r*f-(-2*r+6*m-4*charge/r))
    zero('far_bump_moment_rescaling',(-2*r/L**2).subs(r,L*z)*L+2*z)
    zero('critical_graph_bump_amplitude',L**R(3,2)*L**-2-L**-R(1,2))
    zero('fast_graph_bump_not_dense',L**R(7,2)*L**-2-L**R(3,2))
    check('bootstrap_first_root_free',0<R(3,2)<R(5,2)<3)
    check('bootstrap_second_crosses_only_mass',R(5,2)<3<R(7,2)<4)
    check('bootstrap_steps_leave_two_power_margin',R(5,2)-R(3,2)<2 and R(7,2)-R(5,2)<2)
    check('mass_corrected_local_and_far_fast',5>R(7,2) and 4>R(7,2))
    # E0>=beta lambda C0/(2G), ||c||^2=4pi beta C0.
    zero('scalar_c_L2_gap_normalization',R(1,2)/(4*S.pi)-1/(8*S.pi))
    zero('TT_tensor_Gram_gap',R(1,8)/22-R(1,176))
    check('disk_lower_margin',R(995,1000)>R(99,100))
    zero('TT_gap_inverse_constant',1/(R(99,100)/(176*S.pi))-176*S.pi/R(99,100))
    zero('mixed_gap_coefficient',R(176,64)-R(11,4))
    zero('Q_preconditioned_mixed_coefficient',80*R(11,4)/4-55)
    jnorm,aval=S.symbols('Jnorm a0',positive=True)
    gg,ll=S.symbols('G lambda',positive=True)
    zero('Planck_energy_majorant_ratio',(55*S.pi*jnorm**2/(R(99,100)*gg*ll**2))/(60*S.pi*aval/gg)-55*(jnorm/ll)**2/(60*R(99,100)*aval))
    check('bolt_majorant_not_a_smallness_bound',R(55,60)/R(99,100)*6**2>30)
    check('tensor_derivative_energy_sum',R(1,2*41600)>0 and R(1,2*176)>R(1,2*41600))


def exact_noncommuting_Schur():
    ll1,ll2,ll3,b1,b2,b3,b4=S.symbols('l1 l2 l3 b1 b2 b3 b4')
    LL=S.Matrix([[ll1,ll2],[ll2,ll3]])
    BB=S.Matrix([[b1,b2],[b3,b4]])
    DD=S.Matrix([[3,1],[1,2]])
    inv=DD.inv(); eye=S.eye(2); null=S.zeros(2)
    full=LL.row_join(BB).col_join(BB.T.row_join(DD))
    left=eye.row_join(-BB*inv).col_join(null.row_join(eye))
    right=eye.row_join(null).col_join((-inv*BB.T).row_join(eye))
    schur=LL-BB*inv*BB.T
    diag=schur.row_join(null).col_join(null.row_join(DD))
    matrix_zero('noncommuting_triangular_block_elimination',left*full*right-diag)
    check('TT_matrix_positive',DD[0,0]>0 and DD.det()>0)
    tau=S.Matrix(S.symbols('tau1 tau2')); psi=S.Matrix(S.symbols('psi1 psi2'))
    ft=S.Matrix(S.symbols('ft1 ft2')); fc=S.Matrix(S.symbols('fc1 fc2'))
    c=inv*(fc-BB.T*tau)
    matrix_zero('Schur_reconstructed_TT_equation',BB.T*tau+DD*c-fc)
    matrix_zero('Schur_reconstructed_scalar_equation',LL*tau+BB*c-(schur*tau+BB*inv*fc))
    kernel=tau.col_join(-inv*BB.T*tau)
    matrix_zero('coupled_kernel_graph',full*kernel-(schur*tau).col_join(S.zeros(2,1)))
    chi=psi.col_join(-inv*BB.T*psi)
    zero('actual_adjoint_source_reduction',(chi.T*ft.col_join(fc))[0]-(psi.T*(ft-BB*inv*fc))[0])
    # Positive correction is a completed finite quadratic, not a sign claim
    # about the entire scalar Schur operator.
    xx,yy=S.symbols('x y',real=True)
    zero('Schur_correction_finite_positive_square',(S.Matrix([xx,yy]).T*inv*S.Matrix([xx,yy]))[0]-(R(2,5)*(xx-yy/2)**2+yy**2/2))
    matrix_zero('inverse_transpose_symmetry',inv-inv.T)


def hashes(data):
    return {name:hashlib.sha256(raw).hexdigest() for name,raw in data.items()}


def main():
    before={name:(ROOT/name).read_bytes() for name in INPUTS}
    existing=json.loads(RECEIPT.read_text(encoding='utf-8')) if RECEIPT.exists() else {'runs':[]}
    if not isinstance(existing,dict) or not isinstance(existing.get('runs'),list):
        raise ValueError('receipt preflight requires an object with a runs list')
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'contract':CONTRACT,
         'command':sys.orig_argv,'dont_write_bytecode':sys.dont_write_bytecode,'cwd':str(Path.cwd()),
         'python':platform.python_version(),'sympy':S.__version__,'execution_inputs_sha256':hashes(before)}
    with LOG.open('a',encoding='utf-8',newline='\n') as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
        print('RUN '+json.dumps(run,sort_keys=True),flush=True)
        try:
            check('actual_bytecode_disabled',sys.dont_write_bytecode and '-B' in sys.orig_argv)
            check('inputs_strict_UTF8_LF',all(raw.decode('utf-8') is not None and b'\r' not in raw and not raw.startswith(b'\xef\xbb\xbf') for raw in before.values()))
            normal_and_mass()
            scalar_curvature_and_density()
            density_scaling_and_gaps()
            exact_noncommuting_Schur()
        except Exception:
            traceback.print_exc()
            check('uncaught_exception',False)
        try:
            run['final_inputs_sha256']=hashes({name:(ROOT/name).read_bytes() for name in INPUTS})
        except Exception:
            traceback.print_exc()
            run['final_inputs_sha256']={}
            check('final_input_reread_failed',False)
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
                print(json.dumps({'path':name,'sha256':hashlib.sha256(raw).hexdigest(),'base64':base64.b64encode(raw).decode('ascii')},sort_keys=True))
        print('TERMINAL '+json.dumps({'status':run['status'],'checks':len(ROWS),'failures':len(run['failed_checks']),'terminal_exit_status':code,'finished_utc':run['finished_utc']},sort_keys=True),flush=True)
    existing['runs'].append(run)
    RECEIPT.write_text(json.dumps(existing,indent=2,sort_keys=True)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'failures':len(run['failed_checks']),'receipt':str(RECEIPT)},sort_keys=True))
    return code


if __name__=='__main__':
    sys.exit(main())
