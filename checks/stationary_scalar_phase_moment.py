"""Exact finite identities supporting the actual scalar phase and moment.

Proposed verbatim destination: checks/stationary_scalar_phase_moment.py.
Generalized-field domains, PV/LAP limits, global adjoint identification,
Green boundary limits and source admissibility remain ordinary proofs.
This checker does not evaluate the actual logarithmic phase, the signed
preparation moment, a fast inverse norm or a nonlinear physical solution.
"""
from __future__ import annotations
import base64, contextlib, hashlib, json, platform, sys, traceback
from datetime import datetime, timezone
from pathlib import Path
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/stationary-scalar-phase-moment.json'
LOG=ROOT/'evidence/stationary-scalar-phase-moment.log'
INPUTS=[
    'checks/stationary_scalar_phase_moment.py',
    'notes/stationary-scalar-phase-moment.md',
    'notes/stationary-scalar-outgoing-inverse.md',
    'notes/stationary-weighted-tt-inverse.md',
    'notes/stationary-adjoint-phase.md',
    'notes/stationary-global-weighted-map.md',
    'notes/stationary-global-tt-reduction.md',
    'notes/stationary-disk-and-normal-inverse.md',
    'notes/stationary-weighted-conformal-range.md',
    'notes/stationary-centered-covariant-seed.md',
    'notes/stationary-centered-first-order-geometry.md',
    'notes/stationary-centered-exterior-source.md',
    'notes/stationary-centered-complete-residual.md',
    'notes/uniform-dirac-reference.md',
    'notes/stationary-trace-range-condition.md',
]
CONTRACT={
    'acceptance':'exact residuals zero and stated rational comparisons satisfied',
    'numerical_tolerances':'none',
    'scope':'finite spectral-normalization, noncommuting matrix, rank-one, phase-budget, coordinate, trace-projection and Gaussian-coefficient identities',
    'ordinary_arguments':['regular-bolt generalized scalar fields and their weighted domains',
        'outgoing/incoming limiting absorption and the real principal-value map',
        'rank-one spectral density and source-amplitude limits as operator statements',
        'actual TT/York reconstruction and global adjoint-line identification',
        'convergent seed quadrature and vanishing bolt/outer Green terms',
        'actual source admissibility and limitations of unsigned estimates'],
    'unresolved':['numerical logarithmic phase at the exact retained coefficients',
        'signed preparation or finite-candidate moment and its zero/nonzero value',
        'numerical compatible fast weighted inverse constant',
        'nonlinear control or justified spatial reduction and physical matching'],
    'matrix_scope':'generic symmetric matrix algebra plus a singular mixed-sign fixture; no matrix sample is an operator-domain proof',
    'failure_policy':'append initialized scientific runs and retain failed raw input bytes; unreadable-input or receipt preflight errors abort before scientific execution',
    'hash_policy':'SHA256 of raw UTF8/LF/noBOM bytes without normalization',
}
ROWS=[]
R=S.Rational


def check(label,condition,detail=''):
    if any(row['label']==label for row in ROWS):
        raise ValueError('duplicate check label: '+label)
    row={'label':label,'passed':bool(condition),'detail':str(detail)}
    ROWS.append(row)
    print(json.dumps(row,sort_keys=True),flush=True)


def zero(label,expression):
    expression=S.sympify(expression)
    if expression.atoms(S.Float):
        raise ValueError('inexact Float atom in '+label)
    value=S.simplify(S.cancel(expression))
    check(label,value==0,value)


def matrix_zero(label,matrix):
    if any(entry.atoms(S.Float) for entry in matrix):
        raise ValueError('inexact Float atom in '+label)
    vals=[S.simplify(S.cancel(entry)) for entry in matrix]
    check(label,all(entry==0 for entry in vals),vals)


def spectral_and_physical_normalization():
    r,rho,sigma,beta,F,k,G,la=S.symbols('r rho sigma beta F k G lambda',positive=True)
    delta=S.symbols('delta',real=True)
    measure_factor=4*S.pi*beta*r**2*S.sqrt(F)
    zero('Liouville_physical_measure_factor',measure_factor/S.sqrt(F)-4*S.pi*beta*r**2)
    E=60*S.pi/G+11*la/3
    zero('radial_energy_shift',E-R(9,4)*la-(60*S.pi/G+R(17,12)*la))
    theta=k*rho+delta
    regular=(S.exp(S.I*theta)+S.exp(-S.I*theta))/2
    jost=S.exp(S.I*k*rho)
    wr=regular*S.diff(jost,rho)-S.diff(regular,rho)*jost
    zero('unit_cosine_Jost_Wronskian',wr-S.I*k*S.exp(-S.I*delta))
    zero('Green_derivative_jump_sign',(-1/(S.I*k*S.exp(-S.I*delta)))*(S.I*k*S.exp(-S.I*delta))+1)
    zero('outgoing_Green_coefficient',-1/(S.I*k*S.exp(-S.I*delta))-S.I*S.exp(S.I*delta)/k)
    tsmall,tlarge=S.symbols('theta_small theta_large',real=True)
    green=S.I*S.exp(S.I*tlarge)*S.cos(tsmall)/k
    zero('unit_cosine_spectral_density',S.im(green)-S.cos(tsmall)*S.cos(tlarge)/k)
    zero('principal_value_far_sign',S.re(green)+S.sin(tlarge)*S.cos(tsmall)/k)
    Bnorm=S.sqrt(S.pi*k*rho/2)*S.sqrt(S.pi*k*sigma/2)
    zero('Bessel_spectral_rank_one_coefficient',(S.pi/2)*S.sqrt(rho*sigma)/Bnorm-1/k)
    zero('Nicholson_unit_amplitude_bound',S.pi*k*rho/2*2/(S.pi*k*rho)-1)
    t,e=S.symbols('T E')
    r4=(1/(t-e)-1/(t+4*la))/(e+4*la)
    zero('same_operator_fourth_resolvent',r4-1/((t-e)*(t+4*la)))
    zero('fourth_spectral_density_coefficient',1/(k*(e+4*la))-1/k/(e+4*la))
    v=S.Function('v')(rho); f=S.Function('f')(rho); b=S.Function('b')(rho)
    # D=partial-b and D*=-partial-b after the physical unitary map.
    # This checks the D* adjoint integration identity. Physical div is
    # -D*, so pairing div f with v gives the claimed minus <Dv,f>.
    zero('weak_divergence_pairing_sign',v*(-S.diff(f,rho)-b*f)-(S.diff(v,rho)-b*v)*f+S.diff(v*f,rho))


def real_noncommuting_blocks():
    a,b,d,l1,l2,l3,j1,j2,j3=S.symbols('a b d l1 l2 l3 j1 j2 j3',real=True)
    AA=S.Matrix([[a,b],[b,d]])
    LR=S.Matrix([[l1,l2],[l2,l3]])
    JJ=S.Matrix([[j1,j2],[j2,j3]])
    eye=S.eye(2)
    MR=eye+(AA-LR)*JJ
    matrix_zero('generic_MR_J_symmetry',MR.T*JJ-JJ*MR)
    # Adjugates check the inverse symmetry without dividing by J.
    adj=MR.adjugate(); det=MR.det()
    matrix_zero('generic_MR_adjugate_left',MR*adj-det*eye)
    matrix_zero('generic_MR_adjugate_right',adj*MR-det*eye)
    matrix_zero('generic_J_inverse_symmetry_numerator',JJ*adj-(JJ*adj).T)
    N=(eye+AA*JJ).adjugate()/(eye+AA*JJ).det()
    KR=N*LR*JJ
    matrix_zero('noncommuting_real_Neumann_factorization',(eye+AA*JJ)*(eye-KR)-MR)
    gg=S.Matrix(S.symbols('g1 g2',real=True))
    gR=adj*gg/det
    matrix_zero('real_generalized_field_vertex_reconstruction',gg+(LR-AA)*JJ*gR-gR)
    # This fixture has an actual zero eigenspace and does not commute.
    BB=S.Matrix([[2,1,1],[1,-1,2],[1,2,3]])
    J=S.diag(1,-1,0)
    M=S.eye(3)+BB*J
    Mi=M.inv()
    check('singular_sign_fixture',J.det()==0 and J**2!=S.eye(3) and BB*J!=J*BB)
    check('singular_sign_MR_invertible',M.det()==7)
    matrix_zero('singular_sign_J_inverse_symmetry',J*Mi-(J*Mi).T)
    return J,M,Mi


def rank_one_and_source_algebra(J,MR,MRi):
    c0=S.symbols('c0',positive=True)
    m0=S.symbols('m0',real=True)
    g=S.Matrix(S.symbols('g1 g2 g3',real=True))
    h=S.Matrix(S.symbols('h1 h2 h3',real=True))
    q=MRi*g
    aa=(g.T*J*q)[0]/c0
    zero('real_phase_coupling',S.conjugate(aa)-aa)
    M=MR-S.I*g*g.T*J/c0
    inverse=MRi+S.I*q*(g.T*J*MRi)/(c0*(1-S.I*aa))
    matrix_zero('singular_sign_Sherman_Morrison_left',M*inverse-S.eye(3))
    matrix_zero('singular_sign_Sherman_Morrison_right',inverse*M-S.eye(3))
    matrix_zero('rank_one_incident_inverse',M*(q/(1-S.I*aa))-g)
    zero('rank_one_real_source_cross_symmetry',(g.T*J*MRi*h)[0]-(q.T*J*h)[0])
    hplus=h+S.I*g*m0/c0
    y=inverse*hplus
    numerator=m0+(q.T*J*h)[0]
    zero('rank_one_outgoing_source_elimination',m0+(g.T*J*y)[0]-numerator/(1-S.I*aa))
    # Im L=gg*/c0 is positive, regardless of the signs in J.
    test=S.Matrix(S.symbols('x1 x2 x3',real=True))
    zero('rank_one_spectral_positive_square',(test.T*(g*g.T/c0)*test)[0]-(g.dot(test))**2/c0)
    ps,phi,a=S.symbols('Psi phi a',real=True)
    zero('incident_full_field_normalization',phi+(ps-phi+S.I*a*phi)/(1-S.I*a)-ps/(1-S.I*a))


def phase_amplitude_and_budget():
    a,theta,delta=S.symbols('a theta delta',real=True)
    c0=S.symbols('c0',positive=True)
    mm=S.symbols('moment_numerator',real=True)
    zero('phase_arctangent_cosine',S.cos(S.atan(a))-1/S.sqrt(1+a*a))
    zero('phase_arctangent_sine',S.sin(S.atan(a))-a/S.sqrt(1+a*a))
    zero('standing_wave_phase_sign',S.expand_trig(S.cos(theta+S.atan(a)))-(S.cos(theta)-a*S.sin(theta))/S.sqrt(1+a*a))
    scattering=(1+S.I*a)/(1-S.I*a)
    zero('scattering_ratio_unit_modulus',scattering*S.conjugate(scattering)-1)
    zero('incident_incoming_coefficient',(1-S.I*a)/(2*(1-S.I*a))-R(1,2))
    zero('incident_outgoing_coefficient',R(1,2)+S.I*a/(1-S.I*a)-scattering/2)
    zero('real_unit_cosine_from_incident',((1-S.I*a)/S.sqrt(1+a*a))/(1-S.I*a)-1/S.sqrt(1+a*a))
    outgoing=S.I*S.exp(S.I*delta)*mm/(c0*(1-S.I*a))
    target=S.I*S.exp(S.I*delta)*(1+S.I*a)/S.sqrt(1+a*a)*mm/(c0*S.sqrt(1+a*a))
    zero('normalized_outgoing_moment_coefficient',outgoing-target)
    f=S.symbols('coupled_moment',real=True)
    zero('coupled_source_80pi_squared_normalization',S.I*S.exp(S.I*delta)*(80*S.pi**2*f)/c0-80*S.pi**2*S.I*S.exp(S.I*delta)*f/c0)
    ell=R(4,10**12); normN=2; rho=R(8,10**12); eta=R(2,10**14)
    zero('mixed_Neumann_ratio',normN*ell-rho)
    check('mixed_Neumann_contraction',0<rho<1)
    check('real_auxiliary_inverse_budget',normN/(1-rho)<3)
    check('mixed_phase_budget',ell*normN/(1-rho)<R(81,10**13))
    check('geometric_Born_budget',R(6000,4*10**17)<eta)
    check('geometric_phase_budget',eta/(1-eta)<R(21,10**15))
    zero('arctangent_lipschitz_derivative_identity',1-S.diff(S.atan(a),a)-a*a/(1+a*a))
    check('combined_phase_enclosure',R(81,10**13)+R(21,10**15)<R(9,10**12))
    # A Neumann tail identity, not evaluation of its actual signed terms.
    q=S.symbols('q',positive=True); n=S.symbols('n',integer=True,nonnegative=True)
    zero('inverse_series_tail_geometric_identity',(1-q)*(2*q**(n+1)/(1-q))-2*q**(n+1))
    zero('phase_series_tail_coefficient',2*ell-R(8,10**12))
    gn,hn=S.symbols('g_norm h_norm',nonnegative=True)
    zero('moment_series_tail_coefficient',gn*hn*(2*q**(n+1)/(1-q))-2*gn*hn*q**(n+1)/(1-q))
    check('zeroth_phase_series_error_budget',R(8,10**12)*rho/(1-rho)<R(7,10**23))
    check('first_phase_series_error_budget',R(8,10**12)*rho**2/(1-rho)<R(6,10**34))
    p1,p2,p3,p4=S.symbols('p1 p2 p3 p4')
    K=S.Matrix([[p1,p2],[p3,p4]])
    finite=S.eye(2)+K+K**2+K**3
    matrix_zero('noncommuting_finite_Neumann_telescope',(S.eye(2)-K)*finite-(S.eye(2)-K**4))


def coordinate_and_boundary_factors():
    r,rho,s,z,u,B,cu,cz,kap,la,beta,k,rootF=S.symbols('r rho s z u B c_u c_z kappa lambda beta k rootF',positive=True)
    ell,delta=S.symbols('ell_inf delta',real=True)
    nu=k/S.sqrt(la)
    zero('disk_coordinate_leading_antiderivative',S.diff(kap/(la*r),r)+kap/(la*r**2))
    zero('disk_c_u_factor',(kap/(la*r))/(B/r)-kap/(la*B))
    zero('rescaled_disk_c_u_factor',cz*(kap/(la*B))-cz*kap/(la*B))
    phase_r=k*(ell-S.log(u)/S.sqrt(la))+delta
    Phi_u=k*ell+delta
    zero('proper_to_u_phase_sign',phase_r+(nu*S.log(u)-Phi_u))
    Phi_z=k*ell+nu*S.log(cu)+delta
    zero('u_to_z_phase_sign',S.expand_log((nu*S.log(u)-Phi_u).subs(u,z/cu),force=True)-(nu*S.log(z)-Phi_z))
    amp=(4*S.pi*beta*S.sqrt(la))**(-R(1,2))*r**(-R(3,2))
    target=(4*S.pi*beta*S.sqrt(la))**(-R(1,2))*(cu*B)**(-R(3,2))*z**R(3,2)
    zero('physical_disk_scalar_amplitude',amp.subs(r,cu*B/z)-target)
    zero('bolt_quadrature_regularization',2*s/S.sqrt(2*kap*s**2)-S.sqrt(2/kap))
    # rootF represents sqrt(F); the rationalization is exact.
    free=S.sqrt(la)*r
    zero('geometric_quadrature_reciprocal_identity',1/rootF-1/free+(rootF**2-free**2)/(rootF*free*(rootF+free)))
    check('geometric_quadrature_tail_constant',2<R(3,2)**2 and R(3,2)*3/2<3)
    zero('geometric_quadrature_tail_primitive',S.diff(3*S.sqrt(2)*B**2/(2*S.sqrt(la)*r**2),r)+3*S.sqrt(2)*B**2/(S.sqrt(la)*r**3))
    zero('trace_boundary_phi_Zprime_power',2+2-R(3,2)-6+R(7,2))
    zero('trace_boundary_phiprime_Z_power',2+2-R(5,2)-5+R(7,2))
    check('trace_boundary_negative_power',-R(7,2)<0)
    check('fast_scalar_moment_absolute_integrability',2-R(3,2)-4<-1)
    check('mixed_critical_input_weight_margin',R(1,2)+2==R(5,2)>R(3,2))
    check('fast_source_embeds_weak_weight',2*R(7,2)-6==1>-1)


def trace_projection_and_Gaussian():
    G,la,ep=S.symbols('G lambda epsilon',positive=True)
    R4,Z,TZ,C2,J2,M=S.symbols('R4 Z TZ C2 J2 M')
    atr=1+11*ep*la/3
    E=atr/ep
    l0=la+11*ep*la**2/6
    complete=R4+12*l0+ep*(-TZ-R(11,72)*R4**2+R(7,12)*C2+R(11,6)*J2-20*M)
    projected=(atr*Z-ep*TZ)+ep*(R(7,12)*C2+R(11,6)*J2-20*M-R(11,72)*Z**2)
    zero('full_trace_local_projection_split',complete.subs(R4,-12*la+Z)-projected)
    zero('reference_eigenfield_linear_trace_cancellation',atr-ep*E)
    zero('local_Newton_projection_one_sixth',(80*S.pi**2*ep/(8*S.pi*G)).subs(ep,G/(60*S.pi))-R(1,6))
    tr=S.symbols('trace_defect')
    zero('negative_seed_gradient_trace_sign',-(-tr/(8*S.pi*G))-tr/(8*S.pi*G))
    q,r=S.symbols('q r',positive=True)
    zero('physical_magnetic_projection_coefficient',-20*q**2/(2*r**4)+10*q**2/r**4)
    # Exact Gaussian-Bessel coefficient identity, with positive b,k.
    # Its integration/summation justification is the cited ordinary formula;
    # these checks do not evaluate a transform of the retained source.
    b,k=S.symbols('b k',positive=True)
    n=S.symbols('n',integer=True,nonnegative=True)
    Bessel_coeff=(-1)**n*(k/2)**(2*n)/S.factorial(n)**2
    radial_Gaussian_moment=S.factorial(n)/(2*b**(n+1))
    zero('Gaussian_Bessel_integrated_series_coefficient',Bessel_coeff*radial_Gaussian_moment-(-k*k/(4*b))**n/(2*b*S.factorial(n)))
    zero('Gaussian_Liouville_source_amplitude',S.sqrt(S.pi*k/2)*(S.exp(-k*k/(4*b))/(2*b))-S.sqrt(S.pi*k/2)/(2*b)*S.exp(-k*k/(4*b)))
    x=S.symbols('x',positive=True)
    zero('Gaussian_radial_moment_base',S.diff(-S.exp(-b*x*x)/(2*b),x)-x*S.exp(-b*x*x))
    t=S.symbols('t')
    zero('Gaussian_exponential_series_coefficient',S.diff(S.exp(t),t,4).subs(t,0)/S.factorial(4)-R(1,24))


def hashes(data):
    return {name:hashlib.sha256(raw).hexdigest() for name,raw in data.items()}


def main():
    before={name:(ROOT/name).read_bytes() for name in INPUTS}
    existing=json.loads(RECEIPT.read_text(encoding='utf-8')) if RECEIPT.exists() else {'runs':[]}
    if not isinstance(existing,dict) or not isinstance(existing.get('runs'),list):
        raise ValueError('receipt preflight requires an object with a runs list')
    ROWS.clear()
    run={'started_utc':datetime.now(timezone.utc).isoformat(),'contract':CONTRACT,
         'command':sys.orig_argv,'dont_write_bytecode':sys.dont_write_bytecode,'cwd':str(Path.cwd()),
         'python':platform.python_version(),'sympy':S.__version__,'execution_inputs_sha256':hashes(before)}
    with LOG.open('a',encoding='utf-8',newline='\n') as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
        print('RUN '+json.dumps(run,sort_keys=True),flush=True)
        try:
            check('actual_bytecode_disabled',sys.dont_write_bytecode and '-B' in sys.orig_argv)
            check('inputs_strict_UTF8_LF',all(raw.decode('utf-8') is not None and b'\r' not in raw and not raw.startswith(b'\xef\xbb\xbf') for raw in before.values()))
            spectral_and_physical_normalization()
            J,MR,MRi=real_noncommuting_blocks()
            rank_one_and_source_algebra(J,MR,MRi)
            phase_amplitude_and_budget()
            coordinate_and_boundary_factors()
            trace_projection_and_Gaussian()
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
