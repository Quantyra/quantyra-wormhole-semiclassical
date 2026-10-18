"""Exact finite identities supporting the fixed-seed scalar outgoing inverse.

The bolt realization, spectral form extension, complex Bessel/Born limits,
weak reconstruction and kernel exclusion remain ordinary analytic proofs.
These finite checks do not evaluate an adjoint phase or physical source
moment, or construct a nonlinear exterior/state solution.
"""
from __future__ import annotations
import base64, contextlib, hashlib, json, platform, sys, traceback
from datetime import datetime, timezone
from pathlib import Path
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/stationary-scalar-outgoing-inverse.json'
LOG=ROOT/'evidence/stationary-scalar-outgoing-inverse.log'
INPUTS=[
    'checks/stationary_scalar_outgoing_inverse.py',
    'notes/stationary-scalar-outgoing-inverse.md',
    'notes/stationary-weighted-tt-inverse.md',
    'notes/stationary-curved-quantum-tt.md',
    'notes/stationary-seed-quantum-form.md',
    'notes/stationary-centered-covariant-seed.md',
    'notes/stationary-centered-first-order-geometry.md',
    'notes/stationary-centered-exterior-source.md',
    'notes/stationary-trace-range-condition.md',
    'notes/stationary-global-weighted-map.md',
    'notes/stationary-global-tt-reduction.md',
    'notes/stationary-disk-and-normal-inverse.md',
    'notes/stationary-weighted-conformal-range.md',
    'notes/stationary-adjoint-phase.md',
    'notes/exact-spherical-em-reduction.md',
    'notes/stationary-reference-cold-bolt.md',
    'notes/spherical-anomaly-sources.md',
    'notes/stationary-centered-canonical-response.md',
]
CONTRACT={
    'acceptance':'exact residuals zero and stated rational comparisons satisfied; strict source bounds permit equal majorants',
    'numerical_tolerances':'none',
    'scope':'finite Liouville/curvature, endpoint symbols, Bessel moment, local vertex, resolvent/block and rational budget algebra',
    'ordinary_arguments':['regular-bolt selfadjoint scalar and actual TT realizations',
        'spectral form extension and York return mapping',
        'complex Bessel maximum-modulus and Born limiting-absorption proof',
        'physical derivative kernel and identity-contact estimates',
        'weak outgoing reconstruction and elliptic/domain bootstrap',
        'physical fast kernel exclusion and adjoint-range interpretation'],
    'unresolved':['actual adjoint phase and physical source moment',
        'numerical fast weighted inverse constant',
        'nonlinear candidate transfer and exterior/state/charge matching'],
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
    value=S.simplify(expression)
    check(label,value==0,value)


def matrix_zero(label,matrix):
    if any(entry.atoms(S.Float) for entry in matrix):
        raise ValueError('inexact Float atom in '+label)
    vals=[S.simplify(entry) for entry in matrix]
    check(label,all(entry==0 for entry in vals),vals)


def liouville_and_bolt():
    r,rho,la,kap,rh,f2=S.symbols('r rho lambda kappa rh Frr',positive=True)
    F=S.Function('F')(r)
    dr=lambda z:S.sqrt(F)*S.diff(z,r)
    w=dr(S.log(r**2*S.sqrt(F)))
    potential=dr(w)/2+w**2/4
    target=S.diff(F,r,2)/4+S.diff(F,r)/r-S.diff(F,r)**2/(16*F)
    zero('Liouville_potential_from_physical_measure',potential-target)
    v=S.Function('v')(r)
    D=lambda z:dr(z)-w*z/2
    Da=lambda z:-dr(z)-w*z/2
    zero('physical_gradient_closed_product',Da(D(v))+dr(dr(v))-target*v)
    y=S.Function('y')(rho)
    D0=lambda z:S.diff(z,rho)-z/(2*rho)
    D0a=lambda z:-S.diff(z,rho)-z/(2*rho)
    zero('critical_Bessel_product',D0a(D0(y))+S.diff(y,rho,2)+y/(4*rho**2))
    zero('partner_Bessel_product',D0(D0a(y))+S.diff(y,rho,2)-3*y/(4*rho**2))
    psi=S.Function('psi')(rho)
    zero('critical_form_radial_factor',D0(S.sqrt(rho)*psi)-S.sqrt(rho)*S.diff(psi,rho))
    rb=rh+kap*rho**2/2+kap*f2*rho**4/48
    rbp=S.diff(rb,rho)
    zero('bolt_radius_fourth_jet',S.expand(rbp**2-(2*kap*(rb-rh)+f2*(rb-rh)**2/2)).coeff(rho,4))
    wb=S.series(S.diff(S.log(rb**2*rbp),rho),rho,0,2).removeO()
    vb=S.expand(S.diff(wb,rho)/2+wb**2/4)
    zero('bolt_potential_finite_coefficient',vb.coeff(rho,0)-(f2/6+2*kap/rh))
    zero('bolt_potential_critical_coefficient',vb.coeff(rho,-2)+R(1,4))
    q=S.symbols('q')
    zero('critical_indicial_double_root',-q*(q-1)-R(1,4)+(q-R(1,2))**2)
    zero('partner_indicial_roots',-q*(q-1)+R(3,4)+(q-R(3,2))*(q+R(1,2)))
    check('critical_log_and_regular_L2_powers',2*R(1,2)>-1 and 2*(-R(1,2))==-1)
    check('U_global_absolute_reserve',R(21,4)+22-R(9,4)<30 and R(21,8)+R(9,4)<30)
    check('physical_gradient_shift_reserve',(R(13,4)+R(9,2))/2<4 and R(13,4)**2>R(21,2) and R(9,4)**2>5)


def outer_and_thermal():
    r,la=S.symbols('r lambda',positive=True)
    f=S.Function('f')(r)
    F=la*r**2+f
    zero('outer_curvature_square_cancellation',S.diff(F,r)**2-4*la*F-(4*la*r*S.diff(f,r)+S.diff(f,r)**2-4*la*f))
    t,x=S.symbols('t x',positive=True)
    lead=1-4/t**3+3/t**4
    zero('outer_leading_polynomial_derivative',S.diff(lead,t)-12*(t-1)/t**5)
    zero('outer_leading_value_at_two',lead.subs(t,2)-R(11,16))
    check('outer_F_lower_coefficient',R(11,16)-(R(1,4)+R(1,2))/(1000*4)>R(1,2))
    check('outer_f_coefficient_bound',R(11,4)+R(53,16)/1000<3)
    check('outer_rfprime_coefficient_bound',R(7,2)+R(7,2)/1000<4)
    check('outer_r2fsecond_coefficient_bound',R(17,2)+R(41,4)/1000<9)
    zero('outer_square_ratio_bound',8*4+R(2*4**2,4)+8*3-64)
    check('outer_Liouville_tail_coefficient',R(9,4)+4+4<11)
    # The sharper actual sqrt(2)<16/11 proves 11sqrt(2)/8<2.
    check('outer_integral_below_two_squared',2<R(16,11)**2)
    check('rho2_lower_log_guard',3<R(20,11)**2 and R(5,2)**2>5)
    # Exact exponential-series bounds supporting log(2000)<8, log(16)<3.
    elower=sum((R(1,S.factorial(j)) for j in range(5)),S.Integer(0))
    check('thermal_logarithm_rational_witness',elower**8>2000 and elower**3>16 and 3+10*8<84)
    check('thermal_center_half_margin',R(998,1000)*(1-R(1,10**20))*(1-R(12,10**30))>R(1,2))
    check('rho2_temperature_reserve',2*(2+84)==172)
    check('reference_L1_temperature_reserve',130+60*84<6000)
    kb=2000
    check('reference_Born_ratio',R(6000,4*10**17)<R(2,10**14)<R(1,10**12))
    check('benchmark_sqrt_Glambda_enclosures',R(3,10**17)**2<R(999,1000)*R(1,kb**10) and R(1,kb**10)<R(4,10**17)**2)
    check('benchmark_radial_frequency_lower',180*kb**10>(4*10**17)**2)
    check('weight_L2_reserve',172+1<200)
    check('weighted_zero_derivative_rounding',200/(1-R(1,10**12))<201)
    check('weighted_one_derivative_rounding',3*200==600)
    check('weighted_two_derivative_contact_rounding',600+R(1,4*10**17)<601)


def Bessel_moments_and_resolvents():
    aa=S.symbols('a',positive=True)
    # Mellin moment of K0, obtained by the positive exponential representation.
    moment=lambda j:2**S.Integer(j-1)*S.gamma(R(j+1,2))**2/aa**(j+1)
    zero('K0_moment_degree_zero',moment(0)-S.pi/(2*aa))
    zero('K0_moment_degree_one',moment(1)-1/aa**2)
    zero('K0_moment_degree_two',moment(2)-S.pi/(2*aa**3))
    z=S.symbols('z',positive=True)
    zero('Nicholson_order_zero_coefficient',8/S.pi**2*moment(0).subs(aa,2*z)-2/(S.pi*z))
    zero('Nicholson_order_one_first_majorant',8/S.pi**2*(moment(0)+2*moment(2)).subs(aa,2*z)-(2/(S.pi*z)+1/(S.pi*z**3)))
    zero('Nicholson_order_one_second_majorant',8/S.pi**2*(moment(0)+2*moment(1)).subs(aa,2*z)-(2/(S.pi*z)+4/(S.pi**2*z**2)))
    zero('Bessel_Green_derivative_jump',S.I*S.pi/2*2*S.I/S.pi+1)
    check('Bessel_small_argument_derivative_reserve',(2*R(22,7)+4)/4<4)
    check('Bessel_large_argument_derivative_reserve',R(6,4)<4)
    eta=R(1,10**12); q=R(1,4*10**17)
    check('actual_first_derivative_kernel_reserve',(2+4*q)/(1-eta)<3)
    check('actual_second_derivative_kernel_reserve',2+(4*eta+16*q+16*q*q)/(1-eta)<3)
    T,E,la=S.symbols('T E lambda')
    r4=(1/(T-E)-1/(T+4*la))/(E+4*la)
    zero('fourth_resolvent_exact_same_operator',r4-1/((T-E)*(T+4*la)))
    zero('fourth_resolvent_T_distribution_cancellation',T*r4-(E/(T-E)+4*la/(T+4*la))/(E+4*la))
    zero('homogeneous_component_reconstruction',(T+4*la)/(E+4*la)-(T-E)/(E+4*la)-1)
    zero('homogeneous_components_product',(T-E)*(T+4*la)/(E+4*la)+(T+4*la)*(-(T-E))/(E+4*la))


def seed_trace_and_vertex_budgets():
    u=S.symbols('u',positive=True)
    x,d,delta=S.symbols('x d delta',real=True)
    F=1+x/u**2-2*(1+2*x+delta/2)*u+(1+3*x+d)*u**2-2*d*u**3+d*u**4
    dr=lambda a:-u**2*S.diff(a,u)
    Fr=dr(F); Frr=dr(Fr)
    rr=-Frr-4*u*Fr+2*u**2*(1-F)
    zero('seed_scalar_curvature_from_areal_metric',rr+12*x-d*(4*u**5-6*u**6))
    j=-Frr/4+(F-1)*u**2/2
    zero('seed_tracefree_Ricci_eigenvalue',j+(1+3*x+d)*u**4-5*d*u**5+R(9,2)*d*u**6)
    wc=(Frr/u**2-2*Fr/u+2*(F-1))*u**2/12
    zero('seed_Weyl_amplitude_from_curvature',wc-(-(1+2*x+delta/2)*u**3+(1+3*x+d)*u**4-R(10,3)*d*u**5+R(5,2)*d*u**6))
    # Both estimates are necessary: the factor2 bound has u^4, not u^5.
    zero('Z_weight_four_lower_endpoint_factor',2+u*(4-6*u)-2*(1-u)*(1+3*u))
    # Use extremum of u(4-6u) on [0,1]: max2/3, min-2.
    zero('Z_weight_four_completed_square',u*(4-6*u)-(R(2,3)-6*(u-R(1,3))**2))
    check('J_coefficient_reserve',3+R(23,2)/1000<R(7,2))
    check('J_radial_derivative_reserve',12+R(60,1000)<13 and 13*R(9,4)<30)
    check('J_projector_derivative_reserve',2*30+4*R(7,2)*R(9,4)<100)
    check('Weyl_amplitude_coefficient_reserve',5+R(28,3)/1000<R(501,100))
    check('Weyl_radial_derivative_reserve',(18+R(265,6)/1000)*R(9,4)<43)
    check('Weyl_projector_derivative_reserve',7*43+48*R(9,4)*R(501,100)<1000)
    check('Schouten_orthogonal_norm_reserve',R(49,4)+R(169,144)<25)
    check('Ricci_orthogonal_norm_reserve',49+R(169,4)<196)
    zero('Hessian_scalar_coefficient_sum',(5*372+11*56)+11*56+2*5*56-3652)
    check('Hessian_scalar_lambda_reserve',R(3652,1000)<4)
    check('Weyl_algebra_contact_reserve',2*40*54+8*40*5+4*40**2==12320<13000)
    check('Ricci_algebra_contact_reserve',7*(2*54+2*14)+R(13,2)*7+2*49<1100)
    check('mixed_divergence_derivative_coefficient',7+R(7,12)*6*40+R(11,6)*7<160)
    check('mixed_scalar_derivative_coefficient',R(7,12)*6*1000+R(11,6)*100<3700)
    check('mixed_zero_order_coefficient',R(7,12)*13000+R(11,6)*1100+R(1001,36)+4<10000)
    zero('fixed_flux_vertex_coefficient',40*R(1,2)-20)
    R4,Z,ep,la,M,C2,J2,box=S.symbols('R Z epsilon lambda M C2 J2 boxZ')
    l0=la+11*ep*la**2/6
    trace=R4+12*l0+ep*(box-R(11,72)*R4**2+R(7,12)*C2+R(11,6)*J2-20*M)
    expected=(1+11*ep*la/3)*Z+ep*box-R(11,72)*ep*Z**2+ep*(R(7,12)*C2+R(11,6)*J2-20*M)
    zero('full_trace_vacuum_subtraction',trace.subs(R4,-12*la+Z)-expected)
    check('trace_quantum_curvature_reserve',8+1+R(2800,3)+R(539,6)<1100)
    gl=R(1,2000**10); flux=R(4,10**27)
    estar=R(1001,250)*R(1,10**12)/R(79*10**8)+1100*gl/180+10*flux/180
    check('actual_trace_defect_majorant',estar<R(6,10**22),estar)
    G=S.symbols('G',positive=True)
    zero('mixed_physical_normalization',S.sqrt(80*S.pi**2)/(8*S.pi*G)-S.sqrt(5)/(2*G))
    check('physical_TT_eta_energy_normalization',4*18*2==144)
    check('A_normalizer_rational_guard',R(5,4)/R(99,100)<R(9,8)**2)
    check('A_sqrt_energy_guards',41600*R(22,7)<362**2 and 176*R(22,7)<24**2 and R(22,7)<R(9,5)**2)
    check('A1_quantum_coefficient_guard',R(9,8)*(160*362+100*24)<70000)
    check('A1_return_coefficient_guard',R(9,8)*12*R(9,5)<25)
    # Equality of a numerical majorant is permitted; the actual inputs are strict.
    check('A1_return_smallness',25*R(6,10**22)/R(3,10**17)<=R(5,10**4))
    check('A1_total_derivative_vertex',70000*R(4,10**17)/180+R(5,10**4)<R(1,1000))
    check('A0_classical_coefficient_guard',R(9,8)*7*24<200)
    check('A0_quantum_coefficient_guard',R(9,8)*(3700*362+10000*24)<2*10**6)
    check('A0_magnetic_coefficient_guard',R(9,8)*20*24<600)
    check('A0_total_scalar_vertex',200+2*10**6*gl/180+600*flux/180<201)
    check('local_W_gradient_coefficient',R(11,9)*R(7,2)+R(4,36*1000)<5)
    check('local_W_scalar_vertex',80*R(22,7)*R(1,10**12)/R(79*10**8)<R(1,10**19))
    check('combined_scalar_vertex_rounding',201**2+R(1,10**19)<1000**2)
    check('combined_derivative_vertex_rounding',5+R(1,1000)**2<9)
    check('positive_absorption_local_block',5*gl/180+R(1,10**12)/(6*R(79*10**8))<R(1,10**21))


def signed_noncommuting_blocks_and_final_budget():
    E,la,T=S.symbols('E lambda T')
    lf,af,j,y,src=S.symbols('L A J y source',commutative=False)
    zero('auxiliary_equation_sign',S.expand(y-((lf-af)*j*y+src)-((1+af*j-lf*j)*y-src)))
    AA=S.Matrix([[3,1],[1,2]])
    q1,q2,q3,q4=S.symbols('q1 q2 q3 q4')
    BB=S.Matrix([[q1,q2],[q3,q4]])
    c1,c2,c3,c4,j1,j2,j3,j4=S.symbols('c1 c2 c3 c4 j1 j2 j3 j4')
    CC=S.Matrix([[c1,c2],[c3,c4]])
    JJ=S.Matrix([[j1,j2],[j3,j4]])
    I=S.eye(2); O=S.zeros(2); inv=(I+AA).inv()
    sch=CC-BB.T*inv*BB
    full=(I+AA).row_join(BB*JJ).col_join(BB.T.row_join(I+CC*JJ))
    lower=I.row_join(O).col_join((BB.T*inv).row_join(I))
    diag=(I+AA).row_join(O).col_join(O.row_join(I+sch*JJ))
    upper=I.row_join(inv*BB*JJ).col_join(O.row_join(I))
    matrix_zero('signed_noncommuting_block_factorization',lower*diag*upper-full)
    matrix_zero('positive_Schur_decomposition',sch-(CC-BB.T*AA.inv()*BB+BB.T*(AA.inv()-inv)*BB))
    a=S.symbols('a',nonnegative=True)
    zero('positive_cross_norm_one_half',(a+1)**2-4*a-(a-1)**2)
    b0=R(1,10**21)
    # sqrt(b0)/2<1/10^10 is a deliberately loose exact rational enclosure.
    check('positive_block_inverse_reserve',(1+R(1,10**10))**2/(1-b0)<2 and b0/4<R(1,10**20))
    first=R(201*10**6,180*4*10**17)
    second=R(1200*3000,180)*R(4,10**17)
    third=R(601*9,4*10**17)
    check('outgoing_zero_derivative_budget',first<R(28,10**13),first)
    check('outgoing_one_derivative_budget',second<=R(8,10**13),second)
    check('outgoing_two_derivative_budget',third<R(14,10**15),third)
    check('outgoing_L_budget_sum',R(28,10**13)+R(8,10**13)+R(14,10**15)<R(4,10**12))
    check('absorbed_outgoing_contraction',2*R(4,10**12)==R(8,10**12)<1)
    check('outgoing_inverse_finite_reserve',2/(1-R(8,10**12))<3)
    zero('physical_b_weight_identification',R(3,2)-1-R(1,2))
    check('mixed_critical_reconstruction_margin',R(1,2)+2==R(5,2)>R(3,2))



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
            liouville_and_bolt()
            outer_and_thermal()
            Bessel_moments_and_resolvents()
            seed_trace_and_vertex_budgets()
            signed_noncommuting_blocks_and_final_budget()
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
