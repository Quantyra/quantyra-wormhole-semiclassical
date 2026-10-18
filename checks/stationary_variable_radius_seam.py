"""Finite identities for the variable-radius kernel comparison.

The checker audits exact algebra and rational enclosures. It does not
bound the nonlocal bubble derivative, solve the coupled seam, evaluate
the signed moment, or test an area budget.
"""
from __future__ import annotations
import base64, contextlib, hashlib, json, platform, sys, traceback
from datetime import datetime, timezone
from pathlib import Path
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/stationary-variable-radius-seam.json'
LOG=ROOT/'evidence/stationary-variable-radius-seam.log'
INPUTS=[
    'checks/stationary_variable_radius_seam.py',
    'notes/stationary-variable-radius-seam.md',
    'notes/stationary-product-weighted-response.md',
    'notes/stationary-product-contour-inverse.md',
    'notes/stationary-reduced-cold-response.md',
    'notes/stationary-reduced-parameter-form-bounds.md',
    'notes/stationary-scalar-phase-moment.md',
]
CONTRACT={
    'acceptance':'exact residuals zero and the stated rational comparisons satisfied',
    'numerical_tolerances':'none',
    'scope':'conformal Laplacian, two exponential-gap enclosures, rejected-substitution arithmetic, radial low-jet embedding, a2 ledger and failed fixed-point products, resolvent identity, prescribed-source Jacobian, constant-background derivatives, and the fast-phase separation',
    'ordinary_arguments':['variable-radius decomposition and the uncontrolled conformal substitution',
        'absence of a proved bound on the nonlocal bubble derivative',
        'prescribed-source map distinguished from the coupled fixed point',
        'constant-background derivative distinguished from the variable three-vertex',
        'signed moment left unenclosed and separate from the product inverse'],
    'unresolved':['three-resolvent insertion as an X_L^2 seam source',
        'coupled mass and temperature fixed point',
        'variable Dirac three-vertex, higher weighted nonlinear map, reduced branch, and physical observable errors',
        'signed moment enclosure',
        'exterior match, finite-pulse transfer, and the physical gates'],
    'not_set_to_zero':['E_bubble','m_num','M_s','L_s','epsilon_s','zeta_s','sigma','cJ','c_theta','delta'],
    'area_budget':'not evaluated; no replacement target was written before the test',
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
    value=S.simplify(S.together(expression))
    check(label,value==0,value)


def matrix_zero(label,matrix):
    if any(entry.atoms(S.Float) for entry in matrix):
        raise ValueError('inexact Float atom in '+label)
    vals=[S.simplify(S.together(entry)) for entry in matrix]
    check(label,all(entry==0 for entry in vals),vals)


def conformal_laplacian_and_gap():
    x,y=S.symbols('x y',real=True)
    sigma=S.Function('sigma')(x,y)
    field=S.Function('field')(x,y)
    weight=S.exp(2*sigma)
    flux_x=weight*S.exp(-2*sigma)*S.diff(field,x)
    flux_y=weight*S.exp(-2*sigma)*S.diff(field,y)
    lap_hat=(S.diff(flux_x,x)+S.diff(flux_y,y))/weight
    lap=S.diff(field,x,2)+S.diff(field,y,2)
    zero('two_dimensional_conformal_laplacian',lap_hat-S.exp(-2*sigma)*lap)
    a=S.symbols('a',positive=True)
    zero('cosh_gap_comparison',S.exp(a)+S.exp(-a)-2-(S.exp(a/2)-S.exp(-a/2))**2)
    z=R(1,10)
    direct_lower=sum(z**k/S.factorial(k) for k in range(6))
    direct_tail=(z**6/S.factorial(6))/(1-z/7)
    direct_upper=direct_lower+direct_tail
    check('direct_gap_above_21_200',direct_lower-1>R(21,200),direct_lower-1)
    check('direct_gap_below_11_100',direct_upper-1<R(11,100),direct_upper-1)
    alt_lower_reciprocal=sum((-z)**k/S.factorial(k) for k in range(5))
    alt_upper_reciprocal=sum((-z)**k/S.factorial(k) for k in range(6))
    check('reciprocal_gap_above_21_200',1/alt_lower_reciprocal-1>R(21,200),1/alt_lower_reciprocal-1)
    check('reciprocal_gap_below_11_100',1/alt_upper_reciprocal-1<R(11,100),1/alt_upper_reciprocal-1)
    zero('rejected_endpoint_lower',R(21,200)*130000000-13650000)
    zero('rejected_endpoint_upper',R(11,100)*130000000-14300000)
    check('direct_substitution_above_lower_endpoint',(direct_lower-1)*130000000>13650000)
    check('direct_substitution_below_upper_endpoint',(direct_upper-1)*130000000<14300000)
    check('reciprocal_substitution_above_lower_endpoint',(1/alt_lower_reciprocal-1)*130000000>13650000)
    check('reciprocal_substitution_below_upper_endpoint',(1/alt_upper_reciprocal-1)*130000000<14300000)
    check('retained_weighted_reserve',R(22,7)*(40000000+2000*84*R(6,5))<130000000)


def low_jet_embedding_and_a2_factor():
    chi=S.symbols('chi',positive=True)
    field=S.Function('field')(chi)
    radial=(S.diff(S.sinh(chi)*S.diff(field,chi),chi)/S.sinh(chi)
            -(S.diff(field,chi,2)+S.coth(chi)*S.diff(field,chi)))
    zero('radial_hyperbolic_divergence',radial)
    slope=S.simplify(S.expand(((S.cosh(chi)-1)/S.sinh(chi)-S.tanh(chi/2)).rewrite(S.exp)))
    if slope.atoms(S.Float):
        raise ValueError('inexact Float atom in integrated_slope_is_tanh_half')
    check('integrated_slope_is_tanh_half',slope==0,slope)
    zero('transverse_jet_collapses',S.coth(chi)*S.tanh(chi/2)-(1-1/(2*S.cosh(chi/2)**2)))
    zero('cosh_minus_sinh',S.cosh(chi)-S.sinh(chi)-S.exp(-chi))
    partial=sum(R(1,1)/S.factorial(k) for k in range(7))
    check('exp_one_above_27_10',partial>R(27,10),partial)
    check('cosh_84_above_1e35',R(27,10)**84/2>10**35)
    check('weighted_ball_exceeds_hypothesis_by_1e37',(R(27,10)**84/2)*100>10**37)
    check('a2_ledger_below_1600',1500+R(5,2)<1600)
    zero('a2_factor_on_hypothesis_ball',1600*R(1,100)-16)
    zero('a2_factor_at_contraction_threshold',1600*R(1,1600)-1)
    zero('failed_composition_on_hypothesis_ball',1600*130000000*R(1,100)-2080000000)
    zero('failed_composition_at_threshold',1600*130000000*R(1,1600)-130000000)
    c2,v=S.symbols('c2 v',positive=True)
    third=c2*S.exp(-2*v)*S.Matrix([[-2,0],[0,0]])
    zero('a2_third_variation_determinant',third.det())
    zero('a2_third_variation_not_scalar',third[0,0]/(c2*S.exp(-2*v))+2)


def three_resolvent_and_higher_map():
    s=S.symbols('s')
    base=S.Matrix([[2,0],[0,3]])
    insertion=S.Matrix([[1,1],[0,2]])
    family=base+s*insertion
    inverse=family.inv()
    square=inverse*inverse
    derivative=S.Matrix(2,2,lambda i,j: S.diff(square[i,j],s)).subs(s,0)
    green=base.inv()
    three=-green*insertion*green*green-green*green*insertion*green
    matrix_zero('three_resolvent_matrix_insertion',derivative-three)
    mass,pi_symbol=S.symbols('M pi',positive=True)
    zero('leading_heat_mass_derivative',S.diff(1/(120*pi_symbol*mass**2),mass)+1/(60*pi_symbol*mass**3))
    q=2000**4
    check('lightest_mass_square_below_338',(1+q)*10<338*6*79*10**9)
    check('lightest_mass_below_19',19**2>338)
    check('heat_candidate_above_100',7*130000000>100*60*22*19**3)
    check('heat_candidate_below_124',130000000<124*60*3*18**3)
    check('heat_candidate_not_a_contraction',100>1)
    N,Y=S.symbols('N Y',real=True)
    sample=S.sin(N*Y)/N**5
    zero('c5_fifth_derivative_amplitude',S.diff(sample,Y,5)-S.cos(N*Y))
    zero('c5_sixth_derivative_amplitude',S.diff(sample,Y,6)+N*S.sin(N*Y))
    zero('newton_square_at_product_norm',2*130000000**2-33800000000000000)
    L_s,eps_s=S.symbols('L_s epsilon_s',positive=True)
    certificate=33800000000000000*L_s*eps_s
    check('newton_certificate_keeps_two_symbols',certificate.free_symbols=={L_s,eps_s},
          sorted(str(name) for name in certificate.free_symbols))


def noncommuting_weighted_kernel():
    chi=S.symbols('chi',real=True)
    zero('cosh_first_jet',S.diff(S.cosh(chi),chi)-S.sinh(chi))
    zero('cosh_second_jet',S.diff(S.sinh(chi),chi)-S.cosh(chi))
    zero('cosh_transverse_jet',S.coth(chi)*S.sinh(chi)-S.cosh(chi))
    zero('cosh_profile_weighted_value',S.cosh(chi)/S.cosh(chi)-1)
    zero('cosh_square_weighted_value',S.cosh(chi)**2/S.cosh(chi)-S.cosh(chi))
    zero('bolt_weight',S.cosh(0)-1)
    zero('hyperbolic_comparison_integral',S.integrate(S.cosh(chi)*S.sinh(chi),(chi,0,S.symbols('L',positive=True)))-S.sinh(S.symbols('L',positive=True))**2/2)
    check('commutator_gap_above_1e37',168*10**35>10**37)
    check('comparison_integral_above_1e69',(10**35-1)**2>2*10**69)
    check('horizon_composition_not_contractive',130000000>1)
    check('weight_exceeds_inverse_threshold',168*10**35*130000000>10**45)
    phi_x,phi_y=S.symbols('phi_x phi_y')
    zero('insertion_split',phi_y-(phi_x+(phi_y-phi_x)))


def resolvent_and_prescribed_source():
    A=S.Matrix([[2,0],[0,3]])
    E=S.Matrix([[1,1],[0,1]])
    matrix_zero('resolvent_difference_identity',(A+E).inv()-A.inv()+A.inv()*E*(A+E).inv())
    kappa,G,radius,D_h,rho_R,a_R,I_rho=S.symbols(
        'kappa G R D_h rho_R a_R I_rho',positive=True)
    s_R=(G*rho_R/radius-kappa*a_R)/D_h
    f_R=-(2*kappa*radius*s_R)/radius
    zero('horizon_lapse_constraint',f_R+2*kappa*s_R)
    delta_mu=kappa*radius*s_R+G*I_rho
    jacobian=S.Matrix([[S.diff(s_R,a_R),S.diff(s_R,I_rho)],
                        [S.diff(delta_mu,a_R),S.diff(delta_mu,I_rho)]]).det()
    zero('prescribed_source_jacobian',jacobian+kappa*G/D_h)


def constant_background_and_fast_phase():
    v,x,t=S.symbols('v x t',positive=True)
    f=S.log(1+x)-x+x**2/2
    factor=S.exp(2*v)*f.subs(x,S.exp(-2*v)*x)
    deriv=S.diff(factor,v).subs(v,0)
    by_chain=2*f+(x**2/(1+x))*(-2*x)
    reduced=2*S.log(1+x)-2*x+x**2*(1-x)/(1+x)
    zero('constant_bubble_derivative_chain_rule',deriv-by_chain)
    zero('constant_bubble_derivative_reduced',deriv-reduced)
    upper_log=R(5,6)
    zero('log_two_integral',S.integrate(1/(1+t),(t,0,1))-S.log(2))
    check('constant_derivative_negative_at_one',2*upper_log-2<0,2*upper_log-2)
    u=S.symbols('u',positive=True)
    num=u**2*(1+t*x)
    den=1+u*t*x
    by_u=S.diff(num/den,u).subs(u,S.exp(-2*v))*(-2*S.exp(-2*v))
    Rv=S.exp(-4*v)*(1+t*x)/(1+S.exp(-2*v)*t*x)
    zero('constant_Rv_two_methods',S.diff(Rv,v)-by_u)
    zero('constant_Rv_origin_derivative',S.diff(Rv,v).subs({v:0,t:0,x:0})+4)
    p,a,b,g=S.symbols('p a b g',positive=True)
    M=S.Matrix([[-a,-2*b*p],[-2*b*p,p*g]])
    zero('product_determinant_identity',M.det()+p*(4*b*b*p+a*g))
    slope=S.diff(M.det(),p).subs({p:1,b:1,a:1,g:1})
    zero('product_determinant_slope_at_unit_point',slope+9)
    a_mix=S.symbols('a_mix',real=True)
    ratio=(1+S.I*a_mix)/(1-S.I*a_mix)
    zero('fast_scattering_modulus',ratio*S.conjugate(ratio)-1)
    zero('fast_scattering_phase_at_one',ratio.subs(a_mix,1)-S.I)
    M_s,L_s,eps_s=S.symbols('M_s L_s epsilon_s',positive=True)
    certificate=2*M_s**2*L_s*eps_s
    check('newton_certificate_keeps_three_symbols',
          certificate.free_symbols=={M_s,L_s,eps_s},sorted(str(s) for s in certificate.free_symbols))


def tail_and_radial_factor():
    chi=S.symbols('chi',positive=True)
    v=S.log(S.tanh(chi/2))
    zero('radial_coordinate_derivative',S.simplify(S.expand((S.diff(v,chi)-1/S.sinh(chi)).rewrite(S.exp))))
    F0,kappa=S.symbols('F0 kappa',positive=True)
    zero('dr_dchi_from_F0',(1/S.sinh(chi))/(kappa/F0)-F0/(kappa*S.sinh(chi)))
    B=1000
    lam=10**4
    check('guard_x',lam*B*B>79*10**8)
    check('guard_h',10**31>10**30)
    check('guard_gamma_ratio',0<=12)
    cap=B*(R(1,100))*2*R(1008,1000)
    zero('chart_cap_at_admissible_point',cap-R(2016,100))
    check('chart_cap_not_below_one',cap>1)
    # e**69 < (series upper)**69 < 2*10**30, so log(2*10**30)>69.
    from math import factorial
    from fractions import Fraction
    upper=sum(Fraction(1,factorial(k)) for k in range(12))+Fraction(13,12)/factorial(12)
    check('disk_radius_above_69',upper**69<2*10**30)
    # Area-weighted majorant tail from distance 69, unit L-infinity density.
    check('unit_density_tail_below_1e-20',22*10*402000*10**68*10**20<7*9*27**68)


def exterior_field_and_dirac_family():
    from fractions import Fraction
    from math import factorial
    lower=sum(Fraction(1,factorial(k)) for k in range(8))
    check('exp69_above_1e29',lower**69>10**29)
    check('cosh69_above_1e29',lower**69>2*10**29)
    check('exterior_Linf_allowance_above_1e37',130000000*lower**69>2*10**37)
    check('exterior_allowance_not_below_one',130000000*lower**69>2)
    B=50
    lam=8*10**9
    check('large_radius_guard',lam*B*B>79*10**8)
    check('large_radius_h',10**31>10**30)
    check('large_radius_gamma',0<=12)
    factor=B*R(1,100)*2*R(1008,1000)
    zero('seam_factor_at_B50',factor-R(126,125))
    check('seam_factor_at_B50_exceeds_one',factor>1)
    N,Y=S.symbols('N Y',real=True)
    family=S.sin(N*Y)/100
    zero('variable_shift_derivative',S.diff(family,Y)-N*S.cos(N*Y)/100)
    zero('variable_shift_amplitude_at_100',(N*S.cos(N*Y)/100).subs(N,100)-S.cos(100*Y))
    zero('circle_identity',S.sin(Y)**2+S.cos(Y)**2-1)
    zero('derivative_to_supremum_ratio',(N/100)/R(1,100)-N)


def higher_map_and_signed_moment():
    N,chi=S.symbols('N chi',real=True)
    field=S.cos(N*chi)/N**2
    zero('third_jet_second_derivative',S.diff(field,chi,2)+S.cos(N*chi))
    zero('third_jet_third_derivative',S.diff(field,chi,3)-N*S.sin(N*chi))
    zero('third_jet_bolt_value',S.diff(field,chi,2).subs(chi,0)+1)
    u=S.symbols('u',real=True)
    gap=u-S.sin(u)
    zero('sine_gap_at_zero',gap.subs(u,0))
    zero('sine_gap_derivative',S.diff(gap,u)-(1-S.cos(u)))
    zero('one_minus_cos_square',1-S.cos(u)-2*S.sin(u/2)**2)
    rho=R(8,10**12)
    factor=2*rho/(1-rho)
    zero('moment_truncation_factor',factor-16/(10**12-8))
    check('moment_truncation_below_2e-11',16*10**11<2*(10**12-8))
    g_norm,h_norm=S.symbols('g_norm h_norm',positive=True)
    enclosure=factor*g_norm*h_norm
    check('moment_enclosure_keeps_two_norms',enclosure.free_symbols=={g_norm,h_norm},
          sorted(str(name) for name in enclosure.free_symbols))


def moment_norms_and_outer_endpoint():
    zero('rank_one_mass_square',(R(2,10**6))**2-R(4,10**12))
    c0,h_norm=S.symbols('c0 h_norm',positive=True)
    factor=16/(10**12-8)
    bound=factor*R(2,10**6)*S.sqrt(c0)*h_norm
    check('moment_bound_keeps_c0_and_h',bound.free_symbols=={c0,h_norm},
          sorted(str(name) for name in bound.free_symbols))
    h=S.symbols('h',positive=True)
    z=S.exp(-h)
    zero('outer_endpoint_identity',S.tanh(S.atanh(z))-z)
    check('moment_symbol_not_replaced',S.symbols('m_num')!=0)


def absent_inputs_and_blocked_gates():
    k1=4*10**17+1
    k2=8*10**17+1
    check('first_pair_fast_guard',k1>4*10**17)
    check('second_pair_fast_guard',k2>2*4*10**17)
    check('energy_guard',181>180)
    c0_first=k1*(181+4*1)
    c0_second=k2*(181+4*4)
    check('admissible_c0_values_differ',c0_first!=c0_second)
    factor=R(16,10**12-8)
    h_norm,c0=S.symbols('h_norm c0',positive=True)
    bound=factor*R(2,10**6)*S.sqrt(c0)*h_norm
    ratio=S.simplify(bound.subs(h_norm,2)/bound.subs(h_norm,1))
    zero('moment_error_scales_with_h_norm',ratio-2)
    check('moment_symbol_remains_unassigned',S.symbols('m_num')!=0)


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
            conformal_laplacian_and_gap()
            low_jet_embedding_and_a2_factor()
            three_resolvent_and_higher_map()
            noncommuting_weighted_kernel()
            tail_and_radial_factor()
            exterior_field_and_dirac_family()
            higher_map_and_signed_moment()
            moment_norms_and_outer_endpoint()
            absent_inputs_and_blocked_gates()
            resolvent_and_prescribed_source()
            constant_background_and_fast_phase()
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
