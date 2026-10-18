"""Finite identities and budgets for the product contour and weighted response.

Proposed destination: checks/stationary_product_weighted_inverse.py.
This checks finite algebra and rational reserves, not the analytic pole sum,
functional calculus, weighted kernel theorem or actual exterior seam inverse.
No numerical bubble, spectral solver or nonlinear certificate is implemented.
"""
from __future__ import annotations
import base64, contextlib, hashlib, json, platform, sys, traceback
from datetime import datetime, timezone
from pathlib import Path
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/stationary-product-weighted-inverse.json'
LOG=ROOT/'evidence/stationary-product-weighted-inverse.log'
INPUTS=[
    'checks/stationary_product_weighted_inverse.py',
    'notes/stationary-product-contour-inverse.md',
    'notes/stationary-product-weighted-response.md',
    'notes/stationary-reduced-cold-response.md',
    'notes/stationary-reduced-parameter-form-bounds.md',
    'notes/stationary-covariant-reduced-inverse.md',
    'notes/stationary-centered-canonical-response.md',
    'notes/stationary-centered-covariant-seed.md',
    'notes/stationary-centered-first-order-geometry.md',
    'notes/stationary-centered-exterior-source.md',
    'notes/stationary-centered-complete-residual.md',
    'notes/stationary-reference-cold-bolt.md',
    'notes/stationary-seed-quantum-form.md',
    'notes/stationary-curved-quantum-tt.md',
    'notes/stationary-weighted-tt-inverse.md',
    'notes/stationary-disk-and-normal-inverse.md',
    'notes/stationary-scalar-outgoing-inverse.md',
    'notes/stationary-massive-boundary-layer.md',
    'notes/stationary-paired-dirac-phase.md',
    'notes/stationary-exterior-conformal-domain.md',
]
CONTRACT={
    'acceptance':'exact residuals zero and stated rational comparisons satisfied; a strict source premise may give equality in its numerical majorant',
    'numerical_tolerances':'none',
    'scope':'finite primary-pole normalization/recursion, moment identities, exact restored coefficient, Stieltjes/sectors, coupled inverse, cluster/Dunford algebra, Q identities and weighted geometric budgets',
    'ordinary_arguments':['rank-two massive Green and physical pairing identification with the primary bubble',
        'primary spectral representation, common heat subtraction and locally uniform angular convergence',
        'gamma convexity, total-variation quadrature estimate and continuous contour/pole exclusion',
        'Mehler-Fock inversion, complete-H2 spectral domains and operator-norm Dunford reconstruction',
        'complex-Q oscillatory integral estimate, kernel convergence, Schur and geometric group commutation',
        'bolt-regular C2 source completion and complete-product tail realization'],
    'unresolved':['actual variable-radius quantum operator comparison in the weighted source/response norms',
        'finite-neck nonlocal exterior return and physical mass/temperature matching',
        'nonlinear reduced source-to-observable bound, original signed source condition and all physical gates'],
    'matrix_scope':'finite matrix and scalar identities do not prove infinite-dimensional domains or a nonlocal finite-neck boundary-value inverse',
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


def primary_residue_and_contact():
    m,j,p,tau=S.symbols('M j p tau',positive=True)
    de=m+R(1,2)
    # Primary (D.10), using Gamma(z+1)=z Gamma(z), without evaluating Gamma.
    c_ratio=((4*de+2*j)*(4*de+2*j+1)*(2*de+1+2*j)*(2*de+2+2*j)
             /(16*(2*j+2)*(2*j+3)*(2*de+R(1,2)+2*j)*(2*de+R(3,2)+2*j)))
    primary_ratio=((2*m+j+R(3,2))**2*(m+j+1)**4
       /((j+1)**2*(2*m+2*j+R(3,2))*(2*m+2*j+R(5,2))
                       *(2*m+2*j+2)*(2*m+2*j+3)*c_ratio))
    pole_ratio=(2*m+2*j+R(7,2))/(2*m+2*j+R(3,2))
    target_ratio=(pole_ratio*(j+R(3,2))/(j+1)*(2*m+j+R(3,2))/(2*m+j+1)
                  *((m+j+1)/(m+j+R(3,2)))**2)
    zero('primary_residue_gamma_recursion',primary_ratio*pole_ratio-target_ratio)
    # At j=0 use Gamma(M+3/2)=(M+1/2)Gamma(M+1/2), and
    # Gamma(2M+2)=(2M+1)Gamma(2M+1). W(1)=sqrt(pi)/2.
    ga,gb,gc,gd=S.symbols('Gamma_Mhalf Gamma_Mone Gamma_2Mone Gamma_2Mthreehalf')
    primary0=gd*gb**2/(S.sqrt(S.pi)*ga**2*(2*m+1)*gc*(m+R(1,2)))
    reduced0=(1/S.pi)*(S.sqrt(S.pi)/2)*(gd/gc)*gb**2/((m+R(1,2))*ga)**2
    zero('primary_rank_two_initial_coefficient',primary0-reduced0)
    zero('positive_pole_shift',(2*m+2*j+R(3,2))**2-R(9,4)-(2*m+2*j)*(2*m+2*j+3))
    zero('pole_reference_subtraction',-1/(tau+p)+1/tau-p/(tau*(tau+p)))
    zero('exact_cubic_pole_split',p/(tau*(tau+p))-(p/tau**2-p**2/tau**3+p**3/(tau**3*(tau+p))))
    check('first_massive_pole_support',36*39==1404)
    zero('physical_rank_fold_factor',4*R(1,2)-2)
    G,A,B=S.symbols('G A B',positive=True)
    gamma0=1/(8*S.pi*G*A**2)
    zero('physical_product_pairing_factor',(1/(4*S.pi*B**2*A**2))/gamma0-2*G/B**2)
    # The common a2 polynomial is only the leading quadratic coefficient.
    lam=S.symbols('lambda_H2')
    poly=8*lam**2-32*lam+16
    zero('stationary_common_a2_polynomial',poly.subs(lam,2-p)-poly.subs(lam,2)-8*p**2)
    zero('reduction_p2_coefficient',R(32,5760)-R(1,180))


def moments_and_restored_quadratic():
    x,m,z,n,q,C=S.symbols('x M z n q C',positive=True)
    zero('gamma_ratio_lower_square_reserve',x*x/(x+R(1,2))-(x-R(1,2))-R(1,4)/(x+R(1,2)))
    zero('moment_integral',S.integrate(z**3*S.sqrt(1-z*z),(z,0,1))-R(2,15))
    f2=(x-m*m)/x**6
    zero('moment_maximum_stationary_point',S.diff(f2,x).subs(x,R(6,5)*m*m))
    zero('moment_maximum_square',f2.subs(x,R(6,5)*m*m)-R(125,216)**2/(5*m**10))
    check('moment_total_variation_relative_39',(R(3750,432)**2)/5<R(39,10)**2)
    check('gamma_residue_coarse_2',R(4,3)*(1+R(1,4)/R(37,2))<2)
    check('tau_upper_correction_053',R(1,2)+R(1,2)/18<R(53,100))
    poly=S.expand(((1+R(18,10)*z)*(1-R(53,100)*z)**3-1)/z)
    # Positive constant plus worst negative coefficients at z=1/18.
    bound=poly.coeff(z,0)+sum(min(poly.coeff(z,k),0)*R(1,18)**k for k in range(1,4))
    check('reciprocal_cube_upper_polynomial',bound>R(9,100))
    bern=S.expand(((1-R(3,2)*z)*(1+z/2)**3-1)/z**2)
    check('reciprocal_cube_lower_Bernoulli',all(coef<=0 for coef in S.Poly(bern,z).all_coeffs()))
    check('upper_moment_factor_21',R(25,100)+R(18,10)+R(45,100)/18<R(21,10))
    zero('lower_moment_factor_175',(1-R(3,2)*z)*(1-z/4)-(1-R(7,4)*z)-R(3,8)*z*z)
    check('quadrature_upper_relative_65',R(21,10)+R(39,10)+R(21,10)*R(39,10)/18<R(65,10))
    zero('quadrature_lower_relative_565',R(7,4)+R(39,10)-R(565,100))
    check('quadratic_curvature_reserve_7',R(65,10)<7 and R(565,100)<7)
    zero('third_moment_integral_test',R(1,18)+R(1,6)-R(2,9))
    zero('third_moment_prefactor',R(1,128)*R(19,18)*R(2,9)-R(38,20736))
    check('third_moment_per_mode_500',R(38,20736)<R(1,500))
    zero('angular_fourth_moment_partial_fraction',2*(q+2*n)*C**2/(n*n*(n+q)**2)-(2*C**2/q)*(1/n**2-1/(n+q)**2))
    zero('angular_fourth_moment_reserve',4*R(352,100000)-R(1408,100000))
    check('summed_third_pole_moment',R(1,500)*R(1408,100000)<R(3,100000))
    check('elementary_log_q_lower',(R(11,4))**7<2000)
    check('F_eff_below_minus_51',4-56+R(4,10**21)<-51)
    zero('finite_quadratic_correction',R(7,120)*R(2,3)-R(7,180))
    zero('restored_quadratic_lower_numerator',R(51,120)-R(1,180)-R(7,180)-R(137,360))
    check('restored_quadratic_lower_C_over_10',R(137,360)/R(22,7)>R(1,10))
    check('restored_quadratic_upper_C_over_3',(R(100,120)+R(7,180))/3<R(1,3))
    check('physical_slope_alpha',2*10**24*R(85,10**34)<=R(17,10**9))
    zero('physical_beta_lower',2*R(1,10)-R(1,5))
    zero('physical_beta_upper',2*R(1,3)-R(2,3))
    zero('physical_positive_measure_moment',2*R(3,100000)-R(6,100000))


def Stieltjes_and_contour_sectors():
    x,y,t,r=S.symbols('x y tau r',positive=True)
    re=(x*(t+x)+y*y)/((t+x)**2+y*y)
    zero('Stieltjes_real_part_comparison',re-x/(t+x)-t*y*y/((t+x)*((t+x)**2+y*y)))
    zero('Stieltjes_sector_imaginary_part',S.im((x+S.I*y)/(t+x+S.I*y))-t*y/((t+x)**2+y*y))
    zero('Stieltjes_real_sublinearity',(r/x)*(x/(t+x))-r/(t+r)-r*(r-x)/((t+x)*(t+r)))
    xx=S.symbols('xi',real=True)
    w=S.symbols('w',real=True); delta=w*(w+1)-2; kap=2*w+1
    pp=xx**2-delta-S.I*kap*xx
    zero('parabolic_spectral_change',pp-(2-(w+S.I*xx)*(w+1+S.I*xx)))
    zero('parabola_absolute_square',S.expand(pp*S.conjugate(pp))-(xx**4+(kap**2-2*delta)*xx**2+delta**2))
    check('two_parabola_vertices',delta.subs(w,R(101,100))==R(301,10000) and delta.subs(w,R(99,100))==-R(299,10000))
    check('parabola_minimum_positive',all((kap**2-2*delta).subs(w,v)>0 for v in (R(99,100),R(101,100))))
    check('parabola_distance_from_zero',R(299,10000)>R(1,34))
    check('compact_contour_radius_21',(16+R(301,10000))**2+(4*R(302,100))**2<21**2)
    eg=R(85,10**34); al=R(17,10**9)
    hh=R(2,3)+21/(1-R(21,1404))*R(6,100000)
    check('compact_H_below_etaG',hh<1)
    check('compact_a_below_137',R(61,10)*(21+R(14,10))<137)
    check('compact_alpha_division_5e_minus6',R(61,10)*(1+R(14,10)/R(299,10000))*al<R(5,10**6))
    check('compact_D_margin_399',4*R(999,1000)**2-137*eg-R(5,10**6)>R(399,100))
    check('compact_full_D_majorant',3+R(3,4)*21*eg<R(399,100))
    xlo=16-R(301,10000); ymax=4*R(302,100)
    check('high_frequency_angle_p',ymax/xlo<R(76,100))
    check('high_frequency_angle_a',ymax/(xlo-R(14,10))<R(83,100))
    # For x>=xlo, kappa^2(x+.0301)/x^2 is decreasing.
    zero('right_region_angle_monotonicity',S.diff((x+R(301,10000))/x**2,x)+(x+R(602,10000))/x**3)
    cs=1-R(76,100)*R(83,100)
    check('sector_cosine_one_fifth',cs>0 and 25*cs**2>(1+R(76,100)**2)*(1+R(83,100)**2))
    check('high_a_lower_5r',R(59,10)*(1-R(14,10)/xlo)>5)
    check('high_a_upper_7r',R(61,10)*(1+R(14,10)/xlo)<7)
    check('high_D_positive_margin',4*R(999,1000)**2-7*al>3)
    check('modulus_to_real_ratio_126',1+R(76,100)**2<R(126,100)**2)
    check('real_H_fraction_three_quarters',R(100,126)>R(3,4))
    check('compact_H_modulus_ratio_102',(1404+21)/(R(1404)-R(301,10000))<R(102,100))
    check('global_H_modulus_factor_2',R(126,100)<2 and R(102,100)<2)


def inverse_cluster_and_global_poles():
    a,b,p,g,al,h,FS,FT=S.symbols('a b p g alpha H F_S F_T')
    mat=S.Matrix([[-a,-2*b*p],[-2*b*p,p*g]])
    den=4*b*b*p+a*g
    inv=S.Matrix([[-g/den,-2*b/den],[-2*b/den,a/(p*den)]])
    matrix_zero('coupled_inverse_exact',mat*inv-S.eye(2))
    zero('coupled_determinant_exact',mat.det()+p*den)
    zero('positive_denominator_decomposition',den.subs(g,al+p*h)/p-(4*b*b+a*h+a*al/p))
    check('inverse_SS_constant',R(2,1)/R(3,4)<3)
    check('inverse_ST_constant',(2*R(1001,1000))/R(3,4)<3)
    check('inverse_TT_affine_constants',R(61,10)/R(3,4)<12 and R(61,10)*R(14,10)/R(3,4)<12)
    zero('inverse_compact_TT_entry',12*(1+R(1,34))*34**2-14280)
    check('inverse_uniform_matrix_16000',14280+102<16000)
    check('inverse_high_row_16',3+12*(1+1/(16-R(301,10000)))<16)
    check('inverse_high_frequency_17',16/(1-R(301,10000)/16)<17)
    # Physical real spectrum is a separate estimate, not a contour substitution.
    check('physical_spectrum_a_positive',R(59,10)*(1-R(14,10)/R(9,4))>2)
    check('physical_spectrum_D_margin',4*R(999,1000)**2-R(61,10)*R(17,10**9)>3)
    check('physical_spectrum_P2_graph_bound',3+12*(1+R(4,9))<24)
    circle=R(1,1000); aa=R(61,10)*(circle+R(14,10)); eg=R(85,10**34); alpha=R(17,10**9)
    check('cluster_circle_a_8547',aa<R(8547,1000))
    check('right_region_compact_height',R(302,100)**2*((16-R(301,10000))+R(301,10000))==R(1208,100)**2)
    check('right_region_compact_radius',(16-R(301,10000))**2+R(1208,100)**2<21**2)
    loss=R(61,10)*(1+R(14,10)/circle)*alpha
    check('outside_cluster_alpha_loss',loss<R(15,100000))
    check('outside_cluster_no_zero_budget',4*R(999,1000)**2-loss-137*eg>R(399,100))
    check('cluster_Rouche_budget',aa*(alpha+circle*eg)<4*R(999,1000)**2*circle)
    check('cluster_inverse_row_bound',aa/(R(399,100)*circle**2)+2*R(1001,1000)/(R(399,100)*circle)<2200000)
    check('cluster_nu_square_root',R(2999,1000)**2<9-4*circle)
    check('cluster_nu_displacement',2*circle/(3+R(2999,1000))<R(334,10**6))
    check('cluster_nu_real_part',1-R(334,10**6)>R(999,1000))


def Q_transform_and_Dunford_algebra():
    q=S.symbols('q',positive=True); nu=S.symbols('nu',real=True)
    ff=S.Function('F')(q); aa=(nu+1)/2
    yy=q**aa*ff
    def dd(v): return -2*q*S.diff(v,q)
    radial=dd(dd(yy))+(1+q)/(1-q)*dd(yy)-nu*(nu+1)*yy
    hyper=q*(1-q)*S.diff(ff,q,2)+(nu+R(3,2)-(nu+R(5,2))*q)*S.diff(ff,q)-(nu+1)*ff/2
    reduced=S.simplify(radial/q**aa)*(1-q)/(4*q)
    zero('Q_Euler_hypergeometric_radial_equation',reduced-hyper)
    zero('Q_leading_nu0_beta',S.gamma(1)*S.gamma(R(1,2))/S.gamma(R(3,2))-2)
    zero('Q_leading_nu1_beta',S.gamma(2)*S.gamma(R(1,2))/S.gamma(R(5,2))-R(4,3))
    z=S.symbols('z',positive=True)
    def dz(v): return -z*S.diff(v,z)
    coth=(1+z*z)/(1-z*z); cosh=(1/z+z)/2; sinh=(1/z-z)/2
    q0=S.log((1+z)/(1-z)); q1=cosh*q0-1
    zero('Q0_radial_equation',dz(dz(q0))+coth*dz(q0))
    zero('Q0_delta_flux',sinh*dz(q0)+1)
    zero('Q1_radial_equation',dz(dz(q1))+coth*dz(q1)-2*q1)
    zero('Q1_leading_coefficient',S.limit(q1/z**2,z,0,dir='+')-R(4,3))
    check('Q_high_oscillatory_constant',32<36)
    check('Q_unified_frequency_constant',72<81 and 8<81)
    t,la,a,alpha,b=S.symbols('t lambda a alpha b',nonzero=True)
    basic=1/(t-a)
    zero('Dunford_spectral_residue_sign',S.residue(basic/(la-t),t,la)+1/(la-a))
    zero('Dunford_crossed_simple_pole',S.residue(basic/(la-t),t,a)-1/(la-a))
    zero('Dunford_coalescing_double_pole',S.residue(1/(t*t*(la-t)),t,0)-1/la**2)
    pstar=-a*alpha/(4*b*b)
    cluster=(1/alpha)*(1/la-1/(la-pstar))
    exact=a/(la*(4*b*b*la+a*alpha))
    zero('cluster_cancelled_inverse_without_separate_residue_bounds',cluster-exact)
    zero('cluster_coalescence_limit',S.limit(exact,alpha,0)-a/(4*b*b*la**2))
    zero('cluster_Ward_residue',S.residue(a/(t*(4*b*b*t+a*alpha)),t,0)-1/alpha)
    zero('cluster_physical_residue',S.residue(a/(t*(4*b*b*t+a*alpha)),t,pstar)+1/alpha)
    n=S.symbols('nu_growth')
    zero('cluster_nu_branch_identity',n-1+2*(2-n*(n+1))/(3+(2*n+1)))


def spatial_kernel_and_weighted_geometry():
    check('contour_speed_high',2+R(302,100)/4<=R(2755,1000))
    x=S.symbols('x',positive=True)
    zero('high_frequency_integral',S.integrate(x**(-R(3,2)),(x,4,S.oo))-1)
    integral=2*16000*4*R(1102,100)+2*17*R(2755,1000)
    check('contour_absolute_integral_15e6',integral<1500000)
    check('contour_kernel_constant_4e5',9*1500000/(4*R(3)**2)<400000)
    check('cluster_kernel_constant_1e3',R(1,1000)*2200000*2/(2*R(3))<1000)
    z=S.symbols('z',positive=True)
    sinh=(1/z-z)/2
    zero('weighted_H2_area_cancellation',2*sinh*(1/z)*z*z/S.sqrt(1-z*z)-S.sqrt(1-z*z))
    chi,th=S.symbols('chi theta',real=True)
    zero('cosh_weight_log_derivative',S.diff(S.log(S.cosh(chi)),chi)-S.tanh(chi))
    zero('cosh_weight_Lipschitz_square',1-S.tanh(chi)**2-1/S.cosh(chi)**2)
    tt=R(168,1000)
    exp_upper=1+tt+tt**2/(2*(1-tt/3))
    check('finite_neck_exponential_12',exp_upper<R(6,5))
    check('weighted_row_and_column_13e8',R(22,7)*(40000000+2000*84*R(6,5))<130000000)
    # These are finite coordinate identities; smoothness of boosts at the
    # bolt and commutation with the complete spectral inverse are ordinary.
    v=S.Function('V')(chi)
    def Z1(f): return S.cos(th)*S.diff(f,chi)-S.coth(chi)*S.sin(th)*S.diff(f,th)
    def Z2(f): return S.sin(th)*S.diff(f,chi)+S.coth(chi)*S.cos(th)*S.diff(f,th)
    zero('Killing_first_radial_jet',Z1(v)-S.cos(th)*S.diff(v,chi))
    zero('Killing_second_radial_jet',Z1(Z1(v))-(S.cos(th)**2*S.diff(v,chi,2)+S.sin(th)**2*S.coth(chi)*S.diff(v,chi)))
    zero('Killing_transverse_radial_jet',Z2(Z2(v))-(S.sin(th)**2*S.diff(v,chi,2)+S.cos(th)**2*S.coth(chi)*S.diff(v,chi)))
    zero('Killing_mixed_radial_jet',Z1(Z2(v))-S.sin(th)*S.cos(th)*(S.diff(v,chi,2)-S.coth(chi)*S.diff(v,chi)))
    zero('Killing_reversed_radial_jet',Z2(Z1(v))-Z1(Z2(v)))
    zero('Killing_at_axis_gradient',Z1(v).subs(th,0)-S.diff(v,chi))
    zero('Killing_at_axis_radial_Hessian',Z1(Z1(v)).subs(th,0)-S.diff(v,chi,2))
    zero('Killing_at_axis_transverse_Hessian',Z2(Z2(v)).subs(th,0)-S.coth(chi)*S.diff(v,chi))
    zero('Killing_squared_coefficients_partition',S.cos(th)**2+S.sin(th)**2-1)
    zero('Killing_cross_coefficient_bound_identity',1-4*S.sin(th)**2*S.cos(th)**2-S.cos(2*th)**2)
    # A regular even source has a finite transverse bolt jet, not a free log.
    v0,v2,v4=S.symbols('v0 v2 v4')
    even=v0+v2*chi**2+v4*chi**4
    zero('regular_even_bolt_transverse_jet',S.limit(S.coth(chi)*S.diff(even,chi),chi,0)-2*v2)


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
            primary_residue_and_contact()
            moments_and_restored_quadratic()
            Stieltjes_and_contour_sectors()
            inverse_cluster_and_global_poles()
            Q_transform_and_Dunford_algebra()
            spatial_kernel_and_weighted_geometry()
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
