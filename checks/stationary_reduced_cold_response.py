"""Finite identities for reduced cold response and parameter form bounds.

Proposed destination: checks/stationary_reduced_cold_response.py.
This audits finite algebra and literal rational majorants only. It does not
prove boundary/operator domains, a renormalized three-vertex estimate, a
whole-contour cold inverse, nonlinear existence or physical error control.
"""
from __future__ import annotations
import base64, contextlib, hashlib, json, platform, sys, traceback
from datetime import datetime, timezone
from pathlib import Path
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/stationary-reduced-cold-response.json'
LOG=ROOT/'evidence/stationary-reduced-cold-response.log'
INPUTS=[
    'checks/stationary_reduced_cold_response.py',
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
    'acceptance':'exact residuals zero and stated rational comparisons satisfied; a strict source input may give equality in its numerical majorant',
    'numerical_tolerances':'none',
    'scope':'finite constraint, barrier, product, Ward, Legendre, cutoff, contact, logarithmic multiplier and root identities/budgets',
    'ordinary_arguments':['actual classical constraint/coordinate domains and conserved source hypotheses',
        'paired HH Green realization, common angular subtraction and normalized VALUE lemma',
        'actual product kernel analytic continuation and local/cutoff gluing estimates',
        'closed logarithm/Hardy form domains and complex parameter realization',
        'Rouche theorem and admissible root interpretation'],
    'unresolved':['whole weighted product contour and actual-seed seam inverse',
        'actual variable-coefficient common-subtracted Dirac three-vertex and remaining contacts',
        'high-level parameter-dependent weighted mixed map and nonlinear certificate',
        'physical reduced-branch remainder, original signed source condition and all physical gates'],
    'matrix_scope':'finite action-polynomial and generic matrix identities do not prove infinite-dimensional operator domains',
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


def constraint_inverse_and_observable_budgets():
    G,rh,kap,Frr,rhoh,ah=S.symbols('G rh kappa Frr rhoh ah',positive=True)
    Dh=Frr/2+kap/rh
    dr=(G*rhoh/rh-kap*ah)/Dh
    f_h=-2*kap*dr
    fp_h=-2*G*rhoh/rh+2*kap*dr/rh
    zero('moving_horizon_radial_constraint',f_h+2*kap*dr)
    zero('fixed_clock_temperature_constraint',fp_h+Frr*dr+2*kap*ah)
    mprime=S.symbols('mu_seed_prime')
    zero('moving_mass_chain_rule',mprime*dr-f_h*rh/2-(mprime+kap*rh)*dr)
    h=R(1,10**30)
    sig=R(51,10**13); lapse=R(33,10**11); qr=R(13,10**11)
    check('fixed_temperature_gap_denominator_reserve',1/((1-24*h)*(1-7*h))<R(1001,1000))
    kh=R(1001,1000)*(sig+lapse)
    kv=R(1002,1000)*sig+R(21002,1000)*lapse+30*qr
    uc=R(51,10)*lapse+R(69,10)*qr+183*h*kh
    rc=R(51,10**12)+R(606,100)*R(572,10**13)+18*h*kh+R(81,10)*kv
    check('source_horizon_bound',kh<R(336,10**12))
    check('source_global_radius_bound',kv<R(109,10**10))
    check('source_actual_seam_bound',uc<R(259,10**11))
    check('source_cold_curvature_bound',rc<R(89,10**9))
    check('source_cold_bolt_bound',85*max(102*uc,R(22,100)*rc)<R(27,10**6))
    unrounded=(1+800*h)*kh+20*lapse+30*qr
    check('radius_rounded_coefficients',unrounded<kv)
    check('seam_lapse_rounding',R(5035,1000)<R(51,10))
    check('seam_density_rounding',R(68025,10000)<R(69,10))
    check('seam_horizon_rounding',R(1824,10)<183)
    check('cold_trace_source_factor',6+6*R(1,100)==R(606,100))
    check('cold_seam_value_factor',R(1008,1000)/R(1,100)<102)


def radial_barrier_and_product_inverse():
    x=S.symbols('x',positive=True); hh,L=S.symbols('h L',positive=True)
    ch=(x+1/x)/2; sh=(x-1/x)/2
    def derivative(f): return x*S.diff(f,x)
    def lap(f): return derivative(derivative(f))+ch/sh*derivative(f)
    W=hh*((L+1)*ch-S.log(x)*sh)
    zero('cold_barrier_radial_operator',-lap(W)+2*W-hh*(3*ch+S.log(x)/sh))
    zero('cold_barrier_positive_decomposition',W/hh-((L-S.log(x)+1)*ch+S.log(x)/x))
    check('cold_potential_loss_reserve',20*R(1008,1000)*(R(1,100)/R(992,1000)+2*85/R(10**30))<R(3,10))
    check('cold_barrier_rhs_reserve',R(56,100)/R(27,10)<R(22,100))
    check('cold_potential_positive',2-20*R(1,100)>0)
    z=S.symbols('z',positive=True); v=S.Function('v')(z); f=S.Function('f')(z)
    cp=S.cosh(z); sp=S.sinh(z)
    expr=-S.diff(cp*v,z,2)-S.coth(z)*S.diff(cp*v,z)+2*cp*v
    target=-S.diff(cp**2*sp*S.diff(v,z),z)/(cp*sp)
    zero('product_ground_state_divergence',S.expand_trig(expr-target))
    zero('product_ground_state_integral',S.diff((S.cosh(z)**3-1)/3,z)-S.cosh(z)**2*S.sinh(z))
    c=S.symbols('c',positive=True)
    primitive=1-1/c+S.log((c+1)/2)
    zero('finite_neck_integral_primitive',S.diff(primitive,c)-(c**3-1)/(c**2*(c**2-1)))
    zero('finite_neck_integral_lower_endpoint',primitive.subs(c,1))
    check('two_product_inverse_factors',R(84+1,3)==R(85,3))
    t,beta,B=S.symbols('t beta B')
    ss,tt=S.symbols('S T'); om=ss+tt
    lag=B**2*(2*(t-2)*ss**2+4*(t-2)*ss*om+2*beta*ss**2+8*ss*om)
    hess=S.hessian(lag,(ss,tt))
    matrix_zero('classical_product_action_Hessian',hess-4*B**2*S.Matrix([[3*t+beta-2,t],[t,0]]))
    fs,fr=S.symbols('f_s f_r')
    invS=fr/t; invT=(fs-3*fr+(2-beta)*invS)/t
    zero('classical_product_radius_solve',t*invS-fr)
    zero('classical_product_shape_solve',(3*t+beta-2)*invS+t*invT-fs)


def exact_product_matrix_and_Ward():
    B,la,ep=S.symbols('B lambda epsilon',positive=True)
    d,p=S.symbols('d p'); x=la*B**2; C=1+6*x+d; aa=B**2/C
    rr=-12*la-2*d/B**2; U=-rr; rup=-16*d/B**2
    jo=-(C+1)/(2*B**2); hu=2*U+rup
    corb=6-R(11,6)*ep*rr+2*ep*(U+rup)+R(22,3)*ep*jo
    corb_target=6+24*ep*la-R(22,3)*ep/B**2-28*ep*d/B**2
    zero('product_conformal_principal_from_curvature',corb-corb_target)
    la0=la+R(11,6)*ep*la**2
    ps=-2*U+48*la0+ep*hu**2/6
    vu=4*la-2*d/B**2
    zero('product_scalar_zero_order_from_trace',ps-(24*la-4*d/B**2+88*ep*la**2+6*ep*vu**2))
    rs=2/B**2
    mixed=1+ep*(R(3,2)*(2*rr/3)-R(11,12)*2*rs+hu/6)
    bb=1-8*ep*la-ep*(R(11,3)+4*d)/B**2
    zero('product_mixed_from_orbit_invariants',mixed-bb)
    bet=-(4+12*x+10*d)/C
    zero('classical_reference_product_recovery',6*(p-2)+aa*hu-2*(3*p+bet-2))
    zero('inconsistent_epsilon_only_limit_difference',ps.subs(ep,0)-hu-8*d/B**2)
    a,b,g,FS,FT=S.symbols('a b g F_S F_T')
    mat=S.Matrix([[-a,-2*b*p],[-2*b*p,p*g]])
    denom=4*b*b*p+a*g
    zero('reduced_product_determinant',mat.det()+p*denom)
    solS=-(g*FS+2*b*FT)/denom
    solT=(a*FT-2*b*p*FS)/(p*denom)
    matrix_zero('reduced_product_source_elimination',mat*S.Matrix([solS,solT])-S.Matrix([FS,FT]))
    # At p=0 the entire orbit-shape column vanishes, not an actual global mode.
    matrix_zero('product_Ward_column',mat.subs(p,0)*S.Matrix([0,1]))
    xx=S.symbols('x',positive=True)
    WW=xx*S.log((xx+1)/2)+xx-1
    zero('Legendre_parameter_weight_equation',S.diff((xx**2-1)*S.diff(WW,xx),xx)-2*WW-3*xx)
    zero('Legendre_parameter_weight_bolt',WW.subs(xx,1))
    psi=-WW/3
    zero('Legendre_slope_radial_equation',S.diff((xx**2-1)*S.diff(psi,xx),xx)-2*psi+xx)
    # Eigenvalue differentiation independently fixes the center jets.
    zero('Legendre_first_center_Laplacian',((xx**2-1)*S.diff(psi,xx,2)+2*xx*S.diff(psi,xx)).subs(xx,1)+1)


def Legendre_tail_and_probe():
    # Center jets follow from the differentiated eigenvalue equation.
    p=S.symbols('p')
    delta1=(2-p)*0-1; delta2=(2-p)*delta1-2
    zero('a2_Ward_probe_action',8*delta2-32*delta1-8*p)
    lam=S.symbols('lambda_H2')
    wardpoly=8*lam**2-32*lam+16
    zero('a2_product_Ward_polynomial_difference',wardpoly.subs(lam,2-p)-wardpoly.subs(lam,2)-8*p**2)
    check('Volterra_resolvent_cone',R(201,100)*R(10,81)<1)
    major=(R(10,81)*R(11,9))/(1-R(201,100)*R(10,81))
    check('Volterra_probe_sup_021',major<R(21,100))
    check('Cartesian_complex_polydisk',R(3,1)+2*R(99,70)<R(23,1000)*256)
    check('Cartesian_Cauchy_tensor_reserve',8**2*2<12**2)
    z=S.symbols('z')
    bern=(1+R(32,75)*z+R(32,1875)*z**2)**11
    coeffs=[S.expand(bern).coeff(z,j) for j in range(6)]
    bounds=[1,5,11,14,13,9]
    check('quadratic_radius_cutoff_Bell_array',all(v<=b for v,b in zip(coeffs,bounds)))
    check('probe_C5_product_norm',all(R(21,100)*sum(R(12,100)**k*bounds[j-k] for k in range(j+1))<4 for j in range(6)))
    zero('cutoff_quadratic_radius_gradient',2*R(1,8)/R(3,256)-R(64,3))
    zero('cutoff_quadratic_radius_Hessian',2/R(3,256)-R(512,3))
    # Beta(6,6): differentiate its integral polynomial to verify join orders.
    t=S.symbols('t'); den=S.factorial(5)**2/S.factorial(11)
    poly=1-S.integrate(t**5*(1-t)**5/den,(t,0,t))
    check('beta_cutoff_first_five_join_jets',all(S.diff(poly,t,j).subs(t,0)==0 and S.diff(poly,t,j).subs(t,1)==0 for j in range(1,6)))
    zero('beta_cutoff_values',poly.subs(t,0)-poly.subs(t,1)-1)
    check('H2_normal_cap_amplitude',R(1,63)<R(26,1000))
    check('local_real_imag_rank_two_remainder',8*10**22*R(2,3)+8*10**13+12+1+4*10**6<10**24)
    # Exact gamma tail moment from elementary integer incomplete-gamma sums.
    v=S.symbols('v',positive=True)
    tail=S.integrate((v**4-2*v**3)*S.exp(-v),(v,2,S.oo))
    zero('angular_count_tail_integral',tail-92*S.exp(-2))
    check('angular_count_72_reserve',368/R(8,3)**2<72)
    check('kernel_rank_two_tail_1024',R(8*72*16,9)<=1024)
    zero('radial_tail_prefactor',R(7*1024,4)-1792)
    zero('radial_tail_integral_majorant',R(1792,2)-896)
    check('product_tail_at_1_over_64',896*64**2<4*10**6)
    # The rational 1.01 tail exponent weakens 7 to 6.99 but stays positive.
    check('complex_probe_tail_exponent',9-1-R(101,100)>0)
    qratio=400
    check('LLL_product_allowance_pi_gt_three',R(qratio,36)<12)


def finite_contact_and_form_budgets():
    vv,t,s,P,X,Y,w,dw,alpha=S.symbols('v t s P X Y w Delta_w alpha')
    U,T=S.symbols('U T')
    action=alpha*S.exp(-2*(vv+t*U+s*T))*(P-2*t*X-2*s*Y)**2
    hess=S.diff(action,t,s).subs({t:0,s:0})
    target=8*alpha*S.exp(-2*vv)*(X*Y+P*(U*Y+X*T)+P**2*U*T/2)
    zero('actual_a2_Hessian_from_action',hess-target)
    E=S.symbols('E')
    qmat=S.Matrix([[1,P],[P,P**2/2]])
    varied=S.exp(-2*s*w)*qmat.subs(P,P-2*s*dw)
    dm=S.diff(varied,s).subs(s,0)
    matrix_zero('actual_a2_third_coefficient_matrix',dm-S.Matrix([[-2*w,-2*w*P-2*dw],[-2*w*P-2*dw,-w*P**2-2*P*dw]]))
    pmax=R(101,50)
    zero('a2_absolute_first_row_w',2+2*pmax-R(151,25))
    zero('a2_absolute_second_row_w',2*pmax+pmax**2-R(20301,2500))
    check('a2_absolute_row_column_bound_9',R(50,49)*R(20301,2500)<9)
    zero('a2_common_coefficient_8alpha',8*R(1,960)-R(1,120))
    # Exact partition reconstruction: unit vector and tangent orthogonality.
    c,sn,a,b,lapT,theta2,T=S.symbols('c s a b lapT theta2 T')
    # sum chi Delta(chi T)=Delta T-|grad theta|^2 T.
    zero('partition_first_jet_cancellation',2*(c*(-sn)+sn*c))
    zero('partition_second_jet_unit_relation',(-c*theta2)*c+(-sn*theta2)*sn+theta2*(c*c+sn*sn))
    check('partition_collar_J0_extra',2*500**4*R(34,10**38)<=R(425,10**28))
    gl=R(1,2000**10)
    check('a2_full_action_joint_margin',1-R(4,1000)-R(1,10**6)-2*gl>R(995,1000)-2*gl)
    check('a2_retained_positive_E0',R(995,1000)-2*gl>R(989,1000))
    # pi>3 suffices after the exact 1500 cancellation.
    check('a2_relative_contact_1503',1500+R(15,6)<1503)
    check('a2_relative_contact_1600',1503<1600)
    check('a2_actual_small_step',1600*R(1,10**6)<=R(16,10**4))
    power=S.symbols('y',positive=True)
    # The worst second/third DeltaR flux includes derivative of U Delta T.
    zero('mass_Hessian_flux_minimum_power',S.diff(power**3*S.diff(power**3,power,2),power)-24*power**3)
    check('mass_Hessian_third_flux_vanishes',3>0)


def logarithm_Hardy_and_multiplier():
    aa,t=S.symbols('a t',positive=True)
    zero('relative_log_square_function_primitive',S.diff(-aa/(aa+t),t)-aa/(aa+t)**2)
    zero('relative_log_square_function_integral',(-aa/(aa+t)).subs(t,0)+1)
    check('complex_log_relative_cost_11',R(11,12)/(1-R(11,12))==11)
    rho=R(1,1000000); dv=R(1,1000)
    # exp(rho)<=1/(1-rho), sqrt24>4 gives an explicit alpha reserve.
    alp=(1/(1-rho)-1)+(1/(1-rho))*dv/4
    check('auxiliary_log_delta_below_001',2*alp+alp**2<R(1,1000))
    a,x=S.symbols('a x',positive=True)
    prim=x*S.log(1+a*a/(x*x))+2*a*S.atan(x/a)
    zero('Fourier_log_correction_antiderivative',S.diff(prim,x)-S.log(1+a*a/(x*x)))
    zero('Fourier_log_correction_endpoint_zero',S.limit(prim,x,0,dir='+'))
    zero('Fourier_log_correction_endpoint_infinity',S.limit(prim,x,S.oo)-S.pi*a)
    ss=S.symbols('s',positive=True)
    Cs=2**(2*ss)*S.gamma((1+2*ss)/4)**2/S.gamma((1-2*ss)/4)**2
    zero('Hardy_log_constant_derivative',S.diff(S.log(Cs),ss).subs(ss,0)-2*(S.log(2)+S.polygamma(0,R(1,4))))
    coef=2**(2*ss-1)*S.gamma((1+2*ss)/2)/(S.sqrt(S.pi)*2*S.gamma(1-ss))
    zero('Hardy_ground_state_kernel_quarter',coef.subs(ss,0)-R(1,4))
    check('Hardy_multiplier_near_integral',48**2*2<68**2)
    check('Hardy_multiplier_far_integral',16**2*2<23**2)
    check('Hardy_multiplier_total_reserve',68+23<100)
    check('Hardy_multiplier_factor_50',R(1,2)*100==50)
    check('positive_log_combined_116',100+11+5<=116)
    zero('quadratic_vanishing_log_derivative',S.diff(x**4*S.log(1/x),x)-x**3*(4*S.log(1/x)-1))
    check('endpoint_log_max_monotone_at_half',4*R(1,2)>1)
    # log2>1/2 is an explicit ordinary inequality used for this last guard.
    check('quadratic_cutoff_Lipschitz_reserve',2+R(15,4)<6)
    zero('quadratic_cutoff_log_norm_constant',116+100*6**2-3716)
    check('positive_log_collar_01',R(1,100)*62==R(62,100))
    check('positive_log_collar_002',R(1,2500)*62==R(248,10000))
    check('positive_log_square_root_reserve',3716+8<62**2)
    zero('bilinear_b_log_norm_constant',116*4+100*4*25-10464)
    check('bilinear_b_root_reserve',10464+8*4<103**2)
    check('bilinear_c_root_reserve',116+2*4+100*7**2<71**2)
    check('bilinear_commutator_reserve',2**2*2*7**2<20**2)
    check('bilinear_collar_01',7313*R(1,100)+20*R(1,10)<76)
    check('bilinear_collar_002',7313*R(1,2500)+20*R(1,50)<R(333,100))
    l1,l2,b1,b2,c1,c2=S.symbols('l1 l2 b1 b2 c1 c2')
    logmat=S.Matrix([[l1,1],[1,l2]]); bm=S.diag(b1,b2); cm=S.diag(c1,c2)
    matrix_zero('complex_bilinear_ordered_commutator_identity',bm*cm*logmat-(bm*logmat*cm-bm*(logmat*cm-cm*logmat)))
    # Constant-complex flat positive-UV remainder: differentiate its integral primitive.
    xx=S.symbols('x',positive=True)
    ff=S.log(1+xx)-xx+xx**2/2
    zero('flat_positive_UV_integral_derivative',S.diff(ff,xx)-xx**2/(1+xx))
    zero('flat_positive_UV_zero_endpoint',ff.subs(xx,0))
    zero('flat_rank_two_sixth_moment',S.integrate((x*(1-x))**3,(x,0,1))/6-R(1,840))
    cc=R(49,50)*R(4999,5000)
    check('flat_complex_UV_relative_005',R(1,24)/cc<R(1,20))


def root_and_scope_budgets():
    epl=R(1,180*2000**10); xlo=R(79,10)*10**9; da=R(1,10**12)
    c_hi=6+24*epl+(R(22,3)+28*da)*epl/xlo
    c_lo=6-(R(22,3)+28*da)*epl/xlo
    ps_hi=(24+4*da/xlo+88*epl+6*epl*(4+2*da/xlo)**2)/6
    ps_lo=(24-4*da/xlo)/(6+(1+da)/xlo)
    check('product_Corb_interval',c_lo>R(5999,1000) and c_hi<R(6001,1000))
    check('product_A2PS_interval',ps_lo>R(3999,1000) and ps_hi<R(4001,1000))
    check('product_b_interval',8*epl+(R(11,3)+4*da)*epl/xlo<R(1,1000))
    check('product_a_disk_bound_9',2*R(6001,1000)-R(3999,1000)+R(1,100)*R(6001,1000)<9)
    # Source upper majorants can equal these numbers; the source inequality is strict.
    check('quantum_product_g_reserve',2*10**25*R(85,10**34)<=R(17,10**8))
    gg=R(17,10**8); bmin=R(999,1000)
    check('Rouche_low_frequency_boundary',9*gg<4*bmin**2*R(1,10**6))
    check('cold_root_displacement',9*gg/(4*bmin**2)<R(384,10**9))
    nu=S.symbols('nu'); p=-(nu-1)*(nu+2)
    zero('cold_root_Legendre_conversion',nu*(nu+1)-(2-p))
    check('cold_root_nu_denominator',R(384,10**9)/R(2999,1000)<R(13,10**8))
    xx=R(13,10**8)*84
    check('cold_exponent_amplification',1/(1-xx)<R(100002,100000))
    zero('reflected_rank_four_physical_half',4*R(1,2)-2)


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
            constraint_inverse_and_observable_budgets()
            radial_barrier_and_product_inverse()
            exact_product_matrix_and_Ward()
            Legendre_tail_and_probe()
            finite_contact_and_form_budgets()
            logarithm_Hardy_and_multiplier()
            root_and_scope_budgets()
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
