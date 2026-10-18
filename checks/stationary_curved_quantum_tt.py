"""Focused exact algebra for the curved seed quantum TT lower form.

Analytic FK, finite-part, geodesic-flow, domain and closure arguments are
ordinary proofs.  This checker computes neither a coupled kernel nor a
physical scalar phase, inverse or nonlinear solution.
"""
from __future__ import annotations
import base64, contextlib, hashlib, json, platform, sys, traceback
from datetime import datetime, timezone
from pathlib import Path
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT/'evidence/stationary-curved-quantum-tt.json'
LOG = ROOT/'evidence/stationary-curved-quantum-tt.log'
INPUTS = [
    'checks/stationary_curved_quantum_tt.py',
    'notes/stationary-curved-quantum-tt.md',
    'notes/stationary-seed-quantum-form.md',
    'notes/stationary-centered-canonical-response.md',
    'notes/stationary-reference-cold-bolt.md',
    'notes/stationary-reference-nonbolt-stress.md',
    'notes/stationary-centered-conormal-source.md',
    'notes/stationary-global-tt-reduction.md',
    'notes/stationary-centered-covariant-seed.md',
    'notes/stationary-centered-exterior-source.md',
    'notes/stationary-horizon-source-defect.md',
    'notes/stationary-weyl-radial-sources.md',
    'notes/spherical-anomaly-sources.md',
    'notes/spherical-anomaly-action.md',
    'notes/stationary-paired-dirac-phase.md',
    'notes/stationary-massive-boundary-layer.md',
    'notes/stationary-exterior-conformal-domain.md',
    'notes/stationary-exterior-boundary-weyl.md',
]
CONTRACT = {
    'acceptance': 'exact residuals zero and stated rational comparisons satisfied; strict source inputs permit equality of rational majorants',
    'numerical_tolerances': 'none',
    'scope': 'finite identities and conservative budgets supporting the fixed-seed curved TT form',
    'ordinary_arguments': ['Dirac and paired-domain realization',
        'FK and complete Born estimates', 'common angular finite parts',
        'geodesic-flow and spectral positivity', 'logarithmic Hardy inequality',
        'localization and finite-part decomposition', 'smooth-bolt and outer form closure'],
    'unresolved': ['physical coupled kernel and scalar phase', 'actual source compatibility',
        'quantitative full coupled inverse', 'nonlinear physical branch and matching'],
    'failure_policy': 'append initialized scientific runs and retain all failed raw input bytes; unreadable-input or receipt preflight errors abort before scientific execution',
    'hash_policy': 'SHA256 of raw UTF8/LF/noBOM bytes; no normalization',
}
ROWS = []
R = S.Rational


def check(label, condition, detail=''):
    row = {'label': label, 'passed': bool(condition), 'detail': str(detail)}
    ROWS.append(row)
    print(json.dumps(row, sort_keys=True), flush=True)


def zero(label, expression):
    value = S.simplify(expression)
    check(label, value == 0, value)


def flat_and_common():
    # DLMF product integral, evaluated at the two independently fixed moments.
    def kk(mu, nu):
        return S.sqrt(S.pi)*S.gamma(mu/2)*S.gamma(mu/2+nu)*S.gamma(mu/2-nu)/(4*S.gamma((mu+1)/2))
    zero('Bessel_second_moment', (kk(S.Integer(4), S.Integer(0))-kk(S.Integer(4), S.Integer(1)))/(2*S.pi**2)+1/(6*S.pi**2))
    zero('Bessel_fourth_moment', (kk(S.Integer(6), S.Integer(0))-kk(S.Integer(6), S.Integer(1)))/(2*S.pi**2)+4/(15*S.pi**2))
    zero('Bessel_fifth_absolute_moment', (kk(S.Integer(7), S.Integer(1))-kk(S.Integer(7), S.Integer(0)))/(2*S.pi**2)-R(225,4096))
    sixth=(kk(S.Integer(8),S.Integer(1))-kk(S.Integer(8),S.Integer(0)))/(2*S.pi**2)
    zero('Bessel_sixth_absolute_moment',sixth-48/(35*S.pi**2))
    zero('full_angular_sixth_Taylor_remainder',2*S.pi*sixth/S.factorial(6)-2/(525*S.pi))
    n, vp, vm = S.symbols('n nu_plus nu_minus', positive=True)
    zero('angular_shift_partial_fractions', 2*(2*n+vp+vm)/((n+vp)*(n+vm))-2/(n+vp)-2/(n+vm))
    t, q = S.symbols('t q', positive=True)
    zero('angular_shift_derivative', S.diff(2*(2*n+q)/(n*(n+q)+t), t)+2*(2*n+q)/(n*(n+q)+t)**2)
    z = S.symbols('z', positive=True)
    zero('digamma_recurrence', S.expand_func(S.polygamma(0, z+1))-S.polygamma(0, z)-1/z)
    a = S.symbols('alpha', real=True)
    aa = a*(1-a)
    zero('heat_parameter_second_moment', S.integrate(aa**2, (a, 0, 1))-R(1, 30))
    zero('heat_parameter_third_moment', S.integrate(aa**3, (a, 0, 1))-R(1, 140))
    # Integrate polynomial*log(alpha) termwise; reflection doubles it.
    poly = S.Poly(S.expand(aa**2), a)
    logmoment = 2*sum(-coef/R(power[0]+1)**2 for power, coef in poly.terms())
    zero('Kpp_parameter_log_constant', logmoment+R(47, 900))
    zero('Kpp_constant_107_over_30', 2+R(47, 30)-R(107, 30))
    ll = S.symbols('log_scale', real=True)
    lower = z**2*(S.log(z)/2-ll-R(38, 15))
    zero('Kpp_twice_integration_constants', S.diff(lower, z, 2)-(S.log(z)-2*ll-R(107, 30)))
    zero('K_lower_first_derivative_origin', S.limit(S.diff(lower, z), z, 0, dir='+'))
    zero('K_lower_value_origin', S.limit(lower, z, 0, dir='+'))
    u = S.symbols('u', nonnegative=True)
    rem = S.log(1+u)-u+u**2/2
    zero('flat_remainder_derivative', S.diff(rem, u)-u**2/(1+u))
    zero('flat_third_derivative', S.diff(rem, u, 3)-2/(1+u)**3)
    ss, fg, scale = S.symbols('s Fang log_scale')
    finite = S.exp(2*ss*scale)*(2+fg*ss)
    zero('common_Feff_no_plus_two', S.diff(finite, ss).subs(ss, 0)-fg-4*scale)
    A, du, cf = S.symbols('A Delta_u C_Feff', real=True)
    actual, geod = cf*A*du/(20*S.pi), cf*A*du/(45*S.pi)
    zero('a2_actual_minus_geodesic_second_contact', actual-geod-cf*A*du/(36*S.pi))
    # Independently apply the geodesic fourth polynomial at a normal quadratic.
    zero('geodesic_fourth_operator_at_quadratic', -cf*(3*(-2*A*du)+(-2*A)*du)/(360*S.pi)-geod)
    m, ev = S.symbols('M eigen_A_squared', positive=True)
    zero('constant_probe_resolvent_algebra', m**2*(m**2+ev)/(m**2-ev)**2-m**2/(m**2-ev)-(2*m**4/(m**2-ev)**2-2*m**2/(m**2-ev)))
    b1, b2 = S.symbols('A1 A2', cls=S.Function)
    # The heat weight identity uses the exact termwise relation A1'=-A2.
    zero('p0_heat_weight_derivative', (-2*S.diff(t*b1(t), t)).subs(S.diff(b1(t), t), -b2(t))-(2*t*b2(t)-2*b1(t)))
    v, a0, a2 = S.symbols('v A0 A2', real=True)
    action = -a0*S.exp(2*v)/2+a2*S.exp(-2*v)/2
    zero('p0_common_local_restoration', S.diff(action, v, 2).subs(v, 0)-(-2*a0+2*a2))
    # Independent finite-matrix witness of the noncommutative resolvent
    # commutator used before taking the separated-support diagonal.
    zz=S.symbols('complex_mass',positive=True)
    mat=S.Matrix([[0,1],[-1,0]])
    chi=S.diag(R(1,3),R(2,3))
    inv=(mat+zz*S.eye(2)).inv()
    comm=mat*chi**2-chi**2*mat
    residual=inv*comm*inv-(chi**2*inv-inv*chi**2)
    vals=[S.simplify(x) for x in residual]
    check('first_order_resolvent_commutator_matrix', all(x==0 for x in vals), vals)
    diag=[S.simplify((chi**2*inv-inv*chi**2)[j,j]) for j in range(2)]
    check('first_order_commutator_diagonal_zero', all(x==0 for x in diag), diag)


def geometric_and_cutoff():
    rho = S.symbols('rho', positive=True)
    mu = S.symbols('mu', positive=True)
    # Counting N(E)<=4CE^2 gives these exact polynomial-exponential constants.
    zero('angular_linear_exponential_M2_constant', 4*(2*S.factorial(3)+S.factorial(4))-144)
    zero('angular_Gaussian_M2_constant', 4*(S.factorial(1)+S.factorial(2))-12)
    zero('angular_Gaussian_M4_constant', 4*(2*S.factorial(2)+S.factorial(3))-40)
    check('p0_absolute_heat_weight', 2*(40+12)==104<120)
    check('p0_half_gap_heat_weight', 2*(40*8+12*4)==736<800)
    check('geodesic_tail_J0_J2', R(72,4)==18 and R(72,2)==36)
    check('geodesic_tail_J4_constant', 72*16==1152)
    # E1(9/400)<5 follows from 400/9<e^4 and e^-1<1.
    check('E1_tail_rational_guard', sum(R(4)**j/S.factorial(j) for j in range(9))>R(400,9))
    check('global_flow_quadratic_tail', 7*(R(18*400**4,1024)+36*400**2+5*1152*5)<4*10**9)
    check('flow_quadratic_fourth_jet', 8*3*5==120)
    check('local_geodesic_fourth_ODE', 3*10**6*R(102,100)**4+5*6000*R(102,100)**2*16+2*15*16**2+2*15*R(102,100)*7000<4*10**6)
    check('local_geodesic_quadratic_remainder', R(4*10**6,12)+R(25,4*128)+R(5*4*10**6,24*128**2)+R((4*10**6)**2,24**2*128**3)<4*10**5)
    check('local_flow_fifth_moment_contact', 7*4*10**5*R(1,2)*R(2,3)<2*10**6)
    # Bernstein coefficients of the actual beta(6,6) cutoff, not sampled derivatives.
    x = S.symbols('x', real=True)
    cutoff = 1-S.integrate(2772*x**5*(1-x)**5, (x,0,x))
    zero('beta_cutoff_left_value', cutoff.subs(x,0)-1)
    zero('beta_cutoff_right_value', cutoff.subs(x,1))
    derivative_bounds = [11,220,3960,63360,887040]
    for j, bound in enumerate(derivative_bounds, 1):
        pp = S.Poly(S.diff(cutoff,x,j),x)
        deg = pp.degree()
        power = {i[0]: c for i,c in pp.terms()}
        bern = [sum(power.get(i,0)*S.binomial(k,i)/S.binomial(deg,i) for i in range(k+1)) for k in range(deg+1)]
        check('beta_cutoff_Bernstein_order_'+str(j), max(abs(b) for b in bern)<=bound)
    radial = [1,2,6,36,300]
    Cj = [1]
    for j in range(1,6):
        Cj.append(2**j*sum(derivative_bounds[k-1]*S.bell(j,k,tuple(radial[:j-k+1])) for k in range(1,j+1)))
    check('Cartesian_cutoff_Bell_array', Cj==[1,22,968,42768,1907136,86359680])
    for j in range(1,6):
        num=Cj[j]+2*j*Cj[j-1]+(j*(j-1)*Cj[j-2] if j>=2 else 0)
        check('quadratic_probe_normalized_jet_'+str(j), R(32)**(j-2)*num/(S.factorial(j)*100**j)<3)
    check('normal_chart_H_domain', R(1,16)+R(11,20)*R(1,16)**2<R(7,100))
    check('normal_chart_Hessian', R(118,100)*R(171,160)**2+R(121,100)+R(121,100)/R(149,160)**2<5)
    check('normal_chart_cubic', 14*R(171,160)**3+3*R(118,100)*R(171,160)*R(11,10)+2*R(11,10)**3/R(149,160)**3<30)
    check('normal_chart_v_cubic_remainder', R(101,100)*(30+3*5*R(5,16)+R(5,16)**3)<36)
    check('global_cap_cubic_guard', R(26,1000)*64**3+8*64<10000)
    check('midpoint_Jacobian_bound', R(1,2)*R(1,64)**2<1)
    # Proper-to-normal derivative reserve and the resulting scalar conversion.
    check('proper_Y_fourth_jet', R(109,100)*(14+3*R(11,10)*R(118,100)+R(11,10)**3)<21)
    check('normal_Y_fourth_composition', 21*R(171,160)**4+6*R(26,10)*R(171,160)**2*R(11,10)+3*R(12,10)*R(11,10)**2<52)


def contacts_and_spectral():
    t, V0, L = S.symbols('t V0 LipV', positive=True)
    check('Hardy_aux_potential', 21+36*R(1,6)**2<=22<24)
    check('Hardy_aux_third_jet', 324*10+30==3270)
    check('Hardy_aux_potential_Lipschitz', R(3270,2)+R(12,10)*42/2<2000)
    check('scalar_density_cubic_guard', R(103,100)*(R(3270,2)+R(3,4)*R(12,10)*42+R(1,8)*R(12,10)**3)<1800)
    check('scalar_half_gap_factor', (1-R(22,1296))**-2<R(104,100))
    check('scalar_beta_tail_constant', R(104,100)/5<R(1,4))
    check('scalar_FK_M3_constant', R(8000*3,2)*R(1,32)/(2*3)<100)
    check('scalar_FK_M4_constant', R(484*2,1)*R(1,64)*R(106,100)/(2*3)<3)
    check('scalar_FK_angular_budget', R(200,3)+12*R(352,100000)<70)
    check('scalar_compact_quadratic_tail', 2*64**2<10000)
    check('auxiliary_spectral_escape_tail', R(3*400**4,2*3)<13*10**9)
    h = S.Function('h')(t)
    potential=S.diff(h,t,2)/2+S.diff(h,t)**2/4
    zero('scalar_model_constant_base_mode', -S.diff(S.exp(h/2),t,2)+potential*S.exp(h/2))
    check('Born_one_fifth_moment_budget', (64**3+10000+512)*1000<3*10**8)
    check('Born_two_rank_and_tags', 4*(10000+512)*12288000<6*10**11)
    xx=S.symbols('x', nonnegative=True)
    check('Born_tagged_Neumann_sum', ((3-2*xx)/(1-xx)**2).subs(xx,R(13,100))<4)
    check('Born_higher_six_distance_budget', 2*720*4**6*100*5*8**2*R(1026,1000)*4+2*2*720*4**6*100*8**2<10**12)
    check('Born_angular_second_contact', 3*10**8*R(2,3)+6*10**11*R(2,3)+10**13*4*R(352,100000)+4*10**9<10**13)
    check('complex_UV_changed_ratio', R(13,100)*R(5,3)**2<R(362,1000))
    check('complex_UV_multiplier', R(5,3)**12*R(4,3)**4*R(87,100)/R(638,1000)<R(3*10**14,14*10**10))
    check('cold_constant_probe_remainder', R(48*10**14,1)*R(352,100000)**3<R(21*10**7,1))
    # Buffer repair includes evaluation-radius shrink, not only exponent changes.
    check('local_row_evaluation_buffer', R(19,200)/18<R(1,132))
    check('far_row_evaluation_buffer', R(19,400)/18<R(1,263))
    check('local_buffered_Schur', 1000*(78*10**6)**2*8*135**2<10**24)
    check('far_buffered_Schur', 2*(30000*64*512)**2*8*265**2<2*10**24)
    check('complex_mass_diagonal_Schur', 5*8*(78*10**7)**2*8*130**2<R(33*10**23,1))
    check('complex_mass_Theta_normalization', R(33*10**23,1)*R(352,100000)**2<R(41*10**18,1))
    check('p0_full_budget', R(41*10**18)+R(3*10**16)+R(3*10**13)+R(15*10**9)+R(21*10**7)<10**20)


def localization_and_IMS():
    x=S.symbols('x',real=True)
    P=462*x**6-1980*x**7+3465*x**8-3080*x**9+1386*x**10-252*x**11
    zero('localization_beta_density', S.diff(P,x)-2772*x**5*(1-x)**5)
    check('localization_beta_first_sharp', R(2772,4**5)<3)
    pb=[11,220,3960,63360,887040,10644480]
    for j,bound in enumerate(pb,1):
        pp=S.Poly(S.diff(P,x,j),x); degree=pp.degree()
        powers={i[0]:c for i,c in pp.terms()}
        bern=[sum(powers.get(i,0)*S.binomial(k,i)/S.binomial(degree,i) for i in range(k+1)) for k in range(degree+1)]
        check('localization_Bernstein_derivative_'+str(j),max(abs(b) for b in bern)<=bound)
    th=[5,350,6300,100000,1400000,17000000]
    check('localization_theta_first', R(11,7)*3<5)
    for j in range(2,7):
        check('localization_theta_derivative_'+str(j),R(11,7)*pb[j-1]<th[j-1])
    cb=[5,375,12000,650000,38000000,3000000000]
    for j in range(1,7):
        bound=sum(S.bell(j,k,tuple(th[:j-k+1])) for k in range(1,j+1))
        check('localization_chi_Bell_derivative_'+str(j),bound<=cb[j-1],bound)
    Yj=[1,1,10,100,1000,10000]
    sixth=[S.bell(6,k,tuple(Yj[:7-k])) for k in range(1,7)]
    check('geodesic_sixth_chain_array',sixth==[10000,8500,2115,245,15,1],sixth)
    chi=[1]+[cb[j-1]*100**j for j in range(1,7)]
    check('geodesic_partition_X6',sum(sixth[k-1]*chi[k] for k in range(1,7))<10**22)
    fourth=[S.bell(4,k,tuple(Yj[:5-k])) for k in range(1,5)]
    check('geodesic_partition_X4',sum(fourth[k-1]*chi[k] for k in range(1,5))<10**14)
    # Differentiate Y''=h'(Y)(1-Y'^2) using only the stated real jets.
    check('geodesic_Y5_product_ledger',100+3*12*2+3*2*22+R(1,20)*260<564)
    check('geodesic_Y6_product_ledger',1130+4*100*2+6*12*22+4*2*260+R(1,20)*3400<8994)
    RR=R(1,400)
    # q(s)=s^2 xi(s): ordinary Leibniz, including the sixth cutoff derivative.
    xi=[1]+[pb[j-1]*2**j for j in range(1,7)]
    qcoef=[1,24,970,37100,1280000,40000000,1100000000]
    for j in range(1,7):
        derived=xi[j]+2*j*xi[j-1]+(j*(j-1)*xi[j-2] if j>=2 else 0)
        check('IMS_quadratic_probe_derivative_'+str(j),derived<=qcoef[j])
    qj=[qcoef[j]*RR**(2-j) for j in range(7)]
    product=[sum(S.binomial(j,k)*chi[k]*qj[j-k] for k in range(j+1)) for j in range(7)]
    check('IMS_chi_quadratic_sixth_product',product[6]<6*10**19)
    check('IMS_quadratic_X6',sum(sixth[k-1]*product[k] for k in range(1,7))<10**20)
    check('IMS_quadratic_X4',sum(fourth[k-1]*product[k] for k in range(1,5))<10**12)
    check('IMS_full_angular_point_coefficient',R(8,525*3)*R(352,100000)<R(2,100000))
    check('IMS_p0_point_contact_budget',R(2,100000)*10**22+10**14<10**22)
    check('IMS_p2_point_contact_budget',R(2,100000)*10**20+10**12<10**22)
    check('IMS_principal_kernel_on_anchor',144+2*10**8*RR**2<10000)
    check('IMS_gain_two_separated_kernel',R(25,2)*100**2*10000==1250000000)
    check('IMS_fullflow_cosine_tail',500/RR**4<2*10**13)
    # Exact pointwise weighted fourth IMS identity with k(Y), not only k=1.
    y=S.symbols('Y',real=True)
    theta=S.Function('theta')(y); fun=S.Function('T')(y)
    k=S.Function('weight')(y); warp=S.Function('warp_log_derivative')(y)
    p=S.diff(theta,y); p2=S.diff(theta,y,2)
    f1=S.diff(fun,y); f2=S.diff(fun,y,2)
    lap=f2+warp*f1
    c1=S.cos(theta); c2=S.sin(theta)
    def laplace(g): return S.diff(g,y,2)+warp*S.diff(g,y)
    point=k*(laplace(c1*fun)**2+laplace(c2*fun)**2-lap**2)
    direct=k*(4*p**2*f1**2+((p2+warp*p)**2+p**4)*fun**2+4*p*(p2+warp*p)*fun*f1-2*p**2*fun*lap)
    zero('weighted_fourth_IMS_point_identity',S.trigsimp(S.expand(point-direct)))
    aa=2*k*p**2
    cross=S.diff(aa,y)+4*k*p*(p2+warp*p)
    flux=cross*fun**2/2-aa*fun*f1
    potential=k*((p2+warp*p)**2+p**4)-(S.diff(cross,y)+warp*cross)/2
    integrated=6*k*p**2*f1**2+potential*fun**2+S.diff(flux,y)+warp*flux
    zero('weighted_fourth_IMS_flux_identity',S.expand(direct-integrated))
    zero('second_IMS_localized_gradient',S.trigsimp(S.expand(S.diff(c1*fun,y)**2+S.diff(c2*fun,y)**2-f1**2-p**2*fun**2)))
    # Transition-weighted physical measure and literal SUMS of norms.
    check('transition_weighted_measure',R(52,10**30)*R(16,25)*R(1,100)**4<R(34,10**38))
    # Equality of these rational majorants is intentional; the preceding
    # actual transition measure is strict, hence so is the physical bound.
    check('localized_gradient_extra',R(34,10**38)*500**2<=R(85,10**33))
    check('localized_gradient_sum',R(17,10**28)+R(85,10**33)<R(18,10**28))
    check('literal_J1_aggregate',R(17,10**28)+R(18,10**28)+2*R(18,10**28)+R(1,10**31)<R(1,10**26))
    check('literal_J0_aggregate',3*R(52,10**30)<R(2,10**28))
    check('auxiliary_H2_conversion',2*22**2==968)
    # Literal total loss from the ordinary three-remainder argument is
    # epsilon*(pi/2+2/3), not 3*(epsilon/3)*2.
    loss_upper=R(1,1000)*(R(11,7)+R(2,3))
    check('three_fractional_remainders_H2_loss',loss_upper<R(224,100000))
    check('three_remainder_loss_below_retained_small_margin',loss_upper<R(7,2200))
    # pi<22/7 lowers the initial 33/(600pi); pi>3 raises the target
    # 1/(100pi).  Thus this rational comparison proves the stated
    # surviving fourth margin without evaluating pi numerically.
    check('outer_fourth_margin_after_three_remainders',R(33*7,600*22)-loss_upper>R(1,300))
    check('fractional_split_radius',R(1,4*10**16)<RR**2)
    check('fractional_split_log_guard',sum(R(14)**j/S.factorial(j) for j in range(31))>500000)
    check('fractional_gradient_cost',10**13*(R(1,4)+2*14)<10**15)


def weyl_and_physical():
    v, om, rr, uu, vv, yy, kk = S.symbols('v omega R U V y k',real=True)
    st, tt, de, ds, ww, jj = S.symbols('S T DeltaT DeltaS w j',real=True)
    r_v=S.exp(-2*v*tt)*(rr-2*v*de)
    u_v=S.exp(-2*v*tt)*(uu+v*ds)
    v_v=S.exp(-2*v*tt)*(vv+2*v*ww+v**2*jj)
    o_v=om+v*st
    density=S.exp(2*v*tt)
    lag=o_v*(r_v**2/60-yy*r_v/18+kk)+r_v*(u_v+v_v)/30-yy*v_v/18-u_v**2/10-R(7,90)*u_v*v_v-R(7,180)*v_v**2
    q2=S.expand(S.diff(density*lag,v,2).subs(v,0)/2)
    ff=ds+de/3
    DD=st*(-rr/15+yy/9)+tt*((6*om+1)*rr/45+11*vv/135)-11*ww/135
    LL=tt*(-rr/15+2*uu/5+7*vv/45)-7*ww/45
    JJ=rr/30-7*(uu+vv)/90-yy/18
    WW=-2*rr/15+14*(uu+vv)/45
    ZZ=tt**2*(om*rr**2/30+rr*(uu+vv)/15-uu**2/5-7*uu*vv/45-7*vv**2/90+2*kk*om)+st*tt*(-rr**2/30+2*kk)
    completed=-ff**2/10+(om+R(1,6))*de**2/15+de*DD+ff*LL+JJ*jj+WW*tt*ww-7*ww**2/45+ZZ
    zero('W4_exact_complete_Hessian', S.expand(q2-completed))
    zero('W4_Young_fourth_shift', R(1,6)+R(30,200)-R(19,60))
    check('W4_lower_order_ledger', 25*20000+5*140000+150000+130000+76000<1600000)
    check('W4_T2_coefficient', R(169,30)+R(351,15)+R(22**2,5)+R(7*22*5,45)+R(7*25,90)<145)
    check('W4_Z_integrated_ledger', 145*500+R(169,30)*500<76000)
    # Reference chain rule includes the often-missed second variation of X.
    X = S.symbols('X',positive=True)
    pot=S.Function('DeltaV')
    ref=S.exp(2*v*(st+tt))*pot(X*S.exp(2*v*st))
    ref2=S.diff(ref,v,2).subs(v,0)/4
    expected=(st+tt)**2*pot(X)+2*st*(st+tt)*X*S.diff(pot(X),X)+st**2*(X*S.diff(pot(X),X)+X**2*S.diff(pot(X),X,2))
    zero('reference_second_X_chain_rule', ref2-expected)
    check('reference_complete_ledger', 3341*1100+6682*332+(7000+3341)*100<7*10**6)
    # Strict rational logarithm guards use log2<7/10 and pi<22/7.
    check('Hardy_fourth_reserve', 9-R(38,15)-R(19,60)-4-R(7,20)>R(33,20))
    # The rational majorants meet at equality; the published source bound
    # Gq^2/(lambda B^4)<4e-27 makes the ACTUAL physical result strict.
    check('p0_physical_cost', 1500*10**20*R(4,10**27)<=R(6,10**4))
    check('gluing_row_subtracted_physical_cost', 1000*6*10**24*R(85,10**34)*R(616,100)<R(32,10**5))
    check('reference_physical_cost', 224*10**6*R(22,7)*R(4,10**27)<R(3,10**18))
    # Combined finite local lower-order reserves; analytic localization proof
    # supplies the meaning of each coefficient, not this rational sum.
    check('local_physical_aggregate_below_point004', R(1,1000)+R(2,1000)+R(6,10**4)+R(1,10**12)<R(4,1000))
    check('full_TT_final_comfortable_margin', 1-R(5,1000)-R(1,10**6)>R(99,100))


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
            flat_and_common()
            geometric_and_cutoff()
            contacts_and_spectral()
            localization_and_IMS()
            weyl_and_physical()
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
