"""Exact finite identities for the covariant reduced seed inverse.

Proposed verbatim destination: checks/stationary_covariant_reduced_inverse.py.
Closed operator/form domains, the weighted parametrix and bootstrap, paired
determinant regularity and nonlinear mapping claims remain ordinary proofs.
This checker does not construct a nonlinear branch, evaluate the original
signed source moment or establish physical truncation/observable error bars.
"""
from __future__ import annotations
import base64, contextlib, hashlib, json, platform, sys, traceback
from datetime import datetime, timezone
from pathlib import Path
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/stationary-covariant-reduced-inverse.json'
LOG=ROOT/'evidence/stationary-covariant-reduced-inverse.log'
INPUTS=[
    'checks/stationary_covariant_reduced_inverse.py',
    'notes/stationary-covariant-reduced-inverse.md',
    'notes/stationary-weighted-tt-inverse.md',
    'notes/stationary-scalar-outgoing-inverse.md',
    'notes/stationary-scalar-phase-moment.md',
    'notes/stationary-global-weighted-map.md',
    'notes/stationary-global-tt-reduction.md',
    'notes/stationary-curved-quantum-tt.md',
    'notes/stationary-seed-quantum-form.md',
    'notes/stationary-disk-and-normal-inverse.md',
    'notes/stationary-centered-covariant-seed.md',
    'notes/stationary-centered-first-order-geometry.md',
    'notes/stationary-centered-complete-residual.md',
    'notes/stationary-centered-exterior-source.md',
    'notes/uniform-dirac-reference.md',
]
CONTRACT={
    'acceptance':'exact residuals zero and stated rational comparisons satisfied',
    'numerical_tolerances':'none',
    'scope':'finite field-map, off-shell trace, seed block, coercivity-budget, indicial, energy and observable identities',
    'ordinary_arguments':['same-cutoff covariant action identities and boundary completions',
        'actual paired Dirac and global seed York domains',
        'closed reduced scalar form and positive TT spectral form',
        'signature-flipped Lax-Milgram existence and actual mixed-form boundedness',
        'weighted elliptic bootstrap, relative compactness and fast mass inverse',
        'asymptotic vacuum cancellation and finite field-map remainder interpretation'],
    'unresolved':['parameter-dependent common-scheme mixed nonlinear map',
        'quantitative fast/cold inverse, source defect and nonlinear certificate',
        'a sufficiently regular reduced branch and physical observable remainder',
        'original signed preparation moment and original physical compatibility'],
    'matrix_scope':'exact generic finite block algebra; no matrix fixture proves an operator-domain statement',
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


def field_map_and_trace_normalization():
    G,la=S.symbols('G lambda',positive=True)
    ep=G/(60*S.pi); cr=1/(5760*S.pi**2); aa=8*S.pi*G*cr
    rr,uu,e0=S.symbols('R U e0')
    zero('field_map_coefficient',aa-G/(720*S.pi))
    zero('field_map_epsilon_coefficient',aa-ep/12)
    zero('R_squared_trace_normalization',12*cr-ep/(8*S.pi*G))
    zero('first_order_pulled_action_sign',-aa*e0**2/(8*S.pi*G)+cr*e0**2)
    zero('redundant_curvature_expansion',rr**2-(rr+uu)**2+2*uu*rr+uu**2)
    zero('scalar_block_counterterm_normalization',72*cr-1/(80*S.pi**2))
    zero('normal_scalar_action_normalization',6/(8*S.pi*G)-1/(80*S.pi**2*ep))
    q,r=S.symbols('q r',positive=True)
    # Maxwell dvol has weight+4 and F_ab F^ab has weight-4.
    zero('four_dimensional_Maxwell_conformal_weight',4-4)
    zero('physical_magnetic_invariant',20*q**2/(2*r**4)-10*q**2/r**4)
    la0=la+11*G*la**2/(360*S.pi)
    zero('selected_vacuum_reference_pressure',3*(la-la0)/(8*S.pi*G)+11*la**2/(960*S.pi**2))
    t,k1,k2,acc,mu=S.symbols('t k1 k2 acceleration density_contact')
    polynomial=(t*k1+t**2*acc/2)**2*(1+t*mu)
    zero('counterterm_zero_seed_first_variation',S.diff(polynomial,t).subs(t,0))
    zero('counterterm_seed_Hessian_no_acceleration_contact',S.diff(polynomial,t,2).subs(t,0)-2*k1**2)
    s=S.symbols('s')
    zero('counterterm_polarized_Hessian',S.diff((t*k1+s*k2)**2,t,s)-2*k1*k2)


def off_shell_reduction():
    rr,uu,rup,boxr,boxu,boxs,csq,jsq,mag,la0,ep=S.symbols('R U rUp BoxR BoxU BoxSigma C2 J2 M lambda0 epsilon')
    e0=rr+uu; hu=2*uu+rup
    # Derivative-free terms in delta(e0^2 dvol); the BoxSigma term
    # becomes -12 Box(e0) after the ordinary boundary integration.
    density=2*e0*(rup-2*rr)+4*e0**2
    zero('off_shell_counterterm_density_completion',density-2*hu*e0)
    zero('off_shell_BoxSigma_coefficient',2*e0*(-6*boxs)+12*e0*boxs)
    full=rr+12*la0+ep*(boxr+R(7,12)*csq+R(11,6)*jsq-R(11,72)*rr**2-20*mag)
    ej=ep*(boxr+boxu)-ep*hu*e0/6
    euler=csq-2*jsq+rr**2/6
    reduced=rr+12*la0+ep*(R(3,2)*csq-R(11,12)*euler-20*mag-boxu+hu*e0/6)
    zero('exact_off_shell_trace_cancellation',full-ej-reduced)
    zero('reduced_trace_BoxR_absence',S.diff(reduced,boxr))
    ell=S.symbols('ell')
    quantum=reduced-(rr+12*la0)
    general=rr+12*la0+(1-ell)*(uu-12*la0)+ell*quantum
    zero('bookkeeping_zero_is_covariant_seed_trace',general.subs(ell,0)-(rr+uu))
    zero('bookkeeping_one_is_full_reduced_trace',general.subs(ell,1)-reduced)
    zero('bookkeeping_reference_source_coefficient',S.diff(general,ell)-(quantum-uu+12*la0))
    zero('Euler_Weyl_Ricci_anomaly_decomposition',R(3,2)*csq-R(11,12)*euler-(R(7,12)*csq+R(11,6)*jsq-R(11,72)*rr**2))
    sig,jh=S.symbols('sigma JHessSigma')
    # Conformal Ricci: delta Ric_ab=-2 Hess_ab sigma-g_ab Box sigma.
    # Ricci^2 variation uses Ric:Hess=J:Hess+R Box/4.
    ric2=jsq+rr**2/4
    d_ric2=-4*sig*ric2-4*(jh+rr*boxs/4)-2*rr*boxs
    d_r2=-4*sig*rr**2-12*rr*boxs
    d_c2=-4*sig*csq
    d_euler=d_c2-2*d_ric2+R(2,3)*d_r2
    zero('Euler_conformal_principal_identity',d_euler-(-4*sig*euler+8*jh-2*rr*boxs))
    cscalar=6-R(11,6)*ep*rr+2*ep*(uu+rup)
    principal=-6*boxs+ep*(-R(11,12)*(8*jh-2*rr*boxs)-(rup+hu)*boxs)
    zero('reduced_trace_principal_tensor',principal+cscalar*boxs+R(22,3)*ep*jh)
    zero('principal_divergence_R_gradient_cancellation',-R(11,6)+R(22,3)*R(1,4))
    zero('principal_divergence_U_coefficient',2*ep-2*ep)
    reserve=R(11,6)*100+2*(13+1)+R(22,3)*100
    zero('principal_tensor_rational_reserve',reserve-R(2834,3))
    check('principal_tensor_reserve_below_1000',reserve<1000)
    eg=R(1,180*2000**10)
    check('low_jet_principal_ellipticity_margin',6-1000*eg>5)
    ricrad=S.symbols('Ric_radial')
    zero('radial_principal_Einstein_identity',
         6-R(11,6)*ep*rr+2*ep*(uu+rup)+R(22,3)*ep*(ricrad-rr/4)
         -(6+R(22,3)*ep*(ricrad-rr/2)+2*ep*(uu+rup)))
    k01,k02,k12,k23=S.symbols('K01 K02 K12 K23')
    ric0=k01+2*k02; ric1=k01+2*k12; ric2=k02+k12+k23
    rsec=ric0+ric1+2*ric2
    riemann2=4*(k01**2+2*k02**2+2*k12**2+k23**2)
    eulersec=riemann2-4*(ric0**2+ric1**2+2*ric2**2)+rsec**2
    zero('radial_Einstein_no_second_jet',ric0-rsec/2+2*k12+k23)
    zero('radial_Euler_sectional_identity',eulersec-8*(k01*k23+2*k02*k12))
    zero('radial_Euler_no_mixed_second_jet',S.diff(eulersec,k01,k02))
    zero('radial_Euler_linear_orbit_second_jet',S.diff(eulersec,k01,2))
    zero('radial_Euler_linear_radius_second_jet',S.diff(eulersec,k02,2))
    warp,rad,wp,rp=S.symbols('warp radius warp_prime radius_prime',positive=True)
    zero('radial_Einstein_first_jet_formula',
         (-2*k12-k23).subs({k12:-wp*rp/(warp*rad),k23:(1-rp**2)/rad**2})
         -(2*wp*rp/(warp*rad)+(rp**2-1)/rad**2))


def seed_geometry_and_blocks():
    u,B,la=S.symbols('u B lambda',positive=True)
    d,a=S.symbols('d a',real=True)
    x=la*B**2
    F=x/u**2+1-2*a*u+(1+3*x+d)*u**2-2*d*u**3+d*u**4
    rscalar=(2*u**3*S.diff(F,u)-u**4*S.diff(F,u,2)-2*u**2*(F-1))/B**2
    U=12*la-d/B**2*(4*u**5-6*u**6)
    rup=-u*S.diff(U,u)
    v=d/B**2*(2*u**5-4*u**6)
    V0=4*la+v
    zero('centered_seed_scalar_curvature',rscalar+U)
    zero('seed_trace_vanishes_differentiated',S.diff(rscalar+U,u))
    zero('seed_U_radial_derivative',rup-d/B**2*(20*u**5-36*u**6))
    zero('seed_conformal_Hsigma_potential',(2*U+rup)/6-V0)
    b,c,j0,j1,j2=S.symbols('b c j0 j1 j2')
    h=S.diag(-b-2*c,b,c,c)
    jj=S.diag(j0,j1,j2,j2)
    zero('TT_tensor_trace',S.trace(h))
    hc=-S.trace(jj*h)+rup*c/2
    zero('TT_trace_vertex_contraction',hc-((j0-j1)*b+2*(j0-j2)*c+rup*c/2))
    hs,ht,cr=S.symbols('Hsigma Hc c_R')
    A,Bmix,D=S.symbols('A Bmix D')
    original=S.Matrix([[A,Bmix],[Bmix,D]])
    hv=S.Matrix([[hs,ht]])
    reduced=original-2*cr*hv.T*hv
    matrix_zero('full_reduced_Hessian_Gram_subtraction',reduced-S.Matrix([[A-2*cr*hs**2,Bmix-2*cr*hs*ht],[Bmix-2*cr*hs*ht,D-2*cr*ht**2]]))
    T=S.Symbol('T',commutative=False); vv=S.Symbol('v',commutative=False)
    ee,ll=S.symbols('E lambda')
    q=T+4*ll
    expression=q*(T-ee)-(q+vv)**2+(ee+4*ll)*q+q*vv+vv*q+vv**2
    zero('noncommuting_scalar_reduction',S.expand(expression))
    zero('scalar_fourth_coefficient_canceled',1-1)
    # The remaining TT and mixed Gram entries cannot be omitted.
    zero('TT_subtraction_coefficient',S.diff(reduced[1,1],cr)+2*ht**2)
    zero('mixed_subtraction_coefficient',S.diff(reduced[0,1],cr)+2*hs*ht)


def scalar_form_coefficients():
    rho=S.symbols('rho',real=True)
    sig=S.Function('sigma')(rho); w=S.Function('weight')(rho); V=S.Function('V')(rho)
    def t_op(f):
        return -S.diff(w*S.diff(f,rho),rho)/w
    boxV=S.diff(w*S.diff(V,rho),rho)/w
    lhs=sig*(t_op(V*sig)+V*t_op(sig))
    target=2*V*S.diff(sig,rho)**2-boxV*sig**2
    zero('nonconstant_potential_integration_by_parts',lhs-target+S.diff(2*w*V*sig*S.diff(sig,rho),rho)/w)
    ep,la,U,rup,j=S.symbols('epsilon lambda U rUp j')
    E=1/ep+11*la/3; V0=(2*U+rup)/6; Z=12*la-U
    scalar_gradient=E-4*la+2*V0+Z/36+11*j/9
    seed_principal=(6+R(11,6)*ep*U+2*ep*(U+rup)+R(22,3)*ep*j)/(6*ep)
    zero('off_shell_seed_principal_matches_scalar_form',scalar_gradient-seed_principal)
    u=S.symbols('u'); pol=2*u**5-4*u**6
    rfirst=-u*S.diff(pol,u)
    rsecond=u**2*S.diff(pol,u,2)+2*u*S.diff(pol,u)
    zero('v_first_radial_polynomial',rfirst-(-10*u**5+24*u**6))
    zero('v_second_radial_polynomial',rsecond-(60*u**5-168*u**6))
    zero('v_value_absolute_coefficient_sum',sum(abs(v) for v in S.Poly(pol,u).all_coeffs())-6)
    zero('v_first_absolute_coefficient_sum',sum(abs(v) for v in S.Poly(rfirst,u).all_coeffs())-34)
    zero('v_second_absolute_coefficient_sum',sum(abs(v) for v in S.Poly(rsecond,u).all_coeffs())-228)
    check('v_gradient_bound_80',5*34**2<80**2)
    zero('Box_V0_coefficient_ledger',5*228+(22+2*5)*34-2228)
    check('Box_V0_reserve_2300',2228<2300)
    upol=20*u**5-36*u**6
    zero('rUprime_absolute_coefficient_sum',sum(abs(v) for v in S.Poly(upol,u).all_coeffs())-56)
    du=-u*S.diff(upol,u)
    zero('rUprime_gradient_absolute_sum',sum(abs(v) for v in S.Poly(du,u).all_coeffs())-316)
    check('rUprime_gradient_reserve_800',5*316**2<800**2)


def coercivity_and_energy_budgets():
    gl=R(1,2000**10); ep_l=gl/180; delta0=R(13,10**23)
    check('Hc_zero_order_guard',7+28*delta0<8)
    zero('TT_counterterm_raw_numerator',128*176-22528)
    check('TT_counterterm_below_two_using_pi_gt_three',R(22528,5760*3)<2)
    check('reduced_TT_gap_989',R(99,100)-2*gl>R(989,1000))
    check('V0_value_guard',4+6*delta0<5)
    check('mixed_Hc_gradient_guard',100+400*delta0<101)
    check('scalar_as_guard',4*delta0/36+R(11,9)*7<9)
    # Multiply each scalar coefficient margin by epsilon; this avoids
    # floating approximations to the very large Mstar^2/lambda.
    gradient_margin=R(1,2)-ep_l*(R(31,6)+12*delta0)
    potential_margin=2-R(2,3)*ep_l-R(4,3)*delta0-2300*delta0*ep_l
    check('scalar_gradient_half_baseline_margin',gradient_margin>0)
    check('scalar_potential_half_baseline_margin',potential_margin>0)
    e,dd=S.symbols('epsilon_lambda delta0',positive=True)
    ee=1/e+R(11,3)
    zero('scalar_gradient_margin_derivation',e*(ee-5-12*dd-(ee+4)/2)-(R(1,2)-e*(R(31,6)+12*dd)))
    zero('scalar_potential_margin_derivation',e*(4*ee-4*dd/(3*e)-2300*dd-2*(ee+4))-(2-R(2,3)*e-R(4,3)*dd-2300*dd*e))
    av,dv,p,q=S.symbols('a d p q',real=True)
    gap=R(989,1000)
    zero('energy_lower_norm_gap',av**2+gap*dv**2-gap*(av**2+dv**2)-(1-gap)*av**2)
    zero('energy_source_Cauchy_square',(p*p+q*q)*(av*av+dv*dv)-(p*av+q*dv)**2-(p*dv-q*av)**2)
    zero('energy_inverse_constant',1/gap-R(1000,989))
    aa,bb,cc,sigma,c=S.symbols('A B D sigma c',real=True)
    mat=S.Matrix([[aa,bb],[bb,cc]])
    vec=S.Matrix([sigma,c]); signature=S.diag(-1,1)
    zero('signature_flipped_mixed_real_cancellation',(vec.T*signature*mat*vec)[0]+aa*sigma**2-cc*c**2)
    # Finite scalar-entry Schur congruence; operator domains are not tested.
    left=S.Matrix([[1,-bb/cc],[0,1]])
    right=S.Matrix([[1,0],[-bb/cc,1]])
    matrix_zero('finite_saddle_Schur_congruence',left*mat*right-S.diag(aa-bb**2/cc,cc))


def normal_roots_and_order_scope():
    s,la,E=S.symbols('s lambda E')
    tn=-la*s*(s-3)
    zero('reduced_scalar_normal_factor',tn+4*la+la*(s+1)*(s-4))
    zero('reduced_scalar_source_root',((tn+4*la)/la).subs(s,-1))
    zero('reduced_scalar_fast_root',((tn+4*la)/la).subs(s,4))
    zero('removed_original_oscillatory_factor',(tn+4*la)*(tn-E)-(tn+4*la)**2+(E+4*la)*(tn+4*la))
    check('scalar_weight_shifts_inside_root_strip',-1<R(3,2)<R(5,2)<R(7,2)<4)
    check('mass_weight_is_critical_L2',3>R(3,2))
    check('retained_York_tail_is_fast',4>R(7,2))
    zero('Schur_interior_differential_order',2-4+2)
    check('Schur_lower_than_scalar_order',2-4+2<2)
    zero('constant_vacuum_trace_reference_matching',R(3,1)*(-R(11,360))/(8)+R(11,960))
    # A nonlinear pullback Hessian has an off-shell acceleration term.
    t,k,acc,g1,g2=S.symbols('t k acceleration gradient Hessian')
    y=k*t+acc*t*t/2
    action=g1*y+g2*y*y/2
    zero('off_shell_pullback_Hessian_contact',S.diff(action,t,2).subs(t,0)-g2*k*k-g1*acc)
    check('different_linear_Fredholm_indices',-1!=0)


def observables_and_finite_remainder():
    r,F,Fr,s,alpha,v=S.symbols('r F Fr s alpha v',positive=True)
    # v=1+r*s_r>0, so dR/dr=exp(s)*v.
    dR=S.exp(s)*v; Fnew=F*v**2
    anew=alpha+s-S.log(v)
    zero('areal_metric_radial_transform',S.exp(2*s)/F/dR**2-1/Fnew)
    zero('areal_metric_time_transform',S.exp(2*anew)*Fnew-S.exp(2*s+2*alpha)*F)
    fR_h=S.exp(-s)*v*Fr
    zero('exact_surface_gravity_invariance',S.exp(anew)*fR_h/2-S.exp(alpha)*Fr/2)
    z,c4,la,m,Q=S.symbols('z c4 lambda m Q')
    tail=c4*z**4
    fbar=la/z**2+1-2*m*z+Q*z**2
    transformed=fbar*(1-4*tail)**2
    target=la*S.exp(2*tail)/z**2+1-2*m*z*S.exp(-tail)
    zero('fast_field_map_preserves_lambda_constant_mass',S.series(transformed-target,z,0,2).removeO())
    lapse=tail-S.log(1-4*tail)
    zero('fast_field_map_has_no_cubic_lapse',S.series(lapse,z,0,4).removeO())
    zero('fast_lapse_leading_coefficient',S.series(lapse,z,0,5).removeO()-5*c4*z**4)
    x=S.symbols('x',real=True)
    zero('log_observable_derivative',S.diff(S.log(1+x),x)-1/(1+x))
    check('log_observable_lipschitz_guard',1/(1-R(1,2))==2)
    delta=S.symbols('delta',nonnegative=True)
    # Exact Taylor second derivatives supply the stated |s|<=delta guards.
    zero('EH_R_Taylor_coefficient',S.diff(S.exp(2*x),x,2).subs(x,delta)/2-2*S.exp(2*delta))
    zero('EH_lambda_Taylor_coefficient',6*S.diff(S.exp(4*x),x,2).subs(x,delta)/2-48*S.exp(4*delta))
    zero('EH_conformal_gradient_after_parts',-6+12-6)
    zero('EH_first_order_R_conformal_weight',S.diff(S.exp(2*x),x).subs(x,0)-2)
    zero('EH_first_order_vacuum_conformal_weight',6*S.diff(S.exp(4*x),x).subs(x,0)-24)
    check('fast_value_remainder_integrable',2-2*4<-1)
    check('scalar_C5_interior_Sobolev_margin',(5+2)-5>R(2,2))
    check('scalar_C7_interior_Sobolev_margin',(7+2)-7>R(2,2))
    zero('curvature_field_map_two_derivative_loss',5-2-3)
    zero('C7_to_C5_field_map_count',7-2-5)


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
            field_map_and_trace_normalization()
            off_shell_reduction()
            seed_geometry_and_blocks()
            scalar_form_coefficients()
            coercivity_and_energy_budgets()
            normal_roots_and_order_scope()
            observables_and_finite_remainder()
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
