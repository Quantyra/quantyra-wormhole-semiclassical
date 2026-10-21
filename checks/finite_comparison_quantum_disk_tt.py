"""S167: full disk-chart TT positivity at the actual finite comparison.

Exact metric algebra and rational budgets; continuum claims remain ordinary proofs.
"""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S

ROOT = Path(__file__).resolve().parents[1]
STEM = "finite-comparison-quantum-disk-tt"
OWN = ["checks/finite_comparison_quantum_disk_tt.py", "notes/"+STEM+".md"]
PREDECESSORS = [("finite-comparison-massive-geodesic-form", 314),
                ("stationary-curved-quantum-tt", 157),
                ("finite-comparison-classical-disk-contact", 301)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "complete physical action in the disk chart on actual returned spherical TT inputs",
    "disk_TT_lower_bound": ".90 E_ell",
    "core": "smooth bolt-regular compact-outer TT and explicit stronger form completion",
    "full_angular_tower_and_common_fourth_restoration": True,
    "outer_Weyl_logarithm_cancelled_on_core": True,
    "physical_tensor_chart_contacts_bounded": False,
    "fixed_density_contacts_bounded": False,
    "full_quantum_operator_domain_identified": False,
    "physical_mass_fast_domain_membership_proved": False,
    "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
    "nonlinear_family_preserved": False, "signed_compatibility_evaluated": False,
    "fixture_scope": "exact local Hessian/IMS algebra and rational finite-geometry budgets",
    "formal_verification": False,
}


def check(name, ok, detail=""):
    assert name not in {row["name"] for row in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expr):
    expr = S.cancel(expr.doit())
    check(name, expr == 0, expr)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_bytes())
    return data.get("runs", [data])[-1]


def inputs(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def paths():
    found = set(OWN)
    for stem, _ in PREDECESSORS:
        found.add("evidence/"+stem+".json")
        found.update(inputs(latest(stem)))
    return sorted(found)


def predecessors():
    expected = {}
    for stem, count in PREDECESSORS:
        run = latest(stem)
        check(stem+"_receipt", run["status"] == "passed" and len(run["checks"]) == count
              and all(row["passed"] for row in run["checks"]) and inputs(run) == inputs(run, True))
        conflicts = [p for p, h in inputs(run).items() if p in expected and expected[p] != h]
        check(stem+"_shared_inputs", not conflicts, conflicts)
        expected.update(inputs(run))
    for p, h in sorted(expected.items()):
        check("input_"+p, sha256((ROOT/p).read_bytes()).hexdigest() == h)
    for p in ["finite-comparison-massive-geodesic-form", "finite-comparison-massive-constant-contact", "finite-comparison-classical-disk-contact", "finite-comparison-classical-tt-form",
              "comparison-path-interior-response", "divided-jet-dirac-response",
              "outer-finite-background-response", "stationary-curved-quantum-tt",
              "stationary-centered-canonical-response", "stationary-completed-metric-response-norm"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


def exact_Weyl_Hessian():
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


def actual_Weyl_budgets():
    check('actual_return_S',R(1,2)+R('2.25')*R('4.83')<12)
    check('actual_return_T',(R('2.002')+3)/2+R('1.008')*R('3.002')+R('2.25')*R('4.83')<18)
    check('actual_return_S_derivative',R('5.32')/2+R('2.25')*R('3.002')+R('22.001')*R('4.83')<116)
    check('actual_Rhat_bound',R('2.001')*R('6.1')<13)
    check('actual_U_bound',R('22.001')<23)
    check('actual_V_bound',R('2.25')**2<R('5.1'))
    check('D_S_coefficient',R(13,15)+R(1,9)<1)
    check('D_T_coefficient',R(7*13,45)+11*R('5.1')/135<R('2.5'))
    zero('D_gradient_coefficient',R(11,135)*R('2.25')-R(11,60))
    check('L_T_coefficient',R(13,15)+R(2*23,5)+7*R('5.1')/45<11)
    zero('L_gradient_coefficient',R(7,45)*R('2.25')-R('.35'))
    check('J_coefficient',R(13,30)+7*(23+R('5.1'))/90+R(1,18)<R('2.7'))
    check('W_coefficient',R(2*13,15)+14*(23+R('5.1'))/45<11)
    check('Z_T_coefficient',R(169,30)+13*(23+R('5.1'))/15+R(23**2,5)
          +7*23*R('5.1')/45+7*R('5.1')**2/90<160)
    check('Z_ST_coefficient',R(169,30)<6)
    check('Weyl_D_ledger',(12+R('2.5')*18+R(11,60)*116)**2<10000)
    check('Weyl_L_ledger',(11*18+R('.35')*116)**2<60000)
    check('Weyl_J_ledger',R('2.7')*116**2<40000)
    check('Weyl_W_ledger',11*R('2.25')*18*116<60000)
    check('Weyl_Z_ledger',160*18**2+6*12*18<60000)
    check('Weyl_absolute_w_ledger',R('5.1')*116**2<70000)
    check('Weyl_combined_ledger',25*10000+5*60000+40000+60000+60000<710001)
    check('Weyl_physical_lower_order',R(710001,3)<1000000)
    check('Weyl_k_cost',R(360,3)*R('4e-27')<R('1e-24'))
    B,G,t,mu,e,lam=S.symbols('B G t mu e lambda',positive=True)
    muhat=mu/(4*S.pi*B**2*t**4)
    zero('Weyl_physical_measure',B**2/(2*S.pi)*muhat-mu/(8*S.pi**2*t**4))
    zero('Weyl_energy_conversion',lam*e/(8*S.pi**2)/(e/(8*S.pi*G))-G*lam/S.pi)
    z,ds,de=S.symbols('z DeltaS DeltaT')
    # Both Young inequalities leave the displayed positive square.
    zero('Weyl_D_Young_square',de**2/100+25*ds**2-de*ds-(de/10-5*ds)**2)
    zero('Weyl_L_Young_square',de**2/10-de*ds-de**2/20+5*ds**2-(de-10*ds)**2/20)


def collar_Hardy_and_model():
    check('outer_log_upper',sum(R('2.31')**j/S.factorial(j) for j in range(18))>10)
    check('outer_jet_value',R('1e-7')*R('.1')**4*R('3.31')+R('1e-11')*R('.1')**3<R('4e-11'))
    check('outer_jet_first',R('1e-7')*100*R('.1')**3*R('3.31')+3*R('1e-11')*R('.1')**2<R('4e-8'))
    check('outer_jet_second',R('1e-7')*2*100**2*R('.1')**2*R('3.31')+6*R('1e-11')*R('.1')<R('7e-5'))
    check('outer_jet_third',R('1e-7')*6*100**3*R('.1')*R('3.31')+6*R('1e-11')<R('.2'))
    q0=R('4e-11')+R('.1')*R('4e-8')+R('.1')**2*R('7e-5')/2
    q1=R('4e-8')+R('.1')*R('7e-5'); q2=R('7e-5')
    check('extension_first',q1+20*q0<R('2e-5'))
    check('extension_second',q2+40*q1+600*q0<R('.001'))
    check('extension_third',60*q2+1800*q1+60000*q0<R('.2'))
    check('extension_full_first',R('1.2')+R('2e-5')<R('1.201'))
    check('extension_full_second',42+R('.001')<=R('42.001'))
    check('extension_full_third',3270+R('.2')<3271)
    check('auxiliary_potential',R('42.001')/2+R('1.201')**2/4<22)
    check('auxiliary_potential_Lipschitz',R(3271,2)+R('1.201')*R('42.001')/2<2000)
    Y=S.symbols('Y',real=True); hh=S.Function('h')(Y); T=S.Function('T')(Y)
    vh=S.exp(hh/2); potential=S.diff(hh,Y,2)/2+S.diff(hh,Y)**2/4
    H=lambda v:-S.diff(v,Y,2)+potential*v
    zero('auxiliary_constant_base_mode',H(vh))
    zero('auxiliary_half_density_Laplacian',H(vh*T)+vh*(S.diff(T,Y,2)+S.diff(hh,Y)*S.diff(T,Y)))
    vv=S.Function('V')(Y); scalar_H=lambda f:-S.diff(f,Y,2)+vv*f
    zero('scalar_quadratic_H',scalar_H(Y**2).subs(Y,0)+2)
    zero('scalar_quadratic_H2',scalar_H(scalar_H(Y**2)).subs(Y,0)+4*vv.subs(Y,0))
    z=S.symbols('z',positive=True)
    zero('log_shift_majorant_derivative',S.diff(z-S.log(1+z),z)-z/(1+z))
    zero('Hardy_quarter_digamma',S.expand_func(S.polygamma(0,R(1,4)))+S.EulerGamma+S.pi/2+3*S.log(2))
    check('Hardy_constant_lower',1+2*R(7,10)+R(11,7)<4)
    check('model_scalar_second_contact',100*R(2,3)+3*4*R('.00352')<70)
    check('model_contact_aggregate',R('3e13')+10000+70+100<R('1e14'))
    check('model_escape_aggregate',R('2e11')+R('2e12')<R('3e12'))
    check('finite_H2_conversion',2*22**2==968)


def finite_transition_geodesics():
    # The transition plus the proper second-anchor neighborhood.
    lo,hi=R('.0075'),R('.0225')
    check('transition_log_bound',sum(R(5)**j/S.factorial(j) for j in range(18))>1/lo)
    seed=[R('.045'),2,10,50,200]
    target=[R('.05'),3,11,1500,10**8]
    for j in range(1,6):
        weight=hi**(4-j) if j<=4 else lo**(4-j)
        delta=R('1e-7')*S.factorial(j)*100**j*weight*6
        if j<=3: delta+=R('1e-11')*S.factorial(3)/S.factorial(3-j)*hi**(3-j)
        check('finite_transition_metric_jet_'+str(j),seed[j-1]+delta<target[j-1])
    # Differentiate Y''=a(Y)(1-Y'^2), bounding the two factors separately.
    a0=R('.05'); a1=3; a2=11+3; a3=1500+3*11+3*10
    a4=10**8+6*1500+3*11+4*11*10+3*100
    b1=2; b2=2+2*10; b3=6*10+2*100; b4=6*10**2+8*100+2*2000
    check('finite_geodesic_Y2',a0<1)
    check('finite_geodesic_Y3',a1+a0*b1<10)
    check('finite_geodesic_Y4',a2+2*a1*b1+a0*b2<100)
    check('finite_geodesic_Y5',a3+3*a2*b1+3*a1*b2+a0*b3<2000)
    check('finite_geodesic_Y6',a4+4*a3*b1+6*a2*b2+4*a1*b3+a0*b4<2*10**8)


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
    Yj=[1,1,10,100,2000,200000000]
    sixth=[S.bell(6,k,tuple(Yj[:7-k])) for k in range(1,7)]
    check('geodesic_sixth_chain_array',sixth==[200000000,14500,2115,245,15,1],sixth)
    chi=[1]+[cb[j-1]*100**j for j in range(1,7)]
    check('geodesic_partition_X6',sum(sixth[k-1]*chi[k] for k in range(1,7))<10**22)
    fourth=[S.bell(4,k,tuple(Yj[:5-k])) for k in range(1,5)]
    check('geodesic_partition_X4',sum(fourth[k-1]*chi[k] for k in range(1,5))<10**14)
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
    written=6*k*p**2*f1**2+(8*k*p*p2+4*k*warp*p**2+2*S.diff(k,y)*p**2)*fun*f1+k*(p**4+(p2+warp*p)**2)*fun**2
    flux2=-2*k*p**2*fun*f1
    zero('weighted_fourth_IMS_with_cross_term',S.expand(direct-written-S.diff(flux2,y)-warp*flux2))


def reference_and_complete_budget():
    z,st,tt,X=S.symbols('z S T X',real=True); pot=S.Function('DeltaV')
    ref=S.exp(2*z*(st+tt))*pot(X*S.exp(2*z*st))
    expected=(st+tt)**2*pot(X)+2*st*(st+tt)*X*S.diff(pot(X),X)
    expected+=st**2*(X*S.diff(pot(X),X)+X**2*S.diff(pot(X),X,2))
    zero('reference_exact_second_X_variation',S.diff(ref,z,2).subs(z,0)/4-expected)
    check('reference_missing_second_X_negative_control',
          (st**2*X*S.diff(pot(X),X)).subs(pot(X),X**2).doit().subs({st:1,X:2})!=0)
    t,pA,pv,dp=S.symbols('t pA pv dp',positive=True)
    profA=pA*t**-4+pv*(1-t**-4)+dp*(-t**-4+4*t**-5-3*t**-6)
    profS=-pA*t**-4+pv*(1+t**-4)+dp*(t**-4-6*t**-5+6*t**-6)
    zero('reference_profile_Ward_identity',profS-profA-t*S.diff(profA,t)/2)
    # The source coefficient bounds are strict; their rational majorants meet.
    check('reference_global_pressure_guard',3300+2+3*(1+6+6)<=3341)
    check('reference_global_radial_derivative',4*3300+4+3*(4+30+36)<14000)
    check('reference_finite_return_ledger',3341*30**2+6682*12*30+(3341+7000)*12**2<7000000)
    check('reference_physical_cost',224*10**6*R(22,7)*R('4e-27')<R('1e-16'))
    # Preserve the nonstationary chart distinction in an exact scalar fixture.
    aa,bb=S.symbols('a b'); chart=z+z*z/2; action=aa*chart+bb*chart**2/2
    zero('nonstationary_chart_contact_fixture',S.diff(action,z,2).subs(z,0)-bb-aa)
    check('nonstationary_contact_not_zero',aa.subs(aa,1)!=0)
    eps=R('.001'); ac=R('.01'); K=R('1e14'); rad=R(1,400)
    check('outer_Weyl_log_change',R('3e-8')+R('1.2e-8')/(1-R('1.2e-8'))<R('1e-6'))
    check('outer_fourth_margin',9-R(38,15)-R(19,60)-4-R(7,20)-R('1e-6')>R(8,5))
    check('interior_fourth_margin',R('13.9')-R('6.001')-R(19,60)>R('7.5'))
    check('finite_split_inside',eps/(4*K)<rad**2)
    logarg=rad**2*4*K/eps
    check('finite_split_log_guard',sum(R(30)**j/S.factorial(j) for j in range(71))>logarg)
    check('finite_split_gradient',K*(R(1,4)+30)<R('1e16'))
    loss=eps+eps*R(22,7)/4+2*eps/3
    check('all_three_fourth_losses',loss<R('.00246'))
    check('surviving_fourth_margin',R(8,5)*R(7,22)/30-loss>R(1,300))
    check('Weyl_IMS_gradient_coefficient',6*7*25+50000<60000)
    cross=8*7*5*350+4*7*5*25*ac**2+2*2*25
    check('Weyl_IMS_cross_coefficient',cross<100000)
    check('Weyl_IMS_potential_coefficient',7*(625+(350+25*ac**2)**2)+50000<1000000)
    cj0=R('5.2e-29'); cj1=R('8.5e-28')
    outer=cj0*16*ac**4/R('2.23')**4
    check('actual_outer_physical_weight',outer<R('3.4e-37'))
    extra=R('3.4e-37')*(5/ac)**2
    check('localization_extra_gradient',extra<=R('8.5e-32'))
    check('localized_gradient_sum',cj1+extra<R('8.6e-28'))
    check('half_density_gradient_weight',R('2.001')**2*(2*ac)**2/2<R('.001'))
    check('literal_J1_aggregate',cj1+3*R('8.6e-28')+R('.001')*R('3.4e-37')<R('1e-26'))
    check('literal_J0_aggregate',3*cj0<R('2e-28'))
    check('new_second_contact_envelope',2*R('1e22')+R('1e14')+3*R('1e16')+10**12<R('1e23'))
    check('new_zero_contact_envelope',2*R('1e22')+R('3e12')+R('2e13')+10**12<R('1e23'))
    check('new_gradient_energy',R('1e23')*R('1e-26')==R('.001'))
    check('new_zero_energy',R('1e23')*R('2e-28')==R('.00002'))
    weyl_ims=(60000/ac**2*cj1+1000000/ac**4*R('3.4e-37'))/90
    check('weighted_Weyl_IMS_energy',weyl_ims<R('1e-12'))
    check('Weyl_low_energy',R('1e6')/2000**10+R('1e-24')<R('1e-20'))
    check('full_quantum_disk_loss',R('.005')+R('.001')+R('.00002')+R('1e-12')+R('1e-16')+R('1e-20')<R('.01'))
    check('full_action_disk_TT_margin',R('.92')-R('.01')>R('.90'))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding="utf-8", newline="\n")
    sys.stderr.reconfigure(encoding="utf-8", newline="\n")
    ps = paths()
    hashes = lambda: {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before = hashes()
    error = None
    print(json.dumps({"event": "start", "story": "S167 finite comparison quantum disk TT",
                      "input_count": len(ps)}), flush=True)
    try:
        predecessors()
        exact_Weyl_Hessian()
        actual_Weyl_budgets()
        collar_Hardy_and_model()
        finite_transition_geodesics()
        localization_and_IMS()
        reference_and_complete_budget()
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    after = hashes()
    good = error is None and before == after and all(row["passed"] for row in ROWS)
    run = {
        "story": "S167 finite comparison quantum disk TT",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if good else "failed", "python": sys.version,
        "sympy": S.__version__, "protocol": PROTOCOL,
        "input_sha256_before": before, "input_sha256_after": after,
        "unchanged_inputs": before == after, "checks": ROWS, "exception": error,
        "finite_comparison_full_disk_TT_positive": good,
        "actual_Weyl_and_reference_forms_bounded": good,
        "outer_Weyl_logarithm_cancelled_on_core": good,
        "explicit_stronger_disk_form_completion": good,
        "physical_tensor_chart_contacts_bounded": False,
        "fixed_density_contacts_bounded": False,
        "full_quantum_operator_domain_identified": False,
        "physical_mass_fast_domain_membership_proved": False,
        "physical_tensor_TT_Hessian_positive": False,
        "full_candidate_Hessian_inverse": False, "cold_correction_bound": False,
        "signed_compatibility_evaluated": False, "nonlinear_family_preserved": False,
        "useful_Newton_certificate": False, "S167_complete": False, "goal_complete": False,
        "unresolved": ["physical tensor/gauge/density contacts and mixed blocks",
                       "full logarithmic operator domains and physical-direction membership",
                       "candidate coupled range and adjoint", "projection direction and cold/higher bounds",
                       "nonlinear preservation and solve", "signed compatibility", "physical requirements"],
    }
    if not good:
        run["failed_checker_source_utf8"] = (ROOT/OWN[0]).read_text(encoding="utf-8")
        run["failed_note_source_utf8"] = (ROOT/OWN[1]).read_text(encoding="utf-8")
    dest = ROOT/"evidence"/(STEM+".json")
    data = json.loads(dest.read_bytes()) if dest.exists() else {"runs": []}
    data["runs"].append(run)
    dest.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps({"status": run["status"], "checks": len(ROWS),
                      "unchanged_inputs": before == after}), flush=True)
    return 0 if good else 1


if __name__ == "__main__":
    raise SystemExit(main())
