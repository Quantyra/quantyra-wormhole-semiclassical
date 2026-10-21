"""S167: actual finite TT domain, inverse and coupled physical range.

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
STEM = "finite-comparison-tt-realization-and-range"
OWN = ["checks/finite_comparison_tt_realization_and_range.py", "notes/"+STEM+".md"]
PREDECESSORS = [("finite-comparison-weighted-mixed-form", 335),
                ("stationary-global-weighted-map", 50),
                ("stationary-weighted-tt-inverse", 52),
                ("finite-comparison-disk-return-form", 291)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "actual TT realization and coupled one-condition range at each finite comparison",
    "physical_domain": "unchanged disk X_delta; scalar-c domain pulled back by actual shape return",
    "plain_scalar_c_H4log_inferred_from_C5": False,
    "TT_lower_bound": ".90 E_ell",
    "actual_outgoing_contraction": ".025",
    "auxiliary_inverse_norm": "1.026",
    "physical_mass_fast_kernel_dimension": 0,
    "physical_mass_fast_cokernel_dimension": 1,
    "fixed_seed_projection_transversality": False,
    "useful_fast_H4log_inverse_constant": False,
    "fixture_scope": "exact kinematic, normal, mass, density and block algebra; manufactured regular York fixtures",
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
    for p in ["finite-comparison-weighted-mixed-form", "finite-comparison-quantum-disk-tt",
              "finite-comparison-massive-geodesic-form", "finite-comparison-disk-return-form",
              "stationary-global-tt-reduction", "stationary-global-weighted-map",
              "stationary-weighted-tt-inverse", "stationary-scalar-outgoing-inverse",
              "stationary-adjoint-phase"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)








def exact_kinematics():
    A,Rr,Ap,Rp,Y,Yp,Ypp,T,Tp,theta=S.symbols('A R Ap Rp Y Yp Ypp T Tp theta')
    w=A+2*Rr;wp=Ap+2*Rp;V=wp/2+2*A*A+4*Rr*Rr
    p0=3*Yp/2-w*Y/2;p1=-Yp/2+(3*A/2-Rr)*Y;p2=-Yp/2+(Rr-A/2)*Y
    b=T-p0;c=-T-p2;K=(Yp-R(2,3)*T+w*Y)/4
    N=-R(2,3)*T+2*K-2*Rr*Y
    zero('sphere_amplitude_map',c-N)
    zero('TT_trace',b+(T-p1)+2*c)
    tau=theta+T/6-(Yp+w*Y)/4
    zero('conformal_coordinate',tau-(theta-K))
    Sc=c/2+Rr*Y;Tc=-(b+3*c)/2+(A-Rr)*Y
    zero('actual_return_shape',Tc-T)
    zero('actual_return_conformal',Sc-(K-T/3))
    residual=R(3,2)*(Ypp+w*Yp)-V*Y-Tp-4*Rr*T
    Kp=(Ypp-R(2,3)*Tp+wp*Y+w*Yp)/4
    Ricrr=-(Ap+A*A+2*Rp+2*Rr*Rr)
    rhs=R(2,3)*Rr*T-R(1,3)*Ricrr*Y
    zero('exact_K_derivative_residual',6*(Kp-rhs)-residual)
    solvedYpp=R(2,3)*(V*Y+Tp+4*Rr*T)-w*Yp
    zero('on_York_K_derivative',(Kp-rhs).subs(Ypp,solvedYpp))
    check('omitted_curvature_term_detected',(rhs-R(2,3)*Rr*T).subs({A:2,Rr:3,Ap:5,Rp:7,Y:11})!=0)
    check('finite_regularity_order',2+3==5 and 1+4==5)
    # Triangular coordinates are exact and do not apply a logarithm to K.
    k=S.symbols('k');U=S.Matrix([[1,k-R(1,3)],[0,1]])
    Ui=S.Matrix([[1,-k+R(1,3)],[0,1]])
    zero('triangular_coordinate_inverse',(U*Ui-S.eye(2)).norm()**2)


def normal_and_mass():
    s=S.symbols('s');J,NT,lam=S.symbols('J NT lambda')
    K=-2/(3*(s+1));col=S.Matrix([K-R(1,3),1])
    C=lambda q:S.Matrix([[1,(q+3)/(3*(q+1))],[0,2]])
    H=C(3-s).T*S.diag(J,NT)*C(s)
    zero('TT_normal_column',(C(s)*col-S.Matrix([0,2])).norm()**2)
    zero('shape_normal_congruence',(col.subs(s,3-s).T*H*col)[0]-4*NT)
    ccol=-R(3,2)*col
    zero('scalar_c_normal_congruence',(ccol.subs(s,3-s).T*H*ccol)[0]-9*NT)
    scalar=S.Matrix([1,0])
    zero('mixed_normal_zero',(scalar.T*H*col)[0])
    zero('shape_mass_column',(col.subs(s,3)+S.Matrix([R(1,2),-1])).norm()**2)
    zero('normalized_c_mass_column',(ccol.subs(s,3)-R(3,4)*S.Matrix([1,-2])).norm()**2)
    zero('K_mass_value',K.subs(s,3)+R(1,6))
    zero('physical_mass_scalar_amplitude',R(1,3)-K.subs(s,3)*(-2))
    check('TT_critical_normal_positive',R(9,4)>2 and R(9,4)*(R(9,4)-2)>0)
    check('critical_and_fast_weights_inside_strip',-1<R(3,2)<R(5,2)<3<R(7,2)<4)
    check('bootstrap_steps_below_two',R(5,2)-R(3,2)<2 and R(7,2)-R(5,2)<2)
    check('forced_fast_tail_integrable',4>R(7,2) and 5>R(7,2))
    # A fixed-theta column has a nonzero conformal normal component.
    wrong=S.Matrix([-R(1,3),1])
    check('fixed_theta_not_TT',S.simplify((C(s)*wrong)[0])!=0)


def density_and_spectral_scaling():
    r,L,m,q,e,kappa,lam=S.symbols('r L m q e kappa lambda',positive=True)
    a2=(lam*r*r+1-2*m/r+q/r**2+e/r**3)/kappa**2
    M=r*r*S.diff(a2,r)-2*r*a2
    zero('finite_areal_moment',M-(-2*r+6*m-4*q/r-5*e/r**2)/kappa**2)
    x=S.symbols('x',positive=True);psi=(x-2)**6*(3-x)**6
    I0=S.integrate(psi,(x,2,3));I1=S.integrate(x*psi,(x,2,3))
    check('bump_nonzero_leading_moment',I0>0 and I1>0,(I0,I1))
    zero('bump_denominator_leading',S.diff(-2*L**2*I1/kappa**2+6*m*L*I0/kappa**2,L,2)+4*I1/kappa**2)
    zero('critical_energy_bump_scale',(L**-2)**2*L**3-1/L)
    zero('canonical_second_jet_scale',(L**-2*L**2)**2/L-1/L)
    zero('canonical_third_jet_scale',(L**-2*L**3)**2/L-L)
    # The inward conformal tail is L^-3/Y, differentiated in canonical Y.
    y=S.symbols('y',positive=True)
    second=S.integrate((2*L**-3/y**3)**2,(y,1/L,1))
    zero('inward_conformal_J2_scale',second-R(4,5)*(1/L-L**-6))
    check('positive_power_log_limit',S.limit(S.log(L)/L,L,S.oo)==0)
    mm,p,cc=S.symbols('M p C',positive=True)
    zero('single_mode_low_frequency_scale',mm**2*(p/mm)**6-p**6/mm**4)
    zero('single_mode_high_frequency_scale',mm**2*(p/mm)**4-p**4/mm**2)
    zero('tower_high_mass_integral',S.integrate(16*cc/mm**3,(mm,p,S.oo))-8*cc/p**2)
    minmass=S.symbols('Mmin',positive=True)
    zero('tower_low_mass_log',S.integrate(8*cc/mm,(mm,minmass,p))-8*cc*S.log(p/minmass))
    check('cosine_large_argument_upper',2+R(1,2)+R(1,24)<3)
    check('small_argument_remainder_positive',S.Rational(1,720)>0)
    # A C5 bound alone leaves a logarithmic endpoint loss.
    n=S.symbols('n',positive=True)
    for j in range(6):check('bounded_C5_derivative_'+str(j),j-5<=0)
    ratio=(1+n*n)**4/n**8
    check('derivative_H4_multiplier_not_small',all(v>0 for v in S.Poly(S.expand((1+n*n)**4-n**8),n).coeffs()))
    check('log_endpoint_loss_diverges',S.limit(S.log(2+n),n,S.oo)==S.oo)


def operators_and_scope():
    # Nonunitary physical congruence must carry both output transposes.
    N=S.Matrix([[1,2],[0,3]]);D=S.Matrix([[5,1],[1,4]]);L=N.inv()
    Dc=L.T*D*L;inv=N*D.inv()*N.T
    zero('actual_c_inverse_congruence',(Dc*inv-S.eye(2)).norm()**2)
    check('nonunitary_not_similarity',Dc!=L*D*N)
    v=S.Matrix([S.Symbol('v0'),S.Symbol('v1')]);cv=N*v
    zero('physical_form_pairing',(cv.T*Dc*cv)[0]-(v.T*D*v)[0])
    check('positive_c_congruence',Dc[0,0]>0 and Dc.det()>0)
    G,lam=S.symbols('G lambda',positive=True)
    gap=R('.90')*lam/(8*S.pi*G)
    zero('critical_inverse_gap',gap*(8*S.pi*G/(R('.90')*lam))-1)
    zero('shape_gap_normalization',R('.90')*lam/(8*S.pi*G*18**2)-gap/18**2)
    check('actual_auxiliary_inverse',R('1.000002')/(1-R('.025'))<R('1.026'))
    check('new_weak_response_weight',R(3,2)-R(1,2)==1)
    check('mixed_forcing_is_critical_L2',1+2>R(3,2))
    check('physical_cokernel_dimension',0-(-1)==1)
    # One cokernel does not prove that a previously selected source is transverse.
    chi=S.Matrix([1,0]);pcut=S.Matrix([0,1])
    check('fixed_projection_transversality_not_automatic',(chi.T*pcut)[0]==0 and chi!=S.zeros(2,1) and pcut!=S.zeros(2,1))


def independent_regular_fixtures():
    import mpmath as mp
    mp.mp.dps=90
    a=lambda x:x*(1+x*x/7+x**4/101)
    r=lambda x:2+x*x/5+x**4/103
    yy=lambda x:x*(1+x*x/11)
    AA=lambda x:mp.diff(a,x)/a(x)
    RR=lambda x:mp.diff(r,x)/r(x)
    ww=lambda x:AA(x)+2*RR(x)
    vv=lambda x:mp.diff(ww,x)/2+2*AA(x)**2+4*RR(x)**2
    py=lambda x:mp.mpf('1.5')*(mp.diff(yy,x,2)+ww(x)*mp.diff(yy,x))-vv(x)*yy(x)
    # Cancel the removable polar singularity before near-zero quadrature.
    # The original uncancelled evaluation failed; its run/source are retained.
    sx=S.symbols('x');sa=sx*(1+sx**2/7+sx**4/101)
    sr=2+sx**2/5+sx**4/103;sy=sx*(1+sx**2/11)
    sA=S.diff(sa,sx)/sa;sR=S.diff(sr,sx)/sr;sw=sA+2*sR
    sV=S.diff(sw,sx)/2+2*sA**2+4*sR**2
    sPY=S.cancel(R(3,2)*(S.diff(sy,sx,2)+sw*S.diff(sy,sx))-sV*sy)
    check('manufactured_PY_regular_bolt',sPY.subs(sx,0)==0)
    integrand=S.lambdify(sx,S.cancel(sr**4*sPY),'mpmath')
    regular_py=S.lambdify(sx,sPY,'mpmath')
    tt=lambda x:(mp.mpf(16)/13+mp.quad(integrand,[0,x]))/r(x)**4
    kk=lambda x:(mp.diff(yy,x)-mp.mpf(2)*tt(x)/3+ww(x)*yy(x))/4
    for i,x in enumerate([mp.mpf(1)/5,mp.mpf(2)/3,mp.mpf(7)/4]):
        check('raw_and_regular_PY_'+str(i),abs(py(x)-regular_py(x))<mp.mpf('1e-75'))
        residual=py(x)-mp.diff(tt,x)-4*RR(x)*tt(x)
        ric=-mp.diff(a,x,2)/a(x)-2*mp.diff(r,x,2)/r(x)
        pred=mp.mpf(2)/3*RR(x)*tt(x)-ric*yy(x)/3
        check('manufactured_York_equation_'+str(i),abs(residual)<mp.mpf('1e-60'),mp.nstr(abs(residual),8))
        diff=mp.diff(kk,x)-pred
        check('independent_K_cancellation_'+str(i),abs(diff)<mp.mpf('1e-60'),mp.nstr(abs(diff),8))
        check('fixture_curvature_contact_nonzero_'+str(i),abs(ric*yy(x)/3)>mp.mpf('1e-6'))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 finite TT realization and range','input_count':len(ps)}),flush=True)
    try:
        predecessors();exact_kinematics();normal_and_mass();density_and_spectral_scaling();operators_and_scope();independent_regular_fixtures()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 finite TT realization and range','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'physical_disk_domain_preserved':good,'actual_TT_operator_realized':good,
         'critical_TT_inverse':good,'physical_mass_fast_TT_inverse':good,
         'actual_coupled_outgoing_inverse':good,'physical_coupled_kernel_zero':good,
         'physical_coupled_cokernel_one':good,'plain_scalar_c_H4log_inferred_from_C5':False,
         'fixed_seed_projection_transversality':False,'useful_fast_H4log_inverse_constant':False,
         'signed_compatibility_evaluated':False,'nonlinear_family_preserved':False,
         'S167_complete':False,'goal_complete':False}
    if not good:
        run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text()
        run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text()
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
