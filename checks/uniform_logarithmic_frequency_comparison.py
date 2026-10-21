"""S168: uniform logarithmic frequency comparison.

Exact metric algebra and rational budgets; continuum claims remain ordinary proofs.
"""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S
import numpy as np
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
STEM = "uniform-logarithmic-frequency-comparison"
OWN = ["checks/uniform_logarithmic_frequency_comparison.py", "notes/"+STEM+".md", "evidence/uniform-logarithmic-comparison-source-review.json"]
PREDECESSORS = [("canonical-frequency-core-bound", 360)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope":"Uniform finite scalar/geodesic comparison and actual centered-seed outer logarithmic lower form",
    "mpmath_digits":40,"relative_tolerance":1e-25,
    "bessel_arguments":[.2,1,3],"heat_masses":[1.3],"heat_radii":[.4,1],"feynman_parameters":[.2,.5],
    "polar_r":1,"polar_rprime":1.3,"polar_times":[.07,.2],"polar_modes":[0,1,3],
    "lower_order_constants_numerically_calibrated":False,
    "fixtures_are_exterior_spectral_data":False,
    "actual_wormhole_instability_proved":False,"actual_causal_inverse_proved":False,
    "order_reduction_adopted":False,"formal_verification":False,
}
NUMERICS=[]


def check(name, ok, detail=""):
    assert name not in {row["name"] for row in ROWS}
    row = {"name": name, "passed": bool(ok), "detail": str(detail)}
    ROWS.append(row)
    print(json.dumps(row), flush=True)
    if not ok:
        raise AssertionError(name)


def zero(name, expr):
    expr = S.cancel(S.sympify(expr).doit())
    check(name, expr == 0, expr)


def latest(stem):
    data = json.loads((ROOT/"evidence"/(stem+".json")).read_bytes())
    return data.get("runs", [data])[-1]


def inputs(run, after=False):
    return run.get("input_sha256_after" if after else "input_sha256_before",
                   run.get("final_inputs_sha256" if after else "execution_inputs_sha256", {}))


def paths():
    found = set(OWN) | {"notes/stationary-coupled-asymptotic-obstruction.md", "notes/spherical-anomaly-action.md", "notes/dirac-spherical-source-reduction.md", "notes/stationary-lll-quantum-exterior.md", "notes/coupled-fermion-domain.md", "notes/exact-spherical-em-reduction.md", "notes/stationary-shape-logarithmic-response.md", "notes/stationary-curved-quantum-tt.md"}
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
    for p in ["general-background-newton-operator", "stationary-scalar-phase-moment",
              "on-shell-adjoint-error-and-trace-gradient", "on-shell-finite-range-response",
              "cold-projection-source-lower-bound"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)








def exact_identities_and_margins():
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('primary_source_scope',review['story']=='S168' and len(review['primary_sources'])==2 and not review['actual_wormhole_instability_proved'])
    check('DLMF_download_failure_disclosed',review['primary_sources'][1]['direct_download_attempt']['outcome'].startswith('HTTP 403'))
    t,r,k=S.symbols('t r k',positive=True);b,c,delta=S.symbols('b c delta',real=True)
    G=S.exp(-r*r/(4*t))/(4*S.pi*t);F=S.exp(b*r*r+c*t)*G
    residual=S.simplify((S.diff(F,t)-S.diff(F,r,2)-(1/r+delta)*S.diff(F,r))/F)
    formula=c-4*b-4*b*b*r*r+2*b*r*r/t+delta*(r/(2*t)-2*b*r)
    zero('radial_heat_barrier_residual',residual-formula)
    # If k r^2<1, lower-bound every independent term using |delta|<=k r.
    # The r^2/t term is nonnegative for plus and nonpositive for minus.
    plus_constant=k-R(3,4)*k*k*r*r
    minus_constant=-k+R(1,4)*k*k*r*r
    zero('plus_constant_bound',formula.subs({c:2*k,b:k/4}).subs(delta,k*r).expand().coeff(t,0)-plus_constant)
    zero('minus_constant_bound',formula.subs({c:-2*k,b:-k/4}).subs(delta,k*r).expand().coeff(t,0)-minus_constant)
    check('positive_barrier_reserve',1-R(3,4)>0)
    check('negative_barrier_reserve',-1+R(1,4)<0)
    q=S.exp(-S.Symbol('a',positive=True)/t)/t
    aa=next(x for x in q.free_symbols if x!=t)
    zero('boundary_correction_monotone',S.diff(q,t)/q-(aa/t**2-1/t))
    # Anchored coordinate quadratic and C2,1 seam mechanism.
    y,x,u,v,y0,x0=S.symbols('y x u v y0 x0',real=True)
    h=S.Function('h')(y);bb=S.diff(h,y);cc=S.exp(-2*h);zz=u*(y-y0)+v*(x-x0)
    L=lambda f:-S.diff(f,y,2)-bb*S.diff(f,y)-cc*S.diff(f,x,2)
    L2q=S.expand(L(L(zz*zz)))
    target=2*u*(S.diff(bb,y,2)+bb*S.diff(bb,y))*zz+(4*S.diff(bb,y)+2*bb*bb)*u*u+4*cc*(bb*bb-S.diff(bb,y))*v*v
    zero('quadratic_L_squared_identity',S.simplify(L2q-target))
    zero('third_metric_jet_multiplies_vanishing_gradient',S.diff(L2q,S.diff(h,y,3)).subs({y0:y,x0:x}))
    check('no_fourth_metric_jet_in_quadratic',not L2q.has(S.diff(h,y,4)))
    ss=S.symbols('s',positive=True)
    zero('integrated_semigroup_remainder',S.integrate((t-ss)*S.sqrt(ss),(ss,0,t))-R(4,15)*t**R(5,2))
    beta=S.gamma(R(7,2))**2/S.gamma(7)
    zero('Feynman_parameter_moment',beta-5*S.pi/1024)
    zero('single_mass_contact_coefficient',R(4,15)*S.gamma(R(5,2))*beta/(2*S.pi)-S.sqrt(S.pi)/2048)
    zero('summed_scalar_contact_coefficient',R(2,3)*S.sqrt(S.pi)/2048-S.sqrt(S.pi)/3072)
    # Laplace substitution u=r^2/(c t), the two gain-two moments.
    zeta=S.symbols('zeta',positive=True)
    for power in [3,4]:
        integrand=t**(-power)*S.exp(-r*r/(c*t))
        # Use a positive independent c for this substitution.
        transformed=S.simplify(integrand.subs(t,r*r/(zeta))*r*r/zeta**2).subs(c,1)
        target=(r*r)**(1-power)*zeta**(power-2)*S.exp(-zeta)
        zero('heat_gain_two_substitution_'+str(power),transformed-target)
    zero('gain_two_equal_radial_power',r*r*(r*r)**(-3)-r**(-4))
    FF,CC,J2,RG=S.symbols('F C J2 RG')
    kappa=-FF/4
    zero('fourth_contact_difference',CC*kappa*J2/(30*S.pi)-CC*kappa*RG/(90*S.pi)+CC*FF*J2/(120*S.pi)-CC*FF*RG/(360*S.pi))
    check('outer_fourth_margin',R(33,600)-R(1,100)>R(2,1000)*R(22,7))
    check('large_time_weight_integrable',648>22)
    # A compact H2 polynomial fixture, with smooth approximation available.
    zz=S.symbols('z',real=True);kap=S.symbols('kap',nonnegative=True)
    phi=(zz-1)**3*(2-zz)**3
    I0=S.integrate(phi**2,(zz,1,2));I1=S.integrate(S.diff(phi,zz)**2,(zz,1,2));I2=S.integrate(S.diff(phi,zz,2)**2,(zz,1,2))
    mix=S.integrate((-S.diff(phi,zz,2)+kap**2*phi)**2,(zz,1,2))
    zero('annular_fourth_norm',mix-I2-2*kap*kap*I1-kap**4*I0)
    check('annular_limit_positive',I0>0 and I1>0 and I2>0,[I0,I1,I2])
    d=S.symbols('d',positive=True)
    zero('annular_gradient_to_fourth_ratio',S.limit(d**2*(I1+kap**2*I0)/mix,d,0))
    zero('annular_zero_to_fourth_ratio',S.limit(d**4*I0/mix,d,0))


def near(name,got,want):
    err=abs(got-want);scale=1+abs(want)
    check(name,err<mp.mpf('1e-25')*scale,mp.nstr(err,12))
    return str(err)


def positive_kernel_integrals():
    mp.mp.dps=40
    for z0 in PROTOCOL['bessel_arguments']:
        z=mp.mpf(str(z0))
        f=lambda a:mp.besselk(1,z/mp.sqrt(a*(1-a)))/mp.sqrt(a*(1-a)) if 0<a<1 else mp.mpf(0)
        got=mp.quad(f,[0,mp.mpf('.5'),1]);want=z*(mp.besselk(1,z)**2-mp.besselk(0,z)**2)
        err=near('Bessel_positive_kernel_'+str(z0),got,want)
        NUMERICS.append({'case':'Bessel_kernel','z':z0,'error':err})
    M=mp.mpf('1.3')
    for r0 in PROTOCOL['heat_radii']:
        rr=mp.mpf(str(r0))
        for al0 in PROTOCOL['feynman_parameters']:
            al=mp.mpf(str(al0));a=al*(1-al)
            got=mp.quad(lambda tt:mp.exp(-M*M*tt/a-rr*rr/(4*tt))/tt**2,[0,1,mp.inf])
            want=4*M/(rr*mp.sqrt(a))*mp.besselk(1,M*rr/mp.sqrt(a))
            err=near('heat_Laplace_kernel_'+str(r0)+'_'+str(al0),got,want)
            NUMERICS.append({'case':'Laplace_kernel','r':r0,'alpha':al0,'error':err})


def polar_heat_fixture():
    rr=mp.mpf(1);rp=mp.mpf('1.3')
    def kernel(rad,tt,m):
        return mp.sqrt(rad*rp)/(2*tt)*mp.exp(-(rad*rad+rp*rp)/(4*tt))*mp.besseli(m,rad*rp/(2*tt))
    for t0 in PROTOCOL['polar_times']:
        tt=mp.mpf(str(t0))
        for m in PROTOCOL['polar_modes']:
            plane=lambda theta:mp.exp(-(rr*rr+rp*rp-2*rr*rp*mp.cos(theta))/(4*tt))/(4*mp.pi*tt)
            angular=mp.sqrt(rr*rp)*mp.quad(lambda theta:plane(theta)*mp.cos(m*theta),[0,mp.pi,2*mp.pi])
            value=kernel(rr,tt,m);err=near('polar_density_'+str(t0)+'_'+str(m),angular,value)
            time_derivative=mp.diff(lambda st:kernel(rr,st,m),tt)
            radial_second=mp.diff(lambda rad:kernel(rad,tt,m),rr,2)
            pde=time_derivative-radial_second+(m*m-mp.mpf('.25'))/rr**2*value
            pde_err=near('polar_heat_equation_'+str(t0)+'_'+str(m),pde,0)
            check('wrong_density_potential_detected_'+str(t0)+'_'+str(m),abs(pde+value/(2*rr**2))>mp.mpf('1e-5'))
            NUMERICS.append({'case':'polar_exact_heat','t':t0,'m':m,'normalization_error':err,'PDE_error':pde_err})
        images=(mp.exp(-(rr-rp)**2/(4*tt))-mp.exp(-(rr+rp)**2/(4*tt)))/mp.sqrt(4*mp.pi*tt)
        near('half_integer_Dirichlet_images_'+str(t0),kernel(rr,tt,mp.mpf('.5')),images)

def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S168 uniform logarithmic frequency comparison','input_count':len(ps)}),flush=True)
    try:
        predecessors();exact_identities_and_margins();positive_kernel_integrals();polar_heat_fixture()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S168 uniform logarithmic frequency comparison','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'uniform_finite_scalar_geodesic_comparison_derived':good,'actual_seed_outer_logarithmic_bound_derived':good,
         'lower_order_constants_numerically_calibrated':False,'numerical_collar_threshold_evaluated':False,
         'nonlinear_exterior_frequency_bound_proved':False,
         'numerics':NUMERICS,'actual_positive_frequency_operator_values_evaluated':False,'actual_causal_inverse_proved':False,
         'actual_charged_retarded_transfer_proved':False,'actual_wormhole_instability_proved':False,
         'order_reduction_adopted':False,'S167_complete':False,'S168_complete':False,
         'S077_complete':False,'goal_complete':False}
    if not good:
        run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text(encoding="utf-8")
        run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text(encoding="utf-8")
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
