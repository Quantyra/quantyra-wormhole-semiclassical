"""S168: thermal shape response and frequency-subtracted spectral bounds.

Exact metric algebra and rational budgets; continuum claims remain ordinary proofs.
"""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import sys
import traceback
import sympy as S
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
STEM = "thermal-shape-response-and-subtractions"
OWN = ["checks/thermal_shape_response_and_subtractions.py", "notes/"+STEM+".md", "evidence/thermal-shape-response-source-review.json"]
PREDECESSORS = [("charged-conformal-response-and-constraints", 323)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "Full optical mass vertices; thermal spectral response and three-frequency bound",
    "digits":70,"beta":"13/10","energies":["-2","1/2","3"],
    "imaginary_frequencies":[".4","1.1","2.7"],"commutator_time_cut":256,
    "laplace_scaled_tolerance":"1e-30","matsubara":[0,1,2,3],
    "euclidean_scaled_tolerance":"1e-50","free_energy_steps":["1e-4","5e-5"],
    "difference_scaled_tolerance":"1e-6","refinement_ratio":".35",
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
    found = set(OWN) | {"notes/stationary-coupled-asymptotic-obstruction.md", "notes/spherical-anomaly-action.md", "notes/dirac-spherical-source-reduction.md", "notes/stationary-lll-quantum-exterior.md", "notes/coupled-fermion-domain.md", "notes/stationary-shape-logarithmic-response.md", "notes/stationary-curved-quantum-tt.md"}
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





def exact_identities():
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('source_review_scope',review['story']=='S168' and len(review['primary_sources'])==1 and not review['actual_wormhole_instability_proved'])
    lam,V,B,v,w,a,b=S.symbols('lambda V B v w a b',nonzero=True)
    mass=lam*S.exp(V+a*v+b*w)/B
    zero('optical_first_vertex',S.diff(mass,a).subs({a:0,b:0})-lam*S.exp(V)*v/B)
    zero('optical_second_contact',S.diff(mass,a,b).subs({a:0,b:0})-lam*S.exp(V)*v*w/B)
    zero('LLL_mass_vertex_zero',S.diff(mass,a).subs(lam,0))
    sigma2=S.Matrix([[0,-S.I],[S.I,0]]);sigma1=S.Matrix([[0,1],[1,0]])
    zero('real_radial_derivative_coefficient',S.im(-S.I*sigma2[0,1]))
    check('real_mass_coefficient',all(S.im(x)==0 for x in sigma1))
    X,Y=S.symbols('X Y',positive=True);fa=1/(1+X);fb=1/(1+X*Y)
    zero('KMS_transition_identity',fa-fb-(1-1/Y)*fa*(1-fb))
    delta,om,gamma,s,t,s0=S.symbols('Delta omega gamma s t s0',positive=True)
    zero('retarded_pair_denominator',1/(om-delta)-1/(om+delta)-2*delta/(om**2-delta**2))
    zero('imaginary_frequency_bubble_sign',(2*delta/(om**2-delta**2)).subs(om,S.I*gamma)+2*delta/(delta**2+gamma**2))
    f=-1/(t+s);tay=sum(S.diff(f,s,k).subs(s,s0)*(s-s0)**k/S.factorial(k) for k in range(3))
    zero('third_resolvent_subtraction',f-tay-(s-s0)**3/((t+s)*(t+s0)**3))
    zero('positive_third_derivative',S.diff(f,s,3)-6/(t+s)**4)
    nodes=[S.Rational(1,3),S.Rational(4,3),S.Rational(7,3)]
    interp=sum(f.subs(s,u)*S.prod((s-vv)/(u-vv) for vv in nodes if vv!=u) for u in nodes)
    zero('quadratic_interpolation_remainder',f-interp-S.prod(s-u for u in nodes)/((t+s)*S.prod(t+u for u in nodes)))
    four=[S.Rational(1,3),S.Rational(4,3),S.Rational(7,3),S.Rational(11,3)]
    dd=sum(f.subs(s,u)/S.prod(u-vv for vv in four if vv!=u) for u in four)
    zero('four_frequency_divided_difference',dd-1/S.prod(t+u for u in four))
    c0,c1,c2=S.symbols('c0 c1 c2');poly=c0+c1*s+c2*s**2
    zero('all_local_contacts_cancel',sum(poly.subs(s,u)/S.prod(u-vv for vv in four if vv!=u) for u in four))
    plus=1+s+s*s+R(1,2)-1/(2+s);minus=1-10*s+s*s+R(1,2)-1/(2+s)
    zero('counterexample_equal_static',plus.subs(s,0)-minus.subs(s,0))
    zero('counterexample_positive_at_one',plus.subs(s,1)-R(19,6))
    zero('counterexample_negative_at_one',minus.subs(s,1)+R(47,6))
    zero('counterexample_same_third_derivative',S.diff(plus-minus,s,3))
    zero('counterexample_same_quadratic_contact',S.diff(plus-minus,s,2))
    check('counterexample_gap_and_passivity_insufficient',minus.subs(s,0)>0 and minus.subs(s,1)<0 and plus.subs(s,1)>0)
    beta,e=S.symbols('beta e',positive=True);occ=1/(1+S.exp(beta*e))
    zero('degenerate_thermal_susceptibility',S.diff(occ,e)+beta*occ*(1-occ))
    # UV cumulative growth T^3 log^p T with denominator degree four is integrable.
    zero('UV_power_count',8-3-5)
    zero('positive_measure_cumulative_power',1+5-2*3)
    check('four_denominator_UV_integrable',3-5 < -1)


def finite_thermal_fixture():
    mp.mp.dps=PROTOCOL['digits'];beta=mp.mpf(13)/10
    energies=[mp.mpf(-2),mp.mpf(1)/2,mp.mpf(3)]
    h=mp.diag(energies)
    B=mp.matrix([[mp.mpf(7)/10,mp.mpf(3)/10,mp.mpf(1)/2],[mp.mpf(3)/10,-mp.mpf(1)/5,-mp.mpf(1)/5],[mp.mpf(1)/2,-mp.mpf(1)/5,mp.mpf(2)/5]])
    W=mp.matrix([[mp.mpf(2)/3,mp.mpf(1)/7,0],[mp.mpf(1)/7,mp.mpf(3)/5,mp.mpf(1)/9],[0,mp.mpf(1)/9,-mp.mpf(1)/4]])
    occ=[1/(1+mp.exp(beta*e)) for e in energies];F=mp.diag(occ)
    tr=lambda m:sum(m[i,i] for i in range(3))
    contact=tr(F*W);Q=beta*sum(occ[i]*(1-occ[i])*B[i,i]**2 for i in range(3))
    pairs=[(energies[j]-energies[i],(occ[i]-occ[j])*abs(B[i,j])**2) for i in range(3) for j in range(i+1,3)]
    chi=lambda s:sum(2*d*w/(d*d+s) for d,w in pairs)
    D=lambda s:contact-chi(s)
    check('positive_KMS_pair_weights',all(d>0 and w>0 for d,w in pairs))
    check('finite_degenerate_block_nonzero',Q>0)
    def commutator(time):
        U=mp.diag([mp.exp(1j*e*time) for e in energies]);Ui=mp.diag([mp.exp(-1j*e*time) for e in energies])
        Bt=U*B*Ui
        return -1j*tr(F*(Bt*B-B*Bt))
    for gs in PROTOCOL['imaginary_frequencies']:
        g=mp.mpf(gs)
        integral=mp.quad(lambda time:mp.exp(-g*time)*commutator(time),[0,1,2,4,8,16,32,64,128,256])
        tail=2*sum(w for d,w in pairs)*mp.exp(-g*256)/g
        error=abs(contact+integral-D(g*g));scale=1+abs(D(g*g))
        check('causal_commutator_Laplace_'+gs,error+tail < mp.mpf('1e-30')*scale,mp.nstr(error+tail,15))
        check('wrong_bubble_sign_detected_'+gs,abs(contact+integral-(contact+chi(g*g)))>mp.mpf('1e-4'))
        NUMERICS.append({'case':'Laplace','gamma':gs,'error':mp.nstr(error,25),'tail_bound':mp.nstr(tail,25)})
    def euclidean(tau):
        U=mp.diag([mp.exp(e*tau) for e in energies]);Ui=mp.diag([mp.exp(-e*tau) for e in energies])
        return tr(F*U*B*Ui*(mp.eye(3)-F)*B)
    for n in PROTOCOL['matsubara']:
        nu=2*mp.pi*n/beta
        got=mp.quad(lambda tau:mp.cos(nu*tau)*euclidean(tau),[0,beta/2,beta])
        target=chi(nu*nu)+(Q if n==0 else 0);error=abs(got-target)
        check('independent_Euclidean_covariance_'+str(n),error<mp.mpf('1e-50')*(1+abs(target)),mp.nstr(error,15))
        if n==0:check('static_population_omission_detected',abs(got-chi(0))>mp.mpf('.01'))
        NUMERICS.append({'case':'Matsubara','n':n,'error':mp.nstr(error,25)})
    def free_energy(a):
        eig=mp.eigsy(h+a*B+a*a*W/2,eigvals_only=True)
        return -sum(mp.log(1+mp.exp(-beta*en)) for en in eig)/beta
    target=contact-chi(0)-Q;errors=[]
    for step in PROTOCOL['free_energy_steps']:
        a=mp.mpf(step);got=(free_energy(a)+free_energy(-a)-2*free_energy(0))/(a*a)
        error=abs(got-target);errors.append(error)
        check('eigenvalue_equilibrium_Hessian_'+step,error<mp.mpf('1e-6')*(1+abs(target)),mp.nstr(error,15))
    check('equilibrium_refinement',errors[1]<mp.mpf('.35')*errors[0],[mp.nstr(v,15) for v in errors])
    check('seagull_omission_detected',abs(contact)>mp.mpf('.1'))
    check('causal_vs_equilibrium_population',abs(D(0)-target-Q)<mp.mpf('1e-60'))
    NUMERICS.append({'case':'equilibrium','target':mp.nstr(target,25),'errors':[mp.nstr(v,25) for v in errors]})
    # An actually degenerate block tests its off-diagonal population term.
    de=[mp.mpf(-2),mp.mpf('.5'),mp.mpf('.5')];dh=mp.diag(de)
    df=[1/(1+mp.exp(beta*e)) for e in de]
    dq=beta*sum(df[i]*(1-df[i])*abs(B[i,j])**2 for i in range(3) for j in range(3) if de[i]==de[j])
    dc=sum(2*(df[i]-df[j])*abs(B[i,j])**2/(de[j]-de[i]) for i in range(3) for j in range(i+1,3) if de[i]<de[j])
    dtarg=tr(mp.diag(df)*W)-dc-dq
    def degenerate_energy(a):
        eig=mp.eigsy(dh+a*B+a*a*W/2,eigvals_only=True)
        return -sum(mp.log(1+mp.exp(-beta*en)) for en in eig)/beta
    derr=[]
    for hs in PROTOCOL['free_energy_steps']:
        hstep=mp.mpf(hs)
        got=(degenerate_energy(hstep)+degenerate_energy(-hstep)-2*degenerate_energy(0))/(hstep*hstep)
        derr.append(abs(got-dtarg))
        check('degenerate_block_Hessian_'+hs,derr[-1]<mp.mpf('1e-6')*(1+abs(dtarg)),mp.nstr(derr[-1],15))
    check('degenerate_block_refinement',derr[1]<mp.mpf('.35')*derr[0])
    check('off_diagonal_degenerate_omission_detected',2*beta*df[1]*(1-df[1])*B[1,2]**2>mp.mpf('1e-3'))
    NUMERICS.append({'case':'degenerate_equilibrium','errors':[mp.nstr(v,25) for v in derr]})


def independent_spectral_quadrature():
    mp.mp.dps=70
    mu=lambda t:t*t*mp.exp(-t)
    D=lambda s:2-3*s+5*s*s-mp.quad(lambda t:mu(t)/(t+s),[0,1,mp.inf])
    anchors=[mp.mpf('.3'),mp.mpf('1.2'),mp.mpf('2.4')];s=mp.mpf('3.7')
    interp=sum(D(a)*mp.fprod((s-b)/(a-b) for b in anchors if b!=a) for a in anchors)
    answer=mp.fprod(s-a for a in anchors)*mp.quad(lambda t:mu(t)/((t+s)*mp.fprod(t+a for a in anchors)),[0,1,mp.inf])
    error=abs(D(s)-interp-answer)
    check('continuous_measure_interpolation',error<mp.mpf('1e-60'),mp.nstr(error,20))
    check('continuous_measure_lower_bound',D(s)>interp and answer>0)
    NUMERICS.append({'case':'continuous_positive_measure','error':mp.nstr(error,25),'positive_remainder':mp.nstr(answer,25)})
    # A genuinely UV-divergent raw Stieltjes integral; third subtraction converges.
    a=mp.mpf(1);s=mp.mpf(3)
    kernel=lambda t:(s-a)**3*t*t/((t+s)*(t+a)**3)
    exact=mp.quad(kernel,[1,mp.inf]);previous=mp.mpf(0)
    for cap in [100,1000]:
        val=mp.quad(kernel,[1,10,cap]);tail=exact-val
        check('UV_subtracted_tail_'+str(cap),0<tail<(s-a)**3/cap,mp.nstr(tail,20))
        check('UV_subtracted_monotone_'+str(cap),previous<val<exact)
        previous=val
        raw_difference=s*((cap-1)-s*mp.log((cap+s)/(1+s)))
        independent=mp.quad(lambda t:s*t/(t+s),[1,10,cap])
        check('UV_raw_difference_integral_'+str(cap),abs(raw_difference-independent)<mp.mpf('1e-60'))
        check('UV_raw_difference_linear_lower_bound_'+str(cap),raw_difference>s*(cap-s)/2)

def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S168 thermal shape response and frequency subtractions','input_count':len(ps)}),flush=True)
    try:
        predecessors();exact_identities();finite_thermal_fixture();independent_spectral_quadrature()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S168 thermal shape response and frequency subtractions','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'thermal_shape_spectral_row_derived':good,'third_frequency_form_sign_derived':good,
         'three_positive_frequency_lower_bound_derived':good,'static_passivity_shortcut_refuted':good,
         'numerics':NUMERICS,'actual_positive_frequency_operator_values_evaluated':False,'actual_causal_inverse_proved':False,
         'actual_charged_retarded_transfer_proved':False,'actual_wormhole_instability_proved':False,
         'order_reduction_adopted':False,'S167_complete':False,'S168_complete':False,
         'S077_complete':False,'goal_complete':False}
    if not good:
        run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text()
        run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text()
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
