"""S167: on-shell full finite-comparison value and first-derivative response.

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
STEM = "on-shell-finite-range-response"
OWN = ["checks/on_shell_finite_range_response.py", "notes/"+STEM+".md"]
PREDECESSORS = [("finite-trace-gradient-and-projection-transversality", 357),
                ("on-shell-adjoint-error-and-trace-gradient", 329),
                ("projected-low-correction-and-mass", 294),
                ("projected-first-derivative-bound", 246)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "actual finite Hessian full H1 response; cold/higher neighborhood and nonlinear solve open",
    "source_norm": "F=||t^2 f_tau||_ell+||f_c||_ell; f_tau,f_c scaled by lambda^-2",
    "range_bounds": "||metric||<800 g F; ||gradbar metric||<16000 g F",
    "comparison_bounds": "value/sqrt(gamma0)<1.1e-7 ell^2; gradient/sqrt(gamma0)<2.2e-6 ell^2",
    "operator_split": "V_T=(Tbar-E)C+Q_T; local C term retained",
    "full_return_and_mass_tail_retained": True,
    "finite_mass_coefficient_bound": False,
    "fixture_scope": "exact noncommuting full-block algebra, compatible terminal integration, rational bounds",
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
    for p in ["finite-trace-gradient-and-projection-transversality", "on-shell-adjoint-error-and-trace-gradient",
              "finite-comparison-weighted-mixed-form", "finite-comparison-tt-realization-and-range",
              "finite-comparison-adjoint-and-fixed-projection", "finite-comparison-disk-return-form",
              "global-compatible-scalar-response", "projected-low-correction-and-mass",
              "projected-first-derivative-bound", "endpoint-preserving-compatible-projection"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)


















def positive_elimination_budgets():
    sg=R('4e-17');b0=R('2.6e-13');sb=R('5.1e-7');infl=R('1.001')
    check('local_root_reserve',sb**2>b0)
    check('inverse_local_reserve',1/(1-b0)<infl)
    check('D_square_root',8*R(22,7)/R('.90')<R('5.3')**2)
    check('full_TT_datum_constant',640*R(22,7)**3/R('.90')<150**2)
    check('scalar_datum_constant',80*R(22,7)**2<800)
    N=S.Matrix([[R('1.000001'),R('2.6e-7')],[R('2.6e-7'),R('1.000001')]])
    L=S.Matrix([[R('.02405'),R('2.405e-8')],[R('2.405e-8'),R('2.405e-14')]])
    K=N*L;H=(S.eye(2)-K).inv()
    bound=S.Matrix([[R('1.026'),R('3e-8')],[R('4e-8'),R('1.000001')]])
    check('block_series_converges',max(sum(K[i,j] for j in range(2)) for i in range(2))<R('.025'))
    for i in range(2):
        for j in range(2):check('complex_inverse_block_'+str(i)+str(j),0<=H[i,j]<bound[i,j],H[i,j])
    check('weight_product',32*200==80**2)
    outgoing=S.Matrix([80*(R('1.001')*bb/13+3*aa)*800*sg/180 for aa,bb in [(4*10**6,1000),(4,R('2e-5'))]])
    elliptic=S.Matrix([infl*800/(2*13),infl*sb*800/(2*13)])
    from_scalar=N*outgoing+elliptic
    from_TT=N*S.Matrix([150,0])
    source=S.Matrix([R('150.001'),R('3.91e-5')])
    for i in range(2):
        check('preconditioned_scalar_datum_'+str(i),from_scalar[i]<source[i],from_scalar[i])
        check('preconditioned_TT_datum_'+str(i),from_TT[i]<source[i],from_TT[i])
    final=bound*source
    check('TT_auxiliary_response',final[0]<155,final[0])
    check('local_auxiliary_response',final[1]<R('5e-5'),final[1])


def on_shell_and_terminal_budgets():
    sg=R('4e-17');g=sg**2
    check('complete_trace_pointwise_gradient',R('4e-5')+g/180*(10**28+10**17)+R('1e-24')<R('4.1e-5'))
    check('local_split_value',R('1.2')/180*8*25<R('1.34'))
    check('local_split_gradient',R('1.2')/180*(10**15*25+8*420)<R('1.7e14'))
    check('remainder_divergence_coefficient',R('1.2')/180*200*420<561)
    value=R('1.2')/180*(g*(R('6e15')*420+R('2e20')*25)+R('20.02')*25*R('4e-27'))+R('1.2')*(160*R('3.2e-12')+25*R('4.1e-5'))
    check('remainder_value_coefficient',value<R('.0013'),value)
    moment=R(178**2,2)+R('1.5')/40*(178+R('1.5'))
    check('terminal_weighted_second_moment',moment<16000<127**2,moment)
    check('terminal_scalar_coefficient',R('2.001')*127<300)
    check('terminal_divergence_coefficient',4*127<1000)
    A=800+R('.0013')*155+R('2e-5')*R('5e-5')
    B=561*155*sg+4*R('5e-5')
    check('modified_source_scalar',A<801,A)
    check('modified_source_divergence',B<R('.00021'),B)
    check('modified_outgoing_value',R(300,13)*801+1000*R('.00021')<19000)
    check('outgoing_energy_cross_term',801*19000*sg<1)
    check('modified_outgoing_gradient',15*19000+1+R('.00021')*sg<286000)
    check('restored_outgoing_value',19000+R('1.34')*155*sg<20000)
    check('restored_outgoing_gradient',286000+R('1.7e14')*155*g<290000)
    # Exact operator coefficient split, retaining the complete leading trace term.
    eps,gg,E,T,r0,local,ret=S.symbols('eps g E T r local ret')
    vertex=-S.sqrt(5)/(2*gg)*(-eps*(T-E)*r0+eps*local-ret)
    zero('exact_TT_operator_split',vertex-(S.sqrt(5)*eps/(2*gg)*(T-E)*r0-S.sqrt(5)*eps/(2*gg)*local+S.sqrt(5)/(2*gg)*ret))


def scalar_and_full_response_budgets():
    sg=R('4e-17');infl=R('1.001');sb=R('5.1e-7')
    check('TT_vertex_on_outgoing_field',1000*20000+4*10**6*290000<R('1.2e12'))
    check('local_vertex_on_outgoing_field',R('2e-5')*20000+4*290000<R('1.2e6'))
    check('energy_denominator',180**3>2400**2)
    vertex=infl*(R('1.2e12')+sb*R('1.2e6'))*sg/2400
    value_remainder=20000*sg/180+vertex/2
    check('scalar_value_remainder',value_remainder<R('1e-7'),value_remainder)
    val=max(infl*R(800,4*180),infl*R(150,2*13))+value_remainder
    check('scalar_value_response',val<6,val)
    deriv=R(290000,180)+max(infl*R(800,2*180),infl*R(150,13))+vertex
    check('scalar_gradient_response',deriv<1700,deriv)
    check('scalar_TT_value_response',(R('5.3')*155)**2/(80*3**2)<31**2)
    check('finite_disk_value_response',6+25*31<800)
    check('finite_return_derivative_coefficient',R('.8')*420+73*25==2161)
    check('finite_TT_return_gradient',(2161*155)**2/(R('.90')*80*3**2)<14000**2)
    check('full_disk_gradient_response',1700+14000<16000)
    check('unweighted_TT_source_conversion',26*R('1.000001')<27)
    check('fixed_projection_composition',1+R('1.006')*R('21.3')<R('22.5'))
    forcing=27*R('22.5')*R('2.2e-13')
    check('actual_finite_metric_value',800*forcing<R('1.1e-7'),800*forcing)
    check('actual_finite_metric_gradient',16000*forcing<R('2.2e-6'),16000*forcing)
    check('actual_finite_scalar_value',6*forcing<R('8.1e-10'))
    check('actual_finite_scalar_gradient',1700*forcing<R('2.3e-7'))
    check('actual_finite_TT_value',31*forcing<R('4.2e-9'))


def noncommuting_full_block_fixture():
    T=S.Matrix([[5,1],[1,8]]);H0=T+4*S.eye(2);E=40;e=E+4;scale=2
    D=S.diag(4,9);Droot=S.diag(2,3);Di=Droot.inv()
    VT=S.Matrix([[2,1],[1,3]]);Vs=S.Matrix([R('.01'),R('-.005')]);V=VT.row_join(Vs);J=S.diag(1,1,-1)
    B=VT*Droot/scale;Lscalar=(H0*(T-E*S.eye(2))+Vs*Vs.T)/scale**2
    full=Lscalar.row_join(B).col_join(B.T.row_join(D))
    ft=S.Matrix([R(2,3),R(-5,7)]);fc=S.Matrix([R(11,13),R(3,8)]);source=ft.col_join(fc)
    solution=full.inv()*source;tau=solution[:2,0];c=solution[2:,0]
    p0=scale**2*ft;r=scale*Di*fc;w=V.T*tau-r.col_join(S.zeros(1,1))
    Rp=(T-E*S.eye(2)).inv();R4=(Rp-H0.inv())/e
    A=V.T*H0.inv()*V/e;L=V.T*Rp*V/e;N=(S.eye(3)+A*J).inv();M=S.eye(3)+A*J-L*J
    check('fixture_original_full_equation',full*solution==source)
    check('fixture_auxiliary_equation',M*w==V.T*R4*p0-r.col_join(S.zeros(1,1)))
    check('fixture_positive_elimination',w==(S.eye(3)-N*L*J).inv()*N*(V.T*R4*p0-r.col_join(S.zeros(1,1))))
    check('fixture_TT_reconstruction',c==-Di*w[:2,0]/scale)
    q=p0+V*J*w;z=Rp*q;K=H0+V*J*V.T/e;p=p0-VT*r
    check('fixture_noncommuting_elliptic_operator',K*H0!=H0*K)
    check('fixture_partial_fraction_identity',tau==R4*q and z==H0*tau)
    check('fixture_exact_elliptic_regrouping',tau==z/e-K.inv()*p/e-K.inv()*V*J*V.T*z/e**2)
    C=S.Matrix([[R('.1'),R('.02')],[R('-.03'),R('.05')]]);Q=VT-(T-E*S.eye(2))*C
    z1=z-C*w[:2,0];q1=p0+Q*w[:2,0]-Vs*w[2,0]
    check('fixture_on_shell_source_split',(T-E*S.eye(2))*z1==q1)
    check('fixture_local_split_restoration',Rp*q1+C*w[:2,0]==z)
    check('omitting_local_split_term_detected',Rp*q1!=z)
    check('omitting_elliptic_subtraction_detected',Rp*q/e!=tau)
    bad=(S.eye(3)-L*J).inv()*(V.T*R4*p0-r.col_join(S.zeros(1,1)))
    check('omitting_positive_elimination_detected',M*bad!=M*w)


def independent_terminal_fixture():
    import mpmath as mp
    mp.mp.dps=80
    x=S.symbols('x',positive=True);target=x*S.exp(-x);b=x*x*S.exp(-x)
    q=-S.diff(target,x,2)-9*target;a=q-S.diff(b,x)
    zero('terminal_fixture_exact_compatibility',S.integrate(S.sin(3*x)*q,(x,0,S.oo)))
    r,s=S.symbols('r s',real=True);kernel=S.sin(3*(r-s))/3
    zero('terminal_fixture_diagonal_value',kernel.subs(s,r))
    zero('terminal_fixture_diagonal_derivative',S.diff(kernel,r).subs(s,r)-1)
    aa=S.lambdify(x,a,'mpmath');bb=S.lambdify(x,b,'mpmath')
    for j,rr in enumerate([mp.mpf(0),mp.mpf('.2'),mp.mpf(1),mp.mpf(3)]):
        intervals=[rr,rr+5,rr+20,rr+60,rr+150,mp.inf]
        value=mp.quad(lambda ss:mp.sin(3*(rr-ss))/3*aa(ss)+mp.cos(3*(rr-ss))*bb(ss),intervals)
        expected=rr*mp.exp(-rr)
        check('independent_terminal_integral_'+str(j),abs(value-expected)<mp.mpf('1e-55'),mp.nstr(abs(value-expected),12))
        omitted=mp.quad(lambda ss:mp.sin(3*(rr-ss))/3*aa(ss),intervals)
        check('terminal_divergence_omission_'+str(j),abs(omitted-expected)>mp.mpf('1e-5'))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 on-shell finite range response','input_count':len(ps)}),flush=True)
    try:
        predecessors();positive_elimination_budgets();on_shell_and_terminal_budgets();scalar_and_full_response_budgets();noncommuting_full_block_fixture();independent_terminal_fixture()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 on-shell finite range response','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'finite_full_low_order_response_bounded':good,'finite_projected_value_gradient_bounded':good,
         'on_shell_operator_split_retained':good,'full_return_and_mass_tail_retained':good,
         'finite_mass_coefficient_bound':False,'useful_full_higher_log_inverse_constant':False,
         'cold_pointwise_control':False,'nonlinear_family_preserved':False,'projected_nonlinear_solve':False,
         'signed_physical_compatibility_evaluated':False,'S167_complete':False,'goal_complete':False}
    if not good:
        run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text()
        run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text()
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
