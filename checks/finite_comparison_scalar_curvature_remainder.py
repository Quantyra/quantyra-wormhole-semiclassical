"""S167: scalar curvature cancellation and local Schur absorption.

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
STEM = "finite-comparison-scalar-curvature-remainder"
OWN = ["checks/finite_comparison_scalar_curvature_remainder.py", "notes/"+STEM+".md"]
PREDECESSORS = [("finite-comparison-physical-tt-contacts", 304),
                ("comparison-path-complete-equation-residual", 329),
                ("stationary-scalar-outgoing-inverse", 107)]
R = S.Rational
ROWS = []
PROTOCOL = {
    "scope": "actual scalar curvature remainder and local scalar low-frequency block",
    "curvature_bound": "abs(R4+12lambda)<3.1e-12 lambda/t",
    "local_auxiliary_norm_bound": "2.6e-13",
    "complete_trace_defect_bound": False,
    "mixed_coefficient_estimate_complete": False,
    "full_candidate_operator_realized": False,
    "fixture_scope": "exact algebra and rational budgets; independent warped metric identities",
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
    for p in ["finite-comparison-physical-tt-contacts", "finite-comparison-classical-disk-contact",
              "finite-comparison-classical-tt-form", "stationary-centered-constraint-residual",
              "comparison-path-complete-equation-residual", "stationary-scalar-outgoing-inverse",
              "general-background-newton-operator", "stationary-weighted-tt-inverse"]:
        check("required_pinned_"+p, "notes/"+p+".md" in expected)




def scalar_algebra():
    A,Rr,K,E2,r,ell,u,Sv,up,upp,sp,spp=S.symbols('A R K E2 r ell u S up upp sp spp')
    w=ell*u+Sv
    Ke=S.exp(-2*w)*(K+ell*(upp+A*up)+spp+A*sp)
    E2e=S.exp(-2*w)*(E2+spp+A*sp)
    Are=S.exp(-2*w)*(A+ell*up+sp)*(Rr+sp)
    R2e=S.exp(-2*w)*(Rr+sp)**2
    kr=-Ke; ks=-(E2e-Are+R2e); kt=-Are; kss=S.exp(-2*Sv)/r**2-R2e
    actual=2*(kr+2*ks+2*kt+kss)
    base=-2*K-4*E2-6*Rr**2+2/r**2
    formula=S.exp(-2*w)*(base-2*ell*(upp+A*up)-6*(spp+(A+2*Rr)*sp)-6*sp**2+2*(S.exp(2*ell*u)-1)/r**2)
    check('independent_sectional_scalar_identity',S.simplify(actual-formula)==0)
    x=S.symbols('x'); T=S.Function('T')(x); V=S.Function('V')(x); tau=S.Function('tau')(x)
    drift=S.symbols('drift'); k=S.exp(ell*V); TT=1+tau*(k-1); nu=tau*V; B=k-1-ell*V*k
    L=lambda f:S.diff(f,x,2)+drift*S.diff(f,x)
    zero('exponential_radius_first',S.diff(TT,x)-ell*k*S.diff(nu,x)-B*S.diff(tau,x))
    zero('exponential_radius_second',L(TT)-ell*k*L(nu)-B*L(tau)-ell**2*k*tau*S.diff(V,x)**2)
    zero('logarithm_gradient_square',L(S.log(T))+(S.diff(T,x)/T)**2-L(T)/T)
    tt=S.symbols('t',positive=True)
    zero('return_ratio_cancellation',(k/(1+(1-1/tt)*(k-1))-1)-(k-1)/(tt*(1+(1-1/tt)*(k-1))))
    zero('L4_tau',((2*A*Rr-2*Rr**2)+2*Rr**2)/tt-2*A*Rr/tt)
    q,qq,dd,ln,lt,vy2,n,Bs,ks,Ts=S.symbols('q qq dotR Lnu Ltau Vy2 n B k T')
    e=S.symbols('e'); eu=S.symbols('eu')
    # Use the exact linear trace identity to remove D2u+3L4nu.
    l4s_plus=ell*ks*ln/Ts+Bs*lt/Ts+ell**2*ks*vy2/Ts
    exact=e*(q+ell*dd+2*ell*n*q+6*ell*ln-4*ell*u/r**2-6*l4s_plus+2*(eu-1)/r**2)-q-ell*dd
    rem=q*(e*(1+2*ell*n)-1)+ell*(e-1)*dd-6*e*(ell*(ks/Ts-1)*ln+Bs*lt/Ts+ell**2*ks*vy2/Ts)+2*e*(eu-1-2*ell*u)/r**2
    zero('exact_curvature_remainder',exact-rem)
    # Pointwise omission controls; these values are algebra fixtures only.
    check('missing_gradient_square_detected',(-6*e*ell**2*ks*vy2/Ts).subs({e:1,ell:1,ks:2,vy2:3,Ts:4})!=0)
    check('missing_radius_acceleration_detected',(-6*e*Bs*lt/Ts).subs({e:1,Bs:2,lt:3,Ts:4})!=0)
    s,v,ta=S.symbols('s v tau',real=True)
    logt=S.log(1-ta+ta*S.exp(s*v))
    zero('radius_log_parameter_second',S.diff(logt,s,2)-ta*(1-ta)*v**2*S.exp(s*v)/(1-ta+ta*S.exp(s*v))**2)
    f=S.exp(-2*s)*(1+2*s)-1
    zero('value_remainder_double_zero',f.subs(s,0)+S.diff(f,s).subs(s,0))
    zero('value_remainder_second',S.diff(f,s,2)-(8*s-4)*S.exp(-2*s))
    poly=(tt-1)*(1+2/tt+3/tt**2)
    zero('gradient_coefficient_monotone',S.diff(poly,tt)-(1-1/tt**2+6/tt**3))
    zero('gradient_coefficient_at_four',poly.subs(tt,4)-R(81,16))


def budgets():
    tiny=R('1e-12')/R('7.9e9'); U=R('1.2e-8'); LN=R('2.2e-6')
    check('sharp_seed_curvature',10*tiny<R('1.3e-21'))
    check('sharp_linear_trace',R('6.1e-20')+U*56*tiny<R('6.2e-20'))
    check('radius_L4_bound',R('6.3e-7')+R('2.24')*R('7e-7')<LN)
    # Linear bounds for exponentials, proved by the geometric series.
    check('small_exponential_reserve',1/(1-R('6e-8'))<R('1.001'))
    check('B_second_remainder',R('1.001')/2<R('.501'))
    check('log_radius_second_remainder',R('1.001')**3/2<R('.502'))
    check('value_exponential_majorant',2*R('1.001')**2<R('2.01'))
    check('core_gradient_coefficient',R(81,16)+R(1)/R('7.9e9')+tiny<R('5.064'))
    check('outer_gradient_coefficient',1+R(2,4)+(3+1/R('7.9e9'))/16+tiny/256<R('1.689'))
    core=[R('12.01')*R('2.01')*(4*(2*U)**2+R('.502')*U**2),
          R('6.1e-8')*R('6.2e-20'),
          R('6.02')*U*LN,
          R('6.02')*R('.501')*22*U**2,
          R('6.02')*R('5.064')*R('3e-7')**2,
          R('4.03')*U**2/R('7.9e9')]
    core_rows=list(map(R,['5.75e-14','4e-27','1.60e-13','9.57e-15','2.75e-12','3e-25']))
    outer=[R('12.01')*R('2.01')*(R('1e-8')**2/4**3+R('.502')*R('2e-9')**2/4**4),
           core[1],R('6.02')*R('2e-9')/16*LN,
           R('6.02')*R('.501')*22*(R('2e-9')/16)**2,
           R('6.02')*R('1.689')*R('1.4e-8')**2/4**3,core[-1]]
    outer_rows=list(map(R,['5e-17','4e-27','1.66e-15','1.1e-18','3.12e-17','3e-25']))
    for region,vals,uppers,total in [('core',core,core_rows,R('3e-12')),('outer',outer,outer_rows,R('1.8e-15'))]:
        for i,(v,upper) in enumerate(zip(vals,uppers)):
            check(region+'_remainder_row_'+str(i),v<upper,str(v)+' < '+str(upper))
        check(region+'_rounded_sum',sum(uppers)<total,sum(uppers))
    check('actual_scalar_curvature_bound',R('1.3e-21')+R('6.2e-20')+R('3e-12')<R('3.1e-12'))
    # Actual sectional pullbacks, with separate low derivative inputs.
    exp=R('1.001'); gap=R('6.1e-8')
    K=exp*(R('2.31e-6')+R('3.05e-7')+R('10.5')*gap)
    E2=exp*(R('2.31e-6')+22*gap)
    AR=exp*(R('6e-6')+11*gap)
    R2=exp*(R('3.2e-6')+5*gap)
    check('sectional_K',K<R('3.3e-6'))
    check('sectional_E2',E2<R('3.7e-6'))
    check('sectional_AR',AR<R('7e-6'))
    check('sectional_R2',R2<R('3.6e-6'))
    check('sectional_radial_sphere',E2+AR+R2<R('1.5e-5'))
    check('sectional_sphere_sphere',R2+gap/R('7.9e9')<R('3.6e-6'))
    check('Ricci_Frobenius',2*3*R('1.5e-5')<=R('9e-5'))
    check('finite_J',R('9e-5')+(R('3.1e-12')+R('1.3e-21'))/2<R('.001'))
    check('finite_J_global',7+R('.001')<8)
    check('local_scalar_coefficient',R('3.1e-12')/36+R(11,9)*8<10)
    g=R(1,2000**10); b0=R('2.6e-13')
    check('local_auxiliary_block',10*g/180+R('3.1e-12')/12<b0)
    # sqrt(b0)<5.1e-7, no floating-point root needed.
    check('local_sqrt_bound',b0<R('5.1e-7')**2)
    check('conditional_positive_elimination',(1+R('5.1e-7')/2)**2/(1-b0)<R('1.000002'))
    check('scalar_value_does_not_control_Box',S.diff(S.sin(S.Symbol('z')*S.Symbol('N')),S.Symbol('z'),2).has(S.Symbol('N')**2))


def independent_warped_fixtures():
    import mpmath as mp
    mp.mp.dps=90
    a=lambda x:x*(1+x*x/7+x**4/101)
    r=lambda x:2+x*x/5+x**4/103
    u=lambda x:mp.mpf('.02')*x*x/(1+x*x)
    V=lambda x:mp.mpf('.03')/(1+x*x)
    nu=lambda x:(1-1/r(x))*V(x)
    def geom(x,e):
        rr=lambda y:1+(r(y)-1)*mp.exp(e*V(y))
        ss=lambda y:mp.log(rr(y)/r(y))
        ww=lambda y:e*u(y)+ss(y)
        aa=lambda y:a(y)*mp.exp(ww(y))
        wp=mp.diff(ww,x); lapse2=mp.exp(2*ww(x))
        kr=-(mp.diff(aa,x,2)-wp*mp.diff(aa,x))/(lapse2*aa(x))
        ks=-(mp.diff(rr,x,2)-wp*mp.diff(rr,x))/(lapse2*rr(x))
        kt=-mp.diff(aa,x)*mp.diff(rr,x)/(lapse2*aa(x)*rr(x))
        kss=1/rr(x)**2-mp.diff(rr,x)**2/(lapse2*rr(x)**2)
        return 2*(kr+2*ks+2*kt+kss)
    for ix,x in enumerate([mp.mpf(1)/5,mp.mpf(2)/3,mp.mpf(7)/4]):
        AA=mp.diff(a,x)/a(x); RR=mp.diff(r,x)/r(x)
        L=lambda f:mp.diff(f,x,2)+(AA+2*RR)*mp.diff(f,x)
        curv=geom(x,mp.mpf(0)); dot=mp.diff(lambda e:geom(x,e),mp.mpf(0))
        linear=-2*(u(x)+nu(x))*curv-2*(mp.diff(u,x,2)+AA*mp.diff(u,x))-6*L(nu)+4*u(x)/r(x)**2
        check('fixture_linear_'+str(ix),abs(dot-linear)<mp.mpf('1e-70'),mp.nstr(abs(dot-linear),8))
        for ie,e in enumerate([mp.mpf('.2'),mp.mpf('.7'),mp.mpf(1)]):
            k=mp.exp(e*V(x)); tau=1-1/r(x); TT=1+tau*(k-1); ss=mp.log(TT); w=e*u(x)+ss
            BB=k-1-e*V(x)*k; ee=mp.exp(-2*w)
            lt=L(lambda y:1-1/r(y)); vy2=mp.diff(V,x)**2
            term=-6*ee*e**2*k/TT*tau*vy2
            rem=curv*(ee*(1+2*e*(u(x)+nu(x)))-1)+e*(ee-1)*dot
            rem-=6*ee*(e*(k/TT-1)*L(nu)+(BB/TT)*lt+e**2*k/TT*tau*vy2)
            rem+=2*ee*(mp.exp(2*e*u(x))-1-2*e*u(x))/r(x)**2
            target=geom(x,e)-curv-e*dot
            check('fixture_remainder_'+str(ix)+'_'+str(ie),abs(target-rem)<mp.mpf('1e-70'),mp.nstr(abs(target-rem),8))
            check('fixture_omitted_gradient_'+str(ix)+'_'+str(ie),abs(term)>mp.mpf('1e-12'),mp.nstr(abs(term),8))


def main():
    sys.set_int_max_str_digits(0)
    sys.stdout.reconfigure(encoding='utf-8',newline='\n')
    sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths(); hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes(); error=None
    print(json.dumps({'event':'start','story':'S167 scalar curvature remainder','input_count':len(ps)}),flush=True)
    try:
        predecessors(); scalar_algebra(); budgets(); independent_warped_fixtures()
    except Exception:
        error=traceback.format_exc(); print(error,file=sys.stderr,flush=True)
    after=hashes(); good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 scalar curvature remainder','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','python':sys.version,'sympy':S.__version__,
         'protocol':PROTOCOL,'input_sha256_before':before,'input_sha256_after':after,
         'unchanged_inputs':before==after,'checks':ROWS,'exception':error,
         'scalar_curvature_remainder_bounded':good,'local_scalar_auxiliary_block_bounded':good,
         'complete_trace_defect_bounded':False,'mixed_coefficient_estimate_complete':False,
         'full_candidate_operator_realized':False,'full_candidate_Hessian_inverse':False,
         'signed_compatibility_evaluated':False,'nonlinear_family_preserved':False,
         'S167_complete':False,'goal_complete':False}
    if not good:
        run['failed_checker_source_utf8']=(ROOT/OWN[0]).read_text()
        run['failed_note_source_utf8']=(ROOT/OWN[1]).read_text()
    dest=ROOT/'evidence'/(STEM+'.json')
    data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run); dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':
    raise SystemExit(main())
