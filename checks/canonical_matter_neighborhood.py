"""Frozen canonical matter-domain budgets; no metric Hessian is inverted."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
from math import factorial
import json, traceback
import mpmath as mp
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
STEM='canonical-matter-neighborhood';PREVIOUS='cold-source-feasibility-audit'
OWN=['notes/'+STEM+'.md','checks/canonical_matter_neighborhood.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def exp_upper(x,n=60):
    assert 0<=x<n+2
    return sum(x**j/factorial(j) for j in range(n+1))+x**(n+1)/factorial(n+1)/(1-x/(n+2))

def rational():
    rho=Q('1e-8');ell1=100*rho;e=exp_upper(rho)
    eta=e-1+e*ell1/36;graph=e*(1+ell1/2)
    check('massive_relative_defect',eta<Q('3.78e-8'),str(eta))
    check('matter_transfer_not_metric_target',eta<Q('.25'))
    check('graph_equivalence',graph<Q('1.000000511'))
    check('LLL_relative_inverse',e-1<Q('1.00000001e-8'))
    check('LLL_relative_square_gap',1/exp_upper(2*rho)>1-Q('2e-8'))
    a1=1+50*e/18;a2=1+100*e/18
    check('first_relative_vertex',a1<Q('3.78'))
    check('second_relative_vertex',a2<Q('6.56'))
    check('second_resolvent_budget',2*Q('3.78')**2+Q('6.56')==Q('35.1368')<Q('35.14'))
    complex_eta=exp_upper(Q('.01'))-1+exp_upper(Q('.01'))/36
    check('radial_complex_graph_tube',complex_eta<Q('.0382'))
    check('complex_inverse_budget',1/(1-Q('.0382'))<Q('1.04'))
    outer=Q('1.01e-7');contract=Q('1.5e-7')
    check('stronger_outer_margin',contract-outer==Q('4.9e-8'))
    check('outer_ball_preserved',outer+rho==Q('1.11e-7')<contract)
    check('outer_background_target',contract-rho==Q('1.4e-7'))
    check('conditional_outer_embedding',(contract-outer)/rho==Q('4.9'))
    check('outer_first_source_increment',Q('1.8e-7')*rho==Q('1.8e-15'))
    check('outer_second_source_remainder',Q('.02')*rho**2/2==Q('1e-18'))
    check('older_outer_contract_different',outer>Q('1e-7') and outer<contract)
    check('metric_source_budget_unchanged',Q('1.632e-11')*Q('61.2745098039')+Q('1.36e-9')*Q('1.4705882352')+Q('2e-9')<=Q('5e-9'))
    # A differential, rather than a generic finite-matrix, check of the
    # first-order product rule used in the graph-relative estimate.
    b,bx,u,ux,v,vx,psi,psix=sp.symbols('b bx u ux v vx psi psix')
    A=sp.exp(-b)*(psix-bx*psi/2)
    first=sp.diff(A,b)*u+sp.diff(A,bx)*ux
    second=sp.diff(first,b)*v+sp.diff(first,bx)*vx
    check('local_conformal_first_vertex',sp.simplify(first+u*A+sp.exp(-b)*ux*psi/2)==0)
    check('local_conformal_second_vertex',sp.simplify(second-u*v*A-sp.exp(-b)*(u*vx+v*ux)*psi/2)==0)
    check('local_graph_return',sp.simplify(sp.exp(b)*A+bx*psi/2-psix)==0)
    FIXTURES.append(dict(fixture='exact_graph_and_outer_contract_budgets',eta_upper=str(eta),
        graph_equivalence_upper=str(graph),LLL_relative_inverse_upper=str(e-1),
        first_vertex_upper=str(a1),second_vertex_upper=str(a2),complex_eta_upper=str(complex_eta),
        metric_Hessian_transfer_certified=False,outer_S169_transfer_certified=False))

def adj(a):return a.transpose_conj()
def opnorm2(a):
    s=adj(a)*a
    disc=mp.re((s[0,0]-s[1,1])**2)+4*abs(s[0,1])**2
    return mp.sqrt((mp.re(s[0,0]+s[1,1])+mp.sqrt(disc))/2)

def numerical():
    records=[]
    for digits in [70,110]:
        mp.mp.dps=digits;small=mp.mpf('1e-30');t=mp.mpf('1e-8');ii=mp.j
        A=ii*mp.matrix([[1+small,1-small],[1-small,1+small]])/2
        R0=-ii*mp.matrix([[1+1/small,1-1/small],[1-1/small,1+1/small]])/2
        E=mp.diag([mp.exp(t/2),mp.exp(-t/2)]);Em=mp.diag([mp.exp(-t/2),mp.exp(t/2)])
        Ab=Em*A*Em;Rb=mp.inverse(Ab);closed=E*R0*E
        check('matrix_skew_adjoint_'+str(digits),opnorm2(Ab+adj(Ab))<mp.mpf('1e-60'))
        check('independent_massless_inverse_'+str(digits),opnorm2(Rb-closed)/opnorm2(closed)<mp.mpf('1e-32'))
        relative=opnorm2(Rb-R0)/opnorm2(R0)
        defect=opnorm2((Ab-A)*R0)
        check('small_gap_generic_test_fails_'+str(digits),defect>mp.mpf('1e21'))
        check('exact_inverse_relative_bound_'+str(digits),relative<mp.expm1(t))
        # Finite congruence validates the endpoint derivative algebra only;
        # a generic matrix has no first-order Clifford product rule.
        U=mp.diag([1,-1]);V=mp.diag([2,3]);m=mp.mpf(18)
        Qm=mp.inverse(Ab+m*mp.eye(2))
        check('massive_resolvent_gap_'+str(digits),opnorm2(Qm)<=1/m)
        check('spectral_graph_factor_'+str(digits),opnorm2(Ab*Qm)<=1)
        endpoint=(U*closed+closed*U)/2
        def exact(s):
            Es=mp.diag([mp.exp((t+s)/2),mp.exp(-(t+s)/2)])
            return Es*R0*Es
        numeric=mp.matrix([[mp.diff(lambda ss:exact(ss)[i,j],0) for j in range(2)] for i in range(2)])
        check('LLL_endpoint_first_derivative_'+str(digits),opnorm2(numeric-endpoint)/opnorm2(R0)<mp.mpf('1e-60'))
        endpoint2=(U*V*closed+U*closed*V+V*closed*U+closed*U*V)/4
        def two(s,z):
            Ess=mp.diag([mp.exp((t+s+2*z)/2),mp.exp((-t-s+3*z)/2)])
            return Ess*R0*Ess
        numeric2=mp.matrix([[mp.diff(lambda ss,zz:two(ss,zz)[i,j],(0,0),(1,1)) for j in range(2)] for i in range(2)])
        check('LLL_endpoint_mixed_derivative_'+str(digits),opnorm2(numeric2-endpoint2)/opnorm2(R0)<mp.mpf('1e-60'))
        records.append(dict(digits=digits,relative_inverse_change=mp.nstr(relative,40),
            generic_Neumann_defect=mp.nstr(defect,40),reference_gap_parameter=str(small)))
    for key in ['relative_inverse_change','generic_Neumann_defect']:
        check('cross_precision_'+key,abs(mp.mpf(records[0][key])/mp.mpf(records[1][key])-1)<mp.mpf('1e-30'))
    FIXTURES.append(dict(fixture='noncommuting_congruence_with_small_gap',records=records,
        finite_matrix_is_not_charged_Dirac_kernel=True,actual_LLL_absolute_gap_computed=False,
        continuum_boundary_domain_machine_verified=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_cold_audit_passed',old['status']=='passed' and old['normalized_cold_canonical_source_bound']=='1.36e-9 a_s')
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        rational();numerical()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,
        unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        canonical_matter_neighborhood_proved=error is None,
        massive_matter_relative_defect_bound='3.78e-8' if error is None else None,
        LLL_relative_inverse_change_bound='1.00000001e-8' if error is None else None,
        LLL_absolute_gap_computed=False,stronger_outer_contract_restored=True,
        outer_contract_transferred_to_S169=False,full_metric_Hessian_domain_quantified=False,
        physical_coupled_inverse_proved=False,nonlinear_solution_constructed=False,
        semiclassical_validity_proved=False,physical_stability_proved=False,
        component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
