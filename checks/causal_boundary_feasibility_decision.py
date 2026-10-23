"""Boundary metric probe and feasibility arithmetic, not a quantum solution."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback
import sympy as sp
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='causal-boundary-feasibility-decision'
PREVIOUS='causal-shape-pole-feasibility'
OWN=['notes/'+STEM+'.md','checks/causal_boundary_feasibility_decision.py','evidence/'+STEM+'-review.json']
ROWS=[]; FIXTURES=[]

def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok: raise AssertionError(name)

def geometry():
    v,r,th,ph=sp.symbols('v r theta phi',real=True)
    coords=[v,r,th,ph]; F=sp.Function('F')(v,r)
    g=sp.Matrix([[-F,1,0,0],[1,0,0,0],[0,0,r*r,0],[0,0,0,r*r*sp.sin(th)**2]])
    gi=g.inv(); n=4
    Gamma=[[[sp.simplify(sum(gi[a,d]*(sp.diff(g[d,c],coords[b])+sp.diff(g[d,b],coords[c])-sp.diff(g[b,c],coords[d])) for d in range(n))/2) for c in range(n)] for b in range(n)] for a in range(n)]
    Ric=sp.Matrix(n,n,lambda a,b:sp.simplify(sum(sp.diff(Gamma[c][a][b],coords[c])-sp.diff(Gamma[c][a][c],coords[b])+sum(Gamma[c][c][d]*Gamma[d][a][b]-Gamma[c][b][d]*Gamma[d][a][c] for d in range(n)) for c in range(n))))
    R=sp.simplify(sp.trace(gi*Ric)); Ein=sp.simplify(Ric-g*R/2); mixed=sp.simplify(gi*Ein)
    expectedR=-sp.diff(F,r,2)-4*sp.diff(F,r)/r+2*(1-F)/r**2
    check('four_dimensional_scalar',sp.simplify(R-expectedR)==0)
    check('ingoing_rr_constraint',Ein[1,1]==0)
    check('ingoing_mass_radial_row',sp.simplify(mixed[0,0]-(r*sp.diff(F,r)+F-1)/r**2)==0)
    check('ingoing_mass_flux_row',sp.simplify(mixed[1,0]+sp.diff(F,v)/r)==0)
    G,e,B,lam,cR=sp.symbols('G eta B lambda c_R',positive=True)
    m=sp.Function('m')(v); F0=sp.Function('F0')(r); family=F0-2*G*e*m/r
    vary=lambda expr:sp.simplify(sp.diff(expr.subs(F,family).doit(),e).subs(e,0))
    check('probe_scalar_variation_zero',vary(R)==0)
    check('probe_radial_mass_row_zero',vary(mixed[0,0])==0)
    check('probe_flux_nonzero',sp.simplify(vary(mixed[1,0])-2*G*sp.diff(m,v)/r**2)==0)
    h=2*G*m/r
    check('fixed_conformal_boundary_source',sp.limit(B**2*h/r**2,r,sp.oo)==0)
    check('horizon_regular_metric_tangent',h.subs(r,B)==2*G*m/B)
    check('nonzero_areal_invariant',vary(gi[1,1])==-2*G*m/r)
    R0=-12*lam; Fl=1+32*sp.pi*G*cR*R0
    check('local_R2_flux_multiplier',sp.expand(Fl)==1-384*sp.pi*G*cR*lam)
    A=r*r*sp.diff(F,r,2)-2*r*sp.diff(F,r)+2*F-2
    check('Weyl_density_probe_not_omitted',sp.simplify(vary(A**2/(3*r**4))+8*G*m*A.subs(F,F0).doit()/r**5)==0)
    # This tests the stated invariant formula algebra, not a full quantum Hessian.
    alpha=sp.Matrix([[0,-sp.I],[sp.I,0]]); beta=sp.Matrix([[0,1],[1,0]])
    for sign in [-1,1]:
        P=(sp.eye(2)+sign*beta)/2
        check('projector_flux_zero_'+str(sign),P*alpha*P==sp.zeros(2))
    gamma=sp.symbols('gamma',positive=True); x=sp.symbols('x',real=True)
    check('constraint_boundary_lift_norm_squared',sp.integrate(sp.exp(2*gamma*x)/2,(x,-sp.oo,0))==1/(4*gamma))
    FIXTURES.append(dict(fixture='direct_4d_curvature_and_boundary_probe',scalar=str(R),mass_flux=str(mixed[1,0]),quantum_response_evaluated=False))

def budgets_and_lapse():
    rho=Q('1e-8'); Y=Q('5e-9'); Z0=Q('.25'); Z2=Q('5e7')
    check('preserved_ball_budget',Y+Z0*rho+Z2*rho*rho/2==rho)
    check('contraction_margin',Z0+Z2*rho==Q('.75'))
    check('rounded_down_middle_target',Q('61.2745098039')*Q('1.632e-11')<=Q('1e-9'))
    check('rounded_down_cold_target',Q('1.4705882352')*Q('1.36e-9')<=Q('2e-9'))
    check('scalar_trace_only_error',Q(4,3)*10**18*Q('1.1e-31')<Q('1.467e-13'))
    check('shape_to_trace_allocation',Q(4,3)*10**18*(Q('1.1e-31')+Q('1.25014e-27'))==Q('1.667e-9'))
    check('bolt_reference_plus_difference',Q('3.33e-10')+Q('1.667e-9')==Q('2e-9'))
    check('higher_admission_target',Q('1.3e-4')+rho*87000==Q('1e-3'))
    check('cold_admission_if_unit_embedding',Q('4.9e-7')+rho==Q('5e-7'))
    check('outer_admission_if_unit_embedding',Q('1.4e-7')+rho==Q('1.5e-7'))
    check('old_interior_certificate_fails',Q('4.141944868545865e118')/Q('2.38393e-25')>Q('1.7374e143'))
    check('old_global_certificate_fails',Q('1.069626977892751e201')/Q('2.38393e-25')>Q('4.4868e225'))
    check('local_flux_not_removed_by_R2',1-Q('1.179e-8')>Q('.99999998'))
    z=sp.symbols('z',positive=True); k=sp.symbols('k')
    check('outer_redshift_limit',sp.limit(z*(sp.log(z)+k),z,0,dir='+')==0)
    check('horizon_redshift_term_unbounded',sp.limit(z*sp.log(z),z,sp.oo)==sp.oo)
    # Numerical illustrations of the exact limit, with the retained model constants.
    mp.mp.dps=80; g0=mp.mpf(2000)**-10
    g=2*g0/(1+mp.sqrt(1+11*g0/(90*mp.pi)))
    x0=(mp.sqrt(1+12*mp.pi*mp.mpf(2000)**6)-1)/6
    Fstar=1-384*mp.pi*(mp.mpf('1e22')-1/(5760*mp.pi**2))*g
    a=20*mp.pi*Fstar; zeta=g0/x0*mp.exp(mp.mpf(12)/5)
    om=mp.sqrt(a/mp.lambertw(-a/zeta,-1))
    D=lambda n:1-om**2/(a*n*n)*mp.log(-om**2/(zeta*n*n))
    check('retained_frozen_root_for_lapse_test',abs(D(1))<mp.mpf('1e-70'))
    check('outer_local_comparison_exceeds_transfer_target',abs(D(mp.mpf('1e6')))>mp.mpf('.99'))
    check('near_horizon_local_comparison_large',abs(D(mp.mpf('1e-6')))>mp.mpf('1e10'))
    FIXTURES.append(dict(fixture='rational_budgets_and_lapse_diagnostic',outer_example=mp.nstr(D(mp.mpf('1e6')),35),horizon_example=mp.nstr(D(mp.mpf('1e-6')),35),global_radial_mode_evaluated=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    try:
        check('previous_causal_model_run_passed',old['status']=='passed' and old['frozen_shape_causal_poles_enclosed'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        geometry(); budgets_and_lapse()
    except Exception: error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),kinematic_boundary_implication_rejected=error is None,uniform_frozen_lapse_comparison_rejected=error is None,actual_retarded_characteristic_reduction_proved=False,actual_causal_pole_transfer_bound_proved=False,actual_finite_geometry_growing_mode_proved=False,physical_instability_proved=False,physical_stability_proved=False,semiclassical_validity_proved=False,theory_or_boundary_prescription_changed=False,full_quantum_shape_defect_evaluated=False,full_defect_bolt_pairing_evaluated=False,nonlinear_solution_constructed=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__': main()
