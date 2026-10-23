"""Frozen bulk-equation/symbol audit; no actual metric inverse is certified."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback
import mpmath as mp
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
STEM='physical-metric-block-and-preconditioner'
PREVIOUS='canonical-matter-neighborhood'
OWN=['notes/'+STEM+'.md','checks/physical_metric_block_and_preconditioner.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def algebra():
    P,u,b,C,d,R,Qm=sp.symbols('P u b C d R Qm',real=True)
    fields=[P,u,b]
    grads=[sp.symbols(n+'x '+n+'y',real=True) for n in ['P','u','b']]
    hess=[sp.symbols(n+'xx '+n+'yy',real=True) for n in ['P','u','b']]
    gp,gu,gb=grads;hp,hu,hb=hess
    w=sp.Function('w')(u)
    dot=lambda x,y:sum(a*c for a,c in zip(x,y))
    e=sp.exp(2*P);U=-2*e/C+(w-6*d)*sp.exp(4*P)+Qm
    L=-e*R-4*e*dot(gp,gb)-6*e*dot(gp,gp)+6*e*dot(gu,gu)+sp.exp(2*b)*U
    def total(a,i):
        # The differentiated momenta have no opposite-coordinate jets.
        return sum(sp.diff(a,f)*g[i]+sp.diff(a,g[i])*h[i] for f,g,h in zip(fields,grads,hess))
    EL=[sp.diff(L,f)-sum(total(sp.diff(L,g[i]),i) for i in range(2)) for f,g in zip(fields,grads)]
    Eb=4*e*(sum(hp)+2*dot(gp,gp))+2*sp.exp(2*b)*U
    EP=e*(-2*R+4*sum(hb)+12*sum(hp)+12*dot(gp,gp)+12*dot(gu,gu))+sp.exp(2*b)*sp.diff(U,P)
    Eu=-12*e*(sum(hu)+2*dot(gp,gu))+sp.exp(2*b)*sp.diff(U,u)
    for name,actual,target in zip(['P','u','b'],EL,[EP,Eu,Eb]):
        check('two_dimensional_Euler_row_'+name,sp.simplify(actual-target)==0)
    raw=-e*R+2*e*sum(hb)+6*e*sum(hp)+6*e*dot(gp,gp)+6*e*dot(gu,gu)+sp.exp(2*b)*U
    flux=sum(total(2*e*gb[i]+6*e*gp[i],i) for i in range(2))
    check('reduction_integration_by_parts',sp.simplify(raw-L-flux)==0)
    G,c,lam,A2=sp.symbols('G c lam A2',positive=True)
    F=1-384*sp.pi*G*c*lam
    W=(1-sp.exp(-2*u))**2/(64*sp.pi*G*c)
    phi=(F*sp.exp(2*u)-1)/(32*sp.pi*G*c)
    check('vacuum_potential_identity',sp.simplify(4*W+sp.diff(W,u)-24*lam-2*sp.exp(-2*u)*phi)==0)
    Rg=sp.exp(-2*(P-u))*(sp.exp(-2*b)/A2*(R-2*sum(hb)-6*(sum(hp)-sum(hu))-6*dot([x-y for x,y in zip(gp,gu)],[x-y for x,y in zip(gp,gu)]))+2/(C*A2))
    combined=(EP+Eu).subs({sp.diff(w,u):A2*sp.diff(W,u),w:A2*W,d:lam*A2},simultaneous=True)
    check('exact_curvature_constraint_row',sp.simplify(combined+2*A2*sp.exp(4*P-2*u+2*b)*(Rg-phi))==0)
    f=sp.Function('f')(P-u,b)
    check('complete_quantum_row_cancellation',sp.simplify(sp.diff(f,P)+sp.diff(f,u))==0)
    q,gr,B2=sp.symbols('q gr B2',positive=True)
    magnetic=4*sp.pi*A2*B2*q**2/(8*gr**2*B2**2)
    check('magnetic_flux_normalization',sp.simplify(magnetic/(F*B2/(4*G))-2*sp.pi*G*A2*q**2/(F*gr**2*B2**2))==0)
    T,v,s,ell=sp.symbols('T v s ell',real=True)
    kin=-6*s*s-4*s*b+6*v*v
    check('local_kinetic_diagonalization',sp.expand(kin.subs(s,T-b/3)-(-6*T*T+sp.Rational(2,3)*b*b+6*v*v))==0)
    # Derive (4) from the source invariant formula, rather than only
    # comparing the displayed final scalar-square expressions.
    z=sp.symbols('z')
    metric=sp.diag(2*s+2*b,2*s+2*b,2*s,2*s)
    projector=sp.diag(1,0,1,1)
    t2=sp.trace((metric*projector)**2);t1=sp.trace(metric*projector)
    Ws=sp.exp(-2*z*ell)*((2*z-3)*t2+t1*t1)/(16*sp.pi**2)
    H=(sp.diff(Ws,z).subs(z,0)+sp.Rational(46,15)*Ws.subs(z,0))/120
    target=(s+b/3)**2/(80*sp.pi**2)+b*b*(ell-sp.Rational(6,5))/(120*sp.pi**2)
    check('Dirac_symbol_from_projected_tensor',sp.simplify(H-target)==0)
    dc=1/(5760*sp.pi**2)
    hstar=H.subs(s,T-v-b/3)-2*dc*(6*(T-v))**2
    check('anomaly_complete_fourth_order_cancellation',sp.simplify(hstar-b*b*(ell-sp.Rational(6,5))/(120*sp.pi**2))==0)
    tau,V,zeta=sp.symbols('tau V zeta',positive=True)
    ret=sp.log(1+(sp.exp(zeta)-1)/tau)
    check('radius_return_linearization',sp.diff(ret,zeta).subs(zeta,0)==1/tau)
    t,z0,F0=sp.symbols('t z0 F0',positive=True)
    ratio=t*(sp.log(t/z0)/2-sp.Rational(6,5))/(10*sp.pi*F0)
    check('shape_ratio_minimum_location',sp.simplify(sp.diff(ratio,t).subs(t,z0*sp.exp(sp.Rational(7,5))))==0)
    check('shape_ratio_minimum_value',sp.simplify(ratio.subs(t,z0*sp.exp(sp.Rational(7,5)))+z0*sp.exp(sp.Rational(7,5))/(20*sp.pi*F0))==0)
    check('shape_ratio_strict_convexity',sp.simplify(sp.diff(ratio,t,2)-1/(20*sp.pi*F0*t))==0)
    check('nonlinear_preservation',Q('5e-9')+Q('.25')*Q('1e-8')+Q('5e7')*Q('1e-8')**2/2==Q('1e-8'))
    check('nonlinear_contraction',Q('.25')+Q('5e7')*Q('1e-8')==Q('.75')<1)
    check('regional_return_budget',Q('1.632e-11')*Q('61.2745098039')+Q('1.36e-9')*Q('1.4705882352')+Q('2e-9')<=Q('5e-9'))
    for name,background,E,limit in [('interior','1.3e-4','87000','1e-3'),('cold','4.9e-7','1','5e-7'),('outer','1.4e-7','1','1.5e-7')]:
        check('admission_'+name,Q(background)+Q('1e-8')*Q(E)==Q(limit))
    check('residual_error_target',Q('1.25e-9')/(1-Q('.75'))==Q('5e-9'))
    FIXTURES.append(dict(fixture='bulk_variation_and_local_symbol',symbolic_two_dimensional_jets=True,physical_boundary_domain_certified=False))

def constants(ctx):
    n=ctx.mpf(2000);g0=n**-10;c=ctx.mpf('1e22')-1/(5760*ctx.pi**2)
    g=2*g0/(1+ctx.sqrt(1+4*ctx.mpf(11)/(360*ctx.pi)*g0))
    F=1-384*ctx.pi*c*g;x0=(ctx.sqrt(1+12*ctx.pi*n**6)-1)/6
    return F,g0/x0,1/(96*ctx.pi*c),g

def numerical():
    records=[]
    for digits in [70,110]:
        mp.mp.dps=digits;F,z,m2,g=constants(mp)
        def r(t):return t*(mp.log(t/z)/2-mp.mpf(6)/5)/(10*mp.pi*F)
        closed=20*mp.pi*F/mp.lambertw(20*mp.pi*F/(z*mp.exp(mp.mpf(12)/5)))
        # Independent root evaluation and sign bracket. Interval signs below
        # certify the printed bracket and the complete u_b allowance.
        independent=mp.findroot(lambda t:r(t)-1,(mp.mpf('.6'),mp.mpf('.7')))
        check('independent_crossover_'+str(digits),abs(closed-independent)<mp.mpf('1e-60'))
        check('crossover_equation_'+str(digits),abs(r(closed)-1)<mp.mpf('1e-60'))
        check('Planck_ratio_'+str(digits),mp.mpf('1.5348')<r(1)<mp.mpf('1.5349'))
        check('scalar_mass_ratio_'+str(digits),mp.mpf('3.39e8')<m2/g<mp.mpf('3.40e8'))
        check('scalar_mass_quantum_ratio_'+str(digits),0<r(m2)<mp.mpf('2.13e-25'))
        rec=dict(digits=digits,Fstar=mp.nstr(F,50),Gmu0_squared=mp.nstr(z,50),
            crossover_Gk_squared=mp.nstr(closed,50),crossover_k_sqrtG=mp.nstr(mp.sqrt(closed),50),
            ratio_at_Gk_squared_one=mp.nstr(r(1),50),ratio_at_Gk_squared_one_hundredth=mp.nstr(r(mp.mpf('.01')),50),
            scalar_mass_squared_G=mp.nstr(m2,50),scalar_mass_squared_over_lambda=mp.nstr(m2/g,50),
            ratio_at_scalar_mass=mp.nstr(r(m2),50),negative_model_ratio=mp.nstr(z*mp.exp(mp.mpf(7)/5)/(20*mp.pi*F),50))
        records.append(rec)
    for key in [k for k in records[0] if k!='digits']:
        check('cross_precision_'+key,abs(mp.mpf(records[0][key])/mp.mpf(records[1][key])-1)<mp.mpf('1e-45'))
    iv=mp.iv;iv.dps=70;F,z,m2,g=constants(iv)
    z=z*iv.exp(iv.mpf(['-1.08e-36','1.08e-36']))
    low=iv.mpf('.65442489137407894964960401158');high=iv.mpf('.65442489137407894964960401159')
    r=lambda t:t*(iv.log(t/z)/2-iv.mpf(6)/5)/(10*iv.pi*F)
    check('interval_crossover_lower_all_ub',r(low)<1,str(r(low)))
    check('interval_crossover_upper_all_ub',r(high)>1,str(r(high)))
    negative=z*iv.exp(iv.mpf(7)/5)/(20*iv.pi*F)
    check('interval_model_negative_bound',negative<iv.mpf('7.70e-45'),str(negative))
    check('interval_Fstar_guard',F>iv.mpf('0.99999998821902754'))
    FIXTURES.append(dict(fixture='frozen_symbol_crossover',records=records,
        interval_digits=70,ub_absolute_allowance='5.4e-37',
        actual_low_frequency_operator_computed=False,actual_correction_computed=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_matter_audit_passed',old['status']=='passed' and old['canonical_matter_neighborhood_proved'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        algebra();numerical()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,
        unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        exact_physical_bulk_rows_derived=error is None,mixed_order_principal_symbol_derived=error is None,
        second_order_all_frequency_perturbation_rejected=error is None,
        actual_correction_computed=False,full_metric_Hessian_domain_quantified=False,
        outer_contract_transferred_to_S169=False,physical_coupled_inverse_proved=False,
        nonlinear_solution_constructed=False,semiclassical_validity_proved=False,physical_stability_proved=False,
        component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
