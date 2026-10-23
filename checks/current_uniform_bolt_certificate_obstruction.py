"""Actual source projection and uniform-certificate test; not a metric solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,traceback
import sympy as sp
import mpmath as mp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_geodesic_flow_multiplier as previous
import current_finite_local_action_assembly as local
york=previous.york;prior=previous.prior;endpoint=previous.endpoint
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1]
STEM='current-uniform-bolt-certificate-obstruction';PREVIOUS='current-geodesic-flow-multiplier'
OWN=['notes/'+STEM+'.md','checks/current_uniform_bolt_certificate_obstruction.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];FIXTURES=[];encoded=prior.encoded;Poly=local.Poly;zpow=local.monomial
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)
prior.check=check;york.check=check

def algebra():
    t,z,u,d,delta=sp.symbols('t z u d delta',positive=True)
    f=(t-1)**2*(1+2/t+(3+u)/t**2+d*u/t**4)-delta/t
    pp=1-2*z-5*z*z+12*z**3+u*(-3*z*z+4*z**3)+d*u*(-7*z**4+8*z**5)
    check('actual_thermal_displacement_supersolution',sp.simplify(2*(t-1)*sp.diff(f,t)-3*f-(t-1)**2*pp.subs(z,1/t)-delta*(5*t-2)/t**2)==0)
    w=sp.symbols('w');base=1-2*z-5*z*z+12*z**3;cells=[]
    for i in range(8):
        a=sp.Rational(i,8);b=sp.Rational(i+1,8)
        poly=sp.Poly(base.subs(z,a+(b-a)*w),w)
        bern=[sum(poly.nth(j)*sp.binomial(k,j)/sp.binomial(3,j) for j in range(k+1)) for k in range(4)]
        check('displacement_Bernstein_'+str(i),min(bern)>sp.Rational(1,8))
        cells.append(dict(left=str(a),right=str(b),coefficients=list(map(str,bern))))
    check('finite_seed_displacement_perturbations',sp.Rational(7,8*10**9)+15*sp.Rational(86,10**20)/sp.Integer(8*10**9)<sp.Rational(1,8))
    FIXTURES.append(dict(kind='exact_displacement_Bernstein_cells',cells=cells))
    A=-10*z**4-20*z**5+30*z**6
    check('signed_leading_trace_factorization',sp.expand(A+10*z**4*(1-z)*(1+3*z))==0)
    check('leading_absolute_integral',-sp.integrate(A,(z,0,1))==sp.Rational(22,21))
    check('leading_weighted_absolute_integral',-sp.integrate((1-z)*A,(z,0,1))==sp.Rational(23,84))
    eta,M,x,y=sp.symbols('eta M x y',positive=True)
    f=eta*x*y/M**2+(1-y+eta*y/M)*(1-x+eta*x/M)
    check('all_trial_minimum_corner',sp.simplify(f.subs({x:1,y:1})-eta*(1+eta)/M**2)==0)
    check('one_residual_corner',sp.simplify(f.subs({x:0,y:1})-eta/M)==0)
    check('two_residual_corner',f.subs({x:0,y:0})==1)
    check('bilinear_minimization',sp.diff(f,x,2)==0 and sp.diff(f,y,2)==0)
    check('operator_norm_upper',(sp.Rational(116,100)-sp.Rational(22,10000))**2>sp.Rational(1338,1000))
    check('minimum_corner_order',sp.Rational(116,100)>1+sp.Rational(22,10000))
    # Scalar blocks of the actual symmetric weak equations.
    aa,bb,dd,rs,rt,ls,lt,et,pt=sp.symbols('aa bb dd rs rt ls lt et pt')
    es=(rs+bb*et)/aa;ps=(bb*pt-ls)/aa
    check('Schur_forward_scalar_equation',sp.expand(-aa*es+bb*et+rs)==0)
    check('Schur_adjoint_scalar_equation',sp.expand(-aa*ps+bb*pt-ls)==0)
    check('Schur_forward_TT_equation',sp.expand(bb*es+dd*et+rt-((dd+bb*bb/aa)*et+rt+bb*rs/aa))==0)
    check('Schur_effective_point_functional',sp.expand(ls*es+lt*et-(ls*rs/aa+(lt+bb*ls/aa)*et))==0)

def abs_integral(poly,weighted=False):
    total=arb(0)
    for (k,l),c in poly.c.items():
        assert k>=0 and l==0
        total+=c.abs_upper()/((k+1)*(k+2) if weighted else (k+1))
    return total

def execute():
    p,h,zh,env=york.actual_parameters();NUMERICS.extend(york.NUMERICS)
    Q=Poly({(power,0):env[sym] for power,sym in zip(york.powers,york.qs)})
    C,B2,lam=p['C'],p['B2'],p['g'];q=arb(p['q']);x=lam*B2;d=p['eps']*C;u=1/x;du=d/x
    pi=arb.pi();ee=1/(60*pi);D=96*pi*arb('1e22')-ee;eta=ee*lam;dl=D*lam
    S=ee*q*q/(lam*B2*B2);ff=Q*(C/x)*zpow(-2)
    L=lambda f:zpow(4)*(ff*f.derivative()).derivative()
    r=4*zpow(5)-6*zpow(6);rn=-12+du*r
    tn=(3+u+du)*zpow(4)-du*(5*zpow(5)-arb(9)/2*zpow(6))
    cn=-zpow(3)*ff.derivative()+rn/6
    P=7*cn**2+arb(22)/3*tn**2-arb(11)/72*(rn**2-144)
    Ps=2*(6+du)**2+arb(22)/3*(6+du)*u+arb(22)/3*u*u+22
    Bp=P+Ps*r/2;A=-10*zpow(4)-20*zpow(5)+30*zpow(6)
    remainder=eta*Bp-dl*du*L(r);tau=S*A+remainder
    check('product_root_correlated_interval',du.overlaps(-5*S+eta*Ps/2))
    Kt=arb(11)/36*rn-arb(7)/3*cn+arb(11)/3*tn
    known=dl*L(tau)+eta*Kt*tau;known_value=known.integral(zh)
    I0=S*arb(22)/21+abs_integral(remainder)
    I1=S*arb(23)/84+abs_integral(remainder,True)
    sr=arb('3.33e-10')*h;amax=arb('3.3e-10');Hmax=arb('1.08e-9')
    areal=(490*eta+316*dl*du.abs_upper())*arb('2.11e-6')+(535*eta+56*dl*du.abs_upper())*arb('2.01e-7')
    density=(amax+4*sr)*I0+8*Hmax*I1+sr*tau.value(zh).abs_upper()
    total_error=areal+density
    Ilower=(x/C)*(known_value-total_error);Iupper=(x/C)*(known_value+total_error)
    check('positive_signed_actual_trace_projection',Ilower>arb('1e-38'),encoded(Ilower))
    scalar=zpow(4);di=dict(nu=scalar,u=Poly(0),scalar=True)
    energy=sum(local.scalar_covariant(p,Q,scalar,scalar).values(),Poly(0));A44=energy.integral(zh)
    parts,flux,ident=local.action_bilinear(p,Q,di,di)
    check('independent_complete_scalar_energy',A44.overlaps(sum(parts.values(),Poly(0)).integral(zh)))
    check('actual_scalar_energy_positive',A44>0 and A44<arb('1.334'))
    M=arb('1.16');e0=arb('.0022');Wabs=pi*B2/p['gamma'];ell=(1+h)/h
    EP=Ilower.lower()*ell*zh**4/(M*M*A44.upper())
    Rmin=Wabs.sqrt()*Ilower.lower()/A44.upper().sqrt()
    Lmin=ell*zh**4/(Wabs*A44.upper()).sqrt()
    check('norm_product_identity',EP.overlaps(Rmin*Lmin/(M*M)))
    check('actual_norm_product_lower_test',EP>arb('1e-5'),encoded(EP))
    floor=e0*(1+e0)*EP;target=arb('1.667e-9')
    check('full_V_uniform_certificate_fails_for_every_pair',floor>target,encoded(floor/target))
    def required_eta(budget):return 2*budget/EP/(1+(1+4*budget/EP).sqrt())
    efull=required_eta(target);ehalf=required_eta(target/2)
    check('necessary_uniform_eta_is_smaller',efull<e0 and ehalf<efull)
    check('necessary_eta_threshold_identity',(efull*(1+efull)*EP).overlaps(target))
    Lh=1+(arb('.01')/h).log();scalar_bolt=(1+h)*Lh*(arb(4)/3)*arb('1e18')*arb('1.1e-31')
    check('actual_scalar_Schur_bolt_budget',scalar_bolt<arb('1.2e-11'),encoded(scalar_bolt))
    TTbudget=(target-scalar_bolt)/2;TTproduct=TTbudget/(e0*(1+e0))
    check('positive_remaining_Schur_budget',TTproduct>0 and TTproduct<arb('3.78032e-7'))
    NUMERICS.append(dict(kind='actual_signed_trace_projection',C=encoded(C),B2=encoded(B2),lambda_value=encoded(lam),
        h=encoded(h),zh=encoded(zh),x=encoded(x),d=encoded(d),du=encoded(du),eta=encoded(eta),dl=encoded(dl),S=encoded(S),
        source_remainder_coefficients=remainder.serialized(),known_integrand_coefficients=known.serialized(),
        known_integral=encoded(known_value),known_integral_moments=known.moment_records(zh),
        I0_upper=encoded(I0),I1_upper=encoded(I1),sr_upper=encoded(sr),tau_at_bolt=encoded(tau.value(zh)),
        unknown_areal_terms_upper=encoded(areal),density_and_bolt_terms_upper=encoded(density),
        I_lower=encoded(Ilower),I_upper=encoded(Iupper),physical_action_abs=encoded(Wabs),
        physical_source_projection_sign='negative because W<0 and I>0',full_shape_source_evaluated=False))
    NUMERICS.append(dict(kind='all_trial_uniform_certificate_feasibility',A44=encoded(A44),energy_coefficients=energy.serialized(),
        energy_moments=energy.moment_records(zh),ell_b=encoded(ell),M=encoded(M),eta=encoded(e0),
        source_dual_norm_lower=encoded(Rmin),bolt_dual_norm_lower=encoded(Lmin),
        ideal_forward_norm_lower=encoded(Rmin/M),ideal_adjoint_norm_lower=encoded(Lmin/M),
        unavoidable_EP_lower=encoded(EP),certificate_expression_lower=encoded(floor),failure_factor=encoded(floor/target),
        necessary_uniform_eta_entire_budget=encoded(efull),necessary_uniform_eta_half_budget=encoded(ehalf),
        all_trial_pair_bound=True,actual_operator_difference_lower_bound=False,physical_error_lower_bound=False))
    NUMERICS.append(dict(kind='exact_scalar_Schur_budget',Lh=encoded(Lh),scalar_bolt_upper=encoded(scalar_bolt),
        remaining_TT_half_budget=encoded(TTbudget),ideal_TT_norm_product_target=encoded(TTproduct),
        actual_TT_norms_unknown=True,Schur_solve_executed=False,full_physical_feasibility_test_executed=False))

def main():
    ctx.prec=256;mp.mp.dps=55;review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|
        {'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};deps=endpoint.dependencies();error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',
        protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_flow_multiplier_passed',old['status']=='passed' and old['actual_form_quadrature_error_bounded'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        algebra();execute()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths};depafter=endpoint.dependencies()
    status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,
        finished_utc=datetime.now(timezone.utc).isoformat(),
        actual_signed_trace_projection_certified=status=='passed',full_V_uniform_certificate_rejected=status=='passed',
        all_trial_pair_obstruction=status=='passed',scalar_Schur_bolt_budget_bounded=status=='passed',
        exact_scalar_Schur_reduction_selected=status=='passed',Schur_forward_adjoint_solve_executed=False,
        actual_physical_correction_computed=False,full_solution_feasibility_test_executed=False,
        metric_admission_proved=False,nonlinear_solution_constructed=False,physical_stability_proved=False,
        semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),
                    *ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),exception=error)),flush=True)
    if status!='passed':raise SystemExit(1)

if __name__=='__main__':main()
