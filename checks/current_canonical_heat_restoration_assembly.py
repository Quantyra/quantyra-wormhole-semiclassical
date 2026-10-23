"""Canonical local/LLL rows; not a massive geometric heat Hessian solver."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import sys,json,traceback
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_angular_bolt_cut_certificate as cap
local=cap.local;york=cap.york;endpoint=cap.endpoint;prior=cap.prior
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-canonical-heat-restoration-assembly';PREVIOUS='current-angular-bolt-cut-certificate'
OWN=['notes/'+STEM+'.md','checks/current_canonical_heat_restoration_assembly.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];FIXTURES=[];encoded=prior.encoded;Poly=local.Poly
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)
prior.check=check;york.check=check

def algebra():
    s,L,Z,Zp,F=sp.symbols('s L Z Zprime F')
    aa=sp.exp(2*s*L)*(Z+s*Zp)/(s-1)
    bb=sp.exp(2*s*L)*(2+s*F)
    check('a0_Mellin_finite_part',sp.simplify(sp.diff(aa,s).subs(s,0)+Zp+(1+2*L)*Z)==0)
    check('a2_Mellin_no_extra_two',sp.simplify(sp.diff(bb,s).subs(s,0)-F-4*L)==0)
    t,w,R,vi,vj,di,dj=sp.symbols('t w R vi vj di dj')
    u=t*vi+w*vj;lap=t*di+w*dj
    density=sp.exp(-2*u)*(R-2*lap)**2
    check('a2_gradient_density',sp.expand(sp.diff(density,t).subs({t:0,w:0})+2*R**2*vi+4*R*di)==0)
    expected=4*R**2*vi*vj+8*R*(vi*dj+vj*di)+8*di*dj
    check('a2_Hessian_density',sp.expand(sp.diff(density,t,w).subs({t:0,w:0})-expected)==0)
    N,B,C,G,ga,J=sp.symbols('N B C G gamma J',positive=True)
    W=-N*sp.pi*B**2/(G*ga);time=N*2*sp.pi/ga
    check('vacuum_action_normalization',sp.simplify((-J/(4*sp.pi*C))*time/W-G*J/(2*sp.pi*C*B**2))==0)
    check('curvature_action_normalization',sp.simplify((-C*F/(960*sp.pi))*time/W-G*C*F/(480*sp.pi*B**2))==0)
    check('heat_remainder_normalization',sp.simplify(time/(2*W)+G/B**2)==0)
    # Free noncommutative letters: this checks operator vertices, not domains.
    D,V,U=sp.symbols('D V U',commutative=False)
    # V and U commute as scalar multipliers; their products are symmetrized.
    ei=-(V*D+D*V)/2;ej=-(U*D+D*U)/2
    eij=((V*U+U*V)*D+D*(V*U+U*V))/8+(V*D*U+U*D*V)/4
    left=1-t*V/2-w*U/2+t*w*(V*U+U*V)/8
    transformed=sp.expand(left*D*left)
    check('fixed_Hilbert_D_first_vertex',sp.expand(transformed.coeff(t,1).coeff(w,0)-ei)==0)
    check('fixed_Hilbert_D_mixed_vertex',sp.expand(transformed.coeff(t,1).coeff(w,1)-eij)==0)
    square=sp.expand(-(D+t*ei+w*ej+t*w*eij)**2)
    check('geometric_square_includes_second_vertex',sp.expand(square.coeff(t,1).coeff(w,1)+eij*D+D*eij+ei*ej+ej*ei)==0)
    eta,E,P,D0,Db0=sp.symbols('eta E P D0 Db0',positive=True)
    error=eta*E*P+(Db0+eta*P)*(D0+eta*E)-Db0*D0
    check('operator_error_two_residual_polynomial',sp.expand(error-eta**2*E*P-eta*(E*P+Db0*E+D0*P))==0)
    a,b,M=sp.symbols('a b M',positive=True)
    threshold=2*M/(b+sp.sqrt(b*b+4*a*M))
    check('stable_positive_operator_threshold',sp.simplify(a*threshold**2+b*threshold-M)==0)

def retain_integral(kind,name,density,zh,**extra):
    value=density.integral(zh)
    NUMERICS.append(dict(kind=kind,name=name,value=encoded(value),density_coefficients=density.serialized(),
        full_domain_moments=density.moment_records(zh),**extra))
    return value

def execute():
    p,h,zh,env=york.actual_parameters();NUMERICS.extend(york.NUMERICS)
    Q=Poly({(power,0):env[sym] for power,sym in zip(york.powers,york.qs)})
    R=-Q.derivative().derivative();delta=lambda f:(Q*f.derivative()).derivative()
    C,B2=p['C'],p['B2'];q=arb(p['q']);k=arb(2000)
    C0=(1+12*arb.pi()*k**6).sqrt();b02=k**10*(C0-1)/6
    Lscale=(B2/b02).log()/2
    Hq=(q+1).digamma()+arb.const_euler();Feff=4*arb.const_euler()-2*Hq+4*Lscale
    U0=2*(604*q*q+(1+2*Lscale.abs_upper())*(q*q/6-arb(1)/30))/(arb.pi()*C*B2)
    c0=arb(0,U0.abs_upper());beta2=C*Feff/(480*arb.pi()*B2);lll=q/(6*arb.pi()*B2)
    capfactor=4*C/B2;capthreshold=B2/(4*C)
    check('known_vacuum_enclosure_sufficient',U0<arb('2.39121e-25'))
    check('restored_a2_negative_normalized_coefficient',Feff<0 and Feff>-100)
    check('actual_fixed_scale_bound',Lscale>arb('2.94524e-9') and Lscale<arb('3e-9'),encoded(Lscale))
    check('old_flat_fourth_margin_survives',-Feff/120>arb(1)/1000)
    x=p['g']*B2;outer_scale=x.log()/2-Lscale;kappa4=-Feff/4
    check('old_outer_fourth_margin_survives',x>arb('4e9') and outer_scale>9)
    check('old_interior_fourth_margin_survives',kappa4>arb('13.9'))
    check('coarse_rational_flat_margin_survives',arb(1)/6-4*arb('3e-9')>arb(3)/25)
    check('shared_vacuum_coefficient_not_zeroed',not c0.is_zero() and c0.contains(0))
    NUMERICS.append(dict(kind='common_restoration_coefficients',q=encoded(q),C=encoded(C),B2=encoded(B2),
        Lscale=encoded(Lscale),Feff=encoded(Feff),x=encoded(x),outer_scale_margin=encoded(outer_scale),
        kappa4=encoded(kappa4),flat_fourth_coefficient=encoded(-Feff/120),vacuum_coefficient_upper=encoded(U0),
        shared_signed_vacuum_coefficient=encoded(c0),beta2=encoded(beta2),LLL_Hessian_coefficient=encoded(lll),
        heat_cap_prefactor=encoded(capfactor),heat_cap_threshold_per_unit_error=encoded(capthreshold),
        heat_cap_threshold_for_1e_minus14_entry=encoded(capthreshold*arb('1e-14')),
        differentiated_heat_A_unknown=True,differentiated_heat_alpha_unknown=True,
        heat_cap_requirement='Aij*t0^alpha/alpha <= delta_H*B2/(4GC), alpha>0',
        entry_accuracy_is_physical_residual_target=False))
    trials=local.physical_trials(Q)[2:];matrix=[[None]*3 for _ in range(3)]
    vacuum_gram=[[None]*3 for _ in range(3)];deterministic=[[None]*3 for _ in range(3)]
    for i,trial in enumerate(trials):
        vi=trial['u'];di=delta(vi)
        check('physical_outer_value_and_first_trace_'+str(i),all(power>=3 and log==0 for power,log in vi.c))
        NUMERICS.append(dict(kind='actual_global_shape_trial',name=trial['name'],u=vi.serialized(),
            scalar=trial['nu'].serialized(),smooth_at_actual_bolt=True,actual_solution=False))
        raw=retain_integral('local_restoration_gradient','vacuum_shape_'+str(i),vi/2,zh,i=i,allocation='a0_factor_c0')
        a2=retain_integral('local_restoration_gradient','a2_'+str(i),beta2*(-2*R**2*vi-4*R*di),zh,i=i,allocation='a2')
        ll=retain_integral('local_restoration_gradient','LLL_'+str(i),lll/2*R*vi,zh,i=i,allocation='LLL')
        NUMERICS.append(dict(kind='bounded_restoration_gradient',i=i,known_part=encoded(a2+ll+c0*raw),
            shared_vacuum_factor=encoded(raw),complete_seed_source=False))
        for j in range(i,3):
            vj=trials[j]['u'];dj=delta(vj);tag=str(i)+'_'+str(j)
            gram=retain_integral('local_restoration_Hessian','Gram_'+tag,vi*vj,zh,i=i,j=j,allocation='a0_factor_c0')
            dens2=beta2*(4*R**2*vi*vj+8*R*(vi*dj+vj*di)+8*di*dj)
            hh2=retain_integral('local_restoration_Hessian','a2_'+tag,dens2,zh,i=i,j=j,allocation='a2')
            hll=retain_integral('local_restoration_Hessian','LLL_'+tag,lll*Q*vi.derivative()*vj.derivative(),zh,i=i,j=j,allocation='LLL')
            strong2=(8*beta2*vi*(delta(dj)+R*dj+delta(R*vj)+R**2*vj/2)).integral(zh)
            strongll=(-lll*vi*dj).integral(zh)
            check('a2_independent_strong_pairing_'+tag,hh2.overlaps(strong2))
            check('LLL_independent_strong_pairing_'+tag,hll.overlaps(strongll))
            NUMERICS.append(dict(kind='independent_strong_pairing',i=i,j=j,a2_weak=encoded(hh2),a2_strong=encoded(strong2),
                LLL_weak=encoded(hll),LLL_strong=encoded(strongll),outer_traces_zero=True,
                bolt_flux_zero_by_exact_Q_root=True,full_heat_domain_verified=False))
            pieces,flux,identity=local.action_bilinear(p,Q,trial,trials[j])
            old=retain_integral('retained_local_Hessian','old_'+tag,sum(pieces.values(),Poly(0)),zh,i=i,j=j)
            known=old+hh2+hll+c0*gram;det=old+hh2+hll
            check('known_matrix_entry_width_'+tag,known.upper()-known.lower()<arb('1e-14'))
            if i==j:check('vacuum_Gram_positive_'+str(i),gram>0)
            matrix[i][j]=matrix[j][i]=known;vacuum_gram[i][j]=vacuum_gram[j][i]=gram
            deterministic[i][j]=deterministic[j][i]=det
            NUMERICS.append(dict(kind='assembled_known_TT_entry',i=i,j=j,value=encoded(known),old_local=encoded(old),
                a2=encoded(hh2),LLL=encoded(hll),vacuum_factor=encoded(gram),vacuum_interval=encoded(c0*gram),
                massive_geometric_remainder_evaluated=False))
    NUMERICS.append(dict(kind='known_TT_matrix_with_shared_vacuum',basis=[v['name'] for v in trials],
        matrix=[[encoded(v) for v in row] for row in matrix],
        deterministic_matrix=[[encoded(v) for v in row] for row in deterministic],
        vacuum_Gram_matrix=[[encoded(v) for v in row] for row in vacuum_gram],
        shared_c0=encoded(c0),complete_TT_block=False,finite_basis_completeness_proved=False))
    FIXTURES.append(dict(kind='execution_boundary',canonical_local_and_LLL_assembled=True,
        actual_massive_heat_Hessian_executed=False,actual_forward_bolt_adjoint_residuals_executed=False,
        source_bound_differentiated_as_operator_bound=False,smooth_seam_assumed=False,
        physical_operator_threshold_formula_verified=True,actual_operator_threshold_numerically_evaluated=False,
        unknown_threshold_inputs=['E','P','D0','Db0','Q0','Aop','alpha']))

def main():
    ctx.prec=256;review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|
        {'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};deps=endpoint.dependencies();error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',
        protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_endpoint_certificate_passed',old['status']=='passed' and old['complete_angular_thermal_endpoint_error_bounded'])
        check('previous_scientific_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        algebra();execute()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths};depafter=endpoint.dependencies()
    status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,
        finished_utc=datetime.now(timezone.utc).isoformat(),canonical_local_and_LLL_assembled=status=='passed',
        shared_bounded_vacuum_retained=True,geometric_heat_remainder_identified=True,
        actual_physical_correction_computed=False,full_solution_feasibility_test_executed=False,
        full_adjoint_residual_norm_bounded=False,full_forward_residual_norm_bounded=False,
        complete_renormalized_operator=False,full_angular_curved_remainder_certified=False,
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
