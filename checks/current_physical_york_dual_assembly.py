"""Frozen physical weak transpose and global trial functions; no metric solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import sys,json,traceback
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_all_angular_principal_operator as principal
endpoint=principal.endpoint;prior=principal.prior
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1];STEM='current-physical-york-dual-assembly'
PREVIOUS='current-all-angular-principal-operator'
OWN=['notes/'+STEM+'.md','checks/current_physical_york_dual_assembly.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];FIXTURES=[];encoded=prior.encoded
z=sp.symbols('z',positive=True);qs=sp.symbols('q0 q2 q3 q4 q5 q6');powers=[0,2,3,4,5,6]
Q=sum(c*z**k for c,k in zip(qs,powers));P=Q**2/z**4

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)
prior.check=check

def L(v):return -sp.diff(P*sp.diff(v,z),z)

def trial(m):
    """Global finite-seed TT trial, represented by its polynomial auxiliary J."""
    J=z**(m+4)
    xi=(m+4)*z**(m+3)*Q+(m+4)*z**4*sum(k*c*z**(k+m-1)/sp.Integer(k+m-1) for c,k in zip(qs,powers) if k>=2)
    xi=sp.expand(xi);b=sp.expand(Q*sp.diff(J,z)-xi)
    S=sp.expand(xi/2+Q*J/z);T=sp.expand(-xi-sp.diff(Q*J,z)/2)
    return dict(J=J,xi=xi,b=b,S=S,T=T)

def source(m,row):
    """Manufactured sources from the auxiliary transpose solution w."""
    w=z**(m+5)
    if row=='shape':
        f=sp.Integer(0)
        integral=sum(c*z**(k+m)/sp.Integer(k+m) for c,k in zip(qs,powers)) if m else 0
        g=-4*(m+5)*z**m*Q+2*m*(m+5)*integral
    else:
        f=sp.cancel(z*L(w)/Q);g=sp.Integer(0)
    returned=sp.expand(f/2-g+sp.diff(Q*w,z)/z**4)
    return dict(w=w,f=sp.expand(f),g=sp.expand(g),returned=returned)

def evaluate_polynomial(expr,env):
    poly=sp.Poly(sp.expand(expr),z,*qs);out=arb(0)
    for pp,c in poly.terms():
        value=arb(int(sp.numer(c)))/int(sp.denom(c))
        for var,k in zip([z,*qs],pp):
            if k:value*=env[var]**k
        out+=value
    return out

def integrate_polynomial(expr,env):
    anti=sp.integrate(sp.expand(expr),z)
    return evaluate_polynomial(anti,env)-evaluate_polynomial(anti,{**env,z:arb(0)})

def algebra():
    xi=sp.Function('xi')(z);J=sp.Function('J')(z);b=Q*sp.diff(J,z)-xi
    lhs=sp.diff(Q*b/z**4,z)+sp.diff(Q,z)*xi/z**4
    residual=L(J)-Q/z**4*(4*xi/z-sp.diff(xi,z))
    check('TT_constraint_auxiliary_operator',sp.simplify(lhs+residual)==0)
    h,N,B,G,gam=sp.symbols('h N B G gam',positive=True)
    W=-N*sp.pi*B**2/(G*gam);ell=(1+h)/h
    check('physical_action_area_factor',sp.simplify((-B**2/(2*G))*(2*sp.pi*N/gam)-W)==0)
    check('physical_adjoint_measure_invariant',sp.simplify(W*(ell/W)-ell)==0)
    check('two_bolt_sources_from_radius',sp.diff(sp.Symbol('s')+sp.Symbol('x')/2,sp.Symbol('x'))==sp.Rational(1,2))
    # The constant shape source must have zero auxiliary transpose source.
    check('constant_shape_source_cancellation',sp.simplify(Q*sp.diff(sp.Integer(1),z)/2)==0)
    check('special_zero_return_TT_direction',sp.simplify(4*z**4/z-sp.diff(z**4,z))==0)
    for m in range(3):
        tr=trial(m)
        check('global_trial_TT_'+str(m),sp.simplify(L(tr['J'])-Q/z**4*(4*tr['xi']/z-sp.diff(tr['xi'],z)))==0)
        check('global_trial_bolt_regular_identity_'+str(m),sp.simplify(tr['b']+tr['xi']-Q*sp.diff(tr['J'],z))==0)
        check('global_trial_phase_'+str(m),sp.Poly(tr['xi'],z).coeff_monomial(z**(m+3))==(m+4)*qs[0])
        for row in ['shape','scalar']:
            sr=source(m,row)
            check('manufactured_auxiliary_source_'+row+str(m),sp.simplify(L(sr['w'])-Q*(sr['f']/z+sp.diff(sr['g'],z)/2))==0)
    check('global_mass_cubic_retained',sp.Poly(trial(0)['xi'],z).coeff_monomial(z**3)==4*qs[0])
    check('shape_source_zero_fixture',sp.expand(source(0,'shape')['returned']-25*Q-z*sp.diff(Q,z))==0)
    check('scalar_source_zero_fixture',sp.expand(source(0,'scalar')['returned']-5*Q+4*z*sp.diff(Q,z))==0)

def actual_parameters():
    p=prior.parameters()
    rec=json.loads((ROOT/'evidence/current-bolt-adjoint-domain.json').read_bytes())['runs'][-1]
    old=next(x for x in rec['numerics'] if x['kind']=='actual_current_bolt_normalization')
    h=endpoint.raw.decode(old['h']);v=1+h;zh=1/v
    beta,u,e=p['beta'],p['u'],p['eps']
    tau=v*h*h*(beta*(1+2/v+3/v**2)+u/v**2+e/v**4)
    p={**p,'tau':tau};coeff=[beta,u,-4*beta-2*u-tau,3*beta+u+e,-2*e,e]
    env={z:zh,**dict(zip(qs,coeff))}
    def K(v):return (v+1)*(u+3*beta*(v*v+1))/(2*v**3)-e*(v-3)/(2*v**5)
    check('retained_temperature_root_relation',(h*K(v)).overlaps(p['gamma']))
    check('positive_exact_bolt_slope',2*p['gamma']>0)
    NUMERICS.append(dict(kind='actual_seed_correlated_endpoint',h=encoded(h),z_h=encoded(zh),gamma=encoded(p['gamma']),
        exact_Q_bolt='0',exact_Qprime_bolt='-2 gamma',Q_coefficients=[encoded(x) for x in coeff],
        polynomial_interval_subtraction_not_used_to_prove_bolt_slope=True))
    for N in [1,2]:
        W=-N*arb.pi()*p['B2']/p['gamma'];ell=(1+h)/h
        src=ell/W;stable=-(1+h)*K(v)/(N*arb.pi()*p['B2'])
        check('physical_normalized_adjoint_source_'+str(N),src.overlaps(stable) and src<0)
        check('physical_normalized_source_radius_'+str(N),(-arb.pi()*N*p['B2']*stable/(1+h)).overlaps(K(v)))
        NUMERICS.append(dict(kind='physical_action_and_adjoint_normalization',identical_exteriors=N,W=encoded(W),
            ell_b=encoded(ell),canonical_bolt_delta=encoded(stable),TT_bolt_delta=encoded(stable/2),
            normalization_improvement_claimed=False))
    return p,h,zh,env

def exterior_conversion(p,h,zh,cells):
    A=arb(0);Mf=arb(0);Mg=arb(0);records=[];guards=True
    for i in range(cells):
        left=arb(i)/(2*cells);right=arb(i+1)/(2*cells);width=right-left
        cell=prior.hull(left,right);q,qp=prior.geometry(cell,p)
        guards=guards and q>arb('.11') and q<arb('.17') and qp.abs_upper()<arb('.3')
        inc=(right**5-left**5)/(5*q*q);end=A+inc
        across=prior.hull(A.lower(),end.upper())
        if i==0:
            fm=arb(0,q.abs_upper()*right**5/(25*q.lower()**2))
        else:fm=width*q*across/cell
        gm=(width*abs(qp)*across+(right**5-left**5)/(5*q))/2
        Mf+=fm;Mg+=gm;A=end
        records.append(dict(z=[encoded(left),encoded(right)],Q=encoded(q),Qprime=encoded(qp),
            A_after=encoded(A),f_contribution=encoded(fm),g_envelope_contribution=encoded(gm)))
    factor=2*p['gamma']/(h*zh**4);cf=factor*Mf;cg=factor*Mg
    check('whole_exterior_geometry_'+str(cells),guards)
    check('exterior_f_conversion_'+str(cells),cf.abs_upper()<arb('.04'),encoded(cf))
    check('exterior_g_conversion_'+str(cells),cg.abs_upper()<arb('.08'),encoded(cg))
    NUMERICS.append(dict(kind='actual_exterior_source_conversion',cells=cells,z_cut='.5',
        f_integral=encoded(Mf),g_absolute_envelope_integral=encoded(Mg),factor=encoded(factor),
        f_coefficient=encoded(cf),g_coefficient_upper_enclosure=encoded(cg),cell_records=records,
        gravitational_inverse_or_metric_return=False))
    return cf,cg

def global_pairings(env,p,zh):
    for n in range(3):
        tr=trial(n)
        xi_h=evaluate_polynomial(tr['xi'],env)
        S_h=xi_h/2;T_h=-xi_h+p['gamma']*zh**(n+4)
        check('nonzero_bolt_trial_'+str(n),xi_h.abs_lower()>arb('.01'))
        NUMERICS.append(dict(kind='global_physical_TT_trial',index=n,J_power=n+4,
            fields={k:str(v) for k,v in tr.items()},xi_bolt=encoded(xi_h),S_bolt=encoded(S_h),T_bolt=encoded(T_h),
            mass_cubic=encoded(4*env[qs[0]]) if n==0 else '0',actual_correction=False))
        for m in range(3):
            for row in ['shape','scalar']:
                sr=source(m,row)
                direct=integrate_polynomial(sr['f']*tr['S']+sr['g']*tr['T'],env)
                returned=integrate_polynomial(sr['returned']*tr['xi'],env)
                label=row+str(m)+'_trial'+str(n)
                check('whole_domain_transpose_'+label,direct.overlaps(returned))
                check('whole_domain_pairing_width_'+label,max(direct.upper()-direct.lower(),returned.upper()-returned.lower())<arb('1e-16'))
                NUMERICS.append(dict(kind='global_manufactured_duality',source_row=row,source_index=m,trial_index=n,
                    direct_pairing=encoded(direct),transpose_pairing=encoded(returned),
                    actual_quantum_defect_or_residual=False))
    FIXTURES.append(dict(fixture='global_physical_TT_mass_and_fast_directions',count=3,
        full_trial_space_completeness_proved=False,actual_quantum_pairing=False))

def main():
    ctx.prec=256
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|
        {'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};deps=endpoint.dependencies();error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',
        protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('complete_angular_model_retained',old['status']=='passed' and old['all_angular_principal_model_enclosed'])
        check('previous_scientific_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        algebra();p,h,zh,env=actual_parameters()
        coarse=exterior_conversion(p,h,zh,128);fine=exterior_conversion(p,h,zh,256)
        check('exterior_grid_enclosures_overlap',all(a.overlaps(b) for a,b in zip(coarse,fine)))
        global_pairings(env,p,zh)
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths};depafter=endpoint.dependencies()
    status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,
        finished_utc=datetime.now(timezone.utc).isoformat(),frozen_physical_York_transpose_executed=status=='passed',
        exterior_kinematic_source_conversion_enclosed=status=='passed',global_physical_trial_directions_implemented=status=='passed',
        actual_physical_correction_computed=False,full_solution_feasibility_test_executed=False,
        complete_renormalized_operator=False,full_angular_curved_remainder_certified=False,
        full_adjoint_residual_norm_bounded=False,full_forward_residual_norm_bounded=False,
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
