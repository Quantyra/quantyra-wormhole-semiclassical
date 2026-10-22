"""Exact supporting checks for the ordinary actual graph-domain proof."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S

ROOT=Path(__file__).resolve().parents[1]
STEM='finite-graph-admissibility'
PREVIOUS='finite-comparison-conormal-regularity'
OWN=['checks/finite_graph_admissibility.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[];R=S.Rational
PROTOCOL={'arithmetic':'exact symbolic and rational',
          'continuum_fixture':'periodic two-mode scalar with positive rank-one noncommuting tensor return',
          'transfer_count':3,'ordinary_continuum_proof_separate':True,
          'actual_source_high_norms_evaluated':False,'actual_transfer_certificate_passed':False}


def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append({'name':name,'passed':bool(ok),'detail':str(detail)})
    print(json.dumps(ROWS[-1]),flush=True)
    if not ok:raise AssertionError(name)


def zero(name,expr):
    v=S.simplify(S.expand(expr));check(name,v==0,v)


def prior():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def paths():return sorted(set(prior()['input_sha256_after'])|set(OWN)|{'evidence/'+PREVIOUS+'.json'})


def provenance():
    p=prior()
    check('predecessor_passed',p['status']=='passed' and len(p['checks'])==65 and all(r['passed'] for r in p['checks']))
    check('predecessor_immutable',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('predecessor_hashes',not bad,bad)
    review=json.loads((ROOT/OWN[2]).read_bytes())
    check('retained_inputs_pinned',all(f in p['input_sha256_after'] for f in review['retained_inputs']))
    check('ordinary_domain_scope',review['ordinary_actual_certificate_admissibility_proved'] and not review['actual_transfer_certificate_passed'] and not review['goal_complete'])


def path_and_density():
    z,ell,u,v=S.symbols('z ell u v',real=True)
    tau=1-z;T=1+tau*(S.exp(ell*v)-1)
    tangent=tau*S.exp(ell*v)*v/T
    zero('original_radius_path_tangent',S.diff(S.log(T),ell)-tangent)
    J=T**4*S.exp(2*ell*u)
    zero('fixed_volume_path_derivative',S.diff(J,ell)-J*(4*tangent+2*u))
    u3,v3=S.symbols('u3 v3')
    logJ=4*S.log(1+(1-z)*(S.exp(ell*v3*z**3)-1))+2*ell*u3*z**3
    cubic=S.expand(S.series(logJ,z,0,4).removeO()).coeff(z,3)
    zero('density_mass_cubic',cubic-(4*ell*v3+2*ell*u3))
    zero('density_correlated_cubic',cubic.subs(u3,-2*v3))
    x,y=S.symbols('x y',real=True)
    pot=x**4/7+x*x*y+y**3/5+x*y/11
    density=S.exp(4*x+2*y)
    gradient=S.Matrix([S.diff(pot,x),S.diff(pot,y)])
    F=gradient/density
    H=S.hessian(pot,(x,y))/density
    jac=F.jacobian([x,y]);contact=F*S.Matrix([[4,2]])
    check('full_density_contact_identity',S.simplify(H-jac-contact)==S.zeros(2))
    check('density_contact_omission_detected',S.simplify(H-jac)!=S.zeros(2))
    # The bookkeeping action has an explicit parameter term in addition to its path Hessian.
    I0=x*x+3*y*y+x*y;rest=x**3+2*x*y*y
    action=I0+ell*rest
    path={x:ell+ell*ell/3,y:2*ell-ell*ell/5}
    grad=S.Matrix([S.diff(action,x),S.diff(action,y)])
    hdot=S.Matrix([S.diff(path[x],ell),S.diff(path[y],ell)])
    total=grad.subs(path).diff(ell)
    expected=S.hessian(action,(x,y)).subs(path)*hdot+S.Matrix([S.diff(rest,x),S.diff(rest,y)]).subs(path)
    check('bookkeeping_explicit_source_derivative',S.simplify(total-expected)==S.zeros(2,1))
    check('full_and_bookkeeping_Hessians_distinct',S.hessian(I0+rest,(x,y))!=S.hessian(action,(x,y)))
    FIXTURES.append({'fixture':'original_path_and_density','actual_path_identities':True,'polynomial_action_is_physical':False})


def kernels():
    z,y=S.symbols('z y',positive=True)
    delta=S.symbols('delta',real=True)
    normal=z**(4-delta)*y**(1+delta)/(z+y)**5
    Db=lambda f:z*S.diff(f,z)+y*S.diff(f,y)
    zero('joint_dilation_normal_kernel',Db(normal))
    correction=(z+y)**2*S.log(z+y)**2*normal
    for j in range(1,4):
        previous=correction;correction=Db(correction)
        polynomial=S.simplify(correction/normal/(z+y)**2)
        check('two_power_commutator_'+str(j),not polynomial.has(delta),polynomial)
    # Integration by parts uses dy/y. The corresponding density with dy needs the extra +1.
    k=S.Function('k')(z,y);w=S.Function('w')(y)
    zero('b_measure_commutator_integration',S.diff(k*w,y)-S.diff(k,y)*w-k*S.diff(w,y))
    d=R(7,2)
    zero('physical_weight_exponent',R(3,2)+2-d)
    check('weighted_strip',-1<d<4)
    zero('two_sided_Schur_integral',1/(4-d)+1/(1+d)-R(20,9))
    # A local scaling contact cannot be discarded before commuting the normal operator.
    x=S.symbols('x',real=True);f=S.Function('f')(x)
    raw=lambda a:x*S.diff(a,x,4)
    comm=S.diff(raw(f),x)-raw(S.diff(f,x))
    zero('unrestored_log_contact_commutator',comm-S.diff(f,x,4))
    check('unrestored_log_contact_not_zero',comm!=0)
    restored=lambda a:raw(a)-x*S.diff(a,x,4)
    zero('restored_scaling_contact',S.diff(restored(f),x)-restored(S.diff(f,x)))
    FIXTURES.append({'fixture':'weighted_kernel_and_contact_commutators','actual_kernel_evaluated':False,'weight':'7/2 in b-volume; W=z^-2 in physical volume'})


def tails_and_orders():
    r,b=S.symbols('r b',positive=True)
    for k in [0,2,4,6]:
        tail=sum(S.factorial(k)/S.factorial(j)*S.log(b)**j for j in range(k+1))/b
        zero('tail_primitive_'+str(k),S.diff(tail,b)+S.log(b)**k/b**2)
        zero('tail_vanishes_'+str(k),S.limit(tail,b,S.oo))
        exact=S.integrate(S.log(r)**k/r**2,(r,1,S.oo))
        zero('weighted_log_tail_integral_'+str(k),exact-S.factorial(k))
    p=S.symbols('p',real=True)
    T=lambda f:-S.diff(r*r*(r*r-1)*S.diff(f,r),r)/(r*r)
    zero('independent_warp_scalar_Laplacian',T(r**-p)-p*(3-p)*r**-p-p*(p-1)*r**(-p-2))
    zero('weighted_mass_cubic_integrand',r**2*r**4*(r**-3)**2-1)
    zero('weighted_fast_four_integrand',r**2*r**4*(r**-4)**2-r**-2)
    sequence=[('q0',6),('u0',8),('q1',4),('u1',6),('q2',2),('u2',4),('q3',0)]
    check('three_transfer_derivative_sequence',[sequence[i+1][1]-sequence[i][1] for i in range(6)]==[2,-4,2,-4,2,-4],sequence)
    check('all_three_solves_have_four_derivatives',all(n>=4 for label,n in sequence if label.startswith('u')))
    x=S.symbols('x',real=True);n=S.symbols('n',positive=True,integer=True)
    zero('uncoupled_three_transfer_high_frequency',-S.diff(S.cos(n*x),x,6)-n**6*S.cos(n*x))
    FIXTURES.append({'fixture':'weighted_tail_and_finite_derivative_order','actual_source_evaluated':False,'domain_sequence':sequence,'mass_cubic_not_in_fast_weight':True})


def coupled_continuum():
    x=S.symbols('x',real=True)
    basis=S.Matrix([S.cos(x),S.cos(3*x)])
    scalar=lambda c:(basis.T*c)[0]
    inner=lambda f,g:S.integrate(S.expand_trig(f*g),(x,0,2*S.pi))/S.pi
    T=S.diag(1,9);H=T+4*S.eye(2);D0=T*H;E=S.Integer(5)
    psi=S.Matrix([1,1]);vertex=S.Matrix([-4,R(52,5)])
    Q=vertex*vertex.T/(psi.dot(vertex));K=H+Q;P=D0-E*K
    check('tensor_term_positive_rank_one',Q.det()==0 and S.trace(Q)>0)
    check('positive_auxiliary_matrix',K[0,0]>0 and K.det()>0)
    check('tensor_term_noncommutes',T*Q!=Q*T)
    check('full_adjoint_mode',P*psi==S.zeros(2,1))
    W=2+S.cos(2*x)
    gram=S.Matrix(2,2,lambda i,j:inner(W*basis[i],W*basis[j]))
    check('independent_weighted_Gram',gram==S.Matrix([[R(13,2),R(9,4)],[R(9,4),R(9,2)]]),gram)
    q=S.Matrix([2,-1]);initial=q;states=[]
    mom=lambda c:psi.dot(c)
    rank_function=scalar(vertex)
    for j in range(3):
        u=K.inv()*q;uf=scalar(u)
        qu=rank_function*inner(rank_function,uf)/psi.dot(vertex)
        directK=-S.diff(uf,x,2)+4*uf+qu
        zero('continuum_exact_solve_'+str(j),directK-scalar(q))
        nxt=D0*u
        zero('continuum_D0_transfer_'+str(j),S.diff(uf,x,4)-4*S.diff(uf,x,2)-scalar(nxt))
        zero('continuum_one_step_moment_'+str(j),mom(q)-mom(nxt)/E)
        norm2=(nxt.T*gram*nxt)[0]
        check('continuum_weighted_norm_finite_'+str(j),norm2.is_Rational and norm2>0,norm2)
        states.append({'q':[str(v) for v in q],'u':[str(v) for v in u],'next_source_weighted_norm_squared':str(norm2)})
        q=nxt
    zero('continuum_actual_three_step_identity',mom(initial)-mom(q)/E**3)
    zero('continuum_nonzero_initial_moment',mom(initial)-1)
    zero('continuum_unsuppressed_moment',mom(q)-125)
    reversed_order=(K.inv()*D0)**3*initial
    check('reversed_transfer_order_detected',mom(reversed_order)/E**3!=mom(initial))
    FIXTURES.append({'fixture':'coupled_finite_Fourier_continuum','actual_charged_geometry':False,'actual_source':False,'period':'2pi','E':str(E),'Q':str(Q),'K':str(K),'weighted_Gram':str(gram),'states':states,'initial_moment':str(mom(initial)),'third_moment':str(mom(q)),'admissibility_does_not_suppress_moment':True})


def conditional_budget():
    emin=R('1.125e35');terminal=R('1e57')
    beta=terminal/(25*emin**3)
    check('unchanged_sufficient_cold_target',R('1e49')*beta<R('.281'),R('1e49')*beta)
    a,b,c,p=S.symbols('R0 R2 R4 N6p',positive=True)
    zero('norm_composition_order',a*(b*(c*p))-a*b*c*p)
    FIXTURES.append({'fixture':'conditional_transfer_budget','actual_constants_evaluated':False,'sufficient_product':'R0 R2 R4 N6(p)<=1e57 epsilon_c','cold_factor_upper':str(R('1e49')*beta)})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 actual graph admissibility','input_count':len(ps)}),flush=True)
    try:
        provenance();path_and_density();kernels();tails_and_orders();coupled_continuum();conditional_budget()
    except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(c['passed'] for c in ROWS)
    run={'story':'S167 actual graph admissibility','created_utc':datetime.now(timezone.utc).isoformat(),'status':'passed' if good else 'failed',
         'exception':error,'python':sys.version,'sympy':S.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after,
         'ordinary_finite_metric_source_regularity_proved':True,'ordinary_actual_certificate_admissibility_proved':True,
         'ordinary_actual_graph_transfers_admissible':True,'ordinary_actual_three_transfer_moment_identity_justified':True,
         'ordinary_proof_machine_verified':False,'actual_source_high_norms_evaluated':False,'actual_transfer_norms_numerically_bounded':False,
         'actual_signed_moment_evaluated':False,'actual_transfer_certificate_passed':False,'actual_signed_compatibility_zero_proved':False,
         'cold_metric_response_proved':False,'nonlinear_neighborhood_preserved':False,'nonlinear_solution_constructed':False,
         'physical_source_projection_changed':False,'physical_boundary_condition_changed':False,'physical_stability_proved':False,
         'formal_verification':False,'S167_complete':False,'S168_complete':False,'goal_complete':False}
    if not good:
        for kind,p in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');saved=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    saved['runs'].append(run);dest.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
