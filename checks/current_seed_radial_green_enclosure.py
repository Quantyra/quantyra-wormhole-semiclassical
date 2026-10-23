"""Validated finite-seed radial square resolvents, not a metric solve.

Install python-flint==0.9.0 in an isolated directory and pass --arb-path.
All geometric/ODE enclosures use Arb; no point estimate is promoted to a ball.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import argparse, json, traceback, sys, platform, time
import sympy as sp

parser=argparse.ArgumentParser()
parser.add_argument('--arb-path')
args=parser.parse_args()
if args.arb_path:sys.path.insert(0,str(Path(args.arb_path).resolve()))
import flint
from flint import arb, ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-seed-radial-green-enclosure'
PREVIOUS='cold-defect-simultaneous-feasibility-audit'
OWN=['notes/'+STEM+'.md','checks/current_seed_radial_green_enclosure.py','evidence/'+STEM+'-review.json']
ROWS=[]; FIXTURES=[]; NUMERICS=[]

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {x['name'] for x in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def hull(lo,hi):return lo.union(hi)
def bounds(x):return x.lower(),x.upper()
def encoded(x):
    m,r,e=x.mid_rad_10exp(35)
    return dict(ball=x.str(30),decimal_enclosure=dict(mid=str(m),rad=str(r),exponent=str(e)))
def relative_width(x):return (x.upper()-x.lower())/x.lower()

def flow(v,m,L):
    k=v.sqrt(); t=(k*L).tanh()
    if m is None:return k/t
    return k*(m+k*t)/(k+m*t)

def step(V,initial,L):
    vlo,vhi=bounds(V); lo,hi=initial
    assert vlo>0 and L>0
    return flow(vlo,lo,L).lower(),flow(vhi,hi,L).upper()

def growth(V,initial,L):
    vlo,vhi=bounds(V); lo,hi=initial
    assert hi is not None
    def amplitude(v,m,d):
        k=v.sqrt(); arg=k*d
        return (arg.cosh()+(m/k)*arg.sinh()).log()
    # The cooperative (u,u') system gives monotonicity in all three inputs.
    return hull(amplitude(vlo,lo,L.lower()).lower(),amplitude(vhi,hi,L.upper()).upper())

def parameters():
    k=arb(2000); q=2000**4; pi=arb.pi(); zi=arb('0.265')
    lam0=k**-10; c=11/(360*pi)
    g=2*lam0/(1+(1+4*c*lam0).sqrt())
    Fl=1-384*pi*arb('1e22')*g
    Crn=(1+12*pi*g*k**16/Fl).sqrt(); Brn2=(Crn-1)/(6*g)
    eta=arb(0,'5.4e-23'); v=arb('-1.635e-29','6.5e-31')
    B2=Brn2/(1+eta); C=(Crn*(1+v)+eta)/(1+eta); d=Crn*v/(1+eta)
    u=1/C; eps=d/C; beta=(1-u-eps)/6; tau=arb('5e-66','5e-66')
    C0=(1+12*pi*k**6).sqrt(); b02=k**10*(C0-1)/6
    kap=2*zi*q*C0/(k*b02*b02.sqrt())
    gamma=kap*B2.sqrt()/C; gamma_cl=kap*b02.sqrt()/C0
    check('exact_seed_parameter_guards',g>0 and Fl>arb('0.99999998') and C>arb('4.8e10') and C<arb('5e10') and (q+1)/C>324)
    slope=12*beta*arb(4)/27+u/2+6*eps.abs_upper()+3*tau.abs_upper()
    check('global_positive_radial_potential_guard',slope<arb('.3'),encoded(slope))
    return dict(q=q,u=u,eps=eps,beta=beta,tau=tau,C=C,gamma=gamma,gamma_cl=gamma_cl,g=g,Fl=Fl,B2=B2,kappa=kap)

def geometry(z,par):
    b,u,e,t=par['beta'],par['u'],par['eps'],par['tau']
    omz=1-z
    R=b*(1+2*z+3*z*z)+u*z*z+e*z**4
    Q=omz*omz*R-t*z**3
    Qp=-12*b*z*z*omz+2*u*z*omz*(1-2*z)+2*e*z**3*omz*(2-3*z)-3*t*z*z
    return Q,Qp

def mesh(par,N):
    log2=arb(2).log()
    # The real nodes are defined by this exponential; endpoint balls enclose them.
    nodes=[1-((-10+arb(i)/N)*log2).exp() for i in range(10*N+1)]
    nodes[0]=1-arb(1)/1024; nodes[-1]=arb(0)
    for i in [7*N,9*N]:nodes[i]=1-arb(2)**(-10+i//N)
    cells=[]
    for left,right in zip(nodes,nodes[1:]):
        z=hull(right.lower(),left.upper()); Q,Qp=geometry(z,par)
        assert Q>0
        L=(left-right)/Q.sqrt()
        assert L>0
        cells.append((Q,Qp,L))
    return nodes,cells

def evaluate(par,n,ptarget,N,nodes,cells):
    j=(arb(ptarget)/par['gamma_cl']).floor().unique_fmpz()
    assert j is not None
    omega=arb(int(j))+arb(1)/2; p=omega*par['gamma']; M2=n*(n+par['q'])/par['C']
    assert p>arb('.99')
    potentials=[[],[]]
    for Q,Qp,L in cells:
        base=M2+p*p/Q; shift=p*Qp/(2*Q)
        potentials[0].append(base+shift); potentials[1].append(base-shift)
    obs=[7*N,9*N,10*N]
    left=[[(arb(0),None)],[(arb(0),None)]]
    loggrowth=[arb(0),arb(0)]
    for i,(_,_,L) in enumerate(cells):
        for chir in [0,1]:
            if obs[0]<=i<obs[1]:loggrowth[chir]+=growth(potentials[chir][i],left[chir][-1],L)
            left[chir].append(step(potentials[chir][i],left[chir][-1],L))
    inward=[left[1][-1],left[0][-1]]
    right_at={10*N:list(inward)}
    for i in range(10*N-1,-1,-1):
        L=cells[i][2]
        inward=[step(potentials[chir][i],inward[chir],L) for chir in [0,1]]
        if i in obs:right_at[i]=list(inward)
    diag=[]; raw=[]
    for i in obs:
        parts=[1/(hull(*left[chir][i])+hull(*right_at[i][chir])) for chir in [0,1]]
        trace=sum(parts,arb(0)); Q,_=geometry(nodes[i],par); a=Q.sqrt()/par['gamma']
        raw.append(parts)
        diag.append(dict(z=encoded(nodes[i]),radial_trace=encoded(trace),relative_width=encoded(relative_width(trace)),area_trace_one_thermal_number=encoded(trace/(2*arb.pi()*a))))
    seam=2/(hull(*left[0][-1])+hull(*left[1][-1]))
    assert seam.overlaps(sum(raw[-1],arb(0)))
    separated=sum((raw[1][ch]*(-loggrowth[ch]).exp() for ch in [0,1]),arb(0))
    result=dict(angular_n=n,p_target_label=ptarget,thermal_j=str(j),omega_exact_twice=str(2*int(j)+1),p_actual=encoded(p),M2_actual=encoded(M2),cells=10*N,diagonal=diag,separated_pair_z=['7/8','1/2'],separated_full_double_radial_trace=encoded(separated),separated_relative_width=encoded(relative_width(separated)),outer_seam_identity_enclosed=True,bolt_impedance_initial_enclosure='[0,+infinity], not a substituted boundary condition',full_quantum_sum_evaluated=False)
    return result,[sum(x,arb(0)) for x in raw],separated

def algebra():
    z,b,u,e,t=sp.symbols('z b u e t'); R=b*(1+2*z+3*z*z)+u*z*z+e*z**4
    Q=(1-z)**2*R-t*z**3
    derivative=-12*b*z*z*(1-z)+2*u*z*(1-z)*(1-2*z)+2*e*z**3*(1-z)*(2-3*z)-3*t*z*z
    check('factored_geometry_derivative',sp.expand(sp.diff(Q,z)-derivative)==0)
    check('geodesic_outer_boundary',sp.simplify(sp.diff(Q,z).subs(z,0))==0)
    f=sp.Function('Q')(z); p=sp.symbols('p')
    check('proper_radial_mass_derivative',sp.simplify(-sp.sqrt(f)*sp.diff(p/sp.sqrt(f),z)-p*sp.diff(f,z)/(2*f))==0)
    for v,m,L in [('3','.4','.7'),('325','12','.01'),('8000','100','1/50')]:
        vv,mm,dd=arb(v),arb(m),arb(L); kk=vv.sqrt()
        # Independent exact 2x2 constant-potential transfer matrix.
        c=(kk*dd).cosh(); s=(kk*dd).sinh()
        direct=(kk*s+mm*c)/(c+mm*s/kk)
        check('constant_transfer_'+v,flow(vv,mm,dd).overlaps(direct))
        split=flow(vv,flow(vv,mm,dd/2),dd/2)
        check('constant_transfer_semigroup_'+v,split.overlaps(direct))
    FIXTURES.append(dict(fixture='exact_geometry_and_independent_constant_transfer',boundary_replacement_used=False))

def calculation():
    par=parameters()
    NUMERICS.append(dict(parameters={k:encoded(v) if isinstance(v,arb) else v for k,v in par.items()},preparation='z_i=0.265',quantum_product_enclosure_used=True,pressure_quadrature_used=False))
    previous={}
    for N in [4096,8192]:
        started=time.monotonic(); nodes,cells=mesh(par,N)
        for n,p in [(1,1),(1,8),(2,1)]:
            row,diags,off=evaluate(par,n,p,N,nodes,cells); NUMERICS.append(row)
            prefix='n'+str(n)+'_p'+str(p)+'_N'+str(N)
            check(prefix+'_positive_finite',all(v>0 and v.is_finite() for v in diags+[off]))
            if N==8192:
                coarse_d,coarse_o=previous[(n,p)]
                check(prefix+'_mesh_overlap',all(a.overlaps(b) for a,b in zip(coarse_d+[coarse_o],diags+[off])))
                check(prefix+'_width_refines',all(relative_width(b)<relative_width(a) for a,b in zip(coarse_d+[coarse_o],diags+[off])))
                check(prefix+'_diagonal_precision_target',all(relative_width(v)<arb('.005') for v in diags),[encoded(relative_width(v)) for v in diags])
                check(prefix+'_separated_precision_target',relative_width(off)<arb('.05'),encoded(relative_width(off)))
            previous[(n,p)]=(diags,off)
            print(json.dumps(dict(event='sector_enclosed',angular_n=n,p_target=p,N_oct=N,elapsed_seconds=round(time.monotonic()-started,2))),flush=True)
    FIXTURES.append(dict(fixture='actual_finite_seed_three_mode_enclosures',number_of_sectors=3,preparation_values=['.265'],parameter_errors_enclosed=True,unknown_bolt_data_enclosed=True,renormalized_source_evaluated=False))

def dependency_identity():
    root=Path(flint.__file__).resolve().parent.parent
    payload={str(p.relative_to(root)):sha256(p.read_bytes()).hexdigest() for p in root.rglob('*') if p.is_file() and p.suffix in ['.pyd','.dll','.py'] and ('flint' in str(p.relative_to(root)))}
    return dict(python=sys.version,platform=platform.platform(),python_flint=flint.__version__,binary_and_source_sha256=payload,precision_bits=ctx.prec,dependency_directory=str(root))

def main():
    ctx.prec=256
    review=json.loads((ROOT/OWN[2]).read_bytes()); old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}; error=None; deps=dependency_identity()
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_audit_passed',old['status']=='passed' and old['value_envelope_to_raw_C5_rejected'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',flint.__version__=='0.9.0')
        algebra(); calculation()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    depafter=dependency_identity()
    run.update(status='passed' if error is None and before==after and deps==depafter else 'failed',exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),selected_actual_radial_square_resolvents_enclosed=error is None,bolt_impedance_enclosed_without_boundary_substitution=error is None,full_thermal_angular_sum_evaluated=False,renormalized_quantum_source_evaluated=False,full_quantum_Hessian_evaluated=False,actual_physical_correction_computed=False,raw_C5_admission_certified=False,nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json'); receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run); out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)),flush=True)
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
