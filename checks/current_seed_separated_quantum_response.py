"""Enclose one full separated static shape kernel, not the metric solution.

Run with --arb-path pointing to the isolated python-flint==0.9.0 install.
The preceding immutable evaluator supplies parameters and Arb serialization.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import factorial
import json, traceback, sys, time
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_seed_radial_green_enclosure as prior
from flint import arb, ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-seed-separated-quantum-response'
PREVIOUS='current-seed-radial-green-enclosure'
OWN=['notes/'+STEM+'.md','checks/current_seed_separated_quantum_response.py','evidence/'+STEM+'-review.json']
ROWS=[]; FIXTURES=[]; NUMERICS=[]
encoded=prior.encoded; hull=prior.hull; bounds=prior.bounds
geometry=prior.geometry; dependency_identity=prior.dependency_identity
SCHEDULE=[(1024,256),(2048,512)]

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {x['name'] for x in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)
prior.check=check

def decoded(x):
    z=x['decimal_enclosure']
    return arb(int(z['mid']),int(z['rad']))*arb(10)**int(z['exponent'])

def relwidth(x):return (x.upper()-x.lower())/x.lower()

def flow(W,M,U,L):
    K=(M*M+W*W).sqrt(); T=(K*L).tanh()
    return (W*T+(K-M*T)*U)/(K+M*T+W*T*U)

def step(W,M,initial,L):
    wl,wh=bounds(W); ml,mh=bounds(M); lo,hi=initial
    low=flow(wl,mh,lo,L).lower(); high=flow(wh,ml,hi,L).upper()
    return max(arb(0),low),min(arb(1),high)

def cell_integral(W,M,start,end,L):
    vals=hull(min(start[0],end[0]),max(start[1],end[1]))
    return (M+W*vals)*L

def mesh(par,N):
    log2=arb(2).log()
    nodes=[1-((-4+arb(i)/N)*log2).exp() for i in range(4*N+1)]
    for i in [0,N,3*N,4*N]:nodes[i]=1-arb(2)**(-4+i//N)
    cells=[]
    for left,right in zip(nodes,nodes[1:]):
        Q,_=geometry(hull(right.lower(),left.upper()),par)
        assert Q>0
        invQ=1/Q.sqrt(); L=(left-right)*invQ
        assert L>0
        cells.append((invQ,L))
    return nodes,cells

def evaluate(par,n,p,N,nodes,cells):
    M=(n*(n+par['q'])/par['C']).sqrt()
    A=N; B=3*N; seam=4*N
    left=[(arb(0),arb(1))]; Iab=arb(0); Ias=arb(0)
    for i,(iq,L) in enumerate(cells):
        W=p*iq; end=step(W,M,left[-1],L)
        if i>=A:
            inc=cell_integral(W,M,left[-1],end,L); Ias+=inc
            if i<B:Iab+=inc
        left.append(end)
    right=left[-1]; at={seam:right}; Ibs=arb(0)
    for i in range(seam-1,A-1,-1):
        iq,L=cells[i]; W=p*iq; end=step(W,M,right,L)
        if i>=B:Ibs+=cell_integral(W,M,right,end,L)
        right=end
        if i in [A,B]:at[i]=right
    UA,UB=hull(*left[A]),hull(*left[B]); VB=hull(*at[B]); den=(1+UB*VB)**2
    direct=-4*UA*VB*(-2*Iab).exp()/den
    image=-4*UA*UB*(-2*(Ias+Ibs)).exp()/den
    traces=[]
    for i in [A,B,seam]:
        U,V=hull(*left[i]),hull(*at[i])
        traces.append((1-U*V)/(M*(1+U*V)))
    # Cancel apparent endpoint ratios before evaluating the square trace.
    off=((1-UA)*(1+VB)+(1+UA)*(1-VB))*(-Iab).exp()/(2*M*(1+UB*VB))
    assert direct.is_finite() and image.is_finite()
    return direct,image,traces,off

def symbolic_and_constant_checks():
    U,V,Ub,M=sp.symbols('U V Ub M',positive=True)
    a=M*(1+U)/(1-U); b=M*(1-V)/(1+V)
    gm=(1-U)*(1+V)/(2*M*(1+Ub*V)); gp=gm*a*b/M**2
    SAB=sp.Matrix([[M*gp,-a*gm],[-M**2/a*gp,M*gm]])
    SBA=sp.Matrix([[M*gp,b*gm],[M**2/b*gp,M*gm]])
    check('first_order_matrix_cancellation',sp.factor(sp.trace(SAB*SBA)+4*U*V/(1+Ub*V)**2)==0)
    w=sp.symbols('w',positive=True)
    ua=(a-M)/(a+M)
    adot=M**2+2*w*a-a*a
    check('left_Riccati_change',sp.simplify(2*M*adot/(a+M)**2-(w*(1-U**2)-2*M*U))==0)
    bdot=M**2-2*w*b-b*b
    check('right_Riccati_change',sp.simplify(-2*M*bdot/(M+b)**2-(w*(1-V**2)-2*M*V))==0)
    check('logarithmic_amplitude_identity',sp.simplify(a-w-(M+w*U-(w*(1-U**2)-2*M*U)/(1-U)))==0)
    z=sp.symbols('z')
    check('geometry_tail_polynomial',sp.expand((1-z)**2*(1+2*z+3*z*z)-(1-4*z**3+3*z**4))==0)
    t,c=sp.symbols('t c',positive=True)
    for r in [1,2,3]:
        expression=2*sp.exp(-c*t)*sum(sp.Rational(factorial(2*r+1),factorial(2*r+1-k))*t**(2*r+1-k)/c**(k+1) for k in range(2*r+2))
        check('exact_tail_integral_power_'+str(r),sp.simplify(sp.diff(expression,t)+2*t**(2*r+1)*sp.exp(-c*t))==0)
    for idx,(ww,mm,uu,ll) in enumerate([('0','18','.7','.2'),('3','19','.1','.03'),('100','25','.8','.4')]):
        W,M0,U0,L=map(arb,[ww,mm,uu,ll]); K=(M0*M0+W*W).sqrt()
        cs=(K*L).cosh(); sn=(K*L).sinh()/K
        independent=(W*sn+(cs-M0*sn)*U0)/(cs+M0*sn+W*sn*U0)
        check('independent_constant_transfer_'+str(idx),flow(W,M0,U0,L).overlaps(independent))
        check('constant_semigroup_'+str(idx),flow(W,M0,flow(W,M0,U0,L/2),L/2).overlaps(independent))
        eq=W/(M0+K); stable=-4*eq*eq/(1+eq*eq)**2*(-2*(M0+W*eq)*L).exp()
        direct=-W*W/(M0*M0+W*W)*(-2*K*L).exp()
        check('constant_full_line_bubble_'+str(idx),stable.overlaps(direct))
    FIXTURES.append(dict(fixture='ordinary_identity_checks_and_independent_constant_transfer',continuum_proof_substituted=False))

def tail_series(N,c,q):
    a=18*c; t=arb(N).sqrt(); out=arb(0)
    for r,coef in [(1,arb(800)),(2,arb(2400)/q),(3,arb(1600)/q**2)]:
        assert a*t>2*r
        integ=2*(-a*t).exp()*sum((arb(factorial(2*r+1)//factorial(2*r+1-k))*t**(2*r+1-k)/a**(k+1) for k in range(2*r+2)),arb(0))
        out+=coef*(arb(N)**r*(-a*t).exp()+integ)
    return out

def tail_bounds(par,J):
    scale=arb('1.98')*arb(4).log()
    def constants(alpha,beta):return scale*arb(alpha),scale*arb(beta)*arb(5).sqrt()
    cM,cP=constants('.99','.14')
    angular=8*(1/cP+par['gamma'])*tail_series(3,cM,par['q'])
    cM,cP=constants('.7','.7')
    pmin=par['gamma']*(arb(J)+arb(1)/2)
    freq=8*(-cP*pmin).exp()*(1/cP+par['gamma'])*tail_series(1,cM,par['q'])
    return angular,freq,pmin

def parameter_guards(par):
    q=par['q']; b=par['beta']; u=par['u']; e=par['eps']
    check('global_Q_upper_point_two',b+u/16+e.abs_upper()<arb('.2'))
    check('R_increasing_exact_R_at_one',2*b-4*e.abs_upper()>0 and (6*b+u+e).contains(1))
    check('mass_and_multiplicity_tail_guards',q/par['C']<400 and q/par['C']>324)
    check('buffer_exponential_guard',(arb(18)/4).tanh()>arb('.99') and arb(2).log()>arb('.25'))
    check('angular_split_guard',arb('.99')**2+arb('.14')**2<1)
    check('thermal_split_guard',2*arb('.7')**2<1)

def previous_comparator(par,old):
    N=SCHEDULE[-1][0]; nodes,cells=mesh(par,N)
    oldrows=[x for x in old['numerics'] if x.get('cells')==81920]
    assert len(oldrows)==3
    for row in oldrows:
        n=row['angular_n']; j=int(row['thermal_j']); p=par['gamma']*(arb(j)+arb(1)/2)
        direct,im,ds,off=evaluate(par,n,p,N,nodes,cells)
        refs=[decoded(d['radial_trace']) for d in row['diagonal']]
        refs.append(decoded(row['separated_full_double_radial_trace']))
        name='previous_square_overlap_n'+str(n)+'_p'+str(row['p_target_label'])
        check(name,all(a.overlaps(b) for a,b in zip(ds+[off],refs)))
        NUMERICS.append(dict(kind='independent_previous_square_comparison',n=n,thermal_j=str(j),radial_traces=[encoded(d) for d in ds],separated_square_trace=encoded(off),direct_bubble=encoded(direct),image_bubble=encoded(im)))
    FIXTURES.append(dict(fixture='independent_prior_second_order_square_enclosures',sectors=3,physical_boundary_changed=False))

def full_sum(par,N,Nfreq):
    started=time.monotonic(); nodes,cells=mesh(par,N)
    indices=[]
    for l in range(Nfreq+1):
        pivot=3*Nfreq//4
        nominal=arb(4*l)/pivot if l<=pivot else 4+arb(12*(l-pivot))/(Nfreq-pivot)
        j=(nominal/par['gamma_cl']).floor().unique_fmpz()
        assert j is not None; indices.append(int(j))
    assert indices[0]==0 and all(b>a for a,b in zip(indices,indices[1:]))
    total_direct=arb(0); total_image=arb(0); bins=[]
    for n in [1,2]:
        weight=2*(par['q']+2*n)/arb(par['q'])*(n*(n+par['q'])/par['C'])
        for l,(jl,jr) in enumerate(zip(indices,indices[1:])):
            pl=par['gamma']*(arb(jl)+arb(1)/2); ph=par['gamma']*(arb(jr)-arb(1)/2)
            p=hull(pl.lower(),ph.upper()); multiplicity=par['gamma']*(jr-jl)
            dr,im,_,_=evaluate(par,n,p,N,nodes,cells)
            total_direct+=weight*multiplicity*dr; total_image+=weight*multiplicity*im
            bins.append(dict(n=n,thermal_first=str(jl),thermal_last=str(jr-1),count=str(jr-jl),actual_p=encoded(p),direct=encoded(dr),image=encoded(im)))
            if (l+1)%32==0:print(json.dumps(dict(event='thermal_bins_enclosed',N_oct=N,N_frequency=Nfreq,n=n,bins=l+1,elapsed_seconds=round(time.monotonic()-started,2))),flush=True)
    finite=total_direct+total_image; angular,freq,pmin=tail_bounds(par,indices[-1]); tail=angular+freq
    full=hull(finite.lower()-tail.upper(),finite.upper())
    QA,_=geometry(arb(7)/8,par)
    conversion=-2*par['q']/(arb.pi()*par['B2']*QA.sqrt())
    normalized=conversion*full
    result=dict(kind='full_separated_static_shape_sum',N_oct=N,N_frequency=Nfreq,actual_seed_parameter_box=True,preparation='.265',point_pair_z=['7/8','1/2'],explicit_angular_levels=[1,2],thermal_J_max=str(indices[-1]),frequency_tail_first_p=encoded(pmin),direct_finite_sum=encoded(total_direct),image_finite_sum=encoded(total_image),finite_sum=encoded(finite),angular_tail_absolute_upper=encoded(angular),frequency_tail_absolute_upper=encoded(freq),full_reduced_sum=encoded(full),relative_width=encoded(relwidth(-full)),tail_to_full_lower_magnitude=encoded(tail/(-full).lower()),normalized_disk_row_conversion=encoded(conversion),normalized_disk_row_kernel_per_drho=encoded(normalized),exact_thermal_bins=bins,elapsed_seconds=round(time.monotonic()-started,2))
    return result,full,tail

def calculation(old):
    par=prior.parameters(); parameter_guards(par)
    NUMERICS.append(dict(kind='parameters',values={k:encoded(v) if isinstance(v,arb) else v for k,v in par.items()}))
    previous_comparator(par,old)
    previous=None
    for N,F in SCHEDULE:
        result,full,tail=full_sum(par,N,F); NUMERICS.append(result)
        prefix='full_N'+str(N)+'_F'+str(F)
        check(prefix+'_finite_negative',full.is_finite() and full<0)
        if previous is not None:
            check(prefix+'_mesh_overlap',full.overlaps(previous))
            check(prefix+'_relative_width_refines',relwidth(-full)<relwidth(-previous))
            check(prefix+'_relative_width_target',relwidth(-full)<arb('.20'),encoded(relwidth(-full)))
            check(prefix+'_tail_target',tail/(-full).lower()<arb('.01'),encoded(tail/(-full).lower()))
        previous=full
    FIXTURES.append(dict(fixture='all_angular_and_discrete_thermal_separated_sum',explicit_n=[1,2],angular_tail='every n>=3',frequency_tail='all omitted half-integers, all n',coincident_subtraction_executed=False))

def main():
    ctx.prec=256
    review=json.loads((ROOT/OWN[2]).read_bytes()); old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}; error=None; deps=dependency_identity()
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_radial_enclosure_passed',old['status']=='passed' and old['selected_actual_radial_square_resolvents_enclosed'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        symbolic_and_constant_checks(); calculation(old)
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}; depafter=dependency_identity()
    status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),full_separated_thermal_angular_sum_enclosed=status=='passed',direct_and_paired_image_response_included=status=='passed',coincident_subtraction_executed=False,differentiated_response_operator_bounded=False,renormalized_quantum_source_evaluated=False,full_quantum_Hessian_evaluated=False,actual_physical_correction_computed=False,raw_C5_admission_certified=False,nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json'); receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run); out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)),flush=True)
    if status!='passed':raise SystemExit(1)

if __name__=='__main__':main()
