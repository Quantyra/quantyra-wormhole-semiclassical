"""Complete coincident finite-mode action; NOT a renormalized metric solver."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import sys,json,traceback
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_seed_separated_quantum_response as radial
prior=radial.prior
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]; STEM='current-coincident-dirac-action-assembly'
PREVIOUS='current-vacuum-contact-feasibility-closeout'
OWN=['notes/'+STEM+'.md','checks/current_coincident_dirac_action_assembly.py','evidence/'+STEM+'-review.json']
ROWS=[]; NUMERICS=[]; FIXTURES=[]; encoded=prior.encoded; hull=prior.hull
SCHEDULE=[256,512,1024]; CASES=[(1,'.25'),(1,'1'),(2,'1')]

def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok: raise AssertionError(name)
prior.check=check

def unit_interval(v): return hull(max(arb(0),v.lower()),min(arb(1),v.upper()))

def directions(z):
    if z.upper()<=arb('.5') or z.lower()>=arb('.875'): return [arb(0),arb(0)]
    w=(z-arb('.5'))/arb('.375'); w=unit_interval(w)
    u=unit_interval((4*w*(1-w))**6)
    return [u,unit_interval(w*u)]

def sensitivity_flow(a,f,r,L):
    E=(-a*L).exp()
    return E*r+(1-E)*f/a

def sensitivity_step(W,M,U,u,initial,L):
    damping=2*(M+W*U); forcing=2*M*u*U
    lo=sensitivity_flow(damping.upper(),forcing.lower(),initial[0],L).lower()
    hi=sensitivity_flow(damping.lower(),forcing.upper(),initial[1],L).upper()
    end=(max(arb(0),lo),min(arb(1),hi))
    whole=hull(min(initial[0],end[0]),max(initial[1],end[1]))
    return end,whole

def transfer_cell(W,M,initial,L):
    end=radial.step(W,M,initial,L)
    whole=hull(min(initial[0],end[0]),max(initial[1],end[1]))
    return end,whole

def integrate_matrix(M,cells,basis,initial_sensitivity=0,outer_independent=False):
    # cells hold the actual W and proper-length enclosures.
    left=(arb(0),arb(1)); sensitivities=[(arb(0),arb(initial_sensitivity)) for _ in range(2)]
    forward=[]
    for (W,L),us in zip(cells,basis):
        nxt,U=transfer_cell(W,M,left,L); rr=[]
        for j in range(2): sensitivities[j],whole=sensitivity_step(W,M,U,us[j],sensitivities[j],L); rr.append(whole)
        forward.append((U,rr)); left=nxt
    right=(arb(0),arb(1)) if outer_independent else left
    rs=[(arb(0),arb(1)) for _ in range(2)] if outer_independent else list(sensitivities)
    seam=dict(U=encoded(hull(*right)),sensitivity=[encoded(hull(*s)) for s in rs],encloses_paired_data_without_identifying_a_wall=outer_independent)
    matrix=[[arb(0) for _ in range(2)] for _ in range(2)]
    tadpole=[[arb(0) for _ in range(2)] for _ in range(2)]
    for k in range(len(cells)-1,-1,-1):
        W,L=cells[k]; us=basis[k]; nxt,V=transfer_cell(W,M,right,L); ss=[]
        for j in range(2): rs[j],whole=sensitivity_step(W,M,V,us[j],rs[j],L); ss.append(whole)
        U,rr=forward[k]; den=1+U*V; T=unit_interval((1-U*V)/den)
        for i in range(2):
            for j in range(2):
                local=-2*M*us[i]*us[j]*T*L
                matrix[i][j]+=local-4*M*us[i]*(V*rr[j]+U*ss[j])*L/den**2
                tadpole[i][j]+=local
        right=nxt
    return matrix,tadpole,seam

def continuous_matrix(par,n,p,N):
    nodes,cells=radial.mesh(par,N); M=(n*(n+par['q'])/par['C']).sqrt()
    basis=[directions(hull(right.lower(),left.upper())) for left,right in zip(nodes,nodes[1:])]
    matrix,tadpole,seam=integrate_matrix(M,[(p*iq,L) for iq,L in cells],basis)
    return matrix,tadpole,seam,(nodes,cells)

def decode(x):
    z=x['decimal_enclosure']; e=int(z['exponent'])
    return arb(str(z['mid'])+'e'+str(e),str(z['rad'])+'e'+str(e))

def cold_geometry(par,h,chi):
    beta,u,e=par['beta'],par['u'],par['eps']; v=1+h
    K=(v+1)*(u+3*beta*(v*v+1))/(2*v**3)-e*(v-3)/(2*v**5)
    y=2*(chi/2).sinh()**2; t=1+h*(1+y)
    H1=beta*(t+2*v)+(3*beta+u+e)/(t*v*v)-2*e*(2*t+v)/(t*t*v**3)+e*(3*t*t+2*t*v+v*v)/(t**3*v**4)
    B=2*K*v+y*H1; proper=((y+2)/(t*B)).sqrt()
    # W is singular at chi=0, so compute it only for a strictly positive cell.
    W=K*t*t.sqrt()/(2*(y*B).sqrt()) if chi.lower()>0 else None
    return proper,W

def cold_directions(chi):
    if chi.lower()>=2: return [arb(0),arb(0)]
    x=unit_interval(chi*chi/4); u=unit_interval((1-x)**6)
    return [u,unit_interval(x*u)]

def cold_execution(par):
    old=json.loads((ROOT/'evidence/current-bolt-adjoint-domain.json').read_bytes())['runs'][-1]
    base=next(v for v in old['numerics'] if v['kind']=='actual_current_bolt_normalization')
    h=decode(base['h']); limit=arb('1e-6')
    # A broad direct interval for cosh(chi) is too dependent for a tight
    # ratio, so verify the continuous global guard on a whole-cell cover.
    for i in range(64):
        chi=hull(arb(3*i)/64,arb(3*(i+1))/64); proper,_=cold_geometry(par,h,chi)
        check('factored_bolt_geometry_guard_'+str(i),proper>arb('.5') and proper<2)
    for n in [1,2]:
        M=(n*(n+par['q'])/par['C']).sqrt(); previous=None
        for N in [64,128,256]:
            nodes=[limit*((arb(3)/limit).log()*i/(20*N)).exp() for i in range(20*N+1)]
            nodes[0]=limit; nodes[-1]=arb(3); cells=[]; basis=[]
            for left,right in zip(nodes,nodes[1:]):
                chi=hull(left.lower(),right.upper()); proper,W=cold_geometry(par,h,chi)
                cells.append((W,(right-left)*proper)); basis.append(cold_directions(chi))
            matrix,tadpole,outer=integrate_matrix(M,cells,basis,initial_sensitivity=1,outer_independent=True)
            cap=hull((-20*M*limit).lower(),arb(0))
            matrix=[[v+cap for v in row] for row in matrix]; widths=[]
            for i in range(2):
                for j in range(2):
                    value=matrix[i][j]; widths.append((value.upper()-value.lower())/(-value.mid()))
                    check(f'cold_regulated_entry_{n}_{N}_{i}_{j}',value<0)
                    if previous is not None:
                        check(f'cold_overlapping_narrower_{n}_{N}_{i}_{j}',value.overlaps(previous[i][j]) and value.rad()<previous[i][j].rad())
            check(f'cold_coincident_reciprocity_{n}_{N}',matrix[0][1].overlaps(matrix[1][0]))
            if N==256: check(f'cold_final_width_target_{n}',all(w<arb('.15') for w in widths))
            NUMERICS.append(dict(kind='lowest_thermal_bolt_action_matrix',angular_n=n,thermal_j='0',omega='1/2',subdivisions=N,cells=len(cells),matrix=[[encoded(v) for v in row] for row in matrix],relative_widths=[encoded(w) for w in widths],omitted_tip_interval=encoded(cap),tip_is_included=True,outer_data=outer,complete_renormalized_operator=False))
            previous=matrix
        print(json.dumps(dict(case='regular_bolt',angular_n=n,thermal_j=0,status='finite_mode_assembled',complete_operator=False)),flush=True)

class Jet2:
    """Independent actual derivatives (not factorial-divided coefficients)."""
    def __init__(self,v,d=0,dd=0): self.v,self.d,self.dd=map(arb,[v,d,dd])
    @staticmethod
    def coerce(x): return x if isinstance(x,Jet2) else Jet2(x)
    def __add__(self,x):
        x=self.coerce(x); return Jet2(self.v+x.v,self.d+x.d,self.dd+x.dd)
    __radd__=__add__
    def __neg__(self): return Jet2(-self.v,-self.d,-self.dd)
    def __sub__(self,x): return self+-self.coerce(x)
    def __rsub__(self,x): return self.coerce(x)+-self
    def __mul__(self,x):
        x=self.coerce(x); return Jet2(self.v*x.v,self.d*x.v+self.v*x.d,self.dd*x.v+2*self.d*x.d+self.v*x.dd)
    __rmul__=__mul__
    def inv(self): return Jet2(1/self.v,-self.d/self.v**2,2*self.d**2/self.v**3-self.dd/self.v**2)
    def __truediv__(self,x): return self*self.coerce(x).inv()
    def __rtruediv__(self,x): return self.coerce(x)*self.inv()
    def exp(self):
        v=self.v.exp(); return Jet2(v,v*self.d,v*(self.dd+self.d**2))
    def log(self): return Jet2(self.v.log(),self.d/self.v,self.dd/self.v-(self.d/self.v)**2)
    def sqrt(self):
        v=self.v.sqrt(); return Jet2(v,self.d/(2*v),self.dd/(2*v)-self.d**2/(4*v**3))

def jet_step(W,m,U,L):
    K=(m*m+W*W).sqrt(); E=(-2*K*L).exp(); T=(1-E)/(1+E)
    growth=K*L+((1+E+(m+W*U)*(1-E)/K)/2).log()
    new=(W*T+(K-m*T)*U)/(K+m*T+W*T*U)
    return new,growth

def discrete_characteristic_matrix(par,n,p,mesh):
    nodes,cells=mesh; M=(n*(n+par['q'])/par['C']).sqrt().mid(); values=[]
    for direction in range(3):
        U=Jet2('.5'); amplitude=Jet2(0)
        for k,(iq,L) in enumerate(cells):
            z=((nodes[k]+nodes[k+1])/2).mid(); us=directions(z)
            u=(us[direction] if direction<2 else us[0]+us[1]).mid()
            m=Jet2(M,M*u,M*u*u)
            U,inc=jet_step((p*iq).mid(),m,U,L.mid()); amplitude+=inc
        action=-2*amplitude-(1+U*U).log()
        values.append(action.dd)
    cross=(values[2]-values[0]-values[1])/2
    return [[values[0],cross],[cross,values[1]]]

def identities():
    U,V,M,u,R,S,W=sp.symbols('U V M u R S W',positive=True)
    x,y=sp.symbols('x y',positive=True)
    rate=sp.simplify(((W*x-M*y)*x-y*(M*x+W*y))/x**2).subs(y,U*x)
    check('variable_mass_projective_equation',sp.simplify(rate-(W*(1-U**2)-2*M*U))==0)
    delta=sp.diff(W*(1-U**2)-2*M*U,U)*(-R)+sp.diff(W*(1-U**2)-2*M*U,M)*M*u
    check('complete_projective_sensitivity_equation',sp.simplify(-delta-(-(2*M+2*W*U)*R+2*M*u*U))==0)
    F=-M*(1-U*V)/(1+U*V)
    dF=sp.diff(F,M)*M*u-sp.diff(F,U)*R-sp.diff(F,V)*S
    full=-M*u*(1-U*V)/(1+U*V)-2*M*(V*R+U*S)/(1+U*V)**2
    check('coincident_tadpole_and_bubble_identity',sp.simplify(dF-full)==0)
    q=sp.symbols('q'); bump=(4*q*(1-q))**6
    check('compact_direction_C5',all(sp.diff(bump,q,j).subs(q,endpoint)==0 for j in range(6) for endpoint in [0,1]))
    # Zero W has no separated bubble, and its full Hessian is the tadpole.
    m=Jet2(18,18,18); U0=Jet2(0); U1,amp=jet_step(arb(0),m,U0,arb('.3'))
    action=-2*amp-(1+U1*U1).log()
    check('independent_zero_W_tadpole_fixture',action.dd.contains(arb('-10.8')) and U1.v.contains(0))
    # Compare independent automatic derivatives with exact SymPy derivatives.
    t=sp.symbols('t'); expression=sp.log(1+sp.exp(t))+sp.sqrt(3+sp.exp(2*t))
    jet=(1+Jet2(0,1).exp()).log()+(3+Jet2(0,2).exp()).sqrt()
    check('independent_derivative_arithmetic',sp.simplify(sp.diff(expression,t,2).subs(t,0))==sp.Rational(9,8) and jet.d.contains(1) and jet.dd.contains(arb(9)/8))
    t,v,b,u,e=sp.symbols('t v b u e',positive=True)
    ext=lambda t:(t-1)**2*(b+2*b/t+(3*b+u)/t**2+e/t**4)
    K=(v+1)*(u+3*b*(v*v+1))/(2*v**3)-e*(v-3)/(2*v**5)
    H1=b*(t+2*v)+(3*b+u+e)/(t*v*v)-2*e*(2*t+v)/(t*t*v**3)+e*(3*t*t+2*t*v+v*v)/(t**3*v**4)
    factored=(t-v)*(2*(v-1)*K*v+(t-v)*H1)/t
    check('exact_nonextremal_bolt_factorization',sp.factor(ext(t)-v*ext(v)/t-factored)==0)
    h,y,B=sp.symbols('h y B',positive=True)
    check('exact_cold_proper_chart_cancellation',sp.cancel(h*h*y*(y+2)/(t**4*(h*h*y*B/t**3))-(y+2)/(t*B))==0)
    FIXTURES.append(dict(fixture='exact_first_order_contacts_and_independent_derivatives',continuum_quantum_operator_certified=False))

def execution():
    par=prior.parameters()
    for index,(n,target) in enumerate(CASES):
        j=(arb(target)/par['gamma_cl']).floor().unique_fmpz(); assert j is not None
        j=int(j); p=par['gamma']*(arb(j)+arb(1)/2); previous=None
        for N in SCHEDULE:
            matrix,tadpole,seam,mesh=continuous_matrix(par,n,p,N)
            widths=[]
            for i in range(2):
                for k in range(2):
                    value=matrix[i][k]; width=(value.upper()-value.lower())/(-value.mid())
                    widths.append(width)
                    check(f'negative_regulated_entry_{index}_{N}_{i}_{k}',value<0)
                    if previous is not None:
                        old=previous[i][k]
                        check(f'overlapping_narrower_continuum_entry_{index}_{N}_{i}_{k}',value.overlaps(old) and value.rad()<old.rad())
            check(f'coincident_reciprocity_{index}_{N}',matrix[0][1].overlaps(matrix[1][0]))
            row=dict(kind='continuous_finite_mode_action_matrix',angular_n=n,thermal_j=str(j),p=encoded(p),subdivisions=N,cells=len(mesh[1]),matrix=[[encoded(v) for v in a] for a in matrix],tadpole_matrix=[[encoded(v) for v in a] for a in tadpole],relative_widths=[encoded(w) for w in widths],paired_seam=seam,normalization='one angular eigenstate, both thermal signs, integrated action on one physical half',angular_weight=str(2*(par['q']+2*n)),complete_renormalized_operator=False)
            if N==SCHEDULE[-1]:
                independent=discrete_characteristic_matrix(par,n,p,mesh)
                for i in range(2):
                    for k in range(2):
                        check(f'independent_characteristic_within_continuum_{index}_{i}_{k}',matrix[i][k].contains(independent[i][k]))
                check(f'final_matrix_width_target_{index}',all(w<arb('.15') for w in widths))
                row['independent_piecewise_characteristic_matrix']=[[encoded(v) for v in a] for a in independent]
            NUMERICS.append(row); previous=matrix
        print(json.dumps(dict(case=index,angular_n=n,p_target=target,status='finite_mode_assembled',complete_operator=False)),flush=True)
    FIXTURES.append(dict(fixture='continuous_coincident_matrices_and_independent_characteristic',angular_or_thermal_tail_certified=False))
    cold_execution(par)
    FIXTURES.append(dict(fixture='physical_lowest_thermal_bolt_matrix_with_tip_included',finite_temperature_and_regular_bolt_enclosed=True,complete_mode_sum=False))

def main():
    ctx.prec=256
    review=json.loads((ROOT/OWN[2]).read_bytes()); old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}; deps=prior.dependency_identity(); error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_budget_execution_passed',old['status']=='passed' and old['known_vacuum_contact_budget_substituted'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        identities(); execution()
    except Exception: error=traceback.format_exc()
    after={p:digest(p) for p in paths}; depafter=prior.dependency_identity(); status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),coincident_finite_mode_action_assembly_executed=status=='passed',paired_boundary_sensitivity_included=status=='passed',mass_exponential_tadpole_included=status=='passed',lowest_thermal_physical_bolt_matrix_enclosed=status=='passed',complete_renormalized_operator=False,full_thermal_angular_sum_evaluated=False,coincident_subtraction_executed=False,all_local_renormalized_moments_computed=False,full_bolt_adjoint_computed=False,full_defect_bolt_pairing_evaluated=False,actual_physical_correction_computed=False,differentiated_response_operator_bounded=False,raw_C5_admission_certified=False,nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json'); receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}; receipt['runs'].append(run)
    out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]: log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)),flush=True)
    if status!='passed': raise SystemExit(1)

if __name__=='__main__': main()
