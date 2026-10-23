"""Actual finite-seed local weak action; canonical shape determinant remains."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import sys,json,traceback
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_bolt_contact_parametrix as bolt
york=bolt.york;endpoint=bolt.endpoint;prior=bolt.prior
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-finite-local-action-assembly';PREVIOUS='current-bolt-contact-parametrix'
OWN=['notes/'+STEM+'.md','checks/current_finite_local_action_assembly.py','evidence/'+STEM+'-review.json']
ROWS=[];NUMERICS=[];FIXTURES=[];encoded=prior.encoded
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)
prior.check=check;york.check=check

class Poly:
    """Finite Laurent/log polynomial with interval coefficients."""
    def __init__(self, terms=0):
        if isinstance(terms,Poly):self.c=dict(terms.c)
        elif isinstance(terms,dict):self.c={k:arb(v) for k,v in terms.items() if not arb(v).is_zero()}
        else:self.c={} if arb(terms).is_zero() else {(0,0):arb(terms)}
    def __add__(self,other):
        if isinstance(other,Jet):return NotImplemented
        other=Poly(other);d=dict(self.c)
        for k,v in other.c.items():d[k]=d.get(k,arb(0))+v
        return Poly(d)
    __radd__=__add__
    def __neg__(self):return Poly({k:-v for k,v in self.c.items()})
    def __sub__(self,other):
        if isinstance(other,Jet):return NotImplemented
        return self+-Poly(other)
    def __rsub__(self,other):return Poly(other)+-self
    def __mul__(self,other):
        if isinstance(other,Jet):return NotImplemented
        other=Poly(other);d={}
        for (i,j),a in self.c.items():
            for (k,l),b in other.c.items():
                key=(i+k,j+l);d[key]=d.get(key,arb(0))+a*b
        return Poly(d)
    __rmul__=__mul__
    def __truediv__(self,other):return self*(1/arb(other))
    def __pow__(self,n):
        assert isinstance(n,int) and n>=0
        out=Poly(1)
        for _ in range(n):out=out*self
        return out
    def derivative(self):
        d={}
        for (p,l),v in self.c.items():
            if p:d[(p-1,l)]=d.get((p-1,l),arb(0))+p*v
            if l:d[(p-1,l-1)]=d.get((p-1,l-1),arb(0))+l*v
        return Poly(d)
    def integral(self,zh):
        out=arb(0)
        for (p,l),v in self.c.items():
            assert p>-1 and l in [0,1],(p,l,encoded(v))
            moment=zh**(p+1)/(p+1)
            if l:moment*=zh.log()-arb(1)/(p+1)
            out+=v*moment
        return out
    def moment_records(self,zh):
        out=[]
        for (p,l),v in sorted(self.c.items()):
            assert p>-1 and l in [0,1]
            moment=zh**(p+1)/(p+1)
            if l:moment*=zh.log()-arb(1)/(p+1)
            out.append(dict(power=p,log_power=l,coefficient=encoded(v),moment=encoded(moment),contribution=encoded(v*moment)))
        return out
    def value(self,z):return sum((v*z**p*z.log()**l for (p,l),v in self.c.items()),arb(0))
    def serialized(self):return [dict(power=p,log_power=l,coefficient=encoded(v)) for (p,l),v in sorted(self.c.items())]
    def zero_enclosed(self):return all(v.contains(0) for v in self.c.values())

def monomial(p,l=0):return Poly({(p,l):arb(1)})

class Jet:
    """Value, two first derivatives and their mixed second derivative."""
    def __init__(self,v=0,i=0,j=0,ij=0):self.v,self.i,self.j,self.ij=map(Poly,[v,i,j,ij])
    @staticmethod
    def lift(x):return x if isinstance(x,Jet) else Jet(x)
    def __add__(self,o):
        o=Jet.lift(o);return Jet(self.v+o.v,self.i+o.i,self.j+o.j,self.ij+o.ij)
    __radd__=__add__
    def __neg__(self):return Jet(-self.v,-self.i,-self.j,-self.ij)
    def __sub__(self,o):return self+-Jet.lift(o)
    def __rsub__(self,o):return Jet.lift(o)+-self
    def __mul__(self,o):
        o=Jet.lift(o)
        return Jet(self.v*o.v,self.i*o.v+self.v*o.i,self.j*o.v+self.v*o.j,
            self.ij*o.v+self.i*o.j+self.j*o.i+self.v*o.ij)
    __rmul__=__mul__
    def __truediv__(self,o):return self*(1/arb(o))
    def __pow__(self,n):
        out=Jet(1)
        for _ in range(n):out=out*self
        return out
    def derivative(self):return Jet(*(x.derivative() for x in [self.v,self.i,self.j,self.ij]))

def exp_direction(i,j,k):return Jet(1,k*i,k*j,k*k*i*j)

def action_bilinear(p,Q,di,dj):
    """Weak density for any returned direction jets representable here."""
    C,B2=p['C'],p['B2'];q=p['q'];lam0=arb(2000)**-10
    ni,ui=di['nu'],di['u'];nj,uj=dj['nu'],dj['u']
    phi=Jet(-monomial(0,1),ni,nj);u=Jet(0,ui,uj)
    delta=lambda f:(f.derivative()*Q).derivative()
    L=-Q.derivative().derivative()-2*delta(u)
    M=delta(phi);V=phi.derivative()**2*Q;D=L-6*M-6*V
    en2=exp_direction(ni,nj,2);en4=exp_direction(ni,nj,4)
    eu2=exp_direction(ui,uj,2);eum2=exp_direction(ui,uj,-2)
    LE=en2*D*monomial(-2)/2+en2*eu2*monomial(-2)/C+3*lam0*B2/C*en4*eu2*monomial(-4)
    LEM=-arb.pi()*arb(q)**2*arb(2000)**8/(C*B2)*eu2
    LR=8*arb.pi()*arb('1e22')/(C*B2)*(4*eu2+4*C*D+C*C*eum2*D**2)
    w_bracket=phi*L**2/60+L*M/30-M**2/10+L*V/30-7*M*V/90-7*V**2/180
    LW=(C*C*eum2*w_bracket-C*phi*L/18+(arb(1)/15-arb(q)**2/3)*phi*eu2-C*V/18)/(2*arb.pi()*C*B2)
    flux=(phi*(eum2*L).derivative()-phi.derivative()*(eum2*L))*Q*C/(60*arb.pi()*B2)
    raw_difference=(phi*delta(eum2*L)-(eum2*L)*delta(phi))*C/(60*arb.pi()*B2)
    return dict(Einstein=LE.ij,magnetic=LEM.ij,R2=LR.ij,Weyl=LW.ij),flux.ij,raw_difference.ij-flux.ij.derivative()

def scalar_covariant(p,Q,si,sj):
    C,B2=p['C'],p['B2'];lam0=arb(2000)**-10;x=p['g']*B2;d=p['eps']*C
    R=-12*x+d*(4*monomial(5)-6*monomial(6))
    tau=(1+3*x+d)*monomial(4)-d*(5*monomial(5)-arb(9)/2*monomial(6))
    delta=lambda f:(Q*f.derivative()).derivative()
    Ai=delta(si)-2*Q*si.derivative()*monomial(-1)
    Aj=delta(sj)-2*Q*sj.derivative()*monomial(-1)
    grad=Q*si.derivative()*sj.derivative()*monomial(-2)
    mass=(2*R+48*lam0*B2)*si*sj*monomial(-4)/C
    eps=1/(60*arb.pi());DR=96*arb.pi()*arb('1e22')
    return dict(Einstein=6*grad+mass,magnetic=Poly(0),
        R2=6*DR*C/B2*Ai*Aj-2*DR/B2*R*grad,
        Weyl=-6*eps*C/B2*Ai*Aj+(eps/(6*B2)*R-22*eps/(3*B2)*tau)*grad)

def trace_row_variation(p,Q,direction):
    """Independent Euler trace row, differentiated once in a returned field."""
    C,B2=p['C'],p['B2'];lam0=arb(2000)**-10
    nu,u=direction['nu'],direction['u']
    phi=Jet(-monomial(0,1),nu);uj=Jet(0,u)
    delta=lambda f:(f.derivative()*Q).derivative()
    L=Jet(-Q.derivative().derivative())-2*delta(uj)
    M=delta(phi);V=phi.derivative()**2*Q;D=L-6*M-6*V
    en2=exp_direction(nu,Poly(0),2);en4=exp_direction(nu,Poly(0),4)
    eu2=exp_direction(u,Poly(0),2);eum2=exp_direction(u,Poly(0),-2)
    U=2+C*eum2*D
    FE=en2*D*monomial(-2)+2*en2*eu2*monomial(-2)/C+12*lam0*B2/C*en4*eu2*monomial(-4)
    FR=-96*arb.pi()*arb('1e22')/B2*(delta(U)-2*(U*phi.derivative()*Q).derivative())
    FW=(C*C*(eum2*L**2/60+delta(eum2*L)/30
        +delta(eum2*(-M/5-7*V/90))
        -2*(eum2*(L/30-7*M/90-7*V/90)*phi.derivative()*Q).derivative())
        -C*L/18+(arb(1)/15-arb(p['q'])**2/3)*eu2+C*M/9)/(2*arb.pi()*C*B2)
    return dict(Einstein=FE.i,R2=FR.i,Weyl=FW.i,magnetic=Poly(0))


def physical_trials(Q):
    out=[dict(name='scalar_z4',nu=monomial(4),u=Poly(0),scalar=True),
         dict(name='scalar_z5',nu=monomial(5),u=Poly(0),scalar=True)]
    for m in [0,1]:
        J=monomial(m+4)
        extra=sum((k*c*monomial(k+m-1)/(k+m-1) for (k,l),c in Q.c.items() if k>=2),Poly(0))
        xi=(m+4)*(monomial(m+3)*Q+monomial(4)*extra)
        S=xi/2+Q*J*monomial(-1);T=-xi-(Q*J).derivative()/2
        out.append(dict(name='TT_J'+str(m+4),nu=S,u=T,scalar=False,xi=xi,J=J))
    xi=monomial(4)
    out.append(dict(name='TT_zero_J',nu=xi/2,u=-xi,scalar=False,xi=xi,J=Poly(0)))
    return out

def algebra():
    uncertain=arb(0,'1e-20')
    check('zero_containing_coefficients_retained',len(Poly({(5,0):uncertain}).c)==1 and len(Poly(uncertain).c)==1)
    check('only_exact_zero_coefficients_removed',not Poly({(5,0):arb(0)}).c)
    s,R,la,G=sp.symbols('s R lambda G',positive=True)
    grav=-(sp.exp(2*s)*R+6*la*sp.exp(4*s))/(16*sp.pi*G)
    check('Euclidean_trace_gradient_sign',sp.simplify(sp.diff(grav,s).subs(s,0)+(R+12*la)/(8*sp.pi*G))==0)
    check('Euclidean_scalar_action_second_sign',sp.simplify(sp.diff(grav,s,2).subs(s,0)+(2*R+48*la)/(8*sp.pi*G))==0)
    t,v,ui,uj,D0,di,dj,dij,C=sp.symbols('t v ui uj D0 di dj dij C')
    dd=D0+t*di+v*dj+t*v*dij;uu=t*ui+v*uj
    sq=4*sp.exp(2*uu)+4*C*dd+C*C*sp.exp(-2*uu)*dd*dd
    expected=16*ui*uj+4*C*dij+C*C*(4*ui*uj*D0*D0-4*D0*(ui*dj+uj*di)+2*di*dj+2*D0*dij)
    check('unexpanded_R2_second_variation',sp.expand(sp.diff(sq,t,v).subs({t:0,v:0})-expected)==0)
    check('magnetic_conformal_invariance',sp.diff(sp.exp(4*s)*sp.exp(-4*s),s)==0)
    N,B,ga=sp.symbols('N B gamma',positive=True)
    W=-N*sp.pi*B**2/(G*ga)
    measure=N*8*sp.pi**2*B**4/(C*ga)
    check('normalized_trace_row_is_positive',sp.simplify((-measure/(8*sp.pi*G))/W-B**2/C)==0)

def execute():
    p,h,zh,env=york.actual_parameters();NUMERICS.extend(york.NUMERICS)
    coeff={power:env[sym] for power,sym in zip(york.powers,york.qs)}
    Q=Poly({(power,0):v for power,v in coeff.items()})
    D0=-Q.derivative().derivative()+6*Q.derivative()*monomial(-1)-12*Q*monomial(-2)
    x=p['g']*p['B2'];d=p['eps']*p['C']
    direct=monomial(2)*(2+p['C']*D0)
    expected=-12*x+d*(4*monomial(5)-6*monomial(6))
    check('finite_seed_scalar_curvature_identity',(direct-expected).zero_enclosed(),(direct-expected).serialized())
    trials=physical_trials(Q);matrix=[[None]*5 for _ in range(5)]
    for i,di in enumerate(trials):
        NUMERICS.append(dict(kind='actual_global_trial',name=di['name'],nu=di['nu'].serialized(),u=di['u'].serialized(),
            scalar=di['scalar'],actual_solution=False))
        for j in range(i,5):
            dj=trials[j];parts,flux,identity=action_bilinear(p,Q,di,dj)
            tag=str(i)+'_'+str(j)
            check('raw_Weyl_flux_identity_'+tag,identity.zero_enclosed())
            check('outer_second_variation_flux_vanishes_'+tag,all(k>0 for k,l in flux.c))
            total=sum(parts.values(),Poly(0));integrals={k:v.integral(zh) for k,v in parts.items()}
            combined=total.integral(zh)
            check('local_matrix_width_'+tag,combined.upper()-combined.lower()<arb('1e-14'),encoded(combined))
            check('local_allocations_match_'+tag,combined.overlaps(sum(integrals.values(),arb(0))))
            if di['scalar'] and dj['scalar']:
                cov=scalar_covariant(p,Q,di['nu'],dj['nu'])
                for key,val in integrals.items():
                    independent=cov[key].integral(zh)
                    check('independent_scalar_'+key+'_'+tag,val.overlaps(independent),str((encoded(val),encoded(independent))))
                    NUMERICS.append(dict(kind='independent_complete_scalar_pairing',i=i,j=j,allocation=key,
                        action_value=encoded(val),covariant_trace_value=encoded(independent)))
            if di['scalar'] and not dj['scalar']:
                trace=trace_row_variation(p,Q,dj)
                for key,val in integrals.items():
                    independent=(trace[key]*di['nu']).integral(zh)
                    check('independent_mixed_trace_'+key+'_'+tag,val.overlaps(independent),str((encoded(val),encoded(independent))))
                    NUMERICS.append(dict(kind='independent_complete_mixed_pairing',i=i,j=j,allocation=key,
                        action_value=encoded(val),euler_trace_value=encoded(independent)))
            matrix[i][j]=combined;matrix[j][i]=combined
            NUMERICS.append(dict(kind='actual_finite_local_hessian_pairing',i=i,j=j,
                normalized_value=encoded(combined),allocations={k:encoded(v) for k,v in integrals.items()},
                density_coefficients=total.serialized(),boundary_flux_coefficients=flux.serialized(),
                full_domain_moments=total.moment_records(zh),
                complete_physical_pairing=di['scalar'] or dj['scalar'],canonical_shape_hessian_added=False))
    check('complete_scalar_block_positive',matrix[0][0]>0 and matrix[1][1]>0 and matrix[0][0]*matrix[1][1]-matrix[0][1]**2>0)
    W=-arb.pi()*p['B2']/p['gamma']
    NUMERICS.append(dict(kind='signed_actual_local_matrix',basis=[v['name'] for v in trials],
        matrix=[[encoded(v) for v in row] for row in matrix],physical_action_multiplier=encoded(W),
        complete_scalar_and_mixed_rows=True,complete_TT_block=False,inverse_or_solution_certificate=False))
    FIXTURES.append(dict(kind='full_finite_domain_weak_assembly',directions=5,pairings=15,
        exact_polynomial_log_primitives=True,finite_trial_space_completeness_proved=False))

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
        check('previous_parametrix_passed',old['status']=='passed' and old['contact_complete_local_parametrix_executed'])
        check('previous_scientific_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        algebra();execute()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths};depafter=endpoint.dependencies()
    status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,
        finished_utc=datetime.now(timezone.utc).isoformat(),actual_finite_local_action_assembled=status=='passed',
        complete_scalar_and_mixed_weak_rows_implemented=status=='passed',signed_trace_conversion_corrected=status=='passed',
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
