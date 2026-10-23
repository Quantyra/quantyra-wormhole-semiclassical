"""Current h1 interior modulus and feasibility budgets; not a metric solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as F
from math import factorial, comb
import json, traceback, sys
import sympy as sp
sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_seed_radial_green_enclosure as prior
from flint import arb,ctx

ROOT=Path(__file__).resolve().parents[1]
STEM='current-reference-regularity-feasibility'
PREVIOUS='current-seed-diagonal-restoration-contract'
OWN=['notes/'+STEM+'.md','checks/current_reference_regularity_feasibility.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[];NUMERICS=[]
encoded=prior.encoded
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)
prior.check=check
def A(x):
    if isinstance(x,F):return arb(x.numerator)/x.denominator
    return arb(x)

def cap_source():
    rho=F(1,500);T=[F('.002'),F('.697'),F('1445.21'),F(6143635),F(1000)+4*315*100+6*700000*3+4*3000000000*F('.067')+4000000000000*F('.002')]
    gp=sum(T[j]*rho**j/factorial(j) for j in range(5))
    cp=3*rho**2/2+100*rho**3/6+1000*rho**4/24
    check('fourth_translation_global_value',gp<F('.0205'),gp)
    check('fourth_translation_center',cp<F('7e-6'))
    mass=1/(1-gp-cp)-1
    check('fourth_translation_Born_and_linear_potential',mass<F('.021') and mass*80<2 and 1-gp>F('.979'))
    gradient=3*rho+100*rho**2/2+1000*rho**3/6+(F('.002')*factorial(5)*(96*F('1.415'))**5)*rho**4/24
    check('fourth_translation_inner_gradient',gradient<F('.02'),gradient)
    for k in range(2,6):
        val=F('.002')*F('.96')**k*sum(comb(k+j,j)*(96*F('1.415')*rho)**j for j in range(1,5))
        check('fourth_translation_inner_jet_'+str(k),val<F('.011'),val)
    j=4;bb=600*(j+3);ll=512*(j+3)
    for k in range(j+2):
        v=F(11,10)*(F(ll,bb)+F(21,10000)*sum(F(factorial(k),factorial(k-b))*F(128,bb)**b for b in range(1,k+1)))
        check('fixed_real_cap_row_four_'+str(k),v<1)
    cs=[F(10)]+[F(180*(j+3)*(600*(j+3))**j) for j in range(1,5)]
    pref=8*256**2*88**2*F(31,2500000)
    gj=[pref*sum(comb(j,k)*cs[k]*cs[j-k] for k in range(j+1))/(factorial(j)*100**j) for j in range(5)]
    lj=[F('2.1e12')]+[F(2**((j+1)//2)*j**j)*F('2.1e12')/(factorial(j)*100**j*rho**j) for j in range(1,5)]
    coord=[F('1e14')]+[(gj[j]+lj[j]+F('1e12'))*factorial(j)*100**j for j in range(1,5)]
    raw=[coord[0],coord[1],coord[2],coord[3]+3*coord[1],coord[4]+12*coord[2]+100*coord[1]]
    src=[x*2*F('8.5e-33')*400**2 for x in raw]
    check('proper_fourth_source_below_one',src[4]<1,src[4])
    NUMERICS.append(dict(kind='current_fixed_cap_fourth_source',global_translation=str(gp),gluing_jets=list(map(str,gj)),local_jets=list(map(str,lj)),proper_source_jets=list(map(str,src))))
    FIXTURES.append(dict(fixture='same_cap_fourth_translation_and_actual_domain_gluing',continuum_Born_argument_replaced_by_checks=False))
    return list(map(A,src))

def symbolic_and_local():
    z=sp.symbols('z');a,d,b,c,e=sp.symbols('a d b c e');hm=F(2000)**-10
    vals=[F(1)/F('8e9'),F('8.6e-19'),7*hm**2,F('4.91e-10'),F('3.4e-22')]
    def cb(p):
        return sum(abs(F(int(sp.numer(v)),int(sp.denom(v))))*product(x**j for x,j in zip(vals,exps)) for exps,v in sp.Poly(sp.expand(p),a,d,b,c,e).terms())
    def ac(p):return sum(cb(v) for v in sp.Poly(sp.expand(p),z).all_coeffs())
    dt=lambda f:sp.expand(-z*z*sp.diff(f,z))
    fn=z**-2+a-(4+2*a+b)*z+(3+a+d*a)*z*z-2*d*a*z**3+d*a*z**4
    wr=4*z**5-6*z**6;rn=-12+d*a*wr;tn=(3+a+d*a)*z**4-d*a*(5*z**5-sp.Rational(9,2)*z**6)
    hr=sp.expand(-2*rn*tn+2*(dt(fn)/2+2*z*fn)*dt(rn));ht=sp.expand(hr+2*fn*dt(dt(rn)))
    J=(24*c*(3+a+d*a)+(1-24*c)*e)*z**4-4*d*a*z**5+3*d*a*z**6
    local=F(2)/F('5.9')*ac(sp.expand((J-c*(hr+ht)/2)/z**4))
    check('R2_reference_cancellation_complex_polynomial',local<F('1.16e-22'),local)
    cn=sp.expand(-z**3*sp.diff(fn,z)+rn/6)
    ps=2*(6+d*a)**2+sp.Rational(22,3)*(6+d*a)*a+sp.Rational(22,3)*a*a+22
    bp=sp.expand(7*cn**2+sp.Rational(22,3)*tn**2-sp.Rational(11,72)*(rn**2-144)+ps*wr/2)
    aa=-10*z**4-20*z**5+30*z**6;lr=sp.expand(z**4*sp.diff(fn*sp.diff(wr,z),z))
    weights=[F('1.001')/(60*F('3.14159')*F(2000)**8),hm/(60*F('3.14159')),F('2.946e-9')*vals[0]*vals[1]]
    trace=sum(w*ac(sp.expand(p/z**2)) for w,p in zip(weights,[aa,bp,lr]))/F('11.8')
    check('complete_trace_complex_polynomial',trace<F('1e-24'),trace)
    beta,eps,tau,m,l,ell=sp.symbols('beta eps tau m l ell')
    q=beta*(1-4*z**3+3*z**4)+a*z**2*(1-z)**2+eps*z**4*(1-z)**2-tau*z**3
    f=sp.expand(q/z**2);R=-sp.diff(q,z,2);Tn=-z*sp.diff(f,z);sigma=-l+ell
    dz=lambda x:sp.expand(sp.diff(x,z)+sp.diff(x,l)/z)
    dtt=lambda x:sp.expand(-z*z*dz(x))
    H=sigma*(R**2/60-a*R/18+a*a/15-m*m/3)+R*Tn/30-Tn*Tn/10+11*Tn*f/90-11*f*f/180-a*f/18
    Hr=sigma*(R/30-a/18)+Tn/30;Ht=R/30-Tn/5+11*f/90
    E=H-R*Hr-(z*z*dz(f)+2*z*f)*dz(Hr)-2*z*z*f*dz(dz(Hr))+z*f*dz(Ht)+(z*dz(f)+f)*Ht
    fj=[f]
    for i in range(4):fj.append(dtt(fj[-1]))
    Ef=(sigma*fj[4]/z**2+(4*sigma+2)*fj[3]/z+2*fj[2]-2*z*fj[1])/30
    W=sp.expand(E-f*Ef/z**2)
    vv=[F(1)/F('5.9'),F('1e-10'),F('1e-24'),F('1e-65'),F(400),F(1),F('1e-8')]
    norm=F(0)
    for exps,co in sp.Poly(sp.expand(z**20*W),z,beta,a,eps,tau,m,l,ell).terms():
        power=exps[0]-20;zb=F('.92')**power if power>=0 else F('.45')**power
        norm+=abs(F(int(sp.numer(co)),int(sp.denom(co))))*zb*product(v**j for v,j in zip(vv,exps[1:]))
    weyl=F('8.5e-33')/F('3.14159')*norm
    full=local+F('2e-35')/F('.45')**4+weyl
    check('complete_local_source_complex_disk',full<F('1e-20'),dict(weyl=str(weyl),total=str(full)))
    # Independent conversion of the Euler-alpha formula from t to z.
    f0=sp.Function('f')(z);h=sp.Function('h')(z);k=sp.Function('k')(z)
    direct=-(dt(f0)/z**2+2*f0/z)*dt(h)-2*f0/z**2*dt(dt(h))-f0/z*dt(k)-(dt(f0)/z-f0)*k
    target=-(z*z*sp.diff(f0,z)+2*z*f0)*sp.diff(h,z)-2*z*z*f0*sp.diff(h,z,2)+z*f0*sp.diff(k,z)+(z*sp.diff(f0,z)+f0)*k
    check('Weyl_Euler_coordinate_conversion',sp.simplify(direct-target)==0)
    Q=sp.Function('Q')(z);rp=sp.Function('rpoly')(z);nu,fc,S,T,u=sp.symbols('nu fc S T u')
    curv=-2*T-6*fc+3*z*z*S+nu*z*sp.diff(Q,z,3)
    curv4=-2*T-nu*sp.diff(rp,z)/z
    radius=(curv-curv4)/6-nu*rp/(3*z*z)
    check('radius_trace_cancellation_and_Lie_term',sp.simplify(radius-(-fc+z*z*S/2+nu*(z*sp.diff(Q,z,3)/6+sp.diff(rp,z)/(6*z)-rp/(3*z*z))))==0)
    check('shape_curvature_reconstruction',sp.expand(-curv/2-(T+3*fc-sp.Rational(3,2)*z*z*S-nu*z*sp.diff(Q,z,3)/2))==0)
    NUMERICS.append(dict(kind='complete_local_complex_bounds',R2_reference=str(local),Weyl=str(weyl),complete_shape=str(full),trace=str(trace),disk_radius='.04'))
    FIXTURES.append(dict(fixture='complete_local_source_and_exact_reconstruction_algebra',boundary_realization_proved=False))

def product(xs):
    out=F(1)
    for x in xs:out*=x
    return out

ORDER=4
class Jet:
    def __init__(self,x):self.c=list(x) if isinstance(x,list) else [A(x)]+[arb(0) for _ in range(ORDER)]
    def __add__(self,b):
        b=b if isinstance(b,Jet) else Jet(b)
        return Jet([x+y for x,y in zip(self.c,b.c)])
    __radd__=__add__
    def __neg__(self):return Jet([-x for x in self.c])
    def __sub__(self,b):return self+-tojet(b)
    def __rsub__(self,b):return tojet(b)+-self
    def __mul__(self,b):
        b=tojet(b);return Jet([sum((self.c[j]*b.c[k-j] for j in range(k+1)),arb(0)) for k in range(ORDER+1)])
    __rmul__=__mul__
    def inv(self):
        out=[1/self.c[0]]
        for k in range(1,ORDER+1):out.append(-sum((self.c[j]*out[k-j] for j in range(1,k+1)),arb(0))/self.c[0])
        return Jet(out)
    def __truediv__(self,b):return self*tojet(b).inv()
    def __rtruediv__(self,b):return tojet(b)*self.inv()
    def __pow__(self,n):
        if n<0:return self.inv()**(-n)
        out=Jet(1)
        for j in range(n):out=out*self
        return out
    def sqrt(self):
        out=[self.c[0].sqrt()]
        for k in range(1,ORDER+1):out.append((self.c[k]-sum((out[j]*out[k-j] for j in range(1,k)),arb(0)))/(2*out[0]))
        return Jet(out)
def tojet(x):return x if isinstance(x,Jet) else Jet(x)
def poly(cs,z):
    out=Jet(0)
    for c in reversed(cs):out=out*z+c
    return out
def polyd(cs,n):
    for j in range(n):cs=[(k+1)*cs[k+1] for k in range(len(cs)-1)]
    return cs
def absolute_jets(s):return [(s.c[j]*factorial(j)).abs_upper() for j in range(5)]

def local_proper(zj,M):
    x=[(zj.c[j]*factorial(j)).abs_upper() for j in range(5)]
    f=[arb(M)*factorial(j)/arb('.04')**j for j in range(5)]
    return [f[0],f[1]*x[1],f[2]*x[1]**2+f[1]*x[2],f[3]*x[1]**3+3*f[2]*x[1]*x[2]+f[1]*x[3],f[4]*x[1]**4+6*f[3]*x[1]**2*x[2]+f[2]*(3*x[2]**2+4*x[1]*x[3])+f[1]*x[4]]

def radial_bound(par,canonical,cells):
    beta,a,e,t=par['beta'],par['u'],par['eps'],par['tau']
    qc=[beta,arb(0),a,-4*beta-2*a-t,3*beta+a+e,-2*e,e]
    rc=[-12*beta,arb(0),arb(0),arb(0),arb(0),4*e,-6*e]
    h5=arb(0);h6=arb(0);records=[]
    for index in range(cells):
        lo=arb('.49')+arb('.39')*index/cells;hi=arb('.49')+arb('.39')*(index+1)/cells
        zj=Jet(prior.hull(lo,hi))
        for k in range(4):zj.c[k+1]=-poly(qc,zj).sqrt().c[k]/(k+1)
        Q=poly(qc,zj);root=Q.sqrt();Qp=poly(polyd(qc,1),zj);Qpp=poly(polyd(qc,2),zj);Qppp=poly(polyd(qc,3),zj)
        assert Q.c[0]>0
        H=-Qp/(2*root);L=root/zj;rp=poly(rc,zj);rpp=poly(polyd(rc,1),zj)
        aa=[[Jet(0) for j in range(6)] for i in range(6)];b=[Jet(0) for i in range(6)];d=[Jet(0) for i in range(6)]
        aa[0][1]=Jet(1)
        aa[1][0]=Qpp;aa[1][1]=-H;aa[1][2]=-zj*Qppp/2;aa[1][4]=Jet(3);b[1]=-F(3,2)*zj**2;d[1]=Jet(1)
        aa[2][3]=Jet(1)
        aa[3][0]=2*(2*Q/zj**2-Qp/zj);aa[3][2]=zj*Qppp/6+rpp/(6*zj)-rp/(3*zj**2);aa[3][3]=-H-2*L;aa[3][4]=Jet(-1);b[3]=zj**2/2
        aa[4][4]=-L;aa[4][5]=-2*L;b[4]=L*zj**2
        aa[5][5]=-2*(H+L);b[5]=H*zj**2/2;d[5]=L
        aj=[[absolute_jets(x) for x in row] for row in aa];bj=list(map(absolute_jets,b));dj=list(map(absolute_jets,d))
        extra=local_proper(zj,'1e-20');ss=[x+y for x,y in zip(canonical,extra)];tt=local_proper(zj,'1e-24')
        ws=[[arb(s) for s in ['1.2e-8','7.7e-6','1.2e-8','3e-7','8e-10','1e-10']]]
        for k in range(5):
            ws.append([sum((comb(k,j)*(sum((aj[i][v][j]*ws[k-j][v] for v in range(6)),arb(0))+bj[i][j]*ss[k-j]+dj[i][j]*tt[k-j]) for j in range(k+1)),arb(0)) for i in range(6)])
        h5=arb(max(h5.upper(),ws[4][1].upper()));h6=arb(max(h6.upper(),ws[5][1].upper()))
        records.append(dict(cell=index,H5_upper=encoded(ws[4][1]),H6_upper=encoded(ws[5][1])))
    return h5,h6,records

def execute_geometry(canonical):
    p=prior.parameters()
    check('reference_N_low_bound',arb('.2').sqrt()/arb('.49')*(arb('3e-7')+arb('.88')*arb('1.2e-8'))<arb('3e-7'))
    check('reference_f_low_bound',arb('7.5e-10')*arb('.2')/arb('.49')**2+arb('1e-30')<arb('8e-10'))
    check('reference_pressure_low_bound',arb('.65e-10')/arb('5.9')*(arb('.49')**-2-1)+arb('1e-30')<arb('1e-10'))
    # Independent elementary series identities, including the reciprocal and root.
    j=Jet([arb(2),arb(3),arb(-1),arb(2),arb(1)])
    for name,v,expect in [('reciprocal',j*j.inv(),Jet(1)),('root_square',j.sqrt()**2,j)]:
        check('interval_series_'+name,all((x-y).contains(0) for x,y in zip(v.c,expect.c)))
    kap=4*p['C']/(5*arb.pi()*p['B2']);delta=arb('1e-5');threshold=6*arb('1e-28')/(kap*delta.sqrt())
    for cells in [64,128]:
        h5,h6,records=radial_bound(p,canonical,cells)
        holder=(2*arb('2e8')*h6).sqrt();actual_holder=(2*h5*h6).sqrt()
        err=kap*holder*delta.sqrt()/6;lip=kap*h6*delta/12
        check('current_reference_modulus_target_'+str(cells),holder<threshold,encoded(holder))
        check('principal_derivative_error_'+str(cells),err<arb('1e-28') and lip<arb('1e-28'))
        NUMERICS.append(dict(kind='current_reference_modulus',cells=cells,region_z=['.49','.88'],H5_upper=encoded(h5),H6_upper=encoded(h6),Holder_half_using_old_H5_upper=encoded(holder),Holder_half_using_current_H5_upper=encoded(actual_holder),Holder_half_threshold=encoded(threshold),principal_error_Holder=encoded(err),principal_error_Lipschitz=encoded(lip),threshold_margin_factor=encoded(threshold/holder),cell_bounds=records))
        if cells==128:
            previous=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
            majorant=next(r for r in previous['numerics'] if r['kind']=='all_mode_kernel_majorant')['kappa']
            def decoded(x):
                enc=x['decimal_enclosure'];return arb(enc['mid'],enc['rad'])*arb(10)**int(enc['exponent'])
            ks={int(j):decoded(x) for j,x in majorant.items()}
            for label in ['.00001','.01']:
                cut=arb(label);value=h5/60*sum((ks[j]*cut**(6-j)/(6-j) for j in ks),arb(0));derivative=kap*h6*cut/12
                resolution=next(r for r in previous['numerics'] if r.get('kind')=='near_diagonal_budget' and r['delta']==label)
                check('actual_direction_simultaneous_integration_budget_'+label,value<arb('1e-28') and derivative<arb('1e-28'))
                NUMERICS.append(dict(kind='current_direction_simultaneous_integration_budget',delta=label,all_mode_value_remainder_upper=encoded(value),principal_Lipschitz_derivative_remainder_upper=encoded(derivative),angular_32_attenuation_scale=resolution['angular_index_scale'],thermal_32_attenuation_scale=resolution['thermal_index_scale'],scale_counts_are_truncation_tolerance=False,full_differentiated_remainder_bounded=False))
    FIXTURES.append(dict(fixture='current_reference_radial_jet_and_principal_modulus_test',actual_full_quantum_Hessian_evaluated=False))

def main():
    ctx.prec=256
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};deps=prior.dependency_identity();error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'],dependencies=deps,invocation=sys.argv)
    try:
        check('previous_diagonal_audit_passed',old['status']=='passed' and old['C5_to_C1_without_modulus_rejected'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        check('pinned_interval_library',prior.flint.__version__=='0.9.0')
        canonical=cap_source();symbolic_and_local();execute_geometry(canonical)
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths};depafter=prior.dependency_identity();status='passed' if error is None and before==after and deps==depafter else 'failed'
    run.update(status=status,exception=error,input_sha256_after=after,unchanged_inputs=before==after,unchanged_dependencies=deps==depafter,checks=ROWS,fixtures=FIXTURES,numerics=NUMERICS,finished_utc=datetime.now(timezone.utc).isoformat(),actual_reference_middle_fifth_derivative_modulus_quantified=status=='passed',current_reference_local_sixth_derivative_bounded=status=='passed',principal_derivative_quadrature_target_certified=status=='passed',all_local_renormalized_moments_computed=False,coincident_subtraction_executed=False,differentiated_response_operator_bounded=False,renormalized_quantum_source_evaluated=False,full_quantum_Hessian_evaluated=False,actual_physical_correction_computed=False,raw_C5_admission_certified=False,nonlinear_solution_constructed=False,physical_stability_proved=False,semiclassical_validity_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]};receipt['runs'].append(run)
    out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=status,exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=status,checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)),flush=True)
    if status!='passed':raise SystemExit(1)

if __name__=='__main__':main()
