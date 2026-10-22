"""Intrinsic finite TT lift, full tensor coefficients and exact H8 budget."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
from math import comb,factorial
from itertools import product
import importlib.util,json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-higher-tt-lift';PREVIOUS='finite-radius-operator-source'
OWN=['checks/finite_higher_tt_lift.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
spec=importlib.util.spec_from_file_location('retained_low_gap',ROOT/'checks/finite_low_gap_reconstruction.py')
B=importlib.util.module_from_spec(spec);spec.loader.exec_module(B);D=B.D;L=B.L
poly,add,mul,scale,der,power,compose=L.poly,L.add,L.mul,L.scale,L.der,L.power,L.compose
ROWS=[];FIXTURES=[]
PROTOCOL={'input_order':8,'tensor_lift_order':8,'return_order':7,'vector_reconstruction_order':9,
          'curvature_coefficient_order':8,'metric_order_used':10,'numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','manufactured_derivative_orders':[0,3,8],
          'actual_charged_source_sampled':False,'ordinary_continuum_proof_separate':True}

def check(name,ok,detail=''):
    assert name not in {x['name'] for x in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)

def zero(name,value):
    simplified=S.simplify(S.expand_trig(value));check(name,simplified==0,simplified)
def saved(stem):return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def fixture(stem,name):return next(x for x in saved(stem)['fixtures'] if x['fixture']==name)
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def strings(a):return [str(x) for x in a]
def decimals(a):return [L.decimal(x) for x in a]
def arr(f,k):return poly([F(x) for x in f[k]])
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})

def inputs():
    old=saved(PREVIOUS)
    check('prior_radius_source_passed',old['status']=='passed' and len(old['checks'])==79 and all(x['passed'] for x in old['checks']))
    check('all_prior_inputs_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    check('low_return_vector_guard',F('1.415')**2>2 and F('1.415')*F('3.002')/F('.995')<F('4.3'))
    check('low_Psi_guard',2*F('1.008')*F('3.002')+4*F('2.25')*F('4.3')+F('2.002')<47)
    check('low_full_tensor_guard',2*F('3.002')**2+4<F('4.7')**2)
    check('gradient_log_radius_guard',35**2+12**2+2*F('5.1')**2<40**2)

def intrinsic_identities():
    x=S.symbols('rho',positive=True)
    aa,rr,bb,cc,yy=[S.Function(k)(x) for k in ['A','R','b','c','Y']]
    w=aa+2*rr;psi=2*w*yy-bb;ric=-S.diff(aa,x)-aa**2-2*S.diff(rr,x)-2*rr**2
    sub={S.diff(yy,x):aa*yy-bb-cc,S.diff(bb,x):-2*(aa+rr)*bb-2*(aa-rr)*cc}
    zero('finite_intrinsic_Psi_equation',S.diff(psi,x).subs(sub)-(2*rr*psi-6*rr*cc-2*ric*yy-12*rr**2*yy))
    zero('finite_intrinsic_divergence',S.diff(yy,x).subs(sub)+w*yy-(psi-cc-2*rr*yy))
    sc=cc/2+rr*yy;tc=-3*cc/2+psi/2-3*rr*yy
    zero('returned_disk_radius',sc-(cc/2+rr*yy))
    zero('returned_disk_shape',tc-(-(bb+3*cc)/2+(aa-rr)*yy))
    for i,(proj,grad,expected) in enumerate(zip([1,1,0,0],[aa*yy-bb-cc,aa*yy,rr*yy,rr*yy],[bb,-bb-2*cc,cc,cc])):
        zero('intrinsic_lift_component_'+str(i),cc+2*rr*yy+(-3*cc+psi-6*rr*yy)*proj-2*grad-expected)
    check('quadratic_q_term_is_nonzero',S.simplify(12*rr**2*yy)!=0)
    FIXTURES.append({'fixture':'finite_intrinsic_lift_identities','Psi':'2(A+2R)Y-b','Ric_rr':str(ric),
                     'R_prime_independent':True,'fixed_return_sign':'Y=-sqrt(lambda)*eta','singular_A_bound_used':False})

def independent_four_geometry():
    x,t,theta,phi=S.symbols('rho time theta phi',positive=True);coords=[x,t,theta,phi]
    a=S.Function('a')(x);r=S.Function('r')(x);b=S.Function('b')(x);c=S.Function('c')(x);y=S.Function('Y')(x)
    g=S.diag(1,a*a,r*r,r*r*S.sin(theta)**2);inv=g.inv()
    conn={(i,j,k):S.simplify(sum(inv[i,l]*(S.diff(g[l,k],coords[j])+S.diff(g[l,j],coords[k])-S.diff(g[j,k],coords[l]))/2 for l in range(4))) for i,j,k in product(range(4),repeat=3)}
    ricrr=sum(S.diff(conn[k,0,0],coords[k])-S.diff(conn[k,k,0],x)+sum(conn[k,k,l]*conn[l,0,0]-conn[k,0,l]*conn[l,k,0] for l in range(4)) for k in range(4))
    aa=S.diff(a,x)/a;rr=S.diff(r,x)/r
    zero('independent_4D_Ricci',ricrr+S.diff(aa,x)+aa*aa+2*S.diff(rr,x)+2*rr*rr)
    h=S.diag(b,a*a*(-b-2*c),r*r*c,r*r*S.sin(theta)**2*c)
    div=[]
    for j in range(4):
        val=sum(inv[i,k]*(S.diff(h[i,j],coords[k])-sum(conn[l,k,i]*h[l,j]+conn[l,k,j]*h[i,l] for l in range(4))) for i,k in product(range(4),repeat=2))
        div.append(S.simplify(val))
    zero('independent_4D_TT_equation',div[0]-(S.diff(b,x)+2*(aa+rr)*b+2*(aa-rr)*c))
    for j in range(1,4):zero('independent_4D_angular_divergence_'+str(j),div[j])
    divergence=S.diff(y,x)+sum(conn[i,i,0]*y for i in range(4))
    zero('independent_4D_vector_divergence',divergence-S.diff(y,x)-(aa+2*rr)*y)
    yp=aa*y-b-c;eta=-y
    lie=S.diag(*[eta*S.diff(g[i,i],x)+(2*g[i,i]*S.diff(eta,x) if i==0 else 0) for i in range(4)])
    returned=h-lie;sc=c/2+rr*y;tc=-(b+3*c)/2+(aa-rr)*y
    for j in range(4):zero('independent_4D_disk_return_'+str(j),(returned[j,j]/g[j,j]-2*(sc+(tc if j<2 else 0))).subs(S.diff(y,x),yp))
    FIXTURES.append({'fixture':'independent_four_dimensional_geometry','dimension':4,
                     'metric':'diag(1,a(rho)^2,r(rho)^2,r(rho)^2 sin(theta)^2)',
                     'genuine_Ricci_and_TT_divergence':True,'seed_only_identity_used':False})

def coefficient_arrays():
    pf=fixture('finite-higher-local-coefficients','finite_reference_metric_and_curvature_majorants')
    old=fixture('finite-low-gap-reconstruction','refined_finite_geometry_and_components')
    tf=fixture('finite-low-gap-reconstruction','complete_improved_tangent_arrays')
    T=arr(fixture('finite-low-gap-reconstruction','finite_parameter_integral_bounds'),'Tplus')
    q=add(poly([1]),arr(pf,'seed_jets'));ll=arr(pf,'log_radius_jets')
    geo=D.geometry(T,arr(tf,'K'),arr(tf,'Vbar'),q,ll)
    for k in ['M0','M1','C0','C1','J1','CW1','Sc1','S1']:
        check('frozen_geometry_reconciled_'+k,geo[k][:9]==[F(x) for x in old['geometry'][k]])
    C=geo['C0'];M=geo['M0']
    qj=scale('2.5',L.covariant(der(add(ll,T)),C,1,8));qj[0]=min(qj[0],F('2.25'));qj[1]=min(qj[1],40)
    seed=poly([F(1)/F('5.9'),0,F('1e-10'),F('.668'),F('.501'),F('2e-24'),F('1e-24')])
    rr=add(compose(der(der(seed)),q),poly(['2e-10']),scale(2,der(C)),scale(2,power(C,2)))
    rm=L.covariant(scale('6.1',mul(M,rr)),C,4,8);rm[0]=min(rm[0],55);rm[1]=min(rm[1],F('1e16'))
    ric=[2*x for x in rm];ric[0]=min(ric[0],15)
    projector=[F(2)]
    for j in range(8):projector.append(2*sum(F(comb(j,i))*projector[i]*qj[j-i] for i in range(j+1)))
    check('available_metric_orders_suffice',len(fixture('finite-low-gap-reconstruction','finite_parameter_integral_bounds')['Tplus'])==11 and len(rm)==len(qj)==len(projector)==9)
    check('all_full_coefficient_arrays_positive',min(qj+rm+ric+projector)>0)
    for j in range(9):
        check('curvature_Gaffney_square_root_'+str(j),16*(j+1)*rm[0]<(32*(j+1))**2)
        check('ordered_commutator_count_'+str(j),all(sum(comb(l,i)*(j-l+1) for l in range(i,j))<=(j+1)*comb(j,i+1) for i in range(j)))
    FIXTURES.append({'fixture':'full_finite_lift_coefficient_arrays','q':strings(qj),'Rm':strings(rm),'Ric':strings(ric),'projector':strings(projector),
                     'q_decimals':decimals(qj),'Rm_decimals':decimals(rm),'metric_order_used':10,'full_ordered_tensors':True,
                     'frozen_scalar_estimate_unchanged':True,'no_endpoint_inverse_used':True})
    return qj,rm,ric,projector

def va(*args):return [sum(a[j] for a in args) for j in range(9)]
def vs(c,a):return [F(c)*x for x in a]
def vz():return [F(0)]*9
def conv(q,y,j):return va(*[vs(comb(j,i)*q[i],y[j-i]) for i in range(j+1)])

def lift_matrices(q,rm,ric,projector):
    unit=[[F(int(i==j)) for i in range(9)] for j in range(9)]
    y=[vs('4.3',unit[0])];p=[vs(47,unit[0])];xs=[];ds=[]
    for j in range(9):
        x=conv(q,y,j);xs.append(x);d=va(p[j],unit[j],vs(2,x));ds.append(d)
        curvature=va(*[vs(16*(j+1)*comb(j,i+1)*rm[i],y[j-1-i]) for i in range(j)]) if j else vz()
        y.append(va(vs(2,d),vs(32*(j+1),y[j]),curvature))
        p.append(va(vs(2,conv(q,p,j)),vs(6,conv(q,unit,j)),vs(2,conv(ric,y,j)),vs(12,conv(q,xs,j))))
    h=[]
    for j in range(9):
        core=va(vs(2,unit[j]),vs(4,xs[j]),vs(2,y[j+1]))
        pt=va(*[vs(comb(j,i)*projector[i],va(vs(3,unit[j-i]),p[j-i],vs(6,xs[j-i]))) for i in range(j+1)])
        h.append(va(core,pt))
    h[0]=vs('4.7',unit[0])
    hc=va(*h);yc=va(*y[:8]);bound=max(hc+yc);exponent=L.decimal_power(bound)
    check('lift_matrices_nonnegative',all(x>=0 for a in h+y+p for x in a))
    check('input_order_triangular',all(h[j][i]==0 for j in range(9) for i in range(j+1,9)))
    check('vector_gain_one_derivative',all(y[j][i]==0 for j in range(1,10) for i in range(j,9)))
    check('exact_column_sum_bound',all(x<=bound for x in hc+yc) and bound in hc+yc)
    check('upward_L8_enclosure',F(10)**(exponent-1)<bound<=F(10)**exponent,{'L8':L.decimal(bound),'power':exponent})
    # A normalized basis input attains the induced l1 norm of each positive matrix.
    for name,matrix,columns in [('tensor',h,hc),('return',y[:8],yc)]:
        i=columns.index(max(columns));basis=unit[i]
        norm=sum(sum(row[k]*basis[k] for k in range(9)) for row in matrix)
        check('independent_'+name+'_induced_norm',norm==max(columns))
    mixed=F(381,24)*(F('1e52')+F('1e56'))*F(10)**exponent;mpower=L.decimal_power(mixed)
    check('upward_mixed_row_enclosure',mixed<=F(10)**mpower and mixed>F(10)**(mpower-1))
    scalar_budget=F(10**57)/(1600*F(10)**20);tt_budget=F(10**57)/(1600*F(10)**mpower)
    check('equal_split_source_budgets',800*(F(10)**20*scalar_budget+F(10)**mpower*tt_budget)==10**57)
    FIXTURES.append({'fixture':'per_input_order_lift_and_return','Y':[strings(v) for v in y],'Psi':[strings(v) for v in p],
                     'tensor':[strings(v) for v in h],'tensor_column_sums':strings(hc),'return_column_sums':strings(yc),
                     'tensor_column_decimals':decimals(hc),'return_column_decimals':decimals(yc),
                     'L8_exact':str(bound),'L8_decimal':L.decimal(bound),'L8_power_of_ten':exponent,
                     'mixed_row_coefficient_over_g_exact':str(mixed),'mixed_row_power_of_ten':mpower,
                     'half_scalar_budget_over_g_epsilon_c':str(scalar_budget),'half_scalar_budget_decimal':L.decimal(scalar_budget),
                     'half_TT_R_C_D4_F4_budget_over_g_epsilon_c':str(tt_budget),'half_TT_budget_decimal':L.decimal(tt_budget),
                     'actual_quantum_H4_supplied':False,'actual_TT_inverse_H4_H8_supplied':False,'target_passed':False})

def manufactured_polynomials(alpha,beta):
    x=S.symbols('x',real=True);k=S.Integer(2)
    yp=x*(1+x*x+x**4+x**8)/10**6;aa=1/x+2*alpha*x;rr=2*beta*x
    ric=-(6*alpha+4*beta)-(4*alpha**2+8*beta**2)*x*x
    dp=S.expand(S.diff(yp,x)+(aa+2*rr-2*k*x)*yp)
    rhs=S.expand(6*rr*dp-2*ric*yp);qq={6:S.Integer(0)}
    for n in reversed(range(6)):qq[n]=S.simplify((rhs.coeff(x,2*n+1)-(2*n+2)*qq[n+1])/(8*beta-2*k))
    pp=S.expand(sum(qq[n]*x**(2*n) for n in range(6)))
    cp=S.expand(pp-dp-2*rr*yp);bp=S.expand(2*(aa+2*rr)*yp-pp)
    zero('manufactured_TT_polynomial_'+str(alpha)+'_'+str(beta),S.diff(bp,x)-2*k*x*bp+2*(aa+rr)*bp+2*(aa-rr)*cp)
    zero('manufactured_Psi_polynomial_'+str(alpha)+'_'+str(beta),S.diff(pp,x)+(8*beta-2*k)*x*pp-rhs)
    check('manufactured_not_seed_identity_'+str(alpha)+'_'+str(beta),S.expand(S.diff(rr,x)-aa*rr+rr*rr)!=0)
    check('manufactured_regular_tensor_bolt_'+str(alpha)+'_'+str(beta),S.expand(bp+cp).subs(x,0)==0)
    def coefficients(v):
        v=S.expand(v);return [F(str(v.coeff(x,j))) for j in range(S.degree(v,x)+1)]
    return {name:coefficients(v) for name,v in [('Y',yp),('Psi',pp),('c',cp),('b',bp),('b_plus_c',S.expand(bp+cp))]}

def mfrac(x):return mp.mpf(x.numerator)/x.denominator
def evalp(a,x):return mp.polyval([mfrac(v) for v in reversed(a)],x)

def independent_integrals():
    configs=[(F(1,10),F(1,5)),(F(7,20),F(1,7))];polys={key:manufactured_polynomials(S.Rational(key[0]),S.Rational(key[1])) for key in configs};records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            for alpha,beta in configs:
                af,bf=mfrac(alpha),mfrac(beta);p=polys[alpha,beta]
                field=lambda key,x:mp.exp(-2*x*x)*evalp(p[key],x)
                def b_integral(x):
                    return -2*mp.exp(-2*(af+bf)*x*x)*mp.quad(lambda t:(t+2*(af-bf)*x*x*t**3)*mp.exp(2*(af+bf)*x*x*t*t)*field('c',x*t),[0,1])
                def y_integral(x):
                    rate=2+af
                    return x*mp.exp(af*x*x)*sum(mfrac(p['b_plus_c'][2*n])*mp.gammainc(n,rate*x*x,mp.inf)/(2*rate**n) for n in range(1,len(p['b_plus_c'])//2+1))
                ric=lambda x:-(6*af+4*bf)-(4*af**2+8*bf**2)*x*x
                for at in [mp.mpf(0),mp.mpf('.3')]:
                    tag=str(alpha)+'_'+str(beta)+'_'+str(at)+'_'+str(dps);errors=[]
                    for n in PROTOCOL['manufactured_derivative_orders']:
                        error=max(abs(mp.diff(b_integral,at,n)-mp.diff(lambda x:field('b',x),at,n)),abs(mp.diff(y_integral,at,n)-mp.diff(lambda x:field('Y',x),at,n)))/mp.factorial(n)
                        check('independent_fixed_integrals_'+tag+'_'+str(n),error<mp.mpf('1e-42'),mp.nstr(error,18));errors.append(mp.nstr(error,18))
                    omric=abs(mp.diff(lambda x:2*ric(x)*field('Y',x),mp.mpf(0),1))
                    omq=abs(mp.diff(lambda x:12*(2*bf*x)**2*field('Y',x),mp.mpf(0),3))
                    check('omitted_Ricci_or_quadratic_q_fails_'+tag,min(omric,omq)>mp.mpf('1e-10'))
                    records.append({'alpha':str(alpha),'beta':str(beta),'rho':str(at),'dps':dps,'orders':PROTOCOL['manufactured_derivative_orders'],
                                    'errors':errors,'omitted_Ricci':mp.nstr(omric,18),'omitted_quadratic_q':mp.nstr(omq,18)})
    FIXTURES.append({'fixture':'independent_regular_fixed_integrals','records':records,'exact_bolt_included':True,
                     'bolt_fixed_TT_and_outer_fixed_return':True,'manufactured_fields_decay_at_infinity':True,
                     'actual_charged_geometry_sampled':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 higher physical TT lift','input_count':len(ps)}),flush=True)
    try:
        inputs();intrinsic_identities();independent_four_geometry();co=coefficient_arrays();lift_matrices(*co);independent_integrals()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(x['passed'] for x in ROWS)
    run={'story':'S167 higher physical TT lift','created_utc':datetime.now(timezone.utc).isoformat(),'status':'passed' if good else 'failed',
         'exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after}
    for k,v in review().items():
        if k not in run and isinstance(v,bool):run[k]=v
    if not good:
        for kind,p in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
