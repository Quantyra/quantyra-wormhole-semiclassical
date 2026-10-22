"""Finite middle massive source: actual-domain spatial derivatives through four."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
from math import factorial,comb
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-middle-source-four';PREVIOUS='finite-higher-york-source'
OWN=['checks/finite_middle_source_four.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
PROTOCOL={'spatial_order':4,'metric_order':9,'energy_order':6,'translation_radius':'1e-28',
          'numerical_dps':[60,90],'gauss_nodes':[48,96],'quadrature_tolerance':'1e-25',
          'actual_charged_determinant_sampled':False,'ordinary_continuum_proof_separate':True}

def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)

def zero(name,value):
    value=S.simplify(value);check(name,value==S.zeros(*value.shape) if isinstance(value,S.MatrixBase) else value==0,value)

def saved(stem):return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})
def strings(xs):return [str(x) for x in xs]
def mf(x):return mp.mpf(x.numerator)/x.denominator if isinstance(x,F) else mp.mpf(x)

def exp_coeff(z,n):
    a=[F(1)]+[F(0)]*n
    for j in range(1,n+1):a[j]=sum(k*z[k]*a[j-k] for k in range(1,j+1))/j
    return a

def retained_and_geometry():
    old=saved(PREVIOUS)
    check('previous_York_source_passed',old['status']=='passed' and len(old['checks'])==95 and all(x['passed'] for x in old['checks']))
    check('previous_inputs_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    check('complete_shape_and_solve_remain_open',not old['actual_complete_shape_source_H4_bounded'] and not old['nonlinear_solution_constructed'])
    check('chart_leaves_positive_outer_distance',F('.051')-F('1.028')/40>F('.018'))
    zz=F(1,200)
    lower=1/F('6.1')-F('.668')*zz**3-2*F('1e-24')*zz**5
    check('small_z_seed_gap',lower>F(1,7))
    check('middle_chart_radius_bound',7*zz**2<F('.018')**2)
    gamma=[F('1e4')*factorial(i)*1000**i for i in range(9)]
    for i in range(9):check('seed_connection_Cauchy_'+str(i),F('.216')*factorial(i+1)*640**(i+1)<gamma[i])
    b={(n,0):F('4e24')*factorial(n)*40000**n for n in range(10)}
    for k in range(9):
        for n in range(9-k):
            b[n,k+1]=b[n+1,k]+n*sum(comb(k,i)*gamma[i]*b[n,k-i] for i in range(k+1))
    check('complete_coordinate_conversion_induction',all(v<=F('4e24')*factorial(n+k)*40000**n*F('1e6')**k for (n,k),v in b.items()))
    check('coordinate_nine_jet_reserve',max(b[0,j] for j in range(10))<F('1e85'))
    step=S.symbols('t',real=True)
    eta=S.expand(sum(S.binomial(23,i)*step**i*(1-step)**(23-i) for i in range(12)))
    zero('eta23_exact_first_derivative',S.diff(eta,step)+16224936*step**11*(1-step)**11)
    check('eta23_C11_endpoint_matching',eta.subs(step,0)==1 and eta.subs(step,1)==0 and all(S.diff(eta,step,j).subs(step,k)==0 for j in range(1,12) for k in [0,1]))
    check('fixed_cap_slope',F(16224936,4**11)*160<640)
    check('fixed_cap_global_gradient',F('.067')+F('8e-6')+640*(F('.002')+F('1.2e-8'))<F('1.35'))
    check('cap_branch_complex_value',F(2)**23*F('1.005')**23<F('1e8'))
    cut=[F('1e8')*factorial(j)*200000**j for j in range(6)]
    check('cap_finite_derivatives',max(cut)<F('1e38'))
    check('uncapped_low_derivatives_on_support',max(F('.027')*factorial(j)*640**j+F('.00013')*factorial(j)*150**j for j in range(6))<F('1e16'))
    check('translated_inner_disk_stays_uncapped',F(1,80)+F('.001')<F(3,160))
    check('translated_local_X_bound',max(F('.027')*2**j+F('.00013') for j in range(2,6))<2)
    check('translation_fourth_X_reserve',F('1e86')+32*F('1e38')*F('1e16')<F('1e88'))
    FIXTURES.append({'fixture':'finite_metric_and_fixed_cap','coordinate_u_jets':strings([b[0,j] for j in range(10)]),
                     'cap_derivative_bounds':strings(cut),'X_value_bound':2,'X_fourth_derivative_bound':'1e88',
                     'actual_spatial_analyticity_assumed':False,'cap_fixed_during_observation_derivatives':True})
    return eta,step

def interpolation_and_cauchy():
    vv=S.Matrix([[k**j for j in range(4)] for k in range(4)]);iv=vv.inv()
    zero('Vandermonde_inverse',iv*vv-S.eye(4))
    costs=[sum(abs(iv[j,k])*S.Rational(2*24+k**4,24) for k in range(4))*factorial(j) for j in range(4)]
    check('Banach_interpolation_constants',max(costs)<1000)
    eps=F('1e-6');pchange=1000*sum(eps**j/F(factorial(j)) for j in range(1,5))
    check('finite_polynomial_parameter_cone',pchange<F('.002') and pchange+F('.00013')<F('.011'))
    check('complex_cap_linear_potential',F('1.36')/(1-F('.007'))<2)
    check('complex_cap_absolute_potential',F('.007')/(1-F('.007'))<F('.06'))
    check('interpolation_real_translation_interval',3*F('1e-22')<F('.001'))
    check('high_jet_and_translation_scale',F('1e88')*F('1e-22')**4==1 and F('1e-28')/F('1e-22')==eps)
    etaG=F('8.5e-33');cap=2*etaG*(F(2,3)*F('1e17')+F('1e17'))
    check('complete_cap_value_budget',cap<F('3e-15'))
    raw=[cap]+[2**((j+1)//2)*j**j*F('3e-15')*F('1e28')**j for j in range(1,5)]
    check('complete_polarized_cap_four_jet',max(raw)<F('1e101'))
    FIXTURES.append({'fixture':'finite_translation_interpolation_and_Cauchy','inverse_matrix':[[str(x) for x in iv.row(j)] for j in range(4)],
                     'interpolation_costs':strings(costs),'parameter_change':str(pchange),'cap_value_bound':str(cap),
                     'cap_coordinate_bounds':strings(raw),'finite_polynomial_only':True,'common_fourth_restoration_retained':True})

def finite_rows_and_conversion():
    check('original_row_evaluation_disk_fits',F('.95')*(F(1,160)-F(1,4096))>F(1,176)>F(1,180))
    check('real_mass_coercivity',1-F('.006')-F('.07')/18>F('.98'))
    check('actual_first_mass_jet', (F('.067')+F('8e-6'))/(1-F('.003'))<F('.07'))
    z=[F(0),F('.07')]+[F('.027')*128**j+F('.00013')*150**j for j in range(2,6)]
    ac=exp_coeff(z,5);ec=exp_coeff([v/2 for v in z],5)
    aa=[F(0),F('.07')]+[F('.1')*factorial(j)*300**j for j in range(2,6)]
    check('finite_mass_exponential_jets',all(F('1.01')*ac[j]*factorial(j)<aa[j] for j in range(2,6)))
    check('finite_endpoint_exponential_jets',all(F('1.01')*ec[j]<F('.1')*300**j for j in range(1,6)))
    c=[F(34000)];ratios=[]
    for j in range(5):
        m=j+3;bb=10000*m;ll=8192*m
        r=[F('1.1')*(F(ll,bb)+sum(F(comb(k,b))*aa[b]/bb**b for b in range(1,k+1))) for k in range(j+2)]
        check('mass_retaining_energy_recurrence_'+str(j),max(r)<1);ratios.append(strings(r))
        point=F(33000)*m*(9+4*F(10000,8192)+2*F(10000,8192)**2)
        end=F('1.01')+F('.1')*sum(F(factorial(j),factorial(j-b))*F(300,bb)**b for b in range(1,j+1))
        check('physical_point_and_endpoint_row_'+str(j),point*end<600000*m)
        if j:c.append(F(600000*m)*bb**j)
    gj=[2*F('8.5e-33')*2*320*24*176**3*sum(comb(j,k)*c[k]*c[j-k] for k in range(j+1)) for j in range(5)]
    check('all_angular_gluing_four_jet',max(gj)<F('1e20'))
    x,delta=S.symbols('E delta',positive=True)
    zero('whole_tower_counting_moment',4*delta*S.integrate(x**3*S.exp(-delta*x),(x,0,S.oo))-24/delta**3)
    gam=[F('.3'),F(4000),F('2e6')]
    check('actual_connection_low_jets',4*F('.07')<gam[0] and 4*factorial(2)*z[2]<gam[1] and 4*factorial(3)*z[3]<gam[2])
    cc={(0,k):F(1) for k in range(5)}
    for r in range(4):
        for k in range(4-r):cc[r+1,k]=cc[r,k+1]+r*sum(comb(k,i)*gam[i]*cc[r,k-i] for i in range(k+1)) if r else cc[r,k+1]
    out=[cc[j,0] for j in range(5)]
    check('actual_canonical_tensor_conversion',max(out)<F('1e8'))
    check('finite_middle_massive_C4_reserve',4*F('1e101')*max(out)<F('1e110'))
    FIXTURES.append({'fixture':'finite_coercive_rows_and_complete_gluing','mass_jets':strings(aa),'energy_ratios':ratios,
                     'physical_row_constants':strings(c),'normalized_gluing_bounds':strings(gj),
                     'full_angular_tower_retained':True,'actual_paired_domain_retained':True})
    FIXTURES.append({'fixture':'actual_canonical_tensor_conversion','connection':strings(gam),'tensor_coefficients':strings(out),
                     'final_uniform_C4_bound':'1e110','quantity':'S_M,ell minus its fixed flat-vacuum scalar',
                     'region':'1.001 <= t_seed <= 40','canonical_orbit_dimension':2,'physical_shape_H4_evaluated':False})

def independent_dirac_energy():
    x,y=S.symbols('x y',real=True);r2=x*x+y*y;records=[]
    sx=S.Matrix([[0,1],[1,0]]);sy=S.diag(1,-1)
    for m,alpha,beta in [(18,S.Rational(1,10),S.Rational(1,5)),(37,S.Rational(1,7),S.Rational(1,3))]:
        a=1+alpha*x+beta*r2;w=S.Matrix([1,S.Rational(1,3)])
        dx=-2*x*w;dy=-2*y*w;dw=sx*dx+sy*dy;cross=(w.T*(sx*S.diff(a,x)+sy*S.diff(a,y))*w)[0]
        lhs=((dw+m*a*w).T*(dw+m*a*w))[0]
        rhs=(dx.T*dx+dy.T*dy)[0]+m*m*a*a*(w.T*w)[0]-m*cross
        def integral(poly):
            ans=0
            for (i,j),v in S.Poly(S.expand(poly),x,y).terms():
                if i%2 or j%2:continue
                ans+=v*S.gamma(S.Rational(i+1,2))*S.gamma(S.Rational(j+1,2))/2**S.Rational(i+j+2,2)
            return S.simplify(ans)
        lv,rv=integral(lhs),integral(rhs)
        zero('independent_Gaussian_Dirac_energy_'+str(m),lv-rv)
        check('wrong_cross_sign_detected_'+str(m),integral(2*m*cross)!=0)
        check('nonconstant_mass_positive_'+str(m),1-alpha**2/(4*beta)>0)
        records.append({'M':m,'alpha':str(alpha),'beta':str(beta),'lhs':str(lv),'rhs':str(rv),'wrong_sign_difference':str(integral(2*m*cross))})
    FIXTURES.append({'fixture':'independent_Dirac_energy_identity','records':records,'spatial_dimension':2,'cross_term_retained':True})

def independent_cap_translation(eta,step):
    x=S.symbols('x',real=True);a=S.Rational(3,160);b=S.Rational(1,40);phi=1+x*x+x**4
    taper=eta.subs(step,160*x-3);kap4=16*1000**4*x**4-48*1000**3*x*x+12*1000**2
    records=[]
    pieces=[]
    for lo,hi,chi in [(0,a,S.Integer(1)),(a,b,taper)]:
        field=S.expand(chi*phi)
        polys=[S.expand(kap4*field),S.diff(field,x,4),S.expand(chi*S.diff(phi,x,4))]
        pieces.append((str(lo),str(hi),[S.lambdify(x,p,'mpmath',cse=True) for p in polys]))
    for dps in PROTOCOL['numerical_dps']:
        with mp.workdps(dps):
            previous=None
            for count in PROTOCOL['gauss_nodes']:
                nodes,weights=mp.gauss_quadrature(count,'legendre');vals=[mp.mpf(0)]*3
                for lo,hi,fs in pieces:
                    low=mf(F(lo));high=mf(F(hi));mid=(low+high)/2;half=(high-low)/2
                    for j in range(count):
                        at=mid+half*nodes[j];factor=2*half*weights[j]*mp.exp(-1000*at*at)
                        for k,fun in enumerate(fs):vals[k]+=factor*fun(at)
                err=abs(vals[0]-vals[1]);tag=str(dps)+'_'+str(count)
                check('independent_fixed_cap_translation_'+tag,err<mp.mpf('1e-25'),mp.nstr(err,18))
                check('omitting_cap_jets_fails_'+tag,abs(vals[0]-vals[2])>1,mp.nstr(abs(vals[0]-vals[2]),18))
                if previous is not None:check('cap_quadrature_refinement_'+str(dps),max(abs(v-w) for v,w in zip(vals,previous))<mp.mpf('1e-25'))
                previous=vals
                records.append({'dps':dps,'nodes':count,'kernel_derivative':mp.nstr(vals[0],30),'complete_field_derivative':mp.nstr(vals[1],30),
                                'omitted_cap_derivative':mp.nstr(vals[2],30),'error':mp.nstr(err,18)})
    FIXTURES.append({'fixture':'independent_fixed_cap_translation','records':records,'cap_derivatives_required':True})

def independent_finite_translation():
    x,t=S.symbols('x t',real=True);p=S.Rational(3,1000)
    # At x=p the fractional term |x-p|^(19/2) is C9 and is not C10.
    vals=[S.diff(S.sin(17*x)/10000,x,j).subs(x,p) for j in range(10)]
    nonanalytic_coefficient=S.prod(S.Rational(19,2)-j for j in range(10))/10000
    check('manufactured_metric_C9_not_C10',nonanalytic_coefficient!=0 and S.Rational(19,2)-9>0 and S.Rational(19,2)-10<0)
    # A nonlinear observable of the value and fifth jet uses precisely metric jets through nine.
    q=[vals[j]+S.Rational(1,10**10)*vals[j+5] for j in range(5)]
    poly=sum(q[j]*t**j/S.factorial(j) for j in range(5))
    bell4=q[4]+4*q[1]*q[3]+3*q[2]**2+6*q[1]**2*q[2]+q[1]**4
    zero('finite_translation_fourth_jet_chain_rule',S.diff(S.exp(poly),t,4).subs(t,0)-S.exp(q[0])*bell4)
    check('nonanalytic_spatial_part_not_discarded',nonanalytic_coefficient>0)
    records=[]
    # Direct real finite differences include the nonanalytic term, independently of the polynomial.
    nodes=list(range(-4,5));weights={j:S.finite_diff_weights(j,nodes,0)[-1][-1] for j in [1,2,3,4]}
    for dps in PROTOCOL['numerical_dps']:
        with mp.workdps(dps):
            pp=mp.mpf('.003');coef5=mp.fprod(mp.mpf(19)/2-j for j in range(5))/10000
            def value(s):
                xx=pp+s
                f=mp.sin(17*xx)/10000+abs(s)**mp.mpf('9.5')/10000
                f5=17**5*mp.sin(17*xx+5*mp.pi/2)/10000+coef5*mp.sign(s)*abs(s)**mp.mpf('4.5')
                return mp.exp(f+mp.mpf('1e-10')*f5)
            for hh in ['.001','.0005']:
                h=mp.mpf(hh);errors=[]
                for j in [1,2,3,4]:
                    approx=sum(mf(F(str(w)))*value(k*h) for k,w in zip(nodes,weights[j]))/h**j
                    exact=mp.mpf(str(S.N(S.diff(S.exp(poly),t,j).subs(t,0),dps)))
                    error=abs(approx-exact);errors.append(mp.nstr(error,18))
                    check('real_C9_translation_'+str(dps)+'_'+hh+'_'+str(j),error<mp.mpf('1e-7'),mp.nstr(error,18))
                records.append({'dps':dps,'step':hh,'orders':[1,2,3,4],'errors':errors})
    FIXTURES.append({'fixture':'independent_real_C9_translation','records':records,'field':'sin(17x)/10000+abs(x-.003)^(19/2)/10000',
                     'observable':'exp(field+1e-10 fifth_spatial_derivative)','spatial_holomorphic_extension_used':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 finite middle massive source four jets','input_count':len(ps)}),flush=True)
    try:
        eta,step=retained_and_geometry();interpolation_and_cauchy();finite_rows_and_conversion()
        independent_dirac_energy();independent_cap_translation(eta,step);independent_finite_translation()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
    run={'story':'S167 finite middle massive source four jets','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
         'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after}
    for k,v in review().items():
        if k not in run and isinstance(v,bool):run[k]=v
    for kind,p in zip(['checker','note','review'],OWN):run['frozen_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
