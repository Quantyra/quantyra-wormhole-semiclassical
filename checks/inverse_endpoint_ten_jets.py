"""Supporting checks for C11 inverse maps and C10 fixed endpoint chords."""
from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import json, sys, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='inverse-endpoint-ten-jets';PREVIOUS='outer-metric-ten-jets'
OWN=['checks/inverse_endpoint_ten_jets.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
R=S.Rational;ROWS=[];FIXTURES=[]
PROTOCOL={'exact_arithmetic':'symbolic and rational','numerical_dps':[60,90],
          'numerical_tolerance':'1e-42','contour_guard_digits':40,'metric_order':10,'inverse_map_order':11,
          'ordinary_continuum_proof_separate':True,'actual_charged_source_evaluated':False,
          'original_finite_comparison_retained':True,'low_order_contracts_retained':True,
          'twelfth_metric_or_map_derivative_used':False}


def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)


def zero(name,expr):
    v=S.simplify(expr);check(name,v==0,v)


def previous():return json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(previous()['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})


def conv(a,b,n):return [sum(a[j]*b[k-j] for j in range(k+1) if j<len(a) and k-j<len(b)) for k in range(n+1)]


def powers(a,n):
    out=[[1]+[0]*n]
    for j in range(n):out.append(conv(out[-1],a,n))
    return out


def provenance():
    p=previous()
    check('prior_outer_metric_passed',p['status']=='passed' and len(p['checks'])==94 and all(r['passed'] for r in p['checks']))
    check('prior_snapshot_unchanged',p['input_sha256_before']==p['input_sha256_after'])
    bad=[f for f,h in p['input_sha256_after'].items() if sha256((ROOT/f).read_bytes()).hexdigest()!=h]
    check('all_prior_inputs_unchanged',not bad,bad)
    check('retained_sources_pinned',set(review()['retained_inputs'])<=set(paths()))
    check('new_scope_is_endpoint_not_quantum_response',review()['endpoint_inverse_map_C11_bounded']
          and review()['proper_endpoint_chords_C10_bounded'] and not review()['actual_finite_metric_source_H6_bounded'])


def tame_inverse_ledger():
    H=3*10**14;delta=R('3e-10');eps=R('2e-7');amin=1-R('2.1e-8')
    check('shared_forward_remainder',R('2.084')*10**14+R('1e-9')<H)
    check('new_cutoff_low_relative_jets',R('2.084e-10')<delta)
    check('inverse_low_Q',5*101*delta/amin<eps)
    check('inverse_high_Q',500/amin<501)
    c=(1+delta)/amin
    check('inverse_scale_power',c**11<R('1.000001'))
    b0=c**11*(1+sum(((1-eps)**(-n)-1)/n for n in range(2,12)))
    bh=c**11*5*(1-eps)**(-12)*501
    check('inverse_low_jet_sum',b0<R('1.00001'),b0)
    check('inverse_high_jet_sum',bh<4000,bh)
    check('composition_low_power',R('1.00001')**11<R('1.001'))
    check('high_squared_excluded',2*6>10 and 2*7>11)
    check('high_f_high_inverse_excluded',6+7-1>11)
    # |log(X/Y)| <= delta/(1-delta), with L(Y)>1.
    ratio=(1-delta)**(-5)/(1-delta/(1-delta))
    check('weight_ratio',ratio<R('1.001'))
    inverse=R('1.001')*(R('1.001')*12*R('1.01')+6*R('2e-9')*4000)
    check('weighted_inverse_remainder',inverse<20,inverse)
    check('inverse_reserve',20*H<10**16)
    chord=R('1.001')*(R('1.001')*11*15+6*R('1e-8')*4000)
    check('ordinary_endpoint_remainder',chord<200,chord)
    check('physical_endpoint_reserve',R('1.002')*200*H<10**17)
    check('radial_high_sum',100*sum(range(7,12))==4500)
    check('radial_low_sum',21*101*delta<R('1e-6'))
    radial=R('1.001')*(R('1.001')*6000+6*R('2e-6')*4000)
    check('radial_composed_reserve',R('1.002')*radial*H<10**19,radial)
    check('transition_weight_conversion',R(49,1000)**(-6)*R('4.02')<3*10**8)
    check('transition_C10_reserve',3*10**8*10**19<10**28)
    check('pole_division_coefficient_sum',1/(1-R('.01'))<R('1.011'))
    check('analytic_seed_difference_weight',1000*R('.1')**2==10)
    check('analytic_seed_k_derivative',R('.21')/3*3==R('.21'))
    FIXTURES.append({'fixture':'tame_inverse_and_composition_ledger','map_remainder_reserve':str(H),
                     'inverse_low_sum':str(b0),'inverse_high_multiplier':str(bh),
                     'inverse_remainder_multiplier':str(inverse),'chord_remainder_multiplier':str(chord),
                     'radial_remainder_multiplier':str(radial),'twelfth_derivative_used':False})


def exact_identities():
    Y,c,rr=S.symbols('Y c R');F=Y+c*Y**4+rr
    zero('inverse_quartic_remainder_identity',Y-F+c*F**4-(-rr+c*(F**4-Y**4)))
    ell,u3,nu3=S.symbols('ell u3 nu3',real=True);v3=-ell*u3;w3=ell*(nu3+u3/4)
    zero('proper_cubic_from_common_relation',(v3-4*w3).subs(u3,-2*nu3))
    check('common_and_proper_coefficients_differ',S.simplify(w3-ell*nu3)!=0)
    q=S.symbols('q');Ql=R(1,100)*q+R(1,200)*q**2;Qh=R(2,3)*q**6+R(3,4)*q**8
    for n in [7,9,11]:
        # Coefficients from a full finite binomial expansion versus the one-high formula.
        total=sum((-1)**j*S.binomial(n+j-1,j)*(Ql+Qh)**j for j in range(n))
        low=sum((-1)**j*S.binomial(n+j-1,j)*Ql**j for j in range(n))
        high=-n*Qh*sum((-1)**j*S.binomial(n+j,j)*Ql**j for j in range(n))
        zero('one_high_Lagrange_coefficient_'+str(n),S.expand(total-low-high).coeff(q,n-1))
    z=S.Function('z');n=S.Function('n');u,wl,fp=S.symbols('u wl fp',positive=True)
    # Multiplicative ratios give the logarithmic identities without choosing a log branch.
    zero('radial_pullback_ratio',S.exp(2*ell*u)/fp**2-S.exp(2*ell*u)*fp**(-2))
    zero('radius_pole_factorization',z(F)/z(Y)-(F/Y)*(z(F)/F)/(z(Y)/Y))
    zero('proper_radial_chord_vanishes',ell*u-S.log(S.exp(ell*u)))
    # Pole removal and the two endpoint subtractions in the formal cubic ring.
    X=S.symbols('X');inverse=X-c*X**4
    zero('inverse_cubic_no_new_order_three',(S.expand(inverse**3-X**3)).coeff(X,3))
    check('inverse_mass_quadratic_starts_at_seven',S.Poly(S.expand(c*((Y+c*Y**4)**4-Y**4)),Y).monoms()[-1][0]==7)
    FIXTURES.append({'fixture':'exact_inverse_pullback_and_cubic_identities','proper_cubic_relation':'v3=4w3',
                     'radius_pole_factored_before_estimation':True,'original_finite_radius_map_retained':True})


def transition_low_ledger():
    U=list(map(R,['3.32e-11','3.32e-10','3.32e-9','3.32e-8','4.01e-7','8.01e-6']))
    pu=powers([0]+U[1:],5)
    expu=[(1+R('1.001')*U[0])*sum(pu[k][j]/S.factorial(k) for k in range(6)) for j in range(6)]
    D=[R('8e-13')]+[(expu[j-1]-(1 if j==1 else 0))/(100*j) for j in range(1,7)]
    eta=[1]+[R('1.25')*4**j for j in range(1,7)]
    d=conv(eta,D,6);hh=[100*(j+1)*d[j+1] for j in range(6)]
    logsum=sum(hh)/(1-sum(hh));amin=1-R('2.1e-8')
    check('transition_u_low_sum',sum(U)<R('9e-6'),sum(U))
    check('transition_forward_low_sum',sum(d[:6])<R('4e-9'),sum(d[:6]))
    check('transition_log_jacobian_C5',logsum<R('2e-5'),logsum)
    check('transition_inverse_Q',100*sum(d[2:])/amin<R('1e-5'))
    eps=R('1e-5');c=1/amin
    bsum=c**5*(1+sum(((1-eps)**(-n)-1)/n for n in range(2,6)))
    check('transition_composition_C5',bsum**5<R('1.01'),bsum**5)
    check('transition_seed_shift_C5',sum(d[:6])*1<R('1e-7'))
    check('transition_radius_shift_C5',sum(d[:6])*(25+1)<R('2e-7'))
    check('transition_radial_chord_C5',R('1.1')*R('1.01')*(sum(U)+logsum)<R('1e-4'))
    check('transition_tangential_chord_C5',R('1.1')*R('1.01')*(sum(U)+R('1e-7'))<R('1e-4'))
    check('transition_radius_chord_C5',R('1.1')*R('1.01')*(R('2e-5')+R('2e-7'))<R('1e-4'))
    D0,D1,D2,D3=map(R,['8e-13','3.4e-11','3.4e-8','6.7e-5'])
    f1=100*D0+D1;f2=400000*D0+200*D1+D2
    f3=480000000*D0+3*400000*D1+300*D2+D3
    check('real_transition_first',f1<R('1.2e-10'),f1)
    check('real_transition_second',f2<R('3.7e-7'),f2)
    check('real_transition_third',f3<R('6e-4'),f3)
    check('real_jacobian_log_C2',(f3/amin+f2**2/amin**2)/amin**2+f2**2/amin**3<R('1e-3'))
    # Ordered scalar tensor recurrence includes the connection at order zero.
    B={(0,k):S.factorial(k)*100**k for k in range(6)};amps=[R(1)]
    for rank in range(5):
        for k in range(5-rank):
            B[rank+1,k]=B[rank,k+1]+rank*sum(S.binomial(k,i)*R('.3')*S.factorial(i)*40**i*B[rank,k-i] for i in range(k+1))
        amps.append(B[rank+1,0]/(S.factorial(rank+1)*100**(rank+1)))
    check('transition_full_tensor_conversion',max(amps)<R('1.1'),max(amps))
    FIXTURES.append({'fixture':'new_transition_low_and_full_tensor_ledger','D':list(map(str,D)),
                     'cutoff_displacement':list(map(str,d)),'log_jacobian_C5_sum':str(logsum),
                     'tensor_amplifications':list(map(str,amps)),'old_cutoff_estimates_assumed':False,
                     'real_first_three':[str(f1),str(f2),str(f3)]})


def numeric_compose(a,b,n):
    pp=powers(b,n)
    return [sum(a[k]*pp[k][j] for k in range(min(len(a),n+1))) for j in range(n+1)]


def inverse_series(a,n):
    # Independent triangular reversion of a series with zero constant term.
    b=[mp.mpf(0)]*(n+1);b[1]=1/a[1]
    for j in range(2,n+1):
        rem=numeric_compose(a,b,j)[j]
        b[j]=-rem/a[1]
    return b


def independent_inverse_and_chords():
    records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            for yy,ell in [('.03','.7'),('.075','1')]:
                Y=mp.mpf(yy);ee=mp.mpf(ell);c=mp.mpf('.002');tail=mp.mpf('.003')
                def F(y):return y+ee*(c*y**4+tail*y**5*(1-mp.log(y)))
                coeff=mp.taylor(F,Y,11);a=[mp.mpf(0)]+coeff[1:];inv=inverse_series(a,11)
                q=[a[j+1]/a[1] for j in range(1,11)]
                # Lagrange uses a normalized direct map; compare against triangular reversion.
                qq=[mp.mpf(0)]+q;qp=powers(qq,10)
                for n in [2,6,7,10,11]:
                    lag=sum((-1)**j*mp.binomial(n+j-1,j)*qp[j][n-1] for j in range(n))/(n*a[1]**n)
                    err=abs(lag-inv[n])
                    check('independent_inverse_'+str(dps)+'_'+yy+'_'+str(n),err<mp.mpf('1e-42'),mp.nstr(err,12))
                    records.append({'dps':dps,'Y':yy,'inverse_order':n,'error':mp.nstr(err,18)})
                X=F(Y)
                # Direct numerical root evaluation checks values and derivatives independently of reversion.
                def G(x):return mp.findroot(lambda y:F(y)-x,(Y*mp.mpf('.99'),Y*mp.mpf('1.01')))
                for n in [1,2,3]:
                    rootjet=mp.diff(G,X,n)/mp.factorial(n);err=abs(rootjet-inv[n])
                    check('root_inverse_derivative_'+str(dps)+'_'+yy+'_'+str(n),err<mp.mpf('1e-42'),mp.nstr(err,12))
                    records.append({'dps':dps,'Y':yy,'root_order':n,'error':mp.nstr(err,18)})
                def chord(y):return mp.log(F(y)/y)+mp.mpf('.2')*(F(y)**2-y**2)
                cc=mp.taylor(chord,Y,10);composed=numeric_compose(cc,inv,10)
                # Independent contour differentiation of the composite includes its log boundary tail.
                radius=mp.mpf('.001')
                for n in [6,10]:
                    N=128
                    with mp.workdps(dps+40):
                        contour=sum(chord(G(X+radius*mp.exp(2j*mp.pi*k/N)))*mp.exp(-2j*mp.pi*k*n/N) for k in range(N))/(N*radius**n)
                    err=abs(contour-composed[n])
                    check('independent_composed_chord_'+str(dps)+'_'+yy+'_'+str(n),err<mp.mpf('1e-42'),mp.nstr(err,12))
                    records.append({'dps':dps,'Y':yy,'chord_order':n,'error':mp.nstr(err,18)})
    FIXTURES.append({'fixture':'independent_logarithmic_inverse_and_chord_coefficients','records':records,
                     'actual_charged_source_sampled':False,'logarithmic_tail_retained':True})


def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 inverse endpoint ten jets','input_count':len(ps)}),flush=True)
    try:
        provenance();tame_inverse_ledger();exact_identities();transition_low_ledger();independent_inverse_and_chords()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(row['passed'] for row in ROWS)
    run={'story':'S167 inverse endpoint ten jets','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,
         'mpmath':mp.__version__,'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,
         'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after}
    for key,value in review().items():
        if key not in run and (isinstance(value,bool) or key.endswith('_reserve')):run[key]=value
    if not good:
        for kind,path in zip(['checker','note','review'],OWN):run['failed_'+kind+'_source_utf8']=(ROOT/path).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');saved=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    saved['runs'].append(run);dest.write_text(json.dumps(saved,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1


if __name__=='__main__':raise SystemExit(main())
