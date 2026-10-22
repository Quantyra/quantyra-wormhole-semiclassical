"""Actual physical York adjoint, density conversion and complete shape split."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
from math import factorial,comb
from itertools import product
import importlib.util,json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-higher-york-source';PREVIOUS='finite-refreshed-local-mixed'
OWN=['checks/finite_higher_york_source.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
spec=importlib.util.spec_from_file_location('retained_physical_lift',ROOT/'checks/finite_higher_tt_lift.py')
G=importlib.util.module_from_spec(spec);spec.loader.exec_module(G);L=G.L
ROWS=[];FIXTURES=[]
PROTOCOL={'output_order':4,'density_order':4,'known_scalar_order':3,'shape_source_order':4,
          'numerical_dps':[60,90],'numerical_tolerance':'1e-42','integral_derivative_orders':[0,1,4],
          'actual_charged_source_sampled':False,'ordinary_continuum_proof_separate':True}

def check(name,ok,detail=''):
    assert name not in {x['name'] for x in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)
def zero(name,v):
    v=S.simplify(v);check(name,v==0,v)
def saved(stem):return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def fixture(stem,name):return next(x for x in saved(stem)['fixtures'] if x['fixture']==name)
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})
def strings(a):return [str(x) for x in a]
def decimals(a):return [L.decimal(x) for x in a]
def va(*args):return [sum(a[j] for a in args) for j in range(10)]
def vs(k,a):return [F(k)*v for v in a]
def vz():return [F(0)]*10
def conv(a,b,n):return va(*[vs(comb(n,j)*a[j],b[n-j]) for j in range(n+1)])

def inputs_and_low():
    old=saved(PREVIOUS)
    check('previous_complete_mixed_passed',old['status']=='passed' and len(old['checks'])==99 and all(x['passed'] for x in old['checks']))
    check('all_previous_inputs_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    check('prior_missing_source_and_inverse_retained',not old['actual_TT_source_H4_bounded'] and not old['higher_TT_inverse_quantified'] and not old['actual_transfer_certificate_passed'])
    check('actual_return_inverse_guard',F(2)/(1-F('.005'))**2<F('1.422')**2)
    check('radial_radius_return_guard',F('2.25')*F('1.422')<F('3.2'))
    check('adjoint_v_low_guard',F('1.008')+F('3.2')<F('4.21'))
    check('adjoint_Z_low_guards',F('1.422')*F('3.2')<F('4.6') and F('1.422')*F('4.71')<F('6.7'))
    p=1+2*(F('1.008')+F('3.2'))
    check('adjoint_Psi_low_guards',p*F('3.2')<31 and p*F('4.71')<45)
    check('adjoint_Upsilon_low_guard',45+F(1,2)<46)
    check('physical_L2_transpose_guard',F('4.6')*F('4.7')<22)
    check('intrinsic_kappa_value_guard',F('1.42')**2>2 and F('1.42')*15+2*F('5.1')<32)
    check('density_value_guard',1/(1-6*F('3e-8'))<F('1.001'))
    check('retained_scalar_L2_guard',27*F('2.2e-13')<F('6e-12'))

def exact_adjoint_identities():
    x=S.symbols('rho',positive=True)
    aa,rr,v,z,s,t=[S.Function(n)(x) for n in ['A','R','v','Z','s','t']]
    w=v-t/2;gg=s/2-t;psi=w+2*(aa+rr)*z;up=psi+t/2
    sub={S.diff(v,x):rr*s+(aa-rr)*t-2*(aa+rr)*v,S.diff(z,x):aa*z+w}
    kap=2*(S.diff(aa,x)+aa**2+S.diff(rr,x)+aa*rr)
    zero('adjoint_w_equation',S.diff(w,x).subs(sub)+2*(aa+rr)*w-(2*rr*gg-S.diff(t,x)/2))
    zero('intrinsic_Psi_equation',S.diff(psi,x).subs(sub)-(kap*z+2*rr*gg-S.diff(t,x)/2))
    zero('intrinsic_Upsilon_equation',S.diff(up,x).subs(sub)-(kap*z+2*rr*gg))
    zero('intrinsic_adjoint_divergence',S.diff(z,x).subs(sub)+(aa+2*rr)*z-(up-t/2))
    direct=(s-3*t)/2+v+2*(aa-rr)*z
    zero('intrinsic_adjoint_output',direct-(s/2-3*t/2+up-4*rr*z))
    bs=S.symbols('b_src');zero('scalar_principal_term_cancels',(s/2-3*t/2).subs(t,s/3+2*bs)+3*bs)
    a=S.Function('a')(x);r=S.Function('r')(x);mu=a*r*r;uu=a*a*r*r
    zero('physical_measure_H_adjoint_denominator',a*mu-uu)
    zero('physical_measure_B_adjoint_coefficient',uu/mu-a)
    jj=1+x*x;zz=x*S.exp(-x*x)
    zero('density_return_commutator',S.diff(jj*zz,x)-aa*jj*zz-jj*(S.diff(zz,x)-aa*zz)-S.diff(jj,x)*zz)
    check('density_does_not_commute_with_return',S.diff(jj,x)*zz!=0)
    FIXTURES.append({'fixture':'exact_physical_adjoint_system','physical_radial_measure':'a r^2 d rho',
                     'kappa':'-tr_P Ricbar-2|gradbar log r|^2','scalar_principal_term':'-3 b_src',
                     'density_before_transpose':True,'singular_A_bound_used':False,'R_prime_seed_identity_used':False})

def independent_geometry():
    x,t,th,ph=S.symbols('rho time theta phi',positive=True);coords=[x,t,th,ph]
    a=S.Function('a')(x);r=S.Function('r')(x);g=S.diag(1,a*a,r*r,r*r*S.sin(th)**2);inv=g.inv()
    conn={(i,j,k):S.simplify(sum(inv[i,l]*(S.diff(g[l,k],coords[j])+S.diff(g[j,l],coords[k])-S.diff(g[j,k],coords[l]))/2 for l in range(4))) for i,j,k in product(range(4),repeat=3)}
    def ric(b,d):return S.simplify(sum(S.diff(conn[k,d,b],coords[k])-S.diff(conn[k,k,b],coords[d])+sum(conn[k,k,l]*conn[l,d,b]-conn[k,d,l]*conn[l,k,b] for l in range(4)) for k in range(4)))
    aa=S.diff(a,x)/a;rr=S.diff(r,x)/r;rric=ric(0,0);tric=ric(1,1)/(a*a)
    zero('independent_four_dimensional_orbit_Ricci_trace',-rric-tric-2*rr*rr-2*(S.diff(aa,x)+aa*aa+S.diff(rr,x)+aa*rr))
    zz=S.Function('Z')(x);div=S.diff(zz,x)+sum(conn[i,i,0]*zz for i in range(4))
    zero('independent_four_dimensional_vector_divergence',div-S.diff(zz,x)-(aa+2*rr)*zz)
    FIXTURES.append({'fixture':'independent_warped_adjoint_geometry','dimension':4,'orbit_Ricci_trace':str(S.simplify(rric+tric)),
                     'actual_warped_metric_used':True,'seed_equation_imposed':False})

def coefficient_arrays():
    f=fixture('finite-higher-tt-lift','full_finite_lift_coefficient_arrays')
    q=[F(x) for x in f['q'][:5]];rm=[F(x) for x in f['Rm'][:4]];ric=[F(x) for x in f['Ric'][:4]];pp=[F(x) for x in f['projector'][:4]]
    kap=[sum(F(comb(n,j))*(pp[j]*ric[n-j]+2*q[j]*q[n-j]) for j in range(n+1)) for n in range(4)];kap[0]=min(kap[0],32)
    check('full_coefficient_orders_available',len(q)==5 and len(kap)==len(rm)==4 and min(q+kap+rm)>0)
    for j in range(4):
        check('retained_Gaffney_curvature_'+str(j),16*(j+1)*rm[0]<(32*(j+1))**2)
    FIXTURES.append({'fixture':'full_adjoint_coefficient_arrays','q':strings(q),'Rm':strings(rm),'kappa':strings(kap),
                     'kappa_decimals':decimals(kap),'projector_and_Ricci_derivatives_retained':True})
    return q,rm,kap

def recurrence(q,rm,kap,shape):
    unit=[[F(int(j==k)) for k in range(10)] for j in range(10)];s=unit[:5];t=unit[5:]
    if shape:
        z=[va(vs(F(41,6),s[0]),vs('13.4',t[0]))];u=[va(vs(F(139,3),s[0]),vs(92,t[0]))]
        gg=[va(vs(F(1,6),s[j]),vs(2,t[j])) for j in range(5)]
        extra=[va(vs(F(1,6),s[j]),t[j]) for j in range(5)]
        local=[vs(3,t[j]) for j in range(5)];low=va(vs(F(88,3),s[0]),vs(44,t[0]))
    else:
        z=[va(vs('4.6',s[0]),vs('6.7',t[0]))];u=[va(vs(31,s[0]),vs(46,t[0]))]
        gg=[va(vs(F(1,2),s[j]),t[j]) for j in range(5)];extra=[vs(F(1,2),t[j]) for j in range(5)]
        local=[va(vs(F(1,2),s[j]),vs(F(3,2),t[j])) for j in range(5)];low=va(vs(22,s[0]),vs(22,t[0]))
    for j in range(4):
        curve=va(*[vs(16*(j+1)*comb(j,i+1)*rm[i],z[j-1-i]) for i in range(j)]) if j else vz()
        z.append(va(vs(2,u[j]),vs(2,extra[j]),vs(32*(j+1),z[j]),curve))
        u.append(va(conv(kap,z,j),vs(2,conv(q,gg,j))))
    out=[va(local[j],u[j],vs(4,conv(q,z,j))) for j in range(5)];out[0]=low
    columns=va(*out)
    return {'Z':z,'Upsilon':u,'output':out,'columns':columns}

def conversion_matrices(q,rm,kap):
    disk=recurrence(q,rm,kap,False);shape=recurrence(q,rm,kap,True)
    for name,data in [('disk',disk),('shape',shape)]:
        check(name+'_matrix_nonnegative',all(v>=0 for key in ['Z','Upsilon','output'] for row in data[key] for v in row))
        check(name+'_output_order_triangular',all(data['output'][j][block*5+k]==0 for j in range(5) for block in range(2) for k in range(j+1,5)))
        check(name+'_auxiliary_gain_one_derivative',all(data[key][j][block*5+k]==0 for key in ['Z','Upsilon'] for j in range(1,5) for block in range(2) for k in range(j,5)))
        check(name+'_column_norm_reconciles',all(sum(row[k] for row in data['output'])==data['columns'][k] for k in range(10)))
    check('scalar_fourth_order_not_required',shape['columns'][4]==0 and all(shape['output'][j][j]==0 for j in range(1,5)))
    check('shape_fourth_order_still_required',shape['columns'][9]==3)
    factor=max(disk['columns']);ep=L.decimal_power(factor)
    check('upward_physical_York_H4_enclosure',F(10)**(ep-1)<factor<=F(10)**ep)
    FIXTURES.append({'fixture':'physical_York_H4_matrices','disk':{k:[strings(row) for row in v] if k!='columns' else strings(v) for k,v in disk.items()},
                     'shape':{k:[strings(row) for row in v] if k!='columns' else strings(v) for k,v in shape.items()},
                     'disk_column_decimals':decimals(disk['columns']),'shape_column_decimals':decimals(shape['columns']),
                     'C_Y4_exact':str(factor),'C_Y4_decimal':L.decimal(factor),'C_Y4_power_of_ten':ep,
                     'scalar_order_four_column_zero':True,'shape_order_four_retained':True})
    return disk,shape

def density_and_budget(disk,shape):
    ff=fixture(PREVIOUS,'integrated_finite_field_bounds');gf=fixture(PREVIOUS,'refreshed_finite_component_geometry')
    tp=L.poly([F(x) for x in ff['Tplus']]);cc=L.poly([F(x) for x in gf['geometry']['C']])
    dj=L.covariant(L.scale('1.001',L.exps(L.scale(6,tp))),cc,0,4);dj[0]=min(dj[0],F('1.001'))
    den=[[F(comb(n,l))*dj[n-l] if l<=n else F(0) for l in range(5)] for n in range(5)]
    full=[[sum(disk['output'][n][block*5+k]*den[k][l] for k in range(5)) for block in range(2) for l in range(5)] for n in range(5)]
    cols=va(*full);bound=max(cols);bp=L.decimal_power(bound)
    check('density_product_raw_divided_reconciliation',all(den[n][l]==F(factorial(n),factorial(l)*factorial(n-l))*dj[n-l] for n in range(5) for l in range(n+1)))
    check('density_then_adjoint_full_matrix',all(full[n][block*5+l]==sum(disk['output'][n][block*5+k]*comb(k,l)*dj[k-l] for k in range(l,5)) for n in range(5) for block in range(2) for l in range(5)))
    check('upward_density_York_H4_enclosure',F(10)**(bp-1)<bound<=F(10)**bp)
    sf=fixture('finite-radius-operator-source','complete_refined_scalar_and_anomaly')
    sig=[F(1001,24000)*F(x) for x in sf['physical_covariant_bounds'][:5]];sig[0]=min(sig[0],F('6e-12'))
    scalar_parts=[shape['columns'][j]*sig[j] for j in range(5)];cs=sum(scalar_parts);cb=max(shape['columns'][5:])
    sp=L.decimal_power(cs);cp=L.decimal_power(cb)
    check('known_scalar_H3_allocation',scalar_parts[4]==0 and cs>0 and cb>0)
    check('upward_known_scalar_conversion',F(10)**(sp-1)<cs<=F(10)**sp)
    check('upward_complete_shape_conversion',F(10)**(cp-1)<cb<=F(10)**cp)
    iv=fixture(PREVIOUS,'low_TT_inverse_composition')
    check('full_source_budget_unchanged',iv['Cstart_power_of_ten']==29 and iv['mhi_power_of_ten']==13 and not iv['target_passed'])
    high=F(iv['high_R_C_D4_F4_budget_over_g_epsilon_c']);bs=high/(2*F(10)**sp);bb=high/(2*F(10)**cp)
    check('two_subbudgets_sum_to_retained_high_budget',F(10)**sp*bs+F(10)**cp*bb==high)
    FIXTURES.append({'fixture':'density_conversion_and_complete_shape_budget','density_jets':strings(dj),'density_product_matrix':[strings(row) for row in den],
                     'fixed_density_output_matrix':[strings(row) for row in full],'fixed_density_columns':strings(cols),
                     'C_density_Y4_exact':str(bound),'C_density_Y4_decimal':L.decimal(bound),'C_density_Y4_power_of_ten':bp,
                     'known_scalar_derivative_bounds':strings(sig),'known_scalar_parts':strings(scalar_parts),'known_scalar_part_decimals':decimals(scalar_parts),
                     'C_sigma_exact':str(cs),'C_sigma_decimal':L.decimal(cs),'C_sigma_power_of_ten':sp,
                     'C_b4_exact':str(cb),'C_b4_decimal':L.decimal(cb),'C_b4_power_of_ten':cp,
                     'R_CD4_budget_over_g_squared_epsilon':str(bs),'R_CD4_Fb4_budget_over_g_epsilon':str(bb),
                     'R_CD4_budget_decimal':L.decimal(bs),'R_CD4_Fb4_budget_decimal':L.decimal(bb),
                     'Cstart_exact':iv['Cstart_exact'],'mhi_exact':iv['mhi_exact'],'actual_complete_shape_H4_supplied':False,
                     'actual_TT_H4_supplied':False,'higher_inverse_supplied':False,'target_passed':False})

def manufactured(alpha,beta):
    x=S.symbols('rho',real=True);aa=1/x+2*alpha*x;rr=2*beta*x;k=S.Integer(3)
    zp=x*(1+x*x+2*x**4+x**6)/10**6;tp=(2+x*x+x**4+x**8)/10**6
    dd=lambda p:S.diff(p,x)-2*k*x*p
    wp=S.expand(dd(zp)-aa*zp)
    sp=S.cancel(2*tp+(dd(wp)+2*(aa+rr)*wp+dd(tp)/2)/rr)
    vp=S.expand(wp+tp/2);psip=S.expand(dd(zp)+(aa+2*rr)*zp);up=S.expand(psip+tp/2)
    kap=2*(S.diff(aa,x)+aa*aa+S.diff(rr,x)+aa*rr)
    fc=S.expand(sp/2-3*tp/2+up-4*rr*zp)
    tag=str(alpha)+'_'+str(beta)
    zero('manufactured_v_equation_'+tag,dd(vp)+2*(aa+rr)*vp-rr*sp-(aa-rr)*tp)
    zero('manufactured_Upsilon_equation_'+tag,dd(up)-kap*zp-2*rr*(sp/2-tp))
    check('manufactured_regular_bolt_'+tag,wp.subs(x,0)==0 and zp.subs(x,0)==0 and all(S.Poly(p,x).degree()>=0 for p in [sp,tp,fc]))
    check('manufactured_not_seed_relation_'+tag,S.simplify(S.diff(rr,x)-aa*rr+rr*rr)!=0)
    check('omitted_adjoint_curvature_detected_'+tag,S.expand(kap*zp).coeff(x,1)!=0)
    def coefficients(p):
        p=S.Poly(S.expand(p),x);return [F(str(p.nth(j))) for j in range(p.degree()+1)]
    return {name:coefficients(p) for name,p in [('s',sp),('t',tp),('v',vp),('w',wp),('Z',zp),('Psi',psip),('Upsilon',up),('fc',fc)]}

def exact_pairings(configs,fields):
    x=S.symbols('rho',real=True);G.check=check;records=[]
    for af,bf in configs:
        aa=S.Rational(af);bb=S.Rational(bf);forward=G.manufactured_polynomials(aa,bb)
        def expr(coeff):return sum(S.Rational(v.numerator,v.denominator)*x**j for j,v in enumerate(coeff))
        c,b,y=[expr(forward[key]) for key in ['c','b','Y']]
        s,t,fc=[expr(fields[af,bf][key]) for key in ['s','t','fc']]
        rr=2*bb*x;ar=1/x+2*(aa-bb)*x
        sc=S.expand(c/2+rr*y);tc=S.expand(-(b+3*c)/2+ar*y)
        left=S.expand(s*sc+t*tc);right=S.expand(c*fc);difference=S.Poly(left-right,x)
        rate=5-aa-2*bb
        def physical(p):
            p=S.Poly(p,x)
            assert all(p.nth(j)==0 for j in range(1,p.degree()+1,2))
            return sum(p.nth(2*j)*S.factorial(j)/(2*rate**(j+1)) for j in range(p.degree()//2+1))
        lv=physical(left);rv=physical(right)
        zero('independent_physical_adjoint_pairing_'+str(af)+'_'+str(bf),lv-rv)
        flat=sum(difference.nth(j)*S.gamma(S.Rational(j+1,2))/(2*S.Integer(5)**S.Rational(j+1,2)) for j in range(difference.degree()+1))
        flat=S.simplify(flat)
        check('flat_density_adjoint_pairing_fails_'+str(af)+'_'+str(bf),flat!=0)
        records.append({'alpha':str(af),'beta':str(bf),'physical_left':str(lv),'physical_right':str(rv),
                        'physical_radial_density':'rho exp((alpha+2 beta)rho^2)',
                        'flat_density_discrepancy':str(flat),'forward_decay_rate':2,'source_decay_rate':3})
    FIXTURES.append({'fixture':'independent_physical_adjoint_pairings','records':records,
                     'bolt_and_outer_boundary_retained':True,'actual_charged_geometry_sampled':False})

def mfrac(v):return mp.mpf(v.numerator)/v.denominator
def evalp(a,x):return mp.polyval([mfrac(v) for v in reversed(a)],x)

def independent_integrals():
    configs=[(F(1,10),F(1,5)),(F(7,20),F(1,7))]
    fields={key:manufactured(S.Rational(key[0]),S.Rational(key[1])) for key in configs}
    exact_pairings(configs,fields);records=[]
    for dps in [60,90]:
        with mp.workdps(dps):
            for alpha,beta in configs:
                af,bf=mfrac(alpha),mfrac(beta);p=fields[alpha,beta]
                field=lambda name,x:mp.exp(-3*x*x)*evalp(p[name],x)
                def vi(x):
                    return mp.exp(-2*(af+bf)*x*x)*mp.quad(lambda y:mp.exp(2*(af+bf)*x*x*y*y)*(2*bf*x*x*y**3*field('s',x*y)+(y+2*(af-bf)*x*x*y**3)*field('t',x*y)),[0,1])
                def zj(x):
                    rate=3+af
                    return -mp.exp(af*x*x)*sum(mfrac(p['w'][2*j])*mp.gammainc(j,rate*x*x,mp.inf)/(2*rate**j) for j in range(1,len(p['w'])//2+1))
                def fci(x):return field('s',x)/2-3*field('t',x)/2+vi(x)+2*(1+2*(af-bf)*x*x)*zj(x)
                for at in [mp.mpf(0),mp.mpf('.3')]:
                    errors=[];tag=str(alpha)+'_'+str(beta)+'_'+str(at)+'_'+str(dps)
                    for n in PROTOCOL['integral_derivative_orders']:
                        error=max(abs(mp.diff(vi,at,n)-mp.diff(lambda x:field('v',x),at,n)),
                                  abs(mp.diff(lambda x:x*zj(x),at,n)-mp.diff(lambda x:field('Z',x),at,n)),
                                  abs(mp.diff(fci,at,n)-mp.diff(lambda x:field('fc',x),at,n)))/mp.factorial(n)
                        check('independent_adjoint_integrals_'+tag+'_'+str(n),error<mp.mpf('1e-42'),mp.nstr(error,18));errors.append(mp.nstr(error,18))
                    records.append({'alpha':str(alpha),'beta':str(beta),'rho':str(at),'dps':dps,
                                    'orders':PROTOCOL['integral_derivative_orders'],'errors':errors})
    FIXTURES.append({'fixture':'independent_fixed_adjoint_integrals','records':records,'exact_bolt_included':True,
                     'bolt_fixed_v_and_outer_fixed_Z':True,'all_fields_decay_at_infinity':True,'actual_charged_source_sampled':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 higher physical York source conversion','input_count':len(ps)}),flush=True)
    try:
        inputs_and_low();exact_adjoint_identities();independent_geometry();q,rm,kap=coefficient_arrays()
        disk,shape=conversion_matrices(q,rm,kap);density_and_budget(disk,shape);independent_integrals()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(x['passed'] for x in ROWS)
    run={'story':'S167 higher physical York source conversion','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
         'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after}
    for key,value in review().items():
        if key not in run and isinstance(value,bool):run[key]=value
    for kind,p in zip(['checker','note','review'],OWN):run['frozen_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
