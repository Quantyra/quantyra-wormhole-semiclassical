"""Supporting finite-jet, source and global-norm checks for S167 outer H4."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
from math import factorial,comb
import json,sys,traceback
import sympy as S
from sympy.functions.combinatorial.numbers import stirling
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='finite-outer-shape-source';PREVIOUS='finite-interior-shape-source'
OWN=['checks/finite_outer_shape_source.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
PROTOCOL={'source_order':4,'coefficient_order':10,'dilation_radius':'1e-40',
          'physical_dimension':4,'outer_seed_cut':40,'digits':[60,90],
          'quadrature_nodes':[32,64],'tolerance':'1e-25','actual_quantum_state_sampled':False}

def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)}
    ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)

def saved(stem):return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def fixture(stem,name):return next(x for x in saved(stem)['fixtures'] if x['fixture']==name)
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})
def sf(v):return F(str(v))
def mf(v):return mp.mpf(v.numerator)/v.denominator if isinstance(v,F) else mp.mpf(v)
def dec(v):
    with mp.workdps(30):return mp.nstr(mf(v),20)
def st2(j,k):return int(stirling(j,k,kind=2))
def st1(j,k):return int(stirling(j,k,kind=1,signed=True))
def power10(v):
    n=0
    while F(10)**n<v:n+=1
    return n

def inputs_and_dilation():
    old=saved(PREVIOUS)
    check('interior_receipt_passed',old['status']=='passed' and len(old['checks'])==170 and all(r['passed'] for r in old['checks']))
    check('all_previous_inputs_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    for stem in ['inverse-endpoint-ten-jets','outer-finite-background-response','seed-outer-eight-jets','outer-metric-ten-jets']:
        r=saved(stem)
        check(stem+'_receipt',r['status']=='passed' and all(q['passed'] for q in r['checks']))
        check(stem+'_inputs',r['input_sha256_before']==r['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in r['input_sha256_after'].items()))
    check('outer_coverage',F(3,40)<F('.076')/(1+F('3e-10')))
    check('cap_inside_proper_coordinate_domain',F('.099')<F('.1')*(1-F('3e-10')))
    check('commutator_separation',F('.084')-F('.076')>F('.004'))
    cut=sum(F('4.2')**j/factorial(j) for j in range(11))
    check('cutoff_normalized_jet_sum',cut<67,dec(cut))
    check('capped_low_weighted_contract',F('1e-7')*cut+F('1e-11')*12*2*cut+F('1e-9')<F('1e-5'))
    check('capped_high_weighted_contract',F('1e17')*cut+F('1e10')<F('1e20'))
    K=[]
    for j in range(6):
        row=[]
        for a in range(6):
            v=F('1e20')*sum(F(st2(j,k)*comb(a,r)*factorial(k)*factorial(a+k-r)*100**(k-r),factorial(k-r)*factorial(a)) for k in range(j+1) for r in range(min(a,k)+1))
            row.append(v)
        K.append(max(row))
    rho=F('1e-40');tube=F('1e-5')+sum(rho**j*K[j]/factorial(j) for j in range(1,5))
    check('finite_dilation_cap_tube',tube<F('2e-5'),dec(tube))
    check('fifth_dilation_derivative_uses_ten',5+5==10 and K[5]<F('1e36'),dec(K[5]))
    seed=[]
    for j in range(1,5):
        vals=[]
        for a in range(4):
            vals.append(sum(F(st2(j,k)*comb(a,r)*factorial(k),factorial(k-r))*100**(k-r)*2*factorial(a+k-r)*100**(a+k-r) for k in range(j+1) for r in range(min(a,k)+1)))
        seed.append(max(vals))
    seed_change=sum(rho**j*seed[j-1]/factorial(j) for j in range(1,5))
    check('global_seed_finite_dilation_small',seed_change*F('1e6')<F('1e-20'),dec(seed_change))
    check('longer_seed_cap_heat_budget',F('2.4e24')+F('1.01e24')+F('1.536e22')+F('1e21')<F('5e24'))
    check('cap_scaled_angular_budget',F('.9998')**-6<F('1.002'))
    check('complex_cap_source_budget',2*2*F('8.5e-33')*F('1e25')*F('2e-5')<F('1e-11'))
    y=S.symbols('y',positive=True);b=S.Function('b')(y)
    for j in range(5):
        direct=b
        for _ in range(j):direct=y*S.diff(direct,y)
        for a in range(6):
            rhs=sum(st2(j,k)*comb(a,r)*factorial(k)//factorial(k-r)*y**(k-r)*S.diff(b,y,a+k-r) for k in range(j+1) for r in range(min(a,k)+1))
            check('Euler_ordinary_jet_identity_'+str(j)+'_'+str(a),S.expand(S.diff(direct,y,a)-rhs)==0)
    FIXTURES.append({'fixture':'finite_dilation_and_cap_contract','K':[str(v) for v in K],
                     'rho':str(rho),'tube_bound':str(tube),'seed_relative_change':str(seed_change),
                     'actual_metric_spatial_analyticity_assumed':False,'positive_dilation_defect_omitted':False,
                     'whole_cap_four_jets_match':True,'coefficient_order':10})

def gluing_and_source():
    check('evaluation_ball_fits',F(1,100000)<F('.004')/4)
    check('mass_plus_scale',1+F(10000,18)<600)
    check('point_rows_source_cutoff',2*(F('3e5')**2)*F('1e4')<F('1e16'))
    check('normal_recovery_margin',sum(F(2*factorial(6),factorial(6-r))*(F(1000)/(F('1e8')*18))**r for r in range(1,4))<F('.001'))
    density=F('1e10')
    check('fourth_density_log_budget',2*factorial(4)*100**4+F('1e-7')*factorial(4)*100**4<F('1e10'))
    check('endpoint_four_derivative_budget',2**8*density*5*factorial(4)*F('1e8')**4<F('1e50'))
    check('per_mode_all_four_bound',F('1e16')*F('1e50')*600**6<F('1e85'))
    moment=4*(factorial(7)+5*factorial(6))*1000**7
    defect=F('1e85')*moment*2*F('8.5e-33')*F('1e4')
    check('actual_cap_endpoint_difference',defect<F('1e85'),dec(defect))
    check('actual_domain_graph_cone',10*F('1e-4')<=F('.001'))
    euler=[factorial(j)*F('1e-11')*F('1e40')**j+F('1e85') for j in range(5)]
    ordinary=[sum(abs(st1(j,k))*euler[k] for k in range(j+1)) for j in range(5)]
    check('complete_regular_massive_four_jets',max(ordinary)<F('1e154'),dec(max(ordinary)))
    check('cubic_restore_four_jets',factorial(4)*2*F('8.5e-33')*F('1e-11')/15<F('1e-42'))
    check('forward_map_four_jet_composition',1+sum(F('3e-10')*factorial(j+1)*100**(j+1) for j in range(1,4))<10)
    check('proper_to_common_reserve',1000*F('1e155')*F('1.001')+F('1e80')<F('1e160'))
    product_bound=sum(comb(4,j)*factorial(j)*1000**j for j in range(5))*F('1e160')+F('1e80')
    check('complete_massive_radius_product',product_bound<F('1e176'),dec(product_bound))
    check('massive_bookkeeping_residual',F('6.1')*product_bound<F('1e177'))
    FIXTURES.append({'fixture':'actual_short_cap_and_massive_four_jets','actual_cap_difference':str(defect),
                     'complete_angular_moment_over_C':str(moment),'regular_ordinary_bounds':[str(v) for v in ordinary],
                     'common_complete_scalar_reserve':'1e160 ell y^-1 L^2',
                     'radius_multiplied_reserve':'1e176 ell y^3 L^2','full_angular_tower_retained':True,
                     'short_separation_exponent':'exp(-M/1000)','fixed_vacuum_radius_multiplier_retained':True})

def local_polynomials():
    def mul(a,b):return [sum(a[i]*b[j-i] for i in range(j+1)) for j in range(9)]
    def expjet(a):
        out=[F(1)]+[F(0)]*8
        for j in range(1,9):out[j]=sum(k*a[k]*out[j-k] for k in range(1,j+1))/j
        return out
    low=F(1,10000)*F('3.31')
    check('finite_outer_high_background_jets',all(F('1e14')*100**j*low+2**j<F('1e6')**j for j in range(6,9)))
    check('finite_outer_low_background_jets',all(100**j<F('1e6')**j for j in range(1,6)))
    first=mul(mul([F(2,10**6)**j for j in range(9)],[F(1,1000)**j for j in range(9)]),expjet([F(0)]+[F(2)]*8))
    second=mul(mul([F(j+1)*F(2,10**6)**j for j in range(9)],[F(j+1)*F(1,1000)**j for j in range(9)]),expjet([F(0)]+[F(3)]*8))
    check('original_radius_first_path_jet_reserve',max(first)*F('1.001')*F('1e20')<F('1e25'))
    check('original_radius_second_path_jet_reserve',max(second)*F('1.001')*F('1e40')<F('1e50'))
    u=S.symbols('u0:9');w=S.symbols('w0:9');p=S.symbols('p0:9')
    h,z,e,q=S.symbols('h z e q');variables=(*u[1:],*w[1:],*p,h,z,e,q)
    def dy(f):
        return S.expand(-h*u[1]*S.diff(f,h)-z*p[1]*S.diff(f,z)+sum(S.diff(f,a[j])*a[j+1] for a in (u,w,p) for j in range(8)))
    def dp(f):return S.expand(h*dy(f))
    rr,uu,vv,pp=S.symbols('rr uu vv pp');aa=S.Rational(11,360);bb=S.Rational(1,30);cc=S.Rational(1,20);dd=bb-2*aa/3
    lag=pp*(cc*(rr+2*e)**2/3-4*aa*e*rr-q/3)+bb*rr*uu-3*bb*uu**2+(bb*rr+(2*bb-4*aa)*e)*vv-6*dd*uu*vv-3*dd*vv**2
    lp,lq,lz=[S.diff(lag,v) for v in (rr,uu,vv)]
    R=-2*(dp(dp(w[0]))+dp(w[0])**2);U=dp(dp(p[0]))+dp(w[0])*dp(p[0]);V=dp(p[0])**2
    sub={rr:R,uu:U,vv:V,pp:p[0]};lag,lp,lq,lz=[S.expand(v.subs(sub)) for v in (lag,lp,lq,lz)]
    mt=-lag/2+R*lp/2+dp(dp(lp))+(dp(lq)*dp(p[0])+lq*U)/2
    my=-lag/2+R*lp/2+dp(w[0])*dp(lp)+lz*V-dp(lq)*dp(p[0])/2+lq*U/2
    trace=S.expand(mt+my)
    check('Weyl_magnetic_orbit_trace_retained',S.expand(trace).coeff(q)==p[0]/3)
    index={v:i for i,v in enumerate(variables)}
    orders={v:j for a in (u,w,p) for j,v in enumerate(a) if v in index}
    H=F('1e6');A=F('1e20')
    def budget(expr,j,tag='first'):
        total=F(0);terms=S.Poly(expr,*variables).terms()
        for ex,co in terms:
            deriv=sum(orders.get(v,0)*ex[i] for i,v in enumerate(variables))
            zp=ex[index[z]]
            assert deriv<=j+zp,(j,deriv,zp)
            value=sf(abs(co))*H**deriv*F('1.01')**ex[index[h]]*2**ex[index[p[0]]]*160000**ex[index[q]]
            for v,n in orders.items():value*=factorial(n)**ex[index[v]]
            degree=sum(ex[index[v]] for v in orders if v not in w)+ex[index[h]]+zp
            if tag=='weyl':degree=sum(ex[index[v]] for v in orders)+ex[index[h]]+zp
            if tag=='second':value*=degree**2*F('1e50')+degree*F('1e50')
            elif tag=='mixed':value*=degree*F('1e25')
            else:value*=degree*A
            total+=value
        return total,len(terms)
    weyl=z**4*trace;lll=z**4*R
    kinetic=z**2*h**2*(2*p[2]+2*w[1]*p[1]+4*p[1]**2)
    pot=20*(1+z**2+z**4+z**5+z**6);classical=pot+kinetic;mixed=2*pot
    records=[]
    for j in range(5):
        wb,nw=budget(weyl,j,'weyl');lb,nl=budget(lll,j,'weyl')
        kb,nk=budget(classical,j,'second');mb,nm=budget(mixed,j,'mixed')
        check('Weyl_tagged_complete_order_'+str(j),F('6.1')*F('8.5e-33')*wb/3<F('1e110'),dec(wb))
        check('LLL_tagged_complete_order_'+str(j),F('6.1')*F('1e-30')*lb<F('1e110'),dec(lb))
        check('classical_quadratic_order_'+str(j),F('6.1')*(kb+mb)/2<F('1e105'),dec(kb+mb))
        records.append({'order':j,'Weyl_monomials':nw,'Weyl_tagged_budget':str(wb),'LLL_tagged_budget':str(lb),'classical_second_budget':str(kb+mb)})
        if j<4:weyl=dy(weyl);lll=dy(lll);classical=dy(classical);mixed=dy(mixed)
    a,v,t=S.symbols('a v tau');ss=S.log(1+t*(S.exp(a*v)-1))
    s2=t*(1-t)*v*v*S.exp(a*v)/(1+t*(S.exp(a*v)-1))**2
    check('exact_original_radius_second_path_derivative',S.simplify(S.diff(ss,a,2)-s2)==0)
    check('local_total_below_orbit_reserve',F('1e105')+2*F('1e110')+F('1e177')<F('1e180'))
    FIXTURES.append({'fixture':'complete_outer_local_source_polynomials','records':records,
                     'physical_charge_retained':True,'Weyl_magnetic_trace_retained':True,
                     'original_radius_second_derivative':str(s2),'orbit_ordinary_reserve':'1e180 ell^2 y^3 L^2'})

def tensor_volume_and_global():
    B=F('1e180');gamma=[F('.3'),F(4000),F('2e6')];cc={(0,k):B for k in range(5)}
    for r in range(4):
        for k in range(4-r):cc[r+1,k]=cc[r,k+1]+r*sum(comb(k,i)*gamma[i]*cc[r,k-i] for i in range(k+1)) if r else cc[r,k+1]
    seed=[cc[j,0]*F('.42')**j/factorial(j) for j in range(5)]
    geom=fixture('finite-refreshed-local-mixed','refreshed_finite_component_geometry')
    conn=[F(v) for v in geom['geometry']['C'][:4]];phys={(0,k):seed[k] for k in range(5)}
    for n in range(4):
        for k in range(4-n):phys[n+1,k]=(k+1)*phys[n,k+1]+n*sum(conn[i]*phys[n,k-i] for i in range(k+1))
    P=[2*F('2.5')**n*phys[n,0] for n in range(5)]
    check('actual_full_four_connection_used',conn[0]>27 and conn[1]>216 and len(P)==5)
    volume=36*sum(F(comb(4,k)*factorial(k),3**(k+1)) for k in range(5))
    check('outer_volume_density',F('1.001')*F('.42')*3**4<36)
    check('outer_squared_weight_integral',volume<64,str(volume))
    scalar=[F(1001,24000)*F(v) for v in fixture('finite-radius-operator-source','complete_refined_scalar_and_anomaly')['physical_covariant_bounds'][:5]];scalar[0]=min(scalar[0],F('6e-12'))
    co=(sum(P)+sum(scalar))/6
    ci=F(fixture(PREVIOUS,'complete_interior_shape_budget')['C_I_exact']);cb=co+ci;cf=F('1e4')+F('1e15')*cb
    check('complete_global_shape_H4',cb<F('1e200'),dec(cb))
    check('complete_global_TT_H4_forcing',cf<F('1e215'),dec(cf))
    FIXTURES.append({'fixture':'global_physical_shape_and_TT_source','physical_orbit_coefficients':[str(x) for x in P],
                     'physical_orbit_decimals':[dec(x) for x in P],'actual_connection':[str(x) for x in conn],
                     'volume_integral_reserve':str(volume),'C_I':str(ci),'C_O':str(co),'C_b':str(cb),'C_F4':str(cf),
                     'C_O_decimal':dec(co),'C_b_decimal':dec(cb),'C_F4_decimal':dec(cf),
                     'C_b_power':power10(cb),'C_F4_power':power10(cf),
                     'quadratic_factor':'ell^2/g','full_four_dimensional_tensor_slots_retained':True,
                     'higher_inverse_known':False,'useful_solve_target_passed':False})

def independent_dilation_fixture():
    y=S.symbols('y',positive=True);base=y**3/7+y**4*S.log(y)/9
    derivatives=[S.lambdify(y,S.diff(base,y,a),'mpmath') for a in range(10)]
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            center=mp.mpf('.03');power=mp.mpf('10.5')
            def bd(x,a):
                rough=mp.mpf(0)
                if x>center:rough=10*mp.rf(power-a+1,a)*(x-center)**(power-a)
                return derivatives[a](x)+rough
            for a in [0,5]:
                jets=[]
                for j in range(5):
                    jets.append(sum(st2(j,k)*comb(a,r)*mp.mpf(factorial(k))/factorial(k-r)*center**(k-r)*bd(center,a+k-r) for k in range(j+1) for r in range(min(a,k)+1)))
                for sign in [-1,1]:
                    errs=[]
                    for step in ['.001','.0005']:
                        ss=sign*mp.mpf(step);actual=mp.exp(a*ss)*bd(mp.exp(ss)*center,a)
                        jet=sum(ss**j*jets[j]/factorial(j) for j in range(5));errs.append(abs(actual-jet))
                    check(f'C10_non_C11_dilation_refinement_{dps}_{a}_{sign}',errs[1]<mp.mpf('.04')*errs[0],mp.nstr(errs[1]/errs[0],15))
                    records.append({'dps':dps,'spatial_order':a,'sign':sign,'errors':[mp.nstr(v,20) for v in errs]})
    check('fixture_regular_not_analytic',S.Rational(21,2)>10 and S.Rational(21,2)<11)
    FIXTURES.append({'fixture':'independent_nonanalytic_finite_dilation','records':records,'regularity':'C10 but not C11 at X=.03','actual_geometry_sampled':False})

def independent_heat_and_volume():
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            for ss in [mp.mpf('-.01'),mp.mpf('.02')]:
                tt=mp.mpf('.003');yy=mp.mpf('.04');kk=mp.mpf('.17')
                bridge=lambda r:mp.sin(mp.pi*r)+mp.sin(2*mp.pi*r)/3
                def heat(t,y,scale):
                    f=lambda x:mp.exp(kk*scale*x)
                    ff=mp.quad(lambda r:f(y+mp.sqrt(t)*bridge(r))**2,[0,1])
                    gg=mp.quad(lambda r:kk*scale*f(y+mp.sqrt(t)*bridge(r)),[0,1])
                    return f(y)/mp.sqrt(ff)*mp.exp(t*gg*gg/(4*ff))/(2*mp.pi*t)
                left=heat(mp.exp(2*ss)*tt,mp.exp(ss)*yy,1);right=mp.exp(-2*ss)*heat(tt,yy,mp.exp(ss))
                err=abs(left-right)
                check(f'independent_bridge_heat_scale_{dps}_{ss}',err<mp.mpf('1e-25'),mp.nstr(err,15))
                masses=[mp.mpf(2),mp.mpf(3),mp.mpf(5)];deg=[2,4,6]
                angular=lambda t:sum(d*m*m*mp.exp(-m*m*t) for d,m in zip(deg,masses))
                amp=mp.exp(mp.exp(ss)*yy)
                direct=amp*sum(2*d*m*m/(m*m+1)**3 for d,m in zip(deg,masses))
                scaled=mp.quad(lambda t:angular(mp.exp(2*ss)*t)*mp.exp(6*ss)*t*t*mp.exp(-mp.exp(2*ss)*t)*amp,[0,1,mp.inf])
                omitted=mp.quad(lambda t:angular(t)*mp.exp(6*ss)*t*t*mp.exp(-mp.exp(2*ss)*t)*amp,[0,1,mp.inf])
                check(f'independent_angular_scale_integral_{dps}_{ss}',abs(direct-scaled)<mp.mpf('1e-25'))
                check(f'omitted_angular_scaling_detected_{dps}_{ss}',abs(direct-omitted)>mp.mpf('1e-5'))
                k2=mp.exp(mp.exp(ss)*yy);dilated=mp.exp(4*ss)*k2
                check(f'common_k2_scale_retained_{dps}_{ss}',abs(mp.exp(-4*ss)*dilated-k2)<mp.mpf('1e-25') and abs(dilated-k2)>mp.mpf('.01'))
                records.append({'dps':dps,'s':str(ss),'heat_error':mp.nstr(err,20),'angular_integral_error':mp.nstr(abs(direct-scaled),20),'omitted_angular_error':mp.nstr(abs(direct-omitted),20)})
            exact=36*sum(mp.mpf(comb(4,k)*factorial(k))/3**(k+1) for k in range(5))
            numeric=36*mp.quad(lambda t:mp.exp(-3*t)*(1+t)**4,[0,1,10,mp.inf])
            check('independent_outer_volume_integral_'+str(dps),abs(exact-numeric)<mp.mpf('1e-25'))
    FIXTURES.append({'fixture':'independent_heat_scaling_and_outer_volume','records':records,'actual_quantum_state_sampled':False})

def independent_outer_residual():
    y,a=S.symbols('y a',positive=True);phi=-S.log(y)+a*y**3/10+a*y**4*S.log(y)/20
    u=y**3/30;z2=S.exp(-2*phi);h2=S.exp(-2*a*u)
    K=6+2*z2+3*z2**2-a*z2**2-z2*h2*(2*S.diff(phi,y,2)+4*S.diff(phi,y)**2)
    Q=S.exp(-4*phi)*(1/y+S.log(y)**2)+S.exp(-4*a*y**3/10)
    original=K+a*Q
    identity=K-K.subs(a,0)-a*S.diff(K,a).subs(a,0)+a*(Q-Q.subs(a,0))
    check('original_seed_tangent_residual_identity',S.simplify(original-original.subs(a,0)-a*S.diff(original,a).subs(a,0)-identity)==0)
    records=[]
    # These are manufactured functions; affine subtraction is diagnostic only.
    for j in [0,4]:
        expr=S.diff(original,y,j);fun=S.lambdify((y,a),expr,'mpmath',cse=True)
        d1=S.lambdify((y,a),S.diff(expr,a),'mpmath',cse=True);d2=S.lambdify((y,a),S.diff(expr,a,2),'mpmath',cse=True)
        for dps in PROTOCOL['digits']:
            with mp.workdps(dps):
                for ys in ['.01','.06']:
                    yy=mp.mpf(ys);ell=mp.mpf('.7');direct=fun(yy,ell)-fun(yy,0)-ell*d1(yy,0)
                    scale=yy**j/(ell**2*yy**3*(1-mp.log(yy))**2)
                    for count in PROTOCOL['quadrature_nodes']:
                        nodes,weights=mp.gauss_quadrature(count,'legendre')
                        val=ell**2*sum(weights[k]*(1-(nodes[k]+1)/2)*d2(yy,ell*(nodes[k]+1)/2)/2 for k in range(count))
                        error=abs((direct-val)*scale)
                        check(f'independent_outer_quadratic_{j}_{dps}_{ys}_{count}',error<mp.mpf('1e-25'),mp.nstr(error,15))
                    records.append({'order':j,'dps':dps,'y':ys,'scaled_residual':mp.nstr(direct*scale,20),'error':mp.nstr(error,20)})
    FIXTURES.append({'fixture':'independent_outer_quadratic_residual','records':records,'affine_subtraction_diagnostic_only':True,'actual_quantum_state_sampled':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 complete outer and global physical shape source','input_count':len(ps)}),flush=True)
    try:
        inputs_and_dilation();gluing_and_source();local_polynomials();tensor_volume_and_global()
        independent_dilation_fixture();independent_heat_and_volume();independent_outer_residual()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
    run={'story':'S167 complete outer and global physical shape source','created_utc':datetime.now(timezone.utc).isoformat(),
         'status':'passed' if good else 'failed','exception':error,'python':sys.version,'sympy':S.__version__,'mpmath':mp.__version__,
         'protocol':PROTOCOL,'checks':ROWS,'fixtures':FIXTURES,'input_sha256_before':before,'input_sha256_after':after,'unchanged_inputs':before==after}
    for k,v in review().items():
        if isinstance(v,bool) and k not in run:run[k]=v
    for kind,p in zip(['checker','note','review'],OWN):run['frozen_'+kind+'_source_utf8']=(ROOT/p).read_text(encoding='utf-8')
    dest=ROOT/'evidence'/(STEM+'.json');data=json.loads(dest.read_bytes()) if dest.exists() else {'runs':[]}
    data['runs'].append(run);dest.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({'status':run['status'],'checks':len(ROWS),'unchanged_inputs':before==after}),flush=True)
    return 0 if good else 1

if __name__=='__main__':raise SystemExit(main())
