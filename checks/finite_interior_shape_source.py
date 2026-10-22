"""Complete finite physical shape forcing from the bolt through the middle cut."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
from math import factorial,comb
from itertools import product
import importlib.util,json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='finite-interior-shape-source';PREVIOUS='finite-middle-source-four'
OWN=['checks/finite_interior_shape_source.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[];N=8
PROTOCOL={'source_order':4,'local_metric_order':8,'real_path':[0,1],'complex_path_radius':1,
          'physical_dimension':4,'restricted_cut':40,'numerical_dps':[60,90],
          'quadrature_nodes':[32,64],'contour_nodes':[32,64],'numerical_tolerance':'1e-25',
          'manufactured_radii':['0','.3'],'manufactured_parameters':['0','.4','1'],
          'manufactured_spatial_orders':[0,4],'actual_quantum_state_sampled':False}

def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row);print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)

def zero(name,value):
    value=S.simplify(value);check(name,value==0,value)

def saved(stem):return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def fixture(stem,name):return next(f for f in saved(stem)['fixtures'] if f['fixture']==name)
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})
def strings(xs):return [str(v) for v in xs]
def mf(v):return mp.mpf(v.numerator)/v.denominator if isinstance(v,F) else mp.mpf(v)
def decimal(v):
    with mp.workdps(35):return mp.nstr(mf(v),20)
def poly(a=()):return [F(a[j]) if j<len(a) else F(0) for j in range(N+1)]
def add(*ps):return [sum(p[j] for p in ps) for j in range(N+1)]
def scale(a,p):return [F(a)*v for v in p]
def mul(a,b):return [sum(a[i]*b[j-i] for i in range(j+1)) for j in range(N+1)]
def derivative(p):return [(j+1)*p[j+1] if j<N else F(0) for j in range(N+1)]
def exp0(p):
    assert p[0]==0
    out=poly([1])
    for j in range(1,N+1):out[j]=sum(k*p[k]*out[j-k] for k in range(1,j+1))/j
    return out
def minus_log1(p):
    assert p[0]==0
    out=poly();power=poly([1])
    for j in range(1,N+1):power=mul(power,p);out=add(out,scale(F(1,j),power))
    return out

def inputs_and_cold():
    old=saved(PREVIOUS)
    check('previous_finite_middle_passed',old['status']=='passed' and len(old['checks'])==91 and all(r['passed'] for r in old['checks']))
    check('previous_inputs_unchanged',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    check('global_shape_and_solution_remain_open',not old['actual_complete_shape_source_H4_bounded'] and not old['nonlinear_solution_constructed'])
    radius=F(1,32);ee=F('1.005')
    check('cold_normal_radius_invariant',2*F('.000501')*radius+F('.07')*ee*radius**2/4<F('.00005'))
    check('cold_normal_conformal_invariant',F('2.31')*ee*radius**2/8<F('.002'))
    check('cold_normal_gradient',3*F('2.31')*ee*radius/8<F('.067'))
    check('cold_Picard_row_margins',max(ee*radius**2*(10+2*F('2.31'))/8,ee*radius**2*(F('2.31')+2*F('.07'))/4)<F('.02'))
    check('complex_path_coefficient_cone',2*F('.00013')+F('.002')<F('.011'))
    check('complex_mass_coercivity',F('.989')**2<F('.98') and (F('.989')-F('3e-8'))**2>F('.97'))
    check('retained_energy_factor_for_complex_path',F('.97')*F('1.1')**2>1)
    check('doubled_cap_global_gradient',F('.067')+2*F('8e-6')+640*(F('.002')+2*F('1.2e-8'))<F('1.35'))
    z=poly([0,F('.07')]+[F('.027')*128**j+F('.00026')*150**j for j in range(2,6)])
    aa=exp0(z);endpoint=exp0(scale(F(1,2),z))
    check('complex_mass_derivative_reserves',all(F('1.01')*aa[j]<F('.1')*300**j for j in range(2,6)))
    check('complex_endpoint_derivative_reserves',all(F('1.01')*endpoint[j]<F('.1')*300**j for j in range(1,6)))
    geo=fixture(PREVIOUS,'finite_metric_and_fixed_cap')
    check('doubled_C9_stays_in_translation_reserve',2*max(F(v) for v in geo['coordinate_u_jets'])+F('1e86')<F('1e88'))
    check('fixed_vacuum_value_reserve',F('4.352e-10')+F('1.25e-11')+F('1e-28')<F('1e-8'))
    cap=fixture(PREVIOUS,'finite_translation_interpolation_and_Cauchy');glue=fixture(PREVIOUS,'finite_coercive_rows_and_complete_gluing')
    mass=[F(a)+F(b) for a,b in zip(cap['cap_coordinate_bounds'],glue['normalized_gluing_bounds'])];mass[0]+=F('1e-8')
    FIXTURES.append({'fixture':'cold_bolt_and_complex_path_extension','complex_path_radius':1,'finite_field_size_max':2,
                     'massive_coordinate_bounds':strings(mass),'massive_coordinate_decimals':[decimal(x) for x in mass],
                     'fixed_vacuum_scalar_bound':'1e-8','bolt_included':True,'inverse_thermal_gap_used':False,
                     'full_angular_tower_and_actual_domain_retained':True})
    return mass

def local_field_jets():
    fields=fixture('finite-low-gap-reconstruction','complete_improved_tangent_arrays');records={}
    gamma=[F('1e4')*factorial(i)*1000**i for i in range(8)]
    for name in ['U','Vbar']:
        data=[F(v) for v in fields[name][:9]]
        b={(n,0):F(factorial(n))*400**n*data[n] for n in range(9)}
        for k in range(8):
            for n in range(8-k):b[n,k+1]=b[n+1,k]+n*sum(comb(k,i)*gamma[i]*b[n,k-i] for i in range(k+1))
        arr=[b[0,j]/factorial(j) for j in range(9)]
        check(name+'_coordinate_derivatives_through_eight',all(arr[j]<F('1e7')**j for j in range(1,9)))
        records[name]=strings(arr)
    check('seed_radius_and_normal_factor_jets',all(6*640**j<F('1e7')**j for j in range(1,9)))
    h=poly([0]+[F('1e7')**j for j in range(1,9)]);one=poly([1]);eh=exp0(h);e2h=exp0(scale(2,h))
    a=add(scale('3e-8',add(eh,scale(-1,one))),scale('1.001',mul(eh,add(e2h,scale(-1,one)))))
    pp=add(h,minus_log1(scale(1/(1-F('3e-8')),a)));p=poly([0]+[F('1e9')**j for j in range(1,9)])
    check('complex_nonlinear_radius_jets',all(pp[j]<p[j] and 3*h[j]<p[j] for j in range(1,9)))
    check('complex_radius_does_not_vanish',1-3*F('1e-8')>F('.999999'))
    check('radius_inverse_powers_value',1/(1-F('3e-8'))**6<F('1.001'))
    check('retained_physical_charge',F('1.01')/F('.99')**2<2)
    check('retained_seed_potential_coefficients',6/F('5.9')<2 and 2*(F('1e-10')+3/F('5.9'))<2)
    check('LLL_physical_coefficient',F('8.5e-33')*400/36<F('1e-30'))
    FIXTURES.append({'fixture':'complex_finite_local_metric_jets','coordinate_divided_tangent_bounds':records,
                     'nonlinear_radius_divided_bounds':strings(pp),'common_divided_bound_base':'1e9',
                     'physical_charge_distinct_from_comparison_charge':True,'large_high_jets_used_as_values':False})
    return p

def exact_orbit_normalization():
    r=S.symbols('rho',positive=True);gg=S.symbols('G',positive=True);xx=S.Function('X')(r);ss=S.Function('sigma')(r);vv=S.Function('V')
    lag=(-2*xx*S.diff(ss,r,2)+S.diff(xx,r)**2/(2*xx)+S.exp(2*ss)*vv(xx))/(-4*gg)
    euler=S.diff(lag,ss)-S.diff(S.diff(lag,S.diff(ss,r)),r)+S.diff(S.diff(lag,S.diff(ss,r,2)),r,2)
    zero('independent_reduced_conformal_Euler',euler-(S.diff(xx,r,2)-S.exp(2*ss)*vv(xx))/(2*gg))
    zero('physical_volume_orbit_gradient',euler/(4*S.pi*xx*S.exp(2*ss))-(S.exp(-2*ss)*S.diff(xx,r,2)-vv(xx))/(8*S.pi*gg*xx))
    B,C,lam,rr,theta=S.symbols('B C lambda r Theta',positive=True);g=gg*lam;x=lam*B**2;scan=2*gg*theta/(B**2*C)
    zero('canonical_quantum_physical_normalization',-theta/(4*S.pi*lam**2*rr**4)+(C/x)*(B/rr)**4*scan/(8*S.pi*g))
    o,e,scalar,tt=S.symbols('O E scalar t')
    zero('complete_shape_covector_identity',((tt-scalar/3)/2+o/(16*S.pi*g)+scalar/6).subs(tt,-o/(8*S.pi*g)))
    p=S.Function('phi')(r);sig=S.Function('can_sigma')(r)
    zero('radius_Laplacian_in_canonical_units',S.exp(-2*p)*S.diff(S.exp(2*p),r,2)-(2*S.diff(p,r,2)+4*S.diff(p,r)**2))
    z=S.symbols('z',positive=True);fc=S.Function('F')(z)
    zero('seed_orbit_equation_identity',2*(fc+z*S.diff(fc,z))-2*(-fc+z*S.diff(fc,z))-4*fc)
    FIXTURES.append({'fixture':'exact_complete_orbit_and_shape_normalization',
                     't':'-O_ell/(8 pi g)','b_src':'-O_ell/(16 pi g)-s/6',
                     'orbit_source_uses_trace_only':True,'canonical_pressure_reconstruction_required':False,
                     'fixed_vacuum_radius_multiplier_retained':True})

def independent_four_geometry():
    rho,time,th,ph=S.symbols('rho time theta azimuth',positive=True);coords=[rho,time,th,ph]
    a=S.Function('a')(rho);r=S.Function('r')(rho);metric=S.diag(1,a*a,r*r,r*r*S.sin(th)**2);inv=metric.inv()
    conn={(i,j,k):S.simplify(sum(inv[i,l]*(S.diff(metric[l,k],coords[j])+S.diff(metric[j,l],coords[k])-S.diff(metric[j,k],coords[l]))/2 for l in range(4))) for i,j,k in product(range(4),repeat=3)}
    def ric(b,d):return S.simplify(sum(S.diff(conn[k,d,b],coords[k])-S.diff(conn[k,k,b],coords[d])+sum(conn[k,k,l]*conn[l,d,b]-conn[k,d,l]*conn[l,k,b] for l in range(4)) for k in range(4)))
    rd=[ric(j,j) for j in range(4)];sc=sum(inv[j,j]*rd[j] for j in range(4));orbit=S.simplify(rd[0]+rd[1]/a**2-sc)
    q,l0=S.symbols('Q lambda0');potential=2+6*l0*r*r-2*q/(r*r);X=r*r
    reduced=(S.diff(X,rho,2)+S.diff(a,rho)*S.diff(X,rho)/a-potential)/X
    zero('independent_four_dimensional_orbit_trace',orbit-6*l0+2*q/r**4-reduced)
    for alpha,beta in [(S.Rational(1,10),S.Rational(1,5)),(S.Rational(7,20),S.Rational(1,7))]:
        val=S.simplify(reduced.subs({a:rho*S.exp(alpha*rho*rho),r:S.exp(beta*rho*rho)}).doit())
        bolt=S.simplify(S.limit(val,rho,0));expected=8*beta-2-6*l0+2*q
        zero('manufactured_regular_bolt_orbit_'+str(alpha)+'_'+str(beta),bolt-expected)
    FIXTURES.append({'fixture':'independent_four_dimensional_orbit_trace','dimension':4,
                     'orbit_Einstein_trace':str(orbit),'physical_charge_retained':True,'regular_bolt_examples':2})

def complete_weyl_polynomial():
    spec=importlib.util.spec_from_file_location('retained_complete_Weyl',ROOT/'checks/comparison_path_complete_source_response.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    module.check=lambda name,ok,detail='':check('retained_Weyl_'+name,ok,detail)
    H,e,qm,sig,phi,mt,null,_,_=module.cartesian_polynomials();trace=S.expand(2*mt-null)
    for j in range(6,9):
        for k in range(j+1):sig[j-k,k]=S.Symbol('s'+str(j-k)+str(k));phi[j-k,k]=S.Symbol('p'+str(j-k)+str(k))
    directions=[]
    for axis in range(2):
        d={H:-2*H*sig[(1,0) if axis==0 else (0,1)]}
        for jets in [sig,phi]:
            for (i,j),symbol in jets.items():
                if i+j<8:d[symbol]=jets[(i+1,j) if axis==0 else (i,j+1)]
        directions.append(d)
    def D(expr,axis):return S.expand(sum(S.diff(expr,v)*directions[axis][v] for v in expr.free_symbols if v in directions[axis]))
    variables=[H,e,qm]+list(sig.values())+list(phi.values());orders={v:sum(k) for jets in [sig,phi] for k,v in jets.items()}
    bound={H:F('1.01'),e:F(1),qm:F(160000),sig[0,0]:F(1),phi[0,0]:F(4)}
    for v,order in orders.items():
        if order:bound[v]=F(factorial(order))*F('1e9')**order
    row=[trace];records=[]
    for j in range(5):
        if j:row=[D(expr,0) for expr in row]+[D(row[-1],1)]
        values=[];counts=[]
        for expr in row:
            terms=S.Poly(expr,*variables).terms();value=F(0)
            for powers,co in terms:
                term=F(str(abs(co)))
                assert sum(p*orders.get(v,0) for p,v in zip(powers,variables))<=4+j
                assert powers[0]<=2 and powers[variables.index(phi[0,0])]<=1
                for p,v in zip(powers,variables):
                    if p:term*=bound[v]**p
                value+=term
            values.append(value);counts.append(len(terms))
        tensor=2**((j+1)//2)*max(values);reserve=F('1e10')*factorial(4+j)*F('1e9')**(4+j)
        check('complete_Weyl_orbit_tensor_order_'+str(j),tensor<reserve,decimal(tensor/reserve))
        records.append({'order':j,'component_monomial_counts':counts,'tensor_bound':str(tensor),'reserve':str(reserve)})
    FIXTURES.append({'fixture':'complete_Weyl_orbit_trace_four_derivatives','records':records,
                     'magnetic_trace_monomial':'qbar phi/3','all_spatial_multiindices_retained':True})

def full_interior_budget(mass,p):
    one=poly([1]);H=scale('1.01',exp0(scale(2,p)));em={m:scale('1.001',exp0(scale(m,p))) for m in [2,4,5,6]}
    p0=add(scale('2e-10',em[2]),scale(2,one),scale(2,em[4]),scale('2e-24',add(em[4],scale(4,em[5]),scale(3,em[6]))))
    pp=add(scale('2e-10',em[2]),scale(6,one),scale(4,em[4]));pot=add(scale(3,p0),scale(2,pp))
    dp=derivative(p);ddp=derivative(dp);kin=mul(mul(em[2],H),scale(4,add(ddp,mul(dp,dp))))
    lll=scale('1e-30',mul(scale(4,H),ddp));weyl=poly([F('1e10')*factorial(4+j)*F('1e9')**(4+j)/factorial(j) for j in range(5)])
    massive=poly([mass[j]/factorial(j) for j in range(5)])
    quantum=scale(2,mul(em[4],add(massive,lll,scale(F('8.5e-33')/3,weyl))))
    out=scale('6.1',add(pot,kin,quantum));raw=[factorial(j)*out[j] for j in range(5)]
    check('complete_orbit_bound_positive',all(v>0 for v in raw))
    gamma=[F('.3'),F(4000),F('2e6')];cc={(0,k):raw[k] for k in range(5)}
    for r in range(4):
        for k in range(4-r):cc[r+1,k]=cc[r,k+1]+r*sum(comb(k,i)*gamma[i]*cc[r,k-i] for i in range(k+1)) if r else cc[r,k+1]
    seed=[cc[j,0]/factorial(j) for j in range(5)]
    old=fixture('finite-refreshed-local-mixed','refreshed_finite_component_geometry');conn=[F(v) for v in old['geometry']['C'][:4]]
    phys={(0,k):seed[k] for k in range(5)}
    for n in range(4):
        for k in range(4-n):phys[n+1,k]=(k+1)*phys[n,k+1]+n*sum(conn[i]*phys[n,k-i] for i in range(k+1))
    output=[2*F('2.5')**n*phys[n,0] for n in range(5)]
    check('actual_physical_connection_not_flat',conn[0]>27 and conn[1]>216)
    check('physical_volume_without_inverse_gap',F('1.001')*F(40**3-1,3)<150**2)
    check('real_density_value_reserve',1/(1-6*F('3e-8'))<F('1.001'))
    sf=fixture('finite-radius-operator-source','complete_refined_scalar_and_anomaly');scalar=[F(1001,24000)*F(v) for v in sf['physical_covariant_bounds'][:5]];scalar[0]=min(scalar[0],F('6e-12'))
    bound=F(25,8)*sum(output)+sum(scalar)/6
    check('complete_physical_interior_shape_H4',bound<F('1e110'),decimal(bound))
    power=0
    while F(10)**power<bound:power+=1
    check('upward_interior_shape_enclosure',F(10)**(power-1)<bound<=F(10)**power)
    FIXTURES.append({'fixture':'complete_interior_shape_budget','orbit_coordinate_bounds':strings(raw),
                     'orbit_coordinate_decimals':[decimal(v) for v in raw],
                     'full_physical_orbit_bounds':strings(output),'full_physical_orbit_decimals':[decimal(v) for v in output],
                     'actual_connection_divided_coefficients':strings(conn),'scalar_H4_coefficients':strings(scalar),
                     'C_I_exact':str(bound),'C_I_decimal':decimal(bound),'C_I_power_of_ten':power,'reserve':'1e110',
                     'restricted_region':'regular bolt <= t_seed <= 40','quadratic_path_factor':'ell^2/g',
                     'physical_volume_factor':150,'full_four_dimensional_tensor_slots_retained':True,
                     'complete_global_shape_H4_supplied':False,'target_passed':False})

def independent_quadratic_residual():
    rho,ell=S.symbols('rho ell',real=True);rr=1+rho**2*S.exp(ell*rho**2/50);a=rho*S.exp(rho**2/10+3*ell*rho**2/100);X=rr*rr
    v0=2+S.Rational(3,4)*X-S.Rational(2,3)/X;vp=2+S.Rational(6,7)*X-S.Rational(4,5)/X
    op=S.cancel(((1-ell)*v0+ell*vp-S.diff(X,rho,2)-(1/rho+rho/5+3*ell*rho/50)*S.diff(X,rho))/X+ell*S.exp(ell*rho*rho)/(1+rho*rho)**4)
    derivatives={j:S.diff(op,rho,j) for j in PROTOCOL['manufactured_spatial_orders']};records=[]
    for dps in PROTOCOL['numerical_dps']:
        with mp.workdps(dps):
            for at in PROTOCOL['manufactured_radii']:
                rr0=S.Rational(at);functions={}
                for j,expr in derivatives.items():
                    f=S.simplify(expr.subs(rho,rr0));functions[j]=[S.lambdify(ell,S.diff(f,ell,k),'mpmath',cse=True) for k in range(3)]
                for arg in PROTOCOL['manufactured_parameters']:
                    lam=mp.mpf(arg)
                    for order,ff in functions.items():
                        f,df,ddf=ff;direct=f(lam)-f(0)-lam*df(0);last=None
                        for count in PROTOCOL['quadrature_nodes']:
                            nodes,weights=mp.gauss_quadrature(count,'legendre')
                            integral=lam*lam*sum(weights[k]*(1-(nodes[k]+1)/2)*ddf(lam*(nodes[k]+1)/2)/2 for k in range(count))
                            error=abs(direct-integral);tag=f'{dps}_{at}_{arg}_{order}_{count}'
                            check('independent_real_quadratic_Taylor_'+tag,error<mp.mpf('1e-25'),mp.nstr(error,18))
                            last=integral
                        for count in PROTOCOL['contour_nodes']:
                            radius=mp.mpf('.25');values=[mp.exp(2j*mp.pi*k/count) for k in range(count)]
                            cauchy=2*sum(f(lam+radius*z)/(z*z) for z in values)/(count*radius**2)
                            error=abs(cauchy-ddf(lam));tag=f'{dps}_{at}_{arg}_{order}_{count}'
                            check('independent_complex_Cauchy_'+tag,error<mp.mpf('1e-25'),mp.nstr(error,18))
                        records.append({'dps':dps,'rho':at,'ell':arg,'spatial_order':order,'direct_residual':mp.nstr(direct,24),
                                        'Taylor_integral':mp.nstr(last,24),'last_contour_error':mp.nstr(error,18)})
    FIXTURES.append({'fixture':'independent_regular_metric_quadratic_residual','records':records,
                     'actual_quantum_state_sampled':False,'affine_Taylor_subtraction_is_diagnostic_only':True,
                     'regular_bolt_included':True,'two_dimensional_orbit_trace_used':True})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps};before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 complete finite interior physical shape source','input_count':len(ps)}),flush=True)
    try:
        mass=inputs_and_cold();p=local_field_jets();exact_orbit_normalization();independent_four_geometry()
        complete_weyl_polynomial();full_interior_budget(mass,p);independent_quadratic_residual()
    except Exception:
        error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
    run={'story':'S167 complete finite interior physical shape source','created_utc':datetime.now(timezone.utc).isoformat(),
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
