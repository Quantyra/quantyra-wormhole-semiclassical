"""Fixed-domain distant-seam proof budgets and independent resolvent mechanisms."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as F
import json,sys,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1];STEM='fixed-domain-cold-far-seam'
PREVIOUS='finite-cold-admission-and-response'
OWN=['checks/fixed_domain_cold_far_seam.py','notes/'+STEM+'.md','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
PROTOCOL={'digits':[60,90],'Cauchy_nodes':[32,64],'identity_tolerance':'1e-35','finite_difference_steps':['1e-4','5e-5'],'finite_difference_relative_tolerance':'1e-8','refinement_ratio':'.27','actual_state_sampled':False}

def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    row={'name':name,'passed':bool(ok),'detail':str(detail)};ROWS.append(row)
    print(json.dumps(row),flush=True)
    if not ok:raise AssertionError(name)
def saved(stem):return json.loads((ROOT/'evidence'/(stem+'.json')).read_bytes())['runs'][-1]
def review():return json.loads((ROOT/OWN[2]).read_bytes())
def paths():return sorted(set(saved(PREVIOUS)['input_sha256_after'])|set(OWN)|set(review()['retained_inputs'])|{'evidence/'+PREVIOUS+'.json'})
def mf(x):return mp.mpf(x.numerator)/x.denominator if isinstance(x,F) else mp.mpf(x)



def cone_and_energy():
    old=saved(PREVIOUS)
    check('prior_audit_preserved',old['status']=='passed' and not old['general_cold_family_admitted'] and not old['fixed_domain_differentiated_source_proof'])
    check('prior_input_hashes_preserved',old['input_sha256_before']==old['input_sha256_after'] and all(sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in old['input_sha256_after'].items()))
    ep=F('2e-6');rzprod=F('1e-5')
    check('finite_tube_epsilon',F('1e-6')+1000*F('5e-10')<ep)
    check('actual_value_coefficient',F('.012')*1000>1)
    check('cap_value_coefficient',F('.03')*1000>1)
    check('actual_gradient_coefficient',F('7.7')<8 and 8000>1)
    check('cap_gradient_coefficient',F('.6')*1000>1)
    check('complex_value_increment',F('.03')*rzprod==F('3e-7'))
    check('complex_gradient_increment',8*rzprod==F('8e-5'))
    check('whole_complex_cap_cone',F('.042')+F('.03')*(ep+rzprod)<F('.043'))
    check('local_complex_low_jet',100*(ep+rzprod)*F('.00331')<F('4e-6'))
    check('evaluation_ball_fits_seed_chart',F(4,20*18)<F(1,32))
    check('first_energy',80+F('1.1')<82)
    check('second_energy',82**2+F(10,18)<6730)
    value=F('2.21')+F(10,18)+160*82+11200
    gradient=(F('2.21')+F(10,18))*82+F(22,18)+F(2000,18**2)+160*6730+11200*82
    check('value_cutoff_energy',value<25000)
    check('gradient_cutoff_energy',gradient<F('2.1e6'))
    check('value_geometric_product',4*2*2*9**2*300000**2<F('2e14'))
    check('gradient_geometric_product',4*2*2*2*9**2*300000*F('1e10')<F('1e19'))
    check('density_gradient_reserve',F('1e19')+20*F('2e14')/18<F('1e20'))
    check('Cauchy_value',F('2e14')/rzprod==F('2e19'))
    check('Cauchy_gradient',F('1e20')/rzprod==F('1e25'))
    check('weighted_value_prefactor',120*F('2e19')==F('2.4e21'))
    check('weighted_gradient_prefactor',120*F('1e25')==F('1.2e27'))
    check('counting_integral',12*F(1,18)==F(2,3))
    check('full_value_allocation',F(2,3)*F('1e20')<F('1e21'))
    check('full_bolt_allocation',2*2*F(2,3)*F('1e20')<F('1e21'))
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            inv=mp.exp(mp.mpf('.043'))/(mp.exp(-mp.mpf('.043'))*mp.cos(mp.mpf('3e-7'))-mp.mpf('.25'))
            q=mp.exp(mp.mpf('.025004'))
            check('weighted_cap_inverse_'+str(dps),inv<2)
            check('local_mass_potential_'+str(dps),q<mp.mpf('1.1') and 2*q*mp.mpf('4.000004')<10 and 2*q*(mp.mpf('400.000004')+mp.mpf('4.000004')**2)<2000)
            check('complex_actual_inverse_'+str(dps),mp.exp(mp.mpf('3e-7'))/(1-(mp.exp(mp.mpf('3e-7'))-1))<2)
            check('commutator_'+str(dps),mp.mpf('7.5')*mp.exp(mp.mpf('.043'))<9)
            c0=25000/(2*mp.sqrt(mp.pi))*2*mp.mpf('1.2')*4
            c1=mp.mpf('2.1e6')/(2*mp.sqrt(mp.pi))*2*mp.mpf('1.2')*4+100*mp.mpf('3e5')/18
            check('point_row_'+str(dps),c0<mp.mpf('3e5'))
            check('gradient_row_'+str(dps),c1<mp.mpf('1e10'))
            records.append({'dps':dps,'weighted_inverse_coefficient':mp.nstr(inv,30),'point_row_derived':mp.nstr(c0,30),'gradient_row_derived':mp.nstr(c1,30)})
    FIXTURES.append({'fixture':'fixed_domain_cone_and_point_rows','records':records,'value_cutoff_coefficient':str(value),'gradient_cutoff_coefficient':str(gradient),'T0':'2.4e21','T1':'1.2e27','complete_cold_admission_proved':False})

def distance_and_sum():
    L=S.symbols('L',positive=True)
    logweight=S.log(S.Rational('.01'))+S.Rational('.055')*(1-L)-S.Rational('.25')+S.log(L)
    check('weight_log_concave',S.diff(logweight,L,2)==-1/L**2)
    check('value_mass_extremum_sign',F(3,18)-F('1.9')/2<0)
    check('gradient_mass_extremum_sign',F(4,18)-60/2<0)
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            ends=[1+mp.log(10),81-mp.log(100)]
            vals=[mp.mpf('.01')*mp.exp(mp.mpf('.055')*(1-y)-mp.mpf('.25'))*y for y in ends]
            check('finite_h_weight_bound_'+str(dps),min(vals)>mp.mpf(1)/120)
            d0=mp.mpf('.945')*mp.log(10)-mp.mpf('.25')
            db=mp.mpf('.945')*mp.log(mp.mpf('.01')/(mp.mpf('1.62')*mp.mpf('1e-30')))-mp.mpf('.25')
            check('whole_cold_separation_'+str(dps),d0>mp.mpf('1.9'))
            check('near_bolt_separation_'+str(dps),db>60)
            val=mp.mpf('2.4e21')*18**3*mp.exp(-mp.mpf('15.2'))
            grad=mp.mpf('1.2e27')*18**4*mp.exp(-480)
            check('complete_value_per_mode_budget_'+str(dps),val<mp.mpf('1e20'))
            check('complete_gradient_per_mode_budget_'+str(dps),grad<mp.mpf('1e20'))
            margin=(mp.mpf('1e20')*mp.exp(mp.mpf('15.2'))/18**3)/mp.mpf('2.4e21')
            check('value_prefactor_margin_'+str(dps),margin>28)
            # Independent numerical integration of the full counting envelope.
            integral=mp.quad(lambda y:12/(18+y)**2,[0,18,mp.inf])
            check('Stieltjes_envelope_integral_'+str(dps),abs(integral-mp.mpf(2)/3)<mp.mpf('1e-40'))
            records.append({'dps':dps,'weight_endpoint_values':[mp.nstr(x,30) for x in vals],'D_cold':mp.nstr(d0,30),'D_bolt':mp.nstr(db,30),'value_M3_coefficient':mp.nstr(val,30),'gradient_M3_coefficient':mp.nstr(grad,30),'value_margin':mp.nstr(margin,30)})
    FIXTURES.append({'fixture':'all_mass_and_distance_extrema','records':records,'angular_cutoff_used':False})

def flat_identities():
    x,y,M=S.symbols('x y M',real=True);I=S.I
    g1=S.Matrix([[0,1],[1,0]]);g2=S.Matrix([[0,-I],[I,0]])
    u=S.Matrix([S.exp(-x*x-y*y), (x+I*y)*S.exp(-x*x-y*y)])
    D=lambda q:g1*q.diff(x)+g2*q.diff(y)
    check('constant_principal_Dirac_square',S.simplify(D(D(u))-u.diff(x,2)-u.diff(y,2))==S.zeros(2,1))
    chi=S.exp(-x*x-2*y*y);C=g1*S.diff(chi,x)+g2*S.diff(chi,y)
    check('continuum_Dirac_Leibniz',S.simplify(D(chi*u)-chi*D(u)-C*u)==S.zeros(2,1))
    check('continuum_second_cutoff_commutator_zero',S.simplify(C*(chi*u)-chi*(C*u))==S.zeros(2,1))
    # Fourier norm from a direct radial integral, not an assigned point constant.
    k,m=S.symbols('k m',positive=True)
    integral=S.integrate(k/(2*S.pi*(k*k+m*m)**2),(k,0,S.oo))
    check('Fourier_point_norm_squared',S.simplify(integral-1/(4*S.pi*m*m))==0)
    energy=lambda q:S.expand_complex((S.conjugate(q.T)*q)[0]).expand()
    diff=S.simplify(energy(D(u))-energy(u.diff(x))-energy(u.diff(y)))
    # Pointwise equality is not required; the cross term integrates to zero.
    integrated=S.integrate(diff,(x,-S.oo,S.oo),(y,-S.oo,S.oo))
    check('Dirac_energy_integrated',integrated==0)
    t=S.symbols('t',real=True);eta=1-10*t**3+15*t**4-6*t**5
    check('seam_quintic_slope',S.diff(eta,t).subs(t,S.Rational(1,2))==-S.Rational(15,8))
    check('seam_quintic_extrema',S.factor(S.diff(eta,t,2))==-60*t*(t-1)*(2*t-1))
    FIXTURES.append({'fixture':'flat_point_and_energy_identities','Fourier_point_norm_squared':str(integral),'energy_cross_integral':str(integrated),'continuum_proof_machine_verified':False})

def finite_domain_mechanism():
    records=[]
    for dps in PROTOCOL['digits']:
        with mp.workdps(dps):
            n=9;D=mp.zeros(n)
            for j in range(n-1):D[j,j+1]=1;D[j+1,j]=-1
            chi=mp.diag([0,0,mp.mpf('.5'),1,1,1,mp.mpf('.5'),0,0]);index=4
            f=[mp.mpf('.02')*(1+mp.mpf(j)/n) for j in range(n)]
            fc=list(f);fc[0]*=2;fc[-1]*=mp.mpf('.5')
            sa=[mp.mpf('.01')*j/n for j in range(n)];sc=list(sa);sc[0]+=mp.mpf('.08');sc[-1]-=mp.mpf('.06')
            fa=mp.diag(f);fb=mp.diag(fc)
            for mass in [2,5]:
                def objects(z):
                    ea=mp.diag([mp.exp(-(sa[j]+z*f[j])/2) for j in range(n)])
                    ec=mp.diag([mp.exp(-(sc[j]+z*fc[j])/2) for j in range(n)])
                    aa=ea*D*ea;ac=ec*D*ec
                    A=(aa+mass*mp.eye(n))**-1;B=(ac+mass*mp.eye(n))**-1;C=ac*chi-chi*ac
                    return aa,ac,A,B,C
                z=mp.mpf('.3')+mp.mpf('.07')*1j;aa,ac,A,B,C=objects(z)
                direct=A[index,index]-B[index,index];J=C*chi-chi*C
                glue=(B*J*B/2-B*C*A*C*B)[index,index]
                err=abs(direct-glue)/abs(direct)
                check('two_sided_fixed_domain_identity_'+str(dps)+'_'+str(mass),err<mp.mpf('1e-35'))
                check('lattice_second_commutator_required_'+str(dps)+'_'+str(mass),abs(direct+(B*C*A*C*B)[index,index])/abs(direct)>1)
                check('operators_agree_on_cutoff_'+str(dps)+'_'+str(mass),mp.norm((aa-ac)*chi)<mp.mpf('1e-45') and mp.norm(chi*(aa-ac))<mp.mpf('1e-45'))
                Ap=(A*fa+fa*A)/2-mass*A*fa*A;Bp=(B*fb+fb*B)/2-mass*B*fb*B
                Cp=-(fb*C+C*fb)/2
                Jp=Cp*chi-chi*Cp
                raw=((Bp*J*B+B*Jp*B+B*J*Bp)/2-(Bp*C*A*C*B+B*Cp*A*C*B+B*C*Ap*C*B+B*C*A*Cp*B+B*C*A*C*Bp))[index,index]
                derivative=mp.exp(-2*z*f[index])*mass*(raw-2*f[index]*glue)
                def scalar(w):
                    _,_,Ar,Br,_=objects(w)
                    return mp.exp(-2*w*f[index])*mass*(Ar[index,index]-Br[index,index])
                exact=mp.diff(scalar,z)
                de=abs(derivative-exact)/abs(exact)
                check('all_parameter_insertions_'+str(dps)+'_'+str(mass),de<mp.mpf('1e-35'))
                omitted=mp.exp(-2*z*f[index])*mass*raw
                check('density_omission_detected_'+str(dps)+'_'+str(mass),abs(omitted-exact)/abs(exact)>mp.mpf('.001'))
                contour=[]
                for count in PROTOCOL['Cauchy_nodes']:
                    radius=mp.mpf('.1')
                    value=sum(scalar(z+radius*mp.exp(2j*mp.pi*j/count))*mp.exp(-2j*mp.pi*j/count) for j in range(count))/(count*radius)
                    ce=abs(value-exact)/abs(exact);contour.append(ce)
                    check('parameter_Cauchy_'+str(dps)+'_'+str(mass)+'_'+str(count),ce<mp.mpf('1e-35'))
                errs=[]
                for dt in map(mp.mpf,PROTOCOL['finite_difference_steps']):
                    e=abs((scalar(z+dt)-scalar(z-dt))/(2*dt)-exact)/abs(exact);errs.append(e)
                check('parameter_real_difference_refinement_'+str(dps)+'_'+str(mass),errs[0]<mp.mpf('1e-8') and errs[1]<mp.mpf('.27')*errs[0])
                # All observation changes occur in evaluation rows; A itself is unchanged.
                for j in [3,5]:
                    ee=abs((A[j,j]-B[j,j])-(B*J*B/2-B*C*A*C*B)[j,j])/abs(A[j,j]-B[j,j])
                    check('fixed_middle_neighbor_observation_'+str(dps)+'_'+str(mass)+'_'+str(j),ee<mp.mpf('1e-35'))
                records.append({'dps':dps,'manufactured_mass':mass,'two_sided_relative_error':mp.nstr(err,15),'product_derivative_relative_error':mp.nstr(de,15),'Cauchy_relative_errors':[mp.nstr(e,15) for e in contour],'finite_difference_relative_errors':[mp.nstr(e,15) for e in errs]})
    FIXTURES.append({'fixture':'fixed_domain_resolvent_mechanism','records':records,'actual_quantum_state_sampled':False,'production_mass_gap_tested_by_model':False})

def main():
    sys.stdout.reconfigure(encoding='utf-8',newline='\n');sys.stderr.reconfigure(encoding='utf-8',newline='\n')
    ps=paths();hashes=lambda:{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in ps}
    before=hashes();error=None
    print(json.dumps({'event':'start','story':'S167 fixed-domain cold far seam','input_count':len(ps)}),flush=True)
    try:cone_and_energy();distance_and_sum();flat_identities();finite_domain_mechanism()
    except Exception:error=traceback.format_exc();print(error,file=sys.stderr,flush=True)
    after=hashes();good=error is None and before==after and all(r['passed'] for r in ROWS)
    run={'story':'S167 fixed-domain cold far seam','created_utc':datetime.now(timezone.utc).isoformat(),
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
