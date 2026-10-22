"""Recentered scalar trace response and paired boundary normal jets; supporting checks."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='recentered-r2-scalar-boundary'
PREVIOUS='charged-r2-trace-block'
OWN=['notes/'+STEM+'.md','checks/recentered_r2_scalar_boundary.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def boundary_algebra():
    I=S.eye(2);p1=S.Matrix([[0,1],[1,0]]);p2=S.Matrix([[0,-S.I],[S.I,0]]);p3=S.diag(1,-1)
    tensor=S.kronecker_product;eye=S.eye(8)
    normal=tensor(p3,p1,I)
    tangent=[-tensor(I,p3,I),-tensor(I,p2,p1),-tensor(I,p2,p2)]
    gam=[normal]+tangent;gen=S.I*tensor(p2,p1,I)
    ls=S.symbols('l0 l1 l2 l01 l02 l12',real=True)
    L=S.Matrix([[ls[0],ls[3],ls[4]],[ls[3],ls[1],ls[5]],[ls[4],ls[5],ls[2]]])
    tr=S.trace(L);l2=S.trace(L*L)
    results=[]
    for sign in [1,-1]:
        chi=sign*normal;P=(eye+chi)/2;Pm=eye-P;robin=-tr*P/2
        dchi=[-sign*sum((L[a,b]*tangent[b] for b in range(3)),S.zeros(8)) for a in range(3)]
        sq=sum((d*d for d in dchi),S.zeros(8))
        check('paired_rank_'+str(sign),S.trace(P)==4 and P*P==P and S.trace(chi)==0)
        induced=sum((P*normal*tangent[a]*dchi[a]*P/2 for a in range(3)),S.zeros(8))
        check('induced_Robin_'+str(sign),all(S.simplify(x)==0 for x in induced-robin))
        basis=[eye]+[gam[i]*gam[j] for i in range(4) for j in range(i+1,4)]
        check('curvature_gauge_trace_cancellation_'+str(sign),all(S.trace(chi*E)==0 and gen*E==E*gen for E in basis))
        check('projector_jet_trace_'+str(sign),S.simplify(S.trace(sq)-8*l2)==0)
        fn=S.trace(((84*P-180*Pm)*tr**2+(84*P+60*Pm)*l2)/7+72*robin*tr+240*robin**2-18*sq)/360
        fnn=S.trace(24*tr*eye+120*robin)/360
        check('complete_first_normal_heat_row_'+str(sign),S.simplify(fn-8*(2*tr**2-3*l2)/140)==0)
        check('complete_second_normal_heat_row_'+str(sign),S.simplify(fnn+8*tr/60)==0)
        check('third_normal_heat_row_zero_'+str(sign),30*S.trace(chi)==0)
        restriction=S.simplify(fn+tr*fnn)
        check('finite_Navier_boundary_obstruction_'+str(sign),S.simplify(restriction+8*(tr**2+9*l2)/420)==0)
        check('umbilic_obstruction_nonzero_'+str(sign),restriction.subs(dict(zip(ls,[1,1,1,0,0,0])))==-S.Rational(24,35))
        results.append(dict(domain='Dirac' if sign==1 else 'adjoint',normal_1=str(S.factor(fn)),normal_2=str(S.factor(fnn)),restricted=str(S.factor(restriction))))
    check('partner_average_no_extra_copy',results[0]['normal_1']==results[1]['normal_1'] and results[0]['normal_2']==results[1]['normal_2'])
    FIXTURES.append(dict(fixture='complete_rank_eight_normal_heat_rows',domains=results,angle_derivatives=False,full_physical_boundary_proved=False))

def recenter_algebra():
    r,M,Z,l,ep,D,q,rp=S.symbols('r M Z l ep D q rp',positive=True)
    f=1-2*M/r+Z/r**2+l*r*r
    R=-S.diff(f,r,2)-4*S.diff(f,r)/r+2*(1-f)/r**2
    check('recentered_RN_scalar',S.simplify(R+12*l)==0)
    t=Z/r**4;w=M/r**3-Z/r**4
    H=28*w*w+S.Rational(22,3)*t*t-10*q*q/r**4
    l0=l+S.Rational(11,6)*ep*l*l
    trace=R+12*l0+ep*(H-S.Rational(11,72)*R**2)
    check('exact_vacuum_trace_cancellation',S.simplify(trace-ep*H)==0)
    co=[-10*q*q,28*M*M,-56*M*Z,S.Rational(106,3)*Z*Z];powers=[4,6,7,8]
    check('actual_trace_source_polynomial',S.simplify(H-sum(a/r**p for a,p in zip(co,powers)))==0)
    # Fixed density derivative derived before factorization.
    eta,s,ds,d2s,jh=S.symbols('eta s ds d2s jh')
    EH=S.exp(2*eta*s)*(-12*l-6*eta*ds)+12*l0*S.exp(4*eta*s)
    anomaly=-S.Rational(22,3)*eta*jh-S.Rational(11,72)*(-12*l-6*eta*ds)**2
    deriv=S.diff(EH+ep*anomaly,eta).subs(eta,0)-D*(24*l*ds-6*d2s)
    Cv=1+S.Rational(11,3)*ep*l
    expected=6*D*d2s-(6*Cv+24*D*l)*ds+24*l*Cv*s-S.Rational(22,3)*ep*jh
    check('recentered_density_operator',S.simplify(deriv-expected)==0)
    z=S.symbols('z')
    check('exact_fourth_order_factorization',S.expand(6*(-z+4*l)*(-D*z+Cv)-(6*D*z*z-(6*Cv+24*D*l)*z+24*l*Cv))==0)
    gram=sum(a*b*rp**(3-p-v)/(p+v-3) for a,p in zip(co,powers) for b,v in zip(co,powers))
    y=S.symbols('y',positive=True)
    poly=S.Poly(S.expand(H.subs(r,1/y)**2),y)
    direct=sum(a*rp**(3-power[0])/(power[0]-3) for power,a in poly.terms())
    check('signed_exact_source_L2_integral',S.simplify(gram-direct)==0)
    p=S.symbols('p');P=6*l*(4-p*(p-3))*(Cv-D*l*p*(p-3))
    check('energy_roots_minus_one_four',P.subs(p,-1)==0 and P.subs(p,4)==0)
    check('resonant_derivative',S.simplify(-S.diff(P,p).subs(p,4)-30*l*(Cv-4*D*l))==0)
    a=ep*q*q/(3*l*(Cv-4*D*l))
    check('forced_log_coefficient',S.simplify(-S.diff(P,p).subs(p,4)*a-10*ep*q*q)==0)
    # Independently apply the normal radial differential operator to r^-4 log r.
    lap=lambda fun:l*(r*r*S.diff(fun,r,2)+4*r*S.diff(fun,r))
    tail=r**-4*S.log(r)
    actual=6*(-lap(-D*lap(tail)+Cv*tail)+4*l*(-D*lap(tail)+Cv*tail))
    check('direct_normal_log_substitution',S.simplify(actual-30*l*(Cv-4*D*l)/r**4)==0)
    # Boundary normal rows on the leading AdS geometry, independent of RN expansion.
    fp=l*r*r;sn=-S.sqrt(fp)*S.diff(tail,r)
    snn=fp*S.diff(tail,r,2)+S.diff(fp,r)*S.diff(tail,r)/2
    LL=3*S.sqrt(l);LL2=3*l
    normal_density=S.sqrt(fp)*r*r*((2*LL*LL-3*LL2)*sn/140-LL*snn/60)
    check('normal_boundary_tail_asymptotic',S.simplify(normal_density-l*l/r*(-S.Rational(19,35)*S.log(r)+S.Rational(47,140)))==0)
    FIXTURES.append(dict(fixture='exact_recentered_source_and_tail',source_powers=powers,source_integral=str(S.factor(gram)),leading_log_coefficient=str(a),full_coupled_correction=False))

def rational_bounds():
    pl=Q('3.14159');ph=Q('3.14160');k=2000;g0=Q(1,k**10);gl=(1-Q('1e-34'))*g0;c=Q('1e22')
    check('vacuum_root_interval',gl+11*gl*gl/(360*pl)<g0 and g0+11*g0*g0/(360*ph)>g0)
    flo=1-384*ph*c*g0;Clo=Q(49119000000)
    check('recentered_positive_F',flo>1-Q('1.179e-8'))
    check('recentered_C_lower',1+12*pl*gl*k**16>Clo*Clo)
    beta=Q('1.223e-10')
    check('recentered_inverse_radius',6/(Clo-1)<beta)
    Dl=(96*pl*c-1/(60*pl))*gl;Dh=(96*ph*c-1/(60*ph))*g0
    epg=g0/(60*pl)
    check('recentered_gradient_coercivity',24*Dl-Q(22,3)*epg*beta>Q('7.06e-8'))
    denom=1-4*Dh
    check('separated_positive_indicial_factors',denom>Q('.99999998') and Dl>0)
    check('fast_root_above_18000',1/Dh>18000*17997)
    roots={5:Q('2.23606'),11:Q('3.31662'),13:Q('3.60555'),24:Q('4.89897')}
    check('rational_square_root_lower_bounds',all(v*v<n for n,v in roots.items()))
    U=Q(1,k**8)*(3+beta)/(pl*gl)
    bound=epg/roots[24]*(10*U/roots[5]+28*(2+beta)**2/3+56*(2+beta)*(3+beta)/roots[11]+Q(106,3)*(3+beta)**2/roots[13])
    lead=epg*U/(3*denom)
    check('actual_trace_response_energy_budget',bound<Q('1.82e-29'))
    check('actual_forced_log_budget',lead<Q('6.60e-30'))
    FIXTURES.append(dict(fixture='rational_uniform_response_bounds',g_interval=[str(gl),str(g0)],D_lambda=[str(Dl),str(Dh)],U_upper=str(U),energy_upper=str(bound),leading_log_upper=str(lead),metric_admission_proved=False))

def numerical():
    records=[]
    for prec in [60,90]:
        mp.mp.dps=prec;pi=mp.pi;k=mp.mpf(2000);g0=k**-10
        g=2*g0/(1+mp.sqrt(1+4*(11/(360*pi))*g0))
        # Use lambda=1 units. The dimensionless couplings and physical ratios remain fixed.
        G=g;lam=mp.mpf(1);c=mp.mpf('1e22');q=k**4;gr=1/q
        ep=G/(60*pi);D=96*pi*G*c-ep;Cv=1+mp.mpf(11)/3*ep
        F=1-384*pi*G*c;Z=pi*G*q*q/(gr*gr*F);B=mp.sqrt((mp.sqrt(1+12*Z)-1)/6)
        for ratio in ['1.000001','1.05','2']:
            rp=B*mp.mpf(ratio);M=(rp+Z/rp+rp**3)/2
            powers=[4,6,7,8];co=[-10*q*q,28*M*M,-56*M*Z,mp.mpf(106)/3*Z*Z]
            bp=[a/rp**p for a,p in zip(co,powers)]
            gram=sum(a*b/(p+v-3) for a,p in zip(bp,powers) for b,v in zip(bp,powers))
            integral=mp.quad(lambda u:sum(a*u**(p-2) for a,p in zip(bp,powers))**2,[0,1])
            check('source_integral_'+ratio+'_'+str(prec),abs(gram/integral-1)<mp.mpf('1e-45'))
            nb=ep*mp.sqrt(gram)/mp.sqrt(24);ab=ep*q*q/(3*(Cv-4*D)*rp**4)
            check('numeric_response_budget_'+ratio+'_'+str(prec),nb<mp.mpf('1.82e-29') and ab<mp.mpf('6.60e-30'))
            f=lambda r:1-2*M/r+Z/r**2+r*r
            # Complete normal-smearing polynomial is also evaluated using literal Robin matrices below.
            cuts=[]
            for scale in ['2','10','1000000']:
                r=rp*mp.mpf(scale);fr=f(r);df=mp.diff(f,r)
                Lt=df/(2*mp.sqrt(fr));Ls=mp.sqrt(fr)/r;L=Lt+2*Ls;L2=Lt*Lt+2*Ls*Ls
                residual=-(L*L+9*L2)/420
                check('finite_cut_term_nonzero_'+ratio+'_'+scale+'_'+str(prec),residual<0)
                # Tail fixture is r_+^4 r^-4 log(r/r_+), not an integrated coupled solution.
                tail=lambda rr:(rp/rr)**4*mp.log(rr/rp)
                sn=-mp.sqrt(fr)*mp.diff(tail,r)
                snn=fr*mp.diff(tail,r,2)+df*mp.diff(tail,r)/2
                normal=mp.sqrt(fr)*r*r*((2*L*L-3*L2)*sn/140-L*snn/60)
                scaled=normal*r/rp**4
                leading=-mp.mpf(19)/35*mp.log(r/rp)+mp.mpf(47)/140
                if scale=='1000000':check('normal_tail_limit_'+ratio+'_'+str(prec),abs(scaled-leading)<mp.mpf('1e-12'))
                cuts.append(dict(scale=scale,normal_term=mp.nstr(normal,30),finite_Navier_coefficient_per_rank_without_4pi=mp.nstr(residual,30),tail_asymptotic_error=mp.nstr(abs(scaled-leading),10)))
            # Evaluate the ACTUAL radial operator on its forced leading profile, checking its leading residual.
            lap=lambda fun,r:f(r)*mp.diff(fun,r,2)+(mp.diff(f,r)+2*f(r)/r)*mp.diff(fun,r)
            lead=lambda rr:ab*(rp/rr)**4*mp.log(rr/rp)
            r=rp*mp.mpf('1e6');lp=lambda rr:lap(lead,rr)
            sp=mp.diff(lead,r);spp=mp.diff(lead,r,2);t=Z/r**4
            jh=-t*(f(r)*spp+mp.diff(f,r)*sp)+2*t*f(r)*sp/r
            value=6*D*lap(lp,r)-(6*Cv+24*D)*lp(r)+24*Cv*lead(r)-mp.mpf(22)/3*ep*jh
            leading_source=10*ep*q*q/r**4
            check('actual_radial_log_residual_'+ratio+'_'+str(prec),abs(value/leading_source-1)<mp.mpf('1e-12'))
            records.append(dict(digits=prec,radius_ratio=ratio,g=mp.nstr(g,50),energy_bound=mp.nstr(nb,50),forced_log_amplitude=mp.nstr(ab,50),
                fast_root=mp.nstr(mp.mpf('1.5')+mp.sqrt(mp.mpf('2.25')+Cv/D),50),source_integral_relative_error=mp.nstr(abs(gram/integral-1),10),
                forced_tail_radial_error=mp.nstr(abs(value/leading_source-1),10),cuts=cuts))
    for ratio in ['1.000001','1.05','2']:
        pair=[r for r in records if r['radius_ratio']==ratio]
        for key in ['g','energy_bound','forced_log_amplitude','fast_root']:
            check('precision_'+ratio+'_'+key,abs(mp.mpf(pair[0][key])/mp.mpf(pair[1][key])-1)<mp.mpf('1e-45'))
    FIXTURES.append(dict(fixture='actual_trace_source_integrals_and_boundary_geometry',records=records,quantum_state_samples=False,full_scalar_profile_numerically_solved=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',
        protocol=dict(digits=[60,90],relative_tolerance='1e-45',asymptotic_fixture_tolerance='1e-12',continuum_machine_verified=False))
    error=None
    try:
        check('previous_charged_comparison',old['status']=='passed' and old['local_charged_tensor_lift'] and old['bulk_comparison_inverse_bound']=='1')
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        boundary_algebra();recenter_algebra();rational_bounds();numerical()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,
        input_sha256_after=after,unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,
        finished_utc=datetime.now(timezone.utc).isoformat(),recentered_scalar_source_admissible=error is None,
        scalar_comparison_response_energy_bound='1.82e-29' if error is None else None,
        finite_Navier_Dirac_boundary_term_nonzero=error is None,
        scalar_response_AdS_boundary_limit_zero=error is None,
        physical_boundary_completed=False,physical_coupled_inverse_proved=False,
        coupled_correction_norm_proved=False,nonlinear_solution_constructed=False,
        semiclassical_validity_proved=False,physical_stability_proved=False,
        component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:
            log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
