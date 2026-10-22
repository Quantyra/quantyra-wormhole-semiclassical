"""Frozen independent geometry and budget diagnostics for the mixed form."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as Q
import json,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='recentered-mixed-form';PREVIOUS='recentered-quantum-tt-form'
OWN=['notes/'+STEM+'.md','checks/recentered_mixed_form.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def exact_identities():
    u,v,x,de=S.symbols('u v x delta',positive=True)
    VV=3-12*(1+2*x+de/2)*u**3/x+13*(1+3*x)*u**4/x-S.Rational(7,3)*u*u/x
    tt=(3+1/x)*u**4;PP=3-24*u**3+39*u**4-S.Rational(1791,100)*u**8
    remainder=(-12*u**3-6*de*u**3+13*u**4-S.Rational(7,3)*u*u-S.Rational(597,50)*u**8)/x-S.Rational(199,100)*u**8/x**2
    check('exact_potential_weight_remainder',S.expand(VV-S.Rational(199,100)*tt**2-PP-remainder)==0)
    cert=[]
    for a,b in [(S.Rational(0),S.Rational(1,2)),(S.Rational(1,2),S.Rational(1))]:
        poly=S.Poly(S.expand(PP.subs(u,a+(b-a)*v)),v)
        coeff=[S.simplify(sum(poly.nth(j)*S.binomial(i,j)/S.binomial(8,j) for j in range(i+1))) for i in range(9)]
        rebuilt=sum(coeff[i]*S.binomial(8,i)*v**i*(1-v)**(8-i) for i in range(9))
        check('Bernstein_exact_reconstruction_'+str(a),S.expand(rebuilt-poly.as_expr())==0)
        check('Bernstein_uniform_lower_'+str(a),min(coeff)>=S.Rational(9,100))
        cert.append(dict(interval=[str(a),str(b)],coefficients=[str(z) for z in coeff],minimum=str(min(coeff))))
    # Reconstruct all sectional curvature invariants before varying them.
    A,B,C,E,rr,ww,ttt,lam=S.symbols('A B C E r W t lambda')
    R=2*(A+2*B+2*C+E);Ric2=(A+2*B)**2+(A+2*C)**2+2*(B+C+E)**2
    Rm2=4*(A*A+2*B*B+2*C*C+E*E)
    sub={A:-lam+2*ww-ttt,B:-lam-ww,C:-lam-ww,E:-lam+2*ww+ttt}
    check('RN_curvature_scalar_from_sections',S.expand(R.subs(sub)+12*lam)==0)
    check('RN_Weyl_square_from_sections',S.expand((Rm2-2*Ric2+R*R/3).subs(sub)-48*ww**2)==0)
    check('RN_J_square_from_sections',S.expand((Ric2-R*R/4).subs(sub)-4*ttt**2)==0)
    et,xi,q=S.symbols('e xi q');rad=rr*S.sqrt(1+et*xi)
    check('physical_magnetic_variation',S.simplify(S.diff(q*q/(2*rad**4),et).subs(et,0)+q*q*xi/rr**4)==0)
    sig,Aa,Bb,cc=S.symbols('sigma A b xi');density=S.sqrt((1+et*Aa)*(1+et*Bb))*(1+et*cc)
    check('TT_density_first_variation_zero',S.simplify(S.diff(density,et).subs({et:0,Aa:-Bb-2*cc}))==0)
    T,Y=S.Function('T'),S.Function('Y')
    check('off_shell_density_divergence_identity',S.expand(-Y(rr)*S.diff(T(rr),rr)-T(rr)*(S.diff(Y(rr),rr)+2*Y(rr)/rr)+S.diff(rr**2*T(rr)*Y(rr),rr)/rr**2)==0)
    check('nonzero_trace_return_retained',S.simplify(S.diff(rr**2*rr**-4*rr**-2,rr)/rr**2)!=0)
    FIXTURES.append(dict(fixture='uniform_potential_and_variation_identities',Bernstein=cert,physical_inverse_proved=False))

def rational_bounds():
    pl=Q('3.14159');ph=Q('3.14160');k=2000;g0=Q(1,k**10);gl=(1-Q('1e-34'))*g0;cr=Q('1e22')
    fl=1-384*ph*cr*g0;beta=Q('1.223e-10');epg=g0/(60*pl)
    Dlam=96*ph*cr*g0;Cv=1+11*g0/(180*pl)
    check('reference_Dlambda_budget',Dlam<Q('2.946e-9'))
    check('negative_potential_remainder_bound',12+6+Q(7,3)+6*Q('1.99')<33)
    check('uniform_sharp_potential_margin',Q('.09')-Q(33,8*10**9)-Q(2,(8*10**9)**2)>Q('.08'))
    eta=Q('.005')+Q('7.07e-8')
    check('sharp_TT_gradient_margin',1-eta>Q('.99'))
    check('sharp_TT_potential_margin',Q('.08')-eta>Q('.07'))
    check('sharp_energy_controls_old_energy',Q('.99')>Q('.07'))
    Kmain2=Q(2,1)/(3*Q('1.99')*fl)*(Cv**2+4*Dlam)
    check('main_mixed_weight_bound',Kmain2<Q('.33501')<Q('.579')**2)
    check('tensor_L2_embedding',22/Q('.07')<315)
    check('tensor_gradient_embedding',5200/Q('.07')<75000)
    check('return_vector_embedding',18/Q('.07')<258)
    # Explicit safe bounds for all curvature-variation coefficients.
    check('Riemann_from_sectional_bounds',4*(Q('14.1')**2+4*Q('6.1')**2+Q('8.1')**2)<60**2)
    check('Weyl_Frobenius_bound',48*(5+2*beta)**2<35**2)
    check('J_Frobenius_bound',2*(3+beta)<7)
    check('Schouten_bound',(13+4)/2<9)
    check('W_gradient_bound',5*(3*(2+beta)+4*(3+beta))**2<41**2)
    check('J_gradient_bound',144*5*(3+beta)**2<100**2)
    # sqrt(48)<7, sqrt(10)<3.163; derivative coefficient <=48 sqrt(10)(5+2beta).
    check('Weyl_gradient_bound',7*41+48*Q('3.163')*(5+2*beta)<10**4 and Q('3.163')**2>10)
    check('Weyl_variation_algebraic_coefficient',2*35*60+8*35*9+4*35**2==11620)
    check('J_variation_algebraic_coefficient',7*(2*60+2*13+6)+2*7**2<2000)
    check('combined_anomaly_coefficient',Q(7,12)+Q(11,6)<3)
    check('square_root_upper_reserves',112**2>Q(75000,6) and 56**2>Q(75000,24) and 4**2>Q(315,24) and fl>Q('.999')**2)
    Kcurv=3*10**5*epg*(112+56+4)/Q('.999')
    check('curvature_anomaly_mixed_budget',Kcurv<Q('3e-28'))
    # F_lambda<=1; use the conservative rational square-root lower reserve 6.8.
    check('magnetic_denominator_lower',Q('6.8')**2<24*Q('1.99'))
    Kem=20*Q(1,k**8)/(60*pl**2*Q('6.8'))
    check('magnetic_mixed_budget',Kem<Q('3e-29'))
    U=Q(1,k**8)*(3+beta)/(pl*gl)
    Tbound=epg*(10*U+28*(2+beta)**2+56*(2+beta)*(3+beta)+Q(106,3)*(3+beta)**2)
    check('actual_trace_value_budget',Tbound<Q('2.1e-28'))
    check('return_square_root_guard',Q('6.6')**2>258/(6*fl))
    Kreturn=Q('2.1e-28')*Q('6.6')
    check('off_shell_return_budget',Kreturn<Q('1.5e-27'))
    check('complete_small_remainder',Kcurv+Kem+Kreturn<Q('2e-27'))
    check('relative_mixed_target',(Q('.579')+Q('2e-27'))**2<Q('.336')<Q(3,8))
    check('conditional_Schur_loop',Q(4,3)*Q('.336')==Q('.448')<Q(1,2))
    check('conditional_scalar_inverse',Q(4,3)/(1-Q('.448'))<Q('2.42'))
    FIXTURES.append(dict(fixture='physical_mixed_budget',main_squared_upper=str(Kmain2),curvature_upper=str(Kcurv),magnetic_upper=str(Kem),trace_value_upper=str(Tbound),return_upper=str(Kreturn),relative_form_bound='.336',inverse_product_conditional=True))

def metric_variations():
    records=[]
    for prec in [70,100]:
        mp.mp.dps=prec;kk=mp.mpf(2000);q=kk**4;g0=kk**-10
        G=2*g0/(1+mp.sqrt(1+4*11*g0/(360*mp.pi)));lam=mp.mpf(1);lam0=g0/G
        cr=mp.mpf('1e22');F=1-384*mp.pi*cr*G;ep=G/(60*mp.pi);DD=96*mp.pi*G*cr-ep;Cv=1+11*ep/3
        C=mp.sqrt(1+12*mp.pi*G*kk**16/F);B=mp.sqrt((C-1)/6);Z=B**2+3*B**4
        C0=mp.sqrt(1+12*mp.pi*kk**6);b0=mp.sqrt((C0-1)/6*G/g0)
        kap=mp.mpf('.265')*2*q/kk*G*C0/b0**3;gam=kap*B/C;x=B**2
        K=lambda h:((1+3*x)*((1+h)**-2+(1+h)**-3)+3*x*(1+(1+h)**-1))/(2*C)
        hh=gam
        for _ in range(5):hh-=(hh*K(hh)-gam)/(K(hh)+hh*mp.diff(K,hh))
        rp=B*(1+hh);M=(rp+Z/rp+rp**3)/2
        f=lambda r:1-2*M/r+Z/r**2+r*r;xi=lambda r:(B/r)**3
        def num(r):
            return B**3*(6*M*(r**-2-rp**-2)/(-2)-4*Z*(r**-3-rp**-3)/(-3)-2*(r**-1-rp**-1)/(-1))
        bb=lambda r:-num(r)/(r*r*f(r));aa=lambda r:-bb(r)-2*xi(r)
        def sections(r,e):
            lapse=lambda y:mp.sqrt(f(y)*(1+e*aa(y)))
            radial=lambda y:mp.sqrt((1+e*bb(y))/f(y))
            sphere=lambda y:y*mp.sqrt(1+e*xi(y))
            n=lapse(r);m=radial(r);s=sphere(r)
            np=mp.diff(lapse,r);sp=mp.diff(sphere,r);mp1=mp.diff(radial,r)
            return [-(mp.diff(lapse,r,2)-mp1/m*np)/(m*m*n),-(mp.diff(sphere,r,2)-mp1/m*sp)/(m*m*s),-np*sp/(m*m*n*s),(1-(sp/m)**2)/s**2]
        def invariants(r,e):
            A,Bb,Cc,E=sections(r,e);R=2*(A+2*Bb+2*Cc+E)
            Ric2=(A+2*Bb)**2+(A+2*Cc)**2+2*(Bb+Cc+E)**2
            Rm2=4*(A*A+2*Bb*Bb+2*Cc*Cc+E*E)
            return R,Rm2-2*Ric2+R*R/3,Ric2-R*R/4,q*q/(2*r**4*(1+e*xi(r))**2)
        def trace(r,e):
            RR,CC,JJ,II=invariants(r,e)
            dens=lambda y:y*y*mp.sqrt((1+e*aa(y))*(1+e*bb(y)))*(1+e*xi(y))
            coef=lambda y:dens(y)*f(y)/(1+e*bb(y))
            lap=mp.diff(lambda y:coef(y)*mp.diff(lambda z:invariants(z,e)[0],y),r)/dens(r)
            return RR+12*lam0-DD*lap+ep*(-mp.mpf(11)/72*RR**2+mp.mpf(7)/12*CC+mp.mpf(11)/6*JJ-20*II)
        for ratio in ['1.3','2','7']:
            r=B*mp.mpf(ratio);p=mp.diff(f,r)/(2*f(r));fr=f(r);bv=bb(r);xv=xi(r)
            base=sections(r,0)
            ds=[-bv*base[0]-fr/2*(mp.diff(aa,r,2)+3*p*mp.diff(aa,r)-p*mp.diff(bb,r)),
                -bv*base[1]-fr/2*(mp.diff(xi,r,2)+(2/r+p)*mp.diff(xi,r)-mp.diff(bb,r)/r),
                -bv*base[2]-fr/2*(mp.diff(aa,r)/r+p*mp.diff(xi,r)),
                -xv/r**2+fr*bv/r**2-fr*mp.diff(xi,r)/r]
            for i in range(4):check('section_variation_'+ratio+'_'+str(i)+'_'+str(prec),abs(mp.diff(lambda e:sections(r,e)[i],0)-ds[i])<mp.mpf('1e-45'))
            dR=2*(ds[0]+2*ds[1]+2*ds[2]+ds[3]);directR=-4*Z/r**4*xv
            check('actual_TT_curvature_variation_'+ratio+'_'+str(prec),abs(dR/directR-1)<mp.mpf('1e-45'))
            A,Bb,Cc,E=base;da,db,dc,de=ds;R0=-12
            ric=[A+2*Bb,A+2*Cc,Bb+Cc+E];dric=[da+2*db,da+2*dc,db+dc+de]
            dRic2=2*(ric[0]*dric[0]+ric[1]*dric[1]+2*ric[2]*dric[2]);dRm2=8*(A*da+2*Bb*db+2*Cc*dc+E*de)
            dJ2=dRic2-R0*dR/2;dC2=dRm2-2*dRic2+2*R0*dR/3;dI=-q*q*xv/r**4
            rdot=lambda y:-4*Z/y**4*xi(y)
            lapdot=f(r)*mp.diff(rdot,r,2)+(mp.diff(f,r)+2*f(r)/r)*mp.diff(rdot,r)
            expected=Cv*directR-DD*lapdot+ep*(mp.mpf(7)/12*dC2+mp.mpf(11)/6*dJ2-20*dI)
            density=lambda e:mp.sqrt((1+e*aa(r))*(1+e*bb(r)))*(1+e*xi(r))
            literal=lambda e:density(e)*trace(r,e)
            deriv=mp.diff(literal,0);err=abs(deriv/expected-1)
            check('unexpanded_full_mixed_trace_'+ratio+'_'+str(prec),err<mp.mpf('1e-45'))
            errors=[]
            for ss in ['1e-4','5e-5']:
                step=mp.mpf(ss);fd=(literal(step)-literal(-step))/(2*step);errors.append(abs(fd/expected-1))
            rate=errors[1]/errors[0]
            check('full_metric_difference_error_'+ratio+'_'+str(prec),errors[1]<mp.mpf('1e-6'))
            check('full_metric_quadratic_refinement_'+ratio+'_'+str(prec),mp.mpf('.24')<rate<mp.mpf('.26'))
            records.append(dict(digits=prec,radius_ratio=ratio,full_mixed_row=mp.nstr(expected,60),direct_derivative_error=mp.nstr(err,20),finite_difference_errors=[mp.nstr(z,25) for z in errors],refinement=mp.nstr(rate,20),magnetic_piece=mp.nstr(-20*ep*dI,50),curvature_anomaly_piece=mp.nstr(ep*(mp.mpf(7)/12*dC2+mp.mpf(11)/6*dJ2),50)))
    for ratio in ['1.3','2','7']:
        pair=[r for r in records if r['radius_ratio']==ratio]
        check('mixed_cross_precision_'+ratio,abs(mp.mpf(pair[0]['full_mixed_row'])/mp.mpf(pair[1]['full_mixed_row'])-1)<mp.mpf('1e-45'))
    FIXTURES.append(dict(fixture='unexpanded_actual_RN_metric_and_full_trace_variations',records=records,curved_quantum_kernel_sampled=False,physical_boundary_realization_proved=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_quantum_TT_result',old['status']=='passed' and old['quantum_TT_reference_lower_allowance']=='.005 E_eff')
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        exact_identities();rational_bounds();metric_variations()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        relative_mixed_form_bound='.336' if error is None else None,
        weighted_TT_potential_bound='V>=1.99 t^2/lambda+.08lambda' if error is None else None,
        scope='reference scalar and returned-TT forms on stated cores/completions; inverse product conditional on physical realization',
        physical_boundary_completed=False,actual_form_closedness_in_own_norm_proved=False,physical_coupled_inverse_proved=False,
        coupled_correction_norm_proved=False,nonlinear_solution_constructed=False,semiclassical_validity_proved=False,
        physical_stability_proved=False,component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
