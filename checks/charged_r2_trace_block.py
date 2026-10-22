"""Charged R2 reference and bulk trace block: supporting checks, not a wormhole."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json, traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='charged-r2-trace-block'
PREVIOUS='finite-r2-scalar-window'
OWN=['notes/'+STEM+'.md','checks/charged_r2_trace_block.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]

def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def algebra():
    r,M,Z,l,G,c,eps,D=S.symbols('r M Z l G c eps D',positive=True)
    f=1-2*M/r+Z/r**2+l*r**2
    # Independent orthonormal sectional contractions on the Euclidean metric.
    sect=[-S.diff(f,r,2)/2,-S.diff(f,r)/(2*r),-S.diff(f,r)/(2*r),
          -S.diff(f,r)/(2*r),-S.diff(f,r)/(2*r),(1-f)/r**2]
    pairs=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)]
    ric=[S.simplify(sum(v for v,p in zip(sect,pairs) if i in p)) for i in range(4)]
    R=S.simplify(sum(ric));J=[v-R/4 for v in ric]
    j2=S.simplify(sum(v*v for v in J))
    riem2=4*sum(v*v for v in sect)
    c2=S.simplify(riem2-2*sum(v*v for v in ric)+R*R/3)
    check('RN_scalar_curvature',R==-12*l)
    check('RN_traceless_Ricci',all(S.simplify(v-w)==0 for v,w in zip(J,[-Z/r**4]*2+[Z/r**4]*2)))
    check('RN_curvature_contractions',S.simplify(j2-4*Z**2/r**8)==0 and S.simplify(c2-48*(M/r**3-Z/r**4)**2)==0)
    F0=1-384*S.pi*G*c*l
    # Local tensor equation: algebraic Ricci and Maxwell eigenvalues.
    for i in range(4):
        em=F0*J[i]
        lhs=F0*(ric[i]-R/2)-3*l+8*S.pi*G*c*R**2
        check('local_tensor_component_'+str(i),S.simplify(lhs-em)==0)
    a=16*S.pi*G*c/F0;m0=1/(96*S.pi*G*c)
    check('charged_lift_curvature',S.simplify(6*a*(m0-4*l)-1)==0)
    check('charged_lift_derivative_tensor',S.simplify(-2*F0+32*S.pi*G*c/a)==0)
    check('charged_lift_algebraic_tensor',S.simplify(-12*l+32*S.pi*G*c*R**2-R*F0)==0)
    phi,RR=S.symbols('phi RR')
    check('auxiliary_Ward_identity',S.simplify(-RR*32*S.pi*G*c/2+16*S.pi*G*c*phi-16*S.pi*G*c*(phi-RR))==0)
    b=S.symbols('b',positive=True)
    me=(b+Z/b+l*b**3)/2
    temp=S.diff(f,r).subs({M:me,r:b})/(4*S.pi)
    expected=(1/b-Z/b**3+3*l*b)/(4*S.pi)
    check('temperature_from_bolt',S.simplify(temp-expected)==0)
    check('temperature_derivative_lower_bound',S.simplify(4*S.pi*S.diff(expected,b)-(3*l-1/(12*Z))-3*Z*(1/b**2-1/(6*Z))**2)==0)
    check('extremal_charge_relation',S.simplify(expected.subs(Z,b*b+3*l*b**4))==0)
    # Density differentiation, retaining independent Hessian and trace symbols.
    eta,s,lap,lap2,jh,H=S.symbols('eta s lap lap2 jh H')
    eh=S.exp(2*eta*s)*(R-6*eta*lap)+12*l*S.exp(4*eta*s)
    an=H-eta*S.Rational(22,3)*jh-S.Rational(11,72)*(R-6*eta*lap)**2
    boxvar=-2*R*lap-6*lap2
    dens=S.diff(eh+eps*an,eta).subs(eta,0)-D*boxvar
    target=6*D*lap2-(6+24*D*l+22*eps*l)*lap-S.Rational(22,3)*eps*jh+24*l*s
    check('full_density_trace_operator',S.simplify(dens-target)==0)
    t0=eps*(H-22*l*l)
    dR=-6*lap+24*l*s
    undens=(1+S.Rational(11,3)*eps*l)*dR-D*(24*l*lap-6*lap2)+eps*(-4*s*H-S.Rational(22,3)*jh)
    check('density_contact_cancels_Weyl_gauge_potential',S.simplify(undens+4*s*t0-target)==0)
    u=S.Function('u')(r);t=Z/r**4
    jh_rad=-t*(f*S.diff(u,r,2)+S.diff(f,r)*S.diff(u,r))+2*t*f*S.diff(u,r)/r
    check('traceless_Ricci_divergence_identity',S.simplify(jh_rad-S.diff(-r*r*t*f*S.diff(u,r),r)/r**2)==0)
    kr=6+24*D*l+22*eps*l-S.Rational(22,3)*eps*t
    check('sharp_gradient_lower_bound',S.simplify(kr.subs(t,3*l+1/b**2)-(6+24*D*l-S.Rational(22,3)*eps/b**2))==0)
    # Green flux identity is checked before any quadrature.
    lap_u=S.diff(r*r*f*S.diff(u,r),r)/r**2
    Au=6*D*S.diff(r*r*f*S.diff(lap_u,r),r)/r**2-S.diff(r*r*f*kr*S.diff(u,r),r)/r**2+24*l*u
    energy=6*D*lap_u**2+kr*f*S.diff(u,r)**2+24*l*u**2
    flux=6*D*r*r*f*(u*S.diff(lap_u,r)-S.diff(u,r)*lap_u)-r*r*f*kr*u*S.diff(u,r)
    check('charged_radial_Green_identity',S.simplify(r*r*(u*Au-energy)-S.diff(flux,r))==0)
    FIXTURES.append(dict(fixture='exact_charged_geometry_lift_and_density',bulk_trace_only=True,physical_boundary_proved=False))

def rational_bounds():
    pl=Q('3.14159');ph=Q('3.14160');g=Q(1,2000**10);c=Q('1e22')
    fmin=1-384*ph*c*g;fmax=1-384*pl*c*g
    check('F0_rational_margin',1-Q('1.179e-8')<fmin<fmax<1)
    # Cc exceeds the old C0, whose strict lower bound is checked directly.
    Clo=Q(49119000000)
    check('C0_lower_bound',Clo*Clo<1+12*pl*2000**6)
    beta=6/(Clo-1) # Bc^-2/lambda0 upper bound
    check('inverse_radius_bound',beta<Q('1.223e-10'))
    Dl=(96*pl*c-1/(60*pl))*g;Dh=(96*ph*c-1/(60*ph))*g
    el=g/(60*ph);eh=g/(60*pl)
    gainlo=24*Dl-Q(22,3)*eh*beta;gainhi=24*Dh
    check('coercive_gain_interval',Q('7.06e-8')<gainlo and gainhi<Q('7.08e-8'))
    check('positive_fourth_order_coefficient',Dl>0)
    # lambda0 Qc^2 = pi k^6/F0, so the bound below is conservative.
    qdimlo=pl*2000**6
    check('temperature_strictly_increasing',3-1/(12*qdimlo)>0)
    # Geometric reference shifts, not the quantum product correction.
    check('uniform_curvature_t_bound',3+beta<Q('3.000000001'))
    check('uniform_curvature_W_bound',5+2*beta<Q('5.000000001'))
    check('Schur_design_budget',Q(4,3)*Q(3,8)==Q(1,2))
    FIXTURES.append(dict(fixture='rational_benchmark_margins',F0=[str(fmin),str(fmax)],
        B_inverse_squared_over_lambda_upper=str(beta),D_lambda=[str(Dl),str(Dh)],
        K_gain=[str(gainlo),str(gainhi)],inverse_bound='1 in the fixed E norm',
        physical_domain_transferred=False))

def numerical():
    records=[]
    for prec in [60,90]:
        mp.mp.dps=prec;pi=mp.pi
        for label in ['diagnostic','benchmark']:
            lam=mp.mpf(1)
            if label=='diagnostic':
                G=mp.mpf('.002');c=mp.mpf('.01');q=mp.mpf(3);gr=mp.mpf('.7')
            else:
                G=mp.mpf(2000)**-10;c=mp.mpf('1e22');q=mp.mpf(2000)**4;gr=1/q
            F0=1-384*pi*G*c*lam;Z=pi*G*q*q/(gr*gr*F0)
            Cc=mp.sqrt(1+12*lam*Z);B=mp.sqrt((Cc-1)/(6*lam))
            rp=B*mp.mpf('1.05');M=(rp+Z/rp+lam*rp**3)/2
            ep=G/(60*pi);D=96*pi*G*c-ep
            f=lambda r:1-2*M/r+Z/r**2+lam*r*r
            lap=lambda fun,r:f(r)*mp.diff(fun,r,2)+(mp.diff(f,r)+2*f(r)/r)*mp.diff(fun,r)
            prof=lambda r:mp.exp(-(r/B-mp.mpf('1.3'))**2)
            def density(eta,r):
                s=lambda x:eta*prof(x)
                R=lambda x:mp.exp(-2*s(x))*(-12*lam-6*lap(s,x)-6*f(x)*mp.diff(s,x)**2)
                sp=mp.diff(s,r);spp=mp.diff(s,r,2);hess=[mp.diff(f,r)*sp/2,f(r)*spp+mp.diff(f,r)*sp/2,f(r)*sp/r,f(r)*sp/r]
                dl=sum(hess);vv=f(r)*sp*sp;t=Z/r**4
                jc=[(-t if i<2 else t)-2*hess[i]+dl/2+(2*vv if i==1 else 0)-vv/2 for i in range(4)]
                rt=-12*lam-6*dl-6*vv;weyl2=48*(M/r**3-Z/r**4)**2;I=q*q/(2*r**4)
                return mp.exp(2*s(r))*rt+12*lam*mp.exp(4*s(r))+ep*(7*weyl2/12+mp.mpf(11)/6*sum(v*v for v in jc)-mp.mpf(11)/72*rt*rt-20*I)-D*mp.exp(2*s(r))*(lap(R,r)+2*f(r)*sp*mp.diff(R,r))
            errors=[]
            for scale in ['1.15','1.5','2.4']:
                r=B*mp.mpf(scale);lp=lambda x:lap(prof,x)
                sp=mp.diff(prof,r);spp=mp.diff(prof,r,2);t=Z/r**4
                jh=-t*(f(r)*spp+mp.diff(f,r)*sp)+2*t*f(r)*sp/r
                expected=6*D*lap(lp,r)-(6+24*D*lam+22*ep*lam)*lp(r)-mp.mpf(22)/3*ep*jh+24*lam*prof(r)
                direct=mp.diff(lambda eta:density(eta,r),mp.mpf(0))
                errors.append(abs(direct-expected)/max(1,abs(expected)))
            check('unexpanded_charged_density_'+label+'_'+str(prec),max(errors)<mp.mpf('1e-45'))
            # Direct radial form quadrature with compact endpoint jets.
            lo=B*mp.mpf('1.1');hi=B*mp.mpf('1.7')
            u=lambda r:((r-lo)*(hi-r)/(hi-lo)**2)**3
            lu=lambda r:lap(u,r)
            kr=lambda r:6+24*D*lam+22*ep*lam-mp.mpf(22)/3*ep*Z/r**4
            A=lambda r:6*D*lap(lu,r)-mp.diff(lambda x:x*x*f(x)*kr(x)*mp.diff(u,x),r)/(r*r)+24*lam*u(r)
            v1=mp.quad(lambda r:r*r*u(r)*A(r),[lo,hi])
            v2=mp.quad(lambda r:r*r*(6*D*lu(r)**2+kr(r)*f(r)*mp.diff(u,r)**2+24*lam*u(r)**2),[lo,hi])
            e0=mp.quad(lambda r:r*r*(6*D*lu(r)**2+6*f(r)*mp.diff(u,r)**2+24*lam*u(r)**2),[lo,hi])
            check('charged_form_identity_'+label+'_'+str(prec),abs(v1/v2-1)<mp.mpf('1e-45'))
            if label=='benchmark':check('benchmark_test_form_coercive_'+str(prec),v2>e0)
            records.append(dict(digits=prec,case=label,F0=mp.nstr(F0,50),B_squared=mp.nstr(B*B,50),
                D_lambda=mp.nstr(D*lam,50),K_gain_lower=mp.nstr(24*D*lam-mp.mpf(22)/3*ep/(B*B),50),
                density_relative_errors=[mp.nstr(e,8) for e in errors],form_relative_error=mp.nstr(abs(v1/v2-1),8),
                form_ratio=mp.nstr(v2/e0,50)))
    for label in ['diagnostic','benchmark']:
        pair=[v for v in records if v['case']==label]
        for key in ['F0','B_squared','D_lambda','K_gain_lower','form_ratio']:
            check('precision_'+label+'_'+key,abs(mp.mpf(pair[0][key])/mp.mpf(pair[1][key])-1)<mp.mpf('1e-45'))
    FIXTURES.append(dict(fixture='unexpanded_density_and_actual_radial_forms',records=records,
        scope='charged local reference geometry, complete local trace only; no quantum state samples'))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes())
    old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths}
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
             frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
             protocol=dict(digits=[60,90],relative_tolerance='1e-45',continuum_machine_verified=False),
             previous_goal_turn='progress')
    error=None
    try:
        check('previous_candidate_window',old['status']=='passed' and old['positive_vacuum_scalar'] and old['full_Dirac_homogeneous_theorem_applied'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        algebra();rational_bounds();numerical()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,
        input_sha256_after=after,unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,
        finished_utc=datetime.now(timezone.utc).isoformat(),local_charged_tensor_lift=error is None,
        bulk_comparison_inverse_bound='1' if error is None else None,
        physical_boundary_completed=False,physical_coupled_inverse_proved=False,
        actual_correction_norm_proved=False,nonlinear_solution_constructed=False,
        physical_stability_proved=False,semiclassical_validity_proved=False,
        component_complete=False,goal_complete=False,formal_verification=False)
    output=ROOT/'evidence'/(STEM+'.json')
    receipt=json.loads(output.read_bytes()) if output.exists() else {'runs':[]}
    receipt['runs'].append(run)
    output.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:
            log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
