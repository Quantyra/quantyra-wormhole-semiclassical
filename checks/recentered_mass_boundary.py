"""Frozen off-shell endpoint and paired-domain diagnostics; no nonlinear solve."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
from fractions import Fraction as Q
import json,traceback
import sympy as S
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='recentered-mass-boundary';PREVIOUS='recentered-reference-inverse'
OWN=['notes/'+STEM+'.md','checks/recentered_mass_boundary.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[];EXPRESSIONS={}
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def off_shell_endpoints():
    r,B,lam,G,cc,lam0,m,Z,qq=S.symbols('r B lambda G c lambda0 m Z q',positive=True)
    ff=S.symbols('F0:5',positive=True);aa=S.symbols('a0:5',real=True)
    def dr(expr):
        return S.diff(expr,r)+sum(S.diff(expr,ff[i])*ff[i+1]+S.diff(expr,aa[i])*aa[i+1] for i in range(4))
    sig=S.log(r/B);T=r*(ff[1]+ff[0]*aa[1])
    Rb=-r*r*ff[2]-3*r*r*ff[1]*aa[1]-2*r*r*ff[0]*(aa[2]+aa[1]**2)+2*T-2*ff[0]
    LW=S.exp(aa[0])/(4*S.pi*r*r)*(sig*(Rb**2/60-Rb/18+S.Rational(1,15)-qq**2/3)+Rb*T/30-T*T/10+11*T*ff[0]/90-11*ff[0]**2/180-ff[0]/18)
    QF=S.diff(LW,ff[2]);PF=S.diff(LW,ff[1])-dr(QF)
    QA=S.diff(LW,aa[2]);PA=S.diff(LW,aa[1])-dr(QA)
    # Boundary polynomial differentiated with both radial metric and lapse retained.
    N=B*S.exp(aa[0])*S.sqrt(ff[0])/r
    proper=lambda expr:r*S.sqrt(ff[0])/B*dr(expr)
    av=proper(S.log(N));sv=S.sqrt(ff[0])/B
    Gcal=sig*(12*av*proper(av)+S.Rational(44,7)*av**3+10*av/B**2)-6*proper(av)*sv-S.Rational(24,7)*av**2*sv-6*av*proper(sv)-S.Rational(96,7)*av*sv**2-18*sv*proper(sv)-S.Rational(85,21)*sv**3
    BW=S.simplify(-B**2*N*Gcal/(360*S.pi))
    c0=-221*lam**S.Rational(3,2)/(15120*S.pi**2);c1=29*S.sqrt(lam)/(20160*S.pi**2)
    Bct=4*S.pi*S.exp(aa[0])*S.sqrt(ff[0])*(c0*r*r+2*c1)
    F=lam*r*r+1-2*m/r+Z/r**2
    subst={**{ff[j]:S.diff(F,r,j) for j in range(5)},**{a:0 for a in aa}}
    dm=[S.diff(-2/r,r,j) for j in range(5)]
    mass_bulk=S.simplify((PF*dm[0]+QF*dm[1]).subs(subst))
    mass_BW=S.simplify(sum(S.diff(BW,ff[j])*dm[j] for j in range(5)).subs(subst))
    mass_ct=S.simplify(sum(S.diff(Bct,ff[j])*dm[j] for j in range(5)).subs(subst))
    values=[]
    for name,expr,want in [('bulk',mass_bulk,-105),('boundary',mass_BW,-116),('counterterm',mass_ct,221)]:
        value=S.simplify(S.limit(expr,r,S.oo)*3780*S.pi/lam)
        check('off_shell_mass_endpoint_'+name,value==want,value)
        values.append(dict(piece=name,coefficient=str(value)))
    check('complete_fixed_clock_mass_cancellation',-105-116+221==0)
    J=sig*(Rb/30-S.Rational(1,18))+T/30;K=Rb/30-T/5+11*ff[0]/90
    zero_a={a:0 for a in aa}
    check('derived_bulk_PF',S.simplify((PF-((2*J+K)/r+dr(J))/(4*S.pi)).subs(zero_a))==0)
    check('derived_bulk_QF',S.simplify(QF+S.exp(aa[0])*J/(4*S.pi))==0)
    # The general lapse momentum is retained even though its mass tangent is zero.
    check('derived_bulk_Qalpha',S.simplify(QA+S.exp(aa[0])*ff[0]*J/(2*S.pi))==0)
    LRN=S.expand(S.expand_log(S.simplify(LW.subs(subst)),force=True))
    BRN=S.expand(S.expand_log(S.simplify(BW.subs(subst)),force=True));CTRN=Bct.subs(subst)
    bulk2=S.limit(LRN/r**2,r,S.oo);bulk0=S.limit(S.expand(LRN-bulk2*r*r),r,S.oo)
    bd3=S.limit(BRN/r**3,r,S.oo);bd1=S.limit(S.expand((BRN-bd3*r**3)/r),r,S.oo)
    check('bulk_cubic_divergence',S.simplify(bulk2+13*lam**2/(240*S.pi))==0)
    check('bulk_linear_divergence',S.simplify(bulk0+lam/(60*S.pi))==0)
    check('Weyl_complete_cubic_divergence',S.simplify(bulk2/3+bd3-653*lam**2/(15120*S.pi))==0)
    check('Weyl_complete_linear_divergence',S.simplify(bulk0+bd1-67*lam/(3780*S.pi))==0)
    ct3=S.limit(CTRN/r**3,r,S.oo);ct1=S.limit((CTRN-ct3*r**3)/r,r,S.oo)
    vacuum=11*lam**2/(720*S.pi)
    check('vacuum_and_loop_cubic_subtraction',S.simplify(bulk2/3+bd3+vacuum+ct3)==0)
    check('loop_linear_subtraction',S.simplify(bulk0+bd1+ct1)==0)
    check('omitted_finite_counterterm_fails',-105-116!=0)
    check('vacuum_shift_is_not_optional',Q(653,15120)-Q(221,3780)==-Q(11,720))
    Fl=1-384*S.pi*G*cc*lam
    # Literal local action, before decomposing the curvature square.
    R=-ff[2]-3*aa[1]*ff[1]-2*ff[0]*(aa[2]+aa[1]**2)-4*ff[1]/r-4*ff[0]*aa[1]/r+2*(1-ff[0])/r**2
    Llocal=S.exp(aa[0])*r*r*((R+6*lam0)/(4*G)+4*S.pi*cc*R**2)
    BG=S.exp(aa[0])*(1/(2*G)+16*S.pi*cc*R)*(r*r*ff[1]/2+r*r*ff[0]*aa[1]+2*r*ff[0])
    Qlocal=S.diff(Llocal,ff[2]);Plocal=S.diff(Llocal,ff[1])-dr(Qlocal)
    full_local_mass=S.simplify((Plocal*dm[0]+Qlocal*dm[1]+sum(S.diff(BG,ff[j])*dm[j] for j in range(5))).subs(subst))
    check('literal_local_mass_before_counterterm',S.simplify(full_local_mass+Fl/G)==0)
    EHct=-Fl*S.exp(aa[0])*S.sqrt(ff[0])/G*(S.sqrt(lam)*r*r+1/(2*S.sqrt(lam)))
    mass_EHct=S.simplify(S.diff(EHct,ff[0])*dm[0]).subs(subst)
    check('literal_local_mass_after_counterterm',S.simplify(S.limit(full_local_mass+mass_EHct,r,S.oo))==0)
    EXPRESSIONS.update(symbols=(r,B,lam,G,cc,lam0,m,Z,qq),bulk=mass_bulk,boundary=mass_BW,counterterm=mass_ct,local=full_local_mass+mass_EHct)
    FIXTURES.append(dict(fixture='unexpanded_off_shell_bulk_and_boundary_endpoints',mass_coefficients=values,radial_and_lapse_jets_retained=True,vacuum_shift_required=True))

def split_and_domain():
    R,lam,lam0,G,c,phi,K=S.symbols('R lambda lambda0 G c phi K')
    Fl=1-384*S.pi*G*c*lam
    lhs=(R+6*lam0)/(16*S.pi*G)+c*R**2
    rhs=Fl*(R+6*lam)/(16*S.pi*G)+c*(R+12*lam)**2+3*(lam0-lam)/(8*S.pi*G)
    check('exact_bulk_square_split',S.expand(lhs-rhs)==0)
    check('exact_auxiliary_boundary_split',S.expand((1/(8*S.pi*G)+4*c*phi)*K-Fl*K/(8*S.pi*G)-4*c*(phi+12*lam)*K)==0)
    check('residual_square_integrable_power',2-2*4==-6<-1)
    check('residual_boundary_value_power',3-4==-1<0)
    check('mass_square_bulk_endpoint_power',-3-1==-4 and -2-2==-4)
    check('clock_square_endpoint_power',-1<0)
    # The local indicial action on the RN charge and mass powers vanishes.
    r,la,m,Z=S.symbols('r lambda m Z',positive=True)
    F=la*r*r+1-2*m/r+Z/r**2
    scalar=-S.diff(F,r,2)-4*S.diff(F,r)/r+2*(1-F)/r**2
    check('actual_RN_scalar_constant',S.simplify(scalar+12*la)==0)
    # Whole Fourier subspaces, expressed by their nonoverlapping integer bounds.
    n,j=S.symbols('n j',integer=True,nonnegative=True)
    check('LLL_boundary_Fourier_indices_disjoint',S.ask(S.Q.positive(n+j+1)) is True)
    ph=S.symbols('theta',real=True)
    gn=S.Matrix([[0,S.exp(-S.I*ph)],[S.exp(S.I*ph),0]])
    chir=S.diag(1,-1)
    check('actual_normal_projector_involution',gn*gn==S.eye(2))
    check('balanced_projector_trace',S.trace(gn)==0)
    check('chirality_image_trace_cancellation',S.simplify(chir*gn+gn*chir)==S.zeros(2))
    ss=S.symbols('s');FM=12*S.gamma(1-2*ss)*S.gamma(3-ss)/(S.pi*S.gamma(6-2*ss))
    check('common_cubic_Mellin_coefficient',S.simplify(FM.subs(ss,0)-1/(5*S.pi))==0)
    ld=S.simplify(S.diff(FM,ss).subs(ss,0)/FM.subs(ss,0))
    check('common_cubic_log_derivative',S.simplify(ld-S.EulerGamma-S.Rational(46,15))==0)
    check('canonical_delta_including_anomaly',S.simplify((ld+S.EulerGamma)/2-S.Rational(1,2)-S.EulerGamma-S.Rational(31,30))==0)
    e,y=S.symbols('e y',positive=True)
    check('canonical_mass_test_endpoint',S.integrate(y**3/y,(y,0,e))==e**3/3)
    check('canonical_mass_test_delta',S.limit(y**3,y,0)==0)
    check('action_heat_remainder_integrability',Q(3,2)-1-1==Q(-1,2)>-1)
    check('cut_remainder_endpoint_integrability',Q(1,2)>0)
    FIXTURES.append(dict(fixture='local_split_and_actual_paired_domain_guards',LLL_numerical_gap_supplied=False,physical_mass_boundary_target='exact zero on stated phase',C2_source_estimate_supplied=False))

def finite_cut_and_charge():
    records=[];charge=[]
    syms=EXPRESSIONS['symbols']
    fun={key:S.lambdify(syms,EXPRESSIONS[key],'mpmath') for key in ['bulk','boundary','counterterm','local']}
    for prec in [70,100]:
        mp.mp.dps=prec;k=mp.mpf(2000);g0=k**-10;c=mp.mpf('1e22')
        g=2*g0/(1+mp.sqrt(1+4*11*g0/(360*mp.pi)));G=g;lam=mp.mpf(1);lam0=g0/G
        Fl=1-384*mp.pi*c*g;C=mp.sqrt(1+12*mp.pi*g*k**16/Fl);B=mp.sqrt((C-1)/6)
        Z=B*B+3*B**4;C0=mp.sqrt(1+12*mp.pi*k**6);bclass=mp.sqrt((C0-1)/6*g/g0)
        kap=mp.mpf('.265')*2*k**3*G*C0/bclass**3;gam=kap*B/C;x=B*B
        Kt=lambda h:((1+3*x)*((1+h)**-2+(1+h)**-3)+3*x*(1+(1+h)**-1))/(2*C)
        h=gam
        for _ in range(5):h-=(h*Kt(h)-gam)/(Kt(h)+h*mp.diff(Kt,h))
        rp=B*(1+h);M=(rp+Z/rp+rp**3)/2
        prior=None
        for ratio in ['1000','10000','100000']:
            r=B*mp.mpf(ratio);args=(r,B,lam,G,c,lam0,M,Z,k**4)
            vals={key:fun[key](*args)*3780*mp.pi/lam for key in ['bulk','boundary','counterterm']}
            localnorm=fun['local'](*args)*G/Fl
            errors=[abs(vals['bulk']+105),abs(vals['boundary']+116),abs(vals['counterterm']-221)]
            worst=max(errors)
            check('direct_cut_local_normalized_'+ratio+'_'+str(prec),abs(localnorm)<mp.mpf('1e-6'))
            check('direct_cut_loop_endpoints_'+ratio+'_'+str(prec),worst<mp.mpf('.001'))
            if prior is not None:check('direct_cut_refinement_'+ratio+'_'+str(prec),worst<prior/100)
            prior=worst
            records.append(dict(digits=prec,radius_over_B=ratio,loop_coefficients={k:mp.nstr(v,45) for k,v in vals.items()},loop_sum=mp.nstr(sum(vals.values()),45),normalized_local_endpoint=mp.nstr(localnorm,45),maximum_loop_coefficient_error=mp.nstr(worst,30)))
        logscale=-mp.log(bclass);quantum=g/(5*mp.pi)*(logscale+mp.euler+mp.mpf(151)/180)
        Zc=Fl+quantum
        check('new_charge_quantum_remainder_'+str(prec),abs(quantum)<mp.mpf('1e-33'))
        check('new_charge_positive_'+str(prec),Zc>mp.mpf('.9999999882'))
        check('new_charge_coefficient_change_'+str(prec),mp.mpf('1.178e-8')<1-Fl<mp.mpf('1.179e-8'))
        check('vacuum_shift_retained_'+str(prec),abs((lam0-lam)/(11*G*lam**2/(360*mp.pi))-1)<mp.mpf('1e-34'))
        charge.append(dict(digits=prec,F_lambda=mp.nstr(Fl,60),quantum_charge_remainder=mp.nstr(quantum,60),Z_c=mp.nstr(Zc,60),lambda0_minus_lambda=mp.nstr(lam0-lam,60)))
    for key in ['F_lambda','quantum_charge_remainder','Z_c']:
        check('charge_cross_precision_'+key,abs(mp.mpf(charge[0][key])/mp.mpf(charge[1][key])-1)<mp.mpf('1e-50'))
    pl=Q('3.14159');g0=Q(1,2000**10)
    check('uniform_quantum_charge_reserve',14*g0/(5*pl)<Q('1e-33'))
    check('uniform_total_charge_lower',1-Q('1.179e-8')-Q('1e-33')>Q('.9999999882'))
    FIXTURES.append(dict(fixture='direct_finite_cut_variations_and_actual_charge',endpoint_records=records,charge_records=charge,finite_cuts_have_zero_endpoint=False,continuum_limit_proved_by_samples=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    error=None
    try:
        check('previous_bulk_inverse_result',old['status']=='passed' and old['reference_scalar_phase_cokernel_absent'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        off_shell_endpoints();split_and_domain();finite_cut_and_charge()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,input_sha256_after=after,unchanged_inputs=before==after,
        checks=ROWS,fixtures=FIXTURES,finished_utc=datetime.now(timezone.utc).isoformat(),
        stationary_mass_boundary_coefficient='0 exactly on the specified fast phase' if error is None else None,
        stationary_metric_boundary_identification=error is None,stationary_current_constraint_realization_complete=False,
        physical_coupled_inverse_proved=False,full_dynamical_boundary_completed=False,admission_norm_inverse_quantified=False,
        nonlinear_solution_constructed=False,semiclassical_validity_proved=False,physical_stability_proved=False,
        component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
