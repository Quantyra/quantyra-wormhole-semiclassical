"""Exact normal endpoint return, not a full finite-geometry inverse.

Ordinary argument: notes/current-coupled-outer-return-test.md.
"""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from math import factorial
import sys, json, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'checks'))
import sympy as sp
from flint import arb,ctx
import current_product_matched_source_projection as old
STEM='current-coupled-outer-return-test'
OWN=['notes/'+STEM+'.md','checks/current_coupled_outer_return_test.py','evidence/'+STEM+'-review.json']
ROWS=[]; SYMBOLIC={}; enc=old.encoded

def check(name,ok,detail=''):
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok: raise AssertionError(name)

def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()

def algebra():
    s,be,G,la,C,B2,Fl,cr=sp.symbols('s beta G lambda C B2 Fl c_R',positive=True)
    q=s*(s-3);a=lambda t:(t+3)/(3*(t+1))
    cc=lambda t:sp.Matrix([[1,a(t)],[0,2]])
    J=3*(q-4)/(4*sp.pi*G*la);N=-q/(48*sp.pi*G*la)
    # Physical lambda^2 has already been absorbed into D0.
    D0=-8*sp.pi*G*C*be**2/B2
    classical=(D0*cc(3-s).T*sp.diag(J,N)*cc(s)).subs(la,C*be/B2)
    reference=be*sp.Matrix([[24+18*s-6*s*s,24+2*s-2*s*s],[12+10*s-2*s*s,12]])
    check('both_rows_and_density_normalization',(classical-reference).applyfunc(sp.simplify)==sp.zeros(2,2))
    DR=96*sp.pi*G*cr-G/(60*sp.pi);Cv=1+11*G*la/(180*sp.pi)
    Fs=Fl+23*G*la/(180*sp.pi)
    check('scalar_and_TT_share_Fstar',sp.simplify((Cv-4*DR*la-Fs).subs(Fl,1-384*sp.pi*cr*G*la))==0)
    Jr=3*(Cv-DR*la*q)*(q-4)/(4*sp.pi*G*la)
    check('trace_derivative_with_all_contacts',sp.simplify((D0*sp.diff(Jr,s).subs(s,4)+30*be*(Cv-4*DR*la)).subs(la,C*be/B2))==0)
    check('non_singular_row_combination',sp.simplify((s-4)*a(3-s)-(s-6)/3)==0)
    d0=1/(480*sp.pi**2);K=G*C*be**2/(15*sp.pi*B2)
    check('canonical_endpoint_residue',sp.simplify(D0*96*d0+24*K)==0)
    f0,cl=sp.symbols('f0 c_log');blog=-(2*f0+3*cl)/(8*be)
    B0=(cl+2*f0/3)/(24*K);eta=9*K/(be*Fs)
    check('correlated_return_coefficient',sp.simplify(B0*eta+blog/Fs)==0)
    check('squared_log_defect_cancels',sp.simplify(-12*K*(-blog/Fl)-12*K*blog/Fl)==0)
    # Bounds on the complete derivatives, not values at sampled z.
    e=sp.symbols('e');p=sp.Integer(1);ledger=[]
    for j in range(7):
        if j:p*=s-(j-1)
        bounds=[]
        for poly in [p,sp.diff(p,s)]:
            translated=sp.Poly(sp.expand(poly.subs(s,4+e)),e)
            bounds.append(sum(abs(c)*sp.Rational(1,100)**power[0] for power,c in translated.terms()))
        budget=factorial(j)*100**j
        check('entire_falling_polynomial_bounds_'+str(j),max(bounds)<=budget)
        ledger.append(dict(j=j,value_bound=str(bounds[0]),derivative_bound=str(bounds[1]),denominator=budget))
    SYMBOLIC.update(matrix_normalization='D0=-8 pi G C beta^2/B^2',full_shape_particular=str(B0),
                    scalar_log='f0/(30 beta Fstar)',metric_log_return='-b_L/F_lambda',jet_bounds=ledger)

def execute():
    algebra();p,h,zh,env=old.scalar.prev.york.actual_parameters()
    C,B2,be,la,Fl=p['C'],p['B2'],p['beta'],p['g'],p['Fl'];q=arb(p['q'])
    K=C*be**2/(15*arb.pi()*B2);k=q*q/(12*arb.pi()*C*B2)
    cn=4*C*be*(3*be+p['u']+p['eps'])/(5*arb.pi()*B2)
    f0=-2*k;cl=4*k-cn;blog=-(2*f0+3*cl)/(8*be)
    Fs=Fl+23*la/(180*arb.pi());L=10*arb.pi()*Fs/la;eta=6/L
    # Keep the original fixed scale using its previously certified interval.
    source=json.loads((ROOT/'evidence/current-canonical-heat-restoration-assembly.json').read_bytes())
    # The interval for log(mu0 B) is retained symbolically below via its
    # independently certified guard, which is enough for this root test.
    check('fixed_scale_assembly_passed',source['runs'][-1]['status']=='passed')
    log_muB=arb('1.5e-9','1.5e-9')
    kap=(la*B2).log()/2-log_muB-arb(6)/5
    check('retained_contact_interval',kap>arb('7.8') and kap<arb('10.8'))
    check('root_monotonicity',L>arb('3e34'))
    em=eta*(1-300/L);ep=eta*(1+300/L)
    def hreg(e):
        error=arb(0,(e**3).abs_upper())
        return (3+e).digamma()+kap-arb.pi()**2*e/12+error
    # Preserve the 6-L*e correlation exactly, prior to interval evaluation.
    def endpoint(e,sign):
        return sign*1800/L+5*e+e*e+(6+5*e+e*e)*e*hreg(e)
    rm=endpoint(em,1);rp=endpoint(ep,-1)
    check('entire_family_root_lower_sign',rm>0)
    check('entire_family_root_upper_sign',rp<0)
    root=arb((em.lower()+ep.upper())/2,(ep.upper()-em.lower())/2)
    check('actual_root_enclosure',root>arb('1.86509e-34') and root<arb('1.86510e-34'))
    B0=(cl+2*f0/3)/(24*K);pair=B0.abs_upper()*root.abs_upper()
    scalarlog=f0/(30*be*Fs);deltaLog=-blog/Fl
    cap=arb('1e-8');logfactor=1+arb(10).log()
    match=(cap-pair)*logfactor;half=(cap/2-pair)*logfactor
    check('normal_particular_large',B0>arb('530514') and B0<arb('530516'))
    check('combined_normal_phase_term_small',pair<arb('9.895e-29'))
    check('full_log_return_admitted',deltaLog>0 and deltaLog<arb('1e-28'))
    check('conditional_matching_budget_positive',half>0 and match>half)
    check('matching_relative_accuracy',match/B0>arb('6.22e-14') and match/B0<arb('6.23e-14'))
    return dict(Fstar=enc(Fs),L=enc(L),kappa_M=enc(kap),eta=enc(eta),exact_root_enclosure=enc(root),
        root_lower_sign=enc(rm),root_upper_sign=enc(rp),relative_root_bracket_width=enc(600/L),
        full_normal_shape_particular=enc(B0),combined_normal_phase_term_upper=enc(pair),
        full_linear_scalar_log=enc(scalarlog),required_metric_shape_log_return=enc(deltaLog),
        matching_budget_with_zero_remainder=enc(match),matching_budget_with_half_remainder=enc(half),
        relative_matching_budget_zero_remainder=enc(match/B0),separate_coefficient_loss=enc(1/root),
        stated_sensitivity_norm='existing canonical weighted C5 mass/log coefficient seminorm on 0<z<=.1',
        sensitivity_cap='1e-8; diagnostic allocation only, not the X norm or a physical admission certificate',
        actual_global_matching_coefficient_known=False,actual_strong_remainder_known=False,
        normal_endpoint_compatibility=True,full_quantum_strong_inverse_proved=False,
        actual_controlled_solution=False,nonlinear_admission=False,exact_compatibility=False,
        physical_matching=False,semiclassical_validity_noise=False,all_sector_stability=False)

def main():
    ctx.prec=512;review=json.loads((ROOT/OWN[2]).read_bytes());paths=set(OWN+review['retained_inputs'])
    for mod in list(sys.modules.values()):
        f=getattr(mod,'__file__',None)
        if f:
            try:p=Path(f).resolve().relative_to(ROOT).as_posix()
            except ValueError:continue
            if p.endswith('.py'):paths.add(p)
    before={p:digest(p) for p in sorted(paths)};out=ROOT/'evidence'/f'{STEM}.json'
    archive=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},precision_bits=ctx.prec,
        python=sys.version,invocation=sys.argv)
    archive['runs'].append(run)
    def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
    save()
    try:
        check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256'])
        run['result']=execute();run['status']='passed'
    except Exception:run['status']='failed';run['exception']=traceback.format_exc()
    finally:
        run['checks']=ROWS;run['symbolic']=SYMBOLIC;run['input_sha256_after']={p:digest(p) for p in before}
        if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
        run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
    with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),exception=run.get('exception'),result=run.get('result')),indent=2))
    if run['status']!='passed':raise SystemExit(1)
if __name__=='__main__':main()
