"""All-spectrum constant-rescaling consistency test, not a strong response bound."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import sys, json, traceback
if '--arb-path' in sys.argv:
    sys.path.insert(0, sys.argv[sys.argv.index('--arb-path')+1])
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'checks'))
import sympy as S
from flint import arb, ctx
import current_product_matched_source_projection as old
STEM='current-common-regulator-conversion-test'
OWN=['notes/'+STEM+'.md','checks/current_common_regulator_conversion_test.py',
     'evidence/'+STEM+'-review.json']
ROWS=[]

def check(name,condition):
    assert name not in {v['name'] for v in ROWS}
    ROWS.append(dict(name=name,passed=bool(condition)))
    if not condition: raise AssertionError(name)

def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()

def execute():
    b,s,t,e,A0,A1,A2,C,H0,H1=S.symbols('b s t epsilon A0 A1 A2 C H0 H1', positive=True)
    amin2=A0*2*C*S.exp(-2*b)/(4*S.pi)
    amin1=(A0*H0+A1*2*C*S.exp(-2*b))/(4*S.pi)
    azero=(A0*H1*S.exp(2*b)+A1*H0+2*C*A2*S.exp(-2*b))/(4*S.pi)
    prod=S.expand((A0/t+A1+A2*t)*(2*C*S.exp(-2*b)/t+H0+H1*S.exp(2*b)*t)/(4*S.pi))
    for i,val in [(-2,amin2),(-1,amin1),(0,azero)]:
        check('common_heat_coefficient_'+str(i),S.simplify(prod.coeff(t,i)-val)==0)
    integral=S.integrate((amin2/t**2+amin1/t+azero)/(2*t),(t,e*S.exp(-2*b),e))
    expected=amin2*(S.exp(4*b)-1)/(4*e**2)+amin1*(S.exp(2*b)-1)/(2*e)+b*azero
    check('exact_cutoff_difference',S.simplify(integral-expected)==0)
    first=S.diff(b*azero,b).subs(b,0)
    second=S.diff(b*azero,b,2).subs(b,0)
    check('finite_part_first_derivative',S.simplify(first-(A0*H1+A1*H0+2*C*A2)/(4*S.pi))==0)
    check('finite_part_second_derivative',S.simplify(second-(A0*H1-2*C*A2)/S.pi)==0)
    z0,z1=S.symbols('zeta0 zeta1')
    check('common_zeta_conversion',S.diff((S.exp(2*b*s)-1)*(z0+s*z1)/2,s).subs(s,0)==b*z0)
    pars,h,zh,env=old.scalar.prev.york.actual_parameters()
    be,u,ep,ta=pars['beta'],pars['u'],pars['eps'],pars['tau']
    Cn,B2,q=pars['C'],pars['B2'],arb(pars['q'])
    coeff=[be,arb(0),u,-4*be-2*u-ta,3*be+u+ep,-2*ep,ep]
    evalq=lambda z:sum((v*z**i for i,v in enumerate(coeff)),arb(0))
    evalqp=lambda z:sum((i*coeff[i]*z**(i-1) for i in range(1,len(coeff))),arb(0))
    check('actual_horizon_kept',zh.overlaps(1/(1+h)))
    check('actual_bolt_root',evalq(zh).contains(0))
    check('actual_bolt_slope',(evalqp(zh)+2*pars['gamma']).contains(0))
    r=[(i+2)*(i+1)*coeff[i+2] for i in range(len(coeff)-2)]
    integral_r2=sum((r[i]*r[j]*zh**(i+j+1)/(i+j+1)
                     for i in range(len(r)) for j in range(len(r))),arb(0))
    A0n=2*zh;A1n=evalqp(zh)/6;A2n=-integral_r2/120
    H0n=-q-arb(1)/3;H1n=(q*q/6-arb(1)/30)/Cn
    vol=-2*A0n*H1n/(arb.pi()*B2)
    curv=4*Cn*A2n/(arb.pi()*B2)
    j=vol+curv
    jfirst=-(A0n*H1n+A1n*H0n+2*Cn*A2n)/(2*arb.pi()*B2)
    check('positive_tower_mass_coefficient',H1n>0)
    check('nonzero_integrated_curvature',integral_r2>0)
    check('negative_volume_contact',vol<0)
    check('negative_curvature_contact',curv<0)
    check('nonzero_conversion_Hessian',j<0)
    check('certified_interval_width',j.rad()<arb('1e-48'))
    enc=old.encoded
    return dict(symbolic_cutoff_identity=str(expected),symbolic_common_finite_part=str(b*azero),
      actual_parameters={k:enc(pars[k]) for k in ['C','B2','beta','u','eps','tau','gamma']},
      q=str(pars['q']),z_h=enc(zh),integrated_R_squared=enc(integral_r2),
      A0=enc(A0n),A1=enc(A1n),A2=enc(A2n),H0=enc(H0n),H1=enc(H1n),
      normalized_conversion_first_derivative=enc(jfirst),
      normalized_conversion_second_derivative=enc(j),
      volume_second_derivative=enc(vol),curvature_second_derivative=enc(curv),
      curvature_to_volume_ratio=enc(curv/vol),
      all_spectrum_constant_family_common_limit_proved=True,
      two_dimensional_only_conversion_rejected=True,
      physical_correction_direction=False,variable_conversion_proved=False,
      full_signed_projection_executed=False,strong_response_bound_executed=False,
      selected_remaining_feedback_target='.18368',proposed_heat_target='.18',
      proposed_other_global_target='.00368',
      actual_controlled_solution=False,uniform_nonlinear_admission=False,
      exact_compatibility=False,physical_matching=False,semiclassical_validity_noise=False,
      all_sector_stability=False)

def main():
    ctx.prec=512
    review=json.loads((ROOT/OWN[2]).read_bytes())
    paths=set(OWN+review['retained_inputs'])
    for mod in list(sys.modules.values()):
        f=getattr(mod,'__file__',None)
        if f:
            try:p=Path(f).resolve().relative_to(ROOT).as_posix()
            except ValueError:continue
            if p.endswith('.py'):paths.add(p)
    before={p:digest(p) for p in sorted(paths)}
    out=ROOT/'evidence'/f'{STEM}.json'
    archive=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),status='running',
      input_sha256_before=before,frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},
      precision_bits=ctx.prec,python=sys.version,invocation=sys.argv)
    archive['runs'].append(run)
    def save():out.write_bytes((json.dumps(archive,indent=2)+'\n').encode())
    save()
    try:
        check('ordinary_argument_frozen',digest(OWN[0])==review['ordinary_argument_sha256'])
        run['result']=execute();run['status']='passed'
    except Exception:run['status']='failed';run['exception']=traceback.format_exc()
    finally:
        run['checks']=ROWS;run['input_sha256_after']={p:digest(p) for p in before}
        if before!=run['input_sha256_after']:run['status']='failed';run['input_changed']=True
        run['finished_utc']=datetime.now(timezone.utc).isoformat();save()
    with (ROOT/'evidence'/f'{STEM}.log').open('ab') as f:f.write((json.dumps(run)+'\n').encode())
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),exception=run.get('exception'),result=run.get('result')),indent=2))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
