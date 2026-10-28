"""Support the local-reference rejection; do not certify a numerical F_n."""
from pathlib import Path
from hashlib import sha256
from fractions import Fraction as F
from datetime import datetime, timezone
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-local-reference-range-test'
runtime = json.loads((ROOT/'evidence/current-coupled-mass-adjoint-runtime.json').read_bytes())
sys.path.insert(0, str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb, ctx

INPUTS = [
    'notes/'+STEM+'.md', 'checks/current_local_reference_range_test.py',
    'notes/current-finite-local-action-assembly.md',
    'notes/current-schur-residual-phase-domain.md',
    'notes/recentered-quantum-tt-form.md',
    'notes/current-spatial-inverse-output-test.md',
    'evidence/current-spatial-inverse-output-test.json',
    'evidence/current-coupled-truncation-energy-error.json',
    'evidence/current-coupled-mass-adjoint-runtime.json',
]

def digest(p):
    return sha256(p.read_bytes()).hexdigest()

def dec(v):
    d=v['decimal_enclosure']
    return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])

def enc(v):
    m,r,e=v.mid_rad_10exp(60)
    return {'ball':v.str(40),'decimal_enclosure':
            {'mid':str(m),'rad':str(r),'exponent':str(e)}}

def execute():
    ctx.prec=512
    assert digest(Path(sys.executable))==runtime['python_sha256']
    for pkg in runtime['packages'].values():
        for name,expected in pkg['files'].items():
            assert digest(Path(pkg['root'])/name)==expected
    # Exact highest-jet algebra from the original local action, not a fit.
    L1,M1=F(-2),F(-1,3)
    phi_coefficient=F(2,60)*L1**2
    constant=F(2,30)*L1*M1-F(2,10)*M1**2
    assert phi_coefficient==F(2,15) and constant==F(1,45)
    assert -2-6*M1==0
    q=arb(2000)**4
    kappa=(q+1).digamma()/2-arb.const_euler()/2+arb(0,'1e-8')
    left=kappa+arb('1e-10').log()-arb(1)/6
    right=kappa+arb('2e-10').log()-arb(1)/6
    assert right<arb('-7.58')
    prior=json.loads((ROOT/'evidence/current-spatial-inverse-output-test.json').read_bytes())
    c=prior['result']['constants']
    eps=dec(c['epsilon']);qt=arb('.338')/(1-eps);ft=qt/(1+qt)
    e0=dec(c['geodesic_Et_comparison_error'])
    assert 0<e0<1
    residual=e0*ft/(1-e0)
    assert residual<arb('5.569803e-4')
    W0=dec(c['W0_complete_scalar_physical_upper'])
    W1=dec(c['W1_complete_mixed_physical_upper'])
    zeta=e0/(1-e0)
    debit=W0*(3*zeta+3*zeta**2+zeta**3)+W1*(2*zeta+zeta**2)
    ratio=debit/arb('1.667e-9')
    assert arb('1.0164e23')<ratio<arb('1.0165e23')
    p=json.loads((ROOT/'evidence/current-coupled-truncation-energy-error.json').read_bytes())['runs'][-1]
    assert p['status']=='passed'
    beta,C,B2=[dec(p['result']['parameters'][k]) for k in ['beta','C','B2']]
    pref=C/(15*arb.pi()*B2*beta)
    assert pref*right<arb('-5.65e-33')
    return {
        'constants':{k:enc(v) for k,v in {
            'kappa4_complete_angular':kappa,
            'a4_at_left':left,'a4_at_right':right,
            'p4_at_right_original_energy_units':pref*right,
            'reference_only_residual_upper_exact_applications':residual,
            'numerical_application_error_multiplier':1+qt,
            'reference_only_output_error_upper_expression':debit,
            'reference_only_estimator_over_entire_bolt_allowance':ratio,
        }.items()},
        'numerical_Fn_application_error':None,
        'complete_numerical_residual_bound':None,
        'full_signed_G_evaluated':False,
        'local_positive_reference_rejected':True,
        'physical_requirements_closed':[],
        'construction_obstruction_proved':False,
        'scope':'ordinary actual-domain reference test; requested numerical residual unresolved',
    }

if __name__=='__main__':
    out=ROOT/'evidence'/(STEM+'.json')
    if out.exists(): raise SystemExit('Immutable receipt exists; refusing overwrite.')
    before={p:digest(ROOT/p) for p in INPUTS}
    receipt={'started_utc':datetime.now(timezone.utc).isoformat(),
             'input_sha256_before':before,
             'ordinary_argument_utf8':(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
             'arithmetic_source_utf8':Path(__file__).read_text(encoding='utf-8')}
    try:
        receipt['result']=execute()
        receipt['status']='local_positive_reference_rejected_complete_residual_unresolved'
    except BaseException as exc:
        receipt['status']='failed';receipt['failure']=repr(exc)
        raise
    finally:
        receipt['input_sha256_after']={p:digest(ROOT/p) for p in INPUTS}
        receipt['finished_utc']=datetime.now(timezone.utc).isoformat()
        with out.open('xb') as f:
            f.write((json.dumps(receipt,indent=2)+'\n').encode('utf-8'))
    assert before==receipt['input_sha256_after']
    print(json.dumps({'status':receipt['status'],'constants':{
        k:v['ball'] for k,v in receipt['result']['constants'].items()}},indent=2))
