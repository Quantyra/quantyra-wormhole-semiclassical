"""Conditional block accounting only; no complement operator implemented."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
STEM='current-weak-complement-inverse-attempt'
RUNTIME=ROOT/'evidence/current-coupled-mass-adjoint-runtime.json'
runtime=json.loads(RUNTIME.read_bytes())
sys.path.insert(0,str(Path(runtime['packages']['flint']['root']).parent))
from flint import arb,ctx

def digest(p):return sha256(p.read_bytes()).hexdigest()
def dec(v):
    d=v['decimal_enclosure']
    return arb(int(d['mid']),int(d['rad']))*arb(10)**int(d['exponent'])
def enc(v):
    m,r,e=v.mid_rad_10exp(60)
    return dict(ball=v.str(40),decimal_enclosure=dict(mid=str(m),rad=str(r),exponent=str(e)))

inputs=['notes/'+STEM+'.md','checks/current_weak_complement_inverse_attempt.py',
        'evidence/current-structured-source-graph-test.json',
        'evidence/current-quartic-ritz-rank-test-final.json',
        'checks/current_geodesic_flow_multiplier.py',
        'checks/current_global_geodesic_return.py',
        'evidence/current-coupled-mass-adjoint-runtime.json']
out=ROOT/'evidence'/(STEM+'.json')
assert not out.exists(),'Receipts are immutable.'
protected=json.loads(Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json').read_bytes())
verify=lambda:all(digest(Path(p))==h for p,h in protected.items())
assert verify()
before={p:digest(ROOT/p) for p in inputs}
data=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
          ordinary_argument_utf8=(ROOT/inputs[0]).read_text(encoding='utf-8'),
          arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'))
try:
    ctx.prec=512
    assert digest(Path(sys.executable))==runtime['python_sha256']
    for pkg in runtime['packages'].values():
        for p,h in pkg['files'].items():assert digest(Path(pkg['root'])/p)==h
    source=json.loads((ROOT/inputs[2]).read_bytes())
    assert source['status']=='source_graph_bound_checked_complete_application_unresolved'
    rank=json.loads((ROOT/inputs[3]).read_bytes())
    assert rank['status']=='passed'
    e=dec(source['constants']['actual_reference_comparison'])
    q=dec(rank['result']['constants']['q_t'])
    a,b,d=1-e,q+e,arb('.4')
    c=(1+d)/a
    h=a-b*b*d/a
    assert h>0
    delta=d*((b/h)**2+(1+b*b*c/h)**2).sqrt()
    assert 0<delta<arb('.5')
    target=arb('4.989e-32')
    n=1
    while not q*delta**n<target/2:n+=1
    trunc=q*delta**n
    ceiling=(target-trunc)/(1+q+e)
    constants=dict(e_Q=e,q_t=q,conditional_complement_defect=d,
                   conditional_complement_norm=c,conditional_Schur_inverse_norm=1/h,
                   conditional_full_block_defect=delta,
                   conditional_exact_application_residual=trunc,
                   conditional_complete_application_error_ceiling=ceiling)
    data.update(status='conditional_accounting_only_complete_application_unresolved',
                constants={k:enc(v) for k,v in constants.items()},
                conditional_iteration_count=n,actual_complement_operator=None,
                actual_complement_defect=None,actual_application_error=None,
                actual_full_residual=None,new_consequential_bound=False)
except Exception as error:
    data.update(status='failed',error=repr(error))
    raise
finally:
    data['input_sha256_after']={p:digest(ROOT/p) for p in inputs}
    data['inputs_unchanged']=before==data['input_sha256_after']
    data['protected_paths_unchanged']=verify()
    data['protected_path_count']=len(protected)
    data['finished_utc']=datetime.now(timezone.utc).isoformat()
    out.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
assert data['inputs_unchanged'] and data['protected_paths_unchanged']
print('CONDITIONAL ONLY; no complement constructed or error bound attained.')
print('Conditional iterations:',n)
for k,v in data['constants'].items():print(k+': '+v['ball'])
