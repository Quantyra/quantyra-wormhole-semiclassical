"""Necessary finite-correction test; never a whole-space inverse certificate."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import argparse
import json
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import current_quartic_ritz_rank_test as h

ROOT = Path(__file__).resolve().parents[1]
STEM = 'current-scalar-log-solver-attempt'
BASELINE = Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
INPUTS = ['notes/'+STEM+'.md', 'checks/current_scalar_log_solver_attempt.py',
    'checks/current_quartic_ritz_rank_test.py',
    'notes/current-global-scalar-log-comparison.md',
    'notes/current-physical-log-preconditioner-test.md',
    'evidence/current-physical-log-preconditioner-test.json',
    'evidence/current-quartic-ritz-rank-test-final.json',
    'evidence/current-canonical-heat-restoration-assembly.json',
    'evidence/current-finite-local-action-assembly.json',
    'evidence/current-coupled-mass-adjoint-runtime.json']


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def execute():
    h.ctx.prec=512
    a=h.arb
    assert digest(Path(sys.executable)) == h.runtime['python_sha256']
    for pkg in h.runtime['packages'].values():
        for name, expected in pkg['files'].items():
            assert digest(Path(pkg['root'])/name) == expected
    for name in ['current-physical-log-preconditioner-test', 'current-quartic-ritz-rank-test-final']:
        accepted=json.loads((ROOT/('evidence/'+name+'.json')).read_bytes())
        assert accepted['input_sha256_before']==accepted['input_sha256_after']
        for p,v in accepted['input_sha256_before'].items():
            assert digest(ROOT/p)==v
    r=accepted['result']
    local=h.receipt('current-finite-local-action-assembly')
    ep=h.pick(local,'actual_seed_correlated_endpoint')
    assert ep['exact_Q_bolt']=='0'
    end=h.dec(ep['z_h']); Q=[a(0)]*7
    for k,v in zip([0,2,3,4,5,6],ep['Q_coefficients']): Q[k]=h.dec(v)
    beta=Q[0]; ex=[a(0)]*6
    for k in range(2,7): ex[k-1]=k*Q[k]/(k-1)
    ps=[h.scale(h.add(h.shift(Q,3),h.shift(ex,4)),4), [a(0)]*4+[a(1)]]
    inner=lambda p,q:h.primitive_at(h.shift(h.mul(p,q),-4),end)
    grad=lambda p,q:h.primitive_at(h.shift(h.mul(Q,h.mul(h.diff(p),h.diff(q))),-2),end)/beta
    Lop=lambda p:h.scale(h.add(h.shift(h.mul(Q,h.diff(h.diff(p))),2),
        h.shift(h.mul(h.diff(Q),h.diff(p)),2),
        h.scale(h.shift(h.mul(Q,h.diff(p)),1),-2)),-1/beta)
    known=h.pick(h.receipt('current-canonical-heat-restoration-assembly'),'known_TT_matrix_with_shared_vacuum')
    K=[[-h.dec(known['matrix'][ii][jj])/beta for jj in [0,2]] for ii in [0,2]]
    B=[[h.dec(x) for x in row] for row in r['variational_B_lower_matrix']]
    for row in r['full_quantum_energy_upper_rows']:
        assert 2*h.dec(row['complete_geodesic_nonlocal_upper']) < a('1e-35')
    lo=(1-a('.002784'))*(1-a('1.227e-74'))/(1+a('.0022'))
    hi=(1+a('.002784')+a('.338000000000001'))*(1+a('1.227e-74'))/(1-a('.0022'))
    AL=[[lo*K[i][j]+B[i][j] for j in range(2)] for i in range(2)]
    AU=[[hi*(K[i][j]+(a('1e-35') if i==j else 0)) for j in range(2)] for i in range(2)]
    M1=[[grad(p,q) for q in ps] for p in ps]
    M3=[[grad(Lop(p),Lop(q)) for q in ps] for p in ps]
    f_lower=a('1.5')*(1-a('1.179e-8'))-60*a('1e-33')/a.pi()
    FL=[[f_lower*x for x in row] for row in M1]
    FU=[[a('1.5')*(1+a('1e-30'))*M1[i][j]+3*a('1e-33')/(20*a.pi())*(400*M1[i][j]+14*M3[i][j]) for j in range(2)] for i in range(2)]
    v=[inner(ps[1],p) for p in ps]
    det=lambda M:M[0][0]*M[1][1]-M[0][1]*M[1][0]
    slope=lambda M:M[0][0]*v[1]**2+M[1][1]*v[0]**2-2*M[0][1]*v[0]*v[1]
    M=[[a('1.49')*FU[i][j]-AL[i][j] for j in range(2)] for i in range(2)]
    N=[[AU[i][j]-a('.51')*FL[i][j] for j in range(2)] for i in range(2)]
    assert slope(M)>0 and slope(N)>0
    gmin=-det(M)/(a('1.49')*slope(M)); gmax=det(N)/(a('.51')*slope(N))
    assert a('66.7925')<gmin<gmax<a('90.2951')
    c=v[0]/v[1]
    quad=lambda T:T[0][0]-2*c*T[0][1]+c*c*T[1][1]
    assert quad(FU)>0 and quad(FL)>0
    lower_ratio=quad(AL)/quad(FU)
    upper_envelope_ratio=quad(AU)/quad(FU)
    assert upper_envelope_ratio>a('1.69') and lower_ratio<a('1.31')
    ell_upper=(1+a('.002784'))*(1+a('1.227e-74'))/(1-a('.0022'))*(quad(K)+a('1e-35')*(1+c*c))
    B_upper=a('.338000000000001')*(1+a('1.227e-74'))/(1-a('.0022'))*(quad(K)+a('1e-35')*(1+c*c))
    B_required=a('1.49')*quad(FL)-ell_upper
    assert B_required>0
    conversion=(a('1.341')/a('.997216')).sqrt()
    target=a('.49')/conversion
    constants=dict(energy_to_original_dual_conversion=conversion,
        sufficient_energy_defect_target_NOT_attained=target,
        necessary_gamma_lower=gmin, necessary_gamma_upper=gmax,
        H_orthogonal_coefficient=c, actual_Rayleigh_lower_envelope=lower_ratio,
        separated_upper_envelope_ratio=upper_envelope_ratio,
        scalar_allowance_needed_for_this_row_NOT_attained=B_required,
        current_scalar_upper_allowance=B_upper,
        required_scalar_upper_improvement_ratio=B_upper/B_required)
    return {'precision_bits':h.ctx.prec, 'constants':{k:h.enc(v) for k,v in constants.items()},
        'matrices':{k:[[h.enc(x) for x in row] for row in T] for k,T in
                    dict(A_lower=AL,A_upper=AU,F_lower=FL,F_upper=FU).items()},
        'scope':{'internal_norm_comparison_proved':True,
                 'necessary_window_nonempty':True,'specified_upper_envelope_cannot_certify':True,
                 'actual_candidate_rejected':False,'whole_complement_bound':False,
                 'implemented_inverse':False,'signed_certificate':False,
                 'physical_margin_closed':False,'construction_obstruction':False,'goal_complete':False}}


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',default='evidence/'+STEM+'.json');args=p.parse_args()
    output=ROOT/args.output
    assert not output.exists(), 'Receipts are immutable.'
    protected=json.loads(BASELINE.read_bytes())
    verify=lambda:all(digest(Path(p))==v for p,v in protected.items())
    assert verify()
    before={p:digest(ROOT/p) for p in INPUTS}
    data={'started_utc':datetime.now(timezone.utc).isoformat(),'input_sha256_before':before,
          'ordinary_argument_utf8':(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
          'arithmetic_source_utf8':Path(__file__).read_text(encoding='utf-8')}
    try:
        data['result']=execute();data['status']='internal_transfer_verified_complete_solver_unresolved'
    except Exception as exc:
        data['status']='failed';data['error']=repr(exc);raise
    finally:
        data['input_sha256_after']={p:digest(ROOT/p) for p in INPUTS}
        data['protected_paths_unchanged']=verify();data['protected_path_count']=len(protected)
        data['finished_utc']=datetime.now(timezone.utc).isoformat()
        output.write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8',newline='\n')
    assert before==data['input_sha256_after'] and data['protected_paths_unchanged']
    for k,v in data['result']['constants'].items():print(k+': '+v['ball'])
    print('Receipt: '+str(output))


if __name__=='__main__':main()
