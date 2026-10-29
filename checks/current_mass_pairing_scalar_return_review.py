"""Author review of ordinary equations, full budget and frozen arithmetic."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent))
import current_filtered_test_mass_trace as m
arb,ctx,sp=m.arb,m.ctx,m.sp
ROOT,BASELINE,runtime=m.ROOT,m.BASELINE,m.runtime
STEM='current-mass-pairing-scalar-return'
digest=lambda p:sha256(p.read_bytes()).hexdigest()


def execute():
    ctx.prec=768
    source=ROOT/'evidence'/f'{STEM}-attempt2.json'
    rec=json.loads(source.read_bytes())
    assert rec['status']=='complete_mass_contribution_bounded_full_diagnostic_open'
    assert rec['input_sha256_before']==rec['input_sha256_after']
    assert all(digest(ROOT/p)==h for p,h in rec['input_sha256_after'].items())
    assert digest(Path(sys.executable))==runtime['python_sha256']
    for pkg in runtime['packages'].values():
        assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
    protected=json.loads(BASELINE.read_bytes())
    assert all(digest(Path(p))==h for p,h in protected.items())
    c={k:m.dec(x) for k,x in rec['result']['constants'].items()}
    # Independently differentiate the original inverse York equations.
    z=sp.Symbol('z',positive=True)
    Q=sp.Function('Q')(z);y=sp.Function('y')(z)
    u=sp.Function('u')(z);S=sp.Function('S')(z)
    Py=6/z**2+sp.Rational(2,3)*(sp.diff(Q,z,2)-4*sp.diff(Q,z)/z)/Q
    scalar=u+z**2*sp.diff(y,z)/4-z*y/2
    raw=sp.expand(Q*sp.diff(scalar,z,2)+(sp.diff(Q,z)-2*Q/z)*sp.diff(scalar,z))
    raw=raw.subs(sp.diff(y,z,3),sp.diff(Py,z)*y+Py*sp.diff(y,z)-sp.diff(S,z))
    raw=raw.subs(sp.diff(y,z,2),Py*y-S)
    Cy=-2*Q/z+sp.Rational(5,3)*sp.diff(Q,z)-sp.Rational(2,3)*z*sp.diff(Q,z,2)+z**2*sp.diff(Q,z,3)/6
    By=Q+z**2*sp.diff(Q,z,2)/6-sp.Rational(2,3)*z*sp.diff(Q,z)
    local=Q*sp.diff(u,z,2)+(sp.diff(Q,z)-2*Q/z)*sp.diff(u,z)-z**2*sp.diff(Q*S,z)/4
    assert sp.simplify(raw-local-Cy*y-By*sp.diff(y,z))==0
    ps=sp.Function('psi')(z);pt=sp.Function('psiT')(z)
    Sz=sp.Rational(2,3)*(z*sp.diff(pt,z)-pt)
    localfact=z**2*(3*sp.diff(Q,z)*ps+(4*Q+z*sp.diff(Q,z))*sp.diff(ps,z)
        +z*Q*sp.diff(ps,z,2)-(sp.diff(Q,z)*Sz+sp.Rational(2,3)*z*Q*sp.diff(pt,z,2))/4)
    assert sp.simplify(local.subs({u:z**3*ps,S:Sz}).doit()-localfact)==0
    # Euler scaled derivatives are derived again, independently of recurrence.
    value=z**3*ps
    for k,row in enumerate(rec['result']['exact_euler_coefficients']):
        expected=z**3*sum(row[j]*z**j*sp.diff(ps,z,j) for j in range(k+1))
        assert sp.simplify(value-expected)==0
        value=z*sp.diff(value,z)
    # The weighted profile's physical measure and finite reconstruction
    # factors are present in BOTH weighted and unweighted trial bounds.
    assert c['complete_weighted_TT_reconstruction']>1
    assert c['complete_unweighted_TT_reconstruction']>1
    assert c['complete_scalar_remainder_inverse_ratio']<arb('1e-26')
    assert c['scalar_P_gap']>arb('.1') and c['scalar_D_gap']>arb('.99')
    assert c['complete_weighted_Fref_ratio']<1
    assert c['scalar_bolt_average_upper']<arb('1e-6')
    for values in rec['result']['global_polynomial_Bernstein_guards'].values():
        assert all(m.dec(x)>0 for x in values)
    names=['scalar_observation_mass_debit','scalar_mass_mixed_main_debit',
           'direct_scalar_mass_mixed_remainder_debit','scalar_return_mixed_main_debit',
           'direct_scalar_return_mixed_remainder_debit','complete_local_canonical_pairing_debit']
    assert sum((c[n] for n in names),arb(0)).overlaps(c['complete_local_scalar_mass_block_upper'])
    emass=c['unchanged_massive_and_LLL_mass_block_upper']+c['complete_local_scalar_mass_block_upper']
    assert emass.overlaps(c['complete_mass_vector_norm_upper'])
    target=arb('1.49e34');conv=c['full_E_to_retained_Et_comparison'];eb=c['unchanged_full_response_bridge']
    remain=(target-eb)/conv-emass
    assert remain.overlaps(c['conditional_full_E_account_for_unknown_mass_free_and_errors'])
    assert c['generic_energy_estimator_over_whole_Dj_target']>1
    assert c['direct_mixed_remainder_gain']>1
    assert rec['result']['numerical_contract_pass']==bool(remain>0)
    assert all(rec['result'][k] is None for k in ['actual_mass_free_residual','remaining_application_errors','actual_Dj','actual_Gt'])
    scope=rec['result']['scope']
    assert all(not scope[k] for k in ['complete_signed_response','complete_Dj_diagnostic',
        'full_bolt_allowance','available_physical_margin','original_P32_defect',
        'actual_kappa_m','all_original_coupled_outputs','nonlinear_admission_contraction',
        'nonlinear_neighborhood_transfer','nonlinear_existence','exact_compatibility',
        'physical_matching','semiclassical_validity_noise','all_sector_stability',
        'construction_obstruction','full_goal_complete'])
    first=ROOT/'evidence'/f'{STEM}-attempt1.json'
    return dict(verification='author_review_supporting_ordinary_proof',
        accepted_receipt=str(source.relative_to(ROOT)),accepted_receipt_sha256=digest(source),
        frozen_inputs_checked=len(rec['input_sha256_after']),whole_support_cells=rec['result']['whole_cell_count'],
        protected_paths_unchanged=len(protected),pinned_python_sha256=runtime['python_sha256'],
        rejected_attempts=[dict(path=str(first.relative_to(ROOT)),sha256=digest(first),
            reason='Weighted and unweighted profile bounds omitted the actual TT reconstruction. Numerical pass is not scientifically accepted.')],
        ordinary_hypotheses_reviewed=['same fixed-seed scalar graph closure and original York endpoints',
            'bounded physical weight and uniform truncated auxiliary weights',
            'form factorization retains complete anomaly and off-shell remainder',
            'actual canonical source rows and all cutoff derivatives',
            'paired local endpoint products vanish; accepted quantum block reused without modification',
            'full auxiliary resolvent kernel, complete rational spectrum and Volterra return',
            'direct Laplacian/Bochner norms substitute into complete accepted mixed tensor bound',
            'both residual equations with physical alpha and all mass debits',
            'diagnostic account only; all mass-free and remaining errors unknown'],
        independent_specialist_review=False,Lean_verification=False,
        conclusion='Complete mass contribution bounded; generic energy estimator excessive; full diagnostic and signed response unresolved.')


if __name__=='__main__':
    out=ROOT/'evidence'/f'{STEM}-verification.json'
    assert not out.exists(),'Preserve immutable receipts'
    result=execute()
    result.update(finished_utc=datetime.now(timezone.utc).isoformat(),review_source_sha256=digest(Path(__file__)))
    out.write_bytes((json.dumps(result,indent=2)+'\n').encode('utf-8'))
    print(json.dumps(result,indent=2))
