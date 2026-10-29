"""Uniform full-spectrum positive-time bound; short-time action unresolved."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
import json
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent))
from current_quartic_ritz_rank_test import arb,ctx,dec,enc

ROOT=Path(__file__).resolve().parents[1]
STEM='current-positive-time-quantum-operator'
BASELINE=Path('C:/Users/dfred/AppData/Local/Temp/s169-weak-response-execution-baseline.json')
INPUTS=['notes/'+STEM+'.md','checks/current_positive_time_quantum_operator.py',
        'checks/current_quartic_ritz_rank_test.py',
        'notes/current-structured-weak-response-attempt.md',
        'notes/current-quartic-column-energy-dual.md',
        'notes/current-canonical-heat-restoration-assembly.md',
        'notes/current-coupled-truncation-energy-error.md',
        'notes/current-disk-heat-propagation.md',
        'evidence/current-quartic-column-energy-dual.json',
        'evidence/current-quartic-column-energy-dual-verification.json',
        'evidence/current-coupled-truncation-energy-error.json',
        'evidence/current-complete-polynomial-weak-transfer.json',
        'evidence/current-first-quantum-solve-attempt.json',
        'evidence/current-first-inverse-source-screen.json',
        'evidence/current-scalar-green-source-residual.json',
        'evidence/current-coupled-mass-adjoint-runtime.json']


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def result(name):
    d=json.loads((ROOT/'evidence'/(name+'.json')).read_bytes())
    r=d.get('runs',[d])[-1]
    assert r['status']!='failed'
    return r['result']


def center_and_error(v):
    m,_,exponent=v.mid_rad_10exp(65)
    center=arb(int(m))*arb(10)**int(exponent)
    error=arb(abs(v-center).upper())
    assert (center+arb(0,error)).contains(v)
    return center,error,dict(mantissa=str(m),power_of_ten=int(exponent))


def execute():
    ctx.prec=512
    energy=result('current-coupled-truncation-energy-error')
    C,B2,beta=[dec(energy['parameters'][k]) for k in ['C','B2','beta']]
    epsilon=dec(result('current-complete-polynomial-weak-transfer')['constants']['epsilon'])
    CLt=dec(energy['global_shape_L2_constant'])/(1-epsilon).sqrt()
    CDt2=(100*arb.pi()*B2*beta/C+2*arb(500)**4*dec(energy['global_shape_L2_constant'])**2)/(1-epsilon)
    CDt=CDt2.sqrt()
    col=result('current-quartic-column-energy-dual')['constants']
    M=dec(col['heat_diagonal_bound_M'])
    assert M<734 and CLt.overlaps(dec(col['CLt']))
    t0=arb('.001');q=arb(2000)**4;a=q/C
    assert (q+1)/C>324 and a>324
    c=1/(B2*beta)
    kplus=256*C*M*CLt**2/(B2*beta*324*t0**2)
    N=256;n0=N+1
    I0=arb(0);I2=arb(0)
    for n in range(1,N+1):
        mn=n*(n+q)/C;dn=2*(q+2*n);x=mn*t0
        e2=x.expint(2)
        assert e2>0
        I0+=dn*e2/t0
        I2+=dn*(-x).exp()/mn
    zeta=(-a*t0).exp()
    Stail=zeta**n0*((2*q+4*n0)/(1-zeta)+4*zeta/(1-zeta)**2)
    tail2=Stail/(a*n0);tail0=tail2/t0**2
    I0full=I0+arb(0).union(arb(tail0.upper()))
    I2full=I2+arb(0).union(arb(tail2.upper()))
    a0=-2*c*I0full/arb.pi();a2=c*I2full/(480*arb.pi())
    a0center,da0,a0literal=center_and_error(a0)
    a2center,da2,a2literal=center_and_error(a2)
    H4norm=16*CLt**2+32*CLt*CDt+8*CDt2
    kcontact=da0*CLt**2+da2*H4norm
    local_size=abs(a0)*CLt**2+abs(a2)*H4norm
    tail_contact_error=(2*c*tail0/arb.pi())*CLt**2+(c*tail2/(480*arb.pi()))*H4norm
    k=kplus+kcontact
    assert k<1
    first=result('current-first-quantum-solve-attempt')['constants']
    Ghat=dec(first['proved_source_dual_upper'])
    Kfirst=dec(result('current-first-inverse-source-screen')['constants']['first_residual_to_original_bolt_multiplier'])
    Es=Ghat*k/(1-k);physical=Kfirst*Es
    source=result('current-scalar-green-source-residual')['constants']
    subtotal=dec(source['known_partial_subtotal_not_free_margin'])+physical
    screen=arb('1e33')
    conditional_short_residual=(1-k)*screen-k*Ghat
    constants=dict(CL_t=CLt,CD_t=CDt,actual_heat_diagonal_M=M,
                   positive_time_operator_norm_upper=kplus,
                   full_I0=I0full,full_I2=I2full,
                   angular_I0_tail=tail0,angular_I2_tail=tail2,
                   a0_coefficient=a0,a2_coefficient=a2,
                   a0_center=a0center,a2_center=a2center,
                   a0_absolute_coefficient_error=da0,a2_absolute_coefficient_error=da2,
                   H4_energy_bilinear_norm_upper=H4norm,
                   complete_local_reallocation_size_upper=local_size,
                   contact_coefficient_error_operator_norm=kcontact,
                   omitted_angular_contact_operator_error_upper=tail_contact_error,
                   total_operator_difference_upper=k,
                   first_inverse_energy_error_upper=Es,
                   first_inverse_original_bolt_debit=physical,
                   first_inverse_fraction_of_whole_bolt_allowance=physical/arb('1.667e-9'),
                   known_partial_subtotal_NOT_free_margin=subtotal,
                   conditional_short_residual_ceiling=conditional_short_residual,
                   fraction_of_first_residual_screen_spent=Es/screen)
    assert all(v.is_finite() for v in constants.values())
    return dict(precision_bits=ctx.prec,positive_time_threshold='0.001',
                contact_angular_terms=256,complete_angular_tail_retained=True,
                constants={k:enc(v) for k,v in constants.items()},
                implemented_contact_coefficients={'a0':a0literal,'a2':a2literal},
                actual_short_time_column=None,actual_first_quantum_trial=None,
                actual_complete_signed_residual=None,
                scope=dict(full_energy_space_positive_time_bound=True,
                           complete_angular_and_thermal_sums=True,
                           common_regulator_contact_reallocation_retained=True,
                           contact_coefficient_application_errors_certified=True,
                           exact_first_inverse_perturbation_bound=True,
                           full_first_inverse_numerically_applied=False,
                           complete_signed_Zr_Gt=False,full_bolt_or_coupled_admission=False,
                           neighborhood_or_independent_phase_comparison=False,
                           nonlinear_existence_compatibility_matching_validity_stability=False,
                           construction_obstruction=False,full_goal_complete=False))


if __name__=='__main__':
    out=ROOT/'evidence'/(STEM+'.json')
    if out.exists():
        raise RuntimeError('Preserve immutable scientific evidence.')
    before={p:digest(ROOT/p) for p in INPUTS}
    protected=json.loads(BASELINE.read_bytes())
    assert all(digest(Path(p))==h for p,h in protected.items())
    rt=json.loads((ROOT/INPUTS[-1]).read_bytes());pkg=rt['packages']['flint']
    assert digest(Path(sys.executable))==rt['python_sha256']
    assert all(digest(Path(pkg['root'])/p)==h for p,h in pkg['files'].items())
    rec=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
             ordinary_argument_utf8=(ROOT/INPUTS[0]).read_text(encoding='utf-8'),
             arithmetic_source_utf8=Path(__file__).read_text(encoding='utf-8'),
             python_sha256=digest(Path(sys.executable)))
    try:
        rec['result']=execute()
        after={p:digest(ROOT/p) for p in INPUTS}
        assert after==before
        assert all(digest(Path(p))==h for p,h in protected.items())
        rec.update(status='positive_time_operator_enclosed_short_time_signed_solve_unresolved',
                   input_sha256_after=after,protected_paths_unchanged=len(protected),
                   pinned_flint_files_unchanged=len(pkg['files']))
    except Exception:
        import traceback
        rec.update(status='failed',traceback=traceback.format_exc())
    rec['finished_utc']=datetime.now(timezone.utc).isoformat()
    out.write_bytes((json.dumps(rec,indent=2)+'\n').encode('utf-8'))
    cs=rec.get('result',{}).get('constants',{})
    keys=['positive_time_operator_norm_upper','complete_local_reallocation_size_upper',
          'contact_coefficient_error_operator_norm','total_operator_difference_upper',
          'first_inverse_original_bolt_debit','known_partial_subtotal_NOT_free_margin',
          'conditional_short_residual_ceiling']
    print(json.dumps(dict(status=rec['status'],error=rec.get('traceback'),
                         constants={k:cs[k]['ball'] for k in keys if k in cs}),indent=2))
    if rec['status']=='failed':
        raise SystemExit(1)
