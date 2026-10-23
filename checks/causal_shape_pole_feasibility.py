"""Causal poles of the retained local model; no finite-geometry pole claim."""
from pathlib import Path
from hashlib import sha256
from datetime import datetime, timezone
from fractions import Fraction as Q
import json,traceback
import sympy as sp
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
STEM='causal-shape-pole-feasibility'; PREVIOUS='current-seed-full-trace-defect'
OWN=['notes/'+STEM+'.md','checks/causal_shape_pole_feasibility.py','evidence/'+STEM+'-review.json']
ROWS=[];FIXTURES=[]
def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check(name,ok,detail=''):
    assert name not in {r['name'] for r in ROWS}
    ROWS.append(dict(name=name,passed=bool(ok),detail=str(detail)))
    if not ok:raise AssertionError(name)

def algebra_and_rational():
    om=sp.symbols('omega')
    a,z,y=sp.symbols('a zeta y',positive=True)
    D=1-om**2/a*sp.log(-om**2/z)
    check('causal_model_first_derivative',sp.simplify(sp.diff(D,om)+2*om/a*(sp.log(-om**2/z)+1))==0)
    check('causal_model_second_derivative',sp.simplify(sp.diff(D,om,2)+2/a*(sp.log(-om**2/z)+3))==0)
    check('imaginary_frequency_returns_Euclidean_model',sp.simplify(D.subs(om,sp.I*y)-(1+y**2/a*sp.log(y*y/z)))==0)
    th=sp.symbols('theta',positive=True)
    ft=-th*sp.cot(th)+sp.log(th/sp.sin(th))
    fp=th/sp.sin(th)**2+1/th-2*sp.cot(th)
    check('monotone_branch_equation_derivative',sp.trigsimp(sp.diff(ft,th)-fp)==0)
    U=-th*sp.cot(th)
    check('Lambert_negative_axis_imaginary_part',sp.trigsimp(U*sp.sin(th)+th*sp.cos(th))==0)
    check('Lambert_negative_axis_real_part',sp.trigsimp(U*sp.cos(th)-th*sp.sin(th)+th/sp.sin(th))==0)
    r=Q('.005');d1=2/Q('.83')*(1-Q(1,90));ldiff=2*r/(Q('.80')-r)
    d2=Q(2,62)*(105+ldiff+3)
    check('Rouche_disk_stays_in_upper_half_plane',Q('.012')-r>0)
    check('Rouche_simple_root_derivative_margin',d1>Q('2.38'),float(d1))
    check('Rouche_second_derivative_margin',d2<Q('3.5'),float(d2))
    margin=Q('2.38')*r-Q('3.5')*r*r/2
    check('Rouche_model_boundary_margin',margin==Q('.01185625') and margin>Q('.0118'),str(margin))
    check('unknown_physical_remainder_target',Q('.01')<margin)
    check('unchanged_nonlinear_targets',Q('5e-9')+Q('.25')*Q('1e-8')+Q('5e7')*Q('1e-8')**2/2==Q('1e-8'))
    FIXTURES.append(dict(fixture='causal_model_and_conditional_transfer_identities',disk_radius=str(r),
        model_modulus_lower_bound=str(margin),required_physical_remainder='<.01',actual_remainder_bounded=False))

def constants(ctx):
    k=ctx.mpf(2000);g0=k**-10;cstar=ctx.mpf('1e22')-1/(5760*ctx.pi**2)
    g=2*g0/(1+ctx.sqrt(1+4*ctx.mpf(11)/(360*ctx.pi)*g0))
    F=1-384*ctx.pi*cstar*g;x0=(ctx.sqrt(1+12*ctx.pi*k**6)-1)/6
    a=20*ctx.pi*F;zeta=g0/x0*ctx.exp(ctx.mpf(12)/5)
    return a,zeta,1/(96*ctx.pi*cstar),g

def numerical():
    records=[]
    for digits in [80,120]:
        mp.mp.dps=digits;a,ze,m2,g=constants(mp);loga=mp.log(a/ze)
        f=lambda th:-th/mp.tan(th)+mp.log(th/mp.sin(th))-loga
        lo=mp.pi/2;hi=mp.pi-mp.mpf('1e-5')
        check('monotone_bracket_'+str(digits),f(lo)<0<f(hi))
        for _ in range(4*digits):
            mid=(lo+hi)/2
            if mid==lo or mid==hi:break
            if f(mid)<0:lo=mid
            else:hi=mid
        theta=(lo+hi)/2;U=-theta/mp.tan(theta);W=mp.mpc(U,-theta)
        M=abs(W);om=mp.sqrt(a/W)
        direct=mp.sqrt(a/mp.lambertw(-a/ze,-1))
        residual=1-om**2/a*mp.log(-om**2/ze)
        check('independent_Lambert_branch_'+str(digits),abs(om-direct)<mp.mpf('1e-65'))
        check('causal_sheet_identity_'+str(digits),abs(mp.log(-om**2/ze)-W)<mp.mpf('1e-65'))
        check('growing_pole_residual_'+str(digits),abs(residual)<mp.mpf('1e-65'))
        omminus=-mp.conj(om)
        check('opposite_frequency_growing_pole_'+str(digits),mp.im(omminus)>0 and abs(1-omminus**2/a*mp.log(-omminus**2/ze))<mp.mpf('1e-65'))
        check('positive_real_square_root_formula_'+str(digits),abs(mp.re(om)-mp.sqrt(a*(M+U)/(2*M*M)))<mp.mpf('1e-65'))
        check('positive_imaginary_square_root_formula_'+str(digits),abs(mp.im(om)-theta*mp.sqrt(a/(2*M*M*(M+U))))<mp.mpf('1e-65'))
        records.append(dict(digits=digits,theta=mp.nstr(theta,55),U=mp.nstr(U,55),
            Omega_real=mp.nstr(mp.re(om),55),Omega_imaginary=mp.nstr(mp.im(om),55),
            efold_time_over_sqrtG=mp.nstr(1/mp.im(om),55),scalar_mass_sqrtG=mp.nstr(mp.sqrt(m2),55),
            local_pole_to_scalar_mass_ratio=mp.nstr(abs(om)/mp.sqrt(m2),55),
            model_residual_abs=mp.nstr(abs(residual),8)))
    for key in ['theta','U','Omega_real','Omega_imaginary','efold_time_over_sqrtG','scalar_mass_sqrtG']:
        check('cross_precision_'+key,abs(mp.mpf(records[0][key])-mp.mpf(records[1][key]))<mp.mpf('1e-52'))
    # A deliberately wider rational bracket than the point root error;
    # interval signs certify it uniformly over all current u_b values.
    unit=mp.mpf(10)**28;nn=int(mp.floor(theta*unit))
    ls=mp.nstr(mp.mpf(nn-1)/unit,80);hs=mp.nstr(mp.mpf(nn+2)/unit,80)
    iv=mp.iv;iv.dps=90;ai,zi,mi,gi=constants(iv)
    zi=zi*iv.exp(iv.mpf(['-1.08e-36','1.08e-36']))
    lf=iv.log(ai/zi)
    fiv=lambda th:-th*iv.cos(th)/iv.sin(th)+iv.log(th/iv.sin(th))-lf
    check('interval_root_lower_all_ub',fiv(iv.mpf(ls))<0,str(fiv(iv.mpf(ls))))
    check('interval_root_upper_all_ub',fiv(iv.mpf(hs))>0,str(fiv(iv.mpf(hs))))
    ti=iv.mpf([ls,hs]);ui=-ti*iv.cos(ti)/iv.sin(ti);mm=iv.sqrt(ui**2+ti**2)
    re=iv.sqrt(ai*(mm+ui)/(2*mm**2));im=ti*iv.sqrt(ai/(2*mm**2*(mm+ui)))
    mag=iv.sqrt(ai/mm);ef=1/im
    check('interval_real_frequency',re>iv.mpf('.80') and re<iv.mpf('.82'),str(re))
    check('interval_growing_frequency',im>iv.mpf('.012') and im<iv.mpf('.014'),str(im))
    check('interval_frequency_modulus',mag>iv.mpf('.80') and mag<iv.mpf('.83'),str(mag))
    check('interval_W_modulus',mm>90 and mm<105,str(mm))
    check('interval_a_margin',ai>62,str(ai))
    check('interval_growth_time',ef>71 and ef<84,str(ef))
    negative=zi/(ai*iv.exp(1))
    check('interval_Euclidean_no_negative_zero',negative<iv.mpf('7.70e-45'),str(negative))
    bound=(iv.mpf('1e-6')*(iv.log(iv.mpf('1e-6')/zi)+iv.pi)+zi/iv.exp(1))/ai
    check('complex_model_low_invariant_band',bound<iv.mpf('1.4e-6'),str(bound))
    FIXTURES.append(dict(fixture='causal_pole_interval_and_independent_evaluations',records=records,
        theta_interval=[ls,hs],interval_digits=90,ub_absolute_allowance='5.4e-37',
        real_frequency_interval=str(re),imaginary_frequency_interval=str(im),growth_time_interval=str(ef),
        model_low_invariant_bound=str(bound),actual_boundary_operator_computed=False,
        physical_validity_band_proved=False))

def main():
    review=json.loads((ROOT/OWN[2]).read_bytes());old=json.loads((ROOT/'evidence'/(PREVIOUS+'.json')).read_bytes())['runs'][-1]
    paths=sorted(set(old['input_sha256_after'])|set(OWN)|set(review['retained_inputs'])|{'evidence/'+PREVIOUS+'.json','evidence/'+PREVIOUS+'.log'})
    before={p:digest(p) for p in paths};error=None
    run=dict(started_utc=datetime.now(timezone.utc).isoformat(),input_sha256_before=before,
        frozen_sources_utf8={p:(ROOT/p).read_text(encoding='utf-8') for p in OWN},previous_goal_turn='progress',protocol=review['protocol'])
    try:
        check('previous_full_trace_defect_run_passed',old['status']=='passed' and old['full_disk_trace_defect_bounded'])
        check('previous_inputs_unchanged',all(digest(p)==h for p,h in old['input_sha256_after'].items()))
        algebra_and_rational();numerical()
    except Exception:error=traceback.format_exc()
    after={p:digest(p) for p in paths}
    run.update(status='passed' if error is None and before==after else 'failed',exception=error,
        input_sha256_after=after,unchanged_inputs=before==after,checks=ROWS,fixtures=FIXTURES,
        finished_utc=datetime.now(timezone.utc).isoformat(),
        frozen_shape_causal_poles_enclosed=error is None,conditional_causal_pole_transfer_margin_proved=error is None,
        model_small_invariant_band_checked=error is None,actual_retarded_characteristic_reduction_proved=False,
        actual_causal_pole_transfer_bound_proved=False,actual_finite_geometry_growing_mode_proved=False,
        physical_instability_proved=False,physical_stability_proved=False,semiclassical_validity_proved=False,
        theory_or_boundary_prescription_changed=False,full_quantum_shape_defect_evaluated=False,
        full_defect_bolt_pairing_evaluated=False,nonlinear_solution_constructed=False,
        component_complete=False,goal_complete=False,formal_verification=False)
    out=ROOT/'evidence'/(STEM+'.json');receipt=json.loads(out.read_bytes()) if out.exists() else {'runs':[]}
    receipt['runs'].append(run);out.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    with (ROOT/'evidence'/(STEM+'.log')).open('a',encoding='utf-8',newline='\n') as log:
        for row in [dict(event='start',run=len(receipt['runs']),source_sha256={p:before[p] for p in OWN}),*ROWS,dict(event='finish',status=run['status'],exception=error)]:log.write(json.dumps(row)+'\n')
    print(json.dumps(dict(status=run['status'],checks=len(ROWS),inputs=len(paths),fixtures=len(FIXTURES),exception=error)))
    if run['status']!='passed':raise SystemExit(1)

if __name__=='__main__':main()
