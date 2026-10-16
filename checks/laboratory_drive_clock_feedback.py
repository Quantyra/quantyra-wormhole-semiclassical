"""S135.3 one frozen laboratory waveform; coupled free clock diagnostic."""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_name]='1'
import hashlib
import argparse
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.integrate import solve_ivp,quad
from scipy.optimize import brentq
import sympy as sym

ROOT=Path(__file__).resolve().parents[1]
ETA=1/2000
TAU=1/4
PI=math.pi
TSTAR=PI+TAU
HORIZON=8.
PHASE_MAX=16.
ZSTOP=1/200
REF_INIT=53/200
STATE_TOL=2e-6
IDENTITY_TOL=2e-8
BOUND_ALLOWANCE=1e-8
INPUT_PATHS=['notes/laboratory-drive-clock-feedback.md','checks/laboratory_drive_clock_feedback.py','notes/returned-switch-capture.md','notes/boundary-vertex-matching.md','notes/schwarzian-initial-state.md','notes/jt-prescribed-source-constraint.md','notes/prescribed-source-clock-variation.md']


def S(t):
    x=t/TAU
    if x<=0:return 0.
    if x>=1:return 1.
    return x**3*(10+x*(-15+6*x))


def integral_S(t):
    x=t/TAU
    if x<=0:return 0.
    if x>=1:return t-TAU/2
    return TAU*(2.5*x**4-3*x**5+x**6)


def integral_S_squared(t):
    x=t/TAU
    if x<=0:return 0.
    total=181/462
    if x>=1:return t-TAU+TAU*total
    return TAU*(100*x**7/7-75*x**8/2+115*x**9/3-18*x**10+36*x**11/11)


def theta(t):return 4*ETA*(integral_S(t)-integral_S(t-PI))
def rho0(t):return 4*ETA*(S(t)-S(t-PI))
def d0(t):return 4*ETA*S(t)
def Q0(t):return 16*ETA**2*(integral_S_squared(t)-integral_S_squared(t-PI))
def g0(t):return (theta(t)-Q0(t))/(4*PI*ETA)


def reference(method):
    result=solve_ivp(lambda t,y:[y[1],-y[0]+g0(t)/2],(0,PHASE_MAX),[REF_INIT,0.],method=method,dense_output=True,rtol=2e-13,atol=2e-15,max_step=.005)
    if not result.success:raise AssertionError(result.message)
    return result


def algebra():
    a,b,w,wp,wpp,z,p,q,r,eta=sym.symbols('a b w wp wpp z p q r eta',nonzero=True)
    # Product-rule identities for aligned variables and inverse phase.
    p=w/z
    adot=b+wp*a/z
    bdot=(-1+wpp/z)*a-q/(8*sym.pi*eta)+r
    for value in (adot/p-(b+wp*a/w+a*b/w),bdot/p-((-1+wpp/w)*a-a*a/w+(1+a/w)*(r-q/(8*sym.pi*eta)))):
        assert sym.simplify(value.subs(z,w+a))==0
    x=sym.symbols('x')
    polynomial=10*x**3-15*x**4+6*x**5
    assert sym.integrate(polynomial**2,(x,0,1))==sym.Rational(181,462)
    assert sym.diff(polynomial,x).subs(x,sym.Rational(1,2))==sym.Rational(15,8)
    return {'aligned_and_inverse_phase_identities':'exact symbolic pass','quintic_squared_integral':'181/462','d0_Lipschitz_constant':'30*eta'}


class Evolution:
    def __init__(self,zinit,vinit,ref,aligned=False,refined=False):
        self.ref=ref
        self.aligned=aligned
        self.pieces=[]
        self.initial=(zinit,vinit)
        self.completions=[]
        initial=[zinit-REF_INIT,vinit,0.,0.] if aligned else [zinit,vinit,0.,0.]
        def stop(t,y):return self.unpack(y)[0]-ZSTOP
        stop.terminal=True
        stop.direction=-1
        def phase_stop(t,y):return PHASE_MAX-y[2]
        phase_stop.terminal=True
        phase_stop.direction=-1
        def complete(t,y):return y[2]-TSTAR
        complete.terminal=False
        complete.direction=1
        method='RK45' if aligned else 'DOP853'
        tolerance=(2e-12,2e-14,.005) if refined or aligned else (2e-11,2e-13,.01)
        for start,end in ((0.,PI),(PI,2*PI),(2*PI,HORIZON)):
            result=solve_ivp(self.rhs,(start,end),initial,method=method,dense_output=True,rtol=tolerance[0],atol=tolerance[1],max_step=tolerance[2],events=(stop,phase_stop,complete))
            if not result.success:raise AssertionError(result.message)
            self.pieces.append((start,float(result.t[-1]),result))
            self.completions.extend(map(float,result.t_events[2]))
            initial=result.y[:,-1]
            if result.t[-1]<end:break
        self.end=self.pieces[-1][1]

    def unpack(self,y):
        if not self.aligned:return np.asarray(y)
        a,b,F,Q=y
        w,wp=self.ref.sol(F)
        return np.array([a+w,b+wp,F,Q])

    def state(self,t):
        for start,end,result in reversed(self.pieces):
            if start-1e-10<=t<=end+1e-10:return self.unpack(result.sol(t))
        raise AssertionError(('missing past history',t))

    def pastD(self,t):
        value=0.
        while t>1e-12:
            z,v,F,Q=self.state(t)
            value+=rho0(F)*self.ref.sol(F)[0]/z
            t-=PI
        return value

    def rhs(self,t,y):
        z,v,F,Q=self.unpack(y)
        w,wp=self.ref.sol(F)
        phase_rate=w/z
        rate=rho0(F)*phase_rate
        delayed=self.pastD(t-PI) if t>=PI else 0.
        qprime=rate*rate+2*rate*delayed
        if not self.aligned:return [v,-z+(theta(F)-Q)/(8*PI*ETA),phase_rate,qprime]
        a,b=y[:2]
        wpp=-w+g0(F)/2
        return [b+wp*a/z,(-1+wpp/z)*a-(Q-Q0(F))/(8*PI*ETA),phase_rate,qprime]

    def energy(self,t):
        z,v,F,Q=self.state(t)
        force=(theta(F)-Q)/(4*PI*ETA)
        return v*v+z*z-force*z


def error_checks(evolution,grid):
    # Reconstruct memory independently as a moving-window integral.
    memory_errors=[]
    defect_checks=[]
    for t in np.linspace(.05,evolution.end,33):
        z,v,F,Q=evolution.state(t)
        lo=max(0,t-PI)
        points=[lo,t]
        for j in (TAU,PI,PI+TAU,2*PI):
            if lo<j<t:points.append(j)
        points=sorted(set(points))
        integrated=sum(quad(lambda s:evolution.pastD(s)**2,points[i],points[i+1],epsabs=1e-12,epsrel=1e-10,limit=150)[0] for i in range(len(points)-1))
        memory_errors.append(abs(Q-integrated))
        alpha=lambda s:abs(1-evolution.ref.sol(evolution.state(s)[2])[0]/evolution.state(s)[0]) if s>0 else 0.
        summed=0.
        shifted=t
        while shifted>=0:
            summed+=alpha(shifted)
            shifted-=PI
        integral_alpha=quad(alpha,0,t,epsabs=1e-9,epsrel=1e-8,limit=150)[0]
        bound=4*ETA*summed+30*ETA*integral_alpha
        actual=abs(evolution.pastD(t)-d0(F))
        if actual>bound+BOUND_ALLOWANCE:raise AssertionError('radiation retiming defect inequality failed')
        defect_checks.append({'t':float(t),'actual_D_defect':float(actual),'analytic_bound':float(bound)})
    if max(memory_errors)>IDENTITY_TOL:raise AssertionError(('memory identity failed',max(memory_errors)))
    completion=evolution.completions[0] if evolution.completions else None
    energy_drift=None
    if completion is not None and completion<evolution.end:
        energies=np.array([evolution.energy(t) for t in np.linspace(completion,evolution.end,101)])
        energy_drift=float(np.ptp(energies))
        if energy_drift>IDENTITY_TOL:raise AssertionError(('post-completion energy not constant',energy_drift))
    # Integrate the finite energy/work identity, including returned cross terms.
    def work(t):
        z,v,F,Q=evolution.state(t)
        p=evolution.ref.sol(F)[0]/z
        rate=rho0(F)*p
        delayed=evolution.pastD(t-PI) if t>=PI else 0.
        gprime=rate*(1-rate-2*delayed)/(4*PI*ETA)
        return -gprime*z
    work_integral=quad(work,0,evolution.end,epsabs=1e-10,epsrel=1e-9,limit=250)[0]
    work_error=abs(evolution.energy(evolution.end)-evolution.energy(0)-work_integral)
    if work_error>IDENTITY_TOL:raise AssertionError(('energy/work identity failed',work_error))
    return {'max_memory_integral_error':float(max(memory_errors)),'post_completion_energy_drift':energy_drift,'energy_work_error':float(work_error),'sampled_radiation_defect_checks':defect_checks}


def completion_reporting(record,ref):
    """Reporting from an accepted endpoint; no feedback trajectory is rerun."""
    zinit,vinit=record['initial']
    w,wp=ref.sol(TSTAR)
    zbase=w+(zinit-REF_INIT)*math.cos(TSTAR)+vinit*math.sin(TSTAR)
    vbase=wp-(zinit-REF_INIT)*math.sin(TSTAR)+vinit*math.cos(TSTAR)
    baseline=float(vbase*vbase+zbase*zbase-g0(TSTAR)*zbase)
    delta=float(record['completion_energy']-baseline)
    record.update(same_preparation_global_baseline_completion_energy=baseline,
                  matched_source_completion_energy_difference=delta,
                  old_absolute_energy_correction_one_over_50_budget=bool(abs(delta)<1/50))
    record['comparison_time_windows']={'before_Tstar':'old fixed global-time [0,Tstar] window','after_Tstar':'old fixed global-time [Tstar,8] window; not necessarily the actual post-source interval'}


def summarize(evolution,ref):
    grid=np.linspace(0,evolution.end,8001)
    values=np.array([evolution.state(t) for t in grid]).T
    central=ref.sol(grid)
    zinit,vinit=evolution.initial
    sameprep=central[0]+(zinit-REF_INIT)*np.cos(grid)+vinit*np.sin(grid)
    before=grid<=TSTAR
    after=grid>=TSTAR
    completion=evolution.completions[0] if evolution.completions else None
    output={'initial':[zinit,vinit],'terminal_time':evolution.end,'domain_exit':bool(evolution.end<HORIZON),'sampled_min_z':float(values[0].min()),'source_completion_time':completion,'source_completion_shift':None if completion is None else completion-TSTAR,'max_phase_shift':float(np.max(abs(values[2]-grid))),'common_time_central_error_before_Tstar':float(np.max(abs(values[0,before]-central[0,before]))),'common_time_same_preparation_error_before_Tstar':float(np.max(abs(values[0,before]-sameprep[before]))),'common_time_central_error_after_Tstar':float(np.max(abs(values[0,after]-central[0,after]))),'common_time_same_preparation_error_after_Tstar':float(np.max(abs(values[0,after]-sameprep[after])))}
    output['old_common_time_one_over_200_clock_budget']=bool(output['common_time_same_preparation_error_before_Tstar']<1/200)
    output['old_common_time_post_one_over_50_clock_budget']=bool(output['common_time_same_preparation_error_after_Tstar']<1/50)
    if completion is not None:
        z,v,F,Q=evolution.state(completion)
        gf=(theta(F)-Q)/(4*PI*ETA)
        energy=evolution.energy(completion)
        capture=bool(gf>0 and energy<0)
        future=gf/2-math.sqrt(max(0,gf*gf/4+energy)) if capture else None
        output.update(final_g=float(gf),completion_energy=float(energy),negative_energy_capture_diagnostic=capture,autonomous_future_min_z=None if future is None else float(future),old_minus_one_over_25_energy_margin=bool(energy<-1/25))
        completion_reporting(output,ref)
    return output,grid,values


def main():
    identities=algebra()
    ref=reference('DOP853')
    independent_ref=reference('RK45')
    refgrid=np.linspace(0,PHASE_MAX,16001)
    referror=float(np.max(abs(ref.sol(refgrid)-independent_ref.sol(refgrid))))
    if referror>IDENTITY_TOL:raise AssertionError(('reference disagreement',referror))
    cases=[]
    for zinit in (.25,.28):
        for vinit in (-.0005,0.,.0005):
            primary=Evolution(zinit,vinit,ref)
            refined=Evolution(zinit,vinit,ref,refined=True)
            independent=Evolution(zinit,vinit,independent_ref,aligned=True)
            record,grid,values=summarize(refined,ref)
            comparison_grid=np.linspace(0,min(primary.end,refined.end,independent.end),4001)
            fine=np.array([refined.state(t) for t in comparison_grid])
            coarse=np.array([primary.state(t) for t in comparison_grid])
            alternative=np.array([independent.state(t) for t in comparison_grid])
            refinement=float(np.max(abs(fine-coarse)))
            independent_error=float(np.max(abs(fine-alternative)))
            if max(refinement,independent_error)>STATE_TOL:raise AssertionError(('evolution disagreement',refinement,independent_error))
            completion_errors=[abs(run.completions[0]-refined.completions[0]) for run in (primary,independent)] if refined.completions and primary.completions and independent.completions else []
            energy_errors=[abs(run.energy(run.completions[0])-record['completion_energy']) for run in (primary,independent)] if completion_errors else []
            if completion_errors and max(completion_errors+energy_errors)>STATE_TOL:raise AssertionError('completion or energy disagreement')
            record['validation']={'uniform_refinement_error':refinement,'uniform_independent_aligned_RK45_error':independent_error,'completion_time_errors':completion_errors,'completion_energy_errors':energy_errors,**error_checks(refined,grid)}
            cases.append(record)
            print(json.dumps({key:record[key] for key in ('initial','source_completion_time','completion_energy','old_minus_one_over_25_energy_margin','common_time_same_preparation_error_before_Tstar','common_time_same_preparation_error_after_Tstar')}),flush=True)
    # Identity control is verification detail, not the six-corner scientific test.
    control=Evolution(REF_INIT,0.,ref,refined=True)
    controlgrid=np.linspace(0,HORIZON,4001)
    controlstates=np.array([control.state(t) for t in controlgrid]).T
    expected=np.vstack((ref.sol(controlgrid),controlgrid,np.array([Q0(t) for t in controlgrid])))
    controlerror=float(np.max(abs(controlstates-expected)))
    if controlerror>IDENTITY_TOL:raise AssertionError(('exact-reference control failed',controlerror))
    paths=INPUT_PATHS
    result={'status':'pass','milestone':'E003/S135.3','scientific_verdict':'numerical capture at six corners; old quantitative margins and fixed completion time do not transfer to one laboratory waveform','frozen_model':{'eta':ETA,'reference_initial':[REF_INIT,0.],'horizon':HORIZON,'phase_domain_max':PHASE_MAX,'z_domain_stop':ZSTOP,'global_reference_source_completion':TSTAR,'residual':'r=0; no nonzero residual radius selected'},'symbolic':identities,'reference_independent_error':referror,'same_initial_identity_control_error':controlerror,'cases':cases,'tolerances':{'uniform_state_completion_energy':STATE_TOL,'memory_work_identity':IDENTITY_TOL,'defect_bound_rounding':BOUND_ALLOWANCE},'claim_limits':{'interval_certification':False,'whole_preparation_box_verified':False,'nonzero_residual_range_verified':False,'old_six_corner_quantitative_theorem_transferred':False,'four_dimensional_opening_claimed':False,'S074_complete':False},'inputs_sha256':{path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},'runtime':{'python':sys.version,'platform':platform.platform(),'numpy':np.__version__,'scipy':scipy.__version__,'sympy':sym.__version__}}
    path=ROOT/'evidence/laboratory-drive-clock-feedback.json'
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print('S135.3 one laboratory waveform: PASS numerical verification; inherited quantitative margin transfer FAILS.',flush=True)
    print('No waveform refitting, preparation tuning, interval theorem or four-dimensional opening claim.',flush=True)
    print('Receipt: '+str(path),flush=True)


def refresh_reporting():
    path=ROOT/'evidence/laboratory-drive-clock-feedback.json'
    result=json.loads(path.read_text(encoding='utf-8'))
    assert result['status']=='pass'
    assert result['frozen_model']['reference_initial']==[REF_INIT,0.]
    assert result['frozen_model']['eta']==ETA
    assert result['frozen_model']['horizon']==HORIZON
    original_checker=result.get('reporting_only_refresh',{}).get('evolution_run_checker_sha256',result['inputs_sha256']['checks/laboratory_drive_clock_feedback.py'])
    ref=reference('DOP853')
    for record in result['cases']:
        completion_reporting(record,ref)
        print(json.dumps({key:record[key] for key in ('initial','same_preparation_global_baseline_completion_energy','matched_source_completion_energy_difference','old_absolute_energy_correction_one_over_50_budget')}),flush=True)
    result['reporting_only_refresh']={'evolution_run_checker_sha256':original_checker,'scope':'Retain accepted trajectories/endpoints and verification; independently recompute same-preparation baseline at Tstar and matched-source-completion energy differences; add final reviewed note links/hashes. No evolution or tolerances changed.'}
    result['inputs_sha256']={item:hashlib.sha256((ROOT/item).read_bytes()).hexdigest() for item in INPUT_PATHS}
    path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    failed=sum(not record['old_absolute_energy_correction_one_over_50_budget'] for record in result['cases'])
    print('Reporting-only refresh PASS; '+str(failed)+' of '+str(len(result['cases']))+' matched-completion |delta e|<1/50 budgets fail.',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh-report',action='store_true',help='Recompute endpoint reporting from an existing passing receipt; do not rerun trajectories.')
    args=parser.parse_args()
    refresh_reporting() if args.refresh_report else main()
