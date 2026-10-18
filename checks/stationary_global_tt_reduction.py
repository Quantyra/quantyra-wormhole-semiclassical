"""Exact finite checks for the global TT projection and classical seed gap.

No quantum form bound, actual kernel, adjoint angle or moment is computed.
Analytic mapping, packet, Schur and bolt-closure arguments remain ordinary.
"""
from __future__ import annotations
import base64,contextlib,hashlib,json,platform,sys,traceback
from datetime import datetime,timezone
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/stationary-global-tt-reduction.json'
LOG=ROOT/'evidence/stationary-global-tt-reduction.log'
INPUTS=[
    'checks/stationary_global_tt_reduction.py',
    'notes/stationary-global-tt-reduction.md',
    'notes/stationary-adjoint-phase.md',
    'notes/stationary-global-weighted-map.md',
    'notes/stationary-disk-and-normal-inverse.md',
    'notes/stationary-centered-covariant-seed.md',
    'notes/stationary-centered-exterior-source.md',
    'notes/stationary-shape-logarithmic-response.md',
]
CONTRACT={
    'acceptance':'every symbolic residual exactly zero and rational inequality true',
    'numerical_tolerances':'none',
    'scope':'normal diagonal, York identities, classical TT form and exact gap guards',
    'ordinary_arguments':['escaping packets','global weighted York inverse',
        'finite regularity qualification','Hardy and Schur bounds','smooth bolt closure'],
    'unresolved':['quantum lower-form bound','candidate sharp mixed projection',
        'off-shell transfer','physical kernel','adjoint phase','candidate moments',
        'quantitative coupled inverse','nonlinear and physical matching'],
    'failure_policy':'append every run and retain failed input bytes',
    'hash_policy':'SHA256 of raw UTF8/LF/noBOM bytes, no normalization',
}
ROWS=[]
R=S.Rational


def check(label,condition,detail=''):
    row={'label':label,'passed':bool(condition),'detail':str(detail)}
    ROWS.append(row)
    print(json.dumps(row,sort_keys=True),flush=True)


def zero(label,expression):
    value=S.simplify(S.expand_func(expression))
    check(label,value==0,value)


def normal_and_york():
    rr,g,H=S.symbols('R g H',positive=True)
    q=-rr
    jc=3*(1+11*g/(180*S.pi)+g*q/(60*S.pi))*(q-4)/(4*S.pi*g)
    nt=q*(q+2)*H/(480*S.pi**2)-(1/(48*S.pi*g)+R(23,8640)/S.pi**2)*q
    expected=(rr-4)/(12*S.pi*g)+(35*rr-44)/(2160*S.pi**2)+rr*(rr-2)*H/(120*S.pi**2)
    zero('exact_fixed_theta_normal',4*jc/(9*(rr+4))+4*nt-expected)
    zero('negative_endpoint_formula',expected.subs(rr,R(9,4))-(-7/(48*S.pi*g)+R(139,8640)/S.pi**2+3*H/(640*S.pi**2)))
    zero('positive_endpoint_formula',expected.subs(rr,4)-(R(2,45)+H/15)/S.pi**2)
    check('negative_endpoint_loop_guard',R(139,8640)+R(33,640)<R(1,10) and 70*3>48)
    check('positive_H_guard',R(39,5)-3-2==R(14,5)>0)
    check('root_shift_loop_guard',R(96,2160)+R(8*74,5*120)<2)
    check('benchmark_root_scale',8*R(1,2000**10)<R(1,10**30))
    rho=S.symbols('rho',positive=True)
    aa=S.Function('aa')(rho);rrr=S.Function('rrr')(rho);y=S.Function('y')(rho)
    w=aa+2*rrr
    p0=3*S.diff(y,rho)/2-w*y/2
    p1=-S.diff(y,rho)/2+(3*aa/2-rrr)*y
    p2=-S.diff(y,rho)/2+(rrr-aa/2)*y
    vv=S.diff(w,rho)/2+2*aa**2+4*rrr**2
    zero('York_tracefree',p0+p1+2*p2)
    zero('York_divergence',S.diff(p0,rho)+aa*(p0-p1)+2*rrr*(p0-p2)-(3*(S.diff(y,rho,2)+w*S.diff(y,rho))/2-vv*y))
    zero('York_energy_by_parts',p0**2+p1**2+2*p2**2-3*S.diff(y,rho)**2-2*vv*y**2+S.diff(w*y**2,rho)+w**2*y**2)
    zero('York_kernel_difference_01',p0-p1-2*(S.diff(y,rho)-aa*y))
    zero('York_kernel_difference_02',p0-p2-2*(S.diff(y,rho)-rrr*y))
    ss,lam,tt,sc=S.symbols('s lam T S',positive=True)
    pn=R(3,2)*lam*(ss**2-3*ss)-6*lam
    zero('York_outer_roots',pn-R(3,2)*lam*(ss+1)*(ss-4))
    yn=-2*tt/(3*S.sqrt(lam)*(ss+1))
    zero('York_normal_solution',pn*yn-S.sqrt(lam)*(4-ss)*tt)
    zero('York_normal_conformal',sc+tt/2-(3-ss)*S.sqrt(lam)*yn/4-(sc+(ss+3)*tt/(3*(ss+1))))
    check('York_outer_essential_threshold',R(3,2)*R(9,4)+6==R(75,8)>0)
    check('York_bolt_inverse_square',R(3,2)-R(3,8)==R(3,2)*R(3,4))
    check('York_seed_gap_guard',R(3,2)-R(13,2000)>R(149,100))


def seed_and_tt():
    r=S.symbols('r',positive=True)
    B,lam,m0,q0,d=S.symbols('B lam m0 q0 d',real=True)
    seed=1-2*m0/r+q0/r**2+lam*r**2+d*(B**2/r**2-2*B**3/r**3+B**4/r**4)
    D=lambda v,n=1:S.diff(v,r,n)
    zero('York_seed_potential',D(seed,2)/4+D(seed)/(2*r)-(3*lam/2+q0/(2*r**4)+d*(B**2/(2*r**4)-3*B**3/r**5+3*B**4/r**6)))
    k0=D(seed)+2*seed/r
    zero('TT_Schur_monotonicity',D(k0)-(4*lam-2/r**2+4*m0/r**3+d*(-8*B**3/r**5+10*B**4/r**6)))
    F=S.Function('F')(r);c=S.Function('c')(r);b=S.Function('b')(r);hh=S.Function('H')(r)
    bp=-(2/r+D(F)/F)*b-(D(F)/F-2/r)*c
    df=-F*b+F*(c+r*D(c))-r*D(F)*c/2
    L=3*r*F+r**2*D(F)/2
    K=-F+r*D(F)-r**2*D(F,2)/2
    M=-2*r*F+r**2*D(F)/2+r**3*D(F,2)/2
    force=r**2*F*D(c,2)+L*D(c)+K*c+F*b
    zero('areal_TT_force',(D(r*df)-force).subs(D(b),bp))
    vf=6*F-7*r*D(F)/2+r**2*D(F,2)/2-r**3*D(F,3)/4
    J=-3*r**2*F*c*D(c)-r**3*F*D(c)**2/2-r*F*b*c+M*c**2/2
    zero('complete_TT_integration_certificate',(-(3*c+r*D(c))*force-(3*r**2*F*D(c)**2/2+vf*c**2-4*F*b*c+D(J))).subs(D(b),bp))
    wf=3*F-r*D(F)/2+r**2*D(F,2)/2-r**3*D(F,3)/4
    zero('TT_square_completion',3*r**2*F*D(c)**2/2+vf*c**2-(3*r**2*F*(D(c)+2*c/r)**2/2+wf*c**2)+D(3*r*F*c**2))
    wseed=wf.subs(F,seed).doit()
    zero('TT_seed_potential',wseed-(3-12*m0/r+13*q0/r**2+3*lam*r**2+d*(13*B**2/r**2-51*B**3/r**3+45*B**4/r**4)))
    p=r**4*D(F)-2*r**3*F
    primitive=4*r**2*F*b*hh-2*p*hh**2
    zero('TT_cross_with_bolt_primitive',(-4*F*b*c-2*D(p)*hh**2-D(primitive)).subs({D(b):bp,D(hh):-c/r**2}))
    zero('TT_cross_seed_guard',D(p.subs(F,seed).doit())-(-6*r**2+12*m0*r-4*q0+d*(-4*B**2+6*B**4/r**2)))
    zero('Hardy_boundary_identity',D(r**3*hh**2).subs(D(hh),-c/r**2)-(3*r**2*hh**2-2*r*hh*c))
    u=S.Function('U')(r)
    zero('Schur_first_primitive',D(2*S.sqrt(u))-D(u)/S.sqrt(u))
    zero('Schur_second_primitive',D(-2/S.sqrt(u))-D(u)/u**R(3,2))
    t=S.symbols('t',positive=True)
    zero('F_second_x_minimum',2-8/t**3+18/t**4-R(52,27)-R(2,27)*(t-3)**2*(t**2+6*t+27)/t**4)
    zero('F_second_constant_minimum',-4/t**3+6/t**4+R(1,8)-(t-2)**2*(t**2+4*t+12)/(8*t**4))
    zero('W_constant_minimum',3-12/t+13/t**2-(13*(1/t-R(6,13))**2+R(3,13)))
    check('W_x_calculus_guard',6<2**3 and 18*2<39)
    check('F_second_positive_guard',R(52*1000,27)-R(1,8)-2-50>0)
    check('Schur_monotone_guard',4*1000-2-18>0)
    check('Hardy_cross_guard',8+12*1000-10>0)
    check('W_complete_gap_guard',1000+R(3,13)-6-109-R(16,3)>0)
    check('TT_tensor_norm_constant',4**2+2**2+2==22)
    check('TT_physical_gap_constant',2*22*4==176)
    check('coarse_thermal_guard',16*2000**3*R(1,10**60)<1)
    check('projection_strip_guards',-1<R(3,2)<4 and -1<3<4 and -1<R(7,2)<4)


def hashes(data):
    return {name:hashlib.sha256(raw).hexdigest() for name,raw in data.items()}


def main():
    before={name:(ROOT/name).read_bytes() for name in INPUTS}
    existing=json.loads(RECEIPT.read_text(encoding="utf-8")) if RECEIPT.exists() else {"runs":[]}
    run={"started_utc":datetime.now(timezone.utc).isoformat(),"contract":CONTRACT,
         "command":sys.orig_argv,"dont_write_bytecode":sys.dont_write_bytecode,"cwd":str(Path.cwd()),
         "python":platform.python_version(),"sympy":S.__version__,"execution_inputs_sha256":hashes(before)}
    with LOG.open("a",encoding="utf-8",newline="\n") as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
        print("RUN "+json.dumps(run,sort_keys=True),flush=True)
        try:
            check("actual_bytecode_disabled",sys.dont_write_bytecode and "-B" in sys.orig_argv)
            check("inputs_strict_UTF8_LF",all(raw.decode("utf-8") is not None and b"\r" not in raw
                  and not raw.startswith(b"\xef\xbb\xbf") for raw in before.values()))
            normal_and_york()
            seed_and_tt()
        except Exception:
            traceback.print_exc()
            check("uncaught_exception",False)
        run["final_inputs_sha256"]=hashes({name:(ROOT/name).read_bytes() for name in INPUTS})
        print("FINAL_INPUTS "+json.dumps(run["final_inputs_sha256"],sort_keys=True),flush=True)
        check("frozen_inputs_unchanged",run["execution_inputs_sha256"]==run["final_inputs_sha256"])
        run["checks"]=ROWS
        run["failed_checks"]=[row for row in ROWS if not row["passed"]]
        code=int(bool(run["failed_checks"]))
        run["status"]="failed" if code else "passed"
        run["terminal_exit_status"]=code
        run["finished_utc"]=datetime.now(timezone.utc).isoformat()
        if code:
            run["failed_inputs_base64"]={name:base64.b64encode(raw).decode("ascii") for name,raw in before.items()}
            for name,raw in before.items():
                print(json.dumps({"path":name,"sha256":hashlib.sha256(raw).hexdigest(),
                                  "base64":base64.b64encode(raw).decode("ascii")},sort_keys=True))
        print("TERMINAL "+json.dumps({"status":run["status"],"checks":len(ROWS),
              "failures":len(run["failed_checks"]),"terminal_exit_status":code,
              "finished_utc":run["finished_utc"]},sort_keys=True),flush=True)
    existing["runs"].append(run)
    RECEIPT.write_text(json.dumps(existing,indent=2,sort_keys=True)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps({"status":run["status"],"checks":len(ROWS),"failures":len(run["failed_checks"]),
                      "receipt":str(RECEIPT)},sort_keys=True))
    return code


if __name__=="__main__":
    sys.exit(main())
