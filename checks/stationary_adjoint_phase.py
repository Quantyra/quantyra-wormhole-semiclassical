"""Finite checks for the actual adjoint mass pairing and scalar-phase reduction.

The nonlocal cutoff limit and actual global phase quotient remain ordinary
proofs. No actual kernel, phase angle or candidate moment is computed.
"""
from __future__ import annotations
import base64,contextlib,hashlib,json,platform,sys,traceback
from datetime import datetime,timezone
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/stationary-adjoint-phase.json'
LOG=ROOT/'evidence/stationary-adjoint-phase.log'
INPUTS=[
    'checks/stationary_adjoint_phase.py',
    'notes/stationary-adjoint-phase.md',
    'notes/stationary-global-weighted-map.md',
    'notes/stationary-disk-and-normal-inverse.md',
    'notes/stationary-centered-complete-residual.md',
    'notes/stationary-centered-exterior-source.md',
    'notes/stationary-centered-covariant-seed.md',
    'notes/stationary-fast-phase-variation.md',
    'notes/stationary-weighted-conformal-range.md',
    'notes/stationary-trace-range-condition.md',
]
CONTRACT={
    'acceptance':'every symbolic residual exactly zero and rational inequality true',
    'numerical_tolerances':'none',
    'scope':'normal boundary pairings, signs, Ward/gauge identities and source integrability',
    'ordinary_arguments':['actual nonlocal cutoff limit','finite asymptotic extraction',
        'actual scalar-phase quotient','actual candidate source-tail proof'],
    'unresolved':['actual physical kernel','actual phase angle','candidate source moments',
        'quantitative inverse','nonlinear physical logarithmic phase','physical matching'],
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

def matrix_zero(label,expression):
    values=[S.simplify(x) for x in expression]
    check(label,all(x==0 for x in values),values)


def boundary_pairings():
    s=S.symbols('s')
    at,d0,aa,ee,lam,gg,h0=S.symbols('at d0 aa ee lam gg h0',positive=True)
    q=s*(s-3)
    hf=S.Function('hf')(s)
    nt=d0*q*(q+2)*hf-at*q
    jc=3*(aa+ee*lam*q)*(q-4)/(4*S.pi*gg*lam)
    def C(v):return S.Matrix([[1,(v+3)/(3*(v+1))],[0,2]])
    u0=S.Matrix([1,-1]);um=S.Matrix([1,-2]);es=S.Matrix([1,0])
    matrix_zero('TT_source_vector',C(S.Integer(0))*u0-S.Matrix([0,-2]))
    matrix_zero('TT_mass_vector',C(S.Integer(3))*um-S.Matrix([0,-4]))
    matrix_zero('scalar_vector',C(s)*es-es)
    H=C(3-s).T*S.diag(jc,nt)*C(s)
    pair=(u0.T*S.diff(H,s)*um)[0].subs(s,3)
    zero('reflected_congruence_mass_derivative',pair-8*S.diff(nt,s).subs(s,3))
    slope=S.diff(nt,s).subs(s,3).subs(hf.subs(s,3),h0)
    zero('TT_mass_root_slope',slope-3*(2*d0*h0-at))
    zero('positive_mass_pairing_coefficient',-8*slope-24*(at-2*d0*h0))
    check('loop_coefficient_bound',2*(1+R(54,5))<24)
    zero('relative_loop_bound',24/(30*S.Integer(2000)**10)-R(4,5)/S.Integer(2000)**10)
    check('benchmark_relative_loop_bound',R(4,5)/S.Integer(2000)**10<R(8,10**34))
    z=S.symbols('z',positive=True)
    def D(f):return z*S.diff(f,z)
    cut_mass=z**3*(1-z)**5
    # Smooth through the core for this second-order sign witness; same outer cubic.
    zero('local_outer_flux_sign',S.integrate(S.expand((D(D(cut_mass))-3*D(cut_mass))/z**4),(z,0,1))+3)
    r,beta=S.symbols('r beta',real=True)
    moment=-S.integrate(r*S.exp(-3*r),(r,0,S.oo))-S.integrate(r*S.exp(r),(r,-S.oo,0))
    zero('nonlocal_commutator_first_moment',moment-S.diff(4/(4-beta**2),beta).subs(beta,1))
    omega=S.symbols('omega',positive=True)
    sp=R(3,2)+S.I*omega
    K=3*omega*(aa+4*ee*lam)/(2*S.pi*gg*lam)
    root_sub={aa:ee*lam*(omega**2+R(9,4))}
    zero('scalar_root_derivative',S.expand((S.diff(jc,s).subs(s,sp)+S.I*K).subs(root_sub)))
    zero('scalar_flux_normalization',K.subs(root_sub).subs(ee,gg/(60*S.pi))-omega*(omega**2+R(25,4))/(40*S.pi**2))
    kval=S.symbols('kval',positive=True)
    complex_form=S.Matrix([[0,-S.I*kval],[S.I*kval,0]])
    change=S.Matrix([[R(1,2),1/(2*S.I)],[R(1,2),-1/(2*S.I)]])
    matrix_zero('real_scalar_Green_determinant',change.T*complex_form*change-S.Matrix([[0,kval/2],[-kval/2,0]]))


def covariance_and_source():
    x=S.symbols('x',positive=True)
    a=S.Function('a')(x);r=S.Function('r')(x)
    u=S.Function('u')(x);h=S.Function('h')(x);d=S.Function('d')(x)
    er=(u+d)/2;et=(u-d)/2
    W=a*a*r*r
    ward=S.diff(er,x)+S.diff(a,x)/a*(er-et)+2*S.diff(r,x)/r*(er-h)
    target=S.diff(W*d,x)+W*(S.diff(u,x)+2*S.diff(r,x)/r*u-4*S.diff(r,x)/r*h)
    zero('full_spherical_Ward_reconstruction',2*W*ward-target)
    fs,ft,G=S.symbols('fs ft G',real=True)
    U=8*S.pi*G*ft;H=4*S.pi*G*(fs-ft)
    zero('full_trace_output_reconstruction',(U+2*H)/(8*S.pi*G)-fs)
    zero('orbit_trace_output_reconstruction',U/(8*S.pi*G)-ft)
    eta=S.Function('eta')(x)
    zero('radial_isothermal_gauge_equation',S.diff(eta/x,x)-(S.diff(eta,x)-eta/x)/x)
    c=S.symbols('c',real=True)
    zero('radial_dilation_solution',S.diff(c*x,x)-c*x/x)
    kap,B,lam=S.symbols('kap B lam',positive=True)
    zero('thermal_Phi_source',S.diff(S.log(kap/(B*lam)),kap)-1/kap)
    zero('thermal_V_source',S.diff(S.log(B*S.sqrt(lam)/kap),kap)+1/kap)
    z=S.symbols('z',positive=True)
    for j in range(3):
        zero('candidate_weighted_tail_'+str(j),S.integrate(z**8*z**(-7)*S.log(z)**(2*j)/z,(z,0,1))-S.factorial(2*j))
    zero('scalar_moment_tail',S.integrate(z**R(3,2)*z**4*z**(-4),(z,0,1))-R(2,5))
    check('source_mass_corrections_decay',2+3-3>0 and R(3,2)+3-3>0 and R(7,2)+0-3>0)


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
            boundary_pairings()
            covariance_and_source()
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
