"""Finite checks supporting the actual weighted stationary Fredholm theorem.

Green/contact estimates, Fredholm gluing and finite Mellin range-jump
construction remain ordinary proofs. No actual adjoint mode is computed.
"""
from __future__ import annotations
import base64,contextlib,hashlib,json,platform,sys,traceback
from datetime import datetime,timezone
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
RECEIPT=ROOT/'evidence/stationary-global-weighted-map.json'
LOG=ROOT/'evidence/stationary-global-weighted-map.log'
INPUTS=[
    'checks/stationary_global_weighted_map.py',
    'notes/stationary-global-weighted-map.md',
    'notes/stationary-disk-and-normal-inverse.md',
    'notes/stationary-shape-logarithmic-response.md',
    'notes/stationary-fast-phase-variation.md',
    'notes/stationary-weighted-conformal-range.md',
    'notes/stationary-centered-covariant-seed.md',
    'notes/stationary-paired-dirac-phase.md',
    'notes/stationary-exterior-conformal-domain.md',
    'notes/stationary-exterior-boundary-weyl.md',
]
CONTRACT={
    'acceptance':'every symbolic residual exactly zero and rational inequality true',
    'numerical_tolerances':'none',
    'scope':'finite density, kernel integration, strip-bound and index identities',
    'ordinary_arguments':['actual Green and contact expansion','dyadic operator bounds',
        'Fredholm parametrix','finite Mellin relative-index construction','actual range theorem'],
    'unresolved':['actual kernel and adjoint projections','quantitative inverse norm',
        'nonlinear physical logarithmic phase','exterior/state matching','physical gates'],
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

def weights_and_kernels():
    z,y,a,h,t=S.symbols('z y a h t',positive=True)
    delta=S.symbols('delta',real=True)
    u=S.Function('u')(z)
    zero('physical_Euler_adjoint',-S.diff(z*z**(-4)*u,z)/z**(-4)-(3*u-z*S.diff(u,z)))
    zero('physical_vs_b_weight',z**(-4)*z**(-4)-z**(-7)/z)
    phi,V,ds,dt,eps=S.symbols('phi V ds dt eps',real=True)
    density=S.exp(4*(phi+eps*ds)+2*(V+eps*dt))
    zero('fixed_volume_density_contact',S.diff(density,eps).subs(eps,0)/density.subs(eps,0)-4*ds-2*dt)
    zero('complementary_weight',(delta+(3-delta))-3)
    zero('Green_insertion_degree',S.Integer(-3)+2-1-3+4-(-1))
    zero('stress_kernel_difference_degree',S.Integer(-3)-1-2-(-6))
    zero('radial_kernel_difference_degree',S.Integer(-6)+3-(-3))
    for j in range(3):
        zero('contact_weight_gain_'+str(j),S.Integer(4)-(2-j)-j-2)
    r=S.symbols('r',positive=True)
    zero('tangential_integral_normal',4*S.pi*S.integrate(r**2/(a*a+r*r)**4,(r,0,S.oo))-S.pi**2/(8*a**5))
    zero('tangential_integral_remainder',4*S.pi*S.integrate(r**2/(a*a+r*r)**3,(r,0,S.oo))-S.pi**2/(4*a**3))
    k=z**(4-delta)*y**(1+delta)/(z+y)**5
    scaled=S.powsimp(k.subs(y,z*t))
    zero('weighted_image_kernel_scale',S.powsimp(scaled)-t**(1+delta)/(1+t)**5)
    zero('outer_Schur_integral',S.integrate(t**R(7,2)/(1+t)**5,(t,0,S.oo))-35*S.pi/128)
    zero('dual_Schur_integral',S.integrate(t**(-R(1,2))/(1+t)**5,(t,0,S.oo))-35*S.pi/128)
    zero('mass_far_integral',z**4*S.integrate(S.Integer(1),(y,2*z,h))-z**4*(h-2*z))
    zero('mass_local_integral',z**4*z**(-3)*S.integrate(y**3,(y,0,z/2))-z**5/64)
    zero('mass_output_weight_integral',S.integrate(z**8*z**(-7)/z,(z,0,h))-h)
    rate=S.symbols('rate',positive=True)
    for j in range(3):
        zero('weighted_log_integral_'+str(j),S.integrate(t**(2*j)*S.exp(-rate*t),(t,0,S.oo))-S.factorial(2*j)/rate**(2*j+1))
    check('two_weight_margins',4-R(7,2)==R(1,2) and 1+R(7,2)==R(9,2))

def strip_guards():
    a,t,y,d=S.symbols('a t y d',positive=True)
    summand=t/(y*y+t*t)
    zero('digamma_summand_monotonicity',S.diff(summand,y)+2*t*y/(y*y+t*t)**2)
    zero('digamma_integral_derivative',S.diff(S.pi/2-S.atan(a/t),a)+t/(a*a+t*t))
    zero('digamma_series_guard',1/R(1,2)**2+1/R(1,2)-6)
    check('digamma_value_lower_guard',R(3,5)+2*R(7,10)==2)
    check('quadratic_product_guard',(R(5,2)**2+4)*(R(3,2)**2+4)<81)
    zero('cot_first_pole_guard',2*R(25,16)/R(3,16)-R(50,3))
    zero('cot_tail_guard',2*R(25,16)*R(64,39)-R(200,39))
    check('cot_total_guard',1+R(50,3)+R(200,39)<23)
    check('cot_product_guard',R(13,4)*23<75)
    check('small_strip_guard',9*(15+11)+75<400)
    check('cot_exponential_guard',5*R(2,64)==R(5,32))
    check('large_imaginary_guard',R(5,2)+R(5,32)<3)
    zero('large_real_quadratic_guard',R(7,2)**2-3*R(7,2)+2-4+R(1,4))
    zero('positive_H_lower',R(39,5)-R(5,2)-R(53,10))
    check('large_real_bracket_guard',4*3==12 and 12<400)
    check('benchmark_excludes_strip',30*2000**10>400)
    s=d+S.I*t
    p=100-(s-R(3,2))**2
    zero('holomorphic_reduction_real_part',S.re(S.expand(p))-(100-(d-R(3,2))**2+t**2))
    zero('holomorphic_reduction_strip_floor',100-2**2-96)

def index_and_extension():
    s=S.symbols('s')
    q=s*(s-3)
    zero('source_root_simple_slope',S.diff(q,s).subs(s,0)+3)
    zero('mass_root_simple_slope',S.diff(q,s).subs(s,3)-3)
    ip,im=S.symbols('ip im')
    values=S.solve([ip-im+4,ip+im],[ip,im])
    check('four_roots_adjoint_index',values=={ip:-2,im:2})
    k,r=S.symbols('k r',integer=True)
    kernel=k+1-r;cokernel=k+2-r
    zero('mass_extension_index',kernel-cokernel+1)
    zero('mass_compatibility_dimension',cokernel-kernel-1)
    check('forced_scalar_correction_on_contour',R(3,2)+2==R(7,2))
    check('post_correction_remainder_beyond_contour',min(R(3,2)+4,4)>R(7,2))

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
            weights_and_kernels()
            strip_guards()
            index_and_extension()
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
