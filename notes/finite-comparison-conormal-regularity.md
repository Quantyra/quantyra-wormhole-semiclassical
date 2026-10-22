# All finite conormal orders for the original comparison family

2026-10-21, E003/S167. Ordinary proof, with supporting checks specified
before execution. Keep the fixed prepared benchmark, positive temperature,
full angular tower, HH state, paired domain and common subtraction.

The fixed seed source and the ORIGINAL finite comparison metric have
conormal outer expansions to every finite order. This supplies higher
metric regularity previously unavailable from the quantitative C5 bound.
The new constants are finite but unevaluated. This is not a quantitative
nonlinear neighborhood or fourth-order auxiliary-solve certificate.

## Statement and coordinates

Use [the seed source](stationary-centered-exterior-source.md),
[its quantitative three-derivative estimate](stationary-centered-conormal-source.md)
and [the original completion](stationary-centered-first-order-geometry.md).
Put z=B/r, and distinguish two outer defining functions:

    Yp=sqrt(C) integral_0^z calB(v)^(-1/2) dv,
    Yc=Ddisk-chi,
    calB(z)=x+z^2-2az^3+Q0 z^4+d z^4(1-z)^2,
    Q0=1+3x, C=1+6x+d, x>0.

Both are analytic positive multiples of z near zero. Yp is the proper
product-frame distance used in the source proof; Yc is the common disk
coordinate used in the metric completion. Do not equate their cubic
coefficients. Let Dp=Yp partial_Yp, Dc=Yc partial_Yc, and L(Y)=1+|log Y|.

Here polyhomogeneous means an expansion in integer powers of Y and finite
polynomials in log Y, with a lower bound on the powers and a remainder
estimate after every finite number of Euler derivatives. The expansion
need not converge. Write O_con(Y^m L) when, for each j, its jth Euler
derivative is bounded by a finite constant times Y^m L.

For the actual fixed seed trace Theta=B^2 tau_hat,mass, with
c=-a sqrt(x)/C^(3/2), the result is

    Theta=c C^2/(5pi Yp)+O_con(L(Yp)).                 (1)

For the ACTUAL reconstructed u,V, with the fixed mass m1 and
c_pole=kappa sinh(Ddisk)/lambda_e, it is

    u=u3 Yc^3+O_con(Yc^4 L(Yc)),
    V=V3 Yc^3+O_con(Yc^4 L(Yc)),
    u3=-m1/(lambda_e c_pole^3),
    V3=m1/(2lambda_e c_pole^3), u3+2V3=0.            (2)

The original family sigma_ell=sigma0+ell u and
r_ell-B=(r0-B)exp(ell V), 0<=ell<=1, has polyhomogeneous
compactified metric coefficients to all finite orders. It remains smooth
at the bolt and in the open disk. For every specified finite derivative
order, its conormal constants can be chosen uniformly in ell in [0,1].
They may depend on the fixed benchmark, temperature and derivative order;
no uniform extremal limit, analytic radius or numerical C7/C11 bound is
claimed. All physical model choices are unchanged.

## The seed boundary problem, without a smooth double

calB(0)=x>0. The inverse function theorem applied to Yp(z) makes the
half-side coefficients analytic near the outer boundary. The finite
temperature bolt is a smooth disk center, not another boundary. The
complete Euclidean orbit is a compact smooth disk with smooth outer
boundary and a fixed thermal spin structure.

By [the actual square construction](stationary-massive-boundary-layer.md),
in the outer radial half-density frame its operator is

    L2=-partial_y^2-f(y)^2 partial_t^2
          -gamma_y gamma_t f'(y) partial_t,
    Pi_D v|_0=0, Pi_N partial_y v|_0=0.              (3)

The complementary projectors are constant in this collar. The principal
symbol is scalar positive definite; f is smooth and bounded away from
zero there. L2 is the original nonnegative complementary-adjoint square.
No self-adjoint first-order Dirac replacement is made. The partner trace
equality and per-exterior normalization remain those of the paired-domain
proof. The even double generally has a cubic absolute-value term and is
only C2. We construct the kernel on the smooth HALF instead.

The source-method reference is [Grieser, Notes on Heat Kernel Asymptotics](https://web.math.ku.dk/grubb/notes/heat.pdf),
Theorem 3.1, Proposition 3.8, Lemma 3.10 and Proposition 3.11,
printed pages 16, 22-24. Its displayed theorem is Dirichlet. The following
constant-projector construction supplies the mixed matrix adaptation.

In a boundary patch, freeze the tangential principal coefficient at the
source point y'. Multiply its tangential Gaussian by

    Pi_D [G_tau(y-y')-G_tau(y+y')]
       +Pi_N [G_tau(y-y')+G_tau(y+y')],
    G_tau(s)=(4pi tau)^(-1/2) exp(-s^2/(4tau)).      (4)

This has the delta initial value and satisfies BOTH left boundary
conditions exactly. Boundary cutoffs can have zero normal derivative.
Patch to interior Gaussian approximations in smooth bundle frames,
using no polar-frame singularity at the bolt. The patched kernel K1 has
error R=(partial_tau+L2)K1 of one half parabolic order better than the
uncancelled principal equation: each coefficient difference supplies
a spatial difference and the matrix connection term has only one
derivative. The leading half-space operator annihilates (4).

The convolution estimates apply entry by entry with finite matrix sums.
They use Gaussian decay and the integrable time powers, not a scalar
maximum principle or commuting projectors. The correction
K1-K1*R+K1*R*R-... retains the mixed conditions because the boundary
operator acts on the LEFT K1 in every convolution. After enough factors
the remainder has any prescribed number of derivatives and powers of
sqrt(tau); the time-simplex factorial bound gives convergence at positive
time. Thus the exact matrix kernel has the same uniform diagonal form
as this local construction. Uniqueness for the nonnegative square
identifies it with the physical heat kernel, not an auxiliary state.

After taking the physical spin trace and density conversion, on a fixed
collar this gives

    k2(tau;Y,Y)=tau^-1 [a(tau,Y)+b(sqrt(tau),Y/sqrt(tau))]. (5)

a is smooth down to tau=0,Y=0. b is smooth in its two variables and
rapidly decreasing in the second, with all derivatives. Cutoffs are
understood. In particular its expansion in sqrt(tau) has rapid profiles
b_j(Y/sqrt(tau)) and differentiated remainder estimates to any finite
order. This statement includes the full matrix connections. Taking a
trace of (4) alone would not determine the higher b_j.

## Common angular sum and proper-time integration

Retain the exact mass weight

    Atilde(tau)=sum_(n>=1) 2(q+2n) M_n^2 exp(-M_n^2 tau),
    M_n^2=n(n+q)/C.

It is minus the time derivative of the sphere Dirac-square heat trace;
the zero modes drop out under that derivative. On the smooth closed
sphere the local Gaussian construction gives the ordinary integer heat
expansion with every finite differentiated remainder. Consequently
Atilde=tau^-2 h(tau), with h smooth down to zero. This is one summed
angular heat kernel, not a set of separately renormalized modes.

Multiply (5) by this weight before taking the finite part. The interior
product is tau^-3 times a smooth function of tau,Y. Its first three
Taylor terms are exactly the common bulk subtraction through a2.
Subtracting those terms leaves a smooth-in-Y integrable remainder;
the fixed-scheme restored local coefficients are also smooth here.
This argument retains their values and scale, including the existing
angular finite constants. It does not redefine a finite counterterm.

Each boundary expansion term reduces, after including h, to

    I_j(Y)=integral_0^tau0 tau^(j/2-3) b_j(Y/sqrt(tau)) d tau
          =2Y^(j-4) integral_(Y/sqrt(tau0))^infinity
                           v^(3-j) b_j(v) dv.       (6)

Split the last integral at a fixed positive v and expand b_j at zero.
A Taylor term v^k integrates to a power unless k=j-4, in which case
it integrates to log Y. Rapid decay controls the other endpoint. This
proves integer powers bounded below by -4 and at most one log per term.
For any desired expansion order and number of Euler derivatives, first
retain sufficiently many heat terms and Taylor coefficients. The same
split bounds the remaining integrals and all those derivatives. This
establishes conormal remainder control, rather than differentiating an
unsigned value bound. Smooth interior contributions supply ordinary
nonnegative powers. A different fixed tau0 changes only a smooth term.

The integral over tau>=tau0 is smooth up to the outer boundary to every
finite order. Indeed heat-semigroup smoothing on the compact smooth
boundary problem bounds derivatives at tau0/2, while the remaining
semigroup is bounded in L2. The angular gap M1>18 supplies integrable
exponential decay, including after summing the angular tower. The actual
thermal period, spin structure and global bolt information remain in
these smooth coefficients. They have not been discarded as small.

Thus Theta is polyhomogeneous. The already proved VALUE estimate
Theta-c C^2/(5pi Yp)=O(L(Yp)) now excludes all remaining negative-power
coefficients. Uniqueness of real-power logarithmic expansions makes this
exclusion valid term by term. The remaining constant and log terms, and
all positive-power terms, obey (1) after any finite Euler derivative.
This is the new conclusion; the old value inequality alone would not
have implied it. For example Y^4 sin(Y^-2) has a small value but an Euler
derivative containing -2Y^2 cos(Y^-2), outside O(Y^4 L).

## Exact source and metric reconstruction

Use t=r/B and Bfun(t)=F0(t)/t^2. The actual HH Ward reconstruction is

    I(t)=integral_(1+h)^t Bfun'(s) Theta(s) ds,
    B^2 p2,mass=I/(2F0),
    B^2 rho2,mass=I/(2F0)-Theta/t^2.                (7)

It fixes the integration constant by regularity; no free thermal datum
is added. Analytic changes of a boundary defining function preserve
polyhomogeneity. Products, derivatives and primitives also preserve it:
integrating z^m(log z)^k is elementary, with a higher log power only at
m=-1. At z=0, F0=z^-2 calB is invertible as a Laurent factor. Therefore
(7) is polyhomogeneous. The LLL anomaly and HH integration constant,
local Weyl terms and fixed reference potential have the same property;
their coefficients are rational seed expressions and finite powers of
log z. The angular pressure follows from the retained trace anomaly.

Perform the previously established massive/Weyl cancellation BEFORE
integrating the physical source. The resulting sphere-integrated rho,p
have their known O(z^2 L) tails, now with every conormal derivative.
No divergent massive and Weyl masses are integrated separately. At the
bolt the common-subtracted heat source is smooth and HH regularity
gives rho+p=0 there. All interior derivatives are finite on compact
subsets of the fixed disk; no uniform zero-temperature estimate follows.

The linear constraints give the actual functions

    a(r)=-G integral_r^infinity (rho+p)/(s F0(s)) ds,
    f(r)=-[2kappa R s_R+2G integral_R^r rho(s) ds]/r,
    H=a+(f+s_R F0')/F0, J=integral_r^infinity H/F0,
    xi=s_R-F0 J, V=xi/(r0-B),
    u=a+(f+s_R F0')/(2F0)-s_R/r-(F0'/2-F0/r)J.      (8)

These are products and primitives of the above expansions. At the bolt
the exact fixed-temperature identity H(R)=0 removes the apparent
division singularity; even smooth radiality follows in the disk frame.
At infinity (8) gives f=-2m1/r+O_con(r^-2 log r) and
a=O_con(r^-4 log r), and hence (2).

For clarity, the analytic coordinate relation at the outer endpoint is

    v(z)=-h-kappa B integral_0^z calB(s)^(-1) ds,
    chi(z)=2atanh(exp(v(z))),
    Yc'(0)=kappa B sinh(Ddisk)/x>0.

This proves the coordinate-transfer assertion without identifying Yp
and Yc. It also fixes c_pole and the cubic normalization in (2).

Finally put Z_ell=B/r_ell. Exactly,

    Z_ell=z/[z+(1-z)exp(ell V)],
    log(r_ell/r0)=log[z+(1-z)exp(ell V)].            (9)

The denominator is positive. Polyhomogeneous functions vanishing to
positive order are closed under exp and log(1+...), as follows by a
finite Taylor expansion with differentiated integral remainder. Applying
this to (9) and sigma_ell proves the family assertion, uniformly for
ell in the compact interval. The existing positivity and bolt estimates
remain in force. No new comparison family has replaced the original.

The exact physical lapse identity remains

    alpha_ell-alpha0=2log(r_ell/r0)+2ell u-ell V
                       -log[1+ell(r0-B)V_r].       (10)

Since (r0-B)V_r=-3V3 Yc^3+O_con(Yc^4 L), its cubic coefficient
is 2ell u3+4ell V3=0 at every ell. Thus the higher-order regularity
proof preserves the correlated mass cubic required by the outer source
cancellation. Terms Yc^4 log Yc remain allowed: conormal regularity
does not mean ordinary C-infinity regularity of the compactification.

## Consequence and limits

The earlier quantitative C5 estimate is unchanged. Its continuation to
higher orders is now known to exist for this particular fixed-seed
comparison family. This closes a qualitative metric-regularity gap;
the constants in (1)-(2) are not estimates adequate for the cancellation
targets. A bounded C5 family in general has no such higher uniform
bounds: N^-5 sin(NY) on an interior interval has bounded derivatives
through five and unbounded sixth-derivative suprema as N grows.

The heat construction above used the SMOOTH SEED. Applying it directly
to the finite metrics, whose compactifications may contain Y^4 log Y,
would require a further argument with nonsmooth coefficients. In
particular this note does not prove higher regularity of the FULL
quantum source evaluated on those finite metrics. That source enters
the actual residual p in [the positive solve](finite-positive-weighted-inverse.md).
Nor does the note identify the stronger D_cert closure or bound
D0 K^-1p, D0=T(T+4). Higher metric coefficients alone do not do so.

The checker tests exact boundary models, Mellin identities and independent
proper-time quadrature, Ward reconstruction, coordinate and completion
identities, and controls rejecting invalid differentiations. It does
not evaluate the actual charged heat kernel, quantify the new constants,
or formally verify the continuum heat-parametrix argument. No novelty
is claimed for the general heat method or polyhomogeneous calculus.

Remaining to-do list: differentiated quantum source on the finite
polyhomogeneous comparison and numerical higher-order bounds; fourth-order
certificate admissibility and weighted transfers; signed compatibility,
cold/higher response, nonlinear preservation/existence; physical state,
boundary/validity, stability, matching, charge/return, opening, transport,
preparation, separate creation and full audit.
