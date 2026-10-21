# Full thermal shape operator and its variational inverse

2026-10-21, E003/S168. Ordinary argument before focused verification.
The full research goal remains incomplete.

The [global thermal coercivity theorem](global-thermal-curvature-free-form.md)
has a closed-form realization in its own energy norm. For every positive
Matsubara frequency gamma above that theorem's unchanged threshold, the
complete curvature-free shape core defines a densely defined symmetric
operator in H=L2((0,infinity),dy), with w=r v. Its Friedrichs extension
A_gamma,F satisfies

    A_gamma,F>=gamma^2/(40G),
    ||A_gamma,F^-1||_(H->H)<=40G/gamma^2.                         (1)

This includes the classical action, fixed Weyl action and complete
canonical Dirac response with the original state, domain, contacts and
angular tower. The proof identifies the missing operator-to-H mapping;
it does not infer closability from positivity of an arbitrary form.
It is an auxiliary variational realization at discrete thermal frequencies.
Its identification with the physically admitted gravitational domain,
continuous causal inversion, curvature coupling and nonlinear transfer
remain separate requirements. No numerical threshold is certified.

## 1. A bounded global reconstruction and its adjoint

Use the [global classical notation](global-classical-curvature-free-form.md):

    K=-partial_y^2+gamma^2+V_K,
    B=-partial_y^2-2a partial_y+gamma^2+V_K+c,
    a=F/r, q=r p=T_gamma w=-(1/3)K^-1 B w.

K has its retained positive limit-point realization, V_K>=0. Set

    K0=-partial_y^2+gamma^2+6/y^2, W=V_K-6/y^2,
    b=a-1/y, j=V_K+c-12/y^2,
    epsilon=||W||_infinity/gamma^2,
    delta=2||b||_infinity/gamma
               +(2||b_y||_infinity+||j||_infinity)/gamma^2.

At gamma>=Gamma_C the preceding theorem gives epsilon,delta<=.001.
K=K0+W has the same operator domain as K0. On its smooth compact core,
integration by parts gives the exact identity

    ||K0 z||^2=||z_yy||^2+12||z_y/y||^2
       +2gamma^2||z_y||^2+12gamma^2||z/y||^2+gamma^4||z||^2.       (2)

Weighted Hardy also gives ||z_y/y||^2>=9||z/y^2||^2/4. Both
extend to D(K0) by graph approximation, using its limit-point core.
In particular

    ||z_y/y||<=||K0 z||/sqrt(12),
    ||z/y^2||<=||K0 z||/sqrt(27).

For z=K^-1 f, positivity and Cauchy--Schwarz give
||z||<=||f||/gamma^2, ||z_y||<=||f||/gamma, and
||K0 z||<=(1+epsilon)||f||. The unweighted formal adjoint is

    B^*=-partial_y^2+2a partial_y+gamma^2+V_K+c+2a_y,
    B^*-K=(2/y)partial_y+4/y^2+2b partial_y+2b_y+j-W.            (3)

The divergence contact 2a_y must be retained. Equations (2)-(3) imply

    ||B^* K^-1||<=1+[7/(3sqrt(3))](1+epsilon)+delta+epsilon
                  <12/5,
    ||T_gamma||=||T_gamma^*||<4/5,
    T_gamma^*=-(1/3)B^*K^-1.                                   (4)

To define the first map without assuming its boundedness, initially
use its compact-input formula. Pair it against arbitrary f in H,
move the self-adjoint K inverse, and integrate B against compact w.
The bound for B^*K^-1 then gives its unique bounded extension.
The same pairing proves the stated adjoint identity. This argument
does not commute B with K or replace the charged reconstruction by K0.

## 2. The complete Hessian maps its compact core into H

Let D=C_c^infinity(0,infinity) in w, with its smooth thermal factor.
For each fixed w in D, v=w/r is compact and p=T_gamma w/r has the
retained global tails. We show that the polarized Hessian a_gamma(u,w)
equals <u,A_gamma,0 w> for some A_gamma,0 w in H and every u in D.
All constants in the following point estimates may depend on the
fixed core vector and frequency. No uniform H-to-H bound for the
unbounded shape operator is being asserted.

### Canonical Dirac row

Write Qcan's scalar output in gbar measure as Hcan v, so its w-pairing
output is

    fcan(y)=(C F/r^3)(Hcan v)(y).                               (5)

Here the actual time-mode norm is divided out, as in the preceding
theorem. We prove Hcan v is bounded over the compactified orbit.
This requires the actual contacts as well as the separated kernel.

For the positive geodesic model the retained point estimate is

    |sum d_M P_M v|<=2e-5 C sup|X^6 v|.

The restored fourth contact is bounded by C sup|X^4 v|. These suprema
are finite: v vanishes on neighborhoods of seam and bolt, and its
support lies in smooth geometry. Along the complete geodesic flow it
is constant zero when the flow crosses the nonsmooth seam. Thus no
sixth derivative of the seam metric is required. Full long-flow
returns are already included in this point estimate.

The [global actual/geodesic comparison](stationary-curved-quantum-tt.md),
sections 2-6, specifies the remaining symmetric kernel, its common
finite-part contacts and bounded gluing row mass. Its local separated
remainder is bounded by KC/d^4 in two dimensions. In a compact smooth
neighborhood of supp(v), subtract the actual Taylor polynomial through
degree two. The remainder is O(d^3), whose product with d^-4 and the
two-dimensional volume d dd is integrable. All fixed finite contacts
there have smooth bounded coefficients and finite jets of v. This is
the actual regulated extension; no new diagonal term is selected from
an off-diagonal estimate. Outside that neighborhood the observation
has positive distance from supp(v), so no subtraction is needed.
The bounded long-distance row mass controls the remainder there too.

The potentially unbounded zero contact, bounded in the retained proof
by a constant times (1+1/Y), multiplies v and hence vanishes near the
outer seam. Its distributional seam contact also vanishes. The LLL
Polyakov term is local and has compact output. The canonical-frequency
core audit already extends these full two-dimensional representations
to the present nonzero thermal sector. Consequently Hcan v is bounded.

At the outer end C F/r^3=O(y); at the horizon it is
O(exp(-2kappa_h y)). More directly,

    integral_0^infinity (C F/r^3)^2dy
       =C^2 integral_(r_h)^infinity F/r^6 dr<infinity.           (6)

Thus fcan belongs to H. No optical mass gap, evaluation of the actual
kernel, or interchange with an infinite-time horizon limit is used.

### Local classical and Weyl rows

Polarize the exact local quadratic density in
[the fixed optical action](local-coupled-frequency-polynomial.md).
After integration by parts in its first argument, denote its outputs
against independent w and q by fw and fq. Then its restricted output is

    floc=fw+T_gamma^* fq.                                      (7)

For fixed core w both fw,fq belong to H. Here is an explicit endpoint
argument, which is necessary because p is generally noncompact.

Near y=0, v vanishes. The exact seed has

    e^(2V)=B^2 F/r^2=B^2 lambda+O(y^2),
    Phi=log(r/B)=-log y+constant+O(y^2),
    p=O(y^4), partial_y^j p=O(y^(4-j)), 0<=j<=4.

The regular Frobenius branch supplies the differentiated p bounds;
the seed's reciprocal radius solves an analytic regular ODE at y=0.
In the pure-p and mixed p,v Weyl rows the optical density has fourth
derivatives with bounded leading coefficient, derivative coefficients
with poles of order at most two, and at most a third-order pole after
one derivative of such a coefficient. The unbounded explicit Phi
coefficient multiplies only the v,v fourth term, whose core input is
zero here. Directly in the retained F_W polynomial, the p,p terms are
bounded multiples of (p_yy-gamma^2p)^2, Phi_y p_y(p_yy-gamma^2p),
Phi_yy(p_y^2+gamma^2p^2), Phi_y^2(p_y^2+gamma^2p^2), and p_y^2+gamma^2p^2.
Their Euler row is O(1). The mixed row uses the same poles and at most
four derivatives of p, and is O(1) too. The classical p,p and mixed
rows have coefficients at worst r^2=O(y^-2) in derivative terms and
F r^2=O(y^-4) in potential terms, again giving O(1) outputs against
p,v. Conversion to q,w divides both rows by r, yielding fw,fq=O(y).

At the horizon use rho=exp(-kappa_h y) times a fixed smooth positive
normalization. The nonzero thermal mode is an integer n-th harmonic:
p=rho^n times a smooth function of rho^2. Background scalar fields
are smooth functions of rho^2 and the canonical metric is regular
in Cartesian bolt coordinates. The exact local Euler rows are scalar
densities of a covariant local action; their scalar factors are smooth
n-th harmonics. The optical volume contributes F=O(rho^2), while
1/r is smooth and bounded. Hence fw,fq=O(rho^(n+2)) there.
They are smooth in the compact middle region, proving H membership.

The integrations against p_u have zero outer flux by its O(y^4)
jets. At the bolt the smooth Cartesian calculation has no boundary;
the shrinking-circle flux vanishes. This uses the actual reconstruction,
not an independent cutoff of p. Equation (4) now proves (7) lies in H.

Combining (5)-(7) gives A_gamma,0:D->H. Its symmetry follows from
the actual Euclidean Hessian, the fixed contacts and the adjoint in
(4); no symmetry of a general complex-frequency retarded form is
assumed. The operator contains every term of the preceding coercivity
theorem, so on D

    <w,A_gamma,0 w> >= Gv/(40G)
                  +C J2[v]/(200pi)+C J2[p+v/3]/(40pi)
                  >=gamma^2||w||^2/(40G).                     (8)

## 3. Closure, inverse and response estimates

For clarity, the closure argument can be written without assuming any
logarithmic Sobolev equivalence. Suppose w_j in D is Cauchy in the
positive a_gamma norm and tends to zero in H. For every phi in D,

    a_gamma(w_j,phi)=<w_j,A_gamma,0 phi> ->0,

because A_gamma,0 phi is in H. The core is dense in the abstract form
completion, so the limit vector is orthogonal to a dense set and is
zero. This proves closability. Its closure is a densely defined closed
positive form, and its associated Friedrichs operator proves (1).

This is the standard construction described by
[Teschl, second edition, Theorems 2.13-2.14](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe2.pdf)
and [Ashbaugh et al., equations (2.14)-(2.21)](https://www.mat.univie.ac.at/~gerald/ftp/articles/KreinWeyl.pdf).
The new work is the full charged-core mapping in section 2, with the
bounded global reconstruction in section 1. No new abstract extension
theorem is claimed.

For f in H let w=A_gamma,F^-1 f, v=w/r, p=T_gamma w/r. Then

    ||w||<=40G||f||/gamma^2,
    ||r p||<=32G||f||/gamma^2,
    a_gamma[w]=<w,f><=40G||f||^2/gamma^2,
    Gv<=1600G^2||f||^2/gamma^2,
    J2[v]<=8000pi G||f||^2/(C gamma^2),
    J2[p+v/3]<=1600pi G||f||^2/(C gamma^2).                     (9)

The preceding proof also retained a_gamma[w]>=N_gamma/(33G), since
1/32-1e-24>1/33. Thus the closure embeds in H1_0 for w. The other
norms in (8) pass by closed weak-derivative maps and the bounded
reconstruction, proving (9) on the full form domain. The stronger
completion from the previous theorem embeds continuously in this
closed form domain, since its norm controls a_gamma. Equality of the
two domains is not needed and is not asserted. On the larger form
domain the compensated logarithmic expression is defined by its
combined core limit, without requiring each logarithmic summand to
be finite separately or changing its finite contacts.

The inverse is also the unique weak solution for a continuous functional
on this form space, by Riesz representation. Showing that an actual
curvature or boundary forcing defines such a functional is additional
work. The inverse has not been substituted into the coupled physical
equations without that check.

## 4. Verification and remaining physical problem

The focused checker verifies (2)-(4), the divergence contact, the
spectral-gap and response arithmetic, and exact endpoint jets of the
local Weyl Euler rows. Independent quadrature tests the normal graph
identity. Refined charged boundary-value solves compare the direct
reconstruction with its unweighted adjoint pairing; omission of 2a_y
is detected. These are identity/domain-supporting fixtures, not an
evaluation of the full quantum inverse or a calibration of Gamma_T.
Continuum kernel and closure arguments remain ordinary mathematics;
no Lean run or formal continuum certification is claimed.

All statements concern the charged seed and its fixed matter domain.
The Friedrichs closure of compact gravitational shape inputs is an
auxiliary realization. Its physical boundary/initial interpretation,
the causal operator between thermal frequencies, the curvature cross
block, independent state response and exact nonlinear background are
still unresolved. Thermal invertibility by itself excludes neither
growing modes of the coupled problem nor its dangerous curvature branch.
S168.3-S168.4 and all physical completion gates remain open.

Remaining to-do list: identify the physical shape realization and control
continuous causal frequencies and curvature coupling; prove nonlinear
existence and signed compatibility; complete exterior/state, charge/return,
opening, all-sector stability, transport, preparation, separate creation
and full audit.
