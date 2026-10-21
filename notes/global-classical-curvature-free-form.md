# Global classical curvature-free shape bound at large real frequency

2026-10-21, S168. The classical Einstein--Maxwell--cosmological part of
the retained action has a positive curvature-free shape form on the
whole charged seed exterior at sufficiently large real Laplace frequency.
The conformal reconstruction and both its endpoint tails are retained.
This is a global classical component result. The full shape form also
contains the fixed higher-curvature/Weyl action and the canonical Dirac
response; its global inverse and physical stability remain unresolved.

## 1. Statement and unchanged model

Use the centered nonextremal seed, corrected outer curvature lambda,
physical cosmological parameter lambda0, fixed physical magnetic length
squared Q_phys, prepared clock and the conventions of the
[curvature reconstruction](physical-curvature-reconstruction.md).
Inward tortoise distance obeys r_y=-F, with y=0 at the outer boundary and
y=infinity at the horizon. Let v be any smooth compact exterior optical
shape and let p be its unique curvature-free conformal reconstruction:

    Psi=delta R-r R_r p=0, q=r p, w=r v.

The subscript C below means exactly the classical action in
[the fixed optical action](local-coupled-frequency-polynomial.md),
with the original lambda0 and physical magnetic charge. It is not the
full semiclassical action and is not an order-reduced replacement for it.
For real gamma>=Gamma_C, where the finite coefficient threshold is
specified in (11), we prove

    S_C,gamma[v] >= (1/(40G)) integral r^2(v_y^2+gamma^2 v^2)dy. (1)

The closure in w variables defines a positive classical operator on
L2(dy), with form domain H1_0(0,infinity), and

    ||S_C,gamma,F^-1||_(L2_w -> L2_w) <=40G/gamma^2.       (2)

These are fixed-real-frequency variational realizations. They do not
select physical gravitational boundary data, establish a causal inverse,
or invert the full retarded shape block. Gamma_C is an explicit expression
in suprema of known seed coefficient functions; those suprema have not
been numerically enclosed for the retained benchmark in this milestone.

## 2. Exact whole-exterior representation, including residuals

Define radial coefficients (r derivatives carry the subscript r)

    H_R=2R+r R_r, D=F_r/r-H_R/6, V_K=F D,
    C_v=2F_rr+4F/r^2+8F_r/r,
    a=F/r, c=F(F_rr+4F_r/r+H_R/6),
    K_gamma=-partial_y^2+gamma^2+V_K,
    B_gamma=-partial_y^2-2a partial_y+gamma^2+V_K+c.

Multiplication by r is unitary from L2(r^2dy) to L2(dy). Directly applying
it to the exact curvature identity gives

    q=-(1/3)K_gamma^-1 B_gamma w,
    E_gamma[w]=<B_gamma w,K_gamma^-1 B_gamma w>,
    <q,K_gamma q>=E_gamma[w]/9.                         (3)

The inverse is the already established positive limit-point realization;
V_K>=3lambda F>=0. It is the auxiliary reconstruction inverse and not
the fast scalar factor. For compact w its output has a y^3 outer tail
and a decaying horizon tail. We keep these tails in every integral.

Let T_C=R+12lambda0 and J=1+3lambda0 r^2-Q_phys/r^2. Then the exact
classical congruence along (p,v), with q as in (3), is

    G S_C,gamma[v]=E_gamma[w]/3-2 integral F J w^2/r^2
       -integral F[(2T_C+r R_r/2)q^2+2T_C q w].         (4)

For complex profiles polarize the real form, using real parts in cross
terms. In particular Q_phys is not silently replaced by the seed's
comparison 1/r^2 coefficient.

Here is a derivation which retains the off-shell contacts. The classical
conformal block and direct shape block are

    Q_pp,C=-(1/(2G)) integral [6r^2(p_y^2+gamma^2p^2)
                             +F r^2(2R+48lambda0)p^2],
    Q_vv,C=-(2/G) integral F J v^2.

The fixed-density conformal gradient identity, evaluated at Psi=0, gives

    Q_pp,C+Q_pv,C=-(1/(2G)) integral F r^2[
                         (r R_r+4T_C)p^2+2T_C p v].

Therefore S_C=Q_vv,C-Q_pp,C+2(Q_pp,C+Q_pv,C). Integration by parts in
the positive reconstruction energy converts this expression into (4).
The global p tails justify the endpoint integrations; no seed equation
T_C=0 is assumed. Even the corrected cosmological shift stays present.

## 3. A bounded perturbation of the normal reconstruction problem

On L2(dy) set

    K0_gamma=-partial_y^2+gamma^2+6/y^2,
    B0_gamma=-partial_y^2-(2/y)partial_y+gamma^2+12/y^2,
    E0_gamma=<B0_gamma w,K0_gamma^-1 B0_gamma w>,
    P=integral w^2/y^2,
    N_gamma=integral [w_y^2+(gamma^2+2/y^2)w^2].

The earlier [all-frequency normal estimate](curvature-free-all-frequency-bound.md),
after replacing its shape profile by y w, is exactly

    N_gamma/9 <= E0_gamma-18P <= N_gamma,
    P<=4N_gamma/9, E0_gamma<=9N_gamma.                  (5)

The factor 18 is essential. The first inequality in (5) is the previously
proved positive normal classical block, not an assumption of positivity
of the charged form. The second follows from the half-line Hardy bound
integral w_y^2>=P/4. All profiles in (5) may have arbitrary compact support
in (0,infinity); no annular restriction is made.

The following coefficient functions are bounded on the entire half-line:

    W=V_K-6/y^2, b=a-1/y, b_y,
    j=V_K+c-12/y^2,
    Pbar=F(1/r^2+3lambda-Q_phys/r^4)-3/y^2,
    Wg=F F_r/r-2/y^2.                                  (6)

To verify this globally, write the exact seed as

    F=lambda r^2+1-2m/r+Q_cmp/r^2+c3/r^3+c4/r^4,
    c3=-2d B^3, c4=d B^4.

Its outer tortoise expansion is
r=1/(lambda y)-y/3+(m lambda/2)y^2+O(y^3). Thus the limits in (6) are
respectively 2lambda, 0, 2lambda/3, 4lambda, 2lambda and 2lambda/3;
b=O(y). These statements use the differentiated rational expansion,
not just convergence of r itself. At the simple horizon all undifferenced
coefficients containing F decay exponentially, while 1/y and 1/y^2
tend to zero. The coefficients are smooth on every intermediate compact
interval. This proves finite suprema, without requiring the geometry to
be uniformly close to AdS throughout the exterior.

Let the supremum norm be over (0,infinity) and put

    W*=||W||, A*=2||b||, B*=2||b_y||+||j||,
    P*=||Pbar||, Wg*=||Wg||, D_gamma=A*+B*/gamma,
    epsilon=W*/gamma^2, delta=D_gamma/gamma.

For the positive normal form and every compact test u, integration by
parts gives

    |<u,(-2b partial_y+j)w>|
       <= [2||b|| ||u_y||+(2||b_y||+||j||)||u||]||w||.

Since K0_gamma>=-partial_y^2+gamma^2, this proves the dual-energy estimate

    ||K0_gamma^-1/2(B_gamma-B0_gamma)w||<=D_gamma||w||.  (7)

This is a first-order coefficient estimate; it does not multiply an
unbounded derivative by a mere L2 resolvent norm. Also
(1-epsilon)K0_gamma<=K_gamma<=(1+epsilon)K0_gamma as closed forms.
Their inverse ordering follows either from quadratic minimization or
from the variational identity
<f,K^-1 f>=sup_u[2 Re<f,u>-<u,Ku>]. Combining (5)-(7), for epsilon<1,

    E_gamma >= E0_gamma-(9epsilon+6delta)N_gamma,
    E_gamma <= (3+delta)^2 N_gamma/(1-epsilon).         (8)

For the lower bound expand the squared dual norm, keep its cross term
bounded below by -2sqrt(E0_gamma)D_gamma||w||, and use
gamma||w||<=sqrt(N_gamma). This proves (8) for arbitrary profiles,
including any sequence escaping toward either endpoint.

## 4. Small residuals, finite threshold and coercivity

The retained seed guards give x=lambda B^2>1, |d|<1e-12 and
g=G lambda<(4e-17)^2. They are pinned in the
[centered-source estimates](stationary-centered-exterior-source.md).
Set Delta=lambda0-lambda=11G lambda^2/(360pi)>0. Direct curvature gives

    R+12lambda=(d/B^2)(4z^5-6z^6),
    r R_r=(d/B^2)(-20z^5+36z^6), z=B/r in (0,1).

Consequently the three small quantities

    Delta/lambda,
    ||T_C||/lambda <=12Delta/lambda+10|d|/x,
    ||2T_C+r R_r/2||/lambda <=24Delta/lambda+48|d|/x

are all below sigma=1e-9. The deliberately loose sigma already suffices;
the residuals have not been set to zero.

Choose the explicit finite threshold

    Gamma_C^2=max(lambda,1000W*,(2000A*)^2,2000B*,
                          6000P*,1000Wg*).             (11)

At gamma>=Gamma_C this ensures epsilon<=.001, delta<=.001,
6P*/gamma^2<=.001, Wg*/gamma^2<=.001. Equations (5),(8) imply

    E_gamma-18P >=(1/9-.015)N_gamma,
    E_gamma<10N_gamma.

The bound V_K>=3lambda F and (3) give

    integral F q^2<=10N_gamma/(27lambda),
    integral F w^2<=N_gamma/lambda.                    (12)

For the latter, 3lambda integral Fw^2<=6P+W*||w||^2
<= (8/3+.001)N_gamma<3N_gamma. Split the potential in (4) into
3/y^2+Pbar+3Delta F. Its bounded difference costs at most
.001 N_gamma/3 in (4). The three residual contributions cost at most

    sigma[10/27+2+6]N_gamma,

using Cauchy--Schwarz for the qw term and sqrt(10/27)<1. Altogether

    G S_C,gamma >=[(1/9-.016)/3-sigma(10/27+8)]N_gamma
                  >N_gamma/32.                        (13)

Finally, direct integration by parts gives

    integral r^2(v_y^2+gamma^2v^2)
       =N_gamma+integral Wg w^2 <=1.001 N_gamma.        (14)

The right side of (13) is larger than 1/40 of (14), proving (1).
The same estimates give an upper form bound by a finite multiple of
N_gamma/G. Hardy's inequality makes N_gamma equivalent to the H1_0
norm for each fixed gamma>0. The bilinear form therefore extends to
a bounded coercive closed form there. Its associated positive operator
has spectral lower bound gamma^2/(40G), proving (2). This is a global
form statement for the classical restriction in the w pairing.

## 5. Comparison, verification and remaining operator

[Teschl, second edition, sections 2.3 and 6.5](https://radon.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe2.pdf)
supplies the standard closed-form and relative-form framework, including
Theorem 6.24. The inverse ordering and coefficient estimates used here
were written out above. No novelty is claimed for those abstract results.
The earlier normal form is the strongest applicable retained comparator;
the added result is the whole-exterior charged classical transfer with
its residual terms and separate physical charge retained.

The frozen checker verifies the curvature, unitary operators, asymptotic
coefficient limits, residual budgets and coercivity arithmetic. Independent
boundary-value fixtures evaluate (4) against a direct second variation of
the optical Einstein--Maxwell action and check the reconstruction energy
in two ways. The fixtures vary shape support, frequency, charge and seed
parameters, retain both conformal tails, and refine the left endpoint,
collocation tolerance and quadrature. They do not numerically enclose the
benchmark coefficient suprema or establish a physical stability spectrum.
Some identity fixtures deliberately use larger residuals than the theorem
allows, to detect omission of the residual contacts.

Two pre-acceptance exploratory runs had unsuccessful 1e-10 collocation
solves in their first two fixtures; changing to factored profile evaluation
did not resolve that solver outcome. Both exact development sources and
their observed solver statuses are retained in the source review. Their
original stdout was returned by the tool but was not saved as a standalone
raw file; no acceptance claim uses those unsuccessful solves. Before the
acceptance run, freeze successful solves at tolerances 1e-6 and 1e-8,
left endpoints 1e-4 and 5e-5, and 48/96 quadrature nodes. Require normalized
action/energy identity errors and action refinement below 2e-7. Three
noncommuting matrix fixtures independently test the inverse-form comparison.
All acceptance runs, exact failed sources and raw logs are retained.

The complete shape congruence remains
S_full=S_C+S_Weyl+Q_canonical. A positive S_C is not an inverse for this
sum. The higher-curvature part has higher differential order, so it cannot
be treated as a bounded error in the H1 norm of (1) without further proof.
The existing full logarithmic annular bound remains useful, but its global
extension and the curvature/shape cross block are still required. The
fast scalar factor, physical boundary/initial data, independent state
sources, exact nonlinear background and all physical gates remain open.
No counterterm, magnetic flux, preparation or ultraviolet prescription
is changed. No Lean campaign is selected.

Remaining to-do list: control the full global shape operator and its
curvature coupling, physical domain and charged stability; establish the
nonlinear solve and compatibility; complete exterior/state, charge/return,
opening, all-sector stability, transport, preparation, creation and audit.
