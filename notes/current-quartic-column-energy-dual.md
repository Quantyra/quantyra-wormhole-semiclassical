# A complete massive-response column in the original energy dual

2026-10-27, E003/S169. Ordinary argument, followed by interval arithmetic.
The previous reply drafted a goal and supplied no scientific evidence.
The goal remains the FULL residual T_t-(I+T_t)F_n, including application
errors. This calculation addresses a necessary quantum column of that
residual. It does not replace the goal by a one-column calculation.
Preserve the closed audit, accepted results, P32 and original norms.

## 1. Actual trial, spaces and contribution to the full residual

Use the fixed paired first-order realization, complementary adjoint and
both squares on the actual finite seed. Keep the physical half, mass cubic
and finite quartic logarithms, fixed clock/period/spin/scale, original
cutoff transitions and smooth disk bolt. No smooth reflected domain is
substituted. The reflected cubic coefficient remains C2,1.

Let E_t be the accepted retained TT form completion, with
e_t=ell_t/N, N=|W_action| beta. T_t=ell_t^-1 B, B=b* a^-1 b, and
0<=T_t<=q_t=.338/(1-epsilon). Both original rows are retained in this
Schur reduction. The canonical York maps are

    -(Q^2 J'/z^4)'=Q(4xi/z-xi')/z^4,
    J(0)=0, (Q^2 J'/z^4)(L)=0,
    Uxi=xi/2+QJ/z, Vxi=-xi-(QJ)'/2, L=z_h.

Choose the actual trial psi(z)=z^4. Its York right side vanishes exactly;
the boundary flux condition and J(0)=0 give J=0, Upsi=psi/2, Vpsi=-psi.
Thus the quantum direction v=-z^4 is global, has its first three outer
jets zero and is smooth in the regular bolt disk. There is no added cutoff,
wall or free temperature mode. Other lifts and their cutoffs are unchanged.

For a numerical trial F_n u=sum_j psi_j a_j(u), the COMPLETE residual is

    r_n(u)[w]=B[u,w]/N-(ell_t+B)[F_n u,w]/N,
    ||R_n||=sup_(||u||E_t=1) ||r_n(u)||E_t*.

A column application error d_j in E_t* contributes at most
sum_j ||a_j||E_t* d_j, in addition to source, scalar-solve, coefficient
and arithmetic errors and the exact full residual. The coefficient maps
and the infinite-dimensional residual are NOT determined by this note.
In particular a zero finite Galerkin residual does not bound R_n.

For sensitivity only, the exact one-dimensional Ritz coefficient obeys
||a_psi||<=q_t/||psi||E_t. Since ||Vw||2<=CLt ||w||E_t and L>.99,

    ||psi||E_t >= .99^(9/2)/(3 CLt),
    ||a_psi|| <= alpha1=3 q_t CLt/.99^(9/2).                 (1)

This is a proved envelope for that specified exact one-dimensional map,
not a bound on coefficients in an unassembled multicolumn F_n. Its
denominator, source functional and numerical solve have not been evaluated.

## 2. First missing implication: uniform regulator removal for this column

Put h=Vw. The accepted comparisons give, on the actual energy completion,

    ||h||2<=CLt ||w||E_t, CLt=CL/sqrt(1-epsilon),
    ||h||infinity<=CIt ||w||E_t, CIt=CI/sqrt(1-epsilon).

Here CL<22.474 and CI<7.477e50. These are established seed constants,
not neighborhood bounds. Core approximants converge in both norms. On
the original mass/log core h=h'=h''=0 at the seam; no second-derivative
trace is asserted for arbitrary completion elements.

For fixed smooth v=-z^4 the exact paired heat identity is

    delta_v delta_h k(t)=-2t partial_t delta_v tr_time[h exp(-tL)].

The accepted intertwiner and fixed-Hilbert-space vertex retain both
square partners and the coincidence/second vertex. This is not a raw
mass-only determinant. The seed and the one-parameter smooth v family
satisfy the actual mixed boundary conditions. The extension below uses
continuity in h, not a smooth heat theorem on rough base metrics.

The disk propagation proof supplies a stronger statement than its old
selected smearing: an absolute-diagonal remainder bound before h is
inserted. The complete paired boundary part has L1(dz) norm <=149 t^3.
Its collar polynomial P(V,t), retained in the accepted boundary receipt,
has nonnegative coefficients and monomials V^p t^q of weight p+2q>=7.
On crossing paths Y<=sqrt(t)S, V<=2sqrt(t)S, and E S^p<=2 Gamma(1+p/2).
The fixed-z diagonal density is f(Y)/(pi t) times the differentiated path
difference, with f<=1.06 on the collar. Consequently its collar L-infinity
constant, with power t^(5/2), is at most

    Cb,in=(1.06/pi) sum P_pq 2^p [2 Gamma(1+p/2)]
                                      t0^(q+(p-7)/2), t0=.001.       (2)

This uses the whole differentiated path function, not a heat coefficient.
For exited crossing paths S>1/(4sqrt(t)). The accepted global cap guards
give, per fixed-z density, the bound

    66 exp(t)/(pi t) E[exp(2sqrt(t)S)(1+sqrt(t)S)^2 1_exit].

On this event 2sqrt(t)S<=S^2/8. Also
(1+sqrt(t)S)^2<=2(1+S^2)<=16 exp(S^2/8). Integrating the full bridge tail
P(S>s)<=2exp(-s^2) gives

    Cb,exit <= (2816/pi) exp(t0) t0^(-7/2) exp[-3/(64t0)].   (3)

The right side bounds the ratio to t^(5/2) for every 0<t<=t0; its
logarithmic derivative is positive there. No event indicator is varied.
All nonnearest thermal images use the accepted enlarged-cap Gaussian
argument. Adding the extra t^(-7/2) weight preserves its endpoint maximum:
4 t0^(-7/2) times its full image guard is <1e-100. Keep this debit.

Set Cb=Cb,in+Cb,exit+1e-100. Interpolating the proved L1 and L-infinity
bounds gives a full boundary L2 bound sqrt(149 Cb) t^(11/4). This retains
the rough reflection, all cap excursions and every thermal winding.
The interior, actual-disk and smooth-cap propagation remainders have
POINTWISE bounds in the accepted proof. Their sum, including the
product-gradient contact pointwise allowance 1/(120pi), is <953 t^3.
We use the new cap split directly; there is no discarded old-to-new return.

Let Delta=(Q partial_z)', R=-Q'', H=Delta+R. Define the exact polynomials

    O6[v]=-[H Delta Hv+(Delta R)Hv+Delta(v Delta R)]/105
      -[16R^3v+24R^2 Delta v+24Delta(R^2v)+24Delta(R Delta v)]/1008,
    O6prod[v]=O6[v]+(4/63)v,
    A6=||O6prod[v]||L2(0,1)/(4pi), Cb2=sqrt(149 Cb).

The actual polynomial Q, including finite tau, is used. Extending only
the positive polynomial-square integral from L<1 to 1 is an upper bound;
it does not extend or change the physical operator. All core integrations
by parts have zero outer flux by the mass/log jets and zero bolt flux by
regular-disk matching. Uniform convergence in the proved energy
comparisons then extends the resulting identities without inventing
additional traces. Thus the complete product-centered heat remainder
obeys, for EVERY w in E_t,

    |rprod_hv(t)| <= CLt [A6 t^2+Cb2 t^(11/4)+953 t^3] ||w||E_t. (4)

For the full angular sum Theta(t)=sum_(n>=1)2(q+2n)
exp[-n(n+q)t/C], Theta<=4C/t. In the ORIGINAL normalized energy units,
the massive response is -(G/(B^2 beta)) integral Theta rprod dt/t.
G=1 here. The complete omitted common heat cap at eta<=t0 therefore has

    d_cap(eta) <= 4C CLt/(B^2 beta)
          [A6 eta+(4/7) Cb2 eta^(7/4)+(953/2) eta^2].        (5)

This proves an operator-norm regulator limit E_t->R for this fixed
column. It is not a two-variable rough metric variation or a nonlinear
neighborhood theorem. The exact ultraviolet logarithm and finite local
restorations remain in their common allocations; none is set to zero.

## 3. Complete positive-time estimate without the global supremum loss

The accepted boundary-compatible parametrix has pointwise kernel bound
A/t; its complete Volterra correction is bounded by (E0+Ehigher)t^3.
Use A=max(2 Acap,Ainterior+cutoff*t0^3), retaining the full image guard.
At t*=t0/2, rank two gives

    M=2[A/t*+(E0+Ehigher)t*^3]<734.                         (6)

All quantities come from the accepted actual-disk propagation receipt.
The heat diagonal tr K(t;x,x) decreases for t>0 by its full positive
spectral representation. Hence M bounds it for ALL t>=t*. It is the
kernel in the original Riemannian volume, not a cylindrical density
carried through the bolt. Time reduction leaves precisely dz.

For the paired self-adjoint Dirac A with A^2=L, the exact heat Hessian
coefficient is q_t(a,b)=-t(a+b)^2[g_t(a^2)-g_t(b^2)]/(a^2-b^2),
g_t(x)=x exp(-tx). The accepted scalar proof gives
|q_t(a,b)|<=32(exp(-ta^2/2)+exp(-tb^2/2)). Weighted Cauchy-Schwarz in
the FULL double spectral sum, before time reduction, now gives

    |D2 k(t)[h,v]| <=64 sqrt(tr_time[h^2 exp(-tL/2)])
                          sqrt(tr_time[v^2 exp(-tL/2)])
                    <=64 M ||h||2 ||v||2, t>=t0.            (7)

The factor 64 includes both weighted sums. All angular modes are in
Theta; all thermal/radial modes and finite returns are in this heat
kernel. No finite set of modes defines (7).

Move the derivatives of the completed a2 contact off h on the original
core, retaining the regular-bolt and outer boundary terms (which vanish):

    P4[v]=4R^2v+8R Delta v+8Delta(Rv)+8Delta^2v.

Use ||v||2<=1/3, C0=2 CLt/(3pi), C2=CLt||P4[v]||L2(0,1)/(480pi).
This is an exact identity, not a replacement of the contact by its
principal symbol. It avoids the unnecessary C_D norm of Delta h.

The full product-gradient term remaining after its a0/a2 subtractions is

    cprod(t)v=(8t/pi) integral_0^infinity
       lambda^3(1-exp(-t lambda^2))/(exp(2pi lambda)-1) d lambda *v,
    0<=cprod(t)<=t^2/(63pi).

This preserves the contact's sign and complete spectral definition.
Put Cp=CLt/(189pi) and CH=64 M CLt/3. Since M1^2>324 and
Theta(t)<=4C/t0 exp[-324(t-t0)], the ENTIRE positive-time contribution has

    d_positive <= (4C/(B^2 beta t0*324))
        [CH/t0+C0/t0^2+C2+Cp(t0+1/324)].                    (8)

Thus d_mass=d_cap(t0)+d_positive bounds the complete product-centered
massive remainder column in E_t*, including its common regulator limit.
The separately retained finite a0/a2, Weyl, LLL and classical terms are
NOT omitted from ell_t and are not evaluated by calling (8) complete.
They must still be included in an assembled residual. A center zero with
radius d_mass is a rigorous unsampled envelope for this particular
remainder; it is not an evaluated signed quantum center.

### Match the original canonical assembly and retained operator

The original canonical assembly writes the uncentered heat remainder r,
with the fixed vacuum/local restoration separate. It is related exactly
by rprod=r+cprod(t) int h v. Do not mix these two assemblies. For r use
A6can=||O6[v]||2/(4pi) in (5), and omit Cp from (8). The bulk allowance
953 still dominates the actual propagation remainder without its product
contact. Thus d_can bounds the ORIGINAL canonical massive column; the
product-centered certificate above is retained as a separate equivalent
decomposition with its explicitly accounted local completion.

The accepted known TT matrix has its third trial equal to this psi and
contains the Einstein, magnetic, R2, Weyl, local a0/a2 and LLL terms.
Let hknown be that diagonal entry in I/W units. W is NEGATIVE, so its
contribution in the present energy norm is c=-hknown/beta, approximately
2.26666663. This sign and beta conversion must be retained.

Write x=||psi||E_t and e=epsilon/(1-epsilon). The complete canonical
column bound and the accepted full transport-tail comparison imply

    (1-e)x^2-d_can x-c_upper<=0,
    (1+e)x^2+d_can x-c_lower>=0.                            (9)

The positive roots give proved x_lower,x_upper. This assembles a full
diagonal norm enclosure, not a signed full column. In particular the
retained-operator column transfer costs at most e*x_upper, and the exact
one-dimensional Ritz coefficient has the sharper proved envelope
alpha_exact=q_t/x_lower. This is an actual numerical comparison in the
original energy norm, not a coordinate change. The complete numerical
F_n still has no certified source or coefficient application.

For the one-dimensional sensitivity test, use

    d_cut=alpha_exact*(d_cap,can(t0)+e*x_upper),
    d_unsampled=alpha_exact*(d_can+e*x_upper).                (10)

The first is an error allowance for omitted short time and the full
retained-operator transfer, not for the signed positive-time application.
The second treats that entire unknown signed remainder as a zero-centered
error interval. If it exceeds 1e-30, the missing evaluation must reduce
the ENCLOSURE ERROR; the actual response itself need not be smaller.
All other residual/source/solve/application contributions remain unknown.

To display the original physical sensitivity, the accepted bolt error
function is Gamma(z)=W0(3z+3z^2+z^3)+W1(2z+z^2), with the original radius
conversion inside W0,W1 and z=d/f_t. Evaluate this for (10) only as a
best-case diagnostic with every other residual set to zero. The entire
bolt allowance is 1.667e-9; it is NOT an available remaining margin.
Nonlinear cross terms prevent assigning these values independently when
other residual contributions are unknown. The other original physical,
graph, defect, mass-column and boundary-phase budgets remain open.

## 4. Numerical test and decision boundary

Execute (1)-(8) with outward Arb arithmetic, actual parameter intervals,
the accepted full boundary polynomial and disk kernel constants. Compare
alpha1*d_cap and alpha1*d_mass with 1e-30 only as sensitivity to the
previous sufficient uniform-application target. No positive remaining
margin is assumed, and no multicolumn coefficient is inferred from (1).
Report the improvement over the specified old supremum-based estimate
and any remaining excessive upper allowance; neither is a physical-error
lower bound or a construction obstruction.

An actual F_n still requires its explicit coefficients, complete known
local/scalar/source assembly, the signed retained quantum columns, all
other phase/fast columns, and the infinite-dimensional residual. Its
source, numerical solve and application errors remain unknown. Bounds
(4)-(8) do not control arbitrary pairs or independent boundary variations
and do not transfer the seed estimates to a nonlinear neighborhood.

The method decision follows the executed constants. If the cap is
affordable but the unsampled positive-time envelope is excessive, retain
the paired geometric representation and evaluate that signed column in
the L2-dual topology, rather than refining P32 or the already controlled
cutoff. The next calculation is a certified signed L2 enclosure of the
whole positive-time density for v=-z^4, including its tails and contacts,
as a column of the complete Ritz residual. It must cover ALL test fields;
a few scalar pairings do not give that enclosure.

Remaining to-do list: full residual and actual application errors; signed
G_t and every original coupled/admission and mass/phase budget; nonlinear
existence/contraction and exact compatibility; physical matching,
semiclassical validity/noise, all-sector stability, transport, preparation,
separate creation and the full stable semiclassical wormhole objective.
