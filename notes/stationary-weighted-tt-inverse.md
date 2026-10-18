# The actual weighted TT inverse and remaining scalar equation

2026-10-17 HST, E003/S139 supporting S074/S075. Ordinary argument before
focused verification. **PROGRESS; coupled compatibility remains open.**

At the analytic centered seed the actual TT block is invertible both
at physical L2 weight and on the fast weighted domain with its one
physical mass column. This is an inverse of the complete retained
block, including its massive response and off-shell disk contacts.
It follows from the [curved TT bound](stationary-curved-quantum-tt.md),
the seed York map and the actual weighted parametrix.

The coupled equations therefore have an exact scalar Schur reduction.
Its local fourth-order part and mixed block are given explicitly below.
The reduction fixes the TT response but does not exclude a coupled
kernel or evaluate the physical scalar phase and source moment. The
scalar equation is not positive, and its nonlocal correction cannot
be declared small from the TT L2 gap and small G lambda alone.

## 1. Fixed Hessian, measure and seed coordinates

Keep the analytic [centered seed](stationary-centered-covariant-seed.md),
prepared positive temperature, laboratory clock, flux, paired Dirac
realization and common subtraction. Let A be the complete fixed-volume
disk Hessian of the [weighted map](stationary-global-weighted-map.md):

    I_disk''(u,v)=integral v . A u dmu,
    dmu=4pi beta r^2 dr=gamma(z) z^-3 dz/z,
    gamma(z)=gamma0[1+O(z^2)], gamma0>0.

The outer variable z is proportional to 1/r. All output transposes
below use this fixed measure. Differentiating a gradient per moving
volume would omit its density contact. The seed is stationary for I0,
not for the full retained action.

Use the existing mixed disk domain X_delta. In the end it requires
theta=S+T/3 in z^delta H_b^4 and shape T in z^delta H_b^(4,log),
with measure dz/z and logarithmic multiplier log(2+|D_logz|).
At the bolt use the smooth rotational two-dimensional disk norms.
For the new scalar coordinates put

    Z_delta=z^delta H_b^4,
    W_delta=z^delta H_b^(4,log),
    Y_delta=z^delta L_b^2,

with the same interior/disk patching. At delta=3/2, Y_delta is the
physical scalar L2 space, up to the smooth positive factor gamma.
The complementary weight for the fixed pairing is 3-delta.

Let c be the sphere component of the genuine TT tensor
h_c=diag(-b-2c,b,c,c), with b RADIAL. The exact reconstruction is

    U=r^2 F, M=r^2 F_r-2rF,
    b(r)=-U(r)^-1 integral_(r_h)^r M(y)c(y)dy,
    eta=-a integral_rho^infinity (b+c)/a d rho,
    S_c=c/2-R eta,
    T_c=-(b+3c)/2-(A_orb-R)eta,                      (1)

where a=sqrt(F)/kappa, R=r'/r, A_orb=a'/a, and these last primes
are physical proper-radius derivatives. The actual bolt is
r_h=B(1+h). Its homogeneous radial b~r^-4 tail is fixed by the
integral; it is not an extra adjustable datum.

Define the linear seed reconstruction

    Rcal(tau,c)=(tau+S_c,T_c).                       (2)

It is inverse to the [global York extraction](stationary-global-tt-reduction.md):
solve P Y=T'+4RT and set
tau=S+T/2-[Y'+(A_orb+2R)Y]/4, with c extracted from h_TT.
The two vectors satisfy eta=-Y. Any discrepancy would preserve the
disk gauge and solve eta'=A_orb eta, hence eta=constant*a~z^-1;
the selected outer weight delta>-1 excludes it. Thus these maps are
inverses, rather than only two bounded constructions.

On the analytic seed, for every -1<delta<4,

    Rcal: Z_delta x W_delta -> X_delta              (3)

is a bounded isomorphism. Here are the regularity details required
for using this statement sharply. The vector inverse P^-1 has roots
-1 and 4 and the proved weighted Green bounds in that strip. Its
leading derivative is Y'=2T/3 plus an operator of order minus one,
so tau=theta+order(-1)T. Conversely T_c=-3c/2+order(-1)c and
theta_c=order(-1)c in the interior principal calculus. The seed
coefficients are smooth with bounded conormal jets on rescaled end
annuli. Dyadic frequency blocks of its order-zero maps have rapid
off-diagonal decay. The extra logarithmic ratio (1+j)/(1+k) has
only polynomial growth in |j-k| and is absorbed by this decay.
This proves the logarithmic mapping without claiming it for a merely
C5 nonlinear candidate.

At the bolt the TT integral has leading expression
-2rho^-2 integral_0^rho s c(s)ds, a regular radial averaging operator.
Equivalently use the smooth disk York elliptic inverse and its inverse
return. Both preserve the required radial Sobolev/log norms and
smooth tensor parity. The global z^4 projection tail is allowed
because delta<4. The same maps are bounded at order zero; duality at
weight 3-delta therefore gives the bounded physical output transpose
Rcal* on the same strip.

Use the SCALAR c pairing, not the full TT tensor norm, to define

    Hcal=Rcal* A Rcal = [[L,Bmix],[Bmix*,D]],
    (f_tau,f_c)=Rcal* f.                            (4)

In particular <c,Dc>=I_disk''[(S_c,T_c),(S_c,T_c)]. The preceding
proof's disk estimate gives

    <c,Dc> >= .995 E0[c] >= .99 E0[c]               (5)

on its smooth compact-outer TT core. Here
E0=beta/(2G) integral [3r^2F(c_r+2c/r)^2/2+lambda r^2c^2]dr.
The extra affine-TT return estimate is not needed to deduce (5).
No full-action gauge-annihilation has been assumed.

## 2. Scalar TT parametrix, including the bolt and outer tails

The normal TT reconstruction and the accepted disk change of variables
obey

    Rcal_c(s)=((s+3)/(2(s+1)),-3/2)^T,
    C(s)=[[1,(s+3)/(3(s+1))],[0,2]],
    C(s) Rcal_c(s)=(0,-3)^T.

Using the reflected transpose at 3-s, the actual normal blocks are

    N_D(s)=9 lambda^2 N_T(s),
    N_L(s)=lambda^2 J_c(s), N_Bmix(s)=0.            (6)

The scalar c normalization in (6) is essential. An orthogonally
projected tensor-L2 operator has a different Gram factor.

D has its own elliptic fourth-log principal symbol. In a smooth
Cartesian orbit plane with frozen nonzero covector xi, the TT
constraint and trace give the orbit tensor
-2c(delta_ab-xi_a xi_b/|xi|^2), with sphere components c,c.
The disk return adds the longitudinal orbit component
-2c xi_a xi_b/|xi|^2. Thus at principal level
S_c=c/2, T_c=-3c/2, theta_c=0. The accepted shape coefficient
p^4 log|p|/(120pi^2) becomes

    (3/(160pi^2)) p^4 log|p|                       (7)

on c, in agreement with 9d0, d0=1/(480pi^2). This rotation-covariant
calculation applies through the smooth bolt. It does not infer a
diagonal inverse from invertibility of the coupled matrix. Section 5
also proves directly that the mixed block has at most two local
derivatives, so it cannot supply a missing fourth-log symbol.

The [actual kernel estimates](stationary-global-weighted-map.md)
transfer through the seed York/return maps. On comparable end annuli
their normalized coefficient error is O(z^2 log^N z). The return
Volterra kernels have poles only at the strip edges -1 and 4, and
their actual-minus-normal coefficients have the same two-power gain.
Composition preserves the integrable separated log-coordinate tails
on closed substrips. Separated compact-to-end terms are
O(z^4 log^N z); this includes the TT integral's radial z^4 term.

Use the scalar interior/disk parametrix for (7), and the scalar normal
inverse of (6) in a sufficiently short end. The two-power error there
is small in operator norm and is inverted by a Neumann series.
Cutoff commutators improve differential order and have compact support
or the just-described decaying separated tails, hence are compact.
Full-order errors on a fixed compact annulus are not called compact.
This constructs left and right parametrices for

    D: W_delta -> Y_delta                          (8)

whenever N_T is invertible on that weight. The scalar construction
also upgrades a distributional cokernel null field to the same
mixed-order domain at the complementary weight: use the local
elliptic estimate through four derivatives and the bounded normal
inverse, with the retained finite contact coefficients. This is the
adjoint regularity needed below. The coupled map itself is not
Fredholm at its scalar-phase critical line, and supplies no shortcut
to the scalar TT assertion (8).

## 3. Critical inverse by graph-domain density

On s=3/2+it set R_t=t^2+9/4. The proved normal formula gives

    N_T=A_T R_t+d0 R_t(R_t-2)H_M>0, H_M>14/5.

The [strip root count](stationary-global-weighted-map.md) has only
the simple TT roots 0 and 3 between -1/2 and 7/2. Thus (8) is
Fredholm at delta=3/2. Formal symmetry and the distributional
adjoint regularity just proved identify its kernel and cokernel.
Its critical index is therefore zero.

One more domain argument is needed before applying (5). A general
critical field need not belong to the stronger norm N of the curved
TT proof: its separate positive model and absolute Weyl forms can
diverge, although their combined actual Hessian is defined. Instead
use density in the actual weighted operator domain.

First truncate and smooth c in W_(3/2), preserving smooth bolt
behavior. Each resulting compact sphere component has a radial TT
b with an allowed z^4 tail. For each fixed such c, add a distant bump
a_L psi(r/L), supported in [2L,3L], chosen so

    integral_(r_h)^infinity M[c+a_L psi(r/L)]dr=0.

Since M=-2r+O(1), its coefficient is a_L=O(L^-2). This makes the
entire TT tensor compact at infinity. In log radius the weighted
profile z^-3/2 a_L psi(r/L) is a fixed smooth translate of amplitude
O(L^-1/2); all four b derivatives and the logarithmic multiplier
have the same bound. The analytic coordinate correction has uniformly
bounded decaying jets. Hence the bump tends to zero in W_(3/2).
A diagonal approximation proves exact compact-TT core density there.
It would be incorrect to claim the same density at weight 7/2.

D is already bounded W_(3/2)->Y_(3/2), independently of coercivity.
The quadratic pairing <c,Dc> is therefore continuous in this domain
norm. E0 is continuous as well, since it is a weighted H1 form in
the end and an ordinary regular disk H1 form at the bolt. Thus (5)
extends to ALL of W_(3/2). This is a new graph-domain extension of
the combined Hessian, not an identification with the earlier N domain.

The critical kernel is consequently zero. Index zero gives

    G_TT=D_crit^-1: Y_(3/2) -> W_(3/2).             (9)

The elliptic estimate makes this inverse bounded, and the fixed
pairing makes it symmetric. No arbitrary semibounded quadratic form
has been promoted to a Friedrichs operator. For c=G_TT f_c,

    sqrt(E0[c]) <= (1/.99)||f_c||_(E0-dual),
    <c,Dc> >= .99 lambda/(8piG) ||c||_L2(dmu)^2.    (10)

The full tensor norm is at most 22 times the scalar c norm squared,
so the corresponding tensor gap is .99lambda/(176piG). These are
energy and L2 bounds. A numerical fast H4log inverse constant has
not been established.

## 4. The fast TT inverse includes exactly one mass datum

For f_c in Y_(7/2), let c=G_TT f_c. The two-power actual-minus-normal
remainder improves its weight from 3/2 to 5/2, crossing no TT root.
A second step from 5/2 to 7/2 crosses the simple root 3 and extracts
its residue. Each step is strictly smaller than two, leaving room for
the finite logarithms. The finite coefficient/contact expansion in
the global weighted proof then gives

    c=m_c e_c+v, e_c=chi(z)z^3, v in W_(7/2).       (11)

For the mass term the local correction has degree at least five;
the far tail is O(z^4 log^N z), beyond the final weight. There is no
forced term exactly at 7/2 in this TT bootstrap. The separate radial
z^4 coefficient remains in the fast remainder.

The mass reconstruction is
Rcal_c e_c=(3/4)e_m+X_(7/2), where
e_m=chi z^3(1,-2) is the existing physical disk mass column.
Its bounded image has already been established, and its TT output
under the physical transpose is in Y_(7/2). Thus

    D_m: W_(7/2) plus span{e_c} -> Y_(7/2)          (12)

is onto by (9)-(11), and one-to-one by critical uniqueness. Its
inverse is bounded by the bounded inverse theorem and agrees with
G_TT on fast sources. In particular the TT fast-only map has index
-1 and the TT mass column is transverse. This says nothing automatic
about mass transversality for the full coupled problem.

As a normalization check, its source/mass boundary constant is

    -9gamma0 lambda^2 N_T'(3)
      =27gamma0 lambda^2[A_T-2d0(1-gamma_E+kappa_M)]>0.

It is 9/8 times the earlier disk constant because the scalar c source
and mass vectors are respectively (3/2)(1,-1) and (3/4)(1,-2).
There is still only one correlated physical mass parameter.

## 5. Exact scalar and mixed blocks, with the density contact

The complete local trace, including the fixed magnetic two-form, is
the [retained trace](stationary-weighted-conformal-range.md)

    E=R4+12lambda0+epsilon[Box R4-11R4^2/72
             +7C4^2/12+11J^2/6-20M_em],
    epsilon=G/(60pi), J=Ric-(R4/4)g, M_em=F_em^2,
    Z=R4+12lambda, lambda0=lambda+11epsilon lambda^2/6.

Here E=-g^ab E_ab and the disk conformal gradient is -E/(8piG).
Write T0=-Box>=0 on the regular-bolt physical scalar realization,
a0=1+11epsilon lambda/3 and Mstar^2=a0/epsilon. The exact scalar
diagonal in (4) is

    L=-(L_g+4E)/(8piG)
     =(1/(80pi^2))[(T0+4lambda)(T0-Mstar^2)
       +(1/36)div(Z grad)+(11/9)div(J grad)-Z/(3epsilon)]. (13)

L_g denotes the conformal derivative of E, not the action Hessian.
Its potential is 2R4+48lambda0-4E. The +4E fixed-density contact
cancels that last term exactly; dropping it would change (13).
For the explicit seed, with u=B/r,

    Z=(d/B^2)(4u^5-6u^6).

No constant-scalar-curvature assumption has been made.

For a frozen-background TT tensor h put r_dot=delta_h R4=-J:h.
The exact local derivative of the trace is

    Ttrace h=(1-11epsilon R4/36)r_dot
       +epsilon Box r_dot-epsilon h:Hess R4
       +(7epsilon/12)delta_h C4^2
       +(11epsilon/6)delta_h J^2-20epsilon delta_h M_em,
    delta_h J^2=J:Delta_L h-(R4/2)J:h-2h:J^2,
    delta_h M_em=-2h^ab F_ac F_b^c,                 (14)

where Delta_L h=-Box h-2Riem(h)+Ric o h+h o Ric. Indeed
delta_h Ric=Delta_L h/2 and
delta_h(Box R4)=Box r_dot-h:Hess R4; the derivative of h in the
Laplacian variation vanishes by its TT constraint.

For clarity, the Weyl variation in (14) is a local expression of
order two as well. With Schouten P_ab=(Ric_ab-R4 g_ab/6)/2, set

    K^a_bc=(nabla_b h_c^a+nabla_c h_b^a-nabla^a h_bc)/2,
    Rdot_abcd=h_ae R^e_bcd
               +g_ae(nabla_c K^e_db-nabla_d K^e_cb),
    delta_h C4^2=2C^abcd Rdot_abcd-8C^abcd h_ac P_bd
                         -4h^ae C_abcd C_e^bcd.     (15)

This follows by varying C=Riem-g wedge P and its four raised
indices; the contraction of C with g wedge Pdot is zero. In
particular it gives -4sigma C4^2 on a constant conformal variation
h=2sigma g. Equations (14)-(15) contain at most two derivatives of
h. The Weyl-invariant nonlocal determinant has zero compact
conformal trace; its nonlocal response remains in D.

The actual disk metric direction reconstructed from c is
k=h_c-L_eta g. Scalar covariance and its trace give
E'[L_eta g]=eta.grad E and tr k=-2div eta. Therefore

    Bmix c=-(1/(8piG))[Ttrace h_c-div(E eta_c)].     (16)

This is the exact off-shell return and density contact. It is not
obtained by setting the quantum gradient to zero. If H:c->h_c and
V:c->eta_c denote the frozen maps, its physical scalar transpose is

    Bmix* tau=-(1/(8piG))[H* Ttrace* tau+V*(E grad tau)]. (17)

The plus inside the bracket follows from integration by parts of
-div(E eta). H* here maps to the scalar c pairing; it is not a
silently substituted orthogonal tensor projection.

On the actual seed, in an orthonormal frame, Z=O(r^-5),
J=O(r^-4), C=O(r^-3), M_em=O(r^-4), and E=O(r^-4 log^N r),
with the retained differentiated bounds. Thus the local mixed block
has positive outer decay (at least two powers) and vanishes in the
AdS normal limit. The global York transpose retains its z^4 tails;
one cannot turn the local coefficient decay into an unrestricted
three-weight gain for every input. In the interior it has order at
most two before the order-zero reconstruction maps.

## 6. Actual scalar reduction and the adjoint source moment

Under Rcal the physical mass domain is
Z_(7/2) x [W_(7/2)+span{e_c}]. With the actual inverse (12), define

    S_eff=L-Bmix D_m^-1 Bmix*,
    f_eff=f_tau-Bmix D_m^-1 f_c.                    (18)

The full linear equation is exactly equivalent to

    S_eff tau=f_eff,
    c=D_m^-1(f_c-Bmix* tau).                        (19)

All sources and transposes are those of the retained full action.
On the compact scalar core the inverse is the symmetric critical
inverse, so

    0<=<tau,Bmix G_TT Bmix* tau>
       <=(1/.99)||Bmix* tau||_(E0-dual)^2.          (20)

This is a negative correction to L. It is relatively compact as a
map Z_(7/2)->Y_(7/2). In fact Bmix*:Z_(7/2)->Y_(7/2) is itself
compact: locally its order is at most two, so H4-to-H2 regularity
and Rellich give compactness into L2. In the end, its vanishing
mixed normal coefficients and controlled z^4 York tails give uniform
tail smallness. Compose this compact map with bounded D_m^-1 and
bounded Bmix on the mass-extended TT domain. The finite mass output
is finite rank. This is a statement about the stated domain/target;
it does not assume an unproved L2 norm on a negative Sobolev extension.

The exact triangular block elimination, or the established full index
and (12), gives index S_eff=-1. Its kernel maps isomorphically to
the physical coupled kernel by
tau -> Rcal(tau,-D_m^-1 Bmix* tau). In particular a coupled physical
kernel cannot be purely TT. No value of its dimension is inferred.

The actual adjoint phase also reduces without a vacuum-mode
substitution. A source-free adjoint field in the mass-orthogonal
space has scalar component psi with the accepted oscillatory
z^(3/2) phase and forced corrections. Its TT component is critical-L2,
and Bmix* psi is critical-L2 because of the two-power mixed gain.
Hence it is exactly -G_TT Bmix* psi. If a forcing has a term of
degree 7/2, use this critical inverse and extract its forced term;
do not apply D_m^-1 outside Y_(7/2).

Symmetry and the convergent weighted pairings then give the exact
source identity

    <chi,f>=<psi,f_tau-Bmix G_TT f_c>=<psi,f_eff>.    (21)

The integrations by parts have zero bolt contribution by regularity
and zero outer contribution by the mixed decay and the established
weighted tails. The [adjoint theorem](stationary-adjoint-phase.md)
still supplies k kernel moments and one additional scalar-phase
moment when the coupled kernel has dimension k. Neither psi nor
its actual moment is evaluated by (21).

## 7. Why the gap alone does not finish the scalar problem

Even the classical local piece is
B0 h=J:h/(8piG); the global York projection does not remove this
curvature coupling in general. On the magnetic RN seed with d=0,
J is tracefree and divergence-free because R4 is constant, smooth
at the bolt and decays as r^-4. Thus it is itself a nonzero TT
field. The actual d correction is retained in (14), not set to zero.

Using only the tensor L2 gap from (10), Cauchy--Schwarz gives

    <tau,B0 D^-1 B0* tau>
       <=11/(4*.99*pi*G*lambda) ||J tau||^2.        (22)

This estimate is valid using the TT tensor pairing and its pullback;
it does not alter the scalar normalization of D. Let
Q=(T0+4lambda)^(1/2), so ||Q^-1||<=1/(2sqrt(lambda)). Then (22)
only yields the global upper majorant

    ||80pi^2 Q^-1 B0D^-1B0*Q^-1||
       <=55pi ||J||_infinity^2/(.99 G lambda^2).    (23)

Dividing by Mstar^2=60pi a0/G gives
[55/(60*.99*a0)](||J||_infinity/lambda)^2, which is not small
in G lambda. The seed has two-plus-two eigenvalues +/-j, with

    B^2 j=-(1+3x+d)u^4+5d u^5-(9d/2)u^6,
    ||J||=2|j|.

Near its cold magnetic bolt ||J||/lambda is approximately six;
the majorant is approximately 33. This is a limitation of the
gap-only bound, not a lower bound on the actual Schur correction
and not an instability result.

After preconditioning (18) by Q and multiplying by 80pi^2, the
reference scalar problem is T0-Mstar^2. The remainder consists of
the explicit Z/J terms in (13) and the nonlocal correction in (18).
The next inverse argument must control this actual high-energy
weighted scalar problem, retaining TT derivative/fourth-order
coercivity, spatial coefficient commutators and the York tails.
Replacing Q^-1 by 1/Mstar at all frequencies is invalid. Even an
outgoing scalar inverse would leave the fast physical compatibility
moment to evaluate, since the physical domain excludes both
oscillatory amplitudes.

## 8. Verification and remaining scope

The [focused checker](../checks/stationary_weighted_tt_inverse.py)
audits finite normal congruence, scalar trace/density algebra, mass
normalization, graph-density scaling, block elimination and the
displayed rational bounds. Its immutable
[receipt](../evidence/stationary-weighted-tt-inverse.json) and
[raw log](../evidence/stationary-weighted-tt-inverse.log) record the
execution and exact input hashes. The weighted parametrix, adjoint
regularity, graph-density limit, asymptotic bootstrap and actual
inverse remain ordinary arguments, not outputs of a numerical solve.

This result fixes the TT component of the linear problem at the
analytic seed. It supplies an energy inverse constant, not a numerical
weighted coupled inverse constant. It does not transfer the result
to the finite nonlinear candidate or construct the exterior/state.
All S139.1-S139.2 and S074.2-S074.6 requirements retain their scope.

Remaining to-do list: determine the actual scalar Schur kernel and
adjoint phase, evaluate the retained source moments and quantify the
compatible weighted inverse; transfer regularity and close nonlinear
physical logarithmic phase or justify controlled spatial reduction;
match exterior/state and charge, transfer S138, and complete driven
stress, radius, controller, capture, tail, spatial area and S074-S080.
