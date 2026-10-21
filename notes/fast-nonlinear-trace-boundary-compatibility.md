# Nonlinear trace compatibility of the fast boundary correction

2026-10-21, E003/S168. Ordinary argument before focused verification.
The full stable-wormhole goal remains incomplete.

The second-order conformal tail that cancels the preceding cutoff
oscillation does not satisfy the retained trace equation. Conversely,
the conformal tail required by that equation leaves a strictly nonzero
oscillation in the action for every positive retained coupling. Thus a
conformal-only repair fails already at a necessary bulk equation, within
the explicit asymptotic expansion below. This does not exclude a coupled
metric completion or establish a physical instability or no-go theorem.

We also derive the necessary second-order trace relation including the
independent optical shape. The full quantum shape equation, its boundary
variation, and global matching remain to be solved. No shape correction
is adopted merely because it satisfies this one local equation.

## 1. Unchanged action, regulator and expansion

Keep the action, paired matter domain, magnetic flux, subtraction and
fixed background cut of the [finite-cut test](fast-trace-finite-cut-action.md).
Use z=log(r/B), t=epsilon lambda>0, epsilon=G/(60pi),
k^2=1/t+17/12 and S0=lambda B^3/G. In the optical metric

    g=exp(2Phi)[exp(2V)(d tau^2+dy^2)+B^2 dOmega^2],

write w=Phi-Phi_seed and v=V-V_seed. Consider the second-order outer
expansion

    w=eta exp(-3z/2) f(z) H(tau)
                         +eta^2 exp(-3z) h(z) H(tau)^2,
    v=eta^2 exp(-3z) j(z) H(tau)^2,
    f=a cos(kz+phi), a>0.                                  (1)

H is a fixed nonzero thermal cosine. The displayed coefficients can
also be used for the static trace calculation, but the preceding static
linear action failure already prevents a finite static action on the
simple path. Interpret (1) as an expansion through second amplitude
order with differentiated remainders smaller than its displayed powers.
It is not a theorem that every exact solution has this expansion. The
conformal-only test sets v identically zero. The coupled equation below
only supplies a necessary condition on j in this proposed expansion.

Every radial derivative used below is a z derivative. At fixed frequency,
paired time derivatives gain two radial powers. Charged coefficients
F=lambda r^2[1+O(r^-2)] and the sphere curvature also contribute lower
powers. The normalized second-order trace coefficient therefore has the
same leading equation as the radial normal geometry. Corrections are
O(exp(-2z)) after the second-order exp(-3z) factor is removed, for the
explicit charged seed and smooth profiles used in the fixtures. No
uniform estimate in increasing frequency or t approaching zero follows.

## 2. Derivation from the action and from curvature invariants

Divide the leading static z density by S0 exp(3z). Put
u=1+w'+w'', q=(1+w')^2 and b=v'+v''. The retained local bulk gives

    L_C={exp(2w)[6q+4(1+w')v']
                       +6(1+11t/6)exp(4w+2v)}/4,
    L_W=15t exp(-2v){(z+w)b^2/15-u b/15-u^2/10
                              -b q/15-7u q/90-7q^2/180}.   (2)

These follow from the [full optical density](local-coupled-frequency-polynomial.md),
with no change to the finite-cut boundary action. Intrinsic lower-power
terms are omitted only from this leading coefficient. The independent
canonical determinant has no conformal derivative; it has not been
discarded from the independent shape equation.

For E=exp(-3z) delta integral exp(3z)(L_C+L_W)/delta w, the Euler
operation is partial_w-(D+3)partial_w'+(D+3)^2 partial_w'', D=d/dz.
The background value is zero. Its linear symbols are

    P_w(s)=-(s-1)(s+4)(3t s^2+9t s+11t+3),
    P_v(s)=-(s-3)(s+4)(3t s^2+9t s+11t+3)/3.             (3)

In particular s=-3/2 +/- i k annihilates the first conformal order.

There is an independent local check. The [complete retained trace](charged-conformal-response-and-constraints.md)
is

    T=R+12lambda0+epsilon[7C^2/12+11J^2/6-11R^2/72
                                               -20M+Box R]. (4)

Here J is traceless Ricci, M=F_ab F^ab, and lambda0 includes its
unchanged vacuum shift. In the normal geometry the density identity is
E=exp(4w+2v)T/(2lambda). For lambda=1, an independent warped-metric
calculation gives

    R=2 exp(-2w-2v)[-v'-v''-3(w')^2-9w'-3w''-6],
    C^2=(4/3)exp(-4w-4v)(v'+v'')^2.

Ricci eigenvalues and Box R then recover (3) and the quadratic term
below directly from (4), without using an action variation. The checker
keeps the density contact exp(4w+2v); the seed is not declared on shell.
The underlying conformal Ricci, Paneitz and anomaly identities are
standard [Mazur--Mottola, (38), (40), (44), (46), (54)-(56)](https://arxiv.org/html/hep-th/0106151v1).
Their discarded bulk surface terms do not replace our retained boundary
action or select an ensemble.

## 3. The forced second-order tail

After inserting the first-order f, the normalized quadratic forcing is

    N_t(f,f')=(471t/4+289/4+6/t)f^2
                   +(35t/3-5)f f'-(21t+3)(f')^2.          (5)

Define

    L_h(D)=-(D-4)(D+1)(3t D^2-9t D+11t+3),
    L_j(D)=-(D-6)(D+1)(3t D^2-9t D+11t+3)/3.

The necessary leading equation is

    L_h h+L_j j+N_t(f,f')=0.                              (6)

Take a=1 for the coefficient formulas; every second-order coefficient
below scales by a^2. In the phase basis 1, cos(2kz+2phi), sin(2kz+2phi),
the forcing coefficients are

    n0=(88t^2+47t+3)/(2t),
    n_c=(295t^2+195t+18)/(4t),
    n_s=k(5/2-35t/6).

For j=0 the mean is forced to be

    h0=-(88t^2+47t+3)/[8t(11t+3)].                      (7)

It is not a free choice of action zero. On the doubled-frequency pair,
L_h is the matrix

    L=(5t+3)[[(19t-12)/t,-42k],[42k,(19t-12)/t]],
    det L=4(5t+3)^2(11t+3)(65t+12)/t^2>0.               (8)

The unique periodic trace-compatible pair is h_T=-L^-1(n_c,n_s).
The homogeneous exponents for h are -1,4,3/2 +/- i k. Thus a bounded
coefficient solution of this normal forced equation has this periodic
part, with only the decaying exp(-z) homogeneous addition. There is no
bounded homogeneous freedom that can change its persistent oscillation.
This statement concerns the normal coefficient equation, not global
existence or an estimate for the full nonlinear map.

## 4. Incompatibility with the action cancellation

The preceding exact finite-cut calculation gives the quadratic cutoff
term Q_t+(3t/7)(7h''-33h'+20h), where

    Q_t=-(3001t+546)f^2/504-(71t+21)f f'/14
                                      -(167t/42)(f')^2.

Its doubled-frequency vector and h operator are

    q=(-9(t-9)/56, k(71t+21)/28),
    M=[[-(59t+84)/7,-198kt/7],[198kt/7,-(59t+84)/7]],
    det M=4(65t+12)(227t+147)/49>0.

The action-cancelling vector is h_A=-M^-1 q. Substituting the trace
vector h_T instead leaves R=q+M h_T with cosine coefficient

    R_c=-(78815t^3+74391t^2+16137t+81)
                              /[168(5t+3)(11t+3)]<-3/56. (9)

To check the strict bound, subtract 9(5t+3)(11t+3) from the numerator;
the remainder is 78815t^3+73896t^2+15705t>0. Thus no radial phase
rotation can make the remaining oscillatory vector zero. The quadratic
action on the trace-compatible conformal path still has a nonconvergent
term S0 <H^2> a^2 [R_c cos(2kz+2phi)+R_s sin(2kz+2phi)]. Its amplitude
is greater than (3/56) S0 <H^2> a^2. A constant h or the decaying
homogeneous correction changes no doubled-frequency coefficient.

Conversely, the action-cancelling h_A has trace residual L h_A+n whose
sine coefficient is

    -k(61175t^3+51351t^2+14625t+3969)
                                      /[12t(227t+147)]<0. (10)

No mean correction can cancel it. Independent quantum state changes
have no adjustable bulk trace at fixed metric and action, so they cannot
repair this conformal-only mismatch. The required mean (7) grows as
-1/(8t) when t tends to zero; no small nonlinear neighborhood at the
physical benchmark has been proved by this asymptotic computation.

## 5. The coupled obligation and the next decision

Allowing shape j changes (6). Its doubled-frequency matrix has determinant

    det L_j=4(5t+3)^2(11t+3)(125t+12)/(9t^2)>0.

If h_A were held fixed, (6) would require the nonzero shape pair
j=-L_j^-1(L h_A+n). The mean equation is
(12+44t)h0+(6+22t)j0+n0=0. These are necessary trace conditions only.
Adding j also changes the action and the canonical determinant, so the
preceding conformal action cancellation cannot be reused unchanged.
There is no claim that these algebraic pairs solve the coupled boundary
problem or belong to the physically admitted phase.

The next calculation must match the complete shape equation and boundary
variation, at both the static and doubled thermal frequencies generated
by H^2=(1+cos(2gamma tau))/2. A positive-frequency auxiliary inverse does
not settle the static signed compatibility condition. The full canonical
response, horizon/outer matching and constraints must remain present.

Focused verification derives the leading operators twice, compares the
two forced tails exactly, and checks the predicted mismatch against the
unexpanded trace of charged finite-radius conformal metrics at refined
radii and two precisions. The charged fixtures test the leading trace
limit; they are not an exact solution, calibrated physical benchmark,
full dynamic response or stability calculation. No Lean run is selected.

This closes the proposed conformal-only second-order repair in (1).
It leaves coupled nonlinear boundary completion open. No physical fast
branch exclusion, exact wormhole, instability or impossibility follows.

Remaining to-do list: coupled nonlinear gravitational boundary completion
and static source compatibility; continuous causal response and
constrained stability; nonlinear existence and signed compatibility;
exterior/state, charge/return, opening, all-sector stability, transport,
preparation, separate creation and full audit.
