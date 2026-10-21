# A boundary singular sequence for the actual seed in volume norm

2026-10-20, E003/S168. Ordinary argument before focused verification.
The full wormhole goal remains incomplete.

There is a precise obstruction to one possible global continuation of
the annular theorem. In the spacetime-volume weighted Hilbert space,
the actual seed's full two-row retarded Hessian has a normalized compact
conformal singular sequence at every fixed causal frequency. Its two
output rows tend to zero, although the input norm is one and the inputs
escape toward the outer boundary. Thus no closed realization containing
these compact conformal tests has closed range and finite-dimensional
kernel in this norm. In particular it cannot have a bounded inverse.

This is an actual seed operator statement, not an actual growing mode.
The volume norm differs from the lapse-weighted norm of a wave evolution.
The seed is also off shell. Neither distinction may be dropped when
interpreting this obstruction physically.

## 1. Operator and norm being tested

Keep the action, common contacts, full Dirac tower, fixed state and
neutral spherical sector from the preceding notes. Use outer optical
coordinates (tau,Y), normalized so that the leading orbit is flat.
The actual seed has, with conormal derivatives through the needed order,

    Phi=-log Y+(1/2)log(alpha)+O(Y^2),
    h=log a=O(Y^2).                                      (1)

These follow from u'=sqrt(calB(u)/C), calB(0)=x,
calB'(0)=0, u=B/r, and the retained polynomial seed. The seed has
the required finite jets; no sixth derivative or nonlinear solution
is assumed. In a conformally flat orbit coordinate the same statement
holds after a change x=Y+O(Y^3).

Let sigma=gamma^2 for fixed complex Laplace frequency Re gamma>0.
The fixed-density polarized Hessian has radial density operator D(sigma),
with its two metric rows ordered (p,v). Normalize its output by

    K_gamma=Y^4 D(gamma^2)/C,
    H=L2((0,Y_*),Y^-4 dY; C^2),                         (2)

and extend the positive weight smoothly through the rest of the exterior.
Here D is the operator associated with the full quadratic form, so the
off-diagonal coefficient is half its cross term. The precise extension
of the norm away from Y=0 will not matter. Multiplying the weight by a
positive factor tending to a nonzero constant gives the same conclusion.
This includes the asymptotic actual spacetime volume per unit time.

Only the action of K_gamma on (p,0), p smooth and compact, is needed.
The canonical Weyl-invariant functional depends on v alone. Hence both
rows of this column are LOCAL and of order at most four, at arbitrary
causal frequency. Their support stays within supp p. The full nonlocal
shape-to-shape kernel is neither replaced nor estimated here; it is
multiplied by an exactly zero shape input. The causal conformal-column
identity and the fixed-density Hessian include all off-shell contacts.

## 2. The limiting column in logarithmic radius

Put t=log Y and p(Y)=Y^(3/2)u(t). This is an exact unitary map from
L2(dt) to the conformal component of (2). Write partial=partial_t and

    e=60pi/g+11/3, g=G lambda,
    k^2=e-9/4=60pi/g+17/12>0.

The zero-scaled-frequency normal operators are

    D0=Y^2 partial_Y^2-2Y partial_Y,
    A0=D0-4, B0=Y^2 partial_Y^2-12.

In measure Y^-4 dY the formal transpose is
B0*=Y^2 partial_Y^2-4Y partial_Y-6. Thus the normal conformal
column, after unitary conjugation, is

    P1(partial)=(partial^2-25/4)(partial^2+k^2)/(20pi),
    P2(partial)=(partial^2-2partial-45/4)
                                      (partial^2+k^2)/(60pi). (3)

The order of the mixed transpose matters. In particular, it is not
B0 acting on the conformal input. Both polynomials vanish at partial=ik.
These are precisely the full p,p and v,p columns of the previously
derived coupled normal action, not a scalar equation substituted for
the second metric row.

The actual column obeys

    U^-1 K_gamma U (u,0)
       = (P1(partial)u,P2(partial)u)
              + sum_(j=0..4) b_(j,gamma)(t) partial^j u,
    |b_(j,gamma)(t)| <= M_J exp(2t), t<t_J,             (4)

uniformly for gamma in any fixed compact J of Re gamma>0. Each b is a
two-component coefficient. Constants can depend on the retained model.

To justify the stronger coefficient statement (4), use the exact local
Einstein--Maxwell and Weyl densities in the local-frequency note, not
only convergence of annular forms. Their p,p and v,p Euler derivatives
use background jets through order four. Substitution of (1) into these
finite differential expressions gives (3). Every remaining conormal
coefficient is O(Y^2). The time terms carry gamma^2 Y^2 or gamma^4 Y^4;
sphere terms carry Y^2 and explicit magnetic terms carry Y^4. Deviations
of the radial jets and the normalized leading potential are O(Y^2).
The prescribed lambda0/lambda fixes the constant in e. An undifferentiated
Phi multiplying (Delta V)^2 contributes a logarithm only to the v,v
column; a conformal derivative removes it in the present column.
The smooth coordinate and density changes preserve the estimate.
No background equation is used to drop a residual. This also proves
(4) for any exact background with these same differentiated asymptotics,
conditionally on its existence and the same conformal state prescription.

## 3. A full two-row singular sequence

Choose chi in C_c^infinity(-1,1), ||chi||_L2=1. For L>=1 and R>L set

    u_(L,R)(t)=L^-1/2 chi((t+R)/L) exp(ik t),
    p_(L,R)(Y)=Y^(3/2)u_(L,R)(log Y).                  (5)

The H norm of (p,0) is one. Its support lies in
exp(-R-L)<Y<exp(-R+L), so each member is a compact interior test.
For j<=4, ||partial^j u|| <= sum_(a=0..j) binom(j,a)
k^(j-a) L^-a ||chi^(a)||; these norms stay bounded as L grows.
For either polynomial in (3), exact differentiation gives

    ||P_i(partial)u|| <= sum_(j=1..4)
                   |P_i^(j)(ik)| ||chi^(j)||/(j! L^j). (6)

There is no constant term because P_i(ik)=0. Equation (4) adds at
most M_J exp(-2(R-L)) times the sum of the four derivative norms and
the zeroth norm. Consequently

    sup_(gamma in J) ||K_gamma(p_(L,R),0)||_H
        <= C_(k,chi)/L+C_(J,k,chi) exp(-2(R-L)).        (7)

Take, for example, L_n=n and R_n=n^3+10 for large n. Supports are
eventually disjoint and escape toward zero, so the inputs converge
weakly to zero in H. Equation (7) tends to zero, uniformly on J.
No tail from the nonlocal shape kernel occurs: that kernel has no
conformal input column. This is therefore a singular sequence for
the FULL two-row expression, not a zero of a finite compression.
The constants at the very large benchmark k are not small; this is
an asymptotic theorem and does not certify a physical collar size.

## 4. Consequence for closed realizations

Let T_gamma be ANY closed densely defined realization on H that agrees
with K_gamma on these compact pure conformal tests. If T_gamma had
closed range and finite-dimensional kernel, the closed graph theorem
would give

    dist_H(w,ker T_gamma) <= A ||T_gamma w||_H.         (8)

Indeed the inverse on the range, with values in the orthogonal
complement of the kernel, has closed graph. Apply (8) to (5). Its
projection onto a finite-dimensional kernel tends to zero by weak
convergence, so the left side tends to one and the right side to zero.
This contradiction proves failure of the upper semi-Fredholm property,
and hence of Fredholmness or a bounded everywhere-defined inverse.
No self-adjointness of a retarded operator is asserted or needed.

Any compact bounded correction on H has the same singular sequence,
since it maps this weakly null bounded sequence to a strongly null one.
A finite-rank correction therefore cannot restore (8). Boundary rules
alone cannot restore it while retaining the compact conformal tests
and this norm. For a self-adjoint realization, if one exists, the usual
Weyl criterion says zero lies in its essential spectrum. For the
general retarded realization the direct closed-range statement above
is the claim; no ambiguous essential-spectrum convention is needed.

[Teschl, section 6.4, Lemma 6.17](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe.pdf)
provides the standard self-adjoint singular-sequence criterion and
compact-perturbation context. The two-row column, its actual seed
transfer and the closed-range argument are given here.

## 5. Why this is not yet a physical instability theorem

The spacetime-volume norm is not the evolution norm of a scalar wave.
In the leading AdS metric the Klein--Gordon spatial Hilbert weight is
sqrt(h_spatial)/N proportional to Y^-2 dY, whereas the spacetime
volume per unit time is proportional to Y^-4 dY. Thus for (5),

    ||p||_KG^2 = integral exp(2t)|u(t)|^2 dt
                         <= exp(-2(R-L)) ->0.          (9)

Also the wave pencil is obtained from D0+e-gamma^2Y^2 by multiplying
by Y^-2. That multiplier is unbounded on volume H. A small residual
in (7) supplies no small residual in the corresponding evolution
problem after normalization. Equivalently setting p=Y w yields
the spatial wave expression -partial_Y^2-(e-2)/Y^2 on L2(dY),
whose physical boundary realization remains a separate question.

The actual seed's density Hessian is off shell and still needs its
residual contacts to be related to the physical tensor equations.
The sequences do not supply charged constraint-satisfying initial
data, a pole, a nonlinear exterior, or a finite-time survival failure.
Their proper logarithmic phase scale remains
G k_proper^2=60pi+17g/12>180. No ultraviolet cutoff or order reduction
is adopted. Stronger stationary weighted source/range spaces, such
as those used for S166, are not contradicted by a theorem in H.

## 6. Decision and verification

Do not seek a global bounded inverse of the actual two-row density
operator in the equal volume weights (2), or try to repair it with
finitely many projected modes. That route fails for an explicit
boundary singular sequence. The annular form theorem remains valid;
it neither asserted this inverse nor bounded expanding logarithmic
windows. The next physical transfer must choose the lapse-weighted
evolution/source norms and metric boundary traces, and account for
the constraints and off-shell-to-on-shell change.

Verification differentiates the exact local action independently of
the factored normal column, checks its indicial roots and both mixed
transpose terms, tests variable-coefficient decay on a curved fixture,
and compares integrated compact-envelope residuals with independent
quadrature. Finite fixtures use moderate k and are not the benchmark
spectrum. Polynomial envelopes are H5 compact tests and can be
approximated in H4 by smooth compact tests. The continuum result is
the ordinary proof; no Lean or formal certification is claimed.

Remaining to-do list: physical evolution norms, metric boundary data
and constrained stability transfer or justified exclusion; S167 signed
source, cold/higher response, nonlinear solve and compatibility;
exterior/state, charge/return, opening, full-sector stability, transport,
preparation, separate creation and final audit.
