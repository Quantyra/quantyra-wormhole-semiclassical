# A coercive certificate for the actual signed projection moment

2026-10-21, E003/S167. Ordinary argument before focused verification.

The signed projection coefficient can be bounded without numerically
resolving the fast adjoint oscillations: an exact identity transfers its
moment through a positive operator, retaining the complete TT coupling.
Three transfers give an explicit sufficient cold-source budget, including
errors in each auxiliary solve. This is a certification method, not an
evaluation of the actual coefficient. Its weighted residual inequalities
and domain conditions have not yet been verified for the actual source.

The [cold-source lower bound](cold-projection-source-lower-bound.md)
requires |beta_ell|<epsilon_c/(3e43) if the normalized subtraction is to
have cold-value norm below epsilon_c. The older unsigned bound
5.3e21 ell^2 does not establish that. The
[seed constraint obstruction](seed-residual-constraint-obstruction.md)
also prevents using reduced dynamic seed resonances as physical modes.
We return here to the unchanged static nonlinear problem and its actual
signed coefficient. No projection, physical boundary condition, action
coefficient or quantum state is altered.

## 1. Density-correct source and complete moment

Use the original comparison g_ell, full retained physical Hessian at
that comparison, and its fixed seed-density gradient. The bookkeeping
action used to construct f_ell is distinct from this full Hessian except
at ell=1. The target ell=1 is the physical retained equation residual.

Let t=r_seed/B, W=t^2, gamma0=4pi beta_period B^3,
g=G lambda and Jvol=dmu_ell/dmu0. In actual physical volume set

    (f_tau,f_c)=lambda^-2 Rcal_ell* (Jvol^-1 f_ell),
    p=80pi^2(f_tau-Bbar Dbar^-1 f_c),
    T=-Box_ell/lambda, H=T+4,
    E=60pi/g+11/3, e=E+4,
    V=(V_T,V_s), J=diag(I,-J_s),
    P=H(T-E)-VJV*.

These are the conventions of
[the full finite response](on-shell-finite-range-response.md).
The density conversion precedes the nonlocal York transpose and TT
inverse. It cannot be commuted past either operation. P is the actual
static scalar Schur operator, not the scalar vacuum comparison.

The real adjoint has chi_ell=Rcal_ell(psi,c_psi),
c_psi=-Dbar^-1 Bbar* psi and P psi=0. Its prescribed normalization and
the original projection give

    sigma_ell=L_ell^0(p_cut)>0,
    L_ell^0(f_ell)=lambda^2 <psi,p>_ell/(80pi^2),
    beta_ell=<psi,p>_ell/(80pi^2 sqrt(gamma0) sigma_ell). (1)

Equation (1) follows by the full block transpose, with both source
components retained. Setting f_c=0, replacing psi by phi_ref or removing
the returned mass component would generally change this coefficient.

The moment is initially defined on the retained weak-source space.
For an ordinary weighted L2 source q define

    N(q)=||W q||_L2(dmu_ell)/sqrt(gamma0),
    M(q)=<psi,q>_ell/(80pi^2 sqrt(gamma0) sigma_ell).

The following uniform numerical bound is useful:

    |M(q)| < N(q)/25.                                  (2)

Here and below this applies whenever q is in the stated weighted space;
ordinary L2 membership is not inferred from a weak-source estimate.

To verify the constant, use the two-channel bounds in
[the adjoint error theorem](on-shell-adjoint-error-and-trace-gradient.md)
and the actual d_tr<6e-4 from
[fixed-projection transversality](finite-trace-gradient-and-projection-transversality.md).
In their notation

    xi<30.1(6e-4)+1e-8=.01806001,
    sqrt(g) Gamma_RT<1.026 xi+3e-7(2520),
    sqrt(g) Gamma_Rs<4e-7 xi+1.000002(2520).

The scalar part alone has error at most
1.5 sqrt(g)Gamma_RT+3e-8 sqrt(g)Gamma_Rs+.0013<.031.
The reference D norm is below 4; converting its seed volume to actual
volume costs less than 1.000001. Thus

    ||W^-1 psi||_ell <4.04 lambda^-1/4.

Meanwhile sigma_ell>.048 sigma0 and
sigma0=lambda^-1/4 N0 with N0^2>8.7. Cauchy--Schwarz and

    4.04/[80pi^2(.048)sqrt(8.7)]<1/25

prove (2). The physical period, scalar normalization and lambda factors
are accounted for, rather than reset to one.

## 2. A positive operator with the full quantum coupling

Define the closed form

    K=H+VJV*/E
      =H+V_T V_T*/E-V_s J_s V_s*/E,
    D0=T(T+4).

The earlier positive elimination used the divisor e=E+4. Changing it to
E is harmless only after bounding its effect. The retained local block
obeys ||V_s*H^-1/2||^2/e<=2.6e-13 and e/E<1.000001. Hence

    K >= (1-3e-13)H+V_T V_T*/E,
    ||K^-1|| <=1/[4(1-3e-13)].                         (3)

The positive TT term is not treated as small. Its vertex is bounded on
the H form domain by the established mixed estimate, so adding its
positive form preserves closedness and an equivalent graph norm. The
negative local form has relative bound below 3e-13; the standard
[KLMN theorem, Teschl 6.24](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe2.pdf)
therefore supplies the stated closed positive realization. This is a
static auxiliary inverse, not a stable gravitational evolution.

The exact noncommuting operator identity is

    P=D0-E K.                                         (4)

In particular the source transfer is D0 K^-1. Replacing it by K^-1 D0,
or by T alone, would discard the coupling or change the order of maps.

## 3. Domains and an identity allowing inexact solves

Let C_reg consist of real smooth radial scalar functions that are
regular even fields at the disk bolt and compactly supported away from
the outer endpoint. The full adjoint equation gives

    <psi,P u>=0, u in C_reg.

The TT inverse and the fixed return are included in this Schur pairing.
Bolt parity and the regular measure eliminate the inner boundary term.
The compact scalar input gives no outer scalar boundary term; its TT
response retains the admitted mass and return tails.

The identity extends to any limit of such u in the norm

    ||H^1/2 u|| + ||W K u|| + ||W D0 u||,                (5)

because (2) controls both source pairings. Denote this certificate
domain by D_cert. It restricts test functions used to certify a moment;
it does not impose a new physical metric or matter boundary condition.
Neither K^-1:L2->D_cert nor an iterated weighted graph estimate is
asserted without proof. Smoothness in the open exterior alone does not
supply the endpoint control in (5).

Start from q0=p and choose trial functions u_j in D_cert. Put

    r_j=q_j-K u_j, q_(j+1)=D0 u_j.

Whenever these sources have the required weighted L2 norms, (4) gives

    q_j=E^-1 q_(j+1)+r_j-E^-1 P u_j.

Taking the ACTUAL full moment and iterating proves

    beta_ell=M(p)
      =sum_(j=0)^(n-1) E^-j M(r_j)+E^-n M(q_n),
    |beta_ell| < (1/25)[sum_(j=0)^(n-1) E^-j N(r_j)
                                             +E^-n N(q_n)]. (6)

Thus approximation and truncation errors remain visible. Exact solves
u_j=K^-1 q_j are a special case, provided they belong to D_cert; then
all r_j vanish. The sign of beta_ell is not determined by (6), and a
small nonzero upper bound is not the equation beta_ell=0.

For an approximate operator Ktilde the residual must still be measured
against K. Explicitly r_j=(q_j-Ktilde u_j)+(Ktilde-K)u_j. An error
budget for the actual nonlocal TT block and density/source conversion
is required. Computing a small residual for a scalar approximation or a
truncated angular tower would not satisfy (6).

## 4. A concrete three-step cold-value certificate

The benchmark has g<(4e-17)^2, hence

    E>180/g>E_min=1.125e35.

For epsilon_c>0 the following four actual norm bounds suffice:

| Source term | Required upper bound |
| --- | ---: |
| N(r0) | 2e-43 epsilon_c |
| N(r1) | 2e-8 epsilon_c |
| N(r2) | 2e27 epsilon_c |
| N(q3) | 2e62 epsilon_c |

Indeed the right side of (6), using E_min, is below
2.71e-44 epsilon_c, strictly below epsilon_c/(3e43). The early residual
requires the tightest absolute error because it has no E suppression.
This is a target for controlled computation or analytic bounds, not a
claim that ordinary floating point or the current source bounds meet it.

For three exact admissible solves, the looser single condition

    N(q3)<=1e63 epsilon_c                              (7)

also suffices: N(q3)/(25 E_min^3)<2.81e-44 epsilon_c.
This replaces direct phase sampling at kbar>4e17 by a weighted source
estimate after positive solves. It retains the full quantum operator.
No finite quantitative estimate for that actual N(q3) is supplied here.

The thresholds are for the cold VALUE norm identified in the earlier
necessary-budget theorem. They do not certify a C5 source norm, a cold
metric bound, an exact signed zero, or the nonlinear solution.

## 5. Why the new certificate still requires actual evidence

Positivity (3) alone gives no useful bound for D0 K^-1 in the required
weighted norm. Even in the uncoupled case K=T+4, D0 K^-1=T is unbounded.
On a resonant source Tq=E q, each exact transfer multiplies q by E;
the E^-n factor then gives no improvement at all. This control prevents
misreading (6) as automatic high-frequency suppression.

The existing comparison contract controls five metric derivatives and
the full equation residual in a low weighted norm. It does not by itself
prove membership or a norm bound for three applications of D0 K^-1.
Nor can six source derivatives be inferred by differentiating an unsigned
residual bound. The next task is to construct admissible trials and bound
the four ACTUAL quantities in the table, or prove the corresponding
weighted graph bounds for exact solves. Preserve cold, outer, density,
state and full angular cancellations throughout that calculation.

## 6. Verification contract

Freeze this note, checker and source review and reconcile the preceding
constraint-obstruction receipt and all transitive input hashes. Exact
checks cover the physical source transpose, the scalar-adjoint norm
budget, the E versus e form bound, telescoping with inexact residuals,
and the four numerical thresholds. Noncommuting finite block systems
compare the full adjoint/source pairing against Schur elimination and
the three-step formula at two precisions, including deliberately inexact
positive solves. Independent direct block solves verify K's positive
elimination. Omission controls detect the TT source, incorrect operator
order and discarded residuals. A resonant control exhibits no automatic
suppression, and an increasing-energy sequence detects unboundedness of
the source transfer. These are mechanism fixtures, not the charged
coefficient or a certificate for the actual source.

The abstract positive-form theorem is established mathematics. The new
contribution is the full density-correct charged moment identity and its
goal-specific quantitative budget. Continuum domains remain ordinary
arguments; no Lean campaign or formal continuum certification is claimed.
Every S167 acceptance task and physical gate stays open.

Remaining to-do list: actual weighted residual or graph certificate and
signed compatibility; cold/higher nonlinear response, preserved family
and nonlinear existence; physical constraint/boundary/state and validity
checks; stability, exterior matching, charge/return, opening, transport,
preparation, separate creation and full audit.
