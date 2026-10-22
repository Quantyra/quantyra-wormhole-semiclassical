# Derivative budget for the actual strong Schur source

2026-10-21, E003/S167. Ordinary argument before supporting verification.
The endpoint C10 theorem is retained. It does not yet close the direct
H6 Schur-source estimate: the off-shell disk return contains a seventh
derivative of the full trace, hence an eleventh metric derivative.
There is also a fourth-order logarithm in the nonlocal quantum shape
response. The two issues require different estimates.

This audit identifies the missing derivative and a sufficient split
of the source norms. It supplies no numerical higher source/TT-transfer
constant, no obstruction for the actual correlated Schur source, and
no existence or stability theorem. The original full equations, state,
subtraction, mass column, source and outer-fixed return are unchanged.

## 1. The source can be split before asking for quantum H6

Use the physical finite comparison and the notation of the
[strong source theorem](finite-strong-schur-source.md):

    p=80pi^2(f_tau-Bbar Dbar^{-1} f_c),
    Bbar=Bmix/lambda^2, Dbar=D/lambda^2, g=G lambda,
    ||u||_(H_W^k)=sum_(j=0)^k ||W gradbar^j u||_L2,
    W=t_seed^2, gradbar=lambda^{-1/2} grad.

All norms below use the actual physical volume. The scalar source is
weighted; the TT source remains unweighted. The
[actual-domain theorem](finite-graph-admissibility.md) establishes
qualitatively Dbar^{-1}:H^4 -> H^(8,log) -> H^8 for this realization.
Let C_D,4 denote its as-yet unevaluated H^4-to-H^8 norm. Suppose the
following quantitative row estimate has been proved:

    ||Bbar c||_(H_W^6)<=C_B,6 ||c||_H^8.           (1)

Then the exact source satisfies

    N6(p)<=800/sqrt(gamma0) [||f_tau||_(H_W^6)
                              +C_B,6 C_D,4 ||f_c||_H^4]. (2)

The number 800 uses 80pi^2<800. No commutation of W through the TT
inverse, return or mixed row occurs. Thus the nonlocal shape forcing
need only be bounded in H4 for this route. Asking for its raw H6 norm
would impose an avoidable additional requirement. The local conformal
source still needs H6, including its fourth-order anomaly term.

The map converting the original fixed-volume residual to (f_tau,f_c)
must also be quantified in these norms. Equation (2) is a conditional
estimate for the unchanged source, not a computed upper bound for it.

## 2. A quantitative coefficient contract for the mixed row

The full finite mixed row is exactly

    Bbar c=-(8pi g)^{-1}[Tbar h_c-divbar(Ebar eta_bar_c)],
    Tbar h=lambda^{-1} Ttrace h, Ebar=Etr/lambda,
    eta_bar=sqrt(lambda) eta.                     (3)

Here h_c is the genuine physical TT lift and eta its fixed disk return.
The [weak row](finite-comparison-weighted-mixed-form.md) and
[strong row](finite-strong-schur-source.md) both retain the divergence
term. Tbar h is a local operator of order at most two on a TT tensor:

    Tbar h=sum_(j=0)^2 A_j dot gradbar^j h.        (4)

The contraction operator norms are included in A_j. In particular the
TT curvature identity delta R=-J:h is used before estimating (4).
The Box R variation is Box(-J:h)-h:Hess R; TT trace/divergence terms
vanish, but the disk return in (3) does not vanish off shell.

Define concrete coefficient and lift budgets

    C_loc=max_(j<=2,k<=6) ||W gradbar^k A_j||_infinity,
    C_E=max_(k<=7) ||W gradbar^k Ebar||_infinity,
    ||h_c||_H^8, ||eta_bar_c||_H^7 <= L8 ||c||_H^8. (5)

These numbers are not evaluated here. On the smooth regular core,
the ordered tensor Leibniz rule and sum_(i=0)^6 2^i=127 give

    ||Tbar h||_(H_W^6)<=381 C_loc ||h||_H^8,
    ||divbar(Ebar eta_bar)||_(H_W^6)
                                  <=381 C_E ||eta_bar||_H^7.

For the second estimate first expand into dEbar dot eta_bar and
Ebar divbar eta_bar. The scalar-vector contraction costs one and
the four-dimensional divergence contraction costs at most two.
No derivative of W is introduced: it stays on the coefficient.
Consequently a sufficient version of (1) is

    C_B,6 <= [381/(8pi g)] (C_loc+C_E)L8.          (6)

Completion on the existing core is justified by the actual graph-domain
theorem and the bounds (5); no new boundary realization is selected.
Combining (2),(6) identifies explicit numerical work still needed,
without declaring any of C_D,4,C_loc,C_E,L8 already bounded numerically.

## 3. The return reaches metric order eleven

For j=2,1,0 the coefficients in (4) contain at most two, three and
four metric derivatives, respectively. Six more derivatives therefore
reach order ten. This is consistent with the original metric result.

The full trace, however, contains the nonzero term

    Etr=R+12lambda0+epsilon B_an,
    B_an=(3/2)C4^2-(11/12)Euler4-20M_em+Box R,
    epsilon=G/(60pi)>0.

At the top radial order, the six-times differentiated return term is

    partial_rho^6 div(Etr eta)
                      =eta partial_rho^7 Etr + lower terms. (7)

The displayed contribution contains eleventh derivatives of the
metric through epsilon partial_rho^7 Box R. Its coefficient is not
eliminated by the TT constraint. In the local radial system the exact
equations are

    b'=-2A(b+c)-2Rwarp(b-c),
    eta'=A eta+b+c,
    A=a'/a, Rwarp=r'/r.                            (8)

Differentiating (8) expresses the needed derivatives of b and eta using
lower background jets. It does not express Etr' using only curvature
values or set it to zero. In particular the seed-only relation
Rwarp'=A Rwarp-Rwarp^2 is not an allowed finite-background replacement.

For clarity, the leading conformal symbol gives a local diagnostic
independent of a coordinate-count upper bound. In a flat observation
chart let g_v=exp(2v) delta, with v depending on one coordinate x.
Then exactly

    R_v=-6 exp(-2v)(v''+(v')^2),
    Box_v R_v=exp(-2v)(R_v''+2v' R_v').

At v=0 their variations are -6v'' and -6v''''. Let the local return
vector have nonzero constant coordinate component eta0, and set the
diagnostic trace E_diag=E0+R_v+epsilon Box_v R_v, with constant E0.
The curvature-square variations vanish at the flat basepoint; other
algebraic source variations can add only lower-order terms here.
For this diagnostic the linearized physical divergence is exactly

    delta div_v(E_diag eta0 partial_x)
          =eta0[-6v'''-6epsilon v''''' +4E0 v'].   (9)

The extra 4E0 v' is the physical volume contribution. Six derivatives
of (9), followed by the actual prefactor +1/(8piG) in Bmix, contain

    -eta0 v^(11)/(80pi^2).                        (10)

In particular the G in epsilon cancels this prefactor; its smallness
does not remove the derivative requirement. The tensor-trace part of
(3) has at most tenth background derivatives after six differentiations.
Local coefficients and variations of the return/lift from (8) do not
supply another eleventh-order term to cancel (10). A cancellation in
the full correlated Schur source would be additional work, not a
consequence of TT transversality or the weak divergence notation.

The test v_N=N^-10 cos(Nx) has all derivatives through ten bounded by
one. On a fixed interior interval, the six-times differentiated (9)
is eta0[6/N-6epsilon N+4E0/N^3]sin(Nx). Its L2 norm is unbounded for
fixed epsilon>0. A fixed cutoff equal to one on the interval preserves
the uniform C10 bound, up to a fixed constant. This is a local norm
diagnostic, not the actual comparison, a physical fluctuation, or a
claim that the actual higher derivative fails to exist. The actual
comparison is already known qualitatively to have all finite orders.

Integrating the divergence by parts changes the source norm. For
example a six-derivative bound on Etr eta controls its divergence in
H5, or in H^{-1} after six derivatives. It is not the H6 strong source
used in (2) and the existing three-transfer certificate. Any such
alternative must derive its matching dual-moment bound explicitly.

## 4. The quantum logarithm is a separate issue

The retained flat full-tower multiplier Khat obeys the
[spatial response audit](spatial-response-norm-audit.md):

    |Khat(k)|<=eta_G k^4 log(2+k),
    Khat(k)/(eta_G k^4)>=-|F_eff|/(60pi)
                              +(3/(512pi)) H_N, k>=4M_N. (11)

All unretained tower remainders are positive, so the second inequality
also holds at integer k_N=4 ceil(M_N). On a fixed 2pi torus take
v_N=k_N^-10 cos(k_N x1) and normalize the L2 norm of cosine to one.
The Fourier Sobolev norms satisfy

    ||v_N||_H^10=(1+k_N^-2)^5<=32,
    ||Khat(D)v_N||_H^6>=|Khat(k_N)|/k_N^4.       (12)

Equation (11) makes the latter unbounded. The same inputs have bounded
ordinary C10 norms. A finite common fourth-order contact changes the
constant in (11), not the growth. Thus raw C10/H10 control alone does
not justify a general quantum H6 estimate, even in the flat comparison.
This does not identify the curved operator with the flat multiplier.

The split (2) avoids that particular endpoint. In the flat Fourier
norm define H^(8,log) by the multiplier (1+k^2)^4 log(2+k). Then

    ||Khat(D)v||_H^4<=eta_G ||v||_H^(8,log),
    ||v||_H^(8,log)<=2 ||v||_H^9,                 (13)

using k^4<=(1+k^2)^2 and log(2+k)<=1+k<=2sqrt(1+k^2).
The existing endpoint C10 estimates are sufficient derivative inputs
to seek this lower-output-order curved quantum theorem; (13) alone
does not prove its curved, cold-weighted or actual-domain constants.

## 5. Selected next action and quantitative consequence

The direct strong route requires the following distinct inputs:

| Required estimate | Relevant highest derivative count |
| --- | --- |
| Local conformal source in H6 | Metric ten |
| Quantum shape source in H4 | Four plus a logarithm; seek metric nine or stronger input |
| TT inverse from H4 to H8 | Quantify the existing fourth-log realization |
| Local tensor-trace row in H6 | Metric ten |
| Off-shell return row in H6 | Trace seven, hence metric eleven in the direct ledger |

The [repaired metric audit](repaired-source-metric-transfer.md) stated
that ten was a necessary input, not a sufficient full H6 theorem. This
result identifies the additional off-shell requirement. The preceding
regional C10 and endpoint inverse results remain valid.

Select the direct extension to quantitative original-field C11, using
seed-source order nine in the second-order reconstruction, while
retaining the existing small low-order estimates. This does not
require replacing the quantum H4 task by a raw quantum H6 task, nor
silently declaring the auxiliary C11 coordinate map to be C12.
Then evaluate C_E through seven, C_loc through six, the needed lift
and TT inverse constants, the H6 local source and H4 quantum source.
If a sharper actual-source cancellation is used instead, prove it
before omitting (7). The original sufficient target remains

    R0 R2 R4 N6(p)<=1e57 epsilon_c.

It has not passed. The result of this audit is a corrected derivative
budget and a sufficient conditional norm split, not a finite source
bound or a new acceptance gate completion.

The checker retains exact conformal/divergence identities, omission
controls, symbolic derivative counts, rational product budgets and
independent periodic Fourier and physical-divergence diagnostics.
Continuum claims remain ordinary proofs; no actual charged kernel or
signed moment is sampled. Failed runs retain source snapshots.

Remaining to-do list: quantitative seed C9 and original metric C11;
finite-metric H4 quantum, H6 local/Schur and TT-transfer estimates;
signed compatibility, nonlinear existence, physical state/boundary
and stability; matching, charge/return, opening, transport,
preparation, separate creation and audit.
