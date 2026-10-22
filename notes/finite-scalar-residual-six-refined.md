# Refined actual scalar residual through six derivatives

2026-10-21, E003/S167. Ordinary argument before frozen verification.
Keep the original finite comparison, action, physical state/domain,
magnetic flux, mass and return. Refine the scalar estimate in
[the preceding proof](finite-scalar-residual-six.md). The new exact
rational ledger gives C_tau_new in

    ||f_tau||_(H_W6)/sqrt(gamma0)<=C_tau_new ell^2/g,
    g=G lambda, W=t^2, gamma0=4pi beta B^3.

The scientific receipt reports its size and comparison with the old
coefficient. No useful transfer certificate is assumed in advance.
The quantum/TT source, inverse/lift, nonlinear solve and stability
remain separate requirements.

## 1. Reuse each derivative order before taking maxima

Use the fixed reference gamma, D=z0 nabla_gamma, q=1+Zeta and
L=-log(1-Zeta) from the preceding covariant calculus. All weights
are fixed at the observation; no derivative acts on W or z0.
Use the exact radial arrays U,N from the frozen
[middle](middle-metric-ten-jets.md) and
[cold](cold-metric-ten-jets.md) receipts, through order eight.
These are original tangent fields u,nu, not finite-radius logarithms.

For a radial coefficient array r_j normalized by j!100^j define

    B(0,k)=a_k r_k k!100^k,
    B(n+1,k)=B(n,k+1)
              +n sum_(i=0)^k binom(k,i)c i!s^i B(n,k-i),
    T_j=B(j,0)/(j!100^j).

Middle uses a_k=1,c=6,s=200. Cold away uses a_k=1,c=15,s=200;
cold near uses the retained Cartesian factors a_k=A_k,c=100,s=20.
The radial-to-Cartesian proof applies to each separate r_k: its
Taylor remainder uses only the kth radial derivative. Its cold
weight is monotone on the segment from the bolt to the observation.
The connection recurrence therefore remains valid for nonconstant
arrays. Do not replace every r_k by the tenth-order maximum.

For V=nu/tau, tau=1-z, the earlier analytic inverse-tau disks give
normalized radial arrays

    v_mid,j=1200 sum_(i=0)^j (.2)^i N_mid,j-i,
    v_cold,j=1.2*82 sum_(i=0)^j (.2)^i N_cold,j-i.

The second formula uses a(s)/tau<82 before tensor conversion and
is an unweighted V bound. The regular even scalar V extends through
the bolt. Its near conversion uses the same A_k and connection;
there is no inverse power of the thermal gap in this calculation.

Multiply middle tensor arrays by 50.1^2 to obtain W-weighted
bounds. Multiply cold u,nu arrays by .00331*1.001^2, and cold V
arrays by 1.001^2. Omitting the factors 50.1^2 or 1.001^2 gives
unweighted bounds; cold u,nu still retain the factor .00331.
Since D=z0 nabla_gamma and z0<=1 these are also valid D-jet bounds.

For the outer part use the exact common_u, common_nu and common_V
ordinary remainder arrays from the [outer receipt](outer-metric-ten-jets.md).
At order j, let r be their maximum through j. The retained outer
tensor amplification a_outer,j gives the common weighted allowance

    9*.0331*a_outer,j*r+18e-11*.1.

Here z<Y<3z, Y^2 L(Y)<=.0331, and the final term retains the
mass cubic with coefficient below 1e-11. The same allowance is
valid unweighted. Take the maximum over the cover, then the minimum
with the preceding valid per-order bounds. Let U,K denote the
resulting unweighted and W-weighted arrays for both u and nu;
Vbar,KV similarly bound V. Multiply normalized entries by 100^j
to get divided D jets. Low entries keep the old contracts.

One additional sharp low input is useful. The
[retained gradient proof](finite-trace-gradient-and-projection-transversality.md)
gives |u_y|<3.2e-7/t, |nu_y|<7e-7 and |V_y|<8e-7.
Since |Df|=|f_y|/sqrt(C/x), C/x>5.9, the first divided coefficient
is below 4e-7 unweighted. On t<=50.1 its weighted value is below
.001. On t>=40, the old outer conormal tensor bounds give
W|Du|,W|Dnu|<3e-6, while the independent outer estimate
t^2|V_y|<4.3e-8 gives W|DV|<2e-8. Thus .001 is global. These
are bounds for first derivatives only; no unsigned inequality is
differentiated. Zeroth entries remain U0=3e-8,K0=1e-4,
Vbar0=1.2e-8,KV0=4e-5 whenever the regional minimum is larger.

## 2. Exact scalar-curvature remainder

Use the exact identity (7) of
[the curvature remainder proof](finite-comparison-scalar-curvature-remainder.md),
with its later corrected radius equation. Write

    R=R40/lambda=-12+Z0, n=u+nu, w=ell u+S,
    k=exp(ell V), T=1+tau(k-1), S=log T,
    B=k-1-ell V k, N_R=R_ell/lambda-R-ell dotR.

Then, exactly,

    N_R=R[e^-2w(1+2ell n)-1]+ell(e^-2w-1)dotR
      -6e^-2w[ell(k/T-1)L4nu+(B/T)L4tau
                                  +ell^2(k/T)tau V_y^2]
      +2e^-2w[e^(2ell u)-1-2ell u]/(xt^2).

This retains the cancellation between the Laplacian of log T and
its gradient square. It neither omits the remaining V_y^2 term nor
divides by the bolt speed. The common-coordinate trace equation
still supplies dotR=-Z0-H0+dotA0; use its previous W^2 array Zd.

Let Tp be the OLD finite-field unweighted array with its constant
removed, Up=U-U0 and Vp=Vbar-Vbar0. Those old S jets remain valid.
Define the following positive divided-jet series:

    O=1.001 exp(4Tp), E=exp(Vp+Tp),
    S1=1.001 E K, S1u=1.001 E U,
    S2=1.001 q exp(Vp+2Tp) K KV,
    W1=K+S1,
    J2=O[(4W1^2+2S2)(1+4U)+16K W1].

S1 and S2 have one and two weights, respectively. They follow from
S'=nu exp(ell V-S) and S''=z nu V exp(ell V-2S).
For J=e^-2w(1+2ell n)-1, J(0)=J'(0)=0 and

    J''=e^-2w[4(w')^2(1+2ell n)-8nw'
                                      -2w''(1+2ell n)].

Thus J2 bounds W^2 J'', and J is bounded by ell^2 J2/2.
All exponential operations act on zero-constant spatial jets;
their small actual zeroth values fit the 1.001 factors.

## 3. Seed differential operators and the six-term ledger

The seed physical Laplacian on an orbit scalar is

    L4 f=(C/x)(z/z0)^2[tr_P D^2 f+2Dphi.Df].

The canonical product connection has no mixed orbit/sphere slots.
Since |P|<2, positive arrays for W L4nu and L4tau are

    LN=6.1 q^2[2K''+2L' K'],
    LT=6.1[2q Bseed(q)+q^2 Bseed'(q)].

Bseed is the old positive degree-six polynomial for b. The latter
formula comes from L4tau=(C/x)(2zb-z^2 b_z), before absolute values.
All primes on composed series mean spatial series derivatives;
Bseed' alone is the polynomial derivative before composition.

An optional sharp zeroth replacement is LN0<=.1. On t<=50.1,
W|L4nu|<50.1^2*2.2e-6<.006. On the outer overlap, the old low
tensor remainder reserve 1e-7 and cubic give W|Dnu|<3e-6 and
W|D^2nu|<.0006, so 6.1[2*.0006+2*3e-6]<.008.
Hence .1 is a valid common reserve without changing higher entries.

Put A=1.001 q E KV, Q=.501 E KV^2. These bound
W(k/T-1)/ell and W^2 B/(ell^2 T). The identities
k/T-1=z(k-1)/T and B=-integral_0^ell r V^2 exp(rV)dr
retain the respective zeros. Since 0<=tau<=1, its jets are
majorized by q. Also V_y^2=(C/x)(z/z0)^2 |DV|^2.
The following six positive arrays therefore bound W^2 N_R/ell^2:

    N1=(12+Zabs(q))J2/2,
    N2=2O(U+S1u) Zd,
    N3=6O A LN,
    N4=6O Q LT,
    N5=6O*6.1*1.001 E q^3 (KV')^2,
    N6=4.004*1.3e-10 q^2 exp(4Tp+2Up) K^2.

Here 1/x<1.3e-10. In N2, Zd already carries both weights;
the exponential difference is estimated unweighted. The other terms
carry their two weights on the displayed field factors. This avoids
converting a generic second derivative of the entire curvature tensor.
Set Nnew=sum N_i. Taking coefficientwise minima with the old
NR2/2 remains valid and preserves every earlier improvement.

## 4. Complete residual and actual physical norm

Keep the previous reference-potential remainder RA2, Lie term and
anomaly-difference array DH unchanged. In particular DH still retains
the variation of the Laplacian acting on the nonconstant seed trace,
curvature contraction contacts, magnetic flux and its full degree-twelve
seed polynomial. The new total divided D-jet array through six is

    Pnew=min(Nnew,NR2/2)+RA2/2+Lie+DH.

Use the SAME actual-connection recurrence, six physical derivatives,
Jvol<1.001 and integral_(1+h)^infinity t^-2 dt<1 as before. Thus

    C_tau_new=(1001/24000) sum_(j=0)^6 Hnew_j.

Only a source upper estimate is refined. Smaller L2 and first-order
contracts remain available. At ell=1, using this scalar allowance
alone in the retained sufficient transfer target requires

    R0 R2 R4 <= 1e57 g epsilon_c/(800 C_tau_new),

before the nonnegative TT contribution. This is a condition on using
the allowance, not a lower bound on the residual or a nonexistence
result. Quantify it before choosing the next scientific step.

## Verification scope

Freeze this argument, checker and source review before execution.
Verify inherited hashes, nonconstant tensor recurrences, all six
remainder terms, the two parameter zeros, operator conversion,
per-order improvements and physical volume. Independently compare
sixth spatial coefficients of direct finite warped curvature minus
its seed and tangent with the exact remainder. Use two regular
manufactured base points, two nonzero path parameters and 60/90
digits, contour extraction versus real differentiation. These are
identity tests, not samples of the actual charged source. Preserve
all runs, failures and snapshots; no formal certification is claimed.

Remaining to-do list: assess useful scalar/source response bounds;
actual quantum/TT H4, higher lift/inverse and positive transfer;
nonlinear mapping, preserved neighborhood, signed compatibility
and existence; physical validity, state/boundary, matching, stability,
charge/return, opening, transport, preparation, creation and audit.
