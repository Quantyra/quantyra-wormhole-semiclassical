# Complete cold seed-source derivatives through eight

2026-10-21, E003/S167. Ordinary argument before supporting verification.
Retain the exact centered seed, original paired Dirac realization, HH
state, full angular tower and common physical subtraction scale. With
a_x=s_x[1+log(.01/s_x)], the complete canonical seed source satisfies

    |nabla^j(Theta_seed-Theta_product)|/(j!100^j q^2 a_x)<1e44,
    h<=s_x<=.001, 0<=j<=8.                              (1)

These are full physical covariant tensors, in the same unit-orbit
normalization as the preceding source results. Theta_product is the
homogeneous H2 source with the same sphere, masses and finite terms.
The new contribution is the actual-outer-domain versus global-cap
error E=Theta_seed-Theta_cap. Its corresponding allowance is 1e24.
Combine it with the explicit
[gluing ledger](seed-cold-gluing-eight-jets.md), whose largest entry is
below 2.37e43, and the
[local allowance](seed-cold-local-eight-jets.md) 1e8 to obtain (1).

This completes the cold C8 seed-source input. It does not bound ten
metric derivatives, the finite-metric H6 source or a nonlinear solution.
The smaller original complete C3 and endpoint C5 contracts are unchanged.

## 1. Correct principal Clifford convention before using locality

The preceding gluing note incorrectly calls its principal gamma matrices
anti-Hermitian. For the specified anti-self-adjoint Dsharp they must be
HERMITIAN, with {gamma_i,gamma_j}=2 delta_ij. This is the convention of
the [paired phase derivation](stationary-paired-dirac-phase.md) and
[actual-domain construction](stationary-reference-cold-bolt.md).
For example take gamma1=sigma1, gamma2=sigma2. Then

    Dflat=gamma1 partial_1+gamma2 partial_2,
    Dflat^*=-Dflat, Dflat^*Dflat=-Delta,
    Jspin=gamma1 gamma2, Jspin^*=-Jspin, Jspin^2=-I.      (2)

The canonical expression A=P^a partial_a+div(P)/2 is consequently
formally anti-self-adjoint when the real principal matrices are Hermitian.
Its massive inverse is the original Dsharp inverse, not an inverse
obtained by multiplying a self-adjoint operator by an untracked phase.

For a holomorphic normalizer with DG=lambda R_theta, K=(DG)^-1 and
Q=exp(theta Jspin/2), direct multiplication gives

    Q gamma^i K[a,i] Q^-1=lambda^-1 gamma^a.             (3)

This verifies the physical frame alignment, including its sign. Both
{Jspin,gamma_i}=0 and the Piola identity still give the half-divergence
formula used in the gluing proof. Its earlier algebraic fixture used
anti-Hermitian test matrices; that tested the algebraic cancellation but
did not fix the physical principal-symbol convention or rotation sign.
The new checker tests (2)-(3), density transport and the nonzero omitted
spin term with Hermitian matrices.

Every absolute matrix norm, derivative budget and inverse-composition
ledger in the preceding gluing result is unchanged by this correction.
The inverse estimates there were for the original anti-self-adjoint
physical realization. Read that note's principal matrices with (2);
its claimed norm bounds and its matching locality identity then apply.
The frozen note and run are preserved, with this correction recorded
explicitly rather than silently changing their source hashes.

## 2. Keep the actual inverse fixed

Use exactly the global C10 cap of
[the weighted inverse construction](seed-cold-weighted-inverse-jets.md).
It agrees with the seed to chi1=arcosh(.01/h); its supremum conformal
factor is <.03. Let kappa equal one on the inner side and decrease to
zero over the last chi width 1/4 before chi1. Its physical commutator
C_kappa=[Dsharp,kappa] has norm <9. For an observation x in the cold
region, its H2 distance from the support of C_kappa satisfies

    D_x>=.94 log(.01/s_x)-.3.                           (4)

Indeed dchi/ds=exp(-sigma)/sqrt(P), and the retained centered rational
guards give |sigma|<.051 and sqrt(P)<=1.006s. Thus dchi/ds>.94/s.
Integrate to s=.01 and subtract the cutoff width 1/4; .3 leaves a
strict margin for the centered perturbation. The transfer of these
guards is in the [centered source proof](stationary-centered-exterior-source.md),
section on the cold patch. This bound holds at the regular bolt too.

Write Rcap=(Dcap+M)^-1 and Ract=(Dactual+M)^-1, with M>18. The exact
first-order two-commutator identity is

    E_M(x)=-M tr[Rcap C_kappa Ract C_kappa Rcap] (x,x).  (5)

It has rank two per exterior under the retained paired normalization.
The original domain supplies ||Ract||<=1/M. Differentiating x in (5)
acts only on the two separated cap rows: the cutoff and actual inverse
are fixed. No high derivative of a reflected boundary or of a distant
cap taper is needed. Local counterterms cancel in this difference
because the two metrics agree on a neighborhood of the observation.

## 3. An explicit interior derivative estimate for either row

At each fixed x use its physical conformal normal chart w=(w1,w2),
whose factor sigma_N obeys |sigma_N|<.025 on the independent complex
z,zbar bidisk of radius 1/4. It has sigma_N(0)=d sigma_N(0)=0.
These are exactly the Goursat charts of the local source calculation.
They are regular at the bolt. The cap itself is held fixed when the
observation moves within this chart.

The global H2-cap Green construction bounds each plain separated row
by 10 exp(-M d/4) in its canonical H2 representation when d>=1.
Its proof uses only the small conformal mass-potential amplitude,
so the present cap with |sigma_cap|<.03 meets the earlier .051 guard.
Converting the observation to physical volume costs <exp(.03).
Allow also sqrt(2) to put its two spin components in a Hilbert tensor
norm; use 16 exp(-M d/4) for the physical row. The far argument and
its L2 measure are held fixed. The same bound holds for the adjoint row.

Consider the real disk |w|<=2/M. Its H2 distance from x is at most
2 exp(.025+.03)/M<2.2/M. By (4), the distance from the commutator
support remains >1 there. Remove the normal conformal spin connection
by setting psi=exp(sigma_N/2)r for either physical row (use the adjoint
equation for the incoming row). With X=Mw, it obeys

    (Dflat +/- V(X))psi=0, V(X)=exp[sigma_N(X/M)].       (6)

All estimates are Hilbert-valued: the target includes the other kernel
variable and spin index. For every real |X|<=2, a Cartesian complex
increment with |delta w_i|<=1/16 stays in the original z,zbar bidisk,
since 2/18+2/16<1/4. Tensor Cauchy therefore gives

    ||partial_X^r V||<=1.1 r!(24/M)^r
                         <=1.1 r!(4/3)^r, 0<=r<=10.  (7)

The same bound applies to the endpoint exp(-sigma_N/2).
The factor 24 exceeds 16sqrt(2) and includes all ordered coordinate
indices. It is a bound throughout the real disk, not just at its center.

Let R_n=2-n/12, and U_n=||partial_X^n psi||_(L2(B_Rn)), 0<=n<=10.
Use nested Lipschitz cutoffs between consecutive disks, with gradient
at most 24. Plancherel and (2) give ||grad f||_2=||Dflat f||_2 for a
compactly supported vector f. Differentiate (6) n-1 times and apply
that identity to the cut off derivative. The tensor product rule gives

    U_n<=26 U_(n-1)
       +1.1 sum_(r=1)^(n-1) binom(n-1,r) r!(4/3)^r U_(n-1-r),
    U_n<=n!50^n U_0, 1<=n<=10.                        (8)

For the induction divide by n!50^n. Its coefficient is at most
[26+1.1 sum_(r=1)^(n-1)(2/75)^r]/(50n)<1.
The smaller balls are nested inside those used by every lower order.
Local elliptic regularity, or the same argument with difference
quotients, justifies the derivatives before applying the estimate.
The coefficients are analytic on this interior patch; global cap
analyticity is not assumed.

For point evaluation use eta(X)=(1-|X|^2)^3_+, supported in B1.
It is C2, eta(0)=1, |grad eta|<=6 and |Delta eta|<=12. The Fourier
Cauchy-Schwarz estimate is

    |f(0)|<=(2sqrt(pi))^-1 ||(1-Delta)f||_2.

Apply it to eta partial_X^j psi for j<=8; two extra derivative
orders in (8) are available. In full tensor norms this yields

    |partial_X^j psi(0)|
      <=13U_j+12U_(j+1)+2U_(j+2)
      <5e5 j!50^j U_0.                                (9)

Here R10=7/6>1. The inequality uses
13+600(j+1)+5000(j+1)(j+2)<5e5 for j<=8.

The initial row bound, its distance loss 2.2/M, the normal endpoint,
and sqrt(area(B2))<4 give

    U_0<200 exp(-M D_x/4).

Multiplying back by exp(-sigma_N/2) and summing its product rule
costs at most 1.1/[1-(4/3)/50]. Equations (7)-(9) therefore imply
for each plain row, separately,

    |partial_w^j r(0)|<=2e8 j!(100M)^j exp(-M D_x/4),
    0<=j<=8.                                          (10)

The spin frame and physical density at the observation are included
by the exact conformal equation (6). No derivative of an adjoint with
complexified parameters is used. All observation derivatives here are
real, and the incoming row is estimated using its own adjoint PDE.

## 4. Preserve the cold weight and sum the complete tower

Insert both rows from (10), both commutators and ||Ract||<=1/M into
(5). There are j+1 Leibniz distributions. The trace M cancels the
inverse mass, so

    |partial_w^j E_M(x)|
      <=1e19 (j+1)j!100^j M^j exp(-M D_x/2).           (11)

The constant covers 2*9^2*(2e8)^2=6.48e18. Let t=log(.01/s_x)>=log10.
Since a_x=.01 exp(-t)(1+t), (4) gives

    exp(-M D_x/2)/a_x
      <=100 exp[-(.47M-1)t+.15M]/(1+t)
      <=(1000/3) exp(-.932M), M>=18.                   (12)

The logarithmic derivative of the middle expression in t is
-(.47M-1)-1/(1+t)<0. Its maximum is at t=log10. Use
1+log10>3 and .47log10-.15>.932. Thus the same cold factor survives
at every derivative order, including directly at s_x=h. No division
by the vanishing radial speed at the bolt is needed.

For the exact tower M_n>=18sqrt(n), d_n=2(q+2n)<=6qn. For 0<=j<=8,
M^j exp(-.932M) decreases on M>=18. Put A=.932*18=16.776. The
function n^(1+j/2)exp(-A sqrt(n)) decreases for n>=1. Hence

    sum_n d_n M_n^j exp(-.932M_n)<=q L_j,
    L_j=6*18^j*(6e-8)*[1+2/(16.776-j-3)].              (13)

To prove (13), bound the sum by its first term plus the integral,
put u=sqrt(n), and use u^(j+3)<=exp[(j+3)(u-1)]. The remaining
integral is 2exp(-A)/(A-j-3), with A-j-3>=5.776. The rational
exponential certificate gives exp(-16.776)<6e-8. This includes
the entire infinite angular tail; it is not a truncated simulation.

After the original source conversion C/q<.00352, define

    H_j=1e19 (j+1)(.00352)(1000/3)L_j.                  (14)

Then |partial_w^j E|/(j!100^j q^2 a_x)<=H_j. The absolute bounds
also justify termwise observation differentiation through eight.

## 5. Full tensor conversion and the completed regional input

In the fixed physical normal chart, the conformal Christoffel tensor
has full norm 2|partial sigma_N|. On the complex Cartesian bidisk
|w_i|<=1/16, Cauchy from the radius-1/8 bidisk gives |Gamma|<2.
At the origin Gamma=0 exactly. A second Cauchy estimate gives

    gamma_0=0, |partial^a Gamma(0)|<=gamma_a=2a!25^a,
    a>=1.                                            (15)

Let B(0,k)=k!100^k H_k and, for r+k+1<=8,

    B(r+1,k)=B(r,k+1)
        +r sum_(a=0)^k binom(k,a)gamma_a B(r,k-a),
    J_j=B(j,0)/(j!100^j).                             (16)

This retains all connection actions and nonsymmetric tensor components.
The physical metric at the center is I, so no further frame factor is
needed. Rational evaluation gives max J_j<1e24. These bounds are for
the actual-domain difference alone, in each point's fixed normal chart.

Now combine (16) with the preceding local and gluing PHYSICAL tensor
ledgers. Tensor norms obey the triangle inequality even though their
intermediate coordinate charts differ. At every order through eight,

    J_actual_domain+J_gluing+J_local
       <1e24+2.37e43+1e8<1e44.

This proves (1). There is no new finite subtraction or independently
renormalized angular mode. The
[middle](seed-middle-eight-jets.md) and
[outer](seed-outer-eight-jets.md) C8 seed-source bounds are already
available with their respective weights. All regional eighth-order
source inputs for the proposed metric reconstruction are now quantitative.

## Verification scope and remaining work

The supporting checker freezes the ordinary argument before execution.
It checks the corrected Hermitian symbol, spin alignment and density
formula; disk margins, all ten nested-energy steps, endpoint products
and point evaluation; monotonic cold-weight extraction and the infinite
angular moment estimate; every physical-tensor and combined regional
entry. Independent fixtures use a nonconstant scalar mass-potential
Dirac solution with Hermite derivatives and complex contour integrals,
and direct angular-tail quadrature at 60/90 digits. Neither evaluates
the actual charged kernel. The continuum proof is not machine-verified.
The prior convention mistake and its correction are retained explicitly.

Quantitative metric order ten still requires the radial reconstruction,
Ward pressure derivatives, regional weights and all coefficient products.
Finite-metric quantum response, tensor return and H6 Schur-source
estimates remain required before the target
R0 R2 R4 N6(p)<=1e57 epsilon_c can be assessed. That target has not
passed. Signed compatibility, nonlinear existence, physical state and
boundary admissibility, stability and traversability remain unproved.

Remaining to-do list: metric order ten and finite-metric quantum/Schur
estimates; positive transfer or a sharper actual transfer; cancellation,
signed compatibility and nonlinear existence; physical validity,
state/boundary, stability, matching, charge/return, opening, transport,
preparation, separate creation and audit.
