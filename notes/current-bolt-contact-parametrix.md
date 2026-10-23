# Bolt adjoint with its exact fourth-order contact

2026-10-23, E003/S169. Ordinary argument before frozen execution.
The full stable semiclassical wormhole objective remains incomplete.
The closed feasibility audit and its numerical targets are unchanged.

The full physical adjoint residual cannot be computed as an L2 function
after subtracting only its leading rho^2 log log(1/rho) profile. This
calculation supplies an explicit contact-complete local TT parametrix,
retains the Einstein crossover, and bounds the entire high-frequency
angular-model mismatch. It does not compute the curved/global residual
or the physical correction. The most recent visible continuation only
verified already committed state; classify that continuation as no progress.

## 1. Derive the finite contact in the physical bolt variables

Keep the actual finite seed, c_R=1e22, q=2000^4, original mu0=1/b0,
prepared temperature, paired state, fixed outer clock and mass/log phase.
Use the physical variables and radial measure of
[the bolt-domain argument](current-bolt-adjoint-domain.md):

    L_b=ell_b(s(0)+xi(0)/2), ell_b=(1+h)/h,
    dmu=N (8pi^2 R^4/C) rho d rho [1+O(rho^2)],
    R=B(1+h), K_T=3 N C/20,
    K_s=6 pi N D_R C/G, D_R=96piG c_R-G/(60pi).

The leading high-frequency York return is S=xi/2, T=-3xi/2.
The exact canonical principal model has

    K(p)=C p^4 kappa(p)/(120pi),
    kappa(p)=4log p+2log C-4log(mu0 B)-62/15+o(1).

The angular/time and TT conversion multiplies K by 9pi N/2.
The EXACT Weyl Hessian in section 8 of
[the curved TT proof](stationary-curved-quantum-tt.md) contributes
-C(omega+1/6)(Delta T)^2/(30pi), with omega=log(r/B).
The other fourth-order Weyl square involves S+T/3 and has zero TT
principal part. The number 19/60 in that note is an inequality reserve,
NOT this action coefficient. At the bolt omega_b=log(1+h). Consequently
the complete frozen fourth-order TT symbol is

    K_T p^4(log p+a_b),
    a_b=log C/2-log(mu0 B)-log(1+h)-6/5.                (1)

Indeed 31/30+1/6=6/5. This agrees with the independently retained
physical outer normal contact after replacing its local physical scale
by R/sqrt(C). No scale or finite coupling has changed. The curvature
square has no fourth-order TT symbol: delta R[k]=-4tau xi. Its finite
potential, nonstationary contacts and all lower rows remain in H.

We retain a known positive second-order term in the local approximation:

    K_E=3 pi N F_lambda R^2/(2G), D=K_E/K_T,
    P_b(p)=K_T p^4(log p+a_b)+K_E p^2, p>=1.           (2)

K_E is the effective Einstein TT principal coefficient. This is not a
claim that it is the complete second-order coefficient of H. Other
second-order terms may be left in the residual without creating a
non-L2 point singularity. The actual parameter enclosure verifies a_b>0,
so (2) has no pole on the selected frequency interval.

## 2. A point-source parametrix, not a smooth homogeneous adjoint

Use the radial Hankel convention
u(rho)=int_0^infinity J0(p rho) u_hat(p) p dp. Its delta_rad has
transform one and pairing v(0) in measure rho d rho. Define

    E_T(rho)=(ell_b/2) int_1^infinity J0(p rho) p/P_b(p) dp. (3)

This integral is absolutely convergent for the value. Its derivatives
at the origin must be treated distributionally. Exactly in the frozen
model, P_b E_T=(ell_b/2)[delta_rad-J1(rho)/rho]; the quotient has limit
1/2 at zero. The lower cutoff creates a smooth remainder, not a wall,
a missing point source, or a changed thermal frequency prescription.
Here p is local spatial Fourier momentum, not a Matsubara label.

Multiply E_T by a smooth radial cutoff equal to one near the bolt and
zero before rho=.01. Pair it with the already established scalar
singularity -ell_b rho^2 log(rho)/(4K_s), with the same sort of cutoff.
The two physical delta coefficients are exactly ell_b and ell_b/2.
Apply the original global York lift; the induced outer radial tail is
retained. Compact support of a sphere component is not an extra outer
boundary condition. Both components lie in the local form space.

The exact constant a_b must be kept symbolically/correlatively when
subtracting H E_T. Even a nonzero rounded error Delta a leaves a
Fourier residual asymptotic to (ell_b/2)Delta a/log p. Its squared L2
norm has the divergent integral int p/(log p)^2 dp. Omitting the finite
contact has this same problem. A fixed relative numerical approximation
to the principal symbol does not certify an L2 residual either.
Interval evaluation encloses exact coefficients; it does not redefine
the operator by the midpoint. The prior uniform principal quadrature
remains useful for weak forms but cannot alone perform this cancellation.

## 3. Explicit all-angular asymptotic remainder

This is a remainder for the complete angular PRINCIPAL MODEL, not the
remaining curved/state-dependent Hessian. Put q=2000^4, t0=q^2 and

    F(t)=-2[psi(1+nu_+)+psi(1+nu_-)],
    nu_+/-=(q+/-sqrt(q^2-4t))/2.

For t>=t0 the two arguments are conjugate, have positive real part and
squared modulus t+q+1. DLMF 5.11(ii), with the first omitted Bernoulli
term and |arg z|<=pi/2, gives

    |psi(z)-log z+1/(2z)| <=sqrt(2)/(6|z|^2)<1/(4|z|^2).

It follows, keeping the conjugate reciprocal sum, that

    |-F(t)-2log t| <=(3q+5)/t.                         (4)

For t<t0 use the already proved logarithmic derivative bound
0<=d[-F(exp s)]/ds<=4. The difference in (4) is therefore bounded by
(3q+5)/t0+6log(t0/t). In the exact kappa integral split at
y0=4q^2/(C p^2)<=1. Integrating y and y log(y0/y), and bounding
sqrt(1-y)<=1 on the short interval, gives

    |kappa(p)-[4log p+2log C-4log(mu0 B)-62/15]|
       <= A/p^2+B_tail/p^4,
    A=10(3q+5)/C,
    B_tail=30[3+(3q+5)/q^2]q^4/C^2, p>=2q/sqrt(C).    (5)

The large but finite B_tail is an explicit bound, not an unknown
constant and not a violation of the physical correction target.
Equation (5) is an analytic continuum estimate. Checking frequencies
numerically is supporting verification, not its proof.

Let the exact angular model plus the same Weyl/Einstein terms be P_all.
Its normalized point residual on E_T has transform

    r_hat(p)=[kappa(p)-kappa_asym(p)] /
                         [4(log p+a_b+D/p^2)], p>=1.

For P>=max(1,2q/sqrt(C)), Hankel Plancherel and (5) imply

    int_P^infinity p |r_hat(p)|^2 dp
      <= [A^2/(2P^2)+A B_tail/(2P^4)+B_tail^2/(6P^6)]
                                      /[16(log P+a_b)^2]. (6)

Multiply the square root by ell_b/2 for the physical radial action
source coefficient. This norm uses rho d rho. It is NOT the physical
V* norm or a cold metric error; its measure conversion and the finite
frequency, curved, mixed, cutoff and global pieces remain necessary.
The finite frequency range is continuous and hence locally L2; we do
not assign its integral zero. All angular levels remain in (4)-(6).

## 4. Continuum local realization and its limit

This construction uses the smooth two-dimensional bolt interior, not a
smooth extension across the reflected outer mass cubic. The retained
finite-kernel identification must supply the actual local log symbol and
contacts. In normal Cartesian coordinates its coefficient change has
the structure O(rho^2) order-four-log plus O(rho) order-three-log plus
order-two-log. There are no odd invariant first jets at the smooth
rotational bolt. The frozen York principal lift and its lower symbols
must be included in this statement; it is not inferred from a resolvent
norm alone.

Applied to (3), the order-four-log kernel has off-point degree -2, the
order-three-log kernel degree -1, and lower terms are locally integrable
logarithmic functions. The displayed vanishing coefficients therefore
leave an L2 local remainder; multiplication kills the residual delta
at the point. Mixed rows have order at most two and likewise leave L2.
Cutoff commutators are smooth near the point. This establishes the
local regularity under the retained actual-symbol identification, not
a numerical norm for those curved or global remainders. The remote
operator and the outer mass/log phase still require their original
boundary realization. No closure of that outstanding quantitative step
follows from this argument or from the checker.

## 5. Executed quantitative tests to require before adoption

Substitute every seed parameter into a_b,K_T,K_E,D,ell_b for N=1,2.
Find the unique Einstein/log crossover p_cross by
p_cross^2(log p_cross+a_b)=D. This is a local approximation scale,
not a semiclassical validity cutoff or a physical pole.

At rho=0 formula (3) becomes the positive integral

    E_T(0)=ell_b/(2K_E) int_0^infinity
                         [1+exp(2s)(s+a_b)/D]^-1 ds. (7)

Enclose it by monotone left/right sums on [0,80], on 1024 and 4096
cells. The omitted integral is at most D exp(-160)/[2(80+a_b)].
Store every cell. Require the final relative enclosure width <.001.
This is the value of the LOCAL PARAMETRIX, not the actual adjoint value
and not a correction bound. An independent high-precision quadrature
must lie inside the enclosure. Evaluate (6) at P=1e12,1e16,1e20;
retain both unit-source and physical-source radial norms.

Check (4) directly by complex-ball digamma evaluation at t=q^2 times
1,2,100,1e8. Check the crossover signs using the complete parameter
box. Check both source normalizations and the exact fourth contact
algebra independently of the quadrature. All targets are fixed before
execution, and all attempts and source snapshots are retained.

## 6. Complete frozen-model mismatch, including finite frequencies

After the initial successful frozen run, extend and replay the complete
package to enclose the previously unevaluated finite-frequency MODEL
mismatch. This does not assert that the physical curved remainder is
small. Preserve the initial successful source snapshot and receipt.

Put p0=2q/sqrt(C), kappa_asym(1)=2log C-4log(mu0 B)-62/15,
kappa(0)=-F(0)-4log(mu0 B), and

    M0=|kappa(0)-kappa_asym(1)|+C/q+4log p0.

Indeed -F(t)+F(0)<=t sum d_n/A_n^2, and the telescoping identity
sum d_n/A_n^2=(2/q) sum_(n=1)^q n^-2<4/q gives
kappa(1)-kappa(0)<=C/q; the normalized weight mean y is <1.
For p>=1, 0<=d kappa/dlog p<=8, so
|kappa(p)-kappa_asym(p)|<=M0 for 1<=p<=p0.

Writing r_hat as in section 3, the entire interval p>=1 therefore obeys

    |r_hat| <= max[M0 p0^2/(4D), (A+B_tail/p0^2)/(4D)]. (8)

This is a relative quadratic-form error between the two declared
positive FROZEN radial multipliers. Require (8)<1e-12. It is not Z0
for H, and the all-angular quantum part alone is not being assigned
this relative error. The Einstein term and physical units are essential.

For an explicit complete L2 model residual, split at P1=1e16>p0.
Before multiplication by (ell_b/2)^2, its squared radial norm is bounded by

    M0^2(p0^6-1)/(96 D^2)
      +[A^2(P1^2-p0^2)/2+2A B_tail log(P1/p0)
         +B_tail^2(p0^-2-P1^-2)/2]/(16D^2)
      +the right side of (6) at P1.                   (9)

All terms are finite. This allowance can be large in radial source
units; it is neither an unknown constant nor a failed physical metric
test. The inverse/measure/dual conversion, cutoff, mixed rows and actual
curved response still need their own bounds. Store every contribution.

## 7. Research decision

Continue the physical weak primal-dual route with this analytic singular
enrichment. A leading-log profile or a rounded fourth contact is
insufficient for the proposed strong residual. The full angular model
is retained; (5)-(9) close its high-frequency point-residual integrability
and complete frozen-model mismatch with explicit tails. They do not supply the actual curved
remainder constant, the forward defect, X admission or nonlinear bounds.

The single next calculation is H applied to the cut-off two-component
parametrix together with the forward trial correction, including the
common curved/boundary remainder, to evaluate the actual two residual
norms and signed bolt pairing. The unchanged test is
|Q_app|+||d_b||_V*||d||_V*<1.667e-9, with exterior return <6.67e-10,
Y0<=5e-9 and every admission/nonlinear inequality. The full test remains
unexecuted. No solution or nonexistence theorem follows from this package.

Sources: [DLMF 5.11(ii)](https://dlmf.nist.gov/5.11#ii) supplies the
complex digamma remainder in its stated sector;
[DLMF 10.22(v)](https://dlmf.nist.gov/10.22#v) supplies the Hankel
transform convention after the usual square-root conjugation. Neither
source identifies the actual quantum boundary domain. The explicit
kernel/contact and phase arguments above remain required.

Remaining to-do list: full curved/boundary and finite-frequency residuals,
physical forward/adjoint certificate and X/admission; nonlinear closure,
exact compatibility, matching, causal realization, semiclassical validity,
all-sector stability, charge/return, opening, transport, preparation,
separate creation and the full stable-wormhole objective.
