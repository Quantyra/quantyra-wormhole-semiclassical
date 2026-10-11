# The first return-dependent photon contribution to clock stress

2026-10-11, S053. In the specified projected free-channel model, the
order-g^2 photon vacuum diagram's return-dependent part is finite and
reduces exactly to photon Green functions at the two mouths. Its energy
and first two static clock derivatives are bounded uniformly in the
number and largest length of returns. This is a component of the
renormalized stress calculation, not the full interacting clock kernel.

## 1. Definition and scope of the subtraction

Use q free NS Dirac channels with common throat `I=(-T/2,T/2)` and
lengths `L_j=T+D_j`. The positive exterior lengths D_j are fixed under
clock variation. Keep the Hermitian LLL vertices T_A, `A=(l,a)`, from
the [polarization calculation](propagating-photon-polarization.md), with

    kappa_l=l(l+1), tau_l=q beta_l/(4pi),
    beta_l=(q-1)! q!/[(q-1-l)! (q+l)!], 1<=l<q,
    sum_A tau_A/kappa_A=q(H_q-1)/(4pi).                       (1)

The photon has one of the nonnegative electric dual-scalar extensions
in [the open-load result](open-photon-load.md). At imaginary frequency
`k>0`, let `R_k=(H+k^2)^-1` and `G_k=Gamma R_k Gamma*`, where Gamma
takes the two core endpoint values for each canonical dual scalar.
The extension can mix angular modes. Require the local energy inequality
there, with fixed positive lower core mass and T in a compact positive
interval. Traces and the resolvent are interpreted through their forms;
Gamma* is a boundary functional, not an L2 delta function.

Let Pi_L be the exact free Euclidean current polarization, including its
local vector-current contact. Define Pi_infinity by unwrapping the free
fermion channels, on the *same core and with the same photon operator*.
Set `Delta Pi=Pi_L-Pi_infinity`. The object studied here is precisely

    E_w^(2)=(g^2/(2pi)) int_0^infinity Tr(D_E,k Delta Pi_k) dk. (2)

It is the difference of the first terms in the photon vacuum logarithms
`(1/2)Tr log(1+g^2 D_E Pi)` per unit Euclidean time. The coefficient
already includes g^2. It is not a definition of the entire renormalized
energy. The unwrapped diagram, background photon vacuum, physical
counterterm matching and terms of order g^4 and above are not set to
zero. Both terms in (2) use the same local current contact; it cancels
in their difference. The proof below makes this diagram difference
finite without assigning separate values to divergent traces.

No independent boson is assigned to every noncommuting angular current.
Only the exact free fermion two-point function enters this order.

## 2. Exact unequal-return kernel

For a particle in channel i and a hole in channel j, the right-moving
pair energies are `E_nm=2pi[(n+1/2)/L_i+(m+1/2)/L_j]`. Their Laplace
sum, including normalization, is

    F_ij(z)=1/[4 L_i L_j sinh(pi z/L_i) sinh(pi z/L_j)],
    F_infinity(z)=1/(4pi^2 z^2),
    f_ij(z)=F_ij(z)-F_infinity(z).                           (3)

The pole at zero cancels: `f_ij(0)=-(L_i^-2+L_j^-2)/24`.
The function is even, real on the real axis, analytic for
`|Im z|<min(L_i,L_j)`, and integrable on every fixed horizontal line
in that strip. For k>=0 put

    A_ij(k)=4 int_0^infinity cos(k t) f_ij(t) dt.             (4)

This converges ordinarily, including at k=0. In core potential
components `(A_tau,A_s)` the exact kernel difference for this pair is

    K_ij(k;s,t)=A_ij(k)
       [[cosh(k(s-t)), -i sinh(k(s-t))],
        [-i sinh(k(s-t)), -cosh(k(s-t))]].                   (5)

The full angular kernel is
`Delta Pi_AB=sum_ij (T_A)_ij (T_B)_ji K_ij`. Hermiticity follows
from A_ij=A_ji and the Hermitian vertices.

To derive (5), expand the filled-sea CAR contractions to obtain (3).
The density-density susceptibility is the cosine transform of the real
part of `F_ij(t+i(s-t'))`. Subtract (3) before removing the point split.
Moving its Fourier contour within the analytic strip multiplies the two
chiral branches by `exp(+-k(s-t'))`; their mean is the cosh in (5).
The two chiral current signs and the Euclidean spatial-potential factor
give the off-diagonal sinh entries. Alternatively the local Ward
identity fixes them after the common contact has canceled. The spatial
diagonal has the opposite sign. This also shows why keeping only density
response or omitting the contact before subtraction is incorrect.

For each pair (5) has spatial/potential rank at most two:

    q_+(s)=exp(ks)(1,-i), q_-(s)=exp(-ks)(1,i),
    K_ij=(A_ij/2)(|q_+><q_-|+|q_-><q_+|).                  (6)

The kernel is indefinite; no positivity of this state difference is
assumed. Its potential-component trace vanishes at coincidence. This is
not a claim that the unsubtracted polarization is finite rank.

Three useful independent limits of (4), with continuous k=0 values, are

    A_(L,L)=-k/[pi(exp(kL)-1)],
    A_(L,2L)=k[1-tanh(kL)]/(2pi)-sech(kL)/(4L),
    A_(L,infinity)=-log(1+exp(-kL))/(pi L).                  (7)

For example, the first follows from the partial fractions of csch^2.
For the second set z=pi t/(2L) and use
`csch(2z)csch(z)=-(csch'(z)+sech(z))/2` before Fourier transforming.
The last follows by summing the simple poles at `i n L`. The infinite
length here is a checked limit of the coefficient, not an assumed
normalizable finite-circle vacuum.

## 3. The diagram is an endpoint trace

Write `M=diag(kappa_A)` and

    C_AB(k)=sum_ij (T_A)_ij (T_B)_ji A_ij(k),
    Ctilde=M^(-1/2) C M^(-1/2),
    Q_k=[[1,-cosh(kT)],[-cosh(kT),1]].                       (8)

The Euclidean weighted potential propagator is
`M^(1/2)D_E M^(1/2)=1-L_k R_k L_k*`, with
`L_k u=(u',-ik u)` on the core. This is the Euclidean version of the
previous dual-scalar elimination. In the interior `L_k* q_+=L_k* q_-=0`.
In the full weak sense, however,

    L_k* q_+=exp(kT/2) delta_R-exp(-kT/2) delta_L,
    L_k* q_-=exp(-kT/2) delta_R-exp(kT/2) delta_L.            (9)

The endpoint terms must be retained. The trace of the identity against
Delta Pi is zero, while applying (9) to (6) gives exactly

    E_w^(2)=-(g^2/(2pi)) int_0^infinity
                    Tr_ports[(Ctilde(k) tensor Q_k) G_k] dk. (10)

The tensor is ordered by angular index and then left/right endpoint.
Mixed angular exterior loads are permitted in G. Finite-rank cyclicity
justifies this reduction at each k; no cyclic rearrangement of two
divergent individual vacuum traces is used.

For dual-Dirichlet core ends Gamma vanishes, so this diagram difference
is exactly zero, even with a varying positive core mass profile. This
does not make the full interaction energy or its static curvature zero:
the unwrapped reference diagram is still present. Open ends generally
give a nonzero result. The identity is therefore also a sensitive test
of boundary conventions in a stress calculation.

## 4. Uniform bounds and two static clock derivatives

Choose fixed `0<d0<=min(1,T_min)`, with `D_j>=d0`. Put
`b=T+d0/2`. Moving the whole-line transform of (3) to Im z=b gives

    |A_ij(k)| <= [1/(2d0)+1/(2pi b)] exp(-kb).              (11)

Indeed `|L sinh(pi(t+ib)/L)|^2`
`=L^2[sinh^2(pi t/L)+sin^2(pi b/L)]>=pi^2 t^2+d0^2`
for L>=T+d0. Thus the horizontal L1 norm of F_ij is at most
`1/(4d0)` and that of F_infinity at most `1/(4pi b)`.
Equation (4) is twice the whole-line transform, which proves (11).
The contour shift has no crossed poles because the zero pole is removed.

The endpoint trace inequality, applied to each core interval, gives
`||Gamma u||^2<=C h[u]`, with C independent of angular count and the
exterior. Hence `||G_k||<=C` for all k>0. Also

    ||Ctilde||_1 <= max_ij |A_ij| sum_A tau_A/kappa_A.        (12)

For (12), write C as the sum of rank-one angular matrices
`A_ij v_ij v_ij*`, with `(v_ij)_A=(T_A)_ij`; their weighted squared
norms sum to (1). Combining (10)-(12) with `||Q_k||<=2 exp(kT)`
already proves `|E_w^(2)|<=C g^2 q(H_q-1) d0^-2`.

For clock derivatives D=partial_phi, pull the core to a fixed interval
with endpoint values unchanged, and keep the exterior coordinates fixed.
Assume T is C2, `T_min<=T<=T_max`, and |DT|, |D2T| are uniformly bounded.
Let the positive Euclidean form `A_k(phi)=h_phi+k^2 norm_phi^2` have a
common domain, fixed endpoint map Gamma, and satisfy

    -c_j A_k <= D^j A_k <= c_j A_k, j=1,2,                 (13)

uniformly in k,q and the returns. This is a condition on the *whole*
photon realization, not just a Robin value at one frequency. It holds
for the scaled static core with a fixed nonnegative exterior used in
the checks. Its validity for physical moving mouths is unproved.

Differentiating the resolvent in forms gives
`G'=-Gamma R A' R Gamma*` and
`G''=Gamma(2R A' R A' R-R A'' R)Gamma*`. Equation (13) and the trace
bound imply `||D^j G||<=C_j`, j<=2, without a spectral gap.

Here are the needed differentiated contour estimates. If L=T+D and
z=t+i b(T), direct differentiation of `f_L(z)=1/[L sinh(pi z/L)]`
gives, for a+c<=2,

    |partial_L^a partial_z^c f_L(t+ib)|
      <= C_(a+c) d0^(-a-c) exp(-c_* |t|/L)
                                      /sqrt(t^2+d0^2).     (14)

The constants depend on T_min,T_max, not on D or L. To verify (14),
split |t|<=L and |t|>L. In the first part use the lower bound preceding
(12); each derivative contributes `pi coth(pi z/L)/L`, `1/L`, or
`pi z coth(pi z/L)/L^2`. Their magnitudes are bounded by C/d0 there,
as are the corresponding second-derivative factors by C/d0^2. In the
second part use the exponential lower bound on |sinh|; the at-most-
quadratic powers of |t|/L from L derivatives are absorbed by a smaller
fixed exponential. This proves (14) after increasing C. The reference
`1/z^2` and its first two z derivatives have uniformly integrable bounds
because b>=T_min. The chain and product rules therefore give

    int_R |D^r f_ij(t+i b(T))| dt <= C_r d0^(-1-r), r<=2,
    |D^j A_ij(k)| <= C_j exp(-kb)
                      sum_(r=0)^j k^(j-r) d0^(-1-r).        (15)

The moving contour contributes the derivatives of exp(-kb); it is not
held fixed while differentiating the lengths. The same k factors occur
when differentiating Q. Equations (12)-(15) provide an integrable
majorant for each differentiated trace in (10), including k=0. Thus
for j=0,1,2,

    |D^j E_w^(2)| <= C_j g^2 q(H_q-1) d0^(-2-j).           (16)

This proves differentiability under the integral, rather than inferring
stress control from the earlier current operator norm. It needs a
positive lower *exterior segment length* d0, not a gap or an upper
bound on the largest return. Collapsing D_j to zero is outside this
uniform statement.

For the inverse-tail family, D_min is proportional to delta. At
`g^2=q^-2`, `delta=q^-1/16`, a common bound on all three static
quantities is `O(q^-3/4 log q)`. Relative to the free damping-scale
stiffness q delta this is `O(q^-27/16 log q)`. This comparison concerns
only (2); it does not establish smallness of the omitted reference,
higher-order or nonzero-frequency stress terms.

## 5. Verification, sources and research decision

The [script](../checks/winding_interaction_stress.py) and
[receipt](../evidence/winding-interaction-stress.json) retain checks of
the pair coefficient against the three limits (7), contour shifts,
direct regulated pair sums, the rank and boundary identities, angular
weights, independent continuum photon propagation, and static energy
and clock-derivative refinements. Numerical error criteria are specified
in the script; they do not replace the uniform ordinary proof.

[Blaschke, Carballo-Rubio and Mottola, 1407.8523v2](https://arxiv.org/html/1407.8523v2)
develop the known cylinder-current/Schwinger-term and intermediate-pair
framework, especially sections 2.3-2.4 and Appendix A. Their common
cylinder and Abelian bosonization do not solve the unequal-return,
noncommuting interacting theory here. [Georgi and Warner,
1907.12705v2](https://arxiv.org/html/1907.12705v2) remain the exact
commuting massive-vector baseline. The use of current Ward identities,
imaginary frequency and a shared reference subtraction is established.
The source supplement records reading depth and version identities;
this component is not accepted as an independent novelty endpoint.

**Decision:** continue the physical damping primary and matched
finite-duration transfer fallback. The first return-dependent photon
diagram and its static clock derivatives are now controlled in the
declared model. Next include the renormalized unwrapped reference,
bound higher current/vacuum diagrams and the dynamical stress response,
and close the physical coupled geometry and state. The static derivative
bound alone cannot be inserted as a full retarded clock self-energy.
No physical damping, finite-disturbance survival, global passage or
travel-time advantage is proved. No Lean/Lake or GCP computation occurred.

Remaining to-do list: finish the renormalized reference and higher-order
static/dynamical stress, match the coupled geometry and state, complete
novelty and the consequential result, then independently verify,
GCP-certify suitable endpoints and prepare the reproducible manuscript.
