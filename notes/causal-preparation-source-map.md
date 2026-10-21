# Causal preparation, initial terms and thermal source compatibility

2026-10-21, E003/S168. Ordinary argument before focused verification.
The stable-wormhole goal remains incomplete.

The fixed-HH exponentially prepared response and a finite-time initial
value problem have different source maps. We derive the latter, including
the initial metric amplitude in the moving-domain Dirac representation,
the local fourth-order metric jets, and the affine curvature reconstruction.
For compact initial metric jets the reconstruction offset contributes a
range element, so it changes the reconstructed solution but not the full
coupled cokernel overlaps. Initial constraint data still have to be imposed.

The [preceding conserved-source result](ward-sources-and-fast-compatibility.md)
does not identify arbitrary conserved stresses with quantum preparation
errors. In particular its pure thermal bolt condition is not automatically
a boundary condition on the Laplace transform of a finite-time initial
pulse. This distinction is exhibited explicitly below.

## 1. The actual independent covariance source

Keep the fixed self-adjoint zero-source Dirac Hamiltonian H_0, its paired
domain, magnetic bundle and prepared covariance C=f_beta(H_0). Retain all
angular sectors and the common subtraction. Use the fixed-domain vertices
A_i,D_i,T_ij from the [ADM derivation](adm-dirac-boundary-response.md).
For amplitudes lambda_i(t), the first Hamiltonian variation is

    delta H=sum_i [A_i lambda_i+D_i dot(lambda_i)].

Let c(t)=delta C(t), c(0)=c_0, lambda_i(0)=lambda_i0. The exact linearized
covariance equation and the interior generalized force are

    dot c+i[H_0,c]=-i[delta H,C],
    delta F_i=Tr(A_i c)-Tr(D_i dot c)
                             +sum_j Tr(C T_ij)lambda_j.       (1)

Temporal endpoint work is not discarded by writing this interior force.
For Re s>0, put R_s=(s+i ad_H0)^-1. The Laplace transform of (1) is

    c_hat=R_s{c_eff-i sum_j[A_j+s D_j,C]lambda_hat_j},
    c_eff=c_0+i sum_j[D_j,C]lambda_j0.                        (2)

The PLUS sign in c_eff comes from transforming dot(lambda), not from
an independently changed initial state. The full force becomes

    F_hat_i=sum_j K_ij(s)lambda_hat_j+J_i^C(s),
    K_ij=Tr(C T_ij)
         -i Tr[(A_i-s D_i)R_s[A_j+s D_j,C]],
    J_i^C=Tr[(A_i-s D_i)R_s c_eff]+Tr(D_i c_0).               (3)

The last contact contains c_0, not c_eff. The readout uses A_i-s D_i,
whereas the drive uses A_j+s D_j. These are bilinear Laplace identities;
sesquilinear versions conjugate the test frequency appropriately.
Equations (2)-(3) extend the retained zero-initial-amplitude formula.

At a common finite regulator these traces are literal. In the continuum
they are distribution/form identities where the same renormalized
response and the admitted initial state are defined. They do not prove
that arbitrary c_0 is positive, Hadamard, horizon-regular or compatible
with a changed initial metric. Separate unsmeared traces need not exist.
No new initial-surface subtraction is selected. An ultraviolet mismatch
in an arbitrarily chosen initial state is not cured by (2).

For the optical two-field chart p,v, the shift and boundary transport
vanish: D_i=0. The canonical conformal vertex is zero; the shape vertex
is the retained beta_D m_nu v. Thus an admitted independent c_0 gives
only the canonical shape readout Tr(B_v R_s c_0); the local conformal
anomaly is already in the metric operator. At fixed geometry and domain
the independent bulk trace is zero. Its normalization into (f,g) must
use the original source work and the chosen equation-side sign.

The selected zero-initial-state metric response has c_0=0. Its K_ij term
is part of the operator and must not be reused as an independent source.
This separation agrees with [Anderson, Molina-Paris and Mottola, section I](https://arxiv.org/html/0907.0823v2).
Their calculation does not supply this charged initial-data map or its
continuum estimates. No short-scale branch prescription is adopted.

## 2. No freely growing independent covariance mode

For zero metric drive, the exact solution is

    c(t)=exp(-i H_0 t)c_0 exp(i H_0 t).                      (4)

Unitary conjugation preserves the operator norm and every finite Schatten
norm. Its right-half-plane resolvent is

    R_s c_0=integral_0^infinity exp(-s t)
                        exp(-i H_0 t)c_0 exp(i H_0 t)dt,
    ||R_s c_0|| <= ||c_0||/Re s.                            (5)

For bounded c_0 the integral can be taken weakly; for trace-class or
Hilbert-Schmidt c_0 it converges in that norm. Differentiation on the
generator domain, or integration of weak matrix elements, proves (5).
If a bounded independent covariance has c(t)=exp(s t)c_s with Re s>0,
then norm invariance forces c_s=0. Equivalently s+i ad_H0 has no such
homogeneous null vector. This does not assert that every singular
state tangent is bounded, or that the coupled metric system is stable.

In the fixed-HH exponential response prepared from t=-infinity, therefore,
there is no freely chosen independent growing state source. The covariance
is the driven term in (2). With no additional drive, the nonzero-frequency
metric problem is homogeneous. A finite-time state error instead enters
through R_s c_0 in (3); it is not a freely prescribed exponential stress.
The no-growing statement for the free covariance alone does not exclude
a coupled gravitational mode or establish a physical boundary domain.

## 3. The local metric initial polynomial

In the optical chart use X=(q,w)=(rp,rv). The retained local action and
its fixed contact terms have an even time polynomial

    P(partial_t)=P_0+P_2 partial_t^2+P_4 partial_t^4.        (6)

P_2 and P_4 are the original radial operators obtained from the
[local frequency polynomial](local-coupled-frequency-polynomial.md),
conjugated by the fixed r factors, together with the prescribed canonical
local restoration in the same scheme. These coefficients are not retuned.
The remaining canonical response is a causal memory operator. Its metric
part belongs on the left of the equation; (3) is its independent source.

For initial jets X_j=partial_t^j X(0), Laplace transformation gives

    Laplace[P(partial_t)X]=P(s)X_hat-J_metric(s),
    J_metric=P_2(s X_0+X_1)
             +P_4(s^3 X_0+s^2 X_1+s X_2+X_3).             (7)

In the equation convention P X+F_canonical=J_external, the independent
right side is J_external+J_metric-J^C. This fixes the state-source sign
without assigning a new sign to the stress tensor. If a full residual
equation is being solved around the off-shell seed, its background
residual also belongs on that right side; a constant residual transforms
with 1/s. It is not an independent quantum state source. The seed residual
and all its density contacts remain nonzero wherever the prior notes
retain them.

Equation (7) is exact on compact spatial tests at a common regulator and
as a local differential identity. Combining it with a continuum memory
kernel requires an admitted initial state and a justified initial-value
domain. It does not declare X_0,...,X_3 independent physical data. The
gravitational constraints, gauge, endpoint traces and state matching
still restrict them. In particular third metric jets cannot simply be
varied freely to manufacture an incompatible physical source.

## 4. The affine curvature reconstruction has a range correction

For gamma=s real and admitted by the retained auxiliary inverse, the exact
curvature identity gives

    q_hat=R_gamma U_hat+Q_gamma w_hat+I_gamma,
    I_gamma=(K+gamma^2)^-1 J_kin,
    J_kin=gamma q_0+q_1+(gamma w_0+w_1)/3.                 (8)

To derive (8), transform
U=6F^-1(q_tt+Kq)+(2/F)w_tt+(the retained spatial shape row).
The w initial jets are essential. This agrees with the initial terms
in [physical reconstruction](physical-curvature-reconstruction.md).

If q_0,q_1,w_0,w_1 are smooth and compact in the open exterior, define

    U_I=6 F^-1 J_kin.                                     (9)

It is compact, belongs to D_min, and R_gamma U_I=I_gamma by the exact
positive reconstruction inverse. Thus the affine metric offset in (8)
is the homogeneous reconstruction of the admissible pair (U_I,0).

Write L_0=Atilde+gamma^2+C_gamma-epsilon^-1 P R_gamma,
so L_gamma=L_0+epsilon^-1 D_gamma K_cross^dagger. Here P is the retained
residual multiplier, not (6). The complete auxiliary rows of this offset
have strong trace source f_I=L_0 U_I and shape
Riesz source w_ellI=K_cross^dagger U_I. After shape elimination their
combined source is exactly

    f_I+epsilon^-1 D_gamma w_ellI=L_gamma U_I.             (10)

Every z in ker L_gamma^* annihilates (10). Moving the offset to the right
therefore leaves the full cokernel overlaps unchanged. It must still be
subtracted when reconstructing U and the metric. This proof uses the
full strong/form pair, not separate unproved bounds for every original
offset row. It also explains why dropping the offset from the solution
and claiming source cancellation are different operations.

Apply the [original-source test](ward-sources-and-fast-compatibility.md)
to the normalized original right side from (7), whenever it lies in the
stated source spaces. Compact initial metric reconstruction needs no
additional overlap beyond that test. For noncompact endpoint data,
U_I need not be in D_min; boundary trace terms may survive, and (10)
cannot then be annihilated without checking the domain.

## 5. Initial Ward terms prevent a false thermal identification

For an independent conserved spherical stress on the fixed metric, put
W=N^2 r^2, E=W T^x_x, J=W T^x_t, o=T^t_t+T^x_x and tau=Tr T.
Writing n_0=log N and l=log r, the finite-time transform satisfies

    gamma E_hat-J_hat'=gamma W o_hat+E(0)-W o(0),
    gamma J_hat-E_hat'=-W[(n_0'-l')o_hat+l' tau_hat]+J(0). (11)

These follow by transforming both time derivatives in the exact Ward
system. They differ from the pure-exponential equations used in the
preceding completion theorem. No reinterpretation of E(0),J(0) as free
thermal bulk sources is justified by conservation alone.

There is a useful exact counterexample to imposing thermal bolt behavior
on a general finite-time transform. Choose a nonnegative nonzero smooth
b supported in -2<v<-1 and set, for t>=0,

    T=b(t+x)(dt+dx)^2/(2r^2), E=J=b(t+x)/2, tau=o=0.

This conserved ingoing pulse has zero outer flux for t>=0 and is smooth
at the Lorentzian horizon: in V=exp(kappa_h v) its support stays away
from V=0. In a collar x<-2, however, its Laplace transform is

    E_hat=J_hat=exp(gamma x) A_gamma/2,
    A_gamma=integral_-2^-1 exp(-gamma v)b(v)dv>0.           (12)

At gamma=kappa_h, interpreting that radial transform as a Euclidean
first harmonic would produce the V^-1 tensor coefficient. The original
Lorentzian pulse is regular; its finite-time transform is not the
equilibrium first-harmonic response. This stress example makes no claim
of realization by the Dirac state. It proves that conservation and
Lorentzian smoothness alone cannot justify the identification. The
[first-thermal Dirac theorem](first-thermal-horizon-flux.md) retains its
fixed initial-state, exponentially prepared source hypotheses.

## Decision and verification scope

For the selected fixed-HH exponential preflight, analyze the homogeneous
coupled metric operator and its physically justified boundary domain;
there is no free independent growing covariance to tune. For a finite-time
disturbance, use (2)-(3), (7) and (11), impose the actual initial constraints
and state regularity, then evaluate the original-source overlaps. The
compact reconstruction offset (8) cannot supply an extra cancellation.
No actual adjoint profile, physical initial-data image, nonlinear solution
or all-sector stability theorem is evaluated here. S075 and S077 physical
preparation/disturbance gates remain open.

Focused verification checks nonzero initial metric amplitudes and both
Dirac frequency signs against direct time integrals; the local initial
polynomial; the affine range cancellation; and the initial Ward terms
against a translated pulse. Finite matrices and pulse controls do not
replace the charged continuum kernel. No Lean campaign is selected.

Remaining to-do list: constrained initial metric/state image and its
source overlaps; homogeneous coupled operator with physical gravitational
boundary domain; continuous causal estimates, static signed condition,
nonlinear existence and all-sector stability; exterior/state, charge/return,
opening, transport, preparation, separate creation and full audit.
