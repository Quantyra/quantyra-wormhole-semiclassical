# Coupled high-frequency boundary phase and conditional winding

2026-10-21, E003/S168. Ordinary argument before focused verification.
The stable-wormhole goal remains incomplete.

At each sufficiently high admitted thermal frequency, the full reduced
auxiliary curvature operator has exactly one homogeneous mode on its
maximal graph domain. Exactly one real reflecting phase admits this mode.
Every other separated real phase gives a bijection on the strong source
space. The resonant phase differs from the inverse-square Bessel phase
by O(gamma^-2). This proves a frequency-dependent boundary alternative
for the retained coupled operator, without adopting a physical phase.

Along the Matsubara sequence the resonant phases approach every phase
arbitrarily closely. This does NOT imply an exact resonance at a fixed
phase and a Matsubara frequency. A continuous-frequency resonance theorem
is proved below under an additional real, continuous, uniformly bounded
causal-family hypothesis. That hypothesis is not established for the
actual full shape response. No charged Lorentzian instability is claimed.

## 1. Maximal operator and its single homogeneous vector

Use the actual bounded remainder from the preceding
[high-frequency theorem](high-frequency-zero-fast-injectivity.md):

    A_min,gamma=-partial_y^2-a/y^2+gamma^2 on D_min,
    A_max,gamma=A_min,gamma^*, H=L2(0,infinity),
    L_min,gamma=A_min,gamma+V_gamma,
    L_max,gamma=A_max,gamma+V_gamma,
    ||V_gamma||<=M_*, eta_gamma=C_a M_*/gamma^2<1.          (1)

The maximal domain has both outer traces free and the limit-point
horizon; D_min sets both traces to zero. Both domains are independent
of gamma. The remainder retains the seed potential, reconstruction,
full shape inverse, residual contacts and nonlocal strong feedback.
L_max is not L_min^* when V_gamma is not self-adjoint.

Let b_gamma be the normalized real inverse-square horizon solution,
P_gamma=I-|b_gamma><b_gamma|, and T_gamma the previous compatible
minimal inverse, composed with P_gamma. Its adjoint Q_gamma=T_gamma^*
is used only in this note; it is NOT the physical shape reconstruction
map bearing the same letter in earlier notes. We have

    ||Q_gamma||<=C_a/gamma^2,
    A_max,gamma Q_gamma=I,
    Q_gamma A_max,gamma=P_gamma on D_max,
    Ran Q_gamma subset D_max intersect b_gamma^perp.      (2)

Define

    h_gamma=(I+Q_gamma V_gamma)^-1 b_gamma,
    J_max,gamma=(I+Q_gamma V_gamma)^-1 Q_gamma.            (3)

These Neumann inverses exist by (1). The equations
h=b-Q V h and J_max f=Q(f-V J_max f) show their outputs lie in D_max.
Applying A_max and using (2) gives

    ker L_max,gamma=span{h_gamma}, L_max,gamma J_max,gamma=I,
    <b_gamma,h_gamma>=1, ||h_gamma||<=1/(1-eta_gamma),
    ||h_gamma-b_gamma||<=eta_gamma/(1-eta_gamma),
    ||J_max,gamma||<=C_a/[gamma^2(1-eta_gamma)].           (4)

For completeness, every homogeneous u satisfies
(I+Q V)u=b<b,u>, by applying (2). This proves the whole kernel assertion,
not just construction of one mode. The corresponding homogeneous shape
is w=-K_cross^dagger h. No extra free state source is inserted.

The earlier source-adjoint vector is instead

    z_gamma=(I+Q_gamma V_gamma^*)^-1 b_gamma.             (5)

It annihilates Ran L_min. It generally differs from h_gamma. Interchanging
V and V^* would solve the wrong problem for the nonself-adjoint strong row.

All maps in (3) preserve real radial amplitudes at admitted thermal gamma.
The [thermal spectral representation](thermal-shape-response-and-subtractions.md)
has a real radial Hamiltonian/projector and reciprocal real mass vertices.
Its polarized nonzero-frequency form and its local subtractions are real;
the real full shape form has a conjugation-preserving Riesz inverse.
The reconstruction and mixed maps have real coefficients. Finally the
pure compatible inverse is real by uniqueness on its real compatible
range, despite its representation using complex I_(ik). Thus h is real.

## 2. A graph estimate controls the boundary trace

Fix y0>0 and k=sqrt(a-1/4)>0. Write

    u(y)=sqrt(y)[A_c cos(k log(y/y0))+A_s sin(k log(y/y0))]
                                      +o(y^(3/2)),
    tau(u)=(A_c,A_s).                                   (6)

The derivative remainder is o(y^(1/2)), as in the retained
[metric trace theorem](coupled-fast-metric-boundary-traces.md).
H-norm proximity alone does not control tau. Instead choose a smooth
compact cutoff chi equal to 1 near zero and define v_c and v_s by
multiplying the two leading basis functions in (6) by chi. They belong
to D_max and A_max,1 v_c,A_max,1 v_s are in H. Green's identity is

    <u,A_max,1 v>-<A_max,1 u,v>
             =k[conj(A_c) B_s-conj(A_s) B_c].            (7)

Testing with v_s and v_c gives a finite constant C_tau such that

    ||tau(u)|| <= C_tau (||u||+||A_max,1 u||).            (8)

For example take C_tau=k^-1 times the Euclidean norm of
(||A_max,1 v_s||+||v_s||, ||A_max,1 v_c||+||v_c||).
No actual benchmark value of this constant is evaluated.

Unitary dilation S_gamma u(y)=sqrt(gamma)u(gamma y) multiplies tau
by gamma and rotates it by -k log gamma. Put e=h_gamma-b_gamma and
e_1=S_gamma^-1 e. Equations (1)-(4) imply

    ||e_1|| <= C_a M_*/[gamma^2(1-eta_gamma)],
    ||A_max,1 e_1|| <= M_*/[gamma^2(1-eta_gamma)].

The second bound uses A_max,gamma e=-V_gamma h_gamma and exact dilation;
it does not differentiate the bounded nonlocal remainder. Hence

    ||tau(h_gamma)-tau(b_gamma)||
       <= C_tau(C_a+1) M_*/[gamma(1-eta_gamma)].          (9)

Let N_k=||sqrt(y)K_(ik)(y)||, B_k=|Gamma(ik)|/N_k>0. The
[I-series](https://dlmf.nist.gov/10.25) and
[K connection formula](https://dlmf.nist.gov/10.27) give

    K_(ik)(z)=Re[Gamma(ik)(z/2)^(-ik)]+O(z^2).

This is also the [imaginary-order small-argument formula](https://dlmf.nist.gov/10.45),
with Gamma(1+ik)=ik Gamma(ik). The convergent series justifies its
differentiated remainder. Therefore

    b_gamma(y)=gamma sqrt(y)K_(ik)(gamma y)/N_k,
    tau(b_gamma)=gamma B_k (cos beta_0,gamma,sin beta_0,gamma),
    beta_0,gamma=arg Gamma(ik)-k log(gamma y0/2).          (10)

Set delta_gamma=C_tau(C_a+1)M_*/[B_k gamma^2(1-eta_gamma)].
Above a further finite threshold delta_gamma<1/2. A disk of radius
delta_gamma around a unit vector subtends angle at most arcsin(delta_gamma).
Consequently the real trace angle has the uniquely nearby lift

    beta_gamma=beta_0,gamma+e_gamma,
    |e_gamma|<=arcsin(delta_gamma)<=2 delta_gamma
                                              =O(gamma^-2). (11)

In particular tau(h) is nonzero. Even below this extra trace threshold,
where (1) holds, a zero trace would put h in the already injective D_min
and contradict <b,h>=1.

## 3. Exact alternative for each real reflecting phase

The [boundary symplectic calculation](fast-curvature-boundary-symplectic.md)
identifies separated maximal isotropic real trace lines as

    ell_theta tau(u)=cos(theta) A_c+sin(theta) A_s=0,
    theta modulo pi; D_theta subset D_max.               (12)

The UNIQUE resonant line at each gamma in (1) is

    theta_res(gamma)=beta_gamma+pi/2 modulo pi.           (13)

If theta differs from theta_res, the operator on D_theta is bijective,
with explicit inverse

    J_theta f=J_max f-h [ell_theta tau(J_max f)]/
                                      [ell_theta tau(h)]. (14)

Indeed it solves the full row and its trace lies in D_theta; uniqueness
uses the complete kernel statement (4). At resonance, ker L_theta is
span{h}, and

    Ran L_theta=Ran L_min={f:<z_gamma,f>=0}.              (15)

To see (15), subtract a multiple of h from any u in D_theta to remove
its entire trace; the result is in D_min and has the same image. The
reverse inclusion is immediate. The resonant cokernel is exactly one
dimensional. No source is silently projected onto that range.

One useful quantitative consequence of (8) and dilation is

    ||tau(J_max f)||<=C_tau(C_a+1)||f||/[gamma(1-eta_gamma)].

Here A_max,gamma J_max f=f-V J_max f and eta=C_a M_*/gamma^2.
Combining this with (9)-(14), above delta_gamma<1/2, gives

    ||J_theta|| <= C_a/[gamma^2(1-eta_gamma)]
       + C_tau(C_a+1)/[gamma^2 B_k(1-eta_gamma)^2
           (1-delta_gamma)|sin(theta_res-theta)|].       (16)

The small denominator is retained. This is not a positive energy bound
or a physical stability estimate.

## 4. Discrete winding and a conditional continuous-frequency theorem

At admitted gamma_n=n kappa_h, (11)-(13) give a lifted phase

    theta_res,n=arg Gamma(ik)+pi/2
                   -k log(n kappa_h y0/2)+O(n^-2).       (17)

It tends to minus infinity and eventually decreases: the main increment
is -k log(1+1/n), while the error increment is O(n^-2). For any target
theta, choose the real solutions N_j of the leading phase equal to
theta-j pi and take the nearest integers n_j. The phase error is O(1/n_j),
so the projective phases approach every target as n_j tends to infinity.
This establishes asymptotic density, NOT equality at any Matsubara n.
A countable dense set can miss a specified phase entirely.

Now impose an ADDITIONAL hypothesis: the actual causal shape elimination
defines V_gamma for every sufficiently large real gamma, on the same
strong curvature space and graph domains, with the bound in (1), real
coefficients and operator-norm continuity in gamma. This is not furnished
by the existing Matsubara coercivity theorem.

Under that hypothesis, h_gamma and its trace are continuous. To justify
the pure inverse continuity, fix gamma_0 and let d=gamma^2-gamma_0^2.
For small d, F_gamma=(I+d Q_gamma_0)^-1 Q_gamma_0 is a norm-continuous
right inverse of A_max,gamma. Since b_gamma is H-continuous under dilation,

    Q_gamma=(I-|b_gamma><b_gamma|)F_gamma

is the unique right inverse with range perpendicular to b_gamma, hence
is norm-continuous. Formula (3) gives H-continuity of h. Its differential
equation gives continuity in the fixed maximal graph norm, and (8) then
gives trace continuity. No norm-continuity of the dilation operator on
all of H is assumed.

For each fixed theta define the pure resonant sequence

    gamma_j^0=(2/y0) exp[(arg Gamma(ik)+pi/2-theta+j pi)/k].

The continuous lifted phase differs from the pure one by at most C/gamma^2
above a fixed threshold. At log-frequency distances of order
2C/[k(gamma_j^0)^2] on opposite sides of each sufficiently large pure
root, the phase residual has opposite signs. The intermediate value
theorem supplies a sequence of actual auxiliary homogeneous zeros with

    gamma_j/gamma_j^0=1+O((gamma_j^0)^-2),
    gamma_(j+1)/gamma_j -> exp(pi/k).                     (18)

When C=0 the roots are exactly the pure ones. No uniqueness, simplicity
or analytic pole assertion follows from mere continuity. The intervals
are eventually disjoint for this fixed k>0. Constants are not uniform
in a benchmark limit k->infinity.

## 5. Research consequence and verification scope

The auxiliary boundary alternative and thermal winding are unconditional
within the retained high-frequency hypotheses. Equation (18) is
conditional on a new actual continuous causal-family theorem. The compact-
probe third divided-difference positivity of the thermal spectral measure
does not by itself prove global endpoint-uniform shape coercivity between
Matsubara points: interpolation coefficients have both signs and the
reconstruction and boundary compensation depend on frequency. Arbitrary
interpolation of discrete operators is not the physical causal response.

Neither a fixed physical phase nor continuous causal zeros have been
obtained for the charged seed. The seed is still off shell. All constraints,
initial/state admissibility and the earlier conditional nonlinear
fast-family obstruction remain relevant. No new boundary action, finite
counterterm, source projection or order reduction is adopted. The actual
C_a, C_tau, M_*, B_k and semiclassical frequency window are uncalibrated.

Focused checks compare Bessel values and derivatives with analytic boundary
coefficients at three cutoffs and two precisions. A constant bounded
potential supplies an exactly solvable continuous comparator: its phase
shift is -(k/2)log(1+v0/gamma^2), and its fixed-phase roots are
sqrt((gamma_j^0)^2-v0). A real nonnormal finite row model checks the
maximal right inverse, boundary correction and resonant range identities.
These controls do not approximate the charged kernel or prove the
ordinary infinite-dimensional theorem by finite matrix rank. No Lean
campaign is selected.

Remaining to-do list: actual continuous causal shape realization and
physical gravitational boundary/initial-state domain; quantitative
thresholds and finite frequencies; static signed condition, nonlinear
existence and all-sector stability; exterior/state, charge/return,
opening, transport, preparation, separate creation and full audit.
