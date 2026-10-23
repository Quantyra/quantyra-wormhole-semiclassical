# Cold source calculation and the full solution feasibility audit

2026-10-22, E003/S169. Ordinary argument before frozen verification.
The full objective remains an actual stable semiclassical wormhole.

Decision: retain the specified c_R=10^22 construction and pursue a direct
enclosure of its coupled physical (b,V) correction. Do not extend the
generic higher-inverse chain. The decisive source-loss calculation below
closes a useful cold canonical C3 input, with an allowance 5e17 q^2 a_s,
or 1.36e-9 in canonical source units. The full metric inverse, physical
neighborhood, independent sources, validity and stability remain unproved.
This is a source estimate and research decision, not a controlled solution.

## 1. Complete requirements and available constants

Keep G, flux, coupling, mu0, temperature, paired HH realization and all
finite coefficients of S169 fixed. This is the explicitly changed R^2
theory, not the old zero-additional-c theory. Let X control physical cold
values/gradients and raw C5/cap jets, the middle derivatives, outer weighted
remainder, correlated mass and fixed-temperature return. These norms must
be included with their original physical units; merely renaming an energy
norm X does not establish the required embeddings.

For the complete physical residual F, valid fixed projection Pi, and
preconditioner B on its actual domain, put

    Y0=||B Pi F(0)||_X,
    Z0=||I-B Pi F'(0)||_(X->X),
    ||B Pi [F'(u)-F'(v)]||<=Z2 ||u-v||_X.

The last inequality must hold on EVERY real metric in the radius-rho ball,
with a common physical domain and F continuously differentiable there.
The exact sufficient preservation and contraction inequalities are

    Y0+Z0 rho+(Z2/2)rho^2<=rho,
    Z0+Z2 rho<1.                                      (1)

B must be injective on the projected range. A projected solution is a full
solution only when (I-Pi)F=0 exactly. Use Pi=I for S169 unless its actual
range condition is derived; the old signed projection is not inherited.
At rho=1e-8, the design allocation Y0<=5e-9, Z0<=1/4, Z2<=5e7 gives
preservation 1e-8 and contraction <=3/4. These are targets, not measured
constants. In general the missing Z2 must meet BOTH

    Z2<=2[(1-Z0)rho-Y0]/rho^2,
    Z2<(1-Z0)/rho.                                    (2)

For an approximate center, a certified residual delta in preconditioned
X units gives error <=delta/(1-q), within an admitted invariant ball.
For q<=3/4 and desired error 5e-9, require delta<=1.25e-9. An unbounded
angular, spatial or boundary tail cannot be omitted from that residual.

The numerical admission and block targets are:

| Requirement | Exact target or sufficient allocation at rho=1e-8 | Status |
| --- | --- | --- |
| Every coefficient margin | b_i+E_i rho<=a_i; require positive physical slack where the theorem uses an open ball | Complete b_i,E_i unknown |
| Interior P45 | b_45+1e-8 E_45<=.001; if b_45<=1.3e-4 transfers, E_45<=87000 | Old path bound, conditional transfer |
| Cold intersected contract | b_c+1e-8 E_c<=5e-7; if E_c<=1, b_c<=4.9e-7 | Unknown for full new correction |
| Outer contract | b_o+1e-8 E_o<=1e-7; if E_o<=1, b_o<=9e-8 | Old upper bound 1.01e-7 fails to certify slack; not a lower bound |
| Actual inverse/domain transfer | ||L_ref^-1(L_b-L_ref)||<=.25 in the same full graph after physical identification | Unknown; no division by the tiny thermal gap establishes it |
| Scalar/metric graph block | PC<1 and (M0+P A)/(1-PC)<=5e-9 | A,C,P,M0 unknown in derivative-bearing spaces |
| Sufficient block allocation | PC<=.25, M0<=3e-9, P A<=7.5e-10 | Scalar VALUE constants are not A,C |
| Radius and bolt | r>0, 1+s partial_t V>0, positive bolt coefficient | Must hold on the whole ball |
| Curvature recovery | |phi+12lambda|<=lambda, |R-phi|<=lambda | Algebraic recovery bracket >1-2e-8; its ball embeddings remain needed |
| Boundary and matching | paired adjoints, fixed clock/period, regular bolt, fast phase, correlated mass/return | Exact conditions, not small residual allowances |

For example the missing graph feedback P must satisfy
P<=[5e-9-M0]/[A+5e-9 C], with positive numerator and PC<1.
There is no numerical value to substitute for an unknown complete norm.

## 2. Substitution ledger, without mixing scopes

The latest [middle source result](fixed-cap-seed-source-three.md) is
6e15 q^2 in j!100^j units, on s>=.001,Y>=.05. The canonical conversion is

    eta_G<8.5e-33, q/C<400,
    S_can=2G Theta/(B^2 C),
    2 eta_G(q/C)^2<2.72e-27.                           (3)

Thus the middle allowance is 1.632e-11. Section 4 proves the cold allowance
5e17 q^2 a_s, with a_s=s[1+log(.01/s)], hence 1.36e-9 a_s. These are
seed canonical scalar bounds, not bounds on all four-dimensional metric
sources or Fréchet derivatives over a nonlinear ball.

Let K_m,K_c be norms of the COMPLETE physical preconditioner and return
applied to the respective regional source extensions, including partition
commutators. Let Y_rest include outer, independent Weyl/radius, other
couplings and all residual/tail errors. A sufficient source budget is

    1.632e-11 K_m+1.36e-9 K_c+Y_rest<=5e-9.             (4)

For a specified K_c,Y_rest, the missing K_m must not exceed
[5e-9-1.36e-9 K_c-Y_rest]/1.632e-11, and likewise solve (4) for K_c.
One explicit allocation is K_m<=61.2745098039, K_c<=1.4705882352,
Y_rest<=2e-9 (the displayed downward decimals allocate 1e-9,2e-9,2e-9).
These operator bounds are UNKNOWN. If each source used the whole budget
alone, its thresholds would be 306.3725490196 and 3.6764705882.
If K_m=K_c=1 could be proved, the two sources would leave 3.62368e-9.
This conditional arithmetic is not an assumed unit inverse.

The retained outer Euler C3 allowance is 1e20 q^2 L(Y), for Theta-S,
where S=c C^2/(5pi Y) must be restored with its Weyl partner. It gives
2.72e-7 in those Euler source units and a whole-budget return target
K_o<=.0183823529411. Ordinary weighted derivatives use the larger
2e20 allowance and target .0091911764705. This is an old seed estimate,
not yet the completed product-centered physical source/domain transfer.
Do not silently bound S alone or drop it. Its physical cancellation and
the independent Weyl/radius sources belong in Y_rest before (4) is used.

Keep the latest complete OLD interior and outer allowances visible:

| Quantity in original zero-additional-c proof | Available upper allowance |
| --- | ---: |
| C_I | 9.1591500000000000015e103 |
| C_O | 2.365283519093731664e186 |
| C_tau | 6.87674442436215072e19 |
| York C_sigma, C_b4 | 1421.40132786776365; 4.52219350981899516e14 |
| C_F,I=C_sigma+C_b4 C_I | 4.141944868545865e118 |
| C_F=C_sigma+C_b4(C_I+C_O) | 1.069626977892751e201 |
| Cstart, mhi | 3.123143830205975e28; 2.577600186401461e12 |

With g+(11/(360pi))g^2=2000^-10 and g approximately 9.765625e-34,
the old certificate requires

    R[Cstart+mhi C_D4 C_F/g]<=1e57 g epsilon_c/800.     (5)

Even at epsilon_c=1 its low term requires R<3.90858e-8, whereas the
proved lower bound is R>.97. Removing that allowance entirely still
requires C_F<=2.38393e-25 at R=.97,C_D4=2. The interior and global upper
allowances miss this sufficient threshold by >1.7374e143 and >4.4868e225.
At epsilon_c=1e-4 all thresholds tighten by 1e-4. The old signed target
is |beta|<=1e-49 epsilon_c; it is not a new S169 compatibility condition.
The unknown inverse UPPER constants cannot rescue (5), given its lower
bounds. These are failed certificates, not measured source amplitudes or
nonexistence results. The new middle/cold C3 numbers do not replace C_I:
different source, derivatives, geometry and action are involved.

The actual auxiliary scalar VALUE result is

    ||v||_w<=4.71875e-37+3.125e-36 N2(h_E),
    w=[(s+h)/t]^(2/5), h_min=1.220703125e-34.

With W_c<4.7e18 and W_r<3e20, its scalar-only contribution to the radius
is <1.1e-15 for N2<=1. The exact physical return, however, is

    tau=s/t, zeta=ell_E-v,
    V=log[1+(exp(zeta)-1)/tau],
    delta b=(a-c)/2, delta V=(c/2-v)/tau.               (6)

The auxiliary scalar cancels from b exactly. A generic VALUE certificate
using the sharper W_r=2.940667788792409618e20 requires N2<=approximately
3.400588043500485e-29 to achieve 5e-9; N2<=3.3e-29 suffices. Finite return
requires sup|zeta|/tau<=1-exp(-rho). This severe generic loss is avoidable
only by estimating the actual joint combination c/2-v and its derivatives.
On a compact annulus b_n=n^-2 chi sin(ny) has bounded N2 and unbounded
raw C5 proportional to n^3. Therefore NO finite N2-to-C5 embedding can
close admission. These test fields are not solutions or unstable modes.
The reference energy 1.44e-10 would need a cold embedding <=34.7222,
against its proved lower bound >2.6e15; a .25 transfer tightens that to
26.0417. Reference TT lower bound >.994 and mixed relative bound <.336
remain useful in their reference energy space, not in the missing X norm.

## 3. Decisive calculation: can the cold source loss be closed?

The best retained complete cold C3 allowance was 1e20 q^2 a_s. Blindly
restricting the later full C8 proof to j<=3 actually gives a WORSE
3.816942631332431e20 allowance. The latter ledger, however, keeps every
coordinate, half-density, spin, covariant and actual-domain contribution
explicit. We use that ledger and sharpen its two local rows. The fixed
units, physical state and subtraction are unchanged.

Use exactly the moving normal cap of
[the local theorem](seed-cold-local-eight-jets.md) and
[the gluing theorem](seed-cold-gluing-eight-jets.md), with the corrected
Hermitian Clifford convention of
[the actual-domain theorem](seed-cold-actual-domain-eight-jets.md).
The cap multiplies exp(sigma_N)-1 by the fixed radial cutoff chi, equal
to one through 3/16 and zero past 1/4. Both gluing rows include ONE
commutator, supported on 1/8<=r<=5/32, bounded by 33.

The Goursat equations hold uniformly on the independent parameter
polydisk |u0-1|,|p|,|pbar|<=.01, including the retained thermal auxiliary
disk. They give |sigma_N|<.025, |Bcal''|<2.31 and zero axis values.
Integrating sigma_zw=Bcal'' exp(2sigma)/8 twice along straight segments
to (z,bar z) gives, also for complex parameter values,

    |sigma_N(y)|<=2.31 exp(.05) r^2/8,
    |h(y)|=|chi(exp(sigma_N)-1)|<.32 r^2,
    |h|<.026.                                        (7)

Indeed 2.31 exp(.075)/8<.32. Outside the quarter disk h=0. This keeps
the exact quadratic vanishing at the normal center; it does not assert
that the global cold cap is O(a_s) without a distance weight.

Let R=(D_flat+M(1+h))^-1, S=(D_flat+M)^-1 and theta=M/2, M>18.
Truncated Lipschitz distance weights and bounded perturbation of the free
resolvent give, for both signs of the conjugation,

    ||exp(-theta r) R exp(theta r)||
       <=1/[M(1-.026-.5)]<2.11/M.                    (8)

This is valid for the COMPLEX parameter family: Re(1+h)>=1-|h|, and
the perturbation plus weight has norm <.526M. The adjoint has the same
coercive bound and the range is all L2. The incoming row is estimated
separately; no conjugate-parameter differentiation is used.

The free rank-two kernel has HS square
M^2[K0(Mr)^2+K1(Mr)^2]/(2pi^2). The positive-real integral formulas
[DLMF 10.32.8-9](https://dlmf.nist.gov/10.32) imply
K0(t)<=sqrt(pi/(2t))exp(-t) and
K1(t)<=sqrt(pi/(2t))exp(-t)(1+3/(8t)).
The first follows from cosh v>=1+v^2/2; the second from
sqrt(1+v/2)<=1+v/4 in the nu=1 integral. These are free-kernel bounds,
not a theorem about the charged curved Green function.

For d=M/8>=2.25 the remote weighted free row obeys

    ||1_(r>=1/8) S exp(theta r)||_HS^2
      <=[2+3/(4d)+9/(64d^2)] exp(-d)/2 < .354^2.

Also ||r^2 S exp(theta r)||_HS<5.3/M^2; its squared dimensionless
moment is <=[2*4!+(3/4)*3!+(9/64)*2!]/2<5.3^2.
Using R=S-S M h R before evaluating the row yields

    ||R(0,.) exp(theta r)||_(r>=1/8)
      <.354+2.11(.32)(5.3)/18^2.                     (9)

The point singularity is harmless: the free term is remote, and the
other term contains h=O(r^2). On the commutator annulus the frozen
weight W^-1 is <=exp[2(1.06)(5/32)]<1.4. In the canonical flat-volume
representation Rlocal=exp(sigma_N/2)R exp(sigma_N/2), its observation
factor is one and its other factor has modulus <exp(.025/2). The
commutator already contains its exp(-sigma_N) and slope; 33 includes
them. No additional coordinate or spin endpoint appears: the transported
global operator and local cap coincide there and its center density is
exactly one. Reserve 1.1 for the remaining factor (larger than exp(.025/2)).
Consequently each commutator-weighted row satisfies

    ||W^-1 rL||, ||rR W^-1|| <=12 exp(-M/32),          (10)

because its coefficient is at most

    33(1.4)(1.1) exp(-18/32)
          [.354+2.11(.32)(5.3)/18^2] <12.

We split exp(-M/16) into two exp(-M/32) factors when removing the
remote exponential weight; the remaining coefficient decreases with M.
The proof is uniform on the ORIGINAL analytic normal-parameter polydisk.
It therefore supports the same row Cauchy and cold-difference estimates,
with 400 replaced by 12. It does not import the middle plain-row constant
ten: that row has different cutoffs and does not include a commutator.

## 4. Full low-order cold ledger and numerical consequence

All central inverse and coordinate constants are unchanged:
d=6/5, G=12, u=20, w=25000, k=15000,
D0=360000, D_j=360000*241^(j-1), j>=1. With epsilon=.01,

    T0=362400,
    Tj=sum_(b=1)^j (j-b+1)epsilon^(j-b)D_b
                         +1200(j+1)epsilon^j,
    H_j=4505.6(12/400)^2 150^j Tj.                   (11)

These retain both commutators, spin rank two, all ordered inverse
insertions, one cold weighted vertex and the complete angular majorant
sum d_n exp(-M_n/16)<4q. No inverse or normalization was changed.
For covariant tensors use exactly the existing nonnormal connection:
gamma_a=20 a!300^a, B(0,k)=k!100^k H_k,

    B(r+1,k)=B(r,k+1)+r sum_(a=0)^k binom(k,a)gamma_a B(r,k-a),
    J_j=(6/5) B(j,0)/(j!100^j).                       (12)

The nonzero gamma_0 is retained. The actual-domain error is separately
estimated in its own fixed normal chart by the retained formulas

    L_j=6*18^j*(6e-8)[1+2/(16.776-j-3)],
    H^act_j=1e19(j+1)(.00352)(1000/3)L_j.              (13)

Apply (12) with frame factor one, gamma_0=0 and
gamma_a=2 a!25^a for a>=1. This uses the actual paired-domain inverse,
the seam distance .94 log(.01/s)-.3, M_n>=18sqrt(n) and the whole angular
tail. Different charts may be combined AFTER conversion to physical
tensors, not by equating their ordinary derivatives.

The local moving-cap bound is 1e8 in the same units through eight.
Reserve a further 1e8 for the canonical LLL trace if separated from that
local ledger: its anomaly is bounded by Cq times the curvature variation.
Here R_gamma=-Bcal'', |R_gamma+2|<12a_s on the observation bidisk .01;
Cauchy and (12) through three give <1e4 a_s in normalized physical jets,
and C/q<.00352 makes this allowance <36. The LLL traceless/state response
is an independent physical source, not inferred from this trace bound.

Exact rational evaluation gives the following rounded values:

| j | Gluing with old row 400 | Actual-domain error | New gluing, row 12 |
| --- | ---: | ---: | ---: |
| 0 | 1.95940e9 | 4.83725e12 | 1.76346e6 |
| 1 | 2.91983e11 | 1.75869e14 | 2.62785e8 |
| 2 | 1.05554e16 | 4.80304e15 | 9.49983e12 |
| 3 | 3.81578e20 | 1.16826e17 | 3.43420e17 |

In each row, add 2e8 for the local/LLL reserves. Their maximum is below
4.603e17, hence below 5e17. This proves the stated complete canonical
cold C3 allowance on the actual product-centered SEED. The degree-six
polynomial, thermal and scale guards are the ones checked in the middle
audit; |d|/C<1e-24 and |log(mu0 B)|<1e-8 preserve the strict cold guards.
No finite corrected metric or whole neighborhood is asserted analytic.

The executed feasibility test is informative even before using (10).
Keeping every other constant fixed, the old ledger would need the
commutator-weighted row coefficient <=18.34149 to allocate 2.5e-9 to the
cold source at unit return norm. Alternatively its derivative multiplier
k=15000 would need to fall below approximately 1921.26475. The proof above
achieves row 12 with no change to k. It closes that particular source loss.
With the original row, the normalized actual-domain part alone is
3.17766554814e-10; reducing only local gluing can never remove this
ALLOWANCE. It is not a lower bound on the true domain error. The old
1e20 C3 result remains valid but is weaker than the new bound by 200.

## 5. Hard continuum and physical steps; order of work

The largest source-specific loss left here is the transported inverse
derivative factor 241^(j-1) k^j, followed by the actual-domain row
derivative allowance. More source sharpening may help after the physical
operator is quantified, but neither loss is currently a demonstrated
obstacle to a small actual solution. The generic energy-to-radius loss
and absent C5 implication are earlier obstacles to the proposed solver.

The quantum hypotheses must be checked at their actual scope:

* Spatial seed C3 estimates do not supply first and second METRIC responses.
  The full traceless logarithmic order-four operator, angular tower, common
  subtraction, LLL, Weyl terms, densities and contacts remain required.
  Finite-part cancellation must precede absolute estimates. A scalar trace
  Green function does not determine this tensor response.
* Both partner domains and all endpoint factors remain in the gluing
  identities. Complex auxiliary caps are proof devices; the physical
  neighborhood needs real-metric uniform graph estimates. The cutoff is
  not spatially analytic. The middle proof uses one fixed cap; the present
  cold proof uses the same moving cap for its local and gluing terms.
* The LLL zero-mode argument on the selected disk proves kernel absence,
  not a quantified gap or radius of domain stability. For a fixed domain
  identification a sufficient perturbation target is
  ||delta D_LLL R_LLL,0||<=1/2. If bounded by K_domain rho, it needs
  K_domain<=5e7 at rho=1e-8. K_domain and the base inverse remain unknown
  numerically; no small fixed temperature may be replaced by zero.
* Reference form closedness and normal symbols do not establish the
  actual physical graph on the product-centered seed and every finite
  metric in the ball. The reference kernel construction has a finite but
  unquantified end constant in C h^2(1+|log h|)^N. For a chosen cutoff,
  the remaining transfer allocation d_end requires
  C<=d_end/[h^2(1+|log h|)^N]. Its C,N and appropriate cutoff are not
  numerically supplied. A general boundary logarithm theorem assumes an
  already specified elliptic realization and spectral cut; it does not
  manufacture this determinant-Hessian domain
  ([Gaarde--Grubb](https://arxiv.org/abs/math/0703878)).

Compatibility is not presently a proved obstruction, but belongs BEFORE
further inverse iteration. On the chosen fast phase the stationary mass
momentum and its tangent vanish exactly; the current/reflection and Ward
identities are retained. The scalar curvature-recovery bracket is positive.
None of these authorizes deleting an arbitrary residual or imposing an
extra scalar phase condition. Any remaining cokernel requires an exact
signed equation and a justified physical parameter with nonzero derivative.
The specified temperature and c_R are fixed, not free cancellation knobs.

Likewise |delta B^2/B^2|<5.5e-23, orbit shift <1.71e-29 and horizon
residual <2e-12 are seed/value facts, not finite-domain transfer norms.
The exact cancellation of v from b and joint c/2-v return are promising
because they hold before estimation. Product centering and the correlated
outer singular/Weyl terms are similarly real cancellations. Estimating
them separately can destroy the thermal small factor. Changing units
alone rescales Y0 and rho together and changes no dimensionless test.

Small G lambda<1e-33, positive scalar m^2/lambda approximately 3.3953e8
and F_*>1-1.179e-8 do not prove semiclassical validity or retarded stability.
For a specified smearing and physical tolerance, omitted deterministic
terms must fit their assigned preconditioned residual budget; induced
fluctuations require their own noise-kernel/retarded-response bound.
No noise allowance, omission bound or all-sector retarded spectral margin
has been calculated here. Assigning them zero would be false. The need
to retain induced as well as intrinsic fluctuations is explained by
[Hu, Roura and Verdaguer](https://arxiv.org/abs/gr-qc/0508010).
Euclidean positivity alone does not supply that Lorentzian criterion.

## 6. Research decision and acceptance result

Pursue the sharper method: a source-specific, product-centered nonlocal
physical (b,V) graph solve with the exact auxiliary constraint and joint
return, actual fast boundary/mass phase and both paired matter domains.
First quantify the linear physical realization and residual, including
outer and independent Weyl/radius components, against (4). Then prove
open-ball transfer and all admission margins before using (1). The current
source calculation establishes why this route is still quantitatively
worth testing; it does not settle feasibility of the full construction.

The decisive calculation selected and executed in this audit was whether
the cold row loss can meet a useful source budget while retaining its
vanishing thermal factor and true domain error. Equations (7)-(13) answer
yes. The decisive calculation for EXISTENCE is still an enclosure of the
actual coupled correction and its domain, with certified continuum tails;
that calculation is not completed or replaced by the supporting checker.
The missing boundary/source operators are explicitly listed above. Do not
claim the full goal complete on the strength of this smaller source bound.

Switch construction only upon a documented incompatibility, physical
obstruction, or failure of the actual source-specific method with a reason;
the oversized old upper bounds supply none of those. Preserve S169.2-S169.5
unchecked and S166 as the last completed component. Opening, transport,
preparation and separate creation remain full-objective requirements.

Remaining to-do list: actual coupled physical graph/domain and correction;
outer and complete Weyl/radius sources; positive admission and nonlinear
preservation; exact matching, semiclassical validity, all-sector stability,
charge/return, opening, transport, preparation, separate creation and audit.
