# Uniform admission of the fast graph sector to the mass/log coefficient domain

2026-10-24, E003/S169. Ordinary argument before frozen acceptance.
The full stable-wormhole objective remains ACTIVE and INCOMPLETE; audit CLOSED.

## 1. An admission requirement before the matched operator estimate

The geometric heat-action argument is C2 on bounded ADMITTED weighted C5
coefficient families. Its statement does not automatically include a ball
of metric corrections controlled by the normal graph. The new matched
subtraction must be tested on a specified domain of inputs. This note
supplies the missing continuous embedding for the FAST graph sector,
including the scalar-to-shape triangular change of variables. It does
not bound H_epsilon-P_epsilon, compute the actual residual, or prove that
the nonlinear map preserves the graph ball.

Retain the actual finite seed Q, paired fermion projector and adjoint,
thermal period, regular bolt and one-sided mass/log phase. Use
x=log z, sigma=4.005, delta=sigma-4=.005 and the full normal multipliers
from current-coupled-local-graph-bound.md. Let

 theta_f=z^sigma Jc_sigma^-1 f,
 b_f=z^sigma B_sigma^-1 g,
 nu_f=theta_f-a_op(D)b_f,
 a_op(s)=(s+3)/(3(s+1)),  D=z partial_z,
 max(||f||_H3,||g||_H3)<=E.                         (1)

The comparison functions are defined on the full x line and then
restricted to the original physical interval. This is the already
selected normal estimator, not a new physical boundary condition.
The actual mass cubic, correlated quartic/slow pair, cutoff commutator
and global matching data remain separate components. Their original
coefficients are not set to zero or selected to reduce the source.

The radius E=1e-11 is the EXISTING conditional joint strong-data target.
Neither component of the actual physical data is certified to meet it.
Changing to (1) does not renormalize a source or declare an inverse to
have norm one. We calculate the physical-field embedding constants.

## 2. Whole-frequency H6 constants for both rows

The retained shape estimate is

 ||B_sigma^-1 g||_H6 <= Cb ||g||_H3, Cb=1801.       (2)

For the scalar row put r=(96pi G c_R-G/(60pi))lambda,
Cv=1+11G lambda/(180pi), A=Cv-r sigma(sigma-3)>0. The previously
proved coercivity gives, for every real t,

 |Jc(sigma+it)| >=6beta(A+r t^2)|sigma-4+it||sigma+1+it|. (3)

The H3-to-H6 constant is bounded by the maximum of

 C_lo=2sqrt(2)/[6beta A delta(sigma+1)],
 C_mid=200sqrt(2)/(6beta A),
 C_hi=(1+1e-4)^(3/2)/[12beta sqrt(A r)].            (4)

The ranges are |t|<=1, 1<=|t|<=100, and |t|>=100 respectively.
In the middle range use <t>^3<=2sqrt(2)|t|^3 and |s-4||s+1|>=t^2.
In the last range use t/(A+r t^2)<=1/(2sqrt(A r)). These are entire
frequency inequalities, not sample-based maxima. Define Ctheta=max (4).

The rational triangular multiplier satisfies

 sup_t |a_op(sigma+it)|=(sigma+3)/(3(sigma+1))=:a_sigma,
 Cnu=Ctheta+a_sigma Cb.                            (5)

It commutes with the constant normal Fourier operators. Consequently
the H6 norms of z^-sigma(b_f,theta_f,nu_f) are bounded by
(Cb,Ctheta,Cnu)E. There is no additional normalization of the fields.

## 3. Exact one-dimensional Sobolev constants

Use the retained unitary Fourier convention and H6 multiplier <t>^6.
For w in H6(R), Cauchy--Schwarz gives

 ||partial_x^k w||_infinity<=c_k ||w||_H6, 0<=k<=5,
 (c_0^2,...,c_5^2)=(63,7,3,3,7,63)/512.            (6)

Indeed c_k^2=(1/(2pi)) integral t^(2k)/(1+t^2)^6 dt. Its beta integral
equals (2k)!(10-2k)!/[2*4^5*k!*(5-k)!*5!]. This proof covers the whole
frequency line. H6 also gives continuous derivatives through order five.
No analyticity of the perturbations is imposed.
The beta integral is [DLMF 5.12.1 and 5.12.3](https://dlmf.nist.gov/5.12#E1),
used with both positive parameters; no domain theorem is imported from it.

Let Y(z)=integral_0^z Q(v)^(-1/2)dv be the original seed's proper
distance from the seam. On its physical half, partial_Y=sqrt(Q)partial_z.
For r_f=z^sigma w(log z) define polynomials P_j(z,s), j<=5, by

 partial_Y^j z^s=z^(s-j) Q^((j mod 2)/2) P_j(z,s),
 P_0=1,
 P_(j+1)=(s-j)P_j+z partial_z P_j,                         j even,
 P_(j+1)=Q[(s-j)P_j+z partial_z P_j]+z Q' P_j/2,           j odd. (7)

In the exponent of Q, `(j mod 2)/2` means zero for even j and one half
for odd j. Every P_j is polynomial in z and in the actual Q coefficients.
There is no inverse power of Q at the bolt. Expand
P_j(z,sigma+d)=sum_k p_jk(z)d^k. If M_jk bounds |p_jk| on a given
interval and Q<=Qmax, then

 |partial_Y^j r_f|
 <=z^(sigma-j) L_j ||w||_H6,
 L_j=Qmax^((j mod 2)/2) sum_k M_jk c_k.              (8)

The finite polynomial recurrence (7) is a coordinate identity, not a
truncation of the geometry or perturbation. Enclose its coefficients on
the actual polynomial Q using interval Bernstein coefficients. Retain
the true correlated z_h; do not replace it by an independently rounded 1.

## 4. The actual outer weighted C5 class

Use the unchanged output collar a=.01. Certify Qmin<=Q<=Qmax there,
Qmin>0. Then z/sqrt(Qmax)<=Y(z)<=z/sqrt(Qmin). The admitted coefficient
norm for the fast part is

 A_out(r_f)=max_(0<=j<=5) sup_(0<Y<=Y(a))
 |partial_Y^j r_f|/[j!100^j Y^(4-j)(1+log(1/Y))].    (9)

Y(a)<.1, so this is the ACTUAL outer portion of the retained phase;
the rest of the seed is treated as interior, not erased. From (8),

 A_out(r_f)<=A6 ||w||_H6,
 A6=max_j L_j a^delta T_j /
                 [j!100^j (1+log(sqrt(Qmin)/a))],
 T_j=Qmax^((4-j)/2) for j<=4,
 T_5=Qmin^(-1/2).                                  (10)

This handles the negative power at j=5 correctly. The maximum of
z^delta/[1+log(sqrt(Qmin)/z)] occurs at z=a. No bound on a nonexistent
smooth fourth derivative of the REFLECTED cubic background is used.
The fast correction itself has zero jets through order four at Y=0,
and its even extension is W3,infinity. A pre-existing mass cubic remains
only C2,1 upon reflection; this calculation does not smooth that term.

Evaluate A6 Cb E, A6 Ctheta E and A6 Cnu E at E=1e-11. Preassign the
outer fast-sector target max A_out<1e-8. Also calculate the sufficient
data threshold E_out=1e-8/(A6 Cnu). This is a bound on a graph-neighborhood
of TEST FIELDS, not a measured value of the physical solution.

## 5. Global multiplier bounds and the bolt

The original domain/adjoint proof requires bounded scalar W1,infinity
factors. The actual seed obeys 0<=Q<.2 on [0,z_h]. Since sigma>1,
(6) gives on the whole physical half

 ||r_f||_infinity+||partial_Y r_f||_infinity
 <= W6 ||w||_H6,
 W6=c_0+sqrt(.2)(sigma c_0+c_1).                    (11)

For j<=3, sigma-j>0. Use (7)-(8) on the complete interval to bound
all three proper derivatives, with z^(sigma-j)<=1. This supplies a
uniform W3,infinity bound for the even fast extension needed by the
weak heat-action argument, including complex graph directions. The
small W1 bound keeps the exponential multiplier invertible; it does
not supply a differentiated trace or a strong source estimate.

At the bolt z=z_h, an H6 function of x gives a C5 function of z in a
neighborhood of z_h. The original regular-disk coordinate has z as a
smooth function of squared radius. Composition therefore gives the
required C5 stationary disk field with its original thermal spin/angle
realization. No puncture condition is imposed and no singular polar
connection is used. On the compact middle and bolt atlas the ordinary
C5 constants are finite, depending on the FIXED seed chart derivatives.
This is a uniform bounded-family conclusion, not a numerical smallness
claim for every atlas seminorm or a physical cold-admission certificate.

Execute the conservative W1 guard W6 Cnu E<1e-3 and its sufficient
threshold E_W1=1e-3/(W6 Cnu). For the quantum Dirac multiplier b alone,
the sharper threshold uses Cb. All these estimates preserve the same
first-order paired domain and complementary adjoint by scalar
multiplication. The induced square domain continues to be defined by
composition; it is not frozen independently.

## 6. What this admits, and what stays conditional

Equations (6)-(11) extend by completion from smooth graph inputs to the
whole fast graph ball. They prove a bounded-family embedding into the
actual weighted C5 mass/log domain, with a uniformly controlled even
W3,infinity collar, global W1,infinity multiplication and regular bolt.
Thus the prior weak C2 action theorem applies to this fast family and
its metric directions. This is a continuum domain statement supported
by explicit constants, not an implication drawn from checker passes.

For the COMPLETE correction, write b=b_asym+b_f and nu=nu_asym+nu_f,
where the actual mass/slow pieces and original global returns are retained.
If their cubic coefficients and weighted remainders are bounded in the
existing phase, the sums are in a bounded admitted family. To certify
the selected quantitative guards it suffices, component by component, that

 A_out(asym remainder)+A6 C_component E <=1e-8,
 ||asym||_W1 + W6 C_component E <=1e-3,              (12)

as well as bounded original cubic coefficients and regular bolt data.
These are requirements on the ACTUAL correlated fields. Do not apply
them separately to the large quartic and slow coefficients and discard
their known cancellation. Their actual global coefficients and E remain
unknown. Report the numerical remaining allowances in (12), not a claim
that they are achieved. The particular canonical shape phase guard in
earlier notes is not identical to all of (9)-(12) or to the complete
physical nonlinear X norm.

No uniform embedding is asserted for a V action-form ball, whose C5
counterexample remains valid. Nor is this a proof that the full older X
ball coincides with the graph family (1), that the Newton map maps into
it, or that its derivative satisfies the contraction estimates. The
graph regularity is an explicit stronger hypothesis, with original units.

RESEARCH DECISION: retain the matched geometric-principal method and use
this quantified fast domain for its test inputs. A source/graph target
that meets (10)-(11) can supply the required fast-sector coefficient
regularity; the action-form estimate alone cannot. The unresolved
operator and global matching estimates, rather than fast-sector
regularity under this hypothesis, determine the next feasibility test.

SINGLE NEXT CALCULATION (PROPOSED, NOT EXECUTED): prove and evaluate the
actual matched H_epsilon-P_epsilon strong remainder on this admitted
graph domain, including the symmetric reflected principal kernel,
coefficient derivatives, finite contacts, common regulator, thermal
tails and original global terms. The canonical target remains .18 and
the combined remaining target .18368. The earlier frozen-cap result is
not a bound for this actual operator. No claim about its magnitude is
inferred from the present embedding.

The actual joint affine <=7.49e-12, known-source <=1e-14, independent
Mellin Efar<=1e-10 and mass/slow moments remain open. The prior sufficient
projection enclosure E_rest<2.123e-42<1.287e-40 is unchanged, with an
unmeasured sign, and is not recomputed. At rho=1e-8 the unproved nonlinear
requirements remain Y0+1e-8 Z0+5e-17 Z2<=1e-8 and Z0+1e-8 Z2<1, with
allocations Y0<=5e-9, Z0<=.25, Z2<=5e7. No actual controlled solution,
exact compatibility, physical matching/causal realization, semiclassical
validity/noise or all-sector stability is certified.

Remaining to-do list: actual matched curved response and common limit;
global sources, mass/slow matching and complete nonlinear admission;
preserved ball and contraction, exact compatibility, physical matching,
validity/noise, all-sector stability, charge/return, opening, transport,
preparation, separate creation and the full stable-wormhole objective.
