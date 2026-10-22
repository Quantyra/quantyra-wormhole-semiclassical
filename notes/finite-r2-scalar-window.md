# A specified R-squared candidate with a positive curvature scalar

2026-10-22, E003/S169, following S167/S168. Construction decision and
ordinary comparison theorem; no charged exterior or wormhole is solved.

The original zero-additional-coefficient model is preserved. Its known
vacuum scalar obstruction and unproved physical frequency prescription
make a different, explicitly specified local gravitational action worth
testing before extending its failed uniform inverse certificate. Add the
Lorentzian density c R^2 with c=10^22 at the SAME scale mu0=1/b0. Keep
one Dirac field, its complete determinant and subtraction, G, magnetic
flux, Maxwell coupling, cosmological constant and prepared temperature.
All other additional finite coefficients stay zero. This is a new theory
parameter, not a change of units, renormalization bookkeeping, or an
evaluation of a previously unknown physical coefficient.

The decisive coefficient-window calculation succeeds: this choice admits
the existing full-Dirac homogeneous-product existence theorem and gives
a positive curvature-scalar mass with m^2 G between 3.31e-25 and 3.32e-25,
and m^2/lambda between 3.39e8 and 3.40e8. The full charged response and
quantum fluctuations remain unproved. Select this candidate for the next
physical construction test, retaining S167 as an incomplete alternative.

## 1. Total trace coefficient, with the Lorentzian sign fixed

The [retained anomaly](spherical-anomaly-action.md) has b=1/30 and
epsilon=G/(60pi). Its additional Lorentzian c R^2 changes b by
-192pi^2 c. Equivalently, varying its metric tensor gives

    H_ab=2 R Ric_ab-(1/2)g_ab R^2
                         +2(g_ab Box-nabla_a nabla_b)R,
    tr H=6 Box R.

The Einstein equation gains 16piG c H_ab. The repository's normalized
trace is minus the trace of that tensor equation, so the COMPLETE trace is

    T_c=R+12lambda0
       +epsilon[7 C^2/12+11 J^2/6-11 R^2/72-20 M_em]
       +(epsilon-96piG c) Box R.                       (1)

Here J is traceless Ricci. The nonlocal conformally invariant determinant
still has zero trace; it has not been removed from the shape equations.
On arbitrary compact conformal variations, delta integral R^2 is
-12 integral sigma Box R. The checker also derives this identity from an
unexpanded time-dependent conformal metric, independently of (1).

Do not confuse c with the positive EUCLIDEAN representative
1/(5760pi^2) in the [field-redefinition note](stationary-covariant-reduced-inverse.md).
That representative already describes the retained Dirac Box R term.
Here c is ADDITIONAL in the Lorentzian action. No term is counted twice.
Changing a determinant representative while compensating its finite local
coupling leaves epsilon-96piG c invariant and cannot change the branch.

## 2. Execute the exact vacuum comparison

Let g=G lambda, where lambda is the same positive root
g+11g^2/(360pi)=2000^-10. Constant-curvature R^2 variations vanish on
an Einstein metric, so the outer vacuum relation is unchanged. On the
same exact AdS conformal-vacuum comparison as [S168](ads-fast-scalar-stability-preflight.md),
put

    C_v=1+11g/(180pi), chi=96pi c-1/(60pi).

For psi=delta R, (1) gives

    delta T_c=(C_v-G chi Box)psi,
    (Box-m^2)psi=0, m^2 G=C_v/chi.                     (2)

Thus c_*=1/(5760pi^2) is the sign threshold. At c<c_* the mass is
negative; c=c_* removes this differential trace branch and must be
treated as a different-order problem; c>c_* gives positive m^2.
The original c=0 gives exactly the previously retained negative mass.

The full conformal tensor lift also follows, with its original limited
scope. With Q_ab=nabla_a nabla_b-g_ab Box+3lambda g_ab,

    delta E_ab[2s g]=-2 Q_ab(C_v-G chi Box)s,
    delta R=-6(Box-4lambda)s.                          (3)

Conservation, maximal symmetry and the polynomial trace determine (3),
as in the prior proof. It solves all vacuum tensor components when
(Box-m^2)s=0. For the chosen coefficient m^2/lambda>10^6, so the lift
has nonzero curvature. This is not a lift on the charged non-Einstein
metric, nor a statement about independently varied quantum states.

The radial KG operator in global AdS is now

    A_c=-partial_z^2+(2+m^2/lambda)/sin(z)^2.

For the regular-center, finite-energy realization its quadratic form
is positive. Since m^2/lambda>10^6, the outer endpoint is limit point:
the indicial powers of the curvature are

    Delta_+/-=3/2 +/- sqrt(9/4+m^2/lambda),

with only Delta_+ admissible in the KG norm. There is no two-amplitude
oscillatory scalar boundary phase in this comparison. All higher angular
momenta add a nonnegative centrifugal term. This supplies scalar vacuum
energy control, not positivity of every quantum metric sector.

Independently, u(x)=sin(x)cos(x)^(Delta_+-1) is the actual spherical
vacuum ground radial mode, with frequency squared Delta_+^2. The checker
substitutes it into the differential equation without discretizing a
charged kernel. Standard scalar-tensor equivalence and AdS scalar
theory are established results, not new wormhole results.

## 3. Simultaneous mass and homogeneous-existence window

Choose explicit DESIGN targets

    m^2 G <= 10^-24,        m^2/lambda >= 10^6.         (4)

These separate the candidate scalar frequency from the Planck and
background curvature scales. They are not, by themselves, an EFT error
theorem. Equation (2) makes the exact coefficient requirements

    [10^24 C_v+1/(60pi)]/(96pi) <= c,
    c <= [C_v/(10^6 g)+1/(60pi)]/(96pi).              (5)

The [existing finite-local reference theorem](finite-local-reference-window.md)
uses, for c1=c and other coefficients zero,

    N4=32c, Lambda_ref=8piG/a0^2
             =48pi C0/[2000^10(C0-1)],
    C0=sqrt(1+12pi*2000^6),
    Lambda_ref N4 <= 10^-4.                           (6)

Use different names for Lambda_ref and the physical AdS curvature.
The exact additional upper bound is c<=10^-4/(32 Lambda_ref).
The rounded closed interval

    4e21 <= c <= 3e24                                 (7)

is wholly inside (5)-(6). Rational enclosures of pi, C0 and the original
g root verify the endpoints and monotonicity. This is a nonempty range
of THEORIES, not freedom to vary c while solving the fixed selected one.

At c=10^22 the conservative certificate gives

    Lambda_ref N4 <4.713e-8,
    D_ref=epsilon_ref B_ref+Lambda_ref N4
                                      <7e-23+4.713e-8,
    |b_*^2/b0^2-1|, |a_*^2/a0^2-1| <9.43e-8.        (8)

Here epsilon_ref B_ref<7e-23 is the established full-Dirac pressure
budget; these symbols are not the trace epsilon or corrected radius B.
The contraction theorem includes the ENTIRE determinant, not just a
local-only numerical root. It proves a unique homogeneous product in
the stated one-percent box. Its radii are NEW a_*,b_*. They are not
identified with the old zero-c radii or inserted into an old receipt.

The magnetic lowest massive parameter remains quantitatively separated:
with q=2000^4 and R=9.43e-8,

    sqrt[(q+1)(1-R)/(C0(1+R))] < M1
                < sqrt[(q+1)(1+R)/(C0(1-R))],
    18.04 < M1 <18.06.                                (9)

Thus this first product test does not destroy the original full-tower
mass separation. It does not re-establish the finite exterior source
theorems, which depend on more than M1. The product stress of c R^2 is
pA=4c(y^2-x^2), pS=-pA, x=a_*^-2,y=b_*^-2. It has zero orbit null
contraction, so it supplies no replacement for the quantum opening source.

## 4. An explicit local boundary prescription and its limits

Treat the large c term as a retained scalar-tensor interaction, not a
small perturbation at its own scalar mass. Its exact auxiliary action is

    integral sqrt(-g)[(R-2Lambda)/(16piG)
                          +c(2varphi R-varphi^2)] .   (10)

The varphi equation gives varphi=R. Set F=1+32piG c varphi. For F>0,
the Einstein-frame kinetic density of log F is
-3(grad log F)^2/(32piG). This is the ordinary healthy local scalar
kinetic sign. It is not a ghost analysis of the full Dirac effective
action or a bound on its induced fluctuations.

On the proved product box,

    F >=1-8c Lambda_ref/(1-R)>1-1.20e-8>0.            (11)

The constant-curvature vacuum also obeys F=1-384pi c g>0. Thus the
local scalar-tensor map does not fail on either comparison.

On a smooth non-null finite boundary the local action (10) uses

    S_boundary=epsilon_n integral sqrt(|h|)
                         [1/(8piG)+4c varphi] K,
    delta h_ij=0, delta varphi=0.                     (12)

This cancels the local metric normal-derivative variation; the remaining
scalar boundary variation is zero under the displayed datum. After
eliminating varphi this is a curvature datum, not merely a leading metric
datum. [Dyer--Hinterbichler, (1.4) and section 5](https://arxiv.org/pdf/0809.4033)
derive this established boundary correspondence. No claim that the
complete renormalized boundary problem is now solved follows: the Dirac
paired domain, Weyl boundary terms, AdS counterterms, charge, ensemble
and cutoff limit must all be matched with (12) for the new candidate.

The scale hierarchy in (4) avoids the original scalar's Planck-frequency
problem in this comparison. It does not bound omitted higher operators,
photon/graviton loops or state-dependent noise. For any proposed physical
band and observable one must still bound the omitted response relative
to the retained response, and its propagated metric fluctuations.
[Hu--Verdaguer, sections 3.3 and 6.5](https://arxiv.org/html/0802.0658)
distinguish mean-field response, induced fluctuations and prescriptions
for runaway modes. Their results do not certify this charged candidate.

## 5. Research decision, transfer ledger and verification

Select c=10^22 as an explicitly different candidate for the NEXT complete
construction test. The reason is the executed overlap of a positive
vacuum scalar, a large frequency separation and the full homogeneous
existence window. Neither an oversized old bound nor the vacuum comparison
is used as an impossibility theorem for the original construction.
S167/S168 and their artifacts remain intact and incomplete. Their known
cold radius 2.5e-10 and response constants are not imported as quantitative
bounds for a newly centered exterior.

Next derive the actual charged scalar-tensor equations and mixed response
with (12), recenter on the new product, and test whether the realized
coupled scalar/shape boundary operator has a useful inverse. Positivity
of the vacuum scalar is not an invertibility proof for that operator.
In particular do not set the old signed moment to zero: the old operator,
projection and source belong to a different action. Any source/range
condition for the new physical realization must be derived afresh.

Reuse only identities or parameter-uniform proofs whose hypotheses have
been checked on the new background. Recompute radius/temperature matching,
finite residual, solver norms, full nonlinear remainder, physical charges
and stability. The complete angular tower and negative-energy mechanism
remain required. A positive scalar vacuum mode does not solve a wormhole.

Freeze this note, source review and checker before execution. Reconcile
the previous canonical-admission receipt and every retained input hash.
Verify the nonlinear R^2 variation, total-coefficient sign, vacuum lift,
exact window inequalities and homogeneous pressure formula independently.
At 60/90 digits evaluate the benchmark, compare the explicit AdS radial
mode with its independent differential equation, and test the auxiliary
Legendre elimination. All logs and failure snapshots are retained. These
checks are supporting calculations; no charged state is sampled, no
continuum proof is machine certified and no Lean run is selected.

Remaining to-do list: new charged scalar-tensor/boundary equations and
coupled inverse; actual residual/correction, nonlinear solution and exact
compatibility; full physical matching, validity and stability; charge/
return, opening, transport, preparation, separate creation and full audit.
