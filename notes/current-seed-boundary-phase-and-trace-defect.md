# The current reference correction's physical outer phase and trace defect

2026-10-22, E003/S169. Ordinary argument before supporting verification.
The full stable semiclassical wormhole objective remains incomplete.
This executes the boundary check inside the selected finite-domain BVP.
The closed feasibility audit is not reopened.

## 1. Question, fixed construction and result

Use exactly the current finite seed and evaluated complete reference
correction h1=(f1,a1,s_R) of
[the fixed-clock test](fixed-clock-seed-bolt-feasibility.md). All theory
coefficients, physical flux, scale, paired HH domain, temperature and
outer clock stay fixed. The reference action is a covariant spherical
potential at the CURRENT product radii B,A, held fixed under variation.

The correction has the physical outer mass/fast phase through five areal
derivatives. Its renormalized fixed-clock boundary work vanishes by the
retained phase theorem after the hypotheses below are checked. There is
no new finite-wall condition, freely tuned mass, or two-zero scalar-tail
condition. This closes a candidate-specific boundary test, not the closed
full metric Hessian domain, its inverse or numerical admission norm.

There is a nonzero forced r^-4 curvature coefficient. Its leading full
trace defect can be evaluated exactly: the necessary full linear response
coefficient is larger than the reference coefficient by a relative amount
<1.18e-8, or <2.35e-36 in units lambda B^4. This small tail defect is
correctible at the level of the outer indicial equation. It is not a bound
on its global bolt response or on the rest of the full Hessian defect.

## 2. Transfer the fixed-seed endpoint regularity before using the phase

Let O_j(r^-n L) mean the estimate after every areal derivative of order
at most j, with the corresponding additional inverse powers; L=1+log(r/B).
On the current ANALYTIC one-sided compactified seed, the complete evaluated
sphere-integrated residual source satisfies

    rho,p=O_4(r^-2 L).                                      (1)

This statement is stronger than absolute mass integrability. It is not
obtained by differentiating a previously proved unsigned bound.
Here is the retained parametrix argument and its transfer:

* With z=B/r, the canonical polynomial is
  calB=x+z^2-2a z^3+(1+3x+d)z^4-2d z^5+d z^6.
  It is positive at the face, analytic on its one-sided collar and has
  calB'(0)=0. The boundary, paired projectors, complementary mixed square
  and sphere connection are smooth there. The mass cubic is unchanged in
  form; the new terms start at degree four and higher. We do not require
  a smooth even double of the mass cubic.
* The actual square is positive and has no kernel at this fixed positive
  temperature. The massive inverse estimate and the disk Fourier proof
  for the LLL apply to this seed, as already checked in the canonical
  matter neighborhood. No numerical absolute LLL gap is assigned.
* Apply the full four-dimensional direct/reflected boundary parametrix
  BEFORE the common finite part, as in section 5 of
  [the physical endpoint proof](recentered-physical-endpoint.md).
  Subtract the finite boundary profiles through the needed order, retaining
  their exact finite parts. The leading canonical mass profile and its
  Weyl partner cancel. The physical vacuum is removed by the exact lambda
  relation. The remaining leading powers contain at most one logarithm.
  The separated large-time kernel is smooth at the face by the fixed
  elliptic, kernel-free realization. Its global state coefficients are
  retained and not set to zero.
* Four differentiations of the evaluated source need one more conormal
  derivative than the earlier three-derivative endpoint statement. On this
  fixed analytic seed this uses one more differentiated parametrix term:
  remove finitely many additional terms so that the differentiated remainder
  is still one power beyond the retained expansion. These are added back
  exactly; no new finite subtraction is chosen. Normal differentiation
  consumes one power, and multiplication by z^2 for rho,p converts each
  estimate into (1). Subleading z log^N z terms are absorbed in the
  leading log bound for any fixed N. This is a finite fixed-seed argument,
  not differentiability with respect to a rough unknown metric.
* The local Weyl terms are explicit rational/log polynomials. The added
  R2 and reference-potential sources are rational. Their vacuum parts
  cancel; their remaining sphere-integrated terms start at r^-2. All
  four differentiated versions therefore have (1) as well.

The differential-boundary parameter-ellipticity input is
[Gaarde--Grubb, Theorem 2.4](https://arxiv.org/pdf/math/0703878).
Its differentiated resolvent estimates permit a finite number of normal
derivatives on this smooth fixed seed. The indices in its logarithmic
estimate (2.24) are met by taking the subtracted expansion sufficiently
far. Proposition 2.9 retains a singular logarithmic Green term; it is not
legitimate to call the entire logarithm smoothing. The concrete
subtractions and cancellations remain those of the retained source proof.
The theorem alone does not prove renormalized stress differentiability.
[Vassilevich, sections 3.3 and 5.3](https://arxiv.org/html/hep-th/0306138v3)
provides the mixed heat coefficients; their smoothness/strong-ellipticity
hypotheses are used on the seed, not on a fictitious smooth even double.
Neither reference supplies numerical constants in (1).

## 3. Integrate the actual constraints; no new boundary datum

The preceding fixed-clock test gives an integrable rho and

    delta m=kappa R s_R+G integral_R^infinity rho dr,
    f1(r)=-2delta m/r+(2G/r)integral_r^infinity rho(t)dt,
    a1(r)=-G integral_r^infinity (rho+p)/(t F(t))dt.       (2)

Since F=lambda r^2+1-2m/r+O(r^-2), (1)-(2), including differentiated
fundamental-theorem identities, imply

    f1=-2delta m/r+O_5(r^-2 L),
    a1=O_5(r^-4 L).                                    (3)

For example integral_r^infinity t^-2(1+log(t/B))dt
=r^-1[2+log(r/B)]. Four source derivatives suffice for five derivatives
of the integrated response. The bolt cancellation and period are already
fixed by the reference constraints; (2) does not change either of them.

The current seed itself obeys the SAME phase, with its additional r^-3,
r^-4 metric terms and R4+12lambda=O(r^-5). Equation (3) is therefore an
allowed tangent in the phase of
[the recentered mass-boundary theorem](recentered-mass-boundary.md).
This is a qualitative finite-norm statement. It does not give a numerical
bound in the proper-collar N_(4,5) norm or in the global admission space.
The X^4 log X tail is allowed by that WEIGHTED norm; claiming ordinary
closed-boundary C5 regularity would be false and is unnecessary.

The complete fixed-source boundary momentum has coefficient
-105-116+221=0 in the same lambda delta m/(3780pi) units. The
canonical distribution has no surviving pairing with the mass/fast
tangent, and the residual curvature-square boundary term vanishes as
r^-1 L. Thus the physical boundary one-form is zero here. Its derivative
within this fixed-source phase is also zero, as an identity of that
one-form. This does not promote the C1 action theorem to a quantitative
C2 source theorem. The stationary charge difference, where needed, is
Z_c delta m/G with the already fixed Z_c, not an independently chosen
boundary condition. No finite wall or additional counterterm is added.

This is a boundary compatibility result for h1. A closed graph realization
and inverse for the complete current finite metric Hessian must still be
proved. They do not follow from vanishing boundary work alone.

## 4. Exact curvature of the reference correction

Use the anomaly-complete trace, with epsilon=G/(60pi),
D_R=96piG c_R-epsilon and J=Ric-R4 g/4:

    T[g]=R4+12lambda0
       +epsilon[7 C4^2/12+11 J^2/6-20 I-11 R4^2/72]
       -D_R Box R4,                 I=q^2/(2r^4).      (4)

The evaluated current seed trace T0 is the exact rational expression in
[the product-centered scalar proof](product-centered-scalar-cold-response.md).
Its outer leading term is

    T0=-10epsilon q^2/r^4+O(r^-5).                    (5)

Indeed R4+12lambda starts at r^-5, Box preserves that leading decay,
C4^2 begins at r^-6, J^2 at r^-8, and the constant vacuum cancels.

At fixed areal r, the mixed trace of Gamma_ref is a fixed function of r,
independent of F and alpha. Differentiating the REFERENCE trace equation
therefore gives exactly

    delta R4[h1]=-8piG tr T_rest=-T0.                 (6)

The added R2 source is included in T_rest and T0. Thus this formula
uses the full local anomaly of every angular mode, even though the
traceless state-dependent source is not local. The first reference
curvature is known exactly, without differentiating a source envelope.
At order r^-4, the seed has no coefficient, so

    [R4_seed+delta R4[h1]+12lambda]_(r^-4)=a_ref,
    a_ref=10epsilon q^2>0.                            (7)

An areal-gauge change cannot alter this leading coefficient: the seed
curvature derivative is O(r^-6), and the allowed return cannot create an
r^-4 contribution from it. This is a physical tail coefficient, not a
free normalization.

## 5. Execute the full trace indicial consistency test

Let a_full/r^4 be the leading curvature correction required by the FULL
linearized current equation with the same fixed physical leading data.
For this power Box(r^-4)=4lambda r^-4+O(r^-6). Variation of the local
curvature anomaly contributes (11epsilon lambda/3)delta R4 at leading
order. The magnetic invariant is fixed in areal gauge, and the other
curvature variations, moving-density contact and Box variation on the
seed curvature enter only at higher decay. Consequently

    chi a_full=10epsilon q^2,
    chi=C_v-4D_R lambda,
    C_v=1+11epsilon lambda/3,
    chi=F_lambda+(23/3)epsilon lambda,
    F_lambda=1-384piG c_R lambda.                     (8)

In particular

    a_full=a_ref/chi,
    [T0+D T[g0]h1]_(r^-4)=(chi-1)a_ref.              (9)

This is a nonzero defect, not an assertion that h1 already solves the
full metric equations. The exact benchmark and current product enclosure
give conservative bounds

    1-1.179e-8<chi<1,
    0<a_full/a_ref-1<1.18e-8,
    0<a_ref/(lambda B^4)<1.99e-28,
    0<(a_full-a_ref)/(lambda B^4)<2.35e-36,
    abs((chi-1)a_ref)/(lambda B^4)<2.35e-36.           (10)

These refer to a coefficient at infinity. They are neither a uniform
residual bound nor a cold amplification estimate.

For clarity, if f1 has an r^-2 log(r/B) coefficient F_L and a1 an
r^-4 log(r/B) coefficient A_L, direct curvature differentiation gives

    a_ref=F_L+8lambda A_L.                            (11)

At least one log coefficient is nonzero. Setting all logarithmic tails
to zero would already be incompatible with the actual magnetic anomaly.
The existing physical weighted phase allows them. Equation (11) does
not determine F_L and A_L separately: their traceless row and global
state remain part of the full solve. No arbitrary split is made.

The massive curvature homogeneous roots are
Delta_+/-=3/2 +/- sqrt(9/4+C_v/(D_R lambda)); Delta_+>18000.
The growing branch is excluded by the fixed leading data; the decaying
branch is allowed and fixed by the global solution. No second independent
zero-trace condition is justified. Since chi>0, the forced r^-4
coefficient is nonresonant in the curvature equation. There is no scalar
tail incompatibility in this test.

## 6. Consequence for the selected BVP

The first reference candidate has the correct physical outer phase and
zero fixed-source boundary work. Its leading curvature coefficient needs
the explicit small correction (9)-(10), not a new boundary ensemble.
Thus this boundary diagnostic gives no reason to switch construction.
It also does not close the full finite Hessian domain or its global
inverse. Keep the original shape row and its global kernel moments;
a differentiated equation that loses an r^-4 moment is not a substitute.

The next work remains the SAME full finite-domain matched BVP/defect:
realize the complete Hessian on the physical mass/fast phase and bound
its correction to h1 in the physical weighted return. Required cold value
difference is <1.667e-9, total Y0<=5e-9, Z0<=.25 and Z2<=5e7 on
rho=1e-8 with all admission margins. None follows from the small number
in (10). No additional generic source-estimate chain is selected.

Executed: candidate-specific outer phase and boundary work, exact reference
curvature, and full-trace leading tail defect. Uncompleted: the full BVP,
its numerical weighted norm, preserved nonlinear solution, physical
matching, semiclassical validity and all-sector stability.

Remaining to-do list: full matched finite-domain Hessian realization and
weighted correction; nonlinear control, matching, validity, stability and
the remaining stable-wormhole requirements.
