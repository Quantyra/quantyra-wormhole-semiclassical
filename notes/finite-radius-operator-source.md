# Direct scalar radius operator from the complete source equation

2026-10-22, E003/S167. Ordinary argument before supporting verification.
Keep the original finite comparison, complete angular tower, state,
paired domain, subtraction, mass and fixed return. The preceding
scalar bound is dominated by estimating L4nu through eighth metric
derivatives. Use its exact radius equation and existing full source
derivatives instead. This is a bound on the same scalar residual,
not a new solution, source or boundary realization.

## 1. Combine the seed coefficients before taking bounds

Use z=1/t, b=calB/C, gamma=Gq^2/(CB^2), c_x=C/x<6.1 and
the original common-coordinate equation

    L4nu=c_x z^2[2S u+c_nu nu-F+gamma z^2 T+gamma E_nu],
    S=-b_z/z+2b/z^2,
    c_nu=z b_zzz/6+4(x/C)/z^2+2e z^3-4e z^4,
    E_nu=z^2 E_delta.

On the cold patch T and E_delta include their original homogeneous
product subtractions. On the middle patch they have the original
unsubtracted canonical/local definitions. In both cases

    Rdelta=B^2(p-rho)_rest/q^2=z^2(T+E_delta),
    L4nu=c_x[Au u+Anu nu-z^2 F+gamma z^2 Rdelta],
    Au=-2z b_z+4b,
    Anu=z^3 b_zzz/6+4(x/C)+2e z^5-4e z^6.

The z factors cancel the apparent singular coefficients. For
b=b0+b2 z^2+b3 z^3+b4 z^4-2e z^5+e z^6, exactly

    Au=4b0-2b3 z^3-4b4 z^4+12e z^5-8e z^6,
    Anu=4b0+b3 z^3+4b4 z^4-18e z^5+16e z^6.

Use b0<1/5.9, |b3|<.668, |b4|<.501, |e|<1e-24.
These are regular polynomials on the entire cover. The original
radius-equation anomaly cancellation is retained, as are the Lie
terms in its derivation. No unsigned term is canceled by a bound.

## 2. Fixed reference derivatives and mass conversion

Use D=z0 nabla_gamma, W=t0^2 and q=1+Zeta at a fixed observation,
as in the [finite coefficient proof](finite-higher-local-coefficients.md).
Every polynomial z^m has divided D-jet majorant z0^m q^m.
Thus W z0^2=1 implies

    W jets_D(z^2 F) <= q^2 jets_D(F),
    W jets_D(z^2 Rdelta) <= q^2 jets_D(Rdelta).

The right-hand side arrays are unweighted. There is no attempt to
bound W F on the outer end, where its mass term would diverge.
The polynomial terms Au u,Anu nu use the already established
W-weighted tangent array K; it bounds both u and nu.

For F on the middle patch, convert the frozen per-order mass radial
array using the existing canonical tensor recurrence and drop only
the harmless factor z0^j<=1. In the cold region use the improved
uniform F array from [the divided mass reconstruction](finite-low-gap-reconstruction.md),
with both the regular Cartesian and away conversions. Every tensor
slot and the zeroth connection remain included.

On the outer collar, the exact mass split is

    F=-2(x/C)z(m+mu), |m|<4e-11.

The existing f_tail array bounds its mu part with weight Y^2 L(Y).
The mass part has ordinary divided array A(.14*4e-11)*P(1), with
weight Y, where A(c)_j=c(.02)^j and P(m)_j=|binom(m,j)|/100^j.
Convert both arrays by the outer tensor recurrence. Its connection
insertion gains Y^2 independently of the leading weight m, so the
same coefficients (9/400)j!2^j apply for m=1 and m=2. Since
z/Y<.412, Y<=.1 and Y^2 L<.0331, the unweighted D-jet bounds are

    (.412)^j[.0331 tail_tensor_j+.1 mass_tensor_j]100^j.

Take the maximum of these three regional arrays. This retains the
exact mass mode and the complete mass tail without a moving endpoint.

## 3. Full source tensors and explicit local additions

The [ninth source ledger](seed-nine-metric-eleven.md) supplies m_j and
c_j as FULL canonical source-tensor bounds in units j!100^j q^2,
with cold weight a(s). Its proof retains every ordered connection
slot; those entries are not reconstructed from radial derivatives.
Use m_j 100^j on the middle patch and .00331 c_j 100^j on the
cold patch. The extra D factor z0^j is at most one.

The canonical LLL term is -kappa b_zz, kappa=C/(24pi q)<5e-5.
The middle full-jet bound is kappa B''(q), where B is the old
positive seed polynomial. On the cold patch subtract the product
first. With tau=1-z and delta_bar=delta/C,

    b_zz-2=-6B3 tau+12B4 tau^2-80e tau^3+30e tau^4-6delta_bar z.

Use |B3|<1.334, |B4|<.501, delta_bar<1.001e-60 and the full
D-jet majorant .001+Zeta for tau. This retains the actual local
subtraction and all tensor derivatives. It does not identify the
old radial LLL allowance 2^j with a full tensor allowance.

The explicit E_delta formulas have .005-disk radial bounds
1e16*2^j on the middle patch and 1e22 s0*2^j on the cold patch.
For a uniform cold bound use s0<=.001. Convert these actual radial
jets with the middle and cold full tensor recurrences, including
the Cartesian bolt argument. Denote the resulting D-jet arrays by
Em and Ec. These are derivatives of the local formulas; no unsigned
real value inequality is differentiated.

Since Rdelta=z^2(T+E_delta), safe bounds for the source contribution
W jets_D(z^2 Rdelta) on the two bounded regions are

    q^4 [m_j100^j + LLL_middle + Em],
    q^4 [.00331 c_j100^j + LLL_cold + Ec].

The extra z0 powers have only been enlarged to one. The source
arrays still describe the original complete angular sum and state.

## 4. Outer physical cancellation precedes the tensor bound

Use the [complete outer source proof](outer-metric-ten-jets.md).
At each derivative order j<=8 let K_j=max_(i<=j)o_i, with o_i
the retained outer ordinary source entries. The proof gives

    |Rdelta^(j)|/[j!100^j Y^(2-j)L(Y)] <=2K_j.

The principal 1/Y profile, local Weyl and reference vacuum terms
have already been combined in that physical source. The explicit
LLL/local reserve remains included. Do not differentiate or
integrate the divergent sectors separately.

Apply the outer full tensor recurrence directly to the nonconstant
array 2K_j. The same weight-two calculation gives its unweighted
D-jet bound (.412)^j*.0331*tensor_j*100^j. Multiply by q^2 for
W jets_D(z^2 Rdelta). Take the coefficientwise maximum with the
two bounded-region source arrays in section 3. Call this array Q.
This is a bound on the complete physical pressure difference,
not merely the outer canonical trace remainder.

## 5. Replace only the radius-operator allowance

Let Fu be the full unweighted mass array and let Ap,Np be the
positive Au,Anu polynomials evaluated at q. Then

    LN_direct=6.1[(Ap+Np)K+q^2 Fu+6.8e-28 Q]

bounds W|D^j L4nu|/j! through eight. Take its minimum with the
previous generic LN at each order. Keep the independent bound
W|L4nu|/tau<.005 for the value-times-high-gap product.

Recompute only the scalar curvature remainder term containing
V L4nu, retaining its finite exponential factors and every cold
factor. The other five remainder arrays are unchanged. Take the
minimum with the preceding complete curvature remainder through
eight. Recompute every scalar-trace-dependent anomaly term using
the new DZ, including the Laplacian and trace-square terms. Keep
the old valid curvature-square and magnetic arrays. The finite
geometry, reference-potential remainder and Lie term are unchanged.

Use the complete residual array, actual covariant conversion and
full noncompact physical volume as before. The exact quadratic
path zero and t^-4 decay remain. Report the new scalar H6
coefficient, upward enclosure and scalar-only transfer budget;
no full Schur-source or transfer pass is assumed.

## Verification and independent scope

Freeze this note, checker and review before execution. Check the
combined polynomial identities, the complete physical source
normalization, every regional full-tensor conversion and the old
radius-product reconstruction. Verify that exactly one curvature
remainder term changes and that all unrelated anomaly terms remain.

Independently compare the scalar Laplacian obtained from the full
four-dimensional volume divergence with the complete radius source
equation on nonconstant regular manufactured fields. Use a flat
canonical orbit times a round sphere, z=z0-(x^2+y^2)/4, and a
source chosen from the separate radial equation. Compare every
mixed Cartesian derivative through the selected orders up to eight,
including the center, at 60/90 digits. Dropping the nonzero local
source must fail. These fields test the identities and tensor slots;
they are not samples of the actual charged quantum source.

Remaining to-do list: assess useful complete source and transfer
bounds; quantum/TT H4 and higher lift/inverse; nonlinear mapping,
preserved neighborhood, signed compatibility and existence; physical
validity, state/boundary, matching, stability, charge/return, opening,
transport, preparation, separate creation and audit.
