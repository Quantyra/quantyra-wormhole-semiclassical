# Separate source orders in the scalar reconstruction

2026-10-22, E003/S167. Ordinary argument before supporting verification.
Keep the physical source, complete angular tower, state/domain, common
subtraction, mass, finite path and return. Use the existing source
derivative bounds separately, rather than replacing all of them by
the largest order-eight reserve. Propagate the improvement into the
scalar curvature remainder and its anomaly Laplacian. The receipt
reports the resulting exact rational scalar H6 coefficient and its
quantitative consequence; no transfer pass is presumed.

## 1. Complete per-order canonical sources

Let m_j,c_j,o_j, j<=8, be the frozen middle, cold-total and outer
ordinary entries in the [source ledger](seed-nine-metric-eleven.md).
Their units are j!100^j q^2, with the retained cold a(s) or outer
logarithmic weight. These are the actual seed-source bounds; the
auxiliary C11 step already exists and no new cutoff is introduced.

The [middle reconstruction](middle-metric-ten-jets.md) retains the
explicit canonical LLL term -q C b_zz/(24pi). Its complex modulus
on the .005 proper disk is below q^2, so its normalized jth radial
jet is below 2^j. Use

    Tm_0=1e20, Tm_j=min(2e44,m_j+2^j), 1<=j<=8.

The old value and pressure constants remain valid. For the cold
product-subtracted trace, the extra LLL term is
-C(b_zz-2)/(24pi q). On the cold complex disk |tau|<4s,
the centered seed polynomial bounds |b_zz-2|<33s. Since
C/q<.00352 and 24pi>72, its modulus is below s<a(s).
Cauchy again supplies 2^j a(s). Thus use

    Tc_j=min(2e44,c_j+2^j), j<=8.

Every explicit Weyl/reference local allocation remains as in the
reconstruction proofs. No cancellation is inferred from an unsigned
source bound. In particular E_u and E_nu retain their distinct
middle/cold analytic reserves; the anomaly cancels from the radius
equation only by the original exact identity.

For the outer calculation define K_n=max_(j<=n) o_j. This is at
least 4e19, hence contains all the old 1e8 analytic local reserves.
The [outer source-to-metric proof](outer-metric-ten-jets.md), applied
at each truncation n, bounds the complete rho and p-rho source jets
through n by 2K_n. It keeps the principal 1/Y, Weyl, vacuum and
reference cancellations before integration. The explicit local LLL
allocation is still included in its 1e8 local reserve.

## 2. Reconstruct without a uniform high-order replacement

Use the exact positive recurrences of the three metric proofs, with
the following substitutions only. All coefficients below are in the
same factorial/100 normalization, and gamma<g0=6.8e-28.

Middle: replace its constant source array by Tm in both Ward/mass
and coupled metric recurrences. Keep P0=1e20,F0=1e-8 and
U_j=N_j=1e-4 through five. Each source occurrence is Tm_(n-i),
not Tm_n outside the convolution. The analytic local terms remain
1e16*2^j and 1e18*2^j.

Cold: replace the pressure arrays by

    Qnear_j=5 sum_(i=0)^j (.1)^i Tc_(j-i),
    Qaway_0=Tc_0/2,
    Qaway_(n+1)=sum_(i=0)^n 5*2^i
               [Tc_(n-i)+2Qaway_(n-i)]/[100(n+1)],
    Q_j=max(Qnear_j,Qaway_j).

The regular Ward kernel is unchanged; its differentiated source
weight is at most the observation weight. In the mass and metric
recurrences, replace each old constant K by the corresponding Tc
entry inside its convolution. Keep F0=1e-10, low U,N=1e-4, the
factor two from the regular bolt integral, and the analytic local
reserves 1e22 and 1e24. No extra source derivative is required.

Outer: keep mu0=7e-11 and use

    mu_j=8 gamma_x K_(j-1)/(100j), 1<=j<=9,
    gamma_x<4e-27.

All primitive and analytic-multiplier convolutions remain those of
the original outer proof. In the proper v_(n+2) recurrence replace
6 gamma K by 6 gamma K_n. The drift, Lie and local trace terms
are unchanged. Keep the original proper v,w low bounds 1e-8
through five, return to the original common u,nu,V, and add the
analytic mass remainder reserve 1e-9 before full tensor conversion.
The correlated mass cubic remains separate and unchanged.

These recurrences are triangular and monotone. The new source
arrays lie below the old uniform ones, so they cannot worsen a
retained higher bound. Verify this explicitly, including reproduction
of the former arrays when the old uniform sources are restored.
This is a quantitative refinement of the same tangent fields, not
a new metric family or a new regularity assumption.

## 3. Full tangent tensors and the retained cold factor

Apply the existing middle, cold-near, cold-away and outer tensor
recurrences through ten. Middle carries W<=50.1^2; cold u,nu
carry a(s)<=.00331 and W<1.001^2. For the outer remainder, at
order j use the maximum of its ordinary entries through j plus
1e-9, the frozen outer tensor amplification, and the factor
9*.0331. Add the separate cubic allowance 18e-11*.1. This
explicitly includes both the analytic mass remainder and cubic.

For V=nu/tau use the same relative inverse-tau radial coefficients
1.2(.2)^i/tau. The unweighted middle enlargement is 1200 and the
cold a/tau factor is 82. Construct W tau D^jV separately: middle
uses 1.2 instead of 1200; cold away retains a before multiplying
by tau; cold near uses a/tau<82 along the Taylor segment and
tau<.001 at the observation. These are exactly the estimates of
[the gap-product proof](finite-scalar-cold-gap-product.md).

Take coefficientwise minima with all previously proved tangent
bounds, including the small values, gradients and gap Hessian.
This produces U,K,Vbar,KV and h through ten, where h_j bounds
W tau|D^jV|/j!. The old finite-field S array remains available;
no higher endpoint-dependent coordinate map is differentiated.

## 4. Scalar curvature remainder through eight

Use the complete six-term exact scalar remainder from the preceding
two scalar proofs. Its tensor product rules hold at every available
finite order. Extend the tau-retaining gradient-square formula and
the V L4nu formula through spatial order eight. Terms where two or
more derivatives hit tau keep the ordinary tau-jet bound; they are
not discarded when their other factors contain higher derivatives.

The L4nu array uses two extra tangent derivatives, so order eight
uses the original metric/tangent order ten just reconstructed. The
source order is still eight. Retain every path and exponential
factor and the same W^2 decay. Let Nnew be the resulting bound
for W^2 N_R/ell^2 through eight. Take minima with the old generic
second-parameter remainder through eight and, through six, with
the latest stronger scalar remainder.

The exact tangent trace equation gives the unchanged complete seed
polynomial plus the Lie term. Recompute the latter with the new K
array; call their sum Zd. Then

    DZ=Zd+Nnew

bounds W^2(Z_ell-Z0)/ell through eight, since 0<=ell<=1.
The full degree-twelve seed polynomial remains present.

## 5. Refine the anomaly term without dropping contacts

Use the old valid finite metric/inverse arrays M0,M1 and connection
arrays C0,C1. They need not be sharpened to use a smaller DZ.
With Zabs the old seed trace polynomial, the Laplacian difference is

    BX(DZ)=6.1{M0[DZ''+C0 DZ']
          +M1 Zabs(q)''+(M1 C0+M0 C1)Zabs(q)'}.

Thus its scalar-dependent anomaly contribution is

    A(DZ)=emax[BX(DZ)+(11/3)DZ
                       +(11/72)(ZW0+Zabs(q))DZ],
    emax=2000^-10/180.

The remaining old anomaly array is exactly

    Rest=emax[(7/12)DCsq+(11/6)DJsq]+DMU.

Reconcile this with the frozen original DH-A(DZ_old), using the
retained first-parameter curvature formulas. In particular both
inverse-metric contraction variations and the physical magnetic
flux term remain. Then DHnew=A(DZ)+Rest is a complete anomaly
difference bound through six. Replacing only the Laplacian would
miss the other DZ-dependent terms, so replace all of A together.

The new total scalar array is

    P=Nnew+RA2/2+Lie+DHnew,

where RA2 is recomputed with the new tangent parameter bounds.
Take minima with the previous total array through six. Use the
same actual covariant conversion, Jvol<1.001 and full integral
of t^-2 on [1+h,infinity). The norm coefficient is again
(1001/24000)sum H_j. This still bounds only the actual scalar
source, not the quantum/TT or full Schur source.

## Verification and scientific decision

Freeze note, checker and review before execution. Check the complete
source additions, nonconstant Ward recurrences, uniform-input recovery,
all three reconstructions, mass remainder/cubic and tensor conversion.
Reconcile the old anomaly remainder exactly, retain all new product
indices through eight, and independently compare eighth spatial
coefficients of finite warped curvature with the exact scalar remainder
at 60/90 digits. Manufactured fields do not sample the charged source.
Preserve every run and failure; no formal certification is claimed.

Report the actual scalar coefficient and the budget
R0 R2 R4<=1e57 g epsilon_c/(800 C_tau) for using this scalar
allowance at ell=1, before TT. An oversized allowance is not an
actual residual lower bound or a nonexistence theorem. Reassess the
quantitative consequence before selecting further proof work.

Remaining to-do list: assess useful scalar and quantum/TT source
bounds; higher lift/inverse and positive transfer; nonlinear mapping,
preserved neighborhood, signed compatibility and existence; physical
validity, state/boundary, matching, stability, charge/return, opening,
transport, preparation, separate creation and audit.
