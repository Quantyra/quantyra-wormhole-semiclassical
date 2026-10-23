# Execute the signed resonant cold moment

2026-10-22, E003/S169. Ordinary argument before supporting verification.
The feasibility audit remains closed; this executes a specific component
of its selected signed matched-response calculation. The full stable
semiclassical wormhole objective remains incomplete.

The massive contribution to the dangerous product slope alpha_q is
EXACTLY ZERO in every angular sector. The remaining full slope is

    alpha_q=-Gq/(6pi B^2),
    -1.012546e-31<alpha_q<-1.012545e-31.                 (1)

This replaces the loose 1.7e-8 allowance for that SIGNED MOMENT, not
the complete finite source or an operator norm. It supplies the signed
first cold shape-residual coefficient. The actual finite boundary return
and weighted higher remainder are not evaluated here.

## 1. Evaluate the massive moment exactly

Use the same exact current product, state, paired normalization and
common angular subtraction as [the full product inverse](current-product-full-metric-inverse.md).
For one positive angular mass M, the retained primary bubble residues are

    tau_j=(2M+2j)(2M+2j+3),
    A_j=(2/pi)(2M+2j+3/2)
          W(j+1)W(2M+j+1)/W(M+j+1)^2,
    W(x)=Gamma(x+1/2)/Gamma(x),
    ell_M=M^2 sum_(j>=0) A_j/tau_j^2-1/(12pi).          (2)

The sum is absolutely convergent, with positive terms O(j^-3). Define
t_j=M^2 A_j/tau_j^2 and t0=t_0. Elementary gamma recurrences give

    t_j/t0=(2M+3/2)_j (M+7/4)_j (3/2)_j (M)_j^2
                /[(M+3/4)_j (2M+1)_j (M+5/2)_j^2 j!]. (3)

The cancellation of the intermediate (M+1)_j and (M+3/2)_j factors
is essential. The resulting sum is the very well-poised 5F4 in
[NIST DLMF 16.4.9](https://dlmf.nist.gov/16.4.E9), with

    a=2M+3/2, b=3/2, c=d=M.

Its convergence hypothesis Re(b+c+d-a)<1 holds with left side zero;
every denominator gamma argument is positive for M>0. Rogers--Dougall
therefore evaluates the hypergeometric sum as

    Gamma(2M+1) Gamma(M+5/2)^2
       /[Gamma(2M+5/2) Gamma(M+1)^2 Gamma(5/2)].        (4)

Meanwhile

    t0=(2M+3/2)/(8pi(M+3/2)^2)
       * Gamma(3/2) Gamma(2M+3/2) Gamma(M+1)^2
            /[Gamma(2M+1) Gamma(M+3/2)^2].             (5)

Multiplying (4)-(5), the factors 2M+3/2 and (M+3/2)^2 cancel.
Gamma(3/2)/Gamma(5/2)=2/3 leaves exactly 1/(12pi). Hence

    ell_M=0 for every M>0.                             (6)

This is an exact evaluation of the original pole sum and its already
specified two-dimensional anomaly. It does not add a finite counterterm
per mode. The original common a2 restoration and reduction are proportional
to p^2 and do not contribute to the slope. The retained locally uniform
summability of the subtracted response licenses its derivative at p=0;
there is no summation of unrenormalized operator norms. Thus

    alpha_q=(2G/B^2)[-q/(12pi)+sum_n d_n ell_Mn]
           =-Gq/(6piB^2).

The massless LLL remains as written; it is not obtained by declaring the
zero mass sector to satisfy the massive resolvent hypotheses. All massive
modes and their nonzero higher response remain in H(p).

The numerical enclosure (1) uses the exact product's squared-sphere shift
<5.5e-23 about the current RN reference, not an unvalidated quadrature of
its pressure. The reference formulas are

    g0=2000^-10, g+11g^2/(360pi)=g0,
    F=1-384pi*1e22*g, C_RN=sqrt(1+12pi*g*2000^16/F),
    GY_RN=g/[(C_RN-1)/6],
    Y/Y_RN in [1/(1+5.5e-23),1/(1-5.5e-23)].

No scale, physical field or norm has changed. The exact signed evaluation
is the improvement, not the number of digits in (1).

## 2. The physical cold root with the evaluated slope

Retain the COMPLETE current matrix and coefficients from the previous
product proof, including the local R^2 term and the Stieltjes remainder:

    F(p)=p[B(p-2)+a(p-2)H(p)]+a(p-2)alpha_q,
    a(ell)=9beta ell^2+c ell+d,
    B(ell)=4b^2-beta d+beta ell(12b-c).

On |p|<=2.05e-31 the previous bounds |a|<8.04, |B|>3.99 and
|H|<.667 eta_G still apply. Now |alpha_q|<1.013e-31, so

    8.04*1.013e-31+2.05e-31*8.04*.667*8.5e-33
                                <2.05e-31*3.99.

Rouche's theorem gives exactly one non-Ward root in this disk. It is
real, negative, and therefore has nu>1 on p=2-nu(nu+1). The sign also
follows directly: a(0)<0 and alpha_q<0 make F(0)>0, while F is negative
at the negative real radius by the same strict bound. Consequently

    -2.05e-31<p_root<0, 0<nu-1<6.84e-32,
    exp(84|nu-1|)-1<5.75e-30.                          (7)

This is a product analytic-continuation root, not an allowed finite
exterior zero mode, a retarded spectrum or a proof of a weighted inverse.
The complete-H2 L2 inverse <.51 remains valid. No new global domain is
inferred from (7).

## 3. Apply the signed moment to the seed's first cold jet

Let alpha_g=Y/X and d_g=C-1-6lambda B^2. Retain the exact finite seed
curvature from the product-centered construction:

    B^2 R(t)=-12lambda B^2+d_g(4t^-5-6t^-6), t=r/B.

Its first variation about t=1 is delta R=16d_g Y s, s=t-1.
The leading finite-temperature cold radius is S=s=h cosh chi. In the
product normal model P2 S=0. For its canonical shape T, the exact
linearized curvature identity is

    delta R/X=(6ell-2R/X)S+2P2 T,
    P2 T=chi_g S, chi_g=4+2alpha_g+8d_g alpha_g.         (8)

Thus the growing orbit correction is generally a GENERALIZED Ward field,
not just a multiple of cosh chi: P2 T is nonzero, P2^2 T=0. The first
cold log profile arises from this resonance. Gauge/outer normalization
can add a multiple of cosh chi to T without changing (8).

Apply the full invariant matrix to this formal product first jet. Its
shape row is exactly

    F_T,1=chi_g alpha_q S.                              (9)

Every massive contribution to this coefficient is zero by (6); the
p^2 H and local -beta p^2 terms annihilate this generalized Ward field.
The mixed row -2p j S also vanishes. Analytic functional calculus on
these growing probes is justified by the massive decay with M>18 and
the same cutoff Ward limit as the retained product proof. This is not
an application of the L2 inverse to a growing source.

The current bounds alpha_g<2.1e-11, |d_g alpha_g|<1.71e-29 give

    chi_g<4.000000001, |F_T,1|<4.053e-31 |S|.           (10)

The trace/conformal row is retained too. With epsilon0=(lambda0-lambda)/X,
the exact coefficient after the local cancellations is

    F_S,1/S=-8d_g alpha_g-48epsilon0
       +chi_g GX(6+5alpha_g)/(90pi)-11GX alpha_g/(45pi)
       +48beta d_g alpha_g.                            (11)

Here the apparently O(beta) terms cancel using (8). The actual guards
epsilon0<2e-36, GX<8.5e-33, beta<1.18e-8 imply |F_S,1|<1.4e-28|S|.
Neither row is set to zero. All signs, product contacts and the selected
large local coefficient are retained.

Equations (8)-(11) are the homogeneous normal-model first cold jet. The
fixed-chi expansion is not asserted uniform out to the finite seam
chi approximately log(1/h). It does not bound the second-order geometric
remainder, the actual finite-state kernel difference, or boundary data.
This distinction prevents a formal first jet from becoming a claimed
finite residual certificate.

For scale only, insert (10) into the already executed CLASSICAL positive
neck problem with zero radius seam value. Its constant-envelope gain is
B_L/6<85/6, with the physical bolt factor 1+h<1.001. This ONE SHAPE
component would give radius <5.8e-30. It is an executed comparison of
the evaluated moment, not the complete quantum response; it excludes
the other row, higher remainder and physical boundary contribution.

The previous unsigned envelope gain >301.92 remains correct for its
stated comparison. It is simply much too pessimistic for this particular
signed first-jet component. It remains unsafe to apply the improvement
to the entire finite source or to claim a unit physical cold inverse.

## 4. Boundary priority and what remains of the selected calculation

The exact coefficient (1) removes an unquantified bulk slope from the
selected signed-response calculation. The formal leading cold forcing
no longer provides evidence for the large envelope accumulation. This
supports pursuing the same fixed construction, without a theory switch.

The SINGLE selected calculation remains the complete matched finite-seed
bolt-radius functional and a certified weighted remainder, tested against
the 2e-9 cold allocation. Only its resonant product moment and normal-model
first jet have been executed here. The remaining work concerns the actual
finite boundary map and remainder, not another generic source bound.

In particular the boundary term cannot be chosen to cancel the bulk
moment. Fixed clock/temperature, correlated mass, the auxiliary scalar
fast phase and the actual paired quantum domain must determine it. The
old classical constraint reconstruction applies to a conserved joint
effective source; treating the off-shell quantum variation as separately
conserved would omit metric/density contacts. The current complete
auxiliary rows and their constraint identity must be kept in the coupled
realization. A product Ward root alone supplies no global compatibility
or projection theorem on this finite domain.

The weighted C5 tail and actual boundary-to-bolt return remain unbounded;
an N2 estimate has no generic raw-C5 upgrade. No full metric inverse
defect <=.25, full Y0<=5e-9 or nonlinear Z2<=5e7 has been established.
All preserved-ball, regional, outer-transfer, matching, validity and
stability requirements remain as recorded in the closed audit. In
particular (10) does not certify any of those missing constants.

Decision: use the exact massive cancellation and signed LLL moment in
the already selected boundary-matched response calculation. Prioritize
the prescribed finite boundary realization and its weighted remainder.
Do not reopen the unsigned estimate chain or label this partial execution
as the complete signed response. No existence/nonexistence conclusion
or completed S169 component is claimed.

Remaining to-do list: complete signed matched finite-seed bolt-radius
functional and weighted remainder; physical domain/compatibility and
nonlinear admission; actual solution, matching, semiclassical validity,
all-sector stability, charge/return, opening, transport, preparation,
separate creation and final audit.
