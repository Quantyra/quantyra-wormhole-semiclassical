# Complete physical shape forcing on the finite interior

2026-10-22, E003/S167. Ordinary argument before supporting verification.
Keep the original comparison, bookkeeping action, HH paired domain,
full angular tower, common subtraction, physical charge and fixed return.
Write lambda=lambda_e, x=lambda B^2, g=G lambda and gamma0=4pi beta B^3.
Let I be the fixed region from the regular bolt through t_seed=40.

For the ACTUAL complete shape covector b_src of the
[York source split](finite-higher-york-source.md), define the restricted norm

    F_b4,I=gamma0^(-1/2) sum_(j=0)^4
                 ||1_I nabla_bar,ell^j b_src||_L2(dmu_ell).

The indicator restricts each already differentiated tensor; it is not
differentiated. The finite ledger below gives an exact rational C_I with

    F_b4,I <= C_I ell^2/g, C_I<1e110, 0<=ell<=1.       (1)

This includes the complete canonical massive/LLL, Weyl and classical
contributions, the true physical volume and every four-dimensional
connection slot. The outer region and therefore the global F_b4 remain
unbounded numerically here. No useful solve target or stability gate passes.

## 1. The orbit trace gives the shape row directly

In the common disk let sigma=sigma0+ell u be the canonical orbit factor,
phi=log t_ell and X=B^2 exp(2phi), with

    t_ell=1+(t_seed-1)exp(ell V),
    g2=A^2 exp(2phi+2sigma)(dx^2+dy^2),
    V_ell(X)=(1-ell)V0(X)+ell V_phys(X).

The [exact spherical action split](stationary-centered-covariant-seed.md)
has Euclidean local action

    I_pot=-(4G)^-1 integral sqrt(g2)
                    [X R2+|dX|^2/(2X)+V_ell(X)].

Vary the canonical orbit conformal factor at fixed X. In two dimensions
the kinetic density is conformally invariant, and integration by parts
in the curvature term gives

    F_T,pot=(Box2 X-V_ell)/(8piG X).

The complete canonical Lorentzian scalar is S_can=S_M+S_LLL, normalized
by 2G/(B^2 C). The independent Weyl orbit trace in the
[full source allocation](comparison-path-complete-source-response.md) is
S_W=-eta_G(M_t+M_r)/pi. Its fixed-radius variation has the same scalar
normalization. Dividing by physical four-volume and then by lambda^2 gives

    t=lambda^-2 F_T=-O_ell/(8pi g),
    O_ell=(C/x){z_ell^2[V_ell/C-2 Delta_gamma phi-4|dphi|_gamma^2]
                         +ell z_ell^4(S_can+S_W)},
    z_ell=exp(-phi), eta_G=GC/B^2.                    (2)

Here gamma is the ACTUAL canonical orbit metric exp(2sigma)delta.
Indeed Box2 X=C[2 Delta_gamma phi+4|dphi|_gamma^2]. The canonical
Euclidean action source has the opposite sign to the Lorentzian trace;
this is the sign convention already fixed in the interior source theorem.

For s=f_tau=-E_ell/(8pi g), the exact shape identity is

    b_src=(t-s/3)/2=-O_ell/(16pi g)-s/6.             (3)

No pressure reconstruction is needed for (2). The physical charge Q
remains in V_phys; the comparison charge in V0 is not substituted for it.
The seed and actual complete tangent equations give O_0=dot O_0=0.
The latter is the original first-order reconstruction, not a new
compatibility assumption. These cancellations apply to the entire row.

## 2. Extend the finite massive estimate through the bolt

The [middle finite-source proof](finite-middle-source-four.md) needs only
its stated normal-chart, low-norm, C9 and actual-domain hypotheses.
The cold seed supplies the same normal charts. Its Goursat equations are

    sigma_zw=b_zz(z_radius) exp(2sigma)/8,
    (z_radius)_zw=b_z(z_radius) exp(2sigma)/4.

The [cold normal-chart construction](stationary-reference-cold-bolt.md),
with the [centered coefficient guards](seed-cold-local-eight-jets.md),
has |b_z|<.07, |b_zz|<2.31, |b_zzz|<10 and initial
|p|<=.501s<=.000501. On radius 1/32 the Picard bounds give

    |delta z_radius|<.00005, |sigma|<.002,
    |partial sigma|<.067,

with contraction row sums below .02. These include p=0 at the bolt.
The origin is an ordinary disk point, and no polar 1/rho coefficient
or inverse temperature is used. The finite fixed-chart P5 bound and
global C9 input already cover the cold region. Thus the SAME fixed-cap,
translation and gluing proof applies on all of I.

For the residual estimate we also need a complex path neighborhood.
For each real a in [0,1], allow |ell-a|<=1, hence |ell|<=2. The fixed
metric direction has P5<=1.3e-4, so the doubled background and the
translation polynomial still lie strictly inside the .011 cone.
The cap's fourth translation allowance remains 1e88; its earlier
bound has ample slack after doubling the finite direction.

For complex ell the local massive coefficient differs from its real-a
value by less than 3e-8. The real energy lower constant .98 therefore
becomes at least (sqrt(.98)-3e-8)^2>.97. The old 1.1 energy recurrence
factor still applies. The derivative reserves A_j=.1 j!300^j and the
physical endpoint rows remain valid after replacing .00013 by .00026.
The fixed-domain inverse bounds are the same coefficient-cone bounds
used by S155/S156; the physical boundary is not moved.

Consequently, in each fixed seed chart, the entire massive source is
holomorphic in ell and has coordinate tensor bounds

    |partial^j S_M,ell|<=m_j,
    m_0=1e-8+Kcap+G0,
    m_j=2^(ceil(j/2)) j^j (3e-15)(1e28)^j+G_j, 1<=j<=4, (4)

where Kcap and G_j are the complete cap/gluing bounds in the preceding
receipt. The extra 1e-8 bounds the fixed flat-vacuum scalar using the
retained seed value and nonconstant source bound. It is retained in m0
and in the physical factor z_ell^4; it has not been assigned zero value.
Uniform angular convergence justifies the spatial and path derivatives.

## 3. Finite local metric jets on the same complex disks

Use the improved full D=z0 nabla_gamma,seed divided tangent arrays U,Vbar
from [the low reconstruction](finite-low-gap-reconstruction.md). They
already include 100^j. For either field, start the coordinate-conversion
recurrence of the preceding proof with

    B(n,0)=n!400^n field_n,
    Gamma_i=1e4 i!1000^i.

Through order eight, its exact positive recurrence gives
|partial^j(u,V)|/j! <=(1e7)^j, j>=1. The analytic seed log radius and
normal factor fit the same reserve. The point values retain their
independent small bounds; high coefficients are never used as values.

Here is an explicit finite majorant for the nonlinear radius. Put
h=sum_(j=1)^8(1e7 X)^j. The divided nonconstant coefficients of
tau=1-1/t_seed are bounded by exp(h)-1, with tau0<=1. For |ell|<=2,
|exp(ell V0)-1|<3e-8. Define the nonconstant positive series

    Aplus=3e-8[exp(h)-1]+1.001 exp(h)[exp(2h)-1],
    Phi_plus=h-log[1-Aplus/(1-3e-8)].               (5)

It bounds phi's derivatives. Those of sigma are bounded by 3h.
Exact finite coefficient comparison through eight proves that both
are bounded by P=sum_(j=1)^8 L^j X^j, L=1e9. Also |phi0|<4,
|sigma0,total|<.003 and |z_ell|^m<1.001 for m<=6. All arguments of
the radius logarithm stay away from zero on these complex disks.

## 4. Every local term and physical normalization

All following arrays are divided coordinate tensor jets. Products,
derivatives and exponentials mean finite coefficient operations only.
Set H=1.01 exp(2P), E_m=1.001 exp(mP). The canonical curvature has
majorant 4H P''. Its complete LLL scalar is

    S_LLL=eta_G q R_gamma/(12pi C),

whose coefficient is below 1e-30. This keeps the zero angular level
as well as the massive tower in (4).

The fixed potential guards are 1/C<1e-10, x/C<1/5.9,
lambda0 B^2/C<1, Q/(CB^2)<2 and |d|/C<1e-24; the physical-charge
bound is proved in the [complete compact action](complete-compact-disk-graph.md).
After multiplying by z_ell^2, valid majorants for the two potentials are

    P0=2e-10 E2+2+2E4+2e-24(E4+4E5+3E6),
    Pphys=2e-10 E2+6+4E4,
    Ppot=3P0+2Pphys.                                (6)

The exact kinetic terms in (2) have majorant

    Kkin=E2 H[4P''+4(P')^2].                       (7)

For the Weyl orbit trace use the FULL two-variable polynomials of S160,
M_t+M_r=2M_t-(M_t-M_r). They have total derivative order at most four,
at most one undifferentiated phi and H-power at most two. Differentiate
these actual polynomials through four more spatial derivatives, keeping
the magnetic monomial qbar phi/3 with qbar<160000. With |phi|<4 and
|H|^2<2, their full coordinate-tensor coefficients obey the reserve

    w_j=1e10(4+j)! L^(4+j)/j!, j<=4.                (8)

The checker expands every x/y derivative and verifies this bound from
the actual monomials. It does not differentiate an unsigned source bound.
The Weyl scalar array is (eta_G/3)w, using pi>3.

Let M_j=m_j/j! from (4), and retain all coefficients through four in

    Oplus=6.1{Ppot+Kkin
             +2 E4[M+1e-30(4H P'')+(eta_G/3)w]}.    (9)

This includes the physical potential, canonical massive and LLL,
Weyl trace, radius multipliers and |ell|<=2. The flat vacuum appears
inside M and is multiplied by the actual radius. No physical sector
or density factor is omitted or counted twice.

## 5. Quadratic residual, full physical tensors and volume

Use holomorphy in the unit ell disk around every real a. Cauchy and
O0=dot O0=0 give, before changing coordinates or taking unsigned sums,

    |partial^j O_ell| <=j! [X^j]Oplus *ell^2.       (10)

Indeed the second path derivative is at most twice that bound, and
the real Taylor integral has factor ell^2/2. The metric connection
at the final real ell is applied AFTER this fixed-chart cancellation.

Convert (10) to seed canonical covariant tensors with the raw connection
reserves (.3,4000,2e6) and the full ordered-tensor recurrence of the
preceding proof. The same tensors on the canonical product have zero
sphere components because its factors form a product. Divide order j
by j!, and multiply by z0^j<=1 to obtain divided D arrays H_j.

Now use the actual finite physical connection array C from the
[refreshed local geometry](finite-refreshed-local-mixed.md), with

    H_(0,k)=H_k,
    H_(n+1,k)=(k+1)H_(n,k+1)
                    +n sum_(i=0)^k C_i H_(n,k-i),
    P_n=2(2.5)^n H_(n,0), n<=4.                    (11)

This is the full four-dimensional conversion, including all angular
slots and the zeroth connection term. It uses the actual metric and
sqrt(C/x)<2.5. It is not a radial-derivative proxy for physical H4.

The physical volume in the common seed coordinates is Jvol dmu0,
Jvol=exp(2ell u+4log(t_ell/t_seed))<1.001. Therefore

    Vol_ell(I)/gamma0 <1.001 integral_1^40 t^2 dt <150^2. (12)

The physical period cancels in this ratio. No inverse power of the
small positive thermal gap is introduced at the bolt.

Let C_sigma4 be the retained sum of the scalar residual H4 coefficients
from the [direct radius-source proof](finite-radius-operator-source.md),
using the smaller value allowance 6e-12. Its weighted global norm
bounds its unweighted restriction. Equations (3),(10)-(12) give

    C_I=(25/8)sum_(n=0)^4 P_n+C_sigma4/6 <1e110,    (13)

which proves (1). The receipt retains the exact rational C_I and its
upward decimal enclosure. This is an interior contribution; adding
an unproved outer estimate or declaring global F_b4 finite numerically
would exceed its scope. Cold-weighted nonlinear-family preservation
also remains distinct from this unweighted physical source norm.

## 6. Verification contract and remaining work

Freeze note, checker and review before execution. Verify the cold
Goursat margins, complex path cone and coercivity, all finite local
jet/potential/Weyl bounds, exact source normalization, quadratic
Taylor identity and actual physical tensor/volume conversion. Check
the orbit trace independently from four-dimensional Christoffel/Ricci
components and from the reduced conformal Euler variation. Keep the
original scalar receipt and every historical input unchanged.
Use independent manufactured regular metrics for the orbit normalization
and real/complex Cauchy-Taylor residual diagnostics; these are not
actual quantum-state samples. Preserve every run and source snapshot.

Remaining to-do list: outer complete physical shape H4 source and global
assembly; higher TT inverse and positive transfers; nonlinear mapping,
cold/outer preservation, existence and signed compatibility; physical
matching, stability, charge/return, opening, transport, preparation and creation.
