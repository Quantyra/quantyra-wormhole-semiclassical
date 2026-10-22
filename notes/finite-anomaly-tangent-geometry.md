# Finite geometry from the sharper tangent and complete anomaly terms

2026-10-22, E003/S167. Ordinary argument before supporting verification.
Keep the original physical path, state/domain, full angular tower,
subtraction, mass and fixed return. The preceding scalar bound retained
old finite geometry in its anomaly terms. Refine that geometry using
the already proved tangent arrays, then recompute every term of the
same scalar residual. This does not bound the quantum/TT source.

## 1. A finite improvement of valid spatial derivative bounds

Use the fixed canonical reference D=z0 nabla_gamma and W=t0^2 of
the [local coefficient proof](finite-higher-local-coefficients.md).
The [per-order reconstruction](finite-source-order-scalar.md) supplies
divided ordered D-jet bounds U,K,Vbar,KV,h through ten. U bounds
both original tangents u,nu without W; K bounds them with W. The
last three bound V, W V and W tau D^j V, where nu=tau V.
All weights stay outside differentiation at the observation point.

For 0<=a<=1 let

    S_a=log(1+tau(exp(a V)-1)),
    partial_a S_a=nu exp(a V-S_a), S_0=0.

The old finite bound on a u,S_a has constant 3e-8 and positive
derivative polynomial Tplus. Put Vplus=Vbar-Vbar_0. If Tplus is
any valid bound on the positive spatial derivatives of S_a, then
the ordered tensor chain/product rules give

    S1u=1.001 exp(Vplus+Tplus) U.

This bounds every divided D-jet of partial_a S_a without W,
uniformly in a. The zeroth exponential is below 1.001 by the
retained small values of V and S_a; high derivatives are not
exponentiated as values. Integration from zero to a commutes with
these fixed reference derivatives. Thus for each j>=1,

    Tnext_j=min(Tplus_j,max(U_j,(S1u)_j))

is another simultaneous bound on a u,S_a. Its constant is zero;
the existing finite value bound 3e-8 remains separate. Starting with
the old Tplus, apply this operation exactly six times through ten.
Every step follows from the preceding proved bound. No convergence,
fixed-point solution or implicit regularity premise is needed.
Coefficientwise minima preserve all earlier estimates.

## 2. Rebuild the finite geometry and its first path variation

With the resulting Tplus retain the original exact parameter bounds

    S1=1.001 exp(Vplus+Tplus) K,
    S2=1.001 q exp(Vplus+2Tplus) K KV,
    M0=3.003 exp(2L+4Tplus),
    M1=2 M0(K+S1),
    C0=(3/2)M0 M0',
    C1=(3/2)(M1 M0'+M0 M1').

Here q=1+Zeta and L=-log(1-Zeta) are the old seed majorants.
M0 bounds both metric and inverse jets, M1 their W-weighted first
path derivatives, and C0,C1 the reference connection difference
and its W-weighted first path derivative. Primes are formal spatial
derivatives. All these are the original formulas with sharper inputs.
Use only metric order ten, connection order nine and curvature
order eight from these new arrays.

Let Rref=B''(q)+2e-10. In the same normalization define

    R0=Rref+2C0'+2C0^2, R1=2C1'+4C0 C1,
    Ric0=2R0, Ric1=2R1,
    Sc0=M0 Ric0, Sc1=M1 Ric0+M0 Ric1,
    J1=6.1[ Ric1+(Sc1 M0+Sc0 M1)/4 ],
    CW1=6.1[ M1 R0+M0 R1+2(M1 Ric0+M0 Ric1)
                       +(2/3)(Sc1 M0^2+2Sc0 M0 M1) ].

The reference curvature has zero path derivative. Integrate these
uniform first-variation bounds along the original path, retaining
the explicit seed arrays pJ,pC,pZ from the local coefficient proof:

    JW0=J1+13pJ, CW0=CW1+25pC, ZW0=6.1 Sc1+pZ.

These bound the finite W-weighted covariant components. Take
minima with the previously proved component arrays. A second valid
bound is ZW0<=pZ+DZ when DZ bounds the W^2-weighted scalar
difference divided by a: a<=1 and W>=1 allow that enlargement.
No seed constant curvature is separately assigned a decaying weight.

## 3. Complete scalar and anomaly differences

Recompute the six exact scalar remainder products from the preceding
two scalar refinements, now with the sharper finite Tplus. Retain
the same tangent arrays, the full degree-twelve seed polynomial,
the cold factor on high gap derivatives and every differentiated
tau term. Order eight uses tangent order ten. Take minima with
the preceding scalar remainder Nold through eight. Call the result
N, so DZ=Zd+N bounds W^2(Z_a-Z0)/a, with Zd the full seed
polynomial plus the unchanged Lie bound.

Use the complete anomaly formulas, now with the new finite arrays:

    BX=6.1{M0[DZ''+C0 DZ']
              +M1 Zabs(q)''+(M1 C0+M0 C1)Zabs(q)'},
    DCsq=2 M0^4 CW0 CW1+4 M0^3 M1 CW0^2,
    DJsq=2 M0^2 JW0 J1+2 M0 M1 JW0^2,
    MU=(4e-27/18)*1.001 exp(4Tplus) q^4,
    DMU=4 MU S1,
    DH=emax[BX+(11/3)DZ+(11/72)(ZW0+Zabs(q))DZ
                                 +(7/12)DCsq+(11/6)DJsq]+DMU,
    emax=2000^-10/180.

The four Weyl inverse-metric variations and both Ricci contraction
variations remain. For a W^2 bound, products of two W-weighted
curvature components already carry the two required weights;
an extra W on M1 can be discarded since W>=1. The magnetic
term uses the fixed physical flux. The reference charge is not
substituted for it. No unsigned allowance is subtracted as if it
were a signed physical cancellation.

Recompute the same second-parameter reference-potential array RA2
with the new Tplus,S1,S2. The complete residual array is

    P=N+RA2/2+Lie+DH.

Take minima with the preceding P through six. Convert to actual
covariant derivatives with the new valid C0 using the unchanged
ordered tensor recurrence. Keep Jvol<1.001, the full noncompact
integral and coefficient C_tau=(1001/24000)sum H_j. The exact
quadratic path zero and t^-4 decay are unchanged. This controls
the scalar residual only; all original stronger low contracts remain.

## 4. Independent checks and quantitative decision

Freeze this note, checker and review before execution. Reproduce
the old metric/connection parameter arrays and old retained anomaly
remainder exactly when the old inputs are restored. Check each
finite improvement step, all positive arrays and all contractions.

For an independent finite-log check use polynomial V,nu=tau V,u
on |x|<=1/2. With V=(2-x+3x^2+x^3)/1e5 and tau=tau0+x/10,
tau0=.2,.8, the complex modulus of V is below 3.375e-5 and of
tau below .85. Thus |log(1+tau(exp(a V)-1))|<.001 for a in
[0,1]. Cauchy gives initial divided derivatives .001*2^j.
Apply the same finite improvement and compare all eight actual
derivatives at a=.3,.9 and 60/90 digits. This manufactured test
does not sample the actual charged source.

Independently use a nonconstant warped product

    ds^2=exp(2A)(dx^2+dy^2)+exp(2B) ds_sphere^2,
    A=phi+a u+S_a, B=phi+S_a.

Its sectional curvatures are

    k01=-exp(-2A) A'',
    k02=k03=-exp(-2A)(B''+B'^2-A'B'),
    k12=k13=-exp(-2A)A'B',
    k23=exp(-2B)-exp(-2A)B'^2.

Form Ricci, scalar and Weyl from these sectionals. Compare the
path derivative of their invariant squared norms against the
explicit covariant-component contraction derivative, including
every inverse metric slot, at spatial orders zero, three and six.
Omitting the inverse-metric contributions must fail. Use both
precisions and four nonconstant path/background configurations.
This checks contraction algebra independently of the positive
majorant implementation, not a full physical solution.

Report the resulting C_tau and rounded scalar-only transfer budget
1e57/(800 C_tau). Identify the newly limiting allowance from the
complete ledger before selecting more proof work. No useful
transfer pass, nonlinear existence or stability result is presumed.

Remaining to-do list: useful scalar and quantum/TT source bounds;
higher lift/inverse and positive transfer; nonlinear mapping,
preserved neighborhood, signed compatibility and existence; physical
validity, state/boundary, matching, stability, charge/return, opening,
transport, preparation, separate creation and audit.
