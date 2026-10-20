# Global fast bounds for the scalar component of the coupled inverse

2026-10-20 HST, E003/S164 supporting S139/S074. Ordinary proof and
protocol before execution. **A global component estimate; full
fourth/log inverse and physical completion remain open.**

Let A0 be the retained analytic-seed Hessian, including its mass column,
and let u=G0 f be its unique physical fast solution for f in Ran(A0).
The [existing inverse theorem](stationary-scalar-outgoing-inverse.md)
establishes this compatible inverse qualitatively. Here we quantify
its scalar component while retaining the complete shape coupling.

Use the exact seed York map u=Rcal(tau,c) from the
[TT inverse](stationary-weighted-tt-inverse.md). Put

    g=G lambda, W=(r/B)^2,
    Tbar=(-Box)/lambda, gradbar=lambda^-1/2 grad,
    (fbar_tau,fbar_c)=lambda^-2 Rcal* f,
    Fsrc=||W fbar_tau||_L2(dmu)+||W fbar_c||_L2(dmu),
    dmu=4pi beta r^2 dr.

The main result, uniform on the retained preparation interval, is

    ||W tau||             <=4.4e7 g Fsrc,
    ||W Tbar tau||        <=3.2e9 sqrt(g) Fsrc,
    ||W gradbar tau||     <=4e8 g^(3/4) Fsrc.          (1)

Since g<2000^-10, the three numerical coefficients are respectively
below 4.4e-26, 1.1e-7 and 8e-17. These multiply Fsrc, not an
unconverted comparison-path residual.

The weight is genuinely fast: for z=B/r,
||W v||_dmu=sqrt(4pi beta B^3)||z^(-7/2)v||_L2(dz/z).
The output controls a value, gradient and second-order scalar graph
norm on the entire exterior, including the regular bolt. It is not
the fourth-order norm of tau, the fourth/log norm of the whole metric,
a bound on c's mass coefficient, or a nonlinear solution theorem.
In particular tau is the York scalar; the disk metric also contains
the shape reconstruction. The transformed source norm is explicit,
and no unproved numerical norm for Rcal* is supplied.

## 1. Dimensionless operators and the strengthened vertex weight

Rescale proper distance to rhohat=sqrt(lambda)rho and use Tbar.
Write E=E_sp/lambda and k=sqrt(E-9/4). In this note k is the
dimensionless quantum frequency, not the preparation integer 2000.
The accepted bounds give

    E=60pi/g+11/3, k^2=60pi/g+17/12,
    E>180/g, 13<k sqrt(g)<14, k>4e17,
    rhohat(2B)<172, ||U_geo||_L1(d rhohat)/k<2e-14.  (2)

All scalar and auxiliary Hilbert norms retain the same physical
volume normalization. The rescaling of derivative and vertex
operators below is explicit, so no period factor is dropped.

Use Dbar=D/lambda^2, Bbar=Bmix/lambda^2 and

    V=Fcal/lambda=(sqrt(80pi^2) Bbar Dbar^-1/2, Ubar_s*),
    J=diag(I,-J_s), ||J||<=1.

Then the exact scalar Schur operator is

    Pbar=(Tbar+4)(Tbar-E)-V J V*.                    (3)

The original outgoing proof used w=B/r. Its vertex estimates actually
permit the stronger w=(B/r)^2=W^-1:

    V y=g_y+divbar f_y,
    ||W g_y||<=1000/sqrt(g) ||y||,
    ||W f_y||<=3 ||y||.                            (4)

For the TT vertex, the earlier g/f bounds have a factor (B/r)^3,
so this strengthening costs nothing: the constants are
201/sqrt(g) and .001. For the local scalar vertex,
|a_s|<=5lambda(B/r)^4 and
|Z|<=2|d|lambda(B/r)^4/x, x=lambda B^2. Thus its squared
dimensionless g constant is at most 40pi|d|/(g x), and its f
constant is sqrt(5). Since |d|/x<1.3e-22, combining the two
vertices in their actual auxiliary direct sum gives (4).
No derivative of a square-root coefficient or of Dbar^-1/2 is taken.
All original TT, York-return, trace-density and boundary contacts remain.

## 2. A compatible weak-source Helmholtz bound through the bolt

For a source p=a+divbar b, write A=||W a|| and Bv=||W b||.
The moment pairing with the regular reference field is the physical
weak pairing, including -integral gradbar(phi_ref).b.
If <phi_ref,p>=0, the actual reference outgoing inverse
R+=(Tbar-E-i0)^-1 satisfies

    ||W R+ p|| <=347 A/k+700 Bv.                   (5)

Here is a direct proof with the original bolt realization. Under the
unitary Liouville map the reference operator is

    H_B0+U_geo-k^2, H_B0=-partial_rhohat^2-1/(4rhohat^2).

For 0<rhohat<=sigma define the terminal Bessel kernel

    K_B(rhohat,sigma)=(pi/2)sqrt(rhohat sigma)
      [Y0(k rhohat)J0(k sigma)-J0(k rhohat)Y0(k sigma)].

It has K_B(sigma,sigma)=0 and partial_rhohat K_B(sigma,sigma)=1.
Thus integral_rhohat^infinity K_B p has the correct source sign
for -partial^2-1/(4rhohat^2)-k^2. The sign follows from the
[Bessel Wronskian](https://dlmf.nist.gov/10.5.E2).
The [derivative identities](https://dlmf.nist.gov/10.6) and the
[Nicholson bounds](https://dlmf.nist.gov/10.9.E30) already used
in the outgoing proof imply

    |K_B|<=2/k,
    |(partial_sigma-1/(2sigma))K_B|<3.              (6)

For completeness, normalize sqrt(pi v/2)J0(v), Y0(v); each has
absolute value at most one. For k sigma>=1, the two differentiated
products in (6) each cost at most sqrt(3/2). For k sigma<=1, use
|J0|<=1, |J1(v)|<=v/2 and
|H1(v)|^2<=2/(pi v)+4/(pi^2 v^2). The costs are at most
sqrt(pi/2+1) and sqrt(pi/2)/2, whose sum is below 3.
This proof includes points arbitrarily close to the bolt.

The actual terminal kernel solves the Volterra equation

    K=K_B-integral_rhohat^sigma K_B(rhohat,t)U_geo(t)K(t,sigma)dt.

Its ordered series gives |K|<=2 exp(2eta)/k and the input derivative
bound 3 exp(2eta), where eta=||U_geo||_1/k<2e-14.
There is no differentiated endpoint term because K(sigma,sigma)=0.
The physical input derivative after Liouville conjugation differs
from partial_sigma-1/(2sigma) by a multiplier of norm below 4.
Consequently

    |K|<2.001/k, |D_sigma K|<4.                    (7)

Let j_+ be the actual outgoing field and phi_ref its regular,
unit-amplitude real partner. The usual Green formula contains an
integral from the bolt to rhohat. Precisely when <phi_ref,p>=0,
replace that integral by minus its remaining tail. The result is

    U R+ p(rhohat)=integral_rhohat^infinity
           [K(rhohat,sigma) U a(sigma)
                    -D_sigma K(rhohat,sigma) U b(sigma)]d sigma. (8)

This is the actual outgoing inverse with its regular bolt domain.
For a source with nonzero moment, the terminal formula generally
has the forbidden logarithmic bolt term and is not a physical inverse.
Compatibility is essential, not a free boundary choice.

For r>=2B, F>=lambda r^2/2. Therefore for sigma>=rhohat
the weight ratio W(rhohat)/W(sigma) is at most one, and beyond
rhohat(2B) it decays at least as exp[-sqrt(2)(sigma-rhohat)].
Both Schur integrals of this triangular ratio are at most
172+1/sqrt(2)<173. Applying (7) gives (5).
Also ||w||_L2(d rhohat)^2<173 for w=W^-1.

Initially take smooth g/vector sources supported off the endpoints.
Compatible approximants are obtained by a compact scalar correction
with nonzero reference moment. The weighted weak pairing is continuous.
The triangular bounds and the previously proved outgoing weak-source
limit then extend (5) while retaining the regular bolt realization.
No input or output derivative singularity at the bolt is discarded.

## 3. Positive elliptic factor in the same fast weight

Let R_e=(Tbar+4)^-1. For every p=a+divbar b in the same weighted
weak-source class,

    ||W R_e p||<=80 A+320 Bv.                      (9)

The seed has F_rr>=lambda and r'>=0. The dimensionless orbit
scale satisfies a_orb''/a_orb>=1/2, so
A_orb>=coth(rhohat/sqrt(2))/sqrt(2).
Integrating (a_orb r^2 v^2)' and using Cauchy--Schwarz gives

    ||gradbar v||^2>=||v||^2/8                     (10)

on the regular-bolt scalar form domain. The bolt boundary term
vanishes; no Dirichlet value of v there is imposed.

The seed polynomial also gives F<=lambda r^2+1 for r>=r_h:
with u=B/r<=1 its remainder is
(1-u)^2+x u(-4+3u)-delta u+d u^2(1-u)^2<=1.
Use the Lipschitz weight W_c=max(1,(r/(2B))^2), so
W_c<=W<=4W_c and
|gradbar log W_c|^2<=4+1/x, x>1000.
For v=R_e p set A1=||gradbar(W_c v)|| and B1=||W_c v||.
The weighted form identity, (10) and x>1000 yield

    A1^2+4B1^2-||v gradbar W_c||^2
      >=A1^2-B1^2/x
      >=A1^2/2+B1^2/20
      >=(A1^2+B1^2)/20.

The source pairing is at most
(A+4Bv)sqrt(A1^2+B1^2), since |gradbar log W_c|<3.
Thus sqrt(A1^2+B1^2)<=20(A+4Bv), proving (9).
Apply the calculation first to truncated weights min(W_c,R) and
then take their monotone limit. It holds for weak sources by form
approximation; no unproved weighted spectral inverse is assumed.

## 4. Auxiliary response and the full compatible scalar equation

Set R4+=(R+-R_e)/(E+4). The exact partial-fraction identity is

    Tbar R4+=(E R++4R_e)/(E+4).

For a source satisfying the reference moment, (5) and (9) imply

    ||W R4+ p||<= (81 A+1100 Bv)/E,
    ||W Tbar R4+ p||<=348 A/k+701 Bv.              (11)

The fourth-order bolt domain is the same one used by the earlier
outgoing construction. The distributional copies of p cancel in
the second identity. Weak sources are permitted.

The full Schur solve retains the auxiliary variable y=V* tau:

    M y=V* R4+ p, M=I+Acal J-Lcal J,
    ||M^-1||<=2/(1-8e-12)<2.001,
    tau=R4+(p+VJy).                               (12)

We now quantify y for fast input p, without assuming p itself has
zero reference moment. Use (4), ||w||_2^2<173, the existing outgoing
kernel and derivative estimates, and the identity contact in the
two-derivative kernel. They give

    ||V* R+ p||
      <=[173173/(k sqrt(g))+1557] A
         +[519000/sqrt(g)+3+1557k] Bv.

For the elliptic part, the unweighted form inequalities are
||V*(Tbar+4)^-1/2||<=500/sqrt(g)+3 and
||(Tbar+4)^-1/2p||<=A/2+Bv. Combining both parts in (12),
using E>180/g and 13<k sqrt(g)<14, gives

    ||y||<=3sqrt(g) A+6100sqrt(g) Bv.              (13)

This step uses the outgoing estimate on arbitrary sources and does
not prematurely impose the reference source condition.

Now let p be the actual Schur forcing for a physical compatible
coupled source. The total forcing q=p+VJy has a representation with

    A_q<=3001 A+6100000 Bv,
    B_q<=9sqrt(g) A+1.001 Bv.                     (14)

It is q, not necessarily p, that obeys <phi_ref,q>=0. Indeed
q=(Tbar+4)(Tbar-E)tau, and the physical fast solution has both
regular bolt factors and zero oscillatory tail. Integration by
parts against phi_ref has vanishing endpoint terms. This is also
the exact outgoing-amplitude identity in the
[phase theorem](stationary-scalar-phase-moment.md).
It preserves the full coupled moment rather than substituting
the reference moment for it.

Apply (11) to q. The explicit budgets give

    ||W tau||<=1400 g A+2800000 g Bv,
    ||W Tbar tau||<=100000 sqrt(g) A+1000 Bv.       (15)

This controls the actual coupled scalar component, with the full
nonlocal TT block still inside V and M.

## 5. Physical forcing and the claimed constants

The exact transformed forcing is

    p=80pi^2[fbar_tau-Bbar Dbar^-1 fbar_c].

The actual critical TT gap is
Dbar>=.99/(8pi g), hence
||Dbar^-1/2||<=sqrt(8pi g/.99).
Using the separate TT vertex constants in section 1,

    A<=800||W fbar_tau||+30000||W fbar_c||,
    Bv<.2sqrt(g)||W fbar_c||.                     (16)

Substitution into (15) proves the first two bounds in (1).
For the gradient use the scalar weighted form identity with W:

    ||W gradbar tau||
       <=sqrt(||W Tbar tau|| ||W tau||)+2sqrt(20)||W tau||,

because |gradbar log W|^2<=20. The first two estimates,
g^(1/4)<6e-9 and
sqrt(4.4e7*3.2e9)+2sqrt(20)*4.4e7*6e-9<4e8
give the third. Form cutoffs justify the identity on the same
fast domain. The finite positive temperature is retained throughout.

## Meaning for the remaining inverse problem

The [S163 fourth/log lower bound](weighted-normal-inverse-scale.md)
remains intact. Equation (1) uses fewer output derivatives and an
explicit York-transformed source norm. It supplies a global fast
scalar estimate across the bolt, not a numerical upper bound for
the full fourth/log metric inverse.

The remaining linear work is to control the mass-extended shape
response in its fast fourth/log norm, recover the remaining scalar
derivatives with matched quantitative source estimates, and bound
the source/reconstruction conversions. A small low-order scalar
response alone does not control an unknown metric's higher jets.
The nonlinear family and the signed moment of the actual preparation
source remain unproved. No comparison-path defect is converted into
a physical solution error by this theorem.

## Frozen verification contract

Before execution freeze the note and checker. Revalidate the five
predecessor receipts and their exact current inputs. Verify the
weighted form and partial-fraction identities, all rational budgets,
the source/vertex normalizations and the retained benchmark bounds.

At 80 digits test the terminal Bessel kernel and its input derivative
for k=1,10,1e4 and dimensionless pairs
(1e-8,2e-8),(.1,.5),(.5,2),(2,3),(10,40),(100,100).
Check the derivative jump and the source sign. Identity tolerance is
1e-45; kernel inequalities have no added tolerance.

Use the manufactured regular Bessel field v=sqrt(rho)exp(-rho^2).
Its compatible source is sqrt(rho)(4-4rho^2-k^2)exp(-rho^2).
Compare its direct and weak-source terminal reconstructions at
k=1,3,10 and rho=.01,.3,1,2, using a=-k^2 sqrt(rho)exp(-rho^2)
and b=2rho^(3/2)exp(-rho^2) in the Liouville representation.
Integrate to 12 with the analytic Gaussian tail allowance below 1e-59;
require relative agreement within 1e-40. Check the zero reference
moment as well.

As a failure witness, at k=1 use the source sqrt(rho)exp(-rho^2),
whose moment is nonzero. Its terminal field has logarithmic bolt
slope exp(-1/4)/2. At rho=1e-8 and 1e-10, verify that slope to
relative tolerance 1e-10. This detects the invalid use of a terminal
fast solution without imposing regular-bolt compatibility.

These are manufactured kernel diagnostics, not physical preparation
sources. They do not certify the continuum theorem by sampling.
Preserve all runs and failed sources. Global kernel localization,
weighted form limits and the coupled reconstruction remain ordinary
proofs. No nonlinear solution or physical stability is certified.

Remaining to-do list: full global fourth/log inverse and nonlinear
family; signed compatibility and coupled exterior/state, charge and
return; opening, stability, transport, preparation, creation and
full-result audit.
