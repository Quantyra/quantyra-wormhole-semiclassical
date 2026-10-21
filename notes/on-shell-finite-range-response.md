# Full finite-comparison response in value and first-derivative norms

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary proof
before frozen verification. **Quantitative full low-order response
at the actual finite comparison; no nonlinear solution or cold/higher
neighborhood preservation is asserted.**

Use the original action, finite comparison, physical state/domain,
clock, temperature, flux, mass column and fixed seed projection.
Let B_ell=Jvol A_ell be the actual fixed-density linearized operator,
and let f_ell be the retained complete comparison residual. The unique
physical solution of

    Pi_cut B_ell h_ell=-Pi_cut f_ell

exists by [fixed-projection transversality](finite-trace-gradient-and-projection-transversality.md).
Write h_ell=Rcal_ell(tau,c)=(S,T), gamma0=4pi beta B^3 and
gradbar_ell=lambda^-1/2 grad_ell. Then

    ||(S,T)||_ell/sqrt(gamma0)<1.1e-7 ell^2,
    ||gradbar_ell(S,T)||_ell/sqrt(gamma0)<2.2e-6 ell^2. (1)

These are actual-volume unweighted integral norms of BOTH disk
coordinates, including the returned TT field and its correlated mass
tail. They concern the finite Hessian, not the seed Hessian applied
at a new endpoint. No new numerical bound for the mass coefficient
itself, fast weighted H4/log norm, or pointwise cold error follows.

## 1. Operators, source norm and actual positive elimination

At a fixed ell use dimensionless derivatives and operators

    g=G lambda, Tbar=-Box_ell/lambda, H0=Tbar+4,
    E=60pi/g+11/3, e=E+4, k^2=E-9/4,
    Dbar=D_ell/lambda^2, Bbar=Bmix,ell/lambda^2,
    V=(V_T,V_s)=(sqrt(80pi^2)Bbar Dbar^-1/2,Ubar_s*),
    J=diag(I,-J_s), Pbar=H0(Tbar-E)-VJV*.

All Hilbert norms and transposes use dmu_ell; the unitary Liouville
map for dimensionless proper radius includes its lambda^-1/4 factor.
Keep W=t^2 and omega=t^-1/2. For a compatible full source f define

    (f_tau,f_c)=lambda^-2 Rcal_ell* f,
    F=||W f_tau||_ell+||f_c||_ell.                (2)

The TT source is unweighted in (2); no weighted norm of its transpose
is assumed. The established finite vertex coefficients are

    ||omega^-1 f_i[y]||<=a_i||y||,
    ||omega^-1 g_i[y]||<=b_i g^-1/2||y||,
    V_i y=divbar f_i[y]+g_i[y],
    a=(4e6,4), b=(1000,2e-5).

Also e>180/g, sqrt(E)<15/sqrt(g), 13<k sqrt(g)<14,
sqrt(g)<4e-17. The realized TT form has its .90 gap and
||Dbar^-1/2||<5.3sqrt(g). The positive local auxiliary block
is below b0=2.6e-13.

Let R+=(Tbar-E-i0)^-1, R4=(R+-H0^-1)/e,
A=V*H0^-1V/e, L=V*R+V/e and N=(I+AJ)^-1. The retained
block norm majorants, valid also for complex L, are

    Nmaj=[[1.000001,2.6e-7],[2.6e-7,1.000001]],
    Lmaj=[[.02405,2.405e-8],[2.405e-8,2.405e-14]].

Consequently (I-NLJ)^-1 has the entrywise majorant

    Hmaj=[[1.026,3e-8],[4e-8,1.000001]].          (3)

The matrices bound block operator norms; they do not replace the
actual operators by commuting scalar matrices.

Put p0=80pi^2 f_tau, r=sqrt(80pi^2)Dbar^-1/2 f_c,
y=V*tau and w=y-(r,0). The exact full equations give

    q=p0+VJw, tau=R4 q,
    (I+AJ-LJ)w=V*R4 p0-(r,0),
    c=-Dbar^-1/2 w_T/sqrt(80pi^2).               (4)

In particular ||r||<150sqrt(g)||f_c||. This choice of w retains
the original TT forcing inside the positive elimination.

For completeness set B_T=V_T/sqrt(e), B_s=V_s/sqrt(e) and

    K=H0+B_T B_T*-B_s J_s B_s*
       >=(1-b0)H0+B_T B_T*.

These are the same closed form maps used in the seed low-response
proof; the finite local constant is b0. Factorization gives

    ||K^-1||<=1/[4(1-b0)],
    ||K^-1 B_T||<=1/[2sqrt(1-b0)],
    ||K^-1 B_s||<=sqrt(b0)/[2(1-b0)],
    ||H0^1/2 K^-1||<=1/[2(1-b0)],
    ||H0^1/2 K^-1 B_T||<=1/sqrt(1-b0),
    ||H0^1/2 K^-1 B_s||<=sqrt(b0)/(1-b0).        (5)

The identity N V*H0^-1=V*K^-1 bounds the negative elliptic
part of the right side of (4), after division by e, by
31sqrt(g)||W f_tau|| in the TT channel and
1.6e-5 sqrt(g)||W f_tau|| in the local channel.

For the outgoing part use integral omega^2 dz<200 and
integral W^-2 dz<32, whose L2-weight product is below 80,
and the actual kernels |G|<1.001/k, |D G|<3. Before multiplying
by N its channel coefficients are bounded by

    80[1.001 b_i/13+3a_i](800g/180)||W f_tau||.

After adding -N(r,0), the two N-preconditioned source bounds are
150.001sqrt(g)F and 3.91e-5 sqrt(g)F. Applying (3) gives

    ||w_T||<155sqrt(g)F, ||w_s||<5e-5 sqrt(g)F.  (6)

No smallness of the TT diagonal of A is used.

## 2. Cancel the on-shell TT vertex before bounding the radiation

The [adjoint cancellation](on-shell-adjoint-error-and-trace-gradient.md)
has an operator version. For h the TT lift of Dbar^-1/2 v, let
r_h=-Jbar:h and kappa=sqrt(5)/(120pi). Then

    V_T=(Tbar-E)C+Q_T, C v=kappa r_h.            (7)

The entire (a-epsilon T)r_h term has been factored, including its
classical part. The remainder keeps the curvature, magnetic and
off-shell return terms. Expand the latter as
divbar(Ebar Ybar)=Ebar divbar Ybar+gradbar Ebar.Ybar,
rather than keeping its large derivative coefficient in a generic
weak-source estimate.

The preceding trace-gradient proof implies the pointwise bound
|gradbar Ebar|<4.1e-5 omega: its leading term is 4e-5/t,
and the anomaly derivative costs at most g(1e28+1e17)/180.
Also |Ebar|<3.2e-12/t. The finite TT energy/return bounds give

    ||C v||<1.34sqrt(g)||v||,
    ||gradbar C v||<1.7e14 sqrt(g)||v||,
    Q_T v=divbar F_v+G_v,
    ||omega^-1 F_v||<561sqrt(g)||v||,
    ||omega^-1 G_v||<.0013 g^-1/2||v||.         (8)

Here ||h||,||gradbar h||,||Ybar||,||divbar Ybar|| have
coefficients 25,420,25,160 in units sqrt(g E_TT), and
E_TT[Dbar^-1/2 v]<=||v||^2/.90. The factor
sqrt(5)/(2sqrt(.90)) is below 1.2. Thus the last coefficient,
after multiplication by sqrt(g), is at most

    (1.2/180)[g(6e15*420+2e20*25)+20.02*25*(4e-27)]
       +1.2[160(3.2e-12)+25(4.1e-5)]<.0013.

The derivative-source coefficient is (1.2/180)(200)(420)<561.
The C bounds use |Jbar|<8/t and |gradbar Jbar|<1e15 omega.
The latter large derivative coefficient will be multiplied by g,
not used as an unscaled response error.

Identity (7) is first an exact distributional identity on regular
TT inputs. C maps the TT energy space into scalar H1 with decaying
outer weight. It has zero reference radiation moment under Tbar-E.
On physical fast/mass inputs its bolt parity and outer tails make
the boundary pairings zero; energy approximation extends the identity.
The local term C is retained when applying R+ to (7).

## 3. A compatible terminal-kernel estimate and the full scalar field

For a reference-compatible weak source a+divbar b, put
Aq=||omega^-1 a|| and Bq=||omega^-1 b||. The actual regular
outgoing inverse satisfies

    ||R+(a+divbar b)||<300 Aq/k+1000 Bq.         (9)

Reuse the terminal Bessel/Volterra argument of the
[compatible scalar theorem](global-compatible-scalar-response.md),
now with the finite L1 potential allowance 72000. It gives
|K_terminal|<2.001/k and |D_input K_terminal|<4, since
2||q_geo||1/k<4e-13 and |D-D_Bessel|<5. The input boundary
term is zero because K_terminal(z,z)=0. Compatibility replaces
the bolt-to-z source integral by its remaining tail; it is not a
new boundary choice.

The triangular kernel weighted on its input by omega has squared
Hilbert--Schmidt norm bounded by integral z omega(z)^2 dz.
Using z_40<178 and dz<=1.5dt/t outside forty gives

    integral z/t dz
      <178^2/2+(1.5/40)(178+1.5)<16000<127^2.

Thus 2.001(127)<300 and 4(127)<1000 prove (9), through the
regular bolt and the outer end. For physical compatible full
sources, q in (4) has zero reference moment. Neither p0 alone
nor an uncoupled forcing is assigned that condition.

Let z=R+q=H0 tau and z1=z-Cw_T. By (7), z1 is the reference
outgoing solution for

    q1=p0+Q_T w_T-V_s J_s w_s.

Its reference moment is also zero, since the removed (Tbar-E)Cw_T
has zero moment. Equations (6),(8) give

    Aq1<801F, Bq1<.00021 sqrt(g)F.

Equation (9) therefore gives ||z1||<19000sqrt(g)F.
Testing (Tbar-E)z1=q1 in its regular form realization yields

    ||gradbar z1||<=sqrt(E)||z1||+sqrt(Aq1||z1||)+Bq1
                  <286000F.

This weak-source form estimate needs no derivative of q1. A positive
shifted scalar resolvent identifies z1 in H1 from its L2 solution
and H^-1 forcing, as in the seed low-response proof. Cutoff/form
approximation supplies the integrations by parts. Restoring Cw_T
with (6),(8) gives

    ||z||<20000sqrt(g)F, ||gradbar z||<290000F.   (10)

Now keep BOTH pieces of R4 and use the exact form identity

    p=p0-V_T r,
    tau=z/e-K^-1 p/e-K^-1 VJV*z/e^2.             (11)

It follows from K tau=(H0z-p)/e, without commuting K and H0.
The original finite vertex coefficients and (10) give
||V_T*z||<1.2e12 F and ||V_s*z||<1.2e6 F. Apply (5) to
(11), including the V_T r term, to obtain

    ||tau||<6gF, ||gradbar tau||<1700gF.         (12)

For the value estimate, the p0 coefficient is at most
800/[4(180)(1-b0)] and the r coefficient below
1.001(150)/(2*13); their maximum is below 5.78. The z/vertex
terms cost less than 1e-7 in units gF. For the derivative estimate
the leading coefficient is 290000/180, plus at most
1.001(150)/13 for the TT datum, 1.001(800)/(2*180) for the
scalar datum and a vertex remainder below 1e-7. These fit (12).

## 4. The actual returned TT field and residual composition

Equation (4) and the critical gap give

    ||c||<31gF, E_TT[c]<=||w_T||^2/(.90*80pi^2).

The actual finite disk return satisfies ||Rcal(0,c)||<25||c||.
Its derivative bound is .8||h'_rad||+73||h||, bounded by
2161sqrt(g E_TT[c]). Consequently

    ||Rcal(0,c)||<775gF,
    ||gradbar Rcal(0,c)||<14000gF.

Combining with (12) proves the uniform finite range-response bounds

    ||Rcal(tau,c)||<800gF,
    ||gradbar Rcal(tau,c)||<16000gF.             (13)

The critical solution equals the physical mass-extended solution
by the established uniqueness theorem. No mass column or return
tail was projected away to obtain (13).

Return now to fixed seed-density forcing f. Since Rcal(tau,0)=(tau,0)
and the unweighted finite TT return has norm below 25,

    ||W f_tau||+||f_c||<27lambda^-2 Fdisk,0(f)

after applying Jvol^-1 BEFORE the actual York transpose. The actual
to seed volume cost is below 1.000001, so 26.000026<27 suffices.
This uses no weighted estimate for the TT transpose.

The fixed-projection theorem supplies the range source

    f_comp=f-p_cut L_ell^0(f)/sigma_ell,
    Fdisk,0(f_comp)<22.5 Fdisk,0(f).

Its solve has exactly the original Pi_cut right-hand side. Apply
the [retained complete-residual bound](projected-low-correction-and-mass.md)

    lambda^-2 gamma0^-1/2 Fdisk,0(f_ell)<2.2e-13 ell^2/g.

Equations (13) now give 800(27)(22.5)(2.2e-13)<1.1e-7
and 16000(27)(22.5)(2.2e-13)<2.2e-6, proving (1).
The signed residual is still the actual projection direction times
its unevaluated physical coefficient. A compatible projected solve
does not solve the omitted physical compatibility equation.

## Verification and next decision

Freeze this proof and checker. Reconcile the finite projection,
on-shell, mixed-form and earlier low-response inputs. Check the
operator split, positive elimination, terminal-kernel moment, rational
channel/response budgets and exact source conversion. Include an
independent noncommuting finite-block solve, compatible terminal
integration and omission controls. Preserve every run and failure.
Continuum claims remain ordinary arguments; no Lean run is selected.

This supplies a useful FULL low-order response bound at the finite
comparison. It is not the full higher/logarithmic solver norm. An H1
integral bound does not control a cold point value, a fourth derivative
or its Holder modulus. Next obtain the needed local/cold and higher
response bounds, then preserve the nonlinear neighborhood and certify
the projected nonlinear solve. S167.1-S167.4 remain open.

Remaining to-do list: cold/local and higher response control; nonlinear
mapping, preservation and solve; signed physical compatibility;
exterior/state, charge/return, opening, stability, transport,
preparation, separate creation and full-result audit.
