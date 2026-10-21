# Cold reference overlap and a sufficient full-adjoint error budget

2026-10-20 HST, E003/S167 supporting S139/S074. Ordinary proof before
frozen verification. **The actual scalar reference overlap is positive
uniformly along the finite comparison. Full projection transversality
follows from an explicit adjoint-error budget that remains unproved.**

The large possible common-coordinate phase shift does not destroy
all useful information. On a long cold interval, the retained metric
bound controls the phase difference directly. Its weighted length
is enough to dominate an arbitrary relative phase on the complement.
The actual scalar reference fields, with their regular bolt and full
finite scalar potentials, consequently have a positive overlap.

This is not an overlap theorem for the full coupled adjoints. We
quantify the remaining distinction: a dimensionless weighted error
at most .1 for each full adjoint relative to its own scalar reference
would certify the ORIGINAL fixed projection. Neither error estimate
is supplied by the existing outgoing contraction or asymptotic phase.

## 1. Fixed normalization, geometry and the target norm

Keep the original comparison path, seed projection and unit-amplitude
adjoints of [the finite adjoint construction](finite-comparison-adjoint-and-fixed-projection.md).
Write t=r0/B=1+s, v=1+h, y=sqrt(lambda)rho0 and
z_ell(y)=sqrt(lambda)rho_ell. Thus z_ell'=exp(w), where
w=ell u+S and S=log[1+tau(exp(ell V)-1)], tau=s/t.

Use the common seed-density, two-component weighted inner product

    (f,g)_D=sqrt(lambda) integral t^-4 f dot g dmu0,
    ||f||_D^2=(f,f)_D.                             (1)

It is the scaled dual of the fixed forcing Hilbert norm H0. It
does not impose a new domain or rescale the physical Newton product.
Let phi_ref,ell be the actual scalar reference of the preceding
construction and set Phi_ell=(phi_ref,ell,0) in disk coordinates.
Their physical Liouville amplitudes v_ell satisfy

    v_ell(z_ell)=sqrt(8pi^2 a_ell r_ell^2) phi_ref,ell.

In particular

    (Phi_ell,Phi_0)_D
      =integral t^-4 exp[-(w+2S)/2] v_ell(z_ell) v_0(y)dy,
    ||Phi_ell||_D^2
      =integral t^-4 exp[-w-2S] v_ell(z_ell)^2dy.    (2)

The normalization is the same physical unit cosine amplitude at
each endpoint. No generalized adjoint is asserted to be ordinary L2.

The [metric norm theorem](stationary-completed-metric-response-norm.md)
and [complete equation ledger](comparison-path-complete-equation-residual.md)
give, on h<=s<=.001,

    1/(8*2000^10)<=h<1/2000^10,
    |u|<=4e-6 s Lc(s), |V|<1.2e-8,
    Lc(s)=1+log(.01/s).

The finite completion therefore obeys |w|<4.1e-6 s Lc(s).
Indeed |S|<=1.001 tau |V|<1.202e-8 s. Globally use the already
proved |w|,|S|<3e-8 and |w_y|<.001. These bounds hold on the
original path, not on an assumed nonlinear solver neighborhood.

## 2. Enough weighted cold length and a coherent phase

The exact seed polynomial, with x=lambda B^2>=1000 and |d|<=1, is

    F0(t)=s^2 f(t)-delta/t,
    f(t)=x+2x/t+(3x+1)/t^2+d/t^4,
    delta=v h^2 f(v).

Put P(t)=t s^2 f(t). Its exact derivative is

    P'(t)=s[(t+1)(1+3x(t^2+1))/t^2-d(t-3)/t^4].

For 1<=t<=1.001, the coefficient P'/(2xs) lies between

    a_-=6/(1.001)^2-.001,
    a_+=[2.001(.001+3((1.001)^2+1))+.002]/2.

Since F0=[P(t)-P(v)]/t, a_-/1.001>5.98 and a_+<6.02 give

    5.98(s^2-h^2)<F0/x<6.02(s^2-h^2).             (3)

For s>=.001 the [regional bound](projected-regional-pointwise-bound.md)
gives delta/x<7h^2. Also f/x>=6/t^2: the possible negative d
is covered by the positive 1/t^2 term, and t^2+2t+3>=6.
Thus F0/x>=5.99s^2/t^2 there, with ample thermal margin.

Set s_*=5e-15, y_*=y(1+s_*), a(y)=t(y)^-4, and

    I=integral_0^infinity a(y)dy,
    Ic=integral_0^(y_*) a(y)dy.

Equation (3) yields

    Ic> (1+s_*)^-4 acosh(s_*/h)/sqrt(6.02)>17.5.   (4)

For the last inequality use log(s_* 2000^10)>43,
sqrt(6.02)<2.454 and (1+s_*)^-4>.99999. On the rest of the
exterior, the cold part is below 72/2.445, since
log(2*.001*8*2000^10)<72 and sqrt(5.98)>2.445. The outer part
is below

    (1/sqrt(5.99)) integral_(1.001)^infinity dt/[t^3(t-1)]
      =[log1001-1/1.001-1/(2(1.001)^2)]/sqrt(5.99)<2.22.

Consequently

    17.5<I<31.7.                                  (5)

The exact phase difference is Delta(y)=kbar[z_ell(y)-y], where
4e17<kbar<5e17. On the cold interval use (3) and
|exp(w)-1|<=1.001|w|. Integration by parts gives

    integral_h^s r Lc(r)/sqrt(r^2-h^2)dr
      =Lc(s)sqrt(s^2-h^2)+integral_h^s sqrt(r^2-h^2)/r dr
      <=s[Lc(s)+1].

Since s[Lc(s)+1] is increasing, this proves

    |Delta(y)|
      <(5e17)(1.001)(4.1e-6)s_*[Lc(s_*)+1]/sqrt(5.98)
      <.13, 0<=y<=y_*.                            (6)

Use log(2e12)<28.33 for the last budget. The estimate retains the
bolt-anchored integral; no boundary phase is independently adjusted.
No coherence assumption is made outside this interval.

## 3. Reference fields, averaging and the positive overlap

The actual finite scalar potential has L1 norm below 72000 in
dimensionless proper radius. The same Bessel construction used in
the finite adjoint theorem gives its unit-amplitude reference field
within 2.1e-13 uniformly of

    v_B(z)=sqrt(pi kbar z/2) J0(kbar z), |v_B|<=1.

To see the error, its real resolvent correction is bounded by
eta/(1-eta), eta=72000/kbar<2e-13; unit-amplitude normalization
changes it by at most the square of that allowance. This proof is
for the scalar reference, without a possibly large positive TT block.

For z>0 the far-end Volterra equation for v_B, whose potential is
-1/(4z^2), and |v_B|<=1 give

    |v_B(z)-cos(kbar z-pi/4)|<=1/(4kbar z).         (7)

This follows by integrating the sine kernel from infinity; its
normalization uses the already established Bessel asymptotic. No
asymptotic series is used near the bolt. For y>=1e-8, z_ell>=
exp(-3e-8)y makes the total cosine comparison error below 1e-9.
On the omitted interval 0<=y<=1e-8 use the uniform reference bound
1.000001. Its contribution to each quadratic comparison is below 3e-8.

Here is the needed averaging estimate, including the actual cutoff.
Let b=1, theta or 1-theta, where theta is the ORIGINAL seed quantile
cutoff. It rises monotonically from zero to one, then falls to zero;
its total variation is two. Thus 0<=b<=1 and TV(a b)<=3. For the
sum phase kbar(y+z_ell)-pi/2, its derivative is
kbar(1+exp(w)); for the square phase 2kbar z_ell-pi/2 it is
2kbar exp(w). They are bounded below by 1.99kbar. Integrating
exp(i phase) by parts against a b, the endpoint and variation terms
and integral a|w_y|<.001 I give an absolute bound 10/kbar.
This argument uses bounded variation, not a bound on cutoff derivatives
or numerical values of its quantile radii.

The volume amplitudes in (2) differ from one by less than 1e-7.
Together with (5),(7) and the initial interval, all errors in each
of the following comparisons are less than e=1e-5:

    | ||sqrt(b)Phi_ell||_D^2 - (1/2) integral a b dy |<e,
    | (sqrt(b)Phi_ell,sqrt(b)Phi_0)_D
                         -(1/2) integral a b cos(Delta)dy |<e. (8)

The seed satisfies the same first estimate. In particular every
reference norm is below 4. For b=1 use cos(Delta)>=1-.13^2/2 on
the cold interval and cos(Delta)>=-1 elsewhere. Equations (4),(5)
give the uniform actual reference result

    (Phi_ell,Phi_0)_D
      >[(2-.13^2/2)(17.5)-31.7]/2-e>1.5.         (9)

This bounds reference fields of the actual finite scalar operators,
not just an uncoupled cosine model. The full coupled adjoints are
not substituted into (8) or (9).

## 4. A sufficient budget for the ORIGINAL fixed projection

Define the two currently unbounded errors

    e0=||chi0-Phi0||_D,
    eell=||chi_ell-Phi_ell||_D.                    (10)

The fields here use the actual full coupled unit-amplitude adjoints,
including their returned TT and mass components, in the common seed
disk and density. These errors are finite by the established weighted
domains, but no numerical .1 bound has been proved for either one.

The seed's exact quantile construction gives

    integral (1-theta)t^-4 |chi0|^2 dmu0
                           <=.01 integral t^-4 |chi0|^2 dmu0.

Therefore, with (1) and the triangle inequality,

    ||sqrt(1-theta)Phi0||_D
      <=.1||chi0||_D+e0 <.4+1.1e0.

The first line of (8) for both metrics implies that the squared
finite-reference norm on this same omitted set is at most the seed
one plus 2e. Cauchy--Schwarz and sqrt(A(A+2e))<=A+e then bound
the omitted reference cross-pairing by (.4+1.1e0)^2+e. Combining
this with (9) and the full-field error terms gives

    (theta chi_ell,chi0)_D
      >1.5-(.4+1.1e0)^2-4(e0+eell)-e0 eell-e.     (11)

Thus the concrete sufficient condition

    e0<=.1 AND eell<=.1                            (12)

would imply (theta chi_ell,chi0)_D>.4. This is a conditional theorem,
not a claim that the retained estimates satisfy (12). The more flexible
right side of (11) can be used instead of the symmetric .1 budget.

To translate it back without changing units, put
nu0=||t^-4 chi0||_H0 and p_cut=theta t^-4 chi0/(alpha nu0).
Then sqrt(lambda)nu0^2=||chi0||_D^2 and

    sigma_ell=(theta chi_ell,chi0)_D/[alpha nu0 sqrt(lambda)].

Under (12), ||chi0||_D,||chi_ell||_D<4.1, so

    sigma_ell/sigma0>.4/(4.1)^2>.023,
    ||L_ell^0||_(H0*)/sigma_ell<(4.1)^2/.4<43.     (13)

Here sigma0=nu0. The preceding projected inverse bound would then
be below 44.3 times the still-unquantified physical range inverse norm.
The cutoff, source preparation and physical equation are unchanged.
These conclusions are conditional on (12), unlike (3)-(9).

## 5. What must be bounded next

The full-error target is explicit. With a=a_mix and the actual maps
of the finite adjoint construction,

    chi_ell-Phi_ell
      =[(1/sqrt(1+a^2)-1)phi_ref,ell,0]
       +[(R4_PV VJg_R)/sqrt(1+a^2),0]
       +Rcal_ell(0,-D^-1 B*psi).                  (14)

The first term has norm below 2(.0247)^2<.0013, using
1-1/sqrt(1+a^2)<=a^2/2 and ||Phi_ell||_D<4. A bound .0987
for the norm of the sum of the last two terms would therefore suffice
for (12), at both the seed and the finite comparison.
The interior real-resolvent term and returned TT term require a
weighted pairing estimate. In particular R4_PV includes the negative
elliptic resolvent, and the positive TT elimination in g_R need not
be small. An asymptotic phase or auxiliary series tail alone does not
bound the sum (14) in (1). The same requirement applies at the seed;
its very small asymptotic phase is not silently used as an interior norm.

The next route is to estimate (14), exploiting the scalar reference
equation before absolute bounds on the mixed form. If that budget
cannot be met, evaluate the original signed overlap with controlled
errors. Failure of this sufficient budget would not imply a zero
overlap or a physical nonexistence theorem.

## Frozen verification contract

Freeze the ordinary proof and checker before execution. Reconcile
the finite adjoint, complete-equation, regional metric and scalar
comparison receipts and every shared input hash. Check the exact
thermal polynomial and derivative, rational cold/outer comparisons,
logarithmic bounds, phase integral identity, variation/error budgets,
and the conditional projection constants. Use independent rational
thermal geometries, exact oscillatory integrals with nonconstant
phase, and finite Hilbert-space cutoff examples to check the mechanisms
and distinguish the sufficient condition from its missing hypotheses.
These are supplementary checks, not evaluations of the full adjoint
errors or physical compatibility. The continuum argument remains an
ordinary proof. Retain every scientific run and failed source snapshot.

Remaining to-do list: bound the full scalar/returned-TT adjoint errors
or evaluate the signed overlap; useful full response and cold/higher
control; nonlinear preservation and solve; signed compatibility;
exterior/state, charge/return, opening, stability, transport, preparation,
separate creation and full-result audit.
