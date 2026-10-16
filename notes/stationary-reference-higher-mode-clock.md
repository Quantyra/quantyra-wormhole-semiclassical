# Complete fixed-reference stationary source and clock bound

2026-10-15, E003/S139. Ordinary proof before focused verification.

The retained classical RN/HH reference admits a quantitative bound on the
entire common-renormalized higher angular sector, including the cold bolt
and the actual outer projector. At the fixed benchmark and throughout the
original preparation interval, the result is

```
Theta(t)=b^2 tau_mass(bt), t=r/b, v=r_h/b=1+h,
|Theta'(t)| <= K q^2,             K=10^15,
|Theta(2)| <= K0 q^2,             K0=10^11,
sup |alpha_mass^(1)| < 3.125e-10.
```

The same-background canonical LLL and Weyl contributions give
`sup |alpha_total^(1)| < 3.126e-10`. The superscript means the first
Einstein-constraint response to the complete one-loop source evaluated
on this specified classical geometry. It is not the lapse of a proved
nonlinear semiclassical solution. The previously established failure of
an unrestricted same-regularity Einstein contraction remains intact.

This result includes a trace-value anchor, pointwise individual stress
and radius-source bounds, and an integrated clock consequence. It does
not supply the required metric-response norm, a quantitative finite
vacuum-subtracted mass bound, the full-action Hamiltonian charge,
extremal recentering, or the corrected geometry/state.

## Fixed model and parameters

Retain the single free massless charged four-component Dirac field, the
identical paired stationary clocks and bundles, the actual projectors,
mu0=1/b and zero extra finite curvature coefficients. The sources and
the common four-dimensional subtraction are unchanged. Definitions and
the transformed domain are in
[the conformal-domain note](stationary-exterior-conformal-domain.md).

```
k=2000, q=k^4, eta=1/2000, z_i in [1/4,7/25],
xi=b^2/L4^2, xi(1+3xi)=pi k^6, C=1+6xi,
b^2/G=xi k^10, ell=b/sqrt(C),
F(t)=1+xi t^2-2 a_m/t+Qhat/t^2,
Qhat=1+3xi, a_m=m/b, F(v)=0,
m_n^2=n(n+q)/b^2, d_n=2(q+2n), M_n=m_n ell.
```

In particular `C^2=1+12 pi k^6`, `xi>=k^3`, `xi<2k^3`,
`C/q<.00352`, `M_n>=18 sqrt(n)` and `h<10^-30` at this benchmark.
The elementary bounds used below can be verified with `31/10<pi<22/7`.
These are fixed-k statements; no logarithm is asserted uniform over an
unbounded k family.

The exact temperature relation supplies the indispensable lower bound
on h. Since `kappa0 b/C=2 z_i eta/(xi k^6)`, direct differentiation
of F at v gives

```
kappa0 b/C
 =h [1/v^2+1/v^3+3xi(1+1/v+1/v^2+1/v^3)]/(2C)
 <= h.
```

The bracket decreases from one at v=1. Consequently
`h>=eta/(4k^9)`; the weaker `eta/(8k^9)` already gives
`log(1/h)<80`. The cold proof uses
`L(s)=1+log[(.01+h)/(s+h)]<80`, where `s=t-1`.
An upper bound on h alone would not justify this logarithmic bound.

## Common angular subtraction and the two spatial regions

Let `k2D(tau;x,x)` denote the rank-two orbit heat trace and use physical
proper time in this paragraph. Its local coefficients are

```
k2D=(4pi tau)^-1 [k0+k1 tau+k2 tau^2+...],
k0=2, k1=-R/6, k2=-R^2/120-Box R/60.
H_h(tau)=sum d_n exp(-m_n^2 tau), A(tau)=-H_h'(tau),
zeta_h(z)=sum d_n [n(n+q)]^-z,
H_h(tau)=2b^2/tau+H0+H1 tau+...,
H0=-q-1/3, H1=-zeta_h(-1)/b^2,
Fang=FP_(z=1) zeta_h(z)=4 gamma_E-2 H_q, |Fang|<100.
```

The orbit trace is the massive Green contribution minus the common
four-dimensional heat anomaly, rather than a sum of separately
renormalized two-dimensional determinants. At mu0 b=1 its derivative
and value are exactly

```
partial_y tau_mass
 = integral_0^infty A(tau)
     [partial_y k2D-(partial_y k1+tau partial_y k2)/(4pi)] d tau
   + b^2 Fang partial_y k2/(4pi),

tau_mass
 = integral_0^infty A(tau)
     [k2D-(k0/tau+k1+tau k2)/(4pi)] d tau
   + { (k0/b^2)[zeta_h'(-1)+zeta_h(-1)]
       + b^2 Fang k2 }/(4pi).
```

To track the finite restoration, the regulated mass Green local terms
are `b^(2z-2) k0 zeta_h(z-1)/z`, `b^(2z) k1 zeta_h(z)`, and
`(z+1)b^(2z+2) k2 zeta_h(z+1)`, divided by 4pi and multiplied by
mu0^(2z). The k1 term cancels H0 k1 in the anomaly. The k2 finite
part first contains `b^2(Fang+2) k2`, and the anomaly cancels its
`2b^2 k2`. The a0 anomaly leaves the displayed additional
`zeta_h(-1)` in the value. A different scale would add
`2log(mu b) zeta_h(-1)` to the a0 bracket; none is changed here.

The two detailed proofs are
[the cold-bolt comparison](stationary-reference-cold-bolt.md) and
[the non-bolt stress and anchor estimate](stationary-reference-nonbolt-stress.md).
They retain all Born orders, all angular modes and the actual domain.
Their budgets, established before verification, are:

| Region or contribution | Bound on abs(Theta')/q^2 |
| --- | ---: |
| Cold local analytic remainder | 3.1e6 |
| Cold actual-domain to global-cap localization | 2.3e8 |
| Cold local-flat to global-H2 gluing | 2.2e14 |
| Cold common local restoration | 100 |
| Entire cold patch, h<=s<=.001 | less than 3e14 |
| Entire non-bolt patch, t>=1.001 | less than 3e14 |

Thus K=10^15 is a single conservative bound on the whole exterior.
The cold weighted inverse-variation bounds are estimates for separated
rows with centered weights; they are not small unweighted global
operator derivatives as the observation point approaches the bolt.
Their first/second center variations are used only to differentiate the
fixed-reference localization comparison, not as a Frechet response
bound for a family of unknown metrics.

At the explicit anchor r=2b the same common subtraction gives
`|Theta(2)|<=10^11 q^2`. Its absolute a0 restoration is fixed by

```
zeta_h(-1)=1/30-q^2/6, |zeta_h'(-1)|<1000 q^2.
```

The non-bolt note derives these through Euler--Maclaurin at n=1 and
bounds the remainder, without expanding in q/n outside its convergence
range. Long heat time uses the even RN cap ratio R<=10^10; the larger
ratio of its short-time smooth comparison is not substituted there.

## Stress and radius reconstruction

Use the canonical product allocation of the previous radial-source
note. On the fixed RN reference let `B=F/t^2` and
`I(r)=integral_(r_h)^r B'(R) tau_mass(R) dR`. HH regularity sets the
otherwise independent traceless integration constant to zero. Then

```
p2=I/(2F), rho2=p2-b^2 tau_mass/r^2,
T4=T2/(4pi r^2), J_X,mass=-2G b^2 tau_mass/r^4.
```

Integration by parts gives the particularly useful exact formulas

```
b^2 p2(t) = Theta(t)/(2t^2)
            - [1/(2F(t))] integral_v^t F(s) Theta'(s)/s^2 ds,
b^2 rho2(t) = -Theta(t)/(2t^2)
            - [1/(2F(t))] integral_v^t F(s) Theta'(s)/s^2 ds.
```

F is increasing on this reference: the numerator of F'/2 is
`xi t^4+a_m t-Qhat`, whose derivative is positive and whose horizon
value is positive. The anchor and gradient bound therefore imply

```
|Theta(t)| <= q^2 [K0+K |t-2|],
max(|b^2 rho2|,|b^2 p2|)
 <= q^2 { [K0+K |t-2|]/(2t^2)
          +(K/2)(1/v-1/t) }
 <= (K+K0/2)q^2 <1.1e15 q^2,
|J_X,mass| <= (2G q^2/r^4)[K0+K|t-2|].
```

Here `|t-2|/t^2<=1` and `1/v-1/t<=1` for t>=v>1. The limits
at the bolt are regular: `b^2 rho_h=-Theta(v)/(2v^2)` and
`b^2 p_h=Theta(v)/(2v^2)`. Known Weyl and LLL terms retain their
own exact allocations. The four-dimensional anomaly fixes the
remaining angular stress once these orbit and radius terms are kept.

A derivative bound without the anchor would miss a constant
`tau_mass -> tau_mass+c/b^2`: it changes p2 by `c/(2r^2)`, rho2 by
`-c/(2r^2)` and J_X by `-2Gc/r^4`, while leaving the lapse source
unchanged. This is why the value estimate is part of the theorem. It
does not introduce freedom to add a new coefficient to the fixed model.

The pointwise upper bound above is deliberately coarse at infinity and
cannot be integrated to bound the finite renormalized mass. The earlier
[outer matching theorem](stationary-massive-boundary-layer.md) still
proves the massive/Weyl tail cancellation and finite vacuum-subtracted
geometric mass at each fixed reference. A quantitative bound on that
finite mass, suitable for recentering and the solver, remains separate.

## Integrated first-order clock consequence

The exact static conservation identity gives

```
d_t alpha_mass^(1)
 =-(G/b^2)/(t F(t)^2) integral_v^t F(s) Theta'(s)/s^2 ds,
alpha_mass^(1)(infinity)=0.
```

There is no further 4pi in this formula; the orbit stress already
contains the sphere integral. The horizon root factorization is

```
F(t)=(t-v)/t^2
  [xi t(t^2+vt+v^2)+t-Qhat/v]
 >=xi(t-v)(t-1).
```

For the inequality, subtract the latter polynomial inside the bracket.
The remainder is `xi(v+1)t^2+xi v^2 t+t-Qhat/v`, which is at least
`3xi+1-Qhat=0`. Monotonicity of F now gives

```
|d_t alpha_mass^(1)|
 <= (G K q^2/b^2) (t-v)/(v t^2 F(t))
 <= (G K q^2/(b^2 xi v)) /[t^2(t-1)].
```

The removable horizon limit is controlled by the same inequality.
Integrating from v to infinity yields

```
sup |alpha_mass^(1)|
 <= [G K q^2/(b^2 xi v)] [log(v/(v-1))-1/v]
 <80 K/k^8 =3.125e-10.
```

We used `Gq^2/(b^2 xi)=1/(xi^2 k^2)<=k^-8` and the preceding
temperature lower bound. This includes the cold logarithm rather than
assuming a temperature-independent throat length.

The same classical reference has canonical LLL trace
`Theta_LLL=q[-2+12a_m/t-12Qhat/t^2]/(24pi)`, hence
`|Theta_LLL'|<=qC/4`. Its first-order HH lapse is bounded by
`140/k^9<2.735e-28`. The exact Weyl lapse is

```
alpha_W^(1)=G/(15pi b^2)
  [a_m/t^3+Qhat((3/2)log t-13/8)/t^4],
sup |alpha_W^(1)|<k^-10<10^-33.
```

These are all allocations on the same reference, so their sum gives
`sup |alpha_total^(1)|<3.126e-10`. The isolated nonlinear LLL theorem
is not added to this linear reference sum. Passing the earlier 1e-8
component comparison is not a complete physical clock allowance.

## Verification contract and open scope

The ordinary proof, including the cap interfaces, weighted norms,
angular restoration and value anchor, was reviewed before selecting
focused verification. Verification checks exact algebraic identities,
conservative rational margins, the finite-part and Born arithmetic,
parameter inequalities and source reconstruction. It does not turn
sampling into a uniform heat estimate or compute a corrected metric.
The detailed analytic estimates in the two supporting notes remain
ordinary arguments. No Lean or cloud campaign is selected.

The contract is exact symbolic equality where stated, strict
conservative inequalities for frozen constants, unchanged input hashes
and explicit scope guards. All current source inputs are frozen before
the single checker run. Its receipt and raw log must record the actual
Python -B command, versions, individual outcomes and exit status. Any
failure or later source correction must be retained rather than hidden
by retuning a tolerance.

S139.1-S139.2 and the original S074-S079 gates remain open. In
particular, the compact coordinate variations used in a fixed-reference
localization proof do not resolve the previously proved fourth-order
metric-response obstruction. A higher-order inverse or controlled
spatial order reduction is still required. The reference bound does
not transfer S138 to the unmatched corrected physical state.

Remaining to-do list: quantify the finite mass and horizon recentering;
bound the complete horizon-weighted metric response and justify the
chosen solver; compute the full boundary charge and match the corrected
exterior/state; transfer S138 and complete the driven response, common
capture, tail and renormalized spatial-area work in S074-S079 and S080.
