# Finite-cut LLL boundary map and the static-state normalization test

2026-10-12. The induced proper-time mouth condition and a prescribed
static optical return delay determine the causal characteristic map for the
free LLL state on the finite product-throat model. This converts the
boundary clock into a conserved bulk JT source and reproduces the MMP static
source when the fixed AdS2 Weyl-anomaly term is retained. No curved RN
exterior or physical stability theorem is claimed.

## Exact finite-cut ray and state map

Let the symmetric throat cut be `sigma=+-s(t)` in
`ds2^2=r_e^2 sec^2(sigma)(-d tau^2+d sigma^2)`, with
`tau'=x(t)>0`. The induced proper-time condition is

```text
cos^2(s)=epsilon^2[x(t)^2-sdot(t)^2].                       (1)
```

A ray traveling from the left cut at time `p_epsilon(t)` to the right cut
at time `t` has constant `u=tau-sigma`, hence

```text
tau(t)-s(t)=tau(p_epsilon(t))+s(p_epsilon(t)).               (2)
```

For channel `j` with fixed exterior optical return time `D_j`, put
`r_j(t)=p_epsilon(t)-D_j`. With a static past, define the unwrapped phase
recursively by

```text
R_j(t)=t,                                      t<=0,
R_j(t)=R_j(r_j(t))+L_j,                       t>0,
L_j=T_epsilon+D_j,
T_epsilon=2 acos(epsilon x_*)/x_*,            (3)
```

where `x_*` is the static throat clock rate. Positive bounded `x(t)` makes
the recursion causal and finite on every finite observation interval, as in
the transparent-return construction. This is an exact ray/state map for
the finite product throat joined to fixed optical-delay channels; it does
not represent propagation through a curved RN exterior.

The outgoing right-mouth global null coordinate is `U(t)=tau(t)-s(t)`;
the reflected left-mouth outgoing coordinate has the same form. The
characteristic state map in global throat time is

```text
Y_j(U(t))=x_* R_j(t).                                      (4)
```

In the static past, `U=x_*t-a`, `a=acos(epsilon x_*)`, and therefore
`Y_j(U)=U+a`: the reference map is a translation. For a complex free
fermion channel with chiral central charge `c_j=1`, its circle-vacuum
outgoing stress in exterior optical time is

```text
O_j(t)=-pi c_j/(12 L_j^2) R_j'(t)^2
       -c_j/(24 pi) Sch(R_j,t).                             (5)
```

Transforming to the throat null coordinate gives the chiral state part
plus the local Weyl anomaly,

```text
T_UU,j(U) = Y_j'(U)^2 t_j^0 - c_j/(24 pi) Sch(Y_j,U) + A_U[sigma]
           = (dt/dU)^2 O_j(t) - c_j/(24 pi) Sch(t,U) + A_U[sigma],
t_j^0 = -pi c_j/[12(x_* L_j)^2].                          (6)
```

Here `A_U[sigma]=-(c_j/(12 pi))[(partial_U sigma)^2-partial_U^2 sigma]`
is the geometric anomaly for `ds2^2=-e^(2 sigma)dU dV`. In global AdS2,
`sigma=log(r_e)-log(cos x)`, so `A_U=c_j/(48 pi)`; the same holds for
`A_V`. The equality follows from the Schwarzian chain rule. The other
chirality is obtained at the left mouth by reflection. The state-difference
tensor is traceless and conserved on the fixed product throat; the local
anomaly is also conserved together with its trace component. The total is
the source for the leading JT constraint. Equation (6) includes both the
coordinate-anomaly term and the fixed-background Weyl anomaly needed to
convert the normal-ordered exterior flux into the covariant throat stress.
The mixed component is fixed by the trace anomaly,
`T_UV,j=-(c_j/(12 pi)) partial_U partial_V sigma`; it is common to the
reference and the finite-history state on the fixed background.

## Linear retarded transfer to the bulk source

Write `tau=x_*[t+psi(t)]`, `phi=psi'`, and set
`a=acos(epsilon x_*)`, `b=sqrt(1-epsilon^2 x_*^2)`,
`r=epsilon/b`, `T=2a/x_*`. To first order (1) gives
`delta s=-epsilon x_*/b phi`. Linearizing (2) yields

```text
delta p_epsilon(t)=psi(t)-psi(t-T)
                     +r[phi(t)+phi(t-T)].                  (7)
```

Linearizing the recurrence (3), with prehistory fixed, gives
`delta R_j(t)-delta R_j(t-L_j)=delta p_epsilon(t)`. For the retarded
Fourier convention `exp(-i omega t)`, `Im omega>0`,

```text
delta R_j/phi = [i(1-exp(i omega T))/omega
                  +r(1+exp(i omega T))]
                 /[1-exp(i omega L_j)].                    (8)
```

At fixed bulk coordinate `U`, the moving endpoint changes its boundary
time by `delta t=-psi-r phi`. Thus the variation of the characteristic map
`Y_j(U)` is

```text
eta_j/phi = x_* [ i/omega (Z_j-E)/(1-Z_j)
                   + r(E+Z_j)/(1-Z_j) ],
E=exp(i omega T),  Z_j=exp(i omega L_j).                    (9)
```

The retarded bulk stress follows by linearizing (6). Since its static
reference is constant in `U`,

```text
delta T_UU,j/phi = [-2i Omega t_j^0
                     -i c_j Omega^3/(24 pi)] eta_j/phi,
Omega=omega/x_*.                                           (10)
```

The reflected chirality supplies `delta T_VV,j`; `delta T_UV,j=0` for the
state difference. Equations (8)-(10) are a closed retarded free-LLL source
for the finite-cut fixed-delay model. They are the input to, rather than a
solution of, the constraint `L_ab[delta Phi]=(2G/r_e^2)delta T_ab`.
The transfer identity (9) is checked at 10,000 finite-cut complex-frequency
samples in `checks/finite_cut_lll_boundary_map.ps1`; this checks the algebra
of the linearized map, not the physical assumptions or the whole stress
tensor. The fixed JT constraint and source balance have separate checks in
the linked notes.

## The Weyl anomaly restores the MMP static source

In the sharp-cut, zero-delay limit `epsilon->0`, `D_j=0`, one has
`x_*L_j=pi`, so the chiral circle-vacuum term is `t_j^0=-c_j/(12 pi)`.
Comparing that term alone with the MMP static LLL source would be
inconsistent: MMP's sourced traceless stress includes the fixed throat
conformal-anomaly contribution. Adding `A_U=c_j/(48 pi)` gives

```text
T_hat_UU=T_hat_VV=-c_j/(16 pi),
T_hat_tau tau=T_hat_xx=-c_j/(8 pi),  T_hat_tau x=0.           (11)
```

Here the hat denotes the traceless two-dimensional part; the common mixed
trace-anomaly component is not included in (11). For `q` channels this is
precisely the static stress in
`lll-state-difference-source.md` and MMP equations (5.26)-(5.27), after
integrating the four-dimensional stress over the sphere. Thus the static
source is reproduced only when the circle Casimir state and its Weyl
anomaly are combined. This confirms the existing MMP normalization; it does
not establish that the dynamical channel return map or its local contact is
the full curved-mouth result.

The finite-cut characteristic map therefore supplies a controlled causal
continuation of the reproduced static source inside this reduced model.
The curved-mouth load, massive angular tower, state-dependent electric
current, absolute reference corrections and renormalized joined contact
remain open. No novelty or stability claim follows from the finite-delay
result alone.

Remaining to-do list: solve and independently check all sourced JT
components with the finite-cut retarded source; match the charged matter
flux through curved mouths; include massive-mode and reference corrections;
then reassess the selected stability result and novelty.
