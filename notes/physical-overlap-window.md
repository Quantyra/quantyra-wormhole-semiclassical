# Scale compatibility of the finite-collar estimate with the MMP throat

2026-10-11. This is a parameter-window check, not a matched-geometry proof.
It tests whether the finite flat point-dipole collar estimate can sit outside
an MMP mouth while retaining the large-separation throat regime and its
semiclassical bound.

## Scales

For the collar result in [mouth-collar-scattering.md](mouth-collar-scattering.md),
choose

```text
epsilon_1 = 1/12500, epsilon_0 = 1/25000,
rho_i = epsilon_i d, D = d-2 rho_1,
3 <= q <= 10^12, 0.05 <= Xi=xi d <= 5.
```

The collar theorem then bounds its prepared-input error by `8e-6`. The
exterior comparison window `0.05 <= xi D <= 4` lies inside this `Xi` window,
because `D/d=1-2/12500`.

MMP gives the extremal mouth radius and the short-distance throat scale as

```text
r_e = sqrt(pi) q l_p/g,
d_0 = 16 pi r_e^3/(q G_N)
    = 16 pi^(5/2) q^2 l_p/g^3,
G_N=l_p^2.
```

Their large-separation solution uses `d >> d_0`; semiclassical control of
the throat requires `ell << q^3 l_p`. The geometric-mean choice

```text
d_* = sqrt(d_0 q^3 l_p)
```

balances these two inequalities:

```text
d_*/d_0 = sqrt(q g^3/(16 pi^(5/2))),
d_*/(q^3 l_p) = sqrt(16 pi^(5/2)/(g^3 q)),
rho_0/r_e = (1/25000) 4 pi^(3/4) q^(3/2)/sqrt(g).
```

Thus, for fixed weak coupling, the hierarchy
`d_0 << d_* << q^3 l_p` and `r_e << rho_0 < rho_1 << d_*` is parametrically
compatible once `q` is large compared with `16 pi^(5/2)/g^3`, while still
below the collar theorem's finite cap `10^12`. This removes a scaling
conflict between the collar estimate and the large-separation MMP regime.
It does not establish that a specific value of `q` and the running coupling
meet every desired numerical margin or ultraviolet assumption.

## Scope of the conclusion

The collar estimate is for the flat two-monopole Dirac operator with
independently specified inner and outer incoming data. The scale hierarchy
does not bound the Reissner-Nordström-to-monopole metric and spinor trace
mismatch at `rho_0`, nor does it control the feedback when the outer port is
attached to the full exterior. It therefore supplies no physical return
operator, full stress tensor, or travel-time theorem. The next concrete
physical gate is a quantitative mouth trace/curvature bound at `rho_0`,
followed by feedback composition with the exterior resolvent.

The MMP scale formulas are from [§§2.2, 5.3-5.5](https://arxiv.org/html/1807.04726v3).
