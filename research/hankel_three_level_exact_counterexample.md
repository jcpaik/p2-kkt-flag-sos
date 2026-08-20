# Exact obstruction to the three-level filtered-bosonic relaxation

Let `h` be the normalized monomial basis of `Sym^3(R^3)`, ordered by

`(003),(012),(021),(030),(102),(111),(120),(201),(210),(300)`.

Let `A` be the exact invertible filter from the symmetric cubic feature to
the tangent Pluecker feature, as written in
`hankel_three_level_exact_counterexample.py`.  Thus `G=A H A^T`.

In the standard basis `e_0,...,e_9`, put

```
u0  = (e4+e6)/sqrt(2),
v0  = (u0+(7 sqrt(6)/2)e9)/sqrt(149/2),
u1c = (sqrt(3)e0+e2)/2,       u1s = (sqrt(3)e3+e1)/2,
v1c = (3u1c-e7)/sqrt(10),     v1s = (3u1s-e8)/sqrt(10),
v3c = (e0-sqrt(3)e2)/2,       v3s = (e3-sqrt(3)e1)/2.
```

Define

```
H = (3129/6260) v0 v0^T
  + (7671/31300)(v1c v1c^T+v1s v1s^T)
  + (1/200)(v3c v3c^T+v3s v3s^T).
```

This is manifestly positive semidefinite of rank five.  Its coefficients
were selected so that the two affine normalizations hold exactly:

```
tr H = tr G = 1,        tr(P1 G)=1/5,
```

where `P1` projects onto the spin-one summand of
`wedge^2(H_2)=H_1 direct-sum H_3`.

For the ordinary bosonic reductions `R2=Tr_1 H`, `M=Tr_1 R2`, and
`L=2 P_{H2} R2 P_{H2}`, direct exact calculation gives

```
||G||^2 =     647167327 / 1959380000,
||L||^2 =   13019087509 / 13776890625,
||M||^2 =   19113855963 / 48984500000,
```

and therefore

```
||G||^2 - (1/2)||L||^2 + (1/6)||M||^2 + 1/18
  = -6353145083 / 293907000000 < 0.
```

So PSD at all three levels plus the exact homogeneous block-trace relation
does not prove the desired purity inequality.

The missing condition is genuine full Hankel coherence.  For example,
both `(5,5)` and `(4,6)` encode the total multi-index `(2,2,2)`, so a fully
Hankel matrix must obey `H_55/6=H_46/3`, equivalently `H_55=2H_46`.
Here

```
H_55=0,       H_46=21/6260,
```

so that identity is violated explicitly.  Thus this is a counterexample
only to the larger filtered-bosonic relaxation, not to the original fully
Hankel cone.
