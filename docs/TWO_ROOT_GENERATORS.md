# Two-root generators: the modulated families $\hat C_n$ and $\hat S_n$

Expository companion to [E1 admissible](E1_ADMISSIBLE.md) §3.  Every
displayed identity is machine-checked in exact rational arithmetic by
`python3 solve_e1.py` (54 checks); the relevant check names are quoted
in place.

## 1. What a two-root block is

A **two-root block** is a square in which two points ("roots") are
sampled from $\mu$ and held as parameters, and the squared quantity is a
$\mu$-average over a third point (the "leaf"):

$$\mathcal Q[A]\;=\;\iint_{S^2\times S^2}
\Big(\int_{S^2} A(t_1,t_2,s)\,d\mu(y)\Big)^{2}
\,d\mu(x_1)\,d\mu(x_2)\;\ge\;0,
\qquad
\begin{aligned}
t_1&=x_1\cdot y,\\ t_2&=x_2\cdot y,\\ s&=x_1\cdot x_2,
\end{aligned}$$

for a polynomial $A(t_1,t_2,s)$.  Expanding the square produces a
linear combination of four-point moments of $\mu$ — the labels of the
SDP.  In `sos_search.py` these are the `two_root_*` sectors (basis
monomials $t_1^it_2^js^k$, Gram matrices over them); the
orientation-odd sectors use leaves $\det(x_1,x_2,y)\,B(t_1,t_2,s)$, and
the `_minor` sectors carry the extra root weight $1-s^2$.  The one-root
blocks are the special case where $A$ does not depend on $x_2$.

## 2. The coplanar chart and what (E1) allows

If a sharp certificate exists, every square in it must vanish at every
zero-energy measure $\mu^*=\frac13\delta_{\pm e}+\frac23\nu$
($\nu$ on $C=e^\perp$, $\hat\nu(2)=\hat\nu(6)=0$): this is condition
(E1) of [Limit certificate](LIMIT_CERTIFICATE.md).  For a two-root
block it reads

$$\int A(x_1\cdot y,\;x_2\cdot y,\;x_1\cdot x_2)\,d\mu^*(y)=0
\qquad\text{for all roots } x_1,x_2\in\operatorname{supp}\mu^* .$$

The decisive configuration puts both roots on the equator.  Place $x_2$
at angle $\delta$ from $x_1$ along $C$ and the leaf $y$ at angle
$\theta$ from $x_1$; then

$$t_1=\cos\theta,\qquad t_2=\cos(\theta-\delta),\qquad s=\cos\delta .$$

Because $\nu$ is free except for $\hat\nu(2)=\hat\nu(6)=0$, the
condition forces the restriction of $A$ to this circle to be a
combination of the waves

$$\cos(2\theta-n\delta)\quad\text{and}\quad\cos(6\theta-n\delta),
\qquad n\in\mathbb Z,$$

plus a mean term tied to the value of $A$ at the poles.  The integer
$n$ — how fast the wave precesses with the root separation $\delta$ —
is the direction along which "SDP degree" grows: a degree-$d$ block can
reach modulations $|n|\lesssim d$ only.  The generators below are
canonical polynomial representatives of these waves, corrected so that
the pole-tie also holds.

## 3. The mode-2 family $C_n$

**Definition.**

$$C_0=T_2(t_1)=2t_1^2-1,\qquad C_1=2t_1t_2-s,\qquad
\boxed{\;C_{n+1}=2s\,C_n-C_{n-1}\;}$$

(run backwards for $n<0$).  On the coplanar circle $C_n$ restricts to
$\cos(2\theta-n\delta)$: the seed $C_1$ comes from
$\cos(2\theta-\delta)=\cos\big(\theta+(\theta-\delta)\big)
=t_1t_2-\sin\theta\sin(\theta-\delta)$ together with
$\sin\theta\sin(\theta-\delta)=s-t_1t_2$, and the recurrence is
$2\cos\delta\cos(2\theta-n\delta)=\cos(2\theta-(n{+}1)\delta)
+\cos(2\theta-(n{-}1)\delta)$, kept polynomial because $\cos\delta=s$
is a variable.  Explicitly (check `explicit expansions of C_2, C_3,
C_4, C_-1`):

$$\begin{aligned}
C_{-1}&=4st_1^2-2t_1t_2-s,\\
C_2&=4st_1t_2-2s^2-2t_1^2+1,\\
C_3&=8s^2t_1t_2-4s^3-4st_1^2-2t_1t_2+3s,\\
C_4&=16s^3t_1t_2-8s^4-8s^2t_1^2-8st_1t_2+8s^2+2t_1^2-1 .
\end{aligned}$$

**Lift ambiguity.**  A polynomial representing a given function on the
coplanar variety is unique only modulo the coplanarity relation
$D=1+2t_1t_2s-t_1^2-t_2^2-s^2=\det\operatorname{Gram}(x_1,x_2,y)$,
which vanishes there.  For example $\cos(2\theta-2\delta)=T_2(t_2)$ on
the circle, and indeed (check `C_2 = T_2(t_2) + 2 D`)

$$C_2=T_2(t_2)+2D .$$

The recurrence merely fixes one convenient representative per wave; the
$D$-multiples are accounted separately (the "$(M,Q)$ kernel family" of
[E1 admissible](E1_ADMISSIBLE.md) §3).

**Closed-form solution of the recurrence** (check `closed form C_n =
...`): in the Chebyshev basis of the $s$-recurrence,

$$C_n=C_0\,T_n(s)+2t_1(t_2-st_1)\,U_{n-1}(s),\qquad
S_n=S_0\,T_n(s)+(S_1-sS_0)\,U_{n-1}(s),$$

whence on the Gram body $|\hat C_n|\le\tfrac43+4n$ and
$|\hat S_n|\le\tfrac43+12n$: the theta atoms
$\Theta_q=\sum_nq^{n^2}\hat C_n$ converge for every $|q|<1$
([Wrapper lemmas](WRAPPER_LEMMAS.md) L7).

**One-step lifts in general.**  For every $m\ge2$, exactly on the
coplanar circle (check `one-step lifts ... for m = 2..7`):

$$\cos(m\theta-\delta)=U_{m-1}(t_1)\,t_2-U_{m-2}(t_1)\,s,
\qquad
\cos(m\theta+\delta)=U_m(t_1)\,s-U_{m-1}(t_1)\,t_2,$$

with $U_k$ the Chebyshev polynomials of the second kind
($U_k(\cos\theta)=\sin((k{+}1)\theta)/\sin\theta$).  For $m=2$ these
reproduce $C_1$ and $C_{-1}$.

## 4. The pole-value lemma and the $\tfrac13T_n(s)$ correction

The raw $C_n$ is **not** yet (E1)-admissible.  When the leaf $y$ sits
at a pole, $t_1=t_2=0$, and (check `pole-value lemma`)

$$C_n(0,0,s)\;=\;-\,T_n(s)\qquad(\text{likewise } S_n(0,0,s)=-T_n(s)).$$

The full (E1) condition for the equator-pair configuration is

$$\tfrac13\,\big[\text{value at poles}\big]
+\tfrac23\,\big[\text{mode-0 of the circle restriction}\big]=0 ,$$

and the circle restriction $\cos(2\theta-n\delta)$ has zero mean.
Adding a function $g(s)$ of the root separation alone (mean $g(\cos\delta)$,
pole value $g(s)$) and solving
$\tfrac13\!\left(-T_n+g\right)+\tfrac23\,g=0$ gives $g=\tfrac13T_n$.
Hence the admissible generators:

$$\boxed{\;\hat C_n=C_n+\tfrac13\,T_n(s),\qquad
\hat S_n=S_n+\tfrac13\,T_n(s),\qquad n\in\mathbb Z\;}$$

(checks `modulated families ... admissible for every n` and
`C^_-1, C^_-2, S^_-1 are admissible`; the remaining configuration
classes — one root at a pole, both roots at poles — are verified to
hold as well).  Explicitly:

$$\begin{aligned}
\hat C_0&=2t_1^2-\tfrac23\;=\;(T_2+\tfrac13)(t_1),\\
\hat C_1&=2t_1t_2-\tfrac23\,s,\\
\hat C_2&=4st_1t_2-\tfrac43\,s^2-2t_1^2+\tfrac23,\\
\hat C_3&=8s^2t_1t_2-\tfrac83\,s^3-4st_1^2-2t_1t_2+2s,\\
\hat C_4&=16s^3t_1t_2-\tfrac{16}3\,s^4-8s^2t_1^2-8st_1t_2
+\tfrac{16}3\,s^2+2t_1^2-\tfrac23 .
\end{aligned}$$

$\hat C_0$ is precisely the admissible one-root leaf $T_2+\tfrac13$:
the seven-dimensional one-root layer of
[E1 admissible](E1_ADMISSIBLE.md) §1 is the $n=0$ slice of these
families.

## 5. The mode-6 family $S_n$

Same recurrence, seeded at the surviving circle mode 6 (check
`S_1 = U_5(t1) t2 - U_4(t1) s and S_-1 = ...`):

$$S_0=T_6(t_1),\qquad
S_1=U_5(t_1)\,t_2-U_4(t_1)\,s,\qquad
S_{n+1}=2s\,S_n-S_{n-1},$$

with $U_5(t)=32t^5-32t^3+6t$, $U_4(t)=16t^4-12t^2+1$, and backwards
$S_{-1}=U_6(t_1)\,s-U_5(t_1)\,t_2$, $U_6(t)=64t^6-80t^4+24t^2-1$.
Corrected members:

$$\begin{aligned}
\hat S_0&=T_6(t_1)+\tfrac13,\\
\hat S_1&=U_5(t_1)\,t_2-U_4(t_1)\,s+\tfrac13\,s,\\
\hat S_2&=2s\,U_5(t_1)\,t_2-2s^2\,U_4(t_1)-T_6(t_1)+\tfrac13(2s^2-1).
\end{aligned}$$

## 6. A worked verification: $\hat C_1$

$\hat C_1=2t_1t_2-\tfrac23 s$ against every configuration class of
$\mu^*=\frac13\delta_{\pm e}+\frac23\nu$:

* **Both roots on the equator** (separation $\delta$).  On the circle:
  $2\cos\theta\cos(\theta-\delta)-\tfrac23\cos\delta
  =\cos(2\theta-\delta)+\tfrac13\cos\delta$; the wave integrates to
  zero against every admissible $\nu$ (it is a mode-2 wave and
  $\hat\nu(2)=0$), leaving $\tfrac13\cos\delta$.  At the poles
  $t_1=t_2=0$, so $\hat C_1=-\tfrac23\cos\delta$.  Total:
  $\tfrac13\big({-\tfrac23\cos\delta}\big)
  +\tfrac23\big(\tfrac13\cos\delta\big)=0$ for every $\delta$ and every
  admissible $\nu$.
* **One root at a pole** ($s=0$): $\hat C_1=2t_1t_2$; on the equator
  $t_1=0$, at the poles $t_2=0$ — every term vanishes pointwise.
* **Both roots at poles** ($x_2=\pm x_1$, $s=\pm1$, $t_2=\pm t_1$):
  $\hat C_1=\pm(2t_1^2-\tfrac23)$; equator value $\mp\tfrac23$, pole
  value $\pm\tfrac43$; total
  $\tfrac13(\pm\tfrac43)+\tfrac23(\mp\tfrac23)=0$.

## 7. Degrees, and why these are "the pattern in the degree"

$$\deg\hat C_n=|n|+1\ (n\neq0),\qquad \deg\hat S_n=|n|+5\ (n\neq0),$$

(up to the small offsets of the seeds), so a degree-$d$ relaxation can
only use modulations $|n|\lesssim d$.  This is exactly what raising the
SDP degree buys in the (E1)-admissible cone: nothing in the one-root
layer, new modulation orders here.  The limiting certificate is a
series

$$\sum_{n\in\mathbb Z}\gamma^{(2)}_n\,\hat C_n
+\sum_{n\in\mathbb Z}\gamma^{(6)}_n\,\hat S_n
\qquad\text{inside the Gram blocks},$$

and the observed $\exp(-c\,d^2)$ convergence of the hierarchy bounds
predicts $\gamma_n\sim q^{\,n^2}$.  The proposed finite objects to
adjoin to the SDP are therefore the resummed **theta atoms**

$$\Theta^{(2)}_q=\sum_{n\in\mathbb Z}q^{n^2}\hat C_n,\qquad
\Theta^{(6)}_q=\sum_{n\in\mathbb Z}q^{n^2}\hat S_n,
\qquad q\in\mathbb Q\cap(0,1),$$

whose truncation tails are again sums of admissible squares, so the
one-sided rational cuts of
[Exact zero program](EXACT_ZERO_PROGRAM.md) §4.2 apply.

## 8. Completeness and pointers

These two families do not exhaust the admissible two-root space by
themselves: the full description
([E1 admissible](E1_ADMISSIBLE.md) §3, checks `even sector ...` and
`odd sector ...`) is

* even sector $=$ modulated profiles on circle modes 2 and 6
  (realized by $\hat C_n$, $\hat S_n$)
  $\;\oplus\;$ the kernel family
  $D\,Q-\tfrac13(1-s^2)\,Q(0,0,s)$ with slice conditions on $Q$
  (dimension identity $52=26+26$ at degree 8);
* orientation-odd sector: only the two pole-root slice conditions
  ($B(0,\cdot,0),B(\cdot,0,0)\in\operatorname{span}\{U_1,U_5\}$) —
  nearly free.

`solve_e1.py --export-projection` writes the exact admissible bases per
solver sector; `sos_search.py --e1-project` builds the hierarchy on
them.
