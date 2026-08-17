# Exact structure of the $P_2$ kernel problem

Every identity in this document is verified in exact arithmetic by

```sh
python3 verify_exact_structure.py
```

## 1. Three equivalent forms of the kernel

With $t=\cos\theta$ the angle between two lines,

$$K(t)=32t^6-48t^4+20t^2-\frac43
      \;=\;\underbrace{\cos 6\theta+\cos 2\theta+\frac23}_{\text{Chebyshev}}
      \;=\;\underbrace{2\big(T_3(t)^2+t^2\big)-\frac43}_{\text{squared}} .$$

The Legendre (Gegenbauer) expansion is

$$K=\frac{32}{105}P_0+\frac87P_2-\frac{384}{385}P_4+\frac{512}{231}P_6,$$

so **exactly one coefficient is negative**, the one of degree $4$. Writing
$A_\ell=\mathbb E[P_\ell(X\cdot Y)]\ge0$ (Schoenberg), copositivity is

$$\frac{32}{105}+\frac87A_2-\frac{384}{385}A_4+\frac{512}{231}A_6\ \ge\ 0,$$

which under isotropy ($A_2=0$) is exactly $11-36A_4+80A_6\ge0$.

### The squared form

$$W(t):=T_3(t)^2+t^2=\cos^2 3\theta+\cos^2\theta=16t^6-24t^4+10t^2,\qquad
E(\mu)=2\iint W-\frac43 .$$

So the conjecture is $\iint W\,d\mu\,d\mu\ge\frac23$, equivalently

$$\iint T_3(x\cdot y)^2\,d\mu\,d\mu\;+\;\Big\|\int xx^{\mathsf T}d\mu\Big\|_F^2\ \ge\ \frac23 .$$

The second summand is the frame potential, $\ge\frac13$ with equality iff $\mu$ is
isotropic. The first is $\ge\frac13$ under isotropy — but that statement is
algebraically identical to the target, so the two halves cannot be separated.
Without isotropy $\min\iint T_3^2=0.2678884\ldots<\frac13$, so the coupling is
essential.

In the double angle $u=\cos2\theta=2t^2-1$,

$$W=1-u+2u^3=(1+u)\big(u^2+(1-u)^2\big).$$

## 2. Strata that are settled exactly

### Measures on a great circle

If $\nu$ lives on one great circle, the angle is the planar angle and

$$E(\nu)=|\hat\nu(2)|^2+|\hat\nu(6)|^2+\frac23\ \ge\ \frac23\ >\ 0 .$$

The planar Fourier expansion of $K$ has **no** negative coefficient; all of the
difficulty of the problem is the passage from $S^1$ to $S^2$.

### Pole–equator measures

Let $\mu=w\,\delta_{\pm e}+(1-w)\nu$ with $\nu$ supported on the great circle
$e^{\perp}$. Then

$$\boxed{\;E(\mu)=6\Big(w-\frac13\Big)^{2}
        +(1-w)^2\Big(|\hat\nu(2)|^2+|\hat\nu(6)|^2\Big)\;}$$

an exact sum-of-squares identity. It vanishes precisely when

$$w=\frac13,\qquad \hat\nu(2)=\hat\nu(6)=0 .$$

This is an **infinite-dimensional** (codimension 4) family of zero-energy
measures, closed under $SO(3)$. It contains the orthonormal-basis measure
($\nu$ = two orthogonal lines), the pole + Haar-equator measure, and pole +
regular $m$-gon for every $m\notin\{1,3\}$. These are the "pole–equator
measures" of Bilyk–Matzke–Nathe, arXiv:2409.16508 §8.

Numerically, minimising $E$ over $3$–$30$ atoms from random starts never goes
below $0$, and every one of 100+ converged global minimizers has this form.

### The ONB KKT certificate

For $\mu_{\mathrm{ONB}}$ (mass $\frac13$ on each of three orthogonal lines),
$E=0$ and the rooted potential is exactly

$$U(z)=\int K(z\cdot y)\,d\mu_{\mathrm{ONB}}(y)=32\,z_1^2z_2^2z_3^2\ \ge\ 0 ,$$

which is the first-order KKT inequality, with equality on the support.

## 3. What cannot work

### No two-point LP certificate, at any degree

Suppose $h=\sum_n\hat h_nP_n$ with $\hat h_n\ge0$ for $n\ge1$ and $h\le K$ on
$[-1,1]$. Take $\mu^*=\frac13\delta_{\pm e}+\frac23\sigma_{e^\perp}$, which has
$E(\mu^*)=0$ exactly. Then

$$0=E(\mu^*)\ \ge\ I_h(\mu^*)\ \ge\ \hat h_0\ \ge\ 0,$$

so all inequalities are equalities; $h\le K$ with $I_h(\mu^*)=I_K(\mu^*)$ forces
$h=K$ on the support of the distance distribution of $\mu^*$, which is **all of
$[-1,1]$** (the equator–equator part alone realises every value). Hence
$h\equiv K$, contradicting $\hat h_4=-\frac{384}{385}<0$. $\square$

Imposing isotropy does not rescue this: the same argument applies verbatim with
$\hat h_2$ unrestricted. Any certificate must therefore be genuinely
multi-point.

### The cylindrical decomposition is not enough

Writing $x=(\cos\psi_x)e+(\sin\psi_x)\omega_x$, $a=\cos\psi_x$, $b=\cos\psi_y$,
$p^2=(1-a^2)(1-b^2)$:

$$K(x\cdot y)=\sum_{k=0}^{6}c_k(a,b)\cos\big(k(\varphi_x-\varphi_y)\big),$$

$$\begin{aligned}
c_6&=p^6, & c_5&=12abp^5,\\
c_4&=6p^4\big(11a^2b^2-a^2-b^2\big), &
c_3&=4abp^3\big(55a^2b^2-15a^2-15b^2+3\big).
\end{aligned}$$

If every $c_k$ were a positive semidefinite kernel in $(a,b)$ the conjecture
would follow immediately, since
$\iint c_k\cos(k\Delta\varphi)=\sum_j\lambda_j\big|\int g_j(a)e^{ik\varphi}d\mu\big|^2$.
Only $c_5$ and $c_6$ are (they are rank one). The inner forms of
$c_4,c_3,c_2,c_0$ in the separable features $1,a^2,a^4,\ldots$ have determinants
$-1$, $-60$, $-6480$, $-34992$ respectively. This had to fail: PSD-ness of all
$c_k$ is equivalent to $K$ being a positive-definite kernel on $S^2$, which
$\hat K_4<0$ rules out.


## 4. The SO(3) group reformulation

Lines in $\mathbb R^3$ are exactly the involutions of $SO(3)$: identify $x$ with
$\rho_x$, the $\pi$-rotation about $x$. The product $\rho_x\rho_y$ is a rotation
by $2\theta$ about $x\times y$, so with $\chi_\ell$ the character of the spin-$\ell$
representation ($\chi_\ell(2\theta)=1+2\sum_{m=1}^{\ell}\cos 2m\theta$),

$$\boxed{\;K(\cos\theta)=\tfrac12\big(\chi_3-\chi_2+\chi_1\big)(2\theta)+\tfrac16\;}$$

Let $\omega$ be the push-forward of $\mu$ under $x\mapsto\rho_x$ and set
$A_\ell=\int\pi_\ell(\rho_x)\,d\mu(x)$. Each $\pi_\ell(\rho_x)$ is symmetric
orthogonal and involutive, with $\operatorname{tr}\pi_\ell(\rho_x)=\chi_\ell(\pi)=(-1)^\ell$.
Since $\iint\chi_\ell(\rho_x\rho_y)=\|A_\ell\|_F^2$,

$$E(\mu)=\tfrac12\Big(\|A_3\|_F^2-\|A_2\|_F^2+\|A_1\|_F^2\Big)+\tfrac16 ,$$

so the conjecture is exactly

$$\|A_3\|_F^2+\|A_1\|_F^2+\tfrac13\ \ge\ \|A_2\|_F^2 ,$$

for $A_1,A_2,A_3$ symmetric of sizes $3,5,7$ with traces $-1,1,-1$, realised
simultaneously as averages of $\pi_\ell$ on the involution class. Here
$A_1=2M-I$, so $\|A_1\|_F^2=4p_2-1$. At the ONB the three involutions form the
Klein four-group and $\|A_\ell\|_F^2=\tfrac13,\tfrac73,\tfrac53$, giving $0$.

**Why the zero-energy family looks small here.** Let
$\tau=\text{law of }\rho_x\rho_y$, a probability measure on $SO(3)$ with
$\hat\tau(\ell)=A_\ell^2\succeq0$. For *every* member of the pole–equator family,
$\operatorname{supp}\tau$ lies in

$$O(2)_e=\{\text{rotations about } e\}\ \cup\ \{\pi\text{-rotations about axes}\perp e\},$$

a **one-dimensional subgroup** of the three-dimensional $SO(3)$ (for the ONB it is
the four-element Klein group). So in the group picture the equality set is thin,
unlike on $\mathbb{RP}^2$ where its distance set fills $[-1,1]$.

**But the resulting non-commutative LP collapses.** Seek $h$ on $SO(3)$ with
$h\le f:=\tfrac12(\chi_3-\chi_2+\chi_1)+\tfrac16$ pointwise, matrix Fourier
coefficients $H_\ell\succeq0$ for $\ell\ge1$, and $h_0=0$; then $E\ge0$. If such
an $h$ exists, its conjugation-average $\bar h$ also works ($f$ is a class
function, so $\bar h\le f$; averaging $H_\ell$ gives
$(\operatorname{tr}H_\ell/(2\ell+1))I\succeq0$; $\bar h_0=h_0$). So one may assume
$h$ is a class function, i.e. $h=\sum b_\ell\chi_\ell$ with $b_\ell\ge0$. Each
$\chi_\ell(2\theta)=\langle\pi_\ell(\rho_x),\pi_\ell(\rho_y)\rangle_F$ is a positive
definite kernel on $\mathbb{RP}^2$ with nonnegative Legendre coefficients, so this
cone is *contained* in the two-point cone of §3 — which is already ruled out.
The group formulation therefore does not by itself evade the barrier; its value
is that the non-commutative moment constraints on $(A_1,A_2,A_3)$ — symmetry,
fixed traces, and joint realisability on the involution class — are strictly
stronger than positive definiteness of $\tau$, and are not used by any bound in
the literature.

## 5. Consequences for the search

The minimizer set being infinite-dimensional means any exact certificate must
vanish identically on it. Facial reduction against it is therefore mandatory —
and, as it turns out, **already complete** in the published audit: the
pole–equator faces of regular order $m$ annihilate the target for every
$m\ne1,3$ (verified exactly), yet adding orders $4,5,7$ to the ONB + continuous
reduction leaves the reduced dimension unchanged at $480$.
