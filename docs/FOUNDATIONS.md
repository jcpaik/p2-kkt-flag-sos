# Foundations: problem, certificate framework, and code map

This is the base reference for the $P_2$-kernel copositivity program. Part I states the problem and its equivalent formulations, the exactly solved sub-cases, and the structural barriers.  Part II develops the certificate framework: labels, flags, conditional-expectation squares, the exact certificate theorem, and the wrapper lemmas that turn a verified identity into a proof of the conjecture.  Part III maps the mathematics onto `sos_search.py`, including the complete family inventory of the degree-14 reference command.  The appendix collects the operational conventions (SDPA-GMP build, parameter files, export pipeline, objective shifts).

Companion documents: [SHARP_STRUCTURE.md](SHARP_STRUCTURE.md) (the $E_1$ generator, two-root closed forms, limit/weighted-$E_1$/transport analysis), [ENRICHMENTS.md](ENRICHMENTS.md) (weighted targets, cuts, theta atoms, Jensen/fiber-Toeplitz modules), and [SUBCASES_AND_RECORD.md](SUBCASES_AND_RECORD.md) (cylindrical theorems and the measurement record).

Every identity in Part I marked "verified exactly" is checked in exact arithmetic by `python3 verify_exact_structure.py`.

---

# Part I — The problem

## 1. Lines, energy, copositivity

A line through the origin in $\mathbb R^3$ is a unit vector $x\in S^2$ with $x$ and $-x$ identified; the space of such lines is the real projective plane $\mathbb{RP}^2=S^2/\{x\sim-x\}$.  Because the kernel studied here is even, $K(t)=K(-t)$, its value on two projective lines is well defined by $K(x\cdot y)$.  A probability measure $\mu$ on $\mathbb{RP}^2$ is equivalently an antipodally symmetric probability measure on $S^2$.

For a continuous even kernel $K:[-1,1]\to\mathbb R$, define the energy

$$E(\mu)=\iint K(x\cdot y)\,d\mu(x)\,d\mu(y).$$

The kernel is **copositive** on $\mathbb{RP}^2$ if $E(\mu)\ge0$ for every finite nonnegative measure $\mu$; by homogeneity it suffices to check probability measures.  Copositivity is weaker than positive definiteness: a positive-definite zonal kernel has nonnegative spherical-harmonic coefficients, and the present kernel has a negative degree-four coefficient, so the classical two-point positive-definite-kernel argument cannot apply.

The kernel is

$$K(t)=32t^6-48t^4+20t^2-\frac43, \tag{1.1}$$

arising from the trigonometric expression $\cos(2\theta)+\cos(6\theta)+\frac23$ with $|t|=\cos\theta$.

For an orthonormal basis $e_1,e_2,e_3$, the measure

$$\mu_{\mathrm{ONB}}=\frac13(\delta_{e_1}+\delta_{e_2}+\delta_{e_3})$$

has zero energy.  Any proof of copositivity must therefore be sharp.

## 2. Equivalent forms of the kernel

### 2.1 Three algebraic forms

With $t=\cos\theta$ the angle between two lines,

$$K(t)=32t^6-48t^4+20t^2-\frac43 \;=\;\underbrace{\cos6\theta+\cos2\theta+\frac23}_{\text{Chebyshev}} \;=\;\underbrace{2\big(T_3(t)^2+t^2\big)-\frac43}_{\text{squared}}.$$

### 2.2 Legendre data

The Legendre (Gegenbauer) expansion is

$$K=\frac{32}{105}P_0+\frac87P_2-\frac{384}{385}P_4+\frac{512}{231}P_6,$$

so **exactly one coefficient is negative**, the one of degree $4$. Writing $A_\ell=\mathbb E[P_\ell(X\cdot Y)]\ge0$ (Schoenberg), copositivity is

$$\frac{32}{105}+\frac87A_2-\frac{384}{385}A_4+\frac{512}{231}A_6\ \ge\ 0,$$

which under isotropy ($A_2=0$) is exactly $11-36A_4+80A_6\ \ge\ 0$.

### 2.3 The squared form and the frame potential

$$W(t):=T_3(t)^2+t^2=\cos^23\theta+\cos^2\theta=16t^6-24t^4+10t^2, \qquad E(\mu)=2\iint W-\frac43.$$

So the conjecture is $\iint W\,d\mu\,d\mu\ge\frac23$, equivalently

$$\iint T_3(x\cdot y)^2\,d\mu\,d\mu\;+\;\Big\|\int xx^{\mathsf T}d\mu\Big\|_F^2\ \ge\ \frac23.$$

The second summand is the frame potential, $\ge\frac13$ with equality iff $\mu$ is isotropic.  The first is $\ge\frac13$ under isotropy — but that statement is algebraically identical to the target, so the two halves cannot be separated.  Without isotropy $\min\iint T_3^2=0.2678884\ldots<\frac13$, so the coupling is essential.

### 2.4 Double-angle form

In $u=\cos2\theta=2t^2-1$,

$$W=1-u+2u^3=(1+u)\big(u^2+(1-u)^2\big).$$

## 3. Pair moments, isotropy, and $h_2$

Write $p_j(\mu)=\iint(x\cdot y)^j\,d\mu\,d\mu$.  Then (1.1) gives the target

$$E(\mu)=-\frac43+20p_2-48p_4+32p_6. \tag{3.1}$$

The second-moment matrix is $M_2(\mu)=\int xx^{\mathsf T}d\mu(x)$, and $\mu$ is **isotropic** if $M_2(\mu)=\frac13I$.  Since $p_2=\operatorname{tr}(M_2^2)$ and $\operatorname{tr}M_2=1$, always $p_2\ge\frac13$, with equality precisely at isotropy.  The deficit is a harmonic square:

$$h_2=\frac{3p_2-1}{2}=\sum_m|\hat\mu_{2m}|^2 =\frac32\Big\|M_2-\tfrac13I\Big\|_F^2\ \ge\ 0. \tag{3.2}$$

The current code makes **no isotropy assumption**: it uses target (3.1), keeps $p_2$ and every antipodally even moment as an independent label, and encodes the isotropy deficit through the flag square $h_2\ge0$ and the spin-2 Gram blocks (Part III §17).  In the isotropic branch the target reduces to

$$E=\frac{16}{3}-48p_4+32p_6. \tag{3.3}$$

## 4. Exactly solved strata and the minimizer family

### 4.1 Measures on a great circle

If $\nu$ lives on one great circle, the angle is the planar angle and

$$E(\nu)=|\hat\nu(2)|^2+|\hat\nu(6)|^2+\frac23\ \ge\ \frac23\ >\ 0.$$

The planar Fourier expansion of $K$ has **no** negative coefficient; all of the difficulty of the problem is the passage from $S^1$ to $S^2$.

### 4.2 Pole–equator measures

Let $\mu=w\,\delta_{\pm e}+(1-w)\nu$ with $\nu$ supported on the great circle $e^{\perp}$.  Then (verified exactly)

$$\boxed{\;E(\mu)=6\Big(w-\frac13\Big)^{2} +(1-w)^2\Big(|\hat\nu(2)|^2+|\hat\nu(6)|^2\Big)\;}$$

an exact sum-of-squares identity.  It vanishes precisely when

$$w=\frac13,\qquad \hat\nu(2)=\hat\nu(6)=0.$$

This is an **infinite-dimensional** (codimension 4) family of zero-energy measures, closed under $SO(3)$.  It contains the orthonormal-basis measure ($\nu$ = two orthogonal lines), the pole + Haar-equator measure, and pole + regular $m$-gon for every $m\notin\{1,3\}$.  These are the "pole–equator measures" of Bilyk–Matzke–Nathe, arXiv:2409.16508 §8.

Numerically, minimising $E$ over $3$–$30$ atoms from random starts never goes below $0$, and every one of 100+ converged global minimizers has this form.

### 4.3 The ONB KKT certificate

For $\mu_{\mathrm{ONB}}$ (mass $\frac13$ on each of three orthogonal lines), $E=0$ and the rooted potential is exactly

$$U(z)=\int K(z\cdot y)\,d\mu_{\mathrm{ONB}}(y)=32\,z_1^2z_2^2z_3^2\ \ge\ 0,$$

which is the first-order KKT inequality, with equality on the support.

## 5. Barriers: what cannot work

### 5.1 No two-point LP certificate, at any degree

**Proposition.** There is no $h=\sum_n\hat h_nP_n$ with $\hat h_n\ge0$ for $n\ge1$ and $h\le K$ on $[-1,1]$ certifying the conjecture.

**Proof.** Take $\mu^*=\frac13\delta_{\pm e}+\frac23\sigma_{e^\perp}$, which has $E(\mu^*)=0$ exactly.  Then

$$0=E(\mu^*)\ \ge\ I_h(\mu^*)\ \ge\ \hat h_0\ \ge\ 0,$$

so all inequalities are equalities; $h\le K$ with $I_h(\mu^*)=I_K(\mu^*)$ forces $h=K$ on the support of the distance distribution of $\mu^*$, which is **all of $[-1,1]$** (the equator–equator part alone realises every value).  Hence $h\equiv K$, contradicting $\hat h_4=-\frac{384}{385}<0$. $\square$

Imposing isotropy does not rescue this: the same argument applies verbatim with $\hat h_2$ unrestricted.  Any certificate must therefore be genuinely multi-point: the analogue of a magic function must be a rooted, matrix-valued, higher-point positive kernel rather than a single scalar function of one inner product.  (Equivalently: equality in a pointwise minorant proof against the pole–equator zero measure would force $h(t)=K(t)$ on the entire interval, impossible with $\hat K_4<0$.)

### 5.2 The cylindrical decomposition is not enough

Writing $x=(\cos\psi_x)e+(\sin\psi_x)\omega_x$, $a=\cos\psi_x$, $b=\cos\psi_y$, $p^2=(1-a^2)(1-b^2)$:

$$K(x\cdot y)=\sum_{k=0}^{6}c_k(a,b)\cos\big(k(\varphi_x-\varphi_y)\big),$$

$$\begin{aligned} c_6&=p^6, & c_5&=12abp^5,\\ c_4&=6p^4\big(11a^2b^2-a^2-b^2\big), & c_3&=4abp^3\big(55a^2b^2-15a^2-15b^2+3\big). \end{aligned}$$

If every $c_k$ were a positive semidefinite kernel in $(a,b)$ the conjecture would follow immediately, since $\iint c_k\cos(k\Delta\varphi)=\sum_j\lambda_j\big|\int g_j(a)e^{ik\varphi}d\mu\big|^2$. Only $c_5$ and $c_6$ are (they are rank one).  The inner forms of $c_4,c_3,c_2,c_0$ in the separable features $1,a^2,a^4,\ldots$ have determinants $-1$, $-60$, $-6480$, $-34992$ respectively.  This had to fail: PSD-ness of all $c_k$ is equivalent to $K$ being a positive-definite kernel on $S^2$, which $\hat K_4<0$ rules out. (For what the cylindrical decomposition *does* prove, see [SUBCASES_AND_RECORD.md](SUBCASES_AND_RECORD.md).)

## 6. The $SO(3)$ group reformulation

Lines in $\mathbb R^3$ are exactly the involutions of $SO(3)$: identify $x$ with $\rho_x$, the $\pi$-rotation about $x$.  The product $\rho_x\rho_y$ is a rotation by $2\theta$ about $x\times y$, so with $\chi_\ell$ the character of the spin-$\ell$ representation ($\chi_\ell(2\theta)=1+2\sum_{m=1}^{\ell}\cos2m\theta$),

$$\boxed{\;K(\cos\theta)=\tfrac12\big(\chi_3-\chi_2+\chi_1\big)(2\theta)+\tfrac16\;}$$

Let $\omega$ be the push-forward of $\mu$ under $x\mapsto\rho_x$ and set $A_\ell=\int\pi_\ell(\rho_x)\,d\mu(x)$.  Each $\pi_\ell(\rho_x)$ is symmetric orthogonal and involutive, with $\operatorname{tr}\pi_\ell(\rho_x)=\chi_\ell(\pi)=(-1)^\ell$.  Since $\iint\chi_\ell(\rho_x\rho_y)=\|A_\ell\|_F^2$,

$$E(\mu)=\tfrac12\Big(\|A_3\|_F^2-\|A_2\|_F^2+\|A_1\|_F^2\Big)+\tfrac16,$$

so the conjecture is exactly

$$\|A_3\|_F^2+\|A_1\|_F^2+\tfrac13\ \ge\ \|A_2\|_F^2,$$

for $A_1,A_2,A_3$ symmetric of sizes $3,5,7$ with traces $-1,1,-1$, realised simultaneously as averages of $\pi_\ell$ on the involution class.  Here $A_1=2M-I$, so $\|A_1\|_F^2=4p_2-1$.  At the ONB the three involutions form the Klein four-group and $\|A_\ell\|_F^2=\tfrac13,\tfrac73,\tfrac53$, giving $0$.

**Why the zero-energy family looks small here.**  Let $\tau=\text{law of }\rho_x\rho_y$, a probability measure on $SO(3)$ with $\hat\tau(\ell)=A_\ell^2\succeq0$.  For *every* member of the pole–equator family, $\operatorname{supp}\tau$ lies in

$$O(2)_e=\{\text{rotations about }e\}\ \cup\ \{\pi\text{-rotations about axes}\perp e\},$$

a **one-dimensional subgroup** of the three-dimensional $SO(3)$ (for the ONB it is the four-element Klein group).  So in the group picture the equality set is thin, unlike on $\mathbb{RP}^2$ where its distance set fills $[-1,1]$.

**But the resulting non-commutative LP collapses.**  Seek $h$ on $SO(3)$ with $h\le f:=\tfrac12(\chi_3-\chi_2+\chi_1)+\tfrac16$ pointwise, matrix Fourier coefficients $H_\ell\succeq0$ for $\ell\ge1$, and $h_0=0$; then $E\ge0$.  If such an $h$ exists, its conjugation-average $\bar h$ also works ($f$ is a class function, so $\bar h\le f$; averaging $H_\ell$ gives $(\operatorname{tr}H_\ell/(2\ell+1))I\succeq0$; $\bar h_0=h_0$).  So one may assume $h$ is a class function, i.e. $h=\sum b_\ell\chi_\ell$ with $b_\ell\ge0$.  Each $\chi_\ell(2\theta)=\langle\pi_\ell(\rho_x),\pi_\ell(\rho_y)\rangle_F$ is a positive definite kernel on $\mathbb{RP}^2$ with nonnegative Legendre coefficients, so this cone is *contained* in the two-point cone of §5.1 — which is already ruled out.  The group formulation does not by itself evade the barrier; its value is that the non-commutative moment constraints on $(A_1,A_2,A_3)$ — symmetry, fixed traces, and joint realisability on the involution class — are strictly stronger than positive definiteness of $\tau$, and are not used by any bound in the literature.

## 7. Consequences for the search

The minimizer set being infinite-dimensional means any exact certificate must vanish identically on it.  Facial reduction against it is therefore mandatory — and **already complete** in the published audit: the pole–equator faces of regular order $m$ annihilate the target for every $m\ne1,3$ (verified exactly), yet adding orders $4,5,7$ to the ONB + continuous reduction leaves the reduced dimension unchanged at $480$.

---

# Part II — The certificate framework

## 8. KKT conditions at a minimizer

Define the potential generated by $\mu$:

$$U_\mu(z)=\int K(z\cdot y)\,d\mu(y),\qquad E(\mu)=\int U_\mu(x)\,d\mu(x).$$

Suppose $\mu$ is a global minimizer of $E$ over probability measures. For a trial point $z$ and $\mu_\varepsilon=(1-\varepsilon)\mu+\varepsilon\delta_z$,

$$\left.\frac{d}{d\varepsilon}E(\mu_\varepsilon)\right|_{\varepsilon=0} =2\bigl(U_\mu(z)-E(\mu)\bigr).$$

Minimality therefore implies the global first-variation inequality

$$q_\mu(z):=U_\mu(z)-E(\mu)\ \ge\ 0\qquad\text{for every }z. \tag{8.1}$$

Averaging (8.1) against $\mu$ gives zero, hence

$$q_\mu(x)=0\qquad\text{for }\mu\text{-almost every }x. \tag{8.2}$$

Every support point is thus a minimum of the smooth function $U_\mu$ on the sphere, so on the support

$$g_\mu(x):=\nabla_{S^2}U_\mu(x)=0, \tag{8.3}$$

$$H_\mu(x):=\operatorname{Hess}_{S^2}U_\mu(x)\ \succeq\ 0. \tag{8.4}$$

These are Karush–Kuhn–Tucker (**KKT**) conditions: necessary conditions for a constrained optimum, the constraints being positivity and total mass one for the measure together with the spherical constraint on point positions.  For $s=x\cdot y$ and a tangent vector $v\perp x$,

$$\nabla_{S^2}U_\mu(x)=\int K'(s)\,(y-sx)\,d\mu(y), \tag{8.5}$$

$$\operatorname{Hess}_{S^2}U_\mu(x)[v,v] =\int\left[K''(s)(v\cdot y)^2-sK'(s)\lVert v\rVert^2\right]d\mu(y). \tag{8.6}$$

Equations (8.1)–(8.6) are exact analytic facts; they do not depend on an atomic approximation.  The kernel derivatives used throughout are

$$K'(t)=192t^5-192t^3+40t,\qquad K''(t)=960t^4-576t^2+40.$$

The KKT statements are not valid for an arbitrary measure; the proof architecture applies them to a hypothetical negative global minimizer (§12, §13).

## 9. Moment labels as multigraphs

### 9.1 Gram monomials

For independent samples $X_0,\ldots,X_{n-1}\sim\mu$, a Gram monomial is

$$m_{\mathbf e}(X_0,\ldots,X_{n-1}) =\prod_{0\le i<j<n}(X_i\cdot X_j)^{e_{ij}},$$

with nonnegative integer exponents.  Regard $\mathbf e$ as a loopless multigraph: sampled points are vertices, the multiplicity of edge $ij$ is $e_{ij}$, and the total Gram degree is $\sum_{i<j}e_{ij}$.  The corresponding **moment label** is $M_{\mathbf e}=\mathbb E[m_{\mathbf e}]$. Only a canonical label is stored, after the following exact reductions.

### 9.2 Antipodal parity

Changing the sign of one projective representative $X_i$ leaves its law unchanged and changes the monomial by $(-1)^{d_i}$, $d_i=\sum_{j\ne i}e_{ij}$.  Hence $d_i$ odd for some $i$ implies $M_{\mathbf e}=0$: every surviving vertex has even multigraph degree.

### 9.3 Isotropic contraction (isotropic branch only)

If a sampled vertex has multigraph degree two it can be integrated out exactly under isotropy:

$$\mathbb E_{X_i}[(X_i\cdot u)(X_i\cdot v)] =u^{\mathsf T}\mathbb E[X_iX_i^{\mathsf T}]v=\frac13(u\cdot v),$$

so a degree-two path $u-X_i-v$ contracts to an edge $u-v$ with factor $1/3$, and a doubled edge gives $\mathbb E_{X_i}[(X_i\cdot u)^2]=\frac13$. The isotropic-branch reducer repeats these contractions until no degree-two vertex remains.  **The current non-isotropic code applies no such contraction**: `("pair", 2)` and all degree-two-vertex graphs are independent labels, and the isotropy deficit is carried by the `harmonic_2`, `harmonic_flag_*`, and `spin2_flag` blocks.

### 9.4 Isolated vertices, disconnected components, exchangeability

An isolated sampled vertex contributes $1$ and is removed.  If the graph splits into connected components $H_1,\ldots,H_s$, independence gives $M_H=\prod_\alpha M_{H_\alpha}$; the code represents this by a `product` label.  The samples are iid, so relabeling vertices does not change the moment; each connected graph is replaced by the lexicographically least edge-exponent tuple among all vertex permutations.  Surviving labels have forms such as

```text
("constant",)
("pair", r)
("triangle", e01, e02, e12)
("graph_4", ...)
("graph_5", ...)
("product", label_1, label_2, ...)
```

### 9.5 Exact geometric identities

Besides antipodal parity and (branch-dependent) isotropic contraction, the label algebra uses:

1. **Rank-three Gram identities.**  Four vectors in $\mathbb R^3$ are linearly dependent, so $\det\operatorname{Gram}(X_1,X_2,X_3,X_4)=0$; the identity remains zero after multiplication by any Gram polynomial and averaging.
2. **Exchangeability**, as in §9.4.
3. **Sphere relations** $|x|^2=1$.

## 10. Flags and flag squares

### 10.1 Rooted flags

A **rooted flag** keeps some sampled points fixed as roots and averages over the remaining points.  With roots $X_1,X_2,X_3$,

$$\Phi_{\mathbf e}(X_1,X_2,X_3) =\int\prod_{i=1}^3(X_i\cdot Y)^{e_i}\,d\mu(Y).$$

For any real coefficients $c_{\mathbf e}$,

$$\mathbb E\left[\Big(\sum_{\mathbf e}c_{\mathbf e} \Phi_{\mathbf e}(X_1,X_2,X_3)\Big)^2\right]\ \ge\ 0. \tag{10.1}$$

Expanding (10.1) introduces two independent leaves $Y,Z$ sharing the same roots: it becomes a positive semidefinite constraint on moments of five sampled points.  This is the flag-algebra mechanism: square a conditionally averaged rooted expression, then forget the roots by averaging — analogous to Razborov's flag algebras, with Gram monomials replacing finite graph densities.

### 10.2 The flag moment matrix and gluing

Generally, for shared roots $R=(X_0,\ldots,X_{r-1})$, a leaf $Y$, and polynomial flag functions $f_1,\ldots,f_s$ of the Gram entries, define $F(R)=\mathbb E_Y[(f_1(R,Y),\ldots,f_s(R,Y))^{\mathsf T}]$.  For every $c\in\mathbb R^s$, $\mathbb E_R[(c^{\mathsf T}F(R))^2]\ge0$; equivalently the flag moment matrix

$$\mathcal M=\mathbb E_R[F(R)F(R)^{\mathsf T}]\ \succeq\ 0,$$

whose entries, after gluing a second independent leaf $W$, are $\mathcal M_{ij}=\mathbb E_{R,Y,W}[f_i(R,Y)f_j(R,W)]$.  Each entry reduces to canonical moment labels, giving coefficient matrices $A_L$ with $\mathcal M(y)=\sum_Ly_LA_L$ in the dual variables $y_L$.

### 10.3 Degree, arity, parity sectors

The **arity** is the number of sampled vertices appearing after a square is expanded (one root + two leaves = three-point block; two roots + two leaves = four-point; three roots + two leaves = five-point).  The **degree** is the total polynomial degree in Gram entries.  Increasing degree enriches the functions available on a fixed configuration; increasing arity adds genuinely new multipoint correlations.

Flags are separated according to their parity under $X_i\mapsto-X_i$.  Only flags with the same parity character are paired in a Gram block; this ensures the unrooted square is well defined on projective space.

### 10.4 Root-weighted flags

The strongest current blocks also permit factors among the roots:

$$\Phi_{\mathbf r,\mathbf e}(X_1,\ldots,X_s) =\Big(\prod_{i<j}(X_i\cdot X_j)^{r_{ij}}\Big) \int\prod_i(X_i\cdot Y)^{e_i}\,d\mu(Y). \tag{10.2}$$

## 11. KKT-infused positive and free families

Ordinary flag squares use only expressions such as (10.1).  The KKT conditions provide additional nonnegative expressions:

$$\int\rho_\mu(z)\,q_\mu(z)\,d\sigma(z)\ \ge\ 0 \qquad\text{if }\rho_\mu(z)\ge0, \tag{11.1}$$

$$\int\operatorname{tr}\!\big(B_\mu(x)\operatorname{Hess}_{S^2}U_\mu(x)\big)\,d\mu(x)\ \ge\ 0 \qquad\text{if }B_\mu(x)\succeq0. \tag{11.2}$$

The multipliers $\rho_\mu$ and $B_\mu$ are themselves constructed from rooted flag squares.  The equalities (8.2) and (8.3) may be multiplied by arbitrary signed flag expressions, because their averages remain zero.  A KKT-infused flag certificate therefore combines

$$\text{ordinary flag squares} +\text{PSD multipliers of the Hessian} +\text{nonnegative multipliers of the global gap} +\text{free multipliers of KKT equalities},$$

every term justified by a square, a KKT inequality, or an exact identity.  ("KKT-infused flag algebra" is descriptive, not a standard named theorem.)

## 12. The abstract exact certificate and the certificate theorem

An exact KKT-infused flag certificate is an identity

$$\begin{aligned} E(\mu)={}& \sum_\alpha\left\|\Phi_{\alpha,\mu}\right\|_{L^2(\mu^{r_\alpha})}^2\\ &+\int\operatorname{tr}\!\left(B_\mu(x)H_\mu(x)\right)d\mu(x)\\ &+\int_{\mathbb{RP}^2}\rho_\mu(z)\,q_\mu(z)\,d\sigma(z)\\ &+\int a_\mu(x)\,q_\mu(x)\,d\mu(x)\\ &+\int v_\mu(x)\cdot g_\mu(x)\,d\mu(x)\\ &+\mathcal Z(\mu). \end{aligned} \tag{12.1}$$

The datum in (12.1) must satisfy:

- every $\Phi_{\alpha,\mu}$ is a rooted flag evaluation;
- $B_\mu(x)\succeq0$ pointwise on the tangent plane;
- $\rho_\mu(z)\ge0$;
- $a_\mu$ and $v_\mu$ are arbitrary real equality multipliers;
- $\mathcal Z(\mu)=0$ follows exactly from antipodal symmetry, isotropy, exchangeability, sphere equations, and Gram-rank identities;
- expansion of both sides into moment labels agrees exactly.

At finite degree, the first three types of nonnegative multipliers are encoded by rational positive semidefinite Gram matrices.

**Certificate theorem.**  Suppose (12.1) is an exact identity for every probability measure in a class on which the KKT conditions (8.1)–(8.4) hold.  Then every such measure has $E(\mu)\ge0$.

**Proof.**  Each squared norm in (12.1) is nonnegative.  Since $B_\mu(x)\succeq0$ and $H_\mu(x)\succeq0$, their Frobenius contraction is nonnegative: $\operatorname{tr}(B_\mu(x)H_\mu(x))\ge0$.  Likewise $\rho_\mu(z)q_\mu(z)\ge0$ by (8.1).  The two equality-multiplier integrals vanish by (8.2)–(8.3), and $\mathcal Z(\mu)=0$ identically. Therefore the right side of (12.1) is nonnegative, so $E(\mu)\ge0$. $\square$

**Global reduction theorem.**  Assume:

1. every global minimizer with negative energy lies in the class covered by (12.1); and
2. (12.1) has been verified exactly on that class.

Then $K$ is copositive.

**Proof.**  If copositivity failed, some probability measure would have negative energy.  The space of probability measures on the compact space $\mathbb{RP}^2$ is weak-\* compact and, $K$ being continuous, $E$ is weak-\* continuous; hence a global minimizer $\mu_*$ exists and satisfies $E(\mu_*)<0$.  Global minimality gives the KKT conditions (8.1)–(8.4).  Assumption 1 places $\mu_*$ in the class covered by the certificate; the certificate theorem then gives $E(\mu_*)\ge0$, a contradiction. $\square$

**Coverage.**  The current program makes no isotropy assumption: degree-two sampled vertices are retained in the canonical moment labels, so the certificate algebra covers every antipodally symmetric probability measure satisfying the KKT conditions.  Antipodal symmetry itself is exact ($K$ is even, so symmetrizing $\mu\mapsto\tfrac12(\mu+\mu^-)$ preserves $E$; lemma L0 below), and a negative minimizer may be taken antipodal without loss of generality. Consequently a rationally verified certificate in the current algebra would prove unrestricted copositivity directly — assumption 1 holds for the full class.  The isotropy deficit is represented inside the hierarchy by harmonic flag squares: $h_2=\sum_m|\hat\mu_{2m}|^2\ge0$ and spin-2 Gram blocks containing the deviatoric second moment. However (see [SUBCASES_AND_RECORD.md](SUBCASES_AND_RECORD.md)) the numerical dual value of the non-isotropic hierarchy at degree 14 is approximately $-4.5\times10^{-3}$, far from zero, and the convex relaxation of the spin-2 sector provably leaks at order $\sqrt{h_2}$ against the linear cost of $h_2$.  Coverage is no longer the obstacle; existence of a finite-degree certificate is (the weighted/cut program of [ENRICHMENTS.md](ENRICHMENTS.md) addresses this).

## 13. Wrapper lemmas

These statements turn an exact SDP certificate into a proof of the conjecture.  Notation: $\mathcal M$ is the set of Borel probability measures on $S^2$, $\mathcal M_a\subset\mathcal M$ the antipodally symmetric ones.

### L0. Reduction to antipodal measures

$K$ is even, so for the symmetrization $\bar\mu=\tfrac12(\mu+\mu_-)$ (with $\mu_-$ the pushforward under $x\mapsto-x$) all four sign-combinations in the double integral are equal and $E(\bar\mu)=E(\mu)$.  Hence $\inf_{\mathcal M}E=\inf_{\mathcal M_a}E$ and it suffices to prove $E\ge0$ on $\mathcal M_a$.

### L1. Existence of a minimizer

$\mathcal M_a$ is weak-\* compact (Banach–Alaoglu; antipodality and total mass are weak-\* closed constraints) and $E$ is weak-\* continuous ($K$ is a polynomial, so $E$ is a finite sum of products of moments).  Hence $E$ attains its minimum at some $\mu^\dagger\in\mathcal M_a$.

### L2. First-order conditions at a minimizer

For $\nu\in\mathcal M_a$ and $t\in[0,1]$,

$$E\big((1-t)\mu^\dagger+t\nu\big) =E(\mu^\dagger)+2t\Big(\textstyle\int U_{\mu^\dagger}\,d\nu-E(\mu^\dagger)\Big)+O(t^2),$$

so minimality gives $\int U_{\mu^\dagger}\,d\nu\ge E(\mu^\dagger)$ for every $\nu$; taking $\nu=\tfrac12(\delta_z+\delta_{-z})$:

$$U_{\mu^\dagger}(z)\ \ge\ E(\mu^\dagger)\qquad\text{for every }z\in S^2.$$

Moreover $\int(U_{\mu^\dagger}-E)\,d\mu^\dagger=0$ with a continuous nonnegative integrand, so $U_{\mu^\dagger}=E(\mu^\dagger)$ on $\operatorname{supp}\mu^\dagger$.

### L3. Second-order conditions at support points

By L2 every $x_0\in\operatorname{supp}\mu^\dagger$ is a global minimum of the smooth function $U_{\mu^\dagger}$ on $S^2$, hence $\nabla_{S^2}U_{\mu^\dagger}(x_0)=0$ and $\operatorname{Hess}_{S^2}U_{\mu^\dagger}(x_0)\succeq0$.

### L4. Identities valid for every measure

The following hold identically and may be multiplied by *free-sign* multipliers: $|x|^2=1$ on the sphere; vanishing of odd antipodal moments; exchangeability of sample points; and $\det\operatorname{Gram}(x_1,x_2,x_3,x_4)=0$ for any four points of $\mathbb R^3$ (the `--rank-relations` family).  The identity $\int U_\mu\,d\mu=E(\mu)$ is the constant-multiplier case of the potential relation; nonconstant potential/gradient/Hessian multipliers and the global-gap terms are **not** identities — they are valid only through L2/L3 at a minimizer.

### L5. Validity of a certificate

Suppose rational data give the label-algebra identity

$$E-\lambda=\sum_\beta\langle Q_\beta,\,G_\beta(\mu)\rangle +\sum_j c_j\,R_j(\mu)+\sum_k f_k\,I_k(\mu)$$

with $Q_\beta\succeq0$ (so each pairing $\ge0$ at every measure), $c_j\ge0$ multiplying quantities $R_j\ge0$ **at minimizers** (the L2/L3 families), and $f_k$ free multiplying the L4 identities $I_k\equiv0$.  Evaluating at $\mu^\dagger$ gives $E(\mu^\dagger)\ge\lambda$, hence $E\ge\lambda$ on all of $\mathcal M_a$ by L1.  If no $R_j$ terms appear, the same display evaluated at an arbitrary $\mu$ gives $E\ge\lambda$ directly, with no appeal to L1–L3 (an *all-measures* certificate).

### L6. The weighted target and its composition requirement

With $h_2=\tfrac{3p_2-1}2\ge0$ (zero exactly on isotropic measures), the reduction lemma (see [ENRICHMENTS.md](ENRICHMENTS.md)) states: if $h_2E\ge0$ **for every** $\mu\in\mathcal M_a$, then $E\ge0$ on $\mathcal M_a$ (divide on $\{h_2>0\}$; on $\{h_2=0\}$ perturb by $\delta_{\pm e}$ and use weak-\* continuity).

*Composition caveat.*  The hypothesis must hold for all measures.  A certificate of $h_2E\ge0$ that uses the L2/L3 families proves the inequality only at minimizers of $E$; combined with L1 it yields "$h_2(\mu^\dagger)E(\mu^\dagger)\ge0$", and since a hypothetical counterexample minimizer could be isotropic ($h_2=0$) — indeed the conjectured zero family is isotropic — this is not a contradiction. Valid conclusions require either (i) an all-measures weighted certificate (L5 with no $R_j$ terms), or (ii) L2/L3 re-derived for the functional $W=h_2E$ itself.  By the product rule its first-variation potential is

$$\Phi^W_\mu(z)=h_2(\mu)\,U_\mu(z)+E(\mu)\,u_\mu(z)+c(\mu),\qquad u_\mu(z)=\tfrac32\int(z\cdot x)^2\,d\mu(x),$$

with the $\mu$-dependent constant $c(\mu)$ fixed by $\int\Phi^W_\mu\,d\mu=2W(\mu)$; every term is polynomial in the label algebra.  With this potential, L1–L3 and L5 apply to $W$ verbatim and give $W\ge\lambda$ on all of $\mathcal M_a$.  (Not implemented in the solver.)

### L7. Theta atoms: convergence and truncation cuts

Let $\hat C_n=C_n+\tfrac13T_n(s)$, $\hat S_n=S_n+\tfrac13T_n(s)$ be the modulated two-root generators (see [SHARP_STRUCTURE.md](SHARP_STRUCTURE.md)).  The recurrence solves in closed form (machine-checked exactly in `solve_e1.py`):

$$C_n=C_0\,T_n(s)+2t_1(t_2-st_1)\,U_{n-1}(s),\qquad S_n=S_0\,T_n(s)+(S_1-sS_0)\,U_{n-1}(s).$$

On the Gram body ($|t_1|,|t_2|,|s|\le1$) we have $|T_n|\le1$, $|U_{n-1}|\le n$, $|C_0|,|S_0|\le1$, $|2t_1(t_2-st_1)|\le4$, $|S_1-sS_0|\le12$, so

$$|\hat C_n|\le\tfrac43+4n,\qquad |\hat S_n|\le\tfrac43+12n.$$

Hence for rational $q\in(0,1)$ the theta atoms $\Theta_q=\sum_{n\in\mathbb Z}q^{n^2}\hat C_n$ (and the $\hat S$ analogue) converge absolutely and uniformly, and every Gram pairing of $\Theta_q$-leaves is the absolutely convergent sum of the pairings of its terms.  Since each $\hat C_n$-block pairing is PSD at every measure, the tails are PSD and the one-sided cuts

$$\operatorname{Pair}(\Theta_q)\;\succeq\;\sum_{|n|\le N}q^{n^2}\operatorname{Pair}(\hat C_n) \qquad(N=0,1,2,\dots)$$

are valid at every measure, with exactly rational data ($q^{n^2}$ and the $\hat C_n$ coefficients are rational).  The adjoined atom labels never need closed-form values: they enter the SDP only as variables bounded by these cuts.  (The theta/cut program itself lives in [ENRICHMENTS.md](ENRICHMENTS.md).)

## 14. From squares to an SDP

Choose finite bases of flags.  Give each family an unknown symmetric Gram matrix $Q$; requiring $Q\succeq0$ makes the associated quadratic form a sum of squares.  Expand every block in the common basis of moment labels.  With coefficient matrices $A_L^{(\beta)}$ for each positive family $\beta$, $B_L^{(\gamma)}$ for each matrix equality block $\gamma$, and $R_L^{(r)}$ for each scalar identity $r$, an exact primal certificate is a coefficient identity

$$t_L=\sum_\beta\big\langle A_L^{(\beta)},Q_\beta\big\rangle +\sum_\gamma\big\langle B_L^{(\gamma)},S_\gamma\big\rangle +\sum_r\lambda_rR_L^{(r)}\qquad\text{for every label }L,$$

with $Q_\beta\succeq0$, while $S_\gamma$ and $\lambda_r$ are unrestricted because they multiply identities equal to zero.

The conic **dual** is the formal moment problem

$$\begin{aligned} \text{minimize}\quad& \sum_Lt_Ly_L\quad(\text{the target coefficients}),\\ \text{subject to}\quad& y_{\mathrm{constant}}=1,\\ &\textstyle\sum_Ly_LA_L^{(\beta)}\succeq0 &&\text{for every positive block }\beta,\\ &\textstyle\sum_Ly_LB_L^{(\gamma)}=0 &&\text{for every matrix equality block }\gamma,\\ &\textstyle\sum_Ly_LR_L^{(r)}=0 &&\text{for every scalar identity }r. \end{aligned}$$

Every actual measure satisfying the encoded KKT conditions defines a feasible vector $y_L=M_L(\mu)$; the formal feasible set is larger, so the dual optimum is a lower bound for the encoded branch.  Three standing qualifications:

1. formal moment variables need not come from an actual measure — a negative dual optimum is a relaxation gap, not a counterexample;
2. a numerical value close to zero is **not** an exact proof unless a primal decomposition is recovered and checked exactly;
3. a dual optimum rigorously $\ge0$ would show every actual measure in the encoded class has $E\ge0$; an exact rational primal identity gives the machine-checkable certificate.

## 15. Exactness requirements and verification policy

A publishable finite certificate must contain all of the following.

**Exact Gram matrices.**  Every SOS and positive multiplier block must be a matrix over $\mathbb Q$, or an explicitly described algebraic number field, proved positive semidefinite — e.g. by an exact $LDL^T$ decomposition with nonnegative diagonal or an exact factorization $Q=R^{\mathsf T}R$.

**Exact coefficient identity.**  After expanding every block, the coefficient of every canonical moment label must equal the target coefficient exactly.  A small floating-point residual is not sufficient.

**Exact equality identities.**  All free terms must be explicit combinations of identities that truly vanish: odd-degree antipodal moments; isotropic degree-two contractions (isotropic branch); permutation/exchangeability relations; sphere relations; $\det\operatorname{Gram}_4=0$ and its polynomial multiples; the KKT equalities on the support.

**Correctly signed KKT multipliers.**  Multipliers of $q_\mu(z)\ge0$ must be nonnegative, and multipliers contracted with $H_\mu(x)\succeq0$ must be positive semidefinite tangent matrices. Only multipliers of equality constraints may be unrestricted.

**Coverage of the minimizing class.**  The certificate must apply to every possible negative global minimizer (§12, coverage).

**Infinite-degree alternative.**  It may happen that no finite exact identity exists, while for every $\varepsilon>0$ there is a certificate $E(\mu)+\varepsilon=C_\varepsilon(\mu)$ with $C_\varepsilon$ a valid KKT/SOS expression.  Letting $\varepsilon\downarrow0$ still proves $E(\mu)\ge0$, provided each certificate is valid on the full minimizing class.  Such a family is an asymptotic certificate; its datum is a convergent family of positive rooted kernels and matrix-valued Hessian multipliers — the higher-point analogue of a magic function.

**Verification policy.**  A solver status of `optimal` means only that the numerical termination criteria were satisfied; it is not a proof. An exact certificate should be distributed as: (1) rational Gram matrices and free multipliers; (2) a machine-readable basis/label manifest; (3) an independent exact-arithmetic verifier (`verify_certificate.py`, `verify_exact_structure.py`); (4) a short human-readable derivation of every positive and zero term.

---

# Part III — Mathematics-to-code map

## 16. `sos_search.py` overview

`sos_search.py` constructs either side of the finite semidefinite relaxation of §14:

- the **primal** searches for SOS Gram matrices and KKT multipliers whose expansion equals the target $E$;
- the **dual** (`--dual`) searches for a formal moment functional satisfying every encoded positivity/equality constraint while minimizing $E$.  The dual is numerically much more stable in the current singular problem, but does not itself provide certificate matrices.

Moment labels are stored as multigraphs (§9).  `reduce_graph_matrix` applies antipodal parity, factorization over disconnected components, and canonical relabeling under vertex permutations; the canonical cached entry point is `graph_expectation_label`.  No isotropy contraction is applied: degree-two vertices stay in the label, so `("pair", 2)` and its relatives are genuine moment variables. Unreduced connected moments become labels such as `triangle`, `graph_4`, `graph_5`, or `graph_6`.

## 17. Block-family directory

Positive (PSD-constrained) families:

| family | mathematics |
|---|---|
| `flag_*` | one-root $O(2)$ harmonic flag squares (§19.3) |
| `two_root_*` | two-root even and orientation-odd sectors (§19.5) |
| `star_flag_*` | higher-arity flags with root-to-leaf factors |
| `weighted_flag_*` | higher-arity flags with additional root-root factors (§10.4, §19.6) |
| `harmonic_*` | ordinary nonnegative spherical-harmonic moments, from degree 2 upward (`harmonic_2` encodes $p_2\ge1/3$) |
| `harmonic_flag_*`, `spin2_flag` | spin-$\ell$ Gram blocks of harmonic-weighted unrooted flags; the spin-2 block contains the deviatoric second moment and subsumes the removed isotropy contraction |
| `hessian_*` | scalar and matrix-valued KKT Hessian multipliers (§19.10) |
| `global_*_gap` | nonnegative multipliers of the global first-variation gap (§19.11) |

Free (equality) families, which need not be PSD:

| family | mathematics |
|---|---|
| `potential_flag_*` | $U_\mu(x)-E(\mu)=0$ on the support (§19.8) |
| gradient relations | $\nabla_{S^2}U_\mu(x)=0$ (§19.9) |
| `rank_flag_*`, scalar rank relations | $\det\operatorname{Gram}_4=0$ (§19.12) |
| label reductions | exact exchangeability built into the labels |

## 18. CLI reference

### Degree and arity

- `--degree D` bounds the total Gram-polynomial degree.
- `--max-flag-arity A` adds systematic conditional squares through $A$ sampled vertices after gluing.  Five-point constraints begin at arity 5; six-point at arity 6.
- `--max-root-factor-degree R` permits total root-root degree at most $R$ in the higher-arity flags.  The strongest standard run uses $A=5$, $R=2$.
- `--max-hessian-arity A` adds matrix-valued Hessian multipliers with additional shared conditioning roots.

### Equality-face audit and facial reduction

- `--check-onb` evaluates every PSD block and equality relation on the uniform orthonormal-basis measure.  The target must be zero, PSD blocks must have no negative eigenvalues, equality blocks must vanish.
- `--facial-reduce-onb` parametrizes each primal PSD matrix directly on the nullspace forced by equality at the ONB; valid only for the sharp target with `--target-epsilon 0`.
- `--exact-onb-face` reconstructs that nullspace over $\mathbb Q$, avoiding a numerical eigenvalue cutoff.
- `--pole-equator-faces continuous` applies a second exact equality face; its moment functional is computed by Fourier constant-term arithmetic.  For the underlying measure

$$\mu=\frac13\delta_e+\frac23\nu_{e^\perp},\qquad U_\mu(z)=4q(1-q)^2,\quad q=(z\cdot e)^2,$$

so the global KKT gap and its support Hessian are nonnegative without assuming the desired theorem.
- `--eliminate-free` forms the exact rational quotient by every unrestricted potential, gradient, and Gram-rank multiplier.  For the degree-14, arity-5, root-factor-degree-2 model, 74 generated equality columns have rank 71 among 574 labels, leaving 503 independent coefficient equations.
- `--find-face` and `--numerical-faces` are facial-reduction diagnostics only: a numerically inferred face is never proof data; any certificate found with it must be embedded in the original basis and verified exactly.

### Conditioning and output

- `--scale-constraints` rescales mathematically equivalent constraints before the solve: each PSD or equality matrix is divided by $s=\max_L\|A_L\|_{\max}$, each scalar identity by its largest coefficient.  Since $s>0$ this changes neither PSD signs nor equality solution sets, only conditioning — essential near the singular equality face.  Different scalings produced substantially different floating-point answers in early runs, so no objective is accepted without residual and nesting checks.
- `--summary-only` suppresses the formal-moment dictionary; the JSON still reports solver status, `objective`, number of labels, `minimum_block_eigenvalue`, `maximum_free_residual` (matrix-equality residual), and `maximum_relation_residual`.
- `--tolerance 1e-10` sets MOSEK primal-feasibility, dual-feasibility, and relative-gap tolerances — a numerical termination threshold, not an exact-arithmetic guarantee.

### The reference command

```sh
python sos_search.py \
  --dual --summary-only --scale-constraints --rank-relations \
  --degree 14 --no-pointwise-sos \
  --harmonics --three-point-flags --four-point-flags --two-root-flags \
  --max-flag-arity 5 --max-root-factor-degree 2 \
  --gradient --potential --potential-matrices \
  --hessian --four-point-hessian \
  --global-gap --global-tangent-gaps \
  --tolerance 1e-10
```

## 19. Degree-14 command: complete family inventory

**Branch note.**  The inventory below was derived for the *isotropic* formulation, in which every degree-two sampled vertex is contracted through $\mathbb E[XX^{\mathsf T}]=\tfrac13I$ (§9.3), the target is (3.3), and the degree-two harmonic vanishes.  The current non-isotropic code keeps degree-two vertices as labels and uses target (3.1); label counts, block lists, and the reported bound differ accordingly (see [SUBCASES_AND_RECORD.md](SUBCASES_AND_RECORD.md) for the non-isotropic record).  The derivations below remain correct for the isotropic branch they were written for, and the kernel formulas are branch-independent.

### 19.1 Target and label count

Under isotropy $p_2=\mathbb E_Y[Y^{\mathsf T}\mathbb E[XX^{\mathsf T}]Y]=\frac13$, so the SDP target is

$$\boxed{E=\frac{16}{3}-48p_4+32p_6}$$

with coefficient vector `constant: 16/3`, `pair_4: -48`, `pair_6: 32`.  The command produces **574** distinct canonical labels, one fixed by $y_{\mathrm{constant}}=1$.

### 19.2 Degree budget (`--degree 14`, `--no-pointwise-sos`)

Glued Gram monomials have total degree $\le14$: each half of an ordinary square has degree $\le7$; for a localizing multiplier of degree $d$ the square basis satisfies $2\deg(f)+d\le14$.  Mixed-degree geometric kernels use their conservative maximum degree. `--no-pointwise-sos` disables the ordinary polynomial SOS block $v(a,b,c)^{\mathsf T}Qv(a,b,c)$ in the three Gram variables of three points; the command also omits `--gram-module` (no generic three-point Gram-determinant/principal-minor quadratic module) and `--higher-rank-matrices` (no matrix-valued Gram-determinant identities on the five-point flag spaces).  All conditional flag squares remain active.

### 19.3 `--three-point-flags`: one-root $O(2)$ harmonic flags

Root $X$, leaves $Y,Z$; $a=X\cdot Y$, $b=Y\cdot Z$, $c=Z\cdot X$. Tangent projections $P_XY=Y-aX$, $P_XZ=Z-cX$ give $(P_XY)\cdot(P_XZ)=b-ac$ and $\|P_XY\|^2\|P_XZ\|^2=(1-a^2)(1-c^2)$.  The stabilizer of $X$ is $O(2)$; its weight-$k$ zonal kernel is the polynomial $R_k(a,b,c)=\bigl((1-a^2)(1-c^2)\bigr)^{k/2}\cos(k(\phi_Y-\phi_Z))$, computed by the square-root-free recurrence

$$R_0=1,\qquad R_1=b-ac,\qquad R_{k+1}=2(b-ac)R_k-(1-a^2)(1-c^2)R_{k-1}.$$

For leaf powers $r,s$ with parity $r\equiv k\pmod2$ (correct projective character), the dual imposes $\big(\mathbb E[a^rc^sR_k(a,b,c)]\big)_{r,s}\succeq0$.  Degree-14 blocks:

| block | $k$ | leaf powers | size |
|---|---:|---|---:|
| `flag_0` | 0 | $0,2,4,6$ | 4 |
| `flag_1` | 1 | $1,3,5$ | 3 |
| `flag_2` | 2 | $0,2,4$ | 3 |
| `flag_3` | 3 | $1,3$ | 2 |
| `flag_4` | 4 | $0,2$ | 2 |
| `flag_5` | 5 | $1$ | 1 |
| `flag_6` | 6 | $0$ | 1 |

### 19.4 `--four-point-flags`: the empty-type pair block

For even $r$, the unrooted two-sample flag average is $f_r=\mathbb E_{X,Y}[(X\cdot Y)^r]=p_r$; squaring gives $\mathbb E[(X_0\cdot X_1)^r(X_2\cdot X_3)^s]=p_rp_s$, so

$$\big(\mathbb E[(X_0\cdot X_1)^r(X_2\cdot X_3)^s]\big)_{r,s}\succeq0.$$

At degree 14 the basis is $r\in\{0,2,4,6\}$: `empty_type_flag` is a $4\times4$ PSD block connecting single-pair moments to disconnected `product` labels.

### 19.5 `--two-root-flags`: four-point conditional squares

Roots $X,Z$, one leaf $Y$; basis flag $f_{ijk}(X,Z;Y)=(X\cdot Y)^i(Z\cdot Y)^j(X\cdot Z)^k$, $i+j+k\le7$. Gluing a second leaf $W$:

$$f_{ijk}(X,Z;Y)f_{i'j'k'}(X,Z;W) =(X\cdot Z)^{k+k'}(X\cdot Y)^i(X\cdot W)^{i'}(Z\cdot Y)^j(Z\cdot W)^{j'}.$$

The three parity characters of a flag are $(i+j,\ i+k,\ j+k)\pmod2$, corresponding to leaf $Y$, root $X$, root $Z$ respectively.  Four sectors are used.

**Even-leaf sectors** (direct conditional squares), each with 30 basis flags at degree 14:

| sector | $i+j$ | $i+k$ | $j+k$ |
|---|---:|---:|---:|
| `two_root_even_00` | 0 | 0 | 0 |
| `two_root_even_11` | 0 | 1 | 1 |

**Odd-leaf orientation sectors.**  If $i+j$ is odd the raw conditional leaf integral vanishes by antipodal symmetry; the code pairs it with the oriented tangent factor $\det(X,Z,Y)$.  After gluing the multiplier is $\det(X,Z,Y)\det(X,Z,W)$, independent of orientation. With $r=X\cdot Z$, $p=X\cdot Y$, $q=X\cdot W$, $s=Z\cdot Y$, $t=Z\cdot W$, $u=Y\cdot W$, the exact Gram polynomial is

$$u-r^2u-st-pq+rpt+rqs.$$

| sector | $i+j$ | $i+k$ | $j+k$ |
|---|---:|---:|---:|
| `two_root_odd_01` | 1 | 0 | 1 |
| `two_root_odd_10` | 1 | 1 | 0 |

Each has 30 basis flags before localization.

**Root-pair localizing blocks.**  Since $X,Z\in S^2$, $1-(X\cdot Z)^2\ge0$; a second PSD block in every sector is multiplied by this nonnegative principal minor (odd sectors additionally carry $\det(X,Z,Y)\det(X,Z,W)$), with square basis of degree $\le6$:

| block | size |
|---|---:|
| `two_root_even_00_minor` | 24 |
| `two_root_even_11_minor` | 20 |
| `two_root_odd_01_minor` | 20 |
| `two_root_odd_10_minor` | 20 |

### 19.6 Five-point root-weighted flags (`--max-flag-arity 5 --max-root-factor-degree 2`)

Three shared roots $R=(X_0,X_1,X_2)$ and one leaf $Y$ per half.  A basis flag is indexed by root-edge exponents $\rho=(\rho_{01},\rho_{02},\rho_{12})$ and root-to-leaf exponents $\lambda=(\lambda_0,\lambda_1,\lambda_2)$, representing

$$f_{\rho,\lambda}(R;Y) =\prod_{0\le i<j\le2}(X_i\cdot X_j)^{\rho_{ij}} \prod_{i=0}^2(X_i\cdot Y)^{\lambda_i},$$

with restrictions

$$|\rho|\le2\ (\texttt{--max-root-factor-degree 2}),\qquad |\rho|+|\lambda|\le7\ (\text{glued degree}\le14),\qquad |\lambda|\equiv0\ (\mathrm{mod}\ 2)\ (\text{even leaf}).$$

Parity character at root $X_i$: $\sigma_i=\lambda_i+\sum_{j\ne i}\rho_{\min(i,j),\max(i,j)}\pmod2$; only flags with equal $\sigma$ pair.  Active signatures and sizes:

| block | root signature | size |
|---|---|---:|
| `weighted_flag_5_000` | $(0,0,0)$ | 92 |
| `weighted_flag_5_011` | $(0,1,1)$ | 80 |
| `weighted_flag_5_101` | $(1,0,1)$ | 80 |
| `weighted_flag_5_110` | $(1,1,0)$ | 80 |

Gluing row $(\rho,\lambda)$ against column $(\rho',\lambda')$ with leaves $Y,W$ produces $\prod_{i<j}(X_i\cdot X_j)^{\rho_{ij}+\rho'_{ij}} \prod_i(X_i\cdot Y)^{\lambda_i}\prod_i(X_i\cdot W)^{\lambda'_i}$, reduced as a five-vertex graph.  No six-point block is added at arity 5.

### 19.7 `--harmonics`: two-point spherical-harmonic positivity

The Legendre addition formula gives

$$A_\ell=\mathbb E[P_\ell(X\cdot Y)] =\frac{4\pi}{2\ell+1}\sum_{m=-\ell}^{\ell} \left|\int Y_{\ell m}(x)\,d\mu(x)\right|^2\ \ge\ 0,$$

(normalization convention-dependent, sign not).  Each $A_\ell$ is a $1\times1$ PSD block.  At degree 14 the command adds $\ell=4,6,8,10,12,14$ ($\ell=0$ is the normalization $A_0=1$; under isotropy the degree-two harmonic vanishes and is not added).  The exact polynomials:

$$\begin{aligned} P_4(t)&=\frac{35t^4-30t^2+3}{8},\\ P_6(t)&=\frac{231t^6-315t^4+105t^2-5}{16},\\ P_8(t)&=\frac{6435t^8-12012t^6+6930t^4-1260t^2+35}{128},\\ P_{10}(t)&=\frac{46189t^{10}-109395t^8+90090t^6-30030t^4+3465t^2-63}{256},\\ P_{12}(t)&=\frac{676039t^{12}-1939938t^{10}+2078505t^8-1021020t^6+225225t^4-18018t^2+231}{1024},\\ P_{14}(t)&=\frac{5014575t^{14}-16900975t^{12}+22309287t^{10}-14549535t^8+4849845t^6-765765t^4+45045t^2-429}{2048}. \end{aligned}$$

### 19.8 `--potential` and `--potential-matrices`: support-potential identities

Support stationarity $U_\mu(X)=E(\mu)$ $\mu$-a.s. (8.2), multiplied by $(X\cdot Z)^r$ and averaged over independent $X,Y,Z,W\sim\mu$:

$$\boxed{\mathbb E[K(X\cdot Y)(X\cdot Z)^r] -\mathbb E[K(X\cdot Y)]\,\mathbb E[(Z\cdot W)^r]=0.}$$

Equality constraints; multipliers unrestricted.  At degree 14 the even powers $r=0,2,4,6,8$ are tested; after antipodal and isotropic reduction the $r=0,2$ relations are identically zero, leaving three nonzero relations $r=4,6,8$.  E.g. $r=4$ reduces to

$$0=-48M_{\triangle(0,4,4)}+32M_{\triangle(0,4,6)} +48p_4^2-32p_4p_6.$$

`--potential-matrices` strengthens the scalar equalities by multiplying $U_\mu(X)-E(\mu)$ by the one-root flag Gram matrices of §19.3; a typical entry is

$$0=\mathbb E\!\left[K(X\cdot Y)(X\cdot Z)^r(X\cdot W)^sR_k(X;Z,W)\right] -E(\mu)\,\mathbb E\!\left[(X\cdot Z)^r(X\cdot W)^sR_k(X;Z,W)\right],$$

with $E=\frac{16}{3}-48p_4+32p_6$ substituted in the isotropic branch. These are matrix *equalities* (free multipliers).  Degree-14 bases (shorter than §19.3 because the kernel consumes degree six):

| free block | $k$ | leaf powers | size |
|---|---:|---|---:|
| `potential_flag_0` | 0 | $0,2,4$ | 3 |
| `potential_flag_1` | 1 | $1,3$ | 2 |
| `potential_flag_2` | 2 | $0,2$ | 2 |
| `potential_flag_3` | 3 | $1$ | 1 |
| `potential_flag_4` | 4 | $0$ | 1 |

### 19.9 `--gradient`: first-order stationarity

With $a=X\cdot Y$, $b=Y\cdot Z$, $c=Z\cdot X$, the tangent projection of $Z$ at $X$ satisfies $(Y-aX)\cdot(Z-cX)=b-ac$, so the scalar gradient kernel is

$$\boxed{g(a,b,c)=K'(a)(b-ac).}$$

At a KKT support point $G_\mu(X)=\nabla_{S^2}U_\mu(X)=0$; multiplying by $c^r$ and averaging gives $\mathbb E[c^rK'(a)(b-ac)]=0$.  For $r=0,\ldots,7$, antipodal parity and isotropic contraction annihilate all but $r=3,5,7$: three unrestricted scalar equalities.  E.g. $r=3$:

$$0=-192M_{\triangle(0,4,6)}+192M_{\triangle(1,3,5)} +192M_{\triangle(0,4,4)}-192M_{\triangle(1,3,3)}.$$

### 19.10 `--hessian` and `--four-point-hessian`: Hessian positivity

For tangent vectors $v,w\perp X$ the spherical Hessian is

$$H_X(v,w)=\int\big[K''(a)(v\cdot Y)(w\cdot Y)-K'(a)a(v\cdot w)\big]\,d\mu(Y), \qquad a=X\cdot Y,$$

with $H_X\succeq0$ at a locally minimizing support point.

**Scalarized blocks (`--hessian`).**  For $v_\parallel=Z-cX$ ($c=X\cdot Z$): $v_\parallel\cdot Y=b-ac$, $\|v_\parallel\|^2=1-c^2$, giving the parallel kernel

$$\boxed{h_\parallel(a,b,c)=K''(a)(b-ac)^2-K'(a)a(1-c^2).}$$

For $v_\perp=X\times Z$: $(v_\perp\cdot Y)^2=D(a,b,c)$ where $D(a,b,c)=1+2abc-a^2-b^2-c^2=\det\operatorname{Gram}(X,Y,Z)$, giving

$$\boxed{h_\perp(a,b,c)=K''(a)D(a,b,c)-K'(a)a(1-c^2).}$$

Each nonnegative rooted Hessian quantity is integrated against a polynomial nonnegative on $[-1,1]$, represented as $s(c)^{\mathsf T}Q_0s(c)+(1-c^2)t(c)^{\mathsf T}Q_1t(c)$ with $Q_0,Q_1\succeq0$.  At degree 14: $s(c)=(1,c,c^2,c^3)$ gives the $4\times4$ blocks `hessian_sos` and `perpendicular_hessian_sos`; $t(c)=(1,c,c^2)$ gives the $3\times3$ blocks `hessian_minor` and `perpendicular_hessian_minor` (multiplied by $1-c^2$).

**Bilinear four-point blocks (`--four-point-hessian`).**  With $a=X\cdot Y$, $c=X\cdot Z$, $d=X\cdot W$, $b=Y\cdot Z$, $e=Y\cdot W$, $f=Z\cdot W$ and tangent fields $v_Z=Z-cX$, $v_W=W-dX$:

$$\boxed{H_\parallel=K''(a)(b-ac)(e-ad)-K'(a)a(f-cd).}$$

For the rotated fields $Jv_Z=X\times Z$, $Jv_W=X\times W$, the identity $(Jv_Z\cdot Y)(Jv_W\cdot Y)=(f-cd)(1-a^2)-(b-ac)(e-ad)$ gives

$$\boxed{H_\perp=K''(a)\big[(f-cd)(1-a^2)-(b-ac)(e-ad)\big]-K'(a)a(f-cd).}$$

The vector-field basis uses odd auxiliary powers $c^rv_Z$, $r\in\{1,3\}$ (odd powers make the field invariant under $Z\mapsto-Z$), creating the two $2\times2$ PSD blocks `four_point_parallel_hessian` and `four_point_perpendicular_hessian`. Without `--max-hessian-arity`, no five-point Hessian flag blocks are included.

### 19.11 `--global-gap` and `--global-tangent-gaps`

**Uniform gap.**  Averaging $q_\mu(z)\ge0$ (8.1) uniformly over $S^2$, with degree-zero Legendre coefficient $\int_{S^2}K(z\cdot y)\,d\sigma(z)=\frac{32}{105}$:

$$\frac{32}{105}-E(\mu)\ \ge\ 0, \qquad\text{under isotropy}\quad \frac{32}{105}-E=-\frac{176}{35}+48p_4-32p_6.$$

This is the $1\times1$ PSD block `global_uniform_gap`.

**Tangent trial points.**  With $a=X\cdot Y$, $b=Y\cdot Z$, $c=Z\cdot X$, $d=1-c^2$, evaluate the gap at two trial points built from a support point $X$ and an auxiliary sample $Z$.  Parallel trial $z_\parallel=(Z-cX)/\sqrt d$ gives $z_\parallel\cdot Y=(b-ac)/\sqrt d$ and, using $U_\mu(X)=E$ on the support, the cleared-denominator gap

$$\boxed{G_\parallel(a,b,c)=d^3\left[K\!\left(\frac{b-ac}{\sqrt d}\right)-K(a)\right]} =32(b-ac)^6-48(b-ac)^4d+20(b-ac)^2d^2-\frac43d^3-d^3K(a),$$

a polynomial since $K$ is even of degree six.  Perpendicular trial $z_\perp=(X\times Z)/\sqrt d$ gives $(z_\perp\cdot Y)^2=D(a,b,c)/d$ and

$$\boxed{G_\perp=32D^3-48D^2d+20Dd^2-\frac43d^3-d^3K(a).}$$

Both gaps are nonnegative after averaging over $Y$; the cleared polynomial expressions extend continuously through $d=0$.  The base gap degree budget is 12; at total degree 14 the ordinary multiplier basis is $(1,c)$ and a localizing block uses multiplier $1-c^2$ with constant basis:

| block | size |
|---|---:|
| `global_parallel_tangent_gap` | 2 |
| `global_parallel_tangent_gap_minor` | 1 |
| `global_perpendicular_tangent_gap` | 2 |
| `global_perpendicular_tangent_gap_minor` | 1 |

### 19.12 `--rank-relations`: dimension-three Gram identities

Four vectors in $\mathbb R^3$ are linearly dependent, so

$$\det\begin{pmatrix} 1&g_{01}&g_{02}&g_{03}\\ g_{01}&1&g_{12}&g_{13}\\ g_{02}&g_{12}&1&g_{23}\\ g_{03}&g_{13}&g_{23}&1 \end{pmatrix}=0.$$

The determinant has Gram degree at most four; to stay within total degree 14 it is multiplied by canonical four-vertex Gram monomials of degree at most $14-4=10$.  Only multipliers with even degree at every vertex survive projective parity.  Results are reduced (antipodality, isotropy), canonicalized under all four-vertex permutations, normalized and deduplicated: **61 independent scalar relations** for the stated command, with unrestricted primal coefficients; the dual requires every formal pairing with the determinant relation to vanish. These relations are crucial: the formal moment labels otherwise forget that the sampled vectors live in dimension three.

### 19.13 Complete degree-14 inventory

Positive semidefinite blocks (**37 total**):

| family | blocks | sizes |
|---|---:|---|
| one-root harmonic flags | 7 | $4,3,3,2,2,1,1$ |
| empty-type four-point flag | 1 | $4$ |
| two-root flags | 8 | $30,30,30,30,24,20,20,20$ |
| five-point root-weighted flags | 4 | $92,80,80,80$ |
| Legendre harmonics | 6 | six $1\times1$ |
| uniform global gap | 1 | $1$ |
| tangent global gaps | 4 | $2,1,2,1$ |
| scalar Hessian localizers | 4 | $4,4,3,3$ |
| bilinear four-point Hessians | 2 | $2,2$ |

Matrix equality blocks: the five potential-stationarity flag matrices, sizes $3,2,2,1,1$.

Scalar equalities (**67 total**):

| family | count |
|---|---:|
| gradient stationarity | 3 |
| scalar potential stationarity | 3 |
| four-vector rank identities | 61 |

Formal variables: **574** moment labels, one fixed by $y_{\mathrm{constant}}=1$.

### 19.14 Option-to-implementation map

| option | main implementation |
|---|---|
| `--degree 14` | degree budgets inside `solve` |
| `--no-pointwise-sos` | skips the `sos` module term |
| `--three-point-flags` | `tangent_harmonic_polynomials`, `flag_expectation_matrix` |
| `--four-point-flags` | `empty_type_flag_expectation_matrix` |
| `--two-root-flags` | `two_root_flag_expectation_matrix` |
| `--max-flag-arity 5` | `rooted_weighted_flag_sectors`, `rooted_weighted_flag_expectation_matrix` |
| `--max-root-factor-degree 2` | root-edge degree bound in the weighted flag basis |
| `--harmonics` | `harmonic_pair_vector` |
| `--potential` | `potential_stationarity_relation` |
| `--potential-matrices` | `potential_flag_relation_matrix` |
| `--gradient` | gradient polynomial from `kernel_polynomials` |
| `--hessian` | scalar kernels from `kernel_polynomials` |
| `--four-point-hessian` | `four_point_hessian_polynomials`, `four_point_hessian_expectation_matrix` |
| `--global-gap` | `global_uniform_gap` block in `solve` |
| `--global-tangent-gaps` | `global_tangent_gap_polynomials` |
| `--rank-relations` | `four_point_rank_relations` |
| `--dual` | formal moment branch in `solve` |
| `--scale-constraints` | positive block/relation rescaling in the dual branch |
| `--summary-only` | omits the formal moment dictionary from JSON |
| `--tolerance 1e-10` | MOSEK interior-point tolerance parameters |

### 19.15 Interpreting a reported objective

A small negative dual optimum means the selected finite formal-moment relaxation nearly proves $E\ge0$; the residual may be numerical error, a relaxation gap, or a non-attainment/closure phenomenon.  It is not an actual measure of negative energy, and it is not an exact SOS proof.  A rigorously nonnegative dual optimum would prove $E\ge0$ for every actual measure in the encoded class; a recovered rational primal identity with exact PSD matrices would give a machine-checkable certificate for that class.  For the isotropic branch specifically, turning a branch result into unrestricted copositivity requires either (1) a proof that a hypothetical negative global minimizer can be assumed isotropic, or (2) a hierarchy that retains the full second-moment matrix — which is what the current non-isotropic code does.

---

# Appendix — Operational conventions

## A.1 SDPA-GMP build

Double-precision solvers provably stall on these problems; SDPA-GMP (arbitrary-precision interior-point SDP solver) is required.  Build from source, `nakatamaho/sdpa-gmp`, with GMP via Homebrew.  Two build fixes: `-Wno-int-conversion` for the bundled SPOOLES C code, `-std=c++17` for the C++ part.  The binary, parameter files, and all problem files live in `sdpa_runs/` (gitignored).

## A.2 Invocation and parameter files

Invoke with option-style flags only — positional arguments are misparsed:

```sh
./sdpa_gmp -ds PROBLEM.dat-s -o OUT.result -p param_200bit.sdpa
```

Parameter file conventions:

- `param_200bit.sdpa` — 200-bit precision, `epsilonStar 1e-25`: for bound measurements, where the primal–dual agreement itself is the datum (typically ≥13 digits).
- `param_128bit.sdpa` — 128-bit precision, `epsilonStar 1e-16`: for large problems and selector solves, where full precision is unaffordable and only feasibility/trace behavior is needed.

Verify completion from the `.result` file (`phase.value`, feasibility errors), never from mid-run status.

## A.3 Export pipeline and objective shifts

`sos_search.py --export-sdpa FILE.dat-s` writes exact rational SDPA problem files.  The **reported bound = objValPrimal + objective_shift**, where `objective_shift` is printed in the export JSON:

| target | shift |
|---|---|
| $E$ | $-4/3$ |
| $h_2E$ | $+2/3$ |

**Legacy scale.**  Exports made under the older convention used the $(3/16)$-scaled target, with shifts $-1/4$ (for scaled $E$) and $+1/8$ (for scaled $h_2E$).  Ratios and pole orders are scale-free, but any absolute bound quoted from a legacy export must be rescaled before comparison.

Selector problems (minimal certificate trace subject to a bound on the objective) are exported by `sdpa_selector.py`; exact expansion JSONs and ray data accompany the problems in `sdpa_runs/`.  Certificate candidates are checked by `verify_certificate.py` and structural identities by `verify_exact_structure.py` (§15).

## A.4 Baseline reproduction

The pruned nine-toggle hierarchy at degree 14,

```sh
python3 sos_search.py --export-sdpa deg14_pruned.dat-s \
  --degree 14 --no-pointwise-sos \
  --harmonics --three-point-flags --four-point-flags --two-root-flags \
  --gradient --potential --hessian --global-tangent-gaps --rank-relations
sdpa_gmp -ds deg14_pruned.dat-s -o deg14_pruned.result -p param_200bit.sdpa
```

reproduces the published bound to 10 significant digits: $-4.48560259\times10^{-3}$ (m = 398, 46 iterations, ~6 min, feasibility errors $\sim10^{-26}$).  For the weighted-target and enrichment measurements built on this baseline, see [ENRICHMENTS.md](ENRICHMENTS.md) and [SUBCASES_AND_RECORD.md](SUBCASES_AND_RECORD.md).
