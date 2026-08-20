# Sharp zero-face one-atom chord inequality

Let

\[
K(t)=32t^6-48t^4+20t^2-\frac43
\]

and, for a probability measure \(\mu\) on \(\mathbb {RP}^2\), let
\(E(\mu)=\mathbb E K(X\cdot Y)\).  For three unit vectors put

\[
a=X\cdot Y,\qquad b=X\cdot Z,\qquad c=Y\cdot Z,
\qquad D=1+2abc-a^2-b^2-c^2
\]

and define the symmetric nonnegative kernel

\[
h(X,Y,Z)=\frac83D\sum_{\rm cyc}a^2b^2(c-ab)^2.
\]

Write \(H(\alpha,\beta,\gamma)=\int h\,d\alpha d\beta d\gamma\) and
\(F(\mu)=H(\mu,\mu,\mu)\).  The Gram determinant makes \(h\) vanish
whenever two projective arguments coincide.

## Theorem

Let

\[
\mu_0=\frac13\delta_p+\frac23\nu_{p^\perp},
\qquad \widehat\nu(2)=\widehat\nu(6)=0.
\]

For every projective point \(z\) and every \(0\le t\le1\),

\[
E((1-t)\mu_0+t\delta_z)\ge
24F((1-t)\mu_0+t\delta_z).
\]

The constant 24 is best possible uniformly over this family.

## Proof

Put \(q=(z\cdot p)^2\), \(s=1-q\), and let \(\psi\) be the azimuth of
the equatorial projection of \(z\).  In the frame with this azimuth zero,
set

\[
m=e^{-4i\psi}\widehat\nu(4),\qquad
n=e^{-8i\psi}\widehat\nu(8),\qquad r=\Re m,\qquad u=|m|.
\]

Direct Fourier expansion, using \(\widehat\nu(2)=\widehat\nu(6)=0\), gives

\[
U_0(z):=\int K(z\cdot y)d\mu_0(y)=4qs^2(1-r).
\]

It also gives

\[
H(z,\mu_0,\mu_0)=\frac49(P+Q),
\qquad P=\frac{q^2s^3(1-r)}6,
\qquad Q=\frac{qs^2}{96}L,
\]

where

\[
\begin{aligned}
L={}&2(5q^2-2q+3)-2(2q^2-4q+3)r\\
 &-(7q^2+2q+3)u^2+2\Re(n\bar m)+s^2\Re(m^2).
\end{aligned}
\]

Since \(|n|\le1\),

\[
\Re(n\bar m)\le u,\qquad \Re(m^2)\le u^2.
\]

After these substitutions, the desired inequality
\(72H(z,\mu_0,\mu_0)\le2U_0(z)\) is equivalent to

\[
L\le8(3-2qs)(1-r).
\]

The difference between the right side and the displayed upper bound for
\(L\) is affine decreasing in \(r\), so its minimum subject to \(r\le u\)
occurs at \(r=u\).  There it factors as

\[
2(u-1)\big((3q^2+2q+1)u-3q^2+6q-9\big)\ge0.
\]

Indeed, \(u-1\le0\), while the second factor is increasing in \(u\) and
is at most its value \(8(q-1)\le0\) at \(u=1\).  Hence
\(72H\le2U_0\).

For \(\mu_t=(1-t)\mu_0+t\delta_z\), collision vanishing and
\(E(\mu_0)=F(\mu_0)=0\) give the exact identities

\[
F(\mu_t)=3t(1-t)^2H(z,\mu_0,\mu_0),
\]

\[
E(\mu_t)=2t(1-t)U_0(z)+\frac83t^2.
\]

Therefore

\[
E(\mu_t)-24F(\mu_t)
=t(1-t)\big(2U_0(z)-72(1-t)H(z,\mu_0,\mu_0)\big)
+\frac83t^2\ge0.
\]

For sharpness, split an infinitesimal amount of pole mass toward the
equator and take zero-face measures with \(m\to1\) in the displacement
frame.  At fourth order in the displacement,

\[
\frac{E}{F}\longrightarrow24.
\]

More explicitly, the fourth-order coefficient of the determinant term is

\[
h_4=\frac1{18}-\frac r{108}-\frac{|m|^2}{18}
       +\frac{\Re(n\bar m)}{108}\le\frac{1-r}{9},
\]

and equality is approached as \(r,|m|,\Re(n\bar m)\to1\).

