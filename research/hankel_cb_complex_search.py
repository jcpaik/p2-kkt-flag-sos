"""Search rank-seven ternary-sextic Hankel rays with one complex point pair.

Take a transverse intersection of F=xyz with a real cubic G.  Two coordinate
lines are required to meet G in three real points, while the third contributes
one real point and one complex-conjugate pair.  Cubic evaluations then form a
real hyperplane in R^7 + C.  A quadratic form

    sum_i a_i r_i^2 + 2 Re(a_c z^2)

with a_i>0 has one negative ambient direction.  Its restriction to the CB
hyperplane is PSD and singular exactly when

    sum_i u_i^2/a_i + 2 Re(u_c^2/a_c) = 0.

This produces the missing complex-pair family of rank-seven extreme Hankel
rays from Blekherman's classification.
"""
import math
import numpy as np


def comps(n, k=3):
    if k == 1:
        return [(n,)]
    return [(a,) + rest for a in range(n + 1) for rest in comps(n - a, k - 1)]


def multinomial(a):
    out = math.factorial(sum(a))
    for x in a:
        out //= math.factorial(x)
    return out


d3 = comps(3)
d6 = comps(6)
ix6 = {a: i for i, a in enumerate(d6)}


def moment_map(d):
    basis = comps(d)
    out = np.zeros((len(basis), len(basis), len(d6)))
    for i, a in enumerate(basis):
        for j, b in enumerate(basis):
            scale = np.sqrt(multinomial(a) * multinomial(b))
            for g in comps(3 - d):
                exponent = tuple(a[k] + b[k] + 2 * g[k] for k in range(3))
                out[i, j, ix6[exponent]] += scale * multinomial(g)
    return out


maps = {d: moment_map(d) for d in (1, 2, 3)}
grams = {d: np.einsum("abk,abl->kl", M, M) for d, M in maps.items()}
mass_vector = np.trace(maps[3], axis1=0, axis2=1)
Q = 32 * grams[3] - 48 * grams[2] + 20 * grams[1] - 4 / 3 * np.outer(mass_vector, mass_vector)


def evaluate(points, basis):
    return np.array([[np.prod(p ** np.array(a)) for a in basis] for p in points])


def roots_on_coordinate_lines(coefficients):
    # Coefficient order is d3.  Return roots t for points (0,t,1),
    # (t,0,1), and (t,1,0), respectively.
    coeff = {a: coefficients[i] for i, a in enumerate(d3)}
    polynomials = [
        [coeff[(0, 3, 0)], coeff[(0, 2, 1)], coeff[(0, 1, 2)], coeff[(0, 0, 3)]],
        [coeff[(3, 0, 0)], coeff[(2, 0, 1)], coeff[(1, 0, 2)], coeff[(0, 0, 3)]],
        [coeff[(3, 0, 0)], coeff[(2, 1, 0)], coeff[(1, 2, 0)], coeff[(0, 3, 0)]],
    ]
    return [np.roots(p) for p in polynomials]


def intersection(rng):
    while True:
        coefficients = rng.normal(size=10)
        roots = roots_on_coordinate_lines(coefficients)
        real_counts = [sum(abs(z.imag) < 1e-7 for z in rr) for rr in roots]
        if sorted(real_counts) != [1, 3, 3]:
            continue
        points_real = []
        complex_point = None
        for line, rr in enumerate(roots):
            for z in rr:
                if line == 0:
                    p = np.array([0, z, 1], dtype=complex)
                elif line == 1:
                    p = np.array([z, 0, 1], dtype=complex)
                else:
                    p = np.array([z, 1, 0], dtype=complex)
                # Use real projective normalization for real points and a
                # conjugation-compatible positive real scaling for the pair.
                p /= np.sqrt(np.sum(np.abs(p) ** 2))
                if abs(z.imag) < 1e-7:
                    points_real.append(p.real)
                elif z.imag > 0:
                    complex_point = p
        if len(points_real) != 7 or complex_point is None:
            continue
        points = np.array(points_real + [complex_point, complex_point.conjugate()])
        V = evaluate(points, d3)
        _, singular, Vh = np.linalg.svd(V.T)
        if singular[-1] > 2e-8:
            continue
        u = Vh[-1].conjugate()
        # Fix global phase so the real-point coefficients are real.
        phase = np.angle(np.sum(u[:7] ** 2)) / 2
        u *= np.exp(-1j * phase)
        if np.max(np.abs(u[:7].imag)) > 2e-6:
            continue
        u_real = u[:7].real
        u_complex = (u[7] + u[8].conjugate()) / 2
        if np.min(np.abs(u_real)) < 1e-7 or abs(u_complex) < 1e-7:
            continue
        return points_real, complex_point, u_real, u_complex, coefficients


def ray(rng):
    while True:
        real_points, z, u, uc, coefficients = intersection(rng)
        positive = np.exp(rng.normal(scale=1.8, size=7))
        C = np.sum(u * u / positive)
        tau = C * rng.normal(scale=2.0)
        ac = uc * uc / (-C / 2 + 1j * tau)
        y = positive @ evaluate(real_points, d6) + 2 * np.real(ac * evaluate([z], d6)[0])
        mass = mass_vector @ y
        if mass <= 1e-8:
            continue
        y /= mass
        H = np.einsum("abk,k->ab", maps[3], y)
        eig = np.linalg.eigvalsh(H)
        if eig[0] < -2e-6 or eig[3] < 1e-8:
            continue
        return y, (real_points, z, u, uc, positive, ac, coefficients, eig)


rng = np.random.default_rng(77491)
N = 16000
Y = []
metadata = []
for k in range(N):
    y, meta = ray(rng)
    Y.append(y)
    metadata.append(meta)
    if (k + 1) % 1000 == 0:
        print("sampled", k + 1, flush=True)
Y = np.asarray(Y)
self_q = np.einsum("bi,ij,bj->b", Y, Q, Y)
i0 = int(np.argmin(self_q))
print("self min/max", self_q[i0], self_q.max(), "index", i0)
print("self meta", metadata[i0])

# Pair search includes cross terms between complex-pair rays.
best = (np.inf, None, None, None)
for lo in range(0, N, 250):
    cross = Y[lo : lo + 250] @ Q @ Y.T
    score = cross + np.sqrt(np.maximum(self_q[lo : lo + 250, None] * self_q[None, :], 0))
    for k in range(score.shape[0]):
        score[k, lo + k] = np.inf
    a, b = np.unravel_index(np.argmin(score), score.shape)
    if score[a, b] < best[0]:
        best = (score[a, b], lo + a, b, cross[a, b])
_, i, j, cross = best
print("pair", i, j, "Lorentz score", best[0], "q", self_q[i], self_q[j], "cross", cross)
d = Y[j] - Y[i]
aa = d @ Q @ d
bb = 2 * Y[i] @ Q @ d
t = np.clip(-bb / (2 * aa), 0, 1) if aa > 0 else (1 if self_q[j] < self_q[i] else 0)
mix = Y[i] + t * d
print("mix", t, mix @ Q @ mix)
print("meta1", metadata[i])
print("meta2", metadata[j])
