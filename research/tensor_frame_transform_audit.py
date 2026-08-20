"""Numerical identities and inequalities for the ONB frame transform T_mu."""

import numpy as np

from tensor_pair_circle_projection import polynomial_root_average, energy_data


def frame_transform(points, weights):
    transformed_points = []
    transformed_weights = []
    for x, wx in zip(points, weights):
        for y, wy in zip(points, weights):
            a = x @ y
            s = 1 - a * a
            if s < 1e-14:
                continue
            u = (y - a * x) / np.sqrt(s)
            n = np.cross(x, y) / np.sqrt(s)
            w = a*a*s*s / 4
            transformed_points.extend((x, u, n))
            transformed_weights.extend((8 * wx * wy * w,) * 3)
    return np.asarray(transformed_points), np.asarray(transformed_weights)


def j_cross(first_points, first_weights, second_points, second_weights):
    gram = first_points @ second_points.T
    kernel = 144 * gram**6 - 216 * gram**4 + 87 * gram**2 - 5
    return first_weights @ kernel @ second_weights


def scan(samples=100, seed=20260822):
    rng = np.random.default_rng(seed)
    records = []
    for _ in range(samples):
        points = rng.normal(size=(7, 3))
        points /= np.linalg.norm(points, axis=1)[:, None]
        weights = rng.dirichlet(np.ones(7))
        transformed_points, transformed_weights = frame_transform(points, weights)
        J, F, gap = energy_data(points, weights)
        A = polynomial_root_average(points, weights, 2)
        cross = j_cross(points, weights, transformed_points, transformed_weights)
        self_energy = j_cross(
            transformed_points, transformed_weights,
            transformed_points, transformed_weights,
        )
        H3 = J-108*F-288*A-3*self_energy
        records.append((J, F, A, cross, self_energy, H3))
    records = np.asarray(records)
    print("cross identity", np.max(np.abs(records[:, 3] - 108 * records[:, 1])))
    print("A/JTT ratio", np.min(records[:, 2] / records[:, 4]),
          np.max(records[:, 2] / records[:, 4]))
    design = np.c_[records[:, 4], records[:, 3] ** 2 / records[:, 0]]
    coefficient, *_ = np.linalg.lstsq(design, records[:, 2], rcond=None)
    print("fits", coefficient,
          np.linalg.norm(design @ coefficient - records[:, 2]) / np.linalg.norm(records[:, 2]))
    print("H3 min/max", np.min(records[:, 5]), np.max(records[:, 5]))
    print("(G-288A)/JTT min", np.min((records[:, 0]-records[:, 3]-288*records[:, 2])/records[:, 4]))
    return records


if __name__ == "__main__":
    scan()
