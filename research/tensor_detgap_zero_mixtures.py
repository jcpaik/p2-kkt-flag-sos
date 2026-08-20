"""Audit the determinant-gap conjecture on mixtures of rotated zero measures."""

import numpy as np
import torch
from scipy.spatial.transform import Rotation

from tensor_detgap_opt import Q, values


torch.set_default_dtype(torch.float64)


def pole_equator(m=5):
    points = [[0.0, 0.0, 1.0]]
    weights = [1 / 3]
    for k in range(m):
        theta = np.pi * k / m
        points.append([np.cos(theta), np.sin(theta), 0.0])
        weights.append(2 / (3 * m))
    return np.array(points), np.array(weights)


def onb():
    return np.eye(3), np.ones(3) / 3


def rotation_y(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def evaluate(points, weights):
    raw = torch.tensor(points)
    logits = torch.log(torch.tensor(weights))
    return [float(x) if x.numel() == 1 else x.numpy() for x in values(raw, logits)]


def bhat(points, weights):
    x = torch.tensor(points)
    w = torch.tensor(weights)
    rho = 2 * x[:, :, None] * x[:, None, :] - torch.eye(3)
    action = torch.einsum("aij,nip,bpq,nqj->nab", Q, rho, Q, rho)
    A2 = torch.einsum("n,nab->ab", w, action)
    return ((5 / 4) * (torch.eye(5) - A2)).numpy()


def cross_energy(points1, weights1, points2, weights2):
    t = points1 @ points2.T
    K = 32 * t**6 - 48 * t**4 + 20 * t**2 - 4 / 3
    return np.sum(weights1[:, None] * weights2[None, :] * K)


def scan(base):
    points, weights = base
    B0 = bhat(points, weights)
    best = (1e9, None)
    best_ratio = (1e9, None)
    for theta in np.linspace(0, np.pi / 2, 401):
        rotated = points @ rotation_y(theta).T
        B1 = bhat(rotated, weights)
        cross = cross_energy(points, weights, rotated, weights)
        # Both endpoints have zero self-energy.
        for t in np.linspace(0, 1, 1001):
            energy = 2 * t * (1 - t) * cross
            det = np.linalg.det((1 - t) * B0 + t * B1)
            gap = energy - (32 / 105) * det**2
            if gap < best[0]:
                best = (gap, (theta, t, energy, det))
            if det > 1e-10 and energy / det**2 < best_ratio[0]:
                best_ratio = (energy / det**2, (theta, t, energy, det))
    return best, best_ratio


def scan_random_rotations(base, samples=20000, seed=20260820):
    points, weights = base
    B0 = bhat(points, weights)
    rng = np.random.default_rng(seed)
    best = (1e9, None)
    best_ratio = (1e9, None)
    for R in Rotation.random(samples, random_state=rng).as_matrix():
        rotated = points @ R.T
        B1 = bhat(rotated, weights)
        cross = cross_energy(points, weights, rotated, weights)
        for t in np.linspace(0.005, 0.995, 100):
            energy = 2 * t * (1 - t) * cross
            det = np.linalg.det((1 - t) * B0 + t * B1)
            gap = energy - (32 / 105) * det**2
            if gap < best[0]:
                best = (gap, (R, t, energy, det))
            if det > 1e-10 and energy / det**2 < best_ratio[0]:
                best_ratio = (energy / det**2, (R, t, energy, det))
    return best, best_ratio


if __name__ == "__main__":
    for name, base in [("pole_equator", pole_equator()), ("onb", onb())]:
        best, ratio = scan(base)
        print(name, "best gap", best)
        print(name, "best E/det^2", ratio, "target", 32 / 105)
    best, ratio = scan_random_rotations(onb())
    print("onb random best gap", best)
    print("onb random best E/det^2", ratio, "target", 32 / 105)
