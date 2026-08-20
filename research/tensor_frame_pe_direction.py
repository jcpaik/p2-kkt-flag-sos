"""Directional transform data at a two-frame pole--equator zero measure."""

import numpy as np
import sympy as sp

from tensor_pair_circle_graph_exact import VARIABLES, frame_self_integrand


FRAME_SELF = sp.lambdify(VARIABLES, frame_self_integrand(), "numpy")


def directional_values(lam, phi, z, psi):
    """First variations toward one inserted atom.

    The base measure has a pole of weight 1/3 and two orthogonal equatorial
    pairs, of individual weights ``lam`` and ``1/3-lam``.  It is a mixture
    of two ONB measures sharing the pole, hence all target quantities vanish.
    The direction is ``delta_q-mu0``, where q has squared equatorial height
    z and azimuth psi.
    """
    pole = np.array([0.0, 0.0, 1.0])
    first = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    second = np.array([
        [np.cos(phi), np.sin(phi), 0.0],
        [-np.sin(phi), np.cos(phi), 0.0],
    ])
    inserted = np.array([
        np.sqrt(z)*np.cos(psi), np.sqrt(z)*np.sin(psi), np.sqrt(1-z),
    ])
    points = np.r_[pole[None], first, second, inserted[None]]
    weights = np.array([1/3, lam, lam, 1/3-lam, 1/3-lam, 0.0])
    direction = np.r_[-weights[:-1], 1.0]
    gram = points @ points.T

    kernel = 144*gram**6-216*gram**4+87*gram**2-5
    jprime = direction @ kernel @ weights + weights @ kernel @ direction

    a = gram[:, :, None]
    b = gram[:, None, :]
    c = gram[None, :, :]
    determinant = 1+2*a*b*c-a*a-b*b-c*c
    fkernel = (8/3)*determinant*(
        a*a*b*b*(c-a*b)**2
        + a*a*c*c*(b-a*c)**2
        + b*b*c*c*(a-b*c)**2
    )
    fprime = sum(
        np.einsum(
            "i,j,k,ijk", *(direction if slot == marked else weights
                             for slot in range(3)), fkernel
        )
        for marked in range(3)
    )

    aa = gram[:, :, None, None]
    u = gram[:, None, :, None]
    r = gram[:, None, None, :]
    v = gram[None, :, :, None]
    t = gram[None, :, None, :]
    U = u*u+v*v-2*aa*u*v
    V = r*r+t*t-2*aa*r*t
    L = u*r+v*t-aa*(u*t+v*r)
    UV = U*V
    akernel = aa**4*(
        32*L**6-48*L**4*UV+20*L**2*UV**2-2*UV**3
    )
    aprime = sum(
        np.einsum(
            "i,j,k,l,ijkl", *(direction if slot == marked else weights
                               for slot in range(4)), akernel
        )
        for marked in range(4)
    )
    tkernel = FRAME_SELF(aa, u, r, v, t, gram[None, None, :, :])
    tprime = sum(
        np.einsum(
            "i,j,k,l,ijkl", *(direction if slot == marked else weights
                               for slot in range(4)), tkernel
        )
        for marked in range(4)
    )
    return jprime-108*fprime-288*aprime, tprime, (jprime, fprime, aprime)


if __name__ == "__main__":
    from scipy.optimize import differential_evolution

    def objective(parameters):
        gap, self_energy, _ = directional_values(*parameters)
        return gap/self_energy if self_energy > 1e-10 else 1e3

    result = differential_evolution(
        objective,
        ((0.03, 0.30), (0.05, 1.52), (0.005, 0.95), (0, 1.57)),
        popsize=16, maxiter=100, seed=11, polish=True,
    )
    print(result.fun, result.x, directional_values(*result.x))
