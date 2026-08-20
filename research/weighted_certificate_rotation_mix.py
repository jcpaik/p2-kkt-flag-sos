"""Exact/numeric rotational-copy mixing constraints at a global minimizer."""

import numpy as np
from scipy.special import eval_legendre


COEFF={2:8/7,4:-384/385,6:512/231}


def intensities(x,w):
    g=x@x.T
    return {l:w@eval_legendre(l,g)@w for l in (2,4,6)}


def character_ratio(l,theta):
    return (1+2*sum(np.cos(m*theta) for m in range(1,l+1)))/(2*l+1)


def central_gap(I,theta):
    # I(mu, averaged rotated mu)-E(mu)
    return sum(COEFF[l]*I[l]*(character_ratio(l,theta)-1) for l in I)


def scan(x,w):
    I=intensities(x,w);th=np.linspace(0,np.pi,20001);gap=np.array([central_gap(I,t) for t in th]);i=gap.argmin()
    return I,th[i],gap[i]

