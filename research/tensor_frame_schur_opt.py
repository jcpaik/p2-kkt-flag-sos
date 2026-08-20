"""Stress test the tied two-by-two J-Gram determinant.

The transform T is the weighted ONB transform in tensor_frame_transform_exact.md.
This is diagnostic only; exact identities/counterexamples are handled in
separate symbolic scripts.
"""

import argparse
import numpy as np
import torch

from tensor_pair_circle_gap_batch_opt import objective as moments

torch.set_default_dtype(torch.float64)


def run(n=7, batch=48, steps=6000, seed=734):
    rng = np.random.default_rng(seed+n)
    raw = torch.tensor(rng.normal(size=(batch,n,3)), requires_grad=True)
    logits = torch.tensor(rng.normal(size=(batch,n)), requires_grad=True)
    optimizer = torch.optim.Adam([raw,logits], lr=.02)
    for step in range(steps):
        optimizer.zero_grad()
        _, J, F, _, JTT, *_ = moments(raw, logits)
        cross = 108*F
        # The small positive denominator prevents the optimizer from seeing
        # only the trivial T=0 faces, without changing the sign being tested.
        determinant = J*JTT-cross*cross
        loss = determinant/(JTT+cross+1e-8)
        loss.sum().backward()
        optimizer.step()
        if step in (steps//2, 3*steps//4):
            for group in optimizer.param_groups:
                group["lr"] *= .2
    values = moments(raw, logits)
    _, J, F, A4, JTT, points, weights = values
    cross = 108*F
    determinant = J*JTT-cross*cross
    quotient = determinant/(JTT+cross+1e-8)
    index = int(torch.argmin(quotient))
    print("n,batch,steps,index",n,batch,steps,index)
    print("quotient,det,J,cross,JTT,A4", *[
        float(value[index].detach()) for value in
        (quotient,determinant,J,cross,JTT,A4)
    ])
    print("weights",weights[index].detach().numpy())
    print("points",points[index].detach().numpy())


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--n",type=int,default=7)
    parser.add_argument("--batch",type=int,default=48)
    parser.add_argument("--steps",type=int,default=6000)
    parser.add_argument("--seed",type=int,default=734)
    run(**vars(parser.parse_args()))
