#!/usr/bin/env python3
"""Independent small-system AVBMC-1 and plain MC comparison."""
import argparse
import json
from math import comb, pi
from pathlib import Path
import numpy as np
from scipy.stats import chi2, ks_2samp

L = 6.0
RBIND = 1.5
VIN = 4*pi*RBIND**3/3
VOUT = L**3-VIN
N = 9
EPS = 1.0


def inside(pos, target):
    dr = pos-target
    dr -= L*np.rint(dr/L)
    return np.dot(dr, dr) < RBIND**2


def trial(rng, target, want_inside):
    if want_inside:
        v = rng.normal(size=3)
        v /= np.linalg.norm(v)
        return (target+RBIND*rng.random()**(1/3)*v) % L
    while True:
        p = rng.random(3)*L
        if not inside(p, target):
            return p


def quaternion(rng):
    u = rng.random(3)
    a, b = np.sqrt(1-u[0]), np.sqrt(u[0])
    x, y = 2*pi*u[1], 2*pi*u[2]
    return np.array([b*np.cos(y), a*np.sin(x), a*np.cos(x), b*np.sin(y)])


def particle_energy(pos, j, p, epsilon):
    if epsilon == 0:
        return 0.0
    dr = pos[np.arange(N) != j]-p
    dr -= L*np.rint(dr/L)
    d2 = np.sum(dr*dr, axis=1)
    if np.any(d2 < 0.7**2):
        return np.inf
    return -epsilon*np.count_nonzero(d2 < 1.25**2)


def total_energy(pos, epsilon):
    return sum(particle_energy(pos, j, pos[j], epsilon) for j in range(N))/2


def max_cluster(pos):
    parent = list(range(N))
    def root(x):
        while parent[x] != x:
            x = parent[x]
        return x
    for i in range(N):
        for j in range(i):
            d = pos[i]-pos[j]
            d -= L*np.rint(d/L)
            if np.dot(d,d) < 1.25**2:
                parent[root(i)] = root(j)
    return max(sum(root(j)==root(i) for j in range(N)) for i in range(N))


def avbmc_step(pos, quat, rng, epsilon, k=4, old=False, target_types=None, move_types=None, types=None, p_in=.5):
    # Fixed type pools; in this test all types are eligible.
    if types is None:
        types = np.zeros(N, dtype=int)
    targets = np.flatnonzero(np.isin(types, target_types)) if target_types is not None else np.arange(N)
    t = rng.choice(targets)
    in_move = rng.random() < p_in
    candidates = np.flatnonzero(np.isin(types, move_types)) if move_types is not None else np.arange(N)
    candidates = candidates[candidates != t]
    mask = np.array([inside(pos[j], pos[t]) for j in candidates])
    near, far = candidates[mask], candidates[~mask]
    ni, no = len(near), len(far)
    pool = far if in_move else near
    if len(pool) == 0:
        return
    j = rng.choice(pool)
    e_old = particle_energy(pos, j, pos[j], epsilon)
    old_energies = [e_old]
    for _ in range(k-1):
        p = trial(rng, pos[t], not in_move)
        quaternion(rng)
        old_energies.append(particle_energy(pos, j, p, epsilon))
    new_pos, new_quat, new_energies = [], [], []
    for _ in range(k):
        p = trial(rng, pos[t], in_move)
        q = quaternion(rng)
        new_pos.append(p); new_quat.append(q)
        new_energies.append(particle_energy(pos, j, p, epsilon))
    shift = min(old_energies+new_energies)
    if not np.isfinite(shift):
        return
    old_weights = np.exp(-np.subtract(old_energies, shift))
    weights = np.exp(-np.subtract(new_energies, shift))
    wnew, wold = sum(weights), sum(old_weights)
    if wnew == 0 or wold == 0:
        return
    chosen = rng.choice(k, p=weights/wnew)
    if old:
        # Legacy formulas with n_neigh=1, nbrd=N_in, nnonassoc=N_out.
        # This branch is exercised with K=1 for the ideal-gas diagnostic.
        assert k == 1
        factor = ((ni+1)/no)*(L**3/VIN) if in_move else ((no+1)/ni)*(VIN/L**3)
    else:
        factor = ((1-p_in)/p_in)*(VIN*no)/(VOUT*(ni+1)) if in_move else (p_in/(1-p_in))*(VOUT*ni)/(VIN*(no+1))
    if rng.random() < min(1., factor*wnew/wold):
        pos[j] = new_pos[chosen]
        quat[j] = new_quat[chosen]


def translation_step(pos, rng, epsilon):
    j = rng.integers(N)
    p = rng.random(3)*L
    old = particle_energy(pos,j,pos[j],epsilon)
    new = particle_energy(pos,j,p,epsilon)
    if np.isfinite(new) and rng.random() < min(1.,np.exp(-(new-old))):
        pos[j] = p


def run(seed, epsilon, method, steps=150000, thin=80):
    rng = np.random.default_rng(seed)
    pos = rng.random((N,3))*L
    quat = np.array([quaternion(rng) for _ in range(N)])
    # Equilibrate from the same diffuse start; core overlaps are removed by redraw.
    while not np.isfinite(total_energy(pos,epsilon)):
        pos = rng.random((N,3))*L
    counts, energies, clusters = [], [], []
    for step in range(steps):
        if method == 'plain':
            translation_step(pos,rng,epsilon)
        else:
            avbmc_step(pos,quat,rng,epsilon,k=1 if epsilon==0 else 4,old=method=='old')
            if epsilon == 0 and step%2 == 0:
                translation_step(pos,rng,epsilon)
        if step >= steps//4 and step%thin == 0:
            counts.append(sum(inside(pos[j],pos[0]) for j in range(1,N)))
            energies.append(total_energy(pos,epsilon))
            clusters.append(max_cluster(pos) if epsilon else 1)
    return np.array(counts),np.array(energies),np.array(clusters)


def categorical_p(a,b,lo,hi):
    x=np.bincount(a,minlength=hi+1)[lo:hi+1]
    y=np.bincount(b,minlength=hi+1)[lo:hi+1]
    valid=x+y>0
    obs=np.array([x[valid],y[valid]],dtype=float)
    expected=np.outer(obs.sum(axis=1),obs.sum(axis=0))/obs.sum()
    return float(chi2.sf(((obs-expected)**2/expected).sum(),valid.sum()-1))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--steps',type=int,default=150000)
    args=ap.parse_args()
    ideal={m:run(100+i,0,m,args.steps) for i,m in enumerate(('new','old'))}
    p=VIN/L**3
    expected=np.array([comb(N-1,k)*p**k*(1-p)**(N-1-k) for k in range(N)])
    ideal_stats={}
    for m,(count,_,_) in ideal.items():
        observed=np.bincount(count,minlength=N)
        observed=np.r_[observed[:3],observed[3:].sum()]
        exp_bins=np.r_[expected[:3],expected[3:].sum()]
        chi=np.sum((observed-len(count)*exp_bins)**2/(len(count)*exp_bins))
        ideal_stats[m]={'mean_Nin':float(count.mean()),'binomial_mean':(N-1)*p,'chi2_p':float(chi2.sf(chi,3))}
    plain=run(200,EPS,'plain',args.steps)
    biased=run(201,EPS,'new',args.steps)
    interacting={'energy_KS_p':float(ks_2samp(plain[1],biased[1]).pvalue),
                 'Nin_chi2_p':categorical_p(plain[0],biased[0],0,N-1),
                 'cluster_chi2_p':categorical_p(plain[2],biased[2],1,N),
                 'plain_mean_E':float(plain[1].mean()),'avbmc_mean_E':float(biased[1].mean())}
    out={'ideal':ideal_stats,'interacting':interacting,'steps':args.steps,'samples':len(plain[0])}
    path=Path(__file__).with_name('summary.json')
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':
    main()
