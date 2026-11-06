"""Double deep Q-learning benchmark (as in [5]), PyTorch: one DDQN agent per
storage unit, discrete power levels, shared reward = minus the system cost.

usage: python dqn_train.py env.mat out.mat
env.mat (from netrun 'dqn'): x (training net load, kW), P, E, etac, etad, c0,
a_fuel, qs, qc, e0, mu, sx, dt.
State of agent i: [x/sx, m4/sx, m96/sx, e_i/E_i - 0.5, sin, cos] (the RL
features of the paper). Actions: p_i = level * P_i, levels in LEV.
Environment: the single-bus clipped dispatch of the RL benchmark (util_clip),
7-day episodes drawn from the training window, B episodes in parallel, and a
terminal charge a_fuel * etad_i (e0_i - e_T,i) for energy drawn from the
initial stock. Hyperparameters are chosen by the greedy cost on the training
window only. Weights are saved as W{l}_{i} (in x out) and b{l}_{i}.
"""
import sys, numpy as np, scipy.io as sio, torch, torch.nn as nn

torch.set_num_threads(4)
E_ = sio.loadmat(sys.argv[1], squeeze_me=True)
x = np.asarray(E_['x'], float); n = x.size
P, Emax = np.atleast_1d(E_['P']).astype(float), np.atleast_1d(E_['E']).astype(float)
etac, etad, c0 = [np.atleast_1d(E_[k]).astype(float) for k in ('etac', 'etad', 'c0')]
e0 = np.atleast_1d(E_['e0']).astype(float)
a, qs, qc, mu, sx, dt = [float(E_[k]) for k in ('a_fuel', 'qs', 'qc', 'mu', 'sx', 'dt')]
N = P.size
LEV = np.array([-1, -2/3, -1/3, 0, 1/3, 2/3, 1.0]); nA = LEV.size

xz = x - mu
m4 = np.zeros(n); m96 = np.zeros(n); a4 = b96 = 0.0
for t in range(n):
    a4 = 0.75*a4 + 0.25*xz[t]; b96 = (1 - 1/96)*b96 + xz[t]/96
    m4[t], m96[t] = a4, b96
tod = 2*np.pi*(np.arange(n) % 96)/96
EXO = np.stack([np.clip(xz/sx, -5, 5), np.clip(m4/sx, -5, 5), np.clip(m96/sx, -5, 5), np.sin(tod), np.cos(tod)], 1)
NV = 15*96; NT = n - NV                     # last 15 days: validation for model selection
base = np.mean((a*np.maximum(x, 0) + qs*np.maximum(x, 0)**2 + qc*np.minimum(x, 0)**2)*dt)

def states(t, e):                                  # t: (B,), e: (B,N) -> (B,N,6)
    ex = EXO[t]
    return np.concatenate([np.repeat(ex[:, None, :3], N, 1), (e/Emax - 0.5)[..., None],
                           np.repeat(ex[:, None, 3:], N, 1)], 2)

def step(t, e, act):                               # vectorized over B
    p = LEV[act]*P
    pmax = np.minimum(P, etad*e/dt); pmin = -np.minimum(P, (Emax - e)/(etac*dt))
    p = np.minimum(pmax, np.maximum(pmin, p))
    s = x[t] - p.sum(1)
    en = np.clip(np.where(p > 0, e - p/etad*dt, e - p*etac*dt), 0, Emax)
    cost = (a*np.maximum(s, 0) + qs*np.maximum(s, 0)**2 + qc*np.minimum(s, 0)**2 + (c0*p**2).sum(1))*dt
    return en, cost

def qnet():
    return nn.Sequential(nn.Linear(6, 128), nn.ReLU(), nn.Linear(128, 128), nn.ReLU(), nn.Linear(128, nA))

@torch.no_grad()
def greedy(nets, S):                               # S: (B,N,6) -> (B,N)
    St = torch.as_tensor(S, dtype=torch.float32)
    return np.stack([nets[i](St[:, i]).argmax(1).numpy() for i in range(N)], 1)

def greedy_cost(nets, T0=NT, T1=n):
    e = e0[None].copy(); J = 0.0
    for t in range(T0, T1):
        act = greedy(nets, states(np.array([t]), e))
        e, c = step(np.array([t]), e, act); J += c[0]
    return J + a*(etad*(e0 - e[0])).sum()

import os
STEPS = int(os.environ.get("DQN_STEPS", "100000"))
def run(seed, gamma, steps=STEPS, B=16, L=672):
    rng = np.random.default_rng(seed); torch.manual_seed(seed)
    on = [qnet() for _ in range(N)]; tg = [qnet() for _ in range(N)]
    for i in range(N): tg[i].load_state_dict(on[i].state_dict())
    opt = [torch.optim.Adam(on[i].parameters(), lr=5e-4) for i in range(N)]
    cap = 300000; Sb = np.zeros((cap, N, 6), np.float32); Ab = np.zeros((cap, N), np.int64)
    Rb = np.zeros(cap, np.float32); S2b = np.zeros((cap, N, 6), np.float32); Db = np.zeros(cap, np.float32)
    ptr = full = 0
    t = rng.integers(0, NT - L, B); k = np.zeros(B, int); e = np.repeat(e0[None], B, 0)
    for it in range(steps):
        S = states(t, e)
        eps = max(0.05, 1 - it/(0.6*steps))
        act = greedy(on, S)
        rnd = rng.random((B, N)) < eps
        act[rnd] = rng.integers(0, nA, rnd.sum())
        en, c = step(t, e, act); r = -c/base
        done = k == L - 1
        r = r - done*(a*(etad*(e0 - en)).sum(1))/base
        S2 = states(np.minimum(t + 1, n - 1), en)
        idx = (ptr + np.arange(B)) % cap
        Sb[idx], Ab[idx], Rb[idx], S2b[idx], Db[idx] = S, act, r, S2, done
        ptr = (ptr + B) % cap; full = min(full + B, cap)
        # advance; restart finished episodes
        t = t + 1; k = k + 1; e = en
        if done.any():
            m = done; t[m] = rng.integers(0, NT - L, m.sum()); k[m] = 0; e[m] = e0
        if full >= 5000:
            j = rng.integers(0, full, 256)
            s1 = torch.as_tensor(Sb[j]); s2 = torch.as_tensor(S2b[j])
            rr = torch.as_tensor(Rb[j]); dd = torch.as_tensor(Db[j]); aa = torch.as_tensor(Ab[j])
            for i in range(N):
                with torch.no_grad():
                    an = on[i](s2[:, i]).argmax(1, keepdim=True)                  # double Q
                    y = rr + gamma*(1 - dd)*tg[i](s2[:, i]).gather(1, an).squeeze(1)
                q = on[i](s1[:, i]).gather(1, aa[:, i:i+1]).squeeze(1)
                loss = nn.functional.smooth_l1_loss(q, y)
                opt[i].zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(on[i].parameters(), 10); opt[i].step()
        if it % 500 == 0:
            for i in range(N): tg[i].load_state_dict(on[i].state_dict())
    return on

best = None
CFG = os.environ.get('DQN_CONFIGS')                 # e.g. '1:0.995' fixes the configuration
SEEDS, GAMMAS = ((int(CFG.split(':')[0]),), (float(CFG.split(':')[1]),)) if CFG else ((1, 2), (0.99, 0.995, 0.999))
for seed in SEEDS:
    for gamma in GAMMAS:
        nets = run(seed, gamma)
        J = greedy_cost(nets)
        print(f"seed {seed} gamma {gamma}: validation cost {J:.1f} (no storage {base*NV:.1f})", flush=True)
        if best is None or J < best[0]: best = (J, nets, seed, gamma)
J, nets, seed, gamma = best
out = {'LEV': LEV, 'seed': seed, 'gamma': gamma, 'Jtrain': J}
for i in range(N):
    lin = [m for m in nets[i] if isinstance(m, nn.Linear)]
    for l, m in enumerate(lin):
        out[f'W{l+1}_{i+1}'] = m.weight.detach().numpy().T.astype(float)
        out[f'b{l+1}_{i+1}'] = m.bias.detach().numpy().astype(float)
sio.savemat(sys.argv[2], out)
print('saved', sys.argv[2], 'seed', seed, 'gamma', gamma)
