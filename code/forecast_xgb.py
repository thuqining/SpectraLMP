"""XGBoost net-load forecast for the MPC benchmark (seasonal 14-bus study).

Reads results/netload_s.mat (export_netload_s.m), trains one XGBoost model
per net-load shape on the training steps only, and forecasts each test step
t0 (after x_{t0-1} is observed) for t0-1+h, h = 1..96. Bus net loads are
fixed multiples of three shapes, so the forecast of bus k is
Map[g, k] times that of its shape g. Saves results/forecast_xgb_s.mat
(Fg: T x 96 x G float32, Map: G x nb) and prints out-of-sample errors.
"""
import os, numpy as np, scipy.io as sio, xgboost as xgb

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = sio.loadmat(os.path.join(root, 'results', 'netload_s.mat'))
X = D['X'].astype(float); n0 = int(D['nTrain'][0, 0]); hrs = D['hrs'].ravel()
N, nb = X.shape; T = N - n0; H = 96; spd = 96

# shapes: columns that are multiples of each other
sd = X[:n0].std(0); Z = X / sd
grp = -np.ones(nb, int); reps = []
for k in range(nb):
    for g, r in enumerate(reps):
        if np.allclose(Z[:, k], Z[:, r], atol=1e-9): grp[k] = g; break
    else:
        reps.append(k); grp[k] = len(reps) - 1
G = len(reps); S = X[:, reps]
Map = np.zeros((G, nb))
for k in range(nb):
    r = reps[grp[k]]; Map[grp[k], k] = X[:, k] @ X[:, r] / (X[:, r] @ X[:, r])

dow = ((np.arange(N) * 0.25 / 24) + 3) % 7            # 2016-07-28 is a Thursday (Mon = 0)
def feats(x, t, h):
    tg = t + h
    m16 = np.array([x[max(0, i - 15):i + 1].mean() for i in range(len(x))])
    # level-free features: everything relative to the last observation x_t, so
    # the trees need not extrapolate when the seasonal level leaves the training range
    return np.column_stack([h, np.sin(2*np.pi*hrs[tg]/24), np.cos(2*np.pi*hrs[tg]/24), np.floor(dow[tg]),
                            x[t] - x[t-1], x[t] - x[t-2], x[t] - x[t-3], m16[t] - x[t],
                            x[tg - spd] - x[t], x[tg - 2*spd] - x[t], x[tg - 7*spd] - x[t],
                            x[t - spd] - x[t]])

Fg = np.zeros((T, H, G), np.float32); rep = {}
for g in range(G):
    x = S[:, g]
    t = np.arange(7*spd, n0 - 1, 2); tt, hh = np.meshgrid(t, np.arange(1, H + 1), indexing='ij')
    tt = tt.ravel(); hh = hh.ravel(); ok = tt + hh < n0; tt = tt[ok]; hh = hh[ok]
    # linear forecast of netrun 'mpcfit' (per horizon, training steps only); XGBoost
    # learns its residual, so it keeps the linear accuracy at short horizons
    C = np.zeros((H + 1, 3))
    for h in range(1, H + 1):
        t1 = np.arange(96 - h, n0 - h); t1 = t1[t1 >= 0]
        C[h] = np.linalg.lstsq(np.column_stack([x[t1], x[t1 + h - 96], np.ones(t1.size)]), x[t1 + h], rcond=None)[0]
    lin = lambda t, h: C[h, 0]*x[t] + C[h, 1]*x[t + h - 96] + C[h, 2]
    A = np.column_stack([feats(x, tt, hh), lin(tt, hh) - x[tt]]); y = x[tt + hh] - lin(tt, hh)
    prm = {'max_depth': 6, 'eta': 0.05, 'subsample': 0.8, 'colsample_bytree': 0.8,
           'tree_method': 'hist', 'seed': 1, 'nthread': 8, 'objective': 'reg:squarederror'}
    # per-horizon choice between linear and boosted forecasts on a holdout of the
    # last 15 training days (no test data), then refit on all training steps
    nv = n0 - 15*spd; fit = tt + hh < nv; val = tt >= nv
    bv = xgb.train(prm, xgb.DMatrix(A[fit], label=y[fit]), num_boost_round=400)
    rv = y[val] - bv.predict(xgb.DMatrix(A[val]))
    use = np.zeros(H + 1, bool)
    for h in range(1, H + 1):
        mh = hh[val] == h
        use[h] = np.mean(rv[mh]**2) < np.mean(y[val][mh]**2)
    bst = xgb.train(prm, xgb.DMatrix(A, label=y), num_boost_round=400)
    print('shape', g, 'boosted horizons (h, 15-min steps):', np.flatnonzero(use).tolist()[:5], '... count', use.sum())
    class _M:
        def predict(self, F): return bst.predict(xgb.DMatrix(F))
    mdl = _M()
    oo, h2 = np.meshgrid(n0 - 1 + np.arange(T), np.arange(1, H + 1), indexing='ij')
    oo = oo.ravel(); h2 = h2.ravel(); m = oo + h2 < N           # targets beyond the data are not needed
    P = np.zeros(oo.size); P[m] = lin(oo[m], h2[m]) + use[h2[m]] * mdl.predict(np.column_stack([feats(x, oo[m], h2[m]), lin(oo[m], h2[m]) - x[oo[m]]])); P[~m] = x[oo[~m]]
    Fg[:, :, g] = P.reshape(T, H).astype(np.float32)
    # out-of-sample error at 15 min and 24 h ahead, relative to the training std
    for hz in (1, 96):
        o = n0 - 1 + np.arange(T - hz + 1)
        e_x = Fg[:T - hz + 1, hz - 1, g] - x[o + hz]
        e_p = x[o + hz - spd] - x[o + hz]
        rep[(g, hz)] = (np.sqrt(np.mean(e_x**2)) / sd[reps[g]], np.sqrt(np.mean(e_p**2)) / sd[reps[g]])
cnt = np.bincount(grp, minlength=G)
for hz in (1, 96):
    xg = sum(cnt[g]*rep[(g, hz)][0] for g in range(G)) / nb
    pe = sum(cnt[g]*rep[(g, hz)][1] for g in range(G)) / nb
    print(f'h={hz:3d}: XGBoost RMSE/std {xg:.3f} | day-ahead persistence {pe:.3f}')
sio.savemat(os.path.join(root, 'results', 'forecast_xgb_s.mat'), {'Fg': Fg, 'Map': Map, 'grp': grp + 1})
print('groups', grp + 1, 'saved', Fg.shape)
