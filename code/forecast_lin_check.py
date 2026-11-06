"""Out-of-sample error of the linear MPC forecast (as netrun 'mpcfit'):
x_hat(t+h) = d1_h x_t + d2_h x_{t+h-96} + d3_h, fitted per bus and h on the
training steps, compared with XGBoost (forecast_xgb_s.mat) and day-ahead
persistence, all relative to the training std of each bus."""
import os, numpy as np, scipy.io as sio
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = sio.loadmat(os.path.join(root, 'results', 'netload_s.mat'))
X = D['X'].astype(float); n0 = int(D['nTrain'][0, 0]); N, nb = X.shape; T = N - n0
F = sio.loadmat(os.path.join(root, 'results', 'forecast_xgb_s.mat')); Fg = F['Fg']; Map = F['Map']
sd = X[:n0].std(0)
for h in (1, 4, 16, 48, 96):
    el, ex, ep = [], [], []
    for k in range(nb):
        x = X[:, k]
        tt = np.arange(96 - h, n0 - h); tt = tt[tt >= 0]
        A = np.column_stack([x[tt], x[tt + h - 96], np.ones(tt.size)]); c = np.linalg.lstsq(A, x[tt + h], rcond=None)[0]
        o = n0 - 1 + np.arange(T - h + 1)                       # origin = last observed index
        yl = c[0]*x[o] + c[1]*x[o + h - 96] + c[2]
        yx = Fg[:T - h + 1, h - 1, :] @ Map[:, k]
        y = x[o + h]
        el.append(np.sqrt(np.mean((yl - y)**2))/sd[k]); ex.append(np.sqrt(np.mean((yx - y)**2))/sd[k])
        ep.append(np.sqrt(np.mean((x[o + h - 96] - y)**2))/sd[k])
    print(f'h={h:3d} ({h/4:5.2f} h): linear {np.mean(el):.3f} | XGBoost {np.mean(ex):.3f} | persistence {np.mean(ep):.3f}')
