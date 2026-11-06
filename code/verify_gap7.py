"""Numerical check of the gap identity (7) of Theorem 2 in its own setting.

Single system price, generator cost q s^2, the four units of the base fleet
(c^tau_i, kappa_i from results/pol_s.mat), a small energy leak nu, and a
zero-mean stationary Gaussian net load with seasonal-to-hourly components.
Hourly steps; the leaky SoC deviation is e_{t+1} = r e_t - p_t (r = exp(-nu)),
so the cost spectrum is a_i(w) = c_i + kappa_i / |1 - r e^{-jw}|^2 and all
formulas are the discrete-time versions of (3)-(7).

1. Hindsight optimum H x and its values mu (FFT over the whole sample).
2. Ridge predictor of mu from hourly lags of past net load (as Step 3),
   fitted on the first half, applied causally on the second half.
3. Online dispatch 2 c_i p_i = lambda - mu~_i with the balance, simulated
   step by step; measured gap = its mean cost minus that of H x.
4. Predicted gap (7): (1/2pi) int tr(Pi' M Pi Phi_delta) dw with Phi_delta
   estimated by the periodogram of delta = mu~ - mu.
"""
import os, numpy as np, scipy.io as sio
from scipy.signal import lfilter, csd

rng = np.random.default_rng(1)
pol = sio.loadmat(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "pol_s.mat"),  # from pol_s.m
                  squeeze_me=True, struct_as_record=False)['pol']
c = np.asarray(pol.ceff, float); kap = np.asarray(pol.kappa, float); N = c.size
q, nu = 0.01, 1e-3                      # generator weight, leak (1/h)
r = np.exp(-nu)
n = 2**19                               # hourly steps (~60 years)

def ar1(tau, sd):
    a = np.exp(-1/tau); return lfilter([sd*np.sqrt(1-a*a)], [1, -a], rng.standard_normal(n))
w0, rho = 2*np.pi/24, np.exp(-1/200)    # damped daily resonance
x = (ar1(2, 5) + ar1(48, 5) + ar1(700, 5) +
     lfilter([3.0], [1, -2*rho*np.cos(w0), rho**2], rng.standard_normal(n)))
x -= x.mean()

w = 2*np.pi*np.fft.fftfreq(n)            # rad/step
g = 1/np.abs(1 - r*np.exp(-1j*w))**2     # |E/P|^2
A = c[:, None] + kap[:, None]*g[None, :]  # a_i(w), N x n
den = 1/q + (1/A).sum(0)
X = np.fft.fft(x)
lam_f = 2*X/den                          # hindsight price spectrum
pH = np.real(np.fft.ifft(lam_f[None, :]/(2*A), axis=1))      # N x n
lamH = np.real(np.fft.ifft(lam_f))
mu = lamH[None, :] - 2*c[:, None]*pH                          # hindsight values

# causal ridge predictor on hourly lags (L = 192 h), fitted on the first half
L = 192
def lags(v):
    F = np.zeros((v.size, L))
    for k in range(L): F[k:, k] = v[:v.size-k]
    return F
h = n//2
F = lags(x); Z = np.hstack([np.ones((n, 1)), F])
tr = slice(L, h)
lamR = 0.1*(h-L); R = lamR*np.eye(L+1); R[0, 0] = 0
B = np.linalg.solve(Z[tr].T @ Z[tr] + R, Z[tr].T @ mu[:, tr].T)
mut = (Z @ B).T                                               # N x n

def cost(p, s):                                               # per-hour cost of a trajectory
    e = np.zeros_like(p)
    for i in range(N): e[i] = lfilter([-1.0], [1, -r], p[i])
    return e, (c[:, None]*p**2).sum(0) + (kap[:, None]*e**2).sum(0) + q*s**2

# online dispatch on the second half
D = 2*c
lam = (x + (mut/D[:, None]).sum(0))/((1/D).sum() + 1/(2*q))
p = (lam[None, :] - mut)/D[:, None]; s = lam/(2*q)
sH = lamH/(2*q)
burn = h + 20000                                              # leak and SoC transients
_, Jon = cost(p, s); _, JH = cost(pH, sH)
idx = slice(burn, n - 5000)                                   # drop the FFT wrap-around end
gap_meas = Jon[idx].mean() - JH[idx].mean()

# predicted gap (7)
delta = (mut - mu)[:, idx]
Pi = np.linalg.inv(np.diag(D)) @ (np.ones((N, N)) @ np.linalg.inv(np.diag(D)) /
                                 ((1/D).sum() + 1/(2*q)) - np.eye(N))
# periodogram of delta over the whole test segment (keeps the lowest
# frequencies, where the leaky weight kappa*|E/P|^2 is largest)
m = delta.shape[1]
Dk = np.fft.fft(delta, axis=1)                                 # N x m
wf = 2*np.pi*np.fft.fftfreq(m)
gf = 1/np.abs(1 - r*np.exp(-1j*wf))**2
Zk = Pi @ Dk                                                   # zeta = Pi delta
# (1/2pi) int tr(Pi' M Pi Phi_delta) dw ~ (1/m^2) sum_k zeta_k^* M_k zeta_k
quad = ((c[:, None] + kap[:, None]*gf[None, :])*np.abs(Zk)**2).sum(0) + q*np.abs(Zk.sum(0))**2
gap_pred = quad.sum()/m**2

print(f"hindsight cost J(H)        {JH[idx].mean():10.4f} per h")
print(f"measured gap  J - J(H)     {gap_meas:10.4f} per h ({100*gap_meas/JH[idx].mean():.1f}% of J(H))")
print(f"predicted gap (7)          {gap_pred:10.4f} per h")
print(f"relative difference        {100*(gap_pred-gap_meas)/gap_meas:+.2f}%")
print(f"corr(delta, x) max |rho|   {max(abs(np.corrcoef(delta[i], x[idx])[0,1]) for i in range(N)):.2f}")
