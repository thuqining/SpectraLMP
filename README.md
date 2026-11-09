# SpectraLMP

Code, data pointers, parameters, and full results for the letter

> N. Qi, W. Wei, K. Hou, Z. Liu, H. Jia, D. W. Gao, and M.-Y. Chow, "Online Dispatch of Multi-Duration Storage via Spectral Decomposition of Locational Marginal Prices," submitted to *IEEE Power Engineering Letters*, 2026.

**SpectraLMP** values each storage unit by its bus locational marginal price (LMP), averaged over the past and future with the unit's time constant. Each unit then follows the high-pass part of its own LMP. Online, the method predicts these averages from past net load and solves a single-period network economic dispatch (ED). It needs no LMP forecasts, no value functions, and no multi-period optimization online.

- `paper/` holds the appendix of the letter (proofs and benchmark details, `SpectraLMP_proofs.pdf`).
- This README gives every setting that the letter leaves out for space.

---

## 1. Repository layout

| Folder | Content |
|---|---|
| `code/` | MATLAB code (R2024b) for all results and figures, and the Python XGBoost forecaster of the MPC benchmark |
| `data/` | place the OPSD household data file here (see Section 2) |
| `results/` | paired result tables (CSV) behind Fig. 2(b) and the sensitivity study |
| `paper/` | proofs of Theorem 1, Corollary 1, the noncausality of (5), expansion (6), Theorem 2, and Corollary 2, and the benchmark details |

**Requirements**
- MATLAB R2024b with the Statistics and Machine Learning Toolbox, and Gurobi 13.0 with its MATLAB interface. Set `GUROBI_HOME` if Gurobi is not in `C:\gurobi1302`.
- Python 3.11 with `numpy`, `scipy`, and `xgboost` 3.2.0, for the MPC forecast only.
- All reported tests, including the CPU times in Table I, ran on an Intel i9-13900HX CPU.

---

## 2. Data

- **Source:** the Open Power System Data (OPSD) household data package, 15-min resolution, file `household_data_15min.csv` (CoSSMic project), https://data.open-power-system-data.org/household_data/. Download it and save it as `data/household_data_15min.csv`.
- **Buildings and period:**
  - `letter_data.m` uses buildings industrial 1 and residential 5 and 6 over 425 days, from 2016-07-28 13:15 to 2017-09-26 13:00.
  - These are the CoSSMic buildings DE_KN_industrial1, DE_KN_residential5, and DE_KN_residential6 in Konstanz.
  - Community load is grid import + PV − export. PV is scaled by ×16.09, so that the mean net load is zero.
  - The net load of each building is its scaled net load minus its mean over the 425 days, so **the seasons stay in the net load** (`netrun.m` with `trend = 0`, the default of all sets ending in `14s`).
- **Bus mapping (IEEE 14-bus, MATPOWER case14, slack bus 1):**
  - buses 2–4 follow the total of the three buildings;
  - buses 5, 6, 9, 10 follow industrial 1;
  - buses 11–14 follow residential 6.
  - Each shape is scaled to the training deviation of the total, and each bus takes its IEEE load share of its shape (`net14_loads.m`).
- **Protocol:** the first 60 days train the policies and the next 365 days test them.

---

## 3. System and fleet parameters

**Generators (`net14_setup.m`)**

| Item | Value |
|---|---|
| Linear cost γ_j | 0.30 $/kWh at buses 1, 2; 0.60 $/kWh at buses 3, 6, 8 (case14 ratios) |
| Quadratic cost b_j | 0.005, 0.029, 0.0012, 0.0012, 0.0012 $/(kW²h) at buses 1, 2, 3, 6, 8 (case14 ratios, 0.005 at bus 1) |
| Curtailment generators | at the load buses, cost 0.002 s² $/h; they absorb negative net load |
| Line limits | only line 5–6 is limited in the congested case: 3.3 kW, the 65th percentile of its unconstrained flow in the hindsight optimum of the test year (`res.flowP`); it binds in 41.7% of the test steps |
| Generator bounds and ramps | not limited |

**Base storage fleet 4w (`params_fleet.m`, `synth_policy.m`)**

| Unit | Bus | P̄ (kW) | Duration Ē/P̄ (h) | η^c | η^d | c ($/(kW²h)) | c^τ ($/(kW²h)) | τ (h) |
|---|---|---|---|---|---|---|---|---|
| Battery | 3 | 44.8 | 1 | 0.96 | 0.96 | 1×10⁻³ | 1.78×10⁻³ | 0.89 |
| Battery | 4 | 33.6 | 4 | 0.95 | 0.95 | 3×10⁻³ | 3.98×10⁻³ | 2.9 |
| Flow battery | 13 | 25.2 | 30 | 0.85 | 0.85 | 10×10⁻³ | 12.8×10⁻³ | 81 |
| Hydrogen | 2 | 19.6 | 200 | 0.65 | 0.55 | 30×10⁻³ | 36.4×10⁻³ | 6035 |

- κ_i = c^τ_i/τ_i².
- Initial SoC e_{i,0} = 0.6 Ē_i.
- The other fleets in Fig. 2(b) have 2, 4, or 6 units with durations log-spaced over 4–30 h (fleets 2c, 4c, 6c) or 1–200 h (fleets 2w, 4w, 6w) (`netrun.m`, `netcfg.m`).

**Design and scoring parameters**

| Item | Value |
|---|---|
| Time step | 15 min, which resolves up to the Nyquist frequency of 48 cycles/day, above all crossovers |
| Risk level ρ | 0.32 |
| Predictor memory L | 768 steps (8 days), at most half of the available training data |
| q | 0.0073 $/(kW²h), so that 2qs matches the marginal cost γ₀ + 2q₊s of the producing generators at s = σ_x, with γ₀ = 0.30 and q₊ = 0.005 |
| Terminal value | every policy and the hindsight optimum are charged γ₀ η^d_i per kWh drawn from the initial stock (and credited per kWh added), so no policy gains by emptying its storage (`net5_cost.m`, `termval`) |
| Training hindsight QPs | final energy valued at γ₀ η^d_i |

---

## 4. Design steps of Algorithm 1

**Calibration (Step 1, `synth_policy.m`, `fir_design.m`)**
- κ_i and c^τ_i rise until a causal FIR policy reaches z σ_{p,i} ≤ P̄_i and z σ_{ε,i} ≤ min(e_{i,0}, Ē_i − e_{i,0}). The FIR policy uses L steps, has zero DC gain, and minimizes the discrete cost of the LQ model (3) with g = qs², so its taps solve one KKT system. Here z is the 1 − ρ/2 quantile of a zero-mean Gaussian model.
- These are per-step bounds, and the online QP enforces the hard limits.

**Targets (Step 2, `fit_split.m`)**
- Every month, the training hindsight ED (1) is solved on all past data. Its SoC duals μ_{i,t} are the targets.
- The past part of the average (4) is the exponential average λ̄_{i,t} of the hindsight LMPs at bus n(i) up to t − 1, with time constant τ_i.

**Predictor (Step 3, `price_features.m`, `fit_split.m`)**
- Ridge regression of μ_{i,t} on [1, λ̄_{i,t}, f_t], with ridge weight 0.1 n_s on n_s samples (no penalty on the intercept).
- Features f_t: the total net load averaged over each of the last 192 hours and the sine and cosine of the first two daily harmonics, standardized on the training data.

**Online (Steps 4–5, `method_values_net.m` with `opts.emaTau`)**
- λ̄_{i,t} is updated with the LMPs realized by the online dispatch, and μ̃_{i,t} = max(θ_iᵀ [1, λ̄_{i,t}, f_t], 0).
- A single-period DC-OPF QP values the stored energy of unit i at μ̃_{i,t}.

**Protocol (`netrun.m`)**
- The first 60 days (or 30 or 7 in the history study) are the history at commissioning, and the next 365 days are the test.
- The proposed method, ADP, and DQN are refitted every 30 days on all past data. The other policies are fitted once.
- All policies and the hindsight optimum are scored by (1) with the same e_{i,0} and the terminal value above.

---

## 5. Benchmarks

| Benchmark | Implementation |
|---|---|
| **MPC** (simplified from [3]) | Solves a 24-h DC-OPF QP every 30 min with the terminal penalty κ_i (e_{i,t+24 h} − e_{i,0})². Its forecast is learned from the training data only (`forecast_xgb.py`): the linear predictor d¹ x_{n,t} + d² x_{n,t+k−96} + d³ plus an XGBoost correction from recent lags, the values 1, 2, and 7 days earlier, and the time of day, kept for horizon k only if it lowers the error on the last 15 training days. Test RMSE: 24% and 53% of the net-load standard deviation at 15 min and 24 h ahead, against 24% and 57% for the linear predictor alone and 62% for day-ahead persistence (`forecast_lin_check.py`). |
| **ADP** (simplified from [4]) | Fits v_i(e) = max(v⁰_i + v¹_i e, 0), v¹_i ≤ 0, to μ_{i,t} of the same monthly hindsight QPs, and dispatches as in Algorithm 1. *ADP with features* (`fit_adpf.m`) adds ϑ_iᵀ f_t with the features of Step 3. |
| **DQN** (adapted from [5]) | One double deep Q-network per unit (`dqn_train.py`, PyTorch, 6-128-128-7): state = net load, its 1-h and 24-h exponential averages, SoC, time of day; 7 power levels in [−P̄_i, P̄_i]; shared reward = minus the system cost; 7-day episodes, 16 in parallel, 10⁵ steps. Seed and discount (0.99, 0.995, 0.999) are chosen on the last 15 training days (γ = 0.995); the agents are then retrained at every monthly refit on all past data. A single-period DC-OPF QP tracks the setpoints (`dqn_policy.m`). |
| **Cutoffs** (simplified from [6]) | Fixed low-pass cutoffs split the total net load into storage powers. They are tuned on the last 25 of the 60 training days. The single-period QP completes the dispatch with the first ADP values. |
| **Oracle** (noncausal) | Values e_{i,t} at α_i + β_i μ̄_{i,t}, with the exact averages (4) of the test-year hindsight LMPs and α_i, β_i fitted to the past μ_{i,t} at each refit |

Reference numbers follow the letter.

---

## 6. Additional results

**Diagnostics, base fleet (% of hindsight optimum, uncongested / congested; `seas_table.m`)**

| Variant | Cost |
|---|---|
| Proposed (Algorithm 1) | 109.0 / 109.0 |
| Without the past average λ̄ (same features, targets μ_{i,t}) | 109.2 / 108.9 |
| Earlier design: targets α_i + β_i μ̄_{i,t}, constant value for τ_i > L | 108.9 / 108.7 |
| ADP [4] (Table I) | 113.4 / 113.3 |
| ADP with the features f_t | 109.2 / 108.9 |
| DQN, retrained monthly (Table I) | 145.6 / 147.8 |
| DQN, trained once on 60 days | 153.6 / 156.5 |
| Linear-policy RL, cross-entropy search | 122.7 / 125.1 |
| No storage / myopic single-period ED (uncongested) | 154.1 / 114.4 |
| Oracle (exact future averages) | 103.7 / 103.7 |
| MPC, XGBoost forecast (Table I) | 116.5 / 116.5 |
| MPC, linear forecast alone | 117.6 / 117.5 |
| MPC, perfect 24-h forecasts | 102.9 / 103.1 |

Across the 36 settings below, ADP with the features f_t is within 1.0 point of the proposed method (mean −0.2), and the earlier design costs up to 3.7 points more with 7 days of history.

**Theorem 1 and Corollary 1 (`seas_numbers.m`)**
- The averages (4) explain 96.7%, 89.8%, and 91.8% of the variance of the hindsight μ_{i,t} of the 1-, 4-, and 30-h units without the line limit (96.7%, 89.9%, 94.5% with it).
- The storage shares of the proposed dispatch follow H_i/Σ H_i of (5), with crossovers at 6.9 h, 27 h, and 29 days.

**Fleets, congestion, and history (`results/paired_36.csv`)**
- The proposed method costs less than ADP in all 36 settings, by 7.7–8.4 points on average per history length; the gap ranges from 0.9 (fleet 2w) to 15.0 points.
- With 60 days of history, MPC with the XGBoost forecast costs 1.2–32.0 points more than the proposed method.
- The average cost of the proposed method is 108.5%, 109.3%, and 109.8% with 60, 30, and 7 days of history.
- The columns `adp_features`, `direct`, and `averaged_target` hold ADP with the net-load features f_t, the method without λ̄, and the earlier design.
- ADP with the features f_t is within 0.96 point of the proposed method in all 36 settings and lower in 25, so much of the gain over SoC-only ADP comes from these features.

**Theorem 2 (`verify_gap7.py`)**
- In the single-price LQ model (3) with a leak, the base-fleet weights, and a ridge predictor on hourly lags, the measured gap of the online dispatch matches (7), evaluated from the periodogram of the prediction errors, within 0.6%.

**CPU time (Table I, base fleet without congestion)**
- 1.2 ms per step for the proposed method, ADP, DQN, and the oracle, 1.1 ms for the cutoffs, and 5.8 ms for MPC (excluding the precomputed forecast).

---

## 7. Reproducing the letter

```matlab
cd code
for id = 1:2,  netrun(netcfg(id, 'main14s')); end    % Table I: ADP, cutoffs, oracle, MPC (linear, perfect)
for k = 1:2, c = netcfg(k, 'main14s'); c.methods = {'split', 'prop', 'hindcal'}; c.traj = (k == 1); ...
    c.tag = strrep(c.tag, 'main', 'splitmain'); netrun(c); end  % proposed (Table I) and its trajectories (Figs. 1, 2(c))
export_netload_s                                      % bus net loads for the MPC forecaster
```

```bash
python code/forecast_xgb.py                           # XGBoost-corrected MPC forecast (training data only)
python code/forecast_lin_check.py                     # out-of-sample forecast errors
```

```matlab
for id = 1:12, netrun(netcfg(id, 'mpcxgb14s')); end   % MPC with the XGBoost forecast (Table I, Fig. 2(b))
for id = 1:36, netrun(netcfg(id, 'split14s')); end    % proposed and the variant without the past average (Fig. 2(b))
for id = 1:36, netrun(netcfg(id, 'abl14s')); end      % ADP, ADP with features, earlier design (Fig. 2(b))
% DQN: write the monthly training sets, train each with dqn_train.py, then evaluate
c = netcfg(1, 'main14s'); c.methods = {'dqnm'}; c.traj = false; c.dqnEnvOnly = true; netrun(c)
% shell, per block XX = 01..13: DQN_CONFIGS=1:0.995 python code/dqn_train.py results/dqnm/env_XX.mat results/dqnm/w_XX.mat
% then rerun netrun with c.dqnEnvOnly = false for both line cases (Table I)
fig_idea_s; pol_s; seas_numbers; export_paired_s      % LMP cache, weights, Theorem 1 R^2, frequency split, paired tables
fig_frame1c; fig_results6                             % Figs. 1 and 2 (compile fig_framework.tex for Fig. 1)
seas_table                                            % Table I numbers
```

Each `netrun` call saves `results/net/<tag>.mat`. A base-case run takes about 10–40 minutes on the CPU above, depending on the methods it runs. Sets without the `s` suffix (`main14`, `fleet14`, …) reproduce an earlier version in which the net load excluded its past-week mean; the letter uses only the seasonal sets.

---

## License

Code: MIT. Data: see the OPSD household data package license.
