function p = dqn_policy(Q, ex, e, par)
% DQN_POLICY  Greedy DDQN setpoints (dqn_train.py): ex = [x/sx, m4/sx, m96/sx,
% sin, cos] at step t, e = SoC (1 x N). Returns p (1 x N), positive = discharge.
N = par.N; p = zeros(1, N);
for i = 1:N
    s = [ex(1:3), e(i)/par.E(i) - 0.5, ex(4:5)];
    h = max(s*Q.W1{i} + Q.b1{i}, 0); h = max(h*Q.W2{i} + Q.b2{i}, 0);
    [~, k] = max(h*Q.W3{i} + Q.b3{i});
    p(i) = Q.LEV(k)*par.P(i);
end
end
