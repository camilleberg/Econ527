% Econ 527 Spring 2026
% HW 4 Question 3 - Aiyagari Model
% Professor Greaney
% Written by: Camille Bergeron
% Due: May 22, 2026
%
% Solves the Aiyagari (1994) model:
%

clear; clc; close all;


%% ── Parameters ───────────────────────────────────────────────────────────

Params.gamma   = 2;       % CRRA risk aversion
Params.beta    = 0.96;    % annual discount factor
Params.rho     = 0.9;     % AR(1) persistence of log labor supply
Params.sigma   = 0.2;     % AR(1) innovation std dev
Params.n_z     = 7;       % income grid points
Params.alpha   = 0.36;    % capital share
Params.delta   = 0.10;    % depreciation rate
Params.e_stop  = 1e-6;    % VFI convergence tolerance
Params.max_iter = 2000;   % max VFI iterations

% Asset grid (borrowing limit = 0 per problem set)
Params.a_min  = 0;
Params.a_max  = 60;
Params.n_a    = 150;
Params.curve  = 2;        % polynomial grid curvature

%% ── Income process (Rouwenhorst) ─────────────────────────────────────────
[z_grid, z_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n_z);
z_grid = z_grid(:);   % (n_z x 1) column

% Stationary distribution of z (for computing L)
[V_mc, D_mc] = eig(z_prob');
[~, ii]      = min(abs(diag(D_mc) - 1));
pi_z_stat    = real(V_mc(:, ii));
pi_z_stat    = pi_z_stat / sum(pi_z_stat);   % (n_z x 1)

% Aggregate labour supply: L = E[exp(z)] under stationary distribution
L = sum(exp(z_grid) .* pi_z_stat);
fprintf('Aggregate labour supply L = %.6f\n', L);

%% ── Asset grids ──────────────────────────────────────────────────────────
a_grid      = polynomial_grid(Params.a_min, Params.a_max, Params.n_a, Params.curve);
a_fine_grid = a_grid;   % same grid for distribution (per problem set)

%% ── Firm factor prices given aggregate capital K ─────────────────────────
r_of_K = @(K) Params.alpha     * K^(Params.alpha - 1) * L^(1 - Params.alpha) - Params.delta;
w_of_K = @(K) (1 - Params.alpha) * K^Params.alpha     * L^(-Params.alpha);

%% ── Results directory ────────────────────────────────────────────────────
if ~exist('results', 'dir'); mkdir('results'); end
if ~exist('figs',    'dir'); mkdir('figs');    end

%% ── General Equilibrium: find K* such that K = E[a](K) ──────────────────
%
% At equilibrium: household asset demand equals aggregate capital supply.
% We bisect on K directly (equivalent to bisecting on r since r is
% monotone decreasing in K).
%
% Bounds:
%   K_low  → very small K → very high r → households want to save a lot
%             (net_assets > K_low, so excess demand > 0)
%   K_high → high K → r near zero or negative → little saving
%             (net_assets < K_high, so excess demand < 0)
%
% Excess demand: ED(K) = E[a](r(K), w(K)) - K
% Equilibrium:  ED(K*) = 0

% Reasonable bounds for K
% Steady state without uncertainty: beta*(1+r)=1 => r=1/beta-1-delta
r_max_nouncert = 1/Params.beta - 1;          % max r (no depreciation case)
r_min_bound    = -Params.delta + 0.001;      % r must exceed -delta for finite V

% Translate to K bounds
K_low  = (Params.alpha * L^(1-Params.alpha) / (r_max_nouncert + Params.delta))^(1/(1-Params.alpha)) * 0.5;
K_high = (Params.alpha * L^(1-Params.alpha) / (r_min_bound    + Params.delta))^(1/(1-Params.alpha)) * 1.5;
K_low  = max(K_low, 0.01);

fprintf('K bracket: [%.4f, %.4f]\n', K_low, K_high);
fprintf('r bracket: [%.4f, %.4f]\n', r_of_K(K_low), r_of_K(K_high));

% ── Bisection on K ────────────────────────────────────────────────────────
if isfile('results/aiyagari_Kstar.mat')
    fprintf('Loading cached K*...\n');
    load('results/aiyagari_Kstar.mat', 'K_star');
    fprintf('K* = %.6f  (r* = %.6f)\n', K_star, r_of_K(K_star));
else
    fprintf('\nStarting bisection on K...\n');
    tol_K = 1e-4;

    while (K_high - K_low) > tol_K
        K_mid = (K_low + K_high) / 2;

        % Set prices
        Params.r = r_of_K(K_mid);
        Params.w = w_of_K(K_mid);

        % Household income in this model: labour income = w*exp(z)
        % We pass z_grid as log-labour; f_spe uses exp(z_grid) as income.
        % Wage-adjusted income: household receives w*exp(z), not just exp(z).
        % We re-scale z_grid so that exp(z_scaled) = w*exp(z_original).
        % Equivalently define z_wage_grid = z_grid + log(w).
        z_wage_grid = z_grid + log(Params.w);   % log(w*exp(z)) = z + log(w)

        tic;
        net_assets = f_nad(z_wage_grid, z_prob, @u_fxn, @u_prime_inv, ...
                           Params, a_grid, a_fine_grid);
        t = toc;

        ED = net_assets - K_mid;   % excess demand for capital
        fprintf('K = %.4f | r = %.4f | w = %.4f | E[a] = %.4f | ED = %.4f | t=%.1fs\n', ...
                K_mid, Params.r, Params.w, net_assets, ED, t);

        if ED > 0
            K_low  = K_mid;   % too little capital → raise K
        else
            K_high = K_mid;   % too much capital → lower K
        end
    end

    K_star = (K_low + K_high) / 2;
    fprintf('\nEquilibrium K* = %.6f  r* = %.6f  w* = %.6f\n', ...
            K_star, r_of_K(K_star), w_of_K(K_star));
    save('results/aiyagari_Kstar.mat', 'K_star');
end

%% ── Solve model at equilibrium prices ────────────────────────────────────
Params.r     = r_of_K(K_star);
Params.w     = w_of_K(K_star);
z_wage_grid  = z_grid + log(Params.w);

fprintf('\nSolving full model at K*=%.4f, r*=%.4f, w*=%.4f...\n', ...
        K_star, Params.r, Params.w);

if isfile('results/aiyagari_model.mat')
    fprintf('Loading cached model...\n');
    load('results/aiyagari_model.mat', 'v_fxn','policy_fxn','phi_dist','net_assets');
else
    tic;
    [v_fxn, policy_fxn, phi_dist, net_assets] = f_spe( ...
        z_wage_grid, z_prob, @u_fxn, @u_prime_inv, ...
        Params, a_grid, a_fine_grid);
    fprintf('Solved in %.2f seconds\n', toc);
    fprintf('Market clearing check: E[a] = %.6f, K* = %.6f\n', net_assets, K_star);
    save('results/aiyagari_model.mat', 'v_fxn','policy_fxn','phi_dist','net_assets');
end

%% ── Equilibrium summary ──────────────────────────────────────────────────
r_star = Params.r;
w_star = Params.w;
Y_star = K_star^Params.alpha * L^(1-Params.alpha);
C_star = Y_star - Params.delta * K_star;   % aggregate consumption

fprintf('\n════════════════════════════════════\n');
fprintf('  AIYAGARI EQUILIBRIUM SUMMARY\n');
fprintf('════════════════════════════════════\n');
fprintf('  K*     = %.4f\n', K_star);
fprintf('  r*     = %.4f (%.2f%%)\n', r_star, r_star*100);
fprintf('  w*     = %.4f\n', w_star);
fprintf('  Y*     = %.4f\n', Y_star);
fprintf('  C*     = %.4f\n', C_star);
fprintf('  K*/Y*  = %.4f\n', K_star/Y_star);
fprintf('  E[a]   = %.4f (should = K* = %.4f)\n', net_assets, K_star);
fprintf('════════════════════════════════════\n\n');

%% ── Wealth distribution statistics ──────────────────────────────────────
% Marginal distribution over assets (summing over z states)
phi_marg = sum(phi_dist, 2);   % (n_fine x 1)

% ── Gini coefficient ──────────────────────────────────────────────────────
% Gini = 1 - 2 * integral of Lorenz curve
% Discrete formula: Gini = 1 - sum_{i=1}^{n} (F_i - F_{i-1})*(L_i + L_{i-1})
% where F_i = cumulative population share, L_i = cumulative wealth share

% Sort wealth grid (already sorted by construction)
% Compute cumulative population and wealth shares
cum_pop    = cumsum(phi_marg);                          % F_i
cum_wealth = cumsum(a_fine_grid .* phi_marg);           % numerator of Lorenz
total_wealth = sum(a_fine_grid .* phi_marg);            % = net_assets = K*
lorenz     = cum_wealth / total_wealth;                 % L_i

% Gini via trapezoidal integration of Lorenz curve
F = [0; cum_pop];
L = [0; lorenz];
Gini = 1 - 2 * trapz(F, L);

% ── Top 10% and top 1% wealth shares ─────────────────────────────────────
% Top 10%: households with wealth above the 90th percentile
pct_90 = find(cum_pop >= 0.90, 1, 'first');
pct_99 = find(cum_pop >= 0.99, 1, 'first');

top10_share = 1 - lorenz(pct_90);
top1_share  = 1 - lorenz(pct_99);

fprintf('WEALTH DISTRIBUTION STATISTICS\n');
fprintf('────────────────────────────────────\n');
fprintf('  Gini coefficient : %.4f\n', Gini);
fprintf('  Top 10%% share    : %.4f (%.1f%%)\n', top10_share, top10_share*100);
fprintf('  Top 1%% share     : %.4f (%.1f%%)\n', top1_share,  top1_share*100);
fprintf('────────────────────────────────────\n\n');

%% ── Check wealth grid binding ────────────────────────────────────────────
tol = 1e-6;
lower_bind = sum(policy_fxn(:) <= a_fine_grid(1)   + tol);
upper_bind = sum(policy_fxn(:) >= a_fine_grid(end) - tol);
total_states = Params.n_a * Params.n_z;

fprintf('GRID BINDING CHECK\n');
fprintf('  States at lower bound (a_min=%.2f): %d / %d (%.1f%%)\n', ...
        Params.a_min, lower_bind, total_states, 100*lower_bind/total_states);
fprintf('  States at upper bound (a_max=%.2f): %d / %d (%.1f%%)\n', ...
        Params.a_max, upper_bind, total_states, 100*upper_bind/total_states);
if upper_bind > 0
    warning('Upper wealth grid is binding — consider increasing a_max');
end

%% ── Plots ────────────────────────────────────────────────────────────────

% 1. Value function
fig1 = figure('Name','Value Function');
theme(fig1, "light");
plot(a_fine_grid, v_fxn(:,1),         'b-',  'LineWidth', 2, 'DisplayName', sprintf('z_1 = %.2f', z_grid(1)));
hold on;
plot(a_fine_grid, v_fxn(:,ceil(end/2)), 'g-',  'LineWidth', 2, 'DisplayName', sprintf('z_4 = %.2f', z_grid(ceil(Params.n_z/2))));
plot(a_fine_grid, v_fxn(:,end),        'r-',  'LineWidth', 2, 'DisplayName', sprintf('z_7 = %.2f', z_grid(end)));
xlabel('Assets a');  ylabel('V(a,z)');
title('Value Function — Aiyagari Model');
legend('Location','southeast');
grid on;
saveas(gcf, 'figs/Aiyagari_ValueFunction.png');

% 2. Policy function
fig2 = figure('Name','Policy Function');
theme(fig2, "light");
plot(a_fine_grid, policy_fxn(:,1),          'b-',  'LineWidth', 2, 'DisplayName', 'Low z (z_1)');
hold on;
plot(a_fine_grid, policy_fxn(:,ceil(end/2)), 'g-',  'LineWidth', 2, 'DisplayName', 'Med z (z_4)');
plot(a_fine_grid, policy_fxn(:,end),         'r-',  'LineWidth', 2, 'DisplayName', 'High z (z_7)');
plot(a_fine_grid, a_fine_grid,               'k--', 'LineWidth', 1, 'DisplayName', '45-degree line');
xlabel('Current Assets a');  ylabel('Next Period Assets a''');
title('Savings Policy Function — Aiyagari Model');
legend('Location','northwest');
grid on;
saveas(gcf, 'figs/Aiyagari_PolicyFunction.png');

% 3. Stationary wealth distribution (marginal)
fig3 = figure('Name','Wealth Distribution');
theme(fig3, "light");
% Plot as a smooth density using bar chart of marginal distribution
bar(a_fine_grid, phi_marg, 1, 'FaceColor',[0.2 0.4 0.7], 'EdgeColor','none');
xlabel('Wealth a');  ylabel('Density \phi(a)');
title(sprintf('Stationary Wealth Distribution — Aiyagari Model\n(Gini = %.3f, Top 10%% = %.1f%%, Top 1%% = %.1f%%)', ...
      Gini, top10_share*100, top1_share*100));
xlim([Params.a_min, min(Params.a_max, 30)]);   % zoom in on bulk of distribution
grid on;
saveas(gcf, 'figs/Aiyagari_WealthDist.png');

% 4. Lorenz curve
fig4 = figure('Name','Lorenz Curve');
theme(fig4, "light");
plot(F, L, 'b-', 'LineWidth', 2, 'DisplayName', sprintf('Lorenz (Gini=%.3f)', Gini));
hold on;
plot([0 1],[0 1], 'k--', 'LineWidth', 1, 'DisplayName', '45-degree (perfect equality)');
xlabel('Cumulative population share');
ylabel('Cumulative wealth share');
title('Lorenz Curve — Aiyagari Model');
legend('Location','northwest');
grid on;
saveas(gcf, 'figs/Aiyagari_Lorenz.png');

% 5. Joint distribution (heatmap)
fig5 = figure('Name','Joint Distribution');
theme(fig5, "light");
imagesc(1:Params.n_z, a_fine_grid, phi_dist);
set(gca, 'YDir','normal');
colorbar;
xlabel('Income state z_j');  ylabel('Assets a');
title('Joint Stationary Distribution \phi(a,z)');
xticks(1:Params.n_z);
xticklabels(arrayfun(@(x) sprintf('%.2f',x), z_grid, 'UniformOutput', false));
saveas(gcf, 'figs/Aiyagari_JointDist.png');

fprintf('\nAll figures saved to figs/\n');

%% ═══════════════════════════════════════════════════════════════════════
%  LOCAL FUNCTIONS
%% ═══════════════════════════════════════════════════════════════════════

function u = u_fxn(c, Params)
% CRRA utility: u(c) = c^(1-gamma) / (1-gamma)
    gamma    = Params.gamma;
    u        = -Inf(size(c));
    pos      = c > 0;
    u(pos)   = (c(pos).^(1 - gamma)) / (1 - gamma);
end

function c = u_prime_inv(mu, Params)
% Inverse marginal utility: (u')^{-1}(mu) = mu^{-1/gamma}
    c = mu .^ (-1 / Params.gamma);
end