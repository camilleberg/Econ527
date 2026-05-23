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
Params.n_af   = 150;
Params.curve  = 2;        % polynomial grid curvature

%% Functions

% ── polynomial grid ───────────────────────────────────────────────────────
polynomial_grid = @(a_min, a_max, n_a, curve) ...
    a_min + (a_max - a_min) * linspace(0, 1, n_a) .^ curve;

% ── utility function ──────────────────────────────────────────────────────
u_fxn = @(c) c.^(1-Params.gamma) / (1-Params.gamma);

% ── inverse marginal utility ──────────────────────────────────────────────
u_prime_inv = @(mu) mu .^ (-1 / Params.gamma);

%% ── Grids and income process ─────────────────────────────────────────────

% rouwenhorst for labor 
[l_grid, l_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n_z);

% income levels and stationary distribution
labor_grid = exp(l_grid);                        % (n_z x 1) labor level grid
labor_dist = mc_invdist(l_prob);                 % stationary distribution over z
expected_labor_supply = labor_grid' * labor_dist; % scalar E[l]

% asset grids
a_grid      = polynomial_grid(Params.a_min, Params.a_max, Params.n_a,  Params.curve);
a_fine_grid = polynomial_grid(Params.a_min, Params.a_max, Params.n_af, Params.curve);

%% ── Solve for GE equilibrium ─────────────────────────────────────────────

if ~exist('results', 'dir'); mkdir('results'); end

if isfile('results/aiyagari_model.mat')
    disp('File exists — loading instead of recomputing...');
    load('results/aiyagari_model.mat', 'V_grid', 'policy_fxn', 'phi_dist', ...
         'net_assets', 'r_star', 'w_star', 'k_star');
    disp('Loaded aiyagari model');
else
    disp('File not found — solving Aiyagari model...')

    tic;
    damp = 0.3;   % damping for capital updating
    k    = 5;     % initial guess for aggregate capital
    k_tol = 1e-5; % convergence tolerance for capital market clearing

    for iter = 1:Params.max_iter

        % factor prices from Cobb-Douglas firm FOCs
        r = Params.alpha     * k^(Params.alpha-1) * expected_labor_supply^(1-Params.alpha) - Params.delta;
        w = (1-Params.alpha) * k^Params.alpha      * expected_labor_supply^(-Params.alpha);

        % solve household problem: income = w * l
        [V_grid, policy_fxn, phi_dist, net_assets] = f_spe(r, labor_grid * w, ...
            l_prob, u_fxn, u_prime_inv, Params, a_grid, a_fine_grid);

        % check capital market clearing
        diff = abs(net_assets - k);
        fprintf('iter %4d | k = %.6f | net_assets = %.6f | diff = %.2e\n', ...
                iter, k, net_assets, diff);

        if diff < k_tol
            break
        end

        % damped update
        k = k + damp * (net_assets - k);
    end

    assert(iter < Params.max_iter, 'Aiyagari outer loop did not converge')

    % store equilibrium prices
    r_star = r;
    w_star = w;
    k_star = net_assets;

    time = toc;
    fprintf('\nConverged in %d iterations (%.2f seconds)\n', iter, time);
    fprintf('r* = %.6f,  w* = %.6f,  K* = %.6f\n', r_star, w_star, k_star);

    save('results/aiyagari_model.mat', 'V_grid', 'policy_fxn', 'phi_dist', ...
         'net_assets', 'r_star', 'w_star', 'k_star');
end

%% ── Wealth inequality ────────────────────────────────────────────────────
%
% phi_dist is (n_z x n_af): marginal wealth distribution is sum over z states

% marginal distribution over wealth (sum over income states)
wealth_dist = sum(phi_dist, 1);   % (1 x n_af), sums to 1

% ── Lorenz curve ─────────────────────────────────────────────────────────
% sort by wealth (a_fine_grid already increasing, so no re-sort needed)
cum_pop    = cumsum(wealth_dist);               % cumulative population share
cum_wealth = cumsum(a_fine_grid .* wealth_dist); % cumulative wealth
cum_wealth = cum_wealth / cum_wealth(end);      % normalise to [0,1]

% prepend (0,0) so curve starts at origin
cum_pop    = [0, cum_pop];
cum_wealth = [0, cum_wealth];

% ── Gini coefficient ─────────────────────────────────────────────────────
% Gini = 1 - 2 * area under Lorenz curve  (trapezoid rule)
gini = 1 - 2 * trapz(cum_pop, cum_wealth);

% ── Top 10% and top 1% wealth shares ─────────────────────────────────────
t10_share = 1 - interp1(cum_pop, cum_wealth, 0.90, 'linear', 'extrap');
t1_share  = 1 - interp1(cum_pop, cum_wealth, 0.99, 'linear', 'extrap');

fprintf('\n── Wealth Inequality ─────────────────────\n');
fprintf('Gini coefficient : %.4f\n', gini);
fprintf('Top 10%% share   : %.4f (%.1f%%)\n', t10_share, t10_share*100);
fprintf('Top  1%% share   : %.4f (%.1f%%)\n', t1_share,  t1_share*100);

%% ── Plots ────────────────────────────────────────────────────────────────

if ~exist('figs', 'dir'); mkdir('figs'); end

z_labels = arrayfun(@(i) sprintf('z_%d', i), 1:Params.n_z, 'UniformOutput', false);

% value function
v_fig = figure;
theme(v_fig, "light");
plot(a_fine_grid, V_grid', 'LineWidth', 2);
xlabel('Assets a'); ylabel('Value Function');
title('Value Function — Aiyagari');
legend(z_labels, 'Location', 'best'); grid on;
saveas(gcf, 'figs/aiyagari_value.png');

% policy function
policy_fig = figure;
theme(policy_fig, "light");
plot(a_fine_grid, policy_fxn', 'LineWidth', 2);
hold on;
plot(a_fine_grid, a_fine_grid, '--k', 'LineWidth', 1.5);
xlabel('Current Assets a'); ylabel('Next Period Assets a''');
title('Policy Function — Aiyagari');
legend([z_labels, {'45° line'}], 'Location', 'best'); grid on;
saveas(gcf, 'figs/aiyagari_policy.png');

% stationary distribution
phi_fig = figure;
theme(phi_fig, "light");
plot(a_fine_grid, phi_dist', 'LineWidth', 2);
xlabel('Assets a'); ylabel('Density');
title('Stationary Distribution — Aiyagari');
legend(z_labels, 'Location', 'best'); grid on;
saveas(gcf, 'figs/aiyagari_phi.png');

% Lorenz curve
lorenz_fig = figure;
theme(lorenz_fig, "light");
plot(cum_pop, cum_wealth, 'b-', 'LineWidth', 2); hold on;
plot([0 1], [0 1], '--k', 'LineWidth', 1.5);
xlabel('Cumulative Population Share');
ylabel('Cumulative Wealth Share');
title(sprintf('Lorenz Curve — Gini = %.4f', gini));
legend('Lorenz curve', '45° line (perfect equality)', 'Location', 'best');
grid on;
saveas(gcf, 'figs/aiyagari_lorenz.png');

%% ── Binding constraints check ────────────────────────────────────────────

bind_tol   = 1e-6;
lower_bind = sum(policy_fxn(:) <= Params.a_min + bind_tol);
upper_bind = sum(policy_fxn(:) >= Params.a_max - bind_tol);
fprintf('\nStates at lower bound (a_min = %.2f): %d\n', Params.a_min, lower_bind);
fprintf('States at upper bound (a_max = %.2f): %d\n', Params.a_max, upper_bind);

%% Calculating excess demand
if isfile('results/aiyagari_model_nad.mat')
    disp('File exists and skipping nad calculations...');
    disp('Loading file instead...');
    load('results/aiyagari_model_nad.mat', 'net_assets_grid', 'r_grid');
    disp('Loaded aiyagari nad');
else
    disp('File not found, calculating aiyagari net asset demands...')
    tic;
    r_grid = linspace(-0.02, 1/Params.beta - 1 - 1e-4, 25);
    net_assets_grid = zeros(1, length(r_grid));

    for ir = 1:length(r_grid)
        r_ir = r_grid(ir);
        % wage depends on r via firm FOC: w(r) uses k implied by r
        % At each r, use k_star as the capital stock to get the implied w
        % (i.e. we trace out household demand holding firm side fixed at GE k)
        k_ir = (Params.alpha / (r_ir + Params.delta))^(1/(1-Params.alpha)) ...
               * expected_labor_supply;
        w_ir = (1-Params.alpha) * k_ir^Params.alpha ...
               * expected_labor_supply^(-Params.alpha);

        net_assets_grid(ir) = f_nad(r_ir, labor_grid * w_ir, l_prob, u_fxn, u_prime_inv, ...
                                    Params, a_grid, a_fine_grid);
    end
    time = toc;
    fprintf('Took %.4f seconds to calculate\n', time)
    save('results/aiyagari_model_nad.mat', 'net_assets_grid', 'r_grid');
end

fig_nad = figure;
theme(fig_nad, "light")
plot(r_grid, net_assets_grid, 'b-o', 'LineWidth', 2);
hold on;
yline(0, '--k', 'LineWidth', 1.5);
xlabel('Interest rate r');
ylabel('Net asset demand (= K^d)');
title('Excess Demand Function — Aiyagari');
legend('Household asset demand', 'Market clearing', 'r^*', 'Location', 'best');
grid on;
saveas(gcf, 'figs/aiyagari_excess_demand.png');