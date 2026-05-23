% Econ 527 Spring 2026
% HW 4 Question 2 - Hugget Model
% Professor Greaney
% Written by: Camille Bergeron
% Due: May 22, 2026

clear; clc; close all;

%% Setting up Problem

% ── Inputs──────────────────────────────────────────────────────────────────

% 1. Defining the utility function 
Params.gamma = 2; % risk aversion
Params.rho = 0.9; % persistence of income process
Params.sigma = 0.2; % standard deviation of income shocks
Params.n_z = 7; % number of income grid points


Params.beta   = 0.96;         
Params.e_stop = 1e-6;
Params.max_iter = 10000;

%% Functions
% ── polynomial grid ──────────────────────────────────────────────────────────────────
polynomial_grid = @(a_min, a_max, n_a, curve) ...
    a_min + (a_max - a_min) * linspace(0, 1, n_a) .^ curve;

% ── Utilty function ──────────────────────────────────────────────────────────────────
u_fxn = @(c) c.^(1-Params.gamma)/(1-Params.gamma); 

% ── derivative of the inverse ──────────────────────────────────────────────────────────────────
u_prime_inv = @(mu) mu .^ (-1 / Params.gamma);

% using rouwenhorst method to discretize the AR(1) process for income
[z_grid, z_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n_z); % discretized income grid and transition probabilities using Rouwenhorst's method
disp('Transition probabilities:'); disp(z_prob);

% calculating the income states and distribution
income_grid = exp(z_grid); % income levels (used as z in the model)
earnings_dist = mc_invdist(z_prob); % income distribution                      
expected_earnings = income_grid' * earnings_dist; % expected earnings           


% setting grid parameters
Params.a_min = -expected_earnings; % the mean of income 
Params.a_max = expected_earnings * 50; % 50 times that 
Params.curve = 2; % curvature parameter
Params.n_a = 150; % number of welath grid points 
Params.n_af = 150;

% creating a_next_grid because EGM
a_grid = polynomial_grid(Params.a_min, Params.a_max, Params.n_a, Params.curve);
a_fine_grid = polynomial_grid(Params.a_min, Params.a_max, Params.n_af, Params.curve);

%% Caluclating Equilibrium interest Rate

% we can use the bisection algorithm 
r_low  = -0.05;
r_high = 1/Params.beta - 1 - 1e-4;


if ~exist('results', 'dir')
    mkdir('results');
end

% bisection algortithm 
if isfile('results/hugget_r_star.mat')
    disp('File exists and skipping equilbrium rate calculation...');
    disp('Loading file instead...');
    load('results/hugget_r_star.mat', 'r_star');
    fprintf('\nEquilibrium r* = %.6f\n', r_star);
else
    disp('File not found, calculating equilibrium interest rate distribution...')
    opts = optimset('Display','iter');
    r_star = fminbnd(@(r) f_nad(r, income_grid, z_prob, u_fxn, u_prime_inv, ...
                           Params, a_grid, a_fine_grid)^2, r_low, r_high, opts);
    fprintf('\nEquilibrium r* = %.6f\n', r_star);
    save('results/hugget_r_star.mat', 'r_star')
end



%% Finding the value function, policy function, and density

if isfile('results/hugget_model.mat')
    disp('File exists and skipping hugget model calculation...');
    disp('Loading file instead...');
    load('results/hugget_model.mat', 'V_grid', 'policy_fxn', 'phi_dist', 'net_assets');
    disp('Loaded hugget model');
else
    disp('File not found, calculating hugget model...')
    tic;
    [V_grid, policy_fxn, phi_dist, net_assets] = f_spe(r_star, income_grid, z_prob, u_fxn, u_prime_inv, ...
        Params, a_grid, a_fine_grid);
    time = toc;

    fprintf('Took %.4f seconds to converge\n', time)

    save('results/hugget_model.mat', 'V_grid', 'policy_fxn', 'phi_dist', 'net_assets');
end

%% Plotting the variables 

if ~exist('figs', 'dir')
    mkdir('figs');
end

% value function
% build legend labels for all income states
z_labels = arrayfun(@(i) sprintf('z_%d', i), 1:Params.n_z, 'UniformOutput', false);

% value function — one line per income state
fig_value = figure;
theme(fig_value, "light");
plot(a_fine_grid, V_grid', 'LineWidth', 2);
xlabel('Assets a');
ylabel('Value Function');
title('Value Function');
legend(z_labels, 'Location', 'best');
grid on;
saveas(gcf, 'figs/Value_Function_graph.png');

% policy function — one line per income state plus 45-degree line
fig_policy = figure;
theme(fig_policy, "light");
plot(a_fine_grid, policy_fxn', 'LineWidth', 2);   % (n_af x n_z), one line per z
hold on;
plot(a_fine_grid, a_fine_grid, '--k', 'LineWidth', 1.5); % 45-degree line
xlabel('Current Assets a');
ylabel('Next Period Assets a''');
title('Policy Function');
legend([z_labels, {'45° line'}], 'Location', 'best');
grid on;
saveas(gcf, 'figs/Policy_fxn_graph.png');

% stationary distribution — one line per income state
fig_stat = figure;
theme(fig_stat, "light");
plot(a_fine_grid, phi_dist', 'LineWidth', 2);      % (n_af x n_z), one line per z
xlabel('Assets a');
ylabel('Density');
title('Stationary Distribution');
legend(z_labels, 'Location', 'best');
grid on;
saveas(gcf, 'figs/Phi_Dist_graph.png');

%% Checking if wealth grid is binding 

% to check if welath grid is binding, we check if the policy function is the top of the wealth grid 
% check lower bound -- are any households choosing a' = a_min?
tol        = 1e-6;
lower_bind = sum(policy_fxn(:) <= Params.a_min + tol);
upper_bind = sum(policy_fxn(:) >= Params.a_max - tol);
disp(['States at lower bound (a_min=', num2str(Params.a_min), '): ', num2str(lower_bind)]);
disp(['States at upper bound (a_max=', num2str(Params.a_max), '): ', num2str(upper_bind)]);

%% plottign excess demand 
%  plot the excess demand function evaluated on an
% evenly spaced grid of 25 interest rates over the range [−0.02,0.02].

if isfile('results/hugget_model_nad.mat')
    disp('File exists and skipping nad calculations...');
    disp('Loading file instead...');
    load('results/hugget_model_nad.mat', 'net_assets_grid', 'r_grid');
    disp('Loaded hugget model');
else
    disp('File not found, calculating hugget net asset demands ')
    % plot the excess demand function over 25 interest rates in [-0.02, 0.02]
    tic; 
    r_grid = linspace(-0.02, 0.02, 25);
    net_assets_grid = zeros(1, 25);

    for ir = 1:length(r_grid)
        net_assets_grid(ir) = f_nad(r_grid(ir), income_grid, z_prob, u_fxn, u_prime_inv, ...
                                    Params, a_grid, a_fine_grid);
    end
    time = toc;

    fprintf('Took %.4f seconds to calculate\n', time)
    save('results/hugget_model_nad.mat', 'net_assets_grid', 'r_grid');
end


fig_nad = figure;
theme(fig_nad, "light")
plot(r_grid, net_assets_grid, 'b-o', 'LineWidth', 2);
hold on;
yline(0, '--k', 'LineWidth', 1.5);
xline(r_star, '--r', 'LineWidth', 1.5);
xlabel('Interest rate r');
ylabel('Net asset demand');
title('Excess Demand Function');
legend('Net asset demand', 'Zero line', 'r^*', 'Location', 'best');
grid on;
saveas(gcf, 'figs/huggett_excess_demand.png');

