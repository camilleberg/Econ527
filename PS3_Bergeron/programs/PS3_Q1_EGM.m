% Econ 527 Spring 2026
% HW 3 Question 1 — Endogenous Grid Method
% Professor Greaney
% Written by: Camille Bergeron
% Due: May 08, 2026

clear; clc; close all;

%% Setting up Problem
Params.rho      = 0.9;
Params.sigma    = 0.2;
Params.beta     = 0.96;
Params.gamma    = 1.5;
Params.r        = 0.02;
Params.a_min    = 0;
Params.a_max    = 50;
Params.n_a      = 100;
Params.curve    = 2;
Params.n_z      = 5;
Params.max_iter = 10000;
Params.e_stop   = 1e-4;

a_next_grid      = polynomial_grid(Params.a_min, Params.a_max, Params.n_a, Params.curve);
[z_grid, z_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n_z);

R = 1 + Params.r;

%% ── Initialisation ───────────────────────────────────────────────────────────
% EGM tracks MARGINAL UTILITY mu = c^{-gamma} across iterations, not V.

% Use cash-on-hand consumption as starting c: c = R*a + exp(z) - a_1
% since zeros were throwing errors
a_1    = a_next_grid(1);
c_init = zeros(Params.n_a, Params.n_z);
for iz = 1:Params.n_z
    c_init(:, iz) = max(R * a_next_grid + exp(z_grid(iz)) - a_1, 1e-6);
end
mu_grid = c_init .^ (-Params.gamma);    % initial marginal utility

% Initialise V from cash-on-hand consumption discounted forever
V_grid = zeros(Params.n_a, Params.n_z);
for iz = 1:Params.n_z
    u0 = (c_init(:,iz).^(1-Params.gamma) - 1) / (1-Params.gamma);
    V_grid(:, iz) = u0 / (1 - Params.beta);
end

policy_grid = repmat(a_next_grid, 1, Params.n_z);   % initial policy: a' = a

error_val = Inf;
iteration = 0;

%% EGM loop
tic;
while error_val > Params.e_stop * (1 - Params.beta) && iteration < Params.max_iter

    V_grid_old  = V_grid;
    mu_grid_old = mu_grid;

    %── Step 1: EGM backward step ─────────────────────────────────────────────
    % compute_EGM uses Euler equation to recover c and endogenous a.
    [a_current_grid, c_current_grid] = compute_EGM(Params, z_grid, z_prob, a_next_grid, mu_grid_old);

    %── Step 2: constrained vs unconstrained ──────────────────────────────────
    constrained_mask   = a_current_grid <= Params.a_min;
    unconstrained_mask = ~constrained_mask;

    %── Step 3: value function update ─────────────────────────────────────────
    % Unconstrained: interpolate V from endogenous grid onto fixed grid
    V_grid_unconstrained = calc_unconstrained(Params, a_current_grid, V_grid_old, a_next_grid, unconstrained_mask);

    % Constrained: evaluate Bellman at forced a' = a_1
    [V_grid_constrained, ~] = calc_constrained(Params, z_grid, z_prob, a_next_grid, V_grid_old);

    % Combine: constrained entries override unconstrained interpolation
    V_grid = V_grid_unconstrained;
    V_grid(constrained_mask) = V_grid_constrained(constrained_mask);

    %── Step 4: update policy and marginal utility ────────────────────────────
    % Unconstrained: policy is the recovered endogenous a', consumption from EGM
    c_new = zeros(Params.n_a, Params.n_z);
    for iz = 1:Params.n_z
        unc = unconstrained_mask(:, iz);
        con = constrained_mask(:, iz);

        % unconstrained consumption: interpolate from endogenous grid
        a_endog = a_current_grid(unc, iz);
        c_endog = c_current_grid(unc, iz);
        [a_s, si] = sort(real(a_endog));
        c_s       = c_endog(si);
        [a_u, ui] = unique(a_s);
        c_u       = c_s(ui);
        if length(a_u) >= 2
            c_new(unc, iz) = max(cubic_spline_interpolation(a_u, c_u, a_next_grid(unc)), 1e-10);
        else
            c_new(unc, iz) = c_endog;
        end

        % constrained consumption: c = R*a + exp(z) - a_1
        c_new(con, iz) = max(R * a_next_grid(con) + exp(z_grid(iz)) - a_1, 1e-10);

        % policy
        policy_grid(unc, iz) = a_current_grid(unc, iz);
        policy_grid(con, iz) = a_1;
    end

    % update mu_grid from consumption (not from gradient of V)
    mu_grid = max(c_new, 1e-10) .^ (-Params.gamma);

    %── Step 5: convergence check ─────────────────────────────────────────────
    error_val = max(abs(V_grid(:) - V_grid_old(:)));
    iteration = iteration + 1;

    if mod(iteration, 100) == 0
        fprintf('Iteration %4d,  error = %.4e\n', iteration, error_val);
    end
end
time_egm = toc;

% to create table of policy iterations and time 
if ~exist('results', 'dir')
    mkdir('results');
end

save('results/EGM_policy_grid.mat', 'policy_grid');
save('results/EGM_valuefunction_grid.mat', 'V_grid');
save('results/EGM_a_next_grid.mat', 'a_next_grid');

fprintf('EGM converged after %d iterations in %.2f seconds (error = %.2e)\n', ...
        iteration, time_egm, error_val);



