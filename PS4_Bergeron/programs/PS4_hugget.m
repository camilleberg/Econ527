% Econ 527 Spring 2026
% HW 4 Question 2 - Hugget Model
% Professor Greaney
% Written by: Camille Bergeron
% Due: May 22, 2026

clear; clc; close all;

%% Setting up Problem

% ── Utilty function ──────────────────────────────────────────────────────────────────
function u = u_fxn(c, Params)

    gamma = Params.gamma;
    u        = -Inf(size(c));       % default: -Inf for c <= 0
    pos      = c > 0;
    u(pos)   = (c(pos).^(1 - gamma)) / (1 - gamma);
end

function c = u_prime_inv(mu, Params)
% c = u_prime_inv(mu, Params)
    c = mu .^ (-1 / Params.gamma);
end

function u_prime = u_prime(c, Params) % (Na x 1)  
    u_prime = c.^(-Params.gamma); % inverse of CRRA utility
end


% ── Inputs──────────────────────────────────────────────────────────────────

% 1. Defining the utility function 
Params.gamma = 2; % risk aversion
Params.rho = 0.9; % persistence of income process
Params.sigma = 0.2; % standard deviation of income shocks
Params.n_z = 7; % number of income grid points


Params.beta   = 0.96;         
Params.e_stop = 1e-6;
Params.max_iter = 10000;

% using rouwenhorst method to discretize the AR(1) process for income
[z_grid, z_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n_z); % discretized income grid and transition probabilities using Rouwenhorst's method
disp('Transition probabilities:'); disp(z_prob);

% setting grid parameters
Params.a_min = -mean(exp(z_grid)); % the mean of income 
fprintf('The average income is %d ', Params.a_min);
Params.a_max = mean(exp(z_grid)) * 50; % 50 times that 
Params.curve = 2; % curvature parameter
Params.n_a = 150; % number of welath grid points 

% creating a_next_grid because EGM
a_next_grid = polynomial_grid(Params.a_min, Params.a_max, Params.n_a, Params.curve);
% and using the same for the fine grid
a_fine_grid = polynomial_grid(Params.a_min, Params.a_max, Params.n_a, Params.curve);

%% Caluclating Equilibrium interest Rate

% we can use the bisection algorithm 
r_low  = -0.05;
r_high = 1/Params.beta - 1 - 1e-4;

% bisection algortithm 
while (r_high - r_low) > 1e-4
    r_mid = (r_low + r_high) / 2; 
    Params.r = r_mid; 
    tic;
    net_assets = f_nad(z_grid, z_prob, @u_fxn, @u_prime_inv, Params, a_next_grid, a_fine_grid);
    time_nad = toc; 
    fprintf('Net assets are %.4f at real interest rate %.4f\n and took %.4f to run \n', net_assets, Params.r, time_nad);
    
    if net_assets > 0
        r_low = r_mid;    % excess saving → r too low, raise the floor
    else
        r_high = r_mid;   % excess borrowing → r too high, lower the ceiling
    end
end
fprintf('\nEquilibrium interest rate: r* = %.6f\n', (r_low + r_high)/2);


%% Finding the value function, policy function, and density

% tic;
% [v_fxn, policy_fxn, phi_dist, net_assets] = f_spe(r, markov_chain, u_fxn, u_prime_inv, Params, a_next_grid, a_fine_grid);
% time = toc;

% disp('Took ', num2str(time), " seconds")

%% Checking if wealth grid is binding 