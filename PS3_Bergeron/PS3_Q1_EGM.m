% Econ 527 Spring 2026 
% HW 3 Question 1
% Professor Greaney
% Written by: Camille Bergeron
% Due: May 08, 2026
% Last Edited: May 06, 2026


% Endogenous Grid MEthod

%% Setting up Problem 

clear; clc; close all; % clear workspaces 


% creating parameters, same as in PS2_Q2.m
Params.rho = 0.9; % persistence of income process
Params.sigma = 0.2; % standard deviation of income shocks
Params.beta = 0.96; % discount factor
Params.gamma = 1.5; % risk aversion parameter
Params.r = 0.02; % interest rate
Params.a_min = 0; % minimum income level
Params.a_max = 50; % maximum income level
Params.n_a = 100; % number of grid points for wealth grid
Params.curve = 2; % curvature parameter for polynomial grid
Params.n_z = 5; % number of points in discrete approximation for income process
Params.max_iter = 1000; % maximum number of iterations for value function iteration 
Params.e_stop = 1e-4; % convergence criterion for value function iteration


% creating polynomial wealth grid 
a_next_grid = polynomial_grid(Params.a_min, Params.a_max, Params.n_a, Params.curve); % income grid using polynomial transformation

% using rouwenhorst method to discretize the AR(1) process for income
[z_grid, z_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n_z); % discretized income grid and transition probabilities using Rouwenhorst's method



%% Initialzing values for EGM loop

V_grid = ones(Params.n_a, Params.n_z); % initialize value function grid for EGM loop
% changed initila guess from zeros because the inverse dooesn't liek that 
error = Inf; % initialize error for convergence check in EGM loop
iteration = 0; % initialize iteration count for EGM loop

%% Runing the loop 

tic; % start timer for EGM loop
while error > Params.e_stop * (1-Params.beta) && iteration < Params.max_iter

    V_grid_old = V_grid; % update old value function grid for convergence check in EGM loop

    % Calculating value function for a_min i.e. constrained
    [V_grid_constrained] = calc_constrained(Params, z_grid, z_prob, a_next_grid, V_grid_old); 
    % calculating value function for constrained case where a' = a_min for all (a, iz) pairs using calc_constrained function with V_grid initialized to zeros since we are only calculating the value function for the constrained case where a' = a_min for all (a, iz) pairs    

    % 1. discretize a' and find contnous a
    a_current_grid = compute_EGM(Params, z_grid, z_prob, a_next_grid, V_grid_old); 

    % 2. constarined versus unconstarined
    if any(a_current_grid(:) <= Params.a_min) % check if any are constrained
        disp('Constrained case detected, updating value and policy grids for constrained case');
        constrained_idx = a_current_grid <= Params.a_min; % get indices of constrained cases where a' <= a_min

        V_grid(constrained_idx, :) = V_grid_constrained(constrained_idx, :); % update value grid for values below borrowing limit using boundary bellman value from constrained case
        a_next_grid(constrained_idx, :) = Params.a_min; % update policy grid for values below borrowing limit using boundary bellman value from constrained case
    else
        disp('No constrained case detected, using EGM values for value and policy grids');

        unconstrained_idx = a_current_grid > Params.a_min; % get indices of unconstrained cases where a' >= a_min

        % calculating for unconstrained 
        [V_grid(unconstrained_idx, :), a_next_grid(unconstrained_idx, :)] = ...
            VFI_calc(Params, z_grid, z_prob, a_current_grid(unconstrained_idx), V_grid_old(unconstrained_idx, :)); % calculating value and policy grids for unconstrained cases using VFI_calc function with a_current_grid for unconstrained cases and V_grid initialized to zeros since we are only calculating the value and policy grids for the unconstrained cases where a' >= a_min for all (a, iz) pairs

    end

    error = max(abs(V_grid(:) - V_grid_old(:))); % calculate convergence error for EGM loop
    iteration = iteration + 1; % update iteration count for EGM loop

    if mod(iteration, 50) == 0
        disp(['Iteration: ', num2str(iteration), ', Convergence error: ', num2str(error)]);
    end

end

time = toc; % end timer for EGM loop
disp(['EGM loop completed after ', num2str(iteration), ' iterations with convergence error: ', num2str(error), ' and time: ', num2str(time), ' seconds.']);


