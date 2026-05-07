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
a_grid = polynomial_grid(Params.a_min, Params.a_max, Params.n_a, Params.curve); % income grid using polynomial transformation

% using rouwenhorst method to discretize the AR(1) process for income
[z_grid, z_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n_z); % discretized income grid and transition probabilities using Rouwenhorst's method


%% Starting EGM 

% 1. find a' 

% initializing values
V_grid = zeros(Params.n_a, Params.n_z); % initialize value function grid
policy_grid = zeros(Params.n_a, Params.n_z); % initialize policy function grid
a_grid_egm = zeros(Params.n_a, Params.n_z); % initialize endogenous grid for a choices for each (a', iz) pair
    % this is the new a choices given k'prime 

V_grid_old = V_grid; % initialize old value grid for convergence check

 % VFI loop calculation to get V^{n+1} and a' grifd for each (a, iz) pair
[V_grid, policy_grid] = VFI_calc(Params, z_grid, z_prob, a_grid, V_grid_old, policy_grid);

% 2. given a' at (a, iz) find optimal a that maps to a'
[a_grid_egm] = policy_fxn_egm_calc(Params, z_grid, z_prob, policy_grid, a_grid_egm);

% 3. sub into bellman equation 

% check if new a_grid_egm is within bounds of original a_grid, if not, we need to extrapolate V and policy grid to get values for a_grid_egm outside of original a_grid bounds
%
% psedo code 
% if a_grid_edm(i) > a_min, them VFI calc with a_grid_egm as new grid and get new V_grid and policy_grid
% else, update V_grid(i) to be the budget constraint i.e. using a_grid 

% getting bounday bellman value
a_grid_min = Params.a_min * ones(1, Params.n_a); % next period's chocie have to be borrowing limit
[V_min, policy_min] = solve_bellman(Params, z_grid, z_prob, a_grid_min, Params.a_min, V_grid_old);

below = a_grid_egm < Params.a_min; % craetign index of values below borrowing limti
% updating V_grid and policy_grid for values below borrowing limit using budget constraint
V_grid(below) = V_min; % update value grid for values below borrowing limit using boundary bellman value
policy_grid(below) = policy_min; % update policy grid for values

% updating V_grid and policy_grid for values above borrowing limit using VFI loop with a_grid_egm as new grid
V_grid_old = V_grid; % update old value grid for convergence check
[V_grid, policy_grid] = VFI_calc(Params, z_grid, z_prob, a_grid_egm, V_grid_old, policy_grid); % run VFI loop with a_grid_egm as new grid to get updated value and policy grids for values above borrowing limit    


