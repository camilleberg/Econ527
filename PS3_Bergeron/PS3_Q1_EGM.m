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


% creating polynomial wealth grid EGM
a_prime_grid = polynomial_grid(Params.a_min, Params.a_max, Params.n_a, Params.curve); % income grid using polynomial transformation

% using rouwenhorst method to discretize the AR(1) process for income
[z_grid, z_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n_z); % discretized income grid and transition probabilities using Rouwenhorst's method


%% Starting EGM 

% initializing values
V_grid = zeros(Params.n_a, Params.n_z); % initialize value function grid
a_grid_egm = zeros(Params.n_a, Params.n_z); % initialize policy function grid (a choices for each (a', iz) pair)
    % this will be filled in with updating loops 


%% 1. find a' 

 % VFI loop calculation to get V^{n+1} and a' grifd for each (a, iz) pair
vfi


%% 2. given a' at (a, iz) find optimal a that maps to a'

% backeards induction 
[a_grid_egm] = policy_fxn_egm_calc(Params, z_grid, z_prob, a_prime_grid, a_grid_egm) ;

% 3. sub into bellman equation 
% 4. use interpolation 


