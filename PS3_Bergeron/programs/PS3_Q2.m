% Econ 527 Spring 2026 
% HW 3 Question 2
% Professor Greaney
% Written by: Camille Bergeron
% Due: May 08, 2026
% Last Edited: May 05, 2026


% Income Fluctuations and Stationary Distributon
clear; clc; close all;

%% Setting up Problem 


% creating parameters, same as in PS2_Q2.m
Params.rho = 0.9; % persistence of income process
Params.sigma = 0.2; % standard deviation of income shocks
Params.beta = 0.96; % discount factor
Params.gamma = 1.5; % risk aversion parameter
Params.r = 0.02; % interest rate
Params.a_min = 0; % minimum income level
Params.a_max = 50; % maximum income level
Params.n_a = 1000; % number of grid points for wealth grid
Params.curve = 2; % curvature parameter for polynomial grid
Params.n_z = 5; % number of points in discrete approximation for income process
Params.max_iter = 1000; % maximum number of iterations for value function iteration 
Params.e_stop = 1e-4; % convergence criterion for value function iteration


% fine space grid
a_grid_fine = linspace(Params.a_min, Params.a_max, Params.n_a);

% using rouwenhorst method to discretize the AR(1) process for income
[z_grid, z_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n_z); % discretized income grid and transition probabilities using Rouwenhorst's method
disp('Transition probabilities:'); disp(z_prob);

% importing V_grid and a_grid
load("./results/EGM_policy_grid.mat", 'policy_grid')
load("./results/EGM_valuefunction_grid.mat", 'V_grid')
load("results/EGM_a_next_grid.mat", "a_next_grid")

%% 1. compute policy 

% inialize
policy_fxn  = zeros(length(a_grid_fine), Params.n_z);

% interpolate for every income state 
for iz = 1:Params.n_z
    policy_fxn(:, iz) = cubic_spline_interpolation(a_next_grid, policy_grid(:, iz), a_grid_fine);
end

$$ 2. construct state transition matrix 
