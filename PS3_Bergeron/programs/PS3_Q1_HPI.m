% Econ 527 Spring 2026 
% HW 3 Question 1, pt1
% Professor Greaney
% Written by: Camille Bergeron
% Due: May 08, 2026
% Last Edited: May 06, 2026


% Howard Policy Improvement

% References
% https://github.com/fediskhakov/dcegm
% https://jax.quantecon.org/ifp_egm.html

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
a_current_grid = polynomial_grid(Params.a_min, Params.a_max, Params.n_a, Params.curve); % income grid using polynomial transformation
disp('Wealth grid points:');
disp(a_current_grid(1:10));

% using rouwenhorst method to discretize the AR(1) process for income
[z_grid, z_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n_z); % discretized income grid and transition probabilities using Rouwenhorst's method
disp('Transition probabilities:'); disp(z_prob);

m_iter_bounds = [5, 10, 20]; % number of iterations

%% setting up VFI iteration for outer loop 


% to create table of policy iterations and time 
if ~exist('results', 'dir')
    mkdir('results');
end
results = [transpose(m_iter_bounds), zeros(length(m_iter_bounds), 1), zeros(length(m_iter_bounds), 1)]; % initialize results table
% N x 1 and N x 1 time to make N x 2

for i = 1:length(m_iter_bounds)

    % initializing values
    V_grid = zeros(Params.n_a, Params.n_z); % initialize value function grid
    policy_grid = zeros(Params.n_a, Params.n_z); % initialize policy function grid

    max_policy_iter = m_iter_bounds(i); % get max policy iteration for this loop
    disp(['Running Howard policy improvement with max policy iterations = ', num2str(max_policy_iter), '...']);

    if isfile(['results/policy_grid_m_iter_', num2str(max_policy_iter), '.mat'])
        disp('File exists and skipping Howard policy improvement for this max policy iteration...');
    else
        disp('File does not exist. Running Howard policy improvement...');
    
        % running VFI
        tic; % start timer

        [V_grid, policy_grid, iteration] = howard_policy_improvement(Params, z_grid, z_prob, a_current_grid, max_policy_iter, V_grid, policy_grid); % run VFI loop and get value and policy grids

        % save policy grid for this loop
        save(['results/policy_grid_m_iter_', num2str(max_policy_iter), '.mat'], 'policy_grid');
        save(['results/value_grid_m_iter_', num2str(max_policy_iter), '.mat'], 'V_grid');
        
        time = toc; % saving time 

        % updating results table
        results(i, 2) = time; % store time in results table
        results(i, 3) = iteration; % store number of iterations
        writematrix(results, 'results/policy_improvement_results.csv');

        disp(['Completed Howard policy improvement with max policy iterations = ', num2str(max_policy_iter), ' in ', num2str(time), ' seconds']);
    end
end

if isfile('results/policy_improvement_results.csv')
    disp(readmatrix('results/policy_improvement_results.csv')); % display results table
else
    disp('Results file not found.');
end
 

