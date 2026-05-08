% Econ 527 Spring 2026 
% HW 3 Question 2
% Professor Greaney
% Written by: Camille Bergeron
% Due: May 08, 2026
% Last Edited: May 05, 2026


% Income Fluctuations and Stationary Distributon
clear; clc; close all;

% References
% https://www.mathworks.com/help/econ/count-state-transitions.html
% https://brilliant.org/wiki/stationary-distributions/

% directory set up 

% to create table of policy iterations and time 
if ~exist('results', 'dir')
    mkdir('results');
end

% to create table of policy iterations and time 
if ~exist('figs', 'dir')
    mkdir('figs');
end


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

%%  2. construct state transition matrix 

% total states calc
total_states = Params.n_a * Params.n_z;
disp('Total states:'); disp(total_states);

% orderd 
%          (z1, a1) (z1, a2) ...
% (z1, a1)
% (z2, a2)

% calculating tranisionts
n_a_fine = Params.n_a;
numstates   = Params.n_a * Params.n_z;
policy_P    = zeros(numstates);       

for iz = 1:Params.n_z
    for ia = 1:n_a_fine 
        % definign starting point
        s_from  = (iz - 1) * n_a_fine  + ia;

        % defining next period by policy fxn 
        a_prime = policy_fxn(ia, iz);
        % getting a_prime as either the budget consr
        a_prime = max(a_grid_fine(1), min(a_grid_fine(end), a_prime));

        % finding index immediatley to the left of a_prime
        ik   = min(find(a_grid_fine <= a_prime, 1, 'last'), n_a_fine - 1);

        % interpolation weights 
        w_lo = (a_grid_fine(ik+1) - a_prime)  / (a_grid_fine(ik+1) - a_grid_fine(ik));
        w_hi = (a_prime - a_grid_fine(ik))     / (a_grid_fine(ik+1) - a_grid_fine(ik));

        % going through income states 
        for iz_next = 1:Params.n_z
            pi_z  = z_prob(iz, iz_next);
            s_lo  = (iz_next - 1) * n_a_fine  + ik;
            s_hi  = (iz_next - 1) * n_a_fine  + (ik + 1);
            policy_P(s_lo, s_from) = policy_P(s_lo, s_from) + pi_z * w_lo;
            policy_P(s_hi, s_from) = policy_P(s_hi, s_from) + pi_z * w_hi;
        end
    end
end

%% Calculate stationary dist using eigenvalue method


if isfile('results/stationary_eigen.mat')
    disp('File exists and skipping stationary calculation...');
    disp('Loading file instead...');
    load('results/stationary_eigen.mat', 'stationary_eigen');
else
    disp('File not found, calculating stationary distribution...')
    tic;

    % computing dist
    stationary_eigen = mc_invdist(policy_P);

    time_eigen = toc;

    disp(['Computed stationary distribution usign eigenvalue method in ', num2str(time_eigen), ' seconds']);
    save('results/stationary_eigen.mat', 'stationary_eigen');
end

% Plotting
lambda_grid = reshape(stationary_eigen, Params.n_a, Params.n_z);
lambda_marginal_a = sum(lambda_grid, 2);   % (n_a x 1) marginal over z

fig = figure;
theme(fig, "light");
plot(a_grid_fine, lambda_marginal_a, 'LineWidth', 2);
xlabel('Wealth (a)');
ylabel('Stationary distribution φ(a, z)');
title('Stationary Distribution using Eigen value method');
grid on;    
saveas(gcf, 'figs/stationary_eigen.png');

% net demand for assets (should be 0)
net_demand = sum(lambda_marginal_a - a_grid_fine, 'all');
disp("Net demand for assets is "); disp(net_demand)


%% Using simulation 

Params.n_duration = 2000;
Params.n_sim = 1000;
Params.n_sim_burn = 1000;


% simualting markov chain
if isfile('results/markov_sim_income.mat')
    disp('File exists and skipping simulation...');
    disp('Loading file instead...');
    load('results/markov_sim_income.mat', 'markov_sim_income', 'time_simulate');
else
    disp('File not found, generating simulation...')
    tic;

    % simualting markov, getting invar from mc_indivst i.e eigen
    markov_sim_income = simulate_markov(z_prob, Params.n_duration, Params.n_sim);
    time_simulate = toc;

    disp(['Simulated ', num2str(Params.n_sim ) ' histories with ', num2str(Params.n_duration) ' periods in ', num2str(time_simulate), ' seconds']);
    save('results/markov_sim_income.mat', 'markov_sim_income', 'time_simulate');
end

disp('income draws')
disp(markov_sim_income(1:10))

tic;
% intiialing wealth array
sim_a = simulate_a(policy_fxn, a_grid_fine, markov_sim_income, Params);
time_asset_sim = toc;
disp(sim_a(1:10, 1:10))

% dropping first 1000 sim
% This creates a new version of A without the first 1000 rows
sim_a_reduced = sim_a(Params.n_sim_burn:end, :);

total_time = time_asset_sim + time_simulate;

% computign policy function 
disp('Total time for simulation'); disp(total_time);

% asset demand 
net_asset_demand_sim = mean(sim_a_reduced, 'all');
disp("Net demand for assets is "); disp(net_asset_demand_sim)


