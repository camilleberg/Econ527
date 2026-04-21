% Econ 527 Spring 2026 
% HW 2 Question 2
% Professor Greaney
% Written by: Camille Bergeron
% Due: April 30, 2026
% Last Edited: April 20, 2026


% Income Fluctuations

%% Setting up Problem 

% creating parameters
Params.rho = 0.9; % persistence of income process
Params.sigma = 0.2; % standard deviation of income shocks
Params.beta = 0.96; % discount factor
Params.gamma = 1.5; % risk aversion parameter
Params.r = 0.02; % interest rate
Params.a_min = 0; % minimum income level
Params.a_max = 50; % maximum income level
Params.a_nodes = 100; % number of grid points for wealth grid
Params.curve = 2; % curvature parameter for polynomial grid
Params.n = 5; % number of points in discrete approximation for income process
Params.max_iter = 1000; % maximum number of iterations for value function iteration 
Params.e_stop = 1e-4; % convergence criterion for value function iteration


% creating polynomial wealth grid 
a_grid = polynomial_grid(Params.a_min, Params.a_max, Params.a_nodes, Params.curve); % income grid using polynomial transformation
disp('Wealth grid points:');
disp(a_grid(1:10));

% using rouwenhorst method to discretize the AR(1) process for income
[z_grid, z_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n); % discretized income grid and transition probabilities using Rouwenhorst's method
disp('Transition probabilities:'); disp(z_prob);

% creating state space
state_space = combinations(a_grid, z_grid); % creating state space as the Cartesian product of wealth grid and income grid
disp('State space size:'); disp(size(state_space));
disp('First 10 states:'); disp(state_space(1:10, :));

%% Value Function Iteration
