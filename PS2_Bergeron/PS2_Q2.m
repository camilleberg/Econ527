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

% initialize value function and policy function
V_grid = zeros(size(state_space, 1), 1); % initialize value function grid
policy_grid = zeros(size(state_space, 1), 1); % initialize policy function grid
iteration = 0; % initialize iteration counter
error = Inf; % initialize error for convergence check

while error > Params.e_stop * (1-Params.beta) && iteration < Params.max_iter % continue until convergence or max iterations

    V_grid_old = V_grid; % store old value function grid for convergence check
    V_next = cubic_spline_interpolation(a_grid, V_grid_old, a_grid); % interpolate value function for next period's wealth grid

    for i = 1:size(state_space, 1) % loop over all states
        a = state_space(i, 1); % current wealth
        z = state_space(i, 2); % current income

        % compute expected values
        expected_values = zeros(size(a_grid)); % initialize expected values for each action
        c = consumption(a{1,1}, z{1,1}, a_grid, Params); % compute consumption for each possible next period's wealth
        u = utility(c, Params); % compute utility from consumption for each possible next period's wealth

        % spline interpolation for next period's value function
        

        % compute expected value of next period's value function
        EV_next = 0; % initialize expected value of next period's value function
        for k = 1:length(z_grid) % loop over possible next period's income states
            z_next = z_grid(k); % next period's income
            P = z_prob(j, k); % transition probability from current to next income state
            EV_next = EV_next + P * V_next(state_space == [a_next, z_next]); % sum over all possible next period's states
        end
        
        expected_values = u + Params.beta * EV_next; % total expected value for current action

        [V_grid(i), idx] = max(expected_values); % update value function and policy function for current state
        policy_grid(i) = a_grid(idx); % store optimal next period's wealth as policy function

    end

    error = max(abs(V_grid - V_grid_old)); % compute maximum error between iterations
    iteration = iteration + 1; % increment iteration counter

end



%% Local Functions

% making consumption function 
function c = consumption(a, z, a_next_grid, Params)
    r = Params.r;
    c_possible = (1 + r) * a + z - a_next_grid(:); % broadcasts to (N_a x N_z)
    c = max(0, c_possible);
end

% utilty function
function u = utility(c, Params)
    gamma = Params.gamma; % risk aversion parameter
    if c > 0
        u = (c^(1-gamma)) / (1-gamma); % CRRA utility function
    else    
        u = -Inf; % assign negative infinity utility for non-positive consumption
    end 
end 