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
Params.n_a = 100; % number of grid points for wealth grid
Params.curve = 2; % curvature parameter for polynomial grid
Params.n_z = 5; % number of points in discrete approximation for income process
Params.max_iter = 1000; % maximum number of iterations for value function iteration 
Params.e_stop = 1e-4; % convergence criterion for value function iteration


% creating polynomial wealth grid 
a_grid = polynomial_grid(Params.a_min, Params.a_max, Params.n_a, Params.curve); % income grid using polynomial transformation
disp('Wealth grid points:');
disp(a_grid(1:10));

% using rouwenhorst method to discretize the AR(1) process for income
[z_grid, z_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n_z); % discretized income grid and transition probabilities using Rouwenhorst's method
disp('Transition probabilities:'); disp(z_prob);


%% Value Function Iteration

% initialize value function and policy function
V_grid = zeros(Params.n_a, Params.n_z); % initialize value function grid
policy_grid = zeros(Params.n_a, Params.n_z); % initialize policy function grid
iteration = 0; % initialize iteration counter
error = Inf; % initialize error for convergence check

while error > Params.e_stop * (1-Params.beta) && iteration < Params.max_iter % continue until convergence or max iterations

    V_grid_old = V_grid; % store old value function grid for convergence check

    for ia = 1:Params.n_a
        a = a_grid(ia);
        for iz = 1:Params.n_z
            z = z_grid(iz)

             % consumption for each possible next-period wealth choice
             c_vec = (1 + Params.r) * a + z - a_grid; % (Na x 1)

             % utility 
             u = utility(c, Params); % calcualting current utility

             % interpolating 
             V_next = cubic_spline_interpolation(a_grid, V_grid_old, a_grid); 

             % continuation
             EV = V_next * z_prob(iz, :); % calculating continuation
             total_val = u + Params.beta * EV; % calculating ottla value based on current utility and discounted 

             % optimizing
             [V_grid(ia, iz), policy_grid(ia, iz)] = max(total_val);
        end
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

% cubic spline interpolation for the value function
function V_interp = cubic_spline_interpolation(K_grid, V_grid, K_query)
    % take as input K grid and corresponding V grid, and query points K_query
    % K-query will be new guesses of K and will be updated
    % will return interpolated value function at the query points
    spline_interp = spline(K_grid, V_grid); % create spline interpolation object
    V_interp = ppval(spline_interp, K_query); % evaluate spline interpolation at query points
end 