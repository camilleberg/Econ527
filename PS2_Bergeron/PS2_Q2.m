% Econ 527 Spring 2026 
% HW 2 Question 2
% Professor Greaney
% Written by: Camille Bergeron
% Due: April 30, 2026
% Last Edited: April 20, 2026


% Income Fluctuations

%% Setting up Problem 

tic; 

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
            z = z_grid(iz);

             % consumption for each possible next-period wealth choice
             c_vec = consumption(a, z, a_grid, Params); % (Na x 1)

             % utility 
             u = utility(c_vec, Params); % calcualting current utility

             % continuation
             EV = zeros(Params.n_a, 1); % NA x 1

             % looping over all possible next income values 
            for iz_next = 1:Params.n_z
                % interpolate V(:, iz_next) at a_grid query points
                V_next = cubic_spline_interpolation(a_grid, V_grid_old(:, iz_next), a_grid);
                EV = EV + z_prob(iz, iz_next) * V_next; % (Na x 1)
            end

             total_val = u + Params.beta * EV; % NA x 1 
             total_val = total_val(:);              % force column vector to be safe

             % optimizing
             [V_grid(ia, iz), policy_grid(ia, iz)] = max(total_val);
        end
    end

    error = max(abs(V_grid - V_grid_old), [], 'all'); % compute maximum error between iterations
    iteration = iteration + 1; % increment iteration counter

    if mod(iteration, 50) == 0 
            disp(['Iteration: ', num2str(iteration), ', Error: ', num2str(error)]);
    end
end

time = toc;
disp(['Value function iteration converged in ', num2str(iteration), ' iterations and took ', num2str(time), ' seconds']);


%% Local Functions

% making consumption function 
function c = consumption(a, z, a_grid, Params)
    r = Params.r;
    c_possible= (1 + Params.r) * a + z - a_grid;
    c = max(0, c_possible);
end

% utilty function
function u = utility(c, Params)
    % Vectorized CRRA utility — handles (Na x 1) input
    gamma = Params.gamma;
    u        = -Inf(size(c));       % default: -Inf for c <= 0
    pos      = c > 0;
    u(pos)   = (c(pos).^(1 - gamma)) / (1 - gamma);
end
