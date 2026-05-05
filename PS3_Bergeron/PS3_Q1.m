% Econ 527 Spring 2026 
% HW 3 Question 1
% Professor Greaney
% Written by: Camille Bergeron
% Due: May 08, 2026
% Last Edited: May 05, 2026


% Income Fluctuations

%% Setting up Problem 


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
disp('Wealth grid points:');
disp(a_grid(1:10));

% using rouwenhorst method to discretize the AR(1) process for income
[z_grid, z_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n_z); % discretized income grid and transition probabilities using Rouwenhorst's method
disp('Transition probabilities:'); disp(z_prob);

m_iter_bounds = {5, 10, 20}; % number of iterations

% nested looop]% outer loop -- VFI iteration and save policy fxn 
% inener loop - iterate (to convergence) on BELLMAN with fixed policy function from outer lop

%% setting up vfi iteeration for inner loop

% vfi iteration algorithm from PS2
function [K_grid, V_new, policy_grid_new] = vfi_iter(Params, V_grid, policy_grid, error, max_iter_M)
    % This function will run the vfi iteration for a fixed number of iterations
    % and return the updated value function grid and policy grid.

    m_iter = 0; % initilaize counter 
    while error > Params.e_stop * (1-Params.beta) && iteration < max_iter_M % continue until convergence or max iterations

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
                %[V_grid(ia, iz), policy_grid(ia, iz)] = max(total_val);

                [V_grid(ia, iz), best_idx] = max(total_val);
                if isempty(best_idx) || best_idx < 1 || best_idx > Params.n_a
                    best_idx = 1;  % fallback to minimum savings
                end
                policy_grid(ia, iz) = best_idx;
                K_grid(ia, iz) = a_grid(best_idx); % save the actual capital choice corresponding to the best index
            end
        end

        error = max(abs(V_grid - V_grid_old), [], 'all'); % compute maximum error between iterations
        m_iter = m_iter + 1; % increment iteration counter
    end
    V_new = V_grid; % reasingnconverged output to return value
    policy_grid_new = policy_grid; % reassign converged output to return policy grid
end

%% Setting up bellman for outer loop

function bellman_iter(Params, V_grid, policy_grid, max_iter_M)
    % This function will run the Bellman iteration for a fixed number of iterations
    % and return the updated value function grid and policy grid.
    
    % Initialize error and iteration counter
    error = Inf;
    n_iter = 0; % bellman iterations

    % do bellman iteration

    % check gor convergence and run inner loop for VFI iteration
    if error > Params.e_stop * (1-Params.beta)
        [K_grid, V_grid, policy_grid] = vfi_iter(Params, V_grid, policy_grid, error, max_iter_M); % run VFI iteration for inner loop
    end
    
end

%% setting up to run 
% initilaize value function and policy function


%% Iterating for n = 5, 10, 20 iterations for inner loop
for i = 1:length(m_iter_bounds)
    tic; % start times

    % intialize value function and policy function for each run
    V_grid = zeros(Params.n_a, Params.n_z); % initialize value function grid
    policy_grid = zeros(Params.n_a, Params.n_z); % initialize policy function grid
    m_iter = 0; % initialize iteration counter
    error = Inf; % initialize error for convergence check

    % setting max values
    max_iter_N = m_iter_bounds(i); % number of iterations for inner loop
    [K_grid, V_grid, policy_grid] = vfi_iter(Params, max_iter_N); % run VFI loop with specified number of iterations for inner loop
    time = toc;
    % save results
    save(['PS3_Q2_results_n', num2str(max_iter_N), '.mat'], 'K_grid', 'V_grid', 'policy_grid'); % save results to .mat file
    disp(['Completed Bellman iteration with ', num2str(max_iter_N), ' iterations for inner loop in ', num2str(time), ' seconds.']);

end

%% Local Functions

% making consumption function 
function c = consumption(a, z, a_grid, Params)
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

% cubic spline interpolation for the value function
function V_interp = cubic_spline_interpolation(K_grid, V_grid, K_query)
    % Replace -Inf with a large negative number for spline stability
    V_grid_clean = V_grid;
    V_grid_clean(~isfinite(V_grid_clean)) = -1e10;
    
    spline_interp = spline(K_grid, V_grid_clean);
    V_interp = ppval(spline_interp, K_query);
end

