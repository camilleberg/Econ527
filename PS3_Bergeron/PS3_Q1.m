% Econ 527 Spring 2026 
% HW 3 Question 1
% Professor Greaney
% Written by: Camille Bergeron
% Due: May 08, 2026
% Last Edited: May 05, 2026


% Howar Policy

%% Setting up Problem 

clear; clc; close all; % clear wokrspaces 


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

%% setting up vfi iteeration for outer loop 

% iniitialixing values
V_grid = zeros(Params.n_a, Params.n_z); % initialize value function grid
policy_grid = zeros(Params.n_a, Params.n_z); % initialize policy function grid

for i = 1:length(m_iter_bounds)
    max_policy_iter = m_iter_bounds{i} % get max policy iteration for this loop

    % runnifn vfi 
    [V_grid, policy_grid] = howard_policy_improvement(Params, z_grid, z_prob, a_grid, max_policy_iter, V_grid, policy_grid); % run VFI loop and get value and policy grids

    % save policy grid for this loop
    save(['policy_grid_m_iter_', num2str(max_policy_iter), '.mat'], 'policy_grid');

    % save value grid for this loop
    save(['value_grid_m_iter_', num2str(max_policy_iter), '.mat'], 'V_grid');

    disp(['Completed Howard policy improvement with max policy iterations = ', num2str(max_policy_iter)]);
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


%% VFI loop function 
function [V_grid, policy_grid] = vfi_calc(Params, z_grid, z_prob, a_grid, V_grid_old, policy_grid); % run VFI loop and get value and policy grids
    policy_grid_old = policy_grid; % save old policy grid for convergence check
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
        end
    end
end

%% Policy iteration loop 

function [V_grid] = policy_iter(Params, z_grid, z_prob, a_grid, policy_grid, max_policy_iter, V_grid); % run policy iteration loop and get value grid
    
    % policy iteration loop
    while error_policy > Params.e_stop * (1-Params.beta) && iteration_policy < max_policy_iter % continue until convergence or max iterations
        policy_grid_old = policy_grid; % save old policy grid for convergence check

        % compute policy_fxn values at time t
        % but this is the current policy grid 

        % compute policy values at t+1
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
            end
        end 
        error_policy = max(abs(policy_grid(:) - policy_grid_old(:))); % check convergence of policy function
        iteration_policy = iteration_policy + 1; % increment iteration counter for policy iteration
    end
    disp(['Policy iteration converged in ', num2str(iteration_policy), ' iterations.']);
end

%%  function for howard policy improvement
function [V_grid, policy_grid] = howard_policy_improvement(Params, z_grid, z_prob, a_grid, max_policy_iter, V_grid, policy_grid); % run VFI loop and get value and policy grids
    error = Inf; % initialize error for outer loop convergence check
    iteration = 0; % initialize iteration counter for outer loop
    while error > Params.e_stop * (1-Params.beta) && iteration < Params.max_iter % continue until convergence or max iterations

        % runs reguar VFI loop to get value and policy grids
        V_grid_old = V_grid; % store old value function grid for convergence check
        [V_grid, policy_grid] = vfi_calc(Params, z_grid, z_prob, a_grid, V_grid_old, policy_grid); % run VFI loop and get value and policy grids

        % check convergence of value function
        error = max(abs(V_grid(:) - V_grid_old(:)));
        
        % check if converged, if not run policy iteration loop to update value function with fixed policy function
        if error < Params.e_stop * (1-Params.beta)
            disp('Converged in ', num2str(iteration))
            break; % break if converged
        else
            error_policy = Inf; % initialize error for policy iteration convergence check
            iteration_policy = 0; % initialize iteration counter for policy iteration
            [V_grid, policy_grid] = policy_iter(Params, z_grid, z_prob, a_grid, policy_grid, max_policy_iter, V_grid); % run policy iteration loop and get value grid
        end

        error = max(abs(policy_grid(:) - policy_grid_old(:))); % check convergence of policy function
        iteration = iteration + 1; % increment iteration counter for policy iteration
    end
end