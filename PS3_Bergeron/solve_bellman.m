function [V_value, policy_value] = solve_bellman(Params, z_grid, z_prob, a, V_grid_old)

% [V] = solve_bellman(Params, z_grid, z_prob, a, V_grid_old)
%
% solves bellman equation for given value function grid (V^n) and wealth grid (a)
% returns updated value function grid (V^{n+1}) for each (a, iz) pair
%
% Params: struct of model parameters (e.g. beta, signa, etc.)
% Z_grid: discretized income grid
% z_prob: transition probabikities for income process
% a: current wealth grid
% V_grid_old: current value function (V^n) for each (a, iz) pair

V_value = zeros(Params.n_z); % initialize value function grid
policy_value = zeros(Params.n_z); % initialize policy grid

% consumption for each possible next-period wealth choice
for iz = 1:Params.n_z
    z = z_grid(iz); % current income
    c_vec = consumption(a, z, a_grid, Params); % (Na x 1)

    % utility 
    u = utility(c_vec, Params); % calculating current utility

    % continuation
    EV = zeros(Params.n_a, 1); % NA x 1

    % looping over all possible next income values 
    for iz_next = 1:Params.n_z
    % interpolate V(:, iz_next) at a_grid query points
    V_next = cubic_spline_interpolation(a_grid, V_grid_old(:, iz_next), a_grid);
    EV = EV + z_prob(iz, iz_next) * V_next; % (Na x 1)
    end

    total_val = u + Params.beta * EV; % NA x 1 
    total_val = total_val(:);          % force column vector to be safe
    
    % optimizing
             [V_value(iz), best_idx] = max(total_val);
             if isempty(best_idx) || best_idx < 1 || best_idx > Params.n_a
                best_idx = 1;  % fallback to minimum savings
             end
             policy_value(iz) = best_idx;
end

