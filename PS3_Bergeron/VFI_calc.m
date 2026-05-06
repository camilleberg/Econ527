function [V_grid, policy_grid] = vfi_calc(Params, z_grid, z_prob, a_grid, V_grid_old, policy_grid)

% [V_grid, policy_grid] = vfi_calc(Params, z_grid, z_prob, a_grid, V_grid_old, policy_grid)
%
% value function iteration loop to get update Value grid (V^{n+1}) and policy grid (a') for each (a, iz) pair
% updates both value function and policy function 
%
% Params: struct of model parameters (e.g. beta, signa, etc.)
% Z-grid: discretized income grid
% z_prob: transition probabikities for income process
%a_grid; current wealth grid
% C_grid_old: current value function (V^n) for each (a, iz) pair
% policy_grid: current policy grid (a' choices) for each (a, iz)


    % looping through each (a, iz pair)
    for ia = 1:Params.n_a
        a = a_grid(ia);
        for iz = 1:Params.n_z
            z = z_grid(iz); % curren income 

             % consumption for each possible next-period wealth choice
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
             [V_grid(ia, iz), best_idx] = max(total_val);
             if isempty(best_idx) || best_idx < 1 || best_idx > Params.n_a
                best_idx = 1;  % fallback to minimum savings
             end
             policy_grid(ia, iz) = best_idx;
        end
    end
end