
function [V_grid] = policy_iter(Params, z_grid, z_prob, a_grid, policy_grid, max_policy_iter, V_grid)

% [V_grid] = policy_iter(Params, z_grid, z_prob, a_grid, policy_grid, max_policy_iter, V_grid)
%
% Policy iteration loop to get value grid given a policy grid and returns new value funcion (\Nu^{n+1)})
% Created to impleement Howard policy improvement method
% which iteratively evaluates a fixed policy until convergence before updating the policy grid
%
% Params: struct of model parameters (e.g. beta, signa, etc.)
% Z-grid: discretized income grid
% z_prob: transition probabikities for income process
% a_grid; current wealth grid
% policy_grid: current policy grid (a' choices) for each (a, iz)
% max_policy_iter: maximum number of iterations for policy evaluation loop
% V_grid: current value grid (V^n) for each (a, iz)

    error_policy = Inf;
    iteration_policy = 0;

    while error_policy > Params.e_stop * (1-Params.beta) && iteration_policy < max_policy_iter

        V_grid_old = V_grid;

        % FIX 3: policy EVALUATION — use fixed policy, no max()
        for ia = 1:Params.n_a
            a = a_grid(ia);
            for iz = 1:Params.n_z
                z = z_grid(iz);

                % look up the savings choice from the fixed policy
                a_next_idx = policy_grid(ia, iz);

                % scalar consumption at the chosen next-period asset level
                c_val = (1 + Params.r) * a + z - a_grid(a_next_idx);
                c_val = max(0, c_val);

                % scalar utility
                if c_val > 0
                    u_val = (c_val^(1 - Params.gamma)) / (1 - Params.gamma);
                else
                    u_val = -Inf;
                end

                % expected continuation value at chosen a_next
                EV = 0;
                for iz_next = 1:Params.n_z
                    V_next = cubic_spline_interpolation(a_grid, V_grid_old(:, iz_next), a_grid);
                    EV = EV + z_prob(iz, iz_next) * V_next(a_next_idx);
                end

                V_grid(ia, iz) = u_val + Params.beta * EV;
            end
        end

        error_policy = max(abs(V_grid(:) - V_grid_old(:))); % convergence on value function
        iteration_policy = iteration_policy + 1;
    end
end