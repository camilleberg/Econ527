function [V_grid, policy_grid] = vfi_calc(Params, z_grid, z_prob, a_grid, V_grid_old, policy_grid)

    for ia = 1:Params.n_a
        a = a_grid(ia);
        for iz = 1:Params.n_z
            z = z_grid(iz);

            % consumption for each possible next-period wealth choice (n_a x 1)
            c_vec = consumption(a, z, a_grid, Params);

            % utility (n_a x 1)
            u = utility(c_vec, Params);

            % expected continuation value — (n_a x n_z) * (n_z x 1) = (n_a x 1)
            EV = V_grid_old * z_prob(iz, :)';

            % total value
            total_val = u + Params.beta * EV;
            total_val = total_val(:);

            % optimize
            [V_grid(ia, iz), best_idx] = max(total_val);
            if isempty(best_idx) || best_idx < 1 || best_idx > Params.n_a
                best_idx = 1;
            end
            policy_grid(ia, iz) = best_idx;
        end
    end
end