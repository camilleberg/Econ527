
function [V_grid, policy_grid, iteration] = howard_policy_improvement(Params, z_grid, z_prob, a_grid, max_policy_iter, V_grid, policy_grid)

    conv_error = Inf;
    iteration = 0;

    while conv_error > Params.e_stop * (1-Params.beta) && iteration < Params.max_iter

        V_grid_old = V_grid;

        % run regular VFI loop to get updated value and policy grids
        [V_grid, policy_grid] = vfi_calc(Params, z_grid, z_prob, a_grid, V_grid_old, policy_grid);

        % checking error for convergence
        conv_error = max(abs(V_grid(:) - V_grid_old(:)));

        if conv_error < Params.e_stop * (1-Params.beta)
            disp(['Converged in ', num2str(iteration), ' iterations.']);
            break;
        else
            % run policy iteration to update value function with fixed policy
            [V_grid] = policy_iter(Params, z_grid, z_prob, a_grid, policy_grid, max_policy_iter, V_grid);
        end

        iteration = iteration + 1;
    end
end