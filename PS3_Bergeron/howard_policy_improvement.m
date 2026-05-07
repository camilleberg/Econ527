
function [V_grid, policy_grid, iteration] = howard_policy_improvement(Params, z_grid, z_prob, a_current_grid, max_policy_iter, V_grid, policy_grid)

% function [V_grid, policy_grid, iteration] = howard_policy_improvement(Params, z_grid, z_prob, a_current_grid, max_policy_iter, V_grid, policy_grid)
% 
% runs Howard policy improvement algorithm to solve for value and policy functions
% iteratively updates value function using policy function iteration with fixed policy for max_policy_iter iterations before
% updating policy function using VFI loop to get new policy and value grids
%
% V_grid: initial value grid (can be from previous VFI loop or initialized to zeros)
% z_grid: discretized income grid
% z_prob: transition probabilities for income process
% a_current_grid: current period wealth grid (a choices) for each (a, iz) pair
% policy_grid: initial policy grid (can be from previous VFI loop or initialized to zeros)
% max_policy_iter: maximum number of policy iterations to run before updating policy function within VFI loop
% iteration: number of iterations taken to converge 



    conv_error = Inf;
    iteration = 0;

    while conv_error > Params.e_stop * (1-Params.beta) && iteration < Params.max_iter

        V_grid_old = V_grid;

        % run regular VFI loop to get updated value and policy grids
        [V_grid, policy_grid] = VFI_calc(Params, z_grid, z_prob, a_current_grid, V_grid_old);
        iteration = iteration + 1; % updating vfi iteration count

        % checking error for convergence
        conv_error = max(abs(V_grid(:) - V_grid_old(:)));
        if mod(iteration, 50) == 0
            disp(['Iteration: ', num2str(iteration), ', Convergence error: ', num2str(conv_error)]);
        end
        if(conv_error <= Params.e_stop * (1-Params.beta))
            disp(['Convergence achieved after ', num2str(iteration), ' iterations with convergence error: ', num2str(conv_error)]);
            break;
        end

        policy_iter = 0;
        while policy_iter < max_policy_iter
            % update value grid using policy function iteration with fixed policy
            V_grid = policy_calc(Params, z_grid, z_prob, a_current_grid, V_grid, policy_grid);

            % adding 
            policy_iter = policy_iter + 1;
        end
    end
end
