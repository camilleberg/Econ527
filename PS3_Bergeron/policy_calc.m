function [V_grid] = policy_calc(Params, z_grid, z_prob, a_current_grid, V_grid_old, policy_grid)

% [V_grid_new, policy_grid] = policy_calc(Params, z_grid, z_prob, a_current_grid, V_grid, policy_grid)
%
% policy function iteration loop to get update Value grid (V^{n+1}) and policy grid (a') for each (a, iz) pair
% updates both value function and policy function 
%
% solves policy function iteration for each (a, iz) pair to get updated value and policy
% meant to be inside policy function iteration loop to get updated value and policy grids for each iteratio
%
%
%
% Params: struct of model parameters (e.g. beta, signa, etc.)
% Z_grid: discretized income grid
% z_prob: transition probabikities for income process
% a_current_grid; current period wealth grid (a choices) for each (a, iz)
% V_grid: current value function (V^n) for each (a, iz) pair (from previos vfi iteration)
% policy_grid: current policy grid (a' choices) for each (a, iz)

% saving 
EV = V_grid_old * z_prob';
V_grid = zeros(Params.n_a, Params.n_z);


% unpacking parameters for convenience
beta = Params.beta;
gamma = Params.gamma;
R = 1 + Params.r;


 % looping through each (a, iz pair)
for ia = 1:Params.n_a
    a_curr = a_current_grid(ia);
    for iz = 1:Params.n_z
        z_curr = z_grid(iz); % curren income 

        % solving bellman equation for given (a, iz) and fixed policy to get updated value and policy
        a_next_fixed = policy_grid(ia,iz);
        c_fixed      = R*a_curr + exp(z_curr) - a_next_fixed;
        V_grid(ia,iz) = ((c_fixed)^(1-gamma)-1)/(1-gamma) ...
              + beta * cubic_spline_interpolation(a_current_grid, EV(:,iz), a_next_fixed);

    end
end

end