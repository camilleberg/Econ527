function [V_grid, policy_grid] = VFI_calc(Params, z_grid, z_prob, a_current_grid, V_grid_old)

% [V_grid, policy_grid] = vfi_calc(Params, z_grid, z_prob, a_grid, V_grid_old, policy_grid)
%
% value function iteration loop to get update Value grid (V^{n+1}) and policy grid (a') for each (a, iz) pair
% updates both value function and policy function 
%
% solves bellman equation for each (a, iz) pair to get updated value and policy
% meant to be inside VFI loop to get updated value and policy grids for each iteration
%
% Params: struct of model parameters (e.g. beta, signa, etc.)
% Z_grid: discretized income grid
% z_prob: transition probabikities for income process
% a_current_grid; current period wealth grid (a choices) for each (a, iz)
% V_grid_old: current value function (V^n) for each (a, iz) pair
% policy_grid: current policy grid (a' choices) for each (a, iz)

% saving 
EV          = V_grid_old * z_prob';
V_grid      = zeros(Params.n_a, Params.n_z);
policy_grid = zeros(Params.n_a, Params.n_z);


% unpacking parameters for convenience
beta = Params.beta;
gamma = Params.gamma;
R = 1 + Params.r;
a_min = Params.a_min;
a_max = Params.a_max;

% looping through each (a, iz pair)
for ia = 1:length(a_current_grid)
    a_curr = a_current_grid(ia);
    for iz = 1:length(z_grid)
        z_curr = z_grid(iz); % curren income 

        % solving bellman equation for given (a, iz) pair to get updated value and policy
        obj = @(a_next) - ( (( R*a_curr + exp(z_curr) - a_next)^(1-gamma)-1)/(1-gamma) ...
            + beta * cubic_spline_interpolation(a_current_grid, EV(:,iz), a_next) );

        % Constraints: Cannot borrow past a_min, cannot consume more than wealth
        lower_bound = a_min;
        upper_bound = min(R*a_curr + exp(z_curr) - 1e-8, a_max);

        % Continuous optimization for a_prime
        [a_prime, fval] = fminbnd(obj, lower_bound, upper_bound);

        V_grid(ia,iz) = -fval;
        policy_grid(ia,iz) = a_prime;

    end
end

end