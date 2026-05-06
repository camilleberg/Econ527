function [a_grid_egm] = policy_fxn_egm_calc(Params, z_grid, z_prob, policy_grid, a_grid_egm) 

% [a_grid_egm] = policy_fxn_egm_calc(Params, z_grid, z_prob, policy_grid, a_grid_egm)  
%
% value function iteration loop but for endgenous grid (V^{n+1}) and policy grid (a') for each (a, iz) pair
% takes endogenous grid a' as input and calculates 
% a(a', iz) for each (a', iz) pair
%
% Params: struct of model parameters (e.g. beta, signa, etc.)
% Z_grid: discretized income grid
% z_prob: transition probabikities for income process
% a_prime_grid: current endogenous wealth grid (a' choices) for each (a, iz)
% a_grid_egm: current exogenous wealth grid (a choices) for each (a', iz)

% solving for 
% a(a', iz) = (a'-exp(z) + u'^-1[\beta * E_z(V(a', z') | iz)])/(1+r)

for ia_prime = 1:length(policy_grid(:)) % looping through each (a', iz) pair
    a_prime = policy_grid(ia_prime); % current poly grid
    for iz = 1:Params.n_z
        z= z_grid(iz); % current income 

       % calculate current consumptio
        c = consumption(a_prime, z, a_prime, Params); 
        % consumption for each (a', iz) pair using current endogenous grid a' and income z
        % i.e. next period wealth choice a' and current income z determine current consumption c

        % E_z(V(a', z') | iz)]
        up = utility_prime(c, Params); % marginal utility for each (a', iz) pair
        expected_up = up' * z_prob(iz, :)'; % expected marginal utility for each (a', iz) pair

        % u'^-1[\beta * E_z(V(a', z') | iz)]
        up_inv = max(utility_prime_inv(Params.beta * expected_up, Params), 1e-10); % inverse marginal utility for each 
        
        % (a'-exp(z) + u'^-1[\beta * E_z(V(a', z') | iz)])/(1+r)
        a_grid_egm(ia_prime, iz) = (a_prime - exp(z) + up_inv)/(1 + Params.r); % current wealth a for each (a', iz) pair using EGM formula
        % a is the current wealth level that would lead to the given a' choice and income
    end
end