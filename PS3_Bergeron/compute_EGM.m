function [a_current_grid] = compute_EGM(Params, z_grid, z_prob, a_next_grid, V_grid_old) 

% [a_current_grid] = compute_EGM(Params, z_grid, z_prob, a_next_grid, V_grid_old) 

% value function iteration loop but for endgenous grid (V^{n+1}) and policy grid (a') for each (a, iz) pair
% takes endogenous grid a' as input and calculates 
% a(a', iz) for each (a', iz) pair
%
% Params: struct of model parameters (e.g. beta, signa, etc.)
% Z_grid: discretized income grid
% z_prob: transition probabikities for income process
% a_prime_grid: current endogenous wealth grid (a' choices) for each (a, iz)
% a_next_grid: current exogenous wealth grid (a choices) for each (a', iz)

% solving for 
% a(a', iz) = (a'-exp(z) + u'^-1[\beta * E_z(V(a', z') | iz)])/(1+r)

% saving 
EV = V_grid_old * z_prob';

% initializing values
a_current_grid = zeros(Params.n_a, Params.n_z);

% unpacking parameters for convenience
beta = Params.beta;
R = 1 + Params.r;

 % looping through each (a, iz pair)
for ia = 1:Params.n_a
    a_next = a_next_grid(ia);
    for iz = 1:Params.n_z
        z_curr = z_grid(iz); % current income 

        % solving bellman equation (sort of) for given (a', iz) pair to get current a that maps to a' and iz

        % expected values for next period value function at a' and z' given current iz
        c_curr = max(utility_prime_inv(beta * EV(ia, iz), Params), 0);
        a_curr = (a_next - exp(z_curr) + c_curr) / R;

        a_current_grid(ia, iz) = a_curr;

    end
end

end