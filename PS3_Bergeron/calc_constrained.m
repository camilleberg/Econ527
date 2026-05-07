function [V_grid_constrained] = calc_constrained(Params, z_grid, z_prob, a_current_grid, V_grid_old) 

% [value_fxn_val] = calc_constrained(Params, z_grid, z_prob, a_current_grid, V_grid_old, a_next_grid) 
%
% evaluates the bellman equation for the constrained case where a <= a_min for all possible iz
% to get value function value 
%
% Params: struct of model parameters (e.g. beta, signa, etc.)
% Z_grid: discretized income grid
% z_prob: transition probabilities for income process
% a_current_grid; current period wealth grid (a choices) for each (a, iz)
% V_grid_old: current value function (V^n) for each (a, iz) pair


V_grid_constrained = zeros(Params.n_a, Params.n_z);
EV_constrained = z_prob * V_grid_old(1,:)'; 
% since a' = a_min for all (a, iz) pairs, we only need to look at first 
% row of V_grid_old for expected value calculation since we are only 
% looking at the constrained case where a' = a_min for all (a, iz) pairs

for iz = 1:Params.n_z
    z_curr = z_grid(iz);

    c = consumption(a_current_grid, z_curr, Params.a_min, Params);
    c = c(:);

    % can't use utility fxm cause it said dot indexing wasn't allowed 
    u = (c.^(1-Params.gamma)-1)/(1-Params.gamma);
    u(c <= 0) = -Inf;

    V_grid_constrained(:,iz) = u + Params.beta * EV_constrained(iz);
end


end