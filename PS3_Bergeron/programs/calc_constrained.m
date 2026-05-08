function [V_grid_constrained, policy_grid_constrained] = calc_constrained(Params, z_grid, z_prob, a_next_grid, V_grid_old)

% [V_grid_constrained, policy_grid_constrained] = calc_constrained(...)
%
% Evaluates the Bellman equation for ALL current wealth levels a_i under
% the CONSTRAINED policy a' = a_next_grid(1) = a_1 (first exogenous grid point).
%
% Used for households where the unconstrained optimal a' < a_min, i.e.
% a_current_grid(ia,iz) <= a_min.
%
% Formula:
%   V(a_i, z_j) = u((1+r)*a_i + exp(z_j) - a_1) + beta*E[V(a_1,z')|z_j]
%
% Inputs:
%   a_next_grid - (n_a x 1) exogenous grid; a_1 = a_next_grid(1)
%   V_grid_old  - (n_a x n_z) current value function on exogenous grid

    a_1  = a_next_grid(1);   % BUG FIX 1: constrained choice is a_1, not a_min scalar
    R    = 1 + Params.r;

    % BUG FIX 2: E[V(a_1,z')|z] — continuation value evaluated at a_1 for each iz
    % V_grid_old(1,:) is V at a_1 for each income state (first row of exog grid)
    % z_prob is (n_z x n_z): EV_a1(iz) = sum_{iz'} pi(iz,iz') * V(a_1, iz')
    EV_a1 = z_prob * V_grid_old(1, :)';   % (n_z x 1)

    V_grid_constrained   = zeros(Params.n_a, Params.n_z);
    policy_grid_constrained = ones(Params.n_a, Params.n_z) * a_1;

    for iz = 1:Params.n_z
        z_curr = z_grid(iz);

        % consumption at each a_i when saving a_1: c_i = R*a_i + exp(z) - a_1
        % a_next_grid plays the role of a_current here (we evaluate at all grid pts)
        c = R * a_next_grid + exp(z_curr) - a_1;   % (n_a x 1)

        u      = (c .^ (1 - Params.gamma) - 1) / (1 - Params.gamma);
        u(c <= 0) = -Inf;

        % BUG FIX 3: EV_a1(iz) is a scalar — same for all ia since a'=a_1 fixed
        V_grid_constrained(:, iz) = u + Params.beta * EV_a1(iz);
    end
end