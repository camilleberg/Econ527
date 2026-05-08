function [V_grid_unconstrained] = calc_unconstrained(Params, a_current_grid, V_grid_old, a_next_grid, unconstrained_mask)

% [V_grid_unconstrained] = calc_unconstrained(...)
%
% Unconstrained EGM update: interpolates V from the endogenous grid
% a_current_grid back onto the fixed exogenous grid a_next_grid.
%
% For each income state iz, the endogenous grid points a_current_grid(:,iz)
% are the x-knots and V_grid_old(:,iz) are the corresponding values.
% We fit a spline and evaluate at the fixed a_next_grid query points.
%
% Inputs:
%   a_current_grid   - (n_a x n_z) endogenous current-wealth grid from EGM
%   V_grid_old       - (n_a x n_z) value function on exogenous grid
%   a_next_grid      - (n_a x 1)   fixed exogenous savings grid
%   unconstrained_mask - (n_a x n_z) logical, true where a_current > a_min

    V_grid_unconstrained = zeros(Params.n_a, Params.n_z);

    for iz = 1:Params.n_z
        unc_rows = unconstrained_mask(:, iz);

        a_endog = a_current_grid(unc_rows, iz);
        V_endog = V_grid_old(unc_rows, iz);

        [a_sorted, sort_idx] = sort(a_endog);
        V_sorted             = V_endog(sort_idx);

        [a_unique, unique_idx] = unique(a_sorted);
        V_unique               = V_sorted(unique_idx);

        a_unique = real(a_unique);
        V_unique = real(V_unique);

        % need at least 2 points to interpolate
        if length(a_unique) < 2
            V_grid_unconstrained(:, iz) = V_grid_old(:, iz);
            continue
        end

        % interpolate on full grid
        V_grid_unconstrained(:, iz) = cubic_spline_interpolation(a_unique, V_unique, a_next_grid);
    end
end