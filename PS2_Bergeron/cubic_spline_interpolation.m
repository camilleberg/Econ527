% cubic spline interpolation for the value function
function V_interp = cubic_spline_interpolation(K_grid, V_grid, K_query)
    % take as input K grid and corresponding V grid, and query points K_query
    % K-query will be new guesses of K and will be updated
    % will return interpolated value function at the query points
    spline_interp = spline(K_grid, V_grid); % create spline interpolation object
    V_interp = ppval(spline_interp, K_query); % evaluate spline interpolation at query points
end 