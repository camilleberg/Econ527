% cubic spline interpolation for the value function
function V_interp = cubic_spline_interpolation(K_grid, V_grid, K_query)

% [V_interp = cubic_spline_interpolation(K_grid, V_grid, K_query)
%
%  Cubic spline interpolation between grid points (for value function)
%
% K_grid: existing / known grid points (e.g. asset grid)
% V_grid: value function evaluated at K_grid points / or function to be eval at
% K_query: new grid points where we want to evaluate the value function 
%
% reference: https://www.mathworks.com/help/matlab/ref/spline.html

    V_grid_clean = V_grid;
    V_grid_clean(~isfinite(V_grid_clean)) = -1e10;
    
    spline_interp = spline(K_grid, V_grid_clean);
    V_interp = ppval(spline_interp, K_query);
end