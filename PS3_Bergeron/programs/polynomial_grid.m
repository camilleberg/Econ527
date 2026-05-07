% creating grid points for capital using polynomial transformation
function grid_poly = polynomial_grid(point_a, point_b, nodes, curve)

% grid_poly = polynomial_grid(point_a, point_b, nodes, curve)
%
% function to create a polynomial grid
%
% point_a: minimum grid point
% point_b: maximum grid point
% nodes: no. of grid points
% curve: curvature parameter *(higher vals --> more curved)

    grid_normalized = linspace(0, 1, nodes); % normalized grid in [0,1]
    grid_poly = point_a + (point_b - point_a) * grid_normalized.^curve; % polynomial transformation
end