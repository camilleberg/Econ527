% Econ 527 Spring 2026 Question 2
% Professor Greaney
% Written by: Camille Bergeron
% Due: April 230, 2026
% Last Edited: April 20, 2026

% Neoclassical Growth Model

%% Part 1: AR(1) approximation

% paramters
% for vfi
beta = 0.96; % discount factor
alpha = 0.35; % capital share
delta = 1; % depreciation rate
A = 1; % technology level

% for grid
curve = 2; % curvature parameter for polynomial grid
nodes = 100; % number of grid points
point_a = 0.1; % lower bound of the grid
point_b = 10; % upper bound of the grid

% creating grid points for z
function c_grid_poly = polynomial_grid(point_a, point_b, nodes, curve)
    c_grid_normalized = linspace(0, 1, nodes); % normalized grid in [0,1]
    c_grid_poly = point_a + (point_b - point_a) * c_grid_normalized.^curve; % polynomial transformation
end

% creating capital grid using polynomial transformation
K_grid = polynomial_grid(point_a, point_b, nodes, curve); % capital grid using polynomial transformation
disp('Capital grid points:');
disp(K_grid(1:10));

% cubic spline interpolation for the value function
function V_interp = cubic_spline_interpolation(K_grid, V_grid, K_query)
    spline_interp = spline(K_grid, V_grid); % create spline interpolation object
    V_interp = ppval(spline_interp, K_query); % evaluate spline interpolation at query points
end 

% displaying results
V_interp = cubic_spline_interpolation(K_grid, rand(1, nodes), K_grid); % example usage of cubic spline interpolation
disp('Interpolated value function at capital grid points:');
disp(V_interp(1:10));

% calculating the policy function 