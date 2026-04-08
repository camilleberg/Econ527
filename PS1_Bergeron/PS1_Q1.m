% Econ 527 Spring 2026
% Professor Greaney
% Written by: Camille Bergeron
% Due: April 21, 2026
% Last Edited: April 8, 2026

%% Part 1: AR(1) approximation
n_z = 9; % number of grid points for z
rho = 0.99; % persistence parameter
sigma = 0.2; % standard deviation of the error term
q = 3; % number of standard deviations to cover in the grid
z_grid = linspace(-q*sigma, q*sigma, n_z); % grid for z w 3 st devs
disp('Grid for z:');
disp(z_grid);



