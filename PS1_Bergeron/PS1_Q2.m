% Econ 527 Spring 2026 Question 2
% Professor Greaney
% Written by: Camille Bergeron
% Due: April 21, 2026
% Last Edited: April 11, 2026


%% Part 1: Interpolation

% vector for consumption grid
c_grid = linspace(0.1, 30, 100); % consumption grid from 0.1 to 30 with step of 0.01
disp('Consumption grid:');
disp(c_grid(1:10));

% evaluating utility function at the consumption grid points

function u = crra_utility(c, gamma)
    if gamma == 1
        u = log(c); % logarithmic utility for gamma = 1
    else
        u = c.^(1-gamma) / (1-gamma); % CRRA utility for gamma != 1
    end
end

gamma = 4; % risk aversion parameter
u_grid = crra_utility(c_grid, gamma); % utility grid using CRRA utility function
disp('Original utility values:');
disp(u_grid(1:10));

% spline interpolation
spline_interp = spline(c_grid, u_grid); % create spline interpolation object
c_grid_interp = linspace(0.1, 30, 1000); % new consumption grid for interpolation
u_interp = ppval(spline_interp, c_grid_interp); % evaluate spline interpolation at new dense consumption grid

% calculating error at each grid point 
u_original = crra_utility(c_grid_interp, gamma); % original utility values at the new dense grid
error = abs(u_original - u_interp); % absolute error between original and interpolated utility values
disp('Maximum error between original and interpolated utility values:');
disp(max(error));

% Plotting the original utility values and the interpolated values
figure;
fplot(@(x) crra_utility(x, gamma), [0.1, 30], 'r-', 'DisplayName', 'Original Utility'); % original utility function
hold on;
plot(c_grid_interp, u_interp, 'b-', 'DisplayName', 'Spline Interpolation'); % spline interpolation
xlabel('Consumption');
ylabel('Utility');
title('Utility Function and Spline Interpolation');
legend;
grid on;    
saveas(gcf,'spline_interpolation.fig')
saveas(gcf,'spline_interpolation.png')

%% Part 2: Polynomial grid

% creating polynomilial grid for consumption
curve = 3; % curvature parameter for polynomial grid
c_grid_poly = 0.1 + (30-0.1)* c_grid.^curve; % polynomial grid transformation
disp('Polynomial consumption grid:');
disp(c_grid_poly(1:10));    

u_grid_poly = crra_utility(c_grid_poly, gamma); % utility grid using CRRA utility function
disp('Original utility values:');
disp(u_grid_poly(1:10));

% spline interpolation, again
spline_interp_poly = spline(c_grid_poly, u_grid_poly); % create spline interpolation object
c_grid_interp_poly = linspace(0.1, 30, 1000); % new consumption grid for interpolation
u_interp_poly = ppval(spline_interp_poly, c_grid_interp_poly); % evaluate spline interpolation at new dense consumption grid

% calculating error at each grid point 
u_original_poly = crra_utility(c_grid_interp_poly, gamma); % original utility values at the new dense grid
error_poly = abs(u_original_poly - u_interp_poly); % absolute error between original and interpolated utility values
disp('Maximum error between original and interpolated utility values:');
disp(max(error_poly));

% Plotting the original utility values and the new interpolated values
figure;
fplot(@(x) crra_utility(x, gamma), [30, 40], 'r-', 'DisplayName', 'Original Utility'); % original utility function
hold on;
plot(c_grid_interp_poly, u_interp_poly, 'b-', ...
    'DisplayName', 'Spline Interpolationb Polynomial', Color='#DC94FF'); % spline interpolation   
xlabel('Consumption');
ylabel('Utility');
title('Utility Function and Spline Interpolation with Polynomial Grid, 30-40');
legend;
grid on;    
saveas(gcf,'spline_interpolation_poly1.png')
saveas(gcf,'spline_interpolation_poly.fig')

% Plotting the original utility values and the new interpolated values
figure;
fplot(@(x) crra_utility(x, gamma), [0.01, 0.1], 'r-', 'DisplayName', 'Original Utility'); % original utility function
hold on;
plot(c_grid_interp_poly, u_interp_poly, 'b-', ...
    'DisplayName', 'Spline Interpolationb Polynomial', Color='#DC94FF'); % spline interpolation   
xlabel('Consumption');
ylabel('Utility');
title('Utility Function and Spline Interpolation with Polynomial Grid, 0.01-0.1');
legend;
grid on;    
saveas(gcf,'spline_interpolation_poly2.png')
saveas(gcf,'spline_interpolation_poly2.fig')


%% Part 3: 
