% Econ 527 Spring 2026
% Professor Greaney
% Written by: Camille Bergeron
% Due: April 21, 2026
% Last Edited: April 11, 2026


%% Part 2: Interpolation

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

% error = max(abs(u_grid - u_interp)); % calculate maximum error between original and interpolated utility values
% disp('Maximum error between original and interpolated utility values:');
% disp(error);    
