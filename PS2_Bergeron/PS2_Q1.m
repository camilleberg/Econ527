% Econ 527 Spring 2026 
% HW 2 Question 1
% Professor Greaney
% Written by: Camille Bergeron
% Due: April 30, 2026
% Last Edited: April 20, 2026


% Neoclassical Growth Model

%% Part 1: Neoclassical Growth Model Setup

% paramters
% for vfi
Params.beta = 0.96; % discount factor
Params.alpha = 0.35; % capital share
Params.delta = 1; % depreciation rate
Params.A = 1; % technology level


%% Anlaytical solution for the value function and policy function
b = Params.alpha / (1 - Params.alpha * Params.beta);
a = (1/(1-Params.beta)) * log((Params.A/(1+Params.beta*b)) * ((Params.A*Params.beta*b)/(1+Params.beta*b))^(Params.beta*b));
disp('Analytical value function: V(k) = a + b*log(k)');
disp(['a = ', num2str(a)]);
disp(['b = ', num2str(b)]);

% creating funcion for the analytical value function
function V = analytical_value_function(k, a, b)
    V = a + b * log(k); % analytical value function
end 



%% Numerical solution using value function iteration

% parameters
% for grid
Params.curve = 2; % curvature parameter for polynomial grid
Params.nodes = 100; % number of grid points
Params.point_a = 0.1; % lower bound of the grid
Params.point_b = 10; % upper bound of the grid

% for iteration 
Params.e_stop = 1e-4; % convergence criterion

% creating capital grid using polynomial transformation
K_grid = polynomial_grid(Params.point_a, Params.point_b, Params.nodes, Params.curve); % capital grid using polynomial transformation
disp('Capital grid points:');
disp(K_grid(1:10));

% loop
% evaluate inital guess with interpolatiion object 
% check absolute max error between old and new value function grids
% vaerify convergence and reassign as diff

% value function iteration
V_grid = zeros(size(K_grid)); % initialize value function grid
policy_grid = zeros(size(K_grid)); % initialize policy function grid
iteration = 0; % iteration counter
error = Inf; % initialize error

% making consumption function
function c = consumption(k, k_next, Params)
    c_possible = Params.A * k^Params.alpha - k_next + (1-Params.delta) * k; % possible consumption based on current capital, next period's capital, and production
    c = max(0, c_possible); % ensure consumption is non-negative
end

% utility function
function u = utility(c, Params)
    if c > 0
        u = log(c); % utility from consumption using log utility function
    else
        u = -Inf; % assign negative infinity utility for non-positive consumption
    end
end 

% running through loop 
while error > Params.e_stop * (1-Params.beta) % continue until convergence
    V_grid_old = V_grid; % store old value function grid
    for i = 1:length(K_grid)
        k = K_grid(i); % current capital stock
        % compute the value of consuming all output and investing the rest
        c = consumption(k, K_grid, Params); % consumption for each possible next period capital stock
        u = utility(c, Params); % utility from consumption for each possible next period capital stock

        % calculating new value function for each possible next period capital stock
        V_next = cubic_spline_interpolation(K_grid, V_grid_old, K_grid); % interpolate value function for next period capital stocks
        total_value = u + Params.beta * V_next; % total value for each possible next period capital stock

        % updating value function and policy function
        [V_grid(i), policy_index] = max(total_value); % update value function and policy function
        policy_grid(i) = K_grid(policy_index); % store optimal next period capital stock in policy grid
    end
    error = max(abs(V_grid - V_grid_old)); % compute maximum error between iterations
    iteration = iteration + 1; % increment iteration counter
end

disp(['Value function iteration converged in ', num2str(iteration), ' iterations.']);
disp('Optimal policy function (next period capital stock):');
disp(policy_grid(1:10));

%% Plotting the results
% Plot value function
figure;
fig = figure;
theme(fig, "light");
plot(K_grid, V_grid, 'b-', 'LineWidth', 2);
hold on;
plot(K_grid, analytical_value_function(K_grid, a, b), 'r--', 'LineWidth', 2);
xlabel('Capital Stock (k)');
ylabel('Value Function (V)');
title('Value Function Iteration vs Analytical Solution');
legend('Numerical VFI', 'Analytical Solution');
grid on;    
saveas(gcf, 'figs/Value_Function_Comparison.png');