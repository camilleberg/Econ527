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


% loop
% evaluate inital guess with interpolatiion object 
% check absolute max error between old and new value function grids
% vaerify convergence and reassign as diff


% running through loop 
function[K_grid, V_grid, policy_grid]= vfi_loop(Params, print_iter)

        % creating capital grid using polynomial transformation
    K_grid = polynomial_grid(Params.point_a, Params.point_b, Params.nodes, Params.curve); % capital grid using polynomial transformation
    disp('Capital grid points:');
    disp(K_grid(1:10));

    % value function iteration
    V_grid = zeros(size(K_grid)); % initialize value function grid
    policy_grid = ones(size(K_grid)); % initialize policy function grid
    iteration = 0; % iteration counter
    err = Inf; % initialize error


    % runnign thrgouh loop 
    while err > Params.e_stop * (1 - Params.beta)
        V_grid_old = V_grid; % saving old grid

        % looping oiver every state
        for i = 1:length(K_grid)
            k = K_grid(i);
            c = consumption(k, K_grid, Params); % claculating consumption for all possible next period capital stocks
            u = utility(c); % same for utility 

            % interpolating for next choice of k' 
            V_next = cubic_spline_interpolation(K_grid, V_grid_old, K_grid); 
            total_value = u + Params.beta * V_next; % calculating total value for each possible next period capital stock

            [V_grid(i), idx] = max(total_value); % gets max value for each choice of k' and the index of the optimal choice
            policy_grid(i) = K_grid(idx); % saves policy grid 
        end

        err = max(abs(V_grid - V_grid_old)); % check error 
        iteration = iteration + 1;

        % displaying every 50 
        if mod(iteration, 50) == 0 && print_iter == true
            disp(['Iteration: ', num2str(iteration), ', Error: ', num2str(err)]);
        end
    end


    disp(['Value function iteration converged in ', num2str(iteration), ' iterations.']);
    disp('Optimal policy function (next period capital stock):');
    disp(policy_grid(1:10));

end

[K_grid, V_grid, policy_grid] = vfi_loop(Params, true); % run the value function iteration loop


%% Checking timing 
executionTime = timeit(@() vfi_loop(Params, false));
disp(['Execution time for value function iteration: ', num2str(executionTime), ' seconds.']);

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

%% Local Functions

% making consumption function
function c = consumption(k, k_next, Params)
    c_possible = Params.A * k^Params.alpha - k_next + (1-Params.delta) * k; % possible consumption based on current capital, next period's capital, and production
    c = max(0, c_possible); % ensure consumption is non-negative
end

% utility function
function u = utility(c)
    u = -Inf(size(c));        % initialize all to -Inf
    u(c > 0) = log(c(c > 0)); % only compute log where feasible
end
