% Econ 527 Spring 2026
% HW 4 Question 3 - Aiyagari Model
% Professor Greaney
% Written by: Camille Bergeron
% Due: May 22, 2026
%
% Solves the Aiyagari (1994) model:
%

clear; clc; close all;


%% ── Parameters ───────────────────────────────────────────────────────────

Params.gamma   = 2;       % CRRA risk aversion
Params.beta    = 0.96;    % annual discount factor
Params.rho     = 0.9;     % AR(1) persistence of log labor supply
Params.sigma   = 0.2;     % AR(1) innovation std dev
Params.n_z     = 7;       % income grid points
Params.alpha   = 0.36;    % capital share
Params.delta   = 0.10;    % depreciation rate
Params.e_stop  = 1e-6;    % VFI convergence tolerance
Params.max_iter = 2000;   % max VFI iterations

% Asset grid (borrowing limit = 0 per problem set)
Params.a_min  = 0;
Params.a_max  = 60;
Params.n_a    = 150;
Params.n_af   = 150;
Params.curve  = 2;        % polynomial grid curvature

% rouwenhorst for labor 
[l_grid, l_prob] = rouwenhorst(Params.rho, Params.sigma, Params.n_z);

% calculating the income states and distribution
labor_grid = exp(l_grid); % income states
labor_dist = mc_invdist(l_prob); % income distribution                      
expected_labor_supply = labor_grid' * labor_dist; % expected earnins           


% creating a_next_grid because EGM
a_grid = polynomial_grid(Params.a_min, Params.a_max, Params.n_a, Params.curve);
a_fine_grid = polynomial_grid(Params.a_min, Params.a_max, Params.n_af, Params.curve);

%% Functions
% ── polynomial grid ──────────────────────────────────────────────────────────────────
polynomial_grid = @(a_min, a_max, n_a, curve) ...
    a_min + (a_max - a_min) * linspace(0, 1, n_a) .^ curve;

% ── Utilty function ──────────────────────────────────────────────────────────────────
u_fxn = @(c) c.^(1-Params.gamma)/(1-Params.gamma); 

% ── derivatove of the inverse ──────────────────────────────────────────────────────────────────
u_prime_inv = @(mu) mu .^ (-1 / Params.gamma);

%% Solving model

if isfile('results/aiyagari_model.mat')
    disp('File exists and skipping aiyagari model calculation...');
    disp('Loading file instead...');
    load('results/aiyagari.mat', 'V_grid', 'policy_fxn', 'phi_dist', 'net_assets');
    disp('Loaded aiyagari model');
else
    disp('File not found, calculating aiyagari model...')

    tic;
    damp = 0.3;  % Damping parameter
    k = 5;  % inital guess for capital stock
    dif = inf; 

    % looping to find the correct levels
    for iter = 1:Params.max_iter
        r = Params.alpha * k^(Params.alpha-1) * expected_labor_supply^(1-Params.alpha) - Params.delta;

        % solving for implied wage
        w = (1-Params.alpha) * k^Params.alpha * expected_labor_supply^-Params.alpha;   

        % solving model, adjusting wage (previously was  1)
        [V_grid, policy_fxn, phi_dist, net_assets] = f_spe(r, l_grid*w,...
            l_prob, u_fxn, u_prime_inv, Params, a_grid, a_fine_grid);
        
        % calculating error (should be 0)
        diff = abs(net_assets - k);  
    
        if diff < Params.e_stop
            break
        end
        k = k + damp*(net_assets-k); 
    end
    assert(iter < Params.max_iter) % for iterations

    time = toc;

    fprintf('Took %.4d seconds to converge', time)
    save('results/aiyagari_model.mat', 'V_grid', 'policy_fxn', 'phi_dist', 'net_assets');
end

%% Computing Inequality
[gini, LC] = gini_lc(sum(phi_dist), a_next_grid);
t10_share = 1 - interp1(LC(:,1), LC(:,2), 0.9); 
t1_share = 1 - interp1(LC(:,1), LC(:,2), 0.99); 

