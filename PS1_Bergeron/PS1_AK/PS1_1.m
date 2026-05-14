
clearvars;  clc;  addpath('/Users/bgreaney/Dropbox/Codes');

rng('default');

%--------------------------------------------------------------------------
% Set parameters
%--------------------------------------------------------------------------

rho = 0.99;

sig = 0.2;

uvar = sig^2 / (1-rho^2);               % Unconditional variance

acv = uvar * rho;                       % Autocovariance

cvar = sig^2;                           % Conditional variance

%--------------------------------------------------------------------------
% Discretize AR(1)
%--------------------------------------------------------------------------

N = 9;

m = 3;

[X, PI] = tauchen(N, 0, rho, sig, m);

% [X, PI] = rouwenhorst(rho, sig, N);

%--------------------------------------------------------------------------
% Simulate
%--------------------------------------------------------------------------

mc = dtmc(PI);

P = invdist_markov(PI);

II0 = ceil(1000 * P);

II = simulate(mc, 1000, 'X0', II0);

%--------------------------------------------------------------------------
% Compute statistics. Xt = X(t), Xt1 = X(t+1)
%--------------------------------------------------------------------------

Xt = X(II(1:end-1,:));  Xt = Xt(:);

Xt1 = X(II(2:end,:));   Xt1 = Xt1(:);

rho_sim = regress(Xt1, Xt);

cov_mat = cov(Xt1, Xt);

uvar_sim = cov_mat(1,1);

acv_sim = cov_mat(1,2);

cvar_sim = var(Xt1 - rho_sim*Xt);

RATIOS = [rho_sim/rho; uvar_sim/uvar; acv_sim/acv; cvar_sim/cvar];
