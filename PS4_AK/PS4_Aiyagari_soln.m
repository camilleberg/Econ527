%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ECON 527 Problem Set 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;

clc;

addpath('/Users/bgreaney/Dropbox/Codes');

%--------------------------------------------------------------------------
% Set parameters 
%--------------------------------------------------------------------------

gam = 2;                                        % Risk aversion

bet = 0.96;                                     % Discount factor

rho = 0.9;                                      % AR(1) persistence

sig = 0.2;                                      % AR(1) volatility

alp = 0.36;                                     % Capital share

del = 0.10;                                     % Capital depreciation rate

u = @(c) c.^(1-gam)/(1-gam);                    % Utility function

du_inv = @(c) c.^(-1/gam);                      % Inverse of u'(c)

%--------------------------------------------------------------------------
% Discretize state space
%--------------------------------------------------------------------------

na = 150;                                       % # wealth gridpoints

acur = 2;                                       % Wealth grid curvature

nl = 7;                                         % # income gridpoints

[log_LG, PI] = rouwenhorst(rho, sig, nl);

LG = exp(log_LG);                               % Labor supply states

PL = mc_invdist(PI);                            % Labor supply distribution

el = LG' * PL;                                  % Expected labor supply

amin = 0;

amax = 60;

AG = amin + (amax-amin).*(linspace(0,1,na)).^acur;

AGf = AG;                                       % Fine wealth grid

par = struct('PI',PI,'bet',bet,'u',u,'du_inv',du_inv,'AG',AG,'AGf',AGf);

%--------------------------------------------------------------------------
% Compute equilibrium using fixed point iteration
%--------------------------------------------------------------------------

damp = 0.3;                                     % Damping parameter

maxit = 100;

tol = 1E-5;

k = 6;                                          % Initial guess 

for iter = 1:maxit
    
   r = alp * k^(alp-1) * el^(1-alp) - del;
   
   w = (1-alp) * k^alp * el^-alp;   
    
   [V, A1, PHI, s] = f_spe(r, w*LG, par);
   
   dif = abs(s - k);  % disp(dif);
   
   if dif < tol
      break
   end
   
   k = k + damp*(s-k);
    
end

assert(iter < maxit)

save('PS4_eqbm','V','PHI','k');

%--------------------------------------------------------------------------
% Compute inequality statistics
%--------------------------------------------------------------------------

[gini, LC] = gini_lc(sum(PHI), AG);

t10_share = 1 - interp1(LC(:,1), LC(:,2), 0.9); 

t1_share = 1 - interp1(LC(:,1), LC(:,2), 0.99); 

%--------------------------------------------------------------------------
% Compute excess demand function
%--------------------------------------------------------------------------

nr = 25;

r_min = -0.02;

r_max = 0.04;

RG = linspace(r_min, r_max, nr)';

NAD = NaN(nr,1);                                % Net asset demand

for ir = 1:nr
    
   k = ((RG(ir) + del) / (alp * el^(1-alp))) ^ (1/(alp-1));
   
   w = (1-alp) * k^alp * el^-alp;
    
   [~, ~, ~, s] = f_spe(RG(ir), w*LG, par);
   
   NAD(ir) = s - k;  
   
end

%--------------------------------------------------------------------------
% Display results
%--------------------------------------------------------------------------

subplot(2,2,1);  plot(AG,V,'-o');    title('Value Function');
subplot(2,2,2);  plot(AG,A1,'-o');   title('Policy Function');
subplot(2,2,3);  plot(AG,PHI,'-o');  title('Density Function');
subplot(2,2,4);  plot(RG,NAD,'-o');  title('Excess Demand Function') 

saveas(gcf,'fig2.png') 

i_maxa = max(sum(A1>AG,2));

fprintf('\nEquilibrium interest rate: r = %2.6g \n\n', r)

fprintf('Decumulate assets if wealth > %2.6g (amax = %2.6g)\n\n', ...
   AG(i_maxa), AG(end))

fprintf('Wealth gini, top 10%%, and top 1%% wealth shares: %2.6g, %2.6g, and %2.6g\n',...
   gini, t10_share, t1_share);
