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

u = @(c) c.^(1-gam)/(1-gam);                    % Utility function

du_inv = @(c) c.^(-1/gam);                      % Inverse of u'(c)

%--------------------------------------------------------------------------
% Discretize state space
%--------------------------------------------------------------------------

na = 150;                                       % # wealth gridpoints

acur = 2;                                       % Wealth grid curvature

nz = 7;                                         % # income gridpoints

[log_ZG, PI] = rouwenhorst(rho, sig, nz);

ZG = exp(log_ZG);                               % Earning states

PZ = mc_invdist(PI);                            % Earning distribution

ez = ZG' * PZ;                                  % Expected earnings

amin = -ez;                                     

amax = 50*ez;

AG = amin + (amax-amin).*(linspace(0,1,na)).^acur;

AGf = AG;                                       % Fine wealth grid

par = struct('PI',PI,'bet',bet,'u',u,'du_inv',du_inv,'AG',AG,'AGf',AGf);

%--------------------------------------------------------------------------
% Compute equilibrium 
%--------------------------------------------------------------------------

opts = optimset('Display','iter');

r = fminbnd(@(r) f_nad(r, ZG, par)^2, -0.02, 0.02, opts);

[V, A1, PHI, nad_eqb] = f_spe(r, ZG, par);

%--------------------------------------------------------------------------
% Compute excess demand function
%--------------------------------------------------------------------------

nr = 25;                                        

r_min = -0.02;

r_max = 0.02;

RG = linspace(r_min, r_max, nr);                % Grid of interest rates                                   

NAD = NaN(nr,1);                                % Net asset demands

for ir = 1:nr
   [~, ~, ~, NAD(ir)] = f_spe(RG(ir), ZG, par);
end

%--------------------------------------------------------------------------
% Display results
%--------------------------------------------------------------------------

subplot(2,2,1);  plot(AG,V,'-o');    title('Value Function');
subplot(2,2,2);  plot(AG,A1,'-o');   title('Policy Function');
subplot(2,2,3);  plot(AG,PHI,'-o');  title('Density Function');
subplot(2,2,4);  plot(RG,NAD,'-o');  title('Excess Demand Function') 

saveas(gcf,'fig1.png')

i_maxa = max(sum(A1>AG,2));

fprintf('\nEquilibrium interest rate: r = %2.6g \n\n', r)

fprintf('Decumulate assets if wealth > %2.6g (amax = %2.6g)\n\n', ...
   AG(i_maxa), AG(end))
