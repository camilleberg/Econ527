
clearvars;  clc;  addpath('/Users/bgreaney/Dropbox/Codes');

%--------------------------------------------------------------------------
% Set parameters and make grids
%--------------------------------------------------------------------------

gam = 4;

u = @(c) c.^(1-gam)/(1-gam);

cmin = 0.1;  

cmax = 30;

% Ci = cmin + (cmax-cmin)*linspace(0,1,100);          % Evenly spaced grid 

Ci = cmin + (cmax-cmin)*linspace(0,1,100).^3;     % Unevenly spaced grid 

% Xj = linspace(0.1, 30, 1000);                     % Interpolation tests

% Xj = linspace(30, 40, 100);                         % Extrapolation test 1

Xj = linspace(0.01, 0.1, 100);                      % Extrapolation test 2

%--------------------------------------------------------------------------
% Make cubic spline, compute max error, and plot u(Xj) vs g(Xj)
%--------------------------------------------------------------------------

pp = spline(Ci, u(Ci));

Gj = ppval(pp, Xj);

Fj = abs(u(Xj) - Gj);

fprintf('max error = %2.4g \n\n',max(Fj));

plot(Xj,u(Xj),Xj,Gj,'--','Linewidth',2);  
legend('u(xj)','g(xj)','Location','se','FontSize',12)

% saveas(gcf,'ps1_fig1.png')                      % Even cj, interpolation
% saveas(gcf,'ps1_fig2.png')                      % Uneven cj, interpolation
% saveas(gcf,'ps1_fig3.png')                      % Extrapolation 1
saveas(gcf,'ps1_fig4.png')                      % Extrapolation 2

%--------------------------------------------------------------------------
% Check whether g(Xj) is montonic, strictly concave, and < 0
%--------------------------------------------------------------------------

if min(diff(Gj)) <= 0
   disp('g(xj) is not monotonically increasing')
end

D1 = (Gj(2:end-1) - Gj(1:end-2)) ./ (Xj(2:end-1) - Xj(1:end-2));

D2 = (Gj(3:end) - Gj(2:end-1)) ./ (Xj(3:end) - Xj(2:end-1));

if min(D1 - D2) <= 0
   disp('g(xj) is not strictly concave')
end

if max(Gj) >= 0
   disp('g(xj) is not less than 0')
end

