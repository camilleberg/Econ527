
function [V, A1, PHI, nad] = f_spe(r, ZG, par)

% [V, A1, PHI, nad] = f_spe(r, par)
%
% Given an interest rate and parameters, solves the stationary partial 
% equilibrium of the income fluctuation problem using the endogenous 
% gridpoint method. Returns value function V, policy function A1, 
% stationary density of state variables PHI and net asset demand nad.

%--------------------------------------------------------------------------
% Unpack parameters, set new parameters, initialize objects for later use
%--------------------------------------------------------------------------

PI=par.PI; u=par.u; du_inv=par.du_inv; bet=par.bet; AG=par.AG; AGf=par.AGf;

na = length(AG);

nz = length(ZG);

naf = length(AGf);

maxit = 1000;

tol = (1-bet) * 10^-4;

V = NaN(nz,na);  

dV0 = NaN(nz,na);

A1 = NaN(nz,naf);

pp = cell(nz,1); 

pp_d1 = cell(nz,1); 

%--------------------------------------------------------------------------
% Solve Bellman equation
%--------------------------------------------------------------------------

C_AUX = max(r,0.001) * (AG-AG(1)+0.01) + ZG;

V0 = u(C_AUX) / (1-bet);                            % Initial guess for V

for iter = 1:maxit
    
   % Make cubic splines and compute the derivative of V0 on the grid, dV0
   
   for iz = 1:nz    
      
      pp{iz} = spline(AG, V0(iz,:));
      
      pp_d1{iz} = fnder(pp{iz}, 1);
      
      dV0(iz,:) = ppval(pp_d1{iz}, AG);
      
   end
   
   % Compute expected V', expected dV', and endogenous gridpoints
   
   EV = PI * V0;
   
   EDV = PI * dV0;
   
   AE = (AG - ZG + du_inv(bet*EDV)) / (1+r);
   
   % Update value function 
   
   VE = u((1+r)*AE + ZG - AG) + bet*EV;          % V(endogenous gridpoints)
   
   for iz = 1:nz
       
      V(iz,:) = spline(AE(iz,:), VE(iz,:), AG);  % Unconstrained, exog grid  
      
      IND = AG < AE(iz,1);                       % Constrained, exog grid
      
      V(iz,IND) = u((1+r)*AG(IND) + ZG(iz) - AG(1)) + bet*EV(iz,1);

   end
   
   % Check for convergence
   
   dif = max(max(abs(V - V0)));  % disp(dif)
   
   if dif < tol
      break
   end
   
   V0 = V;
    
end

assert(iter < maxit)

%--------------------------------------------------------------------------
% Compute policy function on fine exogenous wealth grid
%--------------------------------------------------------------------------

for iz = 1:nz
       
   A1(iz,:) = spline(AE(iz,:), AG, AGf);
      
   IND = AG < AE(iz,1);
      
   A1(iz,IND) = AG(1);

end

assert(min(min(A1)) >= AG(1))

% assert(max(max(A1)) <= AG(na))

A1 = min(A1, AG(na));                           % Policy. Impose a' <= amax

%--------------------------------------------------------------------------
% Make transition matrix TM
%--------------------------------------------------------------------------

TM = zeros(naf * nz);

for ia = 1:naf
for iz = 1:nz
    
   j = min(sum(AGf <= A1(iz,ia)), naf-1); 
   
   wgt = (A1(iz,ia) - AGf(j)) / (AGf(j+1) - AGf(j));  
   
   assert(0<=wgt && wgt<=1)
   
   s = (iz-1)*naf + ia;
   
   for iz1 = 1:nz
   
      s1_lo = (iz1-1)*naf + j;
   
      TM(s,s1_lo) = (1-wgt) * PI(iz,iz1); 
      
      TM(s,s1_lo+1) = wgt * PI(iz,iz1); 
      
   end
    
end
end

assert(max(abs(sum(TM,2) - 1)) < 1E-12)

PHI_vec = mc_invdist(TM);

PHI = reshape(PHI_vec,naf,nz)';                     % Invariant dist.

nad = sum(sum(AGf.*PHI));                           % Net asset demand

end

