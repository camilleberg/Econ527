function u_inv = utility_inv(c, Params); % (Na x 1)  

% u_inv = utility_inv(c, Params); 
%
% calcualtes inverse of CRRA utility
% 
% c: consumption level (scalar or vector)
% Params: struc of model

    u_inv = c.^(-Params.gamma); % inverse of CRRA utility
end
