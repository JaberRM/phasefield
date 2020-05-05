function dfdPhi = Calc_dfdPhi(Phi)
% 

% f= Phi.^2 * (5 * Phi.^2 -12*Phi + 8)
% dfdPhi = 2 * Phi .* ( 3 * Phi.^2 -8 * Phi + 6) + Phi.^2 .* ( 6 * Phi - 8);

% %  f  = phi.^2 .* (3-2*phi)
dfdPhi = 2 * Phi .* (3 - 2 * Phi) - 2 * Phi.^2 ;
