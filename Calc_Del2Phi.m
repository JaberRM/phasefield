function [Del2Phi] = Calc_Del2Phi(Phi, dx,dy)


% Periodic
PhiFX = circshift(Phi,1,2);
PhiBX = circshift(Phi,-1,2);
PhiFY = circshift(Phi,1,1);
PhiBY = circshift(Phi,-1,1);

Del2Phi=( PhiFX + PhiFY + PhiBX + PhiBY - 4*Phi) / (dx*dy);

 

% % nonPer
% [dPhidY, dPhidX] = gradient(Phi,hx,hy);
% Del2Phi = del2(Phi,hx,hy);
% Del4Phi = del2(Del2Phi,hx,hy);
% 
% GradPhi2 = dPhidX.^2 + dPhidY.^2;
