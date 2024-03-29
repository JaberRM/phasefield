function GradPhi = Calc_Grad(Phi,dx,dy,Flag)
% Flag = 1 : forward diff
% Flag = -1 : backward diff
% Flag = 0 : central

% Periodic

if Flag == 1
    
    PhiFX = circshift(Phi,1,2);
    PhiFY = circshift(Phi,1,1);
    
    GradPhi(:,:,1)=(PhiFX - Phi) / (dx);
    GradPhi(:,:,2)=(PhiFY - Phi) / (dy);
    
elseif Flag == -1
    
    PhiBX = circshift(Phi,-1,2);
    PhiBY = circshift(Phi,-1,1);
    
    GradPhi(:,:,1)=(Phi - PhiBX) / (dx);
    GradPhi(:,:,2)=(Phi - PhiBY) / (dy);
    
elseif Flag == 0 
    
    PhiFX = circshift(Phi,1,2);
    PhiBX = circshift(Phi,-1,2);
    PhiFY = circshift(Phi,1,1);
    PhiBY = circshift(Phi,-1,1);
    
    GradPhi(:,:,1)=(PhiFX - PhiBX) / (2*dx);
    GradPhi(:,:,2)=(PhiFY - PhiBY) / (2*dy);
    
end