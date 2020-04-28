function GradPhi = Calc_Grad(Phi,dx,dy,Flag)


% Periodic

if Flag == 0
    
    PhiFX = circshift(Phi,1,2);
    % PhiBX = circshift(Phi,-1,2);
    PhiFY = circshift(Phi,1,1);
    % PhiBY = circshift(Phi,-1,1);
    
    GradPhi(:,:,1)=(PhiFX - Phi) / (dx);
    GradPhi(:,:,2)=(PhiFY - Phi) / (dy);
else
    
    %     PhiFX = circshift(Phi,1,2);
    PhiBX = circshift(Phi,-1,2);
    %     PhiFY = circshift(Phi,1,1);
    PhiBY = circshift(Phi,-1,1);
    
    GradPhi(:,:,1)=(Phi - PhiBX) / (dx);
    GradPhi(:,:,2)=(Phi - PhiBY) / (dy);
end