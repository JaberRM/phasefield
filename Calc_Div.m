function DivPhi = Calc_Div(Phi,dx,dy,Flag)

Phi1 = Phi(:,:,1);
Phi2 = Phi(:,:,2);

if Flag==0
    
    % Periodic
    % Phi1FX = circshift(Phi1,1,2);
    Phi1BX = circshift(Phi1,-1,2);
    % Phi2FY = circshift(Phi2,1,1);
    Phi2BY = circshift(Phi2,-1,1);
    
    DivPhi = (Phi1 - Phi1BX) / (dx) + (Phi2 - Phi2BY) / (dy);
else
    % Periodic
    Phi1FX = circshift(Phi1,1,2);
%     Phi1BX = circshift(Phi1,-1,2);
    Phi2FY = circshift(Phi2,1,1);
%     Phi2BY = circshift(Phi2,-1,1);
    
    DivPhi = (Phi1FX - Phi1) / (dx) + (Phi2FY - Phi2) / (dy);
    
end
