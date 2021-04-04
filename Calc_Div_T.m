function DivT = Calc_Div_T(T,dx,dy,Flag)
% Flag = 1 : forward diff
% Flag = -1 : backward diff
% Flag = 0 : central

T1 = T(:,:,1);
T2 = T(:,:,2);
T3 = T(:,:,3);
DivT = zeros( size(T,1) , size(T,2) ,2);

if Flag==-1
    
    % Periodic
    % Phi1FX = circshift(Phi1,1,2);
    T1BX = circshift(T1,-1,2);
    T3BX = circshift(T3,-1,2);

    % Phi2FY = circshift(Phi2,1,1);
    T2BY = circshift(T2,-1,1);
    T3BY = circshift(T3,-1,1);

    
    DivT(:,:,1) = (T1 - T1BX) / (dx) + (T3 - T3BY) / (dy);
    DivT(:,:,2) = (T3 - T3BX) / (dx) + (T2 - T2BY) / (dy);

elseif Flag ==1
%     % Periodic
%     Phi1FX = circshift(T1,1,2);
% %     Phi1BX = circshift(Phi1,-1,2);
%     Phi2FY = circshift(T2,1,1);
% %     Phi2BY = circshift(Phi2,-1,1);
%     
%     DivT = (Phi1FX - T1) / (dx) + (Phi2FY - T2) / (dy);
 
elseif Flag == 0
%     Phi1FX = circshift(T1,1,2);
%     T1BX = circshift(T1,-1,2);
%     Phi2FY = circshift(T2,1,1);
%     T2BY = circshift(T2,-1,1);
%     
%     DivT = (Phi1FX - T1BX) / (2*dx) + (Phi2FY - T2BY) / (2*dy);
end