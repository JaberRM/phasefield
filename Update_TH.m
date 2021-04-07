function [T, Tol, Iter,converged] = Update_TH(T, Phi, Param, TH, Flag_Land)

K   =  TH.K_ice * Phi   +  ( TH.K_land * Flag_Land + TH.K_water * (1-Flag_Land) ) .* (1-Phi);
CP  =  TH.Cp_ice * Phi  +  ( TH.Cp_land * Flag_Land + TH.Cp_water * (1-Flag_Land) ) .* (1-Phi);
Rho =  TH.Rho_ice * Phi +  ( TH.Rho_land * Flag_Land + TH.Rho_water * (1-Flag_Land) ) .* (1-Phi);

% K = TH.K_ice;
% CP = TH.Cp_ice/10;
% Rho = TH.Rho_ice;

f = Calc_Force_TH(T,K,Param.dx,Param.dy) ./ (Rho .* CP);

Ttemp = T + f * Param.dt;
converged = 0;
for Iter = 1:10000
    
    ftemp = Calc_Force_TH(Ttemp,K,Param.dx,Param.dy) ./ (Rho .* CP);
    
    Ttry = T + (f + ftemp) * Param.dt/2;
%     C1 = ApplyBC(C1);
    
    dT = Ttry - Ttemp;
    Tol = max(abs(dT(:))) ;
    if Tol < 1e-6
        T = Ttry;
        converged = 1;
        break
    else
        Ttemp = Ttry;
    end
end


