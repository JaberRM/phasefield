function [T, Tol, Iter,converged] = Update_TH(T, Phi, dPhidt, Param, Mat, Flag_Land)


K   =  Mat.K_ice   * (1-Phi)   +  ( Mat.K_land * Flag_Land + Mat.K_water * (1-Flag_Land) ) .* (Phi);
CP  =  Mat.Cp_ice  * (1-Phi)  +  ( Mat.Cp_land * Flag_Land + Mat.Cp_water * (1-Flag_Land) ) .* (Phi);
Rho =  Mat.Rho_ice * (1-Phi) +  ( Mat.Rho_land * Flag_Land + Mat.Rho_water * (1-Flag_Land) ) .* (Phi);

% K = Mat.K_ice;
% CP = Mat.Cp_ice;
% Rho = Mat.Rho_ice;

f = ( Calc_Force_TH(T,K,Param.dx,Param.dy)  - Mat.Q * dPhidt) ./ (Rho .* CP);

Ttemp = T + f * Param.dt;
converged = 0;
for Iter = 1:10000
    
    ftemp = ( Calc_Force_TH(Ttemp,K,Param.dx,Param.dy) - Mat.Q * dPhidt ) ./ (Rho .* CP);
    
    Ttry = T + (f + ftemp) * Param.dt/2;
%     Ttry = ApplyBC(Ttry);
    
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


