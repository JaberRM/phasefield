function [Phi, Tol, Iter,converged] = Update_AC(Phi,C ,M,EPS,Param)

f = Calc_Force_AC(Phi,C,M,EPS,Param.dx,Param.dy);

Phi0 = Phi + f * Param.dt;
% Phi0 = ApplyBC_AC(Phi0);
converged = 0;
for Iter = 1:10000
    
    f0 = Calc_Force_AC(Phi0,C,M,EPS,Param.dx,Param.dy);
    
    Phi1 = Phi + (f + f0) * Param.dt/2;
%     Phi1 = ApplyBC_AC(Phi1);
    
    dPhi = Phi1- Phi0;
    Tol = max(abs(dPhi(:))) ;
    if Tol < 1e-6
        Phi = Phi1;
        converged = 1;
        break
    else
        Phi0 = Phi1;
    end
end


