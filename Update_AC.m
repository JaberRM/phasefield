function [Phi, Tol, Iter, converged] = Update_AC(Phi,x , Param)

f = Calc_Force_AC(Phi,x,Param);

Phi0 = Phi + f * Param.dt;

converged = 0;
for Iter = 1:400
    
    f0 = Calc_Force_AC(Phi0,x,Param);
    
    Phi1 = Phi + (f + f0) * Param.dt/2;
    
    dPhi = Phi1 - Phi0;
    Tol = max(abs(dPhi(:))) ;
    if Tol < 1e-4
        Phi = Phi1;
        converged = 1;
        break
    else
        Phi0 = Phi1;
    end
end


