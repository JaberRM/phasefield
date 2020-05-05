function [Phi, Tol, Iter] = Update_AC(Phi,C ,M,EPS,dx,dy,dt)

f = Calc_Force_AC(Phi,C,M,EPS,dx,dy);

Phi0 = Phi + f * dt;
% Phi0 = ApplyBC(Phi0);

for Iter = 1:10000
    
    f0 = Calc_Force_AC(Phi0,C,M,EPS,dx,dy);
    
    Phi1 = Phi + (f + f0) * dt/2;
%     Phi1 = ApplyBC(Phi1);
    
    dPhi = Phi1- Phi0;
    Tol = max(abs(dPhi(:))) ;
    if Tol < 1e-6
        Phi = Phi1;
        break
    else
        Phi0 = Phi1;
    end
end


