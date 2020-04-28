function Phi = Update_CH(Phi,M,EPS,dx,dy,dt,MethodFlag)

f = Calc_DrivingFocre(Phi,M,EPS,dx,dy,MethodFlag);

Phi0 = Phi + f * dt;
% Phi0 = ApplyBC(Phi0);

for Iter = 1:10000
    
    f0 = Calc_DrivingFocre(Phi0,M,EPS,dx,dy,MethodFlag);
    
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


