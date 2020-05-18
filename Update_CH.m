function [C, Tol, Iter] = Update_CH(C,Phi,M,EPS,dx,dy,dt,MethodFlag)

f = Calc_Force_CH(C,Phi,M,EPS,dx,dy,MethodFlag);

C0 = C + f * dt;
% C0 = ApplyBC(C0);

for Iter = 1:10000
    
    f0 = Calc_Force_CH(C0,Phi,M,EPS,dx,dy,MethodFlag);
    
    C1 = C + (f + f0) * dt/2;
%     C1 = ApplyBC(C1);
    
    dC = C1 - C0;
    Tol = max(abs(dC(:))) ;
    if Tol < 1e-6
        C = C1;
        break
    else
        C0 = C1;
    end
end


