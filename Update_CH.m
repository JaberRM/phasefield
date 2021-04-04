function [C, Tol, Iter,converged] = Update_CH(C,Phi,M,EPS,Param)

f = Calc_Force_CH(C,Phi,M,EPS,Param.dx,Param.dy);

Ctemp = C + f * Param.dt;
% C0 = ApplyBC(C0);
converged = 0;
for Iter = 1:10000
    
    ftemp = Calc_Force_CH(Ctemp,Phi,M,EPS,Param.dx,Param.dy);
    
    Ctry = C + (f + ftemp) * Param.dt/2;
%     C1 = ApplyBC(C1);
    
    dC = Ctry - Ctemp;
    Tol = max(abs(dC(:))) ;
    if Tol < 1e-6
        C = Ctry;
        converged = 1;
        break
    else
        Ctemp = Ctry;
    end
end


