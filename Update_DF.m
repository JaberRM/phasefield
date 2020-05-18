function [x, Tol, Iter, converged] = Update_DF(x,Phi,dPhidt,Param)


f = Calc_Force_DF(x,Phi,Param) + Param.K * dPhidt; %Eq 11

x0 = x + f * Param.dt;
converged = 0;
for Iter = 1:200
    
    f0 = Calc_Force_DF(x0,Phi,Param) + Param.K * dPhidt;
    
    x1 = x + (f + f0) * Param.dt/2;
    
    dx = x1 - x0;
    Tol = max(abs(dx(:))) ;
    if Tol < 1e-6
        x = x1;
        converged = 1;
        break
    else
        x0 = x1;
    end
end


