function [U, T, Tol, Iter,converged] = Update_Mech(U,Param)


dt = 0.001;


E = Calc_Strain(U,Param);

E(:,1:5,1) = E(:,1:5,1) - 1e-4;

T = Calc_Stress(E,Param);

f = Calc_Div_T(T,Param.dx,Param.dy,-1);

Utemp = U + f * dt;
converged = 0;

for Iter = 1:10000
    
    
    E = Calc_Strain(Utemp,Param);

    E(:,1:5,1) = E(:,1:5,1) - 1e-4;

    T = Calc_Stress(E,Param);
    ftemp = Calc_Div_T(T,Param.dx,Param.dy,-1);
    
    Utry = U + (f + ftemp) * dt/2;
    
    dU = Utry - Utemp;
    
    Tol = max(abs(dU(:))) 
    if Tol < 1e-6
        U = Utry;
        converged = 1;
        break
    else
        Utemp = Utry;
    end
end


