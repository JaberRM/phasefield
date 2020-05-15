function f = Calc_Force_CH(C,Phi,M,EPS,dx,dy,MethodFlag)

% G = G1(c) * f(Phi) + G2(c) * (1-f(Phi))
% G1 = (5*(C-0.0)).^2 + 10;
% G2 = (5*(C-1)).^2 + 10;

dGdC = (20*(C-0.0)) .* Calc_f(Phi) + (20*(C-1)) .*  (1-Calc_f(Phi)) ;

% dFdC = 1e-8*2*(C) ;

Del2C = Calc_Del2(C, dx,dy);

Mu = dGdC - 2 * EPS .* Del2C;

GradMu = Calc_Grad(Mu,dx,dy,MethodFlag);
% GradM  = Calc_Grad(M,dx,dy,MethodFlag);
Div_GradMu = Calc_Div(GradMu,dx,dy,MethodFlag);
% f = M .* Div_GradMu + dot(GradM, GradMu,3) - 2*Reaction_rate;
% f = M .* Div_GradMu + dot(GradM, GradMu,3);
f = M .* Div_GradMu ;
